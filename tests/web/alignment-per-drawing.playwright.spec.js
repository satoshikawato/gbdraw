// OV-3532 (OV-3504, OV-3506): each Linear drawing generates, saves and loads its
// own Similarity alignment. A drawing and its copy start aligned (the Gallery
// Session's plan); the copy's Reset alignment unaligns only the copy, so a
// Generate of the source still draws the source aligned, and after Save and Load
// a Generate of each draws it as it was saved (PD-OI-090).
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');

const GALLERY_SESSION = 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw';

/** @param {import('@playwright/test').Page} page @param {Buffer} bytes @param {string} name */
const importSession = (page, bytes, name) => evaluateWithRetainedPromise(page, async ({ base64, name }) => {
  const binary = atob(base64);
  const data = new Uint8Array(binary.length);
  for (let index = 0; index < binary.length; index += 1) data[index] = binary.charCodeAt(index);
  const result = await window.__GBDRAW_APP__.importSession({ target: { files: [new File([data], name)], value: 'selected' } });
  if (result.status !== 'ok') throw new Error(result.error?.stack || result.status);
  return result.status;
}, { base64: bytes.toString('base64'), name });

// The figure: the x translation each record's feature track is drawn with in
// the shown Result, by record key.
const recordOffsets = (page) => page.evaluate(() => Object.fromEntries([...window.__GBDRAW_APP__.svgContainer
  .querySelectorAll('svg g[data-gbdraw-slot-id="features"][data-record-key]')]
  .map((group) => [group.dataset.recordKey, Math.round(Number(group.dataset.recordTranslationX) * 10) / 10])));
const shownDrawingId = (page) => page.evaluate(async () => (await import('./js/state.js')).state.activeDrawingId.value);
const shownPlanGroup = (page) => page.evaluate(async () => (await import('./js/state.js')).state.similarityAlignmentPlan.value?.groupId ?? null);
const show = async (page, drawingId) => {
  await page.evaluate((id) => window.__GBDRAW_APP__.selectDrawing(id), drawingId);
  await expect.poll(() => shownDrawingId(page)).toBe(drawingId);
};
const generate = async (page) => {
  expect((await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis())).status).toBe('ok');
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing);
};

test('a Linear drawing and its copy each generate, save and load their own Similarity alignment', async ({ page }, testInfo) => {
  test.setTimeout(900000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(GALLERY_SESSION), 'BGC0000708-BGC0000713.gbdraw');
  const source = await shownDrawingId(page);
  const group = await shownPlanGroup(page);
  expect(group).toBeTruthy();
  const aligned = await recordOffsets(page);
  expect(Object.values(aligned).some((x) => x !== 0), 'the Gallery Session draws its records aligned').toBe(true);

  // OV-3532 (a): the copy's Reset alignment leaves the source aligned.
  const copy = await page.evaluate((id) => window.__GBDRAW_APP__.drawingList.duplicateDrawing(id).drawingId, source);
  await expect.poll(() => shownDrawingId(page)).toBe(copy);
  expect((await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.resetSimilarityAlignment('positions'))).status).toBe('ok');
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing);
  const unaligned = await recordOffsets(page);
  expect(Object.keys(unaligned)).toEqual(Object.keys(aligned));
  expect(unaligned).not.toEqual(aligned);
  expect(await shownPlanGroup(page)).toBeNull();

  await show(page, source);
  expect(await shownPlanGroup(page)).toBe(group);
  await generate(page);
  expect(await recordOffsets(page), 'a Generate of the source draws it aligned').toEqual(aligned);

  // OV-3532 (b): saved with the copy shown, loaded, each drawing generates as saved.
  await show(page, copy);
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'two alignments'; });
  const download = page.waitForEvent('download');
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const savedPath = testInfo.outputPath('two-alignments.gbdraw');
  await (await download).saveAs(savedPath);
  await importSession(page, readFileSync(savedPath), 'two-alignments.gbdraw');
  expect(await shownDrawingId(page)).toBe(copy);
  expect(await shownPlanGroup(page), 'the loaded copy has no plan').toBeNull();
  await generate(page);
  expect(await recordOffsets(page), 'a Generate of the loaded copy draws it unaligned').toEqual(unaligned);
  await show(page, source);
  expect(await shownPlanGroup(page)).toBe(group);
  await generate(page);
  expect(await recordOffsets(page), 'a Generate of the loaded source draws it aligned').toEqual(aligned);
});

// D-38 (PD-OI-039): a Similarity alignment review belongs to the drawing it was
// opened in. Over a copy of that drawing (the same records and plan) it stays
// open with its choices and marks nothing on the copy's figure; D-39: Apply
// there is disabled, and its reason names the drawing to show. Back on its own
// drawing, Apply aligns that drawing only. D-40: removing the review's drawing
// ends the review as Cancel does, and the Undo of the removal does not bring it
// back.
const openReview = async (page, proteinId) => {
  await page.evaluate(() => window.__GBDRAW_APP__.closeFeaturePopup());
  const featureId = await page.evaluate((id) => window.__GBDRAW_APP__.extractedFeatures
    .find((item) => item.protein_id === id || item.sourceProteinId === id)?.svg_id, proteinId);
  expect(featureId).toBeTruthy();
  await page.locator(`[data-gbdraw-feature-id="${featureId}"]`).first().dispatchEvent('click', { clientX: -100, clientY: -100 });
  await page.locator('.feature-popup[role="dialog"]').getByRole('button', { name: 'Review alignment options…' }).click();
  const dialog = page.getByRole('dialog', { name: 'Review alignment' });
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
  for (const key of await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDraft.rows
    .filter((row) => !row.choice).map((row) => row.recordKey))) {
    await dialog.locator(`[data-alignment-record-key="${key}"]`).getByRole('radio', { name: /Select .*bp, strand/ }).first().check();
  }
  return dialog;
};
const shownAlignmentState = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  return JSON.stringify([state.similarityAlignmentPlan.value, state.linearRecordTranslations.value]);
});
const applyOutcome = async (page, dialog) => {
  await dialog.getByRole('button', { name: 'Apply', exact: true }).click();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
  await expect(dialog).toBeHidden({ timeout: 180000 });
  return page.evaluate(() => ({ summary: Boolean(window.__GBDRAW_APP__.similarityAlignmentSummary),
    error: window.__GBDRAW_APP__.errorLog?.code ?? null }));
};

test('an open alignment review stays with its drawing and applies only there', async ({ page }) => {
  test.setTimeout(900000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(GALLERY_SESSION), 'BGC0000708-BGC0000713.gbdraw');
  const source = await shownDrawingId(page);
  const copy = await page.evaluate((id) => window.__GBDRAW_APP__.drawingList.duplicateDrawing(id).drawingId, source);
  await expect.poll(() => shownDrawingId(page)).toBe(copy);
  const copyBefore = { figure: await recordOffsets(page), alignment: await shownAlignmentState(page) };
  await show(page, source);
  const sourceBefore = { figure: await recordOffsets(page), alignment: await shownAlignmentState(page) };

  const sourceName = await page.evaluate((id) => window.__GBDRAW_APP__.drawingBar.nameOf(id), source);
  const dialog = await openReview(page, 'CAG38712.1');
  const apply = dialog.getByRole('button', { name: 'Apply', exact: true });
  await show(page, copy);
  await expect(dialog, 'the review stays open over the copy').toBeVisible();
  await expect(page.locator('[data-similarity-alignment-canvas]'), 'and marks nothing on its figure').toHaveCount(0);
  await expect(apply, 'Apply waits for the review\'s drawing').toBeDisabled();
  await expect(page.locator('#similarity-alignment-apply-reason')).toHaveText(`Show ${sourceName} to apply this alignment.`);
  expect({ figure: await recordOffsets(page), alignment: await shownAlignmentState(page) }).toEqual(copyBefore);
  await show(page, source);
  await expect(dialog).toBeVisible();
  await expect(apply).toBeEnabled();
  expect(await applyOutcome(page, dialog), 'Apply on the review\'s own drawing').toMatchObject({ summary: true });
  await show(page, copy);
  expect({ figure: await recordOffsets(page), alignment: await shownAlignmentState(page) }, 'the copy is unchanged')
    .toEqual(copyBefore);

  // D-40: the review's drawing is deleted while the copy is shown.
  await show(page, source);
  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('orthogroups'));
  await openReview(page, 'CAG38712.1');
  await show(page, copy);
  expect((await page.evaluate((id) => window.__GBDRAW_APP__.drawingList.deleteDrawing(id), source)).status).toBe('ok');
  await expect(dialog, 'the review leaves with its drawing').toBeHidden();
  expect(await page.evaluate(() => ({ open: window.__GBDRAW_APP__.showRightDrawer, tab: window.__GBDRAW_APP__.rightDrawerTab })),
    'the Editor is given back on its tab').toEqual({ open: true, tab: 'orthogroups' });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.undoHistory())).toBe(true);
  await show(page, source);
  await expect(dialog, 'the Undo of the Delete brings the drawing back without the review').toBeHidden();
});
