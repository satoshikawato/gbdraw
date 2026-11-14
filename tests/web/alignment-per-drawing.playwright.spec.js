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
