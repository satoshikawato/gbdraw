const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { mkdirSync, readFileSync } = require('node:fs');
const { join } = require('node:path');

const importSession = async (page, bytes, name) => page.evaluate(async ({ bytes, name }) => {
  const file = new File([new Uint8Array(bytes)], name);
  const result = await window.__GBDRAW_APP__.importSession({
    target: { files: [file], value: 'selected' }
  });
  if (result.status !== 'ok') throw new Error(result.error?.stack || result.status);
  return result.status;
}, { bytes: [...bytes], name });

const artifactSnapshot = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const history = window.__GBDRAW_HISTORY__;
  return {
    plan: structuredClone(state.similarityAlignmentPlan.value),
    request: structuredClone(state.lastCommittedRequest?.value || null),
    results: state.results.value.map(({ name, content }) => ({ name, content })),
    history: [history.getUndoCount(), history.getRedoCount(), history.revision.value]
  };
});

test('Similarity alignment UI completes exact-reference, ambiguity, focus, summary, and narrow journeys', async ({ page }, testInfo) => {
  test.setTimeout(300000);
  const directory = testInfo.outputPath('inparalog-fixture');
  mkdirSync(directory, { recursive: true });
  execFileSync('python', ['-c',
    'import sys; from pathlib import Path; from tests.test_api_session import _record_local_collinear_session; _record_local_collinear_session(Path(sys.argv[1]))',
    directory
  ], { cwd: process.cwd(), stdio: 'pipe' });
  const source = readFileSync(join(directory, 'mixed.gbdraw-session.json'), 'utf8');

  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const session = Buffer.from(await page.evaluate(async (raw) => {
    const { projectCanonicalSessionRequest } = await import('./js/services/session-request.js');
    const document = JSON.parse(raw);
    document.renderRequest.comparisons = [];
    document.editorState.featureCatalog.items.forEach((item) => {
      item.orthogroups.forEach((group) => {
        if (group.id === 'og_1') group.orthologEdges = [];
      });
    });
    document.config = projectCanonicalSessionRequest({
      renderRequest: document.renderRequest,
      resources: document.resources,
      webFiles: document.webFiles
    }).config;
    document.config.linearComparisonPlan = {
      mode: 'adjacent', defaultSource: 'losat', edges: []
    };
    document.config.losatProgram = 'blastp';
    document.config.losat = { blastp: { mode: 'orthogroup' } };
    return JSON.stringify(document);
  }, source));
  await importSession(page, session, 'similarity-alignment-inparalog.gbdraw-session.json');

  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) {
    await page.locator('.drawer-toggle').click();
  }
  const drawer = page.locator('.right-drawer');
  await drawer.getByRole('button', { name: 'Similarity groups' }).click();
  await drawer.locator('button').filter({
    has: page.locator('.font-mono', { hasText: /^og_1$/ })
  }).click();

  const align = drawer.getByRole('button', { name: 'Align', exact: true });
  const alignOrient = drawer.getByRole('button', { name: 'Align & orient', exact: true });
  const referenceSelect = drawer.getByLabel('Exact reference record and feature');
  await expect(align).toBeDisabled();
  await expect(alignOrient).toBeDisabled();
  await expect(page.locator('#similarity-alignment-drawer-reason')).toContainText(
    'Select an exact reference record and feature'
  );
  const b0Reference = await page.evaluate(() => (
    window.__GBDRAW_APP__.similarityAlignmentDrawerReferenceOptions('og_1')
      .find(({ anchor }) => anchor.recordKey === 'record_b')?.key || ''
  ));
  expect(b0Reference).not.toBe('');
  await referenceSelect.selectOption(b0Reference);
  await expect(align).toBeEnabled();
  await expect(alignOrient).toBeEnabled();

  const beforeCancel = await artifactSnapshot(page);
  await alignOrient.click();
  const dialog = page.getByRole('dialog', { name: 'Select alignment anchors' });
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await expect(dialog).toHaveAttribute('aria-describedby', 'similarity-alignment-dialog-description');
  await expect(dialog).toContainText('inparalog');
  await expect(dialog).toContainText(/\d+\.\.\d+ · strand [+-]/);
  await expect(dialog).toContainText('Direct evidence: None');
  await expect(dialog).not.toContainText(/score/i);
  const apply = dialog.getByRole('button', { name: 'Apply', exact: true });
  await expect(apply).toBeDisabled();
  await expect(apply).toHaveAttribute('aria-describedby', 'similarity-alignment-apply-reason');
  await expect(page.locator('#similarity-alignment-apply-reason')).toContainText(
    'require Select or Skip'
  );
  expect(await dialog.evaluate((element) => element.contains(document.activeElement))).toBe(true);
  await page.screenshot({
    path: testInfo.outputPath('desktop-ambiguity.png'),
    fullPage: true
  });

  const firstCandidate = dialog.getByRole('radio', { name: /Select feature/ }).first();
  const candidateAnchor = await page.evaluate(() => (
    window.__GBDRAW_APP__.similarityAlignmentDraft.ambiguities[0].candidates[0].anchor
  ));
  const previewState = async () => page.evaluate((anchor) => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((item) => (
      item.recordKey === anchor.recordKey
      && item.biologicalFeatureId === anchor.biologicalFeatureId
    ));
    const candidateElements = Array.from(document.querySelectorAll(
      `[data-gbdraw-feature-id="${CSS.escape(feature.svg_id)}"]`
    ));
    const otherElements = Array.from(document.querySelectorAll('[data-gbdraw-feature-id]'))
      .filter((element) => !candidateElements.includes(element));
    return {
      candidate: candidateElements.map((element) => element.style.opacity),
      otherHighlighted: otherElements.filter((element) => element.style.opacity === '0.7').length
    };
  }, candidateAnchor);
  await firstCandidate.locator('..').hover();
  await expect.poll(previewState).toMatchObject({ otherHighlighted: 0 });
  expect((await previewState()).candidate).toContain('0.7');
  await dialog.getByRole('heading', { name: 'Select alignment anchors' }).hover();
  await expect.poll(async () => (await previewState()).candidate.includes('0.7')).toBe(false);

  const cancel = dialog.getByRole('button', { name: 'Cancel', exact: true });
  await cancel.focus();
  await page.keyboard.press('Tab');
  expect(await dialog.evaluate((element) => element.contains(document.activeElement))).toBe(true);
  await page.keyboard.press('Shift+Tab');
  await dialog.getByRole('radio', { name: /Skip record/ }).check();
  await expect(dialog).toContainText('Selected: Skip');
  await expect(apply).toBeEnabled();
  await dialog.locator('input[type="radio"]:checked').press('Escape');
  await expect(dialog).toBeHidden();
  await expect(alignOrient).toBeFocused();
  expect(await artifactSnapshot(page)).toEqual(beforeCancel);

  await page.setViewportSize({ width: 720, height: 740 });
  await drawer.getByRole('button', { name: 'Similarity groups' }).click();
  const b0FeatureId = await page.evaluate(() => {
    const feature = window.__GBDRAW_APP__.extractedFeatures.find(
      ({ protein_id: proteinId }) => proteinId === 'b0'
    );
    return feature.svg_id;
  });
  await page.locator(`[data-gbdraw-feature-id="${b0FeatureId}"]`).first().dispatchEvent('click', {
    clientX: 240,
    clientY: 260
  });
  const popup = page.locator('.feature-popup[role="dialog"]');
  await expect(popup).toBeVisible();
  const popupAlign = popup.getByRole('button', { name: 'Align', exact: true });
  const popupAlignOrient = popup.getByRole('button', { name: 'Align & orient', exact: true });
  await expect(popupAlign).toHaveAttribute('title', /keep every record orientation unchanged/);
  await expect(popupAlignOrient).toHaveAttribute('title', /reverse targets/);
  await popupAlign.scrollIntoViewIfNeeded();
  await popupAlign.click();
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await dialog.getByRole('radio', { name: /Select feature/ }).first().check();
  await expect(apply).toBeEnabled();
  await apply.scrollIntoViewIfNeeded();
  const applyBox = await apply.boundingBox();
  expect(applyBox).not.toBeNull();
  expect(applyBox.x).toBeGreaterThanOrEqual(0);
  expect(applyBox.x + applyBox.width).toBeLessThanOrEqual(720);
  expect(applyBox.y + applyBox.height).toBeLessThanOrEqual(740);
  await page.screenshot({
    path: testInfo.outputPath('narrow-ambiguity.png'),
    fullPage: true
  });
  await dialog.locator('input[type="radio"]:checked').press('Escape');
  await expect(dialog).toBeHidden();
  await expect(popupAlign).toBeFocused();
  await popup.getByRole('button', { name: 'Close feature popup' }).click();

  await page.evaluate(() => {
    window.__GBDRAW_APP__.linearComparisonPlan.mode = 'none';
    window.__GBDRAW_APP__.losatProgram = 'blastn';
  });
  await drawer.getByRole('button', { name: 'Similarity groups' }).click();
  const undoCountBeforeApply = await page.evaluate(
    () => window.__GBDRAW_HISTORY__.getUndoCount()
  );
  await align.click();
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await dialog.getByRole('radio', { name: /Select feature/ }).first().check();
  await expect(apply).toBeEnabled();
  await apply.click();

  const summary = page.locator('[data-similarity-alignment-summary]');
  await expect(summary).toBeVisible({ timeout: 180000 });
  await expect(summary).toContainText('1 aligned');
  await expect(summary).toContainText('1 unchanged');
  await expect(summary).toContainText('0 explicitly skipped');
  await expect(summary).toContainText('0 no-candidate');
  await expect(summary).toContainText('0 reversed');
  expect(await page.evaluate(
    () => window.__GBDRAW_HISTORY__.getUndoCount()
  )).toBe(undoCountBeforeApply + 1);

  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) {
    await page.locator('.drawer-toggle').click();
  }
  await page.evaluate(() => {
    window.__GBDRAW_APP__.rightDrawerTab = 'orthogroups';
  });
  const inspector = drawer.locator('[data-similarity-alignment-plan-inspector]');
  await expect(inspector).toBeVisible();
  await expect(inspector).toContainText('Exact reference:');
  await expect(inspector).toContainText('Selected by user');
  await expect(inspector).toContainText('Exact reference');
  await page.screenshot({ path: testInfo.outputPath('narrow-final.png'), fullPage: true });
});
