const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { mkdirSync, readFileSync } = require('node:fs');
const { join } = require('node:path');
const { gunzipSync } = require('node:zlib');

const importSession = async (page, bytes, name) => page.evaluate(async ({ bytes, name }) => {
  const file = new File([new Uint8Array(bytes)], name);
  const result = await window.__GBDRAW_APP__.importSession({
    target: { files: [file], value: 'selected' }
  });
  if (result.status !== 'ok') throw new Error(result.error?.stack || result.status);
  return result.status;
}, { bytes: [...bytes], name });

const captureUnhandledRejections = (page) => page.addInitScript(() => {
  window.__GBDRAW_UNHANDLED_REJECTIONS__ = [];
  window.addEventListener('unhandledrejection', (event) => {
    window.__GBDRAW_UNHANDLED_REJECTIONS__.push(String(event.reason?.message || event.reason));
  });
});
const expectNoUnhandledRejections = async (page) => {
  expect(await page.evaluate(() => window.__GBDRAW_UNHANDLED_REJECTIONS__)).toEqual([]);
};

const artifactSnapshot = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const history = window.__GBDRAW_HISTORY__;
  return {
    plan: JSON.parse(JSON.stringify(state.similarityAlignmentPlan.value)),
    request: structuredClone(state.lastCommittedRequest?.value || null),
    results: state.results.value.map(({ name, content }) => ({ name, content })),
    history: [history.getUndoCount(), history.getRedoCount(), history.revision.value]
  };
});

test('Similarity alignment UI completes exact-reference, ambiguity, focus, summary, and narrow journeys', async ({ page, browser }, testInfo) => {
  test.setTimeout(600000);
  const pageErrors = [];
  const consoleErrors = [];
  page.on('pageerror', (error) => pageErrors.push(String(error?.message || error)));
  page.on('console', (message) => {
    if (message.type() === 'error') consoleErrors.push(message.text());
  });
  await captureUnhandledRejections(page);
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
  const b0Feature = page.locator(`[data-gbdraw-feature-id="${b0FeatureId}"]`).first();
  await b0Feature.dispatchEvent('click', { clientX: -100, clientY: -100 });
  const popup = page.locator('.feature-popup[role="dialog"]');
  await expect(popup).toBeVisible();
  await expect(popup).toContainText('b0');
  await expect(popup.getByRole('button', { name: 'Reset', exact: true })).toHaveCount(0);
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
  const similarityTab = drawer.getByRole('button', { name: 'Similarity groups' });
  await expect(similarityTab).toBeEnabled();
  await similarityTab.click();
  const inspector = drawer.locator('[data-similarity-alignment-plan-inspector]');
  await expect(inspector).toBeVisible();
  await expect(inspector).toContainText('Exact reference:');
  await expect(inspector).toContainText('Selected by user');
  await expect(inspector).toContainText('Exact reference');
  await page.screenshot({ path: testInfo.outputPath('narrow-final.png'), fullPage: true });

  const savedState = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return {
      plan: JSON.parse(JSON.stringify(state.similarityAlignmentPlan.value)),
      translations: JSON.parse(JSON.stringify(state.linearRecordTranslations.value)),
      orientations: state.linearSeqs.map(({ uid, region_reverse: reverse }) => ({
        recordKey: uid, reverse: Boolean(reverse)
      })),
      results: state.results.value.map(({ name, content }) => ({ name, content }))
    };
  });
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 's06-active-alignment'; });
  const downloadPromise = page.waitForEvent('download', { timeout: 180000 });
  await page.evaluate(() => window.__GBDRAW_APP__.saveSessionWithTitle());
  const download = await downloadPromise;
  const savedPath = testInfo.outputPath('s06-active-alignment.gbdraw-session.json.gz');
  await download.saveAs(savedPath);

  const freshContext = await browser.newContext();
  const freshPage = await freshContext.newPage();
  await captureUnhandledRejections(freshPage);
  const freshErrors = [];
  freshPage.on('pageerror', (error) => freshErrors.push(String(error?.message || error)));
  await freshPage.goto(new URL('/gbdraw/web/index.html', page.url()).href, {
    waitUntil: 'domcontentloaded'
  });
  await freshPage.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(
    freshPage,
    readFileSync(savedPath),
    's06-active-alignment.gbdraw-session.json.gz'
  );
  const freshState = await freshPage.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return {
      plan: JSON.parse(JSON.stringify(state.similarityAlignmentPlan.value)),
      translations: JSON.parse(JSON.stringify(state.linearRecordTranslations.value)),
      orientations: state.linearSeqs.map(({ uid, region_reverse: reverse }) => ({
        recordKey: uid, reverse: Boolean(reverse)
      })),
      results: state.results.value.map(({ name, content }) => ({ name, content })),
      workers: structuredClone(window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__ || {
        constructions: 0, instances: []
      })
    };
  });
  expect(freshState).toMatchObject({
    plan: savedState.plan,
    translations: savedState.translations,
    orientations: savedState.orientations,
    results: savedState.results
  });
  expect(freshState.workers.constructions).toBe(0);
  expect(freshState.workers.instances).toHaveLength(0);

  await freshPage.evaluate(() => {
    window.__GBDRAW_APP__.showRightDrawer = true;
    window.__GBDRAW_APP__.rightDrawerTab = 'orthogroups';
  });
  const freshInspector = freshPage.locator('[data-similarity-alignment-plan-inspector]');
  await expect(freshInspector).toBeVisible();
  await expect(freshInspector).toContainText('Selected by user');
  const historyBeforeReset = await freshPage.evaluate(
    () => window.__GBDRAW_HISTORY__.getUndoCount()
  );
  await freshPage.locator('[data-similarity-alignment-reset]').click();
  await expect.poll(() => freshPage.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return {
      plan: state.similarityAlignmentPlan.value,
      translations: JSON.parse(JSON.stringify(state.linearRecordTranslations.value)),
      history: window.__GBDRAW_HISTORY__.getUndoCount()
    };
  }), { timeout: 180000 }).toEqual({
    plan: null,
    translations: savedState.translations,
    history: historyBeforeReset + 1
  });
  await freshPage.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(() => freshPage.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return JSON.parse(JSON.stringify(state.similarityAlignmentPlan.value));
  })).toEqual(savedState.plan);
  await freshPage.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(() => freshPage.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return state.similarityAlignmentPlan.value;
  })).toBeNull();
  expect(pageErrors).toEqual([]);
  expect(consoleErrors).toEqual([]);
  expect(freshErrors).toEqual([]);
  await expectNoUnhandledRejections(page);
  await expectNoUnhandledRejections(freshPage);
  await freshContext.close();
});

test('one usable member resolves without a dialog and keeps the orient plan inspectable', async ({ page }, testInfo) => {
  test.setTimeout(180000);
  await captureUnhandledRejections(page);
  const directory = testInfo.outputPath('one-candidate-fixture');
  mkdirSync(directory, { recursive: true });
  execFileSync('python', ['-c',
    'import sys; from pathlib import Path; from tests.test_api_session import _record_local_collinear_session; _record_local_collinear_session(Path(sys.argv[1]))',
    directory
  ], { cwd: process.cwd(), stdio: 'pipe' });
  const source = readFileSync(join(directory, 'mixed.gbdraw-session.json'), 'utf8');
  const errors = [];
  page.on('pageerror', (error) => errors.push(String(error?.message || error)));
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const session = Buffer.from(await page.evaluate(async (raw) => {
    const { projectCanonicalSessionRequest } = await import('./js/services/session-request.js');
    const document = JSON.parse(raw);
    document.renderRequest.comparisons = [];
    document.editorState.featureCatalog.items.forEach((item) => {
      const secondCopy = item.biologicalFeatures.find((feature) => (
        feature.qualifiers?.protein_id?.[0] === 'a1'
      ))?.biologicalFeatureId;
      const group = item.orthogroups.find(({ id }) => id === 'og_1');
      if (!secondCopy || !group) throw new Error('The one-candidate fixture is incomplete.');
      group.members = group.members.filter(({ biologicalFeatureId }) => (
        biologicalFeatureId !== secondCopy
      ));
      group.orthologEdges = [];
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
  await importSession(page, session, 'one-candidate.gbdraw-session.json');

  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) {
    await page.locator('.drawer-toggle').click();
  }
  const drawer = page.locator('.right-drawer');
  await drawer.getByRole('button', { name: 'Similarity groups' }).click();
  await drawer.locator('button').filter({
    has: page.locator('.font-mono', { hasText: /^og_1$/ })
  }).click();
  const reference = await page.evaluate(() => (
    window.__GBDRAW_APP__.similarityAlignmentDrawerReferenceOptions('og_1')
      .find(({ anchor }) => anchor.recordKey === 'record_b')?.key || ''
  ));
  expect(reference).not.toBe('');
  await drawer.getByLabel('Exact reference record and feature').selectOption(reference);
  await drawer.getByRole('button', { name: 'Align & orient', exact: true }).click();
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return state.similarityAlignmentPlan.value?.mode || null;
  }), { timeout: 120000 }).toBe('position_and_orientation');
  await expect(page.getByRole('dialog', { name: 'Select alignment anchors' })).toBeHidden();
  await expect(page.locator('[data-similarity-alignment-summary]')).toContainText('1 aligned');
  const similarityTab = drawer.getByRole('button', { name: 'Similarity groups' });
  await expect(similarityTab).toBeEnabled();
  await similarityTab.click();
  await expect(drawer.locator('[data-similarity-alignment-plan-inspector]')).toBeVisible();
  expect(errors).toEqual([]);
  await expectNoUnhandledRejections(page);
});

test('released Session 44 schema 7 and catalog 4 load and save through the current writer', async ({ page }, testInfo) => {
  await captureUnhandledRejections(page);
  const source = readFileSync(
    'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'
  );
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, source, 'HmmtDNA_basic_circular.gbdraw-session.json');
  const workers = await page.evaluate(() => window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__ || {
    constructions: 0, instances: []
  });
  expect(workers.constructions).toBe(0);
  expect(workers.instances).toHaveLength(0);

  const downloadPromise = page.waitForEvent('download');
  await page.evaluate(() => {
    window.__GBDRAW_APP__.sessionTitle = 'released-session-44-schema-7';
    window.__GBDRAW_APP__.saveSessionWithTitle();
  });
  const download = await downloadPromise;
  const savedPath = testInfo.outputPath('released-session-44-current.gbdraw-session.json.gz');
  await download.saveAs(savedPath);
  const saved = JSON.parse(gunzipSync(readFileSync(savedPath)).toString('utf8'));
  expect(saved.version).toBe(44);
  expect(saved.renderRequest.schema).toBe(8);
  expect(saved.editorState.featureCatalog.schema).toBe(4);
  const previous = JSON.parse(source.toString('utf8'));
  expect(saved.results.map(({ name }) => name)).toEqual(previous.results.map(({ name }) => name));
  const previewEquivalent = await page.evaluate(({ before, after }) => {
    const parser = new DOMParser();
    const normalized = (element) => ({
      tag: element.tagName,
      attributes: [...element.attributes]
        .filter(({ name }) => name !== 'baseProfile' && !name.startsWith('xmlns')) // Sanitizer removes unused declarations.
        .map(({ name, value }) => [name, value]).sort(([a], [b]) => a.localeCompare(b)),
      text: [...element.childNodes]
        .filter(({ nodeType }) => nodeType === Node.TEXT_NODE)
        .map(({ textContent }) => textContent.trim()).filter(Boolean),
      children: [...element.children].map(normalized)
    });
    const original = parser.parseFromString(before, 'image/svg+xml').documentElement;
    const retained = parser.parseFromString(after, 'image/svg+xml').documentElement;
    return JSON.stringify(normalized(original)) === JSON.stringify(normalized(retained));
  }, { before: previous.results[0].content, after: saved.results[0].content });
  expect(previewEquivalent).toBe(true);
  await expectNoUnhandledRejections(page);
});

test('released v40 alignment materializes by stable feature identity without a Worker', async ({ page }) => {
  await captureUnhandledRejections(page);
  const source = readFileSync(
    'tests/fixtures/sessions/BGC0000708-BGC0000713.v40-schema5.json',
    'utf8'
  );
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const observed = await page.evaluate(async (raw) => {
    const session = JSON.parse(raw);
    const { materializeLegacySimilarityAlignment } = await import(
      './js/services/legacy-similarity-alignment.js'
    );
    const { materializeRecordTranslations } = await import(
      './js/app/legend-layout/composition-actions.js'
    );
    const target = session.renderRequest.comparisons.find(
      (comparison) => comparison.kind === 'generatedProteinComparison'
    ).settings.alignOrthogroupFeature;
    const plan = materializeLegacySimilarityAlignment({
      target,
      records: session.renderRequest.records,
      featureCatalog: session.editorState.featureCatalog,
      legacyOrthogroupState: session.orthogroupState
    });
    const svg = new DOMParser().parseFromString(
      session.results[0].content,
      'image/svg+xml'
    ).documentElement;
    const legacyRecordGroups = Array.from(svg.querySelectorAll(
      '[data-gbdraw-composition-role="primary"][data-gbdraw-record-id]'
    )).filter((group) => !group.hasAttribute('data-gbdraw-definition-part'));
    legacyRecordGroups.reverse().forEach((group) => group.parentElement.append(group));
    const recordKeys = session.renderRequest.records.map((record) => record.recordKey);
    const base = recordKeys.map((recordKey) => ({ recordKey, x: 0, y: 0 }));
    const translations = materializeRecordTranslations(svg, base, recordKeys, plan);
    const malformed = structuredClone(plan);
    malformed.records[0].anchor.biologicalFeatureId = 'missing-feature';
    malformed.reference.biologicalFeatureId = 'missing-feature';
    let malformedError = '';
    try {
      materializeRecordTranslations(svg, base, recordKeys, malformed);
    } catch (error) {
      malformedError = String(error?.message || error);
    }
    return {
      translations,
      malformedError,
      workers: structuredClone(window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__ || {
        constructions: 0,
        instances: []
      })
    };
  }, source);
  expect(observed.translations).toEqual([
    { recordKey: 'record-1', x: 577, y: 0 },
    { recordKey: 'record-2', x: 358.65497562715484, y: 0 },
    { recordKey: 'record-3', x: 697.3384456862045, y: 0 },
    { recordKey: 'record-4', x: 369.1372805453176, y: 0 },
    { recordKey: 'record-5', x: 588.9882693298457, y: 0 }
  ]);
  expect(observed.malformedError).toMatch(/cannot bind alignment record "record-1"/);
  expect(observed.workers.constructions).toBe(0);
  expect(observed.workers.instances).toHaveLength(0);
  await expectNoUnhandledRejections(page);
});
