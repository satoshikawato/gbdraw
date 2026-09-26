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
  const { getCommittedCanonicalSession } = await import('./js/services/config.js');
  return {
    session: JSON.parse(JSON.stringify(getCommittedCanonicalSession())),
    plan: JSON.parse(JSON.stringify(state.similarityAlignmentPlan.value)),
    receipt: JSON.parse(JSON.stringify(state.similarityAlignmentResetReceipt.value)),
    request: structuredClone(state.lastCommittedRequest?.value || null),
    results: state.results.value.map(({ name, content }) => ({ name, content })),
    history: [history.getUndoCount(), history.getRedoCount(), history.revision.value]
  };
});

const loadAmbiguousSession = async (page, testInfo) => {
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
};

test('alignment canvas guide and candidates stay transient and share palette choices', async ({ page }, testInfo) => {
  test.setTimeout(180000);
  await loadAmbiguousSession(page, testInfo);
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 's03-canvas-overlay'; });
  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) {
    await page.locator('.drawer-toggle').click();
  }
  const drawer = page.locator('.right-drawer');
  await drawer.getByRole('button', { name: 'Similarity groups' }).click();
  await drawer.locator('button').filter({
    has: page.locator('.font-mono', { hasText: /^og_1$/ })
  }).click();
  const referenceSelect = drawer.getByLabel('Exact reference record and feature');
  const referenceKey = await page.evaluate(() => (
    window.__GBDRAW_APP__.similarityAlignmentDrawerReferenceOptions('og_1')
      .find(({ anchor }) => anchor.recordKey === 'record_b')?.key || ''
  ));
  await referenceSelect.selectOption(referenceKey);
  const align = drawer.getByRole('button', { name: 'Align…', exact: true });
  const before = await artifactSnapshot(page);
  await align.click();
  const dialog = page.getByRole('dialog', { name: 'Select alignment anchors' });
  await expect(dialog).toBeVisible({ timeout: 180000 });
  const overlay = page.locator('[data-similarity-alignment-canvas]');
  const guide = overlay.locator('.gbdraw-alignment-guide');
  const badges = overlay.locator('.gbdraw-alignment-badge');
  await expect(guide).toBeVisible();
  await expect(badges).toHaveCount(2);
  await expect(badges.nth(0)).toHaveText('1');
  await expect(badges.nth(1)).toHaveText('2');
  await expect(dialog.getByRole('radio', { name: /Select .*bp, strand/ }).first().locator('xpath=ancestor::label')).toContainText('1');

  const geometry = () => page.evaluate(async () => {
    const { getFeatureFillElements } = await import('./js/app/feature-dom.js');
    const app = window.__GBDRAW_APP__;
    const svg = document.querySelector('.gbdraw-preview-surface svg');
    const draft = app.similarityAlignmentDraft;
    const center = (anchor) => {
      const feature = app.extractedFeatures.find((item) => (
        item.recordKey === anchor.recordKey
        && item.biologicalFeatureId === anchor.biologicalFeatureId
      ));
      const [element] = getFeatureFillElements(svg, feature.svg_id);
      const rect = element.getBoundingClientRect();
      return { x: (rect.left + rect.right) / 2, y: (rect.top + rect.bottom) / 2 };
    };
    const reference = center(draft.response.reference);
    const line = document.querySelector('.gbdraw-alignment-guide');
    const candidateErrors = draft.rows[0].candidates.map((candidate, index) => {
      const badge = document.querySelectorAll('.gbdraw-alignment-badge')[index];
      if (badge.hidden) return null;
      const point = center(candidate.anchor);
      const rect = badge.getBoundingClientRect();
      return Math.max(Math.abs(point.x - (rect.left + rect.right) / 2),
        Math.abs(point.y - (rect.top + rect.bottom) / 2));
    }).filter((value) => value !== null);
    return {
      guideError: Math.abs(reference.x - line.getBoundingClientRect().left),
      candidateErrors,
      visibleBadges: candidateErrors.length
    };
  });
  await expect.poll(async () => (await geometry()).visibleBadges).toBe(2);
  expect((await geometry()).guideError).toBeLessThan(2);
  expect(Math.max(...(await geometry()).candidateErrors)).toBeLessThan(2);
  const geometryIsAligned = async () => {
    const measured = await geometry();
    return measured.visibleBadges > 0 && measured.guideError < 2
      && measured.candidateErrors.every((error) => error < 2);
  };
  await page.screenshot({ path: testInfo.outputPath('canvas-overlay.png'), fullPage: true });

  const header = await dialog.locator('header').boundingBox();
  await page.mouse.move(header.x + 25, header.y + 20);
  await page.mouse.down();
  await page.mouse.move(50, 45, { steps: 4 });
  await page.mouse.up();
  const firstRadio = dialog.getByRole('radio', { name: /Select .*bp, strand/ }).first();
  const secondRadio = dialog.getByRole('radio', { name: /Select .*bp, strand/ }).nth(1);
  const firstRow = firstRadio.locator('xpath=ancestor::label');
  await badges.first().hover();
  await expect(firstRow).toHaveClass(/ring-2/);
  const firstSvgId = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const anchor = app.similarityAlignmentDraft.rows[0].candidates[0].anchor;
    return app.extractedFeatures.find((item) => item.recordKey === anchor.recordKey
      && item.biologicalFeatureId === anchor.biologicalFeatureId).svg_id;
  });
  const firstFeature = page.locator('[data-gbdraw-feature-id="' + firstSvgId + '"]').first();
  await firstFeature.dispatchEvent('mouseover');
  await expect(firstRow).toHaveClass(/ring-2/);
  await firstFeature.dispatchEvent('click', { clientX: -100, clientY: -100 });
  await expect(firstRadio).toBeChecked();
  await expect(page.locator('.feature-popup[role="dialog"]')).toHaveCount(0);
  await badges.nth(1).click();
  await expect(secondRadio).toBeChecked();
  const secondRow = secondRadio.locator('xpath=ancestor::label');
  await firstFeature.dispatchEvent('mouseover');
  await secondRow.dispatchEvent('mouseenter');
  await expect.poll(() => firstFeature.evaluate((element) => element.style.opacity))
    .not.toBe('0.7');
  await secondRow.dispatchEvent('mouseleave');
  await expect.poll(() => firstFeature.evaluate((element) => element.style.opacity))
    .toBe('0.7');
  await expect(secondRadio).toBeChecked();
  const referenceSvgId = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const anchor = app.similarityAlignmentDraft.response.reference;
    return app.extractedFeatures.find((item) => item.recordKey === anchor.recordKey
      && item.biologicalFeatureId === anchor.biologicalFeatureId).svg_id;
  });
  await page.locator('[data-gbdraw-feature-id="' + referenceSvgId + '"]').first()
    .dispatchEvent('click', { clientX: -100, clientY: -100 });
  await expect(page.locator('.feature-popup[role="dialog"]')).toBeVisible();
  await page.evaluate(() => window.__GBDRAW_APP__.closeFeaturePopup());
  expect(await artifactSnapshot(page)).toEqual(before);

  const canvas = page.locator('.gbdraw-preview-surface').locator('..');
  const zoomBefore = await page.evaluate(() => window.__GBDRAW_APP__.zoom);
  await canvas.dispatchEvent('wheel', { deltaY: -80 });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.zoom)).toBeGreaterThan(zoomBefore);
  await expect.poll(geometryIsAligned).toBe(true);
  const panBefore = await page.evaluate(() => window.__GBDRAW_APP__.canvasPan.x);
  const canvasBox = await canvas.boundingBox();
  await page.mouse.move(canvasBox.x + 240, canvasBox.y + 220);
  await page.mouse.down();
  await page.mouse.move(canvasBox.x + 270, canvasBox.y + 220, { steps: 3 });
  await page.mouse.up();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.canvasPan.x))
    .toBeGreaterThan(panBefore);
  await expect.poll(geometryIsAligned).toBe(true);
  const scrollAfter = await canvas.evaluate((element) => { element.scrollLeft += 40; return element.scrollLeft; });
  expect(scrollAfter).toBeGreaterThan(0);
  await expect.poll(geometryIsAligned).toBe(true);
  await page.setViewportSize({ width: 1600, height: 900 });
  await expect.poll(geometryIsAligned).toBe(true);

  await page.evaluate((svgId) => {
    document.querySelectorAll('[data-gbdraw-feature-id="' + CSS.escape(svgId) + '"]')
      .forEach((element) => element.setAttribute('display', 'none'));
    window.dispatchEvent(new Event('resize'));
  }, firstSvgId);
  await expect(badges.first()).toBeHidden();
  await firstRadio.check();
  await expect(firstRadio).toBeChecked();
  await page.evaluate((svgId) => {
    document.querySelectorAll('[data-gbdraw-feature-id="' + CSS.escape(svgId) + '"]')
      .forEach((element) => element.removeAttribute('display'));
    window.dispatchEvent(new Event('resize'));
  }, firstSvgId);
  expect(JSON.stringify(await artifactSnapshot(page))).not.toContain('gbdraw-alignment-');
  expect(await page.evaluate(() => document.querySelector('.gbdraw-preview-surface svg').outerHTML))
    .not.toContain('gbdraw-alignment-');
  const svgDownloadPromise = page.waitForEvent('download');
  await page.evaluate(() => window.__GBDRAW_APP__.downloadSVG());
  const svgDownload = await svgDownloadPromise;
  const svgPath = testInfo.outputPath('canvas-overlay-export.svg');
  await svgDownload.saveAs(svgPath);
  expect(readFileSync(svgPath, 'utf8')).not.toContain('gbdraw-alignment-');
  const sessionDownloadPromise = page.waitForEvent('download');
  await page.evaluate(() => window.__GBDRAW_APP__.saveSessionWithTitle());
  const sessionDownload = await sessionDownloadPromise;
  const sessionPath = testInfo.outputPath('canvas-overlay-session.json.gz');
  await sessionDownload.saveAs(sessionPath);
  expect(gunzipSync(readFileSync(sessionPath)).toString('utf8'))
    .not.toContain('gbdraw-alignment-');

  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(overlay).toHaveCount(0);
  expect(await artifactSnapshot(page)).toEqual(before);
  await align.click();
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    window.__S03Result = state.results.value[state.selectedResultIndex.value];
    state.results.value = [{ ...window.__S03Result }];
  });
  await expect(overlay).toHaveCount(0);
  await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    state.results.value = [];
  });
  await expect(overlay).toHaveCount(0);
  await page.evaluate(()=>window.__GBDRAW_APP__.cancelSimilarityAlignmentDialog());
  await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    state.results.value = [window.__S03Result];
    delete window.__S03Result;
  });
  await expect(page.locator('.gbdraw-preview-surface svg')).toBeVisible();
  await align.click();
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await dialog.getByRole('radio', { name: /Select .*bp, strand/ }).first().check();
  await dialog.getByRole('button', { name: 'Apply', exact: true }).click();
  await expect(dialog).toBeHidden({ timeout: 180000 });
  await expect(overlay).toHaveCount(0);
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
  await loadAmbiguousSession(page, testInfo);

  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) {
    await page.locator('.drawer-toggle').click();
  }
  const drawer = page.locator('.right-drawer');
  await drawer.getByRole('button', { name: 'Similarity groups' }).click();
  await drawer.locator('button').filter({
    has: page.locator('.font-mono', { hasText: /^og_1$/ })
  }).click();

  const align = drawer.getByRole('button', { name: 'Align…', exact: true });
  const drawerReview = drawer.getByRole('button', { name: 'Review alignment options…' });
  const referenceSelect = drawer.getByLabel('Exact reference record and feature');
  await expect(align).toBeDisabled();
  await expect(drawerReview).toBeDisabled();
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
  await expect(drawerReview).toBeEnabled();

  const beforeCancel = await artifactSnapshot(page);
  const diagramWidthBefore = (await page.locator('.gbdraw-preview-surface').boundingBox()).width;
  const immediate = await page.evaluate(async () => {
    const button = document.querySelector('button[title="Align to the selected exact feature"]');
    button.click();
    await window.Vue.nextTick();
    return {
      text: button.textContent.trim(),
      disabled: button.disabled,
      busy: button.getAttribute('aria-busy'),
      live: document.querySelector('[data-similarity-alignment-status]')?.textContent.trim()
    };
  });
  expect(immediate).toMatchObject({
    text: 'Resolving…', disabled: true, busy: 'true',
    live: 'Resolving similarity alignment…'
  });
  const dialog = page.getByRole('dialog', { name: 'Select alignment anchors' });
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await expect(dialog.getByRole('heading', { name: 'Select alignment anchors' })).toBeFocused();
  await expect(dialog).toHaveAttribute('aria-describedby', 'similarity-alignment-dialog-description');
  const referenceCard = dialog.locator('[data-similarity-alignment-reference]');
  await expect(referenceCard).toContainText('b0');
  await expect(referenceCard).toContainText('record_b');
  await expect(referenceCard).toContainText('121..419 bp');
  await expect(referenceCard).toContainText('Strand: +');
  await expect(dialog.locator('[data-alignment-record-key]')).toHaveCount(1);
  await expect(dialog.locator('[data-similarity-alignment-reason]')).toContainText(
    'The only representative candidate.'
  );
  await dialog.getByRole('radio', {name:/Select .*bp, strand/}).first().check();
  await expect(dialog.getByText('Recommended', { exact: true })).toHaveCount(1);
  await expect(dialog).not.toHaveAttribute('aria-modal', 'true');
  expect(await dialog.evaluate((element) => element.contains(document.activeElement))).toBe(true);
  expect((await page.locator('.gbdraw-preview-surface').boundingBox()).width).toBe(diagramWidthBefore);
  await expect(page.locator('[data-similarity-alignment-count]')).toContainText('0 need a choice');
  const paletteHeader = dialog.locator('header');
  const headerBox = await paletteHeader.boundingBox();
  await page.mouse.move(headerBox.x + 35, headerBox.y + 20);
  await page.mouse.down();
  await page.mouse.move(-200, -200, { steps: 5 });
  await page.mouse.up();
  let paletteBox = await dialog.boundingBox();
  expect(paletteBox.x).toBeGreaterThanOrEqual(11);
  expect(paletteBox.y).toBeGreaterThanOrEqual(11);
  await page.mouse.move(paletteBox.x + 35, paletteBox.y + 20);
  await page.mouse.down();
  await page.mouse.move(2000, 2000, { steps: 5 });
  await page.mouse.up();
  paletteBox = await dialog.boundingBox();
  expect(paletteBox.x + paletteBox.width).toBeLessThanOrEqual(1589);
  expect(paletteBox.y + paletteBox.height).toBeLessThanOrEqual(989);
  await page.setViewportSize({ width: 900, height: 650 });
  paletteBox = await dialog.boundingBox();
  if(await page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentCompact)) {
    await expect(page.locator('[data-similarity-alignment-dock]')).toContainText('Select alignment anchors');
    expect((await page.locator('.preview-viewport').boundingBox()).height).toBeGreaterThanOrEqual(200);
  } else {
    expect(paletteBox.x + paletteBox.width).toBeLessThanOrEqual(889);
    expect(paletteBox.y + paletteBox.height).toBeLessThanOrEqual(639);
  }
  await page.setViewportSize({ width: 1600, height: 1000 });
  const canvasBox = await page.locator('.gbdraw-preview-surface').locator('..').boundingBox();
  const beforeViewport = await page.evaluate(() => ({
    zoom: window.__GBDRAW_APP__.zoom,
    panX: window.__GBDRAW_APP__.canvasPan.x
  }));
  await page.mouse.move(canvasBox.x + 24, canvasBox.y + 100);
  await page.mouse.wheel(0, -100);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.zoom))
    .toBeGreaterThan(beforeViewport.zoom);
  await page.mouse.down();
  await page.mouse.move(canvasBox.x + 64, canvasBox.y + 100, { steps: 4 });
  await page.mouse.up();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.canvasPan.x))
    .toBeGreaterThan(beforeViewport.panX);
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.linearSourceRemovalDialog.sourceUid = app.linearSeqs[0].uid;
    app.linearSourceRemovalDialog.origin = 'card';
    app.linearSourceRemovalDialog.open = true;
  });
  const removalDialog = page.getByRole('dialog', { name: 'Clear or delete File?' });
  await expect(removalDialog).toBeVisible();
  await removalDialog.getByRole('button', { name: 'Cancel' }).focus();
  await page.keyboard.press('Escape');
  await expect(removalDialog).toBeHidden();
  await expect(dialog).toBeVisible();
  expect(await artifactSnapshot(page)).toEqual(beforeCancel);
  await expect(dialog).toContainText(/\d[\d,]*\.\.\d[\d,]* bp · Strand [+-]/);
  await expect(dialog).toContainText('Direct evidence: None');
  await expect(dialog).not.toContainText(/score/i);
  const apply = dialog.getByRole('button', { name: 'Apply', exact: true });
  await expect(apply).toBeEnabled();
  await expect(apply).toHaveAttribute('aria-describedby', 'similarity-alignment-apply-reason');
  await expect(page.locator('#similarity-alignment-apply-reason')).toContainText(
    'complete draft'
  );
  await page.screenshot({
    path: testInfo.outputPath('desktop-ambiguity.png'),
    fullPage: true
  });

  const firstCandidate = dialog.getByRole('radio', { name: /Select .*bp, strand/ }).first();
  const candidateAnchor = await page.evaluate(() => (
    window.__GBDRAW_APP__.similarityAlignmentDraft.rows[0].candidates[0].anchor
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
  const replacementCandidate = dialog.getByRole('radio', { name: /Select .*bp, strand/ }).nth(1);
  await replacementCandidate.check();
  await expect(replacementCandidate).toBeChecked();
  await expect(dialog.getByText('Recommended', { exact: true })).toHaveCount(0);
  await expect(dialog.locator('[data-similarity-alignment-reason]')).toContainText(
    'The only representative candidate.'
  );
  await firstCandidate.check();
  await expect(dialog.getByText('Recommended', { exact: true })).toHaveCount(1);

  await dialog.getByRole('radio', { name: /Skip / }).focus();
  await page.keyboard.press('Space');
  await expect(dialog).toContainText('Unchanged');
  await expect(dialog.getByText('Recommended', { exact: true })).toHaveCount(0);
  await expect(apply).toBeEnabled();
  await apply.focus();
  await page.keyboard.press('Tab');
  expect(await dialog.evaluate((element) => element.contains(document.activeElement))).toBe(false);
  await page.keyboard.press('Escape');
  await expect(dialog).toBeHidden();
  if(!await page.evaluate(()=>window.__GBDRAW_APP__.showRightDrawer)) {
    await expect(page.locator('.drawer-toggle')).toBeFocused();
    expect(await page.evaluate(()=>window.__GBDRAW_APP__.rightDrawerTab)).toBe('orthogroups');
    await page.locator('.drawer-toggle').click();
  } else await expect(align).toBeFocused();
  expect(await artifactSnapshot(page)).toEqual(beforeCancel);

  await align.click();
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await dialog.getByRole('radio', { name: /Select .*bp, strand/ }).first().check();
  const staleOutcome = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    state.selectedOrthogroupId.value = 'another-group';
    const outcome = await window.__GBDRAW_APP__.applySimilarityAlignmentDraft();
    state.selectedOrthogroupId.value = 'og_1';
    return outcome;
  });
  expect(staleOutcome.status).toBe('stale');
  await expect(dialog).toBeHidden();
  expect(await artifactSnapshot(page)).toEqual(beforeCancel);

  await page.setViewportSize({ width: 390, height: 740 });
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
  const popupAlign = popup.getByRole('button', { name: 'Align…', exact: true });
  await expect(popupAlign).toHaveAttribute('title', /Align to this exact feature/);
  const popupReview = popup.getByRole('button', { name: 'Review alignment options…' });
  await expect(popupReview).toBeVisible();
  await popupReview.scrollIntoViewIfNeeded();
  const popupImmediate = await page.evaluate(async () => {
    const review = document.querySelector('[data-similarity-alignment-popup-review]');
    const alignButton = document.querySelector(
      '.feature-popup button[title="Align to this exact feature"]'
    );
    const drawerButton = document.querySelector(
      'button[title="Align to the selected exact feature"]'
    );
    review.click();
    await window.Vue.nextTick();
    return {
      reviewText: review.textContent.trim(),
      reviewDisabled: review.disabled,
      reviewBusy: review.getAttribute('aria-busy'),
      alignDisabled: alignButton.disabled,
      drawerDisabled: drawerButton.disabled,
      drawerBusy: drawerButton.getAttribute('aria-busy')
    };
  });
  expect(popupImmediate).toEqual({
    reviewText: 'Review alignment options…', reviewDisabled: true, reviewBusy: 'true',
    alignDisabled: true, drawerDisabled: true, drawerBusy: 'true'
  });
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await page.setViewportSize({ width: 390, height: 500 });
  const paletteBody = dialog.locator('.custom-scrollbar').first();
  expect(await paletteBody.evaluate((element) => element.scrollHeight > element.clientHeight)).toBe(true);
  expect(await paletteBody.evaluate((element) => element.scrollWidth <= element.clientWidth)).toBe(true);
  expect(await dialog.locator('footer').evaluate(
    (element) => parseFloat(getComputedStyle(element).paddingBottom)
  )).toBeLessThan(24);
  await expect(dialog.locator('[data-similarity-alignment-reference]')).toBeVisible();
  await expect(dialog.locator('[data-similarity-alignment-reason]')).toBeVisible();
  await expect(dialog.getByRole('checkbox', { name: /Match reference direction for/ })).toHaveCount(0);
  await paletteBody.evaluate((element) => { element.scrollTop = element.scrollHeight; });
  expect(await paletteBody.evaluate((element) => element.scrollTop)).toBeGreaterThan(0);
  await page.setViewportSize({ width: 390, height: 740 });
  await dialog.getByRole('radio', { name: /Select .*bp, strand/ }).first().check();
  await expect(apply).toBeEnabled();
  await apply.scrollIntoViewIfNeeded();
  const applyBox = await apply.boundingBox();
  expect(applyBox).not.toBeNull();
  expect(applyBox.x).toBeGreaterThanOrEqual(0);
  expect(applyBox.x + applyBox.width).toBeLessThanOrEqual(390);
  expect(applyBox.y + applyBox.height).toBeLessThanOrEqual(740);
  const narrowPaletteBox = await dialog.boundingBox();
  expect(narrowPaletteBox.x).toBeGreaterThanOrEqual(0);
  expect(narrowPaletteBox.x + narrowPaletteBox.width).toBeLessThanOrEqual(390);
  await page.screenshot({
    path: testInfo.outputPath('narrow-ambiguity.png'),
    fullPage: false
  });
  await dialog.locator('[data-alignment-record-key] input[type="radio"]:checked').press('Escape');
  await expect(dialog).toBeHidden();
  await expect(popupReview).toBeFocused();
  await popup.getByRole('button', { name: 'Close feature popup' }).click();
  expect(await page.evaluate(()=>window.__GBDRAW_APP__.showRightDrawer)).toBe(false);
  await page.locator('.drawer-toggle').click();
  expect(await page.evaluate(()=>window.__GBDRAW_APP__.rightDrawerTab)).toBe('orthogroups');

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
  await dialog.getByRole('radio', { name: /Select .*bp, strand/ }).first().check();
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
  expect(Object.keys(savedState.plan.records.find(({ status }) => status === 'aligned')).sort())
    .toEqual(['anchor', 'rationale', 'recordKey', 'status']);
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
  const regenerate = await freshPage.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const historyBefore = window.__GBDRAW_HISTORY__.getUndoCount();
    const result = await app.runAnalysis();
    return {
      result,
      notice: app.similarityAlignmentNotice,
      historyBefore,
      historyAfter: window.__GBDRAW_HISTORY__.getUndoCount()
    };
  });
  expect(regenerate.result).toEqual({ status: 'ok' });
  expect(regenerate.notice).not.toContain('needs repair');
  expect(regenerate.historyAfter).toBe(regenerate.historyBefore + 1);
  expect(await freshPage.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return {
      plan: JSON.parse(JSON.stringify(state.similarityAlignmentPlan.value)),
      translations: JSON.parse(JSON.stringify(state.linearRecordTranslations.value)),
      results: state.results.value.map(({ name, content }) => ({ name, content }))
    };
  })).toEqual({
    plan: savedState.plan,
    translations: savedState.translations,
    results: savedState.results
  });
  const historyBeforeReset = await freshPage.evaluate(
    () => window.__GBDRAW_HISTORY__.getUndoCount()
  );
  await freshPage.locator('[data-similarity-alignment-reset]').click();
  await freshPage.getByRole('dialog',{name:'Reset alignment',exact:true}).getByRole('button',{name:'Reset',exact:true}).click();
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

const openAlignmentFromDrawer = async (page) => {
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
  await drawer.getByLabel('Exact reference record and feature').selectOption(reference);
  await drawer.getByRole('button', { name: 'Align…', exact: true }).click();
  const dialog = page.getByRole('dialog', { name: 'Select alignment anchors' });
  await expect(dialog).toBeVisible({ timeout: 180000 });
  return dialog;
};

test('exclusive direction modes and Custom follow keyboard choices without Worker jobs', async ({ page }, testInfo) => {
  test.setTimeout(180000);
  await loadAmbiguousSession(page, testInfo);
  const dialog = await openAlignmentFromDrawer(page);
  await page.evaluate(() => {
    window.__s03Jobs = 0;
    const post = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function(message, transfer) {
      window.__s03Jobs++;
      return post.call(this, message, transfer);
    };
  });
  const before = await artifactSnapshot(page);
  const modes = dialog.getByRole('group', { name: 'Alignment direction', exact: true });
  await expect(modes.getByRole('radio')).toHaveCount(4);
  await expect(modes.getByRole('radio', { name: 'Keep current directions', exact: true })).toBeChecked();
  await expect(dialog.getByRole('combobox')).toHaveCount(0);
  const first = dialog.getByRole('radio', { name: /Select .*bp, strand/ }).first();
  const skip = dialog.getByRole('radio', { name: /Skip / });
  await first.focus(); await page.keyboard.press('Space');
  await expect(first).toBeChecked();
  await page.keyboard.press('ArrowDown');
  await expect(dialog.getByRole('radio', { name: /Select .*bp, strand/ }).nth(1)).toBeChecked();
  await page.keyboard.press('ArrowDown'); await expect(skip).toBeChecked();
  await expect(dialog.locator('#alignment-direction-scope')).toContainText('1 selected anchors');
  await modes.getByRole('radio', { name: 'Keep current directions', exact: true }).focus();
  await page.keyboard.press('ArrowRight');
  await expect(modes.getByRole('radio', { name: 'All selected arrows right', exact: true })).toBeChecked();
  await page.keyboard.press('ArrowRight');
  await expect(modes.getByRole('radio', { name: 'All selected arrows left', exact: true })).toBeChecked();
  await page.keyboard.press('ArrowRight');
  await expect(modes.getByRole('radio', { name: 'Custom', exact: true })).toBeChecked();
  await expect(modes.locator('input:checked')).toHaveCount(1);
  const target = dialog.getByLabel('Direction for record_a', { exact: true });
  await expect(target).toBeDisabled();
  await expect(dialog.locator('#alignment-direction-record_a')).toContainText('Skipped by user');
  await first.focus(); await page.keyboard.press('Space');
  await expect(target).toBeEnabled();
  await expect(dialog.locator('#alignment-direction-scope')).toContainText('2 selected anchors');
  const reference = dialog.getByLabel('Direction for record_b', { exact: true });
  await reference.selectOption('left'); await target.selectOption('right');
  await expect(reference).toHaveAccessibleDescription(/exact reference.*Current right-facing; after Align left-facing/);
  await expect(dialog.locator('#alignment-direction-scope')).toHaveAttribute('aria-live', 'polite');
  await page.setViewportSize({ width: 390, height: 740 });
  await reference.scrollIntoViewIfNeeded();
  await expect(reference).toBeVisible(); await expect(target).toBeEnabled();
  await reference.focus(); await page.keyboard.press('Home'); await page.keyboard.press('ArrowDown');
  await expect(reference).toHaveValue('right');
  await target.selectOption('left'); await expect(target).toHaveValue('left');
  expect(await dialog.locator('.custom-scrollbar').evaluate(e => e.scrollWidth <= e.clientWidth)).toBe(true);
  await dialog.screenshot({ path: testInfo.outputPath('custom-390.png') });
  await modes.getByRole('radio', { name: 'Keep current directions', exact: true }).check();
  await expect(dialog.getByRole('combobox')).toHaveCount(0);
  expect(await page.evaluate(() => window.__s03Jobs)).toBe(0);
  expect(await artifactSnapshot(page)).toEqual(before);
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).focus();
  await page.keyboard.press('Enter'); await expect(dialog).toBeHidden();
  await expect(page.locator('.drawer-toggle')).toBeFocused();
});

test('Apply error retains draft, focuses retry guidance, and accepts correction', async ({ page }, testInfo) => {
  test.setTimeout(180000);
  await loadAmbiguousSession(page, testInfo);
  const dialog = await openAlignmentFromDrawer(page);
  await dialog.getByRole('radio',{name:/Select .*bp, strand/}).first().check();
  const before = await artifactSnapshot(page);
  const priorLayout = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const prior = app.form.linear_track_layout;
    const sanitize=window.DOMPurify.sanitize;window.DOMPurify.sanitize=(...args)=>{if(String(args[0]).includes('<svg')){window.DOMPurify.sanitize=sanitize;throw new Error('Forced candidate post-processing failure.');}return sanitize(...args);};
    return prior;
  });
  await dialog.getByRole('button', { name: 'Apply', exact: true }).click();
  const error = dialog.locator('[data-similarity-alignment-error]');
  await expect(error).toBeVisible({ timeout: 180000 });
  await expect(error).toHaveText('Forced candidate post-processing failure.');
  await expect(page.getByRole('alert', { name: 'Generation Error' })
    .locator('.text-sm.font-semibold')).toHaveText('Forced candidate post-processing failure.');
  await expect(error).toBeFocused();
  await expect(dialog.getByRole('button', { name: 'Apply', exact: true })).toBeEnabled();
  await page.setViewportSize({ width: 390, height: 500 });
  await expect.poll(async () => {
    const errorBox = await error.boundingBox();
    return errorBox.x >= 0 && errorBox.x + errorBox.width <= 390;
  }).toBe(true);
  await expect(dialog.getByRole('button', { name: 'Cancel', exact: true })).toBeVisible();
  await expect(dialog.getByRole('button', { name: 'Apply', exact: true })).toBeVisible();
  expect(await dialog.locator('.custom-scrollbar').evaluate(
    (element) => element.scrollWidth <= element.clientWidth
  )).toBe(true);
  const afterFailure = await artifactSnapshot(page);
  expect(afterFailure.results).toEqual(before.results);
  expect(afterFailure.history.slice(0, 2)).toEqual(before.history.slice(0, 2));
  await page.evaluate((layout) => {
    window.__GBDRAW_APP__.form.linear_track_layout = layout;
  }, priorLayout);
  await dialog.getByRole('radio', { name: /Skip / }).check();
  await expect(error).toBeHidden();
  await dialog.getByRole('button', { name: 'Apply', exact: true }).click();
  await expect(dialog).toBeHidden({ timeout: 180000 });
  await expect(page.locator('[data-similarity-alignment-summary]')).toContainText('explicitly skipped');
});

test('one usable member applies directly through one Worker resolve and one History action', async ({ page }, testInfo) => {
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
  const before = await artifactSnapshot(page);
  const drawerReview = drawer.getByRole('button', { name: 'Review alignment options…' });
  const normalAlign = drawer.getByRole('button', { name: 'Align…', exact: true });
  await expect(drawerReview).toBeEnabled();
  await expect(normalAlign).toBeEnabled();
  await drawerReview.focus();
  const immediateReview = await page.evaluate(async () => {
    const review = document.querySelector('[data-similarity-alignment-drawer-review]');
    const align = document.querySelector('button[title="Align to the selected exact feature"]');
    review.click();
    await window.Vue.nextTick();
    return {
      reviewBusy: review.getAttribute('aria-busy'), reviewDisabled: review.disabled,
      alignBusy: align.getAttribute('aria-busy'), alignDisabled: align.disabled,
      live: document.querySelector('[data-similarity-alignment-status]')?.textContent.trim()
    };
  });
  expect(immediateReview).toEqual({
    reviewBusy: 'true', reviewDisabled: true, alignBusy: 'true', alignDisabled: true,
    live: 'Resolving similarity alignment…'
  });
  const resolvedReview = page.getByRole('dialog', { name: 'Select alignment anchors' });
  await expect(resolvedReview).toBeVisible({ timeout: 180000 });
  await expect(resolvedReview.getByRole('heading', { name: 'Select alignment anchors' })).toBeFocused();
  await expect(resolvedReview.locator('[data-similarity-alignment-reference]')).toContainText('b0');
  await expect(resolvedReview.locator('[data-alignment-record-key]')).toHaveCount(1);
  await expect(resolvedReview.getByRole('button', { name: 'Apply', exact: true })).toBeEnabled();
  await expect(resolvedReview.getByRole('checkbox', { name: /Match reference direction for/ }))
    .toHaveCount(0);
  await resolvedReview.getByRole('radio', { name: /Skip / }).check();
  await expect(resolvedReview).toContainText('Unchanged');
  await resolvedReview.getByRole('radio', { name: /Select .*bp, strand/ }).first().check();
  await resolvedReview.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(drawerReview).toBeFocused();
  expect(await artifactSnapshot(page)).toEqual(before);

  await page.setViewportSize({ width: 390, height: 740 });
  const b0FeatureId = await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures
    .find(({ protein_id: proteinId }) => proteinId === 'b0').svg_id);
  await page.locator(`[data-gbdraw-feature-id="${b0FeatureId}"]`).first()
    .dispatchEvent('click', { clientX: -100, clientY: -100 });
  const popup = page.locator('.feature-popup[role="dialog"]');
  await expect(popup).toBeVisible();
  const popupReview = popup.getByRole('button', { name: 'Review alignment options…' });
  await popupReview.scrollIntoViewIfNeeded();
  await popupReview.focus();
  await popupReview.press('Enter');
  await expect(resolvedReview).toBeVisible({ timeout: 180000 });
  await expect(resolvedReview.locator('[data-similarity-alignment-reference]')).toContainText('b0');
  const narrowBox = await resolvedReview.boundingBox();
  expect(narrowBox.x).toBeGreaterThanOrEqual(0);
  expect(narrowBox.x + narrowBox.width).toBeLessThanOrEqual(390);
  await page.keyboard.press('Escape');
  await expect(resolvedReview).toBeHidden();
  await expect(popupReview).toBeFocused();
  await popup.getByRole('button', { name: 'Close feature popup' }).click();
  expect(await artifactSnapshot(page)).toEqual(before);
  await page.setViewportSize({ width: 1600, height: 1000 });
  if(!await page.evaluate(()=>window.__GBDRAW_APP__.showRightDrawer))await page.locator('.drawer-toggle').click();

  await page.evaluate(() => {
    window.__alignmentWorkerResolves = 0;
    window.__alignmentPaletteMounts = 0;
    const postMessage = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function (message, transfer) {
      if (message?.type === 'helper' && message.operation === 'resolveSimilarityAlignment') {
        window.__alignmentWorkerResolves += 1;
      }
      return postMessage.call(this, message, transfer);
    };
    new MutationObserver((changes) => {
      for (const change of changes) {
        for (const node of change.addedNodes) {
          if (node.nodeType === Node.ELEMENT_NODE && (
            node.matches('[data-similarity-alignment-dialog]')
            || node.querySelector('[data-similarity-alignment-dialog]')
          )) window.__alignmentPaletteMounts += 1;
        }
      }
    }).observe(document.body, { childList: true, subtree: true });
  });
  await normalAlign.click();
  await expect(page.getByRole('dialog', { name: 'Select alignment anchors' })).toHaveCount(0);
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return state.similarityAlignmentPlan.value?.schema || null;
  }), { timeout: 120000 }).toBe(2);
  await expect(page.locator('[data-similarity-alignment-summary]')).toContainText('1 aligned');
  const after = await artifactSnapshot(page);
  expect(after.results).not.toEqual(before.results);
  expect(after.history[0]).toBe(before.history[0] + 1);
  expect(after.history[1]).toBe(0);
  expect(await page.evaluate(() => ({
    resolves: window.__alignmentWorkerResolves,
    paletteMounts: window.__alignmentPaletteMounts
  }))).toEqual({ resolves: 1, paletteMounts: 0 });

  await importSession(page, session, 'one-candidate-retry.gbdraw-session.json');
  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) {
    await page.locator('.drawer-toggle').click();
  }
  const similarityTab = drawer.getByRole('button', { name: 'Similarity groups' });
  await similarityTab.click();
  await drawer.locator('button').filter({
    has: page.locator('.font-mono', { hasText: /^og_1$/ })
  }).click();
  await drawer.getByLabel('Exact reference record and feature').selectOption(reference);
  const committed = await artifactSnapshot(page);
  const priorLayout = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const prior = app.form.linear_track_layout;
    const sanitize=window.DOMPurify.sanitize;window.DOMPurify.sanitize=(...args)=>{if(String(args[0]).includes('<svg')){window.DOMPurify.sanitize=sanitize;throw new Error('Forced candidate post-processing failure.');}return sanitize(...args);};
    return prior;
  });
  await normalAlign.click();
  await expect(resolvedReview).toBeVisible({ timeout: 180000 });
  const retryError = resolvedReview.locator('[data-similarity-alignment-error]');
  await expect(retryError).toHaveText('Forced candidate post-processing failure.');
  await expect(page.getByRole('alert', { name: 'Generation Error' })
    .locator('.text-sm.font-semibold')).toHaveText('Forced candidate post-processing failure.');
  await expect(retryError).toBeFocused();
  await expect(resolvedReview.locator('[data-alignment-record-key]')).toHaveCount(1);
  const failed = await artifactSnapshot(page);
  expect(failed.results).toEqual(committed.results);
  expect(failed.history.slice(0, 2)).toEqual(committed.history.slice(0, 2));
  await page.evaluate((layout) => {
    window.__GBDRAW_APP__.form.linear_track_layout = layout;
  }, priorLayout);
  await resolvedReview.getByRole('radio', { name: /Skip / }).check();
  await resolvedReview.getByRole('button', { name: 'Apply', exact: true }).click();
  await expect(resolvedReview).toBeHidden({ timeout: 180000 });
  await expect(page.locator('[data-similarity-alignment-summary]')).toContainText('explicitly skipped');

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
    malformed.records[0].anchor.stableFeatureSvgId = 'missing-feature';
    malformed.reference.biologicalFeatureId = 'missing-feature';
    malformed.reference.stableFeatureSvgId = 'missing-feature';
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

test('Gallery explicit directions own receipts, Reset scopes, fresh Load and ribbon geometry', async ({ page, browser }, testInfo) => {
  test.setTimeout(600000);
  await captureUnhandledRejections(page);
  const pageErrors = [];
  page.on('pageerror', (error) => { pageErrors.push(String(error?.message || error)); console.log('Gallery page error:', String(error)); });
  page.on('console', (message) => { if (message.type() === 'error') console.log('Gallery console error:', message.text()); });
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(
    'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'
  ), 'BGC0000708-BGC0000713.gbdraw-session.json');

  const recordState = (targetPage = page) => targetPage.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return {
      receipt: JSON.parse(JSON.stringify(state.similarityAlignmentResetReceipt.value)),
      plan: JSON.parse(JSON.stringify(state.similarityAlignmentPlan.value)),
      translations: JSON.parse(JSON.stringify(state.linearRecordTranslations.value)),
      orientations: state.linearSeqs.map(({ uid, region_reverse }) => ({
        recordKey: uid, reverseComplement: Boolean(region_reverse)
      })),
      results: state.results.value.map(({ name, content }) => ({ name, content })),
      positions: Array.from(document.querySelectorAll(
        '.gbdraw-preview-surface svg [data-gbdraw-composition-role="primary"][data-gbdraw-record-id]'
      )).filter((element) => !element.hasAttribute('data-gbdraw-definition-part'))
        .map((element) => ({ recordId: element.getAttribute('data-gbdraw-record-id'),
          transform: element.getAttribute('transform') })),
      historyCount: window.__GBDRAW_HISTORY__.getUndoCount()
    };
  });
  const review = async () => {
    const featureId = await page.evaluate(() => {
      const feature = window.__GBDRAW_APP__.extractedFeatures.find(
        (item) => item.protein_id === 'CAG38712.1' || item.sourceProteinId === 'CAG38712.1'
      );
      if (!feature) throw new Error('Gallery reference CAG38712.1 is unavailable.');
      return feature.svg_id;
    });
    await page.locator(`[data-gbdraw-feature-id="${featureId}"]`).first()
      .dispatchEvent('click', { clientX: -100, clientY: -100 });
    await page.locator('.feature-popup[role="dialog"]')
      .getByRole('button', { name: 'Review alignment options…' }).click();
    const dialog = page.getByRole('dialog', { name: 'Select alignment anchors' });
    await expect(dialog).toBeVisible({ timeout: 180000 });
    const unresolved=await page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentDraft.rows.filter(row=>!row.choice).map(row=>row.recordKey));
    for(const key of unresolved)await dialog.locator(`[data-alignment-record-key="${key}"]`).getByRole('radio',{name:/Select .*bp, strand/}).first().check();
    return dialog;
  };
  const geometry = () => page.evaluate(async () => {
    const { getFeatureFillElements } = await import('./js/app/feature-dom.js');
    const { state } = await import('./js/state.js');
    const svg = document.querySelector('.gbdraw-preview-surface svg');
    const plan = state.similarityAlignmentPlan.value;
    const center = (anchor) => {
      const feature = state.extractedFeatures.value.find((item) => (
        item.recordKey === anchor.recordKey && item.biologicalFeatureId === anchor.biologicalFeatureId
      ));
      const bounds = getFeatureFillElements(svg, feature.svg_id).map((element) => element.getBoundingClientRect());
      return (Math.min(...bounds.map((b) => b.left)) + Math.max(...bounds.map((b) => b.right))) / 2;
    };
    const reference = center(plan.reference);
    const offsets = plan.records.filter(({ status }) => status === 'aligned').map(({ recordKey, anchor }) => ({
      recordKey, offset: center(anchor) - reference
    }));
    const ribbons = Array.from(svg.querySelectorAll('[data-gbdraw-pairwise-match-id]'));
    return { offsets, ribbons: ribbons.length,
      targetRibbons: ribbons.filter((element) => element.getAttribute('data-query-record-index') === '1'
        || element.getAttribute('data-subject-record-index') === '1').length };
  });
  const applyReviewed=async()=>{
    const historyBefore=(await recordState()).historyCount;
    await expect(apply).toBeEnabled();
    await apply.click();await page.waitForFunction(()=>!window.__GBDRAW_APP__.similarityAlignmentBusy,null,{timeout:180000});
    const message=await page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentError?.message||'');
    console.log('S02 Gallery: final batch outcome',message,await page.evaluate(()=>{const e=window.__GBDRAW_APP__.errorLog;return {message:e?.message,stack:e?.stack,summary:e?.summary,details:e?.details};}));
    if(message.includes('Validated directions or reference placement changed')){
      expect((await recordState()).historyCount).toBe(historyBefore);
      await apply.click();await page.waitForFunction(()=>!window.__GBDRAW_APP__.similarityAlignmentBusy,null,{timeout:180000});
    }
  };
  const before = await recordState();
  console.log('S02 Gallery: imported baseline');
  // Default Align renders comparisons and preserves every record direction.
  const defaultReference = before.plan.records.find(({ recordKey }) => recordKey === 'record-2').anchor;
  const referenceId = await page.evaluate((anchor) => window.__GBDRAW_APP__.extractedFeatures
    .find((item) => item.recordKey === anchor.recordKey
      && item.biologicalFeatureId === anchor.biologicalFeatureId).svg_id, defaultReference);
  console.log('S02 Gallery: default reference',referenceId);
  await page.locator(`[data-gbdraw-feature-id="${referenceId}"]`).first()
    .dispatchEvent('click', { clientX: -100, clientY: -100 });
  console.log('S02 Gallery: invoking default Align');
  await page.locator('.feature-popup[role="dialog"]')
    .getByRole('button', { name: 'Align…', exact: true }).click();
  console.log('S02 Gallery: default Align invoked');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentBusy),
    { timeout: 180000 }).toBe(false);
  console.log('S02 Gallery: default Align finished',await page.evaluate(()=>({status:window.__GBDRAW_APP__.similarityAlignmentStatus,error:window.__GBDRAW_APP__.similarityAlignmentError?.message})));
  await expect(page.locator('[data-similarity-alignment-summary]')).toContainText('0 reversed');
  expect((await recordState()).orientations).toEqual(before.orientations);
  await expect(page.getByRole('dialog', { name: 'Select alignment anchors' })).toHaveCount(0);
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(recordState).toEqual(before);
  console.log('S02 Gallery: Undo restored baseline');
  const dialog = await review();
  console.log('S02 Gallery: initial review');
  await expect(dialog.locator('[data-similarity-alignment-reference]')).toContainText('Streptomyces lividus CBS 844.73');
  const targetRow = dialog.locator('[data-alignment-record-key="record-2"]');
  await expect(targetRow).toContainText('Streptomyces fradiae');
  await expect(dialog.getByRole('checkbox',{name:'Match reference direction'})).toHaveCount(0);
  const direction = dialog.getByRole('radio', { name: 'All selected arrows left', exact: true });
  await expect(direction).not.toBeChecked();
  await direction.focus();
  await page.keyboard.press('Space');
  await expect(direction).toBeChecked();
  await expect(dialog.locator('[data-similarity-alignment-direction]').filter({hasText:'record-2:'})).toContainText('→ → ←');
  await dialog.screenshot({ path: testInfo.outputPath('direction-review-desktop.png') });
  await page.setViewportSize({ width: 390, height: 740 });
  await direction.scrollIntoViewIfNeeded();
  await expect(direction).toBeVisible();
  const bounds = await dialog.boundingBox();
  expect(bounds.x).toBeGreaterThanOrEqual(0);
  expect(bounds.x + bounds.width).toBeLessThanOrEqual(390);
  const apply = dialog.getByRole('button', { name: 'Apply', exact: true });
  await expect(apply).toBeVisible();
  expect(await dialog.evaluate((element) => element.scrollWidth <= element.clientWidth)).toBe(true);
  for(const height of [740,844]) {
    await page.setViewportSize({width:390,height});
    await expect.poll(()=>page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentCompact)).toBe(true);
    await expect(page.locator('.drawer-toggle')).toBeDisabled();
    await expect(page.locator('.drawer-toggle')).toHaveAttribute('title',/Finish or cancel alignment review/);
    expect(await page.evaluate(()=>window.__GBDRAW_APP__.showRightDrawer)).toBe(false);
    await page.locator('.preview-canvas').evaluate(element => element.scrollIntoView({ block: 'start' }));
    const canvas=await page.locator('.preview-canvas').boundingBox();
    const palette=await dialog.boundingBox();
    expect(canvas.height).toBeGreaterThanOrEqual(200);
    expect(Math.abs(canvas.width-palette.width)).toBeLessThanOrEqual(2);
    const visible = await page.evaluate(box => {
      const header = document.querySelector('.app-header').getBoundingClientRect();
      const generate = document.querySelector('.generate-bar').getBoundingClientRect();
      const overlaps = generate.left < box.x + box.width && generate.right > box.x;
      return Math.min(box.y + box.height, window.innerHeight, overlaps ? generate.top : window.innerHeight)
        - Math.max(box.y, 0, header.bottom);
    }, canvas);
    expect(visible).toBeGreaterThanOrEqual(200);
    await page.screenshot({path:testInfo.outputPath(`canvas-hit-${height}.png`)});
    const hit=await page.evaluate(box=>{
      const element=document.elementFromPoint(box.x+box.width/2,Math.max(0,box.y)+100);
      return {inside:Boolean(element?.closest('.preview-canvas')),tag:element?.tagName,classes:element?.className,
        headerBottom:document.querySelector('.app-header').getBoundingClientRect().bottom,canvas:box};
    },canvas);
    console.log('S02 compact canvas hit',JSON.stringify(hit));
    expect(hit.inside).toBe(true);
    for (const control of [apply, dialog.getByRole('button', { name: 'Cancel', exact: true })]) {
      expect(await control.evaluate(element => {
        const box = element.getBoundingClientRect();
        const hit = document.elementFromPoint(box.x + box.width / 2, box.y + box.height / 2);
        return hit === element || element.contains(hit);
      })).toBe(true);
    }
    console.log('Visible compact canvas height:', height, visible);
    await page.screenshot({ path: testInfo.outputPath(`direction-review-390-${height}.png`) });
  }
  await page.setViewportSize({ width: 1600, height: 1000 });
  // A failed direction Apply retains the checked draft without committing reversals.
  await page.evaluate(() => {
    const sanitize = window.DOMPurify.sanitize;
    window.DOMPurify.sanitize = (...args) => {
      if (String(args[0]).includes('<svg')) {
        window.DOMPurify.sanitize = sanitize;
        throw new Error('Forced direction candidate post-processing failure.');
      }
      return sanitize(...args);
    };
  });
  console.log('S02 Gallery: Apply requested');
  await applyReviewed();
  await expect(dialog.locator('[data-similarity-alignment-error]'))
    .toHaveText('Forced direction candidate post-processing failure.');
  await expect(page.getByRole('alert', { name: 'Generation Error' })
    .locator('.text-sm.font-semibold'))
    .toHaveText('Forced direction candidate post-processing failure.');
  await expect(direction).toBeChecked();
  expect(await recordState()).toEqual(before);
  console.log('S02 Gallery: Apply requested');
  await applyReviewed();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
  const applyStatus = await page.evaluate(() => ({ status: window.__GBDRAW_APP__.similarityAlignmentStatus,
    error: window.__GBDRAW_APP__.similarityAlignmentError?.message, globalError: window.__GBDRAW_APP__.errorLog }));
  expect(applyStatus, JSON.stringify(applyStatus)).toMatchObject({ status: 'idle' });
  await expect(dialog).toBeHidden({ timeout: 180000 });
  const applied = await recordState();
  console.log('S02 Gallery: applied');
  expect(applied.orientations.find(({ recordKey }) => recordKey === 'record-2').reverseComplement).toBe(true);
  expect(applied.historyCount).toBe(before.historyCount + 1);
  const reversedCount = applied.orientations.filter((entry, index) => (
    entry.reverseComplement !== before.orientations[index].reverseComplement
  )).length;
  expect(reversedCount).toBe(4);
  await expect(page.locator('[data-similarity-alignment-summary]')).toContainText(`${reversedCount} reversed`);
  expect(applied.orientations[0]).toEqual(before.orientations[0]);
  expect(applied.translations.map(({ y }) => y)).toEqual(before.translations.map(({ y }) => y));
  const revLabels = () => page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const previous = { open: app.showRightDrawer, tab: app.rightDrawerTab };
    app.showRightDrawer = true;
    app.rightDrawerTab = 'orthogroups';
    await window.Vue.nextTick();
    const labels = [...document.querySelectorAll(
      '[data-similarity-alignment-plan-inspector] [aria-label$=" reversed relative to source"]'
    )].map((element) => element.getAttribute('aria-label'));
    app.showRightDrawer = previous.open;
    app.rightDrawerTab = previous.tab;
    await window.Vue.nextTick();
    return labels;
  });
  expect(await revLabels()).toEqual(applied.orientations.filter(({ reverseComplement }) => reverseComplement)
    .map(({ recordKey }) => `${recordKey} reversed relative to source`));
  const measured = await geometry();
  expect(measured.ribbons).toBeGreaterThan(0);
  expect(measured.targetRibbons).toBeGreaterThan(0);
  expect(Math.max(...measured.offsets.map(({ offset }) => Math.abs(offset)))).toBeLessThanOrEqual(0.5);
  await testInfo.attach('direction-anchor-offsets', { body: JSON.stringify(measured, null, 2), contentType: 'application/json' });
  console.log('Direction Gallery anchor offsets:', JSON.stringify(measured));
  await page.screenshot({ path: testInfo.outputPath('direction-applied.png'), fullPage: true });

  await review();
  await expect(direction).not.toBeChecked();
  await expect(direction).toBeEnabled();
  await page.keyboard.press('Escape');
  await expect(dialog).toBeHidden();
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(recordState).toEqual(before);
  console.log('S02 Gallery: Undo restored baseline');
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(recordState).toEqual(applied);

  console.log('S02 Gallery: positions Reset');
  const reset = await page.evaluate(() => window.__GBDRAW_APP__.resetSimilarityAlignment());
  expect(reset.status).toBe('ok');
  const resetState = await recordState();
  expect(resetState.plan).toBeNull();
  expect(resetState.orientations).toEqual(applied.orientations);
  expect(resetState.translations).toEqual(applied.translations);
  expect(before.positions).toHaveLength(5);
  expect(resetState.positions).toEqual(before.positions);
  await expect(page.locator('[data-similarity-alignment-notice]')).toContainText(
    'Alignment reset: record positions restored; record directions unchanged.'
  );
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(recordState).toEqual(applied);

  await page.evaluate(() => window.__GBDRAW_APP__.closeFeaturePopup());
  await page.getByRole('button', { name: 'Record options for sequence 2', exact: true }).click();
  const reverse = page.getByRole('checkbox', { name: 'Reverse complement for sequence 2', exact: true });
  await expect(reverse).toBeChecked();
  await reverse.uncheck();
  expect((await recordState()).plan).toEqual(applied.plan);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentNotice))
    .not.toContain('Alignment cleared');
  console.log('S02 Gallery: manual Reverse generation');
  const generated = await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis());
  expect(generated.status).toBe('ok');
  const manual = await recordState();
  expect(manual.plan).toEqual(applied.plan);
  expect(manual.receipt).toEqual(applied.receipt);
  await page.evaluate(()=>window.__GBDRAW_APP__.openSimilarityAlignmentReset());
  const resetReview=page.getByRole('dialog',{name:'Reset alignment',exact:true});
  await expect(resetReview).toContainText('This record was manually changed after Align');
  expect(manual.receipt.directions.length).toBeGreaterThan(1);
  const restorationTargets = resetReview.getByRole('list', { name: 'Alignment direction restoration targets' });
  await expect(restorationTargets.getByRole('listitem')).toHaveCount(manual.receipt.directions.length);
  await page.setViewportSize({ width: 390, height: 740 });
  await restorationTargets.getByRole('listitem').last().scrollIntoViewIfNeeded();
  await expect(restorationTargets.getByRole('listitem').last()).toBeVisible();
  expect(await resetReview.evaluate(element => element.scrollWidth <= element.clientWidth)).toBe(true);
  const resetFooter = await resetReview.getByRole('button', { name: 'Reset', exact: true }).boundingBox();
  expect(resetFooter.y + resetFooter.height).toBeLessThanOrEqual(740);
  await resetReview.screenshot({ path: testInfo.outputPath('reset-targets-390.png') });
  await page.setViewportSize({ width: 1600, height: 1000 });
  await resetReview.getByRole('button',{name:'Cancel',exact:true}).click();
  expect(await recordState()).toEqual(manual);
  expect(manual.orientations.find(({ recordKey }) => recordKey === 'record-2').reverseComplement).toBe(false);
  expect(await revLabels()).not.toContain('record-2 reversed relative to source');
  const manualGeometry = await geometry();
  expect(Math.max(...manualGeometry.offsets.map(({ offset }) => Math.abs(offset)))).toBeLessThanOrEqual(0.5);
  expect(manualGeometry.targetRibbons).toBeGreaterThan(0);
  console.log('Manual Reverse anchor offsets:', JSON.stringify(manualGeometry));
  await testInfo.attach('manual-anchor-offsets', { body: JSON.stringify(manualGeometry, null, 2), contentType: 'application/json' });

  // Save the reversed matched result and verify a fresh app loads it without reconstructing Python.
  await review();
  await direction.check();
  console.log('S02 Gallery: Apply requested');
  await applyReviewed();
  await expect(dialog).toBeHidden({ timeout: 180000 });
  const saved = await recordState();
  console.log('S02 Gallery: saving receipt');
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'record-owned-direction'; });
  const downloadPromise = page.waitForEvent('download', { timeout: 180000 });
  await page.evaluate(() => window.__GBDRAW_APP__.saveSessionWithTitle());
  const download = await downloadPromise;
  const savedPath = testInfo.outputPath('record-owned-direction.gbdraw-session.json.gz');
  await download.saveAs(savedPath);
  const freshContext = await browser.newContext();
  const freshPage = await freshContext.newPage();
  await captureUnhandledRejections(freshPage);
  await freshPage.goto(new URL('/gbdraw/web/index.html', page.url()).href, { waitUntil: 'domcontentloaded' });
  await freshPage.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(freshPage, readFileSync(savedPath), 'record-owned-direction.gbdraw-session.json.gz');
  const loaded = await recordState(freshPage);
  expect(loaded).toMatchObject({ plan: saved.plan, translations: saved.translations,
    orientations: saved.orientations, results: saved.results });
  expect(loaded.receipt).toEqual(saved.receipt);
  expect(saved.receipt.directions).toHaveLength(1);
  // Both fresh Load Reset scopes consume the original evidence and start no LOSAT.
  for(const scope of ['positions','positions-and-directions']) {
    const observed=await freshPage.evaluate(async scope=>{
      const {state}=await import('./js/state.js');
      let jobs=0;
      globalThis.__GBDRAW_LOSAT_EXECUTOR__=async batch=>{jobs+=batch.length;throw new Error('Reset must reuse comparisons without a LOSAT job.');};
      const historyBefore=window.__GBDRAW_HISTORY__.getUndoCount();
      const result=await window.__GBDRAW_APP__.resetSimilarityAlignment(scope);
      delete globalThis.__GBDRAW_LOSAT_EXECUTOR__;
      return {result,receipt:state.similarityAlignmentResetReceipt.value,jobs,historyBefore,historyAfter:window.__GBDRAW_HISTORY__.getUndoCount()};
    },scope);
    expect(observed.result.status).toBe('ok');expect(observed.receipt).toBeNull();
    expect(observed.jobs).toBe(0);expect(observed.historyAfter).toBe(observed.historyBefore+1);
    const reset=await recordState(freshPage);
    const expected=structuredClone(saved.orientations);
    if(scope==='positions-and-directions') for(const delta of saved.receipt.directions)expected.find(r=>r.recordKey===delta.recordKey).reverseComplement=delta.before;
    expect(reset.orientations).toEqual(expected);
    await freshPage.evaluate(()=>window.__GBDRAW_HISTORY__.undo());
    expect((await recordState(freshPage)).receipt).toEqual(saved.receipt);
  }
  const corrupt=JSON.parse(gunzipSync(readFileSync(savedPath)).toString('utf8'));
  corrupt.editorState.alignmentResetReceipt.binding='0'.repeat(64);
  const rejected=await freshPage.evaluate(async raw=>window.__GBDRAW_APP__.importSession({target:{files:[new File([raw],'corrupted.json')],value:'selected'}}),JSON.stringify(corrupt));
  expect(rejected.status).toBe('error');expect((await recordState(freshPage)).receipt).toEqual(saved.receipt);
  // Re-render the restored Session and compare the complete SVG tree.
  const regenerated = await freshPage.evaluate(() => window.__GBDRAW_APP__.runAnalysis());
  expect(regenerated.status).toBe('ok');
  const replayed = await recordState(freshPage);
  expect(replayed).toMatchObject({ plan: saved.plan, translations: saved.translations,
    orientations: saved.orientations });
  const svgEquivalent = await freshPage.evaluate(({ before, after }) => {
    const parser = new DOMParser();
    const normalized = (element) => ({
      tag: element.tagName,
      attributes: [...element.attributes]
        .filter(({ name }) => name !== 'baseProfile' && !name.startsWith('xmlns'))
        .map(({ name, value }) => [name, value]).sort(([a], [b]) => a.localeCompare(b)),
      text: [...element.childNodes].filter(({ nodeType }) => nodeType === Node.TEXT_NODE)
        .map(({ textContent }) => textContent.trim()).filter(Boolean),
      children: [...element.children].map(normalized)
    });
    return JSON.stringify(normalized(parser.parseFromString(before, 'image/svg+xml').documentElement))
      === JSON.stringify(normalized(parser.parseFromString(after, 'image/svg+xml').documentElement));
  }, { before: saved.results[0].content, after: replayed.results[0].content });
  expect(svgEquivalent).toBe(true);
  // An open review must never leak guides, badges, or controls into the SVG download.
  await review();
  const exportPromise = page.waitForEvent('download');
  await page.evaluate(() => window.__GBDRAW_APP__.downloadSVG());
  const exported = await exportPromise;
  const exportPath = testInfo.outputPath('record-owned-direction-export.svg');
  await exported.saveAs(exportPath);
  expect(readFileSync(exportPath, 'utf8'))
    .not.toMatch(/gbdraw-alignment-|Select alignment anchors|Match reference direction/);
  const downloadFormat = async (format, suffix) => {
    const pending = page.waitForEvent('download');
    await page.evaluate((method) => window.__GBDRAW_APP__[method](), `download${format}`);
    const file = await pending;
    const path = testInfo.outputPath(`record-owned-direction-${suffix}.${format.toLowerCase()}`);
    await file.saveAs(path);
    const bytes = readFileSync(path);
    return format === 'PNG' ? bytes : bytes.toString('latin1')
      .replace(/\/CreationDate \([^)]*\)/g, '')
      .replace(/\/ID \[\s*<[0-9a-f]+>\s*<[0-9a-f]+>\s*\]/gi, '');
  };
  const reviewPng = await downloadFormat('PNG', 'review');
  const reviewPdf = await downloadFormat('PDF', 'review');
  await page.keyboard.press('Escape');
  expect(await downloadFormat('PNG', 'closed')).toEqual(reviewPng);
  expect(await downloadFormat('PDF', 'closed')).toEqual(reviewPdf);
  await expectNoUnhandledRejections(page);
  await expectNoUnhandledRejections(freshPage);
  expect(pageErrors).toEqual([]);
  await freshContext.close();
});

test('exclusive directions include the minority reference and keep Custom Select/Skip local', async ({ page }, testInfo) => {
  test.setTimeout(180000);
  await page.addInitScript(() => {
    window.__ALIGNMENT_POSTS__ = [];
    const original = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function (message, ...options) {
      window.__ALIGNMENT_POSTS__.push({ type: message.type, operation: message.operation });
      return original.call(this, message, ...options);
    };
  });
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync('gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'), 'directions.gbdraw-session.json');
  const featureId = await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures
    .find(feature => feature.protein_id === 'CAG38712.1' || feature.sourceProteinId === 'CAG38712.1').svg_id);
  await page.locator(`[data-gbdraw-feature-id="${featureId}"]`).first()
    .dispatchEvent('click', { clientX: -100, clientY: -100 });
  await page.locator('.feature-popup[role="dialog"]')
    .getByRole('button', { name: 'Review alignment options…' }).click();
  const dialog = page.getByRole('dialog', { name: 'Select alignment anchors' });
  await expect(dialog).toBeVisible({ timeout: 180000 });
  for (const key of await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDraft.rows
    .filter(row => !row.choice).map(row => row.recordKey))) {
    await dialog.locator(`[data-alignment-record-key="${key}"]`)
      .getByRole('radio', { name: /Select .*bp, strand/ }).first().check();
  }
  const changed = () => page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDirectionPreview.records
    .filter(record => record.beforeReverseComplement !== record.afterReverseComplement)
    .map(record => record.recordKey));
  const before = await artifactSnapshot(page);
  const posts = await page.evaluate(() => window.__ALIGNMENT_POSTS__.length);
  await expect(dialog.getByRole('radio', { name: 'Keep current directions', exact: true })).toBeChecked();
  await expect(dialog.getByRole('combobox')).toHaveCount(0);
  const right = dialog.getByRole('radio', { name: 'All selected arrows right', exact: true });
  await right.focus();
  await page.keyboard.press('Space');
  expect(await changed()).toEqual(['record-1']);
  await dialog.getByRole('radio', { name: 'All selected arrows left', exact: true }).check();
  expect(await changed()).toEqual(['record-2', 'record-3', 'record-4', 'record-5']);
  const custom = dialog.getByRole('radio', { name: 'Custom', exact: true });
  await custom.focus();
  await page.keyboard.press('Space');
  await expect(dialog.getByRole('combobox')).toHaveCount(5);
  expect(await changed()).toEqual([]);
  await dialog.getByLabel('Direction for record-1', { exact: true }).selectOption('right');
  await dialog.getByLabel('Direction for record-2', { exact: true }).selectOption('left');
  expect(await changed()).toEqual(['record-1', 'record-2']);
  const target = dialog.locator('[data-alignment-record-key="record-2"]');
  await target.getByRole('radio', { name: /Skip / }).focus();
  await page.keyboard.press('Space');
  expect(await changed()).toEqual(['record-1']);
  await expect(dialog.getByLabel('Direction for record-2', { exact: true })).toBeDisabled();
  await expect(dialog.locator('[data-similarity-alignment-direction]').filter({ hasText: 'record-2:' }))
    .toContainText('unchanged: Skipped by user.');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDirectionPreview.records.find(record => record.recordKey === 'record-2').exclusion)).toBe('skipped_by_user');
  await target.getByRole('radio', { name: /Select .*bp, strand/ }).first().check();
  expect(await changed()).toEqual(['record-1', 'record-2']);
  await page.setViewportSize({ width: 390, height: 740 });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentCompact)).toBe(true);
  await expect(custom).toBeChecked();
  await expect(dialog.getByLabel('Direction for record-1', { exact: true })).toHaveValue('right');
  await expect(dialog.getByLabel('Direction for record-2', { exact: true })).toHaveValue('left');
  await page.setViewportSize({ width: 1600, height: 1000 });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentCompact)).toBe(false);
  expect(await artifactSnapshot(page)).toEqual(before);
  expect(await page.evaluate(() => window.__ALIGNMENT_POSTS__.length)).toBe(posts);
  await right.check();
  const apply = dialog.getByRole('button', { name: 'Apply', exact: true });
  for (let attempt = 0; attempt < 2; attempt += 1) {
    const batches = await page.evaluate(() => window.__ALIGNMENT_POSTS__.filter(post => post.type === 'helper').length);
    await apply.click();
    await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
    expect(await page.evaluate(() => window.__ALIGNMENT_POSTS__.filter(post => post.type === 'helper').length)).toBe(batches + 1);
    if (!await dialog.isVisible()) break;
    await expect(dialog.locator('[data-similarity-alignment-error]')).toContainText('Validated directions or reference placement changed');
    expect(await artifactSnapshot(page)).toEqual(before);
  }
  await expect(dialog).toBeHidden();
  const committed = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { receipt: state.similarityAlignmentResetReceipt.value,
      history: window.__GBDRAW_HISTORY__.getUndoCount() };
  });
  expect(committed.receipt.directions).toHaveLength(1);
  expect(committed.receipt.directions[0].recordKey).toBe('record-1');
  expect(committed.history).toBe(before.history[0] + 1);
  await page.screenshot({ path: testInfo.outputPath('minority-reference-right.png') });
});

const loadDirectionSession = async (page, testInfo) => {
  const directory = testInfo.outputPath('direction-fixture');
  mkdirSync(directory, { recursive: true });
  execFileSync('python', ['tests/web/helpers/alignment-direction-fixture.py', directory],
    { cwd: process.cwd(), stdio: 'pipe' });
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(join(directory, 'directions.gbdraw-session.json')), 'directions.json');
  // The fixture has admitted source-bound groups and no comparison run is needed.
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.losatProgram = 'blastp'; app.losat.blastp.mode = 'orthogroup';
    await app.refreshLinearRecordSelectors();
    await window.Vue.nextTick();
  });
};

const openDirectionReview = async page => {
  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) await page.locator('.drawer-toggle').click();
  const drawer = page.locator('.right-drawer');
  await drawer.getByRole('button', { name: 'Similarity groups' }).click();
  await drawer.locator('button').filter({ has: page.locator('.font-mono', { hasText: /^group$/ }) }).click();
  const key = await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDrawerReferenceOptions('group')
    .find(({ anchor }) => anchor.recordKey === 'record-1').key);
  await drawer.getByLabel('Exact reference record and feature').selectOption(key);
  await drawer.getByRole('button', { name: 'Review alignment options…' }).click();
  const dialog = page.getByRole('dialog', { name: 'Select alignment anchors' });
  await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentError?.message || '')).toBe('');
  await expect(dialog).toBeVisible({ timeout: 180000 });
  return dialog;
};

const applyAndSettle = async (page, dialog) => {
  await dialog.getByRole('button', { name: 'Apply', exact: true }).press('Enter');
  await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
  return page.evaluate(() => ({ status: window.__GBDRAW_APP__.similarityAlignmentStatus,
    error: window.__GBDRAW_APP__.similarityAlignmentError?.message || '' }));
};

test('minority reference, source strands and Reset scopes have truthful accessible previews', async ({ page }, testInfo) => {
  test.setTimeout(600000);
  await loadDirectionSession(page, testInfo);
  const dialog = await openDirectionReview(page);
  await dialog.getByRole('radio', { name: /Skip skip/ }).press('Space');
  const unknown = dialog.locator('#alignment-direction-record-3');
  const missing = dialog.locator('#alignment-direction-record-5');
  await dialog.getByRole('radio', { name: 'All selected arrows right', exact: true }).check();
  await expect(unknown).toContainText('no known direction; it is never guessed');
  await expect(missing).toContainText('No usable candidate');
  await expect(dialog.locator('#alignment-direction-scope')).toContainText('2 selected anchors');
  const initial = await artifactSnapshot(page);
  const preview = await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDirectionPreview);
  expect(preview.records.filter(r => r.beforeReverseComplement !== r.afterReverseComplement).map(r => r.recordKey)).toEqual(['record-1']);
  expect(preview.reference.deltaX).not.toBe(0);
  await dialog.screenshot({ path: testInfo.outputPath('minority-reference-desktop.png') });
  let outcome = await applyAndSettle(page, dialog);
  if (outcome.error.includes('Validated directions or reference placement changed')) outcome = await applyAndSettle(page, dialog);
  expect(outcome.error).toBe('');
  await expect(dialog).toBeHidden();
  const aligned = await artifactSnapshot(page);
  expect(aligned.receipt.directions).toEqual([{ recordKey: 'record-1', before: false, after: true }]);
  expect(aligned.session.resources).toEqual(initial.session.resources);
  const resetButton = page.locator('[data-similarity-alignment-reset]');
  await resetButton.click();
  const reset = page.getByRole('dialog', { name: 'Reset alignment', exact: true });
  const positions = reset.getByRole('radio', { name: 'Reset positions', exact: true });
  const combined = reset.getByRole('radio', { name: 'Reset positions and alignment direction changes', exact: true });
  await expect(reset.getByRole('heading')).toBeFocused();
  await expect(positions).toBeChecked();
  await expect(reset.getByRole('list', { name: 'Alignment direction restoration targets' })).toContainText('ref');
  await expect(reset.getByRole('list')).toContainText('current reverse; restored forward');
  await positions.focus(); await page.keyboard.press('ArrowDown');
  await expect(combined).toBeChecked(); await expect(positions).not.toBeChecked();
  await expect(reset.locator('#alignment-reset-scope')).toHaveAttribute('aria-live', 'polite');
  await expect(reset.locator('#alignment-reset-scope')).toContainText('later manual direction edits');
  for (const height of [740, 844]) {
    await page.setViewportSize({ width: 390, height });
    await expect(reset.getByRole('button', { name: 'Reset', exact: true })).toBeVisible();
    await reset.getByRole('list').scrollIntoViewIfNeeded();
    expect(await reset.evaluate(e => e.scrollWidth <= e.clientWidth)).toBe(true);
    await reset.screenshot({ path: testInfo.outputPath(`reset-390-${height}.png`) });
  }
  await combined.press('Escape'); await expect(reset).toBeHidden(); await expect(resetButton).toBeFocused();
  await resetButton.press('Enter'); await expect(positions).toBeChecked();
  await reset.getByRole('button', { name: 'Reset', exact: true }).press('Enter');
  await expect(reset).toBeHidden({ timeout: 180000 });
  await expect(resetButton).toHaveCount(0); await expect(page.locator('.drawer-toggle')).toBeFocused();
  const consumed = await artifactSnapshot(page);
  expect(consumed.plan).toBeNull(); expect(consumed.receipt).toBeNull();
  await page.evaluate(() => window.__GBDRAW_APP__.openSimilarityAlignmentReset());
  await expect(reset).toHaveCount(0);
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(() => artifactSnapshot(page)).toEqual({ ...aligned, history: [aligned.history[0], 1, expect.any(Number)] });
  expect((await artifactSnapshot(page)).history[2]).toBeGreaterThan(aligned.history[2]);
  await resetButton.click(); await combined.check();
  await reset.getByRole('button', { name: 'Reset', exact: true }).press('Enter');
  await expect(reset).toBeHidden({ timeout: 180000 });
  const restored = await artifactSnapshot(page);
  expect(restored.session.renderRequest.records).toEqual(initial.session.renderRequest.records);
  expect(restored.plan).toBeNull(); expect(restored.receipt).toBeNull();
  await page.setViewportSize({ width: 1600, height: 1000 });
  const targetOnly = await openDirectionReview(page);
  await targetOnly.getByRole('radio', { name: /Skip skip/ }).check();
  await targetOnly.getByRole('radio', { name: 'All selected arrows left', exact: true }).check();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDirectionPreview.records
    .filter(r => r.beforeReverseComplement !== r.afterReverseComplement).map(r => r.recordKey))).toEqual(['record-2']);
  outcome = await applyAndSettle(page, targetOnly);
  if (outcome.error.includes('Validated directions or reference placement changed')) outcome = await applyAndSettle(page, targetOnly);
  expect(outcome.error).toBe('');
  await expect(targetOnly).toBeHidden();
  expect((await artifactSnapshot(page)).receipt.directions).toEqual([{ recordKey: 'record-2', before: false, after: true }]);
  const keepReview = await openDirectionReview(page);
  await keepReview.getByRole('radio', { name: /Skip skip/ }).check();
  expect((await applyAndSettle(page, keepReview)).error).toBe('');
  await expect(keepReview).toBeHidden();
  expect((await artifactSnapshot(page)).receipt.directions).toEqual([]);
  await resetButton.click();
  await expect(combined).toBeDisabled(); await expect(positions).toBeEnabled();
  await expect(reset.locator('#alignment-reset-direction-reason')).toContainText('latest Align made no direction changes');
  await reset.getByRole('button', { name: 'Cancel', exact: true }).press('Enter');
  await importSession(page, readFileSync('gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'), 'historical-receipt.json');
  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) await page.locator('.drawer-toggle').click();
  await page.locator('.right-drawer').getByRole('button', { name: 'Similarity groups' }).click();
  await resetButton.click();
  await expect(combined).toBeDisabled(); await expect(positions).toBeEnabled();
  await expect(reset.locator('#alignment-reset-direction-reason')).toContainText('no historical alignment direction evidence');
  await positions.focus(); await page.keyboard.press('ArrowDown'); await expect(positions).toBeChecked();
  await positions.press('Escape'); await expect(reset).toBeHidden(); await expect(resetButton).toBeFocused();
});

test('changed final facts require a separate Apply and validation errors preserve Custom for retry', async ({ page }, testInfo) => {
  test.setTimeout(300000);
  await loadDirectionSession(page, testInfo);
  // Perturb only the initial provisional counterfactual; final batches use real Python facts.
  await page.evaluate(() => {
    const BaseWorker = window.Worker;
    window.__s03FinalBatches = 0; window.__s03PerturbPreview = true;
    window.Worker = class extends BaseWorker {
      constructor(...args) {
        super(...args);
        this.requests = new Set();
        this.addEventListener('message', event => {
          if (event.data.type !== 'helper' || !this.requests.delete(event.data.requestId)) return;
          if (window.__s03PerturbPreview && event.data.ok) {
            window.__s03PerturbPreview = false;
            const response = event.data.result;
            const fact = response.projection.records.find(r => r.recordKey === response.reference.recordKey);
            fact.variants[Number(!fact.beforeReverseComplement)].anchors[0].centerX += 5;
          } else if (window.__s03ValidationFailure) {
            window.__s03ValidationFailure = false;
            event.data.ok = false; event.data.error = { message: 'Forced final batch validation failure.' };
          }
        });
      }
      postMessage(message, transfer) {
        if (message.type === 'helper' && message.operation === 'resolveSimilarityAlignment') {
          this.requests.add(message.requestId);
          if (message.payload.projection.orientations !== null) window.__s03FinalBatches++;
        }
        return super.postMessage(message, transfer);
      }
    };
  });
  const dialog = await openDirectionReview(page);
  await page.evaluate(() => { window.__s03FinalBatches = 0; });
  await dialog.getByRole('radio', { name: /Skip skip/ }).check();
  await dialog.getByRole('radio', { name: 'Custom', exact: true }).check();
  await dialog.getByLabel('Direction for record-1', { exact: true }).selectOption('right');
  const before = await artifactSnapshot(page);
  const oldDelta = await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDirectionPreview.reference.deltaX);
  await applyAndSettle(page, dialog);
  await expect(dialog.locator('[data-similarity-alignment-error]')).toContainText('Review the updated preview and Apply again');
  expect(await artifactSnapshot(page)).toEqual(before);
  expect(await page.evaluate(() => window.__s03FinalBatches)).toBe(1);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDirectionPreview.reference.deltaX)).not.toBe(oldDelta);
  await page.evaluate(() => { window.__s03ValidationFailure = true; });
  await applyAndSettle(page, dialog);
  await expect(dialog.locator('[data-similarity-alignment-error]')).toContainText('Forced final batch validation failure');
  await expect(dialog.getByLabel('Direction for record-1', { exact: true })).toHaveValue('right');
  expect(await artifactSnapshot(page)).toEqual(before);
  expect(await page.evaluate(() => window.__s03FinalBatches)).toBe(2);
  await applyAndSettle(page, dialog); await expect(dialog).toBeHidden();
  expect(await page.evaluate(() => window.__s03FinalBatches)).toBe(3);
  expect((await artifactSnapshot(page)).history[0]).toBe(before.history[0] + 1);
});
