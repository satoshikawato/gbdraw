const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { mkdirSync, readFileSync } = require('node:fs');
const { join } = require('node:path');
const { gunzipSync } = require('node:zlib');
const { CURRENT_FEATURE_CATALOG_SCHEMA, CURRENT_REQUEST_SCHEMA, CURRENT_SESSION_VERSION, answerLoadRepairs, evaluateWithRetainedPromise, openLinearFileCards } = require('./helpers/app-lifecycle.cjs');
const { showPhonePane, paneTop } = require('./helpers/drawing-bar.cjs');

const importSession = async (page, bytes, name) => evaluateWithRetainedPromise(page, async ({ bytes, name }) => {
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

// Record keys are binding-only; the visible names come from biological metadata.
const directionSelect = (dialog, recordKey) => dialog.locator(`select[aria-describedby="alignment-direction-${recordKey}"]`);
const directionRow = (dialog, recordKey) => dialog.locator(
  `[data-similarity-alignment-direction]:has(#alignment-direction-${recordKey})`);
const artifactSnapshot = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const history = window.__GBDRAW_HISTORY__;
  const { getCommittedCanonicalSession } = await import('./js/services/config.js');
  return {
    session: JSON.parse(JSON.stringify(getCommittedCanonicalSession())),
    plan: JSON.parse(JSON.stringify(state.similarityAlignmentPlan.value)),
    receipt: JSON.parse(JSON.stringify(state.similarityAlignmentResetReceipt.value)),
    request: structuredClone(state.lastCommittedRequest?.value || null),
    results: window.__GBDRAW_APP__.currentResults().map(({ name, content }) => ({ name, content })),
    history: [history.getUndoCount(), history.getRedoCount(), history.revision.value]
  };
});

// The review is docked in the Result. Rings mark the candidates on the figure, and a
// card lists its candidates in figure order (left to right, undrawn ones last).
const reviewCard = (page, recordKey) => page.evaluate((key) => {
  const card = window.__GBDRAW_APP__.similarityAlignmentReview.groups
    .flatMap((group) => group.cards).find((item) => item.recordKey === key);
  return { chosenKey: card.chosenKey,
    candidates: card.candidates.map(({ key, anchor, name }) => ({ key, anchor, name })) };
}, recordKey);
const candidateRing = (page, candidateKey) => page.locator(
  `button.gbdraw-alignment-ring[data-alignment-candidate-key=${JSON.stringify(candidateKey)}]`);
const previewView = (page) => page.evaluate(() => ({
  zoom: window.__GBDRAW_APP__.zoom, pan: { ...window.__GBDRAW_APP__.canvasPan }
}));
const editorDrawer = (page) => page.evaluate(() => ({
  open: window.__GBDRAW_APP__.showRightDrawer, tab: window.__GBDRAW_APP__.rightDrawerTab
}));
// Whether the view is what Fit gives in the canvas as it is now: the largest whole-percent zoom
// that leaves 8 px around the figure, with the figure centred. `room` is the canvas size.
const figureFit = (page) => page.evaluate(() => {
  const { zoom } = window.__GBDRAW_APP__;
  const canvas = document.querySelector('.preview-canvas');
  const frame = canvas.getBoundingClientRect();
  const figure = document.querySelector('.gbdraw-preview-surface svg').getBoundingClientRect();
  const largest = Math.min((canvas.clientWidth - 16) / (figure.width / zoom), (canvas.clientHeight - 16) / (figure.height / zoom));
  return {
    room: [canvas.offsetWidth, canvas.offsetHeight],
    zoom,
    fitted: Math.abs(zoom - Math.floor(largest * 100 + 1e-9) / 100) < 1e-9
      && Math.abs(figure.left + figure.right - frame.left - frame.right) <= 2
      && Math.abs(figure.top + figure.bottom - frame.top - frame.bottom) <= 2,
    inFrame: figure.left >= frame.left && figure.right <= frame.right
      && figure.top >= frame.top && figure.bottom <= frame.bottom
  };
});
// A change of the canvas size has reached its observers once two frames have been drawn.
const roomSettled = (page) => evaluateWithRetainedPromise(page, () => new Promise((resolve) => {
  requestAnimationFrame(() => requestAnimationFrame(() => resolve(true)));
}));

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
    const document = JSON.parse(raw);
    document.renderRequest.comparisons = [];
    document.editorState.featureCatalog.items.forEach((item) => {
      item.orthogroups.forEach((group) => {
        if (group.id === 'og_1') group.orthologEdges = [];
      });
    });
    // The Python writer has no Web draft: Load takes the request's projection,
    // and the Linear slice (Session 46) holds the comparison draft of this test.
    document.modes = { linear: { config: {
      linearComparisonPlan: { mode: 'adjacent', defaultSource: 'losat', edges: [] },
      losatProgram: 'blastp',
      losat: { blastp: { mode: 'orthogroup' } }
    } } };
    return JSON.stringify(document);
  }, source));
  await importSession(page, session, 'similarity-alignment-inparalog.gbdraw-session.json');
};

test('alignment rings and the reference line stay transient and share the panel choices', async ({ page }, testInfo) => {
  test.setTimeout(180000);
  await loadAmbiguousSession(page, testInfo);
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 's03-canvas-overlay'; });
  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) {
    await page.locator('.drawer-toggle').click();
  }
  const drawer = page.locator('.right-drawer');
  await drawer.getByRole('button', { name: 'Similarity groups' }).click();
  await drawer.locator('[data-similarity-group-row="og_1"]').click();
  const referenceSelect = drawer.getByLabel('Exact reference record and feature');
  const referenceKey = await page.evaluate(() => (
    window.__GBDRAW_APP__.similarityAlignmentDrawerReferenceOptions('og_1')
      .find(({ anchor }) => anchor.recordKey === 'record_b')?.key || ''
  ));
  await referenceSelect.selectOption(referenceKey);
  const align = drawer.getByRole('button', { name: 'Align…', exact: true });
  const before = await artifactSnapshot(page);
  await align.click();
  const dialog = page.getByRole('dialog', { name: 'Review alignment' });
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await expect(dialog.locator('[data-alignment-application-help]')).toContainText('Apply required');
  await expect(page.locator('[data-generation-application-feedback]')).toHaveCount(0);
  const overlay = page.locator('[data-similarity-alignment-canvas]');
  const line = overlay.locator('.gbdraw-alignment-line:not(.casing)');
  const rings = overlay.locator('button.gbdraw-alignment-ring');
  await expect(line).toHaveCount(1);
  await expect(rings).toHaveCount(2);
  await expect(overlay.locator('.gbdraw-alignment-flag')).toHaveText('b0 · reference');
  // A ring and its row share the candidate's name; no number joins them.
  await expect(overlay.locator('.gbdraw-alignment-badge')).toHaveCount(0);
  await expect(dialog.locator('[title="Number on the diagram"]')).toHaveCount(0);
  const [first, second] = (await reviewCard(page, 'record_a')).candidates;
  await expect(candidateRing(page, first.key)).toHaveAttribute('aria-label', /^Select .+ for .+$/);
  await expect(candidateRing(page, first.key)).toHaveText('');

  const geometry = () => page.evaluate(async () => {
    const { getFeatureFillElements } = await import('./js/services/feature-dom.js');
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
    const centerOf = (element) => {
      const rect = element.getBoundingClientRect();
      return { x: (rect.left + rect.right) / 2, y: (rect.top + rect.bottom) / 2 };
    };
    const reference = center(draft.response.reference);
    const referenceLine = document.querySelector('.gbdraw-alignment-line:not(.casing)');
    const ringErrors = draft.rows[0].candidates.map((candidate) => {
      const ring = [...document.querySelectorAll('button.gbdraw-alignment-ring')]
        .find((element) => element.dataset.alignmentCandidateKey === candidate.key);
      if (!ring || ring.hidden) return null;
      const point = center(candidate.anchor);
      const ringCenter = centerOf(ring);
      return Math.max(Math.abs(point.x - ringCenter.x), Math.abs(point.y - ringCenter.y));
    }).filter((value) => value !== null);
    return {
      lineError: referenceLine ? Math.abs(reference.x - centerOf(referenceLine).x) : null,
      ringErrors,
      visibleRings: ringErrors.length
    };
  });
  // The figure was just fitted beside the panel: both rings and the line are in the frame.
  await expect.poll(async () => {
    const measured = await geometry();
    return measured.visibleRings === 2 && measured.lineError !== null && measured.lineError < 2
      && measured.ringErrors.every((error) => error < 2);
  }).toBe(true);
  // A zoom or a pan can move a candidate or the reference out of the frame, which hides its
  // ring or the line; what stays in the frame stays on its feature.
  const geometryIsAligned = async () => {
    const measured = await geometry();
    return measured.visibleRings + (measured.lineError === null ? 0 : 1) > 0
      && (measured.lineError === null || measured.lineError < 2)
      && measured.ringErrors.every((error) => error < 2);
  };
  await page.screenshot({ path: testInfo.outputPath('canvas-overlay.png'), fullPage: true });

  const radios = dialog.getByRole('radio', { name: /Select .*bp, strand/ });
  const firstRadio = radios.first();
  const secondRadio = radios.nth(1);
  const firstRow = firstRadio.locator('xpath=ancestor::label');
  // A ring, its row and the figure name light together.
  await candidateRing(page, first.key).hover();
  await expect(firstRow).toHaveClass(/is-hover/);
  await expect(overlay.locator('.gbdraw-alignment-tag')).toHaveText(first.name);
  await page.mouse.move(5, 5);
  await expect(firstRow).not.toHaveClass(/is-hover/);
  const firstSvgId = await page.evaluate((anchor) => window.__GBDRAW_APP__.extractedFeatures
    .find((item) => item.recordKey === anchor.recordKey
      && item.biologicalFeatureId === anchor.biologicalFeatureId).svg_id, first.anchor);
  const firstFeature = page.locator('[data-gbdraw-feature-id="' + firstSvgId + '"]').first();
  await firstFeature.dispatchEvent('mouseover');
  await expect(firstRow).toHaveClass(/is-hover/);
  await firstFeature.dispatchEvent('click', { clientX: -100, clientY: -100 });
  await expect(firstRadio).toBeChecked();
  await expect(page.locator('.feature-popup[role="dialog"]')).toHaveCount(0);
  await candidateRing(page, second.key).click();
  await expect(secondRadio).toBeChecked();
  await expect(candidateRing(page, second.key)).toHaveAttribute('aria-pressed', 'true');
  await expect(candidateRing(page, first.key)).toHaveAttribute('aria-pressed', 'false');
  await page.mouse.move(5, 5);
  const secondRow = secondRadio.locator('xpath=ancestor::label');
  await firstFeature.dispatchEvent('mouseover');
  await secondRow.dispatchEvent('mouseenter');
  await expect.poll(() => firstFeature.evaluate((element) => element.style.opacity))
    .not.toBe('0.7');
  await secondRow.dispatchEvent('mouseleave');
  await expect.poll(() => firstFeature.evaluate((element) => element.style.opacity))
    .toBe('0.7');
  await expect(secondRadio).toBeChecked();
  // No gesture on the figure opens the feature popup while the review is open.
  const referenceSvgId = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const anchor = app.similarityAlignmentDraft.response.reference;
    return app.extractedFeatures.find((item) => item.recordKey === anchor.recordKey
      && item.biologicalFeatureId === anchor.biologicalFeatureId).svg_id;
  });
  await page.locator('[data-gbdraw-feature-id="' + referenceSvgId + '"]').first()
    .dispatchEvent('click', { clientX: -100, clientY: -100 });
  await expect(page.locator('.feature-popup[role="dialog"]')).toHaveCount(0);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.clickedFeature)).toBeNull();
  expect(await artifactSnapshot(page)).toEqual(before);

  const canvas = page.locator('.gbdraw-preview-surface').locator('..');
  const zoomBefore = await page.evaluate(() => window.__GBDRAW_APP__.zoom);
  await canvas.dispatchEvent('wheel', { deltaY: -80 });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.zoom)).toBeGreaterThan(zoomBefore);
  await expect.poll(geometryIsAligned).toBe(true);
  // The Zoom in button moves the figure through a transition; the overlay follows it.
  await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
  await expect.poll(geometryIsAligned).toBe(true);
  const panBefore = await page.evaluate(() => window.__GBDRAW_APP__.canvasPan.x);
  // The feature shapes and the rings own their pointer input. Pan from a real
  // exposed canvas margin rather than a fixed point that may hit either one.
  const panPoint = await canvas.evaluate((element) => {
    const box = element.getBoundingClientRect();
    for (const y of [box.top + 4, box.bottom - 4, box.top + 20, box.bottom - 20]) {
      for (let x = box.left + 4; x < box.right - 40; x += 12) {
        const hit = document.elementFromPoint(x, y);
        if (hit && element.contains(hit) && !hit.closest('svg')) return { x, y };
      }
    }
    return null;
  });
  expect(panPoint).not.toBeNull();
  await page.mouse.move(panPoint.x, panPoint.y);
  await page.mouse.down();
  await page.mouse.move(panPoint.x + 30, panPoint.y, { steps: 3 });
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
  await expect(candidateRing(page, first.key)).toBeHidden();
  // A candidate the Result does not show is listed last; its row still chooses it.
  const hiddenRadio = radios.and(page.locator(`[aria-label=${JSON.stringify(await firstRadio.getAttribute('aria-label'))}]`));
  await hiddenRadio.check();
  await expect(hiddenRadio).toBeChecked();
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
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.saveSessionWithTitle());
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
  // Another Result under the open review, as a Generate or a Session load draws it: it commits
  // a Session with its own resources, even from the same files. The rings and the line belong
  // to the sources the review was made on, so they go. (A Result of the same sources, as a live
  // setting draws it, keeps the resources and them: OV-1082, OV-1120.)
  await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const { canonicalRenderArtifactOwner } = await import('./js/services/config.js');
    window.__S03Result = window.__GBDRAW_APP__.currentResults()[state.selectedResultIndex.value];
    const canonical = window.__S03Canonical = canonicalRenderArtifactOwner.capture();
    const committed = canonical.committedCanonicalSession;
    if (!committed?.resources) throw new Error('The fixture has no committed resources.');
    canonicalRenderArtifactOwner.restore({ ...canonical,
      committedCanonicalSession: { ...committed, resources: { ...committed.resources } } });
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
    const { canonicalRenderArtifactOwner } = await import('./js/services/config.js');
    canonicalRenderArtifactOwner.restore(window.__S03Canonical);
    state.results.value = [window.__S03Result];
    delete window.__S03Result;
    delete window.__S03Canonical;
  });
  await expect(page.locator('.gbdraw-preview-surface svg')).toBeVisible();
  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) {
    await page.locator('.drawer-toggle').click();
  }
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
  await drawer.locator('[data-similarity-group-row="og_1"]').click();

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
  const dialog = page.getByRole('dialog', { name: 'Review alignment' });
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await expect(dialog.getByRole('heading', { name: 'Review alignment' })).toBeFocused();
  // The description counts what the section headings count, then says what Apply means.
  await expect(dialog).toHaveAttribute('aria-describedby',
    'similarity-alignment-dialog-description similarity-alignment-application-help');
  const sectionCounts = await dialog.locator('[data-alignment-section]').evaluateAll((headings) => Object.fromEntries(
    headings.map((heading) => [heading.dataset.alignmentSection, Number(heading.querySelector('.count').textContent)])));
  const summaryParts = [
    ['choose', 'awaiting a choice'], ['several', 'with several candidates'],
    ['single', 'with one candidate'], ['none', 'without a usable candidate']
  ].filter(([id]) => sectionCounts[id]).map(([id, words]) => `${sectionCounts[id]} ${sectionCounts[id] === 1 ? 'record' : 'records'} ${words}`);
  const recordCount = await dialog.locator('.alignment-review-card').count();
  await expect(dialog.locator('#similarity-alignment-dialog-description')).toHaveText(
    `Similarity group ${await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDraft.response.groupId)}: `
    + [`${recordCount} records`, ...summaryParts].join(', ') + '.');
  const applyChip = dialog.locator('#similarity-alignment-application-help');
  await expect(applyChip).toBeVisible();
  await expect(applyChip).toHaveText('Apply required: nothing changes until Apply succeeds. Cancel keeps the current Result.');
  // The reference gene is placed in text as a candidate row is.
  await expect(dialog.locator('[data-similarity-alignment-reference] .alignment-review-place')).toHaveText(/^[\d,]+\.\.[\d,]+ bp \(.+\)$/);
  // The reference card names its record and gene; it says more only when its record is reversed.
  const referenceCard = dialog.locator('[data-similarity-alignment-reference]');
  await expect(referenceCard).toContainText('Reference');
  await expect(referenceCard).toContainText('b0');
  await expect(referenceCard).toContainText('record_b');
  await expect(referenceCard.locator('[data-similarity-alignment-reference-correction]')).toHaveCount(0);
  await expect(dialog.locator('[data-alignment-record-key]')).toHaveCount(1);
  await expect(dialog.locator('[data-similarity-alignment-reason]')).toContainText(
    'Recommended: a0, the only group representative among the candidates.'
  );
  await dialog.getByRole('radio', {name:/Select .*bp, strand/}).first().check();
  await expect(dialog.getByText('Recommended', { exact: true })).toHaveCount(1);
  await expect(dialog).not.toHaveAttribute('aria-modal', 'true');
  expect(await dialog.evaluate((element) => element.contains(document.activeElement))).toBe(true);
  // The review took the Result: the Editor that opened it is closed until the review ends.
  expect(await editorDrawer(page)).toEqual({ open: false, tab: 'orthogroups' });
  await expect(page.locator('.drawer-toggle')).toBeDisabled();
  await expect(page.locator('.drawer-toggle')).toHaveAttribute('title', /Finish or cancel alignment review/);
  // Below the dock width the panel sits under the figure, which keeps a usable height.
  await page.setViewportSize({ width: 900, height: 650 });
  await expect.poll(async () => (await page.locator('.preview-canvas').boundingBox()).height)
    .toBeGreaterThanOrEqual(200);
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
  const removalDialog = page.getByRole('dialog', { name: /^Remove .*\?$/ });
  await expect(removalDialog).toBeVisible();
  await removalDialog.getByRole('button', { name: 'Cancel' }).focus();
  await page.keyboard.press('Escape');
  await expect(removalDialog).toBeHidden();
  await expect(dialog).toBeVisible();
  expect(await artifactSnapshot(page)).toEqual(beforeCancel);
  await expect(dialog).toContainText(/\d[\d,]*\.\.\d[\d,]* bp \([+-]\)/);
  await expect(dialog).toContainText('no direct hit');
  await expect(dialog).not.toContainText(/score/i);
  const apply = dialog.getByRole('button', { name: 'Apply', exact: true });
  await expect(apply).toBeEnabled();
  await expect(apply).toHaveAttribute('aria-describedby', 'similarity-alignment-apply-reason');
  await expect(page.locator('#similarity-alignment-apply-reason')).toHaveText('Nothing changes until Apply.');
  await page.screenshot({
    path: testInfo.outputPath('desktop-ambiguity.png'),
    fullPage: true
  });

  const firstCandidate = dialog.getByRole('radio', { name: /Select .*bp, strand/ }).first();
  const candidateAnchor = (await reviewCard(page, 'record_a')).candidates[0].anchor;
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
  await dialog.getByRole('heading', { name: 'Review alignment' }).hover();
  await expect.poll(async () => (await previewState()).candidate.includes('0.7')).toBe(false);
  const replacementCandidate = dialog.getByRole('radio', { name: /Select .*bp, strand/ }).nth(1);
  await replacementCandidate.check();
  await expect(replacementCandidate).toBeChecked();
  // The recommendation stays marked on its candidate whatever is selected.
  await expect(dialog.getByText('Recommended', { exact: true })).toHaveCount(1);
  await expect(firstCandidate.locator('xpath=ancestor::label')).toContainText('Recommended');
  await expect(dialog.locator('[data-similarity-alignment-reason]')).toContainText(
    'Recommended: a0, the only group representative among the candidates.'
  );
  await firstCandidate.check();
  await expect(dialog.getByText('Recommended', { exact: true })).toHaveCount(1);

  await dialog.getByLabel(/^Skip /).focus();
  await page.keyboard.press('Space');
  await expect(dialog).toContainText('Unchanged: Skipped by user.');
  await expect(dialog.getByText('Recommended', { exact: true })).toHaveCount(1);
  await expect(apply).toBeEnabled();
  await apply.focus();
  await page.keyboard.press('Tab');
  expect(await dialog.evaluate((element) => element.contains(document.activeElement))).toBe(false);
  await page.keyboard.press('Escape');
  await expect(dialog).toBeHidden();
  // The Editor that was open when the review started is open again on its tab, and the
  // focus is back on the control that opened the review.
  expect(await editorDrawer(page)).toEqual({ open: true, tab: 'orthogroups' });
  await expect(page.locator('.drawer-toggle')).toBeEnabled();
  await expect(align).toBeFocused();
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
  // The popup that opened the review is hidden, not closed, while the review is open.
  await expect(popup).toBeHidden();
  await expect(popup).toHaveCount(1);
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
  await expect(popup).toBeVisible();
  await expect(popupReview).toBeFocused();
  await popup.getByRole('button', { name: 'Close feature popup' }).click();
  expect(await editorDrawer(page)).toEqual({ open: true, tab: 'orthogroups' });

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
  // UX-B06: one line until opened; the reference and each record's rationale are inside.
  await expect(inspector.locator('summary')).toContainText(/^Aligned to .+ in .+$/);
  await expect(inspector.getByText('Exact reference:')).toBeHidden();
  await inspector.locator('summary').click();
  await expect(inspector.getByText('Exact reference:')).toBeVisible();
  await expect(inspector).toContainText('Selected by user');
  await expect(inspector).toContainText('Exact reference');
  await page.screenshot({ path: testInfo.outputPath('narrow-final.png'), fullPage: true });

  const savedState = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return {
      plan: JSON.parse(JSON.stringify(state.similarityAlignmentPlan.value)),
      translations: JSON.parse(JSON.stringify(state.linearRecordTranslations.value)),
      orientations: state.drawingList().find((drawing) => drawing.mode === 'linear').linearSeqs.map(({ uid, region_reverse: reverse }) => ({
        recordKey: uid, reverse: Boolean(reverse)
      })),
      results: window.__GBDRAW_APP__.currentResults().map(({ name, content }) => ({ name, content }))
    };
  });
  expect(Object.keys(savedState.plan.records.find(({ status }) => status === 'aligned')).sort())
    .toEqual(['anchor', 'rationale', 'recordKey', 'status']);
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 's06-active-alignment'; });
  const downloadPromise = page.waitForEvent('download', { timeout: 180000 });
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.saveSessionWithTitle());
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
      orientations: state.drawingList().find((drawing) => drawing.mode === 'linear').linearSeqs.map(({ uid, region_reverse: reverse }) => ({
        recordKey: uid, reverse: Boolean(reverse)
      })),
      results: window.__GBDRAW_APP__.currentResults().map(({ name, content }) => ({ name, content })),
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
  const regenerate = await evaluateWithRetainedPromise(freshPage, async () => {
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
      results: window.__GBDRAW_APP__.currentResults().map(({ name, content }) => ({ name, content }))
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
  await freshPage.getByRole('dialog',{name:'Reset alignment',exact:true}).getByRole('button',{name:'Reset positions',exact:true}).click();
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
  await drawer.locator('[data-similarity-group-row="og_1"]').click();
  const reference = await page.evaluate(() => (
    window.__GBDRAW_APP__.similarityAlignmentDrawerReferenceOptions('og_1')
      .find(({ anchor }) => anchor.recordKey === 'record_b')?.key || ''
  ));
  await drawer.getByLabel('Exact reference record and feature').selectOption(reference);
  await drawer.getByRole('button', { name: 'Align…', exact: true }).click();
  const dialog = page.getByRole('dialog', { name: 'Review alignment' });
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
  const modes = dialog.getByRole('group', { name: 'Arrow direction', exact: true });
  await expect(modes.getByRole('radio')).toHaveCount(4);
  await expect(modes.getByRole('radio', { name: 'Keep as is', exact: true })).toBeChecked();
  await expect(dialog.getByRole('combobox')).toHaveCount(0);
  const first = dialog.getByRole('radio', { name: /Select .*bp, strand/ }).first();
  const skip = dialog.getByLabel(/^Skip /);
  await first.focus(); await page.keyboard.press('Space');
  await expect(first).toBeChecked();
  await page.keyboard.press('ArrowDown');
  await expect(dialog.getByRole('radio', { name: /Select .*bp, strand/ }).nth(1)).toBeChecked();
  await page.keyboard.press('ArrowDown'); await expect(skip).toBeChecked();
  await expect(dialog.locator('#alignment-direction-scope')).toContainText('1 record stays unchanged');
  await modes.getByRole('radio', { name: 'Keep as is', exact: true }).focus();
  await page.keyboard.press('ArrowRight');
  await expect(modes.getByRole('radio', { name: 'All right', exact: true })).toBeChecked();
  await page.keyboard.press('ArrowRight');
  await expect(modes.getByRole('radio', { name: 'All left', exact: true })).toBeChecked();
  await page.keyboard.press('ArrowRight');
  await expect(modes.getByRole('radio', { name: 'Custom', exact: true })).toBeChecked();
  await expect(modes.locator('input:checked')).toHaveCount(1);
  const target = directionSelect(dialog, 'record_a');
  await expect(target).toBeDisabled();
  await expect(dialog.locator('#alignment-direction-record_a')).toContainText('Skipped by user');
  await first.focus(); await page.keyboard.press('Space');
  await expect(target).toBeEnabled();
  await expect(dialog.locator('#alignment-direction-scope')).toContainText('No record is reversed yet');
  const reference = directionSelect(dialog, 'record_b');
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
  await modes.getByRole('radio', { name: 'Keep as is', exact: true }).check();
  await expect(dialog.getByRole('combobox')).toHaveCount(0);
  expect(await page.evaluate(() => window.__s03Jobs)).toBe(0);
  expect(await artifactSnapshot(page)).toEqual(before);
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).focus();
  await page.keyboard.press('Enter'); await expect(dialog).toBeHidden();
  // The review was opened from the Editor, which is open again with the focus on its button.
  await expect(page.locator('.right-drawer').getByRole('button', { name: 'Align…', exact: true })).toBeFocused();
});

test('Apply error retains draft, focuses retry guidance, and accepts correction', async ({ page }, testInfo) => {
  test.setTimeout(180000);
  await loadAmbiguousSession(page, testInfo);
  const dialog = await openAlignmentFromDrawer(page);
  await dialog.getByRole('radio',{name:/Select .*bp, strand/}).first().check();
  const before = await artifactSnapshot(page);
  await expect(dialog.locator('[data-alignment-application-help]')).toContainText('Apply required');
  const application = page.locator('[data-generation-application-feedback]');
  await dialog.getByLabel(/^Skip /).check();
  await expect(application).toHaveCount(0);
  expect(await artifactSnapshot(page)).toEqual(before);
  await dialog.getByRole('radio', { name: /Select .*bp, strand/ }).first().check();
  const priorLayout = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const prior = app.form.linear_track_layout;
    const sanitize=window.DOMPurify.sanitize;window.DOMPurify.sanitize=(...args)=>{if(String(args[0]).includes('<svg')){window.DOMPurify.sanitize=sanitize;throw Object.assign(new Error('PRIVATE_CAUSE_SENTINEL'), { code: 'COMPARISON_IDENTITY', stage: 'result-admission', context: { reason: 'SOURCE_VIEW_CONFLICT' }, cause: new Error('PRIVATE_CAUSE_SENTINEL'), stdout: 'PRIVATE_STDOUT_SENTINEL', secondary: [{code:'CLEANUP_FAILED',stage:'cleanup',message:'PRIVATE_CLEANUP_SENTINEL'}] });}return sanitize(...args);};
    return prior;
  });
  const fitted = await previewView(page);
  await dialog.getByRole('button', { name: 'Apply', exact: true }).click();
  const error = dialog.locator('[data-similarity-alignment-error]');
  await expect(error).toBeVisible({ timeout: 180000 });
  await expect(error).toContainText('The features at the ends of comparison links do not match');
  // The failed run put the Result back; the review goes on with the figure fitted beside it.
  await expect.poll(async () => {
    const view = await previewView(page);
    return view.zoom === fitted.zoom && Math.abs(view.pan.x - fitted.pan.x) < 0.01
      && Math.abs(view.pan.y - fitted.pan.y) < 0.01;
  }).toBe(true);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentError)).toMatchObject({code:'COMPARISON_IDENTITY',operation:'align',stage:'result-admission',secondary:[{code:'CLEANUP_FAILED',stage:'cleanup'}]});
  await expect(page.getByRole('alert', { name: 'Alignment error' })).toHaveCount(1);
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
  // Candidate sanitization failed; the independent settings draft was not changed.
  await expect(application).toHaveCount(0);
  const afterFailure = await artifactSnapshot(page);
  expect(afterFailure.results).toEqual(before.results);
  expect(afterFailure.history.slice(0, 2)).toEqual(before.history.slice(0, 2));
  await page.evaluate((layout) => {
    window.__GBDRAW_APP__.form.linear_track_layout = layout;
  }, priorLayout);
  await dialog.getByLabel(/^Skip /).check();
  await expect(error).toBeHidden();
  await dialog.getByRole('button', { name: 'Apply', exact: true }).click();
  await expect(dialog).toBeHidden({ timeout: 180000 });
  await expect(page.locator('[data-similarity-alignment-summary]')).toContainText('explicitly skipped');
  await expect(application).toHaveCount(0);
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
    // The Python writer has no Web draft: Load takes the request's projection,
    // and the Linear slice (Session 46) holds the comparison draft of this test.
    document.modes = { linear: { config: {
      linearComparisonPlan: { mode: 'adjacent', defaultSource: 'losat', edges: [] },
      losatProgram: 'blastp',
      losat: { blastp: { mode: 'orthogroup' } }
    } } };
    return JSON.stringify(document);
  }, source));
  await importSession(page, session, 'one-candidate.gbdraw-session.json');

  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) {
    await page.locator('.drawer-toggle').click();
  }
  const drawer = page.locator('.right-drawer');
  await drawer.getByRole('button', { name: 'Similarity groups' }).click();
  await drawer.locator('[data-similarity-group-row="og_1"]').click();
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
  const resolvedReview = page.getByRole('dialog', { name: 'Review alignment' });
  await expect(resolvedReview).toBeVisible({ timeout: 180000 });
  await expect(resolvedReview.getByRole('heading', { name: 'Review alignment' })).toBeFocused();
  await expect(resolvedReview.locator('[data-similarity-alignment-reference]')).toContainText('b0');
  await expect(resolvedReview.locator('[data-alignment-record-key]')).toHaveCount(1);
  await expect(resolvedReview.getByRole('button', { name: 'Apply', exact: true })).toBeEnabled();
  await expect(resolvedReview.getByRole('checkbox', { name: /Match reference direction for/ }))
    .toHaveCount(0);
  // One usable candidate: a compact card whose Skip checkbox is the only choice.
  await expect(resolvedReview.getByRole('radio', { name: /Select .*bp, strand/ })).toHaveCount(0);
  await expect(resolvedReview.locator('[data-alignment-section="single"]')).toHaveText('One candidate1 (1 record)');
  await expect(resolvedReview.locator('[data-alignment-record-key] .alignment-review-place')).toHaveText(/^[\d,]+\.\.[\d,]+ bp \(.+\)$/);
  await expect(resolvedReview.locator('[data-alignment-record-key]')).toContainText('a0');
  // The card's second line: the evidence at the left, Skip at the right, level with it.
  expect(await resolvedReview.locator('[data-alignment-record-key] .alignment-review-gene').evaluate((gene) => {
    const box = (selector) => gene.querySelector(selector).getBoundingClientRect();
    const [name, chip, skip] = [box('b'), box('.alignment-review-chip.is-evidence'), box('.alignment-review-skip')];
    return { below: chip.top >= name.bottom - 1, left: Math.abs(chip.left - name.left) <= 1, right: skip.left > chip.right,
      level: Math.abs(chip.top + chip.bottom - skip.top - skip.bottom) <= 4 };
  })).toEqual({ below: true, left: true, right: true, level: true });
  await resolvedReview.getByRole('checkbox', { name: /^Skip / }).check();
  await expect(resolvedReview).toContainText('Unchanged: Skipped by user.');
  await resolvedReview.getByRole('checkbox', { name: /^Skip / }).uncheck();
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
  // A one-click Align shows no review, so nothing is taken over: the Editor is never closed.
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    window.__alignmentDrawerChanges = [];
    window.Vue.watch(() => app.showRightDrawer, (open) => window.__alignmentDrawerChanges.push(open), { flush: 'sync' });
  });
  const viewBeforeAlign = await previewView(page);
  await normalAlign.click();
  await expect(page.getByRole('dialog', { name: 'Review alignment' })).toHaveCount(0);
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
  expect(await page.evaluate(() => window.__alignmentDrawerChanges)).toEqual([]);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)).toBe(true);
  expect(await previewView(page)).toEqual(viewBeforeAlign);

  // Layout editing keeps record drag ahead of pan and consumes the active plan
  // through the existing alignment lifecycle, as one undoable action.
  if (await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) {
    await page.locator('.drawer-toggle').click();
  }
  await page.getByRole('button', { name: 'Toggle layout edit mode' }).click();
  await page.evaluate(() => { window.__GBDRAW_APP__.zoom = 0.6; });
  await expect.poll(() => page.evaluate(() => new DOMMatrix(
    getComputedStyle(window.__GBDRAW_APP__.svgContainer).transform).a)).toBeCloseTo(0.6, 5);
  const recordPoint = await page.evaluate(async () => {
    const a = window.__GBDRAW_APP__;
    const { closestRecordGroup } = await import('./js/app/record-groups.js');
    for (const line of a.svgContainer.querySelectorAll('g[data-gbdraw-record-id] line')) {
      const rect = line.getBoundingClientRect();
      const x = rect.x + rect.width / 2, y = rect.y + rect.height / 2;
      const hit = document.elementFromPoint(x, y);
      const group = closestRecordGroup(hit);
      if (group && !hit.closest('[data-gbdraw-feature-id], text[data-label-editable="true"]')) {
        return { x, y, groupId: group.id };
      }
    }
    throw new Error('No visible record axis drag target');
  });
  const beforeRecordDrag = await artifactSnapshot(page);
  const panBeforeRecordDrag = await page.evaluate(() => ({ ...window.__GBDRAW_APP__.canvasPan }));
  await page.mouse.move(recordPoint.x, recordPoint.y);
  await page.mouse.down();
  await expect(page.locator(`[id="${recordPoint.groupId}"]`)).toHaveCSS('cursor', 'grabbing');
  await page.mouse.move(recordPoint.x + 12, recordPoint.y + 6, { steps: 6 });
  await page.mouse.up();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount()))
    .toBe(beforeRecordDrag.history[0] + 1);
  const afterRecordDrag = await artifactSnapshot(page);
  expect(afterRecordDrag.plan).toBeNull();
  expect(afterRecordDrag.results).not.toEqual(beforeRecordDrag.results);
  expect(await page.evaluate(() => ({ ...window.__GBDRAW_APP__.canvasPan }))).toEqual(panBeforeRecordDrag);
  expect(await page.evaluate(() => window.__alignmentWorkerResolves)).toBe(1);
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  expect((await artifactSnapshot(page)).plan).toEqual(beforeRecordDrag.plan);
  await page.evaluate(() => window.__GBDRAW_APP__.redoHistory());
  expect((await artifactSnapshot(page)).plan).toBeNull();
  await page.getByRole('button', { name: 'Toggle layout edit mode' }).click();

  await importSession(page, session, 'one-candidate-retry.gbdraw-session.json');
  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) {
    await page.locator('.drawer-toggle').click();
  }
  const similarityTab = drawer.getByRole('button', { name: 'Similarity groups' });
  await similarityTab.click();
  await drawer.locator('[data-similarity-group-row="og_1"]').click();
  await drawer.getByLabel('Exact reference record and feature').selectOption(reference);
  const committed = await artifactSnapshot(page);
  const priorLayout = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const prior = app.form.linear_track_layout;
    const sanitize=window.DOMPurify.sanitize;window.DOMPurify.sanitize=(...args)=>{if(String(args[0]).includes('<svg')){window.DOMPurify.sanitize=sanitize;throw Object.assign(new Error('PRIVATE_CAUSE_SENTINEL'), { code: 'COMPARISON_IDENTITY', stage: 'result-admission', context: { reason: 'SOURCE_VIEW_CONFLICT' }, cause: new Error('PRIVATE_CAUSE_SENTINEL'), stdout: 'PRIVATE_STDOUT_SENTINEL', secondary: [{code:'CLEANUP_FAILED',stage:'cleanup',message:'PRIVATE_CLEANUP_SENTINEL'}] });}return sanitize(...args);};
    return prior;
  });
  await normalAlign.click();
  await expect(resolvedReview).toBeVisible({ timeout: 180000 });
  const retryError = resolvedReview.locator('[data-similarity-alignment-error]');
  await expect(retryError).toContainText('The features at the ends of comparison links do not match the displayed features.');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentError)).toMatchObject({code:'COMPARISON_IDENTITY',operation:'align',stage:'result-admission',secondary:[{code:'CLEANUP_FAILED',stage:'cleanup'}]});
  await expect(page.getByRole('alert', { name: 'Alignment error' })).toHaveCount(1);
  await expect(retryError).toBeFocused();
  await expect(resolvedReview.locator('[data-alignment-record-key]')).toHaveCount(1);
  const failed = await artifactSnapshot(page);
  expect(failed.results).toEqual(committed.results);
  expect(failed.history.slice(0, 2)).toEqual(committed.history.slice(0, 2));
  await page.evaluate((layout) => {
    window.__GBDRAW_APP__.form.linear_track_layout = layout;
  }, priorLayout);
  await resolvedReview.getByLabel(/^Skip /).check();
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
  expect(saved.version).toBe(CURRENT_SESSION_VERSION);
  expect(saved.renderRequest.schema).toBe(CURRENT_REQUEST_SCHEMA);
  expect(saved.editorState.featureCatalog.schema).toBe(CURRENT_FEATURE_CATALOG_SCHEMA);
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
      orientations: state.drawingList().find((drawing) => drawing.mode === 'linear').linearSeqs.map(({ uid, region_reverse }) => ({
        recordKey: uid, reverseComplement: Boolean(region_reverse)
      })),
      results: window.__GBDRAW_APP__.currentResults().map(({ name, content }) => ({ name, content })),
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
    const dialog = page.getByRole('dialog', { name: 'Review alignment' });
    await expect(dialog).toBeVisible({ timeout: 180000 });
    const unresolved=await page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentDraft.rows.filter(row=>!row.choice).map(row=>row.recordKey));
    for(const key of unresolved)await dialog.locator(`[data-alignment-record-key="${key}"]`).getByRole('radio',{name:/Select .*bp, strand/}).first().check();
    return dialog;
  };
  const geometry = () => page.evaluate(async () => {
    const { getFeatureFillElements } = await import('./js/services/feature-dom.js');
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
    const ribbons = Array.from(svg.querySelectorAll('[data-gbdraw-match-id]'));
    return { offsets, ribbons: ribbons.length,
      targetRibbons: ribbons.filter((element) => element.getAttribute('data-query-record-index') === '1'
        || element.getAttribute('data-subject-record-index') === '1').length };
  });
  const applyReviewed=async()=>{
    const historyBefore=(await recordState()).historyCount;
    await expect(apply).toBeEnabled();
    await apply.click();await page.waitForFunction(()=>!window.__GBDRAW_APP__.similarityAlignmentBusy,null,{timeout:180000});
    const message=await page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentError?.summary||'');
    console.log('S02 Gallery: final batch outcome',message,await page.evaluate(()=>{const e=window.__GBDRAW_APP__.errorLog;return {message:e?.message,stack:e?.stack,summary:e?.summary,details:e?.details};}));
    if(message.includes('Validated directions or reference placement changed')){
      expect((await recordState()).historyCount).toBe(historyBefore);
      await apply.click();await page.waitForFunction(()=>!window.__GBDRAW_APP__.similarityAlignmentBusy,null,{timeout:180000});
    }
  };
  // Align… and Review alignment options… leave the feature popup open, and the
  // popup opened at the top-left corner can cover header buttons such as Undo.
  const closeFeaturePopup = async () => {
    const popup = page.locator('.feature-popup[role="dialog"]');
    await popup.getByRole('button', { name: 'Close feature popup', exact: true }).click();
    await expect(popup).toHaveCount(0);
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
  console.log('S02 Gallery: default Align finished',await page.evaluate(()=>({status:window.__GBDRAW_APP__.similarityAlignmentStatus,error:window.__GBDRAW_APP__.similarityAlignmentError?.summary})));
  await expect(page.locator('[data-similarity-alignment-summary]')).toContainText('0 reversed');
  expect((await recordState()).orientations).toEqual(before.orientations);
  await expect(page.getByRole('dialog', { name: 'Review alignment' })).toHaveCount(0);
  await closeFeaturePopup();
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(recordState).toEqual(before);
  console.log('S02 Gallery: Undo restored baseline');
  const dialog = await review();
  console.log('S02 Gallery: initial review');
  await expect(dialog.locator('[data-similarity-alignment-reference]')).toContainText('Streptomyces lividus CBS 844.73');
  const targetRow = dialog.locator('[data-alignment-record-key="record-2"]');
  await expect(targetRow).toContainText('Streptomyces fradiae');
  await expect(dialog.getByRole('checkbox',{name:'Match reference direction'})).toHaveCount(0);
  const direction = dialog.getByRole('radio', { name: 'All left', exact: true });
  await expect(direction).not.toBeChecked();
  await direction.focus();
  await page.keyboard.press('Space');
  await expect(direction).toBeChecked();
  await expect(directionRow(dialog, 'record-2')).toContainText('reversed');
  await expect(dialog.locator('#alignment-direction-record-2')).toHaveText(/Current right-facing; after Align left-facing\./);
  await dialog.screenshot({ path: testInfo.outputPath('direction-review-desktop.png') });
  await page.setViewportSize({ width: 390, height: 740 });
  // The review is docked in the Result: a phone shows it in the Preview pane.
  await showPhonePane(page, 'preview');
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
    expect(await page.evaluate(box => [25, 100, 175].map(offset => Boolean(
      document.elementFromPoint(box.x + box.width / 2, box.y + offset)?.closest('.preview-canvas')
    )), canvas)).toEqual([true, true, true]);
    for (const control of [apply, dialog.getByRole('button', { name: 'Cancel', exact: true })]) {
      await exposeReviewControl(page, control);
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
    .toContainText('without recognized diagnostic information');
  await expect(page.getByRole('alert', { name: 'Alignment error' })).toHaveCount(1);
  await expect(direction).toBeChecked();
  expect(await recordState()).toEqual(before);
  console.log('S02 Gallery: Apply requested');
  await applyReviewed();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
  const applyStatus = await page.evaluate(() => ({ status: window.__GBDRAW_APP__.similarityAlignmentStatus,
    error: window.__GBDRAW_APP__.similarityAlignmentError?.summary, globalError: window.__GBDRAW_APP__.errorLog }));
  expect(applyStatus, JSON.stringify(applyStatus)).toMatchObject({ status: 'idle', globalError: null });
  await expect(page.getByRole('alert', { name: 'Generation Error' })).toHaveCount(0);
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
      '[data-similarity-alignment-plan-inspector] [data-alignment-record-key]:has([aria-label$=" reversed relative to source"])'
    )].map((element) => element.dataset.alignmentRecordKey);
    app.showRightDrawer = previous.open;
    app.rightDrawerTab = previous.tab;
    await window.Vue.nextTick();
    return labels;
  });
  expect(await revLabels()).toEqual(applied.orientations.filter(({ reverseComplement }) => reverseComplement)
    .map(({ recordKey }) => recordKey));
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
  await closeFeaturePopup();
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(recordState).toEqual(before);
  console.log('S02 Gallery: Undo restored baseline');
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(recordState).toEqual(applied);

  console.log('S02 Gallery: positions Reset');
  const reset = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.resetSimilarityAlignment());
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
  await openLinearFileCards(page);
  await page.getByRole('button', { name: 'Record options for sequence 2', exact: true }).click();
  const reverse = page.getByRole('checkbox', { name: 'Reverse complement for sequence 2', exact: true });
  await expect(reverse).toBeChecked();
  await reverse.uncheck();
  expect((await recordState()).plan).toEqual(applied.plan);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentNotice))
    .not.toContain('Alignment cleared');
  console.log('S02 Gallery: manual Reverse generation');
  const generated = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  expect(generated.status).toBe('ok');
  const manual = await recordState();
  expect(manual.plan).toEqual(applied.plan);
  expect(manual.receipt).toEqual(applied.receipt);
  await page.evaluate(()=>window.__GBDRAW_APP__.openSimilarityAlignmentReset());
  const resetReview=page.getByRole('dialog',{name:'Reset alignment',exact:true});
  // The targets show only for the directions scope, one line per record.
  await expect(resetReview.getByRole('list', { name: 'Alignment direction restoration targets' })).toHaveCount(0);
  await resetReview.getByRole('radio', { name: /^Positions and directions \(\d+ records\)$/ }).check();
  await expect(resetReview.getByRole('button', { name: 'Reset positions and directions', exact: true })).toBeVisible();
  await expect(resetReview).toContainText('changed after Align');
  await expect(resetReview).toContainText('Direction changes you made after Align on these records are replaced.');
  expect(manual.receipt.directions.length).toBeGreaterThan(1);
  const restorationTargets = resetReview.getByRole('list', { name: 'Alignment direction restoration targets' });
  await expect(restorationTargets.getByRole('listitem')).toHaveCount(manual.receipt.directions.length);
  await page.setViewportSize({ width: 390, height: 740 });
  await restorationTargets.getByRole('listitem').last().scrollIntoViewIfNeeded();
  await expect(restorationTargets.getByRole('listitem').last()).toBeVisible();
  expect(await resetReview.evaluate(element => element.scrollWidth <= element.clientWidth)).toBe(true);
  const resetFooter = await resetReview.getByRole('button', { name: 'Reset positions and directions', exact: true }).boundingBox();
  expect(resetFooter.y + resetFooter.height).toBeLessThanOrEqual(740);
  await resetReview.screenshot({ path: testInfo.outputPath('reset-targets-390.png') });
  await page.setViewportSize({ width: 1600, height: 1000 });
  await resetReview.getByRole('button',{name:'Cancel',exact:true}).click();
  expect(await recordState()).toEqual(manual);
  expect(manual.orientations.find(({ recordKey }) => recordKey === 'record-2').reverseComplement).toBe(false);
  expect(await revLabels()).not.toContain('record-2');
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
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.saveSessionWithTitle());
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
    const observed=await evaluateWithRetainedPromise(freshPage,async scope=>{
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
    await evaluateWithRetainedPromise(freshPage,()=>window.__GBDRAW_HISTORY__.undo());
    expect((await recordState(freshPage)).receipt).toEqual(saved.receipt);
  }
  const corrupt=JSON.parse(gunzipSync(readFileSync(savedPath)).toString('utf8'));
  corrupt.editorState.alignmentResetReceipt.binding='0'.repeat(64);
  // A receipt that no longer binds the request is dropped with a row before
  // anything changes (receipt 〔確認〕, BA-D25); Cancel keeps the restored receipt.
  const pendingCorrupt=evaluateWithRetainedPromise(freshPage,async raw=>window.__GBDRAW_APP__.importSession({target:{files:[new File([raw],'corrupted.json')],value:'selected'}}),JSON.stringify(corrupt));
  await answerLoadRepairs(freshPage,['Session file: the field editorState.alignmentResetReceipt could not be read and was left out.'],'cancel');
  expect((await pendingCorrupt).status).toBe('canceled');expect((await recordState(freshPage)).receipt).toEqual(saved.receipt);
  // Re-render the restored Session and compare the complete SVG tree.
  const regenerated = await evaluateWithRetainedPromise(freshPage, () => window.__GBDRAW_APP__.runAnalysis());
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
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.downloadSVG());
  const exported = await exportPromise;
  const exportPath = testInfo.outputPath('record-owned-direction-export.svg');
  await exported.saveAs(exportPath);
  expect(readFileSync(exportPath, 'utf8'))
    .not.toMatch(/gbdraw-alignment-|Review alignment|Apply required|Nothing changes until Apply|Match reference direction/);
  const downloadFormat = async (format, suffix) => {
    const pending = page.waitForEvent('download');
    await evaluateWithRetainedPromise(page, (method) => window.__GBDRAW_APP__[method](), `download${format}`);
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
  // Reference-only Custom Align C preserves pending form and anchor-center x.
  const beforeC=await recordState(freshPage);
  const referenceCenter=()=>freshPage.evaluate(async anchor=>{
    const {getFeatureFillElements}=await import('./js/services/feature-dom.js');
    const app=window.__GBDRAW_APP__;
    const feature=app.extractedFeatures.find(item=>item.recordKey===anchor.recordKey && item.biologicalFeatureId===anchor.biologicalFeatureId);
    const svg=document.querySelector('.gbdraw-preview-surface svg');
    const bounds=getFeatureFillElements(svg,feature.svg_id).map(element=>element.getBoundingClientRect());
    const screenX=(Math.min(...bounds.map(box=>box.left))+Math.max(...bounds.map(box=>box.right)))/2;
    const screenY=(Math.min(...bounds.map(box=>box.top))+Math.max(...bounds.map(box=>box.bottom)))/2;
    // S01 placement facts exclude automatic composition/viewBox fitting.
    // Remove that admitted global placement as well as viewport pan/zoom.
    const {parseCompositionMetadata}=await import('./js/app/legend-layout/composition-actions.js');
    const automaticX=parseCompositionMetadata(svg).primary.automaticTranslation[0];
    return new DOMPoint(screenX,screenY).matrixTransform(svg.getScreenCTM().inverse()).x-automaticX;
  },saved.plan.reference);
  const centerBeforeC=await referenceCenter();
  const referenceFeature=await freshPage.evaluate(anchor=>window.__GBDRAW_APP__.extractedFeatures
    .find(item=>item.recordKey===anchor.recordKey && item.biologicalFeatureId===anchor.biologicalFeatureId).svg_id,saved.plan.reference);
  await freshPage.locator(`[data-gbdraw-feature-id="${referenceFeature}"]`).first().dispatchEvent('click',{clientX:-100,clientY:-100});
  await freshPage.locator('.feature-popup[role="dialog"]').getByRole('button',{name:'Review alignment options…'}).click();
  const customReview=freshPage.getByRole('dialog',{name:'Review alignment'});
  await expect(customReview).toBeVisible({timeout:180000});
  const unresolvedC=await freshPage.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentDraft.rows.filter(row=>!row.choice).map(row=>row.recordKey));
  for(const key of unresolvedC)await customReview.locator(`[data-alignment-record-key="${key}"]`).getByRole('radio',{name:/Select .*bp, strand/}).first().check();
  await freshPage.evaluate(()=>{
    window.__GBDRAW_APP__.form.prefix='pending-unapplied-prefix';
    window.__S02FinalBatches=0;
    const postMessage=Worker.prototype.postMessage;
    Worker.prototype.postMessage=function(message,transfer){
      if(message.type==='helper' && message.operation==='resolveSimilarityAlignment')window.__S02FinalBatches++;
      return postMessage.call(this,message,transfer);
    };
  });
  console.log('S02 Gallery: Custom reference review ready');
  await customReview.getByRole('radio',{name:'Custom',exact:true}).check();
  await directionSelect(customReview, saved.plan.reference.recordKey).selectOption('right');
  expect(await freshPage.evaluate(()=>window.__S02FinalBatches)).toBe(0);
  const referencePreview=await freshPage.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentDirectionPreview.reference);
  expect(Math.abs(referencePreview.beforeX-centerBeforeC)).toBeLessThanOrEqual(0.5);
  await customReview.getByRole('button',{name:'Apply',exact:true}).click();
  await freshPage.waitForFunction(()=>!window.__GBDRAW_APP__.similarityAlignmentBusy,null,{timeout:180000});
  expect(await freshPage.evaluate(()=>window.__S02FinalBatches)).toBe(1);
  if(await freshPage.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentError?.summary.includes('Validated directions or reference placement changed'))) {
    expect((await recordState(freshPage)).historyCount).toBe(beforeC.historyCount);
    await customReview.getByRole('button',{name:'Apply',exact:true}).click();
    await freshPage.waitForFunction(()=>!window.__GBDRAW_APP__.similarityAlignmentBusy,null,{timeout:180000});
    expect(await freshPage.evaluate(()=>window.__S02FinalBatches)).toBe(2);
  }
  await expect(customReview).toBeHidden({timeout:180000});
  const afterC=await recordState(freshPage);
  expect(afterC.historyCount).toBe(beforeC.historyCount+1);
  expect(afterC.receipt.directions).toEqual([{recordKey:saved.plan.reference.recordKey,before:false,after:true}]);
  expect(afterC.receipt.referenceDeltaX.recordKey).toBe(saved.plan.reference.recordKey);
  const centerAfterC=await referenceCenter();
  console.log('S02 reference logical canvas centers',JSON.stringify({before:centerBeforeC,after:centerAfterC,receipt:afterC.receipt.referenceDeltaX}));
  expect(Math.abs(centerAfterC-centerBeforeC)).toBeLessThanOrEqual(0.5);
  expect(await freshPage.evaluate(()=>window.__GBDRAW_APP__.form.prefix)).toBe('pending-unapplied-prefix');
  expect((await freshPage.evaluate(async()=>{const {getCommittedCanonicalSession}=await import('./js/services/config.js');return getCommittedCanonicalSession().renderRequest.output.prefix;})))
    .not.toBe('pending-unapplied-prefix');
  expect((await evaluateWithRetainedPromise(freshPage,()=>window.__GBDRAW_APP__.resetSimilarityAlignment('positions-and-directions'))).status).toBe('ok');
  expect((await recordState(freshPage)).orientations).toEqual(beforeC.orientations);
  expect(await freshPage.evaluate(()=>window.__GBDRAW_APP__.form.prefix)).toBe('pending-unapplied-prefix');
  await evaluateWithRetainedPromise(freshPage,()=>window.__GBDRAW_HISTORY__.undo());
  expect((await recordState(freshPage)).receipt).toEqual(afterC.receipt);
  await evaluateWithRetainedPromise(freshPage,()=>window.__GBDRAW_HISTORY__.redo());
  expect((await recordState(freshPage)).receipt).toBeNull();
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
  const dialog = page.getByRole('dialog', { name: 'Review alignment' });
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
  await expect(dialog.getByRole('radio', { name: 'Keep as is', exact: true })).toBeChecked();
  await expect(dialog.getByRole('combobox')).toHaveCount(0);
  const right = dialog.getByRole('radio', { name: 'All right', exact: true });
  await right.focus();
  await page.keyboard.press('Space');
  expect(await changed()).toEqual(['record-1']);
  // D-11: reversing the reference record keeps its feature in place; the card
  // says how far the record's left edge moves.
  await expect(dialog.locator('[data-similarity-alignment-reference-correction]'))
    .toContainText(/stays where it is; the record's left edge moves [1-9][\d,]* px \(at 100% zoom\) to the (left|right)\./);
  await dialog.getByRole('radio', { name: 'All left', exact: true }).check();
  expect(await changed()).toEqual(['record-2', 'record-3', 'record-4', 'record-5']);
  const custom = dialog.getByRole('radio', { name: 'Custom', exact: true });
  await custom.focus();
  await page.keyboard.press('Space');
  await expect(dialog.getByRole('combobox')).toHaveCount(5);
  expect(await changed()).toEqual([]);
  await directionSelect(dialog, 'record-1').selectOption('right');
  await directionSelect(dialog, 'record-2').selectOption('left');
  expect(await changed()).toEqual(['record-1', 'record-2']);
  const target = dialog.locator('[data-alignment-record-key="record-2"]');
  await target.getByLabel(/^Skip /).focus();
  await page.keyboard.press('Space');
  expect(await changed()).toEqual(['record-1']);
  await expect(directionSelect(dialog, 'record-2')).toBeDisabled();
  await expect(directionRow(dialog, 'record-2'))
    .toContainText('Unchanged: Skipped by user.');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDirectionPreview.records.find(record => record.recordKey === 'record-2').exclusion)).toBe('skipped_by_user');
  await target.getByRole('checkbox', { name: /^Skip / }).uncheck();
  expect(await changed()).toEqual(['record-1', 'record-2']);
  await page.setViewportSize({ width: 390, height: 740 });
  // At <= 768 px one pane shows; the review is in the Result (Preview pane).
  await showPhonePane(page, 'preview');
  await expect(custom).toBeChecked();
  await expect(directionSelect(dialog, 'record-1')).toHaveValue('right');
  await expect(directionSelect(dialog, 'record-2')).toHaveValue('left');
  await page.setViewportSize({ width: 1600, height: 1000 });
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
    await expect(dialog.locator('[data-similarity-alignment-error]')).toContainText('Review the updated preview and Apply again');
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
  await drawer.locator('[data-similarity-group-row="group"]').click();
  const key = await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDrawerReferenceOptions('group')
    .find(({ anchor }) => anchor.recordKey === 'record-1').key);
  await drawer.getByLabel('Exact reference record and feature').selectOption(key);
  await drawer.getByRole('button', { name: 'Review alignment options…' }).click();
  const dialog = page.getByRole('dialog', { name: 'Review alignment' });
  await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentError?.summary || '')).toBe('');
  await expect(dialog).toBeVisible({ timeout: 180000 });
  return dialog;
};

const applyAndSettle = async (page, dialog) => {
  await dialog.getByRole('button', { name: 'Apply', exact: true }).press('Enter');
  await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
  return page.evaluate(() => ({ status: window.__GBDRAW_APP__.similarityAlignmentStatus,
    error: window.__GBDRAW_APP__.similarityAlignmentError?.summary || '' }));
};

test('minority reference, source strands and Reset scopes have truthful accessible previews', async ({ page }, testInfo) => {
  test.setTimeout(600000);
  await loadDirectionSession(page, testInfo);
  const dialog = await openDirectionReview(page);
  await dialog.getByLabel(/^Skip skip/).press('Space');
  const unknown = dialog.locator('#alignment-direction-record-3');
  const missing = dialog.locator('#alignment-direction-record-5');
  await dialog.getByRole('radio', { name: 'All right', exact: true }).check();
  await expect(unknown).toContainText('no known direction; it is never guessed');
  await expect(missing).toContainText('No usable candidate');
  await expect(dialog.locator('#alignment-direction-scope')).toContainText('all 2 chosen features point right');
  const initial = await artifactSnapshot(page);
  const preview = await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDirectionPreview);
  expect(preview.records.filter(r => r.beforeReverseComplement !== r.afterReverseComplement).map(r => r.recordKey)).toEqual(['record-1']);
  expect(preview.reference.deltaX).not.toBe(0);
  await dialog.screenshot({ path: testInfo.outputPath('minority-reference-desktop.png') });
  let outcome = await applyAndSettle(page, dialog);
  if (outcome.error.includes('Review the updated preview and Apply again')) outcome = await applyAndSettle(page, dialog);
  expect(outcome.error).toBe('');
  await expect(dialog).toBeHidden();
  const aligned = await artifactSnapshot(page);
  expect(aligned.receipt.directions).toEqual([{ recordKey: 'record-1', before: false, after: true }]);
  expect(aligned.session.resources).toEqual(initial.session.resources);
  const resetButton = page.locator('[data-similarity-alignment-reset]');
  await resetButton.click();
  const reset = page.getByRole('dialog', { name: 'Reset alignment', exact: true });
  const positions = reset.getByRole('radio', { name: 'Positions only', exact: true });
  const combined = reset.getByRole('radio', { name: /^Positions and directions/ });
  const resetPositions = reset.getByRole('button', { name: 'Reset positions', exact: true });
  const resetBoth = reset.getByRole('button', { name: 'Reset positions and directions', exact: true });
  await expect(reset.getByRole('heading')).toBeFocused();
  await expect(reset).toHaveAccessibleDescription('Undo the latest Align.');
  await expect(positions).toBeChecked();
  await expect(combined).toHaveAccessibleName('Positions and directions (1 record)');
  // UX-A01: the targets and their disclosure appear with the directions scope only.
  await expect(reset.getByRole('list', { name: 'Alignment direction restoration targets' })).toHaveCount(0);
  await expect(resetPositions).toBeVisible();
  await positions.focus(); await page.keyboard.press('ArrowDown');
  await expect(combined).toBeChecked(); await expect(positions).not.toBeChecked();
  await expect(reset.getByRole('list', { name: 'Alignment direction restoration targets' })).toContainText('ref');
  await expect(reset.getByRole('list')).toContainText('reverse → forward');
  await expect(reset).toContainText('Direction changes you made after Align on these records are replaced.');
  await expect(resetBoth).toBeVisible(); await expect(resetPositions).toHaveCount(0);
  for (const height of [740, 844]) {
    await page.setViewportSize({ width: 390, height });
    await expect(resetBoth).toBeVisible();
    await reset.getByRole('list').scrollIntoViewIfNeeded();
    expect(await reset.evaluate(e => e.scrollWidth <= e.clientWidth)).toBe(true);
    await reset.screenshot({ path: testInfo.outputPath(`reset-390-${height}.png`) });
  }
  await combined.press('Escape'); await expect(reset).toBeHidden(); await expect(resetButton).toBeFocused();
  await resetButton.press('Enter'); await expect(positions).toBeChecked();
  await resetPositions.press('Enter');
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
  await resetBoth.press('Enter');
  await expect(reset).toBeHidden({ timeout: 180000 });
  const restored = await artifactSnapshot(page);
  expect(restored.session.renderRequest.records).toEqual(initial.session.renderRequest.records);
  expect(restored.plan).toBeNull(); expect(restored.receipt).toBeNull();
  await page.setViewportSize({ width: 1600, height: 1000 });
  const targetOnly = await openDirectionReview(page);
  await targetOnly.getByLabel(/^Skip skip/).check();
  await targetOnly.getByRole('radio', { name: 'All left', exact: true }).check();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDirectionPreview.records
    .filter(r => r.beforeReverseComplement !== r.afterReverseComplement).map(r => r.recordKey))).toEqual(['record-2']);
  outcome = await applyAndSettle(page, targetOnly);
  if (outcome.error.includes('Review the updated preview and Apply again')) outcome = await applyAndSettle(page, targetOnly);
  expect(outcome.error).toBe('');
  await expect(targetOnly).toBeHidden();
  expect((await artifactSnapshot(page)).receipt.directions).toEqual([{ recordKey: 'record-2', before: false, after: true }]);
  const keepReview = await openDirectionReview(page);
  await keepReview.getByLabel(/^Skip skip/).check();
  expect((await applyAndSettle(page, keepReview)).error).toBe('');
  await expect(keepReview).toBeHidden();
  expect((await artifactSnapshot(page)).receipt.directions).toEqual([]);
  await resetButton.click();
  await expect(combined).toBeDisabled(); await expect(positions).toBeEnabled();
  await expect(combined).toHaveAccessibleName('Positions and directions');
  await expect(combined).toHaveAccessibleDescription('Not available: the latest Align did not change any record direction.');
  await reset.getByRole('button', { name: 'Cancel', exact: true }).press('Enter');
  // A Session with an Align plan but no Reset receipt (PD-OI-028); the Gallery Session ships an empty one.
  const historical = JSON.parse(readFileSync('gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json', 'utf8'));
  historical.editorState.alignmentResetReceipt = null;
  await importSession(page, Buffer.from(JSON.stringify(historical)), 'historical-receipt.json');
  if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) await page.locator('.drawer-toggle').click();
  await page.locator('.right-drawer').getByRole('button', { name: 'Similarity groups' }).click();
  await resetButton.click();
  await expect(combined).toBeDisabled(); await expect(positions).toBeEnabled();
  await expect(combined).toHaveAccessibleDescription(/does not record which directions the latest Align changed/);
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
  await dialog.getByLabel(/^Skip skip/).check();
  await dialog.getByRole('radio', { name: 'Custom', exact: true }).check();
  await directionSelect(dialog, 'record-1').selectOption('right');
  const before = await artifactSnapshot(page);
  const oldDelta = await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDirectionPreview.reference.deltaX);
  await applyAndSettle(page, dialog);
  await expect(dialog.locator('[data-similarity-alignment-error]')).toContainText('Review the updated preview and Apply again');
  expect(await artifactSnapshot(page)).toEqual(before);
  expect(await page.evaluate(() => window.__s03FinalBatches)).toBe(1);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDirectionPreview.reference.deltaX)).not.toBe(oldDelta);
  await page.evaluate(() => { window.__s03ValidationFailure = true; });
  await applyAndSettle(page, dialog);
  await expect(dialog.locator('[data-similarity-alignment-error]')).toContainText('without recognized diagnostic information');
  await expect(directionSelect(dialog, 'record-1')).toHaveValue('right');
  expect(await artifactSnapshot(page)).toEqual(before);
  expect(await page.evaluate(() => window.__s03FinalBatches)).toBe(2);
  await applyAndSettle(page, dialog); await expect(dialog).toBeHidden();
  expect(await page.evaluate(() => window.__s03FinalBatches)).toBe(3);
  expect((await artifactSnapshot(page)).history[0]).toBe(before.history[0] + 1);
});
// S06 acceptance uses the existing controller, real Gallery sources and native input.
const exposeReviewControl = async (page, control) => {
  await control.evaluate(element => {
    element.scrollIntoView({block:'center',inline:'nearest'});
    const box=element.getBoundingClientRect();
    const header=document.querySelector('.app-header').getBoundingClientRect();
    const bar=document.querySelector('.generate-bar').getBoundingClientRect();
    const bottom=bar.left<box.right&&bar.right>box.left?bar.top:innerHeight;
    if(box.top<header.bottom)window.scrollBy(0,box.top-header.bottom-8);
    else if(box.bottom>bottom)window.scrollBy(0,box.bottom-bottom+8);
  });
  const pointerState=()=>control.evaluate(element=>{
    const box=element.getBoundingClientRect();
    const hit=document.elementFromPoint(box.x+box.width/2,box.y+box.height/2);
    const rect=el=>el?el.getBoundingClientRect().toJSON():null;
    const list=element.closest('.custom-scrollbar');
    return {reachable:hit===element||element.contains(hit),control:element.outerHTML,box:rect(element),hit:hit?.outerHTML.slice(0,800),header:rect(document.querySelector('.app-header')),bar:rect(document.querySelector('.generate-bar')),list:rect(list),scrollTop:list?.scrollTop,scrollHeight:list?.scrollHeight,scrollY,viewport:[innerWidth,innerHeight]};
  });
  try {
    await expect.poll(async()=>(await pointerState()).reachable,{message:'the actual operation center must receive pointer input'}).toBe(true);
  } catch(error) {
    console.log('S06_POINTER_FAILURE',JSON.stringify(await pointerState()));
    throw error;
  }
};
const compactReviewGeometry = (page) => page.evaluate(()=>{
  const canvas=document.querySelector('.preview-canvas');
  const palette=document.querySelector('[data-similarity-alignment-dialog]');
  const rect=el=>{const b=el.getBoundingClientRect();return {x:b.x,y:b.y,width:b.width,height:b.height,left:b.left,right:b.right,top:b.top,bottom:b.bottom}};
  const c=rect(canvas),p=rect(palette),header=rect(document.querySelector('.app-header')),bar=rect(document.querySelector('.generate-bar'));
  const v={left:Math.max(0,c.left),right:Math.min(innerWidth,c.right),top:Math.max(0,c.top,header.bottom),bottom:Math.min(innerHeight,c.bottom,bar.left<c.right&&bar.right>c.left?bar.top:innerHeight)};
  const overlap=Math.max(0,Math.min(p.right,v.right)-Math.max(p.left,v.left))*Math.max(0,Math.min(p.bottom,v.bottom)-Math.max(p.top,v.top));
  const list=palette.querySelector('.custom-scrollbar');
  return {viewport:[innerWidth,innerHeight],canvas:c,review:p,visibleWidth:Math.max(0,v.right-v.left),visibleHeight:Math.max(0,v.bottom-v.top),overlap,
    pageWidth:document.documentElement.scrollWidth,listScrollable:list.scrollHeight>list.clientHeight,
    editorOpen:window.__GBDRAW_APP__.showRightDrawer,editorTab:window.__GBDRAW_APP__.rightDrawerTab};
});

for(const entry of ['resolved explicit','multiple target ambiguity','automatic render retry']) {
  test(`compact review preserves canvas, local intent and recovery: ${entry}`,async({page},info)=>{
    test.setTimeout(600000);
    page.setDefaultTimeout(30000);
    await captureUnhandledRejections(page);
    const pageErrors=[];page.on('pageerror',error=>pageErrors.push(String(error)));
    await page.addInitScript(()=>{
      window.__S06_POSTS__=[];window.__GBDRAW_TEST_HOOKS__={};
      const original=Worker.prototype.postMessage;
      Worker.prototype.postMessage=function(message,...rest){window.__S06_POSTS__.push({type:message.type,operation:message.operation});return original.call(this,message,...rest)};
    });
    await page.setViewportSize({width:1600,height:1000});
    await page.goto('/gbdraw/web/index.html',{waitUntil:'domcontentloaded'});
    await page.waitForFunction(()=>window.__GBDRAW_APP__);
    let seedPath='gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json';
    if(entry==='multiple target ambiguity') {
      // Direct evidence comes from the committed resource; drop og_1's edges there.
      const stripped=info.outputPath('bgc-og1-without-edges.json');
      execFileSync('python',['-c',[
        'import base64, json, sys', 'from dataclasses import replace',
        'from gbdraw.analysis.protein_colinearity import OrthogroupGraphResult, OrthogroupResult',
        'from gbdraw.session_request_codec import decode_canonical_typed_resource, encode_canonical_typed_resource',
        'session = json.load(open(sys.argv[1]))', "entry = session['resources']['comparison-orthogroups']",
        "result = decode_canonical_typed_resource(base64.b64decode(entry['data']), value_kind='orthogroupResult', expected=OrthogroupResult | OrthogroupGraphResult)",
        "drop = lambda values: {key: value for key, value in values.items() if key != 'og_1'}",
        "data = encode_canonical_typed_resource('orthogroupResult', replace(result, ortholog_edges_by_orthogroup_id=drop(result.ortholog_edges_by_orthogroup_id), path_indexes_by_orthogroup_id=drop(result.path_indexes_by_orthogroup_id)))",
        "entry.update(data=base64.b64encode(data).decode(), size=len(data))", "json.dump(session, open(sys.argv[2], 'w'))"
      ].join('\n'),seedPath,stripped],{cwd:process.cwd(),stdio:'pipe'});
      seedPath=stripped;
    }
    const seed=JSON.parse(readFileSync(seedPath,'utf8'));
    if(entry==='multiple target ambiguity') {
      const catalog=seed.editorState.featureCatalog.items[0];
      const group=catalog.orthogroups.find(group=>group.id==='og_1');
      for(const recordKey of ['record-2','record-3']) {
        const existing=new Set(group.members.filter(member=>member.recordKey===recordKey).map(member=>member.biologicalFeatureId));
        const feature=catalog.biologicalFeatures.find(feature=>feature.recordKey===recordKey&&feature.type==='CDS'&&!existing.has(feature.biologicalFeatureId));
        expect(feature).toBeTruthy();
        group.members.push({recordKey,biologicalFeatureId:feature.biologicalFeatureId,representative:false,role:'member',confidence:'low',assignmentReason:'disposable ambiguity fixture'});
      }
      group.member_count=group.members.length;
    }
    await importSession(page,Buffer.from(JSON.stringify(seed)),'s06-gallery.gbdraw-session.json');
    await page.evaluate(async()=>{const app=window.__GBDRAW_APP__;app.openRightDrawerTab('orthogroups');app.selectedOrthogroupId='og_1';app.losatExecution.executionMode='serial';await window.Vue.nextTick()});
    const drawer=page.locator('.right-drawer');
    const reference=drawer.getByLabel('Exact reference record and feature');
    if(!await reference.isVisible())await drawer.locator('[data-similarity-group-row="og_1"]').click();
    const key=await page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentDrawerReferenceOptions('og_1').find(option=>option.anchor.recordKey==='record-1').key);
    await reference.selectOption(key);
    let before=await artifactSnapshot(page);
    if(entry==='automatic render retry')await page.evaluate(()=>{window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse=()=>{throw new Error('S06 forced automatic render failure.')}});
    await page.setViewportSize({width:390,height:844});
    const start=drawer.getByRole('button',{name:entry==='resolved explicit'?'Review alignment options…':'Align…',exact:true});
    await exposeReviewControl(page,start);await start.click();
    const dialog=page.getByRole('dialog',{name:'Review alignment'});
    await expect(dialog).toBeVisible({timeout:180000});
    await expect.poll(()=>page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentBusy),{timeout:180000}).toBe(false);
    if(entry==='multiple target ambiguity')expect(await page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentDraft.rows.filter(row=>row.candidates.length>1).length)).toBeGreaterThanOrEqual(2);
    if(entry==='automatic render retry') {
      await expect(dialog.locator('[data-similarity-alignment-error]')).toContainText('without recognized diagnostic information');
      // Aborting the existing artifact transaction increments an internal
      // invalidation revision; the previous artifact and both History stacks
      // remain unchanged. Resize and local review must retain this settled state.
      const settled=await artifactSnapshot(page);
      expect({...settled,history:[...settled.history.slice(0,2),before.history[2]]}).toEqual(before);
      before=settled;
      for(const viewport of [{width:390,height:740},{width:320,height:740},{width:195,height:422}]) {
        await page.setViewportSize(viewport);
        await exposeReviewControl(page,dialog.locator('[data-alignment-record-key]').first().locator('input').first());
        await exposeReviewControl(page,dialog.getByRole('button',{name:'Cancel',exact:true}));
        await exposeReviewControl(page,dialog.locator('[data-similarity-alignment-error]'));
        await expect(dialog.locator('[data-similarity-alignment-error]')).toContainText('without recognized diagnostic information');
        expect(await artifactSnapshot(page)).toEqual(before);
      }
      await page.setViewportSize({width:390,height:844});
    }
    expect(await artifactSnapshot(page)).toEqual(before);
    await expect(dialog).not.toHaveAttribute('aria-modal','true');
    await expect(page.locator('.drawer-toggle')).toBeDisabled();
    await expect(page.locator('.drawer-toggle')).toHaveAttribute('title',/Finish or cancel alignment review/);
    expect(await page.evaluate(()=>({open:window.__GBDRAW_APP__.showRightDrawer,tab:window.__GBDRAW_APP__.rightDrawerTab}))).toEqual({open:false,tab:'orthogroups'});
    await dialog.getByRole('radio',{name:'Custom',exact:true}).check();
    const target=dialog.locator('[data-alignment-record-key]').first();
    const select=target.getByRole('radio',{name:/Select .*bp, strand/}).first();
    const skip=target.getByLabel(/^Skip /);
    if(await select.count()) {
      await exposeReviewControl(page,select);await select.focus();await page.keyboard.press('Space');
      await exposeReviewControl(page,skip);await skip.focus();await page.keyboard.press('Space');await expect(skip).toBeChecked();
      await exposeReviewControl(page,select);await select.focus();await page.keyboard.press('Space');await expect(select).toBeChecked();
    } else {
      // A one-candidate record: Space on its Skip checkbox skips, and again selects.
      await exposeReviewControl(page,skip);await skip.focus();await page.keyboard.press('Space');await expect(skip).toBeChecked();
      await page.keyboard.press('Space');await expect(skip).not.toBeChecked();
    }
    const custom=directionSelect(dialog, 'record-1');await custom.selectOption('right');
    const local=await page.evaluate(()=>JSON.stringify(window.__GBDRAW_APP__.similarityAlignmentDraft));
    const posts=await page.evaluate(()=>window.__S06_POSTS__.length);
    const geometry=[];
    for(const viewport of [{width:390,height:844},{width:390,height:740},{width:390,height:500},{width:320,height:740},{width:844,height:390},{width:195,height:422},{width:390,height:422},{width:1600,height:1000}]) {
      await page.setViewportSize(viewport);
      if(viewport.width===390&&viewport.height===422)await page.getByRole('searchbox',{name:'Search features',exact:true}).focus();
      if(await select.count())await exposeReviewControl(page,select);
      await exposeReviewControl(page,skip);await exposeReviewControl(page,custom);
      for(const control of [dialog.getByRole('button',{name:'Apply',exact:true}),dialog.getByRole('button',{name:'Cancel',exact:true})])await exposeReviewControl(page,control);
      await page.locator('.preview-canvas').evaluate(el=>{
        el.scrollIntoView({block:'start'});
        // scrollIntoView keeps the canvas's 9rem scroll margin, which in a short phone Preview pane
        // (390x422: 111 px) puts the canvas below the pane's visible end; a user swipes the pane instead.
        const pane=el.closest('.result-pane').getBoundingClientRect();
        if(el.getBoundingClientRect().top>=pane.bottom)el.closest('.result-pane').scrollTop+=el.getBoundingClientRect().top-Math.max(pane.top,document.querySelector('.drawing-bar').getBoundingClientRect().bottom);
      });
      const measured=await compactReviewGeometry(page);geometry.push(measured);
      await info.attach(`review-${viewport.width}x${viewport.height}`,{body:JSON.stringify(measured,null,2),contentType:'application/json'});
      expect(await dialog.evaluate(el=>el.scrollWidth<=el.clientWidth)).toBe(true);
      if(viewport.width===390&&[740,844].includes(viewport.height)) {
        expect(measured.visibleHeight).toBeGreaterThanOrEqual(200);expect(measured.visibleWidth).toBeGreaterThanOrEqual(369);
        expect(measured.overlap).toBe(0);expect(measured.listScrollable).toBe(true);
      }
      // Docked beside the figure or below it, the panel never covers the figure.
      expect(measured.overlap).toBe(0);
      expect(measured.editorOpen).toBe(false);await expect(page.locator('.drawer-toggle')).toBeDisabled();
      // Native pointer input, not a DOM dispatch or a fit-camera assignment.
      const canvas=await page.locator('.preview-canvas').boundingBox();
      const pan=await page.evaluate(()=>({...window.__GBDRAW_APP__.canvasPan}));
      const point={x:canvas.x+5,y:Math.max(canvas.y+5,await page.locator('.app-header').evaluate(el=>el.getBoundingClientRect().bottom+5),await paneTop(page)+5)};
      await page.mouse.move(point.x,point.y);await page.mouse.down();await page.mouse.move(point.x+24,point.y+10,{steps:4});await page.mouse.up();
      await expect.poll(()=>page.evaluate(()=>({...window.__GBDRAW_APP__.canvasPan}))).toEqual({x:pan.x+24,y:pan.y+10});
      const zoom=page.getByRole('button',{name:'Zoom in',exact:true});await exposeReviewControl(page,zoom);const zoomBefore=await page.evaluate(()=>window.__GBDRAW_APP__.zoom);await zoom.click();
      await expect.poll(()=>page.evaluate(()=>window.__GBDRAW_APP__.zoom)).toBeGreaterThan(zoomBefore);
      for(const control of await page.locator('.preview-controls button').all())await exposeReviewControl(page,control);
      expect(await artifactSnapshot(page)).toEqual(before);
      expect(await page.evaluate(()=>JSON.stringify(window.__GBDRAW_APP__.similarityAlignmentDraft))).toBe(local);
      expect(await page.evaluate(()=>window.__S06_POSTS__.length)).toBe(posts);
      // Record readable full-Gallery QA using only the existing camera.
      await page.locator('.preview-canvas').evaluate(el=>el.scrollIntoView({block:'start'}));
      await page.evaluate(async()=>{const app=window.__GBDRAW_APP__,c=document.querySelector('.preview-canvas').getBoundingClientRect(),s=document.querySelector('.gbdraw-preview-surface svg').getBoundingClientRect();app.zoom=Math.min((c.width-16)/(s.width/app.zoom),(c.height-16)/(s.height/app.zoom));await window.Vue.nextTick()});
      await page.waitForTimeout(230);
      await page.evaluate(async()=>{const app=window.__GBDRAW_APP__,c=document.querySelector('.preview-canvas').getBoundingClientRect(),s=document.querySelector('.gbdraw-preview-surface svg').getBoundingClientRect();app.canvasPan.x+=(c.left+c.right-s.left-s.right)/2;app.canvasPan.y+=(c.top+c.bottom-s.top-s.bottom)/2;await window.Vue.nextTick()});
      await page.waitForTimeout(230); // Existing camera transform transition.
      await page.screenshot({path:info.outputPath(`${entry.replaceAll(' ','-')}-${viewport.width}x${viewport.height}.png`)});
    }
    await info.attach('compact-review-geometry',{body:JSON.stringify(geometry,null,2),contentType:'application/json'});
    const apply=dialog.getByRole('button',{name:'Apply',exact:true});await apply.focus();await page.keyboard.press('Tab');
    expect(await dialog.evaluate(el=>el.contains(document.activeElement))).toBe(false);
    await page.setViewportSize({width:390,height:740});
    if(entry==='automatic render retry') {
      await page.evaluate(()=>{delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse});
      for(const key of await page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentDraft.rows.filter(row=>!row.choice).map(row=>row.recordKey))) {
        const choice=dialog.locator(`[data-alignment-record-key="${key}"]`).getByRole('radio',{name:/Select .*bp, strand/}).first();
        await exposeReviewControl(page,choice);await choice.check();
      }
      for(let attempt=0;attempt<2;attempt++) {
        const batches=await page.evaluate(()=>window.__S06_POSTS__.filter(post=>post.type==='helper').length);
        await exposeReviewControl(page,apply);await apply.click();await expect.poll(()=>page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentBusy),{timeout:180000}).toBe(false);
        expect(await page.evaluate(()=>window.__S06_POSTS__.filter(post=>post.type==='helper').length)).toBe(batches+1);
        if(!await dialog.isVisible())break;
        await expect(dialog.locator('[data-similarity-alignment-error]')).toContainText('Review the updated preview and Apply again');
      }
      await expect(dialog).toBeHidden();expect((await artifactSnapshot(page)).history[0]).toBe(before.history[0]+1);
      expect(await page.evaluate(async()=>{const {state}=await import('./js/state.js');return {error:state.errorLog.value,preserved:state.failedGeneratePreservedResult.value}})).toEqual({error:null,preserved:false});
    } else {
      await exposeReviewControl(page,dialog.getByRole('button',{name:'Cancel',exact:true}));
      await dialog.getByRole('button',{name:'Cancel',exact:true}).focus();await page.keyboard.press('Enter');await expect(dialog).toBeHidden();
      expect(await artifactSnapshot(page)).toEqual(before);
    }
    // The review was opened from the Editor: it is open again on its tab, with the
    // focus on the button that opened the review.
    await expect(page.locator('.drawer-toggle')).toBeEnabled();
    expect(await editorDrawer(page)).toEqual({open:true,tab:'orthogroups'});
    await expect(drawer).toHaveAttribute('aria-hidden','false');
    await expect(start).toBeFocused();
    expect(pageErrors).toEqual([]);await expectNoUnhandledRejections(page);
  });
}

test('@pr-smoke native Python engine error reaches Align review and retains choices for Apply retry', async ({page}, info) => {
  test.setTimeout(300000);
  const { failNextNativeRender, inspectSafeDetails } = require('./helpers/operation-error.cjs');
  const logs=[];page.on('console', message=>logs.push(message.text()));page.on('pageerror',error=>logs.push(error.message));
  await page.setViewportSize({width:1600,height:1000});
  await page.goto('/gbdraw/web/index.html',{waitUntil:'domcontentloaded'});
  await page.waitForFunction(()=>window.__GBDRAW_APP__);
  await importSession(page,readFileSync('gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'),'native-align.json');
  await page.evaluate(()=>{const a=window.__GBDRAW_APP__;a.openRightDrawerTab('orthogroups');a.selectedOrthogroupId='og_1';a.losatExecution.executionMode='serial';});
  const drawer=page.locator('.right-drawer');
  const key=await page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentDrawerReferenceOptions('og_1').find(o=>o.anchor.recordKey==='record-1').key);
  await drawer.getByLabel('Exact reference record and feature').selectOption(key);
  await drawer.getByRole('button',{name:'Review alignment options…',exact:true}).click();
  const dialog=page.getByRole('dialog',{name:'Review alignment'});
  await expect(dialog).toBeVisible({timeout:180000});
  await page.waitForFunction(()=>!window.__GBDRAW_APP__.similarityAlignmentBusy);
  for(const key of await page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentDraft.rows.filter(r=>!r.choice).map(r=>r.recordKey))) {
    await dialog.locator(`[data-alignment-record-key="${key}"]`).getByRole('radio',{name:/Select .*bp, strand/}).first().check();
  }
  const before=await artifactSnapshot(page);
  const choices=await page.evaluate(()=>JSON.stringify(window.__GBDRAW_APP__.similarityAlignmentDraft));
  await failNextNativeRender(page);
  await applyAndSettle(page,dialog);
  const alert=dialog.getByRole('alert',{name:'Alignment error'});
  await expect(alert).toContainText('Not a valid regular expression');
  expect(await page.evaluate(()=>window.__GBDRAW_APP__.similarityAlignmentError)).toMatchObject({code:'REGEX_SYNTAX',operation:'align',stage:'rule-validation',context:{position:1,positionUnit:'python-character'}});
  const after=await artifactSnapshot(page);
  // Observer invalidation can advance revision; both History stacks and every
  // artifact/request/orientation field must remain unchanged.
  expect({...after,history:[...after.history.slice(0,2),before.history[2]]}).toEqual(before);
  expect(await page.evaluate(()=>JSON.stringify(window.__GBDRAW_APP__.similarityAlignmentDraft))).toBe(choices);
  await expect(alert).toBeFocused();
  await expect(page.getByRole('alert',{name:'Alignment error'})).toHaveCount(1);
  await page.setViewportSize({width:390,height:900});
  await inspectSafeDetails(page,alert,{code:'REGEX_SYNTAX',operation:'align',stage:'rule-validation'});
  expect(logs.join('\n')).not.toContain('PRIVATE_');
  // Exercise the real run-analysis catch after the native engine-error path.
  // The next render succeeds; its actual admission owner throws once.
  await alert.locator('details summary').press('Enter');
  await page.evaluate(() => {
    const sanitize = window.DOMPurify.sanitize;
    window.DOMPurify.sanitize = (...args) => {
      if (String(args[0]).includes('<svg')) {
        window.DOMPurify.sanitize = sanitize;
        throw Object.assign(new Error('PRIVATE_ADMISSION_SENTINEL'), {
          code: 'COMPARISON_IDENTITY', stage: 'result-admission',
          context: { reason: 'SOURCE_VIEW_CONFLICT' },
          secondary: [{ code: 'CLEANUP_FAILED', stage: 'cleanup', message: 'PRIVATE_CLEANUP_SENTINEL' }]
        });
      }
      return sanitize(...args);
    };
  });
  await exposeReviewControl(page,dialog.getByRole('button',{name:'Apply',exact:true}));
  await applyAndSettle(page,dialog);
  const caught = await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentError);
  expect(caught).toMatchObject({ code: 'COMPARISON_IDENTITY', operation: 'align', stage: 'result-admission',
    context: { reason: 'SOURCE_VIEW_CONFLICT' }, secondary: [{ code: 'CLEANUP_FAILED', stage: 'cleanup' }] });
  const afterCatch = await artifactSnapshot(page);
  expect({ ...afterCatch, history: [...afterCatch.history.slice(0, 2), before.history[2]] }).toEqual(before);
  expect(await page.evaluate(() => JSON.stringify(window.__GBDRAW_APP__.similarityAlignmentDraft))).toBe(choices);
  await expect(alert).toBeFocused();
  await inspectSafeDetails(page, alert, caught);
  expect(logs.join('\n')).not.toContain('PRIVATE_');
  await exposeReviewControl(page,dialog.getByRole('button',{name:'Apply',exact:true}));
  await applyAndSettle(page,dialog);
  await expect(dialog).toBeHidden({timeout:180000});
  expect((await artifactSnapshot(page)).history[0]).toBe(before.history[0]+1);
  expect(await page.evaluate(()=>window.__GBDRAW_APP__.errorLog)).toBeNull();
  expect(logs.join('\n')).not.toContain('PRIVATE_');
});

// OV-380: a comparison edit before Generate is a pending form edit; it keeps the
// committed alignment. OV-381: the drawer offers the active plan's exact
// reference, so Align works from the drawer after a Session load.
test('choosing Upload BLAST TSV without a file keeps a loaded alignment and the drawer offers its reference', async ({ page }) => {
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  // Save another drawer group than the plan's group.
  const seeded = JSON.parse(readFileSync(
    'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json', 'utf8'
  ));
  seeded.orthogroupState.selectedOrthogroupId = 'og_18';
  await importSession(page, Buffer.from(JSON.stringify(seeded)), 'BGC0000708-BGC0000713.gbdraw-session.json');
  const plan = () => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { reference: state.similarityAlignmentPlan.value?.reference ?? null,
      groupId: state.similarityAlignmentPlan.value?.groupId ?? null,
      notice: window.__GBDRAW_APP__.similarityAlignmentNotice };
  });
  const loaded = await plan();
  expect(loaded.groupId).toBeTruthy();
  // The Session saved another drawer group and the plan's group is not the
  // first one, so the drawer must select the plan's group itself.
  const groups = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return state.orthogroups.value.map(({ id }) => id);
  });
  const savedSelection = seeded.orthogroupState.selectedOrthogroupId;
  expect(groups[0]).not.toBe(loaded.groupId);
  expect(savedSelection).not.toBe(loaded.groupId);

  await page.evaluate(() => document.querySelectorAll('details').forEach((details) => { details.open = true; }));
  const upload = page.getByRole('radio', { name: 'Upload BLAST TSV' }).first();
  await upload.check({ force: true });
  await expect(upload).toBeChecked();
  expect(await plan()).toEqual(loaded);
  await expect(page.getByText('Alignment cleared', { exact: false })).toHaveCount(0);

  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('orthogroups'));
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return state.selectedOrthogroupId.value;
  })).toBe(loaded.groupId);
  const reference = page.locator('#similarity-alignment-reference');
  await expect(reference).toHaveValue(JSON.stringify(loaded.reference));
  await expect(page.locator('#similarity-alignment-drawer-reason'))
    .not.toContainText('Select an exact reference');
  await expect(page.locator('button[title="Align to the selected exact feature"]')).toBeEnabled();
});

test('the drawer list scrolls to the group a new alignment plan selects', async ({ page }) => {
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(
    'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'
  ), 'BGC0000708-BGC0000713.gbdraw-session.json');
  const groupId = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return state.similarityAlignmentPlan.value?.groupId || '';
  });
  expect(groupId).toBeTruthy();
  const list = page.locator('[data-similarity-group-list]');
  const planRow = list.locator(`[data-similarity-group-row="${groupId}"]`);
  const rowInView = (row) => row.evaluate((element) => {
    const box = element.getBoundingClientRect();
    const view = element.closest('[data-similarity-group-list]').getBoundingClientRect();
    return box.top >= view.top - 1 && box.bottom <= view.bottom + 1;
  });

  // The Session's plan selected its group before the drawer opened. The row
  // stays in view once the Result card has settled: the Session loaded notice
  // above the card changes its margins, which must not slide the list after it
  // scrolled (OV-1022).
  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('orthogroups'));
  await expect(planRow).toHaveAttribute('aria-current', 'true');
  await page.locator('[aria-label="Result Preview"]').evaluate(
    (card) => Promise.all(card.getAnimations().map((animation) => animation.finished)));
  await expect.poll(() => rowInView(planRow)).toBe(true);

  // With the drawer open, the user picks the first group and the list is back at
  // the top; the next plan (the same plan replaced, as Undo or a reorder does)
  // selects its group again and the list follows it.
  const firstRow = list.locator('[data-similarity-group-row]').first();
  await list.evaluate((element) => { element.scrollTop = 0; });
  await firstRow.click();
  await expect(firstRow).toHaveAttribute('aria-current', 'true');
  await expect.poll(() => rowInView(planRow)).toBe(false);
  await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    state.similarityAlignmentPlan.value = { ...state.similarityAlignmentPlan.value };
  });
  await expect(planRow).toHaveAttribute('aria-current', 'true');
  await expect.poll(() => rowInView(planRow)).toBe(true);
});

// OV-498 found the floating panel ending under the window bottom. The docked panel
// is a part of the Result instead: beside the figure when the Result is wide, below
// it when narrow, and never over it.
test('the Review alignment panel is docked beside or below the figure with Cancel and Apply in the window', async ({ page }) => {
  test.setTimeout(300000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(
    'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'
  ), 'BGC0000708-BGC0000713.gbdraw-session.json');
  await page.locator('[data-gbdraw-feature-id="f94b7f30b_record_2"]').click({ force: true });
  await page.locator('.feature-popup [data-similarity-alignment-popup-review]').click();
  const dialog = page.getByRole('dialog', { name: 'Review alignment' });
  await expect(dialog.locator('[data-alignment-record-key]')).toHaveCount(4, { timeout: 180000 });
  // The panel is in the Result's workspace, not a layer over the page.
  expect(await dialog.evaluate((element) => element.parentElement.classList.contains('preview-workspace'))).toBe(true);
  const edges = () => page.evaluate(() => {
    const box = (selector) => {
      const { left, right, top, bottom, width, height } = document.querySelector(selector).getBoundingClientRect();
      return { left, right, top, bottom, width, height };
    };
    return {
      panel: box('[data-similarity-alignment-dialog]'),
      figure: box('.preview-viewport'),
      applyBottom: box('[data-similarity-alignment-dialog] button[aria-describedby="similarity-alignment-apply-reason"]').bottom,
      headerBottom: box('header.app-header').bottom,
      footerTop: box('[data-app-footer]').top
    };
  });
  const rows = dialog.locator('.custom-scrollbar').first();
  // Wide: a right column 24rem wide over the workspace height; the figure keeps the rest.
  for (const [width, height] of [[1600, 1000], [1600, 640], [1280, 800]]) {
    await page.setViewportSize({ width, height });
    await expect.poll(async () => {
      const { panel, figure, applyBottom, headerBottom, footerTop } = await edges();
      return panel.top >= headerBottom && applyBottom < panel.bottom && panel.bottom <= footerTop
        && Math.round(panel.width) === 384 && Math.abs(panel.left - figure.right) <= 1
        && Math.abs(panel.top - figure.top) <= 1 && figure.width >= 300;
    }, { message: `window ${width}x${height}: ${JSON.stringify(await edges())}` }).toBe(true);
    // The rows scroll inside the panel when they do not fit.
    if (height === 640) expect(await rows.evaluate((element) => element.scrollHeight > element.clientHeight)).toBe(true);
  }
  // Narrow: the panel is below the figure, which keeps a usable height.
  await page.setViewportSize({ width: 390, height: 844 });
  await showPhonePane(page, 'preview');
  await expect.poll(async () => {
    const { panel, figure } = await edges();
    return panel.top >= figure.bottom - 1 && figure.height >= 200
      && Math.abs(panel.left - figure.left) <= 1 && Math.abs(panel.width - figure.width) <= 2
      && panel.right <= 390;
  }, { message: `window 390x844: ${JSON.stringify(await edges())}` }).toBe(true);
  for (const control of [dialog.getByRole('button', { name: 'Apply', exact: true }),
    dialog.getByRole('button', { name: 'Cancel', exact: true })]) {
    await exposeReviewControl(page, control);
  }
  expect(await rows.evaluate((element) => element.scrollWidth <= element.clientWidth)).toBe(true);
});

const galleryAlignmentSession = 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json';
const openPopupReview = async (page, findFeature, argument) => {
  const featureId = await page.evaluate(findFeature, argument);
  expect(featureId).toBeTruthy();
  await page.locator(`[data-gbdraw-feature-id="${featureId}"]`).first()
    .dispatchEvent('click', { clientX: -100, clientY: -100 });
  const popup = page.locator('.feature-popup[role="dialog"]');
  const popupReview = popup.getByRole('button', { name: 'Review alignment options…' });
  await popupReview.click();
  const dialog = page.getByRole('dialog', { name: 'Review alignment' });
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
  return { dialog, popup, popupReview };
};
const featureByProteinId = (proteinId) => window.__GBDRAW_APP__.extractedFeatures
  .find((item) => item.protein_id === proteinId || item.sourceProteinId === proteinId)?.svg_id;
const chooseUnresolved = async (page, dialog) => {
  for (const key of await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDraft.rows
    .filter((row) => !row.choice).map((row) => row.recordKey))) {
    await dialog.locator(`[data-alignment-record-key="${key}"]`)
      .getByRole('radio', { name: /Select .*bp, strand/ }).first().check();
  }
};
// Apply, once more when the final validation asks for a second look at the preview.
const applyReview = async (page, dialog) => {
  const apply = dialog.getByRole('button', { name: 'Apply', exact: true });
  for (let attempt = 0; attempt < 2; attempt += 1) {
    await apply.click();
    await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
    if (!await dialog.isVisible()) break;
    await expect(dialog.locator('[data-similarity-alignment-error]')).toContainText('Apply again');
  }
  await expect(dialog).toBeHidden({ timeout: 180000 });
};

test('the docked review takes the popup, the Editor and the view, and gives them back after Cancel and Apply', async ({ page }) => {
  test.setTimeout(600000);
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(String(error?.message || error)));
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(galleryAlignmentSession), 'BGC0000708-BGC0000713.gbdraw-session.json');
  // The user's own view and Editor tab.
  await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.canvasPan.x += 37;
    app.canvasPan.y += 11;
    app.openRightDrawerTab('orthogroups');
  });
  const before = await previewView(page);
  const { dialog, popup, popupReview } = await openPopupReview(page, featureByProteinId, 'CAG38712.1');
  const clicked = () => page.evaluate(() => window.__GBDRAW_APP__.clickedFeature?.svg_id
    || window.__GBDRAW_APP__.clickedFeature?.feat?.svg_id || null);
  const clickedBefore = await clicked();
  expect(clickedBefore).toBeTruthy();

  // The popup is hidden, not closed; the Editor is closed and its toggle says why.
  await expect(popup).toBeHidden();
  await expect(popup).toHaveCount(1);
  expect(await editorDrawer(page)).toEqual({ open: false, tab: 'orthogroups' });
  await expect(page.locator('.drawer-toggle')).toBeDisabled();
  await expect(page.locator('.drawer-toggle')).toHaveAttribute('title', /Finish or cancel alignment review/);
  await expect(dialog.getByRole('heading', { name: 'Review alignment' })).toBeFocused();
  // The figure is fitted to the room the panel leaves.
  await expect.poll(async () => (await previewView(page)).zoom).toBeLessThan(before.zoom);
  await expect.poll(() => page.evaluate(() => {
    const figure = document.querySelector('.gbdraw-preview-surface svg').getBoundingClientRect();
    const frame = document.querySelector('.preview-viewport').getBoundingClientRect();
    const panel = document.querySelector('[data-similarity-alignment-dialog]').getBoundingClientRect();
    return figure.left >= frame.left && figure.right <= frame.right
      && figure.top >= frame.top && figure.bottom <= frame.bottom && figure.right <= panel.left;
  })).toBe(true);
  const fitted = await previewView(page);

  // A Session saved during the review holds the view from before it.
  const savedView = await evaluateWithRetainedPromise(page, async () => {
    const { exportSession } = (await import('/tests/web/helpers/session-caller.mjs')).callerSession(await import('./js/services/config.js'), (await import('./js/state.js')).state, window.__GBDRAW_APP__.sessionImportRollbackState);
    const result = await exportSession('saved during review');
    const text = await new Response(result.blob.stream().pipeThrough(new DecompressionStream('gzip'))).text();
    const { zoom, canvasPan } = JSON.parse(text).ui;
    return { status: result.status, zoom, pan: canvasPan };
  });
  expect(savedView).toEqual({ status: 'saved', ...before });

  // No gesture on the figure opens, retargets or closes the hidden popup: a click on another
  // feature, a Ctrl-click and a Shift-click (which select), a click on the paper.
  const otherFeatureId = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const draft = app.similarityAlignmentDraft;
    const used = new Set(draft.rows.flatMap((row) => row.candidates.map(({ anchor }) => anchor.biologicalFeatureId))
      .concat([draft.response.reference.biologicalFeatureId]));
    return app.extractedFeatures.find((item) => item.type === 'CDS' && !used.has(item.biologicalFeatureId)).svg_id;
  });
  const otherFeature = page.locator(`[data-gbdraw-feature-id="${otherFeatureId}"]`).first();
  for (const gesture of [{}, { ctrlKey: true }, { shiftKey: true }]) {
    await otherFeature.dispatchEvent('click', { clientX: -100, clientY: -100, ...gesture });
    expect(await clicked(), JSON.stringify(gesture)).toBe(clickedBefore);
  }
  expect(await page.evaluate(() => window.__GBDRAW_APP__.selectedFeatures.length)).toBeGreaterThan(0);
  await page.locator('.gbdraw-preview-surface svg').dispatchEvent('click', { clientX: -100, clientY: -100 });
  expect(await clicked()).toBe(clickedBefore);
  await expect(popup).toBeHidden();
  expect(await previewView(page)).toEqual(fitted);

  // Escape ends the review: the popup, the Editor tab, the view and the focus are back.
  await dialog.getByRole('heading', { name: 'Review alignment' }).focus();
  await page.keyboard.press('Escape');
  await expect(dialog).toBeHidden();
  await expect(popup).toBeVisible();
  await expect(popupReview).toBeFocused();
  await expect.poll(() => previewView(page)).toEqual(before);
  expect(await editorDrawer(page)).toEqual({ open: true, tab: 'orthogroups' });
  await expect(page.locator('.drawer-toggle')).toBeEnabled();
  expect(await clicked()).toBe(clickedBefore);

  // Apply: the view from before the review is the view after it, and Undo keeps it.
  await popupReview.click();
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await expect(popup).toBeHidden();
  await expect.poll(async () => (await previewView(page)).zoom).toBeLessThan(before.zoom);
  await chooseUnresolved(page, dialog);
  await dialog.getByRole('radio', { name: 'All left', exact: true }).check();
  const undoBefore = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  await applyReview(page, dialog);
  await expect.poll(() => previewView(page)).toEqual(before);
  await expect(popup).toBeVisible();
  await expect(popupReview).toBeFocused();
  expect(await editorDrawer(page)).toEqual({ open: true, tab: 'orthogroups' });
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(undoBefore + 1);
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(undoBefore);
  await expect.poll(() => previewView(page)).toEqual(before);
  expect(pageErrors).toEqual([]);
});

// A review opened with the keyboard from a feature popup, whose Apply fails: the review
// stays with the error and its draft while the popup stays hidden, and Cancel gives the
// popup back with focus on Review alignment options…, beside the figure and below it.
for (const [width, height] of [[1280, 800], [390, 844]]) {
  test(`a failed Apply in a review opened by keyboard from a feature popup, then Cancel (${width}x${height})`, async ({ page }) => {
    test.setTimeout(600000);
    await page.setViewportSize({ width, height });
    await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
    await page.waitForFunction(() => window.__GBDRAW_APP__);
    await importSession(page, readFileSync(galleryAlignmentSession), 'BGC0000708-BGC0000713.gbdraw-session.json');
    // A phone user taps Preview to reach the figure (a Load keeps the shown pane).
    await showPhonePane(page, 'preview');
    const featureId = await page.evaluate(featureByProteinId, 'CAG38712.1');
    await page.locator(`[data-gbdraw-feature-id="${featureId}"]`).first()
      .dispatchEvent('click', { clientX: -100, clientY: -100 });
    const popup = page.locator('.feature-popup[role="dialog"]');
    const popupReview = popup.locator('[data-similarity-alignment-popup-review]');
    await popupReview.focus();
    await page.keyboard.press('Enter');
    const dialog = page.locator('[data-similarity-alignment-dialog]');
    await expect(dialog).toBeVisible({ timeout: 180000 });
    await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
    await expect(dialog.locator('#similarity-alignment-dialog-title')).toBeFocused();
    await expect(popup).toBeHidden();
    await chooseUnresolved(page, dialog);
    const draft = () => page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDraft.rows
      .map((row) => [row.recordKey, row.choice?.kind ?? null, row.choice?.candidateKey ?? null]));
    const chosen = await draft();
    const before = await artifactSnapshot(page);
    await page.evaluate(() => {
      const sanitize = window.DOMPurify.sanitize;
      window.DOMPurify.sanitize = (...args) => {
        if (!String(args[0]).includes('<svg')) return sanitize(...args);
        window.DOMPurify.sanitize = sanitize;
        throw Object.assign(new Error('The candidate Result was not admitted.'), {
          code: 'COMPARISON_IDENTITY', stage: 'result-admission', context: { reason: 'SOURCE_VIEW_CONFLICT' }
        });
      };
    });
    const apply = dialog.locator('footer button[aria-describedby="similarity-alignment-apply-reason"]');
    await exposeReviewControl(page, apply);
    await apply.click();
    const error = dialog.locator('[data-similarity-alignment-error]');
    await expect(error).toBeVisible({ timeout: 180000 });
    await expect(error).toBeFocused();
    await expect(apply).toBeEnabled();
    await expect(popup).toBeHidden();
    expect(await draft()).toEqual(chosen);
    expect((await artifactSnapshot(page)).results).toEqual(before.results);

    const cancel = dialog.locator('footer button.btn-secondary');
    await exposeReviewControl(page, cancel);
    await cancel.click();
    await expect(dialog).toHaveCount(0);
    await expect(popup).toBeVisible();
    await expect(popupReview).toBeFocused();
  });
}

// The search row's Open is a command, not a gesture on the figure: it moves the view and
// retargets the hidden popup, which is shown at its feature when the review ends.
test('Open in the feature search row retargets the hidden popup, shown at its feature when the review ends', async ({ page }) => {
  test.setTimeout(600000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(galleryAlignmentSession), 'BGC0000708-BGC0000713.gbdraw-session.json');
  // The whole figure is in the view the review will give back.
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  await expect.poll(async () => (await previewView(page)).zoom).toBeLessThan(1);
  const before = await previewView(page);
  const { dialog, popup } = await openPopupReview(page, featureByProteinId, 'CAG38712.1');
  const clicked = () => page.evaluate(() => window.__GBDRAW_APP__.clickedFeature?.svg_id
    || window.__GBDRAW_APP__.clickedFeature?.feat?.svg_id || null);
  const clickedBefore = await clicked();
  await expect.poll(async () => (await previewView(page)).zoom).toBeLessThan(before.zoom);
  const fitted = await previewView(page);
  const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
  await search.fill('neoA');
  await search.press('Enter');
  const openActive = page.getByRole('button', { name: 'Open active feature', exact: true });
  await expect(openActive).toBeEnabled();
  await openActive.click();
  await expect.poll(() => previewView(page)).not.toEqual(fitted);
  await expect.poll(clicked).not.toBe(clickedBefore);
  const opened = await clicked();
  expect(opened).toBeTruthy();
  await expect(popup).toBeHidden();
  await expect(dialog).toBeVisible();

  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toBeHidden();
  await expect(popup).toBeVisible();
  await expect(popup).toContainText('neoA');
  await expect.poll(() => previewView(page)).toEqual(before);
  // At the feature in the view the review gave back: beside, below or above its centre, 12 px away.
  const offFeature = () => page.evaluate((svgId) => {
    const rects = [...document.querySelectorAll(`.gbdraw-preview-surface [data-gbdraw-feature-id="${svgId}"]`)]
      .map((element) => element.getBoundingClientRect());
    const x = (Math.min(...rects.map((rect) => rect.left)) + Math.max(...rects.map((rect) => rect.right))) / 2;
    const y = (Math.min(...rects.map((rect) => rect.top)) + Math.max(...rects.map((rect) => rect.bottom))) / 2;
    const box = document.querySelector('.feature-popup[role="dialog"]').getBoundingClientRect();
    return [Math.max(box.left - x, x - box.right, 0), Math.max(box.top - y, y - box.bottom, 0)]
      .map(Math.round).sort((a, b) => a - b);
  }, opened);
  await expect.poll(offFeature).toEqual([0, 12]);
  expect(await clicked()).toBe(opened);
});

// A4: what the review took is given back only while the Result is still the review's.
test('a Session loaded under an open review keeps its own view and Editor when the review closes', async ({ page }) => {
  test.setTimeout(600000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const session = JSON.parse(readFileSync(galleryAlignmentSession, 'utf8'));
  await importSession(page, Buffer.from(JSON.stringify(session)), 'BGC0000708-BGC0000713.gbdraw-session.json');
  await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('orthogroups'));
  const before = await previewView(page);
  const { dialog } = await openPopupReview(page, featureByProteinId, 'CAG38712.1');
  await expect.poll(async () => (await previewView(page)).zoom).toBeLessThan(before.zoom);
  // Another Session arrives with its own view.
  session.ui = { ...session.ui, zoom: 0.37, canvasPan: { x: -21, y: 13 } };
  await importSession(page, Buffer.from(JSON.stringify(session)), 'BGC0000708-BGC0000713-other-view.gbdraw-session.json');
  const loaded = await previewView(page);
  expect(loaded).toEqual({ zoom: 0.37, pan: { x: -21, y: 13 } });
  // The review stays open on its draft, without rings: they belonged to the replaced Result.
  await expect(dialog).toBeVisible();
  await expect(page.locator('[data-similarity-alignment-canvas]')).toHaveCount(0);
  // A Session saved now holds the loaded view, not the view from before the review.
  const savedView = await evaluateWithRetainedPromise(page, async () => {
    const { exportSession } = (await import('/tests/web/helpers/session-caller.mjs')).callerSession(await import('./js/services/config.js'), (await import('./js/state.js')).state, window.__GBDRAW_APP__.sessionImportRollbackState);
    const result = await exportSession('saved after load');
    const text = await new Response(result.blob.stream().pipeThrough(new DecompressionStream('gzip'))).text();
    const { zoom, canvasPan } = JSON.parse(text).ui;
    return { zoom, pan: canvasPan };
  });
  expect(savedView).toEqual(loaded);
  // The view is the loaded Session's now: a change of the room does not fit it.
  await page.setViewportSize({ width: 1400, height: 1000 });
  await roomSettled(page);
  expect(await previewView(page)).toEqual(loaded);
  await page.keyboard.press('Escape');
  await expect(dialog).toBeHidden();
  expect(await previewView(page)).toEqual(loaded);
  expect((await editorDrawer(page)).open).toBe(false);
});

// OV-536: while the review shows the fitted figure, the figure follows the room of the Result,
// whatever changed that room: the window, the settings pane's display, a docked section, the
// pane's width.
test('the figure is fitted again whenever the room of the Result changes while Review alignment is open', async ({ page }) => {
  test.setTimeout(600000);
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(String(error?.message || error)));
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(galleryAlignmentSession), 'BGC0000708-BGC0000713.gbdraw-session.json');
  // The user's own view.
  await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
  await page.evaluate(() => { window.__GBDRAW_APP__.canvasPan.x += 37; });
  const before = await previewView(page);
  const { dialog } = await openPopupReview(page, featureByProteinId, 'CAG38712.1');
  const undoAtOpen = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  // The figure is fitted in the canvas as it is now, whose size differs from the one before.
  let room = null;
  const fittedAfterChange = async (state) => {
    const last = String(room);
    await expect.poll(async () => {
      const found = await figureFit(page);
      return found.fitted && found.inFrame && String(found.room) !== last;
    }, { message: `${state}: ${JSON.stringify(await figureFit(page))}` }).toBe(true);
    room = (await figureFit(page)).room;
    return { ...await previewView(page), room };
  };
  const opened = await fittedAfterChange('opened');
  // The same room gives the same fit, whatever the view it replaces.
  const expectOpenedView = (view) => {
    expect(view.room).toEqual(opened.room);
    expect(view.zoom).toBe(opened.zoom);
    expect(Math.abs(view.pan.x - opened.pan.x) + Math.abs(view.pan.y - opened.pan.y)).toBeLessThan(0.01);
  };

  // The window: wider, narrower, and back. The fit has landed once two frames are drawn, which
  // is what the cases that expect no fit wait for.
  await page.setViewportSize({ width: 1900, height: 1000 });
  await roomSettled(page);
  expect((await previewView(page)).zoom).toBeGreaterThan(opened.zoom);
  const wider = await fittedAfterChange('window 1900');
  expect(wider.room[0]).toBeGreaterThan(opened.room[0]);
  await page.setViewportSize({ width: 1340, height: 1000 });
  const narrower = await fittedAfterChange('window 1340');
  expect(narrower.room[0]).toBeLessThan(opened.room[0]);
  expect(narrower.zoom).toBeLessThan(opened.zoom);
  await page.setViewportSize({ width: 1600, height: 1000 });
  expectOpenedView(await fittedAfterChange('window 1600 again'));

  // The settings pane: the icon rail with a section, the rail alone, the stacked pane again.
  await page.getByRole('button', { name: 'Collapse settings to icons', exact: true }).click();
  const railWithSection = await fittedAfterChange('rail with a section');
  expect(railWithSection.room[0]).toBeLessThan(opened.room[0]);
  expect(railWithSection.zoom).toBeLessThan(opened.zoom);
  await page.getByRole('button', { name: 'Close section', exact: true }).click();
  const railAlone = await fittedAfterChange('rail alone');
  expect(railAlone.room[0]).toBeGreaterThan(opened.room[0]);
  expect(railAlone.zoom).toBeGreaterThan(opened.zoom);
  await page.getByRole('button', { name: 'Expand settings', exact: true }).click();
  expectOpenedView(await fittedAfterChange('stacked again'));

  // A zoom and a pan made by hand stay while the room stays: through a hover, and through the
  // canvas hidden and shown again in the same room (the Run info tab) ...
  await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
  await page.evaluate(() => { window.__GBDRAW_APP__.canvasPan.y += 23; });
  const byHand = await previewView(page);
  expect(byHand.zoom).toBeGreaterThan(opened.zoom);
  await dialog.locator('[data-alignment-record-key]').first().hover();
  await roomSettled(page);
  expect(await previewView(page)).toEqual(byHand);
  for (const tab of ['run-info', 'preview']) {
    await page.evaluate((shown) => { window.__GBDRAW_APP__.resultPanelTab = shown; }, tab);
    await roomSettled(page);
  }
  expect(await previewView(page)).toEqual(byHand);
  // ... and the fit replaces them when it changes: here the pane's width drag, which follows
  // the pointer from the first move (OV-620).
  const paneWidth = () => page.evaluate(() => window.__GBDRAW_APP__.sidebarWidth);
  const paneBefore = await paneWidth();
  const handle = await page.locator('.resize-handle').first().boundingBox();
  await page.mouse.move(handle.x + handle.width / 2, handle.y + 200);
  await page.mouse.down();
  await page.mouse.move(handle.x + handle.width / 2 + 1, handle.y + 200);
  expect(await paneWidth()).toBe(paneBefore + 1);
  await page.mouse.move(handle.x + handle.width / 2 + 100, handle.y + 200, { steps: 5 });
  await page.mouse.up();
  expect(await paneWidth()).toBe(paneBefore + 100);
  const dragged = await fittedAfterChange('pane 100 px wider');
  expect(dragged.room[0]).toBe(opened.room[0] - 100);
  expect(dragged.zoom).toBeLessThan(opened.zoom);

  // No fit is a History step, and a Session saved now holds the view from before the review.
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(undoAtOpen);
  const savedView = await evaluateWithRetainedPromise(page, async () => {
    const { exportSession } = (await import('/tests/web/helpers/session-caller.mjs')).callerSession(await import('./js/services/config.js'), (await import('./js/state.js')).state, window.__GBDRAW_APP__.sessionImportRollbackState);
    const result = await exportSession('saved after a fit to the room');
    const text = await new Response(result.blob.stream().pipeThrough(new DecompressionStream('gzip'))).text();
    const { zoom, canvasPan } = JSON.parse(text).ui;
    return { zoom, pan: canvasPan };
  });
  expect(savedView).toEqual(before);

  // Apply records its one History action and gives back the view from before the review.
  await chooseUnresolved(page, dialog);
  await dialog.getByRole('radio', { name: 'All left', exact: true }).check();
  await applyReview(page, dialog);
  await expect.poll(() => previewView(page)).toEqual(before);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(undoAtOpen + 1);
  // The review is over: a change of the room keeps that view.
  await page.setViewportSize({ width: 1900, height: 1000 });
  await roomSettled(page);
  expect(await previewView(page)).toEqual(before);
  expect(pageErrors).toEqual([]);
});

test('outside Review alignment a change of the room of the Result keeps the zoom and the position', async ({ page }) => {
  test.setTimeout(300000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(galleryAlignmentSession), 'BGC0000708-BGC0000713.gbdraw-session.json');
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
  await page.evaluate(() => { window.__GBDRAW_APP__.canvasPan.x += 37; });
  const view = await previewView(page);
  const room = async () => (await figureFit(page)).room;
  const rooms = [await room()];
  for (const change of [
    () => page.setViewportSize({ width: 1900, height: 1000 }),
    () => page.setViewportSize({ width: 1340, height: 800 }),
    () => page.getByRole('button', { name: 'Collapse settings to icons', exact: true }).click(),
    () => page.getByRole('button', { name: 'Close section', exact: true }).click()
  ]) {
    await change();
    await roomSettled(page);
    rooms.push(await room());
    expect(await previewView(page), `room ${rooms.at(-1)}`).toEqual(view);
  }
  // Each step did change the room.
  expect(new Set(rooms.map(String)).size).toBe(rooms.length);
});

// The user's own view and Editor tab, then the review opened from a feature popup.
const openReviewOverOwnView = async (page) => {
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(galleryAlignmentSession), 'BGC0000708-BGC0000713.gbdraw-session.json');
  await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.canvasPan.x += 37;
    app.canvasPan.y += 11;
    app.openRightDrawerTab('orthogroups');
  });
  const before = await previewView(page);
  const { dialog } = await openPopupReview(page, featureByProteinId, 'CAG38712.1');
  await expect.poll(async () => (await figureFit(page)).fitted).toBe(true);
  return { before, dialog };
};

// One live setting change: the CDS color in Colors, drawn on the shown Result.
const cdsColor = (page) => page.locator('input[type="color"][aria-label="CDS feature color"]').first();
const recolorCds = async (page, color) => {
  await cdsColor(page).evaluate((input, value) => {
    input.value = value;
    input.dispatchEvent(new Event('input', { bubbles: true }));
    input.dispatchEvent(new Event('change', { bubbles: true }));
  }, color);
  await expect.poll(() => page.evaluate((value) => window.__GBDRAW_APP__.currentResults()
    .some(({ content }) => String(content).toLowerCase().includes(value)), color), { timeout: 60000 }).toBe(true);
  await roomSettled(page);
};
// Whether each Linear record is drawn reversed.
const recordOrientations = (page) => page.evaluate(async () => (await import('./js/state.js')).state.drawingList().find((drawing) => drawing.mode === 'linear').linearSeqs
  .map(({ uid, region_reverse: reverse }) => [uid, Boolean(reverse)]));
// The Result and the figure draw the color the setting shows (PD-OI-102).
const expectCdsColorDrawn = async (page, color) => {
  await expect(cdsColor(page)).toHaveValue(color);
  expect(await page.evaluate((value) => ({
    results: window.__GBDRAW_APP__.currentResults().every(({ content }) => String(content).toLowerCase().includes(value)),
    figure: document.querySelector('.gbdraw-preview-surface svg').outerHTML.toLowerCase().includes(value)
  }), color)).toEqual({ results: true, figure: true });
};

// OV-1082: a live setting change redraws the Result of the same records (PD-OI-102). The review
// still shows its draft on it: the rings are drawn again, the figure follows the room, and the
// end (Cancel, Apply) gives back the view and the Editor. Apply draws the change (OV-1083).
test('a live setting change under an open review keeps its rings, its fit and the view it gives back', async ({ page }) => {
  test.setTimeout(600000);
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(String(error?.message || error)));
  const { before, dialog } = await openReviewOverOwnView(page);
  const rings = page.locator('button.gbdraw-alignment-ring');
  const ringCount = await rings.count();
  expect(ringCount).toBeGreaterThan(0);
  const fitted = await previewView(page);
  const recolor = (color) => recolorCds(page, color);
  await recolor('#e76f51');
  await expect(dialog).toBeVisible();
  await expect(rings).toHaveCount(ringCount);
  expect(await previewView(page)).toEqual(fitted);
  // The figure is fitted again when the room changes.
  const room = (await figureFit(page)).room;
  await page.setViewportSize({ width: 1500, height: 1000 });
  await expect.poll(async () => {
    const found = await figureFit(page);
    return found.fitted && found.inFrame && String(found.room) !== String(room);
  }).toBe(true);
  // Cancel gives back the view and the Editor from before the review.
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toBeHidden();
  await expect.poll(() => previewView(page)).toEqual(before);
  expect(await editorDrawer(page)).toEqual({ open: true, tab: 'orthogroups' });

  // Apply after a live change: one History action, and the view from before the review.
  const { dialog: again } = await openPopupReview(page, featureByProteinId, 'CAG38712.1');
  await recolor('#2a9d8f');
  await expect(rings).toHaveCount(ringCount);
  const undoBefore = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const preApply = await recordOrientations(page);
  await chooseUnresolved(page, again);
  await again.getByRole('radio', { name: 'All left', exact: true }).check();
  await applyReview(page, again);
  await expect.poll(() => previewView(page)).toEqual(before);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(undoBefore + 1);
  expect(await editorDrawer(page)).toEqual({ open: true, tab: 'orthogroups' });
  await expectCdsColorDrawn(page, '#2a9d8f');
  expect(await recordOrientations(page)).not.toEqual(preApply);
  // Undo puts back the orientation from before Apply, and the live color stays.
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await expect.poll(() => recordOrientations(page)).toEqual(preApply);
  await expectCdsColorDrawn(page, '#2a9d8f');
  expect(pageErrors).toEqual([]);
});

// The font sizes of the Legend the figure draws.
const legendFontSizes = (page) => page.evaluate(() => [...new Set([...document.querySelectorAll(
  '.gbdraw-preview-surface svg g[data-legend-key] text'
)].map((text) => String(parseFloat(getComputedStyle(text).fontSize))))]);
// A live setting row the automatic rerender draws: type the Legend font size and leave the field.
const commitLegendFont = async (page, size) => {
  const input = page.locator('input[aria-label="Legend font size"]');
  for (const details of await input.locator('xpath=ancestor::details').all()) {
    if (await details.getAttribute('open') === null) await details.locator(':scope > summary').press('Enter');
  }
  await input.fill(String(size));
  await input.press('Tab');
  await expect.poll(() => legendFontSizes(page), { timeout: 120000 }).toEqual([String(size)]);
  await roomSettled(page);
};

// OV-1120: the automatic rerender draws a live setting (here the Legend font size) as a new
// Result with a new feature catalog of the same features (PD-OI-102). The review keeps its
// draft, as after a live color, and Apply commits onto the Session that names the new value.
test('a live setting drawn by the automatic rerender under an open review keeps its rings and the view it gives back', async ({ page }) => {
  test.setTimeout(600000);
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(String(error?.message || error)));
  const { before, dialog } = await openReviewOverOwnView(page);
  const rings = page.locator('button.gbdraw-alignment-ring');
  const ringCount = await rings.count();
  expect(ringCount).toBeGreaterThan(0);
  const fontBefore = await legendFontSizes(page);
  await page.evaluate(async () => {
    window.__OV1120Catalog = (await import('./js/state.js')).state.featureCatalog.value;
  });
  await commitLegendFont(page, 21);
  expect(await page.evaluate(async () => (
    (await import('./js/state.js')).state.featureCatalog.value !== window.__OV1120Catalog
  )), 'the rerender installed a new feature catalog').toBe(true);
  await expect(dialog).toBeVisible();
  await expect(rings).toHaveCount(ringCount);
  // Its Undo is a live change too (D-44: B), and so is the Redo.
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await expect.poll(() => legendFontSizes(page), { timeout: 120000 }).toEqual(fontBefore);
  await expect(rings).toHaveCount(ringCount);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await expect.poll(() => legendFontSizes(page), { timeout: 120000 }).toEqual(['21']);
  await roomSettled(page);
  await expect(rings).toHaveCount(ringCount);
  // The figure is fitted again when the room changes.
  const room = (await figureFit(page)).room;
  await page.setViewportSize({ width: 1500, height: 1000 });
  await expect.poll(async () => {
    const found = await figureFit(page);
    return found.fitted && found.inFrame && String(found.room) !== String(room);
  }).toBe(true);
  // Cancel gives back the view and the Editor from before the review.
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toBeHidden();
  await expect.poll(() => previewView(page)).toEqual(before);
  expect(await editorDrawer(page)).toEqual({ open: true, tab: 'orthogroups' });

  // Apply after the change keeps the value the rerender committed.
  const { dialog: again } = await openPopupReview(page, featureByProteinId, 'CAG38712.1');
  await commitLegendFont(page, 23);
  await expect(rings).toHaveCount(ringCount);
  const preApply = await recordOrientations(page);
  await chooseUnresolved(page, again);
  await again.getByRole('radio', { name: 'All left', exact: true }).check();
  await applyReview(page, again);
  await expect.poll(() => previewView(page)).toEqual(before);
  expect(await recordOrientations(page)).not.toEqual(preApply);
  expect(await legendFontSizes(page)).toEqual(['23']);
  await expect(page.locator('input[aria-label="Legend font size"]')).toHaveValue('23');
  expect(pageErrors).toEqual([]);
});

// OV-1083: a live color changed before the review is drawn by Apply too.
test('Apply of a review draws a live color changed before the review', async ({ page }) => {
  test.setTimeout(600000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(galleryAlignmentSession), 'BGC0000708-BGC0000713.gbdraw-session.json');
  await recolorCds(page, '#264653');
  const { dialog } = await openPopupReview(page, featureByProteinId, 'CAG38712.1');
  await chooseUnresolved(page, dialog);
  await dialog.getByRole('radio', { name: 'All left', exact: true }).check();
  await applyReview(page, dialog);
  await expectCdsColorDrawn(page, '#264653');
});

// A Session loaded while the other mode is shown replaces the review's figure: the review gives
// back neither the view nor the Editor (OV-622). This pins the outcome; the load shows the
// review's mode again, so it does not tell which mode the replacement test reads.
test('a Session loaded while the other mode is shown under a review keeps its own view', async ({ page }) => {
  test.setTimeout(600000);
  const { dialog } = await openReviewOverOwnView(page);
  await page.evaluate(() => window.__GBDRAW_APP__.selectDrawing('circular'));
  await expect.poll(() => page.evaluate(async () => (await import('./js/state.js')).state.mode.value)).toBe('circular');
  const session = JSON.parse(readFileSync(galleryAlignmentSession, 'utf8'));
  session.ui = { ...session.ui, zoom: 0.37, canvasPan: { x: -21, y: 13 } };
  await importSession(page, Buffer.from(JSON.stringify(session)), 'BGC0000708-BGC0000713-other-view.gbdraw-session.json');
  const loaded = await previewView(page);
  expect(loaded).toEqual({ zoom: 0.37, pan: { x: -21, y: 13 } });
  await page.keyboard.press('Escape');
  await expect(dialog).toBeHidden();
  expect(await previewView(page)).toEqual(loaded);
  expect((await editorDrawer(page)).open).toBe(false);
});

// OV-622: a mode switch brings no view of its own. After Circular and back, the review still
// holds the view from before it and gives it back; the figure is not fitted again (the draft no
// longer shows on it, as before).
test('a mode switch and back under an open review leaves the view the review gives back', async ({ page }) => {
  test.setTimeout(600000);
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(String(error?.message || error)));
  const { before, dialog } = await openReviewOverOwnView(page);
  const shownMode = () => page.evaluate(async () => (await import('./js/state.js')).state.mode.value);
  for (const target of ['circular', 'linear']) {
    await page.evaluate((next) => window.__GBDRAW_APP__.selectDrawing(next), target);
    await expect.poll(shownMode).toBe(target);
  }
  await expect(dialog).toBeVisible();
  await roomSettled(page);
  const shown = await previewView(page);
  await page.setViewportSize({ width: 1400, height: 1000 });
  await roomSettled(page);
  expect(await previewView(page)).toEqual(shown);
  // Circular has no Similarity tab, so the switch left the Editor on its Features tab; the
  // review gives back the tab it closed (OV-1255).
  expect((await editorDrawer(page)).tab).toBe('features');
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toBeHidden();
  await expect.poll(() => previewView(page)).toEqual(before);
  expect(await editorDrawer(page)).toEqual({ open: true, tab: 'orthogroups' });
  expect(pageErrors).toEqual([]);
});

// OV-1089: the held view is the Linear review's. Escape while Circular is shown ends the
// review there and gives back the Editor (PD-OI-039 rev 4), but not Linear's view from before
// the review: the Circular figure keeps its own.
test('a Linear review ended in Circular gives the Linear view back to no other figure (OV-1089)', async ({ page }) => {
  test.setTimeout(600000);
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(String(error?.message || error)));
  const { before, dialog } = await openReviewOverOwnView(page);
  const shownMode = () => page.evaluate(async () => (await import('./js/state.js')).state.mode.value);
  await page.evaluate(() => window.__GBDRAW_APP__.selectDrawing('circular'));
  await expect.poll(shownMode).toBe('circular');
  // PD-OI-039 rev 4: the review stays open across the switch (Circular has no Result here, so
  // its panel has no Result area to dock in), and Escape ends it there.
  expect(await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDialogOpen)).toBe(true);
  const circular = await previewView(page);
  expect(circular).not.toEqual(before);
  await page.locator('body').focus();
  await page.keyboard.press('Escape');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDialogOpen)).toBe(false);
  await expect(dialog).toHaveCount(0);
  await roomSettled(page);
  expect(await previewView(page)).toEqual(circular);
  expect((await editorDrawer(page)).open).toBe(true);
  expect(pageErrors).toEqual([]);
});

// OV-1255: Apply's new Result leaves the Similarity groups tab without groups for a while; the
// review still gives back the tab it closed, as Cancel does.
test('a successful Apply of a review opened from the Editor reopens the Similarity groups tab', async ({ page }, testInfo) => {
  test.setTimeout(600000);
  await loadAmbiguousSession(page, testInfo);
  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('orthogroups'));
  const drawer = page.locator('.right-drawer');
  await drawer.locator('[data-similarity-group-row="og_1"]').click();
  const reference = await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDrawerReferenceOptions('og_1')
    .find(({ anchor }) => anchor.recordKey === 'record_b').key);
  await drawer.getByLabel('Exact reference record and feature').selectOption(reference);
  await drawer.getByRole('button', { name: 'Review alignment options…' }).click();
  const dialog = page.getByRole('dialog', { name: 'Review alignment' });
  await expect(dialog).toBeVisible({ timeout: 180000 });
  await page.waitForFunction(() => !window.__GBDRAW_APP__.similarityAlignmentBusy, null, { timeout: 180000 });
  expect(await editorDrawer(page)).toEqual({ open: false, tab: 'orthogroups' });
  await chooseUnresolved(page, dialog);
  await applyReview(page, dialog);
  await expect.poll(() => editorDrawer(page)).toEqual({ open: true, tab: 'orthogroups' });
});

// The majanivirus og_13 review opens with records that need a choice.
test('records that need a choice come first; a ring or a row chooses, and the chosen box is on the line', async ({ page }) => {
  test.setTimeout(600000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync('gbdraw/web/gallery/sessions/majanivirus_orthogroup.gbdraw-session.json.gz'),
    'majanivirus_orthogroup.gbdraw-session.json.gz');
  const { dialog } = await openPopupReview(page, (groupId) => window.__GBDRAW_APP__.extractedFeatures
    .find((item) => (item.orthogroup_id || item.orthogroupId || item.similarityGroupId) === groupId)?.svg_id, 'og_13');
  const needKeys = await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDraft.rows
    .filter((row) => row.choice === null).map((row) => row.recordKey));
  expect(needKeys.length).toBeGreaterThan(0);
  // The reference card, then the records that need a choice under their own heading.
  const order = () => dialog.locator('.alignment-review-section, .alignment-review-card').evaluateAll((elements) => (
    elements.map((element) => element.matches('h3') ? '#' + element.dataset.alignmentSection
      : element.hasAttribute('data-similarity-alignment-reference') ? 'reference' : element.dataset.alignmentRecordKey)));
  const opened = await order();
  expect(opened.slice(0, 2)).toEqual(['reference', '#choose']);
  expect(opened.slice(2, 2 + needKeys.length).sort()).toEqual([...needKeys].sort());
  expect(opened[2 + needKeys.length] ?? '#').toMatch(/^#/);
  await expect(dialog.locator('[data-alignment-section="choose"]'))
    .toHaveText(`Choose a gene${needKeys.length} (${needKeys.length} record${needKeys.length === 1 ? '' : 's'})`);
  await expect(dialog.locator('#similarity-alignment-dialog-description'))
    .toContainText(`${needKeys.length} ${needKeys.length === 1 ? 'record' : 'records'} awaiting a choice`);
  await expect(page.locator('#similarity-alignment-apply-reason')).toContainText(/Choose a candidate or Skip for \d+ record/);
  await expect(dialog.getByRole('button', { name: 'Apply', exact: true })).toBeDisabled();
  // No number joins a row to the figure.
  await expect(page.locator('.gbdraw-alignment-badge')).toHaveCount(0);
  await expect(dialog.locator('[title="Number on the diagram"]')).toHaveCount(0);

  const needKey = needKeys[0];
  const card = dialog.locator(`[data-alignment-record-key="${needKey}"]`);
  const strip = card.locator('[data-alignment-strip]');
  await expect(card).toHaveClass(/is-need/);
  const chooseHeading = dialog.locator('[data-alignment-section="choose"]');
  await expect(chooseHeading).toHaveClass(/is-need/);
  await expect(card).toContainText('Choose the gene to line up, or skip this record.');
  await expect(strip).toHaveCount(1);
  await expect(strip.locator('.candidate.is-chosen')).toHaveCount(0);
  const chosenBox = () => strip.evaluate((svg) => {
    const chosen = svg.querySelector('.candidate.is-chosen');
    if (!chosen) return null;
    const centre = Number(chosen.getAttribute('x')) + Number(chosen.getAttribute('width')) / 2;
    return { key: chosen.dataset.alignmentCandidateKey,
      offLine: Math.abs(centre - Number(svg.querySelector('.line').getAttribute('x1'))) };
  });
  const boxKeys = await strip.locator('.candidate').evaluateAll((boxes) => boxes
    .map((box) => ({ key: box.dataset.alignmentCandidateKey, x: Number(box.getAttribute('x')) }))
    .sort((left, right) => left.x - right.x).map(({ key }) => key));
  expect(boxKeys.length).toBeGreaterThan(1);
  const { candidates } = await reviewCard(page, needKey);
  // The rows list the candidates as the strip and the figure show them, left to right.
  expect(candidates.slice(0, boxKeys.length).map(({ key }) => key)).toEqual(boxKeys);
  const ringCentres = [];
  for (const key of boxKeys) {
    const ring = await candidateRing(page, key).boundingBox();
    ringCentres.push(ring.x + ring.width / 2);
  }
  expect([...ringCentres].sort((a, b) => a - b)).toEqual(ringCentres);

  // A ring lights its row and its strip box, names the gene on the figure, and chooses.
  const ring = candidateRing(page, boxKeys[1]);
  await expect(ring).toHaveClass(/is-need/);
  await ring.hover();
  await expect(card.locator('.alignment-review-candidate.is-hover')).toHaveCount(1);
  await expect(strip.locator('.candidate.is-hover')).toHaveAttribute('data-alignment-candidate-key', boxKeys[1]);
  await expect(page.locator('.gbdraw-alignment-tag')).toHaveText(candidates[1].name);
  await ring.click();
  await expect(card.getByRole('radio', { name: /Select .*bp, strand/ }).nth(1)).toBeChecked();
  await expect(candidateRing(page, boxKeys[1])).toHaveAttribute('aria-pressed', 'true');
  await expect.poll(chosenBox).toEqual({ key: boxKeys[1], offLine: expect.any(Number) });
  expect((await chosenBox()).offLine).toBeLessThan(0.01);
  // A row chooses another candidate: its box is the one on the line now.
  await card.getByRole('radio', { name: /Select .*bp, strand/ }).first().check();
  await expect.poll(async () => (await chosenBox())?.key).toBe(boxKeys[0]);
  expect((await chosenBox()).offLine).toBeLessThan(0.01);
  await expect(candidateRing(page, boxKeys[0])).toHaveAttribute('aria-pressed', 'true');
  // The record keeps its place under the heading it opened with.
  expect(await order()).toEqual(opened);
  await expect(card).not.toHaveClass(/is-need/);
  // The heading is amber only while a record under it still needs a choice.
  await expect.poll(() => chooseHeading.evaluate((heading) => heading.classList.contains('is-need')))
    .toBe(needKeys.length > 1);
  // On the Run info tab the figure is not shown and nothing is measured: the strips keep the
  // figure they had, and the rings return with the preview.
  const stripCount = await dialog.locator('[data-alignment-strip]').count();
  await page.evaluate(() => { window.__GBDRAW_APP__.resultPanelTab = 'run-info'; });
  await card.getByRole('radio', { name: /Select .*bp, strand/ }).nth(1).check();
  await expect.poll(async () => (await chosenBox())?.key).toBe(boxKeys[1]);
  expect((await chosenBox()).offLine).toBeLessThan(0.01);
  await expect(dialog.locator('[data-alignment-strip]')).toHaveCount(stripCount);
  await page.evaluate(() => { window.__GBDRAW_APP__.resultPanelTab = 'preview'; });
  await expect(candidateRing(page, boxKeys[1])).toHaveAttribute('aria-pressed', 'true');
  await expect(candidateRing(page, boxKeys[1])).toBeVisible();
  await page.keyboard.press('Escape');
  await expect(dialog).toBeHidden();
});

// C2: every record is drawn in its strip as the figure draws it after Apply. A feature drawn in
// several parts (exons) counts as one span from its first block to its last.
// Each strip box: where its centre lies along its record, as a fraction of the record's extent.
const stripBoxes = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const draft = app.similarityAlignmentDraft;
  const number = (element, name) => Number(element.getAttribute(name));
  return app.similarityAlignmentReview.groups.flatMap((group) => group.cards).flatMap((card) => {
    const selector = card.kind === 'reference' ? '[data-similarity-alignment-reference]'
      : `[data-alignment-record-key="${card.recordKey}"]`;
    const svg = document.querySelector(`[data-similarity-alignment-dialog] ${selector} [data-alignment-strip]`);
    if (!svg) return [];
    const axis = svg.querySelector('.axis');
    const fraction = (x) => (x - number(axis, 'x1')) / (number(axis, 'x2') - number(axis, 'x1'));
    const line = number(svg.querySelector('.line'), 'x1');
    return [...svg.querySelectorAll('.candidate')].map((box) => {
      const key = box.dataset.alignmentCandidateKey;
      const anchor = card.kind === 'reference' ? draft.response.reference
        : card.candidates.find((candidate) => candidate.key === key).anchor;
      const centre = number(box, 'x') + number(box, 'width') / 2;
      return { recordKey: anchor.recordKey, biologicalFeatureId: anchor.biologicalFeatureId,
        chosen: box.classList.contains('is-chosen'), fraction: fraction(centre), side: Math.sign(centre - line),
        onLine: Math.abs(centre - line) };
    });
  });
});
// The drawn span of each feature in the Result, in client pixels, with its block count.
const drawnSpans = (page, anchors) => page.evaluate((anchors) => {
  const app = window.__GBDRAW_APP__;
  const svg = document.querySelector('.gbdraw-preview-surface svg');
  return anchors.map(({ recordKey, biologicalFeatureId }) => {
    const feature = app.extractedFeatures.find((item) => item.recordKey === recordKey
      && item.biologicalFeatureId === biologicalFeatureId);
    const rects = feature ? [...svg.querySelectorAll(`[data-gbdraw-feature-id="${feature.svg_id}"]`)]
      .filter((element) => element.getAttribute('fill') !== 'none').map((element) => element.getBoundingClientRect()) : [];
    if (!rects.length) return { blocks: 0 };
    const left = Math.min(...rects.map((rect) => rect.left));
    const right = Math.max(...rects.map((rect) => rect.right));
    const group = svg.querySelector(`g[data-record-key="${recordKey}"]`);
    const extent = [...group.querySelectorAll(':scope > line:not([data-label-feature-id]), [data-gbdraw-feature-id]')]
      .map((element) => element.getBoundingClientRect()).filter((rect) => rect.width > 0 || rect.height > 0);
    const recordLeft = Math.min(...extent.map((rect) => rect.left));
    const recordRight = Math.max(...extent.map((rect) => rect.right));
    return { blocks: rects.length, left, right, x: (left + right) / 2,
      fraction: ((left + right) / 2 - recordLeft) / (recordRight - recordLeft) };
  });
}, anchors);
const expectFigureAgreesWithStrips = async (page, strips) => {
  const figure = await drawnSpans(page, strips);
  const lineX = figure[strips.findIndex(({ chosen }) => chosen)].x;
  strips.forEach((box, index) => {
    const drawn = { fraction: figure[index].fraction, offLine: figure[index].x - lineX };
    const label = `${box.recordKey} ${box.biologicalFeatureId}: strip ${JSON.stringify(box)}, figure ${JSON.stringify(drawn)}`;
    expect(Math.abs(box.fraction - drawn.fraction), label).toBeLessThan(0.005);
    if (box.chosen) expect(Math.abs(drawn.offLine), label).toBeLessThanOrEqual(1);
    else if (Math.abs(drawn.offLine) > 1) expect(Math.sign(drawn.offLine), label).toBe(box.side);
  });
};

test('the strips of reversed records agree with the figure after Apply in order and side', async ({ page }) => {
  test.setTimeout(600000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync(galleryAlignmentSession), 'BGC0000708-BGC0000713.gbdraw-session.json');
  const { dialog } = await openPopupReview(page, featureByProteinId, 'CAG38712.1');
  await chooseUnresolved(page, dialog);
  await dialog.getByRole('radio', { name: 'All left', exact: true }).check();
  const reversed = await page.evaluate(() => window.__GBDRAW_APP__.similarityAlignmentDirectionPreview.records
    .filter((record) => record.reversed).map((record) => record.recordKey));
  expect(reversed.length).toBeGreaterThan(1);
  const strips = await stripBoxes(page);
  expect(strips.filter(({ recordKey }) => reversed.includes(recordKey)).length).toBeGreaterThan(0);
  for (const box of strips.filter(({ chosen }) => chosen)) expect(box.onLine).toBeLessThan(0.01);
  await applyReview(page, dialog);
  await expectFigureAgreesWithStrips(page, strips);
});

// OV-541: a gene with introns is drawn as several blocks, and Apply centres its whole span.
test('a reference and candidates drawn in several parts get the line, rings and strip boxes around their whole span', async ({ page }) => {
  test.setTimeout(900000);
  await page.setViewportSize({ width: 1600, height: 1000 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, readFileSync('gbdraw/web/gallery/sessions/majanivirus_orthogroup.gbdraw-session.json.gz'),
    'majanivirus_orthogroup.gbdraw-session.json.gz');
  const { dialog } = await openPopupReview(page, (groupId) => window.__GBDRAW_APP__.extractedFeatures
    .find((item) => (item.orthogroup_id || item.orthogroupId || item.similarityGroupId) === groupId)?.svg_id, 'og_3');
  const review = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const cards = app.similarityAlignmentReview.groups.flatMap((group) => group.cards);
    return { reference: app.similarityAlignmentDraft.response.reference, cards: cards.length,
      candidates: cards.flatMap((card) => card.candidates.map(({ key, anchor }) => ({ key, ...anchor }))) };
  });
  expect((await drawnSpans(page, [review.reference]))[0].blocks).toBeGreaterThan(1);
  const centre = (box) => box.x + box.width / 2;
  const line = page.locator('.gbdraw-alignment-line:not(.casing)');
  await expect(line).toHaveCount(1);
  // The view is fitted when the panel opens: the line follows the feature until both rest.
  await page.waitForFunction(() => !document.querySelector('.gbdraw-preview-surface').getAnimations().length);
  const lineOffReference = async () => {
    const lineX = await line.evaluate((element) => {
      const rect = element.getBoundingClientRect();
      return (rect.left + rect.right) / 2;
    });
    return Math.abs(lineX - (await drawnSpans(page, [review.reference]))[0].x);
  };
  await expect.poll(lineOffReference).toBeLessThanOrEqual(1);
  const [reference] = await drawnSpans(page, [review.reference]);
  const referenceRing = await page.locator('.gbdraw-alignment-ring.is-reference').boundingBox();
  expect(Math.abs(centre(referenceRing) - reference.x)).toBeLessThanOrEqual(1);
  expect(referenceRing.width).toBeGreaterThanOrEqual(reference.right - reference.left);
  // The long gene name is cut in the flag; "· reference" stays whole.
  expect(await page.locator('.gbdraw-alignment-flag').evaluate((flag) => {
    const range = document.createRange();
    range.selectNodeContents(flag.lastChild);
    const words = range.getBoundingClientRect();
    const box = flag.getBoundingClientRect();
    return { words: flag.lastChild.textContent, inside: words.left >= box.left && words.right <= box.right,
      nameCut: flag.firstElementChild.scrollWidth > flag.firstElementChild.clientWidth };
  })).toEqual({ words: ' · reference', inside: true, nameCut: true });
  await expect(dialog.locator('[data-alignment-strip]')).toHaveCount(review.cards);
  // A candidate in several parts has one ring around all of them.
  const drawn = await drawnSpans(page, review.candidates);
  const inParts = review.candidates.map((candidate, index) => ({ ...candidate, ...drawn[index] }))
    .filter(({ blocks }) => blocks > 1);
  expect(inParts.length).toBeGreaterThan(0);
  for (const candidate of inParts) {
    const ring = await candidateRing(page, candidate.key).boundingBox();
    expect(Math.abs(centre(ring) - candidate.x), candidate.key).toBeLessThanOrEqual(1);
    expect(ring.width, candidate.key).toBeGreaterThanOrEqual(candidate.right - candidate.left);
  }
  await chooseUnresolved(page, dialog);
  const strips = await stripBoxes(page);
  expect(strips.filter(({ chosen }) => chosen).length).toBeGreaterThan(2);
  for (const box of strips.filter(({ chosen }) => chosen)) expect(box.onLine).toBeLessThan(0.01);
  await applyReview(page, dialog);
  await expectFigureAgreesWithStrips(page, strips);
});
