const fs = require('node:fs');
const { execFile } = require('node:child_process');
const { join } = require('node:path');
const { pathToFileURL } = require('node:url');
const { promisify } = require('node:util');
const { gunzipSync } = require('node:zlib');
const { test, expect } = require('@playwright/test');
const {
  CURRENT_FEATURE_CATALOG_SCHEMA,
  CURRENT_REQUEST_SCHEMA,
  CURRENT_SESSION_VERSION,
  evaluateWithRetainedPromise,
  generateAndWaitForResult,
  openApp,
  waitForAppShell
} = require('./helpers/app-lifecycle.cjs');

const installDiagramRequestObserver = (page) => page.addInitScript(() => {
  window.__GBDRAW_DIAGRAM_RUNS__ = [];
  const NativeWorker = window.Worker;
  window.Worker = new Proxy(NativeWorker, {
    construct(target, args) {
      const worker = Reflect.construct(target, args, target);
      if (!String(args[0] || '').includes('diagram-generation-worker.js')) return worker;
      const nativePostMessage = worker.postMessage.bind(worker);
      worker.postMessage = (message, transfer) => {
        if (message?.type === 'run' && message?.payload?.request) {
          window.__GBDRAW_DIAGRAM_RUNS__.push(
            JSON.parse(JSON.stringify(message.payload.request))
          );
        }
        if (transfer === undefined) return nativePostMessage(message);
        return nativePostMessage(message, transfer);
      };
      return worker;
    }
  });
});

const ROTATE_RECORD = 'Rotate record using this feature';

// Custom position: "Record starts at <reference> shifted by <offset> bp".
const chooseCustomPosition = async (actions, reference, offset) => {
  await actions.getByRole('radio', { name: 'Custom position', exact: true }).check();
  await actions.getByLabel('Record starts at', { exact: true }).selectOption(reference);
  await actions.getByLabel('Shifted by (bp)', { exact: true }).fill(String(offset));
};

const makeCircularRecord = (recordId, gene, start, end) => `LOCUS       ${recordId.padEnd(24)} 360 bp    DNA     circular UNA 01-JAN-2000
DEFINITION  feature popup record rotation acceptance.
ACCESSION   ${recordId}
VERSION     ${recordId}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
     source          1..360
     CDS             ${start}..${end}
                     /gene="${gene}"
                     /locus_tag="${gene}"
                     /product="${gene} protein"
                     /translation="MKKKKKKKKK"
ORIGIN
        1 ${'acgt'.repeat(15)}
       61 ${'acgt'.repeat(15)}
      121 ${'acgt'.repeat(15)}
      181 ${'acgt'.repeat(15)}
      241 ${'acgt'.repeat(15)}
      301 ${'acgt'.repeat(15)}
//
`;

test('feature popup record rotation works by pointer and keyboard in rich and simple layouts', async ({
  page
}) => {
  test.setTimeout(240000);
  page.on('dialog', (dialog) => dialog.dismiss());
  await page.addInitScript(() => {
    window.__ROTATION_HISTORY_COMMITS__ = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      onHistoryDiagnostic: (event) => {
        if (event.type === 'commit' && event.created) {
          window.__ROTATION_HISTORY_COMMITS__.push({ scope: event.scope, label: event.label });
        }
      }
    };
  });
  await installDiagramRequestObserver(page);
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(join(
    process.cwd(), 'tests/test_inputs/HmmtDNA.gbk'
  ));
  await generateAndWaitForResult(page);
  const before = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return {
      svg: window.__GBDRAW_APP__.svgContent,
      prefix: window.__GBDRAW_DIAGRAM_RUNS__.at(-1).output.prefix,
      request: (await import('/gbdraw/web/js/services/config.js'))
        .getCommittedCanonicalRenderRequest(),
      history: window.__GBDRAW_HISTORY__.getUndoCount(),
      biological: state.biologicalFeatures.value.map((feature) => [
        feature.record_key,
        feature.biological_feature_id
      ])
    };
  });
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.richFeaturePopup = true;
    app.form.prefix = 'UNRELATED_PENDING_PREFIX';
  });

  const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
  await search.fill('tRNA');
  await search.press('Enter');
  await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
  const layout = page.getByRole('region', { name: 'Layout · applies on Generate' });
  await expect(layout.getByLabel('Feature placement', { exact: true })).toBeVisible();
  const disclosure = layout.getByRole('button', { name: ROTATE_RECORD, exact: true });
  await expect(disclosure).toHaveAttribute('aria-expanded', 'false');
  await disclosure.click();
  const actions = page.getByRole('region', { name: ROTATE_RECORD, exact: true });
  await expect(actions).toBeVisible();
  await expect(actions.getByRole('radio', { name: 'Start of the record', exact: true })).toBeChecked();
  await expect(actions.getByRole('button', { name: 'About these coordinates', exact: true })).toBeVisible();
  await expect(actions.locator('[data-record-rotation-preview]'))
    .toContainText(/^NC_012920\.1 will start at [\d,]+ · orientation unchanged/);
  await expect(actions.getByLabel('Record starts at', { exact: true })).toHaveCount(0);

  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount()))
    .toBe(before.history);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.svgContent)).toBe(before.svg);
  await chooseCustomPosition(actions, 'midpoint', 2);
  await disclosure.click();
  await expect(actions).toBeHidden();
  await disclosure.click();
  await expect(actions.getByRole('radio', { name: 'Custom position', exact: true })).toBeChecked();
  await expect(actions.getByLabel('Record starts at', { exact: true })).toHaveValue('midpoint');
  await expect(actions.getByLabel('Shifted by (bp)', { exact: true })).toHaveValue('2');
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount()))
    .toBe(before.history);
  expect(await page.evaluate(async () => (
    (await import('/gbdraw/web/js/services/config.js')).getCommittedCanonicalRenderRequest()
  ))).toEqual(before.request);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.svgContent)).toBe(before.svg);
  const expectedStart = await page.evaluate(() => (
    window.__GBDRAW_APP__.featureRecordRotationDraft.startCoordinate
  ));
  expect(expectedStart).toBeGreaterThan(0);
  await actions.getByRole('button', { name: 'Apply and regenerate' }).click();
  await expect(actions.locator('[aria-live="polite"]')).toContainText('Regenerating');
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, {
    timeout: 240000
  });
  await expect(actions.locator('[aria-live="polite"]')).toContainText(
    'Record rotation applied and regenerated.'
  );
  expect(await search.inputValue()).toBe('tRNA');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.form.prefix))
    .toBe('UNRELATED_PENDING_PREFIX');
  expect(await page.evaluate(() => {
    const controls = window.__GBDRAW_APP__.recordDisplayControls;
    const rows = controls.rows?.value || controls.rows || [];
    const row = rows.find((entry) => entry.recordId === 'NC_012920.1') || rows[0];
    const draft = controls.draftFor(row);
    return {
      startCoordinate: draft.startCoordinate,
      anchorFeatureId: draft.anchorIntent?.biologicalFeatureId || ''
    };
  })).toMatchObject({
    startCoordinate: expectedStart
  });
  const committed = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('/gbdraw/web/js/state.js');
    const request = window.__GBDRAW_DIAGRAM_RUNS__.at(-1);
    return {
      schema: request.schema,
      startCoordinate: request.records[0].display.startCoordinate,
      prefix: request.output.prefix,
      history: window.__GBDRAW_HISTORY__.getUndoCount(),
      biological: state.biologicalFeatures.value.map((feature) => [
        feature.record_key,
        feature.biological_feature_id
      ]),
      target: {
        recordKey: app.featureRecordRotationDraft.identity.recordKey,
        biologicalFeatureId: app.featureRecordRotationDraft.identity.biologicalFeatureId
      }
    };
  });
  expect(committed).toMatchObject({
    schema: CURRENT_REQUEST_SCHEMA,
    startCoordinate: expectedStart,
    prefix: before.prefix,
    history: before.history + 1
  });
  expect(await page.evaluate(() => window.__ROTATION_HISTORY_COMMITS__
    .slice(-1))).toEqual([{ scope: 'artifact-replacement', label: 'Rotate record to feature' }]);
  const rotatedRequest = await page.evaluate(async () => (
    (await import('/gbdraw/web/js/services/config.js')).getCommittedCanonicalRenderRequest()
  ));
  expect(rotatedRequest.records[0].display.startCoordinate).toBe(expectedStart);
  expect(rotatedRequest.output.prefix).toBe(before.prefix);
  expect(committed.biological).toEqual(before.biological);
  expect(committed.target.recordKey).toBeTruthy();
  expect(committed.target.biologicalFeatureId).toBeTruthy();
  const transformedSvg = await page.evaluate(() => window.__GBDRAW_APP__.svgContent);
  expect(transformedSvg).not.toBe(before.svg);

  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.svgContent)).toBe(before.svg);
  expect(await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return state.activeDrawing().recordDisplayDrafts.length;
  })).toBe(0);
  expect(await page.evaluate(async () => (
    (await import('/gbdraw/web/js/services/config.js')).getCommittedCanonicalRenderRequest()
  ))).toEqual(before.request);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.form.prefix))
    .toBe('UNRELATED_PENDING_PREFIX');
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount()))
    .toBe(before.history);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.svgContent)).toBe(transformedSvg);
  expect(await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return state.activeDrawing().recordDisplayDrafts[0];
  })).toMatchObject({
    startCoordinate: expectedStart,
    anchorIntent: {
      schema: 1,
      placement: 'anchor',
      anchor: 'midpoint',
      offsetBp: 2
    }
  });

  expect(await page.evaluate(async () => (
    (await import('/gbdraw/web/js/services/config.js')).getCommittedCanonicalRenderRequest()
  ))).toEqual(rotatedRequest);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.form.prefix))
    .toBe('UNRELATED_PENDING_PREFIX');
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount()))
    .toBe(before.history + 1);
  await page.setViewportSize({ width: 390, height: 844 });
  await page.evaluate(() => { window.__GBDRAW_APP__.richFeaturePopup = false; });
  await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
  const simplePopup = page.locator('.feature-popup--simple');
  await expect(simplePopup).toBeVisible();
  const bounds = await simplePopup.boundingBox();
  expect(bounds.x).toBeGreaterThanOrEqual(0);
  expect(bounds.x + bounds.width).toBeLessThanOrEqual(390);

  const simpleDisclosure = simplePopup.getByRole('button', { name: ROTATE_RECORD, exact: true });
  await expect(simpleDisclosure).toHaveAttribute('aria-expanded', 'false');
  await simpleDisclosure.focus();
  await page.keyboard.press('Enter');
  const simpleActions = simplePopup.getByRole('region', { name: ROTATE_RECORD, exact: true });
  // The position radio group: Tab enters at the checked choice, arrows move it.
  const startChoice = simpleActions.getByRole('radio', { name: 'Start of the record', exact: true });
  await expect(startChoice).toBeChecked();
  await simpleDisclosure.focus();
  await page.keyboard.press('Tab');
  await expect(startChoice).toBeFocused();
  await page.keyboard.press('ArrowDown');
  await expect(simpleActions.getByRole('radio', { name: 'End of the record', exact: true })).toBeChecked();
  await page.keyboard.press('ArrowDown');
  await expect(simpleActions.getByRole('radio', { name: 'Custom position', exact: true })).toBeChecked();
  const reference = simpleActions.getByLabel('Record starts at', { exact: true });
  await page.keyboard.press('Tab');
  await expect(reference).toBeFocused();
  await reference.press('Enter');
  await reference.press('ArrowDown');
  await reference.press('Enter');
  await expect(reference).toHaveValue('midpoint');
  const offset = simpleActions.getByLabel('Shifted by (bp)', { exact: true });
  await page.keyboard.press('Tab');
  await expect(offset).toBeFocused();
  await page.keyboard.press('ControlOrMeta+A');
  await page.keyboard.type('-3');
  await expect(offset).toHaveValue('-3');
  for (const control of [simpleDisclosure, reference, offset,
    simpleActions.getByRole('button', { name: 'Apply on Generate', exact: true }),
    simpleActions.getByRole('button', { name: 'Apply and regenerate', exact: true })]) {
    const box = await control.boundingBox();
    expect(box.x).toBeGreaterThanOrEqual(0);
    expect(box.x + box.width).toBeLessThanOrEqual(390);
  }
  const cancel = simpleActions.getByRole('button', { name: 'Cancel', exact: true });
  await cancel.focus();
  await page.keyboard.press('Enter');
  await expect(simplePopup).toBeVisible();
  await expect(simpleDisclosure).toHaveAttribute('aria-expanded', 'false');
  await expect(simpleActions).toBeHidden();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.svgContent)).toBe(transformedSvg);
  await simpleDisclosure.click();
  const staleActions = simplePopup.getByRole('region', { name: ROTATE_RECORD, exact: true });
  const staleBefore = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const target = window.__GBDRAW_APP__.featureRecordRotationDraft.identity;
    const duplicate = (features) => {
      const feature = features.find((entry) => (
        entry.record_key === target.recordKey
        && entry.biological_feature_id === target.biologicalFeatureId
      ));
      features.push(JSON.parse(JSON.stringify(feature)));
    };
    duplicate(state.extractedFeatures.value);
    duplicate(state.biologicalFeatures.value);
    return {
      svg: window.__GBDRAW_APP__.svgContent,
      history: window.__GBDRAW_HISTORY__.getUndoCount(),
      drafts: JSON.stringify(state.activeDrawing().recordDisplayDrafts)
    };
  });
  await staleActions.getByRole('button', { name: 'Apply and regenerate' }).click();
  // One reason source: the reason line states it once; the status line does not repeat it.
  await expect(staleActions.locator('#feature-record-rotation-reason')).toContainText(
    'no longer present'
  );
  await expect(staleActions.locator('[aria-live="polite"]')).not.toContainText(
    'no longer present'
  );
  expect(await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return {
      svg: window.__GBDRAW_APP__.svgContent,
      history: window.__GBDRAW_HISTORY__.getUndoCount(),
      drafts: JSON.stringify(state.activeDrawing().recordDisplayDrafts)
    };
  })).toEqual(staleBefore);
  await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    state.extractedFeatures.value.pop();
    state.biologicalFeatures.value.pop();
  });
  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  await page.setViewportSize({ width: 1280, height: 900 });

  const pendingSave = page.waitForEvent('download');
  await evaluateWithRetainedPromise(page, async () => {
    window.__GBDRAW_APP__.sessionTitle = 'feature-popup-record-rotation';
    await window.__GBDRAW_APP__.saveSessionWithTitle();
  });
  const savedPath = await (await pendingSave).path();
  const saved = JSON.parse(gunzipSync(fs.readFileSync(savedPath)));
  expect(saved.version).toBe(CURRENT_SESSION_VERSION);
  expect(saved.renderRequest.schema).toBe(CURRENT_REQUEST_SCHEMA);
  expect(saved.renderRequest.records[0].display.startCoordinate).toBe(expectedStart);
  expect(saved.modes[saved.ui.mode].config.recordDisplayDrafts[0].anchorIntent).toMatchObject({
    schema: 1,
    placement: 'anchor',
    anchor: 'midpoint',
    offsetBp: 2
  });
  await page.reload({ waitUntil: 'domcontentloaded' });
  await waitForAppShell(page);
  await page.locator('input[accept^=".json,"]').first().setInputFiles(savedPath);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.sessionImportPending), {
    timeout: 180000
  }).toBe(false);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.svgContent)).toBe(transformedSvg);
  expect(await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return state.activeDrawing().recordDisplayDrafts[0];
  })).toMatchObject({
    startCoordinate: expectedStart,
    anchorIntent: {
      schema: 1,
      recordKey: committed.target.recordKey,
      biologicalFeatureId: committed.target.biologicalFeatureId
    }
  });
});

test('same-file record rotation keeps chromosome targets independent', async ({ page }) => {
  test.setTimeout(240000);
  await installDiagramRequestObserver(page);
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({
    name: 'two-chromosomes.gbk',
    mimeType: 'text/plain',
    buffer: Buffer.from(
      makeCircularRecord('chromosome_I', 'dnaA', 21, 105)
      + makeCircularRecord('chromosome_II', 'parB', 151, 225)
    )
  });
  await generateAndWaitForResult(page);
  const baseline = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return {
      svg: window.__GBDRAW_APP__.svgContent,
      biological: state.biologicalFeatures.value.map((feature) => [
        feature.record_key,
        feature.biological_feature_id
      ])
    };
  });

  const rotateFromPopup = async (query, anchor, offset) => {
    const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
    await search.fill(query);
    await search.press('Enter');
    await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
    await page.getByRole('button', { name: ROTATE_RECORD, exact: true }).click();
    const actions = page.getByRole('region', { name: ROTATE_RECORD, exact: true });
    await expect(actions).toBeVisible();
    await chooseCustomPosition(actions, anchor, offset);
    const draft = await page.evaluate(() => ({
      identity: { ...window.__GBDRAW_APP__.featureRecordRotationDraft.identity },
      startCoordinate: window.__GBDRAW_APP__.featureRecordRotationDraft.startCoordinate
    }));
    await actions.getByRole('button', { name: 'Apply and regenerate' }).click();
    await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, {
      timeout: 240000
    });
    await expect(actions.locator('[aria-live="polite"]')).toContainText(
      'Record rotation applied and regenerated.'
    );
    const snapshot = await page.evaluate(async () => {
      const { state } = await import('/gbdraw/web/js/state.js');
      return {
        request: JSON.parse(JSON.stringify(window.__GBDRAW_DIAGRAM_RUNS__.at(-1))),
        svg: window.__GBDRAW_APP__.svgContent,
        biological: state.biologicalFeatures.value.map((feature) => [
          feature.record_key,
          feature.biological_feature_id
        ])
      };
    });
    await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
    return { draft, snapshot };
  };

  const chromosomeI = await rotateFromPopup('dnaA', 'five-prime', -5);
  expect(chromosomeI.snapshot.request.records).toHaveLength(2);
  const firstIndex = chromosomeI.snapshot.request.records.findIndex((record) => (
    record.recordKey === chromosomeI.draft.identity.recordKey
  ));
  expect(firstIndex).toBeGreaterThanOrEqual(0);
  const secondIndex = firstIndex === 0 ? 1 : 0;
  expect(chromosomeI.snapshot.request.records[firstIndex].display.startCoordinate)
    .toBe(chromosomeI.draft.startCoordinate);
  expect(chromosomeI.snapshot.request.records[secondIndex].display.startCoordinate).toBeNull();
  expect(chromosomeI.snapshot.svg).not.toBe(baseline.svg);
  expect(chromosomeI.snapshot.biological).toEqual(baseline.biological);

  const preservedFirstRecord = chromosomeI.snapshot.request.records[firstIndex];
  const chromosomeII = await rotateFromPopup('parB', 'midpoint', 7);
  const secondTargetIndex = chromosomeII.snapshot.request.records.findIndex((record) => (
    record.recordKey === chromosomeII.draft.identity.recordKey
  ));
  expect(secondTargetIndex).toBe(secondIndex);
  expect(chromosomeII.snapshot.request.records[firstIndex]).toEqual(preservedFirstRecord);
  expect(chromosomeII.snapshot.request.records[secondIndex].display.startCoordinate)
    .toBe(chromosomeII.draft.startCoordinate);
  expect(chromosomeII.snapshot.request.records.map((record) => record.presentation.gridRow))
    .toEqual(chromosomeI.snapshot.request.records.map((record) => record.presentation.gridRow));
  expect(chromosomeII.snapshot.request.tracks).toEqual(chromosomeI.snapshot.request.tracks);
  expect(chromosomeII.snapshot.request.comparisons)
    .toEqual(chromosomeI.snapshot.request.comparisons);
  expect(chromosomeII.snapshot.svg).not.toBe(chromosomeI.snapshot.svg);
  expect(chromosomeII.snapshot.biological).toEqual(baseline.biological);
});

// PD-OI-085: Apply on Generate stages a record rotation in the record display
// draft; one Generate Diagram applies every staged record, and Apply and
// regenerate stays target-only (PD-OI-032 item 4).
test('Apply on Generate stages rotations that one Generate applies together', async ({ page }) => {
  test.setTimeout(300000);
  await page.addInitScript(() => {
    window.__ROTATION_HISTORY_COMMITS__ = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      onHistoryDiagnostic: (event) => {
        if (event.type === 'commit' && event.created) {
          window.__ROTATION_HISTORY_COMMITS__.push({ scope: event.scope, label: event.label });
        }
      }
    };
  });
  await installDiagramRequestObserver(page);
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({
    name: 'two-chromosomes.gbk',
    mimeType: 'text/plain',
    buffer: Buffer.from(
      makeCircularRecord('chromosome_I', 'dnaA', 21, 105)
      + makeCircularRecord('chromosome_II', 'parB', 151, 225)
    )
  });
  await generateAndWaitForResult(page);
  const snapshot = () => page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return {
      runs: window.__GBDRAW_DIAGRAM_RUNS__.length,
      svg: window.__GBDRAW_APP__.svgContent,
      history: window.__GBDRAW_HISTORY__.getUndoCount(),
      drafts: Object.fromEntries(state.activeDrawing().recordDisplayDrafts
        .map(({ recordId, startCoordinate }) => [recordId, startCoordinate]))
    };
  });
  const lastStarts = () => page.evaluate(() => window.__GBDRAW_DIAGRAM_RUNS__.at(-1).records
    .map((record) => record.display.startCoordinate)
    .sort((a, b) => a - b));
  const pendingNotice = page.getByText(
    'Record rotation, feature placement, or tolerance has changes pending Generate.'
  );
  const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
  const openRotation = async (query) => {
    await search.fill(query);
    await search.press('Enter');
    await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
    const disclosure = page.getByRole('button', { name: ROTATE_RECORD, exact: true });
    await expect(disclosure).toHaveAttribute('aria-expanded', 'false');
    await disclosure.click();
    const actions = page.getByRole('region', { name: ROTATE_RECORD, exact: true });
    await expect(actions).toBeVisible();
    return actions;
  };
  const closePopup = () => page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  const before = await snapshot();
  const initialStarts = await lastStarts();
  expect(initialStarts).toEqual([null, null]);
  await expect(pendingNotice).toHaveCount(0);

  // Stage chromosome I at the start of dnaA (21) and chromosome II after parB (226).
  let actions = await openRotation('dnaA');
  await expect(actions.locator('[data-record-rotation-preview]'))
    .toHaveText('chromosome_I will start at 21 · orientation unchanged');
  await actions.getByRole('button', { name: 'Apply on Generate', exact: true }).click();
  await expect(actions.locator('[aria-live="polite"]'))
    .toHaveText('Record rotation will apply on the next Generate Diagram.');
  await expect(actions.locator('[data-record-rotation-pending]'))
    .toHaveText('Pending for Generate: chromosome_I will start at 21 · orientation unchanged');
  await expect(actions.locator('[data-record-rotation-preview]')).toHaveCount(0);
  await closePopup();
  actions = await openRotation('parB');
  await expect(actions.locator('[data-record-rotation-pending]')).toHaveCount(0);
  await actions.getByRole('radio', { name: 'End of the record', exact: true }).check();
  await actions.getByRole('button', { name: 'Apply on Generate', exact: true }).click();
  await expect(actions.locator('[data-record-rotation-pending]'))
    .toHaveText('Pending for Generate: chromosome_II will start at 226 · orientation unchanged');
  const staged = await snapshot();
  expect(staged).toMatchObject({ runs: before.runs, svg: before.svg, history: before.history + 2 });
  expect(staged.drafts).toEqual({ chromosome_I: 21, chromosome_II: 226 });
  expect(await page.evaluate(() => window.__ROTATION_HISTORY_COMMITS__.slice(-2))).toEqual([
    { scope: 'intent', label: 'Rotate record to feature on Generate' },
    { scope: 'intent', label: 'Rotate record to feature on Generate' }
  ]);
  await expect(pendingNotice).toBeVisible();

  // One Undo removes only the last staged record; Redo stages it again.
  // History restore closes the popup; reopening shows the current draft.
  await actions.getByRole('button', { name: 'Apply on Generate', exact: true }).focus();
  await page.keyboard.press('ControlOrMeta+z');
  await expect.poll(async () => (await snapshot()).drafts).toEqual({ chromosome_I: 21 });
  actions = await openRotation('parB');
  await expect(actions.locator('[data-record-rotation-pending]')).toHaveCount(0);
  await expect(actions.locator('[data-record-rotation-preview]'))
    .toHaveText('chromosome_II will start at 151 · orientation unchanged');
  await page.keyboard.press('ControlOrMeta+Shift+z');
  await expect.poll(async () => (await snapshot()).drafts).toEqual(staged.drafts);
  expect(await snapshot()).toEqual(staged);
  await page.keyboard.press('Escape');
  await expect(actions).toBeHidden();

  // Reopening shows the staged values.
  for (const [query, expected] of [
    ['dnaA', 'Pending for Generate: chromosome_I will start at 21 · orientation unchanged'],
    ['parB', 'Pending for Generate: chromosome_II will start at 226 · orientation unchanged']
  ]) {
    actions = await openRotation(query);
    await expect(actions.locator('[data-record-rotation-pending]')).toHaveText(expected);
    await closePopup();
  }

  // One Generate Diagram applies both staged records.
  await generateAndWaitForResult(page);
  const generated = await snapshot();
  expect(generated.runs).toBe(before.runs + 1);
  expect(generated.svg).not.toBe(before.svg);
  expect(await lastStarts()).toEqual([21, 226]);
  await expect(pendingNotice).toHaveCount(0);
  actions = await openRotation('dnaA');
  await expect(actions.locator('[data-record-rotation-pending]')).toHaveCount(0);
  await closePopup();

  // Stage chromosome II again, then Apply and regenerate chromosome I only:
  // the staged chromosome II stays staged and is not drawn.
  actions = await openRotation('parB');
  await actions.getByRole('button', { name: 'Apply on Generate', exact: true }).click();
  await expect(actions.locator('[data-record-rotation-pending]'))
    .toHaveText('Pending for Generate: chromosome_II will start at 151 · orientation unchanged');
  await closePopup();
  actions = await openRotation('dnaA');
  await actions.getByRole('radio', { name: 'End of the record', exact: true }).check();
  await actions.getByRole('button', { name: 'Apply and regenerate', exact: true }).click();
  await expect(actions.locator('[aria-live="polite"]')).toHaveText(
    'Record rotation applied and regenerated.', { timeout: 240000 }
  );
  expect(await lastStarts()).toEqual([106, 226]);
  expect((await snapshot()).drafts).toEqual({ chromosome_I: 106, chromosome_II: 151 });
  await expect(pendingNotice).toBeVisible();
  await closePopup();
  actions = await openRotation('parB');
  await expect(actions.locator('[data-record-rotation-pending]'))
    .toHaveText('Pending for Generate: chromosome_II will start at 151 · orientation unchanged');
});

test('both modes record rotation resolves the same circular source anchor', async ({ browser }) => {
  test.setTimeout(300000);
  const source = makeCircularRecord('shared_anchor', 'anchor_gene', 41, 125);
  const rotateInMode = async (mode) => {
    const context = await browser.newContext();
    const page = await context.newPage();
    await installDiagramRequestObserver(page);
    await openApp(page);
    if (mode === 'linear') {
      await page.getByRole('button', { name: 'Linear', exact: true }).click();
    }
    const upload = mode === 'linear'
      ? page.getByTestId('linear-genbank-1')
      : page.getByLabel('GenBank/DDBJ File', { exact: true });
    await upload.setInputFiles({
      name: 'shared-anchor.gbk',
      mimeType: 'text/plain',
      buffer: Buffer.from(source)
    });
    await generateAndWaitForResult(page);
    const beforeSvg = await page.evaluate(() => window.__GBDRAW_APP__.svgContent);
    const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
    await search.fill('anchor_gene');
    await search.press('Enter');
    await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
    await page.getByRole('button', { name: ROTATE_RECORD, exact: true }).click();
    const actions = page.getByRole('region', { name: ROTATE_RECORD, exact: true });
    await chooseCustomPosition(actions, 'midpoint', -11);
    const expected = await page.evaluate(() => ({
      identity: { ...window.__GBDRAW_APP__.featureRecordRotationDraft.identity },
      startCoordinate: window.__GBDRAW_APP__.featureRecordRotationDraft.startCoordinate
    }));
    await actions.getByRole('button', { name: 'Apply and regenerate' }).click();
    await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, {
      timeout: 240000
    });
    await expect(actions.locator('[aria-live="polite"]')).toContainText(
      'Record rotation applied and regenerated.'
    );
    const accepted = await page.evaluate(async (recordKey) => {
      const { state } = await import('/gbdraw/web/js/state.js');
      const request = window.__GBDRAW_DIAGRAM_RUNS__.at(-1);
      const record = request.records.find((entry) => entry.recordKey === recordKey);
      const draft = state.activeDrawing().recordDisplayDrafts.find((entry) => (
        entry.anchorIntent?.recordKey === recordKey
      ));
      return {
        schema: request.schema,
        startCoordinate: record?.display?.startCoordinate,
        anchorIntent: draft?.anchorIntent,
        svg: window.__GBDRAW_APP__.svgContent
      };
    }, expected.identity.recordKey);
    expect(accepted.schema).toBe(CURRENT_REQUEST_SCHEMA);
    expect(accepted.startCoordinate).toBe(expected.startCoordinate);
    expect(accepted.anchorIntent).toMatchObject({
      schema: 1,
      recordKey: expected.identity.recordKey,
      biologicalFeatureId: expected.identity.biologicalFeatureId,
      anchor: 'midpoint',
      offsetBp: -11
    });
    expect(accepted.svg).not.toBe(beforeSvg);
    await context.close();
    return accepted;
  };

  const circular = await rotateInMode('circular');
  const linear = await rotateInMode('linear');
  expect(linear.startCoordinate).toBe(circular.startCoordinate);
  expect(linear.anchorIntent.anchor).toBe(circular.anchorIntent.anchor);
  expect(linear.anchorIntent.offsetBp).toBe(circular.anchorIntent.offsetBp);
});

// Session Load reads no record bytes (776a2f93). Opening Record actions is the
// explicit record action that reads them, so no Generate is needed.
const installRecordReadCounter = (page) => page.addInitScript(() => {
  window.__GBDRAW_RECORD_RESOURCE_READS__ = 0;
  window.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric(metric) {
      if (metric.name === 'resourceByteReadCount') {
        window.__GBDRAW_RECORD_RESOURCE_READS__ += metric.value;
      }
    }
  };
});

const loadSessionFile = async (page, sessionPath) => {
  const dialogPromise = page.waitForEvent('dialog');
  await page.locator('input[accept^=".json,"]').first().setInputFiles(sessionPath);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.sessionImportPending))
    .toBe(false);
};

const openLoadedRecordActions = async (page, { query, recordId, reads = 1 }) => {
  const recordReads = () => page.evaluate(() => window.__GBDRAW_RECORD_RESOURCE_READS__);
  const undoCount = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  expect(await recordReads()).toBe(0);
  await page.evaluate(() => { window.__GBDRAW_APP__.richFeaturePopup = true; });
  const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
  await search.fill(query);
  await search.press('Enter');
  await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
  const disclosure = page.getByRole('button', { name: ROTATE_RECORD, exact: true });
  await expect(disclosure).toHaveAttribute('aria-expanded', 'false');
  expect(await recordReads()).toBe(0);

  await disclosure.click();
  const actions = page.getByRole('region', { name: ROTATE_RECORD, exact: true });
  // Start of the record (default) on a forward record: the feature's leftmost
  // source base becomes display base 1, whatever its strand.
  const expectedStart = await page.evaluate(() => {
    const parts = window.__GBDRAW_APP__.featureRecordRotationDraft.feature.location_parts;
    return parts.length === 1 ? parts[0].start + 1 : null;
  });
  expect(expectedStart).toBeGreaterThan(0);
  await expect(actions.locator('[data-record-rotation-preview]')).toHaveText(
    `${recordId} will start at ${expectedStart.toLocaleString('en-US')} · orientation unchanged`
  );
  await expect(actions.getByRole('button', { name: 'Apply and regenerate' })).toBeEnabled();
  await expect(actions.getByRole('radio', { name: 'End of the record', exact: true })).toBeEnabled();
  await expect(actions).not.toContainText('stale or ambiguous');
  await expect(actions).not.toContainText('Unavailable');
  expect(await recordReads()).toBe(reads);
  expect(await page.evaluate(() => window.__GBDRAW_DIAGRAM_RUNS__.length)).toBe(0);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(undoCount);
  return { actions, expectedStart };
};

test('Circular Record actions read a loaded Gallery Session without Generate', async ({ page }) => {
  test.setTimeout(180000);
  await installRecordReadCounter(page);
  await installDiagramRequestObserver(page);
  await openApp(page);
  await loadSessionFile(page, join(
    process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'
  ));
  await openLoadedRecordActions(page, { query: 'tRNA', recordId: 'NC_012920.1' });
});

test('Linear Record actions rotate a circular record of a loaded Session without Generate', async ({
  page
}) => {
  test.setTimeout(240000);
  page.on('dialog', (dialog) => {
    if (dialog.message() !== 'Session loaded successfully!') dialog.dismiss();
  });
  await installRecordReadCounter(page);
  await installDiagramRequestObserver(page);
  await openApp(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.getByTestId('linear-genbank-1').setInputFiles({
    name: 'loaded-linear.gbk',
    mimeType: 'text/plain',
    buffer: Buffer.from(makeCircularRecord('loaded_linear', 'loaded_gene', 41, 125))
  });
  await generateAndWaitForResult(page);
  const pendingSave = page.waitForEvent('download');
  await evaluateWithRetainedPromise(page, async () => {
    window.__GBDRAW_APP__.sessionTitle = 'record-actions-after-load';
    await window.__GBDRAW_APP__.saveSessionWithTitle();
  });
  const savedPath = await (await pendingSave).path();
  await page.reload({ waitUntil: 'domcontentloaded' });
  await waitForAppShell(page);
  await loadSessionFile(page, savedPath);

  const { actions, expectedStart } = await openLoadedRecordActions(page, {
    query: 'loaded_gene', recordId: 'loaded_linear'
  });
  await actions.getByRole('button', { name: 'Apply and regenerate' }).click();
  await expect(actions.locator('[aria-live="polite"]')).toContainText(
    'Record rotation applied and regenerated.', { timeout: 240000 }
  );
  const request = await page.evaluate(() => window.__GBDRAW_DIAGRAM_RUNS__.at(-1));
  expect(request.mode).toBe('linear');
  expect(request.records.map((record) => record.display.startCoordinate)).toEqual([expectedStart]);
});

// `gbdraw ... --session_output` binds its input files for the Web draft, and
// its request draws each record from those same files, with its crop and
// orientation.
const writeCliSession = async (testInfo, mode, inputs, args = []) => {
  const paths = Object.entries(inputs).map(([name, text]) => {
    const path = testInfo.outputPath(name);
    fs.writeFileSync(path, text);
    return path;
  });
  const prefix = testInfo.outputPath(`cli-${mode}`);
  const session = `${prefix}.gbdraw-session.json`;
  await promisify(execFile)('python', [
    '-m', 'gbdraw.cli', mode, '--gbk', ...paths, ...args, '-o', prefix, '-f', 'svg',
    '--session_output', session
  ], { cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: process.cwd() }, timeout: 300000 });
  return session;
};

test('Circular Record actions read a loaded CLI Session without Generate', async ({
  page
}, testInfo) => {
  test.setTimeout(240000);
  const session = await writeCliSession(testInfo, 'circular', {
    'cli-circular.gbk': makeCircularRecord('cli_circular', 'cli_gene', 41, 125)
  });
  await installRecordReadCounter(page);
  await installDiagramRequestObserver(page);
  await openApp(page);
  await loadSessionFile(page, session);
  await openLoadedRecordActions(page, { query: 'cli_gene', recordId: 'cli_circular' });
});

test('Linear Record actions rotate a multi-record File of a loaded CLI Session without Generate', async ({
  page
}, testInfo) => {
  test.setTimeout(300000);
  page.on('dialog', (dialog) => {
    if (dialog.message() !== 'Session loaded successfully!') dialog.dismiss();
  });
  const session = await writeCliSession(testInfo, 'linear', {
    'cli-single.gbk': makeCircularRecord('cli_single', 'single_gene', 41, 125),
    'cli-multi.gbk': makeCircularRecord('cli_multi_a', 'multi_a_gene', 61, 150)
      + makeCircularRecord('cli_multi_b', 'multi_b_gene', 101, 200)
  });
  await installRecordReadCounter(page);
  await installDiagramRequestObserver(page);
  await openApp(page);
  await loadSessionFile(page, session);

  const { actions, expectedStart } = await openLoadedRecordActions(page, {
    query: 'multi_b_gene', recordId: 'cli_multi_b', reads: 2
  });
  await actions.getByRole('button', { name: 'Apply and regenerate' }).click();
  await expect(actions.locator('[aria-live="polite"]')).toContainText(
    'Record rotation applied and regenerated.', { timeout: 240000 }
  );
  const request = await page.evaluate(() => window.__GBDRAW_DIAGRAM_RUNS__.at(-1));
  expect(request.mode).toBe('linear');
  expect(request.records.map((record) => record.display.startCoordinate))
    .toEqual([null, null, expectedStart]);
});

// One original file stays one File holding its records; record selection
// decides what is drawn. A CLI Session that draws some records of a file
// selects them on the bound file.
test('Linear CLI Session that draws some records of a file loads and generates only those records', async ({
  page
}, testInfo) => {
  test.setTimeout(300000);
  page.on('dialog', (dialog) => {
    if (dialog.message() !== 'Session loaded successfully!') dialog.dismiss();
  });
  const session = await writeCliSession(testInfo, 'linear', {
    'cli-single.gbk': makeCircularRecord('cli_single', 'single_gene', 41, 125),
    'cli-multi.gbk': makeCircularRecord('cli_multi_a', 'multi_a_gene', 61, 150)
      + makeCircularRecord('cli_multi_b', 'multi_b_gene', 101, 200)
  }, ['--record_id', '', '--record_id', 'cli_multi_b']);
  await installRecordReadCounter(page);
  await installDiagramRequestObserver(page);
  await openApp(page);
  await loadSessionFile(page, session);
  const linearRows = () => page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return state.linearSeqs.map((seq) => [seq.uid, seq.gb?.name, seq.region_record_id]);
  });
  await expect(page.locator('[data-linear-source-records]')).toHaveCount(2);
  expect(await linearRows()).toEqual([
    ['record-1', 'cli-single.gbk', ''],
    ['record-2', 'cli-multi.gbk', 'cli_multi_b']
  ]);

  await openLoadedRecordActions(page, { query: 'multi_b_gene', recordId: 'cli_multi_b', reads: 2 });
  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();

  await generateAndWaitForResult(page);
  const request = await page.evaluate(() => window.__GBDRAW_DIAGRAM_RUNS__.at(-1));
  expect(request.records.map((record) => [record.recordKey, record.selector])).toEqual([
    ['record-1', null],
    ['record-2', { kind: 'recordId', value: 'cli_multi_b' }]
  ]);
  expect(await linearRows()).toHaveLength(2);
  const drawnRecordIds = (content) => [...new Set(
    [...content.matchAll(/data-gbdraw-record-id="([^"]+)"/g)].map((match) => match[1])
  )];
  const cliSvg = fs.readFileSync(testInfo.outputPath('cli-linear.svg'), 'utf8');
  expect(drawnRecordIds(cliSvg)).toEqual(['cli_single', 'cli_multi_b']);
  expect(drawnRecordIds(await page.evaluate(() => window.__GBDRAW_APP__.svgContent)))
    .toEqual(drawnRecordIds(cliSvg));
});

test('Linear CLI Session keeps --region and --reverse_complement through Load and Generate', async ({
  page
}, testInfo) => {
  test.setTimeout(300000);
  page.on('dialog', (dialog) => {
    if (dialog.message() !== 'Session loaded successfully!') dialog.dismiss();
  });
  const session = await writeCliSession(testInfo, 'linear', {
    'cli-cropped.gbk': makeCircularRecord('cli_cropped', 'cropped_gene', 61, 150),
    'cli-reversed.gbk': makeCircularRecord('cli_reversed', 'reversed_gene', 41, 125)
  }, ['--region', 'cli_cropped:21-300', '--reverse_complement', '0', '--reverse_complement', '1']);
  await installDiagramRequestObserver(page);
  await openApp(page);
  await loadSessionFile(page, session);
  await page.evaluate(() => { window.__GBDRAW_APP__.richFeaturePopup = true; });
  const openRecordActions = async (query) => {
    const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
    await search.fill(query);
    await search.press('Enter');
    await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
    await page.getByRole('button', { name: ROTATE_RECORD, exact: true }).click();
    return page.getByRole('region', { name: ROTATE_RECORD, exact: true });
  };
  const closePopup = () => page.getByRole('button', { name: 'Close feature popup', exact: true }).click();

  // The reversed record is its bound input file reverse-complemented, so its
  // rotation is available without Generate.
  const reversed = await openRecordActions('reversed_gene');
  await expect(reversed.locator('[data-record-rotation-preview]'))
    .toContainText('cli_reversed will start at');
  await expect(reversed.getByRole('button', { name: 'Apply and regenerate' })).toBeEnabled();
  await closePopup();
  // A cropped record cannot rotate, and says so (PD-OI-032 item 7).
  const cropped = await openRecordActions('cropped_gene');
  await expect(cropped.locator('#feature-record-rotation-reason'))
    .toHaveText('Record rotation is unavailable for a cropped record.');
  await expect(cropped.getByRole('button', { name: 'Apply and regenerate' })).toBeDisabled();
  await closePopup();
  expect(await page.evaluate(() => window.__GBDRAW_DIAGRAM_RUNS__.length)).toBe(0);

  await generateAndWaitForResult(page);
  const request = await page.evaluate(() => window.__GBDRAW_DIAGRAM_RUNS__.at(-1));
  expect(request.records.map((record) => [record.region, record.presentation.reverseComplement]))
    .toEqual([
      [{ selector: null, start: 21, end: 300, reverseComplement: false }, false],
      [null, true]
    ]);
  // Generate draws what the CLI drew, with its source coordinates.
  const semantics = (content) => {
    const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
    const rounded = (value) => String(value || '').replace(/-?\d+\.\d+/g, (number) => Number(number).toFixed(3));
    return {
      records: [...svg.querySelectorAll('[data-gbdraw-record-source-start]')].map((element) => [
        'record-id', 'record-source-start', 'record-source-end', 'record-source-step'
      ].map((name) => element.getAttribute(`data-gbdraw-${name}`))),
      features: [...svg.querySelectorAll('[data-gbdraw-feature-id]')].map((element) => [
        element.getAttribute('data-gbdraw-feature-id'), rounded(element.getAttribute('d'))
      ]),
      text: [...svg.querySelectorAll('text')].map((element) => element.textContent)
    };
  };
  const cliSvg = fs.readFileSync(testInfo.outputPath('cli-linear.svg'), 'utf8');
  const generated = await page.evaluate(semantics, await page.evaluate(() => window.__GBDRAW_APP__.svgContent));
  expect(generated).toEqual(await page.evaluate(semantics, cliSvg));
  expect(generated.records.map((record) => record.slice(1))).toEqual([
    ['21', '300', '1'],
    ['1', '360', '-1']
  ]);
});

test('browser export embeds the exact selected current-schema item and expands references', async ({
  page
}, testInfo) => {
  await page.goto('/');
  const origin = new URL(page.url()).origin;
  const exported = await page.evaluate(async ({ origin, catalogSchema }) => {
    const { enrichSvgWithStandaloneInteractivity } = await import(
      `${origin}/gbdraw/web/js/services/standalone-interactivity.js`
    );
    const fullNote = `${'x'.repeat(49)}😀tail`;
    const catalog = {
      schema: catalogSchema,
      items: [{
        resultIndex: 0,
        resultName: 'diagram.svg',
        recordKeys: ['record-key-a', 'record-key-b'],
        features: [{
          svgId: 'rendered-visible',
          recordKey: 'record-key-a',
          biologicalFeatureId: 'stable-visible',
          fillColor: '#54bcf8'
        }, {
          svgId: 'rendered-collision',
          recordKey: 'record-key-a',
          biologicalFeatureId: 'h_aaaaaaaaaaaaaaaaaaaaaaaaaa',
          fillColor: '#f59e0b'
        }],
        biologicalFeatures: [
          {
            recordKey: 'record-key-a',
            biologicalFeatureId: 'stable-visible',
            record_id: 'rec-visible',
            type: 'CDS',
            start: 0,
            end: 9,
            strand: '+',
            aminoAcidSequence: 'MVISIBLE',
            translationFromAminoAcidSequence: true,
            qualifiers: {
              locus_tag: ['VP_1'],
              protein_id: ['VP_1'],
              old_locus_tag: ['OLD_VP_1'],
              gene: ['visible_gene'],
              product: ['Visible protein'],
              note: [fullNote]
            },
            sequenceSourceIndex: 0
          },
          {
            recordKey: 'record-key-a',
            biologicalFeatureId: 'h_aaaaaaaaaaaaaaaaaaaaaaaaaa',
            type: 'CDS',
            start: 20,
            end: 29,
            strand: '+',
            qualifiers: {
              note: [fullNote]
            },
            sequenceSourceIndex: 0
          },
          {
            recordKey: 'record-key-b',
            biologicalFeatureId: 'stable-hidden',
            record_id: 'rec-hidden',
            type: 'CDS',
            start: 9,
            end: 18,
            strand: '+',
            location_parts: [
              { start: 9, end: 12, strand: '+' },
              { start: 12, end: 18, strand: '+' }
            ],
            product: 'Hidden override',
            amino_acid_sequence: 'MHIDDEN',
            translationFromAminoAcidSequence: true,
            qualifiers: {
              locus_tag: ['HP_1'],
              protein_id: ['HP_1'],
              product: ['Hidden protein']
            },
            sequenceSourceIndex: 1
          }
        ],
        orthogroups: [{
          id: 'og-hidden',
          name: 'hidden-test',
          description: 'Original group description',
          member_count: 2,
          record_coverage_count: 2,
          members: [
            {
              recordKey: 'record-key-a',
              biologicalFeatureId: 'stable-visible',
              representative: true
            },
            {
              recordKey: 'record-key-b',
              biologicalFeatureId: 'stable-hidden'
            }
          ]
        }],
        annotations: [{
          dom_id: 'annotation-review-window',
          id: 'review-window',
          set_id: 'review',
          track_id: 'annotations-1',
          record_id: 'rec-visible',
          record_index: 0,
          segments: [[2, 8]],
          label: 'Review window',
          mark: 'band',
          lane: 0,
          metadata: { reviewer: 'Ada' }
        }],
        comparisonMatches: [],
        sequenceSources: [{
          key: 'linear:record:0',
          origin: 'linear-record',
          recordIndex: 0,
          sequence: `ATGAAATAA${'N'.repeat(11)}ATGCCCTAA`
        }, {
          key: 'linear:record:1',
          origin: 'linear-record',
          recordIndex: 1,
          sequence: `${'N'.repeat(9)}ATGCCCTAA`
        }]
      }, {
        resultIndex: 1,
        resultName: 'other.svg',
        recordKeys: [],
        features: [],
        biologicalFeatures: [],
        orthogroups: [],
        annotations: [],
        comparisonMatches: []
      }]
    };
    for (let index = 0; index < 128; index += 1) {
      catalog.items[0].biologicalFeatures.push({
        recordKey: 'record-key-a',
        biologicalFeatureId: `bulk-feature-${index}`,
        type: 'CDS',
        start: 0,
        end: 3,
        strand: '+',
        qualifiers: {},
        sequenceSourceIndex: 0
      });
    }
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    svg.setAttribute('viewBox', '0 0 120 80');
    svg.innerHTML = `
      <rect id="rendered-visible" data-gbdraw-feature-id="shared-stable"
        data-gbdraw-rendered-feature-id="rendered-visible"
        x="5" y="5" width="25" height="12" fill="#54bcf8" />
      <rect id="rendered-collision" data-gbdraw-feature-id="shared-stable"
        data-gbdraw-rendered-feature-id="rendered-collision"
        x="35" y="5" width="25" height="12" fill="#f59e0b" />
      <g id="annotation-review-window" data-gbdraw-annotation-id="review-window"
        data-gbdraw-annotation-set-id="review"
        data-gbdraw-annotation-track-id="annotations-1"
        data-gbdraw-record-id="rec-visible" data-gbdraw-record-index="0"
        data-gbdraw-annotation-mark="band" data-gbdraw-annotation-label="Review window">
        <rect x="5" y="25" width="20" height="6" fill="#94a3b8" />
      </g>`;
    const enriched = enrichSvgWithStandaloneInteractivity(svg, {
      popupMode: 'rich',
      featureCatalog: catalog,
      catalogResultIndex: 0,
      catalogResultName: 'diagram.svg',
      requireFeatureCatalog: true,
      labelTextFeatureOverrides: {
        'rendered-visible': 'Edited visible label'
      },
      orthogroupNameOverrides: {
        'og-hidden': 'Edited similarity group'
      },
      orthogroupDescriptionOverrides: {
        'og-hidden': 'Edited group description'
      }
    });
    const metadata = svg.querySelector('#gbdraw-interactive-feature-metadata');
    return {
      enriched,
      catalog,
      embedded: JSON.parse(metadata.textContent),
      schema: metadata.getAttribute('data-schema'),
      resultIndex: metadata.getAttribute('data-result-index'),
      resultName: metadata.getAttribute('data-result-name'),
      sourceDisplayLabel: catalog.items[0].features[0].displayLabel,
      sourceGroup: catalog.items[0].orthogroups[0],
      svgText: new XMLSerializer().serializeToString(svg)
    };
  }, { origin, catalogSchema: CURRENT_FEATURE_CATALOG_SCHEMA });

  expect(exported.enriched).toBe(true);
  expect(exported.embedded.schema).toBe(CURRENT_FEATURE_CATALOG_SCHEMA);
  expect(exported.embedded.items).toHaveLength(1);
  expect(exported.embedded.items[0].features[0].displayLabel)
    .toBe('Edited visible label');
  expect(exported.embedded.items[0].orthogroups[0].display_name)
    .toBe('Edited similarity group');
  expect(exported.embedded.items[0].orthogroups[0].description)
    .toBe('Edited group description');
  expect(exported.embedded.items[0].annotations[0].id).toBe('review-window');
  expect(exported.sourceDisplayLabel).toBeUndefined();
  expect(exported.sourceGroup.display_name).toBeUndefined();
  expect(exported.sourceGroup.description).toBe('Original group description');
  expect(exported.schema).toBe(String(CURRENT_FEATURE_CATALOG_SCHEMA));
  expect(exported.resultIndex).toBe('0');
  expect(exported.resultName).toBe('diagram.svg');
  expect(exported.svgText.match(/data-gbdraw-interactive-feature="true"/g)).toHaveLength(2);

  const svgPath = testInfo.outputPath('interactive-v3.svg');
  fs.writeFileSync(svgPath, exported.svgText, 'utf8');
  await page.addInitScript(() => {
    window.__copiedText = '';
    window.__expandedCatalogFeatures = {};
    window.__sourceValidationScans = 0;
    window.__sourceValidationScanCounts = {};
    const sourceSequences = new Set([
      `ATGAAATAA${'N'.repeat(11)}ATGCCCTAA`,
      `${'N'.repeat(9)}ATGCCCTAA`
    ]);
    const nativeRegexTest = RegExp.prototype.test;
    RegExp.prototype.test = function (value) {
      if (this.source === '\\s' && sourceSequences.has(value)) {
        window.__sourceValidationScans += 1;
        window.__sourceValidationScanCounts[value] = (
          window.__sourceValidationScanCounts[value] || 0
        ) + 1;
      }
      return nativeRegexTest.call(this, value);
    };
    const nativeMapSet = Map.prototype.set;
    Map.prototype.set = function (key, value) {
      const biologicalFeatureId = String(
        value && (
          value.biologicalFeatureId || value.biological_feature_id
        ) || ''
      );
      if (
        biologicalFeatureId
        && value.qualifiers
        && Array.isArray(value.location_parts)
      ) {
        window.__expandedCatalogFeatures[biologicalFeatureId] = {
          aminoAcidSequence: (
            value.amino_acid_sequence || value.aminoAcidSequence
          ),
          gene: value.gene,
          locationParts: value.location_parts,
          locusTag: value.locus_tag,
          note: value.note,
          oldLocusTag: value.old_locus_tag,
          product: value.product,
          proteinId: value.protein_id,
          translation: value.qualifiers && value.qualifiers.translation,
          translationMarker: value.translationFromAminoAcidSequence
        };
      }
      return nativeMapSet.call(this, key, value);
    };
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {
        writeText: async (value) => {
          window.__copiedText = String(value);
        }
      }
    });
  });
  await page.goto(pathToFileURL(svgPath).href);
  await page.clock.install();

  const expandedFeatures = await page.evaluate(
    () => window.__expandedCatalogFeatures
  );
  expect(await page.evaluate(() => window.__sourceValidationScans)).toBe(2);
  expect(await page.evaluate(() => Object.values(
    window.__sourceValidationScanCounts
  ))).toEqual([1, 1]);
  expect(expandedFeatures['stable-visible']).toMatchObject({
    aminoAcidSequence: 'MVISIBLE',
    gene: 'visible_gene',
    locusTag: 'VP_1',
    oldLocusTag: 'OLD_VP_1',
    product: 'Visible protein',
    proteinId: 'VP_1',
    note: `${'x'.repeat(49)}😀`,
    translation: ['MVISIBLE'],
    locationParts: [{
      start: 0,
      end: 9,
      strand: '+',
      display: '1..9'
    }]
  });
  expect(expandedFeatures['stable-visible'].translationMarker).toBeUndefined();
  expect(expandedFeatures['stable-hidden'].product).toBe('Hidden override');

  await page.locator('[data-gbdraw-rendered-feature-id="rendered-visible"]').click();
  await expect(page.locator('.gfi-title')).toContainText('Edited visible label');
  await expect(page.locator('#gbdraw-feature-popup')).toContainText(
    'Edited similarity group'
  );
  const memberBlock = page.locator('.gfi-block').filter({
    hasText: 'Similarity-group members'
  }).last();
  await expect(memberBlock.locator('tbody tr')).toHaveCount(2);
  await expect(memberBlock).toContainText('Visible protein');
  await expect(memberBlock).toContainText('Hidden override');
  const groupCopyButtons = memberBlock.locator(
    '.gfi-block-actions [data-copy-feedback-key]'
  );
  const groupNtCopy = groupCopyButtons.nth(0);
  const groupAaCopy = groupCopyButtons.nth(1);
  await expect(groupNtCopy).toHaveText('Copy nt (2)');
  await expect(groupAaCopy).toHaveText('Copy aa (2)');
  await groupAaCopy.click();
  await expect.poll(() => page.evaluate(() => window.__copiedText)).toContain('>VP_1');
  const copied = await page.evaluate(() => window.__copiedText);
  expect(copied).toContain('>HP_1');
  expect(copied).toContain('MHIDDEN');
  await expect(groupAaCopy).toHaveText('Copied!');
  await expect(groupAaCopy).toHaveAccessibleName('Copied!');
  await expect(groupAaCopy).toHaveAttribute('aria-live', 'polite');
  await expect(groupAaCopy).toHaveAttribute('aria-atomic', 'true');
  await expect(groupNtCopy).toHaveText('Copy nt (2)');

  await page.getByRole('button', { name: 'Qualifiers', exact: true }).click();
  await page.getByRole('button', { name: 'Details', exact: true }).click();
  await expect(groupAaCopy).toHaveText('Copied!');
  await page.clock.fastForward(1500);
  await expect(groupAaCopy).toHaveText('Copy aa (2)');

  const firstMemberCopyButtons = memberBlock.locator(
    'tbody tr'
  ).first().locator('[data-copy-feedback-key]');
  const secondMemberAaCopy = memberBlock.locator(
    'tbody tr'
  ).nth(1).locator('[data-copy-feedback-key]').nth(1);
  const firstMemberAaCopy = firstMemberCopyButtons.nth(1);
  await firstMemberAaCopy.click();
  await expect(firstMemberAaCopy).toHaveText('Copied!');
  await expect(secondMemberAaCopy).toHaveText('Copy aa');
  expect(await page.evaluate(() => window.__copiedText)).toContain('>VP_1');
  expect(await page.evaluate(() => window.__copiedText)).not.toContain('>HP_1');
  await page.clock.fastForward(1500);
  await expect(firstMemberAaCopy).toHaveText('Copy aa');

  await groupNtCopy.click();
  await expect(groupNtCopy).toHaveText('Copied!');
  await expect(groupAaCopy).toHaveText('Copy aa (2)');
  const copiedNucleotide = await page.evaluate(() => window.__copiedText);
  expect(copiedNucleotide).toContain('ATGAAATAA');
  expect(copiedNucleotide).toContain('ATGCCCTAA');
  await page.clock.fastForward(1500);

  await page.evaluate(() => {
    window.__manualCopy = null;
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: { writeText: async () => { throw new Error('denied'); } }
    });
    window.prompt = (_title, value) => {
      window.__manualCopy = String(value);
      return null;
    };
  });
  await groupAaCopy.click();
  await expect(groupAaCopy).toHaveText('Copy manually');
  await expect(groupAaCopy).not.toHaveText('Copied!');
  expect(await page.evaluate(() => window.__manualCopy)).toBe(copied);
  await page.clock.fastForward(1500);
  await expect(groupAaCopy).toHaveText('Copy aa (2)');

  await page.evaluate(() => {
    window.__manualCopy = null;
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {}
    });
  });
  await groupAaCopy.click();
  await expect(groupAaCopy).toHaveText('Copy manually');
  expect(await page.evaluate(() => window.__manualCopy)).toBe(copied);
  await page.clock.fastForward(1500);

  await page.evaluate(() => {
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: { writeText: async () => { throw new Error('denied'); } }
    });
    window.prompt = () => { throw new Error('prompt failed'); };
  });
  await groupAaCopy.click();
  await expect(groupAaCopy).toHaveText('Copy failed');
  await expect(groupAaCopy).not.toHaveText('Copied!');
  await page.clock.fastForward(1500);

  await page.evaluate(() => {
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {
        writeText: async (value) => { window.__copiedText = String(value); }
      }
    });
  });
  await groupAaCopy.click();
  await page.clock.fastForward(1000);
  await groupAaCopy.click();
  await page.clock.fastForward(600);
  await expect(groupAaCopy).toHaveText('Copied!');
  await page.clock.fastForward(900);
  await expect(groupAaCopy).toHaveText('Copy aa (2)');

  await page.evaluate(() => {
    window.__copyResolvers = [];
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {
        writeText: (value) => {
          window.__copiedText = String(value);
          return new Promise((resolve, reject) => {
            window.__copyResolvers.push({ reject, resolve });
          });
        }
      }
    });
  });
  await groupAaCopy.click();
  await groupAaCopy.click();
  await page.evaluate(() => window.__copyResolvers[1].resolve());
  await expect(groupAaCopy).toHaveText('Copied!');
  await page.evaluate(() => window.__copyResolvers[0].reject(new Error('stale denial')));
  await expect(groupAaCopy).toHaveText('Copied!');
  await page.clock.fastForward(1500);
  await expect(groupAaCopy).toHaveText('Copy aa (2)');

  await page.evaluate(() => {
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {
        writeText: async (value) => { window.__copiedText = String(value); }
      }
    });
  });
  await groupAaCopy.click();
  await expect(groupAaCopy).toHaveText('Copied!');
  await page.locator('[data-close]').click();
  await page.locator('[data-gbdraw-rendered-feature-id="rendered-visible"]').click();
  const reopenedMemberBlock = page.locator('.gfi-block').filter({
    hasText: 'Similarity-group members'
  }).last();
  await expect(reopenedMemberBlock.locator(
    '.gfi-block-actions [data-copy-feedback-key]'
  ).nth(1)).toHaveText('Copy aa (2)');

  await page.locator('[data-close]').click();
  await page.locator('[data-gbdraw-rendered-feature-id="rendered-collision"]').click();
  await page.getByRole('button', { name: 'Sequence' }).click();
  const nucleotideBlock = page.locator('.gfi-block').filter({
    hasText: 'Nucleotide'
  }).last();
  await nucleotideBlock.getByRole('button', { name: 'Copy', exact: true }).click();
  const nucleotideFasta = await page.evaluate(() => window.__copiedText);
  expect(nucleotideFasta).toContain('>record:21..29');
  expect(nucleotideFasta).toContain('ATGCCCTAA');
  expect(nucleotideFasta).not.toMatch(/h_[a-z2-7]{26}/i);

  await page.locator('[data-close]').click();
  const annotation = page.locator('#annotation-review-window');
  await expect(annotation).toHaveAttribute('data-gbdraw-interactive-annotation', 'true');
  await annotation.click();
  await expect(page.locator('.gfi-title')).toContainText('Review window');
  await expect(page.locator('#gbdraw-feature-popup')).toContainText('3..8');
  const labelRow = page.locator('.gfi-row').filter({ hasText: 'Label' });
  await labelRow.getByRole('button', { name: 'Copy' }).click();
  await expect.poll(() => page.evaluate(() => window.__copiedText))
    .toBe('Review window');

  await page.locator('[data-close]').click();
  await page.getByRole('button', { name: 'Expand feature search' }).click();
  await page.locator('[data-search-query]').fill('Visible protein');
  await page.locator('[data-search-apply]').click();
  await expect(page.locator('.gbdraw-interactive-feature--match')).toHaveCount(1);
});

test('standalone rejects conflicting compact provenance markers', async ({
  page
}, testInfo) => {
  await page.goto('/');
  const origin = new URL(page.url()).origin;
  const variants = await page.evaluate(async ({ origin, catalogSchema }) => {
    const { enrichSvgWithStandaloneInteractivity } = await import(
      `${origin}/gbdraw/web/js/services/standalone-interactivity.js`
    );
    const makeCatalog = () => ({
      schema: catalogSchema,
      items: [{
        resultIndex: 0,
        resultName: 'invalid.svg',
        recordKeys: ['record-key'],
        features: [{
          svgId: 'rendered-invalid',
          recordKey: 'record-key',
          biologicalFeatureId: 'invalid-feature'
        }],
        biologicalFeatures: [{
          recordKey: 'record-key',
          biologicalFeatureId: 'invalid-feature',
          record_id: 'record-id',
          type: 'CDS',
          start: 0,
          end: 3,
          strand: '+',
          sequenceSourceIndex: 0,
          amino_acid_sequence: 'M',
          translationFromAminoAcidSequence: true,
          qualifiers: {}
        }],
        orthogroups: [],
        annotations: [],
        comparisonMatches: [],
        sequenceSources: [{
          origin: 'linear-record',
          recordIndex: 0,
          sequence: 'ATG'
        }]
      }]
    });
    const cases = [];
    for (const conflict of [
      'DIFFERENT', '', ' ', null, [], [null], 0, false
    ]) {
      const catalog = makeCatalog();
      catalog.items[0].biologicalFeatures[0]
        .qualifiers.translation = conflict;
      cases.push({ name: `translation-${cases.length}`, catalog });
    }
    for (const invalidAminoAcid of [0, false, {}, []]) {
      const catalog = makeCatalog();
      catalog.items[0].biologicalFeatures[0]
        .amino_acid_sequence = invalidAminoAcid;
      cases.push({ name: `amino-${cases.length}`, catalog });
    }
    for (const shadowingValue of [null, '']) {
      const catalog = makeCatalog();
      const feature = catalog.items[0].biologicalFeatures[0];
      feature.aminoAcidSequence = feature.amino_acid_sequence;
      feature.amino_acid_sequence = shadowingValue;
      cases.push({ name: `amino-alias-${cases.length}`, catalog });
    }
    for (const invalidSequence of [123, {}, 'AT G']) {
      const catalog = makeCatalog();
      catalog.items[0].sequenceSources[0].sequence = invalidSequence;
      cases.push({ name: `sequence-${cases.length}`, catalog });
    }
    const unreferencedInvalidSource = makeCatalog();
    unreferencedInvalidSource.items[0].sequenceSources.push({
      origin: 'linear-record',
      recordIndex: 0,
      sequence: 'AT G'
    });
    cases.push({
      name: 'sequence-unreferenced',
      catalog: unreferencedInvalidSource
    });
    const coexisting = makeCatalog();
    coexisting.items[0].biologicalFeatures[0].nucleotide_sequence = 'ATG';
    cases.push({ name: 'coexisting-sequence', catalog: coexisting });

    return cases.map(({ name, catalog }) => {
      const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
      svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
      svg.setAttribute('viewBox', '0 0 100 50');
      svg.innerHTML = `
        <rect id="rendered-invalid"
          data-gbdraw-feature-id="invalid-feature"
          data-gbdraw-rendered-feature-id="rendered-invalid"
          x="5" y="5" width="20" height="10" />`;
      enrichSvgWithStandaloneInteractivity(svg, {
        popupMode: 'rich',
        featureCatalog: catalog,
        catalogResultIndex: 0,
        catalogResultName: 'invalid.svg',
        requireFeatureCatalog: true
      });
      return {
        name,
        svgText: new XMLSerializer().serializeToString(svg)
      };
    });
  }, { origin, catalogSchema: CURRENT_FEATURE_CATALOG_SCHEMA });

  await page.addInitScript(() => {
    window.__expandedInvalidFeature = false;
    const nativeMapSet = Map.prototype.set;
    Map.prototype.set = function (key, value) {
      if (
        value
        && (
          value.biologicalFeatureId === 'invalid-feature'
          || value.biological_feature_id === 'invalid-feature'
        )
      ) {
        window.__expandedInvalidFeature = true;
      }
      return nativeMapSet.call(this, key, value);
    };
  });

  for (const variant of variants) {
    const svgPath = testInfo.outputPath(`${variant.name}.svg`);
    fs.writeFileSync(svgPath, variant.svgText, 'utf8');
    await page.goto(pathToFileURL(svgPath).href);
    await expect(page.locator('#gbdraw-feature-search-controls')).toBeAttached();
    expect(await page.evaluate(() => window.__expandedInvalidFeature)).toBe(false);
  }
});

test('Download Interactive SVG forwards live editor overrides without mutating the catalog', async ({
  page
}) => {
  await page.goto('/');
  const origin = new URL(page.url()).origin;
  await page.addScriptTag({
    url: '/gbdraw/web/vendor/vue/vue.global.js'
  });
  await page.addScriptTag({
    url: '/gbdraw/web/vendor/dompurify/purify.min.js'
  });

  const exported = await page.evaluate(async ({ origin, catalogSchema }) => {
    const { state } = await import(`${origin}/gbdraw/web/js/state.js`);
    const { downloadInteractiveSVG } = await import(
      `${origin}/gbdraw/web/js/services/export.js`
    );
    const { captureSvgExport } = await import(
      `${origin}/gbdraw/web/js/services/svg-serialization.js`
    );
    const catalog = {
      schema: catalogSchema,
      items: [{
        resultIndex: 0,
        resultName: 'live.svg',
        recordKeys: ['record-a'],
        features: [{
          svgId: 'rendered-a',
          recordKey: 'record-a',
          biologicalFeatureId: 'biological-a'
        }],
        biologicalFeatures: [{
          recordKey: 'record-a',
          biologicalFeatureId: 'biological-a',
          type: 'CDS',
          start: 0,
          end: 9,
          product: 'Original feature label'
        }],
        orthogroups: [{
          id: 'group-a',
          name: 'Original group',
          description: 'Original description',
          members: [{
            recordKey: 'record-a',
            biologicalFeatureId: 'biological-a'
          }]
        }],
        annotations: [],
        comparisonMatches: []
      }]
    };
    const container = document.createElement('div');
    container.innerHTML = `
      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 50">
        <rect id="rendered-a" data-gbdraw-feature-id="biological-a"
          data-gbdraw-rendered-feature-id="rendered-a"
          x="5" y="5" width="20" height="10" />
      </svg>`;
    document.body.appendChild(container);

    state.results.value = [{ name: 'live.svg', content: container.innerHTML }];
    state.selectedResultIndex.value = 0;
    state.svgContainer.value = container;
    state.featureCatalog.value = catalog;
    state.activeDrawing().featureOverrides[JSON.stringify([state.generatedMode.value, 'record-a', 'biological-a'])] = {
      scope: state.generatedMode.value,
      recordKey: 'record-a',
      biologicalFeatureId: 'biological-a',
      featureVisibility: null,
      labelVisibility: null,
      labelText: 'Live feature label',
      labelSourceText: null
    };
    state.activeDrawing().orthogroupNameOverrides['group-a'] = 'Live group name';
    state.activeDrawing().orthogroupDescriptionOverrides['group-a'] = 'Live description';

    let downloadedBlob = null;
    const originalCreateObjectURL = URL.createObjectURL;
    const originalRevokeObjectURL = URL.revokeObjectURL;
    const originalClick = HTMLAnchorElement.prototype.click;
    URL.createObjectURL = (blob) => {
      downloadedBlob = blob;
      return 'blob:gbdraw-test';
    };
    URL.revokeObjectURL = () => {};
    HTMLAnchorElement.prototype.click = () => {};
    try {
      await downloadInteractiveSVG(captureSvgExport(state, { interactive: true }));
      const svgText = await downloadedBlob.text();
      const doc = new DOMParser().parseFromString(svgText, 'image/svg+xml');
      const metadata = doc.querySelector('#gbdraw-interactive-feature-metadata');
      return {
        embedded: JSON.parse(metadata.textContent),
        sourceFeature: catalog.items[0].features[0],
        sourceGroup: catalog.items[0].orthogroups[0]
      };
    } finally {
      URL.createObjectURL = originalCreateObjectURL;
      URL.revokeObjectURL = originalRevokeObjectURL;
      HTMLAnchorElement.prototype.click = originalClick;
    }
  }, { origin, catalogSchema: CURRENT_FEATURE_CATALOG_SCHEMA });

  expect(exported.embedded.items[0].features[0].displayLabel)
    .toBe('Live feature label');
  expect(exported.embedded.items[0].orthogroups[0].display_name)
    .toBe('Live group name');
  expect(exported.embedded.items[0].orthogroups[0].description)
    .toBe('Live description');
  expect(exported.sourceFeature.displayLabel).toBeUndefined();
  expect(exported.sourceGroup.display_name).toBeUndefined();
  expect(exported.sourceGroup.description).toBe('Original description');
});
