// A Linear card's windows in the app's own state (PD-OI-124). They survive Save and Load: a window across the origin,
// a reversed window, and an empty window (Add window before its Start and End are typed) load as the card saved
// them. Add window and Remove go through the record rows' one transition, which empties the drawing's stored
// comparisons, whose unit indexes the change shifts; a Start or End edit keeps them (F7).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { File } from 'node:buffer';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { gallerySessionDocument, readSavedSession, sessionLoadEvent } from './helpers/session-documents.mjs';
import { validateWebConfigOverrides, pythonLegacySessionUpgrader } from './helpers/python-evaluator.mjs';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }) };
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
installSessionImportWorker();
const { setLegacySessionUpgrader } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
setLegacySessionUpgrader(pythonLegacySessionUpgrader());
const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const { INPUTS_SLICE, bindDrawingCanonicalComparisons } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { createLinearRecordEdits, createLinearSourceEdits } = await import('../../gbdraw/web/js/app/linear-source-edits.js');
const { createLinearComparisonView } = await import('../../gbdraw/web/js/app/linear-comparison-panel.js');
const { exportSession, importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

const LAMBDA = await gallerySessionDocument('lambda_basic_linear');
const lambda = () => structuredClone(LAMBDA);
const load = async (file) => {
  const loaded = await importSession(await sessionLoadEvent(file), { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.ok(['ok', 'legacy'].includes(loaded.status), JSON.stringify(loaded.error));
};
const saved = async () => {
  const result = await exportSession('windows', { askChoice: async () => 'download' });
  assert.equal(result?.status, 'saved', JSON.stringify(result?.error ?? result));
  return readSavedSession(result.blob);
};

test('a card\'s windows load as the card saved them: across the origin, reversed, and empty', async () => {
  await load(lambda());
  const windows = [
    { start: 47001, end: 2000, reverseComplement: false, column: null, acrossOrigin: true },
    { start: 9001, end: 12000, reverseComplement: true, column: null, acrossOrigin: false },
    { start: null, end: null, reverseComplement: false, column: null, acrossOrigin: false }
  ];
  modeDrawing(state, 'linear').linearSeqs[0].windows = structuredClone(windows);
  const document = await saved();
  await load(document);
  assert.deepEqual(modeDrawing(state, 'linear').linearSeqs[0].windows, windows);
});

// The app's own transition of the record rows and its invalidation of the stored comparisons; the ports that only
// other edits need (Depth sources, discovery, the alignment plan) do nothing here.
const recordEditsOfApp = () => {
  const shownDrawingOfMode = () => modeDrawing(state, 'linear');
  const { invalidateLinearComparisonArtifacts } = createLinearComparisonView({
    state, mode: { value: 'linear' }, editLosatRows: () => {}, shownDrawingOfMode
  });
  const through = (_drawing, change) => change();
  const { applyLinearSeqMutation } = createLinearSourceEdits({
    sessionOperationAvailability: () => null, shownAction: () => undefined, retireLegendStylesOfUnnamedCaptions: through,
    changeLinearDepthSources: through, linearDepthSourceLoss: () => { throw new Error('no Depth source changes'); },
    changeRecordSources: through, askingStackPorts: () => { throw new Error('no question asked'); },
    pendingDiscoveryOf: () => ({ expansions: new Map(), uploads: new Set() }), replaceLinearComparisonPlan: () => {},
    invalidateLinearComparisonArtifacts, linearReorderNotice: { value: '' }, retainAlignmentForStableReorder: () => {},
    clearAlignmentForMutation: () => {}, rebindAnnotationRecords: () => {}, discoverDrawingRecords: async () => {}
  });
  return createLinearRecordEdits({
    sessionOperationAvailability: () => null, shownDrawingOfMode, runUndoable: (_label, action) => { action(); return true; },
    clearAlignmentForMutation: () => {}, linearRecordsFor: () => [], linearRecordStatusFor: () => ({}),
    rebindAnnotationRecords: () => {}, applyLinearSeqMutation, shownUnitPlaces: () => new Map()
  });
};
const storedComparisons = () => INPUTS_SLICE.capture(modeDrawing(state, 'linear').slices.inputs).files.linearCanonicalComparisons;
const storeComparison = () => {
  bindDrawingCanonicalComparisons(state.inputFiles, modeDrawing(state, 'linear').slices.inputs,
    [{ kind: 'collinearityResult', file: new File(['pair'], 'pair.tsv') }]);
  assert.equal(storedComparisons().length, 1);
};

test('Add window and Remove empty the stored comparisons; a Start edit keeps them (F7)', async () => {
  await load(lambda());
  const edits = recordEditsOfApp();
  const card = () => modeDrawing(state, 'linear').linearSeqs[0];
  storeComparison();
  edits.addLinearCardWindow(card());
  assert.equal(card().windows.length, 2);
  assert.deepEqual(storedComparisons(), []);
  storeComparison();
  edits.setLinearCardWindow(card(), 1, 'start', 9001);
  assert.equal(storedComparisons().length, 1);
  edits.removeLinearCardWindow(card(), 1);
  assert.equal(card().windows.length, 0);
  assert.deepEqual(storedComparisons(), []);
});
