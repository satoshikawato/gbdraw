import { pythonLegacySessionUpgrader } from './helpers/python-evaluator.mjs';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { readFile } from 'node:fs/promises';
import { gunzipSync } from 'node:zlib';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {};
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};

installSessionImportWorker();
// A Session 27-44 moves to a Session 49 in Python, as the composition root's
// Worker port does it (DESIGN 4.3).
const { setLegacySessionUpgrader } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
setLegacySessionUpgrader(pythonLegacySessionUpgrader());

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const { importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));
const session = await readFile(
  'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json',
  'utf8'
);
const imported = await importSession({
  target: {
    files: [new Blob([session], { type: 'application/json' })],
    value: 'selected'
  }
});

assert.equal(imported.status, 'ok');
const circularGenome = modeDrawing(state, 'circular').files.c_gb;
assert.equal(circularGenome instanceof File, false);
assert.equal(typeof circularGenome.arrayBuffer, 'undefined');
assert.ok(Object.isFrozen(circularGenome));
assert.ok(state.extractedFeatures.value.length > 0);
assert.deepEqual(state.featureEditorStatus.value, { status: 'summary-ready', error: null });
assert.equal(state.featureExtractionPending.value, false);
const loadedFeatureCount = state.extractedFeatures.value.length;
assert.equal(state.featureEditorStatusText.value, `${loadedFeatureCount} features ready`);

// OV-1762: the Editor's feature count is the shown drawing's (its `extractedFeatures`), not the count the
// last Load or Generate in any drawing wrote. A drawing switch swaps the artifact slots.
const { createHistorySnapshotService } = await import('../../gbdraw/web/js/services/history-snapshot.js');
const { createInputFileStore } = await import('../../gbdraw/web/js/services/history-files.js');
const slotSnapshots = createHistorySnapshotService({ state, fileStore: createInputFileStore(), currentResults: () => state.results.value });
// The Load made no Linear drawing (D-15-120); the offer tab makes it.
(await (await import('./helpers/drawing-state.mjs')).modeDrawingsAdder())();
const circularSlot = slotSnapshots.captureArtifactSlot();
slotSnapshots.installArtifactSlot(null, { drawing: modeDrawing(state, 'linear') });
assert.equal(state.featureEditorStatusText.value, 'No features ready', 'the Linear drawing has no Result');
slotSnapshots.installArtifactSlot(circularSlot, { drawing: modeDrawing(state, 'circular') });
assert.equal(state.featureEditorStatusText.value, `${loadedFeatureCount} features ready`);

// The Legend edits are the shown Circular drawing's (Session 49 `drawings[]`, by ID).
const maliciousLegendSession = JSON.parse(session);
const circularSlice = maliciousLegendSession.drawings.find((entry) => entry.id === maliciousLegendSession.ui.activeDrawingId);
circularSlice.editorState = { ...(circularSlice.editorState || {}) };
circularSlice.editorState.legend = {
  entries: [{
  caption: 'unsafe-current',
  color: 'url(javascript:alert(1))',
  showStroke: true,
  featureIds: ['forged-feature']
  }],
  deletedEntries: [{
    caption: 'unsafe-deleted',
    color: 'url(javascript:alert(2))'
  }]
};
// Load lists the entries it drops for their color before it applies the
// file (BA-D25); the port answers [Load].
const maliciousLegendDialogs = [];
const importedMaliciousLegend = await importSession({
  target: {
    files: [new Blob([JSON.stringify(maliciousLegendSession)], { type: 'application/json' })],
    value: 'selected'
  }
}, { confirmLoadRepairs: async (dialog) => { maliciousLegendDialogs.push(dialog); return true; } });
assert.equal(importedMaliciousLegend.status, 'ok');
assert.deepEqual(maliciousLegendDialogs.map(({ items }) => items), [[
  'Legend: entries "unsafe-current" and "unsafe-deleted" (deleted) had a color the app does not accept and were dropped.'
]]);
assert.equal(
  state.activeDrawing().legendEntries.value.some((entry) => entry.caption === 'unsafe-current'),
  false
);
assert.equal(
  state.activeDrawing().deletedLegendEntries.value.some((entry) => entry.caption === 'unsafe-deleted'),
  false
);

const conservationSession = gunzipSync(await readFile(
  'tests/fixtures/sessions/synthetic_conservation.gbdraw-session.json.gz'
));
const importedConservation = await importSession({
  target: {
    files: [new Blob([conservationSession], { type: 'application/json' })],
    value: 'selected'
  }
});
assert.equal(importedConservation.status, 'ok');
const comparisonSources = state.matchSequenceRegistry.values()
  .filter((source) => source.origin === 'homology-comparison');
assert.equal(comparisonSources.length, 3);
assert.deepEqual(
  comparisonSources.map((source) => source.sourceIndex),
  [0, 1, 2]
);
assert.ok(comparisonSources.every((source) => source.sequence.length > 0));

console.log('session definition resource rehydration tests passed');
