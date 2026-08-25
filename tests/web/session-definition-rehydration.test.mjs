import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { readFile } from 'node:fs/promises';
import { gunzipSync } from 'node:zlib';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

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

const { importSession } = await import('../../gbdraw/web/js/services/config.js');
const { state } = await import('../../gbdraw/web/js/state.js');
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
assert.equal(state.files.c_gb instanceof File, false);
assert.equal(typeof state.files.c_gb.arrayBuffer, 'undefined');
assert.ok(Object.isFrozen(state.files.c_gb));
assert.ok(state.extractedFeatures.value.length > 0);
assert.equal(state.featureEditorStatus.status, 'summary-ready');
assert.equal(state.featureEditorStatus.summaryCount, state.extractedFeatures.value.length);
assert.equal(state.featureExtractionPending.value, false);
assert.equal(state.featureExtractionError.value, null);

const maliciousLegendSession = JSON.parse(session);
maliciousLegendSession.editorState.legend = {
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
const importedMaliciousLegend = await importSession({
  target: {
    files: [new Blob([JSON.stringify(maliciousLegendSession)], { type: 'application/json' })],
    value: 'selected'
  }
});
assert.equal(importedMaliciousLegend.status, 'ok');
assert.equal(
  state.legendEntries.value.some((entry) => entry.caption === 'unsafe-current'),
  false
);
assert.equal(
  state.deletedLegendEntries.value.some((entry) => entry.caption === 'unsafe-deleted'),
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
