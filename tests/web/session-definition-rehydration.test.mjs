import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { readFile } from 'node:fs/promises';

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
globalThis.File = File;
globalThis.alert = () => {};

const { importSession } = await import('../../gbdraw/web/js/services/config.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const {
  stageCircularDefinitionSource
} = await import('../../gbdraw/web/js/app/results.js');

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
assert.ok(state.files.c_gb instanceof File);
assert.ok(state.extractedFeatures.value.length > 0);
assert.equal(state.featureEditorStatus.status, 'summary-ready');
assert.equal(state.featureEditorStatus.summaryCount, state.extractedFeatures.value.length);
assert.equal(state.featureExtractionPending.value, false);
assert.equal(state.featureExtractionError.value, null);

let stagedPath = '';
let stagedBytes = null;
const resolvedPath = await stageCircularDefinitionSource({
  inputType: state.cInputType.value,
  inputFile: state.files.c_gb,
  writeFileToFs: async (file, path) => {
    stagedPath = path;
    stagedBytes = new Uint8Array(await file.arrayBuffer());
    return true;
  }
});

assert.equal(resolvedPath, '/definition-input.gb');
assert.equal(stagedPath, resolvedPath);
assert.match(new TextDecoder().decode(stagedBytes), /VERSION\s+NC_012920\.1/);

let gffWriteAttempted = false;
assert.equal(await stageCircularDefinitionSource({
  inputType: 'gff',
  inputFile: state.files.c_gb,
  writeFileToFs: async () => {
    gffWriteAttempted = true;
    return true;
  }
}), null);
assert.equal(gffWriteAttempted, false);

const wssvSession = await readFile(
  'gbdraw/web/gallery/sessions/WSSV_genome_comparison.gbdraw-session.json',
  'utf8'
);
const importedWssv = await importSession({
  target: {
    files: [new Blob([wssvSession], { type: 'application/json' })],
    value: 'selected'
  }
});
assert.equal(importedWssv.status, 'ok');
const comparisonSources = state.matchSequenceRegistry.values()
  .filter((source) => source.origin === 'homology-comparison');
assert.equal(comparisonSources.length, 20);
assert.deepEqual(
  comparisonSources.map((source) => source.sourceIndex),
  Array.from({ length: 20 }, (_, index) => index)
);
assert.ok(comparisonSources.every((source) => source.sequence.length > 0));

console.log('session definition resource rehydration tests passed');
