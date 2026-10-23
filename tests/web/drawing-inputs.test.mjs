// The project's input file store (`state.inputFiles`, DRAWINGS-B design §2):
// one frozen `InputFileRef` per stored file, the bytes behind `fileOf`, and a
// retention that keeps every file a binding names even when no History step
// names it (R11: a file stays while any History or binding names it).
import assert from 'node:assert/strict';
import { test } from 'node:test';

globalThis.window = { Vue: {
  ref: (value) => ({ value }), reactive: (value) => value, toRaw: (value) => value,
  computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}
} };
const { state } = await import('../../gbdraw/web/js/state.js');
const {
  bindInputFile, bindInputFileList, bindLoadedInputs, boundInputFileIds, registerInputFiles, retainHistoryAndBoundFiles
} = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { createHistoryFileStore } = await import('../../gbdraw/web/js/services/history-files.js');

const fileNamed = (name) => new File([`>${name}\nACGT\n`], name, { type: 'text/plain', lastModified: 7 });

test('a bound file that no History step names survives the release and keeps its bytes', async () => {
  const store = state.inputFiles;
  const bound = fileNamed('colors.tsv');
  const boundSeq = fileNamed('record.gb');
  const released = fileNamed('old.tsv');
  const boundRef = store.registerFile(bound);
  const seqRef = store.registerFile(boundSeq);
  const releasedRef = store.registerFile(released);
  const previous = { d_color: state.files.d_color, gb: state.linearSeqs[0].gb };
  try {
    state.files.d_color = boundRef;
    state.linearSeqs[0].gb = seqRef;
    retainHistoryAndBoundFiles(state, store)([]);
    assert.equal(store.fileOf(boundRef), bound, 'the bound card file stays stored');
    assert.equal(store.fileOf(seqRef), boundSeq, 'the bound record file stays stored');
    assert.equal(await store.fileOf(boundRef).text(), '>colors.tsv\nACGT\n');
    assert.equal(store.fileOf(releasedRef), null, 'a file nothing names is released');
    assert.equal(store.has(releasedRef.fileId), false);
  } finally {
    state.files.d_color = previous.d_color;
    state.linearSeqs[0].gb = previous.gb;
  }
});

test('registering one File twice gives one frozen ref, and a ref is not a Blob', () => {
  const store = state.inputFiles;
  const file = fileNamed('genome.fasta');
  const ref = store.registerFile(file);
  assert.equal(store.registerFile(file), ref, 'one ref per stored file');
  assert.ok(Object.isFrozen(ref));
  assert.deepEqual({ ...ref }, {
    fileId: ref.fileId, sourceId: ref.fileId, name: 'genome.fasta', type: 'text/plain', size: file.size, lastModified: 7
  });
  assert.equal(ref instanceof Blob, false);
  assert.equal(typeof ref.arrayBuffer, 'undefined', 'a byte reader must resolve the ref with fileOf');
  assert.equal(store.registerFile(ref), ref, 'a ref of the store stays bound as it is');
  assert.equal(store.registerFile({ ...ref }), ref, 'a ref is never stored as a file; its fileId names the stored one');
  assert.equal(store.registerFile({ ...ref, fileId: 'file-none' }), null, 'a ref the store does not hold names nothing');
  assert.equal(store.fileOf(ref), file);
});

test('the binding writers store each file once and hold refs, dropping places without a file', () => {
  const store = state.inputFiles;
  const files = { d_color: null, c_conservation_blasts: [], c_conservation_fastas: [] };
  const colors = fileNamed('colors.tsv');
  const ref = bindInputFile(store, files, 'd_color', colors);
  assert.equal(files.d_color, ref);
  assert.equal(store.fileOf(files.d_color), colors);
  assert.equal(bindInputFile(store, files, 'd_color', null), null);
  assert.equal(files.d_color, null);
  const first = fileNamed('a.tsv');
  const second = fileNamed('b.tsv');
  const refs = bindInputFileList(store, files, 'c_conservation_blasts', [first, null, second]);
  assert.deepEqual(refs.map((entry) => entry.name), ['a.tsv', 'b.tsv']);
  assert.equal(files.c_conservation_blasts, refs);
  assert.deepEqual(registerInputFiles(store, [first, null]).map((entry) => entry && entry.name), ['a.tsv', null]);
  assert.equal(bindInputFileList(store, files, 'c_conservation_fastas', first)[0], refs[0], 'one ref per File');
});

test('the refs of two views of one stored resource share its source; other files keep their own', () => {
  const descriptor = {};
  const store = createHistoryFileStore({ sourceOf: (file) => file.backing || null });
  const first = store.registerFile({ name: 'two-records.gb', size: 9, backing: descriptor });
  const second = store.registerFile({ name: 'two-records.gb', size: 9, backing: descriptor });
  const other = store.registerFile(fileNamed('other.gb'));
  assert.notEqual(first.fileId, second.fileId, 'each view is stored once, as itself');
  assert.equal(second.sourceId, first.sourceId, 'one File card groups the views of one resource');
  assert.equal(other.sourceId, other.fileId);
});

test('a Load installs its Files as refs at install, and a rolled-back ref stays bound', () => {
  const store = createHistoryFileStore();
  const kept = store.registerFile(fileNamed('kept.gb'));
  const genome = fileNamed('genome.gb');
  const depth = fileNamed('depth.tsv');
  const edgeFile = fileNamed('pair.tsv');
  const view = fileNamed('comparison.tsv');
  const bound = bindLoadedInputs(store, {
    files: {
      c_gb: genome, c_depth: [[depth, null]], c_conservation_blasts: [kept], c_conservation_blasts_source: 'losat-cache',
      linearCanonicalComparisons: [{ kind: 'collinearityResult', file: view }, { kind: 'generatedProteinComparison' }]
    },
    linearSeqs: [{ uid: 'a', gb: genome, gff: null, fasta: null, depth: [depth] }],
    linearComparisonPlan: { mode: 'selected', defaultSource: 'upload', edges: [{ id: 'e', file: edgeFile }] },
    circularRecordDiscovery: { primaryFile: genome, pairedFile: null }
  });
  const genomeRef = store.registerFile(genome);
  assert.equal(bound.files.c_gb, genomeRef, 'the loaded File is registered and the field holds its ref');
  assert.deepEqual(bound.files.c_depth, [[store.registerFile(depth), null]], 'a Depth matrix keeps its places');
  assert.deepEqual(bound.files.c_conservation_blasts, [kept], 'a ref of the rollback stays bound as it is');
  assert.equal(bound.files.c_conservation_blasts_source, 'losat-cache', 'a flag is not a file');
  assert.equal(store.fileOf(bound.files.linearCanonicalComparisons[0].file), view);
  assert.equal(Object.hasOwn(bound.files.linearCanonicalComparisons[1], 'file'), false, 'a generated entry has no file');
  assert.equal(bound.linearSeqs[0].gb, genomeRef, 'one ref per File across fields');
  assert.deepEqual(bound.linearSeqs[0].depth, [store.registerFile(depth)]);
  assert.equal(store.fileOf(bound.linearComparisonPlan.edges[0].file), edgeFile);
  assert.equal(bound.circularRecordDiscovery.primaryFile, genomeRef, 'the discovery names the field\'s ref');
  const ids = boundInputFileIds({
    files: bound.files, linearSeqs: bound.linearSeqs,
    drawingList: () => [{ linearComparisonPlan: bound.linearComparisonPlan }]
  });
  assert.equal(ids.size, 5, 'every bound ref is named (genome, depth, kept, comparison, edge)');
  assert.equal(boundInputFileIds({ files: { c_gb: fileNamed('raw.gb') }, linearSeqs: [], drawingList: () => [] }).size, 0,
    'a binding holds refs only: a File in a field names nothing');
});
