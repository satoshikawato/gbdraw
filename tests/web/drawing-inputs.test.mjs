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
const { createInputFileStore } = await import('../../gbdraw/web/js/services/history-files.js');

const fileNamed = (name) => new File([`>${name}\nACGT\n`], name, { type: 'text/plain', lastModified: 7 });

test('a bound file that no History step names survives the release and keeps its bytes', async () => {
  const store = state.inputFiles;
  const bound = fileNamed('colors.tsv');
  const boundSeq = fileNamed('record.gb');
  const released = fileNamed('old.tsv');
  const boundRef = store.registerFile(bound);
  const seqRef = store.registerFile(boundSeq);
  const releasedRef = store.registerFile(released);
  const rows = state.drawingOfMode('linear').linearSeqs;
  const previous = { d_color: state.files.d_color, gb: rows[0].gb };
  try {
    state.files.d_color = boundRef;
    rows[0].gb = seqRef;
    retainHistoryAndBoundFiles(state, store)([]);
    assert.equal(store.fileOf(boundRef), bound, 'the bound card file stays stored');
    assert.equal(store.fileOf(seqRef), boundSeq, 'the bound record file stays stored');
    assert.equal(await store.fileOf(boundRef).text(), '>colors.tsv\nACGT\n');
    assert.equal(store.fileOf(releasedRef), null, 'a file nothing names is released');
    assert.equal(store.has(releasedRef.fileId), false);
  } finally {
    state.files.d_color = previous.d_color;
    rows[0].gb = previous.gb;
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
  const store = createInputFileStore({ sourceOf: (file) => file.backing || null });
  const first = store.registerFile({ name: 'two-records.gb', size: 9, backing: descriptor });
  const second = store.registerFile({ name: 'two-records.gb', size: 9, backing: descriptor });
  const other = store.registerFile(fileNamed('other.gb'));
  assert.notEqual(first.fileId, second.fileId, 'each view is stored once, as itself');
  assert.equal(second.sourceId, first.sourceId, 'one File card groups the views of one resource');
  assert.equal(other.sourceId, other.fileId);
});

test('a Load installs its Files as refs at install, and a rolled-back ref stays bound', () => {
  const store = createInputFileStore();
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
    files: bound.files,
    drawingList: () => [{ linearSeqs: bound.linearSeqs, linearComparisonPlan: bound.linearComparisonPlan }]
  });
  assert.equal(ids.size, 5, 'every bound ref is named (genome, depth, kept, comparison, edge)');
  assert.equal(boundInputFileIds({ files: { c_gb: fileNamed('raw.gb') }, drawingList: () => [] }).size, 0,
    'a binding holds refs only: a File in a field names nothing');
});

test('a file bound only by a hidden drawing survives History eviction (the union of every drawing\'s bindings)', () => {
  const store = createInputFileStore();
  const genome = fileNamed('hidden.gb');
  const record = fileNamed('hidden-record.gb');
  const genomeRef = store.registerFile(genome);
  const recordRef = store.registerFile(record);
  const shown = { linearComparisonPlan: null, files: { c_gb: null }, linearSeqs: [] };
  const hidden = { linearComparisonPlan: null, files: { c_gb: genomeRef }, linearSeqs: [{ uid: 'r', gb: recordRef }] };
  retainHistoryAndBoundFiles({ files: {}, linearSeqs: [], drawingList: () => [shown, hidden] }, store)([]);
  assert.equal(store.fileOf(genomeRef), genome, 'the hidden drawing\'s input file stays stored');
  assert.equal(store.fileOf(recordRef), record, 'the hidden drawing\'s record file stays stored');
});

const {
  INPUTS_SLICE, addLinearSeq, bindDrawingCanonicalComparisons, bindDrawingInputFile, bindDrawingInputFileList,
  bindDrawingInputFileTree, bindLinearSeqFile, createLinearSeq, moveLinearSeq, removeLinearSeq, replaceLinearSeqs,
  setConservationBlastsSource
} = await import('../../gbdraw/web/js/services/drawing/inputs.js');

// How each member is filled with a value: bound files, rows, a GFF input, and
// a read record catalog; the files it binds, and how a capture's nested list
// is changed to show the capture is a copy. A slice holds only its joined
// members (`INPUTS_SLICE.members`), which the tests below fill.
const MEMBER_FILLS = {
  files: { boundFiles: 4, fill: ({ store, handle, genome }) => {
    bindDrawingInputFile(store, handle, 'c_gb', genome);
    bindDrawingInputFileList(store, handle, 'c_conservation_blasts', [fileNamed('a.tsv'), null]);
    bindDrawingInputFileTree(store, handle, 'c_depth', [[fileNamed('depth.tsv'), null]]);
    setConservationBlastsSource(handle, 'losat-cache');
    bindDrawingCanonicalComparisons(store, handle, [{ kind: 'collinearityResult', file: fileNamed('pair.tsv') }]);
  }, changeCopy: (captured) => captured.files.c_depth[0].push('changed') },
  linearSeqs: { boundFiles: 2, fill: ({ store, handle }) => {
    addLinearSeq(handle, createLinearSeq({ uid: 'second', definition: 'Second' }));
    bindLinearSeqFile(store, handle, 'second', 'gb', fileNamed('second.gb'));
    bindLinearSeqFile(store, handle, 'second', 'depth', [fileNamed('second-depth.tsv')]);
  }, changeCopy: (captured) => captured.linearSeqs[1].depth.push('changed') },
  inputType: { boundFiles: 0, fill: ({ view }) => { view.inputType.value = 'gff'; }, changeCopy: () => {} },
  circularRecordList: { boundFiles: 0, fill: ({ view }) => {
    view.circularRecordList.value = [{ selector: '#1', record_id: 'NC_1' }];
  }, changeCopy: () => {} },
  circularRecordDiscovery: { boundFiles: 0, fill: ({ store, view, genome }) => {
    Object.assign(view.circularRecordDiscovery, {
      status: 'ready', inputType: 'gb', primaryFile: store.registerFile(genome),
      canonicalRecordIdentities: [{ selector: '#1', record_id: 'NC_1', recordKey: 'k1' }]
    });
  }, changeCopy: (captured) => { captured.circularRecordDiscovery.canonicalRecordIdentities[0].recordKey = 'changed'; } }
};
const filledInputs = () => {
  const store = createInputFileStore();
  const { handle, view } = INPUTS_SLICE.create();
  const genome = fileNamed('genome.gb');
  INPUTS_SLICE.members.forEach((member) => MEMBER_FILLS[member].fill({ store, handle, view, genome }));
  return { store, handle, view, genome };
};
const joinedFiles = INPUTS_SLICE.members.includes('files');

test('inputs: install(capture(x)) gives back x, also after reset, and a capture is a copy', () => {
  const { handle, view, store, genome } = filledInputs();
  if (joinedFiles) {
    assert.equal(store.fileOf(view.files.c_gb), genome, 'a file field holds the stored file\'s ref');
    assert.equal(Object.hasOwn(view.files, 'd_color'), false, 'the card files stay project fields (D-15-91 S-2 A)');
  }
  assert.deepEqual(Object.keys(view), [...INPUTS_SLICE.members]);
  const captured = INPUTS_SLICE.capture(handle);
  const json = JSON.stringify(captured);
  INPUTS_SLICE.install(handle, captured);
  assert.equal(JSON.stringify(INPUTS_SLICE.capture(handle)), json, 'install(capture(x))');
  INPUTS_SLICE.reset(handle);
  assert.equal(JSON.stringify(INPUTS_SLICE.capture(handle)), json, 'a reset keeps every input binding (PD-OI-070)');
  const empty = INPUTS_SLICE.create();
  INPUTS_SLICE.install(empty.handle, JSON.parse(json));
  INPUTS_SLICE.install(empty.handle, INPUTS_SLICE.capture(empty.handle));
  assert.equal(JSON.stringify(INPUTS_SLICE.capture(empty.handle)), json, 'another drawing installs the same inputs');
  INPUTS_SLICE.members.forEach((member) => MEMBER_FILLS[member].changeCopy(captured));
  assert.equal(JSON.stringify(INPUTS_SLICE.capture(handle)), json, 'a capture shares no list with the slice');
  if (joinedFiles) assert.equal(INPUTS_SLICE.capture(handle).files.c_gb, view.files.c_gb, 'refs are kept by reference');
});

test('inputs: the row writers find a row by its uid and refuse an unknown one', () => {
  const { handle, view } = filledInputs();
  const first = view.linearSeqs[0].uid;
  moveLinearSeq(handle, 'second', 0);
  assert.deepEqual(view.linearSeqs.map((seq) => seq.uid), ['second', first]);
  assert.throws(() => bindLinearSeqFile(createInputFileStore(), handle, 'none', 'gb', null), /Unknown Linear record row none/);
  assert.throws(() => addLinearSeq(handle, createLinearSeq({ uid: 'second' })), /exists/);
  removeLinearSeq(handle, 'second');
  removeLinearSeq(handle, first);
  assert.equal(view.linearSeqs.length, 1, 'a drawing always holds a row');
  replaceLinearSeqs(handle, [createLinearSeq({ uid: 'only' })]);
  assert.deepEqual(view.linearSeqs.map((seq) => seq.uid), ['only']);
  const other = INPUTS_SLICE.create();
  assert.throws(() => bindLinearSeqFile(createInputFileStore(), /** another slice's handle */ { __slice: 'canvas' }, 'only', 'gb', null));
  assert.deepEqual(other.view.linearSeqs.length, 1, 'a new slice holds one empty row');
});

test('the retention walk names every bound file of a slice\'s drawing', () => {
  const { handle, view, store } = filledInputs();
  const ids = boundInputFileIds({ files: {}, drawingList: () => [{ ...view, linearComparisonPlan: null }] });
  assert.equal(ids.size, INPUTS_SLICE.members.reduce((count, member) => count + MEMBER_FILLS[member].boundFiles, 0),
    'genome, BLAST, Depth, comparison, row GenBank, row Depth, of the joined members');
  retainHistoryAndBoundFiles({ files: {}, drawingList: () => [view] }, store)([]);
  assert.ok(INPUTS_SLICE.capture(handle).linearSeqs.every((seq) => !seq.gb || store.fileOf(seq.gb)));
});
