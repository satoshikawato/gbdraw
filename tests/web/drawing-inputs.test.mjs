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
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const {
  bindInputFile, bindLoadedInputs, boundInputFileIds, registerInputFiles, retainHistoryAndBoundFiles
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
  const rows = state.drawingOf('linear').linearSeqs;
  const hidden = state.drawingOf('linear');
  const previous = { d_color: hidden.files.d_color, gb: rows[0].gb };
  try {
    bindInputFile(store, hidden.slices.inputs, 'd_color', boundRef);
    rows[0].gb = seqRef;
    retainHistoryAndBoundFiles(state, store)([]);
    assert.equal(store.fileOf(boundRef), bound, 'the bound card file stays stored');
    assert.equal(store.fileOf(seqRef), boundSeq, 'the bound record file stays stored');
    assert.equal(await store.fileOf(boundRef).text(), '>colors.tsv\nACGT\n');
    assert.equal(store.fileOf(releasedRef), null, 'a file nothing names is released');
    assert.equal(store.has(releasedRef.fileId), false);
  } finally {
    bindInputFile(store, hidden.slices.inputs, 'd_color', previous.d_color);
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
  const { handle, view } = INPUTS_SLICE.create();
  const colors = fileNamed('colors.tsv');
  const ref = bindInputFile(store, handle, 'd_color', colors);
  assert.equal(view.files.d_color, ref);
  assert.equal(store.fileOf(view.files.d_color), colors);
  assert.equal(bindInputFile(store, handle, 'd_color', null), null);
  assert.equal(view.files.d_color, null);
  const first = fileNamed('a.tsv');
  const second = fileNamed('b.tsv');
  const refs = bindDrawingInputFileList(store, handle, 'c_conservation_blasts', [first, null, second]);
  assert.deepEqual(refs.map((entry) => entry.name), ['a.tsv', 'b.tsv']);
  assert.deepEqual(view.files.c_conservation_blasts, refs);
  assert.deepEqual(registerInputFiles(store, [first, null]).map((entry) => entry && entry.name), ['a.tsv', null]);
  assert.equal(bindDrawingInputFileList(store, handle, 'c_conservation_fastas', first)[0], refs[0], 'one ref per File');
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
    drawingList: () => [{ files: bound.files, linearSeqs: bound.linearSeqs, linearComparisonPlan: bound.linearComparisonPlan }]
  });
  assert.equal(ids.size, 5, 'every bound ref is named (genome, depth, kept, comparison, edge)');
  assert.equal(boundInputFileIds({ drawingList: () => [{ files: { c_gb: fileNamed('raw.gb'), d_color: fileNamed('raw.tsv') }, linearSeqs: [] }] }).size, 0,
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
  retainHistoryAndBoundFiles({ drawingList: () => [shown, hidden] }, store)([]);
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
    assert.equal(view.files.d_color, null, 'each drawing binds its own card files (R16-4)');
  }
  assert.deepEqual(Object.keys(view), [...INPUTS_SLICE.members]);
  const captured = INPUTS_SLICE.capture(handle);
  const json = JSON.stringify(captured);
  INPUTS_SLICE.install(handle, captured);
  assert.equal(JSON.stringify(INPUTS_SLICE.capture(handle)), json, 'install(capture(x))');
  INPUTS_SLICE.reset(handle);
  assert.equal(JSON.stringify(INPUTS_SLICE.capture(handle)), json, 'a reset keeps every input binding (PD-OI-070)');
  INPUTS_SLICE.install(handle, captured);
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
  assert.throws(() => bindLinearSeqFile(createInputFileStore(), handle, 'none', 'gb', null), /UNKNOWN_ROW/);
  assert.throws(() => addLinearSeq(handle, createLinearSeq({ uid: 'second' })), /DUPLICATE_ROW/);
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
  const ids = boundInputFileIds({ drawingList: () => [{ ...view, linearComparisonPlan: null }] });
  assert.equal(ids.size, INPUTS_SLICE.members.reduce((count, member) => count + MEMBER_FILLS[member].boundFiles, 0),
    'genome, BLAST, Depth, comparison, row GenBank, row Depth, of the joined members');
  retainHistoryAndBoundFiles({ drawingList: () => [view] }, store)([]);
  assert.ok(INPUTS_SLICE.capture(handle).linearSeqs.every((seq) => !seq.gb || store.fileOf(seq.gb)));
});

// R16-4 (DESIGN §2, §12): the five card files are each drawing's own. A card
// binds in the drawing it is chosen in; a reset of a drawing removes its card
// files and keeps its other bindings and every other drawing's cards.
test('a card file binds in its drawing only, and a reset removes only that drawing\'s card files', () => {
  const store = createInputFileStore();
  const a = INPUTS_SLICE.create();
  const b = INPUTS_SLICE.create();
  const colors = bindInputFile(store, a.handle, 'd_color', fileNamed('colors.tsv'));
  const rules = bindInputFile(store, b.handle, 't_color', fileNamed('rules.tsv'));
  const genome = bindDrawingInputFile(store, a.handle, 'c_gb', fileNamed('genome.gb'));
  assert.deepEqual([a.view.files.d_color, a.view.files.t_color, b.view.files.d_color, b.view.files.t_color],
    [colors, null, null, rules]);
  assert.deepEqual([...boundInputFileIds({ drawingList: () => [a.view] })].sort(), [colors.fileId, genome.fileId].sort(),
    'a drawing\'s card file is one of its bindings (retention keeps it while any drawing binds it)');
  INPUTS_SLICE.reset(a.handle);
  assert.deepEqual([a.view.files.d_color, a.view.files.c_gb, b.view.files.t_color], [null, genome, rules],
    'the reset drawing loses its card files and keeps its genome; the other drawing keeps its card');
  assert.equal(bindInputFile(store, b.handle, 't_color', null), null, 'Remove unbinds the card');
  assert.equal(b.view.files.t_color, null);
});

// PD-OI-090, R16-4: a card change is a step of the shown drawing. Its Undo,
// from that drawing, gives the card its earlier file back; another drawing's
// card and its steps stay; an Undo while the other drawing is shown undoes
// that drawing's own step, never the first drawing's card step.
test('a card step is the shown drawing\'s: Undo from another drawing never undoes it', async () => {
  const { buildConfigData, applyCapturedConfigData, installShownDrawing } = await import('../../gbdraw/web/js/services/config.js');
  const { createHistoryRouter } = await import('../../gbdraw/web/js/services/history.js');
  const { createHistorySnapshotService } = await import('../../gbdraw/web/js/services/history-snapshot.js');
  const snapshots = createHistorySnapshotService({
    state, fileStore: state.inputFiles, currentResults: () => state.results.value, buildConfigData, applyConfigData: applyCapturedConfigData
  });
  const history = createHistoryRouter({
    drawingIds: ['circular', 'linear'],
    activeDrawingId: () => state.activeDrawingId.value,
    managerOptions: (drawingId) => snapshots.drawingHistory(drawingId),
    sessionHistory: snapshots.sessionHistory()
  });
  const circular = state.drawingOf('circular');
  const linear = state.drawingOf('linear');
  installShownDrawing(circular.id);
  await history.initializeIntentBaseline();
  const rules = fileNamed('rules.tsv');
  const colors = fileNamed('colors.tsv');
  const ref = (file) => state.inputFiles.registerFile(file);
  await history.runUndoable('Change uploaded file', () => { bindInputFile(state.inputFiles, circular.slices.inputs, 't_color', rules); });
  installShownDrawing(linear.id);
  assert.equal(linear.files.t_color, null, 'the other drawing\'s card shows its own file (none)');
  await history.runUndoable('Change uploaded file', () => { bindInputFile(state.inputFiles, linear.slices.inputs, 'd_color', colors); });
  await history.undo();
  assert.deepEqual([linear.files.d_color, circular.files.t_color], [null, ref(rules)],
    'Undo in the second drawing takes back its card only');
  await history.undo();
  assert.equal(circular.files.t_color, ref(rules), 'an Undo while the second drawing is shown never undoes the first drawing\'s card step');
  installShownDrawing(circular.id);
  await history.undo();
  assert.deepEqual([circular.files.t_color, linear.files.d_color], [null, null], 'Undo from the first drawing takes back its card');
  await history.redo();
  assert.deepEqual([circular.files.t_color, linear.files.d_color], [ref(rules), null]);
  bindInputFile(state.inputFiles, circular.slices.inputs, 't_color', null);
});

// S-1 (DESIGN §2): a drawing's History intent holds its inputs slice capture,
// the slice's own copy, not a field list of History's. It leaves out what is
// not an edit: the comparisons a Generate writes (the artifact checkpoint keeps
// them) and the Circular record list and discovery, which are derived status.
// The Linear comparison plan's files stay with the plan (settings slice).
test('a History intent\'s files are the drawing\'s inputs slice capture', async () => {
  const { buildConfigData, applyCapturedConfigData } = await import('../../gbdraw/web/js/services/config.js');
  const { createHistorySnapshotService } = await import('../../gbdraw/web/js/services/history-snapshot.js');
  const { INPUTS_SLICE } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  const snapshots = createHistorySnapshotService({
    state, fileStore: state.inputFiles, currentResults: () => state.results.value, buildConfigData, applyConfigData: applyCapturedConfigData
  });
  const circular = state.drawingOf('circular');
  bindInputFile(state.inputFiles, circular.slices.inputs, 'whitelist', fileNamed('keep.tsv'));
  for (const drawing of [circular, state.drawingOf('linear')]) {
    const { files, linearSeqs, inputType } = INPUTS_SLICE.capture(drawing.slices.inputs);
    const { linearCanonicalComparisons: _generated, ...bound } = files;
    const intent = snapshots.drawingHistory(drawing.id).buildIntent();
    const expected = JSON.parse(JSON.stringify({
      ...bound, inputType, linearSeqs,
      linearComparisons: drawing.linearComparisonPlan.edges.map((edge) => ({ id: edge.id, file: edge.file ?? null }))
    }));
    assert.deepEqual(intent.files, expected, drawing.mode);
  }
  bindInputFile(state.inputFiles, circular.slices.inputs, 'whitelist', null);
});

// S-1 apply side: Undo installs the step's inputs slice capture through the
// slice, so every field the slice holds comes back as captured (no per-mode
// field list; D-15-127 J5), a ref whose file was released binds nothing, and
// the comparisons a Generate wrote stay (an intent leaves them out).
test('Undo of an input-file step installs the step\'s inputs slice capture', async () => {
  const { buildConfigData, applyCapturedConfigData, installShownDrawing } = await import('../../gbdraw/web/js/services/config.js');
  const { createHistoryRouter } = await import('../../gbdraw/web/js/services/history.js');
  const { createHistorySnapshotService } = await import('../../gbdraw/web/js/services/history-snapshot.js');
  const { INPUTS_SLICE, bindDrawingInputFile, bindDrawingCanonicalComparisons, bindLinearSeqFile } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  const store = state.inputFiles;
  const snapshots = createHistorySnapshotService({
    state, fileStore: store, currentResults: () => state.results.value, buildConfigData, applyConfigData: applyCapturedConfigData
  });
  const history = createHistoryRouter({
    drawingIds: ['circular', 'linear'],
    activeDrawingId: () => state.activeDrawingId.value,
    managerOptions: (drawingId) => snapshots.drawingHistory(drawingId),
    sessionHistory: snapshots.sessionHistory()
  });
  const linear = state.drawingOf('linear');
  const slice = linear.slices.inputs;
  const initial = INPUTS_SLICE.capture(slice);
  installShownDrawing(linear.id);
  await history.initializeIntentBaseline();
  const uid = linear.linearSeqs[0].uid;
  await history.runUndoable('Change uploaded file', () => {
    bindLinearSeqFile(store, slice, uid, 'gb', fileNamed('first.gb'));
    bindDrawingInputFile(store, slice, 'c_fasta', fileNamed('kept.fa'));
  });
  const stepCapture = INPUTS_SLICE.capture(slice);
  await history.runUndoable('Change uploaded file', () => {
    bindLinearSeqFile(store, slice, uid, 'gb', fileNamed('second.gb'));
    bindDrawingInputFile(store, slice, 'c_fasta', null);
  });
  const generated = [{ kind: 'generatedProteinComparison', queryRecordIndex: 0, subjectRecordIndex: 1 }];
  bindDrawingCanonicalComparisons(store, slice, generated);
  const released = stepCapture.linearSeqs[0].gb;
  // Only the first row file is released; the step's other file stays stored.
  const stored = [...boundInputFileIds({ drawingList: () => state.drawingList() }), stepCapture.files.c_fasta.fileId];
  store.retainOnly(stored.filter((id) => id !== released.fileId));
  await history.undo();
  const expected = {
    ...stepCapture,
    files: { ...stepCapture.files, linearCanonicalComparisons: generated },
    linearSeqs: stepCapture.linearSeqs.map((seq) => ({ ...seq, gb: null }))
  };
  const { circularRecordList: _list, circularRecordDiscovery: _discovery, ...restored } = INPUTS_SLICE.capture(slice);
  const { circularRecordList: _l, circularRecordDiscovery: _d, ...wanted } = expected;
  assert.deepEqual(restored, wanted, 'the slice holds the step\'s capture: the released row file binds nothing, the generated comparisons stay');
  INPUTS_SLICE.install(slice, initial);
});

// Save writes each drawing's own five card bindings (DESIGN §1
// `drawings[i].webFiles.bindings`).
test('Save writes each drawing\'s own card bindings', async () => {
  const { buildSessionDrawingResources } = await import('../../gbdraw/web/js/services/session-resources.js');
  const store = createInputFileStore();
  const colors = store.registerFile(new File(['CDS\t#123456\n'], 'colors.tsv', { type: 'text/plain', lastModified: 1 }));
  const words = store.registerFile(new File(['hypothetical\n'], 'black.txt', { type: 'text/plain', lastModified: 2 }));
  const drawings = [
    { mode: 'circular', files: { d_color: colors }, linearSeqs: [] },
    { mode: 'linear', files: { blacklist: words }, linearSeqs: [], linearComparisonPlan: null }
  ];
  const { written } = await buildSessionDrawingResources({ inputFiles: store }, drawings.map((drawing) => ({ drawing, committed: null })));
  const cards = written.map(({ webFiles }) => ['d_color', 'blacklist'].map((key) => Boolean(webFiles.bindings[key])));
  assert.deepEqual(cards, [[true, false], [false, true]], 'each drawing names its own card files only');
});

// PD-OI-090 〔入力ファイル〕, T-B4-6: another drawing chooses a file the project
// already holds without an upload. The picker offers the files every drawing
// binds in a slot of the uploader's kind; choosing one binds the same ref (the
// same file ID) through the uploader's own bind path and stores nothing new.
test('a Session file chosen in another drawing binds the same file ID and stores no new file', async () => {
  const { sessionFilesOfKind } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  const store = createInputFileStore();
  const circular = INPUTS_SLICE.create();
  const linear = INPUTS_SLICE.create();
  const genome = bindDrawingInputFile(store, circular.handle, 'c_gb', fileNamed('genome.gb'));
  const colors = bindInputFile(store, circular.handle, 'd_color', fileNamed('colors.tsv'));
  const blast = store.registerFile(fileNamed('pair.tsv'));
  const bindings = { drawingList: () => [
    { ...circular.view, linearComparisonPlan: null },
    { ...linear.view, linearComparisonPlan: { edges: [{ file: blast }] } }
  ] };
  assert.deepEqual(sessionFilesOfKind(bindings, 'genbank'), [genome], 'a GenBank uploader offers the bound GenBank');
  assert.deepEqual(sessionFilesOfKind(bindings, 'd_color'), [colors], 'a card offers its own kind only');
  assert.deepEqual(sessionFilesOfKind(bindings, 'comparison'), [blast]);
  assert.deepEqual(sessionFilesOfKind(bindings, 'fasta'), [], 'no FASTA is bound');
  assert.deepEqual(sessionFilesOfKind(bindings, 'blacklist'), []);
  const [picked] = sessionFilesOfKind(bindings, 'genbank');
  bindLinearSeqFile(store, linear.handle, linear.view.linearSeqs[0].uid, 'gb', picked);
  assert.equal(linear.view.linearSeqs[0].gb, genome, 'the other drawing binds the same ref');
  assert.equal(linear.view.linearSeqs[0].gb.fileId, genome.fileId);
  assert.deepEqual(sessionFilesOfKind(bindings, 'genbank'), [genome], 'one file, bound twice, is offered once');
  assert.equal(store.registerFile(fileNamed('later.gb')).fileId, 'file-4', 'choosing stored no new file');
});
