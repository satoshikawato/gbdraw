// Record discovery runs from the transitions that change a drawing's inputs
// (R10): each one calls the discovery of its own drawing, and the discovery
// writes only that drawing. A drawing switch calls none. Vue's `watch` is a
// no-op here, so a discovery seen below was started by the transition itself.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }), querySelector: () => null };
installFakeSvgDom();
globalThis.File = File;

const { state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');

const genbank = (name, ...records) => new File([records.map(([id, length]) => (
  `LOCUS       ${id} ${length} bp DNA\nACCESSION   ${id}\nORIGIN\n        1 ${'a'.repeat(length)}\n//\n`
)).join('')], name);
const settle = async () => { for (let i = 0; i < 20; i += 1) await new Promise((resolve) => setTimeout(resolve, 0)); };
const app = createAppSetup();
const circular = modeDrawing(state, 'circular');
const linear = modeDrawing(state, 'linear');
// The GenBank file of the drawing's one File card: an upload into the card (a new card while there is none),
// or, with null, a Remove (which removes the card holding no other file).
const setCircularGenbank = (file) => {
  const uid = circular.circularSources[0]?.uid;
  if (uid) app.setCircularSourceFile(uid, 'gb', file);
  else if (file) app.addCircularSourceFiles('gb', [file]);
};
const recordIds = (drawing) => drawing.circularRecordList.value.map((record) => record.record_id);

test('a Circular upload discovers that drawing only; a switch writes no record list', async () => {
  app.selectDrawing('linear');
  app.setLinearSeqPrimaryFile(0, 'gb', genbank('l.gb', ['L1', 30]));
  await settle();
  app.selectDrawing('circular');
  const linearRecordsBefore = linear.linearSeqs.map((seq) => app.linearRecordOptions(seq).map((option) => option.value));
  setCircularGenbank(genbank('c.gb', ['C1', 40], ['C2', 50]));
  await settle();
  assert.deepEqual(recordIds(circular), ['C1', 'C2'], 'the upload discovers its drawing\'s records');
  assert.deepEqual(linear.linearSeqs.map((seq) => app.linearRecordOptions(seq).map((option) => option.value)),
    linearRecordsBefore, 'the other drawing\'s records are as they were');
  const list = circular.circularRecordList.value;
  app.selectDrawing('linear');
  app.selectDrawing('circular');
  await settle();
  assert.equal(circular.circularRecordList.value, list, 'a drawing switch writes no record list');
  assert.equal(circular.circularRecordDiscovery.status, 'ready');
});

test('the input type change and a Circular Remove discover that drawing again', async () => {
  app.selectDrawing('circular');
  app.setCircularInputType('gff');
  await settle();
  assert.equal(circular.circularRecordDiscovery.status, 'idle', 'GFF without its pair: nothing to read');
  assert.deepEqual(recordIds(circular), []);
  app.setCircularInputType('gb');
  await settle();
  assert.deepEqual(recordIds(circular), ['C1', 'C2'], 'back to GenBank: its file is read again');
  setCircularGenbank(null);
  await settle();
  assert.equal(circular.circularRecordDiscovery.status, 'idle');
  assert.deepEqual(recordIds(circular), []);
});

// OV-1669: the idle record has one owner, so a History restore that runs the
// discovery again leaves the record it captured.
test('a discovery with nothing to read leaves the idle record of a new drawing', async () => {
  const { INPUTS_SLICE } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  const recordOf = (slice) => INPUTS_SLICE.capture(slice).circularRecordDiscovery;
  const idle = recordOf(INPUTS_SLICE.create().handle);
  app.selectDrawing('circular');
  setCircularGenbank(genbank('c.gb', ['C1', 40]));
  await settle();
  assert.equal(circular.circularRecordDiscovery.status, 'ready');
  app.setCircularInputType('gff');
  await settle();
  assert.deepEqual(recordOf(circular.slices.inputs), idle, 'GFF without its pair');
  app.setCircularInputType('gb');
  await settle();
  setCircularGenbank(null);
  await settle();
  assert.deepEqual(recordOf(circular.slices.inputs), idle, 'no file');
});

test('a Linear card upload, add, and input type change discover the Linear drawing', async () => {
  app.selectDrawing('linear');
  app.addLinearSeq();
  const added = linear.linearSeqs[linear.linearSeqs.length - 1];
  assert.equal(app.linearRecordOptions(added)[0].value, '', 'an empty card has no records');
  app.setLinearSeqPrimaryFile(linear.linearSeqs.length - 1, 'gb', genbank('l2.gb', ['L2', 20]));
  await settle();
  const uploaded = linear.linearSeqs.find((seq) => app.linearRecordOptions(seq).some((option) => option.value === 'L2'));
  assert.ok(uploaded, 'the upload read the card\'s records without an Inspect');
  app.setLinearInputType('gff');
  await settle();
  assert.ok(linear.linearSeqs.every((seq) => app.linearRecordSelectorDisabled(seq)), 'GFF cards without files read nothing');
  app.setLinearInputType('gb');
  await settle();
  assert.ok(linear.linearSeqs.every((seq) => !app.linearRecordSelectorDisabled(seq)), 'GenBank cards are read again');
  app.selectDrawing('circular');
});

test('no discovery starts while a Load rolls back', async () => {
  app.selectDrawing('circular');
  state.sessionImportRollbackInProgress.value = true;
  try {
    setCircularGenbank(genbank('r.gb', ['R1', 10]));
    await settle();
    assert.notEqual(circular.circularRecordDiscovery.status, 'ready', 'the rollback admits no discovery');
    assert.deepEqual(recordIds(circular), []);
  } finally {
    state.sessionImportRollbackInProgress.value = false;
  }
});

// OV-2489: a Load leaves a Circular record list deferred until Generate reads
// it; the drawing's record identities read it now, as they read a deferred
// Linear card, so New drawing from a region offers the record before a Generate.
test('the record identities read a Circular list still deferred after a Load', async () => {
  app.selectDrawing('circular');
  setCircularGenbank(genbank('d.gb', ['D1', 60]));
  await settle();
  assert.equal(circular.circularRecordDiscovery.status, 'ready');
  // The state a Session Load leaves: the list deferred for the current file.
  circular.circularRecordDiscovery.status = 'deferred';
  circular.circularRecordList.value = [];
  await app.regionDrawings.open();
  assert.deepEqual(app.regionDrawings.dialog.rows.map((row) => [row.caption, row.length]), [['D1', 60]]);
  assert.equal(circular.circularRecordDiscovery.status, 'ready', 'the list was read once, through its owner');
  app.regionDrawings.close();
});

// OV-2492 (D-15-205 A): the region dialog's read of a deferred list writes no
// setting. Record Order (`adv.multi_record_positions`, the request's
// `layout.multiRecordPositions`) stays as the Load left it until Generate or
// Reload Records, which apply the read.
test('the region dialog reads a deferred Circular list without writing Record Order', async () => {
  const { projectDraftRequest } = await import('../../gbdraw/web/js/services/config.js');
  const requestPositions = async () => (await projectDraftRequest(state, circular, {
    linearRecordCatalog: null, recordDisplayRows: [], omittedAnnotationRecordKeys: []
  })).envelope.renderRequest.layout.multiRecordPositions;
  app.selectDrawing('circular');
  setCircularGenbank(genbank('p.gb', ['P1', 70]));
  await settle();
  assert.equal(circular.circularRecordDiscovery.status, 'ready');
  // The state a Session Load without stored Record Order leaves.
  circular.circularRecordDiscovery.status = 'deferred';
  circular.circularRecordList.value = [];
  circular.adv.multi_record_positions.splice(0, circular.adv.multi_record_positions.length);
  const before = await requestPositions();
  await app.regionDrawings.open();
  assert.deepEqual(app.regionDrawings.dialog.rows.map((row) => row.caption), ['P1'], 'the dialog offers the record');
  app.regionDrawings.close();
  assert.deepEqual(circular.adv.multi_record_positions, [], 'the read writes no Record Order');
  assert.deepEqual(await requestPositions(), before, 'the request carries the Record Order the Load left');
  await app.refreshCircularRecordOrder();
  assert.deepEqual(circular.adv.multi_record_positions, [{ selector: circular.circularSources[0].uid, row: 1 }],
    'Reload Records applies the read, under the record\'s identity (the card\'s uid, a one-record file)');
});

// The region dialog's rows come from the record list the deferred read
// writes: the read finishes before the shown request is projected, so a
// Multi-Record Canvas projected after a Load offers every record, not only
// the ones of the empty list the read starts from.
test('the region dialog offers every record of a deferred Multi-Record Canvas', async () => {
  app.selectDrawing('circular');
  setCircularGenbank(genbank('m.gb', ['M1', 80], ['M2', 90]));
  await settle();
  assert.equal(circular.circularRecordDiscovery.status, 'ready');
  circular.form.multi_record_canvas = true;
  // The state a Session Load leaves: the list deferred for the current file.
  circular.circularRecordDiscovery.status = 'deferred';
  circular.circularRecordList.value = [];
  await app.regionDrawings.open();
  assert.deepEqual(app.regionDrawings.dialog.rows.map((row) => [row.caption, row.length]), [['M1', 80], ['M2', 90]]);
  app.regionDrawings.close();
  circular.form.multi_record_canvas = false;
});

// OV-2868: a record-major Depth matrix keeps one row per record of the cards.
// A card edit changes the records; the read after it moves each row with its
// record (a card's records keep their places in its file, D-20), so a Replace
// that changes a card's record count, a Move, and a Move made while a read
// still runs leave every record its own row.
test('a record-major Depth matrix follows its records through a Replace and a Move', async () => {
  const { bindDrawingInputFileTree } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  app.selectDrawing('circular');
  circular.form.multi_record_canvas = false;
  while (circular.circularSources.length > 0) await app.removeCircularSourceCard(circular.circularSources[0].uid);
  await app.addCircularSourceFiles('gb', [genbank('a.gb', ['A1', 40]), genbank('b.gb', ['B1', 50])]);
  await settle();
  assert.deepEqual(recordIds(circular), ['A1', 'B1']);
  const depth = (name) => state.inputFiles.registerFile(new File(['1\t1\n'], name));
  bindDrawingInputFileTree(state.inputFiles, circular.slices.inputs, 'c_depth', [[depth('a.tsv')], [depth('b.tsv')]]);
  const rows = () => circular.files.c_depth.map((row) => row.map((ref) => ref?.name ?? null));
  const [cardA, cardB] = circular.circularSources.map((card) => card.uid);
  await app.setCircularSourceFile(cardA, 'gb', genbank('a2.gb', ['A1', 40], ['A2', 45]));
  await settle();
  assert.deepEqual(recordIds(circular), ['A1', 'A2', 'B1']);
  assert.deepEqual(rows(), [['a.tsv'], [], ['b.tsv']], 'B keeps its row; A\'s first record keeps A\'s');
  // Three Moves, each made while the read of the one before still runs.
  await app.moveCircularSourceCard(1, -1);
  await app.moveCircularSourceCard(0, 1);
  await app.moveCircularSourceCard(1, -1);
  await settle();
  assert.deepEqual(circular.circularSources.map((card) => card.uid), [cardB, cardA]);
  assert.deepEqual(rows(), [['b.tsv'], ['a.tsv'], []]);
  // Undo gives each Move's matrix back with the records it stood for (the
  // test's Replace is no History step of its own: the file input makes it one).
  for (let step = 0; step < 3; step += 1) await app.undoHistory();
  await settle();
  assert.deepEqual(recordIds(circular), ['A1', 'A2', 'B1']);
  assert.deepEqual(rows(), [['a.tsv'], [], ['b.tsv']]);
  bindDrawingInputFileTree(state.inputFiles, circular.slices.inputs, 'c_depth', null);
});
