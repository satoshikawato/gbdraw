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
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');

const genbank = (name, ...records) => new File([records.map(([id, length]) => (
  `LOCUS       ${id} ${length} bp DNA\nACCESSION   ${id}\nORIGIN\n        1 ${'a'.repeat(length)}\n//\n`
)).join('')], name);
const settle = async () => { for (let i = 0; i < 20; i += 1) await new Promise((resolve) => setTimeout(resolve, 0)); };
const app = createAppSetup();
const circular = modeDrawing(state, 'circular');
const linear = modeDrawing(state, 'linear');
const recordIds = (drawing) => drawing.circularRecordList.value.map((record) => record.record_id);

test('a Circular upload discovers that drawing only; a switch writes no record list', async () => {
  app.selectDrawing('linear');
  app.setLinearSeqPrimaryFile(0, 'gb', genbank('l.gb', ['L1', 30]));
  await settle();
  app.selectDrawing('circular');
  const linearRecordsBefore = linear.linearSeqs.map((seq) => app.linearRecordOptions(seq).map((option) => option.value));
  app.setCircularSourceFile('c_gb', genbank('c.gb', ['C1', 40], ['C2', 50]));
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
  app.setCircularSourceFile('c_gb', null);
  await settle();
  assert.equal(circular.circularRecordDiscovery.status, 'idle');
  assert.deepEqual(recordIds(circular), []);
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
    app.setCircularSourceFile('c_gb', genbank('r.gb', ['R1', 10]));
    await settle();
    assert.notEqual(circular.circularRecordDiscovery.status, 'ready', 'the rollback admits no discovery');
    assert.deepEqual(recordIds(circular), []);
  } finally {
    state.sessionImportRollbackInProgress.value = false;
  }
});
