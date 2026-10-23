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
const { retainHistoryAndBoundFiles } = await import('../../gbdraw/web/js/services/drawing/inputs.js');

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
    fileId: ref.fileId, name: 'genome.fasta', type: 'text/plain', size: file.size, lastModified: 7
  });
  assert.equal(ref instanceof Blob, false);
  assert.equal(typeof ref.arrayBuffer, 'undefined', 'a byte reader must resolve the ref with fileOf');
  assert.equal(store.registerFile(ref), null, 'a ref is never stored as a file');
  assert.equal(store.fileOf(ref), file);
});
