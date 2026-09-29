import assert from 'node:assert/strict';
import test from 'node:test';
import { createBoundedJsonReceiver, sendBoundedJson } from '../../gbdraw/web/js/services/bounded-json-transport.js';

test('bounded JSON transfers preserve every value, code unit and own unsafe key', async () => {
  const source = JSON.parse('{"__proto__":{"marker":true}}');
  source.text = 'α😀\ud800\udfff'.repeat(60000);
  source.rows = Array.from({ length: 4000 }, (_, index) => ({ index, text: `row-${index}` }));
  source.empty = [null, {}, [], false, 0, ''];
  const expected = JSON.stringify(source);
  const receiver = createBoundedJsonReceiver();
  let count = 0;
  let transfers = 0;
  await sendBoundedJson(source, async (part, buffers = []) => {
    count++;
    if (buffers.length) {
      assert.equal(buffers.length, 1);
      assert.ok(buffers[0].byteLength <= 256 * 1024);
      transfers++;
    } else assert.ok(Buffer.byteLength(JSON.stringify(part)) < 129 * 1024);
    const reply = structuredClone(part, { transfer: buffers });
    for (const buffer of buffers) assert.equal(buffer.byteLength, 0);
    receiver.receivePart(reply);
    // Acknowledgement completes on a later task, before the next part.
    await new Promise((resolve) => setImmediate(resolve));
  });
  assert.ok(count > 3 && transfers > 1);
  assert.equal(JSON.stringify(receiver.getValue()), expected);
  assert.equal(JSON.stringify(source), expected);
  assert.equal(Object.hasOwn(receiver.getValue(), '__proto__'), true);
  assert.equal(Object.getPrototypeOf(receiver.getValue()), Object.prototype);
});

test('incomplete strings and out of order array batches cannot expose a candidate', () => {
  const receiver = createBoundedJsonReceiver();
  assert.throws(() => receiver.getValue(), /Incomplete/);
  receiver.receivePart({ kind: 'value', path: [], value: [] });
  assert.throws(() => receiver.receivePart({ kind: 'batch', path: [], index: 1, value: [1] }), /batch/);
  receiver.receivePart({ kind: 'string-start', path: [0] });
  assert.throws(() => receiver.getValue(), /Incomplete/);
  assert.throws(() => receiver.receivePart({ kind: 'string-chunk', units: new Uint8Array(2) }), /string/);
});

test('array batches cannot walk inherited transport paths', () => {
  const receiver = createBoundedJsonReceiver();
  receiver.receivePart({ kind: 'value', path: [], value: [] });
  const prototypeLength = Array.prototype.length;
  try {
    assert.throws(() => receiver.receivePart({ kind: 'batch', path: ['__proto__'], index: prototypeLength, value: ['unexpected'] }), /parent/);
  } finally { Array.prototype.length = prototypeLength; }
  assert.equal(Array.prototype.length, prototypeLength);
  assert.deepEqual(receiver.getValue(), []);
});

test('owned Worker results release acknowledged branches while preserving the full received graph', async () => {
  const source = JSON.parse('{"__proto__":{"marker":true}}');
  source.rows = Array.from({ length: 4000 }, (_, index) => ({ index, text: 'α😀'.repeat(80) }));
  source.text = 'α😀\ud800'.repeat(60000);
  const expected = structuredClone(source);
  const receiver = createBoundedJsonReceiver();
  let acknowledgedRows = 0;
  await sendBoundedJson(source, async (part, buffers = []) => {
    if (part.kind === 'batch' && part.path[0] === 'rows') {
      assert.ok(source.rows.slice(0, acknowledgedRows).every(item => item === null));
      assert.deepEqual(source.rows.slice(part.index, part.index + part.value.length), part.value);
      acknowledgedRows = part.index + part.value.length;
    }
    receiver.receivePart(structuredClone(part, { transfer: buffers }));
    await new Promise(resolve => setImmediate(resolve));
  }, [], { consume: true });
  assert.ok(acknowledgedRows > 0);
  assert.deepEqual(receiver.getValue(), expected);
  assert.deepEqual(source, {});
});

test('an unacknowledged owned result branch remains available when transport fails', async () => {
  const source = { rows: Array.from({ length: 4000 }, (_, index) => ({ index, text: 'x'.repeat(200) })) };
  const pending = source.rows[0];
  await assert.rejects(sendBoundedJson(source, async part => {
    if (part.kind === 'batch') throw new Error('ACK failed');
  }, [], { consume: true }), /ACK failed/);
  assert.equal(source.rows[0], pending);
});

test('bounded Worker byte replies transfer exact canonical bytes in ordered chunks', async () => {
  const bytes = Uint8Array.from({ length: 800000 }, (_, index) => index % 251);
  const source = { canonicalResource: { kind: 'collinearity-result', bytes } };
  const receiver = createBoundedJsonReceiver();
  let transferred = 0;
  await sendBoundedJson(source, async (part, buffers = []) => {
    if (part.kind === 'bytes-chunk') {
      assert.equal(buffers.length, 1);
      assert.ok(buffers[0].byteLength <= 256 * 1024);
      transferred += buffers[0].byteLength;
    }
    receiver.receivePart(structuredClone(part, { transfer: buffers }));
  }, [], { consume: true });
  assert.equal(transferred, bytes.byteLength);
  assert.deepEqual(receiver.getValue().canonicalResource.bytes, bytes);
  assert.deepEqual(source, {});
});

test('bounded Worker byte replies reject incomplete or out-of-order chunks', () => {
  const receiver = createBoundedJsonReceiver();
  receiver.receivePart({ kind: 'bytes-start', path: [], length: 3 });
  assert.throws(() => receiver.getValue(), /Incomplete/);
  assert.throws(() => receiver.receivePart({ kind: 'bytes-chunk', offset: 1,
    bytes: new Uint8Array([1]) }), /Invalid bounded bytes chunk/);
  receiver.receivePart({ kind: 'bytes-chunk', offset: 0, bytes: new Uint8Array([1, 2]) });
  assert.throws(() => receiver.receivePart({ kind: 'bytes-end' }), /Incomplete bounded bytes/);
});
