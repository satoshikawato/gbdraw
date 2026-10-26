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
    } else assert.ok(Buffer.byteLength(JSON.stringify(part)) <= 128 * 1024);
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
  receiver.receivePart({ kind: 'string-start', path: [0], encoding: 'utf-16' });
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
  await sendBoundedJson(source, async (message, buffers = []) => {
    // Rows of earlier, acknowledged messages are released (then the whole branch); rows of this message are not yet.
    assert.ok((source.rows || []).slice(0, acknowledgedRows).every(item => item === null));
    for (const part of message.kind === 'parts' ? message.parts : [message]) {
      if (part.kind !== 'batch' || part.path[0] !== 'rows') continue;
      assert.deepEqual(source.rows.slice(part.index, part.index + part.value.length), part.value);
      acknowledgedRows = part.index + part.value.length;
    }
    receiver.receivePart(structuredClone(message, { transfer: buffers }));
    await new Promise(resolve => setImmediate(resolve));
  }, { consume: true });
  assert.ok(acknowledgedRows > 0);
  assert.deepEqual(receiver.getValue(), expected);
  assert.deepEqual(source, {});
});

test('an unacknowledged owned result branch remains available when transport fails', async () => {
  const source = { rows: Array.from({ length: 4000 }, (_, index) => ({ index, text: 'x'.repeat(200) })) };
  const pending = source.rows[0];
  await assert.rejects(sendBoundedJson(source, async message => {
    if ((message.parts || [message]).some(part => part.kind === 'batch')) throw new Error('ACK failed');
  }, { consume: true }), /ACK failed/);
  assert.equal(source.rows[0], pending);
});

test('bounded Worker byte replies transfer exact canonical bytes in ordered chunks', async () => {
  const bytes = Uint8Array.from({ length: 800000 }, (_, index) => index % 251);
  const source = { canonicalResource: { kind: 'collinearity-result', bytes } };
  const receiver = createBoundedJsonReceiver();
  let transferred = 0;
  await sendBoundedJson(source, async (message, buffers = []) => {
    const chunks = (message.parts || [message]).filter(part => part.kind === 'bytes-chunk');
    assert.equal(buffers.length, chunks.length);
    for (const buffer of buffers) {
      assert.ok(buffer.byteLength <= 256 * 1024);
      transferred += buffer.byteLength;
    }
    receiver.receivePart(structuredClone(message, { transfer: buffers }));
  }, { consume: true });
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

/** Sends `value` through a clone of every message and returns the received value and the string chunks sent. */
const roundTrip = async (value) => {
  const receiver = createBoundedJsonReceiver();
  const chunks = [];
  await sendBoundedJson(value, async (message, buffers = []) => {
    for (const part of message.kind === 'parts' ? message.parts : [message]) {
      if (part.kind === 'string-start') chunks.push({ encoding: part.encoding, sizes: [] });
      if (part.kind === 'string-bytes') chunks.at(-1).sizes.push(part.bytes.byteLength);
      if (part.kind === 'string-chunk') chunks.at(-1).sizes.push(part.units.byteLength);
    }
    for (const buffer of buffers) assert.ok(buffer.byteLength <= 256 * 1024);
    receiver.receivePart(structuredClone(message, { transfer: buffers }));
  });
  return { value: receiver.getValue(), chunks };
};

test('long well-formed strings travel as UTF-8 byte chunks and keep every code unit', async () => {
  const cases = {
    ascii: 'ACGTN'.repeat(200000),
    bmp: 'αβγ漢字'.repeat(90000),
    astral: '😀𝔸'.repeat(70000),
    bom: '﻿start'.repeat(40000),
    // Each 3- and 4-byte character crosses a 256 KiB chunk boundary at some offset.
    shifted: `${'a'.repeat(262143)}漢${'b'.repeat(262142)}😀${'c'.repeat(100)}`,
  };
  for (const [name, text] of Object.entries(cases)) {
    const { value, chunks } = await roundTrip({ text });
    assert.equal(value.text, text, name);
    assert.equal(chunks.length, 1, name);
    assert.equal(chunks[0].encoding, 'utf-8', name);
    assert.equal(chunks[0].sizes.reduce((sum, size) => sum + size, 0), Buffer.byteLength(text), name);
    assert.ok(chunks[0].sizes.every(size => size > 0 && size <= 256 * 1024), name);
  }
});

test('long strings with lone surrogates keep the UTF-16 route; short and empty strings stay values', async () => {
  const lone = `${'x'.repeat(200000)}\ud800${'y'.repeat(10)}\udfff`;
  const { value, chunks } = await roundTrip({ lone, empty: '', short: 'é😀\ud800' });
  assert.equal(value.lone, lone);
  assert.equal(value.empty, '');
  assert.equal(value.short, 'é😀\ud800');
  assert.deepEqual(chunks.map(chunk => chunk.encoding), ['utf-16']);
});

test('without String.prototype.isWellFormed a long string travels the UTF-16 route and stays equal', async () => {
  const original = Object.getOwnPropertyDescriptor(String.prototype, 'isWellFormed');
  const text = 'ACGTN'.repeat(200000);
  try {
    delete String.prototype.isWellFormed;
    const { value, chunks } = await roundTrip({ text });
    assert.equal(value.text, text);
    assert.deepEqual(chunks.map(chunk => chunk.encoding), ['utf-16']);
  } finally {
    if (original) Object.defineProperty(String.prototype, 'isWellFormed', original);
  }
});

test('the receiver decodes a multi-byte sequence split across byte chunks and rejects malformed UTF-8', () => {
  const bytes = new TextEncoder().encode('﻿a漢😀');
  const receiver = createBoundedJsonReceiver();
  receiver.receivePart({ kind: 'value', path: [], value: {} });
  receiver.receivePart({ kind: 'string-start', path: ['text'], encoding: 'utf-8' });
  for (const [start, end] of [[0, 2], [2, 5], [5, 9], [9, bytes.length]]) {
    receiver.receivePart({ kind: 'string-bytes', bytes: bytes.subarray(start, end) });
  }
  receiver.receivePart({ kind: 'string-end' });
  assert.equal(receiver.getValue().text, '﻿a漢😀');

  const truncated = createBoundedJsonReceiver();
  truncated.receivePart({ kind: 'string-start', path: [], encoding: 'utf-8' });
  truncated.receivePart({ kind: 'string-bytes', bytes: bytes.subarray(0, 5) });
  assert.throws(() => truncated.receivePart({ kind: 'string-end' }), /string/);

  const mixed = createBoundedJsonReceiver();
  mixed.receivePart({ kind: 'string-start', path: [], encoding: 'utf-8' });
  assert.throws(() => mixed.receivePart({ kind: 'string-chunk', units: new Uint16Array([0x61]) }), /string/);
  assert.throws(() => mixed.receivePart({ kind: 'string-bytes', bytes: new Uint8Array([0xff]) }), /string/);
  const unknown = createBoundedJsonReceiver();
  assert.throws(() => unknown.receivePart({ kind: 'string-start', path: [], encoding: 'latin1' }), /string/);
  const units = createBoundedJsonReceiver();
  units.receivePart({ kind: 'string-start', path: [], encoding: 'utf-16' });
  assert.throws(() => units.receivePart({ kind: 'string-bytes', bytes: new Uint8Array([0x61]) }), /string/);
});
