// The diagram Worker's font registry (custom fonts PR-2, review R1 and R5):
// a Generate that names a font the Worker has not registered fails with
// FONT_NOT_LOADED, a call registers all of its fonts or none, and the replay
// into a new Worker registers font by font and names each font it lost
// (OV-686).
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import test from 'node:test';

import { createWorkerFontRegistry } from '../../gbdraw/web/js/workers/worker-font-registry.js';

const fontOf = (text, label = `${text}.ttf`) => {
  const bytes = new TextEncoder().encode(text).buffer;
  return { identity: createHash('sha256').update(Buffer.from(bytes)).digest('hex'), label, bytes };
};

// A Worker file system and a Python stand-in that refuses files whose bytes start with "bad".
const fakeRuntime = () => {
  const files = new Map();
  const directories = new Set();
  /** @type {string[][]} */
  const calls = [];
  const fs = {
    analyzePath: (path) => ({ exists: directories.has(path) || files.has(path) }),
    mkdir: (path) => { directories.add(path); },
    writeFile: (path, bytes) => { files.set(path, new TextDecoder().decode(bytes)); },
    unlink: (path) => { files.delete(path); }
  };
  const registerInPython = (entries) => {
    calls.push(entries.map((entry) => entry.identity));
    if (entries.some((entry) => files.get(entry.path).startsWith('bad'))) {
      throw { code: 'INPUT_INVALID', stage: 'helper', context: { reason: 'FONT_FILE' } };
    }
    return { fonts: entries.map((entry) => ({ identity: entry.identity })) };
  };
  return { files, calls, registry: createWorkerFontRegistry({ fs, registerInPython }) };
};

test('a Generate that names an unregistered font fails with FONT_NOT_LOADED', async () => {
  const { registry } = fakeRuntime();
  const font = fontOf('regular');
  assert.throws(
    () => registry.requireRegistered([font.identity]),
    (error) => error.code === 'INPUT_INVALID' && error.context.reason === 'FONT_NOT_LOADED'
  );
  await registry.register([font]);
  assert.doesNotThrow(() => registry.requireRegistered([font.identity]));
  assert.doesNotThrow(() => registry.requireRegistered([]));
});

test('a call registers all of its fonts or none, and keeps no file of a refused call', async () => {
  const { registry, files } = fakeRuntime();
  const good = fontOf('good');
  const bad = fontOf('bad font');
  await assert.rejects(registry.register([good, bad]), (error) => error.code === 'INPUT_INVALID');
  assert.equal(files.size, 0);
  assert.throws(() => registry.requireRegistered([good.identity]), (error) => error.context.reason === 'FONT_NOT_LOADED');
});

test('bytes that do not match their identity are refused', async () => {
  const { registry } = fakeRuntime();
  await assert.rejects(
    registry.register([{ ...fontOf('one'), identity: fontOf('other').identity }]),
    (error) => error.code === 'RESOURCE_INVALID'
  );
});

test('OV-686: the replay registers font by font and names each font it lost', async () => {
  const { registry, calls } = fakeRuntime();
  const regular = fontOf('regular', 'Serif Regular.ttf');
  const broken = fontOf('bad bold', 'Serif Bold.ttf');
  const italic = fontOf('italic', 'Serif Italic.ttf');
  const failures = await registry.replay([regular, broken, italic]);
  assert.deepEqual(calls, [[regular.identity], [broken.identity], [italic.identity]]);
  assert.deepEqual(failures, [{ identity: broken.identity, label: 'Serif Bold.ttf' }]);
  assert.doesNotThrow(() => registry.requireRegistered([regular.identity, italic.identity]));
  assert.throws(() => registry.requireRegistered([broken.identity]), (error) => error.context.reason === 'FONT_NOT_LOADED');
});
