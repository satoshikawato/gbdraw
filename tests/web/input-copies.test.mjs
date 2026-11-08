// The temporary copies of the chosen input files (PD-OI-121,
// services/input-copies.js and services/file-content-cache.js): a File the
// input store holds is copied by its first read, every later reading reads the
// copy, a File without a copy is checked against its first read's facts, and a
// copy is deleted when no File uses it. The private file system, the locks and
// the copy Worker are in-memory stand-ins that follow the browser's rules.
import assert from 'node:assert/strict';
import test from 'node:test';
import { webcrypto } from 'node:crypto';
import { gzipSync } from 'node:zlib';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const {
  holdFileCopy, holdReadBytes, keepFileCopy, readFileBytes, readFileContentFacts, releaseFileCopy, useInputCopyStore
} = await import('../../gbdraw/web/js/services/file-content-cache.js');
const { createBrowserInputCopyStore } = await import('../../gbdraw/web/js/services/input-copies.js');

const PIECE = 16 * 1024 * 1024;
const notFound = () => new DOMException('missing', 'NotFoundError');

/**
 * An origin private file system in memory: folders of files, and the copy Worker that writes them. Each write
 * waits for `writeGate()` first, which stands for a slow Worker.
 */
const createFakeFileSystem = ({ writeResult = (buffer) => buffer.byteLength, writeGate = async () => {} } = {}) => {
  const folder = () => ({ entries: new Map() });
  const root = folder();
  const handleOf = (node) => ({
    async getDirectoryHandle(name, { create = false } = {}) {
      if (!node.entries.has(name)) {
        if (!create) throw notFound();
        node.entries.set(name, folder());
      }
      return handleOf(node.entries.get(name));
    },
    async getFileHandle(name, { create = false } = {}) {
      if (!node.entries.has(name)) {
        if (!create) throw notFound();
        node.entries.set(name, { bytes: new Uint8Array(0) });
      }
      const file = node.entries.get(name);
      return { file, getFile: async () => new Blob([file.bytes]) };
    },
    async removeEntry(name, { recursive = false } = {}) {
      const entry = node.entries.get(name);
      if (!entry) throw notFound();
      if (entry.open) throw new DOMException('open', 'NoModificationAllowedError');
      if (entry.entries?.size && !recursive) throw new DOMException('not empty', 'InvalidModificationError');
      node.entries.delete(name);
    },
    async *keys() { yield* [...node.entries.keys()]; }
  });
  const open = new Map();
  const run = async ({ op, folder: path, name, buffer }) => {
    if (op === 'closeAll') {
      for (const file of open.values()) file.open = false;
      open.clear();
      let parent = handleOf(root);
      for (const part of path.slice(0, -1)) parent = await parent.getDirectoryHandle(part);
      await parent.removeEntry(path.at(-1), { recursive: true });
      return 0;
    }
    if (op === 'open') {
      let handle = handleOf(root);
      for (const part of path) handle = await handle.getDirectoryHandle(part);
      const { file } = await handle.getFileHandle(name, { create: true });
      file.open = true;
      open.set(name, file);
      return 0;
    }
    const file = open.get(name);
    if (op === 'close') {
      file.open = false;
      open.delete(name);
      return 0;
    }
    await writeGate();
    if (!open.has(name)) throw new Error('closed');
    const written = writeResult(buffer);
    const bytes = new Uint8Array(buffer).subarray(0, Math.min(written, buffer.byteLength));
    const joined = new Uint8Array(file.bytes.byteLength + bytes.byteLength);
    joined.set(file.bytes);
    joined.set(bytes, file.bytes.byteLength);
    file.bytes = joined;
    return written;
  };
  const createWorker = () => {
    const worker = {
      onmessage: null,
      onerror: null,
      terminate() {},
      postMessage(message, transfer = []) {
        // A posted message detaches the buffers it transfers, as in the browser.
        const request = structuredClone(message, { transfer });
        run(request).then(
          (written) => worker.onmessage({ data: { id: request.id, ok: true, written } }),
          (error) => worker.onmessage({ data: { id: request.id, ok: false, error: String(error) } })
        );
      }
    };
    return worker;
  };
  const tabFolders = () => [...(root.entries.get('gbdraw-input-copies')?.entries.keys() ?? [])];
  const copyFiles = () => [...(root.entries.get('gbdraw-input-copies')?.entries.values() ?? [])]
    .flatMap((tab) => [...tab.entries.values()]);
  return { root, storage: { getDirectory: async () => handleOf(root) }, createWorker, tabFolders, copyFiles };
};

/** Web Locks in one process: a held name refuses `ifAvailable`, and a callback's end releases it. */
const createFakeLocks = () => {
  const held = new Set();
  return {
    held,
    async request(name, options, callback) {
      const [settings, run] = callback ? [options, callback] : [{}, options];
      if (held.has(name)) {
        if (settings.ifAvailable) return run(null);
        throw new Error(`${name} is held`);
      }
      held.add(name);
      try {
        return await run({ name });
      } finally {
        held.delete(name);
      }
    }
  };
};

const flush = async () => {
  for (let round = 0; round < 20; round += 1) await new Promise((resolve) => setImmediate(resolve));
};

/** A chosen File whose reads are counted; `change` stands for the file changed on disk after it was chosen. */
const chosenFile = (bytes, name = 'chosen.gb') => {
  let current = bytes;
  return {
    name,
    reads: 0,
    refused: false,
    async arrayBuffer() {
      this.reads += 1;
      if (this.refused) throw new DOMException('changed', 'NotReadableError');
      return current.slice().buffer;
    },
    change(next) { current = next; }
  };
};

const copyStoreOn = (fileSystem, locks = createFakeLocks(), page = new EventTarget()) => createBrowserInputCopyStore({
  storage: fileSystem.storage, locks, createWorker: fileSystem.createWorker, page
});
/** A gate a test opens: writes wait until then. */
const createGate = () => {
  let open;
  const opened = new Promise((resolve) => { open = resolve; });
  return { wait: () => opened, open };
};
const TIMED_OUT = Symbol('timed out');
/** `promise`'s value, or TIMED_OUT when it has not settled after the event loop ran a while. */
const settledSoon = async (promise) => Promise.race([promise, flush().then(() => TIMED_OUT)]);

/** Chooses `file` (the input store holds it) and waits until its copy is written. */
const choose = async (file) => {
  keepFileCopy(file);
  await readFileContentFacts(file);
  await flush();
};

test.afterEach(() => useInputCopyStore(null));

test('a chosen File is copied by its one read, and later readings read the copy even after the file changed', async () => {
  const fileSystem = createFakeFileSystem();
  useInputCopyStore(copyStoreOn(fileSystem));
  const bytes = new Uint8Array(PIECE + 5).map((_, index) => index % 253);
  const file = chosenFile(bytes);
  await choose(file);
  assert.equal(file.reads, 1, 'the copy is written from the first read');
  assert.equal(fileSystem.copyFiles().length, 1);
  assert.deepEqual(fileSystem.copyFiles()[0].bytes, bytes);
  // Chromium refuses a changed File; Firefox returns its new bytes. Neither is read again.
  file.refused = true;
  assert.deepEqual(await readFileBytes(file), bytes);
  file.refused = false;
  file.change(new Uint8Array([1, 2, 3]));
  await flush();
  assert.deepEqual(await readFileBytes(file), bytes);
  assert.equal(file.reads, 1);
});

test('a gzip input is copied as its expanded content, and the cap still refuses it', async () => {
  const fileSystem = createFakeFileSystem();
  useInputCopyStore(copyStoreOn(fileSystem));
  const plain = new TextEncoder().encode('chr\t1\t5\n'.repeat(512));
  const file = chosenFile(new Uint8Array(gzipSync(plain)), 'depth.tsv.gz');
  await choose(file);
  assert.deepEqual(fileSystem.copyFiles()[0].bytes, plain);
  assert.equal((await readFileContentFacts(file)).byteLength, plain.byteLength, 'the facts count the expanded bytes');
  assert.deepEqual(await readFileBytes(file), plain);
  assert.equal(file.reads, 1);

  globalThis.__GBDRAW_TEST_HOOKS__ = { maxExpandedInputBytes: 1024 };
  try {
    const large = chosenFile(new Uint8Array(gzipSync(plain)), 'large.tsv.gz');
    keepFileCopy(large);
    await assert.rejects(readFileBytes(large), (error) => (
      error.code === 'INPUT_UNREADABLE' && error.context.reason === 'GZIP_EXPANDED_SIZE' && error.context.limit === 0
    ));
    await flush();
    assert.equal(fileSystem.copyFiles().length, 1, 'the refused file leaves no copy');
  } finally {
    delete globalThis.__GBDRAW_TEST_HOOKS__;
  }
});

test('a copy that cannot be written falls back to reading the File, which refuses once it changed', async () => {
  // Chromium's in-memory storage returned 4294967288 from one write().
  const fileSystem = createFakeFileSystem({ writeResult: () => 4294967288 });
  useInputCopyStore(copyStoreOn(fileSystem));
  const bytes = new TextEncoder().encode('LOCUS x\n//\n');
  const file = chosenFile(bytes);
  await choose(file);
  assert.deepEqual(fileSystem.copyFiles(), [], 'the partial copy is deleted');
  assert.deepEqual(await readFileBytes(file), bytes, 'an unchanged File reads as before');
  await flush();
  file.change(new TextEncoder().encode('LOCUS y\n//\n'));
  await assert.rejects(readFileBytes(file), (error) => error.code === 'INPUT_UNREADABLE' && error.stage === 'read');
});

test('without a copy, a File that returns other bytes than its first read is refused (Firefox)', async () => {
  const bytes = new TextEncoder().encode('LOCUS f\n//\n');
  const longer = chosenFile(bytes);
  await readFileBytes(longer);
  await readFileContentFacts(longer);
  await flush();
  longer.change(new TextEncoder().encode('LOCUS f\n//\n\n'));
  await assert.rejects(readFileBytes(longer), (error) => error.code === 'INPUT_UNREADABLE');

  const sameLength = chosenFile(bytes);
  await readFileBytes(sameLength);
  await readFileContentFacts(sameLength);
  await flush();
  sameLength.change(new TextEncoder().encode('LOCUS g\n//\n'));
  await assert.rejects(readFileBytes(sameLength), (error) => error.code === 'INPUT_UNREADABLE');
});

test('a copy is deleted when the input store lets its File go, after the readings that run end', async () => {
  const fileSystem = createFakeFileSystem();
  useInputCopyStore(copyStoreOn(fileSystem));
  const file = chosenFile(new TextEncoder().encode('>x\nACGT\n'));
  await choose(file);
  const release = holdReadBytes();
  releaseFileCopy(file);
  await flush();
  assert.equal(fileSystem.copyFiles().length, 1, 'a Generate or Save that runs keeps the copy');
  release();
  await flush();
  assert.deepEqual(fileSystem.copyFiles(), []);
});

test('Files with the same content share one copy until the last of them is let go', async () => {
  const fileSystem = createFakeFileSystem();
  useInputCopyStore(copyStoreOn(fileSystem));
  const bytes = new TextEncoder().encode('>same\nACGT\n');
  const first = chosenFile(bytes);
  const second = chosenFile(bytes);
  await choose(first);
  await choose(second);
  assert.equal(fileSystem.copyFiles().length, 1);
  releaseFileCopy(first);
  await flush();
  assert.equal(fileSystem.copyFiles().length, 1);
  assert.deepEqual(await readFileBytes(second), bytes);
  releaseFileCopy(second);
  await flush();
  assert.deepEqual(fileSystem.copyFiles(), []);
});

test('the start sweep deletes the folders of tabs that are gone and keeps a live tab\'s', async () => {
  const fileSystem = createFakeFileSystem();
  const locks = createFakeLocks();
  const base = await (await fileSystem.storage.getDirectory()).getDirectoryHandle('gbdraw-input-copies', { create: true });
  await base.getDirectoryHandle('tab-gone', { create: true });
  await (await base.getDirectoryHandle('tab-live', { create: true })).getFileHandle('copy-1.bin', { create: true });
  locks.request('gbdraw-input-copies/tab-live', () => new Promise(() => {}));
  const store = copyStoreOn(fileSystem, locks);
  await store.openWriter();
  await flush();
  const folders = fileSystem.tabFolders();
  assert.ok(!folders.includes('tab-gone'));
  assert.ok(folders.includes('tab-live'));
  const own = folders.find((name) => name !== 'tab-live');
  assert.ok(locks.held.has(`gbdraw-input-copies/${own}`), 'this tab holds its own folder\'s lock');
});

test('a browser without a private file system or locks writes no copies', () => {
  assert.equal(createBrowserInputCopyStore({ storage: {}, locks: createFakeLocks() }), null);
  assert.equal(createBrowserInputCopyStore({ storage: createFakeFileSystem().storage, locks: {} }), null);
});

test('a first reading and its facts resolve before a slow copy is written, for a plain and a gzip file', async () => {
  const gate = createGate();
  const fileSystem = createFakeFileSystem({ writeGate: gate.wait });
  useInputCopyStore(copyStoreOn(fileSystem));
  const plain = new TextEncoder().encode('LOCUS slow\n//\n');
  const expanded = new TextEncoder().encode('LOCUS slow-gzip\n//\n');
  const files = [[chosenFile(plain), plain], [chosenFile(new Uint8Array(gzipSync(expanded)), 'slow.gb.gz'), expanded]];
  for (const [file, content] of files) {
    keepFileCopy(file);
    assert.deepEqual(await settledSoon(readFileBytes(file)), content, `${file.name}: the bytes wait for no copy`);
    assert.equal((await settledSoon(readFileContentFacts(file))).byteLength, content.byteLength, `${file.name}: nor the facts`);
  }
  assert.deepEqual(fileSystem.copyFiles().map(({ bytes }) => bytes.byteLength), [0, 0], 'the copies are still being written');
  gate.open();
  await flush();
  assert.deepEqual(fileSystem.copyFiles().map(({ bytes }) => bytes), [plain, expanded]);
});

test('a gzip file whose copy could not be written keeps its expanded bytes: no second read and no new expansion', async () => {
  const fileSystem = createFakeFileSystem({ writeResult: () => 0 });
  useInputCopyStore(copyStoreOn(fileSystem));
  const plain = new TextEncoder().encode('chr\t1\t9\n'.repeat(64));
  const file = chosenFile(new Uint8Array(gzipSync(plain)), 'unwritten.tsv.gz');
  await choose(file);
  const first = await readFileBytes(file);
  assert.deepEqual(fileSystem.copyFiles(), []);
  await flush();
  assert.equal(await readFileBytes(file), first, 'later readings share the expanded bytes');
  assert.equal(file.reads, 1);
});

test('a tab whose lock is refused makes no folder and writes no copy', async () => {
  const fileSystem = createFakeFileSystem();
  const store = copyStoreOn(fileSystem, { request: async () => { throw new DOMException('refused', 'AbortError'); } });
  assert.equal(await store.openWriter(), null);
  useInputCopyStore(store);
  const file = chosenFile(new TextEncoder().encode('>refused\nACGT\n'));
  await choose(file);
  assert.deepEqual(fileSystem.tabFolders(), []);
});

test('closing the page closes the copies being written and deletes the tab\'s folder', async () => {
  const gate = createGate();
  const fileSystem = createFakeFileSystem({ writeGate: gate.wait });
  const page = new EventTarget();
  useInputCopyStore(copyStoreOn(fileSystem, createFakeLocks(), page));
  const file = chosenFile(new TextEncoder().encode('>closing\nACGT\n'));
  keepFileCopy(file);
  await readFileBytes(file);
  await flush();
  assert.equal(fileSystem.tabFolders().length, 1);
  page.dispatchEvent(new Event('pagehide'));
  await flush();
  assert.deepEqual(fileSystem.tabFolders(), []);
  gate.open();
});

test('a copy that another state still uses outlives the input store\'s release', async () => {
  const fileSystem = createFakeFileSystem();
  useInputCopyStore(copyStoreOn(fileSystem));
  const bytes = new TextEncoder().encode('CDS\tproduct\tx\t#123456\tx\n');
  const file = chosenFile(bytes, 'rules.tsv');
  await choose(file);
  const endUse = holdFileCopy(file);
  releaseFileCopy(file);
  await flush();
  file.refused = true;
  assert.deepEqual(await readFileBytes(file), bytes, 'a Retry reads the copy');
  await flush();
  assert.equal(fileSystem.copyFiles().length, 1);
  endUse();
  await flush();
  assert.deepEqual(fileSystem.copyFiles(), []);
});
