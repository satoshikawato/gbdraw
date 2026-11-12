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
  chosenFileRefusal, contentTokenOfBytes, copyInputStream, holdFileCopy, holdReadBytes, keepFileCopy, readFileBytes, readFileContentFacts,
  releaseFileCopy, storedFileBytes, takeFileBytesForTransfer, useInputCopyStore
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
  /** An open file, or a folder with one inside: OPFS refuses to delete it while a sync access handle is open. */
  const holdsOpen = (entry) => Boolean(entry.open || [...(entry.entries?.values() ?? [])].some(holdsOpen));
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
      return {
        file,
        getFile: async () => {
          copyReads += 1;
          return new Blob([file.bytes]);
        }
      };
    },
    async removeEntry(name, { recursive = false } = {}) {
      const entry = node.entries.get(name);
      if (!entry) throw notFound();
      if (holdsOpen(entry)) throw new DOMException('open', 'NoModificationAllowedError');
      if (entry.entries?.size && !recursive) throw new DOMException('not empty', 'InvalidModificationError');
      node.entries.delete(name);
    },
    async *keys() { yield* [...node.entries.keys()]; }
  });
  const open = new Map();
  let copyReads = 0;
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
          (written) => worker.onmessage(request.op === 'write'
            // A write gives its buffer back, as the copy Worker does.
            ? structuredClone({ data: { id: request.id, ok: true, written, buffer: request.buffer } }, { transfer: [request.buffer] })
            : { data: { id: request.id, ok: true, written } }),
          (error) => worker.onmessage({ data: { id: request.id, ok: false, error: String(error) } })
        );
      }
    };
    return worker;
  };
  const tabFolders = () => [...(root.entries.get('gbdraw-input-copies')?.entries.keys() ?? [])];
  const copyFiles = () => [...(root.entries.get('gbdraw-input-copies')?.entries.values() ?? [])]
    .flatMap((tab) => [...tab.entries.values()]);
  return {
    root, storage: { getDirectory: async () => handleOf(root) }, createWorker, tabFolders, copyFiles, copyReads: () => copyReads
  };
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

/**
 * `promise`'s value; rejects when it has not settled within `ms` of real time. A gzip
 * file expands through a stream, which takes more than a few event-loop rounds.
 */
const within = (promise, ms = 10_000) => {
  let timer;
  const late = new Promise((_, reject) => {
    timer = setTimeout(() => reject(new Error(`not settled within ${ms} ms`)), ms);
  });
  return Promise.race([promise, late]).finally(() => clearTimeout(timer));
};
/** Waits, in real time up to `ms`, until `check()` holds. */
const until = async (check, ms = 10_000) => {
  const end = Date.now() + ms;
  while (!check()) {
    if (Date.now() > end) throw new Error(`not reached within ${ms} ms`);
    await new Promise((resolve) => setTimeout(resolve, 5));
  }
};

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

test('a gzip input is copied as its gzip bytes with their facts, reads expanded, and the cap refuses it when it is chosen', async () => {
  const fileSystem = createFakeFileSystem();
  useInputCopyStore(copyStoreOn(fileSystem));
  const plain = new TextEncoder().encode('chr\t1\t5\n'.repeat(512));
  const gzip = new Uint8Array(gzipSync(plain));
  const file = chosenFile(gzip, 'depth.tsv.gz');
  await choose(file);
  // The copy holds the bytes as chosen, which a Session file stores (receipt session.file-container).
  assert.deepEqual(fileSystem.copyFiles()[0].bytes, gzip);
  assert.deepEqual(await readFileContentFacts(file), { token: await contentTokenOfBytes(gzip), byteLength: gzip.byteLength },
    'the facts are those of the bytes as chosen');
  assert.deepEqual(await readFileBytes(file), plain);
  // A later reading reads the copy and expands it again, even after the file changed on disk.
  file.refused = true;
  await flush();
  assert.deepEqual(await readFileBytes(file), plain);
  assert.equal(file.reads, 1);
  file.refused = false;
  // A plain File of the same content holds other bytes, so it gets its own copy.
  await choose(chosenFile(plain, 'depth.tsv'));
  assert.deepEqual(fileSystem.copyFiles().map(({ bytes }) => bytes.byteLength), [gzip.byteLength, plain.byteLength]);

  const metrics = [];
  globalThis.__GBDRAW_TEST_HOOKS__ = { maxExpandedInputBytes: 1024, onStructuralMetric: (metric) => metrics.push(metric) };
  const counted = (name) => metrics.filter((metric) => metric.name === name).length;
  try {
    // PD-OI-116: a gzip file that expands past the cap is refused right after it is chosen, whatever its kind. Its
    // one read counts the expanded bytes while the copy is written, keeps nothing expanded, and stops at the cap.
    const large = chosenFile(new Uint8Array(gzipSync(Buffer.concat([plain, plain]))), 'large.tsv.gz');
    const refused = (error) => (
      error.code === 'INPUT_UNREADABLE' && error.context.reason === 'GZIP_EXPANDED_SIZE'
      && error.context.inputName === 'large.tsv.gz' && error.context.limit === 0
    );
    keepFileCopy(large);
    assert.ok(refused(await within(chosenFileRefusal(large))), 'the uploader shows the refusal of its read');
    await flush();
    assert.equal(large.reads, 1, 'choosing reads the file once');
    assert.deepEqual([counted('inputExpansionCheckCount'), counted('inputExpansionCount')], [1, 0],
      'it counts the expanded bytes once and keeps none');
    assert.equal(fileSystem.copyFiles().length, 2, 'the refused file keeps no copy');
    // Every later reading of it, as Generate's, gives the same refusal.
    await assert.rejects(readFileContentFacts(large), refused);
    await assert.rejects(readFileBytes(large), refused);
    // A field that asks before it binds the file (another drawing uses the file it
    // replaces) asks for the refusal before the input store holds the file: it is
    // answered once the store holds it and that read ends, right after the answer.
    const asked = chosenFile(new Uint8Array(gzipSync(Buffer.concat([plain, plain]))), 'asked.tsv.gz');
    const answer = chosenFileRefusal(asked);
    await flush();
    assert.equal(asked.reads, 0, 'asking reads nothing before the store holds the file');
    keepFileCopy(asked);
    const refusal = await within(answer);
    assert.deepEqual([refusal.context.reason, refusal.context.inputName], ['GZIP_EXPANDED_SIZE', 'asked.tsv.gz'],
      'the refusal comes once the store holds the file');
    assert.equal(asked.reads, 1, 'the hold reads the file once');
    // A gzip file within the cap is admitted when it is chosen, and copied as chosen.
    const small = chosenFile(new Uint8Array(gzipSync(plain.subarray(0, 1000))), 'small.tsv.gz');
    keepFileCopy(small);
    assert.equal(await within(chosenFileRefusal(small)), null);
    await flush();
    assert.equal(fileSystem.copyFiles().length, 3, 'the admitted file is copied as chosen');
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

test('a first read whose bytes are detached before they are copied keeps its facts\' length, and its copy fails', async () => {
  const fileSystem = createFakeFileSystem();
  useInputCopyStore(copyStoreOn(fileSystem));
  const bytes = new TextEncoder().encode('LOCUS detached\n//\n');
  const file = chosenFile(bytes);
  keepFileCopy(file);
  const read = await readFileBytes(file);
  // The first read's bytes go elsewhere before the copy writes them.
  structuredClone(read.buffer, { transfer: [read.buffer] });
  assert.equal((await readFileContentFacts(file)).byteLength, bytes.byteLength);
  await until(() => fileSystem.copyFiles().every(({ open }) => !open));
  await flush();
  assert.deepEqual(fileSystem.copyFiles(), [], 'the short copy is deleted');
  assert.deepEqual(await readFileBytes(file), bytes, 'the File is read again');
});

test('a copy whose source is detached after its first piece holds fewer bytes than the facts count, and fails', async () => {
  let writes = 0;
  /** @type {Promise<Uint8Array> | null} */
  let first = null;
  /** @type {{ reads: number } | null} */
  let file = null;
  const fileSystem = createFakeFileSystem({
    // While the first piece is written, the facts are hashed and the source bytes go elsewhere.
    writeGate: async () => {
      writes += 1;
      if (writes !== 1) return;
      const bytes = await first;
      await readFileContentFacts(file);
      structuredClone(bytes.buffer, { transfer: [bytes.buffer] });
    }
  });
  useInputCopyStore(copyStoreOn(fileSystem));
  const content = new Uint8Array(PIECE + 1024).map((_, index) => index % 251);
  file = chosenFile(content, 'two-pieces.gb');
  keepFileCopy(file);
  first = readFileBytes(file);
  assert.equal((await within(readFileContentFacts(file))).byteLength, content.byteLength);
  await until(() => writes >= 1 && fileSystem.copyFiles().every(({ open }) => !open));
  await flush();
  assert.equal(writes, 1, 'the loop stopped at the detached source');
  // Lengths and a comparison, not deepEqual: a failure message of 16 MiB arrays would take minutes.
  assert.deepEqual(fileSystem.copyFiles().map(({ bytes }) => bytes.byteLength), [], 'no truncated copy is kept');
  const again = await readFileBytes(file);
  assert.ok(Buffer.from(again).equals(Buffer.from(content)), 'the File is read again');
  assert.equal(file.reads, 2);
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
  // The gate stays closed until both readings resolved: they cannot wait for a copy byte.
  for (const [file, content] of files) {
    keepFileCopy(file);
    assert.deepEqual(await within(readFileBytes(file)), content, `${file.name}: the bytes wait for no copy`);
    assert.ok((await within(readFileContentFacts(file))).byteLength > 0, `${file.name}: nor the facts`);
  }
  assert.deepEqual(fileSystem.copyFiles().map(({ bytes }) => bytes.byteLength), [0, 0], 'the copies are still being written');
  gate.open();
  await until(() => fileSystem.copyFiles().every(({ bytes, open }) => bytes.byteLength && !open));
  assert.deepEqual(fileSystem.copyFiles().map(({ bytes }) => bytes), [plain, new Uint8Array(gzipSync(expanded))]);
});

test('a gzip file whose copy could not be written keeps nothing expanded: a later reading reads the File again, checked', async () => {
  const fileSystem = createFakeFileSystem({ writeResult: () => 0 });
  useInputCopyStore(copyStoreOn(fileSystem));
  const plain = new TextEncoder().encode('chr\t1\t9\n'.repeat(64));
  const gzip = new Uint8Array(gzipSync(plain));
  const file = chosenFile(gzip, 'unwritten.tsv.gz');
  await choose(file);
  assert.deepEqual(await readFileBytes(file), plain);
  assert.deepEqual(fileSystem.copyFiles(), []);
  await flush();
  assert.equal(file.reads, 2, 'with no copy, a later reading reads the File again');
  assert.deepEqual(await readFileBytes(file), plain, 'and expands its bytes again');
  assert.equal(file.reads, 3, 'nothing expanded was kept');
  // The File changed on disk: its gzip bytes no longer match the facts of its first read.
  file.change(new Uint8Array(gzipSync(new TextEncoder().encode('other'))));
  await flush();
  await assert.rejects(readFileBytes(file), (error) => error.code === 'INPUT_UNREADABLE' && !error.context.reason);
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

test('a reading that starts while the copy is written takes the first read\'s bytes, and a transfer then reads no copy', async () => {
  const gate = createGate();
  const fileSystem = createFakeFileSystem({ writeGate: gate.wait });
  useInputCopyStore(copyStoreOn(fileSystem));
  const file = chosenFile(new Uint8Array(PIECE + 9).map((_, index) => index % 249));
  keepFileCopy(file);
  const first = await readFileBytes(file);
  const release = holdReadBytes();
  try {
    assert.equal(await readFileBytes(file), first, 'a Generate reading shares the bytes the page holds');
    const taken = takeFileBytesForTransfer(file);
    assert.equal(await settledSoon(taken), TIMED_OUT, 'a transfer waits for the copy to close');
    gate.open();
    assert.equal(await taken, first.buffer, 'it transfers the bytes the page holds');
    assert.equal(fileSystem.copyReads(), 0, 'and reads no copy');
  } finally {
    release();
  }
  assert.equal(file.reads, 1);
});

test('a copy goes through one buffer, which the Worker gives back after each piece', async () => {
  const fileSystem = createFakeFileSystem();
  const store = copyStoreOn(fileSystem);
  /** @type {{ given: ArrayBuffer[], back: ArrayBuffer[] }} */
  const seen = { given: [], back: [] };
  useInputCopyStore({
    openWriter: async () => {
      const writer = await store.openWriter();
      return writer && {
        ...writer,
        write: async (buffer) => {
          seen.given.push(buffer);
          const back = await writer.write(buffer);
          seen.back.push(back);
          return back;
        }
      };
    }
  });
  const bytes = new Uint8Array(3 * PIECE + 5).map((_, index) => index % 241);
  await choose(chosenFile(bytes));
  assert.deepEqual(fileSystem.copyFiles()[0].bytes, bytes);
  assert.equal(seen.given.length, 4);
  assert.equal(seen.given[1], seen.back[0]);
  assert.equal(seen.given[2], seen.back[1]);
});

test('a failed Feature Edits TSV import keeps a copy, which its Retry reads after the file changed on disk', async () => {
  const { createFeatureEditTableActions } = await import('../../gbdraw/web/js/app/feature-editor/feature-edit-table.js');
  const fileSystem = createFakeFileSystem();
  useInputCopyStore(copyStoreOn(fileSystem));
  const table = new TextEncoder().encode('record_id\tfeature_id\tlabel_text\nr1\tf1\tEdited\n');
  const file = chosenFile(table, 'feature_edits.tsv');
  const sent = [];
  const errorLog = { value: null };
  const actions = createFeatureEditTableActions({
    state: { activeDrawing: () => ({}), errorLog, sessionOperationAvailability: () => null },
    ref: (value) => ({ value }),
    computed: (getter) => ({ get value() { return getter(); } }),
    getCommittedSession: () => ({ renderRequest: { records: [{ recordKey: 'r1' }] }, resources: {} }),
    readResourceRecordCount: async () => 1,
    // The helper fails each time: the test reads what each attempt sent.
    readFeatureOverrideTable: async ({ files }) => {
      sent.push(new Uint8Array(files[0].bytes));
      throw new Error('helper failed');
    },
    projectFeatureEdits: async () => {}
  });
  assert.equal(await actions.loadFeatureEditTable({ target: { files: [file], value: 'feature_edits.tsv' } }), false);
  assert.equal(actions.canRetryFeatureEditTableImport.value, true);
  await until(() => fileSystem.copyFiles().some(({ bytes, open }) => bytes.byteLength && !open));
  file.refused = true;
  assert.equal(await actions.retryFeatureEditTableImport(), false);
  assert.deepEqual(sent, [table, table], 'Retry sends the content first read, from its copy');
  assert.equal(file.reads, 1);
  assert.notEqual(errorLog.value?.code, 'INPUT_UNREADABLE');
});

/** A stream of `bytes` in chunks of `size` bytes (a ZIP entry inflates in small chunks); `fail` ends it with an error after them. */
const chunkedStream = (bytes, size, fail = null) => new ReadableStream({
  start(controller) {
    for (let start = 0; start < bytes.byteLength; start += size) controller.enqueue(bytes.slice(start, start + size));
    if (fail) controller.error(fail);
    else controller.close();
  }
});

test('what a Session file stores for a File streams from its copy, not from memory, with the size and token of its first read', async () => {
  const fileSystem = createFakeFileSystem();
  useInputCopyStore(copyStoreOn(fileSystem));
  const bytes = new Uint8Array(PIECE + 17).map((_, index) => (index * 7) % 251);
  const file = chosenFile(bytes);
  await choose(file);
  const stored = await storedFileBytes(file);
  assert.equal(stored.size, bytes.byteLength);
  assert.equal(await stored.token(), await contentTokenOfBytes(bytes));
  const readsBefore = fileSystem.copyReads();
  assert.ok(Buffer.from(await new Response(await stored.open()).arrayBuffer()).equals(Buffer.from(bytes)));
  assert.equal(fileSystem.copyReads(), readsBefore + 1, 'the copy is streamed');
  assert.equal(file.reads, 1, 'the File is not read again');
});

test('without a copy, what a Session file stores for a File is the File, checked as it streams against its first read', async () => {
  const bytes = new TextEncoder().encode('LOCUS stored\n//\n');
  const gzip = new Uint8Array(gzipSync(bytes));
  /** A real File's shape: it streams and has a size. */
  const streamingFile = (content, name) => {
    let current = content;
    return {
      name, get size() { return current.byteLength; }, async arrayBuffer() { return current.slice().buffer; },
      stream() { return new Blob([current]).stream(); }, change(next) { current = next; }
    };
  };
  // Never read (no copy store): its own size, and no token until a stream of it is hashed.
  const unread = streamingFile(gzip, 'unread.gb.gz');
  const first = await storedFileBytes(unread);
  assert.deepEqual([first.size, first.token], [gzip.byteLength, null]);
  assert.deepEqual(new Uint8Array(await new Response(await first.open()).arrayBuffer()), gzip);
  // Read before: the gzip bytes as chosen, with the facts of that read; a change on disk since refuses the stream.
  const read = streamingFile(gzip, 'read.gb.gz');
  assert.deepEqual(await readFileBytes(read), bytes);
  const stored = await storedFileBytes(read);
  assert.deepEqual([stored.size, await stored.token()], [gzip.byteLength, await contentTokenOfBytes(gzip)]);
  assert.deepEqual(new Uint8Array(await new Response(await stored.open()).arrayBuffer()), gzip);
  const changed = new Uint8Array(gzip);
  changed[changed.length - 1] ^= 1;
  read.change(changed);
  await assert.rejects(new Response(await stored.open()).arrayBuffer(), (error) => error.code === 'INPUT_UNREADABLE');
  read.change(new Uint8Array(3));
  await assert.rejects(stored.open(), (error) => error.code === 'INPUT_UNREADABLE');
});

test('a stream is copied in 16 MiB pieces through one buffer, under the token its bytes have', async () => {
  const fileSystem = createFakeFileSystem();
  const store = copyStoreOn(fileSystem);
  const written = [];
  useInputCopyStore({
    openWriter: async () => {
      const writer = await store.openWriter();
      return writer && { ...writer, write: async (buffer) => { written.push(buffer.byteLength); return writer.write(buffer); } };
    }
  });
  for (const length of [2 * PIECE + 1000, PIECE, 0]) {
    written.length = 0;
    const bytes = new Uint8Array(length).map((_, index) => (index * 13) % 249);
    const copied = await within(copyInputStream(chunkedStream(bytes, 65536 + 3)));
    assert.ok(copied.copy, `copied ${length} bytes`);
    assert.deepEqual(copied.facts, { token: await contentTokenOfBytes(bytes), byteLength: length });
    assert.deepEqual(new Uint8Array(await copied.copy.read()), bytes);
    assert.deepEqual(new Uint8Array(await new Response((await copied.copy.blob()).stream()).arrayBuffer()), bytes);
    assert.deepEqual(written, length === 0 ? [0] : length === PIECE ? [PIECE] : [PIECE, PIECE, 1000]);
    copied.copy.release();
  }
  assert.equal(fileSystem.copyFiles().length, 0, 'a released copy is deleted');
});

test('a stream with the content of a chosen File shares its copy', async () => {
  const fileSystem = createFakeFileSystem();
  useInputCopyStore(copyStoreOn(fileSystem));
  const bytes = new Uint8Array(1000).map((_, index) => index % 7);
  await choose(chosenFile(bytes));
  const copied = await within(copyInputStream(chunkedStream(bytes, 100)));
  assert.ok(copied.copy);
  assert.equal(fileSystem.copyFiles().length, 1);
});

test('a stream that fails, or a copy write that fails or falls short, leaves no copy; no store cancels the stream', async () => {
  const fileSystem = createFakeFileSystem();
  useInputCopyStore(copyStoreOn(fileSystem));
  const bytes = new Uint8Array(PIECE + 10).map((_, index) => index % 5);
  const damaged = Object.assign(new Error('damaged'), { code: 'INPUT_INVALID' });
  await assert.rejects(within(copyInputStream(chunkedStream(bytes, 1 << 20, damaged))), (error) => error === damaged);
  await until(() => fileSystem.copyFiles().length === 0);

  const full = createFakeFileSystem({ writeResult: (buffer) => Math.max(0, buffer.byteLength - 1) });
  useInputCopyStore(copyStoreOn(full));
  assert.deepEqual(await within(copyInputStream(chunkedStream(bytes, 1 << 20))), { copy: null, reason: 'write-failed' });
  await until(() => full.copyFiles().length === 0);

  useInputCopyStore(null);
  let cancelled = false;
  const unread = new ReadableStream({ pull() {}, cancel() { cancelled = true; } });
  assert.deepEqual(await copyInputStream(unread), { copy: null, reason: 'unavailable' });
  assert.equal(cancelled, true);
});
