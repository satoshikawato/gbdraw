// @ts-check
// The temporary copies of the input files the user chose (PD-OI-121): each
// lives in this tab's folder of the browser's private site storage (the origin
// private file system), never leaves the browser, and is deleted when no File
// uses it, when the page closes, or, if the page ended without that, when the
// page next opens. A Worker writes them with `FileSystemSyncAccessHandle`; a
// main-thread `createWritable()` wrote 7 to 11 times slower in Chromium.
import { testFailsInputCopyWrites } from './runtime-test-hooks.js';

const COPY_ROOT = 'gbdraw-input-copies';
/** @param {string} folder */
const lockNameOf = (folder) => `${COPY_ROOT}/${folder}`;

/**
 * One copy file. Files with equal content (one content token) share it.
 * @typedef {object} InputCopy
 * @property {() => Promise<Uint8Array>} read The copied bytes; rejects if the file is gone.
 * @property {() => void} release Ends one File's use of the copy; the last use deletes the file. A second call does nothing.
 */
/**
 * @typedef {object} InputCopyWriter
 * @property {(buffer: ArrayBuffer) => Promise<void>} write Appends `buffer`, which it transfers; rejects unless every byte was written.
 * @property {(token: string) => Promise<InputCopy>} finish Ends the writing; a copy that holds `token` already is shared instead.
 * @property {() => void} abort Ends the writing and deletes the file.
 */
/**
 * Where the copies are written. `openWriter` resolves to null when the browser
 * lets no copy be written (no private file system, no Worker, no lock).
 * @typedef {object} InputCopyStore
 * @property {() => Promise<InputCopyWriter | null>} openWriter
 */
/**
 * A call of the copy Worker (`workers/input-copy-worker.js`): open a copy file in a folder (its path
 * from the private file system's root), append to it, or close it; or, when the page closes, close
 * every copy and delete the folder.
 * @typedef {{ op: 'open', folder: string[], name: string }
 *   | { op: 'write', name: string, buffer: ArrayBuffer }
 *   | { op: 'close', name: string }
 *   | { op: 'closeAll', folder: string[] }} InputCopyWorkerCall
 */
/** @typedef {InputCopyWorkerCall & { id: number }} InputCopyWorkerRequest */
/**
 * The copy Worker's reply: `written` is the byte count `write()` returned.
 * @typedef {{ id: number, ok: true, written: number } | { id: number, ok: false, error: string }} InputCopyWorkerReply
 */
/** @typedef {{ keys: () => AsyncIterable<string> }} ListedFolder */
/** @typedef {{ base: FileSystemDirectoryHandle, folder: FileSystemDirectoryHandle }} CopyFolders */

/**
 * Starts a Worker and returns its call; a crashed Worker rejects every pending call and the next call starts another.
 * @param {() => Worker} createWorker
 */
const workerCaller = (createWorker) => {
  /** @type {Worker | null} */
  let worker = null;
  /** @type {Map<number, { resolve: (written: number) => void, reject: (error: Error) => void }>} */
  const waiting = new Map();
  let nextId = 1;
  /** @param {Error} error */
  const failAll = (error) => {
    worker?.terminate();
    worker = null;
    for (const { reject } of waiting.values()) reject(error);
    waiting.clear();
  };
  const start = () => {
    const started = createWorker();
    started.onmessage = ({ data }) => {
      const reply = /** @type {InputCopyWorkerReply} */ (data);
      const pending = waiting.get(reply.id);
      waiting.delete(reply.id);
      if (reply.ok) pending?.resolve(reply.written);
      else pending?.reject(new Error(reply.error));
    };
    started.onerror = () => failAll(new Error('The input copy Worker stopped.'));
    return started;
  };
  /**
   * @param {InputCopyWorkerCall} request
   * @param {Transferable[]} [transfer]
   * @returns {Promise<number>}
   */
  const call = (request, transfer = []) => new Promise((resolve, reject) => {
    const id = nextId++;
    waiting.set(id, { resolve, reject });
    try {
      worker ??= start();
      worker.postMessage({ ...request, id }, transfer);
    } catch (error) {
      waiting.delete(id);
      reject(error instanceof Error ? error : new Error(String(error)));
    }
  });
  return { call, started: () => worker !== null, stop: () => failAll(new Error('The page closed.')) };
};

/**
 * Deletes the folders of tabs that are gone: a live tab holds its folder's
 * lock, so a folder whose lock is free is a leftover.
 * @param {FileSystemDirectoryHandle} base
 * @param {string} own
 * @param {LockManager} locks
 */
const sweepLeftovers = async (base, own, locks) => {
  /** @type {string[]} */
  const names = [];
  for await (const name of /** @type {ListedFolder} */ (/** @type {unknown} */ (base)).keys()) names.push(name);
  for (const name of names.filter((entry) => entry !== own)) {
    await locks.request(lockNameOf(name), { ifAvailable: true }, async (lock) => {
      if (lock) await base.removeEntry(name, { recursive: true }).catch(() => {});
    });
  }
};

/**
 * Takes this tab's lock, then makes its folder: a folder exists only under a
 * held lock, so another tab's sweep never deletes it. A refused lock makes no
 * folder, and the tab writes no copies.
 * @param {StorageManager} storage
 * @param {LockManager} locks
 * @param {string} own
 * @returns {Promise<CopyFolders>}
 */
const openTabFolder = async (storage, locks, own) => {
  await new Promise((granted, refused) => {
    locks.request(lockNameOf(own), () => {
      granted(undefined);
      return new Promise(() => {});
    }).catch(refused);
  });
  const base = await (await storage.getDirectory()).getDirectoryHandle(COPY_ROOT, { create: true });
  return { base, folder: await base.getDirectoryHandle(own, { create: true }) };
};

/** @typedef {{ addEventListener: (type: 'pagehide', listener: (event: Event) => void) => void }} ClosingPage */

/**
 * The browser's copy store, or null when it has no private file system or no locks.
 * Creating it takes this tab's lock and folder and starts the sweep of leftovers.
 * When the page closes (not into the back-forward cache), the Worker closes the
 * copies and deletes the folder, and the page deletes it too if the Worker did
 * not; the next page's sweep deletes what is left.
 * @param {{ storage?: StorageManager, locks?: LockManager, createWorker?: () => Worker, page?: ClosingPage }} [environment]
 * @returns {InputCopyStore | null}
 */
export const createBrowserInputCopyStore = ({
  storage = globalThis.navigator?.storage,
  locks = globalThis.navigator?.locks,
  createWorker = () => new Worker(new URL('../workers/input-copy-worker.js', import.meta.url), { type: 'module' }),
  page = globalThis
} = {}) => {
  if (typeof storage?.getDirectory !== 'function' || typeof locks?.request !== 'function') return null;
  const own = `tab-${crypto.randomUUID()}`;
  const { call, started, stop } = workerCaller(createWorker);
  /** @type {Promise<CopyFolders | null>} */
  const folders = openTabFolder(storage, locks, own).then((opened) => {
    sweepLeftovers(opened.base, own, locks).catch(() => {});
    page.addEventListener?.('pagehide', (event) => {
      if (/** @type {PageTransitionEvent} */ (event).persisted) return;
      const closed = started() ? call({ op: 'closeAll', folder: [COPY_ROOT, own] }) : Promise.resolve(0);
      closed.catch(() => 0).then(() => {
        stop();
        opened.base.removeEntry(own, { recursive: true }).catch(() => {});
      });
    });
    return opened;
  }, () => null);
  /** @type {Map<string, { name: string, uses: number }>} */
  const byToken = new Map();
  let nextCopy = 1;

  /**
   * @param {FileSystemDirectoryHandle} folder
   * @param {string} token
   * @param {{ name: string, uses: number }} entry
   * @returns {InputCopy}
   */
  const useOf = (folder, token, entry) => {
    let used = true;
    return {
      read: async () => new Uint8Array(await (await (await folder.getFileHandle(entry.name)).getFile()).arrayBuffer()),
      release: () => {
        if (!used) return;
        used = false;
        entry.uses -= 1;
        if (entry.uses > 0) return;
        if (byToken.get(token) === entry) byToken.delete(token);
        folder.removeEntry(entry.name).catch(() => {});
      }
    };
  };

  return {
    openWriter: async () => {
      const opened = await folders;
      if (!opened) return null;
      const name = `copy-${nextCopy++}.bin`;
      const remove = () => { opened.folder.removeEntry(name).catch(() => {}); };
      try {
        await call({ op: 'open', folder: [COPY_ROOT, own], name });
      } catch {
        remove();
        return null;
      }
      return {
        write: async (buffer) => {
          const expected = buffer.byteLength;
          const written = await call({ op: 'write', name, buffer }, [buffer]);
          // A write can return less than it was given, or an error code (Chromium's in-memory storage, 4294967288).
          if (written !== expected || testFailsInputCopyWrites()) throw new Error(`Wrote ${written} of ${expected} bytes.`);
        },
        finish: async (token) => {
          await call({ op: 'close', name });
          const shared = byToken.get(token);
          if (shared) {
            remove();
            shared.uses += 1;
            return useOf(opened.folder, token, shared);
          }
          const entry = { name, uses: 1 };
          byToken.set(token, entry);
          return useOf(opened.folder, token, entry);
        },
        abort: () => { call({ op: 'close', name }).then(remove, remove); }
      };
    }
  };
};
