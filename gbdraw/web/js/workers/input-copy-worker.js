// The Worker that writes the input copies of `services/input-copies.js`:
// `FileSystemSyncAccessHandle` exists only in a Worker. Each copy file is
// open from `open` to `close`, or to `closeAll`, which also deletes the tab's
// folder when the page closes; `write` replies with the byte count `write()`
// returned, which the page checks. Nothing calls `flush()`: a copy is
// temporary, and it only has to outlive the page, not the system.
/** @import { InputCopyWorkerRequest } from '../services/input-copies.js' */

/** @type {Map<string, FileSystemSyncAccessHandle>} */
const handles = new Map();
/** @type {Map<string, number>} */
const sizes = new Map();

/** @param {InputCopyWorkerRequest} request @returns {Promise<number>} */
const run = async (request) => {
  if (request.op === 'open') {
    let folder = await navigator.storage.getDirectory();
    for (const name of request.folder) folder = await folder.getDirectoryHandle(name);
    const file = await folder.getFileHandle(request.name, { create: true });
    handles.set(request.name, await file.createSyncAccessHandle());
    sizes.set(request.name, 0);
    return 0;
  }
  if (request.op === 'closeAll') {
    // The page closes: every copy is closed, so the folder can be deleted.
    for (const open of handles.values()) open.close();
    handles.clear();
    sizes.clear();
    let parent = await navigator.storage.getDirectory();
    for (const name of request.folder.slice(0, -1)) parent = await parent.getDirectoryHandle(name);
    await parent.removeEntry(request.folder[request.folder.length - 1], { recursive: true });
    return 0;
  }
  const handle = handles.get(request.name);
  if (!handle) throw new Error(`No open copy ${request.name}.`);
  if (request.op === 'write') {
    const at = sizes.get(request.name) ?? 0;
    const written = handle.write(new Uint8Array(request.buffer), { at });
    sizes.set(request.name, at + written);
    return written;
  }
  handles.delete(request.name);
  sizes.delete(request.name);
  handle.close();
  return 0;
};

self.onmessage = async ({ data }) => {
  const request = /** @type {InputCopyWorkerRequest} */ (data);
  try {
    self.postMessage({ id: request.id, ok: true, written: await run(request) });
  } catch (error) {
    self.postMessage({ id: request.id, ok: false, error: String(error instanceof Error ? error.message : error) });
  }
};
