// The vendored Pyodide in Node with its heap pushed past 2 GiB, for tests of
// values that Python allocates above 2^31 and that cross to or from JavaScript.
// The heap never shrinks in a Worker, so after one large Generate every later
// value may sit there. The blocks that push the heap are allocated but never
// touched, so the process stays small.
import assert from 'node:assert/strict';
import { fileURLToPath } from 'node:url';
import { PYODIDE_INDEX_URL, PYODIDE_LOCAL_WHEELS } from '../../../gbdraw/web/js/config.js';

export const HIGH_ADDRESS = 2 ** 31;

const webRoot = new URL('../../../gbdraw/web/', import.meta.url);

/** The vendored Pyodide at the folder config.js names (written from gbdraw/_web_assets.py). */
export const loadVendoredPyodide = async () => {
  globalThis.self ??= {};
  const indexURL = new URL(PYODIDE_INDEX_URL, webRoot);
  const { loadPyodide } = await import(new URL('pyodide.mjs', indexURL).href);
  return loadPyodide({ indexURL: fileURLToPath(indexURL) });
};

/**
 * The vendored Pyodide with the wheels the Worker installs and this checkout's gbdraw sources on sys.path
 * (mounted read-only from the repository), for tests that run gbdraw's Python where intp is 32 bits.
 */
export const loadVendoredGbdraw = async () => {
  const pyodide = await loadVendoredPyodide();
  await pyodide.loadPackage(PYODIDE_LOCAL_WHEELS.map((wheel) => fileURLToPath(new URL(wheel, webRoot))), {
    messageCallback: () => {}
  });
  pyodide.FS.mkdirTree('/gbdraw-src');
  pyodide.FS.mount(pyodide.FS.filesystems.NODEFS, { root: fileURLToPath(new URL('../../', webRoot)) }, '/gbdraw-src');
  pyodide.runPython('import sys; sys.path.insert(0, "/gbdraw-src")');
  return pyodide;
};

/**
 * Push the heap past 2 GiB, then take every free block below it of at least
 * 512 bytes, so that each new allocation larger than that lands above 2 GiB.
 * @returns {() => void} Frees the blocks.
 */
export const allocateAboveTwoGiB = (pyodide) => {
  const module = pyodide._module;
  const blocks = [];
  const take = (size) => {
    const block = module._PyMem_RawMalloc(size) >>> 0;
    assert.notEqual(block, 0, 'the Pyodide heap could not grow');
    blocks.push(block);
    return block;
  };
  const tail = pyodide.runPython('lambda: id(b"x" * (1 << 20))');
  try {
    while (tail() < HIGH_ADDRESS) {
      assert.ok(blocks.length < 4, 'the Pyodide heap did not grow past 2 GiB');
      take(2 ** 30 - 2 ** 20);
    }
  } finally {
    tail.destroy();
  }
  for (const size of [2 ** 20, 2 ** 16, 2 ** 12, 2 ** 9]) {
    while (take(size) < HIGH_ADDRESS) { /* a free block below 2 GiB taken */ }
  }
  return () => blocks.forEach((block) => module._PyMem_RawFree(block));
};
