// The shipped LOSAT Wasm file a URL names, as a stubbed `fetch` serves it: the page compiles only the bytes its build
// records for that file (services/losat-runtime.js `compileLosatWasm`), so a test that runs the page's LOSAT loader
// serves the shipped file by its name, whatever the test's Workers then do with the module.
import { readFileSync } from 'node:fs';

const wasmRoot = new URL('../../../gbdraw/web/wasm/losat/', import.meta.url);
/** @type {Map<string, Uint8Array>} */
const files = new Map();

/**
 * The shipped bytes of the file the URL's last path segment names (`losat.wasm` or `losat-threaded.wasm`).
 * @param {string | URL} url
 * @returns {Uint8Array}
 */
export const shippedLosatWasm = (url) => {
  const name = new URL(String(url), 'https://example.test/').pathname.split('/').pop();
  if (name !== 'losat.wasm' && name !== 'losat-threaded.wasm') throw new Error(`no shipped LOSAT file named ${name}`);
  if (!files.has(name)) files.set(name, new Uint8Array(readFileSync(new URL(name, wasmRoot))));
  return /** @type {Uint8Array} */ (files.get(name));
};

/**
 * A `fetch` that answers the shipped LOSAT file a URL names.
 * @param {string | URL} url
 */
export const fetchShippedLosatWasm = async (url) => {
  const bytes = shippedLosatWasm(url);
  return { ok: true, arrayBuffer: async () => bytes.slice().buffer };
};
