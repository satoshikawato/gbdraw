// @ts-check
import {
  GBDRAW_WHEEL_CACHE_BUST,
  GBDRAW_WHEEL_NAME,
  PYODIDE_INDEX_URL,
  PYODIDE_LOCAL_WHEELS,
  PYODIDE_RUNTIME_FILES,
  PYODIDE_RUNTIME_WHEELS,
  PYODIDE_SNAPSHOT_URL
} from '../config.js';

const getDefaultBaseUrl = () => globalThis.window?.location?.href || globalThis.location?.href || './';

/**
 * @param {string} path
 * @param {string} [baseUrl]
 */
export const resolveAssetUrl = (path, baseUrl = getDefaultBaseUrl()) =>
  new URL(path, baseUrl).toString();

export const resolvePyodideIndexUrl = (baseUrl = getDefaultBaseUrl()) =>
  resolveAssetUrl(PYODIDE_INDEX_URL, baseUrl);

export const resolvePyodideModuleUrl = (baseUrl = getDefaultBaseUrl()) =>
  new URL('pyodide.mjs', resolvePyodideIndexUrl(baseUrl)).toString();

export const resolveLocalDependencyWheelUrls = (baseUrl = getDefaultBaseUrl()) =>
  PYODIDE_LOCAL_WHEELS.map((path) => resolveAssetUrl(path, baseUrl));

export const resolveGbdrawWheelUrl = (baseUrl = getDefaultBaseUrl()) => {
  const wheelUrl = new URL(GBDRAW_WHEEL_NAME, baseUrl);
  if (GBDRAW_WHEEL_CACHE_BUST) {
    wheelUrl.searchParams.set('v', GBDRAW_WHEEL_CACHE_BUST);
  }
  return wheelUrl.toString();
};

// Whether this page has asked for the diagram Worker's boot files: by the
// prefetch, or by a Worker boot (its manifest).
let bootAssetsRequested = false;

/**
 * The diagram Worker's boot manifest (its `init` message). The Worker then
 * downloads the boot files itself, so a later prefetch on this page is skipped.
 * @param {string} [baseUrl]
 */
export const buildPyodideAssetManifest = (baseUrl = getDefaultBaseUrl()) => {
  bootAssetsRequested = true;
  return {
    pyodideIndexUrl: resolvePyodideIndexUrl(baseUrl),
    pyodideModuleUrl: resolvePyodideModuleUrl(baseUrl),
    gbdrawWheelUrl: resolveGbdrawWheelUrl(baseUrl),
    localWheelUrls: resolveLocalDependencyWheelUrls(baseUrl),
    pyodideSnapshotUrl: PYODIDE_SNAPSHOT_URL ? resolveAssetUrl(PYODIDE_SNAPSHOT_URL, baseUrl) : null
  };
};

/**
 * Every file the diagram Worker's cold boot downloads: the Pyodide module and
 * runtime files, the distribution wheels it loads by name, the dependency
 * wheels and the gbdraw wheel.
 * @param {string} [baseUrl]
 * @returns {string[]}
 */
export const resolvePyodideBootAssetUrls = (baseUrl = getDefaultBaseUrl()) => {
  const indexUrl = resolvePyodideIndexUrl(baseUrl);
  return [
    resolvePyodideModuleUrl(baseUrl),
    ...PYODIDE_RUNTIME_FILES.map((file) => new URL(file, indexUrl).toString()),
    ...PYODIDE_RUNTIME_WHEELS.map((path) => resolveAssetUrl(path, baseUrl)),
    ...resolveLocalDependencyWheelUrls(baseUrl),
    resolveGbdrawWheelUrl(baseUrl)
  ];
};

/** @typedef {{ connection?: { saveData?: boolean } }} SaveDataNavigator */

// Two downloads at a time leave the browser's other connections to the host
// (six over HTTP/1.1) free for the files the page loads meanwhile.
const BOOT_ASSET_PREFETCH_LANES = 2;

/**
 * Downloads the diagram Worker's boot files into the HTTP cache once per page,
 * without starting the Worker, so a later first Generate reads them from the
 * cache. The requests use the Worker's URLs at low priority, which the HTTP
 * cache does not key on. A lane stops at its first network error; the Worker
 * requests a file the cache lacks again. Nothing is downloaded when the
 * browser asks to save data, or once a diagram Worker has started booting.
 */
export const prefetchPyodideBootAssets = () => {
  if (bootAssetsRequested) return;
  bootAssetsRequested = true;
  if (/** @type {SaveDataNavigator} */ (globalThis.navigator).connection?.saveData === true) return;
  const pending = resolvePyodideBootAssetUrls();
  const next = () => {
    const url = pending.shift();
    if (!url) return;
    fetch(url, { priority: 'low' })
      .then((response) => response.body?.pipeTo(new WritableStream()))
      .then(next, () => {});
  };
  for (let lane = 0; lane < BOOT_ASSET_PREFETCH_LANES; lane += 1) next();
};
