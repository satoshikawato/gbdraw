// @ts-check
import {
  GBDRAW_WHEEL_CACHE_BUST,
  GBDRAW_WHEEL_NAME,
  PYODIDE_INDEX_URL,
  PYODIDE_LOCAL_WHEELS,
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

export const buildPyodideAssetManifest = (baseUrl = getDefaultBaseUrl()) => ({
  pyodideIndexUrl: resolvePyodideIndexUrl(baseUrl),
  pyodideModuleUrl: resolvePyodideModuleUrl(baseUrl),
  gbdrawWheelUrl: resolveGbdrawWheelUrl(baseUrl),
  localWheelUrls: resolveLocalDependencyWheelUrls(baseUrl),
  pyodideSnapshotUrl: PYODIDE_SNAPSHOT_URL ? resolveAssetUrl(PYODIDE_SNAPSHOT_URL, baseUrl) : null
});
