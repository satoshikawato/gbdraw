// The diagram Worker's boot files (services/pyodide-assets.js): the list the
// page downloads into the HTTP cache when a file is chosen, before any Worker
// exists, and that prefetch's rules: once per page, two at a time, never on
// Save-Data, and a failed download is dropped without an unhandled rejection.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import {
  GBDRAW_WHEEL_NAME, PYODIDE_INDEX_URL, PYODIDE_LOCAL_WHEELS, PYODIDE_RUNTIME_FILES, PYODIDE_RUNTIME_WHEELS
} from '../../gbdraw/web/js/config.js';

const BASE = 'http://test.invalid/app/';
const ASSETS = new URL('../../gbdraw/web/js/services/pyodide-assets.js', import.meta.url).href;
// A fresh module instance: a page that has not prefetched yet.
let instance = 0;
const freshAssets = () => import(`${ASSETS}?page=${instance += 1}`);

/** Runs `prefetch` with fetch recorded and Save-Data set as given; returns the requests. */
const prefetchWith = async ({ saveData = false, respond = () => new Response('bytes') } = {}) => {
  const { prefetchPyodideBootAssets } = await freshAssets();
  const requests = [];
  const originals = { fetch: globalThis.fetch, location: globalThis.location, navigator: Object.getOwnPropertyDescriptor(globalThis, 'navigator') };
  globalThis.fetch = async (url, options) => { requests.push([String(url), options]); return respond(); };
  globalThis.location = { href: BASE };
  Object.defineProperty(globalThis, 'navigator', { value: { connection: { saveData } }, configurable: true });
  try {
    prefetchPyodideBootAssets();
    prefetchPyodideBootAssets();
    await new Promise((resolve) => setTimeout(resolve, 10));
  } finally {
    globalThis.fetch = originals.fetch;
    globalThis.location = originals.location;
    Object.defineProperty(globalThis, 'navigator', originals.navigator);
  }
  return requests;
};

test('the boot files are the Worker manifest URLs plus the Pyodide runtime files and distribution wheels', async () => {
  const { buildPyodideAssetManifest, resolvePyodideBootAssetUrls } = await freshAssets();
  const manifest = buildPyodideAssetManifest(BASE);
  const index = new URL(PYODIDE_INDEX_URL, BASE).toString();
  assert.deepEqual(resolvePyodideBootAssetUrls(BASE), [
    manifest.pyodideModuleUrl,
    ...PYODIDE_RUNTIME_FILES.map((file) => `${index}${file}`),
    ...PYODIDE_RUNTIME_WHEELS.map((path) => new URL(path, BASE).toString()),
    ...manifest.localWheelUrls,
    manifest.gbdrawWheelUrl
  ]);
  assert.equal(manifest.localWheelUrls.length, PYODIDE_LOCAL_WHEELS.length);
  assert.match(manifest.gbdrawWheelUrl, new RegExp(`/${GBDRAW_WHEEL_NAME.replaceAll('.', '\\.')}(\\?|$)`));
});

test('a page prefetches every boot file once, in order, at low priority', async () => {
  const { resolvePyodideBootAssetUrls } = await freshAssets();
  const requests = await prefetchWith();
  assert.deepEqual(requests, resolvePyodideBootAssetUrls(BASE).map((url) => [url, { priority: 'low' }]));
});

test('the prefetch keeps two downloads in flight', async () => {
  const { prefetchPyodideBootAssets } = await freshAssets();
  const started = [];
  const finish = [];
  const originals = { fetch: globalThis.fetch, location: globalThis.location };
  globalThis.fetch = (url) => new Promise((resolve) => { started.push(String(url)); finish.push(() => resolve(new Response('bytes'))); });
  globalThis.location = { href: BASE };
  try {
    prefetchPyodideBootAssets();
    assert.equal(started.length, 2);
    finish.shift()();
    await new Promise((resolve) => setTimeout(resolve, 10));
    assert.equal(started.length, 3);
  } finally {
    globalThis.fetch = originals.fetch;
    globalThis.location = originals.location;
  }
});

test('Save-Data prefetches nothing', async () => {
  assert.deepEqual(await prefetchWith({ saveData: true }), []);
});

test('a failed boot file download is dropped without an unhandled rejection', async () => {
  const rejections = [];
  const record = (reason) => rejections.push(reason);
  process.on('unhandledRejection', record);
  try {
    const requests = await prefetchWith({ respond: () => { throw new TypeError('Failed to fetch'); } });
    // Each lane stops at its failure: no retry and no request storm.
    assert.equal(requests.length, 2);
    await new Promise((resolve) => setTimeout(resolve, 10));
  } finally {
    process.off('unhandledRejection', record);
  }
  assert.deepEqual(rejections, []);
});
