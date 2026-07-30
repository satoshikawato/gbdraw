import assert from 'node:assert/strict';

import { createPyodideManager } from '../../gbdraw/web/js/app/pyodide.js';
import { readFileBytes } from '../../gbdraw/web/js/services/file-content-cache.js';
import {
  DIAGRAM_ENGINE_STARTUP_MESSAGE
} from '../../gbdraw/web/js/services/runtime-capabilities.js';

globalThis.location = { href: 'https://example.test/gbdraw/web/' };

const startupError = new Error(
  'Missing packaged asset: gbdraw browser wheel. Rebuild and redeploy.'
);
globalThis.loadPyodide = async () => {
  throw startupError;
};

const ref = (value) => ({ value });
const loadingStatus = ref('');
const state = {
  pyodideReady: ref(false),
  loadingStatus,
  paletteDefinitions: ref({}),
  paletteNames: ref([]),
  selectedPalette: ref('default'),
  currentColors: ref({}),
  appliedPaletteName: ref('default'),
  appliedPaletteColors: ref({}),
  pendingPaletteName: ref(''),
  pendingPaletteColors: ref({}),
  normalizePaletteColors: (value) => value,
  normalizePaletteDefinitions: (value) => value
};
const diagnostics = [];
const originalConsoleError = console.error;
console.error = (...args) => diagnostics.push(args);
try {
  const manager = createPyodideManager({ state });
  assert.equal(await manager.initPyodide(), null);
} finally {
  console.error = originalConsoleError;
}

assert.equal(
  loadingStatus.value,
  `Startup Error: ${DIAGRAM_ENGINE_STARTUP_MESSAGE}`
);
assert.doesNotMatch(loadingStatus.value, /wheel|rebuild|redeploy/i);
assert.equal(diagnostics.length, 1);
assert.equal(diagnostics[0][0], 'Main-thread Pyodide startup failed.');
assert.equal(diagnostics[0][1], startupError);

const writes = [];
const runtime = {
  loadPackage: async () => {},
  pyimport: () => ({ install: async () => {} }),
  runPythonAsync: async () => {},
  runPython: () => JSON.stringify({ default: { CDS: '#cccccc' } }),
  FS: {
    writeFile: (path, bytes) => writes.push({ path, bytes })
  }
};
globalThis.fetch = async () => ({ ok: true, status: 200 });
globalThis.loadPyodide = async () => runtime;
const readyState = {
  ...state,
  pyodideReady: ref(false),
  loadingStatus: ref(''),
  paletteDefinitions: ref({}),
  paletteNames: ref([]),
  currentColors: ref({}),
  appliedPaletteColors: ref({})
};
const readyManager = createPyodideManager({ state: readyState });
assert.equal(await readyManager.initPyodide(), runtime);

let reads = 0;
const sourceBytes = new TextEncoder().encode('cached staging bytes');
const file = {
  async arrayBuffer() {
    reads += 1;
    return sourceBytes.slice().buffer;
  }
};
await Promise.all([
  readyManager.writeFileToFs(file, '/first.gb'),
  readyManager.writeFileToFs(file, '/second.gb')
]);
const cachedBytes = await readFileBytes(file);
assert.equal(reads, 1);
assert.deepEqual(writes.map((entry) => entry.path), ['/first.gb', '/second.gb']);
writes.forEach(({ bytes }) => {
  assert.deepEqual(bytes, cachedBytes);
  assert.notEqual(bytes, cachedBytes);
});

let paletteRuntimeStarts = 0;
globalThis.loadPyodide = async () => {
  paletteRuntimeStarts += 1;
  return runtime;
};
globalThis.fetch = async () => ({
  ok: true,
  status: 200,
  json: async () => ({
    palettes: {
      default: {
        CDS: '#54bcf8',
        repeat_region: '#d3d3d3'
      }
    }
  })
});
const paletteState = {
  ...state,
  pyodideReady: ref(false),
  loadingStatus: ref(''),
  paletteDefinitions: ref({}),
  paletteNames: ref([]),
  currentColors: ref({}),
  appliedPaletteColors: ref({})
};
const paletteManager = createPyodideManager({ state: paletteState });
await paletteManager.loadPaletteAsset();
assert.equal(paletteRuntimeStarts, 0);
assert.equal(paletteState.pyodideReady.value, false);
assert.equal(paletteState.currentColors.value.repeat_region, '#d3d3d3');
