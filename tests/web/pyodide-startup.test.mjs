import assert from 'node:assert/strict';

import { createPyodideManager } from '../../gbdraw/web/js/app/pyodide.js';
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
