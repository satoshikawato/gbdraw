// Builds the diagram Worker's Pyodide snapshot (snapshot.bin, state.json and
// fs.bin) in Node, from the prepared browser wheel, for a test server
// (docs/internal/SELECTIVE_CI.md, "Pyodide snapshot"). The runtime it captures
// is the Worker's cold boot: the dependency wheels config.js lists, the gbdraw
// wheel, then the Worker helpers. Run tools/prepare_browser_wheel.py first.
//
// Usage: node tools/build_pyodide_snapshot.mjs <output directory>
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const root = resolve(fileURLToPath(new URL('..', import.meta.url)));
const web = join(root, 'gbdraw/web');
/** @param {string} path */
const importWeb = (path) => import(pathToFileURL(join(web, path)).href);

const main = async () => {
  if (process.argv.length !== 3) throw new Error('Usage: node tools/build_pyodide_snapshot.mjs <output directory>');
  const outDir = resolve(process.argv[2]);
  const started = performance.now();
  const { GBDRAW_WHEEL_NAME, PYODIDE_INDEX_URL, PYODIDE_LOCAL_WHEELS } = await importWeb('js/config.js');
  const { PYTHON_HELPERS } = await importWeb('js/app/python-helpers.js');
  const {
    PYODIDE_SNAPSHOT_FILES, PYTHON_HELPER_FILE, capturePyodideSnapshot, pyodideSnapshotKey
  } = await importWeb('js/workers/pyodide-snapshot.js');
  const indexURL = `${join(web, PYODIDE_INDEX_URL)}/`;
  const { loadPyodide } = await import(pathToFileURL(join(indexURL, 'pyodide.mjs')).href);
  const localWheels = PYODIDE_LOCAL_WHEELS.map((/** @type {string} */ path) => join(web, path));
  const gbdrawWheel = join(web, GBDRAW_WHEEL_NAME);

  const pyodide = await loadPyodide({ indexURL, _makeSnapshot: true });
  await pyodide.loadPackage('micropip');
  await pyodide.loadPackage(localWheels);
  await pyodide.loadPackage([gbdrawWheel]);
  await pyodide.runPythonAsync(PYTHON_HELPERS);
  pyodide.runPython('import gc; gc.collect()');
  const key = await pyodideSnapshotKey({
    localWheels: PYODIDE_LOCAL_WHEELS, gbdrawWheel: readFileSync(gbdrawWheel), helperFile: readFileSync(PYTHON_HELPER_FILE)
  });
  const { memory, state, files } = capturePyodideSnapshot(pyodide, key);
  mkdirSync(outDir, { recursive: true });
  writeFileSync(join(outDir, PYODIDE_SNAPSHOT_FILES.memory), memory);
  writeFileSync(join(outDir, PYODIDE_SNAPSHOT_FILES.state), JSON.stringify(state));
  writeFileSync(join(outDir, PYODIDE_SNAPSHOT_FILES.files), files);
  console.log(`Pyodide snapshot in ${outDir}: ${memory.length} + ${files.length} bytes, `
    + `${state.libs.length} side modules, ${state.unknownHiwire.length} unnamed JS references, `
    + `${((performance.now() - started) / 1000).toFixed(1)} s`);
};

main().then(
  () => process.exit(0),
  (error) => {
    console.error(error?.stack || error);
    process.exit(1);
  }
);
