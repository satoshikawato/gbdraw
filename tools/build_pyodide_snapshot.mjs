// Builds the diagram Worker's Pyodide snapshot (snapshot.bin, state.json and
// fs.bin) in Node, from the prepared browser wheel, for a test server
// (docs/internal/SELECTIVE_CI.md, "Pyodide snapshot"). The runtime it captures
// is the Worker's cold boot: the same package recipe (installRuntimePackages)
// on the dependency wheels config.js lists and the gbdraw wheel, then the
// Worker helpers. The build fails when an install fails or an installed
// distribution misses a requirement. Run tools/prepare_browser_wheel.py first.
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
  // With WebAssembly JSPI, Pyodide adds a function table entry at startup, so a
  // snapshot made without it does not restore in Chromium, which has it.
  if (typeof (/** @type {Record<string, unknown>} */ (WebAssembly)).Suspending !== 'function') {
    throw new Error(`Node ${process.version} has no WebAssembly JSPI (WebAssembly.Suspending); `
      + 'Chromium cannot restore a snapshot made without it. Use the Node version in .node-version.');
  }
  const outDir = resolve(process.argv[2]);
  const started = performance.now();
  const { GBDRAW_WHEEL_NAME, PYODIDE_INDEX_URL, PYODIDE_LOCAL_WHEELS } = await importWeb('js/config.js');
  const { PYTHON_HELPERS } = await importWeb('js/app/python-helpers.js');
  const {
    PYODIDE_SNAPSHOT_FILES, PYTHON_HELPER_FILE, capturePyodideSnapshot, installRuntimePackages, pyodideSnapshotKey
  } = await importWeb('js/workers/pyodide-snapshot.js');
  const indexURL = `${join(web, PYODIDE_INDEX_URL)}/`;
  const { loadPyodide } = await import(pathToFileURL(join(indexURL, 'pyodide.mjs')).href);
  const wheels = [...PYODIDE_LOCAL_WHEELS.map((/** @type {string} */ path) => join(web, path)), join(web, GBDRAW_WHEEL_NAME)]
    .map((source) => ({ source, bytes: new Uint8Array(readFileSync(source)) }));

  const pyodide = await loadPyodide({ indexURL, _makeSnapshot: true });
  await installRuntimePackages(pyodide, { wheels });
  await pyodide.runPythonAsync(PYTHON_HELPERS);
  // Every requirement of every installed distribution is installed: a
  // dependency the recipe misses fails the build, not a later Generate.
  const missing = JSON.parse(String(pyodide.runPython(`
import importlib.metadata as md, json
from micropip._vendored.packaging.src.packaging.requirements import Requirement
from micropip._vendored.packaging.src.packaging.utils import canonicalize_name
installed = {canonicalize_name(d.metadata['Name']) for d in md.distributions()}
json.dumps(sorted({f"{d.metadata['Name']} requires {r}" for d in md.distributions() for r in (d.requires or [])
  if (q := Requirement(r)) and (not q.marker or q.marker.evaluate({'extra': ''}))
  and canonicalize_name(q.name) not in installed}))
`, { globals: pyodide.toPy({}) })));
  if (missing.length) throw new Error(`The snapshot runtime misses dependencies: ${missing.join('; ')}`);
  pyodide.runPython('import gc; gc.collect()');
  const key = await pyodideSnapshotKey({
    localWheels: wheels.slice(0, -1).map(({ source }) => /** @type {const} */ ([source, readFileSync(source)])),
    gbdrawWheel: readFileSync(wheels[wheels.length - 1].source),
    helperFile: readFileSync(PYTHON_HELPER_FILE)
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
