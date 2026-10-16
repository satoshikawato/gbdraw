// PROTOTYPE (proto/pyodide-snapshot, CI measurement only): build the diagram
// Worker's Pyodide snapshot in Node, without a browser. A Node-built snapshot
// restores in the Chromium Worker with a byte-identical first SVG.
// Usage: node tools/proto-pyodide-snapshot/make_snapshot_node.mjs <out dir>
import { mkdirSync, readdirSync, writeFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = resolve(fileURLToPath(new URL('../..', import.meta.url)));
const outDir = resolve(process.argv[2]);
const web = join(root, 'gbdraw/web');
const indexURL = join(web, 'vendor/pyodide/v0.29.0/full/');
const started = performance.now();
const { loadPyodide } = await import(join(indexURL, 'pyodide.mjs'));
const { GBDRAW_WHEEL_NAME } = await import(join(web, 'js/config.js'));
const { PYTHON_HELPERS } = await import(join(web, 'js/app/python-helpers.js'));
const proto = await import(join(web, 'js/workers/pyodide-snapshot-proto.js'));

const pyodide = await loadPyodide({ indexURL, _makeSnapshot: true });
await pyodide.loadPackage('micropip');
const wheels = join(web, 'vendor/pyodide-wheels');
await pyodide.loadPackage(readdirSync(wheels).filter((f) => f.endsWith('.whl')).map((f) => join(wheels, f)));
await pyodide.loadPackage([join(web, GBDRAW_WHEEL_NAME)]);
await pyodide.runPythonAsync(PYTHON_HELPERS);
const { state, fsBundle } = proto.captureSnapshotState(pyodide);
const { bytes, unknown } = proto.makeSnapshotWithPackages(pyodide);
state.unknownHiwire = unknown;
mkdirSync(outDir, { recursive: true });
writeFileSync(join(outDir, 'snapshot.bin'), bytes);
writeFileSync(join(outDir, 'state.json'), JSON.stringify(state));
writeFileSync(join(outDir, 'fs.bin'), fsBundle);
console.log(`Pyodide snapshot: ${bytes.length} + ${fsBundle.length} bytes, ${state.libs.length} side modules, `
  + `${unknown.length} unnamed JS references, ${((performance.now() - started) / 1000).toFixed(1)} s`);
