// Adds the bytecode of every module to a browser wheel, compiled by the
// vendored Pyodide that loads it (tools/browser_wheel_bytecode.py): a .pyc
// holds its interpreter's magic number, and the Python that builds the wheel
// (3.11 in CI) is not Pyodide's. tools/prepare_browser_wheel.py runs it on the
// wheel it builds, and names the Pyodide directory from gbdraw/_web_assets.py
// (not config.js, which that script writes only afterwards).
//
// Usage: node tools/compile_browser_wheel_bytecode.mjs <wheel> <pyodide-dir>   (wheel rewritten in place)
import { readFileSync, writeFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const root = resolve(fileURLToPath(new URL('..', import.meta.url)));

const main = async () => {
  if (process.argv.length !== 4) {
    throw new Error('Usage: node tools/compile_browser_wheel_bytecode.mjs <wheel> <pyodide-dir>');
  }
  const wheel = resolve(process.argv[2]);
  const indexURL = `${resolve(process.argv[3])}/`;
  const { loadPyodide } = await import(pathToFileURL(join(indexURL, 'pyodide.mjs')).href);
  const pyodide = await loadPyodide({ indexURL, stdout: () => {}, stderr: (line) => console.error(line) });
  pyodide.FS.writeFile('/tmp/source.whl', readFileSync(wheel));
  const count = pyodide.runPython(`${readFileSync(join(root, 'tools/browser_wheel_bytecode.py'), 'utf8')}
add_bytecode('/tmp/source.whl', '/tmp/target.whl')`);
  writeFileSync(wheel, pyodide.FS.readFile('/tmp/target.whl'));
  console.log(`Added ${count} bytecode files (${pyodide.runPython('import sys; sys.implementation.cache_tag')}) to ${wheel}`);
};

main().then(
  () => process.exit(0),
  (error) => {
    console.error(error?.stack || error);
    process.exit(1);
  }
);
