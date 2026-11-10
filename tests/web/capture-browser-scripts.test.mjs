// The browser scripts the documentation captures embed in Python
// (docs/capture/**/*.py, a module-level `_..._SCRIPT = r"""..."""`) must
// parse: one that does not throws a SyntaxError in every capture that runs it
// (OV-3217, a `const` declared twice).
import assert from 'node:assert/strict';
import { readFile, readdir } from 'node:fs/promises';
import { join, relative } from 'node:path';
import test from 'node:test';
import vm from 'node:vm';

const repoRoot = process.cwd();
const captureRoot = join(repoRoot, 'docs', 'capture');
const SCRIPT = /^(_[A-Z0-9_]*SCRIPT)\s*=\s*r"""([\s\S]*?)"""/gm;

const pythonFiles = async (dir) => (await Promise.all(
  (await readdir(dir, { withFileTypes: true })).map((entry) => {
    const path = join(dir, entry.name);
    if (entry.isDirectory()) return pythonFiles(path);
    return entry.name.endsWith('.py') ? [path] : [];
  })
)).flat();

/** @param {string} source */
const scriptsIn = (source) => [...source.matchAll(SCRIPT)].map(([, name, body]) => ({ name, body }));
/** @param {string} body */
const parseError = (body) => {
  try {
    new vm.Script(`(${body})`);
    return null;
  } catch (error) {
    return String(error);
  }
};

test('every browser script a documentation capture embeds parses', async () => {
  const scripts = [];
  for (const path of await pythonFiles(captureRoot)) {
    for (const script of scriptsIn(await readFile(path, 'utf8'))) scripts.push({ file: relative(repoRoot, path), ...script });
  }
  assert.ok(scripts.some(({ name }) => name === '_SVG_INSPECTION_SCRIPT'), 'the SVG inspection script is found');
  assert.deepEqual(scripts.flatMap(({ file, name, body }) => {
    const error = parseError(body);
    return error ? [`${file} ${name}: ${error}`] : [];
  }), []);
});

test('the check finds a script that declares a name twice', () => {
  const [script] = scriptsIn('_X_SCRIPT = r"""\nroot => {\n  const a = 1;\n  const a = 2;\n  return a;\n}\n"""\n');
  assert.match(String(parseError(script.body)), /SyntaxError/);
});
