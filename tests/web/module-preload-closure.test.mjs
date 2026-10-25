// The hosted bundle preloads the static import closure of js/app.js
// (`module_preload_paths` in tools/prepare_cloudflare_pages.py). This test
// derives the closure with V8's own module parser and compares the lists, so
// an import form the bundle's reader misses fails here.
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { test } from 'node:test';
import { fileURLToPath } from 'node:url';

const REPOSITORY_ROOT = fileURLToPath(new URL('../..', import.meta.url));

const V8_CLOSURE = `
const { readFileSync } = require('node:fs');
const { posix } = require('node:path');
const vm = require('node:vm');
const root = process.argv[1];
const seen = new Set();
const pending = ['js/app.js'];
while (pending.length) {
  const path = pending.shift();
  const module = new vm.SourceTextModule(readFileSync(posix.join(root, path), 'utf8'));
  for (const specifier of module.dependencySpecifiers) {
    const target = posix.normalize(posix.join(posix.dirname(path), specifier));
    if (!seen.has(target)) { seen.add(target); pending.push(target); }
  }
}
seen.delete('js/app.js');
process.stdout.write(JSON.stringify([...seen].sort()));
`;

test('the bundle preloads exactly the modules V8 finds statically imported from js/app.js', () => {
  const webRoot = `${REPOSITORY_ROOT}gbdraw/web`;
  const v8Closure = JSON.parse(execFileSync(process.execPath,
    ['--experimental-vm-modules', '--no-warnings', '-e', V8_CLOSURE, webRoot], { encoding: 'utf8' }));
  const bundleList = JSON.parse(execFileSync(process.env.PYTHON || 'python', ['-c', `
import json, runpy, sys
from pathlib import Path
module = runpy.run_path("tools/prepare_cloudflare_pages.py")
print(json.dumps(sorted(module["module_preload_paths"](Path("gbdraw/web")))))
`], { cwd: REPOSITORY_ROOT, encoding: 'utf8' }));
  assert.ok(v8Closure.length > 100, `closure of ${v8Closure.length} modules`);
  assert.deepEqual(bundleList, v8Closure);
});
