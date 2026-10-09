// R14 "Typed boundaries" (gbdraw/web/CLAUDE.md): the noImplicitAny ratchet.
// `noImplicitAny` is not yet on in the guard config. Until it is,
// tests/web/types/tsconfig.no-implicit-any.json is the guard config
// (tests/web/types/tsconfig.json) plus `noImplicitAny`, and
// NO_IMPLICIT_ANY_BASELINE caps how many diagnostics one `tsc` run with it
// reports across gbdraw/web/js/. The guard asserts:
//   1. the ratchet config is exactly the guard config plus `noImplicitAny`;
//   2. every diagnostic of that run is in a module under gbdraw/web/js/;
//   3. the run reports no more diagnostics than NO_IMPLICIT_ANY_BASELINE.
//
// The cap is one number for the whole tree, so moving a module does not change
// it. A pull request may lower it to the new count; it need not. Raising it is
// an authority-only change: NO_IMPLICIT_ANY_BASELINE is a registered
// design-rule allowlist (tools/web-design-rule-guards.json, R14, kind `count`).
// Owner decision of 2026-10-07: this cap stands until typed boundaries phase 4
// (v0.15.0) replaces it with a count per module and then turns the option on.
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { dirname, join } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

const NO_IMPLICIT_ANY_BASELINE = 7540;

const REPOSITORY_ROOT = join(dirname(fileURLToPath(import.meta.url)), '..', '..');
const RATCHET_CONFIG = 'tests/web/types/tsconfig.no-implicit-any.json';
const EXPECTED_RATCHET_CONFIG = { extends: './tsconfig.json', compilerOptions: { noImplicitAny: true } };
const DIAGNOSTIC = /^gbdraw\/web\/js\/(.+?)\((\d+),(\d+)\): error (TS\d+): /;

const tsc = () => {
  const require = createRequire(import.meta.url);
  let packagePath;
  try {
    packagePath = require.resolve('typescript/package.json');
  } catch {
    assert.fail('typescript is not installed: run `npm ci` (R14 runs the pinned devDependency)');
  }
  const pinned = JSON.parse(readFileSync(join(REPOSITORY_ROOT, 'package.json'), 'utf8')).devDependencies.typescript;
  const installed = JSON.parse(readFileSync(packagePath, 'utf8')).version;
  assert.equal(installed, pinned, `typescript ${installed} is installed but package.json pins ${pinned}: run \`npm ci\``);
  const run = spawnSync(
    process.execPath,
    [join(dirname(packagePath), 'bin', 'tsc'), '-p', RATCHET_CONFIG, '--pretty', 'false'],
    { cwd: REPOSITORY_ROOT, encoding: 'utf8', maxBuffer: 1 << 28 }
  );
  assert.ifError(run.error);
  return `${run.stdout}\n${run.stderr}`;
};

test('R14: the noImplicitAny config is the guard config plus noImplicitAny', () => {
  const config = JSON.parse(readFileSync(join(REPOSITORY_ROOT, RATCHET_CONFIG), 'utf8'));
  assert.deepEqual(
    config,
    EXPECTED_RATCHET_CONFIG,
    `${RATCHET_CONFIG} changes more than noImplicitAny; a config change is a rule change (authority-only)`
  );
});

test('R14: noImplicitAny diagnostics stay at or below NO_IMPLICIT_ANY_BASELINE', (t) => {
  assert.ok(Number.isInteger(NO_IMPLICIT_ANY_BASELINE) && NO_IMPLICIT_ANY_BASELINE >= 0);
  const perModule = new Map();
  const unparsed = [];
  let total = 0;
  let inDiagnostic = false;
  tsc().split(/\r?\n/).forEach((line) => {
    if (!line.trim()) return;
    const match = line.match(DIAGNOSTIC);
    if (match) {
      total += 1;
      perModule.set(match[1], (perModule.get(match[1]) || 0) + 1);
      inDiagnostic = true;
    } else if (/^\s/.test(line) && inDiagnostic) {
      // A continuation line of the diagnostic above.
    } else unparsed.push(line);
  });
  assert.deepEqual(unparsed, [], `tsc -p ${RATCHET_CONFIG} printed output outside gbdraw/web/js/`);

  if (total < NO_IMPLICIT_ANY_BASELINE) {
    t.diagnostic(`noImplicitAny: ${total} diagnostics, below the cap ${NO_IMPLICIT_ANY_BASELINE}; `
      + `NO_IMPLICIT_ANY_BASELINE may be lowered to ${total}`);
  }
  const counts = [...perModule].sort(([a], [b]) => a.localeCompare(b)).map(([path, count]) => `${path} ${count}`);
  assert.ok(
    total <= NO_IMPLICIT_ANY_BASELINE,
    `${total} noImplicitAny diagnostics, above NO_IMPLICIT_ANY_BASELINE ${NO_IMPLICIT_ANY_BASELINE}. `
      + 'Declare the types of the new parameters and variables (for example `@param {string} name`), '
      + 'or type the collection a callback iterates; raising the cap is an authority-only change. '
      + `Run \`node node_modules/typescript/bin/tsc -p ${RATCHET_CONFIG} --pretty false\` and compare with dev. `
      + `Diagnostics per module:\n    ${counts.join('\n    ')}`
  );
});
