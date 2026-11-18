import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

import { FRESH_GUARD_FILES, guardCommand } from '../../tools/ci-impact-fresh-guards.mjs';
import { makeTempDirSync } from '../web/helpers/temp-dir.cjs';

const REPOSITORY_ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const PYTHON_STRUCTURE_GUARD = 'tests/test_python_structure_ratchet.py';
const git = (...args) => spawnSync('git', args, { cwd: REPOSITORY_ROOT, encoding: 'utf8', maxBuffer: 1 << 26 });

test('PR / gate merges the pull request into the current dev after its trusted check and runs the fresh guards without a token', () => {
  const workflow = readFileSync(resolve(REPOSITORY_ROOT, '.github/workflows/test.yml'), 'utf8');
  const gate = workflow.match(/\n {2}pr-gate:\n[\s\S]*?(?=\n {2}[a-z0-9-]+:\n|$)/)[0];
  const [trusted, fresh] = gate.split('run: node .ci-trusted-base/tools/ci-impact.mjs gate\n');
  assert.ok(fresh, 'the fresh guards follow the trusted check');
  assert.doesNotMatch(trusted, /\.ci-fresh-merge/);
  assert.match(fresh, /ref: refs\/heads\/\$\{\{ github\.base_ref \}\}\n {10}path: \.ci-fresh-merge\n {10}fetch-depth: 0\n/);
  assert.match(fresh, /persist-credentials: false\n/);
  assert.match(fresh, /HEAD_SHA: \$\{\{ github\.event\.pull_request\.head\.sha \}\}\n/);
  // The merge is built by merge-tree and checked out as a commit (dev first, as in GitHub's merge commit):
  // `git merge` on this blob-less checkout asks GitHub for the merged blob it just wrote (OV-3782).
  assert.match(fresh, /git merge-tree --write-tree --messages "\$dev" "\$HEAD_SHA"/);
  assert.match(fresh, /commit-tree "\$tree" -p "\$dev" -p "\$HEAD_SHA"/);
  assert.doesNotMatch(fresh.replace(/^ *#.*\n/gm, ''), /\bgit\b[^\n]*\smerge\s/);
  // Only exit status 1 of merge-tree is a conflict; any other failure is reported as itself.
  assert.match(fresh, /\[ "\$status" -eq 1 \][\s\S]*title=Conflict with dev[\s\S]*\[ "\$status" -ne 0 \][\s\S]*title=Merge into dev failed/);
  assert.match(fresh, /working-directory: \.ci-fresh-merge\n {8}run: node tools\/ci-impact-fresh-guards\.mjs\n/);
  // Pull request code runs here: no token, no secret, and nothing after a failed trusted check.
  assert.doesNotMatch(fresh.replace(/^ *#.*\n/gm, ''), /token|secrets\.|\n {8}if:/i);
});

test('each fresh guard file exists, and a Python one takes no pytest fixture', () => {
  for (const file of FRESH_GUARD_FILES) {
    const source = readFileSync(resolve(REPOSITORY_ROOT, file), 'utf8');
    if (!file.endsWith('.py')) continue;
    const tests = source.match(/^def test_\w+\(.*$/gm);
    assert.ok(tests?.length, file);
    assert.deepEqual(tests.filter((line) => !/^def test_\w+\(\)/.test(line)), [], file);
  }
});

// #1907 passed CI on a base before #1876, which limited a function to 58 statements. A rerun of
// its gate alone released it after dev turned green on 2026-11-11, and its merge turned dev red.
const HEAD_1907 = '270b5023044b3eacabf4b5579883a4c1a4d8c963';
const MERGE_1907 = '54a6d084a0da39c999edd3f55ee80b57f9e11486';
const MERGE_1876 = '6fa719e3a6e64fc7c471db6dfba0d25ba53dfcaf';
const historyAvailable = [HEAD_1907, MERGE_1907, MERGE_1876]
  .every((commit) => git('cat-file', '-e', `${commit}^{commit}`).status === 0);

test('the Python structure guard passes #1907 as its run tested it and fails it merged into the dev that had #1876', {
  skip: historyAvailable ? false : 'needs the history of dev (a checkout with fetch-depth: 0)'
}, () => {
  assert.ok(FRESH_GUARD_FILES.includes(PYTHON_STRUCTURE_GUARD));
  assert.equal(git('merge-base', '--is-ancestor', MERGE_1876, HEAD_1907).status, 1);
  assert.equal(git('merge-base', '--is-ancestor', MERGE_1876, `${MERGE_1907}^1`).status, 0);
  const root = makeTempDirSync('fresh-guards-');
  const run = (commit) => {
    const directory = join(root, commit);
    for (const file of [PYTHON_STRUCTURE_GUARD, 'gbdraw/diagrams/circular/assemble.py']) {
      const shown = git('show', `${commit}:${file}`);
      assert.equal(shown.status, 0, shown.stderr);
      mkdirSync(dirname(join(directory, file)), { recursive: true });
      writeFileSync(join(directory, file), shown.stdout);
    }
    const [command, args] = guardCommand(PYTHON_STRUCTURE_GUARD);
    return spawnSync(command, args, { cwd: directory, encoding: 'utf8', maxBuffer: 1 << 24 });
  };
  const tested = run(HEAD_1907);
  assert.equal(tested.status, 0, tested.stdout + tested.stderr);
  const fresh = run(MERGE_1907);
  assert.equal(fresh.status, 1, fresh.stdout + fresh.stderr);
  assert.match(fresh.stderr, /_draw_resolved_circular_slot: 59 statements \(limit 58/);
});
