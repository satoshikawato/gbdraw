// OV-110 (D-15-49): the Legend layout owner releases the composition offsets
// of a Result that nothing names, and keeps those of a Result a History step
// names, so its Undo still returns it to the offsets first seen.
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import test from 'node:test';

test('a released Result drops its offsets; a Result History names keeps them', () => {
  const run = spawnSync(process.execPath, [
    '--experimental-test-module-mocks', '--no-warnings',
    new URL('./helpers/composition-delta-release.mjs', import.meta.url).pathname
  ], { encoding: 'utf8' });
  assert.equal(run.status, 0, run.stderr);
  const out = JSON.parse(run.stdout.trim().split('\n').at(-1));
  const moved = { primary: [], legend: [5, 5], title: null };
  assert.deepEqual(out.moved, { R1: moved }, 'the move is recorded');
  assert.deepEqual(out.kept, { R1: moved }, '(b) a Result a step names keeps its move');
  assert.deepEqual(out.restoredTo, { primary: [], legend: [0, 0], title: null }, '(b) and its Undo returns to the first offsets');
  assert.deepEqual(out.restored, ['R1']);
  assert.deepEqual(out.released, {}, '(a) a dropped Result releases its offsets');
  assert.deepEqual(out.restoredAfterRelease, [], '(a) nothing is left to restore');
});
