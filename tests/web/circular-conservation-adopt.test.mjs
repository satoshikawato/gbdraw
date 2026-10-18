// A Generate reconciles the Circular comparison rings only when its adopt
// installed Circular ring files (a Circular LOSAT-cache result). A Linear
// Generate leaves the Circular drawing as Reset Settings left it (OV-1167)
// (A2 review-I2, Medium 1).
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import test from 'node:test';

test('a Linear Generate leaves the Circular rings that Reset turned off', () => {
  const run = spawnSync(process.execPath, [
    '--experimental-test-module-mocks', '--no-warnings',
    new URL('./helpers/rings-after-adopt.mjs', import.meta.url).pathname
  ], { encoding: 'utf8' });
  assert.equal(run.status, 0, run.stderr);
  const { before, after } = JSON.parse(run.stdout.trim().split('\n').at(-1));
  assert.deepEqual(before, { enabled: false, series: 0 });
  assert.deepEqual(after, before, 'the Circular rings stay off with no rows');
});
