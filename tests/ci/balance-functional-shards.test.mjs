import assert from 'node:assert/strict';
import test from 'node:test';

import { measureMinutes, planShards } from '../../tools/balance-functional-shards.mjs';

// A functional report as a CI shard writes it: rootDir is the runner's test directory,
// not the directory the tool runs in.
const ciReport = (cases) => ({
  config: { rootDir: '/home/runner/work/gbdraw/gbdraw/tests/web' },
  suites: [{
    file: 'a.playwright.spec.js',
    specs: cases.map(([file, title, durations]) => ({
      file, line: 1, column: 1, title, tests: [{ results: durations.map((duration) => ({ duration })) }]
    }))
  }]
});

test('measured minutes keep spec files of reports written on another machine', () => {
  const minutes = measureMinutes([
    ciReport([['a.playwright.spec.js', 'one', [60000]], ['contracts/b.playwright.spec.js', 'two', [30000]]]),
    ciReport([['a.playwright.spec.js', 'one', [120000, 60000]]])
  ], 'tests/web');
  // A case's minutes are the median over the runs that report it, retries included.
  assert.deepEqual(Object.fromEntries(minutes), {
    'tests/web/a.playwright.spec.js': 2,
    'tests/web/contracts/b.playwright.spec.js': 0.5
  });
});

test('the plan assigns files longest first and estimates unmeasured files from their cases', () => {
  const plan = planShards({
    casesPerFile: new Map([['a', 2], ['b', 1], ['c', 1], ['d', 4]]),
    fileMinutes: new Map([['a', 6], ['b', 2], ['c', 1]]),
    shardCount: 2
  });
  // 9 measured minutes over 4 measured cases estimate d at 4 * 2.25 = 9 minutes.
  assert.deepEqual(plan.minutes, { a: 6, b: 2, c: 1, d: 9 });
  assert.deepEqual(plan.shards, [['d'], ['a', 'b', 'c']]);
  assert.deepEqual(plan.shardMinutes, [9, 9]);
  assert.deepEqual(plan.unmeasured, ['d']);
});
