import assert from 'node:assert/strict';
import test from 'node:test';

import {
  measureMinutes, measureSecondsPerCase, parseMeasurements, planShards
} from '../../tools/balance-functional-shards.mjs';

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

test('seconds per case are the mean of the median case durations of each file', () => {
  const seconds = measureSecondsPerCase([
    ciReport([['a.playwright.spec.js', 'one', [60000]], ['a.playwright.spec.js', 'two', [30000]]]),
    ciReport([['a.playwright.spec.js', 'one', [120000, 60000]], ['a.playwright.spec.js', 'two', [30000]]])
  ], 'tests/web');
  // one: median of 60 s and 180 s; two: 30 s. Two cases, not four reported ones.
  assert.deepEqual(Object.fromEntries(seconds), { 'tests/web/a.playwright.spec.js': 75 });
});

test('the plan weighs files by their current cases and estimates unmeasured files from the mean case', () => {
  const plan = planShards({
    // a grew from the 2 cases it was measured with to 4; d has no measurement.
    casesPerFile: new Map([['a', 4], ['b', 1], ['c', 1], ['d', 2]]),
    secondsPerCase: new Map([['a', 90], ['b', 120], ['c', 60], ['stale', 600]]),
    shardCount: 2
  });
  // a: 4 * 90 s = 6 min. The 6 measured cases take 540 s, 90 s each, so d is 2 * 90 s = 3 min.
  // A measurement of a file that no longer exists plays no part.
  assert.deepEqual(plan.minutes, { a: 6, b: 2, c: 1, d: 3 });
  assert.deepEqual(plan.shards, [['a'], ['b', 'c', 'd']]);
  assert.deepEqual(plan.shardMinutes, [6, 6]);
  assert.deepEqual(plan.unmeasured, ['d']);
});

test('the measurements name the shard count and the seconds per case of each spec file', () => {
  const parsed = parseMeasurements(JSON.stringify({ shardCount: 8, secondsPerCase: { a: 1.5 } }));
  assert.equal(parsed.shardCount, 8);
  assert.deepEqual(Object.fromEntries(parsed.secondsPerCase), { a: 1.5 });
  for (const data of [
    { secondsPerCase: { a: 1 } },
    { shardCount: 0, secondsPerCase: { a: 1 } },
    { shardCount: 8, secondsPerCase: {} },
    { shardCount: 8, secondsPerCase: { a: -1 } },
    { shardCount: 8, secondsPerCase: { a: '1' } }
  ]) assert.throws(() => parseMeasurements(JSON.stringify(data)), JSON.stringify(data));
});
