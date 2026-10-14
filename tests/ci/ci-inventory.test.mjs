import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import {
  countWorkflowJobs, isoWeek, percentile, readData, renderReport, summarizeFunctional, summarizeRuns
} from '../../tools/ci-inventory.mjs';

const at = (minute) => new Date(Date.UTC(2026, 9, 14, 0, minute)).toISOString();
const job = (name, created, started, completed, extra = {}) => ({
  name, created_at: at(created), started_at: at(started), completed_at: at(completed),
  run_attempt: 1, conclusion: 'success', ...extra
});
const tests = (id, event, headBranch) => ({
  id, name: 'Tests', path: '.github/workflows/test.yml', event, head_branch: headBranch,
  run_attempt: 1, run_started_at: at(0), updated_at: at(60)
});

test('runner-minutes group by workflow, event, and job, and Tests time splits PR and dev runs', () => {
  const runs = [
    tests(1, 'pull_request', 'fix/a'),
    tests(2, 'push', 'dev'),
    { id: 3, name: 'PR #7', path: 'dynamic/github-code-scanning/codeql', event: 'dynamic',
      run_attempt: 1, run_started_at: at(0), updated_at: at(4) }
  ];
  const jobsByRun = new Map([
    [1, [
      job('Playwright functional (shard 1/8)', 0, 10, 40),
      job('Playwright functional (shard 2/8)', 0, 20, 30),
      job('Lint', 0, 0, 0, { conclusion: 'skipped' }),
      // An earlier attempt costs runner-minutes but not the latest wall time.
      job('Lint', 0, 1, 3, { run_attempt: 0 })
    ]],
    [2, [job('Core (Python 3.11)', 0, 2, 12)]]
  ]);
  const summary = summarizeRuns({ runs, jobsByRun });
  assert.equal(summary.runs, 3);
  assert.equal(summary.runDurationOnly, 1);
  assert.equal(summary.runnerMinutes, 56);
  assert.deepEqual(summary.byWorkflow, [
    { name: 'Tests [pull_request]', runs: 1, minutes: 42 },
    { name: 'Tests [push]', runs: 1, minutes: 10 },
    { name: 'codeql [dynamic]', runs: 1, minutes: 4 }
  ]);
  assert.deepEqual(summary.byJob.map((entry) => [entry.name, entry.minutes]), [
    ['Playwright functional', 40], ['Core', 10], ['codeql', 4], ['Lint', 2]
  ]);
  assert.deepEqual(summary.tests.pull_request, { runs: 1, wall: [40, 40], queue: [10, 20] });
  assert.deepEqual(summary.tests['dev push'], { runs: 1, wall: [12, 12], queue: [2, 2] });
});

const shardReport = (cases, durationMinutes) => ({
  config: { rootDir: '/home/runner/work/gbdraw/gbdraw/tests/web' },
  stats: { duration: durationMinutes * 60000 },
  suites: [{
    specs: cases.map(([file, minutes]) => ({
      file, line: 1, column: 1, title: file, tests: [{ results: [{ duration: minutes * 60000 }] }]
    }))
  }]
});

test('functional shards flag a rebalance for a slow shard or an unrecorded spec file', () => {
  const input = (minutes) => ({
    reports: [
      { shard: 2, report: shardReport([['b.spec.js', 2]], 1.5) },
      { shard: 1, report: shardReport([['a.spec.js', 8], ['c.spec.js', 2]], 6) }
    ],
    shardsFile: { minutes, shards: [['tests/web/a.spec.js', 'tests/web/c.spec.js'], ['tests/web/b.spec.js']] },
    casesPerFile: new Map([['tests/web/a.spec.js', 1], ['tests/web/b.spec.js', 1], ['tests/web/c.spec.js', 1]]),
    testDir: 'tests/web'
  });
  const functional = summarizeFunctional(input({ 'tests/web/a.spec.js': 8, 'tests/web/b.spec.js': 2 }));
  assert.deepEqual(functional.shards, [
    { shard: 1, testMinutes: 10, wallMinutes: 6 },
    { shard: 2, testMinutes: 2, wallMinutes: 1.5 }
  ]);
  assert.equal(functional.cases, 3);
  assert.equal(functional.testMinutes, 12);
  assert.equal(functional.slowestToMean, 1.67);
  assert.deepEqual(functional.unmeasured, ['tests/web/c.spec.js']);
  assert.equal(functional.rebalance, true);
  // The balance tool would move c.spec.js next to b.spec.js.
  assert.equal(functional.planChanges, true);
  assert.deepEqual(functional.plannedTestMinutes, [8, 4]);
  assert.deepEqual(functional.topSpecs.map((spec) => spec.file), [
    'tests/web/a.spec.js', 'tests/web/b.spec.js', 'tests/web/c.spec.js'
  ]);

  const even = summarizeFunctional({
    ...input({}),
    reports: [
      { shard: 1, report: shardReport([['a.spec.js', 5]], 3) },
      { shard: 2, report: shardReport([['b.spec.js', 4], ['c.spec.js', 1]], 3) }
    ],
    shardsFile: {
      minutes: { 'tests/web/a.spec.js': 5, 'tests/web/b.spec.js': 4, 'tests/web/c.spec.js': 1 },
      shards: [['tests/web/a.spec.js'], ['tests/web/b.spec.js', 'tests/web/c.spec.js']]
    }
  });
  assert.equal(even.rebalance, false);
  assert.equal(even.planChanges, false);
});

test('the report shows week-over-week changes and carries its numbers for next week', () => {
  const functional = {
    shards: [{ shard: 1, testMinutes: 10, wallMinutes: 6 }], cases: 3, testMinutes: 12,
    slowestToMean: 1.67, unmeasured: ['tests/web/c.spec.js'], rebalance: true, planChanges: true,
    plannedTestMinutes: [8], topSpecs: [{ file: 'tests/web/a.spec.js', cases: 1, minutes: 8 }]
  };
  const runs = summarizeRuns({ runs: [], jobsByRun: new Map() });
  const markdown = renderReport({
    week: '2026-W42', since: new Date(at(0)), until: new Date(at(60)), runs, functional,
    functionalRun: { id: 38490146362, head_sha: '0e66fd7bbfe17adc' }, testJobs: 24,
    previous: { week: '2026-W41', cases: 1000, functionalTestMinutes: 20, runnerMinutes: 0, testJobs: 22 },
    serverUrl: 'https://github.com', repository: 'satoshikawato/gbdraw'
  });
  assert.match(markdown, /\*\*Rebalance the functional shards\*\*: slowest shard 1\.67x the mean \(limit 1\.2x\); 1 spec file\(s\) without recorded minutes: `tests\/web\/c\.spec\.js`\./);
  assert.match(markdown, /\| Functional cases \| 3 \| 1000 \| -997 \|/);
  assert.match(markdown, /\| Jobs in `test\.yml` \| 24 \| 22 \| \+2 \|/);
  assert.match(markdown, /\[38490146362\]\(https:\/\/github\.com\/satoshikawato\/gbdraw\/actions\/runs\/38490146362\)/);
  assert.match(markdown, /\| 1 \| 10 \| 5 \| 6 \| 4 \|/);
  assert.deepEqual(readData(`older text\n${markdown}`), {
    week: '2026-W42', cases: 3, functionalTestMinutes: 12, runnerMinutes: 0, testJobs: 24
  });
  assert.equal(readData('no data here'), null);
});

test('helpers count workflow jobs, ISO weeks, and nearest-rank percentiles', () => {
  assert.equal(countWorkflowJobs('name: x\non: push\njobs:\n  a:\n    steps: []\n  b-c:\n    needs: a\n'), 2);
  assert.ok(countWorkflowJobs(readFileSync('.github/workflows/test.yml', 'utf8')) > 10);
  assert.equal(isoWeek(new Date('2026-10-14T00:00:00Z')), '2026-W42');
  assert.equal(isoWeek(new Date('2027-01-01T00:00:00Z')), '2026-W53');
  assert.equal(percentile([5, 1, 3, 2, 4], 50), 3);
  assert.equal(percentile([5, 1, 3, 2, 4], 90), 5);
  assert.equal(percentile([], 50), null);
});

test('the CI inventory workflow is weekly, read-only except issues, and posts its report', () => {
  const workflow = readFileSync('.github/workflows/ci-inventory.yml', 'utf8');
  assert.match(workflow, /schedule:\n\s+- cron: '[^']+'\n\s+workflow_dispatch:/);
  assert.match(workflow, /permissions:\n {2}contents: read\n {2}actions: read\n {2}issues: write\n/);
  assert.doesNotMatch(workflow, /pull_request|push:|contents: write|pull-requests: write/);
  assert.match(workflow, /node tools\/ci-inventory\.mjs --post >> "\$GITHUB_STEP_SUMMARY"/);
});
