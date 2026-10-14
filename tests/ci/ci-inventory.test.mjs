import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import {
  countWorkflowJobs, isoWeek, latestShardArtifacts, percentile, priorReports, readData, renderReport,
  summarizeFunctional, summarizeRuns
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
  // A run without a job list is one row marked as counted by run duration.
  assert.deepEqual(summary.byWorkflow, [
    { name: 'Tests [pull_request]', runs: 1, minutes: 42 },
    { name: 'Tests [push]', runs: 1, minutes: 10 },
    { name: 'codeql [dynamic] (run duration)', runs: 1, minutes: 4 }
  ]);
  assert.deepEqual(summary.byJob.map((entry) => [entry.name, entry.minutes]), [
    ['Playwright functional', 40], ['Core', 10], ['codeql (run duration)', 4], ['Lint', 2]
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

  // The flag and the printed ratio are the same rounded number: 1.2004 is 1.2, not above it.
  const edge = summarizeFunctional({
    ...input({ 'tests/web/a.spec.js': 6.002, 'tests/web/b.spec.js': 2.5, 'tests/web/c.spec.js': 1.5 }),
    reports: [
      { shard: 1, report: shardReport([['a.spec.js', 6.002]], 3) },
      { shard: 2, report: shardReport([['b.spec.js', 2.5], ['c.spec.js', 1.5]], 2) }
    ]
  });
  assert.equal(edge.slowestToMean, 1.2);
  assert.equal(edge.rebalance, false);
});

test('a re-run counts each shard once, from its latest attempt, and a missing report is left out', () => {
  const artifact = (shard, attempt, extra = {}) => ({
    id: shard * 10 + attempt, name: `playwright-functional-shard-${shard}-traces-38490146362-${attempt}`, ...extra
  });
  assert.deepEqual(latestShardArtifacts([
    artifact(1, 1), artifact(2, 1), artifact(2, 2), artifact(3, 2, { expired: true }),
    { id: 9, name: 'playwright-performance-traces-38490146362-1' }
  ]).map(({ shard, attempt }) => [shard, attempt]), [[1, 1], [2, 2]]);

  const functional = summarizeFunctional({
    reports: [{ shard: 2, report: shardReport([['b.spec.js', 4]], 2) }],
    shardsFile: { minutes: { 'tests/web/a.spec.js': 4, 'tests/web/b.spec.js': 4 },
      shards: [['tests/web/a.spec.js'], ['tests/web/b.spec.js']] },
    casesPerFile: new Map([['tests/web/a.spec.js', 1], ['tests/web/b.spec.js', 1]]),
    testDir: 'tests/web'
  });
  assert.deepEqual(functional.missingShards, [1]);
  assert.equal(functional.cases, 1);
  const markdown = renderReport({
    week: '2026-W42', since: new Date(at(0)), until: new Date(at(60)),
    runs: summarizeRuns({ runs: [], jobsByRun: new Map() }), functional,
    functionalRun: { id: 1, head_sha: '0e66fd7bbfe17adc' }, testJobs: 24, previous: null,
    serverUrl: 'https://github.com', repository: 'satoshikawato/gbdraw'
  });
  // Shard 2's planned minutes come from its own position in the plan, not the row index.
  assert.match(markdown, /\| 2 \| 4 \| 2 \| 2 \| 2 \|/);
  assert.match(markdown, /No functional report for shard 1; the totals leave it out\./);
});

test('the report shows week-over-week changes and carries its numbers for next week', () => {
  const functional = {
    shards: [{ shard: 1, testMinutes: 10, wallMinutes: 6 }], missingShards: [], cases: 3, testMinutes: 12,
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
  const partial = renderReport({
    week: '2026-W42', since: new Date(at(0)), until: new Date(at(60)), runs, functional,
    partial: 'runs created from 2026-10-10T00:00Z on are missing (HTTP 403).',
    functionalRun: { id: 1, head_sha: '0e66fd7bbfe17adc' }, testJobs: 24, previous: null,
    serverUrl: 'https://github.com', repository: 'satoshikawato/gbdraw'
  });
  assert.match(partial, /\*\*Partial report\*\*: runs created from 2026-10-10T00:00Z on are missing \(HTTP 403\)\./);
});

const block = (data) => `<!-- ci-inventory-data ${JSON.stringify(data)} -->`;
const numbers = { cases: 1000, functionalTestMinutes: 400, runnerMinutes: 80000, testJobs: 24 };
const comment = (id, login, data) => ({ id, user: { login }, body: `report\n${typeof data === 'string' ? data : block(data)}` });

test("only the workflow's own comments with well-formed numbers are earlier reports", () => {
  const bot = 'github-actions[bot]';
  const w41 = comment(1, bot, { week: '2026-W41', ...numbers });
  // An outsider's comment is neither this week's report to replace nor the previous one.
  const outsider = comment(2, 'someone', { week: '2026-W42', ...numbers });
  assert.deepEqual(priorReports([w41, outsider], '2026-W42'), {
    thisWeek: null, previous: { week: '2026-W41', ...numbers }
  });
  const w42 = comment(3, bot, { week: '2026-W42', ...numbers, cases: 1068 });
  assert.deepEqual(priorReports([w41, outsider, w42], '2026-W42'), {
    thisWeek: w42, previous: { week: '2026-W41', ...numbers }
  });
  // A block with any malformed field is dropped whole, and only known fields are kept.
  for (const bad of [
    { week: '2026-W41 @someone', ...numbers },
    { week: '2026-W41', ...numbers, cases: '1000 [link](https://example.com)' },
    { week: '2026-W41', ...numbers, runnerMinutes: null },
    { ...numbers }
  ]) assert.equal(readData(block(bad)), null, JSON.stringify(bad));
  assert.equal(readData('<!-- ci-inventory-data {not json} -->'), null);
  assert.deepEqual(readData(block({ week: '2026-W41', ...numbers, note: '<!-- x -->' })), { week: '2026-W41', ...numbers });
  assert.deepEqual(priorReports([comment(4, bot, { week: '2026-W41', ...numbers, testJobs: 'x' })], '2026-W42'),
    { thisWeek: null, previous: null });
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
  assert.match(workflow, /ref: dev\n\s+persist-credentials: false\n/);
});
