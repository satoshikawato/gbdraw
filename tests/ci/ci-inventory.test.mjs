import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import {
  PR_TARGETS, PR_WINDOW, SUPERSEDED_GAP_SECONDS, countWorkflowJobs, createClient, listRuns, isoWeek, latestShardArtifacts, percentile, priorReports, readData,
  renderReport, summarizeFunctional, summarizePullRequestCi, summarizePullRequests, summarizeRuns
} from '../../tools/ci-inventory.mjs';

const noPullRequestCi = summarizePullRequestCi({ runs: [], earlierRuns: [], jobsByRun: new Map(), pulls: [] });
const noShares = { onePassPercent: null, runsPerPullRequest: null, supersededPercent: null, failedPercent: null };
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

// Recorded measurements from minutes per spec file; every fixture file has one case.
const measured = (minutes, shardCount = 2) => ({
  shardCount,
  secondsPerCase: new Map(Object.entries(minutes).map(([file, value]) => [`tests/web/${file}`, value * 60]))
});

test('functional shards ask for fresh measurements for a slow shard or an unmeasured spec file', () => {
  const input = (minutes) => ({
    reports: [
      { shard: 2, report: shardReport([['b.spec.js', 2]], 1.5) },
      { shard: 1, report: shardReport([['a.spec.js', 8], ['c.spec.js', 2]], 6) }
    ],
    measurements: measured(minutes),
    casesPerFile: new Map([['tests/web/a.spec.js', 1], ['tests/web/b.spec.js', 1], ['tests/web/c.spec.js', 1]]),
    testDir: 'tests/web'
  });
  const functional = summarizeFunctional(input({ 'a.spec.js': 8, 'b.spec.js': 2 }));
  assert.deepEqual(functional.shards, [
    { shard: 1, testMinutes: 10, wallMinutes: 6 },
    { shard: 2, testMinutes: 2, wallMinutes: 1.5 }
  ]);
  assert.equal(functional.cases, 3);
  assert.equal(functional.testMinutes, 12);
  assert.equal(functional.slowestToMean, 1.67);
  assert.deepEqual(functional.unmeasured, ['tests/web/c.spec.js']);
  assert.equal(functional.refresh, true);
  // c.spec.js is estimated at the mean case (5 minutes), which already puts it next to b.spec.js.
  assert.equal(functional.planChanges, false);
  assert.deepEqual(functional.plannedTestMinutes, [8, 4]);
  assert.deepEqual(functional.topSpecs.map((spec) => spec.file), [
    'tests/web/a.spec.js', 'tests/web/b.spec.js', 'tests/web/c.spec.js'
  ]);

  const evenReports = [
    { shard: 1, report: shardReport([['a.spec.js', 5]], 3) },
    { shard: 2, report: shardReport([['b.spec.js', 4], ['c.spec.js', 1]], 3) }
  ];
  const even = summarizeFunctional({
    ...input({ 'a.spec.js': 5, 'b.spec.js': 4, 'c.spec.js': 1 }), reports: evenReports
  });
  assert.equal(even.refresh, false);
  assert.equal(even.planChanges, false);
  // Stale measurements (a 3, b 3, c 2 minutes) plan [a, c] and [b]; this run's plan [a] and [b, c].
  const stale = summarizeFunctional({
    ...input({ 'a.spec.js': 3, 'b.spec.js': 3, 'c.spec.js': 2 }), reports: evenReports
  });
  assert.equal(stale.planChanges, true);
  assert.deepEqual(stale.plannedTestMinutes, [5, 5]);

  // The flag and the printed ratio are the same rounded number: 1.2004 is 1.2, not above it.
  const edge = summarizeFunctional({
    ...input({ 'a.spec.js': 6.002, 'b.spec.js': 2.5, 'c.spec.js': 1.5 }),
    reports: [
      { shard: 1, report: shardReport([['a.spec.js', 6.002]], 3) },
      { shard: 2, report: shardReport([['b.spec.js', 2.5], ['c.spec.js', 1.5]], 2) }
    ]
  });
  assert.equal(edge.slowestToMean, 1.2);
  assert.equal(edge.refresh, false);
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
    measurements: measured({ 'a.spec.js': 4, 'b.spec.js': 4 }),
    casesPerFile: new Map([['tests/web/a.spec.js', 1], ['tests/web/b.spec.js', 1]]),
    testDir: 'tests/web'
  });
  assert.deepEqual(functional.missingShards, [1]);
  assert.equal(functional.cases, 1);
  const markdown = renderReport({
    week: '2026-W42', since: new Date(at(0)), until: new Date(at(60)),
    runs: summarizeRuns({ runs: [], jobsByRun: new Map() }), pullRequests: summarizePullRequests([]), functional,
    pullRequestCi: noPullRequestCi, functionalRun: { id: 1, head_sha: '0e66fd7bbfe17adc' }, testJobs: 24, previous: null,
    serverUrl: 'https://github.com', repository: 'satoshikawato/gbdraw'
  });
  // Shard 2's planned minutes come from its own position in the plan, not the row index.
  assert.match(markdown, /\| 2 \| 4 \| 2 \| 2 \| 2 \|/);
  assert.match(markdown, /No functional report for shard 1; the totals leave it out\./);
});

test('the report shows week-over-week changes and carries its numbers for next week', () => {
  const functional = {
    shards: [{ shard: 1, testMinutes: 10, wallMinutes: 6 }], missingShards: [], cases: 3, testMinutes: 12,
    slowestToMean: 1.67, unmeasured: ['tests/web/c.spec.js'], refresh: true, planChanges: true,
    plannedTestMinutes: [8], topSpecs: [{ file: 'tests/web/a.spec.js', cases: 1, minutes: 8 }]
  };
  const runs = summarizeRuns({ runs: [], jobsByRun: new Map() });
  const pullRequests = summarizePullRequests([]);
  const markdown = renderReport({
    week: '2026-W42', since: new Date(at(0)), until: new Date(at(60)), runs, pullRequests, functional,
    pullRequestCi: noPullRequestCi,
    functionalRun: { id: 38490146362, head_sha: '0e66fd7bbfe17adc' }, testJobs: 24,
    previous: { week: '2026-W41', cases: 1000, functionalTestMinutes: 20, runnerMinutes: 0, testJobs: 22 },
    serverUrl: 'https://github.com', repository: 'satoshikawato/gbdraw'
  });
  assert.match(markdown, /\*\*Refresh the functional shard measurements\*\*: slowest shard 1\.67x the mean \(limit 1\.2x\); 1 spec file\(s\) without measurements: `tests\/web\/c\.spec\.js`\./);
  assert.match(markdown, /\| Functional cases \| 3 \| 1000 \| -997 \|/);
  assert.match(markdown, /\| Jobs in `test\.yml` \| 24 \| 22 \| \+2 \|/);
  assert.match(markdown, /\[38490146362\]\(https:\/\/github\.com\/satoshikawato\/gbdraw\/actions\/runs\/38490146362\)/);
  assert.match(markdown, /\| 1 \| 10 \| 5 \| 6 \| 4 \|/);
  assert.deepEqual(readData(`older text\n${markdown}`), {
    week: '2026-W42', cases: 3, functionalTestMinutes: 12, runnerMinutes: 0, testJobs: 24, ...noShares
  });
  assert.equal(readData('no data here'), null);
  const partial = renderReport({
    week: '2026-W42', since: new Date(at(0)), until: new Date(at(60)), runs, pullRequests, functional,
    pullRequestCi: noPullRequestCi,
    partial: 'runs created from 2026-10-10T00:00Z on are missing (HTTP 403).',
    functionalRun: { id: 1, head_sha: '0e66fd7bbfe17adc' }, testJobs: 24, previous: null,
    serverUrl: 'https://github.com', repository: 'satoshikawato/gbdraw'
  });
  assert.match(partial, /\*\*Partial report\*\*: runs created from 2026-10-10T00:00Z on are missing \(HTTP 403\)\./);
});

test('pull request feedback time splits by functional Playwright, takes the median, and warns over target', () => {
  // A green head run: planning at minute 1, its slowest job ends at `end`, and a skipped job
  // stamped later than that does not count.
  const headRun = (end, functional) => [
    job('Plan', 0, 1, 2),
    job(functional ? 'Playwright functional (shard 2/8)' : 'Web contracts', 2, 3, end),
    job('Recipes', 0, 0, 0, { conclusion: 'skipped', started_at: at(end + 30), completed_at: at(end + 30) })
  ];
  const summary = summarizePullRequests([
    headRun(18, false), headRun(26, false), headRun(22, false),
    headRun(41, true), headRun(31, true), headRun(51, true), headRun(35, true),
    // A pull request without a green head run, and one whose jobs were all skipped.
    null, [job('Lint', 0, 0, 0, { conclusion: 'skipped' })]
  ]);
  assert.equal(summary.pullRequests, 9);
  assert.equal(summary.measured, 7);
  assert.deepEqual(summary.kinds.map(({ kind, runs, median, minutes, over }) => [kind, runs, median, minutes, over]), [
    ['without functional Playwright', 3, 21, 20, true],
    // Nearest rank: the second of 30, 34, 40, 50.
    ['with functional Playwright', 4, 34, 35, false]
  ]);

  const report = (pullRequests) => renderReport({
    week: '2026-W42', since: new Date(at(0)), until: new Date(at(60)),
    runs: summarizeRuns({ runs: [], jobsByRun: new Map() }), pullRequests, pullRequestCi: noPullRequestCi,
    functional: summarizeFunctional({
      reports: [], measurements: measured({}, 0), casesPerFile: new Map(), testDir: 'tests/web'
    }),
    functionalRun: { id: 1, head_sha: '0e66fd7bbfe17adc' }, testJobs: 24, previous: null,
    serverUrl: 'https://github.com', repository: 'satoshikawato/gbdraw'
  });
  const over = report(summary);
  assert.match(over, /\*\*Pull request feedback over target\*\*: median 21 minutes without functional Playwright \(target 20\)\./);
  assert.match(over, /last 9 pull requests merged into `dev` \(7 had one\)/);
  assert.match(over, /\| without functional Playwright \| 3 \| 21 \| 20 \| \*\*yes\*\* \|/);
  assert.match(over, /\| with functional Playwright \| 4 \| 34 \| 35 \| no \|/);

  // The warning and the printed median are the same rounded number: 20.03 is 20, not above it.
  const edge = summarizePullRequests([[{
    name: 'Web contracts', conclusion: 'success', started_at: '2026-10-14T00:00:00Z', completed_at: '2026-10-14T00:20:02Z'
  }]]);
  assert.equal(edge.kinds[0].over, false);
  assert.match(report(edge), /Pull request feedback is within target \(median 20 of 20 minutes without functional Playwright; median - of 35 minutes with functional Playwright\)\./);
});

test('the targets, window, and cancelled-run gap the inventory uses are the ones SELECTIVE_CI.md states', () => {
  const doc = readFileSync('docs/internal/SELECTIVE_CI.md', 'utf8').replace(/\s+/g, ' ');
  const [without, withFunctional] = PR_TARGETS;
  assert.equal(without.functional, false);
  assert.equal(withFunctional.functional, true);
  assert.ok(doc.includes(`PR feedback targets ${without.minutes} minutes without functional Playwright`
    + ` and ${withFunctional.minutes} minutes with it`));
  assert.ok(doc.includes(`last ${PR_WINDOW} pull requests merged into \`dev\``));
  assert.ok(doc.includes(`created at most ${SUPERSEDED_GAP_SECONDS} seconds after it ended`));
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
    thisWeek: null, previous: { week: '2026-W41', ...numbers, ...noShares }
  });
  const w42 = comment(3, bot, { week: '2026-W42', ...numbers, cases: 1068 });
  assert.deepEqual(priorReports([w41, outsider, w42], '2026-W42'), {
    thisWeek: w42, previous: { week: '2026-W41', ...numbers, ...noShares }
  });
  // A block with any malformed field is dropped whole, and only known fields are kept.
  for (const bad of [
    { week: '2026-W41 @someone', ...numbers },
    { week: '2026-W41', ...numbers, cases: '1000 [link](https://example.com)' },
    { week: '2026-W41', ...numbers, runnerMinutes: null },
    { ...numbers }
  ]) assert.equal(readData(block(bad)), null, JSON.stringify(bad));
  assert.equal(readData('<!-- ci-inventory-data {not json} -->'), null);
  assert.deepEqual(readData(block({ week: '2026-W41', ...numbers, note: '<!-- x -->' })), {
    week: '2026-W41', ...numbers, ...noShares
  });
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

// Pull request `Tests` runs on 2026-10-14; `at(minute)` stamps them.
const prRun = (id, branch, created, ended, conclusion, extra = {}) => ({
  id, path: '.github/workflows/test.yml', event: 'pull_request', head_branch: branch,
  head_repository: { full_name: 'satoshikawato/gbdraw' }, created_at: at(created), updated_at: at(ended),
  run_attempt: 1, conclusion, ...extra
});
const pull = (number, branch, created, merged) => ({
  number, created_at: at(created), merged_at: at(merged),
  head: { ref: branch, repo: { full_name: 'satoshikawato/gbdraw' } }
});

test('pull request CI counts single-run passes, runs per pull request, and splits cancelled runs', () => {
  const runs = [
    // fix/a: one clean run.
    prRun(1, 'fix/a', 0, 30, 'success'),
    // fix/b: a newer push 30 seconds after the end stops the first run.
    prRun(2, 'fix/b', 0, 10, 'cancelled'),
    { ...prRun(3, 'fix/b', 10, 40, 'success'), created_at: '2026-10-14T00:10:30.000Z' },
    // fix/c: fail-fast cancels the first run; the next push comes 10 minutes later.
    prRun(4, 'fix/c', 0, 20, 'cancelled'),
    prRun(5, 'fix/c', 30, 60, 'success'),
    // fix/d: a failure, then a pass.
    prRun(6, 'fix/d', 0, 15, 'failure'),
    prRun(7, 'fix/d', 20, 50, 'success'),
    // fix/g was opened before the window; its earlier run comes from `earlierRuns`.
    prRun(8, 'fix/g', 5, 35, 'success'),
    // fix/h passed only after a re-run of its one run.
    prRun(9, 'fix/h', 0, 40, 'success', { run_attempt: 2 }),
    // Cancelled with no later run: a failed job other than `PR / gate` makes it failed.
    prRun(10, 'fix/e', 0, 12, 'cancelled'),
    // Without a failed job, or without a job list, it counts as stopped by a newer push.
    prRun(11, 'fix/f', 0, 12, 'cancelled'),
    prRun(12, 'fix/k', 0, 12, 'cancelled'),
    // A fork's branch of the same name is another branch, and a dev push run is not a PR run.
    prRun(13, 'fix/a', 1, 31, 'failure', { head_repository: { full_name: 'someone/gbdraw' } }),
    { ...prRun(14, 'dev', 0, 30, 'failure'), event: 'push' },
    { ...prRun(15, 'fix/a', 0, 30, 'failure'), path: '.github/workflows/web-policy.yml' }
  ];
  const earlierRuns = [prRun(20, 'fix/g', -600, -580, 'failure')];
  const jobsByRun = new Map([
    [10, [job('Core PR (Python 3.11)', 0, 1, 10, { conclusion: 'failure' }), job('Fail fast (Core PR)', 10, 10, 11)]],
    [11, [job('Core PR (Python 3.11)', 0, 1, 10, { conclusion: 'cancelled' }),
      job('PR / gate', 0, 11, 12, { conclusion: 'failure' }),
      // A failure of an earlier attempt does not count.
      job('Lint', 0, 1, 2, { conclusion: 'failure', run_attempt: 0 })]]
  ]);
  const pulls = [
    pull(101, 'fix/a', -1, 35), pull(102, 'fix/b', -1, 45), pull(103, 'fix/c', -1, 65), pull(104, 'fix/d', -1, 55),
    pull(107, 'fix/g', -700, 40), pull(108, 'fix/h', -1, 45),
    // Merged without a run, and a run long before this pull request opened.
    pull(109, 'fix/i', -1, 5), pull(110, 'fix/e', 600, 700)
  ];
  const ci = summarizePullRequestCi({ runs, earlierRuns, jobsByRun, pulls });
  assert.deepEqual(ci, {
    mergedPullRequests: 8,
    withRuns: 6,
    singlePass: 1,
    pullRequestRuns: 6 + 5 + 2,
    runs: { passed: 6, failed: 4, superseded: 3, other: 0 },
    failFast: 2,
    onePassPercent: 16.7,
    runsPerPullRequest: 1.7,
    supersededPercent: 23.1,
    failedPercent: 30.8
  });
  const empty = summarizePullRequestCi({ runs: [], earlierRuns: [], jobsByRun: new Map(), pulls: [] });
  assert.deepEqual([empty.onePassPercent, empty.runsPerPullRequest, empty.supersededPercent, empty.failedPercent],
    [null, null, null, null]);
});

test('the report shows pull request CI next to the previous week and carries it in the data block', () => {
  const functional = summarizeFunctional({
    reports: [], measurements: measured({}, 0), casesPerFile: new Map(), testDir: 'tests/web'
  });
  const pullRequestCi = summarizePullRequestCi({
    runs: [prRun(1, 'fix/a', 0, 30, 'success'), prRun(2, 'fix/b', 0, 15, 'failure')], earlierRuns: [],
    jobsByRun: new Map(), pulls: [pull(101, 'fix/a', -1, 35)]
  });
  const report = (previous) => renderReport({
    week: '2026-W42', since: new Date(at(0)), until: new Date(at(60)),
    runs: summarizeRuns({ runs: [], jobsByRun: new Map() }), pullRequests: summarizePullRequests([]), pullRequestCi,
    functional, functionalRun: { id: 1, head_sha: '0e66fd7bbfe17adc' }, testJobs: 24, previous,
    serverUrl: 'https://github.com', repository: 'satoshikawato/gbdraw'
  });
  const markdown = report({
    week: '2026-W41', ...numbers, onePassPercent: 49, runsPerPullRequest: 2.2, supersededPercent: 29, failedPercent: 20
  });
  assert.match(markdown, /\| Merged pull requests that passed on one run \(%\) \| 100 \| 49 \| \+51 \|/);
  assert.match(markdown, /\| `Tests` runs per merged pull request \| 1 \| 2\.2 \| -1\.2 \|/);
  assert.match(markdown, /\| Pull request runs stopped by a newer push \(%\) \| 0 \| 29 \| -29 \|/);
  assert.match(markdown, /\| Pull request runs that failed \(%\) \| 50 \| 20 \| \+30 \|/);
  assert.match(markdown, /1 pull requests merged into `dev` in the window, 1 with a `Tests` run/);
  assert.match(markdown, /\| failed \| 1 \| 50 \|/);
  assert.deepEqual(readData(markdown), {
    week: '2026-W42', cases: 0, functionalTestMinutes: 0, runnerMinutes: 0, testJobs: 24,
    onePassPercent: 100, runsPerPullRequest: 1, supersededPercent: 0, failedPercent: 50
  });
  // A report from before these numbers still gives the previous week of the others.
  const older = report({ week: '2026-W41', ...numbers });
  assert.match(older, /\| Functional cases \| 0 \| 1000 \| -1000 \|/);
  assert.match(older, /\| Pull request runs that failed \(%\) \| 50 \| - \| - \|/);
  assert.deepEqual(readData(block({ week: '2026-W41', ...numbers })), { week: '2026-W41', ...numbers, ...noShares });
  assert.equal(readData(block({ week: '2026-W41', ...numbers, failedPercent: '20%' })), null);
});

test('a GET that GitHub answers with a server error or a cut-off body is tried again', async (t) => {
  const answers = [];
  t.mock.method(globalThis, 'fetch', async () => answers.shift());
  const reply = (status, body) => new Response(body, { status, headers: { 'x-ratelimit-remaining': '900' } });
  const client = createClient({ token: 't', apiUrl: 'https://api.github.com', repository: 'o/r', retryWaitMs: 0 });
  answers.push(reply(502, '{"message":"Server Error"}'), reply(200, '{"jobs":[{"na'), reply(200, '{"jobs":[]}'));
  assert.deepEqual(await client.json('/actions/runs/1/jobs'), { jobs: [] });
  assert.equal(client.requests(), 3);
  // Three tries at most, and a refusal or a write is never repeated.
  answers.push(reply(502, 'a'), reply(503, 'b'), reply(504, 'c'));
  await assert.rejects(client.json('/actions/runs/2/jobs'), { status: 504 });
  answers.push(reply(403, 'rate limited'));
  await assert.rejects(client.json('/actions/runs/3/jobs'), { status: 403 });
  answers.push(reply(502, 'a'));
  await assert.rejects(client.json('/issues/1/comments', { method: 'POST', body: {} }), { status: 502 });
  assert.equal(client.requests(), 8);
  assert.equal(answers.length, 0);
});

test('a day with more than 1,000 runs is read in parts, since a filtered listing stops at 1,000', async (t) => {
  // 1,500 runs on 2026-10-14, one every 50 seconds from midnight.
  const created = Array.from({ length: 1500 }, (_, index) => Date.UTC(2026, 9, 14) + index * 50_000);
  t.mock.method(globalThis, 'fetch', async (url) => {
    const query = new URL(url).searchParams;
    const [from, to] = query.get('created').split('..').map(Date.parse);
    const matching = created.filter((time) => time >= from && time <= to + 999);
    const perPage = Number(query.get('per_page'));
    const page = Number(query.get('page') || 1);
    const listed = matching.slice(0, 1000).slice((page - 1) * perPage, page * perPage);
    return new Response(JSON.stringify({
      total_count: matching.length,
      workflow_runs: listed.map((time) => ({ id: time, created_at: new Date(time).toISOString() }))
    }), { status: 200 });
  });
  const client = createClient({ token: 't', apiUrl: 'https://api.github.com', repository: 'o/r', retryWaitMs: 0 });
  const { runs, partial } = await listRuns(client, new Date(Date.UTC(2026, 9, 14)), new Date(Date.UTC(2026, 9, 15)));
  assert.equal(partial, null);
  assert.equal(runs.length, 1500);
});
