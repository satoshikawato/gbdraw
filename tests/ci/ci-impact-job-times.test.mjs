import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';
import { PromotionReadinessError, readRunAttemptJobs } from '../../tools/check-promotion-readiness.mjs';
import { WARNING_SHARE, formatJobTimeReport, jobTimes, workflowJobLimits } from '../../tools/ci-impact-job-times.mjs';
import { createImpactPlan, knownJobsFor } from '../../tools/ci-impact-policy.mjs';
import { runCiImpactCli } from '../../tools/ci-impact.mjs';

const REPOSITORY_ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const WORKFLOW_FILE = '.github/workflows/test.yml';
const REPOSITORY = 'satoshikawato/gbdraw';
const RUN_ID = 4242;
const limitsOfTests = () => workflowJobLimits(readFileSync(resolve(REPOSITORY_ROOT, WORKFLOW_FILE), 'utf8'));
const limitOf = (limits, name) => limits.filter(({ pattern }) => pattern.test(name))
  .sort((left, right) => right.literalLength - left.literalLength)[0];
let nextJobId = 1;
const job = (name, minutes, { status = 'completed', conclusion = 'success' } = {}) => {
  const id = nextJobId++;
  return {
    id,
    run_id: RUN_ID,
    name,
    status,
    conclusion,
    started_at: '2026-11-06T01:00:00Z',
    completed_at: new Date(Date.parse('2026-11-06T01:00:00Z') + minutes * 60000).toISOString().replace('.000', ''),
    html_url: `https://github.com/${REPOSITORY}/actions/runs/${RUN_ID}/job/${id}`
  };
};

test('every job of the Tests workflow has its limit, and each run name finds its job', () => {
  const limits = limitsOfTests();
  const workflow = readFileSync(resolve(REPOSITORY_ROOT, WORKFLOW_FILE), 'utf8');
  const ids = [...workflow.split('\njobs:\n')[1].matchAll(/^ {2}([a-z0-9-]+):\n/gm)].map(([, id]) => id);
  assert.deepEqual(limits.map(({ id }) => id), ids);
  for (const [name, id, minutes] of [
    ['Core PR (Python 3.10)', 'core-pr', 25],
    ['Core (Python 3.12)', 'core', 25],
    ['Playwright functional (shard 2/8)', 'playwright-functional', 45],
    ['recipe acceptance (Python 3.10)', 'acceptance-supported-main', 20],
    ['Promotion audit / sweep-labels', 'promotion-audit-sweeps', 120],
    ['Promotion audit / summary', 'promotion-audit-summary', 10],
    ['Web change budget', 'web-change-budget', 10],
    ['Fail fast (Core PR)', 'fail-fast-core-pr', 2],
    ['PR / gate', 'pr-gate', 5],
    // A job without timeout-minutes has GitHub's 360-minute limit.
    ['Lint', 'lint', 360]
  ]) {
    assert.equal(limitOf(limits, name)?.id, id, name);
    assert.equal(limitOf(limits, name).limitMinutes, minutes, name);
  }
});

test('a job without a name matches its ID and its matrix values; step limits do not count', () => {
  const limits = workflowJobLimits([
    'name: Example', 'on: push', 'jobs:',
    '  build:', '    runs-on: ubuntu-latest', '    timeout-minutes: 12',
    '    strategy:', '      matrix:', '        python: ["3.10"]',
    '    steps:', '      - name: Step', '        timeout-minutes: 3', '        run: true',
    '  "quoted":', '  report:', '    name: "Report ${{ matrix.part }} done"', '    steps:', '      - run: true', ''
  ].join('\n'));
  assert.deepEqual(limits.map(({ id, limitMinutes }) => [id, limitMinutes]), [['build', 12], ['report', 360]]);
  assert.equal(limitOf(limits, 'build (3.10)').id, 'build');
  assert.equal(limitOf(limits, 'build').id, 'build');
  assert.equal(limitOf(limits, 'Report 2 done').id, 'report');
  assert.equal(limitOf(limits, 'Report done'), undefined);
  assert.throws(() => workflowJobLimits('name: none\n'), /no jobs section/);
});

test('the report lists finished jobs at or over the share of their limit', () => {
  const limits = limitsOfTests();
  const times = jobTimes({
    limits,
    jobs: [
      job('Core (Python 3.12)', 21),
      job('Core (Python 3.10)', 25 * WARNING_SHARE),
      job('Core (Python 3.11)', 19.9),
      job('Playwright functional (shard 1/8)', 45, { conclusion: 'cancelled' }),
      job('Core PR (Python 3.10)', 0, { conclusion: 'skipped' }),
      job('PR / gate', 0, { status: 'in_progress', conclusion: null }),
      job('Renamed job', 3)
    ]
  });
  assert.deepEqual(times.map(({ name }) => name), [
    'Core (Python 3.12)', 'Core (Python 3.10)', 'Core (Python 3.11)', 'Playwright functional (shard 1/8)', 'Renamed job'
  ]);
  const { summary, annotations } = formatJobTimeReport({ times });
  assert.deepEqual(annotations, [
    '::warning title=Job near its time limit::Playwright functional (shard 1/8) took 45.0 of its 45 minutes (100%25).',
    '::warning title=Job near its time limit::Core (Python 3.12) took 21.0 of its 25 minutes (84%25).',
    '::warning title=Job near its time limit::Core (Python 3.10) took 20.0 of its 25 minutes (80%25).'
  ]);
  assert.match(summary, /^## Job times against their limits\n\nJobs that used 80% or more of their `timeout-minutes`:\n\n- \[Playwright functional \(shard 1\/8\)\]\(https:[^)]+\): 45\.0 of 45 min \(100%\)\n/);
  assert.match(summary, /\n- \[Core \(Python 3\.10\)\]\([^)]+\): 20\.0 of 25 min \(80%\)\n/);
  assert.doesNotMatch(summary, /3\.11/);
  assert.match(summary, /\n- No workflow job matches: Renamed job\n/);

  const quiet = formatJobTimeReport({ times: jobTimes({ limits, jobs: [job('Core (Python 3.11)', 19.9), job('Lint', 1)] }) });
  assert.deepEqual(quiet.annotations, []);
  assert.match(quiet.summary, /- No job used 80% or more of its `timeout-minutes`; the highest share: Core \(Python 3\.11\) took 19\.9 of its 25 minutes \(79%\)\./);
});

const gatePlan = (profile) => createImpactPlan({
  profile,
  impact: 'full',
  decision: 'full',
  basis: profile === 'pr' ? 'FULL_CHANGE' : 'MANUAL_FULL_RUN',
  changeBaseSha: 'a'.repeat(40),
  changeHeadSha: 'b'.repeat(40),
  workflowSha: 'c'.repeat(40),
  changedPathCount: 1,
  inheritedEvidence: null
});
const runGate = async ({ jobs = [], status = 200, needsResult = 'success' } = {}) => {
  const plan = gatePlan('dev');
  const calls = [];
  let stdout = '';
  let stderr = '';
  let summary = '';
  const exitCode = await runCiImpactCli({
    argv: ['gate'],
    cwd: REPOSITORY_ROOT,
    env: {
      CI_IMPACT_PLAN_JSON: JSON.stringify(plan),
      CI_IMPACT_NEEDS_JSON: JSON.stringify(Object.fromEntries([['ci-impact', { result: 'success' }],
        ...knownJobsFor('dev').map((id) => [id, { result: needsResult }])])),
      CI_IMPACT_EXPECTED_PROFILE: 'dev',
      CI_IMPACT_EXPECTED_WORKFLOW_SHA: 'c'.repeat(40),
      CI_IMPACT_REPOSITORY: REPOSITORY,
      CI_IMPACT_WORKFLOW_FILE: WORKFLOW_FILE,
      GITHUB_RUN_ID: String(RUN_ID),
      GITHUB_RUN_ATTEMPT: '2',
      GITHUB_TOKEN: 'gate-token',
      GITHUB_STEP_SUMMARY: '/tmp/ci-impact-gate-summary'
    },
    stdout: { write: (text) => { stdout += text; } },
    stderr: { write: (text) => { stderr += text; } },
    appendFileImpl: (_path, content) => { summary += content; },
    fetchImpl: async (source, options) => {
      calls.push(source);
      assert.equal(options.headers.Authorization, 'Bearer gate-token');
      assert.equal(new URL(source).pathname, `/repos/${REPOSITORY}/actions/runs/${RUN_ID}/attempts/2/jobs`);
      return {
        ok: status === 200,
        status,
        headers: { get: () => null },
        json: async () => (status === 200 ? { total_count: jobs.length, jobs } : { message: 'Server Error gate-token' })
      };
    }
  });
  return { exitCode, stdout, stderr, summary, calls };
};

test('the dev gate reports its run attempt\'s jobs near their limits before it validates them', async () => {
  const passed = await runGate({ jobs: [job('Core (Python 3.12)', 22), job('Lint', 1)] });
  assert.equal(passed.exitCode, 0, passed.stderr);
  assert.equal(passed.calls.length, 1);
  assert.match(passed.stdout, /^::warning title=Job near its time limit::Core \(Python 3\.12\) took 22\.0 of its 25 minutes \(88%25\)\.\n\{/);
  assert.ok(passed.summary.indexOf('## Job times against their limits') < passed.summary.indexOf('## CI impact aggregate validation'));
  // A gate that fails still reports the job that ran out of time.
  const failed = await runGate({ jobs: [job('Core (Python 3.12)', 25, { conclusion: 'cancelled' })], needsResult: 'cancelled' });
  assert.equal(failed.exitCode, 1);
  assert.match(failed.stderr, /REQUIRED_JOB_NOT_SUCCESSFUL/);
  assert.match(failed.summary, /Core \(Python 3\.12\)\]\([^)]+\): 25\.0 of 25 min \(100%\)/);
});

test('a job-time read failure is a warning, not a gate failure, and redacts the token', async () => {
  const result = await runGate({ status: 500 });
  assert.equal(result.exitCode, 0, result.stderr);
  assert.match(result.stdout, /^::warning title=Job times not read::Job times were not read: \[API_REQUEST_FAILED\] /);
  assert.match(result.summary, /- Job times were not read: \\\[API\\_REQUEST\\_FAILED\\\] /);
  assert.doesNotMatch(`${result.stdout}${result.summary}`, /gate-token/);
});

test('the job reader reads one run attempt and refuses another run\'s jobs', async () => {
  const reply = (body) => ({ ok: true, status: 200, headers: { get: () => null }, json: async () => body });
  const read = (jobs) => readRunAttemptJobs({
    repository: REPOSITORY,
    runId: RUN_ID,
    runAttempt: 1,
    token: 'gate-token',
    fetchImpl: async (source) => {
      const url = new URL(source);
      assert.equal(url.pathname, `/repos/${REPOSITORY}/actions/runs/${RUN_ID}/attempts/1/jobs`);
      assert.equal(url.searchParams.get('per_page'), '100');
      return reply({ total_count: jobs.length, jobs });
    }
  });
  assert.deepEqual((await read([job('Lint', 1)])).map(({ name }) => name), ['Lint']);
  await assert.rejects(read([{ ...job('Lint', 1), run_id: RUN_ID + 1 }]),
    (error) => error instanceof PromotionReadinessError && error.code === 'MALFORMED_JOB');
  await assert.rejects(readRunAttemptJobs({ repository: REPOSITORY, runId: Number.NaN, runAttempt: 1, token: 't' }),
    (error) => error.code === 'INVALID_RUN_IDENTITY');
  await assert.rejects(readRunAttemptJobs({ repository: REPOSITORY, runId: RUN_ID, runAttempt: 1, token: '' }),
    (error) => error.code === 'MISSING_TOKEN');
});

test('a gate without a workflow file makes no job-time request', async () => {
  let calls = 0;
  const exitCode = await runCiImpactCli({
    argv: ['gate'],
    cwd: REPOSITORY_ROOT,
    env: {
      CI_IMPACT_PLAN_JSON: JSON.stringify(gatePlan('dev')),
      CI_IMPACT_NEEDS_JSON: JSON.stringify(Object.fromEntries([['ci-impact', { result: 'success' }],
        ...knownJobsFor('dev').map((id) => [id, { result: 'success' }])])),
      CI_IMPACT_EXPECTED_PROFILE: 'dev',
      CI_IMPACT_EXPECTED_WORKFLOW_SHA: 'c'.repeat(40)
    },
    stdout: { write: () => {} },
    stderr: { write: () => {} },
    fetchImpl: async () => { calls += 1; }
  });
  assert.equal(exitCode, 0);
  assert.equal(calls, 0);
});
