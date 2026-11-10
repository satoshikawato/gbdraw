import assert from 'node:assert/strict';
import test from 'node:test';
import {
  createImpactPlan,
  knownJobsFor,
  validateGateResults
} from '../../tools/ci-impact-policy.mjs';
import { runCiImpactCli } from '../../tools/ci-impact.mjs';

const SHA = Object.freeze({
  base: 'a'.repeat(40),
  head: 'b'.repeat(40),
  workflow: 'c'.repeat(40)
});

const plan = (overrides = {}) => createImpactPlan({
  profile: 'pr',
  impact: 'documentation',
  decision: 'selective',
  basis: 'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE',
  changeBaseSha: SHA.base,
  changeHeadSha: SHA.head,
  workflowSha: SHA.workflow,
  changedPathCount: 1,
  inheritedEvidence: {
    workflowPath: '.github/workflows/test.yml',
    aggregateName: 'Dev staging / gate',
    headSha: SHA.base,
    runId: 101,
    aggregateJobId: 202,
    runUrl: 'https://github.com/satoshikawato/gbdraw/actions/runs/101',
    aggregateJobUrl: 'https://github.com/satoshikawato/gbdraw/actions/runs/101/job/202'
  },
  ...overrides
});

const devPlan = (overrides = {}) => plan({
  profile: 'dev',
  basis: 'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE',
  ...overrides
});

const galleryPlan = (overrides = {}) => plan({
  profile: 'gallery',
  basis: 'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE',
  inheritedEvidence: {
    ...plan().inheritedEvidence,
    workflowPath: '.github/workflows/gallery-publication.yml',
    aggregateName: 'Gallery readiness / gate'
  },
  ...overrides
});

const needsFor = (impactPlan = plan()) => Object.fromEntries([
  ['ci-impact', { result: 'success' }],
  ...knownJobsFor(impactPlan.profile).map((jobId) => [
    jobId,
    { result: impactPlan.requiredJobs.includes(jobId) ? 'success' : 'skipped' }
  ])
]);

const validate = (impactPlan, needs) => validateGateResults({
  plan: impactPlan,
  needs,
  knownJobs: knownJobsFor(impactPlan.profile),
  expected: { profile: impactPlan.profile, workflowSha: SHA.workflow }
});

test('gate accepts metadata, documentation, and full PR routes', () => {
  const routes = [
    plan({ impact: 'metadata' }),
    plan(),
    plan({ impact: 'policy-documentation' }),
    plan({
      impact: 'full',
      decision: 'full',
      basis: 'FULL_CHANGE',
      inheritedEvidence: null
    })
  ];
  assert.deepEqual(routes.map(({ requiredJobs }) => requiredJobs), [
    [],
    ['recipes-standard'],
    ['web-change-budget'],
    [
      'web-change-budget',
      'core-pr',
      'recipes-standard',
      'gallery',
      'examples',
      'lint',
      'web-contracts-pr',
      'web-pr-smoke',
      'playwright-functional',
      'playwright-performance'
    ]
  ]);
  routes.forEach((impactPlan) => {
    assert.equal(validate(impactPlan, needsFor(impactPlan)).ok, true);
  });

  const impactPlan = plan();
  const successfulNeeds = needsFor(impactPlan);
  successfulNeeds.gallery.result = 'success';
  const successful = validate(impactPlan, successfulNeeds);
  assert.equal(successful.ok, true);
});

test('gate accepts metadata, documentation, and full dev staging routes', () => {
  const routes = [
    devPlan({ impact: 'metadata' }),
    devPlan(),
    devPlan({ impact: 'policy-documentation' }),
    devPlan({
      impact: 'full',
      decision: 'full',
      basis: 'FULL_CHANGE',
      inheritedEvidence: null
    })
  ];
  assert.deepEqual(routes.map(({ requiredJobs }) => requiredJobs), [
    [],
    ['recipes-standard'],
    ['web-change-budget'],
    [
      'web-change-budget',
      'core',
      'recipes-standard',
      'gallery',
      'browser',
      'playwright-functional',
      'playwright-performance',
      'lint',
      'losat-cache-browser-acceptance'
    ]
  ]);
  routes.forEach((impactPlan) => {
    assert.equal(validate(impactPlan, needsFor(impactPlan)).ok, true);
  });
});

test('Gallery gate accepts both light routes and requires both full jobs', () => {
  const routes = [
    galleryPlan({ impact: 'metadata' }),
    galleryPlan(),
    galleryPlan({ impact: 'policy-documentation' }),
    galleryPlan({
      impact: 'full',
      decision: 'full',
      basis: 'FULL_CHANGE',
      inheritedEvidence: null
    })
  ];
  assert.deepEqual(routes.map(({ requiredJobs }) => requiredJobs), [
    [],
    [],
    [],
    ['browser', 'performance']
  ]);
  routes.forEach((impactPlan) => {
    assert.equal(validate(impactPlan, needsFor(impactPlan)).ok, true);
  });

  const full = routes[3];
  for (const jobId of ['browser', 'performance']) {
    for (const result of ['skipped', 'failure', 'cancelled']) {
      const needs = needsFor(full);
      needs[jobId].result = result;
      assert.throws(
        () => validate(full, needs),
        /required CI job did not succeed/,
        `${jobId}: ${result}`
      );
    }
  }

  const optionalFailure = needsFor(routes[0]);
  optionalFailure.browser.result = 'failure';
  assert.throws(
    () => validate(routes[0], optionalFailure),
    /unrequired CI job failed unexpectedly/
  );
});

test('Gallery gate rejects the wrong profile, workflow SHA, and schema', () => {
  const impactPlan = galleryPlan();
  const needs = needsFor(impactPlan);
  for (const [candidate, expected, message] of [
    [impactPlan, { profile: 'dev', workflowSha: SHA.workflow }, /profile does not match/],
    [impactPlan, { profile: 'gallery', workflowSha: 'd'.repeat(40) }, /workflow SHA does not match/],
    [{ ...impactPlan, schemaVersion: 999 }, {
      profile: 'gallery', workflowSha: SHA.workflow
    }, /schema version is not supported/]
  ]) {
    assert.throws(() => validateGateResults({
      plan: candidate,
      needs,
      knownJobs: knownJobsFor('gallery'),
      expected
    }), message);
  }
});

test('dev gate treats each matrix aggregate as one required job result', () => {
  const full = devPlan({
    impact: 'full',
    decision: 'full',
    basis: 'FULL_CHANGE',
    inheritedEvidence: null
  });
  const fullNeeds = needsFor(full);
  assert.equal(validate(full, fullNeeds).ok, true);
  for (const result of ['skipped', 'failure']) {
    const invalid = needsFor(full);
    invalid['playwright-functional'].result = result;
    assert.throws(
      () => validate(full, invalid),
      /required CI job did not succeed/,
      result
    );
  }

  const metadata = devPlan({ impact: 'metadata' });
  const optionalSuccess = needsFor(metadata);
  optionalSuccess['playwright-functional'].result = 'success';
  assert.equal(validate(metadata, optionalSuccess).ok, true);
  const optionalFailure = needsFor(metadata);
  optionalFailure['playwright-functional'].result = 'failure';
  assert.throws(() => validate(metadata, optionalFailure), /unrequired CI job failed unexpectedly/);
});

test('dev gate rejects the wrong current SHA and malformed inherited evidence', () => {
  const impactPlan = devPlan();
  assert.throws(() => validateGateResults({
    plan: impactPlan,
    needs: needsFor(impactPlan),
    knownJobs: knownJobsFor('dev'),
    expected: { profile: 'dev', workflowSha: 'd'.repeat(40) }
  }), /workflow SHA does not match/);

  const malformed = {
    ...impactPlan,
    inheritedEvidence: { ...impactPlan.inheritedEvidence, runId: '101' }
  };
  assert.throws(() => validateGateResults({
    plan: malformed,
    needs: needsFor(impactPlan),
    knownJobs: knownJobsFor('dev'),
    expected: { profile: 'dev', workflowSha: SHA.workflow }
  }), /Inherited evidence identity is invalid/);
});

test('gate rejects required skipped, failure, and cancellation', () => {
  for (const result of ['skipped', 'failure', 'cancelled']) {
    const impactPlan = plan();
    const needs = needsFor(impactPlan);
    needs['recipes-standard'].result = result;
    assert.throws(
      () => validate(impactPlan, needs),
      /required CI job did not succeed/,
      result
    );
  }
});

test('gate rejects planner failures and unrequired failures', () => {
  const impactPlan = plan();
  const plannerFailure = needsFor(impactPlan);
  plannerFailure['ci-impact'].result = 'failure';
  assert.throws(() => validate(impactPlan, plannerFailure), /planner did not succeed/);

  for (const result of ['failure', 'cancelled', 'timed_out']) {
    const needs = needsFor(impactPlan);
    needs.gallery.result = result;
    assert.throws(
      () => validate(impactPlan, needs),
      /unrequired CI job failed unexpectedly/,
      result
    );
  }
});

test('gate rejects missing expected jobs', () => {
  const impactPlan = plan();
  const needs = needsFor(impactPlan);
  delete needs.lint;
  assert.throws(() => validate(impactPlan, needs), /result is missing or invalid/);
});

test('unknown jobs may succeed or skip but may not fail', () => {
  const impactPlan = plan();
  for (const result of ['success', 'skipped']) {
    const needs = { ...needsFor(impactPlan), future: { result } };
    assert.equal(validate(impactPlan, needs).ok, true);
  }
  const failed = { ...needsFor(impactPlan), future: { result: 'failure' } };
  assert.throws(() => validate(impactPlan, failed), /unknown CI job failed unexpectedly/);
});

test('trusted gate rejects wrong profile, workflow SHA, and helper schema', () => {
  const impactPlan = plan();
  const needs = needsFor(impactPlan);
  assert.throws(() => validateGateResults({
    plan: impactPlan,
    needs,
    knownJobs: knownJobsFor('pr'),
    expected: { profile: 'dev', workflowSha: SHA.workflow }
  }), /profile does not match/);
  assert.throws(() => validateGateResults({
    plan: impactPlan,
    needs,
    knownJobs: knownJobsFor('pr'),
    expected: { profile: 'pr', workflowSha: 'd'.repeat(40) }
  }), /workflow SHA does not match/);
  assert.throws(() => validateGateResults({
    plan: { ...impactPlan, schemaVersion: 999 },
    needs,
    knownJobs: knownJobsFor('pr'),
    expected: { profile: 'pr', workflowSha: SHA.workflow }
  }), /schema version is not supported/);
});

const textWriter = () => {
  let value = '';
  return {
    stream: { write: (chunk) => { value += chunk; } },
    value: () => value
  };
};

// The Actions API as `PR / gate` reads it: `Tests` push runs on dev, and earlier attempts of a
// run being rerun.
const REPOSITORY = 'satoshikawato/gbdraw';
const devRun = (id, conclusion, minute, { status = 'completed', attempt = 1 } = {}) => ({
  id,
  run_number: id,
  run_attempt: attempt,
  workflow_id: 501,
  path: '.github/workflows/test.yml',
  event: 'push',
  head_branch: 'dev',
  head_sha: String(id % 10).repeat(40),
  repository: { full_name: REPOSITORY },
  status,
  conclusion,
  url: `https://api.github.com/repos/${REPOSITORY}/actions/runs/${id}`,
  html_url: `https://github.com/${REPOSITORY}/actions/runs/${id}`,
  created_at: `2026-11-05T10:${String(minute).padStart(2, '0')}:00Z`,
  updated_at: `2026-11-05T11:${String(minute).padStart(2, '0')}:00Z`
});
const reply = (status, body) => ({
  ok: status === 200,
  status,
  headers: { get: () => null },
  json: async () => body
});
// `attempts` maps `<run id>/<attempt>` to that attempt's run object; a missing key is a 404.
const devRunsFetch = (runs, { status = 200, attempts = {} } = {}) => {
  const calls = [];
  const fetchImpl = async (source, options) => {
    calls.push(source);
    assert.equal(options.headers.Authorization, 'Bearer gate-token');
    const url = new URL(source);
    const attempt = url.pathname.match(new RegExp(`^/repos/${REPOSITORY}/actions/runs/(\\d+)/attempts/(\\d+)$`));
    if (attempt) {
      const key = `${attempt[1]}/${attempt[2]}`;
      return Object.hasOwn(attempts, key) ? reply(200, attempts[key]) : reply(404, { message: 'Not Found' });
    }
    assert.equal(url.pathname, `/repos/${REPOSITORY}/actions/workflows/test.yml/runs`);
    assert.equal(url.searchParams.get('branch'), 'dev');
    assert.equal(url.searchParams.get('event'), 'push');
    // Runs still running are listed too: a rerun keeps its previous result.
    assert.equal(url.searchParams.has('status'), false);
    return status === 200
      ? reply(200, { total_count: runs.length, workflow_runs: runs })
      : reply(status, { message: 'Server Error gate-token' });
  };
  return { fetchImpl, calls };
};
const runPrGate = async ({ runs, labels = [], status, attempts, impactPlan = plan() }) => {
  const stdout = textWriter();
  const stderr = textWriter();
  const stub = devRunsFetch(runs, { status, attempts });
  let summary = '';
  const exitCode = await runCiImpactCli({
    argv: ['gate'],
    env: {
      CI_IMPACT_PLAN_JSON: JSON.stringify(impactPlan),
      CI_IMPACT_NEEDS_JSON: JSON.stringify(needsFor(impactPlan)),
      CI_IMPACT_EXPECTED_PROFILE: impactPlan.profile,
      CI_IMPACT_EXPECTED_WORKFLOW_SHA: SHA.workflow,
      CI_IMPACT_REPOSITORY: REPOSITORY,
      CI_IMPACT_PR_LABELS_JSON: JSON.stringify(labels),
      GITHUB_TOKEN: 'gate-token',
      GITHUB_STEP_SUMMARY: '/tmp/ci-impact-gate-summary'
    },
    stdout: stdout.stream,
    stderr: stderr.stream,
    appendFileImpl: (_path, content) => { summary += content; },
    fetchImpl: stub.fetchImpl
  });
  return { exitCode, stdout: stdout.value(), stderr: stderr.value(), summary, calls: stub.calls };
};

test('gate CLI rejects malformed plan and needs JSON', async () => {
  for (const [field, value] of [
    ['CI_IMPACT_PLAN_JSON', '{'],
    ['CI_IMPACT_NEEDS_JSON', '[not-json]']
  ]) {
    const stdout = textWriter();
    const stderr = textWriter();
    const impactPlan = plan();
    const status = await runCiImpactCli({
      argv: ['gate'],
      env: {
        CI_IMPACT_PLAN_JSON: JSON.stringify(impactPlan),
        CI_IMPACT_NEEDS_JSON: JSON.stringify(needsFor(impactPlan)),
        CI_IMPACT_EXPECTED_PROFILE: 'pr',
        CI_IMPACT_EXPECTED_WORKFLOW_SHA: SHA.workflow,
        [field]: value
      },
      stdout: stdout.stream,
      stderr: stderr.stream
    });
    assert.equal(status, 1);
    assert.equal(stdout.value(), '');
    assert.match(stderr.value(), /MALFORMED_JSON/);
  }
});

test('gate CLI rejects a missing plan', async () => {
  const stderr = textWriter();
  const status = await runCiImpactCli({
    argv: ['gate'],
    env: {
      CI_IMPACT_NEEDS_JSON: JSON.stringify(needsFor()),
      CI_IMPACT_EXPECTED_PROFILE: 'pr',
      CI_IMPACT_EXPECTED_WORKFLOW_SHA: SHA.workflow
    },
    stdout: textWriter().stream,
    stderr: stderr.stream
  });
  assert.equal(status, 1);
  assert.match(stderr.value(), /MISSING_ENVIRONMENT/);
});

test('gate CLI emits a passing summary for a valid payload', async () => {
  const { exitCode, stdout, stderr, summary } = await runPrGate({ runs: [devRun(901, 'success', 10)] });
  assert.equal(exitCode, 0, stderr);
  assert.equal(JSON.parse(stdout).ok, true);
  assert.match(summary, /Gate result: pass/);
  assert.match(summary, /`recipes-standard`: `success` \(required\)/);
  assert.match(summary, /`gallery`: `skipped` \(not required\)/);
});

test('PR gate fails while dev is red unless the pull request carries the fix label', async () => {
  for (const conclusion of ['failure', 'timed_out', 'startup_failure']) {
    const runs = [devRun(902, conclusion, 20), devRun(901, 'success', 10)];
    const red = await runPrGate({ runs, labels: ['architecture-change'] });
    assert.equal(red.exitCode, 1, conclusion);
    assert.equal(red.stdout, '');
    assert.match(red.stderr, /\[DEV_RED\]/);
    assert.match(red.stderr, /https:\/\/github\.com\/satoshikawato\/gbdraw\/actions\/runs\/902/);
    assert.match(red.stderr, new RegExp('2'.repeat(40)));
    assert.match(red.stderr, /label dev-red-fix/);
    assert.match(red.stderr, /rerun the PR \/ gate job after dev turns green/);
    assert.match(red.summary, new RegExp(`Dev state: \`red\` \\(\`${conclusion}\`\\)`));
    assert.doesNotMatch(red.summary, /Gate result: pass/);

    const fix = await runPrGate({ runs, labels: ['dev-red-fix'] });
    assert.equal(fix.exitCode, 0, fix.stderr);
    assert.match(fix.summary, /Dev state: `red`[\s\S]*Label `dev-red-fix`: yes[\s\S]*Gate result: pass/);
  }
});

test('PR gate skips dev runs whose conclusion says nothing', async () => {
  const older = [devRun(901, 'success', 10)];
  for (const conclusion of ['cancelled', 'skipped', 'neutral', 'stale', 'action_required']) {
    // Listed oldest first: the reader orders runs by creation time itself.
    const red = await runPrGate({ runs: [devRun(902, 'failure', 20), devRun(903, conclusion, 30)] });
    assert.equal(red.exitCode, 1, conclusion);
    assert.match(red.stderr, /\[DEV_RED\][\s\S]*runs\/902/, conclusion);
    const green = await runPrGate({ runs: [devRun(903, conclusion, 30), ...older] });
    assert.equal(green.exitCode, 0, `${conclusion}: ${green.stderr}`);
    assert.match(green.summary, /Dev state: `green` \(`success`\) — \[run 901\]/);
  }
  const none = await runPrGate({ runs: [devRun(903, 'cancelled', 30)] });
  assert.equal(none.exitCode, 1);
  assert.match(none.stderr, /\[DEV_STATE_UNAVAILABLE\][\s\S]*NO_DEV_STATE/);
  const unknown = await runPrGate({ runs: [devRun(903, 'mystery', 30)] });
  assert.match(unknown.stderr, /\[DEV_STATE_UNAVAILABLE\][\s\S]*UNKNOWN_RUN_CONCLUSION/);
});

test('PR gate fails closed when the dev state cannot be read, and redacts the token', async () => {
  const failed = await runPrGate({ runs: [], status: 500 });
  assert.equal(failed.exitCode, 1);
  assert.equal(failed.stdout, '');
  assert.match(failed.stderr, /\[DEV_STATE_UNAVAILABLE\][\s\S]*"evidenceCode":"API_REQUEST_FAILED"/);
  assert.match(failed.summary, /Dev state: unavailable — `API_REQUEST_FAILED`/);
  assert.doesNotMatch(failed.stderr + failed.summary, /gate-token/);

  const malformedLabels = textWriter();
  const unlabeled = await runCiImpactCli({
    argv: ['gate'],
    env: {
      CI_IMPACT_PLAN_JSON: JSON.stringify(plan()),
      CI_IMPACT_NEEDS_JSON: JSON.stringify(needsFor(plan())),
      CI_IMPACT_EXPECTED_PROFILE: 'pr',
      CI_IMPACT_EXPECTED_WORKFLOW_SHA: SHA.workflow,
      CI_IMPACT_REPOSITORY: REPOSITORY,
      CI_IMPACT_PR_LABELS_JSON: '{"name":"dev-red-fix"}',
      GITHUB_TOKEN: 'gate-token'
    },
    stdout: textWriter().stream,
    stderr: malformedLabels.stream,
    fetchImpl: devRunsFetch([devRun(901, 'success', 10)]).fetchImpl
  });
  assert.equal(unlabeled, 1);
  assert.match(malformedLabels.value(), /\[INVALID_LABELS\]/);
});

test('PR gate counts a run being rerun by its previous attempt', async () => {
  const green = [devRun(901, 'success', 10)];
  const rerun = devRun(902, null, 20, { status: 'in_progress', attempt: 2 });
  const failedAttempt = { ...devRun(902, 'failure', 20), status: 'completed' };
  // A failed run stays red while it reruns, although an older run is green.
  const held = await runPrGate({ runs: [rerun, ...green], attempts: { '902/1': failedAttempt } });
  assert.equal(held.exitCode, 1);
  assert.match(held.stderr, /\[DEV_RED\][\s\S]*runs\/902[\s\S]*attempt 1, and its rerun has not finished/);
  assert.match(held.summary, /Dev state: `red` \(`failure`\) — \[run 902\]\([^)]*\) attempt 1 \(a rerun is in progress\)/);
  assert.ok(held.calls.some((call) => call.endsWith('/actions/runs/902/attempts/1')));
  const fix = await runPrGate({ runs: [rerun, ...green], labels: ['dev-red-fix'], attempts: { '902/1': failedAttempt } });
  assert.equal(fix.exitCode, 0, fix.stderr);

  // A green run being rerun stays green; a cancelled previous attempt says nothing.
  const greenRerun = await runPrGate({ runs: [rerun, devRun(901, 'failure', 10)], attempts: { '902/1': { ...failedAttempt, conclusion: 'success' } } });
  assert.equal(greenRerun.exitCode, 0, greenRerun.stderr);
  const cancelledAttempt = await runPrGate({ runs: [rerun, devRun(901, 'failure', 10)], attempts: { '902/1': { ...failedAttempt, conclusion: 'cancelled' } } });
  assert.match(cancelledAttempt.stderr, /\[DEV_RED\][\s\S]*runs\/901/);

  // A first attempt that is still running has no result yet, and its attempt is not read.
  const first = await runPrGate({ runs: [devRun(903, null, 30, { status: 'queued' }), ...green] });
  assert.equal(first.exitCode, 0, first.stderr);
  assert.equal(first.calls.length, 1);

  // An unreadable or mismatched previous attempt fails closed.
  const missing = await runPrGate({ runs: [rerun, ...green] });
  assert.match(missing.stderr, /\[DEV_STATE_UNAVAILABLE\][\s\S]*API_REQUEST_FAILED/);
  const wrongAttempt = await runPrGate({ runs: [rerun, ...green], attempts: { '902/1': { ...failedAttempt, run_attempt: 2 } } });
  assert.match(wrongAttempt.stderr, /\[DEV_STATE_UNAVAILABLE\][\s\S]*RUN_IDENTITY_MISMATCH/);
  const unfinished = await runPrGate({ runs: [rerun, ...green], attempts: { '902/1': { ...failedAttempt, status: 'in_progress' } } });
  assert.match(unfinished.stderr, /\[DEV_STATE_UNAVAILABLE\][\s\S]*RUN_IDENTITY_MISMATCH/);
});

test('PR gate passes over cancelled attempts to the newest attempt with a result', async () => {
  const green = [devRun(901, 'success', 10)];
  const attemptOf = (id, attempt, conclusion) => ({ ...devRun(id, conclusion, 20, { attempt }), status: 'completed' });
  // (a) A rerun that ended cancelled leaves the run red by its failed first attempt.
  const cancelledRerun = await runPrGate({
    runs: [devRun(903, null, 30, { status: 'in_progress' }), devRun(902, 'cancelled', 20, { attempt: 2 }), ...green],
    attempts: { '902/1': attemptOf(902, 1, 'failure') }
  });
  assert.equal(cancelledRerun.exitCode, 1);
  assert.match(cancelledRerun.stderr, /\[DEV_RED\][\s\S]*runs\/902[\s\S]*ended failure in attempt 1\./);
  assert.match(cancelledRerun.summary, /\[run 902\]\([^)]*\) attempt 1 on/);

  // (b) A third attempt running after a cancelled second: the first attempt's failure decides.
  const thirdAttempt = await runPrGate({
    runs: [devRun(902, null, 20, { status: 'in_progress', attempt: 3 }), ...green],
    attempts: { '902/2': attemptOf(902, 2, 'cancelled'), '902/1': attemptOf(902, 1, 'failure') }
  });
  assert.equal(thirdAttempt.exitCode, 1);
  assert.match(thirdAttempt.stderr, /\[DEV_RED\][\s\S]*runs\/902[\s\S]*attempt 1, and its rerun has not finished/);
  assert.ok(thirdAttempt.calls.some((call) => call.endsWith('/actions/runs/902/attempts/2')));
  assert.ok(thirdAttempt.calls.some((call) => call.endsWith('/actions/runs/902/attempts/1')));

  // Every attempt without a result: the run is skipped and the older run decides.
  const allCancelled = await runPrGate({
    runs: [devRun(902, 'cancelled', 20, { attempt: 3 }), ...green],
    attempts: { '902/2': attemptOf(902, 2, 'skipped'), '902/1': attemptOf(902, 1, 'cancelled') }
  });
  assert.equal(allCancelled.exitCode, 0, allCancelled.stderr);
  assert.match(allCancelled.summary, /Dev state: `green` \(`success`\) — \[run 901\]/);

  // An earlier attempt that cannot be read fails closed.
  const unreadable = await runPrGate({
    runs: [devRun(902, 'cancelled', 20, { attempt: 3 }), ...green],
    attempts: { '902/2': attemptOf(902, 2, 'cancelled') }
  });
  assert.match(unreadable.stderr, /\[DEV_STATE_UNAVAILABLE\][\s\S]*API_REQUEST_FAILED/);
});

test('dev, release, and Gallery gates make no dev-state request', async () => {
  for (const impactPlan of [
    devPlan(),
    plan({ profile: 'release', impact: 'full', decision: 'full', basis: 'MANUAL_FULL_RUN', inheritedEvidence: null }),
    galleryPlan()
  ]) {
    const outcome = await runPrGate({ runs: [devRun(902, 'failure', 20)], impactPlan });
    assert.equal(outcome.exitCode, 0, `${impactPlan.profile}: ${outcome.stderr}`);
    assert.deepEqual(outcome.calls, [], impactPlan.profile);
    assert.doesNotMatch(outcome.summary, /Dev state/, impactPlan.profile);
  }
});

test('release aggregate rejects any skipped or failed exhaustive matrix', () => {
  const release = plan({ profile: 'release', impact: 'full', decision: 'full', basis: 'MANUAL_FULL_RUN', inheritedEvidence: null });
  assert.equal(validate(release, needsFor(release)).ok, true);
  assert.ok(release.requiredJobs.includes('vibrio-generate-release'), 'release evidence must cover Vibrio generation');
  for (const jobId of release.requiredJobs) {
    for (const result of ['skipped', 'failure', 'cancelled']) {
      const needs = needsFor(release);
      needs[jobId].result = result;
      assert.throws(() => validate(release, needs), /required CI job did not succeed/);
    }
  }
});

test('documentation-only PR gate requires selected checks without claiming inherited staging', () => {
  for (const impact of ['documentation', 'policy-documentation']) {
    const documentary = plan({ impact, basis: 'DOCUMENTATION_ONLY_PR', inheritedEvidence: null });
    const needs = needsFor(documentary);
    const outcome = validate(documentary, needs);
    assert.equal(outcome.ok, true);
    assert.equal(outcome.inheritedEvidence, false);
    for (const result of ['failure', 'cancelled', 'skipped']) {
      const job = documentary.requiredJobs[0];
      assert.throws(() => validate(documentary, { ...needs, [job]: { result } }),
        { code: 'REQUIRED_JOB_NOT_SUCCESSFUL' });
    }
  }
});

test('PR gate requires Gallery parity for full and every selective smoke route', () => {
  const routes = [
    plan({ impact: 'full', decision: 'full', basis: 'FULL_CHANGE', inheritedEvidence: null }),
    ...['python-core', 'renderer', 'web-runtime', 'session-persistence', 'gallery', 'losat-integration'].map((impact) => plan({ impact }))
  ];
  for (const impactPlan of routes) {
    assert.ok(impactPlan.requiredJobs.includes('gallery'), impactPlan.impact);
    assert.equal(validate(impactPlan, needsFor(impactPlan)).ok, true);
    for (const result of ['skipped', 'failure', 'cancelled']) {
      const needs = needsFor(impactPlan);
      needs.gallery.result = result;
      assert.throws(() => validate(impactPlan, needs), /required CI job did not succeed/, `${impactPlan.impact}: ${result}`);
    }
    const missing = needsFor(impactPlan);
    delete missing.gallery;
    assert.throws(() => validate(impactPlan, missing), /result is missing or invalid/);
  }
});

test('PR gate requires functional Playwright exactly when the plan selects it', () => {
  const selected = [
    plan({ impact: 'web-runtime' }),
    plan({ impact: 'full', decision: 'full', basis: 'FULL_CHANGE', inheritedEvidence: null }),
    plan({ impact: 'gallery', capabilities: ['renderer', 'gallery'], decision: 'full', basis: 'INHERITED_EVIDENCE_UNAVAILABLE', inheritedEvidence: null })
  ];
  for (const impactPlan of selected) {
    assert.ok(impactPlan.requiredJobs.includes('playwright-functional'), impactPlan.impact);
    assert.equal(validate(impactPlan, needsFor(impactPlan)).ok, true);
    for (const result of ['skipped', 'failure', 'cancelled']) {
      const needs = needsFor(impactPlan);
      needs['playwright-functional'].result = result;
      assert.throws(() => validate(impactPlan, needs), { code: 'REQUIRED_JOB_NOT_SUCCESSFUL' }, result);
    }
    const missing = needsFor(impactPlan);
    delete missing['playwright-functional'];
    assert.throws(() => validate(impactPlan, missing), { code: 'MISSING_OR_INVALID_JOB' });
  }
  const notSelected = [
    plan({ impact: 'python-core' }),
    plan({ impact: 'ci-only', decision: 'full', basis: 'FULL_CHANGE', inheritedEvidence: null })
  ];
  for (const impactPlan of notSelected) {
    assert.ok(!impactPlan.requiredJobs.includes('playwright-functional'), impactPlan.impact);
    const needs = needsFor(impactPlan);
    assert.equal(needs['playwright-functional'].result, 'skipped');
    assert.equal(validate(impactPlan, needs).ok, true);
    needs['playwright-functional'].result = 'failure';
    assert.throws(() => validate(impactPlan, needs), { code: 'UNREQUIRED_JOB_FAILED' });
  }
});

test('dev gate accepts identical-tree and leaf-test plans only when their selected jobs succeed', () => {
  const identical = devPlan({
    impact: 'none',
    capabilities: ['none'],
    basis: 'IDENTICAL_TREE_WITH_DIRECT_PARENT_EVIDENCE',
    changedPathCount: 0
  });
  assert.deepEqual(identical.requiredJobs, []);
  assert.equal(validate(identical, needsFor(identical)).ok, true);

  const leaf = devPlan({
    impact: 'web-runtime',
    capabilities: ['documentation', 'web-runtime'],
    basis: 'LEAF_TEST_CHANGE_WITH_DIRECT_PARENT_EVIDENCE',
    changedPathCount: 2,
    leafTests: [{ path: 'tests/web/right-drawer.playwright.spec.js', kind: 'functional', jobs: ['playwright-functional'] }]
  });
  assert.deepEqual(leaf.requiredJobs, ['recipes-standard', 'playwright-functional']);
  assert.equal(validate(leaf, needsFor(leaf)).ok, true);
  assert.throws(() => validate(leaf, { ...needsFor(leaf), 'playwright-functional': { result: 'skipped' } }),
    { code: 'REQUIRED_JOB_NOT_SUCCESSFUL' });
  const roundTrip = JSON.parse(JSON.stringify(leaf));
  assert.equal(validate(roundTrip, needsFor(leaf)).ok, true);
  assert.throws(() => validate({ ...roundTrip, leafTests: null }, needsFor(leaf)), { code: 'INVALID_LEAF_TESTS' });
});

test('dev gate accepts a merged pull request tree plan only when the jobs it did not inherit succeed', () => {
  const merged = devPlan({
    impact: 'web-runtime',
    capabilities: ['web-runtime'],
    basis: 'IDENTICAL_TREE_WITH_PR_EVIDENCE',
    inheritedEvidence: { ...devPlan().inheritedEvidence, aggregateName: 'PR / gate', headSha: 'f'.repeat(40) }
  });
  assert.deepEqual(merged.requiredJobs, [
    'web-change-budget', 'core', 'recipes-standard', 'browser', 'playwright-performance', 'lint',
    'losat-cache-browser-acceptance'
  ]);
  const needs = needsFor(merged);
  assert.equal(needs['playwright-functional'].result, 'skipped');
  assert.equal(validate(merged, needs).ok, true);
  assert.throws(() => validate(merged, { ...needs, core: { result: 'skipped' } }),
    { code: 'REQUIRED_JOB_NOT_SUCCESSFUL' });
  assert.throws(() => validate(merged, { ...needs, gallery: { result: 'failure' } }),
    { code: 'UNREQUIRED_JOB_FAILED' });
});
