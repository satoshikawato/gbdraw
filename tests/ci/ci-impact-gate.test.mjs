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
    [
      'web-change-budget',
      'core-pr',
      'recipes-standard',
      'gallery',
      'lint',
      'web-pr-smoke'
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
    [
      'web-change-budget',
      'core',
      'recipes-standard',
      'gallery',
      'browser',
      'playwright-functional',
      'playwright-performance',
      'acceptance-supported-main',
      'slow-main',
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
    ['browser', 'performance']
  ]);
  routes.forEach((impactPlan) => {
    assert.equal(validate(impactPlan, needsFor(impactPlan)).ok, true);
  });

  const full = routes[2];
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
    [{ ...impactPlan, schemaVersion: 2 }, {
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
    plan: { ...impactPlan, schemaVersion: 2 },
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
  const stdout = textWriter();
  const stderr = textWriter();
  const impactPlan = plan();
  let summary = '';
  const status = await runCiImpactCli({
    argv: ['gate'],
    env: {
      CI_IMPACT_PLAN_JSON: JSON.stringify(impactPlan),
      CI_IMPACT_NEEDS_JSON: JSON.stringify(needsFor(impactPlan)),
      CI_IMPACT_EXPECTED_PROFILE: 'pr',
      CI_IMPACT_EXPECTED_WORKFLOW_SHA: SHA.workflow,
      GITHUB_STEP_SUMMARY: '/tmp/ci-impact-gate-summary'
    },
    stdout: stdout.stream,
    stderr: stderr.stream,
    appendFileImpl: (_path, content) => { summary += content; }
  });
  assert.equal(status, 0, stderr.value());
  assert.equal(JSON.parse(stdout.value()).ok, true);
  assert.match(summary, /Gate result: pass/);
  assert.match(summary, /`recipes-standard`: `success` \(required\)/);
  assert.match(summary, /`gallery`: `skipped` \(not required\)/);
});
