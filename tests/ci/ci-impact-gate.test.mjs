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

test('gate accepts required success and unrequired skipped or success', () => {
  const impactPlan = plan();
  const skipped = validate(impactPlan, needsFor(impactPlan));
  assert.equal(skipped.ok, true);

  const successfulNeeds = needsFor(impactPlan);
  successfulNeeds.gallery.result = 'success';
  const successful = validate(impactPlan, successfulNeeds);
  assert.equal(successful.ok, true);
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

test('gate rejects wrong profile, workflow SHA, and schema', () => {
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
