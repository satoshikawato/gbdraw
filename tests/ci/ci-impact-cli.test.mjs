import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';
import { PromotionReadinessError } from '../../tools/check-promotion-readiness.mjs';
import {
  buildImpactPlan,
  parseNameStatusZ,
  readPlanConfiguration,
  runCiImpactCli
} from '../../tools/ci-impact.mjs';

const REPOSITORY_ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const SHA = Object.freeze({
  base: 'a'.repeat(40),
  head: 'b'.repeat(40),
  workflow: 'c'.repeat(40)
});

const environment = (overrides = {}) => ({
  CI_IMPACT_PROFILE: 'pr',
  CI_IMPACT_EVENT_NAME: 'pull_request',
  CI_IMPACT_REPOSITORY: 'satoshikawato/gbdraw',
  CI_IMPACT_REPOSITORY_ROOT: REPOSITORY_ROOT,
  CI_IMPACT_CHANGE_BASE_SHA: SHA.base,
  CI_IMPACT_CHANGE_HEAD_SHA: SHA.head,
  CI_IMPACT_WORKFLOW_SHA: SHA.workflow,
  CI_IMPACT_ARCHITECTURE_CHANGE: 'false',
  GITHUB_TOKEN: 'test-token',
  ...overrides
});

const configuration = (overrides = {}) => readPlanConfiguration(
  environment(overrides),
  REPOSITORY_ROOT
);

const gitResult = (...tokens) => ({
  status: 0,
  stdout: Buffer.from(`${tokens.join('\0')}\0`),
  stderr: Buffer.alloc(0)
});

const successfulEvidence = ({
  workflowPath = '.github/workflows/test.yml',
  aggregateName = 'Dev staging / gate'
} = {}) => ({
  workflow: { path: workflowPath },
  run: {
    id: 101,
    headSha: SHA.base,
    url: 'https://github.com/satoshikawato/gbdraw/actions/runs/101'
  },
  aggregateJob: {
    id: 202,
    name: aggregateName,
    url: 'https://github.com/satoshikawato/gbdraw/actions/runs/101/job/202'
  }
});

const textWriter = () => {
  let value = '';
  return {
    stream: { write: (chunk) => { value += chunk; } },
    value: () => value
  };
};

test('NUL name-status parsing preserves newlines and rename endpoints', () => {
  const changes = parseNameStatusZ(Buffer.from(
    'M\0docs/line\nbreak.md\0R100\0docs/old.md\0gbdraw/new.md\0'
  ));
  assert.deepEqual(changes, [
    { status: 'M', paths: ['docs/line\nbreak.md'] },
    { status: 'R100', paths: ['docs/old.md', 'gbdraw/new.md'] }
  ]);
  assert.throws(() => parseNameStatusZ(Buffer.from('M\0docs/no-terminator.md')), /NUL-terminated/);
});

test('PR planning uses a three-dot diff and direct base evidence', async () => {
  let gitArgs;
  let evidenceArguments;
  const outcome = await buildImpactPlan({
    configuration: configuration(),
    token: 'test-token',
    runGitImpl: (_root, args) => {
      gitArgs = args;
      return gitResult('M', 'docs/FAQ.md');
    },
    verifyWorkflowEvidenceImpl: async (args) => {
      evidenceArguments = args;
      return successfulEvidence();
    }
  });
  assert.deepEqual(gitArgs, [
    'diff',
    '--name-status',
    '-z',
    '--find-renames',
    `${SHA.base}...${SHA.head}`,
    '--'
  ]);
  assert.equal(evidenceArguments.expectedHeadSha, SHA.base);
  assert.equal(evidenceArguments.workflowPath, '.github/workflows/test.yml');
  assert.equal(outcome.plan.impact, 'documentation');
  assert.equal(outcome.plan.decision, 'selective');
  assert.equal(outcome.plan.basis, 'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE');
  assert.deepEqual(outcome.plan.requiredJobs, ['recipes-standard']);
});

test('dev and Gallery planning use a two-commit diff and direct parent evidence', async () => {
  for (const [profile, workflowPath, aggregateName] of [
    ['dev', '.github/workflows/test.yml', 'Dev staging / gate'],
    ['gallery', '.github/workflows/gallery-publication.yml', 'Gallery readiness / gate']
  ]) {
    let gitArgs;
    const outcome = await buildImpactPlan({
      configuration: configuration({
        CI_IMPACT_PROFILE: profile,
        CI_IMPACT_EVENT_NAME: 'push'
      }),
      token: 'test-token',
      runGitImpl: (_root, args) => {
        gitArgs = args;
        return gitResult('M', '.gitignore');
      },
      verifyWorkflowEvidenceImpl: async (args) => {
        assert.equal(args.workflowPath, workflowPath);
        assert.equal(args.expectedAggregateName, aggregateName);
        return successfulEvidence({ workflowPath, aggregateName });
      }
    });
    assert.deepEqual(gitArgs.slice(-3), [SHA.base, SHA.head, '--']);
    assert.equal(outcome.plan.basis, 'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE');
    assert.deepEqual(outcome.plan.requiredJobs, []);
  }
});

test('evidence verifier failures fall back to the full profile', async () => {
  const unavailableToken = 'unavailable-token';
  const outcome = await buildImpactPlan({
    configuration: configuration(),
    token: unavailableToken,
    runGitImpl: () => gitResult('M', '.gitignore'),
    verifyWorkflowEvidenceImpl: async () => {
      throw new PromotionReadinessError(
        'API_REQUEST_FAILED',
        `GitHub Actions API request failed for ${unavailableToken}.`
      );
    }
  });
  assert.equal(outcome.plan.impact, 'metadata');
  assert.equal(outcome.plan.decision, 'full');
  assert.equal(outcome.plan.basis, 'INHERITED_EVIDENCE_UNAVAILABLE');
  assert.deepEqual(outcome.plan.requiredJobs, [
    'web-change-budget',
    'core-pr',
    'recipes-standard',
    'gallery',
    'lint',
    'web-pr-smoke'
  ]);
  assert.deepEqual(outcome.evidenceFailure, {
    code: 'API_REQUEST_FAILED',
    reason: 'GitHub Actions API request failed for [REDACTED].'
  });
});

test('full changes never query inherited evidence', async () => {
  let evidenceCalls = 0;
  const outcome = await buildImpactPlan({
    configuration: configuration(),
    token: 'test-token',
    runGitImpl: () => gitResult('M', 'tools/ci-impact.mjs'),
    verifyWorkflowEvidenceImpl: async () => {
      evidenceCalls += 1;
      return successfulEvidence();
    }
  });
  assert.equal(evidenceCalls, 0);
  assert.equal(outcome.plan.impact, 'full');
  assert.equal(outcome.plan.basis, 'FULL_CHANGE');
});

test('invalid and empty diffs fail closed without querying evidence', async () => {
  for (const result of [
    { status: 128, stdout: Buffer.alloc(0), stderr: Buffer.from('missing object') },
    { status: 0, stdout: Buffer.alloc(0), stderr: Buffer.alloc(0) }
  ]) {
    let evidenceCalls = 0;
    const outcome = await buildImpactPlan({
      configuration: configuration(),
      token: 'test-token',
      runGitImpl: () => result,
      verifyWorkflowEvidenceImpl: async () => { evidenceCalls += 1; }
    });
    assert.equal(evidenceCalls, 0);
    assert.equal(outcome.plan.impact, 'full');
    assert.equal(outcome.plan.basis, 'UNKNOWN_OR_INVALID_CHANGE');
  }
});

test('manual runs and architecture-change labels force full execution', async () => {
  let calls = 0;
  const manual = await buildImpactPlan({
    configuration: configuration({
      CI_IMPACT_PROFILE: 'dev',
      CI_IMPACT_EVENT_NAME: 'workflow_dispatch'
    }),
    token: 'test-token',
    runGitImpl: () => { calls += 1; },
    verifyWorkflowEvidenceImpl: async () => { calls += 1; }
  });
  assert.equal(calls, 0);
  assert.equal(manual.plan.basis, 'MANUAL_FULL_RUN');

  const architecture = await buildImpactPlan({
    configuration: configuration({ CI_IMPACT_ARCHITECTURE_CHANGE: 'true' }),
    token: 'test-token',
    runGitImpl: () => gitResult('M', 'docs/FAQ.md'),
    verifyWorkflowEvidenceImpl: async () => { calls += 1; }
  });
  assert.equal(calls, 0);
  assert.equal(architecture.plan.impact, 'documentation');
  assert.equal(architecture.plan.decision, 'full');
  assert.equal(architecture.plan.basis, 'ARCHITECTURE_CHANGE');
});

test('plan command writes one compact output line and escapes summary paths', async () => {
  const stdout = textWriter();
  const stderr = textWriter();
  const writes = new Map();
  const status = await runCiImpactCli({
    argv: ['plan'],
    env: environment({
      GITHUB_OUTPUT: '/tmp/ci-impact-output',
      GITHUB_STEP_SUMMARY: '/tmp/ci-impact-summary'
    }),
    stdout: stdout.stream,
    stderr: stderr.stream,
    appendFileImpl: (path, content) => writes.set(path, (writes.get(path) || '') + content),
    runGitImpl: () => gitResult('M', 'docs/<unsafe>\nname.md'),
    verifyWorkflowEvidenceImpl: async () => successfulEvidence()
  });
  assert.equal(status, 0, stderr.value());
  const output = writes.get('/tmp/ci-impact-output');
  assert.equal(output.split('\n').filter(Boolean).length, 1);
  assert.match(output, /^plan=\{"schemaVersion":1,/);
  const summary = writes.get('/tmp/ci-impact-summary');
  assert.match(summary, /docs\/&lt;unsafe&gt;\\nname\.md/);
  assert.doesNotMatch(summary, /docs\/<unsafe>/);
  assert.match(summary, /Routing: active; pull-request jobs use the trusted-base plan/);
  assert.doesNotMatch(summary, /shadow mode/);
  assert.doesNotMatch(stdout.value(), /test-token/);
});

test('unexpected failures are not converted to full plans and redact tokens', async () => {
  const secret = 'token-with-secret-value';
  const stdout = textWriter();
  const stderr = textWriter();
  const status = await runCiImpactCli({
    argv: ['plan'],
    env: environment({ GITHUB_TOKEN: secret }),
    stdout: stdout.stream,
    stderr: stderr.stream,
    runGitImpl: () => gitResult('M', '.gitignore'),
    verifyWorkflowEvidenceImpl: async () => {
      throw Object.assign(new Error('programming error'), {
        code: 'INJECTED_ERROR',
        details: { diagnostic: secret }
      });
    }
  });
  assert.equal(status, 1);
  assert.equal(stdout.value(), '');
  assert.match(stderr.value(), /INJECTED_ERROR/);
  assert.match(stderr.value(), /\[REDACTED\]/);
  assert.doesNotMatch(stderr.value(), new RegExp(secret));
});

test('CLI rejects unknown arguments and invalid profile/event contracts', async () => {
  const stderr = textWriter();
  const status = await runCiImpactCli({
    argv: ['plan', '--unknown'],
    env: environment(),
    stdout: textWriter().stream,
    stderr: stderr.stream
  });
  assert.equal(status, 1);
  assert.match(stderr.value(), /INVALID_ARGUMENTS/);
  assert.throws(
    () => readPlanConfiguration(environment({ CI_IMPACT_PROFILE: 'gallery' })),
    /PR profile must match/
  );
});

test('workflow routes PR jobs with the trusted base plan and keeps dev routing full', () => {
  const workflow = readFileSync(resolve(REPOSITORY_ROOT, '.github/workflows/test.yml'), 'utf8');
  const workflowJob = (jobId) => workflow.match(
    new RegExp(`\\n  ${jobId}:\\n[\\s\\S]*?(?=\\n  [a-z0-9-]+:\\n|$)`)
  )?.[0];

  assert.match(workflow, /actions: read/);
  assert.doesNotMatch(workflow, /\n    paths(?:-ignore)?:/);

  const planner = workflowJob('ci-impact');
  assert.ok(planner);
  assert.match(planner, /name: CI impact plan/);
  assert.match(planner, /Checkout complete history[\s\S]*fetch-depth: 0/);
  assert.match(planner, /node --test tests\/ci\/\*\.test\.mjs/);
  assert.match(
    planner,
    /ref: \$\{\{ github\.event\.pull_request\.base\.sha \}\}[\s\S]*path: \.ci-trusted-base[\s\S]*persist-credentials: false/
  );
  assert.match(planner, /CI_IMPACT_REPOSITORY_ROOT: \$\{\{ github\.workspace \}\}/);
  assert.match(planner, /run: node \.ci-trusted-base\/tools\/ci-impact\.mjs plan/);
  assert.match(planner, /Build shadow dev CI impact plan[\s\S]*run: node tools\/ci-impact\.mjs plan/);

  for (const jobId of [
    'web-change-budget',
    'core-pr',
    'recipes-standard',
    'gallery',
    'lint',
    'web-pr-smoke'
  ]) {
    const job = workflowJob(jobId);
    assert.ok(job, jobId);
    assert.match(job, /needs: ci-impact/);
    assert.match(job, /needs\.ci-impact\.result == 'success'/);
    assert.match(
      job,
      new RegExp(`contains\\(fromJSON\\(needs\\.ci-impact\\.outputs\\.plan\\)\\.requiredJobs, '${jobId}'\\)`)
    );
    if (['web-change-budget', 'recipes-standard', 'gallery', 'lint'].includes(jobId)) {
      assert.match(job, /!cancelled\(\)/);
      assert.match(job, /github\.event_name == 'push' && github\.ref == 'refs\/heads\/dev'/);
      assert.match(
        job,
        /github\.event_name == 'workflow_dispatch' && github\.ref == 'refs\/heads\/dev'/
      );
    }
  }

  for (const jobId of [
    'core',
    'browser',
    'playwright-functional',
    'playwright-performance',
    'acceptance-supported-main',
    'slow-main',
    'losat-cache-browser-acceptance'
  ]) {
    assert.doesNotMatch(workflowJob(jobId), /requiredJobs/);
  }

  const gate = workflowJob('pr-gate');
  assert.match(gate, /name: PR \/ gate/);
  assert.match(gate, /path: \.ci-trusted-base/);
  assert.match(gate, /CI_IMPACT_PLAN_JSON: \$\{\{ needs\.ci-impact\.outputs\.plan \}\}/);
  assert.match(gate, /CI_IMPACT_NEEDS_JSON: \$\{\{ toJSON\(needs\) \}\}/);
  assert.match(gate, /run: node \.ci-trusted-base\/tools\/ci-impact\.mjs gate/);
  assert.doesNotMatch(gate, /test "\$\{\{ needs\./);
});
