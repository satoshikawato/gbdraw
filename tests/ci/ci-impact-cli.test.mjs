import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';
import {
  GALLERY_PARITY_SPEC,
  IMPACT_PLAN_SCHEMA_VERSION,
  VIBRIO_FULL_GENERATION_SPEC,
  classifyPath,
  createImpactPlan,
  knownJobsFor
} from '../../tools/ci-impact-policy.mjs';
import { PromotionReadinessError } from '../../tools/check-promotion-readiness.mjs';
import {
  buildImpactPlan,
  functionalShardFiles,
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
      return gitResult('M', '.gitignore');
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
  assert.equal(outcome.plan.impact, 'metadata');
  assert.equal(outcome.plan.decision, 'selective');
  assert.equal(outcome.plan.basis, 'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE');
  assert.deepEqual(outcome.plan.requiredJobs, []);
});

test('Product Contract documentation PR selects only Web policy validation', async () => {
  const outcome = await buildImpactPlan({
    configuration: configuration(),
    token: 'test-token',
    runGitImpl: () => gitResult('M', 'docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md'),
    verifyWorkflowEvidenceImpl: async () => assert.fail('Documentation PR must not query staging')
  });
  assert.equal(outcome.plan.impact, 'policy-documentation');
  assert.equal(outcome.plan.decision, 'selective');
  assert.equal(outcome.plan.basis, 'DOCUMENTATION_ONLY_PR');
  assert.equal(outcome.plan.inheritedEvidence, null);
  assert.deepEqual(outcome.plan.requiredJobs, ['web-change-budget']);
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

test('dev metadata and documentation changes select only changed surfaces', async () => {
  for (const [path, impact, requiredJobs] of [
    ['.agents/skills/example/SKILL.md', 'metadata', []],
    ['docs/FAQ.md', 'documentation', ['recipes-standard']],
    ['docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md', 'policy-documentation', ['web-change-budget']]
  ]) {
    let evidenceArguments;
    const outcome = await buildImpactPlan({
      configuration: configuration({
        CI_IMPACT_PROFILE: 'dev',
        CI_IMPACT_EVENT_NAME: 'push'
      }),
      token: 'test-token',
      runGitImpl: () => gitResult('M', path),
      verifyWorkflowEvidenceImpl: async (args) => {
        evidenceArguments = args;
        return successfulEvidence();
      }
    });
    assert.equal(evidenceArguments.expectedHeadSha, SHA.base);
    assert.equal(evidenceArguments.workflowPath, '.github/workflows/test.yml');
    assert.equal(evidenceArguments.expectedAggregateName, 'Dev staging / gate');
    assert.equal(outcome.plan.impact, impact);
    assert.equal(outcome.plan.decision, 'selective');
    assert.equal(outcome.plan.basis, 'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE');
    assert.deepEqual(outcome.plan.requiredJobs, requiredJobs);
  }
});

test('Gallery metadata and documentation changes skip browser and performance', async () => {
  for (const [path, impact] of [
    ['.agents/skills/example/SKILL.md', 'metadata'],
    ['docs/FAQ.md', 'documentation'],
    ['docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md', 'policy-documentation']
  ]) {
    let evidenceArguments;
    const outcome = await buildImpactPlan({
      configuration: configuration({
        CI_IMPACT_PROFILE: 'gallery',
        CI_IMPACT_EVENT_NAME: 'push'
      }),
      token: 'test-token',
      runGitImpl: () => gitResult('M', path),
      verifyWorkflowEvidenceImpl: async (args) => {
        evidenceArguments = args;
        return successfulEvidence({
          workflowPath: '.github/workflows/gallery-publication.yml',
          aggregateName: 'Gallery readiness / gate'
        });
      }
    });
    assert.equal(evidenceArguments.expectedHeadSha, SHA.base);
    assert.equal(evidenceArguments.workflowPath, '.github/workflows/gallery-publication.yml');
    assert.equal(evidenceArguments.expectedAggregateName, 'Gallery readiness / gate');
    assert.equal(outcome.plan.impact, impact);
    assert.equal(outcome.plan.decision, 'selective');
    assert.equal(outcome.plan.basis, 'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE');
    assert.deepEqual(outcome.plan.requiredJobs, []);
  }
});

test('Gallery-impacting surfaces select the full Gallery profile without evidence lookup', async () => {
  for (const path of [
    '.github/workflows/gallery-publication.yml',
    'gbdraw/web/gallery/sessions/example.gbdraw-session.json',
    'gbdraw/web/js/app.js',
    'pyproject.toml',
    'tools/ci-impact.mjs',
    'tests/ci/ci-impact-cli.test.mjs'
  ]) {
    let evidenceCalls = 0;
    const outcome = await buildImpactPlan({
      configuration: configuration({
        CI_IMPACT_PROFILE: 'gallery',
        CI_IMPACT_EVENT_NAME: 'push'
      }),
      token: 'test-token',
      runGitImpl: () => gitResult('M', path),
      verifyWorkflowEvidenceImpl: async () => { evidenceCalls += 1; }
    });
    assert.equal(evidenceCalls, 0, path);
    assert.equal(outcome.plan.impact, classifyPath(path).impact, path);
    assert.equal(outcome.plan.decision, 'full', path);
    assert.equal(outcome.plan.basis, 'FULL_CHANGE', path);
    assert.deepEqual(outcome.plan.requiredJobs, ['browser', 'performance'], path);
  }
});

test('dev control-plane changes run the full profile without inherited evidence', async () => {
  for (const path of ['tools/ci-impact.mjs', '.github/workflows/test.yml']) {
    let evidenceCalls = 0;
    const outcome = await buildImpactPlan({
      configuration: configuration({
        CI_IMPACT_PROFILE: 'dev',
        CI_IMPACT_EVENT_NAME: 'push'
      }),
      token: 'test-token',
      runGitImpl: () => gitResult('M', path),
      verifyWorkflowEvidenceImpl: async () => { evidenceCalls += 1; }
    });
    assert.equal(evidenceCalls, 0, path);
    assert.equal(outcome.plan.impact, classifyPath(path).impact, path);
    assert.equal(outcome.plan.decision, 'full', path);
    assert.equal(outcome.plan.basis, 'FULL_CHANGE', path);
  }
});

test('dev direct-parent staging failures force the current run to full', async () => {
  for (const code of [
    'NO_MATCHING_RUN',
    'RUN_NOT_SUCCESSFUL',
    'AGGREGATE_JOB_NOT_SUCCESSFUL',
    'API_REQUEST_FAILED'
  ]) {
    const outcome = await buildImpactPlan({
      configuration: configuration({
        CI_IMPACT_PROFILE: 'dev',
        CI_IMPACT_EVENT_NAME: 'push'
      }),
      token: 'test-token',
      runGitImpl: () => gitResult('M', '.gitignore'),
      verifyWorkflowEvidenceImpl: async () => {
        throw new PromotionReadinessError(code, `Direct parent evidence unavailable: ${code}`);
      }
    });
    assert.equal(outcome.plan.impact, 'metadata', code);
    assert.equal(outcome.plan.decision, 'full', code);
    assert.equal(outcome.plan.basis, 'INHERITED_EVIDENCE_UNAVAILABLE', code);
    assert.deepEqual(outcome.plan.requiredJobs, [
      'web-change-budget',
      'core',
      'recipes-standard',
      'gallery',
      'browser',
      'playwright-functional',
      'playwright-performance',
      'lint',
      'losat-cache-browser-acceptance'
    ], code);
  }
});

test('documentation-only PRs stay selective without querying any base staging state', async () => {
  for (const code of [
    'NO_MATCHING_RUN', 'RUN_NOT_SUCCESSFUL', 'AGGREGATE_JOB_NOT_SUCCESSFUL', 'API_REQUEST_FAILED'
  ]) {
    for (const [tokens, jobs] of [
      [['M', 'docs/TUTORIALS/1_Intro.md'], ['recipes-standard']],
      [['M', 'docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md'], ['web-change-budget']],
      [['M', '.gitignore', 'M', 'docs/FAQ.md', 'M', 'docs/internal/SELECTIVE_CI.md'],
        ['web-change-budget', 'recipes-standard']],
      [['R100', 'docs/old.md', 'docs/new.md'], ['recipes-standard']],
      [['D', 'docs/FAQ.md'], ['recipes-standard']]
    ]) {
      let calls = 0;
      const outcome = await buildImpactPlan({
        configuration: configuration(),
        token: '',
        runGitImpl: () => gitResult(...tokens),
        verifyWorkflowEvidenceImpl: async () => {
          calls += 1;
          throw new PromotionReadinessError(code, `Base staging unavailable: ${code}`);
        }
      });
      assert.equal(calls, 0);
      assert.equal(outcome.plan.decision, 'selective');
      assert.equal(outcome.plan.basis, 'DOCUMENTATION_ONLY_PR');
      assert.equal(outcome.plan.inheritedEvidence, null);
      assert.equal(outcome.evidenceFailure, undefined);
      assert.deepEqual(outcome.plan.requiredJobs, jobs);
    }
  }
});

test('documentation-only dev pushes without a successful direct parent run the full dev tier', async () => {
  for (const code of [
    'NO_MATCHING_RUN', 'RUN_NOT_SUCCESSFUL', 'AGGREGATE_JOB_NOT_SUCCESSFUL', 'API_REQUEST_FAILED'
  ]) {
    for (const [path, impact] of [
      ['docs/FAQ.md', 'documentation'],
      ['docs/internal/SELECTIVE_CI.md', 'policy-documentation']
    ]) {
      const outcome = await buildImpactPlan({
        configuration: configuration({
          CI_IMPACT_PROFILE: 'dev',
          CI_IMPACT_EVENT_NAME: 'push'
        }),
        token: 'test-token',
        runGitImpl: () => gitResult('M', path),
        verifyWorkflowEvidenceImpl: async () => {
          throw new PromotionReadinessError(code, `Direct parent staging unavailable: ${code}`);
        }
      });
      assert.equal(outcome.plan.impact, impact, code);
      assert.equal(outcome.plan.decision, 'full', code);
      assert.equal(outcome.plan.basis, 'INHERITED_EVIDENCE_UNAVAILABLE', code);
      assert.equal(outcome.plan.inheritedEvidence, null, code);
      assert.deepEqual(outcome.plan.requiredJobs, knownJobsFor('dev'), code);
      assert.equal(outcome.evidenceFailure.code, code);
    }
  }
});

test('documentation-only Gallery pushes fail closed without browser or performance jobs', async () => {
  for (const code of [
    'NO_MATCHING_RUN', 'RUN_NOT_SUCCESSFUL', 'API_REQUEST_FAILED'
  ]) {
    await assert.rejects(buildImpactPlan({
      configuration: configuration({
        CI_IMPACT_PROFILE: 'gallery',
        CI_IMPACT_EVENT_NAME: 'push'
      }),
      token: 'test-token',
      runGitImpl: () => gitResult('M', 'docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md'),
      verifyWorkflowEvidenceImpl: async () => {
        throw new PromotionReadinessError(code, `Direct parent Gallery run failed: ${code}`);
      }
    }), { code: 'DOCUMENTATION_BASE_EVIDENCE_UNAVAILABLE' }, code);
  }
});

test('a successful exact parent aggregate can extend a selective Gallery evidence chain', async () => {
  const outcome = await buildImpactPlan({
    configuration: configuration({
      CI_IMPACT_PROFILE: 'gallery',
      CI_IMPACT_EVENT_NAME: 'push'
    }),
    token: 'test-token',
    runGitImpl: () => gitResult('M', '.gitignore'),
    verifyWorkflowEvidenceImpl: async ({ expectedHeadSha }) => ({
      ...successfulEvidence({
        workflowPath: '.github/workflows/gallery-publication.yml',
        aggregateName: 'Gallery readiness / gate'
      }),
      run: {
        ...successfulEvidence().run,
        headSha: expectedHeadSha,
        composedFromSelectiveEvidence: true
      }
    })
  });
  assert.equal(outcome.plan.decision, 'selective');
  assert.equal(outcome.plan.inheritedEvidence.headSha, SHA.base);
  assert.deepEqual(outcome.plan.requiredJobs, []);
});

test('a zero push before SHA fails closed without querying inherited evidence', async () => {
  for (const profile of ['dev', 'gallery']) {
    let evidenceCalls = 0;
    const outcome = await buildImpactPlan({
      configuration: configuration({
        CI_IMPACT_PROFILE: profile,
        CI_IMPACT_EVENT_NAME: 'push',
        CI_IMPACT_CHANGE_BASE_SHA: '0'.repeat(40)
      }),
      token: 'test-token',
      runGitImpl: () => ({
        status: 128,
        stdout: Buffer.alloc(0),
        stderr: Buffer.from('bad object 0000000000000000000000000000000000000000')
      }),
      verifyWorkflowEvidenceImpl: async () => { evidenceCalls += 1; }
    });
    assert.equal(evidenceCalls, 0, profile);
    assert.equal(outcome.plan.impact, 'full', profile);
    assert.equal(outcome.plan.decision, 'full', profile);
    assert.equal(outcome.plan.basis, 'UNKNOWN_OR_INVALID_CHANGE', profile);
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
    'web-contracts-pr',
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
  assert.equal(outcome.plan.impact, 'ci-only');
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
  for (const [profile, requiredJobs] of [
    ['dev', [
      'web-change-budget',
      'core',
      'recipes-standard',
      'gallery',
      'browser',
      'playwright-functional',
      'playwright-performance',
      'lint',
      'losat-cache-browser-acceptance'
    ]],
    ['gallery', ['browser', 'performance']]
  ]) {
    const manual = await buildImpactPlan({
      configuration: configuration({
        CI_IMPACT_PROFILE: profile,
        CI_IMPACT_EVENT_NAME: 'workflow_dispatch'
      }),
      token: 'test-token',
      runGitImpl: () => { calls += 1; },
      verifyWorkflowEvidenceImpl: async () => { calls += 1; }
    });
    assert.equal(manual.plan.basis, 'MANUAL_FULL_RUN', profile);
    assert.deepEqual(manual.plan.requiredJobs, requiredJobs, profile);
  }
  assert.equal(calls, 0);

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
  assert.match(output, new RegExp(`^plan=\\{"schemaVersion":${IMPACT_PLAN_SCHEMA_VERSION},`));
  const summary = writes.get('/tmp/ci-impact-summary');
  assert.match(summary, /docs\/&lt;unsafe&gt;\\nname\.md/);
  assert.doesNotMatch(summary, /docs\/<unsafe>/);
  assert.match(summary, /Routing: active; pull-request jobs use the trusted-base plan/);
  assert.doesNotMatch(summary, /shadow mode/);
  assert.doesNotMatch(stdout.value(), /test-token/);
});

test('dev plan summary reports active protected-branch routing', async () => {
  const writes = new Map();
  const status = await runCiImpactCli({
    argv: ['plan'],
    env: environment({
      CI_IMPACT_PROFILE: 'dev',
      CI_IMPACT_EVENT_NAME: 'push',
      GITHUB_STEP_SUMMARY: '/tmp/ci-impact-dev-summary'
    }),
    stdout: textWriter().stream,
    stderr: textWriter().stream,
    appendFileImpl: (path, content) => writes.set(path, (writes.get(path) || '') + content),
    runGitImpl: () => gitResult('M', '.gitignore'),
    verifyWorkflowEvidenceImpl: async () => successfulEvidence()
  });
  assert.equal(status, 0);
  const summary = writes.get('/tmp/ci-impact-dev-summary');
  assert.match(summary, /Routing: active; dev staging jobs use the protected-branch plan/);
  assert.doesNotMatch(summary, /shadow mode/);
});

test('Gallery plan summary reports active protected-branch routing', async () => {
  const writes = new Map();
  const status = await runCiImpactCli({
    argv: ['plan'],
    env: environment({
      CI_IMPACT_PROFILE: 'gallery',
      CI_IMPACT_EVENT_NAME: 'push',
      GITHUB_STEP_SUMMARY: '/tmp/ci-impact-gallery-summary'
    }),
    stdout: textWriter().stream,
    stderr: textWriter().stream,
    appendFileImpl: (path, content) => writes.set(path, (writes.get(path) || '') + content),
    runGitImpl: () => gitResult('M', '.gitignore'),
    verifyWorkflowEvidenceImpl: async () => successfulEvidence({
      workflowPath: '.github/workflows/gallery-publication.yml',
      aggregateName: 'Gallery readiness / gate'
    })
  });
  assert.equal(status, 0);
  const summary = writes.get('/tmp/ci-impact-gallery-summary');
  assert.match(
    summary,
    /Routing: active; Gallery readiness jobs use the protected-branch plan/
  );
  assert.doesNotMatch(summary, /shadow mode|observation only/);
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

test('workflow keeps trusted PR routing and activates protected dev routing', () => {
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
  assert.match(planner, /Build dev CI impact plan[\s\S]*run: node tools\/ci-impact\.mjs plan/);

  for (const jobId of [
    'web-change-budget',
    'core-pr',
    'recipes-standard',
    'gallery',
    'lint',
    'web-contracts-pr',
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

  const devJobs = [
    'web-change-budget',
    'core',
    'recipes-standard',
    'gallery',
    'browser',
    'playwright-functional',
    'playwright-performance',
    'lint',
    'losat-cache-browser-acceptance'
  ];
  for (const jobId of devJobs) {
    const job = workflowJob(jobId);
    assert.match(job, /needs: ci-impact/, jobId);
    assert.match(job, /needs\.ci-impact\.result == 'success'/, jobId);
    assert.match(
      job,
      new RegExp(`contains\\(fromJSON\\(needs\\.ci-impact\\.outputs\\.plan\\)\\.requiredJobs, '${jobId}'\\)`),
      jobId
    );
    assert.match(job, /github\.event_name == 'push' && github\.ref == 'refs\/heads\/dev'/, jobId);
    assert.match(
      job,
      /github\.event_name == 'workflow_dispatch' && github\.ref == 'refs\/heads\/dev'/,
      jobId
    );
  }

  const gate = workflowJob('pr-gate');
  assert.match(gate, /name: PR \/ gate/);
  assert.match(gate, /path: \.ci-trusted-base/);
  assert.match(gate, /CI_IMPACT_PLAN_JSON: \$\{\{ needs\.ci-impact\.outputs\.plan \}\}/);
  assert.match(gate, /CI_IMPACT_NEEDS_JSON: \$\{\{ toJSON\(needs\) \}\}/);
  assert.match(gate, /run: node \.ci-trusted-base\/tools\/ci-impact\.mjs gate/);
  assert.doesNotMatch(gate, /test "\$\{\{ needs\./);

  const devGate = workflowJob('dev-staging-gate');
  assert.match(devGate, /name: Dev staging \/ gate/);
  assert.match(devGate, /ref: \$\{\{ github\.sha \}\}/);
  assert.match(devGate, /CI_IMPACT_PLAN_JSON: \$\{\{ needs\.ci-impact\.outputs\.plan \}\}/);
  assert.match(devGate, /CI_IMPACT_NEEDS_JSON: \$\{\{ toJSON\(needs\) \}\}/);
  assert.match(devGate, /CI_IMPACT_EXPECTED_PROFILE: dev/);
  assert.match(devGate, /CI_IMPACT_EXPECTED_WORKFLOW_SHA: \$\{\{ github\.sha \}\}/);
  assert.match(devGate, /run: node tools\/ci-impact\.mjs gate/);
  assert.doesNotMatch(devGate, /test "\$\{\{ needs\./);
});

test('every planned job is a gate dependency that can run for its profile', () => {
  const workflow = readFileSync(resolve(REPOSITORY_ROOT, '.github/workflows/test.yml'), 'utf8');
  const job = (id) => workflow.match(new RegExp(`\\n  ${id}:\\n[\\s\\S]*?(?=\\n  [a-z0-9-]+:\\n|$)`))?.[0] || '';
  const needs = (id) => job(id).match(/\n    needs:\n((?:      - [a-z0-9-]+\n)+)/)?.[1]
    .trim().split('\n').map((line) => line.replace('- ', '').trim()) || [];
  for (const [profile, gate, event] of [
    ['pr', 'pr-gate', /github\.event_name == 'pull_request' && github\.base_ref == 'dev'/],
    ['dev', 'dev-staging-gate', /github\.event_name == 'push' && github\.ref == 'refs\/heads\/dev'/]
  ]) {
    for (const jobId of knownJobsFor(profile)) {
      assert.ok(needs(gate).includes(jobId), `${gate} must need ${jobId}`);
      assert.match(job(jobId), new RegExp(`requiredJobs, '${jobId}'`), jobId);
      assert.match(job(jobId), event, `${jobId} must run for the ${profile} profile`);
    }
  }
});

test('web PR route inherits only exact base evidence and falls back to full on API failure', async () => {
  for (const available of [true, false]) {
    const outcome = await buildImpactPlan({
      configuration: configuration(), token: 'test-token',
      runGitImpl: () => gitResult('M', 'gbdraw/web/js/app/label-editor.js'),
      verifyWorkflowEvidenceImpl: async ({ expectedHeadSha }) => {
        assert.equal(expectedHeadSha, SHA.base);
        if (!available) throw new PromotionReadinessError('API_REQUEST_FAILED', 'unavailable');
        return successfulEvidence();
      }
    });
    assert.equal(outcome.plan.decision, available ? 'selective' : 'full');
    assert.equal(outcome.plan.requiredJobs.includes('core-pr'), !available);
    assert.ok(outcome.plan.requiredJobs.includes('web-contracts-pr'));
    assert.ok(outcome.plan.requiredJobs.includes('web-pr-smoke'));
    assert.ok(outcome.plan.requiredJobs.includes('gallery'));
    assert.ok(outcome.plan.requiredJobs.includes('playwright-functional'));
  }
});

test('release dispatch is exhaustive and cannot be inferred from a routine push', async () => {
  const outcome = await buildImpactPlan({
    configuration: configuration({ CI_IMPACT_PROFILE: 'release', CI_IMPACT_EVENT_NAME: 'workflow_dispatch' }),
    runGitImpl: () => { throw new Error('manual release must not classify a diff'); },
    verifyWorkflowEvidenceImpl: () => { throw new Error('manual release cannot inherit evidence'); }
  });
  assert.equal(outcome.plan.profile, 'release');
  assert.equal(outcome.plan.basis, 'MANUAL_FULL_RUN');
  assert.ok(outcome.plan.requiredJobs.includes('acceptance-supported-main'));
  assert.ok(outcome.plan.requiredJobs.includes('slow-main'));
  assert.ok(outcome.plan.requiredJobs.includes('vibrio-generate-release'));
  assert.throws(() => configuration({ CI_IMPACT_PROFILE: 'release', CI_IMPACT_EVENT_NAME: 'push' }), /explicit dispatch/);
});

test('release workflow binds exhaustive matrices and package/browser contracts to its gate', () => {
  const workflow = readFileSync(resolve(REPOSITORY_ROOT, '.github/workflows/test.yml'), 'utf8');
  const job = (id) => workflow.match(new RegExp(`\\n  ${id}:\\n[\\s\\S]*?(?=\\n  [a-z0-9-]+:\\n|$)`))?.[0] || '';
  const release = job('release-gate');
  assert.match(release, /CI_IMPACT_EXPECTED_PROFILE: release/);
  assert.match(release, /inputs\.tier == 'release'/);
  for (const id of ['core', 'recipes-standard', 'gallery', 'browser', 'playwright-functional', 'playwright-performance', 'losat-cache-browser-acceptance', 'acceptance-supported-main', 'slow-main', 'vibrio-generate-release']) {
    assert.ok(release.includes(`      - ${id}\n`), `release gate missing ${id}`);
  }
  assert.match(job('core'), /python-version: \["3.10", "3.11", "3.12"\]/);
  assert.match(job('acceptance-supported-main'), /python-version: \["3.10", "3.12"\]/);
  assert.match(job('acceptance-supported-main'), /surface: \["recipe", "gallery", "browser"\]/);
  assert.match(job('slow-main'), /python-version: \["3.10", "3.11", "3.12"\]/);
  assert.match(job('slow-main'), /-m "slow and not browser"/);
  assert.match(job('vibrio-generate-release'), /run: npm run test:web:vibrio-generate\n/);
  assert.match(job('vibrio-generate-release'), /requiredJobs, 'vibrio-generate-release'/);
  assert.match(job('browser'), /Run package build integration[\s\S]*-m "slow and not browser"/);
  assert.match(job('browser'), /Run offline GUI browser contracts[\s\S]*-m "slow and browser"/);
  assert.match(job('pr-gate'), /sparse-checkout: tools[\s\S]*node \.ci-trusted-base\/tools\/ci-impact.mjs gate/);
});

test('Gallery alone owns PR parity without expanding dev or release execution', () => {
  const workflow = readFileSync(resolve(REPOSITORY_ROOT, '.github/workflows/test.yml'), 'utf8');
  const job = (id) => workflow.match(new RegExp(`\\n  ${id}:\\n[\\s\\S]*?(?=\\n  [a-z0-9-]+:\\n|$)`))?.[0] || '';
  const gallery = job('gallery');
  const smoke = job('web-pr-smoke');
  assert.equal((workflow.match(/run: npm run test:web:gallery-publication/g) || []).length, 1);
  assert.match(gallery, /\n    timeout-minutes: 25\n/);
  assert.match(gallery, /GBDRAW_GALLERY_PR_PARITY: \$\{\{ github\.event_name == 'pull_request' && github\.base_ref == 'dev' && contains\(fromJSON\(needs\.ci-impact\.outputs\.plan\)\.requiredJobs, 'web-pr-smoke'\) \}\}/);
  assert.match(gallery, /python -m pytest tests\/[\s\S]*-m "gallery and not slow"/);
  for (const name of ['Set up Node.js for Gallery parity', 'Install Gallery parity dependencies', 'Prepare Gallery browser wheel', 'Verify Gallery first-Generate parity']) {
    assert.match(gallery, new RegExp(`name: ${name}\\n        if: env\\.GBDRAW_GALLERY_PR_PARITY == 'true'`));
  }
  assert.match(gallery, /node-version: "20"/);
  assert.match(gallery, /npm ci/);
  assert.match(gallery, /npx playwright install --with-deps chromium/);
  assert.equal((gallery.match(/run: python tools\/prepare_browser_wheel\.py/g) || []).length, 1);
  assert.match(gallery, /if: failure\(\) && env\.GBDRAW_GALLERY_PR_PARITY == 'true'/);
  assert.match(gallery, /path: test-results\//);
  assert.match(smoke, /\n    timeout-minutes: 20\n/);
  assert.match(smoke, /run: npm run test:web:pr-smoke/);
  assert.doesNotMatch(smoke, /test:web:gallery-publication/);
  assert.equal((smoke.match(/run: python tools\/prepare_browser_wheel\.py/g) || []).length, 1);
  assert.match(smoke, /Upload Playwright PR smoke traces/);
});

test('browser jobs seed apt from one verified cache and bound their test steps', () => {
  const workflow = readFileSync(resolve(REPOSITORY_ROOT, '.github/workflows/test.yml'), 'utf8');
  const job = (id) => workflow.match(new RegExp(`\\n  ${id}:\\n[\\s\\S]*?(?=\\n  [a-z0-9-]+:\\n|$)`))?.[0] || '';
  const step = (source, name) => source.match(
    new RegExp(`\\n      - name: ${name}\\n[\\s\\S]*?(?=\\n      - |$)`)
  )?.[0] || '';
  const jobIds = [...workflow.slice(workflow.indexOf('\njobs:\n')).matchAll(/\n  ([a-z0-9-]+):\n/g)]
    .map(([, id]) => id);
  const browserJobs = jobIds.filter((id) => job(id).includes('playwright install --with-deps chromium'));
  assert.deepEqual(browserJobs, [
    'gallery',
    'browser',
    'web-contracts-pr',
    'web-pr-smoke',
    'playwright-functional',
    'playwright-performance',
    'acceptance-supported-main',
    'vibrio-generate-release',
    'losat-cache-browser-acceptance',
    'promotion-audit-sweeps',
    'promotion-audit-journeys',
    'promotion-audit-random-walk'
  ]);

  const cacheSteps = [
    'Restore Playwright system packages',
    'Install Playwright Chromium and system packages',
    'Save Playwright system packages'
  ];
  const prefix = 'playwright-apt-v1-${{ runner.os }}-${{ runner.arch }}-';
  const copies = new Set();
  for (const id of browserJobs) {
    const source = job(id);
    const [restore, install, save] = cacheSteps.map((name) => step(source, name));
    assert.ok(restore && install && save, `${id} must restore, install, and save`);
    assert.ok(source.indexOf(restore) < source.indexOf(install), id);
    assert.ok(source.indexOf(install) < source.indexOf(save), id);
    assert.equal((source.match(/playwright install --with-deps chromium/g) || []).length, 1, id);
    assert.match(restore, /id: playwright-apt\n[\s\S]*uses: actions\/cache\/restore@v4/);
    assert.ok(restore.includes(`key: ${prefix}\n`), id);
    assert.ok(restore.includes(`restore-keys: ${prefix}`), id);
    // apt only receives cached files; the install command and its exit status are unchanged.
    assert.match(
      install,
      /id: playwright-install\n[\s\S]*sudo cp -t \/var\/cache\/apt\/archives\/ ~\/\.cache\/playwright-apt\/\*\.deb [^\n]*\n {10}(?:npx|python -m) playwright install --with-deps chromium\n/
    );
    assert.doesNotMatch(install, /playwright install[^\n]*(?:\|\||continue-on-error)/);
    assert.doesNotMatch(source, /continue-on-error/);
    assert.match(install, /sudo apt-get autoclean/);
    assert.match(install, /sha256sum -- \*\.deb \| LC_ALL=C sort \| sha256sum/);
    assert.match(
      install,
      /echo "apt-cache-key=playwright-apt-v1-\$\{RUNNER_OS\}-\$\{RUNNER_ARCH\}-\$\{sum\}" >> "\$GITHUB_OUTPUT"/
    );
    assert.match(save, /uses: actions\/cache\/save@v4/);
    assert.match(save, /steps\.playwright-install\.outputs\.apt-cache-key != ''/);
    assert.match(
      save,
      /steps\.playwright-install\.outputs\.apt-cache-key != steps\.playwright-apt\.outputs\.cache-matched-key/
    );
    assert.match(save, /key: \$\{\{ steps\.playwright-install\.outputs\.apt-cache-key \}\}/);
    copies.add([restore, install, save].join('')
      .replace(/\n {8}if: (?:env\.GBDRAW_GALLERY_PR_PARITY == 'true'|matrix\.surface == 'browser'|steps\.shard\.outputs\.files != '')(?=\n)/g, '')
      .replace('python -m playwright', 'npx playwright'));
  }
  assert.equal(copies.size, 1, 'every browser job must use the same cache steps');

  const jobTimeouts = Object.fromEntries(browserJobs.map((id) => [
    id,
    Number(job(id).match(/\n    timeout-minutes: (\d+)\n/)?.[1])
  ]));
  assert.deepEqual(jobTimeouts, {
    gallery: 25,
    browser: 20,
    'web-contracts-pr': 20,
    'web-pr-smoke': 20,
    'playwright-functional': 45,
    'playwright-performance': 25,
    'acceptance-supported-main': 20,
    'vibrio-generate-release': 35,
    'losat-cache-browser-acceptance': 20,
    'promotion-audit-sweeps': 120,
    'promotion-audit-journeys': 45,
    'promotion-audit-random-walk': 120
  });
  const stepTimeouts = {
    gallery: { 'Run Gallery tests': 10, 'Verify Gallery first-Generate parity': 10 },
    browser: {
      'Run Web JavaScript tests': 5,
      'Run Python browser tests': 10,
      'Run package build integration': 5,
      'Run offline GUI browser contracts': 5
    },
    'web-contracts-pr': { 'Run fast Web JavaScript contracts': 5, 'Run non-slow Python browser tests': 10 },
    'web-pr-smoke': { 'Run Playwright PR smoke': 10 },
    'promotion-audit-journeys': { 'Run the user journey': 35 },
    'promotion-audit-random-walk': { 'Run the live-vs-Generate random walk': 100 },
    'playwright-performance': { 'Run Playwright performance tests': 10 },
    'vibrio-generate-release': { 'Run Vibrio full generation': 25 }
  };
  for (const [id, limits] of Object.entries(stepTimeouts)) {
    for (const [name, minutes] of Object.entries(limits)) {
      assert.match(step(job(id), name), new RegExp(`\\n {8}timeout-minutes: ${minutes}\\n`), `${id}: ${name}`);
    }
  }

  for (const id of ['ci-impact', 'web-change-budget']) {
    assert.match(
      job(id),
      /uses: actions\/checkout@v4\n {8}with:\n {10}fetch-depth: 0\n(?: {10}#[^\n]*\n)? {10}filter: blob:none\n/,
      id
    );
  }
  assert.match(job('ci-impact'), /node-version: "20"\n {10}cache: npm\n[\s\S]*run: npm ci/);
});

test('aggregate gates and the Gallery planner tolerate a slow checkout', () => {
  const read = (path) => readFileSync(resolve(REPOSITORY_ROOT, path), 'utf8');
  const jobIn = (workflow, id) => workflow.match(
    new RegExp(`\\n  ${id}:\\n[\\s\\S]*?(?=\\n  [a-z0-9-]+:\\n|$)`)
  )?.[0] || '';
  const tests = read('.github/workflows/test.yml');
  const gallery = read('.github/workflows/gallery-publication.yml');
  // A checkout alone has taken 63 s, so a 1-minute gate can cancel a run whose jobs all passed.
  for (const [workflow, id] of [
    [tests, 'pr-gate'],
    [tests, 'dev-staging-gate'],
    [tests, 'release-gate'],
    [gallery, 'readiness-gate']
  ]) {
    assert.match(jobIn(workflow, id), /\n    timeout-minutes: 5\n/, id);
  }
  // A full-history checkout with blobs has taken 4 min 53 s of this 5-minute job.
  const planner = jobIn(gallery, 'ci-impact');
  assert.match(planner, /\n    timeout-minutes: 5\n/);
  assert.match(
    planner,
    /uses: actions\/checkout@v4\n {8}with:\n {10}ref: \$\{\{ github\.sha \}\}\n {10}fetch-depth: 0\n(?: {10}#[^\n]*\n)? {10}filter: blob:none\n/
  );
});


const TREE = Object.freeze({ same: 'd'.repeat(40), other: 'e'.repeat(40) });
const TEST_SHARDS = Object.freeze([
  ['tests/web/a.playwright.spec.js', 'tests/web/right-drawer.playwright.spec.js'],
  ['tests/web/contracts/session-regenerate-intent.playwright.spec.js'],
  [], [], [], [], [], []
]);

const emptyResult = (status = 0) => ({ status, stdout: Buffer.alloc(0), stderr: Buffer.alloc(0) });

const gitRouter = ({
  diff = null,
  trees = [TREE.same, TREE.other],
  grep = [],
  grepStatus,
  blobs = {},
  ancestor = true,
  calls = []
} = {}) => (_root, args) => {
  calls.push(args);
  if (args[0] === 'diff') return diff === null ? emptyResult() : gitResult(...diff);
  if (args[0] === 'rev-parse' && args[1] === '--verify') {
    return { status: 0, stdout: Buffer.from(`${args.at(-1).replace('^{commit}', '')}\n`), stderr: Buffer.alloc(0) };
  }
  if (args[0] === 'rev-parse') {
    return trees === null ? emptyResult(128) : { status: 0, stdout: Buffer.from(`${trees.join('\n')}\n`), stderr: Buffer.alloc(0) };
  }
  if (args[0] === 'merge-base') return emptyResult(ancestor ? 0 : 1);
  if (args[0] === 'grep') {
    const head = args[args.indexOf('--') - 1];
    const output = grep.map(([path, line]) => `${head}:${path}\0${line}\n`).join('');
    return { status: grepStatus ?? (output ? 0 : 1), stdout: Buffer.from(output), stderr: Buffer.alloc(0) };
  }
  if (args[0] === 'cat-file') {
    const spec = args.at(-1);
    const path = spec.slice(spec.indexOf(':') + 1);
    const source = { 'tests/ci/functional-shards.json': JSON.stringify({ minutes: {}, shards: TEST_SHARDS }), ...blobs }[path];
    return source === undefined
      ? { status: 128, stdout: Buffer.alloc(0), stderr: Buffer.from(`missing ${path}`) }
      : { status: 0, stdout: Buffer.from(source), stderr: Buffer.alloc(0) };
  }
  throw new Error(`unexpected git ${args.join(' ')}`);
};

const pushConfiguration = (profile) => configuration({ CI_IMPACT_PROFILE: profile, CI_IMPACT_EVENT_NAME: 'push' });
const parentEvidence = (profile) => successfulEvidence(profile === 'gallery' ? {
  workflowPath: '.github/workflows/gallery-publication.yml',
  aggregateName: 'Gallery readiness / gate'
} : {});

test('identical trees inherit every job on dev and Gallery pushes with direct-parent evidence', async () => {
  for (const profile of ['dev', 'gallery']) {
    const calls = [];
    let evidenceArguments;
    const outcome = await buildImpactPlan({
      configuration: pushConfiguration(profile),
      token: 'test-token',
      runGitImpl: gitRouter({ diff: null, trees: [TREE.same, TREE.same], calls }),
      verifyWorkflowEvidenceImpl: async (args) => {
        evidenceArguments = args;
        return parentEvidence(profile);
      }
    });
    assert.deepEqual(calls.find((args) => args[0] === 'rev-parse'),
      ['rev-parse', `${SHA.base}^{tree}`, `${SHA.head}^{tree}`], profile);
    assert.equal(evidenceArguments.expectedHeadSha, SHA.base, profile);
    assert.equal(outcome.plan.decision, 'selective', profile);
    assert.equal(outcome.plan.basis, 'IDENTICAL_TREE_WITH_DIRECT_PARENT_EVIDENCE', profile);
    assert.deepEqual(outcome.plan.capabilities, ['none'], profile);
    assert.equal(outcome.plan.changedPathCount, 0, profile);
    assert.deepEqual(outcome.plan.requiredJobs, [], profile);
    assert.equal(outcome.plan.inheritedEvidence.headSha, SHA.base, profile);
  }
});

test('identical trees without parent evidence run the complete profile', async () => {
  for (const profile of ['dev', 'gallery']) {
    const outcome = await buildImpactPlan({
      configuration: pushConfiguration(profile),
      token: 'test-token',
      runGitImpl: gitRouter({ diff: null, trees: [TREE.same, TREE.same] }),
      verifyWorkflowEvidenceImpl: async () => {
        throw new PromotionReadinessError('RUN_NOT_SUCCESSFUL', 'parent failed');
      }
    });
    assert.equal(outcome.plan.decision, 'full', profile);
    assert.equal(outcome.plan.basis, 'INHERITED_EVIDENCE_UNAVAILABLE', profile);
    assert.deepEqual(outcome.plan.requiredJobs, knownJobsFor(profile), profile);
  }
});

test('an empty diff without a proven identical tree stays invalid', async () => {
  for (const trees of [[TREE.same, TREE.other], null]) {
    let evidenceCalls = 0;
    const outcome = await buildImpactPlan({
      configuration: pushConfiguration('dev'),
      token: 'test-token',
      runGitImpl: gitRouter({ diff: null, trees }),
      verifyWorkflowEvidenceImpl: async () => { evidenceCalls += 1; }
    });
    assert.equal(evidenceCalls, 0);
    assert.equal(outcome.plan.impact, 'full');
    assert.equal(outcome.plan.basis, 'UNKNOWN_OR_INVALID_CHANGE');
  }
});

const leafOutcome = (profile, gitOptions, evidence = async () => parentEvidence(profile)) => buildImpactPlan({
  configuration: pushConfiguration(profile),
  token: 'test-token',
  runGitImpl: gitRouter(gitOptions),
  verifyWorkflowEvidenceImpl: evidence
});

test('dev leaf-test pushes run only the jobs that execute the changed tests', async () => {
  const allPytest = ['core', 'recipes-standard', 'gallery', 'browser'];
  for (const [diff, blobs, requiredJobs, leafTests] of [
    [['M', 'tests/web/label-editor.test.mjs'], {}, ['browser'],
      [{ path: 'tests/web/label-editor.test.mjs', kind: 'node', jobs: ['browser'] }]],
    [['M', 'tests/web/right-drawer.playwright.spec.js'], {}, ['playwright-functional'],
      [{ path: 'tests/web/right-drawer.playwright.spec.js', kind: 'functional', jobs: ['playwright-functional'] }]],
    [['M', 'tests/web/webapp-performance.playwright.spec.js'], {}, ['playwright-performance'],
      [{ path: 'tests/web/webapp-performance.playwright.spec.js', kind: 'performance', jobs: ['playwright-performance'] }]],
    [['M', GALLERY_PARITY_SPEC], {}, [], [{ path: GALLERY_PARITY_SPEC, kind: 'gallery-parity', jobs: [] }]],
    [['M', VIBRIO_FULL_GENERATION_SPEC], {}, [], [{ path: VIBRIO_FULL_GENERATION_SPEC, kind: 'release-only', jobs: [] }]],
    [['M', 'tests/test_regression.py'], { 'tests/test_regression.py': 'def test_a():\n    assert True\n' }, ['core'],
      [{ path: 'tests/test_regression.py', kind: 'python', jobs: ['core'] }]],
    [['A', 'tests/test_new_browser.py'], { 'tests/test_new_browser.py': 'import pytest\npytestmark = pytest.mark.browser\n' },
      ['core', 'browser'], [{ path: 'tests/test_new_browser.py', kind: 'python', jobs: ['core', 'browser'] }]],
    [['M', 'tests/test_dynamic.py'], { 'tests/test_dynamic.py': 'import pytest\nmarker = getattr(pytest.mark, "browser")\n' },
      allPytest, [{ path: 'tests/test_dynamic.py', kind: 'python', jobs: allPytest }]],
    [['D', 'tests/test_retired.py'], {}, allPytest, [{ path: 'tests/test_retired.py', kind: 'python', jobs: allPytest }]],
    [['M', 'docs/FAQ.md', 'M', 'tests/web/label-editor.test.mjs'], {}, ['recipes-standard', 'browser'],
      [{ path: 'tests/web/label-editor.test.mjs', kind: 'node', jobs: ['browser'] }]],
    [['M', 'docs/internal/SELECTIVE_CI.md', 'M', '.gitignore', 'D', 'tests/web/retired.test.mjs'], {}, ['web-change-budget', 'browser'],
      [{ path: 'tests/web/retired.test.mjs', kind: 'node', jobs: ['browser'] }]]
  ]) {
    const outcome = await leafOutcome('dev', { diff, blobs });
    assert.equal(outcome.plan.decision, 'selective', diff.join(' '));
    assert.equal(outcome.plan.basis, 'LEAF_TEST_CHANGE_WITH_DIRECT_PARENT_EVIDENCE', diff.join(' '));
    assert.deepEqual(outcome.plan.requiredJobs, requiredJobs, diff.join(' '));
    assert.deepEqual(outcome.plan.leafTests, leafTests, diff.join(' '));
  }
});

test('the reference check names the head tree and skips documentation and runner lists', async () => {
  const calls = [];
  const outcome = await leafOutcome('dev', {
    diff: ['M', 'tests/web/request.test.mjs', 'M', 'tests/test_regression.py'],
    blobs: { 'tests/test_regression.py': 'def test_a():\n    pass\n' },
    calls,
    grep: [
      ['tests/web/request.test.mjs', "import './request.test.mjs';"],
      ['tests/test_api_session.py', 'NODE_TEST = "tests/web/session-request.test.mjs"'],
      ['tests/test_other.py', 'def test_regression_handles_inputs():'],
      ['README.md', 'Run tests/web/request.test.mjs and tests/test_regression.py'],
      ['docs/FAQ.md', 'tests/test_regression.py'],
      ['.github/pull_request_template.md', 'tests/web/request.test.mjs'],
      ['.github/workflows/test.yml', 'node tests/web/request.test.mjs'],
      ['tests/ci/functional-shards.json', '"tests/web/request.test.mjs"'],
      ['playwright.config.js', 'request.test.mjs'],
      ['package.json', '"x": "node tests/web/request.test.mjs"']
    ]
  });
  const grep = calls.find((args) => args[0] === 'grep');
  assert.deepEqual(grep.slice(0, 4), ['grep', '-I', '-z', '-F']);
  assert.ok(grep.includes('request.test.mjs'));
  assert.ok(grep.includes('test_regression'));
  assert.equal(grep[grep.indexOf('--') - 1], SHA.head);
  for (const exclusion of [':(exclude)docs/', ':(exclude)tests/ci/', ':(exclude).github/workflows/']) {
    assert.ok(grep.includes(exclusion), exclusion);
  }
  assert.equal(outcome.plan.basis, 'LEAF_TEST_CHANGE_WITH_DIRECT_PARENT_EVIDENCE');
  assert.deepEqual(outcome.plan.requiredJobs, ['core', 'browser']);
});

test('test files that another file uses, shared test code, and unrouted specs keep the complete dev tier', async () => {
  for (const [label, options] of [
    ['run by another test', {
      diff: ['M', 'tests/web/session-request.test.mjs'],
      grep: [['tests/test_api_session.py', 'NODE_TEST = REPOSITORY / "tests/web/session-request.test.mjs"']]
    }],
    ['imported Python module', {
      diff: ['M', 'tests/test_collinearity_units.py'],
      blobs: { 'tests/test_collinearity_units.py': 'def _observe():\n    pass\n' },
      grep: [['tests/test_protein_s078.py', 'from tests.test_collinearity_units import _boundaries']]
    }],
    ['named by a tool', {
      diff: ['M', 'tests/web/feature-catalog.test.mjs'],
      grep: [['tools/web-product-impact-map.json', '"tests/web/feature-catalog.test.mjs"']]
    }],
    ['shared helper', { diff: ['M', 'tests/web/helpers/app-lifecycle.cjs'] }],
    ['shared input', { diff: ['M', 'tests/test_inputs/example.gbk'] }],
    ['spec missing from the shard map', { diff: ['A', 'tests/web/new-feature.playwright.spec.js'] }],
    ['renamed spec without a shard map change', {
      diff: ['R100', 'tests/web/right-drawer.playwright.spec.js', 'tests/web/right-drawer-renamed.playwright.spec.js']
    }],
    ['unrouted contract spec', { diff: ['M', 'tests/web/contracts/unknown.serial.spec.js'] }],
    ['reference scan failure', { diff: ['M', 'tests/web/label-editor.test.mjs'], grepStatus: 2 }],
    ['leaf test with runtime source', {
      diff: ['M', 'tests/web/label-editor.test.mjs', 'M', 'gbdraw/web/js/app/label-editor.js']
    }]
  ]) {
    let evidenceCalls = 0;
    const outcome = await leafOutcome('dev', options, async () => {
      evidenceCalls += 1;
      return parentEvidence('dev');
    });
    assert.equal(evidenceCalls, 0, label);
    assert.equal(outcome.plan.decision, 'full', label);
    assert.equal(outcome.plan.basis, 'FULL_CHANGE', label);
    assert.equal(outcome.plan.leafTests, null, label);
    assert.deepEqual(outcome.plan.requiredJobs, knownJobsFor('dev'), label);
  }
});

test('leaf-test pushes without parent evidence run the complete profile', async () => {
  for (const profile of ['dev', 'gallery']) {
    const outcome = await leafOutcome(profile, { diff: ['M', 'tests/web/label-editor.test.mjs'] }, async () => {
      throw new PromotionReadinessError('NO_MATCHING_RUN', 'parent run replaced');
    });
    assert.equal(outcome.plan.decision, 'full', profile);
    assert.equal(outcome.plan.basis, 'INHERITED_EVIDENCE_UNAVAILABLE', profile);
    assert.equal(outcome.plan.leafTests, null, profile);
    assert.deepEqual(outcome.plan.requiredJobs, knownJobsFor(profile), profile);
  }
});

test('Gallery publication inherits leaf-test changes unless its browser job runs the changed spec', async () => {
  const other = await leafOutcome('gallery', { diff: ['M', 'tests/web/gallery-session-publication.test.mjs'] });
  assert.equal(other.plan.basis, 'LEAF_TEST_CHANGE_WITH_DIRECT_PARENT_EVIDENCE');
  assert.deepEqual(other.plan.requiredJobs, []);
  const parity = await leafOutcome('gallery', { diff: ['M', 'docs/FAQ.md', 'M', GALLERY_PARITY_SPEC] });
  assert.equal(parity.plan.basis, 'LEAF_TEST_CHANGE_WITH_DIRECT_PARENT_EVIDENCE');
  assert.deepEqual(parity.plan.requiredJobs, ['browser']);
  assert.deepEqual(parity.plan.leafTests, [{ path: GALLERY_PARITY_SPEC, kind: 'gallery-parity', jobs: ['browser'] }]);
});

test('pull request planning ignores the leaf-test route', async () => {
  const calls = [];
  const outcome = await buildImpactPlan({
    configuration: configuration(),
    token: 'test-token',
    runGitImpl: gitRouter({ diff: ['M', 'tests/web/label-editor.test.mjs'], calls }),
    verifyWorkflowEvidenceImpl: async () => successfulEvidence()
  });
  assert.equal(calls.some((args) => args[0] === 'grep'), false);
  assert.equal(outcome.plan.basis, 'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE');
  assert.equal(outcome.plan.leafTests, null);
  assert.ok(outcome.plan.requiredJobs.includes('playwright-functional'));
});

const devLeafPlan = (leafTests) => createImpactPlan({
  profile: 'dev',
  impact: 'web-runtime',
  capabilities: ['web-runtime'],
  decision: 'selective',
  basis: 'LEAF_TEST_CHANGE_WITH_DIRECT_PARENT_EVIDENCE',
  changeBaseSha: SHA.base,
  changeHeadSha: SHA.head,
  workflowSha: SHA.workflow,
  changedPathCount: leafTests.length,
  inheritedEvidence: {
    workflowPath: '.github/workflows/test.yml',
    aggregateName: 'Dev staging / gate',
    headSha: SHA.base,
    runId: 101,
    aggregateJobId: 202,
    runUrl: 'https://github.com/satoshikawato/gbdraw/actions/runs/101',
    aggregateJobUrl: 'https://github.com/satoshikawato/gbdraw/actions/runs/101/job/202'
  },
  leafTests
});

const fullDevPlan = () => createImpactPlan({
  profile: 'dev',
  impact: 'web-runtime',
  decision: 'full',
  basis: 'FULL_CHANGE',
  changeBaseSha: SHA.base,
  changeHeadSha: SHA.head,
  workflowSha: SHA.workflow,
  changedPathCount: 1,
  inheritedEvidence: null
});

test('functional shards intersect their assigned specs with the changed leaf specs', () => {
  const plan = devLeafPlan([
    { path: 'tests/web/right-drawer.playwright.spec.js', kind: 'functional', jobs: ['playwright-functional'] }
  ]);
  assert.deepEqual(functionalShardFiles({ plan, shards: TEST_SHARDS, shard: 1 }), ['tests/web/right-drawer.playwright.spec.js']);
  assert.deepEqual(functionalShardFiles({ plan, shards: TEST_SHARDS, shard: 2 }), []);
  assert.deepEqual(functionalShardFiles({ plan: fullDevPlan(), shards: TEST_SHARDS, shard: 1 }), TEST_SHARDS[0]);
  assert.deepEqual(functionalShardFiles({ plan: fullDevPlan(), shards: TEST_SHARDS, shard: 2 }), TEST_SHARDS[1]);
  for (const shard of [0, 9, 1.5]) {
    assert.throws(() => functionalShardFiles({ plan, shards: TEST_SHARDS, shard }), { code: 'INVALID_SHARD' });
  }
  // A trusted base plan of another schema runs the whole shard; a malformed current plan fails.
  const olderSchema = { ...JSON.parse(JSON.stringify(plan)), schemaVersion: IMPACT_PLAN_SCHEMA_VERSION - 1 };
  delete olderSchema.leafTests;
  assert.deepEqual(functionalShardFiles({ plan: olderSchema, shards: TEST_SHARDS, shard: 1 }), TEST_SHARDS[0]);
  assert.throws(() => functionalShardFiles({ plan: { ...plan, leafTests: [] }, shards: TEST_SHARDS, shard: 1 }),
    { code: 'INVALID_LEAF_TESTS' });
  const browserOnly = devLeafPlan([{ path: 'tests/web/label-editor.test.mjs', kind: 'node', jobs: ['browser'] }]);
  assert.throws(() => functionalShardFiles({ plan: browserOnly, shards: TEST_SHARDS, shard: 1 }),
    { code: 'FUNCTIONAL_NOT_REQUIRED' });
});

test('shard-files writes the changed specs of one shard from the real shard map', async () => {
  const shards = JSON.parse(readFileSync(resolve(REPOSITORY_ROOT, 'tests/ci/functional-shards.json'), 'utf8')).shards;
  // devLeafPlan declares the web-runtime capability, so the leaf spec must be
  // one of shard 1's web-runtime specs; the shard order changes on rebalancing.
  const spec = shards[0].find((path) => classifyPath(path).impact === 'web-runtime');
  assert.ok(spec, 'shard 1 holds a web-runtime functional spec');
  const plan = devLeafPlan([{ path: spec, kind: 'functional', jobs: ['playwright-functional'] }]);
  for (const shard of [1, 2]) {
    const writes = new Map();
    const stdout = textWriter();
    const stderr = textWriter();
    const status = await runCiImpactCli({
      argv: ['shard-files', String(shard)],
      env: { CI_IMPACT_PLAN_JSON: JSON.stringify(plan), GITHUB_OUTPUT: '/tmp/shard-output' },
      cwd: REPOSITORY_ROOT,
      stdout: stdout.stream,
      stderr: stderr.stream,
      appendFileImpl: (path, content) => writes.set(path, (writes.get(path) || '') + content)
    });
    assert.equal(status, 0, stderr.value());
    const files = shard === 1 ? spec : '';
    assert.equal(writes.get('/tmp/shard-output'), `files=${files}\n`);
  }
});

const classifyOutcome = async (gitOptions, argv = ['classify', '--base', SHA.base, '--head', SHA.head]) => {
  const stdout = textWriter();
  const stderr = textWriter();
  const status = await runCiImpactCli({
    argv,
    env: {},
    cwd: REPOSITORY_ROOT,
    stdout: stdout.stream,
    stderr: stderr.stream,
    runGitImpl: gitRouter(gitOptions),
    verifyWorkflowEvidenceImpl: async () => assert.fail('classify must not query workflow evidence')
  });
  return { status, stderr: stderr.value(), result: status === 0 ? JSON.parse(stdout.value()) : null };
};

test('classify reports carry-forward verdicts from local Git data only', async () => {
  const all = { releaseEvidenceCarries: true, generatedArtifactChecksCarry: true, localTestEvidenceCarries: true };
  const sessionLog = await classifyOutcome({ diff: ['M', 'docs/internal/web-gui-audit-20260930/SESSION_LOG.md'] });
  assert.equal(sessionLog.status, 0, sessionLog.stderr);
  assert.equal(sessionLog.result.base, SHA.base);
  assert.equal(sessionLog.result.head, SHA.head);
  assert.equal(sessionLog.result.ancestor, true);
  assert.deepEqual(sessionLog.result.capabilities, ['documentation']);
  assert.deepEqual(sessionLog.result.paths, [{
    path: 'docs/internal/web-gui-audit-20260930/SESSION_LOG.md', status: 'M', impact: 'documentation',
    reason: 'DOCUMENTATION_TREE', leaf: null
  }]);
  assert.deepEqual(sessionLog.result.verdicts, all);

  const t13 = await classifyOutcome({ diff: ['M', VIBRIO_FULL_GENERATION_SPEC] });
  assert.deepEqual(t13.result.paths[0].leaf, { kind: 'release-only', markers: null, jobs: [] });
  assert.deepEqual(t13.result.verdicts, { ...all, releaseEvidenceCarries: false });

  const t11 = await classifyOutcome({ diff: [
    'M', 'gbdraw/web/gallery/sessions/example.gbdraw-session.json',
    'M', 'docs/images/tutorials/example.png',
    'M', 'docs/internal/web-gui-audit-20260930/SESSION_LOG.md'
  ] });
  assert.equal(t11.result.verdicts.generatedArtifactChecksCarry, false);

  const python = await classifyOutcome({
    diff: ['M', 'tests/test_recipe_example.py'],
    blobs: { 'tests/test_recipe_example.py': 'import pytest\npytestmark = pytest.mark.recipe\n' }
  });
  assert.deepEqual(python.result.paths[0].leaf, { kind: 'python', markers: ['recipe'], jobs: ['core', 'recipes-standard'] });
  assert.deepEqual(python.result.verdicts, { ...all, releaseEvidenceCarries: false });

  const unrelated = await classifyOutcome({ diff: ['M', 'docs/FAQ.md'], ancestor: false });
  assert.equal(unrelated.result.ancestor, false);
  assert.deepEqual(unrelated.result.verdicts, {
    releaseEvidenceCarries: false, generatedArtifactChecksCarry: false, localTestEvidenceCarries: false
  });

  const identical = await classifyOutcome({ diff: null, trees: [TREE.same, TREE.same] });
  assert.equal(identical.result.identicalTree, true);
  assert.deepEqual(identical.result.verdicts, all);

  for (const argv of [['classify'], ['classify', '--base', SHA.base], ['classify', '--head', SHA.head, '--base'],
    ['classify', '--base', SHA.base, '--head', SHA.head, '--extra', 'x']]) {
    const invalid = await classifyOutcome({}, argv);
    assert.equal(invalid.status, 1, argv.join(' '));
    assert.match(invalid.stderr, /INVALID_ARGUMENTS/);
  }
});

test('functional shards skip setup when no changed spec is assigned to them', () => {
  const workflow = readFileSync(resolve(REPOSITORY_ROOT, '.github/workflows/test.yml'), 'utf8');
  const job = workflow.match(/\n  playwright-functional:\n[\s\S]*?(?=\n  [a-z0-9-]+:\n|$)/)?.[0] || '';
  const step = (name) => job.match(new RegExp(`\\n      - name: ${name}\\n[\\s\\S]*?(?=\\n      - |$)`))?.[0] || '';
  const select = step('Select the changed specs of this shard');
  assert.match(select, /id: shard/);
  assert.match(select, /CI_IMPACT_PLAN_JSON: \$\{\{ needs\.ci-impact\.outputs\.plan \}\}/);
  assert.match(select, /run: node tools\/ci-impact\.mjs shard-files \$\{\{ matrix\.shard \}\}/);
  assert.ok(job.indexOf(step('Set up Node.js')) < job.indexOf(select));
  for (const name of [
    'Set up Python 3.11',
    'Install Playwright functional dependencies',
    'Restore Playwright system packages',
    'Install Playwright Chromium and system packages',
    'Prepare browser wheel',
    'Run Playwright functional tests'
  ]) {
    assert.ok(job.indexOf(select) < job.indexOf(step(name)), name);
    assert.match(step(name), /\n        if: steps\.shard\.outputs\.files != ''\n/, name);
  }
  const run = step('Run Playwright functional tests');
  assert.match(run, /FUNCTIONAL_SPECS: \$\{\{ steps\.shard\.outputs\.files \}\}/);
  assert.match(run, /files="\$FUNCTIONAL_SPECS"/);
  assert.doesNotMatch(run, /require\('\.\/tests\/ci\/functional-shards\.json'\)/);
});
