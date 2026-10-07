import assert from 'node:assert/strict';
import test from 'node:test';
import {
  GALLERY_PARITY_SPEC,
  IMPACT_PLAN_SCHEMA_VERSION,
  VIBRIO_FULL_GENERATION_SPEC,
  carryForwardVerdicts,
  classifyChanges,
  classifyPath,
  createImpactPlan,
  knownJobsFor,
  leafJobsFor,
  leafTestKind,
  requiredJobsFor,
  validateImpactPlan
} from '../../tools/ci-impact-policy.mjs';

const SHA = Object.freeze({
  base: 'a'.repeat(40),
  head: 'b'.repeat(40),
  workflow: 'c'.repeat(40)
});

const evidence = () => ({
  workflowPath: '.github/workflows/test.yml',
  aggregateName: 'Dev staging / gate',
  headSha: SHA.base,
  runId: 101,
  aggregateJobId: 202,
  runUrl: 'https://github.com/satoshikawato/gbdraw/actions/runs/101',
  aggregateJobUrl: 'https://github.com/satoshikawato/gbdraw/actions/runs/101/job/202'
});

const selectivePlan = (overrides = {}) => createImpactPlan({
  profile: 'pr',
  impact: 'documentation',
  decision: 'selective',
  basis: 'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE',
  changeBaseSha: SHA.base,
  changeHeadSha: SHA.head,
  workflowSha: SHA.workflow,
  changedPathCount: 1,
  inheritedEvidence: evidence(),
  ...overrides
});

test('metadata allowlist is deliberately narrow', () => {
  for (const path of [
    '.agents/skills/example/SKILL.md',
    '.claude/settings.json',
    '.codex/config.toml',
    '.cursor/rules/example.mdc',
    '.github/pull_request_template.md',
    '.gitattributes',
    '.gitignore',
    'CITATION.cff',
    'LICENSE.txt',
    'LICENSE_LIBERATION_FONTS.txt'
  ]) {
    assert.equal(classifyPath(path).impact, 'metadata', path);
  }
  assert.equal(classifyPath('.github/ISSUE_TEMPLATE/bug.md').impact, 'full');
  assert.equal(classifyPath('reports/ci.md').impact, 'full');
});

test('root Markdown and docs tree are documentation', () => {
  assert.equal(classifyPath('README.md').impact, 'documentation');
  assert.equal(classifyPath('NEW_PLAN.md').impact, 'documentation');
  assert.equal(classifyPath('docs/internal/policy.md').impact, 'documentation');
  assert.equal(classifyPath('docs/TUTORIALS/1_Intro.md').impact, 'documentation');
  assert.equal(classifyPath('docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md').impact, 'documentation');
  assert.equal(classifyPath('docs/internal/PRODUCT_IMPACT_RATCHET_FINAL_ACCEPTANCE_2026-08-27.md').impact, 'documentation');
  assert.equal(classifyPath('nested/README.md').impact, 'full');
  assert.equal(classifyPath('README.MD').impact, 'full');
});

test('subsystem paths classify without making unknown production paths selective', () => {
  for (const [path, impact] of [
    ['gbdraw/core/sequence.py', 'python-core'],
    ['gbdraw/render/drawers/linear/features.py', 'renderer'],
    ['gbdraw/web/js/app/label-editor.js', 'web-runtime'],
    ['gbdraw/session_request_codec.py', 'session-persistence'],
    ['gbdraw/web/js/services/session-request.js', 'session-persistence'],
    ['gbdraw/web/gallery/examples.json', 'gallery'],
    ['gbdraw/web/js/services/losat.js', 'losat-integration'],
    ['tests/test_regression.py', 'tests-only'],
    ['.github/workflows/test.yml', 'ci-only'],
    ['playwright.config.js', 'ci-only'],
    ['tests/web/playwright/functional.config.js', 'ci-only'],
    ['playwright.functional.config.js', 'full'],
    ['docs/internal/PRODUCT_IMPACT_RATCHET.md', 'policy-documentation'],
    ['pyproject.toml', 'packaging'],
    ['package-lock.json', 'packaging'],
    ['tools/helper.mjs', 'full'],
    ['gbdraw/future/new_engine.py', 'full'],
    ['gbdraw/web/new-runtime.bin', 'full'],
    ['future/unknown.file', 'full']
  ]) assert.equal(classifyPath(path).impact, impact, path);
});

test('only exact normative policy documents receive the policy-documentation route', () => {
  for (const path of [
    'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
    'docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md',
    'docs/internal/PRODUCT_IMPACT_RATCHET.md',
    'docs/internal/SELECTIVE_CI.md',
    'docs/internal/WEB_CHANGE_POLICY.md'
  ]) {
    assert.equal(classifyPath(path).impact, 'policy-documentation', path);
  }
  assert.equal(classifyPath('docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT_COPY.md').impact, 'documentation');
});

test('mixed changes preserve the capability union', () => {
  const metadata = classifyChanges([
    { status: 'M', paths: ['.gitignore'] },
    { status: 'A', paths: ['.agents/skills/new/SKILL.md'] }
  ]);
  assert.equal(metadata.impact, 'metadata');

  const documentation = classifyChanges([
    { status: 'M', paths: ['.gitignore'] },
    { status: 'M', paths: ['docs/FAQ.md'] }
  ]);
  assert.equal(documentation.impact, 'documentation');

  const full = classifyChanges([
    { status: 'M', paths: ['docs/FAQ.md'] },
    { status: 'M', paths: ['gbdraw/cli.py'] }
  ]);
  assert.equal(full.impact, 'python-core');
  assert.deepEqual(full.capabilities, ['documentation', 'python-core']);
});

test('rename, copy, and delete classify every relevant path', () => {
  const rename = classifyChanges([
    { status: 'R100', paths: ['docs/FAQ.md', 'gbdraw/FAQ.md'] }
  ]);
  assert.equal(rename.impact, 'full');
  assert.equal(rename.changedPathCount, 2);

  const copy = classifyChanges([
    { status: 'C75', paths: ['.agents/source.md', 'docs/copied.md'] }
  ]);
  assert.equal(copy.impact, 'documentation');
  assert.equal(copy.changedPathCount, 2);

  const deleted = classifyChanges([{ status: 'D', paths: ['docs/retired.md'] }]);
  assert.equal(deleted.impact, 'documentation');
  assert.equal(deleted.changedPathCount, 1);
});

test('empty, invalid, and unknown Git changes fail closed', () => {
  assert.deepEqual(
    { impact: classifyChanges([]).impact, valid: classifyChanges([]).valid },
    { impact: 'full', valid: false }
  );
  assert.equal(classifyChanges([{ status: 'T', paths: ['docs/FAQ.md'] }]).valid, false);
  assert.equal(classifyChanges([{ status: 'R101', paths: ['a', 'b'] }]).valid, false);
  assert.equal(classifyChanges([{ status: 'M', paths: ['../outside.md'] }]).valid, false);
});

test('profile job registries are exact and centralized', () => {
  assert.deepEqual(requiredJobsFor({
    profile: 'pr', impact: 'metadata', decision: 'selective'
  }), []);
  assert.deepEqual(requiredJobsFor({
    profile: 'pr', impact: 'documentation', decision: 'selective'
  }), ['recipes-standard']);
  assert.deepEqual(requiredJobsFor({
    profile: 'pr', impact: 'policy-documentation', decision: 'selective'
  }), ['web-change-budget']);
  assert.deepEqual(requiredJobsFor({
    profile: 'dev', impact: 'policy-documentation', decision: 'selective'
  }), ['web-change-budget']);
  assert.deepEqual(requiredJobsFor({
    profile: 'pr', impact: 'full', decision: 'full'
  }), [
    'web-change-budget',
    'core-pr',
    'recipes-standard',
    'gallery',
    'lint',
    'web-contracts-pr',
    'web-pr-smoke',
    'playwright-functional'
  ]);
  assert.deepEqual(requiredJobsFor({
    profile: 'pr', impact: 'ci-only', decision: 'full'
  }), [
    'web-change-budget',
    'core-pr',
    'recipes-standard',
    'gallery',
    'lint',
    'web-contracts-pr',
    'web-pr-smoke'
  ]);
  assert.deepEqual(knownJobsFor('pr'), [
    'web-change-budget',
    'core-pr',
    'recipes-standard',
    'gallery',
    'lint',
    'web-contracts-pr',
    'web-pr-smoke',
    'playwright-functional'
  ]);
  assert.deepEqual(knownJobsFor('dev'), [
    'web-change-budget',
    'core',
    'recipes-standard',
    'gallery',
    'browser',
    'playwright-functional',
    'playwright-performance',
    'lint',
    'losat-cache-browser-acceptance'
  ]);
  assert.deepEqual(knownJobsFor('gallery'), ['browser', 'performance']);
  assert.deepEqual(requiredJobsFor({
    profile: 'gallery', impact: 'metadata', decision: 'selective'
  }), []);
  assert.deepEqual(requiredJobsFor({
    profile: 'gallery', impact: 'documentation', decision: 'selective'
  }), []);
  assert.deepEqual(requiredJobsFor({
    profile: 'gallery', impact: 'policy-documentation', decision: 'selective'
  }), []);
  assert.deepEqual(requiredJobsFor({
    profile: 'gallery', impact: 'full', decision: 'full'
  }), ['browser', 'performance']);
});

test('impact plans are strict, policy-derived, and immutable', () => {
  const plan = selectivePlan();
  assert.equal(validateImpactPlan(plan, {
    profile: 'pr', workflowSha: SHA.workflow
  }), true);
  assert.deepEqual(plan.requiredJobs, ['recipes-standard']);
  assert.equal(Object.isFrozen(plan), true);
  assert.equal(Object.isFrozen(plan.requiredJobs), true);
  assert.equal(Object.isFrozen(plan.inheritedEvidence), true);

  const wrongJobs = { ...plan, requiredJobs: [] };
  assert.throws(() => validateImpactPlan(wrongJobs), /required jobs do not match policy/);
  assert.throws(
    () => validateImpactPlan({ ...plan, schemaVersion: 999 }),
    /schema version is not supported/
  );
  assert.throws(
    () => validateImpactPlan({ ...plan, profile: 'future' }),
    /profile is not supported/
  );
  assert.throws(
    () => validateImpactPlan({ ...plan, workflowSha: 'd'.repeat(40) }, {
      workflowSha: SHA.workflow
    }),
    /workflow SHA does not match/
  );
  assert.throws(
    () => validateImpactPlan({ ...plan, unexpected: true }),
    /invalid schema/
  );
});

const jobsForPaths = (changes) => {
  const classification = classifyChanges(changes);
  return requiredJobsFor({ profile: 'pr', ...classification, decision: 'selective' });
};

test('representative PR routes require the changed subsystem and cross-layer smoke', () => {
  for (const [path, expected] of [
    ['README.md', ['recipes-standard']],
    ['gbdraw/render/drawers/linear/features.py', ['web-change-budget', 'core-pr', 'gallery', 'lint', 'web-contracts-pr', 'web-pr-smoke']],
    ['gbdraw/web/js/app/label-editor.js', ['web-change-budget', 'gallery', 'web-contracts-pr', 'web-pr-smoke', 'playwright-functional']],
    ['gbdraw/web/js/services/session-request.js', ['web-change-budget', 'core-pr', 'recipes-standard', 'gallery', 'lint', 'web-contracts-pr', 'web-pr-smoke', 'playwright-functional']],
    ['gbdraw/web/gallery/examples.json', ['web-change-budget', 'gallery', 'web-contracts-pr', 'web-pr-smoke', 'playwright-functional']],
    ['gbdraw/web/js/services/losat.js', ['web-change-budget', 'core-pr', 'gallery', 'lint', 'web-contracts-pr', 'web-pr-smoke', 'playwright-functional']]
  ]) assert.deepEqual(jobsForPaths([{ status: 'M', paths: [path] }]), expected, path);
});

test('mixed capabilities and both rename endpoints contribute independent required jobs', () => {
  const paths = ['gbdraw/render/drawers/linear/features.py', 'gbdraw/web/gallery/examples.json'];
  const expected = ['web-change-budget', 'core-pr', 'gallery', 'lint', 'web-contracts-pr', 'web-pr-smoke', 'playwright-functional'];
  for (const changes of [
    paths.map((path) => ({ status: 'M', paths: [path] })),
    [{ status: 'R100', paths }],
    [{ status: 'C75', paths }]
  ]) assert.deepEqual(jobsForPaths(changes), expected);
  assert.ok(jobsForPaths([{ status: 'D', paths: [paths[0]] }]).includes('core-pr'));
  assert.throws(() => jobsForPaths([{ status: 'R100', paths: [paths[0], 'future/engine.py'] }]), /full coverage/);
  assert.throws(() => jobsForPaths([{ status: 'D', paths: ['gbdraw/unknown.py'] }]), /full coverage/);
});

test('documentation-only edits, copies, renames, and deletes use only documentation jobs', () => {
  const changes = [
    { status: 'M', paths: ['docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md'] },
    { status: 'M', paths: ['docs/TUTORIALS/1_Intro.md'] },
    { status: 'R100', paths: ['docs/old.md', 'docs/new.md'] },
    { status: 'D', paths: ['docs/retired.md'] }
  ];
  assert.deepEqual(jobsForPaths(changes), ['web-change-budget', 'recipes-standard']);
  assert.deepEqual(jobsForPaths([changes[0]]), ['web-change-budget']);
  assert.deepEqual(jobsForPaths([changes[1]]), ['recipes-standard']);
  assert.deepEqual(jobsForPaths([{
    status: 'R100',
    paths: ['docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md', 'docs/internal/retired-contract.md']
  }]), ['web-change-budget', 'recipes-standard']);
});

test('control-plane, dependency, unknown and test-only changes cannot select partial coverage', () => {
  for (const impact of ['ci-only', 'packaging', 'full', 'tests-only']) {
    assert.throws(() => requiredJobsFor({ profile: 'pr', impact, decision: 'selective' }), /full coverage/);
    const fullTier = knownJobsFor('pr').filter((job) => job !== 'playwright-functional');
    assert.deepEqual(
      requiredJobsFor({ profile: 'pr', impact, decision: 'full' }),
      ['full', 'tests-only'].includes(impact) ? knownJobsFor('pr') : fullTier,
      impact
    );
  }
});

test('functional Playwright joins PR plans only through runtime-facing capabilities', () => {
  const functional = ['web-runtime', 'session-persistence', 'gallery', 'losat-integration', 'tests-only', 'full'];
  for (const capability of [
    'metadata', 'documentation', 'policy-documentation', 'tests-only', 'python-core', 'renderer',
    'web-runtime', 'session-persistence', 'gallery', 'losat-integration', 'packaging', 'ci-only', 'full'
  ]) {
    const expected = functional.includes(capability);
    // Evidence fallbacks and the architecture-change label use the full decision.
    assert.equal(
      requiredJobsFor({ profile: 'pr', impact: capability, decision: 'full' }).includes('playwright-functional'),
      expected,
      `${capability} full`
    );
    if (!['tests-only', 'packaging', 'ci-only', 'full'].includes(capability)) {
      assert.equal(
        requiredJobsFor({ profile: 'pr', impact: capability, decision: 'selective' }).includes('playwright-functional'),
        expected,
        `${capability} selective`
      );
    }
  }
  for (const [capabilities, expected] of [
    [['documentation', 'ci-only'], false],
    [['python-core', 'renderer', 'packaging'], false],
    [['web-runtime', 'ci-only'], true],
    [['documentation', 'gallery'], true]
  ]) {
    const plan = { profile: 'pr', impact: capabilities.at(-1), capabilities, decision: 'full' };
    assert.equal(requiredJobsFor(plan).includes('playwright-functional'), expected, capabilities.join('+'));
  }
  for (const profile of ['dev', 'release']) {
    assert.ok(requiredJobsFor({ profile, impact: 'ci-only', decision: 'full' }).includes('playwright-functional'), profile);
  }
});

test('release retains every dev functional job plus supported-version, slow, and Vibrio generation acceptance', () => {
  assert.deepEqual(knownJobsFor('release'), [...knownJobsFor('dev'), 'acceptance-supported-main', 'slow-main', 'vibrio-generate-release']);
  assert.equal(knownJobsFor('dev').includes('vibrio-generate-release'), false, 'Vibrio generation is release-only');
  assert.throws(() => requiredJobsFor({ profile: 'release', impact: 'documentation', decision: 'selective' }), /full coverage/);
  assert.throws(() => requiredJobsFor({ profile: 'dev', impact: 'web-runtime', decision: 'selective' }), /full coverage/);
});

test('capability tampering cannot drop an independent contribution', () => {
  const combined = selectivePlan({ impact: 'web-runtime', capabilities: ['documentation', 'web-runtime'] });
  assert.deepEqual(combined.requiredJobs, ['web-change-budget', 'recipes-standard', 'gallery', 'web-contracts-pr', 'web-pr-smoke', 'playwright-functional']);
  for (const capabilities of [[], ['web-runtime', 'documentation'], ['documentation', 'documentation'], ['future']]) {
    assert.throws(() => validateImpactPlan({ ...combined, capabilities }), /Capabilities/);
  }
  assert.throws(() => validateImpactPlan({ ...combined, requiredJobs: ['web-change-budget', 'gallery', 'web-contracts-pr', 'web-pr-smoke', 'playwright-functional'] }), /required jobs/);
  assert.throws(() => validateImpactPlan({ ...combined, requiredJobs: ['web-change-budget', 'recipes-standard', 'gallery', 'web-contracts-pr', 'web-pr-smoke'] }), /required jobs/);
});


test('ordinary Web changes stay selective when accompanied by their regression tests', () => {
  assert.deepEqual(jobsForPaths([
    { status: 'M', paths: ['gbdraw/web/js/app/label-editor.js'] },
    { status: 'A', paths: ['tests/web/label-editor.test.mjs'] },
    { status: 'M', paths: ['tests/web/right-drawer.playwright.spec.js'] }
  ]), ['web-change-budget', 'gallery', 'web-contracts-pr', 'web-pr-smoke', 'playwright-functional']);
  for (const path of ['tests/web/session-request.test.mjs', 'tests/web/contracts/current-session-lazy-materialization.playwright.spec.js']) {
    assert.equal(classifyPath(path).impact, 'session-persistence', path);
    assert.ok(jobsForPaths([{ status: 'M', paths: [path] }]).includes('core-pr'));
  }
  for (const path of ['tests/conftest.py', 'tests/test_inputs/example.gbk', 'tests/web/architecture-ratchet-fixtures.test.mjs']) {
    assert.throws(() => jobsForPaths([{ status: 'M', paths: [path] }]), /full coverage/, path);
  }
});

test('documentation-only PR basis requires exact documentary scope and no inherited evidence', () => {
  const documentary = selectivePlan({ basis: 'DOCUMENTATION_ONLY_PR', inheritedEvidence: null });
  assert.equal(validateImpactPlan(documentary), true);
  assert.deepEqual(documentary.requiredJobs, ['recipes-standard']);
  for (const impact of ['metadata', 'web-runtime', 'ci-only', 'packaging', 'full']) {
    assert.throws(() => validateImpactPlan({ ...documentary, impact, capabilities: [impact] }),
      { code: 'BASIS_DECISION_MISMATCH' });
  }
  for (const profile of ['dev', 'gallery', 'release']) {
    assert.throws(() => validateImpactPlan({ ...documentary, profile }),
      { code: 'BASIS_DECISION_MISMATCH' });
  }
  assert.throws(() => validateImpactPlan({ ...documentary, decision: 'full' }),
    { code: 'BASIS_DECISION_MISMATCH' });
  assert.throws(() => validateImpactPlan({ ...documentary, inheritedEvidence: evidence() }),
    { code: 'UNEXPECTED_EVIDENCE' });
  assert.throws(() => validateImpactPlan({ ...documentary, requiredJobs: [] }),
    { code: 'REQUIRED_JOBS_MISMATCH' });
  assert.throws(() => selectivePlan({ impact: 'web-runtime', inheritedEvidence: null }),
    { code: 'INVALID_EVIDENCE_SCHEMA' });
});

test('every selective PR smoke route also requires the Gallery parity owner', () => {
  for (const impact of ['python-core', 'renderer', 'web-runtime', 'session-persistence', 'gallery', 'losat-integration']) {
    const requiredJobs = requiredJobsFor({ profile: 'pr', impact, decision: 'selective' });
    assert.ok(requiredJobs.includes('web-pr-smoke'), impact);
    assert.ok(requiredJobs.includes('gallery'), impact);
    const candidate = selectivePlan({ impact });
    assert.throws(() => validateImpactPlan({
      ...candidate, requiredJobs: requiredJobs.filter((job) => job !== 'gallery')
    }), /required jobs/, impact);
  }
});


const parentEvidence = (profile) => ({
  ...evidence(),
  ...(profile === 'gallery' ? {
    workflowPath: '.github/workflows/gallery-publication.yml',
    aggregateName: 'Gallery readiness / gate'
  } : {})
});

const leafPlan = (overrides = {}) => {
  const profile = overrides.profile ?? 'dev';
  return createImpactPlan({
    profile,
    impact: 'web-runtime',
    capabilities: ['web-runtime'],
    decision: 'selective',
    basis: 'LEAF_TEST_CHANGE_WITH_DIRECT_PARENT_EVIDENCE',
    changeBaseSha: SHA.base,
    changeHeadSha: SHA.head,
    workflowSha: SHA.workflow,
    changedPathCount: 1,
    inheritedEvidence: parentEvidence(profile),
    leafTests: [{ path: 'tests/web/label-editor.test.mjs', kind: 'node', jobs: profile === 'dev' ? ['browser'] : [] }],
    ...overrides
  });
};

const identicalPlan = (overrides = {}) => {
  const profile = overrides.profile ?? 'dev';
  return createImpactPlan({
    profile,
    impact: 'none',
    capabilities: ['none'],
    decision: 'selective',
    basis: 'IDENTICAL_TREE_WITH_DIRECT_PARENT_EVIDENCE',
    changeBaseSha: SHA.base,
    changeHeadSha: SHA.head,
    workflowSha: SHA.workflow,
    changedPathCount: 0,
    inheritedEvidence: parentEvidence(profile),
    ...overrides
  });
};

test('leaf test kinds follow the routing table and exclude shared or control-plane test code', () => {
  for (const [path, kind] of [
    ['tests/test_regression.py', 'python'],
    ['tests/web/label-editor.test.mjs', 'node'],
    ['tests/web/contracts/vibrio-sequence-source-coverage.test.mjs', 'node'],
    ['tests/web/right-drawer.playwright.spec.js', 'functional'],
    ['tests/web/contracts/session-regenerate-intent.playwright.spec.js', 'functional'],
    ['tests/web/webapp-performance.playwright.spec.js', 'performance'],
    ['tests/web/vibrio-session-save.performance.playwright.spec.js', 'performance'],
    [GALLERY_PARITY_SPEC, 'gallery-parity'],
    [VIBRIO_FULL_GENERATION_SPEC, 'release-only']
  ]) assert.equal(leafTestKind(path), kind, path);
  for (const path of [
    'tests/conftest.py',
    'tests/utils/svg_compare.py',
    'tests/fixtures/sessions/example.json',
    'tests/test_inputs/example.gbk',
    'tests/reference_outputs/example.svg',
    'tests/run_losat_cache_browser_acceptance.py',
    'tests/web/helpers/app-lifecycle.cjs',
    'tests/web/fixtures/example.json',
    'tests/web/fake-svg-dom.mjs',
    'tests/web/decoration-continuity.playwright.py',
    'tests/web/contracts/unknown.serial.spec.js',
    'tests/web/architecture-contracts.test.mjs',
    'tests/web/product-impact-ratchet-fixtures.test.mjs',
    'tests/web/promotion-readiness.test.mjs',
    'tests/ci/ci-impact-cli.test.mjs',
    'tests/nested/test_example.py',
    'tests/web/playwright/functional.config.js',
    'docs/FAQ.md',
    'gbdraw/web/js/app.js'
  ]) assert.equal(leafTestKind(path), null, path);
});

test('the serial spec constants are the specs their Playwright configurations run', async () => {
  const { readFileSync } = await import('node:fs');
  const config = (name) => readFileSync(new URL(`../../${name}`, import.meta.url), 'utf8');
  assert.match(config('tests/web/playwright/gallery-publication.config.js'),
    new RegExp(`testMatch: '${GALLERY_PARITY_SPEC.split('/').pop().replaceAll('.', '\\.')}'`));
  assert.match(config('tests/web/playwright/vibrio.config.js'),
    new RegExp(`testMatch: '${VIBRIO_FULL_GENERATION_SPEC.split('/').pop().replaceAll('.', '\\.')}'`));
});

test('leaf jobs select only the jobs that run the changed test', () => {
  for (const [kind, dev, gallery] of [
    ['node', ['browser'], []],
    ['functional', ['playwright-functional'], []],
    ['performance', ['playwright-performance'], []],
    ['gallery-parity', [], ['browser']],
    ['release-only', [], []]
  ]) {
    assert.deepEqual(leafJobsFor({ profile: 'dev', kind }), dev, kind);
    assert.deepEqual(leafJobsFor({ profile: 'gallery', kind }), gallery, kind);
  }
  assert.deepEqual(leafJobsFor({ profile: 'dev', kind: 'python', markers: [] }), ['core']);
  assert.deepEqual(leafJobsFor({ profile: 'dev', kind: 'python', markers: ['browser', 'recipe'] }),
    ['core', 'recipes-standard', 'browser']);
  assert.deepEqual(leafJobsFor({ profile: 'dev', kind: 'python', markers: ['slow', 'gallery'] }),
    ['core', 'gallery']);
  assert.deepEqual(leafJobsFor({ profile: 'dev', kind: 'python', markers: null }),
    ['core', 'recipes-standard', 'gallery', 'browser']);
  assert.deepEqual(leafJobsFor({ profile: 'gallery', kind: 'python', markers: null }), []);
});

test('leaf-test plans run documentation jobs plus the jobs of the changed leaf tests', () => {
  const plan = leafPlan({
    impact: 'web-runtime',
    capabilities: ['documentation', 'tests-only', 'web-runtime'],
    changedPathCount: 4,
    leafTests: [
      { path: 'tests/test_regression.py', kind: 'python', jobs: ['core', 'browser'] },
      { path: 'tests/web/label-editor.test.mjs', kind: 'node', jobs: ['browser'] },
      { path: 'tests/web/right-drawer.playwright.spec.js', kind: 'functional', jobs: ['playwright-functional'] }
    ]
  });
  assert.equal(plan.schemaVersion, IMPACT_PLAN_SCHEMA_VERSION);
  assert.deepEqual(plan.requiredJobs, ['core', 'recipes-standard', 'browser', 'playwright-functional']);
  assert.equal(Object.isFrozen(plan.leafTests), true);
  assert.equal(validateImpactPlan(plan), true);

  const releaseOnly = leafPlan({
    impact: 'tests-only',
    capabilities: ['tests-only'],
    leafTests: [{ path: VIBRIO_FULL_GENERATION_SPEC, kind: 'release-only', jobs: [] }]
  });
  assert.deepEqual(releaseOnly.requiredJobs, []);

  const gallery = leafPlan({
    profile: 'gallery',
    impact: 'web-runtime',
    capabilities: ['documentation', 'tests-only', 'web-runtime'],
    leafTests: [
      { path: GALLERY_PARITY_SPEC, kind: 'gallery-parity', jobs: ['browser'] },
      { path: 'tests/web/label-editor.test.mjs', kind: 'node', jobs: [] }
    ]
  });
  assert.deepEqual(gallery.requiredJobs, ['browser']);
  assert.deepEqual(leafPlan({ profile: 'gallery' }).requiredJobs, []);
});

test('the validator accepts test subject capabilities only under a consistent leaf basis', () => {
  const plan = leafPlan();
  assert.throws(() => validateImpactPlan({ ...plan, basis: 'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE', leafTests: null }),
    { code: 'INVALID_SELECTIVE_PLAN' });
  assert.throws(() => validateImpactPlan({ ...plan, basis: 'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE' }),
    { code: 'UNEXPECTED_LEAF_TESTS' });
  for (const leafTests of [
    null,
    [],
    [{ path: 'tests/web/label-editor.test.mjs', kind: 'functional', jobs: ['playwright-functional'] }],
    [{ path: 'tests/web/label-editor.test.mjs', kind: 'node', jobs: [] }],
    [{ path: 'tests/web/label-editor.test.mjs', kind: 'node', jobs: ['browser', 'core'] }],
    [{ path: 'tests/web/label-editor.test.mjs', kind: 'node' }],
    [{ path: 'tests/web/helpers/app-lifecycle.cjs', kind: 'node', jobs: ['browser'] }],
    [{ path: 'tests/test_regression.py', kind: 'python', jobs: ['browser'] }],
    [{ path: 'tests/test_regression.py', kind: 'python', jobs: ['core', 'lint'] }],
    [
      { path: 'tests/web/right-drawer.playwright.spec.js', kind: 'functional', jobs: ['playwright-functional'] },
      { path: 'tests/web/label-editor.test.mjs', kind: 'node', jobs: ['browser'] }
    ]
  ]) {
    assert.throws(() => validateImpactPlan({ ...plan, leafTests }), { code: 'INVALID_LEAF_TESTS' },
      JSON.stringify(leafTests));
  }
  // Every non-documentation capability must come from a listed leaf test.
  assert.throws(() => leafPlan({ impact: 'session-persistence', capabilities: ['web-runtime', 'session-persistence'] }),
    { code: 'INVALID_LEAF_TESTS' });
  assert.throws(() => leafPlan({ impact: 'ci-only', capabilities: ['web-runtime', 'ci-only'] }),
    { code: 'INVALID_LEAF_TESTS' });
  assert.throws(() => leafPlan({ impact: 'documentation', capabilities: ['documentation'] }),
    { code: 'INVALID_LEAF_TESTS' });
  assert.throws(() => leafPlan({ profile: 'pr', inheritedEvidence: evidence() }),
    { code: 'BASIS_DECISION_MISMATCH' });
  assert.throws(() => leafPlan({ decision: 'full', inheritedEvidence: null }),
    { code: 'BASIS_DECISION_MISMATCH' });
  assert.throws(() => leafPlan({ inheritedEvidence: null }), { code: 'INVALID_EVIDENCE_SCHEMA' });
  assert.throws(() => validateImpactPlan({ ...plan, requiredJobs: [] }), { code: 'REQUIRED_JOBS_MISMATCH' });
  const { leafTests: _omitted, ...withoutField } = plan;
  assert.throws(() => validateImpactPlan(withoutField), { code: 'INVALID_PLAN_SCHEMA' });
});

test('identical-tree plans inherit everything only on dev and Gallery pushes', () => {
  for (const profile of ['dev', 'gallery']) {
    const plan = identicalPlan({ profile });
    assert.deepEqual(plan.requiredJobs, [], profile);
    assert.equal(plan.leafTests, null, profile);
    const fallback = createImpactPlan({
      profile,
      impact: 'none',
      capabilities: ['none'],
      decision: 'full',
      basis: 'INHERITED_EVIDENCE_UNAVAILABLE',
      changeBaseSha: SHA.base,
      changeHeadSha: SHA.head,
      workflowSha: SHA.workflow,
      changedPathCount: 0,
      inheritedEvidence: null
    });
    assert.deepEqual(fallback.requiredJobs, knownJobsFor(profile), profile);
  }
  assert.throws(() => identicalPlan({ profile: 'pr' }), { code: 'BASIS_DECISION_MISMATCH' });
  assert.throws(() => identicalPlan({ changedPathCount: 1 }), { code: 'BASIS_IMPACT_MISMATCH' });
  assert.throws(() => identicalPlan({ impact: 'metadata', capabilities: ['metadata'] }), { code: 'BASIS_IMPACT_MISMATCH' });
  assert.throws(() => identicalPlan({ inheritedEvidence: null }), { code: 'INVALID_EVIDENCE_SCHEMA' });
  assert.throws(() => identicalPlan({ decision: 'full', inheritedEvidence: null }), { code: 'BASIS_DECISION_MISMATCH' });
  assert.throws(() => validateImpactPlan({ ...identicalPlan(), leafTests: [] }), { code: 'UNEXPECTED_LEAF_TESTS' });
  // 'none' never mixes with a changed path and never reaches pull requests.
  assert.throws(() => selectivePlan({ impact: 'none', capabilities: ['none'], changedPathCount: 0 }),
    { code: 'BASIS_IMPACT_MISMATCH' });
  assert.throws(() => selectivePlan({
    impact: 'none', capabilities: ['none'], changedPathCount: 0, decision: 'full',
    basis: 'INHERITED_EVIDENCE_UNAVAILABLE', inheritedEvidence: null
  }), { code: 'BASIS_IMPACT_MISMATCH' });
  assert.throws(() => leafPlan({ impact: 'web-runtime', capabilities: ['none', 'web-runtime'] }),
    { code: 'BASIS_IMPACT_MISMATCH' });
  assert.equal(classifyChanges([]).valid, false, 'an empty diff alone stays invalid');
});

const verdictsFor = (paths, overrides = {}) => carryForwardVerdicts({
  ancestor: true,
  valid: true,
  identicalTree: false,
  entries: paths.map((entry) => (typeof entry === 'string'
    ? { path: entry, leaf: null }
    : entry)),
  ...overrides
});

test('carry-forward verdicts follow the allowed path sets', () => {
  const all = { releaseEvidenceCarries: true, generatedArtifactChecksCarry: true, localTestEvidenceCarries: true };
  const none = { releaseEvidenceCarries: false, generatedArtifactChecksCarry: false, localTestEvidenceCarries: false };
  // A SESSION_LOG-only move carries everything.
  assert.deepEqual(verdictsFor(['docs/internal/web-gui-audit-20260930/SESSION_LOG.md']), all);
  assert.deepEqual(verdictsFor(['.gitignore', 'docs/internal/SELECTIVE_CI.md', 'CHANGELOG.md']), all);
  // T13: the release-only Vibrio spec does not carry release evidence.
  assert.deepEqual(verdictsFor([{ path: VIBRIO_FULL_GENERATION_SPEC, leaf: { kind: 'release-only', markers: null } }]), {
    releaseEvidenceCarries: false, generatedArtifactChecksCarry: true, localTestEvidenceCarries: true
  });
  // T11: a Gallery Session and docs screenshots invalidate generated-artifact checks.
  assert.deepEqual(verdictsFor([
    'gbdraw/web/gallery/sessions/example.gbdraw-session.json',
    'docs/images/tutorials/example.png',
    'docs/internal/web-gui-audit-20260930/SESSION_LOG.md'
  ]), none);
  assert.equal(verdictsFor(['docs/images/tutorials/example.png']).generatedArtifactChecksCarry, false);
  assert.deepEqual(verdictsFor(['docs/images/tutorials/example.png']), {
    releaseEvidenceCarries: true, generatedArtifactChecksCarry: false, localTestEvidenceCarries: true
  });
  assert.deepEqual(verdictsFor(['docs/recipes/run_cli_scenarios.py']), none);
  assert.deepEqual(verdictsFor(['docs/capture/run_all.py']), {
    releaseEvidenceCarries: false, generatedArtifactChecksCarry: false, localTestEvidenceCarries: true
  });
  for (const [leaf, expected] of [
    [{ kind: 'node', markers: null }, all],
    [{ kind: 'functional', markers: null }, all],
    [{ kind: 'gallery-parity', markers: null }, all],
    [{ kind: 'python', markers: [] }, all],
    [{ kind: 'python', markers: ['recipe'] }, { ...all, releaseEvidenceCarries: false }],
    [{ kind: 'python', markers: ['slow'] }, { ...all, releaseEvidenceCarries: false, localTestEvidenceCarries: false }],
    [{ kind: 'python', markers: null }, { ...all, releaseEvidenceCarries: false, localTestEvidenceCarries: false }]
  ]) assert.deepEqual(verdictsFor([{ path: 'tests/test_regression.py', leaf }]), expected, JSON.stringify(leaf));
  assert.deepEqual(verdictsFor(['.github/workflows/test.yml']), {
    releaseEvidenceCarries: false, generatedArtifactChecksCarry: true, localTestEvidenceCarries: false
  });
  for (const path of ['tests/web/helpers/app-lifecycle.cjs', 'tests/fixtures/sessions/a.json', 'tests/reference_outputs/a.svg',
    'gbdraw/render/drawers/linear/features.py', 'pyproject.toml', 'tools/audit/parity_replay.py', 'future/unknown.file']) {
    assert.deepEqual(verdictsFor([path]), none, path);
  }
  assert.deepEqual(verdictsFor([], { identicalTree: true }), all);
  assert.deepEqual(verdictsFor(['docs/FAQ.md'], { ancestor: false }), none);
  assert.deepEqual(verdictsFor([], { identicalTree: true, ancestor: false }), none);
  assert.deepEqual(verdictsFor([], { valid: false }), none);
});
