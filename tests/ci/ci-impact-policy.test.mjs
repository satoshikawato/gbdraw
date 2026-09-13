import assert from 'node:assert/strict';
import test from 'node:test';
import {
  classifyChanges,
  classifyPath,
  createImpactPlan,
  knownJobsFor,
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
    '.dockerignore',
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
    ['docs/internal/PRODUCT_IMPACT_RATCHET.md', 'ci-only'],
    ['pyproject.toml', 'packaging'],
    ['package-lock.json', 'packaging'],
    ['tools/helper.mjs', 'full'],
    ['gbdraw/future/new_engine.py', 'full'],
    ['gbdraw/web/new-runtime.bin', 'full'],
    ['future/unknown.file', 'full']
  ]) assert.equal(classifyPath(path).impact, impact, path);
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
    profile: 'pr', impact: 'full', decision: 'full'
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
    'web-pr-smoke'
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
    ['gbdraw/render/drawers/linear/features.py', ['web-change-budget', 'core-pr', 'lint', 'web-contracts-pr', 'web-pr-smoke']],
    ['gbdraw/web/js/app/label-editor.js', ['web-change-budget', 'web-contracts-pr', 'web-pr-smoke']],
    ['gbdraw/web/js/services/session-request.js', ['web-change-budget', 'core-pr', 'recipes-standard', 'lint', 'web-contracts-pr', 'web-pr-smoke']],
    ['gbdraw/web/gallery/examples.json', ['web-change-budget', 'gallery', 'web-contracts-pr', 'web-pr-smoke']]
  ]) assert.deepEqual(jobsForPaths([{ status: 'M', paths: [path] }]), expected, path);
});

test('mixed capabilities and both rename endpoints contribute independent required jobs', () => {
  const paths = ['gbdraw/render/drawers/linear/features.py', 'gbdraw/web/gallery/examples.json'];
  const expected = ['web-change-budget', 'core-pr', 'gallery', 'lint', 'web-contracts-pr', 'web-pr-smoke'];
  for (const changes of [
    paths.map((path) => ({ status: 'M', paths: [path] })),
    [{ status: 'R100', paths }],
    [{ status: 'C75', paths }]
  ]) assert.deepEqual(jobsForPaths(changes), expected);
  assert.ok(jobsForPaths([{ status: 'D', paths: [paths[0]] }]).includes('core-pr'));
  assert.throws(() => jobsForPaths([{ status: 'R100', paths: [paths[0], 'future/engine.py'] }]), /full coverage/);
  assert.throws(() => jobsForPaths([{ status: 'D', paths: ['gbdraw/unknown.py'] }]), /full coverage/);
});

test('control-plane, dependency, unknown and test-only changes cannot select partial coverage', () => {
  for (const impact of ['ci-only', 'packaging', 'full', 'tests-only']) {
    assert.throws(() => requiredJobsFor({ profile: 'pr', impact, decision: 'selective' }), /full coverage/);
    assert.deepEqual(requiredJobsFor({ profile: 'pr', impact, decision: 'full' }), knownJobsFor('pr'));
  }
});

test('release retains every dev functional job plus supported-version and slow acceptance', () => {
  assert.deepEqual(knownJobsFor('release'), [...knownJobsFor('dev'), 'acceptance-supported-main', 'slow-main']);
  assert.throws(() => requiredJobsFor({ profile: 'release', impact: 'documentation', decision: 'selective' }), /full coverage/);
  assert.throws(() => requiredJobsFor({ profile: 'dev', impact: 'web-runtime', decision: 'selective' }), /full coverage/);
});

test('capability tampering cannot drop an independent contribution', () => {
  const combined = selectivePlan({ impact: 'web-runtime', capabilities: ['documentation', 'web-runtime'] });
  assert.deepEqual(combined.requiredJobs, ['web-change-budget', 'recipes-standard', 'web-contracts-pr', 'web-pr-smoke']);
  for (const capabilities of [[], ['web-runtime', 'documentation'], ['documentation', 'documentation'], ['future']]) {
    assert.throws(() => validateImpactPlan({ ...combined, capabilities }), /Capabilities/);
  }
  assert.throws(() => validateImpactPlan({ ...combined, requiredJobs: ['web-change-budget', 'web-contracts-pr', 'web-pr-smoke'] }), /required jobs/);
});


test('ordinary Web changes stay selective when accompanied by their regression tests', () => {
  assert.deepEqual(jobsForPaths([
    { status: 'M', paths: ['gbdraw/web/js/app/label-editor.js'] },
    { status: 'A', paths: ['tests/web/label-editor.test.mjs'] },
    { status: 'M', paths: ['tests/web/right-drawer.playwright.spec.js'] }
  ]), ['web-change-budget', 'web-contracts-pr', 'web-pr-smoke']);
  for (const path of ['tests/web/session-request.test.mjs', 'tests/web/contracts/current-session-lazy-materialization.playwright.spec.js']) {
    assert.equal(classifyPath(path).impact, 'session-persistence', path);
    assert.ok(jobsForPaths([{ status: 'M', paths: [path] }]).includes('core-pr'));
  }
  for (const path of ['tests/conftest.py', 'tests/test_inputs/example.gbk', 'tests/web/architecture-ratchet-fixtures.test.mjs']) {
    assert.throws(() => jobsForPaths([{ status: 'M', paths: [path] }]), /full coverage/, path);
  }
});
