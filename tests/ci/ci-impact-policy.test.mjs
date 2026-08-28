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

test('runtime, tests, workflows, dependencies, and unknown paths are full', () => {
  for (const path of [
    'gbdraw/core/sequence.py',
    'tests/test_regression.py',
    'tools/helper.mjs',
    '.github/workflows/test.yml',
    'pyproject.toml',
    'setup.py',
    'MANIFEST.in',
    'package.json',
    'package-lock.json',
    'playwright.config.js',
    'recipe/meta.yaml',
    'examples/example.gb',
    '_reproduced/result.svg',
    'wrangler.toml',
    'reports/result.json',
    'future/unknown.file'
  ]) {
    assert.equal(classifyPath(path).impact, 'full', path);
  }
});

test('mixed changes use the strongest impact', () => {
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
  assert.equal(full.impact, 'full');
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
  assert.deepEqual(knownJobsFor('pr'), [
    'web-change-budget',
    'core-pr',
    'recipes-standard',
    'gallery',
    'lint',
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
    'acceptance-supported-main',
    'slow-main',
    'lint',
    'losat-cache-browser-acceptance'
  ]);
  assert.deepEqual(knownJobsFor('gallery'), ['browser', 'performance']);
  assert.deepEqual(requiredJobsFor({
    profile: 'gallery', impact: 'documentation', decision: 'selective'
  }), []);
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
    () => validateImpactPlan({ ...plan, schemaVersion: 2 }),
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
