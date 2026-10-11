import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { createRequire } from 'node:module';
import { readFileSync } from 'node:fs';
import path from 'node:path';
import test from 'node:test';

const require = createRequire(import.meta.url);
const cli = require.resolve('@playwright/test/cli');
const list = (...args) => {
  const result = spawnSync(process.execPath, [cli, ...args, '--list', '--reporter=json'], {
    encoding: 'utf8', maxBuffer: 16 * 1024 * 1024
  });
  assert.equal(result.status, 0, result.stderr || result.stdout);
  const report = JSON.parse(result.stdout);
  assert.deepEqual(report.errors, []);
  return report;
};
const specsOf = (report) => {
  const visit = (suites) => suites.flatMap((suite) => [
    ...(suite.specs || []),
    ...visit(suite.suites || [])
  ]);
  return visit(report.suites);
};
const casesOf = (report) => specsOf(report)
  .flatMap((spec) => spec.tests.map(() => `${spec.file}:${spec.title}`));
const collect = (...args) => casesOf(list(...args));

test('expanded PR smoke has 8–19 cases and each remains in full functional acceptance', () => {
  const smoke = collect('test', '--config=tests/web/playwright/pr-smoke.config.js');
  const full = new Set(collect('test', '--config=tests/web/playwright/functional.config.js'));
  assert.ok(smoke.length >= 8 && smoke.length <= 19, `collected ${smoke.length} PR cases`);
  for (const title of smoke) assert.ok(full.has(title), `missing full regression: ${title}`);
  for (const path of [
    'mode-transition-editor-state.playwright.spec.js',
    'mode-record-identity.playwright.spec.js',
    'annotation-download.playwright.spec.js',
    'contracts/active-result-edit-transaction.playwright.spec.js'
  ]) assert.ok([...full].some((title) => title.startsWith(`${path}:`)), `missing moved regression: ${path}`);
});

test('comparison browser contracts run in the required PR contract job and full acceptance', () => {
  const scripts = JSON.parse(readFileSync('package.json', 'utf8')).scripts;
  const [binary, ...args] = scripts['test:web:comparison-contracts'].split(/\s+/);
  assert.equal(binary, 'playwright');
  const contracts = collect(...args);
  const full = new Set(collect('test', '--config=tests/web/playwright/functional.config.js'));
  const smoke = new Set(collect('test', '--config=tests/web/playwright/pr-smoke.config.js'));
  assert.equal(contracts.length, 16);
  assert.ok(contracts.includes('comparison-ui.playwright.spec.js:comparison controls drive appearance and current Session round trips'));
  assert.ok(contracts.includes('linear-multi-record.playwright.spec.js:Collinear inference checkbox skips self searches and reuses matching evidence'));
  assert.ok(contracts.includes('linear-multi-record.playwright.spec.js:protein raw cache survives cancellation and derived options preserve search identity'));
  for (const title of contracts) {
    assert.ok(full.has(title), `missing full regression: ${title}`);
    assert.ok(!smoke.has(title), `duplicated PR execution: ${title}`);
  }
  const workflow = readFileSync('.github/workflows/test.yml', 'utf8');
  const job = workflow.split('  web-contracts-pr:\n')[1].split('  web-pr-smoke:\n')[0];
  assert.match(job, /python -m pytest tests\/[\s\S]*-m "browser and not slow"/);
  const pytestEntry = readFileSync('tests/test_linear_comparison_browser_contracts.py', 'utf8');
  assert.match(pytestEntry, /@pytest\.mark\.browser/);
  assert.match(pytestEntry, /test:web:comparison-contracts/);
  assert.doesNotMatch(pytestEntry, /pytest\.skip|pytest\.mark\.slow/);
  assert.doesNotMatch(job, /continue-on-error/);
  assert.match(workflow.split('  pr-gate:\n')[1], /- web-contracts-pr/);
});

const SHARDS_FILE = 'tests/ci/functional-shards.json';
// Lists why `shards` does not assign each functional spec file to exactly one shard.
const assignmentProblems = (specFiles, shards) => {
  const owners = new Map();
  shards.forEach((files, index) => {
    for (const file of files) owners.set(file, [...(owners.get(file) || []), index + 1]);
  });
  return [
    ...specFiles.filter((file) => !owners.has(file))
      .map((file) => `${file} is not assigned to a shard`),
    ...[...owners].filter(([, shardNumbers]) => shardNumbers.length > 1)
      .map(([file, shardNumbers]) => `${file} is assigned to shards ${shardNumbers.join(', ')}`),
    ...[...owners.keys()].filter((file) => !specFiles.includes(file))
      .map((file) => `${file} is not a functional spec file`)
  ];
};

test('the shard assignment check rejects unassigned, duplicated, and unknown spec files', () => {
  const [a, b] = ['tests/web/a.playwright.spec.js', 'tests/web/b.playwright.spec.js'];
  assert.deepEqual(assignmentProblems([a, b], [[a], [b]]), []);
  assert.deepEqual(assignmentProblems([a, b], [[a], []]), [`${b} is not assigned to a shard`]);
  assert.deepEqual(assignmentProblems([a, b], [[a, b], [b]]), [`${b} is assigned to shards 1, 2`]);
  assert.deepEqual(
    assignmentProblems([a], [[a], ['tests/web/removed.playwright.spec.js']]),
    ['tests/web/removed.playwright.spec.js is not a functional spec file']
  );
});

test('functional CI shards run the checked-in spec file lists and every acceptance case once', () => {
  const workflow = readFileSync('.github/workflows/test.yml', 'utf8');
  const job = workflow.split('  playwright-functional:\n')[1]
    .split('  playwright-performance:\n')[0];
  const { shards } = JSON.parse(readFileSync(SHARDS_FILE, 'utf8'));
  const matrix = JSON.parse(job.match(/shard: (\[[\d, ]+\])/)[1]);
  assert.deepEqual(matrix, shards.map((_, index) => index + 1));
  // tools/ci-impact.mjs shard-files reads this shard map and keeps only changed specs for a
  // leaf-test plan; every other plan runs the shard's whole list.
  assert.match(job, /run: node tools\/ci-impact\.mjs shard-files \$\{\{ matrix\.shard \}\}\n/);
  assert.match(job, /files="\$FUNCTIONAL_SPECS"\n/);
  assert.match(job, /npm run test:web:functional-full -- --reporter=line,github,json \$files\n/);
  assert.doesNotMatch(job, /PWTEST_|--shard/);

  const full = list('test', '--config=tests/web/playwright/functional.config.js');
  // Spec files as the workflow names them: relative to the repository root.
  const root = process.cwd();
  const specFiles = [...new Set(specsOf(full).map((spec) => path
    .relative(root, path.join(full.config.rootDir, spec.file)).split(path.sep).join('/')))];
  assert.deepEqual(
    assignmentProblems(specFiles, shards), [],
    'run tools/balance-functional-shards.mjs; see docs/internal/SELECTIVE_CI.md'
  );
  // Playwright file filters are regular expressions, so list each shard to prove
  // that no filter also selects another shard's file.
  const partition = shards.flatMap((files) => collect(
    'test', '--config=tests/web/playwright/functional.config.js', ...files
  ));
  assert.equal(new Set(partition).size, partition.length, 'a case appears in multiple shards');
  assert.deepEqual(partition.sort(), casesOf(full).sort(), 'the CI matrix loses acceptance cases');
});
