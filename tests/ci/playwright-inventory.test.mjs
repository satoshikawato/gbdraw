import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { createRequire } from 'node:module';
import { readFileSync } from 'node:fs';
import test from 'node:test';

const require = createRequire(import.meta.url);
const cli = require.resolve('@playwright/test/cli');
const collectWith = (env, ...args) => {
  const result = spawnSync(process.execPath, [cli, ...args, '--list', '--reporter=json'], {
    encoding: 'utf8', maxBuffer: 16 * 1024 * 1024, env: { ...process.env, ...env }
  });
  assert.equal(result.status, 0, result.stderr || result.stdout);
  const report = JSON.parse(result.stdout);
  assert.deepEqual(report.errors, []);
  const visit = (suites) => suites.flatMap((suite) => [
    ...(suite.specs || []).flatMap((spec) => spec.tests.map(() => `${spec.file}:${spec.title}`)),
    ...visit(suite.suites || [])
  ]);
  return visit(report.suites);
};
const collect = (...args) => collectWith({}, ...args);

test('expanded PR smoke has 8–19 cases and each remains in full functional acceptance', () => {
  const smoke = collect('test', '--config=playwright.pr-smoke.config.js');
  const full = new Set(collect('test', '--config=playwright.functional.config.js'));
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
  const full = new Set(collect('test', '--config=playwright.functional.config.js'));
  const smoke = new Set(collect('test', '--config=playwright.pr-smoke.config.js'));
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

test('weighted functional CI shards partition the complete acceptance inventory without duplication', () => {
  const workflow = readFileSync('.github/workflows/test.yml', 'utf8');
  const job = workflow.split('  playwright-functional:\n')[1]
    .split('  playwright-performance:\n')[0];
  const shards = JSON.parse(job.match(/shard: (\[[\d, ]+\])/)[1]);
  const total = Number(job.match(/--shard=\$\{\{ matrix\.shard \}\}\/(\d+)/)[1]);
  assert.deepEqual(shards, Array.from({ length: total }, (_, index) => index + 1));
  const weightText = job.match(/\n {10}PWTEST_SHARD_WEIGHTS: "(\d+(?::\d+)*)"\n/)?.[1];
  assert.ok(weightText, 'the functional step must set PWTEST_SHARD_WEIGHTS');
  const weights = weightText.split(':').map(Number);
  assert.equal(weights.length, total);
  const full = collect('test', '--config=playwright.functional.config.js');
  const perShard = shards.map((shard) => collectWith(
    { PWTEST_SHARD_WEIGHTS: weightText },
    'test', '--config=playwright.functional.config.js', `--shard=${shard}/${total}`
  ));
  const partition = perShard.flat();
  assert.equal(new Set(partition).size, partition.length, 'a case appears in multiple shards');
  assert.deepEqual(partition.sort(), [...full].sort(), 'the CI matrix loses acceptance cases');
  // Shard sizes must follow the weights rather than equal counts, so a Playwright
  // release that ignores PWTEST_SHARD_WEIGHTS fails here. Cases that share a
  // serial group or beforeAll hook move between shards together, so sizes are near
  // the weighted split, not equal to it.
  const split = (shares) => {
    const sum = shares.reduce((total, share) => total + share, 0);
    const sizes = shares.map((share) => Math.floor(share * full.length / sum));
    const remainder = full.length - sizes.reduce((total, size) => total + size, 0);
    for (let index = 0; index < remainder; index += 1) sizes[index % sizes.length] += 1;
    return sizes;
  };
  const sizes = perShard.map((cases) => cases.length);
  const distance = (target) => target.reduce((sum, size, index) => sum + Math.abs(size - sizes[index]), 0);
  assert.ok(
    distance(split(weights)) < distance(split(weights.map(() => 1))),
    `shard sizes ${sizes.join(':')} ignore PWTEST_SHARD_WEIGHTS ${weightText}`
  );
});
