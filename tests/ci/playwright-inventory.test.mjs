import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { createRequire } from 'node:module';
import test from 'node:test';

const require = createRequire(import.meta.url);
const cli = require.resolve('@playwright/test/cli');
const collect = (config) => {
  const result = spawnSync(process.execPath, [cli, 'test', `--config=${config}`, '--list', '--reporter=json'], {
    encoding: 'utf8', maxBuffer: 16 * 1024 * 1024
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

test('expanded PR smoke has 8–12 cases and each remains in full functional acceptance', () => {
  const smoke = collect('playwright.pr-smoke.config.js');
  const full = new Set(collect('playwright.functional.config.js'));
  assert.ok(smoke.length >= 8 && smoke.length <= 12, `collected ${smoke.length} PR cases`);
  for (const title of smoke) assert.ok(full.has(title), `missing full regression: ${title}`);
  for (const path of [
    'mode-transition-editor-state.playwright.spec.js',
    'mode-record-identity.playwright.spec.js',
    'annotation-download.playwright.spec.js',
    'contracts/active-result-edit-transaction.playwright.spec.js'
  ]) assert.ok([...full].some((title) => title.startsWith(`${path}:`)), `missing moved regression: ${path}`);
});
