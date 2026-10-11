#!/usr/bin/env node
// Assigns each functional Playwright spec file to one CI shard by measured duration.
//
//   node tools/balance-functional-shards.mjs [functional-report.json ...]
//
// Pass the test-results/functional-report.json of every shard of one or more green
// `Tests` runs to re-measure the minutes of each spec file. Without arguments, the
// minutes recorded in tests/ci/functional-shards.json are reused. A spec file without
// a measurement is estimated from its case count. Files are then assigned longest
// first to the shard with the least work, and tests/ci/functional-shards.json is
// rewritten.
import { spawnSync } from 'node:child_process';
import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import path from 'node:path';

const SHARDS_FILE = 'tests/ci/functional-shards.json';
const WORKERS_PER_SHARD = 2;
const DEFAULT_SHARD_COUNT = 8;

const byText = (a, b) => (a < b ? -1 : a > b ? 1 : 0);
const sum = (values) => values.reduce((total, value) => total + value, 0);
const median = (values) => {
  const sorted = [...values].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
};

// Playwright reports spec files relative to the test directory. The shard file
// lists them relative to the repository root (the working directory), which is
// also how the CLI file filters in the workflow name them.
const specsOf = (report) => {
  const visit = (suites) => suites.flatMap((suite) => [
    ...(suite.specs || []),
    ...visit(suite.suites || [])
  ]);
  const root = process.cwd().split(path.sep).join('/');
  const testDir = report.config.rootDir.split(path.sep).join('/');
  return visit(report.suites).map((spec) => ({
    ...spec, file: path.posix.relative(root, path.posix.join(testDir, spec.file))
  }));
};

const listFunctionalCases = () => {
  const cli = createRequire(import.meta.url).resolve('@playwright/test/cli');
  const result = spawnSync(process.execPath, [
    cli, 'test', '--config=tests/web/playwright/functional.config.js', '--list', '--reporter=json'
  ], { encoding: 'utf8', maxBuffer: 16 * 1024 * 1024 });
  if (result.status !== 0) throw new Error(result.stderr || result.stdout);
  const casesPerFile = new Map();
  for (const spec of specsOf(JSON.parse(result.stdout))) {
    casesPerFile.set(spec.file, (casesPerFile.get(spec.file) || 0) + spec.tests.length);
  }
  return casesPerFile;
};

// Each green run reports a case once, retries included. A file's minutes are the
// sum of the median minutes of its cases over all given reports.
const measureMinutes = (reportFiles) => {
  const caseMinutes = new Map();
  for (const reportFile of reportFiles) {
    for (const spec of specsOf(JSON.parse(readFileSync(reportFile, 'utf8')))) {
      for (const [index, testCase] of spec.tests.entries()) {
        const key = JSON.stringify([spec.file, spec.line, spec.column, spec.title, index]);
        const minutes = sum(testCase.results.map((result) => result.duration)) / 60000;
        caseMinutes.set(key, [...(caseMinutes.get(key) || []), minutes]);
      }
    }
  }
  const fileMinutes = new Map();
  for (const [key, values] of caseMinutes) {
    const [file] = JSON.parse(key);
    fileMinutes.set(file, (fileMinutes.get(file) || 0) + median(values));
  }
  return fileMinutes;
};

const existing = existsSync(SHARDS_FILE) ? JSON.parse(readFileSync(SHARDS_FILE, 'utf8')) : null;
const casesPerFile = listFunctionalCases();
const reports = process.argv.slice(2);
const fileMinutes = reports.length
  ? measureMinutes(reports)
  : new Map(Object.entries(existing?.minutes || {}));
const measuredFiles = [...casesPerFile.keys()].filter((file) => fileMinutes.has(file));
const minutesPerCase = sum(measuredFiles.map((file) => fileMinutes.get(file)))
  / sum(measuredFiles.map((file) => casesPerFile.get(file))) || 1;

const minutes = Object.fromEntries([...casesPerFile]
  .map(([file, cases]) => [
    file, Math.round((fileMinutes.get(file) ?? cases * minutesPerCase) * 100) / 100
  ])
  .sort(([a], [b]) => byText(a, b)));
const shards = Array.from(
  { length: existing?.shards.length || DEFAULT_SHARD_COUNT },
  () => ({ minutes: 0, files: [] })
);
for (const [file, value] of Object.entries(minutes)
  .sort(([a, x], [b, y]) => y - x || byText(a, b))) {
  const lightest = shards.reduce((best, shard) => (shard.minutes < best.minutes ? shard : best));
  lightest.minutes += value;
  lightest.files.push(file);
}
writeFileSync(SHARDS_FILE, `${JSON.stringify({
  minutes,
  shards: shards.map((shard) => shard.files.sort(byText))
}, null, 2)}\n`);
for (const [index, shard] of shards.entries()) {
  const estimate = (shard.minutes / WORKERS_PER_SHARD).toFixed(1);
  console.log(`shard ${index + 1}: ${shard.files.length} files, about ${estimate} min`);
}
