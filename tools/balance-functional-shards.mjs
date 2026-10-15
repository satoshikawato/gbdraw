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
// rewritten. tools/ci-inventory.mjs imports the same functions for its weekly report.
import { spawnSync } from 'node:child_process';
import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import path from 'node:path';
import { pathToFileURL } from 'node:url';

export const SHARDS_FILE = 'tests/ci/functional-shards.json';
export const WORKERS_PER_SHARD = 2;
const DEFAULT_SHARD_COUNT = 8;

const byText = (a, b) => (a < b ? -1 : a > b ? 1 : 0);
const sum = (values) => values.reduce((total, value) => total + value, 0);
const median = (values) => {
  const sorted = [...values].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
};

// Playwright reports spec files relative to its test directory (`config.rootDir`). The
// shard file lists them relative to the repository root, which is also how the CLI file
// filters in the workflow name them. `testDir` is that test directory relative to the
// repository root; a CI report's own rootDir is a path on the runner.
export const specsOf = (report, testDir) => {
  const visit = (suites) => suites.flatMap((suite) => [
    ...(suite.specs || []),
    ...visit(suite.suites || [])
  ]);
  return visit(report.suites).map((spec) => ({ ...spec, file: path.posix.join(testDir, spec.file) }));
};

// Lists the functional cases of the working tree: the test directory relative to the
// repository root (the working directory) and the case count of each spec file.
export const listFunctionalCases = () => {
  const cli = createRequire(import.meta.url).resolve('@playwright/test/cli');
  const result = spawnSync(process.execPath, [
    cli, 'test', '--config=tests/web/playwright/functional.config.js', '--list', '--reporter=json'
  ], { encoding: 'utf8', maxBuffer: 16 * 1024 * 1024 });
  if (result.status !== 0) throw new Error(result.stderr || result.stdout);
  const report = JSON.parse(result.stdout);
  const testDir = path.relative(process.cwd(), report.config.rootDir).split(path.sep).join('/');
  return { testDir, casesPerFile: casesPerFileOf([report], testDir) };
};

// The case count of each spec file in the given reports.
export const casesPerFileOf = (reports, testDir) => {
  const casesPerFile = new Map();
  for (const spec of reports.flatMap((report) => specsOf(report, testDir))) {
    casesPerFile.set(spec.file, (casesPerFile.get(spec.file) || 0) + spec.tests.length);
  }
  return casesPerFile;
};

// Each green run reports a case once, retries included. A file's minutes are the
// sum of the median minutes of its cases over all given reports.
export const measureMinutes = (reports, testDir) => {
  const caseMinutes = new Map();
  for (const report of reports) {
    for (const spec of specsOf(report, testDir)) {
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

// Assigns files longest first to the shard with the least work. Returns the minutes of
// each file (estimated from its case count when unmeasured), the sorted file list and the
// summed test minutes of each shard, and the files without a measurement.
export const planShards = ({ casesPerFile, fileMinutes, shardCount }) => {
  const measuredFiles = [...casesPerFile.keys()].filter((file) => fileMinutes.has(file));
  const minutesPerCase = sum(measuredFiles.map((file) => fileMinutes.get(file)))
    / sum(measuredFiles.map((file) => casesPerFile.get(file))) || 1;
  const minutes = Object.fromEntries([...casesPerFile]
    .map(([file, cases]) => [
      file, Math.round((fileMinutes.get(file) ?? cases * minutesPerCase) * 100) / 100
    ])
    .sort(([a], [b]) => byText(a, b)));
  const shards = Array.from({ length: shardCount }, () => ({ minutes: 0, files: [] }));
  for (const [file, value] of Object.entries(minutes)
    .sort(([a, x], [b, y]) => y - x || byText(a, b))) {
    const lightest = shards.reduce((best, shard) => (shard.minutes < best.minutes ? shard : best));
    lightest.minutes += value;
    lightest.files.push(file);
  }
  return {
    minutes,
    shards: shards.map((shard) => shard.files.sort(byText)),
    shardMinutes: shards.map((shard) => Math.round(shard.minutes * 100) / 100),
    unmeasured: [...casesPerFile.keys()].filter((file) => !fileMinutes.has(file)).sort(byText)
  };
};

const main = (reportFiles) => {
  const existing = existsSync(SHARDS_FILE) ? JSON.parse(readFileSync(SHARDS_FILE, 'utf8')) : null;
  const { testDir, casesPerFile } = listFunctionalCases();
  const plan = planShards({
    casesPerFile,
    fileMinutes: reportFiles.length
      ? measureMinutes(reportFiles.map((file) => JSON.parse(readFileSync(file, 'utf8'))), testDir)
      : new Map(Object.entries(existing?.minutes || {})),
    shardCount: existing?.shards.length || DEFAULT_SHARD_COUNT
  });
  writeFileSync(SHARDS_FILE, `${JSON.stringify({ minutes: plan.minutes, shards: plan.shards }, null, 2)}\n`);
  for (const [index, files] of plan.shards.entries()) {
    const estimate = (plan.shardMinutes[index] / WORKERS_PER_SHARD).toFixed(1);
    console.log(`shard ${index + 1}: ${files.length} files, about ${estimate} min`);
  }
  if (plan.unmeasured.length) console.log(`estimated from case counts: ${plan.unmeasured.join(', ')}`);
};

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  main(process.argv.slice(2));
}
