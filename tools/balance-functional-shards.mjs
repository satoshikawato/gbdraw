#!/usr/bin/env node
// Plans the functional Playwright CI shards from measured case durations
// (docs/internal/SELECTIVE_CI.md "Functional shards").
//
//   node tools/balance-functional-shards.mjs plan
//   node tools/balance-functional-shards.mjs refresh <functional-report.json ...>
//
// tests/ci/functional-shards.json records the shard count and, for each spec file, the
// measured seconds per case. `plan` lists the functional cases of the working tree, weighs
// each spec file by its current case count times its seconds per case (an unmeasured file
// by the mean of the measured cases), and assigns the files longest first to the shard with
// the least work. It prints the plan and, in CI, writes the lists as the `shards` step
// output, which every functional shard job of the run reads. `refresh` rewrites the
// measurements from the test-results/functional-report.json of every shard of one or more
// green `Tests` runs. tools/ci-inventory.mjs imports the same functions for its weekly report.
import { spawnSync } from 'node:child_process';
import { appendFileSync, readFileSync, writeFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import path from 'node:path';
import { pathToFileURL } from 'node:url';

export const SHARDS_FILE = 'tests/ci/functional-shards.json';
export const WORKERS_PER_SHARD = 2;

const byText = (a, b) => (a < b ? -1 : a > b ? 1 : 0);
const sum = (values) => values.reduce((total, value) => total + value, 0);
const round2 = (value) => Math.round(value * 100) / 100;
const median = (values) => {
  const sorted = [...values].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
};

/**
 * The measurements in tests/ci/functional-shards.json.
 * @typedef {object} ShardMeasurements
 * @property {number} shardCount The number of functional shard jobs in the workflow matrix.
 * @property {Map<string, number>} secondsPerCase Measured seconds per case of each spec file,
 *   keyed by its path relative to the repository root.
 */

/**
 * Reads tests/ci/functional-shards.json.
 * @param {string} source The file's text.
 * @returns {ShardMeasurements}
 */
export const parseMeasurements = (source) => {
  const data = JSON.parse(source);
  if (!Number.isInteger(data.shardCount) || data.shardCount < 1) {
    throw new Error(`${SHARDS_FILE}: shardCount must be a positive integer`);
  }
  const entries = Object.entries(data.secondsPerCase ?? {});
  if (!entries.length || entries.some(([, seconds]) => typeof seconds !== 'number' || !(seconds >= 0))) {
    throw new Error(`${SHARDS_FILE}: secondsPerCase must map spec files to non-negative seconds`);
  }
  return { shardCount: data.shardCount, secondsPerCase: new Map(entries) };
};

export const readMeasurements = (file = SHARDS_FILE) => parseMeasurements(readFileSync(file, 'utf8'));

// Playwright reports spec files relative to its test directory (`config.rootDir`). The
// measurements and the shard lists name them relative to the repository root, which is
// also how the CLI file filters in the workflow name them. `testDir` is that test
// directory relative to the repository root; a CI report's own rootDir is a path on the runner.
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
  if (report.errors?.length) throw new Error(JSON.stringify(report.errors));
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

// Each green run reports a case once, retries included. A case's minutes are the median
// of its minutes over all given reports; the result lists them by spec file.
const caseMinutesByFile = (reports, testDir) => {
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
  const byFile = new Map();
  for (const [key, values] of caseMinutes) {
    const [file] = JSON.parse(key);
    byFile.set(file, [...(byFile.get(file) || []), median(values)]);
  }
  return byFile;
};

// A file's minutes are the sum of the median minutes of its cases.
export const measureMinutes = (reports, testDir) => new Map([...caseMinutesByFile(reports, testDir)]
  .map(([file, minutes]) => [file, sum(minutes)]));

// A file's seconds per case are the mean of the median seconds of its cases.
export const measureSecondsPerCase = (reports, testDir) => new Map([...caseMinutesByFile(reports, testDir)]
  .map(([file, minutes]) => [file, round2((sum(minutes) * 60) / minutes.length)]));

// Weighs each listed file by its case count times its seconds per case (an unmeasured file
// by the mean seconds of the measured cases) and assigns the files longest first to the
// shard with the least work. Returns the estimated minutes of each file, the sorted file
// list and the summed test minutes of each shard, and the files without a measurement.
export const planShards = ({ casesPerFile, secondsPerCase, shardCount }) => {
  const measuredFiles = [...casesPerFile.keys()].filter((file) => secondsPerCase.has(file));
  const meanSeconds = sum(measuredFiles.map((file) => casesPerFile.get(file) * secondsPerCase.get(file)))
    / sum(measuredFiles.map((file) => casesPerFile.get(file))) || 60;
  const minutes = Object.fromEntries([...casesPerFile]
    .map(([file, cases]) => [file, round2((cases * (secondsPerCase.get(file) ?? meanSeconds)) / 60)])
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
    shardMinutes: shards.map((shard) => round2(shard.minutes)),
    unmeasured: [...casesPerFile.keys()].filter((file) => !secondsPerCase.has(file)).sort(byText)
  };
};

const formatMeasurements = ({ shardCount, secondsPerCase }) => `${JSON.stringify({
  shardCount,
  secondsPerCase: Object.fromEntries([...secondsPerCase].sort(([a], [b]) => byText(a, b)))
}, null, 2)}\n`;

const planCommand = () => {
  const { shardCount, secondsPerCase } = readMeasurements();
  const { casesPerFile } = listFunctionalCases();
  const plan = planShards({ casesPerFile, secondsPerCase, shardCount });
  for (const [index, files] of plan.shards.entries()) {
    const estimate = (plan.shardMinutes[index] / WORKERS_PER_SHARD).toFixed(1);
    console.log(`shard ${index + 1}: ${files.length} files, about ${estimate} min`);
  }
  if (plan.unmeasured.length) console.log(`estimated from the mean case: ${plan.unmeasured.join(', ')}`);
  if (process.env.GITHUB_OUTPUT) {
    appendFileSync(process.env.GITHUB_OUTPUT, `shards=${JSON.stringify(plan.shards)}\n`, 'utf8');
  }
};

const refreshCommand = (reportFiles) => {
  if (!reportFiles.length) throw new Error('refresh needs the functional-report.json of each shard');
  const { shardCount } = readMeasurements();
  const { testDir } = listFunctionalCases();
  const secondsPerCase = measureSecondsPerCase(
    reportFiles.map((file) => JSON.parse(readFileSync(file, 'utf8'))), testDir
  );
  writeFileSync(SHARDS_FILE, formatMeasurements({ shardCount, secondsPerCase }));
  planCommand();
};

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  const [command, ...rest] = process.argv.slice(2);
  if (command === 'plan' && !rest.length) planCommand();
  else if (command === 'refresh') refreshCommand(rest);
  else {
    console.error('Usage: balance-functional-shards.mjs plan | refresh <functional-report.json ...>');
    process.exitCode = 2;
  }
}
