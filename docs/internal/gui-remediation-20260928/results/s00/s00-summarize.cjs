// Summarize S00 harness samples: node s00-summarize.cjs <samples.jsonl>... > summary.json
// Warm-up samples (declared before measurement) are the only exclusion.
'use strict';

const { readFileSync } = require('node:fs');

const BUDGET = Object.freeze({ visibleP95Ms: 100, visibleMaxMs: 250, compareSettleP95Ms: 200, taskMs: 100 });

const quantile = (values, q) => {
  if (!values.length) return null;
  const sorted = [...values].sort((a, b) => a - b);
  const position = (sorted.length - 1) * q;
  const low = Math.floor(position);
  const high = Math.ceil(position);
  return sorted[low] + (sorted[high] - sorted[low]) * (position - low);
};
const round = (value) => (value === null ? null : Math.round(value * 10) / 10);
const stats = (values) => ({
  n: values.length,
  p50: round(quantile(values, 0.5)),
  p95: round(quantile(values, 0.95)),
  max: round(values.length ? Math.max(...values) : null)
});

const samples = process.argv.slice(2).flatMap((path) => readFileSync(path, 'utf8')
  .split('\n').filter(Boolean).map((line) => JSON.parse(line)))
  .filter((record) => record.kind === 'sample');

const groups = new Map();
// Coverage-pass samples (with call counts) have distorted timing; keep them out.
for (const sample of samples) {
  if (sample.warmup || sample.counts) continue;
  const worker = sample.diagramWorker || {};
  const temperature = sample.op !== 'generate' ? 'n/a'
    : worker.constructionsDelta > 0 ? 'cold' : 'warm';
  // Insert and delete are the same ordinary keystroke; pool them per op.
  const direction = sample.op.startsWith('input-') ? 'any' : sample.direction || sample.sequence || '';
  const key = [sample.target, sample.state, sample.op, direction, temperature].join('|');
  if (!groups.has(key)) groups.set(key, []);
  groups.get(key).push(sample);
}

const rows = [...groups.entries()].map(([key, group]) => {
  const [target, state, op, direction, temperature] = key.split('|');
  const values = (field) => group.map((sample) => sample[field]).filter(Number.isFinite);
  const visible = stats(values('visibleMs'));
  const settle = stats(values('settleMs'));
  const posts = group.map((sample) => sample.diagramWorker || {});
  const row = {
    target, state, op, direction, temperature,
    pages: new Set(group.map((sample) => sample.pageIndex)).size,
    visibleMs: visible,
    settleMs: settle,
    generateTotalMs: op === 'generate' ? stats(group.map((sample) => sample.terminal?.clearedMs).filter(Number.isFinite)) : undefined,
    maxLoafMs: round(Math.max(0, ...values('maxLoafMs'))),
    samplesWithTaskAtLeast100Ms: group.filter((sample) => sample.loafAtLeast100 > 0).length,
    unsettled: group.filter((sample) => !sample.settled).length,
    diagramWorkerRuns: posts.reduce((total, worker) => total + (worker.runsDelta || 0), 0),
    diagramWorkerHelpers: posts.reduce((total, worker) => total + (worker.helpersDelta || 0), 0),
    diagramWorkerConstructions: posts.reduce((total, worker) => total + (worker.constructionsDelta || 0), 0)
  };
  if (op !== 'generate') {
    const failures = [];
    if (visible.p95 > BUDGET.visibleP95Ms) failures.push(`visible p95 ${visible.p95} > ${BUDGET.visibleP95Ms}`);
    if (visible.max > BUDGET.visibleMaxMs) failures.push(`visible max ${visible.max} > ${BUDGET.visibleMaxMs}`);
    if (op === 'compare' && settle.p95 > BUDGET.compareSettleP95Ms) failures.push(`settle p95 ${settle.p95} > ${BUDGET.compareSettleP95Ms}`);
    if (row.samplesWithTaskAtLeast100Ms > 0) failures.push(`${row.samplesWithTaskAtLeast100Ms} samples with a >=${BUDGET.taskMs} ms frame`);
    if (row.diagramWorkerRuns + row.diagramWorkerHelpers > 0) failures.push('diagram Worker dispatch during interaction');
    row.proposedBudget = failures.length ? { verdict: 'FAIL', failures } : { verdict: 'PASS' };
    row.sampleCountMeetsPlan = visible.n >= 20;
  }
  return row;
}).sort((a, b) => [a.target, a.state, a.op, a.direction].join().localeCompare([b.target, b.state, b.op, b.direction].join()));

// Function call counts from the precise-coverage run, per op/state (sum over samples / sample count).
const coverage = {};
for (const sample of samples) {
  if (!sample.counts) continue;
  const key = [sample.target, sample.state, sample.op, sample.direction || sample.sequence || ''].join('|');
  const entry = coverage[key] ||= { samples: 0, perSample: {} };
  entry.samples += 1;
  for (const [name, count] of Object.entries(sample.counts)) {
    entry.perSample[name] = (entry.perSample[name] || 0) + count;
  }
}
for (const entry of Object.values(coverage)) {
  for (const name of Object.keys(entry.perSample)) entry.perSample[name] = round(entry.perSample[name] / entry.samples);
}

process.stdout.write(`${JSON.stringify({ budget: BUDGET, rows, coverage }, null, 1)}\n`);
