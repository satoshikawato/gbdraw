#!/usr/bin/env node
// Weekly CI inventory (docs/internal/SELECTIVE_CI.md "Weekly CI inventory").
//
//   node tools/ci-inventory.mjs [--post]
//
// Reports the last 7 days of GitHub Actions use: runner-minutes by workflow and event
// and by job, the wall and queue time of `Tests` runs, the pull request feedback time
// against its targets, and the functional Playwright shards of the latest green dev
// `Tests` run with the plan that tools/balance-functional-shards.mjs would compute from
// that run's measurements, and the CI of pull requests: how many merged ones passed on one
// run, runs per merged one, and how many of their runs failed or were replaced by a newer
// run. It prints Markdown; with --post it
// also adds the report as this week's comment on the open issue labelled
// `ci-inventory` (creating the label and issue when missing). Week-over-week changes
// come from the data block of the workflow's own previous comment. Job lists are read
// only for `Tests` and `Gallery publication`, and only while RATE_RESERVE requests of the
// hourly budget are left; other runs count by run duration.
//
// Reads GITHUB_TOKEN (or GH_TOKEN), GITHUB_REPOSITORY, GITHUB_API_URL, and
// GITHUB_SERVER_URL. Run it from the repository root: it lists the functional cases
// with @playwright/test and reads the measurements in tests/ci/functional-shards.json and
// .github/workflows/test.yml.
import { spawnSync } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

import {
  WORKERS_PER_SHARD, casesPerFileOf, listFunctionalCases, measureMinutes, measureSecondsPerCase, planShards,
  readMeasurements
} from './balance-functional-shards.mjs';

const LABEL = 'ci-inventory';
// Only comments of the workflow's own token are reports; anyone can comment on the issue.
const REPORT_AUTHOR = 'github-actions[bot]';
const TESTS_WORKFLOW = '.github/workflows/test.yml';
// Workflows whose job lists are read. The others count by run duration: one job list per
// run would need several times the hourly budget of a GITHUB_TOKEN.
const JOB_LIST_WORKFLOWS = new Set([TESTS_WORKFLOW, '.github/workflows/gallery-publication.yml']);
const RUN_DURATION = ' (run duration)';
const DAY_MS = 24 * 60 * 60 * 1000;
const REFRESH_RATIO = 1.2;
const TOP = 15;
// The job lists stop while this many requests of the hourly budget are left. A
// GITHUB_TOKEN has 1,000 requests an hour for the whole repository, shared by every
// workflow (`CI impact plan`, `PR / gate`, and the `Fail fast` jobs among them), so the
// inventory leaves them half; its own issue and comment calls after the job lists use a few.
export const RATE_RESERVE = 500;
const PARALLEL_REQUESTS = 8;
const DATA_BLOCK = /<!-- ci-inventory-data (\{.*?\}) -->/s;
const WEEK = /^\d{4}-W\d{2}$/;
const DATA_NUMBERS = ['cases', 'functionalTestMinutes', 'runnerMinutes', 'testJobs'];
// Pull request CI shares: null for a week without them, and for reports from before them.
const DATA_SHARES = ['onePassPercent', 'runsPerPullRequest', 'supersededPercent', 'failedPercent'];
// Pull request feedback targets (docs/internal/SELECTIVE_CI.md, "Validation tiers"): the
// median wall time, in minutes, of the head `Tests` runs of the last PR_WINDOW pull
// requests merged into dev, by whether a functional Playwright shard ran.
export const PR_WINDOW = 20;
export const PR_TARGETS = [
  { kind: 'without functional Playwright', functional: false, minutes: 20 },
  { kind: 'with functional Playwright', functional: true, minutes: 35 }
];
const FUNCTIONAL_JOB = 'Playwright functional';

const sum = (values) => values.reduce((total, value) => total + value, 0);
const round1 = (value) => Math.round(value * 10) / 10;
const minutesBetween = (from, to) => (from && to ? (Date.parse(to) - Date.parse(from)) / 60000 : 0);
// Nearest-rank percentile; null for no values.
export const percentile = (values, p) => {
  if (!values.length) return null;
  const sorted = [...values].sort((a, b) => a - b);
  return sorted[Math.max(0, Math.ceil((p / 100) * sorted.length) - 1)];
};

// ISO 8601 week of a date, such as 2026-W42.
export const isoWeek = (date) => {
  const day = new Date(Date.UTC(date.getUTCFullYear(), date.getUTCMonth(), date.getUTCDate()));
  day.setUTCDate(day.getUTCDate() + 4 - (day.getUTCDay() || 7));
  const yearStart = Date.UTC(day.getUTCFullYear(), 0, 1);
  const week = Math.ceil(((day - yearStart) / DAY_MS + 1) / 7);
  return `${day.getUTCFullYear()}-W${String(week).padStart(2, '0')}`;
};

// CodeQL default setup runs as `dynamic/github-code-scanning/codeql` under PR-specific names.
const workflowLabel = (run) => (run.path?.startsWith('dynamic/') ? run.path.split('/').pop() : run.name);
// Matrix jobs share the name before their parameters, such as `Playwright functional`.
const jobLabel = (name) => name.replace(/ \(.*\)$/, '');

// Runner-minutes and Tests timing of completed runs. `jobsByRun` maps a run id to its jobs
// of every attempt; a run without a job list counts its own duration under its workflow
// name marked "(run duration)", in both tables.
export const summarizeRuns = ({ runs, jobsByRun }) => {
  const byWorkflow = new Map();
  const byJob = new Map();
  const tests = { pull_request: { walls: [], queues: [] }, 'dev push': { walls: [], queues: [] } };
  let runDurationOnly = 0;
  const add = (map, key, minutes, runs = 0) => {
    const entry = map.get(key) || { runs: 0, minutes: 0 };
    map.set(key, { runs: entry.runs + runs, minutes: entry.minutes + minutes });
  };
  for (const run of runs) {
    const workflow = `${workflowLabel(run)} [${run.event}]`;
    const jobs = jobsByRun.get(run.id);
    if (!jobs) {
      runDurationOnly += 1;
      const minutes = minutesBetween(run.run_started_at, run.updated_at);
      add(byWorkflow, `${workflow}${RUN_DURATION}`, minutes, 1);
      add(byJob, `${workflowLabel(run)}${RUN_DURATION}`, minutes);
      continue;
    }
    add(byWorkflow, workflow, 0, 1);
    for (const job of jobs) {
      const minutes = minutesBetween(job.started_at, job.completed_at);
      add(byWorkflow, workflow, minutes);
      add(byJob, jobLabel(job.name), minutes);
    }
    const kind = run.event === 'pull_request' ? 'pull_request'
      : run.event === 'push' && run.head_branch === 'dev' ? 'dev push' : null;
    if (run.path !== TESTS_WORKFLOW || !kind) continue;
    const ran = jobs.filter((job) => job.run_attempt === run.run_attempt
      && job.conclusion !== 'skipped' && job.started_at && job.completed_at);
    if (!ran.length) continue;
    tests[kind].walls.push(minutesBetween(
      run.run_started_at, ran.map((job) => job.completed_at).sort().pop()
    ));
    tests[kind].queues.push(...ran.map((job) => minutesBetween(job.created_at, job.started_at)));
  }
  const ranked = (map) => [...map].map(([name, entry]) => ({ name, ...entry, minutes: round1(entry.minutes) }))
    .sort((a, b) => b.minutes - a.minutes || (a.name < b.name ? -1 : 1));
  const timing = (kind) => ({
    runs: tests[kind].walls.length,
    wall: [50, 90].map((p) => percentile(tests[kind].walls, p)),
    queue: [50, 90].map((p) => percentile(tests[kind].queues, p))
  });
  return {
    runs: runs.length,
    runDurationOnly,
    runnerMinutes: round1(sum([...byWorkflow.values()].map((entry) => entry.minutes))),
    byWorkflow: ranked(byWorkflow),
    byJob: ranked(byJob),
    tests: { pull_request: timing('pull_request'), 'dev push': timing('dev push') }
  };
};

// The pull request feedback time against PR_TARGETS. `jobLists` holds, for each pull
// request of the window, the jobs of its latest green head `Tests` run, or null when it
// has none. A run's wall time is its first job start to its last job end, so a re-run
// of some jobs counts until the re-run finished. The median is the nearest-rank p50.
export const summarizePullRequests = (jobLists) => {
  const walls = jobLists.flatMap((jobs) => {
    const ran = (jobs || []).filter((job) => job.conclusion !== 'skipped' && job.started_at && job.completed_at);
    if (!ran.length) return [];
    return [{
      functional: ran.some((job) => jobLabel(job.name) === FUNCTIONAL_JOB),
      minutes: minutesBetween(ran.map((job) => job.started_at).sort()[0], ran.map((job) => job.completed_at).sort().pop())
    }];
  });
  return {
    pullRequests: jobLists.length,
    measured: walls.length,
    kinds: PR_TARGETS.map((target) => {
      const minutes = walls.filter((wall) => wall.functional === target.functional).map((wall) => wall.minutes);
      const median = percentile(minutes, 50);
      // One rounded median decides the warning and is the number the report prints.
      return { ...target, runs: minutes.length, median, over: median !== null && round1(median) > target.minutes };
    })
  };
};

/**
 * A workflow run as the REST API lists it; only these fields are read.
 * @typedef {object} WorkflowRun
 * @property {number} id
 * @property {string} path The workflow file, such as `.github/workflows/test.yml`.
 * @property {string} event
 * @property {string} head_branch
 * @property {{ full_name: string } | null} head_repository Null when the head repository was deleted.
 * @property {string} created_at
 * @property {string} updated_at The end of a completed run.
 * @property {number} run_attempt
 * @property {string | null} conclusion
 */
/**
 * A job of one attempt of a run.
 * @typedef {object} WorkflowJob
 * @property {string} name
 * @property {number} run_attempt
 * @property {string | null} conclusion
 */
/**
 * A pull request merged into dev in the window.
 * @typedef {object} MergedPullRequest
 * @property {string} created_at
 * @property {string} merged_at
 * @property {{ ref: string, repo: { full_name: string } | null }} head
 */
/**
 * Pull request CI of the window (docs/internal/SELECTIVE_CI.md, "Weekly CI inventory").
 * The four shares are rounded to one decimal, or null when nothing was counted.
 * @typedef {object} PullRequestCi
 * @property {number} mergedPullRequests Pull requests merged into dev in the window.
 * @property {number} withRuns Those with at least one pull request `Tests` run.
 * @property {number} singlePass Those whose only run passed at its first attempt.
 * @property {number} pullRequestRuns Pull request `Tests` runs created in the window.
 * @property {Record<RunOutcome, number>} runs Those runs by outcome.
 * @property {number} failFast Failed runs that GitHub reports as cancelled (fail-fast).
 * @property {number | null} onePassPercent
 * @property {number | null} runsPerPullRequest
 * @property {number | null} supersededPercent
 * @property {number | null} failedPercent
 */
/** @typedef {'passed' | 'failed' | 'superseded' | 'other'} RunOutcome */

// A run joins a pull request when it ran on its head branch between its opening and its
// merge, give or take this slack: a merged pull request's runs list no pull requests.
const JOIN_SLACK_MS = 5 * 60 * 1000;
// Without a job list, a cancelled run followed by the next run on its branch within this
// time was replaced by it (a push, or a label change); a longer gap means a
// `Fail fast (...)` job cancelled it after a failed job.
export const SUPERSEDED_GAP_SECONDS = 120;
const FAILED = new Set(['failure', 'timed_out', 'startup_failure']);
const GATE_JOB = 'PR / gate';
const branchOf = (/** @type {WorkflowRun} */ run) => `${run.head_repository?.full_name ?? ''}:${run.head_branch}`;
/** @type {(count: number, total: number) => number | null} */
const share = (count, total) => (total ? round1((100 * count) / total) : null);

/**
 * The outcome of a pull request run. A cancelled run with a job list failed when a job of
 * its attempt other than `PR / gate` failed; otherwise a newer run replaced it or someone
 * cancelled it ('superseded'). Without a job list, the gap to the next run on its branch
 * decides (SUPERSEDED_GAP_SECONDS), and a run with no later run counts as replaced.
 * @param {WorkflowRun} run
 * @param {WorkflowRun | undefined} next The next run on the same branch.
 * @param {WorkflowJob[] | undefined} jobs
 * @returns {RunOutcome}
 */
const runOutcome = (run, next, jobs) => {
  if (run.conclusion === 'success') return 'passed';
  if (FAILED.has(run.conclusion ?? '')) return 'failed';
  if (run.conclusion !== 'cancelled') return 'other';
  if (jobs) {
    return jobs.some((job) => job.run_attempt === run.run_attempt && job.name !== GATE_JOB
      && FAILED.has(job.conclusion ?? '')) ? 'failed' : 'superseded';
  }
  return next && (Date.parse(next.created_at) - Date.parse(run.updated_at)) / 1000 > SUPERSEDED_GAP_SECONDS
    ? 'failed' : 'superseded';
};

/**
 * Pull request CI of the window. `runs` are the completed runs created in the window,
 * `earlierRuns` the pull request `Tests` runs before it of pull requests opened before it,
 * and `pulls` the pull requests merged into dev in the window.
 * @param {{ runs: WorkflowRun[], earlierRuns: WorkflowRun[], jobsByRun: Map<number, WorkflowJob[]>,
 *   pulls: MergedPullRequest[] }} input
 * @returns {PullRequestCi}
 */
export const summarizePullRequestCi = ({ runs, earlierRuns, jobsByRun, pulls }) => {
  const isPullRequestRun = (/** @type {WorkflowRun} */ run) => run.path === TESTS_WORKFLOW && run.event === 'pull_request';
  const windowRuns = runs.filter(isPullRequestRun);
  /** @type {Map<string, WorkflowRun[]>} */
  const byBranch = new Map();
  for (const run of [...earlierRuns.filter(isPullRequestRun), ...windowRuns]) {
    const list = byBranch.get(branchOf(run)) || [];
    if (!list.some((other) => other.id === run.id)) list.push(run);
    byBranch.set(branchOf(run), list);
  }
  for (const list of byBranch.values()) list.sort((a, b) => Date.parse(a.created_at) - Date.parse(b.created_at));

  /** @type {Record<RunOutcome, number>} */
  const outcomes = { passed: 0, failed: 0, superseded: 0, other: 0 };
  let failFast = 0;
  for (const run of windowRuns) {
    const list = byBranch.get(branchOf(run)) || [];
    const outcome = runOutcome(run, list[list.findIndex((other) => other.id === run.id) + 1], jobsByRun.get(run.id));
    outcomes[outcome] += 1;
    if (outcome === 'failed' && run.conclusion === 'cancelled') failFast += 1;
  }

  const runsPerPull = pulls.map((pull) => {
    const from = Date.parse(pull.created_at) - JOIN_SLACK_MS;
    const to = Date.parse(pull.merged_at) + JOIN_SLACK_MS;
    return (byBranch.get(`${pull.head.repo?.full_name ?? ''}:${pull.head.ref}`) || [])
      .filter((run) => Date.parse(run.created_at) >= from && Date.parse(run.created_at) <= to);
  }).filter((list) => list.length);
  const singlePass = runsPerPull.filter((list) => list.length === 1
    && list[0].conclusion === 'success' && list[0].run_attempt === 1).length;
  return {
    mergedPullRequests: pulls.length,
    withRuns: runsPerPull.length,
    singlePass,
    pullRequestRuns: windowRuns.length,
    runs: outcomes,
    failFast,
    onePassPercent: share(singlePass, runsPerPull.length),
    runsPerPullRequest: runsPerPull.length ? round1(sum(runsPerPull.map((list) => list.length)) / runsPerPull.length) : null,
    supersededPercent: share(outcomes.superseded, windowRuns.length),
    failedPercent: share(outcomes.failed, windowRuns.length)
  };
};

const SHARD_ARTIFACT = /^playwright-functional-shard-(\d+)-traces-\d+-(\d+)$/;

// The shard artifacts of a run, one per shard: a re-run keeps the earlier attempt's
// artifacts, so only the highest attempt of each shard counts.
export const latestShardArtifacts = (artifacts) => {
  const latest = new Map();
  for (const artifact of artifacts) {
    const match = artifact.name.match(SHARD_ARTIFACT);
    if (!match || artifact.expired) continue;
    const [shard, attempt] = [Number(match[1]), Number(match[2])];
    if (!(latest.get(shard)?.attempt >= attempt)) latest.set(shard, { shard, attempt, artifact });
  }
  return [...latest.values()].sort((a, b) => a.shard - b.shard);
};

// The functional shards of one green run. `reports` are the shard reports with their
// shard numbers (a shard without a report is listed as missing); `measurements` are the
// recorded ones (readMeasurements); `casesPerFile` and `testDir` come from listing the
// functional suite. CI plans the shards from the recorded measurements and the current case
// counts, so an uneven run means that seconds per case moved: refresh the measurements.
export const summarizeFunctional = ({ reports, measurements, casesPerFile, testDir }) => {
  const shardTestMinutes = (report) => sum([...measureMinutes([report], testDir).values()]);
  const shards = reports
    .map(({ shard, report }) => ({
      shard,
      testMinutes: round1(shardTestMinutes(report)),
      wallMinutes: round1((report.stats?.duration || 0) / 60000)
    }))
    .sort((a, b) => a.shard - b.shard);
  const fileMinutes = measureMinutes(reports.map(({ report }) => report), testDir);
  const casesRun = casesPerFileOf(reports.map(({ report }) => report), testDir);
  const { shardCount } = measurements;
  const current = planShards({ casesPerFile, secondsPerCase: measurements.secondsPerCase, shardCount });
  const plan = planShards({
    casesPerFile, secondsPerCase: measureSecondsPerCase(reports.map(({ report }) => report), testDir), shardCount
  });
  const testMinutes = sum(shards.map((shard) => shard.testMinutes));
  const mean = testMinutes / (shards.length || 1);
  const slowest = Math.max(0, ...shards.map((shard) => shard.testMinutes));
  // One rounded ratio decides the flag and is the number the report prints.
  const slowestToMean = mean ? Math.round((slowest / mean) * 100) / 100 : 0;
  return {
    shards,
    missingShards: Array.from({ length: shardCount }, (_, index) => index + 1)
      .filter((number) => !shards.some((shard) => shard.shard === number)),
    cases: sum([...casesRun.values()]),
    testMinutes: round1(testMinutes),
    slowestToMean,
    unmeasured: current.unmeasured,
    refresh: slowestToMean > REFRESH_RATIO || current.unmeasured.length > 0,
    planChanges: JSON.stringify(plan.shards) !== JSON.stringify(current.shards),
    plannedTestMinutes: plan.shardMinutes,
    topSpecs: [...fileMinutes].map(([file, minutes]) => ({ file, cases: casesRun.get(file) || 0, minutes: round1(minutes) }))
      .sort((a, b) => b.minutes - a.minutes || (a.file < b.file ? -1 : 1))
      .slice(0, TOP)
  };
};

// Top-level jobs of a workflow file.
export const countWorkflowJobs = (yaml) => {
  const body = yaml.split(/^jobs:\s*$/m)[1] || '';
  return body.split('\n').filter((line) => /^ {2}[A-Za-z0-9_-]+:\s*$/.test(line)).length;
};

// The numbers of a report's data block, or null unless every field has its type; a missing
// share reads as null.
export const readData = (body) => {
  const match = typeof body === 'string' ? body.match(DATA_BLOCK) : null;
  if (!match) return null;
  let raw;
  try {
    raw = JSON.parse(match[1]);
  } catch {
    return null;
  }
  if (typeof raw?.week !== 'string' || !WEEK.test(raw.week)) return null;
  const data = { week: raw.week };
  for (const field of DATA_NUMBERS) {
    if (typeof raw[field] !== 'number' || !Number.isFinite(raw[field])) return null;
    data[field] = raw[field];
  }
  for (const field of DATA_SHARES) {
    const value = raw[field] ?? null;
    if (value !== null && (typeof value !== 'number' || !Number.isFinite(value))) return null;
    data[field] = value;
  }
  return data;
};

// This week's earlier report comment, to replace, and the previous week's numbers, from
// the reports the workflow itself posted.
export const priorReports = (comments, week) => {
  const reports = comments
    .filter((comment) => comment.user?.login === REPORT_AUTHOR)
    .map((comment) => ({ comment, data: readData(comment.body) }))
    .filter(({ data }) => data);
  const thisWeek = reports.at(-1)?.data.week === week ? reports.at(-1).comment : null;
  return { thisWeek, previous: (thisWeek ? reports.at(-2) : reports.at(-1))?.data ?? null };
};

const SHARE_ROWS = [
  ['onePassPercent', 'Merged pull requests that passed on one run (%)'],
  ['runsPerPullRequest', '`Tests` runs per merged pull request'],
  ['supersededPercent', 'Pull request runs replaced by a newer run or cancelled by hand (%)'],
  ['failedPercent', 'Pull request runs that failed (%)']
];
const RUN_OUTCOMES = { passed: 'passed', failed: 'failed', superseded: 'replaced by a newer run or cancelled by hand', other: 'other' };
const fmt = (value) => (value === null || value === undefined ? '-' : String(round1(value)));
const delta = (now, before) => {
  if (before === null || before === undefined) return '-';
  const change = round1(now - before);
  return change > 0 ? `+${change}` : String(change);
};
const table = (header, rows) => [
  `| ${header.join(' | ')} |`,
  `| ${header.map(() => '---').join(' | ')} |`,
  ...rows.map((row) => `| ${row.join(' | ')} |`)
].join('\n');

// The weekly report as Markdown, ending with its data block.
export const renderReport = ({
  week, since, until, runs, partial = null, pullRequests, pullRequestCi, functional, functionalRun, testJobs, previous,
  serverUrl, repository
}) => {
  const data = {
    week,
    cases: functional.cases,
    functionalTestMinutes: functional.testMinutes,
    runnerMinutes: runs.runnerMinutes,
    testJobs,
    ...Object.fromEntries(DATA_SHARES.map((field) => [field, pullRequestCi[field]]))
  };
  const reasons = [
    functional.slowestToMean > REFRESH_RATIO
      ? `slowest shard ${functional.slowestToMean}x the mean (limit ${REFRESH_RATIO}x)` : null,
    functional.unmeasured.length
      ? `${functional.unmeasured.length} spec file(s) without measurements: ${functional.unmeasured.map((file) => `\`${file}\``).join(', ')}`
      : null
  ].filter(Boolean);
  const overTarget = pullRequests.kinds.filter((kind) => kind.over);
  const runLink =`[${functionalRun.id}](${serverUrl}/${repository}/actions/runs/${functionalRun.id})`;
  const perWorker = (minutes) => fmt(minutes / WORKERS_PER_SHARD);
  return [
    `## CI inventory ${week}`,
    '',
    `Window: ${since.toISOString().slice(0, 16)}Z to ${until.toISOString().slice(0, 16)}Z.`
      + ` What to do with a flagged report: docs/internal/SELECTIVE_CI.md, "Weekly CI inventory".`,
    ...(partial ? ['', `**Partial report**: ${partial}`] : []),
    '',
    functional.refresh
      ? `**Refresh the functional shard measurements**: ${reasons.join('; ')}.`
      : `Functional shards are balanced (slowest ${functional.slowestToMean}x the mean).`,
    '',
    overTarget.length
      ? `**Pull request feedback over target**: ${overTarget.map((kind) =>
        `median ${fmt(kind.median)} minutes ${kind.kind} (target ${kind.minutes})`).join('; ')}.`
      : `Pull request feedback is within target (${pullRequests.kinds.map((kind) =>
        `median ${fmt(kind.median)} of ${kind.minutes} minutes ${kind.kind}`).join('; ')}).`,
    '',
    table(['Measure', 'This week', 'Previous', 'Change'], [
      ['Functional cases', data.cases, previous?.cases ?? '-', delta(data.cases, previous?.cases)],
      ['Functional test-minutes', fmt(data.functionalTestMinutes), fmt(previous?.functionalTestMinutes),
        delta(data.functionalTestMinutes, previous?.functionalTestMinutes)],
      ['Runner-minutes, all workflows', fmt(data.runnerMinutes), fmt(previous?.runnerMinutes),
        delta(data.runnerMinutes, previous?.runnerMinutes)],
      ['Jobs in `test.yml`', data.testJobs, previous?.testJobs ?? '-', delta(data.testJobs, previous?.testJobs)],
      ...SHARE_ROWS.map(([field, measure]) => [
        measure, fmt(data[field]), fmt(previous?.[field]),
        data[field] === null ? '-' : delta(data[field], previous?.[field])
      ])
    ]),
    previous ? `\nPrevious: ${previous.week}.` : '\nNo previous report.',
    '',
    '### Pull request CI',
    '',
    `${pullRequestCi.mergedPullRequests} pull requests merged into \`dev\` in the window,`
      + ` ${pullRequestCi.withRuns} with a \`Tests\` run; ${pullRequestCi.singlePass} passed on their only run at its first attempt.`
      + ` ${pullRequestCi.pullRequestRuns} pull request \`Tests\` runs were created in the window; ${pullRequestCi.failFast}`
      + ' of the failed ones GitHub reports as cancelled, because a `Fail fast (...)` job stopped them after a failed job.'
      + ' Definitions: docs/internal/SELECTIVE_CI.md, "Weekly CI inventory".',
    '',
    table(['Pull request runs', 'Count', 'Share (%)'], Object.entries(RUN_OUTCOMES).map(([outcome, caption]) => [
      caption, pullRequestCi.runs[outcome], fmt(share(pullRequestCi.runs[outcome], pullRequestCi.pullRequestRuns))
    ])),
    '',
    '### Runner-minutes by workflow and event',
    '',
    table(['Workflow [event]', 'Runs', 'Runner-minutes'],
      runs.byWorkflow.map((entry) => [entry.name, entry.runs, fmt(entry.minutes)])),
    '',
    `${runs.runs} completed runs. Rows marked "(run duration)" count ${runs.runDurationOnly} whole run(s), not jobs:`
      + ' job lists are read only for `Tests` and `Gallery publication`, and only while the API request budget lasts.',
    '',
    `### Runner-minutes by job (top ${TOP})`,
    '',
    table(['Job', 'Runner-minutes'], runs.byJob.slice(0, TOP).map((entry) => [entry.name, fmt(entry.minutes)])),
    '',
    '### `Tests` wall and queue time (minutes)',
    '',
    'Wall: run start to its last job. Queue: job created to job started.',
    '',
    table(['Runs', 'Count', 'Wall p50', 'Wall p90', 'Queue p50', 'Queue p90'],
      Object.entries(runs.tests).map(([kind, timing]) => [
        kind, timing.runs, ...timing.wall.map(fmt), ...timing.queue.map(fmt)
      ])),
    '',
    '### Pull request feedback time (minutes)',
    '',
    `The latest green head \`Tests\` run of each of the last ${pullRequests.pullRequests} pull requests merged into`
      + ` \`dev\` (${pullRequests.measured} had one), first job start to last job end.`
      + ' Targets: docs/internal/SELECTIVE_CI.md, "Validation tiers".',
    '',
    table(['Pull requests', 'Count', 'Median', 'Target', 'Over target'],
      pullRequests.kinds.map((kind) => [
        kind.kind, kind.runs, fmt(kind.median), kind.minutes, kind.over ? '**yes**' : 'no'
      ])),
    '',
    `### Functional shards of run ${runLink} (\`${functionalRun.head_sha.slice(0, 8)}\`)`,
    '',
    `Test-minutes add up case durations; a shard runs ${WORKERS_PER_SHARD} workers. `
      + (functional.planChanges
        ? 'The planned column is the plan that measurements refreshed from this run would give.'
        : 'Measurements refreshed from this run would keep the current plan.')
      + (functional.missingShards.length
        ? ` No functional report for shard ${functional.missingShards.join(', ')}; the totals leave it out.` : ''),
    '',
    table(['Shard', 'Test-minutes', 'Per worker', 'Reporter wall', 'Planned per worker'],
      functional.shards.map((shard) => [
        shard.shard, fmt(shard.testMinutes), perWorker(shard.testMinutes), fmt(shard.wallMinutes),
        perWorker(functional.plannedTestMinutes[shard.shard - 1] ?? 0)
      ])),
    '',
    `### Spec files by test-minutes (top ${TOP})`,
    '',
    table(['Spec file', 'Cases', 'Test-minutes'],
      functional.topSpecs.map((spec) => [`\`${spec.file}\``, spec.cases, fmt(spec.minutes)])),
    '',
    `<!-- ci-inventory-data ${JSON.stringify(data)} -->`,
    ''
  ].join('\n');
};

// --- GitHub REST -------------------------------------------------------------------

// GitHub answers some reads with a server error or a body cut short; one such answer would
// end the whole report, so a GET is tried up to READ_TRIES times. Writes are never repeated.
const READ_TRIES = 3;
const RETRY_WAIT_MS = 2000;

/**
 * A failed request. `transient` marks one that a later try of the same read may not
 * repeat: a server error or a body cut short.
 * @typedef {Error & { status?: number, transient?: boolean }} RequestError
 */
/** @type {(error: unknown) => RequestError} */
const transient = (error) => Object.assign(/** @type {Error} */ (error), { transient: true });
/** @type {(message: string, status: number) => RequestError} */
const httpError = (message, status) => Object.assign(new Error(message), { status, transient: status >= 500 });

export const createClient = ({ token, apiUrl, repository, retryWaitMs = RETRY_WAIT_MS }) => {
  let remaining = Infinity;
  let requests = 0;
  const request = async (path, { method = 'GET', body } = {}) => {
    requests += 1;
    const response = await fetch(path.startsWith('http') ? path : `${apiUrl}/repos/${repository}${path}`, {
      method,
      headers: {
        accept: 'application/vnd.github+json',
        authorization: `Bearer ${token}`,
        'x-github-api-version': '2022-11-28',
        ...(body ? { 'content-type': 'application/json' } : {})
      },
      body: body ? JSON.stringify(body) : undefined,
      redirect: 'manual'
    });
    const left = Number(response.headers.get('x-ratelimit-remaining'));
    if (Number.isFinite(left)) remaining = left;
    return response;
  };
  /**
   * Every read goes through here: `once` is tried again after a transient error.
   * @template T
   * @param {() => Promise<T>} once
   * @returns {Promise<T>}
   */
  const retried = async (once) => {
    for (let attempt = 1; ; attempt += 1) {
      try {
        return await once();
      } catch (error) {
        if (!(/** @type {RequestError} */ (error)).transient || attempt >= READ_TRIES) throw error;
      }
      await new Promise((resolve) => setTimeout(resolve, retryWaitMs * attempt));
    }
  };
  const json = async (path, options) => {
    const read = !options || options.method === undefined;
    const once = async () => {
      const response = await request(path, options);
      if (response.status === 404 && read) return null;
      if (!response.ok) {
        throw httpError(`${options?.method || 'GET'} ${path}: ${response.status} ${await response.text()}`, response.status);
      }
      try {
        return JSON.parse(await response.text());
      } catch (error) {
        throw transient(error);
      }
    };
    return read ? retried(once) : once();
  };
  /**
   * An artifact archive: the API answers with a redirect to a signed URL, read without the token.
   * @param {string} url
   * @returns {Promise<Buffer>}
   */
  const download = (url) => retried(async () => {
    const redirect = await request(url);
    await redirect.text();
    const location = redirect.headers.get('location');
    if (!location) throw httpError(`GET ${url}: ${redirect.status}`, redirect.status);
    const response = await fetch(location);
    if (!response.ok) {
      await response.text();
      throw httpError(`GET ${url} (redirected): ${response.status}`, response.status);
    }
    try {
      return Buffer.from(await response.arrayBuffer());
    } catch (error) {
      throw transient(error);
    }
  });
  // Pages until a short page; `key` names the list in each page, or null for a bare list.
  const pages = async (path, key) => {
    const items = [];
    for (let page = 1; ; page += 1) {
      const result = await json(`${path}${path.includes('?') ? '&' : '?'}per_page=100&page=${page}`);
      const list = key ? result?.[key] || [] : result || [];
      items.push(...list);
      if (list.length < 100) return items;
    }
  };
  return { json, download, pages, remaining: () => remaining, requests: () => requests };
};

// A rate-limit or permission refusal ends a listing early; the report says so.
const refused = (error) => [403, 429].includes(error.status);

// A filtered run listing returns at most this many runs.
const RUN_LISTING_LIMIT = 1000;

// The completed runs created from `start` to `end` (milliseconds, end excluded). A part
// with more runs than one listing returns is read in halves.
const listRunsCreated = async (client, start, end) => {
  const created = `${new Date(start).toISOString()}..${new Date(end - 1000).toISOString()}`;
  const path = `/actions/runs?status=completed&created=${encodeURIComponent(created)}`;
  const { total_count: total } = await client.json(`${path}&per_page=1`);
  if (total > RUN_LISTING_LIMIT && end - start > 2000) {
    const middle = start + Math.floor((end - start) / 2000) * 1000;
    return [...await listRunsCreated(client, start, middle), ...await listRunsCreated(client, middle, end)];
  }
  return total ? client.pages(path, 'workflow_runs') : [];
};

// The completed runs created in the window, a day at a time, and why the listing stopped
// early, if it did.
export const listRuns = async (client, since, until) => {
  const runs = [];
  let partial = null;
  for (let start = since.getTime(); start < until.getTime() && !partial; start += DAY_MS) {
    try {
      runs.push(...await listRunsCreated(client, start, Math.min(start + DAY_MS, until.getTime())));
    } catch (error) {
      if (!refused(error)) throw error;
      partial = `runs created from ${new Date(start).toISOString().slice(0, 16)}Z on are missing (HTTP ${error.status}).`;
    }
  }
  return { runs: [...new Map(runs.map((run) => [run.id, run])).values()], partial };
};

// Job lists of every attempt of `Tests` and `Gallery publication` runs, `Tests` first,
// while the request budget lasts.
export const listJobs = async (client, runs) => {
  const jobsByRun = new Map();
  const ordered = runs.filter((run) => JOB_LIST_WORKFLOWS.has(run.path))
    .sort((a, b) => Number(b.path === TESTS_WORKFLOW) - Number(a.path === TESTS_WORKFLOW));
  for (let start = 0; start < ordered.length && client.remaining() > RATE_RESERVE; start += PARALLEL_REQUESTS) {
    const batch = ordered.slice(start, start + PARALLEL_REQUESTS);
    let lists;
    try {
      lists = await Promise.all(batch.map((run) => client.pages(`/actions/runs/${run.id}/jobs?filter=all`, 'jobs')));
    } catch (error) {
      if (!refused(error)) throw error;
      break;
    }
    batch.forEach((run, index) => jobsByRun.set(run.id, lists[index]));
  }
  return jobsByRun;
};

// The pull requests merged into dev, last merged first: every one merged since `since`,
// and at least those among the 100 closed ones updated last. A merge updates a pull
// request, so the listing stops at a page that ends before `since`.
const listMergedPulls = async (client, since) => {
  const closed = [];
  for (let page = 1; ; page += 1) {
    const list = await client.json(`/pulls?state=closed&base=dev&sort=updated&direction=desc&per_page=100&page=${page}`);
    closed.push(...list);
    if (list.length < 100 || Date.parse(list.at(-1).updated_at) < since.getTime()) break;
  }
  return closed.filter((pull) => pull.merged_at).sort((a, b) => (a.merged_at < b.merged_at ? 1 : -1));
};

// The pull request `Tests` runs created before `since` of the pull requests opened before it.
const listEarlierRuns = async (client, pulls, since) => {
  const runs = [];
  for (const pull of pulls.filter((pull) => Date.parse(pull.created_at) < since.getTime())) {
    const from = new Date(Date.parse(pull.created_at) - JOIN_SLACK_MS).toISOString();
    const created = `${from}..${new Date(since.getTime() - 1000).toISOString()}`;
    runs.push(...await client.pages(`/actions/workflows/test.yml/runs?event=pull_request&status=completed`
      + `&branch=${encodeURIComponent(pull.head.ref)}&created=${encodeURIComponent(created)}`, 'workflow_runs'));
  }
  return runs;
};

// The jobs of the latest green head `Tests` run of each of the last PR_WINDOW pull
// requests merged into dev, or null for a pull request without one.
const listPullRequestRuns = async (client, mergedPulls) => {
  const jobLists = [];
  for (const pull of mergedPulls.slice(0, PR_WINDOW)) {
    const { workflow_runs: [run] } = await client.json(
      `/actions/workflows/test.yml/runs?event=pull_request&status=success&head_sha=${pull.head.sha}&per_page=1`
    );
    jobLists.push(run ? await client.pages(`/actions/runs/${run.id}/jobs?filter=latest`, 'jobs') : null);
  }
  return jobLists;
};

// The latest green dev push `Tests` run that kept an artifact for every shard.
const latestFunctionalRun = async (client, shardCount) => {
  const { workflow_runs: runs } = await client.json(
    `/actions/workflows/test.yml/runs?branch=dev&event=push&status=success&per_page=20`
  );
  for (const run of runs) {
    const shards = latestShardArtifacts(await client.pages(`/actions/runs/${run.id}/artifacts`, 'artifacts'));
    if (shards.length === shardCount) return { run, shards };
  }
  throw new Error('no green dev Tests run in the last 20 kept every functional shard artifact');
};

// Downloads each shard artifact and reads its functional-report.json; an artifact
// without one (a job that stopped before its report) is left out.
const readShardReports = async (client, shards) => {
  const dir = mkdtempSync(join(process.env.RUNNER_TEMP || tmpdir(), 'ci-inventory-'));
  try {
    const reports = [];
    for (const { shard, artifact } of shards) {
      const archive = join(dir, `${artifact.id}.zip`);
      writeFileSync(archive, await client.download(artifact.archive_download_url));
      const unzip = spawnSync('unzip', ['-p', archive, 'functional-report.json'], {
        encoding: 'utf8', maxBuffer: 256 * 1024 * 1024
      });
      // unzip exits 11 when the archive has no such file.
      if (unzip.status === 11) continue;
      if (unzip.status !== 0) throw new Error(`artifact ${artifact.name}: ${unzip.stderr}`);
      reports.push({ shard, report: JSON.parse(unzip.stdout) });
    }
    return reports;
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
};

const ISSUE_TITLE = 'Weekly CI inventory';
const ISSUE_BODY = 'The `CI inventory` workflow (`.github/workflows/ci-inventory.yml`) adds one comment'
  + ' a week. docs/internal/SELECTIVE_CI.md, "Weekly CI inventory", says what to do with a flagged report.';

const findIssue = async (client) => (await client.json(`/issues?labels=${LABEL}&state=open&per_page=1`))[0] || null;

const createIssue = async (client) => {
  if (!await client.json(`/labels/${LABEL}`)) {
    await client.json('/labels', {
      method: 'POST', body: { name: LABEL, color: 'c5def5', description: 'Weekly CI inventory reports' }
    });
  }
  return client.json('/issues', { method: 'POST', body: { title: ISSUE_TITLE, body: ISSUE_BODY, labels: [LABEL] } });
};

const main = async (argv) => {
  const post = argv.includes('--post');
  const repository = process.env.GITHUB_REPOSITORY;
  const token = process.env.GITHUB_TOKEN || process.env.GH_TOKEN;
  if (!repository || !token) throw new Error('set GITHUB_REPOSITORY and GITHUB_TOKEN (or GH_TOKEN)');
  const client = createClient({ token, repository, apiUrl: process.env.GITHUB_API_URL || 'https://api.github.com' });
  const until = new Date();
  const since = new Date(until.getTime() - 7 * DAY_MS);
  const week = isoWeek(until);

  const measurements = readMeasurements();
  const { testDir, casesPerFile } = listFunctionalCases();
  const { run: functionalRun, shards } = await latestFunctionalRun(client, measurements.shardCount);
  const functional = summarizeFunctional({
    reports: await readShardReports(client, shards), measurements, casesPerFile, testDir
  });
  // Before the job lists, which use the request budget until RATE_RESERVE is left.
  const mergedPulls = await listMergedPulls(client, since);
  const pullRequests = summarizePullRequests(await listPullRequestRuns(client, mergedPulls));
  const windowPulls = mergedPulls.filter((pull) => Date.parse(pull.merged_at) >= since.getTime()
    && Date.parse(pull.merged_at) < until.getTime());
  const earlierRuns = await listEarlierRuns(client, windowPulls, since);
  const { runs, partial } = await listRuns(client, since, until);
  const jobsByRun = await listJobs(client, runs);
  const summary = summarizeRuns({ runs, jobsByRun });
  const pullRequestCi = summarizePullRequestCi({ runs, earlierRuns, jobsByRun, pulls: windowPulls });

  const issue = await findIssue(client);
  const { thisWeek, previous } = priorReports(
    issue ? await client.pages(`/issues/${issue.number}/comments`, null) : [], week
  );

  const markdown = renderReport({
    week, since, until, runs: summary, partial, pullRequests, pullRequestCi, functional, functionalRun,
    testJobs: countWorkflowJobs(readFileSync(TESTS_WORKFLOW, 'utf8')),
    previous, serverUrl: process.env.GITHUB_SERVER_URL || 'https://github.com', repository
  });
  process.stdout.write(markdown);
  console.error(`GitHub API requests: ${client.requests()}`);
  if (!post) return;
  const target = issue || await createIssue(client);
  // A second run in the same week replaces that week's comment.
  if (thisWeek) await client.json(`/issues/comments/${thisWeek.id}`, { method: 'PATCH', body: { body: markdown } });
  else await client.json(`/issues/${target.number}/comments`, { method: 'POST', body: { body: markdown } });
};

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  main(process.argv.slice(2)).catch((error) => {
    console.error(error.message);
    process.exitCode = 1;
  });
}
