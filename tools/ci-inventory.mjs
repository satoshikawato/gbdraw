#!/usr/bin/env node
// Weekly CI inventory (docs/internal/SELECTIVE_CI.md "Weekly CI inventory").
//
//   node tools/ci-inventory.mjs [--post]
//
// Reports the last 7 days of GitHub Actions use: runner-minutes by workflow and event
// and by job, the wall and queue time of `Tests` runs, and the functional Playwright
// shards of the latest green dev `Tests` run with the plan that
// tools/balance-functional-shards.mjs would write. It prints Markdown; with --post it
// also adds the report as this week's comment on the open issue labelled
// `ci-inventory` (creating the label and issue when missing). Week-over-week changes
// come from the data block of the previous comment.
//
// Reads GITHUB_TOKEN (or GH_TOKEN), GITHUB_REPOSITORY, GITHUB_API_URL, and
// GITHUB_SERVER_URL. Run it from the repository root: it lists the functional cases
// with @playwright/test and reads tests/ci/functional-shards.json and
// .github/workflows/test.yml.
import { spawnSync } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

import {
  SHARDS_FILE, WORKERS_PER_SHARD, listFunctionalCases, measureMinutes, planShards, specsOf
} from './balance-functional-shards.mjs';

export const LABEL = 'ci-inventory';
export const TESTS_WORKFLOW = '.github/workflows/test.yml';
const DAY_MS = 24 * 60 * 60 * 1000;
const REBALANCE_RATIO = 1.2;
const TOP = 15;
// Requests kept for the artifact, issue, and comment calls once the job lists stop.
const RATE_RESERVE = 100;
const DATA_BLOCK = /<!-- ci-inventory-data (\{.*?\}) -->/s;

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
// of every attempt; a run without a job list counts its own duration as one job named
// after the workflow.
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
      add(byWorkflow, workflow, minutes, 1);
      add(byJob, workflowLabel(run), minutes);
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

// The functional shards of one green run. `reports` are the shard reports with their
// shard numbers; `shardsFile` is tests/ci/functional-shards.json; `casesPerFile` and
// `testDir` come from listing the functional suite.
export const summarizeFunctional = ({ reports, shardsFile, casesPerFile, testDir }) => {
  const shardTestMinutes = (report) => sum(measureMinutes([report], testDir).values());
  const shards = reports
    .map(({ shard, report }) => ({
      shard,
      testMinutes: round1(shardTestMinutes(report)),
      wallMinutes: round1((report.stats?.duration || 0) / 60000)
    }))
    .sort((a, b) => a.shard - b.shard);
  const fileMinutes = measureMinutes(reports.map(({ report }) => report), testDir);
  const casesRun = new Map();
  for (const { report } of reports) {
    for (const spec of specsOf(report, testDir)) {
      casesRun.set(spec.file, (casesRun.get(spec.file) || 0) + spec.tests.length);
    }
  }
  const plan = planShards({ casesPerFile, fileMinutes, shardCount: shardsFile.shards.length });
  const testMinutes = sum(shards.map((shard) => shard.testMinutes));
  const mean = testMinutes / (shards.length || 1);
  const slowest = Math.max(0, ...shards.map((shard) => shard.testMinutes));
  const unmeasured = [...casesPerFile.keys()].filter((file) => !(file in shardsFile.minutes)).sort();
  return {
    shards,
    cases: sum([...casesRun.values()]),
    testMinutes: round1(testMinutes),
    slowestToMean: mean ? Math.round((slowest / mean) * 100) / 100 : 0,
    unmeasured,
    rebalance: slowest > REBALANCE_RATIO * mean || unmeasured.length > 0,
    planChanges: JSON.stringify(plan.shards) !== JSON.stringify(shardsFile.shards),
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

export const readData = (body) => {
  const match = body?.match(DATA_BLOCK);
  if (!match) return null;
  try {
    return JSON.parse(match[1]);
  } catch {
    return null;
  }
};

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
export const renderReport = ({ week, since, until, runs, functional, functionalRun, testJobs, previous, serverUrl, repository }) => {
  const data = {
    week,
    cases: functional.cases,
    functionalTestMinutes: functional.testMinutes,
    runnerMinutes: runs.runnerMinutes,
    testJobs
  };
  const reasons = [
    functional.slowestToMean > REBALANCE_RATIO
      ? `slowest shard ${functional.slowestToMean}x the mean (limit ${REBALANCE_RATIO}x)` : null,
    functional.unmeasured.length
      ? `${functional.unmeasured.length} spec file(s) without recorded minutes: ${functional.unmeasured.map((file) => `\`${file}\``).join(', ')}`
      : null
  ].filter(Boolean);
  const runLink = `[${functionalRun.id}](${serverUrl}/${repository}/actions/runs/${functionalRun.id})`;
  const perWorker = (minutes) => fmt(minutes / WORKERS_PER_SHARD);
  return [
    `## CI inventory ${week}`,
    '',
    `Window: ${since.toISOString().slice(0, 16)}Z to ${until.toISOString().slice(0, 16)}Z.`
      + ` What to do with a flagged report: docs/internal/SELECTIVE_CI.md, "Weekly CI inventory".`,
    '',
    functional.rebalance
      ? `**Rebalance the functional shards**: ${reasons.join('; ')}.`
      : `Functional shards are balanced (slowest ${functional.slowestToMean}x the mean).`,
    '',
    table(['Measure', 'This week', 'Previous', 'Change'], [
      ['Functional cases', data.cases, previous?.cases ?? '-', delta(data.cases, previous?.cases)],
      ['Functional test-minutes', fmt(data.functionalTestMinutes), fmt(previous?.functionalTestMinutes),
        delta(data.functionalTestMinutes, previous?.functionalTestMinutes)],
      ['Runner-minutes, all workflows', fmt(data.runnerMinutes), fmt(previous?.runnerMinutes),
        delta(data.runnerMinutes, previous?.runnerMinutes)],
      ['Jobs in `test.yml`', data.testJobs, previous?.testJobs ?? '-', delta(data.testJobs, previous?.testJobs)]
    ]),
    previous ? `\nPrevious: ${previous.week}.` : '\nNo previous report.',
    '',
    '### Runner-minutes by workflow and event',
    '',
    table(['Workflow [event]', 'Runs', 'Runner-minutes'],
      runs.byWorkflow.map((entry) => [entry.name, entry.runs, fmt(entry.minutes)])),
    '',
    `${runs.runs} completed runs.`
      + (runs.runDurationOnly ? ` ${runs.runDurationOnly} run(s) counted by run duration because the API request budget ran out.` : ''),
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
    `### Functional shards of run ${runLink} (\`${functionalRun.head_sha.slice(0, 8)}\`)`,
    '',
    `Test-minutes add up case durations; a shard runs ${WORKERS_PER_SHARD} workers. `
      + (functional.planChanges
        ? 'The planned column is what `tools/balance-functional-shards.mjs` would assign from this run.'
        : '`tools/balance-functional-shards.mjs` would keep the current plan.'),
    '',
    table(['Shard', 'Test-minutes', 'Per worker', 'Reporter wall', 'Planned per worker'],
      functional.shards.map((shard, index) => [
        shard.shard, fmt(shard.testMinutes), perWorker(shard.testMinutes), fmt(shard.wallMinutes),
        perWorker(functional.plannedTestMinutes[index] ?? 0)
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

const createClient = ({ token, apiUrl, repository }) => {
  let remaining = Infinity;
  const request = async (path, { method = 'GET', body } = {}) => {
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
  const json = async (path, options) => {
    const response = await request(path, options);
    if (response.status === 404 && (!options || options.method === undefined)) return null;
    if (!response.ok) throw new Error(`${options?.method || 'GET'} ${path}: ${response.status} ${await response.text()}`);
    return response.json();
  };
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
  return { request, json, pages, remaining: () => remaining };
};

// The completed runs created in the window. A filtered run listing returns at most 1,000
// runs, so the window is read one day at a time.
const listRuns = async (client, since, until) => {
  const runs = [];
  for (let start = since.getTime(); start < until.getTime(); start += DAY_MS) {
    const end = Math.min(start + DAY_MS, until.getTime());
    const created = `${new Date(start).toISOString()}..${new Date(end - 1000).toISOString()}`;
    runs.push(...await client.pages(`/actions/runs?status=completed&created=${encodeURIComponent(created)}`, 'workflow_runs'));
  }
  return [...new Map(runs.map((run) => [run.id, run])).values()];
};

// Job lists of every attempt, `Tests` runs first, until the request budget is spent.
const listJobs = async (client, runs) => {
  const jobsByRun = new Map();
  const ordered = [...runs].sort((a, b) => Number(b.path === TESTS_WORKFLOW) - Number(a.path === TESTS_WORKFLOW));
  for (const run of ordered) {
    if (client.remaining() <= RATE_RESERVE) break;
    jobsByRun.set(run.id, await client.pages(`/actions/runs/${run.id}/jobs?filter=all`, 'jobs'));
  }
  return jobsByRun;
};

const SHARD_ARTIFACT = /^playwright-functional-shard-(\d+)-/;

// The latest green dev push `Tests` run that kept a functional report for every shard.
const latestFunctionalRun = async (client, shardCount) => {
  const { workflow_runs: runs } = await client.json(
    `/actions/workflows/test.yml/runs?branch=dev&event=push&status=success&per_page=20`
  );
  for (const run of runs) {
    const artifacts = (await client.pages(`/actions/runs/${run.id}/artifacts`, 'artifacts'))
      .filter((artifact) => SHARD_ARTIFACT.test(artifact.name) && !artifact.expired);
    if (new Set(artifacts.map((artifact) => artifact.name.match(SHARD_ARTIFACT)[1])).size === shardCount) {
      return { run, artifacts };
    }
  }
  throw new Error('no green dev Tests run in the last 20 kept every functional shard report');
};

// Downloads each shard artifact and reads its functional-report.json.
const readShardReports = async (client, artifacts) => {
  const dir = mkdtempSync(join(process.env.RUNNER_TEMP || tmpdir(), 'ci-inventory-'));
  try {
    const reports = [];
    for (const artifact of artifacts) {
      const redirect = await client.request(artifact.archive_download_url);
      const location = redirect.headers.get('location');
      if (!location) throw new Error(`artifact ${artifact.name}: ${redirect.status}`);
      const archive = join(dir, `${artifact.id}.zip`);
      writeFileSync(archive, Buffer.from(await (await fetch(location)).arrayBuffer()));
      const unzip = spawnSync('unzip', ['-p', archive, 'functional-report.json'], {
        encoding: 'utf8', maxBuffer: 256 * 1024 * 1024
      });
      if (unzip.status !== 0) throw new Error(`artifact ${artifact.name}: ${unzip.stderr}`);
      reports.push({ shard: Number(artifact.name.match(SHARD_ARTIFACT)[1]), report: JSON.parse(unzip.stdout) });
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

  const shardsFile = JSON.parse(readFileSync(SHARDS_FILE, 'utf8'));
  const { testDir, casesPerFile } = listFunctionalCases();
  const { run: functionalRun, artifacts } = await latestFunctionalRun(client, shardsFile.shards.length);
  const functional = summarizeFunctional({
    reports: await readShardReports(client, artifacts), shardsFile, casesPerFile, testDir
  });
  const runs = await listRuns(client, since, until);
  const summary = summarizeRuns({ runs, jobsByRun: await listJobs(client, runs) });

  const issue = await findIssue(client);
  const comments = issue ? (await client.pages(`/issues/${issue.number}/comments`, null))
    .filter((comment) => readData(comment.body)) : [];
  const thisWeek = comments.at(-1) && readData(comments.at(-1).body).week === week ? comments.at(-1) : null;
  const previous = readData((thisWeek ? comments.at(-2) : comments.at(-1))?.body);

  const markdown = renderReport({
    week, since, until, runs: summary, functional, functionalRun,
    testJobs: countWorkflowJobs(readFileSync(TESTS_WORKFLOW, 'utf8')),
    previous, serverUrl: process.env.GITHUB_SERVER_URL || 'https://github.com', repository
  });
  process.stdout.write(markdown);
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
