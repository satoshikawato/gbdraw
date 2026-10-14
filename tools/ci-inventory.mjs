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
// come from the data block of the workflow's own previous comment. Job lists are read
// only for `Tests` and `Gallery publication`; other runs count by run duration.
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
  SHARDS_FILE, WORKERS_PER_SHARD, casesPerFileOf, listFunctionalCases, measureMinutes, planShards
} from './balance-functional-shards.mjs';

const LABEL = 'ci-inventory';
// Only comments of the workflow's own token are reports; anyone can comment on the issue.
const REPORT_AUTHOR = 'github-actions[bot]';
const TESTS_WORKFLOW = '.github/workflows/test.yml';
// Workflows whose job lists are read. The others count by run duration, which keeps a
// normal week within the 1,000 requests an hour of a GITHUB_TOKEN.
const JOB_LIST_WORKFLOWS = new Set([TESTS_WORKFLOW, '.github/workflows/gallery-publication.yml']);
const RUN_DURATION = ' (run duration)';
const DAY_MS = 24 * 60 * 60 * 1000;
const REBALANCE_RATIO = 1.2;
const TOP = 15;
// Requests kept for the artifact, issue, and comment calls once the job lists stop.
const RATE_RESERVE = 100;
const PARALLEL_REQUESTS = 8;
const DATA_BLOCK = /<!-- ci-inventory-data (\{.*?\}) -->/s;
const WEEK = /^\d{4}-W\d{2}$/;
const DATA_NUMBERS = ['cases', 'functionalTestMinutes', 'runnerMinutes', 'testJobs'];

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
// shard numbers (a shard without a report is listed as missing); `shardsFile` is
// tests/ci/functional-shards.json; `casesPerFile` and `testDir` come from listing the
// functional suite.
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
  const casesRun = casesPerFileOf(reports.map(({ report }) => report), testDir);
  const plan = planShards({ casesPerFile, fileMinutes, shardCount: shardsFile.shards.length });
  const testMinutes = sum(shards.map((shard) => shard.testMinutes));
  const mean = testMinutes / (shards.length || 1);
  const slowest = Math.max(0, ...shards.map((shard) => shard.testMinutes));
  const unmeasured = [...casesPerFile.keys()].filter((file) => !(file in shardsFile.minutes)).sort();
  // One rounded ratio decides the flag and is the number the report prints.
  const slowestToMean = mean ? Math.round((slowest / mean) * 100) / 100 : 0;
  return {
    shards,
    missingShards: shardsFile.shards.map((_, index) => index + 1)
      .filter((number) => !shards.some((shard) => shard.shard === number)),
    cases: sum([...casesRun.values()]),
    testMinutes: round1(testMinutes),
    slowestToMean,
    unmeasured,
    rebalance: slowestToMean > REBALANCE_RATIO || unmeasured.length > 0,
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

// The numbers of a report's data block, or null unless every field has its type.
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
  week, since, until, runs, partial = null, functional, functionalRun, testJobs, previous, serverUrl, repository
}) => {
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
    ...(partial ? ['', `**Partial report**: ${partial}`] : []),
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
    `### Functional shards of run ${runLink} (\`${functionalRun.head_sha.slice(0, 8)}\`)`,
    '',
    `Test-minutes add up case durations; a shard runs ${WORKERS_PER_SHARD} workers. `
      + (functional.planChanges
        ? 'The planned column is what `tools/balance-functional-shards.mjs` would assign from this run.'
        : '`tools/balance-functional-shards.mjs` would keep the current plan.')
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

const createClient = ({ token, apiUrl, repository }) => {
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
  const json = async (path, options) => {
    const response = await request(path, options);
    if (response.status === 404 && (!options || options.method === undefined)) return null;
    if (!response.ok) {
      throw Object.assign(
        new Error(`${options?.method || 'GET'} ${path}: ${response.status} ${await response.text()}`),
        { status: response.status }
      );
    }
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
  return { request, json, pages, remaining: () => remaining, requests: () => requests };
};

// A rate-limit or permission refusal ends a listing early; the report says so.
const refused = (error) => [403, 429].includes(error.status);

// The completed runs created in the window, and why the listing stopped early, if it did.
// A filtered run listing returns at most 1,000 runs, so the window is read a day at a time.
const listRuns = async (client, since, until) => {
  const runs = [];
  let partial = null;
  for (let start = since.getTime(); start < until.getTime() && !partial; start += DAY_MS) {
    const end = Math.min(start + DAY_MS, until.getTime());
    const created = `${new Date(start).toISOString()}..${new Date(end - 1000).toISOString()}`;
    try {
      runs.push(...await client.pages(`/actions/runs?status=completed&created=${encodeURIComponent(created)}`, 'workflow_runs'));
    } catch (error) {
      if (!refused(error)) throw error;
      partial = `runs created from ${new Date(start).toISOString().slice(0, 16)}Z on are missing (HTTP ${error.status}).`;
    }
  }
  return { runs: [...new Map(runs.map((run) => [run.id, run])).values()], partial };
};

// Job lists of every attempt of `Tests` and `Gallery publication` runs, `Tests` first,
// while the request budget lasts.
const listJobs = async (client, runs) => {
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
      const redirect = await client.request(artifact.archive_download_url);
      const location = redirect.headers.get('location');
      if (!location) throw new Error(`artifact ${artifact.name}: ${redirect.status}`);
      const archive = join(dir, `${artifact.id}.zip`);
      writeFileSync(archive, Buffer.from(await (await fetch(location)).arrayBuffer()));
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

  const shardsFile = JSON.parse(readFileSync(SHARDS_FILE, 'utf8'));
  const { testDir, casesPerFile } = listFunctionalCases();
  const { run: functionalRun, shards } = await latestFunctionalRun(client, shardsFile.shards.length);
  const functional = summarizeFunctional({
    reports: await readShardReports(client, shards), shardsFile, casesPerFile, testDir
  });
  const { runs, partial } = await listRuns(client, since, until);
  const summary = summarizeRuns({ runs, jobsByRun: await listJobs(client, runs) });

  const issue = await findIssue(client);
  const { thisWeek, previous } = priorReports(
    issue ? await client.pages(`/issues/${issue.number}/comments`, null) : [], week
  );

  const markdown = renderReport({
    week, since, until, runs: summary, partial, functional, functionalRun,
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
