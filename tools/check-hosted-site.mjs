#!/usr/bin/env node
// Hosted site check (docs/internal/WEB_PERIODIC_AUDIT.md, "Hosted site check and rollback").
//
//   node tools/check-hosted-site.mjs --sha <40-hex commit> [--url https://gbdraw.app/]
//     [--build-timeout <s>] [--serve-timeout <s>] [--startup-timeout <s>] [--interval <s>]
//
// Confirms that a commit pushed to main reached the hosted site:
//   1. waits for the `Workers Builds: gbdraw` check run that Cloudflare posts on the commit
//      and requires it to succeed (its summary names the deployed Version ID);
//   2. polls the site until the build stamp of its index.html (tools/stamp_web_build.py)
//      names the commit; a stamp of a later main commit counts as superseded, not failed;
//   3. loads the site once in Chromium and waits for the app's startup predicate (the one
//      tools/verify_gui_offline.py uses) on a cross-origin isolated page without page errors.
// It reports; it gates nothing. The summary goes to stdout and, when set, $GITHUB_STEP_SUMMARY.
//
// Reads GITHUB_TOKEN (or GH_TOKEN), GITHUB_REPOSITORY, and GITHUB_API_URL. Step 3 needs
// @playwright/test and its Chromium (npm ci; npx playwright install chromium).
import { appendFileSync } from 'node:fs';
import { setTimeout as sleep } from 'node:timers/promises';
import { pathToFileURL } from 'node:url';

export const CHECK_NAME = 'Workers Builds: gbdraw';
export const CHECK_APP = 'cloudflare-workers-and-pages';
const DEFAULTS = {
  url: 'https://gbdraw.app/',
  buildTimeout: 1800,
  serveTimeout: 900,
  startupTimeout: 120,
  interval: 30
};
const SHA = /^[0-9a-f]{40}$/;

// The commit named by the build stamp, or null for a page without one.
export const parseBuildStamp = (html) => (
  html.match(/<span class="hosted-build-label" title="Commit ([0-9a-f]{40})">/)?.[1] ?? null
);

// The Version ID line of a Workers Builds check run summary, or null.
export const parseVersionId = (summary) => (
  (summary || '').match(/^Version ID: ([0-9a-f-]{36})\s*$/m)?.[1] ?? null
);

// The newest Workers Builds check run among a commit's check runs, or null.
export const latestBuildRun = (checkRuns) => checkRuns
  .filter((run) => run.name === CHECK_NAME && run.app?.slug === CHECK_APP)
  .reduce((newest, run) => (newest && newest.id > run.id ? newest : run), null);

// What a served stamp means for the expected commit. `compareStatus` is the status of
// GitHub's compare of expected...served: "ahead" means the served commit descends from it.
export const servedState = ({ expected, served, compareStatus }) => {
  if (served === expected) return 'match';
  if (served && compareStatus === 'ahead') return 'superseded';
  return 'pending';
};

// A non-OK answer of the site. A client error (such as a Cloudflare challenge, which
// sets cf-mitigated) does not go away by waiting; a server error may, during a deploy.
export const responseProblem = (status, mitigated) => ({
  retry: status >= 500,
  text: `HTTP ${status}${mitigated ? ` (cf-mitigated: ${mitigated})` : ''}`
});

export const parseArgs = (argv) => {
  const options = { ...DEFAULTS };
  const names = {
    '--sha': 'sha', '--url': 'url', '--build-timeout': 'buildTimeout',
    '--serve-timeout': 'serveTimeout', '--startup-timeout': 'startupTimeout', '--interval': 'interval'
  };
  for (let index = 0; index < argv.length; index += 2) {
    const key = names[argv[index]];
    const value = argv[index + 1];
    if (!key || value === undefined) throw new Error(`unknown or incomplete argument: ${argv[index]}`);
    options[key] = key === 'sha' || key === 'url' ? value : Number(value);
    if (typeof options[key] === 'number' && !(options[key] >= 0)) {
      throw new Error(`${argv[index]} takes a number of seconds`);
    }
  }
  if (!SHA.test(options.sha || '')) throw new Error('--sha takes a full 40-character commit SHA');
  return options;
};

export const renderSummary = (report) => [
  `### Hosted site check: ${report.failure ? 'failed' : 'passed'}`,
  '',
  `- Commit: \`${report.sha}\``,
  ...(report.versionId ? [`- ${CHECK_NAME}: success, Version ID \`${report.versionId}\``] : []),
  ...(report.served === 'match' ? [`- ${report.url} serves this commit`] : []),
  ...(report.served === 'superseded'
    ? [`- ${report.url} serves \`${report.servedSha}\`, a later main commit (superseded)`] : []),
  ...(report.startupSeconds !== undefined
    ? [`- App started in ${report.startupSeconds.toFixed(1)} s, cross-origin isolated`] : []),
  ...(report.failure ? [`- Failure: ${report.failure}`] : []),
  ''
].join('\n');

const createGitHub = () => {
  const repository = process.env.GITHUB_REPOSITORY;
  const token = process.env.GITHUB_TOKEN || process.env.GH_TOKEN;
  if (!repository || !token) throw new Error('set GITHUB_REPOSITORY and GITHUB_TOKEN (or GH_TOKEN)');
  const apiUrl = process.env.GITHUB_API_URL || 'https://api.github.com';
  return async (path) => {
    const response = await fetch(`${apiUrl}/repos/${repository}${path}`, {
      headers: {
        accept: 'application/vnd.github+json',
        authorization: `Bearer ${token}`,
        'x-github-api-version': '2022-11-28'
      }
    });
    if (!response.ok) throw new Error(`GitHub ${path}: HTTP ${response.status}`);
    return response.json();
  };
};

// Calls `probe` every `interval` seconds until it returns a value other than undefined
// or `timeout` seconds pass; then returns undefined.
const poll = async (probe, { timeout, interval }) => {
  const deadline = Date.now() + timeout * 1000;
  for (;;) {
    const value = await probe();
    if (value !== undefined || Date.now() + interval * 1000 > deadline) return value;
    await sleep(interval * 1000);
  }
};

const waitForBuild = async (github, options, report) => {
  const query = `?check_name=${encodeURIComponent(CHECK_NAME)}&per_page=100`;
  const run = await poll(async () => {
    const found = latestBuildRun((await github(`/commits/${options.sha}/check-runs${query}`)).check_runs);
    return found?.status === 'completed' ? found : undefined;
  }, { timeout: options.buildTimeout, interval: options.interval });
  if (!run) throw new Error(`no completed ${CHECK_NAME} check run on the commit after ${options.buildTimeout} s`);
  if (run.conclusion !== 'success') {
    throw new Error(`${CHECK_NAME} concluded ${run.conclusion} (${run.details_url || run.html_url})`);
  }
  report.versionId = parseVersionId(run.output?.summary) || 'not reported';
};

const waitForServedCommit = async (github, options, report) => {
  let lastServed = null;
  let lastSeen = 'nothing';
  const state = await poll(async () => {
    const url = new URL(options.url);
    url.searchParams.set('hosted-check', String(Date.now()));
    let response;
    try {
      response = await fetch(url, { cache: 'no-store' });
    } catch (error) {
      lastSeen = `no response (${error.cause?.code || error.message})`;
      return undefined;
    }
    if (!response.ok) {
      const problem = responseProblem(response.status, response.headers.get('cf-mitigated'));
      if (!problem.retry) throw new Error(`${options.url} answered ${problem.text}`);
      lastSeen = problem.text;
      return undefined;
    }
    lastServed = parseBuildStamp(await response.text());
    lastSeen = lastServed || 'a page without a build stamp';
    if (!lastServed) return undefined;
    const compareStatus = lastServed === options.sha
      ? 'identical'
      : (await github(`/compare/${options.sha}...${lastServed}`)).status;
    const result = servedState({ expected: options.sha, served: lastServed, compareStatus });
    return result === 'pending' ? undefined : result;
  }, { timeout: options.serveTimeout, interval: options.interval });
  if (!state) {
    throw new Error(`${options.url} still serves ${lastSeen} after ${options.serveTimeout} s`);
  }
  report.served = state;
  report.servedSha = lastServed;
};

const checkStartup = async (options, report) => {
  const { chromium } = await import('@playwright/test');
  const browser = await chromium.launch();
  try {
    const page = await browser.newPage();
    const pageErrors = [];
    page.on('pageerror', (error) => pageErrors.push(error.message));
    const started = Date.now();
    await page.goto(options.url, { waitUntil: 'domcontentloaded' });
    try {
      await page.waitForFunction(() => {
        const app = window.__GBDRAW_APP__;
        return Boolean(app) && Object.keys(app.paletteDefinitions || {}).length > 0;
      }, null, { timeout: options.startupTimeout * 1000 });
    } catch (error) {
      throw new Error(`the app did not start within ${options.startupTimeout} s`
        + (pageErrors.length ? `; page errors: ${pageErrors.join(' | ')}` : ''), { cause: error });
    }
    const seconds = (Date.now() - started) / 1000;
    if (pageErrors.length) throw new Error(`page errors during startup: ${pageErrors.join(' | ')}`);
    if (await page.evaluate(() => window.crossOriginIsolated) !== true) {
      throw new Error('the page is not cross-origin isolated');
    }
    report.startupSeconds = seconds;
  } finally {
    await browser.close();
  }
};

const main = async (argv) => {
  const options = parseArgs(argv);
  const report = { sha: options.sha, url: options.url };
  try {
    const github = createGitHub();
    await waitForBuild(github, options, report);
    await waitForServedCommit(github, options, report);
    await checkStartup(options, report);
  } catch (error) {
    report.failure = error.message;
  }
  const summary = renderSummary(report);
  process.stdout.write(summary);
  if (process.env.GITHUB_STEP_SUMMARY) appendFileSync(process.env.GITHUB_STEP_SUMMARY, summary);
  if (report.failure) process.exitCode = 1;
};

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  main(process.argv.slice(2)).catch((error) => {
    console.error(error.message);
    process.exitCode = 1;
  });
}
