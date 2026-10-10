#!/usr/bin/env node
// The promotion audit table of a release-tier `Tests` run
// (.github/workflows/test.yml, `Promotion audit / summary`; docs/internal/WEB_PERIODIC_AUDIT.md).
// It lists each `Promotion audit / <key>` job with its result, the link to its
// `promotion-audit-<key>` artifact, and the findings counted in that artifact:
// failed Playwright tests (`report.json`), failed recipes (`recipes-status.tsv`),
// and DIFF rows of parity replays (`replay-report.json`). Recipes skipped for want
// of a LOSAT runtime are listed to run locally, not counted.
//
//   node tools/audit/release-audit-summary.mjs <downloaded artifacts dir> >> "$GITHUB_STEP_SUMMARY"
//
// Reads GITHUB_TOKEN, GITHUB_REPOSITORY, GITHUB_RUN_ID, GITHUB_RUN_ATTEMPT,
// GITHUB_API_URL, and GITHUB_SERVER_URL.
import { existsSync, readdirSync, readFileSync, statSync } from 'node:fs';
import { basename, join } from 'node:path';
import { pathToFileURL } from 'node:url';

const JOB_PREFIX = 'Promotion audit / ';
const ARTIFACT_PREFIX = 'promotion-audit-';

const files = (dir) => readdirSync(dir).flatMap((name) => {
  const path = join(dir, name);
  return statSync(path).isDirectory() ? files(path) : [path];
});

const readJson = (path) => {
  try {
    return JSON.parse(readFileSync(path, 'utf8'));
  } catch {
    return null;
  }
};

// The findings and notes of one artifact folder, and the journeys' coverage report.
export const readArtifact = (dir) => {
  const result = { findings: 0, notes: [], coverage: null };
  if (!dir || !existsSync(dir)) return { ...result, findings: null, notes: ['no artifact'] };
  for (const path of files(dir)) {
    const name = basename(path);
    if (name === 'report.json') {
      const unexpected = Number(readJson(path)?.stats?.unexpected);
      if (Number.isFinite(unexpected)) {
        result.findings += unexpected;
        if (unexpected) result.notes.push(`${unexpected} failed test(s)`);
      } else {
        result.notes.push('unreadable Playwright report');
      }
    } else if (name === 'recipes-status.tsv') {
      const rows = readFileSync(path, 'utf8').split('\n').filter(Boolean).map((line) => line.split('\t'));
      const skipped = rows.filter(([, , status]) => status === 'skipped-run-locally');
      const failed = rows.filter(([, , status]) => status !== 'pass' && status !== 'skipped-run-locally');
      result.findings += failed.length;
      result.notes.push(`${rows.length - failed.length - skipped.length} of ${rows.length - skipped.length} recipes pass`
        + (failed.length ? `; failed: ${failed.map(([, id]) => id).join(', ')}` : '')
        + (skipped.length ? `; run locally (no LOSAT runtime): ${skipped.map(([, id]) => id).join(', ')}` : ''));
    } else if (name === 'replay-report.json') {
      const rows = readJson(path) || [];
      const diff = rows.filter((row) => row?.outcome === 'DIFF').length;
      const noEffect = rows.filter((row) => row?.outcome === 'NO_EFFECT').length;
      result.findings += diff;
      result.notes.push(`${basename(join(path, '..'))}: ${diff} DIFF, ${noEffect} NO_EFFECT of ${rows.length}`);
    } else if (name === 'losat.txt') {
      const text = readFileSync(path, 'utf8').trim();
      if (!text.startsWith('available')) result.notes.push(`LOSAT runtime unavailable on the runner (${text})`);
    } else if (name === 'seed.txt') {
      result.notes.push(`seed ${readFileSync(path, 'utf8').trim()}`);
    } else if (name === 'coverage.json') {
      result.coverage = readJson(path);
    }
  }
  return result;
};

const cell = (text) => String(text ?? '').replace(/\|/g, '\\|').replace(/\n/g, ' ');

// Markdown: one row per audit job, then the capabilities no journey covers.
export const renderSummary = ({ jobs, artifacts, artifactsDir, runUrl }) => {
  const byName = new Map(artifacts.map((artifact) => [artifact.name, artifact]));
  const rows = [];
  let coverage = null;
  for (const job of jobs.filter((item) => item.name.startsWith(JOB_PREFIX)).sort((left, right) => left.name.localeCompare(right.name))) {
    const key = job.name.slice(JOB_PREFIX.length);
    if (key === 'summary') continue;
    const artifactName = `${ARTIFACT_PREFIX}${key}`;
    const artifact = byName.get(artifactName);
    const read = readArtifact(artifactsDir ? join(artifactsDir, artifactName) : null);
    coverage = coverage || read.coverage;
    rows.push(`| ${cell(key)} | ${cell(job.conclusion || job.status)} | `
      + `${artifact ? `[${cell(artifactName)}](${runUrl}/artifacts/${artifact.id})` : 'none'} | `
      + `${read.findings ?? 'unknown'} | ${cell(read.notes.join('; '))} |`);
  }
  const lines = [
    '## Promotion audit',
    '',
    'Sweeps are report-only: a sweep finding is an audit row, not a failed job (`tools/audit/README.md`).',
    '',
    '| Job | Result | Artifact | Findings | Notes |',
    '| --- | --- | --- | --- | --- |',
    ...(rows.length ? rows : ['| none | | | | |']),
    ''
  ];
  if (coverage?.error) {
    lines.push(`Uncovered capabilities: unknown (${coverage.error})`);
  } else if (coverage) {
    lines.push(`Changed capabilities since \`origin/main\`: ${coverage.changed.join(', ') || 'none'}`);
    lines.push(`Uncovered by the user journeys (hand look with a time limit, or an Owner waiver): ${coverage.uncovered.join(', ') || 'none'}`);
    const uncoveredPaths = coverage.uncoveredPaths || {};
    if (coverage.uncovered.length) {
      lines.push('', '| Uncovered capability | Changed paths |', '| --- | --- |');
      for (const capability of coverage.uncovered) {
        lines.push(`| ${cell(capability)} | ${cell((uncoveredPaths[capability] || []).map((path) => `\`${path}\``).join(', ') || 'unknown')} |`);
      }
    }
  }
  return `${lines.join('\n')}\n`;
};

const listAll = async ({ url, field, token }) => {
  const items = [];
  for (let page = 1; ; page += 1) {
    const response = await fetch(`${url}${url.includes('?') ? '&' : '?'}per_page=100&page=${page}`, {
      headers: { Authorization: `Bearer ${token}`, Accept: 'application/vnd.github+json' }
    });
    if (!response.ok) throw new Error(`${url}: HTTP ${response.status}`);
    const batch = (await response.json())[field] || [];
    items.push(...batch);
    if (batch.length < 100) return items;
  }
};

const main = async () => {
  const artifactsDir = process.argv[2];
  const env = process.env;
  const api = `${env.GITHUB_API_URL || 'https://api.github.com'}/repos/${env.GITHUB_REPOSITORY}/actions/runs/${env.GITHUB_RUN_ID}`;
  const runUrl = `${env.GITHUB_SERVER_URL || 'https://github.com'}/${env.GITHUB_REPOSITORY}/actions/runs/${env.GITHUB_RUN_ID}`;
  const jobs = await listAll({ url: `${api}/attempts/${env.GITHUB_RUN_ATTEMPT || 1}/jobs`, field: 'jobs', token: env.GITHUB_TOKEN });
  const artifacts = await listAll({ url: `${api}/artifacts`, field: 'artifacts', token: env.GITHUB_TOKEN });
  process.stdout.write(renderSummary({ jobs, artifacts, artifactsDir, runUrl }));
};

if (import.meta.url === pathToFileURL(process.argv[1] || '').href) {
  main().catch((error) => {
    console.error(error?.stack || error);
    process.exitCode = 1;
  });
}
