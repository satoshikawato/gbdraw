#!/usr/bin/env node
// Owner-graph report for gbdraw/web/js (implementation plan Phase B1).
//
//   node tools/report-web-owner-graph.mjs --at <revision|worktree> [--json]
//   node tools/report-web-owner-graph.mjs --range <from>..<to> [--first-parent] [--json]
//
// `--at` prints the summary and every subject of each detector at one
// revision. `--range` prints one row per commit in the range that touched
// gbdraw/web/js (first-parent merges only with --first-parent), which is the
// trend table of the audit. Report only: nothing here gates a change.
import { execFileSync } from 'node:child_process';
import { existsSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

import {
  detectWebOwnerGraph,
  summarizeWebOwnerGraph,
  WEB_OWNER_GRAPH_DEFAULTS,
  WEB_OWNER_GRAPH_DETECTOR_IDS
} from './web-owner-graph-detectors.mjs';

const WEB_SOURCE_PREFIX = 'gbdraw/web/js/';
const OWNER_GRAPH_REGISTRY_PATH = 'tools/web-owner-graph.json';

const runGit = (root, args, options = {}) => execFileSync('git', ['-C', root, ...args], {
  encoding: 'utf8',
  maxBuffer: 256 * 1024 * 1024,
  ...options
});

const repositoryRoot = runGit(process.cwd(), ['rev-parse', '--show-toplevel']).trim();

const listSources = (revision) => runGit(repositoryRoot, [
  'ls-tree', '-r', '--name-only', revision, '--', WEB_SOURCE_PREFIX.slice(0, -1)
]).split('\n').filter((path) => /\.[cm]?js$/.test(path));

const readRevisionSources = (revision) => {
  const paths = listSources(revision);
  if (!paths.length) return new Map();
  const input = `${paths.map((path) => `${revision}:${path}`).join('\n')}\n`;
  const output = execFileSync('git', ['-C', repositoryRoot, 'cat-file', '--batch'], { input, maxBuffer: 256 * 1024 * 1024 });
  const sources = new Map();
  let offset = 0;
  paths.forEach((path) => {
    const headerEnd = output.indexOf(10, offset);
    const header = output.subarray(offset, headerEnd).toString('utf8');
    const size = Number(header.split(' ').at(-1));
    if (!Number.isInteger(size)) throw new Error(`Cannot read ${path} at ${revision}: ${header}`);
    const contentStart = headerEnd + 1;
    sources.set(path, output.subarray(contentStart, contentStart + size).toString('utf8'));
    offset = contentStart + size + 1;
  });
  return sources;
};

const readWorktreeSources = () => {
  const tracked = runGit(repositoryRoot, [
    'ls-files', '--cached', '--modified', '--others', '--exclude-standard', '--', WEB_SOURCE_PREFIX.slice(0, -1)
  ]).split('\n').filter(Boolean);
  return new Map([...new Set(tracked)]
    .filter((path) => /\.[cm]?js$/.test(path) && existsSync(join(repositoryRoot, path)))
    .map((path) => [path, readFileSync(join(repositoryRoot, path), 'utf8')]));
};

const readRegistry = (revision) => {
  try {
    const source = revision === 'worktree'
      ? readFileSync(join(repositoryRoot, OWNER_GRAPH_REGISTRY_PATH), 'utf8')
      : runGit(repositoryRoot, ['show', `${revision}:${OWNER_GRAPH_REGISTRY_PATH}`], { stdio: ['ignore', 'pipe', 'ignore'] });
    return JSON.parse(source);
  } catch (_error) {
    return null;
  }
};

// Exported for tests: every detector result and the summary at one revision
// (`'worktree'` reads the working tree).
export const evaluateWebOwnerGraphAt = (revision) => {
  const sources = revision === 'worktree' ? readWorktreeSources() : readRevisionSources(revision);
  const registry = readRegistry(revision);
  const results = detectWebOwnerGraph(sources, registry);
  return { revision, registry: registry ? OWNER_GRAPH_REGISTRY_PATH : 'defaults', results, summary: summarizeWebOwnerGraph(results) };
};
export const revisionExists = (revision) => {
  try {
    runGit(repositoryRoot, ['cat-file', '-e', `${revision}^{commit}`], { stdio: ['ignore', 'pipe', 'ignore'] });
    return true;
  } catch (_error) {
    return false;
  }
};
const evaluate = evaluateWebOwnerGraphAt;

const summaryColumns = (summary) => ({
  injectionEdges: summary.injectionEdges,
  forwardClosures: summary.forwardClosures,
  stateBackdoors: summary.stateBackdoors,
  wholeObjectPorts: summary.wholeObjectPorts,
  triggerSites: summary.triggerSites,
  triggerModules: summary.triggerModules,
  ...Object.fromEntries(Object.entries(summary.projectionShapes).map(([name, count]) => [`shapes:${name}`, count]))
});

const printAt = (revision) => {
  const report = evaluate(revision);
  if (flags.has('--json')) {
    process.stdout.write(`${JSON.stringify(report, null, 2)}\n`);
    return;
  }
  const lines = [
    `# Web owner graph at ${report.revision}`,
    '',
    `- Registry: ${report.registry}`,
    `- Composition roots: ${(readRegistry(revision)?.compositionRoots || WEB_OWNER_GRAPH_DEFAULTS.compositionRoots).join(', ')}`,
    '',
    '## Summary',
    '',
    ...Object.entries(summaryColumns(report.summary)).map(([key, value]) => `- ${key}: ${value}`),
    ''
  ];
  WEB_OWNER_GRAPH_DETECTOR_IDS.forEach((id) => {
    const result = report.results[id];
    lines.push(`## ${id} (${result.subjects.length} subject(s))`, '');
    if (id === 'heavy-derived.trigger-site.v1') {
      Object.entries(result.countsBySubject).forEach(([subject, count]) => lines.push(`- ${subject}: ${count}`));
    } else if (id === 'projection.call-shape.v1') {
      result.observedShapes.forEach(({ domain, path, line, shape }) => lines.push(`- ${domain} ${path}:${line} \`${shape}\``));
    } else {
      result.subjects.forEach((subject) => lines.push(`- ${subject}`));
    }
    if (!result.subjects.length) lines.push('- None');
    lines.push('');
  });
  process.stdout.write(`${lines.join('\n')}\n`);
};

const printRange = (range) => {
  const [from, to] = range.split('..');
  if (!from || !to) throw new Error('--range expects <from>..<to>');
  const logArguments = ['log', '--reverse', '--format=%h%x09%ad%x09%s', '--date=short'];
  if (flags.has('--first-parent')) logArguments.push('--first-parent');
  logArguments.push(`${from}..${to}`, '--', WEB_SOURCE_PREFIX.slice(0, -1));
  const commits = runGit(repositoryRoot, logArguments).trim().split('\n').filter(Boolean)
    .map((line) => { const [sha, date, subject] = line.split('\t'); return { sha, date, subject }; });
  const rows = [{ sha: from, date: '', subject: 'baseline', columns: summaryColumns(evaluate(from).summary) }];
  commits.forEach((commit) => rows.push({ ...commit, columns: summaryColumns(evaluate(commit.sha).summary) }));
  if (flags.has('--json')) {
    process.stdout.write(`${JSON.stringify(rows, null, 2)}\n`);
    return;
  }
  const keys = Object.keys(rows[0].columns);
  const header = ['commit', 'date', ...keys, 'subject'].join('\t');
  const body = rows.map((row, index) => {
    const previous = rows[index - 1];
    const cells = keys.map((key) => {
      const value = row.columns[key];
      return previous && previous.columns[key] !== value ? `${value}*` : String(value);
    });
    const label = /Merge pull request (#\d+)/.exec(row.subject)?.[1] || row.subject.slice(0, 40);
    return [row.sha, row.date, ...cells, label].join('\t');
  });
  process.stdout.write(`${[header, ...body].join('\n')}\n`);
};

const argumentsByName = new Map();
const flags = new Set();
const main = () => {
  for (let index = 2; index < process.argv.length; index += 1) {
    const argument = process.argv[index];
    if (!argument.startsWith('--')) continue;
    const next = process.argv[index + 1];
    if (next && !next.startsWith('--')) {
      argumentsByName.set(argument, next);
      index += 1;
    } else flags.add(argument);
  }
  if (argumentsByName.has('--range')) printRange(argumentsByName.get('--range'));
  else printAt(argumentsByName.get('--at') || 'HEAD');
};

if (process.argv[1] && pathToFileURL(process.argv[1]).href === import.meta.url) main();
