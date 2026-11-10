#!/usr/bin/env node
// Writes docs/internal/NORMS_INDEX.md from the owners of each fact:
//   - the design rules, their guard tests and the declared pull request job
//     (tools/web-design-rule-guards.json), each rule's heading in its document;
//   - the CI jobs, from the CI routing policy (tools/ci-impact-policy.mjs): the
//     jobs a pull request that changes only the guard requires, which must
//     include the declared job, and the dev jobs that run a test of the guard's
//     kind and pytest markers (`leafTestKind`, `leafJobsFor`, with the markers
//     that tools/ci-impact.mjs reads);
//   - the Product Contract records (their headings and fields) and the test
//     files that cite each record id;
//   - every `<document>.md#<anchor>` that a tool or a test cites: the anchor must
//     be a heading of that document when the document exists.
// It fails when a declared job is not required or an anchor of an existing
// document is not one of its headings.
// Run with --check to fail when the written file is stale;
// tests/ci/norms-index.test.mjs runs it in the CI impact job, which every pull
// request and dev push runs.
import { execFileSync } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import {
  classifyChanges,
  leafJobsFor,
  leafTestKind,
  requiredJobsFor,
  requiresFullCoverage
} from './ci-impact-policy.mjs';
import { pythonMarkers } from './ci-impact.mjs';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
export const OUTPUT = 'docs/internal/NORMS_INDEX.md';
const REGISTRY = 'tools/web-design-rule-guards.json';
const CONTRACT = 'docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md';
const TEXT_FILE = /\.(?:cjs|js|json|md|mjs|py|ts|tsv|txt)$/;
const RECORD_ID = /PD-OI-\d+/g;
const RECORD_FIELDS = ['Concern key', 'Scenario revision', 'Supersedes', 'Status'];

// Tracked and new (not ignored) files of the working tree, so that a test written
// before `git add` is indexed as it will be committed.
const repositoryFiles = (root) => execFileSync(
  'git', ['ls-files', '-z', '--cached', '--others', '--exclude-standard'],
  { cwd: root, encoding: 'utf8', maxBuffer: 64 * 1024 * 1024 }
).split('\0').filter((file) => file && fs.existsSync(path.join(root, file))).sort();

// GitHub's heading anchor: the rendered text, lower case, without punctuation,
// with spaces as hyphens; a repeated anchor gets -1, -2, ...
export const slugOf = (text) => text
  .replace(/`([^`]*)`/g, '$1')
  .replace(/!?\[([^\]]*)\]\([^)]*\)/g, '$1')
  .replace(/<[^>]+>/g, '')
  .replace(/(^|[^\p{L}\p{N}])_+([^_]+?)_+(?=[^\p{L}\p{N}]|$)/gu, '$1$2')
  .toLowerCase()
  .replace(/[^\p{L}\p{M}\p{N}\p{Pc}\- ]/gu, '')
  .replace(/ /g, '-');

// Headings outside fenced code, in order: { level, text, slug, line }. A fence
// closes on a bare run of its own character at least as long as the opening one.
export const headingsOf = (source) => {
  const seen = new Map();
  const headings = [];
  let fence = null;
  source.split('\n').forEach((line, index) => {
    const marker = line.match(/^ {0,3}(`{3,}|~{3,})(.*)$/);
    if (fence !== null) {
      if (marker && marker[1][0] === fence[0] && marker[1].length >= fence.length && !marker[2].trim()) fence = null;
      return;
    }
    if (marker && !(marker[1][0] === '`' && marker[2].includes('`'))) {
      fence = marker[1];
      return;
    }
    const heading = fence === null && line.match(/^ {0,3}(#{1,6})[ \t]+(.*?)(?:[ \t]+#+)?[ \t]*$/);
    if (!heading) return;
    const base = slugOf(heading[2]);
    const count = seen.get(base) ?? 0;
    seen.set(base, count + 1);
    headings.push({ level: heading[1].length, text: heading[2], slug: count ? `${base}-${count}` : base, line: index });
  });
  return headings;
};

const code = (value) => `\`${value}\``;
const cell = (value) => String(value).replace(/\|/g, '\\|');
const link = (from, file, slug, text) => `[${text}](${path.posix.relative(path.posix.dirname(from), file)}#${slug})`;
const table = (header, rows) => [
  `| ${header.join(' | ')} |`,
  `| ${header.map(() => '---').join(' | ')} |`,
  ...rows.map((row) => `| ${row.map(cell).join(' | ')} |`)
];

// The jobs the pull request plan requires when a pull request changes only `file`
// (tools/ci-impact.mjs: a full-coverage capability takes the full tier).
const pullRequestJobsFor = (file) => {
  const { impact, capabilities } = classifyChanges([{ status: 'M', paths: [file] }]);
  const decision = requiresFullCoverage('pr', capabilities) ? 'full' : 'selective';
  return requiredJobsFor({ profile: 'pr', impact, decision, capabilities });
};

// The dev jobs that run a test of this file's kind and pytest markers, as the
// planner routes a changed leaf test; null for a fixture or helper.
const devJobsFor = (file, read) => {
  const kind = leafTestKind(file);
  if (kind === null) return null;
  return leafJobsFor({ profile: 'dev', kind, markers: kind === 'python' ? pythonMarkers(read(file)) : null });
};

const ruleSection = ({ registry, headingSlug, read, problems }) => {
  const guardRows = [];
  const allowlistRows = [];
  registry.rules.forEach((rule) => {
    const document = rule.rulesPath ?? registry.rulesPath;
    const slug = headingSlug(document, rule.heading);
    if (slug === null) problems.push(`${rule.id}: ${document} has no heading "${rule.heading}"`);
    const ruleCell = slug === null ? rule.id : link(OUTPUT, document, slug, rule.id);
    rule.guards.forEach((guard) => {
      const required = pullRequestJobsFor(guard.path);
      if (!required.includes(guard.job)) {
        problems.push(`${rule.id} ${guard.path}: the registry declares ${guard.job}, `
          + `but a pull request that changes only this file requires ${required.join(', ')}`);
      }
      const devJobs = devJobsFor(guard.path, read);
      guardRows.push([
        ruleCell, code(guard.path), code(guard.job),
        devJobs === null ? 'not a test' : devJobs.map(code).join(', ') || 'none',
        guard.retireWhen
      ]);
    });
    rule.allowlists.forEach((allowlist) => {
      allowlistRows.push([ruleCell, code(allowlist.path), code(allowlist.symbol), allowlist.kind]);
    });
  });
  return [
    '## Design rules, guard tests, and CI jobs',
    '',
    `Source: ${code(REGISTRY)}. A rule links to its heading.`,
    '',
    '- Pull request job: the job the registry declares. This generator checks that a',
    '  pull request that changes only the guard requires that job',
    `  (${code('tools/ci-impact-policy.mjs')}), and`,
    `  ${code('tests/web/architecture-contracts.test.mjs')} checks that the job's`,
    '  commands run the guard.',
    '- Dev jobs: the jobs of the dev tier that run a test of that kind and its pytest',
    `  markers (${code('leafTestKind')}, ${code('leafJobsFor')}). A dev push may take a job's`,
    '  result from the pull request it merges instead (`docs/internal/SELECTIVE_CI.md`,',
    '  "Merged pull request inheritance"). A fixture or helper runs with the guards that',
    '  read it.',
    '',
    ...table(['Rule', 'Guard', 'Pull request job', 'Dev jobs', 'Retire when'], guardRows),
    '',
    '### Allowlists',
    '',
    ...table(['Rule', 'File', 'Symbol', 'Kind'], allowlistRows),
    ''
  ];
};

// Records start at `### PD-OI-nnn: Title` and end at the next heading of level 1-3.
const contractRecords = (source) => {
  const lines = source.split('\n');
  const headings = headingsOf(source);
  const records = new Map();
  headings.forEach((heading, index) => {
    const id = heading.level === 3 && heading.text.match(/^(PD-OI-\d+): (.+)$/);
    if (!id) return;
    const end = headings.slice(index + 1).find((next) => next.level <= 3)?.line ?? lines.length;
    const fields = {};
    lines.slice(heading.line + 1, end).forEach((line) => {
      const field = line.match(/^- ([A-Z][a-z]+(?: [a-z]+)*): (.*)$/);
      if (field && RECORD_FIELDS.includes(field[1]) && !(field[1] in fields)) fields[field[1]] = field[2];
    });
    records.set(id[1], { title: id[2], slug: heading.slug, body: lines.slice(heading.line + 1, end).join('\n'), fields });
  });
  return records;
};

const supersedes = (id, value) => {
  if (value === undefined || /^none\b/i.test(value)) return '—';
  const earlier = value.match(/^`(PD-OI-\d+)`, scenario revision `(\d+)`/);
  if (!earlier) return value.replace(/[.,]$/, '');
  return earlier[1] === id ? `scenario revision ${earlier[2]}` : `${earlier[1]} scenario revision ${earlier[2]}`;
};
const unquoted = (value) => (value ?? '—').replace(/^`([^`]*)`$/, '$1');

const contractSection = ({ read, files, problems }) => {
  const records = contractRecords(read(CONTRACT));
  if (records.size === 0) problems.push(`${CONTRACT} has no PD-OI record heading`);
  const citing = new Map();
  files.filter((file) => file.startsWith('tests/') && TEXT_FILE.test(file)).forEach((file) => {
    new Set(read(file).match(RECORD_ID) ?? []).forEach((id) => citing.set(id, [...(citing.get(id) ?? []), file]));
  });
  const recordLink = (id) => link(OUTPUT, CONTRACT, records.get(id).slug, id);
  const byNumber = (left, right) => Number(left.slice(6)) - Number(right.slice(6)) || (left < right ? -1 : 1);
  const fileList = (id) => (citing.get(id) ?? []).map((file) => `  - ${code(file)}`);
  // One block per record, ended by a blank line, with no count that a citing
  // test elsewhere changes: two pull requests that change different records edit
  // lines that an unchanged line separates, so their changes to this file merge.
  const recordBlocks = [...records.keys()].sort(byNumber).flatMap((id) => {
    const { title, fields } = records.get(id);
    const earlier = supersedes(id, fields.Supersedes);
    const facts = [
      unquoted(fields.Status), `concern ${code(unquoted(fields['Concern key']))}`,
      `scenario revision ${unquoted(fields['Scenario revision'])}`,
      ...(earlier === '—' ? [] : [`supersedes ${earlier}`])
    ];
    return [`- ${recordLink(id)}: ${title}. ${facts.join(', ')}.`, ...fileList(id), ''];
  });
  const missing = [...citing.keys()].filter((id) => !records.has(id)).sort(byNumber);
  return [
    '## Product Contract records and the tests that cite them',
    '',
    `Source: the \`### PD-OI-nnn\` records of ${code(CONTRACT)}`,
    `and the test files under ${code('tests/')} that name a record id.`,
    '',
    'Each record gives its title, status, concern key, scenario revision, and the',
    'revision it supersedes, followed by the test files that cite it. A record is',
    'revised in place: its Supersedes field names its own earlier scenario revision.',
    'A record replaced in full is removed, so the tests that still cite it are listed',
    'under "Cited ids without a record".',
    '',
    ...recordBlocks,
    '### Cited ids without a record',
    '',
    ...(missing.length ? missing.flatMap((id) => {
      const mentions = [...records.keys()].filter((other) => records.get(other).body.includes(id)).sort(byNumber);
      return [
        `- ${code(id)}: no record. ${mentions.length ? `Records that mention it: ${mentions.map(recordLink).join(', ')}.` : 'No record mentions it.'}`,
        ...fileList(id)
      ];
    }) : ['None.']),
    ''
  ];
};

// `<document>.md#<anchor>` written in a tool or a test, outside a URL.
export const citedAnchors = (source) => [...source.matchAll(/\.md#([\p{L}\p{M}\p{N}\p{Pc}-]+)/gu)].flatMap((match) => {
  let start = match.index;
  while (start > 0 && /[\w./-]/.test(source[start - 1])) start -= 1;
  const document = `${source.slice(start, match.index)}.md`;
  return document.startsWith('//') || document === '.md' ? [] : [{ document, slug: match[1] }];
});

const anchorSection = ({ read, files, headingSlugs, problems }) => {
  const markdown = files.filter((file) => file.endsWith('.md'));
  const resolve = (citer, document) => {
    const candidates = [path.posix.normalize(document), path.posix.join(path.posix.dirname(citer), document)];
    const direct = candidates.find((candidate) => markdown.includes(candidate));
    if (direct) return direct;
    // A link that a generator writes into another document, such as `../page.md`.
    const tail = `/${document.replace(/^(?:\.\.?\/)+/, '')}`;
    const matches = markdown.filter((file) => `/${file}`.endsWith(tail));
    return matches.length === 1 ? matches[0] : null;
  };
  const anchors = new Map();
  const absent = new Map();
  files.filter((file) => /^(?:tests|tools)\//.test(file) && TEXT_FILE.test(file)).forEach((file) => {
    citedAnchors(read(file)).forEach(({ document, slug }) => {
      const target = resolve(file, document);
      if (target === null) {
        // Test data may name a document that does not exist; it is listed, not refused.
        const key = `${document}#${slug}`;
        absent.set(key, [...new Set([...(absent.get(key) ?? []), file])]);
        return;
      }
      const heading = headingSlugs(target).get(slug);
      if (heading === undefined) {
        problems.push(`${file} cites ${document}#${slug}: ${target} has no heading with that anchor`);
        return;
      }
      const key = `${target}#${slug}`;
      const entry = anchors.get(key) ?? { target, slug, heading, citers: new Set() };
      anchors.set(key, entry.citers.add(file) && entry);
    });
  });
  const rows = [...anchors.keys()].sort().map((key) => {
    const { target, slug, heading, citers } = anchors.get(key);
    return [link(OUTPUT, target, slug, `${target}#${slug}`), heading, [...citers].sort().map(code).join(', ')];
  });
  return [
    '## Document anchors cited by tools and tests',
    '',
    'Every `<document>.md#<anchor>` written in a file under `tools/` or `tests/` that',
    'names a document of the repository resolves to one of its headings; generation',
    'fails otherwise. Moving or renaming one of these headings needs the citing files',
    'changed with it.',
    '',
    ...table(['Anchor', 'Heading', 'Cited by'], rows),
    '',
    '### Cited documents that are not in the repository',
    '',
    ...(absent.size ? [...absent.keys()].sort().map((key) => `- ${code(key)}: ${absent.get(key).sort().map(code).join(', ')}`) : ['None.']),
    ''
  ];
};

// The index and the problems that make generation fail.
export const buildNormsIndex = ({
  read = (file) => fs.readFileSync(path.join(ROOT, file), 'utf8'),
  files = repositoryFiles(ROOT),
  registry = JSON.parse(read(REGISTRY))
} = {}) => {
  const problems = [];
  const headingCache = new Map();
  const headingsIn = (document) => headingCache.get(document)
    ?? headingCache.set(document, headingsOf(read(document))).get(document);
  const headingSlugs = (document) => new Map(headingsIn(document).map(({ slug, text }) => [slug, text]));
  const headingSlug = (document, text) => headingsIn(document).find((heading) => heading.text === text)?.slug ?? null;
  const lines = [
    '# Norms index',
    '',
    '<!-- Generated by tools/generate_norms_index.mjs; do not edit manually. -->',
    '',
    'Which test enforces each registered design rule and which CI jobs run it, the',
    'Product Contract records and the tests that cite them, and the document anchors',
    'that tools and tests cite. Run `node tools/generate_norms_index.mjs` after changing',
    'the guard registry, a rule heading, the Product Contract, a test that names a',
    'record id, or a cited heading; `tests/ci/norms-index.test.mjs` fails while this',
    'file is stale.',
    '',
    ...ruleSection({ registry, headingSlug, read, problems }),
    ...contractSection({ read, files, problems }),
    ...anchorSection({ read, files, headingSlugs, problems })
  ];
  return { text: `${lines.join('\n').replace(/\n+$/, '')}\n`, problems };
};

export const main = (argv = process.argv.slice(2)) => {
  const check = argv.includes('--check');
  const { text, problems } = buildNormsIndex();
  if (problems.length) {
    process.stderr.write(`${problems.map((problem) => `norms index: ${problem}`).join('\n')}\n`);
    return 1;
  }
  const output = path.join(ROOT, OUTPUT);
  if (check) {
    const current = fs.existsSync(output) ? fs.readFileSync(output, 'utf8') : null;
    if (current === text) return 0;
    process.stderr.write(`${OUTPUT} is stale: run \`node tools/generate_norms_index.mjs\` and commit the result.\n`);
    return 1;
  }
  fs.writeFileSync(output, text);
  return 0;
};

if (process.argv[1] && import.meta.url === pathToFileURL(path.resolve(process.argv[1])).href) {
  process.exitCode = main();
}
