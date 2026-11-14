import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';
import { headingsOf } from '../../tools/generate_norms_index.mjs';

// The Product Contract names each revision after its pull request and records it
// in the section it changes, and numbers each new record after the pull request
// that adds it, so pull requests that change different sections do not edit the
// same lines or take the same number. This test keeps the numbered revisions and
// the approved decision ID lists of Authority metadata where they stopped, keeps
// the header free of a current revision number, and fails a new record numbered
// in sequence after the last listed one. The two numbers below never rise.
// Record ids are built at run time: tools/generate_norms_index.mjs reads this
// file, and a literal one would enter the index.
const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..', '..');
const CONTRACT = 'docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md';
const RULE = 'docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md#recording-a-revision';
const LAST_NUMBERED_REVISION = 151;
const LAST_LISTED_RECORD = 128;
// Pull request numbers passed this long before records were numbered after them,
// so a record number between LAST_LISTED_RECORD and this one is a sequential one.
const FIRST_PULL_REQUEST_RECORD = 1000;
const REVISIONS_FIELD = [
  '- Revisions: Revisions 1 to 151 are recorded in this section and keep their',
  '  numbers. Each later revision is named after the pull request that makes it,',
  '  as `Revision #<number>`, and is recorded in the section it changes',
  '  ([Recording a revision](#recording-a-revision)).'
].join('\n');
const RECORD = (number) => `PD-OI-${number}`;
const RECORD_HEADING = new RegExp(`^${RECORD('\\d+')}: `);

const read = () => readFileSync(join(ROOT, CONTRACT), 'utf8');

// A top-level list item of Authority metadata with its indented continuation lines.
const metadataField = (lines, start) => {
  let end = start + 1;
  while (end < lines.length && lines[end].startsWith('  ')) end += 1;
  return lines.slice(start, end).join('\n');
};

const revisionProblems = (source) => {
  const lines = source.split('\n');
  const headings = headingsOf(source);
  const metadata = headings.find(({ level, text }) => level === 2 && text === 'Authority metadata');
  if (!metadata) return [`${CONTRACT} has no "## Authority metadata" section`];
  const metadataEnd = headings.find(({ level, line }) => level <= 2 && line > metadata.line)?.line ?? lines.length;
  const inMetadata = (index) => index > metadata.line && index < metadataEnd;
  const problems = [];
  const fields = lines.flatMap((line, index) => (line.startsWith('- Revisions:') ? [index] : []));
  if (fields.length !== 1 || !inMetadata(fields[0]) || metadataField(lines, fields[0]) !== REVISIONS_FIELD) {
    problems.push(`the Revisions field of Authority metadata is missing or changed; it stays as written (${RULE})`);
  }
  const listed = [];
  lines.forEach((line, index) => {
    if (/^- Contract revision:/.test(line)) {
      problems.push(`line ${index + 1} holds a current revision number; name the revision after its pull request (${RULE})`);
    }
    const numbered = line.match(/^(?:- |#+ )Revision (\d+)\b/);
    if (numbered && !(inMetadata(index) && Number(numbered[1]) <= LAST_NUMBERED_REVISION)) {
      problems.push(`line ${index + 1} numbers Revision ${numbered[1]}; name it "Revision #<pull request>" `
        + `and record it in the section it changes (${RULE})`);
    }
    if (inMetadata(index) && /^- (?:Additional )?[Aa]pproved decision IDs:/.test(line)) {
      listed.push(...(metadataField(lines, index).match(new RegExp(RECORD('\\d+'), 'g')) ?? []));
    }
  });
  listed.filter((id) => Number(id.slice(RECORD('').length)) > LAST_LISTED_RECORD).forEach((id) => {
    problems.push(`the approved decision ID lists name ${id}; a record added after them records its approval `
      + `in its own entry (${RULE})`);
  });
  headings.filter(({ level, text }) => level === 3 && RECORD_HEADING.test(text)).forEach(({ line, text }) => {
    const id = text.slice(0, text.indexOf(':'));
    const number = Number(id.slice(RECORD('').length));
    if (number > LAST_LISTED_RECORD && number < FIRST_PULL_REQUEST_RECORD) {
      problems.push(`line ${line + 1}: ${id} is numbered in sequence; number a record added after `
        + `${RECORD(LAST_LISTED_RECORD)} after its pull request (${RULE})`);
    }
  });
  // An entry heading is one level below the section it ends: a record (###) or a
  // `##` section other than the record list, whose `###` headings start records.
  headings.forEach((heading, index) => {
    if (!heading.text.startsWith('Revision #')) return;
    const parent = headings.slice(0, index).reverse().find(({ level }) => level < heading.level);
    const placed = /^Revision #[1-9]\d*$/.test(heading.text) && parent && heading.level === parent.level + 1
      && (parent.level === 3 ? RECORD_HEADING.test(parent.text) : parent.level === 2 && parent.text !== 'Product Decision records');
    if (!placed) {
      problems.push(`line ${heading.line + 1}: "${heading.text}" is not a \`#### Revision #<number>\` heading at the end `
        + `of a record or a \`### Revision #<number>\` heading at the end of a \`##\` section (${RULE})`);
    }
  });
  return problems;
};

// The text of a section that starts at the heading matching `matches`, with
// `lines` inserted before the next heading of its level or above.
const endSection = (source, matches, inserted) => {
  const lines = source.split('\n');
  const headings = headingsOf(source);
  const at = headings.findIndex(matches);
  const end = headings.slice(at + 1).find(({ level }) => level <= headings[at].level)?.line ?? lines.length;
  lines.splice(end, 0, ...inserted, '');
  return lines.join('\n');
};
const recordHeading = (number) => ({ level, text }) => level === 3 && text.startsWith(`${RECORD(number)}: `);
const entry = (heading) => [heading, '', 'What changes, the approval, and the route.'];

test('the Product Contract records no revision by number after the last numbered one', () => {
  assert.deepEqual(revisionProblems(read()), []);
});

test('a numbered revision, a current revision number, or a newly listed record id fails', () => {
  const source = read();
  const last = `- Revision ${LAST_NUMBERED_REVISION} changes:`;
  assert.ok(source.includes(last));
  const numbered = source.replace(last, `- Revision ${LAST_NUMBERED_REVISION + 1} changes: a new record.\n${last}`);
  assert.match(revisionProblems(numbered).join('\n'), new RegExp(`numbers Revision ${LAST_NUMBERED_REVISION + 1};`));
  const current = source.replace(REVISIONS_FIELD, `- Contract revision: \`${LAST_NUMBERED_REVISION + 1}\``);
  assert.deepEqual(revisionProblems(current).map((problem) => problem.split(';')[0]), [
    'the Revisions field of Authority metadata is missing or changed',
    'line 8 holds a current revision number'
  ]);
  const lastListed = `\`${RECORD(LAST_LISTED_RECORD)}\``;
  const listed = source.replace(`and ${lastListed}\n`, `${lastListed}, and \`${RECORD(LAST_LISTED_RECORD + 1)}\`\n`);
  assert.notEqual(listed, source);
  assert.match(revisionProblems(listed).join('\n'), new RegExp(`name ${RECORD(LAST_LISTED_RECORD + 1)};`));
  const numberedEntry = endSection(source, recordHeading('001'), entry(`#### Revision ${LAST_NUMBERED_REVISION + 1}`));
  assert.match(revisionProblems(numberedEntry).join('\n'), new RegExp(`numbers Revision ${LAST_NUMBERED_REVISION + 1};`));
});

test('a new record is numbered after its pull request, not in sequence', () => {
  const source = read();
  const records = ({ level, text }) => level === 2 && text === 'Product Decision records';
  const added = (number) => endSection(source, records, [
    ...entry(`### ${RECORD(number)}: A new record`), '', ...entry(`#### Revision #${number}`)
  ]);
  assert.deepEqual(revisionProblems(added(9001)), []);
  [String(LAST_LISTED_RECORD + 1), `0${LAST_LISTED_RECORD + 1}`].forEach((number) => {
    assert.match(revisionProblems(added(number)).join('\n'), new RegExp(`${RECORD(number)} is numbered in sequence;`));
  });
});

test('an entry heading sits one level below the record or section it ends', () => {
  const source = read();
  assert.deepEqual(revisionProblems(endSection(source, recordHeading('001'), entry('#### Revision #9001'))), []);
  const metadata = ({ level, text }) => level === 2 && text === 'Authority metadata';
  assert.deepEqual(revisionProblems(endSection(source, metadata, entry('### Revision #9001'))), []);
  [
    endSection(source, recordHeading('001'), entry('### Revision #9001')),
    endSection(source, metadata, entry('#### Revision #9001')),
    endSection(source, recordHeading('001'), entry('#### Revision #the pull request'))
  ].forEach((changed) => assert.equal(revisionProblems(changed).length, 1, revisionProblems(changed).join('\n')));
});

test('pull requests that revise two neighbouring records merge without a conflict', () => {
  const source = read();
  const fieldOf = (number) => {
    const record = headingsOf(source).find(recordHeading(number));
    const lines = source.split('\n');
    const index = lines.findIndex((line, at) => at > record.line && line.startsWith('- Scenario revision: '));
    return { index, line: lines[index] };
  };
  // Each pull request raises its record's scenario revision and ends the record with its entry.
  const revise = (number, pullRequest) => {
    const { index, line } = fieldOf(number);
    const lines = endSection(source, recordHeading(number), entry(`#### Revision #${pullRequest}`)).split('\n');
    lines[index] = `${line} Revised by #${pullRequest}.`;
    return lines.join('\n');
  };
  const directory = mkdtempSync(join(tmpdir(), 'contract-revisions-'));
  try {
    const [ours, original, theirs] = ['ours', 'base', 'theirs'].map((name) => join(directory, name));
    writeFileSync(ours, revise('001', 9001));
    writeFileSync(original, source);
    writeFileSync(theirs, revise('002', 9002));
    const merged = spawnSync('git', ['merge-file', '-p', ours, original, theirs], { encoding: 'utf8', maxBuffer: 1 << 26 });
    assert.equal(merged.status, 0, `conflicting hunks: ${merged.status}`);
    assert.deepEqual(revisionProblems(merged.stdout), []);
    assert.equal(merged.stdout.split('\n#### Revision #900').length, 3);
  } finally {
    rmSync(directory, { recursive: true, force: true });
  }
});
