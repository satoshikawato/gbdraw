// R9 / G-J(2) (Web GUI audit 2026-09-30, P08 follow-up): the Web has one
// GenBank header reader, app/genbank-header.js, whose record IDs follow the
// Python loader (IN-02, IN-03). Another module that matches a header keyword
// (LOCUS, ACCESSION, VERSION, DEFINITION, ORGANISM) reads headers a second way
// and can drift from the loader. The baseline of pre-guard sites may only
// shrink; it is now empty.
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import test from 'node:test';

import { countGenBankRecords } from '../../gbdraw/web/js/app/genbank-header.js';

const WEB_ROOT = new URL('../../gbdraw/web/js/', import.meta.url);
const OWNER = 'app/genbank-header.js';
const HEADER_KEYWORD = /\b(?:LOCUS|ACCESSION|VERSION|DEFINITION|ORGANISM)\b/g;

// Shrink-only: file -> header keyword sites outside the owner.
const OUTSIDE_OWNER_BASELINE = {};

// Code without comments, so prose that names a keyword is not a reader.
const stripComments = (source) => source
  .replace(/\/\*[\s\S]*?\*\//g, '')
  .replace(/(^|[^:\\'"`])\/\/.*$/gm, '$1');

export const headerKeywordSites = (source) => stripComments(source).match(HEADER_KEYWORD)?.length ?? 0;

const listModules = (url, prefix = '') => readdirSync(url, { withFileTypes: true }).flatMap((entry) => {
  const path = `${prefix}${entry.name}`;
  if (entry.isDirectory()) return listModules(new URL(`${entry.name}/`, url), `${path}/`);
  return entry.name.endsWith('.js') && !path.endsWith('.generated.js') ? [path] : [];
});

export const headerReaderProblems = (sitesByFile, baseline = OUTSIDE_OWNER_BASELINE) => {
  const problems = [];
  for (const [file, sites] of Object.entries(sitesByFile)) {
    if (file === OWNER || sites === 0) continue;
    const allowed = baseline[file] ?? 0;
    if (sites > allowed) problems.push(`${file}: ${sites} GenBank header keyword site(s), allowed ${allowed}; read headers through ${OWNER}`);
  }
  for (const [file, allowed] of Object.entries(baseline)) {
    const sites = sitesByFile[file] ?? 0;
    if (sites < allowed) problems.push(`${file}: ${sites} site(s) left; lower its baseline to ${sites}`);
  }
  return problems;
};

test('only app/genbank-header.js reads GenBank header keywords (G-J(2), R9)', () => {
  const sitesByFile = Object.fromEntries(listModules(WEB_ROOT).map((file) => [
    file, headerKeywordSites(readFileSync(new URL(file, WEB_ROOT), 'utf8'))
  ]));
  assert.ok(sitesByFile[OWNER] > 0, 'the owner still reads the header keywords');
  assert.deepEqual(headerReaderProblems(sitesByFile), []);
});

test('the reader guard flags a second header reader and a stale baseline', () => {
  const secondReader = "const id = text.match(/^ACCESSION\\s+(\\S+)/m)?.[1];\n";
  const proseOnly = '// ORGANISM lines are read by app/genbank-header.js\nconst DEFINITION_ROLES = 1;\n';
  assert.equal(headerKeywordSites(secondReader), 1);
  assert.equal(headerKeywordSites(proseOnly), 0, 'comments and longer identifiers are not readers');
  assert.match(headerReaderProblems({ 'app/new-reader.js': 1 }).join('\n'), /app\/new-reader\.js: 1 GenBank header keyword site/);
  const baseline = { 'app/old-reader.js': 1 };
  assert.match(headerReaderProblems({ 'app/old-reader.js': 2 }, baseline).join('\n'), /allowed 1/);
  assert.match(headerReaderProblems({ 'app/old-reader.js': 0 }, baseline).join('\n'), /lower its baseline to 0/);
});

test('record counts follow the loader: a record starts at "LOCUS" plus seven spaces', () => {
  const record = 'LOCUS       GOOD1                    10 bp    DNA     linear   UNK 01-JAN-2000\nORIGIN\n        1 acgtacgtac\n//\n';
  // Biopython SeqIO.parse reads only GOOD1 from each of these texts.
  assert.equal(countGenBankRecords(record), 1);
  assert.equal(countGenBankRecords(`${record}LOCUS AB1 100 bp DNA\n//\n`), 1);
  assert.equal(countGenBankRecords(`${record}LOCUS\tAB1 100 bp DNA\n//\n`), 1);
  assert.equal(countGenBankRecords(`${record}LOCUS   AB1                    10 bp    DNA\n//\n`), 1);
  assert.equal(countGenBankRecords(`${record}${record.replace('GOOD1', 'GOOD2')}`.replace(/\n/g, '\r\n')), 2);
  assert.equal(countGenBankRecords(`${record}${record}`), 2);
  assert.equal(countGenBankRecords(''), 0);
});
