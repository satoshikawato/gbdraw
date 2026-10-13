// R9 / G-J(2) (Web GUI audit 2026-09-30, P08 follow-up): the Web has one
// GenBank header reader, services/genbank-header.js, whose record IDs follow the
// Python loader (IN-02, IN-03). Another module that matches a header keyword
// (LOCUS, ACCESSION, VERSION, DEFINITION, ORGANISM) reads headers a second way
// and can drift from the loader. OUTSIDE_OWNER_BASELINE, the pre-guard sites,
// is a registered design-rule allowlist (tools/web-design-rule-guards.json, R9,
// kind `count-map`): a file may not exceed its entry, and adding or raising an
// entry is an authority-only change. It is now empty.
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import test from 'node:test';

import { countGenBankRecords } from '../../gbdraw/web/js/services/genbank-header.js';
import { compareCountBaseline, reportLowerable } from './helpers/count-baseline.mjs';

const WEB_ROOT = new URL('../../gbdraw/web/js/', import.meta.url);
const OWNER = 'services/genbank-header.js';
const HEADER_KEYWORD = /\b(?:LOCUS|ACCESSION|VERSION|DEFINITION|ORGANISM)\b/g;

// file -> header keyword sites outside the owner.
const OUTSIDE_OWNER_BASELINE = {};

// Code without comments, so prose that names a keyword is not a reader.
const stripComments = (source) => source
  .replace(/\/\*[\s\S]*?\*\//g, '')
  .replace(/(^|[^:\\'"`])\/\/.*$/gm, '$1');

export const headerKeywordSites = (source) => stripComments(source).match(HEADER_KEYWORD)?.length ?? 0;

const listModules = (url, prefix = '') => readdirSync(url, { withFileTypes: true }).flatMap((entry) => {
  const path = `${prefix}${entry.name}`;
  if (entry.isDirectory()) return listModules(new URL(`${entry.name}/`, url), `${path}/`);
  return entry.name.endsWith('.js') ? [path] : [];
});

// The files outside the owner with more header keyword sites than their entry.
export const headerReaderProblems = (sitesByFile, baseline = OUTSIDE_OWNER_BASELINE) => {
  const outside = Object.fromEntries(Object.entries(sitesByFile).filter(([file]) => file !== OWNER));
  const { above, lowerable } = compareCountBaseline(outside, baseline);
  return {
    problems: above.map(({ key, observed, recorded }) => (
      `${key}: ${observed} GenBank header keyword site(s), allowed ${recorded}; read headers through ${OWNER}`
    )),
    lowerable
  };
};

test('only services/genbank-header.js reads GenBank header keywords (G-J(2), R9)', (t) => {
  const sitesByFile = Object.fromEntries(listModules(WEB_ROOT).map((file) => [
    file, headerKeywordSites(readFileSync(new URL(file, WEB_ROOT), 'utf8'))
  ]));
  assert.ok(sitesByFile[OWNER] > 0, 'the owner still reads the header keywords');
  const { problems, lowerable } = headerReaderProblems(sitesByFile);
  reportLowerable(t, 'OUTSIDE_OWNER_BASELINE', lowerable);
  assert.deepEqual(problems, []);
});

test('the reader guard flags a second header reader and names a baseline that may be lowered', () => {
  const secondReader = "const id = text.match(/^ACCESSION\\s+(\\S+)/m)?.[1];\n";
  const proseOnly = '// ORGANISM lines are read by services/genbank-header.js\nconst DEFINITION_ROLES = 1;\n';
  assert.equal(headerKeywordSites(secondReader), 1);
  assert.equal(headerKeywordSites(proseOnly), 0, 'comments and longer identifiers are not readers');
  assert.match(headerReaderProblems({ 'app/new-reader.js': 1 }).problems.join('\n'), /app\/new-reader\.js: 1 GenBank header keyword site/);
  const baseline = { 'app/old-reader.js': 1 };
  assert.match(headerReaderProblems({ 'app/old-reader.js': 2 }, baseline).problems.join('\n'), /allowed 1/);
  assert.deepEqual(headerReaderProblems({ 'app/old-reader.js': 0 }, baseline), { problems: [], lowerable: ['app/old-reader.js 1 -> 0'] });
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
