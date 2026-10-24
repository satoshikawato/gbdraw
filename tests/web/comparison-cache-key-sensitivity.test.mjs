// G-F(1) (Web GUI audit 2026-09-30; OIPC-C04, R5): a cached comparison is
// reused only when every scientific input of its identity is unchanged
// (CO-02, CO-03). Each row changes one input and must change the cache key;
// display-only changes must keep it. A builder that drops an input from its
// identity fails the table (see the last test). The derived protein payload
// identity has its Python owner and its table in
// tests/test_derived_losatp_identity.py (OV-1063).
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import test from 'node:test';
import { pathToFileURL } from 'node:url';

import { prepareLosatSourceBatches } from '../../gbdraw/web/js/services/linear-sources.js';

const keyOf = (identity) => createHash('sha256').update(JSON.stringify(identity)).digest('hex');

// Problems of `build` against the table: a row that must change the key and
// does not, or a display-only row that changes it.
export const keySensitivityProblems = (build, baseline, { scientific, displayOnly }) => {
  const baselineKey = keyOf(build(baseline));
  const problems = [];
  for (const [name, change] of scientific) {
    if (keyOf(build(change(baseline))) === baselineKey) problems.push(`reused after a change of ${name}`);
  }
  for (const [name, change] of displayOnly) {
    if (keyOf(build(change(baseline))) !== baselineKey) problems.push(`not reused after display-only ${name}`);
  }
  return problems;
};

const set = (patch) => (input) => ({ ...input, ...patch });
const DISPLAY_ONLY = [
  ['match style', set({ pairwiseMatchStyle: 'curve' })],
  ['comparison height', set({ comparisonHeight: 90 })],
  ['disclosure state', set({ comparisonDisclosureOpen: true })],
  ['scale font size', set({ scaleFontSize: 27 })],
  ['ruler label font size', set({ rulerLabelFontSize: 11 })]
];

// Load a copy of the module so the probe does not share module state.
const loadRawIdentityBuilder = async () => {
  const runAnalysisPath = resolve('gbdraw/web/js/app/run-analysis.js');
  const runAnalysisUrl = pathToFileURL(runAnalysisPath);
  const source = (await readFile(runAnalysisPath, 'utf8'))
    .replace(/from '(\.\.?\/[^']+)'/g, (_match, specifier) => `from '${new URL(specifier, runAnalysisUrl).href}'`);
  assert.match(source, /export const buildLosatCachePayload = \(\{/, 'the raw identity builder is still named buildLosatCachePayload');
  const probePath = join(await mkdtemp(join(tmpdir(), 'gbdraw-key-sensitivity-')), 'run-analysis-probe.mjs');
  await writeFile(probePath, source, 'utf8');
  return (await import(pathToFileURL(probePath))).buildLosatCachePayload;
};

const proteinRawBaseline = {
  identityKind: 'protein', program: 'blastp', outfmt: '6', args: ['--max-target-seqs', '9', '--query-gencode', '11'],
  queryProteinSetHash: 'query-proteins', subjectProteinSetHash: 'subject-proteins',
  queryRuntimeBindingHash: 'query-runtime', subjectRuntimeBindingHash: 'subject-runtime',
  queryRecordInstanceKey: 'query-record', subjectRecordInstanceKey: 'subject-record', searchContext: 'database-a'
};
const PROTEIN_RAW_TABLE = {
  scientific: [
    ['the search arguments (translation table)', set({ args: ['--max-target-seqs', '9', '--query-gencode', '4'] })],
    ['the search arguments (max target seqs)', set({ args: ['--max-target-seqs', '10', '--query-gencode', '11'] })],
    ['the program', set({ program: 'tblastx' })],
    ['the query protein set', set({ queryProteinSetHash: 'query-proteins-2' })],
    ['the subject protein set', set({ subjectProteinSetHash: 'subject-proteins-2' })],
    ['the query runtime binding', set({ queryRuntimeBindingHash: 'query-runtime-2' })],
    ['the subject runtime binding', set({ subjectRuntimeBindingHash: 'subject-runtime-2' })],
    ['the query record instance', set({ queryRecordInstanceKey: 'query-record-2' })],
    ['the subject record instance', set({ subjectRecordInstanceKey: 'subject-record-2' })],
    ['the searched database', set({ searchContext: 'database-b' })]
  ],
  displayOnly: DISPLAY_ONLY
};
const nucleotideRawBaseline = {
  identityKind: 'nucleotide', program: 'blastn', outfmt: '6', args: ['--evalue', '1e-5'],
  queryCanonicalHash: 'query-sequence', subjectCanonicalHash: 'subject-sequence', flow: 'pairwise', searchContext: 'database-a'
};
const NUCLEOTIDE_RAW_TABLE = {
  scientific: [
    ['the search arguments', set({ args: ['--evalue', '1e-6'] })],
    ['the program', set({ program: 'tblastx' })],
    ['the query sequence (region, orientation)', set({ queryCanonicalHash: 'query-sequence-reversed' })],
    ['the subject sequence (region, orientation)', set({ subjectCanonicalHash: 'subject-sequence-cropped' })],
    ['the search flow', set({ flow: 'collinear' })],
    ['the searched database', set({ searchContext: 'database-b' })]
  ],
  displayOnly: DISPLAY_ONLY
};

test('raw LOSAT search identity changes with every scientific input (G-F(1))', async () => {
  const build = await loadRawIdentityBuilder();
  assert.deepEqual(keySensitivityProblems(build, proteinRawBaseline, PROTEIN_RAW_TABLE), [], 'protein');
  assert.deepEqual(keySensitivityProblems(build, nucleotideRawBaseline, NUCLEOTIDE_RAW_TABLE), [], 'nucleotide');
});

test('a source batch searches a different database when a record of it changes (G-F(1))', async () => {
  const sources = { a: { name: 'a.gb' }, b: { name: 'b.gb' } };
  const hashText = async (text) => createHash('sha256').update(text).digest('hex');
  const contexts = async (uids, specs, fastaFor) => (await prepareLosatSourceBatches({
    sequences: uids.map((uid) => ({ uid, gb: sources[uid[0]], gff: null, fasta: null })), specs,
    getEntry: async (index) => ({ fasta: fastaFor(index) }), getDigest: async (index) => hashText(fastaFor(index)), buildArgs: () => ['--task', 'blastn'], hashText, protein: false
  })).batches.map(({ scope, searchContext, query, subject }) => [scope, searchContext, query.hash, subject.hash]);
  const plain = (index) => `>r${index}\nACGT\n`;
  const changed = (changedIndex, text) => (index) => `>r${index}\n${index === changedIndex ? text : 'ACGT'}\n`;

  // Between two source files: both records of the query file search the subject file.
  const between = [{ queryIndex: 0, subjectIndex: 2 }, { queryIndex: 1, subjectIndex: 2 }];
  const betweenBaseline = await contexts(['a-1', 'a-2', 'b-1'], between, plain);
  assert.notDeepEqual(await contexts(['a-1', 'a-2', 'b-1'], between, changed(1, 'TGCA')), betweenBaseline,
    'a reverse-complemented query record changes the batch identity');
  assert.notDeepEqual(await contexts(['a-1', 'a-2', 'b-1'], between, changed(2, 'ACG')), betweenBaseline,
    'a cropped subject record changes the batch identity');
  assert.deepEqual(await contexts(['a-1', 'a-2', 'b-1'], between, plain), betweenBaseline, 'unchanged records keep the batch identity');

  // Within one source file: the database is the file without the query record.
  const within = [{ queryIndex: 0, subjectIndex: 1 }];
  const withinBaseline = await contexts(['a-1', 'a-2', 'a-3'], within, plain);
  assert.equal(withinBaseline[0][0], 'within-source');
  assert.notDeepEqual(await contexts(['a-1', 'a-2', 'a-3'], within, changed(2, 'ACG')), withinBaseline,
    'a cropped record of the searched database changes the batch identity');
  assert.notDeepEqual(await contexts(['a-1', 'a-2', 'a-3'], within, changed(0, 'TGCA')), withinBaseline,
    'a reversed query record changes the batch identity');

  // A requested self comparison searches the record alone, a different database.
  const self = await contexts(['a-1', 'a-2', 'a-3'], [{ queryIndex: 0, subjectIndex: 0 }], plain);
  assert.equal(self[0][0], 'self');
  assert.notEqual(self[0][3], withinBaseline[0][3], 'the self database differs from the within-source database');
});

test('the sensitivity table flags a builder that drops an input or keys on display state', async () => {
  const build = await loadRawIdentityBuilder();
  const dropsSearchContext = (input) => {
    const { searchContext: _searchContext, ...identity } = build(input);
    return identity;
  };
  assert.deepEqual(keySensitivityProblems(dropsSearchContext, proteinRawBaseline, PROTEIN_RAW_TABLE),
    ['reused after a change of the searched database']);
  const keysOnHeight = (input) => ({ ...build(input), height: input.comparisonHeight });
  assert.deepEqual(keySensitivityProblems(keysOnHeight, proteinRawBaseline, PROTEIN_RAW_TABLE),
    ['not reused after display-only comparison height']);
});
