// G-F(1) (Web GUI audit 2026-09-30; OIPC-C04, R5): a cached comparison is
// reused only when every scientific input of its identity is unchanged
// (CO-02, CO-03). Each row changes one input and must change the cache key;
// display-only changes must keep it. A builder that drops an input from its
// identity fails the table (see the last test).
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import test from 'node:test';
import { pathToFileURL } from 'node:url';

import { buildLosatDerivedPayloadCachePayload } from '../../gbdraw/web/js/app/run-analysis.js';
import { prepareLosatSourceBatches } from '../../gbdraw/web/js/app/linear-sources.js';

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

const withRecord = (index, patch) => (input) => ({
  ...input,
  recordPayloads: input.recordPayloads.map((record, position) => (
    position === index ? { ...record, ...patch(record) } : record
  ))
});
const withPair = (patch) => (input) => ({
  ...input,
  pairPayloads: input.pairPayloads.map((pair) => ({ ...pair, ...patch(pair) }))
});
const set = (patch) => (input) => ({ ...input, ...patch });
const DISPLAY_ONLY = [
  ['match style', set({ pairwiseMatchStyle: 'curve' })],
  ['comparison height', set({ comparisonHeight: 90 })],
  ['disclosure state', set({ comparisonDisclosureOpen: true })],
  ['scale font size', set({ scaleFontSize: 27 })],
  ['ruler label font size', set({ rulerLabelFontSize: 11 })]
];

const derivedBaseline = {
  mode: 'orthogroup',
  maxHits: 5,
  bitscore: 50,
  evalue: '1e-5',
  identity: 70,
  alignmentLength: 0,
  orthogroupMembershipMode: 'anchor_core_v1',
  orthogroupMemberMaxHits: 5,
  recordPayloads: [
    { recordIndex: 0, proteinCacheKey: 'protein-a', runtimeBindingHash: 'runtime-a', displayBindingHash: 'display-a', viewTransform: { length: 100, reverse: false } },
    { recordIndex: 1, proteinCacheKey: 'protein-b', runtimeBindingHash: 'runtime-b', displayBindingHash: 'display-b', viewTransform: { length: 200, reverse: true } }
  ],
  pairPayloads: [{ pairIndex: 0, queryIndex: 0, subjectIndex: 1, cacheKey: 'pair-a-b' }]
};
const DERIVED_TABLE = {
  scientific: [
    ['the record set', (input) => ({ ...input, recordPayloads: input.recordPayloads.slice(0, 1) })],
    ['a record region length', withRecord(0, ({ viewTransform }) => ({ viewTransform: { ...viewTransform, length: 90 } }))],
    ['a record orientation', withRecord(0, ({ viewTransform }) => ({ viewTransform: { ...viewTransform, reverse: true } }))],
    ['a record protein set (translation table, Feature visibility)', withRecord(1, () => ({ proteinCacheKey: 'protein-b-gencode-4' }))],
    ['a record runtime binding', withRecord(1, () => ({ runtimeBindingHash: 'runtime-b2' }))],
    ['a record display binding', withRecord(1, () => ({ displayBindingHash: 'display-b2' }))],
    ['the record index of a payload', withRecord(1, () => ({ recordIndex: 2 }))],
    ['a pair endpoint', withPair(() => ({ queryIndex: 1, subjectIndex: 0 }))],
    ['a pair raw search', withPair(() => ({ cacheKey: 'pair-a-b-rerun' }))],
    ['the bit score threshold', set({ bitscore: 60 })],
    ['the E-value threshold', set({ evalue: '1e-6' })],
    ['the identity threshold', set({ identity: 71 })],
    ['the alignment length threshold', set({ alignmentLength: 30 })],
    ['the comparison mode', (input) => ({ ...input, mode: input.mode === 'pairwise' ? 'collinear' : 'pairwise' })],
    ['the orthogroup membership mode', set({ orthogroupMembershipMode: 'all_hits' })],
    ['the orthogroup member max hits', set({ orthogroupMemberMaxHits: 6 })]
  ],
  displayOnly: [
    ...DISPLAY_ONLY,
    ['record payload order', (input) => ({ ...input, recordPayloads: [...input.recordPayloads].reverse() })]
  ]
};

// The raw search identity builder is module-private; load a copy that exports it.
const loadRawIdentityBuilder = async () => {
  const runAnalysisPath = resolve('gbdraw/web/js/app/run-analysis.js');
  const runAnalysisUrl = pathToFileURL(runAnalysisPath);
  const source = (await readFile(runAnalysisPath, 'utf8'))
    .replace('const buildLosatCachePayload = ({', 'export const buildLosatCachePayload = ({')
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

test('derived comparison identity changes with every scientific input (G-F(1))', () => {
  for (const mode of ['orthogroup', 'pairwise']) {
    const table = mode === 'pairwise'
      ? { ...DERIVED_TABLE, scientific: [...DERIVED_TABLE.scientific.filter(([name]) => !name.startsWith('the orthogroup')), ['the Pairwise max hits', set({ maxHits: 6 })]] }
      : DERIVED_TABLE;
    assert.deepEqual(keySensitivityProblems(buildLosatDerivedPayloadCachePayload, { ...derivedBaseline, mode }, table), [], mode);
  }
});

test('raw LOSAT search identity changes with every scientific input (G-F(1))', async () => {
  const build = await loadRawIdentityBuilder();
  assert.deepEqual(keySensitivityProblems(build, proteinRawBaseline, PROTEIN_RAW_TABLE), [], 'protein');
  assert.deepEqual(keySensitivityProblems(build, nucleotideRawBaseline, NUCLEOTIDE_RAW_TABLE), [], 'nucleotide');
});

test('a source batch searches a different database when a record of it changes (G-F(1))', async () => {
  const sources = { a: { name: 'a.gb' }, b: { name: 'b.gb' } };
  const sequences = ['a-1', 'a-2', 'b-1'].map((uid) => ({ uid, gb: sources[uid[0]], gff: null, fasta: null }));
  const hashText = async (text) => createHash('sha256').update(text).digest('hex');
  const plan = (fastaFor) => prepareLosatSourceBatches({
    sequences, specs: [{ queryIndex: 0, subjectIndex: 2 }, { queryIndex: 1, subjectIndex: 2 }],
    getEntry: async (index) => ({ fasta: fastaFor(index) }), buildArgs: () => ['--task', 'blastn'], hashText, protein: false
  });
  const contexts = async (fastaFor) => (await plan(fastaFor)).batches.map(({ searchContext, query, subject }) => [searchContext, query.hash, subject.hash]);
  const baseline = await contexts((index) => `>r${index}\nACGT\n`);
  assert.notDeepEqual(await contexts((index) => `>r${index}\n${index === 1 ? 'TGCA' : 'ACGT'}\n`), baseline,
    'a reverse-complemented query record changes the batch identity');
  assert.notDeepEqual(await contexts((index) => `>r${index}\n${index === 2 ? 'ACG' : 'ACGT'}\n`), baseline,
    'a cropped subject record changes the batch identity');
  assert.deepEqual(await contexts((index) => `>r${index}\nACGT\n`), baseline, 'unchanged records keep the batch identity');
});

test('the sensitivity table flags a builder that drops an input or keys on display state', () => {
  const dropsOrientation = (input) => {
    const identity = buildLosatDerivedPayloadCachePayload(input);
    return { ...identity, records: identity.records.map(({ viewTransform, ...record }) => ({ ...record, length: viewTransform.length })) };
  };
  assert.deepEqual(keySensitivityProblems(dropsOrientation, derivedBaseline, DERIVED_TABLE),
    ['reused after a change of a record orientation']);
  const keysOnHeight = (input) => ({ ...buildLosatDerivedPayloadCachePayload(input), height: input.comparisonHeight });
  assert.deepEqual(keySensitivityProblems(keysOnHeight, derivedBaseline, DERIVED_TABLE),
    ['not reused after display-only comparison height']);
});
