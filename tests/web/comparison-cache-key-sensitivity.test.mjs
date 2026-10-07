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
  pairPayloads: [{ pairIndex: 0, queryIndex: 0, subjectIndex: 1, cacheKey: 'pair-a-b', displayPair: true }]
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
    ['the displayed direction of a pair (displayPair)', withPair(() => ({ displayPair: false }))],
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

const collinearBaseline = {
  ...derivedBaseline,
  mode: 'collinear',
  collinearMinAnchors: 3,
  collinearMaxUnitGap: 2,
  collinearUnitMode: 'auto',
  collinearColorMode: 'orientation',
  collinearAnchorMode: 'rbh',
  collinearMergeOrientation: 'either',
  collinearMaxDiagonalDrift: 5,
  collinearMaxConflictsInMergeGap: 1,
  collinearMaxParalogLinksPerOrthogroup: 2,
  collinearInferOrthogroups: true,
  collinearSearchScope: 'adjacent'
};
const COLLINEAR_TABLE = {
  scientific: [
    ...DERIVED_TABLE.scientific.filter(([name]) => name !== 'the comparison mode'),
    ['the comparison mode', set({ mode: 'orthogroup' })],
    ['the minimum anchors', set({ collinearMinAnchors: 4 })],
    ['the maximum unit gap', set({ collinearMaxUnitGap: 3 })],
    ['the unit mode', set({ collinearUnitMode: 'gene' })],
    ['the color mode', set({ collinearColorMode: 'block' })],
    ['the anchor mode', set({ collinearAnchorMode: 'all' })],
    ['the merge orientation', set({ collinearMergeOrientation: 'same' })],
    ['the maximum diagonal drift', set({ collinearMaxDiagonalDrift: 6 })],
    ['the maximum conflicts in a merge gap', set({ collinearMaxConflictsInMergeGap: 2 })],
    ['the maximum paralog links per orthogroup', set({ collinearMaxParalogLinksPerOrthogroup: 3 })],
    ['orthogroup inference', set({ collinearInferOrthogroups: false })],
    ['the search scope', set({ collinearSearchScope: 'all' })],
    // Only the Collinear converter reads it: with the search scope 'all' it limits the
    // output to the displayed pairs of a CLI grid row layout.
    ['explicit display pairs (CLI grid row layout)', set({ explicitDisplayPairs: true })]
  ],
  displayOnly: DERIVED_TABLE.displayOnly
};

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

test('derived comparison identity changes with every scientific input (G-F(1))', () => {
  for (const mode of ['orthogroup', 'pairwise']) {
    const table = mode === 'pairwise'
      ? { ...DERIVED_TABLE, scientific: [...DERIVED_TABLE.scientific.filter(([name]) => !name.startsWith('the orthogroup')), ['the Pairwise max hits', set({ maxHits: 6 })]] }
      : DERIVED_TABLE;
    assert.deepEqual(keySensitivityProblems(buildLosatDerivedPayloadCachePayload, { ...derivedBaseline, mode }, table), [], mode);
  }
  assert.deepEqual(keySensitivityProblems(buildLosatDerivedPayloadCachePayload, collinearBaseline, COLLINEAR_TABLE), [], 'collinear');
});

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
    getEntry: async (index) => ({ fasta: fastaFor(index) }), buildArgs: () => ['--task', 'blastn'], hashText, protein: false
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
  const dropsDisplayDirection = (input) => buildLosatDerivedPayloadCachePayload({
    ...input,
    pairPayloads: input.pairPayloads.map(({ displayPair, ...pair }) => pair)
  });
  assert.deepEqual(keySensitivityProblems(dropsDisplayDirection, derivedBaseline, DERIVED_TABLE),
    ['reused after a change of the displayed direction of a pair (displayPair)']);
  const dropsExplicitDisplayPairs = (input) => buildLosatDerivedPayloadCachePayload({ ...input, explicitDisplayPairs: undefined });
  assert.deepEqual(keySensitivityProblems(dropsExplicitDisplayPairs, collinearBaseline, COLLINEAR_TABLE),
    ['reused after a change of explicit display pairs (CLI grid row layout)']);
});
