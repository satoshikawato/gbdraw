import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import test from 'node:test';
import { RECORD_SOURCE_FILES } from '../../gbdraw/web/js/contracts/schema.generated.js';
import {
  INPUT_TYPES, inputSourceFiles, inputSourceSlotFiles, isRecordSource, recordSourceOf, recordSourceResourceIds,
  recordSourceFilesOf, recordSourceSlotIds, storedRecordSourceResourceIds
} from '../../gbdraw/web/js/services/input-types.js';
import {
  foldLoadedLinearSources,
  linearCardFileFields,
  linearCardInputComplete,
  groupLinearSourceRecords,
  isPristineLinearSource,
  linearRecordNames,
  linearSourceHasPrimaryInput,
  linearSourceDepthStatus,
  losatSearchSequences,
  moveLinearSourceGroup,
  nucleotideLosatArgs,
  planLinearSourceRemoval,
  planLosatSourceJobs,
  prepareLosatSourceBatches,
  proteinLosatArgs,
  restoreFoldedLinearSources,
  setLinearSourceFolded,
  splitLosatSourceResult
} from '../../gbdraw/web/js/services/linear-sources.js';
import { readFileSync } from 'node:fs';
import { buildLosatJobSpecs } from '../../gbdraw/web/js/services/linear-comparisons.js';
import {
  buildLosatCachePayload,
  extractLosatFastaFast
} from '../../gbdraw/web/js/app/run-analysis.js';
import { createTestInputStore } from './helpers/input-files.mjs';
import {
  adoptCurrentSessionResources,
  createSessionResourceFileView
} from '../../gbdraw/web/js/services/session-resource-backing.js';

// Source grouping keys on a ref's `sourceId`: one registered object is one source.
const inputFiles = createTestInputStore();
const file = (name) => inputFiles.registerFile({ name, size: 0 });
const hashText = async (value) => createHash('sha256').update(value).digest('hex');
/** The entry digest owner of run-analysis.js (`getSeqHash`): `hashText` of the entry's FASTA. */
const digestOf = (getEntry) => async (index) => hashText((await getEntry(index)).fasta);

const sourceRecord = (uid, gb) => ({ uid, gb, gff: null, fasta: null });
const sourceA = file('a.gb');
const sourceB = file('b.gb');
const sourceC = file('c.gb');
const singleSourceRecords = [
  sourceRecord('a-1', sourceA),
  sourceRecord('b-1', sourceB),
  sourceRecord('c-1', sourceC)
];
const singleSourceOrder = moveLinearSourceGroup(singleSourceRecords, 1, -1);
assert.deepEqual(singleSourceOrder.map(({ uid }) => uid), ['b-1', 'a-1', 'c-1']);
assert.strictEqual(singleSourceOrder[0], singleSourceRecords[1], 'record object identity is preserved');
assert.deepEqual(singleSourceRecords.map(({ uid }) => uid), ['a-1', 'b-1', 'c-1'], 'input order is not mutated');

const blockRecords = [
  sourceRecord('a-1', sourceA), sourceRecord('a-2', sourceA),
  sourceRecord('b-1', sourceB), sourceRecord('b-2', sourceB), sourceRecord('b-3', sourceB)
];
assert.deepEqual(
  moveLinearSourceGroup(blockRecords, 1, -1).map(({ uid }) => uid),
  ['b-1', 'b-2', 'b-3', 'a-1', 'a-2'],
  'a multi-record source moves as one block without changing its internal order'
);

for (const [index, direction] of [
  [0, -1], [2, 1], [-1, 1], [3, -1], [1.5, 1], ['1', -1], [1, 0], [1, 2], [1, '-1']
]) {
  assert.deepEqual(
    moveLinearSourceGroup(singleSourceRecords, index, direction),
    singleSourceRecords,
    `invalid source move ${index}:${direction} is a no-op`
  );
}

const duplicateNameA = file('duplicate.gb');
const duplicateNameB = file('duplicate.gb');
const duplicateNameRecords = [sourceRecord('same-a', duplicateNameA), sourceRecord('same-b', duplicateNameB)];
assert.deepEqual(
  moveLinearSourceGroup(duplicateNameRecords, 1, -1).map(({ uid }) => uid),
  ['same-b', 'same-a'],
  'distinct uploads with the same filename remain separate sources'
);

const removalRecords = [
  { ...sourceRecord('a-1', sourceA), region_record_id: 'A1', depth: [file('a.tsv')] },
  { ...sourceRecord('a-2', sourceA), region_record_id: 'A2', depth: [null] },
  sourceRecord('b-1', sourceB)
];
const clearPlan = planLinearSourceRemoval({
  sequences: removalRecords,
  sourceUid: 'a-1',
  intent: 'clear'
});
assert.deepEqual({
  allowed: clearPlan.allowed,
  sourceIndex: clearPlan.sourceIndex,
  insertionIndex: clearPlan.insertionIndex,
  recordCount: clearPlan.recordCount,
  removedUids: clearPlan.removedUids,
  retainedUids: clearPlan.retainedSequences.map(({ uid }) => uid)
}, {
  allowed: true,
  sourceIndex: 0,
  insertionIndex: 0,
  recordCount: 2,
  removedUids: ['a-1', 'a-2'],
  retainedUids: ['b-1']
});
assert.deepEqual(removalRecords.map(({ uid }) => uid), ['a-1', 'a-2', 'b-1'], 'planning does not mutate inputs');
assert.equal(planLinearSourceRemoval({
  sequences: removalRecords,
  sourceUid: 'a-1',
  intent: 'delete'
}).allowed, true);
assert.deepEqual(planLinearSourceRemoval({
  sequences: [sourceRecord('only', sourceA)],
  sourceUid: 'only',
  intent: 'delete'
}), { allowed: false, reason: 'sole-source', sourceIndex: 0 });
assert.deepEqual(planLinearSourceRemoval({
  sequences: removalRecords,
  sourceUid: 'missing',
  intent: 'clear'
}), { allowed: false, reason: 'missing-source' });

const pairedSource = groupLinearSourceRecords([{
  uid: 'paired', gb: null, gff: file('paired.gff3'), fasta: file('paired.fasta'),
  depth: [], definition: '', record_subtitle: '', file_definition: '', file_subtitle: '',
  region_record_id: '', region_start: null, region_end: null, region_reverse: false
}])[0];
assert.equal(linearSourceHasPrimaryInput(pairedSource), true, 'GFF3 and FASTA form one non-empty source');
assert.equal(isPristineLinearSource(pairedSource), false);
const pristineSource = groupLinearSourceRecords([{
  uid: 'blank', gb: null, gff: null, fasta: null, depth: [null, null],
  definition: '', record_subtitle: '', file_definition: '', file_subtitle: '',
  region_record_id: '', region_start: null, region_end: null, region_reverse: false
}])[0];
assert.equal(linearSourceHasPrimaryInput(pristineSource), false);
assert.equal(isPristineLinearSource(pristineSource), true);
pristineSource.sequence.file_subtitle = 'draft';
assert.equal(isPristineLinearSource(pristineSource), false, 'a configured blank slot is not pristine');

const sharedDescriptor = {
  kind: 'web-file', name: 'session.gb', type: 'text/plain', encoding: 'base64',
  data: '', size: 0, lastModified: 0
};
const sessionTable = adoptCurrentSessionResources({ source: sharedDescriptor });
const sessionRecords = [
  sourceRecord('session-1', inputFiles.registerFile(createSessionResourceFileView(sessionTable, 'source'))),
  sourceRecord('session-2', inputFiles.registerFile(createSessionResourceFileView(sessionTable, 'source'))),
  sourceRecord('other', sourceC)
];
assert.deepEqual(groupLinearSourceRecords(sessionRecords).map(({ records }) => records.length), [2, 1]);
assert.deepEqual(
  moveLinearSourceGroup(sessionRecords, 0, 1).map(({ uid }) => uid),
  ['other', 'session-1', 'session-2'],
  'records backed by one Session descriptor move as one source'
);

const sharedDepth = file('depth.tsv');
const sameNamedDepth = file('depth.tsv');
const depthSourceRecords = [
  { ...sourceRecord('depth-1', sourceA), depth: [] },
  { ...sourceRecord('depth-2', sourceA), depth: [] }
];
const depthSource = groupLinearSourceRecords(depthSourceRecords)[0];
assert.deepEqual(linearSourceDepthStatus(depthSource, 0), {
  state: 'empty', file: null, selectedCount: 0, recordCount: 2
});
depthSourceRecords[0].depth = [sharedDepth];
assert.deepEqual(linearSourceDepthStatus(depthSource, 0), {
  state: 'mixed', file: null, selectedCount: 1, recordCount: 2
});
depthSourceRecords[1].depth = [sharedDepth];
assert.deepEqual(linearSourceDepthStatus(depthSource, 0), {
  state: 'common', file: sharedDepth, selectedCount: 2, recordCount: 2
});
depthSourceRecords[1].depth = [sameNamedDepth];
assert.deepEqual(linearSourceDepthStatus(depthSource, 0), {
  state: 'mixed', file: null, selectedCount: 2, recordCount: 2
}, 'same-named independent Depth uploads remain mixed');

const depthDescriptor = {
  kind: 'web-file', name: 'saved-depth.tsv', type: 'text/tab-separated-values', encoding: 'base64',
  data: '', size: 0, lastModified: 0
};
const depthSessionTable = adoptCurrentSessionResources({ depth: depthDescriptor });
depthSourceRecords[0].depth = [inputFiles.registerFile(createSessionResourceFileView(depthSessionTable, 'depth'))];
depthSourceRecords[1].depth = [inputFiles.registerFile(createSessionResourceFileView(depthSessionTable, 'depth'))];
assert.equal(linearSourceDepthStatus(depthSource, 0).state, 'common',
  'views backed by one Session descriptor are one common assignment');
const distinctDepthSessionTable = adoptCurrentSessionResources({
  first: depthDescriptor,
  second: { ...depthDescriptor }
});
depthSourceRecords[0].depth = [inputFiles.registerFile(createSessionResourceFileView(distinctDepthSessionTable, 'first'))];
depthSourceRecords[1].depth = [inputFiles.registerFile(createSessionResourceFileView(distinctDepthSessionTable, 'second'))];
assert.equal(linearSourceDepthStatus(depthSource, 0).state, 'mixed',
  'same-named independent Session resources remain mixed');

const interleavedRecords = [
  sourceRecord('a-1', sourceA), sourceRecord('b-1', sourceB),
  sourceRecord('a-2', sourceA), sourceRecord('b-2', sourceB), sourceRecord('c-1', sourceC)
];
assert.deepEqual(
  moveLinearSourceGroup(interleavedRecords, 1, -1).map(({ uid }) => uid),
  ['b-1', 'b-2', 'a-1', 'a-2', 'c-1'],
  'an explicit move makes legacy interleaved source records contiguous'
);

const files = [file('same.gb'), file('same.gb')];
const sequences = Array.from({ length: 8 }, (_, index) => ({
  uid: `record-${index}`, gb: files[index < 6 ? 0 : 1], gff: null, fasta: null
}));
assert.deepEqual(groupLinearSourceRecords(sequences).map((source) => source.records.length), [6, 2]);
const specs = sequences.flatMap((_, queryIndex) => sequences.map((_, subjectIndex) => ({ queryIndex, subjectIndex })));
const planFor = (records = sequences, jobs = specs, getEntry = async () => ({ fasta: '>duplicate\nACGT\n' })) => prepareLosatSourceBatches({
  sequences: records, specs: jobs, getEntry, getDigest: digestOf(getEntry), buildArgs: () => ['--task', 'blastn'], hashText, protein: false
});
// PD-OI-018 revision 4 (D-19): sources stay batched between files; a record
// searches itself only when that self comparison is requested, and a
// comparison within one source searches the source without the query record.
const plan = await planFor();
const byScope = (batches) => Object.fromEntries(['between-sources', 'within-source', 'self']
  .map((scope) => [scope, batches.filter((batch) => batch.scope === scope).length]));
assert.deepEqual(byScope(plan.batches), { 'between-sources': 2, 'within-source': 8, self: 8 },
  'eight records in two sources: two between-source jobs, one within-source and one self job per record, never 64');
assert.equal(plan.batches.reduce((sum, batch) => sum + batch.specs.length, 0), 64);
assert.equal(
  planLosatSourceJobs({ sequences, specs, buildArgs: () => ['--task', 'blastn'] }).jobs.length,
  plan.batches.length,
  'the job-count estimate and the execution share one plan'
);
for (const batch of plan.batches) {
  assert.equal(batch.query.ids.size, batch.query.indexes.length);
  assert.equal(batch.subject.ids.size, batch.subject.indexes.length);
  if (batch.scope === 'self') {
    assert.deepEqual(batch.query.indexes, batch.subject.indexes, 'a requested self search uses the record alone');
    assert.equal(batch.query.indexes.length, 1);
  } else {
    assert.deepEqual(batch.query.indexes.filter((index) => batch.subject.indexes.includes(index)), [],
      'no unrequested self search');
  }
  if (batch.scope === 'within-source') {
    const [queryIndex] = batch.query.indexes;
    const source = queryIndex < 6 ? [0, 1, 2, 3, 4, 5] : [6, 7];
    assert.deepEqual([...batch.subject.indexes].sort(), source.filter((index) => index !== queryIndex),
      'the within-source database is the source without the query record');
  }
  assert.equal(batch.searchContext === null, batch.query.indexes.length === 1 && batch.subject.indexes.length === 1);
  const jobs = batch.specs.map((spec) => ({ ...spec, cacheKey: `${spec.queryIndex}:${spec.subjectIndex}` }));
  const text = [...batch.query.ids.keys()].flatMap((query) => [...batch.subject.ids.keys()].map(
    (subject) => [query, subject, 100, 4, 0, 0, 1, 4, 1, 4, '1e-20', 50].join('\t')
  )).join('\n');
  const split = splitLosatSourceResult(text, batch, jobs);
  assert.equal(split.length, jobs.length);
  split.forEach((result) => assert.match(result.text, /^duplicate\tduplicate\t100\t/));
  assert.equal(splitLosatSourceResult(text, batch, jobs.slice(0, 1)).length, 1, 'unselected pairs stay omitted');
  assert.throws(() => splitLosatSourceResult('unknown\tunknown', batch, jobs), /unrecognized/);
}
const adjacentSpecs = specs.filter((spec) => spec.queryIndex < 6 && spec.subjectIndex >= 6);
assert.equal((await planFor(sequences, adjacentSpecs)).batches.length, 1);
const intoSecondSource = (batches) => batches.find((batch) => batch.scope === 'between-sources' && batch.subject.indexes.includes(7));
const changed = await planFor(sequences, specs, async (index) => ({ fasta: `>duplicate\n${index === 7 ? 'ACGA' : 'ACGT'}\n` }));
assert.notEqual(intoSecondSource(plan.batches).searchContext, intoSecondSource(changed.batches).searchContext, 'database contents change cache scope');
const reordered = [...sequences].reverse();
const replay = await planFor(reordered);
assert.deepEqual(replay.batches.map((batch) => batch.searchContext).sort(), plan.batches.map((batch) => batch.searchContext).sort(), 'record ordering does not change the searched source sets');

const codePlan = await prepareLosatSourceBatches({
  sequences, specs: adjacentSpecs, getEntry: async () => ({ fasta: '>duplicate\nACGT\n' }),
  getDigest: async () => hashText('>duplicate\nACGT\n'),
  buildArgs: (query) => ['--query-gencode', query === 0 ? '4' : '11', '--db-gencode', '11'],
  hashText, protein: false
});
assert.equal(codePlan.batches.length, 2, 'conflicting explicit translation tables require compatible source batches');
assert.deepEqual(codePlan.batches.map((batch) => batch.query.indexes.length).sort(), [1, 5]);

const noSelf = await prepareLosatSourceBatches({
  sequences, specs: specs.filter(({ queryIndex, subjectIndex }) => queryIndex !== subjectIndex),
  getEntry: async (index) => ({ fasta: `>protein-${index}\nMKK\n` }),
  getDigest: async (index) => hashText(`>protein-${index}\nMKK\n`),
  buildArgs: () => ['--max-target-seqs', '5'], hashText, protein: true
});
for (const batch of noSelf.batches) {
  assert(!batch.query.indexes.some((index) => batch.subject.indexes.includes(index)),
    'no mode submits an unrequested within-record search, including multi-record sources');
}
assert.equal(noSelf.batches.length, 2 + 8);
assert.equal(noSelf.batches.reduce((sum, batch) => sum + batch.specs.length, 0), 56);

// OV-1598: a single-record side searches its entry's FASTA as it is, so its digest is the entry's own (read from
// `getDigest`, the owner run-analysis.js holds per File) and is not computed again. A side whose text differs from
// its entry (renamed records, an added line end) hashes its own text.
test('a single-record LOSAT side reads its entry digest instead of hashing the FASTA again', async () => {
  const fileA = file('a.gbk');
  const cards = ['a0', 'a1', 'b0'].map((uid) => ({ uid, gb: uid === 'b0' ? file('b.gbk') : fileA, gff: null, fasta: null }));
  const fastaFor = (index) => (index === 2 ? '>b0\nACGT' : `>${cards[index].uid}\nACGT\n`);
  const hashed = [];
  const digests = [];
  const plan = await prepareLosatSourceBatches({
    sequences: cards,
    specs: [{ queryIndex: 0, subjectIndex: 0 }, { queryIndex: 0, subjectIndex: 2 }, { queryIndex: 2, subjectIndex: 0 }],
    getEntry: async (index) => ({ fasta: fastaFor(index) }),
    getDigest: async (index) => { digests.push(index); return hashText(fastaFor(index)); },
    buildArgs: () => ['--task', 'blastn'],
    hashText: async (text) => { hashed.push(text); return hashText(text); },
    protein: false
  });
  assert.deepEqual(digests, [0], 'only the record whose side text is its entry text reads the entry digest');
  assert.ok(!hashed.includes(fastaFor(0)), 'the entry FASTA is not hashed again');
  assert.ok(hashed.includes('>b0\nACGT\n'), 'a side with an added line end hashes its own text');
  const self = plan.batches.find((batch) => batch.scope === 'self');
  assert.equal(self.query.hash, await hashText(fastaFor(0)));
  const multi = plan.batches.find((batch) => batch.subject.indexes.length > 1);
  assert.equal(multi.subject.hash, await hashText(multi.subject.fasta), 'a two-record side hashes its joined FASTA');
});

// Record selection D-07: an OFF record is never a query, and its source's
// database keeps it, so turning a record ON or OFF changes no other E-value.
test('the LOSAT plan over the drawn records keeps OFF records in their source database', async () => {
  const fileA = file('a.gbk');
  const fileB = file('b.gbk');
  const cards = [
    { uid: 'a0', gb: fileA }, { uid: 'a1', gb: fileA }, { uid: 'a2', gb: fileA },
    { uid: 'b0', gb: fileB }, { uid: 'b1', gb: fileB }
  ].map((card) => ({ gff: null, fasta: null, ...card }));
  const recordsOff = new Set(['a1', 'b1']);
  const drawn = cards.filter((card) => !recordsOff.has(card.uid));
  const search = losatSearchSequences(drawn, cards);
  assert.deepEqual(search.sequences.map((card) => card.uid), ['a0', 'a2', 'b0', 'a1', 'b1']);
  assert.equal(search.drawnCount, 3);
  // Pairs only between drawn records: a0 -> b0, a2 -> b0, and a0 -> a2 within file A.
  const drawnSpecs = [{ queryIndex: 0, subjectIndex: 2 }, { queryIndex: 1, subjectIndex: 2 }, { queryIndex: 0, subjectIndex: 1 }];
  const getEntry = async (index) => ({ fasta: `>${search.sequences[index].uid}\nACGT\n` });
  const uids = (indexes) => indexes.map((index) => search.sequences[index].uid).sort();
  const planned = await prepareLosatSourceBatches({
    ...search, specs: drawnSpecs, getEntry, getDigest: digestOf(getEntry), buildArgs: () => ['--task', 'blastn'], hashText, protein: false
  });
  const between = planned.batches.find((batch) => batch.scope === 'between-sources');
  assert.deepEqual(uids(between.query.indexes), ['a0', 'a2'], 'an OFF record is not a query');
  assert.deepEqual(uids(between.subject.indexes), ['b0', 'b1'], 'the database is the whole subject file');
  const within = planned.batches.find((batch) => batch.scope === 'within-source');
  assert.deepEqual(uids(within.subject.indexes), ['a1', 'a2'], 'a within-file database keeps the OFF record');
  assert.equal(planLosatSourceJobs({ ...search, specs: drawnSpecs, buildArgs: () => [] }).jobs.length, planned.batches.length);

  // Turning b1 ON keeps the database of file B, so a0 -> b0 keeps its search scope.
  const allOn = losatSearchSequences(cards, cards);
  const allOnPlan = await prepareLosatSourceBatches({
    ...allOn,
    specs: [{ queryIndex: 0, subjectIndex: 3 }, { queryIndex: 2, subjectIndex: 3 }],
    getEntry: async (index) => ({ fasta: `>${allOn.sequences[index].uid}\nACGT\n` }),
    getDigest: async (index) => hashText(`>${allOn.sequences[index].uid}\nACGT\n`),
    buildArgs: () => ['--task', 'blastn'], hashText, protein: false
  });
  const allOnBetween = allOnPlan.batches.find((batch) => batch.scope === 'between-sources');
  assert.equal(allOnBetween.subject.hash, between.subject.hash);
});

// CO-04 (D-19, PD-OI-018 revision 4): two records packaged in one source file
// search the same database as two single-record files, so an unrequested
// self hit cannot fill max-target-seqs.
test('one-file records never search themselves without a request', async () => {
  const planFiles = (packaged) => prepareLosatSourceBatches({
    sequences: [0, 1].map((index) => ({
      uid: `record-${index}`, gb: packaged ? packaged : file(`record-${index}.gbk`), gff: null, fasta: null
    })),
    specs: [{ queryIndex: 0, subjectIndex: 1 }],
    getEntry: async (index) => ({ fasta: `>protein-${index}\nMKK\n` }),
    getDigest: async (index) => hashText(`>protein-${index}\nMKK\n`),
    buildArgs: () => ['--max-target-seqs', '1'],
    hashText,
    protein: true
  });
  const oneFile = await planFiles(file('two-records.gbk'));
  const twoFiles = await planFiles(null);
  for (const batch of oneFile.batches) {
    assert.deepEqual(batch.query.indexes.filter((index) => batch.subject.indexes.includes(index)), [],
      'an unrequested within-record search fills max-target-seqs with self hits');
  }
  const identity = ({ batches }) => batches.map(({ query, subject, searchContext }) => [query.hash, subject.hash, searchContext]);
  assert.deepEqual(identity(oneFile), identity(twoFiles), 'file packaging does not change a two-record search');
});

// Shared vectors (design 3.6, D15): tests/test_losat_job_plan.py runs the same
// cases against gbdraw.comparisons.losat_jobs, so the CLI and the Web search
// the same databases under the same raw keys.
const sharedVectors = (name) => JSON.parse(
  readFileSync(new URL(`../fixtures/${name}`, import.meta.url), 'utf8')
).cases;

test('shared LOSAT FASTA extraction vectors', async () => {
  for (const vector of sharedVectors('losat_fasta_extraction_cases.json')) {
    const entry = await extractLosatFastaFast({
      text: vector.text, fmt: vector.fmt, regionSpec: vector.regionSpec,
      recordSelector: vector.recordSelector
    });
    assert.deepEqual({
      recordId: entry.recordId,
      canonicalLength: entry.canonicalLength,
      fasta: entry.fasta,
      hash: await hashText(entry.fasta)
    }, vector.expected, vector.name);
  }
});

test('shared LOSAT job plan and nucleotide raw-key vectors', async () => {
  for (const vector of sharedVectors('losat_job_plan_cases.json')) {
    const files = new Map();
    const records = vector.records.map((record) => {
      if (!files.has(record.source)) files.set(record.source, file(record.source));
      return { uid: record.uid, gb: files.get(record.source), gff: null, fasta: null, losat_gencode: record.gencode };
    });
    // The args Generate builds for one record pair (`buildLosatArgs` of run-analysis.js), without a genome pair.
    const buildArgs = (query, subject) => nucleotideLosatArgs({
      program: vector.program, task: vector.task, searchSpace: null,
      gencodes: () => [vector.records[query].gencode, vector.records[subject].gencode]
    });
    const specs = vector.specs.map(([queryIndex, subjectIndex]) => ({ queryIndex, subjectIndex }));
    const plan = await prepareLosatSourceBatches({
      sequences: records, specs, getEntry: async (index) => ({ fasta: vector.records[index].fasta }),
      getDigest: async (index) => hashText(vector.records[index].fasta), buildArgs, hashText, protein: false
    });
    assert.equal(planLosatSourceJobs({ sequences: records, specs, buildArgs }).jobs.length, plan.batches.length);
    assert.deepEqual(plan.batches.map((batch) => ({
      scope: batch.scope,
      args: batch.args,
      specs: batch.specs.map((spec) => [spec.queryIndex, spec.subjectIndex]),
      queryIndexes: batch.query.indexes,
      subjectIndexes: batch.subject.indexes,
      queryIds: Object.fromEntries([...batch.query.ids].map(([id, value]) => [id, [value.index, value.originalId]])),
      subjectIds: Object.fromEntries([...batch.subject.ids].map(([id, value]) => [id, [value.index, value.originalId]])),
      queryHash: batch.query.hash,
      subjectHash: batch.subject.hash,
      searchContext: batch.searchContext,
      searchable: batch.searchable
    })), vector.expected.jobs, vector.name);
    const rawKeys = [];
    for (const spec of specs) {
      const payload = buildLosatCachePayload({
        identityKind: 'nucleotide',
        program: vector.program,
        outfmt: '6',
        args: buildArgs(spec.queryIndex, spec.subjectIndex),
        queryCanonicalHash: await hashText(vector.records[spec.queryIndex].fasta),
        subjectCanonicalHash: await hashText(vector.records[spec.subjectIndex].fasta),
        searchContext: plan.bySpec.get(spec).searchContext
      });
      rawKeys.push([spec.queryIndex, spec.subjectIndex, await hashText(JSON.stringify(payload))]);
    }
    assert.deepEqual(rawKeys, vector.expected.rawKeys, vector.name);
  }
});

// LOSATP uses the same plan (design D7): tests/test_losat_job_plan.py runs
// these cases against losatp_job_specs and prepare_losat_batches(protein=True).
// The Web protein raw key comes from build_protein_losat_cache_key through
// the Python helper, so equal args and searchContext give equal keys.
test('shared LOSATP specs and job plan vectors', async () => {
  const cases = JSON.parse(
    readFileSync(new URL('../fixtures/losat_job_plan_cases.json', import.meta.url), 'utf8')
  ).proteinCases;
  assert.ok(cases.length > 0);
  for (const vector of cases) {
    const files = new Map();
    const sequences = vector.records.map((record) => {
      if (!files.has(record.source)) files.set(record.source, file(record.source));
      return { uid: record.uid, gb: files.get(record.source), gff: null, fasta: null };
    });
    const resolution = {
      mode: 'adjacent',
      edges: vector.edges.map(([queryIndex, subjectIndex], ordinal) => ({
        source: 'losat', queryIndex, subjectIndex, ordinal, edgeKey: `edge-${ordinal}`,
        queryUid: vector.records[queryIndex].uid, subjectUid: vector.records[subjectIndex].uid
      }))
    };
    const specs = buildLosatJobSpecs({
      resolution,
      recordCount: vector.records.length,
      recordUids: vector.records.map((record) => record.uid),
      program: 'blastp',
      blastpMode: vector.mode,
      collinearInferOrthogroups: vector.inferOrthogroups,
      collinearSearchScope: vector.searchScope
    });
    assert.deepEqual(specs.map((spec) => [spec.queryIndex, spec.subjectIndex]), vector.specs, vector.name);
    const plan = await prepareLosatSourceBatches({
      sequences, specs, getEntry: async (index) => ({ fasta: vector.records[index].fasta }),
      getDigest: async (index) => hashText(vector.records[index].fasta), buildArgs: () => vector.args, hashText, protein: true
    });
    assert.deepEqual(plan.batches.map((batch) => ({
      scope: batch.scope,
      args: batch.args,
      specs: batch.specs.map((spec) => [spec.queryIndex, spec.subjectIndex]),
      queryIndexes: batch.query.indexes,
      subjectIndexes: batch.subject.indexes,
      queryIds: Object.fromEntries([...batch.query.ids].map(([id, value]) => [id, [value.index, value.originalId]])),
      subjectIds: Object.fromEntries([...batch.subject.ids].map(([id, value]) => [id, [value.index, value.originalId]])),
      queryHash: batch.query.hash,
      subjectHash: batch.subject.hash,
      searchContext: batch.searchContext,
      searchable: batch.searchable
    })), vector.expected.jobs, vector.name);
  }
});

// GX-09: the LOSAT FASTA extraction reports a region outside the record as
// crop_genbank.py does (CI-01): REGION_INVALID with the region correction.
test('LOSAT FASTA extraction reports a region beyond the record as a region failure (GX-09)', async () => {
  const text = '>TINY.1\n' + 'ACGT'.repeat(15) + '\n';
  const extract = (regionSpec) => extractLosatFastaFast({ text, fmt: 'fasta', regionSpec, recordSelector: null });
  await assert.rejects(extract('1000-2000'), (error) => {
    assert.equal(error.code, 'REGION_INVALID');
    assert.deepEqual(error.context, { field: 'region', reason: 'RECORD_BOUNDS' });
    return true;
  });
  // A partial overrun is clamped, and a reversed span is a reverse-complemented crop
  // (parseRegionText, like parse_region_spec), so the ORDER reason has no input here.
  assert.equal((await extract('31-2000')).canonicalLength, 30);
  assert.equal((await extract('2-1')).fasta, '>TINY.1\nGT\n');
});

// A GenBank record without a sequence (a CON record: CONTIG, no ORIGIN) is a
// record of the file for the record list and for Python's SeqIO.parse, so the
// LOSAT reader counts it too: `#2` names that record, whose FASTA only the
// Python reader can report (it reports UndefinedSequenceError), and `#3` names
// the third record, not a fourth.
test('a GenBank record without a sequence keeps its place in the LOSAT records', async () => {
  const record = (id, body) => [
    `LOCUS       ${id}                       8 bp    DNA     linear   UNK 01-JAN-1980`,
    `ACCESSION   ${id}`, 'FEATURES             Location/Qualifiers', ...body, '//'
  ].join('\n');
  const withSequence = (id, bases) => record(id, ['ORIGIN', `        1 ${bases}`]);
  const contigOnly = (id) => record(id, [`CONTIG      join(${id}0:1..8)`]);
  const extract = (text, recordSelector) => extractLosatFastaFast({ text, fmt: 'genbank', regionSpec: null, recordSelector });
  const middle = [withSequence('AAA', 'acgtacgt'), contigOnly('BBB'), withSequence('CCC', 'ttttgggg')].join('\n') + '\n';
  assert.equal((await extract(middle, '')).recordId, 'AAA');
  await assert.rejects(extract(middle, '#2'), /no sequence/);
  await assert.rejects(extract(middle, 'BBB'), /no sequence/);
  assert.equal((await extract(middle, '#3')).fasta, '>CCC\nTTTTGGGG\n');
  // The automatic record is the first record, as the drawing reads it.
  const first = [contigOnly('BBB'), withSequence('AAA', 'acgtacgt')].join('\n') + '\n';
  await assert.rejects(extract(first, ''), /no sequence/);
  assert.equal((await extract(first, '#2')).recordId, 'AAA');
});

// OV-436: two chromosomes of one File were named "Vibrio parahaemolyticus" and
// "NC_004605.1" because only the first drew a Definition.
test('every record of a File is named by the File definition plus its accession', () => {
  const vp = file('vp.gb');
  const va = file('va.gb');
  const card = (uid, gb, definition, region_record_id) => ({ ...sourceRecord(uid, gb), definition, region_record_id });
  const sequences = [
    card('r1', vp, '<i>Vibrio parahaemolyticus</i>', 'NC_004603.1'),
    card('r2', vp, '', 'NC_004605.1'),
    card('r3', va, '', 'NC_022349.1'),
    card('r4', file('single.gb'), '', '')
  ];
  const before = JSON.stringify(sequences);
  const labels = [...linearRecordNames(sequences).values()].map(({ label }) => label);
  assert.deepEqual(labels, [
    'Vibrio parahaemolyticus · NC_004603.1',
    'Vibrio parahaemolyticus · NC_004605.1',
    'NC_022349.1',
    'single.gb'
  ]);
  // The cards (and the diagram headings they draw) are not changed.
  assert.equal(JSON.stringify(sequences), before);
});

// The folded File cards: view state keyed by source uid.
const foldFile = (name) => file(name);
const foldSources = () => {
  const a = foldFile('a.gbk');
  const b = foldFile('b.gbk');
  return groupLinearSourceRecords([
    sourceRecord('a1', a), sourceRecord('a2', a), sourceRecord('b1', b), sourceRecord('empty', null)
  ]);
};

test('a press on a File card header folds and unfolds that card alone', () => {
  const folded = new Set();
  setLinearSourceFolded(folded, 'a1', true);
  setLinearSourceFolded(folded, 'a1', true);
  assert.deepEqual([...folded], ['a1']);
  setLinearSourceFolded(folded, 'b1', true);
  setLinearSourceFolded(folded, 'a1', false);
  setLinearSourceFolded(folded, 'a1', false);
  assert.deepEqual([...folded], ['b1']);
});

test('a loaded Session folds the File cards that hold their input and leaves the others open', () => {
  const folded = new Set(['stale', 'empty']);
  foldLoadedLinearSources(folded, foldSources(), 'gb');
  assert.deepEqual([...folded], ['a1', 'b1'], 'by the uid of each File, not of each record');
  foldLoadedLinearSources(folded, groupLinearSourceRecords([sourceRecord('only', null)]), 'gb');
  assert.deepEqual([...folded], [], 'a Session without files folds nothing');
});

test('a card holds its input when every slot its input type shows has a file', () => {
  const card = (files) => groupLinearSourceRecords([{ uid: 'c', gb: null, gff: null, fasta: null, ...files }])[0];
  const gff = foldFile('g.gff3');
  const fasta = foldFile('g.fasta');
  assert.deepEqual(linearCardFileFields('gb'), ['gb']);
  assert.deepEqual(linearCardFileFields('gff'), ['gff', 'fasta']);
  assert.equal(linearCardInputComplete(card({ gb: foldFile('a.gbk') }), 'gb'), true);
  assert.equal(linearCardInputComplete(card({ gff, fasta }), 'gff'), true);
  assert.equal(linearCardInputComplete(card({ gff }), 'gff'), false, 'GFF3 without FASTA');
  assert.equal(linearCardInputComplete(card({ fasta }), 'gff'), false, 'FASTA without GFF3');
  assert.equal(linearCardInputComplete(card({ gb: foldFile('a.gbk') }), 'gff'), false, 'only a file of the other input type');
  assert.equal(linearCardInputComplete(card({ gff, fasta }), 'gb'), false);
  // The load rule is that answer: an incomplete pair stays open like a card without a file.
  const folded = new Set();
  foldLoadedLinearSources(folded, [card({ gff, fasta }), { ...card({ gff }), uid: 'half' }, { ...card({ gb: foldFile('a.gbk') }), uid: 'other' }], 'gff');
  assert.deepEqual([...folded], ['c']);
});

test('a rolled-back load returns the fold state it replaced', () => {
  const folded = new Set(['a1']);
  const before = [...folded];
  foldLoadedLinearSources(folded, groupLinearSourceRecords([sourceRecord('x', foldFile('x.gbk'))]), 'gb');
  restoreFoldedLinearSources(folded, before);
  assert.deepEqual([...folded], ['a1']);
});

test('the fold state follows a File through Move up and Move down', () => {
  const sequences = foldSources().flatMap(({ records }) => records.map(({ sequence }) => sequence));
  const folded = new Set();
  setLinearSourceFolded(folded, 'b1', true);
  const moved = groupLinearSourceRecords(moveLinearSourceGroup(sequences, 1, -1));
  assert.deepEqual(moved.map(({ uid }) => uid), ['b1', 'a1', 'empty']);
  assert.deepEqual(moved.map(({ uid }) => folded.has(uid)), [true, false, false]);
  assert.equal(moved[0].sequence.gb.name, 'b.gbk');
});

// The input-type table (services/input-types.js): each type's record source files are the generated
// RECORD_SOURCE_FILES of its source kind, the one builder refuses a source that misses a file, and a card holds an
// input only through a slot its type reads.
test('input types read the record source files Python names, and a card holds a file only in its type\'s slots', () => {
  for (const inputType of Object.keys(INPUT_TYPES)) {
    const { sourceKind } = INPUT_TYPES[inputType];
    assert.deepEqual(inputSourceSlotFiles(inputType).map(({ resource, field }) => ({ resource, field })),
      RECORD_SOURCE_FILES[sourceKind], inputType);
  }
  assert.deepEqual(recordSourceOf('gffFasta', ({ slot }) => `${slot}-id`),
    { kind: 'gffFasta', gffResourceId: 'gff-id', fastaResourceId: 'fasta-id' });
  assert.equal(recordSourceOf('gffFasta', ({ slot }) => (slot === 'gff' ? 'gff-id' : null)), null);
  assert.deepEqual(recordSourceOf('fasta', () => 'only'), { kind: 'fasta', fastaResourceId: 'only' });
  assert.equal(isRecordSource({ kind: 'gffFasta', gffResourceId: 'a' }), false);
  assert.equal(isRecordSource({ kind: 'fasta', fastaResourceId: 'a', resourceId: 'b' }), false);
  assert.deepEqual(recordSourceResourceIds({ kind: 'gffFasta', gffResourceId: 'a', fastaResourceId: 'b' }), ['a', 'b']);
  assert.deepEqual(recordSourceSlotIds({ kind: 'fasta', fastaResourceId: 'f', resourceId: 'stale' }),
    { gb: '', gff: '', fasta: 'f' });
  assert.deepEqual(storedRecordSourceResourceIds({ kind: 'unknown', resourceId: 'x' }), []);
  assert.deepEqual(storedRecordSourceResourceIds({ kind: 'fasta', fastaResourceId: ' f ' }), ['f']);
  assert.deepEqual(storedRecordSourceResourceIds({ kind: 'gffFasta', gffResourceId: 'g', fastaResourceId: 'f' }), ['g', 'f']);
  assert.deepEqual(storedRecordSourceResourceIds({ kind: 'gffFasta', gffResourceId: 'g' }), ['g', '']);
  assert.equal(recordSourceFilesOf('fasta'), RECORD_SOURCE_FILES.fasta);
  assert.deepEqual(recordSourceFilesOf('unknown'), []);
  // A card that held GFF3 + FASTA and was switched to FASTA holds nothing once its FASTA is removed.
  assert.equal(inputSourceFiles('fasta', { gff: 'annotation.gff3', fasta: null }).holdsFile, false);
  assert.equal(inputSourceFiles('gff', { gff: 'annotation.gff3', fasta: null }).holdsFile, true);
  assert.equal(inputSourceFiles('gb', { gb: null, fasta: 'genome.fasta' }).holdsFile, false);
});

// PD-OI-018 revision 6: one E-value scale per genome. The Web's args for each unit pair of a figure (a crop, two
// windows, a genome split over two files, and a whole single-record pair, which passes none; LOSATN, TLOSATX and LOSATP)
// are built from the scales Python plans for that figure and equal the command line's raw args
// (tests/test_losat_job_plan.py runs the same figures through the command line and the planner).
test('the Web passes the E-value scale Python plans, as the command line does', async () => {
  const { readGenomeScales } = await import('../../gbdraw/web/js/services/linear-units.js');
  const { searchSpaceCases } = JSON.parse(readFileSync(new URL('../fixtures/losat_job_plan_cases.json', import.meta.url), 'utf8'));
  assert.ok(searchSpaceCases.filter(({ program }) => program === 'blastp').length >= 3);
  for (const vector of searchSpaceCases) {
    const { scaleOf } = readGenomeScales(vector.reply, vector.reply.recordKeys);
    const args = vector.expected.args.map(([query, subject]) => [query, subject, vector.program === 'blastp'
      ? proteinLosatArgs({ pairwise: vector.mode === 'pairwise', candidateLimit: vector.candidateLimit, databaseSize: scaleOf(query, subject) })
      : nucleotideLosatArgs({
        program: vector.program, task: vector.task, gencodes: () => [null, null], searchSpace: scaleOf(query, subject)
      })]);
    assert.deepEqual(args, vector.expected.args, vector.name);
  }
});

// The reply is read against the units Generate searches: another unit list, a scale that is not a positive integer,
// or a unit the reply does not hold is a helper failure, never a search without its scale.
test('a genome scale reply that does not fit the searched units is refused', async () => {
  const { readGenomeScales } = await import('../../gbdraw/web/js/services/linear-units.js');
  const params = { first: 'a.gb', second: 'b.gb', space: '9', limit: '8' };
  const warning = { queryUnit: 0, subjectUnit: 1, text: 'Genomes {first} and {second}', params };
  const good = { recordKeys: ['a', 'b'], scales: [[null, '7'], ['5', null]], warnings: [warning] };
  assert.equal(readGenomeScales(good, ['a', 'b']).scaleOf(0, 1), '7');
  assert.equal(readGenomeScales(good, ['a', 'b']).scaleOf(1, 1), null);
  assert.deepEqual(readGenomeScales(good, ['a', 'b']).warningOf(0, 1), warning);
  assert.equal(readGenomeScales(good, ['a', 'b']).warningOf(1, 0), null);
  const refused = (result, keys = ['a', 'b']) => assert.throws(() => readGenomeScales(result, keys), { code: 'RESULT_INVALID' });
  refused({ ...good, warnings: [{ ...warning, subjectUnit: 2 }] });
  refused({ ...good, warnings: [{ ...warning, text: '' }] });
  refused({ ...good, warnings: [{ ...warning, params: { ...params, space: 9 } }] });
  refused({ ...good, warnings: [{ queryUnit: 0, subjectUnit: 1, message: 'searched per file' }] });
  refused({ recordKeys: good.recordKeys, scales: good.scales });
  refused(good, ['b', 'a']);
  refused(good, ['a']);
  refused({ ...good, scales: [[null, '0'], ['5', null]] });
  refused({ ...good, scales: [[null, '07'], ['5', null]] });
  refused({ ...good, scales: [[null, '1.5'], ['5', null]] });
  refused({ ...good, scales: [[null, 7], ['5', null]] });
  refused({ ...good, scales: [[null, '7']] });
  refused({ recordKeys: ['a', 'b'] });
  refused(null);
  assert.throws(() => readGenomeScales(good, ['a', 'b']).scaleOf(2, 0), { code: 'RESULT_INVALID' });
});

// PD-OI-018 revision 6 (OV-3319): a Circular Generate reads one scale per ring it searches; a reply for other rings, a
// scale that is not a positive integer, or a warning of a ring it does not search is a helper failure.
test('a ring scale reply that does not fit the searched rings is refused', async () => {
  const { readRingScales } = await import('../../gbdraw/web/js/services/linear-units.js');
  const params = { first: 'a.gb', second: 'b.fna', space: '9', limit: '8' };
  const warning = { ring: 1, text: 'Genomes {first} and {second}', params };
  const good = { scales: ['12', null], warnings: [warning] };
  assert.deepEqual(readRingScales(good, 2), good);
  const refused = (result, count = 2) => assert.throws(() => readRingScales(result, count), { code: 'RESULT_INVALID' });
  refused(good, 1);
  refused(good, 3);
  refused({ ...good, scales: ['012', null] });
  refused({ ...good, scales: [12, null] });
  refused({ ...good, warnings: [{ ...warning, ring: 2 }] });
  refused({ ...good, warnings: [{ ...warning, params: { ...params, limit: 8 } }] });
  refused({ scales: good.scales });
  refused(null);
});

// Owner decision 2026-11-10: a genome pair past what LOSAT reads is searched per file with one warning. The banner shows
// the warnings of the pairs Generate searches, each once whatever the direction, translated from the English text
// Python sends with its parameters (the text the command line logs); a pair of the reply no search reads shows none.
test('the banner shows each per-file warning of the searched pairs once, translated', async () => {
  const { readGenomeScales, perFileWarningNotice } = await import('../../gbdraw/web/js/services/linear-units.js');
  const { GENOME_SCALE_WARNING_TEXT } = await import('../../gbdraw/web/js/contracts/card-windows.generated.js');
  const params = (first, second) => ({ first, second, space: '9,223,372,036,854,775,808', limit: '9,223,372,036,854,775,807' });
  const reply = {
    recordKeys: ['a', 'b', 'c'],
    scales: [[null, null, null], [null, null, null], [null, null, null]],
    warnings: [
      { queryUnit: 0, subjectUnit: 1, text: GENOME_SCALE_WARNING_TEXT, params: params('a.gb', 'b.gb') },
      { queryUnit: 1, subjectUnit: 0, text: GENOME_SCALE_WARNING_TEXT, params: params('a.gb', 'b.gb') },
      { queryUnit: 2, subjectUnit: 2, text: GENOME_SCALE_WARNING_TEXT, params: params('c.gb', 'c.gb') }
    ]
  };
  const { warningOf } = readGenomeScales(reply, reply.recordKeys);
  /** @type {{ text: string, tone: string }[]} */
  const shown = [];
  const notice = perFileWarningNotice((entry) => shown.push(entry));
  for (const [query, subject] of [[0, 1], [1, 0], [0, 2], [1, 2]]) notice(warningOf(query, subject));
  const english = 'The search space of genomes a.gb and b.gb (9,223,372,036,854,775,808) is past the largest LOSAT reads '
    + '(9,223,372,036,854,775,807), so their E-values are computed per file, not per genome.';
  assert.deepEqual(shown, [{ text: english, tone: 'warning' }]);

  // The text is a catalog key: a locale that has it shows its translation with the same parameters.
  const { setUiLocale, SOURCE_LOCALE } = await import('../../gbdraw/web/js/utils/i18n.js');
  await setUiLocale('ja', async () => ({ default: { [GENOME_SCALE_WARNING_TEXT]: '{first} / {second}: per file' } }));
  try {
    const translated = [];
    const translatedNotice = perFileWarningNotice((entry) => translated.push(entry.text));
    translatedNotice(warningOf(0, 1));
    translatedNotice(warningOf(2, 2));
    // Two warnings join with the locale's sentence separator (none in Japanese).
    assert.deepEqual(translated, ['a.gb / b.gb: per file', 'a.gb / b.gb: per filec.gb / c.gb: per file']);
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

// Two genomes of 120 Mb and 95 Mb: their search space is past Number.MAX_SAFE_INTEGER, and reaches LOSAT digit for digit.
test('a search space past what a number holds exactly reaches the args unchanged', async () => {
  const { readGenomeScales } = await import('../../gbdraw/web/js/services/linear-units.js');
  const space = String(120000000n * 95000000n);
  assert.ok(BigInt(space) > BigInt(Number.MAX_SAFE_INTEGER));
  const { scaleOf } = readGenomeScales({ recordKeys: ['a', 'b'], scales: [[null, space], [space, null]], warnings: [] }, ['a', 'b']);
  const args = nucleotideLosatArgs({ program: 'blastn', task: 'megablast', gencodes: () => [null, null], searchSpace: scaleOf(0, 1) });
  assert.deepEqual(args, ['--task', 'megablast', '--searchsp', '11400000000000000']);
});

// An OFF record has no genome: it joins the database of the drawn subject it is searched with, under that subject's
// args, so a genome pair's search space does not keep it out (record selection D-07, OV-3143).
test('an OFF record joins a subject database under the args of the drawn subject', () => {
  const a = file('a.gb');
  const b = file('b.gb');
  const sequences = [
    { uid: 'q', gb: a, gff: null, fasta: null }, { uid: 's1', gb: b, gff: null, fasta: null },
    { uid: 's2-off', gb: b, gff: null, fasta: null }
  ];
  const buildArgs = (query, subject, drawnSubject = subject) => nucleotideLosatArgs({
    program: 'blastn', task: 'megablast', gencodes: () => [null, null],
    searchSpace: drawnSubject === 1 ? 600 : null
  });
  const { jobs } = planLosatSourceJobs({ sequences, specs: [{ queryIndex: 0, subjectIndex: 1 }], buildArgs, drawnCount: 2 });
  assert.deepEqual(jobs.map(({ queryIndexes, subjectIndexes, args }) => ({ queryIndexes, subjectIndexes, args })), [
    { queryIndexes: [0], subjectIndexes: [1, 2], args: ['--task', 'megablast', '--searchsp', '600'] }
  ]);
});
