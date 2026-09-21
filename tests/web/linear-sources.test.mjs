import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import {
  groupLinearSourceRecords,
  isPristineLinearSource,
  linearSourceHasPrimaryInput,
  linearSourceDepthStatus,
  moveLinearSourceGroup,
  planLinearSourceRemoval,
  prepareLosatSourceBatches,
  splitLosatSourceResult
} from '../../gbdraw/web/js/app/linear-sources.js';
import {
  adoptCurrentSessionResources,
  createSessionResourceFileView
} from '../../gbdraw/web/js/services/session-resource-backing.js';

const hashText = async (value) => createHash('sha256').update(value).digest('hex');

const sourceRecord = (uid, gb) => ({ uid, gb, gff: null, fasta: null });
const sourceA = { name: 'a.gb' };
const sourceB = { name: 'b.gb' };
const sourceC = { name: 'c.gb' };
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

const duplicateNameA = { name: 'duplicate.gb' };
const duplicateNameB = { name: 'duplicate.gb' };
const duplicateNameRecords = [sourceRecord('same-a', duplicateNameA), sourceRecord('same-b', duplicateNameB)];
assert.deepEqual(
  moveLinearSourceGroup(duplicateNameRecords, 1, -1).map(({ uid }) => uid),
  ['same-b', 'same-a'],
  'distinct uploads with the same filename remain separate sources'
);

const removalRecords = [
  { ...sourceRecord('a-1', sourceA), region_record_id: 'A1', depth: [{ name: 'a.tsv' }] },
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
  uid: 'paired', gb: null, gff: { name: 'paired.gff3' }, fasta: { name: 'paired.fasta' },
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
  sourceRecord('session-1', createSessionResourceFileView(sessionTable, 'source')),
  sourceRecord('session-2', createSessionResourceFileView(sessionTable, 'source')),
  sourceRecord('other', sourceC)
];
assert.deepEqual(groupLinearSourceRecords(sessionRecords).map(({ records }) => records.length), [2, 1]);
assert.deepEqual(
  moveLinearSourceGroup(sessionRecords, 0, 1).map(({ uid }) => uid),
  ['other', 'session-1', 'session-2'],
  'records backed by one Session descriptor move as one source'
);

const sharedDepth = { name: 'depth.tsv' };
const sameNamedDepth = { name: 'depth.tsv' };
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
depthSourceRecords[0].depth = [createSessionResourceFileView(depthSessionTable, 'depth')];
depthSourceRecords[1].depth = [createSessionResourceFileView(depthSessionTable, 'depth')];
assert.equal(linearSourceDepthStatus(depthSource, 0).state, 'common',
  'views backed by one Session descriptor are one common assignment');
const distinctDepthSessionTable = adoptCurrentSessionResources({
  first: depthDescriptor,
  second: { ...depthDescriptor }
});
depthSourceRecords[0].depth = [createSessionResourceFileView(distinctDepthSessionTable, 'first')];
depthSourceRecords[1].depth = [createSessionResourceFileView(distinctDepthSessionTable, 'second')];
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

const files = [{ name: 'same.gb' }, { name: 'same.gb' }];
const sequences = Array.from({ length: 8 }, (_, index) => ({
  uid: `record-${index}`, gb: files[index < 6 ? 0 : 1], gff: null, fasta: null
}));
assert.deepEqual(groupLinearSourceRecords(sequences).map((source) => source.records.length), [6, 2]);
const specs = sequences.flatMap((_, queryIndex) => sequences.map((_, subjectIndex) => ({ queryIndex, subjectIndex })));
const planFor = (records = sequences, jobs = specs, getEntry = async () => ({ fasta: '>duplicate\nACGT\n' })) => prepareLosatSourceBatches({
  sequences: records, specs: jobs, getEntry, buildArgs: () => ['--task', 'blastn'], hashText, protein: false
});
const plan = await planFor();
assert.equal(plan.batches.length, 4, 'eight records in two sources require four LOSAT jobs');
assert.equal(plan.batches.reduce((sum, batch) => sum + batch.specs.length, 0), 64);
for (const batch of plan.batches) {
  assert.equal(batch.query.ids.size, batch.query.indexes.length);
  assert.equal(batch.subject.ids.size, batch.subject.indexes.length);
  assert.match(batch.searchContext, /^[0-9a-f]{64}$/);
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
const changed = await planFor(sequences, specs, async (index) => ({ fasta: `>duplicate\n${index === 7 ? 'ACGA' : 'ACGT'}\n` }));
assert.notEqual(plan.batches[1].searchContext, changed.batches[1].searchContext, 'database contents change cache scope');
const reordered = [...sequences].reverse();
const replay = await planFor(reordered);
assert.deepEqual(replay.batches.map((batch) => batch.searchContext).sort(), plan.batches.map((batch) => batch.searchContext).sort(), 'record ordering does not change the searched source sets');

const codePlan = await prepareLosatSourceBatches({
  sequences, specs: adjacentSpecs, getEntry: async () => ({ fasta: '>duplicate\nACGT\n' }),
  buildArgs: (query) => ['--query-gencode', query === 0 ? '4' : '11', '--db-gencode', '11'],
  hashText, protein: false
});
assert.equal(codePlan.batches.length, 2, 'conflicting explicit translation tables require compatible source batches');
assert.deepEqual(codePlan.batches.map((batch) => batch.query.indexes.length).sort(), [1, 5]);

const noSelf = await prepareLosatSourceBatches({
  sequences, specs: specs.filter(({ queryIndex, subjectIndex }) => queryIndex !== subjectIndex),
  getEntry: async (index) => ({ fasta: `>protein-${index}\nMKK\n` }),
  buildArgs: () => ['--max-target-seqs', '5'], hashText, protein: true,
  excludeSelfComparisons: true
});
for (const batch of noSelf.batches) {
  assert(!batch.query.indexes.some((index) => batch.subject.indexes.includes(index)),
    'Collinear OFF must never submit a within-record search, including multi-record sources');
}
assert.equal(noSelf.batches.length, 34);
assert.equal(noSelf.batches.reduce((sum, batch) => sum + batch.specs.length, 0), 56);
