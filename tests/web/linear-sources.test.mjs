import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { groupLinearSourceRecords, prepareLosatSourceBatches, splitLosatSourceResult } from '../../gbdraw/web/js/app/linear-sources.js';

const hashText = async (value) => createHash('sha256').update(value).digest('hex');
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
