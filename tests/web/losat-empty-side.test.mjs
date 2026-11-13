// A LOSATP source batch whose query or subject side holds no sequence (records without a CDS protein, OV-3493) is not
// searchable: its pairs find no hit and LOSAT, which refuses an empty database, does not run. Python's
// `LosatBatch.searchable` gives the same flag (tests/fixtures/losat_job_plan_cases.json holds both to it), and
// `_search_losatp_record_pairs` answers such a batch the same way.
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import test from 'node:test';

const { prepareLosatSourceBatches, runSearchableLosatSourceJobs } = await import(
  new URL('../../gbdraw/web/js/services/linear-sources.js', import.meta.url).href
);

const hashText = async (/** @type {string} */ text) => createHash('sha256').update(text).digest('hex');

test('a batch with a side of records without proteins is not searchable', async () => {
  const fastas = ['>p1\nMK\n', '', '>p3\nMK\n'];
  // Three single-record files: each record is a side of its own.
  const sequences = fastas.map((_, index) => ({ uid: `record-${index + 1}`, gb: { name: `${index}.gbk` }, gff: null, fasta: null }));
  const specs = [{ queryIndex: 0, subjectIndex: 1 }, { queryIndex: 1, subjectIndex: 2 }, { queryIndex: 0, subjectIndex: 2 }];
  const { batches } = await prepareLosatSourceBatches({
    sequences, specs, getEntry: async (index) => ({ fasta: fastas[index] }),
    getDigest: async (index) => hashText(fastas[index]), buildArgs: () => [], hashText, protein: true
  });
  assert.deepEqual(batches.map((batch) => [batch.query.indexes, batch.subject.indexes, batch.searchable]), [
    [[0], [1], false], [[1], [2], false], [[0], [2], true]
  ]);
});

test('only searchable jobs are searched; the others are answered with no rows', async () => {
  const jobs = [
    { cacheKey: 'a', batch: { searchable: true } },
    { cacheKey: 'b', batch: { searchable: false } },
    { cacheKey: 'c', batch: { searchable: true } }
  ];
  /** @type {string[][]} */
  const ran = [];
  /** @type {unknown[]} */
  const passed = [];
  const results = await runSearchableLosatSourceJobs(jobs, async (searched, options) => {
    ran.push(searched.map(({ cacheKey }) => cacheKey));
    passed.push(options, ...searched.map((job) => 'batch' in job));
    return searched.map(({ cacheKey }) => ({ cacheKey, text: 'p1\tp2\t100\t2\t0\t0\t1\t2\t1\t2\t1e-5\t10\n' }));
  }, { concurrency: 2 });
  assert.deepEqual(ran, [['a', 'c']]);
  // The runner gets the options, and the jobs without their source-job bookkeeping.
  assert.deepEqual(passed, [{ concurrency: 2 }, false, false]);
  assert.deepEqual(Object.fromEntries(results.map(({ cacheKey, text }) => [cacheKey, text !== ''])), { a: true, b: false, c: true });
});

test('jobs that are all unsearchable run no search at all', async () => {
  const results = await runSearchableLosatSourceJobs(
    [{ cacheKey: 'b', batch: { searchable: false } }],
    async () => { throw new Error('SEARCH ATTEMPTED'); },
    {}
  );
  assert.deepEqual(results, [{ cacheKey: 'b', text: '' }]);
});
