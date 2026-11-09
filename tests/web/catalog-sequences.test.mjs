// The catalog-sequence door (D5): the one place the page reads the sequences the stored catalog does not hold. One
// reply per (shown Result, feature or span) is kept; a reply that arrives after the shown Result changed is discarded
// (null); a failed read is forgotten so Retry reads again; a Result without its request (a Result-only Session 27-39)
// gives the reason the sequences are read after the next Generate.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import {
  createCatalogSequences,
  featureSequenceKey,
  featureSequenceReads,
  matchSequenceRecordsOf,
  readRowsFasta
} from '../../gbdraw/web/js/services/catalog-sequences.js';
import { normalizeCaughtError } from '../../gbdraw/web/js/utils/error-normalization.js';

const committedOf = (fastaFiles = []) => ({
  renderRequest: { mode: 'circular', diagramOptions: { conservationFastaFiles: fastaFiles } },
  resources: {}
});
const featureReply = ({ recordKey, biologicalFeatureId }) => ({
  recordKey, biologicalFeatureId, nucleotideSequence: `ATG${biologicalFeatureId}`.toUpperCase().replace(/[^ACGT]/g, 'A'),
  aminoAcidSequence: 'M', sequenceWarnings: []
});

const doorWith = (initialCommitted) => {
  let committed = initialCommitted;
  const calls = { features: [], spans: [] };
  let failNext = null;
  /** @type {Array<() => void>} */
  const holds = [];
  let holding = false;
  const helpers = {
    readFeatureSequences: async (payload) => {
      calls.features.push(payload.features);
      if (holding) await new Promise((resolve) => holds.push(resolve));
      if (failNext) {
        const error = failNext;
        failNext = null;
        throw error;
      }
      return { result: { features: payload.features.map(featureReply) } };
    },
    readMatchSequences: async (payload) => {
      calls.spans.push(payload.spans);
      return { result: { spans: payload.spans.map((span) => ({
        status: 'ok', reason: null, recordId: 'REC', orientation: span.start <= span.end ? '+' : '-',
        length: Math.abs(span.end - span.start) + 1, sequence: 'ACGT'
      })) } };
    },
    readEmbedSequences: async (payload) => ({ result: payload.item })
  };
  const door = createCatalogSequences({ getCommittedSession: () => committed, helpers });
  return {
    door,
    calls,
    setCommitted: (value) => { committed = value; },
    failNextWith: (error) => { failNext = error; },
    hold: () => { holding = true; },
    release: () => { holding = false; holds.splice(0).forEach((resolve) => resolve()); }
  };
};

const kinase = { recordKey: 'rec-1', biologicalFeatureId: 'kinase', bases: 12 };
const ligase = { recordKey: 'rec-1', biologicalFeatureId: 'ligase', bases: 30 };

test('a feature read once per shown Result is kept; another Result reads again', async () => {
  const first = committedOf();
  const harness = doorWith(first);
  const one = await harness.door.readFeatures([kinase, ligase]);
  assert.equal(one.get(featureSequenceKey(kinase)).biologicalFeatureId, 'kinase');
  await harness.door.readFeatures([kinase]);
  assert.equal(harness.calls.features.length, 1, 'the kept reply answers the second read');
  harness.setCommitted(committedOf());
  await harness.door.readFeatures([kinase]);
  assert.equal(harness.calls.features.length, 2, 'a new Result reads again');
});

test('a reply that arrives after the shown Result changed is discarded', async () => {
  const harness = doorWith(committedOf());
  harness.hold();
  const pending = harness.door.readFeatures([kinase]);
  await Promise.resolve();
  harness.setCommitted(committedOf());
  harness.release();
  assert.equal(await pending, null);
  const spans = doorWith(committedOf());
  const pendingSpans = spans.door.readSpans([{ kind: 'record', recordKey: 'rec-1', start: 1, end: 4 }]);
  spans.setCommitted(committedOf());
  assert.equal(await pendingSpans, null);
});

test('a failed read is forgotten, so Retry reads again', async () => {
  const harness = doorWith(committedOf());
  harness.failNextWith(Object.assign(new Error('boom'), { code: 'RESOURCE_INVALID', stage: 'helper' }));
  await assert.rejects(() => harness.door.readFeatures([kinase]));
  const again = await harness.door.readFeatures([kinase]);
  assert.ok(again.get(featureSequenceKey(kinase)));
  assert.equal(harness.calls.features.length, 2);
});

test('a Result without its request reads nothing and says the sequences are read after the next Generate', async () => {
  const harness = doorWith(null);
  const error = await harness.door.readFeatures([kinase]).catch((caught) => caught);
  assert.equal(error.code, 'SEQUENCE_REQUIRES_GENERATE');
  assert.equal(harness.calls.features.length, 0);
  const shown = normalizeCaughtError(error, { operation: 'readFeatureSequences', stage: 'helper' });
  assert.match(shown.summary, /read after the next Generate/);
  assert.equal(matchSequenceRecordsOf({ recordKeys: ['rec-1'], biologicalFeatures: [] }, null), null);
});

test('reads name catalog rows by identity, each once, and rows without an identity are left out', () => {
  const rows = [
    { recordKey: 'rec-1', biologicalFeatureId: 'kinase', start: 0, end: 12 },
    { recordKey: 'rec-1', biologicalFeatureId: 'kinase', start: 0, end: 12 },
    { recordKey: '', biologicalFeatureId: 'orphan', start: 0, end: 5 }
  ];
  assert.deepEqual(featureSequenceReads(rows), [{ recordKey: 'rec-1', biologicalFeatureId: 'kinase', bases: 12 }]);
});

test('the records of a match popup come from the shown catalog item and the committed request', () => {
  const item = {
    recordKeys: ['rec-a', 'rec-b'],
    biologicalFeatures: [
      { recordKey: 'rec-a', record_id: 'NC_1' },
      { recordKey: 'rec-b', record_id: 'NC_2' },
      { recordKey: 'rec-b', record_id: 'NC_2' }
    ]
  };
  const records = matchSequenceRecordsOf(item, committedOf([{ resourceId: 'r1' }, null]));
  assert.deepEqual(records.unitKeys, ['rec-a', 'rec-b']);
  assert.deepEqual(records.unitKeysOfRecordId('NC_2'), ['rec-b']);
  assert.deepEqual(records.unitKeysOfRecordId('NC_9'), []);
  assert.equal(records.comparisonSupplied(0), true);
  assert.equal(records.comparisonSupplied(1), false, 'a never-supplied FASTA is known without reading');
  assert.equal(records.comparisonSupplied(2), false);
});

test('the FASTA of catalog rows reads every row once and leaves out rows without a sequence of the kind', async () => {
  const asked = [];
  const readFeatures = async (reads) => {
    asked.push(reads);
    return new Map(reads.map((read) => [featureSequenceKey(read), {
      ...read, nucleotideSequence: read.biologicalFeatureId === 'kinase' ? 'ATGAAA' : '', aminoAcidSequence: '',
      sequenceWarnings: []
    }]));
  };
  const rows = [
    { read: kinase, feature: { id: 'kinase', record_id: 'REC', start: 0, end: 6, strand: 1, type: 'CDS' } },
    { read: ligase, feature: { id: 'ligase', record_id: 'REC', start: 10, end: 40, strand: 1, type: 'CDS' } },
    { read: null, feature: { id: 'orphan' } }
  ];
  const text = await readRowsFasta(rows, 'nt', readFeatures);
  assert.equal(asked.length, 1);
  assert.equal(asked[0].length, 2);
  assert.match(text, /^>.*\nATGAAA\n$/);
  assert.equal(await readRowsFasta(rows, 'aa', readFeatures), '');
  assert.equal(await readRowsFasta(rows, 'nt', async () => null), null, 'a discarded read gives null');
  assert.equal(await readRowsFasta([{ read: null, feature: {} }], 'nt', readFeatures), '');
});
