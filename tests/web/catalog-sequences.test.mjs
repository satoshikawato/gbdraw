// The catalog-sequence door (D5): the one place the page reads the sequences the stored catalog does not hold. One
// reply per (shown Result, feature or span) is kept; a reply that arrives after the shown Result changed is discarded
// (null); a failed read is forgotten so Retry reads again; a Result without its request (a Result-only Session 27-39)
// gives the reason the sequences are read after the next Generate.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import {
  createCatalogSequences,
  embedReadParts,
  featureSequenceKey,
  featureSequenceReads,
  matchSequenceRecordsOf,
  readRowsFasta
} from '../../gbdraw/web/js/services/catalog-sequences.js';
import { normalizeCaughtError } from '../../gbdraw/web/js/utils/error-normalization.js';
import { embeddedCatalogItemJson } from '../../gbdraw/web/js/services/standalone-interactivity.js';

const committedOf = (fastaFiles = []) => ({
  renderRequest: { mode: 'circular', diagramOptions: { conservationFastaFiles: fastaFiles } },
  resources: {}
});
const featureReply = ({ recordKey, biologicalFeatureId }) => ({
  recordKey, biologicalFeatureId, nucleotideSequence: `ATG${biologicalFeatureId}`.toUpperCase().replace(/[^ACGT]/g, 'A'),
  aminoAcidSequence: 'M', sequenceWarnings: []
});

// A stand-in for Python's export encoder (`embed_catalog_part`): every row gets its sequence; each drawn record has
// one whole-record source (its index the record's place) and a comparison source comes last.
const EMBED_SOURCES = [{ key: 'linear:record:0', recordIndex: 0 }, { key: 'linear:record:1', recordIndex: 1 },
  { key: 'homology:comparison:0:C', sourceIndex: 0 }];
// What the encoder makes of a part, and the reply that carries it as the text the file embeds.
const embedPartObjects = ({ item, part }) => ({
  biologicalFeatures: item.biologicalFeatures.map((row) => ({ ...row, nucleotide_sequence: `N${row.biologicalFeatureId}` })),
  sequenceSources: EMBED_SOURCES.map((source, index) => ({ index, source })).filter(({ source }) => (
    source.recordIndex === undefined ? part.comparisonSources : part.sourcesOf.includes(item.recordKeys[source.recordIndex]))),
  sourceCount: EMBED_SOURCES.length
});
const embedPartReply = (payload) => {
  const { biologicalFeatures, sequenceSources, sourceCount } = embedPartObjects(payload);
  return { biologicalFeatures: biologicalFeatures.map((row) => JSON.stringify(row)).join(','), rowCount: biologicalFeatures.length,
    sequenceSources: sequenceSources.map(({ index, source }) => ({ index, source: JSON.stringify(source) })), sourceCount };
};

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
    readEmbedSequences: async (payload) => ({ result: embedPartReply(payload) })
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

test('the Worker receives plain spans and features even when the popup state holds them as reactive proxies', async () => {
  // The match popup's entries live in Vue state; postMessage cannot clone a Proxy (DataCloneError), so the door
  // sends copies holding only the identity and range fields.
  const sent = [];
  const door = createCatalogSequences({
    getCommittedSession: () => committed,
    helpers: {
      readFeatureSequences: async (payload) => { sent.push(structuredClone(payload)); return { result: { features: payload.features.map(featureReply) } }; },
      readMatchSequences: async (payload) => {
        sent.push(structuredClone(payload));
        return { result: { spans: payload.spans.map(() => ({ status: 'ok', reason: null, recordId: 'REC', orientation: '+', length: 4, sequence: 'ACGT' })) } };
      },
      readEmbedSequences: async (payload) => ({ result: embedPartReply(payload) })
    }
  });
  const committed = committedOf();
  const proxied = (value) => new Proxy(value, {});
  const record = proxied({ kind: 'record', recordKey: 'rec-1', start: 1, end: 4, label: 'Query', available: true });
  const comparison = proxied({ kind: 'comparison', sourceIndex: 0, recordId: 'REC', start: 4, end: 1, label: 'Subject' });
  await door.readSpans([record, comparison]);
  await door.readFeatures([proxied({ ...kinase, label: 'kinase' })]);
  assert.deepEqual(sent.map((payload) => payload.spans ?? payload.features), [
    [{ kind: 'record', recordKey: 'rec-1', start: 1, end: 4 }, { kind: 'comparison', sourceIndex: 0, recordId: 'REC', start: 4, end: 1 }],
    [{ recordKey: 'rec-1', biologicalFeatureId: 'kinase' }]
  ]);
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

const exportItem = {
  resultIndex: 0, resultName: 'd.svg', recordKeys: ['rec-a', 'rec-b'], features: [],
  biologicalFeatures: [
    { recordKey: 'rec-a', biologicalFeatureId: 'a1', start: 0, end: 10 },
    { recordKey: 'rec-a', biologicalFeatureId: 'a2', start: 10, end: 30 },
    { recordKey: 'rec-b', biologicalFeatureId: 'b1', start: 0, end: 5 }
  ],
  orthogroups: [], annotations: [], comparisonMatches: [{ id: 'm1' }]
};

test('an export reads one part per record, never the whole item, with each source once', () => {
  const parts = embedReadParts(exportItem, 25);
  assert.deepEqual(parts.map(({ item, part }) => [item.biologicalFeatures.map((row) => row.biologicalFeatureId), part]), [
    [['a1'], { sourcesOf: ['rec-a'], comparisonSources: false }],
    [['a2'], { sourcesOf: [], comparisonSources: false }],
    [['b1'], { sourcesOf: ['rec-b'], comparisonSources: true }]
  ]);
  assert.deepEqual(Object.keys(parts[0].item), ['recordKeys', 'comparisonMatches', 'biologicalFeatures']);
  // A record without rows gives its sources with the last part; an item without rows is one part.
  const noRows = embedReadParts({ ...exportItem, biologicalFeatures: exportItem.biologicalFeatures.slice(2) }, 25);
  assert.deepEqual(noRows.map(({ part }) => part), [{ sourcesOf: ['rec-b', 'rec-a'], comparisonSources: true }]);
  assert.deepEqual(embedReadParts({ ...exportItem, biologicalFeatures: [] }, 25).map(({ item, part }) => [item.biologicalFeatures, part]),
    [[[], { sourcesOf: ['rec-a', 'rec-b'], comparisonSources: true }]]);
});

test('the parts written as they arrive are the item Python encodes whole, byte for byte', async () => {
  const sent = [];
  const committed = committedOf();
  const door = createCatalogSequences({
    getCommittedSession: () => committed, featureReadBases: 25,
    helpers: { readEmbedSequences: async (payload) => { sent.push(payload.part); return { result: embedPartReply(payload) }; } }
  });
  const sequences = await door.readEmbedJson(exportItem);
  assert.equal(sent.length, 3);
  const whole = embedPartObjects({ item: exportItem, part: { sourcesOf: ['rec-a', 'rec-b'], comparisonSources: true } });
  assert.equal(embeddedCatalogItemJson(exportItem, sequences), JSON.stringify({
    ...exportItem, biologicalFeatures: whole.biologicalFeatures, sequenceSources: whole.sequenceSources.map(({ source }) => source)
  }));
});

test('an export read the shown Result outlived gives null; a source missing or given twice is an invalid reply', async () => {
  let committed = committedOf();
  const outlived = createCatalogSequences({
    getCommittedSession: () => committed, featureReadBases: 25,
    helpers: { readEmbedSequences: async (payload) => { committed = committedOf(); return { result: embedPartReply(payload) }; } }
  });
  assert.equal(await outlived.readEmbedJson(exportItem), null);
  for (const reply of [
    (payload) => ({ ...embedPartReply(payload), sequenceSources: [] }),
    (payload) => ({ ...embedPartReply(payload), sequenceSources: [{ index: 0, source: '{"key":"x"}' }] }),
    (payload) => ({ ...embedPartReply(payload), sequenceSources: [{ index: 0, source: { key: 'x' } }] }),
    (payload) => ({ ...embedPartReply(payload), rowCount: 0 }),
    (payload) => ({ ...embedPartReply(payload), biologicalFeatures: [] })
  ]) {
    const fixed = committedOf();
    const door = createCatalogSequences({
      getCommittedSession: () => fixed, featureReadBases: 25,
      helpers: { readEmbedSequences: async (payload) => ({ result: reply(payload) }) }
    });
    const error = await door.readEmbedJson(exportItem).catch((caught) => caught);
    assert.equal(normalizeCaughtError(error).code, 'RESULT_INVALID');
  }
});
