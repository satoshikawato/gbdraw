import assert from 'node:assert/strict';
import test from 'node:test';

const {
  MAIN_DISPLAY_FRAME_SESSION_VERSION,
  convertMainSessionComparisonFrames
} = await import('../../gbdraw/web/js/services/main-session-comparison-frame.js');

const base64 = (text) => Buffer.from(text, 'utf8').toString('base64');
const TABLE = 'q\ts\t99\t10\t0\t0\t1\t10\t1\t10\t0\t20\n';

// A Linear Session that stored the comparison table of a reversed record (its
// region is reverse-complemented) in the display frame.
const session = ({ version = MAIN_DISPLAY_FRAME_SESSION_VERSION, mode = 'linear', reverse = true, resources } = {}) => ({
  format: 'gbdraw-session',
  version,
  renderRequest: {
    mode,
    records: [
      { source: { kind: 'genbank', resourceId: 'q' }, region: { start: 11, end: 30, reverseComplement: reverse } },
      { source: { kind: 'genbank', resourceId: 's' } }
    ],
    comparisons: [{ kind: 'nucleotideBlast', resourceId: 'table', queryRecordIndex: 0, subjectRecordIndex: 1 }]
  },
  resources: resources ?? {
    q: { encoding: 'base64', data: base64('LOCUS q'), name: 'q.gb', size: 7 },
    s: { encoding: 'base64', data: base64('LOCUS s'), name: 's.gb', size: 7 },
    table: { encoding: 'base64', data: base64(TABLE), name: 'table.tsv', size: TABLE.length, checksum: 'kept' }
  }
});

const recordingConverter = () => {
  const calls = [];
  const convertTable = async (request) => {
    calls.push(request);
    return { result: { tsv: 'converted\n' } };
  };
  return { calls, convertTable };
};
const NO_DISCOVERY = {
  sequence: async () => { throw new Error('no record discovery expected'); },
  gffFasta: async () => { throw new Error('no record discovery expected'); }
};

test('a stored table of a reversed record is rewritten in the search frame', async () => {
  const { calls, convertTable } = recordingConverter();
  const stored = session();
  const read = await convertMainSessionComparisonFrames(stored, { convertTable, discover: NO_DISCOVERY });
  assert.deepEqual(calls, [{ tableText: TABLE, queryFrame: { length: 20, reverse: true }, subjectFrame: { length: 0, reverse: false } }]);
  assert.deepEqual(read.resources.table, {
    encoding: 'base64', data: base64('converted\n'), name: 'table.tsv', size: 10, checksum: 'kept'
  });
  assert.equal(read.resources.q, stored.resources.q);
  assert.equal(stored.resources.table.data, base64(TABLE), 'the stored document is not changed');
});

test('a Session the conversion does not read is returned as it is', async () => {
  const { calls, convertTable } = recordingConverter();
  for (const stored of [
    session({ version: MAIN_DISPLAY_FRAME_SESSION_VERSION + 1 }),
    session({ mode: 'circular' }),
    session({ reverse: false })
  ]) {
    assert.equal(await convertMainSessionComparisonFrames(stored, { convertTable, discover: NO_DISCOVERY }), stored);
  }
  assert.deepEqual(calls, []);
});

// `readSessionDocument` drops a resource whose descriptor cannot be read; the
// request's references to it are its reader's (`readSessionRequests`), which
// reads the request after this conversion.
test('a comparison whose table the Session does not hold is left to the request reader', async () => {
  const { calls, convertTable } = recordingConverter();
  const stored = session();
  delete stored.resources.table;
  assert.equal(await convertMainSessionComparisonFrames(stored, { convertTable, discover: NO_DISCOVERY }), stored);
  assert.deepEqual(calls, []);
});

test('a stored record index is read as a property key, as the stored request was always read', async () => {
  for (const [index, convertedCount] of [['0', 1], [null, 0], [0.5, 0]]) {
    const { calls, convertTable } = recordingConverter();
    const stored = session();
    stored.renderRequest.comparisons[0].queryRecordIndex = index;
    await convertMainSessionComparisonFrames(stored, { convertTable, discover: NO_DISCOVERY });
    assert.equal(calls.length, convertedCount, JSON.stringify(index));
  }
});

// OV-869: a reversed record with no region reads its length from its record
// file; one the Session does not hold (`readSessionDocument` dropped it) leaves
// the comparison to the request reader, which drops it with its row.
test('a comparison whose reversed record file the Session does not hold is left to the request reader (OV-869)', async () => {
  for (const source of [{ kind: 'genbank', resourceId: 'q' }, { kind: 'gffFasta', gffResourceId: 'q', fastaResourceId: 's' }]) {
    const { calls, convertTable } = recordingConverter();
    const stored = session();
    stored.renderRequest.records[0] = { source, presentation: { reverseComplement: true } };
    delete stored.resources.q;
    assert.equal(await convertMainSessionComparisonFrames(stored, { convertTable, discover: NO_DISCOVERY }), stored, source.kind);
    assert.deepEqual(calls, []);
  }
});
