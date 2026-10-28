import assert from 'node:assert/strict';
import test from 'node:test';
import { createHash } from 'node:crypto';

const {
  convertStoredProteinResultFrames, storesDrawnFrameProteinResults, withConvertedResources
} = await import('../../gbdraw/web/js/services/main-session-comparison-frame.js');
const { MAIN_DISPLAY_FRAME_SESSION_VERSION } = await import('../../gbdraw/web/js/contracts/schema.generated.js');

const base64 = (text) => Buffer.from(text, 'utf8').toString('base64');
const TABLE = 'q\ts\t99\t10\t0\t0\t1\t10\t1\t10\t0\t20\n';
const digest = (text) => createHash('sha256').update(text, 'utf8').digest('hex');

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
    table: { encoding: 'base64', data: base64(TABLE), name: 'table.tsv', size: TABLE.length, checksum: digest(TABLE) }
  }
});

// Record 0 is reversed; a precomputed table names its two records.
const protein = (version, reverse = true, kind = 'orthogroupResult') => {
  const data = session({ version, reverse });
  data.renderRequest.comparisons = [kind === 'precomputedProteinComparison'
    ? { kind, resourceId: 'table', encoding: 'canonicalTsv', queryRecordIndex: 0, subjectRecordIndex: 1 }
    : { kind, resourceId: 'table' }];
  return data;
};

test('a Session up to 47 with protein results of a reversed record is named for conversion (OV-1145)', () => {
  assert.equal(storesDrawnFrameProteinResults(protein(47)), true);
  assert.equal(storesDrawnFrameProteinResults(protein(46, true, 'collinearityResult')), true);
  assert.equal(storesDrawnFrameProteinResults(protein(47, true, 'precomputedProteinComparison')), true);
  assert.equal(storesDrawnFrameProteinResults(protein(48)), false);
  assert.equal(storesDrawnFrameProteinResults(protein(47, false)), false);
  assert.equal(storesDrawnFrameProteinResults(protein(47, true, 'nucleotideBlast')), false);
  assert.equal(storesDrawnFrameProteinResults(session({ version: 47, mode: 'circular' })), false);
  // The other mode's committed request is read too.
  const other = { format: 'gbdraw-session', version: 47, renderRequest: { mode: 'circular' },
    otherModeResult: { renderRequest: protein(47).renderRequest }, resources: protein(47).resources };
  assert.equal(storesDrawnFrameProteinResults(other), true);
});

// B-LOAD (D-128): only a result that touches a reversed record is converted,
// so every other Session's preview-only Load starts no Python Worker.
test('a protein result that touches no reversed record names no conversion (OV-1145)', () => {
  const third = (data) => {
    data.renderRequest.records.push({ source: { kind: 'genbank', resourceId: 's' } });
    return data;
  };
  // A table between the two forward records; record 0 is reversed.
  const forwardTable = third(protein(47, true, 'precomputedProteinComparison'));
  Object.assign(forwardTable.renderRequest.comparisons[0], { queryRecordIndex: 1, subjectRecordIndex: 2 });
  assert.equal(storesDrawnFrameProteinResults(forwardTable), false);
  const indexedAsKey = protein(47, true, 'precomputedProteinComparison');
  indexedAsKey.renderRequest.comparisons[0].queryRecordIndex = '0';
  assert.equal(storesDrawnFrameProteinResults(indexedAsKey), true);
  // A typed result touches the records its marker's pairs name, or every record when none does.
  const marked = (pairs) => {
    const data = third(protein(47));
    data.renderRequest.comparisons.push({ kind: 'generatedProteinComparison', mode: 'none', pairs, settings: {} });
    return data;
  };
  assert.equal(storesDrawnFrameProteinResults(marked([{ queryRecordIndex: 1, subjectRecordIndex: 2 }])), false);
  assert.equal(storesDrawnFrameProteinResults(marked([{ queryRecordIndex: 0, subjectRecordIndex: 1 }])), true);
  assert.equal(storesDrawnFrameProteinResults(marked([])), true);
  // A reversed record whose file the Session does not hold is the request reader's to drop (OV-869).
  for (const source of [{ kind: 'genbank', resourceId: 'q' }, { kind: 'gffFasta', gffResourceId: 'q', fastaResourceId: 's' }]) {
    const unheld = protein(47);
    unheld.renderRequest.records[0].source = source;
    delete unheld.resources.q;
    assert.equal(storesDrawnFrameProteinResults(unheld), false, source.kind);
  }
});

const proteinConverter = (resources) => {
  const calls = [];
  const convert = async (request) => {
    calls.push(request);
    return { result: { resources } };
  };
  return { calls, convert };
};

test('Python converts each drawing whose protein results touch a reversed record, once (OV-1145)', async () => {
  const shown = protein(47);
  const kept = protein(47, true, 'precomputedProteinComparison');
  const forward = protein(47, false);
  const { calls, convert } = proteinConverter({ table: { data: base64('F\n'), size: 2 } });
  const converted = await convertStoredProteinResultFrames({ version: 47, inDrawnFrame: true }, [shown, null, forward, kept], { convert });
  assert.deepEqual(calls, [
    { version: 47, canonicalRequest: shown.renderRequest, resources: shown.resources },
    { version: 47, canonicalRequest: kept.renderRequest, resources: kept.resources }
  ]);
  assert.deepEqual([...converted], [['table', { data: base64('F\n'), size: 2 }]]);
  // A Session 48 and a drawing that touches no reversed record start nothing.
  assert.equal((await convertStoredProteinResultFrames({ version: 48, inDrawnFrame: false }, [shown], { convert })).size, 0);
  assert.equal((await convertStoredProteinResultFrames({ version: 47, inDrawnFrame: true }, [forward], { convert: null })).size, 0);
  assert.equal(calls.length, 2);
});

test('a conversion that fails or returns unreadable bytes changes nothing and is reported (OV-1145)', async () => {
  const shown = protein(47);
  const failure = (code) => ({ code, stage: 'helper', operation: 'convertStoredProteinResultFrames' });
  await assert.rejects(convertStoredProteinResultFrames({ version: 47, inDrawnFrame: true }, [shown], { convert: null }), failure('HELPER_PROTOCOL'));
  for (const resources of [null, { table: { data: 1, size: 2 } }, { table: { data: base64('F\n'), size: -1 } }]) {
    await assert.rejects(convertStoredProteinResultFrames({ version: 47, inDrawnFrame: true }, [shown], proteinConverter(resources)),
      failure('RESULT_INVALID'), JSON.stringify(resources));
  }
  // Two drawings that store one resource ID must agree on its bytes.
  let call = 0;
  const convert = async () => ({ result: { resources: { table: { data: base64(call++ ? 'B\n' : 'A\n'), size: 2 } } } });
  await assert.rejects(convertStoredProteinResultFrames({ version: 47, inDrawnFrame: true }, [shown, protein(47)], { convert }), failure('RESULT_INVALID'));
});

test('the converted bytes replace a resource under its stored ID with its size and checksum (OV-1145)', async () => {
  const stored = session().resources;
  const converted = new Map([['table', { data: base64('F\n'), size: 2 }], ['absent', { data: base64('x'), size: 1 }]]);
  const next = await withConvertedResources(stored, converted);
  assert.deepEqual(Object.keys(next), Object.keys(stored));
  assert.equal(next.q, stored.q);
  assert.deepEqual(next.table, {
    encoding: 'base64', data: base64('F\n'), name: 'table.tsv', size: 2,
    checksum: createHash('sha256').update('F\n').digest('hex')
  });
  assert.equal(stored.table.data, base64(TABLE), 'the stored resources are not changed');
  const { checksum: _checksum, ...unchecked } = stored.table;
  assert.equal(Object.hasOwn((await withConvertedResources({ ...stored, table: unchecked }, converted)).table, 'checksum'), false);
  assert.equal(await withConvertedResources(stored, new Map()), stored);
});
