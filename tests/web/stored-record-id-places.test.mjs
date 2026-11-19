// OV-3311: a stored record whose selector is a record ID that a record text reads as a `#N` (`#2`, `#+2`) draws the
// record of that ID after Load, as Python draws it. A Linear card, the Circular Record select and a Circular record
// range hold their record as text, so Load reads the record's place in its file from the record list
// (`readStoredRecordIdPlaces`, services/record-discovery.js) into the stored-leaf readings, through the resources
// it adopted, and the projection gives the text that place (`#k`, the text the Record select gives such a record).
// The card's own read of the file later is the same read. A Session without such a record reads no file; one whose
// record Load cannot place gets a Load repair row.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { File } from 'node:buffer';
import test from 'node:test';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { pythonLegacySessionUpgrader } from './helpers/python-evaluator.mjs';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';
import { readSavedSession, sessionFileOf } from './helpers/session-documents.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }) };
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
installSessionImportWorker();
const metrics = [];
globalThis.__GBDRAW_TEST_HOOKS__ = { onStructuralMetric: (metric) => metrics.push(metric) };
const decodes = () => metrics.filter(({ name }) => name === 'base64DecodeCount').length;

const { setLegacySessionUpgrader } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
setLegacySessionUpgrader(pythonLegacySessionUpgrader());
const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: { overrides: payload.configOverrides } }));
const { CANONICAL_REQUEST_SCHEMA, projectCanonicalSessionRequest } = await import('../../gbdraw/web/js/services/session-request.js');
const { circularDepthRowsAfterRead, committedRecordIdentitySources, readStoredRecordIdPlaces, readOncePerFile } = await import('../../gbdraw/web/js/services/record-discovery.js');
const { committedCardRecords, drawnRecordIdentities, recordKeyMoves, recordSelectorNumber } = await import('../../gbdraw/web/js/services/record-draw-selection.js');
const { buildDisambiguatedRecordEntries } = await import('../../gbdraw/web/js/services/record-options.js');
const { adoptCurrentSessionResources, createSessionResourceFileView } = await import('../../gbdraw/web/js/services/session-resource-backing.js');
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

// Records named `first`, `second`, `#2` and `#+2`, in this order: the ID `#2` is the third record and `#+2` the
// fourth, so reading either text as a `#N` draws another record.
const NAMES = ['first', 'second', '#2', '#+2'];
const genbankOf = (names) => names.map((name, index) => [
  `LOCUS       ${name.padEnd(16)} ${10 + index} bp    DNA     linear   UNK 01-JAN-1980`,
  `DEFINITION  ${name}.`, `ACCESSION   ${name}`, `VERSION     ${name}`, 'FEATURES             Location/Qualifiers',
  'ORIGIN', `        1 ${'a'.repeat(10 + index)}`, '//', ''
].join('\n')).join('');
const genbank = genbankOf(NAMES);
const descriptor = (text, name) => ({
  kind: 'genbank', name, type: 'text/plain', size: Buffer.byteLength(text), lastModified: 0,
  encoding: 'base64', data: Buffer.from(text).toString('base64')
});
const requestRecord = (recordKey, selector, resourceId = 'named') => ({
  recordKey, cardinality: 'exactly_one', source: { kind: 'genbank', resourceId }, selector, windows: [],
  presentation: { label: null, subtitle: null, reverseComplement: false, gridRow: null, gridColumn: null },
  display: { isCircular: null, startCoordinate: null }
});
const request = (mode, records) => ({
  schema: CANONICAL_REQUEST_SCHEMA, mode, grouping: 'single', records,
  diagramOptions: { featurePlacements: [], featureOverrides: [] }, output: { prefix: 'ov3311' }
});
// The record a record text names, as Generate reads it (`selectorPayload`): a `#N` counts records, any other text
// is a record ID.
const drawnName = (text) => {
  const number = recordSelectorNumber(text);
  return number === null ? text : NAMES[number - 1];
};

test('a stored record ID that reads as #N draws the record of that ID after Load, with one decode and one listing', async () => {
  const resources = { named: descriptor(genbank, 'named.gb') };
  const table = adoptCurrentSessionResources(resources);
  const linear = request('linear', [
    requestRecord('record-1', { kind: 'recordId', value: '#2' }),
    requestRecord('record-2', { kind: 'recordId', value: '#+2' }),
    requestRecord('record-3', { kind: 'recordId', value: 'first' })
  ]);
  const circular = request('circular', [requestRecord('record-1', { kind: 'recordId', value: '#+2' })]);
  metrics.length = 0;
  const { places, unresolved } = await readStoredRecordIdPlaces(
    [{ at: 'drawings[0].', renderRequest: linear }, { at: 'drawings[1].', renderRequest: circular }],
    (resourceId) => createSessionResourceFileView(table, resourceId));
  assert.deepEqual(unresolved, []);
  assert.equal(decodes(), 1);
  // The card's own read of its file, through another view of the same resource, is the read Load made.
  let listings = 0;
  await readOncePerFile([createSessionResourceFileView(table, 'named')], 'records:genbank', async () => { listings += 1; return []; });
  assert.deepEqual({ listings, decodes: decodes() }, { listings: 0, decodes: 1 });

  const project = (renderRequest) => projectCanonicalSessionRequest({
    renderRequest, resources, sessionResourceTable: table, initializeCliInputs: true, leafReadings: places
  });
  const cards = project(linear).files.linearSeqs.map((seq) => seq.region_record_id);
  assert.deepEqual(cards, ['#3', '#4', 'first']);
  assert.deepEqual(cards.map(drawnName), ['#2', '#+2', 'first']);
  assert.equal(drawnName(project(circular).config.form.circular_record_selector), '#+2');
});

test('Load reads no file when no stored record ID reads as #N', async () => {
  const read = [];
  const { places, unresolved } = await readStoredRecordIdPlaces([
    { at: '', renderRequest: request('linear', [
      requestRecord('record-1', { kind: 'recordId', value: 'first' }),
      requestRecord('record-2', { kind: 'recordIndex', index: 1 }),
      requestRecord('record-3', null),
      // A `#` text that is no `#N` stays a record ID as it is.
      requestRecord('record-4', { kind: 'recordId', value: '#+ 2' })
    ]) },
    { at: '', renderRequest: null }
  ], (resourceId) => { read.push(resourceId); throw new Error('no read'); });
  assert.deepEqual({ read, places: places.size, unresolved }, { read: [], places: 0, unresolved: [] });
});

test('a stored record ID that no record or several records hold, or whose file cannot be read, has no place', async () => {
  const table = adoptCurrentSessionResources({
    named: descriptor(genbank, 'named.gb'), twice: descriptor(genbankOf([...NAMES, ...NAMES]), 'twice.gb')
  });
  const { places, unresolved } = await readStoredRecordIdPlaces([{ at: 'drawings[0].', renderRequest: request('linear', [
    requestRecord('record-1', { kind: 'recordId', value: '#9' }),
    requestRecord('record-2', { kind: 'recordId', value: '#2' }, 'twice'),
    requestRecord('record-3', { kind: 'recordId', value: '#2' }, 'missing')
  ]) }], (resourceId) => createSessionResourceFileView(table, resourceId));
  assert.equal(places.size, 0);
  assert.deepEqual(unresolved, [1, 2, 3].map((n) => ({
    path: `drawings[0].renderRequest.records[${n - 1}].selector.value`, recordId: n === 1 ? '#9' : '#2'
  })));
});

// The Load itself: a Session 51 a Python writer could write (no Web draft), whose Linear record is stored as the
// record ID `#2` in a file that holds records named `first`, `second`, `#2`, `#+2`.
const pythonWrittenSession = async (recordId) => {
  const session = await readSavedSession(new Blob([readFileSync(new URL('../../gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw', import.meta.url))]));
  session.resources['record-1-genbank'] = descriptor(genbank, 'record-1-genbank-named.gb');
  const [drawing] = session.drawings;
  for (const field of ['config', 'editorState', 'ui', 'webFiles', 'results', 'features', 'orthogroupState']) delete drawing[field];
  delete session.webFiles;
  drawing.renderRequest.records[0].selector = { kind: 'recordId', value: recordId };
  return session;
};
const load = async (session) => {
  const repairs = [];
  metrics.length = 0;
  const loaded = await importSession(
    // A Session 51 is a .gbdraw file; an older Session is its JSON (`sessionFileOf`).
    { target: { files: [new File([await sessionFileOf(session)], 'session')], value: 'selected' } },
    { transformLegacyResultSvg, confirmLoadRepairs: async (dialog) => { repairs.push(...dialog.items); return true; } });
  assert.ok(['ok', 'legacy'].includes(loaded.status), JSON.stringify(loaded.error));
  return { card: modeDrawing(state, 'linear').linearSeqs[0].region_record_id, decodes: decodes(), repairs };
};

test('Load of a Session gives the card the place of its stored record ID and decodes the file once', async () => {
  const plain = await load(await pythonWrittenSession('first'));
  assert.equal(plain.card, 'first');
  const named = await load(await pythonWrittenSession('#2'));
  assert.deepEqual(named, { card: '#3', decodes: plain.decodes, repairs: plain.repairs });
  const missing = await load(await pythonWrittenSession('#9'));
  assert.equal(missing.card, '#9');
  assert.ok(missing.repairs.some((item) => item.includes('"#9"') && item.includes('reads like a record number')), JSON.stringify(missing.repairs));
});

// A Session 31-39 stores its request at the top level; Load reads it the same way.
// A Session 31-39 stores its request at the top level; Load reads Python's Session 51 of it (with its draft, after
// reading the draft over the request's projection), and either path gives the card the record's place.
for (const keepsDraft of [false, true]) {
  test(`Load of a Session 31 ${keepsDraft ? 'with' : 'without'} its draft gives the card the place of its stored record ID`, async () => {
    for (const [recordId, place] of [['#2', '#3'], ['#+2', '#4']]) {
      const session = JSON.parse(gunzipSync(readFileSync(new URL('../fixtures/sessions/lambda_basic_linear.v31-schema1.gbdraw-session.json.gz', import.meta.url))));
      session.resources['record-1-genbank'] = descriptor(genbank, 'record-1-genbank-named.gb');
      for (const field of [...(keepsDraft ? [] : ['config']), 'results', 'orthogroupState']) delete session[field];
      session.renderRequest.records[0].selector = { kind: 'recordId', value: recordId };
      assert.equal((await load(session)).card, place, recordId);
    }
  });
}


// OV-3846: a Replace with the same records in the other order, saved and loaded before its Generate. The Load keeps
// no identities of the committed request, so its records are read from the request's own file, and the next Generate
// moves each record's rows to the record with its ID, as the Depth rows moved at the read after the Replace.
test('after a Load, the committed request\'s records are read from its own file where the card holds another', async () => {
  const resources = { then: descriptor(genbankOf(['first', 'second']), 'ab.gb') };
  const records = [requestRecord('card:1', { kind: 'recordIndex', index: 0 }, 'then'), requestRecord('card:2', { kind: 'recordIndex', index: 1 }, 'then')];
  const now = { mode: 'circular', cards: new Map([['card', buildDisambiguatedRecordEntries([{ recordId: 'second' }, { recordId: 'first' }])]]) };
  const moves = (sources) => Object.fromEntries(recordKeyMoves(drawnRecordIdentities(records, sources), drawnRecordIdentities(records, now)));
  // Read from the card's lists now, the records seem not to have changed: nothing moves (A's rows land on B).
  assert.deepEqual(moves(now), {});
  // The card holds a file the Replace picked; the store gives its bytes by ref.
  const holding = (file) => ({ cards: [{ uid: 'card', gb: { fileId: 'picked' } }], inputType: 'gb', inputFiles: { fileOf: (ref) => (ref ? file : null) } });
  const replaced = new File([genbankOf(['second', 'first'])], 'ba.gb', { type: 'text/plain' });
  const then = await committedRecordIdentitySources({ records, resources, sources: now, ...holding(replaced) });
  assert.deepEqual(then.cards.get('card').map(({ recordId }) => recordId), ['first', 'second']);
  assert.deepEqual(moves(then), { 'card:1': 'card:2', 'card:2': 'card:1' });
  // A card that still holds the request's file keeps its list now, with no read.
  const held = createSessionResourceFileView(adoptCurrentSessionResources(resources), 'then');
  const same = await committedRecordIdentitySources({ records, resources, sources: now, ...holding(held) });
  assert.equal(same.cards.get('card'), now.cards.get('card'));
  // A request file that cannot be read keeps the list now, and the failure is reported.
  const failures = () => metrics.filter(({ name }) => name === 'committedRecordReadFailureCount').map(({ card }) => card);
  const failed = failures().length;
  const unread = await committedRecordIdentitySources({ records: [requestRecord('card:1', null, 'missing')], resources, sources: now, ...holding(replaced) });
  assert.equal(unread.cards.get('card'), now.cards.get('card'));
  assert.deepEqual(failures().slice(failed), ['card']);
});

// OV-3846 review (MEDIUM 1): a Session saved after a Replace (first, first over the committed first, second: a
// repeated ID, so by place) and loaded, then a Replace with (second, first) before the next Generate. The Depth rows
// move at the read after that Replace and the feature edits at the Generate, both from the records the committed
// request drew, read from its own file, so first's Depth file and edits stay with first. Step by step, the rows of
// (first, first) would go by place.
test('after a Load, a second Replace moves the Depth rows where the next Generate moves the feature edits', async () => {
  const resources = { then: descriptor(genbankOf(['first', 'second']), 'ab.gb') };
  const records = [requestRecord('card:1', { kind: 'recordIndex', index: 0 }, 'then'), requestRecord('card:2', { kind: 'recordIndex', index: 1 }, 'then')];
  const now = { mode: 'circular', cards: new Map([['card', buildDisambiguatedRecordEntries([{ recordId: 'second' }, { recordId: 'first' }])]]) };
  const picked = new File([genbankOf(['second', 'first'])], 'ba.gb', { type: 'text/plain' });
  const committed = drawnRecordIdentities(records, await committedRecordIdentitySources({
    records, resources, sources: now, cards: [{ uid: 'card', gb: { fileId: 'picked' } }], inputType: 'gb', inputFiles: { fileOf: (ref) => (ref ? picked : null) }
  }));
  const depth = await circularDepthRowsAfterRead({
    matrix: [['a.tsv'], ['b.tsv']],
    // The cards as the Session was saved: its rows stand for (first, first).
    marks: [{ uid: 'card', primaryFile: 'aa', pairedFile: null, recordIds: ['first', 'first'] }],
    inputs: [{ uid: 'card', spec: null, primaryFile: 'ba', pairedFile: null }],
    records: [[{ record_id: 'second' }, { record_id: 'first' }]],
    fileOf: () => null,
    committed: async () => new Map(committed.map(({ card, cardRecordIds }) => [card, cardRecordIds]))
  });
  assert.deepEqual(depth, { rows: [['b.tsv'], ['a.tsv']] }, 'first\'s Depth file goes to first, at place 2');
  assert.deepEqual(Object.fromEntries(recordKeyMoves(committed, drawnRecordIdentities(records, now))), { 'card:1': 'card:2', 'card:2': 'card:1' },
    'first\'s feature edits go to place 2 too');
});

// OV-3846 review of 91340c36c (MEDIUM 2): a Session saved with (first, second) committed, loaded, and the card's file
// replaced by (second, first), when the committed request's file cannot be read. The Generate then reads the card's
// list now as the records it drew, so first's feature edits stay at place 1; the Depth rows cannot be placed either,
// so they stay where they are too (and Generate reports the matrix): both keep first's rows at place 1.
test('after a Load, a committed file that cannot be read keeps a Replace\'s Depth rows and edits where they are', async () => {
  const resources = { then: descriptor(genbankOf(['first', 'second']), 'ab.gb') };
  const records = [requestRecord('card:1', { kind: 'recordIndex', index: 0 }, 'missing'), requestRecord('card:2', { kind: 'recordIndex', index: 1 }, 'missing')];
  const now = { mode: 'circular', cards: new Map([['card', buildDisambiguatedRecordEntries([{ recordId: 'second' }, { recordId: 'first' }])]]) };
  const picked = new File([genbankOf(['second', 'first'])], 'ba.gb', { type: 'text/plain' });
  const sources = await committedRecordIdentitySources({
    records, resources, sources: now, cards: [{ uid: 'card', gb: { fileId: 'picked' } }], inputType: 'gb', inputFiles: { fileOf: (ref) => (ref ? picked : null) }
  });
  assert.deepEqual([...sources.unread], ['card']);
  const committed = drawnRecordIdentities(records, sources);
  assert.deepEqual(Object.fromEntries(recordKeyMoves(committed, drawnRecordIdentities(records, now))), {}, 'first\'s edits stay at place 1');
  const depth = await circularDepthRowsAfterRead({
    matrix: [['a.tsv'], ['b.tsv']],
    marks: [{ uid: 'card', primaryFile: 'ab', pairedFile: null, recordIds: ['first', 'second'] }],
    inputs: [{ uid: 'card', spec: null, primaryFile: 'ba', pairedFile: null }],
    records: [[{ record_id: 'second' }, { record_id: 'first' }]],
    fileOf: () => null,
    committed: async () => committedCardRecords({ identities: committed, unread: sources.unread })
  });
  assert.deepEqual(depth, { rows: null }, 'first\'s Depth file stays at place 1');
});
