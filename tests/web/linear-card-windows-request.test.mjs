// The request a Linear drawing with card windows writes (PD-OI-124), built by the app's own request owner
// (`projectDraftRequest`, which Save without a Result, New drawing from region and Generate share) and read by the
// real Python codec. Record Layout edits leave stored places on some cells: the request writes a place for every
// cell or for none, as Python reads them (OV-3071), and a card that draws its file whole takes a place as the file:
// Python numbers its records within their row (OV-3072). A drawing whose places can be written as they are writes them
// unchanged. A card that
// chose its record and lists a window reads its record from Generate's record catalog, so an empty window draws the
// whole record and a refusal names the record (B1).
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { evaluatePython } from './helpers/python-evaluator.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {
  createElement: () => ({ click: () => {} }), querySelector: () => null, documentElement: { setAttribute: () => {} }
};
installFakeSvgDom();
globalThis.File = File;

const { state, createDrawingState } = await import('../../gbdraw/web/js/state.js');
const { withNewDrawingDefaults } = await import('../../gbdraw/web/js/services/reset.js');
const { createLinearSeq } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { projectDraftRequest } = await import('../../gbdraw/web/js/services/config.js');
const { annotationSourceKey, buildAnnotationRecordCatalog } = await import('../../gbdraw/web/js/app/annotations/record-catalog.js');
const { createDrawingRecordCatalog } = await import('../../gbdraw/web/js/app/drawing-record-catalog.js');
const { discoverSequenceRecords } = await import('../../gbdraw/web/js/services/record-discovery.js');
const { createDiagramResourceTransport } = await import('../../gbdraw/web/js/services/diagram-resource-staging.js');
const { createLinearRecordEdits } = await import('../../gbdraw/web/js/app/linear-source-edits.js');
const { moveLinearLayoutCell, setLinearRecordRow } = await import('../../gbdraw/web/js/services/linear-record-layout.js');
const { normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
const { liveSettingInputsOf, linearRecordCardIndexes } = await import('../../gbdraw/web/js/services/session-request.js');

/** A linear GenBank record of 1,200 bp, or one file holding several. */
const genbankRecord = (name) => {
  const sequence = 'atgcgtacgt'.repeat(120);
  const lines = [];
  for (let offset = 0; offset < sequence.length; offset += 60) {
    lines.push(`${String(offset + 1).padStart(9)} ${sequence.slice(offset, offset + 60).match(/.{1,10}/g).join(' ')}`);
  }
  return [
    `LOCUS       ${name.padEnd(16)}       1200 bp    DNA     linear BCT 01-JAN-2000`,
    `DEFINITION  ${name} test record.`, `ACCESSION   ${name}`, `VERSION     ${name}.1`,
    'FEATURES             Location/Qualifiers', '     source          1..1200', '     CDS             101..400',
    '                     /locus_tag="x"', 'ORIGIN', ...lines, '//'
  ].join('\n');
};
const fileOf = (...names) => state.inputFiles.registerFile(
  new File([`${names.map(genbankRecord).join('\n')}\n`], `${names[0]}.gb`, { type: 'text/plain', lastModified: 1 }));
const recordsOf = (seq) => discoverSequenceRecords({
  file: state.inputFiles.fileOf(seq.gb), format: 'genbank', runHelperOperation: () => { throw new Error('read by the fast path'); }
});

/** Generate's record catalog of `drawing`, as `prepareLinearRecordCatalog` builds it (null when it builds none). */
const preparedCatalog = async (drawing) => {
  const read = new Map(await Promise.all(drawing.linearSeqs.map(async (seq) => [seq.uid, await recordsOf(seq)])));
  const { prepareLinearRecordCatalog } = createDrawingRecordCatalog({
    state: { activeDrawing: () => drawing }, readInputRecords: async () => { throw new Error('not read privately'); },
    linearRecordCatalogStateFor: (_drawing, seq) => ({ status: 'ready', records: read.get(seq.uid), error: '' }),
    linearRecordErrorModelFor: () => '', refreshLinearRecords: async () => ({})
  });
  const { catalog, error } = await prepareLinearRecordCatalog();
  assert.equal(error, '');
  return catalog;
};
/** The same catalog, built whatever the drawing holds. */
const fullCatalog = async (drawing) => {
  const off = new Set(drawing.recordsOff);
  return buildAnnotationRecordCatalog({ mode: 'linear', linearSources: await Promise.all(drawing.linearSeqs.map(async (seq) => ({
    sourceKey: annotationSourceKey({ scope: 'linear', uid: seq.uid, inputType: 'gb', primaryFile: seq.gb, pairedFile: null }),
    selector: seq.region_record_id, drawn: !off.has(seq.uid), hasInput: true, status: 'ready', records: await recordsOf(seq)
  }))) });
};

const stagedFiles = async (request, resources) => {
  const { resourceManifest, stagedResources } = await createDiagramResourceTransport().prepare({ request, resources });
  const bytesByToken = new Map(await Promise.all(stagedResources.map(async ({ cacheToken, blob }) => [cacheToken, Buffer.from(await blob.arrayBuffer()).toString('base64')])));
  return Object.fromEntries(resourceManifest.map(({ resourceId, cacheToken, name }) => [resourceId, { name, base64: bytesByToken.get(cacheToken) }]));
};
/** 'drawn' when the real Python codec reads `renderRequest` as a Linear request, else its refusal. */
const pythonVerdict = async ({ renderRequest, resources }) => JSON.parse(await evaluatePython(`
import base64, json, os, tempfile
from gbdraw.session_request_codec import decode_canonical_request
directory = tempfile.mkdtemp()
paths = {}
for index, (resource_id, staged) in enumerate(payload['staged'].items()):
    path = os.path.join(directory, f'r{index}_{staged["name"]}')
    with open(path, 'wb') as handle:
        handle.write(base64.b64decode(staged['base64']))
    paths[resource_id] = path
try:
    decode_canonical_request(payload['canonicalRequest'], resource_paths=paths, output_directory=directory)
    result = json.dumps('drawn')
except Exception as exc:
    result = json.dumps(str(exc))
`, { canonicalRequest: renderRequest, staged: await stagedFiles(renderRequest, resources) }));

/** The records Python plans for `renderRequest`, in request order, as `[record key, row, column]`. */
const pythonPlaces = async ({ renderRequest, resources }) => JSON.parse(await evaluatePython(`
import base64, json, os, tempfile
from gbdraw.api.request_render import plan_linear_request
from gbdraw.session_request_codec import decode_canonical_request
directory = tempfile.mkdtemp()
paths = {}
for index, (resource_id, staged) in enumerate(payload['staged'].items()):
    path = os.path.join(directory, f'r{index}_{staged["name"]}')
    with open(path, 'wb') as handle:
        handle.write(base64.b64decode(staged['base64']))
    paths[resource_id] = path
request = decode_canonical_request(payload['canonicalRequest'], resource_paths=paths, output_directory=directory)
result = json.dumps([[record.record_key, record.presentation.grid_row, record.presentation.grid_column]
                     for record in plan_linear_request(request).request.records])
`, { canonicalRequest: renderRequest, staged: await stagedFiles(renderRequest, resources) }));

/**
 * The request of `drawing` (each record as `[key, cardinality, row, column, windows as [start, end, column]]`) and
 * Python's verdict on it.
 */
const requestOf = async (drawing, catalog = undefined) => {
  const linearRecordCatalog = catalog === undefined ? await fullCatalog(drawing) : catalog;
  const { envelope } = await projectDraftRequest(state, drawing, { linearRecordCatalog, recordDisplayRows: [], omittedAnnotationRecordKeys: [] });
  return {
    records: envelope.renderRequest.records.map(({ recordKey, cardinality, presentation, windows }) => [
      recordKey, cardinality, presentation.gridRow, presentation.gridColumn, windows.map(({ start, end, column }) => [start, end, column])
    ]),
    verdict: await pythonVerdict(envelope)
  };
};

let drawings = 0;
/** A Linear drawing with Arrange in rows on, one card per uid, each reading its own one-record file. */
const arrangedDrawing = (uids) => {
  const drawing = withNewDrawingDefaults(state, createDrawingState(`windows-request-${drawings += 1}`, 'linear'));
  drawing.inputType.value = 'gb';
  drawing.linearSeqs.splice(0, Infinity, ...uids.map((uid) => createLinearSeq({ uid, gb: fileOf(`rec${uid}`) })));
  Object.assign(drawing.form, { labels_mode: 'none' });
  drawing.linearRecordLayoutEnabled.value = true;
  return drawing;
};
/** The app's Add window and Remove (the record rows' transition writes the cards it is given). */
const cardEditsOf = (drawing) => createLinearRecordEdits({
  sessionOperationAvailability: () => null, shownDrawingOfMode: () => drawing,
  runUndoable: (_label, action) => { action(); return true; }, clearAlignmentForMutation: () => {},
  linearRecordsFor: () => [], linearRecordStatusFor: () => ({}), rebindAnnotationRecords: () => {}, shownUnitPlaces: () => new Map(),
  applyLinearSeqMutation: (target, items) => { target.linearSeqs.splice(0, Infinity, ...items); }
});
/** Up or Down of a Record Layout cell, by the rule of `moveLinearRecordWithinRow` (app/linear-comparison-edits.js). */
const moveCell = (drawing, uid, direction, window = null) => {
  const next = moveLinearLayoutCell(drawing.linearSeqs, drawing.linearRecordRows, { uid, window }, direction);
  assert.ok(next);
  drawing.linearRecordRows.splice(0, Infinity, ...next);
};

test('a drawing whose stored places can be written as they are writes them unchanged', async () => {
  // As a CLI Session places them: every record has a place, each its own.
  const drawing = arrangedDrawing(['a', 'b', 'c']);
  const one = { canonicalCardinality: 'exactly_one' };
  drawing.linearRecordRows.splice(0, Infinity, { uid: 'a', row: 2, canonicalRow: 2, canonicalColumn: 1, ...one },
    { uid: 'b', row: 1, canonicalRow: 1, canonicalColumn: 1, ...one }, { uid: 'c', row: 1, canonicalRow: 1, canonicalColumn: 2, ...one });
  assert.deepEqual(await requestOf(drawing), { verdict: 'drawn', records: [
    ['a', 'exactly_one', 2, 1, []], ['b', 'exactly_one', 1, 1, []], ['c', 'exactly_one', 1, 2, []]] });
  // No cell has a place: none is written, and a card of a one-record file stays a whole file.
  const plain = arrangedDrawing(['a', 'b']);
  assert.deepEqual(await requestOf(plain), { verdict: 'drawn', records: [['a', 'all', 1, null, []], ['b', 'all', 2, null, []]] });
});

test('a Row change that drops one record\'s place writes a place for every record (OV-3071)', async () => {
  // b and c placed in row 1 and a in row 2, as a CLI Session places them; then a's Row is 1. dev wrote a: none,
  // b: 1, c: 2, and Python refused it ("If one Linear record has grid_column, every record must have it.").
  const drawing = arrangedDrawing(['a', 'b', 'c']);
  const one = { canonicalCardinality: 'exactly_one' };
  drawing.linearRecordRows.splice(0, Infinity, { uid: 'a', row: 2, canonicalRow: 2, canonicalColumn: 1, ...one },
    { uid: 'b', row: 1, canonicalRow: 1, canonicalColumn: 1, ...one }, { uid: 'c', row: 1, canonicalRow: 1, canonicalColumn: 2, ...one });
  setLinearRecordRow(drawing.linearRecordRows, 'a', 1);
  assert.deepEqual(await requestOf(drawing), { verdict: 'drawn', records: [
    ['a', 'exactly_one', 1, 3, []], ['b', 'exactly_one', 1, 1, []], ['c', 'exactly_one', 1, 2, []]] });
  // A card added to that row (a new file, drawn whole) takes the place after them as its file.
  drawing.linearSeqs.push(createLinearSeq({ uid: 'n', gb: fileOf('recn') }));
  drawing.linearRecordRows.push({ uid: 'n', row: 1 });
  const added = await requestOf(drawing);
  assert.equal(added.verdict, 'drawn');
  assert.deepEqual(added.records.at(-1), ['n', 'all', 1, 4, []]);
});

test('Add window and then Remove in an arranged drawing keep a request Python draws (review cases A and B, OV-3072)', async () => {
  // Case A: a and w (a one-window crop) in row 1, c in row 2. Add window, Remove the new window.
  const drawing = arrangedDrawing(['a', 'w', 'c']);
  Object.assign(drawing.linearSeqs[1], { region_start: 1, region_end: 100 });
  drawing.linearRecordRows.splice(0, Infinity, { uid: 'a', row: 1 }, { uid: 'w', row: 1 }, { uid: 'c', row: 2 });
  const edits = cardEditsOf(drawing);
  edits.addLinearCardWindow(drawing.linearSeqs[1]);
  edits.setLinearCardWindow(drawing.linearSeqs[1], 1, 'start', 301);
  edits.setLinearCardWindow(drawing.linearSeqs[1], 1, 'end', 400);
  // a and c draw their files whole: each takes its place as its file.
  assert.deepEqual(await requestOf(drawing), { verdict: 'drawn', records: [
    ['a', 'all', 1, 1, []], ['w', 'exactly_one', 1, null, [[1, 100, 2], [301, 400, 3]]], ['c', 'all', 2, 1, []]] });
  edits.removeLinearCardWindow(drawing.linearSeqs[1], 1);
  assert.deepEqual(await requestOf(drawing), { verdict: 'drawn', records: [
    ['a', 'all', 1, 1, []], ['w', 'exactly_one', 1, 2, [[1, 100, null]]], ['c', 'all', 2, 1, []]] });
  // Case B: a, b, w in row 1, w drawn in two windows. Down on a, then Remove one of w's windows.
  const shared = arrangedDrawing(['a', 'b', 'w']);
  shared.linearSeqs[2].windows = [1, 301].map((start) => ({ start, end: start + 99, reverseComplement: false, column: null, acrossOrigin: false }));
  shared.linearRecordRows.splice(0, Infinity, ...['a', 'b', 'w'].map((uid) => ({ uid, row: 1 })));
  moveCell(shared, 'a', 1);
  cardEditsOf(shared).removeLinearCardWindow(shared.linearSeqs[2], 1);
  const removed = await requestOf(shared);
  assert.equal(removed.verdict, 'drawn');
  assert.deepEqual(removed.records.map(([key, , row, column]) => [key, row, column]), [['b', 1, 1], ['a', 1, 2], ['w', 1, 3]]);
});

test('Add window beside a card that draws a file of three records whole keeps a request Python draws (OV-3072)', async () => {
  // m draws its file of three records whole; w reads another file. Both in row 1; Add window on w writes places.
  const drawing = arrangedDrawing(['w']);
  drawing.linearSeqs.unshift(createLinearSeq({ uid: 'm', gb: fileOf('recx', 'recy', 'recz') }));
  drawing.linearRecordRows.splice(0, Infinity, { uid: 'm', row: 1 }, { uid: 'w', row: 1 });
  const edits = cardEditsOf(drawing);
  edits.addLinearCardWindow(drawing.linearSeqs[1]);
  for (const [window, start, end] of [[0, 1, 100], [1, 301, 400]]) {
    edits.setLinearCardWindow(drawing.linearSeqs[1], window, 'start', start);
    edits.setLinearCardWindow(drawing.linearSeqs[1], window, 'end', end);
  }
  // m is the file at the first place of the row; Python numbers its records 1 to 3 and w's windows after them.
  // dev refused it: "grid_column cannot be set for a Linear RecordCardinality.ALL input".
  assert.deepEqual(await requestOf(drawing), { verdict: 'drawn', records: [
    ['m', 'all', 1, 1, []], ['w', 'exactly_one', 1, null, [[1, 100, 2], [301, 400, 3]]]] });
  const planned = async () => pythonPlaces((await projectDraftRequest(state, drawing, {
    linearRecordCatalog: await fullCatalog(drawing), recordDisplayRows: [], omittedAnnotationRecordKeys: []
  })).envelope);
  assert.deepEqual(await planned(), [['m:1', 1, 1], ['m:2', 1, 2], ['m:3', 1, 3], ['w~w1', 1, 4], ['w~w2', 1, 5]]);
  // Down on m keeps a request Python draws: the file between w's windows, the cards swapped as Down swaps them.
  moveCell(drawing, 'm', 1);
  assert.deepEqual(await requestOf(drawing), { verdict: 'drawn', records: [
    ['w', 'exactly_one', 1, null, [[1, 100, 1], [301, 400, 3]]], ['m', 'all', 1, 2, []]] });
  assert.deepEqual(await planned(), [['w~w1', 1, 1], ['w~w2', 1, 5], ['m:1', 1, 2], ['m:2', 1, 3], ['m:3', 1, 4]]);
});

test('a card drawn in two windows with Draw this record off leaves a request Python draws (review case C)', async () => {
  const drawing = arrangedDrawing(['a', 'b', 'c', 'off']);
  drawing.linearSeqs[3].windows = [1, 301].map((start) => ({ start, end: start + 99, reverseComplement: false, column: null, acrossOrigin: false }));
  drawing.linearRecordRows.splice(0, Infinity, ...['a', 'b', 'c', 'off'].map((uid) => ({ uid, row: 1 })));
  drawing.recordsOff.splice(0, Infinity, 'off');
  moveCell(drawing, 'a', 1);
  drawing.linearSeqs.push(createLinearSeq({ uid: 'n', gb: fileOf('recn') }));
  drawing.linearRecordRows.push({ uid: 'n', row: 1 });
  const written = await requestOf(drawing);
  assert.equal(written.verdict, 'drawn');
  assert.deepEqual(written.records.map(([key, , row, column]) => [key, row, column]), [['b', 1, 1], ['a', 1, 2], ['c', 1, 3], ['n', 1, 4]]);
});

test('a card that chose its record and lists an empty window draws the whole record, without another reader (B1)', async () => {
  // One file of two records; the card chose the second; no region annotation and no card without a chosen record.
  const drawing = arrangedDrawing([]);
  drawing.linearRecordLayoutEnabled.value = false;
  drawing.linearSeqs.splice(0, Infinity, createLinearSeq({ uid: 'm', gb: fileOf('recx', 'recy'), region_record_id: 'recy.1' }));
  // Without a window Generate reads the card's records too: the selector it writes depends on its file's records.
  assert.ok(await preparedCatalog(drawing));
  // Add window on the Region block: the block's own window and the new one, both with Start and End empty.
  cardEditsOf(drawing).addLinearCardWindow(drawing.linearSeqs[0]);
  const catalog = await preparedCatalog(drawing);
  assert.ok(catalog);
  assert.deepEqual(await requestOf(drawing, catalog), { verdict: 'drawn', records: [['m', 'exactly_one', null, null, [[1, 1200, null], [1, 1200, null]]]] });
  // A window across the origin whose Start is not after its End is refused, naming the record.
  Object.assign(drawing.linearSeqs[0].windows[0], { start: 100, end: 900, acrossOrigin: true });
  await assert.rejects(requestOf(drawing, await preparedCatalog(drawing)), (error) => {
    assert.match(normalizeUserFacingError(error, { stage: 'request-validation' }).summary, /Record: recy\.1\./);
    return true;
  });
});

/** The request of `drawing` built with `catalog` (Generate's record catalog), as Generate writes it. */
const envelopeOf = async (drawing, catalog) => (await projectDraftRequest(state, drawing, {
  linearRecordCatalog: catalog, recordDisplayRows: [], omittedAnnotationRecordKeys: []
})).envelope;

test('a card that names the record of a one-record file writes the CLI\'s selector, none, with or without windows', async () => {
  // As the CLI writes it (`_encode_records` in gbdraw/session_request_codec.py): a record of a file that holds one
  // record has no selector; a record of a file of several is named. One rule with Circular (`cliRecordSelector`).
  const drawing = arrangedDrawing([]);
  drawing.linearRecordLayoutEnabled.value = false;
  drawing.linearSeqs.splice(0, Infinity, createLinearSeq({ uid: 'one', gb: fileOf('recone'), region_record_id: 'recone.1' }),
    createLinearSeq({ uid: 'two', gb: fileOf('recx', 'recy'), region_record_id: 'recy.1' }));
  const selectorsOf = async () => {
    const envelope = await envelopeOf(drawing, await preparedCatalog(drawing));
    assert.equal(await pythonVerdict(envelope), 'drawn');
    return envelope.renderRequest.records.map(({ cardinality, selector }) => [cardinality, selector]);
  };
  const expected = [['exactly_one', null], ['exactly_one', { kind: 'recordId', value: 'recy.1' }]];
  assert.deepEqual(await selectorsOf(), expected);
  // Drawn in windows, the record keeps the same selector (the CLI's `--region` form).
  cardEditsOf(drawing).addLinearCardWindow(drawing.linearSeqs[0]);
  drawing.linearSeqs[0].windows.splice(0, Infinity, ...[[1, 100], [201, 300]].map(([start, end]) => (
    { start, end, reverseComplement: false, column: null, acrossOrigin: false })));
  assert.deepEqual(await selectorsOf(), expected);
});

test('a card still draws a committed record of a one-record file, named (a Session saved before) or not', async () => {
  // The live preview applies record names while the cards draw the committed records (`linearRecordCardIndexes`).
  const drawing = arrangedDrawing([]);
  drawing.linearRecordLayoutEnabled.value = false;
  drawing.linearSeqs.splice(0, Infinity, createLinearSeq({ uid: 'one', gb: fileOf('recone'), region_record_id: 'recone.1' }),
    createLinearSeq({ uid: 'two', gb: fileOf('recx', 'recy'), region_record_id: 'recy.1' }));
  const envelope = await envelopeOf(drawing, await preparedCatalog(drawing));
  const inputs = liveSettingInputsOf(drawing, (ref) => state.inputFiles.fileOf(ref));
  // The CLI's form (and Generate's now): the record of the one-record file has no selector.
  const cli = structuredClone(envelope.renderRequest);
  cli.records[0].selector = null;
  assert.deepEqual(linearRecordCardIndexes(cli, inputs, envelope.resources), [0, 1]);
  const named = structuredClone(envelope.renderRequest);
  named.records[0].selector = { kind: 'recordId', value: 'recone.1' };
  assert.deepEqual(linearRecordCardIndexes(named, inputs, envelope.resources), [0, 1]);
  // A card that names another record of a file of several does not draw the committed one.
  const other = structuredClone(envelope.renderRequest);
  other.records[1].selector = { kind: 'recordId', value: 'recx.1' };
  assert.equal(linearRecordCardIndexes(other, inputs, envelope.resources), null);
});
