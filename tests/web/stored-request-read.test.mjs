// The one reader of a stored render request at Load (`readStoredRenderRequest`
// in services/session-request.js, PR-C step C2c): the envelope, the schema
// 5-8 promotion with the Load context, the repair of a request against its
// description with the request policy (references, ring columns, settings of
// the other diagram type), the root (K22), and its Load wiring for the
// top-level request and `otherModeResult` (OV-670).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { File } from 'node:buffer';
import { SESSION_VERSION } from '../../gbdraw/web/js/contracts/schema.generated.js';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { evaluatePython, validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';
import { callerSession } from './helpers/session-caller.mjs';

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

const REPO = new URL('../../', import.meta.url);
const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { setUnmanagedConfigOverrideValidator } = sessionConfig;
const {
  readStoredRenderRequest, readPublishedRenderRequest, projectCommittedEditorIntent
} = await import('../../gbdraw/web/js/services/session-request.js');
const { describeLoadRepairs } = await import('../../gbdraw/web/js/utils/load-repairs.js');
setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { importSession } = callerSession(sessionConfig, state);

const read = (path) => {
  const raw = readFileSync(new URL(path, REPO));
  return JSON.parse((raw[0] === 0x1f ? gunzipSync(raw) : raw).toString('utf8'));
};
const gallery = (name) => read(`gbdraw/web/gallery/sessions/${name}.gbdraw-session.json`);
const clone = (value) => JSON.parse(JSON.stringify(value));
const contextOf = (session, path = 'renderRequest') => ({
  featureCatalog: session.editorState?.featureCatalog ?? null,
  legacyOrthogroupState: session.orthogroupState ?? null,
  resourceIds: new Set(Object.keys(session.resources || {})),
  path
});
const readRequest = (session, change = () => {}) => {
  const request = clone(session.renderRequest);
  change(request);
  return readStoredRenderRequest(request, contextOf(session));
};
const rows = (result) => result.repairs.map(({ code, path, used }) => ({ code, path, ...(used === undefined ? {} : { used }) }));

const input = (session) => ({
  target: { files: [new File([JSON.stringify(session)], 'session.json', { type: 'application/json' })], value: 'selected' }
});
/** Load with [Load] answered; the dialog lines and the result. */
const load = async (session) => {
  const dialogs = [];
  const result = await importSession(input(session), {
    transformLegacyResultSvg, confirmLoadRepairs: async (dialog) => { dialogs.push(dialog.items); return true; }
  });
  return { result, dialogs };
};

test('a request schema the app does not read is a Session-format error (OV-773); a schema that is no integer leaves the request unread (OV-863)', () => {
  for (const schema of [3, 4, 10]) {
    assert.throws(
      () => readRequest(gallery('lambda_basic_linear'), (request) => { request.schema = schema; }),
      (error) => error.code === 'INPUT_INVALID' && error.context?.reason === 'SESSION_FORMAT',
      String(schema)
    );
  }
  for (const schema of [undefined, null, 'nine', '9', 9.5]) {
    assert.deepEqual(readRequest(gallery('lambda_basic_linear'), (request) => { request.schema = schema; }),
      { kind: 'request-dropped', path: 'renderRequest.schema', reason: 'wrong-kind' }, String(schema));
  }
});

test('a request that needs no repair comes back as the stored object (RK-3); schema 1 and 2 are lifted first', () => {
  const session = gallery('lambda_basic_linear');
  const result = readStoredRenderRequest(session.renderRequest, contextOf(session));
  assert.equal(result.value, session.renderRequest);
  assert.deepEqual(result.repairs, []);
  for (const schema of [1, 2]) {
    // A schema 1 or 2 request the lift cannot read leaves nothing to read (K22).
    assert.deepEqual(readStoredRenderRequest({ schema, mode: 'linear', records: [] }, contextOf(session)),
      { kind: 'request-dropped', path: 'renderRequest', reason: 'wrong-kind' });
  }
});

test('a saved Collinear anchor mode that was removed reads as rbh with a row; "either" reads without one (OV-2303, OV-2301)', () => {
  const session = gallery('BGC0000708-BGC0000713');
  const index = session.renderRequest.comparisons.findIndex((item) => item.settings);
  assert.equal(session.renderRequest.comparisons[index].settings.collinearityParams.parameters.mergeOrientation, 'either');
  assert.deepEqual(rows(readRequest(session)), []);
  const result = readRequest(session, (request) => { request.comparisons[index].settings.collinearityAnchorMode = 'one_to_one'; });
  assert.deepEqual(rows(result), [{
    code: 'setting-replaced', path: `renderRequest.comparisons[${index}].settings.collinearityAnchorMode`, used: 'rbh'
  }]);
  assert.equal(result.value.comparisons[index].settings.collinearityAnchorMode, 'rbh');
});

test('a request of schema 5 to 8 is promoted once, with the Load context, without a row (OV-770)', () => {
  for (const path of ['tests/fixtures/sessions/lambda_basic_linear.v40-schema5.json',
    'tests/fixtures/sessions/BGC0000708-BGC0000713.v40-schema5.json',
    'tests/fixtures/sessions/HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz']) {
    const session = read(path);
    const result = readStoredRenderRequest(session.renderRequest, contextOf(session));
    assert.equal(result.value.schema, 9, path);
    assert.deepEqual(result.repairs, [], path);
  }
  // The legacy alignment of BGC schema 5 is read with the Session's catalog (BC-1).
  const bgc = read('tests/fixtures/sessions/BGC0000708-BGC0000713.v40-schema5.json');
  const alignment = readStoredRenderRequest(bgc.renderRequest, contextOf(bgc)).value.layout.similarityAlignment;
  assert.equal(alignment.reference.sourceFeatureIndex, 14);
  // Without the catalog it cannot be read: the alignment is left out with one row.
  const blind = readStoredRenderRequest(bgc.renderRequest, { ...contextOf(bgc), featureCatalog: null, legacyOrthogroupState: null });
  assert.equal(blind.value.layout.similarityAlignment ?? null, null);
  assert.deepEqual(rows(blind), [{ code: 'request-part-dropped', path: 'renderRequest.layout.similarityAlignment' }]);
});

test('an annotation colour the decoder refuses takes its default; a name in any case is kept as stored (C-OD11, C-OD25)', () => {
  const session = gallery('tobacco-chloroplast');
  const result = readRequest(session, (request) => {
    const style = request.diagramOptions.annotations.sets[0].annotations[0].style;
    Object.assign(style, { stroke: 'valueTheDescriptionDoesNotName', fill: '#abcd', labelColor: ' Navy ' });
  });
  const style = result.value.diagramOptions.annotations.sets[0].annotations[0].style;
  assert.deepEqual([style.stroke, style.fill, style.labelColor], ['#404040', null, ' Navy ']);
  const at = 'renderRequest.diagramOptions.annotations.sets[0].annotations[0].style';
  assert.deepEqual(rows(result), [
    { code: 'setting-replaced', path: `${at}.stroke`, used: '#404040' },
    { code: 'setting-replaced', path: `${at}.fill` }
  ]);
});

test('a wrong value in a required field that takes null is read as null, not as its whole part (OV-777)', () => {
  const session = gallery('lambda_basic_linear');
  const table = { resourceId: Object.keys(session.resources)[0], representation: 'canonicalTsv' };
  const result = readRequest(session, (request) => { request.diagramOptions.colors.defaultColorsFile = table; });
  const expected = clone(session.renderRequest.diagramOptions);
  expected.colors.defaultColorsFile = null;
  assert.deepEqual(result.value.diagramOptions, expected);
  assert.deepEqual(rows(result), [{ code: 'request-part-dropped', path: 'renderRequest.diagramOptions.colors.defaultColorsFile' }]);
  // The wording follows the effect: a part replaced by its default says so (N1).
  assert.deepEqual(describeLoadRepairs(result.repairs).items, [
    'Settings: a part (renderRequest.diagramOptions.colors.defaultColorsFile) could not be read; its default is used.'
  ]);
});

test('a Linear layout that omits fields with a default is read as a Linear layout, not as the empty one', () => {
  const result = readRequest(gallery('lambda_basic_linear'), (request) => { request.layout = { recordGapPx: 10 }; });
  assert.deepEqual(result.value.layout, { recordGapPx: 10, multiRecordPositions: null, recordTranslations: [], similarityAlignment: null });
  // Each missing field is listed with its default, as anywhere in a request.
  assert.deepEqual(rows(result), ['multiRecordPositions', 'recordTranslations', 'similarityAlignment']
    .map((key) => ({ code: 'setting-replaced', path: `renderRequest.layout.${key}` })));
});

test('the reader turns only the layout owner\'s refusal into a dropped request; a programming error there is thrown (N2)', () => {
  const session = gallery('lambda_basic_linear');
  const request = clone(session.renderRequest);
  request.layout = { recordGapPx: 10 };
  const [record] = request.records;
  const key = record.recordKey;
  // A defect inside the layout owner, reached through the one read only it makes.
  Object.defineProperty(record, 'recordKey', {
    enumerable: true,
    get: () => {
      if (new Error().stack.includes('canonicalLinearRecordLayout')) throw new TypeError('defect in the layout owner');
      return key;
    }
  });
  assert.throws(() => readStoredRenderRequest(request, contextOf(session)), TypeError);
});

test('a comparison, pair or track axis that names nothing the Session holds is left out with a row (B14, OV-772)', () => {
  const session = gallery('BGC0000708-BGC0000713');
  const result = readRequest(session, (request) => {
    request.comparisons[0].resourceId = 'no-such-resource';
    request.comparisons[1].queryRecordIndex = 99;
  });
  assert.equal(result.value.comparisons.length, session.renderRequest.comparisons.length - 2);
  assert.deepEqual(rows(result), [
    { code: 'request-reference-dropped', path: 'renderRequest.comparisons[0]' },
    { code: 'request-reference-dropped', path: 'renderRequest.comparisons[1]' }
  ]);
  // An axis index takes null, its default, rather than dropping the tracks (OV-777).
  const lambda = gallery('lambda_basic_linear');
  const axis = readRequest(lambda, (request) => { request.diagramOptions.tracks.linearTrackAxisIndex = 3; });
  assert.deepEqual(axis.value.diagramOptions.tracks, lambda.renderRequest.diagramOptions.tracks);
  assert.deepEqual(rows(axis), [{ code: 'request-reference-dropped', path: 'renderRequest.diagramOptions.tracks.linearTrackAxisIndex' }]);
});

test('an unreadable conservation ring drops its column from every ring list, so no ring shifts (B16)', () => {
  const session = gallery('HmmtDNA_basic_circular');
  const ref = (resourceId) => ({ resourceId, representation: 'file' });
  const result = readStoredRenderRequest({
    ...clone(session.renderRequest),
    diagramOptions: {
      ...clone(session.renderRequest.diagramOptions),
      conservationBlastFiles: [ref('ring-a'), null, ref('ring-c')],
      conservationLabels: ['A', 'B', 'C'],
      conservationColors: ['#111111', '#222222', '#333333'],
      conservationFastaFiles: [null, null, null]
    }
  }, { ...contextOf(session), resourceIds: new Set([...Object.keys(session.resources), 'ring-a', 'ring-c']) });
  const options = result.value.diagramOptions;
  assert.deepEqual(options.conservationBlastFiles, [ref('ring-a'), ref('ring-c')]);
  assert.deepEqual(options.conservationLabels, ['A', 'C']);
  assert.deepEqual(options.conservationColors, ['#111111', '#333333']);
  assert.deepEqual(options.conservationFastaFiles, [null, null]);
  assert.deepEqual(rows(result), [{ code: 'setting-entries-dropped', path: 'renderRequest.diagramOptions.conservationBlastFiles' }]);
});

test('a managed conservation slot keeps its ring when a ring column is dropped (B16)', () => {
  // Python binds a slot to the ring whose index is its source_index: the slot
  // of the dropped ring goes, a later one is renumbered (its stored type kept),
  // and the axis boundary follows the slots before it.
  const session = gallery('HmmtDNA_basic_circular');
  const ref = (resourceId) => ({ resourceId, representation: 'file' });
  const slot = (id, renderer, params = {}) => ({
    kind: 'circularTrackSlot', id, renderer, enabled: true, side: 'inside', radius: null, width: null, z: 0,
    params, innerGapPx: null, outerGapPx: null
  });
  const ring = (label, sourceIndex, trackIndex) => slot(`ring-${label}`, 'sequence_conservation', {
    managed: 'conservation', source_index: sourceIndex, track_index: trackIndex, label, color: '#000000'
  });
  const request = clone(session.renderRequest);
  Object.assign(request.diagramOptions, {
    conservationBlastFiles: [ref('ring-a'), null, ref('ring-c')],
    conservationLabels: ['A', 'B', 'C'],
    conservationColors: ['#111111', '#222222', '#333333'],
    conservationFastaFiles: [null, null, null]
  });
  request.diagramOptions.tracks = {
    ...request.diagramOptions.tracks,
    circularTrackSlots: [slot('features', 'features'), ring('A', '0', '1'), ring('B', '1', '2'), ring('C', 2, 3), slot('ticks', 'ticks')],
    circularTrackAxisIndex: 4
  };
  const stored = clone(request);
  const result = readStoredRenderRequest(request, {
    ...contextOf(session), resourceIds: new Set([...Object.keys(session.resources), 'ring-a', 'ring-c'])
  });
  const { tracks } = result.value.diagramOptions;
  assert.deepEqual(tracks.circularTrackSlots, [
    slot('features', 'features'), ring('A', '0', '1'), ring('C', 1, 2), slot('ticks', 'ticks')
  ]);
  assert.equal(tracks.circularTrackAxisIndex, 3);
  assert.deepEqual(result.value.diagramOptions.conservationLabels, ['A', 'C']);
  assert.deepEqual(rows(result), [
    { code: 'setting-entries-dropped', path: 'renderRequest.diagramOptions.conservationBlastFiles' },
    { code: 'setting-entries-dropped', path: 'renderRequest.diagramOptions.tracks.circularTrackSlots' }
  ]);
  assert.deepEqual(request, stored, 'the stored request is not changed');
});

test('a setting of the other diagram type is read at its default without a row (〔出さないもの〕)', () => {
  const result = readRequest(gallery('lambda_basic_linear'), (request) => {
    request.diagramOptions.species = 'Escherichia phage';
    request.diagramOptions.tracks.circularTrackSlots = ['features'];
  });
  assert.equal(result.value.diagramOptions.species, null);
  assert.equal(result.value.diagramOptions.tracks.circularTrackSlots, null);
  assert.deepEqual(result.repairs, []);
});

test('a request that cannot be read at its root loads the Session as settings only, listed (K22, C-OD07)', async () => {
  for (const [name, change] of [
    ['lambda_basic_linear', (request) => { request.records = 'x'; }],
    ['HmmtDNA_basic_circular', (request) => { request.records = []; }],
    ['HmmtDNA_basic_circular', (request) => { request.mode = 'radial'; }]
  ]) {
    const session = gallery(name);
    change(session.renderRequest);
    const { result, dialogs } = await load(session);
    assert.equal(result.status, 'ok', `${name}: ${JSON.stringify(result.error)}`);
    assert.equal(state.results.value.length, 0, name);
    assert.equal(dialogs.length, 1, name);
    assert.match(dialogs[0][0], /^Saved diagram: the settings of the saved diagram \(renderRequest(\.mode)?\) could not be read\. The file loads as settings only/);
    assert.ok(result.repairs.some((repair) => repair.code === 'request-dropped'), name);
  }
});

test('a Session 42 to 47 with no stored request loads as settings only, listed; no writer omits the key (OV-785, C-OD07)', async () => {
  const session = gallery('HmmtDNA_basic_circular');
  delete session.renderRequest;
  const { result, dialogs } = await load(session);
  assert.equal(result.status, 'ok', JSON.stringify(result.error));
  assert.equal(state.results.value.length, 0);
  assert.deepEqual(result.repairs.map(({ code, path }) => [code, path]), [['request-dropped', 'renderRequest']]);
  assert.equal(dialogs.length, 1);
  assert.match(dialogs[0][0], /^Saved diagram: the settings of the saved diagram \(renderRequest\) could not be read\. The file loads as settings only/);
});

test('a Linear layout whose record translations do not cover the drawn records drops the request at its root (vector 027, C-OD20)', async () => {
  const session = gallery('BGC0000708-BGC0000713');
  const change = (request) => { request.layout.recordTranslations[0].recordKey = 'no-such-record'; };
  assert.deepEqual(readRequest(session, change), { kind: 'request-dropped', path: 'renderRequest.layout', reason: 'unknown-record' });
  change(session.renderRequest);
  const { result, dialogs } = await load(session);
  assert.equal(result.status, 'ok', JSON.stringify(result.error));
  assert.equal(state.results.value.length, 0);
  assert.deepEqual(result.repairs.map(({ code, path }) => [code, path]), [['request-dropped', 'renderRequest.layout']]);
  assert.equal(dialogs.length, 1);
  assert.match(dialogs[0][0], /^Saved diagram: the settings of the saved diagram \(renderRequest\.layout\) could not be read\. The file loads as settings only/);
});

test('the other diagram\'s request is read at Load: its repairs are listed, and an unreadable one leaves only that set out (OV-670)', async () => {
  const twoMode = () => evaluatePython('from tests.utils.two_mode_session import two_mode_session\nresult = two_mode_session()');
  const listed = await twoMode();
  // Unknown keys are counted under their object, one row per object (S7).
  Object.assign(listed.otherModeResult.renderRequest.records[0], { bogus: true, extra: 1 });
  const first = await load(listed);
  assert.equal(first.result.status, 'ok', JSON.stringify(first.result.error));
  assert.deepEqual(first.dialogs.flat(), [
    'Settings: 2 entries of otherModeResult.renderRequest.records[0] could not be read and were dropped.'
  ]);
  const dropped = await twoMode();
  dropped.otherModeResult.renderRequest.records = [];
  const second = await load(dropped);
  assert.equal(second.result.status, 'ok', JSON.stringify(second.result.error));
  assert.deepEqual(second.dialogs.flat(), [
    'Saved Linear diagram: its settings (otherModeResult.renderRequest) could not be read; that diagram and its Results are left out.'
  ]);
  assert.ok(state.results.value.length > 0);
});

const sessionConfig2 = await import('../../gbdraw/web/js/services/config.js');
const { getCommittedCanonicalSession } = sessionConfig2;
const { exportSession } = callerSession(sessionConfig2, state);
const frame = () => {
  const intent = modeDrawing(state, 'linear').importedComparisonIntent;
  return [intent.disposition, intent.message];
};

test('a request whose comparisons Load holds is promoted like any other, and the Comparisons notice is as before (C-OD14)', async () => {
  // The frame each input showed before the reader (base b62b9ceb), against
  // what Load now commits: schema 9, no lazy promotion left after Load.
  const missing = 'The saved comparison is missing a required resource.';
  const cases = [
    ['q-frame fixture', () => read('tests/fixtures/sessions/q-frame-main-linear-reverse.v42.gbdraw-session.json.gz'),
      'PRESERVED_READ_ONLY', []],
    ['se06 fixture', () => read('tests/fixtures/sessions/se06-main-linear-blast-cli.v42.gbdraw-session.json.gz'),
      'PRESERVED_READ_ONLY', []],
    // A held comparison's reference to a resource the Session lacks is dropped
    // with its row (OV-899): no Save could write it; the notice is as before.
    ['dangling comparison', () => {
      const session = gallery('BGC0000708-BGC0000713');
      session.renderRequest.comparisons[0].resourceId = 'missing';
      return session;
    }, 'DECISION_REQUIRED', ['request-reference-dropped'], missing],
    ['unknown comparison kind', () => {
      const session = gallery('BGC0000708-BGC0000713');
      session.renderRequest.comparisons[0].kind = 'other';
      return session;
    }, 'DECISION_REQUIRED', ['setting-entries-dropped'], 'The saved comparison kind \'other\' is unsupported.'],
    ['Linear comparison in a Circular request', () => {
      const session = gallery('HmmtDNA_basic_circular');
      session.renderRequest.comparisons = clone(gallery('BGC0000708-BGC0000713').renderRequest.comparisons.slice(0, 1));
      return session;
    }, 'DECISION_REQUIRED', [], 'A Circular session contains unsupported Linear comparison data.']
  ];
  for (const [name, make, disposition, codes, message] of cases) {
    const { result } = await load(make());
    assert.equal(result.status, 'ok', `${name}: ${JSON.stringify(result.error)}`);
    assert.deepEqual(result.repairs.map(({ code }) => code), codes, name);
    assert.equal(frame()[0], disposition, name);
    if (message) assert.equal(frame()[1], message, name);
    assert.equal(getCommittedCanonicalSession().renderRequest.schema, 9, name);
  }
});

test('Save writes the request Load read: the legacy alignment keeps the feature the Load context names (OV-770, BC-1)', async () => {
  const { result } = await load(read('tests/fixtures/sessions/BGC0000708-BGC0000713.v40-schema5.json'));
  assert.equal(result.status, 'ok', JSON.stringify(result.error));
  const saved = await exportSession('probe', { askChoice: async () => 'download' });
  assert.equal(saved.status, 'saved', JSON.stringify(saved.error));
  const session = JSON.parse(gunzipSync(Buffer.from(await saved.blob.arrayBuffer())).toString('utf8'));
  assert.equal(session.renderRequest.schema, 9);
  assert.equal(session.renderRequest.layout.similarityAlignment.reference.sourceFeatureIndex, 14);
});

test('an alignment plan that loses a record decision is left out as one part (vector 002)', () => {
  const session = gallery('BGC0000708-BGC0000713');
  const result = readRequest(session, (request) => { request.layout.similarityAlignment.records[0].recordKey = 'no-such-record'; });
  assert.equal(result.value.layout.similarityAlignment, null);
  assert.deepEqual(rows(result), [{ code: 'request-reference-dropped', path: 'renderRequest.layout.similarityAlignment' }]);
});

test('Gallery publication reads a stored request through the reader and refuses any repair (OV-771, PLAN 12.2)', () => {
  const session = read('tests/fixtures/sessions/lambda_basic_linear.v40-schema5.json');
  const context = { featureCatalog: null, legacyOrthogroupState: null, resourceIds: new Set(Object.keys(session.resources)) };
  // A schema-5 Linear record without a selector reads all records, as Load reads it.
  const published = readPublishedRenderRequest(session.renderRequest, context);
  assert.equal(published.schema, 9);
  assert.deepEqual(published.records.map(({ cardinality }) => cardinality), ['all']);
  const repaired = clone(session.renderRequest);
  repaired.records[0].bogus = true;
  assert.throws(() => readPublishedRenderRequest(repaired, context),
    (error) => error.code === 'INPUT_INVALID' && error.context?.reason === 'FIELDS');
});

test('a saved comparison entry that is neither resource-backed nor a generated comparison is dropped with a row (OV-750)', async () => {
  for (const value of [[null], ['x'], [{ kind: 'generatedProteinComparison' }]]) {
    const session = read('tests/fixtures/sessions/cli-linear-protein.v30.gbdraw-session.json.gz');
    session.files = { ...(session.files || {}), linearCanonicalComparisons: value };
    const { result, dialogs } = await load(session);
    assert.equal(result.status, 'ok', JSON.stringify(result.error));
    assert.deepEqual(state.files.linearCanonicalComparisons, [], JSON.stringify(value));
    assert.ok(result.repairs.some(({ code, path }) => code === 'setting-entries-dropped'
      && path === 'files.linearCanonicalComparisons'), JSON.stringify(value));
    assert.ok(dialogs.flat().some((line) => line.includes('files.linearCanonicalComparisons')));
  }
});

// G3 on the request path: every stored Session the app ships or pins (the
// Gallery and the Session fixtures) loads with no row about its stored
// request and commits it at the current schema (schemas 1 and 2 through the
// lift, PR-D). A row here is a defect of the input's producer.
/** Inputs this Node harness cannot load, by name, with the reason. */
const fixture = (name) => read(`tests/fixtures/sessions/${name}`);

test('every part the reader lists can be dropped alone: a null presentation, translation list or feature row list takes the writer\'s empty form (S1, S4, C-OD25)', () => {
  const lambda = gallery('lambda_basic_linear');
  const presentation = readRequest(lambda, (request) => { request.records[0].presentation = null; });
  assert.deepEqual(presentation.value.records[0].presentation,
    { label: null, subtitle: null, reverseComplement: false, gridRow: null, gridColumn: null });
  assert.deepEqual(rows(presentation), [{ code: 'setting-replaced', path: 'renderRequest.records[0].presentation' }]);
  const translations = readRequest(lambda, (request) => {
    request.layout = { recordGapPx: 12, multiRecordPositions: null, recordTranslations: null, similarityAlignment: null };
  });
  assert.deepEqual(translations.value.layout,
    { recordGapPx: 12, multiRecordPositions: null, recordTranslations: [], similarityAlignment: null });
  assert.deepEqual(rows(translations), [{ code: 'setting-replaced', path: 'renderRequest.layout.recordTranslations' }]);
  for (const key of ['featureOverrides', 'featurePlacements']) {
    const result = readRequest(lambda, (request) => { request.diagramOptions[key] = null; });
    assert.deepEqual(result.value.diagramOptions, { ...lambda.renderRequest.diagramOptions, [key]: [] }, key);
    assert.deepEqual(rows(result), [{ code: 'setting-replaced', path: `renderRequest.diagramOptions.${key}` }], key);
  }
});

test('the schema 5 to 7 promotion leaves a record that is not an object to the repair, never a TypeError (S5)', () => {
  const session = fixture('lambda_basic_linear.v40-schema5.json');
  for (const layout of [undefined, { recordGapPx: 10 }]) {
    const result = readRequest(session, (request) => {
      request.records.push(null);
      if (layout) request.layout = layout;
    });
    assert.equal(result.kind, 'request-dropped', JSON.stringify(layout));
  }
});

test('a CLI or older Session whose request only loses a part loads with its Results, as before the reader (F2, F3)', async () => {
  for (const name of ['cli-linear-tables.v41', 'cli-linear-tables.v42', 'se08-main-linear-cli.v42',
    'cli-linear-tables.v44', 'scale-interval-zero-circular-cli.v44', 'precomputed-ring-unlabelled-cli.v44']) {
    const session = fixture(`${name}.gbdraw-session.json.gz`);
    session.renderRequest.records[0].presentation = null;
    const { result } = await load(session);
    assert.equal(result.status, 'ok', `${name}: ${JSON.stringify(result.error)}`);
    assert.ok(state.results.value.length > 0, name);
    assert.equal(getCommittedCanonicalSession().renderRequest.schema, 9, name);
    const saved = await exportSession('probe', { askChoice: async () => 'download' });
    assert.equal(saved.status, 'saved', `${name}: ${JSON.stringify(saved.error)}`);
  }
});

test('a request that cannot be read at its root in a Session that keeps no settings without it is refused with a typed error, never committed unread (F2, F3, C-OD25)', async () => {
  for (const name of ['cli-linear-tables.v41', 'cli-linear-tables.v42', 'scale-interval-zero-circular-cli.v44']) {
    const session = fixture(`${name}.gbdraw-session.json.gz`);
    session.renderRequest.records = [];
    const { result, dialogs } = await load(session);
    assert.equal(result.status, 'error', name);
    assert.deepEqual([result.error?.code, result.error?.context?.field, result.error?.context?.reason],
      ['INPUT_INVALID', 'schema', 'FIELDS'], `${name}: ${JSON.stringify(result.error)}`);
    assert.deepEqual(dialogs, [], name);
  }
});

// Requests no reader can read at their root: a root the repair drops, none,
// and an object whose schema is no integer (OV-861, OV-863).
const UNREADABLE_REQUESTS = [{ schema: 9, bogus: 1 }, {}, { schema: 'nine' }, { schema: null, mode: 'linear' }];
/** Load refuses with the typed Session-fields error and commits nothing. */
const assertRefusedUnread = async (session, name) => {
  const before = getCommittedCanonicalSession();
  const { result, dialogs } = await load(session);
  assert.deepEqual([result.status, result.error?.code, result.error?.context?.field, result.error?.context?.reason],
    ['error', 'INPUT_INVALID', 'schema', 'FIELDS'], `${name}: ${JSON.stringify(result.error)}`);
  assert.deepEqual(dialogs, [], name);
  assert.equal(getCommittedCanonicalSession(), before, name);
};

test('a Session 31 to 41 whose request cannot be read at its root is refused with the typed Session-fields error (OV-861, OV-863, C-OD17)', async () => {
  for (const name of ['scale-interval-zero-circular-web.v33.gbdraw-session.json.gz', 'BGC0000708-BGC0000713.v39.gbdraw-session.json.gz',
    'HmmtDNA_ATskew.v40-schema5.json', 'cli-linear-tables.v41.gbdraw-session.json.gz']) {
    for (const request of UNREADABLE_REQUESTS) {
      const session = fixture(name);
      session.renderRequest = clone(request);
      await assertRefusedUnread(session, `${name} ${JSON.stringify(request)}`);
    }
  }
});

test('a schema 2 request takes the rows a schema 9 one takes: an empty tracks object, which no writer stored, is a part left out (PR-D)', () => {
  const tracksRows = (session) => readRequest(session, (request) => { request.diagramOptions.tracks = {}; }).repairs
    .filter(({ path }) => path.startsWith('renderRequest.diagramOptions.tracks'))
    .map(({ code, path, effect }) => ({ code, path, effect }));
  const dropped = [{ code: 'request-part-dropped', path: 'renderRequest.diagramOptions.tracks', effect: 'dropped' }];
  assert.deepEqual(tracksRows(gallery('lambda_basic_linear')), dropped);
  assert.deepEqual(tracksRows(fixture('feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz')), dropped);
});

test('a schema 1 or 2 request whose record names a resource the Session lacks is unreadable (K22), and its Session 31-41 is refused (OV-894, E-OD30)', async () => {
  const name = 'scale-interval-zero-circular-web.v33.gbdraw-session.json.gz';
  for (const schema of [1, 2]) {
    const session = fixture(name);
    session.renderRequest.schema = schema;
    delete session.resources['record-1-genbank'];
    assert.deepEqual(readStoredRenderRequest(session.renderRequest, contextOf(session)),
      { kind: 'request-dropped', path: 'renderRequest', reason: 'unknown-resource' }, String(schema));
  }
  const missing = fixture(name);
  delete missing.resources['record-1-genbank'];
  await assertRefusedUnread(missing, 'record resource missing');
  const none = fixture(name);
  none.resources = null;
  await assertRefusedUnread(none, 'resources null');
});

test('a request whose schema is no integer is unreadable (K22): settings only with its Web settings, refused without them; an integer the app does not read stays a Session-format error (OV-863)', async () => {
  for (const request of UNREADABLE_REQUESTS.slice(1)) {
    const session = gallery('lambda_basic_linear');
    session.renderRequest = clone(request);
    const { result, dialogs } = await load(session);
    assert.equal(result.status, 'ok', `${JSON.stringify(request)}: ${JSON.stringify(result.error)}`);
    assert.equal(state.results.value.length, 0);
    assert.deepEqual(result.repairs.map(({ code, path }) => [code, path]), [['request-dropped', 'renderRequest.schema']]);
    assert.match(dialogs[0][0], /^Saved diagram: the settings of the saved diagram \(renderRequest\.schema\) could not be read\. The file loads as settings only/);
    for (const name of ['se08-main-linear-cli.v42.gbdraw-session.json.gz', 'scale-interval-zero-circular-cli.v44.gbdraw-session.json.gz']) {
      const cli = fixture(name);
      cli.renderRequest = clone(request);
      await assertRefusedUnread(cli, `${name} ${JSON.stringify(request)}`);
    }
  }
  for (const schema of [3, 4, 10]) {
    const session = gallery('lambda_basic_linear');
    session.renderRequest.schema = schema;
    const { result } = await load(session);
    assert.deepEqual([result.error?.code, result.error?.context?.reason], ['INPUT_INVALID', 'SESSION_FORMAT'], String(schema));
  }
});

// R4 has one owner (E-OD22, C-OD27): a Session 42, 44 or 47 the CLI or the
// Python API wrote keeps no settings apart from its request, so when Load
// drops that request at its root (K22), `restoresNothing` refuses the file.
test('R4 alone refuses a CLI or Python API Session 42, 44 or 47 whose request is missing, unreadable or of no integer schema: nothing restored, nothing committed (E-OD21, E-OD22)', async () => {
  const directory = mkdtempSync(join(tmpdir(), 'gbdraw-cli-47-'));
  let cliCurrent;
  try {
    const file = join(directory, 'cli.gbdraw-session.json.gz');
    execFileSync('python', ['-m', 'gbdraw.cli', 'linear', '--gbk', fileURLToPath(new URL('tests/test_inputs/NC_001416.gb', REPO)),
      '-o', join(directory, 'cli'), '-f', 'svg', '--session_output', file], {
      cwd: directory, env: { ...process.env, PYTHONPATH: fileURLToPath(REPO) }, stdio: 'pipe'
    });
    cliCurrent = JSON.parse(gunzipSync(readFileSync(file)).toString('utf8'));
  } finally {
    rmSync(directory, { recursive: true, force: true });
  }
  const sessions = [
    ['se08-main-linear-cli.v42', fixture('se08-main-linear-cli.v42.gbdraw-session.json.gz')],
    ['scale-interval-zero-circular-cli.v44', fixture('scale-interval-zero-circular-cli.v44.gbdraw-session.json.gz')],
    ['cli.v47', { ...cliCurrent, version: 47 }],
    ['cli.current', cliCurrent]
  ];
  assert.equal(cliCurrent.version, SESSION_VERSION);
  for (const [name, base] of sessions) {
    assert.ok(!Object.hasOwn(base, 'config') && !Object.hasOwn(base, 'modes'), name);
    const missing = clone(base);
    delete missing.renderRequest;
    await assertRefusedUnread(missing, `${name} missing`);
    for (const request of [...UNREADABLE_REQUESTS, 'x', 5, []]) {
      const session = clone(base);
      session.renderRequest = clone(request);
      await assertRefusedUnread(session, `${name} ${JSON.stringify(request)}`);
    }
  }
});

test('a Session 27 to 30 loads as before whatever its renderRequest key holds: no writer wrote one and Load reads none (OV-862)', async () => {
  const name = 'BGC0000708-BGC0000713.v30.gbdraw-session.json.gz';
  const plain = await load(fixture(name));
  assert.equal(plain.result.status, 'ok', JSON.stringify(plain.result.error));
  const committed = JSON.stringify(getCommittedCanonicalSession());
  for (const request of [...UNREADABLE_REQUESTS, { schema: 10 }, { schema: 3 }]) {
    const session = fixture(name);
    session.renderRequest = clone(request);
    const { result } = await load(session);
    assert.equal(result.status, 'ok', `${JSON.stringify(request)}: ${JSON.stringify(result.error)}`);
    assert.deepEqual(result.repairs, plain.result.repairs, JSON.stringify(request));
    assert.equal(JSON.stringify(getCommittedCanonicalSession()), committed, JSON.stringify(request));
  }
  // Beside `resources`, a readable schema-5 request is still not committed, and
  // Save still asks for Generate first (OV-819).
  const v40 = fixture('BGC0000708-BGC0000713.v40-schema5.json');
  const session = fixture(name);
  Object.assign(session, { resources: v40.resources, renderRequest: v40.renderRequest });
  const { result } = await load(session);
  assert.equal(result.status, 'ok', JSON.stringify(result.error));
  assert.deepEqual(result.repairs, plain.result.repairs);
  assert.equal(JSON.stringify(getCommittedCanonicalSession()), committed);
  await assert.rejects(exportSession('probe', { askChoice: async () => 'download' }), { code: 'SESSION_SAVE_REQUIRES_GENERATE' });
});

const HARNESS_CANNOT_LOAD = new Map([
  // Their Load extracts the features again, which needs the Pyodide Worker (WORKER_INIT here).
  ['feature-edits-circular.v33.gbdraw-session.json.gz', 'feature extraction Worker'],
  ['feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz', 'feature extraction Worker'],
  ['hash-color-rule-linear-rc.v30.gbdraw-session.json.gz', 'feature extraction Worker'],
  // Their comparison tables are converted to the display frame by a Worker
  // (`convertMainSessionComparisonFrames`; OV-775, C-OD19 F2).
  ['q-frame-main-web-losatn.v42.gbdraw-session.json.gz', 'comparison conversion Worker'],
  ['q-frame-main-web-upload.v42.gbdraw-session.json.gz', 'comparison conversion Worker']
]);
const REQUEST_ROW_CODES = new Set(['request-dropped', 'other-result-set-dropped', 'request-part-dropped', 'request-reference-dropped']);
const requestRow = ({ code, path }) => REQUEST_ROW_CODES.has(code) || /^(otherModeResult\.)?renderRequest\b/.test(path);

test('every Gallery Session and Session fixture loads with no stored-request row and commits the current schema (G3)', { timeout: 600000 }, async () => {
  const { readdirSync } = await import('node:fs');
  // Node has no diagram Worker: Python's conversion of the protein results a
  // Session up to 47 stored for a reversed record (OV-1145) answers "nothing changed".
  const { setStoredProteinResultFrameConverter } = await import('../../gbdraw/web/js/services/main-session-comparison-frame.js');
  setStoredProteinResultFrameConverter(async () => ({ result: { resources: {} } }));
  const findings = {};
  let loaded = 0;
  for (const folder of ['gbdraw/web/gallery/sessions/', 'tests/fixtures/sessions/']) {
    for (const name of readdirSync(new URL(folder, REPO)).filter((file) => /\.json(\.gz)?$/.test(file)).sort()) {
      const session = read(folder + name);
      if (session?.format !== 'gbdraw-session' || HARNESS_CANNOT_LOAD.has(name)) continue;
      const stored = session.renderRequest?.schema ?? null;
      const { result } = await load(session);
      const committed = result.status === 'ok' || result.status === 'legacy'
        ? getCommittedCanonicalSession()?.renderRequest?.schema ?? null : 'not loaded';
      const expected = stored === null ? null : 9;
      const rows = (result.repairs || []).filter(requestRow).map(({ code, path }) => `${code} ${path}`);
      if (rows.length || committed !== expected) findings[folder + name] = { status: result.status, rows, committed, expected };
      loaded += 1;
    }
  }
  setStoredProteinResultFrameConverter(null);
  assert.deepEqual(findings, {});
  assert.ok(loaded >= 50, String(loaded));
});

// G2 on the request path (S6): each leaf of a stored request is given a value
// of another kind and an unknown value of its kind, the file is loaded, and the
// committed request may differ from the stored one only where a row names. A
// setting of the other diagram type is read at the one value its type holds,
// without a row (〔出さないもの〕); no mutation is a described spelling.
const { CANONICAL_REQUEST_DESCRIPTION } = await import('../../gbdraw/web/js/contracts/request-shape.generated.js');
const requestType = (name) => {
  const shape = CANONICAL_REQUEST_DESCRIPTION.types[name];
  return shape.kind === 'ref' ? requestType(shape.name) : shape;
};
/** The one value a place of `mode`'s request holds, if its type holds one. */
const onlyValue = (mode, parts) => {
  let shape = requestType(mode === 'linear' ? 'CanonicalLinearRenderRequest' : 'CanonicalCircularRenderRequest');
  for (const part of parts) {
    const next = shape.kind === 'object' ? shape.fields[part] : shape.kind === 'array' ? shape.items
      : shape.kind === 'record' ? shape.values : undefined;
    if (!next) return undefined;
    shape = next.kind === 'ref' ? requestType(next.name) : next;
  }
  if (shape.kind === 'null') return { value: null };
  if (shape.kind === 'tuple' && shape.items.length === 0) return { value: [] };
  return shape.kind === 'literal' && shape.values.length === 1 ? { value: shape.values[0] } : undefined;
};
const REQUEST_MUTATIONS = {
  otherKind: (value) => (value === null ? 'x-mut' : typeof value === 'object' || typeof value === 'string' ? 5 : 'x-mut'),
  sameKind: (value) => (typeof value === 'string' ? 'x-mut' : typeof value === 'number' ? -7.5 : undefined)
};
const placesOf = (value, parts = []) => [
  ...(parts.length ? [parts] : []),
  ...(value && typeof value === 'object' ? Object.entries(value).flatMap(([key, child]) => (
    placesOf(child, [...parts, Array.isArray(value) ? Number(key) : key]))) : [])
];
const differences = (stored, committed, parts = []) => {
  if (JSON.stringify(stored) === JSON.stringify(committed)) return [];
  if (stored && committed && typeof stored === 'object' && typeof committed === 'object'
    && Array.isArray(stored) === Array.isArray(committed)) {
    return [...new Set([...Object.keys(stored), ...Object.keys(committed)])].flatMap((key) => differences(
      stored[key], committed[key], [...parts, Array.isArray(stored) ? Number(key) : key]
    ));
  }
  return [parts];
};
const placeName = (parts) => ['renderRequest', ...parts].map((part) => (typeof part === 'number' ? `[${part}]` : `.${part}`)).join('').slice(1);
const namedBy = (path, repairs) => repairs.some((repair) => repair.path === path || path.startsWith(`${repair.path}.`)
  || path.startsWith(`${repair.path}[`) || repair.path.startsWith(`${path}.`) || repair.path.startsWith(`${path}[`));

test('Load commits a stored request that differs from the stored one only where a row names (G2, S6)', async () => {
  const silent = [];
  let checked = 0;
  for (const name of ['lambda_basic_linear', 'HmmtDNA_basic_circular']) {
    const session = gallery(name);
    const { mode } = session.renderRequest;
    for (const parts of placesOf(session.renderRequest)) {
      for (const [kind, mutate] of Object.entries(REQUEST_MUTATIONS)) {
        const value = mutate(parts.reduce((current, part) => current[part], session.renderRequest));
        if (value === undefined) continue;
        const file = clone(session);
        const owner = parts.slice(0, -1).reduce((current, part) => current[part], file.renderRequest);
        owner[parts.at(-1)] = value;
        const { result } = await load(file);
        if (result.status !== 'ok') continue;
        checked += 1;
        const committed = getCommittedCanonicalSession()?.renderRequest ?? null;
        for (const place of differences(file.renderRequest, committed)) {
          const only = onlyValue(mode, place);
          if (namedBy(placeName(place), result.repairs)) continue;
          if (only && JSON.stringify(place.reduce((current, part) => current?.[part], committed)) === JSON.stringify(only.value)) continue;
          silent.push(`${name} ${kind} ${placeName(parts)} = ${JSON.stringify(value)}: ${placeName(place)}`);
        }
      }
    }
  }
  assert.ok(checked > 200, String(checked));
  assert.deepEqual(silent, [], String(checked));
});

// H2 (PR-D): Load lifts a schema-1 request and commits schema 9, so a live
// edit (the committed request with the drawing's tables) works before any
// Generate, and Save writes schema 9. Save still asks for one Generate while
// the Session holds Results without their feature metadata (D-25, PD-OI-079),
// as it does for a schema-2 Session; the same Session without Results saves.
test('a schema-1 Session takes a live edit and Save before any Generate (H2)', async () => {
  const save = (name) => exportSession(name, { askChoice: async () => 'download' });
  for (const name of ['lambda_basic_linear.v31-schema1', 'HmmtDNA_basic_circular.v32-schema1']) {
    const session = read(`tests/fixtures/sessions/${name}.gbdraw-session.json.gz`);
    assert.equal(session.renderRequest.schema, 1);
    assert.ok(session.results.length > 0);
    assert.equal((await load(session)).result.status, 'ok', name);
    await assert.rejects(save(name), { code: 'SESSION_SAVE_REQUIRES_GENERATE' }, name);
    assert.equal((await load({ ...session, results: [] })).result.status, 'ok', name);
    const committed = getCommittedCanonicalSession();
    assert.equal(committed.renderRequest.schema, 9, name);
    const edited = projectCommittedEditorIntent({ committed, state, drawing: state.activeDrawing() });
    assert.equal(edited.renderRequest.schema, 9, name);
    const saved = await save(name);
    assert.equal(saved.status, 'saved', `${name}: ${saved.error?.code}`);
    const document = JSON.parse(gunzipSync(Buffer.from(await saved.blob.arrayBuffer())).toString('utf8'));
    assert.equal(document.renderRequest.schema, 9, name);
  }
});
