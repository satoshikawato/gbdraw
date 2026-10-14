// The one reader of a stored render request at Load (`readStoredRenderRequest`
// in services/session-request.js, PR-C step C2c): the envelope, the schema
// 5-8 promotion with the Load context, the repair of a request against its
// description with the request policy (references, ring columns, settings of
// the other diagram type), the root (K22), and its Load wiring for the
// top-level request and `otherModeResult` (OV-670).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { File } from 'node:buffer';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

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
const { importSession, setUnmanagedConfigOverrideValidator } = await import('../../gbdraw/web/js/services/config.js');
const { readStoredRenderRequest, readPublishedRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
setUnmanagedConfigOverrideValidator((payload) => ({ result: JSON.parse(execFileSync('python', ['-c', `
import json, sys
from gbdraw.web_support.config_overrides import validate_web_config_overrides_json
p = json.load(sys.stdin)
print(validate_web_config_overrides_json(p['mode'], json.dumps(p['config']),
    json.dumps(p['configOverrides']), json.dumps(p['managedPaths']), p['requireUnmanagedOnly']))
`], { input: JSON.stringify(payload), encoding: 'utf8', cwd: fileURLToPath(REPO) })) }));
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { state } = await import('../../gbdraw/web/js/state.js');

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

test('a request schema the app does not read is a Session-format error (OV-773)', () => {
  for (const schema of [3, 4, 10, undefined]) {
    assert.throws(
      () => readRequest(gallery('lambda_basic_linear'), (request) => { request.schema = schema; }),
      (error) => error.code === 'INPUT_INVALID' && error.context?.reason === 'SESSION_FORMAT',
      String(schema)
    );
  }
});

test('a request that needs no repair comes back as the stored object (RK-3); schema 1 and 2 as stored', () => {
  const session = gallery('lambda_basic_linear');
  const result = readStoredRenderRequest(session.renderRequest, contextOf(session));
  assert.equal(result.value, session.renderRequest);
  assert.deepEqual(result.repairs, []);
  for (const schema of [1, 2]) {
    const stored = { schema, mode: 'linear', records: [] };
    assert.equal(readStoredRenderRequest(stored, contextOf(session)).value, stored);
  }
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

test('an annotation colour the decoder refuses takes its default; another spelling of a name reads as it (C-OD11)', () => {
  const session = gallery('tobacco-chloroplast');
  const result = readRequest(session, (request) => {
    const style = request.diagramOptions.annotations.sets[0].annotations[0].style;
    Object.assign(style, { stroke: 'valueTheDescriptionDoesNotName', fill: '#abcd', labelColor: ' Navy ' });
  });
  const style = result.value.diagramOptions.annotations.sets[0].annotations[0].style;
  assert.deepEqual([style.stroke, style.fill, style.labelColor], ['#404040', null, 'navy']);
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
  const twoMode = () => JSON.parse(execFileSync('python', ['-c',
    'import json; from tests.utils.two_mode_session import two_mode_session; print(json.dumps(two_mode_session()))'
  ], { encoding: 'utf8', cwd: fileURLToPath(REPO), maxBuffer: 1 << 30 }));
  const listed = twoMode();
  listed.otherModeResult.renderRequest.records[0].bogus = true;
  const first = await load(listed);
  assert.equal(first.result.status, 'ok', JSON.stringify(first.result.error));
  assert.deepEqual(first.dialogs.flat(), [
    'Settings: 1 entry of otherModeResult.renderRequest.records[0].bogus could not be read and was dropped.'
  ]);
  const dropped = twoMode();
  dropped.otherModeResult.renderRequest.records = [];
  const second = await load(dropped);
  assert.equal(second.result.status, 'ok', JSON.stringify(second.result.error));
  assert.deepEqual(second.dialogs.flat(), [
    'Saved Linear diagram: its settings (otherModeResult.renderRequest) could not be read; that diagram and its Results are left out.'
  ]);
  assert.ok(state.results.value.length > 0);
});

const { getCommittedCanonicalSession, exportSession } = await import('../../gbdraw/web/js/services/config.js');
const frame = () => {
  const intent = state.drawings.linear.importedComparisonIntent;
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
    ['dangling comparison', () => {
      const session = gallery('BGC0000708-BGC0000713');
      session.renderRequest.comparisons[0].resourceId = 'missing';
      return session;
    }, 'DECISION_REQUIRED', [], missing],
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
// request and commits it at the current schema (schema 1 as stored until its
// reader, PR-D). A row here is a defect of the input's producer.
/** Inputs this Node harness cannot load, by name, with the reason. */
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
      const expected = stored === null ? null : stored === 1 ? 1 : 9;
      const rows = (result.repairs || []).filter(requestRow).map(({ code, path }) => `${code} ${path}`);
      if (rows.length || committed !== expected) findings[folder + name] = { status: result.status, rows, committed, expected };
      loaded += 1;
    }
  }
  assert.deepEqual(findings, {});
  assert.ok(loaded >= 50, String(loaded));
});
