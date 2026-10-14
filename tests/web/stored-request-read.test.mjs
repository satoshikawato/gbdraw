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
const { readStoredRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
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
  const axis = readRequest(gallery('lambda_basic_linear'), (request) => { request.diagramOptions.tracks.linearTrackAxisIndex = 3; });
  assert.equal(axis.value.diagramOptions.tracks, undefined);
  assert.deepEqual(rows(axis), [{ code: 'request-reference-dropped', path: 'renderRequest.diagramOptions.tracks' }]);
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
