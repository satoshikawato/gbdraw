// New drawing from a region (PD-OI-090 〔作成〕(4)): the shown drawing's draft
// request goes to the real Python derivation; the derived drawing is read and
// installed as Load installs a stored drawing, and its request round-trips.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';
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
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');
const { createRegionDrawing } = await import('../../gbdraw/web/js/app/region-drawing.js');
const { bindDrawingInputFile } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const {
  installShownDrawing, installPreparedDrawing, prepareRegionDrawing, projectDraftRequest, serializeActiveRenderFiles
} = await import('../../gbdraw/web/js/services/config.js');
const { buildCanonicalRenderRequest, compareCanonicalRenderRequests } = await import(
  '../../gbdraw/web/js/services/session-request.js'
);

/** A GenBank record of `length` bp with one CDS per `[start, end, tag]`. */
const genbankRecord = (name, length, cds, topology = 'circular') => {
  const features = cds.map(([start, end, tag]) => [
    `     CDS             ${start}..${end}`,
    `                     /locus_tag="${tag}"`,
    `                     /product="${tag} protein"`
  ].join('\n')).join('\n');
  const sequence = 'atgcgtacgt'.repeat(Math.ceil(length / 10)).slice(0, length);
  const lines = [];
  for (let offset = 0; offset < length; offset += 60) {
    const blocks = sequence.slice(offset, offset + 60).match(/.{1,10}/g).join(' ');
    lines.push(`${String(offset + 1).padStart(9)} ${blocks}`);
  }
  return [
    `LOCUS       ${name.padEnd(16)}${String(length).padStart(11)} bp    DNA     ${topology} BCT 01-JAN-2000`,
    `DEFINITION  ${name} test record.`, `ACCESSION   ${name}`, `VERSION     ${name}.1`,
    'FEATURES             Location/Qualifiers', `     source          1..${length}`, features,
    'ORIGIN', ...lines, '//'
  ].join('\n');
};
const genbankFile = (name, records) => new File([`${records.join('\n')}\n`], name, { type: 'text/plain', lastModified: 1 });

/** The Worker helper `deriveRegionDrawing`, run by the real Python owner over staged resources. */
const deriveInPython = async (payload) => ({ result: JSON.parse(await evaluatePython(`
import base64, json, os, tempfile
from gbdraw.web_support.region_drawing import derive_web_region_drawing_json
directory = tempfile.mkdtemp()
paths = {}
for index, (resource_id, descriptor) in enumerate(payload['resources'].items()):
    path = os.path.join(directory, f'r{index}_{descriptor.get("name") or "resource"}')
    data = descriptor['data']
    with open(path, 'wb') as handle:
        handle.write(base64.b64decode(data) if descriptor.get('encoding') == 'base64' else data.encode())
    paths[resource_id] = path
result = derive_web_region_drawing_json(
    json.dumps(payload['canonicalRequest']), json.dumps(paths), directory,
    json.dumps(payload['regions']), json.dumps(payload['adaptSizes']))
`, payload)) });

const facts = { linearRecordCatalog: null, recordDisplayRows: [], omittedAnnotationRecordKeys: [] };

const circularSource = (records, form = {}) => {
  const drawing = installShownDrawing('circular');
  drawing.inputType.value = 'gb';
  bindDrawingInputFile(state.inputFiles, drawing.slices.inputs, 'c_gb', genbankFile('c.gb', records));
  drawing.circularRecordList.value = [];
  Object.assign(drawing.form, { labels_mode: 'none', circular_record_selector: '', multi_record_canvas: false, ...form });
  return drawing;
};

/** The source's request, derived for `regions`, read and installed under `drawingId`. */
const deriveAndInstall = async (source, regions, adaptSizes = true, drawingId = 'region-test') => {
  const sent = await projectDraftRequest(state, source, facts);
  const { renderRequest, resources } = sent.envelope;
  const { result: derived } = await deriveInPython({ canonicalRequest: renderRequest, resources, regions, adaptSizes });
  assert.equal(derived.status, 'derived', JSON.stringify(derived));
  const prepared = await prepareRegionDrawing({ derived, sent });
  const made = createDrawingState(drawingId, prepared.mode);
  installPreparedDrawing(made, prepared);
  return { sent, derived, prepared, made };
};

/** The derived drawing's resources: each one the bytes of the sent resource it reads. */
const derivedResources = ({ derived, sent }) => Object.fromEntries(Object.entries(derived.resourceSources)
  .map(([resourceId, sourceId]) => [resourceId, sent.envelope.resources[sourceId]]));

/** The installed drawing's draft request, as Save projects it, equals the derived request. */
const assertRoundTrip = async (installed) => {
  const { envelope } = await projectDraftRequest(state, installed.made, facts);
  const comparison = await compareCanonicalRenderRequests({
    expectedRequest: installed.derived.renderRequest, expectedResources: derivedResources(installed),
    actualRequest: envelope.renderRequest, actualResources: envelope.resources
  });
  assert.deepEqual(comparison.differences, [], 'the installed drawing writes the derived request');
  assert.equal(comparison.equivalent, true);
};

createAppSetup();
const history = window.__GBDRAW_HISTORY__;
const listView = () => ({
  tabs: state.drawingTabs.value.map(({ id, name }) => ({ id, name })),
  shown: state.activeDrawingId.value,
  timeline: JSON.stringify(history.timeline())
});
const recordA = genbankRecord('recA', 1200, [[101, 400, 'A1'], [601, 900, 'A2']]);

test('(b) a whole Circular record keeps Circular and binds the source drawing\'s file', async () => {
  const source = circularSource([recordA]);
  const installed = await deriveAndInstall(source, [{ recordKey: 'record-1', start: 1, end: 1200 }]);
  assert.equal(installed.derived.mode, 'circular');
  assert.equal(installed.made.mode, 'circular');
  assert.deepEqual(installed.prepared.repairs, []);
  assert.deepEqual(installed.derived.returnedToAuto, []);
  assert.equal(installed.made.files.c_gb.fileId, source.files.c_gb.fileId, 'the source drawing\'s Web file');
  await assertRoundTrip(installed);
});

test('(c) a region that crosses the origin is refused: no drawing, the tabs and History unchanged', async () => {
  circularSource([recordA]);
  const before = listView();
  const created = [];
  const regions = createRegionDrawing({
    shownDrawingId: () => state.activeDrawingId.value,
    projectShownRequest: () => projectDraftRequest(state, state.activeDrawing(), facts),
    deriveRegion: deriveInPython,
    newDrawingFromRegion: async (...args) => { created.push(args); return { status: 'ok', drawingId: 'never' }; }
  });
  const offer = await regions.deriveFromShown([{ recordKey: 'record-1', start: 900, end: 100 }], true);
  assert.deepEqual(offer, { status: 'refused', reason: 'ORIGIN_CROSSING', recordKey: 'record-1', start: 900, end: 100 });
  assert.equal('create' in offer, false, 'a refusal offers no Create');
  assert.deepEqual(created, []);
  assert.deepEqual(listView(), before);
});
