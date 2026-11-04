// New drawing from a region (PD-OI-090 〔作成〕(4)): the shown drawing's draft
// request goes to the real Python derivation; the derived drawing is read and
// installed as Load installs a stored drawing, and its request round-trips.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { createHash } from 'node:crypto';
import { DRAWING_SLICES, modeDrawing } from './helpers/drawing-state.mjs';
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
const { createDrawingList } = await import('../../gbdraw/web/js/app/drawing-list.js');
const { allocateDrawingId, defaultDrawingName } = await import('../../gbdraw/web/js/services/drawing-view.js');
const { withNewDrawingDefaults } = await import('../../gbdraw/web/js/services/reset.js');
const { createAnnotationSet } = await import('../../gbdraw/web/js/services/annotation-state.js');
const { coordinateTarget } = await import('../../gbdraw/web/js/app/annotations/target-actions.js');
const { bindDrawingInputFile } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const {
  installShownDrawing, installPreparedDrawing, prepareRegionDrawing, prepareStoredDrawing, projectDraftRequest
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

const STAGE = `
import base64, json, os, tempfile
directory = tempfile.mkdtemp()
paths = {}
for index, (resource_id, descriptor) in enumerate(payload['resources'].items()):
    path = os.path.join(directory, f'r{index}_{descriptor.get("name") or "resource"}')
    data = descriptor['data']
    with open(path, 'wb') as handle:
        handle.write(base64.b64decode(data) if descriptor.get('encoding') == 'base64' else data.encode())
    paths[resource_id] = path
`;

/** The biological feature IDs of a request's records by start (0-based), as Python reads them. */
const featureIdsInPython = async ({ renderRequest, resources }) => JSON.parse(await evaluatePython(`${STAGE}
from gbdraw.api.request_render import resolve_request_records
from gbdraw.session_request_codec import decode_canonical_request
request = decode_canonical_request(payload['canonicalRequest'], resource_paths=paths, output_directory=directory)
result = json.dumps({
    str(feature.location_parts[0][0]): feature.biological_feature_id
    for item in resolve_request_records(request).provenance for feature in item.source_feature_catalog
})
`, { canonicalRequest: renderRequest, resources }));

/** The Worker helper `deriveRegionDrawing`, run by the real Python owner over staged resources. */
const deriveInPython = async (payload) => ({ result: JSON.parse(await evaluatePython(`${STAGE}
from gbdraw.web_support.region_drawing import derive_web_region_drawing_json
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

// Every slice of a drawing and its input type (as history-per-drawing hashes them).
const hashOf = (drawingId) => {
  const drawing = state.drawingOf(drawingId);
  const slices = Object.fromEntries(Object.entries(DRAWING_SLICES)
    .map(([name, slice]) => [name, slice.capture(drawing.slices[name])]));
  return createHash('sha256').update(JSON.stringify({ slices, inputType: drawing.inputType.value })).digest('hex');
};

// The drawing list over the store and the app's History, with the region
// port as the app builds it; a switch only changes the shown drawing.
const drawingList = createDrawingList({
  drawingTabs: () => state.drawingTabs.value,
  shownDrawingId: () => state.activeDrawingId.value,
  switchBusy: () => null,
  drawingOf: state.drawingOf,
  drawingHasResult: () => false,
  drawingHasEdits: () => false,
  boundFileIdsOf: () => [],
  createDrawing: (drawingId, mode) => withNewDrawingDefaults(state, createDrawingState(drawingId, mode)),
  styleDrawing: () => { throw new Error('not used'); },
  copyDrawing: () => { throw new Error('not used'); },
  regionDrawing: async (derivation) => {
    const prepared = await prepareRegionDrawing(derivation);
    return (drawingId) => {
      const made = createDrawingState(drawingId, prepared.mode);
      installPreparedDrawing(made, prepared);
      return made;
    };
  },
  addDrawing: state.addDrawing,
  removeDrawing: (drawingId) => state.removeDrawing(drawingId, { takeSlot: () => null }),
  putSlot: () => {},
  renameDrawing: state.renameDrawing,
  moveDrawing: state.moveDrawing,
  addDrawingHistory: history.addMember,
  removeDrawingHistory: history.removeMember,
  runListStep: history.runListStep,
  showDrawingNow: installShownDrawing
});
const regionDrawingsOf = () => createRegionDrawing({
  shownDrawingId: () => state.activeDrawingId.value,
  projectShownRequest: () => projectDraftRequest(state, state.activeDrawing(), facts),
  deriveRegion: deriveInPython,
  newDrawingFromRegion: drawingList.newDrawingFromRegion
});

test('(d) Create adds the region drawing as New does; Undo removes it, Redo brings it back', async () => {
  const source = circularSource([recordA]);
  const sourceHash = hashOf(source.id);
  const before = listView();
  const offer = await regionDrawingsOf().deriveFromShown([{ recordKey: 'record-1', start: 301, end: 700 }], true);
  assert.equal(offer.status, 'derived');
  assert.equal(offer.mode, 'linear', 'a part of a record is drawn Linear');
  const drawingId = allocateDrawingId('linear', before.tabs.map(({ id }) => id));
  assert.deepEqual(listView(), before, 'deriving adds nothing');

  const outcome = await offer.create();
  assert.deepEqual(outcome, { status: 'ok', drawingId });
  assert.deepEqual(state.drawingTabs.value.map(({ id, name }) => ({ id, name })),
    [...before.tabs, { id: drawingId, name: defaultDrawingName('linear', drawingId) }], 'the ID and name New gives');
  assert.equal(state.activeDrawingId.value, drawingId, 'the new drawing is shown');
  assert.equal(state.drawingOf(drawingId).mode, 'linear');
  assert.equal(history.timeline().undo[0]?.label, 'New drawing from a region');
  assert.equal(history.managerOf(drawingId).getUndoCount(), 0, 'its History starts empty');
  assert.equal(hashOf(source.id), sourceHash, 'the source drawing is unchanged');

  await history.undo();
  assert.deepEqual(state.drawingTabs.value.map(({ id, name }) => ({ id, name })), before.tabs);
  assert.equal(state.activeDrawingId.value, before.shown, 'the drawing shown before');
  assert.equal(hashOf(source.id), sourceHash);

  await history.redo();
  assert.equal(state.activeDrawingId.value, drawingId);
  assert.equal(state.drawingOf(drawingId).mode, 'linear');
  assert.equal(hashOf(source.id), sourceHash);
  // Leave the store as the next cases expect it.
  await history.undo();
});

// The names of the files registered while `run` runs.
const countingRegistrations = async (run) => {
  const { registerFile } = state.inputFiles;
  const registered = [];
  state.inputFiles.registerFile = (file) => { if (file) registered.push(file.name); return registerFile(file); };
  try {
    return { value: await run(), registered };
  } finally {
    state.inputFiles.registerFile = registerFile;
  }
};

// A Circular record with a per-feature edit and an annotation inside 301..800
// and one of each outside, and a color rule whose caption renames its Legend row.
const recordWithLook = genbankRecord('recA', 1200, [[101, 200, 'A1'], [401, 600, 'A2'], [901, 1100, 'A3']]);
const LOOK_REGION = [{ recordKey: 'record-1', start: 301, end: 800 }];
const RULE = { feat: 'CDS', qual: 'locus_tag', val: 'A2', color: '#ff0000', cap: 'Kinase' };
const circularSourceWithLook = async () => {
  const source = circularSource([recordWithLook]);
  for (const key of Object.keys(source.featureOverrides)) delete source.featureOverrides[key];
  source.annotationSets.splice(0);
  source.manualSpecificRules.splice(0);
  const ids = await featureIdsInPython((await projectDraftRequest(state, source, facts)).envelope);
  for (const [start, labelText] of [['400', 'kept'], ['900', 'dropped']]) {
    source.featureOverrides[JSON.stringify(['record-1', ids[start]])] = {
      recordKey: 'record-1', biologicalFeatureId: ids[start], featureVisibility: null, labelVisibility: null,
      labelText, labelSourceText: null
    };
  }
  source.annotationSets.push(createAnnotationSet({ id: 'marks', annotations: [
    { id: 'in-span', target: coordinateTarget({ start: 450, end: 500 }), label: 'Inside', mark: 'band' },
    { id: 'out-span', target: coordinateTarget({ start: 1000, end: 1050 }), label: 'Outside', mark: 'band' }
  ] }));
  source.manualSpecificRules.push({ ...RULE });
  return source;
};
/** The parts of the look a region drawing carries: inside ones arrive, outside ones do not. */
const assertLookCarried = (drawing) => {
  assert.deepEqual(Object.values(drawing.featureOverrides).map((row) => row.labelText), ['kept']);
  assert.deepEqual(drawing.annotationSets.flatMap((set) => set.annotations.map((item) => item.label)), ['Inside']);
  assert.deepEqual(drawing.manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({ feat, qual, val, color, cap })),
    [RULE], 'the color rule and its Legend name');
};

test('(a) a region of one Circular record is a Linear drawing that carries the look inside the region', async () => {
  const source = await circularSourceWithLook();
  const { value: installed, registered } = await countingRegistrations(() => deriveAndInstall(source, LOOK_REGION));
  assert.equal(installed.made.mode, 'linear');
  assert.deepEqual(installed.prepared.repairs, []);
  // The color table is a table the request writes, which Load holds as its card's file.
  assert.deepEqual(registered, ['specific-colors.tsv'], 'no new GenBank file');
  assert.equal(installed.made.files.c_gb ?? null, null);
  assert.deepEqual(installed.made.linearSeqs.map((row) => row.gb?.fileId), [source.files.c_gb.fileId],
    'the source drawing\'s Web file');
  assertLookCarried(installed.made);
});

test('(f) the same region of a drawing loaded from a Session, without a Result', async () => {
  const sent = await projectDraftRequest(state, await circularSourceWithLook(), facts);
  const prepared = await prepareStoredDrawing({
    mode: 'circular', renderRequest: sent.envelope.renderRequest, resources: sent.envelope.resources, sourceFiles: new Map()
  });
  assert.deepEqual(prepared.repairs, []);
  const loaded = createDrawingState('loaded', 'circular');
  installPreparedDrawing(loaded, prepared);
  const loadedFile = loaded.files.c_gb;
  assert.notEqual(loadedFile.fileId, modeDrawing(state, 'circular').files.c_gb.fileId, 'Load made its own Web file');
  const { value: installed, registered } = await countingRegistrations(() => deriveAndInstall(loaded, LOOK_REGION));
  assert.equal(installed.made.mode, 'linear');
  assert.deepEqual(installed.prepared.repairs, []);
  // The color table is a table the request writes, which Load holds as its card's file.
  assert.deepEqual(registered, ['specific-colors.tsv'], 'no new GenBank file');
  assert.deepEqual(installed.made.linearSeqs.map((row) => row.gb?.fileId), [loadedFile.fileId],
    'the loaded drawing\'s Web file (a Session resource view)');
  assertLookCarried(installed.made);
});

// The round trip of a Linear region drawing needs the Linear record catalog
// (`prepareLinearRecordCatalog`, which reads the records in the Worker); this
// file does not build it yet, so the draft projection misses the GenBank file.
test('(a) round trip of the Linear region drawing of a Circular source', { todo: 'B5-6: build the Linear record catalog in node' }, async () => {
  await assertRoundTrip(await deriveAndInstall(await circularSourceWithLook(), LOOK_REGION));
});
