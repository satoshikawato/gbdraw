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
const { bindDrawingInputFile, createLinearSeq } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const {
  installShownDrawing, installPreparedDrawing, makeRegionDrawing, prepareRegionDrawing, prepareStoredDrawing, projectDraftRequest
} = await import('../../gbdraw/web/js/services/config.js');
const { compareCanonicalRenderRequests } = await import('../../gbdraw/web/js/services/session-request.js');
const { annotationSourceKey, buildAnnotationRecordCatalog } = await import(
  '../../gbdraw/web/js/app/annotations/record-catalog.js'
);
const { discoverSequenceRecords } = await import('../../gbdraw/web/js/services/record-discovery.js');
const { createDiagramResourceTransport } = await import('../../gbdraw/web/js/services/diagram-resource-staging.js');

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

// The bytes the page stages to the Worker for a request (the staging owner,
// `createDiagramResourceTransport`), by resource ID; Python writes each to a file.
/** @param {unknown} request @param {Record<string, unknown>} resources */
const stagedFiles = async (request, resources) => {
  const { resourceManifest, stagedResources } = await createDiagramResourceTransport().prepare({ request, resources });
  const bytesByToken = new Map(stagedResources.map(({ cacheToken, bytes }) => [cacheToken, Buffer.from(bytes).toString('base64')]));
  return Object.fromEntries(resourceManifest.map(({ resourceId, cacheToken, name }) => [resourceId, { name, base64: bytesByToken.get(cacheToken) }]));
};
const STAGE = `
import base64, json, os, tempfile
directory = tempfile.mkdtemp()
paths = {}
for index, (resource_id, staged) in enumerate(payload['staged'].items()):
    path = os.path.join(directory, f'r{index}_{staged["name"]}')
    with open(path, 'wb') as handle:
        handle.write(base64.b64decode(staged['base64']))
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
`, { canonicalRequest: renderRequest, staged: await stagedFiles(renderRequest, resources) }));

/** The Worker helper `deriveRegionDrawing`, run by the real Python owner over staged resources. */
const deriveInPython = async (payload) => ({ result: JSON.parse(await evaluatePython(`${STAGE}
from gbdraw.web_support.region_drawing import derive_web_region_drawing_json
result = derive_web_region_drawing_json(
    json.dumps(payload['canonicalRequest']), json.dumps(paths), directory,
    json.dumps(payload['regions']), json.dumps(payload['adaptSizes']))
`, { ...payload, resources: undefined, staged: await stagedFiles(payload.canonicalRequest, payload.resources) })) });

const noWorker = () => { throw new Error('the records are read by the fast path'); };
/**
 * The facts the app's projection of a drawing reads. A Linear drawing's record
 * catalog is built as `prepareLinearRecordCatalog({ privateCandidate: true })`
 * builds it: the records each card's file holds (read here by the fast path,
 * without the Worker), selected by the card's record selector.
 */
const factsOf = async (drawing) => {
  if (drawing.mode !== 'linear') return { linearRecordCatalog: null, recordDisplayRows: [], omittedAnnotationRecordKeys: [] };
  const inputType = drawing.inputType.value;
  const recordsOff = new Set(drawing.recordsOff);
  const linearSources = await Promise.all(drawing.linearSeqs.map(async (seq) => ({
    sourceKey: annotationSourceKey({ scope: 'linear', uid: seq.uid, inputType, primaryFile: seq.gb, pairedFile: null }),
    selector: seq.region_record_id, drawn: !recordsOff.has(seq.uid), hasInput: Boolean(seq.gb), status: 'ready',
    records: await discoverSequenceRecords({
      file: state.inputFiles.fileOf(seq.gb), format: 'genbank', runHelperOperation: noWorker
    })
  })));
  const catalog = buildAnnotationRecordCatalog({ mode: 'linear', linearSources });
  assert.equal(catalog.status, 'ready', JSON.stringify(catalog.issues));
  return {
    linearRecordCatalog: catalog, recordDisplayRows: [], omittedAnnotationRecordKeys: catalog.omittedRecordKeys || []
  };
};
/** The draft request of `drawing`, shown or not (the request's mode is the drawing's own). */
const projectOf = async (drawing) => projectDraftRequest(state, drawing, await factsOf(drawing));

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
  const sent = await projectOf(source);
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

// A config override's key is a dotted path of its own.
const CONFIG_OVERRIDES = '$.diagramOptions.configOverrides.';
const keysOf = (path) => (path.startsWith(CONFIG_OVERRIDES)
  ? ['diagramOptions', 'configOverrides', path.slice(CONFIG_OVERRIDES.length)] : path.slice(2).split('.'));
/** The value at `path` (`$.a.b`) of `request`, and whether it is there. */
const at = (request, path) => keysOf(path).reduce(
  ({ value }, key) => (value && typeof value === 'object' && Object.hasOwn(value, key)
    ? { value: value[key], found: true } : { value: undefined, found: false }),
  { value: request, found: true });
/** `request` with `value` put at `path`, a copy. */
const putAt = (request, path, value) => {
  const copy = structuredClone(request);
  const keys = keysOf(path);
  keys.slice(0, -1).reduce((parent, key) => (parent[key] ??= {}), copy)[keys.at(-1)] = value;
  return copy;
};

/**
 * The installed drawing's draft request, as Save projects it, writes every
 * value the derived request states. A cross-mode derivation leaves the new
 * mode's settings out (Python carries no Linear setting of a Circular
 * source); the drawing installs them as Load installs a request that leaves
 * them out, and writes them. `leftToLoad` lists the paths the derived request
 * leaves out and the drawing writes with an effect on the request; any other
 * difference fails.
 */
const assertRoundTrip = async (installed, { leftToLoad = [] } = {}) => {
  const { envelope } = await projectOf(installed.made);
  const compare = (expectedRequest) => compareCanonicalRenderRequests({
    expectedRequest, expectedResources: derivedResources(installed),
    actualRequest: envelope.renderRequest, actualResources: envelope.resources
  });
  let expected = installed.derived.renderRequest;
  const taken = [];
  let comparison = await compare(expected);
  while (comparison.differences.length && taken.length <= leftToLoad.length) {
    const { path } = comparison.differences[0];
    assert.equal(at(installed.derived.renderRequest, path).found, false, `the derived request states ${path}`);
    taken.push(path);
    expected = putAt(expected, path, at(envelope.renderRequest, path).value);
    comparison = await compare(expected);
  }
  assert.deepEqual(taken.sort(), [...leftToLoad].sort(), 'the paths the derived request leaves to Load');
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

// The Load readers read what Python derived: a request they would repair
// means Python and the Web readers disagree, so the region is not made.
test('(b) a derived request the Load readers would repair is refused, not installed with a loss', async () => {
  const source = circularSource([recordA]);
  const sent = await projectOf(source);
  const { renderRequest, resources } = sent.envelope;
  const { result: derived } = await deriveInPython({
    canonicalRequest: renderRequest, resources, regions: [{ recordKey: 'record-1', start: 1, end: 1200 }], adaptSizes: true
  });
  assert.equal(derived.status, 'derived', JSON.stringify(derived));
  // A GC content minimum above its maximum: the Web's config reader replaces the maximum.
  const disagreeing = putAt(derived.renderRequest, '$.diagramOptions.configOverrides.objects.gc_content.min_percent', 300);
  // A Load of the same request lists the loss.
  const loaded = await prepareStoredDrawing({ mode: derived.mode, renderRequest: disagreeing, resources, sourceFiles: new Map() });
  assert.notDeepEqual(loaded.repairs, []);
  await assert.rejects(prepareRegionDrawing({ derived: { ...derived, renderRequest: disagreeing }, sent }),
    (error) => error.code === 'RESULT_INVALID' && error.operation === 'deriveRegionDrawing');
});

test('(c) a region across the origin of a circular record is one window across it, kept as such by the card', async () => {
  // D-34: Python draws the span as one window from base 900 through the origin to base 100; the card that holds it
  // writes the same window back (its Start/End fields would read 900-100 as a reversed span, so it is not theirs).
  const source = withNewDrawingDefaults(state, createDrawingState('across-origin', 'linear'));
  source.inputType.value = 'gb';
  source.linearSeqs.splice(0, Infinity, createLinearSeq({ uid: 'mito', gb: state.inputFiles.registerFile(genbankFile('a.gb', [recordA])) }));
  Object.assign(source.form, { labels_mode: 'none' });
  const [{ recordKey }] = (await projectOf(source)).envelope.renderRequest.records;
  const installed = await deriveAndInstall(source, [{ recordKey, start: 900, end: 100 }]);
  assert.equal(installed.derived.mode, 'linear');
  const [record] = installed.derived.renderRequest.records;
  assert.deepEqual(record.windows.map(({ start, end, reverseComplement }) => ({ start, end, reverseComplement })),
    [{ start: 900, end: 100, reverseComplement: false }]);
  const [card] = installed.made.linearSeqs;
  assert.deepEqual([card.region_start, card.region_end, card.region_reverse], [null, null, false]);
  assert.deepEqual(card.windows.map(({ start, end, reverseComplement }) => ({ start, end, reverseComplement })),
    [{ start: 900, end: 100, reverseComplement: false }]);
  await assertRoundTrip(installed);
});

test('(c) a start after the end of a linear record is refused: no drawing, the tabs and History unchanged', async () => {
  circularSource([genbankRecord('recL', 1200, [[101, 400, 'L1']], 'linear')]);
  const before = listView();
  const created = [];
  const regions = createRegionDrawing({
    shownDrawingId: () => state.activeDrawingId.value,
    projectShownRequest: () => projectOf(state.activeDrawing()),
    deriveRegion: deriveInPython,
    newDrawingFromRegion: async (...args) => { created.push(args); return { status: 'ok', drawingId: 'never' }; }
  });
  const offer = await regions.deriveFromShown([{ recordKey: 'record-1', start: 900, end: 100 }], true);
  assert.deepEqual(offer, { status: 'refused', reason: 'START_AFTER_END', recordKey: 'record-1', start: 900, end: 100 });
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

// The drawing list (`drawing-list.js`) over the store and the app's History,
// with the region port the app passes (`makeRegionDrawing`, which reads with
// `prepareRegionDrawing` and makes with `installPreparedDrawing`). The other
// ports are stubs, not the app's wiring: the Playwright case covers that; a
// switch only changes the shown drawing.
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
  regionDrawing: makeRegionDrawing,
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
  projectShownRequest: () => projectOf(state.activeDrawing()),
  deriveRegion: deriveInPython,
  newDrawingFromRegion: drawingList.newDrawingFromRegion
});

test('(d) the drawing list with makeRegionDrawing: Create adds the region drawing as New does; Undo removes it, Redo brings it back', async () => {
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
  const ids = await featureIdsInPython((await projectOf(source)).envelope);
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

// Create pressed twice before the first one adds the drawing (a double click):
// one drawing, one History step.
test('(d) the drawing list with makeRegionDrawing: a second Create while one runs adds one drawing and one History step', async () => {
  circularSource([recordA]);
  const before = listView();
  const undoBefore = history.timeline().undo.length;
  const owner = regionDrawingsOf();
  owner.dialog.rows = [{ recordKey: 'record-1', caption: 'A', length: 1000, checked: true, kind: 'coordinates',
    start: '301', end: '700', featureId: '', features: null, featureReason: '' }];
  await owner.refresh();
  assert.equal(owner.dialog.status, 'derived');
  const outcomes = await Promise.all([owner.confirm(), owner.confirm()]);
  assert.equal(outcomes[1], null, 'the second Create adds nothing');
  assert.equal(state.drawingTabs.value.length, before.tabs.length + 1, 'one drawing');
  assert.equal(history.timeline().undo.length, undoBefore + 1, 'one History step');
  await history.undo();
  assert.deepEqual(state.drawingTabs.value.map(({ id }) => id), before.tabs.map(({ id }) => id));
  assert.equal(state.activeDrawingId.value, before.shown);
});

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
  const sent = await projectOf(await circularSourceWithLook());
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

// The Linear settings a Circular source does not have, which the Web writes
// as Load installs a Linear request that leaves them out (the config defaults
// and the Linear legend and title positions of New Linear drawing).
const LINEAR_SETTINGS_LEFT_TO_LOAD = [
  ...['labels.linear.scope', 'objects.definition.linear.show_replicon', 'objects.definition.linear.show_accession',
    'objects.definition.linear.show_length', 'labels.linear.placement', 'canvas.linear.align_center',
    'canvas.linear.keep_definition_left_aligned', 'canvas.linear.track_layout', 'canvas.linear.ruler_on_axis',
    'objects.scale.style', 'canvas.linear.normalize_length'].map((key) => `$.diagramOptions.configOverrides.${key}`),
  '$.diagramOptions.output', '$.diagramOptions.tracks'
];
// The installed Linear drawing writes the derived request, shown or not.
test('(a) round trip of the Linear region drawing of a Circular source', async () => {
  const installed = await deriveAndInstall(await circularSourceWithLook(), LOOK_REGION);
  await assertRoundTrip(installed, { leftToLoad: LINEAR_SETTINGS_LEFT_TO_LOAD });
  // Owner Q3 (R17-39, A): the derived request states no `output` (Python's
  // region derivation writes none, so Python draws the legend at the right),
  // yet the Web's new Linear drawing starts with legend and title at the
  // bottom and saves them so.
  assert.equal(at(installed.derived.renderRequest, '$.diagramOptions.output').found, false);
  assert.deepEqual({ ...installed.made.layoutPreferences.linear }, { legend: 'bottom', plotTitlePosition: 'bottom' });
  const { envelope } = await projectOf(installed.made);
  assert.deepEqual(at(envelope.renderRequest, '$.diagramOptions.output').value,
    { ...at(envelope.renderRequest, '$.diagramOptions.output').value, legend: 'bottom', plotTitlePosition: 'bottom' });
});

// A Linear drawing with one card whose file holds two records: the card draws both.
const recordB = genbankRecord('recB', 900, [[51, 300, 'B1'], [501, 800, 'B2']], 'linear');
test('(a) regions of two records of one file are one Linear drawing bound to that one file', async () => {
  const source = withNewDrawingDefaults(state, createDrawingState('two-records', 'linear'));
  source.inputType.value = 'gb';
  const fileRef = state.inputFiles.registerFile(genbankFile('two.gb', [recordA, recordB]));
  source.linearSeqs.splice(0, Infinity, createLinearSeq({ uid: 'two', gb: fileRef }));
  Object.assign(source.form, { labels_mode: 'none' });
  const sent = await projectOf(source);
  // The card's records, as the request names them (`<card uid>:<record number>`).
  const regions = [{ recordKey: 'two:1', start: 101, end: 900 }, { recordKey: 'two:2', start: 1, end: 600 }];
  const { value: installed, registered } = await countingRegistrations(() => deriveAndInstall(source, regions));
  assert.equal(installed.made.mode, 'linear');
  assert.deepEqual(installed.prepared.repairs, []);
  assert.deepEqual(registered, [], 'no new entry in the file store');
  assert.deepEqual(installed.made.linearSeqs.map((row) => row.gb?.fileId), [fileRef.fileId, fileRef.fileId],
    'both records bound to the source drawing\'s one Web file');
  assert.equal(new Set(Object.values(installed.derived.resourceSources)).size, 1, 'both read the one sent resource');
  // The source's card drew its two records on one row; the region drawing draws
  // each record on a row of its own.
  assert.deepEqual(installed.derived.renderRequest.records.map((record) => record.presentation.gridRow), [null, null]);
  await assertRoundTrip(installed);
});

// (e) A Linear record of 3 kb with a hand-set tick interval: Auto is 1 kb at
// 3 kb and 100 bp at 1.5 kb (`gbdraw/auto_sizes.py` `linear_tick_interval`),
// so the region drawing draws the interval at Auto and lists it as returned to
// Auto; with "Adapt sizes" off it keeps the hand-set interval.
const sizedRecord = genbankRecord('recS', 3000, [[101, 400, 'S1'], [1601, 2400, 'S2']], 'linear');
const sizedSource = (fileRef) => {
  const source = withNewDrawingDefaults(state, createDrawingState('sized', 'linear'));
  source.inputType.value = 'gb';
  source.linearSeqs.splice(0, Infinity, createLinearSeq({ uid: 'sized', gb: fileRef }));
  Object.assign(source.form, { labels_mode: 'none' });
  source.adv.scale_interval = 500;
  return source;
};
test('(e) a hand-set size whose Auto differs at the region\'s length returns to Auto; Adapt sizes off keeps it', async () => {
  const fileRef = state.inputFiles.registerFile(genbankFile('sized.gb', [sizedRecord]));
  const region = [{ recordKey: 'sized', start: 1, end: 1500 }];
  const adapted = await deriveAndInstall(sizedSource(fileRef), region, true, 'sized-adapted');
  assert.deepEqual(adapted.derived.returnedToAuto, [{ setting: 'objects.scale.interval', explicit: 500, auto: 100 }]);
  assert.equal(adapted.made.adv.scale_interval, null, 'the new drawing draws the interval at Auto');
  assert.equal(at(adapted.derived.renderRequest, '$.diagramOptions.configOverrides.objects.scale.interval').found, false);
  const kept = await deriveAndInstall(sizedSource(fileRef), region, false, 'sized-kept');
  assert.deepEqual(kept.derived.returnedToAuto, []);
  assert.equal(kept.made.adv.scale_interval, 500, 'Adapt sizes off keeps the hand-set interval');
});
