// Circular request writer and Session reader: a table read by free text (a
// label mode, a label scope, a resource id) is read by own key, and an empty
// record list with a region does not dereference null.
import assert from 'node:assert/strict';
import test from 'node:test';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {};

const { state } = await import('../../gbdraw/web/js/state.js');
const { bindDrawingInputFile } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { serializeActiveRenderFiles, installShownDrawing } = await import('../../gbdraw/web/js/services/config.js');
const {
  buildCanonicalRenderRequest, compareCanonicalRenderRequests, decodeCanonicalResourceText,
  projectCanonicalSessionRequest
} = await import('../../gbdraw/web/js/services/session-request.js');
const { adoptCurrentSessionResources, createSessionResourceFileView } = await import(
  '../../gbdraw/web/js/services/session-resource-backing.js'
);

const bytes = new TextEncoder().encode('LOCUS c.gb\n//\n');
const genbank = {
  name: 'c.gb', type: 'text/plain', size: bytes.byteLength, lastModified: 1, arrayBuffer: async () => bytes.slice().buffer
};
const drawing = modeDrawing(state, 'circular');

// Every case starts from one Circular GenBank file and its own form values.
const circularState = (form) => {
  installShownDrawing('circular');
  state.cInputType.value = 'gb';
  bindDrawingInputFile(state.inputFiles, drawing.slices.inputs, 'c_gb', genbank);
  drawing.circularRecordList.value = [];
  Object.assign(drawing.form, {
    labels_mode: 'none', circular_record_selector: '', multi_record_canvas: false,
    circular_region_start: '', circular_region_end: '', ...form
  });
};
const buildWithFiles = async () => buildCanonicalRenderRequest({
  state, drawing, comparisonPlanSnapshot: null, recordDisplayRows: [],
  filesData: await serializeActiveRenderFiles('circular', state, drawing, null)
});
const emptyRecordFiles = () => ({ c_gb: state.inputFiles.fileOf(drawing.files.c_gb), circularRecords: [] });

test('OV-563: a Labels mode outside the table writes the scope none', async () => {
  circularState({ labels_mode: 'toString' });
  const { renderRequest } = await buildWithFiles();
  assert.equal(renderRequest.diagramOptions.configOverrides['labels.circular.scope'], 'none');
});

test('OV-562: a region with an empty offered record list has no selector', () => {
  circularState({ labels_mode: 'out', circular_region_start: 3, circular_region_end: 10 });
  const { renderRequest } = buildCanonicalRenderRequest({ state, drawing, filesData: emptyRecordFiles() });
  assert.equal(renderRequest.records[0].region.selector, null);
  assert.equal(renderRequest.records[0].region.start, 3);
});

test('OV-564: a committed table whose resource id is an Object.prototype member name is written', () => {
  circularState({});
  const committedColors = createSessionResourceFileView(adoptCurrentSessionResources({
    toString: {
      kind: 'colors-default-colors-file', name: 'd.tsv', type: 'text/tab-separated-values',
      size: 12, lastModified: 0, encoding: 'base64', data: Buffer.from('CDS\t#ff0000\n').toString('base64')
    }
  }), 'toString');
  const published = buildCanonicalRenderRequest({
    state: { ...state, canonicalPublicationFiles: { d_color: committedColors } },
    drawing, filesData: emptyRecordFiles()
  });
  assert.equal(published.renderRequest.diagramOptions.colors.defaultColorsFile.resourceId, 'toString');
  assert.equal(Object.hasOwn(published.resources, 'toString'), true);
});

test('OV-565: a committed Circular label scope outside the table shows None', async () => {
  circularState({});
  const built = await buildWithFiles();
  for (const [scope, labelsMode] of [
    ['outer', 'out'], ['both', 'both'], ['none', 'none'], ['xyz', 'none'], ['toString', 'none'], ['valueOf', 'none']
  ]) {
    const session = JSON.parse(JSON.stringify(built));
    session.renderRequest.diagramOptions.configOverrides['labels.circular.scope'] = scope;
    assert.equal(projectCanonicalSessionRequest(session).config.form.labels_mode, labelsMode, scope);
  }
});

test('OV-613: a resource id that names an Object.prototype member and no resource is missing', async () => {
  circularState({});
  const built = await buildWithFiles();
  for (const resourceId of ['__proto__', 'constructor', 'toString']) {
    const session = JSON.parse(JSON.stringify(built));
    session.renderRequest.records[0].source.resourceId = resourceId;
    assert.throws(
      () => projectCanonicalSessionRequest(session),
      { message: `Missing canonical resource: ${resourceId}` }, resourceId
    );
    assert.throws(
      () => decodeCanonicalResourceText(session.resources, resourceId),
      { message: `Missing canonical resource: ${resourceId}` }, resourceId
    );
  }
  // A resource the Session stores under such an id is read (JSON.parse writes it as an own key).
  const stored = JSON.parse('{"__proto__":{"kind":"genbank","encoding":"base64","data":"TE9DVVM="}}');
  assert.equal(decodeCanonicalResourceText(stored, '__proto__'), 'LOCUS');
});

test('OV-614: a request compared with a resource table that lacks its resource names the resource', async () => {
  for (const resourceId of ['absent', '__proto__', 'constructor']) {
    const request = { diagramOptions: { depthTracks: [{ resourceId }] } };
    await assert.rejects(
      compareCanonicalRenderRequests({
        expectedRequest: request, expectedResources: {}, actualRequest: request, actualResources: {}
      }),
      { message: `Canonical request resource '${resourceId}' is missing.` }, resourceId
    );
  }
});
