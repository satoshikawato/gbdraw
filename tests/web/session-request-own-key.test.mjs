// Circular request writer and Session reader: a table read by free text (a
// label mode, a label scope, a resource id) is read by own key, and an empty
// record list with a region does not dereference null.
import assert from 'node:assert/strict';
import test from 'node:test';

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
const { serializeActiveRenderFiles } = await import('../../gbdraw/web/js/services/config.js');
const { buildCanonicalRenderRequest, projectCanonicalSessionRequest } = await import(
  '../../gbdraw/web/js/services/session-request.js'
);
const { adoptCurrentSessionResources, createSessionResourceFileView } = await import(
  '../../gbdraw/web/js/services/session-resource-backing.js'
);

const bytes = new TextEncoder().encode('LOCUS c.gb\n//\n');
const genbank = {
  name: 'c.gb', type: 'text/plain', size: bytes.byteLength, lastModified: 1, arrayBuffer: async () => bytes.slice().buffer
};
const drawing = state.drawings.circular;

// Every case starts from one Circular GenBank file and its own form values.
const circularState = (form) => {
  state.mode.value = 'circular';
  state.cInputType.value = 'gb';
  state.files.c_gb = genbank;
  state.circularRecordList.value = [];
  Object.assign(drawing.form, {
    labels_mode: 'none', circular_record_selector: '', multi_record_canvas: false,
    circular_region_start: '', circular_region_end: '', ...form
  });
};
const buildWithFiles = async () => buildCanonicalRenderRequest({
  state, drawing, comparisonPlanSnapshot: null, recordDisplayRows: [],
  filesData: await serializeActiveRenderFiles('circular', state, drawing, null)
});
const emptyRecordFiles = () => ({ c_gb: state.files.c_gb, circularRecords: [] });

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
