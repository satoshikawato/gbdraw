// Circular request writer: a free-text label mode is read from an own-key
// table, and an empty record list with a region does not dereference null.
import assert from 'node:assert/strict';

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
const { buildCanonicalRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');

const genbank = (name) => {
  const bytes = new TextEncoder().encode(`LOCUS ${name}\n//\n`);
  return { name, type: 'text/plain', size: bytes.byteLength, lastModified: 1, arrayBuffer: async () => bytes.slice().buffer };
};

state.mode.value = 'circular';
state.cInputType.value = 'gb';
state.files.c_gb = genbank('c.gb');
const drawing = state.drawings.circular;

const build = async () => {
  const filesData = await serializeActiveRenderFiles('circular', state, drawing, null);
  return buildCanonicalRenderRequest({ state, drawing, filesData, comparisonPlanSnapshot: null, recordDisplayRows: [] });
};
drawing.form.labels_mode = 'toString';
const { renderRequest } = await build();
assert.equal(renderRequest.diagramOptions.configOverrides['labels.circular.scope'], 'none');

// An empty offered record list with a region: the region has no record to read.
Object.assign(drawing.form, {
  labels_mode: 'out', circular_record_selector: '', multi_record_canvas: false,
  circular_region_start: 3, circular_region_end: 10
});
state.circularRecordList.value = [];
const { renderRequest: cropped } = buildCanonicalRenderRequest({
  state, drawing, filesData: { c_gb: state.files.c_gb, circularRecords: [] }
});
assert.equal(cropped.records[0].region.selector, null);
assert.equal(cropped.records[0].region.start, 3);
