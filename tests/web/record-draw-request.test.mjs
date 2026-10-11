// Record selection (DESIGN 9, review F3): a Linear request with an OFF card
// reads one drawn list through every consumer: the record files, the request
// records and rows, the comparison indexes, the record display rows, and the
// annotations bound to the OFF record, which wait in the draft.
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
const { resolveLinearComparisonPlan } = await import('../../gbdraw/web/js/services/linear-comparisons.js');

const genbank = (name) => {
  const bytes = new TextEncoder().encode(`LOCUS ${name}\n//\n`);
  return { name, type: 'text/plain', size: bytes.byteLength, lastModified: 1, arrayBuffer: async () => bytes.slice().buffer };
};
const blast = genbank('a-to-c.tsv');

state.mode.value = 'linear';
state.lInputType.value = 'gb';
const drawing = state.drawings.linear;
state.linearSeqs.splice(0, Infinity, ...['a', 'b', 'c'].map((uid) => ({
  uid, gb: genbank(`${uid}.gb`), gff: null, fasta: null, depth: null, losat_gencode: 1, definition: '',
  record_subtitle: '', file_definition: '', file_subtitle: '', inferred_definition: '', region_record_id: '',
  region_start: null, region_end: null, region_reverse: false
})));
drawing.recordsOff.splice(0, Infinity, 'b');
drawing.linearRecordLayoutEnabled.value = true;
drawing.linearRecordRows.splice(0, Infinity, { uid: 'a', row: 1 }, { uid: 'b', row: 2 }, { uid: 'c', row: 3 });
Object.assign(drawing.linearComparisonPlan, {
  mode: 'selected',
  defaultSource: 'upload',
  edges: [{ id: 'edge-a-c', queryUid: 'a', subjectUid: 'c', included: true, fileActive: true, losatFilenameActive: false,
    source: 'upload', file: blast, losatFilename: '' }]
});
drawing.recordDisplayDrafts.splice(0, Infinity, {
  sourceUid: 'c', selector: '#1', recordId: 'c', topologyOverride: null, startCoordinate: null,
  reverseComplementOverride: true, anchorIntent: null
});
const annotation = (id, binding) => ({
  id, target: { kind: 'coordinateSpan', recordSelector: null, start: 1, end: 5, strand: null },
  label: id, mark: 'highlight', lane: null, style: null, legendLabel: null,
  metadata: { _gbdraw_web_target_record_key: binding }
});
drawing.annotationSets.splice(0, Infinity, {
  id: 'set', label: 'Set', annotations: [annotation('on-c', 'catalog-c'), annotation('off-b', 'catalog-b')]
});

const comparisonPlanSnapshot = resolveLinearComparisonPlan({
  plan: drawing.linearComparisonPlan,
  sequences: state.linearSeqs,
  recordsOff: drawing.recordsOff,
  layout: drawing.linearRecordRows
});
const filesData = await serializeActiveRenderFiles('linear', state, drawing, comparisonPlanSnapshot);
assert.deepEqual(filesData.linearSeqs.map((seq) => seq.uid), ['a', 'c']);
const recordDisplayRows = ['a', 'b', 'c'].map((uid) => ({
  key: JSON.stringify([uid, '#1']), scope: 'linear', sourceUid: uid, selector: '#1', recordId: uid,
  recordLength: 100, detectedTopology: 'linear', reverse: false, cropped: false
}));
const { renderRequest } = buildCanonicalRenderRequest({
  state, drawing, filesData, comparisonPlanSnapshot, recordDisplayRows,
  omittedAnnotationRecordKeys: ['catalog-b']
});

assert.deepEqual(renderRequest.records.map((record) => record.recordKey), ['a', 'c']);
// The OFF record's row stays in the draft; each ON record keeps its own row.
assert.deepEqual(renderRequest.records.map((record) => record.presentation.gridRow), [1, 3]);
assert.deepEqual(drawing.linearRecordRows.map((row) => row.uid), ['a', 'b', 'c']);
// The Selected edge a -> c compares positions 0 and 1 of the drawn records.
assert.deepEqual(
  renderRequest.comparisons.map((comparison) => [comparison.queryRecordIndex, comparison.subjectRecordIndex]),
  [[0, 1]]
);
// The record display row of c reaches the request on c's record, not on a's.
assert.deepEqual(renderRequest.records.map((record) => record.presentation.reverseComplement), [false, true]);
// The annotation bound to the OFF record waits in the draft.
const requested = renderRequest.diagramOptions.annotations.sets[0].annotations.map((item) => item.id);
assert.deepEqual(requested, ['on-c']);
assert.deepEqual(drawing.annotationSets[0].annotations.map((item) => item.id), ['on-c', 'off-b']);

// Review F7: an error of a drawn card names its card number, OFF cards counted.
// Card 1 (a) is OFF; card 3 (c) reads two records and has a crop without a record.
drawing.recordsOff.splice(0, Infinity, 'a');
const cropped = state.linearSeqs.find((seq) => seq.uid === 'c');
Object.assign(cropped, { region_start: 1, region_end: 10 });
await assert.rejects(
  serializeActiveRenderFiles('linear', state, drawing, {
    comparisonPlan: { mode: 'none', edges: [] },
    linearRecordCatalog: {
      mode: 'linear', status: 'ready', issues: [],
      records: [{ sourceIndex: 0, localIndex: 0 }, { sourceIndex: 1, localIndex: 0 }, { sourceIndex: 1, localIndex: 1 }]
    }
  }),
  { code: 'REGION_INVALID', context: { inputOrdinal: 3, reason: 'SELECT_RECORD_FOR_REGION' } }
);
Object.assign(cropped, { region_start: null, region_end: null });
await assert.rejects(
  serializeActiveRenderFiles('linear', state, drawing, {
    comparisonPlan: { mode: 'none', edges: [] },
    linearRecordCatalog: { mode: 'linear', status: 'ready', issues: [], records: [{ sourceIndex: 0, localIndex: 0 }] }
  }),
  { code: 'NO_RECORDS', context: { inputOrdinal: 3 } }
);
console.log('record draw request tests passed');
