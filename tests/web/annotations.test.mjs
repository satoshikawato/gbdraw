import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-annotations-'));
await cp(join(process.cwd(), 'gbdraw', 'web', 'js'), join(tempRoot, 'js'), { recursive: true });
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
const load = (path) => import(pathToFileURL(join(tempRoot, 'js', path)));

const { createAnnotationSet, normalizeAnnotationSets } = await load('app/annotations/state.js');
const {
  annotationRecordSelector,
  annotationRecordSelectorFromValue,
  annotationRecordSelectorValue,
  coordinateTarget,
  featureTarget,
  parseAnnotationRecordSelectorValue
} = await load('app/annotations/target-actions.js');
const { buildAnnotationRecordCatalog } = await load('app/annotations/record-catalog.js');
const {
  annotationRecordOptions,
  reconcileAnnotationRecordBindings,
  setAnnotationRecordValue
} = await load('app/annotations/record-selector.js');
const { validateAnnotationRecordTargets } = await load('app/annotations/validation.js');
const { encodeAnnotationTable, parseAnnotationTable } = await load('app/annotations/table-codec.js');
const { createAnnotationEditor } = await load('app/annotations.js');
const { buildLinearTrackSlotSpec, normalizeLinearTrackSlots } = await load('app/linear-track-slots.js');
const { buildCircularTrackSlotSpec, normalizeCircularTrackSlots } = await load('app/circular-track-slots.js');
const { hasLinearComparisonIntent } = await load('app/linear-comparisons.js');

const set = createAnnotationSet({
  id: 'review',
  annotations: [
    { id: 'coords', target: coordinateTarget({ start: 10, end: 20 }), label: 'Window', mark: 'band' },
    { id: 'gene', target: featureTarget({ selector: 'locus_tag=ABC_1' }), label: 'Gene', mark: 'bracket' }
  ]
});
const normalized = normalizeAnnotationSets([set]);
assert.equal(normalized[0].annotations[0].target.coordinateSpace, 'source');
assert.deepEqual(normalized[0].annotations[1].target.selectors[0], { key: 'locus_tag', value: 'ABC_1' });

const table = encodeAnnotationTable(normalized);
const restored = parseAnnotationTable(table);
assert.equal(restored[0].id, 'review');
assert.equal(restored[0].annotations[0].target.start, 10);
assert.equal(restored[0].annotations[1].target.selectors[0].value, 'ABC_1');

assert.equal(annotationRecordSelectorValue(annotationRecordSelectorFromValue('')), '');
assert.equal(annotationRecordSelectorValue(annotationRecordSelectorFromValue('RecA')), 'RecA');
assert.equal(annotationRecordSelectorValue(annotationRecordSelectorFromValue('#2')), '#2');
assert.equal(annotationRecordSelectorValue(annotationRecordSelectorFromValue('# 2')), '#2');
assert.equal(annotationRecordSelectorFromValue('NULL'), null);
assert.throws(() => annotationRecordSelectorFromValue('#0'), /must be >= 1/);
assert.throws(() => annotationRecordSelectorFromValue('#record'), /Invalid record selector/);
assert.deepEqual(parseAnnotationRecordSelectorValue('RecA').selector, { kind: 'recordId', value: 'RecA' });
assert.deepEqual(annotationRecordSelector(null, 0), { kind: 'recordIndex', index: 0 });
assert.equal(annotationRecordSelector(null, -1), null);
assert.equal(annotationRecordSelector(null, 'not-a-number'), null);

const state = {
  annotationSets: [],
  selectedFeatures: { value: [] },
  adv: { circular_track_slots: [], linear_track_slots: [] }
};
const editor = createAnnotationEditor({ state });
const created = editor.addAnnotationSet('review');
editor.addCoordinateAnnotation(created, { start: 5, end: 8 });
created.annotations[0].target.record = { kind: 'recordIndex', index: 1 };
editor.setAnnotationTargetKind(created.annotations[0], 'featureSpan');
assert.deepEqual(created.annotations[0].target.record, { kind: 'recordIndex', index: 1 });
editor.setAnnotationTargetKind(created.annotations[0], 'coordinateSpan');
assert.deepEqual(created.annotations[0].target.record, { kind: 'recordIndex', index: 1 });
const copied = editor.duplicateAnnotationSet(created);
assert.equal(state.annotationSets.length, 2);
assert.notEqual(copied.id, created.id);
editor.removeAnnotationSet(copied);
assert.equal(state.annotationSets.length, 1);

const linearCatalog = buildAnnotationRecordCatalog({
  mode: 'linear',
  linearSources: [
    {
      sourceKey: 'source-a',
      hasInput: true,
      status: 'ready',
      selector: '',
      records: [{ selector: '#1', recordId: 'RecA', recordLength: 100 }]
    },
    {
      sourceKey: 'source-b',
      hasInput: true,
      status: 'ready',
      selector: '',
      records: [{ selector: '#1', recordId: 'RecB', recordLength: 200 }]
    }
  ]
});
assert.deepEqual(linearCatalog.records.map((record) => record.value), ['RecA', 'RecB']);
assert.match(linearCatalog.records[0].label, /#1 · RecA · 100 bp/);

const duplicateCatalog = buildAnnotationRecordCatalog({
  mode: 'linear',
  linearSources: [
    { sourceKey: 'dup-a', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'dup', recordLength: 10 }] },
    { sourceKey: 'dup-b', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'dup', recordLength: 20 }] }
  ]
});
assert.deepEqual(duplicateCatalog.records.map((record) => record.value), ['#1', '#2']);

const targetSet = createAnnotationSet({
  id: 'targets',
  annotations: [{ id: 'window', target: coordinateTarget({ start: 1, end: 5 }), label: '', mark: 'band' }]
});
assert.match(validateAnnotationRecordTargets([targetSet], linearCatalog), /targets\/window/);
const featureTargetSet = createAnnotationSet({
  id: 'features',
  annotations: [{ id: 'gene', target: featureTarget({ selector: 'locus_tag=ABC_1' }), label: '', mark: 'bracket' }]
});
assert.match(validateAnnotationRecordTargets([featureTargetSet], linearCatalog), /features\/gene/);
assert.deepEqual(
  annotationRecordOptions(linearCatalog, targetSet.annotations[0]).map((option) => option.value),
  ['', linearCatalog.records[0].key, linearCatalog.records[1].key]
);
setAnnotationRecordValue(linearCatalog, targetSet.annotations[0], linearCatalog.records[1].key);
assert.deepEqual(targetSet.annotations[0].target.record, { kind: 'recordId', value: 'RecB' });
assert.equal(validateAnnotationRecordTargets([targetSet], linearCatalog), '');
assert.equal(encodeAnnotationTable([targetSet]).trim().split('\n')[1].split('\t')[3], 'RecB');
delete targetSet.annotations[0].metadata._gbdraw_web_target_record_key;
targetSet.annotations[0].target.record = { kind: 'recordIndex', index: 2 };
assert.match(validateAnnotationRecordTargets([targetSet], linearCatalog), /out of range/);
targetSet.annotations[0].target.record = { kind: 'recordId', value: 'missing' };
assert.match(validateAnnotationRecordTargets([targetSet], linearCatalog), /not available/);
setAnnotationRecordValue(duplicateCatalog, targetSet.annotations[0], duplicateCatalog.records[0].key);
assert.equal(validateAnnotationRecordTargets([targetSet], duplicateCatalog), '');
targetSet.annotations[0].target.record = { kind: 'recordIndex', index: -1 };
delete targetSet.annotations[0].metadata._gbdraw_web_target_record_key;
assert.match(
  validateAnnotationRecordTargets([targetSet], linearCatalog),
  /valid target record/
);

const automaticSourceCatalog = buildAnnotationRecordCatalog({
  mode: 'linear',
  linearSources: [{
    sourceKey: 'automatic',
    hasInput: true,
    status: 'ready',
    selector: '',
    records: [{ recordId: 'RecA' }, { recordId: 'RecB' }]
  }]
});
assert.deepEqual(automaticSourceCatalog.records.map((record) => record.recordId), ['RecA', 'RecB']);
assert.deepEqual(automaticSourceCatalog.issues, []);

const selectedSourceCatalog = buildAnnotationRecordCatalog({
  mode: 'linear',
  linearSources: [{
    sourceKey: 'selected',
    hasInput: true,
    status: 'ready',
    selector: '#2',
    records: [{ recordId: 'RecA' }, { recordId: 'RecB' }]
  }]
});
assert.deepEqual(selectedSourceCatalog.records.map((record) => record.recordId), ['RecB']);

const gbComparisonCatalog = buildAnnotationRecordCatalog({
  mode: 'linear',
  inputType: 'gb',
  loadComparison: true,
  linearSources: [
    { sourceKey: 'gb-a', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'A1' }, { recordId: 'A2' }] },
    { sourceKey: 'gb-b', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'B1' }, { recordId: 'B2' }] }
  ]
});
assert.deepEqual(gbComparisonCatalog.records.map((record) => record.recordId), ['A1', 'B1']);
const singleGbComparisonCatalog = buildAnnotationRecordCatalog({
  mode: 'linear', inputType: 'gb', loadComparison: true,
  linearSources: [{ sourceKey: 'gb-only', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'A1' }, { recordId: 'A2' }] }]
});
assert.deepEqual(singleGbComparisonCatalog.records.map((record) => record.recordId), ['A1', 'A2']);
const gffComparisonCatalog = buildAnnotationRecordCatalog({
  mode: 'linear', inputType: 'gff', loadComparison: true,
  linearSources: [{ sourceKey: 'gff-only', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'G1' }, { recordId: 'G2' }] }]
});
assert.deepEqual(gffComparisonCatalog.records.map((record) => record.recordId), ['G1']);
assert.equal(hasLinearComparisonIntent({ sequences: [{}, {}], blastSource: 'losat' }), true);
assert.equal(hasLinearComparisonIntent({ sequences: [{}], blastSource: 'losat' }), false);
assert.equal(hasLinearComparisonIntent({ sequences: [{ blast: {} }, {}], blastSource: 'upload' }), true);
assert.equal(hasLinearComparisonIntent({ layoutEnabled: true, comparisons: [{}], sequences: [] }), true);

const duplicateTargetSet = createAnnotationSet({
  id: 'duplicates',
  annotations: [{ id: 'same', target: coordinateTarget({ start: 1, end: 2 }), mark: 'band' }]
});
setAnnotationRecordValue(duplicateCatalog, duplicateTargetSet.annotations[0], duplicateCatalog.records[1].key);
assert.deepEqual(duplicateTargetSet.annotations[0].target.record, { kind: 'recordIndex', index: 1 });
const reorderedDuplicateCatalog = buildAnnotationRecordCatalog({
  mode: 'linear',
  linearSources: [
    { sourceKey: 'dup-b', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'dup', recordLength: 20 }] },
    { sourceKey: 'dup-a', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'dup', recordLength: 10 }] }
  ]
});
reconcileAnnotationRecordBindings([duplicateTargetSet], reorderedDuplicateCatalog);
assert.deepEqual(duplicateTargetSet.annotations[0].target.record, { kind: 'recordIndex', index: 0 });
assert.equal(validateAnnotationRecordTargets([duplicateTargetSet], reorderedDuplicateCatalog), '');
const replacedSourceCatalog = buildAnnotationRecordCatalog({
  mode: 'linear',
  linearSources: [{ sourceKey: 'replacement', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'dup' }] }]
});
assert.match(validateAnnotationRecordTargets([duplicateTargetSet], replacedSourceCatalog), /no longer available/);

targetSet.annotations[0].target.record = { kind: 'recordId', value: 'RecA' };
targetSet.annotations[0].metadata = {};
const circularMultiOutputCatalog = buildAnnotationRecordCatalog({
  mode: 'circular',
  multiRecordCanvas: false,
  circularSource: {
    sourceKey: 'circular', hasInput: true, status: 'ready',
    records: [{ record_id: 'RecA' }, { record_id: 'RecB' }]
  }
});
assert.match(validateAnnotationRecordTargets([targetSet], circularMultiOutputCatalog), /enable Multi-record canvas/);

const circularSingleCatalog = buildAnnotationRecordCatalog({
  mode: 'circular',
  circularSource: { sourceKey: 'single', hasInput: true, status: 'ready', records: [{ record_id: 'RecA' }] }
});
targetSet.annotations[0].target.record = null;
assert.equal(validateAnnotationRecordTargets([targetSet], circularSingleCatalog), '');

const tableLines = table.trimEnd().split('\n');
const tableHeader = tableLines[0].split('\t');
const recordColumn = tableHeader.indexOf('record');
const spacedIndexRow = tableLines[1].split('\t');
spacedIndexRow[recordColumn] = '# 2';
const spacedIndexSets = parseAnnotationTable(`${tableLines[0]}\n${spacedIndexRow.join('\t')}\n`);
assert.deepEqual(spacedIndexSets[0].annotations[0].target.record, { kind: 'recordIndex', index: 1 });
const nullRecordRow = tableLines[1].split('\t');
nullRecordRow[recordColumn] = 'NULL';
const nullRecordSets = parseAnnotationTable(`${tableLines[0]}\n${nullRecordRow.join('\t')}\n`);
assert.equal(nullRecordSets[0].annotations[0].target.record, null);
assert.match(validateAnnotationRecordTargets(nullRecordSets, linearCatalog), /Choose a target record/);
const invalidIndexRow = tableLines[1].split('\t');
invalidIndexRow[recordColumn] = '#0';
assert.throws(
  () => parseAnnotationTable(`${tableLines[0]}\n${invalidIndexRow.join('\t')}\n`),
  /row 2, column 'record'/
);

const linearSlot = normalizeLinearTrackSlots([{
  id: 'notes', renderer: 'annotations', side: 'overlay',
  params: { set_id: 'review', anchor_slot: 'features', layer: 'underlay', show_labels: false }
}])[0];
assert.equal(linearSlot.side, 'overlay');
assert.match(buildLinearTrackSlotSpec(linearSlot), /set_id=review/);
assert.match(buildLinearTrackSlotSpec(linearSlot), /anchor_slot=features/);

const circularSlot = normalizeCircularTrackSlots([{
  id: 'notes', renderer: 'annotations', side: 'outside',
  params: { set_id: 'review', overflow: 'compress' }
}])[0];
assert.equal(circularSlot.side, 'outside');
assert.match(buildCircularTrackSlotSpec(circularSlot), /set_id=review/);
assert.match(buildCircularTrackSlotSpec(circularSlot), /overflow=compress/);
