import assert from 'node:assert/strict';
import test from 'node:test';
import { readFile, cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-annotations-'));
await cp(join(process.cwd(), 'gbdraw', 'web', 'js'), join(tempRoot, 'js'), { recursive: true });
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
const load = (path) => import(pathToFileURL(join(tempRoot, 'js', path)));

const {
  annotationOptionsPayload, createAnnotationSet, draftAnnotationSetsOfRequest, normalizeAnnotationSets
} = await load('services/annotation-state.js');
const {
  annotationRecordSelector,
  annotationRecordSelectorFromValue,
  annotationRecordSelectorValue,
  coordinateTarget,
  featureTarget,
  featureTargetsFromSelection,
  parseAnnotationRecordSelectorValue
} = await load('app/annotations/target-actions.js');
const { buildAnnotationRecordCatalog } = await load('app/annotations/record-catalog.js');
const { resolveCircularRequestRecordSet } = await load('services/record-options.js');
const {
  annotationRecordOptions,
  reconcileAnnotationRecordBindings,
  setAnnotationRecordValue
} = await load('app/annotations/record-selector.js');
const { validateAnnotationRecordTargets } = await load('app/annotations/validation.js');
const {
  encodeAnnotationTable, encodeAnnotationTableWithNotice, parseAnnotationTable, parseAnnotationTableWithNotice
} = await load('app/annotations/table-codec.js');
const { createAnnotationEditor } = await load('app/annotations.js');
const { buildLinearTrackSlotSpec, normalizeLinearTrackSlots } = await load('app/linear-track-slots.js');
const { buildCircularTrackSlotSpec, normalizeCircularTrackSlots } = await load('app/circular-track-slots.js');
const {
  createLinearComparisonEdge,
  hasLinearComparisonIntent,
  resolveLinearComparisonPlan
} = await load('services/linear-comparisons.js');

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

test('malformed imported feature selectors are rejected before replacing the draft', () => {
  for (const selector of [';', 'locus_tag=']) {
    assert.throws(() => parseAnnotationTable(`set_id\tid\tmark\tfeature_selector\ns\ta\thighlight\t${selector}\n`), /feature_selector requires/);
  }
});

test('invalid annotation coordinates cannot be silently clamped or truncated', () => {
  for (const start of ['abc', '0', '-10', '1.5', '', 'Infinity']) {
    assert.throws(() => parseAnnotationTable(`set_id\tid\tmark\tstart\tend\ns\ta\thighlight\t${start}\t30\n`), /positive integers/);
    const invalid = [{ id: 's', annotations: [{ id: 'a', target: coordinateTarget({ start, end: 30 }) }] }];
    assert.deepEqual(validateAnnotationRecordTargets(invalid, { records: [] }), { code: 'ANNOTATION_TARGET', context: { reason: 'POSITIVE_INTEGER' } });
  }
  assert.equal(parseAnnotationTable('set_id\tid\tmark\tstart\tend\ns\ta\thighlight\t10\t30\n')[0].annotations[0].target.start, 10);
});

test('TSV cannot silently replace invalid known annotation choices', () => {
  for (const [column, value] of [['coordinate_space', 'typo'], ['lane', '-1'], ['lane', '1.5']]) {
    assert.throws(() => parseAnnotationTable(`set_id\tid\tmark\tstart\tend\t${column}\ns\ta\tband\t1\t8\t${value}\n`));
  }
  const local = parseAnnotationTable('set_id\tid\tmark\tstart\tend\tcoordinate_space\ns\ta\tBAND\t1\t8\tLOCAL\n')[0].annotations[0];
  assert.equal(local.mark, 'band');
  assert.equal(local.target.coordinateSpace, 'local');
});

test('row style edits do not recolor sibling annotations or their inherited default', () => {
  const state = { annotationSets: [], adv: {} };
  const editor = createAnnotationEditor({ state });
  const set = editor.addAnnotationSet();
  const first = editor.addCoordinateAnnotation(set);
  const second = editor.addCoordinateAnnotation(set);
  const inherited = set.defaultStyle.fill;
  editor.setAnnotationStyle(set, first, 'fill', '#ff0000');
  assert.equal(first.style.fill, '#ff0000');
  assert.equal(second.style, null);
  assert.equal(set.defaultStyle.fill, inherited);
});

test('explicit no-fill survives TSV while omitted fill keeps the Web default', () => {
  const draft = [createAnnotationSet({
    id: 'no-fill', defaultStyle: { fill: null },
    annotations: [
      { id: 'inherited', target: coordinateTarget({ start: 5, end: 8 }), mark: 'band' },
      { id: 'override', target: featureTarget({ selector: 'locus_tag=ABC_1' }),
        mark: 'bracket', style: { fill: null } }
    ]
  })];
  const reimported = parseAnnotationTable(encodeAnnotationTable(draft));
  assert.deepEqual(reimported[0].annotations.map((item) => item.style.fill), [null, null]);
  const omittedFill = parseAnnotationTable('set_id\tid\tmark\tstart\tend\tstroke\ns\ta\tband\t1\t5\t#123456\n');
  assert.equal(omittedFill[0].annotations[0].style.fill, '#94a3b8');
});

test('download reads the current rows without mutating state or resolving a catalog', async () => {
  const draftState = {
    annotationSets: [],
    results: [{ content: '<svg/>', annotations: 'committed, older data' }],
    selectedResultIndex: 0, selectedAnnotation: { setId: 'review', id: 'coords' },
    selectedFeatures: [{ id: 'feature' }], zoom: 2, canvasPan: { x: 5, y: 8 },
    history: ['before', 'after']
  };
  const draftEditor = createAnnotationEditor({
    state: draftState, getRecordCatalog: () => { throw new Error('Download must not resolve records'); }
  });
  const calls = [];
  let blob;
  const original = { document: globalThis.document, create: URL.createObjectURL, revoke: URL.revokeObjectURL };
  globalThis.document = {
    body: {
      appendChild(link) { link.parentNode = this; calls.push('append'); },
      removeChild() { calls.push('remove'); }
    },
    createElement(tag) {
      assert.equal(tag, 'a');
      return { click() { calls.push(['click', this.download, this.href]); } };
    }
  };
  URL.createObjectURL = (value) => { blob = value; calls.push('blob'); return 'blob:annotation'; };
  URL.revokeObjectURL = (url) => calls.push(['revoke', url]);
  try {
    assert.equal(draftEditor.canDownloadAnnotationTable(), false);
    draftEditor.downloadAnnotationTable();
    draftState.annotationSets.push(createAnnotationSet());
    assert.equal(draftEditor.canDownloadAnnotationTable(), false);
    draftEditor.downloadAnnotationTable();
    assert.deepEqual(calls, []);
    draftState.annotationSets.push(createAnnotationSet({ id: 'current', annotations: [
      { id: 'id-bound', target: coordinateTarget({ recordId: 'unique', start: 7, end: 20 }), label: 'Current draft', mark: 'band' },
      { id: 'index-bound', target: featureTarget({ recordIndex: 1, selector: 'locus_tag=ABC_1' }), mark: 'bracket', style: { fill: null } }
    ] }));
    assert.equal(draftEditor.canDownloadAnnotationTable(), true);
    const before = structuredClone(draftState);
    draftEditor.downloadAnnotationTable();
    assert.deepEqual(draftState, before);
    assert.deepEqual(calls, ['blob', 'append', ['click', 'annotations.tsv', 'blob:annotation'], 'remove', ['revoke', 'blob:annotation']]);
    assert.equal(blob.type, 'text/tab-separated-values;charset=utf-8');
    const text = await blob.text();
    assert.equal(text, encodeAnnotationTable(draftState.annotationSets));
    const rows = parseAnnotationTable(text)[0].annotations;
    assert.equal(rows[0].label, 'Current draft');
    assert.deepEqual(rows.map((item) => item.target.record), [
      { kind: 'recordId', value: 'unique' }, { kind: 'recordIndex', index: 1 }
    ]);
  } finally {
    if (original.document === undefined) delete globalThis.document;
    else globalThis.document = original.document;
    URL.createObjectURL = original.create;
    URL.revokeObjectURL = original.revoke;
  }
});

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
test('delete and re-add keeps coordinate and selected-feature annotation IDs unique', () => {
  const state = { annotationSets: [], selectedFeatures: { value: [{
    scope: 'circular', record_key: 'k', biological_feature_id: 'f123'
  }] } };
  const actions = createAnnotationEditor({ state });
  const set = actions.addAnnotationSet();
  for (let i = 0; i < 3; i += 1) actions.addCoordinateAnnotation(set);
  actions.removeAnnotation(set, set.annotations[0]);
  actions.addCoordinateAnnotation(set);
  for (let i = 0; i < 3; i += 1) actions.addSelectedFeatures(set);
  actions.removeAnnotation(set, set.annotations[3]);
  actions.addSelectedFeatures(set);
  const ids = set.annotations.map((item) => item.id);
  assert.equal(new Set(ids).size, ids.length);
  assert.deepEqual(normalizeAnnotationSets([set])[0].annotations.map((item) => item.id), ids);
});

// Design Q4 PR-Q4-5 (OV-03): a selected feature is named by its source
// identity, never by a selector value that a crop or a copy changes.
test('selected annotations name each feature by its source identity', () => {
  const targets = featureTargetsFromSelection([
    { type: 'D-loop', scope: 'circular', record_key: 'mt', biological_feature_id: 'fcf4827e2', selector: { hash: 'fcf4827e2' } },
    { gene: 'duplicated', scope: 'linear', record_key: 'linear-seq-2', biological_feature_id: 'f1234~1',
      selector: { hash: 'f1234' } }
  ]);
  // Each target names the mode of the Result it was selected on (R2, OV-21).
  assert.deepEqual(targets, [
    { kind: 'featureIdentity', scope: 'circular', recordKey: 'mt', biologicalFeatureId: 'fcf4827e2',
      envelope: 'outer_bounds', circularPath: 'shortest' },
    { kind: 'featureIdentity', scope: 'linear', recordKey: 'linear-seq-2', biologicalFeatureId: 'f1234~1',
      envelope: 'outer_bounds', circularPath: 'shortest' }
  ]);
  // A Result without a feature catalog (a Session before 40) has no identities.
  assert.equal(featureTargetsFromSelection([{ selector: { hash: 'f1234' }, locus_tag: 'A_1' }]), null);
  const alerts = [];
  const originalWindow = globalThis.window;
  globalThis.window = { alert: (message) => alerts.push(message) };
  try {
    const legacy = createAnnotationEditor({ state: { annotationSets: [], selectedFeatures: [{ locus_tag: 'A_1' }] } });
    const set = legacy.addAnnotationSet();
    assert.deepEqual(legacy.addSelectedFeatures(set), []);
    assert.deepEqual(set.annotations, []);
    assert.deepEqual(alerts, ['Generate the diagram again to annotate the selected features.']);
  } finally {
    if (originalWindow === undefined) delete globalThis.window;
    else globalThis.window = originalWindow;
  }
});

// Design Q4 6.4: the table writes a selected feature in its current placement.
test('the TSV writes selected features by drawn record position and drawn hash, and counts the rest', () => {
  const identity = (recordKey, biologicalFeatureId) => ({
    kind: 'featureIdentity', scope: 'linear', recordKey, biologicalFeatureId, envelope: 'segments', circularPath: 'reverse'
  });
  const sets = [createAnnotationSet({ id: 'marks', annotations: [
    { id: 'drawn', target: identity('linear-seq-2', 'fb5977f81'), mark: 'band' },
    { id: 'cropped', target: identity('linear-seq-2', 'f0000000a'), mark: 'band' }
  ] })];
  const drawnPlacement = (target) => (target.biologicalFeatureId === 'fb5977f81' ? { recordIndex: 1, hash: 'f3b928d8c' } : null);
  const encoded = encodeAnnotationTableWithNotice(sets, { drawnPlacement });
  assert.equal(encoded.placedFeatureIdentityCount, 1);
  assert.equal(encoded.skippedFeatureIdentityCount, 1);
  const [header, row, ...rest] = encoded.text.trim().split('\n').map((line) => line.split('\t'));
  assert.deepEqual(rest, []);
  const cell = (name) => row[header.indexOf(name)];
  assert.deepEqual(['id', 'record', 'feature_selector', 'envelope', 'circular_path'].map(cell),
    ['drawn', '#2', 'hash=f3b928d8c', 'segments', 'reverse']);
  assert.deepEqual(parseAnnotationTable(encoded.text)[0].annotations[0].target, {
    kind: 'featureSpan', record: { kind: 'recordIndex', index: 1 },
    selectors: [{ key: 'hash', value: 'f3b928d8c' }], envelope: 'segments', circularPath: 'reverse'
  });
  // Without a placement (Run Info has none) no selected-feature row is written.
  assert.equal(encodeAnnotationTableWithNotice(sets).skippedFeatureIdentityCount, 2);
});

// A request carries the selected-feature targets of its own mode and records
// only; the draft keeps the others (design Q4 3.2, R2).
test('the request carries only the selected-feature targets of its records', () => {
  const identity = (recordKey, biologicalFeatureId) => ({
    kind: 'featureIdentity', scope: 'linear', recordKey, biologicalFeatureId
  });
  const sets = [createAnnotationSet({ id: 'marks', annotations: [
    { id: 'here', target: identity('linear-seq-1', 'fa'), mark: 'band' },
    { id: 'expanded', target: identity('linear-seq-2:3', 'fb'), mark: 'band' },
    { id: 'circular', target: identity('circular-A-#1', 'fc'), mark: 'band' },
    { id: 'coordinates', target: coordinateTarget({ start: 1, end: 5 }), mark: 'band' }
  ] })];
  const records = [
    { recordKey: 'linear-seq-1', cardinality: 'exactly_one' },
    { recordKey: 'linear-seq-2', cardinality: 'all' }
  ];
  const payload = annotationOptionsPayload(sets, 'linear', records);
  assert.deepEqual(payload.sets[0].annotations.map((item) => item.id), ['here', 'expanded', 'coordinates']);
  assert.equal(sets[0].annotations.length, 4);
});

// R2, OV-21: both modes can use the same record key for the same feature
// (`record-1` in a Gallery or Python Session). A selected-feature target names
// the mode it was made in, a request carries only the targets of its own mode
// (without the draft-only `scope`, as request schema 9 has it), and the other
// mode's target stays in the draft.
test('a selected-feature target reaches only the requests of its own mode (both ways)', () => {
  const MODES = ['circular', 'linear'];
  const records = [{ recordKey: 'record-1', cardinality: 'exactly_one' }];
  const target = (scope) => ({
    kind: 'featureIdentity', scope, recordKey: 'record-1', biologicalFeatureId: 'f1',
    envelope: 'outer_bounds', circularPath: 'shortest'
  });
  const sets = [createAnnotationSet({ id: 'marks', annotations: MODES.map((scope) => ({
    id: scope, target: target(scope), mark: 'highlight'
  })) })];
  const before = JSON.stringify(sets);
  MODES.forEach((mode) => {
    const payload = annotationOptionsPayload(sets, mode, records);
    assert.deepEqual(payload.sets[0].annotations.map((item) => [item.id, item.target]), [[mode, {
      kind: 'featureIdentity', recordKey: 'record-1', biologicalFeatureId: 'f1',
      envelope: 'outer_bounds', circularPath: 'shortest'
    }]]);
    // A Session written from a request (CLI, Python) gives each target the
    // request's mode, and the draft projects back onto the same request.
    const draft = draftAnnotationSetsOfRequest(payload.sets, mode);
    assert.deepEqual(draft[0].annotations.map((item) => item.target), [target(mode)]);
    assert.deepEqual(annotationOptionsPayload(draft, mode, records), payload);
    const other = mode === 'circular' ? 'linear' : 'circular';
    assert.deepEqual(annotationOptionsPayload(draft, other, records).sets[0].annotations, []);
  });
  assert.equal(JSON.stringify(sets), before);
});

// The draft owner rejects a selected-feature target that names no mode, so a
// path that writes one fails instead of dropping the annotation from every
// request.
test('a selected-feature target without a mode is invalid', () => {
  const set = (target) => [{ id: 'marks', annotations: [{ id: 'one', target }] }];
  const identity = { kind: 'featureIdentity', recordKey: 'record-1', biologicalFeatureId: 'f1' };
  assert.throws(() => normalizeAnnotationSets(set(identity)), /schema|INPUT_INVALID|invalid/i);
  assert.throws(() => normalizeAnnotationSets(set({ ...identity, scope: 'batch' })), /schema|INPUT_INVALID|invalid/i);
  assert.throws(() => normalizeAnnotationSets(set({ ...identity, scope: 'linear', recordKey: '' })),
    /schema|INPUT_INVALID|invalid/i);
  assert.equal(normalizeAnnotationSets(set({ ...identity, scope: 'linear' }))[0].annotations[0].target.scope, 'linear');
});

// The editor names a target's feature from the current Results only in the
// target's own mode: the Circular feature with the same identity is not it.
test('the editor finds a selected-feature target only among features of its mode', () => {
  const feature = { scope: 'circular', record_key: 'record-1', biological_feature_id: 'f1', type: 'CDS',
    locus_tag: 'LT_1', start: 0, end: 30, strand: 1, record_idx: 0, drawnSelector: { hash: 'f1' } };
  const actions = createAnnotationEditor({ state: {
    annotationSets: [], selectedFeatures: { value: [] },
    extractedFeatures: { value: [feature] }, biologicalFeatures: { value: [feature] }
  } });
  const item = (scope) => ({ target: {
    kind: 'featureIdentity', scope, recordKey: 'record-1', biologicalFeatureId: 'f1'
  } });
  assert.doesNotMatch(actions.featureTargetCaption(item('circular')), /not in the current diagram/);
  assert.equal(actions.featureTargetCaption(item('linear')), 'f1 (not in the current diagram)');
});
const created = editor.addAnnotationSet('review');
const addedCoordinate = editor.addCoordinateAnnotation(created, { start: 5, end: 8 });
assert.equal(addedCoordinate.mark, 'highlight');
assert.equal(created.defaultStyle.fill, '#94a3b8');
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
assert.deepEqual(validateAnnotationRecordTargets([targetSet], linearCatalog), { code: 'ANNOTATION_TARGET', context: { reason: 'TARGET_RECORD' } });
const featureTargetSet = createAnnotationSet({
  id: 'features',
  annotations: [{ id: 'gene', target: featureTarget({ selector: 'locus_tag=ABC_1' }), label: '', mark: 'bracket' }]
});
assert.deepEqual(validateAnnotationRecordTargets([featureTargetSet], linearCatalog), { code: 'ANNOTATION_TARGET', context: { reason: 'TARGET_RECORD' } });
assert.deepEqual(
  annotationRecordOptions(linearCatalog, targetSet.annotations[0]).map((option) => option.value),
  ['', linearCatalog.records[0].key, linearCatalog.records[1].key]
);
setAnnotationRecordValue(linearCatalog, targetSet.annotations[0], linearCatalog.records[1].key);
assert.deepEqual(targetSet.annotations[0].target.record, { kind: 'recordId', value: 'RecB' });
assert.equal(validateAnnotationRecordTargets([targetSet], linearCatalog), null);
assert.equal(encodeAnnotationTable([targetSet]).trim().split('\n')[1].split('\t')[3], 'RecB');
delete targetSet.annotations[0].metadata._gbdraw_web_target_record_key;
targetSet.annotations[0].target.record = { kind: 'recordIndex', index: 2 };
assert.deepEqual(validateAnnotationRecordTargets([targetSet], linearCatalog), { code: 'ANNOTATION_TARGET', context: { reason: 'OUT_OF_RANGE' } });
targetSet.annotations[0].target.record = { kind: 'recordId', value: 'missing' };
assert.deepEqual(validateAnnotationRecordTargets([targetSet], linearCatalog), { code: 'ANNOTATION_TARGET', context: { reason: 'NO_MATCH' } });
setAnnotationRecordValue(duplicateCatalog, targetSet.annotations[0], duplicateCatalog.records[0].key);
assert.equal(validateAnnotationRecordTargets([targetSet], duplicateCatalog), null);
targetSet.annotations[0].target.record = { kind: 'recordIndex', index: -1 };
delete targetSet.annotations[0].metadata._gbdraw_web_target_record_key;
assert.deepEqual(validateAnnotationRecordTargets([targetSet], linearCatalog), { code: 'ANNOTATION_TARGET', context: { reason: 'TARGET_RECORD' } });

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
assert.deepEqual(automaticSourceCatalog.records.map((record) => record.sourceIndex), [0, 0]);
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
const emptySourceCatalog = buildAnnotationRecordCatalog({
  mode: 'linear',
  linearSources: [{ sourceKey: 'empty', hasInput: true, status: 'ready', selector: '', records: [] }]
});
assert.equal(emptySourceCatalog.status, 'error');
assert.deepEqual(emptySourceCatalog.issues[0], { code: 'NO_RECORDS', context: { inputOrdinal: 1 } });

const gbComparisonCatalog = buildAnnotationRecordCatalog({
  mode: 'linear',
  linearSources: [
    { sourceKey: 'gb-a', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'A1' }, { recordId: 'A2' }] },
    { sourceKey: 'gb-b', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'B1' }, { recordId: 'B2' }] }
  ]
});
assert.deepEqual(gbComparisonCatalog.records.map((record) => record.recordId), ['A1', 'A2', 'B1', 'B2']);
assert.deepEqual(gbComparisonCatalog.records.map((record) => record.sourceIndex), [0, 0, 1, 1]);
const singleGbComparisonCatalog = buildAnnotationRecordCatalog({
  mode: 'linear',
  linearSources: [{ sourceKey: 'gb-only', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'A1' }, { recordId: 'A2' }] }]
});
assert.deepEqual(singleGbComparisonCatalog.records.map((record) => record.recordId), ['A1', 'A2']);
const gffComparisonCatalog = buildAnnotationRecordCatalog({
  mode: 'linear',
  linearSources: [{ sourceKey: 'gff-only', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'G1' }, { recordId: 'G2' }] }]
});
assert.deepEqual(gffComparisonCatalog.records.map((record) => record.recordId), ['G1', 'G2']);
const comparisonSequences = [{ uid: 'a' }, { uid: 'b' }];
assert.equal(hasLinearComparisonIntent(resolveLinearComparisonPlan({
  plan: { mode: 'adjacent', defaultSource: 'losat', edges: [] },
  sequences: comparisonSequences
})), true);
assert.equal(hasLinearComparisonIntent(resolveLinearComparisonPlan({
  plan: { mode: 'none', defaultSource: 'losat', edges: [] },
  sequences: comparisonSequences
})), false);
assert.equal(hasLinearComparisonIntent(resolveLinearComparisonPlan({
  plan: { mode: 'selected', defaultSource: 'losat', edges: [] },
  sequences: comparisonSequences
})), false);
assert.equal(hasLinearComparisonIntent(resolveLinearComparisonPlan({
  plan: { mode: 'adjacent', defaultSource: 'upload', edges: [] },
  sequences: comparisonSequences
})), false);
assert.equal(hasLinearComparisonIntent(resolveLinearComparisonPlan({
  plan: {
    mode: 'selected',
    defaultSource: 'upload',
    edges: [createLinearComparisonEdge({
      queryUid: 'a', subjectUid: 'b', source: 'upload',
      file: { name: 'a-b.tsv' }, fileActive: true
    })]
  },
  sequences: comparisonSequences
})), true);

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
assert.equal(validateAnnotationRecordTargets([duplicateTargetSet], reorderedDuplicateCatalog), null);
const replacedSourceCatalog = buildAnnotationRecordCatalog({
  mode: 'linear',
  linearSources: [{ sourceKey: 'replacement', hasInput: true, status: 'ready', selector: '', records: [{ recordId: 'dup' }] }]
});
assert.deepEqual(validateAnnotationRecordTargets([duplicateTargetSet], replacedSourceCatalog), { code: 'ANNOTATION_TARGET', context: { reason: 'TARGET_RECORD' } });

// FE-05 (PD-OI-041): a Circular request that draws several records, in one
// figure or one output per record, needs an explicit target record, which
// Python binds to that record alone. A target without a record is ambiguous.
targetSet.annotations[0].target.record = { kind: 'recordId', value: 'RecA' };
targetSet.annotations[0].metadata = {};
const circularMultiOutputCatalog = buildAnnotationRecordCatalog({
  mode: 'circular',
  circularSource: {
    sourceKey: 'circular', hasInput: true, status: 'ready',
    records: [{ record_id: 'RecA' }, { record_id: 'RecB' }]
  }
});
assert.equal(circularMultiOutputCatalog.requiresSelection, true);
assert.equal(validateAnnotationRecordTargets([targetSet], circularMultiOutputCatalog), null);
assert.deepEqual(
  annotationRecordOptions(circularMultiOutputCatalog, targetSet.annotations[0]).map((option) => option.label),
  ['Select target record', '#1 · RecA', '#2 · RecB']
);
targetSet.annotations[0].target.record = null;
assert.deepEqual(validateAnnotationRecordTargets([targetSet], circularMultiOutputCatalog), { code: 'ANNOTATION_TARGET', context: { reason: 'TARGET_RECORD' } });
targetSet.annotations[0].target.record = { kind: 'recordId', value: 'RecA' };
// The catalog offers the records the Circular request draws: the selected
// record of a single presentation, every record of a grid or batch.
const circularRecords = [{ selector: '#1', record_id: 'RecA' }, { selector: '#2', record_id: 'RecB' }];
const drawnIds = (options) => resolveCircularRequestRecordSet({ records: circularRecords, ...options })
  .records.map((record) => record.recordId);
assert.deepEqual(drawnIds({ selector: 'RecB' }), ['RecB']);
assert.deepEqual(drawnIds({ selector: '' }), ['RecA', 'RecB']);
assert.deepEqual(drawnIds({ selector: 'RecB', multiRecordCanvas: true }), ['RecA', 'RecB']);
assert.deepEqual(drawnIds({ selector: 'RecB', groupingIntent: 'batch' }), ['RecA', 'RecB']);
assert.equal(resolveCircularRequestRecordSet({ records: circularRecords, selector: 'missing' }).selectionFailure, 'NO_MATCH');
const circularSelectedCatalog = buildAnnotationRecordCatalog({
  mode: 'circular',
  circularSource: {
    sourceKey: 'circular', hasInput: true, status: 'ready',
    records: resolveCircularRequestRecordSet({ records: circularRecords, selector: 'RecB' }).records
  }
});
assert.deepEqual(circularSelectedCatalog.records.map((record) => record.recordId), ['RecB']);
assert.equal(circularSelectedCatalog.requiresSelection, false);
targetSet.annotations[0].target.record = { kind: 'recordId', value: 'RecB' };
assert.equal(validateAnnotationRecordTargets([targetSet], circularSelectedCatalog), null);
targetSet.annotations[0].target.record = { kind: 'recordId', value: 'RecA' };

const circularSingleCatalog = buildAnnotationRecordCatalog({
  mode: 'circular',
  circularSource: { sourceKey: 'single', hasInput: true, status: 'ready', records: [{ record_id: 'RecA' }] }
});
targetSet.annotations[0].target.record = null;
assert.equal(validateAnnotationRecordTargets([targetSet], circularSingleCatalog), null);

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
assert.deepEqual(validateAnnotationRecordTargets(nullRecordSets, linearCatalog), { code: 'ANNOTATION_TARGET', context: { reason: 'TARGET_RECORD' } });
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

const importCases = JSON.parse(await readFile(join(process.cwd(), 'tests/fixtures/annotations/tsv-import-cases.json'), 'utf8'));
for (const entry of importCases) {
  test(`Annotation TSV parity: ${entry.name}`, () => {
    if (!entry.valid) {
      assert.throws(() => parseAnnotationTable(entry.table));
      assert.throws(() => parseAnnotationTableWithNotice(entry.table));
      return;
    }
    const parsed = parseAnnotationTableWithNotice(entry.table);
    assert.deepEqual(parsed.sets, parseAnnotationTable(entry.control));
    assert.deepEqual(parseAnnotationTable(entry.table), parsed.sets);
    assert.equal(Array.isArray(parsed.sets), true);
    if (entry.ignored.length) {
      for (const name of entry.ignored) assert.ok(parsed.notice.includes(name));
      assert.match(parsed.notice, /not saved in Sessions or TSV re-export/);
    } else assert.equal(parsed.notice, '');
    assert.ok(!JSON.stringify(parsed).includes('PRIVATE-CELL'));
    const encoded = encodeAnnotationTable(parsed.sets);
    for (const name of entry.ignored) assert.ok(!encoded.split('\n')[0].split('\t').includes(name));
  });
}

test('quotes in a label are part of the value through Web write and Web read', () => {
  // The same vector tests/test_annotations.py reads through read_annotation_table (OV-24).
  const vector = importCases.find((entry) => entry.name === 'quotes are part of a label');
  const expected = ['"lead', '"quoted"', 'a"b', 'plain'];
  const sets = parseAnnotationTable(vector.table);
  assert.deepEqual(sets[0].annotations.map((annotation) => annotation.label), expected);
  const encoded = encodeAnnotationTable(sets);
  const [header, ...rows] = encoded.trimEnd().split('\n').map((line) => line.split('\t'));
  assert.deepEqual(rows.map((row) => row[header.indexOf('label')]), expected);
  assert.deepEqual(parseAnnotationTable(encoded), sets);
});

test('file import commits once, separates notices, and preserves draft/Result on failure or stale completion', async () => {
  const state = { annotationSets: [createAnnotationSet({ id: 'before', annotations: [
    { id: 'old', target: coordinateTarget({ start: 1, end: 3 }), mark: 'band' }
  ] })], results: [{ content: '<svg/>', warnings: ['resolved warning'] }] };
  const notices = [];
  const editor = createAnnotationEditor({ state, onImportNotice: (notice) => notices.push(notice) });
  const alerts = [];
  const oldAlert = globalThis.alert;
  globalThis.alert = (message) => alerts.push(message);
  let commits = 0;
  const splice = state.annotationSets.splice.bind(state.annotationSets);
  Object.defineProperty(state.annotationSets, 'splice', { value: (...args) => { commits++; return splice(...args); } });
  const good = importCases[0];
  const input = { files: [{ text: async () => good.table }], value: 'annotations.tsv' };
  try {
    const result = structuredClone(state.results);
    await editor.importAnnotationTableFile({ target: input });
    assert.equal(commits, 1);
    assert.deepEqual(state.annotationSets, parseAnnotationTable(good.control));
    assert.equal(notices.filter(Boolean).length, 1);
    assert.deepEqual(state.results, result);
    for (const entry of importCases.filter((entry) => !entry.valid)) {
      const before = JSON.stringify(state);
      input.files = [{ text: async () => entry.table }];
      await editor.importAnnotationTableFile({ target: input });
      assert.equal(JSON.stringify(state), before);
      assert.equal(commits, 1);
      assert.equal(notices.at(-1), '');
    }
    const before = JSON.stringify(state);
    input.files = [{ text: async () => { throw new Error('read failure'); } }];
    await editor.importAnnotationTableFile({ target: input });
    assert.equal(JSON.stringify(state), before);
    assert.match(alerts.at(-1), /read failure/);
    let finish;
    input.files = [{ text: () => new Promise((resolve) => { finish = resolve; }) }];
    const pending = editor.importAnnotationTableFile({ target: input });
    state.annotationSets[0].annotations[0].label = 'newer draft edit';
    const edited = JSON.stringify(state);
    finish(good.table);
    await pending;
    assert.equal(JSON.stringify(state), edited);
    assert.equal(commits, 1);
    assert.equal(notices.at(-1), '');
    let finishResultRead;
    input.files = [{ text: () => new Promise((resolve) => { finishResultRead = resolve; }) }];
    const staleResult = editor.importAnnotationTableFile({ target: input });
    state.results.splice(0, 1, { content: '<svg>new result</svg>', warnings: [] });
    const replacedResult = JSON.stringify(state);
    finishResultRead(good.table);
    await staleResult;
    assert.equal(JSON.stringify(state), replacedResult);
    assert.equal(commits, 1);
    state.results.splice(0, 1, ...result);
    let finishOld;
    input.files = [{ text: () => new Promise((resolve) => { finishOld = resolve; }) }];
    const obsolete = editor.importAnnotationTableFile({ target: input });
    input.files = [{ text: async () => good.control }];
    await editor.importAnnotationTableFile({ target: input });
    const current = JSON.stringify(state);
    finishOld(good.table);
    await obsolete;
    assert.equal(JSON.stringify(state), current);
    assert.equal(commits, 2);
    assert.equal(notices.at(-1), '');
    assert.deepEqual(state.results, result);
  } finally {
    if (oldAlert === undefined) delete globalThis.alert;
    else globalThis.alert = oldAlert;
  }
});

test('record reconciliation failure cannot partially replace an imported draft', () => {
  const state = { annotationSets: [createAnnotationSet({ id: 'original' })] };
  const before = JSON.stringify(state);
  const editor = createAnnotationEditor({ state, getRecordCatalog: () => { throw new Error('catalog unavailable'); } });
  assert.throws(() => editor.importAnnotationTable(importCases[0].table), /catalog unavailable/);
  assert.equal(JSON.stringify(state), before);
});
