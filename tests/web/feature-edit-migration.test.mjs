// Session 44 -> 45 key migrator (design Q4 4.3): old per-feature edits keyed by
// rendered SVG ID become identity rows through the Session's saved catalog.
// The positive fixtures are Sessions written by first-parent main
// (tests/fixtures/sessions/feature-edits.provenance.json).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import { gunzipSync } from 'node:zlib';
import {
  migrateRenderedIdFeatureEdits,
  migrateSessionAnnotationTargets,
  migrateSessionFeatureEdits,
  migrateSessionFeaturePlacements
} from '../../gbdraw/web/js/services/feature-edit-migration.js';
import { canonicalFeatureOverrides } from '../../gbdraw/web/js/services/feature-placement.js';
import { annotationOptionsPayload, normalizeAnnotationSets } from '../../gbdraw/web/js/services/annotation-state.js';

const fixture = (name) => JSON.parse(gunzipSync(readFileSync(new URL(`../fixtures/sessions/${name}`, import.meta.url))));
// Migrated rows name the mode of the Session's diagram (R2).
const identityRows = (scope) => ({
  key: (recordKey, featureId) => JSON.stringify([scope, recordKey, featureId]),
  row: (recordKey, biologicalFeatureId, fields) => ({
    scope,
    recordKey,
    biologicalFeatureId,
    featureVisibility: null,
    labelVisibility: null,
    labelText: null,
    labelSourceText: null,
    ...fields
  })
});

test('v44 Linear crop and reverse complement: rule 1 maps each drawn ID to its feature', () => {
  const session = fixture('feature-edits-crop-rc.v44.gbdraw-session.json.gz');
  assert.equal(session.version, 44);
  const [testa, testb] = session.renderRequest.records.map((record) => record.recordKey);
  const { key, row } = identityRows('linear');
  const { features, droppedCount, narrowedVisibilityCount } = migrateSessionFeatureEdits({
    features: session.features, mode: session.renderRequest.mode, catalog: session.editorState.featureCatalog
  });
  assert.equal(droppedCount, 0);
  assert.equal(narrowedVisibilityCount, 0);
  assert.deepEqual(features.featureOverrides, {
    // OV-01: the cropped record's hidden misc_feature X, by its source identity.
    [key(testa, 'fb5977f81')]: row(testa, 'fb5977f81', { featureVisibility: 'off' }),
    [key(testa, 'f571eb983')]: row(testa, 'f571eb983', {
      labelText: 'PROBE_A_TRNA', labelSourceText: 'tRNA-Leu'
    }),
    [key(testb, 'f88047061')]: row(testb, 'f88047061', { featureVisibility: 'off', labelVisibility: 'off' })
  });
  for (const field of [
    'featureVisibilityOverrides', 'labelVisibilityOverrides', 'labelTextFeatureOverrides', 'labelTextFeatureOverrideSources'
  ]) assert.equal(Object.hasOwn(features, field), false, field);
  // The saved label table was built from the migrated maps.
  assert.deepEqual(features.labelOverrideRows, []);
  assert.equal(canonicalFeatureOverrides(features.featureOverrides).length, 3);
});

test('v44 Circular canvas with one record twice: each copy keeps its own edits', () => {
  const session = fixture('feature-edits-circular-copies.v44.gbdraw-session.json.gz');
  const { key, row } = identityRows('circular');
  const { features, droppedCount, narrowedVisibilityCount } = migrateSessionFeatureEdits({
    features: session.features, mode: session.renderRequest.mode, catalog: session.editorState.featureCatalog
  });
  assert.equal(droppedCount, 0);
  // Session 44 sent the Feature visibility edit as a `hash` row, which hid the
  // feature in both copies; it now names copy 1 only (Owner decision Q1 = A),
  // so Load says that Generate draws copy 2's feature again.
  assert.equal(narrowedVisibilityCount, 1);
  assert.deepEqual(features.featureOverrides, {
    // Rule 2: the hidden feature is not drawn; `__instance_record_1` names copy 1.
    [key('record-1', 'f3ccacda4')]: row('record-1', 'f3ccacda4', { featureVisibility: 'off' }),
    [key('record-1', 'f01237d96')]: row('record-1', 'f01237d96', {
      labelText: 'COPY1_ONLY', labelSourceText: 'gtg start'
    }),
    [key('record-1', 'fef810304')]: row('record-1', 'fef810304', { labelVisibility: 'off' })
  });
});

test('v33 Session without a catalog maps its edits through its feature metadata', () => {
  const session = fixture('feature-edits-circular.v33.gbdraw-session.json.gz');
  assert.equal(session.version, 33);
  const [recordKey] = session.renderRequest.records.map((record) => record.recordKey);
  const { key, row } = identityRows('circular');
  const { featureOverrides, droppedCount } = migrateRenderedIdFeatureEdits({
    features: session.features,
    mode: 'circular',
    legacy: { features: session.features.extractedFeatures, records: session.renderRequest.records }
  });
  assert.equal(droppedCount, 0);
  assert.deepEqual(featureOverrides, {
    [key(recordKey, 'f406d90f1')]: row(recordKey, 'f406d90f1', { featureVisibility: 'off' }),
    [key(recordKey, 'fbe3a7c0c')]: row(recordKey, 'fbe3a7c0c', {
      labelText: 'V33_LABEL', labelSourceText: 'NADH dehydrogenase subunit 2'
    }),
    [key(recordKey, 'f48b7bf2f')]: row(recordKey, 'f48b7bf2f', { labelVisibility: 'off' })
  });
});

test('rule 3: an edit matching no feature is dropped and counted', () => {
  const session = fixture('feature-edits-crop-rc.v44.gbdraw-session.json.gz');
  const { droppedCount, featureOverrides } = migrateRenderedIdFeatureEdits({
    features: {
      featureVisibilityOverrides: { f0000000a_record_1: 'off', f3b928d8c_record_1: 'off' },
      labelTextFeatureOverrides: { f0000000b: 'gone' }
    },
    mode: 'linear',
    catalog: session.editorState.featureCatalog
  });
  assert.equal(droppedCount, 2);
  assert.equal(Object.keys(featureOverrides).length, 1);
  // Without a catalog or metadata no key resolves.
  assert.equal(migrateRenderedIdFeatureEdits({
    features: { featureVisibilityOverrides: { f3b928d8c_record_1: 'off' } }
  }).droppedCount, 1);
});

test('rule 2 reads record and instance suffixes and needs exactly one candidate', () => {
  const catalog = {
    schema: 4,
    items: [{
      recordKeys: ['seq-a', 'seq-b'],
      features: [],
      biologicalFeatures: [
        { recordKey: 'seq-a', biologicalFeatureId: 'f1~3', stableFeatureId: 'f1', sourceFeatureIndex: 3 },
        { recordKey: 'seq-a', biologicalFeatureId: 'f1~4', stableFeatureId: 'f1', sourceFeatureIndex: 4 },
        { recordKey: 'seq-b', biologicalFeatureId: 'f2', sourceFeatureIndex: 1 }
      ]
    }]
  };
  const { key } = identityRows('linear');
  const { featureOverrides, droppedCount } = migrateRenderedIdFeatureEdits({
    features: {
      featureVisibilityOverrides: {
        'f1_record_1__instance_4_0123456789abcdef': 'off',
        f1_record_1: 'on',
        f2_record_2: 'exclude_matching'
      }
    },
    mode: 'linear',
    catalog
  });
  assert.equal(droppedCount, 1);
  assert.deepEqual(Object.keys(featureOverrides).sort(), [key('seq-a', 'f1~4'), key('seq-b', 'f2')].sort());
  assert.equal(featureOverrides[key('seq-b', 'f2')].featureVisibility, 'exclude_matching');
});

test('a blank label text hid its label and becomes Label visibility Off', () => {
  const { key, row } = identityRows('circular');
  const { featureOverrides } = migrateRenderedIdFeatureEdits({
    features: { labelTextFeatureOverrides: { f2: '  ' } },
    mode: 'circular',
    catalog: { schema: 4, items: [{ recordKeys: ['r'], features: [{ svgId: 'f2', recordKey: 'r', biologicalFeatureId: 'f2' }], biologicalFeatures: [] }] }
  });
  assert.deepEqual(featureOverrides[key('r', 'f2')], row('r', 'f2', { labelVisibility: 'off' }));
});

// Sessions 31-33 Linear: each input is one request record (an ALL record
// without a selector), every input's features have record_idx 0, and the
// rendered ID `<drawn hash>_record_<n>` gives the record's position.
const linearRecords = () => [
  { recordKey: 'seq-a', cardinality: 'exactly_one', region: { selector: null, start: 201, end: 3800, reverseComplement: false },
    presentation: { reverseComplement: false } },
  { recordKey: 'seq-b', cardinality: 'all', region: null, presentation: { reverseComplement: true } },
  { recordKey: 'seq-c', cardinality: 'all', region: null, presentation: { reverseComplement: false } }
];

test('v31-33 Linear: source features read again name each record by its input and drawn hash', () => {
  // As the Session's sources read with its crop and orientation give them:
  // the source hash beside the drawn one, per input (`fileIdx`). Input 3 holds
  // two records, so its record keys are `seq-c:1` and `seq-c:2`.
  const source = (fileIdx, recordIdx, featureIndex, stable, drawn) => ({
    fileIdx, record_idx: recordIdx, feature_index: featureIndex, stable_feature_id: stable,
    svg_id: stable, drawn_selector: { hash: drawn }
  });
  const features = [
    source(0, 0, 2, 'fsrca2', 'fdrwa2'),
    source(1, 0, 2, 'fsrcb2', 'fdrwb2'),
    source(1, 0, 4, 'fsame', 'fdrwb4'),
    source(2, 0, 1, 'fsame', 'fsame'),
    source(2, 1, 1, 'fsame', 'fsame'),
    source(2, 1, 3, 'fdup', 'fdup'),
    source(2, 1, 5, 'fdup', 'fdup')
  ];
  const { key, row } = identityRows('linear');
  const { featureOverrides, droppedCount } = migrateRenderedIdFeatureEdits({
    features: {
      featureVisibilityOverrides: {
        fdrwa2_record_1: 'off',
        fdrwb2_record_2: 'off',
        fsame_record_4: 'exclude_matching',
        'fdup_record_4__instance_5_0123456789abcdef': 'off'
      },
      labelTextFeatureOverrides: { fdrwb4_record_2: 'B_LABEL' },
      labelVisibilityOverrides: { fsame_record_3: 'off' }
    },
    mode: 'linear',
    legacy: { features, records: linearRecords() }
  });
  assert.equal(droppedCount, 0);
  assert.deepEqual(featureOverrides, {
    [key('seq-a', 'fsrca2')]: row('seq-a', 'fsrca2', { featureVisibility: 'off' }),
    [key('seq-b', 'fsrcb2')]: row('seq-b', 'fsrcb2', { featureVisibility: 'off' }),
    [key('seq-b', 'fsame')]: row('seq-b', 'fsame', { labelText: 'B_LABEL' }),
    [key('seq-c:1', 'fsame')]: row('seq-c:1', 'fsame', { labelVisibility: 'off' }),
    [key('seq-c:2', 'fsame')]: row('seq-c:2', 'fsame', { featureVisibility: 'exclude_matching' }),
    [key('seq-c:2', 'fdup~5')]: row('seq-c:2', 'fdup~5', { featureVisibility: 'off' })
  });
});

test('v31-33 Linear: saved metadata names features only of records drawn untransformed', () => {
  // Saved metadata holds the drawn hash only; it equals the source hash where
  // the record was neither cropped nor reverse-complemented.
  const saved = (fileIdx, svgId) => ({
    fileIdx, record_idx: 0, svg_id: svgId, stable_svg_id: svgId.replace(/_record_\d+$/, ''),
    stable_feature_id: svgId.replace(/_record_\d+$/, '')
  });
  const { key, row } = identityRows('linear');
  const { featureOverrides, droppedCount } = migrateRenderedIdFeatureEdits({
    features: {
      featureVisibilityOverrides: { fdrwa2_record_1: 'off', fdrwb2_record_2: 'off', fplain_record_3: 'off' }
    },
    mode: 'linear',
    legacy: {
      features: [saved(0, 'fdrwa2_record_1'), saved(1, 'fdrwb2_record_2'), saved(2, 'fplain_record_3')],
      records: linearRecords()
    }
  });
  assert.equal(droppedCount, 2);
  assert.deepEqual(featureOverrides, {
    [key('seq-c', 'fplain')]: row('seq-c', 'fplain', { featureVisibility: 'off' })
  });
});

// R-7, Owner decision 2026-10-05: a Session before 45 named a selected feature
// in an annotation by `hash=<hash>` (a featureSpan target), which the renderer
// matches in the drawn record. Load moves such a target to the feature's
// source identity only when the figure cannot change: its record is drawn
// without a crop, reverse complement, or rotation, and the hash names exactly
// one feature of the saved catalog, in the record the target binds. Every
// other target stays as saved.
const hashTarget = (hash, record = null, extra = {}) => ({
  kind: 'featureSpan', record, selectors: [{ key: 'hash', value: hash }],
  envelope: 'segments', circularPath: 'reverse', ...extra
});
const annotationSets = (targets) => [{
  id: 'marks', defaultStyle: {}, legendLabel: null,
  annotations: targets.map((target, index) => ({
    id: `a${index + 1}`, target, label: `L${index + 1}`, mark: 'band', lane: null, style: null, legendLabel: null,
    metadata: { kept: String(index) }
  }))
}];
const byId = (recordId) => ({ kind: 'recordId', value: recordId });
const byIndex = (index) => ({ kind: 'recordIndex', index });
const identity = (recordKey, biologicalFeatureId) => ({
  kind: 'featureIdentity', scope: 'linear', recordKey, biologicalFeatureId, envelope: 'segments', circularPath: 'reverse'
});
// Three records drawn in one Result: `plain` untransformed, `copy` the same
// record (ID A) again, and `other` (ID B). `plain` has two CDS at the same
// coordinates (f2~1, f2~2); f3 is in both copies of A.
const synthetic = ({ records: recordOverrides = {}, recordKeys = ['plain', 'copy', 'other'] } = {}) => ({
  catalog: {
    schema: 4,
    items: [{
      recordKeys,
      features: [],
      biologicalFeatures: [
        { recordKey: 'plain', record_id: 'A', biologicalFeatureId: 'f1', sourceFeatureIndex: 1 },
        { recordKey: 'plain', record_id: 'A', biologicalFeatureId: 'f2~2', sourceFeatureIndex: 2 },
        { recordKey: 'plain', record_id: 'A', biologicalFeatureId: 'f2~3', sourceFeatureIndex: 3 },
        { recordKey: 'plain', record_id: 'A', biologicalFeatureId: 'f3', sourceFeatureIndex: 4 },
        { recordKey: 'copy', record_id: 'A', biologicalFeatureId: 'f3', sourceFeatureIndex: 4 },
        { recordKey: 'other', record_id: 'B', biologicalFeatureId: 'f4', sourceFeatureIndex: 1 },
        { recordKey: 'other', record_id: 'B', biologicalFeatureId: 'f5', sourceFeatureIndex: 2 }
      ].filter((feature) => recordKeys.includes(feature.recordKey))
    }]
  },
  records: recordKeys.map((recordKey) => ({
    recordKey, cardinality: 'exactly_one', region: null,
    presentation: { reverseComplement: false }, display: { isCircular: null, startCoordinate: null },
    ...recordOverrides[recordKey]
  }))
});
const migrate = (targets, context = synthetic(), mode = 'linear') => {
  const sets = annotationSets(targets);
  const before = JSON.stringify(sets);
  const result = migrateSessionAnnotationTargets({ annotationSets: sets, mode, ...context });
  assert.equal(JSON.stringify(sets), before, 'the saved sets are not changed in place');
  return { ...result, targets: result.annotationSets[0].annotations.map((item) => item.target) };
};

test('R-7: a hash target of one feature of an untransformed record moves to its identity', () => {
  const { targets, migratedCount, annotationSets: migrated } = migrate([
    hashTarget('f1', byId('B')),
    hashTarget('f4', byId('B')),
    hashTarget('f5', byIndex(2)),
    hashTarget('f5')
  ], synthetic({ recordKeys: ['other'] }));
  // a1 names f1, which record B does not have: the renderer skipped it.
  assert.deepEqual(targets, [hashTarget('f1', byId('B')), identity('other', 'f4'),
    hashTarget('f5', byIndex(2)), identity('other', 'f5')]);
  assert.equal(migratedCount, 2);
  // Only the target changes; the annotation keeps its ID, label, style, and metadata.
  assert.deepEqual({ ...migrated[0].annotations[1], target: null }, {
    id: 'a2', target: null, label: 'L2', mark: 'band', lane: null, style: null, legendLabel: null, metadata: { kept: '1' }
  });
  // The draft owner accepts the target, and the Linear request carries it.
  const records = [{ recordKey: 'other', cardinality: 'exactly_one' }];
  const payload = annotationOptionsPayload(normalizeAnnotationSets(migrated), 'linear', records);
  assert.deepEqual(payload.sets[0].annotations[1].target, {
    kind: 'featureIdentity', recordKey: 'other', biologicalFeatureId: 'f4', envelope: 'segments', circularPath: 'reverse'
  });
  assert.deepEqual(annotationOptionsPayload(normalizeAnnotationSets(migrated), 'circular', records)
    .sets[0].annotations.map((item) => item.id), ['a1', 'a3']);
});

test('R-7: the record the target binds is read as the renderer reads it', () => {
  const kept = [
    // The renderer binds #1 (plain), which lacks f4, or no record (#4).
    hashTarget('f4', byIndex(0)),
    hashTarget('f4', byIndex(3)),
    // Records plain and copy have ID A.
    hashTarget('f1', byId('A')),
    // Without a record selector several records are ambiguous.
    hashTarget('f4')
  ];
  const { targets, migratedCount } = migrate([hashTarget('f4', byIndex(2)), hashTarget('f4', byId('B')), ...kept]);
  assert.deepEqual(targets, [identity('other', 'f4'), identity('other', 'f4'), ...kept]);
  assert.equal(migratedCount, 2);
});

test('R-7: a cropped, reverse-complemented, or rotated record keeps its hash targets', () => {
  const targets = [hashTarget('f4', byId('B'))];
  for (const transform of [
    { region: { selector: null, start: 1, end: 4000, reverseComplement: false } },
    { presentation: { reverseComplement: true } },
    { display: { isCircular: null, startCoordinate: 101 } }
  ]) {
    const result = migrate(targets, synthetic({ records: { other: transform } }));
    assert.deepEqual(result.targets, targets, JSON.stringify(transform));
    assert.equal(result.migratedCount, 0);
  }
  // A transform of another record does not matter.
  assert.equal(migrate(targets, synthetic({ records: { plain: { presentation: { reverseComplement: true } } } }))
    .migratedCount, 1);
});

test('R-7: a hash that names no feature or several features keeps its target', () => {
  const kept = [
    hashTarget('f9', byIndex(0)),
    // Two CDS at the same coordinates share the hash (f2~2, f2~3).
    hashTarget('f2', byIndex(0)),
    hashTarget('f2~2', byIndex(0)),
    // f3 is in both copies of record A.
    hashTarget('f3', byIndex(0))
  ];
  const { targets, migratedCount } = migrate(kept);
  assert.deepEqual(targets, kept);
  assert.equal(migratedCount, 0);
});

test('R-7: only a featureSpan with one hash selector moves', () => {
  const kept = [
    { ...hashTarget('f4', byId('B')), selectors: [{ key: 'locus_tag', value: 'B_1' }] },
    { ...hashTarget('f4', byId('B')), selectors: [{ key: 'hash', value: 'f4' }, { key: 'hash', value: 'f5' }] },
    { ...hashTarget('f4', byId('B')), selectors: [{ key: null, value: 'f4' }] },
    { kind: 'coordinateSpan', record: byId('B'), start: 1, end: 5, coordinateSpace: 'source', wrapsOrigin: false, outOfBounds: 'clip' },
    identity('other', 'f5')
  ];
  const { targets, migratedCount } = migrate(kept);
  assert.deepEqual(targets, kept);
  assert.equal(migratedCount, 0);
  // Without a saved catalog nothing names a source feature.
  assert.equal(migrate([hashTarget('f4', byId('B'))], { catalog: null, records: [] }).migratedCount, 0);
});

test('R-7: a record of an ALL input is read through its request record', () => {
  const context = synthetic({ recordKeys: ['all-input:1', 'all-input:2'] });
  context.catalog.items[0].biologicalFeatures = [
    { recordKey: 'all-input:1', record_id: 'A', biologicalFeatureId: 'f1', sourceFeatureIndex: 1 },
    { recordKey: 'all-input:2', record_id: 'B', biologicalFeatureId: 'f4', sourceFeatureIndex: 1 }
  ];
  context.records = [{ recordKey: 'all-input', cardinality: 'all', region: null,
    presentation: { reverseComplement: false }, display: { isCircular: null, startCoordinate: null } }];
  assert.deepEqual(migrate([hashTarget('f4', byId('B'))], context).targets, [identity('all-input:2', 'f4')]);
  context.records[0].presentation.reverseComplement = true;
  assert.equal(migrate([hashTarget('f4', byId('B'))], context).migratedCount, 0);
});

// The Session 44 positive fixture (selected-feature-annotations.provenance.json):
// feature_1 names one feature of an untransformed record, feature_2 one of two
// CDS at the same coordinates, feature_3 a feature of a reverse-complemented
// record.
test('R-7: the Session 44 fixture moves only the annotation whose feature is certain', () => {
  const session = fixture('selected-feature-annotations.v44.gbdraw-session.json.gz');
  assert.equal(session.version, 44);
  const [testa] = session.renderRequest.records.map((record) => record.recordKey);
  const { annotationSets: migrated, migratedCount } = migrateSessionAnnotationTargets({
    annotationSets: session.config.annotationSets,
    mode: session.renderRequest.mode,
    catalog: session.editorState.featureCatalog,
    records: session.renderRequest.records
  });
  assert.equal(migratedCount, 1);
  const saved = session.config.annotationSets[0].annotations;
  assert.deepEqual(migrated[0].annotations.map((item) => item.target), [
    { kind: 'featureIdentity', scope: 'linear', recordKey: testa, biologicalFeatureId: 'fef810304',
      envelope: 'outer_bounds', circularPath: 'shortest' },
    saved[1].target,
    saved[2].target
  ]);
});

// The same vector pins the Python reader (tests/test_session_compat.py).
test('Session 41-44 placement drafts migrate to the vector both readers share', () => {
  const vector = JSON.parse(readFileSync(new URL('../fixtures/feature-placement-migration.json', import.meta.url), 'utf8'));
  assert.deepEqual(migrateSessionFeaturePlacements(structuredClone(vector.input)), vector.expected);
});
