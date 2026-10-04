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
  migrateSessionFeatureEdits
} from '../../gbdraw/web/js/services/feature-edit-migration.js';
import { canonicalFeatureOverrides } from '../../gbdraw/web/js/services/feature-placement.js';

const fixture = (name) => JSON.parse(gunzipSync(readFileSync(new URL(`../fixtures/sessions/${name}`, import.meta.url))));
const key = (recordKey, featureId) => JSON.stringify([recordKey, featureId]);
const row = (recordKey, biologicalFeatureId, fields) => ({
  recordKey,
  biologicalFeatureId,
  featureVisibility: null,
  labelVisibility: null,
  labelText: null,
  labelSourceText: null,
  ...fields
});

test('v44 Linear crop and reverse complement: rule 1 maps each drawn ID to its feature', () => {
  const session = fixture('feature-edits-crop-rc.v44.gbdraw-session.json.gz');
  assert.equal(session.version, 44);
  const [testa, testb] = session.renderRequest.records.map((record) => record.recordKey);
  const { features, droppedCount } = migrateSessionFeatureEdits({
    features: session.features, catalog: session.editorState.featureCatalog
  });
  assert.equal(droppedCount, 0);
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
  const { features, droppedCount } = migrateSessionFeatureEdits({
    features: session.features, catalog: session.editorState.featureCatalog
  });
  assert.equal(droppedCount, 0);
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
  const { featureOverrides, droppedCount } = migrateRenderedIdFeatureEdits({
    features: session.features,
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
  const { featureOverrides, droppedCount } = migrateRenderedIdFeatureEdits({
    features: {
      featureVisibilityOverrides: {
        'f1_record_1__instance_4_0123456789abcdef': 'off',
        f1_record_1: 'on',
        f2_record_2: 'exclude_matching'
      }
    },
    catalog
  });
  assert.equal(droppedCount, 1);
  assert.deepEqual(Object.keys(featureOverrides).sort(), [key('seq-a', 'f1~4'), key('seq-b', 'f2')].sort());
  assert.equal(featureOverrides[key('seq-b', 'f2')].featureVisibility, 'exclude_matching');
});

test('a blank label text hid its label and becomes Label visibility Off', () => {
  const { featureOverrides } = migrateRenderedIdFeatureEdits({
    features: { labelTextFeatureOverrides: { f2: '  ' } },
    catalog: { schema: 4, items: [{ recordKeys: ['r'], features: [{ svgId: 'f2', recordKey: 'r', biologicalFeatureId: 'f2' }], biologicalFeatures: [] }] }
  });
  assert.deepEqual(featureOverrides[key('r', 'f2')], row('r', 'f2', { labelVisibility: 'off' }));
});
