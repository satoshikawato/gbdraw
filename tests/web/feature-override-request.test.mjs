// Per-feature edits reach the render request as identity rows (request schema 9
// `featureOverrides`, design Q4 3.2, 6.1); the label table carries rules only.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { buildBulkLabelProjection } from '../../gbdraw/web/js/services/label-override-table.js';
import {
  canonicalFeatureOverrides,
  labelTextSettingsVisible,
  requestFeatureOverrides,
  updateFeatureOverride
} from '../../gbdraw/web/js/services/feature-placement.js';

// A drawing's draft key is the identity pair (PR-1: the mode is the drawing's).
const key = (recordKey, featureId) => JSON.stringify([recordKey, featureId]);
const feature = (recordKey, featureId) => ({ scope: 'linear', record_key: recordKey, biological_feature_id: featureId });

test('without label edits the projection is empty and reads no input', () => {
  const result = buildBulkLabelProjection({}, {
    get labelTargets() { throw new Error('irrelevant label traversal'); },
    get featureOverrides() { throw new Error('irrelevant edit traversal'); }
  });
  assert.deepEqual(result, { bulkLabelText: {}, rows: [] });
});

test('the request carries the rows of its records without the Web-only source text (R2)', () => {
  const draft = {};
  updateFeatureOverride(draft, feature('seq-b', 'f9'), { featureVisibility: 'off' });
  updateFeatureOverride(draft, feature('seq-a', 'f1~2'), { labelText: 'Renamed', labelSourceText: 'source' });
  updateFeatureOverride(draft, feature('circular-x', 'f1'), { labelVisibility: 'on' });
  updateFeatureOverride(draft, feature('seq-a', 'f3'), { labelSourceText: 'kept for a bulk edit' });
  assert.deepEqual(requestFeatureOverrides(draft, [{ recordKey: 'seq-a' }, { recordKey: 'seq-b' }]), [
    { recordKey: 'seq-a', biologicalFeatureId: 'f1~2', featureVisibility: null, labelVisibility: null, labelText: 'Renamed' },
    { recordKey: 'seq-b', biologicalFeatureId: 'f9', featureVisibility: 'off', labelVisibility: null, labelText: null }
  ]);
  // An ALL record owns its `<recordKey>:<n>` expansions.
  updateFeatureOverride(draft, feature('all-input:2', 'f4'), { featureVisibility: 'on' });
  assert.equal(requestFeatureOverrides(draft, [{ recordKey: 'all-input', cardinality: 'all' }])[0].recordKey,
    'all-input:2');
});

// B6: a bulk label edit (a global `label` row of a Label TSV) reaches the
// labels of every batch Result: the displayed ones and the recorded sources.
test('a bulk label edit reaches the recorded source of a feature on another Result', () => {
  const draft = {};
  updateFeatureOverride(draft, feature('record-2', 'f6'), { labelSourceText: 'gtg start' });
  updateFeatureOverride(draft, feature('record-3', 'f6'), { labelSourceText: 'gtg start', labelText: 'OWN' });
  const { bulkLabelText, rows } = buildBulkLabelProjection({ 'gtg start': 'BULK' }, {
    labelTargets: [{ identityKey: key('record-1', 'f6'), sourceText: 'gtg start' }],
    featureOverrides: draft
  });
  assert.deepEqual(rows, []);
  const records = ['record-1', 'record-2', 'record-3'].map((recordKey) => ({ recordKey }));
  assert.deepEqual(requestFeatureOverrides(draft, records, { bulkLabelText })
    .map((row) => [row.recordKey, row.labelText]), [
    ['record-1', 'BULK'], ['record-2', 'BULK'], ['record-3', 'OWN']
  ]);
});

test('request rows are one non-blank line and at least one edit', () => {
  assert.throws(() => canonicalFeatureOverrides([{
    recordKey: 'r', biologicalFeatureId: 'f', featureVisibility: null, labelVisibility: null, labelText: null
  }]));
  assert.throws(() => canonicalFeatureOverrides([{
    recordKey: 'r', biologicalFeatureId: 'f', featureVisibility: null, labelVisibility: null, labelText: 'a\tb'
  }]));
  assert.throws(() => canonicalFeatureOverrides([{
    recordKey: 'r', biologicalFeatureId: 'f', featureVisibility: 'hidden', labelVisibility: null, labelText: null
  }]));
  // A draft key encodes its row's identity.
  assert.throws(() => canonicalFeatureOverrides({
    [key('r', 'other')]: {
      recordKey: 'r', biologicalFeatureId: 'f', featureVisibility: 'off',
      labelVisibility: null, labelText: null, labelSourceText: null
    }
  }));
});

test('the label text settings show while the label scope is not None or a feature label is On (OV-222)', () => {
  const labelRow = (labelVisibility) => ({
    recordKey: 'r1', biologicalFeatureId: 'f1', featureVisibility: null, labelVisibility, labelText: null, labelSourceText: null
  });
  assert.equal(labelTextSettingsVisible('none', {}), false);
  assert.equal(labelTextSettingsVisible('none', { [key('r1', 'f1')]: labelRow('off') }), false);
  assert.equal(labelTextSettingsVisible('none', { [key('r1', 'f1')]: labelRow('on') }), true);
  for (const scope of ['all', 'first', 'orthogroup_top', 'out', 'both']) {
    assert.equal(labelTextSettingsVisible(scope, {}), true, scope);
  }
});
