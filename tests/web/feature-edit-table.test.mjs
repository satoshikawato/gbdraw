// Export and Load Feature Edits TSV (design Q4 6.4, Owner Q2 = A and Q3 = A):
// the export reuses the Source recipe's --feature_override_table projection,
// and a load admits Python's rows and replaces the edits of the committed
// records only (R2).
import assert from 'node:assert/strict';
import test from 'node:test';

import { buildFeatureOverrideTable } from '../../gbdraw/web/js/app/run-info.js';
import {
  admitFeatureOverrideTable,
  replaceFeatureEdits
} from '../../gbdraw/web/js/app/feature-editor/feature-edit-table.js';
import { featureIdentityKey } from '../../gbdraw/web/js/services/feature-placement.js';

const record = (recordKey, cardinality = 'exactly_one') => ({
  recordKey, cardinality, source: { kind: 'genbank', resourceId: `source-${recordKey}` }
});
const row = (recordKey, biologicalFeatureId, values = {}) => ({
  recordKey, biologicalFeatureId, featureVisibility: null, labelVisibility: null, labelText: null, ...values
});
const HEADER = 'record\tfeature_selector\tfeature_visibility\tlabel_visibility\tlabel_text';

test('the export writes the Source recipe table: #index of the request record and hash=<identity>', async () => {
  const renderRequest = { records: [record('a'), record('b', 'all'), record('c')] };
  const counts = { 'source-b': 1 };
  const table = await buildFeatureOverrideTable({
    renderRequest,
    resources: {},
    rows: [row('a', 'f1', { labelText: ' Text "quoted" ' }), row('c', 'f2~1', { featureVisibility: 'off', labelVisibility: 'on' })],
    readResourceRecordCount: async (resourceId) => counts[resourceId]
  });
  assert.deepEqual(table, {
    text: `${HEADER}\n#1\thash=f1\t\t\t" Text ""quoted"" "\n#3\thash=f2~1\toff\ton\t\n`,
    reason: ''
  });
  // An input of several records shifts every later #index, so no table carries it.
  counts['source-b'] = 2;
  assert.deepEqual(await buildFeatureOverrideTable({
    renderRequest, resources: {}, rows: [row('c', 'f2', { featureVisibility: 'off' })],
    readResourceRecordCount: async (resourceId) => counts[resourceId]
  }), { text: '', reason: 'feature placements and feature edits require exact materialized records.' });
  assert.deepEqual(await buildFeatureOverrideTable({
    renderRequest: { records: [record('a')] }, resources: {}, rows: [row('b:2', 'f1', { featureVisibility: 'off' })]
  }), { text: '', reason: 'a feature edit names a feature identity that no table row can carry.' });
});

test('a loaded table must return request rows of the committed records', () => {
  const records = [record('a'), record('b', 'all')];
  const rows = [row('b:2', 'f2', { labelText: 'B' }), row('a', 'f1', { featureVisibility: 'off' })];
  assert.deepEqual(admitFeatureOverrideTable({ rows, unmatchedRows: [3] }, records),
    { rows: [rows[1], rows[0]], unmatchedRows: [3] });
  for (const result of [
    null,
    { rows },
    { rows, unmatchedRows: [1] },
    { rows: [row('other', 'f1', { featureVisibility: 'off' })], unmatchedRows: [] },
    { rows: [row('a', 'f1')], unmatchedRows: [] },
    { rows: [{ ...row('a', 'f1', { featureVisibility: 'off' }), labelSourceText: null }], unmatchedRows: [] }
  ]) {
    assert.throws(() => admitFeatureOverrideTable(result, records), (error) => (
      error.code === 'RESULT_INVALID' && error.operation === 'readFeatureOverrideTable'
    ));
  }
});

test('a load replaces the edits of the committed records and keeps the other records and label source texts', () => {
  // The Linear drawing's rows; the Circular drawing's rows are its own map (R2; PR-1).
  const draft = {};
  const put = (value) => {
    draft[featureIdentityKey(value.recordKey, value.biologicalFeatureId)] = value;
  };
  put({ ...row('a', 'old', { featureVisibility: 'off' }), labelSourceText: 'tRNA-Leu' });
  put({ ...row('a', 'gone', { labelVisibility: 'on' }), labelSourceText: null });
  put({ ...row('linear-other', 'f1', { labelText: 'Other record' }), labelSourceText: null });
  const circular = { [featureIdentityKey('a', 'gone')]: { ...row('a', 'gone', { labelText: 'Same key' }), labelSourceText: null } };
  const circularBefore = structuredClone(circular);
  replaceFeatureEdits(draft, [row('a', 'old', { labelText: 'New' }), row('a', 'f9', { featureVisibility: 'exclude_matching' })],
    [record('a')]);
  assert.deepEqual(Object.values(draft).map((value) => [value.recordKey, value.biologicalFeatureId,
    value.featureVisibility, value.labelVisibility, value.labelText, value.labelSourceText]).sort(), [
    ['a', 'f9', 'exclude_matching', null, null, null],
    ['a', 'old', null, null, 'New', 'tRNA-Leu'],
    ['linear-other', 'f1', null, null, 'Other record', null]
  ]);
  assert.deepEqual(circular, circularBefore);
});
