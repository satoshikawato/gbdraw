import assert from 'node:assert/strict';
import { test } from 'node:test';
import { buildLabelOverrideTsv } from '../../gbdraw/web/js/app/feature-editor/label-override-table.js';

test('an empty label edit does not traverse biological metadata to emit an empty table', () => {
  const result = buildLabelOverrideTsv({}, {}, {
    visibilityOverrides: { unknown: 'auto' },
    get extractedFeatures() { throw new Error('irrelevant biological traversal'); },
    get editableLabels() { throw new Error('irrelevant label traversal'); }
  });
  assert.deepEqual(result, { tsv: '', rows: [], skippedFeatureCount: 0,
    skippedFeatureSourceCount: 0, skippedMissingSourceCount: 0, fallbackHashCount: 0 });
});
