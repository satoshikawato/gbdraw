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

// B6: a bulk label edit (a global `label` row of a Label TSV) reaches the
// labels of every batch Result at Generate, not only the displayed ones.
test('a bulk label edit expands to the recorded source of a feature on another Result', () => {
  const features = ['TESTA_0006', 'TESTB_0006'].map((locusTag) => ({
    svg_id: `f_${locusTag}`,
    record_id: locusTag.slice(0, 5),
    type: 'CDS',
    qualifiers: { locus_tag: [locusTag] },
    selector: { hash: `f_${locusTag}`, qualifiers: { locus_tag: [locusTag] } }
  }));
  const { rows } = buildLabelOverrideTsv({}, { 'gtg start': 'BULK' }, {
    extractedFeatures: features,
    editableLabels: [{ featureId: 'f_TESTA_0006', sourceText: 'gtg start', text: 'BULK' }],
    featureOverrideSources: { f_TESTA_0006: 'gtg start', f_TESTB_0006: 'gtg start' }
  });
  assert.deepEqual(rows, [
    'TESTA\tCDS\tlocus_tag\t^TESTA_0006$\tBULK',
    'TESTB\tCDS\tlocus_tag\t^TESTB_0006$\tBULK'
  ]);
});
