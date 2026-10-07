import assert from 'node:assert/strict';
import {
  keepResultNames,
  normalizeLogicalResults
} from '../../gbdraw/web/js/services/result-normalization.js';

assert.deepEqual(
  normalizeLogicalResults([
    { name: 'layout.svg', content: 'plain' },
    { name: 'layout.interactive.svg', content: 'interactive' }
  ]),
  [{ name: 'layout.svg', content: 'plain' }]
);

assert.deepEqual(
  normalizeLogicalResults([
    { name: 'layout.interactive.svg', content: 'interactive' },
    { name: 'layout.svg', content: 'plain' },
    { name: 'layout_2.svg', content: 'batch-2' }
  ]),
  [
    { name: 'layout.svg', content: 'plain' },
    { name: 'layout_2.svg', content: 'batch-2' }
  ]
);

assert.deepEqual(
  normalizeLogicalResults([
    { name: 'interactive-map.svg', content: 'distinct' },
    { name: 'interactive-map.interactive.svg', content: 'paired' },
    { name: 'only.interactive.svg', content: 'legacy-only' }
  ]),
  [
    { name: 'interactive-map.svg', content: 'distinct' },
    { name: 'only.interactive.svg', content: 'legacy-only' }
  ]
);

assert.deepEqual(normalizeLogicalResults(null), []);

// OV-136: a rerender's Worker reply names its Results after the output prefix;
// the Results it draws again keep their names, in every row that names them.
const rerenderReply = () => ({
  results: [{ name: 'out_1.svg', content: 'a' }, { name: 'out_2.svg', content: 'b' }],
  metadata: {
    featureCatalog: { schema: 3, items: [
      { resultIndex: 0, resultName: 'out_1.svg' }, { resultIndex: 1, resultName: 'out_2.svg' }
    ] },
    legendRows: [
      { resultIndex: 0, resultName: 'out_1.svg', drawn: [] }, { resultIndex: 1, resultName: 'out_2.svg', drawn: [] }
    ],
    trackSlotGeometry: { schema: 1, records: [{ resultIndex: 1, resultName: 'out_2.svg', slots: [] }] },
    annotationWarnings: [{ code: 'empty_span', resultIndex: 1, resultName: 'out_2.svg' }],
    comparisonWarnings: [{ code: 'comparison_record_id_unmatched', resultIndex: 0, resultName: 'out_1.svg' }],
    featureIdentityNotices: [{ recordKey: 'record-1', resultIndex: 0 }]
  }
});
const kept = rerenderReply();
keepResultNames(kept, ['HmmtDNA_1', 'HmmtDNA_2']);
assert.deepEqual(kept, {
  results: [{ name: 'HmmtDNA_1', content: 'a' }, { name: 'HmmtDNA_2', content: 'b' }],
  metadata: {
    featureCatalog: { schema: 3, items: [
      { resultIndex: 0, resultName: 'HmmtDNA_1' }, { resultIndex: 1, resultName: 'HmmtDNA_2' }
    ] },
    legendRows: [
      { resultIndex: 0, resultName: 'HmmtDNA_1', drawn: [] }, { resultIndex: 1, resultName: 'HmmtDNA_2', drawn: [] }
    ],
    trackSlotGeometry: { schema: 1, records: [{ resultIndex: 1, resultName: 'HmmtDNA_2', slots: [] }] },
    annotationWarnings: [{ code: 'empty_span', resultIndex: 1, resultName: 'HmmtDNA_2' }],
    comparisonWarnings: [{ code: 'comparison_record_id_unmatched', resultIndex: 0, resultName: 'HmmtDNA_1' }],
    featureIdentityNotices: [{ recordKey: 'record-1', resultIndex: 0 }]
  }
});
// Another number of Results, a blank name, or no names: the engine's names stay.
for (const names of [['HmmtDNA_1'], ['HmmtDNA_1', ' '], null]) {
  const reply = rerenderReply();
  keepResultNames(reply, names);
  assert.deepEqual(reply, rerenderReply());
}
