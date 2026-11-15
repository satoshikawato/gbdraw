import assert from 'node:assert/strict';
import test from 'node:test';
import { RECORD_INDEX_SELECTOR_PATTERN } from '../../gbdraw/web/js/contracts/schema.generated.js';
import { SELECTOR_TEXT, featureSelector, hashSelector, recordIndexSelector } from '../../gbdraw/web/js/utils/selector-syntax.js';

// The Web writes record and feature selectors with the syntax Python reads (gbdraw/io/record_select.py,
// gbdraw/features/selector_values.py), and the help and error texts show the same forms (OV-2790).
test('the selector text the Web writes is the syntax Python reads', () => {
  assert.match(recordIndexSelector(3), new RegExp(`^(?:${RECORD_INDEX_SELECTOR_PATTERN})$`));
  assert.equal(recordIndexSelector(3), '#3');
  assert.equal(hashSelector('a1b2'), 'hash=a1b2');
  assert.equal(featureSelector('locus_tag', 'b0001'), 'locus_tag=b0001');
  assert.deepEqual({ ...SELECTOR_TEXT }, {
    recordNumber: '#<number>', featureSelectorExample: 'gene=COX1', featureSelectorForm: 'locus_tag=<tag>', hashSelector: 'hash='
  });
});
