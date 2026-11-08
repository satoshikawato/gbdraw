// One owner of a feature's first qualifier value (services/feature-utils.js `getQualifierFirstValue`): the
// feature popup, its hover summary, and the match popup's feature rows read the same value. A qualifier
// name that differs only in case matches, as Python's selectors match it (gbdraw/features/selector_values.py);
// the Web's feature rows carry lower-case names (gbdraw/web_support/feature_metadata.py), so for them the
// fallback never applies.
import assert from 'node:assert/strict';
import test from 'node:test';

const { getQualifierFirstValue, resolveDisplayProteinId } = await import('../../gbdraw/web/js/services/feature-utils.js');

test('the first value reads the field, then the qualifier, then a qualifier named in another case', () => {
  assert.equal(getQualifierFirstValue({ gene: 'dnaA', qualifiers: { gene: ['other'] } }, 'gene'), 'dnaA');
  assert.equal(getQualifierFirstValue({ gene: '', qualifiers: { gene: [' ', ' dnaA ', 'x'] } }, ' Gene '), 'dnaA');
  assert.equal(getQualifierFirstValue({ gene: null, qualifiers: { Gene: ['dnaA'] } }, 'gene'), 'dnaA');
  assert.equal(getQualifierFirstValue({ qualifiers: { gene: [], Gene: ['dnaA'] } }, 'gene'), '');
  assert.equal(getQualifierFirstValue({ gene: ['  '], qualifiers: { gene: ['dnaA'] } }, 'gene'), '');
  assert.equal(getQualifierFirstValue({ qualifiers: ['gene'] }, 'gene'), '');
  assert.equal(getQualifierFirstValue({ qualifiers: { gene: 'dnaA' } }, ''), '');
  assert.equal(getQualifierFirstValue(null, 'gene'), '');
  assert.equal(getQualifierFirstValue({ qualifiers: { gene: null, Gene: ['dnaA'] } }, 'gene'), 'dnaA');
  assert.equal(getQualifierFirstValue({ qualifiers: {} }, 'constructor'), '');
});

test('the protein ID display reads its qualifier through the same lookup', () => {
  assert.equal(resolveDisplayProteinId({ qualifiers: { Protein_ID: ['WP_000001.1'] } }), 'WP_000001.1');
  assert.equal(resolveDisplayProteinId({ qualifiers: { protein_id: null, PROTEIN_ID: ['WP_000002.1'] } }), 'WP_000002.1');
});
