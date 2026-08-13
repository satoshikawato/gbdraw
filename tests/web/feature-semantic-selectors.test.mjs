import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

import {
  encodeFeatureSemanticSelector,
  featureSemanticSelectorMatches,
  parseFeatureSemanticSelector
} from '../../gbdraw/web/js/app/feature-editor/semantic-fill-selectors.js';

const vectors = JSON.parse(await readFile(
  fileURLToPath(new URL('../fixtures/feature_semantic_selector_vectors.json', import.meta.url)),
  'utf8'
));

test('semantic selector literals match the shared Python/Web canonical vectors', () => {
  for (const vector of vectors) {
    assert.equal(encodeFeatureSemanticSelector(vector.kind, vector.value), vector.literal);
    assert.deepEqual(parseFeatureSemanticSelector(vector.literal), {
      kind: vector.kind,
      value: vector.value
    });
  }
  assert.equal(parseFeatureSemanticSelector('fs1:rendered-label:non canonical'), null);
  assert.equal(parseFeatureSemanticSelector('fs1:unknown:value'), null);
});

test('semantic matcher supports source/rendered label, caption, and every Similarity membership', () => {
  const feature = {
    type: 'CDS',
    qualifiers: { product: ['Shared annotation'] },
    renderedLabel: 'Shared rendered',
    orthogroupIds: ['og-alpha', 'og-beta']
  };
  assert.equal(featureSemanticSelectorMatches(
    feature,
    encodeFeatureSemanticSelector('source-annotation-label', 'Shared annotation')
  ), true);
  assert.equal(featureSemanticSelectorMatches(
    feature,
    encodeFeatureSemanticSelector('rendered-label', 'Shared rendered')
  ), true);
  for (const groupId of ['og-alpha', 'og-beta']) {
    assert.equal(featureSemanticSelectorMatches(
      feature,
      encodeFeatureSemanticSelector('similarity-group', groupId)
    ), true);
  }
  assert.equal(featureSemanticSelectorMatches(
    feature,
    encodeFeatureSemanticSelector('base-legend-caption', 'Core'),
    { baseLegendCaption: 'core' }
  ), true);
});
