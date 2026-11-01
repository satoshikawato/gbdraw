// The Web side of the feature-start twin (tools/twins.json): the Rotate to
// feature "Start of the record" start base for every case of the shared
// vector, which tests/test_circular_origin_feature.py runs through the Python
// owner of gbdraw circular --origin_feature.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { resolveFeatureAnchor } from '../../gbdraw/web/js/app/record-display/feature-anchor.js';

const VECTORS = JSON.parse(readFileSync(new URL('../fixtures/feature-start-anchor-vectors.json', import.meta.url), 'utf8'));

for (const vector of VECTORS.cases) {
  test(`feature-start vector: ${vector.name}`, () => {
    const result = resolveFeatureAnchor({
      recordLength: vector.recordLength,
      effectiveCircular: true,
      cropped: false,
      currentReverseComplement: vector.reverseComplement,
      identity: { recordKey: 'record-1', biologicalFeatureId: 'feature-1' },
      parts: vector.parts,
      profile: vector.profile,
      intent: { placement: 'feature-start', anchor: null, offsetBp: 0, orientForward: false }
    });
    if (vector.refused) {
      assert.equal(result.eligibility.enabled, false);
      assert.equal(result.eligibility.code, vector.refused);
      assert.equal(result.startCoordinate, null);
    } else {
      assert.equal(result.eligibility.enabled, true, result.eligibility.message);
      assert.equal(result.startCoordinate, vector.startCoordinate);
      assert.equal(result.reverseComplement, vector.reverseComplement);
    }
  });
}
