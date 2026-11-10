// OV-3326: a TLOSATX genetic code that is not one of the NCBI tables LOSAT reads stops LOSAT with a trap. Generate
// refuses it before any search, with the field it was typed in; the tables are Python's (schema.generated.js).
import assert from 'node:assert/strict';
import test from 'node:test';

import { LOSAT_GENETIC_CODES, LOSAT_GENETIC_CODE_TEXT } from '../../gbdraw/web/js/contracts/schema.generated.js';
import { losatRecordGencode, requireLosatGeneticCode } from '../../gbdraw/web/js/services/linear-sources.js';
import { normalizeUserFacingError } from '../../gbdraw/web/js/utils/error-normalization.js';

const refusal = (run) => {
  try {
    run();
  } catch (error) {
    return normalizeUserFacingError(error, { operation: 'generate' });
  }
  assert.fail('the genetic code was not refused');
};

test('the tables are the ones Python names', () => {
  assert.deepEqual([...LOSAT_GENETIC_CODES], [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16,
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]);
  assert.equal(LOSAT_GENETIC_CODE_TEXT, '1–6, 9–16, 21–33');
});

test('a Linear card keeps a table LOSAT reads and its default', () => {
  assert.equal(losatRecordGencode({ losat_gencode: 11 }, 1), 11);
  assert.equal(losatRecordGencode({ losat_gencode: '' }, 1), null);
  assert.equal(losatRecordGencode({}, 1), null);
});

for (const code of [7, 99, 8.5]) test(`a Linear card's Gencode ${code} is refused with its card and field`, () => {
  const model = refusal(() => losatRecordGencode({ losat_gencode: code }, 2));
  assert.equal(model.code, 'COMPARISON_INPUT');
  assert.equal(model.context.reason, 'GENETIC_CODE');
  assert.match(model.summary, /Sequence 2\./);
  assert.doesNotMatch(model.summary, /Comparison sequence/);
  assert.match(model.summary, /Field: Gencode \(this entry\)\./);
  assert.match(model.summary, /Use the number of an NCBI genetic code that LOSAT reads\./);
  assert.match(model.summary, /Genetic codes LOSAT reads: 1–6, 9–16, 21–33\./);
});

test('a comparison ring and the Circular reference name their own fields', () => {
  const ring = refusal(() => requireLosatGeneticCode(99, { field: 'conservation_losat_gencodes', inputOrdinal: 3 }));
  assert.match(ring.summary, /Comparison sequence 3\./);
  assert.match(ring.summary, /Field: Comparison gencode\./);
  const reference = refusal(() => requireLosatGeneticCode(7, { field: 'reference_gencode' }));
  assert.match(reference.summary, /Field: Reference gencode\./);
  assert.equal(requireLosatGeneticCode(4, { field: 'reference_gencode' }), 4);
});
