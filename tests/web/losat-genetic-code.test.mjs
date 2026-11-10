// OV-3326: a TLOSATX genetic code that is not one of the NCBI tables LOSAT reads stops LOSAT with a trap. Generate
// refuses it before any search, with the field it was typed in; the tables are Python's (schema.generated.js).
import assert from 'node:assert/strict';
import test from 'node:test';

import { LOSAT_GENETIC_CODES, LOSAT_GENETIC_CODE_TEXT } from '../../gbdraw/web/js/contracts/schema.generated.js';
import {
  losatReadsGeneticCode, losatRecordGencode, requireLosatGeneticCode, searchedLosatGencode
} from '../../gbdraw/web/js/services/linear-sources.js';
import { comparisonCauseRepeats, normalizeUserFacingError } from '../../gbdraw/web/js/utils/error-normalization.js';

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
  const card = { inputOrdinal: 1, cardUid: 'linear-seq-a' };
  assert.equal(searchedLosatGencode({ losat_gencode: 11 }, card), 11);
  assert.equal(searchedLosatGencode({ losat_gencode: '' }, card), null);
  assert.equal(searchedLosatGencode({}, card), null);
  // The typed value as the job estimate reads it: no refusal, so the estimate counts what Generate would plan.
  assert.equal(losatRecordGencode({ losat_gencode: 7 }), 7);
});

// Review S-3: the File card's "Gencode required" mark and Generate's refusal read one rule.
test('one rule says whether LOSAT reads a code', () => {
  for (const code of LOSAT_GENETIC_CODES) assert.equal(losatReadsGeneticCode(code), true);
  for (const code of [0, 7, 8, 17, 34, 99, 8.5, Number.NaN]) assert.equal(losatReadsGeneticCode(code), false, String(code));
});

for (const code of [7, 99, 8.5]) test(`a Linear card's Gencode ${code} is refused with its card and field`, () => {
  const model = refusal(() => searchedLosatGencode({ losat_gencode: code }, { inputOrdinal: 2, cardUid: 'linear-seq-b' }));
  // The card's identity, whose LOSAT Gencode "Show in Comparison" focuses.
  assert.equal(model.context.cardUid, 'linear-seq-b');
  assert.equal(model.code, 'COMPARISON_INPUT');
  assert.equal(model.context.reason, 'GENETIC_CODE');
  assert.match(model.summary, /Sequence 2\./);
  assert.doesNotMatch(model.summary, /Comparison sequence/);
  assert.match(model.summary, /Field: LOSAT Gencode\./);
  // Only the action that leads to the cause (OV-3317 rule): Retry would refuse the same code again.
  assert.deepEqual(model.actions, ['edit-comparison']);
  assert.equal(comparisonCauseRepeats(model), true);
  assert.match(model.summary, /Use the number of an NCBI genetic code that LOSAT reads\./);
  assert.match(model.summary, /Genetic codes LOSAT reads: 1–6, 9–16, 21–33\./);
});

test('a comparison ring and the Circular reference name their own fields', () => {
  const ring = refusal(() => requireLosatGeneticCode(99, { field: 'conservation_losat_gencodes', inputOrdinal: 3 }));
  assert.match(ring.summary, /Comparison sequence 3\./);
  assert.match(ring.summary, /Field: Comparison gencode\./);
  const reference = refusal(() => requireLosatGeneticCode(7, { field: 'reference_gencode' }));
  assert.match(reference.summary, /Field: Reference gencode\./);
  // Circular shows "Show in Comparison" for both: Pairwise Comparisons holds the two fields.
  for (const model of [ring, reference]) {
    assert.deepEqual(model.actions, ['edit-comparison']);
    assert.equal(comparisonCauseRepeats(model), true);
  }
  assert.equal(requireLosatGeneticCode(4, { field: 'reference_gencode' }), 4);
});
