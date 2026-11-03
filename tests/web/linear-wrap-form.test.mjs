// The Web's Wrap Rows: the words of its refusals and of the export limit (one wording owner,
// utils/error-normalization.js), and the channel of its value.
import assert from 'node:assert/strict';
import test from 'node:test';

import { normalizeUserFacingError, reasonText } from '../../gbdraw/web/js/utils/error-normalization.js';

// What Python sends for each wrap refusal (gbdraw/layout/linear_wrap.py; tests/test_web_wrap_diagnostics.py).
const WRAP_REFUSALS = Object.freeze({
  WRAP: ['COMPARISON_INPUT', 'Remove the comparison inputs (BLAST, comparison pairs, LOSAT, Similarity groups), or turn Wrap Rows off.'],
  WRAP_TRACK_LAYOUT: ['INPUT_INVALID', 'Set Track Layout to Above or Below, move the feature track off the axis in Custom Track Slots, or turn Wrap Rows off.'],
  WRAP_SCALE_STYLE: ['INPUT_INVALID', 'Set Scale Style to Ruler (Ticks), or turn Wrap Rows off.'],
  WRAP_NORMALIZE_LENGTH: ['INPUT_INVALID', 'Turn Normalize Record Lengths off, or turn Wrap Rows off.'],
  WRAP_GRID_ROW: ['INPUT_INVALID', 'Give each record its own row, or turn Wrap Rows off.'],
  WRAP_ANCHOR_FEATURE: ['INPUT_INVALID', 'Reset the alignment, or turn Wrap Rows off.']
});

test('each wrap refusal names what rows cannot draw and the way out, the text its hint shows', () => {
  for (const [reason, [code, way]] of Object.entries(WRAP_REFUSALS)) {
    const text = reasonText(reason);
    assert.ok(text?.endsWith(way), reason);
    assert.doesNotMatch(text, /Preserved session settings|canvas\.linear\.wrap_length/, `${reason}: the retired advice`);
    const model = normalizeUserFacingError({ code, operation: 'generate', stage: 'render', context: { field: 'wrap_length', reason } });
    assert.deepEqual(model.context, { field: 'wrap_length', reason }, reason);
    assert.ok(model.summary.endsWith(`Field: Wrap Rows. ${text}`), `${reason}: ${model.summary}`);
  }
  assert.equal(reasonText('NOT_A_REASON'), null);
  assert.equal(reasonText('toString'), null, 'only the owner\'s own reasons');
});

test('an export over the format limit names the limit, the size and what would fit, and offers SVG', () => {
  const png = normalizeUserFacingError({ code: 'EXPORT_LIMIT', operation: 'export-png', stage: 'export-capture',
    context: { reason: 'PNG_SIZE', limit: 32767, width: 6250, height: 143647, rowLength: 500000, dpi: 68 } });
  assert.equal(png.summary, 'The figure is larger than this file format allows, so no file was made.'
    + ' PNG is limited to 32,767 px per side, and this figure is 6,250 x 143,647 px.'
    + ' A row length of 500,000 bp should fit (estimated from the drawn height).'
    + ' A DPI of 68 or lower would fit.'
    + ' Lower the PNG DPI, or use SVG, which has no size limit.');
  assert.deepEqual(png.actions, ['use-svg']);
  const pdf = normalizeUserFacingError({ code: 'EXPORT_LIMIT', operation: 'export-pdf', stage: 'export-capture',
    context: { reason: 'PDF_SIZE', limit: 14400, width: 1500, height: 34475 } });
  assert.equal(pdf.summary, 'The figure is larger than this file format allows, so no file was made.'
    + ' PDF is limited to 14,400 pt per side, and this figure is 1,500 x 34,475 pt.'
    + ' Set a longer row length in Wrap Rows, or use SVG, which has no size limit.');
  const unbounded = normalizeUserFacingError({ code: 'EXPORT_LIMIT', operation: 'export-png', stage: 'export-capture',
    context: { reason: 'PNG_SIZE', limit: 32767, width: -1, height: 1.5, dpi: 2e7 } });
  assert.deepEqual(unbounded.context, { reason: 'PNG_SIZE', limit: 32767 }, 'only bounded counts are kept');
  assert.equal(unbounded.summary, 'The figure is larger than this file format allows, so no file was made.'
    + ' Lower the PNG DPI, or use SVG, which has no size limit.');
});
