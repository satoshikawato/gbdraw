// @ts-check
// The words of Copy to another drawing (PD-OI-090 (4)): its refusals, keyed by the vocabulary Python declares
// (`MODE_PROFILE_DATA.regionDrawing.copyRefusalReasons`, from gbdraw/api/region_drawing.py); the dialog never reads a
// message's text (R6).
import { msg, tm } from '../../utils/i18n.js';

/** @import { RefusedCopy } from '../../contracts/locate.generated.js' */

/**
 * One line per refusal of a copy (`copyRefusalReasons`, Python's `CopyRefusalReason`).
 * @type {Readonly<Record<RefusedCopy['reason'], string>>}
 */
const COPY_REFUSAL_TEXT = Object.freeze({
  CIRCULAR_WHOLE_RECORDS: msg('A Circular drawing draws several whole records, or one record or range alone.'),
  CIRCULAR_RANGE_TARGET: msg('That Circular drawing draws a range; it takes no other record or range.'),
  SAME_DRAWING: msg('Choose a drawing other than the shown one.')
});

/** @param {RefusedCopy['reason']} reason @returns {string} */
export const copyRefusalLine = (reason) => tm(COPY_REFUSAL_TEXT[reason]);
