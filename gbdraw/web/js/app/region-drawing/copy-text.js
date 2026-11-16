// @ts-check
// The words of Copy to another drawing (PD-OI-090 (4)): its refusals, the spans an existing drawing has already or
// leaves out, and what a paste does not bring. Each table is keyed by a vocabulary Python declares
// (`MODE_PROFILE_DATA.regionDrawing`, from gbdraw/api/region_drawing.py); the dialog never reads a message's text (R6).
import { msg, numberText, t, tm } from '../../utils/i18n.js';

/** @import { AlreadyPresentItem, LeftOutItem, PasteNotCarriedItem, RefusedCopy } from '../../contracts/locate.generated.js' */

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

/**
 * One line per part a paste does not bring (`pasteNotCarriedKinds`, Python's `PasteNotCarriedKind`), with its count.
 * @type {Readonly<Record<PasteNotCarriedItem['kind'], string>>}
 */
const PASTE_NOT_CARRIED_TEXT = Object.freeze({
  ANNOTATIONS: msg('Annotations of these records ({count})'),
  FEATURE_PLACEMENTS: msg('Feature placements of these records ({count})'),
  DEPTH: msg('Depth tracks of these records ({count})'),
  COMPARISON_TABLES: msg('The uploaded comparison tables ({count})'),
  CONSERVATION_RING_TABLES: msg('The uploaded ring tables ({count})'),
  NO_COMPARISON_ROWS: msg("Records the drawing's comparisons have no rows for ({count})")
});

/** @param {PasteNotCarriedItem} item @returns {string} */
export const pasteNotCarriedLine = ({ kind, count }) => tm(PASTE_NOT_CARRIED_TEXT[kind], { count: numberText(count) });

/**
 * A span as the dialog names it: the record's caption, then "whole record" or its source coordinates.
 * @param {string} caption @param {{ start: number, end: number }} span @param {boolean} whole
 */
const spanText = (caption, { start, end }, whole) => (whole
  ? t('{record} (whole record)', { record: caption })
  : t('{record} {start}–{end}', { record: caption, start: numberText(start), end: numberText(end) }));

/**
 * The line saying a chosen span is already in the existing drawing (Skip it / Add it again follows it).
 * @param {AlreadyPresentItem} item @param {(recordKey: string) => string} captionOf @param {string} drawingName
 */
export const alreadyPresentLine = (item, captionOf, drawingName) => t('{span} is already in {drawing}.', {
  span: spanText(captionOf(item.recordKey), item, item.whole), drawing: drawingName
});

/**
 * The line naming a span the existing drawing leaves out, with its reason (`leftOut`; one reason today).
 * @param {LeftOutItem} item @param {(recordKey: string) => string} captionOf
 */
export const leftOutLine = (item, captionOf) => t('{span} is not added: {reason}', {
  span: spanText(captionOf(item.recordKey), item, false), reason: copyRefusalLine(item.reason)
});
