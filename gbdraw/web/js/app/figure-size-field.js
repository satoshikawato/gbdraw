// @ts-check
// The Figure size field, first in the Layout card of both modes (issue #1427). Each side is a number and a unit;
// an empty number is Auto. An edit writes the drawing's `form.figure_width` or `form.figure_height` (the amount as
// typed and the unit, services/figure-size-form.js); the control's History adapter records it as one step, and it
// applies on Generate (the Layout card's note). The unit of each side is view state.
/** @import { DrawingState } from '../state.js' */
/** @import { FigureGrowthCause, FigureLengthFormValue, FigureLengthUnit, FigureSide } from '../services/figure-size-form.js' */
import {
  DEFAULT_FIGURE_LENGTH_UNIT, FIGURE_LENGTH_UNITS, FIGURE_SIZE_FIELDS, figureLengthAmountText, figureLengthInputValue,
  figureLengthUnitOf
} from '../services/figure-size-form.js';
import { listText, msg, t, tm } from '../utils/i18n.js';
import { computed, ref } from '../utils/reactivity.js';

/**
 * @typedef {object} FigureSizeFieldPorts
 * @property {{ activeDrawing: () => DrawingState }} state
 * @property {() => ShownFigureSize | null} shownSize The shown Result's size, null while the drawing shows no
 *   Result.
 * @property {(source: () => unknown, callback: () => void) => unknown} watch Vue `watch`: reads a value from
 *   outside the field (Load, Undo, Reset, a drawing switch).
 */

/**
 * The shown Result's size (`shownFigureSize`), and the requested sides it is larger than with the Web-only
 * additions that made it so (`figureGrowth`; null when it is not larger). `length` names a requested side: `180 mm`.
 * `fitting`: the automatic rerender is fitting the figure, drawn to a requested size, again (after a Legend edit or
 * Canvas padding, PD-OI-114 rev 2); until it ends the Result shows the edit as the Web drew it at once.
 * @typedef {{ first: string, second: string, fitting: boolean,
 *   larger: { sides: { side: FigureSide, length: string }[], causes: FigureGrowthCause[] } | null }} ShownFigureSize
 */

/** Which part of a side an edit comes from. @typedef {'amount' | 'unit'} FigureSizeFieldPart */

const MOVED_ITEMS = msg('moves of the Legend, title, or scale');
/** The causes Python fits a figure to (PD-OI-114 rev 2): when they leave it larger, the fit was refused. */
const FITTED_CAUSES = /** @type {readonly FigureGrowthCause[]} */ (['legend-edits', 'canvas-padding']);

/**
 * Why the shown figure is larger than requested. Python fits a figure drawn to a requested size to its Legend edits
 * and Canvas padding, so either cause left after the automatic rerender means the fit was refused, which the
 * rerender's note above the preview explains (review S3, D-60); a move of the Legend, title, or scale is never fitted.
 * @param {NonNullable<ShownFigureSize['larger']>} larger
 */
const largerSentence = ({ sides, causes }) => {
  const size = listText(sides.map(({ side, length }) => (
    side === 'width' ? t('{length} width', { length }) : t('{length} height', { length })
  )));
  if (causes.some((cause) => FITTED_CAUSES.includes(cause))) {
    return t('Larger than the requested {size}: the edit could not be fitted (see the note above the preview).', { size });
  }
  return causes.includes('moved-items')
    ? t('Larger than the requested {size} after {causes}, which the fit does not include.', { size, causes: tm(MOVED_ITEMS) })
    : t('Larger than the requested {size}.', { size });
};

/** @param {FigureSizeFieldPorts} ports */
export const createFigureSizeField = ({ state, shownSize, watch }) => {
  const form = () => state.activeDrawing().form;
  /** @type {Record<FigureSide, { value: FigureLengthUnit }>} */
  const units = {
    width: ref(/** @type {FigureLengthUnit} */ (DEFAULT_FIGURE_LENGTH_UNIT)),
    height: ref(/** @type {FigureLengthUnit} */ (DEFAULT_FIGURE_LENGTH_UNIT))
  };
  // The value the field just wrote on each side, until the watch below reads it back once. Any other change of a
  // side (a drawing switch, Reset Settings, Load, Undo) shows in its own unit, and Auto in px; so an Auto side whose
  // unit the user picked without an amount returns to px when Reset Settings writes Auto again, while a side the
  // field cleared keeps its unit.
  /** @type {Record<FigureSide, { value: FigureLengthFormValue } | null>} */
  const written = { width: null, height: null };
  for (const side of /** @type {const} */ (['width', 'height'])) {
    watch(() => [state.activeDrawing().id, form()[FIGURE_SIZE_FIELDS[side]]], () => {
      const value = form()[FIGURE_SIZE_FIELDS[side]];
      const own = written[side];
      written[side] = null;
      if (own === null || own.value !== value) units[side].value = figureLengthUnitOf(value);
    });
  }
  const figureWidthAmount = computed(() => figureLengthAmountText(form().figure_width));
  const figureHeightAmount = computed(() => figureLengthAmountText(form().figure_height));
  const amounts = { width: figureWidthAmount, height: figureHeightAmount };
  // One line under the field: the shown Result's size; none while the drawing shows no Result. Under it, that the
  // figure is being fitted while the rerender runs, else why the shown Result is larger than requested, when it is.
  const figureSizeLine = computed(() => {
    const size = shownSize();
    return size === null ? '' : t('Figure: {first} ({second})', { first: size.first, second: size.second });
  });
  const figureSizeFitting = computed(() => (shownSize()?.fitting ? t('Fitting the figure to the requested size.') : ''));
  const figureSizeLarger = computed(() => {
    const size = shownSize();
    return size?.larger && !size.fitting ? largerSentence(size.larger) : '';
  });
  /**
   * @param {Event} event A change of a side's number or unit.
   * @param {FigureSide} side
   * @param {FigureSizeFieldPart} part
   */
  const changeFigureSize = (event, side, part) => {
    const control = /** @type {HTMLInputElement | HTMLSelectElement} */ (event.target);
    if (part === 'unit' && FIGURE_LENGTH_UNITS.includes(control.value)) units[side].value = control.value;
    // A unit change keeps the amount shown (180 mm becomes 180 in); with no amount it writes nothing, and the Auto
    // side keeps the unit chosen.
    const amount = part === 'amount' ? control.value : amounts[side].value;
    if (part === 'unit' && amount === '') return;
    const value = figureLengthInputValue(amount, units[side].value);
    const field = FIGURE_SIZE_FIELDS[side];
    // A write that changes nothing triggers no watch to read it back.
    if (form()[field] === value) return;
    written[side] = { value };
    form()[field] = value;
  };
  return {
    figureLengthUnits: FIGURE_LENGTH_UNITS,
    figureWidthUnit: units.width, figureHeightUnit: units.height, figureWidthAmount, figureHeightAmount,
    figureSizeLine, figureSizeFitting, figureSizeLarger, changeFigureSize
  };
};
