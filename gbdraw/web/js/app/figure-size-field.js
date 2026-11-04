// @ts-check
// The Figure size field, first in the Layout card of both modes (issue #1427). Each side is a number and a unit;
// an empty number is Auto. An edit writes the drawing's `form.figure_width` or `form.figure_height` (the amount as
// typed and the unit, services/figure-size-form.js); the control's History adapter records it as one step, and it
// applies on Generate (the Layout card's note). The unit of each side is view state.
/** @import { DrawingState } from '../state.js' */
/** @import { FigureLengthUnit, FigureSide } from '../services/figure-size-form.js' */
import {
  DEFAULT_FIGURE_LENGTH_UNIT, FIGURE_LENGTH_UNITS, FIGURE_SIZE_FIELDS, figureLengthAmountText, figureLengthInputValue,
  figureLengthUnitOf
} from '../services/figure-size-form.js';
import { t } from '../utils/i18n.js';
import { computed, ref } from '../utils/reactivity.js';

/**
 * @typedef {object} FigureSizeFieldPorts
 * @property {{ activeDrawing: () => DrawingState }} state
 * @property {() => { first: string, second: string } | null} shownSize The shown Result's size
 *   (`shownFigureSize`), null while the drawing shows no Result.
 * @property {(source: () => unknown, callback: () => void) => unknown} watch Vue `watch`: reads a value from
 *   outside the field (Load, Undo, Reset, a drawing switch).
 */

/** Which part of a side an edit comes from. @typedef {'amount' | 'unit'} FigureSizeFieldPart */

/** @param {FigureSizeFieldPorts} ports */
export const createFigureSizeField = ({ state, shownSize, watch }) => {
  const form = () => state.activeDrawing().form;
  /** @type {Record<FigureSide, { value: FigureLengthUnit }>} */
  const units = {
    width: ref(/** @type {FigureLengthUnit} */ (DEFAULT_FIGURE_LENGTH_UNIT)),
    height: ref(/** @type {FigureLengthUnit} */ (DEFAULT_FIGURE_LENGTH_UNIT))
  };
  // The value the field wrote last on each side: a value from anywhere else shows in its own unit.
  /** @type {Record<FigureSide, unknown>} */
  const written = { width: undefined, height: undefined };
  for (const side of /** @type {const} */ (['width', 'height'])) {
    watch(() => [state.activeDrawing().id, form()[FIGURE_SIZE_FIELDS[side]]], () => {
      const value = form()[FIGURE_SIZE_FIELDS[side]];
      if (value !== written[side]) units[side].value = figureLengthUnitOf(value, units[side].value);
      written[side] = undefined;
    });
  }
  const figureWidthAmount = computed(() => figureLengthAmountText(form().figure_width));
  const figureHeightAmount = computed(() => figureLengthAmountText(form().figure_height));
  const amounts = { width: figureWidthAmount, height: figureHeightAmount };
  // One line under the field: the shown Result's size; none while the drawing shows no Result.
  const figureSizeLine = computed(() => {
    const size = shownSize();
    return size === null ? '' : t('Figure: {first} ({second})', size);
  });
  /**
   * @param {Event} event A change of a side's number or unit.
   * @param {FigureSide} side
   * @param {FigureSizeFieldPart} part
   */
  const changeFigureSize = (event, side, part) => {
    const control = /** @type {HTMLInputElement | HTMLSelectElement} */ (event.target);
    if (part === 'unit' && FIGURE_LENGTH_UNITS.includes(control.value)) units[side].value = control.value;
    // A unit change keeps the amount shown (180 mm becomes 180 in); with no amount it writes nothing.
    const amount = part === 'amount' ? control.value : amounts[side].value;
    if (part === 'unit' && amount === '') return;
    const value = figureLengthInputValue(amount, units[side].value);
    written[side] = value;
    state.activeDrawing().form[FIGURE_SIZE_FIELDS[side]] = value;
  };
  return {
    figureLengthUnits: FIGURE_LENGTH_UNITS,
    figureWidthUnit: units.width, figureHeightUnit: units.height, figureWidthAmount, figureHeightAmount,
    figureSizeLine, changeFigureSize
  };
};
