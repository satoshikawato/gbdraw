// @ts-check
// The Figure size field, first in the Layout card of both modes (issue #1427). Each side is a number and a unit;
// an empty number is Auto. An edit writes the drawing's `form.figure_width` or `form.figure_height` (the amount as
// typed and the unit, services/figure-size-form.js); the control's History adapter records it as one step, and it
// applies on Generate (the Layout card's note). The unit of each side is view state.
/** @import { DrawingState } from '../state.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { FigureGrowthCause, FigureLengthUnit, FigureSide } from '../services/figure-size-form.js' */
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

/** @type {Record<FigureGrowthCause, string>} */
const GROWTH_CAUSES = {
  'legend-edits': msg('Legend edits'),
  'canvas-padding': msg('Canvas padding'),
  'moved-items': msg('moves of the Legend, title, or scale')
};

/** @param {NonNullable<ShownFigureSize['larger']>} larger */
const largerSentence = ({ sides, causes }) => {
  const size = listText(sides.map(({ side, length }) => (
    side === 'width' ? t('{length} width', { length }) : t('{length} height', { length })
  )));
  return causes.length === 0
    ? t('Larger than the requested {size}.', { size })
    : t('Larger than the requested {size} after {causes}, which the fit does not include.', {
      size, causes: listText(causes.map((cause) => tm(GROWTH_CAUSES[cause])))
    });
};

/** @param {FigureSizeFieldPorts} ports */
export const createFigureSizeField = ({ state, shownSize, watch }) => {
  const form = () => state.activeDrawing().form;
  /** @type {Record<FigureSide, { value: FigureLengthUnit }>} */
  const units = {
    width: ref(/** @type {FigureLengthUnit} */ (DEFAULT_FIGURE_LENGTH_UNIT)),
    height: ref(/** @type {FigureLengthUnit} */ (DEFAULT_FIGURE_LENGTH_UNIT))
  };
  // The value the field wrote last on each side, and the drawing it wrote it in. Any other value (a drawing switch,
  // Reset Settings, Load, Undo) shows in its own unit, and Auto in px; a side the field cleared keeps its unit.
  /** @type {Record<FigureSide, { drawingId: DrawingId, value: unknown } | null>} */
  const written = { width: null, height: null };
  for (const side of /** @type {const} */ (['width', 'height'])) {
    watch(() => [state.activeDrawing().id, form()[FIGURE_SIZE_FIELDS[side]]], () => {
      const drawing = state.activeDrawing();
      const value = drawing.form[FIGURE_SIZE_FIELDS[side]];
      const own = written[side];
      if (own !== null && own.drawingId === drawing.id && own.value === value) return;
      written[side] = null;
      units[side].value = figureLengthUnitOf(value);
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
    const drawing = state.activeDrawing();
    const value = figureLengthInputValue(amount, units[side].value);
    written[side] = { drawingId: drawing.id, value };
    if (part === 'unit' && amount === '') return;
    drawing.form[FIGURE_SIZE_FIELDS[side]] = value;
  };
  return {
    figureLengthUnits: FIGURE_LENGTH_UNITS,
    figureWidthUnit: units.width, figureHeightUnit: units.height, figureWidthAmount, figureHeightAmount,
    figureSizeLine, figureSizeFitting, figureSizeLarger, changeFigureSize
  };
};
