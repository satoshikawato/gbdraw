// @ts-check
// The figure size of a drawing's form (`form.figure_width`, `form.figure_height`): what the form holds, how a
// stored value reads into it, and how the shown Result's size is written under the field. The request leaves are
// `canvas.figure_width` and `canvas.figure_height` (services/session-request.js compiles them); Python parses and
// checks a length (gbdraw/core/figure_length.py) and fits the figure to it. The Web never converts or validates a
// requested length (R7): the units and their px come from Python through tools/generate_mode_profiles.py.
/** @import { CanonicalRenderRequest } from '../contracts/request.generated.js' */
/** @import { StoredLeafReadings } from './stored-leaf-readings.js' */
import { MODE_PROFILE_DATA } from '../mode-profiles.generated.js';
import { NO_STORED_READINGS, storedLeafKey } from './stored-leaf-readings.js';

/** The request leaves of the figure size (config override paths). */
export const FIGURE_WIDTH_PATH = 'canvas.figure_width';
export const FIGURE_HEIGHT_PATH = 'canvas.figure_height';

/**
 * A unit of a figure length, as Python lists them (`FigureLengthUnit`).
 * @typedef {string} FigureLengthUnit
 */
/**
 * Python's unit table (`figure_length_unit_table` in gbdraw/core/figure_length.py).
 * @typedef {{ units: readonly FigureLengthUnit[], pxPerUnit: Readonly<Record<FigureLengthUnit, number>>,
 *   messageDecimals: Readonly<Record<FigureLengthUnit, number>> }} FigureLengthUnitTable
 */
const UNIT_TABLE = /** @type {FigureLengthUnitTable} */ (MODE_PROFILE_DATA.figureLength);
/** The units of the field's unit selects, in Python's order. */
export const FIGURE_LENGTH_UNITS = UNIT_TABLE.units;
/** The unit an empty side shows (the Owner, 2026-11-04). */
export const DEFAULT_FIGURE_LENGTH_UNIT = /** @type {FigureLengthUnit} */ ('px');

/**
 * `null`: Auto; otherwise a length as the field wrote it, an amount and a unit (`"180mm"`), or Python's canonical
 * text of a stored length.
 * @typedef {null | string} FigureLengthFormValue
 */
/**
 * A stored value read into the form: `lost` when neither the form nor Python can read it, and then the side is Auto.
 * @typedef {{ value: FigureLengthFormValue, lost: boolean }} FigureLengthReading
 */
/** The form fields of the figure size, by side. */
export const FIGURE_SIZE_FIELDS = Object.freeze(/** @type {const} */ ({ width: 'figure_width', height: 'figure_height' }));
/** @typedef {keyof typeof FIGURE_SIZE_FIELDS} FigureSide */
/** The config path of each side. */
export const FIGURE_SIZE_PATHS = Object.freeze(/** @type {const} */ ({ width: FIGURE_WIDTH_PATH, height: FIGURE_HEIGHT_PATH }));

// The text a number field gives (an HTML valid floating-point number), which the field writes before its unit.
const AMOUNT_TEXT = /^-?(?:\d+(?:\.\d+)?|\.\d+)(?:[eE][-+]?\d+)?$/;

/**
 * The amount and unit of a length the form can hold: an amount a number field can show and one of Python's unit
 * names, with nothing between them. Null for any other value, which only Python reads.
 * @param {unknown} value
 * @returns {{ amount: string, unit: FigureLengthUnit } | null}
 */
export const splitFigureLength = (value) => {
  if (typeof value !== 'string') return null;
  const unit = FIGURE_LENGTH_UNITS.find((name) => value.endsWith(name));
  if (unit === undefined) return null;
  const amount = value.slice(0, -unit.length);
  return AMOUNT_TEXT.test(amount) ? { amount, unit } : null;
};

/** @param {unknown} raw */
const isAuto = (raw) => raw === undefined || raw === null || (typeof raw === 'string' && raw.toLowerCase() === 'auto');

/**
 * Reads a stored figure length (a Session's form value, a request or config leaf) as the form holds it: absent,
 * `null` and `"auto"` (any case) are Auto, and a length the form can hold (`splitFigureLength`) reads as it is. Any
 * other value is Python's to read: its reading in `readings`, else lost.
 * @param {FigureSide} side
 * @param {unknown} raw
 * @param {StoredLeafReadings} [readings]
 * @returns {FigureLengthReading}
 */
export const figureLengthFormValue = (side, raw, readings = NO_STORED_READINGS) => {
  if (isAuto(raw)) return { value: null, lost: false };
  if (splitFigureLength(raw) !== null) return { value: /** @type {string} */ (raw), lost: false };
  const canonical = readings.get(storedLeafKey(FIGURE_SIZE_PATHS[side], raw));
  if (isAuto(canonical) && canonical !== undefined) return { value: null, lost: false };
  return splitFigureLength(canonical) === null ? { value: null, lost: true } : { value: /** @type {string} */ (canonical), lost: false };
};

/**
 * The figure size a committed request stores, each side as stored: its override leaf, or else its full config's (a
 * CLI request); null for a side it leaves Auto.
 * @param {CanonicalRenderRequest} request
 * @returns {Record<FigureSide, unknown>}
 */
export const requestFigureSizeLeaves = (request) => {
  const { configOverrides, config } = /** @type {{ configOverrides?: Record<string, unknown> | null,
    config?: { canvas?: Record<string, unknown> } | null }} */ (request.diagramOptions);
  /** @param {FigureSide} side */
  const leaf = (side) => {
    const value = configOverrides && Object.hasOwn(configOverrides, FIGURE_SIZE_PATHS[side])
      ? configOverrides[FIGURE_SIZE_PATHS[side]]
      : config?.canvas?.[FIGURE_SIZE_FIELDS[side]];
    return isAuto(value) ? null : value;
  };
  return { width: leaf('width'), height: leaf('height') };
};

/**
 * The size a committed request asks for, each side as stored, as the export check reads it
 * (gbdraw/web_support/export_limit.py `_FigureSizePayload`); null for a request with both sides Auto.
 * @param {CanonicalRenderRequest} request
 * @returns {Record<FigureSide, unknown> | null}
 */
export const requestedFigureSize = (request) => {
  const leaves = requestFigureSizeLeaves(request);
  return leaves.width === null && leaves.height === null ? null : leaves;
};

/**
 * The value a side's typed amount and unit write: empty is Auto; otherwise the amount as typed and the unit. An
 * amount Python refuses (`0`, `-5`, `1e3`) is sent as typed, and Generate shows Python's refusal (R7).
 * @param {string} text
 * @param {FigureLengthUnit} unit
 * @returns {FigureLengthFormValue}
 */
export const figureLengthInputValue = (text, unit) => {
  const amount = String(text).trim();
  return amount ? `${amount}${unit}` : null;
};

/**
 * The unit a side shows when its value comes from outside the field (Load, Undo, Reset, a drawing switch): the
 * value's unit; Auto and a value the form cannot split keep `current`.
 * @param {unknown} value
 * @param {FigureLengthUnit} current
 */
export const figureLengthUnitOf = (value, current) => splitFigureLength(value)?.unit ?? current;

/**
 * The amount a side's number field shows: empty for Auto.
 * @param {unknown} value
 */
export const figureLengthAmountText = (value) => splitFigureLength(value)?.amount ?? '';

/**
 * `px` in `unit`, rounded as Python's messages round it (`messageDecimals`), without trailing zeros: `180`,
 * `143.2`, `7.09`. Display only: never a request value, an export attribute, or error text.
 * @param {number} px
 * @param {FigureLengthUnit} unit
 */
const lengthText = (px, unit) => String(Number((px / UNIT_TABLE.pxPerUnit[unit]).toFixed(UNIT_TABLE.messageDecimals[unit])));

/**
 * The size of a shown Result, as the line under the field shows it: first in the unit its committed request asked
 * for (its width's unit, else its height's), else in px; then in px, or in mm when the first is px (the Owner,
 * 2026-11-04): `180 × 143.2 mm` and `680 × 541 px`, or `1400 × 1450 px` and `370.4 × 383.6 mm`.
 * @param {{ widthPx: number, heightPx: number, requested: Record<FigureSide, unknown> }} inputs The Result's root
 *   size in px and its committed request's leaves (`requestFigureSizeLeaves`).
 * @returns {{ first: string, second: string }}
 */
export const shownFigureSize = ({ widthPx, heightPx, requested }) => {
  const unit = splitFigureLength(requested.width)?.unit ?? splitFigureLength(requested.height)?.unit ?? 'px';
  /** @param {FigureLengthUnit} shownUnit */
  const size = (shownUnit) => `${lengthText(widthPx, shownUnit)} × ${lengthText(heightPx, shownUnit)} ${shownUnit}`;
  return { first: size(unit), second: size(unit === 'px' ? 'mm' : 'px') };
};

// A side counts as larger than requested from half a px on: the line shows px without decimals.
const LARGER_PX = 0.5;

/**
 * Web-only additions after Python's fit, which the fit never sees (Owner D-26: until they reach Python): Legend
 * edits (the replan grew the canvas past the box), Canvas padding, and a moved Legend, title, or scale.
 * @typedef {'legend-edits' | 'canvas-padding' | 'moved-items'} FigureGrowthCause
 * @typedef {{ sides: FigureSide[], causes: FigureGrowthCause[] }} FigureGrowth
 */

/**
 * Why a shown Result drawn to a requested size is larger than the canvas it asked for, or null when it is not.
 * @param {{ shown: { width: number, height: number }, base: { width: number, height: number },
 *   requested: import('./legend-layout.js').RequestedCanvas, padded: boolean, moved: boolean }} inputs `shown` is the
 *   root's px size; `base` the canvas before padding and moved items grew it (the recorded automatic canvas);
 *   `requested` the composition metadata's `requestedCanvas`.
 * @returns {FigureGrowth | null}
 */
export const figureGrowth = ({ shown, base, requested, padded, moved }) => {
  /** @type {FigureSide[]} */
  const sides = [];
  /** @type {Set<FigureGrowthCause>} */
  const causes = new Set();
  for (const [side, px] of /** @type {const} */ ([['width', requested.widthPx], ['height', requested.heightPx]])) {
    if (px === null || shown[side] < px + LARGER_PX) continue;
    sides.push(side);
    if (base[side] >= px + LARGER_PX) causes.add('legend-edits');
    if (shown[side] >= base[side] + LARGER_PX) {
      if (padded) causes.add('canvas-padding');
      if (moved) causes.add('moved-items');
    }
  }
  return sides.length === 0 ? null : { sides, causes: [...causes] };
};

/**
 * A requested side as the size line names it: `180 mm`; null for a value the form cannot split.
 * @param {unknown} value
 */
export const requestedLengthText = (value) => {
  const split = splitFigureLength(value);
  return split === null ? null : `${split.amount} ${split.unit}`;
};
