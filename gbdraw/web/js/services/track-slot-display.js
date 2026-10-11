// @ts-check
export const normalizeOptionalText = (value) => {
  const text = String(value ?? '').trim();
  return text.length > 0 ? text : null;
};

const roundDisplayNumber = (value, digits = 1) => {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return '';
  const fixed = numeric.toFixed(digits);
  return fixed.replace(/\.0+$/, '').replace(/(\.\d*?)0+$/, '$1');
};

export const isManualSlotValue = (value) => normalizeOptionalText(value) !== null;

/**
 * The note an Auto field shows: a value measured in the last render, or,
 * before one, an estimate marked as such (TK-15).
 * @param {string} text
 * @param {string} unit
 * @param {boolean} estimate
 */
const autoNote = (text, unit, estimate) => {
  if (!text) return '';
  return estimate ? `≈ ${text} ${unit} (estimate)` : `${text} ${unit} (auto)`;
};

/** @param {unknown} value @param {{ estimate?: boolean }} [options] */
export const formatPxAuto = (value, { estimate = false } = {}) => (
  autoNote(roundDisplayNumber(value, 1), 'px', estimate)
);

/**
 * `awayFromAxis` rounds a ticks row's anchor to two decimals away from the
 * Axis (down inside it, up outside it), so that typing the shown value back
 * does not move the ticks into the neighbouring row (GX-18).
 * @param {unknown} value
 * @param {{ estimate?: boolean, awayFromAxis?: boolean }} [options]
 */
export const formatRadiusFactorAuto = (value, { estimate = false, awayFromAxis = false } = {}) => {
  const numeric = Number(value);
  const rounded = awayFromAxis && Number.isFinite(numeric)
    ? (numeric < 1 ? Math.floor((numeric * 100) + 1e-9) : Math.ceil((numeric * 100) - 1e-9)) / 100
    : numeric;
  return autoNote(roundDisplayNumber(rounded, 2), 'R', estimate);
};

// tick_sides_for_tick_label_layout (gbdraw/tracks/circular.py): the side of
// the anchor the tick marks grow to. tick_only follows the row's side;
// label_only draws no marks.
const TICK_SIDE_BY_LAYOUT = Object.freeze({ label_out_tick_in: 'inside', label_in_tick_out: 'outside' });

/**
 * The anchor of a rendered ticks row, as a factor of the Axis radius. The
 * payload's radiusFactor is the centre of the tick band (its persisted
 * meaning); `r` pins the anchor the marks grow from (GX-18). An explicit
 * width is the tick length; an Auto row uses the default length
 * max(6 px, 0.025 R) (_default_tick_length_px, gbdraw/svg/circular_ticks.py).
 * @param {{ radiusFactor?: unknown, widthPx?: unknown, side?: unknown } | null | undefined} slotGeometry
 * @param {unknown} axisRadiusPx
 * @param {unknown} tickLabelLayout
 * @returns {number | null}
 */
export const tickAnchorRadiusFactor = (slotGeometry, axisRadiusPx, tickLabelLayout) => {
  const centre = Number(slotGeometry?.radiusFactor);
  const axis = Number(axisRadiusPx);
  if (!Number.isFinite(centre) || !Number.isFinite(axis) || axis <= 0) return null;
  const width = Number(slotGeometry?.widthPx);
  const lengthPx = Number.isFinite(width) && width > 0 ? width : Math.max(6, 0.025 * axis);
  const layout = String(tickLabelLayout || 'label_out_tick_in').trim().toLowerCase();
  const tickSide = layout === 'tick_only'
    ? (slotGeometry?.side === 'outside' ? 'outside' : 'inside')
    : TICK_SIDE_BY_LAYOUT[layout];
  const sign = tickSide === 'outside' ? 1 : (tickSide === 'inside' ? -1 : 0);
  return centre - ((sign * lengthPx) / (2 * axis));
};

/**
 * @typedef {object} FindTrackSlotGeometryOptions
 * @property {{ records?: any[] } | null} [geometry] The Python track geometry.
 * @property {number} [resultIndex]
 * @property {number} [recordIndex]
 * @property {string | null} [slotId]
 */

// Python emits geometry only for rendered rows, so geometry belongs to a row
// by slot ID alone; a row without rendered geometry has none.
/**
 * The geometry record of the shown Result and record, which also carries the
 * record's `axisRadiusPx`.
 * @param {FindTrackSlotGeometryOptions} [options]
 */
export const findTrackSlotGeometryRecord = ({ geometry, resultIndex = 0, recordIndex = 0 } = {}) => {
  if (!geometry || typeof geometry !== 'object' || !Array.isArray(geometry.records)) return null;
  const wantedResult = Number(resultIndex) || 0;
  const wantedRecord = Number(recordIndex) || 0;
  return geometry.records.find((record) => (
    Number(record?.resultIndex ?? 0) === wantedResult &&
    Number(record?.recordIndex ?? 0) === wantedRecord
  )) || geometry.records.find((record) => Number(record?.resultIndex ?? 0) === wantedResult) || null;
};

/** @param {FindTrackSlotGeometryOptions} [options] */
export const findTrackSlotGeometry = (options = {}) => {
  const wantedSlotId = normalizeOptionalText(options.slotId);
  if (!wantedSlotId) return null;
  const matchingRecord = findTrackSlotGeometryRecord(options);
  if (!matchingRecord || !Array.isArray(matchingRecord.slots)) return null;
  return matchingRecord.slots.find(
    (slot) => normalizeOptionalText(slot?.slotId) === wantedSlotId
  ) || null;
};
