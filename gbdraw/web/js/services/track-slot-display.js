// @ts-check
import { MODE_PROFILE_DATA } from '../mode-profiles.generated.js';
import { t } from '../utils/i18n.js';

/**
 * @param {unknown} value
 * @returns {string | null}
 */
export const normalizeOptionalText = (value) => {
  const text = String(value ?? '').trim();
  return text.length > 0 ? text : null;
};

/**
 * @param {unknown} value
 * @param {number} [digits]
 */
const roundDisplayNumber = (value, digits = 1) => {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return '';
  const fixed = numeric.toFixed(digits);
  return fixed.replace(/\.0+$/, '').replace(/(\.\d*?)0+$/, '$1');
};

/** @param {unknown} value */
export const isManualSlotValue = (value) => normalizeOptionalText(value) !== null;

/**
 * The packaged labels.length_threshold of gbdraw/data/config.toml, by diagram mode.
 * @type {Readonly<Record<'circular' | 'linear', number>>}
 */
const SIZE_CLASS_THRESHOLD_BP = MODE_PROFILE_DATA.autoSizes.sizeClassThresholdBp;

/**
 * What the size class reads from a drawing: the configuration values its
 * Session stored that the Web does not manage, by dotted config path.
 * @typedef {object} SizeClassDrawing
 * @property {Readonly<Record<string, unknown>>} unmanagedConfigOverrides
 */

/**
 * The short/long boundary the next Generate of a drawing uses: the
 * labels.length_threshold its Session stored, otherwise the packaged one.
 * @param {'circular' | 'linear'} mode
 * @param {SizeClassDrawing} drawing the drawing of that mode
 * @returns {number}
 */
const sizeClassThresholdBp = (mode, drawing) => {
  const stored = drawing.unmanagedConfigOverrides[`labels.length_threshold.${mode}`];
  return typeof stored === 'number' ? stored : SIZE_CLASS_THRESHOLD_BP[mode];
};

/**
 * The size class that picks the short or long Auto sizes of a drawing: that
 * of the longest record.
 * @param {'circular' | 'linear'} mode
 * @param {number[]} recordLengths
 * @param {SizeClassDrawing} drawing the drawing of that mode
 * @returns {'short' | 'long' | null} null when no length is known
 */
export const sizeClassForRecordLengths = (mode, recordLengths, drawing) => {
  if (recordLengths.length === 0) return null;
  return Math.max(...recordLengths) < sizeClassThresholdBp(mode, drawing) ? 'short' : 'long';
};

/**
 * The size class every diagram of a Generate takes, by the record lengths of
 * each diagram.
 * @param {'circular' | 'linear'} mode
 * @param {number[][]} diagramRecordLengths
 * @param {SizeClassDrawing} drawing the drawing of that mode
 * @returns {'short' | 'long' | null} null when no length is known or the diagrams differ
 */
export const sharedSizeClass = (mode, diagramRecordLengths, drawing) => {
  const sizeClasses = new Set(diagramRecordLengths.map((lengths) => sizeClassForRecordLengths(mode, lengths, drawing)));
  return sizeClasses.size === 1 ? [...sizeClasses][0] : null;
};

/**
 * Which records set the Auto Linear label font size, by label display scope
 * (linear_label_size_records in gbdraw/auto_sizes.py).
 * @type {Readonly<Record<string, string>>}
 */
const LINEAR_LABEL_SIZE_RECORDS = MODE_PROFILE_DATA.autoSizes.linearLabelSizeRecords;

/**
 * The lengths of the records that set the Auto Linear label font size.
 * @param {number[]} recordLengths the drawn records, in drawing order
 * @param {string} labelScope the label display scope of the drawing (Show Labels)
 * @returns {number[]}
 */
export const linearLabelSizeLengths = (recordLengths, labelScope) => (
  LINEAR_LABEL_SIZE_RECORDS[labelScope] === 'first' ? recordLengths.slice(0, 1) : recordLengths
);

/**
 * The note an Auto field shows: a value measured in the last render, or,
 * before one, an estimate marked as such (TK-15). `variesByRecord` says the
 * measured value differs between records (a Linear features row).
 * @param {string} text
 * @param {string} unit
 * @param {boolean} estimate
 * @param {boolean} [variesByRecord]
 */
const autoNote = (text, unit, estimate, variesByRecord = false) => {
  if (!text) return '';
  const params = { value: text, unit };
  if (estimate) return t('≈ {value} {unit} (estimate)', params);
  return variesByRecord ? t('{value} {unit} (auto; varies by record)', params) : t('{value} {unit} (auto)', params);
};

/**
 * The size note of an Overlay annotation row that no Generate has drawn yet:
 * the rows it follows, by slot id (`follows features → gc_content`), because
 * its size comes from those rows' resolved bands (OV-1300, LP-D99). Empty for
 * any other row. Circular Width/Radius and Linear height read it.
 * @param {{ renderer?: unknown, side?: unknown, params?: unknown } | null | undefined} slot
 * @returns {string}
 */
export const overlayFollowsNote = (slot) => {
  if (!slot || slot.renderer !== 'annotations' || slot.side !== 'overlay') return '';
  const params = slot.params && typeof slot.params === 'object'
    ? /** @type {{ anchor_slot?: unknown, to_slot?: unknown }} */ (slot.params)
    : {};
  const anchor = normalizeOptionalText(params.anchor_slot);
  if (!anchor) return '';
  const target = normalizeOptionalText(params.to_slot);
  return target ? t('follows {anchor} → {target}', { anchor, target }) : t('follows {anchor}', { anchor });
};

/** @param {unknown} value @param {{ estimate?: boolean, variesByRecord?: boolean }} [options] */
export const formatPxAuto = (value, { estimate = false, variesByRecord = false } = {}) => (
  autoNote(roundDisplayNumber(value, 1), 'px', estimate, variesByRecord)
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
/** @type {Readonly<Partial<Record<string, 'inside' | 'outside'>>>} */
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
 * One rendered slot of the Python track geometry.
 * @typedef {object} TrackSlotGeometry
 * @property {unknown} [slotId]
 * @property {unknown} [radiusFactor]
 * @property {unknown} [widthPx]
 * @property {unknown} [side]
 * @property {unknown} [innerGapPx]
 * @property {unknown} [outerGapPx]
 * @property {unknown} [heightPx]
 * @property {unknown} [spacingAfterPx]
 */

/**
 * One rendered record of the Python track geometry.
 * @typedef {object} TrackSlotGeometryRecord
 * @property {unknown} [resultIndex]
 * @property {unknown} [recordIndex]
 * @property {unknown} [axisRadiusPx]
 * @property {TrackSlotGeometry[]} [slots]
 */

/**
 * @typedef {object} FindTrackSlotGeometryOptions
 * @property {{ records?: TrackSlotGeometryRecord[] } | null} [geometry] The Python track geometry.
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
