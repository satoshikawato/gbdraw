// @ts-check
// Load's reading of the saved color and stroke edits (a Legend row's
// `colorOverrides` and `strokeOverrides`, a feature's stroke overrides): keys
// trimmed, values as a Session holds them, and an entry without a usable value
// left out. Session Load (`normalizeEditorStateData`, services/config.js) and
// the Session 46 Legend reading (`legendEditsByRowKey`,
// services/session-drawings-upgrade.js) read through it. Python twins:
// `_legend_color_overrides`, `_stroke_override_map`, `_stroke_width` in
// gbdraw/session_io.py.
import { normalizeOptionalHexColor } from '../utils/color-utils.js';
import { isUnsetValue, noteReplacedLeaf } from '../utils/field-reading.js';
/** @import { UnvalidatedFields } from './session-authority.js' */
/** @import { FieldLeaf } from '../utils/field-reading.js' */

/**
 * A saved stroke override as it is read: unvalidated.
 * @typedef {NonNullable<UnvalidatedFields<'strokeColor' | 'strokeWidth' | 'originalStrokeColor'
 *   | 'originalStrokeWidth'>>} StrokeOverrideFields
 */
/**
 * What Load could not read of a stroke override map: the fields replaced in a
 * kept entry (`at` is `<key>.<field>`) and the number of entries dropped
 * whole. An entry or field that holds nothing is neither.
 * @typedef {{ leaves: FieldLeaf[], dropped: number }} StrokeOverrideLosses
 */
/**
 * What Load could not read of a Legend's edits: the caption of each color
 * dropped, and the stroke overrides' losses.
 * @typedef {{ colors: string[], strokes: StrokeOverrideLosses }} LegendEditsUnread
 */
/**
 * A stroke override as Load reads it.
 * @typedef {{ strokeColor?: string, strokeWidth?: number, originalStrokeColor?: string | null,
 *   originalStrokeWidth?: number | null }} StrokeOverride
 */

/**
 * Whether a stored original stroke color is the `none` the app saves for a
 * figure drawn without a stroke: Load reads it as no color, as before, and
 * lists nothing.
 * @param {unknown} value
 */
export const isNoStrokeColor = (value) => typeof value === 'string' && value.trim().toLowerCase() === 'none';

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isPlainObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

/** @param {unknown} value */
export const normalizeStrokeWidth = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric >= 0 ? numeric : null;
};

/**
 * @param {unknown} source
 * @param {string[] | null} [dropped] Written: the caption of each color that holds a value the app does not accept.
 */
export const normalizeLegendColorOverrides = (source, dropped = null) => {
  /** @type {Record<string, string>} */
  const normalized = {};
  if (!isPlainObject(source)) return normalized;
  Object.entries(source).forEach(([key, value]) => {
    const caption = String(key || '').trim();
    const color = normalizeOptionalHexColor(value);
    if (caption && !color && !isUnsetValue(value)) dropped?.push(caption);
    if (!caption || !color) return;
    normalized[caption] = color;
  });
  return normalized;
};

/**
 * @param {unknown} source
 * @param {{ requireOverride?: boolean, leaves?: FieldLeaf[] | null }} [options] `leaves`: written, each field
 *   replaced (`at` is the field).
 * @returns {StrokeOverride | null}
 */
const normalizeStrokeOverride = (source, { requireOverride = false, leaves = null } = {}) => {
  if (!isPlainObject(source)) return null;
  /** @type {StrokeOverride} */
  const normalized = {};
  const fields = /** @type {StrokeOverrideFields} */ (source);
  const strokeColor = normalizeOptionalHexColor(fields.strokeColor);
  const strokeWidth = normalizeStrokeWidth(fields.strokeWidth);
  if (strokeColor !== null) normalized.strokeColor = strokeColor;
  if (strokeWidth !== null) normalized.strokeWidth = strokeWidth;
  if (Object.hasOwn(source, 'originalStrokeColor')) {
    normalized.originalStrokeColor = normalizeOptionalHexColor(fields.originalStrokeColor);
  }
  if (Object.hasOwn(source, 'originalStrokeWidth')) {
    normalized.originalStrokeWidth = normalizeStrokeWidth(fields.originalStrokeWidth);
  }
  for (const field of /** @type {const} */ (['strokeColor', 'strokeWidth', 'originalStrokeColor', 'originalStrokeWidth'])) {
    const used = normalized[field] ?? undefined;
    // A color read as its hex value keeps its meaning (OV-160).
    if (field.endsWith('Color') && (used !== undefined || (field === 'originalStrokeColor' && isNoStrokeColor(fields[field])))) continue;
    noteReplacedLeaf(leaves, field, fields[field], used);
  }
  if (requireOverride && !Object.hasOwn(normalized, 'strokeColor') && !Object.hasOwn(normalized, 'strokeWidth')) {
    return null;
  }
  return Object.keys(normalized).length > 0 ? normalized : null;
};

/** @param {unknown} value Whether a stored entry holds anything (an empty object holds nothing). */
const holdsValue = (value) => (isPlainObject(value)
  ? Object.values(value).some((field) => !isUnsetValue(field)) : !isUnsetValue(value));

/**
 * @param {unknown} source
 * @param {{ requireOverride?: boolean, losses?: StrokeOverrideLosses | null }} [options] `losses`: written.
 */
export const normalizeStrokeOverrideMap = (source, { requireOverride = false, losses = null } = {}) => {
  /** @type {Record<string, StrokeOverride>} */
  const normalized = {};
  if (!isPlainObject(source)) return normalized;
  Object.entries(source).forEach(([key, value]) => {
    const normalizedKey = String(key || '').trim();
    if (!normalizedKey) return;
    /** @type {FieldLeaf[]} */
    const leaves = [];
    const override = normalizeStrokeOverride(value, { requireOverride, leaves });
    if (!override) {
      if (losses && holdsValue(value)) losses.dropped += 1;
      return;
    }
    losses?.leaves.push(...leaves.map((leaf) => ({ ...leaf, at: `${normalizedKey}.${leaf.at}` })));
    normalized[normalizedKey] = override;
  });
  return normalized;
};
