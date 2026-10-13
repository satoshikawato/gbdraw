// @ts-check
// Load's reading of the saved color and stroke edits (a Legend row's
// `colorOverrides` and `strokeOverrides`, a feature's stroke overrides): keys
// trimmed, values as a Session holds them, and an entry without a usable value
// left out. Session Load (`normalizeEditorStateData`, services/config.js) and
// the Session 46 Legend reading (`legendEditsByRowKey`,
// services/mode-scoped-migration.js) read through it. Python twins:
// `_legend_color_overrides`, `_stroke_override_map`, `_stroke_width` in
// gbdraw/session_io.py.
import { normalizeOptionalHexColor } from '../utils/color-utils.js';
/** @import { UnvalidatedFields } from './session-authority.js' */

/**
 * A saved stroke override as it is read: unvalidated.
 * @typedef {NonNullable<UnvalidatedFields<'strokeColor' | 'strokeWidth' | 'originalStrokeColor'
 *   | 'originalStrokeWidth'>>} StrokeOverrideFields
 */
/**
 * A stroke override as Load reads it.
 * @typedef {{ strokeColor?: string, strokeWidth?: number, originalStrokeColor?: string | null,
 *   originalStrokeWidth?: number | null }} StrokeOverride
 */

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isPlainObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

/** @param {unknown} value */
export const normalizeStrokeWidth = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric >= 0 ? numeric : null;
};

/** @param {unknown} source */
export const normalizeLegendColorOverrides = (source) => {
  /** @type {Record<string, string>} */
  const normalized = {};
  if (!isPlainObject(source)) return normalized;
  Object.entries(source).forEach(([key, value]) => {
    const caption = String(key || '').trim();
    const color = normalizeOptionalHexColor(value);
    if (!caption || !color) return;
    normalized[caption] = color;
  });
  return normalized;
};

/**
 * @param {unknown} source
 * @param {{ requireOverride?: boolean }} [options]
 * @returns {StrokeOverride | null}
 */
const normalizeStrokeOverride = (source, { requireOverride = false } = {}) => {
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
  if (requireOverride && !Object.hasOwn(normalized, 'strokeColor') && !Object.hasOwn(normalized, 'strokeWidth')) {
    return null;
  }
  return Object.keys(normalized).length > 0 ? normalized : null;
};

/**
 * @param {unknown} source
 * @param {{ requireOverride?: boolean }} [options]
 */
export const normalizeStrokeOverrideMap = (source, { requireOverride = false } = {}) => {
  /** @type {Record<string, StrokeOverride>} */
  const normalized = {};
  if (!isPlainObject(source)) return normalized;
  Object.entries(source).forEach(([key, value]) => {
    const normalizedKey = String(key || '').trim();
    if (!normalizedKey) return;
    const override = normalizeStrokeOverride(value, { requireOverride });
    if (!override) return;
    normalized[normalizedKey] = override;
  });
  return normalized;
};
