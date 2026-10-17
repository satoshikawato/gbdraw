// @ts-check
import { MODE_PROFILE_DATA } from './mode-profiles.generated.js';
import { diagnosticError } from './utils/error-normalization.js';
import { DECIMAL_NUMBER_PATTERN } from './utils/optional-positive-number.js';

/**
 * The shape of the generated `MODE_PROFILE_DATA` (tools/generate_mode_profiles.py)
 * that this module reads.
 * @typedef {'circular' | 'linear'} ModeProfileName
 * @typedef {{ integer: boolean, maximum: number | null, minimum: number, reason: string }} ComparisonDomain
 * @typedef {{
 *   comparison: { alignmentLength: number, bitscore: number, evalue: number, identity: number },
 *   linearAxisColor: string | null,
 *   linearRulerAxisColor: string | null,
 *   tracks: { gc: boolean, skew: boolean }
 * }} ModeProfile
 * @typedef {{
 *   comparisonDomains: Record<string, ComparisonDomain>,
 *   defaultFontFamily: string,
 *   featureTypes: string[],
 *   modes: Record<ModeProfileName, ModeProfile>,
 *   version: number
 * }} ModeProfileData
 * @typedef {{ bitscore: number, evalue: string, identity: number, alignment_length: number }} ComparisonFilters
 * @typedef {{ min_bitscore: number, evalue: string, identity: number, alignment_length: number }} ComparisonState
 */

/** @type {readonly string[]} */
const MODE_NAMES = Object.freeze(['circular', 'linear']);
/**
 * @param {unknown} mode
 * @returns {ModeProfileName}
 */
const normalizeMode = (mode) => {
  const normalized = String(mode || '').trim().toLowerCase();
  if (!MODE_NAMES.includes(normalized)) {
    throw new TypeError(`Unsupported diagram mode: ${String(mode)}`);
  }
  return /** @type {ModeProfileName} */ (normalized); // checked against MODE_NAMES above
};

/** @param {unknown} value */
const formatEvalue = (value) => {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return String(value);
  if (numeric === 0) return '0';
  return numeric
    .toExponential()
    .replace(/\.?0+e/, 'e')
    .replace(/e\+?(-?)0+(\d+)/, 'e$1$2');
};

/**
 * @param {unknown} left
 * @param {unknown} right
 */
const valuesEquivalent = (left, right) => {
  const leftBlank = left === null || left === undefined ||
    (typeof left === 'string' && left.trim() === '');
  const rightBlank = right === null || right === undefined ||
    (typeof right === 'string' && right.trim() === '');
  if (leftBlank || rightBlank) return leftBlank && rightBlank;

  const leftNumber = Number(left);
  const rightNumber = Number(right);
  if (Number.isFinite(leftNumber) && Number.isFinite(rightNumber)) {
    const tolerance = Number.EPSILON *
      Math.max(1, Math.abs(leftNumber), Math.abs(rightNumber)) * 8;
    return Math.abs(leftNumber - rightNumber) <= tolerance;
  }
  return String(left ?? '').trim() === String(right ?? '').trim();
};

export const MODE_PROFILE_VERSION = /** @type {ModeProfileData} */ (MODE_PROFILE_DATA).version;
export const MODE_DEFAULT_FEATURE_TYPES = Object.freeze([.../** @type {ModeProfileData} */ (MODE_PROFILE_DATA).featureTypes]);
/** The packaged text font list (config.toml `objects.text.font_family`). */
export const PACKAGED_FONT_FAMILY = /** @type {ModeProfileData} */ (MODE_PROFILE_DATA).defaultFontFamily;

/** @param {unknown} mode */
export const modeProfile = (mode) => /** @type {ModeProfileData} */ (MODE_PROFILE_DATA).modes[normalizeMode(mode)];

/** @param {unknown} mode */
export const trackDefaultsForMode = (mode) => {
  const tracks = modeProfile(mode).tracks;
  return {
    gc: Boolean(tracks.gc),
    skew: Boolean(tracks.skew)
  };
};

/**
 * @param {unknown} mode
 * @returns {ComparisonFilters}
 */
export const comparisonFiltersForMode = (mode) => {
  const comparison = modeProfile(mode).comparison;
  return {
    bitscore: comparison.bitscore,
    evalue: formatEvalue(comparison.evalue),
    identity: comparison.identity,
    alignment_length: comparison.alignmentLength
  };
};

/**
 * @param {unknown} mode
 * @returns {ComparisonState}
 */
export const comparisonStateForMode = (mode) => {
  const filters = comparisonFiltersForMode(mode);
  return {
    min_bitscore: filters.bitscore,
    evalue: filters.evalue,
    identity: filters.identity,
    alignment_length: filters.alignment_length
  };
};

// [draft key, generated domain key, diagnostic field]
/** @type {ReadonlyArray<readonly [string, string, keyof ComparisonFilters]>} */
const COMPARISON_THRESHOLDS = Object.freeze([
  ['min_bitscore', 'bitscore', 'bitscore'],
  ['evalue', 'evalue', 'evalue'],
  ['identity', 'identity', 'identity'],
  ['alignment_length', 'alignmentLength', 'alignment_length']
]);

/**
 * The comparison thresholds of one Generate: a blank field takes the mode
 * default, and every value is evaluated on the Python-owned domains in
 * mode-profiles.generated.js. A violation is a typed INPUT_INVALID with the
 * same field and reason Python reports; the draft is never rewritten. The
 * e-value keeps its trimmed text so derived cache keys keep their shape.
 * @param {Record<string, unknown> | null | undefined} adv
 * @param {unknown} mode
 * @returns {Record<string, string | number>} by generated domain key
 */
export const resolveComparisonThresholds = (adv, mode) => {
  const defaults = comparisonFiltersForMode(mode);
  /** @type {Record<string, string | number>} */
  const resolved = {};
  for (const [draftKey, domainKey, field] of COMPARISON_THRESHOLDS) {
    const domain = /** @type {ModeProfileData} */ (MODE_PROFILE_DATA).comparisonDomains[domainKey];
    const raw = adv?.[draftKey];
    const blank = raw === null || raw === undefined || (typeof raw === 'string' && raw.trim() === '');
    const value = blank ? defaults[field] : raw;
    const text = typeof value === 'string' ? value.trim() : '';
    const numeric = typeof value === 'number' ? value
      : DECIMAL_NUMBER_PATTERN.test(text) ? Number(text) : NaN;
    if (
      !Number.isFinite(numeric)
      || numeric < domain.minimum
      || (domain.maximum !== null && numeric > domain.maximum)
      || (domain.integer && !Number.isInteger(numeric))
    ) {
      throw diagnosticError('INPUT_INVALID', { field, reason: domain.reason });
    }
    resolved[domainKey] = field === 'evalue' ? (text || String(numeric)) : numeric;
  }
  return resolved;
};

/**
 * @param {unknown} mode
 * @param {string} field
 */
export const comparisonProfileDefault = (mode, field) => {
  const defaults = comparisonStateForMode(mode);
  if (!Object.prototype.hasOwnProperty.call(defaults, field)) {
    throw new TypeError(`Unsupported comparison field: ${String(field)}`);
  }
  return defaults[/** @type {keyof ComparisonState} */ (field)]; // an own key, checked above
};

/** @param {unknown} mode */
export const managedAdvStateForMode = (mode) => ({
  ...comparisonStateForMode(mode),
  axis_stroke_color: modeProfile(mode).linearAxisColor,
  pairwise_match_style: normalizeMode(mode) === 'linear' ? 'curve' : 'ribbon',
  plot_title_font_size: null,
  def_font_size: null
});

/** @param {unknown} mode */
const managedStateForMode = (mode) => ({ ...managedAdvStateForMode(mode), plot_title: '' });

// Whether Ruler on Axis can draw; otherwise the flag is ignored, as in
// `linear_axis_ruler_applies` (gbdraw/config/models/canvas.py).
/** @param {{ show_scale?: boolean, scale_style: string, linear_track_layout: string }} form */
export const linearAxisRulerApplies = (form) => (
  form.show_scale !== false &&
  form.scale_style === 'ruler' &&
  ['above', 'below'].includes(form.linear_track_layout)
);

/**
 * @param {{ axisColor?: unknown, rulerOnAxis?: boolean, managed?: boolean }} [options]
 * @returns {string | null}
 */
export const effectiveLinearAxisColor = ({
  axisColor = null,
  rulerOnAxis = false,
  managed = false
} = {}) => {
  const profile = modeProfile('linear');
  const normalized = axisColor === null || axisColor === undefined
    ? ''
    : String(axisColor).trim();
  if (!managed && normalized) return normalized;
  return (
    rulerOnAxis
      ? profile.linearRulerAxisColor
      : profile.linearAxisColor
  ) || normalized || null;
};

// Whether a drawing's profile field still holds its mode's default: the
// Linear axis color then follows the ruler (`effectiveLinearAxisColor`).
/**
 * @param {string} mode
 * @param {string} field
 * @param {unknown} value
 */
export const isModeProfileDefault = (mode, field, value) => {
  /** @type {Record<string, unknown>} */
  const defaults = managedStateForMode(mode);
  return Object.prototype.hasOwnProperty.call(defaults, field) && valuesEquivalent(value, defaults[field]);
};

export { MODE_PROFILE_DATA };
