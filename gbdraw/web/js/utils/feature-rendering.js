// @ts-check
export const FEATURE_RENDERING_VALUES = Object.freeze([
  'arrow',
  'rectangle',
  'underlay'
]);

export const DEFAULT_ARROW_HEAD_LENGTH_RATIO = 'auto';
export const DEFAULT_ARROW_SHAFT_WIDTH_RATIO = 1.0;

export const DEFAULT_FEATURE_RENDERINGS = Object.freeze({
  CDS: 'arrow',
  rRNA: 'arrow',
  tRNA: 'arrow',
  tmRNA: 'arrow',
  ncRNA: 'arrow',
  misc_RNA: 'arrow',
  repeat_region: 'underlay'
});

const FEATURE_RENDERING_VALUE_SET = new Set(FEATURE_RENDERING_VALUES);

/** @param {unknown} featureType */
export const defaultFeatureRendering = (featureType) => {
  // The key is any feature type, not only the defaulted ones: read own keys only.
  const key = String(featureType || '').trim();
  return Object.hasOwn(DEFAULT_FEATURE_RENDERINGS, key)
    ? /** @type {Readonly<Record<string, string>>} */ (DEFAULT_FEATURE_RENDERINGS)[key]
    : 'rectangle';
};

/** @param {unknown} value */
export const normalizeFeatureRendering = (value) => {
  const normalized = String(value ?? '').trim().toLowerCase();
  if (!FEATURE_RENDERING_VALUE_SET.has(normalized)) {
    throw new Error(
      `Unsupported feature rendering '${String(value)}'. Expected arrow, rectangle, or underlay.`
    );
  }
  return normalized;
};

/** @param {unknown} value Unset or Auto. */
const isAutoArrowHeadLengthRatio = (value) => value === null || value === undefined
  || (typeof value === 'string' && ['', 'auto'].includes(value.trim().toLowerCase()));
/** @param {unknown} value @returns {value is number} */
const isArrowHeadLengthRatio = (value) => typeof value === 'number' && Number.isFinite(value) && value > 0;
/** @param {unknown} value Unset. */
const isUnsetArrowShaftWidthRatio = (value) => value === null || value === undefined
  || (typeof value === 'string' && value.trim() === '');
/** @param {unknown} value @returns {value is number} */
const isArrowShaftWidthRatio = (value) => typeof value === 'number' && Number.isFinite(value) && value > 0 && value <= 1;

/** @param {unknown} value */
export const normalizeArrowHeadLengthRatio = (value) => {
  if (isAutoArrowHeadLengthRatio(value)) return DEFAULT_ARROW_HEAD_LENGTH_RATIO;
  if (!isArrowHeadLengthRatio(value)) {
    throw new Error('Arrow head length ratio must be Auto or a positive finite number.');
  }
  return value;
};

/** @param {unknown} value */
export const normalizeArrowShaftWidthRatio = (value) => {
  if (isUnsetArrowShaftWidthRatio(value)) return DEFAULT_ARROW_SHAFT_WIDTH_RATIO;
  if (!isArrowShaftWidthRatio(value)) {
    throw new Error('Arrow shaft width ratio must be a finite number greater than 0 and at most 1.');
  }
  return value;
};

/** @param {unknown} value */
export const arrowHeadLengthRatioForState = (value) => {
  const normalized = normalizeArrowHeadLengthRatio(value);
  return normalized === DEFAULT_ARROW_HEAD_LENGTH_RATIO ? null : normalized;
};

/** @param {{ arrow_head_length_ratio?: unknown, arrow_shaft_width_ratio?: unknown } | null} [adv] */
export const normalizedArrowGeometryState = (adv = {}) => ({
  arrow_head_length_ratio: arrowHeadLengthRatioForState(
    adv?.arrow_head_length_ratio
  ),
  arrow_shaft_width_ratio: normalizeArrowShaftWidthRatio(
    adv?.arrow_shaft_width_ratio
  )
});

/** A stored ratio saved as its text reads as its number. @param {unknown} value */
const persistedArrowGeometryValue = (value) => {
  if (typeof value !== 'string' || value.trim() === '') return value;
  const numeric = Number(value.trim());
  return Number.isFinite(numeric) ? numeric : value;
};

/** @param {{ arrow_head_length_ratio?: unknown, arrow_shaft_width_ratio?: unknown } | null} [adv] */
export const normalizedPersistedArrowGeometryState = (adv = {}) =>
  normalizedArrowGeometryState({
    arrow_head_length_ratio: persistedArrowGeometryValue(
      adv?.arrow_head_length_ratio
    ),
    arrow_shaft_width_ratio: persistedArrowGeometryValue(
      adv?.arrow_shaft_width_ratio
    )
  });

/**
 * A stored arrow head length ratio as the drawing holds it (null for Auto);
 * a value `normalizeArrowHeadLengthRatio` refuses is Auto, and lost.
 * @param {unknown} value
 * @returns {{ value: number | null, lost: boolean }}
 */
export const readPersistedArrowHeadLengthRatio = (value) => {
  const persisted = persistedArrowGeometryValue(value);
  if (isAutoArrowHeadLengthRatio(persisted)) return { value: null, lost: false };
  return isArrowHeadLengthRatio(persisted) ? { value: persisted, lost: false } : { value: null, lost: true };
};
/**
 * A stored arrow shaft width ratio; a value `normalizeArrowShaftWidthRatio`
 * refuses is the default, and lost.
 * @param {unknown} value
 * @returns {{ value: number, lost: boolean }}
 */
export const readPersistedArrowShaftWidthRatio = (value) => {
  const persisted = persistedArrowGeometryValue(value);
  if (isUnsetArrowShaftWidthRatio(persisted)) return { value: DEFAULT_ARROW_SHAFT_WIDTH_RATIO, lost: false };
  return isArrowShaftWidthRatio(persisted) ? { value: persisted, lost: false } : { value: DEFAULT_ARROW_SHAFT_WIDTH_RATIO, lost: true };
};

/** @param {unknown} source */
export const normalizeFeatureRenderingMap = (source) => {
  if (source === null || source === undefined) return {};
  if (typeof source !== 'object' || Array.isArray(source)) {
    throw new Error('Feature renderings must be an object.');
  }
  /** @type {Record<string, string>} */
  const normalized = {};
  Object.entries(source).forEach(([featureTypeRaw, rendering]) => {
    const featureType = String(featureTypeRaw || '').trim();
    if (!featureType) throw new Error('Feature rendering type must not be empty.');
    normalized[featureType] = normalizeFeatureRendering(rendering);
  });
  return normalized;
};

export const createDefaultFeatureRenderings = () => ({
  ...DEFAULT_FEATURE_RENDERINGS
});
