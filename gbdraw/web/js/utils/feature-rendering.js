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

export const defaultFeatureRendering = (featureType) => (
  DEFAULT_FEATURE_RENDERINGS[String(featureType || '').trim()] || 'rectangle'
);

export const normalizeFeatureRendering = (value) => {
  const normalized = String(value ?? '').trim().toLowerCase();
  if (!FEATURE_RENDERING_VALUE_SET.has(normalized)) {
    throw new Error(
      `Unsupported feature rendering '${String(value)}'. Expected arrow, rectangle, or underlay.`
    );
  }
  return normalized;
};

export const normalizeArrowHeadLengthRatio = (value) => {
  if (
    value === null ||
    value === undefined ||
    (
      typeof value === 'string' &&
      ['', 'auto'].includes(value.trim().toLowerCase())
    )
  ) {
    return DEFAULT_ARROW_HEAD_LENGTH_RATIO;
  }
  if (typeof value !== 'number' || !Number.isFinite(value) || value <= 0) {
    throw new Error('Arrow head length ratio must be Auto or a positive finite number.');
  }
  return value;
};

export const normalizeArrowShaftWidthRatio = (value) => {
  if (
    value === null ||
    value === undefined ||
    (typeof value === 'string' && value.trim() === '')
  ) {
    return DEFAULT_ARROW_SHAFT_WIDTH_RATIO;
  }
  if (
    typeof value !== 'number' ||
    !Number.isFinite(value) ||
    value <= 0 ||
    value > 1
  ) {
    throw new Error('Arrow shaft width ratio must be a finite number greater than 0 and at most 1.');
  }
  return value;
};

export const arrowHeadLengthRatioForState = (value) => {
  const normalized = normalizeArrowHeadLengthRatio(value);
  return normalized === DEFAULT_ARROW_HEAD_LENGTH_RATIO ? null : normalized;
};

export const normalizeFeatureRenderingMap = (source) => {
  if (source === null || source === undefined) return {};
  if (typeof source !== 'object' || Array.isArray(source)) {
    throw new Error('Feature renderings must be an object.');
  }
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
