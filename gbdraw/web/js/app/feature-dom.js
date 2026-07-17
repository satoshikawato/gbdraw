export const FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id';
export const FEATURE_PART_ATTRIBUTE = 'data-gbdraw-feature-part';
export const FEATURE_PART_BLOCK = 'block';
export const FEATURE_PART_CONNECTOR = 'connector';

export const FEATURE_SELECTOR = [
  `path[${FEATURE_ID_ATTRIBUTE}]`,
  `polygon[${FEATURE_ID_ATTRIBUTE}]`,
  `rect[${FEATURE_ID_ATTRIBUTE}]`,
  'path[id^="f"]',
  'polygon[id^="f"]',
  'rect[id^="f"]'
].join(', ');

const FEATURE_ELEMENT_SUFFIX_RE = /__(?:part|line)\d+$/;
const CONNECTOR_ID_SUFFIX_RE = /__line\d+$/;

export const normalizeFeatureIdentity = (value) =>
  String(value || '').trim().replace(FEATURE_ELEMENT_SUFFIX_RE, '');

export const getFeatureIdentity = (element) =>
  normalizeFeatureIdentity(
    element?.getAttribute?.(FEATURE_ID_ATTRIBUTE) ||
    element?.getAttribute?.('id') ||
    element?.id ||
    ''
  );

export const getFeaturePart = (element) => {
  const explicitPart = String(element?.getAttribute?.(FEATURE_PART_ATTRIBUTE) || '').trim();
  if (explicitPart === FEATURE_PART_BLOCK || explicitPart === FEATURE_PART_CONNECTOR) {
    return explicitPart;
  }

  const elementId = String(element?.getAttribute?.('id') || element?.id || '').trim();
  return CONNECTOR_ID_SUFFIX_RE.test(elementId) ? FEATURE_PART_CONNECTOR : FEATURE_PART_BLOCK;
};

export const isFeatureFillTarget = (element) => getFeaturePart(element) === FEATURE_PART_BLOCK;

export const filterFeatureFillTargets = (elements) =>
  Array.from(elements || []).filter(isFeatureFillTarget);
