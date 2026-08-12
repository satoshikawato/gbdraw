export const FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id';
export const RENDERED_FEATURE_ID_ATTRIBUTE = 'data-gbdraw-rendered-feature-id';
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

const featureElementIndexCache = new WeakMap();

const FEATURE_ELEMENT_SUFFIX_RE = /__(?:part|line)\d+$/;
const CONNECTOR_ID_SUFFIX_RE = /__line\d+$/;

export const normalizeFeatureIdentity = (value) =>
  String(value || '').trim().replace(FEATURE_ELEMENT_SUFFIX_RE, '');

export const getFeatureIdentity = (element) =>
  normalizeFeatureIdentity(
    element?.getAttribute?.(RENDERED_FEATURE_ID_ATTRIBUTE) ||
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

export const buildFeatureElementIndex = (svg, { markCursor = false } = {}) => {
  const indexed = new Map();
  if (!svg) return indexed;
  Array.from(svg.querySelectorAll(FEATURE_SELECTOR)).forEach((element) => {
    const id = getFeatureIdentity(element);
    if (!id) return;
    if (!indexed.has(id)) indexed.set(id, []);
    indexed.get(id).push(element);
    if (markCursor && element?.style) element.style.cursor = 'pointer';
  });
  featureElementIndexCache.set(svg, indexed);
  return indexed;
};

export const getFeatureElementIndex = (svg, { rebuild = false, markCursor = false } = {}) => {
  if (!svg) return new Map();
  const indexed = rebuild || !featureElementIndexCache.has(svg)
    ? buildFeatureElementIndex(svg, { markCursor })
    : (featureElementIndexCache.get(svg) || new Map());
  if (markCursor) {
    indexed.forEach((elements) => elements.forEach((element) => {
      if (element?.style) element.style.cursor = 'pointer';
    }));
  }
  return indexed;
};

export const clearFeatureElementIndex = (svg) => {
  if (svg) featureElementIndexCache.delete(svg);
};

export const getFeatureElements = (svg, featureId, featureIndex = null) => {
  const normalizedId = String(featureId || '').trim();
  if (!svg || !normalizedId) return [];

  const indexed = featureIndex || featureElementIndexCache.get(svg);
  const indexedElements = indexed?.get?.(normalizedId);
  if (indexedElements?.length) return indexedElements;

  const byId = svg.getElementById?.(normalizedId) || svg.querySelector?.(`#${CSS.escape(normalizedId)}`);
  return byId ? [byId] : [];
};

export const getFeatureFillElements = (svg, featureId, featureIndex = null) =>
  filterFeatureFillTargets(getFeatureElements(svg, featureId, featureIndex));
