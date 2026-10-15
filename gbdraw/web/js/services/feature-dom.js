// @ts-check
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

/** @typedef {Map<string, Element[]>} FeatureElementIndex Rendered elements by normalized feature identity. */

/** An SVG root, or any element that holds the rendered features; the id lookup is optional. */
/** @typedef {Element & { getElementById?: (id: string) => Element | null }} FeatureSvgRoot */

/** @type {WeakMap<Element, FeatureElementIndex>} */
const featureElementIndexCache = new WeakMap();

const FEATURE_ELEMENT_SUFFIX_RE = /__(?:part|line)\d+$/;
const CONNECTOR_ID_SUFFIX_RE = /__line\d+$/;

/** @param {unknown} value */
export const normalizeFeatureIdentity = (value) =>
  String(value || '').trim().replace(FEATURE_ELEMENT_SUFFIX_RE, '');

/** @param {Element | null | undefined} element */
export const getFeatureIdentity = (element) =>
  normalizeFeatureIdentity(
    element?.getAttribute?.(RENDERED_FEATURE_ID_ATTRIBUTE) ||
    element?.getAttribute?.(FEATURE_ID_ATTRIBUTE) ||
    element?.getAttribute?.('id') ||
    element?.id ||
    ''
  );

/** @param {Element | null | undefined} element */
export const getFeaturePart = (element) => {
  const explicitPart = String(element?.getAttribute?.(FEATURE_PART_ATTRIBUTE) || '').trim();
  if (explicitPart === FEATURE_PART_BLOCK || explicitPart === FEATURE_PART_CONNECTOR) {
    return explicitPart;
  }

  const elementId = String(element?.getAttribute?.('id') || element?.id || '').trim();
  return CONNECTOR_ID_SUFFIX_RE.test(elementId) ? FEATURE_PART_CONNECTOR : FEATURE_PART_BLOCK;
};

// An automatic feature underlay (gbdraw/annotations/feature_underlays.py): a
// feature-owned block the renderer draws without a stroke (`stroke="none"`,
// `stroke-width="0"`, gbdraw/render/drawers/*/annotations.py), not with the
// feature block stroke.
export const AUTO_FEATURE_UNDERLAY_ATTRIBUTE = 'data-gbdraw-auto-feature-underlay';
export const AUTO_FEATURE_UNDERLAY_STROKE = Object.freeze({ color: 'none', width: 0 });
/** @param {Element | null | undefined} element */
export const isAutoFeatureUnderlay = (element) => element?.getAttribute?.(AUTO_FEATURE_UNDERLAY_ATTRIBUTE) === 'true';

/** @param {Element | null | undefined} element */
const isFeatureFillTarget = (element) => (
  getFeaturePart(element) === FEATURE_PART_BLOCK
  && !String(element?.getAttribute?.('id') || element?.id || '').endsWith('__outline')
);

/** @param {Iterable<Element> | ArrayLike<Element> | null | undefined} elements */
export const filterFeatureFillTargets = (elements) =>
  Array.from(elements || []).filter(isFeatureFillTarget);

/**
 * @param {FeatureSvgRoot | null | undefined} svg
 * @returns {FeatureElementIndex}
 */
export const buildFeatureElementIndex = (svg) => {
  /** @type {FeatureElementIndex} */
  const indexed = new Map();
  if (!svg) return indexed;
  Array.from(svg.querySelectorAll(FEATURE_SELECTOR)).forEach((element) => {
    const id = getFeatureIdentity(element);
    if (!id) return;
    if (!indexed.has(id)) indexed.set(id, []);
    // `indexed` holds `id` after the line above.
    /** @type {Element[]} */ (indexed.get(id)).push(element);
  });
  featureElementIndexCache.set(svg, indexed);
  return indexed;
};

/**
 * @param {FeatureSvgRoot | null | undefined} svg
 * @param {{ rebuild?: boolean }} [options]
 * @returns {FeatureElementIndex}
 */
export const getFeatureElementIndex = (svg, { rebuild = false } = {}) => {
  if (!svg) return new Map();
  const indexed = rebuild || !featureElementIndexCache.has(svg)
    ? buildFeatureElementIndex(svg)
    : (featureElementIndexCache.get(svg) || new Map());
  return indexed;
};

/** @param {FeatureSvgRoot | null | undefined} svg */
export const clearFeatureElementIndex = (svg) => {
  if (svg) featureElementIndexCache.delete(svg);
};

/**
 * @param {FeatureSvgRoot | null | undefined} svg
 * @param {unknown} featureId
 * @param {FeatureElementIndex | null} [featureIndex]
 * @returns {Element[]}
 */
export const getFeatureElements = (svg, featureId, featureIndex = null) => {
  const normalizedId = String(featureId || '').trim();
  if (!svg || !normalizedId) return [];

  const indexed = featureIndex || getFeatureElementIndex(svg);
  const indexedElements = indexed?.get?.(normalizedId);
  if (indexedElements?.length) return indexedElements;

  const byId = svg.getElementById?.(normalizedId) || svg.querySelector?.(`#${CSS.escape(normalizedId)}`);
  return byId ? [byId] : [];
};

/**
 * @param {FeatureSvgRoot | null | undefined} svg
 * @param {unknown} featureId
 * @param {FeatureElementIndex | null} [featureIndex]
 */
export const getFeatureFillElements = (svg, featureId, featureIndex = null) =>
  filterFeatureFillTargets(getFeatureElements(svg, featureId, featureIndex));
