// @ts-check
import {
  getFeatureElementIndex,
  getFeatureElements
} from '../feature-editor/svg-actions.js';
import {
  PREVIEW_FEATURE_SEARCH_ACTIVE_CLASS,
  PREVIEW_FEATURE_SEARCH_MATCH_CLASS,
  PREVIEW_FEATURE_SEARCH_ROOT_ACTIVE_CLASS,
  PREVIEW_FEATURE_SEARCH_ROOT_UPDATING_CLASS,
  setClassToken,
  stripPreviewFeatureSearchClasses
} from '../../services/svg-serialization.js';

export const resolvePreviewSvg = (root) => {
  if (!root) return null;
  if (root.matches?.('svg')) return root;
  return root.querySelector?.('svg') || null;
};

export const getPreviewFeatureElementIndex = (svg) => getFeatureElementIndex(svg);

/**
 * The DOM state that the search classes last applied to one preview SVG.
 * @typedef {object} PreviewFeatureSearchDomState
 * @property {SVGSVGElement | null} svg
 * @property {boolean} queryActive
 * @property {Set<string>} matchedIds
 * @property {string} activeId
 * @property {number} updateGeneration
 * @property {any} firstFrame handle of the first animation frame (or timeout)
 * @property {any} secondFrame handle of the second animation frame (or timeout)
 */

/**
 * Rendered feature elements by SVG id.
 * @typedef {Map<string, Element[]>} FeatureElementIndex
 */

/** @returns {PreviewFeatureSearchDomState} */
export const createPreviewFeatureSearchDomState = () => ({
  svg: null,
  queryActive: false,
  matchedIds: new Set(),
  activeId: '',
  updateGeneration: 0,
  firstFrame: null,
  secondFrame: null
});

const frameApi = () => {
  const request = typeof requestAnimationFrame === 'function'
    ? requestAnimationFrame
    : (callback) => setTimeout(callback, 0);
  const cancel = typeof cancelAnimationFrame === 'function'
    ? cancelAnimationFrame
    : clearTimeout;
  return { request, cancel };
};

const cancelPreviewFeatureSearchFrames = (appliedState) => {
  const { cancel } = frameApi();
  if (appliedState.firstFrame != null) cancel(appliedState.firstFrame);
  if (appliedState.secondFrame != null) cancel(appliedState.secondFrame);
  appliedState.firstFrame = null;
  appliedState.secondFrame = null;
};

const setFeatureClass = (featureIndex, svgId, className, enabled) => {
  (featureIndex.get(String(svgId || '')) || []).forEach((element) => {
    setClassToken(element, className, enabled);
  });
};

/**
 * @typedef {object} ApplyPreviewActiveSearchMatchOptions
 * @property {FeatureElementIndex} [featureIndex]
 * @property {PreviewFeatureSearchDomState} [appliedState]
 * @property {string} [activeId]
 */

/**
 * @param {ApplyPreviewActiveSearchMatchOptions} [options]
 */
export const applyPreviewActiveSearchMatch = ({
  featureIndex,
  appliedState,
  activeId = ''
} = {}) => {
  if (!featureIndex || !appliedState) return;
  const previousActiveId = appliedState.activeId;
  const nextActiveId = String(activeId || '').trim();
  if (previousActiveId === nextActiveId) return;
  if (previousActiveId) {
    setFeatureClass(featureIndex, previousActiveId, PREVIEW_FEATURE_SEARCH_ACTIVE_CLASS, false);
  }
  if (nextActiveId && appliedState.matchedIds.has(nextActiveId)) {
    setFeatureClass(featureIndex, nextActiveId, PREVIEW_FEATURE_SEARCH_ACTIVE_CLASS, true);
    appliedState.activeId = nextActiveId;
  } else {
    appliedState.activeId = '';
  }
};

/**
 * @typedef {object} ApplyPreviewFeatureSearchClassesOptions
 * @property {SVGSVGElement | null} [svg]
 * @property {string[]} [matches]
 * @property {string} [activeId]
 * @property {boolean} [queryActive]
 * @property {FeatureElementIndex | null} [featureIndex]
 * @property {PreviewFeatureSearchDomState} [appliedState]
 */

/**
 * @param {ApplyPreviewFeatureSearchClassesOptions} [options]
 */
export const applyPreviewFeatureSearchClasses = ({
  svg,
  matches = [],
  activeId = '',
  queryActive = false,
  featureIndex = null,
  appliedState = createPreviewFeatureSearchDomState()
} = {}) => {
  if (!svg) return appliedState;
  const index = featureIndex || getPreviewFeatureElementIndex(svg);
  if (appliedState.svg && appliedState.svg !== svg) {
    stripPreviewFeatureSearchClasses(appliedState.svg);
    appliedState.matchedIds = new Set();
    appliedState.activeId = '';
    appliedState.queryActive = false;
  }
  appliedState.svg = svg;
  const matchedIds = new Set(
    (queryActive && Array.isArray(matches) ? matches : []).map((id) => String(id || '').trim()).filter(Boolean)
  );
  const normalizedActiveId = String(activeId || '').trim();
  appliedState.matchedIds.forEach((svgId) => {
    if (!matchedIds.has(svgId)) setFeatureClass(index, svgId, PREVIEW_FEATURE_SEARCH_MATCH_CLASS, false);
  });
  matchedIds.forEach((svgId) => {
    if (!appliedState.matchedIds.has(svgId)) setFeatureClass(index, svgId, PREVIEW_FEATURE_SEARCH_MATCH_CLASS, true);
  });
  setClassToken(svg, PREVIEW_FEATURE_SEARCH_ROOT_ACTIVE_CLASS, Boolean(queryActive));
  appliedState.matchedIds = matchedIds;
  appliedState.queryActive = Boolean(queryActive);
  applyPreviewActiveSearchMatch({ featureIndex: index, appliedState, activeId: normalizedActiveId });
  return appliedState;
};

/**
 * @param {ApplyPreviewFeatureSearchClassesOptions} [options]
 */
export const schedulePreviewFeatureSearchClasses = (options = {}) => {
  const { svg, appliedState } = options;
  if (!svg || !appliedState) return appliedState;
  cancelPreviewFeatureSearchFrames(appliedState);
  const generation = ++appliedState.updateGeneration;
  const { request } = frameApi();
  setClassToken(svg, PREVIEW_FEATURE_SEARCH_ROOT_UPDATING_CLASS, true);
  appliedState.firstFrame = request(() => {
    if (generation !== appliedState.updateGeneration) return;
    appliedState.firstFrame = null;
    applyPreviewFeatureSearchClasses(options);
    appliedState.secondFrame = request(() => {
      if (generation !== appliedState.updateGeneration) return;
      appliedState.secondFrame = null;
      setClassToken(svg, PREVIEW_FEATURE_SEARCH_ROOT_UPDATING_CLASS, false);
    });
  });
  return appliedState;
};

export const disposePreviewFeatureSearchDomState = (appliedState) => {
  if (!appliedState) return;
  cancelPreviewFeatureSearchFrames(appliedState);
  appliedState.updateGeneration += 1;
  stripPreviewFeatureSearchClasses(appliedState.svg);
  appliedState.svg = null;
  appliedState.queryActive = false;
  appliedState.matchedIds = new Set();
  appliedState.activeId = '';
};

/** @param {FeatureElementIndex | null} [featureIndex] */
const getFeatureRect = (svg, featureId, featureIndex = null) => {
  const elements = getFeatureElements(svg, featureId, featureIndex);
  let left = Infinity;
  let top = Infinity;
  let right = -Infinity;
  let bottom = -Infinity;
  elements.forEach((element) => {
    const rect = typeof element.getBoundingClientRect === 'function'
      ? element.getBoundingClientRect()
      : null;
    if (!rect || rect.right <= rect.left || rect.bottom <= rect.top) return;
    left = Math.min(left, rect.left);
    top = Math.min(top, rect.top);
    right = Math.max(right, rect.right);
    bottom = Math.max(bottom, rect.bottom);
  });
  if (!Number.isFinite(left) || !Number.isFinite(top) || right <= left || bottom <= top) return null;
  return { left, top, right, bottom, width: right - left, height: bottom - top };
};

/** @param {FeatureElementIndex | null} [featureIndex] */
export const getFeatureScreenCenter = (svg, featureId, featureIndex = null) => {
  const rect = getFeatureRect(svg, featureId, featureIndex);
  if (!rect) return null;
  return {
    clientX: rect.left + (rect.width / 2),
    clientY: rect.top + (rect.height / 2)
  };
};

/**
 * @typedef {object} CenterPreviewFeatureOptions
 * @property {SVGSVGElement | null} [svg]
 * @property {string} [featureId]
 * @property {FeatureElementIndex | null} [featureIndex]
 * @property {HTMLElement | null} [canvasContainer]
 * @property {{ x: number, y: number } | null} [canvasPan]
 * @property {number} [maxDelta]
 */

/**
 * @param {CenterPreviewFeatureOptions} [options]
 */
export const centerPreviewFeature = ({
  svg,
  featureId,
  featureIndex = null,
  canvasContainer,
  canvasPan,
  maxDelta = 1600
} = {}) => {
  if (!svg || !canvasContainer || !canvasPan) return false;
  const featureRect = getFeatureRect(svg, featureId, featureIndex);
  const containerRect = typeof canvasContainer.getBoundingClientRect === 'function'
    ? canvasContainer.getBoundingClientRect()
    : null;
  if (!featureRect || !containerRect || containerRect.width <= 0 || containerRect.height <= 0) return false;

  const featureCenterX = featureRect.left + (featureRect.width / 2);
  const featureCenterY = featureRect.top + (featureRect.height / 2);
  const containerCenterX = containerRect.left + (containerRect.width / 2);
  const containerCenterY = containerRect.top + (containerRect.height / 2);
  const clampDelta = (value) => Math.max(-maxDelta, Math.min(maxDelta, value));

  canvasPan.x += clampDelta(containerCenterX - featureCenterX);
  canvasPan.y += clampDelta(containerCenterY - featureCenterY);
  return true;
};
