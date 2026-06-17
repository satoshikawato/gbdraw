import {
  buildFeatureElementIndex,
  getFeatureElements
} from '../feature-editor/svg-actions.js';

export const PREVIEW_FEATURE_SEARCH_MATCH_CLASS = 'gbdraw-preview-feature-search-match';
export const PREVIEW_FEATURE_SEARCH_ACTIVE_CLASS = 'gbdraw-preview-feature-search-active-match';
export const PREVIEW_FEATURE_SEARCH_DIMMED_CLASS = 'gbdraw-preview-feature-search-dimmed';

export const PREVIEW_FEATURE_SEARCH_CLASSES = Object.freeze([
  PREVIEW_FEATURE_SEARCH_MATCH_CLASS,
  PREVIEW_FEATURE_SEARCH_ACTIVE_CLASS,
  PREVIEW_FEATURE_SEARCH_DIMMED_CLASS
]);

export const resolvePreviewSvg = (root) => {
  if (!root) return null;
  if (root.matches?.('svg')) return root;
  return root.querySelector?.('svg') || null;
};

const setClassToken = (element, token, enabled) => {
  if (!element) return;
  if (element.classList?.toggle) {
    element.classList.toggle(token, Boolean(enabled));
    return;
  }
  const existing = String(element.getAttribute('class') || '').trim();
  const tokens = existing ? existing.split(/\s+/).filter(Boolean) : [];
  const nextTokens = tokens.filter((entry) => entry !== token);
  if (enabled) nextTokens.push(token);
  if (nextTokens.length) {
    element.setAttribute('class', nextTokens.join(' '));
  } else {
    element.removeAttribute('class');
  }
};

export const stripPreviewFeatureSearchClasses = (svg) => {
  if (!svg) return;
  PREVIEW_FEATURE_SEARCH_CLASSES.forEach((className) => {
    svg.querySelectorAll(`.${className}`).forEach((element) => {
      setClassToken(element, className, false);
    });
  });
};

export const getPreviewFeatureElementIndex = (svg) => buildFeatureElementIndex(svg, { rebuild: true });

export const applyPreviewFeatureSearchClasses = ({
  svg,
  matches = [],
  activeId = '',
  queryActive = false,
  featureIndex = null
} = {}) => {
  if (!svg) return new Map();
  const index = featureIndex || getPreviewFeatureElementIndex(svg);
  const matchedIds = new Set((Array.isArray(matches) ? matches : []).map((id) => String(id || '').trim()));
  const normalizedActiveId = String(activeId || '').trim();
  index.forEach((elements, svgId) => {
    const isMatch = queryActive && matchedIds.has(svgId);
    const isActive = isMatch && svgId === normalizedActiveId;
    elements.forEach((element) => {
      setClassToken(element, PREVIEW_FEATURE_SEARCH_MATCH_CLASS, isMatch);
      setClassToken(element, PREVIEW_FEATURE_SEARCH_ACTIVE_CLASS, isActive);
      setClassToken(element, PREVIEW_FEATURE_SEARCH_DIMMED_CLASS, queryActive && !isMatch);
    });
  });
  return index;
};

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

export const getFeatureScreenCenter = (svg, featureId, featureIndex = null) => {
  const rect = getFeatureRect(svg, featureId, featureIndex);
  if (!rect) return null;
  return {
    clientX: rect.left + (rect.width / 2),
    clientY: rect.top + (rect.height / 2)
  };
};

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
