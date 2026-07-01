const { computed } = window.Vue;

export const FEATURE_SELECTION_ID_ATTRIBUTE = 'data-gbdraw-feature-id';
export const SELECTABLE_FEATURE_SELECTOR = `[${FEATURE_SELECTION_ID_ATTRIBUTE}]`;

export const FEATURE_SELECTED_CLASS = 'gbdraw-feature-selected';
export const FEATURE_ANCHOR_CLASS = 'gbdraw-feature-selection-anchor';
export const FEATURE_CANDIDATE_CLASS = 'gbdraw-feature-selection-candidate';

export const FEATURE_SELECTION_CLASSES = Object.freeze([
  FEATURE_SELECTED_CLASS,
  FEATURE_ANCHOR_CLASS,
  FEATURE_CANDIDATE_CLASS,
  'feature-selection-marquee',
  'feature-selection-status'
]);

const FEATURE_PART_SUFFIX_RE = /__part\d+$/;
const TWO_PI = Math.PI * 2;
const MARQUEE_THRESHOLD_PX = 4;

export const normalizeFeatureSelectionId = (value) =>
  String(value || '').trim().replace(FEATURE_PART_SUFFIX_RE, '');

export const getFeatureSelectionScope = (feature) => [
  feature?.fileIdx ?? '',
  feature?.record_idx ?? '',
  feature?.displayRecordId || feature?.record_id || ''
].join('::');

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

export const stripFeatureSelectionClasses = (svg) => {
  if (!svg) return;
  FEATURE_SELECTION_CLASSES.forEach((className) => {
    svg.querySelectorAll(`.${className}`).forEach((element) => {
      setClassToken(element, className, false);
    });
  });
};

const getElementSelectionId = (element) =>
  normalizeFeatureSelectionId(element?.getAttribute?.(FEATURE_SELECTION_ID_ATTRIBUTE));

const getSelectableTargetFromPoint = (eventLike, svg) => {
  if (!svg || !Number.isFinite(eventLike?.clientX) || !Number.isFinite(eventLike?.clientY)) return null;
  const stack = typeof document.elementsFromPoint === 'function'
    ? document.elementsFromPoint(eventLike.clientX, eventLike.clientY)
    : [];
  for (const element of stack) {
    if (!element || !svg.contains(element)) continue;
    const target = element.matches?.(SELECTABLE_FEATURE_SELECTOR)
      ? element
      : element.closest?.(SELECTABLE_FEATURE_SELECTOR);
    if (target && svg.contains(target)) return target;
  }
  return null;
};

export const getSelectableFeatureTarget = (eventLike, svg) => {
  const pointTarget = getSelectableTargetFromPoint(eventLike, svg);
  if (pointTarget) return pointTarget;
  const target = eventLike?.target?.closest?.(SELECTABLE_FEATURE_SELECTOR) || null;
  return target && svg?.contains?.(target) ? target : null;
};

export const applyFeatureSelectionClasses = (svg, state) => {
  if (!svg || !state) return;
  const selected = state.selectedFeatureIds?.value || new Set();
  const anchorId = normalizeFeatureSelectionId(state.selectedFeatureAnchorId?.value);
  stripFeatureSelectionClasses(svg);
  svg.querySelectorAll(SELECTABLE_FEATURE_SELECTOR).forEach((element) => {
    const id = getElementSelectionId(element);
    if (!id) return;
    const isSelected = selected.has(id);
    setClassToken(element, FEATURE_SELECTED_CLASS, isSelected);
    setClassToken(element, FEATURE_ANCHOR_CLASS, Boolean(isSelected && anchorId && id === anchorId));
  });
};

const normalizePointerRect = (drag) => {
  const left = Math.min(drag.startX, drag.currentX);
  const top = Math.min(drag.startY, drag.currentY);
  const right = Math.max(drag.startX, drag.currentX);
  const bottom = Math.max(drag.startY, drag.currentY);
  return {
    left,
    top,
    right,
    bottom,
    width: right - left,
    height: bottom - top
  };
};

const rectsIntersect = (left, right) => (
  left.left <= right.right &&
  left.right >= right.left &&
  left.top <= right.bottom &&
  left.bottom >= right.top
);

const getElementScreenRect = (element) => {
  if (!element || element.getAttribute?.('display') === 'none') return null;
  const rect = typeof element.getBoundingClientRect === 'function'
    ? element.getBoundingClientRect()
    : null;
  if (!rect || rect.right <= rect.left || rect.bottom <= rect.top) return null;
  return {
    left: rect.left,
    top: rect.top,
    right: rect.right,
    bottom: rect.bottom,
    width: rect.right - rect.left,
    height: rect.bottom - rect.top
  };
};

const unionRects = (rects) => {
  const validRects = rects.filter(Boolean);
  if (validRects.length === 0) return null;
  const left = Math.min(...validRects.map((rect) => rect.left));
  const top = Math.min(...validRects.map((rect) => rect.top));
  const right = Math.max(...validRects.map((rect) => rect.right));
  const bottom = Math.max(...validRects.map((rect) => rect.bottom));
  if (right <= left || bottom <= top) return null;
  return {
    left,
    top,
    right,
    bottom,
    width: right - left,
    height: bottom - top,
    centerX: left + ((right - left) / 2),
    centerY: top + ((bottom - top) / 2)
  };
};

const normalizeAngle = (angle) => {
  const normalized = angle % TWO_PI;
  return normalized < 0 ? normalized + TWO_PI : normalized;
};

const angleOnClockwiseInterval = (startAngle, endAngle, angle) => {
  const span = normalizeAngle(endAngle - startAngle);
  const delta = normalizeAngle(angle - startAngle);
  return delta <= span + 1e-9;
};

export const createFeatureSelection = ({ state, onMounted = null, onUnmounted = null } = {}) => {
  const {
    selectedFeatureIds,
    selectedFeatureAnchorId,
    featureSelectionStatus,
    featureSelectionSuppressNextClick,
    featureSelectionDrag,
    extractedFeatures,
    featuresBySvgId,
    selectedFeatures,
    svgContainer,
    clickedFeature,
    clickedPairwiseMatch,
    clickedLabel,
    mode
  } = state;

  let statusTimeoutId = null;
  let pendingDrag = null;

  const getCurrentSvg = () => svgContainer.value?.querySelector?.('svg') || null;

  const setStatus = (message, { timeoutMs = 2200 } = {}) => {
    if (statusTimeoutId !== null) {
      window.clearTimeout(statusTimeoutId);
      statusTimeoutId = null;
    }
    featureSelectionStatus.value = String(message || '');
    if (featureSelectionStatus.value && timeoutMs > 0) {
      statusTimeoutId = window.setTimeout(() => {
        featureSelectionStatus.value = '';
        statusTimeoutId = null;
      }, timeoutMs);
    }
  };

  const closePreviewPopups = () => {
    if (clickedFeature?.value) clickedFeature.value = null;
    if (clickedPairwiseMatch?.value) clickedPairwiseMatch.value = null;
    if (clickedLabel?.value) clickedLabel.value = null;
  };

  const syncFeatureSelectionClasses = (svg = getCurrentSvg()) => {
    applyFeatureSelectionClasses(svg, state);
  };

  const replaceSelectedFeatureIds = (ids, { anchorId = selectedFeatureAnchorId.value, status = '' } = {}) => {
    const nextIds = new Set(
      Array.from(ids || [])
        .map(normalizeFeatureSelectionId)
        .filter(Boolean)
    );
    selectedFeatureIds.value = nextIds;
    selectedFeatureAnchorId.value = normalizeFeatureSelectionId(anchorId);
    if (status) setStatus(status);
    syncFeatureSelectionClasses();
  };

  const clearFeatureSelection = ({ keepAnchorId = '', clearStatus = false } = {}) => {
    selectedFeatureIds.value = new Set();
    selectedFeatureAnchorId.value = normalizeFeatureSelectionId(keepAnchorId);
    featureSelectionDrag.active = false;
    featureSelectionDrag.committed = false;
    pendingDrag = null;
    if (clearStatus) setStatus('', { timeoutMs: 0 });
    syncFeatureSelectionClasses();
  };

  const markPlainFeatureClick = (featureId) => {
    clearFeatureSelection({ keepAnchorId: featureId, clearStatus: true });
  };

  const getFeatureById = (featureId) => {
    const id = normalizeFeatureSelectionId(featureId);
    if (!id) return null;
    const bySvgId = featuresBySvgId?.value instanceof Map ? featuresBySvgId.value : new Map();
    return bySvgId.get(id) ||
      (Array.isArray(extractedFeatures.value)
        ? extractedFeatures.value.find((feature) => normalizeFeatureSelectionId(feature?.svg_id) === id)
        : null) ||
      null;
  };

  const getSelectableElementsById = (svg) => {
    const indexed = new Map();
    if (!svg) return indexed;
    svg.querySelectorAll(SELECTABLE_FEATURE_SELECTOR).forEach((element) => {
      const id = getElementSelectionId(element);
      if (!id) return;
      if (!indexed.has(id)) indexed.set(id, []);
      indexed.get(id).push(element);
    });
    return indexed;
  };

  const buildRenderedFeatureGeometry = (svg, scope = '') => {
    const elementIndex = getSelectableElementsById(svg);
    const geometries = [];
    (Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []).forEach((feature) => {
      const id = normalizeFeatureSelectionId(feature?.svg_id);
      if (!id) return;
      if (scope && getFeatureSelectionScope(feature) !== scope) return;
      const elements = elementIndex.get(id) || [];
      if (elements.length === 0) return;
      const rect = unionRects(elements.map(getElementScreenRect));
      if (!rect) return;
      geometries.push({
        id,
        feature,
        scope: getFeatureSelectionScope(feature),
        rect
      });
    });
    return geometries;
  };

  const toggleFeatureSelection = (featureId) => {
    const id = normalizeFeatureSelectionId(featureId);
    if (!id || !getFeatureById(id)) return false;
    const nextIds = new Set(selectedFeatureIds.value || []);
    if (nextIds.has(id)) {
      nextIds.delete(id);
      replaceSelectedFeatureIds(nextIds, {
        anchorId: selectedFeatureAnchorId.value === id ? '' : selectedFeatureAnchorId.value
      });
    } else {
      nextIds.add(id);
      replaceSelectedFeatureIds(nextIds, { anchorId: id });
    }
    closePreviewPopups();
    return true;
  };

  const selectFeatureRange = (featureId, { additive = false } = {}) => {
    const clickedId = normalizeFeatureSelectionId(featureId);
    const clicked = getFeatureById(clickedId);
    if (!clickedId || !clicked) return false;

    const anchorId = normalizeFeatureSelectionId(selectedFeatureAnchorId.value);
    const anchor = getFeatureById(anchorId);
    if (!anchorId || !anchor) {
      const nextIds = additive ? new Set(selectedFeatureIds.value || []) : new Set();
      nextIds.add(clickedId);
      replaceSelectedFeatureIds(nextIds, { anchorId: clickedId });
      closePreviewPopups();
      return true;
    }

    const anchorScope = getFeatureSelectionScope(anchor);
    const clickedScope = getFeatureSelectionScope(clicked);
    if (anchorScope !== clickedScope) {
      setStatus('Range selection is limited to one record.');
      closePreviewPopups();
      return false;
    }

    const svg = getCurrentSvg();
    const geometries = buildRenderedFeatureGeometry(svg, anchorScope);
    const anchorGeometry = geometries.find((entry) => entry.id === anchorId);
    const clickedGeometry = geometries.find((entry) => entry.id === clickedId);
    if (!anchorGeometry || !clickedGeometry) return false;

    let rangeIds = [];
    if (mode.value === 'circular') {
      const recordBounds = unionRects(geometries.map((entry) => entry.rect));
      if (!recordBounds) return false;
      const centerX = recordBounds.centerX;
      const centerY = recordBounds.centerY;
      const angleFor = (geometry) => normalizeAngle(
        Math.atan2(geometry.rect.centerY - centerY, geometry.rect.centerX - centerX)
      );
      const startAngle = angleFor(anchorGeometry);
      const endAngle = angleFor(clickedGeometry);
      rangeIds = geometries
        .filter((entry) => angleOnClockwiseInterval(startAngle, endAngle, angleFor(entry)))
        .map((entry) => entry.id);
    } else {
      const startX = anchorGeometry.rect.centerX;
      const endX = clickedGeometry.rect.centerX;
      const minX = Math.min(startX, endX);
      const maxX = Math.max(startX, endX);
      rangeIds = geometries
        .filter((entry) => entry.rect.centerX >= minX && entry.rect.centerX <= maxX)
        .sort((left, right) =>
          left.rect.centerX - right.rect.centerX ||
          left.rect.centerY - right.rect.centerY ||
          String(left.feature?.type || '').localeCompare(String(right.feature?.type || '')) ||
          left.id.localeCompare(right.id)
        )
        .map((entry) => entry.id);
    }

    const nextIds = additive ? new Set(selectedFeatureIds.value || []) : new Set();
    rangeIds.forEach((id) => nextIds.add(id));
    replaceSelectedFeatureIds(nextIds, { anchorId });
    closePreviewPopups();
    return true;
  };

  const selectIntersectingFeatures = (svg, marqueeRect, additive = false) => {
    const nextIds = additive ? new Set(selectedFeatureIds.value || []) : new Set();
    const geometries = buildRenderedFeatureGeometry(svg);
    geometries.forEach((entry) => {
      if (rectsIntersect(entry.rect, marqueeRect)) {
        nextIds.add(entry.id);
      }
    });
    const firstNewId = Array.from(nextIds).find(Boolean) || '';
    replaceSelectedFeatureIds(nextIds, {
      anchorId: firstNewId || selectedFeatureAnchorId.value,
      status: nextIds.size > 0 ? `${nextIds.size.toLocaleString()} features selected.` : 'No rendered features in selection.'
    });
    closePreviewPopups();
  };

  const startMarqueePointer = (event, svg) => {
    if (!event?.shiftKey || event.button !== 0 || !svg) return false;
    pendingDrag = {
      svg,
      pointerId: event.pointerId,
      startX: event.clientX,
      startY: event.clientY,
      additive: Boolean(event.ctrlKey || event.metaKey)
    };
    featureSelectionDrag.active = false;
    featureSelectionDrag.committed = false;
    featureSelectionDrag.startX = event.clientX;
    featureSelectionDrag.startY = event.clientY;
    featureSelectionDrag.currentX = event.clientX;
    featureSelectionDrag.currentY = event.clientY;
    featureSelectionDrag.additive = pendingDrag.additive;
    if (typeof svg.setPointerCapture === 'function' && Number.isFinite(event.pointerId)) {
      try {
        svg.setPointerCapture(event.pointerId);
      } catch (_err) {
        // Pointer capture is best-effort for cross-browser SVG support.
      }
    }
    return true;
  };

  const moveMarqueePointer = (event) => {
    if (!pendingDrag || event.pointerId !== pendingDrag.pointerId) return false;
    const dx = event.clientX - pendingDrag.startX;
    const dy = event.clientY - pendingDrag.startY;
    const moved = Math.hypot(dx, dy);
    if (!featureSelectionDrag.active && moved < MARQUEE_THRESHOLD_PX) return false;
    if (!featureSelectionDrag.active) {
      closePreviewPopups();
      featureSelectionDrag.active = true;
      featureSelectionDrag.committed = false;
      featureSelectionDrag.additive = pendingDrag.additive;
    }
    featureSelectionDrag.currentX = event.clientX;
    featureSelectionDrag.currentY = event.clientY;
    return true;
  };

  const commitMarqueePointer = (event) => {
    if (!pendingDrag || event.pointerId !== pendingDrag.pointerId) return false;
    const wasActive = featureSelectionDrag.active;
    featureSelectionDrag.currentX = Number.isFinite(event.clientX) ? event.clientX : featureSelectionDrag.currentX;
    featureSelectionDrag.currentY = Number.isFinite(event.clientY) ? event.clientY : featureSelectionDrag.currentY;
    if (typeof pendingDrag.svg?.releasePointerCapture === 'function' && Number.isFinite(event.pointerId)) {
      try {
        pendingDrag.svg.releasePointerCapture(event.pointerId);
      } catch (_err) {
        // Pointer capture is best-effort for cross-browser SVG support.
      }
    }
    if (wasActive) {
      selectIntersectingFeatures(
        pendingDrag.svg,
        normalizePointerRect(featureSelectionDrag),
        Boolean(pendingDrag.additive)
      );
      featureSelectionSuppressNextClick.value = true;
      featureSelectionDrag.committed = true;
    }
    featureSelectionDrag.active = false;
    pendingDrag = null;
    return wasActive;
  };

  const cancelMarqueePointer = () => {
    if (pendingDrag?.svg && typeof pendingDrag.svg.releasePointerCapture === 'function') {
      try {
        pendingDrag.svg.releasePointerCapture(pendingDrag.pointerId);
      } catch (_err) {
        // Pointer capture is best-effort for cross-browser SVG support.
      }
    }
    pendingDrag = null;
    featureSelectionDrag.active = false;
    featureSelectionDrag.committed = false;
  };

  const consumeSuppressNextClick = () => {
    if (!featureSelectionSuppressNextClick.value) return false;
    featureSelectionSuppressNextClick.value = false;
    return true;
  };

  const featureSelectionMarqueeStyle = computed(() => {
    if (!featureSelectionDrag.active) return {};
    const rect = normalizePointerRect(featureSelectionDrag);
    return {
      left: `${Math.round(rect.left)}px`,
      top: `${Math.round(rect.top)}px`,
      width: `${Math.round(rect.width)}px`,
      height: `${Math.round(rect.height)}px`
    };
  });

  const handleEscape = (event) => {
    if (event.key !== 'Escape') return;
    if (featureSelectionDrag.active || pendingDrag) {
      cancelMarqueePointer();
      event.preventDefault();
      return;
    }
    const target = event.target;
    if (target?.closest?.('.feature-popup, .pairwise-match-popup, .label-popup, input, textarea, select, [contenteditable="true"]')) {
      return;
    }
    if ((selectedFeatureIds.value || new Set()).size > 0) {
      clearFeatureSelection({ clearStatus: true });
      event.preventDefault();
    }
  };

  const handleWindowBlur = () => {
    if (featureSelectionDrag.active || pendingDrag) cancelMarqueePointer();
  };

  if (typeof onMounted === 'function' && typeof onUnmounted === 'function') {
    onMounted(() => {
      document.addEventListener('keydown', handleEscape);
      window.addEventListener('blur', handleWindowBlur);
    });
    onUnmounted(() => {
      document.removeEventListener('keydown', handleEscape);
      window.removeEventListener('blur', handleWindowBlur);
      if (statusTimeoutId !== null) window.clearTimeout(statusTimeoutId);
    });
  }

  return {
    clearFeatureSelection,
    consumeSuppressNextClick,
    featureSelectionMarqueeStyle,
    getFeatureById,
    getSelectableFeatureTarget,
    markPlainFeatureClick,
    replaceSelectedFeatureIds,
    selectFeatureRange,
    startMarqueePointer,
    moveMarqueePointer,
    commitMarqueePointer,
    cancelMarqueePointer,
    selectedFeatures,
    setStatus,
    syncFeatureSelectionClasses,
    toggleFeatureSelection
  };
};
