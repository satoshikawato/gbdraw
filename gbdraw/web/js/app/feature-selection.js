// @ts-check
import { setClassToken } from '../services/preview-affordances.js';
import { featureParts, onTheWay, partsMidpoint } from '../services/region-axis.js';
import { numberText } from '../utils/i18n.js';
import { capturePointer, holdTextSelection } from '../utils/pointer-gesture.js';

// createFeatureSelection reads `window.Vue` when it runs, not when the module loads: other modules import this
// module's pure helpers (feature-editor/svg-actions.js imports getSelectableFeatureTarget), and loading them must not
// need Vue.

export const FEATURE_SELECTION_ID_ATTRIBUTE = 'data-gbdraw-feature-id';
export const RENDERED_FEATURE_SELECTION_ID_ATTRIBUTE =
  'data-gbdraw-rendered-feature-id';
export const SELECTABLE_FEATURE_SELECTOR = [
  `[${RENDERED_FEATURE_SELECTION_ID_ATTRIBUTE}]`,
  `[${FEATURE_SELECTION_ID_ATTRIBUTE}]`
].join(', ');

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
const FEATURE_SELECTION_CLASS_SELECTOR = FEATURE_SELECTION_CLASSES.map((className) => `.${className}`).join(', ');

const FEATURE_PART_SUFFIX_RE = /__part\d+$/;
const TWO_PI = Math.PI * 2;
const MARQUEE_THRESHOLD_PX = 4;
const TOOLBAR_EDGE_MARGIN_PX = 8;

/**
 * @import { ExtractedFeature } from '../services/feature-metadata-extraction.js'
 * @import { RegionAxisRecord } from '../services/region-axis.js'
 * @import { WebState } from '../state.js'
 */

/**
 * A listed feature as the selection reads it: a row of `state.js`'s
 * `extractedFeatures`, typed by its owner.
 * @typedef {ExtractedFeature} SelectionFeature
 */

/** @typedef {{ left: number, top: number, right: number, bottom: number }} Edges */
/** @typedef {Edges & { width: number, height: number }} ScreenRect */
/** @typedef {ScreenRect & { centerX: number, centerY: number }} UnionRect */
/** @typedef {{ id: string, feature: SelectionFeature, scope: string, rect: UnionRect }} RenderedFeatureGeometry */

/**
 * The marquee `featureSelectionDrag` of `state.js`.
 * @typedef {{
 *   active: boolean, committed: boolean, startX: number, startY: number,
 *   currentX: number, currentY: number, additive: boolean
 * }} FeatureSelectionDrag
 */

/** @typedef {{ clientX: number, clientY: number }} PointerPosition */

/**
 * The pointer position and target of an event the selection reads.
 * @typedef {{ clientX?: number, clientY?: number, target?: EventTarget | null }} SelectionEventLike
 */

/**
 * The members of the `state.js` state that `applyFeatureSelectionClasses` reads.
 * @typedef {Pick<FeatureSelectionState, 'selectedFeatureIds' | 'selectedFeatureAnchorId'>} FeatureSelectionClassState
 */

/** @param {unknown} value */
export const normalizeFeatureSelectionId = (value) =>
  String(value || '').trim().replace(FEATURE_PART_SUFFIX_RE, '');

/** @param {SelectionFeature | null | undefined} feature */
export const getFeatureSelectionScope = (feature) => [
  feature?.fileIdx ?? '',
  feature?.record_idx ?? '',
  feature?.displayRecordId || feature?.record_id || ''
].join('::');

/** @param {Element | null | undefined} svg */
export const stripFeatureSelectionClasses = (svg) => {
  if (!svg) return;
  FEATURE_SELECTION_CLASSES.forEach((className) => {
    svg.querySelectorAll(`.${className}`).forEach((element) => {
      setClassToken(element, className, false);
    });
  });
};

/** @param {Element | null | undefined} element */
const getElementSelectionId = (element) =>
  normalizeFeatureSelectionId(
    element?.getAttribute?.(RENDERED_FEATURE_SELECTION_ID_ATTRIBUTE) ||
    element?.getAttribute?.(FEATURE_SELECTION_ID_ATTRIBUTE)
  );

/**
 * @param {SelectionEventLike | null | undefined} eventLike
 * @param {Element | null | undefined} svg
 * @returns {Element | null}
 */
const getSelectableTargetFromPoint = (eventLike, svg) => {
  if (!svg || !Number.isFinite(eventLike?.clientX) || !Number.isFinite(eventLike?.clientY)) return null;
  const stack = typeof document.elementsFromPoint === 'function'
    // The Number.isFinite checks above are true only for an eventLike with both numbers.
    ? document.elementsFromPoint(/** @type {PointerPosition} */ (eventLike).clientX, /** @type {PointerPosition} */ (eventLike).clientY)
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

/**
 * @param {SelectionEventLike | null | undefined} eventLike
 * @param {Element | null | undefined} svg
 * @returns {Element | null}
 */
export const getSelectableFeatureTarget = (eventLike, svg) => {
  const pointTarget = getSelectableTargetFromPoint(eventLike, svg);
  if (pointTarget) return pointTarget;
  // A pointer event on the preview is dispatched to an Element; `closest?.` covers a target without it.
  const target = /** @type {Element | null | undefined} */ (eventLike?.target)?.closest?.(SELECTABLE_FEATURE_SELECTOR) || null;
  return target && svg?.contains?.(target) ? target : null;
};

/**
 * @param {Element | null | undefined} svg
 * @param {FeatureSelectionClassState | null | undefined} state
 */
export const applyFeatureSelectionClasses = (svg, state) => {
  if (!svg || !state) return;
  const selected = state.selectedFeatureIds.value;
  const anchorId = normalizeFeatureSelectionId(state.selectedFeatureAnchorId.value);
  stripFeatureSelectionClasses(svg);
  svg.querySelectorAll(SELECTABLE_FEATURE_SELECTOR).forEach((element) => {
    const id = getElementSelectionId(element);
    if (!id) return;
    const isSelected = selected.has(id);
    setClassToken(element, FEATURE_SELECTED_CLASS, isSelected);
    setClassToken(element, FEATURE_ANCHOR_CLASS, Boolean(isSelected && anchorId && id === anchorId));
  });
};

/**
 * @param {Pick<FeatureSelectionDrag, 'startX' | 'startY' | 'currentX' | 'currentY'>} drag
 * @returns {ScreenRect}
 */
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

/**
 * @param {Edges} left
 * @param {Edges} right
 */
const rectsIntersect = (left, right) => (
  left.left <= right.right &&
  left.right >= right.left &&
  left.top <= right.bottom &&
  left.bottom >= right.top
);

/**
 * @param {number} value
 * @param {number | undefined} min
 * @param {number | undefined} max
 */
const clampNumber = (value, min, max) => {
  if (!Number.isFinite(value)) return 0;
  // Number.isFinite is true only for a number, so past these checks min and max are numbers.
  if (!Number.isFinite(min) || !Number.isFinite(max) || /** @type {number} */ (min) > /** @type {number} */ (max)) return value;
  return Math.min(/** @type {number} */ (max), Math.max(/** @type {number} */ (min), value));
};

/**
 * @param {Element | null | undefined} element
 * @returns {ScreenRect | null}
 */
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

/**
 * @param {(ScreenRect | null)[]} rects
 * @returns {UnionRect | null}
 */
const unionRects = (rects) => {
  // filter(Boolean) removes the nulls.
  const validRects = /** @type {ScreenRect[]} */ (rects.filter(Boolean));
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

/** @param {number} angle */
const normalizeAngle = (angle) => {
  const normalized = angle % TWO_PI;
  return normalized < 0 ? normalized + TWO_PI : normalized;
};

/**
 * @param {number} startAngle
 * @param {number} endAngle
 * @param {number} angle
 */
const angleOnClockwiseInterval = (startAngle, endAngle, angle) => {
  const span = normalizeAngle(endAngle - startAngle);
  const delta = normalizeAngle(angle - startAngle);
  return delta <= span + 1e-9;
};

/**
 * The members of the `state.js` state that the selection reads and writes; a
 * ref is `{ value }` because Vue comes from `window.Vue`.
 * @typedef {object} FeatureSelectionState
 * @property {{ value: Set<string> }} selectedFeatureIds
 * @property {{ value: string }} selectedFeatureAnchorId
 * @property {{ value: string }} featureSelectionStatus
 * @property {{ value: boolean }} featureSelectionSuppressNextClick
 * @property {FeatureSelectionDrag} featureSelectionDrag
 * @property {{ value: SelectionFeature[] }} extractedFeatures
 * @property {{ value: Map<string, SelectionFeature> }} featuresBySvgId
 * @property {{ value: HTMLElement | null }} svgContainer
 * @property {WebState['clickedFeature']} clickedFeature The open popup, closed by setting null.
 * @property {WebState['clickedPairwiseMatch']} clickedPairwiseMatch The open popup, closed by setting null.
 * @property {{ value: string }} mode
 */

/**
 * @typedef {object} FeatureSelectionOptions
 * @property {FeatureSelectionState} state
 * @property {((callback: () => void) => void) | null} [onMounted] Vue `onMounted`.
 * @property {((callback: () => void) => void) | null} [onUnmounted] Vue `onUnmounted`.
 * @property {(recordKey: string) => RegionAxisRecord | null} [axisRecordOf] The axis Python wrote for a record of
 *   the shown Result, else null (a Result saved before Python wrote it).
 * @property {() => boolean} [previewPopupsHeld] While it holds, a selection gesture leaves the feature
 *   and match popups as they are.
 */

/** @param {FeatureSelectionOptions} options */
export const createFeatureSelection = ({
  state, onMounted = null, onUnmounted = null, previewPopupsHeld = () => false, axisRecordOf = () => null
}) => {
  const { computed, reactive } = window.Vue;
  const {
    selectedFeatureIds,
    selectedFeatureAnchorId,
    featureSelectionStatus,
    featureSelectionSuppressNextClick,
    featureSelectionDrag,
    extractedFeatures,
    featuresBySvgId,
    svgContainer,
    clickedFeature,
    clickedPairwiseMatch,
    mode
  } = state;

  /** @type {number | null} */
  let statusTimeoutId = null;
  /** @type {{ svg: Element, pointerId: number, startX: number, startY: number, additive: boolean } | null} */
  let pendingDrag = null;
  /** Ends the marquee's hold on the pointer: its capture and the hold on text selection (OV-2478). */
  let releaseGesture = () => {};
  /** @type {{ startX: number, startY: number, originX: number, originY: number, minX?: number, maxX?: number, minY?: number, maxY?: number } | null} */
  let activeToolbarDrag = null;
  const toolbarOffset = reactive({ x: 0, y: 0 });

  const getCurrentSvg = () => svgContainer.value?.querySelector('svg') || null;

  /**
   * @param {unknown} message
   * @param {{ timeoutMs?: number }} [options]
   */
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
    if (previewPopupsHeld()) return;
    if (clickedFeature.value) clickedFeature.value = null;
    if (clickedPairwiseMatch.value) clickedPairwiseMatch.value = null;
  };

  const syncFeatureSelectionClasses = (svg = getCurrentSvg()) => {
    applyFeatureSelectionClasses(svg, state);
  };

  /** @param {HTMLElement | null} toolbar */
  const getToolbarDragBounds = (toolbar) => {
    if (!toolbar || typeof toolbar.getBoundingClientRect !== 'function') return null;
    const toolbarRect = toolbar.getBoundingClientRect();
    const boundsElement = toolbar.offsetParent || toolbar.parentElement || svgContainer.value?.parentElement || null;
    const boundsRect = boundsElement?.getBoundingClientRect?.() || {
      left: 0,
      top: 0,
      right: window.innerWidth,
      bottom: window.innerHeight
    };
    return {
      minX: toolbarOffset.x + boundsRect.left + TOOLBAR_EDGE_MARGIN_PX - toolbarRect.left,
      maxX: toolbarOffset.x + boundsRect.right - TOOLBAR_EDGE_MARGIN_PX - toolbarRect.right,
      minY: toolbarOffset.y + boundsRect.top + TOOLBAR_EDGE_MARGIN_PX - toolbarRect.top,
      maxY: toolbarOffset.y + boundsRect.bottom - TOOLBAR_EDGE_MARGIN_PX - toolbarRect.bottom
    };
  };

  const featureSelectionToolbarStyle = computed(() => ({
    transform: `translate(${Math.round(toolbarOffset.x)}px, ${Math.round(toolbarOffset.y)}px)`
  }));

  /** @param {PointerEvent} event */
  const moveToolbarDrag = (event) => {
    if (!activeToolbarDrag) return;
    const nextX = activeToolbarDrag.originX + ((Number(event.clientX) || 0) - activeToolbarDrag.startX);
    const nextY = activeToolbarDrag.originY + ((Number(event.clientY) || 0) - activeToolbarDrag.startY);
    toolbarOffset.x = clampNumber(nextX, activeToolbarDrag.minX, activeToolbarDrag.maxX);
    toolbarOffset.y = clampNumber(nextY, activeToolbarDrag.minY, activeToolbarDrag.maxY);
    event.preventDefault();
  };

  const stopToolbarDrag = () => {
    if (!activeToolbarDrag) return;
    activeToolbarDrag = null;
    document.removeEventListener('pointermove', moveToolbarDrag, true);
    document.removeEventListener('pointerup', stopToolbarDrag, true);
    document.removeEventListener('pointercancel', stopToolbarDrag, true);
    window.removeEventListener('blur', stopToolbarDrag);
  };

  /** @param {PointerEvent} event */
  const startToolbarDrag = (event) => {
    if (event.button != null && event.button !== 0) return;
    // The toolbar and its children are Elements; `closest?.` covers a target without it.
    const fromHandle = Boolean(/** @type {Element | null} */ (event.target)?.closest?.('.feature-selection-toolbar-drag-handle'));
    if (!fromHandle && /** @type {Element | null} */ (event.target)?.closest?.('input, select, button, label, textarea, a')) return;
    stopToolbarDrag();
    // The selector matches the toolbar, an HTML element.
    const toolbar = /** @type {HTMLElement | null} */ (/** @type {Element | null} */ (event.currentTarget)?.closest?.('.feature-selection-toolbar') || null);
    /** @type {Partial<NonNullable<ReturnType<typeof getToolbarDragBounds>>>} */
    const bounds = getToolbarDragBounds(toolbar) || {};
    activeToolbarDrag = {
      startX: Number(event.clientX) || 0,
      startY: Number(event.clientY) || 0,
      originX: Number(toolbarOffset.x) || 0,
      originY: Number(toolbarOffset.y) || 0,
      minX: bounds.minX,
      maxX: bounds.maxX,
      minY: bounds.minY,
      maxY: bounds.maxY
    };
    document.addEventListener('pointermove', moveToolbarDrag, true);
    document.addEventListener('pointerup', stopToolbarDrag, true);
    document.addEventListener('pointercancel', stopToolbarDrag, true);
    window.addEventListener('blur', stopToolbarDrag);
    event.preventDefault();
    event.stopPropagation();
  };

  /**
   * @param {Iterable<unknown> | null | undefined} ids
   * @param {{ anchorId?: unknown, status?: string }} [options]
   */
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

  /** @param {{ keepAnchorId?: unknown, clearStatus?: boolean, syncDom?: boolean }} [options] */
  const clearFeatureSelection = ({
    keepAnchorId = '',
    clearStatus = false,
    syncDom = true
  } = {}) => {
    // An empty selection stays the same Set, so a plain click or a background
    // click re-renders nothing, and the Result is restyled only when it still
    // shows a selection class (OV-913).
    const alreadyEmpty = !selectedFeatureIds.value.size;
    if (!alreadyEmpty) selectedFeatureIds.value = new Set();
    selectedFeatureAnchorId.value = normalizeFeatureSelectionId(keepAnchorId);
    featureSelectionDrag.active = false;
    featureSelectionDrag.committed = false;
    pendingDrag = null;
    if (clearStatus) setStatus('', { timeoutMs: 0 });
    if (!syncDom) return;
    const svg = getCurrentSvg();
    if (alreadyEmpty && !svg?.querySelector(FEATURE_SELECTION_CLASS_SELECTOR)) return;
    syncFeatureSelectionClasses(svg);
  };

  /** @param {unknown} featureId */
  const markPlainFeatureClick = (featureId) => {
    clearFeatureSelection({ keepAnchorId: featureId, clearStatus: true });
  };

  /**
   * @param {unknown} featureId
   * @returns {SelectionFeature | null}
   */
  const getFeatureById = (featureId) => {
    const id = normalizeFeatureSelectionId(featureId);
    if (!id) return null;
    const bySvgId = featuresBySvgId.value;
    return bySvgId.get(id) ||
      extractedFeatures.value.find((feature) => normalizeFeatureSelectionId(feature?.svg_id) === id) ||
      null;
  };

  /** @param {Element | null | undefined} svg */
  const getSelectableElementsById = (svg) => {
    /** @type {Map<string, Element[]>} */
    const indexed = new Map();
    if (!svg) return indexed;
    svg.querySelectorAll(SELECTABLE_FEATURE_SELECTOR).forEach((element) => {
      const id = getElementSelectionId(element);
      if (!id) return;
      if (!indexed.has(id)) indexed.set(id, []);
      // Set just above.
      /** @type {Element[]} */ (indexed.get(id)).push(element);
    });
    return indexed;
  };

  /**
   * @param {Element | null | undefined} svg
   * @param {string} [scope]
   */
  const buildRenderedFeatureGeometry = (svg, scope = '') => {
    const elementIndex = getSelectableElementsById(svg);
    /** @type {RenderedFeatureGeometry[]} */
    const geometries = [];
    extractedFeatures.value.forEach((feature) => {
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

  /** @param {unknown} featureId */
  const toggleFeatureSelection = (featureId) => {
    const id = normalizeFeatureSelectionId(featureId);
    if (!id || !getFeatureById(id)) return false;
    const nextIds = new Set(selectedFeatureIds.value);
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

  /**
   * @param {unknown} featureId
   * @param {{ additive?: boolean }} [options]
   */
  const selectFeatureRange = (featureId, { additive = false } = {}) => {
    const clickedId = normalizeFeatureSelectionId(featureId);
    const clicked = getFeatureById(clickedId);
    if (!clickedId || !clicked) return false;

    const anchorId = normalizeFeatureSelectionId(selectedFeatureAnchorId.value);
    const anchor = getFeatureById(anchorId);
    if (!anchorId || !anchor) {
      const nextIds = additive ? new Set(selectedFeatureIds.value) : new Set();
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

    // The record's axis orders the features as their bases are drawn (OV-2481). A Result saved before Python wrote
    // the axis keeps the estimate from the drawn shapes until its next Generate.
    const axisRecord = typeof anchor.record_key === 'string' ? axisRecordOf(anchor.record_key) : null;
    /** @type {string[]} */
    let rangeIds = [];
    if (axisRecord) {
      const { axis } = axisRecord;
      const placed = geometries.flatMap((entry) => {
        const position = partsMidpoint(axis.source, featureParts(entry.feature));
        return position === null ? [] : [{ id: entry.id, position }];
      });
      const from = placed.find((entry) => entry.id === anchorId);
      const to = placed.find((entry) => entry.id === clickedId);
      if (!from || !to) return false;
      rangeIds = placed
        .filter((entry) => onTheWay(axis, from.position, to.position, entry.position))
        .sort((left, right) => left.position - right.position || left.id.localeCompare(right.id))
        .map((entry) => entry.id);
    } else if (mode.value === 'circular') {
      const recordBounds = unionRects(geometries.map((entry) => entry.rect));
      if (!recordBounds) return false;
      const centerX = recordBounds.centerX;
      const centerY = recordBounds.centerY;
      /** @param {RenderedFeatureGeometry} geometry */
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

    const nextIds = additive ? new Set(selectedFeatureIds.value) : new Set();
    rangeIds.forEach((id) => nextIds.add(id));
    replaceSelectedFeatureIds(nextIds, { anchorId });
    closePreviewPopups();
    return true;
  };

  // A click on a feature with Ctrl or Cmd alone toggles it; otherwise the range from the anchor, added with Ctrl or Cmd.
  /** @param {string} featureId @param {{ ctrlKey: boolean, metaKey: boolean, shiftKey: boolean }} modifiers */
  const selectFeatureByClick = (featureId, { ctrlKey, metaKey, shiftKey }) => {
    const additive = Boolean(ctrlKey || metaKey);
    return additive && !shiftKey ? toggleFeatureSelection(featureId) : selectFeatureRange(featureId, { additive });
  };

  /**
   * @param {Element} svg
   * @param {Edges} marqueeRect
   * @param {boolean} [additive]
   */
  const selectIntersectingFeatures = (svg, marqueeRect, additive = false) => {
    const nextIds = additive ? new Set(selectedFeatureIds.value) : new Set();
    const geometries = buildRenderedFeatureGeometry(svg);
    geometries.forEach((entry) => {
      if (rectsIntersect(entry.rect, marqueeRect)) {
        nextIds.add(entry.id);
      }
    });
    const firstNewId = Array.from(nextIds).find(Boolean) || '';
    replaceSelectedFeatureIds(nextIds, {
      anchorId: firstNewId || selectedFeatureAnchorId.value,
      status: nextIds.size > 0 ? `${numberText(nextIds.size)} features selected.` : 'No rendered features in selection.'
    });
    closePreviewPopups();
  };

  /**
   * @param {PointerEvent | null | undefined} event
   * @param {Element | null | undefined} svg
   */
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
    const releaseText = holdTextSelection();
    const releaseCapture = capturePointer(svg, event.pointerId);
    releaseGesture = () => {
      releaseCapture();
      releaseText();
      releaseGesture = () => {};
    };
    return true;
  };

  /** @param {PointerEvent} event */
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

  /** @param {PointerEvent} event */
  const commitMarqueePointer = (event) => {
    if (!pendingDrag || event.pointerId !== pendingDrag.pointerId) return false;
    const wasActive = featureSelectionDrag.active;
    featureSelectionDrag.currentX = Number.isFinite(event.clientX) ? event.clientX : featureSelectionDrag.currentX;
    featureSelectionDrag.currentY = Number.isFinite(event.clientY) ? event.clientY : featureSelectionDrag.currentY;
    releaseGesture();
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
    releaseGesture();
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

  /** @param {KeyboardEvent} event */
  const handleEscape = (event) => {
    if (event.key !== 'Escape') return;
    if (activeToolbarDrag) {
      stopToolbarDrag();
      event.preventDefault();
      return;
    }
    if (featureSelectionDrag.active || pendingDrag) {
      cancelMarqueePointer();
      event.preventDefault();
      return;
    }
    // The key event of an element; `closest?.` covers a target without it.
    const target = /** @type {Element | null} */ (event.target);
    if (target?.closest?.('.feature-popup, .pairwise-match-popup, input, textarea, select, [contenteditable="true"]')) {
      return;
    }
    if (selectedFeatureIds.value.size > 0) {
      clearFeatureSelection({ clearStatus: true });
      event.preventDefault();
    }
  };

  const handleWindowBlur = () => {
    if (featureSelectionDrag.active || pendingDrag) cancelMarqueePointer();
    if (activeToolbarDrag) stopToolbarDrag();
  };

  if (typeof onMounted === 'function' && typeof onUnmounted === 'function') {
    onMounted(() => {
      document.addEventListener('keydown', handleEscape);
      window.addEventListener('blur', handleWindowBlur);
    });
    onUnmounted(() => {
      document.removeEventListener('keydown', handleEscape);
      window.removeEventListener('blur', handleWindowBlur);
      stopToolbarDrag();
      if (statusTimeoutId !== null) window.clearTimeout(statusTimeoutId);
    });
  }

  return {
    clearFeatureSelection,
    consumeSuppressNextClick,
    featureSelectionMarqueeStyle,
    featureSelectionToolbarStyle,
    markPlainFeatureClick,
    replaceSelectedFeatureIds,
    selectFeatureByClick,
    startToolbarDrag,
    startMarqueePointer,
    moveMarqueePointer,
    commitMarqueePointer,
    cancelMarqueePointer,
    syncFeatureSelectionClasses
  };
};
