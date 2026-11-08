// @ts-check
// The popups of the Result: the feature popup, the match popup, and the
// settings section shown beside the Result. Each is placed again once
// rendered, dragged by its handle, and kept inside the window; the feature and
// match popups are resized from their corner. The composition root passes the
// popup state of state.js and the reactions of its owners as ports.
import { placePopupNearPoint, placePopupWithin, popupViewportBounds } from '../utils/popup-bounds.js';
/** @import { PopupPlacement } from '../utils/popup-bounds.js' */
/** @import { PopupPosition, WebState } from '../state.js' */
/** @import { FeatureListRow } from './pairwise-match-popup.js' */
/** @import { MatchFramePopup } from './preview-item-popup.js' */
import { matchPopupOf } from './preview-item-popup.js';

const { computed, onMounted, onUnmounted, watch } = window.Vue;

/** @typedef {{ value: HTMLElement | null }} PopupElementRef */
/** @typedef {{ active: boolean, offsetX: number, offsetY: number }} PopupDragState */
/** @typedef {{ width: number, height: number }} PopupSize */
/** @typedef {{ active: boolean, startX: number, startY: number, startWidth: number, startHeight: number }} PopupResizeState */
/** @typedef {{ minWidth: number, minHeight: number, maxWidth: number, maxHeight: number }} PopupConstraints */
/** @typedef {{ start: (event: MouseEvent & { target: Element }) => void, end: () => void, keepInWindow: () => void }} PopupDrag */
/** @typedef {(opened: MatchFramePopup | null) => void} MarkSelectedPairwiseMatch The SVG owner's mark of the open match's bands. */
/** @typedef {(feature: object, event: MouseEvent) => void} OpenFeatureEditorForFeature */

/** @param {Element | null} target */
const isInteractiveTarget = (target) => {
  if (!target) return false;
  return Boolean(target.closest('input, textarea, select, button, label, a, [data-nodrag="true"]'));
};

const FEATURE_POPUP_MARGIN = 12;
const FEATURE_POPUP_RICH_MIN_WIDTH = 360;
const FEATURE_POPUP_SIMPLE_MIN_WIDTH = 300;
const FEATURE_POPUP_MIN_HEIGHT = 220;
const PAIRWISE_MATCH_POPUP_MIN_WIDTH = 280;
const PAIRWISE_MATCH_POPUP_MIN_HEIGHT = 180;

/** @param {unknown} value @param {number} min @param {number} max */
const clampNumber = (value, min, max) => {
  const safeMin = Number.isFinite(min) ? min : 0;
  const safeMax = Number.isFinite(max) ? Math.max(safeMin, max) : safeMin;
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return safeMin;
  return Math.min(Math.max(numeric, safeMin), safeMax);
};

// A popup's size limits at its position, inside the popup bounds.
/**
 * @param {number} left
 * @param {number} top
 * @param {number} desiredMinWidth
 * @param {number} desiredMinHeight
 */
const popupConstraints = (left, top, desiredMinWidth, desiredMinHeight) => {
  const bounds = popupViewportBounds();
  const availableWidth = Math.max(1, bounds.right - (FEATURE_POPUP_MARGIN * 2));
  const availableHeight = Math.max(1, bounds.bottom - bounds.top - (FEATURE_POPUP_MARGIN * 2));
  const minWidth = Math.min(desiredMinWidth, availableWidth);
  const minHeight = Math.min(desiredMinHeight, availableHeight);
  return {
    minWidth,
    minHeight,
    maxWidth: Math.max(minWidth, bounds.right - left - FEATURE_POPUP_MARGIN),
    maxHeight: Math.max(minHeight, bounds.bottom - top - FEATURE_POPUP_MARGIN)
  };
};

/**
 * @param {PopupPlacement} pos
 * @param {{ width: number, height: number }} size
 * @param {{ maxHeight: number }} constraints
 */
const popupStyle = (pos, size, constraints) => ({
  top: `${pos.y}px`,
  left: `${pos.x}px`,
  maxHeight: `${pos.maxHeight > 0 ? Math.min(pos.maxHeight, constraints.maxHeight) : constraints.maxHeight}px`,
  ...(size.width > 0 ? { width: `${size.width}px` } : {}),
  ...(size.height > 0 ? { height: `${size.height}px` } : {})
});

// OV-371: a popup opened at a clicked point was placed with an expected
// size; once rendered it is placed again with its own size and content height.
/**
 * @param {PopupPosition} pos
 * @param {HTMLElement | null} element
 * @param {string} bodySelector
 */
const placeRenderedPopup = (pos, element, bodySelector) => {
  const point = pos.anchor;
  pos.anchor = null;
  if (!point || !element) return;
  const body = element.querySelector(bodySelector);
  const height = element.offsetHeight + (body ? body.scrollHeight - body.clientHeight : 0);
  Object.assign(pos, placePopupNearPoint({ point, width: element.offsetWidth, height, bounds: popupViewportBounds() }));
};

// A press on the handle of an open popup drags it; every popup of the Result moves this way.
// A drag keeps the popup inside the popup bounds (OV-384: below the header), and a window that
// changed size moves an open popup the least distance back inside them (OV-626).
/**
 * @param {object} popup
 * @param {{ active: boolean, offsetX: number, offsetY: number }} popup.drag
 * @param {() => boolean} popup.isOpen
 * @param {(placement: PopupPlacement) => void} popup.moveTo The popup owner's move to a placement.
 * @param {() => HTMLElement | null} popup.element
 * @param {(target: Element) => boolean} popup.isHandle Whether a press on `target` drags the popup.
 */
const createPopupDrag = ({ drag, isOpen, moveTo, element, isHandle }) => {
  /** @param {MouseEvent} event */
  const move = (event) => {
    const popup = element();
    if (!drag.active || !isOpen()) return;
    moveTo(placePopupWithin({
      x: event.clientX - drag.offsetX, y: event.clientY - drag.offsetY,
      width: popup?.offsetWidth || 0, height: popup?.offsetHeight || 0, bounds: popupViewportBounds()
    }));
  };
  const end = () => {
    drag.active = false;
    document.removeEventListener('mousemove', move);
    document.removeEventListener('mouseup', end);
  };
  /** @param {MouseEvent & { target: Element }} event A press on the popup; its target is an element. */
  const start = (event) => {
    const popup = element();
    if (event.button !== 0 || !popup || !isOpen() || !isHandle(event.target)) return;
    const rect = popup.getBoundingClientRect();
    drag.active = true;
    drag.offsetX = event.clientX - rect.left;
    drag.offsetY = event.clientY - rect.top;
    document.addEventListener('mousemove', move);
    document.addEventListener('mouseup', end);
    event.preventDefault();
  };
  const keepInWindow = () => {
    const popup = element();
    if (!popup || !isOpen()) return;
    const { left: x, top: y } = popup.getBoundingClientRect();
    const within = placePopupWithin({ x, y, width: popup.offsetWidth, height: popup.offsetHeight, bounds: popupViewportBounds() });
    if (within.x !== x || within.y !== y) moveTo(within);
  };
  return { start, end, keepInWindow };
};

// The heading of the section in the settings popup is the popup's handle. It is the section's
// own disclosure heading, which does not fold the section there.
const SETTINGS_POPUP_HEADING = '#settings-popup > [data-settings-section] > details > summary';
/** @param {Element} target */
const isSettingsPopupHeading = (target) => !isInteractiveTarget(target) && Boolean(target.closest(SETTINGS_POPUP_HEADING));
/** @param {MouseEvent & { target: Element }} event A click in the popup; its target is an element. */
const keepSettingsPopupSectionOpen = (event) => {
  if (isSettingsPopupHeading(event.target)) event.preventDefault();
};

/**
 * The similarity group the match popup lists members of, and the gene pair of a collinear block picked in it.
 * @param {object} options
 * @param {{ value: MatchFramePopup | null }} options.clickedPairwiseMatch
 * @param {{ value: string }} options.selectedPairwiseBlockOrthogroupId
 * @param {OpenFeatureEditorForFeature} options.openFeatureEditorForFeature
 * @param {MarkSelectedPairwiseMatch} options.markSelectedPairwiseMatch
 */
const createPairwiseGroupView = ({
  clickedPairwiseMatch, selectedPairwiseBlockOrthogroupId, openFeatureEditorForFeature, markSelectedPairwiseMatch
}) => {
  // A collinear block carries the `blockOrthogroups` of its gene pairs, a similarity-group ribbon its `group`.
  const pairwiseBlockOrthogroups = computed(() => matchPopupOf(clickedPairwiseMatch.value)?.blockOrthogroups ?? []);
  const selectedPairwiseBlockOrthogroup = computed(() => {
    const selectedId = String(selectedPairwiseBlockOrthogroupId.value || '').trim();
    if (!selectedId) return null;
    return pairwiseBlockOrthogroups.value.find((/** @type {{ id?: unknown }} */ group) => String(group?.id || '').trim() === selectedId) || null;
  });
  // The group the match popup lists members of: a similarity-group ribbon's
  // own group, or the group picked in a collinear block's gene-pair table.
  const pairwiseGroupView = computed(() => {
    const match = matchPopupOf(clickedPairwiseMatch.value);
    return match?.block ? selectedPairwiseBlockOrthogroup.value : match?.group || null;
  });
  /** @param {{ id?: unknown } | null | undefined} group */
  const selectPairwiseBlockOrthogroup = (group) => {
    const id = String(group?.id || '').trim();
    selectedPairwiseBlockOrthogroupId.value = selectedPairwiseBlockOrthogroupId.value === id ? '' : id;
  };
  /** @param {FeatureListRow} row @param {MouseEvent} event */
  const openPairwiseFeatureRow = (row, event) => {
    if (!row?.canOpen || !row?.feature) return;
    openFeatureEditorForFeature(row.feature, event);
  };

  // The bands of the open match popup carry `gbdraw-match-selected`; the SVG
  // owner marks them from its band indexes (OV-937).
  watch(clickedPairwiseMatch, (/** @type {MatchFramePopup | null} */ match) => {
    selectedPairwiseBlockOrthogroupId.value = '';
    markSelectedPairwiseMatch(match);
  });
  return { pairwiseGroupView, selectPairwiseBlockOrthogroup, openPairwiseFeatureRow };
};

/**
 * The match popup's resize from its corner.
 * @param {object} options
 * @param {{ value: MatchFramePopup | null }} options.clickedPairwiseMatch
 * @param {PopupElementRef} options.pairwiseMatchPopupRef
 * @param {PopupPosition} options.clickedPairwiseMatchPos
 * @param {PopupResizeState} options.pairwiseMatchPopupResize
 * @param {PopupSize} options.pairwiseMatchPopupSize
 * @param {(left?: number, top?: number) => PopupConstraints} options.getPairwiseMatchPopupConstraints
 * @param {PopupDrag} options.pairwiseMatchPopupDragging
 */
const createPairwiseMatchPopupResize = ({
  clickedPairwiseMatch, pairwiseMatchPopupRef, clickedPairwiseMatchPos, pairwiseMatchPopupResize, pairwiseMatchPopupSize,
  getPairwiseMatchPopupConstraints, pairwiseMatchPopupDragging
}) => {
  /** @param {MouseEvent} event */
  const onPairwiseMatchPopupResize = (event) => {
    if (!pairwiseMatchPopupResize.active) return;
    const constraints = getPairwiseMatchPopupConstraints(clickedPairwiseMatchPos.x, clickedPairwiseMatchPos.y);
    const nextWidth = pairwiseMatchPopupResize.startWidth + (event.clientX - pairwiseMatchPopupResize.startX);
    const nextHeight = pairwiseMatchPopupResize.startHeight + (event.clientY - pairwiseMatchPopupResize.startY);
    pairwiseMatchPopupSize.width = clampNumber(nextWidth, constraints.minWidth, constraints.maxWidth);
    pairwiseMatchPopupSize.height = clampNumber(nextHeight, constraints.minHeight, constraints.maxHeight);
    event.preventDefault();
  };

  const endPairwiseMatchPopupResize = () => {
    if (!pairwiseMatchPopupResize.active) return;
    pairwiseMatchPopupResize.active = false;
    document.removeEventListener('mousemove', onPairwiseMatchPopupResize);
    document.removeEventListener('mouseup', endPairwiseMatchPopupResize);
  };

  /** @param {MouseEvent & { target: Element }} event A press on the popup; its target is an element. */
  const startPairwiseMatchPopupResize = (event) => {
    if (event.button !== 0) return;
    if (!clickedPairwiseMatch.value) return;
    const popup = pairwiseMatchPopupRef.value;
    if (!popup) return;
    const rect = popup.getBoundingClientRect();
    const constraints = getPairwiseMatchPopupConstraints(rect.left, rect.top);

    pairwiseMatchPopupDragging.end();

    clickedPairwiseMatchPos.maxHeight = 0;
    pairwiseMatchPopupResize.active = true;
    pairwiseMatchPopupResize.startX = event.clientX;
    pairwiseMatchPopupResize.startY = event.clientY;
    pairwiseMatchPopupResize.startWidth = clampNumber(rect.width, constraints.minWidth, constraints.maxWidth);
    pairwiseMatchPopupResize.startHeight = clampNumber(rect.height, constraints.minHeight, constraints.maxHeight);
    pairwiseMatchPopupSize.width = pairwiseMatchPopupResize.startWidth;
    pairwiseMatchPopupSize.height = pairwiseMatchPopupResize.startHeight;
    document.addEventListener('mousemove', onPairwiseMatchPopupResize);
    document.addEventListener('mouseup', endPairwiseMatchPopupResize);
    event.preventDefault();
    event.stopPropagation();
  };
  return startPairwiseMatchPopupResize;
};

/**
 * The feature popup's resize from its corner.
 * @param {object} options
 * @param {WebState['clickedFeature']} options.clickedFeature
 * @param {PopupElementRef} options.featurePopupRef
 * @param {PopupPosition} options.clickedFeaturePos
 * @param {PopupResizeState} options.featurePopupResize
 * @param {PopupSize} options.featurePopupSize
 * @param {(left?: number, top?: number) => PopupConstraints} options.getFeaturePopupConstraints
 * @param {PopupDrag} options.featurePopupDragging
 */
const createFeaturePopupResize = ({
  clickedFeature, featurePopupRef, clickedFeaturePos, featurePopupResize, featurePopupSize,
  getFeaturePopupConstraints, featurePopupDragging
}) => {
  /** @param {MouseEvent} event */
  const onFeaturePopupResize = (event) => {
    if (!featurePopupResize.active) return;
    const constraints = getFeaturePopupConstraints(clickedFeaturePos.x, clickedFeaturePos.y);
    const nextWidth = featurePopupResize.startWidth + (event.clientX - featurePopupResize.startX);
    const nextHeight = featurePopupResize.startHeight + (event.clientY - featurePopupResize.startY);
    featurePopupSize.width = clampNumber(nextWidth, constraints.minWidth, constraints.maxWidth);
    featurePopupSize.height = clampNumber(nextHeight, constraints.minHeight, constraints.maxHeight);
    event.preventDefault();
  };

  const endFeaturePopupResize = () => {
    if (!featurePopupResize.active) return;
    featurePopupResize.active = false;
    document.removeEventListener('mousemove', onFeaturePopupResize);
    document.removeEventListener('mouseup', endFeaturePopupResize);
  };

  /** @param {MouseEvent & { target: Element }} event A press on the popup; its target is an element. */
  const startFeaturePopupResize = (event) => {
    if (event.button !== 0) return;
    if (!clickedFeature.value) return;
    const popup = featurePopupRef.value;
    if (!popup) return;
    const rect = popup.getBoundingClientRect();
    const constraints = getFeaturePopupConstraints(rect.left, rect.top);

    featurePopupDragging.end();

    clickedFeaturePos.maxHeight = 0;
    featurePopupResize.active = true;
    featurePopupResize.startX = event.clientX;
    featurePopupResize.startY = event.clientY;
    featurePopupResize.startWidth = clampNumber(rect.width, constraints.minWidth, constraints.maxWidth);
    featurePopupResize.startHeight = clampNumber(rect.height, constraints.minHeight, constraints.maxHeight);
    featurePopupSize.width = featurePopupResize.startWidth;
    featurePopupSize.height = featurePopupResize.startHeight;
    document.addEventListener('mousemove', onFeaturePopupResize);
    document.addEventListener('mouseup', endFeaturePopupResize);
    event.preventDefault();
    event.stopPropagation();
  };
  return startFeaturePopupResize;
};

/**
 * The popup state of state.js and the owners' reactions the Result popups use.
 * @typedef {object} ResultPopupsPorts
 * @property {WebState} state The app state; the settings popup's place is `state.settingsPopup`.
 * @property {WebState['clickedFeature']} clickedFeature The open feature popup's payload, or null.
 * @property {PopupPosition} clickedFeaturePos
 * @property {{ value: MatchFramePopup | null }} clickedPairwiseMatch The open match popup's payload, or null.
 * @property {PopupPosition} clickedPairwiseMatchPos
 * @property {{ value: boolean }} richFeaturePopup The app-level rich popup preference.
 * @property {PopupElementRef} featurePopupRef
 * @property {PopupDragState} featurePopupDrag
 * @property {PopupSize} featurePopupSize
 * @property {PopupResizeState} featurePopupResize
 * @property {PopupElementRef} pairwiseMatchPopupRef
 * @property {PopupDragState} pairwiseMatchPopupDrag
 * @property {PopupSize} pairwiseMatchPopupSize
 * @property {PopupResizeState} pairwiseMatchPopupResize
 * @property {{ value: string }} selectedPairwiseBlockOrthogroupId The gene pair picked in a collinear block's table.
 * @property {(placement: PopupPlacement) => void} moveSettingsPopup The settings pane's one move of its popup.
 * @property {OpenFeatureEditorForFeature} openFeatureEditorForFeature
 * @property {MarkSelectedPairwiseMatch} markSelectedPairwiseMatch The SVG owner's mark of the open match's bands.
 */

/**
 * @param {ResultPopupsPorts} ports
 */
export const createResultPopups = ({
  state, clickedFeature, clickedFeaturePos, clickedPairwiseMatch, clickedPairwiseMatchPos, richFeaturePopup,
  featurePopupRef, featurePopupDrag, featurePopupSize, featurePopupResize,
  pairwiseMatchPopupRef, pairwiseMatchPopupDrag, pairwiseMatchPopupSize, pairwiseMatchPopupResize,
  selectedPairwiseBlockOrthogroupId, moveSettingsPopup, openFeatureEditorForFeature, markSelectedPairwiseMatch
}) => {
  // The popup's size follows the app-level rich popup preference.
  const getFeaturePopupConstraints = (left = clickedFeaturePos.x, top = clickedFeaturePos.y) => popupConstraints(
    left, top,
    richFeaturePopup.value === false ? FEATURE_POPUP_SIMPLE_MIN_WIDTH : FEATURE_POPUP_RICH_MIN_WIDTH,
    FEATURE_POPUP_MIN_HEIGHT
  );
  const getPairwiseMatchPopupConstraints = (left = clickedPairwiseMatchPos.x, top = clickedPairwiseMatchPos.y) => popupConstraints(
    left, top, PAIRWISE_MATCH_POPUP_MIN_WIDTH, PAIRWISE_MATCH_POPUP_MIN_HEIGHT
  );

  const featurePopupStyle = computed(() => popupStyle(clickedFeaturePos, featurePopupSize, getFeaturePopupConstraints()));
  const pairwiseMatchPopupStyle = computed(() => popupStyle(clickedPairwiseMatchPos, pairwiseMatchPopupSize, getPairwiseMatchPopupConstraints()));
  // The settings popup has the width of its stylesheet and the height of its section.
  const settingsPopupStyle = computed(() => {
    const pos = state.settingsPopup.value;
    return pos ? popupStyle(pos, { width: 0, height: 0 }, popupConstraints(pos.x, pos.y, 0, 0)) : null;
  });

  watch(clickedFeature, () => placeRenderedPopup(clickedFeaturePos, featurePopupRef.value, '.feature-popup-body'), { flush: 'post' });
  watch(clickedPairwiseMatch, () => placeRenderedPopup(clickedPairwiseMatchPos, pairwiseMatchPopupRef.value, '.pairwise-match-popup-body'), { flush: 'post' });

  const featurePopupDragging = createPopupDrag({
    drag: featurePopupDrag,
    isOpen: () => Boolean(clickedFeature.value),
    moveTo: (placement) => { Object.assign(clickedFeaturePos, placement); },
    element: () => featurePopupRef.value,
    isHandle: (target) => !featurePopupResize.active && !isInteractiveTarget(target)
  });
  const pairwiseMatchPopupDragging = createPopupDrag({
    drag: pairwiseMatchPopupDrag,
    isOpen: () => Boolean(clickedPairwiseMatch.value),
    moveTo: (placement) => { Object.assign(clickedPairwiseMatchPos, placement); },
    element: () => pairwiseMatchPopupRef.value,
    isHandle: (target) => !pairwiseMatchPopupResize.active && !isInteractiveTarget(target)
  });
  const settingsPopupDrag = createPopupDrag({
    drag: { active: false, offsetX: 0, offsetY: 0 },
    isOpen: () => Boolean(state.settingsPopup.value),
    moveTo: moveSettingsPopup,
    element: () => document.getElementById('settings-popup'),
    isHandle: isSettingsPopupHeading
  });
  const keepPopupsInWindow = () => [featurePopupDragging, pairwiseMatchPopupDragging, settingsPopupDrag].forEach((popup) => popup.keepInWindow());
  onMounted(() => window.addEventListener('resize', keepPopupsInWindow));
  onUnmounted(() => window.removeEventListener('resize', keepPopupsInWindow));

  const { pairwiseGroupView, selectPairwiseBlockOrthogroup, openPairwiseFeatureRow } = createPairwiseGroupView({
    clickedPairwiseMatch, selectedPairwiseBlockOrthogroupId, openFeatureEditorForFeature, markSelectedPairwiseMatch
  });
  const startPairwiseMatchPopupResize = createPairwiseMatchPopupResize({
    clickedPairwiseMatch, pairwiseMatchPopupRef, clickedPairwiseMatchPos, pairwiseMatchPopupResize, pairwiseMatchPopupSize,
    getPairwiseMatchPopupConstraints, pairwiseMatchPopupDragging
  });
  const startFeaturePopupResize = createFeaturePopupResize({
    clickedFeature, featurePopupRef, clickedFeaturePos, featurePopupResize, featurePopupSize,
    getFeaturePopupConstraints, featurePopupDragging
  });
  return {
    featurePopupStyle,
    pairwiseMatchPopupStyle,
    settingsPopupStyle,
    featurePopupDragging,
    pairwiseMatchPopupDragging,
    settingsPopupDrag,
    keepSettingsPopupSectionOpen,
    pairwiseGroupView,
    selectPairwiseBlockOrthogroup,
    openPairwiseFeatureRow,
    startFeaturePopupResize,
    startPairwiseMatchPopupResize
  };
};
