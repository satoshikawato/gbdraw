// @ts-check
import { recordStructuralMetric } from '../services/runtime-test-hooks.js';

// Keep one burst alive through the existing 0.2 s transform transition. This
// also absorbs browser delivery delay between nominally 80 ms wheel inputs.
const WHEEL_BURST_QUIET_MS = 220;
const WHEEL_TRANSITION_FALLBACK_MS = 260;
// The preview zoom range; the zoom-in and zoom-out buttons repeat it in index.html.
const MIN_ZOOM = 0.1;
const MAX_ZOOM = 5;
// Fit leaves the canvas padding (p-2) free around the diagram.
const FIT_MARGIN_PX = 8;

/** @param {number} value The zoom range and the wheel's 0.1 steps. */
const clampZoom = (value) => Math.round(Math.max(MIN_ZOOM, Math.min(MAX_ZOOM, value)) * 10) / 10;

/**
 * The width the open Editor drawer covers at the right of the canvas: the
 * stylesheet's `--preview-editor-reserve`, which the search and controls keep
 * free too (0 when the drawer is closed or docked below the canvas).
 * @param {HTMLElement} container The preview canvas.
 */
const editorReserveWidth = (container) => {
  const probe = container.ownerDocument?.createElement?.('div');
  if (!probe) return 0;
  probe.style.cssText = 'position:absolute;visibility:hidden;height:0;width:var(--preview-editor-reserve,0px)';
  container.appendChild(probe);
  const { width } = probe.getBoundingClientRect();
  probe.remove();
  return width;
};

/** @param {Record<string, any>} state Shape owned by state.js. */
export const createPanZoom = (state) => {
  const { zoom, layoutRepositionMode, isPanning, panStart, canvasPan, canvasContainerRef, svgContainer } = state;
  /** @type {number | null} */
  let panFrameId = null;
  /** @type {{ x: number, y: number } | null} */
  let pendingPanPointer = null;
  const previewTransformInteractionSources = new Set();
  /** @type {number | null} */
  let wheelBurstTimerId = null;
  /** @type {number | null} */
  let wheelFallbackTimerId = null;
  /** @type {HTMLElement | null} */
  let wheelTransitionTarget = null;
  let wheelBurstComplete = false;
  let wheelTransitionComplete = false;
  const previewTransformInteractionListeners = new Set();

  const previewTransformInteraction = Object.freeze({
    isActive: () => previewTransformInteractionSources.size > 0,
    subscribe(listener) {
      if (typeof listener !== 'function') return () => {};
      previewTransformInteractionListeners.add(listener);
      return () => previewTransformInteractionListeners.delete(listener);
    }
  });

  const notifyPreviewTransformInteraction = (active, kind, event, reconcile) => {
    previewTransformInteractionListeners.forEach((listener) => {
      listener({ active, kind, event, reconcile });
    });
  };

  const beginPreviewTransformInteraction = (kind, event) => {
    const wasActive = previewTransformInteractionSources.size > 0;
    previewTransformInteractionSources.add(kind);
    if (wasActive) return false;
    recordStructuralMetric('previewTransformInteractionStartCount', 1, { kind });
    notifyPreviewTransformInteraction(true, kind, event, false);
    return true;
  };

  const endPreviewTransformInteraction = (kind, { reconcile = true } = {}) => {
    if (!previewTransformInteractionSources.delete(kind)) return false;
    if (previewTransformInteractionSources.size > 0) return false;
    recordStructuralMetric('previewTransformInteractionEndCount', 1, { kind });
    notifyPreviewTransformInteraction(false, kind, null, reconcile);
    return true;
  };

  const clearWheelTimers = () => {
    if (wheelBurstTimerId !== null) {
      window.clearTimeout(wheelBurstTimerId);
      wheelBurstTimerId = null;
    }
    if (wheelFallbackTimerId !== null) {
      window.clearTimeout(wheelFallbackTimerId);
      wheelFallbackTimerId = null;
    }
  };

  const detachWheelTransitionListener = () => {
    wheelTransitionTarget?.removeEventListener?.('transitionend', handleWheelTransitionEnd);
    wheelTransitionTarget = null;
  };

  const finishWheelInteraction = ({ reconcile = true } = {}) => {
    if (!previewTransformInteractionSources.has('wheel')) return false;
    clearWheelTimers();
    detachWheelTransitionListener();
    wheelBurstComplete = false;
    wheelTransitionComplete = false;
    return endPreviewTransformInteraction('wheel', { reconcile });
  };

  function handleWheelTransitionEnd(event) {
    if (
      !previewTransformInteractionSources.has('wheel')
      || event?.target !== wheelTransitionTarget
      || event?.propertyName !== 'transform'
    ) return;
    wheelTransitionComplete = true;
    if (wheelBurstComplete) finishWheelInteraction();
  }

  const scheduleWheelInteractionEnd = () => {
    clearWheelTimers();
    wheelBurstComplete = false;
    wheelTransitionComplete = false;
    wheelBurstTimerId = window.setTimeout(() => {
      wheelBurstTimerId = null;
      wheelBurstComplete = true;
      if (wheelTransitionComplete) finishWheelInteraction();
    }, WHEEL_BURST_QUIET_MS);
    // A clamped zoom may not produce transitionend. Bound that path just past
    // the existing 0.2 s transform transition without changing the transition.
    wheelFallbackTimerId = window.setTimeout(() => {
      wheelFallbackTimerId = null;
      finishWheelInteraction();
    }, WHEEL_TRANSITION_FALLBACK_MS);
  };

  const beginWheelInteraction = (event) => {
    const isNewWheelInteraction = !previewTransformInteractionSources.has('wheel');
    beginPreviewTransformInteraction('wheel', event);
    if (isNewWheelInteraction) {
      wheelTransitionTarget = svgContainer.value;
      wheelTransitionTarget?.addEventListener?.('transitionend', handleWheelTransitionEnd);
    }
    scheduleWheelInteractionEnd();
  };

  const cancelPreviewTransformInteraction = ({ reconcile = false } = {}) => {
    const wasActive = previewTransformInteractionSources.size > 0;
    cancelPanFrame();
    pendingPanPointer = null;
    isPanning.value = false;
    if (canvasContainerRef.value) canvasContainerRef.value.style.cursor = 'grab';
    clearWheelTimers();
    detachWheelTransitionListener();
    wheelBurstComplete = false;
    wheelTransitionComplete = false;
    previewTransformInteractionSources.clear();
    if (!wasActive) return false;
    recordStructuralMetric('previewTransformInteractionEndCount', 1, { kind: 'cancel' });
    notifyPreviewTransformInteraction(false, 'cancel', null, reconcile);
    return true;
  };

  const isLayoutRepositionModeEnabled = () => Boolean(layoutRepositionMode?.value);

  const isFormControlTarget = (target) =>
    Boolean(target?.closest?.('button, input, textarea, select, a, [role="button"]'));

  const isSvgEditingTarget = (target) =>
    Boolean(
      target?.closest?.(
        [
          'text[data-label-editable="true"]',
          '[data-gbdraw-feature-id]',
          'path[id^="f"]',
          'polygon[id^="f"]',
          'rect[id^="f"]',
          '[data-gbdraw-pairwise-match-id]',
          '[data-match-kind]',
          '[data-pairwise-match-style]',
          '[data-collinearity-block-id]',
          '[data-collinear-group-scope]'
        ].join(', ')
      )
    );

  const cancelPanFrame = () => {
    if (panFrameId !== null) {
      cancelAnimationFrame(panFrameId);
      panFrameId = null;
    }
  };

  const applyPreviewTransform = (panX, panY, zoomLevel, disableTransition = isPanning.value) => {
    if (!svgContainer.value) return;
    svgContainer.value.style.transform = `translate(${panX}px, ${panY}px) scale(${zoomLevel})`;
    svgContainer.value.style.transformOrigin = 'top center';
    svgContainer.value.style.transition = disableTransition ? 'none' : 'transform 0.2s';
    svgContainer.value.style.willChange = disableTransition ? 'transform' : '';
  };

  const getPanPosition = (clientX, clientY) => {
    const dx = clientX - panStart.x;
    const dy = clientY - panStart.y;
    return {
      x: panStart.panX + dx,
      y: panStart.panY + dy
    };
  };

  const flushPanUpdate = (clientX, clientY) => {
    const nextPan = getPanPosition(clientX, clientY);
    applyPreviewTransform(nextPan.x, nextPan.y, zoom.value, true);
    return nextPan;
  };

  /**
   * Replace the whole preview viewport; Reset and Fit both end here.
   * @param {{ x: number, y: number }} pan
   * @param {number} zoomLevel
   */
  const setPreviewViewport = (pan, zoomLevel) => {
    cancelPreviewTransformInteraction({ reconcile: false });
    panStart.x = 0;
    panStart.y = 0;
    panStart.panX = 0;
    panStart.panY = 0;
    canvasPan.x = pan.x;
    canvasPan.y = pan.y;
    zoom.value = zoomLevel;
    applyPreviewTransform(canvasPan.x, canvasPan.y, zoom.value, false);
  };

  /** @param {{ resetZoom?: boolean, pan?: { x?: number, y?: number } | null }} [options] */
  const resetPreviewViewport = ({ resetZoom = false, pan = null } = {}) => setPreviewViewport(
    { x: Number(pan?.x) || 0, y: Number(pan?.y) || 0 },
    resetZoom ? 1.0 : zoom.value
  );

  // Fit (UI-11): the largest whole-percent zoom that shows the whole diagram in
  // the visible preview frame (the canvas left of an open Editor drawer), and
  // the pan that centres it under the top-center origin.
  const fitPreviewToViewport = () => {
    const container = canvasContainerRef.value;
    const surface = svgContainer.value;
    const svg = surface?.querySelector?.('svg');
    if (!container || !svg) return;
    // Measure the current transform itself, not a frame of its transition.
    applyPreviewTransform(canvasPan.x, canvasPan.y, zoom.value, true);
    const surfaceBox = surface.getBoundingClientRect();
    const svgBox = svg.getBoundingClientRect();
    const frameBox = container.getBoundingClientRect();
    const frameWidth = container.clientWidth - editorReserveWidth(container);
    const width = svgBox.width / zoom.value;
    const height = svgBox.height / zoom.value;
    if (!(width > 0 && height > 0)) return;
    const largest = Math.min(
      (frameWidth - 2 * FIT_MARGIN_PX) / width,
      (container.clientHeight - 2 * FIT_MARGIN_PX) / height
    );
    // Round down to a whole percent so the diagram fills the frame without
    // overflowing it; the wheel and buttons step by 0.1 from there.
    const nextZoom = Math.max(MIN_ZOOM, Math.min(MAX_ZOOM, Math.floor(largest * 100 + 1e-9) / 100));
    // The transform origin without the current pan, and the diagram centre's
    // unscaled offset from it.
    const originX = surfaceBox.left + surfaceBox.width / 2 - canvasPan.x;
    const originY = surfaceBox.top - canvasPan.y;
    const centerX = (svgBox.left + svgBox.width / 2 - surfaceBox.left - surfaceBox.width / 2) / zoom.value;
    const centerY = (svgBox.top + svgBox.height / 2 - surfaceBox.top) / zoom.value;
    setPreviewViewport({
      x: frameBox.left + container.clientLeft + frameWidth / 2 - originX - nextZoom * centerX,
      y: frameBox.top + container.clientTop + container.clientHeight / 2 - originY - nextZoom * centerY
    }, nextZoom);
  };

  const handleWheel = (event) => {
    beginWheelInteraction(event);
    zoom.value = clampZoom(zoom.value + (event.deltaY > 0 ? -0.1 : 0.1));
    applyPreviewTransform(canvasPan.x, canvasPan.y, zoom.value, isPanning.value);
  };

  const startPan = (event) => {
    if (event.button !== 0) return;
    if (event.shiftKey && svgContainer.value?.querySelector?.('svg')) return;
    const container = canvasContainerRef.value;
    if (!container) return;

    const target = event.target;
    if (isFormControlTarget(target) || isSvgEditingTarget(target)) return;

    const closestGroup = target.closest?.('g[id]');
    if (closestGroup) {
      const groupId = closestGroup.id;
      if (groupId.startsWith('f') && !target.closest('.gbdraw-preview-layout-target')) {
        return;
      }
      if (isLayoutRepositionModeEnabled() && target.closest('svg')) return;
    }
    if (isLayoutRepositionModeEnabled() && target.tagName === 'path' && target.closest('svg')) {
      return;
    }

    // Own the accepted pan instead of starting native SVG text selection/drag.
    event.preventDefault?.();
    cancelPanFrame();
    pendingPanPointer = null;
    beginPreviewTransformInteraction('pan', event);
    isPanning.value = true;
    panStart.x = event.clientX;
    panStart.y = event.clientY;
    panStart.panX = canvasPan.x;
    panStart.panY = canvasPan.y;
    container.style.cursor = 'grabbing';
    applyPreviewTransform(canvasPan.x, canvasPan.y, zoom.value, true);
  };

  const doPan = (event) => {
    if (!isPanning.value) return;
    pendingPanPointer = { x: event.clientX, y: event.clientY };
    if (panFrameId !== null) return;
    panFrameId = requestAnimationFrame(() => {
      panFrameId = null;
      if (!isPanning.value || !pendingPanPointer) return;
      flushPanUpdate(pendingPanPointer.x, pendingPanPointer.y);
    });
  };

  const endPan = (event) => {
    const wasPanning = isPanning.value;
    // Cancellation and capture loss may carry (0, 0), not a pointer position.
    const finalPointer =
      event?.type !== 'pointercancel' && event?.type !== 'lostpointercapture'
      && Number.isFinite(event?.clientX) && Number.isFinite(event?.clientY)
        ? { x: event.clientX, y: event.clientY }
        : pendingPanPointer;
    cancelPanFrame();

    if (isPanning.value && finalPointer) {
      const nextPan = flushPanUpdate(finalPointer.x, finalPointer.y);
      canvasPan.x = nextPan.x;
      canvasPan.y = nextPan.y;
    }

    pendingPanPointer = null;
    isPanning.value = false;
    const container = canvasContainerRef.value;
    if (container) {
      container.style.cursor = 'grab';
    }
    applyPreviewTransform(canvasPan.x, canvasPan.y, zoom.value, false);
    if (wasPanning) endPreviewTransformInteraction('pan');
  };

  const disposePanZoom = () => {
    cancelPreviewTransformInteraction({ reconcile: false });
    previewTransformInteractionListeners.clear();
  };

  return {
    handleWheel,
    startPan,
    doPan,
    endPan,
    resetPreviewViewport,
    fitPreviewToViewport,
    cancelPreviewTransformInteraction,
    disposePanZoom,
    previewTransformInteraction
  };
};

/** @param {Record<string, any>} state Shape owned by state.js. */
export const createSidebarResize = (state) => {
  const { sidebarWidth, isResizing } = state;

  const doResize = (event) => {
    if (!isResizing.value) return;
    const newWidth = event.clientX - 16;
    sidebarWidth.value = Math.max(240, Math.min(500, newWidth));
  };

  const stopResizing = () => {
    isResizing.value = false;
    document.removeEventListener('mousemove', doResize);
    document.removeEventListener('mouseup', stopResizing);
  };

  const startResizing = () => {
    isResizing.value = true;
    document.addEventListener('mousemove', doResize);
    document.addEventListener('mouseup', stopResizing);
  };

  return { startResizing };
};

/**
 * @typedef {Object} GlobalUiEventsOptions
 * @property {Record<string, any>} state Shape owned by state.js.
 * @property {(callback: () => void) => void} onMounted
 * @property {(callback: () => void) => void} onUnmounted
 * @property {() => void} closeRightDrawer
 */

/** @param {GlobalUiEventsOptions} options */
export const setupGlobalUiEvents = ({
  state,
  onMounted,
  onUnmounted,
  closeRightDrawer
}) => {
  const {
    clickedFeature,
    clickedPairwiseMatch,
    clickedLabel,
    showCanvasControls
  } = state;

  // UI-02: a modal dialog answers Escape and its own clicks (its Cancel), so
  // nothing behind it closes. The target is read through closest() because
  // a choice unmounts the dialog before this document listener runs.
  const modalDialogOpen = () => Boolean(document.querySelector('[role="dialog"][aria-modal="true"]'));

  const closeFeaturePopup = (e) => {
    if (
      !e.target.closest('[data-modal-overlay], [aria-modal="true"]')
      && !e.target.closest('.feature-popup')
      && !e.target.closest('.pairwise-match-popup')
      && !e.target.closest('.label-popup')
      && !e.target.closest('[data-similarity-alignment-overlay]')
    ) {
      if (clickedFeature.value) clickedFeature.value = null;
      if (clickedPairwiseMatch?.value) clickedPairwiseMatch.value = null;
      if (clickedLabel.value) clickedLabel.value = null;
    }
  };

  const handleEscapeKey = (e) => {
    if (e.key === 'Escape' && !modalDialogOpen()) {
      if (clickedFeature.value) clickedFeature.value = null;
      if (clickedPairwiseMatch?.value) clickedPairwiseMatch.value = null;
      if (clickedLabel.value) clickedLabel.value = null;
      if (showCanvasControls.value) showCanvasControls.value = false;
      closeRightDrawer();
    }
  };

  onMounted(() => {
    document.addEventListener('click', closeFeaturePopup);
    document.addEventListener('keydown', handleEscapeKey);
  });

  onUnmounted(() => {
    document.removeEventListener('click', closeFeaturePopup);
    document.removeEventListener('keydown', handleEscapeKey);
  });
};
