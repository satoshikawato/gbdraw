import { recordStructuralMetric } from '../services/runtime-test-hooks.js';

// Keep one burst alive through the existing 0.2 s transform transition. This
// also absorbs browser delivery delay between nominally 80 ms wheel inputs.
const WHEEL_BURST_QUIET_MS = 220;
const WHEEL_TRANSITION_FALLBACK_MS = 260;

export const createPanZoom = (state) => {
  const { zoom, layoutRepositionMode, isPanning, panStart, canvasPan, canvasContainerRef, svgContainer } = state;
  let panFrameId = null;
  let pendingPanPointer = null;
  const previewTransformInteractionSources = new Set();
  let wheelBurstTimerId = null;
  let wheelFallbackTimerId = null;
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

  const resetPreviewViewport = ({ resetZoom = false, pan = null } = {}) => {
    cancelPreviewTransformInteraction({ reconcile: false });
    panStart.x = 0;
    panStart.y = 0;
    panStart.panX = 0;
    panStart.panY = 0;
    canvasPan.x = Number(pan?.x) || 0;
    canvasPan.y = Number(pan?.y) || 0;
    if (resetZoom) {
      zoom.value = 1.0;
    }

    const container = canvasContainerRef.value;
    if (container) {
      container.style.cursor = 'grab';
    }
    applyPreviewTransform(canvasPan.x, canvasPan.y, zoom.value, false);
  };

  const handleWheel = (event) => {
    beginWheelInteraction(event);
    const delta = event.deltaY > 0 ? -0.1 : 0.1;
    const newZoom = Math.max(0.1, Math.min(5, zoom.value + delta));
    zoom.value = Math.round(newZoom * 10) / 10;
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
      if (groupId.startsWith('f')) {
        return;
      }
      if (isLayoutRepositionModeEnabled() && target.closest('svg')) return;
    }
    if (isLayoutRepositionModeEnabled() && target.tagName === 'path' && target.closest('svg')) {
      return;
    }

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
    const finalPointer =
      typeof event?.clientX === 'number' && typeof event?.clientY === 'number'
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
    cancelPreviewTransformInteraction,
    disposePanZoom,
    previewTransformInteraction
  };
};

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

  const closeFeaturePopup = (e) => {
    if (!e.target.closest('.feature-popup') && !e.target.closest('.pairwise-match-popup') && !e.target.closest('.label-popup')) {
      if (clickedFeature.value) clickedFeature.value = null;
      if (clickedPairwiseMatch?.value) clickedPairwiseMatch.value = null;
      if (clickedLabel.value) clickedLabel.value = null;
    }
  };

  const handleEscapeKey = (e) => {
    if (e.key === 'Escape') {
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
