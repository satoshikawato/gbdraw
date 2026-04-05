export const createPanZoom = (state) => {
  const { zoom, isPanning, panStart, canvasPan, canvasContainerRef, svgContainer } = state;
  let panFrameId = null;
  let pendingPanPointer = null;

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

  const resetPreviewViewport = ({ resetZoom = false } = {}) => {
    cancelPanFrame();
    pendingPanPointer = null;
    isPanning.value = false;
    panStart.x = 0;
    panStart.y = 0;
    panStart.panX = 0;
    panStart.panY = 0;
    canvasPan.x = 0;
    canvasPan.y = 0;
    if (resetZoom) {
      zoom.value = 1.0;
    }

    const container = canvasContainerRef.value;
    if (container) {
      container.style.cursor = 'grab';
    }
    applyPreviewTransform(0, 0, zoom.value, false);
  };

  const handleWheel = (event) => {
    const delta = event.deltaY > 0 ? -0.1 : 0.1;
    const newZoom = Math.max(0.1, Math.min(5, zoom.value + delta));
    zoom.value = Math.round(newZoom * 10) / 10;
    applyPreviewTransform(canvasPan.x, canvasPan.y, zoom.value, isPanning.value);
  };

  const startPan = (event) => {
    if (event.button !== 0) return;
    const container = canvasContainerRef.value;
    if (!container) return;

    const target = event.target;
    const closestGroup = target.closest?.('g[id]');
    if (closestGroup) {
      const groupId = closestGroup.id;
      if (
        groupId === 'legend' ||
        groupId === 'feature_legend' ||
        groupId === 'pairwise_legend' ||
        groupId === 'horizontal_legend' ||
        groupId === 'vertical_legend' ||
        groupId.startsWith('comparison') ||
        groupId.startsWith('f') ||
        target.closest('svg')
      ) {
        return;
      }
    }
    if (target.tagName === 'path' && target.closest('svg')) {
      return;
    }

    cancelPanFrame();
    pendingPanPointer = null;
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
  };

  return { handleWheel, startPan, doPan, endPan, resetPreviewViewport };
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

export const setupGlobalUiEvents = ({ state, onMounted, onUnmounted }) => {
  const {
    clickedFeature,
    clickedLabel,
    showCanvasControls,
    showLegendPanel,
    showFeaturePanel
  } = state;

  const closeFeaturePopup = (e) => {
    if (!e.target.closest('.feature-popup') && !e.target.closest('.label-popup')) {
      if (clickedFeature.value) clickedFeature.value = null;
      if (clickedLabel.value) clickedLabel.value = null;
    }
  };

  const handleEscapeKey = (e) => {
    if (e.key === 'Escape') {
      if (clickedFeature.value) clickedFeature.value = null;
      if (clickedLabel.value) clickedLabel.value = null;
      if (showCanvasControls.value) showCanvasControls.value = false;
      if (showLegendPanel.value) showLegendPanel.value = false;
      if (showFeaturePanel.value) showFeaturePanel.value = false;
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
