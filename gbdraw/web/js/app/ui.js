export const createPanZoom = (state) => {
  const { zoom, isPanning, panStart, canvasPan, canvasContainerRef } = state;

  const handleWheel = (event) => {
    const delta = event.deltaY > 0 ? -0.1 : 0.1;
    const newZoom = Math.max(0.1, Math.min(5, zoom.value + delta));
    zoom.value = Math.round(newZoom * 10) / 10;
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

    isPanning.value = true;
    panStart.x = event.clientX;
    panStart.y = event.clientY;
    panStart.panX = canvasPan.x;
    panStart.panY = canvasPan.y;
    container.style.cursor = 'grabbing';
  };

  const doPan = (event) => {
    if (!isPanning.value) return;
    const container = canvasContainerRef.value;
    if (!container) return;

    const dx = event.clientX - panStart.x;
    const dy = event.clientY - panStart.y;
    canvasPan.x = panStart.panX + dx;
    canvasPan.y = panStart.panY + dy;
  };

  const endPan = () => {
    isPanning.value = false;
    const container = canvasContainerRef.value;
    if (container) {
      container.style.cursor = 'grab';
    }
  };

  return { handleWheel, startPan, doPan, endPan };
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
    showCanvasControls,
    showLegendPanel,
    showFeaturePanel
  } = state;

  const closeFeaturePopup = (e) => {
    if (!e.target.closest('.fixed.z-50')) {
      if (clickedFeature.value) clickedFeature.value = null;
    }
  };

  const handleEscapeKey = (e) => {
    if (e.key === 'Escape') {
      if (clickedFeature.value) clickedFeature.value = null;
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
