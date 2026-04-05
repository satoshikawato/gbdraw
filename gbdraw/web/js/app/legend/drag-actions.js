import { parseTransform } from './utils.js';

export const createLegendDragActions = ({ state, extractLegendEntries }) => {
  const {
    results,
    selectedResultIndex,
    svgContainer,
    legendDragging,
    legendDragStart,
    legendOriginalTransform,
    legendInitialTransform,
    legendCurrentOffset,
    zoom,
    skipCaptureBaseConfig
  } = state;
  let legendDragFrameId = null;
  let pendingLegendPointer = null;

  const cancelLegendDragFrame = () => {
    if (legendDragFrameId !== null) {
      cancelAnimationFrame(legendDragFrameId);
      legendDragFrameId = null;
    }
  };

  const applyLegendDragPosition = (clientX, clientY) => {
    if (!legendDragging.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const deltaX = (clientX - legendDragStart.x) / zoom.value;
    const deltaY = (clientY - legendDragStart.y) / zoom.value;

    const newX = legendOriginalTransform.value.x + deltaX;
    const newY = legendOriginalTransform.value.y + deltaY;

    legendGroup.setAttribute('transform', `translate(${newX}, ${newY})`);
    legendCurrentOffset.x = deltaX;
    legendCurrentOffset.y = deltaY;
  };

  const startLegendDrag = (e) => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    e.preventDefault();
    e.stopPropagation();

    cancelLegendDragFrame();
    pendingLegendPointer = null;
    legendDragging.value = true;
    legendDragStart.x = e.clientX;
    legendDragStart.y = e.clientY;

    const currentTransform = parseTransform(legendGroup.getAttribute('transform'));
    legendOriginalTransform.value = { ...currentTransform };
    legendGroup.style.willChange = 'transform';
  };

  const onLegendDrag = (e) => {
    if (!legendDragging.value) return;
    pendingLegendPointer = { x: e.clientX, y: e.clientY };
    if (legendDragFrameId !== null) return;
    legendDragFrameId = requestAnimationFrame(() => {
      legendDragFrameId = null;
      if (!pendingLegendPointer) return;
      applyLegendDragPosition(pendingLegendPointer.x, pendingLegendPointer.y);
    });
  };

  const endLegendDrag = (e) => {
    if (!legendDragging.value) return;
    const finalPointer =
      typeof e?.clientX === 'number' && typeof e?.clientY === 'number'
        ? { x: e.clientX, y: e.clientY }
        : pendingLegendPointer;
    cancelLegendDragFrame();
    if (finalPointer) {
      applyLegendDragPosition(finalPointer.x, finalPointer.y);
    }

    if (svgContainer.value) {
      const svg = svgContainer.value.querySelector('svg');
      const legendGroup = svg?.getElementById('legend');
      if (legendGroup) {
        legendGroup.style.willChange = '';
      }
    }

    pendingLegendPointer = null;
    legendDragging.value = false;
  };

  const resetLegendPositionOnly = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const initial = legendInitialTransform.value;
    legendGroup.setAttribute('transform', `translate(${initial.x}, ${initial.y})`);
    legendCurrentOffset.x = 0;
    legendCurrentOffset.y = 0;

    skipCaptureBaseConfig.value = true;
    const idx = selectedResultIndex.value;
    if (idx >= 0 && results.value.length > idx) {
      const serializer = new XMLSerializer();
      results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svg) };
    }
  };

  const resetLegendPosition = () => {
    resetLegendPositionOnly();
    extractLegendEntries();
  };

  const setupLegendDrag = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const initialTransform = parseTransform(legendGroup.getAttribute('transform'));
    legendInitialTransform.value = { ...initialTransform };

    legendGroup.style.cursor = 'move';
    legendGroup.onmousedown = startLegendDrag;

    svg.onmousemove = onLegendDrag;
    svg.onmouseup = endLegendDrag;
    svg.onmouseleave = endLegendDrag;
  };

  return {
    endLegendDrag,
    onLegendDrag,
    resetLegendPosition,
    resetLegendPositionOnly,
    setupLegendDrag,
    startLegendDrag
  };
};
