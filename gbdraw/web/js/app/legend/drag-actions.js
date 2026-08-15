import { parseTransform } from './utils.js';
import {
  bindCompositionMetadata,
  COMPOSITION_SCHEMA_ATTRIBUTE,
  compositionUserDeltas
} from '../legend-layout/composition-actions.js';
import { replaceLeadingTranslate } from '../legend-layout/transform-utils.js';

export const createLegendDragActions = ({
  state,
  extractLegendEntries,
  history = null,
  previewRuntime
}) => {
  if (typeof previewRuntime?.commitDomEdit !== 'function') {
    throw new Error('createLegendDragActions requires the preview runtime edit protocol.');
  }
  const {
    svgContainer,
    legendDragging,
    legendDragStart,
    legendOriginalTransform,
    legendInitialTransform,
    legendCurrentOffset,
    layoutRepositionMode,
    zoom
  } = state;
  let legendDragFrameId = null;
  let pendingLegendPointer = null;
  let legendDragTxPromise = null;
  let legendDragContext = null;

  const isLayoutRepositionModeEnabled = () => Boolean(layoutRepositionMode?.value);

  const setElementCursor = (element, cursor) => {
    if (!element?.style) return;
    if (cursor) {
      element.style.cursor = cursor;
    } else {
      element.style.removeProperty('cursor');
    }
  };

  const cancelLegendDragFrame = () => {
    if (legendDragFrameId !== null) {
      cancelAnimationFrame(legendDragFrameId);
      legendDragFrameId = null;
    }
  };

  const applyLegendDragPosition = (clientX, clientY) => {
    if (!legendDragging.value) return;
    const legendGroup = legendDragContext?.binding.legend.targets[0] || null;
    if (!legendGroup) return;

    const deltaX = (clientX - legendDragStart.x) / zoom.value;
    const deltaY = (clientY - legendDragStart.y) / zoom.value;

    const newX = legendOriginalTransform.value.x + deltaX;
    const newY = legendOriginalTransform.value.y + deltaY;

    legendGroup.setAttribute(
      'transform',
      replaceLeadingTranslate(legendGroup.getAttribute('transform'), newX, newY)
    );
    legendCurrentOffset.x = newX - legendInitialTransform.value.x;
    legendCurrentOffset.y = newY - legendInitialTransform.value.y;
  };

  const startLegendDrag = (e) => {
    if (!isLayoutRepositionModeEnabled()) return;
    if (e.shiftKey) return;
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    const binding = bindCompositionMetadata(svg);
    const legendGroup = binding.legend.targets[0] || null;
    if (!legendGroup) return;

    e.preventDefault();
    e.stopPropagation();

    cancelLegendDragFrame();
    pendingLegendPointer = null;
    legendDragContext = { binding, svg };
    legendDragTxPromise = history?.begin
      ? history.begin('Move legend', { source: 'legend-drag' })
      : null;
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

  const endLegendDrag = async (e) => {
    if (!legendDragging.value) return;
    const finalPointer =
      typeof e?.clientX === 'number' && typeof e?.clientY === 'number'
        ? { x: e.clientX, y: e.clientY }
        : pendingLegendPointer;
    cancelLegendDragFrame();
    if (finalPointer) {
      applyLegendDragPosition(finalPointer.x, finalPointer.y);
    }

    const completedDragContext = legendDragContext;
    const completedLegendGroup = completedDragContext?.binding.legend.targets[0] || null;
    if (completedLegendGroup) {
      completedLegendGroup.style.willChange = '';
    }

    pendingLegendPointer = null;
    legendDragging.value = false;
    legendDragContext = null;

    if (completedDragContext?.svg) {
      const svg = completedDragContext.svg;
      previewRuntime.commitDomEdit({
        reason: 'legend-position',
        invalidateIndexes: ['legend'],
        mutate: () => Boolean(svg)
      });
    }

    const tx = legendDragTxPromise ? await legendDragTxPromise : null;
    legendDragTxPromise = null;
    if (tx && history?.commit) await history.commit(tx);
  };

  const refreshLegendDragAffordances = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    if (svg.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) !== '1') return;
    const legendGroup = bindCompositionMetadata(svg).legend.targets[0] || null;
    if (!legendGroup) return;

    setElementCursor(legendGroup, isLayoutRepositionModeEnabled() ? 'grab' : '');
  };

  const resetLegendPositionOnly = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    const binding = bindCompositionMetadata(svg);
    const legendGroup = binding.legend.targets[0] || null;
    if (!legendGroup) return;

    const automatic = binding.metadata.legend?.automaticTranslation || [0, 0];
    const initial = { x: automatic[0], y: automatic[1] };
    legendInitialTransform.value = initial;
    legendGroup.setAttribute(
      'transform',
      replaceLeadingTranslate(legendGroup.getAttribute('transform'), initial.x, initial.y)
    );
    legendCurrentOffset.x = 0;
    legendCurrentOffset.y = 0;

    previewRuntime.commitDomEdit({
      reason: 'legend-position-reset',
      invalidateIndexes: ['legend'],
      mutate: () => true
    });
  };

  const resetLegendPosition = () => {
    resetLegendPositionOnly();
    extractLegendEntries();
  };

  const setupLegendDrag = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    if (svg.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) !== '1') return;
    const binding = bindCompositionMetadata(svg);
    const legendGroup = binding.legend.targets[0] || null;
    if (!legendGroup) return;

    const automatic = binding.metadata.legend.automaticTranslation;
    const offsets = compositionUserDeltas(svg).legend || [0, 0];
    legendInitialTransform.value = { x: automatic[0], y: automatic[1] };
    legendCurrentOffset.x = offsets[0];
    legendCurrentOffset.y = offsets[1];

    legendGroup.onmousedown = startLegendDrag;
    refreshLegendDragAffordances();

    svg.onmousemove = onLegendDrag;
    svg.onmouseup = endLegendDrag;
    svg.onmouseleave = endLegendDrag;
  };

  return {
    endLegendDrag,
    onLegendDrag,
    refreshLegendDragAffordances,
    resetLegendPosition,
    resetLegendPositionOnly,
    setupLegendDrag,
    startLegendDrag
  };
};
