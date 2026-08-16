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
  if (
    typeof previewRuntime?.commitDomEdit !== 'function'
    || typeof previewRuntime?.beginDomEditLease !== 'function'
  ) {
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
  let legendDragLease = null;

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

    const changed = legendDragLease?.mutate?.(({ svg, mutation }) => {
      if (svg !== legendDragContext?.svg) return false;
      return mutation.setAttribute(
        legendGroup,
        'transform',
        replaceLeadingTranslate(legendGroup.getAttribute('transform'), newX, newY)
      );
    });
    if (!changed) return;
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
    const lease = previewRuntime.beginDomEditLease({
      reason: 'legend-position',
      invalidateIndexes: ['legend']
    });
    if (!lease) return;
    lease.mutate(({ mutation }) => {
      mutation.captureState(legendCurrentOffset);
      mutation.captureProperty(legendInitialTransform, 'value', { deep: true });
      return false;
    });

    e.preventDefault();
    e.stopPropagation();

    cancelLegendDragFrame();
    pendingLegendPointer = null;
    legendDragContext = { binding, svg };
    legendDragLease = lease;
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

    const completedLease = legendDragLease;
    const completedTxPromise = legendDragTxPromise;
    legendDragLease = null;
    legendDragTxPromise = null;

    try {
      if (completedDragContext?.svg) completedLease?.commit?.();
      else completedLease?.cancel?.();
      const tx = completedTxPromise ? await completedTxPromise : null;
      if (tx && history?.commit) await history.commit(tx);
    } catch (error) {
      try {
        completedLease?.cancel?.();
      } catch (rollbackError) {
        if (error instanceof Error) {
          Object.defineProperty(error, 'rollbackErrors', {
            configurable: true,
            value: Object.freeze([
              ...(Array.isArray(error.rollbackErrors) ? error.rollbackErrors : []),
              rollbackError
            ])
          });
        }
      }
      try {
        const tx = completedTxPromise ? await completedTxPromise : null;
        if (tx && history?.cancel) history.cancel(tx);
      } catch (_historyError) {
        // The live-edit failure remains primary; History has no committed entry.
      }
      throw error;
    }
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
    const committed = previewRuntime.commitDomEdit({
      reason: 'legend-position-reset',
      invalidateIndexes: ['legend'],
      mutate: ({ svg: targetSvg, mutation }) => {
        if (targetSvg !== svg) return false;
        const binding = bindCompositionMetadata(targetSvg);
        const legendGroup = binding.legend.targets[0] || null;
        if (!legendGroup) return false;
        const automatic = binding.metadata.legend?.automaticTranslation || [0, 0];
        const initial = { x: automatic[0], y: automatic[1] };
        mutation.captureProperty(legendInitialTransform, 'value', { deep: true });
        mutation.captureState(legendCurrentOffset);
        const changed = mutation.setAttribute(
          legendGroup,
          'transform',
          replaceLeadingTranslate(legendGroup.getAttribute('transform'), initial.x, initial.y)
        );
        if (changed) {
          legendInitialTransform.value = initial;
          legendCurrentOffset.x = 0;
          legendCurrentOffset.y = 0;
        }
        return changed;
      }
    });
    return committed.changed;
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
