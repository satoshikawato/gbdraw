import { parseTransform, replaceLeadingTranslate } from './transform-utils.js';
import {
  bindCompositionMetadata,
  COMPOSITION_SCHEMA_ATTRIBUTE,
  compositionUserDeltas
} from './composition-actions.js';
import {
  closestRecordGroup,
  isMultiRecordCanvasSvg,
  isRecordGroup
} from '../record-groups.js';

export const createDiagramDragActions = ({ state, history = null, previewRuntime }) => {
  if (
    typeof previewRuntime?.commitDomEdit !== 'function'
    || typeof previewRuntime?.beginDomEditLease !== 'function'
  ) {
    throw new Error('createDiagramDragActions requires the preview runtime edit protocol.');
  }
  const {
    svgContainer,
    diagramElements,
    diagramElementIds,
    diagramElementOriginalTransforms,
    diagramOffset,
    diagramDragging,
    diagramDragStart,
    lengthBarElement,
    lengthBarOriginalTransform,
    lengthBarUserOffset,
    plotTitleElement,
    plotTitleDragging,
    plotTitleDragStart,
    plotTitleAutoTransform,
    plotTitleUserOffset,
    layoutRepositionMode,
    zoom
  } = state;

  const LEGEND_GROUP_IDS = new Set([
    'legend',
    'feature_legend',
    'pairwise_legend',
    'horizontal_legend',
    'vertical_legend'
  ]);
  const isLegendOwnedGroup = (group) =>
    Boolean(
      group &&
      (
        LEGEND_GROUP_IDS.has(String(group.id || '')) ||
        group.getAttribute?.('data-gbdraw-role') === 'comparison-legend' ||
        group.closest?.('#legend')
      )
    );

  let activeDragElements = [];
  let activeDragMode = 'group'; // 'group' | 'record' | 'length_bar' | 'plot_title'
  let activeDragOriginalTransforms = new Map();
  let activeLengthBarOffsetStart = { x: 0, y: 0 };
  let activePlotTitleOffsetStart = { x: 0, y: 0 };
  let diagramDragFrameId = null;
  let pendingDiagramPointer = null;
  let diagramDragTxPromise = null;
  let diagramDragLease = null;
  let activeDragSvg = null;

  const isLayoutRepositionModeEnabled = () => Boolean(layoutRepositionMode?.value);

  const setElementCursor = (element, cursor) => {
    if (!element?.style) return;
    if (cursor) {
      element.style.cursor = cursor;
    } else {
      element.style.removeProperty('cursor');
    }
  };

  const cancelDiagramDragFrame = () => {
    if (diagramDragFrameId !== null) {
      cancelAnimationFrame(diagramDragFrameId);
      diagramDragFrameId = null;
    }
  };

  const beginDragTransaction = (label) => {
    diagramDragTxPromise = history?.begin
      ? history.begin(label, { source: 'diagram-drag' })
      : null;
  };

  const beginPositionLease = (svg) => {
    const lease = previewRuntime.beginDomEditLease({
      reason: 'diagram-position',
      invalidateIndexes: ['legend']
    });
    if (!lease) return false;
    diagramDragLease?.cancel?.();
    diagramDragLease = lease;
    activeDragSvg = svg;
    return true;
  };

  const isLengthBarGroup = (group) => (group?.id || '') === 'length_bar';
  const isPlotTitleGroup = (group) => (group?.id || '') === 'plot_title';

  const setTranslate = (el, x, y, mutation = null) => {
    if (!el) return;
    const value = replaceLeadingTranslate(el.getAttribute('transform'), x, y);
    if (mutation) return mutation.setAttribute(el, 'transform', value);
    el.setAttribute('transform', value);
    return true;
  };

  const normalizeTransform = (transform) => {
    const x = Number(transform?.x);
    const y = Number(transform?.y);
    return {
      x: Number.isFinite(x) ? x : 0,
      y: Number.isFinite(y) ? y : 0
    };
  };

  const assignPlotTitleElement = (group) => {
    if (plotTitleElement.value && plotTitleElement.value !== group) {
      plotTitleElement.value.style.opacity = '1';
      setElementCursor(plotTitleElement.value, '');
    }
    plotTitleElement.value = group || null;
    if (plotTitleElement.value) {
      plotTitleElement.value.style.opacity = '1';
      setElementCursor(plotTitleElement.value, isLayoutRepositionModeEnabled() ? 'grab' : '');
    }
  };

  const applyPlotTitleTransform = (group = plotTitleElement.value, mutation = null) => {
    if (!group) return;
    const autoTransform = normalizeTransform(plotTitleAutoTransform.value);
    const nextX = autoTransform.x + plotTitleUserOffset.x;
    const nextY = autoTransform.y + plotTitleUserOffset.y;
    return setTranslate(group, nextX, nextY, mutation);
  };

  const clearPlotTitleState = () => {
    assignPlotTitleElement(null);
    plotTitleDragging.value = false;
    plotTitleDragStart.x = 0;
    plotTitleDragStart.y = 0;
    plotTitleAutoTransform.value = { x: 0, y: 0 };
    plotTitleUserOffset.x = 0;
    plotTitleUserOffset.y = 0;
    activePlotTitleOffsetStart = { x: 0, y: 0 };
  };

  const applyLengthBarTransform = (mutation = null) => {
    const group = lengthBarElement.value;
    if (!group) return;
    const base = normalizeTransform(lengthBarOriginalTransform.value);
    const nextX = base.x + lengthBarUserOffset.x;
    const nextY = base.y + lengthBarUserOffset.y;
    return setTranslate(group, nextX, nextY, mutation);
  };

  const resetLengthBarPosition = (mutation = null) => {
    lengthBarUserOffset.x = 0;
    lengthBarUserOffset.y = 0;
    if (lengthBarElement.value) {
      lengthBarElement.value.style.opacity = '1';
      applyLengthBarTransform(mutation);
    }
  };

  const clearLengthBarState = () => {
    if (lengthBarElement.value) {
      lengthBarElement.value.style.opacity = '1';
      setElementCursor(lengthBarElement.value, '');
    }
    lengthBarElement.value = null;
    lengthBarOriginalTransform.value = { x: 0, y: 0 };
    lengthBarUserOffset.x = 0;
    lengthBarUserOffset.y = 0;
    activeLengthBarOffsetStart = { x: 0, y: 0 };
  };

  const syncLengthBarElement = (svg, preserveOffset = false) => {
    const nextLengthBar = svg?.getElementById('length_bar') || null;
    if (!nextLengthBar) {
      clearLengthBarState();
      return;
    }

    const hadLengthBarState = !!lengthBarElement.value;
    lengthBarElement.value = nextLengthBar;
    lengthBarElement.value.style.opacity = '1';
    setElementCursor(lengthBarElement.value, isLayoutRepositionModeEnabled() ? 'grab' : '');

    if (!preserveOffset || !hadLengthBarState) {
      lengthBarOriginalTransform.value = parseTransform(nextLengthBar.getAttribute('transform'));
    }
    if (!preserveOffset) {
      lengthBarUserOffset.x = 0;
      lengthBarUserOffset.y = 0;
    }
    applyLengthBarTransform();
  };

  const syncPlotTitleElement = (svg) => {
    if (svg?.getAttribute?.(COMPOSITION_SCHEMA_ATTRIBUTE) !== '1') {
      clearPlotTitleState();
      return;
    }
    const binding = bindCompositionMetadata(svg);
    const nextPlotTitleGroup = binding.title.targets[0] || null;
    if (!nextPlotTitleGroup) {
      clearPlotTitleState();
      return;
    }

    const generatedTransform = {
      x: binding.metadata.title.automaticTranslation[0],
      y: binding.metadata.title.automaticTranslation[1]
    };
    const titleDelta = compositionUserDeltas(svg).title || [0, 0];
    assignPlotTitleElement(nextPlotTitleGroup);
    plotTitleAutoTransform.value = generatedTransform;
    plotTitleUserOffset.x = titleDelta[0];
    plotTitleUserOffset.y = titleDelta[1];
    applyPlotTitleTransform(nextPlotTitleGroup);
  };

  const refreshDiagramDragAffordances = () => {
    const enabled = isLayoutRepositionModeEnabled();
    const elements = Array.isArray(diagramElements.value) ? diagramElements.value : [];
    const hasRecordGroups = elements.some((el) => isRecordGroup(el));

    elements.forEach((el) => {
      if (!enabled) {
        setElementCursor(el, '');
        return;
      }
      if (hasRecordGroups) {
        setElementCursor(el, isRecordGroup(el) ? 'grab' : '');
      } else {
        setElementCursor(el, 'grab');
      }
    });

    setElementCursor(lengthBarElement.value, enabled ? 'grab' : '');
    setElementCursor(plotTitleElement.value, enabled ? 'grab' : '');
  };

  const startLengthBarDrag = (e, group, svg) => {
    if (!isLayoutRepositionModeEnabled()) return;
    if (!group) return;
    if (!beginPositionLease(svg)) return;

    e.preventDefault();
    cancelDiagramDragFrame();
    pendingDiagramPointer = null;
    beginDragTransaction('Move length bar');
    lengthBarElement.value = group;
    diagramDragging.value = true;
    plotTitleDragging.value = false;
    diagramDragStart.x = e.clientX;
    diagramDragStart.y = e.clientY;
    activeLengthBarOffsetStart = { x: lengthBarUserOffset.x, y: lengthBarUserOffset.y };
    activeDragMode = 'length_bar';
    activeDragElements = [group];
    activeDragOriginalTransforms = new Map([[group, parseTransform(group.getAttribute('transform'))]]);
    group.style.opacity = '0.8';
    group.style.willChange = 'transform';

    document.addEventListener('mousemove', onDiagramDrag);
    document.addEventListener('mouseup', endDiagramDrag);
  };

  const startPlotTitleDrag = (e, group, svg) => {
    if (!isLayoutRepositionModeEnabled()) return;
    if (!group) return;
    if (!beginPositionLease(svg)) return;

    e.preventDefault();
    cancelDiagramDragFrame();
    pendingDiagramPointer = null;
    beginDragTransaction('Move plot title');
    assignPlotTitleElement(group);
    diagramDragging.value = true;
    plotTitleDragging.value = true;
    plotTitleDragStart.x = e.clientX;
    plotTitleDragStart.y = e.clientY;
    activePlotTitleOffsetStart = { x: plotTitleUserOffset.x, y: plotTitleUserOffset.y };
    activeDragMode = 'plot_title';
    activeDragElements = [group];
    activeDragOriginalTransforms = new Map([[group, parseTransform(group.getAttribute('transform'))]]);
    group.style.opacity = '0.8';
    group.style.willChange = 'transform';

    document.addEventListener('mousemove', onDiagramDrag);
    document.addEventListener('mouseup', endDiagramDrag);
  };

  const startDiagramDrag = (e) => {
    if (!isLayoutRepositionModeEnabled()) return;
    if (e.shiftKey) return;
    if (
      e.target.closest('#legend') ||
      e.target.closest('text[data-label-editable="true"]') ||
      e.target.closest('path[id^="f"], polygon[id^="f"], rect[id^="f"]')
    ) {
      return;
    }
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const clickedLengthBar = e.target.closest('#length_bar');
    if (clickedLengthBar) {
      startLengthBarDrag(e, clickedLengthBar, svg);
      return;
    }

    const clickedPlotTitle = e.target.closest('#plot_title');
    if (clickedPlotTitle) {
      startPlotTitleDrag(e, clickedPlotTitle, svg);
      return;
    }

    const isMultiRecordCanvas = isMultiRecordCanvasSvg(svg);
    let dragTargets = [];

    if (isMultiRecordCanvas) {
      // Multi-record mode: allow dragging only the clicked record group.
      const clickedRecordGroup = closestRecordGroup(e.target);
      if (!clickedRecordGroup) return;
      dragTargets = [clickedRecordGroup];
      activeDragMode = 'record';
    } else {
      const clickedGroup = e.target.closest('g[id]');
      if (!clickedGroup) return;

      if (isLegendOwnedGroup(clickedGroup)) return;
      if (isLengthBarGroup(clickedGroup)) return;
      if (isPlotTitleGroup(clickedGroup)) return;

      dragTargets = diagramElements.value;
      if (dragTargets.length === 0) return;
      activeDragMode = 'group';
    }

    e.preventDefault();
    if (!beginPositionLease(svg)) return;
    cancelDiagramDragFrame();
    pendingDiagramPointer = null;
    beginDragTransaction(activeDragMode === 'record' ? 'Move record' : 'Move diagram');
    diagramDragging.value = true;
    plotTitleDragging.value = false;
    diagramDragStart.x = e.clientX;
    diagramDragStart.y = e.clientY;

    activeDragElements = dragTargets;
    activeDragOriginalTransforms = new Map();

    activeDragElements.forEach((el) => {
      activeDragOriginalTransforms.set(el, parseTransform(el.getAttribute('transform')));
      el.style.opacity = '0.8';
      el.style.willChange = 'transform';
    });

    document.addEventListener('mousemove', onDiagramDrag);
    document.addEventListener('mouseup', endDiagramDrag);
  };

  const getActiveDragStart = () => {
    return activeDragMode === 'plot_title' ? plotTitleDragStart : diagramDragStart;
  };

  const applyDiagramDragPosition = (clientX, clientY) => {
    if (!diagramDragging.value || activeDragElements.length === 0) return;

    const dragStart = getActiveDragStart();
    const deltaX = (clientX - dragStart.x) / zoom.value;
    const deltaY = (clientY - dragStart.y) / zoom.value;

    return diagramDragLease?.mutate?.(({ svg, mutation }) => {
      if (svg !== activeDragSvg) return false;
      if (activeDragMode === 'plot_title') {
        plotTitleUserOffset.x = activePlotTitleOffsetStart.x + deltaX;
        plotTitleUserOffset.y = activePlotTitleOffsetStart.y + deltaY;
        return applyPlotTitleTransform(activeDragElements[0] || plotTitleElement.value, mutation);
      }

      if (activeDragMode === 'length_bar') {
        const group = activeDragElements[0];
        const original = activeDragOriginalTransforms.get(group) || { x: 0, y: 0 };
        return setTranslate(group, original.x + deltaX, original.y + deltaY, mutation);
      }

      let changed = false;
      activeDragElements.forEach((el) => {
        const original = activeDragOriginalTransforms.get(el) || { x: 0, y: 0 };
        const newX = original.x + deltaX;
        const newY = original.y + deltaY;
        changed = setTranslate(el, newX, newY, mutation) || changed;
      });
      return changed;
    });
  };

  const onDiagramDrag = (e) => {
    if (!diagramDragging.value || activeDragElements.length === 0) return;
    pendingDiagramPointer = { x: e.clientX, y: e.clientY };
    if (diagramDragFrameId !== null) return;
    diagramDragFrameId = requestAnimationFrame(() => {
      diagramDragFrameId = null;
      if (!pendingDiagramPointer) return;
      applyDiagramDragPosition(pendingDiagramPointer.x, pendingDiagramPointer.y);
    });
  };

  const endDiagramDrag = async (e) => {
    if (!diagramDragging.value) return;

    const dragStart = getActiveDragStart();
    const currentX = typeof e?.clientX === 'number' ? e.clientX : dragStart.x;
    const currentY = typeof e?.clientY === 'number' ? e.clientY : dragStart.y;
    cancelDiagramDragFrame();
    applyDiagramDragPosition(currentX, currentY);
    const deltaX = (currentX - dragStart.x) / zoom.value;
    const deltaY = (currentY - dragStart.y) / zoom.value;
    if (activeDragMode === 'group') {
      const svg = svgContainer.value?.querySelector?.('svg') || null;
      const deltas = svg ? compositionUserDeltas(svg).primary : [];
      diagramOffset.x = deltas[0]?.[0] || 0;
      diagramOffset.y = deltas[0]?.[1] || 0;
    } else if (activeDragMode === 'length_bar') {
      lengthBarUserOffset.x = activeLengthBarOffsetStart.x + deltaX;
      lengthBarUserOffset.y = activeLengthBarOffsetStart.y + deltaY;
    }
    diagramDragLease?.commit?.();
    diagramDragLease = null;
    activeDragSvg = null;

    diagramDragging.value = false;
    plotTitleDragging.value = false;
    document.removeEventListener('mousemove', onDiagramDrag);
    document.removeEventListener('mouseup', endDiagramDrag);

    activeDragElements.forEach((el) => {
      el.style.opacity = '1';
      el.style.willChange = '';
    });
    activeDragElements = [];
    activeDragOriginalTransforms = new Map();
    activeDragMode = 'group';
    activeLengthBarOffsetStart = { x: lengthBarUserOffset.x, y: lengthBarUserOffset.y };
    activePlotTitleOffsetStart = { x: plotTitleUserOffset.x, y: plotTitleUserOffset.y };
    pendingDiagramPointer = null;
    const tx = diagramDragTxPromise ? await diagramDragTxPromise : null;
    diagramDragTxPromise = null;
    if (tx && history?.commit) await history.commit(tx);
  };

  const setupDiagramDrag = (preserveOffset = false) => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    if (svg.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) !== '1') {
      diagramElements.value = [];
      diagramElementIds.value = [];
      return;
    }
    const binding = bindCompositionMetadata(svg);
    if (!preserveOffset) {
      diagramOffset.x = 0;
      diagramOffset.y = 0;
    }

    const foundElements = [...binding.primary.targets];
    const foundIds = foundElements.map((element) => element.id || '').filter(Boolean);

    diagramElements.value = foundElements;
    diagramElementIds.value = foundIds;
    syncLengthBarElement(svg, preserveOffset);
    syncPlotTitleElement(svg);

    const originalTransforms = new Map();
    foundElements.forEach((el) => {
      const transform = {
        x: binding.metadata.primary.automaticTranslation[0],
        y: binding.metadata.primary.automaticTranslation[1]
      };
      originalTransforms.set(el, transform);
    });
    diagramElementOriginalTransforms.value = originalTransforms;
    const primaryDeltas = compositionUserDeltas(svg).primary;
    diagramOffset.x = primaryDeltas[0]?.[0] || 0;
    diagramOffset.y = primaryDeltas[0]?.[1] || 0;
    refreshDiagramDragAffordances();

    svg.removeEventListener('mousedown', startDiagramDrag);
    svg.addEventListener('mousedown', startDiagramDrag);
  };

  return {
    endDiagramDrag,
    onDiagramDrag,
    refreshDiagramDragAffordances,
    resetLengthBarPosition,
    setupDiagramDrag,
    startDiagramDrag
  };
};
