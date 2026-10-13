// @ts-check
import { parseTransform, releasePressedTransform, replaceLeadingTranslate } from '../../services/svg-transform.js';
import {
  bindCompositionMetadata,
  COMPOSITION_SCHEMA_ATTRIBUTE,
  compositionUserDeltas
} from './composition-actions.js';
import { closestRecordGroup, isMultiRecordCanvasSvg } from '../record-groups.js';
import { setClassToken } from '../../services/svg-serialization.js';

/**
 * @typedef {object} DiagramDragActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {((label?: string, options?: { source?: string, owner?: unknown }) => Promise<any>) | null} [beginHistoryTransaction]
 *   History's begin of one step (R11); resolves to the transaction, or null when History is busy.
 * @property {((transaction: any, options?: Record<string, any>) => Promise<any>) | null} [commitHistoryTransaction]
 *   History's commit of the step that `beginHistoryTransaction` opened.
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 * @property {LegendLayoutRecordDragPort | null} [similarityAlignmentLifecycle]
 *   The alignment owner's record-drag hooks, registered by the root (R13).
 */

/**
 * @typedef {object} LegendLayoutRecordDragPort
 * @property {() => void} [beforeRecordDrag]
 * @property {(options: { moved: boolean }) => void} [afterRecordDrag]
 */

/** @param {DiagramDragActionsOptions} options */
export const createDiagramDragActions = ({
  state,
  beginHistoryTransaction = null,
  commitHistoryTransaction = null,
  commitActiveResultEdit = null,
  similarityAlignmentLifecycle = null
}) => {
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
  /** @param {Element | null | undefined} group */
  const isLegendOwnedGroup = (group) =>
    Boolean(
      group &&
      (
        LEGEND_GROUP_IDS.has(String(group.id || '')) ||
        group.getAttribute?.('data-gbdraw-role') === 'comparison-legend' ||
        group.closest?.('#legend')
      )
    );

  /** @type {SVGGraphicsElement[]} */
  let activeDragElements = [];
  let activeDragMode = 'group'; // 'group' | 'record' | 'length_bar' | 'plot_title'
  /** @type {Map<SVGGraphicsElement, { x: number, y: number }>} */
  let activeDragOriginalTransforms = new Map();
  /**
   * What each pressed target was before the press. A target whose release
   * never arrived keeps its first record.
   * @type {Map<SVGGraphicsElement, { transform: string | null, opacity: string }>}
   */
  let pressedDragTargets = new Map();
  let activeLengthBarOffsetStart = { x: 0, y: 0 };
  let activePlotTitleOffsetStart = { x: 0, y: 0 };
  /** @type {number | null} */
  let diagramDragFrameId = null;
  /** @type {{ x: number, y: number } | null} */
  let pendingDiagramPointer = null;
  /** @type {Promise<any> | null} */
  let diagramDragTxPromise = null;

  const isLayoutRepositionModeEnabled = () => Boolean(layoutRepositionMode?.value);

  /**
   * @param {{ style?: CSSStyleDeclaration } | null | undefined} element
   * @param {string | null | undefined} cursor
   */
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

  /** @param {string} label */
  const beginDragTransaction = (label) => {
    // Each drag gesture owns its transaction and settles a focused control's (N-18).
    diagramDragTxPromise = beginHistoryTransaction
      ? beginHistoryTransaction(label, { source: 'diagram-drag', owner: Symbol(label) })
      : null;
  };

  /** @param {{ id?: string } | null | undefined} group */
  const isLengthBarGroup = (group) => (group?.id || '') === 'length_bar';
  /** @param {{ id?: string } | null | undefined} group */
  const isPlotTitleGroup = (group) => (group?.id || '') === 'plot_title';

  /**
   * @param {Element | null | undefined} el
   * @param {number} x
   * @param {number} y
   */
  const setTranslate = (el, x, y) => {
    if (!el) return;
    el.setAttribute(
      'transform',
      replaceLeadingTranslate(el.getAttribute('transform'), x, y)
    );
  };

  /** @param {{ x?: unknown, y?: unknown } | null | undefined} transform */
  const normalizeTransform = (transform) => {
    const x = Number(transform?.x);
    const y = Number(transform?.y);
    return {
      x: Number.isFinite(x) ? x : 0,
      y: Number.isFinite(y) ? y : 0
    };
  };

  /** @param {SVGGraphicsElement | null | undefined} group */
  const assignPlotTitleElement = (group) => {
    if (plotTitleElement.value && plotTitleElement.value !== group) {
      setElementCursor(plotTitleElement.value, '');
    }
    plotTitleElement.value = group || null;
    if (plotTitleElement.value) {
      setElementCursor(plotTitleElement.value, isLayoutRepositionModeEnabled() ? 'grab' : '');
    }
  };

  const applyPlotTitleTransform = (group = plotTitleElement.value) => {
    if (!group) return;
    const autoTransform = normalizeTransform(plotTitleAutoTransform.value);
    const nextX = autoTransform.x + plotTitleUserOffset.x;
    const nextY = autoTransform.y + plotTitleUserOffset.y;
    setTranslate(group, nextX, nextY);
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

  const applyLengthBarTransform = () => {
    const group = lengthBarElement.value;
    if (!group) return;
    const base = normalizeTransform(lengthBarOriginalTransform.value);
    const nextX = base.x + lengthBarUserOffset.x;
    const nextY = base.y + lengthBarUserOffset.y;
    setTranslate(group, nextX, nextY);
  };

  const resetLengthBarPosition = () => {
    lengthBarUserOffset.x = 0;
    lengthBarUserOffset.y = 0;
    if (lengthBarElement.value) applyLengthBarTransform();
  };

  const clearLengthBarState = () => {
    if (lengthBarElement.value) setElementCursor(lengthBarElement.value, '');
    lengthBarElement.value = null;
    lengthBarOriginalTransform.value = { x: 0, y: 0 };
    lengthBarUserOffset.x = 0;
    lengthBarUserOffset.y = 0;
    activeLengthBarOffsetStart = { x: 0, y: 0 };
  };

  /**
   * @param {SVGSVGElement | null | undefined} svg
   * @param {boolean} [preserveOffset]
   */
  const syncLengthBarElement = (svg, preserveOffset = false) => {
    // The Python-drawn scale bar is a `g` group.
    const nextLengthBar = /** @type {SVGGraphicsElement | null} */ (svg?.getElementById('length_bar') || null);
    if (!nextLengthBar) {
      clearLengthBarState();
      return;
    }

    const sameElement = lengthBarElement.value === nextLengthBar;
    lengthBarElement.value = nextLengthBar;
    setElementCursor(lengthBarElement.value, isLayoutRepositionModeEnabled() ? 'grab' : '');

    if (!preserveOffset) {
      lengthBarUserOffset.x = 0;
      lengthBarUserOffset.y = 0;
    }
    if (!preserveOffset || !sameElement) {
      // A replacement SVG already contains its own scale translation. Binding
      // adopts that geometry; only a later drag applies an additional offset.
      const current = parseTransform(nextLengthBar.getAttribute('transform'));
      lengthBarOriginalTransform.value = {
        x: current.x - lengthBarUserOffset.x,
        y: current.y - lengthBarUserOffset.y
      };
    }
    applyLengthBarTransform();
  };

  /** @param {Element | null | undefined} svg */
  const syncPlotTitleElement = (svg) => {
    if (svg?.getAttribute?.(COMPOSITION_SCHEMA_ATTRIBUTE) !== '1') {
      clearPlotTitleState();
      return;
    }
    const binding = bindCompositionMetadata(svg);
    const nextPlotTitleGroup = binding.title.targets[0] || null;
    if (!nextPlotTitleGroup || !binding.metadata.title) {
      // bindRole fails on a title target without metadata, so the metadata check never differs.
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

  // Whether a press drags one record (a Multi-Record Canvas) or the whole diagram.
  const dragsOneRecord = () => isMultiRecordCanvasSvg(svgContainer.value?.querySelector('svg'));

  // What a press on `target` drags: its record on a Multi-Record Canvas, else
  // every bound primary group; null when the press drags nothing.
  /**
   * @param {Element} target
   * @param {boolean} [oneRecord] `dragsOneRecord()`, read once by a caller that asks for many targets.
   * @returns {{ mode: 'record' | 'group', targets: SVGGraphicsElement[] } | null}
   */
  const diagramDragOf = (target, oneRecord = dragsOneRecord()) => {
    if (oneRecord) {
      const record = closestRecordGroup(target);
      // `closestRecordGroup` returns a `g` element.
      return record ? { mode: 'record', targets: [/** @type {SVGGraphicsElement} */ (record)] } : null;
    }
    const group = target.closest('g[id]');
    if (!group || isLegendOwnedGroup(group) || isLengthBarGroup(group) || isPlotTitleGroup(group)) return null;
    // setupDiagramDrag stores the bound primary groups (or []) in `diagramElements`.
    /** @type {SVGGraphicsElement[]} */
    const targets = diagramElements.value;
    return targets.length > 0 ? { mode: 'group', targets } : null;
  };

  // A primary group shows `grab` only when a press on it drags it; one that no
  // press drags (a Linear record definition) shows no drag cursor (PD-OI-053).
  const refreshDiagramDragAffordances = () => {
    const enabled = isLayoutRepositionModeEnabled();
    /** @type {SVGGraphicsElement[]} */
    const elements = Array.isArray(diagramElements.value) ? diagramElements.value : [];
    const oneRecord = enabled && dragsOneRecord();

    elements.forEach((el) => {
      setElementCursor(el, !enabled ? '' : diagramDragOf(el, oneRecord)?.targets.includes(el) ? 'grab' : 'default');
    });

    [lengthBarElement.value, plotTitleElement.value].forEach((element) => {
      setClassToken(element, 'gbdraw-preview-layout-target', true);
      setElementCursor(element, activeDragElements.includes(element) && diagramDragging.value
        ? 'grabbing' : enabled ? 'grab' : 'help');
    });
  };

  // The press of every drag: its targets show the drag and keep what they
  // were, so the release can put it back (`endDiagramDrag`).
  /** @param {SVGGraphicsElement[]} targets */
  const pressDragTargets = (targets) => {
    activeDragElements = targets;
    activeDragOriginalTransforms = new Map();
    targets.forEach((el) => {
      const transform = el.getAttribute('transform');
      activeDragOriginalTransforms.set(el, parseTransform(transform));
      if (!pressedDragTargets.has(el)) pressedDragTargets.set(el, { transform, opacity: el.style.opacity });
      el.style.opacity = '0.8';
      el.style.willChange = 'transform';
      setElementCursor(el, 'grabbing');
    });
    document.addEventListener('mousemove', onDiagramDrag);
    document.addEventListener('mouseup', endDiagramDrag);
  };

  /**
   * @param {MouseEvent} e
   * @param {SVGGraphicsElement | null | undefined} group
   */
  const startLengthBarDrag = (e, group) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!isLayoutRepositionModeEnabled()) return;
    if (!group) return;

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
    pressDragTargets([group]);
  };

  /**
   * @param {MouseEvent} e
   * @param {SVGGraphicsElement | null | undefined} group
   */
  const startPlotTitleDrag = (e, group) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!isLayoutRepositionModeEnabled()) return;
    if (!group) return;

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
    pressDragTargets([group]);
  };

  /** @param {MouseEvent & { target: Element }} e A mouse event on the mounted SVG, whose targets are elements. */
  const startDiagramDrag = (e) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!isLayoutRepositionModeEnabled()) return;
    if (e.shiftKey) return;
    if (
      e.target.closest('#legend') ||
      e.target.closest('text[data-label-editable="true"]') ||
      e.target.closest('path[id^="f"], polygon[id^="f"], rect[id^="f"]')
    ) {
      return;
    }
    if (!svgContainer.value?.querySelector('svg')) return;

    const clickedLengthBar = e.target.closest('#length_bar');
    if (clickedLengthBar) {
      startLengthBarDrag(e, /** @type {SVGGraphicsElement} */ (clickedLengthBar)); // the `g#length_bar` group
      return;
    }

    const clickedPlotTitle = e.target.closest('#plot_title');
    if (clickedPlotTitle) {
      startPlotTitleDrag(e, /** @type {SVGGraphicsElement} */ (clickedPlotTitle)); // the `g#plot_title` group
      return;
    }

    const drag = diagramDragOf(e.target);
    if (!drag) return;
    activeDragMode = drag.mode;

    e.preventDefault();
    cancelDiagramDragFrame();
    pendingDiagramPointer = null;
    beginDragTransaction(activeDragMode === 'record' ? 'Move record' : 'Move diagram');
    if (activeDragMode === 'record') {
      similarityAlignmentLifecycle?.beforeRecordDrag?.();
    }
    diagramDragging.value = true;
    plotTitleDragging.value = false;
    diagramDragStart.x = e.clientX;
    diagramDragStart.y = e.clientY;
    pressDragTargets(drag.targets);
  };

  const getActiveDragStart = () => {
    return activeDragMode === 'plot_title' ? plotTitleDragStart : diagramDragStart;
  };

  /**
   * @param {number} clientX
   * @param {number} clientY
   */
  const applyDiagramDragPosition = (clientX, clientY) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!diagramDragging.value || activeDragElements.length === 0) return;

    const dragStart = getActiveDragStart();
    const deltaX = (clientX - dragStart.x) / zoom.value;
    const deltaY = (clientY - dragStart.y) / zoom.value;

    if (activeDragMode === 'plot_title') {
      plotTitleUserOffset.x = activePlotTitleOffsetStart.x + deltaX;
      plotTitleUserOffset.y = activePlotTitleOffsetStart.y + deltaY;
      applyPlotTitleTransform(activeDragElements[0] || plotTitleElement.value);
      return;
    }

    if (activeDragMode === 'length_bar') {
      const group = activeDragElements[0];
      const original = activeDragOriginalTransforms.get(group) || { x: 0, y: 0 };
      setTranslate(group, original.x + deltaX, original.y + deltaY);
      return;
    }

    activeDragElements.forEach((el) => {
      const original = activeDragOriginalTransforms.get(el) || { x: 0, y: 0 };
      const newX = original.x + deltaX;
      const newY = original.y + deltaY;
      setTranslate(el, newX, newY);
    });
  };

  /** @param {MouseEvent} e */
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

  /** @param {MouseEvent} [e] */
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

    diagramDragging.value = false;
    plotTitleDragging.value = false;
    document.removeEventListener('mousemove', onDiagramDrag);
    document.removeEventListener('mouseup', endDiagramDrag);

    // Only a drag that moved a target is an edit (PD-OI-053): the release puts
    // the drag feedback back, and an unmoved target's transform as written.
    const moved = [...pressedDragTargets].map(([el, pressed]) => {
      el.style.opacity = pressed.opacity;
      el.style.willChange = '';
      return releasePressedTransform(el, pressed.transform);
    }).some(Boolean);
    pressedDragTargets = new Map();
    const completedDragMode = activeDragMode;
    activeDragElements = [];
    activeDragOriginalTransforms = new Map();
    activeDragMode = 'group';
    activeLengthBarOffsetStart = { x: lengthBarUserOffset.x, y: lengthBarUserOffset.y };
    activePlotTitleOffsetStart = { x: plotTitleUserOffset.x, y: plotTitleUserOffset.y };
    pendingDiagramPointer = null;
    refreshDiagramDragAffordances();
    if (completedDragMode === 'record') {
      similarityAlignmentLifecycle?.afterRecordDrag?.({ moved });
    }
    if (moved) commitActiveResultEdit?.('diagram-drag');
    const tx = diagramDragTxPromise ? await diagramDragTxPromise : null;
    diagramDragTxPromise = null;
    if (tx && commitHistoryTransaction) await commitHistoryTransaction(tx);
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
