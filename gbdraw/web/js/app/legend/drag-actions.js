// @ts-check
import { parseTransform } from '../../services/legend-svg.js';
import { setClassToken } from '../../services/svg-serialization.js';
import {
  bindCompositionMetadata,
  COMPOSITION_SCHEMA_ATTRIBUTE,
  compositionUserDeltas
} from '../legend-layout/composition-actions.js';
import { replaceLeadingTranslate } from '../../services/svg-transform.js';

/**
 * @typedef {object} LegendDragActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(options?: { replaceGeneratedInventory?: boolean }) => void} extractLegendEntries
 *   The Legend entry owner's re-read of the mounted Legend.
 * @property {((label?: string, options?: { source?: string, owner?: unknown }) => Promise<any>) | null} [beginHistoryTransaction]
 *   History's begin of one step (R11); resolves to the transaction, or null when History is busy.
 * @property {((transaction: any, options?: Record<string, any>) => Promise<any>) | null} [commitHistoryTransaction]
 *   History's commit of the step that `beginHistoryTransaction` opened.
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 */

/** @param {LegendDragActionsOptions} options */
export const createLegendDragActions = ({
  state,
  extractLegendEntries,
  beginHistoryTransaction = null,
  commitHistoryTransaction = null,
  commitActiveResultEdit = null
}) => {
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
  /** @type {number | null} */
  let legendDragFrameId = null;
  /** @type {{ x: number, y: number } | null} */
  let pendingLegendPointer = null;
  /** @type {Promise<any> | null} */
  let legendDragTxPromise = null;
  /** @type {{ binding: ReturnType<typeof bindCompositionMetadata>, svg: any } | null} */
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
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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
    // The drag gesture owns its transaction and settles a focused control's (N-18).
    legendDragTxPromise = beginHistoryTransaction
      ? beginHistoryTransaction('Move legend', { source: 'legend-drag', owner: Symbol('Move legend') })
      : null;
    legendDragging.value = true;
    legendDragStart.x = e.clientX;
    legendDragStart.y = e.clientY;

    const currentTransform = parseTransform(legendGroup.getAttribute('transform'));
    legendOriginalTransform.value = { ...currentTransform };
    legendGroup.style.willChange = 'transform';
    setElementCursor(legendGroup, 'grabbing');
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
      setElementCursor(completedLegendGroup, isLayoutRepositionModeEnabled() ? 'grab' : 'help');
    }

    pendingLegendPointer = null;
    legendDragging.value = false;
    legendDragContext = null;

    if (completedDragContext?.svg) commitActiveResultEdit?.('legend-drag');

    const tx = legendDragTxPromise ? await legendDragTxPromise : null;
    legendDragTxPromise = null;
    if (tx && commitHistoryTransaction) await commitHistoryTransaction(tx);
  };

  const refreshLegendDragAffordances = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    if (svg.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) !== '1') return;
    const legendGroup = bindCompositionMetadata(svg).legend.targets[0] || null;
    if (!legendGroup) return;

    setClassToken(legendGroup, 'gbdraw-preview-layout-target', true);
    setElementCursor(legendGroup, legendDragging.value ? 'grabbing' : isLayoutRepositionModeEnabled() ? 'grab' : 'help');
  };

  const resetLegendPositionOnly = () => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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

    commitActiveResultEdit?.('legend-position-reset');
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
    // bindRole fails on a legend target without metadata, so the metadata check never differs.
    if (!legendGroup || !binding.metadata.legend) return;

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
