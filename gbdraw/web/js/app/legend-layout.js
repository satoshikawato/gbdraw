// @ts-check
import { captureDecorationContinuity } from './legend-layout/decoration-continuity.js';
import { createLegendCanvasActions } from './legend-layout/canvas-actions.js';
import { createDiagramDragActions } from './legend-layout/diagram-drag.js';
import { createLegendRepositionActions } from './legend-layout/reposition-actions.js';
import {
  applyCompositionUserDeltas,
  compositionUserDeltas,
  resetCompositionUserDeltas
} from './legend-layout/composition-actions.js';

/**
 * The preview owner's ports that this root reads (R13); only this root holds
 * the preview owner.
 * @typedef {object} LegendLayoutPreviewPort
 * @property {(reason: string) => boolean} commitActiveResultEdit
 * @property {(result: any) => string} getResultIdentity
 * @property {() => { svg?: Element | null, resultIdentity?: string } | null} getActiveRuntime
 * @property {(
 *   resultIndex: number,
 *   edit: (svg: SVGSVGElement, context: { mounted: boolean }) => boolean,
 *   reason?: string
 * ) => boolean} commitResultEdit
 */

/**
 * @typedef {object} LegendLayoutOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(svg: SVGSVGElement, options?: { side?: string }) => import('../services/legend-layout.js').LayoutBox | null} layOutLegend
 *   The Legend manager's layout of the Legend for a side as Python lays it out.
 * @property {((label?: string, options?: { source?: string, owner?: unknown }) => Promise<any>) | null} [beginHistoryTransaction]
 *   History's begin of one step (R11); resolves to the transaction, or null when History is busy.
 * @property {((transaction: any, options?: Record<string, any>) => Promise<any>) | null} [commitHistoryTransaction]
 *   History's commit of the step that `beginHistoryTransaction` opened.
 * @property {LegendLayoutPreviewPort | null} [previewRuntime]
 * @property {{ beforeRecordDrag?: () => void, afterRecordDrag?: (options: { moved: boolean }) => void } | null} [similarityAlignmentLifecycle]
 *   The alignment owner's record-drag hooks, registered by the root (R13).
 */

/** @param {LegendLayoutOptions} options */
export const createLegendLayout = ({
  state,
  layOutLegend,
  beginHistoryTransaction = null,
  commitHistoryTransaction = null,
  previewRuntime = null,
  similarityAlignmentLifecycle = null
}) => {
  // R13: the drag, canvas, and position owners commit their edit through the
  // preview owner's port; only this root holds the preview owner.
  const commitActiveResultEdit = previewRuntime?.commitActiveResultEdit;
  const diagramActions = createDiagramDragActions({
    state,
    beginHistoryTransaction,
    commitHistoryTransaction,
    commitActiveResultEdit,
    similarityAlignmentLifecycle
  });
  const canvasActions = createLegendCanvasActions({ state, commitActiveResultEdit });
  const repositionActions = createLegendRepositionActions({
    state,
    layOutLegend,
    commitActiveResultEdit
  });

  const resetAllPositions = () => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const svg = state.svgContainer.value?.querySelector?.('svg') || null;
    if (!svg) return;
    diagramActions.resetLengthBarPosition();
    const binding = resetCompositionUserDeltas(svg);
    repositionActions.syncStateFromComposition(svg, binding);
    previewRuntime?.commitActiveResultEdit('layout-position-reset');
  };

  // B17 (R11, D-07): History records composition offsets per Result, keyed by
  // the Result's committed identity. A Result appears in the record only while
  // its offsets differ from the offsets History first saw it with, so showing
  // a Result records nothing, and Undo and Redo restore the Result a step was
  // made on. Only the displayed Result is read, and only when the mounted root
  // is bound to it; another Result keeps the offsets last seen or restored.
  const firstSeenDeltas = new Map();
  const movedDeltas = new Map();
  const recordSeenDeltas = (identity, deltas) => {
    if (!firstSeenDeltas.has(identity)) firstSeenDeltas.set(identity, deltas);
    if (JSON.stringify(deltas) === JSON.stringify(firstSeenDeltas.get(identity))) {
      movedDeltas.delete(identity);
    } else {
      movedDeltas.set(identity, deltas);
    }
  };
  const resultIdentity = (result) => previewRuntime?.getResultIdentity?.(result) || '';
  const captureCompositionIntent = () => {
    const runtime = previewRuntime?.getActiveRuntime?.() || null;
    // state is untyped; the container holds our own Result, so this is an SVG or null.
    const svg = /** @type {SVGSVGElement | null} */ (state.svgContainer.value?.querySelector?.('svg') || null);
    if (svg && runtime?.svg === svg && runtime.resultIdentity) {
      /** @type {Record<string, any> | null} */
      let deltas = null;
      try {
        deltas = compositionUserDeltas(svg);
      } catch (_error) {
        deltas = null;
      }
      if (deltas) recordSeenDeltas(runtime.resultIdentity, deltas);
    }
    const record = {};
    if (movedDeltas.size === 0) return record;
    state.results.value.forEach((result) => {
      const identity = resultIdentity(result);
      if (movedDeltas.has(identity)) record[identity] = movedDeltas.get(identity);
    });
    return record;
  };

  // Restores the recorded offsets of the named Results; a Result absent from
  // the record returns to the offsets History first saw it with.
  const reconcileCompositionUserDeltas = (record, identities) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    let reconciled = false;
    identities.forEach((identity) => {
      const deltas = record?.[identity] || firstSeenDeltas.get(identity);
      const index = state.results.value.findIndex((result) => resultIdentity(result) === identity);
      if (!deltas || index < 0) return;
      const changed = previewRuntime?.commitResultEdit(index, (svg, { mounted }) => {
        const applied = applyCompositionUserDeltas(svg, deltas);
        if (applied.changed && mounted) repositionActions.syncStateFromComposition(svg, applied.binding);
        return applied.changed;
      }, 'layout-composition-reconcile');
      recordSeenDeltas(identity, deltas);
      reconciled = reconciled || Boolean(changed);
    });
    return reconciled;
  };

  return {
    /** @param {Record<string, any>} drawing The drawing of the run that keeps the decorations. */
    captureDecorationContinuity: (canonical, projectRecordIdentity, drawing) => captureDecorationContinuity({
      canonical, projectRecordIdentity, results: state.results.value, catalog: state.featureCatalog.value,
      mountedSvg: state.svgContainer.value?.querySelector?.('svg') || null,
      selectedResultIndex: state.selectedResultIndex.value,
      canvasPadding: drawing.canvasPadding
    }),
    ...canvasActions,
    ...repositionActions,
    refreshDiagramDragAffordances: diagramActions.refreshDiagramDragAffordances,
    captureCompositionIntent,
    reconcileCompositionUserDeltas,
    resetAllPositions,
    setupDiagramDrag: diagramActions.setupDiagramDrag
  };
};
