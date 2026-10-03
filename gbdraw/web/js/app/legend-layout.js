import { captureDecorationContinuity } from './legend-layout/decoration-continuity.js';
import { createLegendCanvasActions } from './legend-layout/canvas-actions.js';
import { createDiagramDragActions } from './legend-layout/diagram-drag.js';
import { createLegendRepositionActions } from './legend-layout/reposition-actions.js';
import {
  applyCompositionUserDeltas,
  resetCompositionUserDeltas
} from './legend-layout/composition-actions.js';

export const createLegendLayout = ({
  state,
  legendActions,
  history = null,
  previewRuntime = null,
  similarityAlignmentLifecycle = null
}) => {
  const diagramActions = createDiagramDragActions({
    state,
    history,
    previewRuntime,
    similarityAlignmentLifecycle
  });
  const canvasActions = createLegendCanvasActions({ state, previewRuntime });
  const repositionActions = createLegendRepositionActions({
    state,
    legendActions,
    previewRuntime
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

  const reconcileCompositionUserDeltas = (deltas) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const svg = state.svgContainer.value?.querySelector?.('svg') || null;
    if (!svg || !deltas) return false;
    const { binding, changed } = applyCompositionUserDeltas(svg, deltas);
    if (!changed) return false;
    repositionActions.syncStateFromComposition(svg, binding);
    previewRuntime?.commitActiveResultEdit('layout-composition-reconcile');
    return true;
  };

  return {
    captureDecorationContinuity: (canonical, projectRecordIdentity) => captureDecorationContinuity({
      canonical, projectRecordIdentity, results: state.results.value, catalog: state.featureCatalog.value,
      mountedSvg: state.svgContainer.value?.querySelector?.('svg') || null,
      selectedResultIndex: state.selectedResultIndex.value,
      canvasPadding: state.canvasPadding
    }),
    ...canvasActions,
    ...repositionActions,
    refreshDiagramDragAffordances: diagramActions.refreshDiagramDragAffordances,
    reconcileCompositionUserDeltas,
    resetAllPositions,
    setupDiagramDrag: diagramActions.setupDiagramDrag
  };
};
