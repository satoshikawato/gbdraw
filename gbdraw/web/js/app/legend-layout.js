import { createLegendCanvasActions } from './legend-layout/canvas-actions.js';
import { createDiagramDragActions } from './legend-layout/diagram-drag.js';
import { createLegendRepositionActions } from './legend-layout/reposition-actions.js';
import {
  applyCompositionUserDeltas,
  resetCompositionUserDeltas
} from './legend-layout/composition-actions.js';

export const createLegendLayout = ({ state, legendActions, history = null, previewRuntime }) => {
  const diagramActions = createDiagramDragActions({ state, history, previewRuntime });
  const canvasActions = createLegendCanvasActions({ state, previewRuntime });
  const repositionActions = createLegendRepositionActions({
    state,
    legendActions,
    diagramActions,
    previewRuntime
  });

  const resetAllPositions = () => {
    const svg = state.svgContainer.value?.querySelector?.('svg') || null;
    if (!svg) return;
    diagramActions.resetLengthBarPosition();
    const binding = resetCompositionUserDeltas(svg);
    repositionActions.syncStateFromComposition(svg, binding);
    canvasActions.persistCurrentSvg(svg);
  };

  const reconcileCompositionUserDeltas = (deltas) => {
    const svg = state.svgContainer.value?.querySelector?.('svg') || null;
    if (!svg || !deltas) return false;
    const { binding, changed } = applyCompositionUserDeltas(svg, deltas);
    if (!changed) return false;
    repositionActions.syncStateFromComposition(svg, binding);
    canvasActions.persistCurrentSvg(svg);
    return true;
  };

  return {
    ...canvasActions,
    ...repositionActions,
    refreshDiagramDragAffordances: diagramActions.refreshDiagramDragAffordances,
    reconcileCompositionUserDeltas,
    resetAllPositions,
    setupDiagramDrag: diagramActions.setupDiagramDrag
  };
};
