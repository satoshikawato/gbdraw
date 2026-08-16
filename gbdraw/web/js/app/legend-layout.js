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
    return canvasActions.runCanvasAction('composition-position-reset', () => (
      canvasActions.persistCurrentSvg(svg, ({ svg: targetSvg, mutation }) => {
        diagramActions.resetLengthBarPosition(mutation);
        const binding = resetCompositionUserDeltas(targetSvg, { mutation });
        repositionActions.syncStateFromComposition(targetSvg, binding);
        return true;
      }, 'composition-position-reset')
    ));
  };

  const reconcileCompositionUserDeltas = (deltas, { lease = null } = {}) => {
    const svg = state.svgContainer.value?.querySelector?.('svg') || null;
    if (!svg || !deltas) return false;
    return canvasActions.persistCurrentSvg(svg, ({ svg: targetSvg, mutation }) => {
      const { binding, changed } = applyCompositionUserDeltas(targetSvg, deltas, { mutation });
      if (!changed) return false;
      repositionActions.syncStateFromComposition(targetSvg, binding);
      return true;
    }, 'composition-position-replay', { lease });
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
