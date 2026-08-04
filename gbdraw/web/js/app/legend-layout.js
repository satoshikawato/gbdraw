import { createLegendCanvasActions } from './legend-layout/canvas-actions.js';
import { createDiagramDragActions } from './legend-layout/diagram-drag.js';
import { createLegendRepositionActions } from './legend-layout/reposition-actions.js';
import { resetCompositionUserDeltas } from './legend-layout/composition-actions.js';

export const createLegendLayout = ({ state, debugLog, legendActions, svgActions, history = null }) => {
  const diagramActions = createDiagramDragActions({ state, debugLog, history });
  const canvasActions = createLegendCanvasActions({ state });
  const repositionActions = createLegendRepositionActions({
    state,
    debugLog,
    legendActions,
    svgActions,
    diagramActions
  });

  const resetAllPositions = () => {
    const svg = state.svgContainer.value?.querySelector?.('svg') || null;
    if (!svg) return;
    const binding = resetCompositionUserDeltas(svg);
    repositionActions.syncStateFromComposition(svg, binding);
    diagramActions.resetLengthBarPosition();
    canvasActions.persistCurrentSvg(svg);
  };

  return {
    ...canvasActions,
    ...repositionActions,
    refreshDiagramDragAffordances: diagramActions.refreshDiagramDragAffordances,
    resetAllPositions,
    setupDiagramDrag: diagramActions.setupDiagramDrag
  };
};
