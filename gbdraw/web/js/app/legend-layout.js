import { createLegendCanvasActions } from './legend-layout/canvas-actions.js';
import { createDiagramDragActions } from './legend-layout/diagram-drag.js';
import { createLegendRepositionActions } from './legend-layout/reposition-actions.js';

export const createLegendLayout = ({ state, debugLog, legendActions, svgActions }) => {
  const diagramActions = createDiagramDragActions({ state });
  const canvasActions = createLegendCanvasActions({ state });
  const repositionActions = createLegendRepositionActions({
    state,
    debugLog,
    legendActions,
    svgActions,
    diagramActions
  });

  const resetAllPositions = () => {
    diagramActions.resetDiagramPosition();
    legendActions.resetLegendPositionOnly();
  };

  return {
    ...canvasActions,
    ...repositionActions,
    resetAllPositions,
    setupDiagramDrag: diagramActions.setupDiagramDrag
  };
};
