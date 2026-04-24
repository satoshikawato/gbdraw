import { createLegendCanvasActions } from './legend-layout/canvas-actions.js';
import { createDiagramDragActions } from './legend-layout/diagram-drag.js';
import { createLegendRepositionActions } from './legend-layout/reposition-actions.js';

export const createLegendLayout = ({
  state,
  debugLog,
  legendActions,
  svgActions,
  onDiagramDragCommitted = null
}) => {
  const diagramActions = createDiagramDragActions({ state, onDiagramDragCommitted });
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
    diagramActions.resetPlotTitlePosition();
    legendActions.resetLegendPositionOnly();
  };

  return {
    ...canvasActions,
    ...repositionActions,
    clearPlotTitleState: diagramActions.clearPlotTitleState,
    resetAllPositions,
    resetPlotTitlePosition: diagramActions.resetPlotTitlePosition,
    setPlotTitleAutoTransform: diagramActions.setPlotTitleAutoTransform,
    setupDiagramDrag: diagramActions.setupDiagramDrag
  };
};
