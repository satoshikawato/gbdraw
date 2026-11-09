// @ts-check
import { applyCanvasToSvg, bindCompositionMetadata } from './composition-actions.js';
import { CANVAS_SLICE, setCanvasPaddingSide as setDrawingCanvasPaddingSide } from '../../services/drawing/canvas.js';
/** @import { WebState } from '../../state.js' */

/**
 * @typedef {object} LegendCanvasActionsOptions
 * @property {Pick<WebState, 'activeDrawing' | 'sessionOperationAvailability' | 'svgContainer'>} state
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 */

/** @param {LegendCanvasActionsOptions} options */
export const createLegendCanvasActions = ({ state, commitActiveResultEdit = null }) => {
  const { svgContainer } = state;

  const currentSvg = () => svgContainer.value?.querySelector?.('svg') || null;

  // One canvas padding applies to every Result (D-09, PD-OI-064): Generate
  // applies it before the candidate is published, and a displayed Result
  // receives the current value. Applying the same padding twice is a no-op.
  // Showing a Result or editing the padding places no item, so the Result
  // keeps the canvas it was saved with (PD-OI-052).
  const applyCanvasPadding = () => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    const svg = currentSvg();
    if (!svg) return false;
    const binding = bindCompositionMetadata(svg);
    if (!applyCanvasToSvg(svg, drawing.canvasPadding, { binding, refit: false })) return false;
    return Boolean(commitActiveResultEdit?.('canvas-padding'));
  };

  const resetCanvasPadding = () => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    CANVAS_SLICE.reset(drawing.slices.canvas);
    return applyCanvasPadding();
  };

  // A padding field's edit applies to the displayed Result at once.
  /**
   * @param {'top' | 'right' | 'bottom' | 'left'} side
   * @param {unknown} value
   */
  const setCanvasPaddingSide = (side, value) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    setDrawingCanvasPaddingSide(state.activeDrawing().slices.canvas, side, value);
    return applyCanvasPadding();
  };

  return {
    applyCanvasPadding,
    resetCanvasPadding,
    setCanvasPaddingSide
  };
};
