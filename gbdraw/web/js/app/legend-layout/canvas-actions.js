// @ts-check
import {
  applyCanvasPaddingToSvg,
  bindCompositionMetadata,
  compositionUserDeltas
} from './composition-actions.js';

/**
 * @typedef {object} LegendCanvasActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 */

/** @param {LegendCanvasActionsOptions} options */
export const createLegendCanvasActions = ({ state, commitActiveResultEdit = null }) => {
  const {
    svgContainer,
    diagramElements,
    diagramElementOriginalTransforms,
    diagramOffset,
    legendInitialTransform,
    legendCurrentOffset,
    plotTitleAutoTransform,
    plotTitleUserOffset,
    generatedLegendPosition
  } = state;

  const currentSvg = () => svgContainer.value?.querySelector?.('svg') || null;

  // One canvas padding applies to every Result (D-09, PD-OI-064): Generate
  // applies it before the candidate is published, and a displayed Result
  // receives the current value. Applying the same padding twice is a no-op.
  const applyCanvasPadding = () => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    const svg = currentSvg();
    if (!svg) return false;
    bindCompositionMetadata(svg);
    if (!applyCanvasPaddingToSvg(svg, drawing.canvasPadding)) return false;
    return Boolean(commitActiveResultEdit?.('canvas-padding'));
  };

  const resetCanvasPadding = () => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    drawing.canvasPadding.top = 0;
    drawing.canvasPadding.right = 0;
    drawing.canvasPadding.bottom = 0;
    drawing.canvasPadding.left = 0;
    return applyCanvasPadding();
  };

  const captureBaseConfig = () => {
    const svg = currentSvg();
    if (!svg) return null;
    const binding = bindCompositionMetadata(svg);
    const { metadata } = binding;
    const deltas = compositionUserDeltas(svg);

    diagramElements.value = [...binding.primary.targets];
    diagramElementOriginalTransforms.value = new Map(
      binding.primary.targets.map((target) => [
        target,
        {
          x: metadata.primary.automaticTranslation[0],
          y: metadata.primary.automaticTranslation[1]
        }
      ])
    );
    diagramOffset.x = deltas.primary[0]?.[0] || 0;
    diagramOffset.y = deltas.primary[0]?.[1] || 0;
    generatedLegendPosition.value = metadata.legendSide;

    if (metadata.legend) {
      legendInitialTransform.value = {
        x: metadata.legend.automaticTranslation[0],
        y: metadata.legend.automaticTranslation[1]
      };
      legendCurrentOffset.x = deltas.legend?.[0] || 0;
      legendCurrentOffset.y = deltas.legend?.[1] || 0;
    } else {
      legendInitialTransform.value = { x: 0, y: 0 };
      legendCurrentOffset.x = 0;
      legendCurrentOffset.y = 0;
    }
    if (metadata.title) {
      plotTitleAutoTransform.value = {
        x: metadata.title.automaticTranslation[0],
        y: metadata.title.automaticTranslation[1]
      };
      plotTitleUserOffset.x = deltas.title?.[0] || 0;
      plotTitleUserOffset.y = deltas.title?.[1] || 0;
    } else {
      plotTitleAutoTransform.value = { x: 0, y: 0 };
      plotTitleUserOffset.x = 0;
      plotTitleUserOffset.y = 0;
    }

    return binding;
  };

  return {
    applyCanvasPadding,
    captureBaseConfig,
    resetCanvasPadding
  };
};
