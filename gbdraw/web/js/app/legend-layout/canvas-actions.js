import { serializeCleanSvg } from '../../services/svg-serialization.js';
import {
  applyCanvasPaddingToSvg,
  bindCompositionMetadata,
  compositionUserDeltas
} from './composition-actions.js';

export const createLegendCanvasActions = ({ state }) => {
  const {
    svgContainer,
    canvasPadding,
    originalSvgStroke,
    diagramElements,
    diagramElementOriginalTransforms,
    diagramOffset,
    legendInitialTransform,
    legendCurrentOffset,
    plotTitleAutoTransform,
    plotTitleUserOffset,
    generatedLegendPosition,
    selectedResultIndex,
    results,
    skipCaptureBaseConfig,
    skipPositionReapply
  } = state;

  const currentSvg = () => svgContainer.value?.querySelector?.('svg') || null;

  const persistCurrentSvg = (svg = currentSvg()) => {
    const index = selectedResultIndex.value;
    if (!svg || index < 0 || index >= results.value.length) return false;
    skipCaptureBaseConfig.value = true;
    skipPositionReapply.value = true;
    const nextResults = [...results.value];
    nextResults[index] = {
      ...results.value[index],
      content: serializeCleanSvg(svg)
    };
    results.value = nextResults;
    return true;
  };

  // One canvas padding applies to every Result (D-09, PD-OI-064): Generate
  // applies it before the candidate is published, and a displayed Result
  // receives the current value. Applying the same padding twice is a no-op.
  const applyCanvasPadding = () => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    const svg = currentSvg();
    if (!svg) return false;
    bindCompositionMetadata(svg);
    if (!applyCanvasPaddingToSvg(svg, canvasPadding)) return false;
    return persistCurrentSvg(svg);
  };

  const resetCanvasPadding = () => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    canvasPadding.top = 0;
    canvasPadding.right = 0;
    canvasPadding.bottom = 0;
    canvasPadding.left = 0;
    return applyCanvasPadding();
  };

  const captureOriginalStroke = () => {
    const svg = currentSvg();
    if (!svg) return;
    const firstFeaturePath = svg.querySelector('path[id^="f"]');
    if (!firstFeaturePath) return;
    const strokeWidth = Number.parseFloat(firstFeaturePath.getAttribute('stroke-width'));
    originalSvgStroke.value = {
      color: firstFeaturePath.getAttribute('stroke'),
      width: Number.isFinite(strokeWidth) ? strokeWidth : null
    };
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
    captureOriginalStroke,
    persistCurrentSvg,
    resetCanvasPadding
  };
};
