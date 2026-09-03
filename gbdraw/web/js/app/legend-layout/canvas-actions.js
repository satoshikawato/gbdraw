import { serializeCleanSvg } from '../../services/svg-serialization.js';
import {
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

  const applyCanvasPadding = () => {
    const svg = currentSvg();
    if (!svg) return;
    bindCompositionMetadata(svg);
    const currentViewBox = svg.dataset.originalViewBox || svg.getAttribute('viewBox');
    if (!svg.dataset.originalViewBox) svg.dataset.originalViewBox = currentViewBox;
    const [baseX, baseY, baseWidth, baseHeight] = String(svg.dataset.originalViewBox)
      .trim().split(/\s+/).map(Number);
    if ([baseX, baseY, baseWidth, baseHeight].some((value) => !Number.isFinite(value))) {
      throw new Error('The SVG base viewBox is invalid. Regenerate the diagram.');
    }
    const originalWidth = Number.parseFloat(svg.dataset.originalWidth || svg.getAttribute('width')) || baseWidth;
    const originalHeight = Number.parseFloat(svg.dataset.originalHeight || svg.getAttribute('height')) || baseHeight;
    if (!svg.dataset.originalWidth) {
      svg.dataset.originalWidth = String(originalWidth);
      svg.dataset.originalHeight = String(originalHeight);
    }
    const width = baseWidth + canvasPadding.left + canvasPadding.right;
    const height = baseHeight + canvasPadding.top + canvasPadding.bottom;
    svg.setAttribute('viewBox', `${baseX - canvasPadding.left} ${baseY - canvasPadding.top} ${width} ${height}`);
    svg.setAttribute('width', `${originalWidth * (width / baseWidth)}px`);
    svg.setAttribute('height', `${originalHeight * (height / baseHeight)}px`);
    persistCurrentSvg(svg);
  };

  const resetCanvasPadding = () => {
    canvasPadding.top = 0;
    canvasPadding.right = 0;
    canvasPadding.bottom = 0;
    canvasPadding.left = 0;
    const svg = currentSvg();
    if (!svg) return;
    bindCompositionMetadata(svg);
    if (svg.dataset.originalViewBox) svg.setAttribute('viewBox', svg.dataset.originalViewBox);
    if (svg.dataset.originalWidth) {
      svg.setAttribute('width', `${svg.dataset.originalWidth}px`);
      svg.setAttribute('height', `${svg.dataset.originalHeight}px`);
    }
    persistCurrentSvg(svg);
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
