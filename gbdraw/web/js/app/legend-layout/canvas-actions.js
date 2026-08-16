import {
  bindCompositionMetadata,
  compositionUserDeltas
} from './composition-actions.js';

export const createLegendCanvasActions = ({ state, previewRuntime }) => {
  if (
    typeof previewRuntime?.commitDomEdit !== 'function'
    || typeof previewRuntime?.runDomEditSync !== 'function'
  ) {
    throw new Error('createLegendCanvasActions requires the preview runtime edit protocol.');
  }
  const {
    svgContainer,
    canvasPadding,
    originalSvgStroke,
    diagramElements,
    diagramElementOriginalTransforms,
    diagramOffset,
    lengthBarUserOffset,
    legendInitialTransform,
    legendCurrentOffset,
    plotTitleAutoTransform,
    plotTitleUserOffset,
    generatedLegendPosition,
    skipPositionReapply
  } = state;

  const currentSvg = () => svgContainer.value?.querySelector?.('svg') || null;
  const canonicalCanvasState = () => [
    { target: canvasPadding },
    { target: diagramElements, key: 'value' },
    { target: diagramElementOriginalTransforms, key: 'value', deep: true },
    { target: diagramOffset },
    { target: lengthBarUserOffset },
    { target: legendInitialTransform, key: 'value', deep: true },
    { target: legendCurrentOffset },
    { target: plotTitleAutoTransform, key: 'value', deep: true },
    { target: plotTitleUserOffset },
    { target: generatedLegendPosition, key: 'value' },
    { target: skipPositionReapply, key: 'value' }
  ];
  const runCanvasAction = (reason, action) => previewRuntime.runDomEditSync({
    reason,
    canonicalState: canonicalCanvasState(),
    action
  });

  const persistCurrentSvg = (
    svg = currentSvg(),
    mutate,
    reason = 'canvas-layout',
    { lease = null } = {}
  ) => {
    if (!svg) return false;
    if (typeof mutate !== 'function') {
      throw new TypeError('persistCurrentSvg requires the SVG mutation callback.');
    }
    const apply = ({ svg: targetSvg, mutation }) => {
        if (targetSvg !== svg) return false;
        mutation.captureProperty(skipPositionReapply, 'value');
        return mutate({ svg: targetSvg, mutation });
      };
    const changed = lease
      ? lease.mutate(apply)
      : previewRuntime.commitDomEdit({
          reason,
          invalidateIndexes: ['legend'],
          mutate: apply
        }).changed;
    if (changed) skipPositionReapply.value = true;
    return changed;
  };

  const applyCanvasPadding = () => {
    const svg = currentSvg();
    if (!svg) return;
    return persistCurrentSvg(svg, ({ svg: targetSvg, mutation }) => {
      bindCompositionMetadata(targetSvg);
      const currentViewBox = targetSvg.dataset.originalViewBox || targetSvg.getAttribute('viewBox');
      const originalViewBox = targetSvg.dataset.originalViewBox || currentViewBox;
      const [baseX, baseY, baseWidth, baseHeight] = String(originalViewBox)
        .trim().split(/\s+/).map(Number);
      if ([baseX, baseY, baseWidth, baseHeight].some((value) => !Number.isFinite(value))) {
        throw new Error('The SVG base viewBox is invalid. Regenerate the diagram.');
      }
      const originalWidth = Number.parseFloat(targetSvg.dataset.originalWidth || targetSvg.getAttribute('width')) || baseWidth;
      const originalHeight = Number.parseFloat(targetSvg.dataset.originalHeight || targetSvg.getAttribute('height')) || baseHeight;
      let changed = false;
      if (!targetSvg.dataset.originalViewBox) {
        changed = mutation.setAttribute(targetSvg, 'data-original-view-box', currentViewBox) || changed;
      }
      if (!targetSvg.dataset.originalWidth) {
        changed = mutation.setAttribute(targetSvg, 'data-original-width', originalWidth) || changed;
        changed = mutation.setAttribute(targetSvg, 'data-original-height', originalHeight) || changed;
      }
      const width = baseWidth + canvasPadding.left + canvasPadding.right;
      const height = baseHeight + canvasPadding.top + canvasPadding.bottom;
      changed = mutation.setAttribute(
        targetSvg,
        'viewBox',
        `${baseX - canvasPadding.left} ${baseY - canvasPadding.top} ${width} ${height}`
      ) || changed;
      changed = mutation.setAttribute(targetSvg, 'width', `${originalWidth * (width / baseWidth)}px`) || changed;
      changed = mutation.setAttribute(targetSvg, 'height', `${originalHeight * (height / baseHeight)}px`) || changed;
      return changed;
    });
  };

  const resetCanvasPadding = () => runCanvasAction('canvas-padding-reset', () => {
    canvasPadding.top = 0;
    canvasPadding.right = 0;
    canvasPadding.bottom = 0;
    canvasPadding.left = 0;
    return applyCanvasPadding();
  });

  const updateCanvasPadding = (side, value) => {
    if (!['top', 'right', 'bottom', 'left'].includes(side)) return false;
    const nextValue = Number(value);
    if (!Number.isFinite(nextValue) || nextValue < 0 || canvasPadding[side] === nextValue) {
      return false;
    }
    return runCanvasAction('canvas-padding-update', () => {
      canvasPadding[side] = nextValue;
      return applyCanvasPadding();
    });
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
    resetCanvasPadding,
    runCanvasAction,
    updateCanvasPadding
  };
};
