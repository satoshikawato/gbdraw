import { serializeCleanSvg } from '../../services/svg-serialization.js';
import {
  applyCompositionEdit,
  bindCompositionMetadata,
  compositionUserDeltas,
  reconcileCompositionTitle
} from './composition-actions.js';

const isHorizontalSide = (side) => side === 'top' || side === 'bottom';

const setLegendVariant = (legendGroup, side) => {
  const horizontal = legendGroup?.querySelector?.('#legend_horizontal') || null;
  const vertical = legendGroup?.querySelector?.('#legend_vertical') || null;
  if (!horizontal || !vertical) return false;
  if (isHorizontalSide(side)) {
    horizontal.removeAttribute('display');
    vertical.setAttribute('display', 'none');
  } else {
    horizontal.setAttribute('display', 'none');
    vertical.removeAttribute('display');
  }
  return true;
};

export const createLegendRepositionActions = ({
  state,
  legendActions
}) => {
  const {
    svgContent,
    svgContainer,
    generatedLegendPosition,
    diagramElements,
    diagramElementOriginalTransforms,
    diagramOffset,
    legendInitialTransform,
    legendCurrentOffset,
    plotTitleAutoTransform,
    plotTitleUserOffset,
    canvasPadding,
    selectedResultIndex,
    results,
    skipCaptureBaseConfig,
    skipPositionReapply
  } = state;
  const {
    reflowDualLegendLayout,
    reflowSingleLegendLayout
  } = legendActions;
  const persist = (svg) => {
    skipCaptureBaseConfig.value = true;
    skipPositionReapply.value = true;
    const index = selectedResultIndex.value;
    if (index >= 0 && index < results.value.length) {
      const nextResults = [...results.value];
      nextResults[index] = {
        ...results.value[index],
        content: serializeCleanSvg(svg)
      };
      results.value = nextResults;
    }
  };

  const syncStateFromComposition = (svg, binding = bindCompositionMetadata(svg)) => {
    const { metadata } = binding;
    const deltas = compositionUserDeltas(svg);
    generatedLegendPosition.value = metadata.legendSide;
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
  };

  const repositionForLegendChange = (newPosition, _oldPosition, _options = {}) => {
    if (!svgContainer.value || !svgContent.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const binding = bindCompositionMetadata(svg);
    const legendGroup = binding.legend.targets[0] || null;
    if (!legendGroup && newPosition !== 'none') {
      throw new Error('This diagram has no legend composition target. Regenerate it with a legend before changing its side.');
    }

    if (legendGroup && newPosition !== 'none') {
      legendGroup.removeAttribute('display');
      const hasDualLegend = setLegendVariant(legendGroup, newPosition);
      if (hasDualLegend) {
        reflowDualLegendLayout(svg);
      } else {
        const widthHint = isHorizontalSide(newPosition)
          ? binding.metadata.primary.finalBounds.width
          : null;
        reflowSingleLegendLayout(
          svg,
          isHorizontalSide(newPosition) ? 'horizontal' : 'vertical',
          widthHint
        );
      }
    }

    const nextBinding = applyCompositionEdit(svg, { legendSide: newPosition, canvasPadding });
    syncStateFromComposition(svg, nextBinding);
    persist(svg);
    return true;
  };

  const refreshLegendGeometry = () => {
    if (!svgContainer.value || !svgContent.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;
    const binding = bindCompositionMetadata(svg);
    if (!binding.legend.metadata || binding.metadata.legendSide === 'none') return false;
    return repositionForLegendChange(binding.metadata.legendSide, binding.metadata.legendSide, {
      preserveManualOffsets: true
    });
  };

  const refreshCompositionGeometry = ({ titleSide = null, titleTarget = undefined } = {}) => {
    if (!svgContainer.value || !svgContent.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;
    const binding = titleTarget === undefined
      ? applyCompositionEdit(svg, {
          titleSide: titleSide ?? bindCompositionMetadata(svg).metadata.titleSide,
          canvasPadding
        })
      : reconcileCompositionTitle(
          svg,
          titleTarget,
          titleSide ?? (titleTarget ? 'bottom' : 'none'),
          { canvasPadding }
        );
    syncStateFromComposition(svg, binding);
    persist(svg);
    return true;
  };

  return {
    refreshCompositionGeometry,
    refreshLegendGeometry,
    repositionForLegendChange,
    syncStateFromComposition
  };
};
