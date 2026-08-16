import {
  applyCompositionEdit,
  bindCompositionMetadata,
  compositionUserDeltas,
  reconcileCompositionTitle
} from './composition-actions.js';

const isHorizontalSide = (side) => side === 'top' || side === 'bottom';

const setLegendVariant = (legendGroup, side, mutation = null) => {
  const horizontal = legendGroup?.querySelector?.('#legend_horizontal') || null;
  const vertical = legendGroup?.querySelector?.('#legend_vertical') || null;
  if (!horizontal || !vertical) return false;
  if (isHorizontalSide(side)) {
    if (mutation) {
      mutation.removeAttribute(horizontal, 'display');
      mutation.setAttribute(vertical, 'display', 'none');
    } else {
      horizontal.removeAttribute('display');
      vertical.setAttribute('display', 'none');
    }
  } else {
    if (mutation) {
      mutation.setAttribute(horizontal, 'display', 'none');
      mutation.removeAttribute(vertical, 'display');
    } else {
      horizontal.setAttribute('display', 'none');
      vertical.removeAttribute('display');
    }
  }
  return true;
};

export const createLegendRepositionActions = ({
  state,
  legendActions,
  previewRuntime
}) => {
  if (typeof previewRuntime?.commitDomEdit !== 'function') {
    throw new Error('createLegendRepositionActions requires the preview runtime edit protocol.');
  }
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
    skipPositionReapply
  } = state;
  const {
    reflowDualLegendLayout,
    reflowSingleLegendLayout
  } = legendActions;
  const captureRepositionState = (mutation) => {
    mutation.captureProperty(generatedLegendPosition, 'value');
    mutation.captureProperty(diagramElements, 'value');
    mutation.captureProperty(diagramElementOriginalTransforms, 'value', { deep: true });
    mutation.captureState(diagramOffset);
    mutation.captureProperty(legendInitialTransform, 'value', { deep: true });
    mutation.captureState(legendCurrentOffset);
    mutation.captureProperty(plotTitleAutoTransform, 'value', { deep: true });
    mutation.captureState(plotTitleUserOffset);
    mutation.captureProperty(skipPositionReapply, 'value');
  };
  const commitReposition = (svg, mutate, { lease = null } = {}) => {
    const apply = ({ svg: targetSvg, mutation }) => {
      if (targetSvg !== svg) return false;
      captureRepositionState(mutation);
      return mutate({ svg: targetSvg, mutation });
    };
    const changed = lease
      ? lease.mutate(apply)
      : previewRuntime.commitDomEdit({
      reason: 'legend-reposition',
      invalidateIndexes: ['legend'],
      mutate: apply
    }).changed;
    if (changed) skipPositionReapply.value = true;
    return changed;
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

  const repositionForLegendChange = (newPosition, _oldPosition, options = {}) => {
    if (!svgContainer.value || !svgContent.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const binding = bindCompositionMetadata(svg);
    const legendGroup = binding.legend.targets[0] || null;
    if (!legendGroup && newPosition !== 'none') {
      throw new Error('This diagram has no legend composition target. Regenerate it with a legend before changing its side.');
    }

    return commitReposition(svg, ({ svg: targetSvg, mutation }) => {
      if (legendGroup && newPosition !== 'none') {
        mutation.removeAttribute(legendGroup, 'display');
        const hasDualLegend = setLegendVariant(legendGroup, newPosition, mutation);
        if (hasDualLegend) {
          reflowDualLegendLayout(targetSvg, mutation);
        } else {
          const widthHint = isHorizontalSide(newPosition)
            ? binding.metadata.primary.finalBounds.width
            : null;
          reflowSingleLegendLayout(
            targetSvg,
            isHorizontalSide(newPosition) ? 'horizontal' : 'vertical',
            widthHint,
            mutation
          );
        }
      }

      const nextBinding = applyCompositionEdit(targetSvg, {
        legendSide: newPosition,
        canvasPadding,
        mutation
      });
      syncStateFromComposition(targetSvg, nextBinding);
      return true;
    }, options);
  };

  const refreshLegendGeometry = (options = {}) => {
    if (!svgContainer.value || !svgContent.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;
    const binding = bindCompositionMetadata(svg);
    if (!binding.legend.metadata || binding.metadata.legendSide === 'none') return false;
    return repositionForLegendChange(binding.metadata.legendSide, binding.metadata.legendSide, {
      preserveManualOffsets: true,
      ...options
    });
  };

  const refreshCompositionGeometry = ({ titleSide = null, titleTarget = undefined } = {}) => {
    if (!svgContainer.value || !svgContent.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;
    return commitReposition(svg, ({ svg: targetSvg, mutation }) => {
      const binding = titleTarget === undefined
        ? applyCompositionEdit(targetSvg, {
            titleSide: titleSide ?? bindCompositionMetadata(targetSvg).metadata.titleSide,
            canvasPadding,
            mutation
          })
        : reconcileCompositionTitle(
            targetSvg,
            titleTarget,
            titleSide ?? (titleTarget ? 'bottom' : 'none'),
            { canvasPadding, mutation }
          );
      syncStateFromComposition(targetSvg, binding);
      return true;
    });
  };

  return {
    refreshCompositionGeometry,
    refreshLegendGeometry,
    repositionForLegendChange,
    syncStateFromComposition
  };
};
