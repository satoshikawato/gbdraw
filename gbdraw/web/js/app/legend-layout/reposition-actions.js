// @ts-check
/** @import { DrawingState } from '../../state.js' */
import {
  applyCompositionEdit,
  bindCompositionMetadata,
  compositionUserDeltas
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

/**
 * @typedef {object} LegendRepositionActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(svg: SVGSVGElement, options?: { side?: string }) => import('../../services/legend-layout.js').LayoutBox | null} layOutLegend
 *   The Legend manager's layout of the Legend for a side as Python lays it out;
 *   returns its local bounds, or null when it laid nothing out.
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 */

/** @param {LegendRepositionActionsOptions} options */
export const createLegendRepositionActions = ({
  state,
  layOutLegend,
  commitActiveResultEdit = null
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
    plotTitleUserOffset
  } = state;
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

  /**
   * @param {DrawingState} drawing
   * @param {string} newPosition
   * @param {string} [_oldPosition]
   * @param {{ preserveManualOffsets?: boolean, commit?: boolean }} [options]
   *   `commit: false` leaves the commit to the caller (the mount binder).
   */
  const repositionForLegendChange = (drawing, newPosition, _oldPosition, options = {}) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!svgContainer.value || !svgContent.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const binding = bindCompositionMetadata(svg);
    const legendGroup = binding.legend.targets[0] || null;
    // A diagram generated without a legend cannot show one in place; the
    // Result stays unchanged and the next Generate applies the side (GE-07).
    if (!legendGroup && newPosition !== 'none') return false;

    // The Legend is laid out for the side as Python lays it out (zero shift),
    // and docked with the bounds that layout gives.
    /** @type {import('../../services/legend-layout.js').LayoutBox | null} */
    let legendLocalBox = null;
    if (legendGroup && newPosition !== 'none') {
      legendGroup.removeAttribute('display');
      setLegendVariant(legendGroup, newPosition);
      legendLocalBox = layOutLegend(svg, { side: newPosition });
    }

    const nextBinding = applyCompositionEdit(svg, { legendSide: newPosition, canvasPadding: drawing.canvasPadding, legendLocalBox });
    syncStateFromComposition(svg, nextBinding);
    if (options.commit !== false) commitActiveResultEdit?.('legend-position');
    return true;
  };

  /** @param {{ commit?: boolean }} [options] */
  const refreshLegendGeometry = ({ commit = true } = {}) => {
    const drawing = state.activeDrawing();
    if (!svgContainer.value || !svgContent.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;
    const binding = bindCompositionMetadata(svg);
    if (!binding.legend.metadata || binding.metadata.legendSide === 'none') return false;
    return repositionForLegendChange(drawing, binding.metadata.legendSide, binding.metadata.legendSide, {
      preserveManualOffsets: true,
      commit
    });
  };

  return {
    refreshLegendGeometry,
    syncStateFromComposition
  };
};
