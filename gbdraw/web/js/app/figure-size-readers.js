// @ts-check
// The readers that the Layout card's Figure size and Wrap Rows fields (app/figure-size-field.js,
// app/linear-wrap-field.js) take from the composition root: the shown Result's size beside the
// requested sides, and what the shown Linear drawing holds that rows cannot draw with. The root
// creates both fields with these readers. It passes the state, the shown Result's edit and mount
// revisions, the committed request of the Result a drawing shows, the moved decorations, and the
// Linear labels' shared row as ports.
/** @import { DrawingState, WebState } from '../state.js' */
/** @import { CanonicalRenderRequest } from '../contracts/request.generated.js' */
/** @import { LinearWrapBlockerInputs } from '../services/linear-wrap-form.js' */
/** @import { ShownFigureSize } from './figure-size-field.js' */
import { figureGrowth, requestFigureSizeLeaves, requestedLengthText, shownFigureSize } from '../services/figure-size-form.js';
import { effectiveLinearSlotPlacement } from '../services/linear-track-slot-model.js';
import { movesARecordSideways } from '../services/record-draw-selection.js';
import { getSvgDimensions } from '../services/svg-dimensions.js';
import { draftFeatureSlot } from './feature-editor/placement-actions.js';
import { requestedCanvasOf } from './legend-layout/composition-actions.js';

const { shallowRef, watch } = window.Vue;

/**
 * @typedef {object} FigureSizeReaderPorts
 * @property {WebState} state
 * @property {{ readonly value: number }} resultEditRevision Counts the edits of the shown Result's root.
 * @property {{ readonly value: number }} mountedResultRevision Counts the Results the preview mounted.
 * @property {(drawing: DrawingState) => CanonicalRenderRequest | null} shownCommittedRequest The committed request
 *   of the Result the drawing shows; null while it shows no Result of its own.
 * @property {() => boolean} decorationsMoved Whether a Web-only change moved the diagram, Legend, title, or scale.
 * @property {{ readonly value: boolean }} linearLabelHasSharedRow Whether a Linear label row is shared by records.
 */

/**
 * @typedef {object} FigureSizeReaders
 * @property {() => ShownFigureSize | null} shownSize The Figure size field's reading of the shown Result.
 * @property {() => LinearWrapBlockerInputs} blockerInputs The Wrap Rows field's blocker inputs.
 */

/**
 * Registers the watch that reads the shown Result's root; the root calls it where that watch stood.
 * @param {FigureSizeReaderPorts} ports
 * @returns {FigureSizeReaders}
 */
export const createFigureSizeReaders = ({
  state, resultEditRevision, mountedResultRevision, shownCommittedRequest, decorationsMoved, linearLabelHasSharedRow
}) => {
  const { svgContent, svgContainer, labelReflowProcessing } = state;
  // Figure size (issue #1427): the line under the field reads the shown Result's root (px) after each change of it,
  // and its committed request's sides; none while the drawing shows no Result of its own.
  // A root drawn to a requested size also names its requested canvas and the canvas before Web-only additions grew it.
  // While the automatic rerender runs on such a root (the shown drawing: a switch waits for it), the line says the
  // figure is being fitted (PD-OI-114 rev 2): a Legend edit or Canvas padding shows at once, then Python fits it.
  /** @type {{ value: { px: { width: number, height: number }, sized: ReturnType<typeof requestedCanvasOf> } | null }} */
  const shownFigureRoot = shallowRef(null);
  watch(() => [svgContent.value, resultEditRevision.value, mountedResultRevision.value, state.activeDrawingId.value], () => {
    const root = svgContainer.value?.querySelector?.('svg') ?? null;
    const px = getSvgDimensions(root);
    shownFigureRoot.value = root && px ? { px, sized: requestedCanvasOf(root) } : null;
  }, { flush: 'post' });
  const shownSize = () => {
    const drawing = state.activeDrawing();
    const request = shownCommittedRequest(drawing);
    const root = shownFigureRoot.value;
    if (!request || !root) return null;
    const requested = requestFigureSizeLeaves(request);
    const growth = root.sized && figureGrowth({
      shown: root.px, base: root.sized.base, requested: root.sized.requested,
      padded: root.sized.addsPadding(drawing.canvasPadding), moved: decorationsMoved()
    });
    return {
      ...shownFigureSize({ widthPx: root.px.width, heightPx: root.px.height, requested }),
      fitting: Boolean(root.sized) && labelReflowProcessing.value,
      larger: growth ? {
        sides: growth.sides.map((side) => ({ side, length: requestedLengthText(requested[side]) ?? `${Math.round(root.px[side])} px` })),
        causes: growth.causes
      } : null
    };
  };
  // Wrap Rows (PD-OI-113): its blockers read the owners the request reads. The comparison blocker reads the comparison
  // inputs only; a drawn record moved sideways is read from the translations the request carries, whatever wrote them
  // (Align at feature counts only through the records its Apply moved), and a Similarity alignment plan counts as
  // moved records, also when its comparison inputs were removed (Q5 A).
  const blockerInputs = () => {
    const drawing = state.activeDrawing();
    return {
      comparisons: drawing.linearComparisonResolution.value.hasComparisonIntent,
      sharedRecordRow: linearLabelHasSharedRow.value,
      movedRecords: movesARecordSideways(state.linearRecordTranslations.value, drawing.drawnLinearSeqs.value),
      alignmentPlan: state.similarityAlignmentPlan.value !== null,
      featureTrackOnAxis: Boolean(drawing.adv.linear_track_slots_enabled)
        && effectiveLinearSlotPlacement(draftFeatureSlot({ mode: drawing.mode, form: drawing.form, adv: drawing.adv })) === 'overlay'
    };
  };
  return { shownSize, blockerInputs };
};
