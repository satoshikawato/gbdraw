// @ts-check
// The Review alignment panel of a Similarity alignment, docked in the Result workspace, and the
// review's takeover of the preview: the figure overlay and the strips, the measures taken again
// after each change of the figure, the fitted view and the closed drawer the review holds while it
// is open and gives back at its end, the focus it returns, Escape, Reset, and Apply. The alignment
// owner (app/similarity-alignment.js) keeps the draft and the plan; the composition root passes its
// members, the pan/zoom owner's hold, the preview runtime's Result identity, the artifact slots'
// read, and the Editor drawer, the settings popup and the feature popup as ports, and fills
// `similarityAlignmentPorts` from what the review returns.
/** @import { WebState } from '../state.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { ArtifactSlotStore } from '../services/artifact-slot.js' */
/** @import { AlignmentFigureGeometry, AlignmentOverlayRecord } from '../services/alignment-review-view.js' */
/** @import { AliasSource } from '../services/feature-identity.js' */
/** @import { SvgResultEntry } from '../services/svg-result-ingestion.js' */
/** @import { DirectionPreviewRecord } from './similarity-alignment.js' */
import { ALIGNMENT_REFERENCE_KEY, alignmentOverlayRecords, buildAlignmentReviewView } from '../services/alignment-review-view.js';
import { t } from '../utils/i18n.js';

const { computed, nextTick, onMounted, onUnmounted, ref, shallowRef, watch } = window.Vue;

/**
 * One record of the alignment draft, as the figure overlay and the Review panel read it.
 * @typedef {{
 *   recordKey: string, recordLabel: string, choice: { kind: string, candidateKey?: string | null } | null,
 *   candidates: readonly { key: string, anchor: Record<string, unknown>, name: string, label: string, coordinates: string }[]
 * }} AlignmentDraftRow
 */
/**
 * The alignment owner's draft, as the review reads it.
 * @typedef {{
 *   rows: readonly AlignmentDraftRow[], reference: { recordKey: string, name: string },
 *   response: { reference: AliasSource }, openedWithoutChoice: readonly string[]
 * }} AlignmentDraft
 */
/** @typedef {Record<string, unknown>} AlignmentOutcome What the alignment owner's Apply and Reset return; the review reads its `status`. */
/**
 * What the review took when it opened: the Editor drawer and its tab, the feature popup, and the Result of the
 * review's drawing it last showed its draft on.
 * @typedef {{ drawerOpen: boolean, drawerTab: string, popup: unknown, drawingId: DrawingId, figure: string }} AlignmentTakeover
 */
/** @typedef {{ value: AlignmentTakeover | null }} AlignmentTakeoverRef The root's takeover ref; the preview's gestures read it. */
/**
 * What the figure overlay draws: the reference and each record's candidates.
 * @typedef {{
 *   reference: { key: string, recordKey: string, anchor: AliasSource, name: string }, records: AlignmentOverlayRecord[],
 *   onSelect: (recordKey: string, anchor: AliasSource) => void, onHover: (recordKey: string, candidateKey: string) => void
 * }} AlignmentOverlayRequest
 */

/**
 * @typedef {object} SimilarityAlignmentTakeoverPorts
 * @property {WebState} state
 * @property {AlignmentTakeoverRef} similarityAlignmentTakeover
 * @property {() => boolean} similarityAlignmentHoldShown Whether the shown drawing is the review's.
 * @property {(options: { restore: boolean }) => void} releasePreviewViewport The pan/zoom owner's release of its hold.
 * @property {() => void} fitPreviewToViewport
 * @property {(result: SvgResultEntry | null | undefined) => string} getResultIdentity The preview runtime's identity of
 *   a Result.
 * @property {ArtifactSlotStore['readOfDrawing']} readOfDrawing
 * @property {(tab: string) => unknown} openRightDrawerTab
 * @property {() => Promise<unknown>} showFeaturePopupAtFeature The feature popup opened again beside its feature.
 * @property {{ readonly value: boolean }} alignmentDialogOpen
 * @property {() => Promise<AlignmentOutcome>} applyAlignmentReset
 * @property {() => Promise<AlignmentOutcome>} applyAlignmentDraft
 */

/**
 * @typedef {object} SimilarityAlignmentTakeover
 * @property {(reviewDrawingId: DrawingId) => string} reviewDrawingFigure
 * @property {() => void} dropSimilarityAlignmentView
 * @property {() => void} fitPreviewForSimilarityAlignment
 * @property {() => void} endSimilarityAlignmentTakeover
 * @property {(event: Event | null | undefined) => void} rememberSimilarityAlignmentInvoker
 * @property {() => Promise<void>} restoreSimilarityAlignmentFocus
 * @property {() => Promise<void>} focusSimilarityAlignmentDialog
 * @property {(reset?: boolean) => Promise<unknown>} applySimilarityAlignmentDialog The review's Apply or
 *   Reset; the session's busy answer while a Session operation runs
 */

/**
 * The review's takeover and its end, the focus it returns, and Apply. Registers no watch, computed or hook.
 * @param {SimilarityAlignmentTakeoverPorts} ports
 * @returns {SimilarityAlignmentTakeover}
 */
export const createSimilarityAlignmentTakeover = ({
  state, similarityAlignmentTakeover, similarityAlignmentHoldShown, releasePreviewViewport, fitPreviewToViewport,
  getResultIdentity, readOfDrawing, openRightDrawerTab, showFeaturePopupAtFeature, alignmentDialogOpen,
  applyAlignmentReset, applyAlignmentDraft
}) => {
  const { results, selectedResultIndex, clickedFeature, showRightDrawer } = state;
  // The takeover ends with the review (Cancel, the close button, Escape, Apply). The pan/zoom
  // owner holds the view from before the fit, so a Session saved or a Result checkpointed
  // meanwhile keeps that view. While Apply runs, the ending waits for its outcome, and its new
  // Result is the review's own.
  let similarityAlignmentApplying = false;
  // The identity of the selected Result of the review's drawing, shown or kept in its slot.
  /** @param {DrawingId} reviewDrawingId */
  const reviewDrawingFigure = (reviewDrawingId) => getResultIdentity(readOfDrawing(reviewDrawingId,
    state.activeDrawingId.value, () => results.value[Number(selectedResultIndex.value) || 0],
    (slot) => slot.values.results[slot.values.selectedResultIndex]));
  // The draft no longer matches the figure. Another Result of the review's mode replaced it (a
  // Session load, History, a new Generate), whichever mode is shown: the view and the drawer are
  // that change's, and the review gives neither back. A mode switch brings no view: the review's
  // own Result, kept in its slot or shown again, leaves the hold, so the end gives the view back.
  const dropSimilarityAlignmentView = () => {
    const takeover = similarityAlignmentTakeover.value;
    if (!takeover || similarityAlignmentApplying || reviewDrawingFigure(takeover.drawingId) === takeover.figure) return;
    takeover.drawerOpen = false;
    releasePreviewViewport({ restore: false });
  };
  /** @type {HTMLElement | null} */
  let similarityAlignmentReturnFocus = null;
  // The figure is fitted to the room the docked panel leaves, once the panel is in the page. The
  // pan/zoom owner fits it again when that room changes, for as long as it holds the view.
  const fitPreviewForSimilarityAlignment = () => {
    if (alignmentDialogOpen.value) fitPreviewToViewport();
  };
  const endSimilarityAlignmentTakeover = () => {
    const takeover = similarityAlignmentTakeover.value;
    if (!takeover || similarityAlignmentApplying) return;
    // The held view is the review's mode's: ended while the other mode is shown, the review
    // does not write it onto that mode's figure (OV-1089).
    const restore = similarityAlignmentHoldShown();
    similarityAlignmentTakeover.value = null;
    if (takeover.drawerOpen) openRightDrawerTab(takeover.drawerTab);
    // The drawer narrows the visible frame, so the view returns after it.
    releasePreviewViewport({ restore });
    // Open in the feature search row retargeted the hidden popup for the fitted view.
    if (clickedFeature.value && clickedFeature.value !== takeover.popup) void showFeaturePopupAtFeature();
  };
  /** @param {Event | null | undefined} event */
  const rememberSimilarityAlignmentInvoker = (event) => {
    // `currentTarget` is the element the handler is bound to; the declared type of the variable already holds the focused element.
    similarityAlignmentReturnFocus = /** @type {HTMLElement | null} */ (event?.currentTarget || document.activeElement || null);
  };
  const restoreSimilarityAlignmentFocus = async () => {
    await nextTick();
    const target = similarityAlignmentReturnFocus;
    if (!target) return;
    similarityAlignmentReturnFocus = null;
    const hiddenEditorInvoker = target?.closest?.('.right-drawer') && !showRightDrawer.value;
    if (!hiddenEditorInvoker && target?.isConnected && target.getClientRects().length && typeof target.focus === 'function') {
      target.focus();
      return;
    }
    /** @type {HTMLElement | null} */ (document.querySelector('.drawer-toggle'))?.focus();
  };
  const focusSimilarityAlignmentDialog = async () => {
    await nextTick();
    document.getElementById('similarity-alignment-dialog-title')?.focus();
  };
  const applySimilarityAlignmentDialog = async (reset = false) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!reset) similarityAlignmentApplying = true;
    let outcome;
    try {
      outcome = await (reset ? applyAlignmentReset() : applyAlignmentDraft());
    } finally {
      if (!reset) similarityAlignmentApplying = false;
    }
    if (!reset) {
      if (!alignmentDialogOpen.value) {
        endSimilarityAlignmentTakeover();
        void restoreSimilarityAlignmentFocus();
      } else if (outcome.status === 'error' || outcome.status === 'canceled') {
        // The failed run put the Result and its saved view back; the review goes on fitted.
        await nextTick();
        fitPreviewForSimilarityAlignment();
      }
    }
    if (outcome.status === 'error' || outcome.status === 'reviewing') {
      await nextTick();
      /** @type {HTMLElement | null} */ (document.querySelector('[data-similarity-alignment-error], [data-similarity-alignment-reset-error]'))?.focus();
    }
    return outcome;
  };
  return {
    reviewDrawingFigure, dropSimilarityAlignmentView, fitPreviewForSimilarityAlignment, endSimilarityAlignmentTakeover,
    rememberSimilarityAlignmentInvoker, restoreSimilarityAlignmentFocus, focusSimilarityAlignmentDialog,
    applySimilarityAlignmentDialog
  };
};

/**
 * @typedef {object} SimilarityAlignmentReviewPorts
 * @property {WebState} state
 * @property {{ readonly value: number }} resultEditRevision Counts the edits of the shown Result's root.
 * @property {AlignmentTakeoverRef} similarityAlignmentTakeover
 * @property {(figureShown: () => boolean) => void} holdPreviewViewport The pan/zoom owner's hold of the view.
 * @property {() => void} clearAlignmentOverlay
 * @property {(request: AlignmentOverlayRequest) => AlignmentFigureGeometry | null} showAlignmentOverlay
 * @property {() => void} closeRightDrawer
 * @property {() => unknown} closeSettingsPopup
 * @property {{ readonly value: boolean }} alignmentDialogOpen
 * @property {{ readonly value: AlignmentDraft | null }} alignmentDraft
 * @property {{ readonly value: string }} alignmentStatus
 * @property {() => boolean} isDraftArtifactCurrent Whether the figure still shows the draft's sources.
 * @property {(recordKey: string, anchor: AliasSource) => unknown} selectAlignmentCandidate
 * @property {{ readonly value: { records: readonly DirectionPreviewRecord[] } | null }} alignmentDirectionPreview
 * @property {{ readonly value: boolean }} alignmentResetDialogOpen
 * @property {() => void} openAlignmentReset
 * @property {{ readonly value: unknown }} alignmentError
 * @property {(options: { groupId?: unknown, mode?: string }) => Promise<unknown>} startAlignmentFromDrawer
 * @property {() => void} cancelAlignment
 * @property {SimilarityAlignmentTakeover['reviewDrawingFigure']} reviewDrawingFigure
 * @property {SimilarityAlignmentTakeover['dropSimilarityAlignmentView']} dropSimilarityAlignmentView
 * @property {SimilarityAlignmentTakeover['fitPreviewForSimilarityAlignment']} fitPreviewForSimilarityAlignment
 * @property {SimilarityAlignmentTakeover['endSimilarityAlignmentTakeover']} endSimilarityAlignmentTakeover
 * @property {SimilarityAlignmentTakeover['rememberSimilarityAlignmentInvoker']} rememberSimilarityAlignmentInvoker
 * @property {SimilarityAlignmentTakeover['restoreSimilarityAlignmentFocus']} restoreSimilarityAlignmentFocus
 * @property {SimilarityAlignmentTakeover['focusSimilarityAlignmentDialog']} focusSimilarityAlignmentDialog
 */

/**
 * The Review alignment panel, docked in the Result workspace. The figure overlay and the
 * panel's strips show one draft: the review passes the draft to the preview's owner, which
 * returns the figure geometry the strips are laid out from.
 *
 * Registers the review's watches, its two computeds and its keydown listener's hooks in their order;
 * the root calls it where they stood.
 * @param {SimilarityAlignmentReviewPorts} ports
 */
export const createSimilarityAlignmentReview = ({
  state, resultEditRevision, similarityAlignmentTakeover, holdPreviewViewport, clearAlignmentOverlay,
  showAlignmentOverlay, closeRightDrawer, closeSettingsPopup, alignmentDialogOpen, alignmentDraft, alignmentStatus,
  isDraftArtifactCurrent, selectAlignmentCandidate, alignmentDirectionPreview, alignmentResetDialogOpen,
  openAlignmentReset, alignmentError, startAlignmentFromDrawer, cancelAlignment, reviewDrawingFigure,
  dropSimilarityAlignmentView, fitPreviewForSimilarityAlignment, endSimilarityAlignmentTakeover,
  rememberSimilarityAlignmentInvoker, restoreSimilarityAlignmentFocus, focusSimilarityAlignmentDialog
}) => {
  const { resultPanelTab, clickedFeature, showRightDrawer, rightDrawerTab } = state;
  const similarityAlignmentCanvasHover = ref(/** @type {{ recordKey: string, candidateKey: string } | null} */ (null));
  const similarityAlignmentGeometry = shallowRef(/** @type {AlignmentFigureGeometry | null} */ (null));
  const refreshSimilarityAlignmentCanvas = () => {
    similarityAlignmentCanvasHover.value = null;
    // Read only after the checks below pass: an open dialog has a draft.
    const current = /** @type {AlignmentDraft} */ (alignmentDraft.value);
    if (!alignmentDialogOpen.value) {
      clearAlignmentOverlay();
      similarityAlignmentGeometry.value = null;
      return;
    }
    if (alignmentStatus.value !== 'reviewing') {
      // While Apply runs, the strips keep the figure they were measured from.
      clearAlignmentOverlay();
      return;
    }
    if (!isDraftArtifactCurrent()) {
      clearAlignmentOverlay();
      similarityAlignmentGeometry.value = null;
      dropSimilarityAlignmentView();
      return;
    }
    // A live rerender of the draft's records draws a new Result of the same figure (PD-OI-102).
    if (similarityAlignmentTakeover.value) similarityAlignmentTakeover.value.figure = reviewDrawingFigure(similarityAlignmentTakeover.value.drawingId);
    // Without a mounted preview (the Run info tab) nothing is measured: the strips keep the
    // geometry of the same Result.
    similarityAlignmentGeometry.value = showAlignmentOverlay({
      reference: { key: ALIGNMENT_REFERENCE_KEY, recordKey: current.reference.recordKey,
        anchor: current.response.reference, name: current.reference.name },
      records: alignmentOverlayRecords(current.rows),
      onSelect: selectAlignmentCandidate,
      onHover: (recordKey, candidateKey) => {
        similarityAlignmentCanvasHover.value = recordKey ? { recordKey, candidateKey } : null;
      }
    }) ?? similarityAlignmentGeometry.value;
  };
  watch([
    alignmentDialogOpen,
    () => alignmentDraft.value?.rows,
    alignmentStatus
  ], refreshSimilarityAlignmentCanvas, { flush: 'post' });
  // The preview is shown again after this flush's watchers have run: measure on the next tick.
  watch(resultPanelTab, () => nextTick(refreshSimilarityAlignmentCanvas));
  // An edit of the shown root (a drag, Reset Layout, an Undo of a record move)
  // leaves the Result object in place (R1), so the review measures again on it.
  watch(resultEditRevision, refreshSimilarityAlignmentCanvas, { flush: 'post' });
  const similarityAlignmentReview = computed(() => {
    const draft = alignmentDraft.value;
    const directions = alignmentDirectionPreview.value?.records;
    return draft && directions ? buildAlignmentReviewView({
      rows: draft.rows, directions, geometry: similarityAlignmentGeometry.value,
      openedWithoutChoice: draft.openedWithoutChoice
    }) : null;
  });
  const similarityAlignmentEditorDisabledReason = computed(() => (
    alignmentDialogOpen.value ? t('Finish or cancel alignment review before opening Editor.') : ''
  ));
  watch(alignmentDialogOpen, fitPreviewForSimilarityAlignment, { flush: 'post' });
  watch(alignmentDialogOpen, (/** @type {boolean} */ open) => {
    if (!open) {
      endSimilarityAlignmentTakeover();
      return;
    }
    // A review shown again while an Apply is pending keeps what the first one took.
    const reviewDrawingId = state.activeDrawingId.value;
    similarityAlignmentTakeover.value ??= { drawerOpen: showRightDrawer.value, drawerTab: rightDrawerTab.value, popup: clickedFeature.value,
      drawingId: reviewDrawingId, figure: reviewDrawingFigure(reviewDrawingId) };
    // The figure is fitted again only while it shows the draft.
    holdPreviewViewport(isDraftArtifactCurrent);
    closeRightDrawer();
    closeSettingsPopup();
  }, { flush: 'sync' });
  /** @param {KeyboardEvent} event */
  const handleSimilarityAlignmentEscape = (event) => {
    if (event.key !== 'Escape' || !alignmentDialogOpen.value) return;
    // The target of a key press is the focused element.
    if (document.querySelector('[data-linear-source-removal-dialog]')
      || /** @type {Element | null} */ (event.target)?.closest?.('.fixed.inset-0')) return;
    event.preventDefault();
    event.stopImmediatePropagation();
    void cancelSimilarityAlignmentDialog();
  };
  /** @param {Event | null | undefined} event */
  const openSimilarityAlignmentReset = (event) => {
    rememberSimilarityAlignmentInvoker(event);
    openAlignmentReset();
  };
  watch([alignmentDialogOpen, alignmentResetDialogOpen], async (/** @type {[unknown, unknown]} */ [open, resetOpen]) => {
    if (resetOpen) {
      await nextTick();
      document.getElementById('similarity-alignment-reset-title')?.focus();
      return;
    }
    // After Apply the focus returns once the outcome has given the drawer back.
    if (open || similarityAlignmentTakeover.value) return;
    void restoreSimilarityAlignmentFocus();
  });
  onMounted(() => {
    document.addEventListener('keydown', handleSimilarityAlignmentEscape, true);
  });
  onUnmounted(() => {
    document.removeEventListener('keydown', handleSimilarityAlignmentEscape, true);
  });
  /** @template T @param {T} outcome @returns {Promise<T>} */
  const finishSimilarityAlignmentStart = async (outcome) => {
    if (alignmentDialogOpen.value) {
      if (alignmentError.value) {
        await nextTick();
        /** @type {HTMLElement | null} */ (document.querySelector('[data-similarity-alignment-error]'))?.focus();
      } else await focusSimilarityAlignmentDialog();
    } else await restoreSimilarityAlignmentFocus();
    return outcome;
  };
  /** @param {string} groupId @param {Event | null} [event] @param {string} [mode] */
  const startSimilarityAlignmentFromDrawer = async (groupId, event = null, mode = 'align') => {
    rememberSimilarityAlignmentInvoker(event);
    return finishSimilarityAlignmentStart(
      await startAlignmentFromDrawer({ groupId, mode })
    );
  };
  const cancelSimilarityAlignmentDialog = async () => {
    cancelAlignment();
  };
  return {
    similarityAlignmentCanvasHover, refreshSimilarityAlignmentCanvas, similarityAlignmentReview,
    similarityAlignmentEditorDisabledReason, openSimilarityAlignmentReset, finishSimilarityAlignmentStart,
    startSimilarityAlignmentFromDrawer, cancelSimilarityAlignmentDialog
  };
};
