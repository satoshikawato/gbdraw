// @ts-check
// The drawing switch (R15): the one switch between the drawings
// (`selectDrawing`), its work without the gate (`showDrawingNow`, which the
// drawing list's steps call), why a switch waits (`drawingSwitchBusy`), the
// drawing an operation error names, and the swap of the artifact slots with
// the install of the arriving slot, which Session Load also uses for the slot
// it opens on (`installArrivingSlot`). The composition root passes the state, its
// refs, the members of History, the artifact slots, the preview runtime, the
// Legend and feature editors, the feature selection, the Preview search and
// the notice banner, and the install of the shown drawing
// (`installShownDrawing`) as ports.
import { normalizeUserFacingError } from '../utils/error-normalization.js';
/** @import { UserFacingError } from '../utils/error-normalization.js' */
import { msg } from '../utils/i18n.js';
import { createArtifactSlot } from '../services/artifact-slot.js';
/** @import { DrawingState, SessionBusy, WebState } from '../state.js' */
/** @import { DrawingId, DrawingIdentity } from '../services/drawing-view.js' */
/** @import { ArtifactSlot } from '../services/artifact-slot.js' */
/** @import { PythonLegendKey } from '../services/legend-svg.js' */
/** @import { SvgResultEntry } from '../services/svg-result-ingestion.js' */

const { computed } = window.Vue;

const HISTORY_RESTORE_BUSY = Object.freeze({ status: 'busy', reason: msg('Undo or Redo in progress. Retry after it finishes.') });
const RULE_WRITE_BUSY = Object.freeze({ status: 'busy', reason: 'A color rule edit is being matched. Retry after it finishes.' });

/**
 * Why History cannot take a step now, as History words it, or null.
 * @typedef {{ readonly status: string, readonly reason: string }} HistoryBusy
 */

/**
 * @typedef {object} DrawingSwitchPorts
 * @property {WebState} state
 * @property {{ readonly value: SvgResultEntry[] }} results The displayed drawing's Results.
 * @property {{ readonly value: number | string }} selectedResultIndex
 * @property {WebState['errorLog']} errorLog The operation error shown, or null.
 * @property {{ value: boolean }} failedGeneratePreservedResult Whether "Showing the last successful result" shows.
 * @property {(operation?: string) => SessionBusy | null} sessionOperationAvailability
 * @property {() => HistoryBusy | null} historyAvailability Why History's open step owns the state, or null.
 * @property {() => boolean} historyRestoring Whether the shown drawing's History restores a step now.
 * @property {() => boolean} traversalPending Whether an Undo or Redo waits to run.
 * @property {() => Promise<unknown>} settlePendingIntent Ends the open History step of the shown drawing.
 * @property {() => Promise<unknown>} refreshCurrentIntent History's capture of the shown drawing's intent.
 * @property {(reset: () => void) => void} adoptSessionReset Resets Session-level transient values without a step.
 * @property {() => boolean} ruleWritePending Whether a color rule table edit is being matched.
 * @property {() => unknown} suspendSpecificRulePatternDrafts
 * @property {(options?: { clearStatus?: boolean, syncDom?: boolean }) => unknown} clearFeatureSelection
 * @property {(drawing: DrawingState) => void} rememberDepartingResultProjection
 * @property {(departing: DrawingState, arriving: DrawingState) => UserFacingError | null} auxiliaryImportFailureOnArrival
 *   The arriving drawing's card failure error, or null.
 * @property {() => void} resetModeTransientUi
 * @property {(id: DrawingId) => DrawingState} installShownDrawing The one writer of the shown drawing.
 * @property {() => void} dismissNotice
 * @property {(departing: DrawingId, arriving: DrawingId) => void} switchDrawingSearch
 * @property {() => void} rebindAnnotationRecords Rebinds the shown drawing's annotations to its record catalog (R10).
 * @property {() => void} clearActiveRuntime Releases the displayed Result.
 * @property {(result: SvgResultEntry) => string} getResultIdentity
 * @property {(bindingOptions?: { modeArrival?: boolean }) => unknown} presentSelectedResult
 * @property {(identity: string | undefined, inventory: readonly PythonLegendKey[]) => void} rememberResultInventory
 * @property {(departing: Readonly<ArtifactSlot>, arrivingId: DrawingId) => Readonly<ArtifactSlot> | null} artifactSlotsSwap
 *   Keeps the departing slot and takes the arriving drawing's, or null.
 * @property {() => Readonly<ArtifactSlot>} captureArtifactSlot The displayed artifact as a slot.
 * @property {(slot: Readonly<ArtifactSlot>, options: { drawing: DrawingIdentity }) => void} installArtifactSlot
 *   Installs the slot, with its Legend inventory and colors.
 */

/**
 * Registers one computed (`drawingSwitchAvailable`). The root creates it
 * before the drawing list, whose steps and busy reason read what it returns.
 * @param {DrawingSwitchPorts} ports
 */
export const createDrawingSwitch = ({
  state, results, selectedResultIndex, errorLog, failedGeneratePreservedResult, sessionOperationAvailability,
  historyAvailability, historyRestoring, traversalPending, settlePendingIntent, refreshCurrentIntent, adoptSessionReset,
  ruleWritePending, suspendSpecificRulePatternDrafts, clearFeatureSelection, rememberDepartingResultProjection,
  auxiliaryImportFailureOnArrival, resetModeTransientUi, installShownDrawing, dismissNotice, switchDrawingSearch,
  rebindAnnotationRecords, clearActiveRuntime, getResultIdentity, presentSelectedResult, rememberResultInventory,
  artifactSlotsSwap, captureArtifactSlot, installArtifactSlot
}) => {
  // E1: installs the arriving slot as the displayed artifact, with its Legend
  // inventory and colors. The Legend owner also keeps each Result's inventory
  // (OV-47): the arriving Result's goes to it, and the binder adopts it when
  // the Result is displayed, as for a batch Result.
  /** @param {Readonly<ArtifactSlot>} slot */
  const installArrivingSlot = (slot) => {
    installArtifactSlot(slot, { drawing: state.drawingOf(slot.drawingId) });
    const shown = results.value[Number(selectedResultIndex.value) || 0];
    if (shown) rememberResultInventory(getResultIdentity(shown), [...slot.legend.originalOrder]);
  };
  /**
   * Shows the artifact of `arrivingDrawing`: keeps the displayed one under the
   * drawing it belongs to, which is the arriving drawing when the shown
   * drawing displayed the arriving one's, and installs the arriving drawing's
   * slot, or its empty slot when it has no artifact yet.
   * @param {DrawingIdentity} arrivingDrawing
   */
  const swapArtifactSlots = (arrivingDrawing) => {
    clearActiveRuntime();
    installArrivingSlot(artifactSlotsSwap(captureArtifactSlot(), arrivingDrawing.id)
      ?? createArtifactSlot({ drawing: arrivingDrawing, values: {} }));
  };

  // PD-OI-090 (R10, R15): the one switch between the drawings, in one task.
  // It writes no setting and records no History step: each drawing has its
  // own History, and Undo and Redo walk the shown drawing's steps. No watcher
  // writes state on a switch. Steps, each one call:
  //   1. Settle the departing drawing: its open History step ends here, in its
  //      own History; pending rule-pattern work; the feature selection; and the
  //      displayed Result's projection state (it followed every live edit).
  //   2. Swap the artifact slots (`swapArtifactSlots`).
  //   3. Show the arriving drawing (`installShownDrawing`); the template now
  //      binds it. Its History takes the shared selections in its intent
  //      that changed while it was not shown.
  //   4. Reset the departing drawing's transient UI; the Preview search is
  //      kept per drawing. The feature panel tab it resets is no Session
  //      change (D-15-53).
  //   5. Show the arriving Result as a selection.
  // A switch waits for Generate, a label rerender, Save and Load, so each
  // Result lands in its own drawing's slot (OIPC-C07), and for a History
  // restore or an open History transaction.
  // OIC-030: also while a color rule table edit is matched, so that edit is
  // one step of its drawing.
  const drawingSwitchBusy = () => sessionOperationAvailability('history')
    || historyAvailability()
    || (historyRestoring() || traversalPending() ? HISTORY_RESTORE_BUSY : null)
    || (ruleWritePending() ? RULE_WRITE_BUSY : null);
  const drawingSwitchAvailable = computed(() => !drawingSwitchBusy());
  // The drawing an operation error names (PD-OI-079 `context.drawingId`), if
  // any: read by the switch, the drawing removal, and E1's Generate. Its
  // producers write a `DrawingId`; the normalizer, a leaf, keeps it as text.
  const errorDrawingId = () => /** @type {DrawingId | undefined} */ (
    normalizeUserFacingError(errorLog.value)?.context.drawingId);
  // The switch work, steps 1-5, without the gate: `selectDrawing` and a step
  // of the drawing list (`app/drawing-list.js`, which shows a drawing inside
  // its project step, R16-2) call it. Showing the shown drawing does nothing.
  /** @param {DrawingId} drawingId */
  const showDrawingNow = (drawingId) => {
    const arriving = state.drawingOf(drawingId);
    const departing = state.activeDrawing();
    if (arriving === departing) return;
    // 1. Settle the departing drawing; its capture runs in this call.
    void settlePendingIntent();
    // A card's failed import is its drawing's card state (D-15-126): read
    // before the departing drawing's message is cleared below.
    const arrivingCardFailure = auxiliaryImportFailureOnArrival(departing, arriving);
    suspendSpecificRulePatternDrafts();
    clearFeatureSelection({ clearStatus: true, syncDom: false });
    rememberDepartingResultProjection(departing);
    // 2. Artifact slots.
    swapArtifactSlots(arriving);
    // 3. The shown drawing.
    installShownDrawing(arriving.id);
    void refreshCurrentIntent();
    // 4. Transient UI. "Showing the last successful result" named the
    // departing drawing's Result, and the notice the departing drawing's work.
    adoptSessionReset(resetModeTransientUi);
    dismissNotice();
    switchDrawingSearch(departing.id, arriving.id);
    failedGeneratePreservedResult.value = false;
    // OV-372: an operation error belongs to the drawing it was raised in,
    // except one that names the arriving drawing's Result (PD-OI-079).
    if (errorLog.value && errorDrawingId() !== arriving.id) errorLog.value = null;
    if (arrivingCardFailure !== null && errorLog.value === null) errorLog.value = arrivingCardFailure;
    // The arriving drawing's annotations, over its mode's record catalog.
    rebindAnnotationRecords();
    // 5. Presentation: the arriving Result is shown as a selection.
    presentSelectedResult({ modeArrival: true });
  };
  // `false` when the drawing is already shown, the busy reason when a switch
  // waits, else `{ status: 'ok' }`. An unknown ID throws (`state.drawingOf`):
  // the tabs and the drawing menu offer only existing drawings (OV-2151).
  /** @param {DrawingId} drawingId */
  const selectDrawing = (drawingId) => {
    const arriving = state.drawingOf(drawingId);
    if (arriving === state.activeDrawing()) return false;
    const busy = drawingSwitchBusy();
    if (busy) return busy;
    showDrawingNow(arriving.id);
    return { status: 'ok' };
  };

  return { installArrivingSlot, drawingSwitchBusy, drawingSwitchAvailable, errorDrawingId, showDrawingNow, selectDrawing };
};
