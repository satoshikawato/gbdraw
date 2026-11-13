// @ts-check
// Reset Settings and Reset Layout: the reset of one drawing to a new
// drawing's settings of its mode with what its mode resets beyond its slices
// (`resetDrawing`), the question Reset Settings asks
// (`requestResetSettings`), Reset this drawing and Reset all drawings, each one
// History step that shows the reset draft on the displayed Result, the
// Instant Preview toggle's step, and Reset Layout. The settings reset itself
// is services/reset.js's. The composition root passes the state, its refs,
// the Linear File edit, the editor intent projection, and the members of
// History, the artifact store, the Legend and feature editors, the settings
// pane, the notice banner, the Results owner and the match sequences as
// ports. It keeps whether Reset Layout has something to reset
// (`layoutResetAvailable`).
import { LIVE_EDIT_DOMAINS } from './candidate-render.js';
import { buildDefaultCircularRecordPositions } from './circular-record-panel.js';
import { ensureDepthTrackConfigCount, nameDepthSeriesAfterFiles } from './depth-tracks.js';
import { paletteColorsEqual } from './svg-styles.js';
import { resetLayoutState, resetProject, resetSettings as resetDrawingSettings } from '../services/reset.js';
import { msg, t } from '../utils/i18n.js';
import { wholeRecordCardWindows } from '../services/record-draw-selection.js';
/** @import { DiagramMode, DrawingState, SessionBusy, WebState } from '../state.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { CheckpointHold } from '../services/history.js' */
/** @import { AskChoice } from './app-messages.js' */
/** @import { ApplyLinearSeqMutation } from './linear-source-edits.js' */

/**
 * @typedef {object} DrawingResetPorts
 * @property {WebState} state
 * @property {{ value: string }} circularTrackNewRenderer The Circular track menu's new-track choice.
 * @property {{ value: string }} linearTrackNewRenderer The Linear track menu's new-track choice.
 * @property {(drawing: DrawingState) => void} rederiveCircularConservation
 *   Derives the comparison rings again from the kept comparison files, without the LOSAT result of the old settings.
 * @property {ApplyLinearSeqMutation} applyLinearSeqMutation
 * @property {(drawing: DrawingState, options?: { preserveLosatCacheInfo?: boolean }) => void} invalidateLinearComparisonArtifacts
 * @property {(drawingId: DrawingId, values: Readonly<Record<string, unknown>>) => void} writeArtifactValues
 *   Writes slot values of a drawing: into the displayed artifact while it is shown, else into its slot.
 */

/**
 * Registers nothing. The root creates it once the Linear File edits exist.
 * @param {DrawingResetPorts} ports
 */
export const createDrawingReset = ({
  state, circularTrackNewRenderer, linearTrackNewRenderer, rederiveCircularConservation, applyLinearSeqMutation,
  invalidateLinearComparisonArtifacts, writeArtifactValues
}) => {
  // What a drawing's reset changes beyond its slices, by its mode.
  /** @type {Readonly<Record<DiagramMode, (drawing: DrawingState) => void>>} */
  const RESET_FOLLOW_UPS = Object.freeze({
    circular: (drawing) => {
      // The track menu's new-track choice is one per mode; it returns to a
      // new drawing's with the drawing.
      circularTrackNewRenderer.value = 'dinucleotide_skew';
      // OV-1249: the loaded Circular records keep a Record Order row each, at
      // their defaults; Reset resets settings, not the records.
      drawing.adv.multi_record_positions.push(...buildDefaultCircularRecordPositions(drawing));
      // OV-1167 (Q-DR-7 B): the comparison rings are derived again from the
      // kept ring files, as a new drawing with those files uploaded has them;
      // a ring search of the old settings goes with them (OV-2758).
      rederiveCircularConservation(drawing);
    },
    linear: (drawing) => {
      linearTrackNewRenderer.value = 'dinucleotide_skew';
      // Linear records return to their File defaults and inferred
      // definitions; Files, record selections, File defaults, and depth stay.
      // The mutation also clears the alignment plan (D-15).
      applyLinearSeqMutation(drawing, drawing.linearSeqs.map((sequence) => ({
        ...sequence,
        definition: '',
        record_subtitle: '',
        ...wholeRecordCardWindows()
      })), { alignmentMutation: msg('settings reset.') });
      invalidateLinearComparisonArtifacts(drawing);
    }
  });
  // One drawing returns to a new drawing's settings of its mode; the
  // one-per-Session values stay (R17-5). OV-1168, OV-1195: its artifact takes
  // the default palette applied, and a drawing not shown loses its alignment
  // Reset receipt (the shown one's plan is cleared by its follow-up).
  /** @param {DrawingState} drawing */
  const resetDrawing = (drawing) => {
    const artifact = resetDrawingSettings(state, drawing);
    writeArtifactValues(drawing.id, drawing.id === state.activeDrawingId.value
      ? artifact
      : { ...artifact, similarityAlignmentResetReceipt: null });
    RESET_FOLLOW_UPS[drawing.mode](drawing);
    // The reset Depth series fit the kept Depth files (OV-1164), and each
    // kept file names its series (OV-2822).
    ensureDepthTrackConfigCount(drawing);
    nameDepthSeriesAfterFiles(drawing);
  };

  return { resetDrawing };
};

/**
 * @typedef {(label: string, fn: () => unknown) => unknown} RunResetStep One History step that runs `fn`.
 * @typedef {object} ResetSettingsPorts
 * @property {WebState} state
 * @property {(operation?: string) => SessionBusy | null} sessionOperationAvailability
 * @property {AskChoice} askChoice
 * @property {(drawing: DrawingState) => void} resetDrawing One drawing's reset (`createDrawingReset`).
 * @property {CheckpointHold} holdAutomaticRerender
 *   Holds the automatic rerender requests until `run` ends; one rerender follows.
 * @property {<T>(drawing: DrawingState, write: () => Promise<T>) => Promise<T>} followRuleEdits
 *   Runs `write`; the rule owner follows the rules it changed.
 * @property {{ readonly value: Element | null }} svgContainer The element that holds the mounted SVG.
 * @property {{ readonly value: Record<string, string> }} appliedPaletteColors The palette the displayed Result shows.
 * @property {{ readonly value: string[] }} originalLegendOrder The captions of the rows Python drew.
 * @property {(options: { colors?: boolean, visibility?: boolean, rerender?: boolean, labels?: boolean,
 *   domains?: readonly string[] }) => Promise<string[] | false>} projectMountedEditorIntent
 * @property {(keys: string[]) => void} drawReturningLegendRows
 * @property {readonly string[]} STROKE_DOMAINS The paint domains a stroke edit shows.
 * @property {() => void} resetAllPositions
 * @property {() => void} resetCanvasPadding
 * @property {(options?: { resetZoom?: boolean }) => void} resetPreviewViewport
 * @property {() => void} dismissNotice
 * @property {() => void} closeSettingsPopup
 * @property {() => void} clearSpecificRulePatternDrafts
 * @property {() => unknown} extractLegendEntries The Legend list follows the displayed Result.
 * @property {RunResetStep} runProjectStep A project step: one Undo from any drawing (R11).
 * @property {RunResetStep} runUndoable A step of the shown drawing.
 * @property {(enabled: boolean) => void} storePaletteInstantPreview The Results owner's Instant Preview setter.
 */

/**
 * Registers nothing. The root creates it after `createDrawingReset`; only the
 * root's return reads what it returns.
 * @param {ResetSettingsPorts} ports
 */
export const createResetSettings = ({
  state, sessionOperationAvailability, askChoice, resetDrawing, holdAutomaticRerender, followRuleEdits, svgContainer,
  appliedPaletteColors, originalLegendOrder, projectMountedEditorIntent, drawReturningLegendRows, STROKE_DOMAINS,
  resetAllPositions, resetCanvasPadding, resetPreviewViewport, dismissNotice, closeSettingsPopup,
  clearSpecificRulePatternDrafts, extractLegendEntries, runProjectStep, runUndoable,
  storePaletteInstantPreview
}) => {
  // Reset Settings asks first (PD-OI-090 〔Reset〕): Reset this drawing resets
  // the shown drawing as its step, Reset all drawings every drawing and the
  // one-per-Session values as one project step (R17-5), and Cancel changes
  // nothing and records no step. The body is D-04's two lines and one line
  // for the one-per-Session values only Reset all drawings resets (D-15-191).
  // `resetShownDrawing` and `resetSettings` do not ask.
  const requestResetSettings = async () => {
    const sessionBusy = sessionOperationAvailability('not-during-run');
    if (sessionBusy) return sessionBusy;
    const choice = await askChoice({
      title: t('Reset settings?'),
      lines: [
        { icon: 'ph-warning', text: t('Your edits are removed too') },
        { icon: 'ph-arrow-counter-clockwise', text: t('Undo restores them') },
        { icon: 'ph-stack', text: t('Reset all drawings also resets the settings every drawing shares') }
      ],
      actions: [
        { id: 'reset-drawing', label: t('Reset this drawing'), tone: 'danger' },
        { id: 'reset', label: t('Reset all drawings'), tone: 'danger' }
      ]
    });
    if (choice === 'reset-drawing') return resetShownDrawing();
    if (choice === 'reset') return resetSettings();
    return false;
  };
  // Both Reset buttons reset the shown drawing: `resetDrawings` resets it (and
  // for Reset all drawings, every other drawing and the Session values) inside
  // the step `runStep` opens, and the shown Result shows the reset draft.
  /**
   * @param {(reset: () => Promise<boolean>) => unknown} runStep The History step
   * @param {() => void} resetDrawings
   */
  const runReset = (runStep, resetDrawings) => {
    const sessionBusy = sessionOperationAvailability('not-during-run');
    if (sessionBusy) return sessionBusy;
    // A notice (a loaded Session's warnings, a table import) described the
    // settings the reset replaces. Generate keeps it: it describes the Session.
    dismissNotice();
    closeSettingsPopup();

    // OV-287, OV-289, OV-290 (PD-OI-066): the displayed Result shows the
    // reset draft as the next Generate draws it. One projection of the editor
    // intent shows the fills, Legend colors, visibility and labels, and the
    // executor shows the strokes and the Legend rows the reset intent leaves
    // (a deleted row Python drew returns, and one the Result's bytes lack asks
    // for the rerender, O-2; a row added in the editor and deleted stays
    // gone); the rule owner follows the removed rules (a
    // changed Legend source asks for the rerender, as their Undo does), and
    // the Legend list follows the Result. The reset itself runs at once, so a
    // caller reads the reset settings right after resetSettings(). One
    // rerender follows the step: the rerender requests wait until it ends
    // (OV-1076).
    return holdAutomaticRerender(() => runStep(() => followRuleEdits(state.activeDrawing(), async () => {
      clearSpecificRulePatternDrafts();
      const shown = Boolean(svgContainer.value?.querySelector?.('svg'));
      const palette = { ...appliedPaletteColors.value };
      // The deleted rows Generate drew return (a row added in the editor and
      // deleted stays gone).
      const returning = state.activeDrawing().deletedLegendEntries.value
        .map((/** @type {{ caption?: string, originalCaption?: string }} */ entry) => String(entry.originalCaption || entry.caption || ''))
        .filter((/** @type {string} */ key) => originalLegendOrder.value.includes(key));
      resetDrawings();
      if (shown) {
        // A Reset that changes the palette leaves the fills to the palette watcher.
        await projectMountedEditorIntent({
          colors: paletteColorsEqual(appliedPaletteColors.value, palette),
          visibility: true,
          rerender: true,
          labels: true,
          domains: [...STROKE_DOMAINS, ...LIVE_EDIT_DOMAINS.legendRows]
        });
        drawReturningLegendRows(returning);
        extractLegendEntries();
      }
      return true;
    })));
  };
  // Reset all drawings is one project step (D-15-51): one Undo restores every
  // drawing and the Session values. OV-1166 (Q-DR-6 A): every drawing's
  // multi-record layout takes a new drawing's defaults (its sides unset),
  // whatever the canvas was before.
  const resetSettings = () => runReset((reset) => runProjectStep(msg('Reset settings'), reset), () => {
    resetProject(state);
    state.drawingList().forEach((drawing) => resetDrawing(drawing));
  });
  // Reset this drawing is a step of the shown drawing (R17-5): the other
  // drawings, their Histories, and the Session values stay.
  const resetShownDrawing = () => runReset((reset) => runUndoable(msg('Reset this drawing'), reset),
    () => resetDrawing(state.activeDrawing()));

  // OV-1160 (Owner R15-23 Q-DR-5 A): turned on, the Instant Preview toggle
  // applies every drawing's queued palette, so it is one project step whose
  // parts hold each drawing's palette (D-15-51); turned off, it changes the
  // Session-level flag alone.
  /** @param {boolean} enabled */
  const setPaletteInstantPreview = (enabled) => (enabled
    ? runProjectStep('Change setting', () => storePaletteInstantPreview(true))
    : runUndoable('Change setting', () => storePaletteInstantPreview(false)));

  const resetLayout = () => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    resetAllPositions();
    resetCanvasPadding();
    resetLayoutState(state);
    resetPreviewViewport({ resetZoom: true });
  };

  return { requestResetSettings, resetSettings, resetShownDrawing, setPaletteInstantPreview, resetLayout };
};
