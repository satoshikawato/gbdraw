// @ts-check
// The steps of the Result binder that the composition root wires into the
// preview runtime (`configureMountedResultBinder`): the Legend a mounted
// Result adopts, its composition, drag, stroke and canvas bindings, the label
// editor's sync, the feature selection, and what runs once the Result is
// ready. The root keeps the delegated interactions step, which reads the
// similarity alignment port the root fills later, and the Result paint record,
// whose methods arrive as ports. The editor intent (app/editor-intent.js)
// supplies the projection the steps call.
import { COMPOSITION_METADATA_ATTRIBUTE, COMPOSITION_SCHEMA_ATTRIBUTE } from './legend-layout/composition-actions.js';
import { drawnBlockStroke, legendStructureEdited } from '../services/legend-svg.js';
import { recordStructuralMetric } from '../services/runtime-test-hooks.js';
import { normalizeOptionalHexColor } from '../utils/color-utils.js';
import { normalizeUserFacingError } from '../utils/error-normalization.js';
/** @import { DrawingState, WebState } from '../state.js' */
/** @import { SpecificColorRule } from '../services/specific-color-rules.js' */
/** @import { LiveRedrawLegend } from './legend/entry-actions.js' */
/** @import { EditorPaintState } from './result-paint-record.js' */

/**
 * The binding options the steps read (the preview runtime's
 * `PreviewBindingOptions`).
 * @typedef {object} BoundResultOptions
 * @property {boolean} [isIncrementalEdit]
 * @property {boolean} [skipLegendExtraction]
 * @property {boolean} [trustedRestore]
 * @property {boolean} [replaceGeneratedLegend]
 * @property {boolean} [modeArrival]
 * @property {readonly string[]} [requiredLabelFeatureIds]
 * @property {readonly string[]} [optionalLabelFeatureIds]
 * @property {{ featureIds: readonly string[], report: (error: unknown) => void }} [reportedLabelBinding]
 * @property {EditorPaintState} [drawnPaint]
 * @property {number} [historyRequestedAt]
 * @property {LiveRedrawLegend} [liveRedrawLegend]
 */

/**
 * What a binder step reads of one observed mount (the preview runtime's
 * `MountedResultContext`).
 * @typedef {object} BoundResultContext
 * @property {SVGSVGElement} root
 * @property {string} sourceClass
 * @property {number} resultIndex
 * @property {string} resultIdentity
 * @property {number} rootGeneration
 * @property {string} phase
 * @property {Readonly<BoundResultOptions>} bindingOptions
 */

/**
 * @typedef {object} MountedLegendBindingPorts
 * @property {WebState} state
 * @property {{ readonly value: Record<string, unknown>[] }} results The displayed drawing's Results.
 * @property {() => string[]} liveResultIdentities The Results of all drawings.
 * @property {(drawing: DrawingState, context: BoundResultContext, resultLegendOrder: string[]) => Promise<void>} projectEditorIntentOnDisplay
 * @property {(drawing: DrawingState, context: BoundResultContext) => void} rememberCommittedEditorState
 *   The editor intent's (app/editor-intent.js).
 * @property {() => boolean} isLegendLayoutReady
 * @property {() => Promise<unknown>} prepareLegendLayout
 * @property {(svg: SVGSVGElement, options?: { resultIdentity?: string, liveResultIdentities?: string[] }) => string[]} captureResultInventory
 * @property {(identity: string | undefined, options?: { restored?: boolean }) => unknown} adoptResultInventory
 * @property {(options?: { commit?: boolean }) => unknown} onLegendGeometryChanged
 * @property {(options?: { replaceGeneratedInventory?: boolean, liveRedraw?: LiveRedrawLegend | null, liveResultIdentities?: string[] }) => unknown} extractLegendEntries
 * @property {() => unknown} setupLegendDrag
 *   The Legend owner's methods.
 * @property {(result: Record<string, unknown>) => string} getResultIdentity
 * @property {(reason: string) => unknown} commitActiveResultEdit
 *   The preview runtime's methods.
 * @property {(write: () => void, request: { requestedAt: number }) => unknown} adoptDrawnWrite History's adoption of a drawn write.
 * @property {() => unknown} captureBaseConfig
 * @property {(preserveOffset?: boolean) => unknown} setupDiagramDrag
 * @property {() => unknown} applyCanvasPadding
 *   The Legend layout's methods.
 */

/**
 * Registers nothing. The root creates it after the editor intent, whose
 * projection its Legend step calls, and wires what it returns into the
 * preview runtime.
 * @param {MountedLegendBindingPorts} ports
 */
export const createMountedLegendBinding = ({
  state, results, liveResultIdentities, projectEditorIntentOnDisplay, rememberCommittedEditorState,
  isLegendLayoutReady, prepareLegendLayout, captureResultInventory, adoptResultInventory, onLegendGeometryChanged,
  extractLegendEntries, setupLegendDrag, getResultIdentity, commitActiveResultEdit, adoptDrawnWrite, captureBaseConfig,
  setupDiagramDrag, applyCanvasPadding
}) => {
  // A legacy imported SVG without composition metadata stays unbound: it has
  // no composition to capture and no canvas to pad.
  // Mounted Results whose Legend `adoptLegend` laid out; the commit waits for
  // the later binding steps (`initializeStrokeAndCanvas`).
  /** @type {WeakSet<Element>} */
  const mountedLegendLayouts = new WeakSet();
  /** @param {BoundResultContext} context */
  const shouldBindComposition = (context) => (
    context.root.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) !== null
    || context.root.getAttribute(COMPOSITION_METADATA_ATTRIBUTE) !== null
    || (!context.bindingOptions.isIncrementalEdit && context.sourceClass !== 'legacy-import')
  );
  return {
    /** @param {BoundResultContext} context */
    async adoptLegend(context) {
      const drawing = state.activeDrawing();
      // The Legend layout port (zero shift) reads the bundled-font metrics;
      // they load with the first Result that has a Legend, so the Legend
      // edits, which lay the Legend out synchronously, find them loaded.
      // Once loaded, the binder does not wait: an edit waiting for this
      // binding (a rule commit after a rerender) reads its Legend unchanged.
      if (context.root.getElementById?.('legend') && !isLegendLayoutReady()) {
        await prepareLegendLayout().catch((error) => {
          console.error('The Legend layout could not load its font metrics.', normalizeUserFacingError(error));
        });
      }
      // OV-47: each Result has its own default Legend order. A Result being
      // displayed is read before its editor intent is projected; it then shows
      // its own inventory, so a Generate or rerender made while another Result
      // was displayed does not replace it.
      // A Result shown again by a drawing switch is a selection in its own
      // drawing, whose Legend edits it was last shown with (PD-OI-090).
      const selecting = context.phase === 'result-selection' && !context.bindingOptions.trustedRestore;
      if (selecting) {
        const resultLegendOrder = captureResultInventory(context.root, {
          resultIdentity: context.resultIdentity,
          liveResultIdentities: liveResultIdentities()
        });
        // A Result shown again by a mode switch arrives in its own drawing,
        // whose Legend edits name that Result's inventory: the inventory is
        // adopted first, so the departing mode's order is never compared with
        // the arriving drawing's rows. A batch Result of the same drawing is
        // compared with the order shown until now (B20, OV-47).
        const arriving = Boolean(context.bindingOptions.modeArrival);
        if (arriving) adoptResultInventory(context.resultIdentity);
        await projectEditorIntentOnDisplay(drawing, context, resultLegendOrder);
        if (!arriving) adoptResultInventory(context.resultIdentity);
      } else {
        rememberCommittedEditorState(drawing, context);
      }
      if (context.bindingOptions.trustedRestore) {
        adoptResultInventory(context.resultIdentity, { restored: true });
        return;
      }
      // A Result the renderer drew: generated, newly displayed, or drawn again
      // by an automatic rerender. An incremental edit, a History restore, and
      // a loaded Session show bytes already laid out.
      const drawn = !context.bindingOptions.isIncrementalEdit
        || Boolean(context.bindingOptions.replaceGeneratedLegend);
      const adopt = () => {
        // It shows the Legend editor's edits laid out as Python lays the edited
        // rows out, before its entries are read (zero shift; OV-122, OV-124,
        // OV-126, OV-127): a Result whose records show a Legend structure edit.
        if (
          drawn && context.root.getAttribute(COMPOSITION_METADATA_ATTRIBUTE) !== null
          && legendStructureEdited(context.root)
        ) {
          onLegendGeometryChanged({ commit: false });
          mountedLegendLayouts.add(context.root);
        }
        if (context.bindingOptions.skipLegendExtraction) return;
        recordStructuralMetric('legendDomFullScanCount', 1, {
          phase: context.phase,
          rootGeneration: context.rootGeneration
        });
        // A draw that replays an edited Legend order replays it on the displayed
        // mode's Results only.
        extractLegendEntries({
          replaceGeneratedInventory: !selecting && drawn,
          liveRedraw: context.bindingOptions.liveRedrawLegend ?? null,
          liveResultIdentities: results.value.map(getResultIdentity)
        });
      };
      // OV-872: the Legend a label rerender draws is no edit of a step begun
      // after it was requested (an action waiting for it, D-40).
      const { historyRequestedAt } = context.bindingOptions;
      if (historyRequestedAt === undefined) adopt();
      else adoptDrawnWrite(adopt, { requestedAt: historyRequestedAt });
    },
    /** @param {BoundResultContext} context */
    bindComposition(context) {
      if (context.bindingOptions.trustedRestore) return;
      if (shouldBindComposition(context)) captureBaseConfig();
    },
    /** @param {BoundResultContext} context */
    setupDragAffordances(context) {
      setupLegendDrag();
      setupDiagramDrag(Boolean(context.bindingOptions.isIncrementalEdit));
    },
    /** @param {BoundResultContext} context */
    initializeStrokeAndCanvas(context) {
      const legendLaidOut = mountedLegendLayouts.delete(context.root);
      if (!context.bindingOptions.trustedRestore && !context.bindingOptions.isIncrementalEdit) {
        // The block stroke Python drew on the shown Result: the Result's
        // block-stroke default, which the Session format keeps as
        // `originalSvgStroke`.
        const blockStroke = drawnBlockStroke(context.root);
        if (blockStroke) state.originalSvgStroke.value = { ...blockStroke, color: normalizeOptionalHexColor(blockStroke.color) };
        // Generate already padded its candidates; another batch Result shows
        // the current canvas padding when it is displayed (D-09).
        if (shouldBindComposition(context)) applyCanvasPadding();
      }
      // The Legend laid out at mount is committed with the bindings the
      // steps since added, so the Result's content is its mounted SVG (R1).
      if (legendLaidOut) commitActiveResultEdit('legend-position');
    }
  };
};

/**
 * @typedef {object} MountedResultReadinessPorts
 * @property {WebState} state
 * @property {{ readonly value: boolean }} isFeatureDrawerMounted
 * @property {{ readonly value: unknown }} clickedFeature The open feature popup's feature, or null.
 * @property {{ readonly show: boolean }} labelTextScopeDialog
 * @property {{ readonly value: boolean }} labelReflowProcessing Whether a label rerender runs.
 * @property {{ value: number }} mountedResultRevision Counts the Results made ready.
 * @property {() => string} takeLabelProjectionResultIdentity
 * @property {() => { identity: string, rules: SpecificColorRule[] | null, labelIntentGone: boolean } | null} takeDisplayAsks
 * @property {(drawing: DrawingState) => EditorPaintState} currentEditorProjectionState
 * @property {(context: BoundResultContext) => void} showEditsMadeWhileDrawn
 *   The editor intent's (app/editor-intent.js).
 * @property {(options: { requiredFeatureIds?: readonly string[], optionalFeatureIds?: readonly string[], reportedLabelBinding?: { featureIds: readonly string[], report: (error: unknown) => void } }) => unknown} syncLabelEditor
 * @property {(rules: SpecificColorRule[]) => unknown} followRestoredSpecificRules
 * @property {() => unknown} requestAutomaticRerender
 * @property {(drawing: DrawingState) => void} redrawShown
 *   The feature editor's methods.
 * @property {(options?: { clearStatus?: boolean, syncDom?: boolean }) => unknown} clearFeatureSelection The feature selection's clear.
 * @property {(drawn: EditorPaintState | null | undefined, current: EditorPaintState) => boolean} resultPaintDrewPaint
 *   Whether a label rerender shows the current paint (Result paint record).
 * @property {() => unknown} handleMountedResultReady The preview feature search's reader of a ready Result.
 * @property {() => Promise<unknown>} keepRedrawnResult History keeps the redrawn figure.
 */

/**
 * Registers nothing. The root creates it after the editor intent and wires
 * what it returns into the preview runtime.
 * @param {MountedResultReadinessPorts} ports
 */
export const createMountedResultReadiness = ({
  state, isFeatureDrawerMounted, clickedFeature, labelTextScopeDialog, labelReflowProcessing, mountedResultRevision,
  takeLabelProjectionResultIdentity, takeDisplayAsks, currentEditorProjectionState, showEditsMadeWhileDrawn,
  syncLabelEditor, followRestoredSpecificRules, requestAutomaticRerender, redrawShown, clearFeatureSelection,
  resultPaintDrewPaint, handleMountedResultReady, keepRedrawnResult
}) => {
  /** @param {DrawingState} drawing */
  const shouldSyncMountedLabelEditor = (drawing) => (
    isFeatureDrawerMounted.value
    || Boolean(clickedFeature.value)
    || labelTextScopeDialog.show
    || Object.values(drawing.featureOverrides).some((row) => row.labelText !== null || row.labelVisibility !== null)
    || Object.keys(drawing.labelTextBulkOverrides).length > 0
  );
  return {
    /** @param {BoundResultContext} context */
    synchronizeLabelEditor(context) {
      const drawing = state.activeDrawing();
      const labelsChanged = takeLabelProjectionResultIdentity() === context.resultIdentity;
      if (!context.bindingOptions.trustedRestore && (labelsChanged || shouldSyncMountedLabelEditor(drawing))) {
        syncLabelEditor({
          requiredFeatureIds: context.bindingOptions.requiredLabelFeatureIds,
          optionalFeatureIds: context.bindingOptions.optionalLabelFeatureIds,
          reportedLabelBinding: context.bindingOptions.reportedLabelBinding
        });
      }
    },
    /** @param {BoundResultContext} context */
    reconcileSelection(context) {
      if (
        !context.bindingOptions.trustedRestore
        && !context.bindingOptions.isIncrementalEdit
      ) {
        clearFeatureSelection({ clearStatus: true, syncDom: false });
      }
    },
    /** @param {BoundResultContext} context */
    afterReady(context) {
      mountedResultRevision.value += 1;
      // OV-1280: a Result a mode switch shows again, whose drawing changed its
      // rules while it was not shown (a file card's Remove, an Undo or Redo
      // there), is drawn again when its Legend cannot show them, as after a
      // History restore of the rules. OV-2293: any Result displayed is drawn
      // again when its bytes hold a label edit that left its drawing while it
      // was not shown (a Reset Settings, an Undo). The display records what it
      // shows, so a failed rerender does not repeat here.
      const display = takeDisplayAsks();
      if (display?.identity === context.resultIdentity) {
        if (context.bindingOptions.modeArrival && display.rules) followRestoredSpecificRules(display.rules);
        if (display.labelIntentGone) requestAutomaticRerender();
      }
      // D-36: a label rerender pays the redraw its drawing owes when, as its
      // Result is ready, no later rerender runs and its figure shows the
      // drawing's current paint, also after an edit made since its paint was
      // committed (read before that edit is shown on it, OV-346).
      // Its drawing is the shown one: a mode switch waits for a rerender
      // (`drawingSwitchBusy`). D-24 (OV-731): History keeps the figure
      // as the last step's.
      const drawing = state.activeDrawing();
      const drewCurrent = context.phase === 'label-reflow' && !labelReflowProcessing.value
        && resultPaintDrewPaint(context.bindingOptions.drawnPaint, currentEditorProjectionState(drawing));
      handleMountedResultReady();
      showEditsMadeWhileDrawn(context);
      if (!drewCurrent) return;
      redrawShown(drawing);
      keepRedrawnResult().catch((error) => {
        console.error('History could not keep the redrawn figure.', normalizeUserFacingError(error));
      });
    }
  };
};
