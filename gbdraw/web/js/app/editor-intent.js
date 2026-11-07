// @ts-check
// The editor intent shown on the displayed Result (design §5.4, R3): the
// compile and show of the active drawing's paint and Legend edits
// (`showEditorIntent`), the live edit transition (`editEditorIntent`), the
// projection of the canonical editor intent onto the mounted Result
// (`projectMountedEditorIntent`) with the Legend restore that returns rows,
// and what each Result showed: the intent a departing Result left with, the
// projection onto a Result when it is displayed, and the display of fills and
// visibility a blocker declined. The composition root keeps the History apply
// registration (`setAfterApplyHistoryIntent`), which calls the projection, and
// the Result paint record, whose methods arrive as ports. The Result binder
// (app/mounted-result-binder.js) takes the two values kept for its steps
// through `takeLabelProjectionResultIdentity` and `takeDisplayAsks`.
import { compileDirectEditorMutationPlan, LIVE_EDIT_DOMAINS } from './candidate-render.js';
import { drawnLabelIntentGone, editorPaintState, paintRuleTable, sameColors } from './result-paint-record.js';
import { admitFeatureCatalog, biologicalFeatureKey, resultCatalogFeatures } from '../services/feature-catalog.js';
import { displayedFeatureAddressing } from '../services/feature-override-identity.js';
import { displayedLegendRowContext, featureDrawnContext } from '../services/feature-visibility.js';
import { drawsPythonLegendRow, legendStructureEdited } from '../services/legend-svg.js';
import { recordStructuralMetric } from '../services/runtime-test-hooks.js';
import { appliedFeatureColors } from '../utils/color-utils.js';
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { toRaw } from '../utils/vue-raw.js';
/** @import { DrawingState, WebState } from '../state.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { KeptLegendRow, LegendRow } from '../services/drawing/paint.js' */
/** @import { SvgMutationOperations } from '../services/svg-result-ingestion.js' */
/** @import { DisplayedLegendRowContext } from '../services/feature-visibility.js' */
/** @import { SpecificColorRule } from '../services/specific-color-rules.js' */
/** @import { UnvalidatedFields } from '../services/session-authority.js' */
/** @import { DrawnTrack } from '../services/track-legend-rows.js' */
/** @import { AddressedFeature } from '../services/feature-override-identity.js' */
/** @import { EditorPaintState } from './result-paint-record.js' */

const { watch } = window.Vue;

/**
 * What the editor intent reads of a mounted Result's binding context (the
 * preview runtime's `MountedResultContext`).
 * @typedef {object} DisplayedResultContext
 * @property {SVGSVGElement} root
 * @property {number} resultIndex
 * @property {string} resultIdentity
 * @property {number} rootGeneration
 * @property {string} phase
 * @property {Readonly<{ trustedRestore?: boolean, drawnPaint?: EditorPaintState }>} bindingOptions
 */

/**
 * @typedef {object} EditorIntentShowPorts
 * @property {WebState} state
 * @property {{ readonly value: Record<string, unknown>[] }} results The displayed drawing's Results.
 * @property {{ readonly value: readonly AddressedFeature[] | null }} extractedFeatures The features read from a Result without a catalog.
 * @property {{ readonly value: string[] }} originalLegendOrder The captions of the rows Python drew.
 * @property {{ readonly value: Element | null }} svgContainer The element that holds the mounted SVG.
 * @property {{ readonly value: number | string }} selectedResultIndex
 * @property {() => readonly DrawnTrack[]} readDrawnTracks The tracks the displayed Result draws.
 * @property {() => { diagramOptions?: unknown } | null} getCommittedCanonicalRenderRequest The committed render request.
 * @property {() => { resultIndex: number } | null} getActiveRuntime The mounted Result's runtime (preview runtime).
 * @property {(operations: SvgMutationOperations | null | undefined, options?: { domains?: readonly string[], afterApply?: ((svg: SVGSVGElement, reconciled: boolean) => void) | null }) => boolean} applyEditorOperations
 *   Reconciles the displayed Result with the operations and commits it once (preview runtime).
 * @property {() => unknown} readLegendRowNames Reads the Legend rows' names (Legend owner).
 * @property {(options?: { commit?: boolean }) => unknown} onLegendGeometryChanged Lays the Legend out (Legend owner).
 * @property {(label: string, fn: () => unknown) => Promise<unknown>} runUndoable One History step (History).
 */

/**
 * Registers nothing. The root creates it after the feature editor's actions,
 * and every earlier reader of what it returns is a closure that runs after
 * setup.
 * @param {EditorIntentShowPorts} ports
 */
export const createEditorIntentShow = ({
  state, results, extractedFeatures, originalLegendOrder, svgContainer, selectedResultIndex, readDrawnTracks,
  getCommittedCanonicalRenderRequest, getActiveRuntime, applyEditorOperations, readLegendRowNames, onLegendGeometryChanged,
  runUndoable
}) => {
  // What the live compile reads of the displayed Result's Legend state.
  /** @param {DisplayedLegendRowContext} context */
  const legendRowsOf = ({ pythonRows, features, sharedLegendKeys }) => ({ pythonRows, features, sharedLegendKeys });
  /**
   * @param {DrawingState} drawing
   * @param {number} resultIndex
   * @param {{
   *   replayDefaultLegendOrder?: string[] | null, domains?: readonly string[] | null,
   *   ruleRows?: import('./candidate-render.js').RuleLegendRows | null
   * }} [options]
   *   The operation domains the displayed Result shows, previewed live and compiled alone; null
   *   compiles the edits Python does not draw, as Generate does. `ruleRows`: a rule commit's rows.
   * @returns {SvgMutationOperations | null}
   */
  const compileDisplayedResultOperations = (drawing, resultIndex, { replayDefaultLegendOrder = null, domains = null, ruleRows = null } = {}) => {
    const catalog = toRaw(state.featureCatalog.value);
    const svg = svgContainer.value?.querySelector?.('svg') || null;
    // A Result without a catalog (a Session older than 40) is reached through
    // the features read from it and the fills it draws.
    const plan = compileDirectEditorMutationPlan({
      catalogAdmission: catalog
        ? admitFeatureCatalog(catalog, toRaw(results.value), { mode: drawing.mode })
        : displayedFeatureAddressing(
          toRaw(extractedFeatures.value) || [], /** @type {string[]} */ (results.value.map((/** @type {{ name?: unknown } | null | undefined} */ result) => result?.name)), resultIndex, svg
        ),
      featureColorOverrides: drawing.featureColorOverrides,
      featureStrokeOverrides: drawing.featureStrokeOverrides,
      featureOverrides: drawing.featureOverrides,
      legendEntries: drawing.legendEntries.value,
      deletedLegendEntries: drawing.deletedLegendEntries.value,
      dormantLegendEntries: drawing.dormantLegendEntries.value,
      originalLegendOrder: originalLegendOrder.value,
      addedLegendCaptions: drawing.addedLegendCaptions.value,
      legendColorOverrides: drawing.legendColorOverrides,
      legendStrokeOverrides: drawing.legendStrokeOverrides,
      manualSpecificRules: drawing.manualSpecificRules,
      replayDefaultLegendOrder,
      livePreview: domains ? {
        domains,
        paletteColors: appliedFeatureColors(state),
        drawnTracks: readDrawnTracks(),
        ...legendRowsOf(displayedLegendRowContext(state, drawing)),
        drawnContext: domains.includes('featureVisibility')
          // The committed request may be one a Session stored (unchecked): its options are read as unvalidated.
          ? featureDrawnContext(drawing, { diagramOptions: /** @type {UnvalidatedFields<'selectedFeaturesSet'>} */ (getCommittedCanonicalRenderRequest()?.diagramOptions) })
          : null,
        ruleRows
      } : null
    });
    return plan.operationsByResult[resultIndex] || null;
  };
  // The port of the editor intent onto the displayed Result (design §5.4):
  // the operations of `domains` Generate compiles from the active drawing
  // (the compile runs only the stages they need) unless given, reconciled by
  // the executor Generate uses, which returns each attribute of `domains` no
  // operation sets to Python's value; then one commit of the Result. Live
  // edits, History, and the display of a batch Result call it. The Legend is
  // laid out once when the reconcile changed its rows (zero shift).
  /**
   * @param {{
   *   domains: readonly string[],
   *   operations?: SvgMutationOperations | null,
   *   ruleRows?: import('./candidate-render.js').RuleLegendRows | null
   * }} options
   */
  const showEditorIntent = ({ domains, operations = undefined, ruleRows = null }) => {
    let shown = operations;
    if (shown === undefined) {
      const resultIndex = getActiveRuntime()?.resultIndex ?? (Number(selectedResultIndex.value) || 0);
      try {
        shown = compileDisplayedResultOperations(state.activeDrawing(), resultIndex, { domains, ruleRows });
      } catch (error) {
        console.error('Editor edits could not be compiled for the displayed Result.', normalizeUserFacingError(error));
        return false;
      }
    }
    return applyEditorOperations(shown, {
      domains,
      // The Legend is laid out without its own commit: the display commits
      // the Result once (S4).
      afterApply: (_svg, reconciled) => {
        if (!reconciled) return;
        readLegendRowNames();
        onLegendGeometryChanged({ commit: false });
      }
    });
  };
  // The paint domains each edit kind shows, live and on Undo and Redo.
  const {
    strokes: STROKE_DOMAINS, legendFills: LEGEND_FILL_DOMAINS, visibility: VISIBILITY_DOMAINS
  } = LIVE_EDIT_DOMAINS;
  /**
   * The one transition of a live edit of the paint intent (design §5.8): the
   * action writes the active drawing's intent and returns true when it
   * changed, the displayed Result shows it, and both are one History step (a
   * control's open step is joined).
   * @template {unknown[]} A
   * @param {string} label
   * @param {(...args: A) => unknown} mutate
   * @param {{ domains?: readonly string[] }} [options]
   * @returns {(...args: A) => Promise<unknown>}
   */
  const editEditorIntent = (label, mutate, { domains = STROKE_DOMAINS } = {}) => (
    ...args
  ) => runUndoable(label, () => showingEditorIntent(mutate, domains)(...args));
  /**
   * The action that writes the intent and, when it changed, shows `domains`.
   * @template {unknown[]} A
   * @param {(...args: A) => unknown} mutate
   * @param {readonly string[]} domains
   * @returns {(...args: A) => Promise<unknown>}
   */
  const showingEditorIntent = (mutate, domains) => async (...args) => {
    const changed = await mutate(...args);
    if (changed === true) showEditorIntent({ domains });
    return changed;
  };
  // A displayed Result shows the Legend structure edits and the paint domains
  // whose intent changed since it was last shown; each paint domain it shows
  // returns to Python's values where no operation sets it (OV-144).
  const LEGEND_STRUCTURE_DOMAINS = LIVE_EDIT_DOMAINS.legendStructure;
  return {
    compileDisplayedResultOperations, showEditorIntent, STROKE_DOMAINS, LEGEND_FILL_DOMAINS, VISIBILITY_DOMAINS,
    editEditorIntent, showingEditorIntent, LEGEND_STRUCTURE_DOMAINS
  };
};

/**
 * @typedef {object} MountedEditorIntentProjectionPorts
 * @property {WebState} state
 * @property {{ readonly value: Element | null }} svgContainer The element that holds the mounted SVG.
 * @property {(options?: { prepareRules?: boolean, show?: boolean, saving?: boolean }) => Promise<boolean> | boolean} projectPaletteAndRules
 *   Projects the palette and the rules; false when they did not project.
 * @property {(options: { rerender?: boolean, reflow?: boolean, labels?: boolean, show?: boolean }) => unknown} projectFeatureVisibility
 *   Projects the feature visibility (feature editor).
 * @property {() => unknown} reconcileLabelOverrides Reconciles the label intent (feature editor).
 * @property {(drawing: DrawingState) => readonly string[]} fillDomainsOf The domains a display of the feature fills shows.
 * @property {readonly string[]} VISIBILITY_DOMAINS
 * @property {(options: { domains: readonly string[] }) => unknown} showEditorIntent
 * @property {() => unknown} requestAutomaticRerender Asks Python to draw the Result again (feature editor).
 * @property {(label: string, fn: () => unknown) => Promise<unknown>} runUndoable One History step (History).
 */

/**
 * Registers nothing. The root creates it after `createEditorIntentShow` and
 * registers History's apply, which calls the projection, right after it.
 * @param {MountedEditorIntentProjectionPorts} ports
 */
export const createMountedEditorIntentProjection = ({
  state, svgContainer, projectPaletteAndRules, projectFeatureVisibility, reconcileLabelOverrides, fillDomainsOf,
  VISIBILITY_DOMAINS, showEditorIntent, requestAutomaticRerender, runUndoable
}) => {
  // One projection of the canonical editor intent onto the mounted Result,
  // shared by History apply, the display of another batch Result, and Load
  // Feature Edits TSV (D-07, R3). A History step of the Legend rows and a newly
  // displayed Result show the Legend operations Generate applies (`domains`).
  // A loaded table (`reflow`) also places the
  // labels, as a visibility edit does. The palette and the rules (`colors`)
  // and the visibility prepare their matches through their projections; the
  // fills, the visibility, and the other paint `domains` of the step then show
  // in one compile through `showEditorIntent`, or in the caller's (`show:
  // false`). Resolves to the domains shown, or false when the palette and
  // rules did not project (rules changed meanwhile, or a Session operation).
  /**
   * @param {{
   *   colors?: boolean, prepareRules?: boolean, visibility?: boolean, rerender?: boolean,
   *   reflow?: boolean, labels?: boolean, domains?: readonly string[], show?: boolean, saving?: boolean
   * }} [options]
   * @returns {Promise<string[] | false>}
   */
  const projectMountedEditorIntent = async ({
    colors = false,
    prepareRules = colors,
    visibility = false,
    rerender = false,
    reflow = false,
    labels = false,
    domains = [],
    show = true,
    saving = false
  } = {}) => {
    if (colors && !await projectPaletteAndRules({ prepareRules, show: false, saving })) return false;
    // The visibility projection follows with the labels, so one edit queues
    // one label request.
    if (visibility) await projectFeatureVisibility({ rerender, reflow, labels, show: false });
    if (labels && !visibility) reconcileLabelOverrides();
    const shown = [...new Set([...(colors ? fillDomainsOf(state.activeDrawing()) : []), ...(visibility ? VISIBILITY_DOMAINS : []), ...domains])];
    if (show && shown.length > 0) showEditorIntent({ domains: shown });
    return shown;
  };

  // O-2 (D-15-6 (5)): a Result saved before the executor kept a deleted row
  // lacks Python's row, so a Restore or Reset that returns the row asks Python
  // to draw it again. A Result Python drew in this session lacks only a row
  // Python did not draw there (its Legend row facts; U3a review M2); a loaded
  // Result has no facts.
  /** @param {string[]} keys The Python keys of the returning rows. */
  const drawReturningLegendRows = (keys) => {
    const svg = svgContainer.value?.querySelector?.('svg');
    const drawn = state.displayedResultMetadata()?.legendRowFacts?.drawn;
    if (keys.some((key) => (!drawn || drawn.has(key)) && !drawsPythonLegendRow(svg, key))) {
      requestAutomaticRerender();
    }
  };

  // OV-154: a Restore is one History step of the Legend intent, shown with the
  // rows' fills and strokes (OV-293) through the port.
  /** @param {() => unknown} restore Returns the Python keys of the restored rows, as `restoreDeletedLegendEntries` does. */
  const restoreLegendRows = (restore) => {
    const restored = restore();
    if (!Array.isArray(restored)) return restored;
    showEditorIntent({ domains: LIVE_EDIT_DOMAINS.deletedRows });
    drawReturningLegendRows(restored);
    return true;
  };
  /** @param {string} label @param {() => unknown} restore */
  const restoreLegendItems = (label, restore) => runUndoable(label, () => restoreLegendRows(restore));
  return { projectMountedEditorIntent, drawReturningLegendRows, restoreLegendRows, restoreLegendItems };
};

/**
 * @typedef {object} DepartedResultIntentPorts
 * @property {WebState} state
 * @property {{ readonly value: string[] }} originalLegendOrder The captions of the rows Python drew.
 * @property {(current: EditorPaintState) => string} resultPaintDepart
 *   The displayed Result leaves; its identity, or '' (Result paint record).
 */

/**
 * Registers nothing.
 * @param {DepartedResultIntentPorts} ports
 */
export const createDepartedResultIntent = ({ state, originalLegendOrder, resultPaintDepart }) => {
  /** @param {DrawingState} drawing @returns {EditorPaintState} */
  const currentEditorProjectionState = (drawing) => editorPaintState(state, drawing);
  // E1: the Result shown until a mode switch followed every live edit. When
  // it is shown again with the same editor intent, nothing is projected and
  // its bytes stay; otherwise the edits made since are projected onto it.
  /** @type {Map<string, string>} */
  const departedResultIntent = new Map();
  /** @param {DrawingState} drawing */
  const displayedIntentSignature = (drawing) => JSON.stringify([
    currentEditorProjectionState(drawing),
    drawing.deletedLegendEntries.value.map((/** @type {KeptLegendRow} */ entry) => entry.originalCaption || entry.caption),
    [...(drawing.addedLegendCaptions.value || [])],
    drawing.legendEntries.value.filter((/** @type {LegendRow} */ entry) => /** @type {string} */ (entry.caption) !== entry.originalCaption)
      .map((/** @type {LegendRow} */ entry) => [entry.originalCaption, entry.caption]),
    // The direct additions, which no Result's inventory lists.
    drawing.legendEntries.value.filter((/** @type {LegendRow} */ entry) => !originalLegendOrder.value.includes(entry.originalCaption || entry.caption))
      .map((/** @type {LegendRow} */ entry) => entry.caption)
  ]);
  /** @param {DrawingState} drawing */
  const rememberDepartingResultProjection = (drawing) => {
    const departed = resultPaintDepart(currentEditorProjectionState(drawing));
    if (departed) departedResultIntent.set(departed, displayedIntentSignature(drawing));
  };
  // Whether a Result draws a feature; an unknown catalog answers yes.
  /** @param {DrawingState} drawing @param {number} resultIndex @returns {(recordKey: string, featureId: string) => boolean} */
  const resultDrawsFeature = (drawing, resultIndex) => {
    const rendered = resultCatalogFeatures(state, drawing.mode, resultIndex)?.renderedByIdentity;
    return (recordKey, featureId) => !rendered || rendered.has(biologicalFeatureKey(recordKey, featureId));
  };
  return {
    currentEditorProjectionState, departedResultIntent, displayedIntentSignature, rememberDepartingResultProjection,
    resultDrawsFeature
  };
};

/**
 * @typedef {object} DeclinedPaintDisplayPorts
 * @property {WebState} state
 * @property {(saving: boolean) => boolean} sessionBlocks Whether a Session operation blocks the projection.
 * @property {{ readonly value: boolean }} ruleMatchingPending Whether the rule preparation runs.
 * @property {(drawing: DrawingState) => EditorPaintState} currentEditorProjectionState
 * @property {(options?: { colors?: boolean, visibility?: boolean, saving?: boolean }) => Promise<string[] | false>} projectMountedEditorIntent
 * @property {(identity: string) => Partial<EditorPaintState> | null} resultPaintLacking
 *   What the displayed Result still shows of the state before a declined display (Result paint record).
 * @property {(identity: string, current: EditorPaintState, previous: EditorPaintState) => void} resultPaintShown
 *   Records what the displayed Result shows (Result paint record).
 */

/**
 * Registers nothing at setup; a declined display watches its blocker at run
 * time.
 * @param {DeclinedPaintDisplayPorts} ports
 */
export const createDeclinedPaintDisplay = ({
  state, sessionBlocks, ruleMatchingPending, currentEditorProjectionState, projectMountedEditorIntent,
  resultPaintLacking, resultPaintShown
}) => {
  // A display that a Session operation or the rule preparation declined shows
  // the fills and visibility once the blocker clears, while its Result stays
  // displayed, so the next display reads the Result with them. A Save shows
  // them before it writes the Result (`saving`: the Save does not block). A
  // projection declined again without a blocker waits for the next display.
  const displayBlocked = (saving = false) => sessionBlocks(saving) || ruleMatchingPending.value;
  /** @type {(() => void) | null} */
  let stopDeclinedPaintWait = null;
  /**
   * @param {string} identity
   * @param {{ saving?: boolean }} [options]
   */
  const showDeclinedPaint = async (identity, { saving = false } = {}) => {
    stopDeclinedPaintWait?.();
    stopDeclinedPaintWait = null;
    if (displayBlocked(saving)) {
      stopDeclinedPaintWait = watch(() => displayBlocked(), (/** @type {unknown} */ blocked) => { if (!blocked) void showDeclinedPaint(identity); });
      return;
    }
    const lacking = resultPaintLacking(identity);
    if (!lacking?.colors || lacking.visibility === undefined) return;
    const current = currentEditorProjectionState(state.activeDrawing());
    try {
      const shown = await projectMountedEditorIntent({
        colors: !sameColors(lacking.colors, current.colors),
        visibility: lacking.visibility !== current.visibility,
        saving
      });
      if (shown === false) {
        if (displayBlocked()) void showDeclinedPaint(identity);
        return;
      }
      if (resultPaintLacking(identity) === lacking) resultPaintShown(identity, current, current);
    } catch (error) {
      console.error('Editor edits could not be shown on the displayed Result.', normalizeUserFacingError(error));
    }
  };

  // A label rerender's displayed Result, once bound, shows the fills and
  // visibility of the edits made while it ran (OV-346).
  /** @param {Pick<DisplayedResultContext, 'resultIdentity' | 'bindingOptions'>} context */
  const showEditsMadeWhileDrawn = (context) => {
    if (context.bindingOptions.drawnPaint && resultPaintLacking(context.resultIdentity)) {
      void showDeclinedPaint(context.resultIdentity);
    }
  };
  return { displayBlocked, showDeclinedPaint, showEditsMadeWhileDrawn };
};

/**
 * @typedef {object} DisplayedResultProjectionPorts
 * @property {WebState} state
 * @property {() => string[]} liveResultIdentities The Results of all drawings.
 * @property {(drawingId: DrawingId) => string[]} resultIdentitiesOf The Results of one drawing.
 * @property {(drawing: DrawingState) => EditorPaintState} currentEditorProjectionState
 * @property {Map<string, string>} departedResultIntent The intent each departed Result left with.
 * @property {(drawing: DrawingState) => string} displayedIntentSignature
 * @property {(drawing: DrawingState, resultIndex: number) => (recordKey: string, featureId: string) => boolean} resultDrawsFeature
 * @property {(drawing: DrawingState, resultIndex: number, options?: { replayDefaultLegendOrder?: string[] | null, domains?: readonly string[] | null }) => SvgMutationOperations | null} compileDisplayedResultOperations
 * @property {(options: { colors?: boolean, visibility?: boolean, rerender?: boolean, show?: boolean }) => Promise<string[] | false>} projectMountedEditorIntent
 * @property {(options: { domains: readonly string[], operations?: SvgMutationOperations | null }) => unknown} showEditorIntent
 * @property {(saving?: boolean) => boolean} displayBlocked
 * @property {(identity: string) => Promise<void>} showDeclinedPaint
 * @property {readonly string[]} LEGEND_STRUCTURE_DOMAINS
 * @property {readonly string[]} STROKE_DOMAINS
 * @property {(identities: readonly string[], displayed: string, current: EditorPaintState, options?: { restored?: boolean, drawn?: EditorPaintState, drawingStateOf?: (identity: string) => EditorPaintState }) => void} resultPaintCommit
 * @property {(identity: string) => void} resultPaintRebind
 * @property {(identity: string, current: EditorPaintState) => EditorPaintState} resultPaintDisplay
 * @property {(identity: string) => EditorPaintState | undefined} resultPaintDrawn
 * @property {(identity: string, current: EditorPaintState, previous: EditorPaintState, options?: { declined?: boolean, failed?: boolean }) => void} resultPaintShown
 *   The Result paint record's methods.
 */

/**
 * Registers nothing. The root creates it last of the editor intent, before
 * the Result binder, whose steps call what it returns.
 * @param {DisplayedResultProjectionPorts} ports
 */
export const createDisplayedResultProjection = ({
  state, liveResultIdentities, resultIdentitiesOf, currentEditorProjectionState, departedResultIntent,
  displayedIntentSignature, resultDrawsFeature, compileDisplayedResultOperations, projectMountedEditorIntent,
  showEditorIntent, displayBlocked, showDeclinedPaint, LEGEND_STRUCTURE_DOMAINS, STROKE_DOMAINS, resultPaintCommit,
  resultPaintRebind, resultPaintDisplay, resultPaintDrawn, resultPaintShown
}) => {
  // A displayed Result whose label intent changed since it was last shown
  // receives the label projection in the binder's label step, also when no
  // label intent remains (an undone or replaced Label TSV import).
  let labelProjectionResultIdentity = '';
  /** @param {DrawingState} drawing @param {DisplayedResultContext} context */
  const rememberCommittedEditorState = (drawing, context) => {
    const identities = liveResultIdentities();
    const current = currentEditorProjectionState(drawing);
    // A label rerender's Results show the editor state it drew from; the
    // edits made while it ran are shown on the displayed Result once it is
    // bound (`afterReady`, OV-346), and on another Result when displayed.
    const drawn = context.bindingOptions.drawnPaint || current;
    if (drawn.labels !== current.labels) labelProjectionResultIdentity = context.resultIdentity;
    // A loaded Session and a History restore keep each Result's bytes as
    // they were saved or kept, a batch Result not displayed since the last
    // edits without them (also a Session older than 40, whose strokes
    // reached only the mounted Result).
    resultPaintCommit(identities, context.resultIdentity, current, {
      restored: context.phase === 'session-load' || Boolean(context.bindingOptions.trustedRestore),
      drawn,
      drawingStateOf: (identity) => {
        const holder = state.drawingList().find(({ id }) => resultIdentitiesOf(id).includes(identity));
        return holder ? currentEditorProjectionState(holder) : current;
      }
    });
    [...departedResultIntent.keys()].forEach((identity) => {
      if (!identities.includes(identity)) departedResultIntent.delete(identity);
    });
  };
  // What a displayed Result asks of the rerender, read once it is ready
  // (`afterReady`): `rules`, the rule table it was last shown with, or the one
  // it was drawn with when that is unknown (a restored Result), for OV-1280;
  // `labelIntentGone`, whether its bytes hold a label edit that left its drawing
  // since it was last shown (OV-2293).
  /** @type {{ identity: string, rules: SpecificColorRule[] | null, labelIntentGone: boolean } | null} */
  let displayAsks = null;
  // D-07 (PD-OI-062): a batch Result shows the canonical color, visibility,
  // Legend, and label edits when it is displayed. Labels follow in the
  // binder's label step.
  /**
   * @param {DrawingState} drawing
   * @param {DisplayedResultContext} context
   * @param {string[]} resultLegendOrder
   */
  const projectEditorIntentOnDisplay = async (drawing, context, resultLegendOrder) => {
    const identity = context.resultIdentity;
    const current = currentEditorProjectionState(drawing);
    displayAsks = null;
    const departedIntent = departedResultIntent.get(identity);
    departedResultIntent.delete(identity);
    if (departedIntent !== undefined && departedIntent === displayedIntentSignature(drawing)) {
      resultPaintRebind(identity);
      labelProjectionResultIdentity = '';
      recordStructuralMetric('displayedResultEditorProjectionCount', 0, {
        phase: context.phase,
        rootGeneration: context.rootGeneration
      });
      return;
    }
    const previous = resultPaintDisplay(identity, current);
    const drawn = resultPaintDrawn(identity) ?? current;
    displayAsks = {
      identity,
      rules: paintRuleTable(previous) ?? paintRuleTable(drawn),
      labelIntentGone: drawn.labels !== current.labels && drawnLabelIntentGone(
        drawn.labels, previous.labels, current.labels, resultDrawsFeature(drawing, context.resultIndex)
      )
    };
    const colors = !sameColors(previous.colors, current.colors);
    const visibility = previous.visibility !== current.visibility;
    const strokes = previous.strokes !== current.strokes;
    labelProjectionResultIdentity = previous.labels !== current.labels ? identity : '';
    // B20: a Result last shown with another Legend order receives the current
    // order, also the default order, which is its own generated order (OV-47).
    const replayDefaultLegendOrder = previous.legendOrder !== current.legendOrder
      ? resultLegendOrder : null;
    /** @type {SvgMutationOperations | null} */
    let operations = null;
    /** @type {string[]} */
    let domains = [...LEGEND_STRUCTURE_DOMAINS];
    // The palette and rules declined (rules changed meanwhile, or a Session
    // operation): the Result keeps the fills and visibility it shows until
    // the blocker clears. A failed projection or compile showed none.
    let declined = false;
    let failed = false;
    try {
      // The palette, rules, and visibility prepare their matches when their
      // intent changed (`colors`, `visibility`); then one compile of the
      // domains the Result shows.
      // OV-2293: a Result whose bytes were drawn with other visibility, or
      // restored and not known to show any (a Reset before a Save and Load), asks
      // for the rerender when they lack a feature.
      const preparation = {
        colors, visibility, rerender: previous.visibility === '' || drawn.visibility !== current.visibility, show: false
      };
      let prepared = await projectMountedEditorIntent(preparation);
      // Declined with nothing blocking (the rules changed meanwhile): prepared
      // once more here, so the display still compiles once (R14-3).
      if (prepared === false && !displayBlocked()) prepared = await projectMountedEditorIntent(preparation);
      declined = prepared === false;
      domains = [...domains, ...(prepared || []), ...(strokes ? STROKE_DOMAINS : [])];
      operations = compileDisplayedResultOperations(drawing, context.resultIndex, { replayDefaultLegendOrder, domains });
    } catch (error) {
      failed = true;
      console.error('Editor edits could not be compiled for the displayed Result.', normalizeUserFacingError(error));
    }
    const shown = operations;
    const hasOperations = shown !== null
      && domains.some((domain) => shown[/** @type {keyof SvgMutationOperations} */ (domain)].length > 0);
    // A Result that shows a Legend structure edit no operation keeps returns
    // to Python's rows (the reconcile's reverse).
    const projects = colors || visibility || strokes || hasOperations || legendStructureEdited(context.root);
    recordStructuralMetric('displayedResultEditorProjectionCount', projects ? 1 : 0, {
      phase: context.phase,
      rootGeneration: context.rootGeneration
    });
    if (projects) {
      try {
        showEditorIntent({ domains, operations });
      } catch (error) {
        failed = true;
        console.error('Editor edits could not be shown on the displayed Result.', normalizeUserFacingError(error));
      }
    }
    resultPaintShown(identity, current, previous, { declined, failed });
    if (declined && !failed && displayBlocked()) await showDeclinedPaint(identity);
  };
  // The binder's label step and its last step each read once what the
  // projection kept for them (`synchronizeLabelEditor`, `afterReady`).
  const takeLabelProjectionResultIdentity = () => {
    const identity = labelProjectionResultIdentity;
    labelProjectionResultIdentity = '';
    return identity;
  };
  const takeDisplayAsks = () => {
    const asks = displayAsks;
    displayAsks = null;
    return asks;
  };
  return { rememberCommittedEditorState, projectEditorIntentOnDisplay, takeLabelProjectionResultIdentity, takeDisplayAsks };
};
