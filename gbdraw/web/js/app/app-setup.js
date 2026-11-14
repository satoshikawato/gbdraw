// @ts-check
/** @import { FeaturePopup } from '../services/feature-popup.js' */
/** @import { RecordPlacementRun } from './align-at-feature.js' */
/** @import { GenerationPayload } from '../services/diagram-generation.js' */
/** @import { RecordSpanRead } from '../services/diagram-worker-protocol.js' */
/** @import { LinearUnitPlanPayload } from '../services/linear-units.js' */
/** @import { PlanGenomeScales } from './run-analysis.js' */
/** @import { CardWindowRecordInput } from '../contracts/card-windows.generated.js' */
/** @import { ResultLabelFit } from '../services/label-fit.js' */
/** @import { DiagramMode, DrawingState, LinearSeq } from '../state.js' */
/** @import { NoticeLink } from './app-messages.js' */
/** @import { RulePreparation } from './rule-matching.js' */
/** @import { CircularDiagramRecords } from '../services/record-options.js' */
/** @import { FeatureEditorOptions } from './feature-editor.js' */
/** @import { CommitSpecificRulesPort } from './feature-editor/color-actions.js' */
/** @import { LegendRowColorOf } from './candidate-render.js' */
/** @import { LegacyResultSvgData, LegacyResultSvgTransform, SessionFileInput } from '../services/config.js' */
/** @import { RecordKeyedFeature } from '../services/feature-identity.js' */
/** @import { ArtifactSlot, ArtifactSlotSnapshot, ArtifactSlotStore, LoadedArtifactSlots } from '../services/artifact-slot.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { GalleryExample } from '../services/gallery-examples.js' */
/** @import { LinearSourceGroup } from '../services/linear-sources.js' */
/** @import { RecordIdentitySources } from '../services/record-draw-selection.js' */
/** @import { RegionSourceFeature } from './region-drawing.js' */
import { createRulePreparation } from './rule-matching.js';
import { displayedLegendRowColors, editorPaintDomains, LIVE_EDIT_DOMAINS, RECORDED_LEGEND_ROW_COLORS } from './candidate-render.js';
import { createDeclinedPaintDisplay, createDepartedResultIntent, createDisplayedResultProjection, createEditorIntentShow, createMountedEditorIntentProjection } from './editor-intent.js';
import { createMountedLegendBinding, createMountedResultReadiness } from './mounted-result-binder.js';
import { createResultPaintRecord, visibilityPaintState } from './result-paint-record.js';
import {
  displayedLegendRowContext, requestFeatureVisibilityRules
} from '../services/feature-visibility.js';
import { legendEntryKey } from '../services/drawing/paint.js';
import { DEFAULT_LINEAR_OUTPUT_PREFIX } from '../utils/download-names.js';
import { ONE_CIRCLE_LIMITS } from '../services/one-circle-layout.js';
import { catalogNameIndex } from '../services/feature-catalog.js';
import { labelCapacityNoticeText, labelCapacityText } from '../services/label-fit.js';
import { prefetchPyodideBootAssets } from '../services/pyodide-assets.js';
import { PNG_DPI_CHOICES } from '../services/session-active-config-contract.js';
import { createSettingsPaneController } from './settings-pane.js';
import {
  createLinearComparisonFocus, createLinearComparisonTimeline, createLinearComparisonView, createLinearLayoutReaders,
  linearRecordRowFor
} from './linear-comparison-panel.js';
import {
  createLinearComparisonGapEdits, createLinearComparisonModes, createLinearComparisonPairEdits,
  createLinearComparisonPlanWriter
} from './linear-comparison-edits.js';
import { createLinearFileCards } from './linear-file-cards.js';
import { createCircularFileCards } from './circular-file-cards.js';
import { createInputTypeEdits, createInputUploads } from './input-edits.js';
import {
  createLinearRecordEdits,
  createLinearSourceEdits,
  createLinearSourceFileReplacement,
  createLinearSourceMoves,
  createLinearSourceRemovalDialog
} from './linear-source-edits.js';
import { createRecordDisplayControls, RECORD_TARGET_NOT_DISCOVERED } from './record-display-options.js';
import { createPreviewItemPopup, defaultColorScopeLine } from './preview-item-popup.js';
import { previewItemPopupOf } from '../services/match-frame-popup.js';
import { legendRowRules, specificRuleHasQualifier } from '../services/specific-color-rules.js';
import { drawnTracks, hiddenTrackLegendKeys, loadedTracks, trackLegendRows } from '../services/track-legend-rows.js';
import { isCurrentFeature } from '../services/feature-identity.js';
import {
  createFeatureRecordRotationAction,
  createFeatureRecordRotationWorkflow
} from './record-display/feature-record-rotation.js';
import { createDrawingState, drawingHasEdits, state, sessionOperationAvailability } from '../state.js';
import { createDrawingBar, createDrawingList, createInputReplacement, moveTabFocus } from './drawing-list.js';
import { createDrawingSwitch } from './drawing-switch.js';
import { createInputReplacementPorts } from './input-replacement-ports.js';
import { createDrawingReset, createResetSettings } from './reset-settings.js';
import { createRegionDrawing } from './region-drawing.js';
import { createFindInOtherFiles } from './find-in-other-files.js';
import { createRegionSelect } from './region-select.js';
import { readRegionAxisRecords, shownAxisRecords as shownAxisRecordsOf } from '../services/region-axis.js';
import { requestRecordSpanSequences } from '../services/record-span-sequences.js';
import { catalogSequenceHelpers, createCatalogSequences } from '../services/catalog-sequences.js';
import { createSequenceReads } from './sequence-reads.js';
import { downloadTextFile } from '../services/text-download.js';
import {
  adoptCanonicalRenderArtifacts,
  applyCapturedConfigData,
  applyEditorStateData,
  applyFeatureStateData,
  applyOrthogroupStateData,
  applyResultsData,
  applyRunStateData,
  applyUiStateData,
  applySessionUiData,
  buildConfigData,
  buildEditorStateData,
  buildFeatureStateData,
  buildOrthogroupStateData,
  buildRunStateData,
  buildUiStateData,
  canonicalRenderArtifactOwner,
  exportSession,
  disposeSessionOperations,
  getCommittedCanonicalSession,
  getCommittedCanonicalRenderRequest,
  committedLinearRecordCards,
  committedSettingNotePredicate,
  committedOneCircleChange,
  committedTrackHideDrawing,
  commitLiveSettingValues,
  readCommittedCanonicalRevision,
  getLoadedTrackSettings,
  getCommittedLinearDefinitionVisibility,
  readCommittedResourceRecordCount,
  assertDrawingInputs,
  importSession as importSessionFromFile,
  inputFilePayload,
  serializeActiveRenderFiles,
  serializeFile,
  writeSessionDocumentFile,
  projectDraftRequest,
  makeRegionDrawing,
  heldDepthFiles,
  setUnmanagedConfigOverrideValidator,
  setPreviewRuntime,
  applyLabelFitValues,
  setHeldPreviewViewportReader,
  installShownDrawing,
  duplicateDrawingState,
  applyDrawingStyle,
  writeLosatCacheInfo,
  CIRCULAR_RECORD_FORM_DEFAULTS,
  linearCardWindowPlanRecords,
  linearDrawingUnitPlaces,
  linearDrawingWindowsOf
} from '../services/config.js';
import { artifactSlotAs } from '../services/artifact-slot.js';
import {
  setStoredProteinResultFrameConverter
} from '../services/main-session-comparison-frame.js';
import { setLegacySessionUpgrader } from '../services/legacy-session-upgrade.js';
import { fetchGalleryExampleSession, readGalleryExamples } from '../services/gallery-examples.js';
import { createHistoryRouter } from '../services/history.js';
import { bindDrawingInputFileTree, createLinearSeq, retainHistoryAndBoundFiles, boundInputFileIds } from '../services/drawing/inputs.js';
import { createHistorySnapshotService } from '../services/history-snapshot.js';
import {
  linearFileCardRecord,
  getLinearSourceDefaultDefinition,
  inferredDefinitionForRecord,
  setLinearSourceDefaultDefinition,
  getLinearSourceDefaultSubtitle,
  setLinearSourceDefaultSubtitle,
  resolveLinearRecordEffectiveDefinition,
  resolveLinearRecordEffectiveSubtitle,
  losatGencodeProblem
} from '../services/linear-sources.js';
import { LOSAT_GENETIC_CODE_TEXT } from '../contracts/schema.generated.js';
import { captureSvgExport, serializeCleanSvg } from '../services/svg-serialization.js';
import { copyPendingTextToClipboard, copyTextToClipboard } from '../utils/clipboard.js';
import { toRaw } from '../utils/vue-raw.js';
import { readHostedBuild } from '../utils/hosted-build.js';
import { withNewDrawingDefaults } from '../services/reset.js';
import { allLiveSettingRequestPaths, liveSettingRequestPaths, SETTING_NOTE_APPLIES, SETTING_NOTE_ONE_CIRCLE_CHANGE, changedTrackHides } from '../services/setting-application.js';
import {
  DIAGRAM_HELPER_OPERATIONS,
  disposeDiagramGenerationWorker,
  runDiagramHelperOperation,
  setDiagramResourceRetainSource,
  setDiagramWorkerFontSource
} from '../services/diagram-generation.js';
import { createFontSettings } from './font-settings.js';
import { customFontSaveWarnings } from '../services/session-fonts.js';
import { requestFontFamily, webFontFamilyNames } from '../services/font-store.js';
import {
  recordSessionLifecycleEvent
} from '../services/runtime-test-hooks.js';
import { createPanZoom, createSidebarResize, setupGlobalUiEvents } from './ui.js';
import { createFeatureEditor } from './feature-editor.js';
import { COLLINEAR_COLOR_MODE_TEXT } from '../services/pairwise-match-popup.js';
import { createFeatureSelection } from './feature-selection.js';
import { createPreviewFeatureSearch } from './feature-search/preview-actions.js';
import { createSvgStyles } from './svg-styles.js';
import { createLegendManager } from './legend.js';
import { createPaletteLoader } from './palettes.js';
import { afterPaint, createRunAnalysis } from './run-analysis.js';
import { createSimilarityAlignmentActions, featureNameParts } from './similarity-alignment.js';
import { createRecordAutoOrientation } from './record-auto-orientation.js';
import { createFeatureAnchorAction } from './align-at-feature.js';
import { createRecordPositionReset } from './record-position-reset.js';
import { createLabelFitDialog } from './label-fit-dialog.js';
import { createResultPopups } from './result-popups.js';
import {
  activeDepthTrackCount, circularDepthRepresentatives, circularDepthRows, circularDepthSources, createDepthFiles,
  createDepthSeriesEdits, createDepthTrackPanels, depthTrackCountLabel, ensureDepthTrackConfigCount, getLinearDepthFile,
  hasLinearDepthFiles, linearDepthLogicalWidth, linearDepthRows, linearDepthSources, normalizeDepthSlotTrackIndex,
  refreshDepthTrackLabelsAfterRemoval, syncDepthTrackSlotLabelsForTrack
} from './depth-tracks.js';
import { createLinearRecordLabels, definitionLineStyleRows } from './linear-record-labels.js';
import {
  createCircularConservationComparisonFiles, createCircularConservationSeriesPanel, createCircularConservationSources
} from './circular-conservation.js';
import { createCircularRecordOrder, createCircularRecordPresentation } from './circular-record-panel.js';
import {
  createClickedFeatureCopy, createImportedComparisonChoice, createMissingBlastTsvChoice, createResultExport, createRunInfoCopy
} from './result-actions.js';
import {
  diagnosticError, losatGeneticCodeHint, normalizeCaughtError, normalizeUserFacingError
} from '../utils/error-normalization.js';
import { shownRunInfoNoteParts } from './run-info.js';
import { createLegendLayout } from './legend-layout.js';
import { pythonLegendRowKeys, stepRefitsDrawnLegend } from '../services/legend-draft-rows.js';
/** @import { HistoryChangeLike } from '../services/setting-application.js' */
import {
  COMPOSITION_METADATA_ATTRIBUTE,
  COMPOSITION_SCHEMA_ATTRIBUTE,
  normalizeLegacyComposition,
  normalizeStoredLegendSide,
  requestedCanvasOf,
  storedLegendSideIsOfMode
} from './legend-layout/composition-actions.js';
import { createResultsManager } from './results.js';
import { createAuxiliaryFileImports, setupWatchers } from './watchers.js';
import { createDialogChoice, setupHistoryInputs } from './history-inputs.js';
import { createHistoryControls } from './history-list.js';
import { createChoicePrompt, createNoticeBanner } from './app-messages.js';
import { setupHistoryShortcuts } from './history-shortcuts.js';
import { createPreviewRuntime } from './preview-runtime.js';
import {
  createOrthogroupEditor
} from './orthogroups.js';
import { createRightDrawerController } from './right-drawer.js';
import { createCircularTrackSlotEditor } from './circular-track-slots.js';
import { createLinearTrackSlotEditor } from './linear-track-slots.js';
import { askTrackSourceRemoval, askingAnnotationTrackAdd, createAskingStackPorts } from './track-slot-edits.js';
import { createFeatureDensityTrackEditor } from './feature-density-tracks.js';
import { createAnnotationEditor } from './annotations.js';
import { createDrawingLegendRetirement, trackDataLegendCaptions } from './legend/track-data-styles.js';
import {
  reconcileAnnotationRecordBindings, resultRecordBindings
} from './annotations/record-selector.js';
import { validateAnnotationRecordTargets } from './annotations/validation.js';
import { createLosatSettings } from './losat-settings.js';
import { createAutoValueDisplay } from './auto-value-display.js';
import { trapDialogFocus } from './dialog-focus-trap.js';
import { createFeatureEditHistory, createLegendEditHistory, createSelectedFeatureHistory } from './history-wrappers.js';
import { createPreservedConfigOverrides } from './preserved-config-overrides.js';
import { createSessionStatus } from './session-status.js';
import { createFeatureChoiceReaders, createScaleStyleAvailability, createShownDrawingMembers } from './shown-drawing-settings.js';
import { createCircularMultiRecordCanvas, createCircularRecordDisplay, createCircularSourceEdits } from './circular-record-controls.js';
import { createAnnotationRecordRebind, createDrawingRecordCatalog } from './drawing-record-catalog.js';
import { createRecordDiscoveryRequests, createRecordDiscoveryRuns } from './record-discovery-requests.js';
import { createRecordSettingsSources, createRecordSettingsTransitions } from './record-settings.js';
import { createShownDrawingReaders } from './shown-drawing-readers.js';
import { createLinearWrapField } from './linear-wrap-field.js';
import { createFigureSizeField } from './figure-size-field.js';
import { createFigureSizeReaders } from './figure-size-readers.js';
import { createOrthogroupPopup } from './orthogroup-popup.js';
import { createSimilarityAlignmentReview, createSimilarityAlignmentTakeover } from './similarity-alignment-review.js';
import { createLinearRecordSelector } from './linear-record-selector.js';
import { createLinearTypographyController } from './linear-typography.js';
import {
  buildDisambiguatedRecordEntries,
  circularEstimateLengths,
  circularRequestDiagrams,
  formatRecordLength
} from '../services/record-options.js';
import {
  setLinearRecordRow as updateLinearRecordRow
} from '../services/linear-record-layout.js';
import {
  missingUploadPairsToResolve
} from '../services/linear-comparisons.js';
import { catalogUnitPlaces } from '../services/linear-units.js';
import { addAnnotationSet as addSearchSetToDrawing } from '../services/drawing/annotation-sets.js';
import { createSequenceSearch } from './sequence-search.js';
import { resultRecordBinding } from '../services/sequence-search.js';
import { createRecordSelection } from './record-selection.js';
import {
  comparisonProfileDefault, IDENTITY_PERCENT_MAXIMUM, identityScaleWithEnd
} from '../mode-profiles.js';
import { joinSentences, msg, t, tm, tx } from '../utils/i18n.js';
import { errorText } from '../utils/error-text.js';
import { importedComparisonExecution } from '../services/imported-comparison-intent.js';
import {
  activeDepthTrackIndices,
  compactDepthFileSlots,
  depthSlotTrackIndex,
  isDefaultManagedDepthSlot,
  representativeDepthFiles,
  reindexDepthSlots,
  removeDepthTrackColumnAt,
  depthSeriesRemovalLoss
} from '../services/depth-track-state.js';
/** @import { HistoryManagerOptions } from '../services/history.js' */
/** @import { RulePreparationOptions } from './rule-matching.js' */
/** @import { SpecificColorRule } from '../services/specific-color-rules.js' */
/** @import { ReindexDepthSlotsOptions } from '../services/depth-track-state.js' */
/** @import { CircularTrackSlot } from '../services/circular-track-slot-model.js' */
/** @import { LinearTrackSlot } from '../services/linear-track-slot-model.js' */
/** @import { SequenceRecord } from '../services/record-discovery.js' */
/** @import { AlignmentOrientationCheckpoint, TargetDraftCheckpoint } from './record-display-options.js' */

const { onMounted, onUnmounted, watch, nextTick, computed, provide, ref, shallowRef, reactive } = window.Vue;

/** @type {Promise<typeof import('../services/export.js')> | null} */
let exportServicePromise = null;

const loadExportService = () => {
  exportServicePromise ??= import('../services/export.js');
  return exportServicePromise;
};

/** @typedef {ReturnType<typeof createHistoryRouter>} HistoryRouter */
/**
 * @typedef {object} SessionImportRollbackOptions
 * @property {ArtifactSlotStore} artifactSlots
 *   The kept generated artifacts of the drawings not shown (E1), held by reference.
 * @property {HistoryRouter['captureMembers']} captureHistoryMembers
 *   The History router's members: a Load replaces them at its store commit (OV-1769, OV-2131).
 * @property {HistoryRouter['replaceMembers']} replaceHistoryMembers
 * @property {HistoryRouter['restoreMembers']} restoreHistoryMembers
 * @property {() => Readonly<ArtifactSlot>} captureDisplayedArtifact The History snapshot service's slot capture
 *   of the shown drawing.
 * @property {(slot: Readonly<ArtifactSlot>) => void} installDisplayedArtifact
 *   Installs a captured slot after the import rollback restored the document.
 * @property {{ value: number }} featureListScrollTop
 * @property {{ value: HTMLElement | null }} featureListScrollRef
 * @property {{ value: string }} selectedPairwiseBlockOrthogroupId
 * @property {(() => any) | null} [captureSpecificRulePatternDrafts] The feature editor's capture of the unsaved rule pattern drafts.
 * @property {((drafts: any) => void) | null} [restoreSpecificRulePatternDrafts]
 */

/**
 * What `capture` holds for `restore`.
 * @typedef {object} SessionImportRollbackSnapshot
 * @property {ArtifactSlotSnapshot} artifactSlots
 * @property {ReturnType<HistoryRouter['captureMembers']>} historyMembers
 * @property {Readonly<ArtifactSlot>} displayedArtifact
 * @property {number} featureListScrollTop
 * @property {string} selectedPairwiseBlockOrthogroupId
 * @property {unknown} [specificRulePatternDrafts] The unsaved rule pattern drafts, as the feature editor captured them; handed back as they are.
 */

/** @param {SessionImportRollbackOptions} options */
export const createSessionImportRollbackState = ({
  artifactSlots,
  captureHistoryMembers,
  replaceHistoryMembers,
  restoreHistoryMembers,
  captureDisplayedArtifact,
  installDisplayedArtifact,
  featureListScrollTop,
  featureListScrollRef,
  selectedPairwiseBlockOrthogroupId,
  captureSpecificRulePatternDrafts = null,
  restoreSpecificRulePatternDrafts = null
}) => ({
  capture: () => ({
    artifactSlots: artifactSlots.snapshot(),
    historyMembers: captureHistoryMembers(),
    // The displayed artifact's parts outside state: its transport identity,
    // match-sequence owner, committed Session and CLI helper files.
    displayedArtifact: captureDisplayedArtifact(),
    featureListScrollTop: featureListScrollTop.value,
    selectedPairwiseBlockOrthogroupId: selectedPairwiseBlockOrthogroupId.value,
    ...(captureSpecificRulePatternDrafts ? { specificRulePatternDrafts: captureSpecificRulePatternDrafts() } : {})
  }),
  // A Load's commit: each drawing it made starts with an empty History.
  replaceHistoryMembers,
  // The previous drawings' Histories, which the Load replaced (`replaceMembers`): the rollback gives them back in the
  // task that restores the store, before `restore` (OV-2132).
  restoreHistoryMembers: (/** @type {SessionImportRollbackSnapshot} */ snapshot) => {
    restoreHistoryMembers(snapshot.historyMembers);
  },
  restore: async (/** @type {SessionImportRollbackSnapshot} */ snapshot) => {
    artifactSlots.restore(snapshot.artifactSlots);
    installDisplayedArtifact(snapshot.displayedArtifact);
    // The restored document holds each drawing's Depth series, and so its Depth panel count (OV-1165, OV-2200).
    await nextTick();
    featureListScrollTop.value = snapshot.featureListScrollTop;
    if (featureListScrollRef.value) {
      featureListScrollRef.value.scrollTop = snapshot.featureListScrollTop;
    }
    selectedPairwiseBlockOrthogroupId.value =
      snapshot.selectedPairwiseBlockOrthogroupId;
    if (Object.hasOwn(snapshot, 'specificRulePatternDrafts')) {
      restoreSpecificRulePatternDrafts?.(snapshot.specificRulePatternDrafts);
    }
  }
});

// A Result's root tag names the composition metadata attributes.
const COMPOSITION_ROOT_ATTRIBUTE = new RegExp(`\\s(?:${COMPOSITION_SCHEMA_ATTRIBUTE}|${COMPOSITION_METADATA_ATTRIBUTE})=`);
// The transform services/config.js applies to each loaded Result of a Session
// before it commits (R13 port): a Result without composition metadata gets the
// legacy composition when its set has one, and a Result with metadata a legend
// side of its mode (OV-672). The saved bytes carry the saved edits. A Session
// 40+ Result takes it only when `appliesToContent` (OV-273).
/** @type {LegacyResultSvgTransform} */
export const transformLegacyResultSvg = Object.assign((/** @type {SVGSVGElement} */ svg, /** @type {LegacyResultSvgData} */ { composition, mode, repairs, path }) => {
  if (
    svg.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) === null
    && svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE) === null
  ) {
    if (!composition) return false;
    normalizeLegacyComposition(svg, composition);
    return true;
  }
  return normalizeStoredLegendSide(svg, mode, { repairs, path });
}, {
  /** @param {unknown} content @param {'circular' | 'linear'} mode @param {boolean} composes */
  appliesToContent: (content, mode, composes) => {
    const rootTag = /<svg\b[^>]*>/i.exec(String(content ?? ''))?.[0] || '';
    return COMPOSITION_ROOT_ATTRIBUTE.test(rootTag) ? !storedLegendSideIsOfMode(rootTag, mode) : composes;
  }
});

// The shown Result as the region tools (Select region, New drawing from region) and Sequence search read it.
/** The mounted root of the shown Result, which the band of Select region and the outline of search hits follow. */
const shownResultRoot = () => state.svgContainer.value?.querySelector('svg') ?? null;
// The Worker helper that plans a Linear request's drawn units (`planLinearUnitSearch`).
/** @param {LinearUnitPlanPayload} payload */
const planLinearUnits = (payload) => runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.PLAN_LINEAR_UNITS, payload);
/**
 * The E-value scales of a drawing's LOSAT searches: Python plans the drawing's draft request, as Save writes it, staged
 * as Generate stages its request (PD-OI-018 revision 6).
 * @type {PlanGenomeScales}
 */
const planDraftGenomeScales = async (drawing, { program, ...facts }) => {
  const { envelope: { renderRequest, resources } } = await projectDraftRequest(state, drawing, { ...facts, planLinearUnits });
  return (await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.PLAN_GENOME_SCALES, {
    canonicalRequest: renderRequest, resources: /** @type {GenerationPayload['resources']} */ (resources), program
  })).result;
};
// The Worker helper that plans a Linear card's window numbers and omitted lengths (PD-OI-124).
/** @param {{ records: CardWindowRecordInput[] }} payload */
const planCardWindows = (payload) => runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.PLAN_CARD_WINDOWS, payload);
/** The records of every Result of the shown drawing whose axis Python wrote (the subjects of Sequence search). */
const shownResultsAxisRecords = () => Array.from({ length: state.results.value.length },
  (_result, index) => readRegionAxisRecords(state.trackSlotResolvedGeometry.value, index)).flat();
const closeRegionPopup = () => {
  if (previewItemPopupOf(state.clickedPairwiseMatch.value)?.previewItem.kind === 'region') state.clickedPairwiseMatch.value = null;
};
// The bases of spans of the shown drawing's records, read by Python from the committed request, the one the shown
// Result was made from: Copy FASTA and Use as search query of Select region, and the subjects of Sequence search.
/** @param {readonly RecordSpanRead[]} spans */
const readShownRecordSpans = (spans) => requestRecordSpanSequences(
  getCommittedCanonicalSession(), spans,
  (payload) => runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.READ_RECORD_SPAN_SEQUENCES, { ...payload,
    resources: /** @type {GenerationPayload['resources']} */ (payload.resources) })
);
/** A rule synthesis helper of the diagram Worker: the `result` it answers for one `request`. @param {string} operation */
const ruleSynthesisHelper = (operation) => async (/** @type {object} */ request) => (await runDiagramHelperOperation(operation, { request })).result;
// The sequences the catalog does not store, read by Python from the committed request the shown Result was made from
// (D5): the Sequence tab, member Copy and Download, the sequence fields of Feature search, the match popup's Copy and
// FASTA, and Interactive SVG export. The one door keeps each reply while the same Result is shown.
/** @type {Parameters<typeof catalogSequenceHelpers>[0]} */
const runCatalogSequenceHelper = (operation, payload) => runDiagramHelperOperation(operation, { ...payload, resources: /** @type {GenerationPayload['resources']} */ (payload.resources) });
const catalogSequences = createCatalogSequences({ getCommittedSession: getCommittedCanonicalSession, helpers: catalogSequenceHelpers(runCatalogSequenceHelper, DIAGRAM_HELPER_OPERATIONS) });
const sequenceReads = createSequenceReads({ readFeatures: catalogSequences.readFeatures, readSpans: catalogSequences.readSpans,
  copyPending: copyPendingTextToClipboard, downloadText: downloadTextFile, isShown: (popup) => state.clickedFeature.value === popup });
// The annotation sets of the committed request the shown Result was made from (Sequence search names a new set apart).
const shownCommittedSets = () => getCommittedCanonicalSession()?.renderRequest.diagramOptions.annotations?.sets ?? [];
// The features (identity and name) the shown Result's admitted catalog holds under a committed record key, or null.
/** @param {string} recordKey @returns {RegionSourceFeature[] | null} */
const shownResultFeatures = (recordKey) => {
  const catalog = toRaw(state.featureCatalog.value);
  const rows = catalog ? catalogNameIndex(catalog).rowsByRecord.get(recordKey) : undefined;
  if (!rows) return null;
  const seen = new Set();
  return rows.flatMap((row) => {
    const featureId = String(row.biologicalFeatureId ?? '');
    if (!featureId || seen.has(featureId)) return [];
    seen.add(featureId);
    return [{ featureId, label: featureNameParts({ biologicalFeatureId: featureId }, row).label }];
  });
};

// The two setting-note readers that the Live setting notes inject: whether a note applies and whether the pending
// layout change is to One circle.
const provideSettingNoteReaders = () => {
  provide(SETTING_NOTE_APPLIES, committedSettingNotePredicate);
  provide(SETTING_NOTE_ONE_CIRCLE_CHANGE, committedOneCircleChange);
};

export const createAppSetup = () => {
  provideSettingNoteReaders();
  setUnmanagedConfigOverrideValidator((payload) => runDiagramHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.VALIDATE_CONFIG_OVERRIDES,
    payload
  ));
  setStoredProteinResultFrameConverter((payload) => runDiagramHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.CONVERT_STORED_PROTEIN_RESULT_FRAMES,
    payload
  ));
  setLegacySessionUpgrader((payload) => runDiagramHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.UPGRADE_LEGACY_SESSION,
    payload
  ));
  const {
    processing,
    processingStatus,
    labelCapacityNotices,
    sessionSavePending,
    sessionImportPending,
    generationCancelRequested,
    errorLog,
    sessionTitle,
    results,
    selectedResultIndex,
    failedGeneratePreservedResult,
    generationFailureRecovery,
    resultPanelTab,
    lastRunInfo,
    annotationWarnings,
    comparisonWarnings,
    labelFitReports,
    svgContent,
    svgResultIdentity,
    zoom,
    layoutRepositionMode,
    isPanning,
    canvasPan,
    canvasContainerRef,
    mode,
    selectedAnnotation,
    losatCacheInfo,
    losatThreadingStatus,
    orthogroups,
    featureOrthogroupIndex,
    selectedOrthogroupId,
    orthogroupSearch,
    orthogroupSortMode,
    showRightDrawer,
    rightDrawerTab,
    linearReorderNotice,
    linearFoldedSources,
    paletteDefinitions,
    paletteNames,
    paletteInstantPreviewEnabled,
    losatExecution,
    richFeaturePopup,
    appliedPaletteName,
    appliedPaletteColors,
    newSpecRule,
    specificRulePresets,
    specificRuleQualifierSuggestions,
    selectedSpecificPreset,
    specificRulePresetLoading,
    downloadDpi,
    downloadWhiteBackground,
    extractedFeatures,
    biologicalFeatures,
    selectedFeatureIds,
    selectedFeatureAnchorId,
    featureSelectionStatus,
    featureSelectionDrag,
    selectedFeatureCount,
    selectedFeatures,
    hasFeatureSelection,
    featureEditorStatus,
    featureEditorStatusText,
    featureExtractionPending,
    featureRecordIds,
    selectedFeatureRecordIdx,
    featurePanelTab,
    featureSearchInput,
    featureSearch,
    previewFeatureSearchInput,
    previewFeatureSearchQuery,
    previewFeatureSearchField,
    previewFeatureSearchQualifierKey,
    previewFeatureSearchUseRegex,
    previewFeatureSearchMatches,
    previewFeatureSearchMatchDetails,
    previewFeatureSearchActiveIndex,
    previewFeatureSearchError,
    previewFeatureSearchRenderedCount,
    featureListScrollTop,
    featureListViewportHeight,
    isFeatureDrawerMounted,
    visibleFeatureRows,
    featureRecordPickerVisible,
    featureListTopSpacerPx,
    featureListExpandedId,
    featureListBottomSpacerPx,
    labelSearch,
    editableLabels,
    filteredEditableLabels,
    autoLabelReflowEnabled,
    labelReflowProcessing,
    labelReflowLastError,
    svgContainer,
    clickedFeature,
    clickedFeaturePos,
    clickedPairwiseMatch,
    clickedPairwiseMatchPos,
    pairwiseMatchPopupRef,
    pairwiseMatchPopupDrag,
    pairwiseMatchPopupSize,
    pairwiseMatchPopupResize,
    featurePopupRef,
    featurePopupDrag,
    featurePopupSize,
    featurePopupResize,
    featureStyleScopeDialog,
    featureVisibilityScopeDialog,
    legendRenameDialog,
    paletteColorsDialog,
    resetColorDialog,
    labelTextScopeDialog,
    hiddenLabelTextDialog,
    labelOnDialog,
    sidebarWidth,
    originalLegendOrder,
    showCanvasControls,
    skipCaptureBaseConfig,
    featureKeys,
    defaultColorKeys,
    newColorFeat,
    newColorVal,
    newPriorityRule,
    newFeatureToAdd,
    filteredFeatures,
    featureListState
  } = state;
  /** @type {ReturnType<typeof createSimilarityAlignmentActions> | null} */
  let similarityAlignmentActions = null;
  // R13 port: the drawer, the preview binder, and the result watchers read the
  // alignment owner through these ports; the root assigns them once the
  // alignment owner exists.
  /** @type {{ refreshCanvas: () => void, editorDisabledReason: () => string }} */
  const similarityAlignmentPorts = {
    refreshCanvas: () => {},
    editorDisabledReason: () => ''
  };
  /**
   * The shown drawing when it is of `drawingMode`, else null: a mode panel's
   * handler does nothing for another mode's drawing. Resolved at each call,
   * never cached.
   * @param {DiagramMode} drawingMode
   * @returns {DrawingState | null}
   */
  const shownDrawingOfMode = (drawingMode) => {
    const drawing = state.activeDrawing();
    return drawing.mode === drawingMode ? drawing : null;
  };
  // The record rows of the shown Linear drawing, none when another drawing is
  // shown: the File cards and the record names its panel lists.
  /** @returns {readonly LinearSeq[]} */
  const shownLinearSeqs = () => {
    const drawing = shownDrawingOfMode('linear');
    return drawing === null ? [] : drawing.linearSeqs;
  };
  /**
   * A mode panel's mutating handler: the session operation's busy result
   * first, then nothing when the shown drawing is not of `drawingMode`, else
   * `action` on the shown drawing. One owner of the admission order.
   * @template T
   * @param {DiagramMode} drawingMode
   * @param {(drawing: DrawingState) => T} action
   * @returns {T | import('../state.js').SessionBusy | undefined}
   */
  const shownAction = (drawingMode, action) => {
    const busy = sessionOperationAvailability();
    if (busy) return busy;
    const drawing = shownDrawingOfMode(drawingMode);
    return drawing ? action(drawing) : undefined;
  };
  /**
   * The drawing the displayed Result belongs to when it is the Linear one,
   * else null (the group names are the Linear drawing's). Resolved at each call.
   * @returns {DrawingState | null}
   */
  const displayedLinearDrawing = () => {
    const drawing = historySnapshots.displayedArtifactDrawing();
    return drawing.mode === 'linear' ? drawing : null;
  };
  const linearTypography = createLinearTypographyController({
    shownDrawing: () => shownDrawingOfMode('linear'),
    mutationAvailability: sessionOperationAvailability
  });

  const {
    comparisonHeightValidationError, invalidateLinearComparisonArtifacts, linearComparisonUi,
    linearComparisonGlobalAction, canRunLinearLosat
  } = createLinearComparisonView({ state, mode, losatCacheInfo, writeLosatCacheInfo, shownDrawingOfMode });
  const {
    linearSourceGroups, linearSourceRemovalDialog, linearSourceRemovalReturnFocus, linearSourceRemovalTarget,
    linearSourceRemovalCanDelete, linearSourceFileName, linearCardFileLine, shownInputComplete, linearCardHoldsInput,
    isLinearSourceFolded, syncLinearSourceFold, linearSourceRemovalTargetName
  } = createLinearFileCards({ shownDrawingOfMode, shownLinearSeqs, linearFoldedSources });
  const {
    linearComparisonTimeline, linearComparisonPairForEdgeKey, linearComparisonRecordLabel, linearComparisonRecordName
  } = createLinearComparisonTimeline({ state, shownLinearSeqs });
  // The settings pane's display, and the one way to bring a setting into
  // view: every "show this setting" path below ends in `revealSetting`, so a
  // setting of a section the icon rail does not show is docked first.
  const settingsPaneActions = createSettingsPaneController({ state });
  const {
    openLinearComparisonDisclosure, focusLinearComparisonPair, focusLinearComparisonIssue, chooseLinearComparisonPairFile,
    comparisonCauseShown, showComparisonCause
  } = createLinearComparisonFocus({
    state, mode, errorLog, linearComparisonUi, revealSetting: settingsPaneActions.revealSetting,
    revealSettingsSection: settingsPaneActions.revealSettingsSection
  });
  const { shownLinearRecordRowFor, linearLayoutTokens, linearLosatCacheInfoByEdgeKey, linearLosatWindowPairs, linearWindowPairSide } = createLinearLayoutReaders({
    state, shownDrawingOfMode, losatCacheInfo
  });

  const {
    pendingDiscoveryOf, requestRecordListAfterUpload, listedDrawingOf, readInputRecords
  } = createRecordDiscoveryRequests({ state });
  /**
   * @param {DrawingState} drawing
   * @param {{ uid: string, records: SequenceRecord[], inferDefinitions?: boolean }} discovered
   */
  const expandDiscoveredLinearRecords = (drawing, { uid, records, inferDefinitions = false }) => {
    const { expansions } = pendingDiscoveryOf(drawing);
    const namesOfThisFile = expansions.get(uid);
    expansions.delete(uid);
    const index = drawing.linearSeqs.findIndex((seq) => seq.uid === uid);
    if (index < 0) return;
    const source = drawing.linearSeqs[index];
    if (namesOfThisFile === undefined || records.length < 2 || source.region_record_id) {
      if (inferDefinitions) {
        source.inferred_definition = inferredDefinitionForRecord(records, source.region_record_id);
      }
      return;
    }
    const row = linearRecordRowFor(drawing, uid, index + 1);
    // Expanded records keep only File-level values. A crop or record display
    // value of the replaced File does not apply to the new records (IN-04).
    // A name entered on the card of this file before it had several records
    // becomes the File's name where the File has none (OV-486).
    const typed = namesOfThisFile ? source : null;
    const expanded = buildDisambiguatedRecordEntries(records).map((record, recordIndex) => createLinearSeq({
      uid: recordIndex === 0 ? uid : undefined,
      gb: source.gb,
      gff: source.gff,
      fasta: source.fasta,
      depth: source.depth,
      losat_gencode: source.losat_gencode,
      file_definition: source.file_definition || typed?.definition || '',
      file_subtitle: source.file_subtitle || typed?.record_subtitle || '',
      inferred_definition: inferDefinitions ? record.inferredDefinition || '' : '',
      region_record_id: record.value
    }));
    applyLinearSeqMutation(drawing, [
      ...drawing.linearSeqs.slice(0, index), ...expanded, ...drawing.linearSeqs.slice(index + 1)
    ], { alignmentMutation: msg('record selector changed.') });
    expanded.forEach((seq) => updateLinearRecordRow(drawing.linearRecordRows, seq.uid, row));
    return true;
  };
  /** @param {{ drawingId: DrawingId, uid: string, records: SequenceRecord[] }} discovered the records of one card, for the drawing the discovery ran for */
  const handleLinearRecordsDiscovered = ({ drawingId, uid, records }) => {
    const drawing = listedDrawingOf(drawingId);
    if (!drawing) return false;
    // Only an upload infers record definitions and opens the record list; a
    // loaded Session keeps its own.
    const fromUpload = !state.sessionImportRollbackInProgress.value && pendingDiscoveryOf(drawing).uploads.delete(uid);
    const expanded = expandDiscoveredLinearRecords(drawing, {
      uid, records, inferDefinitions: fromUpload && drawing.inputType.value === 'gb'
    });
    if (fromUpload) requestRecordListAfterUpload(drawing, 'linear', uid, records.length);
    // After the expansion, which gives records 2..n of the file new keys: an
    // expansion rebinds in its source edit (`applyLinearSeqMutation`).
    if (!expanded) rebindAnnotationRecords();
    return expanded;
  };
  const paletteLoader = createPaletteLoader({ state });
  const linearRecordSelector = createLinearRecordSelector({
    state,
    reactive,
    onRecordsDiscovered: handleLinearRecordsDiscovered,
    recordReader: ({ inputType, primaryFile, pairedFile }) => readInputRecords(inputType, primaryFile, pairedFile)
  });
  const {
    getCircularRecordDiscoveryState, circularAnnotationSourceKey, linearAnnotationSource, getAnnotationRecordCatalog,
    recordCatalogError, prepareLinearRecordCatalog
  } = createDrawingRecordCatalog({
    state, readInputRecords, linearRecordCatalogStateFor: linearRecordSelector.catalogStateFor,
    linearRecordErrorModelFor: linearRecordSelector.errorModelFor, refreshLinearRecords: linearRecordSelector.refresh
  });
  const { rebindAnnotationRecords } = createAnnotationRecordRebind({
    state, getAnnotationRecordCatalog, sessionOperationAvailability
  });
  const {
    handleCircularRecordsDiscovered, materializeAutomaticLinearRecords, prepareLinearGenerate
  } = createRecordDiscoveryRuns({
    pendingDiscoveryOf, requestRecordListAfterUpload, listedDrawingOf, rebindAnnotationRecords,
    refreshLinearRecords: linearRecordSelector.refresh
  });
  // Bumped by each edit of the displayed Result's root (R1).
  const resultEditRevision = shallowRef(0);
  // Counts the Result roots the binder made ready: a History restore can mount a root whose content equals the
  // content the preview already held (its live edits were never serialized), which neither ref above shows.
  const mountedResultRevision = shallowRef(0);
  const previewRuntime = createPreviewRuntime({ state, serializeSvg: serializeCleanSvg, editRevision: resultEditRevision });
  setPreviewRuntime({
    clearActiveRuntime: previewRuntime.clearActiveRuntime,
    currentResults: previewRuntime.currentResults
  });

  // Before the History snapshot service, which clears its copy feedback (D-37).
  const orthogroupActions = createOrthogroupEditor({ state, getCommittedRequest: getCommittedCanonicalRenderRequest, memberAction: sequenceReads.memberAction });
  const historySnapshots = createHistorySnapshotService({
    state,
    fileStore: state.inputFiles,
    nextTick,
    // The discovery owner exists once the runner does; a restore runs after setup.
    discoverDrawingRecords: (drawing) => discoverDrawingRecords(drawing),
    buildConfigData,
    applyConfigData: applyCapturedConfigData,
    buildUiStateData,
    applyUiStateData,
    applySessionUiData,
    buildFeatureStateData,
    applyFeatureStateData,
    buildEditorStateData,
    applyEditorStateData,
    buildOrthogroupStateData,
    applyOrthogroupStateData,
    currentResults: previewRuntime.currentResults,
    applyResultsData,
    buildRunStateData,
    applyRunStateData,
    clearSimilarityGroupCopyFeedback: orthogroupActions.clearCopyFeedback
  });
  // R13: a History restore keeps the open specific-color pattern drafts, and
  // Undo and Redo of a rule edit follow the restored rules as the edit did
  // (OV-43): the rule owner asks for the rerender when they change a Legend
  // source. A Generate restore brings the Results drawn with its rules. The
  // feature editor takes History's undoable runs, so the root registers its
  // ports once both owners exist.
  /** @type {{
   *   captureSpecificRulePatternDrafts: ReturnType<typeof createFeatureEditor>['captureSpecificRulePatternDrafts'],
   *   restoreSpecificRulePatternDrafts: ReturnType<typeof createFeatureEditor>['restoreSpecificRulePatternDrafts'],
   *   followRestoredSpecificRules: ReturnType<typeof createFeatureEditor>['followRestoredSpecificRules'],
   *   retainRulesForRestore: RulePreparation['retain']
   * }} */
  // Each field is null only until the root registers it below (R13); no restore
  // runs before then, so the declared function types hold at every read.
  const specificRuleRestorePorts = {
    captureSpecificRulePatternDrafts: /** @type {never} */ (null),
    restoreSpecificRulePatternDrafts: /** @type {never} */ (null),
    followRestoredSpecificRules: /** @type {never} */ (null),
    retainRulesForRestore: /** @type {never} */ (null)
  };
  /**
   * @template {unknown[]} A
   * @template R
   * @param {(...args: A) => R} restore
   * @param {A} args
   */
  const restoreWithSpecificRuleDrafts = async (restore, ...args) => {
    const drafts = specificRuleRestorePorts.captureSpecificRulePatternDrafts();
    try {
      return await restore(...args);
    } finally {
      if (drafts) specificRuleRestorePorts.restoreSpecificRulePatternDrafts(drafts);
    }
  };
  // A restore, or Reset Settings, that rewrites the rules: the rule owner
  // follows the rules `write` replaced. A checkpoint figure History recorded
  // as drawn shows its rules as Python draws them, so its restore asks for no
  // rerender (D-24 (OV-731)); one recorded as not drawn is drawn again by
  // History's own one request (D-174), so the rule owner adds none.
  /**
   * @template T
   * @param {DrawingState} drawing
   * @param {() => Promise<T>} write
   * @param {{ drawn?: boolean | undefined }} [restored] `drawn`: History's record of a checkpoint figure; unset when unknown.
   * @returns {Promise<T>}
   */
  const followRuleEdits = async (drawing, write, { drawn } = {}) => {
    // The rule owner follows the shown drawing; a drawing that is not shown
    // follows its rules when it is shown again (D-15-51).
    if (drawing !== state.activeDrawing()) return write();
    const previousRules = drawing.manualSpecificRules.map((/** @type {SpecificColorRule} */ rule) => ({ ...rule }));
    specificRuleRestorePorts.retainRulesForRestore(previousRules);
    try {
      const written = await write();
      if (drawn === undefined) specificRuleRestorePorts.followRestoredSpecificRules(previousRules);
      return written;
    } finally {
      specificRuleRestorePorts.retainRulesForRestore([]);
    }
  };
  // An operation's outcome shows in the notice banner (no native alert).
  const noticeBanner = createNoticeBanner({ reactive });
  const showNotice = noticeBanner.show;
  // PD-OI-090: one History per drawing behind the router, which keeps the
  // manager API for the shown drawing and the project stream (the
  // Session-level values, Reset all drawings); the drawings share one byte
  // budget and one file retention (OV-112). The router owns its members: these
  // are the drawings at start, and `managerOptions` makes the History of any
  // drawing the store holds (R11).
  const history = createHistoryRouter({
    drawingIds: state.drawingList().map((drawing) => drawing.id),
    activeDrawingId: () => state.activeDrawingId.value,
    sessionHistory: historySnapshots.sessionHistory(),
    // OV-1260: dropped older steps are told, not silent.
    onHistoryLimit: (text) => showNotice({ text: tm(text), tone: 'warning' }),
    makeRef: ref,
    estimateFileBytes: (fileIds) => state.inputFiles.estimateBytes(fileIds, inputFilePayload),
    // A file stays stored while any History step or binding names it (R11).
    retainOnlyFiles: retainHistoryAndBoundFiles(state, state.inputFiles),
    // A History finds its drawing by ID at each call, never a held object: a
    // legacy configuration Load commits a copy of the shown drawing under its
    // ID and keeps that drawing's History (OV-1790, D-15-114).
    managerOptions: (drawingId) => {
      const drawing = () => state.drawingOf(drawingId);
      const ports = historySnapshots.drawingHistory(drawingId);
      return /** @type {HistoryManagerOptions} */ ({
        buildIntent: ports.buildIntent,
        applyIntent: (...args) => followRuleEdits(
          drawing(), () => restoreWithSpecificRuleDrafts(ports.applyIntent, ...args)
        ),
        buildCheckpoint: ports.buildCheckpoint,
        applyCheckpoint: (checkpoint, restored) => followRuleEdits(
          drawing(), () => restoreWithSpecificRuleDrafts(ports.applyCheckpoint, checkpoint), restored
        ),
        captureGeneratedArtifactHandle: () => historySnapshots.captureGeneratedArtifactHandle(drawing()),
        restoreGeneratedArtifactHandle: (...args) => restoreWithSpecificRuleDrafts(historySnapshots.restoreGeneratedArtifactHandle, ...args),
        signatureFor: historySnapshots.snapshotSignature,
        collectCurrentFileIds: historySnapshots.collectCurrentFileIds,
        identityMapPath: ports.identityMapPath,
        editedPartOf: ports.editedPartOf,
        makeRef: ref,
        mutationAvailability: sessionOperationAvailability
      });
    }
  });
  // The choice in flight of a popup choice dialog: its busy state and its close (D-12).
  const dialogChoice = createDialogChoice({ mutationPending: history.mutationPending, runUndoable: history.runUndoable, ref });
  // Confirmations ask with the app's choice dialog (no native confirm).
  const choicePrompt = createChoicePrompt({ reactive });
  const askChoice = choicePrompt.ask;
  const askText = choicePrompt.askText;
  // The Font section: the drawing's text font and the fonts loaded in this page.
  const fontSettings = createFontSettings({
    state,
    runUndoable: history.runUndoable,
    registerFontsInWorker: async (fonts) => (await runDiagramHelperOperation(
      DIAGRAM_HELPER_OPERATIONS.REGISTER_FONTS, { fonts }
    )).result,
    askChoice,
    committedFontFamily: () => requestFontFamily(getCommittedCanonicalRenderRequest())
  });
  // A new diagram Worker registers the page's fonts, and a Generate names
  // those it uses; a font the new Worker could not register again is named
  // (OV-686), and loading it again registers it again (OV-740).
  setDiagramWorkerFontSource({
    readFonts: () => state.loadedFonts.value,
    onFontsLost: (lost) => {
      fontSettings.markFontsLost(lost.map((font) => font.identity));
      showNotice({
        text: t('Fonts loaded in this page could not be prepared again for drawing. Load them again in Font, or choose another font.'),
        items: lost.map((font) => font.label)
      });
    }
  });
  /**
   * A reader of one card of the shown drawing, for the template.
   * @template T
   * @param {(drawing: DrawingState, seq: LinearSeq) => T} read
   * @returns {(seq: LinearSeq) => T}
   */
  const readingShownLinearCard = (read) => (seq) => read(state.activeDrawing(), seq);
  const recordDisplayControls = createRecordDisplayControls({ state, computed, watch,
    linearRecordsFor: linearRecordSelector.recordsFor, linearWindowsOf: linearDrawingWindowsOf,
    linearRecordStatusFor: linearRecordSelector.statusFor,
    linearRecordErrorFor: linearRecordSelector.errorFor,
    runUndoable: history.runUndoable,
    getCommittedRequest: getCommittedCanonicalRenderRequest, getCommittedSession: getCommittedCanonicalSession });
  // Every write of a record source (a Circular file, an input type, the Linear
  // cards) runs through it, so a replaced source's display drafts go.
  const { changeRecordSources } = recordDisplayControls;
  // Orient records from comparisons: Python decides from the last Generate, the record display writes the flags.
  const recordAutoOrientation = createRecordAutoOrientation({ computed, ref,
    getCommittedSession: getCommittedCanonicalSession,
    readCommittedRevision: readCommittedCanonicalRevision,
    busy: () => sessionOperationAvailability(),
    // The committed Session's resource table holds the descriptors its writer wrote.
    decideOrientation: (payload) => runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.DECIDE_AUTO_ORIENTATION, { ...payload,
      resources: /** @type {GenerationPayload['resources']} */ (payload.resources) }),
    writeDecidedOrientations: recordDisplayControls.writeDecidedOrientations,
    errorMessage: (error) => errorText(normalizeCaughtError(error, { operation: 'decideAutoOrientation', stage: 'helper' }).summaryParts) });
  window.__GBDRAW_HISTORY__ = history;
  // R13: owners receive this port, not the record display controls. The check
  // is pure (services/feature-identity.js) over the record display's binding.
  /** @param {RecordKeyedFeature} feature */
  const isCurrentResultFeature = (feature) => isCurrentFeature(feature, recordDisplayControls.sourceBinding());
  // Undo, Redo, and the History list (PD-OI-112).
  const { canUndoHistory, canRedoHistory, undoHistoryTitle, redoHistoryTitle, undoHistory, redoHistory, historyList } = createHistoryControls({
    revision: () => history.revision.value, canUndo: history.canUndo, canRedo: history.canRedo, undo: history.undo, redo: history.redo,
    timeline: history.timeline, travelTo: history.travelTo, computed, ref, nextTick
  });
  // D-03: the Rich Feature Popup preference records no Undo step, but a change
  // of it is unsaved work that Load Session asks about (D-07 (b)).
  const markRichFeaturePopupChanged = () => history.markChangeOutsideHistory();
  const { unmanagedConfigOverrideEntries, resetUnmanagedConfigOverride } = createPreservedConfigOverrides({
    state, runUndoable: history.runUndoable
  });
  /** @param {unknown} index */
  const selectResult = (index) => previewRuntime.selectResult(index);

  const {
    handleWheel,
    startPan,
    doPan,
    endPan,
    resetPreviewViewport,
    fitPreviewToViewport,
    cancelPreviewTransformInteraction,
    disposePanZoom,
    isPreviewTransformActive, subscribePreviewTransformChanges,
    holdPreviewViewport,
    releasePreviewViewport,
    heldPreviewViewport, hidePreviewUntilShown, showHiddenPreview, storedFigurePreview
  } = createPanZoom(state);
  const { startResizing } = createSidebarResize(state);

  const { ruleMatchingPending } = state;
  // An edit's rule preparation: Save and Load wait for it, not for a popup's (OV-377).
  const ruleMatchingBlocking = ref(false);
  // Python's rule evaluation (R7): the rule preparation's, and the Label TSV
  // import's own stateless one (`evaluateLabelRules`).
  /** @type {RulePreparationOptions['evaluate']} */
  const evaluateRules = async (payload, options) => (await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES, payload, options)).result;
  const rulePreparation = createRulePreparation({
    state,
    pending: ruleMatchingPending,
    blocking: ruleMatchingBlocking,
    evaluate: evaluateRules,
    synthesize: ruleSynthesisHelper(DIAGRAM_HELPER_OPERATIONS.SYNTHESIZE_COLOR_RULES), synthesisQualifiers: ruleSynthesisHelper(DIAGRAM_HELPER_OPERATIONS.SYNTHESIS_QUALIFIERS),
    visibilityRules: () => requestFeatureVisibilityRules(state.activeDrawing().featureVisibilityManualRules)
  });
  specificRuleRestorePorts.retainRulesForRestore = rulePreparation.retain;
  // R13: a Legend row a specific-color rule draws commits its edit through the
  // rule owner; the root registers the port once the feature editor exists.
  /** @type {{ commitLegendRowRules: CommitSpecificRulesPort }} */
  // Null only until the root registers the port below; the Legend cannot commit a row before that.
  const legendRowRulePorts = { commitLegendRowRules: /** @type {never} */ (null) };
  const legendActions = createLegendManager({
    state,
    commitLegendRowRules: (...args) => legendRowRulePorts.commitLegendRowRules(...args),
    beginHistoryTransaction: history.begin,
    commitHistoryTransaction: history.commit,
    commitActiveResultEdit: previewRuntime.commitActiveResultEdit,
    readActiveResultIdentity: () => previewRuntime.getActiveRuntime()?.resultIdentity,
    readLegendRowColors: (drawing) => readLegendRowColors(drawing),
    readHiddenTrackRows: (drawing) => hiddenTrackLegendRows(drawing),
    legendRowNames: shallowRef(new Map()),
    showLegendStructure: () => showEditorIntent({ domains: LIVE_EDIT_DOMAINS.legendStructure })
  });
  // R13: the palette watcher reacts through the root's palette and rules
  // projection, registered once the style owner it applies through exists.
  /** @type {{ projectPaletteAndRules: FeatureEditorOptions['projectPaletteAndRules'] }} */
  // Null only until the root registers the port below; the palette watcher reads it after setup.
  const paletteRulePorts = { projectPaletteAndRules: /** @type {never} */ (null) };
  // OV-484: the tracks of the request that drew the displayed Result, the one
  // reader of a track's slot and colors for the Legend panel, the live compile,
  // the feature popup and the palette's track painter. A Result loaded without
  // a request (a Session 27 to 30) is described by the track settings Load recorded.
  const readDrawnTracks = () => {
    const request = getCommittedCanonicalRenderRequest();
    return request ? drawnTracks(request) : loadedTracks(getLoadedTrackSettings());
  };
  // The Legend rows Python drew for the displayed Result's tracks, by key.
  /** @param {DrawingState} drawing */
  const displayedTrackLegendRows = (drawing) => {
    const { pythonRows, sharedLegendKeys } = displayedLegendRowContext(state, drawing);
    return trackLegendRows(readDrawnTracks(), pythonRows, sharedLegendKeys);
  };
  // The keys of the displayed Result's Legend rows whose drawn track the draft
  // hides (OV-1198, OV-1215), by the hide answer the live track projection
  // reads: the Legend list leaves them out, so a hidden track's rows leave and
  // return with it. A Result without a request, and one whose live track apply
  // is refused, keeps every row until Generate.
  /** @param {DrawingState} drawing @returns {string[]} */
  const hiddenTrackLegendRows = (drawing) => hiddenTrackLegendKeys(
    displayedTrackLegendRows(drawing), getCommittedCanonicalRenderRequest(), committedTrackHideDrawing(drawing)
  );
  const svgActions = createSvgStyles({
    state,
    watch,
    nextTick,
    commitActiveResultEdit: previewRuntime.commitActiveResultEdit,
    projectPaletteAndRules: (...args) => paletteRulePorts.projectPaletteAndRules(...args),
    readDrawnTracks
  });
  // OV-276, OV-282 (D-26): a Legend panel row shows the color the palette
  // gives its swatch on the displayed Result (the live compile's, by Python's
  // rows), derived and never written into the rows. A new Result brings new
  // Python rows. The Legend and rule owners read the same cached colors
  // through `readLegendRowColors`; a drawing that is not shown has no
  // displayed Result, and its rows show the colors they record.
  const legendRowColors = computed(() => displayedLegendRowColors(state, state.activeDrawing(), readDrawnTracks()));
  /** @type {import('./candidate-render.js').ReadLegendRowColorsPort} */
  const readLegendRowColors = (drawing) => (drawing === state.activeDrawing() ? legendRowColors.value : RECORDED_LEGEND_ROW_COLORS);
  /** @param {Parameters<LegendRowColorOf>[0]} entry */
  const legendEntryColor = (entry) => legendRowColors.value.listed(entry);
  /** @param {Parameters<LegendRowColorOf>[0]} entry */
  const deletedLegendEntryColor = (entry) => legendRowColors.value.deleted(entry);
  // The palette and the specific-color rules on the mounted Result (R3): the
  // one call of their projection, shared by `projectMountedEditorIntent`, the
  // palette watcher, and a rule commit. It prepares the rule matches first
  // unless the caller prepared them (`prepareRules: false`, which applies at
  // once), and resolves to false when the rules changed meanwhile or a Session
  // operation is in progress (the Result stays as it is, tracks and features
  // alike). The style owner paints the tracks; feature and Legend row fills
  // show through the port, compiled with Generate's precedence, unless the
  // caller shows them in its own compile (`show: false`). The Save that shows
  // what a display declined because of it (`saving`) does not count itself.
  // A rule commit's Legend rows (`legendRows`) show in the same compile.
  /** @param {boolean} saving */
  const sessionBlocks = (saving) => Boolean(sessionOperationAvailability()) && !(saving && sessionSavePending.value);
  // The domains a display of the feature fills shows. A Legend row's stroke
  // reaches the features drawn in the row's color (`legendRowFeatureIds`), so a
  // fill edit can change which features it reaches, and a rule commit can move
  // the row to another key: while the drawing holds a Legend row stroke, a fill
  // display also shows the strokes (OV-1985, OV-1986).
  /** @param {DrawingState} drawing @returns {readonly string[]} */
  const fillDomainsOf = (drawing) => (Object.keys(drawing.legendStrokeOverrides).length > 0
    ? [...LIVE_EDIT_DOMAINS.fills, ...LIVE_EDIT_DOMAINS.strokes]
    : LIVE_EDIT_DOMAINS.fills);
  /**
   * @param {{
   *   recolor?: { recolorPairwise?: boolean, recolorCollinear?: boolean }, prepareRules?: boolean, show?: boolean, saving?: boolean,
   *   legendRows?: import('./candidate-render.js').RuleLegendRows | null
   * }} [options]
   */
  const projectPaletteAndRules = ({ recolor = {}, prepareRules = true, show = true, saving = false, legendRows = null } = {}) => {
    const blocked = () => sessionBlocks(saving);
    const project = () => {
      if (blocked()) return false;
      svgActions.applyPaletteToSvg(recolor);
      if (show) {
        const fills = fillDomainsOf(state.activeDrawing());
        showEditorIntent(legendRows
          ? { domains: [...fills, ...LIVE_EDIT_DOMAINS.legendStructure], ruleRows: legendRows }
          : { domains: fills });
      }
      return true;
    };
    return prepareRules
      ? Promise.resolve(rulePreparation.prepare(undefined, { blocked })).then((prepared) => (prepared ? project() : false))
      : project();
  };
  paletteRulePorts.projectPaletteAndRules = projectPaletteAndRules;
  // What the Review alignment panel took over while it is shown, given back when it ends: the
  // Editor drawer reopens on its tab, the popups (hidden meanwhile) show again, and the view
  // from before the fit returns. While it is set, no gesture on the figure and no click
  // outside opens, retargets or closes a popup. The one predicate of that hold: the figure's
  // gestures read it, and the preview's `preview-reviewing` class (the hint and the hand cursor
  // of a drawn element) is bound to the same computed.
  // `figure` is the identity of the Result of the review's drawing (`drawingId`) that the review last
  // showed its draft on.
  // `drawerTab` is the Editor tab the review closed: Apply's new Result can leave a tab
  // unavailable for a while, which moves the shown tab (OV-1255).
  const similarityAlignmentTakeover = shallowRef(/** @type {{ drawerOpen: boolean, drawerTab: string, popup: unknown, drawingId: DrawingId, figure: string } | null} */ (null));
  // The hold belongs to the review's drawing (OV-1089): while another drawing is shown, its figure
  // is not the held one, so the end of the review does not give it the held view and its
  // captures keep its own view.
  const similarityAlignmentHoldShown = () => {
    const takeover = similarityAlignmentTakeover.value;
    return !takeover || takeover.drawingId === state.activeDrawingId.value;
  };
  // The Session writer and the artifact checkpoint save the view the pan/zoom owner holds for
  // the shown mode. The departing mode's slot is captured before a switch changes `mode`.
  setHeldPreviewViewportReader(() => (similarityAlignmentHoldShown() ? heldPreviewViewport() : null));
  const previewPopupsHeldNow = computed(() => similarityAlignmentTakeover.value !== null);
  const previewPopupsHeld = () => previewPopupsHeldNow.value;
  // The records of the shown Result whose axis Python wrote: Select region and Shift+click range selection read them.
  const shownAxisRecords = () => shownAxisRecordsOf(state);
  const {
    consumeSuppressNextClick, selectFeatureByClick, markPlainFeatureClick, clearFeatureSelection, startMarqueePointer, moveMarqueePointer, commitMarqueePointer, cancelMarqueePointer,
    syncFeatureSelectionClasses, replaceSelectedFeatureIds, featureSelectionMarqueeStyle, featureSelectionToolbarStyle, startToolbarDrag
  } = createFeatureSelection({
    state, onMounted, onUnmounted, previewPopupsHeld,
    axisRecordOf: (recordKey) => shownAxisRecords().find((record) => record.recordKey === recordKey) ?? null
  });
  // The palette owner comes before the feature editor, whose popup sets a
  // type's default color through it (D-15).
  const resultsManager = createResultsManager({
    state, closeAfterDialogChoice: dialogChoice.closeAfterChoice, writeArtifactValues: historySnapshots.writeArtifactValues
  });
  const featureActions = createFeatureEditor({
    state,
    previewPopupsHeld, shownAxisRecords,
    rulePreparation,
    evaluateLabelRules: evaluateRules,
    runUndoable: history.runUndoable,
    runUndoableCheckpoint: history.runUndoableCheckpoint,
    closeAfterDialogChoice: dialogChoice.closeAfterChoice,
    getCommittedRequest: getCommittedCanonicalRenderRequest,
    readDrawnTracks,
    getCommittedSession: getCommittedCanonicalSession,
    readResourceRecordCount: readCommittedResourceRecordCount,
    readFeatureOverrideTable: (payload) => runDiagramHelperOperation(
      DIAGRAM_HELPER_OPERATIONS.READ_FEATURE_OVERRIDE_TABLE, payload
    ),
    isCurrentFeature: isCurrentResultFeature,
    isPatternEditAvailable: () => !sessionImportPending.value,
    nextTick,
    prepareFileLegendEntries: /** @type {any} */ (legendActions.prepareFileLegendEntries),
    moveLegendRecords: legendActions.moveLegendRecords,
    readLegendRowColors,
    consumeSuppressNextClick, selectFeatureByClick, markPlainFeatureClick, clearFeatureSelection, startMarqueePointer, moveMarqueePointer, commitMarqueePointer, cancelMarqueePointer,
    commitActiveResultEdit: previewRuntime.commitActiveResultEdit,
    currentResults: previewRuntime.currentResults,
    getResultIdentity: previewRuntime.getResultIdentity,
    showEditorIntent: (options) => showEditorIntent(options),
    selectResult,
    isPreviewTransformActive, subscribePreviewTransformChanges,
    projectPaletteAndRules,
    projectFeatureEdits: () => projectMountedEditorIntent({ visibility: true, rerender: true, reflow: true, labels: true }),
    readUserDefaultColor: resultsManager.readUserDefaultColor,
    setDefaultColor: resultsManager.setDefaultColor,
    editDefaultColors: resultsManager.editDefaultColors,
    legendRowIsNamed: legendActions.legendRowIsNamed,
    legendRowName: legendActions.legendRowName,
    findLegendRowNamed: legendActions.findLegendRowNamed,
    keepBothNamesOf: legendActions.keepBothNamesOf,
    legendEntryName: legendActions.legendEntryName,
    writeLegendRowColor: legendActions.writeLegendRowColor,
    shownLegendName: legendActions.shownLegendName,
    askChoice,
    showNotice,
    revealSetting: settingsPaneActions.revealSetting,
    showElementSettings: settingsPaneActions.showSettingsBesideResult,
    closeElementSettings: settingsPaneActions.closeSettingsPopup
  });
  legendRowRulePorts.commitLegendRowRules = featureActions.commitSpecificRules;
  specificRuleRestorePorts.captureSpecificRulePatternDrafts = featureActions.captureSpecificRulePatternDrafts;
  specificRuleRestorePorts.restoreSpecificRulePatternDrafts = featureActions.restoreSpecificRulePatternDrafts;
  specificRuleRestorePorts.followRestoredSpecificRules = featureActions.followRestoredSpecificRules;
  // D-24 (OV-731): History's one reader of a figure not drawn yet (OV-2560), and what a restored side brings back (D-174).
  history.setFigurePort({
    drawn: (drawingId) => !featureActions.isRedrawOwed(state.drawingOf(drawingId)),
    restored: (drawingId, drawn) => featureActions.redrawRestored(state.drawingOf(drawingId), drawn)
  });
  // R13: the drawer and the feature search come after the owners they react
  // through, so each receives its ports directly.
  const rightDrawerActions = createRightDrawerController({ state, watch, computed,
    onClose: featureActions.suspendSpecificRulePatternDrafts,
    onOpen: () => settingsPaneActions.showPhonePane('preview'),
    focusReturn: {
      isFocusInDrawer: () => Boolean(document.querySelector('.right-drawer')?.contains(document.activeElement)),
      focusToggle: () => /** @type {HTMLElement | null} */ (document.querySelector('.drawer-toggle'))?.focus()
    },
    getOpenDisabledReason: () => similarityAlignmentPorts.editorDisabledReason() });
  const previewFeatureSearch = createPreviewFeatureSearch({
    state,
    watch,
    nextTick,
    computed,
    isActiveResultReady: previewRuntime.isActiveResultReady,
    resolveOrthogroups: () => {
      const drawing = displayedLinearDrawing();
      return orthogroups.value.map((group) => ({
        ...group,
        display_name: orthogroupActions.resolveOrthogroupName(group, drawing),
        description: orthogroupActions.resolveOrthogroupDescription(group, drawing)
      }));
    },
    openFeatureEditorForFeature: featureActions.openFeatureEditorForFeature, readFeatureSequences: catalogSequences.readFeatures,
    focusSearchQuery: () => /** @type {HTMLElement | null} */ (document.querySelector('.preview-feature-search .gfs-query'))?.focus(),
    focusSearchToggle: () => /** @type {HTMLElement | null} */ (document.querySelector('[data-preview-search-toggle]'))?.focus()
  });

  // A Result switch closes the feature popup through its close owner: the
  // popup belongs to the Result it was opened on (OV-1269).
  watch(selectedResultIndex, () => {
    clearFeatureSelection({ clearStatus: true, syncDom: false });
    closeFeaturePopup();
  });
  watch(svgContent, () => {
    similarityAlignmentPorts.refreshCanvas();
    if (!skipCaptureBaseConfig.value) {
      clearFeatureSelection({ clearStatus: true, syncDom: false });
    }
  });
  watch(() => results.value[selectedResultIndex.value], () => similarityAlignmentPorts.refreshCanvas(), { flush: 'post' });

  /** @type {ReturnType<typeof setTimeout> | null} */
  let featureSearchDebounceId = null;
  const featureListScrollRef = ref(null);
  const selectedPairwiseBlockOrthogroupId = ref('');
  const resetFeatureListScroll = () => {
    featureListScrollTop.value = 0;
    if (featureListScrollRef.value) featureListScrollRef.value.scrollTop = 0;
  };
  /** @param {Event} event */
  const handleFeatureListScroll = (event) => {
    // A scroll event of the feature list element.
    const target = /** @type {HTMLElement | null} */ (event?.currentTarget || event?.target);
    featureListScrollTop.value = Number(target?.scrollTop || 0);
    const nextHeight = Number(target?.clientHeight || 0);
    if (nextHeight > 0) featureListViewportHeight.value = nextHeight;
  };
  watch(featureSearchInput, (/** @type {unknown} */ value) => {
    if (featureSearchDebounceId !== null) {
      clearTimeout(featureSearchDebounceId);
      featureSearchDebounceId = null;
    }
    const delay = String(value || '').trim() ? 120 : 0;
    featureSearchDebounceId = setTimeout(() => {
      featureSearch.value = String(value || '');
      resetFeatureListScroll();
      featureSearchDebounceId = null;
    }, delay);
  });
  watch(
    () => [
      selectedFeatureRecordIdx.value, selectedResultIndex.value, showRightDrawer.value,
      rightDrawerTab.value, extractedFeatures.value.length
    ],
    resetFeatureListScroll
  );

  /** @type {(() => void) | null} */
  let disposeHistoryInputs = null;
  setupGlobalUiEvents({
    state,
    onMounted,
    onUnmounted,
    closeRightDrawer: rightDrawerActions.closeRightDrawer,
    closeSettingsPopup: settingsPaneActions.closeSettingsPopup,
    previewPopupsHeld
  });
  setupHistoryShortcuts({ history, onMounted, onUnmounted });
  onMounted(async () => {
    disposeHistoryInputs = setupHistoryInputs({
      root: document.getElementById('app'),
      history,
      nextTick
    });
    await history.captureBaseline(msg('Initial state'));
  });
  onUnmounted(() => {
    if (featureSearchDebounceId !== null) clearTimeout(featureSearchDebounceId);
    if (typeof disposeHistoryInputs === 'function') disposeHistoryInputs();
    if (window.__GBDRAW_HISTORY__ === history) delete window.__GBDRAW_HISTORY__;
    previewFeatureSearch.dispose();
    featureActions.dispose();
    disposePanZoom();
    setUnmanagedConfigOverrideValidator(null);
    setStoredProteinResultFrameConverter(null);
    setLegacySessionUpgrader(null);
    setDiagramWorkerFontSource(null);
    setDiagramResourceRetainSource(null);
    disposeSessionOperations();
    disposeDiagramGenerationWorker();
  });

  const circularTrackNewRenderer = ref('dinucleotide_skew');
  const linearTrackNewRenderer = ref('dinucleotide_skew');
  const circularTrackSlotsPanelOpen = ref(false);
  const linearTrackSlotsPanelOpen = ref(false);
  const toggleCircularTrackSlotsPanel = () => {
    circularTrackSlotsPanelOpen.value = !circularTrackSlotsPanelOpen.value;
  };
  const toggleLinearTrackSlotsPanel = () => {
    linearTrackSlotsPanelOpen.value = !linearTrackSlotsPanelOpen.value;
  };
  const circularConservationFastaInput = ref(null);
  // The stack editors receive the feature placement transition as one port (R10, R13).
  const { changeTrackLayout } = featureActions.placementActions;
  // Legend styles and names follow the captions that track data names (OV-65,
  // OV-87): the annotation editor receives this one port, and the Depth source
  // and label transitions below run through it. A retired rename and the row of
  // a hidden track reach the displayed Result through the one Legend projection (R3).
  // Each writer passes the drawing it changes, so a replacement that reaches a
  // drawing not shown measures and retires that drawing's captions (review S2).
  /** @type {<T>(drawing: DrawingState, change: () => T) => T} */
  const retireLegendStylesOfUnnamedCaptions = createDrawingLegendRetirement({
    shownDrawingId: () => state.activeDrawingId.value,
    projectLegendEntries: () => { showEditorIntent({ domains: [...LIVE_EDIT_DOMAINS.legendRows, ...LIVE_EDIT_DOMAINS.strokes] }); },
    hiddenRows: (/** @type {DrawingState} */ drawing) => hiddenTrackLegendRows(drawing),
    relistLegendEntries: () => { legendActions.extractLegendEntries(); },
    namedCaptions: (/** @type {DrawingState} */ drawing) => trackDataLegendCaptions({
      annotationSets: drawing.annotationSets,
      depthTracks: drawing.adv.depth_tracks,
      ...(drawing.mode === 'linear'
        ? { depthSlots: drawing.adv.linear_track_slots, sourcedDepthTrackIndexes: activeDepthTrackIndices(linearDepthRows(drawing)) }
        : {
            depthSlots: drawing.adv.circular_track_slots,
            sourcedDepthTrackIndexes: circularDepthRepresentatives(drawing)
              .flatMap((file, index) => (file ? [index] : []))
          })
    })
  });
  // The stack and annotation editors edit the shown drawing.
  /** @type {<T>(change: () => T) => T} */
  const retireShownLegendStyles = (change) => retireLegendStylesOfUnnamedCaptions(state.activeDrawing(), change);
  // The diagrams the next Circular Generate draws, from the record set of the
  // request, read once: the Auto notes take the record lengths of each, and the
  // track slot estimates those of the diagram that draws the shown Result's
  // record. These estimates read the shown Circular drawing's settings; none
  // while another drawing is shown (its panel shows no Circular estimate).
  const circularDiagrams = computed(() => {
    const drawing = shownDrawingOfMode('circular');
    return drawing ? circularRequestDiagrams(drawing.circularRecordList.value, drawing) : [];
  });
  const circularEstimateRecordLengths = computed(() => circularEstimateLengths(
    circularDiagrams.value, getCommittedCanonicalRenderRequest(), state.selectedResultIndex.value
  ));
  const readCircularEstimateLengths = () => circularEstimateRecordLengths.value;
  const circularTrackSlotEditor = createCircularTrackSlotEditor({
    state, changeTrackLayout, askChoice, runUndoable: history.runUndoable, estimateRecordLengths: readCircularEstimateLengths,
    changeTrackVisibility: retireShownLegendStyles
  });
  const { chooseCircularRecordLayout } = createCircularMultiRecordCanvas({ state, rebindAnnotationRecords, projectCircularSkewHidden: circularTrackSlotEditor.projectCircularSkewHidden });
  const linearTrackSlotEditor = createLinearTrackSlotEditor({
    state, changeTrackLayout, changeTrackVisibility: retireShownLegendStyles
  });
  // PD-OI-104 scenario revisions 6 and 7: the shown stacks' asking ports (track-slot-edits.js).
  const { askingStackPorts, askingRowSwitch, askingRowRemoval } = createAskingStackPorts({
    askChoice, runUndoable: history.runUndoable, busy: sessionOperationAvailability, shownDrawing: () => state.activeDrawing()
  });
  // The Feature density section (app/feature-density-tracks.js); a track's stack rows change through the feature
  // placement transition (R10) and the mode's normalizer.
  const featureDensityTrackEditor = createFeatureDensityTrackEditor({
    shownDrawing: () => state.activeDrawing(), busy: sessionOperationAvailability, askingStackPorts, changeTrackLayout,
    normalizeStack: (drawing) => (drawing.mode === 'circular'
      ? circularTrackSlotEditor.normalizeCircularTrackSlots(drawing) : linearTrackSlotEditor.normalizeLinearTrackSlots())
  });
  const annotationImportNotice = ref('');
  const annotationEditor = createAnnotationEditor({
    state, getRecordCatalog: () => getAnnotationRecordCatalog(state.activeDrawing()), retireLegendStylesOfUnnamedCaptions: retireShownLegendStyles, askChoice,
    onImportNotice: (notice) => { annotationImportNotice.value = notice; }
  });
  // The Show Depth, and Linear Show GC Content and Show GC Skew, checkboxes:
  // the drawn tracks and their Legend rows leave and return together through
  // the Legend port (OV-1198, OV-1215, R10).
  /** @param {'show_depth' | 'show_gc' | 'show_skew'} key @param {unknown} checked */
  const setTrackShown = (key, checked) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const drawing = state.activeDrawing();
    const shown = Boolean(checked);
    retireLegendStylesOfUnnamedCaptions(drawing, () => {
      if (key === 'show_gc') drawing.form.show_gc = shown;
      else if (key === 'show_skew') drawing.form.show_skew = shown;
      else drawing.form.show_depth = shown;
    });
  };
  const {
    linearRegionNote, linearCardWindowView, linearLayoutCellRows, linearUnitPlacesOf, linearCardRecord, circularRegionNote, circularWindowRows, shownInputSummary,
    drawnLinearRecordLengths, committedLinearRecordsOf, awaitLinearRecordDiscovery, decorationsMoved, canvasPadded, shownCommittedRequest, shownRecordIds, shownLinearWrapLeaf, shownLinearWrapReading
  } = createShownDrawingReaders({
    state, shownDrawingOfMode, linearRecordsFor: linearRecordSelector.recordsFor,
    linearRecordStatusFor: linearRecordSelector.statusFor, refreshLinearRecords: linearRecordSelector.refresh,
    committedLinearRecordCards, getCommittedCanonicalRenderRequest, windowsOf: linearDrawingWindowsOf,
    unitPlacesOf: linearDrawingUnitPlaces, planCardWindows,
    cardWindowPlanRecords: linearCardWindowPlanRecords, recordDisplayRows: () => recordDisplayControls.allRows.value
  });
  const losatSettings = createLosatSettings({ state, planLinearUnits, shownPlacesOf: () => linearUnitPlacesOf(state.activeDrawing()) });
  // Auto values follow the records the next Generate draws; until they are
  // known the values name their rule.
  const autoValueDisplay = createAutoValueDisplay(state, {
    circularDiagramLengths: () => circularDiagrams.value.map((/** @type {CircularDiagramRecords} */ diagram) => diagram.lengths),
    linearRecordLengths: drawnLinearRecordLengths,
    linearWrappedAxisLength: () => {
      const reading = shownLinearWrapReading.value;
      if (reading.kind === 'off' || reading.kind === 'one-row') return undefined;
      return reading.kind === 'rows' ? reading.drawnAxisLength : null;
    }
  });
  const {
    hasCircularDepthFiles, canShowDepthTrack, circularDepthTrackRows, linearDepthTrackRows, linearSourceDepthRows,
    linearSourceDepthSummary, depthTrackRows, linearDepthTrackCoverageLabel, linearDepthTrackIndexOptions,
    getDepthTrackLabel, depthSeriesDefaultName, depthSeriesName, getDepthTrackColor, getDepthTrackLegendLabelForSlot,
    getCircularDepthFile
  } = createDepthTrackPanels({ state, shownDrawingOfMode });
  const enabledOptionClass = 'text-slate-700 cursor-pointer';
  const disabledOptionClass = 'text-slate-400 cursor-not-allowed opacity-60';
  const depthToggleOptionClass = computed(() => (
    canShowDepthTrack.value
      ? enabledOptionClass
      : disabledOptionClass
  ));
  /** @param {unknown} value */
  const optionalNumberInputValue = (value) => value ?? '';
  /** @param {Record<string, unknown> | null | undefined} target @param {string} key @param {unknown} value @param {boolean} [numeric] */
  const setOptionalNumberInputValue = (target, key, value, numeric = false) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!target || typeof target !== 'object') return;
    const text = String(value ?? '').trim();
    target[key] = text === '' ? null : numeric ? Number(text) : text;
  };
  // The problems a File card marks, for the open card and its folded line:
  // a file whose inspection failed, and a record's LOSAT Gencode under
  // TLOSATX (`losatGencodeProblem`): "Gencode required" when it is empty, the
  // tables LOSAT reads when the typed one is not one of them (OV-3326).
  /** @param {LinearSourceGroup} source a File card of the shown drawing */
  const linearSourceInspectionFailed = (source) => Boolean(linearRecordSelector.errorFor(state.activeDrawing(), source.sequence));
  /** @param {LinearSeq} seq */
  const linearRecordGencodeProblem = (seq) => {
    const drawing = shownDrawingOfMode('linear');
    return drawing?.hasActiveLinearLosatIntent.value && drawing.losatProgram.value === 'tblastx' ? losatGencodeProblem(seq) : null;
  };
  /** A File card's mark: a typed table LOSAT does not read before an empty one. @param {LinearSourceGroup} source */
  const linearSourceGencodeProblem = (source) => {
    const problems = source.records.map(({ sequence }) => linearRecordGencodeProblem(sequence));
    return problems.includes('unread') ? 'unread' : problems.includes('empty') ? 'empty' : null;
  };
  /** @param {'empty' | 'unread' | null} problem */
  const linearGencodeMark = (problem) => (problem === 'unread' ? errorText(losatGeneticCodeHint(LOSAT_GENETIC_CODE_TEXT))
    : problem === 'empty' ? t('Gencode required') : '');
  const {
    linearLabelHasSharedRow, linearLabelVisibilitySummary, linearLabelAutoFields, linearLabelAutoDisclosure,
    focusLinearLabelVisibility, getDefinitionLineStyleSize, setDefinitionLineStyleSize, getDefinitionLineStyleWeight,
    setDefinitionLineStyleWeight, getDefinitionLineStyleFill, setDefinitionLineStyleColor, getDefinitionLineStyleColorMode,
    getDefinitionLineStyleSwatchValue, isDefinitionLineStyleMuted
  } = createLinearRecordLabels({
    state, mode, sessionOperationAvailability, getCommittedLinearDefinitionVisibility,
    revealSetting: settingsPaneActions.revealSetting, optionalNumberInputValue, setOptionalNumberInputValue
  });
  /** @param {unknown} position */
  const legendPositionLabel = (position) => ({
    right: t('Right'),
    left: t('Left'),
    top: t('Top'),
    bottom: t('Bottom'),
    upper_left: t('Upper Left'),
    upper_right: t('Upper Right'),
    lower_left: t('Lower Left'),
    lower_right: t('Lower Right'),
    none: t('None')
  })[String(position || '').trim().toLowerCase()] || t('None');
  const {
    setDepthTrackColor, setDepthTrackLabel, setDepthTrackLegendLabelForSlot, syncDepthTrackSlotLabel
  } = createDepthSeriesEdits({
    state, mode, sessionOperationAvailability, retireLegendStylesOfUnnamedCaptions
  });
  const {
    replaceCircularDepthFile, setCircularDepthFile, setLinearDepthFiles, setLinearDepthFile, setLinearSourceDepthFile,
    clearLinearSourceDepthFile, addCircularDepthTrack, addLinearDepthTrack
  } = createDepthFiles({
    state, shownAction, retireLegendStylesOfUnnamedCaptions, askingStackPorts,
    changeCircularDepthSources: circularTrackSlotEditor.changeCircularDepthSources,
    circularDepthSourceLoss: circularTrackSlotEditor.circularDepthSourceLoss,
    changeLinearDepthSources: linearTrackSlotEditor.changeLinearDepthSources,
    linearDepthSourceLoss: linearTrackSlotEditor.linearDepthSourceLoss,
    runUndoable: history.runUndoable
  });
  // The Depth file rows, and the series that stay, after the removal of a
  // series: the removal writes them and its question reads them.
  /** @param {DrawingState} drawing @param {number} idx */
  const circularSeriesRemoval = (drawing, idx) => {
    const rows = removeDepthTrackColumnAt(circularDepthRows(drawing), idx).map((row) => compactDepthFileSlots(row));
    return { rows, keptSeries: representativeDepthFiles(rows).length };
  };
  /** @param {DrawingState} drawing */
  const keptLinearSeriesCount = (drawing) => Math.max(0, linearDepthLogicalWidth(drawing) - 1);
  /** @param {DrawingState} drawing @param {number} index @param {() => void} [turnOffReaching] */
  const removeCircularDepthTrack = (drawing, index, turnOffReaching) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0) return;
    // Removing a series retires the Legend styles of its rows (OV-65).
    retireLegendStylesOfUnnamedCaptions(drawing, () => {
      circularTrackSlotEditor.changeCircularDepthSources(drawing, () => {
        const previousFiles = circularDepthRepresentatives(drawing);
        const { rows, keptSeries } = circularSeriesRemoval(drawing, idx);
        bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_depth', rows);
        if (idx < drawing.adv.depth_tracks.length) drawing.adv.depth_tracks.splice(idx, 1);
        refreshDepthTrackLabelsAfterRemoval(drawing, previousFiles, circularDepthRepresentatives(drawing), idx);
        ensureDepthTrackConfigCount(drawing, activeDepthTrackCount(drawing));
        // The removed series' rows before the Axis lower its index, as in Linear,
        // so no other row crosses the Axis (R10).
        const axis = drawing.adv.circular_track_slots_axis_index;
        const removedBeforeAxis = Number.isInteger(axis)
          ? drawing.adv.circular_track_slots.slice(0, axis)
            .filter((/** @type {CircularTrackSlot} */ slot) => isDefaultManagedDepthSlot(slot) && depthSlotTrackIndex(slot) === idx).length
          : 0;
        drawing.adv.circular_track_slots.splice(
          0,
          drawing.adv.circular_track_slots.length,
          ...reindexDepthSlots(/** @type {ReindexDepthSlotsOptions} */ ({
            slots: drawing.adv.circular_track_slots,
            removedIndex: idx,
            activeCount: keptSeries
          }))
        );
        if (removedBeforeAxis) drawing.adv.circular_track_slots_axis_index = axis - removedBeforeAxis;
        turnOffReaching?.();
        syncDepthTrackSlotLabelsForTrack(drawing, idx);
        circularTrackSlotEditor.normalizeCircularTrackSlots(drawing);
      }, { setsFile: false, series: 'removed' });
    });
  };
  // A series removal asks first when it leaves track rows that Generate refuses.
  /** @param {number} index */
  const askingRemoveCircularDepthTrack = (index) => shownAction('circular', (drawing) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0) return removeCircularDepthTrack(drawing, index);
    const ports = askingStackPorts('circular_track_slots');
    return askTrackSourceRemoval(ports, {
      kind: 'depth-series',
      series: idx,
      loss: depthSeriesRemovalLoss(ports.rows(), idx, circularSeriesRemoval(drawing, idx).keptSeries),
      sources: () => circularDepthSources(drawing),
      label: msg('Remove Depth series'),
      apply: (turnOffReaching) => removeCircularDepthTrack(drawing, index, turnOffReaching)
    }) ?? removeCircularDepthTrack(drawing, index);
  });
  /** @param {DrawingState} drawing @param {number} index @param {() => void} [turnOffReaching] */
  const removeLinearDepthTrack = (drawing, index, turnOffReaching) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0) return;
    const logicalWidth = linearDepthLogicalWidth(drawing);
    if (idx >= logicalWidth) return;
    const keptSeries = keptLinearSeriesCount(drawing);
    // Removing a series retires the Legend styles of its rows (OV-65).
    retireLegendStylesOfUnnamedCaptions(drawing, () => {
      linearTrackSlotEditor.changeLinearDepthSources(drawing, () => {
        const nextRows = removeDepthTrackColumnAt(linearDepthRows(drawing), idx);
        drawing.linearSeqs.forEach((seq, recordIndex) => {
          seq.depth = nextRows[recordIndex] || [];
        });
        if (idx < drawing.adv.depth_tracks.length) drawing.adv.depth_tracks.splice(idx, 1);
        ensureDepthTrackConfigCount(drawing, activeDepthTrackCount(drawing));
        const previousAxisIndex = Number(drawing.adv.linear_track_slots_axis_index);
        const removedManagedSlotCountBeforeAxis = Number.isInteger(previousAxisIndex)
          ? drawing.adv.linear_track_slots.reduce((/** @type {number} */ count, /** @type {LinearTrackSlot} */ slot, /** @type {number} */ slotIndex) => {
              if (slotIndex >= previousAxisIndex || !isDefaultManagedDepthSlot(slot)) return count;
              return depthSlotTrackIndex(slot) === idx ? count + 1 : count;
            }, 0)
          : 0;
        drawing.adv.linear_track_slots.splice(
          0,
          drawing.adv.linear_track_slots.length,
          ...reindexDepthSlots(/** @type {ReindexDepthSlotsOptions} */ ({
            slots: drawing.adv.linear_track_slots,
            removedIndex: idx,
            activeCount: keptSeries
          }))
        );
        if (Number.isInteger(previousAxisIndex)) {
          drawing.adv.linear_track_slots_axis_index = Math.max(
            0,
            previousAxisIndex - removedManagedSlotCountBeforeAxis
          );
        }
        turnOffReaching?.();
        syncDepthTrackSlotLabelsForTrack(drawing, 0);
        linearTrackSlotEditor.syncLinearDepthSlotHeightsFromDepthTracks();
        linearTrackSlotEditor.normalizeLinearTrackSlots();
      }, { setsFile: false, series: 'removed' });
    });
  };
  /** @param {number} index */
  const askingRemoveLinearDepthTrack = (index) => shownAction('linear', (drawing) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= linearDepthLogicalWidth(drawing)) return removeLinearDepthTrack(drawing, index);
    const ports = askingStackPorts('linear_track_slots');
    return askTrackSourceRemoval(ports, {
      kind: 'depth-series',
      series: idx,
      loss: depthSeriesRemovalLoss(ports.rows(), idx, keptLinearSeriesCount(drawing)),
      sources: () => linearDepthSources(drawing),
      label: msg('Remove Depth series'),
      apply: (turnOffReaching) => removeLinearDepthTrack(drawing, index, turnOffReaching)
    }) ?? removeLinearDepthTrack(drawing, index);
  });
  const {
    setRingGenomeLists, setCircularConservationSourceFiles, replaceCircularConservationBlasts, adoptCircularRingResult,
    askRingSourceRemoval, reconcileCircularConservation, rederiveCircularConservation, reconcileLoadedCircularConservation,
    refreshCircularConservationSlots
  } = createCircularConservationSources({
    state, writeLosatCacheInfo, askingStackPorts,
    syncCircularConservationSlots: circularTrackSlotEditor.syncCircularConservationSlots,
    rekeyCircularConservationRows: circularTrackSlotEditor.rekeyCircularConservationRows,
    circularConservationRowsLeaving: circularTrackSlotEditor.circularConservationRowsLeaving
  });
  const circularConservationSeriesPanel = createCircularConservationSeriesPanel({
    state, shownAction, shownDrawingOfMode, writeLosatCacheInfo, circularConservationFastaInput,
    reconcileCircularConservation, refreshCircularConservationSlots, askRingSourceRemoval
  });
  const {
    settleComparisonRecordLabels, addCircularConservationComparisonFile, addRingRowFiles,
    removeRingRowFile, removeCircularConservationSource
  } = createCircularConservationComparisonFiles({
    state, shownAction, writeLosatCacheInfo, runUndoable: history.runUndoable, setRingGenomeLists,
    setCircularConservationSourceFiles, reconcileCircularConservation, refreshCircularConservationSlots, askRingSourceRemoval
  });
  // The auxiliary file imports still applying: the watchers add them, Save and
  // Load read them, and Generate waits for them (OV-1193).
  const auxiliaryFileImports = createAuxiliaryFileImports(reactive, state.inputFiles);
  // E1: the artifacts of the drawings not shown, in the slot store of the
  // History snapshot service; only its named operations change it.
  const { artifactSlots } = historySnapshots;
  // The diagram Worker keeps the resources of every drawing's committed
  // request, the shown drawing's and each kept slot's, and the files each
  // drawing's draft holds now while its request omits them.
  setDiagramResourceRetainSource({
    readCommittedResourceSets: () => [
      getCommittedCanonicalSession(),
      ...artifactSlots.kept().map((slot) => /** @type {ReturnType<typeof getCommittedCanonicalSession> | undefined} */ (
        slot.runtimeState?.canonical?.committedCanonicalSession))
    ].flatMap((committed) => (committed ? [{ request: committed.renderRequest, resources: committed.resources }] : [])),
    readDraftHeldFiles: () => state.drawingList().map((drawing) => ({
      drawingId: drawing.id, files: heldDepthFiles(state, drawing)
    }))
  });
  // The Results of all drawings: the shown one's and every kept one's. The
  // per-Result caches (Legend inventories, retired Legend entries, projected
  // editor state) keep a kept Result's entries while another drawing is shown.
  const liveResultIdentities = () => [
    ...results.value,
    ...artifactSlots.kept().flatMap((slot) => slot.values.results || [])
  ].map(previewRuntime.getResultIdentity);
  // The Results of one drawing: the shown one's, or those its slot keeps.
  /** @param {DrawingId} drawingId */
  const resultsOf = (drawingId) => artifactSlots.readOfDrawing(drawingId, state.activeDrawingId.value,
    () => results.value, (slot) => slot.values.results) || [];
  /** @param {DrawingId} drawingId */
  const resultIdentitiesOf = (drawingId) => resultsOf(drawingId).map(previewRuntime.getResultIdentity);
  // R13: a moved record clears the alignment plan and keeps the aligned record
  // positions. The alignment owner runs its candidate through Generate, which
  // reads this layout owner, so the root registers both ports once the
  // alignment owner exists.
  /** @type {{
   *   beforeRecordDrag: ReturnType<typeof createSimilarityAlignmentActions>['beforeRecordDrag'],
   *   afterRecordDrag: ReturnType<typeof createSimilarityAlignmentActions>['afterRecordDrag']
   * }} */
  // Null only until the root registers both ports once the alignment owner exists; no drag runs before then.
  const recordDragAlignmentPorts = {
    beforeRecordDrag: /** @type {never} */ (null),
    afterRecordDrag: /** @type {never} */ (null)
  };
  const legendLayout = createLegendLayout({
    state,
    layOutLegend: legendActions.layOutLegend,
    beginHistoryTransaction: history.begin,
    commitHistoryTransaction: history.commit,
    previewRuntime,
    similarityAlignmentLifecycle: {
      beforeRecordDrag: () => recordDragAlignmentPorts.beforeRecordDrag(),
      afterRecordDrag: (options) => recordDragAlignmentPorts.afterRecordDrag(options)
    }
  });
  legendActions.setLegendGeometryChangedHandler(legendLayout.refreshLegendGeometry);
  // History's capture registers once its owner exists (R13).
  // OV-110: each capture first releases the offsets of the Results that no
  // drawing shows or keeps and no History step names.
  historySnapshots.registerCapture('composition', (drawingId) => {
    legendLayout.releaseCompositionDeltas(new Set([...liveResultIdentities(), ...history.heldIdentities()]));
    return legendLayout.captureCompositionIntent(resultIdentitiesOf(drawingId));
  });
  // Each Result's bytes reflect the editor state it was committed or last
  // shown with (`createResultPaintRecord`); the binder below reads it (R13).
  const resultPaintRecord = createResultPaintRecord();
  // DESIGN A4: the record discovery owner's lists, by mode; a Linear card
  // whose list is still deferred after a Load is read here. Callers read them
  // through the runner's `readDrawingRecordIdentities`, which reads a deferred
  // Circular list first (OV-2489).
  /** @param {DrawingState} drawing @returns {Promise<RecordIdentitySources>} */
  const readRecordIdentitySources = async (drawing) => {
    if (drawing.mode === 'circular') return circularIdentitySources(drawing);
    await awaitLinearRecordDiscovery(drawing);
    return {
      mode: drawing.mode,
      cards: new Map(drawing.linearSeqs.map((seq) => [String(seq.uid), buildDisambiguatedRecordEntries(linearRecordSelector.recordsFor(drawing, seq))]))
    };
  };
  const {
    runInLosatOverlay, runLosatJobsInOverlay,
    runAnalysis: runGeneratedDiagramAnalysis,
    runCommittedCanonicalCandidate,
    projectCommittedRecordTransform,
    projectCommittedSimilarityAlignment,
    cancelRunAnalysis,
    runLabelReflow,
    holdAutomaticRerender,
    requestLiveSettingApply,
    discardLiveSettingRequests,
    refreshCircularRecordOrder: refreshRecordOrderOfDrawing,
    readDrawingRecordIdentities,
    releaseSourceInputFailure,
    downloadCliHelperFiles,
    losatResultRowsOf, clearableLosatResultCount, swapLosatResultRows, losatRowsOf,
    downloadLosatCache: downloadLosatCacheOfDrawing,
    downloadLosatPair: downloadLosatPairOfDrawing,
    clearLosatCache,
    getLosatPairDefaultName: losatPairDefaultNameOfDrawing,
    captureGeneratedArtifactRuntimeState,
    restoreGeneratedArtifactRuntimeState
  } = createRunAnalysis({
    state,
    rulePreparation,
    settlePendingInputs: async () => {
      await settleComparisonRecordLabels();
      await auxiliaryFileImports.settle();
    },
    isCurrentFeature: isCurrentResultFeature, writeSessionDocumentFile, serializeInputFile: serializeFile,
    serializeCanonicalFiles: (comparisonPlanSnapshot, linearRecordCatalog, drawing) => (
      serializeActiveRenderFiles(state, drawing, {
        comparisonPlan: comparisonPlanSnapshot,
        linearRecordCatalog,
        hiddenDepth: true
      })
    ),
    prepareLinearRecordCatalog, planGenomeScales: planDraftGenomeScales,
    readRecordIdentitySources,
    linearCardRecord,
    readLinearCardRecords: linearRecordSelector.refresh,
    placeCircularRecordChoices: (drawing) => placeCircularRecordChoices(drawing), // created below, with the record settings
    recordDisplayRows: recordDisplayControls.allRows,
    assertDrawingInputs,
    readDraftSignature: () => history.getCurrentIntentSignature(),
    readHistoryStepSerial: () => history.intentStepSerial(),
    adoptDrawnWrite: (write, request) => history.adoptDrawnWrite(write, request),
    closeLabelTextScopeDialog: featureActions.closeLabelTextScopeDialog,
    clearLabelBuildNotices: featureActions.clearLabelBuildNotices,
    moveDrawnLabelVisibility: featureActions.moveDrawnLabelVisibility,
    askChoice,
    showNotice,
    // Only a Generate that searched its drawing's rings installs them, and the
    // rings' reconcile pairs them with the rings it searched
    // (`adoptCircularRingResult`); a Linear Generate, a Generate that searched no
    // ring and a candidate of the committed request leave them.
    adoptCanonicalRenderArtifacts: (canonical, { ringAdoption = null, ...options }) => adoptCircularRingResult(
      options.drawing, (installRings) => adoptCanonicalRenderArtifacts(canonical, { ...options, installRings }), ringAdoption
    ),
    getCommittedCanonicalSession,
    readCommittedCanonicalRevision,
    commitLiveSettingValues,
    readLiveRedrawLegend: legendActions.liveRedrawLegend,
    captureDecorationContinuity: legendLayout.captureDecorationContinuity,
    captureGeneratedArtifactHandle: historySnapshots.captureGeneratedArtifactHandle,
    displayedArtifactDrawing: historySnapshots.displayedArtifactDrawing,
    captureGeneratedArtifactOwnerSet: historySnapshots.captureGeneratedArtifactOwnerSet,
    installGeneratedArtifactOwnerSet: historySnapshots.installGeneratedArtifactOwnerSet,
    installRedrawnArtifact: historySnapshots.installRedrawnArtifact,
    restoreGeneratedArtifactHandle: (...args) => restoreWithSpecificRuleDrafts(historySnapshots.restoreGeneratedArtifactHandle, ...args),
    emptyKeptProteinIdentityManifests: historySnapshots.emptyKeptProteinIdentityManifests,
    setGeneratedArtifactIdentity: historySnapshots.setGeneratedArtifactIdentity,
    runGeneratedArtifactReplacement: (...args) => history.runUndoableArtifactReplacement(...args),
    previewRuntime,
    nextTick,
    onGeneratedArtifactCheckpointCapture: ({ phase, diagnostics } = {}) => {
      recordSessionLifecycleEvent(`generate-history-${String(phase || '')}`, {
        diagnostics: diagnostics && typeof diagnostics === 'object'
          ? diagnostics
          : {}
      });
    },
    resetPreviewViewport,
    validateAnnotationTargets: () => {
      const drawing = state.activeDrawing();
      const catalog = getAnnotationRecordCatalog(drawing);
      reconcileAnnotationRecordBindings(drawing, catalog);
      return validateAnnotationRecordTargets(drawing.annotationSets, catalog);
    },
    omittedAnnotationRecordKeys: () => getAnnotationRecordCatalog(state.activeDrawing()).omittedRecordKeys || [],
    searchRecordSelector: (bindingKey) => shownResultRecordBindings().selectorOfBinding(bindingKey),
    onCircularRecordsDiscovered: handleCircularRecordsDiscovered
  });
  // OV-2678: the rerender a checkpoint step asks for, its own and its listener's, starts once the step ends.
  history.setCheckpointHold(holdAutomaticRerender);
  // D-15-89: the panel actions (the template, `window.__GBDRAW_APP__`, the record editor, and the file watcher) act
  // for the shown drawing; this binding resolves it once per action and the runner reads the one it is given.
  /** @param {{ suppress?: boolean, automatic?: boolean } | undefined} [options] */
  const refreshCircularRecordOrder = (options) => refreshRecordOrderOfDrawing(options, state.activeDrawing());
  // R10: the record discovery of one drawing. Each transition that changes the drawing's discovery inputs calls it
  // for that drawing: a source file, the input type, the Linear cards, the end of a Load or a Save, and a History
  // restore of inputs. A drawing switch calls nothing. It is the one admission of a discovery: none starts while a
  // Load rolls back or a Load or History restore installs (the semantic file watchers are suppressed); a Save or Load
  // still pending makes the owners return busy, and its end calls the discovery again. It writes only `drawing`.
  /** @param {DrawingState} drawing */
  const discoverDrawingRecords = (drawing) => {
    if (state.sessionImportRollbackInProgress.value || state.semanticFileWatchersSuppressed.value) return undefined;
    if (drawing.mode === 'circular') return refreshRecordOrderOfDrawing({ automatic: true }, drawing);
    // UJ-08: a discovery that reads every Linear card ends a Generate failure
    // about those inputs, as the Circular discovery does.
    return Promise.resolve(linearRecordSelector.refresh(drawing)).then((outcome) => {
      if (drawing.linearSeqs.every((seq) => ['ready', 'error'].includes(linearRecordSelector.statusFor(drawing, seq)))) {
        releaseSourceInputFailure(drawing);
      }
      return outcome;
    });
  };
  /** @param {string} edgeKey @param {unknown} customName */
  const downloadLosatPair = (edgeKey, customName) => downloadLosatPairOfDrawing(state.activeDrawing(), edgeKey, customName);
  const losatResultRows = computed(() => losatResultRowsOf(state.activeDrawing()));
  const downloadLosatCache = () => downloadLosatCacheOfDrawing(state.activeDrawing());
  /** @param {string} edgeKey */
  const getLosatPairDefaultName = (edgeKey) => losatPairDefaultNameOfDrawing(state.activeDrawing(), edgeKey);
  // R1 (c), PD-OI-102: a History step, a checkpoint (Reset Settings, a rule
  // preset), and their Undo and Redo that change live setting rows of the
  // shown drawing ask the automatic rerender to apply them (the router passes
  // the shown drawing's steps only).
  // Figure size (PD-OI-114 rev 2, Q-D1 A): a step that changes what Python draws of the shown drawing's Legend rows or
  // its Canvas padding (its commit, Undo, or Redo) shows at once, and a figure drawn to a requested size is then fitted
  // again by the automatic rerender, which draws them (`projectCommittedEditorIntent`). An Auto figure is not redrawn.
  // History marks a checkpoint side undrawn while that fit is owed and keeps the rerender's figure as the last step's;
  // a restore of an undrawn side asks the label owner to draw it again (`redrawRestored`, D-24).
  /**
   * @param {readonly HistoryChangeLike[] | null} changes
   * @param {'undo' | 'redo' | undefined} direction
   * @param {boolean | undefined} undrawn
   */
  const refitSizedFigure = (changes, direction, undrawn) => {
    const root = svgContainer.value?.querySelector?.('svg');
    if (!root || !changes || !requestedCanvasOf(root)) return;
    const drawing = state.activeDrawing();
    // A commit whose caller asks for the automatic rerender itself (a color rule that adds a Legend row) is fitted
    // by that rerender; a second request would draw the figure twice.
    if (undrawn && !direction) return;
    if (!stepRefitsDrawnLegend(changes, direction, {
      legendEntries: drawing.legendEntries.value, deletedLegendEntries: drawing.deletedLegendEntries.value,
      rules: drawing.manualSpecificRules, canvasPadding: drawing.canvasPadding, originalLegendOrder: originalLegendOrder.value,
      // The rows the request names (`requestLegendEditLeaves`): those Python drew in every Result of the drawing.
      pythonRowKeys: pythonLegendRowKeys(originalLegendOrder.value, results.value)
    })) return;
    if (featureActions.requestAutomaticRerender() && !direction) history.markLastStepUndrawn();
  };
  /** @param {{ recordKey: string, biologicalFeatureId: string }} identity */
  const resolvePopupRotationFeature = ({ recordKey, biologicalFeatureId }) => {
    /** @param {{ record_key?: unknown, recordKey?: unknown, biological_feature_id?: unknown, biologicalFeatureId?: unknown } | null | undefined} feature */
    const matchesIdentity = (feature) => (
      String(feature?.record_key ?? feature?.recordKey ?? '') === recordKey
      && String(
        feature?.biological_feature_id ?? feature?.biologicalFeatureId ?? ''
      ) === biologicalFeatureId
    );
    const renderedMatches = extractedFeatures.value.filter(matchesIdentity);
    if (renderedMatches.length === 1) return renderedMatches[0];
    const biologicalMatches = biologicalFeatures.value.filter(matchesIdentity);
    return biologicalMatches.length === 1 ? biologicalMatches[0] : null;
  };
  const featureRecordRotationAction = createFeatureRecordRotationAction({
    targetForFeature: recordDisplayControls.targetForFeature,
    setResolvedTransform: recordDisplayControls.setResolvedTransform,
    getCommittedSession: getCommittedCanonicalSession,
    projectCommittedRecordTransform,
    // R13: the root runs the candidate with the record display's target-draft
    // checkpoint; the rotation owner holds neither owner.
    runRecordRotation: ({ canonical, row, transform }) => runCommittedCanonicalCandidate({
      canonical,
      captureIntentCheckpoint: () => recordDisplayControls.captureTargetDraft(row),
      restoreIntentCheckpoint: (/** @type {TargetDraftCheckpoint} */ checkpoint) => recordDisplayControls.restoreTargetDraft(checkpoint),
      commitIntent: () => recordDisplayControls.commitResolvedTransform(row, transform)
    }),
    resolveCurrentFeature: resolvePopupRotationFeature,
    isCurrentFeature: isCurrentResultFeature,
    // The same discovery that Generate and the File card run for the mode.
    readRecords: () => (mode.value === 'linear'
      ? linearRecordSelector.refresh(state.activeDrawing())
      : refreshCircularRecordOrder())
  });
  const featureRecordRotation = createFeatureRecordRotationWorkflow({
    action: featureRecordRotationAction,
    makeReactive: reactive,
    onRebind: (feature, identity) => {
      const popup = clickedFeature.value;
      if (!popup) return;
      const popupFeature = popup.feat;
      const popupRecordKey = String(
        popupFeature.record_key ?? popupFeature.recordKey ?? ''
      );
      const popupFeatureId = String(
        popupFeature.biological_feature_id ?? popupFeature.biologicalFeatureId ?? ''
      );
      if (popupRecordKey !== identity.recordKey
        || popupFeatureId !== identity.biologicalFeatureId) return;
      popup.feat = feature;
      const renderedId = String(
        feature?.rendered_feature_svg_id
        ?? feature?.renderedFeatureSvgId
        ?? feature?.svg_id
        ?? ''
      );
      if (renderedId) popup.svg_id = renderedId;
    }
  });
  const recordActionsExpanded = ref(false);
  watch(clickedFeature, (/** @type {FeaturePopup | null} */ popup) => {
    recordActionsExpanded.value = false;
    if (popup) {
      featureRecordRotation.open({ feature: popup.feat });
    } else {
      featureRecordRotation.close();
    }
  }, { flush: 'sync' });
  const closeFeaturePopup = () => {
    featureRecordRotation.cancel();
    recordActionsExpanded.value = false;
    clickedFeature.value = null;
  };
  const toggleRecordActions = () => {
    if (!recordActionsExpanded.value && !featureRecordRotation.draft.active && clickedFeature.value?.feat) {
      featureRecordRotation.open({ feature: clickedFeature.value.feat });
    }
    recordActionsExpanded.value = !recordActionsExpanded.value;
    if (recordActionsExpanded.value) featureRecordRotation.readRecords();
  };
  // A record that rotation cannot use (a linear or cropped record)
  // shows the rotation row disabled with this reason; a record a Session Load
  // left unread is read when the row opens.
  const featureRecordRotationUnavailable = computed(() => {
    const { targetFailure, disabledReason } = featureRecordRotation.draft;
    return targetFailure && targetFailure !== RECORD_TARGET_NOT_DISCOVERED && !recordActionsExpanded.value
      ? disabledReason : '';
  });
  // Cancel discards the draft and starts a fresh one for the open popup, as
  // opening the popup does, so a record found unusable stays explained.
  const cancelRecordActions = () => {
    featureRecordRotation.cancel();
    recordActionsExpanded.value = false;
    if (clickedFeature.value?.feat) featureRecordRotation.open({ feature: clickedFeature.value.feat });
  };
  historySnapshots.setGeneratedArtifactRuntimeOwner({
    capture: () => ({
      ...captureGeneratedArtifactRuntimeState(),
      canonical: canonicalRenderArtifactOwner.capture()
    }),
    // `null` is an empty mode slot (E1): no committed Session, no CLI helper files.
    restore: (snapshot, options) => {
      canonicalRenderArtifactOwner.restore(snapshot?.canonical ?? null);
      return restoreGeneratedArtifactRuntimeState(snapshot ?? {}, options);
    }
  });
  const {
    importAuxiliaryTable,
    canRetryAuxiliaryImportFailure, retryAuxiliaryImportFailure, auxiliaryImportFailureOnArrival, resetModeTransientUi
  } = setupWatchers({
    state,
    isRuleInputCurrent: rulePreparation.isCurrent, snapshotRuleInputs: rulePreparation.snapshot,
    trackAuxiliaryFileImport: auxiliaryFileImports.track,
    ref, computed, watch,
    nextTick,
    onMounted,
    legendActions,
    featureActions,
    legendLayout,
    resultsManager,
    runLabelReflow,
    resetPreviewViewport,
    resetRightDrawer: rightDrawerActions.resetRightDrawer,
    closeLabelTextScopeDialog: featureActions.closeLabelTextScopeDialog,
    closeHiddenLabelTextDialog: featureActions.closeHiddenLabelTextDialog,
    clearLabelBuildNotices: featureActions.clearLabelBuildNotices,
    previewRuntime,
    // The shown drawing may take its palette's colors: its palette draft syncs.
    preparePaletteDefinitions: async () => {
      await paletteLoader.loadPaletteAsset();
      resultsManager.syncPaletteDraftState();
    }
  });

  // Save and Load also wait for an edit still applying. Its owners exist only
  // here, so the root composes the availability that Save, Load, and their
  // controls read (R13).
  const sessionPreparationBusyReason = () => (
    history.mutationPending() || ruleMatchingBlocking.value || auxiliaryFileImports.isPending()
      ? msg('Applying an edit. Retry after the edit finishes.')
      : ''
  );
  /** @param {string} operation */
  const sessionSaveLoadAvailability = (operation) => (
    sessionOperationAvailability(operation, sessionPreparationBusyReason)
  );
  const semanticMutationAvailable = computed(() => !sessionOperationAvailability());
  const dialogChoicePending = dialogChoice.pending;
  const sessionSaveAvailable = computed(() => !sessionSaveLoadAvailability('save'));
  const sessionLoadAvailable = computed(() => !sessionSaveLoadAvailability('load'));
  // The reason is English (its owners mark it with `msg`); it is translated where it is shown.
  const sessionBusyReason = computed(() => tm(sessionSaveLoadAvailability('save')?.reason || ''));
  // OV-378: the first feature popup after a Session load waits for Python to
  // match the rules; say so, since Save and Load stay available meanwhile.
  const featureDetailsPending = computed(() => ruleMatchingPending.value && !ruleMatchingBlocking.value);
  const circularRecordPresentationPanel = ref(null);
  // UJ-09 (Owner 2026-10-05): loading a Session replaces the work and clears
  // History, so every route that loads one asks first when History changed
  // since the last Save or Load. Cancel loads nothing.
  /** @param {string} actionLabel @param {() => void} proceed */
  const confirmSessionReplacement = async (actionLabel, proceed) => {
    if (await history.hasChangesSinceSavePoint() && await askChoice({
      title: t('Replace the current work?'),
      description: joinSentences(t('Loading a Session replaces the current settings and Results and clears Undo history.'),
        t('Changes since the last Save Session or Load Session are lost.')),
      actions: [{ id: 'load', label: actionLabel, tone: 'primary' }]
    }) !== 'load') return;
    // `proceed` may open the file picker: the user activation of the click that
    // asked, or of the click that confirmed, still holds after these awaits.
    proceed();
  };
  // The template's hidden file input (`ref="sessionInput"`).
  const sessionInput = ref(null);
  const openSessionFilePicker = () => confirmSessionReplacement(
    t('Load Session'), () => sessionInput.value?.click()
  );
  let nextSessionPreviewToken = 1;
  const circularDrawings = () => state.drawingList().filter((drawing) => drawing.mode === 'circular');
  // The Load rollback port: the one the Load options carry and `window.__GBDRAW_APP__` hands to page callers.
  const sessionImportRollbackState = createSessionImportRollbackState({
    artifactSlots,
    captureHistoryMembers: () => history.captureMembers(),
    replaceHistoryMembers: (drawingIds) => history.replaceMembers(drawingIds),
    restoreHistoryMembers: (members) => history.restoreMembers(members),
    captureDisplayedArtifact: () => historySnapshots.captureArtifactSlot(),
    // The rollback restored the pan already; the install keeps it.
    installDisplayedArtifact: (slot) => historySnapshots.installArtifactSlot(slot, {
      drawing: state.drawingOf(slot.drawingId), ui: { canvasPan: { x: canvasPan.x, y: canvasPan.y } }
    }),
    featureListScrollTop,
    featureListScrollRef,
    selectedPairwiseBlockOrthogroupId,
    captureSpecificRulePatternDrafts: featureActions.captureSpecificRulePatternDrafts,
    restoreSpecificRulePatternDrafts: featureActions.restoreSpecificRulePatternDrafts
  });
  // A loaded Session says so in the notice banner, with what the load changed
  // without a loss as a list; what it could not restore was listed in the
  // Load dialog before. A Gallery example's notice also offers its tutorial
  // (`noticeLink`, #1383 PD-OI-108).
  /** @param {SessionFileInput} event @param {{ noticeLink?: NoticeLink | null }} [options] */
  const importSession = async (event, { noticeLink = null } = {}) => {
    const name = String(event?.target?.files?.[0]?.name || '');
    let fitLoaded = false;
    try {
      const result = await loadSessionFromFile(event);
      // R10: a Load that installed a Session, or rolled back to the previous
      // one, ends with each drawing's discovery, unless the Session defers it.
      const recordDiscoveryDeferred = Boolean(result && 'recordDiscoveryDeferred' in result && result.recordDiscoveryDeferred);
      if (['ok', 'legacy', 'error', 'stale'].includes(String(result?.status)) && !recordDiscoveryDeferred) {
        state.drawingList().forEach((drawing) => { void discoverDrawingRecords(drawing); });
      }
      // The notice of an earlier operation does not stand for this load; a
      // Load the user canceled changed nothing, the notice included.
      if (name && !['canceled', 'skipped', 'stale'].includes(String(result?.status))) noticeBanner.dismiss();
      if (result?.status === 'ok') {
        showNotice({ text: name ? t('Session loaded: {file}.', { file: name }) : t('Session loaded.'), items: result.notices || [] });
      } else if (result?.status === 'legacy') {
        showNotice({ text: t('Configuration loaded from an older gbdraw. Save Session to keep it in the current format.') });
      }
      if (noticeLink && (result?.status === 'ok' || result?.status === 'legacy')) noticeBanner.addLink(noticeLink);
      // The phone layout shows a loaded Result in the Preview, as Generate does; a Load without one, Settings.
      if (result?.status === 'ok' || result?.status === 'legacy') settingsPaneActions.showPhonePane(svgContent.value ? 'preview' : 'settings');
      // R17-25: a file Session Load with a Result fits the figure to the shown Preview at every width, after the
      // phone pane switch and the whole notice (its height sizes the canvas); the figure it mounted shows only then,
      // fitted, and also if this end throws (OV-2702). Generate resets to 100 %; a legacy configuration Load keeps the zoom.
      fitLoaded = result?.status === 'ok' && Boolean(svgContent.value);
      if (fitLoaded) await nextTick();
      // A loaded Session makes a Generate likely: download the Worker's boot
      // files now, after the Load, so they never slow its own requests.
      if (result?.status === 'ok' || result?.status === 'legacy') prefetchPyodideBootAssets();
      return result;
    } finally { showHiddenPreview({ fit: fitLoaded }); }
  };
  // The Sessions this page has installed: a Load replaces the document, so
  // an owner that holds what an action did to the previous one (Align at
  // feature's applied state) reads this count instead of the drawn values.
  const sessionReplacements = ref(0);
  /** @param {SessionFileInput} event */
  const loadSessionFromFile = (event) => importSessionFromFile(event, {
    shownDrawing: state.activeDrawing,
    availability: sessionSaveLoadAvailability,
    // Before the Session replaces the current one: a font it names that this
    // page cannot draw asks first (receipt 〔読み込むとき〕); Cancel loads nothing.
    confirmSessionFonts: fontSettings.confirmSessionFonts,
    // BA-D25: one question before a Load that cannot restore everything.
    confirmLoadRepairs: async (dialog) => {
      const load = await askChoice({ ...dialog, actions: [{ id: 'load', label: t('Load'), tone: 'primary' }] }) === 'load';
      if (load) noticeBanner.dismiss();
      return load;
    },
    transformLegacyResultSvg,
    beforeImport: async () => {
      sequenceSearch.loadStarts();
      await nextTick();
      await afterPaint();
      recordSessionLifecycleEvent('session-import-paint-opportunity-completed');
    },
    installLoadedArtifactSlots: (/** @type {LoadedArtifactSlots} */ slots) => installLoadedArtifactSlots(slots),
    beforePreviewMount: (/** @type {{ results: any[], resultIndex: number }} */ {
      results: importedResults, resultIndex
    }) => {
      hidePreviewUntilShown();
      const selectedResult = importedResults[resultIndex] || null;
      if (!selectedResult) return null;
      const token = `session-load:${nextSessionPreviewToken++}`;
      return previewRuntime.registerReadinessExpectation(/** @type {any} */ ({
        result: selectedResult,
        resultIndex,
        artifactIdentity: token,
        generationToken: token,
        catalogState: state.featureCatalog.value,
        phase: 'session-load',
        bindingOptions: { isIncrementalEdit: true },
        isCurrent: () => (
          results.value[resultIndex] === selectedResult
          && Number(selectedResultIndex.value) === resultIndex
        )
      }));
    },
    rollbackState: sessionImportRollbackState, storedFigurePreview,
    afterImport: async (/** @type {{ status: string, ringSettingStored: boolean, decompressedCharacters?: number, isCurrent: () => boolean }} */ result) => {
      // Only a Session Load comes here; a legacy configuration Load returns
      // before it (OV-1792).
      if (result?.status === 'ok') {
        // The comparison rings read from the loaded files and series (an older
        // Session's series names and rings setting included), before the History baseline.
        for (const drawing of circularDrawings()) reconcileLoadedCircularConservation(drawing, result.ringSettingStored);
        // The loaded draft's live values wait for their next change or Generate (6.7).
        discardLiveSettingRequests();
        historySnapshots.clearGeneratedArtifactIdentity({
          retainedBytes: Number(result.decompressedCharacters || 0) * 2
        });
        // A replaced document resets transient UI: the selection named features
        // of the previous Session.
        clearFeatureSelection({ clearStatus: true, syncDom: false });
        previewFeatureSearch.resetModeSearches();
        sequenceSearch.reset();
        settingsPaneActions.closeSettingsPopup();
        await nextTick();
        // The Load's commit gave each loaded drawing an empty History
        // (`replaceHistoryMembers`); a legacy configuration Load does not
        // come here: its drawing keeps its History (OV-1790).
        recordSessionLifecycleEvent('history-baseline-start');
        if (!await history.initializeIntentBaseline(msg('Loaded session'), { isCurrent: result.isCurrent })) {
          throw new Error('Session loading canceled.');
        }
        featureActions.clearSpecificRulePatternDrafts();
        annotationImportNotice.value = '';
        recordSessionLifecycleEvent('history-baseline-end');
        if (circularRecordPresentationPanel.value) circularRecordPresentationPanel.value.open = false;
        closeLegendStrokeOptions();
        sessionReplacements.value += 1;
      }
    }
  });
  // UJ-06 (Owner 2026-10-08): Open Gallery example shows only where this origin
  // serves the Gallery catalog (gbdraw.app, not the pip GUI), and its chooser
  // lists every catalog example. The choice closes the chooser, asks the UJ-09
  // question, and loads the Session through the Load Session path. The catalog
  // read does not hold up startup. A loaded Session, an example
  // included, opens fitted (importSession).
  const galleryExamples = ref(/** @type {GalleryExample[]} */ ([]));
  void readGalleryExamples().then((examples) => { galleryExamples.value = examples; });
  const exampleChooserOpen = ref(false);
  const exampleSessionPending = ref(false);
  const openExampleChooser = () => { exampleChooserOpen.value = true; };
  const closeExampleChooser = () => { exampleChooserOpen.value = false; };
  /** @param {GalleryExample} example */
  const loadExampleSession = async (example) => {
    closeExampleChooser();
    // The chooser returns focus to its button before the question opens.
    await nextTick();
    await confirmSessionReplacement(t('Load example'), async () => {
      exampleSessionPending.value = true;
      try {
        const file = await fetchGalleryExampleSession(example);
        // The import resolves after the shown Result's preview is mounted. #1383 (PD-OI-108):
        // the notice of the loaded example offers its tutorial, found by the example's Gallery id.
        await importSession({ target: { files: [file], value: '' } }, {
          noticeLink: example.tutorialUrl ? {
            text: t('This diagram has a step-by-step tutorial.'),
            label: tx('English-only Gallery tutorial', 'Open the tutorial'),
            href: example.tutorialUrl
          } : null
        });
      } catch (error) {
        errorLog.value = normalizeUserFacingError(error);
      } finally {
        exampleSessionPending.value = false;
      }
    });
  };

  const {
    legendRowHasRules,
    updateLegendEntryColor,
    deleteLegendEntry,
    moveLegendEntry,
    sortLegendEntries,
    sortLegendEntriesByDefault,
    resetLegendPosition,
    getLegendEntryStrokeColor,
    getLegendEntryStrokeWidth,
    isLegendStrokeOptionsOpen,
    toggleLegendStrokeOptions,
    closeLegendStrokeOptions,
    setLegendEntryStrokeColorValue,
    updateLegendEntryStrokeColor,
    updateLegendEntryStrokeWidth,
    resetLegendEntryStroke,
    resetAllStrokes,
    restoreDeletedLegendEntries
  } = legendActions;

  const {
    addCustomColor,
    addPriorityRule,
    setLabelFilterMode,
    setBlacklistText,
    addWhitelistRule,
    setWhitelistRuleField,
    removeWhitelistRule,
    removePriorityRule,
    addFeature,
    removeFeature,
    getFeatureShape,
    setFeatureShape,
    addSpecificRule,
    applySpecificRulePreset,
    clearAllSpecificRules,
    downloadSpecificRulesTsv,
    moveSpecificRuleDown,
    moveSpecificRuleUp,
    removeSpecificRule,
    setSpecificRuleField,
    specificRulePattern, specificRulePatternDraft, specificRulePatternFieldId,
    editSpecificRulePattern, retrySpecificRulePattern, revertSpecificRulePattern,
    addFeatureVisibilityRule,
    downloadFeatureVisibilityRulesTsv,
    featureVisibilityQualifierSuggestions,
    featureVisibilityRuleDetail,
    getFeatureColor,
    getFeatureColorValue,
    handleFeatureVisibilityScopeChoice,
    moveFeatureVisibilityRuleDown,
    moveFeatureVisibilityRuleUp,
    removeFeatureVisibilityRule,
    setFeatureVisibility,
    setFeatureVisibilityRuleField,
    updateClickedFeatureVisibility,
    requestFeatureColorChange,
    setFeatureColorValue,
    updateClickedFeatureColor,
    cancelFeatureStyleScope,
    cancelLegendRename,
    cancelResetColor,
    handleColorScopeChoice,
    handleFeatureStyleScopeChoice,
    handleLegendNameCommit,
    handleLegendRenameChoice,
    joinLegendGroup,
    renameLegendEntry,
    handleResetColorChoice,
    resetClickedFeatureFillColor,
    getFeatureStrokeColorValue,
    setClickedFeatureStrokeColorValue,
    setClickedFeatureStrokeWidthValue,
    updateClickedFeatureStroke,
    resetClickedFeatureStroke,
    applyColorToSelectedFeatures,
    applyStrokeToSelectedFeatures,
    buildSelectedFeaturesVisibilityCommand,
    setFeatureColor,
    openFeatureEditorForFeature,
    getEditableLabelByFeatureId,
    syncLabelEditor,
    downloadLabelOverrideTable,
    loadLabelOverrideTable,
    updateClickedFeatureLabelText,
    handleLabelTextScopeChoice,
    handleHiddenLabelTextChoice,
    handleLabelOnChoice,
    clickedFeatureLabelHint,
    hiddenLabelTextMessage,
    requestLabelTextChangeByFeatureId,
    requestLabelTextChangeByKey,
    projectFeatureVisibility,
    reconcileLabelOverrides,
    resetAllLabelTextOverrides
  } = featureActions;

  const {
    compileDisplayedResultOperations, showEditorIntent, STROKE_DOMAINS, LEGEND_FILL_DOMAINS, VISIBILITY_DOMAINS,
    editEditorIntent, showingEditorIntent, LEGEND_STRUCTURE_DOMAINS
  } = createEditorIntentShow({
    state,
    results,
    extractedFeatures,
    originalLegendOrder,
    svgContainer,
    selectedResultIndex,
    readDrawnTracks,
    getCommittedCanonicalRenderRequest,
    getActiveRuntime: previewRuntime.getActiveRuntime,
    applyEditorOperations: previewRuntime.applyEditorOperations,
    readLegendRowNames: legendActions.readLegendRowNames,
    onLegendGeometryChanged: legendActions.onLegendGeometryChanged,
    runUndoable: history.runUndoable
  });
  const {
    projectMountedEditorIntent, drawReturningLegendRows, restoreLegendRows, restoreLegendItems
  } = createMountedEditorIntentProjection({
    state,
    svgContainer,
    projectPaletteAndRules,
    projectFeatureVisibility,
    reconcileLabelOverrides,
    fillDomainsOf,
    VISIBILITY_DOMAINS,
    showEditorIntent,
    requestAutomaticRerender: featureActions.requestAutomaticRerender,
    runUndoable: history.runUndoable
  });

  historySnapshots.setAfterApplyHistoryIntent(async (_intent, /** @type {{ domains?: Set<string>, changes?: Record<string, any>, direction?: string }} */ { domains, changes, direction } = {}) => {
    if (!svgContainer.value?.querySelector?.('svg')) return;
    const changedDomains = domains instanceof Set ? domains : new Set();
    // B17: restore the offsets of each Result whose composition this step changed.
    const compositionResults = new Set();
    (Array.isArray(changes) ? changes : []).forEach(({ path, before, after } = {}) => {
      if (path?.[0] !== 'ui' || path[1] !== 'compositionUserDeltas') return;
      if (path.length > 2) compositionResults.add(path[2]);
      else [before, after].forEach((record) => Object.keys(record || {}).forEach((key) => compositionResults.add(key)));
    });
    if (compositionResults.size > 0) {
      legendLayout.reconcileCompositionUserDeltas(_intent?.ui?.compositionUserDeltas, compositionResults);
    }
    const colors = changedDomains.has('config') || changedDomains.has('features');
    const rulesChanged = Array.isArray(changes) && changes.some((change) => (
      change?.path?.[0] === 'config' && change.path[1] === 'rules'
    ));
    const editorState = changedDomains.has('editorState');
    // The Legend side of this step: `from`, the list the step leaves. A step
    // made on another batch Result restores that Result's list, so the Legend
    // owner writes the displayed Result's list before it shows (B19).
    const legendChange = (Array.isArray(changes) ? changes : []).find(({ path } = {}) => (
      path?.length === 3 && path[0] === 'editorState' && path[1] === 'legend' && path[2] === 'entries'
    ));
    if (legendChange) legendActions.adoptRestoredLegend({ from: legendChange[direction === 'undo' ? 'after' : 'before'] });
    const projected = await projectMountedEditorIntent({
      colors,
      prepareRules: changedDomains.has('features') || rulesChanged || !rulePreparation.isPrepared(),
      visibility: changedDomains.has('features'),
      rerender: true,
      labels: changedDomains.has('features') || editorState,
      domains: editorPaintDomains(changes)
    });
    if (!projected) return;
    await nextTick();
  });

  const {
    currentEditorProjectionState, departedResultIntent, displayedIntentSignature, rememberDepartingResultProjection,
    resultDrawsFeature
  } = createDepartedResultIntent({ state, originalLegendOrder, resultPaintDepart: resultPaintRecord.depart });
  const { displayBlocked, showDeclinedPaint, showEditsMadeWhileDrawn } = createDeclinedPaintDisplay({
    state,
    sessionBlocks,
    ruleMatchingPending,
    currentEditorProjectionState,
    projectMountedEditorIntent,
    resultPaintLacking: resultPaintRecord.lacking,
    resultPaintShown: resultPaintRecord.shown
  });
  const {
    rememberCommittedEditorState, projectEditorIntentOnDisplay, takeLabelProjectionResultIdentity, takeDisplayAsks
  } = createDisplayedResultProjection({
    state,
    liveResultIdentities,
    resultIdentitiesOf,
    currentEditorProjectionState,
    departedResultIntent,
    displayedIntentSignature,
    resultDrawsFeature,
    compileDisplayedResultOperations,
    projectMountedEditorIntent,
    showEditorIntent,
    displayBlocked,
    showDeclinedPaint,
    LEGEND_STRUCTURE_DOMAINS,
    STROKE_DOMAINS,
    resultPaintCommit: resultPaintRecord.commit,
    resultPaintRebind: resultPaintRecord.rebind,
    resultPaintDisplay: resultPaintRecord.display,
    resultPaintDrawn: resultPaintRecord.drawn,
    resultPaintShown: resultPaintRecord.shown
  });
  // The Result binder's steps call the editor intent above. No Result mounts
  // during setup: the preview runtime reads the steps only when a mounted
  // Result binds or is released, so they are configured here.
  const {
    adoptLegend, bindComposition, setupDragAffordances, initializeStrokeAndCanvas
  } = createMountedLegendBinding({
    state,
    results,
    liveResultIdentities,
    projectEditorIntentOnDisplay,
    rememberCommittedEditorState,
    isLegendLayoutReady: legendActions.isLegendLayoutReady,
    prepareLegendLayout: legendActions.prepareLegendLayout,
    captureResultInventory: legendActions.captureResultInventory,
    adoptResultInventory: legendActions.adoptResultInventory,
    onLegendGeometryChanged: legendActions.onLegendGeometryChanged,
    extractLegendEntries: legendActions.extractLegendEntries,
    setupLegendDrag: legendActions.setupLegendDrag,
    getResultIdentity: previewRuntime.getResultIdentity,
    commitActiveResultEdit: previewRuntime.commitActiveResultEdit,
    adoptDrawnWrite: history.adoptDrawnWrite,
    captureBaseConfig: legendLayout.captureBaseConfig,
    setupDiagramDrag: legendLayout.setupDiagramDrag,
    applyCanvasPadding: legendLayout.applyCanvasPadding
  });
  const { synchronizeLabelEditor, reconcileSelection, afterReady } = createMountedResultReadiness({
    state,
    isFeatureDrawerMounted,
    clickedFeature,
    labelTextScopeDialog,
    labelReflowProcessing,
    mountedResultRevision,
    takeLabelProjectionResultIdentity,
    takeDisplayAsks,
    currentEditorProjectionState,
    showEditsMadeWhileDrawn,
    syncLabelEditor: featureActions.syncLabelEditor,
    followRestoredSpecificRules: featureActions.followRestoredSpecificRules,
    requestAutomaticRerender: featureActions.requestAutomaticRerender,
    redrawShown: featureActions.redrawShown,
    clearFeatureSelection,
    resultPaintDrewPaint: resultPaintRecord.drewPaint,
    resultPaintDrawn: resultPaintRecord.drawn,
    handleMountedResultReady: previewFeatureSearch.handleMountedResultReady,
    keepRedrawnResult: history.keepRedrawnResult
  });
  previewRuntime.configureMountedResultBinder({
    adoptLegend,
    bindComposition,
    setupDragAffordances,
    installDelegatedInteractions(context) {
      cancelPreviewTransformInteraction({ reconcile: false });
      featureActions.attachSvgFeatureHandlers({
        root: context.root,
        phase: context.phase,
        rootGeneration: context.rootGeneration
      });
      similarityAlignmentPorts.refreshCanvas();
      featureActions.preparePairwiseInteractionAffordances({
        root: context.root,
        phase: context.phase,
        rootGeneration: context.rootGeneration
      });
      featureActions.preparePreviewItemAffordances({ root: context.root });
    },
    synchronizeLabelEditor,
    initializeStrokeAndCanvas,
    reconcileSelection,
    afterReady
  });

  const {
    installArrivingSlot, drawingSwitchBusy, drawingSwitchAvailable, errorDrawingId, showDrawingNow, selectDrawing
  } = createDrawingSwitch({
    state, results, selectedResultIndex, errorLog, failedGeneratePreservedResult, sessionOperationAvailability,
    historyAvailability: history.historyAvailability, historyRestoring: () => history.restoring.value,
    traversalPending: history.traversalPending, settlePendingIntent: history.settlePendingIntent,
    refreshCurrentIntent: history.refreshCurrentIntent, adoptSessionReset: history.adoptSessionReset,
    ruleWritePending: featureActions.ruleWritePending,
    suspendSpecificRulePatternDrafts: featureActions.suspendSpecificRulePatternDrafts,
    clearFeatureSelection, rememberDepartingResultProjection,
    auxiliaryImportFailureOnArrival, resetModeTransientUi, installShownDrawing, dismissNotice: noticeBanner.dismiss,
    switchDrawingSearch: previewFeatureSearch.switchDrawingSearch, rebindAnnotationRecords,
    clearActiveRuntime: previewRuntime.clearActiveRuntime, getResultIdentity: previewRuntime.getResultIdentity,
    presentSelectedResult: previewRuntime.presentSelectedResult,
    rememberResultInventory: legendActions.rememberResultInventory, artifactSlotsSwap: artifactSlots.swap,
    captureArtifactSlot: historySnapshots.captureArtifactSlot, installArtifactSlot: historySnapshots.installArtifactSlot
  });
  // The operations on the drawing list, each one project step (R11, R16-2).
  const drawingList = createDrawingList({
    drawingTabs: () => state.drawingTabs.value,
    shownDrawingId: () => state.activeDrawingId.value,
    switchBusy: drawingSwitchBusy,
    drawingOf: state.drawingOf,
    drawingHasResult: (drawingId) => resultsOf(drawingId).length > 0,
    drawingHasEdits: (drawingId) => drawingHasEdits(state.drawingOf(drawingId)),
    boundFileIdsOf: (drawingId) => boundInputFileIds({ drawingList: () => [state.drawingOf(drawingId)] }),
    createDrawing: (drawingId, mode) => withNewDrawingDefaults(state, createDrawingState(drawingId, mode)),
    styleDrawing: (sourceId, drawingId, mode) => {
      const source = state.drawingOf(sourceId);
      return applyDrawingStyle(withNewDrawingDefaults(state, createDrawingState(drawingId, mode)), source,
        resultsManager.paletteColorsOf(source));
    },
    // The copy's Result is its source's slot under its ID (the shown drawing's through `captureArtifactSlot`, R1).
    copyDrawing: (sourceId, copyId) => {
      const drawing = duplicateDrawingState(state.drawingOf(sourceId), copyId);
      const slot = artifactSlots.readOfDrawing(sourceId, state.activeDrawingId.value,
        () => historySnapshots.captureArtifactSlot(), (kept) => kept);
      return { drawing, held: { slot: slot ? artifactSlotAs(slot, drawing) : null, losatRows: losatRowsOf(sourceId), losatClears: historySnapshots.losatClearCount() } };
    },
    // A region drawing is read as Load reads a stored drawing, then made under its new ID.
    regionDrawing: makeRegionDrawing,
    addDrawing: state.addDrawing,
    // An error that names a drawing leaves with it (R17-29), as on a switch away (OV-372): every Delete, its Redo,
    // and the Undo of New and Duplicate remove a drawing here, which hands over its slot and LOSAT rows (OV-3209).
    removeDrawing: (drawingId, left = []) => {
      if (errorDrawingId() === drawingId) errorLog.value = null;
      return { slot: state.removeDrawing(drawingId, artifactSlots), losatRows: swapLosatResultRows(drawingId, left), losatClears: historySnapshots.losatClearCount() };
    },
    putHeld: (drawingId, { slot, losatRows, losatClears }) => {
      const cleared = historySnapshots.putHeldSlot(drawingId, slot, losatClears);
      return swapLosatResultRows(drawingId, cleared ? [] : losatRows);
    },
    renameDrawing: state.renameDrawing,
    moveDrawing: state.moveDrawing,
    addDrawingHistory: history.addMember,
    removeDrawingHistory: history.removeMember,
    runListStep: history.runListStep,
    showDrawingNow
  });
  // The drawing bar: its tabs with the Linear offer tab, and the tab menu's
  // New Circular, New Linear, Duplicate, New with this style, Rename, Move,
  // and Delete (B4a-5, B4b).
  const drawingBar = createDrawingBar({
    newDrawing: drawingList.newDrawing,
    duplicateDrawing: drawingList.duplicateDrawing,
    newDrawingWithStyle: drawingList.newDrawingWithStyle,
    renameDrawing: drawingList.renameDrawing,
    moveDrawing: drawingList.moveDrawing,
    deleteDrawing: drawingList.deleteDrawing,
    needsDeleteConfirmation: drawingList.needsDeleteConfirmation,
    drawingTabs: () => state.drawingTabs.value,
    shownDrawingId: () => state.activeDrawingId.value,
    modeOf: (drawingId) => state.drawingOf(drawingId).mode,
    askText,
    askChoice
  });
  const drawingBarTabs = computed(() => drawingBar.tabs());
  // Why the tabs wait: the switch's own busy reason (R15).
  const drawingSwitchBusyReason = computed(() => tm(drawingSwitchBusy()?.reason || ''));
  // Session Load: the loaded Session's artifact slots, before its preview
  // mounts. `opening` is the shown drawing's slot when its set is not the
  // committed one (Load installed the committed set itself); each `stashed`
  // slot waits in its drawing. The stash starts again with the loaded Session.
  // The opening slot installs as a switch installs the arriving one.
  /** @param {LoadedArtifactSlots} slots */
  const installLoadedArtifactSlots = ({ opening, stashed }) => {
    artifactSlots.replaceLoaded(stashed);
    if (opening) installArrivingSlot(opening);
  };

  const { requestPaletteChange, selectPalette, requestResetColors } = resultsManager;
  const {
    undoableAction, featureIdentityNoticeSummary, removeUnmatchedFeatureEdits, addFeatureVisibilityRuleWithHistory,
    moveFeatureVisibilityRuleDownWithHistory, moveFeatureVisibilityRuleUpWithHistory, removeFeatureVisibilityRuleWithHistory,
    setFeatureVisibilityRuleFieldWithHistory, setFeatureVisibilityWithHistory, updateClickedFeatureVisibilityWithHistory,
    handleFeatureVisibilityScopeChoiceWithHistory, requestFeatureColorChangeWithHistory, setFeatureColorValueWithHistory,
    updateClickedFeatureColorWithHistory, updateLegendEntryColorShown, updateLegendEntryColorInPopup,
    setLegendEntryStrokeColorValueWithHistory, scopeChoiceWithHistory, handleColorScopeChoiceWithHistory,
    handleFeatureStyleScopeChoiceWithHistory, featureStyleScopeLabelOptions
  } = createFeatureEditHistory({
    state, runUndoable: history.runUndoable, withHistory: dialogChoice.withHistory, editEditorIntent, showingEditorIntent,
    LEGEND_FILL_DOMAINS, STROKE_DOMAINS, addFeatureVisibilityRule, moveFeatureVisibilityRuleDown, moveFeatureVisibilityRuleUp,
    removeFeatureVisibilityRule, setFeatureVisibilityRuleField, setFeatureVisibility, updateClickedFeatureVisibility,
    handleFeatureVisibilityScopeChoice, requestFeatureColorChange, setFeatureColorValue, updateClickedFeatureColor,
    updateLegendEntryColor, legendRowHasRules, setTypeDefaultColor: featureActions.setTypeDefaultColor,
    preparedLegendRowDefaultColorType: featureActions.preparedLegendRowDefaultColorType, setLegendEntryStrokeColorValue,
    cancelFeatureStyleScope, handleColorScopeChoice, handleFeatureStyleScopeChoice
  });
  const {
    handleLegendNameCommitWithHistory, joinLegendGroupWithHistory, legendReorder, handleLegendReorderKey,
    startLegendReorder, overLegendReorder, dropLegendReorder, endLegendReorder, handleLegendRenameChoiceWithHistory,
    handlePaletteColorsChoiceWithHistory, handleResetColorChoiceWithHistory
  } = createLegendEditHistory({
    state, undoableAction, scopeChoiceWithHistory, runUndoableCheckpoint: history.runUndoableCheckpoint,
    handleLegendNameCommit, joinLegendGroup, moveLegendEntry, handleLegendRenameChoice, restoreLegendRows,
    restoreDeletedLegendEntries, cancelLegendRename, handlePaletteColorsChoice: resultsManager.handlePaletteColorsChoice,
    cancelPaletteColorsDialog: resultsManager.cancelPaletteColorsDialog, handleResetColorChoice, cancelResetColor
  });
  const {
    resetClickedFeatureFillColorWithHistory, updateClickedFeatureStrokeWithHistory,
    setClickedFeatureStrokeColorValueWithHistory, setClickedFeatureStrokeWidthValueWithHistory,
    resetClickedFeatureStrokeWithHistory, setFeatureColorWithHistory, selectedFeatureBulkColor, selectedFeatureBulkCaption,
    selectedFeatureBulkVisibility, selectedFeatureBulkStrokeColor, selectedFeatureBulkStrokeWidth, applySelectedFeatureColor,
    applySelectedFeatureVisibility, applySelectedFeatureStroke, openFirstSelectedFeature,
    updateClickedFeatureLabelTextWithHistory, handleLabelTextScopeChoiceWithHistory, handleHiddenLabelTextChoiceWithHistory,
    requestLabelTextChangeByFeatureIdWithHistory, requestLabelTextChangeByKeyWithHistory, commitListLabelDraft,
    resetAllLabelTextOverridesWithHistory, loadLabelOverrideTableWithHistory
  } = createSelectedFeatureHistory({
    state, undoableAction, editEditorIntent, runUndoable: history.runUndoable,
    runUndoableCommand: history.runUndoableCommand, resetClickedFeatureFillColor, updateClickedFeatureStroke,
    setClickedFeatureStrokeColorValue, setClickedFeatureStrokeWidthValue, resetClickedFeatureStroke, setFeatureColor,
    applyColorToSelectedFeatures, syncFeatureSelectionClasses,
    buildSelectedFeaturesVisibilityCommand, clearFeatureSelection,
    applyStrokeToSelectedFeatures, openFeatureEditorForFeature, updateClickedFeatureLabelText, handleLabelTextScopeChoice,
    handleHiddenLabelTextChoice, requestLabelTextChangeByFeatureId, requestLabelTextChangeByKey,
    getEditableLabelByFeatureId: featureActions.getEditableLabelByFeatureId, resetAllLabelTextOverrides,
    loadLabelOverrideTable
  });
  const {
    runInfoCopyStatus, exactReplayCopyStatus, runInfoElapsedText, runInfoReproducibilityText, runInfoHasCliHelperFiles,
    copyRunCommand, copyExactReplayCommand
  } = createRunInfoCopy({ lastRunInfo });

  const runAnalysis = () => {
    const drawing = state.activeDrawing();
    // The phone layout shows the Preview, where the run's progress, Result or error appears; a new
    // Result comes in with its workspace at the pane's top while the Preview is still shown.
    settingsPaneActions.showPhonePane('preview', { workspace: false });
    const patternDrafts = featureActions.captureSpecificRulePatternDrafts();
    return runGeneratedDiagramAnalysis(null, null, null, {
      prepareGenerate: async () => {
        if (mode.value === 'linear') await prepareLinearGenerate(drawing);
        if (similarityAlignmentActions?.validateBeforeGenerate) {
          const validation = await similarityAlignmentActions.validateBeforeGenerate();
          if (validation?.status !== 'ok') {
            failedGeneratePreservedResult.value = results.value.length > 0;
            if (similarityAlignmentActions.dialogOpen.value) {
              await focusSimilarityAlignmentDialog();
            }
            return validation;
          }
        }
        const comparisonPlanSnapshot = mode.value === 'linear'
          ? drawing.linearComparisonResolution.value
          : null;
        const comparisonExecution = importedComparisonExecution({
          intent: drawing.importedComparisonIntent,
          draftResolution: comparisonPlanSnapshot
        });
        if (!comparisonExecution.ok) {
          errorLog.value = normalizeUserFacingError(comparisonExecution.message, {
            operation: 'generate', stage: 'request-validation',
            // Only the Linear Comparison panel offers the choice the code names.
            code: mode.value === 'linear' ? comparisonExecution.fallbackCode : undefined
          });
          generationFailureRecovery.value = results.value.length ? 'preserved' : 'no-result';
          failedGeneratePreservedResult.value = results.value.length > 0;
          if (mode.value === 'linear') await focusLinearComparisonIssue();
          return { status: 'error', error: errorLog.value };
        }
        const missingBlastTsvPairs = comparisonExecution.mode === 'draft'
          ? missingUploadPairsToResolve(comparisonPlanSnapshot)
          : [];
        if (missingBlastTsvPairs.length) {
          missingBlastTsvDialog.pairs = missingBlastTsvPairs;
          missingBlastTsvDialog.drawingId = state.activeDrawingId.value;
          missingBlastTsvDialog.show = true;
          return { status: 'blocked', reason: 'missing-blast-tsv' };
        }
        return { status: 'ready', comparisonPlanSnapshot, comparisonExecution };
      },
      afterGenerate: async (result, runDrawing) => {
        if (result?.status === 'error' && mode.value === 'linear') {
          await focusLinearComparisonIssue();
        }
        // Labels that do not fit drew nothing: the dialog computes the rows the bounded search did not reach.
        if (result?.status === 'error' && result.labelFit) {
          void openLabelFitDialog(result.labelFit.report, {
            rowToken: result.labelFit.rowToken, drawing: runDrawing
          });
        }
        if (result?.status === 'ok') {
          featureActions.clearSpecificRulePatternDrafts();
          await rulePreparation.prepare();
          clearFeatureSelection({ clearStatus: true });
          if (settingsPaneActions.settingsPane.value.pane === 'preview') settingsPaneActions.showPhonePane('preview');
        } else {
          featureActions.restoreSpecificRulePatternDrafts(patternDrafts);
        }
      }
    });
  };
  // E1 (PD-OI-079): an error about another drawing's Result (Save of a
  // legacy Result kept in that drawing) offers Generate there: the action
  // shows that drawing first. The drawing exists: its removal took the error
  // with it (R17-29), and a Load clears every error.
  const generateFromError = () => {
    const target = errorDrawingId();
    if (target !== undefined && target !== state.activeDrawingId.value) {
      const switched = selectDrawing(target);
      if (switched && switched.status === 'busy') return switched;
    }
    return runAnalysis();
  };

  const cancelGeneration = () => cancelRunAnalysis();

  similarityAlignmentActions = createSimilarityAlignmentActions({
    state,
    getOrthogroupById: orthogroupActions.getOrthogroupById,
    getEnrichedOrthogroupMembers: orthogroupActions.getEnrichedOrthogroupMembers,
    selectedOrthogroupDetail: orthogroupActions.selectedOrthogroupDetail,
    getCommittedRequest: getCommittedCanonicalRenderRequest,
    getCommittedSession: getCommittedCanonicalSession,
    projectCommittedAlignment: /** @type {any} */ (projectCommittedSimilarityAlignment),
    // R13: the root runs the candidate with the record display's orientation
    // checkpoint; the alignment owner holds neither owner.
    // The Linear drawing is pinned here, before the run's first await.
    runRecordAlignment: ({ orientations, ...run }) => {
      const drawing = shownDrawingOfMode('linear');
      if (!drawing) throw diagnosticError('INPUT_INVALID', { field: 'drawing', reason: 'LINEAR_DRAWING_REQUIRED' }, { stage: 'record-alignment' });
      return runCommittedCanonicalCandidate({ ...run,
        captureIntentCheckpoint: () => recordDisplayControls.captureAlignmentOrientationIntent(drawing, orientations),
        restoreIntentCheckpoint: (/** @type {AlignmentOrientationCheckpoint[]} */ checkpoint) => recordDisplayControls.restoreAlignmentOrientationIntent(drawing, checkpoint),
        commitIntent: () => recordDisplayControls.commitAlignmentOrientations(drawing, orientations) });
    },
    cancelRunAnalysis,
    runHelperOperation: runDiagramHelperOperation,
    resolveOperation: DIAGRAM_HELPER_OPERATIONS.RESOLVE_SIMILARITY_ALIGNMENT,
    getCurrentSvg: () => svgContainer.value?.querySelector?.('svg') || null,
    previewCandidate: featureActions.previewAlignmentCandidate,
    clearCandidatePreview: featureActions.clearAlignmentCandidatePreview,
    onError: (error) => { errorLog.value = error; }
  });
  recordDragAlignmentPorts.beforeRecordDrag = similarityAlignmentActions.beforeRecordDrag;
  recordDragAlignmentPorts.afterRecordDrag = similarityAlignmentActions.afterRecordDrag;
  // The comparison plan's edits clear the alignment candidate: created once that owner exists.
  const {
    replaceLinearComparisonPlan, mutateLinearComparisonPlan, syncLinearRecordLayout, setLinearRecordRow,
    setLinearRecordLayoutEnabled, moveLinearRecordWithinRow
  } = createLinearComparisonPlanWriter({
    shownAction, shownDrawingOfMode, runUndoable: history.runUndoable,
    clearAlignmentForMutation: similarityAlignmentActions.clearForMutation, invalidateLinearComparisonArtifacts,
    materializeAutomaticLinearRecords
  });
  const {
    runLinearComparisonIssueAction, setLinearComparisonGlobalAction, setLinearComparisonLosatMode,
    setLinearComparisonLosatpMode
  } = createLinearComparisonModes({
    shownAction, runUndoable: history.runUndoable, mutateLinearComparisonPlan, openLinearComparisonDisclosure,
    chooseLinearComparisonPairFile, focusLinearComparisonIssue,
    clearAlignmentForMutation: similarityAlignmentActions.clearForMutation, invalidateLinearComparisonArtifacts
  });
  const {
    addLinearComparison, omitLinearComparisonPairs, omitLinearComparison, clearSelectedLinearComparisons,
    setLinearComparisonEndpoint, setLinearComparisonSource, setLinearComparisonFile, swapLinearComparisonFile,
    reuseLinearComparisonFile, deactivateLinearComparisonFile
  } = createLinearComparisonPairEdits({
    state, shownAction, replaceLinearComparisonPlan, syncLinearRecordLayout, focusLinearComparisonPair,
    linearComparisonPairForEdgeKey
  });
  const {
    setLinearComparisonLosatFilename, reuseLinearComparisonLosatFilename, deactivateLinearComparisonLosatFilename,
    setResolvedLinearComparisonLosatFilename, reuseResolvedLinearComparisonLosatFilename,
    deactivateResolvedLinearComparisonLosatFilename, setWindowPairLosatFilename, addLinearComparisonBatch,
    setLinearComparisonGapAction, setLinearComparisonCardFile
  } = createLinearComparisonGapEdits({
    sessionOperationAvailability, shownAction, replaceLinearComparisonPlan, syncLinearRecordLayout,
    linearComparisonPairForEdgeKey, omitLinearComparisonPairs, setLinearComparisonFile,
    windowPairEdgeKeysOf: (cardEdgeKey) => (linearLosatWindowPairs.value[cardEdgeKey] ?? []).map(({ edgeKey }) => edgeKey)
  });
  // The saved comparison's choice replaces a plan, which clears the alignment candidate: created once that owner exists.
  const {
    importedComparisonActionLabels, inheritImportedComparison, replaceImportedComparison, clearImportedComparison,
    importedComparisonNeedsResolution, importedComparisonActionText, importedComparisonCanInherit, importedComparisonCanReplace
  } = createImportedComparisonChoice({ state, runUndoable: history.runUndoable, errorLog, replaceLinearComparisonPlan });
  // Align at feature and Reset record positions write the record translations as one step without the
  // alignment plan, as a manual record movement does (PD-OI-028). The display starts Align at feature writes for the
  // records it rotated go into those records' drafts in the same step (R13: the root runs the candidate with the
  // record display's draft checkpoint), so Generate draws them and one Undo restores the drafts with the Result. The
  // records are read first, as the Rotate popup reads them, when a Session Load left them unread.
  /** @type {RecordPlacementRun} */
  const runRecordPlacement = async ({ canonical, label, starts }) => {
    const drawing = shownDrawingOfMode('linear');
    if (!drawing) throw diagnosticError('INPUT_INVALID', { field: 'drawing', reason: 'LINEAR_DRAWING_REQUIRED' }, { stage: 'record-alignment' });
    if (starts.length) await linearRecordSelector.refresh(drawing);
    const targets = recordDisplayControls.recordStartTargets(drawing, starts);
    return runCommittedCanonicalCandidate({ canonical, label, alignmentResetReceipt: null,
      captureIntentCheckpoint: () => recordDisplayControls.captureRecordStarts(targets),
      restoreIntentCheckpoint: (/** @type {TargetDraftCheckpoint[]} */ checkpoint) => recordDisplayControls.restoreRecordStarts(checkpoint),
      commitIntent: () => recordDisplayControls.commitRecordStarts(drawing, targets) });
  };
  const featureAnchor = createFeatureAnchorAction({
    materializedSession: similarityAlignmentActions.materializedSession,
    hasAlignmentPlan: () => Boolean(state.similarityAlignmentPlan.value),
    runHelperOperation: runDiagramHelperOperation,
    resolveOperation: DIAGRAM_HELPER_OPERATIONS.RESOLVE_FEATURE_ANCHOR,
    runPlacement: runRecordPlacement,
    notePlanReplaced: similarityAlignmentActions.notePlanReplaced,
    projectPlacement: projectCommittedSimilarityAlignment,
    committedRequest: getCommittedCanonicalRenderRequest,
    sessionReplacements: () => sessionReplacements.value
  });
  // Reset record positions writes x 0 for every record through the channel Align at feature's Apply writes through.
  const recordPositionReset = createRecordPositionReset({
    materializedSession: similarityAlignmentActions.materializedSession,
    hasAlignmentPlan: () => Boolean(state.similarityAlignmentPlan.value),
    notePlanReplaced: similarityAlignmentActions.notePlanReplaced,
    runPlacement: runRecordPlacement,
    projectPlacement: projectCommittedSimilarityAlignment,
    onError: (error) => { errorLog.value = error; }
  });
  const { shownSize, blockerInputs } = createFigureSizeReaders({
    state, resultEditRevision, mountedResultRevision, shownCommittedRequest, decorationsMoved, linearLabelHasSharedRow
  });
  const figureSizeField = createFigureSizeField({
    state,
    watch,
    shownSize
  });
  const linearWrapField = createLinearWrapField({
    state,
    changeLayoutSetting: featureActions.placementActions.changeLayoutSetting,
    shownWrap: () => shownLinearWrapReading.value,
    blockerInputs,
    watch
  });
  const { missingBlastTsvDialog, missingBlastTsvRows, resolveMissingBlastTsv } = createMissingBlastTsvChoice({
    state, withHistory: dialogChoice.withHistory, closeAfterChoice: dialogChoice.closeAfterChoice, choicePending: dialogChoice.pending,
    omitLinearComparisonPairs, linearComparisonPairForEdgeKey, runAnalysis
  });
  const {
    reviewDrawingFigure, dropSimilarityAlignmentView, fitPreviewForSimilarityAlignment, endSimilarityAlignmentTakeover,
    rememberSimilarityAlignmentInvoker, restoreSimilarityAlignmentFocus, focusSimilarityAlignmentDialog,
    applySimilarityAlignmentDialog
  } = createSimilarityAlignmentTakeover({
    state, similarityAlignmentTakeover, similarityAlignmentHoldShown, releasePreviewViewport, fitPreviewToViewport,
    getResultIdentity: previewRuntime.getResultIdentity, readOfDrawing: artifactSlots.readOfDrawing,
    openRightDrawerTab: rightDrawerActions.openRightDrawerTab,
    showFeaturePopupAtFeature: featureActions.showFeaturePopupAtFeature,
    alignmentDialogOpen: similarityAlignmentActions.dialogOpen, applyAlignmentReset: similarityAlignmentActions.applyReset,
    applyAlignmentDraft: similarityAlignmentActions.applyDraft
  });
  const {
    similarityAlignmentCanvasHover, refreshSimilarityAlignmentCanvas, similarityAlignmentReview,
    similarityAlignmentEditorDisabledReason, openSimilarityAlignmentReset, finishSimilarityAlignmentStart,
    startSimilarityAlignmentFromDrawer, cancelSimilarityAlignmentDialog
  } = createSimilarityAlignmentReview({
    state, resultEditRevision, similarityAlignmentTakeover, holdPreviewViewport,
    clearAlignmentOverlay: featureActions.clearAlignmentOverlay, showAlignmentOverlay: featureActions.showAlignmentOverlay,
    closeRightDrawer: rightDrawerActions.closeRightDrawer, closeSettingsPopup: settingsPaneActions.closeSettingsPopup,
    alignmentDialogOpen: similarityAlignmentActions.dialogOpen, alignmentDraft: similarityAlignmentActions.draft,
    alignmentStatus: similarityAlignmentActions.status,
    isDraftArtifactCurrent: similarityAlignmentActions.isDraftArtifactCurrent,
    selectAlignmentCandidate: similarityAlignmentActions.selectCandidate,
    alignmentDirectionPreview: similarityAlignmentActions.directionPreview,
    alignmentResetDialogOpen: similarityAlignmentActions.resetDialogOpen,
    openAlignmentReset: similarityAlignmentActions.openReset, alignmentError: similarityAlignmentActions.error,
    startAlignmentFromDrawer: similarityAlignmentActions.startFromDrawer, cancelAlignment: similarityAlignmentActions.cancel,
    reviewDrawingFigure, dropSimilarityAlignmentView, fitPreviewForSimilarityAlignment, endSimilarityAlignmentTakeover,
    rememberSimilarityAlignmentInvoker, restoreSimilarityAlignmentFocus, focusSimilarityAlignmentDialog
  });
  similarityAlignmentPorts.refreshCanvas = refreshSimilarityAlignmentCanvas;
  similarityAlignmentPorts.editorDisabledReason = () => similarityAlignmentEditorDisabledReason.value;
  const {
    canUseClickedOrthogroupActions, clickedOrthogroupDetail, alignByClickedOrthogroup, highlightClickedOrthogroup,
    clearOrthogroupHighlight, openOrthogroupInDrawer, reselectSimilarityAlignmentReference, openClickedOrthogroupInEditor
  } = createOrthogroupPopup({
    state, displayedLinearDrawing, getOrthogroupById: orthogroupActions.getOrthogroupById,
    getEnrichedOrthogroupMembers: orthogroupActions.getEnrichedOrthogroupMembers,
    groupOrthogroupMembersByRecord: orthogroupActions.groupOrthogroupMembersByRecord,
    resolveOrthogroupName: orthogroupActions.resolveOrthogroupName,
    getOrthogroupSequenceCount: orthogroupActions.getOrthogroupSequenceCount,
    orthogroupScope: orthogroupActions.orthogroupScope, orthogroupScopeLabel: orthogroupActions.orthogroupScopeLabel,
    highlightOrthogroupById: orthogroupActions.highlightOrthogroupById,
    removeOrthogroupHighlight: orthogroupActions.clearOrthogroupHighlight, selectOrthogroup: orthogroupActions.selectOrthogroup,
    startAlignmentFromPopup: similarityAlignmentActions.startFromPopup,
    setAlignmentDrawerReference: similarityAlignmentActions.setDrawerReference,
    alignmentRepair: similarityAlignmentActions.repair, openRightDrawerTab: rightDrawerActions.openRightDrawerTab,
    rememberSimilarityAlignmentInvoker, finishSimilarityAlignmentStart
  });

  const { resetAllPositions, resetCanvasPadding } = legendLayout;

  // Reset Layout has something to reset once an element moved, the
  // canvas padding changed, or the view was panned or zoomed.
  const layoutResetAvailable = computed(() => (
    decorationsMoved() || Boolean(Number(canvasPan.x) || Number(canvasPan.y)) || zoom.value !== 1 || canvasPadded(state.activeDrawing())
  ));

  const {
    featurePopupStyle, pairwiseMatchPopupStyle, settingsPopupStyle, featurePopupDragging, pairwiseMatchPopupDragging,
    settingsPopupDrag, keepSettingsPopupSectionOpen, pairwiseGroupView, selectPairwiseBlockOrthogroup, openPairwiseFeatureRow,
    startFeaturePopupResize, startPairwiseMatchPopupResize
  } = createResultPopups({
    state, clickedFeature, clickedFeaturePos, clickedPairwiseMatch, clickedPairwiseMatchPos, richFeaturePopup,
    featurePopupRef, featurePopupDrag, featurePopupSize, featurePopupResize,
    pairwiseMatchPopupRef, pairwiseMatchPopupDrag, pairwiseMatchPopupSize, pairwiseMatchPopupResize,
    selectedPairwiseBlockOrthogroupId, moveSettingsPopup: settingsPaneActions.moveSettingsPopup,
    openFeatureEditorForFeature, markSelectedPairwiseMatch: featureActions.markSelectedPairwiseMatch
  });

  const {
    normalizeSessionTitle, errorDisplay, sessionSaveNeedsGenerate, reloadAfterOperationError, sessionTitleLabel
  } = createSessionStatus({ state });
  const { canUseLinearRulerOnAxis, canUseCircularScaleStyling } = createScaleStyleAvailability({
    state, circularTrackSlotEffectiveEnabled: circularTrackSlotEditor.circularTrackSlotEffectiveEnabled
  });

  const {
    clickedFeatureHeading, clickedFeatureSummary, copyClickedFeatureDetails, copyClickedFeatureQualifiers
  } = createClickedFeatureCopy({ clickedFeature });

  const {
    downloadText, canRetryInteractiveSvgExport, downloadSVG, downloadInteractiveSVG, downloadPNG, downloadPDF
  } = createResultExport({
    state, errorLog, showNotice, captureSvgExport, loadExportService, runDiagramHelperOperation, DIAGRAM_HELPER_OPERATIONS,
    displayedArtifactDrawing: historySnapshots.displayedArtifactDrawing, embedInExports: fontSettings.embedInExports,
    shownInputSummary, shownLinearWrapLeaf, shownCommittedRequest, awaitLinearRecordDiscovery, committedLinearRecordsOf,
    readEmbedJson: catalogSequences.readEmbedJson
  });

  const { popupLegendGroups, featureVisibilityFeatureSuggestions } = createFeatureChoiceReaders({
    state, shownLegendName: legendActions.shownLegendName
  });

  // The Session title is asked in the app's dialog; Save asks for it when the Session has none.
  /** @param {{ title: string, actionLabel: string }} request */
  const askSessionTitle = ({ title, actionLabel }) => askText({
    title, label: t('Session title'), value: normalizeSessionTitle(sessionTitle.value), actionLabel
  });

  const editSessionTitle = async () => {
    const busy = sessionOperationAvailability();
    if (busy) return busy;
    const input = await askSessionTitle({ title: t('Rename Session'), actionLabel: t('Rename') });
    if (input === null) return;
    // The title is part of the History snapshot: a rename is its own step.
    const title = normalizeSessionTitle(input);
    if (title === normalizeSessionTitle(sessionTitle.value)) return;
    return history.runUndoable(msg('Rename Session'), () => { sessionTitle.value = title; }, { detail: () => ({ value: title }) });
  };

  const saveSessionWithTitle = () => exportSession(null, {
    shownDrawing: state.activeDrawing, planLinearUnits,
    // A record read that finished while the Save was pending was deferred: the
    // end of the Save calls each drawing's discovery again (R10).
    onSettled: () => state.drawingList().forEach((drawing) => { void discoverDrawingRecords(drawing); }),
    availability: sessionSaveLoadAvailability,
    askChoice,
    recordDisplayRows: recordDisplayControls.allRows,
    // Save writes every Result: the kept artifacts of the drawings not shown go beside the shown one (E1).
    readKeptArtifacts: () => artifactSlots.kept().filter((slot) => slot.drawingId !== state.activeDrawingId.value),
    resolveTitle: async () => {
      let title = normalizeSessionTitle(sessionTitle.value);
      if (!title) {
        const input = await askSessionTitle({ title: t('Save Session'), actionLabel: t('Save') });
        if (input === null) {
          recordSessionLifecycleEvent('session-save-title-canceled');
          return null;
        }
        title = normalizeSessionTitle(input);
        sessionTitle.value = title;
      }
      return title;
    },
    beforeExport: async (/** @type {{ draftRequest: boolean }} */ { draftRequest }) => {
      // The drawing Save writes, pinned before the first await.
      const drawing = state.activeDrawing();
      await nextTick();
      await afterPaint();
      recordSessionLifecycleEvent('session-save-paint-opportunity-completed');
      // Save writes the displayed Result as Generate draws it: it waits for
      // the Result being displayed, then shows what its display declined.
      await previewRuntime.pendingReadiness()?.catch(() => {});
      const displayed = previewRuntime.getActiveRuntime()?.resultIdentity;
      if (displayed) await showDeclinedPaint(displayed, { saving: true });
      recordSessionLifecycleEvent('session-save-catalog-preparation-start');
      /** @type {Awaited<ReturnType<typeof prepareLinearRecordCatalog>>['catalog']} */
      let catalog = null;
      /** @type {Awaited<ReturnType<typeof prepareLinearRecordCatalog>>['error']} */
      let error = '';
      if (draftRequest) {
        ({ catalog, error } = await prepareLinearRecordCatalog({ privateCandidate: true }));
        await afterPaint();
      }
      recordSessionLifecycleEvent('session-save-catalog-preparation-end', {
        reusedCommittedSession: !draftRequest
      });
      if (error) throw error;
      return {
        linearRecordCatalog: catalog,
        omittedAnnotationRecordKeys: draftRequest ? getAnnotationRecordCatalog(drawing).omittedRecordKeys || [] : []
      };
    },
    onError: (error) => { errorLog.value = normalizeUserFacingError(error); },
    // UJ-09: a saved download is the History position Load Session compares with.
    // The notice names the file the download started.
    onSaved: (/** @type {{ filename: string }} */ { filename }) => {
      history.markSavePoint();
      // D-02: a font gbdraw does not bundle is saved by name only (receipt 〔警告〕).
      showNotice({
        text: t('Session saved: {file}.', { file: filename }),
        // Each drawing's font, and the font of the shown Result, whose committed
        // request the Session saves too (OV-742).
        items: customFontSaveWarnings(
          [
            ...state.drawingList().map((drawing) => drawing.adv.font_family),
            requestFontFamily(getCommittedCanonicalRenderRequest())
          ],
          webFontFamilyNames(state.loadedFonts.value)
        )
      });
    }
  });

  // OV-503: Edit in the Editor drawer's list opens the popup beside the drawer
  // and the toggle on the drawer's edge (OV-505). Without an event the popup
  // is centred.
  /** @param {Record<string, unknown>} feat @param {{ clientX?: number, clientY?: number, currentTarget?: EventTarget | null } | null} [event] */
  const openFeatureEditorFromList = (feat, event = null) => {
    const drawer = event?.currentTarget instanceof Element ? event.currentTarget.closest('.right-drawer') : null;
    const rects = [drawer, drawer?.parentElement?.querySelector(':scope > .drawer-toggle')]
      .flatMap((element) => (element ? [element.getBoundingClientRect()] : []));
    return openFeatureEditorForFeature(feat, event && {
      clientX: event.clientX,
      clientY: event.clientY,
      panel: rects.length ? {
        left: Math.min(...rects.map((rect) => rect.left)), top: Math.min(...rects.map((rect) => rect.top)),
        right: Math.max(...rects.map((rect) => rect.right)), bottom: Math.max(...rects.map((rect) => rect.bottom))
      } : null
    });
  };

  const {
    getCircularRecordOrderLabel, circularRecordDiscoveryState, circularRecordInspectionEnabled, inspectCircularSourceRecords,
    circularDrawnRecordCount, circularRecordPresentationOptions, circularRecordPresentationError, circularOwnRecordBlock,
    circularSingleRecordPresentationEnabled, circularRecordLayout
  } = createCircularRecordPresentation({
    state, mode, semanticMutationAvailable, sessionOperationAvailability, shownDrawingOfMode, refreshCircularRecordOrder, shownRecordIds
  });

  // The legend-entry and record-definition popups (match popup frame): views
  // of the displayed Result and the settings; their edits are the Legend
  // panel's and the sidebar's.
  const previewItemPopup = createPreviewItemPopup({
    state,
    computed,
    revealSetting: settingsPaneActions.revealSetting,
    readLegendEntryColor: legendEntryColor,
    readShownLegendName: legendActions.shownLegendName,
    readLegendRowRules: (entry) => legendRowRules(legendEntryKey(entry), displayedLegendRowContext(state, state.activeDrawing())),
    countLegendCaptionFeatures: featureActions.countLegendCaptionFeatures,
    readLegendTrackKind: (entry) => displayedTrackLegendRows(state.activeDrawing()).get(legendEntryKey(entry))?.kind ?? null,
    readLegendRowDefaultColorType: (entry) => featureActions.legendRowDefaultColorType(legendEntryKey(entry)),
    readUserDefaultColor: (type) => resultsManager.readUserDefaultColor(state.activeDrawing(), type),
    getCommittedRequest: getCommittedCanonicalRenderRequest,
    readAlignmentInspector: () => similarityAlignmentActions.activePlanInspector.value,
    readCommittedRecordRows: recordDisplayControls.committedRecordRows,
    readCircularOwnRecordBlock: () => circularOwnRecordBlock.value,
    runUndoable: history.runUndoable,
    sessionBusy: sessionOperationAvailability,
    // The popup sits outside the app root that records the Legend panel's text
    // edits, so its rename records its own step.
    renameLegendRow: undoableAction(msg('Rename legend item'), renameLegendEntry),
    legendRenameSettled: () => (legendRenameDialog.show ? new Promise((resolve) => {
      const stop = watch(() => legendRenameDialog.show, (/** @type {unknown} */ shown) => { if (!shown) { stop(); resolve(undefined); } });
    }) : Promise.resolve())
  });

  // Labels that do not fit (LP-D43, LP-D48): the dialog over a Generate that
  // drew nothing, or over crossing leaders from the notice bar.
  // The dialog's state and the actions its template calls go to the template as they are.
  const { openLabelFitDialog, ...labelFitTemplate } = createLabelFitDialog({
    reactive,
    computeRow: async (rowToken, kind) => (await runDiagramHelperOperation(
      DIAGRAM_HELPER_OPERATIONS.COMPUTE_LABEL_FIT_ROW, { rowToken, kind }
    )).result,
    writeRowSettings: ({ values, trackOptions }, drawing) => applyLabelFitValues(drawing, values, drawing.mode, trackOptions),
    withHistory: dialogChoice.withHistory,
    closeAfterChoice: dialogChoice.closeAfterChoice,
    choicePending: () => dialogChoice.pending.value,
    runAnalysis: () => runAnalysis(),
    errorSummary: (error) => errorText(normalizeCaughtError(error, { operation: 'computeLabelFitRow', stage: 'helper' }).summaryParts),
    // The record-definition popup (created above) owns which field writes a record's name lines.
    showRecordText: (recordIndex) => previewItemPopup.showRecordTextSettings(recordIndex)
  });
  // The notice bar over the shown Result drawn with crossing leaders or labels
  // over the definition text (WARN-NOTICE, DEF-WARN).
  const labelFitNotice = computed(() => {
    const entry = /** @type {readonly ResultLabelFit[]} */ (labelFitReports.value).find((item) => item.resultIndex === selectedResultIndex.value
      && item.resultName === results.value[selectedResultIndex.value]?.name);
    return entry && entry.report.drawn ? entry : null;
  });
  // A live Feature visibility edit no render drew since: the summary counts what the render drew,
  // so it is compared with the visibility of the drawing whose Result is shown.
  const shownInputSummaryStale = computed(() => {
    const summary = shownInputSummary.value;
    return summary !== null
      && summary.drawnVisibility !== visibilityPaintState(historySnapshots.displayedArtifactDrawing());
  });
  const openLabelFitNotice = () => {
    const entry = labelFitNotice.value;
    if (!entry) return undefined;
    return openLabelFitDialog(entry.report, { rowToken: entry.rowToken, drawing: historySnapshots.displayedArtifactDrawing() });
  };

  const {
    circularRecordSelectionEnabled, setCircularRecordPresentationSelector, showCircularCanvasSetting,
    showCircularTicksSlot, waitForCircularRecordScrollAnchor
  } = createCircularRecordDisplay({
    state, semanticMutationAvailable, circularRecordDiscoveryState, circularSingleRecordPresentationEnabled,
    sessionOperationAvailability, shownDrawingOfMode, rebindAnnotationRecords,
    revealSetting: settingsPaneActions.revealSetting, circularTrackSlotsPanelOpen, circularRecordPresentationPanel,
    afterPaint
  });
  const { sourceReplacementRetiresPresentation, replaceCircularSourceFile, setCircularSourceFile, addCircularSourceFiles,
    removeCircularSourceOf, removeCircularSourceCard, canMoveCircularSource, moveCircularSourceCard } = createCircularSourceEdits({
    state, CIRCULAR_RECORD_FORM_DEFAULTS, changeRecordSources, pendingDiscoveryOf, discoverDrawingRecords, shownAction, runUndoable: history.runUndoable
  });
  const circularFileCards = createCircularFileCards({ shownDrawingOfMode, circularRecordDiscoveryState, addCircularSourceFiles, removeCircularSourceCard });

  const {
    getCircularRecordRow, getCircularRecordRowOptions, setCircularRecordRow, canMoveCircularRecordOrderUp,
    canMoveCircularRecordOrderDown, moveCircularRecordOrderUp, moveCircularRecordOrderDown, resetCircularRecordOrder
  } = createCircularRecordOrder({ shownAction, shownDrawingOfMode });

  // The Linear File edits: created once the alignment owner and the comparison plan writer exist.
  const { applyLinearSeqMutation, addLinearSeq, removeLinearSource, askLinearFileRemoval } = createLinearSourceEdits({
    sessionOperationAvailability, shownAction, retireLegendStylesOfUnnamedCaptions,
    changeLinearDepthSources: linearTrackSlotEditor.changeLinearDepthSources,
    linearDepthSourceLoss: linearTrackSlotEditor.linearDepthSourceLoss,
    changeRecordSources, askingStackPorts, pendingDiscoveryOf, replaceLinearComparisonPlan,
    invalidateLinearComparisonArtifacts, linearReorderNotice,
    retainAlignmentForStableReorder: similarityAlignmentActions.retainForStableReorder,
    clearAlignmentForMutation: similarityAlignmentActions.clearForMutation, rebindAnnotationRecords, discoverDrawingRecords
  });
  // Reset Settings and Reset Layout: created once the Linear File edits exist (a Linear
  // drawing's reset applies `applyLinearSeqMutation`); only the return reads them.
  const { resetDrawing } = createDrawingReset({
    state, circularTrackNewRenderer, linearTrackNewRenderer, rederiveCircularConservation, applyLinearSeqMutation,
    invalidateLinearComparisonArtifacts, writeArtifactValues: historySnapshots.writeArtifactValues
  });
  const {
    requestResetSettings, resetSettings, resetShownDrawing, setPaletteInstantPreview, resetLayout
  } = createResetSettings({
    state, sessionOperationAvailability, askChoice, resetDrawing, holdAutomaticRerender, followRuleEdits,
    svgContainer, appliedPaletteColors, originalLegendOrder, projectMountedEditorIntent, drawReturningLegendRows,
    STROKE_DOMAINS, resetAllPositions, resetCanvasPadding, resetPreviewViewport, dismissNotice: noticeBanner.dismiss,
    closeSettingsPopup: settingsPaneActions.closeSettingsPopup,
    clearSpecificRulePatternDrafts: featureActions.clearSpecificRulePatternDrafts,
    extractLegendEntries: legendActions.extractLegendEntries, runProjectStep: history.runProjectStep,
    runUndoable: history.runUndoable,
    storePaletteInstantPreview: resultsManager.setPaletteInstantPreview
  });
  const {
    cancelLinearSourceRemoval, applyLinearSourceRemoval, requestLinearSourceRemoval, removeLastLinearSeq
  } = createLinearSourceRemovalDialog({
    linearSourceGroups, linearSourceRemovalDialog, linearSourceRemovalReturnFocus, sessionOperationAvailability,
    shownDrawingOfMode, runUndoable: history.runUndoable, removeLinearSource, askLinearFileRemoval
  });
  const { replaceLinearSourceFile, setLinearSeqPrimaryFile } = createLinearSourceFileReplacement({
    state, shownAction, pendingDiscoveryOf, applyLinearSeqMutation, askLinearFileRemoval, sourceReplacementRetiresPresentation
  });
  const { applyInputReplacementToDrawings, askScope } = createInputReplacementPorts({
    state, replaceCircularSourceFile, replaceCircularDepthFile, replaceCircularConservationBlasts, replaceLinearSourceFile,
    setLinearDepthFiles, swapLinearComparisonFile, barTabs: drawingBar.tabs, askChoice
  });
  // The replacement of an input file for the shown drawing or every drawing
  // that uses it (PD-OI-090 〔入力ファイル〕): its owner is `app/drawing-list.js`,
  // created after every writer it reaches.
  const inputReplacement = createInputReplacement({
    drawingTabs: () => state.drawingTabs.value,
    shownDrawingId: () => state.activeDrawingId.value,
    switchBusy: drawingSwitchBusy,
    drawingOf: state.drawingOf,
    runProjectStep: history.runProjectStep,
    runDrawingStep: history.runUndoable,
    applyInputReplacementToDrawings,
    askScope,
    // The busy outcome has no control to disable: the notice banner words the
    // drawing list's own busy reason, as the tabs' title does (R15).
    onBusy: (busy) => showNotice({ text: tm(busy.reason), tone: 'warning' })
  });
  // The uploaders of the shown drawing's bound input files: each names the
  // file its control binds and its own write; `replaceFromShownDrawing`
  // decides whether to ask (PD-OI-090 〔入力ファイル〕). The five cards and
  // the similarity-ring file list are not replacements of one bound file.
  const { replaceFromShownDrawing } = inputReplacement;
  const {
    uploadCircularSourceFile, uploadCircularDepthFile, uploadLinearSeqPrimaryFile, uploadLinearDepthFile,
    uploadLinearSourceDepthFile, uploadLinearComparisonCardFile
  } = createInputUploads({
    shownAction, replaceFromShownDrawing, replaceCircularSourceFile, replaceCircularDepthFile, replaceLinearSourceFile,
    setLinearDepthFiles, setLinearSourceDepthFile, linearComparisonPairForEdgeKey, setLinearComparisonCardFile
  });

  const { circularIdentitySources, placeCircularRecordChoices, recordSettingsOwners, recordEditData, recordSourceOf } = createRecordSettingsSources({
    shownDrawingOfMode, linearSourceGroups, linearSourceFileName, linearCardRecord,
    linearAnnotationSource, getCircularRecordDiscoveryState, circularAnnotationSourceKey,
    takeBackRuleRows: featureActions.takeBackRuleRows
  });
  const { deleteRecordSettings, removeRecordSource, askRecordFileRemoval, readShownRecords } = createRecordSettingsTransitions({
    state, shownDrawingOfMode, recordSettingsOwners, recordEditData, retireLegendStylesOfUnnamedCaptions,
    projectPaletteAndRules, applyLinearSeqMutation, removeCircularSourceOf, removeLinearSource, askLinearFileRemoval,
    dialogChoicePending: dialogChoice.pending, mutationPending: history.mutationPending,
    refreshLinearRecords: linearRecordSelector.refresh, getCircularRecordDiscoveryState, inspectCircularSourceRecords
  });
  const recordSelection = createRecordSelection({
    reactive,
    computed,
    // The OFF list of the shown drawing when it is of `listMode` (every port
    // here acts on the shown drawing; the lists and checkboxes are in its panel).
    recordsOff: (listMode) => {
      const drawing = shownDrawingOfMode(listMode);
      return drawing ? drawing.recordsOff : [];
    },
    source: recordSourceOf,
    runUndoable: history.runUndoable,
    withDialogChoice: dialogChoice.withHistory,
    closeAfterDialogChoice: dialogChoice.closeAfterChoice,
    // A toggle invalidates comparison results and the alignment plan until Generate (P-9).
    afterRecordSetChange: (listMode) => {
      rebindAnnotationRecords();
      if (listMode !== 'linear') return;
      const drawing = shownDrawingOfMode('linear');
      if (drawing) invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo: true });
      similarityAlignmentActions?.clearForMutation?.(msg('records to draw changed.'));
    },
    removeSource: removeRecordSource,
    askRemoveSource: askRecordFileRemoval,
    deleteRecordSettings,
    prepareRecords: readShownRecords,
    autoOpenRequests: () => pendingDiscoveryOf(state.activeDrawing()).listRequests
  });

  const { linearSourceMoveBlockedReason, canMoveLinearSource, moveLinearSource } = createLinearSourceMoves({
    state, linearSourceGroups, shownAction, shownDrawingOfMode, runUndoable: history.runUndoable, applyLinearSeqMutation
  });
  const {
    setLinearRecordSelector, setLinearCardWindow, addLinearCardWindow, removeLinearCardWindow, resetLinearRecordDefinition, resetLinearRecordSubtitle
  } = createLinearRecordEdits({
    sessionOperationAvailability, shownDrawingOfMode, runUndoable: history.runUndoable,
    clearAlignmentForMutation: similarityAlignmentActions.clearForMutation,
    linearRecordsFor: linearRecordSelector.recordsFor, linearRecordStatusFor: linearRecordSelector.statusFor,
    rebindAnnotationRecords, applyLinearSeqMutation,
    shownUnitPlaces: () => catalogUnitPlaces(state.featureCatalog.value)
  });
  const { setCircularInputType, setLinearInputType } = createInputTypeEdits({
    sessionOperationAvailability, shownDrawingOfMode, changeRecordSources, discoverDrawingRecords,
    clearAlignmentForMutation: similarityAlignmentActions.clearForMutation
  });

  const { labelSettingsShown, drawingMember, defaultColorGroups } = createShownDrawingMembers({
    state, linearComparisonGlobalAction, userDefaultColorRows: resultsManager.userDefaultColorRows
  });

  // New drawing from a region and Find in other files: the shown drawing's draft request, as Save stores it, goes to Python.
  const projectShownRequest = async () => {
    const drawing = state.activeDrawing();
    const { catalog, error } = drawing.mode === 'linear'
      ? await prepareLinearRecordCatalog({ privateCandidate: true })
      : { catalog: null, error: '' };
    if (error) throw error;
    return projectDraftRequest(state, drawing, {
      linearRecordCatalog: catalog, planLinearUnits,
      recordDisplayRows: recordDisplayControls.allRows.value,
      omittedAnnotationRecordKeys: getAnnotationRecordCatalog(drawing).omittedRecordKeys || []
    });
  };
  const regionDrawings = createRegionDrawing({
    shownDrawingId: () => state.activeDrawingId.value,
    projectShownRequest,
    // The projection's resource table holds the descriptors Save would write.
    deriveRegion: (payload) => runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.DERIVE_REGION_DRAWING, { ...payload,
      resources: /** @type {GenerationPayload['resources']} */ (payload.resources) }),
    newDrawingFromRegion: drawingList.newDrawingFromRegion,
    readRecordSources: () => readDrawingRecordIdentities(state.activeDrawing()),
    resultFeatures: shownResultFeatures
  });
  const findInOtherFiles = createFindInOtherFiles({
    shownDrawingId: () => state.activeDrawingId.value, projectShownRequest, runInOverlay: runInLosatOverlay,
    runLocateStep: (payload) => runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.FIND_IN_OTHER_FILES, payload),
    shownDrawing: () => state.activeDrawing(), newDrawingFromSearch: drawingList.newDrawingFromRegion
  });
  /** Reads the record catalog of `drawing` when it is unread (OV-373); one that cannot be read throws (OV-2562). @param {DrawingState} drawing */
  const readRecordCatalog = async (drawing) => {
    if (getAnnotationRecordCatalog(drawing).status !== 'ready') await readShownRecords(drawing.mode);
    const catalog = getAnnotationRecordCatalog(drawing);
    if (catalog.status !== 'ready') throw recordCatalogError(drawing, catalog);
  };
  // How a row binds a record of the shown Result (Add to annotation set, Sequence search, its live draw; PD-OI-117).
  const shownResultRecordBindings = () => resultRecordBindings(shownResultsAxisRecords(), recordSettingsOwners(state.activeDrawing()));
  const sequenceSearch = createSequenceSearch({
    state,
    readAxisRecords: shownResultsAxisRecords,
    readRecordSpans: readShownRecordSpans,
    committedSets: shownCommittedSets,
    readRecordBindings: shownResultRecordBindings,
    runLosatJobs: runLosatJobsInOverlay,
    runHelper: runDiagramHelperOperation,
    helperOperations: { readQueries: DIAGRAM_HELPER_OPERATIONS.READ_SEARCH_QUERIES, buildHitSet: DIAGRAM_HELPER_OPERATIONS.SEARCH_HIT_ANNOTATION_SET },
    addSearchSet: (drawingId, set) => history.runUndoable(msg('Sequence search'),
      () => addSearchSetToDrawing(state.drawingOf(drawingId).slices.annotationSets, set)),
    shownRoot: shownResultRoot,
    closeRegionPopup,
    revealSection: () => settingsPaneActions.revealSettingsSection('sequence-search')
  });
  // Select region reads the axis of the shown Result from its track-slot geometry (Python's), and its popup is the
  // match frame's (the one Result popup).
  const regionSelect = createRegionSelect({
    readAxisRecords: shownAxisRecords,
    shownRoot: shownResultRoot,
    layoutEdit: layoutRepositionMode,
    openSelection: () => {
      const item = previewItemPopupOf(clickedPairwiseMatch.value)?.previewItem;
      return item?.kind === 'region' ? item : null;
    },
    openRegionPopup: (selection, element, anchor) => featureActions.openPreviewItemPopup({ element, item: selection }, anchor),
    closeRegionPopup,
    drawnFeatures: () => toRaw(extractedFeatures.value) ?? [],
    selectFeatures: (featureIds) => replaceSelectedFeatureIds(new Set(featureIds), { anchorId: featureIds[0] ?? '' }),
    // Binds as a search hit does; the popup is teleported out of the app root (no History adapter): its own step.
    addRegionAnnotation: async ({ recordKey, start, end }) => {
      const drawing = state.activeDrawing();
      await readRecordCatalog(drawing);
      if (drawing !== state.activeDrawing()) return null;
      return history.runUndoable(msg('Add to annotation set'), () => annotationEditor.addRegionAnnotation({ record: resultRecordBinding(shownResultRecordBindings(), recordKey), start, end }));
    },
    newDrawingFromRegion: regionDrawings.openWithRegion, findInOtherFiles: findInOtherFiles.open,
    // Python reads the bases from the committed request, the one the shown Result was made from.
    readRegionFasta: async ({ recordKey, start, end }) => (await readShownRecordSpans([{ recordKey, start, end }])).fasta,
    useAsSearchQuery: sequenceSearch.useSpanAsQuery,
    // While the tool is on, a touch that misses the tracks pans as a mouse drag does.
    panStart: startPan,
    panMove: doPan,
    panEnd: endPan
  });
  // Registered once sequenceSearch exists: a search's step names the sets it added (the live annotation-sets row).
  history.setSettingsStepListener((changes, direction, undrawn) => {
    requestLiveSettingApply(changes ? liveSettingRequestPaths(changes) : allLiveSettingRequestPaths(), mode.value,
      changes ? changedTrackHides(changes) : undefined, direction === 'undo', sequenceSearch.searchSetIdsOf(changes));
    refitSizedFigure(changes, direction, undrawn);
  });
  return {
    sequenceSearch,
    hostedBuild: readHostedBuild(document),
    fontSettings,
    defaultColorGroups,
    recordDisplayControls,
    featureRecordRotationDraft: featureRecordRotation.draft,
    recordActionsExpanded,
    toggleRecordActions,
    featureRecordRotationUnavailable,
    cancelRecordActions,
    setFeatureRecordRotationPosition: featureRecordRotation.setPosition,
    setFeatureRecordRotationReference: featureRecordRotation.setReference,
    setFeatureRecordRotationOffset: featureRecordRotation.setOffset,
    setFeatureRecordRotationOrientForward: featureRecordRotation.setOrientForward,
    applyFeatureRecordRotation: featureRecordRotation.apply,
    stageFeatureRecordRotation: featureRecordRotation.stage,
    closeFeaturePopup,
    featurePlacementActions: featureActions.placementActions,
    sessionSaveNeedsGenerate,
    processing,
    processingStatus,
    labelCapacityNotices,
    labelCapacityNoticeText,
    sessionImportPending,
    sessionSavePending,
    semanticMutationAvailable,
    sessionOperationAvailability,
    dialogChoicePending,
    missingBlastTsvDialog,
    missingBlastTsvRows,
    chooseLinearComparisonPairFile,
    comparisonCauseShown, showComparisonCause,
    runLinearComparisonIssueAction,
    resolveMissingBlastTsv,
    ...labelFitTemplate,
    labelCapacityText,
    labelFitNotice,
    shownInputSummary,
    shownInputSummaryStale,
    openLabelFitNotice,
    sessionSaveAvailable,
    sessionLoadAvailable,
    sessionBusyReason,
    generationCancelRequested,
    errorLog,
    errorDisplay,
    canRetryInteractiveSvgExport,
    reloadAfterOperationError,
    selectDrawing, drawingList, regionDrawings, findInOtherFiles,
    regionSelect,
    regionSelectOn: regionSelect.active,
    regionSelectAvailable: regionSelect.available,
    regionSelectUnavailableReason: regionSelect.unavailableReason,
    regionOverlay: regionSelect.overlay,
    regionOverlayMatrix: regionSelect.overlayMatrix,
    regionPopup: regionSelect.popup,
    drawingBar, drawingBarTabs, moveTabFocus, drawingSwitchAvailable,
    drawingSwitchBusyReason,
    canRetrySpecificRuleFailure: featureActions.canRetrySpecificRuleFailure,
    canEditSpecificRuleFailure: featureActions.canEditSpecificRuleFailure,
    canRetryLabelImportFailure: featureActions.canRetryLabelImportFailure,
    retryLabelImportFailure: undoableAction(msg('Load label edits'), featureActions.retryLabelImportFailure),
    editLabelImportFailure: featureActions.editLabelImportFailure,
    canRetryFeatureEditTableImport: featureActions.canRetryFeatureEditTableImport,
    retryFeatureEditTableImport: undoableAction(msg('Load feature edits'), featureActions.retryFeatureEditTableImport),
    reselectFeatureEditTable: featureActions.reselectFeatureEditTable,
    retrySpecificRuleFailure: featureActions.retrySpecificRuleFailure,
    editSpecificRuleFailure: featureActions.editSpecificRuleFailure,
    sessionTitle,
    sessionTitleLabel,
    results,
    // The Results with the displayed Result's edits written (R1); the test hooks read them.
    currentResults: previewRuntime.currentResults,
    selectedResultIndex,
    failedGeneratePreservedResult,
    generationFailureRecovery,
    collinearColorModeOptions: Object.entries(COLLINEAR_COLOR_MODE_TEXT),
    importedComparisonIntent: drawingMember('importedComparisonIntent'),
    importedComparisonNeedsResolution,
    importedComparisonActionLabels,
    importedComparisonActionText,
    importedComparisonCanInherit,
    importedComparisonCanReplace,
    inheritImportedComparison,
    replaceImportedComparison,
    clearImportedComparison,
    selectResult,
    resultPanelTab,
    lastRunInfo,
    annotationWarnings,
    featureIdentityNoticeSummary,
    removeUnmatchedFeatureEdits,
    comparisonWarnings,
    runInfoCopyStatus,
    exactReplayCopyStatus,
    svgContent,
    svgResultIdentity,
    zoom,
    layoutRepositionMode,
    isPanning,
    handleWheel,
    fitPreviewToViewport,
    canvasPan,
    canvasContainerRef,
    startPan,
    doPan,
    endPan,
    sidebarWidth,
    startResizing,
    settingsPane: settingsPaneActions.settingsPane,
    showPhonePane: settingsPaneActions.showPhonePane,
    revealSetting: settingsPaneActions.revealSetting,
    revealSettingsSection: settingsPaneActions.revealSettingsSection,
    closeSettingsPopup: settingsPaneActions.closeSettingsPopup,
    settingsPopupStyle,
    startSettingsPopupDrag: settingsPopupDrag.start,
    keepSettingsPopupSectionOpen,
    showSettingsPaneRail: settingsPaneActions.showSettingsPaneRail,
    showSettingsPaneStacked: settingsPaneActions.showSettingsPaneStacked,
    closeDockedSettings: settingsPaneActions.closeDockedSettings,
    mode,
    setCircularInputType,
    // The shown drawing's input type (each mode panel is the shown drawing's).
    inputType: drawingMember('inputType'),
    losatProgram: drawingMember('losatProgram'),
    circularConservation: drawingMember('circularConservation'),
    annotationSets: drawingMember('annotationSets'),
    selectedAnnotation,
    ...featureDensityTrackEditor,
    addAnnotationSet: annotationEditor.addAnnotationSet,
    renameAnnotationSet: annotationEditor.renameAnnotationSet,
    duplicateAnnotationSet: annotationEditor.duplicateAnnotationSet,
    removeAnnotationSet: annotationEditor.removeAnnotationSet,
    setAnnotationSetLegendLabel: annotationEditor.setAnnotationSetLegendLabel,
    addCoordinateAnnotation: annotationEditor.addCoordinateAnnotation,
    addSelectedFeatureAnnotations: annotationEditor.addSelectedFeatures,
    removeAnnotation: annotationEditor.removeAnnotation,
    setAnnotationTargetKind: annotationEditor.setAnnotationTargetKind,
    setAnnotationField: annotationEditor.setAnnotationField,
    annotationImportNotice,
    importAnnotationTableFile: undoableAction(msg('Import annotations'), annotationEditor.importAnnotationTableFile),
    renameAnnotation: annotationEditor.renameAnnotation,
    setAnnotationStyle: annotationEditor.setAnnotationStyle,
    canDownloadAnnotationTable: annotationEditor.canDownloadAnnotationTable,
    downloadAnnotationTable: annotationEditor.downloadAnnotationTable,
    annotationFeatureCaption: annotationEditor.featureTargetCaption,
    annotationRecordOptions: annotationEditor.recordOptionsFor,
    annotationRecordValue: annotationEditor.recordValueFor,
    setAnnotationRecord: annotationEditor.setRecordValue,
    annotationRecordRequired: annotationEditor.recordIsRequired,
    annotationRecordMissing: annotationEditor.recordIsMissing,
    annotationRecordMissingMessage: annotationEditor.recordMissingMessage,
    annotationRecordShowsSelect: annotationEditor.recordShowsSelect,
    annotationRecordCaption: annotationEditor.recordCaptionFor,
    circularConservationFastaInput,
    ...circularConservationSeriesPanel,
    addRingRowFiles, removeRingRowFile,
    addCircularConservationComparisonFile,
    removeCircularConservationSource,
    depthTrackRows,
    circularDepthTrackRows,
    linearDepthTrackRows,
    linearSourceDepthRows,
    linearSourceDepthSummary,
    linearDepthTrackCoverageLabel,
    linearDepthTrackIndexOptions,
    hasCircularDepthFiles,
    hasLinearDepthFiles,
    canShowDepthTrack,
    enabledOptionClass,
    depthToggleOptionClass,
    depthTrackCountLabel,
    depthSeriesDefaultName,
    depthSeriesName,
    normalizeDepthSlotTrackIndex,
    getDepthTrackLabel,
    setDepthTrackLabel,
    getDepthTrackColor,
    setDepthTrackColor,
    getDepthTrackLegendLabelForSlot,
    setDepthTrackLegendLabelForSlot,
    syncDepthTrackSlotLabel,
    addCircularDepthTrack,
    addLinearDepthTrack,
    removeCircularDepthTrack: askingRemoveCircularDepthTrack,
    removeLinearDepthTrack: askingRemoveLinearDepthTrack,
    getCircularDepthFile,
    setCircularDepthFile,
    getLinearDepthFile,
    setLinearDepthFile,
    setLinearSourceDepthFile,
    clearLinearSourceDepthFile,
    linearSeqs: drawingMember('linearSeqs'),
    // The shown drawing's input fields (its five card files among them) and
    // Circular record list.
    files: drawingMember('files'),
    circularRecordList: drawingMember('circularRecordList'),
    linearSourceGroups,
    linearCardFileLine,
    linearCardHoldsInput,
    isLinearSourceFolded,
    syncLinearSourceFold,
    linearSourceInspectionFailed,
    linearRecordGencodeProblem,
    linearSourceGencodeProblem,
    linearGencodeMark,
    linearRecordLayoutEnabled: drawingMember('linearRecordLayoutEnabled'),
    linearRecordGap: drawingMember('linearRecordGap'),
    linearRecordRows: drawingMember('linearRecordRows'),
    linearComparisonPlan: drawingMember('linearComparisonPlan'),
    linearComparisonResolution: drawingMember('linearComparisonResolution'),
    linearComparisonGlobalAction,
    linearComparisonUi,
    canRunLinearLosat,
    hasLinearComparisonIntent: drawingMember('hasLinearComparisonIntent'),
    hasActiveLinearLosatIntent: drawingMember('hasActiveLinearLosatIntent'),
    hasActiveLinearUploadIntent: drawingMember('hasActiveLinearUploadIntent'),
    linearComparisonTimeline,
    linearComparisonRecordLabel,
    linearComparisonRecordName,
    linearLosatCacheInfoByEdgeKey, linearLosatWindowPairs, linearWindowPairSide,
    linearLayoutTokens,
    syncLinearRecordLayout,
    setLinearRecordLayoutEnabled,
    setLinearRecordRow,
    moveLinearRecordWithinRow,
    setLinearComparisonGlobalAction,
    setLinearComparisonLosatMode,
    setLinearComparisonLosatpMode,
    setLinearComparisonGapAction,
    addLinearComparison,
    omitLinearComparison,
    clearSelectedLinearComparisons,
    setLinearComparisonEndpoint,
    setLinearComparisonSource,
    setLinearComparisonFile,
    setLinearComparisonCardFile,
    reuseLinearComparisonFile,
    deactivateLinearComparisonFile,
    setLinearComparisonLosatFilename,
    reuseLinearComparisonLosatFilename,
    deactivateLinearComparisonLosatFilename,
    setResolvedLinearComparisonLosatFilename,
    reuseResolvedLinearComparisonLosatFilename,
    deactivateResolvedLinearComparisonLosatFilename,
    setWindowPairLosatFilename,
    addLinearComparisonBatch,
    linearRecordRowFor: shownLinearRecordRowFor,
    linearReorderNotice,
    addLinearSeq,
    removeLastLinearSeq,
    requestLinearSourceRemoval,
    applyLinearSourceRemoval,
    cancelLinearSourceRemoval,
    trapDialogFocus,
    linearSourceRemovalDialog,
    linearSourceRemovalTarget,
    linearSourceRemovalTargetName,
    linearSourceRemovalCanDelete,
    setLinearSeqPrimaryFile,
    setLinearInputType,
    setLinearRecordSelector,
    setLinearCardWindow, addLinearCardWindow, removeLinearCardWindow,
    linearSourceMoveBlockedReason,
    canMoveLinearSource,
    moveLinearSource,
    linearFileCardRecord,
    getLinearSourceDefaultDefinition,
    setLinearSourceDefaultDefinition,
    getLinearSourceDefaultSubtitle,
    setLinearSourceDefaultSubtitle,
    resolveLinearRecordEffectiveDefinition,
    resolveLinearRecordEffectiveSubtitle,
    resetLinearRecordDefinition,
    resetLinearRecordSubtitle,
    // The template lists the shown drawing's cards; its selector readers read
    // that drawing (OV-2149: a card's uid names it only inside its drawing).
    linearRecordOptions: readingShownLinearCard(linearRecordSelector.optionsFor),
    // The Inspect records button and the test API read the shown drawing's cards.
    refreshLinearRecordSelectors: () => linearRecordSelector.refresh(state.activeDrawing()),
    linearRecordSelectorDisabled: readingShownLinearCard(linearRecordSelector.isDisabled),
    linearRecordSelectorError: readingShownLinearCard(linearRecordSelector.errorFor),
    linearRecordSelectorWarning: readingShownLinearCard(linearRecordSelector.warningFor),
    form: drawingMember('form'),
    adv: drawingMember('adv'),
    linearTypographyLinked: drawingMember('linearTypographyLinked'),
    setLinearTypographyLinked: linearTypography.setLinked,
    setLinearRulerLabelFontSize: linearTypography.setRulerLabelFontSize,
    setLinearScaleFontSize: linearTypography.setScaleFontSize,
    unmanagedConfigOverrideEntries,
    resetUnmanagedConfigOverride,
    comparisonProfileDefault,
    recordAutoOrientation,
    identityPercentMaximum: IDENTITY_PERCENT_MAXIMUM,
    /** @param {0 | 1} end @param {string} text */
    setIdentityScaleEnd: (end, text) => {
      const drawing = state.activeDrawing();
      drawing.adv.identity_scale = identityScaleWithEnd(drawing.adv.identity_scale, end, text);
    },
    comparisonHeightValidationError,
    optionalNumberInputValue,
    setOptionalNumberInputValue,
    autoValueText: autoValueDisplay.autoValueText,
    autoValueVisible: autoValueDisplay.autoValueVisible,
    canUseLinearRulerOnAxis,
    ...linearWrapField,
    ...figureSizeField,
    canUseCircularScaleStyling,
    circularTrackNewRenderer,
    linearTrackNewRenderer,
    circularTrackSlotsPanelOpen,
    toggleCircularTrackSlotsPanel,
    linearTrackSlotsPanelOpen,
    toggleLinearTrackSlotsPanel,
    circularTrackRenderers: circularTrackSlotEditor.circularTrackRenderers,
    circularTrackSlotEditorKey: circularTrackSlotEditor.circularTrackSlotEditorKey,
    circularTrackSlotOptionsOpen: circularTrackSlotEditor.circularTrackSlotOptionsOpen,
    setCircularTrackSlotOptionsOpen: circularTrackSlotEditor.setCircularTrackSlotOptionsOpen,
    circularTrackSlotOptionsSummary: circularTrackSlotEditor.circularTrackSlotOptionsSummary,
    circularDrawingDinucleotide: circularTrackSlotEditor.circularDrawingDinucleotide,
    updateCircularTrackSlotMeasure: circularTrackSlotEditor.updateCircularTrackSlotMeasure,
    updateCircularTrackSlotGap: circularTrackSlotEditor.updateCircularTrackSlotGap,
    circularTrackRendererLabel: circularTrackSlotEditor.circularTrackRendererLabel,
    resetCircularTrackSlotsFromSimpleControls: circularTrackSlotEditor.resetCircularTrackSlotsFromSimpleControls,
    resetCircularTrackSlotsToPreset: circularTrackSlotEditor.resetCircularTrackSlotsToPreset,
    setCircularTrackSlotsEnabled: circularTrackSlotEditor.setCircularTrackSlotsEnabled,
    setCircularGcSuppressed: circularTrackSlotEditor.setCircularGcSuppressed,
    setTrackShown,
    setCircularSkewSuppressed: circularTrackSlotEditor.setCircularSkewSuppressed,
    setCircularDinucleotide: circularTrackSlotEditor.setCircularDinucleotide,
    setCircularTrackSlotNt: circularTrackSlotEditor.setCircularTrackSlotNt,
    addCircularTrackSlot: askingAnnotationTrackAdd(askingStackPorts('circular_track_slots'), circularTrackSlotEditor.addCircularTrackSlot),
    circularTrackAddBlockReason: circularTrackSlotEditor.circularTrackAddBlockReason,
    duplicateCircularTrackSlot: circularTrackSlotEditor.duplicateCircularTrackSlot,
    canDuplicateCircularTrackSlot: circularTrackSlotEditor.canDuplicateCircularTrackSlot,
    removeCircularTrackSlot: askingRowRemoval(circularTrackSlotEditor.removeCircularTrackSlot, 'circular_track_slots'),
    setCircularTrackSlotEnabled: askingRowSwitch(circularTrackSlotEditor.setCircularTrackSlotEnabled, 'circular_track_slots'),
    circularTrackSlotEffectiveEnabled: circularTrackSlotEditor.circularTrackSlotEffectiveEnabled,
    circularTrackSlotHiddenBySuppress: circularTrackSlotEditor.circularTrackSlotHiddenBySuppress,
    circularTrackSlotSuppressMessage: circularTrackSlotEditor.circularTrackSlotSuppressMessage,
    moveCircularTrackSlot: circularTrackSlotEditor.moveCircularTrackSlot,
    canMoveCircularTrackSlot: circularTrackSlotEditor.canMoveCircularTrackSlot,
    moveCircularTrackSlotOutside: circularTrackSlotEditor.moveCircularTrackSlotOutside,
    moveCircularTrackSlotInside: circularTrackSlotEditor.moveCircularTrackSlotInside,
    moveCircularTrackSlotToAxis: circularTrackSlotEditor.moveCircularTrackSlotToAxis,
    canMoveCircularTrackSlotOutside: circularTrackSlotEditor.canMoveCircularTrackSlotOutside,
    canMoveCircularTrackSlotInside: circularTrackSlotEditor.canMoveCircularTrackSlotInside,
    canMoveCircularTrackSlotToAxis: circularTrackSlotEditor.canMoveCircularTrackSlotToAxis,
    updateCircularTrackSlotRenderer: circularTrackSlotEditor.updateCircularTrackSlotRenderer,
    updateCircularTrackSlotPlacement: circularTrackSlotEditor.updateCircularTrackSlotPlacement,
    updateCircularTrackFeatureLane: circularTrackSlotEditor.updateCircularTrackFeatureLane,
    circularTrackSlotIssue: circularTrackSlotEditor.circularTrackSlotIssue,
    circularTrackGlobalIssues: circularTrackSlotEditor.circularTrackGlobalIssues,
    circularAnnotationAnchorOptions: circularTrackSlotEditor.circularAnnotationAnchorOptions,
    circularAnnotationAnchorIsKnown: circularTrackSlotEditor.circularAnnotationAnchorIsKnown,
    annotationTrackMarkOptions: circularTrackSlotEditor.annotationTrackMarkOptions,
    circularAnnotationMarkSelected: circularTrackSlotEditor.circularAnnotationMarkSelected,
    setCircularAnnotationMarkSelected: circularTrackSlotEditor.setCircularAnnotationMarkSelected,
    circularAnnotationLaneGapValue: circularTrackSlotEditor.circularAnnotationLaneGapValue,
    setCircularAnnotationLaneGap: circularTrackSlotEditor.setCircularAnnotationLaneGap,
    circularAnnotationPaddingValue: circularTrackSlotEditor.circularAnnotationPaddingValue,
    setCircularAnnotationPadding: circularTrackSlotEditor.setCircularAnnotationPadding,
    circularAnnotationCoverAnchor: circularTrackSlotEditor.circularAnnotationCoverAnchor,
    setCircularAnnotationCoverAnchor: circularTrackSlotEditor.setCircularAnnotationCoverAnchor,
    setCircularOverlayAnnotationParam: circularTrackSlotEditor.setCircularOverlayAnnotationParam,
    supportsCircularTrackSlotPlacement: circularTrackSlotEditor.supportsCircularTrackSlotPlacement,
    circularTrackSlots: circularTrackSlotEditor.circularTrackSlots,
    circularTrackStackEntries: circularTrackSlotEditor.circularTrackStackEntries,
    circularTrackSlotCliSpec: circularTrackSlotEditor.circularTrackSlotCliSpec,
    circularTrackSlotDisplayLabel: circularTrackSlotEditor.circularTrackSlotDisplayLabel,
    circularTrackSlotDisplayMeta: circularTrackSlotEditor.circularTrackSlotDisplayMeta,
    circularTrackSlotLegendLabelPlaceholder: circularTrackSlotEditor.circularTrackSlotLegendLabelPlaceholder,
    circularTrackSlotColor: circularTrackSlotEditor.circularTrackSlotColor,
    circularTrackSlotHasSkewColorOverride: circularTrackSlotEditor.circularTrackSlotHasSkewColorOverride,
    circularTrackSlotSkewColorValue: circularTrackSlotEditor.circularTrackSlotSkewColorValue,
    circularTrackSlotGeometryAutoText: circularTrackSlotEditor.circularTrackSlotGeometryAutoText,
    circularTrackSlotGeometryHasManual: circularTrackSlotEditor.circularTrackSlotGeometryHasManual,
    circularTrackSlotGeometryUnitSuffix: circularTrackSlotEditor.circularTrackSlotGeometryUnitSuffix,
    setCircularTrackSlotSkewColor: circularTrackSlotEditor.setCircularTrackSlotSkewColor,
    clearCircularTrackSlotSkewColor: circularTrackSlotEditor.clearCircularTrackSlotSkewColor,
    isManagedCircularConservationSlot: circularTrackSlotEditor.isManagedCircularConservationSlot,
    circularTrackPresetSummary: circularTrackSlotEditor.circularTrackPresetSummary,
    circularTrackSlotUsesPresetGeometry: circularTrackSlotEditor.circularTrackSlotUsesPresetGeometry,
    linearTrackRenderers: linearTrackSlotEditor.linearTrackRenderers,
    linearTrackSlotEditorKey: linearTrackSlotEditor.linearTrackSlotEditorKey,
    linearTrackRendererLabel: linearTrackSlotEditor.linearTrackRendererLabel,
    resetLinearTrackSlotsFromSimpleControls: linearTrackSlotEditor.resetLinearTrackSlotsFromSimpleControls,
    setLinearTrackSlotsEnabled: linearTrackSlotEditor.setLinearTrackSlotsEnabled,
    addLinearTrackSlot: askingAnnotationTrackAdd(askingStackPorts('linear_track_slots'), linearTrackSlotEditor.addLinearTrackSlot),
    linearTrackAddBlockReason: linearTrackSlotEditor.linearTrackAddBlockReason,
    duplicateLinearTrackSlot: linearTrackSlotEditor.duplicateLinearTrackSlot,
    canDuplicateLinearTrackSlot: linearTrackSlotEditor.canDuplicateLinearTrackSlot,
    removeLinearTrackSlot: askingRowRemoval(linearTrackSlotEditor.removeLinearTrackSlot, 'linear_track_slots'),
    setLinearTrackSlotEnabled: askingRowSwitch(linearTrackSlotEditor.setLinearTrackSlotEnabled, 'linear_track_slots'),
    moveLinearTrackSlot: linearTrackSlotEditor.moveLinearTrackSlot,
    canMoveLinearTrackSlot: linearTrackSlotEditor.canMoveLinearTrackSlot,
    moveLinearTrackSlotAbove: linearTrackSlotEditor.moveLinearTrackSlotAbove,
    moveLinearTrackSlotBelow: linearTrackSlotEditor.moveLinearTrackSlotBelow,
    moveLinearTrackSlotToAxis: linearTrackSlotEditor.moveLinearTrackSlotToAxis,
    canMoveLinearTrackSlotAbove: linearTrackSlotEditor.canMoveLinearTrackSlotAbove,
    canMoveLinearTrackSlotBelow: linearTrackSlotEditor.canMoveLinearTrackSlotBelow,
    canMoveLinearTrackSlotToAxis: linearTrackSlotEditor.canMoveLinearTrackSlotToAxis,
    updateLinearTrackSlotRenderer: linearTrackSlotEditor.updateLinearTrackSlotRenderer,
    updateLinearTrackSlotPlacement: linearTrackSlotEditor.updateLinearTrackSlotPlacement,
    linearTrackSlotIssue: linearTrackSlotEditor.linearTrackSlotIssue,
    linearTrackGlobalIssues: linearTrackSlotEditor.linearTrackGlobalIssues,
    linearAnnotationAnchorOptions: linearTrackSlotEditor.linearAnnotationAnchorOptions,
    linearAnnotationAnchorIsKnown: linearTrackSlotEditor.linearAnnotationAnchorIsKnown,
    linearAnnotationMarkSelected: linearTrackSlotEditor.linearAnnotationMarkSelected,
    setLinearAnnotationMarkSelected: linearTrackSlotEditor.setLinearAnnotationMarkSelected,
    linearAnnotationLaneGapValue: linearTrackSlotEditor.linearAnnotationLaneGapValue,
    setLinearAnnotationLaneGap: linearTrackSlotEditor.setLinearAnnotationLaneGap,
    linearAnnotationPaddingValue: linearTrackSlotEditor.linearAnnotationPaddingValue,
    setLinearAnnotationPadding: linearTrackSlotEditor.setLinearAnnotationPadding,
    linearAnnotationCoverAnchor: linearTrackSlotEditor.linearAnnotationCoverAnchor,
    setLinearAnnotationCoverAnchor: linearTrackSlotEditor.setLinearAnnotationCoverAnchor,
    setLinearOverlayAnnotationParam: linearTrackSlotEditor.setLinearOverlayAnnotationParam,
    linearTrackSlotHeightValue: linearTrackSlotEditor.linearTrackSlotHeightValue,
    linearTrackSlotGeometryAutoText: linearTrackSlotEditor.linearTrackSlotGeometryAutoText,
    linearTrackSlotGeometryHasManual: linearTrackSlotEditor.linearTrackSlotGeometryHasManual,
    linearTrackSlotGeometryUnitSuffix: linearTrackSlotEditor.linearTrackSlotGeometryUnitSuffix,
    setLinearTrackSlotHeight: linearTrackSlotEditor.setLinearTrackSlotHeight,
    linearTrackSlotHasSkewColorOverride: linearTrackSlotEditor.linearTrackSlotHasSkewColorOverride,
    linearTrackSlotSkewColorValue: linearTrackSlotEditor.linearTrackSlotSkewColorValue,
    setLinearTrackSlotSkewColor: linearTrackSlotEditor.setLinearTrackSlotSkewColor,
    clearLinearTrackSlotSkewColor: linearTrackSlotEditor.clearLinearTrackSlotSkewColor,
    syncLinearDepthSlotHeightsFromDepthTracks: linearTrackSlotEditor.syncLinearDepthSlotHeightsFromDepthTracks,
    linearTrackSlots: linearTrackSlotEditor.linearTrackSlots,
    linearTrackStackEntries: linearTrackSlotEditor.linearTrackStackEntries,
    linearTrackSlotCliSpec: linearTrackSlotEditor.linearTrackSlotCliSpec,
    linearTrackSlotDisplayLabel: linearTrackSlotEditor.linearTrackSlotDisplayLabel,
    linearTrackSlotDisplayMeta: linearTrackSlotEditor.linearTrackSlotDisplayMeta,
    linearTrackSlotLegendLabelPlaceholder: linearTrackSlotEditor.linearTrackSlotLegendLabelPlaceholder,
    linearDrawingDinucleotide: linearTrackSlotEditor.linearDrawingDinucleotide,
    linearTrackSlotPlacementLabel: linearTrackSlotEditor.linearTrackSlotPlacementLabel,
    linearTrackSlotUsesPresetGeometry: linearTrackSlotEditor.linearTrackSlotUsesPresetGeometry,
    losat: drawingMember('losat'),
    losatExecution,
    richFeaturePopup,
    markRichFeaturePopupChanged,
    ...losatSettings,
    losatCacheInfo,
    losatThreadingStatus,
    orthogroups,
    featureOrthogroupIndex,
    // Similarity groups are the shown drawing's.
    orthogroupNameOverrides: computed(() => state.activeDrawing().orthogroupNameOverrides),
    orthogroupDescriptionOverrides: computed(() => state.activeDrawing().orthogroupDescriptionOverrides),
    selectedOrthogroupId,
    orthogroupSearch,
    orthogroupSortMode,
    showRightDrawer,
    rightDrawerTab,
    orthogroupCount: orthogroupActions.orthogroupCount,
    filteredOrthogroups: orthogroupActions.filteredOrthogroups,
    selectedOrthogroup: orthogroupActions.selectedOrthogroup,
    selectedOrthogroupDetail: orthogroupActions.selectedOrthogroupDetail,
    resolveOrthogroupName: (/** @type {import('./orthogroups.js').GroupOrId} */ groupOrId) => (
      orthogroupActions.resolveOrthogroupName(groupOrId, shownDrawingOfMode('linear'))
    ),
    orthogroupDescriptionFieldText: orthogroupActions.orthogroupDescriptionFieldText,
    orthogroupScope: orthogroupActions.orthogroupScope,
    orthogroupScopeLabel: orthogroupActions.orthogroupScopeLabel,
    orthogroupRows: orthogroupActions.orthogroupRows,
    isOrthogroupRenamed: orthogroupActions.isOrthogroupRenamed,
    hasOrthogroupMemberSequence: orthogroupActions.hasOrthogroupMemberSequence,
    orthogroupCopyFeedbackLabel: orthogroupActions.orthogroupCopyFeedbackLabel, orthogroupMemberRead: orthogroupActions.orthogroupMemberRead,
    orthogroupMemberRangeText: orthogroupActions.memberRangeText,
    copyOrthogroupSequences: orthogroupActions.copyOrthogroupSequences,
    downloadOrthogroupSequences: orthogroupActions.downloadOrthogroupSequences,
    copyOrthogroupMemberSequence: orthogroupActions.copyOrthogroupMemberSequence,
    downloadOrthogroupMemberSequence: orthogroupActions.downloadOrthogroupMemberSequence,
    selectOrthogroup: orthogroupActions.selectOrthogroup,
    setOrthogroupNameOverride: orthogroupActions.setOrthogroupNameOverride,
    setOrthogroupDescriptionOverride: orthogroupActions.setOrthogroupDescriptionOverride,
    resetOrthogroupRename: orthogroupActions.resetOrthogroupRename,
    orthogroupDormantNames: orthogroupActions.orthogroupDormantNames,
    clearOrthogroupDormantOverrides: orthogroupActions.clearOrthogroupDormantOverrides,
    highlightOrthogroupById: orthogroupActions.highlightOrthogroupById,
    similarityAlignmentDraft: similarityAlignmentActions.draft,
    featureAnchorSelector: featureAnchor.selector,
    featureAnchorBusy: featureAnchor.busy,
    featureAnchorPreview: featureAnchor.preview,
    featureAnchorPreviewLines: featureAnchor.previewLines,
    featureAnchorError: featureAnchor.error,
    featureAnchorNotice: featureAnchor.notice,
    resolveFeatureAnchor: featureAnchor.resolve,
    applyFeatureAnchor: featureAnchor.apply,
    cancelFeatureAnchor: featureAnchor.cancel,
    recordPositionResetBusy: recordPositionReset.busy,
    resetRecordPositions: recordPositionReset.reset,
    similarityAlignmentStatus: similarityAlignmentActions.status,
    similarityAlignmentBusy: similarityAlignmentActions.busy,
    similarityAlignmentError: similarityAlignmentActions.error,
    similarityAlignmentSummary: similarityAlignmentActions.summary,
    similarityAlignmentNotice: similarityAlignmentActions.notice,
    similarityAlignmentRepair: similarityAlignmentActions.repair,
    similarityAlignmentDialogOpen: similarityAlignmentActions.dialogOpen,
    previewPopupsHeld: previewPopupsHeldNow,
    similarityAlignmentApplyDisabledReason: similarityAlignmentActions.applyDisabledReason,
    similarityAlignmentDrawerReferenceKey: similarityAlignmentActions.drawerReferenceKey,
    featureDetailsPending,
    similarityAlignmentPlanInspector: similarityAlignmentActions.activePlanInspector,
    canApplySimilarityAlignment: similarityAlignmentActions.canApply,
    resetSimilarityAlignment: similarityAlignmentActions.resetAlignment,
    reselectSimilarityAlignmentReference,
    startSimilarityAlignmentFromPopup: similarityAlignmentActions.startFromPopup,
    startSimilarityAlignmentFromDrawer,
    similarityAlignmentDrawerReferenceOptions: similarityAlignmentActions.drawerReferenceOptions,
    similarityAlignmentDrawer: similarityAlignmentActions.drawer,
    setSimilarityAlignmentDrawerReference: similarityAlignmentActions.setDrawerReference,
    selectSimilarityAlignmentCandidate: similarityAlignmentActions.selectCandidate,
    skipSimilarityAlignmentRecord: similarityAlignmentActions.skipRecord,
    similarityAlignmentDirectionPreview: similarityAlignmentActions.directionPreview,
    setSimilarityAlignmentDirectionMode: similarityAlignmentActions.setDirectionMode,
    setSimilarityAlignmentCustomDirection: similarityAlignmentActions.setCustomDirection,
    similarityAlignmentCustomDirection: similarityAlignmentActions.customDirection,
    similarityAlignmentResetPreview: similarityAlignmentActions.resetPreview,
    similarityAlignmentResetDialogOpen: similarityAlignmentActions.resetDialogOpen,
    similarityAlignmentResetScope: similarityAlignmentActions.resetScope,
    openSimilarityAlignmentReset,
    cancelSimilarityAlignmentReset: similarityAlignmentActions.cancelReset,
    applySimilarityAlignmentReset: () => applySimilarityAlignmentDialog(true),
    applySimilarityAlignmentDraft: () => applySimilarityAlignmentDialog(),
    cancelSimilarityAlignmentDraft: cancelSimilarityAlignmentDialog,
    cancelSimilarityAlignmentDialog,
    similarityAlignmentReview,
    similarityAlignmentEditorDisabledReason,
    similarityAlignmentCanvasHover,
    previewSimilarityAlignmentCandidate: similarityAlignmentActions.previewCandidate,
    clearSimilarityAlignmentCandidatePreview: similarityAlignmentActions.clearCandidatePreview,
    isRightDrawerTabAvailable: rightDrawerActions.isRightDrawerTabAvailable,
    openRightDrawerTab: rightDrawerActions.openRightDrawerTab,
    toggleRightDrawer: rightDrawerActions.toggleRightDrawer,
    closeRightDrawer: rightDrawerActions.closeRightDrawer,
    openOrthogroupInDrawer,
    refreshCircularRecordOrder,
    importAuxiliaryTable, waitForAuxiliaryFileImport: auxiliaryFileImports.wait, canRetryAuxiliaryImportFailure,
    retryAuxiliaryImportFailure: () => history.runUndoableCheckpoint(msg('Change uploaded file'), retryAuxiliaryImportFailure, { shouldCommit: result => result !== false }),
    circularRecordPresentationOptions,
    circularRecordPresentationError,
    circularSingleRecordPresentationEnabled,
    circularRecordDiscoveryState,
    circularRecordInspectionEnabled,
    inspectCircularSourceRecords,
    circularRecordSelectionEnabled,
    showCircularCanvasSetting,
    showCircularTicksSlot,
    waitForCircularRecordScrollAnchor,
    DEFAULT_LINEAR_OUTPUT_PREFIX,
    setCircularRecordPresentationSelector,
    chooseCircularRecordLayout, circularRecordLayout, ONE_CIRCLE_LIMITS,
    setCircularSourceFile, addCircularSourceFiles, removeCircularSourceCard, canMoveCircularSource, moveCircularSourceCard, circularFileCards,
    replaceInputFile: inputReplacement.replaceInputFile,
    uploadCircularSourceFile,
    uploadCircularDepthFile,
    uploadLinearSeqPrimaryFile,
    shownInputComplete,
    uploadLinearDepthFile,
    uploadLinearSourceDepthFile,
    uploadLinearComparisonCardFile,
    recordSelection,
    circularDrawnRecordCount,
    formatRecordLength,
    paletteDefinitions,
    paletteNames,
    selectedPalette: drawingMember('selectedPalette'),
    currentColors: drawingMember('currentColors'),
    editDefaultColor: resultsManager.editDefaultColor,
    // Read-only here: the toggle's transition is its one writer.
    paletteInstantPreviewEnabled: computed(() => paletteInstantPreviewEnabled.value),
    setPaletteInstantPreview,
    appliedPaletteName,
    appliedPaletteColors,
    pendingPaletteName: drawingMember('pendingPaletteName'),
    pendingPaletteColors: drawingMember('pendingPaletteColors'),
    hasPendingPaletteDraft: drawingMember('hasPendingPaletteDraft'),
    requestPaletteChange,
    selectPalette,
    requestResetColors,
    losatResultRows, clearableLosatResultCount: computed(clearableLosatResultCount),
    downloadLosatCache,
    downloadLosatPair,
    clearLosatCache,
    getLosatPairDefaultName,
    readShownRecords,
    getCircularRecordOrderLabel,
    getCircularRecordRow,
    getCircularRecordRowOptions,
    setCircularRecordRow,
    canMoveCircularRecordOrderUp,
    canMoveCircularRecordOrderDown,
    moveCircularRecordOrderUp,
    moveCircularRecordOrderDown,
    resetCircularRecordOrder,
    filterMode: drawingMember('filterMode'),
    manualBlacklist: drawingMember('manualBlacklist'),
    manualWhitelist: drawingMember('manualWhitelist'),
    setLabelFilterMode,
    setBlacklistText,
    addWhitelistRule,
    setWhitelistRuleField,
    removeWhitelistRule,
    removePriorityRule,
    featureKeys,
    defaultColorKeys,
    newColorFeat,
    newColorVal,
    addCustomColor,
    newFeatureToAdd,
    addFeature,
    removeFeature,
    getFeatureShape,
    setFeatureShape,
    manualSpecificRules: drawingMember('manualSpecificRules'),
    ruleMatchingPending,
    newSpecRule,
    specificRuleHasQualifier,
    specificRulePresets,
    specificRuleQualifierSuggestions,
    selectedSpecificPreset,
    specificRulePresetLoading,
    addSpecificRule,
    applySpecificRulePreset,
    clearAllSpecificRules,
    downloadSpecificRulesTsv,
    moveSpecificRuleDown,
    moveSpecificRuleUp,
    removeSpecificRule,
    setSpecificRuleField,
    specificRulePattern, specificRulePatternDraft, specificRulePatternFieldId,
    editSpecificRulePattern, retrySpecificRulePattern, revertSpecificRulePattern,
    extractedFeatures,
    featureEditorStatus,
    featureEditorStatusText,
    featureExtractionPending,
    featureRecordIds,
    selectedFeatureRecordIdx,
    featurePanelTab,
    featureSearchInput,
    featureSearch,
    previewFeatureSearchInput,
    previewFeatureSearchQuery,
    previewFeatureSearchField,
    previewFeatureSearchQualifierKey,
    previewFeatureSearchUseRegex,
    previewFeatureSearchMatches,
    previewFeatureSearchMatchDetails,
    previewFeatureSearchActiveIndex,
    previewFeatureSearchError,
    previewFeatureSearchRenderedCount,
    previewFeatureSearchFieldOptions: previewFeatureSearch.previewFeatureSearchFieldOptions,
    previewFeatureSearchQualifierEnabled: previewFeatureSearch.previewFeatureSearchQualifierEnabled,
    previewFeatureSearchHasMatches: previewFeatureSearch.previewFeatureSearchHasMatches,
    previewFeatureSearchCanOpenActive: previewFeatureSearch.previewFeatureSearchCanOpenActive,
    previewFeatureSearchCanSearch: previewFeatureSearch.previewFeatureSearchCanSearch,
    previewFeatureSearchStatusText: previewFeatureSearch.previewFeatureSearchStatusText,
    previewFeatureSearchActiveDetail: previewFeatureSearch.previewFeatureSearchActiveDetail,
    applyPreviewFeatureSearch: previewFeatureSearch.applySearch,
    goToNextPreviewFeatureSearchMatch: previewFeatureSearch.goToNext,
    goToPreviousPreviewFeatureSearchMatch: previewFeatureSearch.goToPrevious,
    clearPreviewFeatureSearch: previewFeatureSearch.clearSearch,
    openPreviewFeatureSearchActiveMatch: previewFeatureSearch.openActiveMatch,
    previewFeatureSearchVisible: state.previewFeatureSearchVisible,
    previewFeatureSearchHiddenMatchCount: previewFeatureSearch.previewFeatureSearchHiddenMatchCount,
    previewFeatureSearchSequenceFailure: previewFeatureSearch.previewFeatureSearchSequenceFailure,
    hidePreviewFeatureSearchBar: previewFeatureSearch.hideSearchBar,
    togglePreviewFeatureSearchBar: previewFeatureSearch.toggleSearchBar,
    selectedFeatureIds,
    selectedFeatureAnchorId,
    featureSelectionStatus,
    featureSelectionDrag,
    selectedFeatureCount,
    selectedFeatures,
    hasFeatureSelection,
    featureSelectionMarqueeStyle, featureSelectionToolbarStyle, startFeatureSelectionToolbarDrag: startToolbarDrag, clearFeatureSelection,
    openFirstSelectedFeature,
    selectedFeatureBulkColor,
    selectedFeatureBulkCaption,
    selectedFeatureBulkVisibility,
    selectedFeatureBulkStrokeColor,
    selectedFeatureBulkStrokeWidth,
    applySelectedFeatureColor,
    applySelectedFeatureVisibility,
    applySelectedFeatureStroke,
    visibleFeatureRows,
    featureRecordPickerVisible,
    featureListTopSpacerPx,
    featureListExpandedId,
    featureListBottomSpacerPx,
    isFeatureDrawerMounted,
    featureListScrollRef,
    handleFeatureListScroll,
    labelSearch,
    editableLabels,
    filteredEditableLabels,
    labelTextBulkOverrides: drawingMember('labelTextBulkOverrides'),
    autoLabelReflowEnabled,
    labelReflowProcessing,
    labelReflowLastError,
    filteredFeatures,
    featureListState,
    featureColorOverrides: drawingMember('featureColorOverrides'),
    featureVisibilityManualRules: drawingMember('featureVisibilityManualRules'),
    featureVisibilityRules: drawingMember('featureVisibilityRules'),
    featureOverrides: drawingMember('featureOverrides'),
    labelSettingsShown,
    featureStrokeOverrides: drawingMember('featureStrokeOverrides'),
    addFeatureVisibilityRule: addFeatureVisibilityRuleWithHistory,
    downloadFeatureVisibilityRulesTsv,
    featureVisibilityFeatureSuggestions,
    featureVisibilityQualifierSuggestions,
    featureVisibilityRuleDetail,
    getFeatureColor,
    getFeatureColorValue,
    moveFeatureVisibilityRuleDown: moveFeatureVisibilityRuleDownWithHistory,
    moveFeatureVisibilityRuleUp: moveFeatureVisibilityRuleUpWithHistory,
    removeFeatureVisibilityRule: removeFeatureVisibilityRuleWithHistory,
    setFeatureVisibility: setFeatureVisibilityWithHistory,
    setFeatureVisibilityRuleField: setFeatureVisibilityRuleFieldWithHistory,
    requestFeatureColorChange: requestFeatureColorChangeWithHistory,
    setFeatureColorValue: setFeatureColorValueWithHistory,
    setFeatureColor: setFeatureColorWithHistory,
    getEditableLabelByFeatureId,
    svgContainer,
    clickedFeature,
    clickedFeaturePos,
    clickedPairwiseMatch,
    clickedPairwiseMatchPos,
    featurePopupRef,
    featurePopupStyle,
    startFeaturePopupDrag: featurePopupDragging.start,
    startFeaturePopupResize,
    pairwiseMatchPopupRef,
    pairwiseMatchPopupStyle,
    startPairwiseMatchPopupDrag: pairwiseMatchPopupDragging.start,
    startPairwiseMatchPopupResize,
    selectedPairwiseBlockOrthogroupId,
    pairwiseGroupView,
    selectPairwiseBlockOrthogroup,
    openPairwiseFeatureRow,
    clickedFeatureSummary,
    clickedFeatureHeading,
    copyClickedFeatureDetails,
    copyClickedFeatureQualifiers,
      copyText: copyTextToClipboard,
    downloadText, sequenceReads,
    canUseClickedOrthogroupActions,
    clickedOrthogroupDetail,
    alignByClickedOrthogroup,
    highlightClickedOrthogroup,
    clearOrthogroupHighlight,
    openClickedOrthogroupInEditor,
    popupLegendGroups,
    updateClickedFeatureColor: updateClickedFeatureColorWithHistory,
    updateClickedFeatureVisibility: updateClickedFeatureVisibilityWithHistory,
    handleFeatureVisibilityScopeChoice: handleFeatureVisibilityScopeChoiceWithHistory,
    handleLegendNameCommit: handleLegendNameCommitWithHistory,
    joinLegendGroup: joinLegendGroupWithHistory,
    resetClickedFeatureFillColor: resetClickedFeatureFillColorWithHistory,
    updateClickedFeatureStroke: updateClickedFeatureStrokeWithHistory,
    getFeatureStrokeColorValue,
    setClickedFeatureStrokeColorValue: setClickedFeatureStrokeColorValueWithHistory,
    setClickedFeatureStrokeWidthValue: setClickedFeatureStrokeWidthValueWithHistory,
    resetClickedFeatureStroke: resetClickedFeatureStrokeWithHistory,
    featureStyleScopeDialog,
    featureStyleScopeLabelOptions,
    featureVisibilityScopeDialog,
    handleColorScopeChoice: handleColorScopeChoiceWithHistory,
    handleFeatureStyleScopeChoice: handleFeatureStyleScopeChoiceWithHistory,
    legendRenameDialog,
    handleLegendRenameChoice: handleLegendRenameChoiceWithHistory,
    paletteColorsDialog,
    handlePaletteColorsChoice: handlePaletteColorsChoiceWithHistory,
    resetColorDialog,
    handleResetColorChoice: handleResetColorChoiceWithHistory,
    labelTextScopeDialog,
    hiddenLabelTextDialog,
    labelOnDialog,
    // The dialog answers the Apply whose History step is still open.
    handleLabelOnChoice,
    clickedFeatureLabelHint,
    hiddenLabelTextMessage,
    updateClickedFeatureLabelText: updateClickedFeatureLabelTextWithHistory,
    handleLabelTextScopeChoice: handleLabelTextScopeChoiceWithHistory,
    handleHiddenLabelTextChoice: handleHiddenLabelTextChoiceWithHistory,
    requestLabelTextChangeByFeatureId: requestLabelTextChangeByFeatureIdWithHistory,
    commitListLabelDraft,
    requestLabelTextChangeByKey: requestLabelTextChangeByKeyWithHistory,
    resetAllLabelTextOverrides: resetAllLabelTextOverridesWithHistory,
    downloadLabelOverrideTable,
    loadLabelOverrideTable: loadLabelOverrideTableWithHistory,
    downloadFeatureEditTable: featureActions.downloadFeatureEditTable,
    loadFeatureEditTable: undoableAction(msg('Load feature edits'), featureActions.loadFeatureEditTable),
    syncLabelEditor,
    openFeatureEditorFromList,
    legendEntries: drawingMember('legendEntries'),
    legendEntryColor,
    previewItemView: previewItemPopup.view,
    setPreviewRecordField: previewItemPopup.setRecordField,
    showPreviewRecordSettings: previewItemPopup.showRecordSettings,
    renameLegendEntryInPopup: previewItemPopup.renameLegendRow,
    commitPreviewItemField: previewItemPopup.commitField,
    deletedLegendEntryColor,
    updateLegendEntryColor: updateLegendEntryColorShown,
    updateLegendEntryColorInPopup,
    defaultColorScopeLine,
    renameLegendEntry,
    ...legendActions.legendListView,
    // A Legend row delete is one History step of the intent, which Undo and
    // Redo show through the port and lay out as Python does, canvas included
    // (OV-125). The deleted row's stroke leaves its features (OV-293).
    deleteLegendEntry: editEditorIntent(msg('Delete legend item'), deleteLegendEntry, { domains: [...LEGEND_STRUCTURE_DOMAINS, ...STROKE_DOMAINS] }),
    // OV-154: a Restore returns deleted rows in one checkpoint step, as a delete
    // removes them; the palette then reaches the returned rows.
    deletedLegendEntries: drawingMember('deletedLegendEntries'),
    restoreDeletedLegendEntry: /** @param {number} index */ (index) => restoreLegendItems(
      msg('Restore legend item'), () => restoreDeletedLegendEntries([index])
    ),
    restoreAllDeletedLegendEntries: () => restoreLegendItems(msg('Restore legend items'), () => restoreDeletedLegendEntries(null)),
    moveLegendEntry,
    legendReorder,
    handleLegendReorderKey,
    startLegendReorder,
    overLegendReorder,
    dropLegendReorder,
    endLegendReorder,
    sortLegendEntries,
    sortLegendEntriesByDefault,
    resetLegendPosition: undoableAction(msg('Reset legend position'), resetLegendPosition),
    getLegendEntryStrokeColor,
    getLegendEntryStrokeWidth,
    isLegendStrokeOptionsOpen,
    toggleLegendStrokeOptions,
    setLegendEntryStrokeColorValue: setLegendEntryStrokeColorValueWithHistory,
    updateLegendEntryStrokeColor: editEditorIntent(msg('Change legend stroke color'), updateLegendEntryStrokeColor),
    updateLegendEntryStrokeWidth: editEditorIntent(msg('Change legend stroke width'), updateLegendEntryStrokeWidth),
    resetLegendEntryStroke: editEditorIntent(msg('Reset legend stroke'), resetLegendEntryStroke),
    resetAllStrokes: editEditorIntent(msg('Reset all strokes'), resetAllStrokes),
    resetAllPositions: undoableAction(msg('Reset positions'), resetAllPositions),
    resetLayout: undoableAction(msg('Reset layout'), resetLayout),
    layoutResetAvailable,
    canvasPadding: drawingMember('canvasPadding'),
    showCanvasControls,
    resetCanvasPadding,
    setCanvasPaddingSide: legendLayout.setCanvasPaddingSide,
    definitionLineStyleRows,
    linearLabelVisibilitySummary,
    linearLabelAutoFields,
    linearLabelAutoDisclosure,
    focusLinearLabelVisibility,
    legendPositionLabel,
    getDefinitionLineStyleSize,
    setDefinitionLineStyleSize,
    getDefinitionLineStyleWeight,
    setDefinitionLineStyleWeight,
    getDefinitionLineStyleFill,
    setDefinitionLineStyleColor,
    getDefinitionLineStyleColorMode,
    getDefinitionLineStyleSwatchValue,
    isDefinitionLineStyleMuted,
    downloadDpi,
    downloadWhiteBackground,
    pngDpiChoices: PNG_DPI_CHOICES,
    runAnalysis,
    generateFromError,
    cancelGeneration,
    downloadSVG,
    downloadInteractiveSVG,
    downloadPNG,
    downloadPDF,
    copyRunCommand,
    copyExactReplayCommand,
    downloadCliHelperFiles,
    runInfoElapsedText,
    runInfoReproducibilityText,
    shownRunInfoNoteParts,
    runInfoHasCliHelperFiles,
    resetSettings,
    resetShownDrawing,
    requestResetSettings,
    saveSessionWithTitle,
    editSessionTitle,
    importSession,
    sessionImportRollbackState,
    sessionInput,
    openSessionFilePicker,
    galleryExamples,
    exampleChooserOpen,
    openExampleChooser,
    closeExampleChooser,
    loadExampleSession,
    exampleSessionPending,
    choicePrompt: choicePrompt.prompt,
    resolveChoicePrompt: choicePrompt.resolve,
    submitChoicePromptField: choicePrompt.submitField,
    circularRegionNote,
    circularWindowRows,
    linearRegionNote,
    linearCardWindowView, linearLayoutCellRows,
    notice: noticeBanner.notice,
    dismissNotice: noticeBanner.dismiss,
    circularRecordPresentationPanel,
    canUndoHistory,
    canRedoHistory,
    undoHistoryTitle,
    redoHistoryTitle,
    undoHistory,
    redoHistory,
    historyList,
    manualPriorityRules: drawingMember('manualPriorityRules'),
    newPriorityRule,
    addPriorityRule
  };
};
