// @ts-check
/** @import { RecordPlacementRun } from './align-at-feature.js' */
/** @import { GenerationPayload } from '../services/diagram-generation.js' */
/** @import { RecordSpanRead } from '../services/diagram-worker-protocol.js' */
/** @import { AskingStackPorts, SwitchedTrackRow } from './track-slot-edits.js' */
/** @import { SvgMutationOperations } from '../services/svg-result-ingestion.js' */
/** @import { ResultLabelFit } from '../services/label-fit.js' */
/** @import { DiagramMode, DrawingState, LinearSeq } from '../state.js' */
/** @import { KeptLegendRow, LegendRow } from '../services/drawing/paint.js' */
/** @import { NoticeLink } from './app-messages.js' */
/** @import { InputFileCandidate, InputFileRef } from '../services/history-files.js' */
/** @import { DrawingInputSite } from './drawing-list.js' */
/** @import { ReplaceableInputSlot } from '../services/drawing/inputs.js' */
/** @import { RulePreparation } from './rule-matching.js' */
/** @import { CircularDiagramRecords } from '../services/record-options.js' */
/** @import { EditorPaintState } from './result-paint-record.js' */
/** @import { LinearDraftRecord } from '../services/linear-wrap-form.js' */
/** @import { FeatureEditorOptions } from './feature-editor.js' */
/** @import { ColorActionsRuleActions } from './feature-editor/color-actions.js' */
/** @import { LegendRowColorOf } from './candidate-render.js' */
/** @import { UserFacingError } from '../utils/error-normalization.js' */
/** @import { AnnotationCatalogSource } from './annotations/record-catalog.js' */
/** @import { LegacyResultSvgData, LegacyResultSvgTransform, SessionFileInput } from '../services/config.js' */
/** @import { RecordKeyedFeature } from '../services/feature-identity.js' */
/** @import { ArtifactSlot, ArtifactSlotSnapshot, ArtifactSlotStore, LoadedArtifactSlots } from '../services/artifact-slot.js' */
/** @import { DrawingId, DrawingIdentity } from '../services/drawing-view.js' */
/** @import { DisplayedLegendRowContext } from '../services/feature-visibility.js' */
/** @import { GalleryExample } from '../services/gallery-examples.js' */
/** @import { LinearSourceGroup } from '../services/linear-sources.js' */
/** @import { RecordDrawKey, RecordIdentitySources, RecordSettingsOwner, RecordEditData } from '../services/record-draw-selection.js' */
/** @import { RegionSourceFeature } from './region-drawing.js' */
/** @import { RecordListRequest, RecordSelectionMode, RecordSource } from './record-selection.js' */
/** @import { CircularRecordsDiscovered } from './run-analysis.js' */
import { createRulePreparation } from './rule-matching.js';
import { compileDirectEditorMutationPlan, displayedLegendRowColors, editorPaintDomains, LIVE_EDIT_DOMAINS, RECORDED_LEGEND_ROW_COLORS } from './candidate-render.js';
import { createResultPaintRecord, drawnLabelIntentGone, editorPaintState, paintRuleTable, sameColors, visibilityPaintState } from './result-paint-record.js';
import {
  countUnresolvedFeatureEdits, displayedLegendRowContext, featureDrawnContext, removeUnresolvedFeatureEdits,
  requestFeatureVisibilityRules
} from '../services/feature-visibility.js';
import { drawnBlockStroke, drawsPythonLegendRow, legendStructureEdited } from '../services/legend-svg.js';
import { legendEntryKey } from '../services/drawing/paint.js';
import { SETTINGS_SLICE } from '../services/drawing/settings.js';
import { DEFAULT_LINEAR_OUTPUT_PREFIX } from '../utils/download-names.js';
import { admitFeatureCatalog, biologicalFeatureKey, catalogNameIndex, resultCatalogFeatures } from '../services/feature-catalog.js';
import { labelCapacityNoticeText, labelCapacityText } from '../services/label-fit.js';
import { inputSummaryOf } from '../services/input-summary.js';
import { prefetchPyodideBootAssets } from '../services/pyodide-assets.js';
import { labelSettingsVisible } from '../services/feature-placement.js';
import { displayedFeatureAddressing } from '../services/feature-override-identity.js';
import { PNG_DPI_CHOICES } from '../services/session-active-config-contract.js';
import { settingsScrollOwner } from '../services/settings-sections.js';
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
import { colorRuleLegendGroups, groupLegendRowKeys, legendRowRules, specificRuleHasQualifier } from '../services/specific-color-rules.js';
import { namedFeatureInstance } from '../services/feature-selector.js';
import { drawnTracks, hiddenTrackLegendKeys, loadedTracks, trackLegendRows } from '../services/track-legend-rows.js';
import { isCurrentFeature } from '../services/feature-identity.js';
import {
  createFeatureRecordRotationAction,
  createFeatureRecordRotationWorkflow
} from './record-display/feature-record-rotation.js';
import { createDrawingState, drawingHasEdits, state, sessionOperationAvailability } from '../state.js';
import { createDrawingBar, createDrawingList, createInputReplacement, moveTabFocus } from './drawing-list.js';
import { createRegionDrawing } from './region-drawing.js';
import { createRegionSelect } from './region-select.js';
import { readRegionAxisRecords } from '../services/region-axis.js';
import { requestRecordSpanSequences } from '../services/record-span-sequences.js';
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
  committedTrackHideDrawing,
  commitLiveSettingValues,
  readCommittedCanonicalRevision,
  getLoadedTrackSettings,
  getCommittedLinearDefinitionVisibility,
  readCommittedResourceRecordCount,
  assertDrawingInputs,
  importSession as importSessionFromFile,
  serializeActiveRenderFiles,
  encodeSessionResources,
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
  CIRCULAR_RECORD_FORM_DEFAULTS
} from '../services/config.js';
import { artifactSlotAs } from '../services/artifact-slot.js';
import {
  setStoredProteinResultFrameConverter
} from '../services/main-session-comparison-frame.js';
import { setLegacySessionUpgrader } from '../services/legacy-session-upgrade.js';
import { fetchGalleryExampleSession, readGalleryExamples } from '../services/gallery-examples.js';
import { createHistoryRouter } from '../services/history.js';
import { bindDrawingInputFile, bindDrawingInputFileTree, createLinearSeq, retainHistoryAndBoundFiles, boundInputFileIds } from '../services/drawing/inputs.js';
import { CIRCULAR_INPUT_FIELDS, circularInputSlots, inputSourceFiles, inputTypeSpec } from '../services/input-types.js';
import { createHistorySnapshotService } from '../services/history-snapshot.js';
import {
  groupLinearSourceRecords,
  linearFileCardRecord,
  planLinearSourceRemoval,
  getLinearSourceDefaultDefinition,
  inferredDefinitionForRecord,
  setLinearSourceDefaultDefinition,
  getLinearSourceDefaultSubtitle,
  setLinearSourceDefaultSubtitle,
  resolveLinearRecordEffectiveDefinition,
  resolveLinearRecordEffectiveSubtitle
} from '../services/linear-sources.js';
import { captureSvgExport, serializeCleanSvg } from '../services/svg-serialization.js';
import { copyTextToClipboard } from '../utils/clipboard.js';
import { toRaw } from '../utils/vue-raw.js';
import { readHostedBuild } from '../utils/hosted-build.js';
import { resetLayoutState, resetProject, resetSettings as resetDrawingSettings, withNewDrawingDefaults } from '../services/reset.js';
import { allLiveSettingRequestPaths, liveSettingRequestPaths, SETTING_NOTE_APPLIES, changedTrackHides } from '../services/setting-application.js';
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
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from '../services/runtime-test-hooks.js';
import { createPanZoom, createSidebarResize, setupGlobalUiEvents } from './ui.js';
import { appliedFeatureColors, COMPARISON_COLOR_KEYS, normalizeOptionalHexColor } from '../utils/color-utils.js';
import { createFeatureEditor } from './feature-editor.js';
import { COLLINEAR_COLOR_MODE_TEXT } from './pairwise-match-popup.js';
import { createFeatureSelection } from './feature-selection.js';
import { createPreviewFeatureSearch } from './feature-search/preview-actions.js';
import { createSvgStyles, paletteColorsEqual } from './svg-styles.js';
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
import {
  buildDefaultCircularRecordPositions, createCircularRecordOrder, createCircularRecordPresentation
} from './circular-record-panel.js';
import {
  createClickedFeatureCopy, createImportedComparisonChoice, createMissingBlastTsvChoice, createResultExport, createRunInfoCopy
} from './result-actions.js';
import { diagnosticError, normalizeCaughtError, normalizeUserFacingError } from '../utils/error-normalization.js';
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
import { askTrackSourceRemoval, askingTrackRowRemoval, askingTrackRowSwitch } from './track-slot-edits.js';
import { createAnnotationEditor } from './annotations.js';
import { createDrawingLegendRetirement, trackDataLegendCaptions } from './legend/track-data-styles.js';
import {
  annotationRecordKeysOf,
  annotationSourceKey,
  buildAnnotationRecordCatalog
} from './annotations/record-catalog.js';
import {
  reconcileAnnotationRecordBindings, resultRecordBindings
} from './annotations/record-selector.js';
import { validateAnnotationRecordTargets } from './annotations/validation.js';
import { createLosatSettings } from './losat-settings.js';
import { createAutoValueDisplay } from './auto-value-display.js';
import { createLinearWrapField } from './linear-wrap-field.js';
import { createFigureSizeField } from './figure-size-field.js';
import { createFigureSizeReaders } from './figure-size-readers.js';
import { createOrthogroupPopup } from './orthogroup-popup.js';
import { createSimilarityAlignmentReview, createSimilarityAlignmentTakeover } from './similarity-alignment-review.js';
import { requestLinearWrapLeaf, shownLinearWrap, wrapFormValue } from '../services/linear-wrap-form.js';
import { createLinearRecordSelector } from './linear-record-selector.js';
import { createLinearTypographyController } from './linear-typography.js';
import {
  buildDisambiguatedRecordEntries,
  circularDrawingRecordSet,
  circularDrawingRegionNote,
  circularEstimateLengths,
  circularRequestDiagrams,
  formatRecordLength,
  regionClampNote,
  linearCardRegion,
  resolveDisambiguatedRecordSelection
} from '../services/record-options.js';
import {
  setLinearRecordRow as updateLinearRecordRow
} from '../services/linear-record-layout.js';
import {
  missingUploadPairsToResolve
} from '../services/linear-comparisons.js';
import {
  RECORD_SETTINGS_DEFAULTS,
  circularRecordDrawKey,
  circularRecordRequestKeys,
  linearCardHasSettings,
  recordsWithEdits,
  removeRecordEdits
} from '../services/record-draw-selection.js';
import { ANNOTATION_RECORD_BINDING_KEY } from '../services/annotation-state.js';
import { circularMultiRecordDefaults, updateActiveLayoutPreference } from '../services/layout-preferences.js';
import { addAnnotationSet as addSearchSetToDrawing, removeAnnotations } from '../services/drawing/annotation-sets.js';
import { createSequenceSearch } from './sequence-search.js';
import { resultRecordBinding } from '../services/sequence-search.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';
import { createRecordSelection } from './record-selection.js';
import {
  circularDiscoveryForInput,
  discoverInputRecords
} from '../services/record-discovery.js';
import { normalizeFileList } from '../services/conservation-series.js';
import {
  comparisonProfileDefault, IDENTITY_PERCENT_MAXIMUM, identityScaleWithEnd, linearAxisRulerApplies
} from '../mode-profiles.js';
import { joinSentences, listText, msg, t, tm, tn, tx } from '../utils/i18n.js';
import { errorText } from '../utils/error-text.js';
import {
  IMPORTED_COMPARISON_DISPOSITIONS,
  importedComparisonExecution
} from '../services/imported-comparison-intent.js';
import {
  activeDepthTrackIndices,
  compactDepthFileSlots,
  depthFileSlotsFromValue,
  depthSlotTrackIndex,
  isDefaultManagedDepthSlot,
  representativeDepthFiles,
  reindexDepthSlots,
  removeDepthTrackColumnAt,
  depthSeriesRemovalLoss
} from '../services/depth-track-state.js';
/** @import { HistoryManagerOptions } from '../services/history.js' */
/** @import { RulePreparationOptions } from './rule-matching.js' */
/** @import { AnnotationRecordCatalog } from './annotations/record-catalog.js' */
/** @import { SpecificColorRule } from '../services/specific-color-rules.js' */
/** @import { MountedResultContext } from './preview-runtime.js' */
/** @import { UnresolvedNotice } from '../services/feature-visibility.js' */
/** @import { OrthogroupGroup } from './orthogroups.js' */
/** @import { ReindexDepthSlotsOptions } from '../services/depth-track-state.js' */
/** @import { CircularTrackSlot } from '../services/circular-track-slot-model.js' */
/** @import { LinearTrackSlot } from '../services/linear-track-slot-model.js' */
/** @import { SequenceRecord } from '../services/record-discovery.js' */
/** @import { AlignmentOrientationCheckpoint, TargetDraftCheckpoint } from './record-display-options.js' */
/** @import { UnvalidatedFields } from '../services/session-authority.js' */

const { onMounted, onUnmounted, watch, nextTick, computed, provide, ref, shallowRef, reactive, isRef } = window.Vue;

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
/** The records of every Result of the shown drawing whose axis Python wrote (the subjects of Sequence search). */
const shownResultsAxisRecords = () => Array.from({ length: state.results.value.length },
  (_result, index) => readRegionAxisRecords(state.trackSlotResolvedGeometry.value, index)).flat();
const closeRegionPopup = () => {
  if (state.clickedPairwiseMatch.value?.previewItem?.kind === 'region') state.clickedPairwiseMatch.value = null;
};
// The bases of spans of the shown drawing's records, read by Python from the committed request, the one the shown
// Result was made from: Copy FASTA and Use as search query of Select region, and the subjects of Sequence search.
/** @param {readonly RecordSpanRead[]} spans */
const readShownRecordSpans = (spans) => requestRecordSpanSequences(
  getCommittedCanonicalSession(), spans,
  (payload) => runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.READ_RECORD_SPAN_SEQUENCES, { ...payload,
    resources: /** @type {GenerationPayload['resources']} */ (payload.resources) })
);
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

const HISTORY_RESTORE_BUSY = Object.freeze({ status: 'busy', reason: msg('Undo or Redo in progress. Retry after it finishes.') });
const RULE_WRITE_BUSY = Object.freeze({ status: 'busy', reason: 'A color rule edit is being matched. Retry after it finishes.' });
export const createAppSetup = () => {
  provide(SETTING_NOTE_APPLIES, committedSettingNotePredicate);
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
    featureIdentityNotices,
    featureEditRemovalCount,
    comparisonWarnings,
    labelFitReports,
    inputSummaries,
    matchSequenceRegistry,
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
    clickedLabel,
    clickedLabelPos,
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
    showLinearComparisonCause
  } = createLinearComparisonFocus({ state, mode, errorLog, linearComparisonUi, revealSetting: settingsPaneActions.revealSetting });
  const { shownLinearRecordRowFor, linearLayoutTokens, linearLosatCacheInfoByEdgeKey } = createLinearLayoutReaders({
    state, shownDrawingOfMode, losatCacheInfo
  });

  /**
   * What a drawing's next record discoveries do, for the rows and input of
   * that drawing only (OV-2148: a Duplicate keeps its row uids, so a uid alone
   * names a row of several drawings).
   * - `expansions`: the Linear rows whose next discovery may expand them, each
   *   with whether the row's own names were entered for the file being read.
   *   They were not when the file replaced another: those names belong to it.
   * - `uploads`, `circularUpload`: the Linear Files and the Circular input just
   *   uploaded: their first discovery infers record definitions (Linear
   *   GenBank) and asks for the record list of a file with many records (D-04).
   *   A Session load, a rollback, and a History step never add one.
   * - `listRequests`: those asks, in upload order; the drawing's panel opens them.
   * @typedef {{ expansions: Map<string, boolean>, uploads: Set<string>, circularUpload: InputFileRef | null,
   *   listRequests: RecordListRequest[] }} PendingRecordDiscovery
   */
  /** @type {WeakMap<DrawingState, PendingRecordDiscovery>} */
  const pendingDiscoveries = new WeakMap();
  /** @param {DrawingState} drawing @returns {PendingRecordDiscovery} */
  const pendingDiscoveryOf = (drawing) => {
    const known = pendingDiscoveries.get(drawing);
    if (known) return known;
    /** @type {PendingRecordDiscovery} */
    const pending = { expansions: new Map(), uploads: new Set(), circularUpload: null, listRequests: reactive([]) };
    pendingDiscoveries.set(drawing, pending);
    return pending;
  };
  /**
   * @param {DrawingState} drawing
   * @param {RecordSelectionMode} listMode @param {string} sourceKey @param {number} recordCount
   */
  const requestRecordListAfterUpload = (drawing, listMode, sourceKey, recordCount) => {
    if (recordCount > WEB_UX_PROFILE.recordList.autoOpenAbove) {
      pendingDiscoveryOf(drawing).listRequests.push({ mode: listMode, sourceKey });
    }
  };
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
  /**
   * The drawing a discovery ran for, or null when it was deleted while the
   * discovery ran: such a discovery changes nothing (OV-2148).
   * @param {DrawingId} drawingId @returns {DrawingState | null}
   */
  const listedDrawingOf = (drawingId) => (
    state.drawingTabs.value.some((tab) => tab.id === drawingId) ? state.drawingOf(drawingId) : null
  );
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
  /**
   * A Circular drawing's discovery settled (success or error): the first after
   * its input's upload asks for the record list.
   * @param {CircularRecordsDiscovered} discovered
   */
  const handleCircularRecordsDiscovered = ({ drawingId, primaryFile, pairedFile, recordCount }) => {
    rebindAnnotationRecords();
    const drawing = listedDrawingOf(drawingId);
    if (!drawing) return;
    const pending = pendingDiscoveryOf(drawing);
    const upload = pending.circularUpload;
    if (!upload || (upload !== primaryFile && upload !== pairedFile)) return;
    pending.circularUpload = null;
    requestRecordListAfterUpload(drawing, 'circular', 'circular', recordCount);
  };
  /** @param {DrawingState} drawing */
  const materializeAutomaticLinearRecords = async (drawing) => {
    if (drawing.mode !== 'linear') return;
    // A row that already waits for its discovery keeps what its writer said
    // about its names (a replacing file does not own the row's record name).
    const { expansions } = pendingDiscoveryOf(drawing);
    drawing.linearSeqs.forEach((seq) => {
      if (!seq.region_record_id && !expansions.has(seq.uid)) expansions.set(seq.uid, true);
    });
    await linearRecordSelector.refresh(drawing);
  };
  /**
   * The records of one input: a GenBank file, a GFF3 file with its FASTA, or a FASTA alone.
   * The bindings name the files by ref; their bytes come from the store.
   * @param {string} inputType
   * @param {InputFileRef | null} primaryFile
   * @param {InputFileRef | null} pairedFile
   */
  const readInputRecords = (inputType, primaryFile, pairedFile) => discoverInputRecords(
    inputTypeSpec(inputType), state.inputFiles.fileOf(primaryFile), state.inputFiles.fileOf(pairedFile)
  );
  const paletteLoader = createPaletteLoader({ state });
  const linearRecordSelector = createLinearRecordSelector({
    state,
    reactive,
    onRecordsDiscovered: handleLinearRecordsDiscovered,
    recordReader: ({ inputType, primaryFile, pairedFile }) => readInputRecords(inputType, primaryFile, pairedFile)
  });
  /** @param {DrawingState} drawing */
  const getCircularRecordDiscoveryState = (drawing) => circularDiscoveryForInput(drawing);
  /** @param {DrawingState} drawing */
  const circularAnnotationSourceKey = (drawing) => annotationSourceKey({
    scope: 'circular',
    inputType: drawing.inputType.value,
    ...inputSourceFiles(drawing.inputType.value, circularInputSlots(drawing.files))
  });
  /**
   * A Linear record card as the annotation record catalog reads it.
   * @param {DrawingState} drawing the drawing that holds the card
   * @param {LinearSeq} seq
   * @param {boolean} drawn
   * @returns {AnnotationCatalogSource}
   */
  const linearAnnotationSource = (drawing, seq, drawn) => {
    const inputType = drawing.inputType.value;
    const { primaryFile, pairedFile, complete } = inputSourceFiles(inputType, seq);
    return {
      sourceKey: annotationSourceKey({
        scope: 'linear',
        uid: seq.uid,
        inputType,
        primaryFile,
        pairedFile
      }),
      selector: seq.region_record_id,
      drawn,
      hasInput: complete,
      ...linearRecordSelector.catalogStateFor(drawing, seq)
    };
  };
  /**
   * The records one drawing's annotations can name: a Circular drawing's
   * record set, or a Linear drawing's drawn cards (`linearSourcesOverride`
   * replaces the cards a Linear drawing reads).
   * @param {DrawingState} drawing
   * @param {AnnotationCatalogSource[] | null} [linearSourcesOverride]
   */
  const getAnnotationRecordCatalog = (drawing, linearSourcesOverride = null) => {
    const circularSource = () => {
      const circularDiscovery = getCircularRecordDiscoveryState(drawing);
      const circularRecordSet = circularDrawingRecordSet(circularDiscovery.records, drawing);
      return {
        sourceKey: circularAnnotationSourceKey(drawing),
        hasInput: circularDiscovery.hasInput,
        status: circularDiscovery.status,
        error: circularDiscovery.error,
        records: circularRecordSet.records,
        omittedRecords: circularRecordSet.omittedRecords
      };
    };
    const linearSources = () => {
      const recordsOff = new Set(drawing.recordsOff);
      return linearSourcesOverride || drawing.linearSeqs.map((seq) => linearAnnotationSource(drawing, seq, !recordsOff.has(seq.uid)));
    };
    return buildAnnotationRecordCatalog(drawing.mode === 'linear'
      ? { mode: 'linear', linearSources: linearSources() }
      : { mode: 'circular', circularSource: circularSource() });
  };
  // R10 (D-15-46): each transition that changes the annotation record catalog
  // (a discovery that ends, a mode switch, record selection, a Linear source
  // edit, a record selector, Multi-Record Canvas) rebinds the shown drawing's
  // annotations to the records the catalog resolves, in that transition's
  // History step. Load and its rollback are busy and end in a discovery; a
  // History restore installs the bindings it captured; Generate rebinds in
  // `validateAnnotationTargets`. A transition that leaves the catalog's
  // status and records as the last rebind saw them rebinds nothing (a refresh
  // after Save or at the end of Load), as the former catalog watcher did.
  let reboundCatalog = '';
  const rebindAnnotationRecords = () => {
    const drawing = state.activeDrawing();
    const catalog = getAnnotationRecordCatalog(drawing);
    const seen = `${catalog.status}:${catalog.signature}`;
    if (seen === reboundCatalog) return;
    reboundCatalog = seen;
    if (catalog.status === 'ready' && !sessionOperationAvailability()) {
      reconcileAnnotationRecordBindings(drawing, catalog);
    }
  };
  // A Circular drawing whose Multi-Record Canvas is on and has no multi-record
  // sides stored gets the toggle's smart defaults. Load, its rollback, and a
  // History restore install the sides they read or captured (R11, D-15-48).
  /** @param {DrawingState} drawing */
  const storeCircularMultiRecordDefaults = (drawing) => {
    const defaults = drawing.mode === 'circular' && drawing.form.multi_record_canvas
      ? circularMultiRecordDefaults(drawing.layoutPreferences) : null;
    if (defaults) updateActiveLayoutPreference(drawing.layoutPreferences, 'circular', true, defaults);
  };
  // The Multi-Record Canvas checkbox, after its `v-model` wrote the shown
  // drawing: in the control's History step.
  const changeCircularMultiRecordCanvas = () => {
    storeCircularMultiRecordDefaults(state.activeDrawing());
    rebindAnnotationRecords();
  };
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
  const orthogroupActions = createOrthogroupEditor({ state, getCommittedRequest: getCommittedCanonicalRenderRequest });
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
  // rerender (D-24 (OV-731)).
  /**
   * @template T
   * @param {DrawingState} drawing
   * @param {() => Promise<T>} write
   * @param {{ drawn: boolean }} [restored]
   * @returns {Promise<T>}
   */
  const followRuleEdits = async (drawing, write, { drawn } = { drawn: false }) => {
    // The rule owner follows the shown drawing; a drawing that is not shown
    // follows its rules when it is shown again (D-15-51).
    if (drawing !== state.activeDrawing()) return write();
    const previousRules = drawing.manualSpecificRules.map((/** @type {SpecificColorRule} */ rule) => ({ ...rule }));
    specificRuleRestorePorts.retainRulesForRestore(previousRules);
    try {
      const written = await write();
      if (!drawn) specificRuleRestorePorts.followRestoredSpecificRules(previousRules);
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
    estimateFileBytes: state.inputFiles.estimateBytes,
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
    linearRecordsFor: linearRecordSelector.recordsFor,
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
  /** @param {unknown} value */
  const displayPreservedConfigValue = (value) => {
    const serialized = JSON.stringify(value) ?? String(value);
    return serialized.length <= 240
      ? serialized
      : `${serialized.slice(0, 237)}...`;
  };
  const unmanagedConfigOverrideEntries = computed(() => (
    Object.entries(state.activeDrawing().unmanagedConfigOverrides)
      .sort(([left], [right]) => left.localeCompare(right))
      .map(([path, value]) => ({
        path,
        value,
        displayValue: displayPreservedConfigValue(value)
      }))
  ));
  /** @param {string} path */
  const resetUnmanagedConfigOverride = (path) => {
    const drawing = state.activeDrawing();
    return history.runUndoable(
      msg('Reset preserved setting'),
      () => {
        if (!Object.prototype.hasOwnProperty.call(drawing.unmanagedConfigOverrides, path)) {
          return false;
        }
        delete drawing.unmanagedConfigOverrides[path];
        return true;
      }
    );
  };
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
    previewTransformInteraction,
    holdPreviewViewport,
    releasePreviewViewport,
    heldPreviewViewport
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
    synthesize: async (request) => (await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.SYNTHESIZE_COLOR_RULES, { request })).result,
    visibilityRules: () => requestFeatureVisibilityRules(state.activeDrawing().featureVisibilityManualRules)
  });
  specificRuleRestorePorts.retainRulesForRestore = rulePreparation.retain;
  // R13: a Legend row a specific-color rule draws commits its edit through the
  // rule owner; the root registers the port once the feature editor exists.
  /** @type {{ commitLegendRowRules: ColorActionsRuleActions['commitSpecificRules'] }} */
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
  const shownAxisRecords = () => readRegionAxisRecords(state.trackSlotResolvedGeometry.value, Number(selectedResultIndex.value) || 0);
  const featureSelection = createFeatureSelection({
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
    previewPopupsHeld,
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
    featureSelection,
    commitActiveResultEdit: previewRuntime.commitActiveResultEdit,
    currentResults: previewRuntime.currentResults,
    getResultIdentity: previewRuntime.getResultIdentity,
    showEditorIntent: (options) => showEditorIntent(options),
    selectResult,
    previewTransformInteraction,
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
      return orthogroups.value.map((/** @type {OrthogroupGroup} */ group) => ({
        ...group,
        display_name: orthogroupActions.resolveOrthogroupName(group, drawing),
        description: orthogroupActions.resolveOrthogroupDescription(group, drawing)
      }));
    },
    openFeatureEditorForFeature: featureActions.openFeatureEditorForFeature,
    focusSearchQuery: () => /** @type {HTMLElement | null} */ (document.querySelector('.preview-feature-search .gfs-query'))?.focus(),
    focusSearchToggle: () => /** @type {HTMLElement | null} */ (document.querySelector('[data-preview-search-toggle]'))?.focus()
  });

  // A Result switch closes the feature popup through its close owner: the
  // popup belongs to the Result it was opened on (OV-1269).
  watch(selectedResultIndex, () => {
    featureSelection.clearFeatureSelection({ clearStatus: true, syncDom: false });
    closeFeaturePopup();
  });
  watch(svgContent, () => {
    similarityAlignmentPorts.refreshCanvas();
    if (!skipCaptureBaseConfig.value) {
      featureSelection.clearFeatureSelection({ clearStatus: true, syncDom: false });
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
  const linearTrackSlotEditor = createLinearTrackSlotEditor({
    state, changeTrackLayout, changeTrackVisibility: retireShownLegendStyles
  });
  // PD-OI-104 scenario revisions 6 and 7: a row's own switch, a row's Remove, and
  // the removal of a Depth file or series ask before they take out of the shown
  // stack a row that other rows start at or extend to.
  /** @param {'circular_track_slots' | 'linear_track_slots'} slotsKey */
  const askingStackPorts = (slotsKey) => ({
    askChoice, runUndoable: history.runUndoable, busy: sessionOperationAvailability,
    rows: () => /** @type {readonly SwitchedTrackRow[]} */ (state.activeDrawing().adv[slotsKey]),
    mode: /** @type {'circular' | 'linear'} */ (slotsKey === 'circular_track_slots' ? 'circular' : 'linear')
  });
  /** @template {SwitchedTrackRow} Row @param {(slot: Row | null | undefined, enabled: unknown, turnOffReaching?: (() => void) | Event) => unknown} setEnabled @param {'circular_track_slots' | 'linear_track_slots'} slotsKey */
  const askingRowSwitch = (setEnabled, slotsKey) => askingTrackRowSwitch(/** @type {AskingStackPorts<Row>} */ (askingStackPorts(slotsKey)), setEnabled);
  /** @param {(index: unknown, turnOffReaching?: (() => void) | Event) => unknown} removeRow @param {'circular_track_slots' | 'linear_track_slots'} slotsKey */
  const askingRowRemoval = (removeRow, slotsKey) => askingTrackRowRemoval(askingStackPorts(slotsKey), removeRow);
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
  const losatSettings = createLosatSettings({ state });
  /**
   * The record a Linear entry draws, once its records are read: the one its
   * selector names, or the only record of its file.
   * @param {DrawingState} drawing the drawing that holds the card
   * @param {LinearSeq} seq
   */
  const linearSeqRecord = (drawing, seq) => {
    const records = linearRecordSelector.recordsFor(drawing, seq);
    const selection = resolveDisambiguatedRecordSelection(records, seq.region_record_id);
    return selection.record || (selection.status === 'unspecified' && records.length === 1 ? records[0] : null);
  };
  // The note beside a Linear entry's Region fields (`regionClampNote`). The
  // request draws the larger coordinate as the region's end (reverse when typed first).
  /** @param {LinearSeq} seq a card of the shown drawing */
  const linearRegionNote = (seq) => {
    const region = linearCardRegion(seq, 1, linearSeqRecord(state.activeDrawing(), seq)?.recordLength);
    return region && !('error' in region) ? regionClampNote(region.drawn) : '';
  };
  const circularRegionNote = computed(() => {
    const drawing = shownDrawingOfMode('circular');
    return drawing ? circularDrawingRegionNote(drawing.circularRecordList.value, drawing) : '';
  });
  // The input summary of the shown Result: its warnings in the status area, its lines in Run info, its
  // wrapped rows.
  const shownInputSummary = computed(() => inputSummaryOf(
    inputSummaries.value, selectedResultIndex.value, results.value[selectedResultIndex.value]?.name
  ));
  // The records the next Linear Generate draws, as their cards ask for them:
  // each source record's length and the region its card states. Known once
  // every drawn record is read and every region can be drawn, else none.
  /** @param {DrawingState | null | undefined} drawing @returns {LinearDraftRecord[]} */
  const linearDraftRecordsOf = (drawing) => {
    if (!drawing || drawing.mode !== 'linear') return [];
    const records = drawing.drawnLinearSeqs.value.map((seq, index) => {
      const sourceLength = Number(linearSeqRecord(drawing, seq)?.recordLength);
      const region = linearCardRegion(seq, index + 1);
      if (!(sourceLength > 0) || (region !== null && 'error' in region)) return null;
      return { sourceLength,
        region: region === null ? null : { start: Math.min(region.start, region.end), end: Math.max(region.start, region.end) } };
    });
    return records.length && records.every((record) => record !== null) ? /** @type {LinearDraftRecord[]} */ (records) : [];
  };
  // Their lengths, while none is cropped: the Web predicts no cropped length.
  /** @param {DrawingState | null | undefined} drawing @returns {number[]} */
  const linearRecordLengthsOf = (drawing) => {
    const records = linearDraftRecordsOf(drawing);
    return records.every((record) => record.region === null) ? records.map((record) => record.sourceLength) : [];
  };
  const drawnLinearRecordLengths = () => linearRecordLengthsOf(state.activeDrawing());
  // The records a shown Result without an input summary (a loaded one) drew:
  // each committed record's region as its request asked for it, with the
  // length of the source record its card reads, matched by record key, so a
  // card edited since Load changes neither (D-73). [] while a record has no
  // card that still draws it, or its length is unknown.
  /** @param {DrawingState} drawing @returns {LinearDraftRecord[]} */
  const committedLinearRecordsOf = (drawing) => {
    const records = (committedLinearRecordCards(drawing) ?? []).map(({ record, seq }) => {
      const sourceLength = Number(linearSeqRecord(drawing, seq)?.recordLength);
      const { region } = record;
      return sourceLength > 0 ? { sourceLength, region: region && { start: region.start, end: region.end } } : null;
    });
    return records.length && records.every((record) => record !== null) ? /** @type {LinearDraftRecord[]} */ (records) : [];
  };
  // A Linear card whose record list is still deferred after a Load is read
  // now, through the record selector, as Generate reads it (DESIGN A4).
  /** @param {DrawingState} drawing */
  const awaitLinearRecordDiscovery = async (drawing) => {
    if (drawing.linearSeqs.some((seq) => linearRecordSelector.statusFor(drawing, seq) === 'deferred')) {
      const outcome = linearRecordSelector.refresh(drawing);
      if (outcome instanceof Promise) await outcome;
    }
  };
  // Web-only changes of the shown canvas: a moved diagram, Legend, title, or scale, and Canvas padding.
  const decorationsMoved = () => [state.diagramOffset, state.legendCurrentOffset, state.plotTitleUserOffset, state.lengthBarUserOffset]
    .some((offset) => Boolean(Number(offset.x) || Number(offset.y)));
  /** @param {{ canvasPadding: Record<string, unknown> }} drawing */
  const canvasPadded = (drawing) => (/** @type {const} */ (['top', 'right', 'bottom', 'left'])).some((side) => Number(drawing.canvasPadding[side]));
  // The committed request of the Result a drawing shows; null while it shows no Result of its own.
  /** @param {{ id: string }} drawing */
  const shownCommittedRequest = (drawing) => (
    state.generatedDrawingId.value === drawing.id && results.value.length > 0 ? getCommittedCanonicalRenderRequest() : null
  );
  // Wrap Rows as the shown Result drew it (PD-OI-113), from Python's input
  // summary and the committed request of the drawing whose Result is shown:
  // the field's Auto text and the Auto hints of sizes and ticks read it.
  /** @param {DrawingState} drawing The committed request's stored wrap leaf while this drawing's Result is shown. */
  const shownLinearWrapLeaf = (drawing) => {
    const request = drawing.mode === 'linear' ? shownCommittedRequest(drawing) : null;
    return request?.mode === 'linear' ? requestLinearWrapLeaf(request) : undefined;
  };
  const shownLinearWrapReading = computed(() => {
    const drawing = state.activeDrawing();
    const leaf = shownLinearWrapLeaf(drawing);
    return shownLinearWrap({
      value: drawing.form.linear_wrap_length,
      requested: leaf === undefined ? undefined : wrapFormValue(leaf).value,
      summary: shownInputSummary.value ?? undefined,
      records: linearDraftRecordsOf(drawing)
    });
  });
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
  // a file whose inspection failed, and a record that TLOSATX cannot search
  // without a translation table.
  /** @param {LinearSourceGroup} source a File card of the shown drawing */
  const linearSourceInspectionFailed = (source) => Boolean(linearRecordSelector.errorFor(state.activeDrawing(), source.sequence));
  /** @param {LinearSeq} seq */
  const linearRecordGencodeRequired = (seq) => {
    const drawing = shownDrawingOfMode('linear');
    if (!drawing) return false;
    return drawing.hasActiveLinearLosatIntent.value && drawing.losatProgram.value === 'tblastx'
      && !(Number(seq.losat_gencode) >= 1);
  };
  /** @param {LinearSourceGroup} source */
  const linearSourceGencodeRequired = (source) => source.records.some(({ sequence }) => linearRecordGencodeRequired(sequence));
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
    clearDerivedCircularConservationBlasts, setCircularConservationSourceFiles, replaceCircularConservationBlasts,
    askRingSourceRemoval, reconcileCircularConservation, setCircularConservationEnabled, refreshCircularConservationSlots
  } = createCircularConservationSources({
    state, writeLosatCacheInfo,
    syncCircularConservationSlots: circularTrackSlotEditor.syncCircularConservationSlots,
    circularConservationRowsLeaving: circularTrackSlotEditor.circularConservationRowsLeaving,
    askingStackPorts
  });
  const {
    setCircularConservationUploadFiles, circularConservationSeriesRows, canMoveCircularConservationSeries,
    moveCircularConservationSeries, openCircularConservationComparisonFilePicker, changeCircularConservationSource,
    refreshShownCircularConservationSlots
  } = createCircularConservationSeriesPanel({
    state, shownAction, shownDrawingOfMode, writeLosatCacheInfo, circularConservationFastaInput,
    reconcileCircularConservation, refreshCircularConservationSlots, askRingSourceRemoval
  });
  const {
    settleComparisonRecordLabels, addCircularConservationComparisonFile, setCircularConservationCompanionFile,
    removeCircularConservationSource
  } = createCircularConservationComparisonFiles({
    state, shownAction, writeLosatCacheInfo, runUndoable: history.runUndoable, clearDerivedCircularConservationBlasts,
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
  /** @param {DrawingState} drawing */
  const shouldSyncMountedLabelEditor = (drawing) => (
    isFeatureDrawerMounted.value
    || Boolean(clickedFeature.value)
    || labelTextScopeDialog.show
    || Object.values(drawing.featureOverrides).some((row) => row.labelText !== null || row.labelVisibility !== null)
    || Object.keys(drawing.labelTextBulkOverrides).length > 0
  );
  // A legacy imported SVG without composition metadata stays unbound: it has
  // no composition to capture and no canvas to pad.
  // Mounted Results whose Legend `adoptLegend` laid out; the commit waits for
  // the later binding steps (`initializeStrokeAndCanvas`).
  /** @type {WeakSet<Element>} */
  const mountedLegendLayouts = new WeakSet();
  /** @param {MountedResultContext} context */
  const shouldBindComposition = (context) => (
    context.root.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) !== null
    || context.root.getAttribute(COMPOSITION_METADATA_ATTRIBUTE) !== null
    || (!context.bindingOptions.isIncrementalEdit && context.sourceClass !== 'legacy-import')
  );
  // Each Result's bytes reflect the editor state it was committed or last
  // shown with (`createResultPaintRecord`); the binder below reads it (R13).
  const resultPaintRecord = createResultPaintRecord();
  previewRuntime.configureMountedResultBinder({
    async adoptLegend(context) {
      const drawing = state.activeDrawing();
      // The Legend layout port (zero shift) reads the bundled-font metrics;
      // they load with the first Result that has a Legend, so the Legend
      // edits, which lay the Legend out synchronously, find them loaded.
      // Once loaded, the binder does not wait: an edit waiting for this
      // binding (a rule commit after a rerender) reads its Legend unchanged.
      if (context.root.getElementById?.('legend') && !legendActions.isLegendLayoutReady()) {
        await legendActions.prepareLegendLayout().catch((error) => {
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
        const resultLegendOrder = legendActions.captureResultInventory(context.root, {
          resultIdentity: context.resultIdentity,
          liveResultIdentities: liveResultIdentities()
        });
        // A Result shown again by a mode switch arrives in its own drawing,
        // whose Legend edits name that Result's inventory: the inventory is
        // adopted first, so the departing mode's order is never compared with
        // the arriving drawing's rows. A batch Result of the same drawing is
        // compared with the order shown until now (B20, OV-47).
        const arriving = Boolean(context.bindingOptions.modeArrival);
        if (arriving) legendActions.adoptResultInventory(context.resultIdentity);
        await projectEditorIntentOnDisplay(drawing, context, resultLegendOrder);
        if (!arriving) legendActions.adoptResultInventory(context.resultIdentity);
      } else {
        rememberCommittedEditorState(drawing, context);
      }
      if (context.bindingOptions.trustedRestore) {
        legendActions.adoptResultInventory(context.resultIdentity, { restored: true });
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
          legendActions.onLegendGeometryChanged({ commit: false });
          mountedLegendLayouts.add(context.root);
        }
        if (context.bindingOptions.skipLegendExtraction) return;
        recordStructuralMetric('legendDomFullScanCount', 1, {
          phase: context.phase,
          rootGeneration: context.rootGeneration
        });
        // A draw that replays an edited Legend order replays it on the displayed
        // mode's Results only.
        legendActions.extractLegendEntries({
          replaceGeneratedInventory: !selecting && drawn,
          liveRedraw: context.bindingOptions.liveRedrawLegend ?? null,
          liveResultIdentities: results.value.map(previewRuntime.getResultIdentity)
        });
      };
      // OV-872: the Legend a label rerender draws is no edit of a step begun
      // after it was requested (an action waiting for it, D-40).
      const { historyRequestedAt } = context.bindingOptions;
      if (historyRequestedAt === undefined) adopt();
      else history.adoptDrawnWrite(adopt, { requestedAt: historyRequestedAt });
    },
    bindComposition(context) {
      if (context.bindingOptions.trustedRestore) return;
      if (shouldBindComposition(context)) legendLayout.captureBaseConfig();
    },
    setupDragAffordances(context) {
      legendActions.setupLegendDrag();
      legendLayout.setupDiagramDrag(Boolean(context.bindingOptions.isIncrementalEdit));
    },
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
    synchronizeLabelEditor(context) {
      const drawing = state.activeDrawing();
      const labelsChanged = labelProjectionResultIdentity === context.resultIdentity;
      labelProjectionResultIdentity = '';
      if (!context.bindingOptions.trustedRestore && (labelsChanged || shouldSyncMountedLabelEditor(drawing))) {
        featureActions.syncLabelEditor({
          requiredFeatureIds: context.bindingOptions.requiredLabelFeatureIds,
          optionalFeatureIds: context.bindingOptions.optionalLabelFeatureIds,
          reportedLabelBinding: context.bindingOptions.reportedLabelBinding
        });
      }
    },
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
        if (shouldBindComposition(context)) legendLayout.applyCanvasPadding();
      }
      // The Legend laid out at mount is committed with the bindings the
      // steps since added, so the Result's content is its mounted SVG (R1).
      if (legendLaidOut) previewRuntime.commitActiveResultEdit('legend-position');
    },
    reconcileSelection(context) {
      if (
        !context.bindingOptions.trustedRestore
        && !context.bindingOptions.isIncrementalEdit
      ) {
        featureSelection.clearFeatureSelection({ clearStatus: true, syncDom: false });
      }
    },
    afterReady(context) {
      mountedResultRevision.value += 1;
      // OV-1280: a Result a mode switch shows again, whose drawing changed its
      // rules while it was not shown (a file card's Remove, an Undo or Redo
      // there), is drawn again when its Legend cannot show them, as after a
      // History restore of the rules. OV-2293: any Result displayed is drawn
      // again when its bytes hold a label edit that left its drawing while it
      // was not shown (a Reset Settings, an Undo). The display records what it
      // shows, so a failed rerender does not repeat here.
      const display = displayAsks;
      displayAsks = null;
      if (display?.identity === context.resultIdentity) {
        if (context.bindingOptions.modeArrival && display.rules) featureActions.followRestoredSpecificRules(display.rules);
        if (display.labelIntentGone) featureActions.requestAutomaticRerender();
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
        && resultPaintRecord.drewPaint(context.bindingOptions.drawnPaint, currentEditorProjectionState(drawing));
      previewFeatureSearch.handleMountedResultReady();
      showEditsMadeWhileDrawn(context);
      if (!drewCurrent) return;
      featureActions.redrawShown(drawing);
      history.keepRedrawnResult().catch((error) => {
        console.error('History could not keep the redrawn figure.', normalizeUserFacingError(error));
      });
    }
  });
  // DESIGN A4: the record discovery owner's lists, by mode; a Linear card
  // whose list is still deferred after a Load is read here. Callers read them
  // through the runner's `readDrawingRecordIdentities`, which reads a deferred
  // Circular list first (OV-2489).
  /** @param {DrawingState} drawing @returns {Promise<RecordIdentitySources>} */
  const readRecordIdentitySources = async (drawing) => {
    if (drawing.mode === 'circular') return { mode: drawing.mode, circular: circularRecordEntries(drawing) };
    await awaitLinearRecordDiscovery(drawing);
    return {
      mode: drawing.mode,
      linear: new Map(drawing.linearSeqs.map((seq) => [String(seq.uid), buildDisambiguatedRecordEntries(linearRecordSelector.recordsFor(drawing, seq))]))
    };
  };
  const {
    runLosatJobsInOverlay,
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
    downloadLosatCache,
    downloadLosatPair: downloadLosatPairOfDrawing,
    setLosatPairFilename: setLosatPairFilenameOfDrawing,
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
    isCurrentFeature: isCurrentResultFeature, encodeSessionResources,
    serializeCanonicalFiles: (comparisonPlanSnapshot, linearRecordCatalog, drawing) => (
      serializeActiveRenderFiles(state, drawing, {
        comparisonPlan: comparisonPlanSnapshot,
        linearRecordCatalog,
        hiddenDepth: true
      })
    ),
    prepareLinearRecordCatalog,
    readRecordIdentitySources,
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
    // A Generate that adopts a Circular LOSAT result reconciles its drawing's
    // rings; a Linear Generate and a Circular one without such a result leave
    // them. The ring rows move to its keys only if the draft's rings are those
    // the request was built from (`keepCircularConservationRowsThrough` states why).
    adoptCanonicalRenderArtifacts: (canonical, { builtFromDraft = false, ...options }) => {
      const { drawing } = options;
      const adopt = () => adoptCanonicalRenderArtifacts(canonical, options);
      const adopted = builtFromDraft
        ? circularTrackSlotEditor.keepCircularConservationRowsThrough(drawing, adopt)
        : adopt();
      if (adopted) reconcileCircularConservation(drawing);
    },
    setCircularConservationEnabled,
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
    setGeneratedArtifactIdentity: historySnapshots.setGeneratedArtifactIdentity,
    runGeneratedArtifactReplacement: (...args) => (
      history.runUndoableArtifactReplacement(...args)
    ),
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
  // D-15-89: the panel actions (the template, `window.__GBDRAW_APP__`, the
  // record editor, and the file watcher) act for the shown drawing; this
  // binding resolves it once per action and the runner reads the one it is given.
  /** @param {{ suppress?: boolean, automatic?: boolean } | undefined} [options] */
  const refreshCircularRecordOrder = (options) => refreshRecordOrderOfDrawing(options, state.activeDrawing());
  // R10: the record discovery of one drawing. Each transition that changes
  // the drawing's discovery inputs calls it for that drawing: a source file,
  // the input type, the Linear cards, the end of a Load or a Save, and a
  // History restore of inputs. A drawing switch calls nothing. It is the one
  // admission of a discovery: none starts while a Load rolls back or a Load
  // or History restore installs (the semantic file watchers are suppressed);
  // a Save or Load still pending makes the owners return busy, and its end
  // calls the discovery again. It writes only `drawing`.
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
  /** @param {string} edgeKey @param {unknown} customName */
  const setLosatPairFilename = (edgeKey, customName) => setLosatPairFilenameOfDrawing(state.activeDrawing(), edgeKey, customName);
  /** @param {string} edgeKey */
  const getLosatPairDefaultName = (edgeKey) => losatPairDefaultNameOfDrawing(state.activeDrawing(), edgeKey);
  // R1 (c), PD-OI-102: a History step, a checkpoint (Reset Settings, a rule
  // preset), and their Undo and Redo that change live setting rows of the
  // shown drawing ask the automatic rerender to apply them (the router passes
  // the shown drawing's steps only).
  // Figure size (PD-OI-114 rev 2, Q-D1 A): a step that changes what Python
  // draws of the shown drawing's Legend rows or its Canvas padding (its
  // commit, Undo, or Redo) shows at once, and a figure drawn to a requested
  // size is then fitted again by the automatic rerender, which draws them
  // (`projectCommittedEditorIntent`). An Auto figure is not redrawn. A
  // committed checkpoint step's figure is then not the one its state draws,
  // so History keeps the rerender's figure as the step's, and the restore of a
  // checkpoint figure History knows is not drawn fits it again (D-24).
  /**
   * @param {readonly HistoryChangeLike[] | null} changes
   * @param {'undo' | 'redo' | undefined} direction
   * @param {boolean | undefined} undrawn
   */
  const refitSizedFigure = (changes, direction, undrawn) => {
    const root = svgContainer.value?.querySelector?.('svg');
    if (!root || !requestedCanvasOf(root)) return;
    if (!changes) {
      if (undrawn) featureActions.requestAutomaticRerender();
      return;
    }
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
      const popupFeature = clickedFeature.value?.feat;
      if (!popupFeature) return;
      const popupRecordKey = String(
        popupFeature.record_key ?? popupFeature.recordKey ?? ''
      );
      const popupFeatureId = String(
        popupFeature.biological_feature_id ?? popupFeature.biologicalFeatureId ?? ''
      );
      if (popupRecordKey !== identity.recordKey
        || popupFeatureId !== identity.biologicalFeatureId) return;
      clickedFeature.value.feat = feature;
      const renderedId = String(
        feature?.rendered_feature_svg_id
        ?? feature?.renderedFeatureSvgId
        ?? feature?.svg_id
        ?? ''
      );
      if (renderedId) clickedFeature.value.svg_id = renderedId;
    }
  });
  const recordActionsExpanded = ref(false);
  watch(clickedFeature, (/** @type {{ feat?: Record<string, unknown> } | null} */ popup) => {
    recordActionsExpanded.value = false;
    if (popup?.feat) {
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
    rulePreparation,
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
    // The phone layout shows a loaded Result in the Preview, as Generate does; a Load
    // without a Result shows Settings.
    if (result?.status === 'ok' || result?.status === 'legacy') settingsPaneActions.showPhonePane(svgContent.value ? 'preview' : 'settings');
    // R17-25: a file Session Load with a Result fits the figure to the shown Preview, at every
    // width, after the phone pane switch and the whole notice (its height sizes the canvas), so
    // the measured canvas is the one shown; Generate resets to 100 %. A legacy configuration
    // Load brings no Result, so the figure already shown keeps its zoom.
    if (result?.status === 'ok' && svgContent.value) {
      await nextTick();
      fitPreviewToViewport();
    }
    // A loaded Session makes a Generate likely: download the Worker's boot
    // files now, after the Load, so they never slow its own requests.
    if (result?.status === 'ok' || result?.status === 'legacy') prefetchPyodideBootAssets();
    return result;
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
      const selectedResult = importedResults[resultIndex] || null;
      if (!selectedResult) return null;
      const token = `session-load:${nextSessionPreviewToken++}`;
      return previewRuntime.registerReadinessExpectation(/** @type {any} */ ({
        result: selectedResult,
        resultIndex,
        artifactIdentity: token,
        generationToken: token,
        catalogState: state.featureCatalog?.value || null,
        phase: 'session-load',
        bindingOptions: { isIncrementalEdit: true },
        isCurrent: () => (
          results.value[resultIndex] === selectedResult
          && Number(selectedResultIndex.value) === resultIndex
        )
      }));
    },
    rollbackState: sessionImportRollbackState,
    afterImport: async (/** @type {{ status: string, decompressedCharacters?: number, isCurrent: () => boolean }} */ result) => {
      // Only a Session Load comes here; a legacy configuration Load returns
      // before it (OV-1792).
      if (result?.status === 'ok') {
        // The comparison rings read from the loaded files and series (an older
        // Session's series names included), before the History baseline.
        for (const drawing of circularDrawings()) reconcileCircularConservation(drawing);
        // The loaded draft's live values wait for their next change or Generate (6.7).
        discardLiveSettingRequests();
        historySnapshots.clearGeneratedArtifactIdentity({
          retainedBytes: Number(result.decompressedCharacters || 0) * 2
        });
        // A replaced document resets transient UI: the selection named features
        // of the previous Session.
        featureSelection.clearFeatureSelection({ clearStatus: true, syncDom: false });
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
      featureActions.requestAutomaticRerender();
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
  const restoreLegendItems = (label, restore) => history.runUndoable(label, () => restoreLegendRows(restore));

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

  // A displayed Result whose label intent changed since it was last shown
  // receives the label projection in the binder's label step, also when no
  // label intent remains (an undone or replaced Label TSV import).
  let labelProjectionResultIdentity = '';
  /** @param {DrawingState} drawing @returns {EditorPaintState} */
  const currentEditorProjectionState = (drawing) => editorPaintState(state, drawing);
  /** @param {DrawingState} drawing @param {MountedResultContext} context */
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
    resultPaintRecord.commit(identities, context.resultIdentity, current, {
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
    const departed = resultPaintRecord.depart(currentEditorProjectionState(drawing));
    if (departed) departedResultIntent.set(departed, displayedIntentSignature(drawing));
  };
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
          toRaw(extractedFeatures.value) || [], results.value.map((/** @type {{ name?: unknown } | null | undefined} */ result) => result?.name), resultIndex, svg
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
      const resultIndex = previewRuntime.getActiveRuntime()?.resultIndex ?? (Number(selectedResultIndex.value) || 0);
      try {
        shown = compileDisplayedResultOperations(state.activeDrawing(), resultIndex, { domains, ruleRows });
      } catch (error) {
        console.error('Editor edits could not be compiled for the displayed Result.', normalizeUserFacingError(error));
        return false;
      }
    }
    return previewRuntime.applyEditorOperations(shown, {
      domains,
      // The Legend is laid out without its own commit: the display commits
      // the Result once (S4).
      afterApply: (_svg, reconciled) => {
        if (!reconciled) return;
        legendActions.readLegendRowNames();
        legendActions.onLegendGeometryChanged({ commit: false });
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
  ) => history.runUndoable(label, () => showingEditorIntent(mutate, domains)(...args));
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
  // What a displayed Result asks of the rerender, read once it is ready
  // (`afterReady`): `rules`, the rule table it was last shown with, or the one
  // it was drawn with when that is unknown (a restored Result), for OV-1280;
  // `labelIntentGone`, whether its bytes hold a label edit that left its drawing
  // since it was last shown (OV-2293).
  /** @type {{ identity: string, rules: SpecificColorRule[] | null, labelIntentGone: boolean } | null} */
  let displayAsks = null;
  // Whether a Result draws a feature; an unknown catalog answers yes.
  /** @param {DrawingState} drawing @param {number} resultIndex @returns {(recordKey: string, featureId: string) => boolean} */
  const resultDrawsFeature = (drawing, resultIndex) => {
    const rendered = resultCatalogFeatures(state, drawing.mode, resultIndex)?.renderedByIdentity;
    return (recordKey, featureId) => !rendered || rendered.has(biologicalFeatureKey(recordKey, featureId));
  };
  // D-07 (PD-OI-062): a batch Result shows the canonical color, visibility,
  // Legend, and label edits when it is displayed. Labels follow in the
  // binder's label step.
  /**
   * @param {DrawingState} drawing
   * @param {any} context
   * @param {string[]} resultLegendOrder
   */
  const projectEditorIntentOnDisplay = async (drawing, context, resultLegendOrder) => {
    const identity = context.resultIdentity;
    const current = currentEditorProjectionState(drawing);
    displayAsks = null;
    const departedIntent = departedResultIntent.get(identity);
    departedResultIntent.delete(identity);
    if (departedIntent !== undefined && departedIntent === displayedIntentSignature(drawing)) {
      resultPaintRecord.rebind(identity);
      labelProjectionResultIdentity = '';
      recordStructuralMetric('displayedResultEditorProjectionCount', 0, {
        phase: context.phase,
        rootGeneration: context.rootGeneration
      });
      return;
    }
    const previous = resultPaintRecord.display(identity, current);
    const drawn = resultPaintRecord.drawn(identity) ?? current;
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
    /** @type {ReturnType<typeof compileDisplayedResultOperations>} */
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
    resultPaintRecord.shown(identity, current, previous, { declined, failed });
    if (declined && !failed && displayBlocked()) await showDeclinedPaint(identity);
  };
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
    const lacking = resultPaintRecord.lacking(identity);
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
      if (resultPaintRecord.lacking(identity) === lacking) resultPaintRecord.shown(identity, current, current);
    } catch (error) {
      console.error('Editor edits could not be shown on the displayed Result.', normalizeUserFacingError(error));
    }
  };

  // A label rerender's displayed Result, once bound, shows the fills and
  // visibility of the edits made while it ran (OV-346).
  /** @param {{ resultIdentity: string, bindingOptions: Readonly<import('./preview-runtime.js').PreviewBindingOptions> }} context */
  const showEditsMadeWhileDrawn = (context) => {
    if (context.bindingOptions.drawnPaint && resultPaintRecord.lacking(context.resultIdentity)) {
      void showDeclinedPaint(context.resultIdentity);
    }
  };

  // E1: shows the artifact of drawing `nextId`. Keeps the displayed one in
  // the departing drawing's slot and installs the arriving slot, or an empty
  // one. The Legend owner
  // keeps each Result's inventory (OV-47): the arriving Result's goes to it,
  // and the binder adopts it when the Result is displayed, as for a batch
  // Result.
  /** @param {Readonly<ArtifactSlot>} slot An installed slot. */
  const rememberShownResultInventory = (slot) => {
    const shown = results.value[Number(selectedResultIndex.value) || 0];
    if (shown) legendActions.rememberResultInventory(previewRuntime.getResultIdentity(shown), [...slot.legendInventory]);
  };
  /**
   * The displayed artifact is kept under the drawing it belongs to, which is
   * the arriving drawing when the shown drawing displayed the arriving one's.
   * @param {DrawingIdentity} arrivingDrawing
   * @returns {Readonly<ArtifactSlot> | null} The installed slot (null: the drawing has no artifact yet).
   */
  const swapArtifactSlots = (arrivingDrawing) => {
    previewRuntime.clearActiveRuntime();
    const arriving = artifactSlots.swap(historySnapshots.captureArtifactSlot(), arrivingDrawing.id);
    historySnapshots.installArtifactSlot(arriving, { drawing: arrivingDrawing });
    if (arriving) rememberShownResultInventory(arriving);
    return arriving;
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
    || history.historyAvailability()
    || (history.restoring.value || history.traversalPending() ? HISTORY_RESTORE_BUSY : null)
    || (featureActions.ruleWritePending() ? RULE_WRITE_BUSY : null);
  const drawingSwitchAvailable = computed(() => !drawingSwitchBusy());
  // The drawing an operation error names (PD-OI-079 `context.drawingId`), if
  // any: read by the switch, the drawing removal, and E1's Generate. Its
  // producers write a `DrawingId`; the normalizer, a leaf, keeps it as text.
  const errorDrawingId = () => /** @type {DrawingId | undefined} */ (normalizeUserFacingError(errorLog.value)?.context.drawingId);
  // The switch work, steps 1-5, without the gate: `selectDrawing` and a step
  // of the drawing list (`app/drawing-list.js`, which shows a drawing inside
  // its project step, R16-2) call it. Showing the shown drawing does nothing.
  /** @param {DrawingId} drawingId */
  const showDrawingNow = (drawingId) => {
    const arriving = state.drawingOf(drawingId);
    const departing = state.activeDrawing();
    if (arriving === departing) return;
    // 1. Settle the departing drawing; its capture runs in this call.
    void history.settlePendingIntent();
    // A card's failed import is its drawing's card state (D-15-126): read
    // before the departing drawing's message is cleared below.
    const arrivingCardFailure = auxiliaryImportFailureOnArrival(departing, arriving);
    featureActions.suspendSpecificRulePatternDrafts();
    featureSelection.clearFeatureSelection({ clearStatus: true, syncDom: false });
    rememberDepartingResultProjection(departing);
    // 2. Artifact slots.
    swapArtifactSlots(arriving);
    // 3. The shown drawing.
    installShownDrawing(arriving.id);
    void history.refreshCurrentIntent();
    // 4. Transient UI. "Showing the last successful result" named the
    // departing drawing's Result, and the notice the departing drawing's work.
    history.adoptSessionReset(resetModeTransientUi);
    noticeBanner.dismiss();
    previewFeatureSearch.switchDrawingSearch(departing.id, arriving.id);
    failedGeneratePreservedResult.value = false;
    // OV-372: an operation error belongs to the drawing it was raised in,
    // except one that names the arriving drawing's Result (PD-OI-079).
    if (errorLog.value && errorDrawingId() !== arriving.id) errorLog.value = null;
    if (arrivingCardFailure !== null && errorLog.value === null) errorLog.value = arrivingCardFailure;
    // The arriving drawing's annotations, over its mode's record catalog.
    rebindAnnotationRecords();
    // 5. Presentation: the arriving Result is shown as a selection.
    previewRuntime.presentSelectedResult({ modeArrival: true });
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
    // The copy's Result is its source's slot under the copy's ID: the shown
    // drawing's through `captureArtifactSlot` (it reads `currentResults`, R1).
    copyDrawing: (sourceId, copyId) => {
      const drawing = duplicateDrawingState(state.drawingOf(sourceId), copyId);
      const slot = artifactSlots.readOfDrawing(sourceId, state.activeDrawingId.value,
        () => historySnapshots.captureArtifactSlot(), (kept) => kept);
      return { drawing, slot: slot ? artifactSlotAs(slot, drawing) : null };
    },
    // A region drawing is read as Load reads a stored drawing, then made under its new ID.
    regionDrawing: makeRegionDrawing,
    addDrawing: state.addDrawing,
    // An error that names a drawing leaves with it (R17-29), as it leaves on
    // a switch away (OV-372): every Delete, its Redo, and the Undo of New and
    // Duplicate remove a drawing here.
    removeDrawing: (drawingId) => {
      if (errorDrawingId() === drawingId) errorLog.value = null;
      return state.removeDrawing(drawingId, artifactSlots);
    },
    putSlot: artifactSlots.putSlot,
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
  // committed one; each `stashed` slot waits in its drawing. The stash starts
  // again with the loaded Session.
  /** @param {LoadedArtifactSlots} slots */
  const installLoadedArtifactSlots = ({ opening, stashed }) => {
    artifactSlots.replaceLoaded(stashed);
    if (!opening) return;
    historySnapshots.installArtifactSlot(opening, { drawing: state.drawingOf(opening.drawingId) });
    rememberShownResultInventory(opening);
  };

  const { requestPaletteChange, selectPalette, requestResetColors } = resultsManager;
  /**
   * @template {unknown[]} A
   * @template R
   * @param {string} label
   * @param {(...args: A) => R} fn
   * @param {{ joinsOpenStep?: boolean }} [options] History's options for the step.
   */
  const undoableAction = (label, fn, options) => (/** @type {A} */ ...args) => history.runUndoable(label, () => fn(...args), options);
  // Python's report of the per-feature edits the Results do not draw (design
  // Q4 3.4): dormant edits outside the crop or display, edits whose feature
  // the source does not have, and the edits a source-replacing Generate removed.
  const featureIdentityNoticeSummary = computed(() => {
    const drawing = state.activeDrawing();
    const notices = Array.isArray(featureIdentityNotices.value) ? featureIdentityNotices.value : [];
    const summary = {
      dormant: notices.filter((/** @type {UnresolvedNotice} */ notice) => notice.status !== 'unresolved').length,
      unmatched: countUnresolvedFeatureEdits({ featureEdits: drawing.slices.featureEdits, notices }),
      removed: Number(featureEditRemovalCount.value) || 0
    };
    return summary.dormant || summary.unmatched || summary.removed ? summary : null;
  });
  // An explicit removal of the edits Python reported unresolved (R2).
  const removeUnmatchedFeatureEdits = undoableAction(msg('Remove unmatched feature edits'), () => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    return removeUnresolvedFeatureEdits({
      featureEdits: drawing.slices.featureEdits, notices: featureIdentityNotices.value
    }) > 0;
  });
  const addFeatureVisibilityRuleWithHistory = undoableAction(msg('Add feature visibility rule'), addFeatureVisibilityRule);
  const moveFeatureVisibilityRuleDownWithHistory = undoableAction(msg('Move feature visibility rule'), moveFeatureVisibilityRuleDown);
  const moveFeatureVisibilityRuleUpWithHistory = undoableAction(msg('Move feature visibility rule'), moveFeatureVisibilityRuleUp);
  const removeFeatureVisibilityRuleWithHistory = undoableAction(msg('Remove feature visibility rule'), removeFeatureVisibilityRule);
  const setFeatureVisibilityRuleFieldWithHistory = undoableAction(
    msg('Edit feature visibility rule'),
    setFeatureVisibilityRuleField
  );
  const setFeatureVisibilityWithHistory = undoableAction(msg('Change feature visibility'), setFeatureVisibility);
  const updateClickedFeatureVisibilityWithHistory = undoableAction(
    msg('Change feature visibility'),
    updateClickedFeatureVisibility
  );
  const handleFeatureVisibilityScopeChoiceWithHistory = undoableAction(
    msg('Change feature visibility'),
    handleFeatureVisibilityScopeChoice
  );
  const requestFeatureColorChangeWithHistory = undoableAction(msg('Change feature color'), requestFeatureColorChange);
  const setFeatureColorValueWithHistory = undoableAction(msg('Change feature color'), setFeatureColorValue);
  const updateClickedFeatureColorWithHistory = undoableAction(msg('Change feature color'), updateClickedFeatureColor);
  const editLegendEntryColor = editEditorIntent(msg('Change legend color'), updateLegendEntryColor, { domains: LEGEND_FILL_DOMAINS });
  /** @param {number} index @param {string} color */
  const updateLegendEntryColorShown = (index, color) => (
    legendRowHasRules(index) ? updateLegendEntryColor(index, color) : editLegendEntryColor(index, color)
  );
  const setTypeDefaultColorWithHistory = undoableAction(msg('Change feature color'), featureActions.setTypeDefaultColor);
  // OV-1281: the Legend item popup's Color on a row a type's default color draws sets that color (OIC-029).
  // `input` shows the row's color again afterwards (`shownColor()`): a refused edit leaves the view as it was.
  /** @param {number} index @param {string} color @param {HTMLInputElement} input @param {() => string} shownColor */
  const updateLegendEntryColorInPopup = async (index, color, input, shownColor) => {
    try {
      return await requestLegendEntryColorInPopup(index, color);
    } finally {
      await nextTick();
      if (input.isConnected && input.value !== shownColor()) input.value = shownColor();
    }
  };
  /** @param {number} index @param {string} color */
  const requestLegendEntryColorInPopup = async (index, color) => {
    const entry = state.activeDrawing().legendEntries.value[index];
    const type = entry ? await featureActions.preparedLegendRowDefaultColorType(legendEntryKey(entry)) : null;
    if (typeof type === 'string') return setTypeDefaultColorWithHistory(type, color);
    // null: no default color draws the row (the swatch only); undefined: the matches went stale, so do nothing.
    return type === null ? updateLegendEntryColorShown(index, color) : type;
  };
  const setLegendEntryStrokeColorValueWithHistory = editEditorIntent(
    msg('Change legend stroke color'),
    setLegendEntryStrokeColorValue
  );
  // OV-161, D-12: a dialog choice is one History step, Cancel records none, and
  // neither starts while History's step is open (`createDialogChoice`).
  /**
   * @param {() => string} label
   * @param {(choice: string, ...rest: any[]) => any} handler
   * @param {() => void} [cancel]
   */
  const scopeChoiceWithHistory = (label, handler, cancel = cancelFeatureStyleScope) => (
    dialogChoice.withHistory(label, handler, cancel)
  );
  const handleColorScopeChoiceWithHistory = scopeChoiceWithHistory(() => msg('Change feature color'), handleColorScopeChoice);
  const handleFeatureColorScopeChoiceWithHistory = scopeChoiceWithHistory(() => msg('Change feature color'), handleFeatureStyleScopeChoice);
  const handleFeatureStrokeScopeChoiceWithHistory = scopeChoiceWithHistory(
    () => msg('Change feature stroke'),
    showingEditorIntent(handleFeatureStyleScopeChoice, STROKE_DOMAINS)
  );
  /** @param {string} choice @param {unknown[]} rest */
  const handleFeatureStyleScopeChoiceWithHistory = (choice, ...rest) => (
    featureStyleScopeDialog.kind === 'stroke'
      ? handleFeatureStrokeScopeChoiceWithHistory(choice, ...rest)
      : handleFeatureColorScopeChoiceWithHistory(choice, ...rest)
  );
  // The label scopes the Color/Stroke Change Scope dialog offers; a
  // label that repeats the Legend name, or the other label, is not offered.
  // OV-648: the color owner decides which label scopes select other features.
  const featureStyleScopeLabelOptions = computed(() => ({
    displayLabel: featureStyleScopeDialog.displayLabelOffered === true,
    annotationLabel: featureStyleScopeDialog.annotationLabelOffered === true
  }));
  const handleLegendNameCommitWithHistory = undoableAction(msg('Rename legend item'), handleLegendNameCommit);
  const joinLegendGroupWithHistory = undoableAction(msg('Rename legend item'), joinLegendGroup);
  // The Legend tab's grip reorders a row by drag or by the Up and Down arrow
  // keys; neither is a click, so each move records its own History step.
  const moveLegendEntryWithHistory = undoableAction(msg('Move legend item'), moveLegendEntry);
  const legendReorder = reactive({ from: /** @type {number | null} */ (null), over: /** @type {number | null} */ (null) });
  /** @param {number} from @param {number} to @param {boolean} [refocus] */
  const reorderLegendEntry = (from, to, refocus = false) => {
    const count = state.activeDrawing().legendEntries.value.length;
    if (from === to || to < 0 || to >= count) return;
    moveLegendEntryWithHistory(from, to);
    if (refocus) {
      nextTick(() => {
        const handle = document.querySelector(`[data-legend-reorder-handle="${to}"]`);
        if (handle instanceof HTMLElement) handle.focus();
      });
    }
  };
  /** @param {KeyboardEvent} event @param {number} index */
  const handleLegendReorderKey = (event, index) => {
    const step = event.key === 'ArrowUp' ? -1 : (event.key === 'ArrowDown' ? 1 : 0);
    if (!step) return;
    event.preventDefault();
    reorderLegendEntry(index, index + step, true);
  };
  /** @param {DragEvent} event @param {number} index */
  const startLegendReorder = (event, index) => {
    legendReorder.from = index;
    legendReorder.over = index;
    const row = event.target instanceof Element ? event.target.closest('[data-legend-drawer-row]') : null;
    if (event.dataTransfer) {
      event.dataTransfer.effectAllowed = 'move';
      event.dataTransfer.setData('text/plain', String(index));
      if (row) event.dataTransfer.setDragImage(row, 12, 12);
    }
  };
  /** @param {DragEvent} event @param {number} index */
  const overLegendReorder = (event, index) => {
    if (legendReorder.from === null) return;
    event.preventDefault();
    legendReorder.over = index;
  };
  /** @param {DragEvent} event @param {number} index */
  const dropLegendReorder = (event, index) => {
    const from = legendReorder.from;
    if (from === null) return;
    event.preventDefault();
    legendReorder.from = null;
    legendReorder.over = null;
    reorderLegendEntry(from, index);
  };
  const endLegendReorder = () => {
    legendReorder.from = null;
    legendReorder.over = null;
  };
  // R15-3 (OV-285): the Merge of a rename onto a deleted row's caption is the
  // Restore of that row (by its Python key), then the merge, in one History
  // checkpoint: the merge's rule commit joins it, so one Undo returns both.
  // Keep both onto a deleted row in the renamed row's color restores it the
  // same way: the two are one row (D-28).
  /** @param {string} choice */
  const chooseLegendRename = (choice) => {
    const key = legendRenameDialog.deletedTargetKey;
    const restores = choice === 'merge' || (choice === 'keep-both' && legendRenameDialog.keepBothRestores);
    if (!restores || !key) return handleLegendRenameChoice(choice);
    return history.runUndoableCheckpoint(msg('Rename legend item'), async () => {
      const deleted = state.activeDrawing().deletedLegendEntries.value;
      const index = deleted.findIndex((/** @type {KeptLegendRow} */ entry) => (entry.originalCaption || entry.caption) === key);
      if (index >= 0) restoreLegendRows(() => restoreDeletedLegendEntries([index]));
      return handleLegendRenameChoice(choice);
    });
  };
  const handleLegendRenameChoiceWithHistory = scopeChoiceWithHistory(
    () => msg('Rename legend item'),
    chooseLegendRename,
    cancelLegendRename
  );
  // D-15: the palette dialog's choice is one History step too.
  const handlePaletteColorsChoiceWithHistory = scopeChoiceWithHistory(
    () => msg('Change setting'),
    resultsManager.handlePaletteColorsChoice,
    resultsManager.cancelPaletteColorsDialog
  );
  const handleResetColorChoiceWithHistory = scopeChoiceWithHistory(
    () => msg('Reset feature color'),
    handleResetColorChoice,
    cancelResetColor
  );
  const resetClickedFeatureFillColorWithHistory = undoableAction(msg('Reset feature color'), resetClickedFeatureFillColor);
  const updateClickedFeatureStrokeWithHistory = editEditorIntent(msg('Change feature stroke'), updateClickedFeatureStroke);
  const setClickedFeatureStrokeColorValueWithHistory = editEditorIntent(
    msg('Change feature stroke'),
    setClickedFeatureStrokeColorValue
  );
  const setClickedFeatureStrokeWidthValueWithHistory = editEditorIntent(
    msg('Change feature stroke'),
    setClickedFeatureStrokeWidthValue
  );
  const resetClickedFeatureStrokeWithHistory = editEditorIntent(msg('Reset feature stroke'), resetClickedFeatureStroke);
  const setFeatureColorWithHistory = undoableAction(msg('Change feature color'), setFeatureColor);
  const selectedFeatureBulkColor = ref('#2563eb');
  const selectedFeatureBulkCaption = ref('Selected features');
  const selectedFeatureBulkVisibility = ref('off');
  const selectedFeatureBulkStrokeColor = ref('#1f2937');
  const selectedFeatureBulkStrokeWidth = ref(1.5);
  const applySelectedFeatureColor = () => history.runUndoable(msg('Change selected feature color'), async () => {
    const changed = await applyColorToSelectedFeatures(
      selectedFeatures.value,
      selectedFeatureBulkColor.value,
      selectedFeatureBulkCaption.value
    );
    if (changed) featureSelection.syncFeatureSelectionClasses();
    return changed;
  });
  const applySelectedFeatureVisibility = async () => {
    const changed = await history.runUndoableCommand(msg('Change selected feature visibility'), () =>
      buildSelectedFeaturesVisibilityCommand(selectedFeatures.value, selectedFeatureBulkVisibility.value)
    );
    if (changed) featureSelection.clearFeatureSelection({ clearStatus: true });
    return changed;
  };
  const editSelectedFeatureStroke = editEditorIntent(msg('Change selected feature stroke'), () => applyStrokeToSelectedFeatures(
    selectedFeatures.value,
    selectedFeatureBulkStrokeColor.value,
    selectedFeatureBulkStrokeWidth.value
  ));
  const applySelectedFeatureStroke = async () => {
    const changed = await editSelectedFeatureStroke();
    if (changed === true) featureSelection.syncFeatureSelectionClasses();
    return changed;
  };
  const openFirstSelectedFeature = (event = null) => {
    const first = selectedFeatures.value[0] || null;
    if (!first) return null;
    return openFeatureEditorForFeature(first, event);
  };
  const updateClickedFeatureLabelTextWithHistory = undoableAction(msg('Change label text'), updateClickedFeatureLabelText);
  const handleLabelTextScopeChoiceWithHistory = undoableAction(msg('Change label text'), handleLabelTextScopeChoice);
  // The Label Not Shown choice answers the label edit whose step is open: one step, under that step's name.
  // A registered test reads the label from this line (feature-label-visual-unit.test.mjs).
  const handleHiddenLabelTextChoiceWithHistory = undoableAction(msg('Change label visibility'), handleHiddenLabelTextChoice, { joinsOpenStep: true });
  const requestLabelTextChangeByFeatureIdWithHistory = undoableAction(
    msg('Change label text'),
    requestLabelTextChangeByFeatureId
  );
  const requestLabelTextChangeByKeyWithHistory = undoableAction(msg('Change label text'), requestLabelTextChangeByKey);
  // The Features list's label input commits what was typed; its draft then follows the label.
  /** @param {string} featureId @param {string} text */
  const commitListLabelDraft = (featureId, text) => {
    const committed = requestLabelTextChangeByFeatureIdWithHistory(featureId, text);
    const entry = featureActions.getEditableLabelByFeatureId(featureId);
    if (entry) entry.draftText = entry.text;
    return committed;
  };
  const resetAllLabelTextOverridesWithHistory = undoableAction(msg('Reset label edits'), resetAllLabelTextOverrides);
  const loadLabelOverrideTableWithHistory = undoableAction(msg('Load label edits'), loadLabelOverrideTable);
  const {
    runInfoCopyStatus, exactReplayCopyStatus, runInfoElapsedText, runInfoReproducibilityText, runInfoHasCliHelperFiles,
    copyRunCommand, copyExactReplayCommand
  } = createRunInfoCopy({ lastRunInfo });

  /** @param {Pick<AnnotationRecordCatalog, 'issues'>} catalog */
  const catalogIssueError = (catalog) => {
    const issue = catalog.issues[0];
    return issue ? diagnosticError(issue.code, issue.context) : diagnosticError('INPUT_UNREADABLE');
  };
  /**
   * Why a drawing's record catalog is not ready, as its record list says it: for Linear, the first drawn card's own
   * read error (its record select's), else the catalog's first issue. Generate and Add to annotation set read it here.
   * @param {DrawingState} drawing
   * @param {Pick<AnnotationRecordCatalog, 'issues'>} catalog
   */
  const recordCatalogError = (drawing, catalog) => (drawing.mode === 'linear'
    // A stored '' (no error) has no `code`: the read is undefined, as for null.
    && drawing.drawnLinearSeqs.value.map(seq => linearRecordSelector.errorModelFor(drawing, seq)).find(error => /** @type {UserFacingError} */ (error)?.code))
    || catalogIssueError(catalog);
  async function prepareLinearRecordCatalog({ privateCandidate = false } = {}) {
    const drawing = state.activeDrawing();
    if (drawing.mode !== 'linear') return { catalog: null, error: '' };
    const hasAutomaticSequence = drawing.linearSeqs.some((sequence) => {
      if (String(sequence?.region_record_id || '').trim()) return false;
      return inputSourceFiles(drawing.inputType.value, sequence).complete;
    });
    const hasRegionAnnotations = drawing.annotationSets.some((/** @type {{ annotations?: unknown } | null | undefined} */ set) => (
      Array.isArray(set?.annotations) && set.annotations.length > 0
    ));
    if (!hasAutomaticSequence && !hasRegionAnnotations) {
      return { catalog: null, error: '' };
    }
    if (privateCandidate) {
      const linearRecordsOff = new Set(drawing.recordsOff);
      const sources = await Promise.all(drawing.linearSeqs.map(async (seq) => {
        const drawn = !linearRecordsOff.has(seq.uid);
        const inputType = drawing.inputType.value;
        const { primaryFile, pairedFile } = inputSourceFiles(inputType, seq);
        const sourceKey = annotationSourceKey({ scope: 'linear', uid: seq.uid, inputType, primaryFile, pairedFile });
        try {
          const records = await readInputRecords(inputType, primaryFile, pairedFile);
          return { sourceKey, selector: seq.region_record_id, drawn, hasInput: Boolean(primaryFile), status: 'ready', records };
        } catch (error) {
          /** @type {UserFacingError | ''} The model the card's record select would store for this failure. */
          const model = normalizeUserFacingError(error, { operation: 'listSequenceRecords', stage: 'helper' }) || '';
          return { sourceKey, selector: seq.region_record_id, drawn, hasInput: Boolean(primaryFile), status: 'error', error: model, records: [] };
        }
      }));
      const catalog = getAnnotationRecordCatalog(drawing, sources);
      return catalog.status === 'ready' ? { catalog, error: '' }
        : { catalog: null, error: catalogIssueError(catalog) };
    }
    let catalog = getAnnotationRecordCatalog(drawing);
    if (catalog.status !== 'ready') {
      try {
        await linearRecordSelector.refresh(drawing);
      } catch (error) {
        return { catalog: null, error: normalizeUserFacingError(error, { operation: 'listSequenceRecords', stage: 'helper' }) };
      }
      catalog = getAnnotationRecordCatalog(drawing);
    }
    return catalog.status === 'ready' ? { catalog, error: '' } : { catalog: null, error: recordCatalogError(drawing, catalog) };
  }

  const runAnalysis = () => {
    const drawing = state.activeDrawing();
    // The phone layout shows the Preview, where the run's progress, Result or error appears; a new
    // Result comes in with its workspace at the pane's top while the Preview is still shown.
    settingsPaneActions.showPhonePane('preview', { workspace: false });
    const patternDrafts = featureActions.captureSpecificRulePatternDrafts();
    return runGeneratedDiagramAnalysis(null, null, null, {
      prepareGenerate: async () => {
        if (mode.value === 'linear') {
          if (drawing.importedComparisonIntent.disposition === IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE) {
            await materializeAutomaticLinearRecords(drawing);
          } else {
            await linearRecordSelector.refresh(drawing);
          }
        }
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
          featureSelection.clearFeatureSelection({ clearStatus: true });
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
    getRecordCatalog: () => getAnnotationRecordCatalog(state.activeDrawing()),
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
    deactivateResolvedLinearComparisonLosatFilename, addLinearComparisonBatch, setLinearComparisonGapAction,
    setLinearComparisonCardFile
  } = createLinearComparisonGapEdits({
    sessionOperationAvailability, shownAction, replaceLinearComparisonPlan, syncLinearRecordLayout,
    linearComparisonPairForEdgeKey, omitLinearComparisonPairs, setLinearComparisonFile
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
      // kept ring files, as a new drawing with those files uploaded has them.
      reconcileCircularConservation(drawing);
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
        region_start: null,
        region_end: null,
        region_reverse: false
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
    historySnapshots.writeArtifactValues(drawing.id, drawing.id === state.activeDrawingId.value
      ? artifact
      : { ...artifact, similarityAlignmentResetReceipt: null });
    RESET_FOLLOW_UPS[drawing.mode](drawing);
    // The reset Depth series fit the kept Depth files (OV-1164).
    ensureDepthTrackConfigCount(drawing);
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
    noticeBanner.dismiss();
    settingsPaneActions.closeSettingsPopup();

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
      featureActions.clearSpecificRulePatternDrafts();
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
        legendActions.extractLegendEntries();
      }
      return true;
    })));
  };
  // Reset all drawings is one project step (D-15-51): one Undo restores every
  // drawing and the Session values. OV-1166 (Q-DR-6 A): every drawing's
  // multi-record layout takes a new drawing's defaults (its sides unset),
  // whatever the canvas was before.
  const resetSettings = () => runReset((reset) => history.runProjectStep(msg('Reset settings'), reset), () => {
    resetProject(state);
    state.drawingList().forEach((drawing) => resetDrawing(drawing));
    // As 0.14.0 Reset settings did; Reset this drawing keeps the shown
    // Result's match sequences with the Result (OV-2203).
    matchSequenceRegistry.reset();
  });
  // Reset this drawing is a step of the shown drawing (R17-5): the other
  // drawings, their Histories, and the Session values stay.
  const resetShownDrawing = () => runReset((reset) => history.runUndoable(msg('Reset this drawing'), reset),
    () => resetDrawing(state.activeDrawing()));

  // OV-1160 (Owner R15-23 Q-DR-5 A): turned on, the Instant Preview toggle
  // applies every drawing's queued palette, so it is one project step whose
  // parts hold each drawing's palette (D-15-51); turned off, it changes the
  // Session-level flag alone.
  /** @param {boolean} enabled */
  const setPaletteInstantPreview = (enabled) => (enabled
    ? history.runProjectStep('Change setting', () => resultsManager.setPaletteInstantPreview(true))
    : history.runUndoable('Change setting', () => resultsManager.setPaletteInstantPreview(false)));

  const resetLayout = () => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    resetAllPositions();
    resetCanvasPadding();
    resetLayoutState(state);
    resetPreviewViewport({ resetZoom: true });
  };

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

  /** @param {unknown} value */
  const normalizeSessionTitle = (value) => {
    if (value === null || value === undefined) return '';
    return String(value).trim();
  };

  const errorDisplay = computed(() => normalizeUserFacingError(errorLog.value));
  // D-25: a Result without current feature metadata (a Session older than 40)
  // is saved only after one Generate; the notice comes from the same state.
  const sessionSaveNeedsGenerate = computed(() => results.value.length > 0 && !state.featureCatalog?.value);
  const reloadAfterOperationError = () => window.location.reload();

  const sessionTitleLabel = computed(() => {
    const title = normalizeSessionTitle(sessionTitle.value);
    return title || t('Untitled session');
  });

  const canUseLinearRulerOnAxis = computed(
    () => linearAxisRulerApplies(state.activeDrawing().form)
  );

  const canUseCircularScaleStyling = computed(
    () => {
      const drawing = state.activeDrawing();
      return (
        drawing.adv.circular_track_slots_enabled
          ? drawing.adv.circular_track_slots.some((/** @type {CircularTrackSlot} */ slot) => (
              slot?.renderer === 'ticks' &&
              circularTrackSlotEditor.circularTrackSlotEffectiveEnabled(slot)
            ))
          : drawing.form.show_scale !== false
      );
    }
  );

  const {
    clickedFeatureHeading, clickedFeatureSummary, copyClickedFeatureDetails, copyClickedFeatureQualifiers
  } = createClickedFeatureCopy({ clickedFeature });

  const {
    downloadText, canRetryInteractiveSvgExport, downloadSVG, downloadInteractiveSVG, downloadPNG, downloadPDF
  } = createResultExport({
    state, errorLog, showNotice, captureSvgExport, loadExportService, runDiagramHelperOperation, DIAGRAM_HELPER_OPERATIONS,
    displayedArtifactDrawing: historySnapshots.displayedArtifactDrawing, embedInExports: fontSettings.embedInExports,
    shownInputSummary, shownLinearWrapLeaf, shownCommittedRequest, awaitLinearRecordDiscovery, committedLinearRecordsOf
  });

  // The list under the popup's Legend name: the color-rule Legend groups of
  // the clicked feature's type, by name (PD-OI-061 revision 4). Each is shown
  // as the Legend shows its row: by the row's drawn name, and a group that
  // draws no row by its name as written, never by a key.
  const popupLegendGroups = computed(() => {
    const type = clickedFeature.value?.feat?.type;
    if (type === undefined) return [];
    const drawing = state.activeDrawing();
    const rowKeyOf = groupLegendRowKeys(displayedLegendRowContext(state, drawing));
    return colorRuleLegendGroups(drawing.manualSpecificRules, String(type))
      .map((group) => ({ group, shownName: legendActions.shownLegendName(rowKeyOf(group)) }))
      .sort((left, right) => left.group.name.localeCompare(right.group.name) || left.shownName.localeCompare(right.shownName));
  });

  const featureVisibilityFeatureSuggestions = computed(() => {
    const values = new Set(['*', ...featureKeys]);
    for (const feat of Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []) {
      const type = String(feat?.type || '').trim();
      if (type) values.add(type);
    }
    return Array.from(values);
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
    shownDrawing: state.activeDrawing,
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
    circularSingleRecordPresentationEnabled
  } = createCircularRecordPresentation({
    state, mode, semanticMutationAvailable, sessionOperationAvailability, shownDrawingOfMode, refreshCircularRecordOrder
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
  const { labelFitDialog, openLabelFitDialog, resolveLabelFitChoice, editLabelFitRecordText } = createLabelFitDialog({
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

  const circularRecordSelectionEnabled = computed(() => semanticMutationAvailable.value
    && mode.value === 'circular'
    && !state.activeDrawing().form.multi_record_canvas && circularRecordDiscoveryState.value.status === 'ready');
  // CI-03: the per-record presentation a source file carries in each mode.
  // Replacing or removing the file resets it in the replacement's History step.
  // The OFF records name records of the old file too (record selection).
  /** @param {DrawingState} drawing */
  const retireCircularRecordSelector = ({ form, adv, recordsOff }) => {
    form.circular_record_selector = CIRCULAR_RECORD_FORM_DEFAULTS.circular_record_selector;
    if (adv.circular_grouping_intent === 'single') adv.circular_grouping_intent = 'auto';
    if (recordsOff.length) recordsOff.splice(0);
  };
  // The whole record presentation: the record form fields return to their
  // defaults (the one list a style copy leaves out too). Each field is written
  // by name, so the track-layout writer guard reads no whole-form write; a
  // node case pins that every field of the list is written.
  /** @param {DrawingState} drawing */
  const retireCircularRecordPresentation = (drawing) => {
    const { form } = drawing;
    retireCircularRecordSelector(drawing);
    form.circular_region_start = CIRCULAR_RECORD_FORM_DEFAULTS.circular_region_start;
    form.circular_region_end = CIRCULAR_RECORD_FORM_DEFAULTS.circular_region_end;
    form.circular_reverse = CIRCULAR_RECORD_FORM_DEFAULTS.circular_reverse;
    form.circular_record_label = CIRCULAR_RECORD_FORM_DEFAULTS.circular_record_label;
    form.circular_record_subtitle = CIRCULAR_RECORD_FORM_DEFAULTS.circular_record_subtitle;
  };
  // D-32: one record replacing one record (a new version of the same genome)
  // keeps its crop and titles; a removed file or a replaced multi-record file
  // retires them. The record selector names a record of the old file, so any
  // replacement retires it.
  /** @param {{ removed: boolean, previousRecordCount: number }} replacement */
  const sourceReplacementRetiresPresentation = ({ removed, previousRecordCount }) => removed || previousRecordCount > 1;
  /**
   * A Circular drawing's source file (upload, replace, Remove), with the
   * drawing's reactions to its new source: the shown drawing's controls and
   * the replacement for every drawing that uses a file write through it.
   * @param {DrawingState} drawing
   * @param {'c_gb' | 'c_gff' | 'c_fasta'} field
   * @param {InputFileCandidate | InputFileRef} value
   */
  const replaceCircularSource = (drawing, field, value) => {
    const nextValue = value ?? null;
    const previous = drawing.files[field];
    const nextRef = state.inputFiles.registerFile(nextValue);
    if (previous === nextRef) return;
    const previousRecordCount = drawing.circularRecordList.value.length;
    changeRecordSources(drawing, () => { bindDrawingInputFile(state.inputFiles, drawing.slices.inputs, field, nextValue); });
    if (nextRef) pendingDiscoveryOf(drawing).circularUpload = nextRef;
    void discoverDrawingRecords(drawing);
    // An upload into an empty slot starts reading another input (after an
    // Input switch), whose records the OFF list does not name.
    if (!previous) {
      if (drawing.recordsOff.length) drawing.recordsOff.splice(0);
      return;
    }
    if (sourceReplacementRetiresPresentation({
      removed: !nextValue, previousRecordCount
    })) retireCircularRecordPresentation(drawing);
    else retireCircularRecordSelector(drawing);
  };
  /**
   * The Circular source file controls of the shown drawing.
   * @param {'c_gb' | 'c_gff' | 'c_fasta'} field
   * @param {File | null} value
   */
  const setCircularSourceFile = (field, value) => shownAction('circular', (drawing) => replaceCircularSource(drawing, field, value));
  /** @param {unknown} value */
  const setCircularRecordPresentationSelector = (value) => {
    const busy = sessionOperationAvailability();
    if (busy) return busy;
    const drawing = shownDrawingOfMode('circular');
    if (!drawing) return { status: 'unavailable' };
    if (!circularRecordSelectionEnabled.value) return { status: 'unavailable' };
    const normalized = String(value || '').trim();
    drawing.form.circular_record_selector = normalized;
    drawing.adv.circular_grouping_intent = normalized ? 'single' : 'auto';
    rebindAnnotationRecords();
    return { status: 'ok' };
  };
  const showCircularCanvasSetting = async () => {
    const drawing = shownDrawingOfMode('circular');
    const control = /** @type {HTMLElement | null} */ (document.querySelector('[data-circular-canvas-setting]'));
    if (!control || !drawing || !drawing.form.multi_record_canvas) {
      return { status: 'unavailable' };
    }
    await settingsPaneActions.revealSetting(control);
    return { status: 'ok' };
  };
  // With Custom Track Slots on, a Ticks row owns the coordinate scale: the link beside Show Coordinate Scale
  // reveals that row (the first enabled one, else the first), or the Custom Track Slots panel when the stack
  // has none. It changes no setting.
  const showCircularTicksSlot = async () => {
    const drawing = state.activeDrawing();
    if (mode.value !== 'circular' || !drawing.adv.circular_track_slots_enabled) return { status: 'unavailable' };
    circularTrackSlotsPanelOpen.value = true;
    await nextTick();
    /** @type {CircularTrackSlot[]} */
    const slots = Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots : [];
    const ticks = slots.filter((slot) => slot.renderer === 'ticks');
    const slot = ticks.find((candidate) => candidate.enabled) || ticks[0];
    const row = slot ? document.querySelector(`[data-capture="circular-track-slot-${CSS.escape(slot.id)}"]`) : null;
    const panel = document.getElementById('circular-custom-track-slots-panel');
    const toggle = /** @type {HTMLElement | null} */ (document.querySelector('[data-track-slots-toggle]'));
    const control = /** @type {HTMLElement | null} */ (row?.querySelector('input, select, button') ?? null);
    if (row instanceof HTMLElement && control) {
      await settingsPaneActions.revealSetting(row, { focus: control });
    } else if (panel && toggle) {
      await settingsPaneActions.revealSetting(panel, { align: 'start', focus: toggle });
    } else {
      return { status: 'unavailable' };
    }
    return { status: 'ok' };
  };
  watch(() => [mode.value, circularSingleRecordPresentationEnabled.value,
    circularRecordDiscoveryState.value.primaryFile, circularRecordDiscoveryState.value.pairedFile,
    state.activeDrawing().form.circular_record_selector], async (/** @type {unknown[]} */ current, /** @type {unknown[]} */ previous = []) => {
    if (!current[1] || current.every((value, index) => Object.is(value, previous[index]))) return;
    const origin = document.activeElement;
    const pane = document.querySelector('.settings-scroll');
    const scrollOwner = pane ? settingsScrollOwner(pane) : document.scrollingElement;
    const anchored = origin && pane?.contains(origin) && origin.getClientRects().length;
    const top = anchored ? origin.getBoundingClientRect().top : null;
    const scrollTop = scrollOwner?.scrollTop;
    await nextTick();
    if (!circularSingleRecordPresentationEnabled.value
      || circularRecordDiscoveryState.value.primaryFile !== current[2]
      || circularRecordDiscoveryState.value.pairedFile !== current[3]) return;
    if (circularRecordPresentationPanel.value) circularRecordPresentationPanel.value.open = true;
    await afterPaint();
    if (scrollOwner && document.activeElement === origin) {
      // `top` is a number whenever `anchored` is truthy; `scrollTop` was read from
      // the same non-null `scrollOwner`.
      scrollOwner.scrollTop = anchored && origin.isConnected
        ? scrollOwner.scrollTop + origin.getBoundingClientRect().top - /** @type {number} */ (top)
        : /** @type {number} */ (scrollTop);
    }
  }, { flush: 'sync' });

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
  const {
    cancelLinearSourceRemoval, applyLinearSourceRemoval, requestLinearSourceRemoval, removeLastLinearSeq
  } = createLinearSourceRemovalDialog({
    linearSourceGroups, linearSourceRemovalDialog, linearSourceRemovalReturnFocus, sessionOperationAvailability,
    shownDrawingOfMode, runUndoable: history.runUndoable, removeLinearSource, askLinearFileRemoval
  });
  // Keeps Tab inside the modal dialog of the overlay that handles the keydown.
  /** @param {KeyboardEvent & { currentTarget: Element }} event Tab on the overlay that holds the dialog. */
  const trapDialogFocus = (event) => {
    const dialog = /** @type {HTMLElement | null} */ (event.currentTarget?.querySelector('[role="dialog"]'));
    const controls = /** @type {HTMLElement[]} */ (Array.from(dialog?.querySelectorAll('button:not(:disabled)') || []));
    if (!controls.length) {
      // A busy choice disables every button: Tab stays on the dialog panel.
      if (dialog) {
        event.preventDefault();
        dialog.focus({ preventScroll: true });
      }
      return;
    }
    // `controls` is not empty, so `dialog` is an element and `last` is a button.
    const first = controls[0];
    const last = /** @type {HTMLElement} */ (controls.at(-1));
    if (event.shiftKey && (document.activeElement === first || !(/** @type {HTMLElement} */ (dialog)).contains(document.activeElement))) {
      event.preventDefault();
      last.focus();
    } else if (!event.shiftKey && (document.activeElement === last || !(/** @type {HTMLElement} */ (dialog)).contains(document.activeElement))) {
      event.preventDefault();
      first.focus();
    }
  };
  const { replaceLinearSourceFile, setLinearSeqPrimaryFile } = createLinearSourceFileReplacement({
    state, shownAction, pendingDiscoveryOf, applyLinearSeqMutation, askLinearFileRemoval, sourceReplacementRetiresPresentation
  });
  /** @param {readonly (InputFileRef | null)[]} slots @param {string} fileId @returns {number[]} */
  const slotIndexesOfFile = (slots, fileId) => slots.flatMap((ref, idx) => (ref && ref.fileId === fileId ? [idx] : []));
  /**
   * The writer of each slot an input file replacement names: every
   * `ReplaceableInputSlot` has one (a missing one fails the type check), and
   * each is the body the shown drawing's control of that slot runs.
   * @type {Readonly<Record<ReplaceableInputSlot, (drawing: DrawingState, site: DrawingInputSite, fileId: string, file: InputFileCandidate | InputFileRef) => void>>}
   */
  const INPUT_REPLACEMENT_WRITERS = Object.freeze({
    c_gb: (drawing, _site, _fileId, file) => replaceCircularSource(drawing, 'c_gb', file),
    c_gff: (drawing, _site, _fileId, file) => replaceCircularSource(drawing, 'c_gff', file),
    c_fasta: (drawing, _site, _fileId, file) => replaceCircularSource(drawing, 'c_fasta', file),
    c_depth: (drawing, _site, fileId, file) => {
      const indexes = new Set(circularDepthRows(drawing).flatMap((row) => slotIndexesOfFile(row, fileId)));
      indexes.forEach((idx) => replaceCircularDepthFile(drawing, idx, file));
    },
    c_conservation_blasts: (drawing, _site, fileId, file) => {
      const ref = state.inputFiles.registerFile(file);
      replaceCircularConservationBlasts(drawing, normalizeFileList(drawing.files.c_conservation_blasts)
        .map((entry) => (entry.fileId === fileId ? ref : entry)));
    },
    'row.gb': (drawing, site, _fileId, file) => replaceLinearSourceFile(drawing, String(site.rowUid), 'gb', file),
    'row.gff': (drawing, site, _fileId, file) => replaceLinearSourceFile(drawing, String(site.rowUid), 'gff', file),
    'row.fasta': (drawing, site, _fileId, file) => replaceLinearSourceFile(drawing, String(site.rowUid), 'fasta', file),
    'row.depth': (drawing, site, fileId, file) => {
      const seq = drawing.linearSeqs.find((entry) => entry.uid === site.rowUid);
      if (!seq) return;
      slotIndexesOfFile(depthFileSlotsFromValue(seq.depth), fileId).forEach((idx) => setLinearDepthFiles(drawing, [seq], idx, file));
    },
    edge: (drawing, site, _fileId, file) => swapLinearComparisonFile(drawing, String(site.edgeId), file)
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
    applyInputReplacementToDrawings: (sites, fileId, file) => sites.forEach((site) => {
      INPUT_REPLACEMENT_WRITERS[site.slot](state.drawingOf(site.drawingId), site, fileId, file);
    }),
    askScope: async (previous, otherDrawingIds) => {
      // Named as their tabs show them: a default name in the UI language.
      const others = drawingBar.tabs().flatMap((tab) => (
        tab.kind === 'drawing' && otherDrawingIds.includes(tab.id) ? [tab.label] : []));
      const choice = await askChoice({
        title: t('Replace input file?'),
        description: tn(others.length, '{name} is also used in {drawings}.', '{name} is also used in {count} drawings: {drawings}.', {
          name: previous.name, drawings: listText(others)
        }),
        actions: [
          { id: 'drawing', label: t('This drawing only'), tone: 'secondary' },
          { id: 'all', label: t('Every drawing that uses this file'), tone: 'primary' }
        ]
      });
      return choice === 'drawing' || choice === 'all' ? choice : null;
    },
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
    shownAction, replaceFromShownDrawing, replaceCircularSource, replaceCircularDepthFile, replaceLinearSourceFile,
    setLinearDepthFiles, setLinearSourceDepthFile, linearComparisonPairForEdgeKey, setLinearComparisonCardFile
  });

  // Record selection (D-01..D-10): the record sources and edits its owner
  // reads, and the transitions it asks the root for.
  /** @param {DrawingState} drawing the drawing that holds the card @param {LinearSeq} seq The discovered record a Linear card draws. */
  const linearCardRecord = (drawing, seq) => {
    const selection = resolveDisambiguatedRecordSelection(linearRecordSelector.recordsFor(drawing, seq), seq.region_record_id);
    return selection.record || (selection.entries.length === 1 ? selection.entries[0] : null);
  };
  /** @param {DrawingState} drawing the drawing that holds the card @param {LinearSeq} seq @returns {RecordSettingsOwner} */
  const linearRecordSettingsOwner = (drawing, seq) => ({
    key: seq.uid, requestKeys: [seq.uid], ownsExpansions: true,
    bindingKeys: annotationRecordKeysOf(linearAnnotationSource(drawing, seq, false), seq.uid),
    displaySource: seq.uid, displaySelector: null
  });
  /** @param {DrawingState} drawing */
  const circularRecordEntries = (drawing) => buildDisambiguatedRecordEntries(
    getCircularRecordDiscoveryState(drawing).records.map((record) => ({
      ...record, recordId: record?.record_id ?? record?.recordId, recordLength: record?.record_length ?? record?.recordLength
    }))
  );
  /** @returns {RecordSettingsOwner[]} */
  /** @param {DrawingState} drawing */
  const circularRecordSettingsOwners = (drawing) => {
    const discovery = getCircularRecordDiscoveryState(drawing);
    const bindingKeys = annotationRecordKeysOf(
      { sourceKey: circularAnnotationSourceKey(drawing), hasInput: discovery.hasInput, status: discovery.status, records: discovery.records },
      'circular-source'
    );
    return circularRecordEntries(drawing).map((entry, index) => ({
      key: circularRecordDrawKey(entry),
      requestKeys: circularRecordRequestKeys(entry),
      bindingKeys: bindingKeys[index] ? [bindingKeys[index]] : [],
      displaySource: 'circular', displaySelector: circularRecordDrawKey(entry)
    }));
  };
  /** The record settings owners of a drawing's records. @param {DrawingState} drawing */
  const recordSettingsOwners = (drawing) => (drawing.mode === 'circular' ? circularRecordSettingsOwners(drawing) : drawing.linearSeqs.map((seq) => linearRecordSettingsOwner(drawing, seq)));
  /** @param {DrawingState} drawing @returns {RecordEditData} */
  const recordEditData = (drawing) => ({
    featureEdits: drawing.slices.featureEdits,
    paint: drawing.slices.paint,
    annotationSets: drawing.annotationSets,
    removeAnnotation: (/** @type {object} */ annotation) => removeAnnotations(drawing.slices.annotationSets, [annotation]),
    colorRuleEdits: drawing.manualSpecificRules.flatMap((row) => {
      const instance = namedFeatureInstance(row);
      return instance ? [{
        recordKey: instance.recordKey,
        remove: () => featureActions.takeBackRuleRows(drawing, (rules) => rules.filter((rule) => rule !== row))
      }] : [];
    }),
    recordDisplayDrafts: drawing.recordDisplayDrafts,
    comparisonEdges: drawing.mode === 'linear' ? drawing.linearComparisonPlan.edges : [],
    annotationBindingField: ANNOTATION_RECORD_BINDING_KEY
  });
  /** @param {RecordSelectionMode} listMode @param {string} sourceKey @returns {RecordSource | null} */
  const recordSourceOf = (listMode, sourceKey) => {
    const drawing = shownDrawingOfMode(listMode);
    if (!drawing) return null;
    if (listMode === 'circular') {
      if (sourceKey !== 'circular' || getCircularRecordDiscoveryState(drawing).status !== 'ready') return null;
      const owners = circularRecordSettingsOwners(drawing);
      if (owners.length === 0) return null;
      const edited = recordsWithEdits(owners, recordEditData(drawing));
      const entries = circularRecordEntries(drawing);
      return {
        key: 'circular',
        name: (({ primaryFile, pairedFile }) => [primaryFile, pairedFile])(inputSourceFiles(drawing.inputType.value, circularInputSlots(drawing.files)))
          .filter(Boolean).map((file) => String(file?.name || 'Unnamed file')).join(' + '),
        records: owners.map((owner, index) => ({
          key: owner.key, recordId: entries[index].recordId, length: entries[index].recordLength ?? null,
          hasSettings: edited.has(owner.key)
        }))
      };
    }
    const group = linearSourceGroups.value.find((entry) => entry.uid === sourceKey);
    if (!group) return null;
    const cards = group.records.map(({ sequence }) => sequence);
    const edited = recordsWithEdits(cards.map((seq) => linearRecordSettingsOwner(drawing, seq)), recordEditData(drawing));
    return {
      key: group.uid,
      name: linearSourceFileName(group),
      records: cards.map((seq) => {
        const record = linearCardRecord(drawing, seq);
        return {
          key: seq.uid, recordId: String(record?.recordId || seq.region_record_id || ''),
          length: record?.recordLength ?? null, hasSettings: edited.has(seq.uid) || linearCardHasSettings(seq)
        };
      })
    };
  };
  /**
   * Delete settings (D-08), inside the owner's History step: the records' card
   * settings and rotation, feature edits and placements, annotations, and
   * comparison pairs.
   * @param {RecordSelectionMode} listMode
   * @param {RecordDrawKey[]} keys
   */
  const deleteRecordSettings = (listMode, keys) => {
    const drawing = shownDrawingOfMode(listMode);
    if (!drawing) return false;
    const deleted = new Set(keys);
    const owners = recordSettingsOwners(drawing).filter((owner) => deleted.has(owner.key));
    const ruleCount = drawing.manualSpecificRules.length;
    retireLegendStylesOfUnnamedCaptions(drawing, () => removeRecordEdits(owners, recordEditData(drawing)));
    // The color rows it took back color nothing drawn; the rule matches follow the table.
    if (drawing.manualSpecificRules.length !== ruleCount) void projectPaletteAndRules();
    if (listMode === 'linear') {
      applyLinearSeqMutation(drawing, drawing.linearSeqs.map((seq) => (
        deleted.has(seq.uid) ? { ...seq, ...RECORD_SETTINGS_DEFAULTS } : seq
      )), { alignmentMutation: msg('record settings deleted.') });
    }
    return true;
  };
  /** @param {DrawingState} drawing @param {string} sourceKey The Linear File's key in the record list. */
  const recordFileRemovalPlan = (drawing, sourceKey) => {
    const intent = /** @type {'clear' | 'delete'} */ (groupLinearSourceRecords(drawing.linearSeqs).length > 1 ? 'delete' : 'clear');
    return { intent, plan: planLinearSourceRemoval({ sequences: drawing.linearSeqs, sourceUid: sourceKey, intent }) };
  };
  /**
   * D-06 Remove File, inside the dialog's History step: a Linear File is
   * deleted, or cleared when it is the only File; the Circular input loses its
   * GenBank file or its GFF3 and FASTA pair.
   * @param {RecordSelectionMode} listMode
   * @param {string} sourceKey
   */
  const removeRecordSource = (listMode, sourceKey) => {
    if (listMode === 'circular') {
      // The shown Circular drawing's input, which its File card removes from.
      const drawing = shownDrawingOfMode('circular');
      if (drawing === null) return undefined;
      const { primary, paired } = inputTypeSpec(drawing.inputType.value);
      if (paired) setCircularSourceFile(CIRCULAR_INPUT_FIELDS[primary], null);
      return setCircularSourceFile(CIRCULAR_INPUT_FIELDS[paired ?? primary], null);
    }
    const drawing = shownDrawingOfMode('linear');
    if (!drawing) return false;
    const { intent, plan } = recordFileRemovalPlan(drawing, sourceKey);
    return plan.allowed ? removeLinearSource(drawing, plan, intent) : false;
  };
  /** @param {RecordSelectionMode} listMode @param {string} sourceKey @param {string} label */
  const askRecordFileRemoval = (listMode, sourceKey, label) => {
    const drawing = shownDrawingOfMode('linear');
    if (listMode !== 'linear' || !drawing || dialogChoice.pending.value || history.mutationPending()) return null;
    const { intent, plan } = recordFileRemovalPlan(drawing, sourceKey);
    return plan.allowed
      ? askLinearFileRemoval(drawing, plan.retainedSequences, label, (turnOffReaching) => removeLinearSource(drawing, plan, intent, turnOffReaching))
      : null;
  };
  // After a Session load the records of the shown drawing's inputs are read on demand, as Generate reads them: the
  // record list and Select region's Add to annotation set ask for them.
  /** @param {RecordSelectionMode} listMode */
  const readShownRecords = async (listMode) => {
    if (listMode === 'linear') return linearRecordSelector.refresh(state.activeDrawing());
    const drawing = shownDrawingOfMode('circular');
    if (!drawing) return null;
    return getCircularRecordDiscoveryState(drawing).status === 'ready' ? null : inspectCircularSourceRecords();
  };
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
    setLinearRecordSelector, setLinearRecordCrop, resetLinearRecordDefinition, resetLinearRecordSubtitle
  } = createLinearRecordEdits({
    sessionOperationAvailability, shownDrawingOfMode, runUndoable: history.runUndoable,
    clearAlignmentForMutation: similarityAlignmentActions.clearForMutation,
    linearRecordsFor: linearRecordSelector.recordsFor, linearRecordStatusFor: linearRecordSelector.statusFor,
    rebindAnnotationRecords
  });
  const { setCircularInputType, setLinearInputType } = createInputTypeEdits({
    sessionOperationAvailability, shownDrawingOfMode, changeRecordSources, discoverDrawingRecords,
    clearAlignmentForMutation: similarityAlignmentActions.clearForMutation
  });

  const labelSettingsShown = computed(() => {
    const drawing = state.activeDrawing();
    const scope = mode.value === 'linear' ? drawing.form.show_labels_linear : drawing.form.labels_mode;
    return labelSettingsVisible(scope, drawing.featureOverrides);
  });

  // The template's names for the members of the shown mode's drawing: a ref
  // member is its value, any other member the object itself. Only a setting's
  // ref takes a write, its control's (plan-r2 §13 item 3); an editor member's
  // writers are ports, and a derived value is read-only.
  /** @param {keyof DrawingState} key */
  const drawingMember = (key) => {
    if (!isRef(state.activeDrawing()[key])) return computed(() => state.activeDrawing()[key]);
    if (!(/** @type {readonly string[]} */ (SETTINGS_SLICE.members)).includes(key)) return computed(() => state.activeDrawing()[key].value);
    return computed({
      get: () => state.activeDrawing()[key].value,
      /** @param {unknown} value */
      set: (value) => SETTINGS_SLICE.install(state.activeDrawing().slices.settings, { [key]: value })
    });
  };

  // The default colors of the feature types the shown figure draws,
  // and the user's own default colors, come first; the others wait under More.
  // The values are read untracked, so editing a color never moves it to the
  // other grid; the grouping follows the drawn types, the palette, the keys,
  // and a replaced color table (Reset, palette change, Session load, an
  // imported Default colors file).
  const defaultColorGroups = computed(() => {
    const drawing = state.activeDrawing();
    /** @type {Record<string, string>} */
    const colors = drawing.currentColors.value || {};
    const keys = Object.keys(colors);
    const drawn = new Set((state.extractedFeatures.value || [])
      .map((/** @type {{ type?: string }} */ feature) => String(feature?.type || '')));
    if (!drawn.size) return { primary: keys, more: /** @type {string[]} */ ([]) };
    const circular = mode.value === 'circular';
    if (circular ? !drawing.form.suppress_gc : drawing.form.show_gc) drawn.add('gc_content');
    if (circular ? !drawing.form.suppress_skew : drawing.form.show_skew) ['skew_high', 'skew_low'].forEach((key) => drawn.add(key));
    if (!circular && drawing.linearSeqs.length > 1 && linearComparisonGlobalAction.value !== 'none') {
      COMPARISON_COLOR_KEYS.forEach((key) => drawn.add(key));
    }
    resultsManager.userDefaultColorRows(drawing, toRaw(colors)).forEach(([key]) => drawn.add(key));
    const primary = keys.filter((key) => key === 'default' || drawn.has(key));
    return { primary, more: keys.filter((key) => !drawn.has(key) && key !== 'default') };
  });

  // New drawing from a region: the shown drawing's draft request, as Save stores it, goes to Python.
  const projectShownRequest = async () => {
    const drawing = state.activeDrawing();
    const { catalog, error } = drawing.mode === 'linear'
      ? await prepareLinearRecordCatalog({ privateCandidate: true })
      : { catalog: null, error: '' };
    if (error) throw error;
    return projectDraftRequest(state, drawing, {
      linearRecordCatalog: catalog,
      recordDisplayRows: recordDisplayControls.allRows.value || [],
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
      const item = clickedPairwiseMatch.value?.previewItem;
      return item?.kind === 'region' ? item : null;
    },
    openRegionPopup: (selection, element, anchor) => featureActions.openPreviewItemPopup({ element, item: selection }, anchor),
    closeRegionPopup,
    drawnFeatures: () => toRaw(extractedFeatures.value) ?? [],
    selectFeatures: (featureIds) => featureSelection.replaceSelectedFeatureIds(new Set(featureIds), { anchorId: featureIds[0] ?? '' }),
    // Binds as a search hit does; the popup is teleported out of the app root (no History adapter): its own step.
    addRegionAnnotation: async ({ recordKey, start, end }) => {
      const drawing = state.activeDrawing();
      await readRecordCatalog(drawing);
      if (drawing !== state.activeDrawing()) return null;
      return history.runUndoable(msg('Add to annotation set'), () => annotationEditor.addRegionAnnotation({ record: resultRecordBinding(shownResultRecordBindings(), recordKey), start, end }));
    },
    newDrawingFromRegion: regionDrawings.openWithRegion,
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
    showLinearComparisonCause,
    runLinearComparisonIssueAction,
    resolveMissingBlastTsv,
    labelFitDialog,
    labelCapacityText,
    labelFitNotice,
    shownInputSummary,
    shownInputSummaryStale,
    openLabelFitNotice,
    resolveLabelFitChoice,
    editLabelFitRecordText,
    sessionSaveAvailable,
    sessionLoadAvailable,
    sessionBusyReason,
    generationCancelRequested,
    errorLog,
    errorDisplay,
    canRetryInteractiveSvgExport,
    reloadAfterOperationError,
    selectDrawing,
    drawingList,
    regionDrawings,
    regionSelect,
    regionSelectOn: regionSelect.active,
    regionSelectAvailable: regionSelect.available,
    regionSelectUnavailableReason: regionSelect.unavailableReason,
    regionOverlay: regionSelect.overlay,
    regionOverlayMatrix: regionSelect.overlayMatrix,
    regionPopup: regionSelect.popup,
    drawingBar,
    drawingBarTabs,
    moveTabFocus,
    drawingSwitchAvailable,
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
    circularConservationSeriesRows,
    setCircularConservationUploadFiles,
    setCircularConservationCompanionFile,
    canMoveCircularConservationSeries,
    moveCircularConservationSeries,
    openCircularConservationComparisonFilePicker,
    addCircularConservationComparisonFile,
    removeCircularConservationSource,
    changeCircularConservationSource,
    refreshCircularConservationSlots: refreshShownCircularConservationSlots,
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
    linearRecordGencodeRequired,
    linearSourceGencodeRequired,
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
    linearLosatCacheInfoByEdgeKey,
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
    setLinearRecordCrop,
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
    addCircularTrackSlot: circularTrackSlotEditor.addCircularTrackSlot,
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
    addLinearTrackSlot: linearTrackSlotEditor.addLinearTrackSlot,
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
    orthogroupCopyFeedbackLabel: orthogroupActions.orthogroupCopyFeedbackLabel,
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
    DEFAULT_LINEAR_OUTPUT_PREFIX,
    setCircularRecordPresentationSelector,
    changeCircularMultiRecordCanvas,
    setCircularSourceFile,
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
    downloadLosatCache,
    downloadLosatPair,
    setLosatPairFilename,
    clearLosatCache,
    getLosatPairDefaultName,
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
    hidePreviewFeatureSearchBar: previewFeatureSearch.hideSearchBar,
    togglePreviewFeatureSearchBar: previewFeatureSearch.toggleSearchBar,
    selectedFeatureIds,
    selectedFeatureAnchorId,
    featureSelectionStatus,
    featureSelectionDrag,
    selectedFeatureCount,
    selectedFeatures,
    hasFeatureSelection,
    featureSelectionMarqueeStyle: featureSelection.featureSelectionMarqueeStyle,
    featureSelectionToolbarStyle: featureSelection.featureSelectionToolbarStyle,
    startFeatureSelectionToolbarDrag: featureSelection.startToolbarDrag,
    clearFeatureSelection: featureSelection.clearFeatureSelection,
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
    clickedLabel,
    clickedLabelPos,
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
    downloadText,
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
    // A Legend row delete is one History step of the intent, which Undo and
    // Redo show through the port and lay out as Python does, canvas included
    // (OV-125). The deleted row's stroke leaves its features (OV-293).
    shownLegendName: legendActions.shownLegendName,
    shownLegendEntryName: legendActions.shownLegendEntryName,
    legendEntryMerged: legendActions.legendEntryMerged,
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
    linearRegionNote,
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
