// @ts-check
/** @import { SvgMutationOperations } from '../services/svg-result-ingestion.js' */
/** @import { FileLike } from '../services/session-resource-backing.js' */
/** @import { ResultLabelFit } from '../services/label-fit.js' */
/** @import { DiagramMode, DrawingState, LegendEntry, LinearSeq, PopupPosition } from '../state.js' */
/** @import { RulePreparation } from './rule-matching.js' */
/** @import { CircularDiagramRecords } from '../services/record-options.js' */
/** @import { EditorPaintState } from './result-paint-record.js' */
/** @import { FeatureEditorOptions } from './feature-editor.js' */
/** @import { ColorActionsRuleActions } from './feature-editor/color-actions.js' */
/** @import { LegendRowColorOf } from './candidate-render.js' */
/** @import { UserFacingError } from '../utils/error-normalization.js' */
/** @import { AnnotationCatalogSource } from './annotations/record-catalog.js' */
/** @import { LegacyResultSvgData, LegacyResultSvgTransform, SessionFileInput } from '../services/config.js' */
/** @import { RecordKeyedFeature } from '../services/feature-identity.js' */
/** @import { ArtifactSlot } from '../services/artifact-slot.js' */
/** @import { DisplayedLegendRowContext } from '../services/feature-visibility.js' */
/** @import { GalleryExample } from '../services/gallery-examples.js' */
/** @import { LinearComparisonEdge, LinearComparisonIssue, LinearComparisonPair, LinearComparisonPlan, LinearComparisonResolvedEdge, LinearComparisonSource } from '../services/linear-comparisons.js' */
/** @import { LinearComparisonUiIssue } from './comparison-ui.js' */
/** @import { LinearSourceGroup, LinearSourceRemovalPlan } from '../services/linear-sources.js' */
/** @import { RecordDrawKey, RecordSettingsOwner, RecordEditData } from '../services/record-draw-selection.js' */
/** @import { RecordListRequest, RecordSelectionMode, RecordSource } from './record-selection.js' */
/** @import { AlignmentFigureGeometry } from '../services/alignment-review-view.js' */
/** @typedef {{ opening: Readonly<ArtifactSlot> | null, stashed: Readonly<ArtifactSlot> | null }} LoadedArtifactSlots */
import { createRulePreparation } from './rule-matching.js';
import { compileDirectEditorMutationPlan, displayedLegendRowColors, editorPaintDomains, LIVE_EDIT_DOMAINS } from './candidate-render.js';
import { createResultPaintRecord, editorPaintState, paintRuleTable, sameColors } from './result-paint-record.js';
import {
  countUnresolvedFeatureEdits, displayedLegendRowContext, featureDrawnContext, removeUnresolvedFeatureEdits,
  requestFeatureVisibilityRules
} from '../services/feature-visibility.js';
import { drawnBlockStroke, drawsPythonLegendRow, legendEntryKey, legendStructureEdited } from '../services/legend-svg.js';
import { admitFeatureCatalog } from '../services/feature-catalog.js';
import { readArtifactOfMode } from '../services/artifact-slot.js';
import { labelSettingsVisible } from '../services/feature-placement.js';
import { displayedFeatureAddressing } from '../services/feature-override-identity.js';
import { PNG_DPI_CHOICES, createDefaultLosatpHitLimits } from '../services/session-active-config-contract.js';
import { settingsScrollOwner } from '../services/settings-sections.js';
import { createSettingsPaneController } from './settings-pane.js';
import { createRecordDisplayControls, RECORD_TARGET_NOT_DISCOVERED } from './record-display-options.js';
import { createPreviewItemPopup, defaultColorScopeLine } from './preview-item-popup.js';
import { colorRuleLegendGroups, groupLegendRowKeys, legendRowRules } from '../services/specific-color-rules.js';
import { drawnTracks, hiddenTrackLegendKeys, loadedTracks, trackLegendRows } from '../services/track-legend-rows.js';
import { isCurrentFeature } from '../services/feature-identity.js';
import {
  createFeatureRecordRotationAction,
  createFeatureRecordRotationWorkflow
} from './record-display/feature-record-rotation.js';
import { state, sessionOperationAvailability, createLinearSeq, normalizeLinearSeqList } from '../state.js';
import {
  adoptCanonicalRenderArtifacts,
  applyCapturedConfigData,
  applyEditorStateData,
  applyFeatureStateData,
  applyOrthogroupStateData,
  applyResultsData,
  applyRunStateData,
  applyUiStateData,
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
  committedLiveSettingApplies,
  commitLiveSettingValues,
  readCommittedCanonicalRevision,
  getLoadedTrackSettings,
  getCommittedLinearDefinitionVisibility,
  readCommittedResourceRecordCount,
  assertActiveModeInputs,
  importSession as importSessionFromFile,
  SESSION_VERSION,
  serializeActiveRenderFiles,
  serializeResults,
  setUnmanagedConfigOverrideValidator,
  setPreviewRuntime,
  applyLabelFitValues,
  setHeldPreviewViewportReader
} from '../services/config.js';
import {
  ALIGNMENT_REFERENCE_KEY, alignmentOverlayRecords, buildAlignmentReviewView
} from '../services/alignment-review-view.js';
import {
  setMainSessionComparisonFrameConverter, setStoredProteinResultFrameConverter
} from '../services/main-session-comparison-frame.js';
import { fetchGalleryExampleSession, readGalleryExamples } from '../services/gallery-examples.js';
import { createHistoryManager } from '../services/history.js';
import { createHistoryFileStore } from '../services/history-files.js';
import { createHistorySnapshotService } from '../services/history-snapshot.js';
import { cloneJsonData } from '../services/json-clone.js';
import {
  groupLinearSourceRecords,
  linearCardFileFields,
  linearCardInputComplete,
  linearRecordNames,
  linearFileCardRecord,
  isPristineLinearSource,
  linearSourceHasPrimaryInput,
  linearSourceDepthStatus,
  moveLinearSourceGroup,
  planLinearSourceRemoval,
  setLinearSourceFolded,
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
import { downloadTextFile } from '../services/text-download.js';
import { resetLayoutState, resetSettings as resetSettingsState } from '../services/reset.js';
import { allLiveSettingRequestPaths, liveSettingRequestPaths, SETTING_NOTE_APPLIES } from '../services/setting-application.js';
import {
  DIAGRAM_HELPER_OPERATIONS,
  disposeDiagramGenerationWorker,
  runDiagramHelperOperation,
  setDiagramWorkerFontSource
} from '../services/diagram-generation.js';
import { createFontSettings } from './font-settings.js';
import { customFontSaveWarnings } from '../services/session-fonts.js';
import { requestFontFamily, webFontFamilyNames } from '../services/font-store.js';
import { SLICE_MODES } from '../services/mode-scoped-migration.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from '../services/runtime-test-hooks.js';
import { createPanZoom, createSidebarResize, setupGlobalUiEvents } from './ui.js';
import { appliedFeatureColors, COMPARISON_COLOR_KEYS, colorValueMode, normalizeOptionalHexColor, toNativeColorInputValue } from '../utils/color-utils.js';
import { createFeatureEditor } from './feature-editor.js';
import { groupCountText } from './pairwise-match-popup.js';
import { createFeatureSelection } from './feature-selection.js';
import { createPreviewFeatureSearch } from './feature-search/preview-actions.js';
import { createSvgStyles, paletteColorsEqual } from './svg-styles.js';
import { createLegendManager } from './legend.js';
import { createPaletteLoader } from './palettes.js';
import { afterPaint, buildDefaultMultiRecordPositions, createRunAnalysis } from './run-analysis.js';
import { createSimilarityAlignmentActions } from './similarity-alignment.js';
import { createLabelFitDialog } from './label-fit-dialog.js';
import { diagnosticError, normalizeCaughtError, normalizeUserFacingError } from '../utils/error-normalization.js';
import { placePopupNearPoint, placePopupWithin, popupViewportBounds } from '../utils/popup-bounds.js';
/** @import { PopupPlacement } from '../utils/popup-bounds.js' */
import { formatElapsedMs, reproducibilityLabel } from './run-info.js';
import { createLegendLayout } from './legend-layout.js';
import {
  COMPOSITION_METADATA_ATTRIBUTE,
  COMPOSITION_SCHEMA_ATTRIBUTE,
  normalizeLegacyComposition,
  normalizeStoredLegendSide,
  storedLegendSideIsOfMode
} from './legend-layout/composition-actions.js';
import { createResultsManager } from './results.js';
import { createAuxiliaryFileImports, setupWatchers } from './watchers.js';
import { createDialogChoice, setupHistoryInputs } from './history-inputs.js';
import { createChoicePrompt, createNoticeBanner } from './app-messages.js';
import { setupHistoryShortcuts } from './history-shortcuts.js';
import { createPreviewRuntime } from './preview-runtime.js';
import {
  createOrthogroupEditor,
  resolveUniqueOrthogroupMemberForFeature
} from './orthogroups.js';
import { createRightDrawerController } from './right-drawer.js';
import {
  createCircularTrackSlotEditor,
  estimateCircularConservationLayoutWarning
} from './circular-track-slots.js';
import { createLinearTrackSlotEditor } from './linear-track-slots.js';
import { createAnnotationEditor } from './annotations.js';
import { buildLegendStyleRetirement, trackDataLegendCaptions } from './legend/track-data-styles.js';
import {
  annotationRecordKeysOf,
  annotationSourceKey,
  buildAnnotationRecordCatalog
} from './annotations/record-catalog.js';
import {
  reconcileAnnotationRecordBindings
} from './annotations/record-selector.js';
import { validateAnnotationRecordTargets } from './annotations/validation.js';
import { classifyOptionalPositiveNumber } from '../utils/optional-positive-number.js';
import { createLosatSettings } from './losat-settings.js';
import { createAutoValueDisplay } from './auto-value-display.js';
import { createLinearRecordSelector, RECORDS_NOT_INSPECTED_LABEL } from './linear-record-selector.js';
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
  linearRecordLayoutHasSharedRow,
  linearRecordPositionTokens,
  moveLinearRecordInRow,
  planLinearSourceRowMove,
  reconcileLinearRecordLayout,
  setLinearRecordRow as updateLinearRecordRow
} from '../services/linear-record-layout.js';
import {
  describeLinearLabelVisibility,
  requireLinearLabelVisibilityMode,
  resolveLinearLabelVisibility
} from '../services/linear-label-visibility.js';
import {
  LINEAR_COMPARISON_MODES,
  LINEAR_COMPARISON_SOURCES,
  adjacentRowPairs,
  buildLinearComparisonTimeline,
  createLinearComparisonEdge,
  setComparisonEdgeFile,
  linearComparisonEdgeKey,
  materializeResolvedEdgesAsSelectedPlan,
  missingUploadPairsToResolve,
  normalizeLinearComparisonPlan,
  reconcileLinearComparisonPlan,
  resolveLinearComparisonPlan
} from '../services/linear-comparisons.js';
import {
  RECORD_SETTINGS_DEFAULTS,
  circularRecordDrawKey,
  circularRecordRequestKeys,
  circularRecordsToDraw,
  linearCardHasSettings,
  pruneRecordsOff,
  recordsWithEdits,
  removeRecordEdits
} from '../services/record-draw-selection.js';
import { ANNOTATION_RECORD_BINDING_KEY } from '../services/annotation-state.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';
import { createRecordSelection } from './record-selection.js';
import {
  projectLinearComparisonLosatModeSelection,
  projectLinearComparisonLosatpModeSelection,
  projectLinearComparisonUi
} from './comparison-ui.js';
import {
  circularDiscoveryForInput,
  discoverComparisonSequenceRecordLabel,
  discoverGffFastaRecords,
  discoverSequenceRecords
} from '../services/record-discovery.js';
import {
  applyComparisonSequenceRecordLabel,
  conservationSourceDescriptors,
  defaultConservationSeriesLabel,
  moveConservationSeriesEntry,
  normalizeFileList,
  orderedConservationSources,
  orderedOptionalConservationFiles,
  parseConservationLabelText,
  reconcileConservationSeries,
  withoutConservationSeries
} from '../services/conservation-series.js';
import {
  getDepthTrackFallbackLabel,
  getDepthTrackLabelFromFile,
  isDepthTrackAutoLabel
} from './depth-tracks.js';
import { comparisonProfileDefault, linearAxisRulerApplies } from '../mode-profiles.js';
import { msg, t, tm, tn } from '../utils/i18n.js';
import {
  IMPORTED_COMPARISON_ACTIONS,
  IMPORTED_COMPARISON_DISPOSITIONS,
  importedComparisonExecution,
  resolveImportedComparisonAction
} from '../services/imported-comparison-intent.js';
import {
  activeDepthTrackIndices,
  clearDepthTrackSourceAt,
  compactDepthFileSlots,
  depthFileSlotsFromValue,
  depthSlotTrackIndex,
  depthTrackCoverageCount,
  depthTrackMatrixWidth,
  ensureDepthTrackConfigCount as ensureDepthTrackConfigCountEntries,
  depthTrackConfigForEdit,
  isDefaultManagedDepthSlot,
  isRecordMajorDepthFileMatrix,
  normalizeDepthTrackConfig as normalizeDepthTrackConfigEntry,
  normalizeRecordMajorDepthFileRows,
  padDepthFileSlots,
  representativeDepthFiles,
  reindexDepthSlots,
  removeDepthTrackColumnAt,
  syncDepthSlotLabels,
  uploadedDepthFileCount
} from '../services/depth-track-state.js';
/** @import { HistoryManagerOptions } from '../services/history.js' */
/** @import { RulePreparationOptions } from './rule-matching.js' */
/** @import { AnnotationRecordCatalog } from './annotations/record-catalog.js' */
/** @import { SpecificColorRule } from '../services/specific-color-rules.js' */
/** @import { MountedResultContext } from './preview-runtime.js' */
/** @import { UnresolvedNotice } from '../services/feature-visibility.js' */
/** @import { FeatureListRow } from './pairwise-match-popup.js' */
/** @import { OrthogroupGroup } from './orthogroups.js' */
/** @import { SequenceRecordEntry } from '../services/record-discovery.js' */
/** @import { DepthFileSlots, RawDepthTrackConfig, ReindexDepthSlotsOptions, SyncDepthSlotLabelsOptions } from '../services/depth-track-state.js' */
/** @import { CircularTrackSlot } from '../services/circular-track-slot-model.js' */
/** @import { LinearTrackSlot } from '../services/linear-track-slot-model.js' */
/** @import { RawConservationSeriesEntry } from '../services/conservation-series.js' */
/** @import { LinearRecordRow } from '../services/linear-record-layout.js' */
/** @import { SequenceRecord } from '../services/record-discovery.js' */
/** @import { AlignmentOrientationCheckpoint, TargetDraftCheckpoint } from './record-display-options.js' */
/** @import { LosatCacheInfoEntry, MultiRecordPosition } from '../services/config.js' */
/** @import { UnvalidatedFields } from '../services/session-authority.js' */

const { onMounted, onUnmounted, watch, nextTick, computed, provide, ref, shallowRef, reactive, isRef } = window.Vue;

/** @type {Promise<typeof import('../services/export.js')> | null} */
let exportServicePromise = null;

const loadExportService = () => {
  exportServicePromise ??= import('../services/export.js');
  return exportServicePromise;
};

/**
 * @typedef {object} SessionImportRollbackOptions
 * @property {Record<'circular' | 'linear', Readonly<ArtifactSlot> | null>} artifactSlots
 *   The modes' stashed generated artifacts (E1), held by reference.
 * @property {() => Readonly<ArtifactSlot>} captureDisplayedArtifact The History snapshot service's slot capture.
 * @property {(slot: Readonly<ArtifactSlot>) => void} installDisplayedArtifact
 *   Installs a captured slot after the import rollback restored the document.
 * @property {{ circular: number }} depthTrackUiCounts
 * @property {Record<string, any>[]} depthTracks
 * @property {{ value: number }} featureListScrollTop
 * @property {{ value: HTMLElement | null }} featureListScrollRef
 * @property {{ value: string }} selectedPairwiseBlockOrthogroupId
 * @property {(() => any) | null} [captureSpecificRulePatternDrafts] The feature editor's capture of the unsaved rule pattern drafts.
 * @property {((drafts: any) => void) | null} [restoreSpecificRulePatternDrafts]
 */

/**
 * What `capture` holds for `restore`.
 * @typedef {object} SessionImportRollbackSnapshot
 * @property {SessionImportRollbackOptions['artifactSlots']} artifactSlots
 * @property {Readonly<ArtifactSlot>} displayedArtifact
 * @property {number} circularDepthTrackUiCount
 * @property {SessionImportRollbackOptions['depthTracks']} depthTracks
 * @property {number} featureListScrollTop
 * @property {string} selectedPairwiseBlockOrthogroupId
 * @property {unknown} [specificRulePatternDrafts] The unsaved rule pattern drafts, as the feature editor captured them; handed back as they are.
 */

/** @param {SessionImportRollbackOptions} options */
export const createSessionImportRollbackState = ({
  artifactSlots,
  captureDisplayedArtifact,
  installDisplayedArtifact,
  depthTrackUiCounts,
  depthTracks,
  featureListScrollTop,
  featureListScrollRef,
  selectedPairwiseBlockOrthogroupId,
  captureSpecificRulePatternDrafts = null,
  restoreSpecificRulePatternDrafts = null
}) => ({
  capture: () => ({
    artifactSlots: { ...artifactSlots },
    // The displayed artifact's parts outside state: its transport identity,
    // match-sequence owner, committed Session and CLI helper files.
    displayedArtifact: captureDisplayedArtifact(),
    circularDepthTrackUiCount: depthTrackUiCounts.circular,
    depthTracks: cloneJsonData(depthTracks),
    featureListScrollTop: featureListScrollTop.value,
    selectedPairwiseBlockOrthogroupId: selectedPairwiseBlockOrthogroupId.value,
    ...(captureSpecificRulePatternDrafts ? { specificRulePatternDrafts: captureSpecificRulePatternDrafts() } : {})
  }),
  restore: async (/** @type {SessionImportRollbackSnapshot} */ snapshot) => {
    Object.assign(artifactSlots, snapshot.artifactSlots);
    installDisplayedArtifact(snapshot.displayedArtifact);
    depthTrackUiCounts.circular = snapshot.circularDepthTrackUiCount;
    depthTracks.splice(
      0,
      depthTracks.length,
      ...cloneJsonData(snapshot.depthTracks)
    );
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

const HISTORY_RESTORE_BUSY = Object.freeze({ status: 'busy', reason: msg('Undo or Redo in progress. Retry after it finishes.') });
// E1: a History step that switches the diagram mode; its changes name `ui.mode`.
/** @param {unknown} changes */
const historyStepSwitchesMode = (changes) => (Array.isArray(changes) ? changes : [])
  .some((/** @type {{ path?: unknown[] }} */ { path } = {}) => path?.[0] === 'ui' && path[1] === 'mode');

export const createAppSetup = () => {
  provide(SETTING_NOTE_APPLIES, committedLiveSettingApplies);
  setUnmanagedConfigOverrideValidator((payload) => runDiagramHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.VALIDATE_CONFIG_OVERRIDES,
    payload
  ));
  setMainSessionComparisonFrameConverter((payload) => runDiagramHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.CONVERT_MAIN_SESSION_COMPARISON_FRAME,
    payload
  ));
  setStoredProteinResultFrameConverter((payload) => runDiagramHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.CONVERT_STORED_PROTEIN_RESULT_FRAMES,
    payload
  ));
  const {
    processing,
    processingStatus,
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
    matchSequenceRegistry,
    svgContent,
    svgResultIdentity,
    zoom,
    layoutRepositionMode,
    isPanning,
    canvasPan,
    canvasContainerRef,
    mode,
    cInputType,
    lInputType,
    files,
    selectedAnnotation,
    linearSeqs,
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
    circularRecordList,
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
    featureExtractionError,
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
  const linearTypography = createLinearTypographyController({
    adv: state.drawings.linear.adv,
    linked: state.drawings.linear.linearTypographyLinked,
    mutationAvailability: sessionOperationAvailability
  });

  const comparisonHeightValidationError = computed(() => {
    const drawing = state.activeDrawing();
    if (
      mode.value !== 'linear' ||
      drawing.linearComparisonResolution.value?.hasComparisonIntent !== true
    ) return '';
    return classifyOptionalPositiveNumber(drawing.adv.comparison_height).status === 'invalid'
      // A diagnostic error object is truthy, so it always normalizes to a model.
      ? /** @type {UserFacingError} */ (normalizeUserFacingError(diagnosticError('INPUT_INVALID', { field: 'match_height', reason: 'POSITIVE_OR_AUTO' }))).summary
      : '';
  });

  /**
   * @param {LinearComparisonEdge | undefined} left
   * @param {LinearComparisonEdge | undefined} right
   */
  const sameLinearComparisonEdge = (left, right) => (
    left?.id === right?.id &&
    left?.queryUid === right?.queryUid &&
    left?.subjectUid === right?.subjectUid &&
    left?.included === right?.included &&
    left?.fileActive === right?.fileActive &&
    left?.losatFilenameActive === right?.losatFilenameActive &&
    left?.source === right?.source &&
    left?.file === right?.file &&
    left?.losatFilename === right?.losatFilename
  );

  /** @param {DrawingState} drawing */
  const reindexLinearLosatCacheInfo = (drawing) => {
    if (!Array.isArray(losatCacheInfo.value)) return;
    // Resolution indexes are positions among the drawn records.
    const indexByUid = new Map(
      drawing.drawnLinearSeqs.value.map((sequence, index) => [String(sequence?.uid || ''), index])
    );
    const resolvedByEdgeKey = new Map(
      drawing.linearComparisonResolution.value.edges.map((/** @type {LinearComparisonResolvedEdge} */ edge) => [edge.edgeKey, edge])
    );
    losatCacheInfo.value = losatCacheInfo.value.flatMap((/** @type {LosatCacheInfoEntry} */ entry) => {
      const edgeKey = String(entry?.edgeKey || '');
      if (!edgeKey) return [entry];
      const resolved = resolvedByEdgeKey.get(edgeKey);
      const queryUid = String(resolved?.queryUid || entry?.queryUid || '');
      const subjectUid = String(resolved?.subjectUid || entry?.subjectUid || '');
      const queryIndex = indexByUid.get(queryUid);
      const subjectIndex = indexByUid.get(subjectUid);
      if (!Number.isInteger(queryIndex) || !Number.isInteger(subjectIndex)) return [];
      return [{
        ...entry,
        edgeKey,
        queryUid,
        subjectUid,
        queryIndex,
        subjectIndex,
        ordinal: Number.isInteger(Number(resolved?.ordinal))
          ? Number(resolved.ordinal)
          : entry.ordinal
      }];
    });
  };

  /** @param {DrawingState} drawing */
  const invalidateLinearComparisonArtifacts = (drawing, { preserveLosatCacheInfo = false } = {}) => {
    files.linearCanonicalComparisons = [];
    if (Array.isArray(losatCacheInfo.value)) {
      if (preserveLosatCacheInfo) reindexLinearLosatCacheInfo(drawing);
      else losatCacheInfo.value = losatCacheInfo.value.filter((/** @type {LosatCacheInfoEntry} */ entry) => !entry?.edgeKey);
    }
  };

  // OV-380: a pair set to a source that has no data yet (Upload BLAST TSV
  // without a file) compares nothing new, so it keeps the alignment. Every
  // other pair must keep its current comparison; the pending pairs count as
  // their current one.
  /** @typedef {import('../services/linear-comparisons.js').LinearComparisonResolvedEdge} ResolvedEdge */
  /** @param {DrawingState} drawing @param {import('../services/linear-comparisons.js').LinearComparisonPlan} normalized */
  const onlyPendingSourcesChanged = (drawing, normalized) => {
    /** @param {ResolvedEdge} edge */
    const pending = (edge) => edge.source === LINEAR_COMPARISON_SOURCES.UPLOAD && !(edge.fileActive && edge.file);
    /** @type {readonly ResolvedEdge[]} */
    const before = drawing.linearComparisonResolution.value.edges;
    /** @type {readonly ResolvedEdge[]} */
    const after = resolveLinearComparisonPlan({
      plan: normalized, sequences: linearSeqs, layout: effectiveLinearComparisonLayout(drawing),
      losatProgram: drawing.losatProgram.value, blastpMode: drawing.losat.blastp?.mode
    }).edges;
    if (!after.some(pending)) return false;
    const beforeByKey = new Map(before.map((edge) => [edge.edgeKey, edge]));
    const effective = after.map((edge) => (pending(edge) ? beforeByKey.get(edge.edgeKey) : edge));
    return effective.length === before.length && effective.every((edge, index) => (
      edge?.edgeKey === before[index].edgeKey && edge.source === before[index].source
      && edge.file === before[index].file && edge.losatFilename === before[index].losatFilename
    ));
  };

  /**
   * @param {DrawingState} drawing
   * @param {unknown} nextPlan An unvalidated plan; `normalizeLinearComparisonPlan` validates it.
   */
  const replaceLinearComparisonPlan = (drawing, nextPlan, { invalidate = true } = {}) => {
    const normalized = normalizeLinearComparisonPlan(nextPlan);
    const changed = normalized.mode !== drawing.linearComparisonPlan.mode
      || normalized.defaultSource !== drawing.linearComparisonPlan.defaultSource
      || normalized.edges.length !== drawing.linearComparisonPlan.edges.length
      || normalized.edges.some((edge, index) => (
        !sameLinearComparisonEdge(edge, drawing.linearComparisonPlan.edges[index])
      ));
    if (invalidate && changed && !onlyPendingSourcesChanged(drawing, normalized)) {
      similarityAlignmentActions?.clearForMutation?.(msg('comparison configuration changed.'));
    }
    drawing.linearComparisonPlan.mode = normalized.mode;
    drawing.linearComparisonPlan.defaultSource = normalized.defaultSource;
    drawing.linearComparisonPlan.edges.splice(
      0,
      drawing.linearComparisonPlan.edges.length,
      ...normalized.edges
    );
    if (invalidate) invalidateLinearComparisonArtifacts(drawing);
  };

  /**
   * @param {DrawingState} drawing
   * @param {(next: LinearComparisonPlan) => void} mutator
   */
  const mutateLinearComparisonPlan = (drawing, mutator) => history.runUndoable(msg('Change comparisons'), () => {
    const next = normalizeLinearComparisonPlan(drawing.linearComparisonPlan);
    mutator(next);
    replaceLinearComparisonPlan(drawing, next);
  });

  /** @param {DrawingState} drawing */
  const effectiveLinearComparisonLayout = (drawing) => (
    drawing.linearRecordLayoutEnabled.value ? drawing.linearRecordRows : []
  );

  /** @param {DrawingState} drawing */
  const syncLinearComparisonRecords = (drawing, { invalidate = true } = {}) => {
    const next = reconcileLinearComparisonPlan(drawing.linearComparisonPlan, linearSeqs);
    const currentEdges = drawing.linearComparisonPlan.edges;
    const unchanged = (
      next.mode === drawing.linearComparisonPlan.mode &&
      next.defaultSource === drawing.linearComparisonPlan.defaultSource &&
      next.edges.length === currentEdges.length &&
      next.edges.every((edge, index) => sameLinearComparisonEdge(edge, currentEdges[index]))
    );
    if (!unchanged) replaceLinearComparisonPlan(drawing, next, { invalidate });
    return !unchanged;
  };

  const syncLinearRecordLayout = ({ preserveLosatCacheInfo = false } = {}) => {
    const drawing = state.drawings.linear;
    const next = reconcileLinearRecordLayout(linearSeqs, drawing.linearRecordRows);
    const rowsUnchanged = next.length === drawing.linearRecordRows.length && next.every((entry, index) => (
      entry.uid === drawing.linearRecordRows[index]?.uid && entry.row === drawing.linearRecordRows[index]?.row
    ));
    if (!rowsUnchanged) drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...next);
    const comparisonsChanged = syncLinearComparisonRecords(drawing, { invalidate: false });
    if (!rowsUnchanged || comparisonsChanged) {
      invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo });
    }
    return !rowsUnchanged || comparisonsChanged;
  };
  /** @param {string} uid @param {unknown} row */
  const setLinearRecordRow = (uid, row) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    syncLinearRecordLayout({ preserveLosatCacheInfo: true });
    const previous = drawing.linearRecordRows.find((/** @type {LinearRecordRow} */ entry) => entry.uid === uid)?.row;
    updateLinearRecordRow(drawing.linearRecordRows, uid, row);
    if (drawing.linearRecordRows.find((/** @type {LinearRecordRow} */ entry) => entry.uid === uid)?.row !== previous) {
      invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo: true });
    }
  };
  /** @param {boolean} enabled */
  const setLinearRecordLayoutEnabled = (enabled) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const nextEnabled = Boolean(enabled);
    if (drawing.linearRecordLayoutEnabled.value === nextEnabled) return;
    drawing.linearRecordLayoutEnabled.value = nextEnabled;
    syncLinearRecordLayout({ preserveLosatCacheInfo: true });
    invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo: true });
    if (nextEnabled) return materializeAutomaticLinearRecords();
  };
  /** @param {string} uid @param {number} direction */
  const moveLinearRecordWithinRow = (uid, direction) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = moveLinearRecordInRow(linearSeqs, drawing.linearRecordRows, uid, direction);
    drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...next);
    syncLinearComparisonRecords(drawing, { invalidate: false });
    invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo: true });
  };

  const linearComparisonUi = computed(() => {
    const drawing = state.activeDrawing();
    return projectLinearComparisonUi({
      plan: drawing.linearComparisonPlan,
      resolution: drawing.linearComparisonResolution.value,
      adjacentEdgeKeys: adjacentRowPairs(drawing.drawnLinearSeqs.value, effectiveLinearComparisonLayout(drawing), true)
        .map(([queryUid, subjectUid]) => linearComparisonEdgeKey(queryUid, subjectUid)),
      losatProgram: drawing.losatProgram.value,
      blastpMode: drawing.losat.blastp?.mode,
      inferOrthogroups: drawing.losat.blastp?.collinearInferOrthogroups === true,
      filters: {
        min_bitscore: drawing.adv.min_bitscore,
        evalue: drawing.adv.evalue,
        identity: drawing.adv.identity,
        alignment_length: drawing.adv.alignment_length
      }
    });
  });
  // The pressed global action follows the projected intent (UJ-05).
  const linearComparisonGlobalAction = computed(() => {
    const { intentKey } = linearComparisonUi.value;
    return intentKey === 'custom' ? 'selected' : intentKey;
  });
  const canRunLinearLosat = computed(() => state.drawings.linear.drawnLinearSeqs.value.filter((sequence) => (
    lInputType.value === 'gff'
      ? sequence.gff && sequence.fasta
      : sequence.gb
  )).length >= 2);

  /** @type {{ readonly value: LinearSourceGroup[] }} */
  const linearSourceGroups = computed(() => groupLinearSourceRecords(linearSeqs));
  const linearSourceRemovalDialog = reactive({ open: false, sourceUid: '', origin: '' });
  const linearSourceRemovalReturnFocus = ref(null);
  const linearSourceRemovalTarget = computed(() => (
    linearSourceGroups.value.find((source) => source.uid === linearSourceRemovalDialog.sourceUid) || null
  ));
  const linearSourceRemovalCanDelete = computed(() => linearSourceGroups.value.length > 1);
  /** @param {{ name?: string }} file */
  const chosenFileName = (file) => String(file.name || t('Unnamed file'));
  /**
   * The file names of a Linear File, as its dialogs name it: every file the
   * row holds, which a removal takes, also one of the other input type.
   * @param {LinearSourceGroup} source
   */
  const linearSourceFileName = (source) => {
    const names = [source.sequence.gb, source.sequence.gff, source.sequence.fasta].filter(Boolean).map(chosenFileName);
    return names.length ? names.join(' + ') : t('File {n}', { n: linearSourceGroups.value.indexOf(source) + 1 });
  };
  const LINEAR_FILE_SLOT_NAMES = { gb: 'GenBank', gff: 'GFF3', fasta: 'FASTA' };
  /**
   * What a File card's header line says of the slots the card shows: each
   * file's name, "no FASTA" for a slot still empty, "No file" for none.
   * @param {LinearSourceGroup} source
   */
  const linearCardFileLine = (source) => {
    const fields = linearCardFileFields(lInputType.value);
    if (!fields.some((field) => source.sequence[field])) return t('No file');
    return fields.map((field) => (
      source.sequence[field] ? chosenFileName(source.sequence[field]) : t('no {file}', { file: LINEAR_FILE_SLOT_NAMES[field] })
    )).join(' + ');
  };
  /** @param {LinearSourceGroup} source */
  const linearCardHoldsInput = (source) => linearCardInputComplete(source, lInputType.value);
  // A File card folds to its header line. The disclosure reports every change
  // of its own state: the user's press, and a path that opened the card to
  // show a setting inside it.
  /** @param {LinearSourceGroup} source */
  const isLinearSourceFolded = (source) => linearFoldedSources.has(source.uid);
  /** @param {LinearSourceGroup} source @param {Event} event */
  const syncLinearSourceFold = (source, event) => setLinearSourceFolded(
    linearFoldedSources, source.uid, !(/** @type {HTMLDetailsElement} */ (event.currentTarget)).open
  );
  const linearSourceRemovalTargetName = computed(() => {
    const source = linearSourceRemovalTarget.value;
    return source ? linearSourceFileName(source) : t('Unavailable File');
  });
  const linearComparisonTimeline = computed(() => {
    const drawing = state.activeDrawing();
    return buildLinearComparisonTimeline({
      sequences: linearSeqs,
      recordsOff: drawing.recordsOff,
      layout: effectiveLinearComparisonLayout(drawing),
      plan: drawing.linearComparisonPlan,
      resolution: drawing.linearComparisonResolution.value
    });
  });
  /** @param {string} edgeKey */
  const linearComparisonPairForEdgeKey = (edgeKey) => {
    for (const row of linearComparisonTimeline.value.rows) {
      const pair = row.boundaryAfter?.pairs.find((/** @type {LinearComparisonPair} */ entry) => entry.edgeKey === edgeKey);
      if (pair) return pair;
    }
    return null;
  };
  // Every panel names a record the same way (OV-436).
  const linearRecordNameMap = computed(() => linearRecordNames(linearSeqs));
  /** @param {string} uid */
  const linearComparisonRecordLabel = (uid) => linearRecordNameMap.value.get(uid)?.label || t('Record');
  // The parts <record-name> lays out: the name shrinks, the accession stays whole.
  /** @param {string} uid */
  const linearComparisonRecordName = (uid) => {
    const entry = linearRecordNameMap.value.get(uid);
    return { name: entry?.name || t('Record'), detail: entry?.accession || '' };
  };
  // The settings pane's display, and the one way to bring a setting into
  // view: every "show this setting" path below ends in `revealSetting`, so a
  // setting of a section the icon rail does not show is docked first.
  const settingsPaneActions = createSettingsPaneController({ state });
  /** @param {string} disclosureKey */
  const openLinearComparisonDisclosure = async (disclosureKey) => {
    if (!disclosureKey) return null;
    const details = /** @type {HTMLDetailsElement[]} */ ([...document.querySelectorAll('[data-linear-comparison-disclosure]')])
      .find((element) => element.dataset.linearComparisonDisclosure === disclosureKey);
    if (!details) return null;
    details.open = true;
    await nextTick();
    return details;
  };
  watch([mode, () => state.drawings.linear.hasActiveLinearLosatIntent.value], (/** @type {[DiagramMode, boolean]} */ [activeMode, activeLosat]) => {
    if (activeMode === 'linear' && activeLosat) openLinearComparisonDisclosure('settings');
  }, { flush: 'post' });

  /** @param {string} edgeKey */
  const focusLinearComparisonPair = async (edgeKey) => {
    await openLinearComparisonDisclosure('selected-pairs');
    const container = /** @type {HTMLElement[]} */ ([...document.querySelectorAll('[data-edge-key]')])
      .find((element) => element.dataset.edgeKey === edgeKey);
    const control = /** @type {HTMLElement | null | undefined} */ (container?.querySelector('input, select, button'));
    if (control) await settingsPaneActions.revealSetting(control);
  };

  /** @param {{ disclosureKey: string, focusTargetKey: string, edgeKey?: string, edgeId?: string } | null | undefined} target */
  const focusLinearComparisonIssue = async (target = linearComparisonUi.value.issues[0]) => {
    if (!target) return false;
    const details = await openLinearComparisonDisclosure(target.disclosureKey);
    /** @type {HTMLElement | null} */
    let container = details;
    if (target.edgeKey) {
      container = /** @type {HTMLElement[]} */ ([...document.querySelectorAll('[data-edge-key]')])
        .find((element) => element.dataset.edgeKey === target.edgeKey) || container;
    }
    if (target.edgeId && target.focusTargetKey === 'pair-row') {
      container = /** @type {HTMLElement[]} */ ([...document.querySelectorAll('[data-linear-unplaced-draft]')])
        .find((element) => element.dataset.linearUnplacedDraft === target.edgeId) || container;
    }
    const marked = container?.querySelector(
      `[data-comparison-focus="${target.focusTargetKey}"]`
    );
    const focusable = marked?.matches?.('input, select, button, [tabindex]')
      ? marked
      : marked?.querySelector?.('input, select, button, [tabindex]');
    const control = /** @type {HTMLElement | null | undefined} */ (focusable || container?.querySelector?.('input, select, button, [tabindex]'));
    if (control) await settingsPaneActions.revealSetting(control);
    return true;
  };

  // The one way to choose a pair's BLAST TSV outside its row:
  // the pair's own uploader opens its file chooser inside the click, so the
  // file takes the uploader's availability check and History step.
  /** @param {string} edgeKey @param {{ focus?: boolean }} [options] */
  const chooseLinearComparisonPairFile = (edgeKey, { focus = false } = {}) => {
    void openLinearComparisonDisclosure('selected-pairs');
    const uploader = /** @type {HTMLElement[]} */ ([...document.querySelectorAll('[data-linear-comparison-pair-upload]')])
      .find((element) => element.dataset.linearComparisonPairUpload === edgeKey);
    const button = /** @type {HTMLElement | null | undefined} */ (uploader?.querySelector('[role="button"]'));
    if (focus && button) void settingsPaneActions.revealSetting(button);
    /** @type {HTMLInputElement | null | undefined} */ (uploader?.querySelector('input[type="file"]'))?.click();
    return Boolean(uploader);
  };

  // An action of an entry in the Comparison card's issue list; a
  // setting it changes is one Undo step.
  /** @param {LinearComparisonUiIssue} issue @param {string} actionKey */
  const runLinearComparisonIssueAction = (issue, actionKey) => {
    if (actionKey === 'all-losat') return setLinearComparisonGlobalAction('losat');
    if (actionKey === 'pairwise') {
      return history.runUndoable(msg('Change setting'), () => setLinearComparisonLosatpMode('pairwise'));
    }
    if (actionKey === 'choose-file' && chooseLinearComparisonPairFile(issue.edgeKey, { focus: true })) return true;
    return focusLinearComparisonIssue(issue);
  };

  // The Generate error's "Show in Comparison": the first issue's
  // control, the Execution setting for a LOSAT runtime failure, or the
  // Comparison card when the cause has no control.
  const showLinearComparisonCause = async () => {
    const code = normalizeUserFacingError(errorLog.value)?.code;
    const runtimeCause = code === 'LOSAT_THREADING_UNAVAILABLE' || code === 'LOSAT_RUNTIME'
      ? { disclosureKey: 'settings', focusTargetKey: 'losat-execution' } : null;
    if (await focusLinearComparisonIssue(linearComparisonUi.value.issues[0] || runtimeCause)) return;
    const card = /** @type {HTMLElement | null} */ (document.querySelector('[data-linear-comparison-card]'));
    if (!card) return;
    await settingsPaneActions.revealSetting(card, {
      align: 'start', focus: /** @type {HTMLElement | null} */ (card.querySelector('button:not([disabled])')) ?? undefined
    });
  };

  /** @param {unknown} action */
  const setLinearComparisonGlobalAction = async (action) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const normalized = String(action || '').trim().toLowerCase();
    const result = await mutateLinearComparisonPlan(drawing, (next) => {
      if (normalized === 'none') {
        next.mode = LINEAR_COMPARISON_MODES.NONE;
        return;
      }
      next.mode = LINEAR_COMPARISON_MODES.ADJACENT;
      next.defaultSource = normalized === LINEAR_COMPARISON_SOURCES.UPLOAD
        ? LINEAR_COMPARISON_SOURCES.UPLOAD
        : LINEAR_COMPARISON_SOURCES.LOSAT;
    });
    if (normalized === LINEAR_COMPARISON_SOURCES.LOSAT) {
      await openLinearComparisonDisclosure('settings');
    }
    return result;
  };

  /** @param {string} modeKey */
  const setLinearComparisonLosatMode = (modeKey) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const selection = projectLinearComparisonLosatModeSelection({ modeKey });
    if (!selection.selectable || !selection.patch) return false;
    const nextProgram = selection.patch.losatProgram;
    if (drawing.losatProgram.value === nextProgram) return true;
    similarityAlignmentActions?.clearForMutation?.(msg('comparison program changed.'));
    drawing.losatProgram.value = nextProgram;
    invalidateLinearComparisonArtifacts(drawing);
    return true;
  };

  /** @param {string} modeKey */
  const setLinearComparisonLosatpMode = (modeKey) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const selection = projectLinearComparisonLosatpModeSelection({
      plan: drawing.linearComparisonPlan,
      modeKey
    });
    if (!selection.selectable || !selection.patch) return false;
    const nextBlastpMode = selection.patch.blastpMode;
    if (drawing.losat.blastp?.mode === nextBlastpMode) return true;
    similarityAlignmentActions?.clearForMutation?.(msg('comparison mode changed.'));
    const hitLimits = drawing.losat.blastp.hitLimitsByMode ||= createDefaultLosatpHitLimits();
    hitLimits[drawing.losat.blastp.mode] = {
      candidateLimit: drawing.losat.blastp.candidateLimit,
      orthogroupMemberMaxHits: drawing.losat.blastp.orthogroupMemberMaxHits
    };
    Object.assign(drawing.losat.blastp, hitLimits[nextBlastpMode]);
    drawing.losat.blastp.mode = nextBlastpMode;
    invalidateLinearComparisonArtifacts(drawing);
    return true;
  };

  /** @param {DrawingState} drawing */
  const selectedPlanForEdit = (drawing) => {
    if (drawing.linearComparisonPlan.mode === LINEAR_COMPARISON_MODES.SELECTED) {
      return normalizeLinearComparisonPlan(drawing.linearComparisonPlan);
    }
    return materializeResolvedEdgesAsSelectedPlan(
      drawing.linearComparisonPlan,
      drawing.linearComparisonResolution.value
    );
  };

  /**
   * @param {readonly LinearComparisonEdge[]} edges
   * @param {string} id
   */
  const findEdgeIndex = (edges, id) => edges.findIndex((edge) => edge.id === id);
  /**
   * @param {readonly LinearComparisonEdge[]} edges
   * @param {Pick<LinearComparisonPair, 'edgeKey' | 'edgeId'> & Partial<Pick<LinearComparisonPair, 'draft' | 'resolved'>>} pair
   */
  const findEdgeIndexForPair = (edges, pair) => {
    const ownerId = pair?.draft?.id || pair?.resolved?.id || pair?.edgeId || '';
    const ownerIndex = ownerId ? findEdgeIndex(edges, ownerId) : -1;
    if (ownerIndex >= 0) return ownerIndex;
    const matchingIndexes = edges
      .map((edge, index) => (
        linearComparisonEdgeKey(edge.queryUid, edge.subjectUid) === pair?.edgeKey ? index : -1
      ))
      .filter((index) => index >= 0);
    return matchingIndexes.length === 1 ? matchingIndexes[0] : -1;
  };

  /**
   * @param {LinearComparisonPlan} next
   * @param {{ id?: string, queryUid: string, subjectUid: string, source?: LinearComparisonSource }} edge
   */
  const upsertSelectedComparison = (next, {
    id = '',
    queryUid,
    subjectUid,
    source = next.defaultSource
  }) => {
    const edgeKey = linearComparisonEdgeKey(queryUid, subjectUid);
    let index = id ? findEdgeIndex(next.edges, id) : -1;
    if (index < 0) {
      const matchingIndexes = next.edges
        .map((edge, edgeIndex) => (
          linearComparisonEdgeKey(edge.queryUid, edge.subjectUid) === edgeKey ? edgeIndex : -1
        ))
        .filter((edgeIndex) => edgeIndex >= 0);
      if (!id || matchingIndexes.length === 1) index = matchingIndexes[0] ?? -1;
    }
    if (index < 0) {
      next.edges.push(createLinearComparisonEdge({
        queryUid,
        subjectUid,
        included: true,
        source
      }));
      return next.edges[next.edges.length - 1];
    }
    const edge = next.edges[index];
    edge.queryUid = String(queryUid || '');
    edge.subjectUid = String(subjectUid || '');
    edge.included = true;
    edge.source = source;
    return edge;
  };

  const addLinearComparison = async () => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const drawn = drawing.drawnLinearSeqs.value;
    if (drawn.length < 2) return;
    syncLinearRecordLayout();
    const [firstPair] = adjacentRowPairs(drawn, effectiveLinearComparisonLayout(drawing));
    const [queryUid, subjectUid] = firstPair || [drawn[0].uid, drawn[1].uid];
    const next = selectedPlanForEdit(drawing);
    upsertSelectedComparison(next, { queryUid, subjectUid });
    replaceLinearComparisonPlan(drawing, next);
    await focusLinearComparisonPair(linearComparisonEdgeKey(queryUid, subjectUid));
  };
  /**
   * An omitted edge stays in the plan only to retain its file or LOSAT filename.
   * @param {LinearComparisonPlan} next
   * @param {number} index
   */
  const omitPlanEdge = (next, index) => {
    if (index < 0) return;
    const edge = next.edges[index];
    edge.included = false;
    if (!edge.file && !String(edge.losatFilename || '').trim()) next.edges.splice(index, 1);
  };
  /**
   * Sets the pairs to No comparison as one plan replacement: a pair's radio,
   * and the pairs of the BLAST TSV missing dialog (D-54).
   * @param {DrawingState} drawing
   * @param {string[]} edgeKeys
   */
  const omitLinearComparisonPairs = (drawing, edgeKeys) => {
    const next = selectedPlanForEdit(drawing);
    edgeKeys.forEach((edgeKey) => {
      const pair = linearComparisonPairForEdgeKey(edgeKey);
      if (pair) omitPlanEdge(next, findEdgeIndexForPair(next.edges, pair));
    });
    replaceLinearComparisonPlan(drawing, next);
  };
  /** @param {string} id */
  const omitLinearComparison = (id) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = selectedPlanForEdit(drawing);
    const index = findEdgeIndex(next.edges, id);
    if (index < 0) return;
    omitPlanEdge(next, index);
    replaceLinearComparisonPlan(drawing, next);
  };
  const clearSelectedLinearComparisons = () => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = selectedPlanForEdit(drawing);
    next.edges = next.edges
      .filter((edge) => edge.file || String(edge.losatFilename || '').trim())
      .map((edge) => ({ ...edge, included: false }));
    replaceLinearComparisonPlan(drawing, next);
  };
  /** @param {string} id @param {'queryUid' | 'subjectUid'} endpoint @param {unknown} uid */
  const setLinearComparisonEndpoint = (id, endpoint, uid) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!['queryUid', 'subjectUid'].includes(endpoint)) return;
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge[endpoint] = String(uid || '');
    edge.included = true;
    replaceLinearComparisonPlan(drawing, next);
  };
  /** @param {string} id @param {unknown} source */
  const setLinearComparisonSource = (id, source) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const normalized = source === LINEAR_COMPARISON_SOURCES.LOSAT
      ? LINEAR_COMPARISON_SOURCES.LOSAT
      : LINEAR_COMPARISON_SOURCES.UPLOAD;
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge.source = normalized;
    edge.included = true;
    replaceLinearComparisonPlan(drawing, next);
  };
  /** @param {string} id @param {File | null} file */
  const setLinearComparisonFile = (id, file) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    setComparisonEdgeFile(edge, file);
    edge.source = LINEAR_COMPARISON_SOURCES.UPLOAD;
    edge.included = Boolean(file) || edge.included;
    replaceLinearComparisonPlan(drawing, next);
  };
  /** @param {string} id */
  const reuseLinearComparisonFile = (id) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge?.file) return;
    edge.fileActive = true;
    edge.source = LINEAR_COMPARISON_SOURCES.UPLOAD;
    edge.included = true;
    replaceLinearComparisonPlan(drawing, next);
  };
  /** @param {string} id */
  const deactivateLinearComparisonFile = (id) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge?.file) return;
    edge.fileActive = false;
    if (edge.source === LINEAR_COMPARISON_SOURCES.UPLOAD) edge.included = false;
    replaceLinearComparisonPlan(drawing, next);
  };
  /** @param {string} id @param {unknown} value */
  const setLinearComparisonLosatFilename = (id, value) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge.losatFilename = String(value || '');
    edge.losatFilenameActive = Boolean(edge.losatFilename.trim());
    replaceLinearComparisonPlan(drawing, next, { invalidate: false });
  };
  /** @param {string} id */
  const reuseLinearComparisonLosatFilename = (id) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge || !String(edge.losatFilename || '').trim()) return;
    edge.losatFilenameActive = true;
    edge.source = LINEAR_COMPARISON_SOURCES.LOSAT;
    edge.included = true;
    replaceLinearComparisonPlan(drawing, next);
  };
  /** @param {string} id */
  const deactivateLinearComparisonLosatFilename = (id) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge.losatFilenameActive = false;
    replaceLinearComparisonPlan(drawing, next, { invalidate: false });
  };
  /**
   * @param {DrawingState} drawing
   * @param {string} edgeKey
   * @param {(edge: LinearComparisonEdge) => void} updater
   */
  const updateResolvedLosatFilenameDraft = (drawing, edgeKey, updater) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const resolved = drawing.linearComparisonResolution.value.edges.find((/** @type {LinearComparisonResolvedEdge} */ edge) => edge.edgeKey === edgeKey);
    if (!resolved) return;
    const next = normalizeLinearComparisonPlan(drawing.linearComparisonPlan);
    const pair = linearComparisonPairForEdgeKey(edgeKey) || {
      edgeKey,
      edgeId: resolved.id,
      resolved
    };
    let index = findEdgeIndexForPair(next.edges, pair);
    if (index < 0) {
      next.edges.push(createLinearComparisonEdge({
        queryUid: resolved.queryUid,
        subjectUid: resolved.subjectUid,
        included: false,
        source: /** @type {any} */ (LINEAR_COMPARISON_SOURCES.LOSAT)
      }));
      index = next.edges.length - 1;
    }
    updater(next.edges[index]);
    replaceLinearComparisonPlan(drawing, next, { invalidate: false });
  };
  /** @param {string} edgeKey @param {unknown} value */
  const setResolvedLinearComparisonLosatFilename = (edgeKey, value) => {
    const drawing = state.drawings.linear;
    updateResolvedLosatFilenameDraft(drawing, edgeKey, (edge) => {
      edge.losatFilename = String(value || '');
      edge.losatFilenameActive = Boolean(edge.losatFilename.trim());
    });
  };
  /** @param {string} edgeKey */
  const reuseResolvedLinearComparisonLosatFilename = (edgeKey) => {
    const drawing = state.drawings.linear;
    updateResolvedLosatFilenameDraft(drawing, edgeKey, (edge) => {
      if (String(edge.losatFilename || '').trim()) edge.losatFilenameActive = true;
    });
  };
  /** @param {string} edgeKey */
  const deactivateResolvedLinearComparisonLosatFilename = (edgeKey) => {
    const drawing = state.drawings.linear;
    updateResolvedLosatFilenameDraft(drawing, edgeKey, (edge) => {
      edge.losatFilenameActive = false;
    });
  };
  const addLinearComparisonBatch = (allPairs = false) => {
    const drawing = state.drawings.linear;
    syncLinearRecordLayout();
    const next = selectedPlanForEdit(drawing);
    adjacentRowPairs(drawing.drawnLinearSeqs.value, effectiveLinearComparisonLayout(drawing), allPairs).forEach(([queryUid, subjectUid]) => {
      upsertSelectedComparison(next, { queryUid, subjectUid });
    });
    replaceLinearComparisonPlan(drawing, next);
  };
  /** @param {string} edgeKey @param {string} action */
  const setLinearComparisonGapAction = (edgeKey, action) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const pair = linearComparisonPairForEdgeKey(edgeKey);
    if (!pair) return;
    if (action === 'none') {
      omitLinearComparisonPairs(drawing, [edgeKey]);
      return;
    }
    const next = selectedPlanForEdit(drawing);
    upsertSelectedComparison(next, {
      id: pair.edgeId,
      queryUid: pair.queryUid,
      subjectUid: pair.subjectUid,
      source: action === LINEAR_COMPARISON_SOURCES.UPLOAD
        ? LINEAR_COMPARISON_SOURCES.UPLOAD
        : LINEAR_COMPARISON_SOURCES.LOSAT
    });
    replaceLinearComparisonPlan(drawing, next);
  };
  /** @param {string} edgeKey @param {File | null} file */
  const setLinearComparisonCardFile = (edgeKey, file) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    let pair = linearComparisonPairForEdgeKey(edgeKey);
    let draft = pair?.draft || null;
    if (!draft) {
      setLinearComparisonGapAction(edgeKey, LINEAR_COMPARISON_SOURCES.UPLOAD);
      pair = linearComparisonPairForEdgeKey(edgeKey);
      draft = pair?.draft || null;
    }
    if (draft) setLinearComparisonFile(draft.id, file);
  };
  /** @param {string} uid @param {number} fallback */
  const linearRecordRowFor = (uid, fallback) => {
    const drawing = state.drawings.linear;
    return drawing.linearRecordRows.find((/** @type {LinearRecordRow} */ entry) => entry.uid === uid)?.row || fallback;
  };
  const linearLayoutTokens = computed(() => {
    const drawing = state.activeDrawing();
    return (
      drawing.linearRecordLayoutEnabled.value
        ? linearRecordPositionTokens(linearSeqs, drawing.linearRecordRows)
        : []
    );
  });
  const linearLosatCacheInfoByEdgeKey = computed(() => Object.fromEntries(
    (Array.isArray(losatCacheInfo.value) ? losatCacheInfo.value : [])
      .filter((/** @type {LosatCacheInfoEntry} */ entry) => String(entry?.edgeKey || ''))
      .map((/** @type {LosatCacheInfoEntry} */ entry) => [String(entry.edgeKey), entry])
  ));

  // The Linear rows whose next record discovery may expand them, each with
  // whether the row's own names were entered for the file being read. They
  // were not when the file replaced another: those names belong to that file.
  /** @type {Map<string, boolean>} */
  const pendingLinearRecordExpansions = new Map();
  // The Linear Files and the Circular input just uploaded: their first record
  // discovery infers record definitions (Linear GenBank) and opens the record
  // list of a file with many records (D-04). A Session load, a rollback, and a
  // History step never add one.
  /** @type {Set<string>} */
  const pendingLinearUploads = new Set();
  /** @type {File | null} */
  let pendingCircularUpload = null;
  /** @type {RecordListRequest[]} */
  const recordListRequests = reactive([]);
  /** @param {RecordSelectionMode} listMode @param {string} sourceKey @param {number} recordCount */
  const requestRecordListAfterUpload = (listMode, sourceKey, recordCount) => {
    if (recordCount > WEB_UX_PROFILE.recordList.autoOpenAbove) recordListRequests.push({ mode: listMode, sourceKey });
  };
  /**
   * @param {DrawingState} drawing
   * @param {{ uid: string, records: SequenceRecord[], inferDefinitions?: boolean }} discovered
   */
  const expandDiscoveredLinearRecords = (drawing, { uid, records, inferDefinitions = false }) => {
    const namesOfThisFile = pendingLinearRecordExpansions.get(uid);
    pendingLinearRecordExpansions.delete(uid);
    const index = linearSeqs.findIndex((seq) => seq.uid === uid);
    if (index < 0) return;
    const source = linearSeqs[index];
    if (namesOfThisFile === undefined || records.length < 2 || source.region_record_id) {
      if (inferDefinitions) {
        source.inferred_definition = inferredDefinitionForRecord(records, source.region_record_id);
      }
      return;
    }
    const row = linearRecordRowFor(uid, index + 1);
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
      ...linearSeqs.slice(0, index), ...expanded, ...linearSeqs.slice(index + 1)
    ], { alignmentMutation: msg('record selector changed.') });
    expanded.forEach((seq) => updateLinearRecordRow(drawing.linearRecordRows, seq.uid, row));
    return true;
  };
  /** @param {{ uid: string, records: SequenceRecord[] }} discovered */
  const handleLinearRecordsDiscovered = ({ uid, records }) => {
    const drawing = state.drawings.linear;
    const isRollbackOrSessionLoad = Boolean(
      state.sessionImportRollbackInProgress?.value ||
      state.sessionResourceDiscoveryDeferred?.value
    );
    // Only an upload infers record definitions and opens the record list; a
    // loaded Session keeps its own.
    const fromUpload = !isRollbackOrSessionLoad && pendingLinearUploads.delete(uid);
    const expanded = expandDiscoveredLinearRecords(drawing, {
      uid, records, inferDefinitions: fromUpload && lInputType.value === 'gb'
    });
    if (fromUpload) requestRecordListAfterUpload('linear', uid, records.length);
    return expanded;
  };
  /**
   * The Circular input's first discovery after its upload (success or error).
   * @param {{ primaryFile: File | null, pairedFile: File | null, recordCount: number }} discovered
   */
  const handleCircularRecordsDiscovered = ({ primaryFile, pairedFile, recordCount }) => {
    if (!pendingCircularUpload || (pendingCircularUpload !== primaryFile && pendingCircularUpload !== pairedFile)) return;
    pendingCircularUpload = null;
    requestRecordListAfterUpload('circular', 'circular', recordCount);
  };
  const materializeAutomaticLinearRecords = async () => {
    if (mode.value !== 'linear') return;
    // A row that already waits for its discovery keeps what its writer said
    // about its names (a replacing file does not own the row's record name).
    linearSeqs.forEach((seq) => {
      if (!seq.region_record_id && !pendingLinearRecordExpansions.has(seq.uid)) {
        pendingLinearRecordExpansions.set(seq.uid, true);
      }
    });
    await linearRecordSelector.refresh();
  };
  const paletteLoader = createPaletteLoader({ state });
  const linearRecordSelector = createLinearRecordSelector({
    state,
    reactive,
    onRecordsDiscovered: handleLinearRecordsDiscovered,
    recordReader: ({ inputType, primaryFile, pairedFile }) => (
      inputType === 'gff'
        ? discoverGffFastaRecords({ gffFile: primaryFile, fastaFile: pairedFile })
        : discoverSequenceRecords({ file: primaryFile, format: 'genbank' })
    )
  });
  const getCircularRecordDiscoveryState = () => circularDiscoveryForInput(state);
  const circularAnnotationSourceKey = () => annotationSourceKey({
    scope: 'circular',
    inputType: cInputType.value,
    primaryFile: cInputType.value === 'gff' ? files.c_gff : files.c_gb,
    pairedFile: cInputType.value === 'gff' ? files.c_fasta : null
  });
  /**
   * A Linear record card as the annotation record catalog reads it.
   * @param {LinearSeq} seq
   * @param {boolean} drawn
   * @returns {AnnotationCatalogSource}
   */
  const linearAnnotationSource = (seq, drawn) => {
    const primaryFile = lInputType.value === 'gff' ? seq.gff : seq.gb;
    const pairedFile = lInputType.value === 'gff' ? seq.fasta : null;
    return {
      sourceKey: annotationSourceKey({
        scope: 'linear',
        uid: seq.uid,
        inputType: lInputType.value,
        primaryFile,
        pairedFile
      }),
      selector: seq.region_record_id,
      drawn,
      hasInput: Boolean(primaryFile && (lInputType.value !== 'gff' || pairedFile)),
      status: linearRecordSelector.statusFor(seq),
      error: linearRecordSelector.errorFor(seq),
      records: linearRecordSelector.recordsFor(seq)
    };
  };
  /** @param {AnnotationCatalogSource[] | null} [linearSourcesOverride] */
  const getAnnotationRecordCatalog = (linearSourcesOverride = null) => {
    const drawing = state.drawings.circular;
    const circularDiscovery = getCircularRecordDiscoveryState();
    const circularRecordSet = circularDrawingRecordSet(circularDiscovery.records, drawing);
    const linearRecordsOff = new Set(state.drawings.linear.recordsOff);
    return buildAnnotationRecordCatalog(/** @type {any} */ ({
      mode: mode.value,
      circularSource: {
        sourceKey: circularAnnotationSourceKey(),
        hasInput: circularDiscovery.hasInput,
        status: circularDiscovery.status,
        error: circularDiscovery.error,
        records: circularRecordSet.records,
        omittedRecords: circularRecordSet.omittedRecords
      },
      linearSources: linearSourcesOverride
        || linearSeqs.map((seq) => linearAnnotationSource(seq, !linearRecordsOff.has(seq.uid)))
    }));
  };
  const previewRuntime = createPreviewRuntime({ state, serializeSvg: serializeCleanSvg });
  setPreviewRuntime(previewRuntime);

  const historyFileStore = createHistoryFileStore();
  // Before the History snapshot service, which clears its copy feedback (D-37).
  const orthogroupActions = createOrthogroupEditor({ state, getCommittedRequest: getCommittedCanonicalRenderRequest });
  const historySnapshots = createHistorySnapshotService({
    state,
    fileStore: historyFileStore,
    nextTick,
    normalizeLinearSeqList,
    buildConfigData,
    applyConfigData: applyCapturedConfigData,
    buildUiStateData,
    applyUiStateData,
    buildFeatureStateData,
    applyFeatureStateData,
    buildEditorStateData,
    applyEditorStateData,
    buildOrthogroupStateData,
    applyOrthogroupStateData,
    serializeResults,
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
  const history = createHistoryManager(/** @type {HistoryManagerOptions} */ ({
    buildIntent: historySnapshots.buildHistoryIntent,
    applyIntent: (...args) => followRuleEdits(
      state.activeDrawing(), () => restoreWithSpecificRuleDrafts(historySnapshots.applyHistoryIntent, ...args)
    ),
    buildCheckpoint: historySnapshots.buildArtifactCheckpoint,
    applyCheckpoint: (checkpoint, restored) => followRuleEdits(
      state.activeDrawing(), () => restoreWithSpecificRuleDrafts(historySnapshots.applyArtifactCheckpoint, checkpoint), restored
    ),
    captureGeneratedArtifactHandle: historySnapshots.captureGeneratedArtifactHandle,
    restoreGeneratedArtifactHandle: (...args) => restoreWithSpecificRuleDrafts(historySnapshots.restoreGeneratedArtifactHandle, ...args),
    compareGeneratedArtifactHandles: historySnapshots.compareGeneratedArtifactHandles,
    signatureFor: historySnapshots.snapshotSignature,
    fileStore: historyFileStore,
    collectCurrentFileIds: historySnapshots.collectCurrentFileIds,
    // OV-1260: dropped older steps are told, not silent.
    onHistoryLimit: (text) => showNotice({ text: tm(text), tone: 'warning' }),
    makeRef: ref,
    mutationAvailability: sessionOperationAvailability,
    // An Undo or Redo of a mode switch waits as the mode buttons do (E1).
    stepAvailability: (/** @type {{ changes?: unknown }} */ { changes }) => (
      historyStepSwitchesMode(changes) ? diagramModeOperationBusy() : null
    )
  }));
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
  const recordDisplayControls = createRecordDisplayControls({ state, computed, watch,
    linearRecordsFor: linearRecordSelector.recordsFor,
    linearRecordStatusFor: linearRecordSelector.statusFor,
    linearRecordErrorFor: linearRecordSelector.errorFor,
    runUndoable: history.runUndoable,
    getCommittedRequest: getCommittedCanonicalRenderRequest, getCommittedSession: getCommittedCanonicalSession });
  window.__GBDRAW_HISTORY__ = history;
  // R13: owners receive this port, not the record display controls. The check
  // is pure (services/feature-identity.js) over the record display's binding.
  /** @param {RecordKeyedFeature} feature */
  const isCurrentResultFeature = (feature) => isCurrentFeature(feature, recordDisplayControls.sourceBinding());
  const canUndoHistory = computed(() => {
    void history.revision.value;
    return history.canUndo();
  });
  const canRedoHistory = computed(() => {
    void history.revision.value;
    return history.canRedo();
  });
  const undoHistoryTitle = computed(() => {
    void history.revision.value;
    // A step label is History's identity (services/history.js compares it):
    // it stays English and is translated here, where it is shown.
    const label = history.undoLabel();
    return label ? t('Undo {action}', { action: tm(label) }) : t('Undo');
  });
  const redoHistoryTitle = computed(() => {
    void history.revision.value;
    const label = history.redoLabel();
    return label ? t('Redo {action}', { action: tm(label) }) : t('Redo');
  });
  const undoHistory = () => history.undo();
  const redoHistory = () => history.redo();
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
      `Reset preserved setting ${path}`,
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
    readShownLegendColor: (entry) => legendEntryColor(entry),
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
  // hides (OV-1198, OV-1215): the live compile deletes them, and the Legend
  // list leaves them out, so a hidden track's rows leave and return with it.
  /** @param {DrawingState} drawing @returns {string[]} */
  const hiddenTrackLegendRows = (drawing) => hiddenTrackLegendKeys(displayedTrackLegendRows(drawing),
    drawing === state.drawings.linear ? 'linear' : 'circular', drawing);
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
  // Python rows.
  const legendRowColors = computed(() => displayedLegendRowColors(state, state.activeDrawing(), readDrawnTracks()));
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
        showEditorIntent(legendRows
          ? { domains: [...FILL_DOMAINS, ...LIVE_EDIT_DOMAINS.legendStructure], ruleRows: legendRows }
          : { domains: FILL_DOMAINS });
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
  // `figure` is the identity of the Result of the review's mode (`mode`) that the review last
  // showed its draft on.
  // `drawerTab` is the Editor tab the review closed: Apply's new Result can leave a tab
  // unavailable for a while, which moves the shown tab (OV-1255).
  const similarityAlignmentTakeover = shallowRef(/** @type {{ drawerOpen: boolean, drawerTab: string, popup: unknown, mode: 'circular' | 'linear', figure: string } | null} */ (null));
  // The hold belongs to the review's mode (OV-1089): while the other mode is shown, its figure
  // is not the held one, so the end of the review does not give it the held view and its
  // captures keep its own view.
  const similarityAlignmentHoldShown = () => {
    const takeover = similarityAlignmentTakeover.value;
    return !takeover || takeover.mode === (mode.value === 'circular' ? 'circular' : 'linear');
  };
  // The Session writer and the artifact checkpoint save the view the pan/zoom owner holds for
  // the shown mode. The departing mode's slot is captured before a switch changes `mode`.
  setHeldPreviewViewportReader(() => (similarityAlignmentHoldShown() ? heldPreviewViewport() : null));
  const previewPopupsHeldNow = computed(() => similarityAlignmentTakeover.value !== null);
  const previewPopupsHeld = () => previewPopupsHeldNow.value;
  const featureSelection = createFeatureSelection({ state, onMounted, onUnmounted, previewPopupsHeld });
  // The palette owner comes before the feature editor, whose popup sets a
  // type's default color through it (D-15).
  const resultsManager = createResultsManager({ state, closeAfterDialogChoice: dialogChoice.closeAfterChoice });
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
    featureSelection,
    commitActiveResultEdit: previewRuntime.commitActiveResultEdit,
    showEditorIntent: (options) => showEditorIntent(options),
    selectResult,
    previewTransformInteraction,
    projectPaletteAndRules,
    projectFeatureEdits: () => projectMountedEditorIntent({ visibility: true, rerender: true, reflow: true, labels: true }),
    readUserDefaultColor: resultsManager.readUserDefaultColor,
    setDefaultColor: resultsManager.setDefaultColor,
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
    resolveOrthogroups: () => orthogroups.value.map((/** @type {OrthogroupGroup} */ group) => ({
      ...group,
      display_name: orthogroupActions.resolveOrthogroupName(group),
      description: orthogroupActions.resolveOrthogroupDescription(group)
    })),
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
    await history.captureBaseline('Initial state');
  });
  onUnmounted(() => {
    if (featureSearchDebounceId !== null) clearTimeout(featureSearchDebounceId);
    if (typeof disposeHistoryInputs === 'function') disposeHistoryInputs();
    if (window.__GBDRAW_HISTORY__ === history) delete window.__GBDRAW_HISTORY__;
    previewFeatureSearch.dispose();
    featureActions.dispose();
    disposePanZoom();
    setUnmanagedConfigOverrideValidator(null);
    setMainSessionComparisonFrameConverter(null);
    setStoredProteinResultFrameConverter(null);
    setDiagramWorkerFontSource(null);
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
  /**
   * @template T
   * @param {() => T} change
   * @returns {T}
   */
  const retireLegendStylesOfUnnamedCaptions = (change) => {
    const drawing = state.activeDrawing();
    return buildLegendStyleRetirement({
      legendColorOverrides: drawing.legendColorOverrides,
      legendStrokeOverrides: drawing.legendStrokeOverrides,
      legendEntries: drawing.legendEntries,
      dormantLegendEntries: drawing.dormantLegendEntries,
      projectLegendEntries: () => { showEditorIntent({ domains: [...LIVE_EDIT_DOMAINS.legendRows, ...LIVE_EDIT_DOMAINS.strokes] }); },
      hiddenRows: () => hiddenTrackLegendRows(drawing),
      relistLegendEntries: () => { legendActions.extractLegendEntries(); },
      namedCaptions: () => trackDataLegendCaptions({
        annotationSets: drawing.annotationSets,
        depthTracks: drawing.adv.depth_tracks,
        ...(mode.value === 'linear'
          ? { depthSlots: drawing.adv.linear_track_slots, sourcedDepthTrackIndexes: activeDepthTrackIndices(linearDepthRows()) }
          : {
              depthSlots: drawing.adv.circular_track_slots,
              sourcedDepthTrackIndexes: circularDepthRepresentatives()
                .flatMap((file, index) => (file ? [index] : []))
            })
      })
    })(change);
  };
  // The diagrams the next Circular Generate draws, from the record set of the
  // request, read once: the Auto notes take the record lengths of each, and the
  // track slot estimates those of the diagram that draws the shown Result's
  // record.
  const circularDiagrams = computed(() => circularRequestDiagrams(state.circularRecordList.value, state.drawings.circular));
  const circularEstimateRecordLengths = computed(() => circularEstimateLengths(
    circularDiagrams.value, getCommittedCanonicalRenderRequest(), state.selectedResultIndex.value
  ));
  const readCircularEstimateLengths = () => circularEstimateRecordLengths.value;
  const circularTrackSlotEditor = createCircularTrackSlotEditor({
    state, changeTrackLayout, askChoice, runUndoable: history.runUndoable, estimateRecordLengths: readCircularEstimateLengths,
    changeTrackVisibility: retireLegendStylesOfUnnamedCaptions
  });
  const linearTrackSlotEditor = createLinearTrackSlotEditor({
    state, changeTrackLayout, changeTrackVisibility: retireLegendStylesOfUnnamedCaptions
  });
  const annotationImportNotice = ref('');
  const annotationEditor = createAnnotationEditor({
    state, getRecordCatalog: getAnnotationRecordCatalog, retireLegendStylesOfUnnamedCaptions, askChoice,
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
    retireLegendStylesOfUnnamedCaptions(() => {
      if (key === 'show_gc') drawing.form.show_gc = shown;
      else if (key === 'show_skew') drawing.form.show_skew = shown;
      else drawing.form.show_depth = shown;
    });
  };
  watch(
    () => {
      const catalog = getAnnotationRecordCatalog();
      return `${catalog.status}:${catalog.signature}`;
    },
    () => {
      const drawing = state.activeDrawing();
      const catalog = getAnnotationRecordCatalog();
      if (catalog.status === 'ready' && !sessionOperationAvailability()) {
        reconcileAnnotationRecordBindings(drawing.annotationSets, catalog);
      }
    },
    { immediate: true }
  );
  const circularConservationLayoutWarning = computed(() => estimateCircularConservationLayoutWarning(
    state, state.drawings.circular, circularEstimateRecordLengths.value
  ));
  const losatSettings = createLosatSettings({ state });
  /**
   * The record a Linear entry draws, once its records are read: the one its
   * selector names, or the only record of its file.
   * @param {LinearSeq} seq
   */
  const linearSeqRecord = (seq) => {
    const records = linearRecordSelector.recordsFor(seq);
    const selection = resolveDisambiguatedRecordSelection(records, seq.region_record_id);
    return selection.record || (selection.status === 'unspecified' && records.length === 1 ? records[0] : null);
  };
  // The note beside a Linear entry's Region fields (`regionClampNote`). The
  // request draws the larger coordinate as the region's end (reverse when typed first).
  /** @param {LinearSeq} seq */
  const linearRegionNote = (seq) => {
    const region = linearCardRegion(seq, 1, linearSeqRecord(seq)?.recordLength);
    return region && !('error' in region) ? regionClampNote(region.drawn) : '';
  };
  const circularRegionNote = computed(() => circularDrawingRegionNote(state.circularRecordList.value, state.drawings.circular));
  // Auto values follow the records the next Generate draws. The Linear ones are
  // known once every drawn record is read and none is cropped; until then the
  // values name their rule.
  const autoValueDisplay = createAutoValueDisplay(state, {
    circularDiagramLengths: () => circularDiagrams.value.map((/** @type {CircularDiagramRecords} */ diagram) => diagram.lengths),
    linearRecordLengths: () => {
      const lengths = state.drawings.linear.drawnLinearSeqs.value.map((seq) => {
        if (String(seq.region_start ?? '').trim() || String(seq.region_end ?? '').trim()) return NaN;
        return Number(linearSeqRecord(seq)?.recordLength);
      });
      return lengths.length && lengths.every((length) => length > 0) ? lengths : [];
    }
  });
  const depthTrackDefaultColors = [
    '#4A90E2',
    '#E45756',
    '#2CA02C',
    '#F28E2B',
    '#9467BD',
    '#8C564B',
    '#17BECF',
    '#7F7F7F'
  ];
  /** @param {DepthFileSlots | DepthFileSlots[]} value */
  const depthFileCount = (value) => (
    isRecordMajorDepthFileMatrix(value)
      ? representativeDepthFiles(value).filter(Boolean).length
      : uploadedDepthFileCount(value)
  );
  /** @param {DepthFileSlots | DepthFileSlots[]} value */
  const depthTrackCountLabel = (value) => {
    const count = depthFileCount(value);
    return tn(count, '{count} TSV', '{count} TSVs');
  };
  const hasCircularDepthFiles = computed(() => depthFileCount(files.c_depth) > 0);
  const hasAnyLinearDepthFiles = computed(() => (
    linearSeqs.some((seq) => depthFileCount(seq?.depth) > 0)
  ));
  const canShowDepthTrack = computed(() => (
    mode.value === 'linear' ? hasAnyLinearDepthFiles.value : hasCircularDepthFiles.value
  ));
  const enabledOptionClass = 'text-slate-700 cursor-pointer';
  const disabledOptionClass = 'text-slate-400 cursor-not-allowed opacity-60';
  const depthToggleOptionClass = computed(() => (
    canShowDepthTrack.value
      ? enabledOptionClass
      : disabledOptionClass
  ));
  /** @param {LinearSeq | null | undefined} seq */
  const hasLinearDepthFiles = (seq) => depthFileCount(seq?.depth) > 0;
  const depthTrackUiCounts = reactive({
    circular: 1
  });
  const circularDepthRecordCount = () => {
    const discoveredCount = Array.isArray(circularRecordList.value)
      ? circularRecordList.value.length
      : 0;
    if (discoveredCount > 0) return discoveredCount;
    return isRecordMajorDepthFileMatrix(files.c_depth)
      ? Math.max(1, files.c_depth.length)
      : 1;
  };
  const circularDepthRows = () => normalizeRecordMajorDepthFileRows(
    files.c_depth,
    circularDepthRecordCount()
  );
  const circularDepthRepresentatives = () => representativeDepthFiles(circularDepthRows());
  /** @param {DepthFileSlots | DepthFileSlots[]} slots @param {number} [uiCount] */
  const sourceDepthTrackCount = (slots, uiCount = 1) => Math.max(
    1,
    Number(uiCount) || 1,
    representativeDepthFiles(slots).length
  );
  // The Depth panels read the drawing's series and never write them (OV-109):
  // the Depth transitions and the Session reconcile grow, trim and label the
  // list, and a series it does not hold yet shows its defaults.
  /** @param {DrawingState} drawing @param {number} count */
  const rowsForDepthTrackCount = (drawing, count) => {
    const normalizedCount = Math.max(1, Number(count) || 1);
    return Array.from({ length: normalizedCount }, (_, index) => ({
      index,
      key: `depth-track-${index}`,
      config: drawing.adv.depth_tracks[index] || normalizeDepthTrackConfig(drawing, null, index)
    }));
  };
  const linearDepthRows = () => linearSeqs.map((seq) => depthFileSlotsFromValue(seq?.depth));
  const linearDepthLogicalWidth = () => depthTrackMatrixWidth(linearDepthRows());
  const linearDepthTrackUiCount = () => Math.max(1, linearDepthLogicalWidth());
  /** @param {number} width */
  const padLinearDepthRows = (width) => {
    const targetWidth = Math.max(0, Number(width) || 0);
    linearSeqs.forEach((seq) => {
      seq.depth = padDepthFileSlots(seq.depth, targetWidth);
    });
  };
  /** @param {number} index */
  const depthTrackFallbackColor = (index) => depthTrackDefaultColors[index % depthTrackDefaultColors.length];
  /** @param {DrawingState} drawing */
  const depthTrackConfigDefaults = (drawing) => ({
    labelForIndex: getDepthTrackFallbackLabel,
    colorForIndex: depthTrackFallbackColor,
    depthColor: drawing.adv.depth_color,
    depthHeight: drawing.adv.depth_height,
    largeTickInterval: null,
    smallTickInterval: null,
    tickFontSize: null
  });
  /** @param {DrawingState} drawing @param {RawDepthTrackConfig | null | undefined} entry @param {number} index */
  const normalizeDepthTrackConfig = (drawing, entry, index) => (
    normalizeDepthTrackConfigEntry(entry, index, depthTrackConfigDefaults(drawing))
  );
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
  /** @param {'circular' | 'linear'} [drawingMode] */
  const activeDepthTrackCount = (drawingMode = mode.value) => {
    if (drawingMode === 'linear') {
      return linearDepthTrackUiCount();
    }
    return sourceDepthTrackCount(files.c_depth, depthTrackUiCounts.circular);
  };
  /** @param {DrawingState} drawing */
  const ensureDepthTrackConfigCount = (drawing, count = activeDepthTrackCount()) => {
    const targetCount = Math.max(1, Number(count) || 1);
    const normalized = ensureDepthTrackConfigCountEntries(
      drawing.adv.depth_tracks,
      targetCount,
      depthTrackConfigDefaults(drawing)
    );
    drawing.adv.depth_tracks.splice(0, drawing.adv.depth_tracks.length, ...normalized);
  };
  const circularDepthTrackRows = computed(() => {
    const drawing = state.drawings.circular;
    return rowsForDepthTrackCount(
      drawing,
      sourceDepthTrackCount(files.c_depth, depthTrackUiCounts.circular)
    );
  });
  const linearDepthTrackRows = () => {
    const drawing = state.drawings.linear;
    return rowsForDepthTrackCount(drawing, linearDepthTrackUiCount());
  };
  /** @param {LinearSourceGroup} source */
  const linearSourceDepthRows = (source) => linearDepthTrackRows().map((track) => ({
    ...track,
    status: linearSourceDepthStatus(source, track.index)
  }));
  /** @param {LinearSourceGroup} source */
  const linearSourceDepthSummary = (source) => {
    const tracks = linearSourceDepthRows(source);
    const attachedCount = tracks.filter((track) => track.status.selectedCount > 0).length;
    if (attachedCount === 0) return t('No depth track attached');
    return tn(attachedCount, '{count} depth track attached', '{count} depth tracks attached');
  };
  // The problems a File card marks, for the open card and its folded line:
  // a file whose inspection failed, and a record that TLOSATX cannot search
  // without a translation table.
  /** @param {LinearSourceGroup} source */
  const linearSourceInspectionFailed = (source) => Boolean(linearRecordSelector.errorFor(source.sequence));
  /** @param {LinearSeq} seq */
  const linearRecordGencodeRequired = (seq) => {
    const drawing = state.drawings.linear;
    return drawing.hasActiveLinearLosatIntent.value && drawing.losatProgram.value === 'tblastx'
      && !(Number(seq.losat_gencode) >= 1);
  };
  /** @param {LinearSourceGroup} source */
  const linearSourceGencodeRequired = (source) => source.records.some(({ sequence }) => linearRecordGencodeRequired(sequence));
  const depthTrackRows = computed(() => {
    const drawing = state.activeDrawing();
    return rowsForDepthTrackCount(drawing, activeDepthTrackCount());
  });
  /** @param {number} trackIndex */
  const linearDepthTrackCoverageLabel = (trackIndex) => {
    const covered = depthTrackCoverageCount(linearDepthRows(), trackIndex);
    const total = linearSeqs.length;
    return tn(total, '{covered}/{count} record', '{covered}/{count} records', { covered });
  };
  const linearDepthTrackIndexOptions = () => {
    const rows = linearDepthRows();
    const active = new Set(activeDepthTrackIndices(rows));
    return Array.from({ length: linearDepthLogicalWidth() }, (_, trackIndex) => ({
      trackIndex,
      label: getDepthTrackLabel(trackIndex),
      coverage: linearDepthTrackCoverageLabel(trackIndex),
      disabled: !active.has(trackIndex)
    }));
  };
  const definitionLineStyleRows = Object.freeze([
    { key: 'name', get label() { return t('Name / Species'); } },
    { key: 'subtitle', get label() { return t('Subtitle'); } },
    { key: 'replicon', get label() { return t('Replicon'); }, visibilityKey: 'linear_show_replicon', visibilityType: 'boolean' },
    { key: 'accession', get label() { return t('Accession'); }, visibilityKey: 'linear_accession_visibility', visibilityType: 'mode' },
    { key: 'length', get label() { return t('Length / Coordinates'); }, visibilityKey: 'linear_length_visibility', visibilityType: 'mode' }
  ]);
  const linearLabelHasSharedRow = computed(() => {
    const drawing = state.activeDrawing();
    return linearRecordLayoutHasSharedRow(
      drawing.drawnLinearSeqs.value,
      drawing.linearRecordRows,
      { enabled: Boolean(drawing.linearRecordLayoutEnabled.value) }
    );
  });
  /** @param {unknown} mode */
  const linearLabelVisibilitySummary = (mode) => describeLinearLabelVisibility(mode, {
    hasSharedRow: linearLabelHasSharedRow.value
  });
  const linearLabelAutoFields = computed(() => {
    const drawing = state.activeDrawing();
    return definitionLineStyleRows.filter((row) => (
      row.visibilityType === 'mode'
      && requireLinearLabelVisibilityMode(drawing.adv[row.visibilityKey]) === 'auto'
    ));
  });
  const linearLabelAutoDisclosure = computed(() => {
    if (!linearLabelAutoFields.value.length) return '';
    const fields = linearLabelAutoFields.value.map((/** @type {{ label: string }} */ row) => row.label).join(' and ');
    const shown = resolveLinearLabelVisibility('auto', {
      hasSharedRow: linearLabelHasSharedRow.value
    });
    // The shown Result already has Auto's outcome: say what it shows.
    const applied = getCommittedLinearDefinitionVisibility();
    const effect = linearLabelAutoFields.value
      .every((/** @type {{ visibilityKey: string }} */ row) => applied[row.visibilityKey] === shown)
      ? (shown ? t('Auto shows these fields throughout the diagram') : t('Auto hides these fields throughout the diagram'))
      : (shown
        ? t('Auto will show these fields throughout the diagram on the next successful Generate')
        : t('Auto will hide these fields throughout the diagram on the next successful Generate'));
    return shown
      ? t('{fields}: {effect} because no rendered row contains multiple records.', { fields, effect })
      : t('{fields}: {effect} because at least one rendered row contains multiple records.', { fields, effect })
        + ' ' + t('Choose Show in Record Labels to keep a field visible.');
  });
  /** @param {string} key */
  const focusLinearLabelVisibility = async (key) => {
    if (mode.value !== 'linear') return;
    const select = document.getElementById(`linear-label-visibility-${key}`);
    if (select) await settingsPaneActions.revealSetting(select);
  };
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
  /** @param {DrawingState} drawing @param {string} kind */
  const ensureDefinitionLineStyle = (drawing, kind) => {
    const key = String(kind || '');
    if (
      !drawing.adv.linear_definition_line_styles ||
      typeof drawing.adv.linear_definition_line_styles !== 'object' ||
      Array.isArray(drawing.adv.linear_definition_line_styles)
    ) {
      drawing.adv.linear_definition_line_styles = {};
    }
    const existing = drawing.adv.linear_definition_line_styles[key];
    if (!existing || typeof existing !== 'object' || Array.isArray(existing)) {
      drawing.adv.linear_definition_line_styles[key] = {
        font_size: null,
        font_weight: null,
        fill: null
      };
    }
    return drawing.adv.linear_definition_line_styles[key];
  };
  /** @param {string} kind */
  const getDefinitionLineStyleSize = (kind) => {
    const drawing = state.activeDrawing();
    return optionalNumberInputValue(
      ensureDefinitionLineStyle(drawing, kind).font_size
    );
  };
  /** @param {string} kind @param {unknown} value */
  const setDefinitionLineStyleSize = (kind, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    setOptionalNumberInputValue(ensureDefinitionLineStyle(drawing, kind), 'font_size', value);
  };
  /** @param {string} kind */
  const getDefinitionLineStyleWeight = (kind) => {
    const drawing = state.activeDrawing();
    return ensureDefinitionLineStyle(drawing, kind).font_weight ?? '';
  };
  /** @param {string} kind @param {unknown} value */
  const setDefinitionLineStyleWeight = (kind, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const normalized = String(value || '').trim().toLowerCase();
    ensureDefinitionLineStyle(drawing, kind).font_weight = normalized === 'bold' ? 'bold' : null;
  };
  /** @param {string} kind */
  const getDefinitionLineStyleFill = (kind) => {
    const drawing = state.activeDrawing();
    return ensureDefinitionLineStyle(drawing, kind).fill ?? '';
  };
  /** @param {string} kind @param {unknown} value */
  const setDefinitionLineStyleColor = (kind, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const normalized = String(value || '').trim();
    ensureDefinitionLineStyle(drawing, kind).fill = normalized || null;
  };
  /** @param {string} kind */
  const getDefinitionLineStyleColorMode = (kind) => (
    colorValueMode(getDefinitionLineStyleFill(kind))
  );
  /** @param {string} kind */
  const getDefinitionLineStyleSwatchValue = (kind) => {
    return toNativeColorInputValue(getDefinitionLineStyleFill(kind));
  };
  /** @param {{ visibilityKey?: string, visibilityType?: string }} row */
  const isDefinitionLineStyleMuted = (row) => {
    const drawing = state.activeDrawing();
    const key = row?.visibilityKey;
    if (!key) return false;
    if (row.visibilityType === 'mode') {
      return !resolveLinearLabelVisibility(drawing.adv[key], {
        hasSharedRow: linearLabelHasSharedRow.value
      });
    }
    return drawing.adv[key] === false;
  };
  /** @param {CircularTrackSlot | LinearTrackSlot | null | undefined} slot */
  const normalizeDepthSlotTrackIndex = (slot) => {
    const rawTrackIndex = Number(slot?.params?.track_index);
    return Number.isInteger(rawTrackIndex) && rawTrackIndex >= 0 ? rawTrackIndex : 0;
  };
  // A series' settings for an edit of the shown drawing, or null when the index
  // names no series of its mode: the edit is dropped and `depth_tracks` does not
  // grow (TK-03).
  /** @param {DrawingState} drawing @param {number} index */
  const depthTrackConfigForIndex = (drawing, index) => {
    if (!Array.isArray(drawing.adv.depth_tracks)) drawing.adv.depth_tracks = [];
    return depthTrackConfigForEdit(
      drawing.adv.depth_tracks, index, activeDepthTrackCount(), depthTrackConfigDefaults(drawing)
    );
  };
  // A series' saved settings, or its defaults while the drawing holds none (read only, OV-109).
  /** @param {DrawingState} drawing @param {number} index */
  const readDepthTrackConfig = (drawing, index) => (
    drawing.adv.depth_tracks[index] || normalizeDepthTrackConfig(drawing, null, index)
  );
  /** @param {number} index */
  const getDepthTrackLabel = (index) => {
    const config = readDepthTrackConfig(state.activeDrawing(), Math.max(0, Number(index) || 0));
    return String(config?.label ?? '');
  };
  /** @param {number} index */
  const getDepthTrackColor = (index) => {
    const idx = Math.max(0, Number(index) || 0);
    const config = readDepthTrackConfig(state.activeDrawing(), idx);
    return String(config?.color || depthTrackFallbackColor(idx));
  };
  /** @param {number} index @param {unknown} value */
  const setDepthTrackColor = (index, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const idx = Math.max(0, Number(index) || 0);
    const color = String(value ?? '').trim();
    const config = depthTrackConfigForIndex(drawing, idx);
    if (!config) return;
    config.color = color || depthTrackFallbackColor(idx);
  };
  /** @param {DrawingState} drawing */
  const depthTrackSlotCollections = (drawing) => [
    Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots : [],
    Array.isArray(drawing.adv.linear_track_slots) ? drawing.adv.linear_track_slots : []
  ];
  /** @param {DrawingState} drawing @param {number} index */
  const syncDepthTrackSlotLabelsForTrack = (drawing, index) => {
    depthTrackSlotCollections(drawing).forEach((slots) => {
      syncDepthSlotLabels(/** @type {SyncDepthSlotLabelsOptions} */ ({
        slots,
        depthTracks: drawing.adv.depth_tracks,
        activeCount: drawing.adv.depth_tracks.length
      }));
    });
    void index;
  };
  /** @param {number} index @param {unknown} value */
  const setDepthTrackLabel = (index, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const idx = Math.max(0, Number(index) || 0);
    const config = depthTrackConfigForIndex(drawing, idx);
    if (!config) return;
    retireLegendStylesOfUnnamedCaptions(() => {
      config.label = String(value ?? '');
      syncDepthTrackSlotLabelsForTrack(drawing, idx);
    });
  };
  /** @param {CircularTrackSlot | LinearTrackSlot | null | undefined} slot */
  const getDepthTrackLegendLabelForSlot = (slot) => (
    getDepthTrackLabel(normalizeDepthSlotTrackIndex(slot))
  );
  /** @param {CircularTrackSlot | LinearTrackSlot | null | undefined} slot @param {unknown} value */
  const setDepthTrackLegendLabelForSlot = (slot, value) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!slot) return;
    const idx = normalizeDepthSlotTrackIndex(slot);
    // A slot whose index names no series has no series label to edit (TK-03).
    if (!depthTrackConfigForIndex(state.activeDrawing(), idx)) return;
    retireLegendStylesOfUnnamedCaptions(() => {
      slot.params = slot.params && typeof slot.params === 'object' ? { ...slot.params } : {};
      const label = String(value ?? '');
      if (label.trim()) {
        slot.params.legend_label = label;
      } else {
        delete slot.params.legend_label;
      }
      setDepthTrackLabel(idx, label);
    });
  };
  /** @param {LinearTrackSlot | null | undefined} slot */
  const syncDepthTrackSlotLabel = (slot) => {
    const drawing = state.activeDrawing();
    if (!slot || slot.renderer !== 'depth') return;
    const trackIndex = normalizeDepthSlotTrackIndex(slot);
    const hasSource = mode.value === 'linear'
      ? activeDepthTrackIndices(linearDepthRows()).includes(trackIndex)
      : Boolean(circularDepthRepresentatives()[trackIndex]);
    if (hasSource) delete slot.depth_binding_error;
    syncDepthTrackSlotLabelsForTrack(drawing, trackIndex);
  };
  /** @type {string[]} */
  const depthTrackAutoLabels = [];
  /** @param {DrawingState} drawing @param {DepthFileSlots} previousFiles @param {DepthFileSlots} nextFiles @param {number} removedIndex */
  const refreshDepthTrackLabelsAfterRemoval = (drawing, previousFiles, nextFiles, removedIndex) => {
    depthFileSlotsFromValue(nextFiles).forEach((file, newIndex) => {
      const oldIndex = newIndex >= removedIndex ? newIndex + 1 : newIndex;
      const config = drawing.adv.depth_tracks[newIndex];
      if (!config) return;
      const oldFile = depthFileSlotsFromValue(previousFiles)[oldIndex] || null;
      const currentLabel = String(config.label ?? '').trim();
      if (
        isDepthTrackAutoLabel(currentLabel, oldIndex, oldFile) ||
        currentLabel === depthTrackAutoLabels[oldIndex]
      ) {
        const nextLabel = file
          ? getDepthTrackLabelFromFile(file, newIndex)
          : getDepthTrackFallbackLabel(newIndex);
        config.label = nextLabel;
        depthTrackAutoLabels[newIndex] = nextLabel;
      } else {
        depthTrackAutoLabels[newIndex] = currentLabel;
      }
    });
    depthTrackAutoLabels.length = depthFileSlotsFromValue(nextFiles).length;
  };
  /**
   * @param {DrawingState} drawing
   * @param {number} index
   * @param {FileLike | null} file
   * @param {FileLike | null} [previousFile]
   */
  const updateDepthTrackLabelFromFile = (drawing, index, file, previousFile = null) => {
    if (!file) return;
    const config = drawing.adv.depth_tracks[index];
    if (!config) return;
    const currentLabel = String(config.label ?? '').trim();
    if (isDepthTrackAutoLabel(currentLabel, index, previousFile) || currentLabel === depthTrackAutoLabels[index]) {
      const nextLabel = getDepthTrackLabelFromFile(file, index);
      config.label = nextLabel;
      depthTrackAutoLabels[index] = nextLabel;
      syncDepthTrackSlotLabelsForTrack(drawing, index);
    }
  };
  /** @param {number} index */
  const getCircularDepthFile = (index) => circularDepthRepresentatives()[Number(index)] || null;
  /** @param {number} index @param {FileLike | null} file */
  const setCircularDepthFile = (index, file) => {
    const drawing = state.drawings.circular;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const idx = Math.max(0, Number(index) || 0);
    ensureDepthTrackConfigCount(drawing, idx + 1);
    depthTrackUiCounts.circular = Math.max(depthTrackUiCounts.circular, idx + 1);
    const rows = circularDepthRows();
    const previousFile = circularDepthRepresentatives()[idx] || null;
    // A new file can rename an auto-labeled series, so the label follows in the same retirement.
    retireLegendStylesOfUnnamedCaptions(() => {
      circularTrackSlotEditor.changeCircularDepthSources(() => {
        rows.forEach((row) => {
          row[idx] = file || null;
        });
        files.c_depth = rows.map((row) => compactDepthFileSlots(row));
      });
      if (file) {
        updateDepthTrackLabelFromFile(drawing, idx, file, previousFile);
        drawing.form.show_depth = true;
      }
    });
  };
  /** @param {LinearSeq | null | undefined} seq @param {number} index */
  const getLinearDepthFile = (seq, index) => depthFileSlotsFromValue(seq?.depth)[Number(index)] || null;
  /** @param {DrawingState} drawing @param {LinearSeq[]} sequences @param {number} index @param {FileLike | null} file */
  const setLinearDepthFiles = (drawing, sequences, index, file) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const targets = Array.from(sequences || []).filter(Boolean);
    if (!targets.length) return;
    const idx = Math.max(0, Number(index) || 0);
    const logicalWidth = Math.max(linearDepthLogicalWidth(), idx + 1);
    padLinearDepthRows(logicalWidth);
    ensureDepthTrackConfigCount(drawing, logicalWidth);
    const previousFile = getLinearDepthFile(targets[0], idx);
    retireLegendStylesOfUnnamedCaptions(() => {
      linearTrackSlotEditor.changeLinearDepthSources(() => {
        targets.forEach((seq) => {
          const slots = depthFileSlotsFromValue(seq.depth);
          if (file) {
            slots[idx] = file;
            seq.depth = slots;
          } else {
            seq.depth = clearDepthTrackSourceAt(slots, idx, logicalWidth);
          }
        });
      });
      if (file) {
        updateDepthTrackLabelFromFile(drawing, idx, file, previousFile);
        drawing.form.show_depth = true;
      }
    });
  };
  /** @param {LinearSeq} seq @param {number} index @param {FileLike | null} file */
  const setLinearDepthFile = (seq, index, file) => {
    const drawing = state.drawings.linear;
    return setLinearDepthFiles(drawing, [seq], index, file);
  };
  /** @param {LinearSourceGroup} source @param {number} index @param {FileLike | null} file */
  const setLinearSourceDepthFile = (source, index, file) => {
    const drawing = state.drawings.linear;
    return setLinearDepthFiles(
      drawing,
      (source?.records || []).map(({ sequence }) => sequence),
      index,
      file
    );
  };
  /** @param {LinearSourceGroup} source @param {number} index */
  const clearLinearSourceDepthFile = (source, index) => history.runUndoable(
    msg('Clear File Depth TSV'),
    () => setLinearSourceDepthFile(source, index, null)
  );
  const addCircularDepthTrack = () => {
    const drawing = state.drawings.circular;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    depthTrackUiCounts.circular = sourceDepthTrackCount(files.c_depth, depthTrackUiCounts.circular) + 1;
    ensureDepthTrackConfigCount(drawing, depthTrackUiCounts.circular);
    retireLegendStylesOfUnnamedCaptions(() => { if (hasCircularDepthFiles.value) drawing.form.show_depth = true; });
  };
  const addLinearDepthTrack = () => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const nextCount = linearDepthTrackUiCount() + 1;
    padLinearDepthRows(nextCount);
    ensureDepthTrackConfigCount(drawing, nextCount);
    retireLegendStylesOfUnnamedCaptions(() => { if (hasAnyLinearDepthFiles.value) drawing.form.show_depth = true; });
  };
  /** @param {number} index */
  const removeCircularDepthTrack = (index) => {
    const drawing = state.drawings.circular;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0) return;
    // Removing a series retires the Legend styles of its rows (OV-65).
    retireLegendStylesOfUnnamedCaptions(() => {
      const count = sourceDepthTrackCount(files.c_depth, depthTrackUiCounts.circular);
      const previousFiles = circularDepthRepresentatives();
      files.c_depth = removeDepthTrackColumnAt(circularDepthRows(), idx)
        .map((row) => compactDepthFileSlots(row));
      if (idx < drawing.adv.depth_tracks.length) drawing.adv.depth_tracks.splice(idx, 1);
      // Show Depth goes off with the drawing's last Depth source (OV-82).
      if (!hasCircularDepthFiles.value) drawing.form.show_depth = false;
      depthTrackUiCounts.circular = count <= 1 ? 1 : Math.max(1, count - 1);
      refreshDepthTrackLabelsAfterRemoval(drawing, previousFiles, circularDepthRepresentatives(), idx);
      ensureDepthTrackConfigCount(drawing, activeDepthTrackCount('circular'));
      const activeFileCount = circularDepthRepresentatives().length;
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
          activeCount: activeFileCount,
          managedPredicate: isDefaultManagedDepthSlot
        }))
      );
      if (removedBeforeAxis) drawing.adv.circular_track_slots_axis_index = axis - removedBeforeAxis;
      syncDepthTrackSlotLabelsForTrack(drawing, idx);
      circularTrackSlotEditor.normalizeCircularTrackSlots();
    });
  };
  /** @param {number} index */
  const removeLinearDepthTrack = (index) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0) return;
    const logicalWidth = linearDepthLogicalWidth();
    if (idx >= logicalWidth) return;
    // Removing a series retires the Legend styles of its rows (OV-65).
    retireLegendStylesOfUnnamedCaptions(() => {
      const nextRows = removeDepthTrackColumnAt(linearDepthRows(), idx);
      linearSeqs.forEach((seq, recordIndex) => {
        seq.depth = nextRows[recordIndex] || [];
      });
      if (idx < drawing.adv.depth_tracks.length) drawing.adv.depth_tracks.splice(idx, 1);
      // Show Depth goes off with the drawing's last Depth source (OV-82).
      if (!hasAnyLinearDepthFiles.value) drawing.form.show_depth = false;
      depthTrackAutoLabels.splice(idx, 1);
      ensureDepthTrackConfigCount(drawing, activeDepthTrackCount('linear'));
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
          activeCount: Math.max(0, logicalWidth - 1),
          managedPredicate: isDefaultManagedDepthSlot
        }))
      );
      if (Number.isInteger(previousAxisIndex)) {
        drawing.adv.linear_track_slots_axis_index = Math.max(
          0,
          previousAxisIndex - removedManagedSlotCountBeforeAxis
        );
      }
      syncDepthTrackSlotLabelsForTrack(drawing, 0);
      linearTrackSlotEditor.syncLinearDepthSlotHeightsFromDepthTracks();
      linearTrackSlotEditor.normalizeLinearTrackSlots();
    });
  };
  watch(
    () => [
      files.c_depth,
      linearSeqs.map((seq) => depthFileSlotsFromValue(seq.depth).length).join(','),
      linearSeqs.map((seq) => seq.uid).join(','),
      depthTrackUiCounts.circular
    ],
    () => {
      // OV-1164: a Load fits both drawings' series to the loaded sources and
      // then sets the row count; the previous document's count does not pad them.
      if (sessionImportPending.value) return;
      depthTrackUiCounts.circular = Math.max(
        depthTrackUiCounts.circular,
        sourceDepthTrackCount(files.c_depth, 1)
      );
      // Each drawing holds a series for each of its own mode's Depth sources.
      ensureDepthTrackConfigCount(state.drawings.circular, activeDepthTrackCount('circular'));
      ensureDepthTrackConfigCount(state.drawings.linear, activeDepthTrackCount('linear'));
    },
    { deep: true, immediate: true }
  );
  /** @param {DrawingState} drawing */
  const isCircularConservationUploadSource = (drawing) => (
    String(drawing.circularConservation.source || '').trim().toLowerCase() === 'upload'
  );
  /** @param {DrawingState} drawing */
  const isDerivedCircularConservationReplay = (drawing) => (
    !isCircularConservationUploadSource(drawing) &&
    files.c_conservation_blasts_source === 'losat-cache' &&
    normalizeFileList(files.c_conservation_blasts).length > 0
  );
  /** @param {DrawingState} drawing */
  const getCircularConservationSourceFiles = (drawing) => (
    isCircularConservationUploadSource(drawing) || isDerivedCircularConservationReplay(drawing)
      ? normalizeFileList(files.c_conservation_blasts)
      : normalizeFileList(files.c_conservation_fastas)
  );
  /** @param {DrawingState} drawing */
  const syncCircularConservationEnabled = (drawing, sourceFiles = getCircularConservationSourceFiles(drawing)) => {
    drawing.circularConservation.enabled = normalizeFileList(sourceFiles).length > 0;
  };
  const clearDerivedCircularConservationBlasts = () => {
    if (files.c_conservation_blasts_source !== 'losat-cache') return;
    files.c_conservation_blasts = [];
    files.c_conservation_blasts_source = null;
  };
  /** @param {DrawingState} drawing @param {FileLike | FileLike[] | null | undefined} nextFiles */
  const setCircularConservationSourceFiles = (drawing, nextFiles) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const normalized = normalizeFileList(nextFiles);
    if (isCircularConservationUploadSource(drawing)) {
      files.c_conservation_blasts = normalized;
      files.c_conservation_blasts_source = null;
    } else {
      clearDerivedCircularConservationBlasts();
      files.c_conservation_fastas = normalized;
    }
    losatCacheInfo.value = [];
    syncCircularConservationSeries();
  };
  /** @param {FileLike | FileLike[] | null | undefined} nextFiles */
  const setCircularConservationUploadFiles = (nextFiles) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    files.c_conservation_blasts = normalizeFileList(nextFiles);
    files.c_conservation_blasts_source = null;
    files.c_conservation_sequence_sources = [];
    losatCacheInfo.value = [];
    syncCircularConservationSeries();
  };
  const syncCircularConservationSeries = () => {
    const drawing = state.drawings.circular;
    const sourceFiles = getCircularConservationSourceFiles(drawing);
    if (isDerivedCircularConservationReplay(drawing)) {
      drawing.circularConservation.enabled = true;
      if (drawing.adv.circular_track_slots_enabled === true) {
        circularTrackSlotEditor.syncCircularConservationSlots();
      }
      return;
    }
    syncCircularConservationEnabled(drawing, sourceFiles);
    const legacyLabels = parseConservationLabelText(drawing.circularConservation.labels);
    const nextSeries = reconcileConservationSeries({
      sourceFiles,
      previousSeries: drawing.circularConservation.series,
      legacyLabels
    });
    drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, ...nextSeries);
    if (drawing.adv.circular_track_slots_enabled === true) {
      circularTrackSlotEditor.syncCircularConservationSlots();
    }
  };
  const circularConservationSeriesRows = computed(() => {
    const drawing = state.activeDrawing();
    return (Array.isArray(drawing.circularConservation.series) ? drawing.circularConservation.series : []).map((/** @type {RawConservationSeriesEntry} */ entry, /** @type {number} */ index) => ({
      index,
      filename: String(entry?.fileName || `source_${Number(index) + 1}`).trim(),
      sourceLabel: isCircularConservationUploadSource(drawing)
        ? t('BLAST {n}', { n: Number(index) + 1 })
        : t('Comparison {n}', { n: Number(index) + 1 }),
      sourceIndex: Number.isInteger(Number(entry?.sourceIndex)) ? Number(entry.sourceIndex) : index,
      comparisonSequenceFilename: String(
        files.c_conservation_sequence_sources?.[
          Number.isInteger(Number(entry?.sourceIndex)) ? Number(entry.sourceIndex) : index
        ]?.name || ''
      ),
      defaultLabel: defaultConservationSeriesLabel(
        { name: entry?.fileName },
        Number.isInteger(Number(entry?.sourceIndex)) ? Number(entry.sourceIndex) : index
      )
    }));
  });
  /** @param {number} index @param {number} direction */
  const canMoveCircularConservationSeries = (index, direction) => {
    const drawing = state.drawings.circular;
    const idx = Number(index);
    const target = idx + Math.sign(Number(direction));
    return (
      Array.isArray(drawing.circularConservation.series) &&
      Number.isInteger(idx) &&
      idx >= 0 &&
      idx < drawing.circularConservation.series.length &&
      target >= 0 &&
      target < drawing.circularConservation.series.length
    );
  };
  /** @param {number} index @param {number} direction */
  const moveCircularConservationSeries = (index, direction) => {
    const drawing = state.drawings.circular;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (moveConservationSeriesEntry(drawing.circularConservation.series, index, direction)) {
      if (drawing.adv.circular_track_slots_enabled === true) {
        circularTrackSlotEditor.syncCircularConservationSlots();
      }
    }
  };
  const openCircularConservationComparisonFilePicker = () => {
    circularConservationFastaInput.value?.click();
  };
  const pendingComparisonRecordLabels = new Set();
  // Also waits for a read that starts meanwhile (another Add Seq).
  const settleComparisonRecordLabels = async () => {
    while (pendingComparisonRecordLabels.size > 0) {
      await Promise.allSettled([...pendingComparisonRecordLabels]);
    }
  };
  // The auxiliary file imports still applying: the watchers add them, Save and
  // Load read them, and Generate waits for them (OV-1193).
  const auxiliaryFileImports = createAuxiliaryFileImports(reactive);
  /** @param {Event & { target: HTMLInputElement }} event */
  const addCircularConservationComparisonFile = (event) => {
    const drawing = state.drawings.circular;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const target = event?.target || null;
    const selectedFile = Array.from(target?.files || []).filter(Boolean)[0] || null;
    if (!selectedFile) return;
    // D12: the Python ring reader names a GenBank or DDBJ row after it is added
    // (a read failure keeps the file-name default; Generate reports it).
    // Generate waits for these reads before it reads the row labels.
    const pending = discoverComparisonSequenceRecordLabel({ file: selectedFile })
      .then((recordLabel) => applyComparisonSequenceRecordLabel({
        series: drawing.circularConservation.series,
        sourceFiles: getCircularConservationSourceFiles(drawing),
        file: selectedFile,
        recordLabel
      }))
      .catch(() => false)
      .finally(() => pendingComparisonRecordLabels.delete(pending));
    pendingComparisonRecordLabels.add(pending);
    // B22 (R11): the input is data-history-managed. This one Add Seq step
    // captures its before now, the row appears at once, and the step commits
    // after the label read, so Undo and Redo restore the label.
    const step = history.runUndoable(msg('Change uploaded file'), settleComparisonRecordLabels);
    clearDerivedCircularConservationBlasts();
    files.c_conservation_fastas = [...normalizeFileList(files.c_conservation_fastas), selectedFile];
    losatCacheInfo.value = [];
    syncCircularConservationSeries();
    if (target) target.value = '';
    return step;
  };
  /** @param {number} sourceIndex @param {Event & { target: HTMLInputElement }} event */
  const setCircularConservationCompanionFile = (sourceIndex, event) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const index = Number(sourceIndex);
    if (!Number.isInteger(index) || index < 0) return;
    const selectedFile = Array.from(event?.target?.files || []).filter(Boolean)[0] || null;
    const next = Array.isArray(files.c_conservation_sequence_sources)
      ? [...files.c_conservation_sequence_sources]
      : [];
    while (next.length <= index) next.push(null);
    next[index] = selectedFile;
    files.c_conservation_sequence_sources = next;
    if (event?.target) event.target.value = '';
  };
  /** @param {number} index */
  const removeCircularConservationSource = (index) => {
    const drawing = state.drawings.circular;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.circularConservation.series.length) return;
    if (isDerivedCircularConservationReplay(drawing)) {
      const orderedBlasts = orderedConservationSources(
        files.c_conservation_blasts,
        drawing.circularConservation
      ).map((entry) => entry.file);
      const orderedFastas = orderedOptionalConservationFiles(
        files.c_conservation_fastas,
        drawing.circularConservation
      );
      const orderedSequenceSources = orderedOptionalConservationFiles(
        files.c_conservation_sequence_sources,
        drawing.circularConservation
      );
      const kept = withoutConservationSeries({
        blasts: orderedBlasts,
        fastas: orderedFastas,
        sequenceSources: orderedSequenceSources,
        series: drawing.circularConservation.series
      }, idx);
      drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, ...kept.series);
      drawing.circularConservation.labels = kept.labels;
      files.c_conservation_blasts = kept.blasts;
      files.c_conservation_fastas = kept.fastas;
      files.c_conservation_sequence_sources = kept.sequenceSources;
      files.c_conservation_blasts_source = kept.blasts.length > 0
        ? 'losat-cache'
        : null;
      drawing.circularConservation.enabled = kept.enabled;
      losatCacheInfo.value = [];
      if (drawing.adv.circular_track_slots_enabled === true) {
        circularTrackSlotEditor.syncCircularConservationSlots();
      }
      return;
    }
    const entry = drawing.circularConservation.series[idx];
    const sourceFiles = getCircularConservationSourceFiles(drawing);
    const descriptors = conservationSourceDescriptors(sourceFiles);
    let sourceIndex = descriptors.findIndex((descriptor) => descriptor.sourceKey === String(entry?.sourceKey || ''));
    if (sourceIndex < 0) {
      const fileName = String(entry?.fileName || '').trim();
      sourceIndex = descriptors.findIndex((descriptor) => descriptor.fileName === fileName);
    }
    if (sourceIndex < 0 && idx < sourceFiles.length) sourceIndex = idx;
    if (sourceIndex < 0 || sourceIndex >= sourceFiles.length) return;
    if (isCircularConservationUploadSource(drawing)) {
      files.c_conservation_sequence_sources = (Array.isArray(files.c_conservation_sequence_sources)
        ? files.c_conservation_sequence_sources
        : [])
        .filter((/** @type {unknown} */ _, /** @type {number} */ fileIndex) => fileIndex !== sourceIndex);
    }
    setCircularConservationSourceFiles(drawing, sourceFiles.filter((_, fileIndex) => fileIndex !== sourceIndex));
  };
  watch(
    () => {
      const drawing = state.drawings.circular;
      return [
        drawing.circularConservation.source,
        files.c_conservation_blasts,
        files.c_conservation_blasts_source,
        files.c_conservation_fastas,
        drawing.circularConservation.labels
      ];
    },
    syncCircularConservationSeries,
    { deep: true, immediate: true }
  );
  watch(
    () => {
      const drawing = state.drawings.circular;
      return [
        drawing.adv.circular_track_slots_enabled,
        drawing.circularConservation.enabled,
        drawing.circularConservation.source,
        drawing.circularConservation.series.map((/** @type {RawConservationSeriesEntry} */ entry) => `${entry?.sourceKey || ''}:${entry?.label || ''}:${entry?.color || ''}`).join('|')
      ];
    },
    (/** @type {[unknown, ...unknown[]]} */ [slotsEnabled]) => {
      if (slotsEnabled) circularTrackSlotEditor.syncCircularConservationSlots();
    }
  );
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
  historySnapshots.registerCapture('composition', legendLayout.captureCompositionIntent);
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
  // E1: each mode's generated artifact. The displayed mode's artifact is the
  // installed state; the other mode's waits here, by reference, until the mode
  // transition installs it. Only the transition and Session Load (with its
  // rollback) write it.
  /** @type {Record<'circular' | 'linear', Readonly<ArtifactSlot> | null>} */
  const artifactSlots = { circular: null, linear: null };
  // The Results of both modes. The per-Result caches (Legend inventories,
  // retired Legend entries, projected editor state) keep a stashed Result's
  // entries while the other mode is shown.
  const liveResultIdentities = () => [
    ...results.value,
    ...(artifactSlots.circular?.values.results || []),
    ...(artifactSlots.linear?.values.results || [])
  ].map(previewRuntime.getResultIdentity);
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
      // A Result shown again by a mode switch is a selection in its own
      // mode's drawing, whose Legend edits it was last shown with (PD-OI-086).
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
      // An automatic rerender and a record transform draw the committed
      // request again, so they draw a track the draft hid since (Show Depth
      // off, its last source removed, a Reset Settings): the track and its
      // Legend rows leave the new Result at once, as they left the one it
      // replaces (OV-1198 H2), also when the track claims no row (R3 H-B). A
      // selected Result was drawn before the draft hid it; its rows left with
      // the projection above. Generate draws none, so its Result stays Python's.
      const redrawn = context.phase === 'label-reflow' || context.phase === 'target-record-transform';
      if (selecting || redrawn) svgActions.applyTrackVisibility();
      if (redrawn && hiddenTrackLegendRows(drawing).length > 0) showEditorIntent({ domains: LIVE_EDIT_DOMAINS.legendRows });
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
      // OV-1280: a Result a mode switch shows again, whose drawing changed its
      // rules while it was not shown (a file card's Remove, an Undo or Redo
      // there), is drawn again when its Legend cannot show them, as after a
      // History restore of the rules.
      const arrival = arrivalShownRules;
      arrivalShownRules = null;
      if (context.bindingOptions.modeArrival && arrival?.rules && arrival.identity === context.resultIdentity) {
        featureActions.followRestoredSpecificRules(arrival.rules);
      }
      // D-36: a label rerender pays the redraw its drawing owes when, as its
      // Result is ready, no later rerender runs and its figure shows the
      // drawing's current paint, also after an edit made since its paint was
      // committed (read before that edit is shown on it, OV-346).
      // Its drawing is the shown one: a mode switch waits for a rerender
      // (`diagramModeSwitchBusy`). D-24 (OV-731): History keeps the figure
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
  const {
    runAnalysis: runGeneratedDiagramAnalysis,
    runCommittedCanonicalCandidate,
    projectCommittedRecordTransform,
    projectCommittedSimilarityAlignment,
    cancelRunAnalysis,
    runLabelReflow,
    holdAutomaticRerender,
    requestLiveSettingApply,
    discardLiveSettingRequests,
    refreshCircularRecordOrder,
    releaseSourceInputFailure,
    downloadCliHelperFiles,
    downloadLosatCache,
    downloadLosatPair,
    setLosatPairFilename,
    clearLosatCache,
    getLosatPairDefaultName,
    captureGeneratedArtifactRuntimeState,
    restoreGeneratedArtifactRuntimeState
  } = createRunAnalysis({
    state,
    rulePreparation,
    settlePendingInputs: async () => {
      await settleComparisonRecordLabels();
      await auxiliaryFileImports.settle();
    },
    isCurrentFeature: isCurrentResultFeature,
    serializeCanonicalFiles: (comparisonPlanSnapshot, linearRecordCatalog, drawing) => (
      serializeActiveRenderFiles(state.mode.value, state, drawing, {
        comparisonPlan: comparisonPlanSnapshot,
        linearRecordCatalog
      })
    ),
    prepareLinearRecordCatalog,
    recordDisplayRows: recordDisplayControls.allRows,
    assertActiveModeInputs,
    readDraftSignature: () => history.getCurrentIntentSignature(),
    readHistoryStepSerial: () => history.intentStepSerial(),
    adoptDrawnWrite: (write, request) => history.adoptDrawnWrite(write, request),
    closeLabelTextScopeDialog: featureActions.closeLabelTextScopeDialog,
    clearLabelBuildNotices: featureActions.clearLabelBuildNotices,
    askChoice,
    showNotice,
    canonicalSessionVersion: SESSION_VERSION,
    adoptCanonicalRenderArtifacts,
    getCommittedCanonicalSession,
    readCommittedCanonicalRevision,
    commitLiveSettingValues,
    captureDecorationContinuity: legendLayout.captureDecorationContinuity,
    captureGeneratedArtifactHandle: historySnapshots.captureGeneratedArtifactHandle,
    captureGeneratedArtifactOwnerSet: historySnapshots.captureGeneratedArtifactOwnerSet,
    installGeneratedArtifactOwnerSet: historySnapshots.installGeneratedArtifactOwnerSet,
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
      const catalog = getAnnotationRecordCatalog();
      reconcileAnnotationRecordBindings(drawing.annotationSets, catalog);
      return validateAnnotationRecordTargets(drawing.annotationSets, catalog);
    },
    omittedAnnotationRecordKeys: () => getAnnotationRecordCatalog().omittedRecordKeys || [],
    onCircularRecordsDiscovered: handleCircularRecordsDiscovered
  });
  // R1 (c), PD-OI-102: a History step, a checkpoint (Reset Settings, a rule
  // preset), and their Undo and Redo that change live setting rows of the
  // shown drawing ask the automatic rerender to apply them.
  history.setSettingsStepListener((changes) => {
    requestLiveSettingApply(changes ? liveSettingRequestPaths(changes, mode.value) : allLiveSettingRequestPaths(), mode.value);
  });
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
      ? linearRecordSelector.refresh()
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
    canRetryAuxiliaryImportFailure, retryAuxiliaryImportFailure, resetModeTransientUi
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
    refreshCircularRecordOrder,
    // UJ-08: a discovery that reads every Linear source ends a Generate
    // failure about those inputs, as the Circular discovery does.
    refreshLinearRecordSelectors: async (/** @type {{ suppress?: boolean } | undefined} */ options) => {
      const outcome = await linearRecordSelector.refresh(options);
      const read = linearSeqs.every((/** @type {Record<string, any>} */ seq) => ['ready', 'error'].includes(linearRecordSelector.statusFor(seq)));
      if (read) releaseSourceInputFailure('linear');
      return outcome;
    },
    resetPreviewViewport,
    resetRightDrawer: rightDrawerActions.resetRightDrawer,
    closeLabelTextScopeDialog: featureActions.closeLabelTextScopeDialog,
    closeHiddenLabelTextDialog: featureActions.closeHiddenLabelTextDialog,
    clearLabelBuildNotices: featureActions.clearLabelBuildNotices,
    previewRuntime,
    preparePaletteDefinitions: paletteLoader.loadPaletteAsset
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
      description: t('Loading a Session replaces the current settings and Results and clears Undo history.')
        + ' ' + t('Changes since the last Save Session or Load Session are lost.'),
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
  // A loaded Session says so in the notice banner, with what the load changed
  // without a loss as a list; what it could not restore was listed in the
  // Load dialog before.
  /** @param {SessionFileInput} event */
  const importSession = async (event) => {
    const name = String(event?.target?.files?.[0]?.name || '');
    const result = await loadSessionFromFile(event);
    // The notice of an earlier operation does not stand for this load; a
    // Load the user canceled changed nothing, the notice included.
    if (name && !['canceled', 'skipped', 'stale'].includes(String(result?.status))) noticeBanner.dismiss();
    if (result?.status === 'ok') {
      // The Depth panel shows each loaded Circular series, with or without a Depth file (BA-OD94).
      depthTrackUiCounts.circular = Math.max(1, state.drawings.circular.adv.depth_tracks.length);
      showNotice({ text: name ? t('Session loaded: {file}', { file: name }) : t('Session loaded.'), items: result.notices || [] });
    } else if (result?.status === 'legacy') {
      showNotice({ text: t('Configuration loaded from an older gbdraw. Save Session to keep it in the current format.') });
    }
    return result;
  };
  /** @param {SessionFileInput} event */
  const loadSessionFromFile = (event) => importSessionFromFile(event, {
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
    rollbackState: createSessionImportRollbackState({
      artifactSlots,
      captureDisplayedArtifact: historySnapshots.captureArtifactSlot,
      // The rollback restored the pan already; the install keeps it.
      installDisplayedArtifact: (slot) => historySnapshots.installArtifactSlot(slot, {
        mode: slot.mode, ui: { canvasPan: { x: canvasPan.x, y: canvasPan.y } }
      }),
      depthTrackUiCounts,
      depthTracks: state.activeDrawing().adv.depth_tracks,
      featureListScrollTop,
      featureListScrollRef,
      selectedPairwiseBlockOrthogroupId,
      captureSpecificRulePatternDrafts: featureActions.captureSpecificRulePatternDrafts,
      restoreSpecificRulePatternDrafts: featureActions.restoreSpecificRulePatternDrafts
    }),
    afterImport: async (/** @type {{ status: string, decompressedCharacters?: number, isCurrent: () => boolean }} */ result) => {
      if (result?.status === 'ok' || result?.status === 'legacy') {
        // The loaded draft's live values wait for their next change or Generate (6.7).
        discardLiveSettingRequests();
        historySnapshots.clearGeneratedArtifactIdentity({
          retainedBytes: result?.status === 'ok'
            ? Number(result.decompressedCharacters || 0) * 2
            : 0
        });
        // A replaced document resets transient UI: the selection named features
        // of the previous Session.
        featureSelection.clearFeatureSelection({ clearStatus: true, syncDom: false });
        previewFeatureSearch.resetModeSearches();
        settingsPaneActions.closeSettingsPopup();
        await nextTick();
        recordSessionLifecycleEvent('history-baseline-start');
        if (!await history.initializeIntentBaseline('Loaded session', { isCurrent: result.isCurrent })) {
          throw new Error('Session loading canceled.');
        }
        featureActions.clearSpecificRulePatternDrafts();
        annotationImportNotice.value = '';
        recordSessionLifecycleEvent('history-baseline-end');
        if (circularRecordPresentationPanel.value) circularRecordPresentationPanel.value.open = false;
        closeLegendStrokeOptions();
      }
    }
  });
  // UJ-06 (Owner 2026-10-08): Load an example shows only where this origin
  // serves the Gallery catalog (gbdraw.app, not the pip GUI), and its chooser
  // lists every catalog example. The choice closes the chooser, asks the UJ-09
  // question, and loads the Session through the Load Session path. The catalog
  // read does not hold up startup. A loaded example opens fitted, as the Fit
  // button leaves it (Owner 2026-10-08); Load Session keeps the saved zoom.
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
        // The import resolves after the shown Result's preview is mounted.
        const { status } = await importSession({ target: { files: [file], value: '' } }) || {};
        if (status === 'ok' || status === 'legacy') {
          await nextTick();
          fitPreviewToViewport();
        }
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
    addWhitelistRule,
    removeWhitelistRule,
    editWhitelistRule,
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
    canEditFeatureColor,
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
    const shown = [...new Set([...(colors ? FILL_DOMAINS : []), ...(visibility ? VISIBILITY_DOMAINS : []), ...domains])];
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
    // E1: a step that switched the mode showed that mode's own Result through
    // the transition, which projects the editor intent onto it on display.
    if (historyStepSwitchesMode(changes)) return;
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
    // The restored form decides which track groups the mounted Result shows.
    // The visibility watcher is suppressed while a step restores files (a Depth
    // source), so the step applies it here once the restored container is
    // mounted. A step that switched the mode leaves it to the mode's own
    // restore: the mounted Result must be of the current mode (OV-66, R3).
    if (changedDomains.has('files') && changedDomains.has('config')
      && state.generatedMode.value === mode.value) {
      svgActions.applyTrackVisibility();
    }
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
      drawn
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
    drawing.deletedLegendEntries.value.map((/** @type {LegendEntry} */ entry) => entry.originalCaption || entry.caption),
    [...(drawing.addedLegendCaptions.value || [])],
    drawing.legendEntries.value.filter((/** @type {LegendEntry} */ entry) => /** @type {string} */ (entry.caption) !== entry.originalCaption)
      .map((/** @type {LegendEntry} */ entry) => [entry.originalCaption, entry.caption]),
    // The direct additions, which no Result's inventory lists.
    drawing.legendEntries.value.filter((/** @type {LegendEntry} */ entry) => !originalLegendOrder.value.includes(entry.originalCaption || entry.caption))
      .map((/** @type {LegendEntry} */ entry) => entry.caption)
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
        ? admitFeatureCatalog(catalog, toRaw(results.value), { mode: state.generatedMode.value })
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
        ruleRows,
        hiddenTrackRows: hiddenTrackLegendRows(drawing)
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
    strokes: STROKE_DOMAINS, legendFills: LEGEND_FILL_DOMAINS, fills: FILL_DOMAINS, visibility: VISIBILITY_DOMAINS
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
  // OV-1280: the rule table the Result a mode switch shows again was last
  // shown with (`resultPaintRecord`), read once it is ready (`afterReady`);
  // null rules: unknown (a restored Result).
  /** @type {{ identity: string, rules: SpecificColorRule[] | null } | null} */
  let arrivalShownRules = null;
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
    arrivalShownRules = null;
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
    if (context.bindingOptions.modeArrival) arrivalShownRules = { identity, rules: paintRuleTable(previous) };
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
      let prepared = await projectMountedEditorIntent({ colors, visibility, show: false });
      // Declined with nothing blocking (the rules changed meanwhile): prepared
      // once more here, so the display still compiles once (R14-3).
      if (prepared === false && !displayBlocked()) prepared = await projectMountedEditorIntent({ colors, visibility, show: false });
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
    if (!projects && !failed) return;
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

  // E1: shows `nextMode`'s artifact. Keeps the displayed one in its mode's
  // slot and installs the arriving slot, or an empty one. The Legend owner
  // keeps each Result's inventory (OV-47): the arriving Result's goes to it,
  // and the binder adopts it when the Result is displayed, as for a batch
  // Result.
  /** @param {Readonly<ArtifactSlot>} slot An installed slot. */
  const rememberShownResultInventory = (slot) => {
    const shown = results.value[Number(selectedResultIndex.value) || 0];
    if (shown) legendActions.rememberResultInventory(previewRuntime.getResultIdentity(shown), [...slot.legendInventory]);
  };
  /**
   * @param {'circular' | 'linear'} nextMode
   * @returns {Readonly<ArtifactSlot> | null} The installed slot (null: the mode has no artifact yet).
   */
  const swapArtifactSlots = (nextMode) => {
    previewRuntime.clearActiveRuntime();
    const departing = historySnapshots.captureArtifactSlot();
    artifactSlots[departing.mode] = departing;
    const arriving = artifactSlots[nextMode];
    historySnapshots.installArtifactSlot(arriving, { mode: nextMode });
    artifactSlots[nextMode] = null;
    if (arriving) rememberShownResultInventory(arriving);
    return arriving;
  };

  // E1 (R10, R13): the one transition between the diagram modes, in one task
  // so a switch is one History step. Undo and Redo of a switch run it through
  // the History port; no watcher writes state on a mode change. Callers check
  // availability first (`setDiagramMode`); the History port calls it inside
  // its own operation. Steps, each one call:
  //   1. Settle the departing mode: pending rule-pattern work, the feature
  //      selection, and the displayed Result's projection state (it followed
  //      every live edit until now).
  //   2. Swap the artifact slots (`swapArtifactSlots`).
  //   3. Set `mode`. Each mode has its own drawing (PD-OI-086), so the
  //      template now binds the arriving drawing and no setting is written.
  //   4. Reset the departing mode's transient UI; the Preview search is kept
  //      per mode.
  //   5. Show the arriving Result as a selection.
  /**
   * @param {'circular' | 'linear'} nextMode
   * @returns {boolean} Whether the mode changed.
   */
  const transitionDiagramMode = (nextMode) => {
    const drawing = state.activeDrawing();
    const previousMode = mode.value;
    if (nextMode === previousMode) return false;
    // 1. Settle the departing mode; `drawing` is its drawing, resolved before
    // the mode changes.
    featureActions.suspendSpecificRulePatternDrafts();
    featureSelection.clearFeatureSelection({ clearStatus: true, syncDom: false });
    rememberDepartingResultProjection(drawing);
    // 2. Artifact slots.
    swapArtifactSlots(nextMode);
    // 3. Mode.
    mode.value = nextMode;
    // 4. Transient UI. "Showing the last successful result" named the
    // departing mode's Result, and the notice the departing mode's work.
    resetModeTransientUi();
    noticeBanner.dismiss();
    previewFeatureSearch.switchModeSearch(previousMode, nextMode);
    failedGeneratePreservedResult.value = false;
    // OV-372: an operation error belongs to the mode it was raised in, except
    // one that names the arriving mode's Result (PD-OI-079).
    if (errorLog.value && normalizeUserFacingError(errorLog.value)?.context?.diagramMode !== nextMode) {
      errorLog.value = null;
    }
    // 5. Presentation: the arriving Result is shown as a selection.
    previewRuntime.presentSelectedResult({ modeArrival: true });
    return true;
  };
  historySnapshots.registerModeTransition(transitionDiagramMode);
  // The committed Session of a mode's artifact, shown or kept in its slot.
  historySnapshots.registerModeCommittedSession((targetMode) => readArtifactOfMode(targetMode, mode.value,
    artifactSlots, getCommittedCanonicalSession,
    (slot) => slot.runtimeState?.canonical?.committedCanonicalSession) ?? null);
  // A switch waits for Generate, a label rerender, Save and Load, so each
  // Result lands in its own mode's slot (OIPC-C07): the mode buttons and an
  // Undo or Redo of a switch read this one predicate. The buttons also wait
  // for a History restore or open History transaction, so a restore finds the
  // mode it was captured in.
  const diagramModeOperationBusy = () => sessionOperationAvailability('history');
  const diagramModeSwitchBusy = () => diagramModeOperationBusy()
    || history.historyAvailability()
    || (history.restoring.value || history.traversalPending() ? HISTORY_RESTORE_BUSY : null);
  const diagramModeSwitchAvailable = computed(() => !diagramModeSwitchBusy());
  /** @param {string} nextMode */
  const setDiagramMode = (nextMode) => {
    if ((nextMode !== 'circular' && nextMode !== 'linear') || nextMode === mode.value) return false;
    const busy = diagramModeSwitchBusy();
    if (busy) return busy;
    transitionDiagramMode(nextMode);
    return { status: 'ok' };
  };
  // Session Load (E1): the loaded Session's artifact slots, before its preview
  // mounts. `opening` is the other Result set's slot when the Session opens on
  // that set's mode; `stashed` waits in its mode. The stash starts again with
  // the loaded Session.
  /** @param {LoadedArtifactSlots} slots */
  const installLoadedArtifactSlots = ({ opening, stashed }) => {
    artifactSlots.circular = null;
    artifactSlots.linear = null;
    if (stashed) artifactSlots[stashed.mode] = stashed;
    if (!opening) return;
    historySnapshots.installArtifactSlot(opening, { mode: opening.mode });
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
      unmatched: countUnresolvedFeatureEdits({
        featureOverrides: drawing.featureOverrides,
        featurePlacementOverrides: drawing.featurePlacementOverrides,
        notices
      }),
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
      featureOverrides: drawing.featureOverrides,
      featurePlacementOverrides: drawing.featurePlacementOverrides,
      notices: featureIdentityNotices.value
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
      const index = deleted.findIndex((/** @type {LegendEntry} */ entry) => (entry.originalCaption || entry.caption) === key);
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
  // The label stays an unmarked literal (its Undo title stays English): a
  // registered test reads it from this line (feature-label-visual-unit.test.mjs).
  const handleHiddenLabelTextChoiceWithHistory = undoableAction('Change label visibility', handleHiddenLabelTextChoice, { joinsOpenStep: true });
  const requestLabelTextChangeByFeatureIdWithHistory = undoableAction(
    msg('Change label text'),
    requestLabelTextChangeByFeatureId
  );
  const requestLabelTextChangeByKeyWithHistory = undoableAction(msg('Change label text'), requestLabelTextChangeByKey);
  const resetAllLabelTextOverridesWithHistory = undoableAction(msg('Reset label edits'), resetAllLabelTextOverrides);
  const loadLabelOverrideTableWithHistory = undoableAction(msg('Load label edits'), loadLabelOverrideTable);
  const runInfoCopyStatus = ref('');
  const exactReplayCopyStatus = ref('');

  /** @param {{ elapsedMs?: unknown } | null | undefined} info */
  const runInfoElapsedText = (info) => formatElapsedMs(info?.elapsedMs);
  /** @param {{ reproducibility?: { level?: unknown } | null } | null | undefined} info */
  const runInfoReproducibilityText = (info) => reproducibilityLabel(info?.reproducibility?.level);
  const runInfoHasCliHelperFiles = computed(() =>
    Array.isArray(lastRunInfo.value?.helperFiles) && lastRunInfo.value.helperFiles.length > 0
  );
  /** @param {unknown} commandValue @param {{ value: string }} status */
  const copyRunInfoCommand = async (commandValue, status) => {
    const command = String(commandValue || '');
    if (!command) return;
    try {
      await copyTextToClipboard(command);
      const copied = t('Copied');
      status.value = copied;
      setTimeout(() => {
        if (status.value === copied) status.value = '';
      }, 1600);
    } catch (error) {
      console.warn('Failed to copy the requested command.', normalizeUserFacingError(error));
      const failed = t('Copy failed');
      status.value = failed;
      setTimeout(() => {
        if (status.value === failed) status.value = '';
      }, 2200);
    }
  };
  const copyRunCommand = () => copyRunInfoCommand(
    lastRunInfo.value?.sourceRecipe?.command || lastRunInfo.value?.command,
    runInfoCopyStatus
  );
  const copyExactReplayCommand = () => copyRunInfoCommand(
    lastRunInfo.value?.exactReplay?.command || lastRunInfo.value?.sessionCommand,
    exactReplayCopyStatus
  );

  /** @param {Pick<AnnotationRecordCatalog, 'issues'>} catalog */
  const catalogIssueError = (catalog) => {
    const issue = catalog.issues[0];
    return issue ? diagnosticError(issue.code, issue.context) : diagnosticError('INPUT_UNREADABLE');
  };
  async function prepareLinearRecordCatalog({ privateCandidate = false } = {}) {
    const drawing = state.drawings.linear;
    if (mode.value !== 'linear') return { catalog: null, error: '' };
    const hasAutomaticSequence = linearSeqs.some((sequence) => {
      if (String(sequence?.region_record_id || '').trim()) return false;
      const primary = lInputType.value === 'gff' ? sequence?.gff : sequence?.gb;
      return Boolean(primary && (lInputType.value !== 'gff' || sequence?.fasta));
    });
    const hasRegionAnnotations = drawing.annotationSets.some((/** @type {{ annotations?: unknown } | null | undefined} */ set) => (
      Array.isArray(set?.annotations) && set.annotations.length > 0
    ));
    if (!hasAutomaticSequence && !hasRegionAnnotations) {
      return { catalog: null, error: '' };
    }
    if (privateCandidate) {
      const linearRecordsOff = new Set(drawing.recordsOff);
      const sources = await Promise.all(linearSeqs.map(async (seq) => {
        const drawn = !linearRecordsOff.has(seq.uid);
        const inputType = lInputType.value;
        const primaryFile = inputType === 'gff' ? seq.gff : seq.gb;
        const pairedFile = inputType === 'gff' ? seq.fasta : null;
        const sourceKey = annotationSourceKey({ scope: 'linear', uid: seq.uid, inputType, primaryFile, pairedFile });
        try {
          const records = inputType === 'gff'
            ? await discoverGffFastaRecords({ gffFile: primaryFile, fastaFile: pairedFile })
            : await discoverSequenceRecords({ file: primaryFile, format: 'genbank' });
          return { sourceKey, selector: seq.region_record_id, drawn, hasInput: Boolean(primaryFile), status: 'ready', records };
        } catch (error) {
          return { sourceKey, selector: seq.region_record_id, drawn, hasInput: Boolean(primaryFile), status: 'error', error: error.message, records: [] };
        }
      }));
      const catalog = getAnnotationRecordCatalog(sources);
      return catalog.status === 'ready' ? { catalog, error: '' }
        : { catalog: null, error: catalogIssueError(catalog) };
    }
    let catalog = getAnnotationRecordCatalog();
    if (catalog.status !== 'ready') {
      try {
        await linearRecordSelector.refresh();
      } catch (error) {
        return { catalog: null, error: normalizeUserFacingError(error, { operation: 'listSequenceRecords', stage: 'helper' }) };
      }
      catalog = getAnnotationRecordCatalog();
    }
    return catalog.status === 'ready'
      ? { catalog, error: '' }
      : {
          catalog: null,
          // A stored '' (no error) has no `code`: the read is undefined, as for null.
          error: drawing.drawnLinearSeqs.value.map(seq => linearRecordSelector.errorModelFor(seq)).find(error => /** @type {UserFacingError} */ (error)?.code)
            || catalogIssueError(catalog)
        };
  }

  const runAnalysis = () => {
    const drawing = state.activeDrawing();
    const patternDrafts = featureActions.captureSpecificRulePatternDrafts();
    return runGeneratedDiagramAnalysis(null, null, null, {
      prepareGenerate: async () => {
        if (mode.value === 'linear') {
          if (drawing.importedComparisonIntent.disposition === IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE) {
            await materializeAutomaticLinearRecords();
          } else {
            await linearRecordSelector.refresh();
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
          missingBlastTsvDialog.show = true;
          return { status: 'blocked', reason: 'missing-blast-tsv' };
        }
        return { status: 'ready', comparisonPlanSnapshot, comparisonExecution };
      },
      afterGenerate: async (result) => {
        if (result?.status === 'error' && mode.value === 'linear') {
          await focusLinearComparisonIssue();
        }
        // Labels that do not fit drew nothing: the dialog computes the rows the bounded search did not reach.
        if (result?.status === 'error' && result.labelFit) {
          void openLabelFitDialog(result.labelFit.report, {
            rowToken: result.labelFit.rowToken, mode: mode.value === 'linear' ? 'linear' : 'circular'
          });
        }
        if (result?.status === 'ok') {
          featureActions.clearSpecificRulePatternDrafts();
          await rulePreparation.prepare();
          featureSelection.clearFeatureSelection({ clearStatus: true });
        } else {
          featureActions.restoreSpecificRulePatternDrafts(patternDrafts);
        }
      }
    });
  };
  // E1 (PD-OI-079): an error about another diagram mode's Result (Save of a
  // legacy Result kept in that mode) offers Generate there: the action shows
  // that mode first.
  const generateFromError = () => {
    const target = errorLog.value?.context?.diagramMode;
    if ((target === 'circular' || target === 'linear') && target !== mode.value) {
      const switched = setDiagramMode(target);
      if (switched && switched.status === 'busy') return switched;
    }
    return runAnalysis();
  };

  /** @param {DrawingState} drawing @param {unknown} action */
  const chooseImportedComparisonAction = (drawing, action) => history.runUndoable(
    `${String(action || '').toLowerCase()} imported comparison`,
    () => {
      const outcome = resolveImportedComparisonAction({
        intent: drawing.importedComparisonIntent,
        action,
        draftResolution: drawing.linearComparisonResolution.value
      });
      if (!outcome.ok) {
        errorLog.value = normalizeUserFacingError(outcome.message, { operation: 'generate', stage: 'request-validation' });
        return false;
      }
      if (outcome.action === IMPORTED_COMPARISON_ACTIONS.CLEAR) {
        replaceLinearComparisonPlan(drawing, { mode: 'none', defaultSource: 'losat', edges: [] });
      }
      errorLog.value = null;
      return true;
    }
  );
  const inheritImportedComparison = () => {
    const drawing = state.activeDrawing();
    return chooseImportedComparisonAction(
      drawing,
      IMPORTED_COMPARISON_ACTIONS.INHERIT
    );
  };
  const replaceImportedComparison = () => {
    const drawing = state.activeDrawing();
    return chooseImportedComparisonAction(
      drawing,
      IMPORTED_COMPARISON_ACTIONS.REPLACE
    );
  };
  const clearImportedComparison = () => {
    const drawing = state.activeDrawing();
    return chooseImportedComparisonAction(
      drawing,
      IMPORTED_COMPARISON_ACTIONS.CLEAR
    );
  };
  const importedComparisonNeedsResolution = computed(() => (
    state.activeDrawing().importedComparisonIntent.disposition !== IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE
  ));
  const importedComparisonCanInherit = computed(() => (
    state.activeDrawing().importedComparisonIntent.disposition
      === IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY
  ));
  const importedComparisonCanReplace = computed(() => {
    const drawing = state.activeDrawing();
    return (
      drawing.linearComparisonResolution.value.valid
      && drawing.linearComparisonResolution.value.hasComparisonIntent
    );
  });

  const cancelGeneration = () => cancelRunAnalysis();

  similarityAlignmentActions = createSimilarityAlignmentActions({
    state,
    getOrthogroupById: orthogroupActions.getOrthogroupById,
    getEnrichedOrthogroupMembers: orthogroupActions.getEnrichedOrthogroupMembers,
    selectedOrthogroupDetail: orthogroupActions.selectedOrthogroupDetail,
    getRecordCatalog: getAnnotationRecordCatalog,
    getCommittedRequest: getCommittedCanonicalRenderRequest,
    getCommittedSession: getCommittedCanonicalSession,
    projectCommittedAlignment: /** @type {any} */ (projectCommittedSimilarityAlignment),
    // R13: the root runs the candidate with the record display's orientation
    // checkpoint; the alignment owner holds neither owner.
    runRecordAlignment: ({ orientations, ...run }) => runCommittedCanonicalCandidate({ ...run,
      captureIntentCheckpoint: () => recordDisplayControls.captureAlignmentOrientationIntent(orientations),
      restoreIntentCheckpoint: (/** @type {AlignmentOrientationCheckpoint[]} */ checkpoint) => recordDisplayControls.restoreAlignmentOrientationIntent(checkpoint),
      commitIntent: () => recordDisplayControls.commitAlignmentOrientations(orientations) }),
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
  // Generate asks about Upload pairs without a file: each pair's row
  // takes its own file; Skip sets the pairs still without one to No comparison
  // in one History step, then Generates.
  /** @type {{ show: boolean, pairs: ReturnType<typeof missingUploadPairsToResolve> }} */
  const missingBlastTsvDialog = reactive({ show: false, pairs: [] });
  const closeMissingBlastTsvDialog = () => { missingBlastTsvDialog.show = false; };
  const setMissingBlastTsvPairsToNone = dialogChoice.withHistory(
    () => msg('Set pairs to No comparison'),
    () => {
      // The pairs as they are now: Undo may have changed them while the dialog was open.
      const drawing = state.drawings.linear;
      const pairs = missingUploadPairsToResolve(drawing.linearComparisonResolution.value);
      if (pairs.length) omitLinearComparisonPairs(drawing, pairs.map((entry) => entry.edgeKey));
      dialogChoice.closeAfterChoice(closeMissingBlastTsvDialog);
      return true;
    },
    closeMissingBlastTsvDialog
  );
  // The dialog's rows: the pairs it opened with, each with its file as it is now.
  const missingBlastTsvRows = computed(() => {
    const missing = new Set(/** @type {readonly LinearComparisonIssue[]} */ (state.drawings.linear.linearComparisonResolution.value?.errors || [])
      .filter((issue) => issue.code === 'missing-upload').map((issue) => issue.edgeKey));
    return missingBlastTsvDialog.pairs.map((pair) => {
      const found = linearComparisonPairForEdgeKey(pair.edgeKey);
      return { ...pair, queryUid: found?.queryUid || '', subjectUid: found?.subjectUid || '',
        fileName: found?.draft?.fileActive ? String(found.draft.file?.name || '') : '', missing: missing.has(pair.edgeKey) };
    });
  });
  // Labels that do not fit (LP-D43, LP-D48): the dialog over a Generate that
  // drew nothing, or over crossing leaders from the notice bar.
  const { labelFitDialog, openLabelFitDialog, resolveLabelFitChoice } = createLabelFitDialog({
    reactive,
    computeRow: async (rowToken, kind) => (await runDiagramHelperOperation(
      DIAGRAM_HELPER_OPERATIONS.COMPUTE_LABEL_FIT_ROW, { rowToken, kind }
    )).result,
    writeRowSettings: (values, rowMode) => applyLabelFitValues(state.drawings[rowMode], values, rowMode),
    withHistory: dialogChoice.withHistory,
    closeAfterChoice: dialogChoice.closeAfterChoice,
    choicePending: () => dialogChoice.pending.value,
    runAnalysis: () => runAnalysis(),
    errorSummary: (error) => normalizeCaughtError(error, { operation: 'computeLabelFitRow', stage: 'helper' }).summary
  });
  // The notice bar over the shown Result drawn with crossing leaders or labels
  // over the definition text (WARN-NOTICE, DEF-WARN).
  const labelFitNotice = computed(() => {
    const entry = /** @type {readonly ResultLabelFit[]} */ (labelFitReports.value).find((item) => item.resultIndex === selectedResultIndex.value
      && item.resultName === results.value[selectedResultIndex.value]?.name);
    return entry && entry.report.drawn ? entry : null;
  });
  const openLabelFitNotice = () => {
    const entry = labelFitNotice.value;
    if (!entry) return undefined;
    return openLabelFitDialog(entry.report, { rowToken: entry.rowToken, mode: mode.value === 'linear' ? 'linear' : 'circular' });
  };
  /** @param {'generate' | 'none' | 'cancel'} choice */
  const resolveMissingBlastTsv = async (choice) => {
    if (choice === 'generate') {
      if (dialogChoice.pending.value) return undefined;
      closeMissingBlastTsvDialog();
      return runAnalysis();
    }
    return await setMissingBlastTsvPairsToNone(choice) === true ? runAnalysis() : undefined;
  };
  // The Review alignment panel, docked in the Result workspace. The figure overlay and the
  // panel's strips show one draft: the root passes the draft to the preview's owner, which
  // returns the figure geometry the strips are laid out from.
  const similarityAlignmentCanvasHover = ref(/** @type {{ recordKey: string, candidateKey: string } | null} */ (null));
  const similarityAlignmentGeometry = shallowRef(/** @type {AlignmentFigureGeometry | null} */ (null));
  // The takeover ends with the review (Cancel, the close button, Escape, Apply). The pan/zoom
  // owner holds the view from before the fit, so a Session saved or a Result checkpointed
  // meanwhile keeps that view. While Apply runs, the ending waits for its outcome, and its new
  // Result is the review's own.
  let similarityAlignmentApplying = false;
  // The identity of the selected Result of the review's mode, shown or kept in its slot.
  /** @param {'circular' | 'linear'} reviewMode */
  const reviewModeFigure = (reviewMode) => previewRuntime.getResultIdentity(readArtifactOfMode(reviewMode,
    mode.value, artifactSlots, () => results.value[Number(selectedResultIndex.value) || 0],
    (slot) => slot.values.results[slot.values.selectedResultIndex]));
  // The draft no longer matches the figure. Another Result of the review's mode replaced it (a
  // Session load, History, a new Generate), whichever mode is shown: the view and the drawer are
  // that change's, and the review gives neither back. A mode switch brings no view: the review's
  // own Result, kept in its slot or shown again, leaves the hold, so the end gives the view back.
  const dropSimilarityAlignmentView = () => {
    const takeover = similarityAlignmentTakeover.value;
    if (!takeover || similarityAlignmentApplying || reviewModeFigure(takeover.mode) === takeover.figure) return;
    takeover.drawerOpen = false;
    releasePreviewViewport({ restore: false });
  };
  similarityAlignmentPorts.refreshCanvas = () => {
    similarityAlignmentCanvasHover.value = null;
    const current = similarityAlignmentActions.draft.value;
    if (!similarityAlignmentActions.dialogOpen.value) {
      featureActions.clearAlignmentOverlay();
      similarityAlignmentGeometry.value = null;
      return;
    }
    if (similarityAlignmentActions.status.value !== 'reviewing') {
      // While Apply runs, the strips keep the figure they were measured from.
      featureActions.clearAlignmentOverlay();
      return;
    }
    if (!similarityAlignmentActions.isDraftArtifactCurrent()) {
      featureActions.clearAlignmentOverlay();
      similarityAlignmentGeometry.value = null;
      dropSimilarityAlignmentView();
      return;
    }
    // A live rerender of the draft's records draws a new Result of the same figure (PD-OI-102).
    if (similarityAlignmentTakeover.value) similarityAlignmentTakeover.value.figure = reviewModeFigure(similarityAlignmentTakeover.value.mode);
    // Without a mounted preview (the Run info tab) nothing is measured: the strips keep the
    // geometry of the same Result.
    similarityAlignmentGeometry.value = featureActions.showAlignmentOverlay({
      reference: { key: ALIGNMENT_REFERENCE_KEY, recordKey: current.reference.recordKey,
        anchor: current.response.reference, name: current.reference.name },
      records: alignmentOverlayRecords(current.rows),
      onSelect: similarityAlignmentActions.selectCandidate,
      onHover: (recordKey, candidateKey) => {
        similarityAlignmentCanvasHover.value = recordKey ? { recordKey, candidateKey } : null;
      }
    }) ?? similarityAlignmentGeometry.value;
  };
  watch([
    similarityAlignmentActions.dialogOpen,
    () => similarityAlignmentActions.draft.value?.rows,
    similarityAlignmentActions.status
  ], similarityAlignmentPorts.refreshCanvas, { flush: 'post' });
  // The preview is shown again after this flush's watchers have run: measure on the next tick.
  watch(resultPanelTab, () => nextTick(similarityAlignmentPorts.refreshCanvas));
  const similarityAlignmentReview = computed(() => {
    const draft = similarityAlignmentActions.draft.value;
    const directions = similarityAlignmentActions.directionPreview.value?.records;
    return draft && directions ? buildAlignmentReviewView({
      rows: draft.rows, directions, geometry: similarityAlignmentGeometry.value,
      openedWithoutChoice: draft.openedWithoutChoice
    }) : null;
  });
  /** @type {HTMLElement | null} */
  let similarityAlignmentReturnFocus = null;
  const similarityAlignmentEditorDisabledReason = computed(() => (
    similarityAlignmentActions.dialogOpen.value ? t('Finish or cancel alignment review before opening Editor.') : ''
  ));
  similarityAlignmentPorts.editorDisabledReason = () => similarityAlignmentEditorDisabledReason.value;
  // The figure is fitted to the room the docked panel leaves, once the panel is in the page. The
  // pan/zoom owner fits it again when that room changes, for as long as it holds the view.
  const fitPreviewForSimilarityAlignment = () => {
    if (similarityAlignmentActions.dialogOpen.value) fitPreviewToViewport();
  };
  watch(similarityAlignmentActions.dialogOpen, fitPreviewForSimilarityAlignment, { flush: 'post' });
  const endSimilarityAlignmentTakeover = () => {
    const takeover = similarityAlignmentTakeover.value;
    if (!takeover || similarityAlignmentApplying) return;
    // The held view is the review's mode's: ended while the other mode is shown, the review
    // does not write it onto that mode's figure (OV-1089).
    const restore = similarityAlignmentHoldShown();
    similarityAlignmentTakeover.value = null;
    if (takeover.drawerOpen) rightDrawerActions.openRightDrawerTab(takeover.drawerTab);
    // The drawer narrows the visible frame, so the view returns after it.
    releasePreviewViewport({ restore });
    // Open in the feature search row retargeted the hidden popup for the fitted view.
    if (clickedFeature.value && clickedFeature.value !== takeover.popup) void featureActions.showFeaturePopupAtFeature();
  };
  watch(similarityAlignmentActions.dialogOpen, (/** @type {boolean} */ open) => {
    if (!open) {
      endSimilarityAlignmentTakeover();
      return;
    }
    // A review shown again while an Apply is pending keeps what the first one took.
    const reviewMode = mode.value === 'circular' ? 'circular' : 'linear';
    similarityAlignmentTakeover.value ??= { drawerOpen: showRightDrawer.value, drawerTab: rightDrawerTab.value, popup: clickedFeature.value,
      mode: reviewMode, figure: reviewModeFigure(reviewMode) };
    // The figure is fitted again only while it shows the draft.
    holdPreviewViewport(similarityAlignmentActions.isDraftArtifactCurrent);
    rightDrawerActions.closeRightDrawer();
    settingsPaneActions.closeSettingsPopup();
  }, { flush: 'sync' });
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
  /** @param {KeyboardEvent} event */
  const handleSimilarityAlignmentEscape = (event) => {
    if (event.key !== 'Escape' || !similarityAlignmentActions.dialogOpen.value) return;
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
    similarityAlignmentActions.openReset();
  };
  watch([similarityAlignmentActions.dialogOpen, similarityAlignmentActions.resetDialogOpen], async (/** @type {[unknown, unknown]} */ [open, resetOpen]) => {
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
    if (similarityAlignmentActions.dialogOpen.value) {
      if (similarityAlignmentActions.error.value) {
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
      await similarityAlignmentActions.startFromDrawer({ groupId, mode })
    );
  };
  const cancelSimilarityAlignmentDialog = async () => {
    similarityAlignmentActions.cancel();
  };
  const applySimilarityAlignmentDialog = async (reset = false) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!reset) similarityAlignmentApplying = true;
    let outcome;
    try {
      outcome = await (reset ? similarityAlignmentActions.applyReset() : similarityAlignmentActions.applyDraft());
    } finally {
      if (!reset) similarityAlignmentApplying = false;
    }
    if (!reset) {
      if (!similarityAlignmentActions.dialogOpen.value) {
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
  const canUseClickedOrthogroupActions = computed(() => {
    const drawing = state.activeDrawing();
    const cf = clickedFeature.value;
    return Boolean(
      cf &&
      mode.value === 'linear' &&
      drawing.hasActiveLinearLosatIntent.value &&
      drawing.linearComparisonResolution.value.valid &&
      drawing.losatProgram.value === 'blastp' &&
      drawing.losat.blastp?.mode === 'orthogroup' &&
      lInputType.value === 'gb' &&
      cf.feat?.type === 'CDS' &&
      cf.orthogroupId
    );
  });

  const clickedOrthogroupDetail = computed(() => {
    const cf = clickedFeature.value;
    const orthogroupId = String(cf?.orthogroupId || '').trim();
    if (!orthogroupId) return null;
    const group = orthogroupActions.getOrthogroupById(orthogroupId);
    if (!group) return null;
    const members = orthogroupActions.getEnrichedOrthogroupMembers(group);
    const currentMember = resolveUniqueOrthogroupMemberForFeature(cf?.feat, members);
    const membersByRecord = orthogroupActions.groupOrthogroupMembersByRecord(members);
    const memberCount = Number(group.member_count || members.length || 0);
    const recordCoverage = Number(group.record_coverage_count || membersByRecord.length || 0);
    return {
      id: orthogroupId,
      displayName: orthogroupActions.resolveOrthogroupName(group),
      memberCount,
      recordCoverage,
      countText: groupCountText(memberCount, recordCoverage),
      genesText: groupCountText(memberCount, 0),
      ntSequenceCount: orthogroupActions.getOrthogroupSequenceCount(group, 'nt'),
      aaSequenceCount: orthogroupActions.getOrthogroupSequenceCount(group, 'aa'),
      currentMember,
      membersByRecord,
      scopeLabel: orthogroupActions.orthogroupScope(group) === 'global'
        ? ''
        : orthogroupActions.orthogroupScopeLabel(group)
    };
  });

  const alignByClickedOrthogroup = async (event = null, mode = 'align') => {
    const detail = clickedOrthogroupDetail.value;
    if (!detail?.id) return { status: 'rejected' };
    rememberSimilarityAlignmentInvoker(event);
    const outcome = await similarityAlignmentActions.startFromPopup({
      groupId: detail.id,
      reference: detail.currentMember,
      mode
    });
    return finishSimilarityAlignmentStart(outcome);
  };

  const highlightClickedOrthogroup = () => {
    const cf = clickedFeature.value;
    const orthogroupId = String(cf?.orthogroupId || '').trim();
    if (!orthogroupId) return;
    orthogroupActions.highlightOrthogroupById(orthogroupId);
  };

  const clearOrthogroupHighlight = () => {
    orthogroupActions.clearOrthogroupHighlight();
  };

  /** @param {string} orthogroupId */
  const openOrthogroupInDrawer = (orthogroupId) => {
    if (!orthogroupActions.selectOrthogroup(orthogroupId)) return false;
    rightDrawerActions.openRightDrawerTab('orthogroups');
    return true;
  };

  const reselectSimilarityAlignmentReference = () => {
    const groupId = similarityAlignmentActions.repair.value?.groupId
      || state.similarityAlignmentPlan.value?.groupId
      || '';
    similarityAlignmentActions.setDrawerReference(groupId, '');
    return openOrthogroupInDrawer(groupId);
  };

  const openClickedOrthogroupInEditor = () => {
    const orthogroupId = String(clickedFeature.value?.orthogroupId || '').trim();
    if (!openOrthogroupInDrawer(orthogroupId)) return;
    clickedFeature.value = null;
  };

  const { resetAllPositions, resetCanvasPadding } = legendLayout;

  // Reset Settings resets both drawings (PD-OI-070); the Linear records and
  // comparisons are the Linear drawing's. The button asks first and says what
  // the reset keeps; `resetSettings` itself does not ask.
  const requestResetSettings = async () => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const choice = await askChoice({
      title: t('Reset all settings?'),
      description: t('Circular and Linear both return to their default settings.') + ' ' + t('This also removes your edits:'),
      items: [
        t('Colors, color rules, and the palette'),
        t('Labels, Feature visibility, and Feature placements'),
        t('Legend names, deleted Legend rows, Similarity group names, and Region Annotations'),
        t('Linear record definitions, subtitles, regions, reverse complement, and alignment'),
        t('Linear comparisons, which turn off; their uploaded tables stay'),
        t('The Default colors, Specific rules, Blacklist, Whitelist, and Priority files, with the settings they loaded')
      ],
      note: t('Input files (genome, Depth, and Circular conservation files), the records you chose, and File defaults stay.') + ' ' + t('Undo brings the settings back.'),
      actions: [{ id: 'reset', label: t('Reset all settings'), tone: 'danger' }],
      cancelLabel: t('Keep settings')
    });
    if (choice !== 'reset') return false;
    return resetSettings();
  };
  const resetSettings = () => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
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
    return holdAutomaticRerender(() => history.runUndoableCheckpoint(msg('Reset settings'), () => followRuleEdits(state.activeDrawing(), async () => {
      featureActions.clearSpecificRulePatternDrafts();
      const shown = Boolean(svgContainer.value?.querySelector?.('svg'));
      const palette = { ...appliedPaletteColors.value };
      // The deleted rows Generate drew return (a row added in the editor and
      // deleted stays gone).
      const returning = state.activeDrawing().deletedLegendEntries.value
        .map((/** @type {{ caption?: string, originalCaption?: string }} */ entry) => String(entry.originalCaption || entry.caption || ''))
        .filter((/** @type {string} */ key) => originalLegendOrder.value.includes(key));
      resetSettingsState(state);
      // OV-1249: the loaded Circular records keep a Record Order row each, at
      // their defaults; Reset resets settings, not the records.
      state.drawings.circular.adv.multi_record_positions.push(...buildDefaultCircularRecordPositions());
      // Linear records return to their File defaults and inferred definitions;
      // Files, record selections, File defaults, and depth stay. The mutation
      // also clears the alignment plan (D-15).
      applyLinearSeqMutation(drawing, linearSeqs.map((sequence) => ({
        ...sequence,
        definition: '',
        record_subtitle: '',
        region_start: null,
        region_end: null,
        region_reverse: false
      })), { alignmentMutation: msg('settings reset.') });
      invalidateLinearComparisonArtifacts(drawing);
      matchSequenceRegistry?.reset?.();
      circularTrackNewRenderer.value = 'dinucleotide_skew';
      linearTrackNewRenderer.value = 'dinucleotide_skew';
      depthTrackUiCounts.circular = 1;
      ensureDepthTrackConfigCount(state.drawings.circular, activeDepthTrackCount('circular'));
      ensureDepthTrackConfigCount(drawing, activeDepthTrackCount('linear'));
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
  const layoutResetAvailable = computed(() => {
    const moved = (/** @type {{ x: number, y: number }} */ offset) => Boolean(Number(offset.x) || Number(offset.y));
    const padding = state.activeDrawing().canvasPadding;
    return [state.diagramOffset, state.legendCurrentOffset, state.plotTitleUserOffset, state.lengthBarUserOffset, canvasPan].some(moved)
      || zoom.value !== 1
      || ['top', 'right', 'bottom', 'left'].some((side) => Number(padding[side]));
  });

  /** @param {Element | null} target */
  const isInteractiveTarget = (target) => {
    if (!target) return false;
    return Boolean(target.closest('input, textarea, select, button, label, a, [data-nodrag="true"]'));
  };

  const FEATURE_POPUP_MARGIN = 12;
  const FEATURE_POPUP_RICH_MIN_WIDTH = 360;
  const FEATURE_POPUP_SIMPLE_MIN_WIDTH = 300;
  const FEATURE_POPUP_MIN_HEIGHT = 220;
  const PAIRWISE_MATCH_POPUP_MIN_WIDTH = 280;
  const PAIRWISE_MATCH_POPUP_MIN_HEIGHT = 180;

  /** @param {unknown} value @param {number} min @param {number} max */
  const clampNumber = (value, min, max) => {
    const safeMin = Number.isFinite(min) ? min : 0;
    const safeMax = Number.isFinite(max) ? Math.max(safeMin, max) : safeMin;
    const numeric = Number(value);
    if (!Number.isFinite(numeric)) return safeMin;
    return Math.min(Math.max(numeric, safeMin), safeMax);
  };

  // A popup's size limits at its position, inside the popup bounds.
  /**
   * @param {number} left
   * @param {number} top
   * @param {number} desiredMinWidth
   * @param {number} desiredMinHeight
   */
  const popupConstraints = (left, top, desiredMinWidth, desiredMinHeight) => {
    const bounds = popupViewportBounds();
    const availableWidth = Math.max(1, bounds.right - (FEATURE_POPUP_MARGIN * 2));
    const availableHeight = Math.max(1, bounds.bottom - bounds.top - (FEATURE_POPUP_MARGIN * 2));
    const minWidth = Math.min(desiredMinWidth, availableWidth);
    const minHeight = Math.min(desiredMinHeight, availableHeight);
    return {
      minWidth,
      minHeight,
      maxWidth: Math.max(minWidth, bounds.right - left - FEATURE_POPUP_MARGIN),
      maxHeight: Math.max(minHeight, bounds.bottom - top - FEATURE_POPUP_MARGIN)
    };
  };
  // The popup's size follows the app-level rich popup preference.
  const getFeaturePopupConstraints = (left = clickedFeaturePos.x, top = clickedFeaturePos.y) => popupConstraints(
    left, top,
    richFeaturePopup.value === false ? FEATURE_POPUP_SIMPLE_MIN_WIDTH : FEATURE_POPUP_RICH_MIN_WIDTH,
    FEATURE_POPUP_MIN_HEIGHT
  );
  const getPairwiseMatchPopupConstraints = (left = clickedPairwiseMatchPos.x, top = clickedPairwiseMatchPos.y) => popupConstraints(
    left, top, PAIRWISE_MATCH_POPUP_MIN_WIDTH, PAIRWISE_MATCH_POPUP_MIN_HEIGHT
  );

  /**
   * @param {PopupPlacement} pos
   * @param {{ width: number, height: number }} size
   * @param {{ maxHeight: number }} constraints
   */
  const popupStyle = (pos, size, constraints) => ({
    top: `${pos.y}px`,
    left: `${pos.x}px`,
    maxHeight: `${pos.maxHeight > 0 ? Math.min(pos.maxHeight, constraints.maxHeight) : constraints.maxHeight}px`,
    ...(size.width > 0 ? { width: `${size.width}px` } : {}),
    ...(size.height > 0 ? { height: `${size.height}px` } : {})
  });
  const featurePopupStyle = computed(() => popupStyle(clickedFeaturePos, featurePopupSize, getFeaturePopupConstraints()));
  const pairwiseMatchPopupStyle = computed(() => popupStyle(clickedPairwiseMatchPos, pairwiseMatchPopupSize, getPairwiseMatchPopupConstraints()));
  // The settings popup has the width of its stylesheet and the height of its section.
  const settingsPopupStyle = computed(() => {
    const pos = state.settingsPopup.value;
    return pos ? popupStyle(pos, { width: 0, height: 0 }, popupConstraints(pos.x, pos.y, 0, 0)) : null;
  });

  // OV-371: a popup opened at a clicked point was placed with an expected
  // size; once rendered it is placed again with its own size and content height.
  /**
   * @param {PopupPosition} pos
   * @param {HTMLElement | null} element
   * @param {string} bodySelector
   */
  const placeRenderedPopup = (pos, element, bodySelector) => {
    const point = pos.anchor;
    pos.anchor = null;
    if (!point || !element) return;
    const body = element.querySelector(bodySelector);
    const height = element.offsetHeight + (body ? body.scrollHeight - body.clientHeight : 0);
    Object.assign(pos, placePopupNearPoint({ point, width: element.offsetWidth, height, bounds: popupViewportBounds() }));
  };
  watch(clickedFeature, () => placeRenderedPopup(clickedFeaturePos, featurePopupRef.value, '.feature-popup-body'), { flush: 'post' });
  watch(clickedPairwiseMatch, () => placeRenderedPopup(clickedPairwiseMatchPos, pairwiseMatchPopupRef.value, '.pairwise-match-popup-body'), { flush: 'post' });

  // A press on the handle of an open popup drags it; every popup of the Result moves this way.
  // A drag keeps the popup inside the popup bounds (OV-384: below the header), and a window that
  // changed size moves an open popup the least distance back inside them (OV-626).
  /**
   * @param {object} popup
   * @param {{ active: boolean, offsetX: number, offsetY: number }} popup.drag
   * @param {() => boolean} popup.isOpen
   * @param {(placement: PopupPlacement) => void} popup.moveTo The popup owner's move to a placement.
   * @param {() => HTMLElement | null} popup.element
   * @param {(target: Element) => boolean} popup.isHandle Whether a press on `target` drags the popup.
   */
  const createPopupDrag = ({ drag, isOpen, moveTo, element, isHandle }) => {
    /** @param {MouseEvent} event */
    const move = (event) => {
      const popup = element();
      if (!drag.active || !isOpen()) return;
      moveTo(placePopupWithin({
        x: event.clientX - drag.offsetX, y: event.clientY - drag.offsetY,
        width: popup?.offsetWidth || 0, height: popup?.offsetHeight || 0, bounds: popupViewportBounds()
      }));
    };
    const end = () => {
      drag.active = false;
      document.removeEventListener('mousemove', move);
      document.removeEventListener('mouseup', end);
    };
    /** @param {MouseEvent & { target: Element }} event A press on the popup; its target is an element. */
    const start = (event) => {
      const popup = element();
      if (event.button !== 0 || !popup || !isOpen() || !isHandle(event.target)) return;
      const rect = popup.getBoundingClientRect();
      drag.active = true;
      drag.offsetX = event.clientX - rect.left;
      drag.offsetY = event.clientY - rect.top;
      document.addEventListener('mousemove', move);
      document.addEventListener('mouseup', end);
      event.preventDefault();
    };
    const keepInWindow = () => {
      const popup = element();
      if (!popup || !isOpen()) return;
      const { left: x, top: y } = popup.getBoundingClientRect();
      const within = placePopupWithin({ x, y, width: popup.offsetWidth, height: popup.offsetHeight, bounds: popupViewportBounds() });
      if (within.x !== x || within.y !== y) moveTo(within);
    };
    return { start, end, keepInWindow };
  };
  const featurePopupDragging = createPopupDrag({
    drag: featurePopupDrag,
    isOpen: () => Boolean(clickedFeature.value),
    moveTo: (placement) => { Object.assign(clickedFeaturePos, placement); },
    element: () => featurePopupRef.value,
    isHandle: (target) => !featurePopupResize.active && !isInteractiveTarget(target)
  });
  const pairwiseMatchPopupDragging = createPopupDrag({
    drag: pairwiseMatchPopupDrag,
    isOpen: () => Boolean(clickedPairwiseMatch.value),
    moveTo: (placement) => { Object.assign(clickedPairwiseMatchPos, placement); },
    element: () => pairwiseMatchPopupRef.value,
    isHandle: (target) => !pairwiseMatchPopupResize.active && !isInteractiveTarget(target)
  });
  // The heading of the section in the settings popup is the popup's handle. It is the section's
  // own disclosure heading, which does not fold the section there.
  const SETTINGS_POPUP_HEADING = '#settings-popup > [data-settings-section] > details > summary';
  /** @param {Element} target */
  const isSettingsPopupHeading = (target) => !isInteractiveTarget(target) && Boolean(target.closest(SETTINGS_POPUP_HEADING));
  const settingsPopupDrag = createPopupDrag({
    drag: { active: false, offsetX: 0, offsetY: 0 },
    isOpen: () => Boolean(state.settingsPopup.value),
    moveTo: settingsPaneActions.moveSettingsPopup,
    element: () => document.getElementById('settings-popup'),
    isHandle: isSettingsPopupHeading
  });
  const keepPopupsInWindow = () => [featurePopupDragging, pairwiseMatchPopupDragging, settingsPopupDrag].forEach((popup) => popup.keepInWindow());
  onMounted(() => window.addEventListener('resize', keepPopupsInWindow));
  onUnmounted(() => window.removeEventListener('resize', keepPopupsInWindow));
  /** @param {MouseEvent & { target: Element }} event A click in the popup; its target is an element. */
  const keepSettingsPopupSectionOpen = (event) => {
    if (isSettingsPopupHeading(event.target)) event.preventDefault();
  };

  const pairwiseBlockOrthogroups = computed(() => (
    Array.isArray(clickedPairwiseMatch.value?.blockOrthogroups)
      ? clickedPairwiseMatch.value.blockOrthogroups
      : []
  ));
  const selectedPairwiseBlockOrthogroup = computed(() => {
    const selectedId = String(selectedPairwiseBlockOrthogroupId.value || '').trim();
    if (!selectedId) return null;
    return pairwiseBlockOrthogroups.value.find((/** @type {{ id?: unknown }} */ group) => String(group?.id || '').trim() === selectedId) || null;
  });
  // The group the match popup lists members of: a similarity-group ribbon's
  // own group, or the group picked in a collinear block's gene-pair table.
  const pairwiseGroupView = computed(() => (
    clickedPairwiseMatch.value?.block
      ? selectedPairwiseBlockOrthogroup.value
      : clickedPairwiseMatch.value?.group || null
  ));
  /** @param {{ id?: unknown } | null | undefined} group */
  const selectPairwiseBlockOrthogroup = (group) => {
    const id = String(group?.id || '').trim();
    selectedPairwiseBlockOrthogroupId.value = selectedPairwiseBlockOrthogroupId.value === id ? '' : id;
  };
  /** @param {FeatureListRow} row @param {MouseEvent} event */
  const openPairwiseFeatureRow = (row, event) => {
    if (!row?.canOpen || !row?.feature) return;
    openFeatureEditorForFeature(row.feature, event);
  };

  // The bands of the open match popup carry `gbdraw-match-selected`; the SVG
  // owner marks them from its band indexes (OV-937).
  watch(clickedPairwiseMatch, (/** @type {{ id?: unknown, matchKind?: unknown, orthogroupId?: unknown } | null} */ match) => {
    selectedPairwiseBlockOrthogroupId.value = '';
    featureActions.markSelectedPairwiseMatch(match);
  });

  /** @param {MouseEvent} event */
  const onPairwiseMatchPopupResize = (event) => {
    if (!pairwiseMatchPopupResize.active) return;
    const constraints = getPairwiseMatchPopupConstraints(clickedPairwiseMatchPos.x, clickedPairwiseMatchPos.y);
    const nextWidth = pairwiseMatchPopupResize.startWidth + (event.clientX - pairwiseMatchPopupResize.startX);
    const nextHeight = pairwiseMatchPopupResize.startHeight + (event.clientY - pairwiseMatchPopupResize.startY);
    pairwiseMatchPopupSize.width = clampNumber(nextWidth, constraints.minWidth, constraints.maxWidth);
    pairwiseMatchPopupSize.height = clampNumber(nextHeight, constraints.minHeight, constraints.maxHeight);
    event.preventDefault();
  };

  const endPairwiseMatchPopupResize = () => {
    if (!pairwiseMatchPopupResize.active) return;
    pairwiseMatchPopupResize.active = false;
    document.removeEventListener('mousemove', onPairwiseMatchPopupResize);
    document.removeEventListener('mouseup', endPairwiseMatchPopupResize);
  };

  /** @param {MouseEvent & { target: Element }} event A press on the popup; its target is an element. */
  const startPairwiseMatchPopupResize = (event) => {
    if (event.button !== 0) return;
    if (!clickedPairwiseMatch.value) return;
    const popup = pairwiseMatchPopupRef.value;
    if (!popup) return;
    const rect = popup.getBoundingClientRect();
    const constraints = getPairwiseMatchPopupConstraints(rect.left, rect.top);

    pairwiseMatchPopupDragging.end();

    clickedPairwiseMatchPos.maxHeight = 0;
    pairwiseMatchPopupResize.active = true;
    pairwiseMatchPopupResize.startX = event.clientX;
    pairwiseMatchPopupResize.startY = event.clientY;
    pairwiseMatchPopupResize.startWidth = clampNumber(rect.width, constraints.minWidth, constraints.maxWidth);
    pairwiseMatchPopupResize.startHeight = clampNumber(rect.height, constraints.minHeight, constraints.maxHeight);
    pairwiseMatchPopupSize.width = pairwiseMatchPopupResize.startWidth;
    pairwiseMatchPopupSize.height = pairwiseMatchPopupResize.startHeight;
    document.addEventListener('mousemove', onPairwiseMatchPopupResize);
    document.addEventListener('mouseup', endPairwiseMatchPopupResize);
    event.preventDefault();
    event.stopPropagation();
  };

  /** @param {MouseEvent} event */
  const onFeaturePopupResize = (event) => {
    if (!featurePopupResize.active) return;
    const constraints = getFeaturePopupConstraints(clickedFeaturePos.x, clickedFeaturePos.y);
    const nextWidth = featurePopupResize.startWidth + (event.clientX - featurePopupResize.startX);
    const nextHeight = featurePopupResize.startHeight + (event.clientY - featurePopupResize.startY);
    featurePopupSize.width = clampNumber(nextWidth, constraints.minWidth, constraints.maxWidth);
    featurePopupSize.height = clampNumber(nextHeight, constraints.minHeight, constraints.maxHeight);
    event.preventDefault();
  };

  const endFeaturePopupResize = () => {
    if (!featurePopupResize.active) return;
    featurePopupResize.active = false;
    document.removeEventListener('mousemove', onFeaturePopupResize);
    document.removeEventListener('mouseup', endFeaturePopupResize);
  };

  /** @param {MouseEvent & { target: Element }} event A press on the popup; its target is an element. */
  const startFeaturePopupResize = (event) => {
    if (event.button !== 0) return;
    if (!clickedFeature.value) return;
    const popup = featurePopupRef.value;
    if (!popup) return;
    const rect = popup.getBoundingClientRect();
    const constraints = getFeaturePopupConstraints(rect.left, rect.top);

    featurePopupDragging.end();

    clickedFeaturePos.maxHeight = 0;
    featurePopupResize.active = true;
    featurePopupResize.startX = event.clientX;
    featurePopupResize.startY = event.clientY;
    featurePopupResize.startWidth = clampNumber(rect.width, constraints.minWidth, constraints.maxWidth);
    featurePopupResize.startHeight = clampNumber(rect.height, constraints.minHeight, constraints.maxHeight);
    featurePopupSize.width = featurePopupResize.startWidth;
    featurePopupSize.height = featurePopupResize.startHeight;
    document.addEventListener('mousemove', onFeaturePopupResize);
    document.addEventListener('mouseup', endFeaturePopupResize);
    event.preventDefault();
    event.stopPropagation();
  };

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

  // The popup header names the gene (else the feature label), then
  // its product, then "<record ID> · <location>" with grouped digits (parts
  // joined by "; ", since a grouped number holds a comma), all from payload
  // fields only.
  const clickedFeatureHeading = computed(() => {
    const cf = clickedFeature.value;
    const title = String(cf?.gene || cf?.label || '');
    return { title, product: cf?.product && cf.product !== title ? cf.product : '' };
  });
  const clickedFeatureSummary = computed(() => {
    const cf = clickedFeature.value;
    const location = String(cf?.location || '').replaceAll(', ', '; ')
      .replace(/\d{4,}/g, (digits) => Number(digits).toLocaleString('en-US'));
    return [cf?.recordId, location].filter(Boolean).join(' · ');
  });

  // One Copy per tab: the header lines and the rows, or every qualifier.
  const copyClickedFeatureDetails = () => {
    const cf = clickedFeature.value;
    if (!cf) return;
    const { title, product } = clickedFeatureHeading.value;
    copyTextToClipboard([title, product, clickedFeatureSummary.value, ...cf.detailRows.map((/** @type {{ label: string, value: string }} */ row) => `${row.label}: ${row.value}`)]
      .filter(Boolean).join('\n'));
  };
  const copyClickedFeatureQualifiers = () => {
    const rows = clickedFeature.value?.qualifierRows || [];
    copyTextToClipboard(rows.flatMap((/** @type {{ key: string, values: string[] }} */ row) => row.values.map((value) => `/${row.key}="${value}"`)).join('\n'));
  };

  /** @param {unknown} filename @param {unknown} text @param {string} [type] */
  const downloadText = (filename, text, type = 'text/plain;charset=utf-8') => {
    const value = String(text ?? '');
    if (!value) return;
    downloadTextFile(String(filename || 'gbdraw.txt'), value, type);
  };

  let latestExportOperation = 0;
  const failedInteractiveSvgExport = ref(null);
  const canRetryInteractiveSvgExport = computed(() => Boolean(failedInteractiveSvgExport.value
    && errorLog.value === failedInteractiveSvgExport.value));
  /** @param {'downloadSVG' | 'downloadInteractiveSVG' | 'downloadPNG' | 'downloadPDF'} methodName @param {string} operation */
  const runExportAction = async (methodName, operation) => {
    const operationId = ++latestExportOperation;
    const previousError = errorLog.value;
    try {
      const snapshot = captureSvgExport(state, { interactive: methodName === 'downloadInteractiveSVG' });
      const exportService = await loadExportService();
      const exportMethod = exportService?.[methodName];
      if (typeof exportMethod !== 'function') {
        throw new Error('The export service did not provide the requested action.');
      }
      // Custom fonts in exports (D-05, D-13): the Worker holds the fonts of
      // this page, so only a page with a loaded font embeds. SVG and
      // Interactive SVG follow Embed font in exported SVG; a PNG always draws
      // with the fonts, and a PDF always registers the faces its text uses.
      /** @type {import('../services/export.js').ExportFontPorts} */
      const fontPorts = {
        embedFonts: state.loadedFonts.value.length > 0
          ? async (svg, kind) => (await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.EMBED_FONTS, { svg, kind })).result
          : null,
        embedInSvg: fontSettings.embedInExports.value,
        resolvePdfFonts: async (runs, known) => (
          await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.PDF_FONTS, { runs, known })
        ).result
      };
      /** @type {import('../services/export.js').ExportOutcome} */
      const result = await exportMethod(snapshot, fontPorts);
      if (operationId === latestExportOperation && errorLog.value === previousError && previousError?.operation?.startsWith('export-')) errorLog.value = null;
      if (result.fontWarnings.length > 0) {
        showNotice({ text: t('Exported {file}.', { file: result.filename }) + ' ' + t('Find it in your browser\'s downloads.'), items: result.fontWarnings });
      }
      return result;
    } catch (error) {
      const normalized = normalizeUserFacingError(error, { operation, stage: 'export-capture' });
      if (operationId !== latestExportOperation || errorLog.value !== previousError) return { status: 'stale' };
      errorLog.value = normalized;
      failedInteractiveSvgExport.value = methodName === 'downloadInteractiveSVG' ? normalized : null;
      return { status: 'error', error: normalized };
    }
  };

  const downloadSVG = () => runExportAction('downloadSVG', 'export-svg');
  const downloadInteractiveSVG = () => (
    runExportAction('downloadInteractiveSVG', 'export-svg')
  );
  const downloadPNG = () => runExportAction('downloadPNG', 'export-png');
  const downloadPDF = () => runExportAction('downloadPDF', 'export-pdf');

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
    return history.runUndoable(msg('Rename Session'), () => { sessionTitle.value = title; });
  };

  const saveSessionWithTitle = () => exportSession(null, {
    availability: sessionSaveLoadAvailability,
    askChoice,
    recordDisplayRows: recordDisplayControls.allRows,
    // Save writes every Result: the other mode's slot goes beside the shown one (E1).
    readOtherModeArtifact: () => artifactSlots[mode.value === 'linear' ? 'circular' : 'linear'],
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
        omittedAnnotationRecordKeys: draftRequest ? getAnnotationRecordCatalog().omittedRecordKeys || [] : []
      };
    },
    onError: (error) => { errorLog.value = normalizeUserFacingError(error); },
    // UJ-09: a saved download is the History position Load Session compares with.
    // The notice names the file the download started.
    onSaved: (/** @type {{ filename: string }} */ { filename }) => {
      history.markSavePoint();
      // D-02: a font gbdraw does not bundle is saved by name only (receipt 〔警告〕).
      showNotice({
        text: t('Session saved: {file}.', { file: filename }) + ' ' + t('Find it in your browser\'s downloads.'),
        // Each drawing's font, and the font of the shown Result, whose committed
        // request the Session saves too (OV-742).
        items: customFontSaveWarnings(
          [
            ...SLICE_MODES.map((drawingMode) => state.drawings[drawingMode].adv.font_family),
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

  /** @param {unknown} selector */
  const getCircularRecordOrderLabel = (selector) => {
    const normalized = String(selector || '').trim();
    if (!normalized) return '';
    const records = Array.isArray(circularRecordList.value) ? circularRecordList.value : [];
    const matched = records.find((/** @type {{ selector?: unknown }} */ entry) => String(entry?.selector || '').trim() === normalized);
    if (!matched) return normalized;
    return `${normalized} (${String(matched.record_id || '').trim() || t('Unknown')})`;
  };

  const circularRecordDiscoveryState = computed(getCircularRecordDiscoveryState);
  const circularRecordInspectionEnabled = computed(() => semanticMutationAvailable.value
    && mode.value === 'circular'
    && circularRecordDiscoveryState.value.hasInput
    && circularRecordDiscoveryState.value.status !== 'loading'
    && !state.semanticFileWatchersSuppressed.value);
  const inspectCircularSourceRecords = () => sessionOperationAvailability() || (circularRecordInspectionEnabled.value
    ? refreshCircularRecordOrder()
    : Promise.resolve({ status: 'unavailable', reason: msg('Source inspection is unavailable.') }));
  const circularRecordPresentationEntries = () => buildDisambiguatedRecordEntries(
    circularRecordDiscoveryState.value.records.map(
      (/** @type {SequenceRecordEntry} */ record) => ({
        ...record,
        recordId: record?.record_id ?? record?.recordId,
        recordLength: record?.record_length ?? record?.recordLength
      })
    )
  );

  const circularDrawnRecordCount = computed(() => (
    circularRecordsToDraw(circularRecordPresentationEntries(), state.drawings.circular.recordsOff).length
  ));
  const circularRecordPresentationOptions = computed(() => {
    const drawing = state.activeDrawing();
    const entries = circularRecordPresentationEntries();
    const current = String(drawing.form.circular_record_selector || '').trim();
    const selection = resolveDisambiguatedRecordSelection(entries, current);
    // Without an inspected catalog the saved selector is unverified, not missing.
    const inspected = circularRecordDiscoveryState.value.status === 'ready';
    const drawnCount = circularDrawnRecordCount.value;
    // OV-379: a source whose records are not read says so, as the Linear select does.
    const automaticLabel = circularRecordDiscoveryState.value.status === 'deferred'
      ? tm(RECORDS_NOT_INSPECTED_LABEL)
      : entries.length > 1 || drawing.adv.circular_grouping_intent === 'batch'
        ? (drawnCount < entries.length ? t('{drawn} of {count} records (separate diagrams)', { drawn: drawnCount, count: entries.length }) : t('All records (separate diagrams)'))
        : t('Automatic (only record)');
    return [
      { value: '', label: automaticLabel, synthetic: false },
      ...(current && selection.status !== 'resolved'
        ? [{
            value: current,
            label: !inspected ? t('{selector} (not inspected)', { selector: current })
              : selection.status === 'ambiguous' ? t('{selector} (ambiguous)', { selector: current })
                : t('{selector} (not found)', { selector: current }),
            synthetic: true
          }]
        : []),
      ...entries.map((record) => ({
        value: record.value,
        label: `${record.recordId} (${formatRecordLength(record.recordLength)})${record.usesIndex ? ` [${record.selector}]` : ''}`,
        synthetic: false
      }))
    ];
  });

  const circularRecordPresentationError = computed(() => {
    const drawing = state.activeDrawing();
    const current = String(drawing.form.circular_record_selector || '').trim();
    if (!current || circularRecordDiscoveryState.value.status !== 'ready') return '';
    const selection = resolveDisambiguatedRecordSelection(
      circularRecordPresentationEntries(),
      current
    );
    if (selection.status === 'ambiguous') {
      return t('Record selector \'{selector}\' is ambiguous in the current input.', { selector: current });
    }
    if (selection.status === 'missing') {
      return t('Record selector \'{selector}\' was not found in the current input.', { selector: current });
    }
    return '';
  });

  // Why a Circular record's own fields (crop, orientation, title lines) are
  // unavailable, in the order the sidebar explains it; '' when they are not.
  const circularOwnRecordBlock = computed(() => {
    const drawing = state.activeDrawing();
    if (drawing.form.multi_record_canvas) return 'multi-record-canvas';
    if (circularRecordDiscoveryState.value.status !== 'ready') return 'not-inspected';
    const entries = circularRecordPresentationEntries();
    const selection = resolveDisambiguatedRecordSelection(entries, drawing.form.circular_record_selector);
    return drawing.adv.circular_grouping_intent !== 'batch'
      && (selection.status === 'resolved' || (selection.status === 'unspecified' && entries.length === 1))
      ? '' : 'no-single-record';
  });
  const circularSingleRecordPresentationEnabled = computed(() => mode.value === 'circular' && !circularOwnRecordBlock.value);

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

  const circularRecordSelectionEnabled = computed(() => semanticMutationAvailable.value
    && mode.value === 'circular'
    && !state.activeDrawing().form.multi_record_canvas && circularRecordDiscoveryState.value.status === 'ready');
  // CI-03: the per-record presentation a source file carries in each mode.
  // Replacing or removing the file resets it in the replacement's History step.
  const RETIRED_LINEAR_RECORD_PRESENTATION = Object.freeze({ region_record_id: '', region_start: null,
    region_end: null, region_reverse: false, definition: '', record_subtitle: '' });
  // The record reference names a record of the old file, so any replacement retires it.
  const RETIRED_LINEAR_RECORD_REFERENCE = Object.freeze({ region_record_id: '' });
  // The OFF records name records of the old file too (record selection).
  /** @param {DrawingState} drawing */
  const retireCircularRecordSelector = ({ form, adv, recordsOff }) => {
    form.circular_record_selector = '';
    if (adv.circular_grouping_intent === 'single') adv.circular_grouping_intent = 'auto';
    if (recordsOff.length) recordsOff.splice(0);
  };
  /** @param {DrawingState} drawing */
  const retireCircularRecordPresentation = (drawing) => {
    const { form } = drawing;
    retireCircularRecordSelector(drawing);
    form.circular_region_start = null;
    form.circular_region_end = null;
    form.circular_reverse = false;
    form.circular_record_label = '';
    form.circular_record_subtitle = '';
  };
  // D-32: one record replacing one record (a new version of the same genome)
  // keeps its crop and titles; a removed file or a replaced multi-record file
  // retires them. The record selector names a record of the old file, so any
  // replacement retires it.
  /** @param {{ removed: boolean, previousRecordCount: number }} replacement */
  const sourceReplacementRetiresPresentation = ({ removed, previousRecordCount }) => removed || previousRecordCount > 1;
  /**
   * The Circular source file controls (upload, replace, Remove).
   * @param {'c_gb' | 'c_gff' | 'c_fasta'} field
   * @param {File | null} value
   */
  const setCircularSourceFile = (field, value) => {
    const drawing = state.drawings.circular;
    const busy = sessionOperationAvailability();
    if (busy) return busy;
    const nextValue = value ?? null;
    const previous = files[field];
    if (previous === nextValue) return;
    files[field] = nextValue;
    if (nextValue) pendingCircularUpload = nextValue;
    // An upload into an empty slot starts reading another input (after an
    // Input switch), whose records the OFF list does not name.
    if (!previous) {
      if (drawing.recordsOff.length) drawing.recordsOff.splice(0);
      return;
    }
    if (sourceReplacementRetiresPresentation({
      removed: !nextValue, previousRecordCount: circularRecordList.value.length
    })) retireCircularRecordPresentation(drawing);
    else retireCircularRecordSelector(drawing);
  };
  /** @param {unknown} value */
  const setCircularRecordPresentationSelector = (value) => {
    const drawing = state.drawings.circular;
    const busy = sessionOperationAvailability();
    if (busy) return busy;
    if (!circularRecordSelectionEnabled.value) return { status: 'unavailable' };
    const normalized = String(value || '').trim();
    drawing.form.circular_record_selector = normalized;
    drawing.adv.circular_grouping_intent = normalized ? 'single' : 'auto';
    return { status: 'ok' };
  };
  const showCircularCanvasSetting = async () => {
    const drawing = state.drawings.circular;
    const control = /** @type {HTMLElement | null} */ (document.querySelector('[data-circular-canvas-setting]'));
    if (!control || mode.value !== 'circular' || !drawing.form.multi_record_canvas) {
      return { status: 'unavailable' };
    }
    await settingsPaneActions.revealSetting(control);
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

  const buildDefaultCircularRecordPositions = () => buildDefaultMultiRecordPositions(
    (Array.isArray(circularRecordList.value) ? circularRecordList.value : [])
      .map((/** @type {{ selector?: unknown }} */ entry) => entry?.selector)
  );

  /** @param {MultiRecordPosition | null | undefined} position */
  const getCircularRecordRow = (position) => {
    const drawing = state.drawings.circular;
    const rowValue = Number(position?.row);
    const maxRow = Math.max(1, Array.isArray(drawing.adv.multi_record_positions) ? drawing.adv.multi_record_positions.length : 1);
    if (!Number.isInteger(rowValue) || rowValue <= 0) return 1;
    return Math.min(rowValue, maxRow);
  };

  const getCircularRecordRowOptions = () => {
    const drawing = state.drawings.circular;
    const count = Array.isArray(drawing.adv.multi_record_positions) ? drawing.adv.multi_record_positions.length : 0;
    const maxRow = Math.max(1, count);
    return Array.from({ length: maxRow }, (_unused, index) => index + 1);
  };

  /** @param {DrawingState} drawing */
  const sortCircularRecordPositionsByRow = (drawing) => {
    if (!Array.isArray(drawing.adv.multi_record_positions) || drawing.adv.multi_record_positions.length <= 1) return;
    const sorted = drawing.adv.multi_record_positions
      .map((/** @type {MultiRecordPosition} */ entry, /** @type {number} */ index) => ({ ...entry, __index: index }))
      .sort((/** @type {MultiRecordPosition & { __index: number }} */ left, /** @type {MultiRecordPosition & { __index: number }} */ right) => {
        const leftRow = Number(left.row);
        const rightRow = Number(right.row);
        if (leftRow !== rightRow) return leftRow - rightRow;
        return left.__index - right.__index;
      })
      .map((/** @type {MultiRecordPosition & { __index: number }} */ { __index, ...entry }) => entry);
    drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length, ...sorted);
  };

  /** @param {number} index @param {unknown} rowValue */
  const setCircularRecordRow = (index, rowValue) => {
    const drawing = state.drawings.circular;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.adv.multi_record_positions.length) return;
    const target = drawing.adv.multi_record_positions[idx];
    if (!target || typeof target !== 'object' || Array.isArray(target)) return;
    const maxRow = Math.max(1, drawing.adv.multi_record_positions.length);
    const parsedRow = Number(rowValue);
    const normalizedRow = Number.isInteger(parsedRow) && parsedRow > 0 ? Math.min(parsedRow, maxRow) : 1;
    target.row = normalizedRow;
    sortCircularRecordPositionsByRow(drawing);
  };

  /** @param {number} index */
  const canMoveCircularRecordOrderUp = (index) => {
    const drawing = state.drawings.circular;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx <= 0 || idx >= drawing.adv.multi_record_positions.length) return false;
    const currentRow = getCircularRecordRow(drawing.adv.multi_record_positions[idx]);
    const prevRow = getCircularRecordRow(drawing.adv.multi_record_positions[idx - 1]);
    return currentRow === prevRow;
  };

  /** @param {number} index */
  const canMoveCircularRecordOrderDown = (index) => {
    const drawing = state.drawings.circular;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.adv.multi_record_positions.length - 1) return false;
    const currentRow = getCircularRecordRow(drawing.adv.multi_record_positions[idx]);
    const nextRow = getCircularRecordRow(drawing.adv.multi_record_positions[idx + 1]);
    return currentRow === nextRow;
  };

  /** @param {number} index */
  const moveCircularRecordOrderUp = (index) => {
    const drawing = state.drawings.circular;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const idx = Number(index);
    if (!canMoveCircularRecordOrderUp(idx)) return;
    const next = [...drawing.adv.multi_record_positions];
    const temp = next[idx - 1];
    next[idx - 1] = next[idx];
    next[idx] = temp;
    drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length, ...next);
  };

  /** @param {number} index */
  const moveCircularRecordOrderDown = (index) => {
    const drawing = state.drawings.circular;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const idx = Number(index);
    if (!canMoveCircularRecordOrderDown(idx)) return;
    const next = [...drawing.adv.multi_record_positions];
    const temp = next[idx + 1];
    next[idx + 1] = next[idx];
    next[idx] = temp;
    drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length, ...next);
  };

  const resetCircularRecordOrder = () => {
    const drawing = state.drawings.circular;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const defaults = buildDefaultCircularRecordPositions();
    drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length, ...defaults);
  };

  /** @param {DrawingState} drawing @param {readonly Partial<LinearSeq>[]} items */
  const applyLinearSeqMutation = (
    drawing, items,
    {
      preserveLosatCacheInfo = false,
      layoutEntries = drawing.linearRecordRows,
      alignmentMutation = msg('source set changed.'),
      retiredRecordKeys = /** @type {string[]} */ ([])
    } = {}
  ) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const depthWidth = linearDepthLogicalWidth();
    const next = normalizeLinearSeqList(items);
    if (depthWidth > 0) {
      next.forEach((seq) => {
        seq.depth = padDepthFileSlots(seq.depth, depthWidth);
      });
    }
    // Removing a File can remove the last source of a Depth series.
    retireLegendStylesOfUnnamedCaptions(() => linearTrackSlotEditor.changeLinearDepthSources(() => {
      linearSeqs.splice(0, linearSeqs.length, ...next);
    }));
    const activeUids = new Set(next.map((seq) => seq.uid));
    pendingLinearRecordExpansions.forEach((_, uid) => {
      if (!activeUids.has(uid)) pendingLinearRecordExpansions.delete(uid);
    });
    pendingLinearUploads.forEach((uid) => {
      if (!activeUids.has(uid)) pendingLinearUploads.delete(uid);
    });
    const nextRows = reconcileLinearRecordLayout(linearSeqs, layoutEntries);
    drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...nextRows);
    // A card that stays but now reads another file (`retiredRecordKeys`) is ON.
    const retired = new Set(retiredRecordKeys);
    const recordsOff = pruneRecordsOff(drawing.recordsOff, new Set([...activeUids].filter((uid) => !retired.has(uid))));
    if (recordsOff.length !== drawing.recordsOff.length) drawing.recordsOff.splice(0, Infinity, ...recordsOff);
    replaceLinearComparisonPlan(
      drawing,
      reconcileLinearComparisonPlan(drawing.linearComparisonPlan, linearSeqs),
      { invalidate: false }
    );
    invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo });
    linearReorderNotice.value = '';
    if (alignmentMutation === 'stable-reorder') {
      similarityAlignmentActions?.retainForStableReorder?.(
        linearSeqs.map((sequence) => sequence.uid)
      );
    } else {
      similarityAlignmentActions?.clearForMutation?.(alignmentMutation);
    }
  };

  const addLinearSeq = () => {
    const drawing = state.drawings.linear;
    return applyLinearSeqMutation(drawing, [...linearSeqs, createLinearSeq()]);
  };

  const restoreLinearSourceRemovalFocus = async () => {
    await nextTick();
    const target = linearSourceRemovalReturnFocus.value;
    linearSourceRemovalReturnFocus.value = null;
    if (target?.isConnected && typeof target.focus === 'function') target.focus();
  };
  const closeLinearSourceRemovalDialog = ({ restoreFocus = true } = {}) => {
    linearSourceRemovalDialog.open = false;
    linearSourceRemovalDialog.sourceUid = '';
    linearSourceRemovalDialog.origin = '';
    if (restoreFocus) void restoreLinearSourceRemovalFocus();
    else linearSourceRemovalReturnFocus.value = null;
  };
  const focusLinearSourceRemovalDialog = async () => {
    await nextTick();
    /** @type {HTMLElement | null} */ (document.querySelector('[data-linear-source-removal-primary]'))?.focus();
  };
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
  // A folded card takes the focus on its header line and stays folded.
  /** @param {number} sourceIndex */
  const focusLinearSourceAfterRemoval = async (sourceIndex) => {
    await nextTick();
    const cards = document.querySelectorAll('[data-linear-source-card]');
    const targetIndex = Math.max(0, Math.min(Number(sourceIndex) || 0, cards.length - 1));
    const target = /** @type {HTMLElement | null} */ (cards[targetIndex]?.querySelector(
      '[data-linear-source-fold]:not([open]) > summary, .upload-zone[tabindex="0"]'
    ) || document.querySelector('[data-linear-file-add]'));
    if (typeof target?.focus === 'function') target.focus();
  };
  /** @param {LinearSourceGroup} source @param {string} origin @param {EventTarget | null} [returnFocus] */
  const openLinearSourceRemovalDialog = (source, origin, returnFocus = null) => {
    if (!source?.uid) return false;
    linearSourceRemovalDialog.sourceUid = source.uid;
    linearSourceRemovalDialog.origin = origin;
    linearSourceRemovalDialog.open = true;
    linearSourceRemovalReturnFocus.value = returnFocus;
    void focusLinearSourceRemovalDialog();
    return true;
  };
  const cancelLinearSourceRemoval = () => closeLinearSourceRemovalDialog();
  /**
   * Deletes a Linear File, or clears it to one blank File (`clear`), inside
   * the caller's History step: the removal dialog and Remove File (D-06).
   * @param {Extract<LinearSourceRemovalPlan, { allowed: true }>} plan
   * @param {'clear' | 'delete'} intent
   */
  const removeLinearSource = (plan, intent) => {
    const next = [...plan.retainedSequences];
    if (intent === 'clear') next.splice(plan.insertionIndex, 0, createLinearSeq());
    return applyLinearSeqMutation(state.drawings.linear, next);
  };
  /** @param {'clear' | 'delete'} intent */
  const applyLinearSourceRemoval = async (intent) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (linearSourceRemovalDialog.origin === 'global' && intent !== 'delete') return false;
    const plan = planLinearSourceRemoval({
      sequences: linearSeqs,
      sourceUid: linearSourceRemovalDialog.sourceUid,
      intent
    });
    if (!plan.allowed) {
      closeLinearSourceRemovalDialog();
      return false;
    }
    const operation = await history.runUndoable(
      intent === 'clear' ? msg('Clear Linear File') : msg('Delete Linear File'),
      () => removeLinearSource(plan, intent)
    );
    closeLinearSourceRemovalDialog({ restoreFocus: false });
    await focusLinearSourceAfterRemoval(plan.sourceIndex);
    return operation;
  };
  /** @param {LinearSourceGroup} source @param {EventTarget | null} [returnFocus] */
  const requestLinearSourceRemoval = (source, returnFocus = null) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!source || !linearSourceHasPrimaryInput(source)) return false;
    return openLinearSourceRemovalDialog(source, 'card', returnFocus);
  };
  /** @param {{ currentTarget?: EventTarget | null } | null} [event] */
  const removeLastLinearSeq = (event = null) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const source = linearSourceGroups.value.at(-1);
    if (!source || linearSourceGroups.value.length <= 1) return false;
    if (isPristineLinearSource(source)) {
      linearSourceRemovalDialog.sourceUid = source.uid;
      linearSourceRemovalDialog.origin = 'global';
      return applyLinearSourceRemoval('delete');
    }
    return openLinearSourceRemovalDialog(source, 'global', event?.currentTarget || null);
  };

  /** @param {number} index @param {'gb' | 'gff' | 'fasta'} field @param {FileLike | null | undefined} value */
  const setLinearSeqPrimaryFile = (index, field, value) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= linearSeqs.length) return;
    if (!['gb', 'gff', 'fasta'].includes(field)) return;

    const nextValue = value ?? null;
    const seq = linearSeqs[idx];
    if (seq[field] === nextValue) return;
    // `groupLinearSourceRecords` puts every record card of `linearSeqs` in a group, and `seq` is one of them.
    const group = /** @type {LinearSourceGroup} */ (linearSourceGroups.value.find((entry) => (
      entry.records.some(({ sequence }) => sequence.uid === seq.uid)
    )));
    const members = new Set(group.records.map(({ sequence }) => sequence.uid));
    const replacement = createLinearSeq({
      ...group.sequence,
      [field]: nextValue,
      ...(field === 'gb' && nextValue ? { file_definition: '', file_subtitle: '' } : {}),
      ...(!seq[field] ? {} : sourceReplacementRetiresPresentation({
        removed: !nextValue, previousRecordCount: group.records.length
      }) ? RETIRED_LINEAR_RECORD_PRESENTATION : RETIRED_LINEAR_RECORD_REFERENCE)
    });
    const keepSource = field === 'gb' ? Boolean(nextValue) : Boolean(replacement.gff || replacement.fasta);
    applyLinearSeqMutation(drawing, linearSeqs.flatMap((entry) => (
      entry.uid === group.uid ? (keepSource ? [replacement] : [])
        : members.has(entry.uid) ? [] : [entry]
    )), { alignmentMutation: msg('source replaced.'), retiredRecordKeys: [...members] });
    if (keepSource) {
      pendingLinearRecordExpansions.set(replacement.uid, !seq[field]);
      pendingLinearUploads.add(replacement.uid);
    }
  };

  // Record selection (D-01..D-10): the record sources and edits its owner
  // reads, and the transitions it asks the root for.
  /** @param {LinearSeq} seq The discovered record a Linear card draws. */
  const linearCardRecord = (seq) => {
    const selection = resolveDisambiguatedRecordSelection(linearRecordSelector.recordsFor(seq), seq.region_record_id);
    return selection.record || (selection.entries.length === 1 ? selection.entries[0] : null);
  };
  /** @param {LinearSeq} seq @returns {RecordSettingsOwner} */
  const linearRecordSettingsOwner = (seq) => ({
    key: seq.uid, requestKeys: [seq.uid], ownsExpansions: true,
    bindingKeys: annotationRecordKeysOf(linearAnnotationSource(seq, false), seq.uid),
    displaySource: seq.uid, displaySelector: null
  });
  const circularRecordEntries = () => buildDisambiguatedRecordEntries(
    getCircularRecordDiscoveryState().records.map((record) => ({
      ...record, recordId: record?.record_id ?? record?.recordId, recordLength: record?.record_length ?? record?.recordLength
    }))
  );
  /** @returns {RecordSettingsOwner[]} */
  const circularRecordSettingsOwners = () => {
    const discovery = getCircularRecordDiscoveryState();
    const bindingKeys = annotationRecordKeysOf(
      { sourceKey: circularAnnotationSourceKey(), hasInput: discovery.hasInput, status: discovery.status, records: discovery.records },
      'circular-source'
    );
    return circularRecordEntries().map((entry, index) => ({
      key: circularRecordDrawKey(entry),
      requestKeys: circularRecordRequestKeys(entry),
      bindingKeys: bindingKeys[index] ? [bindingKeys[index]] : [],
      displaySource: 'circular', displaySelector: circularRecordDrawKey(entry)
    }));
  };
  /** @param {DrawingState} drawing @returns {RecordEditData} */
  const recordEditData = (drawing) => ({
    featureOverrides: drawing.featureOverrides,
    featurePlacementOverrides: drawing.featurePlacementOverrides,
    featureStrokeOverrides: drawing.featureStrokeOverrides,
    annotationSets: drawing.annotationSets,
    recordDisplayDrafts: drawing.recordDisplayDrafts,
    comparisonEdges: drawing === state.drawings.linear ? drawing.linearComparisonPlan.edges : [],
    annotationBindingField: ANNOTATION_RECORD_BINDING_KEY
  });
  /** @param {RecordSelectionMode} listMode @param {string} sourceKey @returns {RecordSource | null} */
  const recordSourceOf = (listMode, sourceKey) => {
    if (listMode === 'circular') {
      if (sourceKey !== 'circular' || getCircularRecordDiscoveryState().status !== 'ready') return null;
      const owners = circularRecordSettingsOwners();
      if (owners.length === 0) return null;
      const edited = recordsWithEdits(owners, recordEditData(state.drawings.circular));
      const entries = circularRecordEntries();
      return {
        key: 'circular',
        name: [cInputType.value === 'gff' ? files.c_gff : files.c_gb, cInputType.value === 'gff' ? files.c_fasta : null]
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
    const edited = recordsWithEdits(cards.map(linearRecordSettingsOwner), recordEditData(state.drawings.linear));
    return {
      key: group.uid,
      name: linearSourceFileName(group),
      records: cards.map((seq) => {
        const record = linearCardRecord(seq);
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
    const drawing = state.drawings[listMode];
    const deleted = new Set(keys);
    const owners = listMode === 'circular'
      ? circularRecordSettingsOwners().filter((owner) => deleted.has(owner.key))
      : linearSeqs.filter((seq) => deleted.has(seq.uid)).map(linearRecordSettingsOwner);
    retireLegendStylesOfUnnamedCaptions(() => removeRecordEdits(owners, recordEditData(drawing)));
    if (listMode === 'linear') {
      applyLinearSeqMutation(drawing, linearSeqs.map((seq) => (
        deleted.has(seq.uid) ? { ...seq, ...RECORD_SETTINGS_DEFAULTS } : seq
      )), { alignmentMutation: msg('record settings deleted.') });
    }
    return true;
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
      if (cInputType.value === 'gff') {
        setCircularSourceFile('c_gff', null);
        return setCircularSourceFile('c_fasta', null);
      }
      return setCircularSourceFile('c_gb', null);
    }
    const intent = linearSourceGroups.value.length > 1 ? 'delete' : 'clear';
    const plan = planLinearSourceRemoval({ sequences: linearSeqs, sourceUid: sourceKey, intent });
    return plan.allowed ? removeLinearSource(plan, intent) : false;
  };
  const recordSelection = createRecordSelection({
    reactive,
    computed,
    recordsOff: (listMode) => state.drawings[listMode].recordsOff,
    source: recordSourceOf,
    runUndoable: history.runUndoable,
    withDialogChoice: dialogChoice.withHistory,
    closeAfterDialogChoice: dialogChoice.closeAfterChoice,
    // A toggle invalidates comparison results and the alignment plan until Generate (P-9).
    afterRecordSetChange: (listMode) => {
      if (listMode !== 'linear') return;
      invalidateLinearComparisonArtifacts(state.drawings.linear, { preserveLosatCacheInfo: true });
      similarityAlignmentActions?.clearForMutation?.(msg('record drawing changed.'));
    },
    removeSource: removeRecordSource,
    deleteRecordSettings,
    // After a Session load the records are read on demand, as Generate reads them.
    prepareRecords: async (listMode) => {
      if (listMode === 'linear') return linearRecordSelector.refresh();
      return getCircularRecordDiscoveryState().status === 'ready' ? null : inspectCircularSourceRecords();
    },
    autoOpenRequests: recordListRequests
  });

  /** @param {DrawingState} drawing @param {number} sourceIndex @param {number} direction */
  const linearSourceMovePlan = (drawing, sourceIndex, direction) => planLinearSourceRowMove({
    sourceGroups: linearSourceGroups.value,
    entries: drawing.linearRecordRows,
    sourceIndex,
    direction
  });
  const linearSourceMoveBlockedReason = computed(() => {
    const drawing = state.activeDrawing();
    return (
      linearSourceMovePlan(drawing, 0, 1).reason === 'custom-layout'
        ? t('File order is unavailable because Record Layout is custom. Use Advanced comparison and layout → Record Layout to give each File its own consecutive rows.')
        : ''
    );
  });
  /** @param {number} sourceIndex @param {number} direction */
  const canMoveLinearSource = (sourceIndex, direction) => {
    const drawing = state.drawings.linear;
    return (
      linearSourceMovePlan(drawing, sourceIndex, direction).allowed
    );
  };

  /** @param {number} sourceIndex @param {number} direction */
  const moveLinearSource = (sourceIndex, direction) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const plan = linearSourceMovePlan(drawing, sourceIndex, direction);
    if (!plan.allowed) return;
    const next = moveLinearSourceGroup(linearSeqs, sourceIndex, direction);
    return history.runUndoable(msg('Move File'), () => {
      applyLinearSeqMutation(drawing, next, {
        preserveLosatCacheInfo: true,
        layoutEntries: plan.rows,
        alignmentMutation: 'stable-reorder'
      });
    });
  };

  /** @param {LinearSeq} sequence @param {unknown} value */
  const setLinearRecordSelector = (sequence, value) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!sequence) return false;
    const next = String(value || '');
    if (String(sequence.region_record_id || '') === next) return false;
    similarityAlignmentActions?.clearForMutation?.(msg('record selector changed.'));
    sequence.region_record_id = next;
    // The inferred definition follows the record the row now selects (D-12).
    if (lInputType.value === 'gb' && linearRecordSelector.statusFor(sequence) === 'ready') {
      sequence.inferred_definition = inferredDefinitionForRecord(linearRecordSelector.recordsFor(sequence), next);
    }
    return true;
  };

  /** @param {unknown} value */
  const setCircularInputType = (value) => {
    const busy = sessionOperationAvailability();
    if (busy) return busy;
    const next = value === 'gff' ? 'gff' : 'gb';
    if (cInputType.value === next) return { status: 'ok' };
    cInputType.value = next;
    // The other input's records are not the ones the OFF list names.
    const drawing = state.drawings.circular;
    if (drawing.recordsOff.length) drawing.recordsOff.splice(0);
    return { status: 'ok' };
  };

  /** @param {unknown} value */
  const setLinearInputType = (value) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = value === 'gff' ? 'gff' : 'gb';
    if (lInputType.value === next) return false;
    similarityAlignmentActions?.clearForMutation?.(msg('source type changed.'));
    lInputType.value = next;
    return true;
  };

  /** @param {LinearSeq} sequence @param {'region_start' | 'region_end'} field @param {unknown} value */
  const setLinearRecordCrop = (sequence, field, value) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!sequence || !['region_start', 'region_end'].includes(field)) return false;
    const next = value === '' || value === null || value === undefined
      ? null
      : Number(value);
    if (Object.is(sequence[field], next)) return false;
    similarityAlignmentActions?.clearForMutation?.(msg('record crop changed.'));
    sequence[field] = next;
    return true;
  };

  /** @param {LinearSeq} seq */
  const resetLinearRecordDefinition = (seq) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!seq) return;
    history.runUndoable(msg('Reset record definition'), () => {
      seq.definition = '';
    });
  };

  /** @param {LinearSeq} seq */
  const resetLinearRecordSubtitle = (seq) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!seq) return;
    history.runUndoable(msg('Reset record subtitle'), () => {
      seq.record_subtitle = '';
    });
  };

  const labelSettingsShown = computed(() => {
    const drawing = state.activeDrawing();
    const scope = mode.value === 'linear' ? drawing.form.show_labels_linear : drawing.form.labels_mode;
    return labelSettingsVisible(scope, drawing.featureOverrides);
  });

  // The template's names for the members of the shown mode's drawing: a ref
  // member reads and writes its value, any other member is the object itself.
  /** @param {keyof DrawingState} key */
  const drawingMember = (key) => (isRef(state.activeDrawing()[key])
    ? computed({
      get: () => state.activeDrawing()[key].value,
      /** @param {unknown} value */
      set: (value) => { state.activeDrawing()[key].value = value; }
    })
    : computed(() => state.activeDrawing()[key]));

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
    if (!circular && state.linearSeqs.length > 1 && linearComparisonGlobalAction.value !== 'none') {
      COMPARISON_COLOR_KEYS.forEach((key) => drawn.add(key));
    }
    resultsManager.userDefaultColorRows(drawing, toRaw(colors)).forEach(([key]) => drawn.add(key));
    const primary = keys.filter((key) => key === 'default' || drawn.has(key));
    return { primary, more: keys.filter((key) => !drawn.has(key) && key !== 'default') };
  });

  return {
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
    sessionImportPending,
    sessionSavePending,
    semanticMutationAvailable,
    dialogChoicePending,
    missingBlastTsvDialog,
    missingBlastTsvRows,
    chooseLinearComparisonPairFile,
    showLinearComparisonCause,
    runLinearComparisonIssueAction,
    resolveMissingBlastTsv,
    labelFitDialog,
    labelFitNotice,
    openLabelFitNotice,
    resolveLabelFitChoice,
    sessionSaveAvailable,
    sessionLoadAvailable,
    sessionBusyReason,
    generationCancelRequested,
    errorLog,
    errorDisplay,
    canRetryInteractiveSvgExport,
    reloadAfterOperationError,
    setDiagramMode,
    diagramModeSwitchAvailable,
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
    selectedResultIndex,
    failedGeneratePreservedResult,
    generationFailureRecovery,
    importedComparisonIntent: drawingMember('importedComparisonIntent'),
    importedComparisonNeedsResolution,
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
    cInputType,
    lInputType,
    losatProgram: drawingMember('losatProgram'),
    files,
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
    circularConservationLayoutWarning,
    circularConservationFastaInput,
    circularConservationSeriesRows,
    setCircularConservationUploadFiles,
    setCircularConservationCompanionFile,
    canMoveCircularConservationSeries,
    moveCircularConservationSeries,
    openCircularConservationComparisonFilePicker,
    addCircularConservationComparisonFile,
    removeCircularConservationSource,
    syncCircularConservationSeries,
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
    getDepthTrackLabel,
    setDepthTrackLabel,
    getDepthTrackColor,
    setDepthTrackColor,
    getDepthTrackLegendLabelForSlot,
    setDepthTrackLegendLabelForSlot,
    syncDepthTrackSlotLabel,
    addCircularDepthTrack,
    addLinearDepthTrack,
    removeCircularDepthTrack,
    removeLinearDepthTrack,
    getCircularDepthFile,
    setCircularDepthFile,
    getLinearDepthFile,
    setLinearDepthFile,
    setLinearSourceDepthFile,
    clearLinearSourceDepthFile,
    linearSeqs,
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
    linearRecordRowFor,
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
    linearRecordOptions: linearRecordSelector.optionsFor,
    refreshLinearRecordSelectors: linearRecordSelector.refresh,
    linearRecordSelectorDisabled: linearRecordSelector.isDisabled,
    linearRecordSelectorError: linearRecordSelector.errorFor,
    linearRecordSelectorWarning: linearRecordSelector.warningFor,
    form: drawingMember('form'),
    adv: drawingMember('adv'),
    linearTypographyLinked: drawingMember('linearTypographyLinked'),
    setLinearTypographyLinked: linearTypography.setLinked,
    setLinearRulerLabelFontSize: linearTypography.setRulerLabelFontSize,
    setLinearScaleFontSize: linearTypography.setScaleFontSize,
    unmanagedConfigOverrideEntries,
    resetUnmanagedConfigOverride,
    comparisonProfileDefault,
    comparisonHeightValidationError,
    optionalNumberInputValue,
    setOptionalNumberInputValue,
    autoValueText: autoValueDisplay.autoValueText,
    autoValueVisible: autoValueDisplay.autoValueVisible,
    canUseLinearRulerOnAxis,
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
    removeCircularTrackSlot: circularTrackSlotEditor.removeCircularTrackSlot,
    setCircularTrackSlotEnabled: circularTrackSlotEditor.setCircularTrackSlotEnabled,
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
    removeLinearTrackSlot: linearTrackSlotEditor.removeLinearTrackSlot,
    setLinearTrackSlotEnabled: linearTrackSlotEditor.setLinearTrackSlotEnabled,
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
    // Similarity groups are the Linear drawing's.
    orthogroupNameOverrides: state.drawings.linear.orthogroupNameOverrides,
    orthogroupDescriptionOverrides: state.drawings.linear.orthogroupDescriptionOverrides,
    selectedOrthogroupId,
    orthogroupSearch,
    orthogroupSortMode,
    showRightDrawer,
    rightDrawerTab,
    orthogroupCount: orthogroupActions.orthogroupCount,
    filteredOrthogroups: orthogroupActions.filteredOrthogroups,
    selectedOrthogroup: orthogroupActions.selectedOrthogroup,
    selectedOrthogroupDetail: orthogroupActions.selectedOrthogroupDetail,
    resolveOrthogroupName: orthogroupActions.resolveOrthogroupName,
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
    circularRecordList,
    refreshCircularRecordOrder,
    waitForAuxiliaryFileImport: auxiliaryFileImports.wait, canRetryAuxiliaryImportFailure,
    retryAuxiliaryImportFailure: () => history.runUndoableCheckpoint(msg('Change uploaded file'), retryAuxiliaryImportFailure, { shouldCommit: result => result !== false }),
    circularRecordPresentationOptions,
    circularRecordPresentationError,
    circularSingleRecordPresentationEnabled,
    circularRecordDiscoveryState,
    circularRecordInspectionEnabled,
    inspectCircularSourceRecords,
    circularRecordSelectionEnabled,
    showCircularCanvasSetting,
    setCircularRecordPresentationSelector,
    setCircularSourceFile,
    recordSelection,
    circularDrawnRecordCount,
    formatRecordLength,
    paletteDefinitions,
    paletteNames,
    selectedPalette: drawingMember('selectedPalette'),
    currentColors: drawingMember('currentColors'),
    paletteInstantPreviewEnabled,
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
    addWhitelistRule,
    removeWhitelistRule,
    editWhitelistRule,
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
    featureExtractionError,
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
    canEditFeatureColor,
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
    runInfoHasCliHelperFiles,
    resetSettings,
    requestResetSettings,
    saveSessionWithTitle,
    editSessionTitle,
    importSession,
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
    manualPriorityRules: drawingMember('manualPriorityRules'),
    newPriorityRule,
    addPriorityRule
  };
};
