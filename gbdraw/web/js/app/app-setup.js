// @ts-check
/** @import { DrawingState } from '../state.js' */
/** @import { RulePreparation } from './rule-matching.js' */
/** @import { FeatureEditorOptions } from './feature-editor.js' */
/** @import { ColorActionsRuleActions } from './feature-editor/color-actions.js' */
/** @import { UserFacingError } from '../utils/error-normalization.js' */
/** @import { AnnotationCatalogSource } from './annotations/record-catalog.js' */
/** @import { LegacyResultSvgTransform } from '../services/config.js' */
/** @import { ArtifactSlot } from '../services/artifact-slot.js' */
/** @typedef {{ opening: Readonly<ArtifactSlot> | null, stashed: Readonly<ArtifactSlot> | null }} LoadedArtifactSlots */
import { createRulePreparation } from './rule-matching.js';
import { compileDirectEditorMutationPlan } from './candidate-render.js';
import {
  countUnresolvedFeatureEdits, removeUnresolvedFeatureEdits, requestFeatureVisibilityRules
} from '../services/feature-visibility.js';
import { isLegendOrderEdited } from '../services/legend-svg.js';
import { admitFeatureCatalog } from '../services/feature-catalog.js';
import { createDefaultLosatpHitLimits } from '../services/session-active-config-contract.js';
import { createRecordDisplayControls } from './record-display-options.js';
import { isCurrentFeature } from '../services/feature-identity.js';
import {
  createFeatureRecordRotationAction,
  createFeatureRecordRotationWorkflow
} from './record-display/feature-record-rotation.js';
import { state, sessionOperationAvailability, createLinearSeq, normalizeLinearSeqList } from '../state.js';
import {
  adoptCanonicalRenderArtifacts,
  applyConfigData,
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
  readCommittedResourceRecordCount,
  assertActiveModeInputs,
  importSession as importSessionFromFile,
  SESSION_VERSION,
  serializeActiveRenderFiles,
  serializeResults,
  setUnmanagedConfigOverrideValidator,
  setPreviewRuntime
} from '../services/config.js';
import { setMainSessionComparisonFrameConverter } from '../services/main-session-comparison-frame.js';
import { createHistoryManager } from '../services/history.js';
import { createHistoryFileStore } from '../services/history-files.js';
import { createHistorySnapshotService } from '../services/history-snapshot.js';
import { cloneJsonData } from '../services/json-clone.js';
import {
  groupLinearSourceRecords,
  isPristineLinearSource,
  linearSourceHasPrimaryInput,
  linearSourceDepthStatus,
  moveLinearSourceGroup,
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
import { downloadTextFile } from '../services/text-download.js';
import { resetLayoutState, resetSettings as resetSettingsState } from '../services/reset.js';
import {
  DIAGRAM_HELPER_OPERATIONS,
  disposeDiagramGenerationWorker,
  runDiagramHelperOperation
} from '../services/diagram-generation.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from '../services/runtime-test-hooks.js';
import { createPanZoom, createSidebarResize, setupGlobalUiEvents } from './ui.js';
import { colorValueMode, toNativeColorInputValue } from '../utils/color-utils.js';
import { createFeatureEditor } from './feature-editor.js';
import { PAIRWISE_MATCH_SELECTOR } from './pairwise-match-popup.js';
import { createFeatureSelection } from './feature-selection.js';
import { createPreviewFeatureSearch } from './feature-search/preview-actions.js';
import { createSvgStyles } from './svg-styles.js';
import { createLegendManager } from './legend.js';
import { createPaletteLoader } from './palettes.js';
import { afterPaint, createRunAnalysis } from './run-analysis.js';
import { createSimilarityAlignmentActions } from './similarity-alignment.js';
import { diagnosticError, normalizeUserFacingError } from '../utils/error-normalization.js';
import { formatElapsedMs, reproducibilityLabel } from './run-info.js';
import { createLegendLayout } from './legend-layout.js';
import {
  COMPOSITION_METADATA_ATTRIBUTE,
  COMPOSITION_SCHEMA_ATTRIBUTE,
  normalizeLegacyComposition
} from './legend-layout/composition-actions.js';
import { applyStrokeOverridesToSvg } from './legend/stroke-actions.js';
import { createResultsManager } from './results.js';
import { setupWatchers } from './watchers.js';
import { setupHistoryInputs } from './history-inputs.js';
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
import { createLinearRecordSelector } from './linear-record-selector.js';
import { createLinearTypographyController } from './linear-typography.js';
import {
  buildDisambiguatedRecordEntries,
  formatRecordLength,
  resolveCircularRequestRecordSet,
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
  linearComparisonEdgeKey,
  materializeResolvedEdgesAsSelectedPlan,
  normalizeLinearComparisonPlan,
  plainTextLinearRecordLabel,
  reconcileLinearComparisonPlan
} from '../services/linear-comparisons.js';
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
  reconcileConservationSeries
} from '../services/conservation-series.js';
import {
  getDepthTrackFallbackLabel,
  getDepthTrackLabelFromFile,
  isDepthTrackAutoLabel
} from './depth-tracks.js';
import { comparisonProfileDefault } from '../mode-profiles.js';
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

const { onMounted, onUnmounted, watch, nextTick, computed, ref, reactive, isRef } = window.Vue;
const toRaw = window.Vue.toRaw || ((value) => value);

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
  restore: async (snapshot) => {
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

// The transform services/config.js applies to each Result of an older Session
// before it commits (R13 port): a Result without composition metadata gets the
// legacy composition, and the saved strokes are projected into it.
/** @type {LegacyResultSvgTransform} */
export const transformLegacyResultSvg = (svg, { composition, strokes }) => {
  let compositionChanged = false;
  if (
    svg.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) === null
    && svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE) === null
  ) {
    normalizeLegacyComposition(svg, composition);
    compositionChanged = true;
  }
  const strokeCount = applyStrokeOverridesToSvg({ svg, ...strokes });
  return compositionChanged || strokeCount > 0;
};

const HISTORY_RESTORE_BUSY = Object.freeze({ status: 'busy', reason: 'Undo or Redo in progress. Retry after it finishes.' });
// E1: a History step that switches the diagram mode; its changes name `ui.mode`.
/** @param {unknown} changes */
const historyStepSwitchesMode = (changes) => (Array.isArray(changes) ? changes : [])
  .some((/** @type {{ path?: unknown[] }} */ { path } = {}) => path?.[0] === 'ui' && path[1] === 'mode');

export const createAppSetup = () => {
  setUnmanagedConfigOverrideValidator((payload) => runDiagramHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.VALIDATE_CONFIG_OVERRIDES,
    payload
  ));
  setMainSessionComparisonFrameConverter((payload) => runDiagramHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.CONVERT_MAIN_SESSION_COMPARISON_FRAME,
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
    resetColorDialog,
    labelTextScopeDialog,
    hiddenLabelTextDialog,
    labelOnDialog,
    sidebarWidth,
    originalLegendOrder,
    newLegendCaption,
    newLegendColor,
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
  /** @type {{ refreshCanvas: () => void, reviewBlocksEditor: () => boolean }} */
  const similarityAlignmentPorts = {
    refreshCanvas: () => {},
    reviewBlocksEditor: () => false
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
    const indexByUid = new Map(
      linearSeqs.map((sequence, index) => [String(sequence?.uid || ''), index])
    );
    const resolvedByEdgeKey = new Map(
      drawing.linearComparisonResolution.value.edges.map((edge) => [edge.edgeKey, edge])
    );
    losatCacheInfo.value = losatCacheInfo.value.flatMap((entry) => {
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
      else losatCacheInfo.value = losatCacheInfo.value.filter((entry) => !entry?.edgeKey);
    }
  };

  /** @param {DrawingState} drawing */
  const replaceLinearComparisonPlan = (drawing, nextPlan, { invalidate = true } = {}) => {
    const normalized = normalizeLinearComparisonPlan(nextPlan);
    const changed = normalized.mode !== drawing.linearComparisonPlan.mode
      || normalized.defaultSource !== drawing.linearComparisonPlan.defaultSource
      || normalized.edges.length !== drawing.linearComparisonPlan.edges.length
      || normalized.edges.some((edge, index) => (
        !sameLinearComparisonEdge(edge, drawing.linearComparisonPlan.edges[index])
      ));
    if (invalidate && changed) {
      similarityAlignmentActions?.clearForMutation?.('comparison configuration changed.');
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

  /** @param {DrawingState} drawing */
  const mutateLinearComparisonPlan = (drawing, mutator) => history.runUndoable('Change comparisons', () => {
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
  const setLinearRecordRow = (uid, row) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    syncLinearRecordLayout({ preserveLosatCacheInfo: true });
    const previous = drawing.linearRecordRows.find((entry) => entry.uid === uid)?.row;
    updateLinearRecordRow(drawing.linearRecordRows, uid, row);
    if (drawing.linearRecordRows.find((entry) => entry.uid === uid)?.row !== previous) {
      invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo: true });
    }
  };
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
      adjacentEdgeKeys: adjacentRowPairs(linearSeqs, effectiveLinearComparisonLayout(drawing), true)
        .map(([queryUid, subjectUid]) => linearComparisonEdgeKey(queryUid, subjectUid)),
      losatProgram: drawing.losatProgram.value,
      blastpMode: drawing.losat.blastp?.mode,
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
  const canRunLinearLosat = computed(() => linearSeqs.filter((sequence) => (
    lInputType.value === 'gff'
      ? sequence.gff && sequence.fasta
      : sequence.gb
  )).length >= 2);

  const linearSourceGroups = computed(() => groupLinearSourceRecords(linearSeqs));
  const linearSourceRemovalDialog = reactive({ open: false, sourceUid: '', origin: '' });
  const linearSourceRemovalReturnFocus = ref(null);
  const linearSourceRemovalTarget = computed(() => (
    linearSourceGroups.value.find((source) => source.uid === linearSourceRemovalDialog.sourceUid) || null
  ));
  const linearSourceRemovalCanDelete = computed(() => linearSourceGroups.value.length > 1);
  const linearSourceRemovalTargetName = computed(() => {
    const source = linearSourceRemovalTarget.value;
    if (!source) return 'Unavailable File';
    const sequence = source.sequence || source.records?.[0]?.sequence || {};
    const names = [sequence.gb, sequence.gff, sequence.fasta]
      .filter(Boolean)
      .map((file) => String(file?.name || 'Unnamed file'));
    return names.length ? names.join(' + ') : `File ${linearSourceGroups.value.indexOf(source) + 1}`;
  });
  const linearComparisonTimeline = computed(() => {
    const drawing = state.activeDrawing();
    return buildLinearComparisonTimeline({
      sequences: linearSeqs,
      layout: effectiveLinearComparisonLayout(drawing),
      plan: drawing.linearComparisonPlan,
      resolution: drawing.linearComparisonResolution.value
    });
  });
  const linearComparisonPairForEdgeKey = (edgeKey) => {
    for (const row of linearComparisonTimeline.value.rows) {
      const pair = row.boundaryAfter?.pairs.find((entry) => entry.edgeKey === edgeKey);
      if (pair) return pair;
    }
    return null;
  };
  const linearComparisonRecordLabel = (uid) => {
    const sequence = linearSeqs.find((entry) => entry.uid === uid);
    return plainTextLinearRecordLabel(
      sequence?.definition ||
      sequence?.region_record_id ||
      sequence?.gb?.name ||
      sequence?.gff?.name ||
      sequence?.fasta?.name ||
      'Record'
    );
  };
  const openLinearComparisonDisclosure = async (disclosureKey) => {
    if (!disclosureKey) return null;
    const details = /** @type {HTMLDetailsElement[]} */ ([...document.querySelectorAll('[data-linear-comparison-disclosure]')])
      .find((element) => element.dataset.linearComparisonDisclosure === disclosureKey);
    if (!details) return null;
    details.open = true;
    await nextTick();
    return details;
  };
  watch([mode, () => state.drawings.linear.hasActiveLinearLosatIntent.value], ([activeMode, activeLosat]) => {
    if (activeMode === 'linear' && activeLosat) openLinearComparisonDisclosure('settings');
  }, { flush: 'post' });

  const focusLinearComparisonPair = async (edgeKey) => {
    await openLinearComparisonDisclosure('selected-pairs');
    const container = /** @type {HTMLElement[]} */ ([...document.querySelectorAll('[data-edge-key]')])
      .find((element) => element.dataset.edgeKey === edgeKey);
    /** @type {HTMLElement | null | undefined} */ (container?.querySelector('input, select, button'))?.focus();
  };

  const focusLinearComparisonIssue = async () => {
    const target = linearComparisonUi.value.errorTargets[0];
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
    /** @type {HTMLElement | null | undefined} */ (focusable || container?.querySelector?.('input, select, button, [tabindex]'))?.focus();
    return true;
  };

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

  const setLinearComparisonLosatMode = (modeKey) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const selection = projectLinearComparisonLosatModeSelection({ modeKey });
    if (!selection.selectable || !selection.patch) return false;
    const nextProgram = selection.patch.losatProgram;
    if (drawing.losatProgram.value === nextProgram) return true;
    similarityAlignmentActions?.clearForMutation?.('comparison program changed.');
    drawing.losatProgram.value = nextProgram;
    invalidateLinearComparisonArtifacts(drawing);
    return true;
  };

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
    similarityAlignmentActions?.clearForMutation?.('comparison mode changed.');
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

  const findEdgeIndex = (edges, id) => edges.findIndex((edge) => edge.id === id);
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
    if (linearSeqs.length < 2) return;
    syncLinearRecordLayout();
    const [firstPair] = adjacentRowPairs(
      linearSeqs,
      effectiveLinearComparisonLayout(drawing)
    );
    const [queryUid, subjectUid] = firstPair || [linearSeqs[0].uid, linearSeqs[1].uid];
    const next = selectedPlanForEdit(drawing);
    upsertSelectedComparison(next, { queryUid, subjectUid });
    replaceLinearComparisonPlan(drawing, next);
    await focusLinearComparisonPair(linearComparisonEdgeKey(queryUid, subjectUid));
  };
  const omitLinearComparison = (id) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = selectedPlanForEdit(drawing);
    const index = findEdgeIndex(next.edges, id);
    if (index < 0) return;
    const edge = next.edges[index];
    edge.included = false;
    if (!edge.file && !String(edge.losatFilename || '').trim()) next.edges.splice(index, 1);
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
  const setLinearComparisonFile = (id, file) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = selectedPlanForEdit(drawing);
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge.file = file || null;
    edge.fileActive = Boolean(file);
    edge.source = LINEAR_COMPARISON_SOURCES.UPLOAD;
    edge.included = Boolean(file) || edge.included;
    replaceLinearComparisonPlan(drawing, next);
  };
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
  /** @param {DrawingState} drawing */
  const updateResolvedLosatFilenameDraft = (drawing, edgeKey, updater) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const resolved = drawing.linearComparisonResolution.value.edges.find((edge) => edge.edgeKey === edgeKey);
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
  const setResolvedLinearComparisonLosatFilename = (edgeKey, value) => {
    const drawing = state.drawings.linear;
    updateResolvedLosatFilenameDraft(drawing, edgeKey, (edge) => {
      edge.losatFilename = String(value || '');
      edge.losatFilenameActive = Boolean(edge.losatFilename.trim());
    });
  };
  const reuseResolvedLinearComparisonLosatFilename = (edgeKey) => {
    const drawing = state.drawings.linear;
    updateResolvedLosatFilenameDraft(drawing, edgeKey, (edge) => {
      if (String(edge.losatFilename || '').trim()) edge.losatFilenameActive = true;
    });
  };
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
    adjacentRowPairs(linearSeqs, effectiveLinearComparisonLayout(drawing), allPairs).forEach(([queryUid, subjectUid]) => {
      upsertSelectedComparison(next, { queryUid, subjectUid });
    });
    replaceLinearComparisonPlan(drawing, next);
  };
  const setLinearComparisonGapAction = (edgeKey, action) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const pair = linearComparisonPairForEdgeKey(edgeKey);
    if (!pair) return;
    const next = selectedPlanForEdit(drawing);
    const index = findEdgeIndexForPair(next.edges, pair);
    if (action === 'none') {
      if (index >= 0) {
        next.edges[index].included = false;
        if (!next.edges[index].file && !String(next.edges[index].losatFilename || '').trim()) {
          next.edges.splice(index, 1);
        }
      }
      replaceLinearComparisonPlan(drawing, next);
      return;
    }
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
  const linearRecordRowFor = (uid, fallback) => {
    const drawing = state.drawings.linear;
    return drawing.linearRecordRows.find((entry) => entry.uid === uid)?.row || fallback;
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
      .filter((entry) => String(entry?.edgeKey || ''))
      .map((entry) => [String(entry.edgeKey), entry])
  ));

  const pendingLinearRecordExpansions = new Set();
  const pendingLinearMetadataInference = new Set();
  /** @param {DrawingState} drawing */
  const expandDiscoveredLinearRecords = (drawing, { uid, records, inferDefinitions = false }) => {
    const expanding = pendingLinearRecordExpansions.delete(uid);
    const index = linearSeqs.findIndex((seq) => seq.uid === uid);
    if (index < 0) return;
    const source = linearSeqs[index];
    if (!expanding || records.length < 2 || source.region_record_id) {
      if (inferDefinitions) {
        source.inferred_definition = inferredDefinitionForRecord(records, source.region_record_id);
      }
      return;
    }
    const row = linearRecordRowFor(uid, index + 1);
    // Expanded records keep only File-level values. A crop or record display
    // value of the replaced File does not apply to the new records (IN-04).
    const expanded = buildDisambiguatedRecordEntries(records).map((record, recordIndex) => createLinearSeq({
      uid: recordIndex === 0 ? uid : undefined,
      gb: source.gb,
      gff: source.gff,
      fasta: source.fasta,
      depth: source.depth,
      losat_gencode: source.losat_gencode,
      file_definition: source.file_definition,
      file_subtitle: source.file_subtitle,
      inferred_definition: inferDefinitions ? record.inferredDefinition || '' : '',
      region_record_id: record.value
    }));
    applyLinearSeqMutation(drawing, [
      ...linearSeqs.slice(0, index), ...expanded, ...linearSeqs.slice(index + 1)
    ], { alignmentMutation: 'record selector changed.' });
    expanded.forEach((seq) => updateLinearRecordRow(drawing.linearRecordRows, seq.uid, row));
    return true;
  };
  const handleLinearRecordsDiscovered = ({ uid, records }) => {
    const drawing = state.drawings.linear;
    const isRollbackOrSessionLoad = Boolean(
      state.sessionImportRollbackInProgress?.value ||
      state.sessionResourceDiscoveryDeferred?.value
    );
    // Only an upload infers record definitions; a loaded Session keeps its own.
    const inferDefinitions = !isRollbackOrSessionLoad && pendingLinearMetadataInference.delete(uid);
    return expandDiscoveredLinearRecords(drawing, { uid, records, inferDefinitions });
  };
  const materializeAutomaticLinearRecords = async () => {
    if (mode.value !== 'linear') return;
    linearSeqs.forEach((seq) => {
      if (!seq.region_record_id) pendingLinearRecordExpansions.add(seq.uid);
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
  /** @param {AnnotationCatalogSource[] | null} [linearSourcesOverride] */
  const getAnnotationRecordCatalog = (linearSourcesOverride = null) => {
    const drawing = state.drawings.circular;
    const circularPrimaryFile = cInputType.value === 'gff' ? files.c_gff : files.c_gb;
    const circularPairedFile = cInputType.value === 'gff' ? files.c_fasta : null;
    const circularDiscovery = getCircularRecordDiscoveryState();
    return buildAnnotationRecordCatalog(/** @type {any} */ ({
      mode: mode.value,
      circularSource: {
        sourceKey: annotationSourceKey({
          scope: 'circular',
          inputType: cInputType.value,
          primaryFile: circularPrimaryFile,
          pairedFile: circularPairedFile
        }),
        hasInput: circularDiscovery.hasInput,
        status: circularDiscovery.status,
        error: circularDiscovery.error,
        records: resolveCircularRequestRecordSet(/** @type {any} */ ({
          records: circularDiscovery.records,
          selector: drawing.form.circular_record_selector,
          multiRecordCanvas: drawing.form.multi_record_canvas,
          groupingIntent: drawing.adv.circular_grouping_intent
        })).records
      },
      linearSources: linearSourcesOverride || linearSeqs.map((seq) => {
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
          hasInput: Boolean(primaryFile && (lInputType.value !== 'gff' || pairedFile)),
          status: linearRecordSelector.statusFor(seq),
          error: linearRecordSelector.errorFor(seq),
          records: linearRecordSelector.recordsFor(seq)
        };
      })
    }));
  };
  const previewRuntime = createPreviewRuntime({ state, serializeSvg: serializeCleanSvg });
  setPreviewRuntime(previewRuntime);

  const historyFileStore = createHistoryFileStore();
  const historySnapshots = createHistorySnapshotService({
    state,
    fileStore: historyFileStore,
    nextTick,
    normalizeLinearSeqList,
    buildConfigData,
    applyConfigData,
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
    applyRunStateData
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
  const restoreWithSpecificRuleDrafts = async (restore, ...args) => {
    const drafts = specificRuleRestorePorts.captureSpecificRulePatternDrafts();
    try {
      return await restore(...args);
    } finally {
      if (drafts) specificRuleRestorePorts.restoreSpecificRulePatternDrafts(drafts);
    }
  };
  /** @param {DrawingState} drawing */
  const restoreRuleEdits = async (drawing, restore, ...args) => {
    const previousRules = drawing.manualSpecificRules.map((rule) => ({ ...rule }));
    specificRuleRestorePorts.retainRulesForRestore(previousRules);
    try {
      const restored = await restoreWithSpecificRuleDrafts(restore, ...args);
      specificRuleRestorePorts.followRestoredSpecificRules(previousRules);
      return restored;
    } finally {
      specificRuleRestorePorts.retainRulesForRestore([]);
    }
  };
  const history = createHistoryManager(/** @type {any} */ ({
    buildIntent: historySnapshots.buildHistoryIntent,
    applyIntent: (...args) => {
      const drawing = state.activeDrawing();
      return restoreRuleEdits(drawing, historySnapshots.applyHistoryIntent, ...args);
    },
    buildCheckpoint: historySnapshots.buildArtifactCheckpoint,
    applyCheckpoint: (...args) => {
      const drawing = state.activeDrawing();
      return restoreRuleEdits(drawing, historySnapshots.applyArtifactCheckpoint, ...args);
    },
    captureGeneratedArtifactHandle: historySnapshots.captureGeneratedArtifactHandle,
    restoreGeneratedArtifactHandle: (...args) => restoreWithSpecificRuleDrafts(historySnapshots.restoreGeneratedArtifactHandle, ...args),
    compareGeneratedArtifactHandles: historySnapshots.compareGeneratedArtifactHandles,
    signatureFor: historySnapshots.snapshotSignature,
    fileStore: historyFileStore,
    collectCurrentFileIds: historySnapshots.collectCurrentFileIds,
    makeRef: ref,
    mutationAvailability: sessionOperationAvailability,
    // An Undo or Redo of a mode switch waits as the mode buttons do (E1).
    stepAvailability: (/** @type {{ changes?: unknown }} */ { changes }) => (
      historyStepSwitchesMode(changes) ? diagramModeOperationBusy() : null
    )
  }));
  const recordDisplayControls = createRecordDisplayControls({ state, computed, watch,
    linearRecordsFor: linearRecordSelector.recordsFor,
    linearRecordStatusFor: linearRecordSelector.statusFor,
    linearRecordErrorFor: linearRecordSelector.errorFor,
    runUndoable: history.runUndoable,
    getCommittedRequest: getCommittedCanonicalRenderRequest, getCommittedSession: getCommittedCanonicalSession });
  window.__GBDRAW_HISTORY__ = history;
  // R13: owners receive this port, not the record display controls. The check
  // is pure (services/feature-identity.js) over the record display's binding.
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
    const label = history.undoLabel();
    return label ? `Undo ${label}` : 'Undo';
  });
  const redoHistoryTitle = computed(() => {
    void history.revision.value;
    const label = history.redoLabel();
    return label ? `Redo ${label}` : 'Redo';
  });
  const undoHistory = () => history.undo();
  const redoHistory = () => history.redo();
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
    previewTransformInteraction
  } = createPanZoom(state);
  const { startResizing } = createSidebarResize(state);

  const specificRuleNotice = ref('');
  const ruleMatchingPending = ref(false);
  // Python's rule evaluation (R7): the rule preparation's, and the Label TSV
  // import's own stateless one (`evaluateLabelRules`).
  const evaluateRules = async (payload, options) => (await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES, payload, options)).result;
  const rulePreparation = createRulePreparation({
    state,
    pending: ruleMatchingPending,
    notify: notice => { specificRuleNotice.value = notice; },
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
    readActiveResultIdentity: () => previewRuntime.getActiveRuntime()?.resultIdentity
  });
  // History captures register once their owner exists (R13).
  historySnapshots.registerCapture('legend', legendActions.captureLegendEntryOwners);
  // R13: the palette watcher reacts through the root's palette and rules
  // projection, registered once the style owner it applies through exists.
  /** @type {{ projectPaletteAndRules: FeatureEditorOptions['projectPaletteAndRules'] }} */
  // Null only until the root registers the port below; the palette watcher reads it after setup.
  const paletteRulePorts = { projectPaletteAndRules: /** @type {never} */ (null) };
  const svgActions = createSvgStyles({
    state,
    watch,
    nextTick,
    commitActiveResultEdit: previewRuntime.commitActiveResultEdit,
    projectPaletteAndRules: (...args) => paletteRulePorts.projectPaletteAndRules(...args)
  });
  // The palette and the specific-color rules on the mounted Result (R3): the
  // one call of their projection, shared by `projectMountedEditorIntent`, the
  // palette watcher, and a rule commit. It prepares the rule matches first
  // unless the caller prepared them (`prepareRules: false`, which applies at
  // once), and resolves to false when the rules changed meanwhile.
  const projectPaletteAndRules = ({ recolor = {}, prepareRules = true } = {}) => {
    const project = () => {
      svgActions.applyPaletteToSvg(recolor);
      svgActions.applySpecificRulesToSvg();
      return true;
    };
    return prepareRules
      ? Promise.resolve(rulePreparation.prepare()).then((prepared) => (prepared ? project() : false))
      : project();
  };
  paletteRulePorts.projectPaletteAndRules = projectPaletteAndRules;
  const featureSelection = createFeatureSelection(/** @type {any} */ ({ state, onMounted, onUnmounted }));
  const featureActions = createFeatureEditor({
    state,
    rulePreparation,
    evaluateLabelRules: evaluateRules,
    runUndoable: history.runUndoable,
    runUndoableCheckpoint: history.runUndoableCheckpoint,
    getCommittedRequest: getCommittedCanonicalRenderRequest,
    getCommittedSession: getCommittedCanonicalSession,
    readResourceRecordCount: readCommittedResourceRecordCount,
    readFeatureOverrideTable: (payload) => runDiagramHelperOperation(
      DIAGRAM_HELPER_OPERATIONS.READ_FEATURE_OVERRIDE_TABLE, payload
    ),
    isCurrentFeature: isCurrentResultFeature,
    isPatternEditAvailable: () => !sessionImportPending.value,
    nextTick,
    prepareFileLegendEntries: /** @type {any} */ (legendActions.prepareFileLegendEntries),
    extractLegendEntries: legendActions.extractLegendEntries,
    onLegendGeometryChanged: legendActions.onLegendGeometryChanged,
    featureSelection,
    commitActiveResultEdit: previewRuntime.commitActiveResultEdit,
    applyFeatureVisibilityChanges: previewRuntime.applyFeatureVisibilityChanges,
    selectResult,
    previewTransformInteraction,
    projectPaletteAndRules,
    projectFeatureEdits: () => projectMountedEditorIntent({ visibility: true, rerender: true, reflow: true, labels: true })
  });
  legendRowRulePorts.commitLegendRowRules = featureActions.commitSpecificRules;
  specificRuleRestorePorts.captureSpecificRulePatternDrafts = featureActions.captureSpecificRulePatternDrafts;
  specificRuleRestorePorts.restoreSpecificRulePatternDrafts = featureActions.restoreSpecificRulePatternDrafts;
  specificRuleRestorePorts.followRestoredSpecificRules = featureActions.followRestoredSpecificRules;
  // R13: the drawer and the feature search come after the owners they react
  // through, so each receives its ports directly.
  const rightDrawerActions = createRightDrawerController({ state, watch,
    onClose: featureActions.suspendSpecificRulePatternDrafts,
    focusReturn: {
      isFocusInDrawer: () => Boolean(document.querySelector('.right-drawer')?.contains(document.activeElement)),
      focusToggle: () => /** @type {HTMLElement | null} */ (document.querySelector('.drawer-toggle'))?.focus()
    },
    getOpenDisabledReason: () => similarityAlignmentPorts.reviewBlocksEditor()
      ? 'Finish or cancel alignment review before opening Editor.' : '' });
  const orthogroupActions = createOrthogroupEditor({ state });
  const previewFeatureSearch = createPreviewFeatureSearch({
    state,
    watch,
    nextTick,
    computed,
    isActiveResultReady: previewRuntime.isActiveResultReady,
    resolveOrthogroups: () => orthogroups.value.map((group) => ({
      ...group,
      display_name: orthogroupActions.resolveOrthogroupName(group),
      description: orthogroupActions.resolveOrthogroupDescription(group)
    })),
    openFeatureEditorForFeature: featureActions.openFeatureEditorForFeature
  });

  watch(selectedResultIndex, () => {
    featureSelection.clearFeatureSelection({ clearStatus: true, syncDom: false });
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
  const handleFeatureListScroll = (event) => {
    const target = event?.currentTarget || event?.target;
    featureListScrollTop.value = Number(target?.scrollTop || 0);
    const nextHeight = Number(target?.clientHeight || 0);
    if (nextHeight > 0) featureListViewportHeight.value = nextHeight;
  };
  watch(featureSearchInput, (value) => {
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
    closeRightDrawer: rightDrawerActions.closeRightDrawer
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
  // and label transitions below run through it. A retired rename reaches the
  // displayed Result through the one Legend projection (R3).
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
      projectLegendEntries: () => { void projectMountedEditorIntent({ legend: {} }); },
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
  const circularTrackSlotEditor = createCircularTrackSlotEditor({ state, changeTrackLayout });
  const linearTrackSlotEditor = createLinearTrackSlotEditor({ state, changeTrackLayout });
  const annotationImportNotice = ref('');
  const annotationEditor = createAnnotationEditor({
    state, getRecordCatalog: getAnnotationRecordCatalog, retireLegendStylesOfUnnamedCaptions,
    onImportNotice: (notice) => { annotationImportNotice.value = notice; }
  });
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
    state, state.drawings.circular
  ));
  const losatSettings = createLosatSettings({ state });
  const autoValueDisplay = createAutoValueDisplay(state);
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
  const depthFileCount = (value) => (
    isRecordMajorDepthFileMatrix(value)
      ? representativeDepthFiles(value).filter(Boolean).length
      : uploadedDepthFileCount(value)
  );
  const depthTrackCountLabel = (value) => {
    const count = depthFileCount(value);
    return count === 1 ? '1 TSV' : `${count} TSVs`;
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
  const sourceDepthTrackCount = (slots, uiCount = 1) => Math.max(
    1,
    Number(uiCount) || 1,
    representativeDepthFiles(slots).length
  );
  // The Depth panels read the drawing's series and never write them (OV-109):
  // the Depth transitions and the Session reconcile grow, trim and label the
  // list, and a series it does not hold yet shows its defaults.
  /** @param {DrawingState} drawing */
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
  const padLinearDepthRows = (width) => {
    const targetWidth = Math.max(0, Number(width) || 0);
    linearSeqs.forEach((seq) => {
      seq.depth = padDepthFileSlots(seq.depth, targetWidth);
    });
  };
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
  /** @param {DrawingState} drawing */
  const normalizeDepthTrackConfig = (drawing, entry, index) => (
    normalizeDepthTrackConfigEntry(entry, index, depthTrackConfigDefaults(drawing))
  );
  const optionalNumberInputValue = (value) => value ?? '';
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
  const linearSourceDepthRows = (source) => linearDepthTrackRows().map((track) => ({
    ...track,
    status: linearSourceDepthStatus(source, track.index)
  }));
  const linearSourceDepthSummary = (source) => {
    const tracks = linearSourceDepthRows(source);
    const attachedCount = tracks.filter((track) => track.status.selectedCount > 0).length;
    if (attachedCount === 0) return 'No depth track attached';
    return `${attachedCount} depth track${attachedCount === 1 ? '' : 's'} attached`;
  };
  const depthTrackRows = computed(() => {
    const drawing = state.activeDrawing();
    return rowsForDepthTrackCount(drawing, activeDepthTrackCount());
  });
  const linearDepthTrackCoverageLabel = (trackIndex) => {
    const covered = depthTrackCoverageCount(linearDepthRows(), trackIndex);
    const total = linearSeqs.length;
    return `${covered}/${total} record${total === 1 ? '' : 's'}`;
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
    { key: 'name', label: 'Name / Species' },
    { key: 'subtitle', label: 'Subtitle' },
    { key: 'replicon', label: 'Replicon', visibilityKey: 'linear_show_replicon', visibilityType: 'boolean' },
    { key: 'accession', label: 'Accession', visibilityKey: 'linear_accession_visibility', visibilityType: 'mode' },
    { key: 'length', label: 'Length / Coordinates', visibilityKey: 'linear_length_visibility', visibilityType: 'mode' }
  ]);
  const linearLabelHasSharedRow = computed(() => {
    const drawing = state.activeDrawing();
    return linearRecordLayoutHasSharedRow(
      linearSeqs,
      drawing.linearRecordRows,
      { enabled: Boolean(drawing.linearRecordLayoutEnabled.value) }
    );
  });
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
    const fields = linearLabelAutoFields.value.map((row) => row.label).join(' and ');
    const shown = resolveLinearLabelVisibility('auto', {
      hasSharedRow: linearLabelHasSharedRow.value
    });
    return shown
      ? `${fields}: Auto will show these fields throughout the diagram on the next successful Generate because no rendered row contains multiple records.`
      : `${fields}: Auto will hide these fields throughout the diagram on the next successful Generate because at least one rendered row contains multiple records. Choose Show in Record Labels to keep a field visible.`;
  });
  const focusLinearLabelVisibility = async (key) => {
    if (mode.value !== 'linear') return;
    const select = document.getElementById(`linear-label-visibility-${key}`);
    if (!select) return;
    // The select is rendered inside the Record Labels <details> of index.html.
    /** @type {HTMLDetailsElement} */ (select.closest('details')).open = true;
    await nextTick();
    select.scrollIntoView({ block: 'center' });
    select.focus({ preventScroll: true });
  };
  const legendPositionLabel = (position) => ({
    right: 'Right',
    left: 'Left',
    top: 'Top',
    bottom: 'Bottom',
    upper_left: 'Upper Left',
    upper_right: 'Upper Right',
    lower_left: 'Lower Left',
    lower_right: 'Lower Right',
    none: 'None'
  })[String(position || '').trim().toLowerCase()] || 'None';
  /** @param {DrawingState} drawing */
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
  const getDefinitionLineStyleSize = (kind) => {
    const drawing = state.activeDrawing();
    return optionalNumberInputValue(
      ensureDefinitionLineStyle(drawing, kind).font_size
    );
  };
  const setDefinitionLineStyleSize = (kind, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    setOptionalNumberInputValue(ensureDefinitionLineStyle(drawing, kind), 'font_size', value);
  };
  const getDefinitionLineStyleWeight = (kind) => {
    const drawing = state.activeDrawing();
    return ensureDefinitionLineStyle(drawing, kind).font_weight ?? '';
  };
  const setDefinitionLineStyleWeight = (kind, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const normalized = String(value || '').trim().toLowerCase();
    ensureDefinitionLineStyle(drawing, kind).font_weight = normalized === 'bold' ? 'bold' : null;
  };
  const getDefinitionLineStyleFill = (kind) => {
    const drawing = state.activeDrawing();
    return ensureDefinitionLineStyle(drawing, kind).fill ?? '';
  };
  const setDefinitionLineStyleColor = (kind, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const normalized = String(value || '').trim();
    ensureDefinitionLineStyle(drawing, kind).fill = normalized || null;
  };
  const getDefinitionLineStyleColorMode = (kind) => (
    colorValueMode(getDefinitionLineStyleFill(kind))
  );
  const getDefinitionLineStyleSwatchValue = (kind) => {
    return toNativeColorInputValue(getDefinitionLineStyleFill(kind));
  };
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
  const getDepthTrackLabel = (index) => {
    const config = readDepthTrackConfig(state.activeDrawing(), Math.max(0, Number(index) || 0));
    return String(config?.label ?? '');
  };
  const getDepthTrackColor = (index) => {
    const idx = Math.max(0, Number(index) || 0);
    const config = readDepthTrackConfig(state.activeDrawing(), idx);
    return String(config?.color || depthTrackFallbackColor(idx));
  };
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
  /** @param {DrawingState} drawing */
  const syncDepthTrackSlotLabelsForTrack = (drawing, index) => {
    depthTrackSlotCollections(drawing).forEach((slots) => {
      syncDepthSlotLabels(/** @type {any} */ ({
        slots,
        depthTracks: drawing.adv.depth_tracks,
        activeCount: drawing.adv.depth_tracks.length
      }));
    });
    void index;
  };
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
  const getDepthTrackLegendLabelForSlot = (slot) => (
    getDepthTrackLabel(normalizeDepthSlotTrackIndex(slot))
  );
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
  const depthTrackAutoLabels = [];
  /** @param {DrawingState} drawing */
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
  /** @param {DrawingState} drawing */
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
  const getCircularDepthFile = (index) => circularDepthRepresentatives()[Number(index)] || null;
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
  const getLinearDepthFile = (seq, index) => depthFileSlotsFromValue(seq?.depth)[Number(index)] || null;
  /** @param {DrawingState} drawing */
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
  const setLinearDepthFile = (seq, index, file) => {
    const drawing = state.drawings.linear;
    return setLinearDepthFiles(drawing, [seq], index, file);
  };
  const setLinearSourceDepthFile = (source, index, file) => {
    const drawing = state.drawings.linear;
    return setLinearDepthFiles(
      drawing,
      (source?.records || []).map(({ sequence }) => sequence),
      index,
      file
    );
  };
  const clearLinearSourceDepthFile = (source, index) => history.runUndoable(
    'Clear File Depth TSV',
    () => setLinearSourceDepthFile(source, index, null)
  );
  const addCircularDepthTrack = () => {
    const drawing = state.drawings.circular;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    depthTrackUiCounts.circular = sourceDepthTrackCount(files.c_depth, depthTrackUiCounts.circular) + 1;
    ensureDepthTrackConfigCount(drawing, depthTrackUiCounts.circular);
    if (hasCircularDepthFiles.value) drawing.form.show_depth = true;
  };
  const addLinearDepthTrack = () => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const nextCount = linearDepthTrackUiCount() + 1;
    padLinearDepthRows(nextCount);
    ensureDepthTrackConfigCount(drawing, nextCount);
    if (hasAnyLinearDepthFiles.value) drawing.form.show_depth = true;
  };
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
          .filter((slot) => isDefaultManagedDepthSlot(slot) && depthSlotTrackIndex(slot) === idx).length
        : 0;
      drawing.adv.circular_track_slots.splice(
        0,
        drawing.adv.circular_track_slots.length,
        ...reindexDepthSlots(/** @type {any} */ ({
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
        ? drawing.adv.linear_track_slots.reduce((count, slot, slotIndex) => {
            if (slotIndex >= previousAxisIndex || !isDefaultManagedDepthSlot(slot)) return count;
            return depthSlotTrackIndex(slot) === idx ? count + 1 : count;
          }, 0)
        : 0;
      drawing.adv.linear_track_slots.splice(
        0,
        drawing.adv.linear_track_slots.length,
        ...reindexDepthSlots(/** @type {any} */ ({
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
  /** @param {DrawingState} drawing */
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
    return (Array.isArray(drawing.circularConservation.series) ? drawing.circularConservation.series : []).map((entry, index) => ({
      index,
      filename: String(entry?.fileName || `source_${Number(index) + 1}`).trim(),
      sourceLabel: `${isCircularConservationUploadSource(drawing) ? 'BLAST' : 'Comparison'} ${Number(index) + 1}`,
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
    const step = history.runUndoable('Change uploaded file', settleComparisonRecordLabels);
    clearDerivedCircularConservationBlasts();
    files.c_conservation_fastas = [...normalizeFileList(files.c_conservation_fastas), selectedFile];
    losatCacheInfo.value = [];
    syncCircularConservationSeries();
    if (target) target.value = '';
    return step;
  };
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
      orderedBlasts.splice(idx, 1);
      orderedFastas.splice(idx, 1);
      orderedSequenceSources.splice(idx, 1);
      drawing.circularConservation.series.splice(idx, 1);
      drawing.circularConservation.series.forEach((seriesEntry, seriesIndex) => {
        seriesEntry.sourceIndex = seriesIndex;
      });
      drawing.circularConservation.labels = drawing.circularConservation.series
        .map((seriesEntry) => seriesEntry.label)
        .join(',');
      files.c_conservation_blasts = orderedBlasts;
      files.c_conservation_fastas = orderedFastas;
      files.c_conservation_sequence_sources = orderedSequenceSources;
      files.c_conservation_blasts_source = orderedBlasts.length > 0
        ? 'losat-cache'
        : null;
      drawing.circularConservation.enabled = orderedBlasts.length > 0;
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
        .filter((_, fileIndex) => fileIndex !== sourceIndex);
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
        drawing.circularConservation.series.map((entry) => `${entry?.sourceKey || ''}:${entry?.label || ''}:${entry?.color || ''}`).join('|')
      ];
    },
    ([slotsEnabled]) => {
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
      // It shows the Legend editor's edits laid out as Python lays the edited
      // rows out, before its entries are read (zero shift; OV-122, OV-124,
      // OV-126, OV-127).
      if (
        drawn && context.root.getAttribute(COMPOSITION_METADATA_ATTRIBUTE) !== null
        && legendActions.layOutMountedLegendEdits(context.root)
      ) {
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
        legendActions.captureOriginalStroke();
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
    afterReady() {
      previewFeatureSearch.handleMountedResultReady();
    }
  });
  const {
    runAnalysis: runGeneratedDiagramAnalysis,
    runCommittedCanonicalCandidate,
    projectCommittedRecordTransform,
    projectCommittedSimilarityAlignment,
    cancelRunAnalysis,
    runLabelReflow,
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
    settleComparisonRecordLabels,
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
    closeLabelTextScopeDialog: featureActions.closeLabelTextScopeDialog,
    clearLabelBuildNotices: featureActions.clearLabelBuildNotices,
    canonicalSessionVersion: SESSION_VERSION,
    adoptCanonicalRenderArtifacts,
    getCommittedCanonicalSession,
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
    }
  });
  const resolvePopupRotationFeature = ({ recordKey, biologicalFeatureId }) => {
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
      restoreIntentCheckpoint: (checkpoint) => recordDisplayControls.restoreTargetDraft(checkpoint),
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
  watch(clickedFeature, (popup) => {
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
  const cancelRecordActions = () => {
    featureRecordRotation.cancel();
    recordActionsExpanded.value = false;
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
  const resultsManager = createResultsManager({ state });

  const {
    waitForAuxiliaryFileImport, auxiliaryFileImportPending, canRetryAuxiliaryImportFailure, retryAuxiliaryImportFailure,
    resetModeTransientUi
  } = setupWatchers({
    state,
    rulePreparation,
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
      const ready = linearSeqs.every((/** @type {Record<string, any>} */ seq) => linearRecordSelector.statusFor(seq) === 'ready');
      if (ready) releaseSourceInputFailure('linear');
      return outcome;
    },
    resetPreviewViewport,
    resetRightDrawer: rightDrawerActions.resetRightDrawer,
    closeLabelTextScopeDialog: featureActions.closeLabelTextScopeDialog,
    clearLabelBuildNotices: featureActions.clearLabelBuildNotices,
    previewRuntime,
    preparePaletteDefinitions: paletteLoader.loadPaletteAsset
  });

  // Save and Load also wait for an edit still applying. Its owners exist only
  // here, so the root composes the availability that Save, Load, and their
  // controls read (R13).
  const sessionPreparationBusyReason = () => (
    history.mutationPending() || ruleMatchingPending.value || auxiliaryFileImportPending()
      ? 'Applying an edit. Retry after the edit finishes.'
      : ''
  );
  const sessionSaveLoadAvailability = (operation) => (
    sessionOperationAvailability(operation, sessionPreparationBusyReason)
  );
  const semanticMutationAvailable = computed(() => !sessionOperationAvailability());
  const sessionSaveAvailable = computed(() => !sessionSaveLoadAvailability('save'));
  const sessionLoadAvailable = computed(() => !sessionSaveLoadAvailability('load'));
  const sessionBusyReason = computed(() => sessionSaveLoadAvailability('save')?.reason || '');
  const circularRecordPresentationPanel = ref(null);
  let nextSessionPreviewToken = 1;
  const importSession = (event) => importSessionFromFile(event, {
    availability: sessionSaveLoadAvailability,
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
    afterImport: async (result) => {
      if (result?.status === 'ok' || result?.status === 'legacy') {
        historySnapshots.clearGeneratedArtifactIdentity({
          retainedBytes: result?.status === 'ok'
            ? Number(result.decompressedCharacters || 0) * 2
            : 0
        });
        // A replaced document resets transient UI: the selection named features
        // of the previous Session.
        featureSelection.clearFeatureSelection({ clearStatus: true, syncDom: false });
        await nextTick();
        recordSessionLifecycleEvent('history-baseline-start');
        if (!await history.initializeIntentBaseline('Loaded session', { isCurrent: result.isCurrent })) {
          throw new Error('Session loading canceled.');
        }
        featureActions.clearSpecificRulePatternDrafts();
        annotationImportNotice.value = '';
        specificRuleNotice.value = '';
        recordSessionLifecycleEvent('history-baseline-end');
        if (circularRecordPresentationPanel.value) circularRecordPresentationPanel.value.open = false;
        closeLegendStrokeOptions();
      }
    }
  });

  const {
    addNewLegendEntry,
    updateLegendEntryColor,
    deleteLegendEntry,
    moveLegendEntryUp,
    moveLegendEntryDown,
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
    reconcileLegendEntries,
    reconcileStrokeOverrides,
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
    handleColorScopeChoice,
    handleFeatureStyleScopeChoice,
    handleLegendNameCommit,
    handleLegendRenameChoice,
    selectLegendNameOption,
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
  // Feature Edits TSV (D-07, R3). History restores the mounted Legend
  // inventory; a newly displayed Result receives the diagram-wide Legend
  // operations Generate applies. A loaded table (`reflow`) also places the
  // labels, as a visibility edit does. The palette and the rules (`colors`)
  // project through `projectPaletteAndRules`.
  /**
   * @param {{
   *   colors?: boolean, prepareRules?: boolean, visibility?: boolean, rerender?: boolean,
   *   reflow?: boolean, labels?: boolean,
   *   legend?: Parameters<typeof reconcileLegendEntries>[0] | null,
   *   strokes?: Parameters<typeof reconcileStrokeOverrides>[0] | null
   * }} [options]
   */
  const projectMountedEditorIntent = async ({
    colors = false,
    prepareRules = colors,
    visibility = false,
    rerender = false,
    reflow = false,
    legend = null,
    strokes = null,
    labels = false
  } = {}) => {
    if (colors && !await projectPaletteAndRules({ prepareRules })) return false;
    if (visibility) await projectFeatureVisibility({ rerender, reflow });
    if (legend) reconcileLegendEntries(legend);
    if (strokes) reconcileStrokeOverrides(strokes);
    if (labels) reconcileLabelOverrides();
    return true;
  };

  /**
   * @param {string} label
   * @param {number[] | null} indexes
   */
  const restoreLegendItems = (label, indexes) => history.runUndoableCheckpoint(label, async () => {
    const restored = await restoreDeletedLegendEntries(indexes);
    if (restored === true) await projectMountedEditorIntent({ colors: true, prepareRules: false });
    return restored;
  });

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
    // The Legend side of this step: the restored entry owners and `from`, the
    // list the step leaves (the restored list when the step kept it). The
    // Legend owner decides what it describes (B19).
    const legendChange = (Array.isArray(changes) ? changes : []).find(({ path } = {}) => (
      path?.length === 3 && path[0] === 'editorState' && path[1] === 'legend' && path[2] === 'entries'
    ));
    const projected = await projectMountedEditorIntent({
      colors,
      prepareRules: changedDomains.has('features') || rulesChanged || !rulePreparation.isPrepared(),
      visibility: changedDomains.has('features'),
      rerender: true,
      legend: editorState
        ? {
            entryOwners: _intent.editorState.legend.entryOwners,
            from: legendChange
              ? legendChange[direction === 'undo' ? 'after' : 'before']
              : _intent.editorState.legend.entries
          }
        : null,
      strokes: editorState ? { changes } : null,
      labels: changedDomains.has('features') || editorState
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

  // Each Result's bytes reflect the editor state it was committed or last
  // shown with. A displayed Result receives a domain only when that state
  // changed since, so a Result without new edits gets no projection work and
  // an Undo reaches a Result that is displayed again.
  const projectedEditorStateByResult = new Map();
  let lastBoundResultIdentity = '';
  // A displayed Result whose label intent changed since it was last shown
  // receives the label projection in the binder's label step, also when no
  // label intent remains (an undone or replaced Label TSV import).
  let labelProjectionResultIdentity = '';
  /** @param {DrawingState} drawing */
  const currentEditorProjectionState = (drawing) => ({
    colors: [
      toRaw(appliedPaletteColors.value),
      JSON.stringify([drawing.manualSpecificRules, drawing.featureColorOverrides, drawing.legendColorOverrides])
    ],
    visibility: JSON.stringify([
      Object.values(drawing.featureOverrides).map((row) => [row.recordKey, row.biologicalFeatureId, row.featureVisibility]),
      drawing.featureVisibilityManualRules
    ]),
    labels: JSON.stringify([
      Object.values(drawing.featureOverrides).map((row) => [
        row.recordKey, row.biologicalFeatureId, row.labelVisibility, row.labelText, row.labelSourceText
      ]),
      drawing.labelTextBulkOverrides
    ]),
    // An edited Legend order, or '' for the default order (D-08).
    legendOrder: isLegendOrderEdited(drawing.legendEntries.value, originalLegendOrder.value)
      ? JSON.stringify(drawing.legendEntries.value.map((entry) => entry.caption))
      : ''
  });
  const sameColors = (left, right) => left[0] === right[0] && left[1] === right[1];
  /** @param {DrawingState} drawing */
  const rememberCommittedEditorState = (drawing, context) => {
    const current = currentEditorProjectionState(drawing);
    const identities = new Set(liveResultIdentities().filter(Boolean));
    identities.forEach((identity) => {
      if (!projectedEditorStateByResult.has(identity)) projectedEditorStateByResult.set(identity, current);
    });
    [...projectedEditorStateByResult.keys()].forEach((identity) => {
      if (!identities.has(identity)) projectedEditorStateByResult.delete(identity);
    });
    [...departedResultIntent.keys()].forEach((identity) => {
      if (!identities.has(identity)) departedResultIntent.delete(identity);
    });
    lastBoundResultIdentity = context.resultIdentity;
  };
  // E1: the Result shown until a mode switch followed every live edit. When
  // it is shown again with the same editor intent, nothing is projected and
  // its bytes stay; otherwise the edits made since are projected onto it.
  /** @type {Map<string, string>} */
  const departedResultIntent = new Map();
  /** @param {DrawingState} drawing */
  const displayedIntentSignature = (drawing) => JSON.stringify([
    currentEditorProjectionState(drawing),
    drawing.featureStrokeOverrides,
    drawing.legendStrokeOverrides,
    drawing.deletedLegendEntries.value.map((entry) => entry.originalCaption || entry.caption),
    [...(drawing.addedLegendCaptions.value || [])],
    drawing.legendEntries.value.filter((entry) => entry.caption !== entry.originalCaption)
      .map((entry) => [entry.originalCaption, entry.caption]),
    // The direct additions, which no Result's inventory lists.
    drawing.legendEntries.value.filter((entry) => !originalLegendOrder.value.includes(entry.originalCaption || entry.caption))
      .map((entry) => entry.caption)
  ]);
  /** @param {DrawingState} drawing */
  const rememberDepartingResultProjection = (drawing) => {
    if (lastBoundResultIdentity && projectedEditorStateByResult.has(lastBoundResultIdentity)) {
      projectedEditorStateByResult.set(lastBoundResultIdentity, currentEditorProjectionState(drawing));
      departedResultIntent.set(lastBoundResultIdentity, displayedIntentSignature(drawing));
    }
    lastBoundResultIdentity = '';
  };
  /**
   * @param {DrawingState} drawing
   * @param {number} resultIndex
   * @param {{ replayDefaultLegendOrder?: string[] | null }} [options]
   */
  const compileDisplayedResultOperations = (drawing, resultIndex, { replayDefaultLegendOrder = null } = {}) => {
    const catalog = toRaw(state.featureCatalog.value);
    if (!catalog) return null;
    const plan = compileDirectEditorMutationPlan({
      catalogAdmission: admitFeatureCatalog(catalog, toRaw(results.value), { mode: state.generatedMode.value }),
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
      replayDefaultLegendOrder
    });
    return plan.operationsByResult[resultIndex] || null;
  };
  const DISPLAY_PROJECTED_DOMAINS = Object.freeze([
    'featureFills', 'featureStrokes', 'featureVisibility',
    'legendFills', 'legendStrokes', 'legendRenames', 'legendDeletes', 'legendAdds', 'legendOrder'
  ]);
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
    const departedIntent = departedResultIntent.get(identity);
    departedResultIntent.delete(identity);
    if (departedIntent !== undefined && departedIntent === displayedIntentSignature(drawing)) {
      lastBoundResultIdentity = identity;
      labelProjectionResultIdentity = '';
      recordStructuralMetric('displayedResultEditorProjectionCount', 0, {
        phase: context.phase,
        rootGeneration: context.rootGeneration
      });
      return;
    }
    // The Result shown until now followed every live edit.
    if (lastBoundResultIdentity && lastBoundResultIdentity !== identity
      && projectedEditorStateByResult.has(lastBoundResultIdentity)) {
      projectedEditorStateByResult.set(lastBoundResultIdentity, current);
    }
    lastBoundResultIdentity = identity;
    const previous = projectedEditorStateByResult.get(identity) || current;
    const colors = !sameColors(previous.colors, current.colors);
    const visibility = previous.visibility !== current.visibility;
    labelProjectionResultIdentity = previous.labels !== current.labels ? identity : '';
    // B20: a Result last shown with another Legend order receives the current
    // order, also the default order, which is its own generated order (OV-47).
    const replayDefaultLegendOrder = previous.legendOrder !== current.legendOrder
      ? resultLegendOrder : null;
    /** @type {ReturnType<typeof compileDisplayedResultOperations>} */
    let operations = null;
    try {
      operations = compileDisplayedResultOperations(drawing, context.resultIndex, { replayDefaultLegendOrder });
    } catch (error) {
      console.error('Editor edits could not be compiled for the displayed Result.', normalizeUserFacingError(error));
    }
    const hasOperations = Boolean(operations)
      && DISPLAY_PROJECTED_DOMAINS.some((domain) => operations[domain].length > 0);
    const legend = {
      resultIdentity: identity,
      liveResultIdentities: liveResultIdentities(),
      deletedCaptions: (operations?.legendDeletes || []).map(({ caption }) => caption)
    };
    const restoresLegend = legendActions.hasRetiredResultLegend(legend);
    const projects = colors || visibility || hasOperations || restoresLegend;
    recordStructuralMetric('displayedResultEditorProjectionCount', projects ? 1 : 0, {
      phase: context.phase,
      rootGeneration: context.rootGeneration
    });
    if (!projects) return;
    try {
      await projectMountedEditorIntent({ colors, visibility });
      const legendChanged = legendActions.prepareDisplayedResultLegend(context.root, legend);
      previewRuntime.applyEditorOperations(hasOperations ? operations : null, {
        afterApply: () => { if (legendChanged) legendActions.onLegendGeometryChanged(); }
      });
      projectedEditorStateByResult.set(identity, current);
    } catch (error) {
      console.error('Editor edits could not be shown on the displayed Result.', normalizeUserFacingError(error));
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
  //   4. Reset the departing mode's transient UI.
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
    // departing mode's Result.
    resetModeTransientUi();
    failedGeneratePreservedResult.value = false;
    // 5. Presentation: the arriving Result is shown as a selection.
    previewRuntime.presentSelectedResult({ modeArrival: true });
    return true;
  };
  historySnapshots.registerModeTransition(transitionDiagramMode);
  // The committed Session of a mode's artifact, shown or kept in its slot.
  historySnapshots.registerModeCommittedSession((targetMode) => (
    targetMode === mode.value
      ? getCommittedCanonicalSession()
      : artifactSlots[targetMode]?.runtimeState?.canonical?.committedCanonicalSession ?? null
  ));
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

  const { updatePalette, resetColors } = resultsManager;
  const undoableAction = (label, fn) => (...args) => history.runUndoable(label, () => fn(...args));
  // Python's report of the per-feature edits the Results do not draw (design
  // Q4 3.4): dormant edits outside the crop or display, edits whose feature
  // the source does not have, and the edits a source-replacing Generate removed.
  const featureIdentityNoticeSummary = computed(() => {
    const drawing = state.activeDrawing();
    const notices = Array.isArray(featureIdentityNotices.value) ? featureIdentityNotices.value : [];
    const summary = {
      dormant: notices.filter((notice) => notice.status !== 'unresolved').length,
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
  const removeUnmatchedFeatureEdits = undoableAction('Remove unmatched feature edits', () => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    return removeUnresolvedFeatureEdits({
      featureOverrides: drawing.featureOverrides,
      featurePlacementOverrides: drawing.featurePlacementOverrides,
      notices: featureIdentityNotices.value
    }) > 0;
  });
  const addFeatureVisibilityRuleWithHistory = undoableAction('Add feature visibility rule', addFeatureVisibilityRule);
  const moveFeatureVisibilityRuleDownWithHistory = undoableAction('Move feature visibility rule', moveFeatureVisibilityRuleDown);
  const moveFeatureVisibilityRuleUpWithHistory = undoableAction('Move feature visibility rule', moveFeatureVisibilityRuleUp);
  const removeFeatureVisibilityRuleWithHistory = undoableAction('Remove feature visibility rule', removeFeatureVisibilityRule);
  const setFeatureVisibilityRuleFieldWithHistory = undoableAction(
    'Edit feature visibility rule',
    setFeatureVisibilityRuleField
  );
  const setFeatureVisibilityWithHistory = undoableAction('Change feature visibility', setFeatureVisibility);
  const updateClickedFeatureVisibilityWithHistory = undoableAction(
    'Change feature visibility',
    updateClickedFeatureVisibility
  );
  const handleFeatureVisibilityScopeChoiceWithHistory = undoableAction(
    'Change feature visibility',
    handleFeatureVisibilityScopeChoice
  );
  const requestFeatureColorChangeWithHistory = undoableAction('Change feature color', requestFeatureColorChange);
  const setFeatureColorValueWithHistory = undoableAction('Change feature color', setFeatureColorValue);
  const updateClickedFeatureColorWithHistory = undoableAction('Change feature color', updateClickedFeatureColor);
  const setLegendEntryStrokeColorValueWithHistory = undoableAction(
    'Change legend stroke color',
    setLegendEntryStrokeColorValue
  );
  // OV-161: Cancel closes a scope dialog at once and records no History step;
  // only a real choice is one undoable step, as with askLabelOn and
  // handleLabelOnChoice. Cancel used to wait for the step's intent capture,
  // which after a Session load took seconds.
  /**
   * @param {() => string} label
   * @param {(choice: string, ...rest: any[]) => any} handler
   */
  const scopeChoiceWithHistory = (label, handler) => (/** @type {string} */ choice, /** @type {any[]} */ ...rest) => (
    choice === 'cancel'
      ? cancelFeatureStyleScope()
      : history.runUndoable(label(), () => handler(choice, ...rest))
  );
  const handleColorScopeChoiceWithHistory = scopeChoiceWithHistory(() => 'Change feature color', handleColorScopeChoice);
  const handleFeatureStyleScopeChoiceWithHistory = scopeChoiceWithHistory(
    () => (featureStyleScopeDialog.kind === 'stroke' ? 'Change feature stroke' : 'Change feature color'),
    handleFeatureStyleScopeChoice
  );
  const handleLegendNameCommitWithHistory = undoableAction('Rename legend item', handleLegendNameCommit);
  const handleLegendRenameChoiceWithHistory = undoableAction('Rename legend item', handleLegendRenameChoice);
  const handleResetColorChoiceWithHistory = undoableAction('Reset feature color', handleResetColorChoice);
  const resetClickedFeatureFillColorWithHistory = undoableAction('Reset feature color', resetClickedFeatureFillColor);
  const updateClickedFeatureStrokeWithHistory = undoableAction('Change feature stroke', updateClickedFeatureStroke);
  const setClickedFeatureStrokeColorValueWithHistory = undoableAction(
    'Change feature stroke',
    setClickedFeatureStrokeColorValue
  );
  const setClickedFeatureStrokeWidthValueWithHistory = undoableAction(
    'Change feature stroke',
    setClickedFeatureStrokeWidthValue
  );
  const resetClickedFeatureStrokeWithHistory = undoableAction('Reset feature stroke', resetClickedFeatureStroke);
  const setFeatureColorWithHistory = undoableAction('Change feature color', setFeatureColor);
  const selectedFeatureBulkColor = ref('#2563eb');
  const selectedFeatureBulkCaption = ref('Selected features');
  const selectedFeatureBulkVisibility = ref('off');
  const selectedFeatureBulkStrokeColor = ref('#1f2937');
  const selectedFeatureBulkStrokeWidth = ref(1.5);
  const applySelectedFeatureColor = () => history.runUndoable('Change selected feature color', async () => {
    const changed = await applyColorToSelectedFeatures(
      selectedFeatures.value,
      selectedFeatureBulkColor.value,
      selectedFeatureBulkCaption.value
    );
    if (changed) featureSelection.syncFeatureSelectionClasses();
    return changed;
  });
  const applySelectedFeatureVisibility = async () => {
    const changed = await history.runUndoableCommand('Change selected feature visibility', () =>
      buildSelectedFeaturesVisibilityCommand(selectedFeatures.value, selectedFeatureBulkVisibility.value)
    );
    if (changed) featureSelection.clearFeatureSelection({ clearStatus: true });
    return changed;
  };
  const applySelectedFeatureStroke = () => history.runUndoable('Change selected feature stroke', async () => {
    const changed = applyStrokeToSelectedFeatures(
      selectedFeatures.value,
      selectedFeatureBulkStrokeColor.value,
      selectedFeatureBulkStrokeWidth.value
    );
    if (changed) featureSelection.syncFeatureSelectionClasses();
    return changed;
  });
  const openFirstSelectedFeature = (event = null) => {
    const first = selectedFeatures.value[0] || null;
    if (!first) return null;
    return openFeatureEditorForFeature(first, event);
  };
  const updateClickedFeatureLabelTextWithHistory = undoableAction('Change label text', updateClickedFeatureLabelText);
  const handleLabelTextScopeChoiceWithHistory = undoableAction('Change label text', handleLabelTextScopeChoice);
  const handleHiddenLabelTextChoiceWithHistory = undoableAction('Change label visibility', handleHiddenLabelTextChoice);
  const requestLabelTextChangeByFeatureIdWithHistory = undoableAction(
    'Change label text',
    requestLabelTextChangeByFeatureId
  );
  const requestLabelTextChangeByKeyWithHistory = undoableAction('Change label text', requestLabelTextChangeByKey);
  const resetAllLabelTextOverridesWithHistory = undoableAction('Reset label edits', resetAllLabelTextOverrides);
  const loadLabelOverrideTableWithHistory = undoableAction('Load label edits', loadLabelOverrideTable);
  const runInfoCopyStatus = ref('');
  const exactReplayCopyStatus = ref('');

  const runInfoElapsedText = (info) => formatElapsedMs(info?.elapsedMs);
  const runInfoReproducibilityText = (info) => reproducibilityLabel(info?.reproducibility?.level);
  const runInfoHasCliHelperFiles = computed(() =>
    Array.isArray(lastRunInfo.value?.helperFiles) && lastRunInfo.value.helperFiles.length > 0
  );
  const copyRunInfoCommand = async (commandValue, status) => {
    const command = String(commandValue || '');
    if (!command) return;
    try {
      await copyTextToClipboard(command);
      status.value = 'Copied';
      setTimeout(() => {
        if (status.value === 'Copied') status.value = '';
      }, 1600);
    } catch (error) {
      console.warn('Failed to copy the requested command.', normalizeUserFacingError(error));
      status.value = 'Copy failed';
      setTimeout(() => {
        if (status.value === 'Copy failed') status.value = '';
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
    const hasRegionAnnotations = drawing.annotationSets.some((set) => (
      Array.isArray(set?.annotations) && set.annotations.length > 0
    ));
    if (!hasAutomaticSequence && !hasRegionAnnotations) {
      return { catalog: null, error: '' };
    }
    if (privateCandidate) {
      const sources = await Promise.all(linearSeqs.map(async (seq) => {
        const inputType = lInputType.value;
        const primaryFile = inputType === 'gff' ? seq.gff : seq.gb;
        const pairedFile = inputType === 'gff' ? seq.fasta : null;
        const sourceKey = annotationSourceKey({ scope: 'linear', uid: seq.uid, inputType, primaryFile, pairedFile });
        try {
          const records = inputType === 'gff'
            ? await discoverGffFastaRecords({ gffFile: primaryFile, fastaFile: pairedFile })
            : await discoverSequenceRecords({ file: primaryFile, format: 'genbank' });
          return { sourceKey, selector: seq.region_record_id, hasInput: Boolean(primaryFile), status: 'ready', records };
        } catch (error) {
          return { sourceKey, selector: seq.region_record_id, hasInput: Boolean(primaryFile), status: 'error', error: error.message, records: [] };
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
          error: linearSeqs.map(seq => linearRecordSelector.errorModelFor(seq)).find(error => error?.code)
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
          errorLog.value = normalizeUserFacingError(comparisonExecution.message, { operation: 'generate', stage: 'request-validation' });
          generationFailureRecovery.value = results.value.length ? 'preserved' : 'no-result';
          failedGeneratePreservedResult.value = results.value.length > 0;
          if (mode.value === 'linear') await focusLinearComparisonIssue();
          return { status: 'error', error: errorLog.value };
        }
        return { status: 'ready', comparisonPlanSnapshot, comparisonExecution };
      },
      afterGenerate: async (result) => {
        if (result?.status === 'error' && mode.value === 'linear') {
          await focusLinearComparisonIssue();
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

  /** @param {DrawingState} drawing */
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
    getRecordCatalog: getAnnotationRecordCatalog,
    getCommittedRequest: getCommittedCanonicalRenderRequest,
    getCommittedSession: getCommittedCanonicalSession,
    projectCommittedAlignment: /** @type {any} */ (projectCommittedSimilarityAlignment),
    // R13: the root runs the candidate with the record display's orientation
    // checkpoint; the alignment owner holds neither owner.
    runRecordAlignment: ({ orientations, ...run }) => runCommittedCanonicalCandidate({ ...run,
      captureIntentCheckpoint: () => recordDisplayControls.captureAlignmentOrientationIntent(orientations),
      restoreIntentCheckpoint: (checkpoint) => recordDisplayControls.restoreAlignmentOrientationIntent(checkpoint),
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
  const similarityAlignmentCanvasHover = ref(null);
  similarityAlignmentPorts.refreshCanvas = () => {
    if (!similarityAlignmentActions.dialogOpen.value
      || similarityAlignmentActions.status.value !== 'reviewing'
      || !similarityAlignmentActions.isDraftArtifactCurrent()) {
      featureActions.clearAlignmentOverlay();
      similarityAlignmentCanvasHover.value = null;
      return;
    }
    const current = similarityAlignmentActions.draft.value;
    similarityAlignmentCanvasHover.value = null;
    featureActions.showAlignmentOverlay({
      reference: current.response.reference,
      ambiguities: current.rows.filter(({ candidates }) => candidates.length > 0),
      onSelect: similarityAlignmentActions.selectCandidate,
      onHover: (recordKey, candidateKey) => {
        similarityAlignmentCanvasHover.value = recordKey ? { recordKey, candidateKey } : null;
      }
    });
  };
  watch([
    similarityAlignmentActions.dialogOpen,
    () => similarityAlignmentActions.draft.value?.rows,
    similarityAlignmentActions.status
  ], similarityAlignmentPorts.refreshCanvas, { flush: 'post' });
  /** @type {HTMLElement | null} */
  let similarityAlignmentReturnFocus = null;
  const similarityAlignmentPaletteRef = ref(null);
  /** @type {{ x: number | null, y: number | null }} */
  const similarityAlignmentPalettePosition = reactive({ x: null, y: null });
  const similarityAlignmentCompact = ref(false);
  const syncSimilarityAlignmentCompact = () => {
    const preview = document.querySelector('[aria-label="Result Preview"]');
    similarityAlignmentCompact.value = Boolean(preview
      && getComputedStyle(preview).getPropertyValue('--alignment-review-compact').trim() === '1');
  };
  const similarityAlignmentEditorDisabledReason = computed(() => (
    similarityAlignmentCompact.value && similarityAlignmentActions.dialogOpen.value
      ? 'Finish or cancel alignment review before opening Editor.' : ''
  ));
  similarityAlignmentPorts.reviewBlocksEditor = () => Boolean(similarityAlignmentEditorDisabledReason.value);
  /** @type {ResizeObserver | null} */
  let similarityAlignmentPreviewObserver = null;
  watch([similarityAlignmentActions.dialogOpen, similarityAlignmentCompact], ([open, compact]) => {
    if (open && compact) {
      stopSimilarityAlignmentPaletteDrag();
      rightDrawerActions.closeRightDrawer();
    }
  }, { flush: 'sync' });
  /** @type {{ x: number, y: number } | null} */
  let similarityAlignmentPaletteDrag = null;
  const clampSimilarityAlignmentPalette = () => {
    const palette = similarityAlignmentPaletteRef.value;
    if (!palette || similarityAlignmentCompact.value) return;
    const margin = 12;
    const maxX = Math.max(margin, window.innerWidth - palette.offsetWidth - margin);
    const maxY = Math.max(margin, window.innerHeight - palette.offsetHeight - margin);
    similarityAlignmentPalettePosition.x = Math.min(
      Math.max(similarityAlignmentPalettePosition.x ?? maxX, margin), maxX
    );
    similarityAlignmentPalettePosition.y = Math.min(
      Math.max(similarityAlignmentPalettePosition.y ?? margin, margin), maxY
    );
  };
  const similarityAlignmentPaletteStyle = computed(() => similarityAlignmentCompact.value ? {} : ({
    left: similarityAlignmentPalettePosition.x === null ? undefined : `${similarityAlignmentPalettePosition.x}px`,
    right: similarityAlignmentPalettePosition.x === null ? undefined : 'auto',
    top: similarityAlignmentPalettePosition.y === null ? undefined : `${similarityAlignmentPalettePosition.y}px`
  }));
  const moveSimilarityAlignmentPalette = (event) => {
    if (!similarityAlignmentPaletteDrag) return;
    similarityAlignmentPalettePosition.x = event.clientX - similarityAlignmentPaletteDrag.x;
    similarityAlignmentPalettePosition.y = event.clientY - similarityAlignmentPaletteDrag.y;
    clampSimilarityAlignmentPalette();
    event.preventDefault();
  };
  const stopSimilarityAlignmentPaletteDrag = () => {
    similarityAlignmentPaletteDrag = null;
    document.removeEventListener('pointermove', moveSimilarityAlignmentPalette);
    document.removeEventListener('pointerup', stopSimilarityAlignmentPaletteDrag);
    document.removeEventListener('pointercancel', stopSimilarityAlignmentPaletteDrag);
  };
  const startSimilarityAlignmentPaletteDrag = (event) => {
    if (similarityAlignmentCompact.value || event.button !== 0
      || event.target.closest('button, a, input, select, textarea')) return;
    const rect = similarityAlignmentPaletteRef.value?.getBoundingClientRect();
    if (!rect) return;
    similarityAlignmentPaletteDrag = { x: event.clientX - rect.left, y: event.clientY - rect.top };
    document.addEventListener('pointermove', moveSimilarityAlignmentPalette);
    document.addEventListener('pointerup', stopSimilarityAlignmentPaletteDrag);
    document.addEventListener('pointercancel', stopSimilarityAlignmentPaletteDrag);
    event.preventDefault();
    event.stopPropagation();
  };
  const rememberSimilarityAlignmentInvoker = (event) => {
    similarityAlignmentReturnFocus = event?.currentTarget || document.activeElement || null;
  };
  const restoreSimilarityAlignmentFocus = async () => {
    await nextTick();
    const target = similarityAlignmentReturnFocus;
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
    clampSimilarityAlignmentPalette();
    document.getElementById('similarity-alignment-dialog-title')?.focus();
  };
  const handleSimilarityAlignmentEscape = (event) => {
    if (event.key !== 'Escape' || !similarityAlignmentActions.dialogOpen.value) return;
    if (document.querySelector('[data-linear-source-removal-dialog]')
      || event.target?.closest?.('.fixed.inset-0')) return;
    event.preventDefault();
    event.stopImmediatePropagation();
    void cancelSimilarityAlignmentDialog();
  };
  const openSimilarityAlignmentReset = (event) => {
    rememberSimilarityAlignmentInvoker(event);
    similarityAlignmentActions.openReset();
  };
  watch([similarityAlignmentActions.dialogOpen, similarityAlignmentActions.resetDialogOpen], async ([open, resetOpen]) => {
    if (resetOpen) {
      await nextTick();
      document.getElementById('similarity-alignment-reset-title')?.focus();
      return;
    }
    if (open) return;
    stopSimilarityAlignmentPaletteDrag();
    similarityAlignmentPalettePosition.x = null;
    similarityAlignmentPalettePosition.y = null;
    void restoreSimilarityAlignmentFocus();
  });
  onMounted(() => {
    similarityAlignmentPreviewObserver = new ResizeObserver(syncSimilarityAlignmentCompact);
    const pane = document.querySelector('.result-pane');
    if (pane) similarityAlignmentPreviewObserver.observe(pane);
    document.addEventListener('keydown', handleSimilarityAlignmentEscape, true);
    window.addEventListener('resize', clampSimilarityAlignmentPalette);
  });
  onUnmounted(() => {
    similarityAlignmentPreviewObserver?.disconnect();
    stopSimilarityAlignmentPaletteDrag();
    document.removeEventListener('keydown', handleSimilarityAlignmentEscape, true);
    window.removeEventListener('resize', clampSimilarityAlignmentPalette);
  });
  watch(() => results.value.length, syncSimilarityAlignmentCompact, { flush: 'post' });
  const finishSimilarityAlignmentStart = async (outcome) => {
    if (similarityAlignmentActions.dialogOpen.value) {
      if (similarityAlignmentActions.error.value) {
        await nextTick();
        /** @type {HTMLElement | null} */ (document.querySelector('[data-similarity-alignment-error]'))?.focus();
      } else await focusSimilarityAlignmentDialog();
    } else await restoreSimilarityAlignmentFocus();
    return outcome;
  };
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
    const outcome = await (reset ? similarityAlignmentActions.applyReset() : similarityAlignmentActions.applyDraft());
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
    return {
      id: orthogroupId,
      displayName: orthogroupActions.resolveOrthogroupName(group),
      description: orthogroupActions.resolveOrthogroupDescription(group),
      scope: orthogroupActions.orthogroupScope(group),
      scopeLabel: orthogroupActions.orthogroupScopeLabel(group),
      candidates: Array.isArray(group.nameCandidates) ? group.nameCandidates : [],
      memberCount: Number(group.member_count || members.length || 0),
      recordCoverage: Number(group.record_coverage_count || membersByRecord.length || 0),
      ntSequenceCount: orthogroupActions.getOrthogroupSequenceCount(group, 'nt'),
      aaSequenceCount: orthogroupActions.getOrthogroupSequenceCount(group, 'aa'),
      currentMember,
      membersByRecord
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

  const openOrthogroupInDrawer = (orthogroupId) => {
    if (!orthogroupActions.selectOrthogroup(orthogroupId)) return false;
    rightDrawerActions.openRightDrawerTab('orthogroups');
    return true;
  };

  const reselectSimilarityAlignmentReference = () => {
    const groupId = similarityAlignmentActions.repair.value?.groupId
      || state.similarityAlignmentPlan.value?.groupId
      || '';
    similarityAlignmentActions.drawerReferenceKey.value = '';
    return openOrthogroupInDrawer(groupId);
  };

  const openClickedOrthogroupInEditor = () => {
    const orthogroupId = String(clickedFeature.value?.orthogroupId || '').trim();
    if (!openOrthogroupInDrawer(orthogroupId)) return;
    clickedFeature.value = null;
  };

  const { resetAllPositions, resetCanvasPadding } = legendLayout;

  // Reset Settings resets both drawings (PD-OI-070); the Linear records and
  // comparisons are the Linear drawing's.
  const resetSettings = () => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const proceed = window.confirm(
      'Reset all settings to the webapp defaults?\n\nUploaded files and current results will be kept.'
    );
    if (!proceed) return false;

    return history.runUndoableCheckpoint('Reset settings', async () => {
      featureActions.clearSpecificRulePatternDrafts();
      resetSettingsState(state);
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
      })), { alignmentMutation: 'settings reset.' });
      invalidateLinearComparisonArtifacts(drawing);
      matchSequenceRegistry?.reset?.();
      circularTrackNewRenderer.value = 'dinucleotide_skew';
      linearTrackNewRenderer.value = 'dinucleotide_skew';
      depthTrackUiCounts.circular = 1;
      ensureDepthTrackConfigCount(state.drawings.circular, activeDepthTrackCount('circular'));
      ensureDepthTrackConfigCount(drawing, activeDepthTrackCount('linear'));
      return true;
    });
  };

  const resetLayout = () => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    resetAllPositions();
    resetCanvasPadding();
    resetLayoutState(state);
    resetPreviewViewport({ resetZoom: true });
  };

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

  const clampNumber = (value, min, max) => {
    const safeMin = Number.isFinite(min) ? min : 0;
    const safeMax = Number.isFinite(max) ? Math.max(safeMin, max) : safeMin;
    const numeric = Number(value);
    if (!Number.isFinite(numeric)) return safeMin;
    return Math.min(Math.max(numeric, safeMin), safeMax);
  };

  // The popup's size follows the app-level rich popup preference.
  const getFeaturePopupConstraints = (left = clickedFeaturePos.x, top = clickedFeaturePos.y) => {
    const viewportWidth = Math.max(1, window.innerWidth || 1);
    const viewportHeight = Math.max(1, window.innerHeight || 1);
    const availableWidth = Math.max(1, viewportWidth - (FEATURE_POPUP_MARGIN * 2));
    const availableHeight = Math.max(1, viewportHeight - (FEATURE_POPUP_MARGIN * 2));
    const desiredMinWidth =
      richFeaturePopup.value === false ? FEATURE_POPUP_SIMPLE_MIN_WIDTH : FEATURE_POPUP_RICH_MIN_WIDTH;
    const minWidth = Math.min(desiredMinWidth, availableWidth);
    const minHeight = Math.min(FEATURE_POPUP_MIN_HEIGHT, availableHeight);
    return {
      minWidth,
      minHeight,
      maxWidth: Math.max(minWidth, viewportWidth - left - FEATURE_POPUP_MARGIN),
      maxHeight: Math.max(minHeight, viewportHeight - top - FEATURE_POPUP_MARGIN)
    };
  };

  const featurePopupStyle = computed(() => {
    const style = {
      top: `${clickedFeaturePos.y}px`,
      left: `${clickedFeaturePos.x}px`,
      maxHeight: `${getFeaturePopupConstraints().maxHeight}px`
    };
    if (featurePopupSize.width > 0) {
      style.width = `${featurePopupSize.width}px`;
    }
    if (featurePopupSize.height > 0) {
      style.height = `${featurePopupSize.height}px`;
    }
    return style;
  });

  const getPairwiseMatchPopupConstraints = (left = clickedPairwiseMatchPos.x, top = clickedPairwiseMatchPos.y) => {
    const viewportWidth = Math.max(1, window.innerWidth || 1);
    const viewportHeight = Math.max(1, window.innerHeight || 1);
    const availableWidth = Math.max(1, viewportWidth - (FEATURE_POPUP_MARGIN * 2));
    const availableHeight = Math.max(1, viewportHeight - (FEATURE_POPUP_MARGIN * 2));
    const minWidth = Math.min(PAIRWISE_MATCH_POPUP_MIN_WIDTH, availableWidth);
    const minHeight = Math.min(PAIRWISE_MATCH_POPUP_MIN_HEIGHT, availableHeight);
    return {
      minWidth,
      minHeight,
      maxWidth: Math.max(minWidth, viewportWidth - left - FEATURE_POPUP_MARGIN),
      maxHeight: Math.max(minHeight, viewportHeight - top - FEATURE_POPUP_MARGIN)
    };
  };

  const pairwiseMatchPopupStyle = computed(() => {
    const style = {
      top: `${clickedPairwiseMatchPos.y}px`,
      left: `${clickedPairwiseMatchPos.x}px`,
      maxHeight: `${getPairwiseMatchPopupConstraints().maxHeight}px`
    };
    if (pairwiseMatchPopupSize.width > 0) {
      style.width = `${pairwiseMatchPopupSize.width}px`;
    }
    if (pairwiseMatchPopupSize.height > 0) {
      style.height = `${pairwiseMatchPopupSize.height}px`;
    }
    return style;
  });

  const pairwiseBlockOrthogroups = computed(() => (
    Array.isArray(clickedPairwiseMatch.value?.blockOrthogroups)
      ? clickedPairwiseMatch.value.blockOrthogroups
      : []
  ));
  const selectedPairwiseBlockOrthogroup = computed(() => {
    const selectedId = String(selectedPairwiseBlockOrthogroupId.value || '').trim();
    if (!selectedId) return null;
    return pairwiseBlockOrthogroups.value.find((group) => String(group?.id || '').trim() === selectedId) || null;
  });
  const renderedPairwiseMatchSections = computed(() => {
    const sections = Array.isArray(clickedPairwiseMatch.value?.sections)
      ? clickedPairwiseMatch.value.sections
      : [];
    const selectedGroup = selectedPairwiseBlockOrthogroup.value;
    if (!selectedGroup) return sections;
    const selectedSection = {
      title: 'Selected similarity group',
      rows: Array.isArray(selectedGroup.detailRows) ? selectedGroup.detailRows : [],
      memberRows: Array.isArray(selectedGroup.memberRows) ? selectedGroup.memberRows : [],
      memberCopyText: selectedGroup.memberCopyText || '',
      memberNtFasta: selectedGroup.memberNtFasta || '',
      memberAaFasta: selectedGroup.memberAaFasta || '',
      memberNtFilename: selectedGroup.memberNtFilename || '',
      memberAaFilename: selectedGroup.memberAaFilename || ''
    };
    const output = [];
    sections.forEach((section) => {
      output.push(section);
      if (Array.isArray(section?.blockOrthogroups)) output.push(selectedSection);
    });
    return output;
  });
  const selectPairwiseBlockOrthogroup = (group) => {
    selectedPairwiseBlockOrthogroupId.value = String(group?.id || '').trim();
  };
  const openPairwiseFeatureRow = (row, event) => {
    if (!row?.canOpen || !row?.feature) return;
    openFeatureEditorForFeature(row.feature, event);
  };

  watch(clickedPairwiseMatch, (match) => {
    selectedPairwiseBlockOrthogroupId.value = '';
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return;
    const matchId = String(match?.id || '').trim();
    const orthogroupId = String(
      match?.matchKind === 'orthogroup' ? match?.orthogroupId : ''
    ).trim();
    svg.querySelectorAll(PAIRWISE_MATCH_SELECTOR).forEach((element) => {
      const elementMatchId = String(
        element.getAttribute('data-gbdraw-match-id') ||
        element.getAttribute('data-gbdraw-pairwise-match-id') ||
        ''
      ).trim();
      const elementOrthogroupIds = String(element.getAttribute('data-orthogroup-id') || '')
        .split(';')
        .map((value) => value.trim())
        .filter(Boolean);
      const elementMatchKind = String(element.getAttribute('data-match-kind') || '').trim().toLowerCase();
      const isOrthogroupMatch = elementMatchKind === 'orthogroup' || (
        !elementMatchKind &&
        !element.hasAttribute('data-collinearity-block-id') &&
        elementOrthogroupIds.length > 0
      );
      const selected = orthogroupId
        ? isOrthogroupMatch && elementOrthogroupIds.includes(orthogroupId)
        : Boolean(matchId) && elementMatchId === matchId;
      element.classList.toggle('gbdraw-match-selected', selected);
    });
  });

  const onPairwiseMatchPopupDrag = (event) => {
    if (!pairwiseMatchPopupDrag.active || pairwiseMatchPopupResize.active) return;
    const popup = pairwiseMatchPopupRef.value;
    const width = popup?.offsetWidth || 420;
    const height = popup?.offsetHeight || 360;
    const margin = FEATURE_POPUP_MARGIN;
    const maxX = Math.max(margin, window.innerWidth - width - margin);
    const maxY = Math.max(margin, window.innerHeight - height - margin);
    const nextX = event.clientX - pairwiseMatchPopupDrag.offsetX;
    const nextY = event.clientY - pairwiseMatchPopupDrag.offsetY;
    clickedPairwiseMatchPos.x = Math.min(Math.max(nextX, margin), maxX);
    clickedPairwiseMatchPos.y = Math.min(Math.max(nextY, margin), maxY);
  };

  const endPairwiseMatchPopupDrag = () => {
    if (!pairwiseMatchPopupDrag.active) return;
    pairwiseMatchPopupDrag.active = false;
    document.removeEventListener('mousemove', onPairwiseMatchPopupDrag);
    document.removeEventListener('mouseup', endPairwiseMatchPopupDrag);
  };

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

  const startPairwiseMatchPopupDrag = (event) => {
    if (event.button !== 0) return;
    if (!clickedPairwiseMatch.value) return;
    if (pairwiseMatchPopupResize.active) return;
    if (isInteractiveTarget(event.target)) return;
    const popup = pairwiseMatchPopupRef.value;
    if (!popup) return;
    const rect = popup.getBoundingClientRect();
    pairwiseMatchPopupDrag.active = true;
    pairwiseMatchPopupDrag.offsetX = event.clientX - rect.left;
    pairwiseMatchPopupDrag.offsetY = event.clientY - rect.top;
    document.addEventListener('mousemove', onPairwiseMatchPopupDrag);
    document.addEventListener('mouseup', endPairwiseMatchPopupDrag);
    event.preventDefault();
  };

  const startPairwiseMatchPopupResize = (event) => {
    if (event.button !== 0) return;
    if (!clickedPairwiseMatch.value) return;
    const popup = pairwiseMatchPopupRef.value;
    if (!popup) return;
    const rect = popup.getBoundingClientRect();
    const constraints = getPairwiseMatchPopupConstraints(rect.left, rect.top);

    pairwiseMatchPopupDrag.active = false;
    document.removeEventListener('mousemove', onPairwiseMatchPopupDrag);
    document.removeEventListener('mouseup', endPairwiseMatchPopupDrag);

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

  const onFeaturePopupDrag = (event) => {
    if (!featurePopupDrag.active || featurePopupResize.active) return;
    const popup = featurePopupRef.value;
    const width = popup?.offsetWidth || 360;
    const height = popup?.offsetHeight || 260;
    const margin = 12;
    const maxX = Math.max(margin, window.innerWidth - width - margin);
    const maxY = Math.max(margin, window.innerHeight - height - margin);
    const nextX = event.clientX - featurePopupDrag.offsetX;
    const nextY = event.clientY - featurePopupDrag.offsetY;
    clickedFeaturePos.x = Math.min(Math.max(nextX, margin), maxX);
    clickedFeaturePos.y = Math.min(Math.max(nextY, margin), maxY);
  };

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

  const endFeaturePopupDrag = () => {
    if (!featurePopupDrag.active) return;
    featurePopupDrag.active = false;
    document.removeEventListener('mousemove', onFeaturePopupDrag);
    document.removeEventListener('mouseup', endFeaturePopupDrag);
  };

  const startFeaturePopupDrag = (event) => {
    if (event.button !== 0) return;
    if (!clickedFeature.value) return;
    if (featurePopupResize.active) return;
    if (isInteractiveTarget(event.target)) return;
    const popup = featurePopupRef.value;
    if (!popup) return;
    const rect = popup.getBoundingClientRect();
    featurePopupDrag.active = true;
    featurePopupDrag.offsetX = event.clientX - rect.left;
    featurePopupDrag.offsetY = event.clientY - rect.top;
    document.addEventListener('mousemove', onFeaturePopupDrag);
    document.addEventListener('mouseup', endFeaturePopupDrag);
    event.preventDefault();
  };

  const startFeaturePopupResize = (event) => {
    if (event.button !== 0) return;
    if (!clickedFeature.value) return;
    const popup = featurePopupRef.value;
    if (!popup) return;
    const rect = popup.getBoundingClientRect();
    const constraints = getFeaturePopupConstraints(rect.left, rect.top);

    featurePopupDrag.active = false;
    document.removeEventListener('mousemove', onFeaturePopupDrag);
    document.removeEventListener('mouseup', endFeaturePopupDrag);

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
    return title || 'Untitled session';
  });

  const canUseLinearRulerOnAxis = computed(
    () => {
      const drawing = state.activeDrawing();
      return (
        drawing.form.show_scale !== false &&
        drawing.form.scale_style === 'ruler' &&
        ['above', 'below'].includes(drawing.form.linear_track_layout)
      );
    }
  );

  const canUseCircularScaleStyling = computed(
    () => {
      const drawing = state.activeDrawing();
      return (
        drawing.adv.circular_track_slots_enabled
          ? drawing.adv.circular_track_slots.some((slot) => (
              slot?.renderer === 'ticks' &&
              circularTrackSlotEditor.circularTrackSlotEffectiveEnabled(slot)
            ))
          : drawing.form.show_scale !== false
      );
    }
  );

  // Popup header line "<record ID>: <location> · <gene>", from payload fields
  // only; the gene is shown when it differs from the label above it.
  const clickedFeatureSummary = computed(() => {
    const cf = clickedFeature.value;
    const located = [cf?.recordId, cf?.location || ''].filter(Boolean).join(': ');
    return cf?.gene && cf.gene !== cf.label ? `${located} · ${cf.gene}` : located;
  });

  const downloadText = (filename, text, type = 'text/plain;charset=utf-8') => {
    const value = String(text ?? '');
    if (!value) return;
    downloadTextFile(String(filename || 'gbdraw.txt'), value, type);
  };

  let latestExportOperation = 0;
  const failedInteractiveSvgExport = ref(null);
  const canRetryInteractiveSvgExport = computed(() => Boolean(failedInteractiveSvgExport.value
    && errorLog.value === failedInteractiveSvgExport.value));
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
      const result = await exportMethod(snapshot, {
        loadPdfFont: async (filename) => {
          const result = await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.READ_PDF_FONT, { filename });
          return result.result.base64;
        }
      });
      if (operationId === latestExportOperation && errorLog.value === previousError && previousError?.operation?.startsWith('export-')) errorLog.value = null;
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

  const specificRuleLegendOptions = computed(() => {
    const drawing = state.activeDrawing();
    const byCaption = new Map();
    for (const rule of drawing.manualSpecificRules) {
      if (!rule) continue;
      const caption = String(rule.cap || '').trim();
      if (!caption) continue;
      const key = caption.toLowerCase();
      const isHashRule = String(rule.qual || '').toLowerCase() === 'hash';
      const existing = byCaption.get(key);
      if (!existing || (existing.isHashRule && !isHashRule)) {
        byCaption.set(key, {
          caption,
          color: String(rule.color || '').trim(),
          isHashRule
        });
      }
    }
    return Array.from(byCaption.values())
      .sort((a, b) => a.caption.localeCompare(b.caption))
      .map(({ caption, color }) => ({ caption, color }));
  });

  const featureVisibilityFeatureSuggestions = computed(() => {
    const values = new Set(['*', ...featureKeys]);
    for (const feat of Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []) {
      const type = String(feat?.type || '').trim();
      if (type) values.add(type);
    }
    return Array.from(values);
  });

  const editSessionTitle = () => {
    const busy = sessionOperationAvailability();
    if (busy) return busy;
    const current = normalizeSessionTitle(sessionTitle.value);
    const input = prompt('Session title', current);
    if (input === null) return;
    sessionTitle.value = normalizeSessionTitle(input);
  };

  const saveSessionWithTitle = () => exportSession(null, {
    availability: sessionSaveLoadAvailability,
    recordDisplayRows: recordDisplayControls.allRows,
    // Save writes every Result: the other mode's slot goes beside the shown one (E1).
    readOtherModeArtifact: () => artifactSlots[mode.value === 'linear' ? 'circular' : 'linear'],
    resolveTitle: () => {
      let title = normalizeSessionTitle(sessionTitle.value);
      if (!title) {
        const input = prompt('Session title', '');
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
      recordSessionLifecycleEvent('session-save-catalog-preparation-start');
      /** @type {Awaited<ReturnType<typeof prepareLinearRecordCatalog>>['catalog']} */
      let catalog = null;
      let error = '';
      if (draftRequest) {
        ({ catalog, error } = await prepareLinearRecordCatalog({ privateCandidate: true }));
        await afterPaint();
      }
      recordSessionLifecycleEvent('session-save-catalog-preparation-end', {
        reusedCommittedSession: !draftRequest
      });
      if (error) throw error;
      return { linearRecordCatalog: catalog };
    },
    onError: (error) => { errorLog.value = normalizeUserFacingError(error); }
  });

  const openFeatureEditorFromList = (feat, event) => {
    return openFeatureEditorForFeature(feat, event);
  };

  const getCircularRecordOrderLabel = (selector) => {
    const normalized = String(selector || '').trim();
    if (!normalized) return '';
    const records = Array.isArray(circularRecordList.value) ? circularRecordList.value : [];
    const matched = records.find((entry) => String(entry?.selector || '').trim() === normalized);
    if (!matched) return normalized;
    return `${normalized} (${String(matched.record_id || '').trim() || 'Unknown'})`;
  };

  const circularRecordDiscoveryState = computed(getCircularRecordDiscoveryState);
  const circularRecordInspectionEnabled = computed(() => semanticMutationAvailable.value
    && mode.value === 'circular'
    && circularRecordDiscoveryState.value.hasInput
    && circularRecordDiscoveryState.value.status !== 'loading'
    && !state.semanticFileWatchersSuppressed.value);
  const inspectCircularSourceRecords = () => sessionOperationAvailability() || (circularRecordInspectionEnabled.value
    ? refreshCircularRecordOrder()
    : Promise.resolve({ status: 'unavailable', reason: 'Source inspection is unavailable.' }));
  const circularRecordPresentationEntries = () => buildDisambiguatedRecordEntries(
    circularRecordDiscoveryState.value.records.map(
      (record) => ({
        ...record,
        recordId: record?.record_id ?? record?.recordId,
        recordLength: record?.record_length ?? record?.recordLength
      })
    )
  );

  const circularRecordPresentationOptions = computed(() => {
    const drawing = state.activeDrawing();
    const entries = circularRecordPresentationEntries();
    const current = String(drawing.form.circular_record_selector || '').trim();
    const selection = resolveDisambiguatedRecordSelection(entries, current);
    // Without an inspected catalog the saved selector is unverified, not missing.
    const inspected = circularRecordDiscoveryState.value.status === 'ready';
    const automaticLabel = entries.length > 1 || drawing.adv.circular_grouping_intent === 'batch'
      ? 'All records (separate diagrams)'
      : 'Automatic (only record)';
    return [
      { value: '', label: automaticLabel, synthetic: false },
      ...(current && selection.status !== 'resolved'
        ? [{
            value: current,
            label: `${current} (${!inspected ? 'not inspected' : selection.status === 'ambiguous' ? 'ambiguous' : 'not found'})`,
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
      return `Record selector '${current}' is ambiguous in the current input.`;
    }
    if (selection.status === 'missing') {
      return `Record selector '${current}' was not found in the current input.`;
    }
    return '';
  });

  const circularSingleRecordPresentationEnabled = computed(() => {
    const drawing = state.activeDrawing();
    if (mode.value !== 'circular' || drawing.form.multi_record_canvas
      || drawing.adv.circular_grouping_intent === 'batch'
      || circularRecordDiscoveryState.value.status !== 'ready') return false;
    const entries = circularRecordPresentationEntries();
    const selection = resolveDisambiguatedRecordSelection(entries, drawing.form.circular_record_selector);
    return selection.status === 'resolved'
      || (selection.status === 'unspecified' && entries.length === 1);
  });

  const circularRecordSelectionEnabled = computed(() => semanticMutationAvailable.value
    && mode.value === 'circular'
    && !state.activeDrawing().form.multi_record_canvas && circularRecordDiscoveryState.value.status === 'ready');
  // CI-03: the per-record presentation a source file carries in each mode.
  // Replacing or removing the file resets it in the replacement's History step.
  const RETIRED_LINEAR_RECORD_PRESENTATION = Object.freeze({ region_record_id: '', region_start: null,
    region_end: null, region_reverse: false, definition: '', record_subtitle: '' });
  /** @param {DrawingState} drawing */
  const retireCircularRecordSelector = ({ form, adv }) => {
    form.circular_record_selector = '';
    if (adv.circular_grouping_intent === 'single') adv.circular_grouping_intent = 'auto';
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
    if (!previous) return;
    if (sourceReplacementRetiresPresentation({
      removed: !nextValue, previousRecordCount: circularRecordList.value.length
    })) retireCircularRecordPresentation(drawing);
    else retireCircularRecordSelector(drawing);
  };
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
    for (let section = control.closest('details'); section; section = section.parentElement?.closest('details') ?? null) {
      section.open = true;
    }
    await nextTick();
    control.focus();
    control.scrollIntoView({ block: 'nearest' });
    return { status: 'ok' };
  };
  watch(() => [mode.value, circularSingleRecordPresentationEnabled.value,
    circularRecordDiscoveryState.value.primaryFile, circularRecordDiscoveryState.value.pairedFile,
    state.activeDrawing().form.circular_record_selector], async (current, previous = []) => {
    if (!current[1] || current.every((value, index) => Object.is(value, previous[index]))) return;
    const origin = document.activeElement;
    const pane = document.querySelector('.settings-scroll');
    const scrollOwner = pane && getComputedStyle(pane).overflowY !== 'visible'
      ? pane : document.scrollingElement;
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

  const buildDefaultCircularRecordPositions = () => {
    const selectors = Array.isArray(circularRecordList.value)
      ? circularRecordList.value
          .map((entry) => String(entry?.selector || '').trim())
          .filter(Boolean)
      : [];
    if (selectors.length === 0) return [];
    const cols = Math.ceil(Math.sqrt(selectors.length));
    return selectors.map((selector, index) => ({
      selector,
      row: Math.floor(index / cols) + 1
    }));
  };

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
      .map((entry, index) => ({ ...entry, __index: index }))
      .sort((left, right) => {
        const leftRow = Number(left.row);
        const rightRow = Number(right.row);
        if (leftRow !== rightRow) return leftRow - rightRow;
        return left.__index - right.__index;
      })
      .map(({ __index, ...entry }) => entry);
    drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length, ...sorted);
  };

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

  const canMoveCircularRecordOrderUp = (index) => {
    const drawing = state.drawings.circular;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx <= 0 || idx >= drawing.adv.multi_record_positions.length) return false;
    const currentRow = getCircularRecordRow(drawing.adv.multi_record_positions[idx]);
    const prevRow = getCircularRecordRow(drawing.adv.multi_record_positions[idx - 1]);
    return currentRow === prevRow;
  };

  const canMoveCircularRecordOrderDown = (index) => {
    const drawing = state.drawings.circular;
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.adv.multi_record_positions.length - 1) return false;
    const currentRow = getCircularRecordRow(drawing.adv.multi_record_positions[idx]);
    const nextRow = getCircularRecordRow(drawing.adv.multi_record_positions[idx + 1]);
    return currentRow === nextRow;
  };

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

  /** @param {DrawingState} drawing */
  const applyLinearSeqMutation = (
    drawing, items,
    {
      preserveLosatCacheInfo = false,
      layoutEntries = drawing.linearRecordRows,
      alignmentMutation = 'source set changed.'
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
    pendingLinearRecordExpansions.forEach((uid) => {
      if (!activeUids.has(uid)) pendingLinearRecordExpansions.delete(uid);
    });
    pendingLinearMetadataInference.forEach((uid) => {
      if (!activeUids.has(uid)) pendingLinearMetadataInference.delete(uid);
    });
    const nextRows = reconcileLinearRecordLayout(linearSeqs, layoutEntries);
    drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...nextRows);
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
  const trapDialogFocus = (event) => {
    const dialog = event.currentTarget?.querySelector('[role="dialog"]');
    const controls = Array.from(dialog?.querySelectorAll('button:not(:disabled)') || []);
    if (!controls.length) return;
    const first = controls[0];
    const last = controls.at(-1);
    if (event.shiftKey && (document.activeElement === first || !dialog.contains(document.activeElement))) {
      event.preventDefault();
      last.focus();
    } else if (!event.shiftKey && (document.activeElement === last || !dialog.contains(document.activeElement))) {
      event.preventDefault();
      first.focus();
    }
  };
  const focusLinearSourceAfterRemoval = async (sourceIndex) => {
    await nextTick();
    const cards = document.querySelectorAll('[data-linear-source-card]');
    const targetIndex = Math.max(0, Math.min(Number(sourceIndex) || 0, cards.length - 1));
    const target = /** @type {HTMLElement | null} */ (cards[targetIndex]?.querySelector('.upload-zone[tabindex="0"]')
      || document.querySelector('[data-linear-file-add]'));
    if (typeof target?.focus === 'function') target.focus();
  };
  /** @param {EventTarget | null} [returnFocus] */
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
  const applyLinearSourceRemoval = async (intent) => {
    const drawing = state.drawings.linear;
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
    const next = [...plan.retainedSequences];
    if (intent === 'clear') next.splice(plan.insertionIndex, 0, createLinearSeq());
    const operation = await history.runUndoable(
      intent === 'clear' ? 'Clear Linear File' : 'Delete Linear File',
      () => applyLinearSeqMutation(drawing, next)
    );
    closeLinearSourceRemovalDialog({ restoreFocus: false });
    await focusLinearSourceAfterRemoval(plan.sourceIndex);
    return operation;
  };
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
    const group = linearSourceGroups.value.find((entry) => (
      entry.records.some(({ sequence }) => sequence.uid === seq.uid)
    ));
    const members = new Set(group.records.map(({ sequence }) => sequence.uid));
    const replacement = createLinearSeq({
      ...group.sequence,
      [field]: nextValue,
      ...(field === 'gb' && nextValue ? { file_definition: '', file_subtitle: '' } : {}),
      ...(seq[field] && sourceReplacementRetiresPresentation({
        removed: !nextValue, previousRecordCount: group.records.length
      }) ? RETIRED_LINEAR_RECORD_PRESENTATION : {})
    });
    const keepSource = field === 'gb' ? Boolean(nextValue) : Boolean(replacement.gff || replacement.fasta);
    applyLinearSeqMutation(drawing, linearSeqs.flatMap((entry) => (
      entry.uid === group.uid ? (keepSource ? [replacement] : [])
        : members.has(entry.uid) ? [] : [entry]
    )), { alignmentMutation: 'source replaced.' });
    if (keepSource) pendingLinearRecordExpansions.add(replacement.uid);
    if (keepSource && field === 'gb') pendingLinearMetadataInference.add(replacement.uid);
  };

  /** @param {DrawingState} drawing */
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
        ? 'File order is unavailable because Record Layout is custom. Use Advanced comparison and layout → Record Layout to give each File its own consecutive rows.'
        : ''
    );
  });
  const canMoveLinearSource = (sourceIndex, direction) => {
    const drawing = state.drawings.linear;
    return (
      linearSourceMovePlan(drawing, sourceIndex, direction).allowed
    );
  };

  const moveLinearSource = (sourceIndex, direction) => {
    const drawing = state.drawings.linear;
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const plan = linearSourceMovePlan(drawing, sourceIndex, direction);
    if (!plan.allowed) return;
    const next = moveLinearSourceGroup(linearSeqs, sourceIndex, direction);
    return history.runUndoable('Move File', () => {
      applyLinearSeqMutation(drawing, next, {
        preserveLosatCacheInfo: true,
        layoutEntries: plan.rows,
        alignmentMutation: 'stable-reorder'
      });
    });
  };

  const setLinearRecordSelector = (sequence, value) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!sequence) return false;
    const next = String(value || '');
    if (String(sequence.region_record_id || '') === next) return false;
    similarityAlignmentActions?.clearForMutation?.('record selector changed.');
    sequence.region_record_id = next;
    // The inferred definition follows the record the row now selects (D-12).
    if (lInputType.value === 'gb' && linearRecordSelector.statusFor(sequence) === 'ready') {
      sequence.inferred_definition = inferredDefinitionForRecord(linearRecordSelector.recordsFor(sequence), next);
    }
    return true;
  };

  const setCircularInputType = (value) => {
    const busy = sessionOperationAvailability();
    if (busy) return busy;
    cInputType.value = value;
    return { status: 'ok' };
  };

  const setLinearInputType = (value) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const next = value === 'gff' ? 'gff' : 'gb';
    if (lInputType.value === next) return false;
    similarityAlignmentActions?.clearForMutation?.('source type changed.');
    lInputType.value = next;
    return true;
  };

  const setLinearRecordCrop = (sequence, field, value) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!sequence || !['region_start', 'region_end'].includes(field)) return false;
    const next = value === '' || value === null || value === undefined
      ? null
      : Number(value);
    if (Object.is(sequence[field], next)) return false;
    similarityAlignmentActions?.clearForMutation?.('record crop changed.');
    sequence[field] = next;
    return true;
  };

  const resetLinearRecordDefinition = (seq) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!seq) return;
    history.runUndoable('Reset record definition', () => {
      seq.definition = '';
    });
  };

  const resetLinearRecordSubtitle = (seq) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!seq) return;
    history.runUndoable('Reset record subtitle', () => {
      seq.record_subtitle = '';
    });
  };

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

  return {
    recordDisplayControls,
    featureRecordRotationDraft: featureRecordRotation.draft,
    recordActionsExpanded,
    toggleRecordActions,
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
    retryLabelImportFailure: undoableAction('Load label edits', featureActions.retryLabelImportFailure),
    editLabelImportFailure: featureActions.editLabelImportFailure,
    canRetryFeatureEditTableImport: featureActions.canRetryFeatureEditTableImport,
    retryFeatureEditTableImport: undoableAction('Load feature edits', featureActions.retryFeatureEditTableImport),
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
    specificRuleNotice,
    importAnnotationTableFile: undoableAction('Import annotations', annotationEditor.importAnnotationTableFile),
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
    annotationRecordDisabled: annotationEditor.recordIsDisabled,
    annotationRecordMissingMessage: annotationEditor.recordMissingMessage,
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
    updateCircularTrackSlotMeasure: circularTrackSlotEditor.updateCircularTrackSlotMeasure,
    circularTrackRendererLabel: circularTrackSlotEditor.circularTrackRendererLabel,
    resetCircularTrackSlotsFromSimpleControls: circularTrackSlotEditor.resetCircularTrackSlotsFromSimpleControls,
    resetCircularTrackSlotsToPreset: circularTrackSlotEditor.resetCircularTrackSlotsToPreset,
    setCircularTrackSlotsEnabled: circularTrackSlotEditor.setCircularTrackSlotsEnabled,
    setCircularGcSuppressed: circularTrackSlotEditor.setCircularGcSuppressed,
    setCircularSkewSuppressed: circularTrackSlotEditor.setCircularSkewSuppressed,
    addCircularTrackSlot: circularTrackSlotEditor.addCircularTrackSlot,
    canAddCircularTrackRenderer: circularTrackSlotEditor.canAddCircularTrackRenderer,
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
    canAddLinearTrackRenderer: linearTrackSlotEditor.canAddLinearTrackRenderer,
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
    linearTrackSlotPlacementLabel: linearTrackSlotEditor.linearTrackSlotPlacementLabel,
    linearTrackSlotUsesPresetGeometry: linearTrackSlotEditor.linearTrackSlotUsesPresetGeometry,
    losat: drawingMember('losat'),
    losatExecution,
    richFeaturePopup,
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
    selectedOrthogroupMembersByRecord: orthogroupActions.selectedOrthogroupMembersByRecord,
    resolveOrthogroupName: orthogroupActions.resolveOrthogroupName,
    resolveOrthogroupDescription: orthogroupActions.resolveOrthogroupDescription,
    orthogroupScope: orthogroupActions.orthogroupScope,
    orthogroupScopeLabel: orthogroupActions.orthogroupScopeLabel,
    orthogroupRows: orthogroupActions.orthogroupRows,
    isOrthogroupRenamed: orthogroupActions.isOrthogroupRenamed,
    getOrthogroupSequenceCount: orthogroupActions.getOrthogroupSequenceCount,
    hasOrthogroupSequence: orthogroupActions.hasOrthogroupSequence,
    hasOrthogroupMemberSequence: orthogroupActions.hasOrthogroupMemberSequence,
    orthogroupCopyFeedbackLabel: orthogroupActions.orthogroupCopyFeedbackLabel,
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
    similarityAlignmentUnresolvedCount: similarityAlignmentActions.unresolvedCount,
    similarityAlignmentApplyDisabledReason: similarityAlignmentActions.applyDisabledReason,
    similarityAlignmentDrawerReferenceKey: similarityAlignmentActions.drawerReferenceKey,
    similarityAlignmentPlanInspector: similarityAlignmentActions.activePlanInspector,
    canApplySimilarityAlignment: similarityAlignmentActions.canApply,
    resetSimilarityAlignment: similarityAlignmentActions.resetAlignment,
    reselectSimilarityAlignmentReference,
    startSimilarityAlignmentFromPopup: similarityAlignmentActions.startFromPopup,
    startSimilarityAlignmentFromDrawer,
    similarityAlignmentDrawerReferenceOptions: similarityAlignmentActions.drawerReferenceOptions,
    setSimilarityAlignmentDrawerReference: similarityAlignmentActions.setDrawerReference,
    similarityAlignmentDrawerDisabledReason: similarityAlignmentActions.drawerDisabledReason,
    selectSimilarityAlignmentCandidate: similarityAlignmentActions.selectCandidate,
    skipSimilarityAlignmentRecord: similarityAlignmentActions.skipRecord,
    similarityAlignmentDirectionPreview: similarityAlignmentActions.directionPreview,
    setSimilarityAlignmentDirectionMode: similarityAlignmentActions.setDirectionMode,
    setSimilarityAlignmentCustomDirection: similarityAlignmentActions.setCustomDirection,
    similarityAlignmentResetPreview: similarityAlignmentActions.resetPreview,
    similarityAlignmentResetDialogOpen: similarityAlignmentActions.resetDialogOpen,
    similarityAlignmentResetScope: similarityAlignmentActions.resetScope,
    openSimilarityAlignmentReset,
    cancelSimilarityAlignmentReset: similarityAlignmentActions.cancelReset,
    applySimilarityAlignmentReset: () => applySimilarityAlignmentDialog(true),
    applySimilarityAlignmentDraft: () => applySimilarityAlignmentDialog(),
    cancelSimilarityAlignmentDraft: cancelSimilarityAlignmentDialog,
    cancelSimilarityAlignmentDialog,
    similarityAlignmentPaletteRef,
    similarityAlignmentPaletteStyle,
    similarityAlignmentCompact,
    similarityAlignmentEditorDisabledReason,
    similarityAlignmentCanvasHover,
    startSimilarityAlignmentPaletteDrag,
    previewSimilarityAlignmentCandidate: similarityAlignmentActions.previewCandidate,
    clearSimilarityAlignmentCandidatePreview: similarityAlignmentActions.clearCandidatePreview,
    isRightDrawerTabAvailable: rightDrawerActions.isRightDrawerTabAvailable,
    openRightDrawerTab: rightDrawerActions.openRightDrawerTab,
    toggleRightDrawer: rightDrawerActions.toggleRightDrawer,
    closeRightDrawer: rightDrawerActions.closeRightDrawer,
    openOrthogroupInDrawer,
    circularRecordList,
    refreshCircularRecordOrder,
    waitForAuxiliaryFileImport, canRetryAuxiliaryImportFailure,
    retryAuxiliaryImportFailure: () => history.runUndoableCheckpoint('Change uploaded file', retryAuxiliaryImportFailure, { shouldCommit: result => result !== false }),
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
    updatePalette,
    resetColors,
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
    startFeaturePopupDrag,
    startFeaturePopupResize,
    pairwiseMatchPopupRef,
    pairwiseMatchPopupStyle,
    startPairwiseMatchPopupDrag,
    startPairwiseMatchPopupResize,
    selectedPairwiseBlockOrthogroupId,
    renderedPairwiseMatchSections,
    selectPairwiseBlockOrthogroup,
    openPairwiseFeatureRow,
    clickedFeatureSummary,
      copyText: copyTextToClipboard,
    downloadText,
    canUseClickedOrthogroupActions,
    clickedOrthogroupDetail,
    alignByClickedOrthogroup,
    highlightClickedOrthogroup,
    clearOrthogroupHighlight,
    openClickedOrthogroupInEditor,
    specificRuleLegendOptions,
    updateClickedFeatureColor: updateClickedFeatureColorWithHistory,
    updateClickedFeatureVisibility: updateClickedFeatureVisibilityWithHistory,
    handleFeatureVisibilityScopeChoice: handleFeatureVisibilityScopeChoiceWithHistory,
    handleLegendNameCommit: handleLegendNameCommitWithHistory,
    selectLegendNameOption,
    resetClickedFeatureFillColor: resetClickedFeatureFillColorWithHistory,
    updateClickedFeatureStroke: updateClickedFeatureStrokeWithHistory,
    getFeatureStrokeColorValue,
    setClickedFeatureStrokeColorValue: setClickedFeatureStrokeColorValueWithHistory,
    setClickedFeatureStrokeWidthValue: setClickedFeatureStrokeWidthValueWithHistory,
    resetClickedFeatureStroke: resetClickedFeatureStrokeWithHistory,
    featureStyleScopeDialog,
    featureVisibilityScopeDialog,
    handleColorScopeChoice: handleColorScopeChoiceWithHistory,
    handleFeatureStyleScopeChoice: handleFeatureStyleScopeChoiceWithHistory,
    legendRenameDialog,
    handleLegendRenameChoice: handleLegendRenameChoiceWithHistory,
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
    loadFeatureEditTable: undoableAction('Load feature edits', featureActions.loadFeatureEditTable),
    syncLabelEditor,
    openFeatureEditorFromList,
    legendEntries: drawingMember('legendEntries'),
    newLegendCaption,
    newLegendColor,
    updateLegendEntryColor,
    renameLegendEntry,
    // A Legend that gains or loses a row is one checkpoint step, so Undo and
    // Redo return the Legend as it was laid out, canvas included (OV-125).
    deleteLegendEntry: /** @param {number} index */ (index) => history.runUndoableCheckpoint(
      'Delete legend item',
      () => deleteLegendEntry(index)
    ),
    addNewLegendEntry: () => history.runUndoableCheckpoint('Add legend item', addNewLegendEntry),
    // OV-154: a Restore returns deleted rows in one checkpoint step, as a delete
    // removes them; the palette then reaches the returned rows.
    deletedLegendEntries: drawingMember('deletedLegendEntries'),
    restoreDeletedLegendEntry: /** @param {number} index */ (index) => restoreLegendItems('Restore legend item', [index]),
    restoreAllDeletedLegendEntries: () => restoreLegendItems('Restore legend items', null),
    moveLegendEntryUp,
    moveLegendEntryDown,
    sortLegendEntries,
    sortLegendEntriesByDefault,
    resetLegendPosition: undoableAction('Reset legend position', resetLegendPosition),
    getLegendEntryStrokeColor,
    getLegendEntryStrokeWidth,
    isLegendStrokeOptionsOpen,
    toggleLegendStrokeOptions,
    setLegendEntryStrokeColorValue: setLegendEntryStrokeColorValueWithHistory,
    updateLegendEntryStrokeColor,
    updateLegendEntryStrokeWidth: undoableAction('Change legend stroke width', updateLegendEntryStrokeWidth),
    resetLegendEntryStroke: undoableAction('Reset legend stroke', resetLegendEntryStroke),
    resetAllStrokes,
    resetAllPositions: undoableAction('Reset positions', resetAllPositions),
    resetLayout: undoableAction('Reset layout', resetLayout),
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
    saveSessionWithTitle,
    editSessionTitle,
    importSession,
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
