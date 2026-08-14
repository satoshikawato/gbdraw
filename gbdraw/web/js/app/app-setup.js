import { state, createLinearSeq, normalizeLinearSeqList } from '../state.js';
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
  exportSession,
  importSession as importSessionFromFile,
  SESSION_VERSION,
  serializeActiveRenderFiles,
  serializeResults,
  setPreviewRuntime
} from '../services/config.js';
import { createHistoryManager } from '../services/history.js';
import { createHistoryFileStore } from '../services/history-files.js';
import { createHistorySnapshotService } from '../services/history-snapshot.js';
import { cloneJsonData } from '../services/json-clone.js';
import { readFileText } from '../services/file-content-cache.js';
import { serializeCleanSvg } from '../services/svg-serialization.js';
import { copyTextToClipboard } from '../utils/clipboard.js';
import { downloadTextFile } from '../services/text-download.js';
import { resetLayoutState, resetSettings as resetSettingsState } from '../services/reset.js';
import { disposeDiagramGenerationWorker } from '../services/diagram-generation.js';
import { recordSessionLifecycleEvent } from '../services/runtime-test-hooks.js';
import { createPanZoom, createSidebarResize, setupGlobalUiEvents } from './ui.js';
import { colorValueMode, toNativeColorInputValue } from './color-utils.js';
import { createFeatureEditor } from './feature-editor.js';
import { PAIRWISE_MATCH_SELECTOR } from './pairwise-match-popup.js';
import { createFeatureSelection } from './feature-selection.js';
import { createPreviewFeatureSearch } from './feature-search/preview-actions.js';
import { createSvgStyles } from './svg-styles.js';
import { createLegendManager } from './legend.js';
import { createPaletteLoader } from './palettes.js';
import { createRunAnalysis } from './run-analysis.js';
import { normalizeUserFacingError } from '../services/error-normalization.js';
import { formatElapsedMs, reproducibilityLabel } from './run-info.js';
import { createLegendLayout } from './legend-layout.js';
import { compositionUserDeltas } from './legend-layout/composition-actions.js';
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
import {
  linearRecordPositionTokens,
  moveLinearRecordInRow,
  reconcileLinearRecordLayout,
  setLinearRecordRow as updateLinearRecordRow
} from './linear-record-layout.js';
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
} from './linear-comparisons.js';
import {
  projectLinearComparisonLosatModeSelection,
  projectLinearComparisonLosatpModeSelection,
  projectLinearComparisonUi
} from './comparison-ui.js';
import { discoverGffFastaRecords, discoverSequenceRecords } from './record-discovery.js';
import {
  conservationSourceDescriptors,
  defaultConservationSeriesLabel,
  moveConservationSeriesEntry,
  normalizeFileList,
  orderedConservationSources,
  orderedOptionalConservationFiles,
  parseConservationLabelText,
  reconcileConservationSeries
} from './conservation-series.js';
import {
  getDepthTrackFallbackLabel,
  getDepthTrackLabelFromFile,
  isDepthTrackAutoLabel
} from './depth-tracks.js';
import { comparisonProfileDefault } from '../mode-profiles.js';
import {
  activeDepthTrackIndices,
  clearDepthTrackSourceAt,
  compactDepthFileSlots,
  depthFileSlotsFromValue,
  depthSlotTrackIndex,
  depthTrackCoverageCount,
  depthTrackMatrixWidth,
  ensureDepthTrackConfigCount as ensureDepthTrackConfigCountEntries,
  ensureDepthTrackConfigShape,
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
} from './depth-track-state.js';

const { onMounted, onUnmounted, watch, nextTick, computed, ref, reactive } = window.Vue;

let exportServicePromise = null;

const loadExportService = () => {
  exportServicePromise ??= import('../services/export.js');
  return exportServicePromise;
};

export const createSessionImportRollbackState = ({
  depthTrackUiCounts,
  depthTracks,
  featureListScrollTop,
  featureListScrollRef,
  selectedPairwiseBlockOrthogroupId
}) => ({
  capture: () => ({
    circularDepthTrackUiCount: depthTrackUiCounts.circular,
    depthTracks: cloneJsonData(depthTracks),
    featureListScrollTop: featureListScrollTop.value,
    selectedPairwiseBlockOrthogroupId: selectedPairwiseBlockOrthogroupId.value
  }),
  restore: async (snapshot) => {
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
  }
});

export const createAppSetup = () => {
  const {
    processing,
    processingStatus,
    generationCancelRequested,
    errorLog,
    sessionTitle,
    results,
    selectedResultIndex,
    resultPanelTab,
    lastRunInfo,
    pairwiseMatchFactors,
    matchSequenceRegistry,
    svgContent,
    zoom,
    layoutRepositionMode,
    isPanning,
    panStart,
    canvasPan,
    canvasContainerRef,
    mode,
    cInputType,
    lInputType,
    losatProgram,
    files,
    circularConservation,
    annotationSets,
    selectedAnnotation,
    linearSeqs,
    linearRecordLayoutEnabled,
    linearRecordGap,
    linearRecordRows,
    linearComparisonPlan,
    linearComparisonResolution,
    hasLinearComparisonIntent,
    hasActiveLinearLosatIntent,
    hasActiveLinearUploadIntent,
    form,
    adv,
    losat,
    losatCacheInfo,
    losatThreadingStatus,
    orthogroups,
    featureOrthogroupIndex,
    selectedOrthogroupAlignmentFeature,
    orthogroupNameOverrides,
    orthogroupDescriptionOverrides,
    selectedOrthogroupId,
    orthogroupSearch,
    orthogroupSortMode,
    showRightDrawer,
    rightDrawerTab,
    linearReorderNotice,
    circularRecordList,
    circularRecordDiscovery,
    paletteDefinitions,
    paletteNames,
    selectedPalette,
    currentColors,
    paletteInstantPreviewEnabled,
    appliedPaletteName,
    appliedPaletteColors,
    pendingPaletteName,
    pendingPaletteColors,
    hasPendingPaletteDraft,
    filterMode,
    manualBlacklist,
    manualWhitelist,
    manualSpecificRules,
    newSpecRule,
    specificRulePresets,
    specificRuleQualifierSuggestions,
    selectedSpecificPreset,
    specificRulePresetLoading,
    downloadDpi,
    extractedFeatures,
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
    featureListTopSpacerPx,
    featureListBottomSpacerPx,
    labelSearch,
    editableLabels,
    filteredEditableLabels,
    labelTextFeatureOverrides,
    labelTextBulkOverrides,
    labelTextFeatureOverrideSources,
    labelVisibilityOverrides,
    labelOverrideContextKey,
    labelOverrideBuildWarning,
    autoLabelReflowEnabled,
    labelReflowProcessing,
    labelReflowLastError,
    featureColorOverrides,
    featureVisibilityManualRules,
    featureVisibilityRules,
    featureVisibilityOverrides,
    featureVisibilitySelectorCache,
    featureStrokeOverrides,
    labelLayoutDirtyReason,
    resultGenerationKey,
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
    globalLabelModeDialog,
    sidebarWidth,
    isResizing,
    legendEntries,
    deletedLegendEntries,
    originalLegendOrder,
    originalLegendColors,
    newLegendCaption,
    newLegendColor,
    legendStrokeOverrides,
    legendColorOverrides,
    originalSvgStroke,
    legendDragging,
    legendDragStart,
    legendOriginalTransform,
    legendInitialTransform,
    legendCurrentOffset,
    diagramDragging,
    diagramDragStart,
    diagramOffset,
    diagramElementIds,
    diagramElementOriginalTransforms,
    diagramElements,
    canvasPadding,
    showCanvasControls,
    generatedLegendPosition,
    skipCaptureBaseConfig,
    skipPositionReapply,
    skipExtractOnSvgChange,
    featureKeys,
    defaultColorKeys,
    newColorFeat,
    newColorVal,
    manualPriorityRules,
    newPriorityRule,
    newFeatureToAdd,
    addedLegendCaptions,
    fileLegendCaptions,
    filteredFeatures
  } = state;
  const rightDrawerActions = createRightDrawerController({ state, watch });

  const comparisonHeightValidationError = computed(() => {
    if (
      mode.value !== 'linear' ||
      linearComparisonResolution.value?.hasComparisonIntent !== true
    ) return '';
    return classifyOptionalPositiveNumber(adv.comparison_height).status === 'invalid'
      ? 'Pairwise Match Height must be Auto or a positive finite number.'
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

  const reindexLinearLosatCacheInfo = () => {
    if (!Array.isArray(losatCacheInfo.value)) return;
    const indexByUid = new Map(
      linearSeqs.map((sequence, index) => [String(sequence?.uid || ''), index])
    );
    const resolvedByEdgeKey = new Map(
      linearComparisonResolution.value.edges.map((edge) => [edge.edgeKey, edge])
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

  const invalidateLinearComparisonArtifacts = ({ preserveLosatCacheInfo = false } = {}) => {
    files.linearCanonicalComparisons = [];
    if (Array.isArray(losatCacheInfo.value)) {
      if (preserveLosatCacheInfo) reindexLinearLosatCacheInfo();
      else losatCacheInfo.value = losatCacheInfo.value.filter((entry) => !entry?.edgeKey);
    }
  };

  const replaceLinearComparisonPlan = (nextPlan, { invalidate = true } = {}) => {
    const normalized = normalizeLinearComparisonPlan(nextPlan);
    linearComparisonPlan.mode = normalized.mode;
    linearComparisonPlan.defaultSource = normalized.defaultSource;
    linearComparisonPlan.edges.splice(
      0,
      linearComparisonPlan.edges.length,
      ...normalized.edges
    );
    if (invalidate) invalidateLinearComparisonArtifacts();
  };

  const mutateLinearComparisonPlan = (mutator) => {
    const next = normalizeLinearComparisonPlan(linearComparisonPlan);
    mutator(next);
    replaceLinearComparisonPlan(next);
  };

  const effectiveLinearComparisonLayout = () => (
    linearRecordLayoutEnabled.value ? linearRecordRows : []
  );

  const syncLinearComparisonRecords = ({ invalidate = true } = {}) => {
    const next = reconcileLinearComparisonPlan(linearComparisonPlan, linearSeqs);
    const currentEdges = linearComparisonPlan.edges;
    const unchanged = (
      next.mode === linearComparisonPlan.mode &&
      next.defaultSource === linearComparisonPlan.defaultSource &&
      next.edges.length === currentEdges.length &&
      next.edges.every((edge, index) => sameLinearComparisonEdge(edge, currentEdges[index]))
    );
    if (!unchanged) replaceLinearComparisonPlan(next, { invalidate });
    return !unchanged;
  };

  const syncLinearRecordLayout = ({ preserveLosatCacheInfo = false } = {}) => {
    const next = reconcileLinearRecordLayout(linearSeqs, linearRecordRows);
    const rowsUnchanged = next.length === linearRecordRows.length && next.every((entry, index) => (
      entry.uid === linearRecordRows[index]?.uid && entry.row === linearRecordRows[index]?.row
    ));
    if (!rowsUnchanged) linearRecordRows.splice(0, linearRecordRows.length, ...next);
    const comparisonsChanged = syncLinearComparisonRecords({ invalidate: false });
    if (!rowsUnchanged || comparisonsChanged) {
      invalidateLinearComparisonArtifacts({ preserveLosatCacheInfo });
    }
    return !rowsUnchanged || comparisonsChanged;
  };
  const setLinearRecordRow = (uid, row) => {
    syncLinearRecordLayout({ preserveLosatCacheInfo: true });
    const previous = linearRecordRows.find((entry) => entry.uid === uid)?.row;
    updateLinearRecordRow(linearRecordRows, uid, row);
    if (linearRecordRows.find((entry) => entry.uid === uid)?.row !== previous) {
      invalidateLinearComparisonArtifacts({ preserveLosatCacheInfo: true });
    }
  };
  const setLinearRecordLayoutEnabled = (enabled) => {
    const nextEnabled = Boolean(enabled);
    if (linearRecordLayoutEnabled.value === nextEnabled) return;
    linearRecordLayoutEnabled.value = nextEnabled;
    syncLinearRecordLayout({ preserveLosatCacheInfo: true });
    invalidateLinearComparisonArtifacts({ preserveLosatCacheInfo: true });
  };
  const moveLinearRecordWithinRow = (uid, direction) => {
    const next = moveLinearRecordInRow(linearSeqs, linearRecordRows, uid, direction);
    linearRecordRows.splice(0, linearRecordRows.length, ...next);
    syncLinearComparisonRecords({ invalidate: false });
    invalidateLinearComparisonArtifacts({ preserveLosatCacheInfo: true });
  };

  const linearComparisonGlobalAction = computed(() => {
    if (linearComparisonPlan.mode === LINEAR_COMPARISON_MODES.NONE) return 'none';
    if (linearComparisonPlan.mode !== LINEAR_COMPARISON_MODES.ADJACENT) return 'selected';
    return linearComparisonPlan.defaultSource;
  });

  const linearComparisonUi = computed(() => projectLinearComparisonUi({
    plan: linearComparisonPlan,
    resolution: linearComparisonResolution.value,
    losatProgram: losatProgram.value,
    blastpMode: losat.blastp?.mode,
    filters: {
      min_bitscore: adv.min_bitscore,
      evalue: adv.evalue,
      identity: adv.identity,
      alignment_length: adv.alignment_length
    }
  }));

  const linearComparisonTimeline = computed(() => buildLinearComparisonTimeline({
    sequences: linearSeqs,
    layout: effectiveLinearComparisonLayout(),
    plan: linearComparisonPlan,
    resolution: linearComparisonResolution.value
  }));
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
      sequence?.gb?.name ||
      sequence?.gff?.name ||
      sequence?.fasta?.name ||
      'Record'
    );
  };
  const openLinearComparisonDisclosure = async (disclosureKey) => {
    if (!disclosureKey) return null;
    const details = [...document.querySelectorAll('[data-linear-comparison-disclosure]')]
      .find((element) => element.dataset.linearComparisonDisclosure === disclosureKey);
    if (!details) return null;
    details.open = true;
    await nextTick();
    return details;
  };
  const focusLinearComparisonPair = async (edgeKey) => {
    await openLinearComparisonDisclosure('selected-pairs');
    const container = [...document.querySelectorAll('[data-edge-key]')]
      .find((element) => element.dataset.edgeKey === edgeKey);
    container?.querySelector('input, select, button')?.focus();
  };

  const focusLinearComparisonIssue = async () => {
    const target = linearComparisonUi.value.errorTargets[0];
    if (!target) return false;
    const details = await openLinearComparisonDisclosure(target.disclosureKey);
    let container = details;
    if (target.edgeKey) {
      container = [...document.querySelectorAll('[data-edge-key]')]
        .find((element) => element.dataset.edgeKey === target.edgeKey) || container;
    }
    if (target.edgeId && target.focusTargetKey === 'pair-row') {
      container = [...document.querySelectorAll('[data-linear-unplaced-draft]')]
        .find((element) => element.dataset.linearUnplacedDraft === target.edgeId) || container;
    }
    const marked = container?.querySelector(
      `[data-comparison-focus="${target.focusTargetKey}"]`
    );
    const focusable = marked?.matches?.('input, select, button, [tabindex]')
      ? marked
      : marked?.querySelector?.('input, select, button, [tabindex]');
    (focusable || container?.querySelector?.('input, select, button, [tabindex]'))?.focus();
    return true;
  };

  const setLinearComparisonGlobalAction = (action) => {
    const normalized = String(action || '').trim().toLowerCase();
    mutateLinearComparisonPlan((next) => {
      if (normalized === 'none') {
        next.mode = LINEAR_COMPARISON_MODES.NONE;
        return;
      }
      next.mode = LINEAR_COMPARISON_MODES.ADJACENT;
      next.defaultSource = normalized === LINEAR_COMPARISON_SOURCES.UPLOAD
        ? LINEAR_COMPARISON_SOURCES.UPLOAD
        : LINEAR_COMPARISON_SOURCES.LOSAT;
    });
  };

  const setLinearComparisonLosatMode = (modeKey) => {
    const selection = projectLinearComparisonLosatModeSelection({ modeKey });
    if (!selection.selectable || !selection.patch) return false;
    const nextProgram = selection.patch.losatProgram;
    if (losatProgram.value === nextProgram) return true;
    losatProgram.value = nextProgram;
    invalidateLinearComparisonArtifacts();
    return true;
  };

  const setLinearComparisonLosatpMode = (modeKey) => {
    const selection = projectLinearComparisonLosatpModeSelection({
      plan: linearComparisonPlan,
      modeKey
    });
    if (!selection.selectable || !selection.patch) return false;
    const nextBlastpMode = selection.patch.blastpMode;
    if (losat.blastp?.mode === nextBlastpMode) return true;
    losat.blastp.mode = nextBlastpMode;
    invalidateLinearComparisonArtifacts();
    return true;
  };

  const selectedPlanForEdit = () => {
    if (linearComparisonPlan.mode === LINEAR_COMPARISON_MODES.SELECTED) {
      return normalizeLinearComparisonPlan(linearComparisonPlan);
    }
    return materializeResolvedEdgesAsSelectedPlan(
      linearComparisonPlan,
      linearComparisonResolution.value
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
    if (linearSeqs.length < 2) return;
    syncLinearRecordLayout();
    const [firstPair] = adjacentRowPairs(
      linearSeqs,
      effectiveLinearComparisonLayout()
    );
    const [queryUid, subjectUid] = firstPair || [linearSeqs[0].uid, linearSeqs[1].uid];
    const next = selectedPlanForEdit();
    upsertSelectedComparison(next, { queryUid, subjectUid });
    replaceLinearComparisonPlan(next);
    await focusLinearComparisonPair(linearComparisonEdgeKey(queryUid, subjectUid));
  };
  const omitLinearComparison = (id) => {
    const next = selectedPlanForEdit();
    const index = findEdgeIndex(next.edges, id);
    if (index < 0) return;
    const edge = next.edges[index];
    edge.included = false;
    if (!edge.file && !String(edge.losatFilename || '').trim()) next.edges.splice(index, 1);
    replaceLinearComparisonPlan(next);
  };
  const clearSelectedLinearComparisons = () => {
    const next = selectedPlanForEdit();
    next.edges = next.edges
      .filter((edge) => edge.file || String(edge.losatFilename || '').trim())
      .map((edge) => ({ ...edge, included: false }));
    replaceLinearComparisonPlan(next);
  };
  const setLinearComparisonEndpoint = (id, endpoint, uid) => {
    if (!['queryUid', 'subjectUid'].includes(endpoint)) return;
    const next = selectedPlanForEdit();
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge[endpoint] = String(uid || '');
    edge.included = true;
    replaceLinearComparisonPlan(next);
  };
  const setLinearComparisonSource = (id, source) => {
    const normalized = source === LINEAR_COMPARISON_SOURCES.LOSAT
      ? LINEAR_COMPARISON_SOURCES.LOSAT
      : LINEAR_COMPARISON_SOURCES.UPLOAD;
    const next = selectedPlanForEdit();
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge.source = normalized;
    edge.included = true;
    replaceLinearComparisonPlan(next);
  };
  const setLinearComparisonFile = (id, file) => {
    const next = selectedPlanForEdit();
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge.file = file || null;
    edge.fileActive = Boolean(file);
    edge.source = LINEAR_COMPARISON_SOURCES.UPLOAD;
    edge.included = Boolean(file) || edge.included;
    replaceLinearComparisonPlan(next);
  };
  const reuseLinearComparisonFile = (id) => {
    const next = selectedPlanForEdit();
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge?.file) return;
    edge.fileActive = true;
    edge.source = LINEAR_COMPARISON_SOURCES.UPLOAD;
    edge.included = true;
    replaceLinearComparisonPlan(next);
  };
  const deactivateLinearComparisonFile = (id) => {
    const next = selectedPlanForEdit();
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge?.file) return;
    edge.fileActive = false;
    if (edge.source === LINEAR_COMPARISON_SOURCES.UPLOAD) edge.included = false;
    replaceLinearComparisonPlan(next);
  };
  const setLinearComparisonLosatFilename = (id, value) => {
    const next = selectedPlanForEdit();
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge.losatFilename = String(value || '');
    edge.losatFilenameActive = Boolean(edge.losatFilename.trim());
    replaceLinearComparisonPlan(next, { invalidate: false });
  };
  const reuseLinearComparisonLosatFilename = (id) => {
    const next = selectedPlanForEdit();
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge || !String(edge.losatFilename || '').trim()) return;
    edge.losatFilenameActive = true;
    edge.source = LINEAR_COMPARISON_SOURCES.LOSAT;
    edge.included = true;
    replaceLinearComparisonPlan(next);
  };
  const deactivateLinearComparisonLosatFilename = (id) => {
    const next = selectedPlanForEdit();
    const edge = next.edges.find((entry) => entry.id === id);
    if (!edge) return;
    edge.losatFilenameActive = false;
    replaceLinearComparisonPlan(next, { invalidate: false });
  };
  const updateResolvedLosatFilenameDraft = (edgeKey, updater) => {
    const resolved = linearComparisonResolution.value.edges.find((edge) => edge.edgeKey === edgeKey);
    if (!resolved) return;
    const next = normalizeLinearComparisonPlan(linearComparisonPlan);
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
        source: LINEAR_COMPARISON_SOURCES.LOSAT
      }));
      index = next.edges.length - 1;
    }
    updater(next.edges[index]);
    replaceLinearComparisonPlan(next, { invalidate: false });
  };
  const setResolvedLinearComparisonLosatFilename = (edgeKey, value) => {
    updateResolvedLosatFilenameDraft(edgeKey, (edge) => {
      edge.losatFilename = String(value || '');
      edge.losatFilenameActive = Boolean(edge.losatFilename.trim());
    });
  };
  const reuseResolvedLinearComparisonLosatFilename = (edgeKey) => {
    updateResolvedLosatFilenameDraft(edgeKey, (edge) => {
      if (String(edge.losatFilename || '').trim()) edge.losatFilenameActive = true;
    });
  };
  const deactivateResolvedLinearComparisonLosatFilename = (edgeKey) => {
    updateResolvedLosatFilenameDraft(edgeKey, (edge) => {
      edge.losatFilenameActive = false;
    });
  };
  const addLinearComparisonBatch = (allPairs = false) => {
    syncLinearRecordLayout();
    const next = selectedPlanForEdit();
    adjacentRowPairs(linearSeqs, effectiveLinearComparisonLayout(), allPairs).forEach(([queryUid, subjectUid]) => {
      upsertSelectedComparison(next, { queryUid, subjectUid });
    });
    replaceLinearComparisonPlan(next);
  };
  const setLinearComparisonGapAction = (edgeKey, action) => {
    const pair = linearComparisonPairForEdgeKey(edgeKey);
    if (!pair) return;
    const next = selectedPlanForEdit();
    const index = findEdgeIndexForPair(next.edges, pair);
    if (action === 'none') {
      if (index >= 0) {
        next.edges[index].included = false;
        if (!next.edges[index].file && !String(next.edges[index].losatFilename || '').trim()) {
          next.edges.splice(index, 1);
        }
      }
      replaceLinearComparisonPlan(next);
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
    replaceLinearComparisonPlan(next);
  };
  const setLinearComparisonCardFile = (edgeKey, file) => {
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
    return linearRecordRows.find((entry) => entry.uid === uid)?.row || fallback;
  };
  const linearLayoutTokens = computed(() => (
    linearRecordLayoutEnabled.value
      ? linearRecordPositionTokens(linearSeqs, linearRecordRows)
      : []
  ));
  const linearLosatCacheInfoByEdgeKey = computed(() => Object.fromEntries(
    (Array.isArray(losatCacheInfo.value) ? losatCacheInfo.value : [])
      .filter((entry) => String(entry?.edgeKey || ''))
      .map((entry) => [String(entry.edgeKey), entry])
  ));

  const paletteLoader = createPaletteLoader({ state });
  const linearRecordSelector = createLinearRecordSelector({
    state,
    reactive,
    recordReader: ({ inputType, primaryFile, pairedFile }) => (
      inputType === 'gff'
        ? discoverGffFastaRecords({
            gffFile: primaryFile,
            fastaFile: pairedFile,
            readText: readFileText
          })
        : discoverSequenceRecords({
            file: primaryFile,
            format: 'genbank',
            readText: readFileText
          })
    )
  });
  const getAnnotationRecordCatalog = (loadComparisonOverride = null) => {
    const inputType = mode.value === 'linear' ? lInputType.value : cInputType.value;
    const circularPrimaryFile = cInputType.value === 'gff' ? files.c_gff : files.c_gb;
    const circularPairedFile = cInputType.value === 'gff' ? files.c_fasta : null;
    const circularHasInput = Boolean(
      circularPrimaryFile && (cInputType.value !== 'gff' || circularPairedFile)
    );
    const circularIsCurrent =
      circularRecordDiscovery.inputType === cInputType.value &&
      circularRecordDiscovery.primaryFile === circularPrimaryFile &&
      circularRecordDiscovery.pairedFile === circularPairedFile;
    const loadComparison = loadComparisonOverride == null
      ? mode.value === 'linear' && hasLinearComparisonIntent.value
      : Boolean(loadComparisonOverride);
    return buildAnnotationRecordCatalog({
      mode: mode.value,
      inputType,
      loadComparison,
      multiRecordCanvas: form.multi_record_canvas,
      circularSource: {
        sourceKey: annotationSourceKey({
          scope: 'circular',
          inputType: cInputType.value,
          primaryFile: circularPrimaryFile,
          pairedFile: circularPairedFile
        }),
        hasInput: circularHasInput,
        status: circularIsCurrent
          ? circularRecordDiscovery.status
          : (circularHasInput ? 'loading' : 'idle'),
        error: circularIsCurrent ? circularRecordDiscovery.error : '',
        records: circularIsCurrent ? circularRecordList.value : []
      },
      linearSources: linearSeqs.map((seq) => {
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
    });
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
    buildCompositionIntent: () => {
      const svg = svgContainer.value?.querySelector?.('svg') || null;
      if (!svg) return null;
      try {
        return compositionUserDeltas(svg);
      } catch (_error) {
        return null;
      }
    },
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
  const history = createHistoryManager({
    buildIntent: historySnapshots.buildHistoryIntent,
    applyIntent: historySnapshots.applyHistoryIntent,
    buildCheckpoint: historySnapshots.buildArtifactCheckpoint,
    applyCheckpoint: historySnapshots.applyArtifactCheckpoint,
    signatureFor: historySnapshots.snapshotSignature,
    fileStore: historyFileStore,
    collectCurrentFileIds: historySnapshots.collectCurrentFileIds,
    makeRef: ref
  });
  window.__GBDRAW_HISTORY__ = history;
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
  const selectResult = (index) => previewRuntime.selectResult(index);

  const { handleWheel, startPan, doPan, endPan, resetPreviewViewport } = createPanZoom(state);
  const { startResizing } = createSidebarResize(state);

  const legendActions = createLegendManager({
    state,
    history,
    previewRuntime
  });
  const svgActions = createSvgStyles({
    state,
    watch,
    nextTick,
    legendActions,
    previewRuntime
  });
  const featureSelection = createFeatureSelection({ state, onMounted, onUnmounted });
  const featureActions = createFeatureEditor({
    state,
    nextTick,
    legendActions,
    svgActions,
    featureSelection,
    previewRuntime
  });
  const previewFeatureSearch = createPreviewFeatureSearch({
    state,
    watch,
    nextTick,
    computed,
    reactive,
    openFeatureEditorForFeature: featureActions.openFeatureEditorForFeature
  });

  watch(selectedResultIndex, (newIndex, oldIndex) => {
    if (newIndex !== oldIndex) previewRuntime.flushActiveResult({ markIncremental: false });
    featureSelection.clearFeatureSelection({ clearStatus: true });
  });
  watch(mode, () => {
    featureSelection.clearFeatureSelection({ clearStatus: true });
  });
  watch(svgContent, () => {
    if (!skipCaptureBaseConfig.value) {
      featureSelection.clearFeatureSelection({ clearStatus: true });
    }
  });

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
    () => [selectedFeatureRecordIdx.value, showRightDrawer.value, rightDrawerTab.value, extractedFeatures.value.length],
    resetFeatureListScroll
  );

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
  const circularTrackSlotEditor = createCircularTrackSlotEditor({ state });
  const linearTrackSlotEditor = createLinearTrackSlotEditor({ state });
  const annotationEditor = createAnnotationEditor({ state, getRecordCatalog: getAnnotationRecordCatalog });
  watch(
    () => {
      const catalog = getAnnotationRecordCatalog();
      return `${catalog.status}:${catalog.signature}`;
    },
    () => {
      const catalog = getAnnotationRecordCatalog();
      if (catalog.status === 'ready') {
        reconcileAnnotationRecordBindings(annotationSets, catalog);
      }
    },
    { immediate: true }
  );
  const circularConservationLayoutWarning = computed(() => estimateCircularConservationLayoutWarning(state));
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
  const circularSeparateStrandsDisabled = computed(() => (
    mode.value === 'circular' && Boolean(adv.resolve_overlaps)
  ));
  const circularResolveOverlapsDisabled = computed(() => (
    mode.value === 'circular' && Boolean(form.separate_strands)
  ));
  const circularSeparateStrandsOptionClass = computed(() => (
    circularSeparateStrandsDisabled.value ? disabledOptionClass : enabledOptionClass
  ));
  const circularResolveOverlapsOptionClass = computed(() => (
    circularResolveOverlapsDisabled.value ? disabledOptionClass : enabledOptionClass
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
  const rowsForDepthTrackCount = (count) => {
    const normalizedCount = Math.max(1, Number(count) || 1);
    ensureDepthTrackEditableConfigCount(normalizedCount);
    return Array.from({ length: normalizedCount }, (_, index) => ({
      index,
      key: `depth-track-${index}`,
      config: adv.depth_tracks[index] || normalizeDepthTrackConfig(null, index)
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
  const depthTrackConfigDefaults = () => ({
    labelForIndex: getDepthTrackFallbackLabel,
    colorForIndex: depthTrackFallbackColor,
    depthColor: adv.depth_color,
    depthHeight: adv.depth_height,
    largeTickInterval: null,
    smallTickInterval: null,
    tickFontSize: null
  });
  const normalizeDepthTrackConfig = (entry, index) => (
    normalizeDepthTrackConfigEntry(entry, index, depthTrackConfigDefaults())
  );
  const optionalNumberInputValue = (value) => value ?? '';
  const setOptionalNumberInputValue = (target, key, value) => {
    if (!target || typeof target !== 'object') return;
    const text = String(value ?? '').trim();
    target[key] = text === '' ? null : text;
  };
  const activeDepthTrackCount = () => {
    if (mode.value === 'linear') {
      return linearDepthTrackUiCount();
    }
    return sourceDepthTrackCount(files.c_depth, depthTrackUiCounts.circular);
  };
  const ensureDepthTrackConfigCount = (count = activeDepthTrackCount()) => {
    const targetCount = Math.max(1, Number(count) || 1);
    const normalized = ensureDepthTrackConfigCountEntries(
      adv.depth_tracks,
      targetCount,
      depthTrackConfigDefaults()
    );
    adv.depth_tracks.splice(0, adv.depth_tracks.length, ...normalized);
  };
  const ensureDepthTrackEditableConfigCount = (count = activeDepthTrackCount()) => {
    const targetCount = Math.max(1, Number(count) || 1);
    if (!Array.isArray(adv.depth_tracks)) adv.depth_tracks = [];
    ensureDepthTrackConfigShape(adv.depth_tracks, targetCount, depthTrackConfigDefaults());
  };
  const circularDepthTrackRows = computed(() => rowsForDepthTrackCount(
    sourceDepthTrackCount(files.c_depth, depthTrackUiCounts.circular)
  ));
  const linearDepthTrackRows = () => rowsForDepthTrackCount(linearDepthTrackUiCount());
  const depthTrackRows = computed(() => rowsForDepthTrackCount(activeDepthTrackCount()));
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
    { key: 'replicon', label: 'Replicon', visibilityKey: 'linear_show_replicon' },
    { key: 'accession', label: 'Accession', visibilityKey: 'linear_show_accession' },
    { key: 'length', label: 'Length / Coord.', visibilityKey: 'linear_show_length' }
  ]);
  const ensureDefinitionLineStyle = (kind) => {
    const key = String(kind || '');
    if (
      !adv.linear_definition_line_styles ||
      typeof adv.linear_definition_line_styles !== 'object' ||
      Array.isArray(adv.linear_definition_line_styles)
    ) {
      adv.linear_definition_line_styles = {};
    }
    const existing = adv.linear_definition_line_styles[key];
    if (!existing || typeof existing !== 'object' || Array.isArray(existing)) {
      adv.linear_definition_line_styles[key] = {
        font_size: null,
        font_weight: null,
        fill: null
      };
    }
    return adv.linear_definition_line_styles[key];
  };
  const getDefinitionLineStyleSize = (kind) => optionalNumberInputValue(
    ensureDefinitionLineStyle(kind).font_size
  );
  const setDefinitionLineStyleSize = (kind, value) => {
    setOptionalNumberInputValue(ensureDefinitionLineStyle(kind), 'font_size', value);
  };
  const getDefinitionLineStyleWeight = (kind) => ensureDefinitionLineStyle(kind).font_weight ?? '';
  const setDefinitionLineStyleWeight = (kind, value) => {
    const normalized = String(value || '').trim().toLowerCase();
    ensureDefinitionLineStyle(kind).font_weight = normalized === 'bold' ? 'bold' : null;
  };
  const getDefinitionLineStyleFill = (kind) => ensureDefinitionLineStyle(kind).fill ?? '';
  const setDefinitionLineStyleColor = (kind, value) => {
    const normalized = String(value || '').trim();
    ensureDefinitionLineStyle(kind).fill = normalized || null;
  };
  const getDefinitionLineStyleColorMode = (kind) => (
    colorValueMode(getDefinitionLineStyleFill(kind))
  );
  const getDefinitionLineStyleSwatchValue = (kind) => {
    return toNativeColorInputValue(getDefinitionLineStyleFill(kind));
  };
  const isDefinitionLineStyleMuted = (row) => {
    const key = row?.visibilityKey;
    return key ? adv[key] === false : false;
  };
  const normalizeDepthSlotTrackIndex = (slot) => {
    const rawTrackIndex = Number(slot?.params?.track_index);
    return Number.isInteger(rawTrackIndex) && rawTrackIndex >= 0 ? rawTrackIndex : 0;
  };
  const depthTrackConfigForIndex = (index) => {
    const idx = Math.max(0, Number(index) || 0);
    ensureDepthTrackEditableConfigCount(idx + 1);
    return adv.depth_tracks[idx];
  };
  const getDepthTrackLabel = (index) => {
    const config = depthTrackConfigForIndex(index);
    return String(config?.label ?? '');
  };
  const getDepthTrackColor = (index) => {
    const idx = Math.max(0, Number(index) || 0);
    const config = depthTrackConfigForIndex(idx);
    return String(config?.color || depthTrackFallbackColor(idx));
  };
  const setDepthTrackColor = (index, value) => {
    const idx = Math.max(0, Number(index) || 0);
    const color = String(value ?? '').trim();
    const config = depthTrackConfigForIndex(idx);
    config.color = color || depthTrackFallbackColor(idx);
  };
  const depthTrackSlotCollections = () => [
    Array.isArray(adv.circular_track_slots) ? adv.circular_track_slots : [],
    Array.isArray(adv.linear_track_slots) ? adv.linear_track_slots : []
  ];
  const syncDepthTrackSlotLabelsForTrack = (index) => {
    depthTrackSlotCollections().forEach((slots) => {
      syncDepthSlotLabels({
        slots,
        depthTracks: adv.depth_tracks,
        activeCount: adv.depth_tracks.length
      });
    });
    void index;
  };
  const setDepthTrackLabel = (index, value) => {
    const idx = Math.max(0, Number(index) || 0);
    const config = depthTrackConfigForIndex(idx);
    config.label = String(value ?? '');
    syncDepthTrackSlotLabelsForTrack(idx);
  };
  const getDepthTrackLegendLabelForSlot = (slot) => (
    getDepthTrackLabel(normalizeDepthSlotTrackIndex(slot))
  );
  const setDepthTrackLegendLabelForSlot = (slot, value) => {
    if (!slot) return;
    const idx = normalizeDepthSlotTrackIndex(slot);
    slot.params = slot.params && typeof slot.params === 'object' ? { ...slot.params } : {};
    const label = String(value ?? '');
    if (label.trim()) {
      slot.params.legend_label = label;
    } else {
      delete slot.params.legend_label;
    }
    setDepthTrackLabel(idx, label);
  };
  const syncDepthTrackSlotLabel = (slot) => {
    if (!slot || slot.renderer !== 'depth') return;
    const trackIndex = normalizeDepthSlotTrackIndex(slot);
    const hasSource = mode.value === 'linear'
      ? activeDepthTrackIndices(linearDepthRows()).includes(trackIndex)
      : Boolean(circularDepthRepresentatives()[trackIndex]);
    if (hasSource) delete slot.depth_binding_error;
    syncDepthTrackSlotLabelsForTrack(trackIndex);
  };
  const depthTrackAutoLabels = [];
  const refreshDepthTrackLabelsAfterRemoval = (previousFiles, nextFiles, removedIndex) => {
    depthFileSlotsFromValue(nextFiles).forEach((file, newIndex) => {
      const oldIndex = newIndex >= removedIndex ? newIndex + 1 : newIndex;
      const config = adv.depth_tracks[newIndex];
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
  const updateDepthTrackLabelFromFile = (index, file, previousFile = null) => {
    if (!file) return;
    const config = adv.depth_tracks[index];
    if (!config) return;
    const currentLabel = String(config.label ?? '').trim();
    if (isDepthTrackAutoLabel(currentLabel, index, previousFile) || currentLabel === depthTrackAutoLabels[index]) {
      const nextLabel = getDepthTrackLabelFromFile(file, index);
      config.label = nextLabel;
      depthTrackAutoLabels[index] = nextLabel;
      syncDepthTrackSlotLabelsForTrack(index);
    }
  };
  const getCircularDepthFile = (index) => circularDepthRepresentatives()[Number(index)] || null;
  const setCircularDepthFile = (index, file) => {
    const idx = Math.max(0, Number(index) || 0);
    ensureDepthTrackConfigCount(idx + 1);
    depthTrackUiCounts.circular = Math.max(depthTrackUiCounts.circular, idx + 1);
    const rows = circularDepthRows();
    const previousFile = circularDepthRepresentatives()[idx] || null;
    rows.forEach((row) => {
      row[idx] = file || null;
    });
    files.c_depth = rows.map((row) => compactDepthFileSlots(row));
    if (file) {
      updateDepthTrackLabelFromFile(idx, file, previousFile);
      form.show_depth = true;
    }
  };
  const getLinearDepthFile = (seq, index) => depthFileSlotsFromValue(seq?.depth)[Number(index)] || null;
  const setLinearDepthFile = (seq, index, file) => {
    if (!seq) return;
    const idx = Math.max(0, Number(index) || 0);
    const logicalWidth = Math.max(linearDepthLogicalWidth(), idx + 1);
    padLinearDepthRows(logicalWidth);
    ensureDepthTrackConfigCount(logicalWidth);
    const slots = depthFileSlotsFromValue(seq.depth);
    const previousFile = slots[idx] || null;
    if (file) {
      slots[idx] = file;
      seq.depth = slots;
    } else {
      seq.depth = clearDepthTrackSourceAt(slots, idx, logicalWidth);
    }
    if (file) {
      updateDepthTrackLabelFromFile(idx, file, previousFile);
      form.show_depth = true;
    }
  };
  const addCircularDepthTrack = () => {
    depthTrackUiCounts.circular = sourceDepthTrackCount(files.c_depth, depthTrackUiCounts.circular) + 1;
    ensureDepthTrackConfigCount(depthTrackUiCounts.circular);
    if (canShowDepthTrack.value) form.show_depth = true;
  };
  const addLinearDepthTrack = () => {
    const nextCount = linearDepthTrackUiCount() + 1;
    padLinearDepthRows(nextCount);
    ensureDepthTrackConfigCount(nextCount);
    if (canShowDepthTrack.value) form.show_depth = true;
    if (adv.linear_track_slots_enabled && form.show_depth) {
      linearTrackSlotEditor.ensureLinearTrackDepthSlots();
    }
  };
  const removeCircularDepthTrack = (index) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0) return;
    const count = sourceDepthTrackCount(files.c_depth, depthTrackUiCounts.circular);
    const previousFiles = circularDepthRepresentatives();
    files.c_depth = removeDepthTrackColumnAt(circularDepthRows(), idx)
      .map((row) => compactDepthFileSlots(row));
    if (idx < adv.depth_tracks.length) adv.depth_tracks.splice(idx, 1);
    depthTrackUiCounts.circular = count <= 1 ? 1 : Math.max(1, count - 1);
    refreshDepthTrackLabelsAfterRemoval(previousFiles, circularDepthRepresentatives(), idx);
    ensureDepthTrackConfigCount(activeDepthTrackCount());
    const activeFileCount = circularDepthRepresentatives().length;
    adv.circular_track_slots.splice(
      0,
      adv.circular_track_slots.length,
      ...reindexDepthSlots({
        slots: adv.circular_track_slots,
        removedIndex: idx,
        activeCount: activeFileCount,
        managedPredicate: isDefaultManagedDepthSlot
      })
    );
    syncDepthTrackSlotLabelsForTrack(idx);
    circularTrackSlotEditor.normalizeCircularTrackSlots();
    if (adv.circular_track_slots_enabled && form.show_depth && activeFileCount > 0) {
      circularTrackSlotEditor.ensureCircularTrackDepthSlot();
    }
  };
  const removeLinearDepthTrack = (_seq, index) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0) return;
    const logicalWidth = linearDepthLogicalWidth();
    if (idx >= logicalWidth) return;
    const nextRows = removeDepthTrackColumnAt(linearDepthRows(), idx);
    linearSeqs.forEach((seq, recordIndex) => {
      seq.depth = nextRows[recordIndex] || [];
    });
    if (idx < adv.depth_tracks.length) adv.depth_tracks.splice(idx, 1);
    depthTrackAutoLabels.splice(idx, 1);
    ensureDepthTrackConfigCount(activeDepthTrackCount());
    const activeFileCount = activeDepthTrackIndices(linearDepthRows()).length;
    const previousAxisIndex = Number(adv.linear_track_slots_axis_index);
    const removedManagedSlotCountBeforeAxis = Number.isInteger(previousAxisIndex)
      ? adv.linear_track_slots.reduce((count, slot, slotIndex) => {
          if (slotIndex >= previousAxisIndex || !isDefaultManagedDepthSlot(slot)) return count;
          return depthSlotTrackIndex(slot, slotIndex) === idx ? count + 1 : count;
        }, 0)
      : 0;
    adv.linear_track_slots.splice(
      0,
      adv.linear_track_slots.length,
      ...reindexDepthSlots({
        slots: adv.linear_track_slots,
        removedIndex: idx,
        activeCount: Math.max(0, logicalWidth - 1),
        managedPredicate: isDefaultManagedDepthSlot
      })
    );
    if (Number.isInteger(previousAxisIndex)) {
      adv.linear_track_slots_axis_index = Math.max(
        0,
        previousAxisIndex - removedManagedSlotCountBeforeAxis
      );
    }
    syncDepthTrackSlotLabelsForTrack(0);
    linearTrackSlotEditor.syncLinearDepthSlotHeightsFromDepthTracks();
    linearTrackSlotEditor.normalizeLinearTrackSlots();
    if (adv.linear_track_slots_enabled && form.show_depth && activeFileCount > 0) {
      linearTrackSlotEditor.ensureLinearTrackDepthSlots();
    }
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
      ensureDepthTrackConfigCount(activeDepthTrackCount());
    },
    { deep: true, immediate: true }
  );
  watch(
    () => [canShowDepthTrack.value, form.show_depth],
    ([available, showDepth]) => {
      if (!available && showDepth) form.show_depth = false;
    },
    { immediate: true }
  );
  watch(
    () => [mode.value, form.separate_strands, adv.resolve_overlaps],
    ([currentMode, separateStrands, resolveOverlaps]) => {
      if (currentMode === 'circular' && separateStrands && resolveOverlaps) {
        adv.resolve_overlaps = false;
      }
    },
    { immediate: true }
  );
  watch(mode, (nextMode, previousMode) => {
    if (nextMode === previousMode) return;
    if (state.semanticFileWatchersSuppressed.value) {
      state.modeProfileStateManager.invalidate(nextMode);
    } else {
      state.modeProfileStateManager.transition(adv, previousMode, nextMode);
    }
    matchSequenceRegistry?.reset?.();
    clickedPairwiseMatch.value = null;
  });
  const isCircularConservationUploadSource = () => (
    String(circularConservation.source || '').trim().toLowerCase() === 'upload'
  );
  const isDerivedCircularConservationReplay = () => (
    !isCircularConservationUploadSource() &&
    files.c_conservation_blasts_source === 'losat-cache' &&
    normalizeFileList(files.c_conservation_blasts).length > 0
  );
  const getCircularConservationSourceFiles = () => (
    isCircularConservationUploadSource() || isDerivedCircularConservationReplay()
      ? normalizeFileList(files.c_conservation_blasts)
      : normalizeFileList(files.c_conservation_fastas)
  );
  const syncCircularConservationEnabled = (sourceFiles = getCircularConservationSourceFiles()) => {
    circularConservation.enabled = normalizeFileList(sourceFiles).length > 0;
  };
  const clearDerivedCircularConservationBlasts = () => {
    if (files.c_conservation_blasts_source !== 'losat-cache') return;
    files.c_conservation_blasts = [];
    files.c_conservation_blasts_source = null;
  };
  const setCircularConservationSourceFiles = (nextFiles) => {
    const normalized = normalizeFileList(nextFiles);
    if (isCircularConservationUploadSource()) {
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
    files.c_conservation_blasts = normalizeFileList(nextFiles);
    files.c_conservation_blasts_source = null;
    files.c_conservation_sequence_sources = [];
    losatCacheInfo.value = [];
    syncCircularConservationSeries();
  };
  const syncCircularConservationSeries = () => {
    const sourceFiles = getCircularConservationSourceFiles();
    if (isDerivedCircularConservationReplay()) {
      circularConservation.enabled = true;
      if (adv.circular_track_slots_enabled === true) {
        circularTrackSlotEditor.syncCircularConservationSlots();
      }
      return;
    }
    syncCircularConservationEnabled(sourceFiles);
    const legacyLabels = parseConservationLabelText(circularConservation.labels);
    const nextSeries = reconcileConservationSeries({
      sourceFiles,
      previousSeries: circularConservation.series,
      legacyLabels
    });
    circularConservation.series.splice(0, circularConservation.series.length, ...nextSeries);
    if (adv.circular_track_slots_enabled === true) {
      circularTrackSlotEditor.syncCircularConservationSlots();
    }
  };
  const circularConservationSeriesRows = computed(() => {
    return (Array.isArray(circularConservation.series) ? circularConservation.series : []).map((entry, index) => ({
      index,
      filename: String(entry?.fileName || `source_${Number(index) + 1}`).trim(),
      sourceLabel: `${isCircularConservationUploadSource() ? 'BLAST' : 'Comparison'} ${Number(index) + 1}`,
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
    const idx = Number(index);
    const target = idx + Math.sign(Number(direction));
    return (
      Array.isArray(circularConservation.series) &&
      Number.isInteger(idx) &&
      idx >= 0 &&
      idx < circularConservation.series.length &&
      target >= 0 &&
      target < circularConservation.series.length
    );
  };
  const moveCircularConservationSeries = (index, direction) => {
    if (moveConservationSeriesEntry(circularConservation.series, index, direction)) {
      if (adv.circular_track_slots_enabled === true) {
        circularTrackSlotEditor.syncCircularConservationSlots();
      }
    }
  };
  const openCircularConservationComparisonFilePicker = () => {
    circularConservationFastaInput.value?.click();
  };
  const addCircularConservationComparisonFile = (event) => {
    const target = event?.target || null;
    const selectedFile = Array.from(target?.files || []).filter(Boolean)[0] || null;
    if (!selectedFile) return;
    clearDerivedCircularConservationBlasts();
    files.c_conservation_fastas = [...normalizeFileList(files.c_conservation_fastas), selectedFile];
    losatCacheInfo.value = [];
    syncCircularConservationSeries();
    if (target) target.value = '';
  };
  const setCircularConservationCompanionFile = (sourceIndex, event) => {
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
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= circularConservation.series.length) return;
    if (isDerivedCircularConservationReplay()) {
      const orderedBlasts = orderedConservationSources(
        files.c_conservation_blasts,
        circularConservation
      ).map((entry) => entry.file);
      const orderedFastas = orderedOptionalConservationFiles(
        files.c_conservation_fastas,
        circularConservation
      );
      const orderedSequenceSources = orderedOptionalConservationFiles(
        files.c_conservation_sequence_sources,
        circularConservation
      );
      orderedBlasts.splice(idx, 1);
      orderedFastas.splice(idx, 1);
      orderedSequenceSources.splice(idx, 1);
      circularConservation.series.splice(idx, 1);
      circularConservation.series.forEach((seriesEntry, seriesIndex) => {
        seriesEntry.sourceIndex = seriesIndex;
      });
      circularConservation.labels = circularConservation.series
        .map((seriesEntry) => seriesEntry.label)
        .join(',');
      files.c_conservation_blasts = orderedBlasts;
      files.c_conservation_fastas = orderedFastas;
      files.c_conservation_sequence_sources = orderedSequenceSources;
      files.c_conservation_blasts_source = orderedBlasts.length > 0
        ? 'losat-cache'
        : null;
      circularConservation.enabled = orderedBlasts.length > 0;
      losatCacheInfo.value = [];
      if (adv.circular_track_slots_enabled === true) {
        circularTrackSlotEditor.syncCircularConservationSlots();
      }
      return;
    }
    const entry = circularConservation.series[idx];
    const sourceFiles = getCircularConservationSourceFiles();
    const descriptors = conservationSourceDescriptors(sourceFiles);
    let sourceIndex = descriptors.findIndex((descriptor) => descriptor.sourceKey === String(entry?.sourceKey || ''));
    if (sourceIndex < 0) {
      const fileName = String(entry?.fileName || '').trim();
      sourceIndex = descriptors.findIndex((descriptor) => descriptor.fileName === fileName);
    }
    if (sourceIndex < 0 && idx < sourceFiles.length) sourceIndex = idx;
    if (sourceIndex < 0 || sourceIndex >= sourceFiles.length) return;
    if (isCircularConservationUploadSource()) {
      files.c_conservation_sequence_sources = (Array.isArray(files.c_conservation_sequence_sources)
        ? files.c_conservation_sequence_sources
        : [])
        .filter((_, fileIndex) => fileIndex !== sourceIndex);
    }
    setCircularConservationSourceFiles(sourceFiles.filter((_, fileIndex) => fileIndex !== sourceIndex));
  };
  watch(
    () => [
      circularConservation.source,
      files.c_conservation_blasts,
      files.c_conservation_blasts_source,
      files.c_conservation_fastas,
      circularConservation.labels
    ],
    syncCircularConservationSeries,
    { deep: true, immediate: true }
  );
  watch(
    () => [
      adv.circular_track_slots_enabled,
      form.show_depth,
      form.suppress_gc,
      form.suppress_skew,
      representativeDepthFiles(files.c_depth).length,
      depthTrackUiCounts.circular
    ],
    ([slotsEnabled, showDepth]) => {
      if (slotsEnabled) {
        circularTrackSlotEditor.normalizeCircularTrackSlots();
        circularTrackSlotEditor.syncCircularConservationSlots();
      }
      if (slotsEnabled && showDepth) {
        circularTrackSlotEditor.ensureCircularTrackDepthSlot();
      }
      if (slotsEnabled) {
        adv.depth_tracks.forEach((_track, index) => syncDepthTrackSlotLabelsForTrack(index));
      }
    }
  );
  watch(
    () => [
      adv.circular_track_slots_enabled,
      circularConservation.enabled,
      circularConservation.source,
      circularConservation.series.map((entry) => `${entry?.sourceKey || ''}:${entry?.label || ''}:${entry?.color || ''}`).join('|')
    ],
    ([slotsEnabled]) => {
      if (slotsEnabled) circularTrackSlotEditor.syncCircularConservationSlots();
    }
  );
  const legendLayout = createLegendLayout({ state, legendActions, history });
  legendActions.setLegendGeometryChangedHandler(legendLayout.refreshLegendGeometry);
  const {
    runAnalysis: runGeneratedDiagramAnalysis,
    cancelRunAnalysis,
    runLabelReflow,
    refreshCircularRecordOrder,
    downloadCliHelperFiles,
    downloadLosatCache,
    downloadLosatPair,
    setLosatPairFilename,
    clearLosatCache,
    getLosatPairDefaultName
  } = createRunAnalysis({
    state,
    serializeCanonicalFiles: (comparisonPlanSnapshot, linearRecordCatalog = null) => (
      serializeActiveRenderFiles(state.mode.value, state, {
        comparisonPlan: comparisonPlanSnapshot,
        linearRecordCatalog
      })
    ),
    prepareLinearRecordCatalog,
    canonicalSessionVersion: SESSION_VERSION,
    adoptCanonicalRenderArtifacts,
    buildGeneratedArtifactSnapshot: historySnapshots.buildGeneratedArtifactSnapshot,
    applyGeneratedArtifactSnapshot: historySnapshots.applyGeneratedArtifactSnapshot,
    resetPreviewViewport,
    validateAnnotationTargets: ({ loadComparison }) => {
      const catalog = getAnnotationRecordCatalog(loadComparison);
      reconcileAnnotationRecordBindings(annotationSets, catalog);
      return validateAnnotationRecordTargets(annotationSets, catalog);
    }
  });
  const resultsManager = createResultsManager({
    state,
    legendLayout,
    rerenderLinearDefinitions: runLabelReflow
  });

  setupWatchers({
    state,
    watch,
    nextTick,
    onMounted,
    legendActions,
    svgActions,
    featureActions,
    legendLayout,
    resultsManager,
    runLabelReflow,
    refreshCircularRecordOrder,
    refreshLinearRecordSelectors: linearRecordSelector.refresh,
    resetPreviewViewport,
    resetRightDrawer: rightDrawerActions.resetRightDrawer,
    previewRuntime,
    preparePaletteDefinitions: paletteLoader.loadPaletteAsset
  });

  const refreshLoadedSessionSvgLayout = async () => {
    await nextTick();
    await new Promise((resolve) => {
      if (typeof window.requestAnimationFrame === 'function') {
        window.requestAnimationFrame(() => resolve());
      } else {
        setTimeout(resolve, 0);
      }
    });

    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return;
    try {
      legendLayout.captureBaseConfig();
      legendActions.setupLegendDrag();
      legendLayout.setupDiagramDrag(false);
    } catch (error) {
      errorLog.value = {
        summary: error?.message || 'The saved SVG composition metadata is invalid.',
        details: []
      };
      throw error;
    }
  };

  const synchronizeLoadedSessionLegendEntries = async () => {
    if (!svgContent.value) return;

    let svgReady = false;
    for (let attempt = 0; attempt < 20; attempt += 1) {
      await nextTick();
      if (svgContainer.value?.querySelector?.('svg')) {
        svgReady = true;
        break;
      }
      await new Promise((resolve) => setTimeout(resolve, 50));
    }
    if (!svgReady) return;

    // The admitted, sanitized Result is the sole visual authority. Extraction
    // only carries nonvisual metadata when both its caption and color match.
    legendActions.extractLegendEntries();
  };

  const importSession = async (event) => {
    const result = await importSessionFromFile(event, {
      afterLoad: refreshLoadedSessionSvgLayout,
      rollbackState: createSessionImportRollbackState({
        depthTrackUiCounts,
        depthTracks: adv.depth_tracks,
        featureListScrollTop,
        featureListScrollRef,
        selectedPairwiseBlockOrthogroupId
      })
    });
    if (result?.status === 'ok' || result?.status === 'legacy') {
      await nextTick();
      if (result?.status === 'ok') {
        await synchronizeLoadedSessionLegendEntries();
      }
      recordSessionLifecycleEvent('history-baseline-start');
      await history.initializeIntentBaseline('Loaded session');
      recordSessionLifecycleEvent('history-baseline-end');
    }
    return result;
  };

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
    setLegendEntryStrokeColorValue,
    updateLegendEntryStrokeColor,
    updateLegendEntryStrokeWidth,
    reconcileLegendEntries,
    reconcileStrokeOverrides,
    resetLegendEntryStroke,
    resetAllStrokes
  } = legendActions;

  const {
    addCustomColor,
    addPriorityRule,
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
    addFeatureVisibilityRule,
    downloadFeatureVisibilityRulesTsv,
    featureVisibilityQualifierSuggestions,
    featureVisibilityRuleDetail,
    getFeatureColor,
    getFeatureColorValue,
    canEditFeatureColor,
    getFeatureVisibility,
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
    handleGlobalLabelModeChoice,
    requestLabelTextChangeByFeatureId,
    requestLabelTextChangeByKey,
    reconcileFeatureVisibility,
    reconcileLabelOverrides,
    resetAllLabelTextOverrides
  } = featureActions;

  historySnapshots.setAfterApplyHistoryIntent(async (_intent, { domains, changes } = {}) => {
    if (!svgContainer.value?.querySelector?.('svg')) return;
    const changedDomains = domains instanceof Set ? domains : new Set();
    if (changedDomains.has('ui')) {
      legendLayout.reconcileCompositionUserDeltas(_intent?.ui?.compositionUserDeltas);
    }
    if (changedDomains.has('config') || changedDomains.has('features')) {
      svgActions.applyPaletteToSvg();
      svgActions.applySpecificRulesToSvg();
    }
    if (changedDomains.has('features')) {
      reconcileFeatureVisibility();
      reconcileLabelOverrides();
    }
    if (changedDomains.has('editorState')) {
      reconcileLegendEntries({ restoreColorState: true });
      reconcileStrokeOverrides({ changes });
      reconcileLabelOverrides();
    }
    await nextTick();
  });

  const { updatePalette, resetColors, cancelDefinitionUpdate } = resultsManager;
  const undoableAction = (label, fn) => (...args) => history.runUndoable(label, () => fn(...args));
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
  const handleColorScopeChoiceWithHistory = undoableAction('Change feature color', handleColorScopeChoice);
  const handleFeatureStyleScopeChoiceWithHistory = (...args) => history.runUndoable(
    featureStyleScopeDialog.kind === 'stroke' ? 'Change feature stroke' : 'Change feature color',
    () => handleFeatureStyleScopeChoice(...args)
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
  const handleGlobalLabelModeChoiceWithHistory = undoableAction('Change label visibility', handleGlobalLabelModeChoice);
  const requestLabelTextChangeByFeatureIdWithHistory = undoableAction(
    'Change label text',
    requestLabelTextChangeByFeatureId
  );
  const requestLabelTextChangeByKeyWithHistory = undoableAction('Change label text', requestLabelTextChangeByKey);
  const resetAllLabelTextOverridesWithHistory = undoableAction('Reset label edits', resetAllLabelTextOverrides);
  const loadLabelOverrideTableWithHistory = undoableAction('Load label edits', loadLabelOverrideTable);
  const runInfoCopyStatus = ref('');

  const runInfoElapsedText = (info) => formatElapsedMs(info?.elapsedMs);
  const runInfoReproducibilityText = (info) => reproducibilityLabel(info?.reproducibility?.level);
  const runInfoHasCliHelperFiles = computed(() =>
    Array.isArray(lastRunInfo.value?.helperFiles) && lastRunInfo.value.helperFiles.length > 0
  );
  const copyRunCommand = async () => {
    const command = String(lastRunInfo.value?.command || '');
    if (!command) return;
    try {
      await copyTextToClipboard(command);
      runInfoCopyStatus.value = 'Copied';
      setTimeout(() => {
        if (runInfoCopyStatus.value === 'Copied') runInfoCopyStatus.value = '';
      }, 1600);
    } catch (error) {
      console.warn('Failed to copy run command:', error);
      runInfoCopyStatus.value = 'Copy failed';
      setTimeout(() => {
        if (runInfoCopyStatus.value === 'Copy failed') runInfoCopyStatus.value = '';
      }, 2200);
    }
  };

  async function prepareLinearRecordCatalog(loadComparison = false) {
    if (mode.value !== 'linear') return { catalog: null, error: '' };
    const hasAutomaticSequence = linearSeqs.some((sequence) => {
      if (String(sequence?.region_record_id || '').trim()) return false;
      const primary = lInputType.value === 'gff' ? sequence?.gff : sequence?.gb;
      return Boolean(primary && (lInputType.value !== 'gff' || sequence?.fasta));
    });
    const hasRegionAnnotations = annotationSets.some((set) => (
      Array.isArray(set?.annotations) && set.annotations.length > 0
    ));
    if (!hasAutomaticSequence && !hasRegionAnnotations) {
      return { catalog: null, error: '' };
    }
    let catalog = getAnnotationRecordCatalog(loadComparison);
    if (catalog.status !== 'ready') {
      try {
        await linearRecordSelector.refresh();
      } catch (error) {
        console.warn('Failed to start Linear record discovery:', error);
      }
      catalog = getAnnotationRecordCatalog(loadComparison);
    }
    return catalog.status === 'ready'
      ? { catalog, error: '' }
      : {
          catalog: null,
          error: catalog.issues[0] || 'Could not read records from the Linear input file(s).'
        };
  }

  const runAnalysis = async () => history.runUndoableCheckpoint('Generate diagram', async () => {
    cancelDefinitionUpdate();
    const comparisonPlanSnapshot = mode.value === 'linear'
      ? linearComparisonResolution.value
      : null;

    const result = await runGeneratedDiagramAnalysis(comparisonPlanSnapshot);
    if (result?.status === 'error' && mode.value === 'linear') {
      await focusLinearComparisonIssue();
    }
    if (result?.status === 'ok') {
      featureSelection.clearFeatureSelection({ clearStatus: true });
    }
    return result;
  }, {
    shouldCommit: (result) => result?.status === 'ok'
  });

  const cancelGeneration = () => {
    cancelRunAnalysis();
  };

  const orthogroupActions = createOrthogroupEditor({
    state,
    runAnalysis
  });
  const canUseClickedOrthogroupActions = computed(() => {
    const cf = clickedFeature.value;
    return Boolean(
      cf &&
      mode.value === 'linear' &&
      hasActiveLinearLosatIntent.value &&
      linearComparisonResolution.value.valid &&
      losatProgram.value === 'blastp' &&
      losat.blastp?.mode === 'orthogroup' &&
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

  const alignByClickedOrthogroup = async () => {
    const cf = clickedFeature.value;
    if (!cf?.orthogroupId) return;
    selectedOrthogroupAlignmentFeature.value = String(cf.orthogroupId || '').trim();
    clickedFeature.value = null;
    await runAnalysis();
  };

  const resetOrthogroupAlignment = async () => {
    clickedFeature.value = null;
    await orthogroupActions.resetOrthogroupAlignment();
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

  const openClickedOrthogroupInEditor = () => {
    const orthogroupId = String(clickedFeature.value?.orthogroupId || '').trim();
    if (!openOrthogroupInDrawer(orthogroupId)) return;
    clickedFeature.value = null;
  };

  const { resetAllPositions, resetCanvasPadding } = legendLayout;

  const resetSettings = () => {
    const proceed = window.confirm(
      'Reset all settings to the webapp defaults?\n\nUploaded files and current results will be kept.'
    );
    if (!proceed) return false;

    return history.runUndoableCheckpoint('Reset settings', async () => {
      cancelDefinitionUpdate();
      resetSettingsState(state);
      invalidateLinearComparisonArtifacts();
      matchSequenceRegistry?.reset?.();
      circularTrackNewRenderer.value = 'dinucleotide_skew';
      linearTrackNewRenderer.value = 'dinucleotide_skew';
      depthTrackUiCounts.circular = 1;
      ensureDepthTrackConfigCount(activeDepthTrackCount());
      return true;
    });
  };

  const resetLayout = () => {
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

  const getFeaturePopupConstraints = (left = clickedFeaturePos.x, top = clickedFeaturePos.y) => {
    const viewportWidth = Math.max(1, window.innerWidth || 1);
    const viewportHeight = Math.max(1, window.innerHeight || 1);
    const availableWidth = Math.max(1, viewportWidth - (FEATURE_POPUP_MARGIN * 2));
    const availableHeight = Math.max(1, viewportHeight - (FEATURE_POPUP_MARGIN * 2));
    const desiredMinWidth =
      adv.rich_feature_popup === false ? FEATURE_POPUP_SIMPLE_MIN_WIDTH : FEATURE_POPUP_RICH_MIN_WIDTH;
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
      left: `${clickedFeaturePos.x}px`
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
      left: `${clickedPairwiseMatchPos.x}px`
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

  const sessionTitleLabel = computed(() => {
    const title = normalizeSessionTitle(sessionTitle.value);
    return title || 'Untitled session';
  });

  const canUseLinearRulerOnAxis = computed(
    () => (
      form.show_scale !== false &&
      form.scale_style === 'ruler' &&
      ['above', 'below'].includes(form.linear_track_layout)
    )
  );

  const canUseCircularScaleStyling = computed(
    () => (
      adv.circular_track_slots_enabled
        ? adv.circular_track_slots.some((slot) => (
            slot?.renderer === 'ticks' &&
            circularTrackSlotEditor.circularTrackSlotEffectiveEnabled(slot)
          ))
        : form.show_scale !== false
    )
  );

  const clickedFeatureLocation = computed(() => {
    const cf = clickedFeature.value;
    if (!cf) return '';
    if (cf.location) return cf.location;
    const feat = cf.feat;
    if (!feat) return '';
    const startVal = Number(feat.start);
    const endVal = Number(feat.end);
    const startPos = Number.isFinite(startVal) ? startVal + 1 : feat.start;
    const endPos = Number.isFinite(endVal) ? endVal : feat.end;
    if (startPos === undefined || endPos === undefined || startPos === null || endPos === null) return '';
    const strand = feat.strand ? ` (${feat.strand})` : '';
    return `${startPos}..${endPos}${strand}`;
  });

  const downloadText = (filename, text, type = 'text/plain;charset=utf-8') => {
    const value = String(text ?? '');
    if (!value) return;
    downloadTextFile(String(filename || 'gbdraw.txt'), value, type);
  };

  const runExportAction = async (methodName, label) => {
    try {
      const exportService = await loadExportService();
      const exportMethod = exportService?.[methodName];
      if (typeof exportMethod !== 'function') {
        throw new Error('The export service did not provide the requested action.');
      }
      return await exportMethod();
    } catch (error) {
      const normalized = normalizeUserFacingError(error);
      errorLog.value = {
        type: 'Export error',
        message: `${label} export failed: ${normalized?.summary || 'Unknown export error.'}`,
        details: normalized?.details || []
      };
      return { status: 'error' };
    }
  };

  const downloadSVG = () => runExportAction('downloadSVG', 'SVG');
  const downloadInteractiveSVG = () => (
    runExportAction('downloadInteractiveSVG', 'Interactive SVG')
  );
  const downloadPNG = () => runExportAction('downloadPNG', 'PNG');
  const downloadPDF = () => runExportAction('downloadPDF', 'PDF');

  const specificRuleLegendOptions = computed(() => {
    const byCaption = new Map();
    for (const rule of manualSpecificRules) {
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
    const current = normalizeSessionTitle(sessionTitle.value);
    const input = prompt('Session title', current);
    if (input === null) return;
    sessionTitle.value = normalizeSessionTitle(input);
  };

  const saveSessionWithTitle = async () => {
    let title = normalizeSessionTitle(sessionTitle.value);
    if (!title) {
      const input = prompt('Session title', '');
      if (input === null) return;
      title = normalizeSessionTitle(input);
      sessionTitle.value = title;
    }
    try {
      const comparisonPlanSnapshot = mode.value === 'linear'
        ? linearComparisonResolution.value
        : null;
      const { catalog, error } = await prepareLinearRecordCatalog(
        comparisonPlanSnapshot?.hasComparisonIntent
      );
      if (error) throw new Error(error);
      return await exportSession(title, { linearRecordCatalog: catalog });
    } catch (error) {
      errorLog.value = normalizeUserFacingError(error);
      return { status: 'error' };
    }
  };

  const openFeatureEditorFromList = (feat, event) => {
    openFeatureEditorForFeature(feat, event);
  };

  const getCircularRecordOrderLabel = (selector) => {
    const normalized = String(selector || '').trim();
    if (!normalized) return '';
    const records = Array.isArray(circularRecordList.value) ? circularRecordList.value : [];
    const matched = records.find((entry) => String(entry?.selector || '').trim() === normalized);
    if (!matched) return normalized;
    return `${normalized} (${String(matched.record_id || '').trim() || 'Unknown'})`;
  };

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
    const rowValue = Number(position?.row);
    const maxRow = Math.max(1, Array.isArray(adv.multi_record_positions) ? adv.multi_record_positions.length : 1);
    if (!Number.isInteger(rowValue) || rowValue <= 0) return 1;
    return Math.min(rowValue, maxRow);
  };

  const getCircularRecordRowOptions = () => {
    const count = Array.isArray(adv.multi_record_positions) ? adv.multi_record_positions.length : 0;
    const maxRow = Math.max(1, count);
    return Array.from({ length: maxRow }, (_unused, index) => index + 1);
  };

  const sortCircularRecordPositionsByRow = () => {
    if (!Array.isArray(adv.multi_record_positions) || adv.multi_record_positions.length <= 1) return;
    const sorted = adv.multi_record_positions
      .map((entry, index) => ({ ...entry, __index: index }))
      .sort((left, right) => {
        const leftRow = Number(left.row);
        const rightRow = Number(right.row);
        if (leftRow !== rightRow) return leftRow - rightRow;
        return left.__index - right.__index;
      })
      .map(({ __index, ...entry }) => entry);
    adv.multi_record_positions.splice(0, adv.multi_record_positions.length, ...sorted);
  };

  const setCircularRecordRow = (index, rowValue) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= adv.multi_record_positions.length) return;
    const target = adv.multi_record_positions[idx];
    if (!target || typeof target !== 'object' || Array.isArray(target)) return;
    const maxRow = Math.max(1, adv.multi_record_positions.length);
    const parsedRow = Number(rowValue);
    const normalizedRow = Number.isInteger(parsedRow) && parsedRow > 0 ? Math.min(parsedRow, maxRow) : 1;
    target.row = normalizedRow;
    sortCircularRecordPositionsByRow();
  };

  const canMoveCircularRecordOrderUp = (index) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx <= 0 || idx >= adv.multi_record_positions.length) return false;
    const currentRow = getCircularRecordRow(adv.multi_record_positions[idx]);
    const prevRow = getCircularRecordRow(adv.multi_record_positions[idx - 1]);
    return currentRow === prevRow;
  };

  const canMoveCircularRecordOrderDown = (index) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= adv.multi_record_positions.length - 1) return false;
    const currentRow = getCircularRecordRow(adv.multi_record_positions[idx]);
    const nextRow = getCircularRecordRow(adv.multi_record_positions[idx + 1]);
    return currentRow === nextRow;
  };

  const moveCircularRecordOrderUp = (index) => {
    const idx = Number(index);
    if (!canMoveCircularRecordOrderUp(idx)) return;
    const next = [...adv.multi_record_positions];
    const temp = next[idx - 1];
    next[idx - 1] = next[idx];
    next[idx] = temp;
    adv.multi_record_positions.splice(0, adv.multi_record_positions.length, ...next);
  };

  const moveCircularRecordOrderDown = (index) => {
    const idx = Number(index);
    if (!canMoveCircularRecordOrderDown(idx)) return;
    const next = [...adv.multi_record_positions];
    const temp = next[idx + 1];
    next[idx + 1] = next[idx];
    next[idx] = temp;
    adv.multi_record_positions.splice(0, adv.multi_record_positions.length, ...next);
  };

  const resetCircularRecordOrder = () => {
    const defaults = buildDefaultCircularRecordPositions();
    adv.multi_record_positions.splice(0, adv.multi_record_positions.length, ...defaults);
  };

  const applyLinearSeqMutation = (items, { preserveLosatCacheInfo = false } = {}) => {
    const depthWidth = linearDepthLogicalWidth();
    const next = normalizeLinearSeqList(items);
    if (depthWidth > 0) {
      next.forEach((seq) => {
        seq.depth = padDepthFileSlots(seq.depth, depthWidth);
      });
    }
    linearSeqs.splice(0, linearSeqs.length, ...next);
    const nextRows = reconcileLinearRecordLayout(linearSeqs, linearRecordRows);
    linearRecordRows.splice(0, linearRecordRows.length, ...nextRows);
    replaceLinearComparisonPlan(
      reconcileLinearComparisonPlan(linearComparisonPlan, linearSeqs),
      { invalidate: false }
    );
    invalidateLinearComparisonArtifacts({ preserveLosatCacheInfo });
    linearReorderNotice.value = '';
  };

  const addLinearSeq = () => {
    applyLinearSeqMutation([...linearSeqs, createLinearSeq()]);
  };

  const removeLinearSeqAt = (index) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= linearSeqs.length) return;
    const current = Array.from(linearSeqs);
    const next = current.filter((_, currentIndex) => currentIndex !== idx);
    applyLinearSeqMutation(next);
  };

  const removeLastLinearSeq = () => {
    if (linearSeqs.length <= 1) return;
    removeLinearSeqAt(linearSeqs.length - 1);
  };

  const setLinearSeqPrimaryFile = (index, field, value) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= linearSeqs.length) return;
    if (!['gb', 'gff', 'fasta'].includes(field)) return;

    const nextValue = value ?? null;
    const seq = linearSeqs[idx];

    if (field === 'gb') {
      if (!nextValue) {
        removeLinearSeqAt(idx);
        return;
      }
      seq.gb = nextValue;
      invalidateLinearComparisonArtifacts();
      linearReorderNotice.value = '';
      return;
    }

    const otherField = field === 'gff' ? 'fasta' : 'gff';
    if (!nextValue && !seq[otherField]) {
      removeLinearSeqAt(idx);
      return;
    }

    seq[field] = nextValue;
    invalidateLinearComparisonArtifacts();
    linearReorderNotice.value = '';
  };

  const canMoveLinearSeqUp = (index) => {
    const idx = Number(index);
    return Number.isInteger(idx) && idx > 0 && idx < linearSeqs.length;
  };

  const canMoveLinearSeqDown = (index) => {
    const idx = Number(index);
    return Number.isInteger(idx) && idx >= 0 && idx < linearSeqs.length - 1;
  };

  const reorderLinearSeqs = (fromIndex, toIndex) => {
    const from = Number(fromIndex);
    const to = Number(toIndex);
    if (!Number.isInteger(from) || !Number.isInteger(to)) return;
    if (from < 0 || to < 0 || from >= linearSeqs.length || to >= linearSeqs.length || from === to) return;

    const current = Array.from(linearSeqs);
    const [moved] = current.splice(from, 1);
    current.splice(to, 0, moved);
    applyLinearSeqMutation(current, { preserveLosatCacheInfo: true });
  };

  const moveLinearSeqUp = (index) => {
    if (!canMoveLinearSeqUp(index)) return;
    reorderLinearSeqs(index, Number(index) - 1);
  };

  const moveLinearSeqDown = (index) => {
    if (!canMoveLinearSeqDown(index)) return;
    reorderLinearSeqs(index, Number(index) + 1);
  };

  return {
    processing,
    processingStatus,
    generationCancelRequested,
    errorLog,
    errorDisplay,
    sessionTitle,
    sessionTitleLabel,
    results,
    selectedResultIndex,
    selectResult,
    resultPanelTab,
    lastRunInfo,
    runInfoCopyStatus,
    svgContent,
    zoom,
    layoutRepositionMode,
    isPanning,
    handleWheel,
    canvasPan,
    canvasContainerRef,
    startPan,
    doPan,
    endPan,
    sidebarWidth,
    startResizing,
    mode,
    cInputType,
    lInputType,
    losatProgram,
    files,
    circularConservation,
    annotationSets,
    selectedAnnotation,
    addAnnotationSet: annotationEditor.addAnnotationSet,
    renameAnnotationSet: annotationEditor.renameAnnotationSet,
    duplicateAnnotationSet: annotationEditor.duplicateAnnotationSet,
    removeAnnotationSet: annotationEditor.removeAnnotationSet,
    addCoordinateAnnotation: annotationEditor.addCoordinateAnnotation,
    addSelectedFeatureAnnotations: annotationEditor.addSelectedFeatures,
    removeAnnotation: annotationEditor.removeAnnotation,
    setAnnotationTargetKind: annotationEditor.setAnnotationTargetKind,
    importAnnotationTableFile: annotationEditor.importAnnotationTableFile,
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
    linearDepthTrackCoverageLabel,
    linearDepthTrackIndexOptions,
    hasCircularDepthFiles,
    hasLinearDepthFiles,
    canShowDepthTrack,
    depthToggleOptionClass,
    circularSeparateStrandsDisabled,
    circularResolveOverlapsDisabled,
    circularSeparateStrandsOptionClass,
    circularResolveOverlapsOptionClass,
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
    linearSeqs,
    linearRecordLayoutEnabled,
    linearRecordGap,
    linearRecordRows,
    linearComparisonPlan,
    linearComparisonResolution,
    linearComparisonGlobalAction,
    linearComparisonUi,
    hasLinearComparisonIntent,
    hasActiveLinearLosatIntent,
    hasActiveLinearUploadIntent,
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
    setLinearSeqPrimaryFile,
    canMoveLinearSeqUp,
    canMoveLinearSeqDown,
    moveLinearSeqUp,
    moveLinearSeqDown,
    linearRecordOptions: linearRecordSelector.optionsFor,
    linearRecordSelectorDisabled: linearRecordSelector.isDisabled,
    linearRecordSelectorError: linearRecordSelector.errorFor,
    linearRecordSelectorWarning: linearRecordSelector.warningFor,
    form,
    adv,
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
    circularTrackRendererLabel: circularTrackSlotEditor.circularTrackRendererLabel,
    resetCircularTrackSlotsFromSimpleControls: circularTrackSlotEditor.resetCircularTrackSlotsFromSimpleControls,
    resetCircularTrackSlotsToPreset: circularTrackSlotEditor.resetCircularTrackSlotsToPreset,
    applyCircularTrackPreset: circularTrackSlotEditor.applyCircularTrackPreset,
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
    ensureLinearTrackDepthSlots: linearTrackSlotEditor.ensureLinearTrackDepthSlots,
    setLinearTrackSlotsEnabled: linearTrackSlotEditor.setLinearTrackSlotsEnabled,
    addLinearTrackSlot: linearTrackSlotEditor.addLinearTrackSlot,
    canAddLinearTrackRenderer: linearTrackSlotEditor.canAddLinearTrackRenderer,
    duplicateLinearTrackSlot: linearTrackSlotEditor.duplicateLinearTrackSlot,
    canDuplicateLinearTrackSlot: linearTrackSlotEditor.canDuplicateLinearTrackSlot,
    removeLinearTrackSlot: linearTrackSlotEditor.removeLinearTrackSlot,
    moveLinearTrackSlot: linearTrackSlotEditor.moveLinearTrackSlot,
    canMoveLinearTrackSlot: linearTrackSlotEditor.canMoveLinearTrackSlot,
    moveLinearTrackSlotAbove: linearTrackSlotEditor.moveLinearTrackSlotAbove,
    moveLinearTrackSlotBelow: linearTrackSlotEditor.moveLinearTrackSlotBelow,
    moveLinearTrackSlotToAxis: linearTrackSlotEditor.moveLinearTrackSlotToAxis,
    moveLinearTrackSlotToPlacement: linearTrackSlotEditor.moveLinearTrackSlotToPlacement,
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
    losat,
    ...losatSettings,
    losatCacheInfo,
    losatThreadingStatus,
    orthogroups,
    featureOrthogroupIndex,
    selectedOrthogroupAlignmentFeature,
    orthogroupNameOverrides,
    orthogroupDescriptionOverrides,
    selectedOrthogroupId,
    orthogroupSearch,
    orthogroupSortMode,
    showRightDrawer,
    rightDrawerTab,
    orthogroupCount: orthogroupActions.orthogroupCount,
    selectedAlignmentTargetLabel: orthogroupActions.selectedAlignmentTargetLabel,
    filteredOrthogroups: orthogroupActions.filteredOrthogroups,
    selectedOrthogroup: orthogroupActions.selectedOrthogroup,
    selectedOrthogroupMembersByRecord: orthogroupActions.selectedOrthogroupMembersByRecord,
    resolveOrthogroupName: orthogroupActions.resolveOrthogroupName,
    resolveOrthogroupDescription: orthogroupActions.resolveOrthogroupDescription,
    orthogroupScope: orthogroupActions.orthogroupScope,
    orthogroupScopeLabel: orthogroupActions.orthogroupScopeLabel,
    isOrthogroupRenamed: orthogroupActions.isOrthogroupRenamed,
    getOrthogroupSequenceCount: orthogroupActions.getOrthogroupSequenceCount,
    hasOrthogroupSequence: orthogroupActions.hasOrthogroupSequence,
    hasOrthogroupMemberSequence: orthogroupActions.hasOrthogroupMemberSequence,
    copyOrthogroupSequences: orthogroupActions.copyOrthogroupSequences,
    downloadOrthogroupSequences: orthogroupActions.downloadOrthogroupSequences,
    copyOrthogroupMemberSequence: orthogroupActions.copyOrthogroupMemberSequence,
    downloadOrthogroupMemberSequence: orthogroupActions.downloadOrthogroupMemberSequence,
    selectOrthogroup: orthogroupActions.selectOrthogroup,
    setOrthogroupNameOverride: orthogroupActions.setOrthogroupNameOverride,
    setOrthogroupDescriptionOverride: orthogroupActions.setOrthogroupDescriptionOverride,
    resetOrthogroupRename: orthogroupActions.resetOrthogroupRename,
    highlightOrthogroupById: orthogroupActions.highlightOrthogroupById,
    alignOrthogroupById: orthogroupActions.alignOrthogroupById,
    isRightDrawerTabAvailable: rightDrawerActions.isRightDrawerTabAvailable,
    openRightDrawerTab: rightDrawerActions.openRightDrawerTab,
    toggleRightDrawer: rightDrawerActions.toggleRightDrawer,
    closeRightDrawer: rightDrawerActions.closeRightDrawer,
    openOrthogroupInDrawer,
    circularRecordList,
    paletteDefinitions,
    paletteNames,
    selectedPalette,
    currentColors,
    paletteInstantPreviewEnabled,
    appliedPaletteName,
    appliedPaletteColors,
    pendingPaletteName,
    pendingPaletteColors,
    hasPendingPaletteDraft,
    updatePalette,
    resetColors,
    downloadLosatCache,
    downloadLosatPair,
    setLosatPairFilename,
    clearLosatCache,
    getLosatPairDefaultName,
    refreshCircularRecordOrder,
    getCircularRecordOrderLabel,
    getCircularRecordRow,
    getCircularRecordRowOptions,
    setCircularRecordRow,
    canMoveCircularRecordOrderUp,
    canMoveCircularRecordOrderDown,
    moveCircularRecordOrderUp,
    moveCircularRecordOrderDown,
    resetCircularRecordOrder,
    filterMode,
    manualBlacklist,
    manualWhitelist,
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
    manualSpecificRules,
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
    previewFeatureSearchStyle: previewFeatureSearch.previewFeatureSearchStyle,
    startPreviewFeatureSearchDrag: previewFeatureSearch.startDrag,
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
    featureListTopSpacerPx,
    featureListBottomSpacerPx,
    isFeatureDrawerMounted,
    featureListScrollRef,
    handleFeatureListScroll,
    labelSearch,
    editableLabels,
    filteredEditableLabels,
    labelTextFeatureOverrides,
    labelTextBulkOverrides,
    labelTextFeatureOverrideSources,
    labelVisibilityOverrides,
    labelOverrideContextKey,
    labelOverrideBuildWarning,
    autoLabelReflowEnabled,
    labelReflowProcessing,
    labelReflowLastError,
    filteredFeatures,
    featureColorOverrides,
    featureVisibilityManualRules,
    featureVisibilityRules,
    featureVisibilityOverrides,
    featureVisibilitySelectorCache,
    featureStrokeOverrides,
    labelLayoutDirtyReason,
    addFeatureVisibilityRule: addFeatureVisibilityRuleWithHistory,
    downloadFeatureVisibilityRulesTsv,
    featureVisibilityFeatureSuggestions,
    featureVisibilityQualifierSuggestions,
    featureVisibilityRuleDetail,
    getFeatureColor,
    getFeatureColorValue,
    getFeatureVisibility,
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
    clickedFeatureLocation,
      copyText: copyTextToClipboard,
    downloadText,
    canUseClickedOrthogroupActions,
    clickedOrthogroupDetail,
    alignByClickedOrthogroup,
    highlightClickedOrthogroup,
    clearOrthogroupHighlight,
    resetOrthogroupAlignment,
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
    globalLabelModeDialog,
    updateClickedFeatureLabelText: updateClickedFeatureLabelTextWithHistory,
    handleLabelTextScopeChoice: handleLabelTextScopeChoiceWithHistory,
    handleGlobalLabelModeChoice: handleGlobalLabelModeChoiceWithHistory,
    requestLabelTextChangeByFeatureId: requestLabelTextChangeByFeatureIdWithHistory,
    requestLabelTextChangeByKey: requestLabelTextChangeByKeyWithHistory,
    resetAllLabelTextOverrides: resetAllLabelTextOverridesWithHistory,
    downloadLabelOverrideTable,
    loadLabelOverrideTable: loadLabelOverrideTableWithHistory,
    syncLabelEditor,
    openFeatureEditorFromList,
    legendEntries,
    newLegendCaption,
    newLegendColor,
    updateLegendEntryColor,
    renameLegendEntry,
    deleteLegendEntry,
    addNewLegendEntry,
    moveLegendEntryUp,
    moveLegendEntryDown,
    sortLegendEntries,
    sortLegendEntriesByDefault,
    resetLegendPosition,
    getLegendEntryStrokeColor,
    getLegendEntryStrokeWidth,
    setLegendEntryStrokeColorValue: setLegendEntryStrokeColorValueWithHistory,
    updateLegendEntryStrokeColor,
    updateLegendEntryStrokeWidth,
    resetLegendEntryStroke,
    resetAllStrokes,
    resetAllPositions,
    resetLayout,
    canvasPadding,
    showCanvasControls,
    resetCanvasPadding,
    definitionLineStyleRows,
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
    cancelGeneration,
    downloadSVG,
    downloadInteractiveSVG,
    downloadPNG,
    downloadPDF,
    copyRunCommand,
    downloadCliHelperFiles,
    runInfoElapsedText,
    runInfoReproducibilityText,
    runInfoHasCliHelperFiles,
    resetSettings,
    saveSessionWithTitle,
    editSessionTitle,
    importSession,
    canUndoHistory,
    canRedoHistory,
    undoHistoryTitle,
    redoHistoryTitle,
    undoHistory,
    redoHistory,
    manualPriorityRules,
    newPriorityRule,
    addPriorityRule
  };
};
