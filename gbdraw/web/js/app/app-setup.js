import { state, createLinearSeq, normalizeLinearSeqList, reconcileLinearSeqPairData } from '../state.js';
import { debugLog } from '../config.js';
import { downloadSVG, downloadInteractiveSVG, downloadPNG, downloadPDF } from '../services/export.js';
import {
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
  serializeResults,
  setPreviewRuntime
} from '../services/config.js';
import { createHistoryManager } from '../services/history.js';
import { createHistoryFileStore } from '../services/history-files.js';
import { createHistorySnapshotService } from '../services/history-snapshot.js';
import { serializeCleanSvg } from '../services/svg-serialization.js';
import { downloadTextFile } from '../services/text-download.js';
import { resetLayoutState, resetSettings as resetSettingsState } from '../services/reset.js';
import {
  disposeDiagramGenerationWorker,
  preinitializeDiagramGenerationWorker
} from '../services/diagram-generation.js';
import { createPanZoom, createSidebarResize, setupGlobalUiEvents } from './ui.js';
import { createFeatureEditor } from './feature-editor.js';
import { PAIRWISE_MATCH_SELECTOR } from './pairwise-match-popup.js';
import { createFeatureSelection } from './feature-selection.js';
import { createPreviewFeatureSearch } from './feature-search/preview-actions.js';
import { createSvgStyles } from './svg-styles.js';
import { createLegendManager } from './legend.js';
import { createPyodideManager } from './pyodide.js';
import { createRunAnalysis } from './run-analysis.js';
import { formatElapsedMs, reproducibilityLabel } from './run-info.js';
import { createLegendLayout } from './legend-layout.js';
import { createResultsManager } from './results.js';
import { setupWatchers } from './watchers.js';
import { setupHistoryInputs } from './history-inputs.js';
import { setupHistoryShortcuts } from './history-shortcuts.js';
import { createPreviewRuntime } from './preview-runtime.js';
import { createOrthogroupEditor } from './orthogroups.js';
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
  addLinearComparison as appendLinearComparison,
  adjacentRowPairs,
  hasLinearComparisonIntent,
  reconcileLinearComparisons
} from './linear-comparisons.js';
import { discoverGffFastaRecords, discoverSequenceRecords } from './record-discovery.js';
import {
  conservationSourceDescriptors,
  defaultConservationSeriesLabel,
  moveConservationSeriesEntry,
  normalizeFileList,
  parseConservationLabelText,
  reconcileConservationSeries
} from './conservation-series.js';
import {
  getDepthTrackFallbackLabel,
  getDepthTrackLabelFromFile,
  isDepthTrackAutoLabel
} from './depth-tracks.js';
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

export const createAppSetup = () => {
  const {
    pyodideReady,
    diagramGenerationWorkerReady,
    diagramGenerationWorkerStatus,
    diagramGenerationWorkerError,
    processing,
    processingStatus,
    generationCancelRequested,
    loadingStatus,
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
    circularLegendPosition,
    linearLegendPosition,
    cInputType,
    lInputType,
    blastSource,
    losatProgram,
    files,
    circularConservation,
    annotationSets,
    selectedAnnotation,
    linearSeqs,
    linearRecordLayoutEnabled,
    linearRecordGap,
    linearRecordRows,
    linearComparisons,
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
    showFeaturePanel,
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
    colorScopeDialog,
    featureVisibilityScopeDialog,
    legendRenameDialog,
    resetColorDialog,
    labelTextScopeDialog,
    globalLabelModeDialog,
    sidebarWidth,
    isResizing,
    showLegendPanel,
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
    circularBaseConfig,
    linearBaseConfig,
    diagramElementBaseTransforms,
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

  const syncLinearRecordLayout = () => {
    const next = reconcileLinearRecordLayout(linearSeqs, linearRecordRows);
    linearRecordRows.splice(0, linearRecordRows.length, ...next);
    const nextComparisons = reconcileLinearComparisons(linearSeqs, linearComparisons);
    linearComparisons.splice(0, linearComparisons.length, ...nextComparisons);
  };
  const setLinearRecordRow = (uid, row) => {
    syncLinearRecordLayout();
    updateLinearRecordRow(linearRecordRows, uid, row);
  };
  const moveLinearRecordWithinRow = (uid, direction) => {
    const next = moveLinearRecordInRow(linearSeqs, linearRecordRows, uid, direction);
    linearRecordRows.splice(0, linearRecordRows.length, ...next);
  };
  const addLinearComparison = () => {
    if (linearSeqs.length < 2) return;
    syncLinearRecordLayout();
    appendLinearComparison(linearComparisons, linearSeqs[0].uid, linearSeqs[1].uid);
  };
  const removeLinearComparison = (index) => {
    if (Number.isInteger(index) && index >= 0) linearComparisons.splice(index, 1);
  };
  const setLinearComparisonFile = (index, file) => {
    if (linearComparisons[index]) linearComparisons[index].file = file || null;
  };
  const addLinearComparisonBatch = (allPairs = false) => {
    syncLinearRecordLayout();
    adjacentRowPairs(linearSeqs, linearRecordRows, allPairs).forEach(([queryUid, subjectUid]) => {
      appendLinearComparison(linearComparisons, queryUid, subjectUid);
    });
  };
  const linearRecordRowFor = (uid, fallback) => {
    syncLinearRecordLayout();
    return linearRecordRows.find((entry) => entry.uid === uid)?.row || fallback;
  };
  const linearLayoutTokens = computed(() => (
    linearRecordLayoutEnabled.value
      ? linearRecordPositionTokens(linearSeqs, linearRecordRows)
      : []
  ));

  const pyodideManager = createPyodideManager({ state });
  const getPyodide = pyodideManager.getPyodide;
  const linearRecordSelector = createLinearRecordSelector({
    state,
    reactive,
    recordReader: ({ inputType, primaryFile, pairedFile, temporaryPathPrefix }) => (
      inputType === 'gff'
        ? discoverGffFastaRecords({
            gffFile: primaryFile,
            fastaFile: pairedFile,
            pyodide: getPyodide(),
            writeFileToFs: pyodideManager.writeFileToFs,
            gffTemporaryPath: `${temporaryPathPrefix}.gff`,
            fastaTemporaryPath: `${temporaryPathPrefix}.fasta`
          })
        : discoverSequenceRecords({
            file: primaryFile,
            format: 'genbank',
            pyodide: getPyodide(),
            writeFileToFs: pyodideManager.writeFileToFs,
            temporaryPath: `${temporaryPathPrefix}.gb`
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
      ? hasLinearComparisonIntent({
          layoutEnabled: linearRecordLayoutEnabled.value,
          comparisons: linearComparisons,
          sequences: linearSeqs,
          blastSource: blastSource.value
        })
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
    buildSnapshot: historySnapshots.buildHistorySnapshot,
    applySnapshot: historySnapshots.applyHistorySnapshot,
    snapshotSignature: historySnapshots.snapshotSignature,
    fileStore: historyFileStore,
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

  const legendActions = createLegendManager({ state, getPyodide, debugLog, history });
  const svgActions = createSvgStyles({ state, watch, legendActions });
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
  setupGlobalUiEvents({ state, onMounted, onUnmounted });
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
  const linearTrackSlotsPanelOpen = ref(false);
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
  const getDefinitionLineStyleSwatchValue = (kind) => {
    const fill = getDefinitionLineStyleFill(kind);
    return /^#[0-9a-fA-F]{6}$/.test(fill) ? fill : '#000000';
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
    matchSequenceRegistry?.reset?.();
    clickedPairwiseMatch.value = null;
  });
  const isCircularConservationUploadSource = () => (
    String(circularConservation.source || '').trim().toLowerCase() === 'upload'
  );
  const getCircularConservationSourceFiles = () => (
    isCircularConservationUploadSource()
      ? normalizeFileList(files.c_conservation_blasts)
      : normalizeFileList(files.c_conservation_fastas)
  );
  const syncCircularConservationEnabled = (sourceFiles = getCircularConservationSourceFiles()) => {
    circularConservation.enabled = normalizeFileList(sourceFiles).length > 0;
  };
  const setCircularConservationSourceFiles = (nextFiles) => {
    const normalized = normalizeFileList(nextFiles);
    if (isCircularConservationUploadSource()) {
      files.c_conservation_blasts = normalized;
    } else {
      files.c_conservation_fastas = normalized;
    }
    losatCacheInfo.value = [];
    syncCircularConservationSeries();
  };
  const setCircularConservationUploadFiles = (nextFiles) => {
    files.c_conservation_blasts = normalizeFileList(nextFiles);
    files.c_conservation_sequence_sources = [];
    losatCacheInfo.value = [];
    syncCircularConservationSeries();
  };
  const syncCircularConservationSeries = () => {
    const sourceFiles = getCircularConservationSourceFiles();
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
  const legendLayout = createLegendLayout({ state, debugLog, legendActions, svgActions, history });
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
    getPyodide,
    writeFileToFs: pyodideManager.writeFileToFs,
    refreshFeatureOverrides: featureActions.refreshFeatureOverrides,
    resetPreviewViewport,
    validateAnnotationTargets: ({ loadComparison }) => {
      const catalog = getAnnotationRecordCatalog(loadComparison);
      reconcileAnnotationRecordBindings(annotationSets, catalog);
      return validateAnnotationRecordTargets(annotationSets, catalog);
    }
  });
  const resultsManager = createResultsManager({
    state,
    getPyodide,
    legendLayout,
    rerenderLinearDefinitions: runLabelReflow
  });

  setupWatchers({
    state,
    watch,
    nextTick,
    onMounted,
    debugLog,
    pyodideManager,
    legendActions,
    svgActions,
    featureActions,
    legendLayout,
    resultsManager,
    runLabelReflow,
    refreshCircularRecordOrder,
    refreshLinearRecordSelectors: linearRecordSelector.refresh,
    resetPreviewViewport,
    previewRuntime,
    prepareDiagramGenerationWorker: async () => {
      diagramGenerationWorkerReady.value = false;
      diagramGenerationWorkerError.value = null;
      diagramGenerationWorkerStatus.value = 'Preparing diagram engine...';
      try {
        await preinitializeDiagramGenerationWorker();
        diagramGenerationWorkerReady.value = true;
        diagramGenerationWorkerStatus.value = 'Diagram engine ready.';
      } catch (error) {
        const message = error instanceof Error ? error.message : String(error);
        diagramGenerationWorkerError.value = message;
        diagramGenerationWorkerStatus.value = `Diagram engine failed to start: ${message}`;
        console.error(error);
      }
    }
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

    if (mode.value !== 'linear') return;
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return;
    const legendGroup = svg.getElementById?.('legend');
    if (!legendGroup || legendGroup.getAttribute('display') === 'none') return;

    const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
    const verticalLegend = legendGroup.querySelector('#legend_vertical');
    if (horizontalLegend && verticalLegend) {
      legendActions.reflowDualLegendLayout(svg);
    } else {
      legendActions.updatePairwiseLegendPositions(svg);
    }
    legendActions.recenterCurrentLegendRoot(svg);

    const idx = selectedResultIndex.value;
    if (idx >= 0 && results.value.length > idx) {
      results.value[idx] = { ...results.value[idx], content: serializeCleanSvg(svg) };
    }
  };

  const restoreLoadedSessionLegendEntries = async (editorState) => {
    const entries = Array.isArray(editorState?.legend?.entries)
      ? editorState.legend.entries
          .map((entry) => ({
            ...entry,
            caption: String(entry?.caption || '').trim(),
            color: String(entry?.color || '').trim()
          }))
          .filter((entry) => entry.caption && entry.color)
      : [];
    if (!entries.length || !svgContent.value) return;

    let pyodideReadyForLegend = pyodideReady.value;
    for (let attempt = 0; !pyodideReadyForLegend && attempt < 150; attempt += 1) {
      await new Promise((resolve) => setTimeout(resolve, 100));
      pyodideReadyForLegend = pyodideReady.value;
    }
    if (!pyodideReadyForLegend) return;

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

    legendActions.extractLegendEntries();

    const expectedCaptions = new Set(entries.map((entry) => entry.caption));
    for (const entry of [...legendEntries.value]) {
      const caption = String(entry?.caption || '').trim();
      if (caption && !expectedCaptions.has(caption)) {
        legendActions.removeLegendEntry(caption);
      }
    }

    for (const entry of entries) {
      await legendActions.addLegendEntry(entry.caption, entry.color);
    }

    legendEntries.value = entries.map((entry) => ({
      ...entry,
      showStroke: Boolean(entry.showStroke),
      featureIds: Array.isArray(entry.featureIds) ? entry.featureIds : []
    }));
  };

  const importSession = async (event) => {
    const result = await importSessionFromFile(event, { afterLoad: refreshLoadedSessionSvgLayout });
    if (result?.status === 'ok' || result?.status === 'legacy') {
      await nextTick();
      if (result?.status === 'ok') {
        await restoreLoadedSessionLegendEntries(result.data?.editorState);
      }
      await history.captureBaseline('Loaded session');
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
    updateLegendEntryStrokeColor,
    updateLegendEntryStrokeWidth,
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
    updateClickedFeatureColor,
    handleColorScopeChoice,
    handleLegendNameCommit,
    handleLegendRenameChoice,
    selectLegendNameOption,
    renameLegendEntry,
    handleResetColorChoice,
    resetClickedFeatureFillColor,
    updateClickedFeatureStroke,
    resetClickedFeatureStroke,
    applyStrokeToAllSiblings,
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
    resetAllLabelTextOverrides
  } = featureActions;

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
  const updateClickedFeatureColorWithHistory = undoableAction('Change feature color', updateClickedFeatureColor);
  const handleColorScopeChoiceWithHistory = undoableAction('Change feature color', handleColorScopeChoice);
  const handleLegendNameCommitWithHistory = undoableAction('Rename legend item', handleLegendNameCommit);
  const handleLegendRenameChoiceWithHistory = undoableAction('Rename legend item', handleLegendRenameChoice);
  const handleResetColorChoiceWithHistory = undoableAction('Reset feature color', handleResetColorChoice);
  const resetClickedFeatureFillColorWithHistory = undoableAction('Reset feature color', resetClickedFeatureFillColor);
  const updateClickedFeatureStrokeWithHistory = undoableAction('Change feature stroke', updateClickedFeatureStroke);
  const resetClickedFeatureStrokeWithHistory = undoableAction('Reset feature stroke', resetClickedFeatureStroke);
  const applyStrokeToAllSiblingsWithHistory = undoableAction('Change feature stroke', applyStrokeToAllSiblings);
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
  const copyTextFallback = (text) => {
    const textarea = document.createElement('textarea');
    textarea.value = text;
    textarea.setAttribute('readonly', '');
    textarea.style.position = 'fixed';
    textarea.style.left = '-9999px';
    document.body.appendChild(textarea);
    textarea.select();
    try {
      document.execCommand('copy');
    } finally {
      textarea.remove();
    }
  };
  const copyRunCommand = async () => {
    const command = String(lastRunInfo.value?.command || '');
    if (!command) return;
    try {
      if (navigator.clipboard?.writeText) {
        await navigator.clipboard.writeText(command);
      } else {
        copyTextFallback(command);
      }
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

  const runAnalysis = async () => history.runUndoable('Generate diagram', async () => {
    cancelDefinitionUpdate();
    if (!pyodideReady.value) {
      if (processing.value) return { status: 'skipped' };
      generationCancelRequested.value = false;
      processing.value = true;
      processingStatus.value = loadingStatus.value || 'Preparing Python environment...';
      await pyodideManager.initPyodide();
      if (generationCancelRequested.value || !pyodideReady.value) {
        const wasCanceled = generationCancelRequested.value;
        processing.value = false;
        processingStatus.value = '';
        generationCancelRequested.value = false;
        return { status: wasCanceled ? 'canceled' : 'error' };
      }
    }

    const hasRegionAnnotations = annotationSets.some((set) => (
      Array.isArray(set?.annotations) && set.annotations.length > 0
    ));
    if (hasRegionAnnotations) {
      let catalog = getAnnotationRecordCatalog();
      if (catalog.status !== 'ready') {
        if (mode.value === 'linear') await linearRecordSelector.refresh();
        else await refreshCircularRecordOrder();
        catalog = getAnnotationRecordCatalog();
      }
      reconcileAnnotationRecordBindings(annotationSets, catalog);
      const annotationTargetError = validateAnnotationRecordTargets(annotationSets, catalog);
      if (annotationTargetError) {
        errorLog.value = { summary: annotationTargetError, details: [] };
        processing.value = false;
        processingStatus.value = '';
        return { status: 'error' };
      }
    }

    if (diagramGenerationWorkerError.value) {
      diagramGenerationWorkerReady.value = false;
      diagramGenerationWorkerError.value = null;
      diagramGenerationWorkerStatus.value = 'Preparing diagram engine...';
    }

    featureSelection.clearFeatureSelection({ clearStatus: true });
    const result = await runGeneratedDiagramAnalysis();
    featureSelection.clearFeatureSelection({ clearStatus: true });
    if (result?.status === 'ok' && !diagramGenerationWorkerReady.value) {
      diagramGenerationWorkerReady.value = true;
      diagramGenerationWorkerStatus.value = 'Diagram engine ready.';
      diagramGenerationWorkerError.value = null;
    }
    return result;
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
      blastSource.value === 'losat' &&
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
    const group = (Array.isArray(orthogroups.value) ? orthogroups.value : [])
      .find((entry) => String(entry?.id || '').trim() === orthogroupId);
    if (!group) return null;
    const members = orthogroupActions.getEnrichedOrthogroupMembers(group);
    const currentSvgId = String(cf?.svg_id || '').trim();
    const currentRecordIndex = Number(cf?.orthogroupMember?.recordIndex);
    const currentMember = members.find((member) => (
      String(member?.featureSvgId || '').trim() === currentSvgId &&
      (!Number.isInteger(currentRecordIndex) || Number(member?.recordIndex) === currentRecordIndex)
    )) || cf.orthogroupMember || null;
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

  const openClickedOrthogroupInEditor = () => {
    const orthogroupId = String(clickedFeature.value?.orthogroupId || '').trim();
    if (!orthogroupId) return;
    orthogroupActions.openOrthogroupInDrawer(orthogroupId);
    clickedFeature.value = null;
  };

  const { resetAllPositions, resetCanvasPadding } = legendLayout;

  const resetSettings = () => {
    const proceed = window.confirm(
      'Reset all settings to the webapp defaults?\n\nUploaded files and current results will be kept.'
    );
    if (!proceed) return;

    cancelDefinitionUpdate();
    resetSettingsState(state);
    matchSequenceRegistry?.reset?.();
    circularTrackNewRenderer.value = 'dinucleotide_skew';
    linearTrackNewRenderer.value = 'dinucleotide_skew';
    depthTrackUiCounts.circular = 1;
    ensureDepthTrackConfigCount(activeDepthTrackCount());
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

  const selectedPairwiseBlockOrthogroupId = ref('');
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
      title: 'Selected orthogroup',
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
    if (!row?.feature?.svg_id) return;
    openFeatureEditorForFeature(row.feature, event);
  };

  watch(clickedPairwiseMatch, (match) => {
    selectedPairwiseBlockOrthogroupId.value = '';
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return;
    const matchId = String(match?.matchId || '').trim();
    svg.querySelectorAll(PAIRWISE_MATCH_SELECTOR).forEach((element) => {
      const elementMatchId = String(
        element.getAttribute('data-gbdraw-match-id') ||
        element.getAttribute('data-gbdraw-pairwise-match-id') ||
        ''
      ).trim();
      element.classList.toggle('gbdraw-match-selected', Boolean(matchId) && elementMatchId === matchId);
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

  const getLastLine = (text) => {
    const trimmed = String(text || '').trim();
    if (!trimmed) return '';
    const lines = trimmed.split(/\r?\n/);
    return lines[lines.length - 1] || '';
  };

  const normalizeErrorSections = (sections) =>
    sections.filter((section) => section && typeof section.text === 'string' && section.text.trim() !== '');

  const buildErrorDisplay = (value) => {
    if (!value) return null;
    if (typeof value === 'object' && value.summary) {
      return {
        summary: value.summary,
        details: normalizeErrorSections(value.details || [])
      };
    }
    if (typeof value === 'string') {
      const trimmed = value.trim();
      if (trimmed.startsWith('{') && trimmed.endsWith('}')) {
        try {
          const parsed = JSON.parse(trimmed);
          if (parsed && typeof parsed === 'object' && parsed.summary) {
            return {
              summary: parsed.summary,
              details: normalizeErrorSections(parsed.details || [])
            };
          }
        } catch {
          // Fall through to plain string handling.
        }
      }
      const summary = getLastLine(trimmed) || 'Unknown error';
      const details = trimmed ? normalizeErrorSections([{ label: 'Details', text: trimmed }]) : [];
      return { summary, details };
    }
    return { summary: String(value), details: [] };
  };

  const errorDisplay = computed(() => buildErrorDisplay(errorLog.value));

  const sessionTitleLabel = computed(() => {
    const title = normalizeSessionTitle(sessionTitle.value);
    return title || 'Untitled session';
  });

  const canUseLinearRulerOnAxis = computed(
    () => form.scale_style === 'ruler' && ['above', 'below'].includes(form.linear_track_layout)
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

  const copyText = async (text) => {
    const value = String(text ?? '');
    if (navigator.clipboard?.writeText) {
      await navigator.clipboard.writeText(value);
      return;
    }

    const textarea = document.createElement('textarea');
    textarea.value = value;
    textarea.setAttribute('readonly', '');
    textarea.style.position = 'fixed';
    textarea.style.left = '-9999px';
    textarea.style.top = '0';
    document.body.appendChild(textarea);
    textarea.select();
    try {
      document.execCommand('copy');
    } finally {
      document.body.removeChild(textarea);
    }
  };

  const downloadText = (filename, text, type = 'text/plain;charset=utf-8') => {
    const value = String(text ?? '');
    if (!value) return;
    downloadTextFile(String(filename || 'gbdraw.txt'), value, type);
  };

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
    await exportSession(title);
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

  const formatLinearReorderCount = (count, label) => {
    const numeric = Number(count);
    return `${numeric} ${label}${numeric === 1 ? '' : 's'}`;
  };

  const buildLinearReorderNotice = ({ clearedBlastSlots = 0, clearedLosatNames = 0, actionLabel = 'Reordering' } = {}) => {
    const parts = [];
    if (clearedBlastSlots > 0) {
      parts.push(formatLinearReorderCount(clearedBlastSlots, 'BLAST TSV slot'));
    }
    if (clearedLosatNames > 0) {
      parts.push(formatLinearReorderCount(clearedLosatNames, 'LOSAT filename'));
    }
    if (parts.length === 0) return '';
    return `${actionLabel} cleared ${parts.join(' and ')} because adjacent pairs changed.`;
  };

  const applyLinearSeqMutation = (items, { actionLabel = 'Updating sequences' } = {}) => {
    const depthWidth = linearDepthLogicalWidth();
    const { linearSeqs: next, clearedBlastSlots, clearedLosatNames } = reconcileLinearSeqPairData(Array.from(linearSeqs), items);
    if (depthWidth > 0) {
      next.forEach((seq) => {
        seq.depth = padDepthFileSlots(seq.depth, depthWidth);
      });
    }
    linearSeqs.splice(0, linearSeqs.length, ...next);
    losatCacheInfo.value = [];
    linearReorderNotice.value = buildLinearReorderNotice({ clearedBlastSlots, clearedLosatNames, actionLabel });
  };

  const addLinearSeq = () => {
    applyLinearSeqMutation([...linearSeqs, createLinearSeq()], { actionLabel: 'Adding a sequence' });
  };

  const removeLinearSeqAt = (index) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= linearSeqs.length) return;
    const current = Array.from(linearSeqs);
    const next = current.filter((_, currentIndex) => currentIndex !== idx);
    applyLinearSeqMutation(next, { actionLabel: 'Removing a sequence' });
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
      losatCacheInfo.value = [];
      linearReorderNotice.value = '';
      return;
    }

    const otherField = field === 'gff' ? 'fasta' : 'gff';
    if (!nextValue && !seq[otherField]) {
      removeLinearSeqAt(idx);
      return;
    }

    seq[field] = nextValue;
    losatCacheInfo.value = [];
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
    applyLinearSeqMutation(current, { actionLabel: 'Reordering' });
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
    pyodideReady,
    diagramGenerationWorkerReady,
    diagramGenerationWorkerStatus,
    diagramGenerationWorkerError,
    processing,
    processingStatus,
    generationCancelRequested,
    loadingStatus,
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
    blastSource,
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
    linearComparisons,
    linearLayoutTokens,
    syncLinearRecordLayout,
    setLinearRecordRow,
    moveLinearRecordWithinRow,
    addLinearComparison,
    removeLinearComparison,
    setLinearComparisonFile,
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
    optionalNumberInputValue,
    setOptionalNumberInputValue,
    autoValueText: autoValueDisplay.autoValueText,
    autoValueVisible: autoValueDisplay.autoValueVisible,
    canUseLinearRulerOnAxis,
    circularTrackNewRenderer,
    linearTrackNewRenderer,
    linearTrackSlotsPanelOpen,
    toggleLinearTrackSlotsPanel,
    circularTrackRenderers: circularTrackSlotEditor.circularTrackRenderers,
    circularTrackRendererLabel: circularTrackSlotEditor.circularTrackRendererLabel,
    resetCircularTrackSlotsFromSimpleControls: circularTrackSlotEditor.resetCircularTrackSlotsFromSimpleControls,
    resetCircularTrackSlotsToPreset: circularTrackSlotEditor.resetCircularTrackSlotsToPreset,
    applyCircularTrackPreset: circularTrackSlotEditor.applyCircularTrackPreset,
    setCircularTrackSlotsEnabled: circularTrackSlotEditor.setCircularTrackSlotsEnabled,
    setCircularGcSuppressed: circularTrackSlotEditor.setCircularGcSuppressed,
    setCircularSkewSuppressed: circularTrackSlotEditor.setCircularSkewSuppressed,
    addCircularTrackSlot: circularTrackSlotEditor.addCircularTrackSlot,
    duplicateCircularTrackSlot: circularTrackSlotEditor.duplicateCircularTrackSlot,
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
    linearTrackRendererLabel: linearTrackSlotEditor.linearTrackRendererLabel,
    resetLinearTrackSlotsFromSimpleControls: linearTrackSlotEditor.resetLinearTrackSlotsFromSimpleControls,
    ensureLinearTrackDepthSlots: linearTrackSlotEditor.ensureLinearTrackDepthSlots,
    setLinearTrackSlotsEnabled: linearTrackSlotEditor.setLinearTrackSlotsEnabled,
    addLinearTrackSlot: linearTrackSlotEditor.addLinearTrackSlot,
    duplicateLinearTrackSlot: linearTrackSlotEditor.duplicateLinearTrackSlot,
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
    openRightDrawerTab: orthogroupActions.openRightDrawerTab,
    closeRightDrawer: orthogroupActions.closeRightDrawer,
    openOrthogroupInDrawer: orthogroupActions.openOrthogroupInDrawer,
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
    showFeaturePanel,
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
    getFeatureVisibility,
    moveFeatureVisibilityRuleDown: moveFeatureVisibilityRuleDownWithHistory,
    moveFeatureVisibilityRuleUp: moveFeatureVisibilityRuleUpWithHistory,
    removeFeatureVisibilityRule: removeFeatureVisibilityRuleWithHistory,
    setFeatureVisibility: setFeatureVisibilityWithHistory,
    setFeatureVisibilityRuleField: setFeatureVisibilityRuleFieldWithHistory,
    requestFeatureColorChange: requestFeatureColorChangeWithHistory,
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
    copyText,
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
    resetClickedFeatureStroke: resetClickedFeatureStrokeWithHistory,
    applyStrokeToAllSiblings: applyStrokeToAllSiblingsWithHistory,
    colorScopeDialog,
    featureVisibilityScopeDialog,
    handleColorScopeChoice: handleColorScopeChoiceWithHistory,
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
    showLegendPanel,
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
