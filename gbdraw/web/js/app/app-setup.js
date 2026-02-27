import { state } from '../state.js';
import { debugLog } from '../config.js';
import { downloadSVG, downloadPNG, downloadPDF } from '../services/export.js';
import { exportConfig, exportSession, importConfig, importSession } from '../services/config.js';
import { createPanZoom, createSidebarResize, setupGlobalUiEvents } from './ui.js';
import { createFeatureEditor } from './feature-editor.js';
import { createSvgStyles } from './svg-styles.js';
import { createLegendManager } from './legend.js';
import { createPyodideManager } from './pyodide.js';
import { createRunAnalysis } from './run-analysis.js';
import { createLegendLayout } from './legend-layout.js';
import { createResultsManager } from './results.js';
import { setupWatchers } from './watchers.js';

const { onMounted, onUnmounted, watch, nextTick, computed } = window.Vue;

export const createAppSetup = () => {
  const {
    pyodideReady,
    processing,
    loadingStatus,
    errorLog,
    sessionTitle,
    results,
    selectedResultIndex,
    pairwiseMatchFactors,
    svgContent,
    zoom,
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
    linearSeqs,
    form,
    adv,
    losat,
    losatCacheInfo,
    paletteNames,
    selectedPalette,
    currentColors,
    filterMode,
    manualBlacklist,
    manualWhitelist,
    manualSpecificRules,
    newSpecRule,
    specificRulePresets,
    selectedSpecificPreset,
    specificRulePresetLoading,
    downloadDpi,
    extractedFeatures,
    featureRecordIds,
    selectedFeatureRecordIdx,
    showFeaturePanel,
    featureSearch,
    featureColorOverrides,
    svgContainer,
    clickedFeature,
    clickedFeaturePos,
    featurePopupRef,
    featurePopupDrag,
    colorScopeDialog,
    resetColorDialog,
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
    newColorFeat,
    newColorVal,
    manualPriorityRules,
    newPriorityRule,
    newFeatureToAdd,
    addedLegendCaptions,
    fileLegendCaptions,
    filteredFeatures
  } = state;

  const pyodideManager = createPyodideManager({ state });
  const getPyodide = pyodideManager.getPyodide;

  const { handleWheel, startPan, doPan, endPan } = createPanZoom(state);
  const { startResizing } = createSidebarResize(state);

  const legendActions = createLegendManager({ state, getPyodide, debugLog });
  const svgActions = createSvgStyles({ state, watch, legendActions });
  const featureActions = createFeatureEditor({
    state,
    nextTick,
    legendActions,
    svgActions
  });

  setupGlobalUiEvents({ state, onMounted, onUnmounted });

  const legendLayout = createLegendLayout({ state, debugLog, legendActions, svgActions });
  const resultsManager = createResultsManager({ state, getPyodide });
  const {
    runAnalysis,
    downloadLosatCache,
    downloadLosatPair,
    setLosatPairFilename,
    clearLosatCache,
    getLosatPairDefaultName
  } = createRunAnalysis({
    state,
    getPyodide,
    writeFileToFs: pyodideManager.writeFileToFs,
    refreshFeatureOverrides: featureActions.refreshFeatureOverrides
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
    resultsManager
  });

  const {
    addNewLegendEntry,
    updateLegendEntryColor,
    updateLegendEntryCaption,
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
    getFeatureColor,
    canEditFeatureColor,
    requestFeatureColorChange,
    updateClickedFeatureColor,
    handleColorScopeChoice,
    handleLegendNameCommit,
    selectLegendNameOption,
    handleResetColorChoice,
    resetClickedFeatureFillColor,
    updateClickedFeatureStroke,
    resetClickedFeatureStroke,
    applyStrokeToAllSiblings,
    setFeatureColor
  } = featureActions;

  const { updatePalette, resetColors } = resultsManager;

  const { resetAllPositions, resetCanvasPadding } = legendLayout;

  const isInteractiveTarget = (target) => {
    if (!target) return false;
    return Boolean(target.closest('input, textarea, select, button, label, a, [data-nodrag="true"]'));
  };

  const onFeaturePopupDrag = (event) => {
    if (!featurePopupDrag.active) return;
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

  const endFeaturePopupDrag = () => {
    if (!featurePopupDrag.active) return;
    featurePopupDrag.active = false;
    document.removeEventListener('mousemove', onFeaturePopupDrag);
    document.removeEventListener('mouseup', endFeaturePopupDrag);
  };

  const startFeaturePopupDrag = (event) => {
    if (event.button !== 0) return;
    if (!clickedFeature.value) return;
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

  return {
    pyodideReady,
    processing,
    loadingStatus,
    errorLog,
    errorDisplay,
    sessionTitle,
    sessionTitleLabel,
    results,
    selectedResultIndex,
    svgContent,
    zoom,
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
    linearSeqs,
    form,
    adv,
    losat,
    losatCacheInfo,
    paletteNames,
    selectedPalette,
    currentColors,
    updatePalette,
    resetColors,
    downloadLosatCache,
    downloadLosatPair,
    setLosatPairFilename,
    clearLosatCache,
    getLosatPairDefaultName,
    filterMode,
    manualBlacklist,
    manualWhitelist,
    featureKeys,
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
    selectedSpecificPreset,
    specificRulePresetLoading,
    addSpecificRule,
    applySpecificRulePreset,
    clearAllSpecificRules,
    extractedFeatures,
    featureRecordIds,
    selectedFeatureRecordIdx,
    showFeaturePanel,
    featureSearch,
    filteredFeatures,
    featureColorOverrides,
    getFeatureColor,
    requestFeatureColorChange,
    setFeatureColor,
    canEditFeatureColor,
    svgContainer,
    clickedFeature,
    clickedFeaturePos,
    featurePopupRef,
    startFeaturePopupDrag,
    clickedFeatureLocation,
    specificRuleLegendOptions,
    updateClickedFeatureColor,
    handleLegendNameCommit,
    selectLegendNameOption,
    resetClickedFeatureFillColor,
    updateClickedFeatureStroke,
    resetClickedFeatureStroke,
    applyStrokeToAllSiblings,
    colorScopeDialog,
    handleColorScopeChoice,
    resetColorDialog,
    handleResetColorChoice,
    showLegendPanel,
    legendEntries,
    newLegendCaption,
    newLegendColor,
    updateLegendEntryColor,
    updateLegendEntryCaption,
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
    canvasPadding,
    showCanvasControls,
    resetCanvasPadding,
    downloadDpi,
    runAnalysis,
    downloadSVG,
    downloadPNG,
    downloadPDF,
    exportConfig,
    saveSessionWithTitle,
    editSessionTitle,
    importConfig,
    importSession,
    manualPriorityRules,
    newPriorityRule,
    addPriorityRule
  };
};
