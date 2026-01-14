import { state } from '../state.js';
import { debugLog } from '../config.js';
import { downloadSVG, downloadPNG, downloadPDF } from '../services/export.js';
import { exportConfig, importConfig } from '../services/config.js';
import { createPanZoom, createSidebarResize, setupGlobalUiEvents } from './ui.js';
import { createFeatureEditor } from './feature-editor.js';
import { createSvgStyles } from './svg-styles.js';
import { createLegendManager } from './legend.js';
import { createPyodideManager } from './pyodide.js';
import { createRunAnalysis } from './run-analysis.js';
import { createLegendLayout } from './legend-layout.js';
import { createResultsManager } from './results.js';
import { setupWatchers } from './watchers.js';

const { onMounted, onUnmounted, watch, nextTick } = window.Vue;

export const createAppSetup = () => {
  const {
    pyodideReady,
    processing,
    loadingStatus,
    errorLog,
    results,
    selectedResultIndex,
    pairwiseMatchFactors,
    svgContent,
    zoom,
    isPanning,
    panStart,
    canvasContainerRef,
    mode,
    circularLegendPosition,
    linearLegendPosition,
    cInputType,
    lInputType,
    files,
    linearSeqs,
    form,
    adv,
    paletteNames,
    selectedPalette,
    currentColors,
    filterMode,
    manualBlacklist,
    manualWhitelist,
    manualSpecificRules,
    newSpecRule,
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
  const { runAnalysis } = createRunAnalysis({
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
    addSpecificRule,
    clearAllSpecificRules,
    getFeatureColor,
    canEditFeatureColor,
    updateClickedFeatureColor,
    handleColorScopeChoice,
    handleResetColorChoice,
    resetClickedFeatureFillColor,
    updateClickedFeatureStroke,
    resetClickedFeatureStroke,
    applyStrokeToAllSiblings,
    setFeatureColor
  } = featureActions;

  const { updatePalette, resetColors } = resultsManager;

  const { resetAllPositions, resetCanvasPadding } = legendLayout;

  return {
    pyodideReady,
    processing,
    loadingStatus,
    errorLog,
    results,
    selectedResultIndex,
    svgContent,
    zoom,
    handleWheel,
    canvasContainerRef,
    startPan,
    doPan,
    endPan,
    sidebarWidth,
    startResizing,
    mode,
    cInputType,
    lInputType,
    files,
    linearSeqs,
    form,
    adv,
    paletteNames,
    selectedPalette,
    currentColors,
    updatePalette,
    resetColors,
    filterMode,
    manualBlacklist,
    manualWhitelist,
    featureKeys,
    newColorFeat,
    newColorVal,
    addCustomColor,
    newFeatureToAdd,
    addFeature,
    manualSpecificRules,
    newSpecRule,
    addSpecificRule,
    clearAllSpecificRules,
    extractedFeatures,
    featureRecordIds,
    selectedFeatureRecordIdx,
    showFeaturePanel,
    featureSearch,
    filteredFeatures,
    featureColorOverrides,
    getFeatureColor,
    setFeatureColor,
    canEditFeatureColor,
    svgContainer,
    clickedFeature,
    clickedFeaturePos,
    updateClickedFeatureColor,
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
    importConfig,
    manualPriorityRules,
    newPriorityRule,
    addPriorityRule
  };
};
