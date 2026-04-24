import { state, createLinearSeq, reconcileLinearSeqPairData } from '../state.js';
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
import { createHistoryManager } from './history.js';
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
    linearReorderNotice,
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
    featurePanelTab,
    featureSearch,
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
    featureVisibilityOverrides,
    svgContainer,
    clickedFeature,
    clickedFeaturePos,
    featurePopupRef,
    featurePopupDrag,
    clickedLabel,
    clickedLabelPos,
    colorScopeDialog,
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

  const { handleWheel, startPan, doPan, endPan, resetPreviewViewport } = createPanZoom(state);
  const { startResizing } = createSidebarResize(state);

  const historyManager = createHistoryManager({ state, nextTick });
  const settleForHistoryCapture = async ({ waitForLabelReflow = false } = {}) => {
    await nextTick();
    await Promise.resolve();

    if (waitForLabelReflow) {
      const deadline = Date.now() + 5000;
      let sawActiveReflow = labelReflowProcessing.value;
      while (Date.now() < deadline) {
        if (labelReflowProcessing.value) {
          sawActiveReflow = true;
          await new Promise((resolve) => window.setTimeout(resolve, 25));
          continue;
        }
        if (sawActiveReflow) break;
        await new Promise((resolve) => window.setTimeout(resolve, 25));
        if (!labelReflowProcessing.value) break;
      }
    }

    await nextTick();
    await Promise.resolve();
  };

  const captureHistoryCheckpoint = async (label, options = {}) => {
    await settleForHistoryCapture(options);
    await historyManager.pushCheckpoint(label);
  };

  const wrapHistoryAction = (action, label, options = {}) => {
    return async (...args) => {
      const result = await action(...args);
      await captureHistoryCheckpoint(label, options);
      return result;
    };
  };

  const legendActions = createLegendManager({
    state,
    getPyodide,
    debugLog,
    onLegendDragCommitted: () => {
      captureHistoryCheckpoint('legend-drag');
    }
  });
  const svgActions = createSvgStyles({ state, watch, legendActions });
  const featureActions = createFeatureEditor({
    state,
    nextTick,
    legendActions,
    svgActions
  });

  const legendLayout = createLegendLayout({
    state,
    debugLog,
    legendActions,
    svgActions,
    onDiagramDragCommitted: () => {
      captureHistoryCheckpoint('diagram-drag');
    }
  });
  const {
    runAnalysis: runGeneratedDiagramAnalysis,
    runLabelReflow,
    refreshCircularRecordOrder,
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
    resetPreviewViewport
  });
  const resultsManager = createResultsManager({
    state,
    getPyodide,
    legendLayout,
    rerenderLinearDefinitions: runLabelReflow
  });

  setupGlobalUiEvents({
    state,
    onMounted,
    onUnmounted,
    historyActions: historyManager
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
    resetPreviewViewport
  });

  const {
    addNewLegendEntry: addNewLegendEntryRaw,
    updateLegendEntryColor: updateLegendEntryColorRaw,
    deleteLegendEntry: deleteLegendEntryRaw,
    moveLegendEntryUp: moveLegendEntryUpRaw,
    moveLegendEntryDown: moveLegendEntryDownRaw,
    sortLegendEntries: sortLegendEntriesRaw,
    sortLegendEntriesByDefault: sortLegendEntriesByDefaultRaw,
    resetLegendPosition: resetLegendPositionRaw,
    getLegendEntryStrokeColor,
    getLegendEntryStrokeWidth,
    updateLegendEntryStrokeColor: updateLegendEntryStrokeColorRaw,
    updateLegendEntryStrokeWidth: updateLegendEntryStrokeWidthRaw,
    resetLegendEntryStroke: resetLegendEntryStrokeRaw,
    resetAllStrokes: resetAllStrokesRaw
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
    getFeatureVisibility,
    setFeatureVisibility: setFeatureVisibilityRaw,
    updateClickedFeatureVisibility: updateClickedFeatureVisibilityRaw,
    requestFeatureColorChange: requestFeatureColorChangeRaw,
    updateClickedFeatureColor: updateClickedFeatureColorRaw,
    handleColorScopeChoice: handleColorScopeChoiceRaw,
    handleLegendNameCommit: handleLegendNameCommitRaw,
    handleLegendRenameChoice: handleLegendRenameChoiceRaw,
    selectLegendNameOption,
    renameLegendEntry: renameLegendEntryRaw,
    handleResetColorChoice: handleResetColorChoiceRaw,
    resetClickedFeatureFillColor: resetClickedFeatureFillColorRaw,
    updateClickedFeatureStroke: updateClickedFeatureStrokeRaw,
    resetClickedFeatureStroke: resetClickedFeatureStrokeRaw,
    applyStrokeToAllSiblings: applyStrokeToAllSiblingsRaw,
    setFeatureColor: setFeatureColorRaw,
    openFeatureEditorForFeature,
    getEditableLabelByFeatureId,
    syncLabelEditor,
    downloadLabelOverrideTable,
    loadLabelOverrideTable,
    updateClickedFeatureLabelText: updateClickedFeatureLabelTextRaw,
    handleLabelTextScopeChoice: handleLabelTextScopeChoiceRaw,
    handleGlobalLabelModeChoice: handleGlobalLabelModeChoiceRaw,
    requestLabelTextChangeByFeatureId,
    requestLabelTextChangeByKey,
    resetAllLabelTextOverrides: resetAllLabelTextOverridesRaw
  } = featureActions;

  const { updatePalette, resetColors, cancelDefinitionUpdate } = resultsManager;
  const { resetAllPositions: resetAllPositionsRaw, resetCanvasPadding: resetCanvasPaddingRaw } = legendLayout;

  const addSpecificRuleWithHistory = wrapHistoryAction(addSpecificRule, 'specific-rule-add');
  const applySpecificRulePresetWithHistory = wrapHistoryAction(
    applySpecificRulePreset,
    'specific-rule-preset'
  );
  const clearAllSpecificRulesWithHistory = wrapHistoryAction(
    clearAllSpecificRules,
    'specific-rule-clear'
  );
  const setFeatureVisibility = wrapHistoryAction(setFeatureVisibilityRaw, 'feature-visibility');
  const updateClickedFeatureVisibility = wrapHistoryAction(
    updateClickedFeatureVisibilityRaw,
    'feature-visibility'
  );
  const requestFeatureColorChange = wrapHistoryAction(requestFeatureColorChangeRaw, 'feature-color');
  const updateClickedFeatureColor = wrapHistoryAction(updateClickedFeatureColorRaw, 'feature-color');
  const handleColorScopeChoice = wrapHistoryAction(handleColorScopeChoiceRaw, 'feature-color');
  const handleLegendNameCommit = wrapHistoryAction(handleLegendNameCommitRaw, 'legend-name');
  const handleLegendRenameChoice = wrapHistoryAction(handleLegendRenameChoiceRaw, 'legend-rename');
  const renameLegendEntry = wrapHistoryAction(renameLegendEntryRaw, 'legend-rename');
  const handleResetColorChoice = wrapHistoryAction(handleResetColorChoiceRaw, 'feature-color-reset');
  const resetClickedFeatureFillColor = wrapHistoryAction(
    resetClickedFeatureFillColorRaw,
    'feature-color-reset'
  );
  const updateClickedFeatureStroke = wrapHistoryAction(updateClickedFeatureStrokeRaw, 'feature-stroke');
  const resetClickedFeatureStroke = wrapHistoryAction(resetClickedFeatureStrokeRaw, 'feature-stroke-reset');
  const applyStrokeToAllSiblings = wrapHistoryAction(applyStrokeToAllSiblingsRaw, 'feature-stroke');
  const setFeatureColor = wrapHistoryAction(setFeatureColorRaw, 'feature-color');
  const updateClickedFeatureLabelText = wrapHistoryAction(
    updateClickedFeatureLabelTextRaw,
    'label-edit',
    { waitForLabelReflow: true }
  );
  const handleLabelTextScopeChoice = wrapHistoryAction(
    handleLabelTextScopeChoiceRaw,
    'label-edit',
    { waitForLabelReflow: true }
  );
  const handleGlobalLabelModeChoice = wrapHistoryAction(
    handleGlobalLabelModeChoiceRaw,
    'label-visibility',
    { waitForLabelReflow: true }
  );
  const resetAllLabelTextOverrides = wrapHistoryAction(
    resetAllLabelTextOverridesRaw,
    'label-reset',
    { waitForLabelReflow: true }
  );

  const addNewLegendEntry = wrapHistoryAction(addNewLegendEntryRaw, 'legend-add');
  const updateLegendEntryColor = wrapHistoryAction(updateLegendEntryColorRaw, 'legend-color');
  const deleteLegendEntry = wrapHistoryAction(deleteLegendEntryRaw, 'legend-delete');
  const moveLegendEntryUp = wrapHistoryAction(moveLegendEntryUpRaw, 'legend-reorder');
  const moveLegendEntryDown = wrapHistoryAction(moveLegendEntryDownRaw, 'legend-reorder');
  const sortLegendEntries = wrapHistoryAction(sortLegendEntriesRaw, 'legend-sort');
  const sortLegendEntriesByDefault = wrapHistoryAction(sortLegendEntriesByDefaultRaw, 'legend-sort');
  const resetLegendPosition = wrapHistoryAction(resetLegendPositionRaw, 'legend-reset');
  const updateLegendEntryStrokeColor = wrapHistoryAction(
    updateLegendEntryStrokeColorRaw,
    'legend-stroke'
  );
  const updateLegendEntryStrokeWidth = wrapHistoryAction(
    updateLegendEntryStrokeWidthRaw,
    'legend-stroke'
  );
  const resetLegendEntryStroke = wrapHistoryAction(resetLegendEntryStrokeRaw, 'legend-stroke-reset');
  const resetAllStrokes = wrapHistoryAction(resetAllStrokesRaw, 'legend-stroke-reset');
  const resetAllPositions = wrapHistoryAction(resetAllPositionsRaw, 'position-reset');
  const resetCanvasPadding = wrapHistoryAction(resetCanvasPaddingRaw, 'canvas-padding');

  const setCanvasPadding = async (side, value) => {
    if (!['top', 'right', 'bottom', 'left'].includes(side)) return;
    const normalized = value === '' || value === null || value === undefined ? 0 : Number(value);
    const nextValue = Number.isFinite(normalized) ? normalized : 0;
    if (canvasPadding[side] === nextValue) return;
    canvasPadding[side] = nextValue;
    await captureHistoryCheckpoint('canvas-padding');
  };

  const undo = async () => historyManager.undo();
  const redo = async () => historyManager.redo();

  const runAnalysis = async () => {
    cancelDefinitionUpdate();
    historyManager.clearHistory();
    const result = await runGeneratedDiagramAnalysis();
    if (result?.status === 'ok') {
      await historyManager.resetHistory('generate');
    } else {
      historyManager.clearHistory();
    }
    return result;
  };

  const importSessionWithHistory = async (e) => {
    historyManager.clearHistory();
    await importSession(e);
    if (results.value.length > 0) {
      await historyManager.resetHistory('session-import');
    }
  };

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
    const { linearSeqs: next, clearedBlastSlots, clearedLosatNames } = reconcileLinearSeqPairData(Array.from(linearSeqs), items);
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
    linearReorderNotice,
    addLinearSeq,
    removeLastLinearSeq,
    setLinearSeqPrimaryFile,
    canMoveLinearSeqUp,
    canMoveLinearSeqDown,
    moveLinearSeqUp,
    moveLinearSeqDown,
    form,
    adv,
    canUseLinearRulerOnAxis,
    losat,
    losatCacheInfo,
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
    addSpecificRule: addSpecificRuleWithHistory,
    applySpecificRulePreset: applySpecificRulePresetWithHistory,
    clearAllSpecificRules: clearAllSpecificRulesWithHistory,
    extractedFeatures,
    featureRecordIds,
    selectedFeatureRecordIdx,
    showFeaturePanel,
    featurePanelTab,
    featureSearch,
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
    featureVisibilityOverrides,
    getFeatureColor,
    getFeatureVisibility,
    setFeatureVisibility,
    canUndo: historyManager.canUndo,
    canRedo: historyManager.canRedo,
    undo,
    redo,
    requestFeatureColorChange,
    setFeatureColor,
    canEditFeatureColor,
    getEditableLabelByFeatureId,
    svgContainer,
    clickedFeature,
    clickedFeaturePos,
    clickedLabel,
    clickedLabelPos,
    featurePopupRef,
    startFeaturePopupDrag,
    clickedFeatureLocation,
    specificRuleLegendOptions,
    updateClickedFeatureColor,
    updateClickedFeatureVisibility,
    handleLegendNameCommit,
    selectLegendNameOption,
    resetClickedFeatureFillColor,
    updateClickedFeatureStroke,
    resetClickedFeatureStroke,
    applyStrokeToAllSiblings,
    colorScopeDialog,
    handleColorScopeChoice,
    legendRenameDialog,
    handleLegendRenameChoice,
    resetColorDialog,
    handleResetColorChoice,
    labelTextScopeDialog,
    globalLabelModeDialog,
    updateClickedFeatureLabelText,
    handleLabelTextScopeChoice,
    handleGlobalLabelModeChoice,
    requestLabelTextChangeByFeatureId,
    requestLabelTextChangeByKey,
    resetAllLabelTextOverrides,
    downloadLabelOverrideTable,
    loadLabelOverrideTable,
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
    canvasPadding,
    setCanvasPadding,
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
    importSession: importSessionWithHistory,
    manualPriorityRules,
    newPriorityRule,
    addPriorityRule
  };
};
