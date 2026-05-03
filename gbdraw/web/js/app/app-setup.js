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
import { createResultsManager } from './results.js';
import { setupWatchers } from './watchers.js';
import { createOrthogroupEditor } from './orthogroups.js';

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
    featureExtractionPending,
    featureExtractionError,
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
    getFeatureColor,
    canEditFeatureColor,
    getFeatureVisibility,
    setFeatureVisibility,
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

  const runAnalysis = async () => {
    cancelDefinitionUpdate();
    return runGeneratedDiagramAnalysis();
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
    const members = Array.isArray(group.members) ? group.members : [];
    const currentSvgId = String(cf?.svg_id || '').trim();
    const currentRecordIndex = Number(cf?.orthogroupMember?.recordIndex);
    const currentMember = members.find((member) => (
      String(member?.featureSvgId || '').trim() === currentSvgId &&
      (!Number.isInteger(currentRecordIndex) || Number(member?.recordIndex) === currentRecordIndex)
    )) || cf.orthogroupMember || null;
    const grouped = new Map();
    members.forEach((member) => {
      const recordIndex = Number(member?.recordIndex);
      const key = Number.isInteger(recordIndex) ? recordIndex : -1;
      if (!grouped.has(key)) grouped.set(key, []);
      grouped.get(key).push(member);
    });
    const membersByRecord = Array.from(grouped.entries())
      .sort((left, right) => left[0] - right[0])
      .map(([recordIndex, recordMembers]) => ({
        recordIndex,
        recordLabel: recordIndex >= 0
          ? (linearSeqs[recordIndex]?.name || linearSeqs[recordIndex]?.gb?.name || linearSeqs[recordIndex]?.gff?.name || `Record ${recordIndex + 1}`)
          : 'Record',
        members: recordMembers
      }));
    return {
      id: orthogroupId,
      displayName: orthogroupActions.resolveOrthogroupName(group),
      description: orthogroupActions.resolveOrthogroupDescription(group),
      candidates: Array.isArray(group.nameCandidates) ? group.nameCandidates : [],
      memberCount: Number(group.member_count || members.length || 0),
      recordCoverage: Number(group.record_coverage_count || membersByRecord.length || 0),
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
    isOrthogroupRenamed: orthogroupActions.isOrthogroupRenamed,
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
    featureExtractionPending,
    featureExtractionError,
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
    canUseClickedOrthogroupActions,
    clickedOrthogroupDetail,
    alignByClickedOrthogroup,
    highlightClickedOrthogroup,
    clearOrthogroupHighlight,
    resetOrthogroupAlignment,
    openClickedOrthogroupInEditor,
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
