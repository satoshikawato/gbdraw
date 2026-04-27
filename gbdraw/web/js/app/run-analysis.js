import { runLosatPairsParallel } from '../services/losat.js';
import { buildLabelOverrideTsv } from './feature-editor/label-override-table.js';

const downloadTextFile = (filename, text) => {
  const safeName = filename || 'losat.tsv';
  const blob = new Blob([text], { type: 'text/tab-separated-values' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = safeName;
  link.click();
  URL.revokeObjectURL(url);
};

const hashText = async (text) => {
  if (globalThis.crypto?.subtle) {
    const buffer = await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text));
    return Array.from(new Uint8Array(buffer))
      .map((b) => b.toString(16).padStart(2, '0'))
      .join('');
  }
  let hash = 2166136261;
  for (let i = 0; i < text.length; i++) {
    hash ^= text.charCodeAt(i);
    hash = Math.imul(hash, 16777619);
  }
  return `fnv1a-${(hash >>> 0).toString(16)}`;
};

const getNow = () => (globalThis.performance?.now ? performance.now() : Date.now());
const formatDuration = (ms) => `${(ms / 1000).toFixed(2)}s`;

const makeSafeFilename = (name) => {
  const cleaned = String(name || '').replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '');
  return cleaned || 'losat';
};
const escapeRegexLiteral = (value) => String(value ?? '').replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
const normalizeFeatureVisibilityMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return normalized === 'on' || normalized === 'off' ? normalized : 'default';
};
const normalizeMultiRecordSizeMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  if (normalized === 'sqrt') return 'auto';
  return ['auto', 'linear', 'equal'].includes(normalized) ? normalized : 'auto';
};
const normalizeMultiRecordMinRadiusRatio = (value) => {
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric > 0 && numeric <= 1 ? numeric : 0.55;
};
const normalizeMultiRecordColumnGapRatio = (value) => {
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric >= 0 ? numeric : 0.10;
};
const normalizeMultiRecordRowGapRatio = (value) => {
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric >= 0 ? numeric : 0.05;
};
const DEFAULT_LINEAR_BLAST_FILTERS = Object.freeze({
  bitscore: 50,
  evalue: '1e-2',
  identity: 0,
  alignment_length: 0
});
const normalizeBlastThresholdNumber = (value, defaultValue, { integer = false } = {}) => {
  if (value === null || value === undefined || value === '') return defaultValue;
  const numeric = Number(value);
  if (!Number.isFinite(numeric) || numeric < 0) return defaultValue;
  if (integer && !Number.isInteger(numeric)) return defaultValue;
  return numeric;
};
const normalizeBlastThresholdText = (value, defaultValue) => {
  const normalized = String(value ?? '').trim();
  return normalized === '' ? defaultValue : normalized;
};
const normalizeMultiRecordPositions = (value, { maxRow = Number.POSITIVE_INFINITY } = {}) => {
  if (!Array.isArray(value)) return [];
  const deduped = [];
  const seen = new Set();
  value.forEach((item) => {
    let selector = '';
    let row = 1;
    if (item && typeof item === 'object' && !Array.isArray(item)) {
      selector = String(item.selector ?? '').trim();
      row = Number(item.row);
    } else if (typeof item === 'string') {
      const raw = String(item || '').trim();
      if (!raw || !raw.includes('@')) return;
      const parts = raw.split('@');
      if (parts.length < 2) return;
      selector = parts.slice(0, -1).join('@').trim();
      row = Number(parts[parts.length - 1]);
    }
    if (!selector || seen.has(selector)) return;
    const normalizedMaxRow = Number.isInteger(maxRow) && maxRow > 0 ? maxRow : Number.POSITIVE_INFINITY;
    const normalizedRowRaw = Number.isInteger(row) && row > 0 ? row : 1;
    const normalizedRow = Number.isFinite(normalizedMaxRow)
      ? Math.min(normalizedRowRaw, normalizedMaxRow)
      : normalizedRowRaw;
    seen.add(selector);
    deduped.push({ selector, row: normalizedRow });
  });
  return deduped;
};
const sortMultiRecordPositionsByRow = (positions) => {
  if (!Array.isArray(positions)) return [];
  return positions
    .map((entry, index) => ({ ...entry, __index: index }))
    .sort((left, right) => {
      const leftRow = Number(left.row);
      const rightRow = Number(right.row);
      if (leftRow !== rightRow) return leftRow - rightRow;
      return left.__index - right.__index;
    })
    .map(({ __index, ...entry }) => entry);
};
const buildDefaultMultiRecordPositions = (selectors) => {
  const normalizedSelectors = Array.isArray(selectors)
    ? selectors.map((value) => String(value ?? '').trim()).filter(Boolean)
    : [];
  if (normalizedSelectors.length === 0) return [];
  const cols = Math.ceil(Math.sqrt(normalizedSelectors.length));
  return normalizedSelectors.map((selector, index) => ({
    selector,
    row: Math.floor(index / cols) + 1
  }));
};
const mergeCircularRecordPositions = (records, currentPositions) => {
  const availableSelectors = Array.isArray(records)
    ? records.map((entry) => String(entry?.selector || '').trim()).filter(Boolean)
    : [];
  if (availableSelectors.length === 0) return [];
  const availableSet = new Set(availableSelectors);
  const defaultPositions = buildDefaultMultiRecordPositions(availableSelectors);
  const defaultRowBySelector = new Map(defaultPositions.map((entry) => [entry.selector, entry.row]));
  const normalizedCurrent = normalizeMultiRecordPositions(currentPositions, { maxRow: availableSelectors.length });
  const nextPositions = [];
  const seen = new Set();

  normalizedCurrent.forEach((entry) => {
    if (!availableSet.has(entry.selector) || seen.has(entry.selector)) return;
    seen.add(entry.selector);
    nextPositions.push({
      selector: entry.selector,
      row: Number.isInteger(entry.row) && entry.row > 0 ? entry.row : (defaultRowBySelector.get(entry.selector) || 1)
    });
  });
  availableSelectors.forEach((selector) => {
    if (seen.has(selector)) return;
    seen.add(selector);
    nextPositions.push({
      selector,
      row: defaultRowBySelector.get(selector) || 1
    });
  });
  return sortMultiRecordPositionsByRow(
    normalizeMultiRecordPositions(nextPositions, { maxRow: availableSelectors.length })
  );
};
const buildMultiRecordPositionToken = (entry) => {
  if (!entry || typeof entry !== 'object') return '';
  const selector = String(entry.selector || '').trim();
  const row = Number(entry.row);
  if (!selector || !Number.isInteger(row) || row <= 0) return '';
  return `${selector}@${row}`;
};
const normalizeCircularPlotTitlePosition = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['none', 'top', 'bottom'].includes(normalized) ? normalized : 'none';
};
const normalizeLinearPlotTitlePosition = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['center', 'top', 'bottom'].includes(normalized) ? normalized : 'bottom';
};

export const createRunAnalysis = ({
  state,
  getPyodide,
  writeFileToFs,
  refreshFeatureOverrides,
  resetPreviewViewport
}) => {
  const {
    pyodideReady,
    processing,
    results,
    selectedResultIndex,
    errorLog,
    zoom,
    skipCaptureBaseConfig,
    skipPositionReapply,
    pairwiseMatchFactors,
    addedLegendCaptions,
    fileLegendCaptions,
    featureColorOverrides,
    featureVisibilityOverrides,
    legendEntries,
    deletedLegendEntries,
    legendColorOverrides,
    originalLegendOrder,
    originalLegendColors,
    selectedPalette,
    currentColors,
    appliedPaletteName,
    appliedPaletteColors,
    pendingPaletteName,
    pendingPaletteColors,
    filterMode,
    manualSpecificRules,
    manualWhitelist,
    manualBlacklist,
    manualPriorityRules,
    form,
    adv,
    mode,
    cInputType,
    lInputType,
    blastSource,
    losatProgram,
    losat,
    losatCacheInfo,
    losatCache,
    circularRecordList,
    files,
    linearSeqs,
    generatedLegendPosition,
    generatedMode,
    generatedMultiRecordCanvas,
    generatedCircularPlotTitlePosition,
    shouldDeferCircularPreviewUpdates,
    extractedFeatures,
    featureRecordIds,
    selectedFeatureRecordIdx,
    editableLabels,
    clickedLabel,
    labelTextScopeDialog,
    labelTextFeatureOverrides,
    labelTextBulkOverrides,
    labelTextFeatureOverrideSources,
    labelVisibilityOverrides,
    labelOverrideBuildWarning,
    labelReflowProcessing,
    labelReflowLastError
  } = state;
  let linearLabelSupportCache = null;
  let featureShapeSupportCache = null;
  let circularMultiRecordCanvasSupportCache = null;
  let pendingReflowRequestId = 0;
  let activeReflowRequestId = 0;
  let pendingReflowReason = 'label-edit';

  const getLastLine = (text) => {
    const trimmed = String(text || '').trim();
    if (!trimmed) return '';
    const lines = trimmed.split(/\r?\n/);
    return lines[lines.length - 1] || '';
  };

  const normalizeSections = (sections) =>
    sections.filter((section) => section && typeof section.text === 'string' && section.text.trim() !== '');

  const formatPythonError = (err) => {
    if (err && typeof err === 'object') {
      const details = normalizeSections([
        { label: 'STDERR', text: err.stderr || '' },
        { label: 'STDOUT', text: err.stdout || '' },
        { label: 'Traceback', text: err.traceback || '' }
      ]);
      let summary =
        (err.type === 'SystemExit' && err.stderr ? getLastLine(err.stderr) : '') ||
        err.message ||
        getLastLine(err.stderr || err.stdout || err.traceback) ||
        'Unknown error';
      if (err.type && summary && !summary.startsWith(err.type)) {
        summary = `${err.type}: ${summary}`;
      }
      return { summary, details };
    }

    const text = String(err || '').trim();
    const summary = getLastLine(text) || 'Unknown error';
    const details = normalizeSections([{ label: 'Details', text }]);
    return { summary, details };
  };

  const formatJsError = (err) => {
    const message = err?.message ? String(err.message) : String(err || 'Unknown error');
    return { summary: message, details: [] };
  };

  const getSeqLabel = (seq, fallback) => {
    const definition = String(seq?.definition || '').trim();
    if (definition) return definition;
    if (fallback) return fallback;
    const file = seq?.gb || seq?.fasta || seq?.gff;
    if (file?.name) {
      return String(file.name).replace(/\.[^.]+$/, '');
    }
    return '';
  };

  const normalizeLabel = (label, fallback) => {
    const base = String(label || '').trim() || String(fallback || '');
    const dotted = base.replace(/[\\s/]+/g, '.').replace(/\.+/g, '.').replace(/^\.|\.$/g, '');
    const safe = makeSafeFilename(dotted);
    return safe || makeSafeFilename(String(fallback || 'losat'));
  };

  const buildLosatSuffix = () => (losatProgram.value === 'blastn' ? 'losatn' : 'tlosatx');

  const buildLosatFilename = (leftLabel, rightLabel) => {
    const left = normalizeLabel(leftLabel, 'seq_1');
    const right = normalizeLabel(rightLabel, 'seq_2');
    return `${left}.${right}.${buildLosatSuffix()}.tsv`;
  };

  const getLosatPairDefaultName = (pairIndex, queryEntry = null, subjectEntry = null) => {
    const leftLabel = getSeqLabel(linearSeqs[pairIndex], queryEntry?.recordId || `seq_${pairIndex + 1}`);
    const rightLabel = getSeqLabel(linearSeqs[pairIndex + 1], subjectEntry?.recordId || `seq_${pairIndex + 2}`);
    return buildLosatFilename(leftLabel, rightLabel);
  };

  const normalizeLosatFilename = (name, fallback) => {
    const raw = String(name || '').trim() || String(fallback || '');
    const withExt = raw.toLowerCase().endsWith('.tsv') ? raw : `${raw}.tsv`;
    return makeSafeFilename(withExt);
  };

  const getLosatParallelWorkers = () => {
    const raw = String(losat.parallelWorkers || 'auto').trim().toLowerCase();
    if (raw === 'auto') return undefined;
    const parsed = Number(raw);
    return Number.isInteger(parsed) && parsed >= 1 && parsed <= 4 ? parsed : undefined;
  };

  const getLinearLabelOptionSupport = () => {
    if (linearLabelSupportCache) return linearLabelSupportCache;
    const pyodide = getPyodide();
    if (!pyodide) {
      linearLabelSupportCache = {
        placement: false,
        rotation: false,
        linear_label_spacing: false,
        track_layout: false,
        track_axis_gap: false,
        ruler_on_axis: false,
        scale_font_size: true,
        ruler_label_font_size: false,
        ruler_label_color: false,
        plot_title: false,
        plot_title_position: false,
        plot_title_font_size: false,
        show_replicon: false,
        hide_accession: false,
        hide_length: false
      };
      return linearLabelSupportCache;
    }
    try {
      const raw = pyodide.runPython(`
import inspect, json
import gbdraw.linear as _gbdraw_linear
_source = inspect.getsource(_gbdraw_linear._get_args)
json.dumps({
  "placement": "--label_placement" in _source,
  "rotation": "--label_rotation" in _source,
  "linear_label_spacing": "--linear_label_spacing" in _source,
  "track_layout": "--track_layout" in _source,
  "track_axis_gap": "--track_axis_gap" in _source,
  "ruler_on_axis": "--ruler_on_axis" in _source,
  "scale_font_size": "--scale_font_size" in _source,
  "ruler_label_font_size": "--ruler_label_font_size" in _source,
  "ruler_label_color": "--ruler_label_color" in _source,
  "plot_title": "--plot_title" in _source,
  "plot_title_position": "--plot_title_position" in _source,
  "plot_title_font_size": "--plot_title_font_size" in _source,
  "show_replicon": "--show_replicon" in _source,
  "hide_accession": "--hide_accession" in _source,
  "hide_length": "--hide_length" in _source,
})
      `);
      linearLabelSupportCache = JSON.parse(String(raw));
    } catch (_err) {
      linearLabelSupportCache = {
        placement: false,
        rotation: false,
        linear_label_spacing: false,
        track_layout: false,
        track_axis_gap: false,
        ruler_on_axis: false,
        scale_font_size: true,
        ruler_label_font_size: false,
        ruler_label_color: false,
        plot_title: false,
        plot_title_position: false,
        plot_title_font_size: false,
        show_replicon: false,
        hide_accession: false,
        hide_length: false
      };
    }
    return linearLabelSupportCache;
  };

  const normalizeFeatureShape = (value) => (String(value || '').trim().toLowerCase() === 'arrow' ? 'arrow' : 'rectangle');

  const getFeatureShapeOptionSupport = () => {
    if (featureShapeSupportCache) return featureShapeSupportCache;
    const pyodide = getPyodide();
    if (!pyodide) {
      featureShapeSupportCache = { circular: false, linear: false };
      return featureShapeSupportCache;
    }
    try {
      const raw = pyodide.runPython(`
import inspect, json
import gbdraw.circular as _gbdraw_circular
import gbdraw.linear as _gbdraw_linear
_circular_source = inspect.getsource(_gbdraw_circular._get_args)
_linear_source = inspect.getsource(_gbdraw_linear._get_args)
json.dumps({
  "circular": "--feature_shape" in _circular_source,
  "linear": "--feature_shape" in _linear_source,
})
      `);
      featureShapeSupportCache = JSON.parse(String(raw));
    } catch (_err) {
      featureShapeSupportCache = { circular: false, linear: false };
    }
    return featureShapeSupportCache;
  };

  const getCircularMultiRecordCanvasOptionSupport = () => {
    if (circularMultiRecordCanvasSupportCache) return circularMultiRecordCanvasSupportCache;
    const pyodide = getPyodide();
    if (!pyodide) {
      circularMultiRecordCanvasSupportCache = {
        circular: false,
        multi_record_size_mode: false,
        multi_record_min_radius_ratio: false,
        multi_record_column_gap_ratio: false,
        multi_record_row_gap_ratio: false,
        multi_record_position: false,
        plot_title: false,
        plot_title_position: false,
        plot_title_font_size: false,
        keep_full_definition_with_plot_title: false,
        tick_label_font_size: false,
        circular_label_spacing: false
      };
      return circularMultiRecordCanvasSupportCache;
    }
    try {
      const raw = pyodide.runPython(`
import inspect, json
import gbdraw.circular as _gbdraw_circular
_source = inspect.getsource(_gbdraw_circular._get_args)
json.dumps({
  "circular": "--multi_record_canvas" in _source,
  "multi_record_size_mode": "--multi_record_size_mode" in _source,
  "multi_record_min_radius_ratio": "--multi_record_min_radius_ratio" in _source,
  "multi_record_column_gap_ratio": "--multi_record_column_gap_ratio" in _source,
  "multi_record_row_gap_ratio": "--multi_record_row_gap_ratio" in _source,
  "multi_record_position": "--multi_record_position" in _source,
  "plot_title": "--plot_title" in _source,
  "plot_title_position": "--plot_title_position" in _source,
  "plot_title_font_size": "--plot_title_font_size" in _source,
  "keep_full_definition_with_plot_title": "--keep_full_definition_with_plot_title" in _source,
  "tick_label_font_size": "--tick_label_font_size" in _source,
  "circular_label_spacing": "--circular_label_spacing" in _source,
})
      `);
      circularMultiRecordCanvasSupportCache = JSON.parse(String(raw));
    } catch (_err) {
      circularMultiRecordCanvasSupportCache = {
        circular: false,
        multi_record_size_mode: false,
        multi_record_min_radius_ratio: false,
        multi_record_column_gap_ratio: false,
        multi_record_row_gap_ratio: false,
        multi_record_position: false,
        plot_title: false,
        plot_title_position: false,
        plot_title_font_size: false,
        keep_full_definition_with_plot_title: false,
        tick_label_font_size: false,
        circular_label_spacing: false
      };
    }
    return circularMultiRecordCanvasSupportCache;
  };

  const downloadLosatPair = async (pairIndex, customName) => {
    const entry = losatCacheInfo.value?.[pairIndex];
    const cacheMap = losatCache.value;
    if (!entry || !cacheMap) return;
    const cached = cacheMap.get(entry.key);
    if (!cached || typeof cached.text !== 'string') return;
    const defaultName = getLosatPairDefaultName(pairIndex);
    const filename = normalizeLosatFilename(
      customName,
      entry.filename || defaultName || `losat_pair_${pairIndex + 1}.tsv`
    );
    entry.filename = filename;
    downloadTextFile(filename, cached.text);
  };

  const setLosatPairFilename = (pairIndex, customName) => {
    const entry = losatCacheInfo.value?.[pairIndex];
    if (!entry) return;
    const defaultName = getLosatPairDefaultName(pairIndex);
    entry.filename = normalizeLosatFilename(
      customName,
      entry.filename || defaultName || `losat_pair_${pairIndex + 1}.tsv`
    );
  };

  const resetLabelScopeDialogState = () => {
    labelTextScopeDialog.show = false;
    labelTextScopeDialog.labelKey = '';
    labelTextScopeDialog.newText = '';
    labelTextScopeDialog.sourceText = '';
    labelTextScopeDialog.featureId = '';
    labelTextScopeDialog.matchingCount = 0;
  };

  const refreshCircularRecordOrder = async () => {
    if (!Array.isArray(adv.multi_record_positions)) {
      adv.multi_record_positions = [];
    }
    const pyodide = getPyodide();
    if (
      mode.value !== 'circular' ||
      cInputType.value !== 'gb' ||
      !files.c_gb ||
      !pyodideReady.value ||
      !pyodide
    ) {
      circularRecordList.value = [];
      if (!files.c_gb || cInputType.value !== 'gb') {
        adv.multi_record_positions.splice(0, adv.multi_record_positions.length);
      }
      return;
    }

    try {
      await writeFileToFs(files.c_gb, '/input.gb');
      const payloadRaw = pyodide.globals.get('list_genbank_records')('/input.gb');
      const payload = JSON.parse(String(payloadRaw || '{}'));
      if (payload?.error) {
        console.warn('Failed to read circular record list:', payload.error);
        circularRecordList.value = [];
        adv.multi_record_positions.splice(0, adv.multi_record_positions.length);
        return;
      }

      const nextRecords = [];
      const seenSelectors = new Set();
      (Array.isArray(payload?.records) ? payload.records : []).forEach((entry, index) => {
        const selector = String(entry?.selector ?? `#${index + 1}`).trim();
        if (!selector || seenSelectors.has(selector)) return;
        seenSelectors.add(selector);
        const recordId = String(entry?.record_id ?? '').trim() || `Record_${index + 1}`;
        nextRecords.push({ selector, record_id: recordId });
      });
      circularRecordList.value = nextRecords;
      const nextPositions = mergeCircularRecordPositions(nextRecords, adv.multi_record_positions);
      adv.multi_record_positions.splice(0, adv.multi_record_positions.length, ...nextPositions);
    } catch (error) {
      console.warn('Failed to refresh circular record order:', error);
      circularRecordList.value = [];
      adv.multi_record_positions.splice(0, adv.multi_record_positions.length);
    }
  };

  const runAnalysisInternal = async ({ runMode = 'manual', requestId = 0 } = {}) => {
    if (!pyodideReady.value) return { status: 'skipped' };
    const pyodide = getPyodide();
    if (!pyodide) return { status: 'skipped' };

    const isReflow = runMode === 'reflow';
    if (isReflow && mode.value === 'circular' && shouldDeferCircularPreviewUpdates.value) {
      return { status: 'skipped' };
    }
    const previousSelectedResultIndex = selectedResultIndex.value;
    const editableLabelsSnapshot = Array.isArray(editableLabels.value)
      ? editableLabels.value.map((entry) => ({ ...entry }))
      : [];
    const featureOverrideSourcesSnapshot = Object.fromEntries(
      Object.entries(labelTextFeatureOverrideSources || {}).map(([featureId, sourceText]) => [
        String(featureId || ''),
        String(sourceText ?? '')
      ])
    );
    const visibilityOverridesSnapshot = Object.fromEntries(
      Object.entries(labelVisibilityOverrides || {}).map(([featureId, modeValue]) => [
        String(featureId || ''),
        String(modeValue || '')
      ])
    );
    const activeRunColors = isReflow ? appliedPaletteColors.value : currentColors.value;

    if (isReflow) {
      labelReflowProcessing.value = true;
      labelReflowLastError.value = null;
      skipCaptureBaseConfig.value = true;
      skipPositionReapply.value = true;
    } else {
      processing.value = true;
      results.value = [];
      selectedResultIndex.value = 0;
      errorLog.value = null;
      if (typeof resetPreviewViewport === 'function') {
        resetPreviewViewport({ resetZoom: true });
      } else {
        zoom.value = 1.0;
      }
      skipCaptureBaseConfig.value = false;
      skipPositionReapply.value = false;
      pairwiseMatchFactors.value = {};
      clickedLabel.value = null;
      resetLabelScopeDialogState();
      addedLegendCaptions.value = new Set();
      fileLegendCaptions.value = new Set();
      Object.keys(featureColorOverrides).forEach((k) => delete featureColorOverrides[k]);
      legendEntries.value = [];
      deletedLegendEntries.value = [];
      Object.keys(legendColorOverrides).forEach((k) => delete legendColorOverrides[k]);
      originalLegendOrder.value = [];
      originalLegendColors.value = {};
      window._origPairwiseMin = activeRunColors.pairwise_match_min || '#FFE7E7';
      window._origPairwiseMax = activeRunColors.pairwise_match_max || '#FF7272';
    }
    labelOverrideBuildWarning.value = '';

    try {
      let args = [];
      let regionSpecs = [];
      let recordSelectors = [];
      let reverseFlags = [];
      let virtualBlastFiles = [];

      if (form.prefix && form.prefix.trim() !== '') args.push('-o', form.prefix.trim());
      if (form.species) args.push('--species', form.species);
      if (form.strain) args.push('--strain', form.strain);
      if (form.separate_strands) args.push('--separate_strands');

      if (adv.features.length) args.push('-k', adv.features.join(','));
      if (adv.window_size) args.push('--window', adv.window_size);
      if (adv.step_size) args.push('--step', adv.step_size);
      if (adv.nt && adv.nt !== 'GC') args.push('--nt', adv.nt);

      if (adv.def_font_size) args.push('--definition_font_size', adv.def_font_size);
      if (adv.label_font_size) args.push('--label_font_size', adv.label_font_size);

      if (adv.block_stroke_width !== null) args.push('--block_stroke_width', adv.block_stroke_width);
      if (adv.block_stroke_color) args.push('--block_stroke_color', adv.block_stroke_color);
      if (adv.line_stroke_width !== null) args.push('--line_stroke_width', adv.line_stroke_width);
      if (adv.line_stroke_color) args.push('--line_stroke_color', adv.line_stroke_color);
      if (adv.axis_stroke_width !== null) args.push('--axis_stroke_width', adv.axis_stroke_width);
      if (adv.axis_stroke_color) args.push('--axis_stroke_color', adv.axis_stroke_color);

      if (adv.legend_box_size) args.push('--legend_box_size', adv.legend_box_size);
      if (adv.legend_font_size) args.push('--legend_font_size', adv.legend_font_size);
      if (adv.resolve_overlaps) args.push('--resolve_overlaps');

      let dContent = '';
      for (const [k, v] of Object.entries(activeRunColors)) dContent += `${k}\t${v}\n`;
      pyodide.FS.writeFile('/combined_d.tsv', dContent);
      args.push('-d', '/combined_d.tsv');

      let tContent = '';
      manualSpecificRules.forEach((r) => {
        tContent += `${r.feat}\t${r.qual}\t${r.val}\t${r.color}\t${r.cap}\n`;
      });
      if (tContent.trim() !== '') {
        pyodide.FS.writeFile('/combined_t.tsv', tContent);
        args.push('-t', '/combined_t.tsv');
      }

      if (filterMode.value === 'Blacklist') {
        if (manualBlacklist.value) {
          args.push('--label_blacklist', manualBlacklist.value.replace(/\n/g, ','));
        }
      } else if (filterMode.value === 'Whitelist') {
        if (manualWhitelist.length > 0) {
          let wlContent = '';
          manualWhitelist.forEach((r) => {
            if (r.feat && r.qual) wlContent += `${r.feat}\t${r.qual}\t${r.key}\n`;
          });
          pyodide.FS.writeFile('/manual_wl.tsv', wlContent);
          args.push('--label_whitelist', '/manual_wl.tsv');
        }
      }

      let pContent = '';
      manualPriorityRules.forEach((r) => {
        pContent += `${r.feat}\t${r.order}\n`;
      });
      if (pContent.trim() !== '') {
        pyodide.FS.writeFile('/priority.tsv', pContent);
        args.push('--qualifier_priority', '/priority.tsv');
      }

      const labelOverride = buildLabelOverrideTsv(labelTextFeatureOverrides, labelTextBulkOverrides, {
        editableLabels: editableLabelsSnapshot,
        extractedFeatures: extractedFeatures.value,
        featureOverrideSources: featureOverrideSourcesSnapshot,
        visibilityOverrides: visibilityOverridesSnapshot
      });
      if (labelOverride.skippedMissingSourceCount > 0) {
        labelOverrideBuildWarning.value = `${labelOverride.skippedMissingSourceCount} feature override row(s) were skipped due to missing source label context.`;
      }
      if (labelOverride.tsv) {
        pyodide.FS.writeFile('/web_label_table.tsv', labelOverride.tsv);
        args.push('--label_table', '/web_label_table.tsv');
      }
      const featureVisibilityRows = [];
      Object.entries(featureVisibilityOverrides || {}).forEach(([featureIdRaw, modeRaw]) => {
        const featureId = String(featureIdRaw || '').trim();
        if (!featureId) return;
        const mode = normalizeFeatureVisibilityMode(modeRaw);
        if (mode === 'default') return;
        const action = mode === 'on' ? 'show' : 'hide';
        featureVisibilityRows.push(`*\t*\thash\t^${escapeRegexLiteral(featureId)}$\t${action}`);
      });
      if (featureVisibilityRows.length > 0) {
        pyodide.FS.writeFile('/web_feature_table.tsv', `${featureVisibilityRows.join('\n')}\n`);
        args.push('--feature_table', '/web_feature_table.tsv');
      }
      if (!isReflow) {
        editableLabels.value = [];
      }

      const selectedFeatureShapes = Array.isArray(adv.features)
        ? adv.features
            .map((featureTypeRaw) => {
              const featureType = String(featureTypeRaw || '').trim();
              if (!featureType) return null;
              const shape = normalizeFeatureShape(adv.feature_shapes?.[featureType]);
              return `${featureType}=${shape}`;
            })
            .filter((assignment) => typeof assignment === 'string' && assignment.length > 0)
        : [];

      if (mode.value === 'circular') {
        const multiCanvasSupport = getCircularMultiRecordCanvasOptionSupport();
        const normalizedCircularPlotTitle = String(form.plot_title || '').trim();
        const normalizedPlotTitlePosition = normalizeCircularPlotTitlePosition(adv.plot_title_position);
        const hasPlotTitleFontSize =
          adv.plot_title_font_size !== null &&
          adv.plot_title_font_size !== undefined &&
          adv.plot_title_font_size !== '';
        const parsedPlotTitleFontSize = hasPlotTitleFontSize
          ? Number(adv.plot_title_font_size)
          : null;
        const normalizedPlotTitleFontSize =
          parsedPlotTitleFontSize !== null &&
          Number.isFinite(parsedPlotTitleFontSize) &&
          parsedPlotTitleFontSize > 0
            ? parsedPlotTitleFontSize
            : null;
        const keepFullDefinitionWithPlotTitle = Boolean(adv.keep_full_definition_with_plot_title);
        form.plot_title = normalizedCircularPlotTitle;
        adv.plot_title_position = normalizedPlotTitlePosition;
        adv.plot_title_font_size = normalizedPlotTitleFontSize;
        adv.keep_full_definition_with_plot_title = keepFullDefinitionWithPlotTitle;

        if (selectedFeatureShapes.length > 0) {
          const shapeOptionSupport = getFeatureShapeOptionSupport();
          if (!shapeOptionSupport.circular) {
            throw new Error(
              'Current gbdraw wheel does not support --feature_shape. Rebuild and redeploy the web wheel.'
            );
          }
          selectedFeatureShapes.forEach((assignment) => {
            args.push('--feature_shape', assignment);
          });
        }
        args.push('--track_type', form.track_type, '-l', form.legend);
        const wantsCircularPlotTitleOption = normalizedCircularPlotTitle.length > 0;
        if (wantsCircularPlotTitleOption) {
          if (!multiCanvasSupport.plot_title) {
            throw new Error(
              'Current gbdraw wheel does not support --plot_title for circular diagrams. Rebuild and redeploy the web wheel.'
            );
          }
          args.push('--plot_title', normalizedCircularPlotTitle);
        }
        if (!multiCanvasSupport.plot_title_position) {
          throw new Error(
            'Current gbdraw wheel does not support circular plot title layout options. Rebuild and redeploy the web wheel.'
          );
        }
        args.push('--plot_title_position', normalizedPlotTitlePosition);
        if (
          normalizedPlotTitlePosition !== 'none' &&
          normalizedPlotTitleFontSize !== null &&
          Number.isFinite(normalizedPlotTitleFontSize)
        ) {
          if (!multiCanvasSupport.plot_title_font_size) {
            throw new Error(
              'Current gbdraw wheel does not support --plot_title_font_size. Rebuild and redeploy the web wheel.'
            );
          }
          args.push('--plot_title_font_size', String(normalizedPlotTitleFontSize));
        }
        if (keepFullDefinitionWithPlotTitle) {
          if (!multiCanvasSupport.keep_full_definition_with_plot_title) {
            throw new Error(
              'Current gbdraw wheel does not support --keep_full_definition_with_plot_title. Rebuild and redeploy the web wheel.'
            );
          }
          args.push('--keep_full_definition_with_plot_title');
        }
        const labelsModeRaw =
          typeof form.labels_mode === 'string'
            ? form.labels_mode
            : (form.allow_inner_labels ? 'both' : (form.show_labels ? 'out' : 'none'));
        const labelsMode = String(labelsModeRaw || 'none').trim().toLowerCase();
        if (labelsMode === 'out') args.push('--labels');
        if (labelsMode === 'both') args.push('--labels', 'both');
        if (form.suppress_gc) args.push('--suppress_gc');
        if (form.suppress_skew) args.push('--suppress_skew');
        if (form.multi_record_canvas) {
          if (!multiCanvasSupport.circular) {
            throw new Error(
              'Current gbdraw wheel does not support --multi_record_canvas. Rebuild and redeploy the web wheel.'
            );
          }
          if (
            !multiCanvasSupport.multi_record_size_mode ||
            !multiCanvasSupport.multi_record_min_radius_ratio ||
            !multiCanvasSupport.multi_record_column_gap_ratio ||
            !multiCanvasSupport.multi_record_row_gap_ratio
          ) {
            throw new Error(
              'Current gbdraw wheel does not support multi-record size/grid spacing options. Rebuild and redeploy the web wheel.'
            );
          }
          if (!multiCanvasSupport.plot_title_position) {
            throw new Error(
              'Current gbdraw wheel does not support multi-record plot-title layout options. Rebuild and redeploy the web wheel.'
            );
          }
          const effectiveRecordPositions = mergeCircularRecordPositions(
            circularRecordList.value,
            adv.multi_record_positions
          );
          const shouldPassRecordPositions = effectiveRecordPositions.length > 0;
          if (shouldPassRecordPositions && !multiCanvasSupport.multi_record_position) {
            throw new Error(
              'Current gbdraw wheel does not support --multi_record_position. Rebuild and redeploy the web wheel.'
            );
          }
          const normalizedSizeMode = normalizeMultiRecordSizeMode(adv.multi_record_size_mode);
          const normalizedMinRatio = normalizeMultiRecordMinRadiusRatio(adv.multi_record_min_radius_ratio);
          const normalizedColumnGapRatio = normalizeMultiRecordColumnGapRatio(adv.multi_record_column_gap_ratio);
          const normalizedRowGapRatio = normalizeMultiRecordRowGapRatio(adv.multi_record_row_gap_ratio);
          adv.multi_record_size_mode = normalizedSizeMode;
          adv.multi_record_min_radius_ratio = normalizedMinRatio;
          adv.multi_record_column_gap_ratio = normalizedColumnGapRatio;
          adv.multi_record_row_gap_ratio = normalizedRowGapRatio;
          adv.multi_record_positions.splice(
            0,
            adv.multi_record_positions.length,
            ...effectiveRecordPositions
          );
          adv.plot_title_position = normalizedPlotTitlePosition;
          adv.plot_title_font_size = normalizedPlotTitleFontSize;
          args.push('--multi_record_canvas');
          args.push('--multi_record_size_mode', normalizedSizeMode);
          args.push('--multi_record_min_radius_ratio', String(normalizedMinRatio));
          args.push('--multi_record_column_gap_ratio', String(normalizedColumnGapRatio));
          args.push('--multi_record_row_gap_ratio', String(normalizedRowGapRatio));
          if (shouldPassRecordPositions) {
            effectiveRecordPositions.forEach((entry) => {
              const token = buildMultiRecordPositionToken(entry);
              if (!token) return;
              args.push('--multi_record_position', token);
            });
          }
        }

        if (adv.outer_label_x_offset) args.push('--outer_label_x_radius_offset', adv.outer_label_x_offset);
        if (adv.outer_label_y_offset) args.push('--outer_label_y_radius_offset', adv.outer_label_y_offset);
        if (adv.inner_label_x_offset) args.push('--inner_label_x_radius_offset', adv.inner_label_x_offset);
        if (adv.inner_label_y_offset) args.push('--inner_label_y_radius_offset', adv.inner_label_y_offset);
        if (
          adv.circular_label_spacing !== null &&
          adv.circular_label_spacing !== undefined &&
          adv.circular_label_spacing !== ''
        ) {
          if (!multiCanvasSupport.circular_label_spacing) {
            throw new Error(
              'Current gbdraw wheel does not support --circular_label_spacing. Rebuild and redeploy the web wheel.'
            );
          }
          args.push('--circular_label_spacing', adv.circular_label_spacing);
        }
        if (
          adv.feature_width_circular !== null &&
          adv.feature_width_circular !== undefined &&
          adv.feature_width_circular !== '' &&
          Number(adv.feature_width_circular) > 0
        ) {
          args.push('--feature_width', adv.feature_width_circular);
        }
        if (!form.suppress_gc) {
          if (
            adv.gc_content_width_circular !== null &&
            adv.gc_content_width_circular !== undefined &&
            adv.gc_content_width_circular !== '' &&
            Number(adv.gc_content_width_circular) > 0
          ) {
            args.push('--gc_content_width', adv.gc_content_width_circular);
          }
          if (
            adv.gc_content_radius_circular !== null &&
            adv.gc_content_radius_circular !== undefined &&
            adv.gc_content_radius_circular !== '' &&
            Number(adv.gc_content_radius_circular) > 0
          ) {
            args.push('--gc_content_radius', adv.gc_content_radius_circular);
          }
        }
        if (!form.suppress_skew) {
          if (
            adv.gc_skew_width_circular !== null &&
            adv.gc_skew_width_circular !== undefined &&
            adv.gc_skew_width_circular !== '' &&
            Number(adv.gc_skew_width_circular) > 0
          ) {
            args.push('--gc_skew_width', adv.gc_skew_width_circular);
          }
          if (
            adv.gc_skew_radius_circular !== null &&
            adv.gc_skew_radius_circular !== undefined &&
            adv.gc_skew_radius_circular !== '' &&
            Number(adv.gc_skew_radius_circular) > 0
          ) {
            args.push('--gc_skew_radius', adv.gc_skew_radius_circular);
          }
        }
        if (adv.scale_interval) args.push('--scale_interval', adv.scale_interval);
        if (
          adv.tick_label_font_size !== null &&
          adv.tick_label_font_size !== undefined &&
          adv.tick_label_font_size !== ''
        ) {
          if (!multiCanvasSupport.tick_label_font_size) {
            throw new Error(
              'Current gbdraw wheel does not support --tick_label_font_size. Rebuild and redeploy the web wheel.'
            );
          }
          args.push('--tick_label_font_size', adv.tick_label_font_size);
        }

        if (cInputType.value === 'gb') {
          if (!files.c_gb) throw new Error('Please upload a GenBank file.');
          await writeFileToFs(files.c_gb, '/input.gb');
          args.push('--gbk', '/input.gb');
        } else {
          if (!files.c_gff || !files.c_fasta) throw new Error('GFF3 and FASTA are required.');
          await writeFileToFs(files.c_gff, '/input.gff');
          await writeFileToFs(files.c_fasta, '/input.fasta');
          args.push('--gff', '/input.gff', '--fasta', '/input.fasta');
        }
      } else {
        if (selectedFeatureShapes.length > 0) {
          const shapeOptionSupport = getFeatureShapeOptionSupport();
          if (!shapeOptionSupport.linear) {
            throw new Error(
              'Current gbdraw wheel does not support --feature_shape. Rebuild and redeploy the web wheel.'
            );
          }
          selectedFeatureShapes.forEach((assignment) => {
            args.push('--feature_shape', assignment);
          });
        }
        args.push('--scale_style', form.scale_style);
        const normalizedTrackLayout =
          form.linear_track_layout === 'spreadout'
            ? 'above'
            : form.linear_track_layout === 'tuckin'
              ? 'below'
              : (form.linear_track_layout || 'middle');
        if (form.align_center) args.push('--align_center');
        if (form.show_gc) args.push('--show_gc');
        if (form.show_skew) args.push('--show_skew');
        if (form.normalize_length) args.push('--normalize_length');
        if (form.legend !== 'right') args.push('-l', form.legend);
        adv.min_bitscore = normalizeBlastThresholdNumber(
          adv.min_bitscore,
          DEFAULT_LINEAR_BLAST_FILTERS.bitscore
        );
        adv.evalue = normalizeBlastThresholdText(adv.evalue, DEFAULT_LINEAR_BLAST_FILTERS.evalue);
        adv.identity = normalizeBlastThresholdNumber(
          adv.identity,
          DEFAULT_LINEAR_BLAST_FILTERS.identity
        );
        adv.alignment_length = normalizeBlastThresholdNumber(
          adv.alignment_length,
          DEFAULT_LINEAR_BLAST_FILTERS.alignment_length,
          { integer: true }
        );
        args.push(
          '--bitscore',
          adv.min_bitscore,
          '--evalue',
          adv.evalue,
          '--identity',
          adv.identity,
          '--alignment_length',
          adv.alignment_length
        );

        const normalizedPlotTitle = String(form.plot_title || '').trim();
        const normalizedPlotTitlePosition = normalizeLinearPlotTitlePosition(adv.plot_title_position);
        const hasPlotTitleFontSize =
          adv.plot_title_font_size !== null &&
          adv.plot_title_font_size !== undefined &&
          adv.plot_title_font_size !== '';
        const parsedPlotTitleFontSize = hasPlotTitleFontSize ? Number(adv.plot_title_font_size) : null;
        const normalizedPlotTitleFontSize =
          parsedPlotTitleFontSize !== null &&
          Number.isFinite(parsedPlotTitleFontSize) &&
          parsedPlotTitleFontSize > 0
            ? parsedPlotTitleFontSize
            : null;
        adv.linear_show_replicon = adv.linear_show_replicon === true;
        adv.linear_show_accession = adv.linear_show_accession !== false;
        adv.linear_show_length = adv.linear_show_length !== false;
        form.plot_title = normalizedPlotTitle;
        adv.plot_title_position = normalizedPlotTitlePosition;
        adv.plot_title_font_size = normalizedPlotTitleFontSize;

        if (form.show_labels_linear !== 'none') {
          args.push('--show_labels');
          if (form.show_labels_linear === 'first') args.push('first');
        }
        const normalizedLabelPlacement = adv.label_placement === 'on_feature' ? 'above_feature' : adv.label_placement;
        const wantsPlacementOption = normalizedLabelPlacement && normalizedLabelPlacement !== 'auto';
        const wantsRotationOption = adv.label_rotation !== null && adv.label_rotation !== undefined && adv.label_rotation !== '';
        const wantsLinearLabelSpacingOption =
          adv.linear_label_spacing !== null && adv.linear_label_spacing !== undefined && adv.linear_label_spacing !== '';
        const wantsTrackLayoutOption = normalizedTrackLayout !== 'middle';
        const wantsTrackAxisGapOption =
          adv.track_axis_gap !== null && adv.track_axis_gap !== undefined && adv.track_axis_gap !== '';
        const wantsRulerLabelFontOption =
          adv.scale_font_size !== null && adv.scale_font_size !== undefined && adv.scale_font_size !== '';
        const wantsRulerLabelColorOption =
          adv.ruler_label_color !== null &&
          adv.ruler_label_color !== undefined &&
          String(adv.ruler_label_color).trim() !== '';
        const wantsPlotTitleOption = normalizedPlotTitle !== '';
        const wantsPlotTitlePositionOption = normalizedPlotTitlePosition !== 'bottom';
        const wantsPlotTitleFontSizeOption = normalizedPlotTitleFontSize !== null;
        const wantsShowRepliconOption = adv.linear_show_replicon === true;
        const wantsHideAccessionOption = adv.linear_show_accession === false;
        const wantsHideLengthOption = adv.linear_show_length === false;
        const wantsRulerOnAxisOption =
          Boolean(form.linear_ruler_on_axis) &&
          form.scale_style === 'ruler' &&
          (normalizedTrackLayout === 'above' || normalizedTrackLayout === 'below');
        const linearLabelSupport = getLinearLabelOptionSupport();
        if (
          wantsPlacementOption ||
          wantsRotationOption ||
          wantsLinearLabelSpacingOption ||
          wantsTrackLayoutOption ||
          wantsTrackAxisGapOption ||
          wantsRulerOnAxisOption ||
          wantsRulerLabelColorOption ||
          wantsPlotTitleOption ||
          wantsPlotTitlePositionOption ||
          wantsPlotTitleFontSizeOption ||
          wantsShowRepliconOption ||
          wantsHideAccessionOption ||
          wantsHideLengthOption
        ) {
          if (wantsPlacementOption && !linearLabelSupport.placement) {
            throw new Error("Current gbdraw wheel does not support --label_placement. Rebuild and redeploy the web wheel.");
          }
          if (wantsRotationOption && !linearLabelSupport.rotation) {
            throw new Error("Current gbdraw wheel does not support --label_rotation. Rebuild and redeploy the web wheel.");
          }
          if (wantsLinearLabelSpacingOption && !linearLabelSupport.linear_label_spacing) {
            throw new Error("Current gbdraw wheel does not support --linear_label_spacing. Rebuild and redeploy the web wheel.");
          }
          if (wantsTrackLayoutOption && !linearLabelSupport.track_layout) {
            throw new Error("Current gbdraw wheel does not support --track_layout. Rebuild and redeploy the web wheel.");
          }
          if (wantsTrackAxisGapOption && !linearLabelSupport.track_axis_gap) {
            throw new Error("Current gbdraw wheel does not support --track_axis_gap. Rebuild and redeploy the web wheel.");
          }
          if (wantsRulerOnAxisOption && !linearLabelSupport.ruler_on_axis) {
            throw new Error("Current gbdraw wheel does not support --ruler_on_axis. Rebuild and redeploy the web wheel.");
          }
          if (wantsRulerLabelColorOption && !linearLabelSupport.ruler_label_color) {
            throw new Error("Current gbdraw wheel does not support --ruler_label_color. Rebuild and redeploy the web wheel.");
          }
          if (wantsPlotTitleOption && !linearLabelSupport.plot_title) {
            throw new Error("Current gbdraw wheel does not support --plot_title. Rebuild and redeploy the web wheel.");
          }
          if (wantsPlotTitlePositionOption && !linearLabelSupport.plot_title_position) {
            throw new Error("Current gbdraw wheel does not support --plot_title_position. Rebuild and redeploy the web wheel.");
          }
          if (wantsPlotTitleFontSizeOption && !linearLabelSupport.plot_title_font_size) {
            throw new Error("Current gbdraw wheel does not support --plot_title_font_size. Rebuild and redeploy the web wheel.");
          }
          if (wantsShowRepliconOption && !linearLabelSupport.show_replicon) {
            throw new Error("Current gbdraw wheel does not support --show_replicon. Rebuild and redeploy the web wheel.");
          }
          if (wantsHideAccessionOption && !linearLabelSupport.hide_accession) {
            throw new Error("Current gbdraw wheel does not support --hide_accession. Rebuild and redeploy the web wheel.");
          }
          if (wantsHideLengthOption && !linearLabelSupport.hide_length) {
            throw new Error("Current gbdraw wheel does not support --hide_length. Rebuild and redeploy the web wheel.");
          }
        }
        if (wantsPlotTitleOption) args.push('--plot_title', normalizedPlotTitle);
        if (wantsPlotTitlePositionOption) args.push('--plot_title_position', normalizedPlotTitlePosition);
        if (wantsPlotTitleFontSizeOption) args.push('--plot_title_font_size', String(normalizedPlotTitleFontSize));
        if (wantsShowRepliconOption) args.push('--show_replicon');
        if (wantsHideAccessionOption) args.push('--hide_accession');
        if (wantsHideLengthOption) args.push('--hide_length');
        if (normalizedLabelPlacement && normalizedLabelPlacement !== 'auto') {
          args.push('--label_placement', normalizedLabelPlacement);
        }
        if (adv.label_rotation !== null && adv.label_rotation !== undefined && adv.label_rotation !== '') {
          args.push('--label_rotation', adv.label_rotation);
        }
        if (adv.linear_label_spacing !== null && adv.linear_label_spacing !== undefined && adv.linear_label_spacing !== '') {
          args.push('--linear_label_spacing', adv.linear_label_spacing);
        }
        if (normalizedTrackLayout !== 'middle') {
          args.push('--track_layout', normalizedTrackLayout);
        }
        if (adv.track_axis_gap !== null && adv.track_axis_gap !== undefined && adv.track_axis_gap !== '') {
          args.push('--track_axis_gap', adv.track_axis_gap);
        }
        if (
          Boolean(form.linear_ruler_on_axis) &&
          form.scale_style === 'ruler' &&
          (normalizedTrackLayout === 'above' || normalizedTrackLayout === 'below')
        ) {
          args.push('--ruler_on_axis');
        }

        if (adv.feature_height) args.push('--feature_height', adv.feature_height);
        if (adv.gc_height) args.push('--gc_height', adv.gc_height);
        if (adv.comparison_height) args.push('--comparison_height', adv.comparison_height);

        if (adv.scale_interval) args.push('--scale_interval', adv.scale_interval);
        if (wantsRulerLabelFontOption) {
          if (form.scale_style === 'ruler') {
            if (linearLabelSupport.ruler_label_font_size) {
              args.push('--ruler_label_font_size', adv.scale_font_size);
            } else if (linearLabelSupport.scale_font_size) {
              args.push('--scale_font_size', adv.scale_font_size);
            } else {
              throw new Error("Current gbdraw wheel does not support ruler label font options. Rebuild and redeploy the web wheel.");
            }
          } else if (linearLabelSupport.scale_font_size) {
            args.push('--scale_font_size', adv.scale_font_size);
          } else {
            throw new Error("Current gbdraw wheel does not support --scale_font_size. Rebuild and redeploy the web wheel.");
          }
        }
        if (wantsRulerLabelColorOption) args.push('--ruler_label_color', adv.ruler_label_color);
        if (adv.scale_stroke_width) args.push('--scale_stroke_width', adv.scale_stroke_width);
        if (adv.scale_stroke_color) args.push('--scale_stroke_color', adv.scale_stroke_color);

        const recordLabels = linearSeqs.map((seq) => (seq.definition ?? '').toString());
        const hasRecordLabels = recordLabels.some((label) => label.trim() !== '');
        if (hasRecordLabels) {
          recordLabels.forEach((label) => {
            args.push('--record_label', label);
          });
        }

        const buildRegionSpec = (seq, idx) => {
          const hasStart = seq.region_start !== null && seq.region_start !== undefined && seq.region_start !== '';
          const hasEnd = seq.region_end !== null && seq.region_end !== undefined && seq.region_end !== '';
          const recordIdRaw = seq.region_record_id ? String(seq.region_record_id).trim() : '';
          const wantsReverse = Boolean(seq.region_reverse);
          if (hasStart !== hasEnd) {
            throw new Error(`Sequence #${idx + 1}: Provide both Region start and end, or leave both empty.`);
          }

          recordSelectors.push(recordIdRaw || '');

          let reverseFlag = wantsReverse;

          if (hasStart && hasEnd) {
            const start = Number(seq.region_start);
            const end = Number(seq.region_end);
            if (!Number.isFinite(start) || !Number.isFinite(end)) {
              throw new Error(`Sequence #${idx + 1}: Region start/end must be numbers.`);
            }
            if (!Number.isInteger(start) || !Number.isInteger(end)) {
              throw new Error(`Sequence #${idx + 1}: Region start/end must be integers.`);
            }
            if (start < 1 || end < 1) {
              throw new Error(`Sequence #${idx + 1}: Region start/end must be >= 1.`);
            }
            const specBody = `${start}-${end}${wantsReverse ? ':rc' : ''}`;
            const fileLabel = lInputType.value === 'gb' ? `seq_${idx}.gb` : `seq_${idx}.gff`;
            const cliSpec = recordIdRaw ? `${fileLabel}:${recordIdRaw}:${specBody}` : `#${idx + 1}:${specBody}`;
            const fileSpec = specBody;
            reverseFlag = false;
            reverseFlags.push(reverseFlag);
            return { cli: cliSpec, file: fileSpec };
          }

          reverseFlags.push(reverseFlag);
          return null;
        };

        let inputArgs = [];
        let blastArgs = [];
        const useLosat = blastSource.value === 'losat';
        const fastaCache = new Map();
        const fastaHashCache = new Map();
        let extractFirstFasta = null;
        let cacheInfo = [];
        const cacheMap = losatCache.value || new Map();
        const losatTiming = useLosat
          ? {
              inputWriteMs: 0,
              fastaExtractionMs: 0,
              cacheHashMs: 0,
              jobBuildMs: 0,
              executionMs: 0,
              blastWriteMs: 0,
              totalFastaChars: 0,
              cacheHits: 0,
              cacheMisses: 0,
              totalPairs: 0,
              uniqueJobs: 0
            }
          : null;

        if (useLosat) {
          extractFirstFasta = pyodide.globals.get('extract_first_fasta');
        } else {
          losatCacheInfo.value = [];
        }

        const getSeqEntry = (idx) => {
          if (fastaCache.has(idx)) return fastaCache.get(idx);
          const startedAt = getNow();
          const path = lInputType.value === 'gb' ? `/seq_${idx}.gb` : `/seq_${idx}.fasta`;
          const fmt = lInputType.value === 'gb' ? 'genbank' : 'fasta';
          const regionSpec = regionSpecs[idx]?.file || null;
          const recordSelector = recordSelectors[idx] ?? '';
          const reverseFlag = reverseFlags[idx] ? '1' : '0';
          const res = JSON.parse(extractFirstFasta(path, fmt, regionSpec, recordSelector, reverseFlag));
          if (res.error) throw new Error(res.error);
          const entry = {
            fasta: res.fasta,
            recordId: res.record_id || `seq_${idx + 1}`
          };
          fastaCache.set(idx, entry);
          if (losatTiming) {
            losatTiming.fastaExtractionMs += getNow() - startedAt;
            losatTiming.totalFastaChars += entry.fasta.length;
          }
          return entry;
        };

        const getSeqHash = async (idx) => {
          if (fastaHashCache.has(idx)) return fastaHashCache.get(idx);
          const entry = getSeqEntry(idx);
          const startedAt = getNow();
          const hash = await hashText(entry.fasta);
          if (losatTiming) losatTiming.cacheHashMs += getNow() - startedAt;
          fastaHashCache.set(idx, hash);
          return hash;
        };

        const buildCacheKey = async (argsKey, queryIdx, subjectIdx) => {
          const queryHash = await getSeqHash(queryIdx);
          const subjectHash = await getSeqHash(subjectIdx);
          const payload = JSON.stringify({
            program: losatProgram.value,
            outfmt: String(losat.outfmt || '6'),
            args: argsKey,
            queryHash,
            subjectHash
          });
          return hashText(payload);
        };

        const buildCacheFilename = (pairIndex, queryEntry, subjectEntry) =>
          getLosatPairDefaultName(pairIndex, queryEntry, subjectEntry);

        const pushArg = (arr, flag, value) => {
          if (value === null || value === undefined || value === '') return;
          if (typeof value === 'number' && !Number.isFinite(value)) return;
          const valueStr = String(value);
          if (valueStr.startsWith('-')) {
            arr.push(`${flag}=${valueStr}`);
          } else {
            arr.push(flag, valueStr);
          }
        };

        const getGencode = (idx) => {
          const raw = linearSeqs[idx]?.losat_gencode;
          if (raw === null || raw === undefined || raw === '') return null;
          const num = Number(raw);
          if (!Number.isFinite(num)) return null;
          return num;
        };

        const buildLosatArgs = (queryIdx, subjectIdx) => {
          const args = [];
          if (losatProgram.value === 'blastn') {
            pushArg(args, '--task', losat.blastn.task);
          } else {
            pushArg(args, '--query-gencode', getGencode(queryIdx));
            pushArg(args, '--db-gencode', getGencode(subjectIdx));
          }
          return args;
        };

        {
          const inputWriteStartedAt = getNow();
          for (let i = 0; i < linearSeqs.length; i++) {
            const seq = linearSeqs[i];
            if (lInputType.value === 'gb') {
              if (!seq.gb) throw new Error(`Sequence #${i + 1}: Missing GenBank file.`);
              await writeFileToFs(seq.gb, `/seq_${i}.gb`);
              inputArgs.push(`/seq_${i}.gb`);
            } else {
              if (!seq.gff || !seq.fasta) throw new Error(`Sequence #${i + 1}: GFF3 and FASTA are required.`);
              await writeFileToFs(seq.gff, `/seq_${i}.gff`);
              await writeFileToFs(seq.fasta, `/seq_${i}.fasta`);
            }
          }
          if (losatTiming) losatTiming.inputWriteMs += getNow() - inputWriteStartedAt;
        }

        regionSpecs = linearSeqs.map((seq, idx) => buildRegionSpec(seq, idx));
        regionSpecs.forEach((spec) => {
          if (spec?.cli) args.push('--region', spec.cli);
        });

        recordSelectors.forEach((selector) => {
          args.push('--record_id', selector);
        });
        reverseFlags.forEach((flag) => {
          args.push('--reverse_complement', flag ? '1' : '0');
        });

        if (useLosat) {
          const losatPairs = [];
          const losatJobs = [];
          const pendingJobKeys = new Set();
          const jobBuildStartedAt = getNow();

          for (let i = 0; i < linearSeqs.length - 1; i++) {
            const queryEntry = getSeqEntry(i);
            const subjectEntry = getSeqEntry(i + 1);
            const losatArgs = buildLosatArgs(i, i + 1);
            const cacheKey = await buildCacheKey(losatArgs, i, i + 1);
            const cached = cacheMap.get(cacheKey);
            const hasCachedText = typeof cached?.text === 'string';
            losatTiming.totalPairs += 1;
            if (hasCachedText) losatTiming.cacheHits += 1;
            else losatTiming.cacheMisses += 1;
            const pair = {
              pairIndex: i,
              cacheKey,
              filename: buildCacheFilename(i, queryEntry, subjectEntry)
            };
            losatPairs.push(pair);
            cacheInfo.push({
              key: cacheKey,
              filename: pair.filename
            });

            if (!hasCachedText && !pendingJobKeys.has(cacheKey)) {
              pendingJobKeys.add(cacheKey);
              losatJobs.push({
                pairIndex: i,
                cacheKey,
                program: losatProgram.value,
                queryFasta: queryEntry.fasta,
                subjectFasta: subjectEntry.fasta,
                outfmt: losat.outfmt || '6',
                extraArgs: losatArgs
              });
            }
          }
          losatTiming.uniqueJobs = losatJobs.length;
          losatTiming.jobBuildMs += getNow() - jobBuildStartedAt;

          if (losatJobs.length > 0) {
            const executionStartedAt = getNow();
            const losatResults = await runLosatPairsParallel(losatJobs, {
              concurrency: getLosatParallelWorkers()
            });
            losatTiming.executionMs += getNow() - executionStartedAt;
            losatResults.forEach((result) => {
              cacheMap.set(result.cacheKey, { text: result.text });
            });
          }

          const blastWriteStartedAt = getNow();
          losatPairs.forEach((pair) => {
            const cached = cacheMap.get(pair.cacheKey);
            const blastText = typeof cached?.text === 'string' ? cached.text : '';
            const blastPath = `/blast_${pair.pairIndex}.txt`;
            virtualBlastFiles.push({ path: blastPath, text: blastText });
            blastArgs.push(blastPath);
          });
          losatTiming.blastWriteMs += getNow() - blastWriteStartedAt;
          console.info(
            [
              `LOSAT timing: pairs=${losatTiming.totalPairs}`,
              `cache hits=${losatTiming.cacheHits}`,
              `misses=${losatTiming.cacheMisses}`,
              `unique jobs=${losatTiming.uniqueJobs}`,
              `input FS write=${formatDuration(losatTiming.inputWriteMs)}`,
              `FASTA extraction=${formatDuration(losatTiming.fastaExtractionMs)}`,
              `cache hashing=${formatDuration(losatTiming.cacheHashMs)}`,
              `job build=${formatDuration(losatTiming.jobBuildMs)}`,
              `execution=${formatDuration(losatTiming.executionMs)}`,
              `BLAST FS write=${formatDuration(losatTiming.blastWriteMs)}`,
              `FASTA chars=${losatTiming.totalFastaChars.toLocaleString()}`
            ].join(', ')
          );
        } else {
          for (let i = 0; i < linearSeqs.length - 1; i++) {
            const seq = linearSeqs[i];
            if (seq.blast) {
              await writeFileToFs(seq.blast, `/blast_${i}.txt`);
              blastArgs.push(`/blast_${i}.txt`);
            }
          }
        }
        if (extractFirstFasta) {
          extractFirstFasta.destroy();
        }
        if (useLosat) {
          losatCacheInfo.value = cacheInfo;
          losatCache.value = cacheMap;
        }
        if (lInputType.value === 'gb') args.push('--gbk', ...inputArgs);
        else {
          let gffs = [];
          let fastas = [];
          for (let i = 0; i < linearSeqs.length; i++) {
            gffs.push(`/seq_${i}.gff`);
            fastas.push(`/seq_${i}.fasta`);
          }
          args.push('--gff', ...gffs, '--fasta', ...fastas);
        }
        if (blastArgs.length) args.push('-b', ...blastArgs);
      }

      console.log('CMD:', args.join(' '));
      const gbdrawStartedAt = getNow();
      const jsonResult = pyodide
        .globals
        .get('run_gbdraw_wrapper')(
          mode.value,
          pyodide.toPy(args.map(String)),
          virtualBlastFiles.length ? JSON.stringify(virtualBlastFiles) : null
        );
      console.info(`gbdraw ${mode.value} wrapper execution: ${formatDuration(getNow() - gbdrawStartedAt)}.`);
      const res = JSON.parse(jsonResult);
      if (res.error) {
        if (isReflow) {
          labelReflowLastError.value = formatPythonError(res.error)?.summary || 'Auto reflow failed';
          return { status: 'error' };
        }
        errorLog.value = formatPythonError(res.error);
        return { status: 'error' };
      }

      if (isReflow && requestId !== pendingReflowRequestId) {
        return { status: 'stale' };
      }

      if (isReflow) {
        skipCaptureBaseConfig.value = true;
        skipPositionReapply.value = true;
      }

      results.value = res;
      if (isReflow) {
        if (res.length > 0) {
          const safeIndex = Math.max(0, Math.min(previousSelectedResultIndex, res.length - 1));
          selectedResultIndex.value = safeIndex;
        }
      }

      generatedLegendPosition.value = form.legend;
      generatedMode.value = mode.value;
      generatedMultiRecordCanvas.value = mode.value === 'circular' ? Boolean(form.multi_record_canvas) : false;
      if (mode.value === 'circular') {
        generatedCircularPlotTitlePosition.value = normalizeCircularPlotTitlePosition(adv.plot_title_position);
      }

      extractedFeatures.value = [];

      if (mode.value === 'circular' && cInputType.value === 'gb') {
        try {
          const featJson = pyodide.globals.get('extract_features_from_genbank')('/input.gb', null, null, null, null);
          const featData = JSON.parse(featJson);
          if (featData?.error) {
            console.warn('Feature extraction failed for circular input:', featData.error);
          } else if (Array.isArray(featData?.features)) {
            extractedFeatures.value = featData.features;
            featureRecordIds.value = featData.record_ids || [];
            selectedFeatureRecordIdx.value = 0;
            refreshFeatureOverrides(featData.features);
            console.log(
              `Extracted ${featData.features.length} features from ${featData.record_ids.length} record(s) for color editor.`
            );
          } else {
            console.warn('Feature extraction returned an unexpected payload for circular input.', featData);
          }
        } catch (e) {
          console.log('Could not extract features:', e);
        }
      } else if (mode.value === 'linear' && lInputType.value === 'gb' && linearSeqs.length > 0) {
        try {
          let allFeatures = [];
          let allRecordLabels = [];
          for (let i = 0; i < linearSeqs.length; i++) {
            const regionSpec = regionSpecs[i]?.file || null;
            const recordSelector = recordSelectors[i] ?? '';
            const reverseFlag = reverseFlags[i] ? '1' : '0';
            const featJson = pyodide
              .globals
              .get('extract_features_from_genbank')(
                `/seq_${i}.gb`,
                regionSpec,
                recordSelector,
                reverseFlag,
                null
              );
            const featData = JSON.parse(featJson);
            if (featData?.error) {
              console.warn(`Feature extraction failed for linear input #${i + 1}:`, featData.error);
            } else if (Array.isArray(featData?.features)) {
              featData.features.forEach((f) => {
                f.fileIdx = i;
                f.displayRecordId = `File ${i + 1}: ${f.record_id}`;
                f.id = `file${i}_${f.id}`;
              });
              allFeatures = allFeatures.concat(featData.features);
              featData.record_ids.forEach((rid, ridx) => {
                allRecordLabels.push({ label: `File ${i + 1}: ${rid}`, fileIdx: i, recordIdx: ridx });
              });
            } else {
              console.warn(`Feature extraction returned an unexpected payload for linear input #${i + 1}.`, featData);
            }
          }
          extractedFeatures.value = allFeatures;
          featureRecordIds.value = allRecordLabels.map((r) => r.label);
          selectedFeatureRecordIdx.value = 0;
          refreshFeatureOverrides(allFeatures);
          console.log(
            `Extracted ${allFeatures.length} features from ${linearSeqs.length} file(s) for color editor.`
          );
        } catch (e) {
          console.log('Could not extract features:', e);
        }
      }
      if (!isReflow) {
        appliedPaletteName.value = String(selectedPalette?.value || appliedPaletteName.value || 'default');
        appliedPaletteColors.value = { ...currentColors.value };
        pendingPaletteName.value = '';
        pendingPaletteColors.value = {};
      }
      return { status: 'ok' };
    } catch (e) {
      if (isReflow) {
        labelReflowLastError.value = formatJsError(e)?.summary || 'Auto reflow failed';
        return { status: 'error' };
      }
      errorLog.value = formatJsError(e);
      return { status: 'error' };
    } finally {
      if (isReflow) {
        labelReflowProcessing.value = false;
      } else {
        processing.value = false;
      }
    }
  };

  const runAnalysis = async () => runAnalysisInternal({ runMode: 'manual' });

  const runLabelReflow = async (reason = 'label-edit') => {
    pendingReflowRequestId += 1;
    pendingReflowReason = String(reason || 'label-edit');
    if (activeReflowRequestId !== 0) return;

    while (activeReflowRequestId < pendingReflowRequestId) {
      activeReflowRequestId = pendingReflowRequestId;
      await runAnalysisInternal({
        runMode: 'reflow',
        requestId: activeReflowRequestId,
        reason: pendingReflowReason
      });
    }

    activeReflowRequestId = 0;
  };

  const downloadLosatCache = async () => {
    if (!losatCacheInfo.value || losatCacheInfo.value.length === 0) return;
    const cacheMap = losatCache.value;
    if (!cacheMap || cacheMap.size === 0) return;

    const totalChars = losatCacheInfo.value.reduce((sum, entry) => {
      const cached = cacheMap.get(entry.key);
      return sum + (cached?.text ? cached.text.length : 0);
    }, 0);

    if (totalChars > 50 * 1024 * 1024) {
      const proceed = confirm(
        `LOSAT TSV export will download about ${(totalChars / (1024 * 1024)).toFixed(1)} MB. Continue?`
      );
      if (!proceed) return;
    }

    for (let idx = 0; idx < losatCacheInfo.value.length; idx += 1) {
      const entry = losatCacheInfo.value[idx];
      const cached = cacheMap.get(entry.key);
      if (!cached) continue;
      const filename = entry.filename || `losat_pair_${idx + 1}.tsv`;
      downloadTextFile(filename, cached.text);
      await new Promise((resolve) => setTimeout(resolve, 0));
    }
  };

  const clearLosatCache = () => {
    if (losatCache.value) {
      losatCache.value.clear();
    }
    losatCacheInfo.value = [];
  };

  return {
    runAnalysis,
    runLabelReflow,
    refreshCircularRecordOrder,
    downloadLosatCache,
    downloadLosatPair,
    setLosatPairFilename,
    clearLosatCache,
    getLosatPairDefaultName
  };
};
