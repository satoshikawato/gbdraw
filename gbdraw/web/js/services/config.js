import { state, normalizeLinearSeqList, collapseEmptyLinearSeqList } from '../state.js';
import { resolveColorToHex } from '../app/color-utils.js';

const SESSION_VERSION = 7;

const safeDeepMerge = (target, source) => {
  if (!source || typeof source !== 'object') return;

  Object.keys(source).forEach((key) => {
    // 1. Prevent prototype pollution
    if (['__proto__', 'constructor', 'prototype'].includes(key)) return;

    // 2. Ignore keys not present in target (whitelisting effect)
    if (!Object.prototype.hasOwnProperty.call(target, key)) return;

    const targetValue = target[key];
    const sourceValue = source[key];

    // 3. Recursive merge for objects
    if (
      targetValue &&
      typeof targetValue === 'object' &&
      !Array.isArray(targetValue) &&
      sourceValue &&
      typeof sourceValue === 'object' &&
      !Array.isArray(sourceValue)
    ) {
      safeDeepMerge(targetValue, sourceValue);
      return;
    }

    // 4. For arrays, intentionally overwrite (replacing lists of settings is natural)
    if (Array.isArray(targetValue) && Array.isArray(sourceValue)) {
      target[key].splice(0, target[key].length, ...sourceValue);
      return;
    }

    // 5. Update value only if types match or initial value is null
    if (typeof targetValue === typeof sourceValue || targetValue === null) {
      target[key] = sourceValue;
    }
  });
};

const downloadJson = (data, filename) => {
  const blob = new Blob([JSON.stringify(data, null, 2)], {
    type: 'application/json'
  });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
};

const makeSafeFilename = (name) => {
  const cleaned = String(name || '')
    .replace(/[^\w.-]+/g, '_')
    .replace(/^_+|_+$/g, '');
  return cleaned || 'gbdraw_session';
};

const buildSessionFilename = (title) => {
  const base = String(title || '').trim();
  if (!base) return 'gbdraw_session.json';
  const safe = makeSafeFilename(base);
  return `${safe}.gbdraw-session.json`;
};

const normalizeCircularPlotTitlePosition = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['none', 'top', 'bottom'].includes(normalized) ? normalized : 'none';
};

const normalizeLinearPlotTitlePosition = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['center', 'top', 'bottom'].includes(normalized) ? normalized : 'bottom';
};

const normalizeFeatureShape = (value) => (String(value || '').trim().toLowerCase() === 'arrow' ? 'arrow' : 'rectangle');

const normalizeFeatureShapes = (featureShapes) => {
  const normalized = {};
  if (!featureShapes || typeof featureShapes !== 'object' || Array.isArray(featureShapes)) {
    return normalized;
  }
  Object.entries(featureShapes).forEach(([featureTypeRaw, shape]) => {
    const featureType = String(featureTypeRaw || '').trim();
    if (!featureType) return;
    normalized[featureType] = normalizeFeatureShape(shape);
  });
  return normalized;
};

let lastSessionFilename = null;

const buildConfigData = () => ({
  form: state.form,
  adv: state.adv,
  losat: state.losat,
  colors: state.currentColors.value,
  palette: state.selectedPalette.value,
  rules: state.manualSpecificRules,
  filterMode: state.filterMode.value,
  whitelist: state.manualWhitelist,
  blacklistText: state.manualBlacklist.value,
  blastSource: state.blastSource.value,
  losatProgram: state.losatProgram.value
});

const shouldSuppressCircularMultiRecordDefaults = (incomingForm) => {
  if (state.mode.value !== 'circular') return false;
  if (!incomingForm || typeof incomingForm !== 'object' || Array.isArray(incomingForm)) return false;
  if (!Object.prototype.hasOwnProperty.call(incomingForm, 'multi_record_canvas')) return false;
  return state.form.multi_record_canvas === false && incomingForm.multi_record_canvas === true;
};

const syncActiveModePlotTitlePosition = () => {
  if (state.mode.value === 'linear') {
    state.linearPlotTitlePosition.value = normalizeLinearPlotTitlePosition(state.adv.plot_title_position);
    state.adv.plot_title_position = state.linearPlotTitlePosition.value;
    return;
  }

  state.circularPlotTitlePosition.value = normalizeCircularPlotTitlePosition(state.adv.plot_title_position);
  state.adv.plot_title_position = state.circularPlotTitlePosition.value;
};

const restoreSessionPlotTitlePositions = (ui = {}) => {
  const activeMode = state.mode.value === 'linear' ? 'linear' : 'circular';
  const activePosition =
    activeMode === 'linear'
      ? normalizeLinearPlotTitlePosition(state.adv.plot_title_position)
      : normalizeCircularPlotTitlePosition(state.adv.plot_title_position);
  const hasCircularPlotTitlePosition =
    typeof ui.circularPlotTitlePosition === 'string' && ui.circularPlotTitlePosition.trim() !== '';
  const hasLinearPlotTitlePosition =
    typeof ui.linearPlotTitlePosition === 'string' && ui.linearPlotTitlePosition.trim() !== '';

  state.circularPlotTitlePosition.value = hasCircularPlotTitlePosition
    ? normalizeCircularPlotTitlePosition(ui.circularPlotTitlePosition)
    : activeMode === 'circular'
      ? activePosition
      : 'none';
  state.linearPlotTitlePosition.value = hasLinearPlotTitlePosition
    ? normalizeLinearPlotTitlePosition(ui.linearPlotTitlePosition)
    : activeMode === 'linear'
      ? activePosition
      : 'bottom';
  state.adv.plot_title_position =
    activeMode === 'linear'
      ? state.linearPlotTitlePosition.value
      : state.circularPlotTitlePosition.value;
};

const applyConfigData = (data) => {
  if (data.form) safeDeepMerge(state.form, data.form);
  if (data.adv) safeDeepMerge(state.adv, data.adv);
  if (state.adv.label_placement === 'on_feature') {
    state.adv.label_placement = 'above_feature';
  }
  const rawTrackAxisGap = state.adv.track_axis_gap;
  if (
    rawTrackAxisGap === null ||
    rawTrackAxisGap === undefined ||
    rawTrackAxisGap === '' ||
    String(rawTrackAxisGap).trim().toLowerCase() === 'auto'
  ) {
    state.adv.track_axis_gap = null;
  } else {
    const numericTrackAxisGap = Number(rawTrackAxisGap);
    state.adv.track_axis_gap = Number.isFinite(numericTrackAxisGap) && numericTrackAxisGap >= 0
      ? numericTrackAxisGap
      : null;
  }
  if (state.form.linear_track_layout === 'spreadout') {
    state.form.linear_track_layout = 'above';
  } else if (state.form.linear_track_layout === 'tuckin') {
    state.form.linear_track_layout = 'below';
  } else if (!['above', 'middle', 'below'].includes(state.form.linear_track_layout)) {
    state.form.linear_track_layout = 'middle';
  }
  state.form.plot_title = String(state.form.plot_title || '');
  state.adv.feature_shapes = normalizeFeatureShapes(state.adv.feature_shapes);
  const normalizedMultiRecordSizeMode = String(state.adv.multi_record_size_mode || '').trim().toLowerCase();
  if (normalizedMultiRecordSizeMode === 'sqrt') {
    state.adv.multi_record_size_mode = 'auto';
  } else {
    state.adv.multi_record_size_mode = ['auto', 'linear', 'equal'].includes(normalizedMultiRecordSizeMode)
      ? normalizedMultiRecordSizeMode
      : 'auto';
  }
  const numericMinRadiusRatio = Number(state.adv.multi_record_min_radius_ratio);
  state.adv.multi_record_min_radius_ratio =
    Number.isFinite(numericMinRadiusRatio) && numericMinRadiusRatio > 0 && numericMinRadiusRatio <= 1
      ? numericMinRadiusRatio
      : 0.55;
  const numericColumnGapRatio = Number(state.adv.multi_record_column_gap_ratio);
  state.adv.multi_record_column_gap_ratio =
    Number.isFinite(numericColumnGapRatio) && numericColumnGapRatio >= 0
      ? numericColumnGapRatio
      : 0.10;
  const numericRowGapRatio = Number(state.adv.multi_record_row_gap_ratio);
  state.adv.multi_record_row_gap_ratio =
    Number.isFinite(numericRowGapRatio) && numericRowGapRatio >= 0
      ? numericRowGapRatio
      : 0.05;
  const rawMultiRecordPositions = Array.isArray(state.adv.multi_record_positions)
    ? state.adv.multi_record_positions
    : [];
  const dedupedMultiRecordPositions = [];
  const seenMultiRecordSelectors = new Set();
  rawMultiRecordPositions.forEach((entry) => {
    if (!entry || typeof entry !== 'object' || Array.isArray(entry)) return;
    const selector = String(entry.selector ?? '').trim();
    if (!selector || seenMultiRecordSelectors.has(selector)) return;
    const rowValue = Number(entry.row);
    const normalizedRow = Number.isInteger(rowValue) && rowValue > 0 ? rowValue : 1;
    seenMultiRecordSelectors.add(selector);
    dedupedMultiRecordPositions.push({ selector, row: normalizedRow });
  });
  state.adv.multi_record_positions = dedupedMultiRecordPositions
    .map((entry, index) => ({ ...entry, __index: index }))
    .sort((left, right) => {
      if (left.row !== right.row) return left.row - right.row;
      return left.__index - right.__index;
    })
    .map(({ __index, ...entry }) => entry);
  if (state.mode.value === 'linear') {
    state.adv.plot_title_position = normalizeLinearPlotTitlePosition(state.adv.plot_title_position);
  } else {
    state.adv.plot_title_position = normalizeCircularPlotTitlePosition(state.adv.plot_title_position);
  }
  syncActiveModePlotTitlePosition();
  const rawPlotTitleFontSize = state.adv.plot_title_font_size;
  if (
    rawPlotTitleFontSize === null ||
    rawPlotTitleFontSize === undefined ||
    rawPlotTitleFontSize === ''
  ) {
    state.adv.plot_title_font_size = null;
  } else {
    const numericPlotTitleFontSize = Number(rawPlotTitleFontSize);
    state.adv.plot_title_font_size =
      Number.isFinite(numericPlotTitleFontSize) && numericPlotTitleFontSize > 0
        ? numericPlotTitleFontSize
        : null;
  }
  state.adv.keep_full_definition_with_plot_title =
    state.adv.keep_full_definition_with_plot_title === true;
  state.adv.linear_show_replicon = state.adv.linear_show_replicon === true;
  state.adv.linear_show_accession = state.adv.linear_show_accession !== false;
  state.adv.linear_show_length = state.adv.linear_show_length !== false;
  if (data.losat) safeDeepMerge(state.losat, data.losat);
  if (data.colors) {
    const normalized = {};
    Object.entries(data.colors).forEach(([key, value]) => {
      normalized[key] = resolveColorToHex(String(value || '').trim());
    });
    state.currentColors.value = normalized;
  }
  if (data.palette) state.selectedPalette.value = data.palette;

  if (data.rules && Array.isArray(data.rules)) {
    state.manualSpecificRules.length = 0;
    data.rules.forEach((r) => {
      state.manualSpecificRules.push({
        feat: String(r.feat || ''),
        qual: String(r.qual || ''),
        val: String(r.val || ''),
        color: resolveColorToHex(String(r.color || '#000000')),
        cap: String(r.cap || ''),
        fromFile: !!r.fromFile
      });
    });
  }
  if (data.filterMode) state.filterMode.value = data.filterMode;
  if (data.whitelist && Array.isArray(data.whitelist)) {
    state.manualWhitelist.length = 0;
    data.whitelist.forEach((w) => {
      state.manualWhitelist.push({
        feat: String(w.feat || ''),
        qual: String(w.qual || ''),
        key: String(w.key || '')
      });
    });
  }
  if (data.blacklistText !== undefined) state.manualBlacklist.value = String(data.blacklistText || '');
  if (data.blastSource) state.blastSource.value = String(data.blastSource);
  if (data.losatProgram) state.losatProgram.value = String(data.losatProgram);
};

const bufferToBase64 = (buffer) => {
  let binary = '';
  const bytes = new Uint8Array(buffer);
  const chunkSize = 0x8000;
  for (let i = 0; i < bytes.length; i += chunkSize) {
    binary += String.fromCharCode(...bytes.subarray(i, i + chunkSize));
  }
  return btoa(binary);
};

const base64ToUint8 = (base64) => {
  const binary = atob(base64);
  const len = binary.length;
  const bytes = new Uint8Array(len);
  for (let i = 0; i < len; i += 1) {
    bytes[i] = binary.charCodeAt(i);
  }
  return bytes;
};

const serializeFile = async (file) => {
  if (!file) return null;
  const buffer = await file.arrayBuffer();
  return {
    name: file.name || 'file',
    type: file.type || '',
    size: file.size || buffer.byteLength,
    lastModified: file.lastModified || Date.now(),
    data: bufferToBase64(buffer)
  };
};

const deserializeFile = (entry) => {
  if (!entry || !entry.data) return null;
  const bytes = base64ToUint8(entry.data);
  return new File([bytes], entry.name || 'file', {
    type: entry.type || 'application/octet-stream',
    lastModified: entry.lastModified || Date.now()
  });
};

const serializeResults = () => {
  const currentSvg = (() => {
    if (!state.svgContainer.value) return null;
    const svg = state.svgContainer.value.querySelector('svg');
    if (!svg) return null;
    const serializer = new XMLSerializer();
    return serializer.serializeToString(svg);
  })();

  return state.results.value.map((res, idx) => ({
    name: res.name || `Result ${idx + 1}`,
    content: idx === state.selectedResultIndex.value && currentSvg ? currentSvg : res.content
  }));
};

const serializeLosatCache = () => {
  const cacheMap = state.losatCache?.value;
  if (!cacheMap || cacheMap.size === 0) return [];
  const info = Array.isArray(state.losatCacheInfo.value) ? state.losatCacheInfo.value : [];
  const entries = [];
  const seen = new Set();

  info.forEach((entry, idx) => {
    if (!entry || !entry.key) return;
    const cached = cacheMap.get(entry.key);
    if (!cached || typeof cached.text !== 'string') return;
    entries.push({
      key: entry.key,
      filename: entry.filename || `losat_pair_${idx + 1}.tsv`,
      text: cached.text
    });
    seen.add(entry.key);
  });

  cacheMap.forEach((value, key) => {
    if (seen.has(key)) return;
    if (!value || typeof value.text !== 'string') return;
    entries.push({ key, filename: '', text: value.text });
  });

  return entries;
};

const applyLosatCache = (entries) => {
  const map = new Map();
  const info = [];

  if (Array.isArray(entries)) {
    entries.forEach((entry, idx) => {
      if (!entry || !entry.key || typeof entry.text !== 'string') return;
      map.set(entry.key, { text: entry.text });
      info.push({
        key: entry.key,
        filename: entry.filename || `losat_pair_${idx + 1}.tsv`
      });
    });
  }

  state.losatCache.value = map;
  state.losatCacheInfo.value = info;
};

const serializeFiles = async () => {
  const normalizedLinearSeqs = normalizeLinearSeqList(state.linearSeqs);
  const linearSeqs = await Promise.all(
    normalizedLinearSeqs.map(async (seq) => ({
      uid: seq.uid,
      gb: await serializeFile(seq.gb),
      gff: await serializeFile(seq.gff),
      fasta: await serializeFile(seq.fasta),
      blast: await serializeFile(seq.blast),
      losat_gencode: seq.losat_gencode ?? 1,
      losat_filename: seq.losat_filename ?? '',
      definition: seq.definition ?? '',
      region_record_id: seq.region_record_id ?? '',
      region_start: seq.region_start ?? null,
      region_end: seq.region_end ?? null,
      region_reverse: !!seq.region_reverse
    }))
  );

  return {
    c_gb: await serializeFile(state.files.c_gb),
    c_gff: await serializeFile(state.files.c_gff),
    c_fasta: await serializeFile(state.files.c_fasta),
    d_color: await serializeFile(state.files.d_color),
    t_color: await serializeFile(state.files.t_color),
    blacklist: await serializeFile(state.files.blacklist),
    whitelist: await serializeFile(state.files.whitelist),
    qualifier_priority: await serializeFile(state.files.qualifier_priority),
    linearSeqs
  };
};

const applyFiles = (filesData) => {
  state.files.c_gb = null;
  state.files.c_gff = null;
  state.files.c_fasta = null;
  state.files.d_color = null;
  state.files.t_color = null;
  state.files.blacklist = null;
  state.files.whitelist = null;
  state.files.qualifier_priority = null;
  state.linearReorderNotice.value = '';

  if (!filesData) {
    state.linearSeqs.splice(0, state.linearSeqs.length, ...normalizeLinearSeqList([]));
    return { collapsedLinearSeqs: false };
  }

  state.files.c_gb = deserializeFile(filesData.c_gb);
  state.files.c_gff = deserializeFile(filesData.c_gff);
  state.files.c_fasta = deserializeFile(filesData.c_fasta);
  state.files.d_color = deserializeFile(filesData.d_color);
  state.files.t_color = deserializeFile(filesData.t_color);
  state.files.blacklist = deserializeFile(filesData.blacklist);
  state.files.whitelist = deserializeFile(filesData.whitelist);
  state.files.qualifier_priority = deserializeFile(filesData.qualifier_priority);

  if (Array.isArray(filesData.linearSeqs)) {
    const loadedLinearSeqs = filesData.linearSeqs.map((seq) => ({
      uid: seq.uid,
      gb: deserializeFile(seq.gb),
      gff: deserializeFile(seq.gff),
      fasta: deserializeFile(seq.fasta),
      blast: deserializeFile(seq.blast),
      losat_gencode: seq.losat_gencode ?? 1,
      losat_filename: seq.losat_filename ?? '',
      definition: seq.definition ?? '',
      region_record_id: seq.region_record_id ?? '',
      region_start: seq.region_start ?? null,
      region_end: seq.region_end ?? null,
      region_reverse: !!seq.region_reverse
    }));
    const normalized = normalizeLinearSeqList(loadedLinearSeqs);
    const collapsed = collapseEmptyLinearSeqList(loadedLinearSeqs);
    const collapsedLinearSeqs = collapsed.length !== normalized.length;
    state.linearSeqs.splice(0, state.linearSeqs.length, ...collapsed);
    return { collapsedLinearSeqs };
  }

  state.linearSeqs.splice(0, state.linearSeqs.length, ...normalizeLinearSeqList([]));
  return { collapsedLinearSeqs: false };
};

export const exportConfig = () => {
  const configData = buildConfigData();
  downloadJson(configData, 'gbdraw_config.json');
};

export const exportSession = async (titleOverride = null) => {
  const losatEntries = serializeLosatCache();
  const losatBytes = losatEntries.reduce((sum, entry) => sum + (entry.text ? entry.text.length : 0), 0);
  const currentLegend = state.form.legend;
  const isLinear = state.mode.value === 'linear';
  const savedCircularLegend = isLinear ? state.circularLegendPosition.value : currentLegend;
  const savedLinearLegend = isLinear ? currentLegend : state.linearLegendPosition.value;
  const currentPlotTitlePosition = state.adv.plot_title_position;
  const savedCircularPlotTitlePosition = isLinear
    ? normalizeCircularPlotTitlePosition(state.circularPlotTitlePosition.value)
    : normalizeCircularPlotTitlePosition(currentPlotTitlePosition);
  const savedLinearPlotTitlePosition = isLinear
    ? normalizeLinearPlotTitlePosition(currentPlotTitlePosition)
    : normalizeLinearPlotTitlePosition(state.linearPlotTitlePosition.value);
  const resolvedTitle =
    typeof titleOverride === 'string'
      ? titleOverride.trim()
      : typeof state.sessionTitle?.value === 'string'
        ? state.sessionTitle.value.trim()
        : '';
  const sessionFilename = buildSessionFilename(resolvedTitle);
  const totalBytes =
    (state.files.c_gb?.size || 0) +
    (state.files.c_gff?.size || 0) +
    (state.files.c_fasta?.size || 0) +
    (state.files.d_color?.size || 0) +
    (state.files.t_color?.size || 0) +
    (state.files.blacklist?.size || 0) +
    (state.files.whitelist?.size || 0) +
    (state.files.qualifier_priority?.size || 0) +
    state.linearSeqs.reduce((sum, seq) => {
      return (
        sum +
        (seq.gb?.size || 0) +
        (seq.gff?.size || 0) +
        (seq.fasta?.size || 0) +
        (seq.blast?.size || 0)
      );
    }, 0) +
    losatBytes;

  if (totalBytes > 50 * 1024 * 1024) {
    const proceed = confirm(
      `Session file will include ${(totalBytes / (1024 * 1024)).toFixed(
        1
      )} MB of input data. Continue?`
    );
    if (!proceed) return;
  }

  const sessionData = {
    format: 'gbdraw-session',
    version: SESSION_VERSION,
    createdAt: new Date().toISOString(),
    title: resolvedTitle || undefined,
    config: buildConfigData(),
    ui: {
      mode: state.mode.value,
      zoom: state.zoom.value,
      canvasPan: { x: state.canvasPan.x, y: state.canvasPan.y },
      canvasPadding: { ...state.canvasPadding },
      selectedResultIndex: state.selectedResultIndex.value,
      generatedLegendPosition: state.generatedLegendPosition.value,
      legend: currentLegend,
      circularLegendPosition: savedCircularLegend,
      linearLegendPosition: savedLinearLegend,
      circularPlotTitlePosition: savedCircularPlotTitlePosition,
      linearPlotTitlePosition: savedLinearPlotTitlePosition,
      featurePanelTab: state.featurePanelTab.value,
      cInputType: state.cInputType.value,
      lInputType: state.lInputType.value,
      downloadDpi: state.downloadDpi.value,
      autoLabelReflow: Boolean(state.autoLabelReflowEnabled.value)
    },
    files: await serializeFiles(),
    results: serializeResults(),
    features: {
      extractedFeatures: state.extractedFeatures.value,
      featureRecordIds: state.featureRecordIds.value,
      selectedFeatureRecordIdx: state.selectedFeatureRecordIdx.value,
      featureColorOverrides: JSON.parse(JSON.stringify(state.featureColorOverrides)),
      featureVisibilityOverrides: JSON.parse(JSON.stringify(state.featureVisibilityOverrides)),
      labelTextFeatureOverrides: JSON.parse(JSON.stringify(state.labelTextFeatureOverrides)),
      labelTextBulkOverrides: JSON.parse(JSON.stringify(state.labelTextBulkOverrides)),
      labelTextFeatureOverrideSources: JSON.parse(JSON.stringify(state.labelTextFeatureOverrideSources)),
      labelVisibilityOverrides: JSON.parse(JSON.stringify(state.labelVisibilityOverrides)),
      labelOverrideContextKey: String(state.labelOverrideContextKey.value || '')
    },
    losatCache: {
      entries: losatEntries
    }
  };

  if (lastSessionFilename && lastSessionFilename === sessionFilename) {
    const proceed = confirm(`Download "${sessionFilename}" again? Your browser may overwrite or rename the file.`);
    if (!proceed) return;
  }
  lastSessionFilename = sessionFilename;
  downloadJson(sessionData, sessionFilename);
};

export const importConfig = async (e) => {
  const file = e.target.files[0];
  if (!file) return;

  if (file.size > 10 * 1024 * 1024) {
    alert('Config file is too large.');
    return;
  }

  try {
    const text = await file.text();
    const data = JSON.parse(text, (key, value) => {
      if (key === '__proto__' || key === 'constructor' || key === 'prototype') {
        return undefined;
      }
      return value;
    });
    state.suppressCircularMultiRecordDefaults.value = shouldSuppressCircularMultiRecordDefaults(data.form);
    applyConfigData(data);
    alert('Configuration loaded successfully!');
  } catch (err) {
    console.error(err);
    state.suppressCircularMultiRecordDefaults.value = false;
    alert('Failed to load config: Invalid JSON structure.');
  } finally {
    e.target.value = '';
  }
};

export const importSession = async (e) => {
  const file = e.target.files[0];
  if (!file) return;

  if (file.size > 200 * 1024 * 1024) {
    alert('Session file is too large.');
    return;
  }

  try {
    const text = await file.text();
    const data = JSON.parse(text, (key, value) => {
      if (key === '__proto__' || key === 'constructor' || key === 'prototype') {
        return undefined;
      }
      return value;
    });

    if (!data || data.format !== 'gbdraw-session') {
      alert('Invalid session file.');
      return;
    }

    const ui = data.ui || {};
    state.sessionTitle.value = typeof data.title === 'string' ? data.title : '';
    if (ui.mode) state.mode.value = ui.mode;
    if (ui.cInputType) state.cInputType.value = ui.cInputType;
    if (ui.lInputType) state.lInputType.value = ui.lInputType;
    if (ui.downloadDpi) state.downloadDpi.value = ui.downloadDpi;
    state.autoLabelReflowEnabled.value = Boolean(ui.autoLabelReflow);
    state.labelOverrideBuildWarning.value = '';
    if (ui.featurePanelTab === 'labels' || ui.featurePanelTab === 'colors') {
      state.featurePanelTab.value = ui.featurePanelTab;
    } else {
      state.featurePanelTab.value = 'colors';
    }
    state.generatedMode.value = ui.mode === 'linear' ? 'linear' : 'circular';
    if (ui.circularLegendPosition) state.circularLegendPosition.value = ui.circularLegendPosition;
    if (ui.linearLegendPosition) state.linearLegendPosition.value = ui.linearLegendPosition;
    if (ui.generatedLegendPosition) state.generatedLegendPosition.value = ui.generatedLegendPosition;

    if (data.config) {
      state.suppressCircularMultiRecordDefaults.value = shouldSuppressCircularMultiRecordDefaults(data.config.form);
      applyConfigData(data.config);
    }
    restoreSessionPlotTitlePositions(ui);

    const { collapsedLinearSeqs } = applyFiles(data.files);
    applyLosatCache(data.losatCache?.entries);
    if (collapsedLinearSeqs) {
      state.losatCacheInfo.value = [];
    }

    state.skipCaptureBaseConfig.value = false;
    state.skipPositionReapply.value = false;

    if (Array.isArray(data.results)) {
      state.results.value = data.results.map((res, idx) => ({
        name: res.name || `Result ${idx + 1}`,
        content: res.content || ''
      }));
    } else {
      state.results.value = [];
    }

    const features = data.features || {};
    if (Array.isArray(features.extractedFeatures)) {
      state.extractedFeatures.value = features.extractedFeatures;
    } else {
      state.extractedFeatures.value = [];
    }
    if (Array.isArray(features.featureRecordIds)) {
      state.featureRecordIds.value = features.featureRecordIds;
    } else {
      state.featureRecordIds.value = [];
    }
    if (Number.isInteger(features.selectedFeatureRecordIdx)) {
      state.selectedFeatureRecordIdx.value = features.selectedFeatureRecordIdx;
    } else {
      state.selectedFeatureRecordIdx.value = 0;
    }
    if (features.featureColorOverrides && typeof features.featureColorOverrides === 'object') {
      Object.keys(state.featureColorOverrides).forEach((k) => delete state.featureColorOverrides[k]);
      Object.entries(features.featureColorOverrides).forEach(([key, value]) => {
        state.featureColorOverrides[key] = value;
      });
    } else {
      Object.keys(state.featureColorOverrides).forEach((k) => delete state.featureColorOverrides[k]);
    }
    if (features.featureVisibilityOverrides && typeof features.featureVisibilityOverrides === 'object') {
      Object.keys(state.featureVisibilityOverrides).forEach((k) => delete state.featureVisibilityOverrides[k]);
      Object.entries(features.featureVisibilityOverrides).forEach(([key, value]) => {
        const mode = String(value || '').trim().toLowerCase();
        if (mode !== 'on' && mode !== 'off') return;
        state.featureVisibilityOverrides[String(key || '')] = mode;
      });
    } else {
      Object.keys(state.featureVisibilityOverrides).forEach((k) => delete state.featureVisibilityOverrides[k]);
    }

    if (features.labelTextFeatureOverrides && typeof features.labelTextFeatureOverrides === 'object') {
      Object.keys(state.labelTextFeatureOverrides).forEach((k) => delete state.labelTextFeatureOverrides[k]);
      Object.entries(features.labelTextFeatureOverrides).forEach(([key, value]) => {
        state.labelTextFeatureOverrides[key] = String(value ?? '');
      });
    } else {
      Object.keys(state.labelTextFeatureOverrides).forEach((k) => delete state.labelTextFeatureOverrides[k]);
    }

    if (features.labelTextBulkOverrides && typeof features.labelTextBulkOverrides === 'object') {
      Object.keys(state.labelTextBulkOverrides).forEach((k) => delete state.labelTextBulkOverrides[k]);
      Object.entries(features.labelTextBulkOverrides).forEach(([key, value]) => {
        state.labelTextBulkOverrides[key] = String(value ?? '');
      });
    } else {
      Object.keys(state.labelTextBulkOverrides).forEach((k) => delete state.labelTextBulkOverrides[k]);
    }
    if (features.labelTextFeatureOverrideSources && typeof features.labelTextFeatureOverrideSources === 'object') {
      Object.keys(state.labelTextFeatureOverrideSources).forEach((k) => delete state.labelTextFeatureOverrideSources[k]);
      Object.entries(features.labelTextFeatureOverrideSources).forEach(([key, value]) => {
        state.labelTextFeatureOverrideSources[key] = String(value ?? '');
      });
    } else {
      Object.keys(state.labelTextFeatureOverrideSources).forEach((k) => delete state.labelTextFeatureOverrideSources[k]);
    }
    if (features.labelVisibilityOverrides && typeof features.labelVisibilityOverrides === 'object') {
      Object.keys(state.labelVisibilityOverrides).forEach((k) => delete state.labelVisibilityOverrides[k]);
      Object.entries(features.labelVisibilityOverrides).forEach(([key, value]) => {
        const normalized = String(value || '').trim().toLowerCase();
        if (normalized !== 'on' && normalized !== 'off') return;
        state.labelVisibilityOverrides[String(key || '')] = normalized;
      });
    } else {
      Object.keys(state.labelVisibilityOverrides).forEach((k) => delete state.labelVisibilityOverrides[k]);
    }
    state.labelOverrideContextKey.value = String(features.labelOverrideContextKey || '');

    const resultCount = state.results.value.length;
    if (resultCount > 0) {
      const desiredIndex =
        Number.isInteger(ui.selectedResultIndex) && ui.selectedResultIndex >= 0
          ? ui.selectedResultIndex
          : 0;
      state.selectedResultIndex.value = Math.min(desiredIndex, resultCount - 1);
    } else {
      state.selectedResultIndex.value = 0;
    }

    if (ui.canvasPadding) {
      state.canvasPadding.top = ui.canvasPadding.top || 0;
      state.canvasPadding.right = ui.canvasPadding.right || 0;
      state.canvasPadding.bottom = ui.canvasPadding.bottom || 0;
      state.canvasPadding.left = ui.canvasPadding.left || 0;
    }
    if (ui.canvasPan) {
      state.canvasPan.x = ui.canvasPan.x || 0;
      state.canvasPan.y = ui.canvasPan.y || 0;
    }
    if (typeof ui.zoom === 'number') {
      state.zoom.value = ui.zoom;
    }

    if (ui.legend) {
      state.form.legend = ui.legend;
    } else if (state.mode.value === 'linear' && ui.linearLegendPosition) {
      state.form.legend = ui.linearLegendPosition;
    } else if (state.mode.value === 'circular' && ui.circularLegendPosition) {
      state.form.legend = ui.circularLegendPosition;
    }

    if (ui.circularLegendPosition && state.mode.value !== 'circular') {
      state.circularLegendPosition.value = ui.circularLegendPosition;
    } else if (state.mode.value === 'circular') {
      state.circularLegendPosition.value = state.form.legend;
    }
    if (ui.linearLegendPosition && state.mode.value !== 'linear') {
      state.linearLegendPosition.value = ui.linearLegendPosition;
    } else if (state.mode.value === 'linear') {
      state.linearLegendPosition.value = state.form.legend;
    }

    alert('Session loaded successfully!');
  } catch (err) {
    console.error(err);
    state.suppressCircularMultiRecordDefaults.value = false;
    alert('Failed to load session: Invalid JSON structure.');
  } finally {
    e.target.value = '';
  }
};
