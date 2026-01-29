import { state } from '../state.js';
import { resolveColorToHex } from '../app/color-utils.js';

const SESSION_VERSION = 1;

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

const applyConfigData = (data) => {
  if (data.form) safeDeepMerge(state.form, data.form);
  if (data.adv) safeDeepMerge(state.adv, data.adv);
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

const createEmptyLinearSeq = () => ({
  gb: null,
  gff: null,
  fasta: null,
  blast: null,
  losat_gencode: 1,
  losat_filename: '',
  definition: '',
  region_record_id: '',
  region_start: null,
  region_end: null,
  region_reverse: false
});

const serializeFiles = async () => {
  const linearSeqs = await Promise.all(
    state.linearSeqs.map(async (seq) => ({
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

  if (!filesData) {
    state.linearSeqs.splice(0, state.linearSeqs.length, createEmptyLinearSeq());
    return;
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
    const normalized = filesData.linearSeqs.map((seq) => ({
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
    if (normalized.length > 0) {
      state.linearSeqs.splice(0, state.linearSeqs.length, ...normalized);
    } else {
      state.linearSeqs.splice(0, state.linearSeqs.length, createEmptyLinearSeq());
    }
  } else {
    state.linearSeqs.splice(0, state.linearSeqs.length, createEmptyLinearSeq());
  }
};

export const exportConfig = () => {
  const configData = buildConfigData();
  downloadJson(configData, 'gbdraw_config.json');
};

export const exportSession = async () => {
  const losatEntries = serializeLosatCache();
  const losatBytes = losatEntries.reduce((sum, entry) => sum + (entry.text ? entry.text.length : 0), 0);
  const currentLegend = state.form.legend;
  const isLinear = state.mode.value === 'linear';
  const savedCircularLegend = isLinear ? state.circularLegendPosition.value : currentLegend;
  const savedLinearLegend = isLinear ? currentLegend : state.linearLegendPosition.value;
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
      cInputType: state.cInputType.value,
      lInputType: state.lInputType.value,
      downloadDpi: state.downloadDpi.value
    },
    files: await serializeFiles(),
    results: serializeResults(),
    features: {
      extractedFeatures: state.extractedFeatures.value,
      featureRecordIds: state.featureRecordIds.value,
      selectedFeatureRecordIdx: state.selectedFeatureRecordIdx.value,
      featureColorOverrides: JSON.parse(JSON.stringify(state.featureColorOverrides))
    },
    losatCache: {
      entries: losatEntries
    }
  };

  downloadJson(sessionData, 'gbdraw_session.json');
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
    applyConfigData(data);
    alert('Configuration loaded successfully!');
  } catch (err) {
    console.error(err);
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
    if (ui.mode) state.mode.value = ui.mode;
    if (ui.cInputType) state.cInputType.value = ui.cInputType;
    if (ui.lInputType) state.lInputType.value = ui.lInputType;
    if (ui.downloadDpi) state.downloadDpi.value = ui.downloadDpi;
    if (ui.circularLegendPosition) state.circularLegendPosition.value = ui.circularLegendPosition;
    if (ui.linearLegendPosition) state.linearLegendPosition.value = ui.linearLegendPosition;
    if (ui.generatedLegendPosition) state.generatedLegendPosition.value = ui.generatedLegendPosition;

    if (data.config) {
      applyConfigData(data.config);
    }

    applyFiles(data.files);
    applyLosatCache(data.losatCache?.entries);

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
    alert('Failed to load session: Invalid JSON structure.');
  } finally {
    e.target.value = '';
  }
};
