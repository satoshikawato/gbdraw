import { prepareLosatRuntime, runLosatPairsParallel } from '../services/losat.js';
import {
  cancelDiagramGeneration,
  DiagramGenerationCanceledError,
  isDiagramGenerationCanceled,
  runDiagramGeneration
} from '../services/diagram-generation.js';
import { buildLabelOverrideTsv } from './feature-editor/label-override-table.js';
import {
  applyCircularTrackOrderPlacements,
  buildCircularTrackSlotSpec,
  clampCircularTrackAxisIndex,
  hasEnabledCircularTrackRenderer,
  inferLegacyAxisIndexFromFeature
} from './circular-track-slots.js';
import {
  normalizeFileList,
  orderedConservationSources
} from './conservation-series.js';

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
const fastaExtractionCache = new WeakMap();
const FASTA_EXTRACTION_CACHE_LIMIT = 12;
const proteinExtractionCache = new WeakMap();
const PROTEIN_EXTRACTION_CACHE_LIMIT = 16;
const featureExtractionCache = new WeakMap();
const FEATURE_EXTRACTION_CACHE_LIMIT = 16;
const LOSAT_CACHE_SCHEMA = 2;

const isRawLosatCacheEntry = (entry) =>
  Boolean(entry) &&
  entry.schema === LOSAT_CACHE_SCHEMA &&
  entry.kind === 'raw-losat' &&
  typeof entry.text === 'string';

const makeSafeFilename = (name) => {
  const cleaned = String(name || '').replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '');
  return cleaned || 'losat';
};
const normalizeRecordSelectorText = (value) => {
  const normalized = String(value ?? '').trim();
  if (!normalized || ['none', 'null', 'jsnull', 'undefined', 'jsundefined', '-'].includes(normalized.toLowerCase())) {
    return '';
  }
  return normalized;
};
const parseRegionText = (value) => {
  const text = String(value || '').trim();
  if (!text) return null;
  const match = text.match(/^(\d+)(?:\.\.|-)(\d+)(?::(rc|rev|reverse|minus|-))?$/i);
  if (!match) throw new Error(`Invalid region spec for LOSAT FASTA extraction: ${text}`);
  let start = Number(match[1]);
  let end = Number(match[2]);
  let reverse = Boolean(match[3]);
  if (!Number.isInteger(start) || !Number.isInteger(end) || start < 1 || end < 1) {
    throw new Error(`Invalid region coordinates for LOSAT FASTA extraction: ${text}`);
  }
  if (start > end) {
    [start, end] = [end, start];
    reverse = true;
  }
  return { start, end, reverse };
};
const reverseComplementSequence = (sequence) => {
  const complements = {
    A: 'T',
    C: 'G',
    G: 'C',
    T: 'A',
    U: 'A',
    R: 'Y',
    Y: 'R',
    S: 'S',
    W: 'W',
    K: 'M',
    M: 'K',
    B: 'V',
    D: 'H',
    H: 'D',
    V: 'B',
    N: 'N'
  };
  let out = '';
  const upper = String(sequence || '').toUpperCase();
  for (let i = upper.length - 1; i >= 0; i -= 1) {
    out += complements[upper[i]] || 'N';
  }
  return out;
};
const wrapFastaSequence = (sequence) => {
  const lines = [];
  for (let i = 0; i < sequence.length; i += 60) lines.push(sequence.slice(i, i + 60));
  return lines.join('\n');
};
const buildFastaText = (record) => `>${record.id}\n${wrapFastaSequence(record.sequence)}\n`;
const getFastaSequenceLength = (fasta) =>
  String(fasta || '')
    .split(/\r?\n/)
    .filter((line) => line && !line.startsWith('>'))
    .join('')
    .replace(/\s+/g, '').length;
const selectParsedRecord = (records, selectorRaw) => {
  if (!records.length) throw new Error('No records found');
  const selector = normalizeRecordSelectorText(selectorRaw);
  if (!selector) return records[0];
  if (selector.startsWith('#')) {
    const idx = Number(selector.slice(1).trim()) - 1;
    if (!Number.isInteger(idx) || idx < 0 || idx >= records.length) {
      throw new Error(`Record selector ${selector} is out of range (loaded ${records.length} record(s)).`);
    }
    return records[idx];
  }
  const matches = records.filter((record) => record.id === selector);
  if (matches.length === 0) throw new Error(`Record selector '${selector}' did not match any record ID.`);
  if (matches.length > 1) throw new Error(`Record selector '${selector}' matched multiple records. Use #index to disambiguate.`);
  return matches[0];
};
const parseFastaRecordsFast = (text) => {
  const records = [];
  let current = null;
  String(text || '').split(/\r?\n/).forEach((line) => {
    if (line.startsWith('>')) {
      if (current) records.push({ ...current, sequence: current.parts.join('').toUpperCase() });
      const header = line.slice(1).trim();
      current = { id: header.split(/\s+/)[0] || `record_${records.length + 1}`, parts: [] };
      return;
    }
    if (current && line.trim()) current.parts.push(line.replace(/\s+/g, ''));
  });
  if (current) records.push({ ...current, sequence: current.parts.join('').toUpperCase() });
  return records;
};
const parseGenbankRecordsFast = (text) => {
  const records = [];
  const recordChunks = String(text || '').split(/^\/\/\s*$/m);
  recordChunks.forEach((chunk) => {
    const originMatch = chunk.match(/\nORIGIN\b([\s\S]*)$/i);
    if (!originMatch) return;
    const locusMatch = chunk.match(/^LOCUS\s+(\S+)/m);
    const accessionMatch = chunk.match(/^ACCESSION\s+(\S+)/m);
    const versionMatch = chunk.match(/^VERSION\s+(\S+)/m);
    const id = versionMatch?.[1] || accessionMatch?.[1] || locusMatch?.[1] || `record_${records.length + 1}`;
    const sequence = originMatch[1].replace(/[^A-Za-z]/g, '').toUpperCase();
    if (sequence) records.push({ id, sequence });
  });
  return records;
};
const applyLosatSequenceTransforms = (record, regionSpec, reverseFlag) => {
  const region = parseRegionText(regionSpec);
  let sequence = record.sequence;
  if (String(reverseFlag).trim().toLowerCase() === '1') sequence = reverseComplementSequence(sequence);
  if (region) {
    const start = Math.max(0, region.start - 1);
    const end = Math.min(sequence.length, region.end);
    if (start >= end && (region.start !== 1 || region.end !== sequence.length)) {
      throw new Error(`Start position (${region.start}) must be less than end position (${region.end}).`);
    }
    sequence = sequence.slice(start, end);
    if (region.reverse) sequence = reverseComplementSequence(sequence);
  }
  return { id: record.id, sequence };
};
const getCachedFastaExtraction = (file, key) => {
  const byKey = fastaExtractionCache.get(file);
  return byKey?.get(key) || null;
};
const setCachedFastaExtraction = (file, key, value) => {
  let byKey = fastaExtractionCache.get(file);
  if (!byKey) {
    byKey = new Map();
    fastaExtractionCache.set(file, byKey);
  }
  if (byKey.size >= FASTA_EXTRACTION_CACHE_LIMIT) byKey.delete(byKey.keys().next().value);
  byKey.set(key, value);
};
const getFileFingerprint = (file) => {
  if (!file) return null;
  return {
    name: String(file.name || ''),
    size: Number(file.size || 0),
    lastModified: Number(file.lastModified || 0)
  };
};
const getCachedProteinExtraction = (file, key) => {
  if (!file) return null;
  const byKey = proteinExtractionCache.get(file);
  return byKey?.get(key) || null;
};
const setCachedProteinExtraction = (file, key, value) => {
  if (!file || value?.error) return;
  let byKey = proteinExtractionCache.get(file);
  if (!byKey) {
    byKey = new Map();
    proteinExtractionCache.set(file, byKey);
  }
  if (byKey.size >= PROTEIN_EXTRACTION_CACHE_LIMIT) byKey.delete(byKey.keys().next().value);
  byKey.set(key, value);
};
const cloneFeatureExtractionData = (data) => ({
  features: Array.isArray(data?.features)
    ? data.features.map((feature) => ({
        ...feature,
        qualifiers:
          feature?.qualifiers && typeof feature.qualifiers === 'object'
            ? Object.fromEntries(
                Object.entries(feature.qualifiers).map(([key, value]) => [
                  key,
                  Array.isArray(value) ? [...value] : value
                ])
              )
            : feature?.qualifiers
      }))
    : [],
  record_ids: Array.isArray(data?.record_ids) ? [...data.record_ids] : []
});
const getCachedFeatureExtraction = (file, key) => {
  if (!file) return null;
  const byKey = featureExtractionCache.get(file);
  const cached = byKey?.get(key) || null;
  return cached ? cloneFeatureExtractionData(cached) : null;
};
const setCachedFeatureExtraction = (file, key, value) => {
  if (!file || value?.error) return;
  let byKey = featureExtractionCache.get(file);
  if (!byKey) {
    byKey = new Map();
    featureExtractionCache.set(file, byKey);
  }
  if (byKey.size >= FEATURE_EXTRACTION_CACHE_LIMIT) byKey.delete(byKey.keys().next().value);
  byKey.set(key, cloneFeatureExtractionData(value));
};
const measureTiming = (entries, label, fn) => {
  const startedAt = getNow();
  const result = fn();
  entries.push({ label, ms: getNow() - startedAt });
  return result;
};
const logPostGbdrawTimings = (entries) => {
  if (!entries || entries.length === 0) return;
  console.groupCollapsed('post-gbdraw timing');
  entries.forEach(({ label, ms, details }) => {
    console.info(`${label}: ${formatDuration(ms)}${details ? ` (${details})` : ''}`);
  });
  console.groupEnd();
};
const extractLosatFastaFast = async ({ file, text, fmt, regionSpec, recordSelector, reverseFlag }) => {
  if (typeof text !== 'string' && !file?.text) {
    throw new Error('Input file is not available for browser FASTA extraction.');
  }
  const sourceText = typeof text === 'string' ? text : await file.text();
  const records = fmt === 'genbank' ? parseGenbankRecordsFast(sourceText) : parseFastaRecordsFast(sourceText);
  const selected = selectParsedRecord(records, recordSelector);
  const transformed = applyLosatSequenceTransforms(selected, regionSpec, reverseFlag);
  return {
    fasta: buildFastaText(transformed),
    recordId: transformed.id,
    canonicalLength: transformed.sequence.length
  };
};
const extractAllLosatFastaFast = async ({ file, text, fmt }) => {
  if (typeof text !== 'string' && !file?.text) {
    throw new Error('Input file is not available for browser FASTA extraction.');
  }
  const sourceText = typeof text === 'string' ? text : await file.text();
  const records = fmt === 'genbank' ? parseGenbankRecordsFast(sourceText) : parseFastaRecordsFast(sourceText);
  if (!records.length) throw new Error('No records found for circular conservation reference.');
  return {
    fasta: records.map((record) => buildFastaText(record)).join(''),
    recordIds: records.map((record) => record.id),
    canonicalLength: records.reduce((sum, record) => sum + String(record.sequence || '').length, 0)
  };
};
const escapeRegexLiteral = (value) => String(value ?? '').replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
const buildConservationSeries = (sourceFiles, circularConservation) => {
  return orderedConservationSources(sourceFiles, circularConservation).map((entry) => ({
    label: entry.label,
    color: entry.color
  }));
};
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
const normalizeBlastpMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['pairwise', 'orthogroup', 'collinear'].includes(normalized) ? normalized : 'orthogroup';
};
const normalizeCircularConservationLosatProgram = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return normalized === 'tblastx' ? 'tblastx' : 'blastn';
};
const normalizePositiveInteger = (value, fallback = 1) => {
  const parsed = Number(value);
  return Number.isInteger(parsed) && parsed > 0 ? parsed : fallback;
};
const normalizeCollinearColorMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
  if (normalized === 'identity') return 'average_identity';
  return ['average_identity', 'orientation', 'orientation_identity'].includes(normalized) ? normalized : 'orientation';
};
const normalizeCollinearAnchorMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
  const aliases = {
    raw: 'all',
    all_hits: 'all',
    top_n: 'all',
    topn: 'all',
    one2one: 'one_to_one',
    mutual_best: 'one_to_one',
    top1: 'one_to_one',
    top_1: 'one_to_one',
    reciprocal_best: 'rbh',
    strict_rbh: 'rbh'
  };
  const resolved = aliases[normalized] || normalized;
  return ['all', 'one_to_one', 'rbh'].includes(resolved) ? resolved : 'rbh';
};
const normalizeCollinearSearchScope = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
  return ['adjacent', 'all'].includes(normalized) ? normalized : 'adjacent';
};
const normalizePairwiseMatchStyle = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['ribbon', 'curve'].includes(normalized) ? normalized : 'ribbon';
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
    diagramGenerationWorkerReady,
    diagramGenerationWorkerStatus,
    diagramGenerationWorkerError,
    processing,
    processingStatus,
    generationCancelRequested,
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
    losatThreadingStatus,
    losatCache,
    circularConservation,
    orthogroups,
    featureOrthogroupIndex,
    selectedOrthogroupAlignmentFeature,
    orthogroupNameOverrides,
    orthogroupDescriptionOverrides,
    selectedOrthogroupId,
    circularRecordList,
    files,
    linearSeqs,
    generatedLegendPosition,
    generatedMode,
    generatedMultiRecordCanvas,
    generatedCircularPlotTitlePosition,
    shouldDeferCircularPreviewUpdates,
    extractedFeatures,
    featureExtractionPending,
    featureExtractionError,
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
  let featureExtractionRequestId = 0;
  let latestGenerationToken = 0;
  let activeLosatAbortController = null;

  const getLastLine = (text) => {
    const trimmed = String(text || '').trim();
    if (!trimmed) return '';
    const lines = trimmed.split(/\r?\n/);
    return lines[lines.length - 1] || '';
  };

  const normalizeSections = (sections) =>
    sections.filter((section) => section && typeof section.text === 'string' && section.text.trim() !== '');

  const extractCircularTrackSlotError = (err) => {
    const texts = [
      err?.message,
      err?.stderr,
      err?.stdout,
      err?.traceback
    ];
    for (const text of texts) {
      const lines = String(text || '').split(/\r?\n/).map((line) => line.trim()).filter(Boolean);
      for (const line of lines) {
        const cleaned = line.replace(/^(ValueError|RuntimeError|ValidationError):\s*/, '');
        if (/^Circular track slot '.+' cannot fit inside\b/.test(cleaned)) {
          return cleaned;
        }
      }
    }
    return '';
  };

  const formatPythonError = (err) => {
    if (err && typeof err === 'object') {
      const details = normalizeSections([
        { label: 'STDERR', text: err.stderr || '' },
        { label: 'STDOUT', text: err.stdout || '' },
        { label: 'Traceback', text: err.traceback || '' }
      ]);
      const circularTrackSlotError = extractCircularTrackSlotError(err);
      let summary =
        circularTrackSlotError ||
        (err.type === 'SystemExit' && err.stderr ? getLastLine(err.stderr) : '') ||
        err.message ||
        getLastLine(err.stderr || err.stdout || err.traceback) ||
        'Unknown error';
      if (!circularTrackSlotError && err.type && summary && !summary.startsWith(err.type)) {
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

  const getGenerationCancelReason = (signal) =>
    signal?.reason instanceof Error ? signal.reason : new DiagramGenerationCanceledError();

  const waitForCancelablePromise = (promise, signal) => {
    if (!signal) return promise;
    if (signal.aborted) return Promise.reject(getGenerationCancelReason(signal));
    return new Promise((resolve, reject) => {
      const cleanup = () => signal.removeEventListener('abort', handleAbort);
      const handleAbort = () => {
        cleanup();
        reject(getGenerationCancelReason(signal));
      };
      signal.addEventListener('abort', handleAbort, { once: true });
      Promise.resolve(promise).then(
        (value) => {
          cleanup();
          if (signal.aborted) {
            reject(getGenerationCancelReason(signal));
            return;
          }
          resolve(value);
        },
        (error) => {
          cleanup();
          reject(error);
        }
      );
    });
  };

  const buildFeatureExtractionCacheKey = ({ regionSpec, recordSelector, reverseFlag, selectedFeatures }) =>
    JSON.stringify({
      regionSpec: String(regionSpec || ''),
      recordSelector: normalizeRecordSelectorText(recordSelector),
      reverseFlag: reverseFlag ? '1' : '0',
      selectedFeatures: Array.isArray(selectedFeatures) && selectedFeatures.length ? selectedFeatures : 'all'
    });

  const readFeatureExtractionData = ({
    pyodide,
    path,
    file,
    regionSpec,
    recordSelector,
    reverseFlag,
    selectedFeatures = null,
    timingEntries,
    timingLabel
  }) => {
    const cacheKey = buildFeatureExtractionCacheKey({ regionSpec, recordSelector, reverseFlag, selectedFeatures });
    const cached = getCachedFeatureExtraction(file, cacheKey);
    if (cached) {
      timingEntries.push({ label: timingLabel, ms: 0, details: 'cache hit' });
      return cached;
    }

    const selectedFeaturesArg =
      Array.isArray(selectedFeatures) && selectedFeatures.length ? JSON.stringify(selectedFeatures) : null;
    const featData = measureTiming(timingEntries, timingLabel, () => {
      const featJson = pyodide
        .globals
        .get('extract_features_from_genbank')(
          path,
          regionSpec || null,
          recordSelector || null,
          reverseFlag ? '1' : '0',
          selectedFeaturesArg
        );
      return JSON.parse(featJson);
    });
    setCachedFeatureExtraction(file, cacheKey, featData);
    return featData;
  };

  const buildOrthogroupIndexKey = (recordIndex, svgId) => `${Number(recordIndex)}:${String(svgId || '').trim()}`;

  const pruneOrthogroupOverrides = (groupIds, { clearAll = false } = {}) => {
    const validIds = new Set(Array.isArray(groupIds) ? groupIds.map((id) => String(id || '').trim()).filter(Boolean) : []);
    const pruneMap = (overrideMap) => {
      Object.keys(overrideMap).forEach((id) => {
        if (clearAll || !validIds.has(id)) delete overrideMap[id];
      });
    };
    pruneMap(orthogroupNameOverrides);
    pruneMap(orthogroupDescriptionOverrides);
  };

  const setOrthogroupMetadata = (orthogroupPayload) => {
    const groups = Array.isArray(orthogroupPayload) ? orthogroupPayload : [];
    const index = new Map();
    const groupIds = [];
    groups.forEach((group) => {
      const orthogroupId = String(group?.id || '').trim();
      if (orthogroupId) groupIds.push(orthogroupId);
      const members = Array.isArray(group?.members) ? group.members : [];
      const memberCount = Number(group?.member_count || members.length || 0);
      const recordCoverage = Number(group?.record_coverage_count || new Set(
        members.map((member) => Number(member?.recordIndex)).filter((recordIndex) => Number.isInteger(recordIndex))
      ).size || 0);
      members.forEach((member) => {
        const featureSvgId = String(member?.featureSvgId || '').trim();
        const recordIndex = Number(member?.recordIndex);
        if (!featureSvgId || !Number.isInteger(recordIndex)) return;
        const entry = {
          orthogroupId,
          orthogroupMemberCount: memberCount,
          orthogroupRecordCoverage: recordCoverage,
          proteinId: String(member?.proteinId || '').trim(),
          sourceProteinId: String(member?.sourceProteinId || '').trim(),
          orthogroupRepresentative: Boolean(member?.representative),
          orthogroupMember: member
        };
        index.set(buildOrthogroupIndexKey(recordIndex, featureSvgId), entry);
        if (!index.has(featureSvgId)) index.set(featureSvgId, entry);
      });
    });
    orthogroups.value = groups;
    featureOrthogroupIndex.value = index;
    pruneOrthogroupOverrides(groupIds);
    if (!selectedOrthogroupId.value || !groupIds.includes(String(selectedOrthogroupId.value || '').trim())) {
      selectedOrthogroupId.value = groupIds[0] || '';
    }
  };

  const clearOrthogroupMetadata = ({ clearSelection = false, clearOverrides = clearSelection } = {}) => {
    orthogroups.value = [];
    featureOrthogroupIndex.value = new Map();
    if (clearSelection) {
      selectedOrthogroupId.value = '';
      selectedOrthogroupAlignmentFeature.value = '';
    }
    if (clearOverrides) pruneOrthogroupOverrides([], { clearAll: true });
  };

  const enrichFeatureWithOrthogroup = (feature, recordIndex) => {
    const svgId = String(feature?.svg_id || '').trim();
    if (!svgId) return feature;
    const index = featureOrthogroupIndex.value instanceof Map ? featureOrthogroupIndex.value : new Map();
    const entry = index.get(buildOrthogroupIndexKey(recordIndex, svgId)) || index.get(svgId);
    if (!entry) return feature;
    return {
      ...feature,
      proteinId: entry.proteinId,
      sourceProteinId: entry.sourceProteinId,
      orthogroupId: entry.orthogroupId,
      orthogroupMemberCount: entry.orthogroupMemberCount,
      orthogroupRecordCoverage: entry.orthogroupRecordCoverage,
      orthogroupRepresentative: entry.orthogroupRepresentative,
      orthogroupMember: entry.orthogroupMember
    };
  };

  const isCurrentFeatureExtractionContext = (context) =>
    Boolean(context) &&
    context.requestId === featureExtractionRequestId &&
    context.mode === mode.value &&
    context.cInputType === cInputType.value &&
    context.lInputType === lInputType.value;

  const extractFeaturesForColorEditor = (context) => {
    if (!context) return;
    if (!isCurrentFeatureExtractionContext(context)) {
      if (context.requestId === featureExtractionRequestId) {
        featureExtractionPending.value = false;
      }
      return;
    }
    const pyodide = getPyodide();
    if (!pyodide) {
      featureExtractionPending.value = false;
      return;
    }

    const timingEntries = [];
    try {
      if (context.mode === 'circular' && context.cInputType === 'gb') {
        const featData = readFeatureExtractionData({
          pyodide,
          path: '/input.gb',
          file: context.circularFile,
          regionSpec: null,
          recordSelector: null,
          reverseFlag: false,
          timingEntries,
          timingLabel: 'feature extraction circular input'
        });
        if (!isCurrentFeatureExtractionContext(context)) return;
        if (featData?.error) {
          console.warn('Feature extraction failed for circular input:', featData.error);
          featureExtractionError.value = { summary: String(featData.error), details: [] };
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
      } else if (context.mode === 'linear' && context.lInputType === 'gb' && context.linearSeqs.length > 0) {
        let allFeatures = [];
        const allRecordLabels = [];

        for (let i = 0; i < context.linearSeqs.length; i++) {
          if (!isCurrentFeatureExtractionContext(context)) return;
          const seq = context.linearSeqs[i] || {};
          const regionSpec = context.regionSpecs[i]?.file || null;
          const recordSelector = context.recordSelectors[i] ?? '';
          const reverseFlag = Boolean(context.reverseFlags[i]);
          const featData = readFeatureExtractionData({
            pyodide,
            path: `/seq_${i}.gb`,
            file: seq.gb || null,
            regionSpec,
            recordSelector,
            reverseFlag,
            timingEntries,
            timingLabel: `feature extraction linear input #${i + 1}`
          });
          if (featData?.error) {
            console.warn(`Feature extraction failed for linear input #${i + 1}:`, featData.error);
          } else if (Array.isArray(featData?.features)) {
            const features = featData.features.map((feature) => ({
              ...enrichFeatureWithOrthogroup(feature, i),
              fileIdx: i,
              displayRecordId: `File ${i + 1}: ${feature.record_id}`,
              id: `file${i}_${feature.id}`
            }));
            allFeatures = allFeatures.concat(features);
            (featData.record_ids || []).forEach((rid, ridx) => {
              allRecordLabels.push({ label: `File ${i + 1}: ${rid}`, fileIdx: i, recordIdx: ridx });
            });
          } else {
            console.warn(`Feature extraction returned an unexpected payload for linear input #${i + 1}.`, featData);
          }
        }

        if (!isCurrentFeatureExtractionContext(context)) return;
        extractedFeatures.value = allFeatures;
        featureRecordIds.value = allRecordLabels.map((r) => r.label);
        selectedFeatureRecordIdx.value = 0;
        refreshFeatureOverrides(allFeatures);
        console.log(
          `Extracted ${allFeatures.length} features from ${context.linearSeqs.length} file(s) for color editor.`
        );
      }
    } catch (e) {
      if (context.requestId === featureExtractionRequestId) {
        featureExtractionError.value = formatJsError(e);
      }
      console.log('Could not extract features:', e);
    } finally {
      if (context.requestId === featureExtractionRequestId) {
        featureExtractionPending.value = false;
        logPostGbdrawTimings(timingEntries);
      }
    }
  };

  const scheduleFeatureExtraction = (context) => {
    featureExtractionPending.value = true;
    featureExtractionError.value = null;
    const run = () => extractFeaturesForColorEditor(context);
    if (typeof globalThis.requestIdleCallback === 'function') {
      globalThis.requestIdleCallback(run, { timeout: 500 });
    } else {
      setTimeout(run, 0);
    }
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

  const buildLosatSuffix = () => {
    if (losatProgram.value === 'blastn') return 'losatn';
    if (losatProgram.value === 'blastp') return 'losatp';
    return 'tlosatx';
  };

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
    return Number.isInteger(parsed) && parsed >= 1 ? parsed : undefined;
  };

  const getLosatExecutionMode = () => {
    const raw = String(losat.executionMode || 'auto').trim().toLowerCase();
    return ['auto', 'serial', 'threaded'].includes(raw) ? raw : 'auto';
  };

  const getLosatThreadsPerJob = () => {
    if (losatProgram.value !== 'blastp') return 1;
    const raw = String(losat.threadsPerJob || 'auto').trim().toLowerCase();
    if (raw === 'auto') return undefined;
    const parsed = Number(raw);
    return Number.isInteger(parsed) && parsed >= 1 ? parsed : undefined;
  };

  const getLosatTotalThreadBudget = () => {
    const raw = String(losat.totalThreadBudget || 'safe').trim().toLowerCase();
    if (raw === 'safe' || raw === 'auto') return undefined;
    if (raw === 'available') {
      return Math.max(1, Number(globalThis.navigator?.hardwareConcurrency || 4) || 4);
    }
    const parsed = Number(raw);
    if (!Number.isInteger(parsed) || parsed < 1) return undefined;
    const hardwareBudget = Math.max(1, Number(globalThis.navigator?.hardwareConcurrency || 4) || 4);
    return Math.min(parsed, hardwareBudget);
  };

  const getLinearLabelOptionSupport = () => {
    if (linearLabelSupportCache) return linearLabelSupportCache;
    const pyodide = getPyodide();
    if (!pyodide) {
      linearLabelSupportCache = {
        placement: false,
        rotation: false,
        linear_label_spacing: false,
        label_rendering: false,
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
        hide_length: false,
        orthogroup_alignment: false,
        pairwise_match_style: false,
        keep_definition_left_aligned: false
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
  "label_rendering": "--label_rendering" in _source,
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
  "orthogroup_alignment": "--align_orthogroup_feature" in _source,
  "pairwise_match_style": "--pairwise_match_style" in _source,
  "keep_definition_left_aligned": "--keep_definition_left_aligned" in _source,
})
      `);
      linearLabelSupportCache = JSON.parse(String(raw));
    } catch (_err) {
      linearLabelSupportCache = {
        placement: false,
        rotation: false,
        linear_label_spacing: false,
        label_rendering: false,
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
        hide_length: false,
        orthogroup_alignment: false,
        pairwise_match_style: false,
        keep_definition_left_aligned: false
      };
    }
    return linearLabelSupportCache;
  };

  const normalizeLabelRendering = (value) => {
    const normalized = String(value || '').trim().toLowerCase();
    return ['embedded_only', 'external_only'].includes(normalized) ? normalized : 'auto';
  };

  const normalizePositiveNumberOrNull = (value) => {
    if (value === null || value === undefined || value === '') return null;
    const numeric = Number(value);
    return Number.isFinite(numeric) && numeric > 0 ? numeric : null;
  };

  const normalizeNonNegativeNumberOrNull = (value) => {
    if (value === null || value === undefined || value === '') return null;
    const numeric = Number(value);
    return Number.isFinite(numeric) && numeric >= 0 ? numeric : null;
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
        center_reserved_radius: false,
        tick_label_font_size: false,
        circular_label_spacing: false,
        label_rendering: false,
        circular_track_slot: false,
        circular_track_axis_index: false,
        conservation_blast: false,
        conservation_reference: false,
        conservation_labels: false,
        conservation_colors: false,
        conservation_ring_width: false,
        conservation_ring_gap: false
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
  "center_reserved_radius": "--center_reserved_radius" in _source,
  "tick_label_font_size": "--tick_label_font_size" in _source,
  "circular_label_spacing": "--circular_label_spacing" in _source,
  "label_rendering": "--label_rendering" in _source,
  "circular_track_slot": "--circular_track_slot" in _source,
  "circular_track_axis_index": "--circular_track_axis_index" in _source,
  "conservation_blast": "--conservation_blast" in _source,
  "conservation_reference": "--conservation_reference" in _source,
  "conservation_labels": "--conservation_labels" in _source,
  "conservation_colors": "--conservation_colors" in _source,
  "conservation_ring_width": "--conservation_ring_width" in _source,
  "conservation_ring_gap": "--conservation_ring_gap" in _source,
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
        center_reserved_radius: false,
        tick_label_font_size: false,
        circular_label_spacing: false,
        label_rendering: false,
        circular_track_slot: false,
        circular_track_axis_index: false,
        conservation_blast: false,
        conservation_reference: false,
        conservation_labels: false,
        conservation_colors: false,
        conservation_ring_width: false,
        conservation_ring_gap: false
      };
    }
    return circularMultiRecordCanvasSupportCache;
  };

  const runCircularLayoutPreflight = (preflightArgs) => {
    const pyodide = getPyodide();
    if (!pyodide) return null;
    const runWrapper = pyodide.globals.get('run_gbdraw_wrapper');
    const pyArgs = pyodide.toPy((preflightArgs || []).map((arg) => String(arg)));
    try {
      const resultJson = runWrapper('circular', pyArgs, null);
      return JSON.parse(String(resultJson || 'null'));
    } finally {
      pyArgs.destroy?.();
      runWrapper.destroy?.();
    }
  };

  const downloadLosatPair = async (pairIndex, customName) => {
    const entry = losatCacheInfo.value?.[pairIndex];
    const cacheMap = losatCache.value;
    if (!entry || !cacheMap) return;
    const cached = cacheMap.get(entry.key);
    if (!isRawLosatCacheEntry(cached)) return;
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

  const validateDepthInputPresence = () => {
    if (!form.show_depth) return '';
    if (mode.value === 'circular') {
      return files.c_depth ? '' : 'Please upload a Depth TSV file or disable Show depth track.';
    }

    const depthCount = linearSeqs.filter((seq) => Boolean(seq.depth)).length;
    if (depthCount === 0) {
      return 'Please upload at least one Depth TSV file or disable Show depth track.';
    }
    if (depthCount !== 1 && depthCount !== linearSeqs.length) {
      return 'Upload one Depth TSV for all records, or one Depth TSV per sequence.';
    }
    return '';
  };

  const shouldSuppressPairwiseIdentityLegend = () => {
    return (
      mode.value === 'linear' &&
      blastSource.value === 'losat' &&
      losatProgram.value === 'blastp' &&
      normalizeBlastpMode(losat.blastp?.mode) === 'collinear' &&
      normalizeCollinearColorMode(losat.blastp?.collinearColorMode) === 'orientation'
    );
  };

  const removePairwiseIdentityLegendFromSvgContent = (content) => {
    if (typeof content !== 'string' || !content.includes('pairwise_legend')) return content;
    const doc = new DOMParser().parseFromString(content, 'image/svg+xml');
    const parseError = doc.querySelector('parsererror');
    if (parseError) return content;
    let removed = false;
    doc.querySelectorAll('#pairwise_legend, [id="pairwise_legend"]').forEach((legend) => {
      legend.remove();
      removed = true;
    });
    if (!removed) return content;
    return new XMLSerializer().serializeToString(doc.documentElement);
  };

  const stripPairwiseIdentityLegendsFromResults = (items) => {
    if (!Array.isArray(items)) return items;
    return items.map((item) => {
      if (!item || typeof item.content !== 'string') return item;
      return {
        ...item,
        content: removePairwiseIdentityLegendFromSvgContent(item.content)
      };
    });
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
        const recordLength = Number(entry?.record_length ?? 0);
        nextRecords.push({
          selector,
          record_id: recordId,
          record_length: Number.isFinite(recordLength) && recordLength > 0 ? recordLength : null
        });
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
    const isReflow = runMode === 'reflow';
    if (!pyodideReady.value) return { status: 'skipped' };
    if (!diagramGenerationWorkerReady.value) {
      if (!isReflow) {
        const message = diagramGenerationWorkerError.value
          ? `Diagram engine is not ready: ${diagramGenerationWorkerError.value}`
          : diagramGenerationWorkerStatus.value || 'Diagram engine is still preparing.';
        errorLog.value = formatJsError(new Error(message));
      }
      return { status: 'skipped' };
    }
    const pyodide = getPyodide();
    if (!pyodide) return { status: 'skipped' };

    const generationToken = ++latestGenerationToken;
    let keepProcessingStatus = false;
    let generationAbortController = null;
    let generationAbortSignal = null;
    const setProcessingStatus = (message) => {
      if (!isReflow) processingStatus.value = String(message || '');
    };
    const throwIfGenerationCanceled = () => {
      if (!isReflow && generationCancelRequested.value) {
        throw new DiagramGenerationCanceledError();
      }
    };
    if (isReflow && mode.value === 'circular' && shouldDeferCircularPreviewUpdates.value) {
      return { status: 'skipped' };
    }
    const depthInputError = validateDepthInputPresence();
    if (depthInputError) {
      if (isReflow) {
        labelReflowLastError.value = depthInputError;
      } else {
        errorLog.value = formatJsError(new Error(depthInputError));
      }
      return { status: 'error' };
    }
    generationAbortController =
      !isReflow && typeof AbortController === 'function' ? new AbortController() : null;
    generationAbortSignal = generationAbortController?.signal || null;
    if (!isReflow) activeLosatAbortController = generationAbortController;
    const previousSelectedResultIndex = selectedResultIndex.value;
    const cloneJsonSafe = (value, fallback) => {
      try {
        return JSON.parse(JSON.stringify(value));
      } catch (_err) {
        return fallback;
      }
    };
    const manualCancelSnapshot = isReflow
      ? null
      : {
          results: Array.isArray(results.value) ? results.value.map((entry) => ({ ...entry })) : [],
          selectedResultIndex: selectedResultIndex.value,
          zoom: zoom.value,
          pairwiseMatchFactors: { ...(pairwiseMatchFactors.value || {}) },
          addedLegendCaptions: new Set(addedLegendCaptions.value || []),
          fileLegendCaptions: new Set(fileLegendCaptions.value || []),
          featureColorOverrides: cloneJsonSafe(featureColorOverrides, {}),
          legendEntries: cloneJsonSafe(legendEntries.value || [], []),
          deletedLegendEntries: cloneJsonSafe(deletedLegendEntries.value || [], []),
          legendColorOverrides: cloneJsonSafe(legendColorOverrides, {}),
          originalLegendOrder: cloneJsonSafe(originalLegendOrder.value || [], []),
          originalLegendColors: cloneJsonSafe(originalLegendColors.value || {}, {}),
          extractedFeatures: cloneJsonSafe(extractedFeatures.value || [], []),
          editableLabels: cloneJsonSafe(editableLabels.value || [], []),
          featureExtractionPending: featureExtractionPending.value,
          featureExtractionError: featureExtractionError.value,
          featureRecordIds: cloneJsonSafe(featureRecordIds.value || [], []),
          selectedFeatureRecordIdx: selectedFeatureRecordIdx.value,
          labelOverrideBuildWarning: labelOverrideBuildWarning.value
        };
    const restoreManualCancelSnapshot = () => {
      if (!manualCancelSnapshot) return;
      results.value = manualCancelSnapshot.results;
      selectedResultIndex.value = Math.max(
        0,
        Math.min(manualCancelSnapshot.selectedResultIndex, Math.max(0, manualCancelSnapshot.results.length - 1))
      );
      zoom.value = manualCancelSnapshot.zoom;
      pairwiseMatchFactors.value = { ...manualCancelSnapshot.pairwiseMatchFactors };
      addedLegendCaptions.value = new Set(manualCancelSnapshot.addedLegendCaptions);
      fileLegendCaptions.value = new Set(manualCancelSnapshot.fileLegendCaptions);
      Object.keys(featureColorOverrides).forEach((k) => delete featureColorOverrides[k]);
      Object.assign(featureColorOverrides, cloneJsonSafe(manualCancelSnapshot.featureColorOverrides, {}));
      legendEntries.value = cloneJsonSafe(manualCancelSnapshot.legendEntries, []);
      deletedLegendEntries.value = cloneJsonSafe(manualCancelSnapshot.deletedLegendEntries, []);
      Object.keys(legendColorOverrides).forEach((k) => delete legendColorOverrides[k]);
      Object.assign(legendColorOverrides, cloneJsonSafe(manualCancelSnapshot.legendColorOverrides, {}));
      originalLegendOrder.value = cloneJsonSafe(manualCancelSnapshot.originalLegendOrder, []);
      originalLegendColors.value = cloneJsonSafe(manualCancelSnapshot.originalLegendColors, {});
      extractedFeatures.value = cloneJsonSafe(manualCancelSnapshot.extractedFeatures, []);
      editableLabels.value = cloneJsonSafe(manualCancelSnapshot.editableLabels, []);
      featureExtractionPending.value = manualCancelSnapshot.featureExtractionPending;
      featureExtractionError.value = manualCancelSnapshot.featureExtractionError;
      featureRecordIds.value = cloneJsonSafe(manualCancelSnapshot.featureRecordIds, []);
      selectedFeatureRecordIdx.value = manualCancelSnapshot.selectedFeatureRecordIdx;
      labelOverrideBuildWarning.value = manualCancelSnapshot.labelOverrideBuildWarning;
    };
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
      featureExtractionRequestId += 1;
      generationCancelRequested.value = false;
      featureExtractionPending.value = false;
      featureExtractionError.value = null;
      processing.value = true;
      processingStatus.value = 'Preparing input files...';
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
      const generationFileMap = new Map();
      const textEncoder = new TextEncoder();
      const textDecoder = new TextDecoder();
      const getPayloadName = (path, fallback = 'input') => {
        const name = String(path || '').split('/').filter(Boolean).pop();
        return name || fallback;
      };
      const toTransferableBuffer = (bytes) => {
        if (bytes instanceof ArrayBuffer) return bytes;
        if (ArrayBuffer.isView(bytes)) {
          const { buffer, byteOffset, byteLength } = bytes;
          if (byteOffset === 0 && byteLength === buffer.byteLength) return buffer;
          return buffer.slice(byteOffset, byteOffset + byteLength);
        }
        return new Uint8Array(bytes || []).buffer;
      };
      const stageGenerationBytes = (path, name, bytes) => {
        const normalizedPath = String(path || '').trim();
        if (!normalizedPath) return;
        generationFileMap.set(normalizedPath, {
          path: normalizedPath,
          name: name || getPayloadName(normalizedPath),
          bytes: toTransferableBuffer(bytes)
        });
      };
      const stageTextFile = (path, text) => {
        throwIfGenerationCanceled();
        const bytes = textEncoder.encode(String(text ?? ''));
        pyodide.FS.writeFile(path, bytes);
        stageGenerationBytes(path, getPayloadName(path), bytes);
      };
      const stageUploadedFile = async (fileObj, path, { cacheText = false, textCache = null } = {}) => {
        if (!fileObj) return false;
        throwIfGenerationCanceled();
        const buffer = await fileObj.arrayBuffer();
        throwIfGenerationCanceled();
        const bytes = new Uint8Array(buffer);
        pyodide.FS.writeFile(path, bytes);
        if (cacheText && textCache) {
          textCache.set(fileObj, textDecoder.decode(bytes));
        }
        stageGenerationBytes(path, fileObj.name || getPayloadName(path), buffer);
        return true;
      };

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
      stageTextFile('/combined_d.tsv', dContent);
      args.push('-d', '/combined_d.tsv');

      let tContent = '';
      manualSpecificRules.forEach((r) => {
        tContent += `${r.feat}\t${r.qual}\t${r.val}\t${r.color}\t${r.cap}\n`;
      });
      if (tContent.trim() !== '') {
        stageTextFile('/combined_t.tsv', tContent);
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
          stageTextFile('/manual_wl.tsv', wlContent);
          args.push('--label_whitelist', '/manual_wl.tsv');
        }
      }

      let pContent = '';
      manualPriorityRules.forEach((r) => {
        pContent += `${r.feat}\t${r.order}\n`;
      });
      if (pContent.trim() !== '') {
        stageTextFile('/priority.tsv', pContent);
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
        stageTextFile('/web_label_table.tsv', labelOverride.tsv);
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
        stageTextFile('/web_feature_table.tsv', `${featureVisibilityRows.join('\n')}\n`);
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
      const appendDepthStyleArgs = () => {
        if (adv.depth_color) args.push('--depth_color', adv.depth_color);
        if (adv.depth_min !== null && adv.depth_min !== undefined && adv.depth_min !== '') {
          args.push('--depth_min', adv.depth_min);
        }
        if (adv.depth_max !== null && adv.depth_max !== undefined && adv.depth_max !== '') {
          args.push('--depth_max', adv.depth_max);
        }
        if (
          adv.depth_min !== null &&
          adv.depth_min !== undefined &&
          adv.depth_min !== '' &&
          adv.depth_max !== null &&
          adv.depth_max !== undefined &&
          adv.depth_max !== '' &&
          Number(adv.depth_min) > Number(adv.depth_max)
        ) {
          throw new Error('Depth minimum must be less than or equal to depth maximum.');
        }
        if (adv.depth_normalize === true) args.push('--depth_log_scale');
        if (adv.depth_share_axis === true) args.push('--share_depth_axis');
        if (adv.depth_show_axis === false) args.push('--hide_depth_axis');
        if (
          adv.depth_window_size !== null &&
          adv.depth_window_size !== undefined &&
          adv.depth_window_size !== ''
        ) {
          if (Number(adv.depth_window_size) <= 0) throw new Error('Depth window must be greater than 0.');
          args.push('--depth_window', adv.depth_window_size);
        }
        if (
          adv.depth_step_size !== null &&
          adv.depth_step_size !== undefined &&
          adv.depth_step_size !== ''
        ) {
          if (Number(adv.depth_step_size) <= 0) throw new Error('Depth step must be greater than 0.');
          args.push('--depth_step', adv.depth_step_size);
        }
        if (adv.depth_show_ticks === false) args.push('--hide_depth_ticks');
        if (
          adv.depth_tick_interval !== null &&
          adv.depth_tick_interval !== undefined &&
          adv.depth_tick_interval !== ''
        ) {
          if (Number(adv.depth_tick_interval) <= 0) throw new Error('Depth large tick interval must be greater than 0.');
          args.push('--depth_large_tick_interval', adv.depth_tick_interval);
        }
        if (
          mode.value === 'circular' &&
          adv.depth_small_tick_interval !== null &&
          adv.depth_small_tick_interval !== undefined &&
          adv.depth_small_tick_interval !== ''
        ) {
          if (Number(adv.depth_small_tick_interval) <= 0) throw new Error('Depth small tick interval must be greater than 0.');
          args.push('--depth_small_tick_interval', adv.depth_small_tick_interval);
        }
        if (
          adv.depth_tick_font_size !== null &&
          adv.depth_tick_font_size !== undefined &&
          adv.depth_tick_font_size !== ''
        ) {
          if (Number(adv.depth_tick_font_size) <= 0) throw new Error('Depth tick font size must be greater than 0.');
          args.push('--depth_tick_font_size', adv.depth_tick_font_size);
        }
      };
      const appendGcContentPercentArgs = () => {
        adv.gc_content_mode = String(adv.gc_content_mode || '').trim().toLowerCase() === 'percent'
          ? 'percent'
          : 'deviation';
        if (adv.gc_content_mode !== 'percent') return;

        args.push('--gc_content_mode', 'percent');
        const minPercent = Number(adv.gc_content_min_percent);
        const maxPercent = Number(adv.gc_content_max_percent);
        if (!Number.isFinite(minPercent)) {
          throw new Error('GC content minimum percent must be a finite number.');
        }
        if (!Number.isFinite(maxPercent)) {
          throw new Error('GC content maximum percent must be a finite number.');
        }
        if (minPercent > maxPercent) {
          throw new Error('GC content minimum percent must be less than or equal to maximum percent.');
        }
        adv.gc_content_min_percent = minPercent;
        adv.gc_content_max_percent = maxPercent;
        args.push('--gc_content_min_percent', String(minPercent));
        args.push('--gc_content_max_percent', String(maxPercent));

        if (adv.gc_content_show_axis === false) args.push('--hide_gc_content_axis');
        if (adv.gc_content_show_ticks === false) args.push('--hide_gc_content_ticks');
        if (
          adv.gc_content_tick_interval !== null &&
          adv.gc_content_tick_interval !== undefined &&
          adv.gc_content_tick_interval !== ''
        ) {
          if (Number(adv.gc_content_tick_interval) <= 0) throw new Error('GC content large tick interval must be greater than 0.');
          args.push('--gc_content_large_tick_interval', adv.gc_content_tick_interval);
        }
        if (
          adv.gc_content_small_tick_interval !== null &&
          adv.gc_content_small_tick_interval !== undefined &&
          adv.gc_content_small_tick_interval !== ''
        ) {
          if (Number(adv.gc_content_small_tick_interval) <= 0) throw new Error('GC content small tick interval must be greater than 0.');
          args.push('--gc_content_small_tick_interval', adv.gc_content_small_tick_interval);
        }
        if (
          adv.gc_content_tick_font_size !== null &&
          adv.gc_content_tick_font_size !== undefined &&
          adv.gc_content_tick_font_size !== ''
        ) {
          if (Number(adv.gc_content_tick_font_size) <= 0) throw new Error('GC content tick font size must be greater than 0.');
          args.push('--gc_content_tick_font_size', adv.gc_content_tick_font_size);
        }
      };
      appendGcContentPercentArgs();

      if (mode.value === 'circular') {
        const multiCanvasSupport = getCircularMultiRecordCanvasOptionSupport();
        const normalizedCircularPlotTitle = String(form.plot_title || '').trim();
        const normalizedPlotTitlePosition = normalizeCircularPlotTitlePosition(adv.plot_title_position);
        const useCircularTrackSlots = adv.circular_track_slots_enabled === true;
        const circularTrackAxisIndex = clampCircularTrackAxisIndex(
          adv.circular_track_slots_axis_index,
          Array.isArray(adv.circular_track_slots) ? adv.circular_track_slots.length : 0
        );
        const circularTrackSlots = useCircularTrackSlots
          ? applyCircularTrackOrderPlacements(
              adv.circular_track_slots,
              adv.nt,
              form.track_type,
              circularTrackAxisIndex
            )
          : [];
        if (useCircularTrackSlots) {
          if (!multiCanvasSupport.circular_track_slot || !multiCanvasSupport.circular_track_axis_index) {
            throw new Error(
              'Current gbdraw wheel does not support --circular_track_slot and --circular_track_axis_index. Rebuild and redeploy the web wheel.'
            );
          }
          adv.circular_track_slots.splice(0, adv.circular_track_slots.length, ...circularTrackSlots);
          const circularTrackOnAxisIndex = circularTrackSlots.findIndex((slot) => slot?.side === 'overlay');
          const normalizedCircularTrackAxisIndex = circularTrackOnAxisIndex >= 0
            ? circularTrackOnAxisIndex
            : (
                circularTrackAxisIndex === null
                  ? inferLegacyAxisIndexFromFeature(circularTrackSlots, form.track_type)
                  : circularTrackAxisIndex
              );
          adv.circular_track_slots_axis_index = clampCircularTrackAxisIndex(
            normalizedCircularTrackAxisIndex,
            circularTrackSlots.length
          );
        }
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
        const normalizedCenterReservedRadius = normalizeNonNegativeNumberOrNull(adv.center_reserved_radius);
        form.plot_title = normalizedCircularPlotTitle;
        adv.plot_title_position = normalizedPlotTitlePosition;
        adv.plot_title_font_size = normalizedPlotTitleFontSize;
        adv.keep_full_definition_with_plot_title = keepFullDefinitionWithPlotTitle;
        adv.center_reserved_radius = normalizedCenterReservedRadius;

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
        args.push('--track_type', form.track_type);
        args.push('-l', form.legend);
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
        if (useCircularTrackSlots && normalizedCenterReservedRadius !== null) {
          if (!multiCanvasSupport.center_reserved_radius) {
            throw new Error(
              'Current gbdraw wheel does not support --center_reserved_radius. Rebuild and redeploy the web wheel.'
            );
          }
          args.push('--center_reserved_radius', String(normalizedCenterReservedRadius));
        }
        const labelsModeRaw =
          typeof form.labels_mode === 'string'
            ? form.labels_mode
            : (form.allow_inner_labels ? 'both' : (form.show_labels ? 'out' : 'none'));
        const labelsMode = String(labelsModeRaw || 'none').trim().toLowerCase();
        if (labelsMode === 'out') args.push('--labels');
        if (labelsMode === 'both') args.push('--labels', 'both');
        const normalizedLabelRendering = labelsMode === 'none' ? 'auto' : normalizeLabelRendering(adv.label_rendering);
        adv.label_rendering = normalizedLabelRendering;
        if (normalizedLabelRendering !== 'auto') {
          if (!multiCanvasSupport.label_rendering) {
            throw new Error(
              'Current gbdraw wheel does not support --label_rendering. Rebuild and redeploy the web wheel.'
            );
          }
          args.push('--label_rendering', normalizedLabelRendering);
        }
        if (!useCircularTrackSlots) {
          if (form.suppress_gc) args.push('--suppress_gc');
          if (form.suppress_skew) args.push('--suppress_skew');
        }
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
        const normalizedCircularLabelSpacing = normalizePositiveNumberOrNull(adv.circular_label_spacing);
        adv.circular_label_spacing = normalizedCircularLabelSpacing;
        if (normalizedCircularLabelSpacing !== null) {
          if (!multiCanvasSupport.circular_label_spacing) {
            throw new Error(
              'Current gbdraw wheel does not support --circular_label_spacing. Rebuild and redeploy the web wheel.'
            );
          }
          args.push('--circular_label_spacing', normalizedCircularLabelSpacing);
        }
        if (
          !useCircularTrackSlots &&
          adv.feature_width_circular !== null &&
          adv.feature_width_circular !== undefined &&
          adv.feature_width_circular !== '' &&
          Number(adv.feature_width_circular) > 0
        ) {
          args.push('--feature_width', adv.feature_width_circular);
        }
        if (!useCircularTrackSlots && !form.suppress_gc) {
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
        if (!useCircularTrackSlots && !form.suppress_skew) {
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
        if (useCircularTrackSlots) {
          args.push('--circular_track_axis_index', String(adv.circular_track_slots_axis_index));
          circularTrackSlots.forEach((slot) => {
            args.push(
              '--circular_track_slot',
              buildCircularTrackSlotSpec(slot, adv.nt, form.track_type, {
                includeSide: false,
                forceSplitLane: true
              })
            );
          });
        }
        const hasCircularDepthFile = Boolean(files.c_depth);
        const circularSlotNeedsDepth = useCircularTrackSlots && hasEnabledCircularTrackRenderer(circularTrackSlots, 'depth');
        if (form.show_depth || circularSlotNeedsDepth) {
          if (!hasCircularDepthFile) throw new Error('Please upload a Depth TSV file or disable Show depth track.');
          await stageUploadedFile(files.c_depth, '/depth.tsv');
          args.push('--depth', '/depth.tsv');
          args.push('--show_depth');
          appendDepthStyleArgs();
          if (
            !useCircularTrackSlots &&
            adv.depth_width_circular !== null &&
            adv.depth_width_circular !== undefined &&
            adv.depth_width_circular !== '' &&
            Number(adv.depth_width_circular) > 0
          ) {
            args.push('--depth_width', adv.depth_width_circular);
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
          await stageUploadedFile(files.c_gb, '/input.gb');
          args.push('--gbk', '/input.gb');
        } else {
          if (!files.c_gff || !files.c_fasta) throw new Error('GFF3 and FASTA are required.');
          await stageUploadedFile(files.c_gff, '/input.gff');
          await stageUploadedFile(files.c_fasta, '/input.fasta');
          args.push('--gff', '/input.gff', '--fasta', '/input.fasta');
        }

        const sourceMode = String(circularConservation.source || '').trim().toLowerCase() === 'upload'
          ? 'upload'
          : 'losat';
        const circularConservationSourceFiles = sourceMode === 'upload'
          ? normalizeFileList(files.c_conservation_blasts)
          : normalizeFileList(files.c_conservation_fastas);
        const shouldDrawCircularPairwiseComparisons = circularConservationSourceFiles.length > 0;
        circularConservation.enabled = shouldDrawCircularPairwiseComparisons;

        if (shouldDrawCircularPairwiseComparisons) {
          if (!multiCanvasSupport.conservation_blast || !multiCanvasSupport.conservation_reference) {
            throw new Error(
              'Current gbdraw wheel does not support circular pairwise comparison rings. Rebuild and redeploy the web wheel.'
            );
          }
          const width = normalizePositiveNumberOrNull(circularConservation.ring_width);
          const gap = normalizePositiveNumberOrNull(circularConservation.ring_gap);
          const writeEmptyConservationPreflightFiles = (count) => {
            const paths = [];
            for (let index = 0; index < count; index += 1) {
              const path = `/conservation_preflight_${index}.txt`;
              pyodide.FS.writeFile(path, new Uint8Array());
              paths.push(path);
            }
            return paths;
          };
          const runConservationLayoutPreflight = (blastPaths, reference) => {
            if (isReflow || !Array.isArray(blastPaths) || blastPaths.length === 0) return true;
            throwIfGenerationCanceled();
            setProcessingStatus('Checking circular track layout...');
            const preflightArgs = [
              ...args,
              '--conservation_blast',
              ...blastPaths,
              '--conservation_reference',
              reference
            ];
            const preflight = runCircularLayoutPreflight(preflightArgs);
            if (!preflight?.error) return true;
            const formatted = formatPythonError(preflight.error);
            if (isReflow) {
              labelReflowLastError.value = formatted?.summary || 'Auto reflow failed';
            } else {
              errorLog.value = formatted;
            }
            return false;
          };
          const appendConservationStyleArgs = (series) => {
            const labels = series.map((entry) => entry.label);
            const colors = series.map((entry) => entry.color);
            if (labels.length > 0) {
              if (!multiCanvasSupport.conservation_labels) {
                throw new Error(
                  'Current gbdraw wheel does not support --conservation_labels. Rebuild and redeploy the web wheel.'
                );
              }
              args.push('--conservation_labels', ...labels);
            }
            if (colors.length > 0) {
              if (!multiCanvasSupport.conservation_colors) {
                throw new Error(
                  'Current gbdraw wheel does not support --conservation_colors. Rebuild and redeploy the web wheel.'
                );
              }
              args.push('--conservation_colors', ...colors);
            }
            if (width !== null) {
              if (!multiCanvasSupport.conservation_ring_width) {
                throw new Error(
                  'Current gbdraw wheel does not support --conservation_ring_width. Rebuild and redeploy the web wheel.'
                );
              }
              args.push('--conservation_ring_width', String(width));
            }
            if (gap !== null) {
              if (!multiCanvasSupport.conservation_ring_gap) {
                throw new Error(
                  'Current gbdraw wheel does not support --conservation_ring_gap. Rebuild and redeploy the web wheel.'
                );
              }
              args.push('--conservation_ring_gap', String(gap));
            }
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
          };

          const runCircularLosatConservation = async (comparisonFiles) => {
            const circularLosatProgram = normalizeCircularConservationLosatProgram(
              circularConservation.losat_program
            );
            circularConservation.losat_program = circularLosatProgram;
            const queryGencode = normalizePositiveInteger(circularConservation.query_gencode, 1);
            const subjectGencode = normalizePositiveInteger(circularConservation.subject_gencode, 1);
            circularConservation.query_gencode = queryGencode;
            circularConservation.subject_gencode = subjectGencode;
            const extraArgs = [];
            if (circularLosatProgram === 'tblastx') {
              extraArgs.push('--query-gencode', String(queryGencode));
              extraArgs.push('--db-gencode', String(subjectGencode));
            } else {
              const normalizedTask = String(losat.blastn?.task || 'megablast').trim() || 'megablast';
              extraArgs.push('--task', normalizedTask);
            }
            const circularLosatSuffix = circularLosatProgram === 'tblastx' ? 'tlosatx' : 'losatn';
            const subjectFile = cInputType.value === 'gb' ? files.c_gb : files.c_fasta;
            const subjectFmt = cInputType.value === 'gb' ? 'genbank' : 'fasta';
            const subjectEntry = await extractAllLosatFastaFast({
              file: subjectFile,
              fmt: subjectFmt
            });
            const subjectHash = await hashText(subjectEntry.fasta);
            const subjectSequenceKey = `circular-subject:${subjectHash}`;
            const sequenceEntriesByKey = new Map([[subjectSequenceKey, subjectEntry.fasta]]);
            const cacheMap = losatCache.value || new Map();
            const cacheInfo = [];
            const losatPairs = [];
            const losatJobs = [];
            const pendingJobKeys = new Set();
            const executionMode = getLosatExecutionMode();
            const runtimeCompatibility = executionMode === 'serial'
              ? 'serial-v1'
              : 'threaded-compatible-v1';

            for (let index = 0; index < comparisonFiles.length; index += 1) {
              throwIfGenerationCanceled();
              const fileObj = comparisonFiles[index];
              const queryFasta = await fileObj.text();
              if (getFastaSequenceLength(queryFasta) <= 0) {
                throw new Error(`Pairwise comparison FASTA #${index + 1} has no sequence data.`);
              }
              const queryHash = await hashText(queryFasta);
              const querySequenceKey = `circular-query:${queryHash}`;
              sequenceEntriesByKey.set(querySequenceKey, queryFasta);
              const cacheKey = await hashText(JSON.stringify({
                cacheSchema: LOSAT_CACHE_SCHEMA,
                flow: 'circular-conservation',
                losatRuntimeCompatibility: runtimeCompatibility,
                losatThreadsPerJob: 1,
                program: circularLosatProgram,
                outfmt: String(losat.outfmt || '6'),
                args: extraArgs,
                queryCanonicalHash: queryHash,
                subjectCanonicalHash: subjectHash
              }));
              const fallbackName = makeSafeFilename(
                `${String(fileObj?.name || `comparison_${index + 1}`).replace(/\.[^.]+$/, '')}.circular_conservation.${circularLosatSuffix}.tsv`
              );
              const pair = {
                sourceIndex: index,
                cacheKey,
                filename: fallbackName
              };
              losatPairs.push(pair);
              cacheInfo.push({
                key: cacheKey,
                filename: fallbackName,
                display: true
              });
              if (!isRawLosatCacheEntry(cacheMap.get(cacheKey)) && !pendingJobKeys.has(cacheKey)) {
                pendingJobKeys.add(cacheKey);
                losatJobs.push({
                  pairIndex: index,
                  cacheKey,
                  program: circularLosatProgram,
                  querySequenceKey,
                  subjectSequenceKey,
                  queryCanonicalHash: queryHash,
                  subjectCanonicalHash: subjectHash,
                  outfmt: losat.outfmt || '6',
                  extraArgs
                });
              }
            }

            if (losatJobs.length > 0) {
              setProcessingStatus(`Running ${circularLosatSuffix.toUpperCase()} conservation: 0/${losatJobs.length} jobs complete`);
              const runtime = await prepareLosatRuntime({ includeThreaded: executionMode !== 'serial' }).catch((error) => {
                console.warn('LOSAT runtime warmup failed; execution will report the error.', error);
                return null;
              });
              if (runtime?.threaded && losatThreadingStatus) {
                const { wasmModule: _wasmModule, ...threadedStatus } = runtime.threaded;
                losatThreadingStatus.value = threadedStatus;
              }
              const losatResults = await runLosatPairsParallel(losatJobs, {
                concurrency: getLosatParallelWorkers(),
                executionMode,
                totalThreadBudget: getLosatTotalThreadBudget(),
                threadsPerJob: 1,
                sequences: sequenceEntriesByKey,
                signal: generationAbortSignal,
                onRuntimeStatus: (status) => {
                  losatThreadingStatus.value = status;
                },
                onProgress: ({ completed, total }) => {
                  if (generationAbortSignal?.aborted || generationCancelRequested.value) return;
                  setProcessingStatus(`Running ${circularLosatSuffix.toUpperCase()} conservation: ${completed}/${total} jobs complete`);
                }
              });
              losatResults.forEach((result) => {
                const job = losatJobs.find((item) => item.cacheKey === result.cacheKey);
                cacheMap.set(result.cacheKey, {
                  schema: LOSAT_CACHE_SCHEMA,
                  kind: 'raw-losat',
                  text: result.text,
                  program: circularLosatProgram,
                  queryCanonicalHash: job?.queryCanonicalHash || '',
                  subjectCanonicalHash: job?.subjectCanonicalHash || ''
                });
              });
            } else {
              setProcessingStatus('Using cached LOSAT conservation results...');
            }

            const paths = [];
            for (const pair of losatPairs) {
              const cached = cacheMap.get(pair.cacheKey);
              const blastText = isRawLosatCacheEntry(cached) ? cached.text : '';
              const blastPath = `/conservation_blast_${pair.sourceIndex}.txt`;
              stageTextFile(blastPath, blastText);
              paths.push(blastPath);
            }
            losatCacheInfo.value = cacheInfo;
            losatCache.value = cacheMap;
            return paths;
          };

          let conservationBlastPaths = [];
          let conservationReference = 'auto';
          let conservationSeries = [];
          if (sourceMode === 'upload') {
            const blastFiles = circularConservationSourceFiles;
            if (blastFiles.length === 0) {
              throw new Error('Please upload at least one BLAST outfmt 6/7 file for Pairwise Comparisons.');
            }
            const conservationEntries = orderedConservationSources(blastFiles, circularConservation);
            conservationSeries = buildConservationSeries(blastFiles, circularConservation);
            for (let index = 0; index < conservationEntries.length; index += 1) {
              const blastPath = `/conservation_blast_${index}.txt`;
              await stageUploadedFile(conservationEntries[index].file, blastPath);
              conservationBlastPaths.push(blastPath);
            }
            const rawReference = String(circularConservation.reference || 'auto').trim().toLowerCase();
            conservationReference = ['query', 'subject'].includes(rawReference) ? rawReference : 'auto';
            losatCacheInfo.value = [];
            appendConservationStyleArgs(conservationSeries);
            if (!runConservationLayoutPreflight(conservationBlastPaths, conservationReference)) {
              return { status: 'error' };
            }
          } else {
            const comparisonFiles = circularConservationSourceFiles;
            if (comparisonFiles.length === 0) {
              throw new Error('Please upload at least one comparison FASTA file for Pairwise Comparisons.');
            }
            const conservationEntries = orderedConservationSources(comparisonFiles, circularConservation);
            const orderedComparisonFiles = conservationEntries.map((entry) => entry.file);
            conservationSeries = buildConservationSeries(comparisonFiles, circularConservation);
            appendConservationStyleArgs(conservationSeries);
            const preflightBlastPaths = writeEmptyConservationPreflightFiles(orderedComparisonFiles.length);
            if (!runConservationLayoutPreflight(preflightBlastPaths, 'subject')) {
              return { status: 'error' };
            }
            conservationBlastPaths = await runCircularLosatConservation(orderedComparisonFiles);
            conservationReference = 'subject';
          }

          args.push('--conservation_blast', ...conservationBlastPaths);
          args.push('--conservation_reference', conservationReference);
        } else {
          losatCacheInfo.value = [];
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
        const normalizedPairwiseMatchStyle = normalizePairwiseMatchStyle(adv.pairwise_match_style);
        adv.pairwise_match_style = normalizedPairwiseMatchStyle;
        if (form.align_center) args.push('--align_center');
        if (form.show_gc) args.push('--show_gc');
        if (form.show_skew) args.push('--show_skew');
        if (form.normalize_length) args.push('--normalize_length');
        if (form.legend !== 'right') args.push('-l', form.legend);
        const useLosat = blastSource.value === 'losat';
        const useProteinBlastp = useLosat && losatProgram.value === 'blastp';
        const blastpMode = normalizeBlastpMode(losat.blastp?.mode);
        const useOrthogroupBlastp = useProteinBlastp && blastpMode === 'orthogroup';
        const useCollinearBlastp = useProteinBlastp && blastpMode === 'collinear';
        const selectedOrthogroupTarget = String(selectedOrthogroupAlignmentFeature.value || '').trim();
        const wantsOrthogroupAlignmentOption =
          useOrthogroupBlastp && lInputType.value === 'gb' && selectedOrthogroupTarget !== '';
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
        const blastpDisplayMaxHits = normalizeBlastThresholdNumber(
          losat.blastp?.maxHits,
          5,
          { integer: true }
        );
        losat.blastp.mode = blastpMode;
        losat.blastp.maxHits = Math.max(1, blastpDisplayMaxHits);
        losat.blastp.candidateLimit = null;
        losat.blastp.collinearMinAnchors = Math.max(
          1,
          normalizeBlastThresholdNumber(losat.blastp?.collinearMinAnchors, 1, { integer: true })
        );
        losat.blastp.collinearMaxGeneGap = Math.max(
          0,
          normalizeBlastThresholdNumber(losat.blastp?.collinearMaxGeneGap, 0, { integer: true })
        );
        losat.blastp.collinearMaxDiagonalDrift = Math.max(
          0,
          normalizeBlastThresholdNumber(losat.blastp?.collinearMaxDiagonalDrift, 0, { integer: true })
        );
        losat.blastp.collinearColorMode = normalizeCollinearColorMode(losat.blastp?.collinearColorMode);
        losat.blastp.collinearAnchorMode = normalizeCollinearAnchorMode(losat.blastp?.collinearAnchorMode);
        losat.blastp.collinearSearchScope = normalizeCollinearSearchScope(losat.blastp?.collinearSearchScope);
        const collinearSearchScope = losat.blastp.collinearSearchScope;
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
        let normalizedLabelRendering = form.show_labels_linear === 'none' ? 'auto' : normalizeLabelRendering(adv.label_rendering);
        if (normalizedLabelPlacement === 'above_feature') {
          normalizedLabelRendering = 'auto';
        }
        adv.label_rendering = normalizedLabelRendering;
        const normalizedLinearLabelSpacing = normalizePositiveNumberOrNull(adv.linear_label_spacing);
        adv.linear_label_spacing = normalizedLinearLabelSpacing;
        const wantsPlacementOption = normalizedLabelPlacement && normalizedLabelPlacement !== 'auto';
        const wantsLabelRenderingOption = normalizedLabelRendering !== 'auto';
        const wantsRotationOption = adv.label_rotation !== null && adv.label_rotation !== undefined && adv.label_rotation !== '';
        const wantsLinearLabelSpacingOption = normalizedLinearLabelSpacing !== null;
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
        const wantsPairwiseMatchStyleOption = normalizedPairwiseMatchStyle !== 'ribbon';
        const wantsKeepDefinitionLeftAlignedOption = form.keep_definition_left_aligned === true;
        const wantsRulerOnAxisOption =
          Boolean(form.linear_ruler_on_axis) &&
          form.scale_style === 'ruler' &&
          (normalizedTrackLayout === 'above' || normalizedTrackLayout === 'below');
        const linearLabelSupport = getLinearLabelOptionSupport();
        if (
          wantsPlacementOption ||
          wantsLabelRenderingOption ||
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
          wantsHideLengthOption ||
          wantsPairwiseMatchStyleOption ||
          wantsOrthogroupAlignmentOption ||
          wantsKeepDefinitionLeftAlignedOption
        ) {
          if (wantsPlacementOption && !linearLabelSupport.placement) {
            throw new Error("Current gbdraw wheel does not support --label_placement. Rebuild and redeploy the web wheel.");
          }
          if (wantsLabelRenderingOption && !linearLabelSupport.label_rendering) {
            throw new Error("Current gbdraw wheel does not support --label_rendering. Rebuild and redeploy the web wheel.");
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
          if (wantsPairwiseMatchStyleOption && !linearLabelSupport.pairwise_match_style) {
            throw new Error("Current gbdraw wheel does not support --pairwise_match_style. Rebuild and redeploy the web wheel.");
          }
          if (wantsOrthogroupAlignmentOption && !linearLabelSupport.orthogroup_alignment) {
            throw new Error("Current gbdraw wheel does not support --align_orthogroup_feature. Rebuild and redeploy the web wheel.");
          }
          if (wantsKeepDefinitionLeftAlignedOption && !linearLabelSupport.keep_definition_left_aligned) {
            throw new Error("Current gbdraw wheel does not support --keep_definition_left_aligned. Rebuild and redeploy the web wheel.");
          }
        }
        if (wantsPlotTitleOption) args.push('--plot_title', normalizedPlotTitle);
        if (wantsPlotTitlePositionOption) args.push('--plot_title_position', normalizedPlotTitlePosition);
        if (wantsPlotTitleFontSizeOption) args.push('--plot_title_font_size', String(normalizedPlotTitleFontSize));
        if (wantsShowRepliconOption) args.push('--show_replicon');
        if (wantsHideAccessionOption) args.push('--hide_accession');
        if (wantsHideLengthOption) args.push('--hide_length');
        if (wantsPairwiseMatchStyleOption) args.push('--pairwise_match_style', normalizedPairwiseMatchStyle);
        if (wantsOrthogroupAlignmentOption) {
          args.push('--align_orthogroup_feature', selectedOrthogroupTarget);
        }
        if (wantsKeepDefinitionLeftAlignedOption) {
          args.push('--keep_definition_left_aligned');
        }
        if (normalizedLabelPlacement && normalizedLabelPlacement !== 'auto') {
          args.push('--label_placement', normalizedLabelPlacement);
        }
        if (normalizedLabelRendering !== 'auto') {
          args.push('--label_rendering', normalizedLabelRendering);
        }
        if (adv.label_rotation !== null && adv.label_rotation !== undefined && adv.label_rotation !== '') {
          args.push('--label_rotation', adv.label_rotation);
        }
        if (normalizedLinearLabelSpacing !== null) {
          args.push('--linear_label_spacing', normalizedLinearLabelSpacing);
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

        const viewTransformSpecs = [];
        const buildRegionSpec = (seq, idx) => {
          const hasStart = seq.region_start !== null && seq.region_start !== undefined && seq.region_start !== '';
          const hasEnd = seq.region_end !== null && seq.region_end !== undefined && seq.region_end !== '';
          const recordIdRaw = seq.region_record_id ? String(seq.region_record_id).trim() : '';
          const wantsReverse = Boolean(seq.region_reverse);
          if (hasStart !== hasEnd) {
            throw new Error(`Sequence #${idx + 1}: Provide both Region start and end, or leave both empty.`);
          }

          recordSelectors.push(recordIdRaw || '');

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
            const canonicalStart = Math.min(start, end);
            const canonicalEnd = Math.max(start, end);
            const coordinateReverse = start > end;
            const displayReverse = wantsReverse || coordinateReverse;
            const specBody = `${start}-${end}${wantsReverse ? ':rc' : ''}`;
            const canonicalSpecBody = `${canonicalStart}-${canonicalEnd}`;
            const fileLabel = lInputType.value === 'gb' ? `seq_${idx}.gb` : `seq_${idx}.gff`;
            const cliSpec = recordIdRaw ? `${fileLabel}:${recordIdRaw}:${specBody}` : `#${idx + 1}:${specBody}`;
            const fileSpec = canonicalSpecBody;
            reverseFlags.push(false);
            viewTransformSpecs.push({ reverse: displayReverse });
            return { cli: cliSpec, file: fileSpec, displayFile: specBody };
          }

          reverseFlags.push(wantsReverse);
          viewTransformSpecs.push({ reverse: wantsReverse });
          return null;
        };

        let inputArgs = [];
        let blastArgs = [];
        const fastaCache = new Map();
        const fastaHashCache = new Map();
        const sequenceEntriesByKey = new Map();
        const linearFileTextCache = new WeakMap();
        let extractFirstFasta = null;
        let extractProteinFasta = null;
        let convertProteinBlast = null;
        let convertNucleotideBlast = null;
        if (useOrthogroupBlastp) {
          clearOrthogroupMetadata();
        } else {
          clearOrthogroupMetadata({ clearSelection: true });
        }
        let cacheInfo = [];
        const cacheMap = losatCache.value || new Map();
        const losatTiming = useLosat
          ? {
              inputWriteMs: 0,
              fastaExtractionMs: 0,
              cacheHashMs: 0,
              jobBuildMs: 0,
              jobBuildWallMs: 0,
              runtimeWaitMs: 0,
              executionMs: 0,
              blastWriteMs: 0,
              totalFastaChars: 0,
              fastaCacheHits: 0,
              proteinExtractionCacheHits: 0,
              fastaJsExtractions: 0,
              fastaPyodideFallbacks: 0,
              proteinConversionCacheHits: 0,
              proteinFilteredHitCacheHits: 0,
              proteinFilteredHitCacheMisses: 0,
              cacheHashHits: 0,
              cacheHits: 0,
              cacheMisses: 0,
              totalPairs: 0,
              uniqueJobs: 0
            }
          : null;
        const losatExecutionMode = useLosat ? getLosatExecutionMode() : 'serial';
        const losatRequestedThreadsPerJob = useLosat ? getLosatThreadsPerJob() : undefined;
        const losatRequestedTotalThreadBudget = useLosat ? getLosatTotalThreadBudget() : undefined;
        const losatRuntimeWarmup = useLosat
          ? prepareLosatRuntime({ includeThreaded: losatExecutionMode !== 'serial' }).then((runtime) => {
              if (runtime?.threaded && losatThreadingStatus) {
                const { wasmModule: _wasmModule, ...threadedStatus } = runtime.threaded;
                losatThreadingStatus.value = threadedStatus;
              }
              return runtime;
            }).catch((error) => {
              console.warn('LOSAT runtime warmup failed; execution will report the error if LOSAT is used.', error);
              return null;
            })
          : null;

        if (useLosat) {
          if (useProteinBlastp) {
            extractProteinFasta = pyodide.globals.get('extract_cds_protein_fasta');
            convertProteinBlast = pyodide.globals.get('convert_losatp_blastp_pairs_to_genomic_payload');
          } else {
            extractFirstFasta = pyodide.globals.get('extract_first_fasta');
            convertNucleotideBlast = pyodide.globals.get('convert_losat_nucleotide_to_display_tsv');
          }
        } else {
          losatCacheInfo.value = [];
        }

        let proteinRecordInstanceKeys = [];

        const buildProteinRecordInstanceKeys = async () => {
          const occurrences = new Map();
          const keys = [];
          for (let idx = 0; idx < linearSeqs.length; idx += 1) {
            const sourceFile = lInputType.value === 'gb' ? linearSeqs[idx]?.gb : linearSeqs[idx]?.gff;
            const pairedFastaFile = lInputType.value === 'gff' ? linearSeqs[idx]?.fasta : null;
            const descriptor = {
              inputFormat: lInputType.value,
              primaryFile: getFileFingerprint(sourceFile),
              pairedFastaFile: getFileFingerprint(pairedFastaFile),
              recordSelector: recordSelectors[idx] ?? '',
              canonicalRegionSpec: regionSpecs[idx]?.file || null
            };
            const base = (await hashText(JSON.stringify(descriptor))).slice(0, 16);
            const occurrence = (occurrences.get(base) || 0) + 1;
            occurrences.set(base, occurrence);
            keys.push(occurrence > 1 ? `r_${base}_o${occurrence}` : `r_${base}`);
          }
          return keys;
        };

        const getSeqEntry = async (idx) => {
          if (fastaCache.has(idx)) return fastaCache.get(idx);
          const startedAt = getNow();
          const path = lInputType.value === 'gb'
            ? `/seq_${idx}.gb`
            : (useProteinBlastp ? `/seq_${idx}.gff` : `/seq_${idx}.fasta`);
          const fmt = lInputType.value === 'gb'
            ? 'genbank'
            : (useProteinBlastp ? 'gff' : 'fasta');
          const pairedFastaPath = lInputType.value === 'gff' && useProteinBlastp
            ? `/seq_${idx}.fasta`
            : null;
          const regionSpec = regionSpecs[idx]?.file || null;
          const recordSelector = recordSelectors[idx] ?? '';
          const reverseFlag = '0';
          const sourceFile = lInputType.value === 'gb'
            ? linearSeqs[idx]?.gb
            : (useProteinBlastp ? linearSeqs[idx]?.gff : linearSeqs[idx]?.fasta);
          const sourceText = sourceFile ? linearFileTextCache.get(sourceFile) : null;
          const pairedFastaFile = lInputType.value === 'gff' && useProteinBlastp
            ? linearSeqs[idx]?.fasta
            : null;
          const recordInstanceKey = useProteinBlastp
            ? (proteinRecordInstanceKeys[idx] || `r_${idx + 1}`)
            : '';
          const persistentCacheKey = useProteinBlastp
            ? JSON.stringify({
                inputFormat: lInputType.value,
                primaryFile: getFileFingerprint(sourceFile),
                pairedFastaFile: getFileFingerprint(pairedFastaFile),
                regionSpec,
                recordSelector,
                recordInstanceKey,
                recordIndex: idx
              })
            : JSON.stringify({ fmt, regionSpec, recordSelector, reverseFlag });
          const usePersistentFastaCache = !useProteinBlastp;
          const cachedEntry = sourceFile
            ? (
                useProteinBlastp
                  ? getCachedProteinExtraction(sourceFile, persistentCacheKey)
                  : (usePersistentFastaCache ? getCachedFastaExtraction(sourceFile, persistentCacheKey) : null)
              )
            : null;
          let entry = cachedEntry;

          if (entry) {
            if (losatTiming) {
              if (useProteinBlastp) losatTiming.proteinExtractionCacheHits += 1;
              else losatTiming.fastaCacheHits += 1;
            }
          } else {
            if (useProteinBlastp) {
              const res = JSON.parse(
                extractProteinFasta(path, fmt, pairedFastaPath, regionSpec, recordSelector, reverseFlag, idx, recordInstanceKey)
              );
              if (res.error) throw new Error(res.error);
              const fastaHash = await hashText(res.fasta || '');
              const proteinCacheKey = await hashText(JSON.stringify({ extraction: persistentCacheKey, fastaHash }));
              entry = {
                fasta: res.fasta,
                recordId: res.record_id || `seq_${idx + 1}`,
                canonicalLength: Number(res.record_length || 0) || getFastaSequenceLength(res.fasta),
                proteinMap: res.protein_map || {},
                proteinCount: res.protein_count || 0,
                proteinCacheKey,
                sequenceKey: `protein:${fastaHash}`,
                hash: fastaHash
              };
              if (losatTiming) losatTiming.fastaPyodideFallbacks += 1;
            } else {
              try {
                entry = await extractLosatFastaFast({
                  file: sourceFile,
                  text: sourceText,
                  fmt,
                  regionSpec,
                  recordSelector,
                  reverseFlag
                });
                if (losatTiming) losatTiming.fastaJsExtractions += 1;
              } catch (fastError) {
                const res = JSON.parse(extractFirstFasta(path, fmt, regionSpec, recordSelector, reverseFlag));
                if (res.error) throw new Error(res.error);
                entry = {
                  fasta: res.fasta,
                  recordId: res.record_id || `seq_${idx + 1}`,
                  canonicalLength: Number(res.record_length || 0) || getFastaSequenceLength(res.fasta)
                };
                if (losatTiming) {
                  losatTiming.fastaPyodideFallbacks += 1;
                  console.warn('LOSAT browser FASTA extraction fell back to Pyodide:', fastError);
                }
              }
            }
            if (sourceFile && useProteinBlastp) {
              setCachedProteinExtraction(sourceFile, persistentCacheKey, entry);
            } else if (sourceFile && usePersistentFastaCache) {
              setCachedFastaExtraction(sourceFile, persistentCacheKey, entry);
            }
          }
          fastaCache.set(idx, entry);
          if (losatTiming) {
            losatTiming.fastaExtractionMs += getNow() - startedAt;
            losatTiming.totalFastaChars += entry.fasta.length;
          }
          return entry;
        };

        const writeLinearFileToFs = async (fileObj, path, { cacheText = false } = {}) => {
          return stageUploadedFile(fileObj, path, {
            cacheText,
            textCache: linearFileTextCache
          });
        };

        const getSeqHash = async (idx) => {
          if (fastaHashCache.has(idx)) return fastaHashCache.get(idx);
          const entry = await getSeqEntry(idx);
          if (entry.hash) {
            if (losatTiming) losatTiming.cacheHashHits += 1;
            if (!entry.sequenceKey) entry.sequenceKey = `seq:${entry.hash}`;
            fastaHashCache.set(idx, entry.hash);
            return entry.hash;
          }
          const startedAt = getNow();
          const hash = await hashText(entry.fasta);
          if (losatTiming) losatTiming.cacheHashMs += getNow() - startedAt;
          entry.hash = hash;
          if (!entry.sequenceKey) entry.sequenceKey = `seq:${hash}`;
          fastaHashCache.set(idx, hash);
          return hash;
        };

        const getViewTransform = async (idx) => {
          const entry = await getSeqEntry(idx);
          const length = Number(entry.canonicalLength || 0) || getFastaSequenceLength(entry.fasta);
          return {
            length,
            reverse: Boolean(viewTransformSpecs[idx]?.reverse)
          };
        };

        const buildCacheKey = async (argsKey, queryIdx, subjectIdx) => {
          const queryHash = await getSeqHash(queryIdx);
          const subjectHash = await getSeqHash(subjectIdx);
          const payload = JSON.stringify({
            cacheSchema: LOSAT_CACHE_SCHEMA,
            losatRuntimeCompatibility: losatExecutionMode === 'serial'
              ? 'serial-v1'
              : 'threaded-compatible-v1',
            losatThreadsPerJob: losatExecutionMode === 'serial'
              ? 1
              : (losatRequestedThreadsPerJob || 'auto'),
            program: losatProgram.value,
            outfmt: String(losat.outfmt || '6'),
            args: argsKey,
            queryCanonicalHash: queryHash,
            subjectCanonicalHash: subjectHash
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

        const getBlastpCandidateLimit = () => {
          if (!useProteinBlastp) return null;
          if (useOrthogroupBlastp) return 1;
          if (useCollinearBlastp) {
            return losat.blastp.collinearAnchorMode === 'all' ? null : 1;
          }
          return Math.max(1, losat.blastp.maxHits);
        };

        const buildLosatArgs = (queryIdx, subjectIdx) => {
          const args = [];
          if (losatProgram.value === 'blastn') {
            pushArg(args, '--task', losat.blastn.task);
          } else if (losatProgram.value === 'tblastx') {
            pushArg(args, '--query-gencode', getGencode(queryIdx));
            pushArg(args, '--db-gencode', getGencode(subjectIdx));
          } else {
            pushArg(args, '--max-hsps-per-subject', 1);
            pushArg(args, '--max-target-seqs', getBlastpCandidateLimit());
          }
          return args;
        };

        {
          const inputWriteStartedAt = getNow();
          for (let i = 0; i < linearSeqs.length; i++) {
            const seq = linearSeqs[i];
            if (lInputType.value === 'gb') {
              if (!seq.gb) throw new Error(`Sequence #${i + 1}: Missing GenBank file.`);
              await writeLinearFileToFs(seq.gb, `/seq_${i}.gb`, { cacheText: useLosat });
              inputArgs.push(`/seq_${i}.gb`);
            } else {
              if (!seq.gff || !seq.fasta) throw new Error(`Sequence #${i + 1}: GFF3 and FASTA are required.`);
              await writeLinearFileToFs(seq.gff, `/seq_${i}.gff`);
              await writeLinearFileToFs(seq.fasta, `/seq_${i}.fasta`, { cacheText: useLosat });
            }
          }
          if (losatTiming) losatTiming.inputWriteMs += getNow() - inputWriteStartedAt;
        }

        regionSpecs = linearSeqs.map((seq, idx) => buildRegionSpec(seq, idx));
        if (useProteinBlastp) {
          proteinRecordInstanceKeys = await buildProteinRecordInstanceKeys();
        }
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
          setProcessingStatus('Preparing LOSAT jobs...');
          const losatPairs = [];
          const losatJobs = [];
          const pendingJobKeys = new Set();
          const jobBuildStartedAt = getNow();
          const fastaExtractionBeforeJobBuild = losatTiming.fastaExtractionMs;
          const cacheHashBeforeJobBuild = losatTiming.cacheHashMs;

          const jobSpecs = [];
          if (useOrthogroupBlastp) {
            for (let i = 0; i < linearSeqs.length; i++) {
              for (let j = i + 1; j < linearSeqs.length; j++) {
                jobSpecs.push({ queryIndex: i, subjectIndex: j, pairIndex: i });
                jobSpecs.push({ queryIndex: j, subjectIndex: i, pairIndex: i });
              }
            }
          } else if (useCollinearBlastp) {
            for (let i = 0; i < linearSeqs.length - 1; i++) {
              const subjectEnd = collinearSearchScope === 'all' ? linearSeqs.length : i + 2;
              for (let j = i + 1; j < subjectEnd; j++) {
                jobSpecs.push({ queryIndex: i, subjectIndex: j, pairIndex: i });
                if (losat.blastp.collinearAnchorMode === 'rbh') {
                  jobSpecs.push({ queryIndex: j, subjectIndex: i, pairIndex: i });
                }
              }
            }
          } else {
            for (let i = 0; i < linearSeqs.length - 1; i++) {
              jobSpecs.push({ queryIndex: i, subjectIndex: i + 1, pairIndex: i });
            }
          }

          for (const spec of jobSpecs) {
            throwIfGenerationCanceled();
            const queryEntry = await getSeqEntry(spec.queryIndex);
            throwIfGenerationCanceled();
            const subjectEntry = await getSeqEntry(spec.subjectIndex);
            throwIfGenerationCanceled();
            const losatArgs = buildLosatArgs(spec.queryIndex, spec.subjectIndex);
            const cacheKey = await buildCacheKey(losatArgs, spec.queryIndex, spec.subjectIndex);
            throwIfGenerationCanceled();
            const queryCanonicalHash = await getSeqHash(spec.queryIndex);
            throwIfGenerationCanceled();
            const subjectCanonicalHash = await getSeqHash(spec.subjectIndex);
            throwIfGenerationCanceled();
            if (!queryEntry.sequenceKey || !subjectEntry.sequenceKey) {
              throw new Error('LOSAT sequence cache key was not prepared.');
            }
            sequenceEntriesByKey.set(queryEntry.sequenceKey, queryEntry.fasta);
            sequenceEntriesByKey.set(subjectEntry.sequenceKey, subjectEntry.fasta);
            const cached = cacheMap.get(cacheKey);
            const hasCachedText = isRawLosatCacheEntry(cached);
            losatTiming.totalPairs += 1;
            if (hasCachedText) losatTiming.cacheHits += 1;
            else losatTiming.cacheMisses += 1;
            const isAdjacentForwardDisplayPair = spec.subjectIndex === spec.queryIndex + 1;
            const pair = {
              pairIndex: spec.pairIndex,
              queryIndex: spec.queryIndex,
              subjectIndex: spec.subjectIndex,
              cacheKey,
              filename: buildCacheFilename(spec.pairIndex, queryEntry, subjectEntry),
              displayPair: isAdjacentForwardDisplayPair
            };
            losatPairs.push(pair);
            if (isAdjacentForwardDisplayPair) {
              cacheInfo.push({
                key: cacheKey,
                filename: pair.filename,
                display: true
              });
            }

            if (!hasCachedText && !pendingJobKeys.has(cacheKey)) {
              pendingJobKeys.add(cacheKey);
              losatJobs.push({
                pairIndex: spec.pairIndex,
                cacheKey,
                program: losatProgram.value,
                querySequenceKey: queryEntry.sequenceKey,
                subjectSequenceKey: subjectEntry.sequenceKey,
                queryCanonicalHash,
                subjectCanonicalHash,
                outfmt: losat.outfmt || '6',
                extraArgs: losatArgs
              });
            }
          }
          losatTiming.uniqueJobs = losatJobs.length;
          const jobBuildWallMs = getNow() - jobBuildStartedAt;
          const nestedFastaMs = losatTiming.fastaExtractionMs - fastaExtractionBeforeJobBuild;
          const nestedHashMs = losatTiming.cacheHashMs - cacheHashBeforeJobBuild;
          losatTiming.jobBuildWallMs += jobBuildWallMs;
          losatTiming.jobBuildMs += Math.max(0, jobBuildWallMs - nestedFastaMs - nestedHashMs);

          if (losatJobs.length > 0) {
            setProcessingStatus(`Running LOSAT: 0/${losatJobs.length} LOSAT jobs complete`);
            const runtimeWaitStartedAt = getNow();
            await waitForCancelablePromise(losatRuntimeWarmup, generationAbortSignal);
            throwIfGenerationCanceled();
            losatTiming.runtimeWaitMs += getNow() - runtimeWaitStartedAt;
            const executionStartedAt = getNow();
            const losatResults = await runLosatPairsParallel(losatJobs, {
              concurrency: getLosatParallelWorkers(),
              executionMode: losatExecutionMode,
              totalThreadBudget: losatRequestedTotalThreadBudget,
              threadsPerJob: losatRequestedThreadsPerJob,
              sequences: sequenceEntriesByKey,
              signal: generationAbortSignal,
              onRuntimeStatus: (status) => {
                losatThreadingStatus.value = status;
              },
              onProgress: ({ completed, total }) => {
                if (generationAbortSignal?.aborted || generationCancelRequested.value) return;
                setProcessingStatus(`Running LOSAT: ${completed}/${total} LOSAT jobs complete`);
              }
            });
            throwIfGenerationCanceled();
            losatTiming.executionMs += getNow() - executionStartedAt;
            losatResults.forEach((result) => {
              const job = losatJobs.find((item) => item.cacheKey === result.cacheKey);
              cacheMap.set(result.cacheKey, {
                schema: LOSAT_CACHE_SCHEMA,
                kind: 'raw-losat',
                text: result.text,
                program: losatProgram.value,
                queryCanonicalHash: job?.queryCanonicalHash || '',
                subjectCanonicalHash: job?.subjectCanonicalHash || ''
              });
            });
          } else {
            setProcessingStatus('Using cached LOSAT results...');
          }

          const blastWriteStartedAt = getNow();
          throwIfGenerationCanceled();
          if (useProteinBlastp) {
            if (!convertProteinBlast) {
              throw new Error('Current gbdraw wheel does not support LOSATP orthogroup metadata. Rebuild and redeploy the web wheel.');
            }
            const recordPayloads = [];
            for (let i = 0; i < linearSeqs.length; i += 1) {
              const entry = await getSeqEntry(i);
              recordPayloads.push({
                recordIndex: i,
                recordId: entry.recordId || `seq_${i + 1}`,
                proteinMap: entry.proteinMap || {},
                proteinCacheKey: entry.proteinCacheKey || entry.sequenceKey || `record:${i}`,
                viewTransform: await getViewTransform(i)
              });
            }
            const pairPayloads = [];
            for (const pair of losatPairs) {
              throwIfGenerationCanceled();
              const cached = cacheMap.get(pair.cacheKey);
              const losatText = isRawLosatCacheEntry(cached) ? cached.text : '';
              pairPayloads.push({
                pairIndex: pair.pairIndex,
                queryIndex: pair.queryIndex,
                subjectIndex: pair.subjectIndex,
                cacheKey: pair.cacheKey,
                blastText: losatText
              });
            }
            setProcessingStatus(
              useCollinearBlastp
                ? 'Converting LOSAT protein links to collinear ribbons...'
                : 'Converting LOSAT protein hits...'
            );
            const convertedPayload = JSON.parse(
              convertProteinBlast(
                JSON.stringify({ records: recordPayloads, pairs: pairPayloads }),
                blastpMode,
                Math.max(1, losat.blastp.maxHits),
                adv.min_bitscore,
                adv.evalue,
                adv.identity,
                adv.alignment_length,
                losat.blastp.collinearMinAnchors,
                losat.blastp.collinearMaxGeneGap,
                'cds',
                losat.blastp.collinearColorMode,
                losat.blastp.collinearAnchorMode,
                50,
                25,
                losat.blastp.collinearMaxDiagonalDrift,
                0,
                0,
                collinearSearchScope
              )
            );
            if (convertedPayload.error) throw new Error(convertedPayload.error);
            const conversionCache = convertedPayload.cache || {};
            if (conversionCache.convertedPayloadHit) losatTiming.proteinConversionCacheHits += 1;
            losatTiming.proteinFilteredHitCacheHits += Number(conversionCache.filteredHitCacheHits || 0);
            losatTiming.proteinFilteredHitCacheMisses += Number(conversionCache.filteredHitCacheMisses || 0);
            if (useOrthogroupBlastp || useCollinearBlastp) {
              setOrthogroupMetadata(convertedPayload.orthogroups || []);
            } else {
              clearOrthogroupMetadata({ clearSelection: true });
            }
            const convertedPairs = Array.isArray(convertedPayload.pairs) ? convertedPayload.pairs : [];
            for (const converted of convertedPairs) {
              const pairIndex = Number(converted?.pair_index);
              if (!Number.isInteger(pairIndex)) continue;
              const blastPath = `/blast_${pairIndex}.txt`;
              virtualBlastFiles.push({
                path: blastPath,
                text: converted.tsv || '',
                rows: Array.isArray(converted.rows) ? converted.rows : []
              });
              blastArgs.push(blastPath);
            }
          } else {
            if (!convertNucleotideBlast) {
              throw new Error('Current gbdraw wheel does not support LOSAT display coordinate conversion. Rebuild and redeploy the web wheel.');
            }
            for (const pair of losatPairs) {
              throwIfGenerationCanceled();
              const cached = cacheMap.get(pair.cacheKey);
              const blastText = isRawLosatCacheEntry(cached) ? cached.text : '';
              const converted = JSON.parse(
                convertNucleotideBlast(
                  blastText,
                  JSON.stringify(await getViewTransform(pair.queryIndex)),
                  JSON.stringify(await getViewTransform(pair.subjectIndex))
                )
              );
              if (converted.error) throw new Error(converted.error);
              const blastPath = `/blast_${pair.pairIndex}.txt`;
              virtualBlastFiles.push({
                path: blastPath,
                text: converted.tsv || '',
                rows: Array.isArray(converted.rows) ? converted.rows : []
              });
              blastArgs.push(blastPath);
            }
          }
          losatTiming.blastWriteMs += getNow() - blastWriteStartedAt;
          console.info(
            [
              `LOSAT timing: pairs=${losatTiming.totalPairs}`,
              `cache hits=${losatTiming.cacheHits}`,
              `misses=${losatTiming.cacheMisses}`,
              `unique jobs=${losatTiming.uniqueJobs}`,
              `input FS write=${formatDuration(losatTiming.inputWriteMs)}`,
              `FASTA extraction=${formatDuration(losatTiming.fastaExtractionMs)}`,
              `FASTA cache hits=${losatTiming.fastaCacheHits}`,
              `protein extraction cache hits=${losatTiming.proteinExtractionCacheHits}`,
              `JS FASTA=${losatTiming.fastaJsExtractions}`,
              `Pyodide FASTA=${losatTiming.fastaPyodideFallbacks}`,
              `protein conversion cache hits=${losatTiming.proteinConversionCacheHits}`,
              `filtered hit cache hits=${losatTiming.proteinFilteredHitCacheHits}`,
              `filtered hit cache misses=${losatTiming.proteinFilteredHitCacheMisses}`,
              `cache hashing=${formatDuration(losatTiming.cacheHashMs)}`,
              `cache hash hits=${losatTiming.cacheHashHits}`,
              `job build=${formatDuration(losatTiming.jobBuildMs)}`,
              `job build wall=${formatDuration(losatTiming.jobBuildWallMs)}`,
              `runtime wait=${formatDuration(losatTiming.runtimeWaitMs)}`,
              `execution=${formatDuration(losatTiming.executionMs)}`,
              `BLAST FS write=${formatDuration(losatTiming.blastWriteMs)}`,
              `FASTA chars=${losatTiming.totalFastaChars.toLocaleString()}`
            ].join(', ')
          );
        } else {
          for (let i = 0; i < linearSeqs.length - 1; i++) {
            const seq = linearSeqs[i];
            if (seq.blast) {
              await stageUploadedFile(seq.blast, `/blast_${i}.txt`);
              blastArgs.push(`/blast_${i}.txt`);
            }
          }
        }
        if (extractFirstFasta) {
          extractFirstFasta.destroy();
        }
        if (extractProteinFasta) {
          extractProteinFasta.destroy();
        }
        if (convertProteinBlast) {
          convertProteinBlast.destroy();
        }
        if (convertNucleotideBlast) {
          convertNucleotideBlast.destroy();
        }
        if (useLosat) {
          losatCacheInfo.value = cacheInfo;
          losatCache.value = cacheMap;
        }
        const depthEntries = linearSeqs
          .map((seq, idx) => ({ file: seq.depth, idx }))
          .filter((entry) => Boolean(entry.file));
        if (form.show_depth) {
          if (depthEntries.length === 0) {
            throw new Error('Please upload at least one Depth TSV file or disable Show depth track.');
          }
          if (depthEntries.length !== 1 && depthEntries.length !== linearSeqs.length) {
            throw new Error('Upload one Depth TSV for all records, or one Depth TSV per sequence.');
          }
          const depthPaths = [];
          if (depthEntries.length === 1 && linearSeqs.length > 1) {
            await stageUploadedFile(depthEntries[0].file, '/depth.tsv');
            depthPaths.push('/depth.tsv');
          } else {
            for (const entry of depthEntries) {
              const depthPath = `/seq_${entry.idx}.depth.tsv`;
              await stageUploadedFile(entry.file, depthPath);
              depthPaths.push(depthPath);
            }
          }
          args.push('--depth', ...depthPaths);
          args.push('--show_depth');
          appendDepthStyleArgs();
          if (
            adv.depth_height !== null &&
            adv.depth_height !== undefined &&
            adv.depth_height !== '' &&
            Number(adv.depth_height) > 0
          ) {
            args.push('--depth_height', adv.depth_height);
          }
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
      throwIfGenerationCanceled();
      setProcessingStatus('Rendering SVG...');
      const gbdrawStartedAt = getNow();
      const generationResponse = await runDiagramGeneration({
        mode: mode.value,
        args: args.map(String),
        files: Array.from(generationFileMap.values()),
        virtualBlastFiles
      });
      console.info(`gbdraw ${mode.value} wrapper execution: ${formatDuration(getNow() - gbdrawStartedAt)}.`);
      const postGbdrawTimingEntries = [];
      const res = generationResponse.results;
      if (res.error) {
        logPostGbdrawTimings(postGbdrawTimingEntries);
        if (isReflow) {
          labelReflowLastError.value = formatPythonError(res.error)?.summary || 'Auto reflow failed';
          return { status: 'error' };
        }
        errorLog.value = formatPythonError(res.error);
        return { status: 'error' };
      }

      if (generationToken !== latestGenerationToken) {
        return { status: 'stale' };
      }

      if (isReflow && requestId !== pendingReflowRequestId) {
        return { status: 'stale' };
      }

      if (isReflow) {
        skipCaptureBaseConfig.value = true;
        skipPositionReapply.value = true;
      }

      measureTiming(postGbdrawTimingEntries, 'run-analysis assign results', () => {
        results.value = shouldSuppressPairwiseIdentityLegend()
          ? stripPairwiseIdentityLegendsFromResults(res)
          : res;
      });
      logPostGbdrawTimings(postGbdrawTimingEntries);
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

      if (!isReflow) {
        extractedFeatures.value = [];
        featureRecordIds.value = [];
        selectedFeatureRecordIdx.value = 0;

        const shouldExtractFeatures =
          (mode.value === 'circular' && cInputType.value === 'gb') ||
          (mode.value === 'linear' && lInputType.value === 'gb' && linearSeqs.length > 0);

        if (shouldExtractFeatures) {
          scheduleFeatureExtraction({
            requestId: featureExtractionRequestId,
            mode: mode.value,
            cInputType: cInputType.value,
            lInputType: lInputType.value,
            circularFile: files.c_gb || null,
            linearSeqs: linearSeqs.map((seq) => ({ gb: seq.gb || null })),
            regionSpecs: regionSpecs.map((spec) => ({ file: spec?.displayFile || spec?.file || null })),
            recordSelectors: [...recordSelectors],
            reverseFlags: [...reverseFlags]
          });
        } else {
          featureExtractionPending.value = false;
          featureExtractionError.value = null;
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
      if (isDiagramGenerationCanceled(e)) {
        if (isReflow) {
          labelReflowLastError.value = null;
        } else {
          restoreManualCancelSnapshot();
          errorLog.value = null;
          processingStatus.value = 'Canceled.';
          keepProcessingStatus = true;
        }
        return { status: 'canceled' };
      }
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
        if (activeLosatAbortController === generationAbortController) {
          activeLosatAbortController = null;
        }
        if (!keepProcessingStatus) processingStatus.value = '';
        generationCancelRequested.value = false;
        processing.value = false;
      }
    }
  };

  const runAnalysis = async () => runAnalysisInternal({ runMode: 'manual' });

  const cancelRunAnalysis = () => {
    latestGenerationToken += 1;
    generationCancelRequested.value = true;
    if (activeLosatAbortController && !activeLosatAbortController.signal.aborted) {
      activeLosatAbortController.abort(new DiagramGenerationCanceledError());
    }
    if (processing.value) {
      processingStatus.value = 'Canceling generation...';
    }
    return cancelDiagramGeneration();
  };

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
      return sum + (isRawLosatCacheEntry(cached) ? cached.text.length : 0);
    }, 0);

    if (totalChars > 50 * 1024 * 1024) {
      const proceed = confirm(
        `Raw LOSAT TSV export will download about ${(totalChars / (1024 * 1024)).toFixed(1)} MB. Continue?`
      );
      if (!proceed) return;
    }

    for (let idx = 0; idx < losatCacheInfo.value.length; idx += 1) {
      const entry = losatCacheInfo.value[idx];
      const cached = cacheMap.get(entry.key);
      if (!isRawLosatCacheEntry(cached)) continue;
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
    cancelRunAnalysis,
    runLabelReflow,
    refreshCircularRecordOrder,
    downloadLosatCache,
    downloadLosatPair,
    setLosatPairFilename,
    clearLosatCache,
    getLosatPairDefaultName
  };
};
