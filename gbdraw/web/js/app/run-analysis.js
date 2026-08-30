import { prepareLosatRuntime, runLosatPairsParallel } from '../services/losat.js';
import {
  cancelDiagramGeneration,
  DIAGRAM_HELPER_OPERATIONS,
  DiagramGenerationCanceledError,
  isDiagramGenerationCanceled,
  runDiagramGeneration,
  runDiagramHelperOperation
} from '../services/diagram-generation.js';
import { buildCanonicalRenderRequest } from '../services/session-request.js';
import {
  buildLabelOverrideTsv,
  serializeLabelOverrideRows
} from './feature-editor/label-override-table.js';
import {
  applyCircularSuppressControlsToSlots,
  applyCircularTrackOrderPlacements,
  clampCircularTrackAxisIndex,
  hasEnabledCircularTrackRenderer,
  inferLegacyAxisIndexFromFeature
} from './circular-track-slots.js';
import {
  normalizeFileList,
  orderedConservationSources
} from './conservation-series.js';
import {
  applyLinearTrackOrderPlacements,
  clampLinearTrackAxisIndex,
  normalizeLinearTrackSlots,
  resolveLinearTrackAxisIndex
} from './linear-track-slots.js';
import { getDepthTrackFallbackLabel } from './depth-tracks.js';
import {
  depthFileSlotsFromValue,
  depthSlotTrackIndex,
  depthTrackCoverageCount,
  depthTrackMatrixWidth,
  isRecordMajorDepthFileMatrix,
  normalizeRecordMajorDepthFileRows,
  representativeDepthFiles,
  syncDepthSlotLabels
} from './depth-track-state.js';
import { encodeAnnotationTable } from './annotations/table-codec.js';
import {
  normalizeCollinearSearchScope
} from './losat-normalization.js';
import { buildRunInfo } from './run-info.js';
import {
  buildPairwiseLosatJobSpecs,
  resolveLinearComparisonPlan
} from './linear-comparisons.js';
import {
  buildDefaultColorOverrideTsv,
  normalizePaletteColors
} from './color-utils.js';
import { serializeSpecificRules } from './file-imports.js';
import { serializeFeatureVisibilityRules } from './feature-visibility.js';
import {
  normalizeDefinitionLineStyleState
} from './definition-line-style-state.js';
import { createZipBlob } from '../utils/zip.js';
import { classifyOptionalPositiveNumber } from '../utils/optional-positive-number.js';
import { cloneJsonData, cloneJsonValue } from '../services/json-clone.js';
import { downloadBlob, downloadTextFile } from '../services/text-download.js';
import {
  normalizeCircularPlotTitlePosition,
  normalizeLinearPlotTitlePosition
} from './plot-title-position.js';
import {
  requireCurrentCircularMultiRecordSizeMode,
  requireCurrentCollinearAnchorMode,
  requireCurrentCollinearColorMode,
  requireCurrentCollinearMaxConflicts,
  requireCurrentCollinearMaxDiagonalDrift,
  requireCurrentCollinearMaxParalogLinks,
  requireCurrentCollinearMaxUnitGap,
  requireCurrentCollinearMergeOrientation,
  requireCurrentCollinearMinAnchors,
  requireCurrentCollinearSearchScope,
  requireCurrentCollinearUnitMode,
  requireCurrentLinearLabelPlacement,
  requireCurrentLinearTrackLayout,
  requireCurrentOrthogroupMemberMaxHits,
  requireCurrentOrthogroupMembershipMode,
  requireCurrentProteinBlastpCandidateLimit,
  requireCurrentProteinBlastpMaxHits,
  requireCurrentProteinBlastpMode
} from './current-option-values.js';
import {
  discoverGffFastaRecords,
  discoverSequenceRecords
} from './record-discovery.js';
import {
  LOSAT_DERIVED_CACHE_SCHEMA,
  NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
  PROTEIN_LOSAT_CACHE_SCHEMA,
  classifyRawLosatCacheEntry,
  emptyProteinIdentityManifest,
  getCurrentRawLosatCacheEntry,
  isCurrentRawLosatCacheEntry,
  isLosatDerivedCacheEntry,
  mergeProteinIdentityManifests,
  normalizeLosatArgs,
  sameLosatArgs,
  transitionLegacyProteinCandidate,
  validateDerivedProteinReferences,
  validateProteinIdentityManifest
} from './losat-cache.js';
import { comparisonFiltersForMode } from '../mode-profiles.js';
import { normalizeUserFacingError } from '../services/error-normalization.js';
import {
  cloneFileBytesForTransfer,
  readFileBytes,
  readFileText
} from '../services/file-content-cache.js';
import {
  validateFeatureCatalog
} from '../services/feature-catalog.js';
import {
  buildOrthogroupFeatureIndex,
  enrichFeaturesWithOrthogroups
} from '../services/orthogroup-feature-metadata.js';
import {
  prepareCandidateRenderCommit,
  prepareReflowResultCommit
} from './candidate-render.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from '../services/runtime-test-hooks.js';
import {
  IMPORTED_COMPARISON_DISPOSITIONS,
  createImportedComparisonIntentState,
  inheritCommittedComparisonIntent
} from '../services/imported-comparison-intent.js';

const DEFAULT_CIRCULAR_CONSERVATION_BLAST_FILTERS = Object.freeze(
  comparisonFiltersForMode('circular')
);
const DEFAULT_LINEAR_BLAST_FILTERS = Object.freeze(
  comparisonFiltersForMode('linear')
);

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
const LOSAT_DERIVED_CACHE_LIMIT = 16;
export const LOSAT_EXPORT_CONFIRM_THRESHOLD_BYTES = 50 * 1024 * 1024;

export const totalHydratedLosatExportBytes = (hydratedResults) => (
  (Array.isArray(hydratedResults) ? hydratedResults : []).reduce((sum, result) => {
    const reportedBytes = Number(result?.utf8Bytes);
    return sum + (
      Number.isFinite(reportedBytes) && reportedBytes >= 0
        ? reportedBytes
        : new TextEncoder().encode(String(result?.text || '')).byteLength
    );
  }, 0)
);

export const losatExportBytesExceedLimit = (
  totalBytes,
  limitBytes = LOSAT_EXPORT_CONFIRM_THRESHOLD_BYTES
) => Number(totalBytes) > Number(limitBytes);

export const confirmHydratedLosatExport = (
  totalBytes,
  confirmDownload = globalThis.confirm,
  limitBytes = LOSAT_EXPORT_CONFIRM_THRESHOLD_BYTES
) => (
  !losatExportBytesExceedLimit(totalBytes, limitBytes) ||
  typeof confirmDownload !== 'function' ||
  Boolean(confirmDownload(
    `Raw LOSAT TSV export will download about ${(totalBytes / (1024 * 1024)).toFixed(1)} MB. Continue?`
  ))
);

const buildLosatCachePayload = ({
  identityKind,
  flow,
  program,
  outfmt,
  args,
  queryCanonicalHash,
  subjectCanonicalHash,
  queryProteinSetHash,
  subjectProteinSetHash,
  queryRuntimeBindingHash,
  subjectRuntimeBindingHash,
  queryRecordInstanceKey,
  subjectRecordInstanceKey
}) => {
  if (identityKind === 'protein') {
    return {
      cacheSchema: PROTEIN_LOSAT_CACHE_SCHEMA,
      identityKind: 'protein',
      program: String(program || 'blastp'),
      outfmt: String(outfmt || '6'),
      args: normalizeLosatArgs(args),
      idEncoding: 'runtime-handle-v1',
      queryProteinSetHash: String(queryProteinSetHash || ''),
      subjectProteinSetHash: String(subjectProteinSetHash || ''),
      queryRuntimeBindingHash: String(queryRuntimeBindingHash || ''),
      subjectRuntimeBindingHash: String(subjectRuntimeBindingHash || ''),
      queryRecordInstanceKey: String(queryRecordInstanceKey || ''),
      subjectRecordInstanceKey: String(subjectRecordInstanceKey || '')
    };
  }
  const payload = {
    cacheSchema: NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    program,
    outfmt: String(outfmt || '6'),
    args: normalizeLosatArgs(args),
    queryCanonicalHash,
    subjectCanonicalHash
  };
  if (flow) payload.flow = flow;
  return payload;
};

const getRawLosatCacheEntry = (cacheMap, cacheKey, metadata, manifest = null) => {
  if (!cacheMap) return null;
  const direct = getCurrentRawLosatCacheEntry(cacheMap, cacheKey, metadata, manifest);
  if (direct) return direct;
  if (metadata?.identityKind === 'protein') return null;

  const expectedProgram = String(metadata?.program || '');
  const expectedOutfmt = String(metadata?.outfmt || '6');
  const expectedFlow = String(metadata?.flow || '');
  const expectedQuery = String(metadata?.queryCanonicalHash || '');
  const expectedSubject = String(metadata?.subjectCanonicalHash || '');

  for (const [key, entry] of cacheMap.entries()) {
    if (classifyRawLosatCacheEntry(entry) !== 'nucleotide-current') continue;
    if (String(entry.queryCanonicalHash || '') !== expectedQuery) continue;
    if (String(entry.subjectCanonicalHash || '') !== expectedSubject) continue;
    const entryProgram = String(entry.program || expectedProgram);
    if (entryProgram && expectedProgram && entryProgram !== expectedProgram) continue;
    if (entry.outfmt && String(entry.outfmt) !== expectedOutfmt) continue;
    if (entry.flow && String(entry.flow) !== expectedFlow) continue;
    if (Array.isArray(entry.args) && !sameLosatArgs(entry.args, metadata?.args)) continue;
    return { key, entry };
  }

  return null;
};

const promoteRawLosatCacheEntry = (cacheMap, cacheKey, found, metadata) => {
  if (!cacheMap || !found?.entry) return found?.entry || null;
  if (classifyRawLosatCacheEntry(found.entry) !== 'nucleotide-current') {
    return found.entry;
  }
  const promoted = {
    ...found.entry,
    program: metadata.program || found.entry.program || '',
    outfmt: String(metadata.outfmt || found.entry.outfmt || '6'),
    args: normalizeLosatArgs(metadata.args || found.entry.args),
    queryCanonicalHash: metadata.queryCanonicalHash || found.entry.queryCanonicalHash || '',
    subjectCanonicalHash: metadata.subjectCanonicalHash || found.entry.subjectCanonicalHash || ''
  };
  if (metadata.flow) promoted.flow = metadata.flow;
  if (found.key && found.key !== cacheKey) cacheMap.delete(found.key);
  cacheMap.set(cacheKey, promoted);
  return promoted;
};

const pruneLosatDerivedCache = (cacheMap) => {
  if (!cacheMap || typeof cacheMap.delete !== 'function') return;
  while (cacheMap.size > LOSAT_DERIVED_CACHE_LIMIT) {
    const oldestKey = cacheMap.keys().next().value;
    if (oldestKey === undefined) break;
    cacheMap.delete(oldestKey);
  }
};

export const buildLosatDerivedPayloadCachePayload = ({
  mode,
  maxHits,
  bitscore,
  evalue,
  identity,
  alignmentLength,
  collinearMinAnchors,
  collinearMaxUnitGap,
  collinearUnitMode,
  collinearColorMode,
  collinearAnchorMode,
  collinearMergeOrientation,
  collinearMaxDiagonalDrift,
  collinearMaxConflictsInMergeGap,
  collinearMaxParalogLinksPerOrthogroup,
  collinearSearchScope,
  orthogroupMembershipMode,
  orthogroupMemberMaxHits,
  recordPayloads,
  pairPayloads
}) => {
  const normalizedMode = String(mode || 'pairwise');
  const payload = {
    cacheSchema: LOSAT_DERIVED_CACHE_SCHEMA,
    semantics: 'derived-option-conformance-v1',
    idEncoding: 'runtime-handle-v1',
    converter: 'convert_losatp_blastp_pairs_to_genomic_payload',
    featureIdentity: 'stable-source-rendered-display-v1',
    mode: normalizedMode,
    thresholds: {
      bitscore: String(bitscore),
      evalue: String(evalue),
      identity: String(identity),
      alignmentLength: String(alignmentLength)
    },
    records: (Array.isArray(recordPayloads) ? recordPayloads : [])
      .map((record) => ({
        recordIndex: Number(record?.recordIndex),
        proteinCacheKey: String(record?.proteinCacheKey || ''),
        runtimeBindingHash: String(record?.runtimeBindingHash || ''),
        displayBindingHash: String(record?.displayBindingHash || ''),
        viewTransform: {
          length: Number(record?.viewTransform?.length || 0),
          reverse: Boolean(record?.viewTransform?.reverse)
        }
      }))
      .sort((left, right) => left.recordIndex - right.recordIndex),
    pairs: (Array.isArray(pairPayloads) ? pairPayloads : [])
      .map((pair) => ({
        pairIndex: Number(pair?.pairIndex),
        queryIndex: Number(pair?.queryIndex),
        subjectIndex: Number(pair?.subjectIndex),
        cacheKey: String(pair?.cacheKey || '')
      }))
  };
  if (normalizedMode === 'pairwise') {
    payload.pairwise = { maxHits: Number(maxHits) || 5 };
  }
  if (['orthogroup', 'collinear'].includes(normalizedMode)) {
    payload.orthogroup = {
      membershipMode: String(orthogroupMembershipMode || 'anchor_core_v1'),
      memberMaxHits: Number(orthogroupMemberMaxHits) || 5
    };
  }
  if (normalizedMode === 'collinear') {
    payload.collinear = {
      minAnchors: String(collinearMinAnchors),
      // Persisted LOSAT derived-cache schema 3 uses this historical identity key.
      maxGeneGap: String(collinearMaxUnitGap),
      unitMode: String(collinearUnitMode || 'auto'),
      colorMode: String(collinearColorMode || 'orientation'),
      anchorMode: String(collinearAnchorMode || 'rbh'),
      mergeOrientation: String(collinearMergeOrientation || 'either'),
      maxDiagonalDrift: String(collinearMaxDiagonalDrift),
      maxConflictsInMergeGap: String(collinearMaxConflictsInMergeGap),
      maxParalogLinksPerOrthogroup: String(collinearMaxParalogLinksPerOrthogroup),
      searchScope: String(collinearSearchScope || 'adjacent')
    };
  }
  return payload;
};

const getLosatDerivedCacheEntry = (cacheMap, key, manifest) => {
  if (!cacheMap || !key) return null;
  const entry = cacheMap.get(key);
  if (
    !isLosatDerivedCacheEntry(entry, { allowLegacy: false }) ||
    !validateDerivedProteinReferences(entry, manifest) ||
    !hasRequiredCanonicalAnalysisResource(entry.mode, entry.payload)
  ) return null;
  cacheMap.delete(key);
  cacheMap.set(key, entry);
  return entry.payload;
};

export const hasRequiredCanonicalAnalysisResource = (mode, payload) => {
  const normalizedMode = String(mode || '').trim().toLowerCase();
  if (!['orthogroup', 'collinear'].includes(normalizedMode)) return true;
  const resource = normalizedMode === 'collinear'
    ? payload?.collinearityResult
    : payload?.orthogroupResult;
  const expectedKind = normalizedMode === 'collinear' ? 'result' : 'orthogroupResult';
  const expectedType = normalizedMode === 'collinear' ? 'CollinearityResult' : 'OrthogroupResult';
  return Boolean(
    resource
    && typeof resource === 'object'
    && !Array.isArray(resource)
    && [1, 2].includes(resource.schema)
    && resource.kind === expectedKind
    && resource.value?.type === expectedType
    && resource.value.fields
    && typeof resource.value.fields === 'object'
    && !Array.isArray(resource.value.fields)
  );
};

export const resolveProteinBlastpCandidateLimit = (candidateLimit) => (
  requireCurrentProteinBlastpCandidateLimit(candidateLimit)
);

const inferredResolvedProteinMode = (comparisons) => {
  const values = Array.isArray(comparisons) ? comparisons : [];
  if (values.some((comparison) => comparison?.kind === 'collinearityResult')) {
    return 'collinear';
  }
  if (values.some((comparison) => comparison?.kind === 'orthogroupResult')) {
    return 'orthogroup';
  }
  if (values.some((comparison) => comparison?.kind === 'precomputedProteinComparison')) {
    return 'pairwise';
  }
  return '';
};

const sameNumber = (left, right) => Number(left) === Number(right);

const canReuseResolvedProteinArtifacts = ({
  canonicalComparisons,
  committedRequest,
  active
}) => {
  const persisted = Array.isArray(canonicalComparisons) ? canonicalComparisons : [];
  const committed = Array.isArray(committedRequest?.comparisons)
    ? committedRequest.comparisons
    : [];
  const persistedMarker = persisted.find((comparison) => (
    comparison?.kind === 'generatedProteinComparison' && comparison.mode === 'none'
  ));
  const committedMarker = committed.find((comparison) => (
    comparison?.kind === 'generatedProteinComparison' && comparison.mode === 'none'
  ));
  if (!persistedMarker || !committedMarker || !active) return false;

  const persistedMode = inferredResolvedProteinMode(persisted);
  const committedMode = inferredResolvedProteinMode(committed);
  if (!persistedMode || persistedMode !== committedMode || committedMode !== active.mode) {
    return false;
  }

  const settings = committedMarker.settings || {};
  const filters = committedRequest.diagramOptions || {};
  if (
    !sameNumber(filters.bitscore ?? DEFAULT_LINEAR_BLAST_FILTERS.bitscore, active.bitscore)
    || !sameNumber(filters.evalue ?? DEFAULT_LINEAR_BLAST_FILTERS.evalue, active.evalue)
    || !sameNumber(filters.identity ?? DEFAULT_LINEAR_BLAST_FILTERS.identity, active.identity)
    || !sameNumber(
      filters.alignmentLength ?? DEFAULT_LINEAR_BLAST_FILTERS.alignment_length,
      active.alignmentLength
    )
    || (settings.proteinBlastpCandidateLimit ?? null) !== active.candidateLimit
  ) {
    return false;
  }

  if (active.mode === 'pairwise') {
    return sameNumber(settings.proteinBlastpMaxHits ?? 5, active.maxHits);
  }
  if (
    String(settings.orthogroupMembershipMode || 'anchor_core_v1')
      !== active.orthogroupMembershipMode
    || !sameNumber(settings.orthogroupMemberMaxHits ?? 5, active.memberMaxHits)
  ) {
    return false;
  }
  if (active.mode === 'orthogroup') return true;

  const parameters = settings.collinearityParams?.parameters || {};
  return (
    sameNumber(parameters.minAnchors ?? 1, active.minAnchors)
    && sameNumber(parameters.maxUnitGap ?? 0, active.maxUnitGap)
    && sameNumber(parameters.maxDiagonalDrift ?? 0, active.maxDiagonalDrift)
    && sameNumber(parameters.maxConflicts ?? 1, active.maxConflicts)
    && String(parameters.mergeOrientation || 'either') === active.mergeOrientation
    && String(settings.collinearityUnitMode || 'auto') === active.unitMode
    && String(settings.collinearityAnchorMode || 'rbh') === active.anchorMode
    && String(settings.collinearitySearchScope || 'adjacent') === active.searchScope
    && String(settings.collinearityColorMode || 'orientation') === active.colorMode
    && sameNumber(
      settings.collinearMaxParalogLinksPerOrthogroup ?? 2,
      active.maxParalogLinks
    )
  );
};

const stripRuntimeCacheStats = (payload) => {
  const cloned = cloneJsonData(payload);
  if (cloned && typeof cloned === 'object' && !Array.isArray(cloned)) {
    delete cloned.cache;
  }
  return cloned;
};

const LEGACY_PROTEIN_REFERENCE_RE = /p_[A-Za-z0-9._%+-]+?_\d+_\d+_(?:-1|0|1)_[0-9a-f]{12}(?:_[2-9][0-9]*)?/g;

const collectLegacyProteinReferences = (...values) => {
  const references = new Set();
  const visited = new WeakSet();
  const visit = (value) => {
    if (typeof value === 'string') {
      for (const match of value.matchAll(LEGACY_PROTEIN_REFERENCE_RE)) {
        references.add(match[0]);
      }
      return;
    }
    if (Array.isArray(value)) {
      value.forEach(visit);
      return;
    }
    if (!value || typeof value !== 'object' || visited.has(value)) return;
    visited.add(value);
    Object.entries(value).forEach(([key, item]) => {
      visit(key);
      visit(item);
    });
  };
  values.forEach(visit);
  return Array.from(references).sort();
};

const rewriteMappedProteinReferences = (value, idMap) => {
  if (!idMap || typeof idMap !== 'object' || Array.isArray(idMap)) return cloneJsonData(value);
  if (typeof value === 'string') {
    if (Object.prototype.hasOwnProperty.call(idMap, value)) return String(idMap[value]);
    return value.replace(
      LEGACY_PROTEIN_REFERENCE_RE,
      (reference) => (
        Object.prototype.hasOwnProperty.call(idMap, reference)
          ? String(idMap[reference])
          : reference
      )
    );
  }
  if (Array.isArray(value)) {
    return value.map((item) => rewriteMappedProteinReferences(item, idMap));
  }
  if (!value || typeof value !== 'object') return value;
  const rewritten = {};
  Object.entries(value).forEach(([key, item]) => {
    const nextKey = Object.prototype.hasOwnProperty.call(idMap, key) ? String(idMap[key]) : key;
    if (Object.prototype.hasOwnProperty.call(rewritten, nextKey)) {
      throw new Error(`Protein reference migration produced duplicate key '${nextKey}'.`);
    }
    rewritten[nextKey] = rewriteMappedProteinReferences(item, idMap);
  });
  return rewritten;
};

const setLosatDerivedCacheEntry = (cacheMap, key, { mode, payload, manifest }) => {
  if (!cacheMap || !key || !payload || typeof payload !== 'object' || Array.isArray(payload)) {
    return null;
  }
  const entry = {
    schema: LOSAT_DERIVED_CACHE_SCHEMA,
    kind: 'derived-losatp-payload',
    idEncoding: 'runtime-handle-v1',
    key,
    mode: String(mode || ''),
    payload: stripRuntimeCacheStats(payload)
  };
  if (
    !isLosatDerivedCacheEntry(entry, { allowLegacy: false }) ||
    !validateDerivedProteinReferences(entry, manifest) ||
    !hasRequiredCanonicalAnalysisResource(entry.mode, entry.payload)
  ) return null;
  if (cacheMap.has(key)) cacheMap.delete(key);
  cacheMap.set(key, entry);
  pruneLosatDerivedCache(cacheMap);
  return entry.payload;
};

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
  const sourceText = typeof text === 'string' ? text : await readFileText(file);
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
  const sourceText = typeof text === 'string' ? text : await readFileText(file);
  const records = fmt === 'genbank' ? parseGenbankRecordsFast(sourceText) : parseFastaRecordsFast(sourceText);
  if (!records.length) throw new Error('No records found for circular conservation reference.');
  return {
    fasta: records.map((record) => buildFastaText(record)).join(''),
    recordIds: records.map((record) => record.id),
    canonicalLength: records.reduce((sum, record) => sum + String(record.sequence || '').length, 0)
  };
};
const buildConservationSeries = (sourceFiles, circularConservation) => {
  return orderedConservationSources(sourceFiles, circularConservation).map((entry) => ({
    label: entry.label,
    color: entry.color
  }));
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
export const createRunAnalysis = ({
  state,
  serializeCanonicalFiles,
  canonicalSessionVersion,
  adoptCanonicalRenderArtifacts,
  getCommittedCanonicalRenderRequest = null,
  getCommittedCanonicalSession = null,
  captureGeneratedArtifactHandle,
  restoreGeneratedArtifactHandle,
  setGeneratedArtifactIdentity = null,
  resetPreviewViewport,
  validateAnnotationTargets = null,
  prepareLinearRecordCatalog = null,
  losatExecutor = runLosatPairsParallel,
  prepareCandidateCommit = prepareCandidateRenderCommit,
  prepareReflowCommit = prepareReflowResultCommit
}) => {
  const {
    processing,
    processingStatus,
    generationCancelRequested,
    results,
    selectedResultIndex,
    failedGeneratePreservedResult,
    resultGenerationKey,
    resultPanelTab,
    lastRunInfo,
    trackSlotResolvedGeometry,
    errorLog,
    semanticFileWatchersSuppressed,
    sessionImportRollbackInProgress,
    zoom,
    canvasPan,
    skipCaptureBaseConfig,
    skipPositionReapply,
    matchSequenceRegistry,
    featureColorOverrides,
    featureVisibilityRules,
    featureStrokeOverrides,
    legendColorOverrides,
    legendStrokeOverrides,
    selectedPalette,
    currentColors,
    paletteDefinitions,
    appliedPaletteName,
    appliedPaletteColors,
    pendingPaletteName,
    pendingPaletteColors,
    filterMode,
    manualSpecificRules,
    manualWhitelist,
    manualPriorityRules,
    form,
    adv,
    mode,
    cInputType,
    lInputType,
    losatProgram,
    losat,
    losatCacheInfo,
    losatThreadingStatus,
    losatCache,
    losatDerivedCache,
    proteinIdentityManifest,
    legacyProteinRawCandidates,
    legacyProteinDerivedEvidence,
    circularConservation,
    annotationSets,
    orthogroups,
    collinearGroups,
    featureOrthogroupIndex,
    selectedOrthogroupAlignmentFeature,
    orthogroupNameOverrides,
    orthogroupDescriptionOverrides,
    selectedOrthogroupId,
    circularRecordList,
    circularRecordDiscovery,
    files,
    linearSeqs,
    linearRecordLayoutEnabled,
    linearRecordRows,
    linearComparisonPlan,
    linearComparisonResolution,
    importedComparisonIntent,
    generatedLegendPosition,
    generatedMode,
    generatedMultiRecordCanvas,
    generatedCircularPlotTitlePosition,
    shouldDeferCircularPreviewUpdates,
    extractedFeatures,
    biologicalFeatures,
    featureCatalog,
    featureSelectorSafetyScope,
    featureEditorStatus,
    featureExtractionPending,
    featureExtractionError,
    featureRecordIds,
    selectedFeatureRecordIdx,
    editableLabels,
    labelTextScopeDialog,
    labelTextFeatureOverrides,
    canonicalLabelOverrideRows,
    labelTextBulkOverrides,
    labelTextFeatureOverrideSources,
    labelVisibilityOverrides,
    labelOverrideBuildWarning,
    labelReflowProcessing,
    labelReflowLastError
  } = state;
  if (
    typeof captureGeneratedArtifactHandle !== 'function'
    || typeof restoreGeneratedArtifactHandle !== 'function'
  ) {
    throw new Error('createRunAnalysis requires generated artifact handle handlers.');
  }
  const executeLosatJobs = (...args) => {
    const override = globalThis.__GBDRAW_LOSAT_EXECUTOR__;
    return (typeof override === 'function' ? override : losatExecutor)(...args);
  };
  let pendingReflowRequestId = 0;
  let activeReflowRequestId = 0;
  let pendingReflowReason = 'label-edit';
  let featureExtractionRequestId = 0;
  let latestGenerationToken = 0;
  let circularRecordRefreshGeneration = 0;
  let activeCircularRecordRefresh = null;
  let activeLosatAbortController = null;
  let latestCliHelperFiles = [];
  let latestCliHelperArchiveName = 'out-cli-files.zip';
  const cloneCliHelperFiles = (files) => (
    (Array.isArray(files) ? files : []).map((file) => ({ ...file }))
  );
  const captureGeneratedArtifactRuntimeState = () => {
    const helperFiles = cloneCliHelperFiles(latestCliHelperFiles);
    return {
      latestCliHelperFiles: helperFiles,
      latestCliHelperArchiveName,
      losatTelemetry: cloneJsonValue(globalThis.__GBDRAW_LAST_LOSAT_TELEMETRY__, null),
      retainedBytes: helperFiles.reduce((total, file) => (
        total
        + String(file?.name || '').length * 2
        + String(file?.data || '').length * 2
        + Math.max(0, Number(file?.retainedBytes) || 0)
        + 256
      ), String(latestCliHelperArchiveName || '').length * 2 + 16_384)
    };
  };
  const restoreGeneratedArtifactRuntimeState = (runtimeState = {}, { ui = {} } = {}) => {
    latestCliHelperFiles = cloneCliHelperFiles(runtimeState?.latestCliHelperFiles);
    latestCliHelperArchiveName = String(
      runtimeState?.latestCliHelperArchiveName || 'out-cli-files.zip'
    );
    globalThis.__GBDRAW_LAST_LOSAT_TELEMETRY__ = cloneJsonValue(
      runtimeState?.losatTelemetry,
      null
    );
    if (typeof resetPreviewViewport === 'function') {
      resetPreviewViewport({ pan: ui?.canvasPan });
    }
  };
  const recordDiscoverySuppressed = () => Boolean(
    semanticFileWatchersSuppressed?.value ||
    sessionImportRollbackInProgress?.value
  );

  const buildLatestCliHelperFiles = (runInfo, generatedCliFileMap, archiveBaseName) => {
    const helperFiles = Array.isArray(runInfo?.helperFiles) ? runInfo.helperFiles : [];
    if (helperFiles.length === 0) {
      return {
        files: [],
        archiveName: 'out-cli-files.zip'
      };
    }

    const bySlot = new Map();
    generatedCliFileMap.forEach((entry) => {
      const slot = String(entry?.slot || '').trim();
      if (slot && !bySlot.has(slot)) bySlot.set(slot, entry);
    });

    const files = helperFiles
      .map((helper) => {
        const path = String(helper?.path || '').trim();
        const slot = String(helper?.slot || '').trim();
        const entry = generatedCliFileMap.get(path) || bySlot.get(slot);
        if (!entry) return null;
        return {
          name: String(helper?.name || entry.name || 'helper.tsv'),
          retainedBytes: Math.max(0, Number(entry.retainedBytes) || 0),
          ...(typeof entry.buildData === 'function'
            ? { buildData: entry.buildData }
            : { data: entry.data })
        };
      })
      .filter(Boolean);
    const archiveStem = makeSafeFilename(`${archiveBaseName || form.prefix || 'out'}-cli-files`);
    return {
      files,
      archiveName: `${archiveStem}.zip`
    };
  };

  const downloadCliHelperFiles = () => {
    if (!latestCliHelperFiles.length) {
      alert('No CLI helper files are available for the latest run.');
      return;
    }
    const materializedFiles = latestCliHelperFiles.map((file) => ({
      name: file.name,
      data: typeof file.buildData === 'function' ? file.buildData() : file.data
    }));
    const totalChars = materializedFiles.reduce(
      (sum, file) => sum + String(file.data ?? '').length,
      0
    );
    if (totalChars > 50 * 1024 * 1024) {
      const proceed = confirm(
        `CLI helper file export will download about ${(totalChars / (1024 * 1024)).toFixed(1)} MB. Continue?`
      );
      if (!proceed) return;
    }
    downloadBlob(createZipBlob(materializedFiles), latestCliHelperArchiveName);
  };

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
    const circularTrackSlotError = extractCircularTrackSlotError(err);
    return normalizeUserFacingError(
      circularTrackSlotError
        ? { type: err?.type || 'ValidationError', message: circularTrackSlotError, notes: err?.notes }
        : err
    );
  };

  const formatJsError = (err) => normalizeUserFacingError(err);

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
    const index = buildOrthogroupFeatureIndex(groups);
    const groupIds = groups
      .map((group) => String(group?.id || '').trim())
      .filter(Boolean);
    orthogroups.value = groups;
    featureOrthogroupIndex.value = index;
    extractedFeatures.value = enrichFeaturesWithOrthogroups(extractedFeatures.value, index);
    if (biologicalFeatures) {
      biologicalFeatures.value = enrichFeaturesWithOrthogroups(
        biologicalFeatures.value,
        index
      );
    }
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

  const setFeatureEditorStatus = (updates = {}) => {
    if (!featureEditorStatus || typeof featureEditorStatus !== 'object') return;
    Object.assign(featureEditorStatus, {
      status: updates.status ?? featureEditorStatus.status,
      generationId: updates.generationId ?? featureEditorStatus.generationId,
      error: updates.error === undefined ? featureEditorStatus.error : updates.error,
      summaryCount: updates.summaryCount ?? featureEditorStatus.summaryCount,
      detailsCacheSize: updates.detailsCacheSize ?? featureEditorStatus.detailsCacheSize
    });
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

  const getResolvedLinearEdge = (edgeKey) => {
    const normalizedKey = String(edgeKey || '').trim();
    if (!normalizedKey) return null;
    const edges = linearComparisonResolution?.value?.edges;
    return (Array.isArray(edges) ? edges : []).find((edge) => edge.edgeKey === normalizedKey) || null;
  };

  const getLosatCacheInfoEntry = (edgeKey) => {
    const entries = Array.isArray(losatCacheInfo.value) ? losatCacheInfo.value : [];
    if (Number.isInteger(edgeKey)) return entries[edgeKey] || null;
    const normalizedKey = String(edgeKey || '').trim();
    return entries.find((entry) => String(entry?.edgeKey || '') === normalizedKey) || null;
  };

  const getLosatPairDefaultName = (edgeKey, queryEntry = null, subjectEntry = null) => {
    const cacheEntry = getLosatCacheInfoEntry(edgeKey);
    const edge = getResolvedLinearEdge(edgeKey) || cacheEntry;
    const queryIndex = Number(edge?.queryIndex);
    const subjectIndex = Number(edge?.subjectIndex);
    const leftLabel = getSeqLabel(
      linearSeqs[queryIndex],
      queryEntry?.recordId || `seq_${Number.isInteger(queryIndex) ? queryIndex + 1 : 1}`
    );
    const rightLabel = getSeqLabel(
      linearSeqs[subjectIndex],
      subjectEntry?.recordId || `seq_${Number.isInteger(subjectIndex) ? subjectIndex + 1 : 2}`
    );
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

  const hydrateLosatDownloadText = async (cacheKey, cached) => {
    if (classifyRawLosatCacheEntry(cached) !== 'protein-current') {
      return {
        text: String(cached?.text || ''),
        utf8Bytes: new TextEncoder().encode(String(cached?.text || '')).byteLength
      };
    }
    const response = await runDiagramHelperOperation(
      DIAGRAM_HELPER_OPERATIONS.HYDRATE_PROTEIN_LOSAT_TSV,
      {
        entry: cloneJsonData({ ...cached, key: String(cacheKey || '') }),
        identityManifest: cloneJsonData(proteinIdentityManifest.value)
      }
    );
    const result = response.result;
    if (result.status !== 'ok' || typeof result.text !== 'string') {
      throw new Error(
        result.error ||
        'Protein raw TSV export contains an unresolved internal reference.'
      );
    }
    return {
      text: result.text,
      utf8Bytes: Number(result.utf8Bytes) || new TextEncoder().encode(result.text).byteLength
    };
  };

  const downloadLosatPair = async (edgeKey, customName) => {
    const entry = getLosatCacheInfoEntry(edgeKey);
    const cacheMap = losatCache.value;
    if (!entry || !cacheMap) return;
    const cached = cacheMap.get(entry.key);
    if (!isCurrentRawLosatCacheEntry(cached)) return;
    const defaultName = getLosatPairDefaultName(entry.edgeKey || edgeKey);
    const fallbackOrdinal = Number.isInteger(Number(entry.ordinal))
      ? Number(entry.ordinal)
      : 0;
    const filename = normalizeLosatFilename(
      customName,
      entry.filename || defaultName || `losat_pair_${fallbackOrdinal + 1}.tsv`
    );
    entry.filename = filename;
    const hydrated = await hydrateLosatDownloadText(entry.key, cached);
    downloadTextFile(
      filename || 'losat.tsv',
      hydrated.text,
      'text/tab-separated-values'
    );
  };

  const setLosatPairFilename = (edgeKey, customName) => {
    const entry = getLosatCacheInfoEntry(edgeKey);
    if (!entry) return;
    const defaultName = getLosatPairDefaultName(entry.edgeKey || edgeKey);
    const fallbackOrdinal = Number.isInteger(Number(entry.ordinal))
      ? Number(entry.ordinal)
      : 0;
    entry.filename = normalizeLosatFilename(
      customName,
      entry.filename || defaultName || `losat_pair_${fallbackOrdinal + 1}.tsv`
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

  const circularDiscoveryTargetsCurrentInput = () => {
    const inputType = cInputType.value;
    const primaryFile = inputType === 'gff' ? files.c_gff : files.c_gb;
    const pairedFile = inputType === 'gff' ? files.c_fasta : null;
    return (
      circularRecordDiscovery.inputType === inputType &&
      circularRecordDiscovery.primaryFile === primaryFile &&
      circularRecordDiscovery.pairedFile === pairedFile
    );
  };

  const circularDiscoveryMatchesCurrentInput = () => (
    circularRecordDiscovery.status === 'ready' &&
    circularDiscoveryTargetsCurrentInput()
  );

  const validateDepthInputPresence = () => {
    const customDepthRequested = mode.value === 'linear'
      ? (
          adv.linear_track_slots_enabled === true &&
          (Array.isArray(adv.linear_track_slots) ? adv.linear_track_slots : []).some((slot) => (
            slot?.enabled !== false && String(slot?.renderer || '') === 'depth'
          ))
        )
      : (
          adv.circular_track_slots_enabled === true &&
          (Array.isArray(adv.circular_track_slots) ? adv.circular_track_slots : []).some((slot) => (
            slot?.enabled !== false && String(slot?.renderer || '') === 'depth'
          ))
    );
    if (!form.show_depth && !customDepthRequested) return '';
    if (mode.value === 'circular') {
      const discoveredCount = (
        circularDiscoveryMatchesCurrentInput() &&
        Array.isArray(circularRecordList.value)
      )
        ? circularRecordList.value.length
        : 0;
      const recordCount = discoveredCount > 0
        ? discoveredCount
        : (
            isRecordMajorDepthFileMatrix(files.c_depth)
              ? Math.max(1, files.c_depth.length)
              : 1
          );
      const rows = normalizeRecordMajorDepthFileRows(files.c_depth, recordCount);
      if (isRecordMajorDepthFileMatrix(files.c_depth) && files.c_depth.length !== recordCount) {
        return `Circular Depth matrix has ${files.c_depth.length} record rows; expected ${recordCount}.`;
      }
      if (!rows.some((row) => row.some(Boolean))) {
        return 'Please upload a Depth TSV file or disable Show depth track.';
      }
      const logicalWidth = depthTrackMatrixWidth(rows);
      for (let trackIndex = 0; trackIndex < logicalWidth; trackIndex += 1) {
        if (depthTrackCoverageCount(rows, trackIndex) > 0) continue;
        return `Depth series #${trackIndex + 1} (logical track index ${trackIndex}) has no TSV source in any record. Add a TSV or remove the series.`;
      }
      return '';
    }

    const depthFileSlots = (value) => (Array.isArray(value) ? value.slice() : (value ? [value] : []));
    const perRecordDepthFiles = linearSeqs.map((seq) => depthFileSlots(seq.depth));
    const depthCount = perRecordDepthFiles.reduce((sum, items) => sum + items.filter(Boolean).length, 0);
    if (depthCount === 0) {
      return 'Please upload at least one Depth TSV file or disable Show depth track.';
    }
    const logicalWidth = depthTrackMatrixWidth(perRecordDepthFiles);
    for (let trackIndex = 0; trackIndex < logicalWidth; trackIndex += 1) {
      if (depthTrackCoverageCount(perRecordDepthFiles, trackIndex) > 0) continue;
      return `Depth series #${trackIndex + 1} (logical track index ${trackIndex}) has no TSV source in any record. Add a TSV or remove the series.`;
    }
    return '';
  };

  const shouldSuppressPairwiseIdentityLegend = (comparisonPlanSnapshot = null) => {
    return (
      mode.value === 'linear' &&
      comparisonPlanSnapshot?.hasLosatIntent === true &&
      losatProgram.value === 'blastp' &&
      normalizeBlastpMode(losat.blastp?.mode) === 'collinear' &&
      normalizeCollinearColorMode(losat.blastp?.collinearColorMode) === 'orientation'
    );
  };

  const runCircularRecordRefresh = async ({ suppress = false } = {}) => {
    const refreshGeneration = ++circularRecordRefreshGeneration;
    if (suppress || recordDiscoverySuppressed()) return;
    if (!Array.isArray(adv.multi_record_positions)) {
      adv.multi_record_positions = [];
    }
    const inputType = cInputType.value;
    const primaryFile = inputType === 'gff' ? files.c_gff : files.c_gb;
    const pairedFile = inputType === 'gff' ? files.c_fasta : null;
    const hasCompleteInput = Boolean(primaryFile && (inputType !== 'gff' || pairedFile));
    const hasActiveInput = mode.value === 'circular' && hasCompleteInput;
    Object.assign(circularRecordDiscovery, {
      status: hasActiveInput ? 'loading' : 'idle',
      error: '',
      inputType,
      primaryFile: primaryFile || null,
      pairedFile: pairedFile || null
    });
    if (
      !hasActiveInput
    ) {
      circularRecordList.value = [];
      adv.multi_record_positions.splice(0, adv.multi_record_positions.length);
      return;
    }
    if (
      refreshGeneration !== circularRecordRefreshGeneration ||
      recordDiscoverySuppressed() ||
      mode.value !== 'circular' ||
      cInputType.value !== inputType ||
      (inputType === 'gff' ? files.c_gff : files.c_gb) !== primaryFile ||
      (inputType === 'gff' ? files.c_fasta : null) !== pairedFile
    ) return;

    try {
      const records = inputType === 'gff'
        ? await discoverGffFastaRecords({
            gffFile: primaryFile,
            fastaFile: pairedFile,
            readText: readFileText
          })
        : await discoverSequenceRecords({
            file: primaryFile,
            format: 'genbank',
            readText: readFileText
          });
      if (
        refreshGeneration !== circularRecordRefreshGeneration ||
        recordDiscoverySuppressed() ||
        mode.value !== 'circular' ||
        cInputType.value !== inputType ||
        (inputType === 'gff' ? files.c_gff : files.c_gb) !== primaryFile ||
        (inputType === 'gff' ? files.c_fasta : null) !== pairedFile
      ) return;
      const nextRecords = records.map((entry) => ({
        selector: entry.selector,
        record_id: entry.recordId,
        record_length: entry.recordLength
      }));
      circularRecordList.value = nextRecords;
      circularRecordDiscovery.status = 'ready';
      const nextPositions = mergeCircularRecordPositions(nextRecords, adv.multi_record_positions);
      adv.multi_record_positions.splice(0, adv.multi_record_positions.length, ...nextPositions);
    } catch (error) {
      if (
        refreshGeneration !== circularRecordRefreshGeneration ||
        recordDiscoverySuppressed() ||
        mode.value !== 'circular' ||
        cInputType.value !== inputType ||
        (inputType === 'gff' ? files.c_gff : files.c_gb) !== primaryFile ||
        (inputType === 'gff' ? files.c_fasta : null) !== pairedFile
      ) return;
      console.warn('Failed to refresh circular record order:', error);
      circularRecordList.value = [];
      circularRecordDiscovery.status = 'error';
      circularRecordDiscovery.error = 'Could not read records from the circular input file(s).';
      adv.multi_record_positions.splice(0, adv.multi_record_positions.length);
    }
  };

  const refreshCircularRecordOrder = (options = {}) => {
    const inputType = cInputType.value;
    const fingerprint = [
      Boolean(options.suppress || recordDiscoverySuppressed()),
      mode.value,
      inputType,
      inputType === 'gff' ? files.c_gff : files.c_gb,
      inputType === 'gff' ? files.c_fasta : null
    ];
    if (
      activeCircularRecordRefresh &&
      fingerprint.length === activeCircularRecordRefresh.fingerprint.length &&
      fingerprint.every((value, index) => (
        Object.is(value, activeCircularRecordRefresh.fingerprint[index])
      ))
    ) {
      return activeCircularRecordRefresh.promise;
    }
    const entry = { fingerprint, promise: null };
    entry.promise = runCircularRecordRefresh(options).finally(() => {
      if (activeCircularRecordRefresh === entry) activeCircularRecordRefresh = null;
    });
    activeCircularRecordRefresh = entry;
    return entry.promise;
  };

  const runAnalysisInternal = async ({
    runMode = 'manual',
    requestId = 0,
    comparisonPlanSnapshot = null,
    generatedArtifactHandle = null,
    comparisonExecution = null
  } = {}) => {
    const isReflow = runMode === 'reflow';
    if (!isReflow) {
      recordSessionLifecycleEvent('generate-start');
      recordSessionLifecycleEvent('generation-input-resolution-start');
    }
    const useCommittedComparison = comparisonExecution?.mode === 'inherit';
    const forceEmptyComparison = useCommittedComparison || comparisonExecution?.mode === 'clear';
    const activeComparisonPlanSnapshot = mode.value === 'linear'
      ? (
          forceEmptyComparison
            ? resolveLinearComparisonPlan({
                plan: { mode: 'none', defaultSource: 'losat', edges: [] },
                sequences: linearSeqs,
                layout: linearRecordLayoutEnabled.value ? linearRecordRows : [],
                losatProgram: losatProgram.value,
                blastpMode: normalizeBlastpMode(losat.blastp?.mode)
              })
            : comparisonPlanSnapshot || resolveLinearComparisonPlan({
                plan: linearComparisonPlan,
                sequences: linearSeqs,
                layout: linearRecordLayoutEnabled.value ? linearRecordRows : [],
                losatProgram: losatProgram.value,
                blastpMode: normalizeBlastpMode(losat.blastp?.mode)
              })
        )
      : null;
    let linearRecordCatalog = null;

    const generationToken = ++latestGenerationToken;
    let keepProcessingStatus = false;
    let canceledAttemptOwnsPresentation = false;
    let generationAbortController = null;
    let generationAbortSignal = null;
    const committedArtifactHandle = isReflow
      ? null
      : (generatedArtifactHandle || await captureGeneratedArtifactHandle());
    const restoreCommittedArtifact = async () => {
      if (!committedArtifactHandle) return;
      await restoreGeneratedArtifactHandle(committedArtifactHandle);
    };
    const finishCanceledManualRun = async () => {
      if (latestGenerationToken !== generationToken + 1) {
        return { status: 'stale' };
      }
      canceledAttemptOwnsPresentation = true;
      await restoreCommittedArtifact();
      errorLog.value = null;
      processingStatus.value = 'Canceled.';
      generationCancelRequested.value = false;
      processing.value = false;
      keepProcessingStatus = true;
      return { status: 'canceled' };
    };
    const setProcessingStatus = (message) => {
      if (!isReflow && generationToken === latestGenerationToken) {
        processingStatus.value = String(message || '');
      }
    };
    const throwIfGenerationCanceled = () => {
      if (!isReflow && generationCancelRequested.value) {
        throw new DiagramGenerationCanceledError();
      }
    };
    if (isReflow && mode.value === 'circular' && shouldDeferCircularPreviewUpdates.value) {
      return { status: 'skipped' };
    }
    generationAbortController =
      !isReflow && typeof AbortController === 'function' ? new AbortController() : null;
    generationAbortSignal = generationAbortController?.signal || null;
    if (!isReflow) activeLosatAbortController = generationAbortController;
    if (mode.value === 'linear' && typeof prepareLinearRecordCatalog === 'function') {
      if (!isReflow) {
        processing.value = true;
        processingStatus.value = 'Reading input records...';
      }
      let prepared;
      try {
        prepared = await prepareLinearRecordCatalog(
          activeComparisonPlanSnapshot?.hasComparisonIntent
        );
      } catch (error) {
        prepared = {
          catalog: null,
          error: error?.message || 'Could not read records from the Linear input file(s).'
        };
      } finally {
        if (!isReflow && generationToken === latestGenerationToken) {
          processing.value = false;
          processingStatus.value = '';
        }
      }
      if (!isReflow && generationToken !== latestGenerationToken) {
        if (activeLosatAbortController === generationAbortController) {
          activeLosatAbortController = null;
        }
        return generationAbortSignal?.aborted
          ? finishCanceledManualRun()
          : { status: 'stale' };
      }
      if (prepared?.error) {
        const message = String(prepared.error);
        if (isReflow) labelReflowLastError.value = message;
        else {
          await restoreCommittedArtifact();
          errorLog.value = formatJsError(new Error(message));
        }
        if (activeLosatAbortController === generationAbortController) {
          activeLosatAbortController = null;
        }
        return { status: 'error' };
      }
      linearRecordCatalog = prepared?.catalog || null;
    }
    if (!isReflow && mode.value === 'circular') {
      const inputType = cInputType.value;
      const primaryFile = inputType === 'gff' ? files.c_gff : files.c_gb;
      const pairedFile = inputType === 'gff' ? files.c_fasta : null;
      const hasCompleteInput = Boolean(primaryFile && (inputType !== 'gff' || pairedFile));
      if (
        hasCompleteInput &&
        !circularDiscoveryMatchesCurrentInput()
      ) {
        processing.value = true;
        processingStatus.value = 'Reading input records...';
        try {
          await refreshCircularRecordOrder();
        } finally {
          if (generationToken === latestGenerationToken) {
            processing.value = false;
            processingStatus.value = '';
          }
        }
        if (generationToken !== latestGenerationToken) {
          if (activeLosatAbortController === generationAbortController) {
            activeLosatAbortController = null;
          }
          return generationAbortSignal?.aborted
            ? finishCanceledManualRun()
            : { status: 'stale' };
        }
        if (!circularDiscoveryMatchesCurrentInput()) {
          const message = circularRecordDiscovery.error || 'Could not read records from the circular input file(s).';
          await restoreCommittedArtifact();
          errorLog.value = formatJsError(new Error(message));
          if (activeLosatAbortController === generationAbortController) {
            activeLosatAbortController = null;
          }
          return { status: 'error' };
        }
      }
    }
    const depthInputError = validateDepthInputPresence();
    if (depthInputError) {
      if (isReflow) {
        labelReflowLastError.value = depthInputError;
      } else {
        await restoreCommittedArtifact();
        errorLog.value = formatJsError(new Error(depthInputError));
      }
      if (activeLosatAbortController === generationAbortController) {
        activeLosatAbortController = null;
      }
      return { status: 'error' };
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
    const manualRunStartedAt = isReflow ? null : getNow();
    const manualRunStartedAtIso = isReflow ? null : new Date().toISOString();
    let structuredLosatTelemetry = null;
    if (!isReflow) globalThis.__GBDRAW_LAST_LOSAT_TELEMETRY__ = null;
    const legacyPromotionTransaction = [];
    let legacyPromotionCommitted = false;
    let commitProteinMigration = null;
    let pendingLosatCacheCommit = null;

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
      setFeatureEditorStatus({
        status: 'idle',
        generationId: featureExtractionRequestId,
        error: null,
        summaryCount: 0,
        detailsCacheSize: 0
      });
      processing.value = true;
      processingStatus.value = 'Preparing input files...';
      resultPanelTab.value = 'preview';
      errorLog.value = null;
      skipCaptureBaseConfig.value = false;
      skipPositionReapply.value = false;
      resetLabelScopeDialogState();
      window._origPairwiseMin = activeRunColors.pairwise_match_min || '#FFE7E7';
      window._origPairwiseMax = activeRunColors.pairwise_match_max || '#FF7272';
    }
    labelOverrideBuildWarning.value = '';

    try {
      if (mode.value === 'linear') {
        if (!activeComparisonPlanSnapshot || !Array.isArray(activeComparisonPlanSnapshot.edges)) {
          throw new Error('A resolved Linear comparison plan is required.');
        }
        if (activeComparisonPlanSnapshot.error) {
          throw new Error(activeComparisonPlanSnapshot.error);
        }
      }
      const annotationLoadComparison = (
        mode.value === 'linear' && activeComparisonPlanSnapshot?.hasComparisonIntent === true
      );
      if (typeof validateAnnotationTargets === 'function' && annotationSets.length > 0) {
        const annotationError = validateAnnotationTargets({
          loadComparison: annotationLoadComparison
        });
        if (annotationError) throw new Error(annotationError);
      }
      let regionSpecs = [];
      let recordSelectors = [];
      let reverseFlags = [];
      const resolvedComparisons = [];
      let resolvedCircularConservation = [];
      const runInfoFileMap = new Map();
      const generatedCliFileMap = new Map();
      const textEncoder = new TextEncoder();
      const textDecoder = new TextDecoder();
      const getPayloadName = (path, fallback = 'input') => {
        const name = String(path || '').split('/').filter(Boolean).pop();
        return name || fallback;
      };
      const generatedSlotForPath = (path) => {
        const normalizedPath = String(path || '').trim();
        if (normalizedPath === '/combined_d.tsv') return 'generatedFiles.combined_d';
        if (normalizedPath === '/combined_t.tsv') return 'generatedFiles.combined_t';
        if (normalizedPath === '/manual_wl.tsv') return 'generatedFiles.manual_wl';
        if (normalizedPath === '/priority.tsv') return 'generatedFiles.priority';
        if (normalizedPath === '/web_label_table.tsv') return 'generatedFiles.web_label_table';
        if (normalizedPath === '/web_feature_visibility_table.tsv') return 'generatedFiles.web_feature_visibility_table';
        if (normalizedPath === '/web_feature_table.tsv') return 'generatedFiles.web_feature_visibility_table';
        if (normalizedPath === '/web_annotations.tsv') return 'generatedFiles.web_annotations';
        const conservationMatch = normalizedPath.match(/^\/conservation_blast_(\d+)\.txt$/);
        if (conservationMatch) return `generatedFiles.circular_conservation_blasts[${Number(conservationMatch[1])}]`;
        const blastMatch = normalizedPath.match(/^\/blast_(\d+)\.txt$/);
        if (blastMatch) return `generatedFiles.losat_blasts[${Number(blastMatch[1])}]`;
        return normalizedPath ? `generatedFiles.${getPayloadName(normalizedPath).replace(/[^\w]+/g, '_')}` : '';
      };
      const registerRunInfoFile = (path, { name = '', slot = '', kind = 'uploaded' } = {}) => {
        const normalizedPath = String(path || '').trim();
        if (!normalizedPath) return;
        const displayName = String(name || getPayloadName(normalizedPath)).trim() || getPayloadName(normalizedPath);
        runInfoFileMap.set(normalizedPath, {
          path: normalizedPath,
          name: displayName,
          slot: String(slot || '').trim(),
          kind: kind === 'generated' ? 'generated' : 'uploaded'
        });
      };
      const recordGeneratedCliFile = (path, data, { name = '', slot = '' } = {}) => {
        const normalizedPath = String(path || '').trim();
        if (!normalizedPath) return;
        const displayName = String(name || getPayloadName(normalizedPath)).trim() || getPayloadName(normalizedPath);
        generatedCliFileMap.set(normalizedPath, {
          path: normalizedPath,
          name: displayName,
          slot: String(slot || '').trim(),
          data: String(data ?? '')
        });
      };
      const recordDeferredGeneratedCliFile = (
        path,
        buildData,
        { name = '', slot = '', retainedBytes = 0 } = {}
      ) => {
        const normalizedPath = String(path || '').trim();
        if (!normalizedPath) return;
        if (typeof buildData !== 'function') {
          throw new TypeError('A deferred CLI helper file requires a builder.');
        }
        const displayName = String(name || getPayloadName(normalizedPath)).trim()
          || getPayloadName(normalizedPath);
        generatedCliFileMap.set(normalizedPath, {
          path: normalizedPath,
          name: displayName,
          slot: String(slot || '').trim(),
          retainedBytes: Math.max(0, Number(retainedBytes) || 0),
          buildData
        });
      };
      const stageTextFile = (
        path,
        text,
        { name = '', slot = '' } = {}
      ) => {
        throwIfGenerationCanceled();
        const displayName = name || getPayloadName(path);
        const resolvedSlot = slot || generatedSlotForPath(path);
        registerRunInfoFile(path, {
          name: displayName,
          slot: resolvedSlot,
          kind: 'generated'
        });
        recordGeneratedCliFile(path, text, {
          name: displayName,
          slot: resolvedSlot
        });
      };
      const stageUploadedFile = async (fileObj, path, {
        cacheText = false,
        textCache = null,
        slot = ''
      } = {}) => {
        if (!fileObj) return false;
        throwIfGenerationCanceled();
        const bytes = await readFileBytes(fileObj);
        throwIfGenerationCanceled();
        if (cacheText && textCache) {
          textCache.set(fileObj, textDecoder.decode(bytes));
        }
        registerRunInfoFile(path, {
          name: fileObj.name || getPayloadName(path),
          slot,
          kind: 'uploaded'
        });
        return true;
      };

      const normalizedOutputPrefix = String(form.prefix || '').trim();

      const activePaletteName = String(
        isReflow ? appliedPaletteName.value : (selectedPalette?.value || appliedPaletteName.value || 'default')
      ).trim() || 'default';
      const paletteBaseColors = normalizePaletteColors(
        paletteDefinitions.value?.[activePaletteName] ||
        paletteDefinitions.value?.default ||
        {}
      );
      const dContent = buildDefaultColorOverrideTsv({
        colors: activeRunColors,
        paletteColors: paletteBaseColors
      });
      if (dContent.trim() !== '') {
        stageTextFile('/combined_d.tsv', `${dContent}\n`);
      }

      const tContent = serializeSpecificRules(manualSpecificRules);
      if (tContent.trim() !== '') {
        stageTextFile('/combined_t.tsv', tContent);
      }

      if (Array.isArray(annotationSets) && annotationSets.length > 0) {
        stageTextFile('/web_annotations.tsv', encodeAnnotationTable(annotationSets), {
          name: 'annotations.tsv',
          slot: 'generatedFiles.web_annotations'
        });
      }

      if (filterMode.value === 'Whitelist') {
        if (manualWhitelist.length > 0) {
          let wlContent = '';
          manualWhitelist.forEach((r) => {
            if (r.feat && r.qual) wlContent += `${r.feat}\t${r.qual}\t${r.key}\n`;
          });
          stageTextFile('/manual_wl.tsv', wlContent);
        }
      }

      let pContent = '';
      manualPriorityRules.forEach((r) => {
        pContent += `${r.feat}\t${r.order}\n`;
      });
      if (pContent.trim() !== '') {
        stageTextFile('/priority.tsv', pContent);
      }

      const labelOverride = buildLabelOverrideTsv(labelTextFeatureOverrides, labelTextBulkOverrides, {
        editableLabels: editableLabelsSnapshot,
        extractedFeatures: extractedFeatures.value,
        featureOverrideSources: featureOverrideSourcesSnapshot,
        visibilityOverrides: visibilityOverridesSnapshot
      });
      const effectiveLabelOverrideTsv = serializeLabelOverrideRows(canonicalLabelOverrideRows.value) || labelOverride.tsv;
      if (labelOverride.skippedMissingSourceCount > 0) {
        labelOverrideBuildWarning.value = `${labelOverride.skippedMissingSourceCount} feature override row(s) were skipped due to missing source label context.`;
      }
      if (effectiveLabelOverrideTsv) {
        stageTextFile('/web_label_table.tsv', effectiveLabelOverrideTsv);
      }
      let featureVisibilityTablePath = null;
      let featureVisibilityCacheKey = '';
      const featureVisibilityTsv = serializeFeatureVisibilityRules(featureVisibilityRules?.value || []);
      if (featureVisibilityTsv.trim()) {
        featureVisibilityTablePath = '/web_feature_visibility_table.tsv';
        featureVisibilityCacheKey = featureVisibilityTsv;
        stageTextFile(featureVisibilityTablePath, featureVisibilityCacheKey);
      }
      const validateDepthStyleSettings = () => {
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
        if (
          adv.depth_window_size !== null &&
          adv.depth_window_size !== undefined &&
          adv.depth_window_size !== ''
        ) {
          if (Number(adv.depth_window_size) <= 0) throw new Error('Depth window must be greater than 0.');
        }
        if (
          adv.depth_step_size !== null &&
          adv.depth_step_size !== undefined &&
          adv.depth_step_size !== ''
        ) {
          if (Number(adv.depth_step_size) <= 0) throw new Error('Depth step must be greater than 0.');
        }
        if (
          adv.depth_large_tick_interval !== null &&
          adv.depth_large_tick_interval !== undefined &&
          adv.depth_large_tick_interval !== ''
        ) {
          if (Number(adv.depth_large_tick_interval) <= 0) throw new Error('Depth large tick interval must be greater than 0.');
        }
        if (
          mode.value === 'circular' &&
          adv.depth_small_tick_interval !== null &&
          adv.depth_small_tick_interval !== undefined &&
          adv.depth_small_tick_interval !== ''
        ) {
          if (Number(adv.depth_small_tick_interval) <= 0) throw new Error('Depth small tick interval must be greater than 0.');
        }
        if (
          adv.depth_tick_font_size !== null &&
          adv.depth_tick_font_size !== undefined &&
          adv.depth_tick_font_size !== ''
        ) {
          if (Number(adv.depth_tick_font_size) <= 0) throw new Error('Depth tick font size must be greater than 0.');
        }
      };
      const normalizeDepthTrackValue = (value) => {
        if (value === null || value === undefined || value === '') return null;
        const numeric = Number(value);
        return Number.isFinite(numeric) && numeric > 0 ? String(numeric) : null;
      };
      const depthTrackEntriesFromRows = (rows) => Array.from(
        { length: depthTrackMatrixWidth(rows) },
        (_, index) => {
          const filesForRecords = rows.map((row) => row[index] || null);
          const file = filesForRecords.find(Boolean) || null;
          return {
            file,
            files: filesForRecords,
            index
          };
        }
      );
      const ensureDepthTrackConfigAt = (index) => {
        const idx = Math.max(0, Number(index) || 0);
        if (!Array.isArray(adv.depth_tracks)) adv.depth_tracks = [];
        while (adv.depth_tracks.length <= idx) {
          const nextIndex = adv.depth_tracks.length;
          adv.depth_tracks.push({
            label: getDepthTrackFallbackLabel(nextIndex),
            color: nextIndex === 0 ? String(adv.depth_color || '#4A90E2') : '',
            height: null,
            large_tick_interval: null,
            small_tick_interval: null,
            tick_font_size: null
          });
        }
        const config = adv.depth_tracks[idx];
        if (!config || typeof config !== 'object' || Array.isArray(config)) {
          adv.depth_tracks[idx] = {
            label: getDepthTrackFallbackLabel(idx),
            color: idx === 0 ? String(adv.depth_color || '#4A90E2') : '',
            height: null,
            large_tick_interval: null,
            small_tick_interval: null,
            tick_font_size: null
          };
        }
        return adv.depth_tracks[idx];
      };
      const syncDepthSlotLegendLabelsFromTrackConfigs = (slots, trackFiles = []) => {
        if (!Array.isArray(slots)) return;
        const trackCount = Array.isArray(trackFiles)
          ? trackFiles.length
          : 0;
        for (let trackIndex = 0; trackIndex < trackCount; trackIndex += 1) {
          ensureDepthTrackConfigAt(trackIndex);
        }
        syncDepthSlotLabels({
          slots,
          depthTracks: adv.depth_tracks,
          activeCount: trackCount
        });
      };
      const validateEnabledDepthSlots = (slots, activeCount, kind) => {
        const count = Math.max(0, Number(activeCount) || 0);
        let fallbackIndex = 0;
        (Array.isArray(slots) ? slots : []).forEach((slot) => {
          if (!slot || slot.enabled === false || String(slot.renderer || '') !== 'depth') return;
          const trackIndex = depthSlotTrackIndex(slot, fallbackIndex);
          fallbackIndex += 1;
          if (trackIndex !== null && trackIndex >= 0 && trackIndex < count) return;
          const slotId = String(slot.id || 'depth').trim() || 'depth';
          const indexText = trackIndex === null ? 'unknown' : String(trackIndex);
          throw new Error(
            `${kind} depth slot '${slotId}' references removed depth track index ${indexText}. Select an existing Depth TSV or remove the slot.`
          );
        });
      };
      const linearDepthRepresentativeFiles = () => {
        const rows = linearSeqs.map((seq) => depthFileSlotsFromValue(seq.depth));
        const maxDepthTracks = depthTrackMatrixWidth(rows);
        return Array.from({ length: maxDepthTracks }, (_, trackIndex) => (
          rows.map((row) => row[trackIndex]).find(Boolean) || null
        ));
      };
      const syncLinearDepthSlotHeightFromTrackConfig = (slot) => {
        if (!slot || String(slot.renderer || '') !== 'depth') return;
        const rawTrackIndex = Number(slot.params?.track_index);
        const trackIndex = Number.isInteger(rawTrackIndex) && rawTrackIndex >= 0 ? rawTrackIndex : 0;
        const config = ensureDepthTrackConfigAt(trackIndex);
        const configuredHeight = normalizeDepthTrackValue(config.height);
        if (configuredHeight) {
          slot.height = configuredHeight;
          return;
        }
        const slotHeight = normalizeDepthTrackValue(slot.height);
        if (slotHeight) {
          config.height = Number(slotHeight);
        }
      };
      const normalizeGcContentPercentState = () => {
        adv.gc_content_mode = String(adv.gc_content_mode || '').trim().toLowerCase() === 'percent'
          ? 'percent'
          : 'deviation';
        if (adv.gc_content_mode !== 'percent') return;

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

        if (
          adv.gc_content_tick_interval !== null &&
          adv.gc_content_tick_interval !== undefined &&
          adv.gc_content_tick_interval !== ''
        ) {
          if (Number(adv.gc_content_tick_interval) <= 0) throw new Error('GC content large tick interval must be greater than 0.');
        }
        if (
          adv.gc_content_small_tick_interval !== null &&
          adv.gc_content_small_tick_interval !== undefined &&
          adv.gc_content_small_tick_interval !== ''
        ) {
          if (Number(adv.gc_content_small_tick_interval) <= 0) throw new Error('GC content small tick interval must be greater than 0.');
        }
        if (
          adv.gc_content_tick_font_size !== null &&
          adv.gc_content_tick_font_size !== undefined &&
          adv.gc_content_tick_font_size !== ''
        ) {
          if (Number(adv.gc_content_tick_font_size) <= 0) throw new Error('GC content tick font size must be greater than 0.');
        }
      };
      normalizeGcContentPercentState();

      if (mode.value === 'circular') {
        const normalizedCircularPlotTitle = String(form.plot_title || '').trim();
        const normalizedPlotTitlePosition = normalizeCircularPlotTitlePosition(adv.plot_title_position);
        const useCircularTrackSlots = adv.circular_track_slots_enabled === true;
        const circularTrackAxisIndex = clampCircularTrackAxisIndex(
          adv.circular_track_slots_axis_index,
          Array.isArray(adv.circular_track_slots) ? adv.circular_track_slots.length : 0
        );
        const circularTrackSlots = useCircularTrackSlots
          ? applyCircularSuppressControlsToSlots(
              applyCircularTrackOrderPlacements(
                adv.circular_track_slots,
                adv.nt,
                form.track_type,
                circularTrackAxisIndex
              ),
              form
            )
          : [];
        if (useCircularTrackSlots) {
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

        const labelsModeRaw =
          typeof form.labels_mode === 'string'
            ? form.labels_mode
            : (form.allow_inner_labels ? 'both' : (form.show_labels ? 'out' : 'none'));
        const labelsMode = String(labelsModeRaw || 'none').trim().toLowerCase();
        const normalizedLabelRendering = labelsMode === 'none' ? 'auto' : normalizeLabelRendering(adv.label_rendering);
        adv.label_rendering = normalizedLabelRendering;
        const normalizedCircularLabelPlacement =
          String(adv.circular_label_placement || '').trim().toLowerCase() === 'radial'
            ? 'radial'
            : 'horizontal';
        adv.circular_label_placement = normalizedCircularLabelPlacement;
        if (form.multi_record_canvas) {
          const effectiveRecordPositions = mergeCircularRecordPositions(
            circularRecordList.value,
            adv.multi_record_positions
          );
          const normalizedSizeMode = requireCurrentCircularMultiRecordSizeMode(
            adv.multi_record_size_mode
          );
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
        }

        const normalizedCircularLabelSpacing = normalizePositiveNumberOrNull(adv.circular_label_spacing);
        adv.circular_label_spacing = normalizedCircularLabelSpacing;
        const discoveredCircularRecordCount = (
          circularDiscoveryMatchesCurrentInput() &&
          Array.isArray(circularRecordList.value)
        )
          ? circularRecordList.value.length
          : 0;
        const circularDepthRecordCount = discoveredCircularRecordCount > 0
          ? discoveredCircularRecordCount
          : (
              isRecordMajorDepthFileMatrix(files.c_depth)
                ? Math.max(1, files.c_depth.length)
                : 1
            );
        const circularDepthRows = normalizeRecordMajorDepthFileRows(
          files.c_depth,
          circularDepthRecordCount
        );
        const circularDepthEntries = depthTrackEntriesFromRows(circularDepthRows);
        if (useCircularTrackSlots) {
          let circularDepthSlotOrdinal = 0;
          circularTrackSlots.forEach((slot) => {
            if (slot?.enabled === false || String(slot?.renderer || '') !== 'depth') return;
            const parsedTrackIndex = Number(slot?.params?.track_index);
            const trackIndex = Number.isInteger(parsedTrackIndex) && parsedTrackIndex >= 0
              ? parsedTrackIndex
              : circularDepthSlotOrdinal;
            if (!slot.params || !Number.isInteger(parsedTrackIndex) || parsedTrackIndex < 0) {
              slot.params = { ...(slot.params || {}), track_index: trackIndex };
            }
            circularDepthSlotOrdinal += 1;
          });
          syncDepthSlotLegendLabelsFromTrackConfigs(
            circularTrackSlots,
            representativeDepthFiles(circularDepthRows)
          );
          validateEnabledDepthSlots(circularTrackSlots, circularDepthEntries.length, 'Circular');
        }
        const hasCircularDepthFile = circularDepthEntries.length > 0;
        const circularSlotNeedsDepth = useCircularTrackSlots && hasEnabledCircularTrackRenderer(circularTrackSlots, 'depth');
        if (form.show_depth || circularSlotNeedsDepth) {
          if (!hasCircularDepthFile) throw new Error('Please upload a Depth TSV file or disable Show depth track.');
          for (let depthIndex = 0; depthIndex < circularDepthEntries.length; depthIndex += 1) {
            const entry = circularDepthEntries[depthIndex];
            if (!entry.file) {
              throw new Error(
                `Depth series #${depthIndex + 1} (logical track index ${depthIndex}) has no TSV source in any record.`
              );
            }
          }
          validateDepthStyleSettings();
        }

        if (cInputType.value === 'gb') {
          if (!files.c_gb) throw new Error('Please upload a GenBank file.');
        } else {
          if (!files.c_gff || !files.c_fasta) throw new Error('GFF3 and FASTA are required.');
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
          circularConservation.ring_width = normalizePositiveNumberOrNull(
            circularConservation.ring_width
          );
          circularConservation.ring_gap = normalizePositiveNumberOrNull(
            circularConservation.ring_gap
          );
          const normalizeConservationState = () => {
            adv.min_bitscore = normalizeBlastThresholdNumber(
              adv.min_bitscore,
              DEFAULT_CIRCULAR_CONSERVATION_BLAST_FILTERS.bitscore
            );
            adv.evalue = normalizeBlastThresholdText(
              adv.evalue,
              DEFAULT_CIRCULAR_CONSERVATION_BLAST_FILTERS.evalue
            );
            adv.identity = normalizeBlastThresholdNumber(
              adv.identity,
              DEFAULT_CIRCULAR_CONSERVATION_BLAST_FILTERS.identity
            );
            adv.alignment_length = normalizeBlastThresholdNumber(
              adv.alignment_length,
              DEFAULT_CIRCULAR_CONSERVATION_BLAST_FILTERS.alignment_length,
              { integer: true }
            );
          };

          const runCircularLosatConservation = async (comparisonEntries) => {
            const circularLosatProgram = normalizeCircularConservationLosatProgram(
              circularConservation.losat_program
            );
            circularConservation.losat_program = circularLosatProgram;
            const subjectGencode = normalizePositiveInteger(circularConservation.subject_gencode, 1);
            circularConservation.subject_gencode = subjectGencode;
            const buildExtraArgs = (comparisonGencode) => {
              if (circularLosatProgram === 'tblastx') {
                return [
                  '--query-gencode',
                  String(normalizePositiveInteger(comparisonGencode, 1)),
                  '--db-gencode',
                  String(subjectGencode)
                ];
              }
              const normalizedTask = String(losat.blastn?.task || 'megablast').trim() || 'megablast';
              return ['--task', normalizedTask];
            };
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
            const cacheMap = new Map(losatCache.value || []);
            const cacheInfo = [];
            const losatPairs = [];
            const losatJobs = [];
            const pendingJobKeys = new Set();
            const executionMode = getLosatExecutionMode();

            for (let index = 0; index < comparisonEntries.length; index += 1) {
              throwIfGenerationCanceled();
              const comparisonEntry = comparisonEntries[index];
              const fileObj = comparisonEntry?.file || comparisonEntry;
              const queryFasta = await readFileText(fileObj);
              if (getFastaSequenceLength(queryFasta) <= 0) {
                throw new Error(`Pairwise comparison FASTA #${index + 1} has no sequence data.`);
              }
              const queryHash = await hashText(queryFasta);
              const querySequenceKey = `circular-query:${queryHash}`;
              sequenceEntriesByKey.set(querySequenceKey, queryFasta);
              const extraArgs = buildExtraArgs(comparisonEntry?.losat_gencode);
              const cacheMetadata = {
                flow: 'circular-conservation',
                program: circularLosatProgram,
                outfmt: String(losat.outfmt || '6'),
                args: extraArgs,
                queryCanonicalHash: queryHash,
                subjectCanonicalHash: subjectHash
              };
              const cacheKey = await hashText(JSON.stringify(buildLosatCachePayload(cacheMetadata)));
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
              const cached = getRawLosatCacheEntry(cacheMap, cacheKey, cacheMetadata);
              const hasCachedText = Boolean(cached);
              if (cached) promoteRawLosatCacheEntry(cacheMap, cacheKey, cached, cacheMetadata);
              if (!hasCachedText && !pendingJobKeys.has(cacheKey)) {
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
              const losatResults = await executeLosatJobs(losatJobs, {
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
                  schema: NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
                  kind: 'raw-losat',
                  identityKind: 'nucleotide',
                  text: result.text,
                  program: circularLosatProgram,
                  flow: 'circular-conservation',
                  outfmt: String(losat.outfmt || '6'),
                  args: job?.extraArgs || [],
                  queryCanonicalHash: job?.queryCanonicalHash || '',
                  subjectCanonicalHash: job?.subjectCanonicalHash || ''
                });
              });
            } else {
              setProcessingStatus('Using cached LOSAT conservation results...');
            }

            const resolved = [];
            for (const pair of losatPairs) {
              const cached = cacheMap.get(pair.cacheKey);
              const blastText = isCurrentRawLosatCacheEntry(cached) ? cached.text : '';
              const blastPath = `/conservation_blast_${pair.sourceIndex}.txt`;
              stageTextFile(blastPath, blastText, {
                name: pair.filename,
                slot: `generatedFiles.circular_conservation_blasts[${pair.sourceIndex}]`
              });
              resolved.push({
                path: blastPath,
                name: pair.filename,
                text: blastText
              });
            }
            pendingLosatCacheCommit = { cacheInfo, cacheMap, derivedCacheMap: null };
            return resolved;
          };

          normalizeConservationState();
          if (sourceMode === 'upload') {
            const blastFiles = circularConservationSourceFiles;
            if (blastFiles.length === 0) {
              throw new Error('Please upload at least one BLAST outfmt 6/7 file for Pairwise Comparisons.');
            }
            const rawReference = String(circularConservation.reference || 'auto').trim().toLowerCase();
            circularConservation.reference = ['query', 'subject'].includes(rawReference)
              ? rawReference
              : 'auto';
            losatCacheInfo.value = [];
          } else {
            const comparisonFiles = circularConservationSourceFiles;
            if (comparisonFiles.length === 0) {
              throw new Error('Please upload at least one comparison FASTA file for Pairwise Comparisons.');
            }
            const conservationEntries = orderedConservationSources(comparisonFiles, circularConservation);
            const conservationSeries = buildConservationSeries(comparisonFiles, circularConservation);
            const conservationResults = await runCircularLosatConservation(conservationEntries);
            circularConservation.reference = 'subject';
            resolvedCircularConservation = conservationResults.map((result, index) => {
              const source = conservationEntries[index] || {};
              const style = conservationSeries[index] || {};
              return {
                name: result.name || getPayloadName(result.path),
                text: result.text,
                sourceIndex: Number(source.sourceIndex ?? index),
                label: String(style.label || `Comparison ${index + 1}`),
                color: String(style.color || '#D9EAF7')
              };
            });
          }
        } else {
          losatCacheInfo.value = [];
        }
      } else {
        requireCurrentLinearTrackLayout(form.linear_track_layout);
        const useLinearTrackSlots = adv.linear_track_slots_enabled === true;
        let linearTrackSlots = [];
        let linearTrackSlotAxisIndex = null;
        let linearSlotNeedsDepth = false;
        if (useLinearTrackSlots) {
          linearTrackSlots = normalizeLinearTrackSlots(
            adv.linear_track_slots,
            adv.nt,
            form.linear_track_layout
          );
          linearSlotNeedsDepth = linearTrackSlots.some(
            (slot) => slot.enabled !== false && slot.renderer === 'depth'
          );
          const depthTrackCount = depthTrackMatrixWidth(
            linearSeqs.map((seq) => depthFileSlotsFromValue(seq.depth))
          );
          let nextDepthIndex = 0;
          linearTrackSlots.forEach((slot) => {
            if (slot.renderer !== 'depth') return;
            slot.params = slot.params && typeof slot.params === 'object' ? { ...slot.params } : {};
            if (slot.enabled === false && slot.depth_binding_error) {
              delete slot.params.track_index;
              return;
            }
            const rawTrackIndex = Number(slot.params.track_index);
            if (!Number.isInteger(rawTrackIndex) || rawTrackIndex < 0) {
              slot.params.track_index = nextDepthIndex;
            }
            syncLinearDepthSlotHeightFromTrackConfig(slot);
            nextDepthIndex = Number(slot.params.track_index) + 1;
          });
          syncDepthSlotLegendLabelsFromTrackConfigs(
            linearTrackSlots,
            linearDepthRepresentativeFiles()
          );
          validateEnabledDepthSlots(linearTrackSlots, depthTrackCount, 'Linear');
          linearTrackSlotAxisIndex = clampLinearTrackAxisIndex(
            adv.linear_track_slots_axis_index,
            linearTrackSlots.length
          );
          linearTrackSlotAxisIndex = resolveLinearTrackAxisIndex(
            linearTrackSlots,
            linearTrackSlotAxisIndex
          );
          linearTrackSlots = applyLinearTrackOrderPlacements(
            linearTrackSlots,
            linearTrackSlotAxisIndex,
            adv.nt,
            form.linear_track_layout
          );
          adv.linear_track_slots_axis_index = linearTrackSlotAxisIndex;
          adv.linear_track_slots.splice(0, adv.linear_track_slots.length, ...linearTrackSlots);
        }
        const comparisonResolution = activeComparisonPlanSnapshot;
        const hasComparisonIntent = comparisonResolution.hasComparisonIntent === true;
        const hasLosatIntent = comparisonResolution.hasLosatIntent === true;
        const useProteinBlastp = hasLosatIntent && losatProgram.value === 'blastp';
        const blastpMode = useProteinBlastp
          ? requireCurrentProteinBlastpMode(losat.blastp?.mode)
          : String(losat.blastp?.mode ?? '');
        const usePairwiseBlastp = useProteinBlastp && blastpMode === 'pairwise';
        const useOrthogroupBlastp = useProteinBlastp && blastpMode === 'orthogroup';
        const useCollinearBlastp = useProteinBlastp && blastpMode === 'collinear';
        if (hasComparisonIntent) {
          adv.pairwise_match_style = normalizePairwiseMatchStyle(adv.pairwise_match_style);
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
        }

        const blastpMaxHits = usePairwiseBlastp
          ? requireCurrentProteinBlastpMaxHits(losat.blastp?.maxHits)
          : 5;
        const blastpCandidateLimit = useProteinBlastp
          ? resolveProteinBlastpCandidateLimit(losat.blastp?.candidateLimit)
          : null;
        const orthogroupMembershipMode = useOrthogroupBlastp || useCollinearBlastp
          ? requireCurrentOrthogroupMembershipMode(
              losat.blastp?.orthogroupMembershipMode
            )
          : 'anchor_core_v1';
        const orthogroupMemberMaxHits = useOrthogroupBlastp || useCollinearBlastp
          ? requireCurrentOrthogroupMemberMaxHits(
              losat.blastp?.orthogroupMemberMaxHits
            )
          : 5;
        const collinearMinAnchors = useCollinearBlastp
          ? requireCurrentCollinearMinAnchors(losat.blastp?.collinearMinAnchors)
          : 1;
        const collinearMaxUnitGap = useCollinearBlastp
          ? requireCurrentCollinearMaxUnitGap(losat.blastp?.collinearMaxUnitGap)
          : 0;
        const collinearMaxDiagonalDrift = useCollinearBlastp
          ? requireCurrentCollinearMaxDiagonalDrift(
              losat.blastp?.collinearMaxDiagonalDrift
            )
          : 0;
        const collinearMaxConflictsInMergeGap = useCollinearBlastp
          ? requireCurrentCollinearMaxConflicts(
              losat.blastp?.collinearMaxConflictsInMergeGap
            )
          : 1;
        const collinearMaxParalogLinksPerOrthogroup = useCollinearBlastp
          ? requireCurrentCollinearMaxParalogLinks(
              losat.blastp?.collinearMaxParalogLinksPerOrthogroup
            )
          : 2;
        const collinearColorMode = useCollinearBlastp
          ? requireCurrentCollinearColorMode(losat.blastp?.collinearColorMode)
          : 'orientation';
        const collinearUnitMode = useCollinearBlastp
          ? requireCurrentCollinearUnitMode(losat.blastp?.collinearUnitMode)
          : 'auto';
        const collinearAnchorMode = useCollinearBlastp
          ? requireCurrentCollinearAnchorMode(losat.blastp?.collinearAnchorMode)
          : 'rbh';
        const collinearMergeOrientation = useCollinearBlastp
          ? requireCurrentCollinearMergeOrientation(
              losat.blastp?.collinearMergeOrientation
            )
          : 'either';
        const collinearSearchScope = useCollinearBlastp
          ? requireCurrentCollinearSearchScope(losat.blastp?.collinearSearchScope)
          : 'adjacent';

        const reuseResolvedProteinArtifacts = useProteinBlastp
          && !legacyProteinRawCandidates.value?.entries?.length
          && canReuseResolvedProteinArtifacts({
            canonicalComparisons: files.linearCanonicalComparisons,
            committedRequest: typeof getCommittedCanonicalRenderRequest === 'function'
              ? getCommittedCanonicalRenderRequest()
              : null,
            active: {
              mode: blastpMode,
              candidateLimit: blastpCandidateLimit,
              bitscore: adv.min_bitscore,
              evalue: adv.evalue,
              identity: adv.identity,
              alignmentLength: adv.alignment_length,
              maxHits: blastpMaxHits,
              orthogroupMembershipMode,
              memberMaxHits: orthogroupMemberMaxHits,
              minAnchors: collinearMinAnchors,
              maxUnitGap: collinearMaxUnitGap,
              maxDiagonalDrift: collinearMaxDiagonalDrift,
              maxConflicts: collinearMaxConflictsInMergeGap,
              maxParalogLinks: collinearMaxParalogLinksPerOrthogroup,
              colorMode: collinearColorMode,
              unitMode: collinearUnitMode,
              anchorMode: collinearAnchorMode,
              mergeOrientation: collinearMergeOrientation,
              searchScope: collinearSearchScope
            }
          });
        const useLosat = hasLosatIntent && !reuseResolvedProteinArtifacts;
        if (useLosat) recordSessionLifecycleEvent('losat-cache-preparation-start');

        if (useProteinBlastp) {
          losat.blastp.mode = blastpMode;
          losat.blastp.candidateLimit = blastpCandidateLimit;
        }
        if (usePairwiseBlastp) {
          losat.blastp.maxHits = blastpMaxHits;
        } else if (useOrthogroupBlastp) {
          losat.blastp.orthogroupMembershipMode = orthogroupMembershipMode;
          losat.blastp.orthogroupMemberMaxHits = orthogroupMemberMaxHits;
        } else if (useCollinearBlastp) {
          losat.blastp.orthogroupMembershipMode = orthogroupMembershipMode;
          losat.blastp.orthogroupMemberMaxHits = orthogroupMemberMaxHits;
          losat.blastp.collinearMinAnchors = collinearMinAnchors;
          losat.blastp.collinearMaxUnitGap = collinearMaxUnitGap;
          losat.blastp.collinearMaxDiagonalDrift = collinearMaxDiagonalDrift;
          losat.blastp.collinearMaxConflictsInMergeGap = collinearMaxConflictsInMergeGap;
          losat.blastp.collinearMaxParalogLinksPerOrthogroup =
            collinearMaxParalogLinksPerOrthogroup;
          losat.blastp.collinearColorMode = collinearColorMode;
          losat.blastp.collinearUnitMode = collinearUnitMode;
          losat.blastp.collinearAnchorMode = collinearAnchorMode;
          losat.blastp.collinearMergeOrientation = collinearMergeOrientation;
          losat.blastp.collinearSearchScope = collinearSearchScope;
        }

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
        adv.linear_definition_line_styles = normalizeDefinitionLineStyleState(adv.linear_definition_line_styles);
        form.plot_title = normalizedPlotTitle;
        adv.plot_title_position = normalizedPlotTitlePosition;
        adv.plot_title_font_size = normalizedPlotTitleFontSize;

        const normalizedLabelPlacement = requireCurrentLinearLabelPlacement(
          adv.label_placement
        );
        let normalizedLabelRendering = form.show_labels_linear === 'none' ? 'auto' : normalizeLabelRendering(adv.label_rendering);
        if (normalizedLabelPlacement === 'above_feature') {
          normalizedLabelRendering = 'auto';
        }
        adv.label_rendering = normalizedLabelRendering;
        const normalizedLinearLabelSpacing = normalizePositiveNumberOrNull(adv.linear_label_spacing);
        adv.linear_label_spacing = normalizedLinearLabelSpacing;
        if (hasComparisonIntent) {
          const comparisonHeight = classifyOptionalPositiveNumber(adv.comparison_height);
          if (comparisonHeight.status === 'invalid') {
            throw new Error('Pairwise Match Height must be Auto or a positive finite number.');
          }
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
            reverseFlags.push(false);
            viewTransformSpecs.push({ reverse: displayReverse });
            return { file: canonicalSpecBody, displayFile: specBody };
          }

          reverseFlags.push(wantsReverse);
          viewTransformSpecs.push({ reverse: wantsReverse });
          return null;
        };

        const fastaCache = new Map();
        const fastaHashCache = new Map();
        const sequenceEntriesByKey = new Map();
        const linearFileTextCache = new WeakMap();
        if (!reuseResolvedProteinArtifacts) {
          collinearGroups.value = [];
          if (useOrthogroupBlastp) {
            clearOrthogroupMetadata();
          } else {
            clearOrthogroupMetadata({ clearSelection: true });
          }
        }
        let cacheInfo = [];
        const cacheMap = new Map(losatCache.value || []);
        const derivedCacheMap = useProteinBlastp
          ? new Map(losatDerivedCache.value || [])
          : null;
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
              fastaWorkerFallbacks: 0,
              proteinDerivedPayloadCacheHits: 0,
              proteinDerivedPayloadCacheMisses: 0,
              proteinConversionCacheHits: 0,
              proteinFilteredHitCacheHits: 0,
              proteinFilteredHitCacheMisses: 0,
              rawTsvEntryCount: 0,
              rawTsvBytes: 0,
              rawTsvLargestEntryBytes: 0,
              simultaneousParsedTables: 0,
              helperRequestMetadataBytes: 0,
              helperRequestRawTransferBytes: 0,
              helperRequestFileCount: 0,
              rawJobs: [],
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

        if (!hasLosatIntent) {
          losatCacheInfo.value = [];
        }

        let proteinRecordInstanceKeys = [];

        const buildProteinRecordInstanceKeys = async () => {
          const used = new Set();
          return linearSeqs.map((sequence, index) => {
            const base = String(sequence?.uid || `record-${index + 1}`).trim() || `record-${index + 1}`;
            let key = base;
            let suffix = 2;
            while (used.has(key)) {
              key = `${base}-${suffix}`;
              suffix += 1;
            }
            used.add(key);
            return key;
          });
        };

        const getSeqEntry = async (idx) => {
          if (fastaCache.has(idx)) return fastaCache.get(idx);
          const startedAt = getNow();
          const fmt = lInputType.value === 'gb'
            ? 'genbank'
            : (useProteinBlastp ? 'gff' : 'fasta');
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
                pairedFastaFile: getFileFingerprint(pairedFastaFile),
                regionSpec,
                recordSelector,
                recordInstanceKey,
                recordIndex: idx,
                featureVisibility: featureVisibilityCacheKey,
                proteinMapSchema: 4
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
              const helperFiles = [{
                role: 'source',
                bytes: await cloneFileBytesForTransfer(sourceFile)
              }];
              if (pairedFastaFile) {
                helperFiles.push({
                  role: 'fasta',
                  bytes: await cloneFileBytesForTransfer(pairedFastaFile)
                });
              }
              if (featureVisibilityTablePath && featureVisibilityCacheKey) {
                helperFiles.push({
                  role: 'visibility',
                  bytes: textEncoder.encode(featureVisibilityCacheKey).buffer
                });
              }
              const response = await runDiagramHelperOperation(
                DIAGRAM_HELPER_OPERATIONS.EXTRACT_CDS_PROTEIN_FASTA,
                {
                  files: helperFiles,
                  format: fmt,
                  regionSpec,
                  recordSelector,
                  reverseFlag: reverseFlag === '1',
                  recordIndex: idx,
                  recordInstanceKey
                }
              );
              const res = response.result;
              if (res.error) throw new Error(res.error);
              const fastaHash = await hashText(res.fasta || '');
              const proteinCacheKey = String(res.display_binding_hash || '');
              entry = {
                fasta: res.fasta,
                recordId: res.record_id || `seq_${idx + 1}`,
                canonicalLength: Number(res.record_length || 0) || getFastaSequenceLength(res.fasta),
                proteinMap: res.protein_map || {},
                proteinCount: res.protein_count || 0,
                proteinCacheKey,
                proteinSetHash: String(res.protein_set_hash || ''),
                recordAnalysisId: String(res.record_analysis_id || ''),
                recordInstanceKey: String(res.record_instance_key || recordInstanceKey),
                runtimeBindingHash: String(res.runtime_binding_hash || ''),
                displayBindingHash: String(res.display_binding_hash || ''),
                identityManifest: res.identity_manifest || null,
                sequenceKey: `protein:${fastaHash}`,
                hash: fastaHash
              };
              if (losatTiming) losatTiming.fastaWorkerFallbacks += 1;
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
                const response = await runDiagramHelperOperation(
                  DIAGRAM_HELPER_OPERATIONS.EXTRACT_FIRST_FASTA,
                  {
                    files: [{
                      role: 'source',
                      bytes: await cloneFileBytesForTransfer(sourceFile)
                    }],
                    format: fmt,
                    regionSpec,
                    recordSelector,
                    reverseFlag: reverseFlag === '1'
                  }
                );
                const res = response.result;
                if (res.error) throw new Error(res.error);
                entry = {
                  fasta: res.fasta,
                  recordId: res.record_id || `seq_${idx + 1}`,
                  canonicalLength: Number(res.record_length || 0) || getFastaSequenceLength(res.fasta)
                };
                if (losatTiming) {
                  losatTiming.fastaWorkerFallbacks += 1;
                  console.warn('LOSAT browser FASTA extraction fell back to the diagram Worker:', fastError);
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

        const prepareLinearFile = async (fileObj, path, { cacheText = false, slot = '' } = {}) => {
          return stageUploadedFile(fileObj, path, {
            cacheText,
            textCache: linearFileTextCache,
            slot
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

        const buildCacheMetadata = async (argsKey, queryIdx, subjectIdx) => {
          if (useProteinBlastp) {
            const queryEntry = await getSeqEntry(queryIdx);
            const subjectEntry = await getSeqEntry(subjectIdx);
            return {
              identityKind: 'protein',
              program: 'blastp',
              outfmt: String(losat.outfmt || '6'),
              args: argsKey,
              queryProteinSetHash: queryEntry.proteinSetHash,
              subjectProteinSetHash: subjectEntry.proteinSetHash,
              queryRuntimeBindingHash: queryEntry.runtimeBindingHash,
              subjectRuntimeBindingHash: subjectEntry.runtimeBindingHash,
              queryRecordInstanceKey: queryEntry.recordInstanceKey,
              subjectRecordInstanceKey: subjectEntry.recordInstanceKey
            };
          }
          const queryHash = await getSeqHash(queryIdx);
          const subjectHash = await getSeqHash(subjectIdx);
          return {
            identityKind: 'nucleotide',
            program: losatProgram.value,
            outfmt: String(losat.outfmt || '6'),
            args: argsKey,
            queryCanonicalHash: queryHash,
            subjectCanonicalHash: subjectHash
          };
        };

        const buildCacheKey = async (metadata) => {
          if (metadata?.identityKind === 'protein') {
            const response = await runDiagramHelperOperation(
              DIAGRAM_HELPER_OPERATIONS.BUILD_PROTEIN_LOSAT_CACHE_KEY,
              {
                identityManifest: cloneJsonData(proteinIdentityManifest.value),
                queryRecordInstanceKey: metadata.queryRecordInstanceKey,
                subjectRecordInstanceKey: metadata.subjectRecordInstanceKey,
                expectedOptions: {
                  program: metadata.program,
                  outfmt: metadata.outfmt,
                  args: normalizeLosatArgs(metadata.args)
                }
              }
            );
            const result = response.result;
            if (result.error || !result.key) {
              throw new Error(result.error || 'Protein cache key generation failed.');
            }
            return String(result.key);
          }
          return hashText(JSON.stringify(buildLosatCachePayload(metadata)));
        };

        const tryPromoteLegacyProteinEntry = async ({
          cacheKey,
          metadata,
          queryEntry,
          subjectEntry
        }) => {
          if (
            metadata?.identityKind !== 'protein' ||
            !legacyProteinRawCandidates?.value
          ) return null;
          const envelope = legacyProteinRawCandidates.value;
          const pending = (Array.isArray(envelope?.entries) ? envelope.entries : [])
            .map((candidate, index) => ({ candidate, index }))
            .filter(({ candidate }) => candidate?.state === 'pending' && candidate?.originalEntry)
            .map(({ candidate, index }) => ({
              candidateIndex: index,
              entry: candidate.originalEntry
            }));
          if (pending.length === 0) return null;

          const response = await runDiagramHelperOperation(
            DIAGRAM_HELPER_OPERATIONS.PROMOTE_LEGACY_LOSATP_CACHE,
            {
              candidates: cloneJsonData(pending),
              queryFasta: String(queryEntry?.fasta || ''),
              subjectFasta: String(subjectEntry?.fasta || ''),
              queryProteinMap: cloneJsonData(queryEntry?.proteinMap || {}),
              subjectProteinMap: cloneJsonData(subjectEntry?.proteinMap || {}),
              identityManifest: cloneJsonData(proteinIdentityManifest.value),
              expectedOptions: {
                program: metadata.program,
                outfmt: metadata.outfmt,
                args: normalizeLosatArgs(metadata.args)
              }
            }
          );
          const result = response.result;
          if (result.status === 'error') {
            throw new Error(result.error || 'Legacy protein cache migration failed.');
          }
          const candidateIndex = Number(result.candidateIndex);
          if (
            result.status !== 'promoted' ||
            !Number.isInteger(candidateIndex) ||
            candidateIndex < 0 ||
            typeof result.text !== 'string' ||
            !result.entry ||
            typeof result.entry !== 'object'
          ) return null;
          if (String(result.entry.key || '') !== cacheKey) {
            legacyProteinRawCandidates.value = transitionLegacyProteinCandidate(
              envelope,
              candidateIndex,
              'rejected',
              'Promoted legacy key does not match the current directional cache key.'
            );
            return null;
          }

          const promoted = {
            ...result.entry,
            text: result.text,
            migratedFromSchema: 2
          };
          delete promoted.key;
          delete promoted.filename;
          delete promoted.display;
          cacheMap.set(cacheKey, promoted);
          const verified = getCurrentRawLosatCacheEntry(
            cacheMap,
            cacheKey,
            metadata,
            proteinIdentityManifest.value
          );
          if (!verified) {
            cacheMap.delete(cacheKey);
            legacyProteinRawCandidates.value = transitionLegacyProteinCandidate(
              envelope,
              candidateIndex,
              'rejected',
              'Rewritten legacy TSV does not resolve through the current protein manifest.'
            );
            return null;
          }
          legacyPromotionTransaction.push({
            cacheMap,
            cacheKey,
            candidateIndex,
            proteinIdMap: result.proteinIdMap || {}
          });
          return verified;
        };

        const buildCacheFilename = (spec, queryEntry, subjectEntry) => {
          const edge = comparisonResolution.edges.find(
            (candidate) => candidate.edgeKey === spec.edgeKey
          );
          const fallback = getLosatPairDefaultName(
            spec.edgeKey,
            queryEntry,
            subjectEntry
          );
          return normalizeLosatFilename(
            edge?.losatFilenameActive ? edge.losatFilename : '',
            fallback
          );
        };

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
          return blastpCandidateLimit;
        };

        const buildLosatArgs = (queryIdx, subjectIdx) => {
          const args = [];
          if (losatProgram.value === 'blastn') {
            pushArg(args, '--task', losat.blastn.task);
          } else if (losatProgram.value === 'tblastx') {
            pushArg(args, '--query-gencode', getGencode(queryIdx));
            pushArg(args, '--db-gencode', getGencode(subjectIdx));
          } else {
            if (!useOrthogroupBlastp && !useCollinearBlastp) {
              pushArg(args, '--max-hsps-per-subject', 1);
            }
            pushArg(args, '--max-target-seqs', getBlastpCandidateLimit());
          }
          return args;
        };

        {
          const inputWriteStartedAt = getNow();
          const losatRecordIndexes = (useOrthogroupBlastp || useCollinearBlastp)
            ? new Set(linearSeqs.map((_, index) => index))
            : new Set(comparisonResolution.edges
                .filter((edge) => edge.source === 'losat')
                .flatMap((edge) => [edge.queryIndex, edge.subjectIndex]));
          for (let i = 0; i < linearSeqs.length; i++) {
            const seq = linearSeqs[i];
            if (lInputType.value === 'gb') {
              if (!seq.gb) throw new Error(`Sequence #${i + 1}: Missing GenBank file.`);
              if (useLosat && losatRecordIndexes.has(i)) {
                await prepareLinearFile(seq.gb, `/seq_${i}.gb`, {
                  cacheText: true,
                  slot: `files.linearSeqs[${i}].gb`
                });
              }
            } else {
              if (!seq.gff || !seq.fasta) throw new Error(`Sequence #${i + 1}: GFF3 and FASTA are required.`);
              if (useLosat && losatRecordIndexes.has(i)) {
                await prepareLinearFile(seq.gff, `/seq_${i}.gff`, {
                  slot: `files.linearSeqs[${i}].gff`
                });
                await prepareLinearFile(seq.fasta, `/seq_${i}.fasta`, {
                  cacheText: true,
                  slot: `files.linearSeqs[${i}].fasta`
                });
              }
            }
          }
          if (losatTiming) losatTiming.inputWriteMs += getNow() - inputWriteStartedAt;
        }

        regionSpecs = linearSeqs.map((seq, idx) => buildRegionSpec(seq, idx));
        if (useLosat && useProteinBlastp) {
          proteinRecordInstanceKeys = await buildProteinRecordInstanceKeys();
          const proteinRecordIndexes = (useOrthogroupBlastp || useCollinearBlastp)
            ? linearSeqs.map((_, index) => index)
            : Array.from(new Set(comparisonResolution.edges
                .filter((edge) => edge.source === 'losat')
                .flatMap((edge) => [edge.queryIndex, edge.subjectIndex])))
                .sort((left, right) => left - right);
          const proteinEntries = [];
          for (const index of proteinRecordIndexes) {
            proteinEntries.push(await getSeqEntry(index));
          }
          const manifests = proteinEntries.map((entry) => entry.identityManifest);
          if (manifests.some((manifest) => !validateProteinIdentityManifest(manifest))) {
            throw new Error(
              'Protein comparison metadata could not be validated. Reload the page and try again.'
            );
          }
          proteinIdentityManifest.value = mergeProteinIdentityManifests(manifests);
          const legacyReferenceIds = collectLegacyProteinReferences(
            orthogroups.value,
            selectedOrthogroupAlignmentFeature.value,
            extractedFeatures.value,
            biologicalFeatures?.value
          );
          if (legacyReferenceIds.length > 0) {
            const response = await runDiagramHelperOperation(
              DIAGRAM_HELPER_OPERATIONS.RESOLVE_LEGACY_PROTEIN_REFERENCES,
              {
                proteinRecords: cloneJsonData(proteinEntries.map((entry) => ({
                  proteinMap: entry.proteinMap || {},
                  fasta: entry.fasta || ''
                }))),
                identityManifest: cloneJsonData(proteinIdentityManifest.value),
                referenceIds: cloneJsonData(legacyReferenceIds)
              }
            );
            const result = response.result;
            if (result.status !== 'resolved' || !result.proteinIdMap) {
              throw new Error(
                result.error || 'Legacy protein UI reference migration failed.'
              );
            }
            const proteinIdMap = result.proteinIdMap;
            orthogroups.value = rewriteMappedProteinReferences(
              orthogroups.value,
              proteinIdMap
            );
            selectedOrthogroupAlignmentFeature.value = rewriteMappedProteinReferences(
              selectedOrthogroupAlignmentFeature.value,
              proteinIdMap
            );
            extractedFeatures.value = rewriteMappedProteinReferences(
              extractedFeatures.value,
              proteinIdMap
            );
            if (biologicalFeatures) {
              biologicalFeatures.value = rewriteMappedProteinReferences(
                biologicalFeatures.value,
                proteinIdMap
              );
            }
            if (collectLegacyProteinReferences(
              orthogroups.value,
              selectedOrthogroupAlignmentFeature.value,
              extractedFeatures.value,
              biologicalFeatures?.value
            ).length > 0) {
              throw new Error('Legacy protein UI references remain after migration.');
            }
          }
        }

        if (useLosat) {
          setProcessingStatus('Preparing LOSAT jobs...');
          const losatPairs = [];
          const losatJobs = [];
          const pendingJobKeys = new Set();
          const jobBuildStartedAt = getNow();
          const fastaExtractionBeforeJobBuild = losatTiming.fastaExtractionMs;
          const cacheHashBeforeJobBuild = losatTiming.cacheHashMs;

          const jobSpecs = [];
          const resolvedLosatEdges = comparisonResolution.edges.filter(
            (edge) => edge.source === 'losat'
          );
          const edgeForOrdinal = (ordinal) => (
            resolvedLosatEdges.find((edge) => edge.ordinal === ordinal) ||
            resolvedLosatEdges[Math.max(0, Math.min(ordinal, resolvedLosatEdges.length - 1))]
          );
          const pushExpandedJobSpec = (queryIndex, subjectIndex, ordinal) => {
            const edge = edgeForOrdinal(ordinal);
            if (!edge) return;
            jobSpecs.push({
              edgeKey: edge.edgeKey,
              ordinal: edge.ordinal,
              queryUid: edge.queryUid,
              subjectUid: edge.subjectUid,
              queryIndex,
              subjectIndex,
              program: losatProgram.value
            });
          };
          if (useOrthogroupBlastp) {
            for (let i = 0; i < linearSeqs.length; i++) {
              pushExpandedJobSpec(i, i, Math.min(i, Math.max(0, resolvedLosatEdges.length - 1)));
              for (let j = i + 1; j < linearSeqs.length; j++) {
                pushExpandedJobSpec(i, j, Math.min(i, Math.max(0, resolvedLosatEdges.length - 1)));
                pushExpandedJobSpec(j, i, Math.min(i, Math.max(0, resolvedLosatEdges.length - 1)));
              }
            }
          } else if (useCollinearBlastp) {
            for (let i = 0; i < linearSeqs.length; i++) {
              pushExpandedJobSpec(i, i, Math.min(i, Math.max(0, resolvedLosatEdges.length - 1)));
            }
            if (collinearSearchScope === 'all') {
              for (let i = 0; i < linearSeqs.length - 1; i++) {
                for (let j = i + 1; j < linearSeqs.length; j++) {
                  pushExpandedJobSpec(i, j, Math.min(i, Math.max(0, resolvedLosatEdges.length - 1)));
                  pushExpandedJobSpec(j, i, Math.min(i, Math.max(0, resolvedLosatEdges.length - 1)));
                }
              }
            } else {
              resolvedLosatEdges.forEach((edge) => {
                pushExpandedJobSpec(edge.queryIndex, edge.subjectIndex, edge.ordinal);
                pushExpandedJobSpec(edge.subjectIndex, edge.queryIndex, edge.ordinal);
              });
            }
          } else {
            jobSpecs.push(...buildPairwiseLosatJobSpecs({
              resolution: comparisonResolution,
              program: losatProgram.value,
              blastpMode
            }));
          }

          for (const spec of jobSpecs) {
            throwIfGenerationCanceled();
            const queryEntry = await getSeqEntry(spec.queryIndex);
            throwIfGenerationCanceled();
            const subjectEntry = await getSeqEntry(spec.subjectIndex);
            throwIfGenerationCanceled();
            const losatArgs = buildLosatArgs(spec.queryIndex, spec.subjectIndex);
            const cacheMetadata = await buildCacheMetadata(losatArgs, spec.queryIndex, spec.subjectIndex);
            const cacheKey = await buildCacheKey(cacheMetadata);
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
            let cached = getRawLosatCacheEntry(
              cacheMap,
              cacheKey,
              cacheMetadata,
              proteinIdentityManifest.value
            );
            if (!cached && useProteinBlastp) {
              cached = await tryPromoteLegacyProteinEntry({
                cacheKey,
                metadata: cacheMetadata,
                queryEntry,
                subjectEntry
              });
              throwIfGenerationCanceled();
            }
            const hasCachedText = Boolean(cached);
            if (cached) promoteRawLosatCacheEntry(cacheMap, cacheKey, cached, cacheMetadata);
            losatTiming.totalPairs += 1;
            if (hasCachedText) losatTiming.cacheHits += 1;
            else losatTiming.cacheMisses += 1;
            const resolvedEdge = comparisonResolution.edges.find(
              (edge) => edge.edgeKey === spec.edgeKey
            );
            const isResolvedDisplayPair = Boolean(
              resolvedEdge &&
              spec.queryIndex === resolvedEdge.queryIndex &&
              spec.subjectIndex === resolvedEdge.subjectIndex
            );
            const pair = {
              pairIndex: spec.ordinal,
              ordinal: spec.ordinal,
              edgeKey: spec.edgeKey,
              queryUid: spec.queryUid,
              subjectUid: spec.subjectUid,
              queryIndex: spec.queryIndex,
              subjectIndex: spec.subjectIndex,
              cacheKey,
              filename: buildCacheFilename(spec, queryEntry, subjectEntry),
              displayPair: isResolvedDisplayPair
            };
            losatPairs.push(pair);
            losatTiming.rawJobs.push({
              queryRecordIndex: spec.queryIndex,
              subjectRecordIndex: spec.subjectIndex,
              cacheKey,
              args: [...losatArgs]
            });
            if (
              isResolvedDisplayPair &&
              !cacheInfo.some((entry) => entry.edgeKey === spec.edgeKey)
            ) {
              cacheInfo.push({
                key: cacheKey,
                filename: pair.filename,
                display: true,
                edgeKey: spec.edgeKey,
                ordinal: spec.ordinal,
                queryUid: resolvedEdge.queryUid,
                subjectUid: resolvedEdge.subjectUid,
                queryIndex: resolvedEdge.queryIndex,
                subjectIndex: resolvedEdge.subjectIndex
              });
            }

            if (!hasCachedText && !pendingJobKeys.has(cacheKey)) {
              pendingJobKeys.add(cacheKey);
              losatJobs.push({
                pairIndex: spec.ordinal,
                ordinal: spec.ordinal,
                edgeKey: spec.edgeKey,
                queryUid: spec.queryUid,
                subjectUid: spec.subjectUid,
                queryIndex: spec.queryIndex,
                subjectIndex: spec.subjectIndex,
                cacheKey,
                program: losatProgram.value,
                querySequenceKey: queryEntry.sequenceKey,
                subjectSequenceKey: subjectEntry.sequenceKey,
                queryCanonicalHash,
                subjectCanonicalHash,
                outfmt: losat.outfmt || '6',
                extraArgs: losatArgs,
                cacheMetadata
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
            const losatResults = await executeLosatJobs(losatJobs, {
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
              const cacheMetadata = job?.cacheMetadata || {};
              const isProteinEntry = cacheMetadata.identityKind === 'protein';
              const rawEntry = {
                schema: isProteinEntry
                  ? PROTEIN_LOSAT_CACHE_SCHEMA
                  : NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
                kind: 'raw-losat',
                identityKind: isProteinEntry ? 'protein' : 'nucleotide',
                ...(isProteinEntry ? { idEncoding: 'runtime-handle-v1' } : {}),
                text: result.text,
                program: losatProgram.value,
                outfmt: String(losat.outfmt || '6'),
                args: job?.extraArgs || [],
                ...(isProteinEntry
                  ? {
                      queryProteinSetHash: cacheMetadata.queryProteinSetHash,
                      subjectProteinSetHash: cacheMetadata.subjectProteinSetHash,
                      queryRuntimeBindingHash: cacheMetadata.queryRuntimeBindingHash,
                      subjectRuntimeBindingHash: cacheMetadata.subjectRuntimeBindingHash,
                      queryRecordInstanceKey: cacheMetadata.queryRecordInstanceKey,
                      subjectRecordInstanceKey: cacheMetadata.subjectRecordInstanceKey,
                      queryCanonicalHash: job?.queryCanonicalHash || '',
                      subjectCanonicalHash: job?.subjectCanonicalHash || ''
                    }
                  : {
                      queryCanonicalHash: job?.queryCanonicalHash || '',
                      subjectCanonicalHash: job?.subjectCanonicalHash || ''
                    })
              };
              cacheMap.set(result.cacheKey, rawEntry);
            });
          } else {
            setProcessingStatus('Using cached LOSAT results...');
          }

          const blastWriteStartedAt = getNow();
          throwIfGenerationCanceled();
          if (useProteinBlastp) {
            const recordPayloads = [];
            const recordIndexes = (useOrthogroupBlastp || useCollinearBlastp)
              ? linearSeqs.map((_, index) => index)
              : Array.from(new Set(losatPairs.flatMap((pair) => [
                  pair.queryIndex,
                  pair.subjectIndex
                ]))).sort((left, right) => left - right);
            for (const i of recordIndexes) {
              const entry = await getSeqEntry(i);
              recordPayloads.push({
                recordIndex: i,
                recordId: entry.recordId || `seq_${i + 1}`,
                proteinMap: entry.proteinMap || {},
                proteinCacheKey: entry.proteinCacheKey || entry.sequenceKey || `record:${i}`,
                runtimeBindingHash: entry.runtimeBindingHash || '',
                displayBindingHash: entry.displayBindingHash || '',
                viewTransform: await getViewTransform(i)
              });
            }
            const pairPayloads = [];
            const rawTsvParts = [];
            let rawTsvOffset = 0;
            let rawTsvLargestEntryBytes = 0;
            for (const pair of losatPairs) {
              throwIfGenerationCanceled();
              const cached = cacheMap.get(pair.cacheKey);
              const losatText = isCurrentRawLosatCacheEntry(cached) ? cached.text : '';
              const rawTsvBytes = new Blob([losatText]).size;
              pairPayloads.push({
                pairIndex: pair.pairIndex,
                ordinal: pair.ordinal,
                edgeKey: pair.edgeKey,
                queryIndex: pair.queryIndex,
                subjectIndex: pair.subjectIndex,
                displayPair: pair.displayPair,
                cacheKey: pair.cacheKey,
                rawTsvOffset,
                rawTsvBytes
              });
              rawTsvParts.push(losatText);
              rawTsvOffset += rawTsvBytes;
              rawTsvLargestEntryBytes = Math.max(rawTsvLargestEntryBytes, rawTsvBytes);
            }
            losatTiming.rawTsvEntryCount = pairPayloads.length;
            losatTiming.rawTsvBytes = rawTsvOffset;
            losatTiming.rawTsvLargestEntryBytes = rawTsvLargestEntryBytes;
            setProcessingStatus(
              useCollinearBlastp
                ? 'Converting LOSAT protein links to collinear ribbons...'
                : 'Converting LOSAT protein hits...'
            );
            const useDerivedProteinPayloadCache = useOrthogroupBlastp || useCollinearBlastp;
            let derivedCacheKey = '';
            let convertedPayload = null;
            if (useDerivedProteinPayloadCache) {
              const derivedCachePayload = buildLosatDerivedPayloadCachePayload({
                mode: blastpMode,
                maxHits: blastpMaxHits,
                bitscore: adv.min_bitscore,
                evalue: adv.evalue,
                identity: adv.identity,
                alignmentLength: adv.alignment_length,
                collinearMinAnchors,
                collinearMaxUnitGap,
                collinearUnitMode,
                collinearColorMode,
                collinearAnchorMode,
                collinearMergeOrientation,
                collinearMaxDiagonalDrift,
                collinearMaxConflictsInMergeGap,
                collinearMaxParalogLinksPerOrthogroup,
                collinearSearchScope,
                orthogroupMembershipMode,
                orthogroupMemberMaxHits,
                recordPayloads,
                pairPayloads
              });
              derivedCacheKey = await hashText(JSON.stringify(derivedCachePayload));
              convertedPayload = getLosatDerivedCacheEntry(
                derivedCacheMap,
                derivedCacheKey,
                proteinIdentityManifest.value
              );
              if (convertedPayload) {
                convertedPayload = {
                  ...convertedPayload,
                  cache: {
                    ...(convertedPayload.cache || {}),
                    derivedPayloadHit: true
                  }
                };
                losatTiming.proteinDerivedPayloadCacheHits += 1;
              } else {
                losatTiming.proteinDerivedPayloadCacheMisses += 1;
              }
            }
            if (!convertedPayload) {
              const pairManifestBytes = textEncoder.encode(JSON.stringify({
                records: recordPayloads,
                pairs: pairPayloads
              }));
              const rawTsvBuffer = await new Blob(rawTsvParts, {
                type: 'text/tab-separated-values'
              }).arrayBuffer();
              losatTiming.helperRequestMetadataBytes = pairManifestBytes.byteLength;
              losatTiming.helperRequestRawTransferBytes = rawTsvBuffer.byteLength;
              losatTiming.helperRequestFileCount = 2;
              const response = await runDiagramHelperOperation(
                DIAGRAM_HELPER_OPERATIONS.CONVERT_LOSATP_PAIRS_TO_GENOMIC_PAYLOAD,
                {
                  files: [
                    { role: 'pairs', bytes: pairManifestBytes.buffer },
                    { role: 'rawTsv', bytes: rawTsvBuffer }
                  ],
                  mode: blastpMode,
                  maxHits: blastpMaxHits,
                  bitscore: adv.min_bitscore,
                  evalue: adv.evalue,
                  identity: adv.identity,
                  alignmentLength: adv.alignment_length,
                  collinearMinAnchors,
                  collinearMaxUnitGap,
                  collinearUnitMode,
                  collinearColorMode,
                  collinearAnchorMode,
                  collinearMergeOrientation,
                  collinearMaxDiagonalDrift,
                  collinearMaxConflictsInMergeGap,
                  collinearMaxParalogLinksPerOrthogroup,
                  collinearSearchScope,
                  orthogroupMembershipMode,
                  orthogroupMemberMaxHits
                }
              );
              convertedPayload = response.result;
              if (useDerivedProteinPayloadCache && !convertedPayload?.error) {
                setLosatDerivedCacheEntry(derivedCacheMap, derivedCacheKey, {
                  mode: blastpMode,
                  payload: convertedPayload,
                  manifest: proteinIdentityManifest.value
                });
              }
            }
            if (convertedPayload.error) throw new Error(convertedPayload.error);
            if (!hasRequiredCanonicalAnalysisResource(blastpMode, convertedPayload)) {
              throw new Error(
                'Protein comparison analysis did not return its canonical typed result.'
              );
            }
            const conversionCache = convertedPayload.cache || {};
            if (conversionCache.convertedPayloadHit) losatTiming.proteinConversionCacheHits += 1;
            losatTiming.proteinFilteredHitCacheHits += Number(conversionCache.filteredHitCacheHits || 0);
            losatTiming.proteinFilteredHitCacheMisses += Number(conversionCache.filteredHitCacheMisses || 0);
            losatTiming.simultaneousParsedTables = Number(
              conversionCache.simultaneousParsedTables || 0
            );
            if (useCollinearBlastp) {
              resolvedComparisons.push({
                kind: 'collinearityResult',
                typedResource: convertedPayload.collinearityResult
              });
            } else if (useOrthogroupBlastp) {
              resolvedComparisons.push({
                kind: 'orthogroupResult',
                typedResource: convertedPayload.orthogroupResult
              });
            }
            const convertedPairs = Array.isArray(convertedPayload.pairs) ? convertedPayload.pairs : [];
            for (const converted of convertedPairs) {
              const pairIndex = Number(converted?.pair_index);
              if (!Number.isInteger(pairIndex)) continue;
              const blastPath = `/blast_${pairIndex}.txt`;
              const sourcePair =
                losatPairs.find((pair) => pair.pairIndex === pairIndex && pair.displayPair) ||
                losatPairs.find((pair) => pair.pairIndex === pairIndex);
              const blastName = sourcePair?.filename || getPayloadName(blastPath);
              const blastSlot = `generatedFiles.losat_blasts[${pairIndex}]`;
              registerRunInfoFile(blastPath, {
                name: blastName,
                slot: blastSlot,
                kind: 'generated'
              });
              recordGeneratedCliFile(blastPath, converted.tsv || '', {
                name: blastName,
                slot: blastSlot
              });
              if (!useCollinearBlastp) {
                resolvedComparisons.push({
                  kind: 'precomputedProteinComparison',
                  edgeKey: sourcePair?.edgeKey || '',
                  ordinal: Number(sourcePair?.ordinal),
                  queryRecordIndex: Number(sourcePair?.queryIndex),
                  subjectRecordIndex: Number(sourcePair?.subjectIndex),
                  rows: Array.isArray(converted.rows) ? converted.rows : []
                });
              }
            }
          } else {
            for (const pair of losatPairs) {
              throwIfGenerationCanceled();
              const cached = cacheMap.get(pair.cacheKey);
              const blastText = isCurrentRawLosatCacheEntry(cached) ? cached.text : '';
              const response = await runDiagramHelperOperation(
                DIAGRAM_HELPER_OPERATIONS.CONVERT_LOSAT_NUCLEOTIDE_TO_DISPLAY_TSV,
                {
                  blastText,
                  queryViewTransform: await getViewTransform(pair.queryIndex),
                  subjectViewTransform: await getViewTransform(pair.subjectIndex)
                }
              );
              const converted = response.result;
              if (converted.error) throw new Error(converted.error);
              const blastPath = `/blast_${pair.pairIndex}.txt`;
              const blastName = pair.filename || getPayloadName(blastPath);
              const blastSlot = `generatedFiles.losat_blasts[${pair.pairIndex}]`;
              registerRunInfoFile(blastPath, {
                name: blastName,
                slot: blastSlot,
                kind: 'generated'
              });
              recordGeneratedCliFile(blastPath, converted.tsv || '', {
                name: blastName,
                slot: blastSlot
              });
              resolvedComparisons.push({
                kind: 'nucleotideBlast',
                edgeKey: pair.edgeKey,
                ordinal: pair.ordinal,
                queryRecordIndex: pair.queryIndex,
                subjectRecordIndex: pair.subjectIndex,
                text: converted.tsv || '',
              });
            }
          }
          losatTiming.blastWriteMs += getNow() - blastWriteStartedAt;
          if (legacyPromotionTransaction.length > 0) {
            commitProteinMigration = () => {
              let nextLegacyEnvelope = legacyProteinRawCandidates.value;
              const promotedProteinIdMap = {};
              legacyPromotionTransaction.forEach(({
                candidateIndex,
                proteinIdMap
              }) => {
                nextLegacyEnvelope = transitionLegacyProteinCandidate(
                  nextLegacyEnvelope,
                  candidateIndex,
                  'promoted'
                );
                Object.entries(proteinIdMap || {}).forEach(([oldId, runtimeHandle]) => {
                  const previous = promotedProteinIdMap[oldId];
                  if (previous && previous !== runtimeHandle) {
                    throw new Error(`Protein reference '${oldId}' migrated ambiguously.`);
                  }
                  promotedProteinIdMap[oldId] = runtimeHandle;
                });
              });
              legacyProteinRawCandidates.value = nextLegacyEnvelope;
              if (!nextLegacyEnvelope.entries.some((candidate) => candidate.state === 'pending')) {
                legacyProteinDerivedEvidence.value = { schema: 1, entries: [] };
              }
              if (Object.keys(promotedProteinIdMap).length > 0) {
                orthogroups.value = rewriteMappedProteinReferences(
                  orthogroups.value,
                  promotedProteinIdMap
                );
                selectedOrthogroupAlignmentFeature.value = rewriteMappedProteinReferences(
                  selectedOrthogroupAlignmentFeature.value,
                  promotedProteinIdMap
                );
                extractedFeatures.value = rewriteMappedProteinReferences(
                  extractedFeatures.value,
                  promotedProteinIdMap
                );
                if (biologicalFeatures) {
                  biologicalFeatures.value = rewriteMappedProteinReferences(
                    biologicalFeatures.value,
                    promotedProteinIdMap
                  );
                }
              }
              legacyPromotionCommitted = true;
            };
          }
          structuredLosatTelemetry = {
            schema: 1,
            totalPairs: losatTiming.totalPairs,
            cacheHits: losatTiming.cacheHits,
            cacheMisses: losatTiming.cacheMisses,
            uniqueJobs: losatTiming.uniqueJobs,
            workerCalls: losatTiming.uniqueJobs,
            proteinDerivedPayloadCacheHits: losatTiming.proteinDerivedPayloadCacheHits,
            proteinDerivedPayloadCacheMisses: losatTiming.proteinDerivedPayloadCacheMisses,
            mode: useProteinBlastp ? blastpMode : null,
            candidateLimitRequested: useProteinBlastp ? blastpCandidateLimit : null,
            candidateLimitEffective: useProteinBlastp ? blastpCandidateLimit : null,
            collinearSearchScope: useCollinearBlastp ? collinearSearchScope : null,
            program: useProteinBlastp ? 'blastp' : losatProgram.value,
            outfmt: String(losat.outfmt || '6'),
            rawTsvEntryCount: losatTiming.rawTsvEntryCount,
            rawTsvBytes: losatTiming.rawTsvBytes,
            rawTsvLargestEntryBytes: losatTiming.rawTsvLargestEntryBytes,
            simultaneousParsedTables: losatTiming.simultaneousParsedTables,
            helperRequestMetadataBytes: losatTiming.helperRequestMetadataBytes,
            helperRequestRawTransferBytes: losatTiming.helperRequestRawTransferBytes,
            helperRequestFileCount: losatTiming.helperRequestFileCount,
            rawIdentitySchema: useProteinBlastp
              ? PROTEIN_LOSAT_CACHE_SCHEMA
              : NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
            rawJobs: cloneJsonData(losatTiming.rawJobs)
          };
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
              `Worker FASTA=${losatTiming.fastaWorkerFallbacks}`,
              `derived payload cache hits=${losatTiming.proteinDerivedPayloadCacheHits}`,
              `derived payload cache misses=${losatTiming.proteinDerivedPayloadCacheMisses}`,
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
        }
        if (useLosat) {
          pendingLosatCacheCommit = {
            cacheInfo: cacheInfo.sort(
              (left, right) => Number(left?.ordinal) - Number(right?.ordinal)
            ),
            cacheMap,
            derivedCacheMap
          };
          recordSessionLifecycleEvent('losat-cache-preparation-end', {
            cacheHits: Number(losatTiming?.cacheHits || 0),
            cacheMisses: Number(losatTiming?.cacheMisses || 0),
            derivedPayloadCacheHits: Number(
              losatTiming?.proteinDerivedPayloadCacheHits || 0
            ),
            derivedPayloadCacheMisses: Number(
              losatTiming?.proteinDerivedPayloadCacheMisses || 0
            )
          });
        }
        if ((!useLinearTrackSlots && form.show_depth) || linearSlotNeedsDepth) {
          const depthRows = linearSeqs.map((seq) => depthFileSlotsFromValue(seq.depth));
          const totalDepthFiles = depthRows.reduce((sum, row) => sum + row.filter(Boolean).length, 0);
          if (totalDepthFiles === 0) {
            throw new Error('Please upload at least one Depth TSV file or disable the depth track.');
          }
          const maxDepthTracks = depthTrackMatrixWidth(depthRows);
          for (let depthTrackIndex = 0; depthTrackIndex < maxDepthTracks; depthTrackIndex += 1) {
            const entries = depthRows.map((row, idx) => ({ file: row[depthTrackIndex] || null, idx }));
            const presentEntries = entries.filter((entry) => Boolean(entry.file));
            if (presentEntries.length === 0) {
              throw new Error(
                `Depth series #${depthTrackIndex + 1} (logical track index ${depthTrackIndex}) has no TSV source in any record. Add a TSV or remove the series.`
              );
            }
          }
          validateDepthStyleSettings();
          adv.depth_height = normalizePositiveNumberOrNull(adv.depth_height);
        }
      }

      if (annotationSets.length > 0) {
        stageTextFile('/web_annotations.tsv', encodeAnnotationTable(annotationSets), {
          name: 'annotations.tsv',
          slot: 'generatedFiles.web_annotations'
        });
      }

      throwIfGenerationCanceled();
      setProcessingStatus('Rendering SVG...');
      if (typeof serializeCanonicalFiles !== 'function') {
        throw new Error('Canonical input serialization is unavailable.');
      }
      if (!isReflow) recordSessionLifecycleEvent('generation-input-resolution-end');
      recordSessionLifecycleEvent('serialize-canonical-files-start');
      const serializedFiles = await serializeCanonicalFiles(
        activeComparisonPlanSnapshot,
        linearRecordCatalog
      );
      recordSessionLifecycleEvent('serialize-canonical-files-end');
      throwIfGenerationCanceled();
      const candidateFiles = forceEmptyComparison
        ? { ...serializedFiles, linearCanonicalComparisons: [] }
        : serializedFiles;
      const canonicalCircularConservation = resolvedCircularConservation.map((entry) => ({
        ...entry,
        fasta: candidateFiles.c_conservation_fastas?.[entry.sourceIndex] || null
      }));
      recordSessionLifecycleEvent('canonical-request-construction-start');
      const canonical = buildCanonicalRenderRequest({
        state,
        filesData: candidateFiles,
        comparisonPlanSnapshot: activeComparisonPlanSnapshot,
        resolvedComparisons,
        resolvedCircularConservation: canonicalCircularConservation
      });
      if (useCommittedComparison) {
        if (typeof getCommittedCanonicalSession !== 'function') {
          throw new Error('The preserved comparison owner is unavailable.');
        }
        inheritCommittedComparisonIntent({
          candidate: canonical,
          committed: getCommittedCanonicalSession()
        });
      }
      recordSessionLifecycleEvent('canonical-request-construction-end');
      const canonicalResourceEntries = Object.values(canonical.resources || {});
      const canonicalResourceCount = canonicalResourceEntries.length;
      const canonicalResourceDeclaredBytes = canonicalResourceEntries.reduce(
        (total, resource) => total + (Number(resource?.size) || 0),
        0
      );
      const canonicalResourceBase64Characters = canonicalResourceEntries.reduce(
        (total, resource) => total + (
          resource?.encoding === 'base64' && typeof resource?.data === 'string'
            ? resource.data.length
            : 0
        ),
        0
      );
      recordSessionLifecycleEvent('canonical-request-built');
      recordSessionLifecycleEvent('canonical-resource-count', {
        value: canonicalResourceCount
      });
      recordSessionLifecycleEvent('canonical-resource-declared-bytes', {
        value: canonicalResourceDeclaredBytes
      });
      recordSessionLifecycleEvent('canonical-resource-base64-characters', {
        value: canonicalResourceBase64Characters
      });
      recordStructuralMetric('canonicalResourceCount', canonicalResourceCount);
      recordStructuralMetric('canonicalResourceDeclaredBytes', canonicalResourceDeclaredBytes);
      recordStructuralMetric(
        'canonicalResourceBase64Characters',
        canonicalResourceBase64Characters
      );
      if (!Number.isInteger(canonicalSessionVersion)) {
        throw new Error('Canonical session version is unavailable.');
      }
      const canonicalReplayPath = '/canonical-render-session.gbdraw-session.json';
      const canonicalReplayName = makeSafeFilename(
        `${normalizedOutputPrefix || 'out'}.gbdraw-session.json`
      );
      const buildCanonicalReplayText = () => {
        recordSessionLifecycleEvent('canonical-replay-json-start');
        const replayText = JSON.stringify({
          format: 'gbdraw-session',
          version: canonicalSessionVersion,
          createdAt: manualRunStartedAtIso || new Date().toISOString(),
          renderRequest: canonical.renderRequest,
          resources: canonical.resources,
          losatCache: { entries: [] },
          losatDerivedCache: { entries: [] },
          proteinIdentityManifest: emptyProteinIdentityManifest()
        });
        recordStructuralMetric('canonicalReplayFullSerializationCount');
        recordSessionLifecycleEvent('canonical-replay-json-end');
        recordSessionLifecycleEvent('canonical-replay-json-characters', {
          value: replayText.length
        });
        return replayText;
      };
      registerRunInfoFile(canonicalReplayPath, {
        name: canonicalReplayName,
        slot: 'generatedFiles.canonical_render_session',
        kind: 'generated'
      });
      recordDeferredGeneratedCliFile(canonicalReplayPath, buildCanonicalReplayText, {
        name: canonicalReplayName,
        slot: 'generatedFiles.canonical_render_session',
        retainedBytes:
          canonicalResourceDeclaredBytes
          + canonicalResourceBase64Characters * 2
          + 65_536
      });
      const gbdrawStartedAt = getNow();
      const generationResponse = await runDiagramGeneration({
        request: canonical.renderRequest,
        resources: canonical.resources
      });
      console.info(`gbdraw ${mode.value} typed request render: ${formatDuration(getNow() - gbdrawStartedAt)}.`);
      const postGbdrawTimingEntries = [];
      const res = generationResponse.results;
      const generationMetadata = (
        generationResponse.metadata &&
        typeof generationResponse.metadata === 'object' &&
        !Array.isArray(generationResponse.metadata)
      )
        ? generationResponse.metadata
        : {};
      if (generationToken !== latestGenerationToken) {
        if (!isReflow && generationAbortSignal?.aborted) {
          return finishCanceledManualRun();
        }
        return { status: 'stale' };
      }
      if (res?.error) {
        logPostGbdrawTimings(postGbdrawTimingEntries);
        if (isReflow) {
          labelReflowLastError.value = formatPythonError(res.error)?.summary || 'Auto reflow failed';
          return { status: 'error' };
        }
        await restoreCommittedArtifact();
        errorLog.value = formatPythonError(res.error);
        return { status: 'error' };
      }
      if (!Array.isArray(res)) {
        throw new Error('The diagram engine returned an invalid Result list.');
      }

      if (isReflow && requestId !== pendingReflowRequestId) {
        return { status: 'stale' };
      }

      if (isReflow) {
        skipCaptureBaseConfig.value = true;
        skipPositionReapply.value = true;
      }

      const suppressPairwiseIdentityLegend = shouldSuppressPairwiseIdentityLegend(
        activeComparisonPlanSnapshot
      );
      if (!isReflow) recordSessionLifecycleEvent('candidate-result-validation-start');
      const candidateCatalog = isReflow
        ? null
        : validateFeatureCatalog(generationMetadata.featureCatalog, res, { adopt: true });
      if (!isReflow) recordSessionLifecycleEvent('candidate-result-validation-end');
      recordSessionLifecycleEvent('result-admission-start');
      const candidateCommit = isReflow
        ? measureTiming(
            postGbdrawTimingEntries,
            'run-analysis commit sanitized reflow results',
            () => prepareReflowCommit({
              results: res,
              suppressPairwiseIdentityLegend,
              features: extractedFeatures.value,
              featureStrokeOverrides,
              legendColorOverrides,
              legendStrokeOverrides
            })
          )
        : measureTiming(
            postGbdrawTimingEntries,
            'run-analysis sanitize and reapply editor overrides',
            () => prepareCandidateCommit({
              results: res,
              catalog: candidateCatalog,
              mode: mode.value,
              featureColorOverrides,
              featureStrokeOverrides,
              legendColorOverrides,
              legendStrokeOverrides,
              manualSpecificRules,
              legacyFeatures: committedArtifactHandle?.features?.extractedFeatures || [],
              suppressPairwiseIdentityLegend
            })
          );
      recordSessionLifecycleEvent('result-admission-end');

      if (generationToken !== latestGenerationToken) {
        if (!isReflow && generationAbortSignal?.aborted) {
          return finishCanceledManualRun();
        }
        return { status: 'stale' };
      }

      if (isReflow && requestId !== pendingReflowRequestId) {
        return { status: 'stale' };
      }

      measureTiming(postGbdrawTimingEntries, 'run-analysis assign results', () => {
        if (isReflow) results.value = candidateCommit.results;
      });
      let candidateRunInfo = null;
      let candidateCliHelpers = null;
      if (!isReflow && manualRunStartedAt !== null) {
        candidateRunInfo = buildRunInfo({
          mode: mode.value,
          args: ['--session', canonicalReplayPath],
          fileMetadata: runInfoFileMap,
          elapsedMs: getNow() - manualRunStartedAt,
          resultCount: candidateCommit.results.length,
          startedAtIso: manualRunStartedAtIso
        });
        if (structuredLosatTelemetry) {
          candidateRunInfo.losatTelemetry = cloneJsonData(structuredLosatTelemetry);
        }
        candidateCliHelpers = buildLatestCliHelperFiles(
          candidateRunInfo,
          generatedCliFileMap,
          normalizedOutputPrefix || 'out'
        );
      }
      logPostGbdrawTimings(postGbdrawTimingEntries);
      if (isReflow) {
        if (res.length > 0) {
          const safeIndex = Math.max(0, Math.min(previousSelectedResultIndex, res.length - 1));
          selectedResultIndex.value = safeIndex;
        }
      }

      if (!isReflow) {
        recordSessionLifecycleEvent('preview-result-commit-start');
        results.value = candidateCommit.results;
        trackSlotResolvedGeometry.value = generationMetadata.trackSlotGeometry || null;
        if (resultGenerationKey) resultGenerationKey.value += 1;
        selectedResultIndex.value = Math.max(
          0,
          Math.min(previousSelectedResultIndex, Math.max(0, candidateCommit.results.length - 1))
        );
        featureCatalog.value = candidateCatalog;
        extractedFeatures.value = candidateCommit.featureState.extractedFeatures;
        if (biologicalFeatures) {
          biologicalFeatures.value = candidateCommit.featureState.biologicalFeatures;
        }
        featureSelectorSafetyScope.value =
          candidateCommit.featureState.featureSelectorSafetyScope;
        featureRecordIds.value = candidateCommit.featureState.featureRecordIds;
        selectedFeatureRecordIdx.value = 0;
        editableLabels.value = [];
        setOrthogroupMetadata(candidateCommit.featureState.orthogroups);
        collinearGroups.value = Array.isArray(candidateCommit.featureState.collinearGroups)
          ? candidateCommit.featureState.collinearGroups
          : [];
        matchSequenceRegistry?.reset?.(candidateCommit.featureState.sequenceSources);
        featureExtractionPending.value = false;
        featureExtractionError.value = null;
        Object.keys(featureColorOverrides).forEach((key) => delete featureColorOverrides[key]);
        Object.assign(
          featureColorOverrides,
          cloneJsonValue(candidateCommit.featureColorOverrides, {})
        );
        Object.keys(featureStrokeOverrides).forEach((key) => delete featureStrokeOverrides[key]);
        Object.assign(
          featureStrokeOverrides,
          cloneJsonValue(candidateCommit.featureStrokeOverrides, {})
        );
        setFeatureEditorStatus({
          status: extractedFeatures.value.length ? 'summary-ready' : 'idle',
          generationId: featureExtractionRequestId,
          error: null,
          summaryCount: extractedFeatures.value.length,
          detailsCacheSize: 0
        });
        generatedLegendPosition.value = form.legend;
        generatedMode.value = mode.value;
        generatedMultiRecordCanvas.value =
          mode.value === 'circular' ? Boolean(form.multi_record_canvas) : false;
        if (mode.value === 'circular') {
          generatedCircularPlotTitlePosition.value =
            normalizeCircularPlotTitlePosition(adv.plot_title_position);
        }
        lastRunInfo.value = candidateRunInfo;
        latestCliHelperFiles = cloneCliHelperFiles(candidateCliHelpers?.files);
        latestCliHelperArchiveName =
          candidateCliHelpers?.archiveName || 'out-cli-files.zip';
        resultPanelTab.value = 'preview';
        if (typeof resetPreviewViewport === 'function') {
          resetPreviewViewport({ resetZoom: true });
        } else {
          zoom.value = 1.0;
        }
        recordSessionLifecycleEvent('preview-result-commit-end');
      }
      if (!isReflow) {
        appliedPaletteName.value = String(selectedPalette?.value || appliedPaletteName.value || 'default');
        appliedPaletteColors.value = { ...currentColors.value };
        pendingPaletteName.value = '';
        pendingPaletteColors.value = {};
      }
      commitProteinMigration?.();
      if (!isReflow && typeof setGeneratedArtifactIdentity === 'function') {
        setGeneratedArtifactIdentity(generationResponse.artifactIdentity, {
          results: candidateCommit.results
        });
      }
      if (!isReflow && pendingLosatCacheCommit) {
        losatCacheInfo.value = pendingLosatCacheCommit.cacheInfo;
        losatCache.value = pendingLosatCacheCommit.cacheMap;
        if (pendingLosatCacheCommit.derivedCacheMap) {
          losatDerivedCache.value = pendingLosatCacheCommit.derivedCacheMap;
        }
      }
      if (!isReflow) {
        globalThis.__GBDRAW_LAST_LOSAT_TELEMETRY__ = cloneJsonData(
          structuredLosatTelemetry
        );
      }
      if (!isReflow && typeof adoptCanonicalRenderArtifacts === 'function') {
        adoptCanonicalRenderArtifacts(canonical, { adoptOwnedRequest: true });
      }
      if (!isReflow && !useCommittedComparison && importedComparisonIntent) {
        Object.assign(
          importedComparisonIntent,
          createImportedComparisonIntentState(),
          { disposition: IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE }
        );
      }
      return { status: 'ok' };
    } catch (e) {
      if (!legacyPromotionCommitted) {
        legacyPromotionTransaction.forEach(({ cacheMap, cacheKey }) => {
          cacheMap.delete(cacheKey);
        });
      }
      if (isDiagramGenerationCanceled(e)) {
        if (isReflow) {
          labelReflowLastError.value = null;
        } else {
          return finishCanceledManualRun();
        }
        return { status: 'canceled' };
      }
      if (!isReflow && generationToken !== latestGenerationToken) {
        return { status: 'stale' };
      }
      if (isReflow) {
        labelReflowLastError.value = formatJsError(e)?.summary || 'Auto reflow failed';
        return { status: 'error' };
      }
      await restoreCommittedArtifact();
      errorLog.value = formatJsError(e);
      return { status: 'error' };
    } finally {
      if (isReflow) {
        labelReflowProcessing.value = false;
      } else {
        if (activeLosatAbortController === generationAbortController) {
          activeLosatAbortController = null;
        }
        if (generationToken === latestGenerationToken || canceledAttemptOwnsPresentation) {
          if (!keepProcessingStatus) processingStatus.value = '';
          generationCancelRequested.value = false;
          processing.value = false;
        }
      }
    }
  };

  const runAnalysis = async (
    comparisonPlanSnapshot = null,
    generatedArtifactHandle = null,
    comparisonExecution = null
  ) => {
    const outcome = await runAnalysisInternal({
      runMode: 'manual',
      comparisonPlanSnapshot,
      generatedArtifactHandle,
      comparisonExecution
    });
    if (['error', 'canceled'].includes(outcome?.status)) {
      failedGeneratePreservedResult.value = results.value.length > 0;
    } else if (outcome?.status === 'ok') {
      failedGeneratePreservedResult.value = false;
    }
    return outcome;
  };

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
    const hydratedEntries = (await Promise.all(losatCacheInfo.value.map(async (entry, idx) => {
      const cached = cacheMap.get(entry.key);
      if (!isCurrentRawLosatCacheEntry(cached)) return null;
      return {
        entry,
        cached,
        idx,
        hydrated: await hydrateLosatDownloadText(entry.key, cached)
      };
    }))).filter(Boolean);
    const totalBytes = totalHydratedLosatExportBytes(
      hydratedEntries.map((item) => item.hydrated)
    );

    if (!confirmHydratedLosatExport(totalBytes)) return;

    for (const { entry, idx, hydrated } of hydratedEntries) {
      const filename = entry.filename || `losat_pair_${idx + 1}.tsv`;
      downloadTextFile(
        filename || 'losat.tsv',
        hydrated.text,
        'text/tab-separated-values'
      );
      await new Promise((resolve) => setTimeout(resolve, 0));
    }
  };

  const clearLosatCache = () => {
    if (losatCache.value) {
      losatCache.value.clear();
    }
    if (losatDerivedCache.value) {
      losatDerivedCache.value.clear();
    }
    proteinIdentityManifest.value = emptyProteinIdentityManifest();
    legacyProteinRawCandidates.value = { schema: 1, entries: [] };
    legacyProteinDerivedEvidence.value = { schema: 1, entries: [] };
    losatCacheInfo.value = [];
  };

  return {
    runAnalysis,
    cancelRunAnalysis,
    captureGeneratedArtifactRuntimeState,
    restoreGeneratedArtifactRuntimeState,
    runLabelReflow,
    refreshCircularRecordOrder,
    downloadCliHelperFiles,
    downloadLosatCache,
    downloadLosatPair,
    setLosatPairFilename,
    clearLosatCache,
    getLosatPairDefaultName
  };
};
