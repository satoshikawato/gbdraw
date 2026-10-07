// @ts-check
/** @import { DrawingState } from '../state.js' */
/** @import { RulePreparation } from './rule-matching.js' */
/** @import { PreviousResultRestoreOptions, ReadinessExpectation, ReadinessExpectationOptions, ReadyReceipt } from './preview-runtime.js' */
import { validateAnnotationWarnings } from '../services/session-feature-metadata.js';
import { validateComparisonWarnings } from '../services/comparison-warnings.js';
import { keepResultNames } from '../services/result-normalization.js';
import {
  nameFeaturePlacementFailure,
  requestFeatureOverrides,
  restorePlacements,
  validateFeatureIdentityNotices
} from '../services/feature-placement.js';
import { rekeyOrthogroupOverrides } from '../services/orthogroup-feature-metadata.js';
import { resolveLinearRegionBounds } from './feature-metadata-extraction.js';
import { buildSimilarityAlignmentResetReceipt, validateSimilarityAlignmentResetReceipt } from '../services/session-active-config-contract.js';
import { prepareLosatRuntime, runLosatPairsParallel } from '../services/losat.js';
import { losatRecordGencode, prepareLosatSourceBatches, splitLosatSourceResult } from './linear-sources.js';
import {
  cancelDiagramGeneration,
  DIAGRAM_HELPER_OPERATIONS,
  DiagramGenerationCanceledError,
  isDiagramGenerationCanceled,
  runDiagramGeneration,
  runDiagramHelperOperation
} from '../services/diagram-generation.js';
import {
  projectCompositionRecordIdentity,
  buildCanonicalRenderRequest,
  bindCanonicalTypedResource,
  committedFeatureVisibilityMatches,
  projectCommittedEditorIntent,
  projectCommittedRecordTransform,
  projectCommittedSimilarityAlignment,
  promoteCanonicalRenderRequestToCurrent,
  readCanonicalResourceRecordCount,
  requestLabelProjection,
  requestLabelTableTsv
} from '../services/session-request.js';
import { labelDrawingBlocker } from './feature-editor/label-actions.js';
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
} from '../services/conservation-series.js';
import {
  applyLinearTrackOrderPlacements,
  clampLinearTrackAxisIndex,
  normalizeLinearTrackSlots,
  resolveLinearTrackAxisIndex
} from './linear-track-slots.js';
import { getDepthTrackFallbackLabel } from './depth-tracks.js';
import {
  activeDepthTrackIndices,
  depthFileSlotsFromValue,
  depthTrackCoverageCount,
  depthTrackMatrixWidth,
  isRecordMajorDepthFileMatrix,
  normalizeRecordMajorDepthFileRows,
  depthSeriesLegendCaptions,
  representativeDepthFiles,
  syncDepthSlotLabels
} from '../services/depth-track-state.js';
import { encodeAnnotationTable } from './annotations/table-codec.js';
import {
  CustomTrackPlanValidationError,
  customTrackPlanIssues,
  validateCustomTrackPlan
} from '../services/track-slot-validation.js';
import { buildRunInfo, buildSourceRecipe, summarizeLosatRuntimes } from './run-info.js';
import {
  buildLosatJobSpecs,
  resolveLinearComparisonPlan
} from '../services/linear-comparisons.js';
import {
  buildDefaultColorOverrideTsv,
  normalizePaletteColors
} from '../utils/color-utils.js';
import {
  serializeLabelWhitelistRules,
  serializeQualifierPriorityRules,
  serializeSpecificRules
} from '../services/file-imports.js';
import { rebindRuleColorOverrides } from './rule-matching.js';
import {
  pruneUnmatchedFeatureOverrides,
  serializeFeatureVisibilityRules
} from '../services/feature-visibility.js';
import {
  normalizeDefinitionLineStyleState
} from '../services/definition-line-style-state.js';
import { requireLinearLabelVisibilityMode } from '../services/linear-label-visibility.js';
import { createZipBlob } from '../utils/zip.js';
import { classifyOptionalPositiveNumber } from '../utils/optional-positive-number.js';
import { cloneJsonData, cloneJsonValue } from '../services/json-clone.js';
import { bytesToBase64, bytesToText } from '../services/byte-utils.js';
import { downloadBlob, downloadTextFile } from '../services/text-download.js';
import {
  normalizeCircularPlotTitlePosition,
  normalizeLinearPlotTitlePosition
} from '../services/layout-preferences.js';
import {
  normalizeCurrentPairwiseMatchStyle,
  requireCurrentCircularMultiRecordSizeMode,
  requireCurrentCollinearAnchorMode,
  requireCurrentCollinearColorMode,
  requireCurrentCollinearMaxConflicts,
  requireCurrentCollinearMaxDiagonalDrift,
  requireCurrentCollinearMaxParalogLinks,
  requireCurrentCollinearMaxUnitGap,
  requireCurrentCollinearMergeOrientation,
  requireCurrentCollinearMinAnchors,
  requireCurrentCollinearInferOrthogroups,
  requireCurrentCollinearSearchScope,
  requireCurrentCollinearUnitMode,
  requireCurrentLinearLabelPlacement,
  requireCurrentLinearTrackLayout,
  requireCurrentOrthogroupMemberMaxHits,
  requireCurrentOrthogroupMembershipMode,
  requireCurrentProteinBlastpCandidateLimit,
  requireCurrentProteinBlastpMaxHits,
  requireCurrentProteinBlastpMode
} from '../services/current-option-values.js';
import {
  circularDiscoveryForInput,
  discoverGffFastaRecords,
  discoverSequenceRecords,
  discoveryErrorIsFinal
} from './record-discovery.js';
import { genbankHeaderIds } from '../services/genbank-header.js';
import {
  LOSAT_DERIVED_CACHE_SCHEMA,
  NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
  PROTEIN_LOSAT_CACHE_SCHEMA,
  buildValidatedProteinIdentityIndex,
  classifyRawLosatCacheEntry,
  emptyProteinIdentityManifest,
  getCurrentRawLosatCacheEntry,
  isCurrentRawLosatCacheEntry,
  isLosatDerivedCacheEntry,
  mergeProteinIdentityManifests,
  normalizeLosatArgs,
  releaseValidatedProteinIdentityIndex,
  sameLosatArgs,
  transitionLegacyProteinCandidate,
  losatEdgeFilename,
  validateDerivedProteinReferences,
  webLosatRuntimeRecord
} from './losat-cache.js';
import { comparisonFiltersForMode, resolveComparisonThresholds } from '../mode-profiles.js';
import { diagnosticError, liveEditFailure, normalizeUserFacingError } from '../utils/error-normalization.js';
import {
  cloneFileBytesForTransfer,
  readFileBytes,
  readFileText
} from '../services/file-content-cache.js';
import {
  admitFeatureCatalog
} from '../services/feature-catalog.js';
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

const DEFAULT_LINEAR_BLAST_FILTERS = Object.freeze(
  comparisonFiltersForMode('linear')
);

export const afterPaint = ({ requestFrame = globalThis.requestAnimationFrame } = {}) => {
  if (typeof requestFrame !== 'function') return Promise.resolve();
  return new Promise((resolve) => {
    requestFrame(() => requestFrame(resolve));
  });
};

export const afterFrame = ({ requestFrame = globalThis.requestAnimationFrame } = {}) => {
  if (typeof requestFrame !== 'function') return Promise.resolve();
  return new Promise((resolve) => requestFrame(resolve));
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

// Label visibility On draws the label whatever Show Labels and the label
// filters select, so it must bind, unless the request cannot draw it: a feature
// drawn as underlay, or Label Rendering = Embedded Only. That On takes effect
// when the label can be drawn (Owner Q1, Q2). The operations reach only the
// features the request's catalog lists as rendered, which Python draws
// (should_render_feature), so a hidden feature has none. Label text alone and
// Off follow those settings, which may leave the feature unlabeled.
const labelOperationIds = (operations) => (operations || [])
  .map((operation) => String(operation?.renderedId || '').trim())
  .filter(Boolean);
const forcedLabelFeatureIds = (operations, { features, diagramOptions }) => {
  const featuresById = new Map((features || []).map((feature) => [String(feature?.svg_id || '').trim(), feature]));
  return Object.freeze([...new Set(labelOperationIds(
    (operations?.labelVisibility || []).filter((operation) => operation?.mode === 'on')
  ))].filter((featureId) => !labelDrawingBlocker(featuresById.get(featureId), diagramOptions)));
};

const getNow = () => (typeof globalThis.performance?.now === 'function' ? performance.now() : Date.now());
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

/**
 * @typedef {object} LosatCacheMetadata
 * @property {string} [identityKind] `protein` for a BLASTP search, otherwise nucleotide.
 * @property {string} flow
 * @property {string} program
 * @property {string} [outfmt]
 * @property {string[]} args
 * @property {string} [queryCanonicalHash]
 * @property {string} [subjectCanonicalHash]
 * @property {string} [queryProteinSetHash]
 * @property {string} [subjectProteinSetHash]
 * @property {string} [queryRuntimeBindingHash]
 * @property {string} [subjectRuntimeBindingHash]
 * @property {string} [queryRecordInstanceKey]
 * @property {string} [subjectRecordInstanceKey]
 * @property {Record<string, any>} [searchContext]
 */

/** @param {LosatCacheMetadata} metadata */
export const buildLosatCachePayload = ({
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
  subjectRecordInstanceKey,
  searchContext
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
      subjectRecordInstanceKey: String(subjectRecordInstanceKey || ''),
      ...(searchContext ? { searchContext } : {})
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
  if (searchContext) payload.searchContext = searchContext;
  return payload;
};

/**
 * @param {Record<string, any> | null} [manifest]
 * @param {object | null} [identityIndex]
 */
const getRawLosatCacheEntry = (cacheMap, cacheKey, metadata, manifest = null, identityIndex = null) => {
  if (!cacheMap) return null;
  const direct = getCurrentRawLosatCacheEntry(cacheMap, cacheKey, metadata, manifest, { identityIndex });
  if (direct) return direct;
  if (metadata?.identityKind === 'protein') return null;

  const expectedProgram = String(metadata?.program || '');
  const expectedOutfmt = String(metadata?.outfmt || '6');
  const expectedFlow = String(metadata?.flow || '');
  const expectedQuery = String(metadata?.queryCanonicalHash || '');
  const expectedSubject = String(metadata?.subjectCanonicalHash || '');

  for (const [key, entry] of cacheMap.entries()) {
    if (classifyRawLosatCacheEntry(entry) !== 'nucleotide-current') continue;
    if ((entry.searchContext ?? null) !== (metadata.searchContext ?? null)) continue;
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
    cacheMap.set(cacheKey, found.entry);
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

/**
 * The option values a derived protein payload depends on, as the options form holds them.
 * @typedef {object} LosatDerivedCacheMetadata
 * @property {any} [mode]
 * @property {any} [maxHits]
 * @property {any} [bitscore]
 * @property {any} [evalue]
 * @property {any} [identity]
 * @property {any} [alignmentLength]
 * @property {any} [collinearMinAnchors]
 * @property {any} [collinearMaxUnitGap]
 * @property {any} [collinearUnitMode]
 * @property {any} [collinearColorMode]
 * @property {any} [collinearAnchorMode]
 * @property {any} [collinearMergeOrientation]
 * @property {any} [collinearMaxDiagonalDrift]
 * @property {any} [collinearMaxConflictsInMergeGap]
 * @property {any} [collinearMaxParalogLinksPerOrthogroup]
 * @property {any} [collinearSearchScope]
 * @property {any} [orthogroupMembershipMode]
 * @property {any} [orthogroupMemberMaxHits]
 * @property {boolean} [collinearInferOrthogroups]
 * @property {boolean} [explicitDisplayPairs]
 * @property {Record<string, any>[]} [recordPayloads]
 * @property {Record<string, any>[]} [pairPayloads]
 */

/** @param {LosatDerivedCacheMetadata} metadata */
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
  collinearInferOrthogroups = true,
  orthogroupMembershipMode,
  orthogroupMemberMaxHits,
  explicitDisplayPairs = false,
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
        cacheKey: String(pair?.cacheKey || ''),
        // The converter returns the pairs flagged here, in this direction.
        displayPair: pair?.displayPair === true
      }))
  };
  if (normalizedMode === 'pairwise') {
    payload.pairwise = { maxHits: Number(maxHits) || 5 };
  }
  if (['orthogroup', 'collinear'].includes(normalizedMode)) {
    payload.pathRepresentation = 'lossless-graph-v1';
    payload.orthogroup = {
      membershipMode: String(orthogroupMembershipMode || 'anchor_core_v1'),
      memberMaxHits: requireCurrentOrthogroupMemberMaxHits(orthogroupMemberMaxHits)
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
      inferOrthogroups: requireCurrentCollinearInferOrthogroups(collinearInferOrthogroups),
      searchScope: String(collinearSearchScope || 'adjacent'),
      // With the search scope 'all' the CLI grid row layout limits the output to the displayed pairs.
      explicitDisplayPairs: explicitDisplayPairs === true
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
  const expectedType = normalizedMode === 'collinear' ? 'CollinearityResult' : 'OrthogroupGraphResult';
  return Boolean(
    resource
    && typeof resource === 'object'
    && !Array.isArray(resource)
    && resource.schema === 3
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
  committedSession,
  sequences,
  active
}) => {
  const committedRequest = committedSession?.renderRequest || null;
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
  // Selected proteins follow the Feature visibility rules (CO-02).
  if (!committedFeatureVisibilityMatches(committedSession, active.featureVisibility, active.featureOverrides)) {
    return false;
  }

  // Derived rows carry view coordinates and feature IDs. Raw LOSATP evidence
  // remains reusable, but a reversed view needs these rows to be projected again.
  if (committedRequest.records?.length !== sequences.length || sequences.some((seq, index) => {
    const record = committedRequest.records[index];
    const start = Number(seq.region_start) || null;
    const end = Number(seq.region_end) || null;
    const reverse = Boolean(seq.region_reverse)
      || (start !== null && end !== null && start > end);
    return record.recordKey !== seq.uid
      || String((record.region?.selector || record.selector)?.value || '') !== String(seq.region_record_id || '').trim()
      || (record.region?.start ?? null) !== (start !== null && end !== null ? Math.min(start, end) : start)
      || (record.region?.end ?? null) !== (start !== null && end !== null ? Math.max(start, end) : end)
      || reverse !== Boolean(record.region?.reverseComplement || record.presentation?.reverseComplement);
  })) return false;

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
    || (settings.orthogroupMemberMaxHits ?? null) !== active.memberMaxHits
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
    && (settings.collinearInferOrthogroups ?? true) === active.inferOrthogroups
    && String(settings.collinearitySearchScope || 'adjacent') === active.searchScope
    && String(settings.collinearityColorMode || 'orientation') === active.colorMode
    && sameNumber(
      settings.collinearMaxParalogLinksPerOrthogroup ?? 2,
      active.maxParalogLinks
    )
  );
};

const stripRuntimeCacheStats = (payload) => {
  // Canonical comparison data is read-only across rendering and cache reuse.
  // Only the transient statistics owner differs; do not duplicate the full result.
  const { cache: _cache, ...stored } = payload;
  return stored;
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
// The deferred replay builder outlives its Generate through the CLI helper
// files. Build it outside that scope so it retains only the published artifact,
// not the pre-Generate rollback handle and, through it, every earlier artifact.
const createCanonicalReplayTextBuilder = ({
  version,
  startedAtIso,
  renderRequest,
  resources,
  publishedArtifact
}) => () => {
  recordSessionLifecycleEvent('canonical-replay-json-start');
  const replayText = JSON.stringify({
    format: 'gbdraw-session',
    version,
    createdAt: startedAtIso || new Date().toISOString(),
    renderRequest,
    resources,
    // Export only after successful generation; reuse the published artifact.
    results: publishedArtifact.results,
    editorState: { featureCatalog: publishedArtifact.featureCatalog },
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
  let current = /** @type {{ id: string, parts: string[] } | null} */ (null);
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
    const header = originMatch ? genbankHeaderIds(chunk) : null;
    if (!header || !originMatch) return;
    const sequence = originMatch[1].replace(/[^A-Za-z]/g, '').toUpperCase();
    if (sequence) records.push({ id: header.recordId, sequence });
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
export const extractLosatFastaFast = async ({ file, text, fmt, regionSpec, recordSelector, reverseFlag }) => {
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
/** @param {{ file?: any, text?: string | null, fmt: string }} options */
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

/**
 * @typedef {object} CanonicalRenderCandidateExecutionOptions
 * @property {Record<string, any>} canonical The candidate's request and resources.
 * @property {string} mode
 * @property {string} [kind]
 * @property {((progress: any) => void) | null} [onProgress]
 * @property {() => boolean} [shouldAdmit]
 * @property {((payload: Record<string, any>, options?: { onProgress?: ((progress: any) => void) | null }) => Promise<any>) | null} [generationExecutor] Test seam.
 * @property {typeof admitFeatureCatalog} [catalogAdmission]
 * @property {typeof prepareCandidateRenderCommit} [prepareCommit]
 * @property {Record<string, any>} [prepareCommitInput]
 * @property {((canonical: Record<string, any>, catalogState: any) => any) | null} [decorationContinuity]
 * @property {any[]} [timingEntries]
 * @property {readonly string[] | null} [resultNames] The names of the Results a rerender draws again,
 *   which they keep (OV-136); Generate's Results take the names the engine gives them.
 */

/** @param {CanonicalRenderCandidateExecutionOptions} options */
export const executeCanonicalRenderCandidate = async ({
  canonical,
  mode,
  kind = 'generate',
  onProgress = null,
  shouldAdmit = () => true,
  generationExecutor = null,
  catalogAdmission = admitFeatureCatalog,
  prepareCommit = prepareCandidateRenderCommit,
  prepareCommitInput = {},
  decorationContinuity = null,
  timingEntries = [],
  resultNames = null
}) => {
  if (!canonical?.renderRequest || !canonical?.resources) {
    throw new Error('Canonical candidate execution requires a request and resources.');
  }
  recordStructuralMetric('canonicalCandidateExecutionCount', 1, { kind });
  const startedAt = getNow();
  const generationResponse = generationExecutor
    ? await generationExecutor({
        request: canonical.renderRequest,
        resources: canonical.resources
      }, { onProgress })
    : await runDiagramGeneration({
        request: canonical.renderRequest,
        resources: canonical.resources
      }, { onProgress: onProgress });
  const results = generationResponse.results;
  if (!shouldAdmit()) {
    return { status: /** @type {const} */ ('superseded'), generationResponse, elapsedMs: getNow() - startedAt };
  }
  if (results?.error) {
    return {
      status: /** @type {const} */ ('engine-error'),
      generationResponse,
      engineError: results.error,
      elapsedMs: getNow() - startedAt
    };
  }
  if (!Array.isArray(results)) {
    throw new Error('The diagram engine returned an invalid Result list.');
  }
  keepResultNames(generationResponse, resultNames);
  const metadata = generationResponse.metadata
    && typeof generationResponse.metadata === 'object'
    && !Array.isArray(generationResponse.metadata)
    ? generationResponse.metadata
    : {};
  const annotationWarnings = validateAnnotationWarnings(metadata.annotationWarnings, results);
  const comparisonWarnings = validateComparisonWarnings(metadata.comparisonWarnings, results);
  // Edits by source identity that a Result does not draw (design Q4 3.4).
  const featureIdentityNotices = Array.isArray(metadata.featureIdentityNotices)
    ? metadata.featureIdentityNotices : [];
  validateFeatureIdentityNotices(featureIdentityNotices, results);
  recordSessionLifecycleEvent('candidate-result-validation-start');
  const catalogState = catalogAdmission(metadata.featureCatalog, results, {
    adopt: true,
    mode
  });
  recordSessionLifecycleEvent('candidate-result-validation-end');
  recordSessionLifecycleEvent('result-admission-start');
  const commit = measureTiming(
    timingEntries,
    kind === 'reflow'
      ? 'run-analysis commit sanitized reflow results'
      : 'run-analysis sanitize and reapply editor overrides',
    () => prepareCommit({
      generationResponse,
      catalogAdmission: catalogState,
      // The Results carry the feature types of the request that drew them (R13).
      selectedFeatureTypes: canonical.renderRequest.diagramOptions?.selectedFeaturesSet ?? null,
      ...prepareCommitInput,
      resultTransforms: decorationContinuity?.(canonical, catalogState) || []
    })
  );
  recordSessionLifecycleEvent('result-admission-end');
  return {
    status: /** @type {const} */ ('ok'),
    generationResponse,
    generationMetadata: metadata,
    annotationWarnings,
    featureIdentityNotices,
    comparisonWarnings,
    results,
    catalogAdmission: catalogState,
    catalog: catalogState.catalog,
    commit,
    elapsedMs: getNow() - startedAt
  };
};

/**
 * The members of the preview owner's runtime (app/preview-runtime.js) that
 * Generate reads: it registers the readiness a candidate must meet, and rolls
 * the selection back when the candidate is rejected.
 * @typedef {object} RunAnalysisPreviewRuntime
 * @property {(expectation: ReadinessExpectationOptions) => ReadinessExpectation} registerReadinessExpectation
 * @property {(generationToken: string, reason: Error) => void} invalidateReadinessExpectation
 * @property {(receipt: any, reason: string) => void} invalidateReadyReceipt
 * @property {(options: PreviousResultRestoreOptions) => Promise<any>} restorePreviousSelectedResult
 * @property {(index: number) => any} selectResult
 * @property {(result: any) => any} getResultIdentity
 * @property {() => ({ readyReceipt?: any } | null)} [getActiveRuntime]
 */

/**
 * @typedef {object} RunAnalysisOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {RulePreparation | null} [rulePreparation] The rule owner's preparation, which Generate runs before it draws.
 * @property {(feature: Record<string, any>) => boolean} isCurrentFeature Whether the feature belongs to the displayed Result.
 * @property {(comparisonPlanSnapshot: Record<string, any> | null, linearRecordCatalog: any, drawing: DrawingState) => any} serializeCanonicalFiles
 *   The Session owner's serialization of the active render files, chosen by the run's drawing.
 * @property {number} canonicalSessionVersion `SESSION_VERSION` of the Session owner.
 * @property {(canonical: Record<string, any>, options: { adoptOwnedRequest?: boolean }) => void} adoptCanonicalRenderArtifacts
 *   The Session owner's adoption of a candidate's request, resources, and files.
 * @property {(() => (Record<string, any> | null)) | null} [getCommittedCanonicalSession] The committed canonical Session.
 * @property {(session: Record<string, any> | null | undefined, projectIdentity: (...args: any[]) => any, drawing: Record<string, any>) => any} [captureDecorationContinuity]
 *   The Legend layout owner's capture of the decoration positions a rerender of the drawing keeps.
 * @property {() => Record<string, any> | Promise<Record<string, any>>} captureGeneratedArtifactHandle History's capture of the current generated artifact.
 * @property {() => ({ results: any[] } & Record<string, any>)} captureGeneratedArtifactOwnerSet
 * @property {(ownerSet: Record<string, any>, options: { selectedResultIndex: number, installResults: (results: any[]) => void }) => void} installGeneratedArtifactOwnerSet
 * @property {(handle: Record<string, any>) => any} restoreGeneratedArtifactHandle
 * @property {((artifactIdentity: any, options: { results: any[] }) => void) | null} [setGeneratedArtifactIdentity]
 * @property {((label: string, execute: (handle?: Record<string, any>) => any, options: Record<string, any>) => Promise<any>) | null} [runGeneratedArtifactReplacement]
 *   History's undoable replacement of the generated artifact.
 * @property {RunAnalysisPreviewRuntime | null} [previewRuntime]
 * @property {() => Promise<void>} [nextTick] Vue `nextTick`.
 * @property {() => Promise<void>} [waitForAfterPaint]
 * @property {() => Promise<void>} [waitForPostBindFrame]
 * @property {((capture: { phase?: string, diagnostics?: any }) => void) | null} [onGeneratedArtifactCheckpointCapture]
 * @property {(options?: { pan?: any, resetZoom?: boolean }) => void} resetPreviewViewport The preview owner's viewport reset.
 * @property {(() => ({ code: string, context?: any } | null | undefined)) | null} [validateAnnotationTargets]
 * @property {(() => Promise<Record<string, any>>) | null} [prepareLinearRecordCatalog]
 * @property {{ value: Record<string, any>[] } | null} [recordDisplayRows]
 *   The draft record display rows the request reads (`recordDisplayControls.allRows`, R13).
 * @property {() => Promise<void>} [settleComparisonRecordLabels] D12: ring rows added just before Generate are named before it reads them.
 * @property {((mode: string, state: any) => void) | null} [assertActiveModeInputs]
 *   services/config.js owns the active-mode input check shared with Save.
 * @property {() => void} closeLabelTextScopeDialog The label owner's port (app/feature-editor/label-actions.js).
 * @property {(options?: { rerender?: boolean }) => void} clearLabelBuildNotices The label owner's port.
 * @property {typeof prepareCandidateRenderCommit} [prepareCandidateCommit] Test seam.
 */

/** @param {RunAnalysisOptions} options */
export const createRunAnalysis = ({
  state,
  rulePreparation = null,
  isCurrentFeature,
  serializeCanonicalFiles,
  canonicalSessionVersion,
  adoptCanonicalRenderArtifacts,
  getCommittedCanonicalSession = null,
  captureDecorationContinuity = () => null,
  captureGeneratedArtifactHandle,
  captureGeneratedArtifactOwnerSet,
  installGeneratedArtifactOwnerSet,
  restoreGeneratedArtifactHandle,
  setGeneratedArtifactIdentity = null,
  runGeneratedArtifactReplacement = null,
  previewRuntime = null,
  nextTick = () => Promise.resolve(),
  waitForAfterPaint = afterPaint,
  waitForPostBindFrame = afterFrame,
  onGeneratedArtifactCheckpointCapture = null,
  resetPreviewViewport,
  validateAnnotationTargets = null,
  prepareLinearRecordCatalog = null,
  // The draft record display rows the request reads (`recordDisplayControls.allRows`, R13).
  recordDisplayRows = null,
  // D12: ring rows added just before Generate are named before it reads them.
  settleComparisonRecordLabels = async () => {},
  // services/config.js owns the active-mode input check shared with Save.
  assertActiveModeInputs = null,
  // R13: the label owner's ports (app/feature-editor/label-actions.js).
  closeLabelTextScopeDialog,
  clearLabelBuildNotices,
  prepareCandidateCommit = prepareCandidateRenderCommit
}) => {
  const {
    processing,
    processingStatus,
    generationCancelRequested,
    results,
    selectedResultIndex,
    failedGeneratePreservedResult,
    generationFailureRecovery,
    resultGenerationKey,
    resultPanelTab,
    pairwiseMatchFactors,
    errorLog,
    semanticFileWatchersSuppressed,
    sessionImportRollbackInProgress,
    zoom,
    skipCaptureBaseConfig,
    matchSequenceRegistry,
    paletteDefinitions,
    appliedPaletteName,
    appliedPaletteColors,
    mode,
    cInputType,
    lInputType,
    losatCacheInfo,
    losatThreadingStatus,
    losatCache,
    losatDerivedCache,
    proteinIdentityManifest,
    legacyProteinRawCandidates,
    legacyProteinDerivedEvidence,
    orthogroups,
    selectedOrthogroupAlignmentFeature,
    selectedOrthogroupId,
    circularRecordList,
    circularRecordDiscovery,
    files,
    linearSeqs,
    shouldDeferCircularPreviewUpdates,
    extractedFeatures,
    biologicalFeatures,
    featureCatalog,
    featureEditRemovalCount,
    featureEditorStatus,
    featureExtractionPending,
    featureExtractionError,
    selectedFeatureRecordIdx,
    labelReflowProcessing,
    labelReflowLastError,
    originalLegendOrder
  } = state;
  if (
    typeof captureGeneratedArtifactHandle !== 'function'
    || typeof captureGeneratedArtifactOwnerSet !== 'function'
    || typeof installGeneratedArtifactOwnerSet !== 'function'
    || typeof restoreGeneratedArtifactHandle !== 'function'
  ) {
    throw new Error('createRunAnalysis requires generated artifact transaction handlers.');
  }
  if (
    !previewRuntime
    || typeof previewRuntime.registerReadinessExpectation !== 'function'
    || typeof previewRuntime.restorePreviousSelectedResult !== 'function'
    || typeof previewRuntime.selectResult !== 'function'
  ) {
    throw new Error('createRunAnalysis requires PreviewRuntime readiness and selection handlers.');
  }
  const executeLosatJobs = (...args) => {
    const override = globalThis.__GBDRAW_LOSAT_EXECUTOR__;
    return (typeof override === 'function' ? override : runLosatPairsParallel)(...args);
  };
  // Show Depth off leaves the Depth sources out of the request, so Python cannot name
  // the Depth rows a Legend style may still address; the draft names them (OV-81).
  /** @param {DrawingState} drawing */
  const unrequestedDepthCaptions = (drawing, canonical) => (
    canonical.renderRequest.diagramOptions?.depthTracks
      ? []
      : depthSeriesLegendCaptions({
          depthTracks: drawing.adv.depth_tracks,
          slots: canonical.renderRequest.mode === 'linear' ? drawing.adv.linear_track_slots : drawing.adv.circular_track_slots
        })
  );
  // A target-record transform and a label reflow draw the committed Results
  // again, with the settings and edits of the committed request's mode.
  /** @param {Record<string, any> | null | undefined} committed */
  const committedDrawing = (committed) => state.drawings[committed?.renderRequest?.mode === 'linear' ? 'linear' : 'circular'];
  let pendingReflowRequestId = 0;
  let activeReflowRequestId = 0;
  // A rerender requested while Generate (or a committed-candidate run) is in
  // progress waits for it: the rerender would take the newer generation token
  // and silently supersede the run (OV-48). It runs once after the run settles.
  let reflowDeferredByProcessing = false;
  let featureExtractionRequestId = 0;
  let latestGenerationToken = 0;
  let latestOperationId = 0;
  let circularRecordRefreshGeneration = 0;
  /** @type {{ fingerprint: any[], promise: Promise<any> | null } | null} */
  let activeCircularRecordRefresh = null;
  /** @type {AbortController | null} */
  let activeLosatAbortController = null;
  /** @type {readonly Record<string, any>[]} */
  let latestCliHelperFiles = Object.freeze([]);
  let latestCliHelperArchiveName = 'out-cli-files.zip';
  let latestCliHelperRetainedBytes = 0;
  const cloneCliHelperFiles = (files) => (
    (Array.isArray(files) ? files : []).map((file) => ({ ...file }))
  );
  const publishGeneratedArtifactRuntimeState = ({ files, archiveName, losatTelemetry }) => {
    latestCliHelperFiles = Object.freeze(cloneCliHelperFiles(files));
    latestCliHelperArchiveName = String(archiveName || 'out-cli-files.zip');
    latestCliHelperRetainedBytes = latestCliHelperFiles.reduce((total, file) => (
      total
      + String(file?.name || '').length * 2
      + String(file?.data || '').length * 2
      + Math.max(0, Number(file?.retainedBytes) || 0)
      + 256
    ), latestCliHelperArchiveName.length * 2 + 16_384);
    globalThis.__GBDRAW_LAST_LOSAT_TELEMETRY__ = losatTelemetry;
  };
  const captureGeneratedArtifactRuntimeState = () => (
    Object.freeze({
      latestCliHelperFiles,
      latestCliHelperArchiveName,
      losatTelemetry: globalThis.__GBDRAW_LAST_LOSAT_TELEMETRY__,
      retainedBytes: latestCliHelperRetainedBytes
    })
  );
  /**
   * @param {Record<string, any>} [runtimeState]
   * @param {{ ui?: { canvasPan?: any } }} [options]
   */
  const restoreGeneratedArtifactRuntimeState = (runtimeState = {}, { ui = {} } = {}) => {
    latestCliHelperFiles = runtimeState?.latestCliHelperFiles || Object.freeze([]);
    latestCliHelperArchiveName = String(
      runtimeState?.latestCliHelperArchiveName || 'out-cli-files.zip'
    );
    latestCliHelperRetainedBytes = Math.max(0, Number(runtimeState?.retainedBytes) || 0);
    globalThis.__GBDRAW_LAST_LOSAT_TELEMETRY__ = runtimeState?.losatTelemetry ?? null;
    if (typeof resetPreviewViewport === 'function') {
      resetPreviewViewport({ pan: ui?.canvasPan });
    }
  };
  const generatedArtifactTransactionOwner = Object.freeze({
    /**
     * @param {Record<string, any>} ownerSet
     * @param {{ runtimeState?: Record<string, any> | null }} [options]
     */
    build(ownerSet, { runtimeState = null } = {}) {
      recordSessionLifecycleEvent('artifact.candidate-completed');
      recordStructuralMetric('generatedArtifactCandidateBuildCount', 1);
      return Object.freeze({
        ownerSet: Object.freeze(ownerSet),
        runtimeState
      });
    },
    activate(candidate, { selectedResultIndex = 0 } = {}) {
      recordSessionLifecycleEvent('artifact.activation-started', { selectedResultIndex });
      installGeneratedArtifactOwnerSet(candidate.ownerSet, {
        selectedResultIndex,
        installResults(nextResults) {
          results.value = nextResults;
        }
      });
      publishGeneratedArtifactRuntimeState(candidate.runtimeState || {});
      recordStructuralMetric('generatedArtifactActivationCount', 1);
      recordSessionLifecycleEvent('artifact.activation-completed', { selectedResultIndex });
    },
    finalize() {
      recordSessionLifecycleEvent('artifact.finalization-started');
      recordStructuralMetric('generatedArtifactFinalizeCount', 1);
      recordSessionLifecycleEvent('artifact.finalization-completed');
    },
    async restore(handle) {
      recordStructuralMetric('generatedArtifactRollbackCount', 1);
      const currentAlert = errorLog.value;
      try {
        await restoreGeneratedArtifactHandle(handle);
      } finally {
        // Snapshot restoration clones presentation. Retain the notification
        // that owned the rollback boundary, without deriving its cause here.
        const restoredAlert = handle.mutableIntent?.presentation?.errorLog;
        if (JSON.stringify(normalizeUserFacingError(errorLog.value))
          === JSON.stringify(normalizeUserFacingError(restoredAlert))) errorLog.value = currentAlert;
      }
    }
  });
  const recordDiscoverySuppressed = () => Boolean(
    semanticFileWatchersSuppressed?.value ||
    sessionImportRollbackInProgress?.value
  );

  /** @param {DrawingState} drawing */
  const buildLatestCliHelperFiles = (drawing, runInfo, generatedCliFileMap, archiveBaseName) => {
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
    const archiveStem = makeSafeFilename(`${archiveBaseName || drawing.form.prefix || 'out'}-cli-files`);
    return {
      files,
      archiveName: `${archiveStem}.zip`
    };
  };

  const downloadCliHelperFiles = () => {
    if (!latestCliHelperFiles.length) {
      alert('No reproducibility files are available for the latest run.');
      return;
    }
    const materializedFiles = latestCliHelperFiles.map((file) => ({
      name: file.name,
      data: typeof file.buildData === 'function' ? file.buildData() : file.data
    }));
    const archiveNames = new Set();
    materializedFiles.forEach((file) => {
      const name = String(file.name || '');
      const key = name.replace(/[A-Z]/g, (letter) => letter.toLowerCase());
      if (!name || archiveNames.has(key)) {
        throw new Error('Reproducibility bundle filenames must be unique before ZIP creation.');
      }
      archiveNames.add(key);
    });
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

  const formatError = (cause, operation = 'generate', stage = 'request-validation') =>
    normalizeUserFacingError(cause || { code: 'UNKNOWN' }, { operation, stage });
  /**
   * @param {any} cause
   * @param {{ handle?: Record<string, any> | null, restore?: (() => any) | null, operation?: string, stage?: string,
   *   isCurrent?: () => boolean, isCurrentOperation?: () => boolean, recovery?: string | null }} [options]
   */
  const failOperation = async (cause, { handle, restore = null, operation = 'generate',
    stage = 'render', isCurrent = () => true, isCurrentOperation = () => true, recovery = null } = {}) => {
    if (!isCurrentOperation()) return { status: 'stale' };
    const error = formatError(cause, operation, stage);
    const previousResults = handle?.ownerSet?.results || [];
    if (!recovery) {
      try {
        const restored = restore ? await restore() : false;
        const currentResults = captureGeneratedArtifactOwnerSet().results;
        const unchanged = currentResults.length === previousResults.length
          && currentResults.every((result, index) => result === previousResults[index]);
        recovery = !handle ? (results.value.length ? 'restore-failed' : 'no-result')
          : previousResults.length === 0 ? 'no-result'
          : restored ? 'restored' : unchanged ? 'preserved' : 'restore-failed';
      } catch (_) { recovery = 'restore-failed'; }
    }
    if (!isCurrent()) return { status: 'stale' };
    errorLog.value = error;
    if (generationFailureRecovery) generationFailureRecovery.value = recovery;
    failedGeneratePreservedResult.value = ['preserved', 'restored'].includes(recovery);
    return { status: 'error', error, recovery };
  };

  // Raw searches finish before the artifact transaction. Keep only the latest
  // search's entries for retry, without changing the saved Result. Cache owner
  // replacement (Clear Cache, Session load, or History) invalidates the retry.
  /** @type {{ owner: any, entries: Map<string, any> } | null} */
  let completedLosatSearch = null;
  /**
   * @param {Record<string, any> | null} [manifest]
   * @param {object | null} [identityIndex]
   */
  const getReusableLosatCacheEntry = (cacheMap, cacheKey, metadata, manifest = null, identityIndex = null) => {
    if (completedLosatSearch?.owner !== losatCache.value) completedLosatSearch = null;
    return getRawLosatCacheEntry(cacheMap, cacheKey, metadata, manifest, identityIndex)
      || getRawLosatCacheEntry(completedLosatSearch?.entries, cacheKey, metadata, manifest, identityIndex);
  };
  const retainCompletedLosatSearch = (cacheMap, pairs) => {
    completedLosatSearch = {
      owner: losatCache.value,
      entries: new Map(pairs.map(({ cacheKey }) => [cacheKey, cacheMap.get(cacheKey)]))
    };
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

  // D-21: names follow the exact member set; the rest stay dormant.
  /** @param {DrawingState} drawing */
  const rekeyCommittedOrthogroupOverrides = (drawing, previousGroups, candidateGroups) => {
    const next = rekeyOrthogroupOverrides({
      previousGroups,
      candidateGroups,
      names: drawing.orthogroupNameOverrides,
      descriptions: drawing.orthogroupDescriptionOverrides,
      dormant: drawing.orthogroupDormantOverrides
    });
    [[drawing.orthogroupNameOverrides, next.names], [drawing.orthogroupDescriptionOverrides, next.descriptions],
      [drawing.orthogroupDormantOverrides, next.dormant]].forEach(([target, values]) => {
      Object.keys(target).forEach((key) => delete target[key]);
      Object.assign(target, values);
    });
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

  /** @param {DrawingState} drawing */
  const buildLosatSuffix = (drawing) => {
    if (drawing.losatProgram.value === 'blastn') return 'losatn';
    if (drawing.losatProgram.value === 'blastp') return 'losatp';
    return 'tlosatx';
  };

  /** @param {DrawingState} drawing */
  const buildLosatFilename = (drawing, leftLabel, rightLabel) => (
    losatEdgeFilename(leftLabel, rightLabel, buildLosatSuffix(drawing))
  );

  /** @param {DrawingState} drawing */
  const getResolvedLinearEdge = (drawing, edgeKey) => {
    const normalizedKey = String(edgeKey || '').trim();
    if (!normalizedKey) return null;
    const edges = drawing.linearComparisonResolution?.value?.edges;
    return (Array.isArray(edges) ? edges : []).find((edge) => edge.edgeKey === normalizedKey) || null;
  };

  const getLosatCacheInfoEntry = (edgeKey) => {
    const entries = Array.isArray(losatCacheInfo.value) ? losatCacheInfo.value : [];
    if (Number.isInteger(edgeKey)) return entries[edgeKey] || null;
    const normalizedKey = String(edgeKey || '').trim();
    return entries.find((entry) => String(entry?.edgeKey || '') === normalizedKey) || null;
  };

  /**
   * @param {DrawingState} drawing
   * @param {string} edgeKey
   * @param {{ recordId?: string } | null} [queryEntry]
   * @param {{ recordId?: string } | null} [subjectEntry]
   */
  const losatPairDefaultName = (drawing, edgeKey, queryEntry = null, subjectEntry = null) => {
    const cacheEntry = getLosatCacheInfoEntry(edgeKey);
    const edge = getResolvedLinearEdge(drawing, edgeKey) || cacheEntry;
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
    return buildLosatFilename(drawing, leftLabel, rightLabel);
  };
  // LOSAT pairs are the Linear comparisons: their names follow the Linear drawing.
  /**
   * @param {string} edgeKey
   * @param {{ recordId?: string } | null} [queryEntry]
   * @param {{ recordId?: string } | null} [subjectEntry]
   */
  const getLosatPairDefaultName = (edgeKey, queryEntry = null, subjectEntry = null) => {
    const drawing = state.drawings.linear;
    return losatPairDefaultName(drawing, edgeKey, queryEntry, subjectEntry);
  };

  const normalizeLosatFilename = (name, fallback) => {
    const raw = String(name || '').trim() || String(fallback || '');
    const withExt = raw.toLowerCase().endsWith('.tsv') ? raw : `${raw}.tsv`;
    return makeSafeFilename(withExt);
  };

  // How LOSAT runs is one app-level setting (`state.losatExecution`) for both
  // drawings: the Linear comparisons and the Circular conservation series.
  const getLosatParallelWorkers = () => {
    const raw = String(state.losatExecution.parallelWorkers || 'auto').trim().toLowerCase();
    if (raw === 'auto') return undefined;
    const parsed = Number(raw);
    return Number.isInteger(parsed) && parsed >= 1 ? parsed : undefined;
  };

  const getLosatExecutionMode = () => {
    const raw = String(state.losatExecution.executionMode || 'auto').trim().toLowerCase();
    return ['auto', 'serial', 'threaded'].includes(raw) ? raw : 'auto';
  };

  /** @param {DrawingState} drawing */
  const getLosatThreadsPerJob = (drawing) => {
    if (drawing.losatProgram.value !== 'blastp') return 1;
    const raw = String(state.losatExecution.threadsPerJob || 'auto').trim().toLowerCase();
    if (raw === 'auto') return undefined;
    const parsed = Number(raw);
    return Number.isInteger(parsed) && parsed >= 1 ? parsed : undefined;
  };

  const getLosatTotalThreadBudget = () => {
    const raw = String(state.losatExecution.totalThreadBudget || 'safe').trim().toLowerCase();
    if (raw === 'safe' || raw === 'auto') return undefined;
    if (raw === 'available') {
      return Math.max(1, Number(globalThis.navigator?.hardwareConcurrency || 4) || 4);
    }
    const parsed = Number(raw);
    if (!Number.isInteger(parsed) || parsed < 1) return undefined;
    const hardwareBudget = Math.max(1, Number(globalThis.navigator?.hardwareConcurrency || 4) || 4);
    return Math.min(parsed, hardwareBudget);
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
    const drawing = state.drawings.linear;
    const entry = getLosatCacheInfoEntry(edgeKey);
    const cacheMap = losatCache.value;
    if (!entry || !cacheMap) return;
    const cached = cacheMap.get(entry.key);
    if (!isCurrentRawLosatCacheEntry(cached)) return;
    const defaultName = losatPairDefaultName(drawing, entry.edgeKey || edgeKey);
    const fallbackOrdinal = Number.isInteger(Number(entry.ordinal))
      ? Number(entry.ordinal)
      : 0;
    const filename = normalizeLosatFilename(
      customName,
      entry.filename || defaultName || `losat_pair_${fallbackOrdinal + 1}.tsv`
    );
    if (!state.sessionOperationAvailability?.()) {
      losatCacheInfo.value = losatCacheInfo.value.map((candidate) => (
        candidate === entry ? { ...candidate, filename } : candidate
      ));
    }
    const hydrated = await hydrateLosatDownloadText(entry.key, cached);
    downloadTextFile(
      filename || 'losat.tsv',
      hydrated.text,
      'text/tab-separated-values'
    );
  };

  const setLosatPairFilename = (edgeKey, customName) => {
    const drawing = state.drawings.linear;
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = getLosatCacheInfoEntry(edgeKey);
    if (!entry) return;
    const defaultName = losatPairDefaultName(drawing, entry.edgeKey || edgeKey);
    const fallbackOrdinal = Number.isInteger(Number(entry.ordinal))
      ? Number(entry.ordinal)
      : 0;
    const filename = normalizeLosatFilename(
      customName,
      entry.filename || defaultName || `losat_pair_${fallbackOrdinal + 1}.tsv`
    );
    losatCacheInfo.value = losatCacheInfo.value.map((candidate) => (
      candidate === entry ? { ...candidate, filename } : candidate
    ));
  };

  const circularDiscoveryTargetsCurrentInput = () => circularDiscoveryForInput(state).current;

  const circularDiscoveryMatchesCurrentInput = () => (
    circularRecordDiscovery.status === 'ready' &&
    circularDiscoveryTargetsCurrentInput()
  );

  /** @param {DrawingState} drawing */
  const validateDepthInputPresence = (drawing) => {
    const linear = mode.value === 'linear';
    const slots = linear ? drawing.adv.linear_track_slots : drawing.adv.circular_track_slots;
    const customDepthRequested = (
      (linear ? drawing.adv.linear_track_slots_enabled : drawing.adv.circular_track_slots_enabled) === true &&
      (Array.isArray(slots) ? slots : []).some((slot) => (
        slot?.enabled !== false && String(slot?.renderer || '') === 'depth'
      ))
    );
    if (!drawing.form.show_depth && !customDepthRequested) return '';
    let rows;
    if (linear) {
      rows = linearSeqs.map((seq) => depthFileSlotsFromValue(seq.depth));
    } else {
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
      if (isRecordMajorDepthFileMatrix(files.c_depth) && files.c_depth.length !== recordCount) {
        return `Circular Depth matrix has ${files.c_depth.length} record rows; expected ${recordCount}.`;
      }
      rows = normalizeRecordMajorDepthFileRows(files.c_depth, recordCount);
    }
    if (customDepthRequested) {
      // An enabled Depth row must reference a logical series with a source in
      // some record; Generate stops on the row issue the track editor shows
      // (PD-OI-083).
      const issues = customTrackPlanIssues(validateCustomTrackPlan({
        mode: mode.value,
        slots,
        depthTrackCount: depthTrackMatrixWidth(rows),
        depthSourcedTrackIndexes: activeDepthTrackIndices(rows)
      })).filter((issue) => String(issue.code || '').startsWith('depth_'));
      if (issues.length > 0) return new CustomTrackPlanValidationError(issues);
    }
    if (!rows.some((row) => row.some(Boolean))) {
      return linear
        ? 'Please upload at least one Depth TSV file or disable Show depth track.'
        : 'Please upload a Depth TSV file or disable Show depth track.';
    }
    const logicalWidth = depthTrackMatrixWidth(rows);
    for (let trackIndex = 0; trackIndex < logicalWidth; trackIndex += 1) {
      if (depthTrackCoverageCount(rows, trackIndex) > 0) continue;
      return `Depth series #${trackIndex + 1} (logical track index ${trackIndex}) has no TSV source in any record. Add a TSV or remove the series.`;
    }
    return '';
  };

  const runCircularRecordRefresh = async ({ suppress = false, automatic = false, reuseFinalError = false } = {}) => {
    const drawing = state.drawings.circular;
    // An inactive Circular source keeps its records and Multi-Record Canvas
    // order; only a read still running for it is superseded.
    if (mode.value !== 'circular') {
      if (circularRecordDiscovery.status === 'loading') {
        circularRecordRefreshGeneration += 1;
        circularRecordDiscovery.status = 'idle';
      }
      return;
    }
    if (automatic && circularDiscoveryTargetsCurrentInput()
      && ['deferred', 'ready', 'error'].includes(circularRecordDiscovery.status)) return;
    // Generate reports a final error of these exact files without reading them
    // again; Retry source inspection still reads.
    if (reuseFinalError && circularDiscoveryTargetsCurrentInput()
      && circularRecordDiscovery.status === 'error' && discoveryErrorIsFinal(circularRecordDiscovery.error)) return;
    const refreshGeneration = ++circularRecordRefreshGeneration;
    if (suppress || recordDiscoverySuppressed()) return;
    if (!Array.isArray(drawing.adv.multi_record_positions)) {
      drawing.adv.multi_record_positions = [];
    }
    const inputType = cInputType.value;
    const primaryFile = inputType === 'gff' ? files.c_gff : files.c_gb;
    const pairedFile = inputType === 'gff' ? files.c_fasta : null;
    const hasCompleteInput = Boolean(primaryFile && (inputType !== 'gff' || pairedFile));
    const preserveCanonicalRecordKeys = circularDiscoveryTargetsCurrentInput();
    const preservedRecordKeys = new Map();
    if (preserveCanonicalRecordKeys) {
      [
        ...(Array.isArray(circularRecordDiscovery.canonicalRecordIdentities)
          ? circularRecordDiscovery.canonicalRecordIdentities
          : []),
        ...(Array.isArray(circularRecordList.value) ? circularRecordList.value : [])
      ]
        .forEach((record) => {
          const recordKey = String(record?.recordKey || '').trim();
          if (!recordKey) return;
          [record?.selector, record?.record_id, record?.recordId]
            .map((value) => String(value || '').trim())
            .filter(Boolean)
            .forEach((value) => preservedRecordKeys.set(value, recordKey));
        });
    }
    Object.assign(circularRecordDiscovery, {
      status: hasCompleteInput ? 'loading' : 'idle',
      error: '',
      inputType,
      primaryFile: primaryFile || null,
      pairedFile: pairedFile || null,
      canonicalRecordIdentities: preserveCanonicalRecordKeys
        ? circularRecordDiscovery.canonicalRecordIdentities
        : []
    });
    circularRecordList.value = [];
    if (!hasCompleteInput) {
      drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length);
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
        ? await discoverGffFastaRecords({ gffFile: primaryFile, fastaFile: pairedFile })
        : await discoverSequenceRecords({ file: primaryFile, format: 'genbank' });
      if (
        refreshGeneration !== circularRecordRefreshGeneration ||
        recordDiscoverySuppressed() ||
        mode.value !== 'circular' ||
        cInputType.value !== inputType ||
        (inputType === 'gff' ? files.c_gff : files.c_gb) !== primaryFile ||
        (inputType === 'gff' ? files.c_fasta : null) !== pairedFile
      ) return;
      const busy = state.sessionOperationAvailability?.();
      if (busy) {
        circularRecordDiscovery.status = 'deferred';
        circularRecordDiscovery.primaryFile = null;
        return busy;
      }
      const nextRecords = records.map((entry) => {
        const recordKey = preservedRecordKeys.get(String(entry.selector || '').trim())
          || preservedRecordKeys.get(String(entry.recordId || '').trim())
          || '';
        return {
          selector: entry.selector,
          record_id: entry.recordId,
          record_length: entry.recordLength,
          detectedTopology: entry.detectedTopology,
          ...(recordKey ? { recordKey } : {})
        };
      });
      circularRecordList.value = nextRecords;
      circularRecordDiscovery.status = 'ready';
      // Identity belongs to the source, including while its mode is inactive.
      circularRecordDiscovery.canonicalRecordIdentities = nextRecords
        .filter((record) => record.recordKey)
        .map(({ selector, record_id, recordKey }) => ({ selector, record_id, recordKey }));
      const nextPositions = mergeCircularRecordPositions(nextRecords, drawing.adv.multi_record_positions);
      drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length, ...nextPositions);
    } catch (error) {
      if (
        refreshGeneration !== circularRecordRefreshGeneration ||
        recordDiscoverySuppressed() ||
        mode.value !== 'circular' ||
        cInputType.value !== inputType ||
        (inputType === 'gff' ? files.c_gff : files.c_gb) !== primaryFile ||
        (inputType === 'gff' ? files.c_fasta : null) !== pairedFile
      ) return;
      const busy = state.sessionOperationAvailability?.();
      if (busy) {
        circularRecordDiscovery.status = 'deferred';
        circularRecordDiscovery.primaryFile = null;
        return busy;
      }

      circularRecordList.value = [];
      circularRecordDiscovery.status = 'error';
      circularRecordDiscovery.error = formatError(error, inputType === 'gff' ? 'listGffFastaRecords' : 'listSequenceRecords', 'helper');
      drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length);
    }
  };

  const refreshCircularRecordOrder = (options = {}) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const inputType = cInputType.value;
    const fingerprint = [
      Boolean(options.suppress || recordDiscoverySuppressed()),
      mode.value,
      inputType,
      inputType === 'gff' ? files.c_gff : files.c_gb,
      inputType === 'gff' ? files.c_fasta : null
    ];
    const inflightRefresh = activeCircularRecordRefresh;
    if (
      inflightRefresh &&
      fingerprint.length === inflightRefresh.fingerprint.length &&
      fingerprint.every((value, index) => (
        Object.is(value, inflightRefresh.fingerprint[index])
      ))
    ) {
      return inflightRefresh.promise;
    }
    const entry = { fingerprint, promise: /** @type {Promise<any> | null} */ (null) };
    entry.promise = runCircularRecordRefresh(options).finally(() => {
      if (activeCircularRecordRefresh === entry) activeCircularRecordRefresh = null;
    });
    activeCircularRecordRefresh = entry;
    return entry.promise;
  };

  // R13: Generate and the label reflow prepare the specific-color rules they
  // render here, the Generate compiler's one trigger of the rule preparation.
  // `isCurrent` is the run's own staleness: the candidate is null when the run
  // went stale during the preparation, and is admitted (`shouldAdmit`) while
  // the run and the rule inputs it was prepared from stay current. The run
  // commits the feature color overrides rebound to the normalized captions.
  /**
   * @param {DrawingState} drawing
   * @param {() => boolean} isCurrent
   * @param {Record<string, any>} [options]
   */
  const prepareAndAdmitCandidate = async (drawing, isCurrent, options) => {
    if (!rulePreparation) {
      return { rules: drawing.manualSpecificRules, featureColorOverrides: drawing.featureColorOverrides, shouldAdmit: isCurrent, notifyChanges: () => {} };
    }
    const candidate = await rulePreparation.prepareCandidate(drawing.manualSpecificRules, options);
    if (!candidate || !isCurrent()) return null;
    return {
      rules: candidate.rules,
      featureColorOverrides: rebindRuleColorOverrides(drawing.featureColorOverrides, drawing.manualSpecificRules, candidate.rules),
      shouldAdmit: () => rulePreparation.isCurrent(candidate.snapshot) && isCurrent(),
      notifyChanges: () => rulePreparation.notifyChanges(candidate)
    };
  };

  /**
   * @param {DrawingState} drawing
   * @param {{ decorationContinuity?: any, comparisonPlanSnapshot?: Record<string, any> | null,
   *   generatedArtifactHandle?: Record<string, any> | null, comparisonExecution?: Record<string, any> | null,
   *   isCurrentOperation?: () => boolean, isCurrentAlert?: () => boolean }} [options]
   */
  const runAnalysisInternal = async (drawing, {
    decorationContinuity = null,
    comparisonPlanSnapshot = null,
    generatedArtifactHandle = null,
    comparisonExecution = null,
    isCurrentOperation = () => true,
    isCurrentAlert = () => true
  } = {}) => {
    /** @type {Awaited<ReturnType<typeof prepareAndAdmitCandidate>>} */
    let colorCandidate = null;
    let candidateRules = drawing.manualSpecificRules;
    // The request reads the run's drawing with the color rules the run
    // admits (`candidateRules`), and the label projection the run built.
    let requestDrawing = drawing;
    /** @type {ReturnType<typeof requestLabelProjection> | null} */
    let generatedLabelProjection = null;
    let failureStage = 'request-validation';
    recordSessionLifecycleEvent('generate-start');
    recordSessionLifecycleEvent('generation-input-resolution-start');
    const useCommittedComparison = comparisonExecution?.mode === 'inherit';
    const forceEmptyComparison = useCommittedComparison || comparisonExecution?.mode === 'clear';
    const activeComparisonPlanSnapshot = mode.value === 'linear'
      ? (
          forceEmptyComparison
            ? resolveLinearComparisonPlan({
                plan: { mode: 'none', defaultSource: 'losat', edges: [] },
                sequences: linearSeqs,
                layout: drawing.linearRecordLayoutEnabled.value ? drawing.linearRecordRows : [],
                losatProgram: drawing.losatProgram.value,
                blastpMode: normalizeBlastpMode(drawing.losat.blastp?.mode)
              })
            : comparisonPlanSnapshot || resolveLinearComparisonPlan({
                plan: drawing.linearComparisonPlan,
                sequences: linearSeqs,
                layout: drawing.linearRecordLayoutEnabled.value ? drawing.linearRecordRows : [],
                losatProgram: drawing.losatProgram.value,
                blastpMode: normalizeBlastpMode(drawing.losat.blastp?.mode)
              })
        )
      : null;
    let linearRecordCatalog = null;

    const generationToken = ++latestGenerationToken;
    let canceledAttemptOwnsPresentation = false;
    /** @type {AbortController | null} */
    let generationAbortController = null;
    /** @type {AbortSignal | null} */
    let generationAbortSignal = null;
    const committedArtifactHandle = generatedArtifactHandle || await captureGeneratedArtifactHandle();
    /** @type {ReturnType<typeof generatedArtifactTransactionOwner.build> | null} */
    let activatedGeneratedArtifactCandidate = null;
    /** @type {ReadyReceipt | null} */
    let acceptedCandidateReadyReceipt = null;
    let workingProteinIdentityManifest = committedArtifactHandle?.ownerSet
      ?.proteinIdentityManifest ?? proteinIdentityManifest.value;
    let workingLegacyProteinRawCandidates = committedArtifactHandle?.ownerSet
      ?.legacyProteinRawCandidates ?? legacyProteinRawCandidates.value;
    let workingLegacyProteinDerivedEvidence = committedArtifactHandle?.ownerSet
      ?.legacyProteinDerivedEvidence ?? legacyProteinDerivedEvidence.value;
    let workingOrthogroups = committedArtifactHandle?.ownerSet?.orthogroups
      ?? orthogroups.value;
    let workingExtractedFeatures = committedArtifactHandle?.ownerSet?.extractedFeatures
      ?? extractedFeatures.value;
    let workingBiologicalFeatures = committedArtifactHandle?.ownerSet?.biologicalFeatures
      ?? biologicalFeatures?.value;
    const hasSourceBoundEditorIntent = Object.keys(drawing.featureOverrides).length > 0
      || Object.keys(drawing.featurePlacementOverrides || {}).length > 0
      || Object.keys(drawing.featureStrokeOverrides).length > 0
      || Object.keys(drawing.legendColorOverrides).length > 0 || Object.keys(drawing.legendStrokeOverrides).length > 0
      || drawing.legendEntries.value.some(entry => entry.originalCaption && entry.originalCaption !== entry.caption)
      || drawing.dormantLegendEntries.value.length > 0;
    // The committed records whose source this Generate replaces or drops. The
    // metadata of a Session without a feature catalog has no record keys.
    const replacedFeatures = hasSourceBoundEditorIntent
      ? [...new Map((workingBiologicalFeatures || []).filter(feature => feature.record_key)
        .map(feature => [feature.record_key, feature])).values()]
        .filter(feature => !isCurrentFeature(feature))
      : [];
    const sourceReplaced = replacedFeatures.length > 0;
    const previousCommittedRequest = typeof getCommittedCanonicalSession === 'function'
      ? getCommittedCanonicalSession()?.renderRequest : null;
    const previousRequestRecords = previousCommittedRequest?.mode === mode.value
      ? previousCommittedRequest.records || [] : [];
    let workingLosatCacheInfo = committedArtifactHandle?.ownerSet?.losatCacheInfo
      ?? losatCacheInfo.value;
    let workingSelectedOrthogroupId = selectedOrthogroupId.value;
    let workingSelectedOrthogroupAlignmentFeature = selectedOrthogroupAlignmentFeature.value;
    const restoreCommittedArtifact = async () => {
      if (!committedArtifactHandle || !activatedGeneratedArtifactCandidate) return false;
      if (captureGeneratedArtifactOwnerSet().results !== activatedGeneratedArtifactCandidate.ownerSet.results) return false;
      if (acceptedCandidateReadyReceipt) {
        previewRuntime.invalidateReadyReceipt(
          acceptedCandidateReadyReceipt,
          'The activated candidate entered rollback.'
        );
      } else {
        previewRuntime.invalidateReadinessExpectation(
          String(generationToken),
          new Error('The activated candidate entered rollback.')
        );
      }
      recordSessionLifecycleEvent('artifact.rollback-started');
      await previewRuntime.restorePreviousSelectedResult({
        handle: committedArtifactHandle,
        restore: () => generatedArtifactTransactionOwner.restore(committedArtifactHandle)
      });
      activatedGeneratedArtifactCandidate = null;
      acceptedCandidateReadyReceipt = null;
      recordSessionLifecycleEvent('artifact.rollback-completed');
      return true;
    };
    const finishCanceledManualRun = async () => {
      if (!isCurrentOperation() || latestGenerationToken !== generationToken + 1) {
        return { status: 'stale' };
      }
      canceledAttemptOwnsPresentation = true;
      await restoreCommittedArtifact();
      if (!isCurrentOperation()) return { status: 'stale' };
      if (isCurrentAlert()) errorLog.value = null;
      processingStatus.value = 'Canceled.';
      generationCancelRequested.value = false;
      return { status: 'canceled' };
    };
    const setProcessingStatus = (message) => {
      if (generationToken === latestGenerationToken) {
        processingStatus.value = String(message || '');
      }
    };
    const onDiagramProgress = ({ stage }) => {
      const message = {
        'preparing-runtime': 'Preparing diagram runtime (first use)...',
        'preparing-resources': 'Preparing diagram input resources...',
        rendering: 'Rendering diagram...',
        finalizing: 'Finalizing diagram results...'
      }[stage];
      if (message) setProcessingStatus(message);
    };
    const throwIfGenerationCanceled = () => {
      if (generationCancelRequested.value) {
        throw new DiagramGenerationCanceledError();
      }
    };
    generationAbortController = typeof AbortController === 'function' ? new AbortController() : null;
    generationAbortSignal = generationAbortController?.signal || null;
    activeLosatAbortController = generationAbortController;
    if (mode.value === 'linear' && typeof prepareLinearRecordCatalog === 'function') {
      processingStatus.value = 'Reading input records...';
      let prepared;
      try {
        prepared = await prepareLinearRecordCatalog();
      } catch (error) {
        prepared = {
          catalog: null,
          error
        };
      } finally {
        if (generationToken === latestGenerationToken) {
          processingStatus.value = 'Preparing input files...';
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
      if (prepared?.error) {
        const error = formatError(prepared.error);
        if (activeLosatAbortController === generationAbortController) {
          activeLosatAbortController = null;
        }
        return failOperation(error, { handle: committedArtifactHandle,
          restore: restoreCommittedArtifact, isCurrent: isCurrentAlert, isCurrentOperation });
      }
      linearRecordCatalog = prepared?.catalog || null;
    }
    if (mode.value === 'circular') {
      const inputType = cInputType.value;
      const primaryFile = inputType === 'gff' ? files.c_gff : files.c_gb;
      const pairedFile = inputType === 'gff' ? files.c_fasta : null;
      const hasCompleteInput = Boolean(primaryFile && (inputType !== 'gff' || pairedFile));
      if (
        hasCompleteInput &&
        !circularDiscoveryMatchesCurrentInput()
      ) {
        processingStatus.value = 'Reading input records...';
        try {
          await refreshCircularRecordOrder({ reuseFinalError: true });
        } finally {
          if (generationToken === latestGenerationToken) {
            processingStatus.value = 'Preparing input files...';
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
          const message = circularRecordDiscovery.error || diagnosticError('INPUT_UNREADABLE');
          const outcome = await failOperation(message, { handle: committedArtifactHandle,
            restore: restoreCommittedArtifact, isCurrent: isCurrentAlert, isCurrentOperation });
          if (activeLosatAbortController === generationAbortController) {
            activeLosatAbortController = null;
          }
          return outcome;
        }
      }
    }
    const depthInputError = validateDepthInputPresence(drawing);
    if (depthInputError) {
      const error = formatError(depthInputError);
      if (activeLosatAbortController === generationAbortController) {
        activeLosatAbortController = null;
      }
      return failOperation(error, { handle: committedArtifactHandle,
        restore: restoreCommittedArtifact, isCurrent: isCurrentAlert, isCurrentOperation });
    }
    const previousSelectedResultIndex = selectedResultIndex.value;
    const activeRunColors = drawing.currentColors.value;
    const manualRunStartedAt = getNow();
    const manualRunStartedAtIso = new Date().toISOString();
    /** @type {Record<string, any> | null} */
    let structuredLosatTelemetry = null;
    const legacyPromotionTransaction = [];
    let legacyPromotionCommitted = false;
    /** @type {((candidateOwnerSet: Record<string, any>) => { ownerSet: Record<string, any>, selectedOrthogroupAlignmentFeature: any }) | null} */
    let commitProteinMigration = null;
    /** @type {{ cacheInfo: any[], cacheMap: Map<any, any>, derivedCacheMap: Map<any, any> | null } | null} */
    let pendingLosatCacheCommit = null;

    featureExtractionRequestId += 1;
    processingStatus.value = 'Preparing input files...';
    resultPanelTab.value = 'preview';
    if (isCurrentAlert()) errorLog.value = null;
    skipCaptureBaseConfig.value = false;
    closeLabelTextScopeDialog();
    window._origPairwiseMin = activeRunColors.pairwise_match_min || '#FFE7E7';
    window._origPairwiseMax = activeRunColors.pairwise_match_max || '#FF7272';
    clearLabelBuildNotices();

    try {
      colorCandidate = await prepareAndAdmitCandidate(
        drawing,
        () => generationToken === latestGenerationToken && !generationCancelRequested.value,
        { onProgress: onDiagramProgress }
      );
      if (!colorCandidate) {
        throwIfGenerationCanceled();
        return { status: 'stale' };
      }
      candidateRules = colorCandidate.rules;
      requestDrawing = { ...drawing, manualSpecificRules: candidateRules };
      if (mode.value === 'linear') {
        if (!activeComparisonPlanSnapshot || !Array.isArray(activeComparisonPlanSnapshot.edges)) {
          throw new Error('A resolved Linear comparison plan is required.');
        }
        if (activeComparisonPlanSnapshot.error) {
          throw new Error(activeComparisonPlanSnapshot.error);
        }
      }
      if (typeof validateAnnotationTargets === 'function' && drawing.annotationSets.length > 0) {
        const annotationError = validateAnnotationTargets();
        if (annotationError) throw diagnosticError(annotationError.code, annotationError.context);
      }
      let regionSpecs = [];
      let recordSelectors = [];
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
      /**
       * @param {File | null | undefined} fileObj
       * @param {string} path
       * @param {{ cacheText?: boolean, textCache?: WeakMap<WeakKey, any> | null, slot?: string }} [options]
       */
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

      const normalizedOutputPrefix = String(drawing.form.prefix || '').trim();

      const activePaletteName = String(
        drawing.selectedPalette?.value || appliedPaletteName.value || 'default'
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

      const tContent = serializeSpecificRules(candidateRules);
      if (tContent.trim() !== '') {
        stageTextFile('/combined_t.tsv', tContent);
      }

      if (Array.isArray(drawing.annotationSets) && drawing.annotationSets.length > 0) {
        stageTextFile('/web_annotations.tsv', encodeAnnotationTable(drawing.annotationSets), {
          name: 'annotations.tsv',
          slot: 'generatedFiles.web_annotations'
        });
      }

      if (drawing.filterMode.value === 'Whitelist') {
        const wlContent = serializeLabelWhitelistRules(drawing.manualWhitelist);
        if (wlContent) {
          stageTextFile('/manual_wl.tsv', wlContent);
        }
      }

      const pContent = serializeQualifierPriorityRules(drawing.manualPriorityRules);
      if (pContent.trim() !== '') {
        stageTextFile('/priority.tsv', pContent);
      }

      // The staged label table is the request's: label rules only; per-feature
      // label edits travel as featureOverrides rows (design Q4). One projection
      // per Generate serves the staged copy and the request (CW-02).
      generatedLabelProjection = requestLabelProjection(state, requestDrawing);
      const labelTableTsv = requestLabelTableTsv(state, requestDrawing, generatedLabelProjection);
      if (labelTableTsv) {
        stageTextFile('/web_label_table.tsv', labelTableTsv);
      }
      // Per-feature Feature visibility decides LOSATP's proteins as the rules do;
      // the helper resolves one record's rows as Generate does (R4).
      const proteinVisibilityRows = (recordKeys) => requestFeatureOverrides(
        drawing.featureOverrides, recordKeys.map((recordKey) => ({ recordKey }))
      ).filter((row) => row.featureVisibility !== null)
        .map((row) => ({ ...row, labelVisibility: null, labelText: null }));
      /** @type {string | null} */
      let featureVisibilityTablePath = null;
      let featureVisibilityCacheKey = '';
      const featureVisibilityTsv = serializeFeatureVisibilityRules(drawing.featureVisibilityRules?.value || []);
      if (featureVisibilityTsv.trim()) {
        featureVisibilityTablePath = '/web_feature_visibility_table.tsv';
        featureVisibilityCacheKey = featureVisibilityTsv;
        stageTextFile(featureVisibilityTablePath, featureVisibilityCacheKey);
      }
      const validateDepthStyleSettings = () => {
        if (
          drawing.adv.depth_min !== null &&
          drawing.adv.depth_min !== undefined &&
          drawing.adv.depth_min !== '' &&
          drawing.adv.depth_max !== null &&
          drawing.adv.depth_max !== undefined &&
          drawing.adv.depth_max !== '' &&
          Number(drawing.adv.depth_min) > Number(drawing.adv.depth_max)
        ) {
          throw new Error('Depth minimum must be less than or equal to depth maximum.');
        }
        if (
          drawing.adv.depth_large_tick_interval !== null &&
          drawing.adv.depth_large_tick_interval !== undefined &&
          drawing.adv.depth_large_tick_interval !== ''
        ) {
          if (Number(drawing.adv.depth_large_tick_interval) <= 0) throw new Error('Depth large tick interval must be greater than 0.');
        }
        if (
          mode.value === 'circular' &&
          drawing.adv.depth_small_tick_interval !== null &&
          drawing.adv.depth_small_tick_interval !== undefined &&
          drawing.adv.depth_small_tick_interval !== ''
        ) {
          if (Number(drawing.adv.depth_small_tick_interval) <= 0) throw new Error('Depth small tick interval must be greater than 0.');
        }
        if (
          drawing.adv.depth_tick_font_size !== null &&
          drawing.adv.depth_tick_font_size !== undefined &&
          drawing.adv.depth_tick_font_size !== ''
        ) {
          if (Number(drawing.adv.depth_tick_font_size) <= 0) throw new Error('Depth tick font size must be greater than 0.');
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
        if (!Array.isArray(drawing.adv.depth_tracks)) drawing.adv.depth_tracks = [];
        while (drawing.adv.depth_tracks.length <= idx) {
          const nextIndex = drawing.adv.depth_tracks.length;
          drawing.adv.depth_tracks.push({
            label: getDepthTrackFallbackLabel(nextIndex),
            color: nextIndex === 0 ? String(drawing.adv.depth_color || '#4A90E2') : '',
            height: null,
            large_tick_interval: null,
            small_tick_interval: null,
            tick_font_size: null
          });
        }
        const config = drawing.adv.depth_tracks[idx];
        if (!config || typeof config !== 'object' || Array.isArray(config)) {
          drawing.adv.depth_tracks[idx] = {
            label: getDepthTrackFallbackLabel(idx),
            color: idx === 0 ? String(drawing.adv.depth_color || '#4A90E2') : '',
            height: null,
            large_tick_interval: null,
            small_tick_interval: null,
            tick_font_size: null
          };
        }
        return drawing.adv.depth_tracks[idx];
      };
      const syncDepthSlotLegendLabelsFromTrackConfigs = (slots, trackFiles = []) => {
        if (!Array.isArray(slots)) return;
        const trackCount = Array.isArray(trackFiles)
          ? trackFiles.length
          : 0;
        for (let trackIndex = 0; trackIndex < trackCount; trackIndex += 1) {
          ensureDepthTrackConfigAt(trackIndex);
        }
        syncDepthSlotLabels(/** @type {any} */ ({
          slots,
          depthTracks: drawing.adv.depth_tracks,
          activeCount: trackCount
        }));
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
        drawing.adv.gc_content_mode = String(drawing.adv.gc_content_mode || '').trim().toLowerCase() === 'percent'
          ? 'percent'
          : 'deviation';
        if (drawing.adv.gc_content_mode !== 'percent') return;

        const minPercent = Number(drawing.adv.gc_content_min_percent);
        const maxPercent = Number(drawing.adv.gc_content_max_percent);
        if (!Number.isFinite(minPercent)) {
          throw new Error('GC content minimum percent must be a finite number.');
        }
        if (!Number.isFinite(maxPercent)) {
          throw new Error('GC content maximum percent must be a finite number.');
        }
        if (minPercent > maxPercent) {
          throw new Error('GC content minimum percent must be less than or equal to maximum percent.');
        }
        drawing.adv.gc_content_min_percent = minPercent;
        drawing.adv.gc_content_max_percent = maxPercent;

        if (
          drawing.adv.gc_content_tick_interval !== null &&
          drawing.adv.gc_content_tick_interval !== undefined &&
          drawing.adv.gc_content_tick_interval !== ''
        ) {
          if (Number(drawing.adv.gc_content_tick_interval) <= 0) throw new Error('GC content large tick interval must be greater than 0.');
        }
        if (
          drawing.adv.gc_content_small_tick_interval !== null &&
          drawing.adv.gc_content_small_tick_interval !== undefined &&
          drawing.adv.gc_content_small_tick_interval !== ''
        ) {
          if (Number(drawing.adv.gc_content_small_tick_interval) <= 0) throw new Error('GC content small tick interval must be greater than 0.');
        }
        if (
          drawing.adv.gc_content_tick_font_size !== null &&
          drawing.adv.gc_content_tick_font_size !== undefined &&
          drawing.adv.gc_content_tick_font_size !== ''
        ) {
          if (Number(drawing.adv.gc_content_tick_font_size) <= 0) throw new Error('GC content tick font size must be greater than 0.');
        }
      };
      normalizeGcContentPercentState();

      if (mode.value === 'circular') {
        const normalizedCircularPlotTitle = String(drawing.form.plot_title || '').trim();
        const normalizedPlotTitlePosition = normalizeCircularPlotTitlePosition(drawing.adv.plot_title_position);
        const useCircularTrackSlots = drawing.adv.circular_track_slots_enabled === true;
        const circularTrackAxisIndex = clampCircularTrackAxisIndex(
          drawing.adv.circular_track_slots_axis_index,
          Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots.length : 0
        );
        const circularTrackSlots = useCircularTrackSlots
          ? applyCircularSuppressControlsToSlots(
              applyCircularTrackOrderPlacements(
                drawing.adv.circular_track_slots,
                drawing.adv.nt,
                drawing.form.track_type,
                circularTrackAxisIndex
              ),
              drawing.form
            )
          : [];
        if (useCircularTrackSlots) {
          drawing.adv.circular_track_slots.splice(0, drawing.adv.circular_track_slots.length, ...circularTrackSlots);
          const circularTrackOnAxisIndex = circularTrackSlots.findIndex((slot) => slot?.side === 'overlay');
          const normalizedCircularTrackAxisIndex = circularTrackOnAxisIndex >= 0
            ? circularTrackOnAxisIndex
            : (
                circularTrackAxisIndex === null
                  ? inferLegacyAxisIndexFromFeature(circularTrackSlots, drawing.form.track_type)
                  : circularTrackAxisIndex
              );
          drawing.adv.circular_track_slots_axis_index = clampCircularTrackAxisIndex(
            normalizedCircularTrackAxisIndex,
            circularTrackSlots.length
          );
        }
        // Numeric settings (plot title font size, center radius, label
        // spacing, multi-record ratios, ring geometry) are projected literally
        // by the request and validated by Python; Generate keeps the draft.
        const keepFullDefinitionWithPlotTitle = Boolean(drawing.adv.keep_full_definition_with_plot_title);
        drawing.form.plot_title = normalizedCircularPlotTitle;
        drawing.adv.plot_title_position = normalizedPlotTitlePosition;
        drawing.adv.keep_full_definition_with_plot_title = keepFullDefinitionWithPlotTitle;

        const normalizedCircularLabelPlacement =
          String(drawing.adv.circular_label_placement || '').trim().toLowerCase() === 'radial'
            ? 'radial'
            : 'horizontal';
        drawing.adv.circular_label_placement = normalizedCircularLabelPlacement;
        if (drawing.form.multi_record_canvas) {
          const effectiveRecordPositions = mergeCircularRecordPositions(
            circularRecordList.value,
            drawing.adv.multi_record_positions
          );
          requireCurrentCircularMultiRecordSizeMode(drawing.adv.multi_record_size_mode);
          drawing.adv.multi_record_positions.splice(
            0,
            drawing.adv.multi_record_positions.length,
            ...effectiveRecordPositions
          );
        }

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
        }
        const hasCircularDepthFile = circularDepthEntries.length > 0;
        const circularSlotNeedsDepth = useCircularTrackSlots && hasEnabledCircularTrackRenderer(circularTrackSlots, 'depth');
        if (drawing.form.show_depth || circularSlotNeedsDepth) {
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

        assertActiveModeInputs?.('circular', state);

        const sourceMode = String(drawing.circularConservation.source || '').trim().toLowerCase() === 'upload'
          ? 'upload'
          : 'losat';
        const circularConservationSourceFiles = sourceMode === 'upload'
          ? normalizeFileList(files.c_conservation_blasts)
          : normalizeFileList(files.c_conservation_fastas);
        const shouldDrawCircularPairwiseComparisons = circularConservationSourceFiles.length > 0;
        drawing.circularConservation.enabled = shouldDrawCircularPairwiseComparisons;

        if (shouldDrawCircularPairwiseComparisons) {
          setProcessingStatus('Preparing conservation comparisons...');
          // Thresholds are evaluated before conservation work; the request
          // projects the same resolution, and the draft keeps what was typed.
          resolveComparisonThresholds(drawing.adv, 'circular');

          const runCircularLosatConservation = async (comparisonEntries) => {
            const circularLosatProgram = normalizeCircularConservationLosatProgram(
              drawing.circularConservation.losat_program
            );
            drawing.circularConservation.losat_program = circularLosatProgram;
            const subjectGencode = normalizePositiveInteger(drawing.circularConservation.subject_gencode, 1);
            const buildExtraArgs = (comparisonGencode) => {
              if (circularLosatProgram === 'tblastx') {
                return [
                  '--query-gencode',
                  String(normalizePositiveInteger(comparisonGencode, 1)),
                  '--db-gencode',
                  String(subjectGencode)
                ];
              }
              const normalizedTask = String(drawing.losat.blastn?.task || 'megablast').trim() || 'megablast';
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
              // One sequence reader (D12): the Python helper reads FASTA, GenBank
              // or DDBJ and returns the query FASTA that the CLI ring hashes.
              const sequenceResponse = await runDiagramHelperOperation(
                DIAGRAM_HELPER_OPERATIONS.READ_COMPARISON_SEQUENCE,
                { files: [{ role: 'source', bytes: await cloneFileBytesForTransfer(fileObj) }] }
              );
              if (sequenceResponse.result?.error) throw sequenceResponse.result.error;
              const queryFasta = String(sequenceResponse.result?.fasta || '');
              const queryHash = await hashText(queryFasta);
              const querySequenceKey = `circular-query:${queryHash}`;
              sequenceEntriesByKey.set(querySequenceKey, queryFasta);
              const extraArgs = buildExtraArgs(comparisonEntry?.losat_gencode);
              const cacheMetadata = {
                flow: 'circular-conservation',
                program: circularLosatProgram,
                outfmt: String(drawing.losat.outfmt || '6'),
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
              const cached = getReusableLosatCacheEntry(cacheMap, cacheKey, cacheMetadata);
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
                  outfmt: drawing.losat.outfmt || '6',
                  extraArgs
                });
              }
            }

            if (losatJobs.length > 0) {
              setProcessingStatus('Preparing comparison search runtime...');
              const runtime = await prepareLosatRuntime({ includeThreaded: executionMode !== 'serial' }).catch((error) => {
                console.warn('LOSAT runtime warmup failed; execution will report the error.', formatError(error));
                return null;
              });
              if (runtime?.threaded && losatThreadingStatus) {
                const { wasmModule: _wasmModule, ...threadedStatus } = runtime.threaded;
                losatThreadingStatus.value = threadedStatus;
              }
              setProcessingStatus(`Running ${circularLosatSuffix.toUpperCase()} conservation: 0/${losatJobs.length} jobs complete`);
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
                  outfmt: String(drawing.losat.outfmt || '6'),
                  args: job?.extraArgs || [],
                  queryCanonicalHash: job?.queryCanonicalHash || '',
                  subjectCanonicalHash: job?.subjectCanonicalHash || '',
                  runtime: webLosatRuntimeRecord(circularLosatProgram)
                });
              });
            } else {
              setProcessingStatus('Using cached LOSAT conservation results...');
            }
            throwIfGenerationCanceled();
            retainCompletedLosatSearch(cacheMap, losatPairs);

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

          if (sourceMode === 'upload') {
            const blastFiles = circularConservationSourceFiles;
            if (blastFiles.length === 0) {
              throw new Error('Please upload at least one BLAST outfmt 6/7 file for Pairwise Comparisons.');
            }
            const rawReference = String(drawing.circularConservation.reference || 'auto').trim().toLowerCase();
            drawing.circularConservation.reference = ['query', 'subject'].includes(rawReference)
              ? rawReference
              : 'auto';
            workingLosatCacheInfo = [];
          } else {
            const comparisonFiles = circularConservationSourceFiles;
            if (comparisonFiles.length === 0) {
              throw new Error('Please upload at least one comparison sequence file for Pairwise Comparisons.');
            }
            const conservationEntries = orderedConservationSources(comparisonFiles, drawing.circularConservation);
            const conservationSeries = buildConservationSeries(comparisonFiles, drawing.circularConservation);
            const conservationResults = await runCircularLosatConservation(conservationEntries);
            drawing.circularConservation.reference = 'subject';
            resolvedCircularConservation = conservationResults.map((result, index) => {
              /** @type {Record<string, any>} */
              const source = conservationEntries[index] || {};
              /** @type {Record<string, any>} */
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
          workingLosatCacheInfo = [];
        }
      } else {
        requireCurrentLinearTrackLayout(drawing.form.linear_track_layout);
        const useLinearTrackSlots = drawing.adv.linear_track_slots_enabled === true;
        let linearTrackSlots = [];
        /** @type {number | null} */
        let linearTrackSlotAxisIndex = null;
        let linearSlotNeedsDepth = false;
        if (useLinearTrackSlots) {
          linearTrackSlots = normalizeLinearTrackSlots(
            drawing.adv.linear_track_slots,
            drawing.adv.nt,
            drawing.form.linear_track_layout
          );
          linearSlotNeedsDepth = linearTrackSlots.some(
            (slot) => slot.enabled !== false && slot.renderer === 'depth'
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
          linearTrackSlotAxisIndex = clampLinearTrackAxisIndex(
            drawing.adv.linear_track_slots_axis_index,
            linearTrackSlots.length
          );
          linearTrackSlotAxisIndex = resolveLinearTrackAxisIndex(
            linearTrackSlots,
            linearTrackSlotAxisIndex
          );
          linearTrackSlots = applyLinearTrackOrderPlacements(
            linearTrackSlots,
            linearTrackSlotAxisIndex,
            drawing.adv.nt,
            drawing.form.linear_track_layout
          );
          drawing.adv.linear_track_slots_axis_index = linearTrackSlotAxisIndex;
          drawing.adv.linear_track_slots.splice(0, drawing.adv.linear_track_slots.length, ...linearTrackSlots);
        }
        // This branch is Linear mode, where the 'A resolved Linear comparison plan is required.' check above
        // throws when the snapshot is null.
        const comparisonResolution = /** @type {NonNullable<typeof activeComparisonPlanSnapshot>} */ (activeComparisonPlanSnapshot);
        const hasComparisonIntent = comparisonResolution.hasComparisonIntent === true;
        const hasLosatIntent = comparisonResolution.hasLosatIntent === true;
        const useProteinBlastp = hasLosatIntent && drawing.losatProgram.value === 'blastp';
        const blastpMode = useProteinBlastp
          ? requireCurrentProteinBlastpMode(drawing.losat.blastp?.mode)
          : String(drawing.losat.blastp?.mode ?? '');
        const usePairwiseBlastp = useProteinBlastp && blastpMode === 'pairwise';
        const useOrthogroupBlastp = useProteinBlastp && blastpMode === 'orthogroup';
        const useCollinearBlastp = useProteinBlastp && blastpMode === 'collinear';
        if (hasComparisonIntent) {
          setProcessingStatus('Preparing comparisons...');
          drawing.adv.pairwise_match_style = normalizeCurrentPairwiseMatchStyle(
            drawing.adv.pairwise_match_style,
            'ribbon'
          );
        }
        // One resolution per Generate feeds LOSAT post-processing and caches.
        const comparisonThresholds = hasComparisonIntent
          ? resolveComparisonThresholds(drawing.adv, 'linear')
          : null;

        const blastpMaxHits = usePairwiseBlastp
          ? requireCurrentProteinBlastpMaxHits(drawing.losat.blastp?.maxHits)
          : 5;
        const blastpCandidateLimit = useProteinBlastp
          ? resolveProteinBlastpCandidateLimit(drawing.losat.blastp?.candidateLimit)
          : null;
        const orthogroupMembershipMode = useOrthogroupBlastp || useCollinearBlastp
          ? requireCurrentOrthogroupMembershipMode(
              drawing.losat.blastp?.orthogroupMembershipMode
            )
          : 'anchor_core_v1';
        const orthogroupMemberMaxHits = useOrthogroupBlastp || useCollinearBlastp
          ? requireCurrentOrthogroupMemberMaxHits(
              drawing.losat.blastp?.orthogroupMemberMaxHits
            )
          : null;
        const collinearMinAnchors = useCollinearBlastp
          ? requireCurrentCollinearMinAnchors(drawing.losat.blastp?.collinearMinAnchors)
          : 1;
        const collinearMaxUnitGap = useCollinearBlastp
          ? requireCurrentCollinearMaxUnitGap(drawing.losat.blastp?.collinearMaxUnitGap)
          : 0;
        const collinearMaxDiagonalDrift = useCollinearBlastp
          ? requireCurrentCollinearMaxDiagonalDrift(
              drawing.losat.blastp?.collinearMaxDiagonalDrift
            )
          : 0;
        const collinearMaxConflictsInMergeGap = useCollinearBlastp
          ? requireCurrentCollinearMaxConflicts(
              drawing.losat.blastp?.collinearMaxConflictsInMergeGap
            )
          : 1;
        const collinearMaxParalogLinksPerOrthogroup = useCollinearBlastp
          ? requireCurrentCollinearMaxParalogLinks(
              drawing.losat.blastp?.collinearMaxParalogLinksPerOrthogroup
            )
          : 2;
        const collinearColorMode = useCollinearBlastp
          ? requireCurrentCollinearColorMode(drawing.losat.blastp?.collinearColorMode)
          : 'orientation';
        const collinearUnitMode = useCollinearBlastp
          ? requireCurrentCollinearUnitMode(drawing.losat.blastp?.collinearUnitMode)
          : 'auto';
        const collinearAnchorMode = useCollinearBlastp
          ? requireCurrentCollinearAnchorMode(drawing.losat.blastp?.collinearAnchorMode)
          : 'rbh';
        const collinearMergeOrientation = useCollinearBlastp
          ? requireCurrentCollinearMergeOrientation(
              drawing.losat.blastp?.collinearMergeOrientation
            )
          : 'either';
        const collinearInferOrthogroups = requireCurrentCollinearInferOrthogroups(drawing.losat.blastp?.collinearInferOrthogroups);
        const collinearSearchScope = useCollinearBlastp
          ? requireCurrentCollinearSearchScope(drawing.losat.blastp?.collinearSearchScope)
          : 'adjacent';

        const reuseResolvedProteinArtifacts = useProteinBlastp
          && !workingLegacyProteinRawCandidates?.entries?.length
          && canReuseResolvedProteinArtifacts({
            canonicalComparisons: files.linearCanonicalComparisons,
            sequences: linearSeqs,
            committedSession: typeof getCommittedCanonicalSession === 'function'
              ? getCommittedCanonicalSession()
              : null,
            active: {
              featureVisibility: featureVisibilityCacheKey,
              featureOverrides: proteinVisibilityRows(linearSeqs.map((sequence) => String(sequence?.uid || ''))),
              mode: blastpMode,
              candidateLimit: blastpCandidateLimit,
              ...comparisonThresholds,
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
              inferOrthogroups: collinearInferOrthogroups,
              searchScope: collinearSearchScope
            }
          });
        const useLosat = hasLosatIntent && !reuseResolvedProteinArtifacts;
        if (useLosat) recordSessionLifecycleEvent('losat-cache-preparation-start');

        if (useProteinBlastp) {
          drawing.losat.blastp.mode = blastpMode;
          drawing.losat.blastp.candidateLimit = blastpCandidateLimit;
        }
        if (usePairwiseBlastp) {
          drawing.losat.blastp.maxHits = blastpMaxHits;
        } else if (useOrthogroupBlastp) {
          drawing.losat.blastp.orthogroupMembershipMode = orthogroupMembershipMode;
          drawing.losat.blastp.orthogroupMemberMaxHits = orthogroupMemberMaxHits;
        } else if (useCollinearBlastp) {
          drawing.losat.blastp.orthogroupMembershipMode = orthogroupMembershipMode;
          drawing.losat.blastp.orthogroupMemberMaxHits = orthogroupMemberMaxHits;
          drawing.losat.blastp.collinearMinAnchors = collinearMinAnchors;
          drawing.losat.blastp.collinearMaxUnitGap = collinearMaxUnitGap;
          drawing.losat.blastp.collinearMaxDiagonalDrift = collinearMaxDiagonalDrift;
          drawing.losat.blastp.collinearMaxConflictsInMergeGap = collinearMaxConflictsInMergeGap;
          drawing.losat.blastp.collinearMaxParalogLinksPerOrthogroup =
            collinearMaxParalogLinksPerOrthogroup;
          drawing.losat.blastp.collinearColorMode = collinearColorMode;
          drawing.losat.blastp.collinearUnitMode = collinearUnitMode;
          drawing.losat.blastp.collinearAnchorMode = collinearAnchorMode;
          drawing.losat.blastp.collinearMergeOrientation = collinearMergeOrientation;
          drawing.losat.blastp.collinearSearchScope = collinearSearchScope;
        }

        const normalizedPlotTitle = String(drawing.form.plot_title || '').trim();
        const normalizedPlotTitlePosition = normalizeLinearPlotTitlePosition(drawing.adv.plot_title_position);
        drawing.adv.linear_show_replicon = drawing.adv.linear_show_replicon === true;
        drawing.adv.linear_accession_visibility = requireLinearLabelVisibilityMode(
          drawing.adv.linear_accession_visibility,
          'Linear Accession visibility'
        );
        drawing.adv.linear_length_visibility = requireLinearLabelVisibilityMode(
          drawing.adv.linear_length_visibility,
          'Linear Length / Coordinates visibility'
        );
        drawing.adv.linear_definition_line_styles = normalizeDefinitionLineStyleState(drawing.adv.linear_definition_line_styles);
        drawing.form.plot_title = normalizedPlotTitle;
        drawing.adv.plot_title_position = normalizedPlotTitlePosition;

        requireCurrentLinearLabelPlacement(drawing.adv.label_placement);
        if (hasComparisonIntent) {
          const comparisonHeight = classifyOptionalPositiveNumber(drawing.adv.comparison_height);
          if (comparisonHeight.status === 'invalid') {
            throw diagnosticError('INPUT_INVALID', { field: 'match_height', reason: 'POSITIVE_OR_AUTO' });
          }
        }

        const viewTransformSpecs = [];
        const buildRegionSpec = (seq, idx) => {
          const recordIdRaw = seq.region_record_id ? String(seq.region_record_id).trim() : '';
          const wantsReverse = Boolean(seq.region_reverse);
          const bounds = resolveLinearRegionBounds(seq, idx);

          recordSelectors.push(recordIdRaw || '');

          if (bounds) {
            const { start, end } = bounds;
            const canonicalStart = Math.min(start, end);
            const canonicalEnd = Math.max(start, end);
            const coordinateReverse = start > end;
            const displayReverse = wantsReverse || coordinateReverse;
            const specBody = `${start}-${end}${wantsReverse ? ':rc' : ''}`;
            const canonicalSpecBody = `${canonicalStart}-${canonicalEnd}`;
            viewTransformSpecs.push({ reverse: displayReverse });
            return { file: canonicalSpecBody, displayFile: specBody };
          }

          viewTransformSpecs.push({ reverse: wantsReverse });
          return null;
        };

        const fastaCache = new Map();
        const fastaHashCache = new Map();
        const sequenceEntriesByKey = new Map();
        const linearFileTextCache = new WeakMap();
        if (!reuseResolvedProteinArtifacts) {
          if (useOrthogroupBlastp) {
            workingOrthogroups = [];
          } else {
            workingOrthogroups = [];
            workingSelectedOrthogroupId = '';
            workingSelectedOrthogroupAlignmentFeature = '';
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
              rawJobs: /** @type {Record<string, any>[]} */ ([]),
              cacheHashHits: 0,
              cacheHits: 0,
              cacheMisses: 0,
              totalPairs: 0,
              uniqueJobs: 0
            }
          : null;
        const losatExecutionMode = useLosat ? getLosatExecutionMode() : 'serial';
        const losatRequestedThreadsPerJob = useLosat ? getLosatThreadsPerJob(drawing) : undefined;
        const losatRequestedTotalThreadBudget = useLosat ? getLosatTotalThreadBudget() : undefined;
        const losatRuntimeWarmup = useLosat
          ? prepareLosatRuntime({ includeThreaded: losatExecutionMode !== 'serial' }).then((runtime) => {
              if (runtime?.threaded && losatThreadingStatus) {
                const { wasmModule: _wasmModule, ...threadedStatus } = runtime.threaded;
                losatThreadingStatus.value = threadedStatus;
              }
              return runtime;
            }).catch((error) => {
              console.warn('LOSAT runtime warmup failed; execution will report the error if LOSAT is used.', formatError(error));
              return null;
            })
          : null;

        if (!hasLosatIntent) {
          workingLosatCacheInfo = [];
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
          throwIfGenerationCanceled();
          if (fastaCache.has(idx)) return fastaCache.get(idx);
          const startedAt = getNow();
          const fmt = lInputType.value === 'gb'
            ? 'genbank'
            : (useProteinBlastp ? 'gff' : 'fasta');
          const regionSpec = regionSpecs[idx]?.file || null;
          const recordSelector = recordSelectors[idx] ?? '';
          const reverseFlag = /** @type {string} */ ('0');
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
          const recordVisibilityRows = useProteinBlastp ? proteinVisibilityRows([recordInstanceKey]) : [];
          const persistentCacheKey = useProteinBlastp
            ? JSON.stringify({
                inputFormat: lInputType.value,
                pairedFastaFile: getFileFingerprint(pairedFastaFile),
                regionSpec,
                recordSelector,
                recordInstanceKey,
                recordIndex: idx,
                featureVisibility: featureVisibilityCacheKey,
                featureOverrides: recordVisibilityRows,
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
                  recordInstanceKey,
                  featureOverrides: recordVisibilityRows.length ? recordVisibilityRows : null
                }
              );
              const res = response.result;
              if (res.error) throw res.error;
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
                if (res.error) throw res.error;
                entry = {
                  fasta: res.fasta,
                  recordId: res.record_id || `seq_${idx + 1}`,
                  canonicalLength: Number(res.record_length || 0) || getFastaSequenceLength(res.fasta)
                };
                if (losatTiming) {
                  losatTiming.fastaWorkerFallbacks += 1;
                  console.warn('LOSAT browser FASTA extraction fell back to the diagram Worker.', formatError(fastError));
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
              outfmt: String(drawing.losat.outfmt || '6'),
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
            program: drawing.losatProgram.value,
            outfmt: String(drawing.losat.outfmt || '6'),
            args: argsKey,
            queryCanonicalHash: queryHash,
            subjectCanonicalHash: subjectHash
          };
        };

        const tryPromoteLegacyProteinEntry = async ({
          cacheKey,
          metadata,
          queryEntry,
          subjectEntry,
          identityIndex
        }) => {
          if (
            metadata?.identityKind !== 'protein' ||
            !workingLegacyProteinRawCandidates
          ) return null;
          const envelope = workingLegacyProteinRawCandidates;
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
              identityManifest: cloneJsonData(workingProteinIdentityManifest),
              expectedOptions: {
                program: metadata.program,
                outfmt: metadata.outfmt,
                args: normalizeLosatArgs(metadata.args)
              }
            }
          );
          const result = response.result;
          if (result.status === 'error') {
            throw result.error || new Error('Legacy protein cache migration failed.');
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
            workingLegacyProteinRawCandidates = transitionLegacyProteinCandidate(
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
            workingProteinIdentityManifest,
            { identityIndex }
          );
          if (!verified) {
            cacheMap.delete(cacheKey);
            workingLegacyProteinRawCandidates = transitionLegacyProteinCandidate(
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
          const fallback = losatPairDefaultName(
            drawing,
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

        const getBlastpCandidateLimit = () => {
          if (!useProteinBlastp) return null;
          return blastpCandidateLimit;
        };

        const buildLosatArgs = (queryIdx, subjectIdx) => {
          const args = [];
          if (drawing.losatProgram.value === 'blastn') {
            pushArg(args, '--task', drawing.losat.blastn.task);
          } else if (drawing.losatProgram.value === 'tblastx') {
            pushArg(args, '--query-gencode', losatRecordGencode(linearSeqs[queryIdx]));
            pushArg(args, '--db-gencode', losatRecordGencode(linearSeqs[subjectIdx]));
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
          assertActiveModeInputs?.('linear', state);
          for (let i = 0; i < linearSeqs.length; i++) {
            const seq = linearSeqs[i];
            if (lInputType.value === 'gb') {
              if (useLosat && losatRecordIndexes.has(i)) {
                await prepareLinearFile(seq.gb, `/seq_${i}.gb`, {
                  cacheText: true,
                  slot: `files.linearSeqs[${i}].gb`
                });
              }
            } else {
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
          workingProteinIdentityManifest = mergeProteinIdentityManifests(manifests, {
            invalidInputMessage:
              'Protein comparison metadata could not be validated. Reload the page and try again.'
          });
          const legacyReferenceIds = collectLegacyProteinReferences(
            workingOrthogroups,
            workingSelectedOrthogroupAlignmentFeature,
            workingExtractedFeatures,
            workingBiologicalFeatures
          );
          if (legacyReferenceIds.length > 0) {
            const response = await runDiagramHelperOperation(
              DIAGRAM_HELPER_OPERATIONS.RESOLVE_LEGACY_PROTEIN_REFERENCES,
              {
                proteinRecords: cloneJsonData(proteinEntries.map((entry) => ({
                  proteinMap: entry.proteinMap || {},
                  fasta: entry.fasta || ''
                }))),
                identityManifest: cloneJsonData(workingProteinIdentityManifest),
                referenceIds: cloneJsonData(legacyReferenceIds)
              }
            );
            const result = response.result;
            if (result.status !== 'resolved' || !result.proteinIdMap) {
              throw result.error || new Error('Legacy protein UI reference migration failed.');
            }
            const proteinIdMap = result.proteinIdMap;
            workingOrthogroups = rewriteMappedProteinReferences(
              workingOrthogroups,
              proteinIdMap
            );
            workingSelectedOrthogroupAlignmentFeature = rewriteMappedProteinReferences(
              workingSelectedOrthogroupAlignmentFeature,
              proteinIdMap
            );
            workingExtractedFeatures = rewriteMappedProteinReferences(
              workingExtractedFeatures,
              proteinIdMap
            );
            if (workingBiologicalFeatures) {
              workingBiologicalFeatures = rewriteMappedProteinReferences(
                workingBiologicalFeatures,
                proteinIdMap
              );
            }
            if (collectLegacyProteinReferences(
              workingOrthogroups,
              workingSelectedOrthogroupAlignmentFeature,
              workingExtractedFeatures,
              workingBiologicalFeatures
            ).length > 0) {
              throw new Error('Legacy protein UI references remain after migration.');
            }
          }
        }

        // `losatTiming` is built exactly when `useLosat` holds; the second test narrows it for this block.
        if (useLosat && losatTiming) {
          setProcessingStatus('Preparing LOSAT jobs...');
          const losatPairs = [];
          const losatJobs = [];
          const pendingJobKeys = new Set();
          const jobBuildStartedAt = getNow();
          const fastaExtractionBeforeJobBuild = losatTiming.fastaExtractionMs;
          const cacheHashBeforeJobBuild = losatTiming.cacheHashMs;

          const jobSpecs = buildLosatJobSpecs(/** @type {any} */ ({
            resolution: comparisonResolution,
            recordCount: linearSeqs.length,
            recordUids: linearSeqs.map((seq) => seq.uid),
            program: drawing.losatProgram.value,
            blastpMode,
            collinearInferOrthogroups,
            collinearSearchScope
          }));

          const sourcePlan = await prepareLosatSourceBatches({
            sequences: linearSeqs,
            specs: jobSpecs,
            getEntry: getSeqEntry,
            buildArgs: buildLosatArgs,
            hashText,
            protein: useProteinBlastp
          });
          const preparedJobs = [];
          for (const spec of jobSpecs) {
            throwIfGenerationCanceled();
            const losatArgs = buildLosatArgs(spec.queryIndex, spec.subjectIndex);
            const cacheMetadata = await buildCacheMetadata(
              losatArgs, spec.queryIndex, spec.subjectIndex
            );
            const batch = sourcePlan.bySpec.get(spec);
            if (batch.searchContext) cacheMetadata.searchContext = batch.searchContext;
            preparedJobs.push({ spec, losatArgs, cacheMetadata, batch });
          }
          /** @type {string[] | null} */
          let proteinCacheKeys = null;
          if (useProteinBlastp && preparedJobs.length > 0) {
            throwIfGenerationCanceled();
            const response = await runDiagramHelperOperation(
              DIAGRAM_HELPER_OPERATIONS.BUILD_PROTEIN_LOSAT_CACHE_KEYS,
              {
                identityManifest: cloneJsonData(workingProteinIdentityManifest),
                pairs: preparedJobs.map(({ cacheMetadata }) => ({
                  queryRecordInstanceKey: cacheMetadata.queryRecordInstanceKey,
                  subjectRecordInstanceKey: cacheMetadata.subjectRecordInstanceKey,
                  expectedOptions: {
                    program: cacheMetadata.program,
                    outfmt: cacheMetadata.outfmt,
                    args: normalizeLosatArgs(cacheMetadata.args),
                    ...(cacheMetadata.searchContext ? { searchContext: cacheMetadata.searchContext } : {})
                  }
                }))
              }
            );
            throwIfGenerationCanceled();
            const result = response.result;
            if (
              result.error || !Array.isArray(result.keys)
              || result.keys.length !== preparedJobs.length
              || result.keys.some((key) => !/^[0-9a-f]{64}$/.test(key))
            ) {
              throw result.error || new Error('Protein cache key generation failed.');
            }
            proteinCacheKeys = result.keys;
          }

          // mergeProteinIdentityManifests made a private deep copy for this run.
          // Nothing mutates or publishes it during this loop; helper calls receive
          // clones, including legacy promotion. Release before search/render and
          // publication so no index can outlive input, Session or History changes.
          const identityIndex = useProteinBlastp && preparedJobs.length > 0
            ? buildValidatedProteinIdentityIndex(workingProteinIdentityManifest)
            : null;
          if (useProteinBlastp && preparedJobs.length > 0 && !identityIndex) {
            throw new Error('Protein comparison identity manifest is invalid.');
          }
          // CLI Sessions can specify a complete two-dimensional row layout.
          // For collinearity, every pair across adjacent rows is displayed;
          // the Web comparison-plan edges alone do not represent those pairs.
          const canonicalGridRows = useCollinearBlastp && drawing.linearRecordLayoutEnabled.value
            ? linearSeqs.map((seq) => drawing.linearRecordRows.find((entry) => entry.uid === seq.uid))
            : null;
          const hasCanonicalGridRows = canonicalGridRows?.length > 1
            && canonicalGridRows.every((entry) => entry?.canonicalRow === entry.row
              && Number.isInteger(entry.canonicalColumn) && entry.canonicalColumn > 0
              && entry.canonicalCardinality === 'exactly_one');
          try {
            for (const [jobIndex, { spec, losatArgs, cacheMetadata, batch }] of preparedJobs.entries()) {
              throwIfGenerationCanceled();
              const queryEntry = await getSeqEntry(spec.queryIndex);
              throwIfGenerationCanceled();
              const subjectEntry = await getSeqEntry(spec.subjectIndex);
              throwIfGenerationCanceled();
              const cacheKey = useProteinBlastp
                // The loop runs only for a non-empty `preparedJobs`, and then `proteinCacheKeys` was set above for a protein run.
                ? /** @type {string[]} */ (proteinCacheKeys)[jobIndex]
                : await hashText(JSON.stringify(buildLosatCachePayload(cacheMetadata)));
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
              let cached = getReusableLosatCacheEntry(
                cacheMap,
                cacheKey,
                cacheMetadata,
                workingProteinIdentityManifest,
                identityIndex
              );
              if (!cached && useProteinBlastp && !cacheMetadata.searchContext) {
                cached = await tryPromoteLegacyProteinEntry({
                  cacheKey,
                  metadata: cacheMetadata,
                  queryEntry,
                  subjectEntry,
                  identityIndex
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
              const isResolvedDisplayPair = hasCanonicalGridRows
                ? spec.queryIndex < spec.subjectIndex && Math.abs(
                    canonicalGridRows[spec.queryIndex].row
                    - canonicalGridRows[spec.subjectIndex].row
                  ) === 1
                : Boolean(
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
                  program: drawing.losatProgram.value,
                  querySequenceKey: queryEntry.sequenceKey,
                  subjectSequenceKey: subjectEntry.sequenceKey,
                  queryCanonicalHash,
                  subjectCanonicalHash,
                  outfmt: drawing.losat.outfmt || '6',
                  extraArgs: losatArgs,
                  cacheMetadata,
                  batch
                });
              }
            }
          } finally {
            releaseValidatedProteinIdentityIndex(identityIndex);
          }
          const sourceJobs = [];
          for (const batch of sourcePlan.batches) {
            const members = losatJobs.filter((job) => job.batch === batch);
            if (members.length === 0) continue;
            sequenceEntriesByKey.set(batch.query.sequenceKey, batch.query.fasta);
            sequenceEntriesByKey.set(batch.subject.sequenceKey, batch.subject.fasta);
            sourceJobs.push({
              ...members[0],
              cacheKey: await hashText(JSON.stringify([
                drawing.losatProgram.value, drawing.losat.outfmt || '6', batch.args,
                batch.query.hash, batch.subject.hash
              ])),
              querySequenceKey: batch.query.sequenceKey,
              subjectSequenceKey: batch.subject.sequenceKey,
              queryRecordIndexes: batch.query.indexes,
              subjectRecordIndexes: batch.subject.indexes,
              recordPairs: members.map((job) => [job.queryIndex, job.subjectIndex]),
              members,
              batch
            });
          }
          losatTiming.rawJobs = sourceJobs.map((job) => ({
            queryRecordIndexes: job.queryRecordIndexes,
            subjectRecordIndexes: job.subjectRecordIndexes,
            recordPairs: job.recordPairs,
            cacheKey: job.cacheKey,
            args: [...job.extraArgs]
          }));
          losatTiming.uniqueJobs = sourceJobs.length;
          const jobBuildWallMs = getNow() - jobBuildStartedAt;
          const nestedFastaMs = losatTiming.fastaExtractionMs - fastaExtractionBeforeJobBuild;
          const nestedHashMs = losatTiming.cacheHashMs - cacheHashBeforeJobBuild;
          losatTiming.jobBuildWallMs += jobBuildWallMs;
          losatTiming.jobBuildMs += Math.max(0, jobBuildWallMs - nestedFastaMs - nestedHashMs);

          if (sourceJobs.length > 0) {
            // A failure while searching is reported as a LOSAT failure (CO-01).
            failureStage = 'losat';
            setProcessingStatus('Preparing comparison search runtime...');
            const runtimeWaitStartedAt = getNow();
            await waitForCancelablePromise(losatRuntimeWarmup, generationAbortSignal);
            throwIfGenerationCanceled();
            losatTiming.runtimeWaitMs += getNow() - runtimeWaitStartedAt;
            setProcessingStatus(`Running LOSAT: 0/${sourceJobs.length} source jobs complete`);
            const executionStartedAt = getNow();
            const sourceResults = await executeLosatJobs(sourceJobs.map(({ members, batch, ...job }) => job), {
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
                setProcessingStatus(`Running LOSAT: ${completed}/${total} source jobs complete`);
              }
            });
            throwIfGenerationCanceled();
            failureStage = 'request-validation';
            losatTiming.executionMs += getNow() - executionStartedAt;
            const losatResults = sourceResults.flatMap((result) => {
              const job = sourceJobs.find((item) => item.cacheKey === result.cacheKey);
              if (!job) throw new Error('LOSAT returned an unknown source job.');
              return splitLosatSourceResult(result.text, job.batch, job.members);
            });
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
                program: drawing.losatProgram.value,
                outfmt: String(drawing.losat.outfmt || '6'),
                args: job?.extraArgs || [],
                ...(cacheMetadata.searchContext ? { searchContext: cacheMetadata.searchContext } : {}),
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
                    }),
                runtime: webLosatRuntimeRecord(drawing.losatProgram.value)
              };
              cacheMap.set(result.cacheKey, rawEntry);
            });
          } else {
            setProcessingStatus('Using cached LOSAT results...');
          }
          throwIfGenerationCanceled();
          // Legacy promotions must repeat their transaction until a Result commits.
          retainCompletedLosatSearch(cacheMap, losatPairs.filter(({ cacheKey }) => (
            !legacyPromotionTransaction.some((promotion) => promotion.cacheKey === cacheKey)
          )));

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
            /** @type {Record<string, any> | null} */
            let convertedPayload = null;
            if (useDerivedProteinPayloadCache) {
              const derivedCachePayload = buildLosatDerivedPayloadCachePayload({
                mode: blastpMode,
                maxHits: blastpMaxHits,
                ...comparisonThresholds,
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
                collinearInferOrthogroups,
                orthogroupMembershipMode,
                orthogroupMemberMaxHits,
                explicitDisplayPairs: Boolean(hasCanonicalGridRows),
                recordPayloads,
                pairPayloads
              });
              derivedCacheKey = await hashText(JSON.stringify(derivedCachePayload));
              convertedPayload = getLosatDerivedCacheEntry(
                derivedCacheMap,
                derivedCacheKey,
                workingProteinIdentityManifest
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
                  ...comparisonThresholds,
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
                  collinearInferOrthogroups,
                  orthogroupMembershipMode,
                  orthogroupMemberMaxHits,
                  explicitDisplayPairs: Boolean(hasCanonicalGridRows)
                }
              );
              // The converter returns a JSON object (convert_losatp_blastp_pairs_to_genomic_payload in python-helpers.js).
              convertedPayload = /** @type {Record<string, any>} */ (response.result);
              if (["orthogroup", "collinear"].includes(blastpMode)) {
                const resourceKey = blastpMode === "collinear" ? "collinearityResult" : "orthogroupResult";
                const canonical = convertedPayload.canonicalResource;
                if (!(canonical?.bytes instanceof Uint8Array) || canonical.bytes.byteLength !== canonical.size) {
                  throw new Error("Analysis canonical resource bytes are missing.");
                }
                convertedPayload[resourceKey] = JSON.parse(bytesToText(canonical.bytes));
                bindCanonicalTypedResource(convertedPayload[resourceKey], {
                  kind: canonical.kind, type: "application/json", size: canonical.size,
                  lastModified: 0, encoding: "base64", data: bytesToBase64(canonical.bytes)
                });
                delete convertedPayload.canonicalResource;
              }
              if (useDerivedProteinPayloadCache && !convertedPayload?.error) {
                setLosatDerivedCacheEntry(derivedCacheMap, derivedCacheKey, {
                  mode: blastpMode,
                  payload: convertedPayload,
                  manifest: workingProteinIdentityManifest
                });
              }
            }
            if (convertedPayload.error) throw convertedPayload.error;
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
            setProcessingStatus('Preparing nucleotide comparison results...');
            for (const pair of losatPairs) {
              throwIfGenerationCanceled();
              const cached = cacheMap.get(pair.cacheKey);
              // Raw search-frame rows; the Python planner projects orientation (PD-OI-073).
              const blastText = isCurrentRawLosatCacheEntry(cached) ? String(cached.text || '') : '';
              const blastPath = `/blast_${pair.pairIndex}.txt`;
              const blastName = pair.filename || getPayloadName(blastPath);
              const blastSlot = `generatedFiles.losat_blasts[${pair.pairIndex}]`;
              registerRunInfoFile(blastPath, {
                name: blastName,
                slot: blastSlot,
                kind: 'generated'
              });
              recordGeneratedCliFile(blastPath, blastText, {
                name: blastName,
                slot: blastSlot
              });
              resolvedComparisons.push({
                kind: 'nucleotideBlast',
                edgeKey: pair.edgeKey,
                ordinal: pair.ordinal,
                queryRecordIndex: pair.queryIndex,
                subjectRecordIndex: pair.subjectIndex,
                text: blastText,
              });
            }
          }
          losatTiming.blastWriteMs += getNow() - blastWriteStartedAt;
          if (legacyPromotionTransaction.length > 0) {
            commitProteinMigration = (candidateOwnerSet) => {
              let nextLegacyEnvelope = candidateOwnerSet.legacyProteinRawCandidates;
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
              let nextLegacyEvidence = candidateOwnerSet.legacyProteinDerivedEvidence;
              if (!nextLegacyEnvelope.entries.some((candidate) => candidate.state === 'pending')) {
                nextLegacyEvidence = { schema: 1, entries: [] };
              }
              let nextOrthogroups = candidateOwnerSet.orthogroups;
              let nextExtractedFeatures = candidateOwnerSet.extractedFeatures;
              let nextBiologicalFeatures = candidateOwnerSet.biologicalFeatures;
              let nextSelectedAlignmentFeature = workingSelectedOrthogroupAlignmentFeature;
              if (Object.keys(promotedProteinIdMap).length > 0) {
                nextOrthogroups = rewriteMappedProteinReferences(
                  nextOrthogroups,
                  promotedProteinIdMap
                );
                nextSelectedAlignmentFeature = rewriteMappedProteinReferences(
                  nextSelectedAlignmentFeature,
                  promotedProteinIdMap
                );
                nextExtractedFeatures = rewriteMappedProteinReferences(
                  nextExtractedFeatures,
                  promotedProteinIdMap
                );
                if (nextBiologicalFeatures) {
                  nextBiologicalFeatures = rewriteMappedProteinReferences(
                    nextBiologicalFeatures,
                    promotedProteinIdMap
                  );
                }
              }
              legacyPromotionCommitted = true;
              return {
                ownerSet: {
                  ...candidateOwnerSet,
                  legacyProteinRawCandidates: nextLegacyEnvelope,
                  legacyProteinDerivedEvidence: nextLegacyEvidence,
                  orthogroups: nextOrthogroups,
                  extractedFeatures: nextExtractedFeatures,
                  biologicalFeatures: nextBiologicalFeatures
                },
                selectedOrthogroupAlignmentFeature: nextSelectedAlignmentFeature
              };
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
            program: useProteinBlastp ? 'blastp' : drawing.losatProgram.value,
            outfmt: String(drawing.losat.outfmt || '6'),
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
        if ((!useLinearTrackSlots && drawing.form.show_depth) || linearSlotNeedsDepth) {
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
        }
      }

      if (drawing.annotationSets.length > 0) {
        stageTextFile('/web_annotations.tsv', encodeAnnotationTable(drawing.annotationSets), {
          name: 'annotations.tsv',
          slot: 'generatedFiles.web_annotations'
        });
      }

      throwIfGenerationCanceled();
      setProcessingStatus('Preparing render inputs and session...');
      await nextTick();
      await waitForAfterPaint();
      throwIfGenerationCanceled();
      if (typeof serializeCanonicalFiles !== 'function') {
        throw new Error('Canonical input serialization is unavailable.');
      }
      recordSessionLifecycleEvent('generation-input-resolution-end');
      recordSessionLifecycleEvent('serialize-canonical-files-start');
      const serializedFiles = await serializeCanonicalFiles(
        activeComparisonPlanSnapshot,
        linearRecordCatalog,
        requestDrawing
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
        drawing: requestDrawing,
        generatedLabelProjection,
        filesData: candidateFiles,
        recordDisplayRows: recordDisplayRows?.value || [],
        comparisonPlanSnapshot: activeComparisonPlanSnapshot,
        resolvedComparisons,
        resolvedCircularConservation: canonicalCircularConservation
      });
      if (useCommittedComparison) {
        if (typeof getCommittedCanonicalSession !== 'function') {
          throw new Error('The preserved comparison owner is unavailable.');
        }
        // The candidate is a current-schema request; an older committed
        // request (for example a v42 CLI sidecar) is promoted before reuse.
        const committedSession = getCommittedCanonicalSession();
        inheritCommittedComparisonIntent({
          candidate: canonical,
          committed: committedSession?.renderRequest
            ? {
                ...committedSession,
                renderRequest: promoteCanonicalRenderRequestToCurrent(
                  committedSession.renderRequest,
                  { featureCatalog: featureCatalog?.value || null }
                )
              }
            : committedSession
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
      /** @type {{ results: any[] | null, featureCatalog: Record<string, any> | null | undefined }} */
      const publishedReplayArtifact = { results: null, featureCatalog: null };
      const buildCanonicalReplayText = createCanonicalReplayTextBuilder({
        version: canonicalSessionVersion,
        startedAtIso: manualRunStartedAtIso,
        renderRequest: canonical.renderRequest,
        resources: canonical.resources,
        publishedArtifact: publishedReplayArtifact
      });
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
      /** @type {Awaited<ReturnType<typeof buildSourceRecipe>> | null} */
      let sourceRecipe = null;
      if (manualRunStartedAt !== null) {
        const generatedFileNameHints = new Map();
        generatedCliFileMap.forEach((file) => {
          const slot = String(file?.slot || '').trim();
          const name = String(file?.name || '').trim();
          if (slot && name && !generatedFileNameHints.has(slot)) {
            generatedFileNameHints.set(slot, name);
          }
        });
        // Establish source counts before Worker transfer releases the file-content cache.
        sourceRecipe = await buildSourceRecipe(/** @type {any} */ ({
          ...canonical,
          generatedFileNameHints,
          readResourceRecordCount: (resourceId, kind) => (
            readCanonicalResourceRecordCount(canonical.resources, resourceId, kind)
          )
        }));
        throwIfGenerationCanceled();
      }
      const postGbdrawTimingEntries = [];
      failureStage = 'render';
      const canonicalExecution = await executeCanonicalRenderCandidate({
        canonical,
        decorationContinuity,
        mode: mode.value,
        kind: 'generate',
        shouldAdmit: colorCandidate.shouldAdmit,
        onProgress: onDiagramProgress,
        prepareCommit: prepareCandidateCommit,
        prepareCommitInput: {
          sourceReplaced,
          featureColorOverrides: colorCandidate.featureColorOverrides,
          featureStrokeOverrides: drawing.featureStrokeOverrides,
          featureOverrides: drawing.featureOverrides,
          legendEntries: drawing.legendEntries.value,
          deletedLegendEntries: drawing.deletedLegendEntries.value,
          dormantLegendEntries: drawing.dormantLegendEntries.value,
          originalLegendOrder: originalLegendOrder.value,
          addedLegendCaptions: drawing.addedLegendCaptions.value,
          unrequestedDepthCaptions: unrequestedDepthCaptions(drawing, canonical),
          legendColorOverrides: drawing.legendColorOverrides,
          legendStrokeOverrides: drawing.legendStrokeOverrides,
          manualSpecificRules: candidateRules
        },
        timingEntries: postGbdrawTimingEntries
      });
      console.info(`gbdraw ${mode.value} typed request render: ${formatDuration(canonicalExecution.elapsedMs)}.`);
      setProcessingStatus('Preparing preview...');
      await nextTick();
      await waitForAfterPaint();
      throwIfGenerationCanceled();
      if (canonicalExecution.status === 'superseded'
        || generationToken !== latestGenerationToken) {
        if (generationAbortSignal?.aborted) {
          return finishCanceledManualRun();
        }
        return { status: 'stale' };
      }
      if (canonicalExecution.status === 'engine-error') {
        logPostGbdrawTimings(postGbdrawTimingEntries);
        // R6: Python names the failed placement row; the Result still shows its feature.
        return await failOperation(nameFeaturePlacementFailure(canonicalExecution.engineError,
          canonical.renderRequest, extractedFeatures.value), { handle: committedArtifactHandle,
          restore: restoreCommittedArtifact, isCurrent: isCurrentAlert, isCurrentOperation });
      }
      const {
        generationResponse,
        generationMetadata,
        catalogAdmission: candidateCatalogAdmission,
        catalog: candidateCatalog,
        commit: candidateCommit
      } = canonicalExecution;
      publishedReplayArtifact.results = candidateCommit.results;
      publishedReplayArtifact.featureCatalog = candidateCatalog;

      if (generationToken !== latestGenerationToken) {
        if (generationAbortSignal?.aborted) {
          return finishCanceledManualRun();
        }
        return { status: 'stale' };
      }

      /** @type {ReturnType<typeof buildRunInfo> | null} */
      let candidateRunInfo = null;
      /** @type {ReturnType<typeof buildLatestCliHelperFiles> | null} */
      let candidateCliHelpers = null;
      if (manualRunStartedAt !== null) {
        candidateRunInfo = buildRunInfo(/** @type {any} */ ({
          mode: mode.value,
          sourceRecipe,
          exactReplayArgs: ['--session', canonicalReplayPath],
          fileMetadata: runInfoFileMap,
          elapsedMs: getNow() - manualRunStartedAt,
          resultCount: candidateCommit.results.length,
          startedAtIso: manualRunStartedAtIso,
          losatComparisons: mode.value === 'linear' && activeComparisonPlanSnapshot?.hasLosatIntent === true,
          losatRuntimes: summarizeLosatRuntimes(
            pendingLosatCacheCommit?.cacheInfo || workingLosatCacheInfo,
            pendingLosatCacheCommit?.cacheMap || losatCache.value
          ),
          featureIdentityNotices: generationMetadata.featureIdentityNotices
        }));
        // Assigned above under the same `manualRunStartedAt !== null` test, which always holds: it is a number.
        /** @type {NonNullable<typeof sourceRecipe>} */ (sourceRecipe).generatedFiles.forEach((file) => {
          recordGeneratedCliFile(
            file.path,
            file.data instanceof Uint8Array ? textDecoder.decode(file.data) : file.data,
            { name: file.name, slot: file.slot }
          );
        });
        if (structuredLosatTelemetry) {
          candidateRunInfo.losatTelemetry = cloneJsonData(structuredLosatTelemetry);
        }
        candidateCliHelpers = buildLatestCliHelperFiles(
          drawing,
          candidateRunInfo,
          generatedCliFileMap,
          normalizedOutputPrefix || 'out'
        );
      }
      logPostGbdrawTimings(postGbdrawTimingEntries);
      await validateSimilarityAlignmentResetReceipt(
        state.similarityAlignmentResetReceipt?.value, canonical
      );
      recordSessionLifecycleEvent('preview-result-commit-start');
      const candidateGroups = Array.isArray(candidateCommit.featureState.orthogroups)
        ? candidateCommit.featureState.orthogroups
        : [];
      const candidateOrthogroupIndex = candidateCommit.featureState.featureOrthogroupIndex;
      const candidateExtractedFeatures = candidateCommit.featureState.extractedFeatures;
      const candidateBiologicalFeatures = candidateCommit.featureState.biologicalFeatures;
      const currentOwnerSet = captureGeneratedArtifactOwnerSet();
      /** @type {Record<string, any>} */
      let candidateOwnerSet = {
        ...currentOwnerSet,
        results: candidateCommit.results,
        featureCatalog: candidateCatalog,
        extractedFeatures: candidateExtractedFeatures,
        biologicalFeatures: candidateBiologicalFeatures,
        featureRecordIds: candidateCommit.featureState.featureRecordIds,
        orthogroups: candidateGroups,
        featureOrthogroupIndex: candidateOrthogroupIndex,
        collinearGroups: Array.isArray(candidateCommit.featureState.collinearGroups)
          ? candidateCommit.featureState.collinearGroups
          : [],
        trackSlotResolvedGeometry: generationMetadata.trackSlotGeometry || null,
        annotationWarnings: canonicalExecution.annotationWarnings,
        featureIdentityNotices: canonicalExecution.featureIdentityNotices,
        comparisonWarnings: canonicalExecution.comparisonWarnings,
        specificRules: candidateRules,
        fileLegendCaptions: new Set(candidateRules.filter(rule => rule.fromFile && rule.cap).map(rule => rule.cap)),
        proteinIdentityManifest: workingProteinIdentityManifest,
        legacyProteinRawCandidates: workingLegacyProteinRawCandidates,
        legacyProteinDerivedEvidence: workingLegacyProteinDerivedEvidence,
        losatCache: pendingLosatCacheCommit?.cacheMap || losatCache.value,
        losatDerivedCache:
          pendingLosatCacheCommit?.derivedCacheMap || losatDerivedCache.value,
        losatCacheInfo: pendingLosatCacheCommit?.cacheInfo || workingLosatCacheInfo,
        matchSequenceOwner: matchSequenceRegistry?.buildTrustedOwner?.(
          candidateCommit.featureState.sequenceSources
        ) || currentOwnerSet.matchSequenceOwner,
        lastRunInfo: candidateRunInfo,
        pairwiseMatchFactors: { ...(pairwiseMatchFactors?.value || {}) },
        editableLabels: [],
        generatedLegendPosition: drawing.form.legend,
        generatedMode: mode.value,
        generatedMultiRecordCanvas:
          mode.value === 'circular' ? Boolean(drawing.form.multi_record_canvas) : false,
        generatedCircularPlotTitlePosition: mode.value === 'circular'
          ? normalizeCircularPlotTitlePosition(drawing.adv.plot_title_position)
          : currentOwnerSet.generatedCircularPlotTitlePosition,
        appliedPaletteName: String(
          drawing.selectedPalette?.value || appliedPaletteName.value || 'default'
        ),
        appliedPaletteColors: { ...drawing.currentColors.value },
        pendingPaletteName: '',
        pendingPaletteColors: {}
      };
      let migratedSelectedAlignmentFeature = workingSelectedOrthogroupAlignmentFeature;
      if (commitProteinMigration) {
        const migrated = commitProteinMigration(candidateOwnerSet);
        candidateOwnerSet = migrated.ownerSet;
        migratedSelectedAlignmentFeature = migrated.selectedOrthogroupAlignmentFeature;
      }
      const generatedArtifactCandidate = generatedArtifactTransactionOwner.build(
        candidateOwnerSet,
        {
          runtimeState: {
            files: candidateCliHelpers?.files,
            archiveName: candidateCliHelpers?.archiveName,
            losatTelemetry: cloneJsonData(structuredLosatTelemetry)
          }
        }
      );
      const nextSelectedResultIndex = Math.max(
        0,
        Math.min(previousSelectedResultIndex, Math.max(0, candidateCommit.results.length - 1))
      );
      const selectedCandidateResult = candidateCommit.results[nextSelectedResultIndex] || null;
      if (!selectedCandidateResult) {
        throw new Error('The generated artifact has no selected preview Result.');
      }
      const selectedMutationOperations = (
        candidateCommit.mutationPlan?.operationsByResult?.[nextSelectedResultIndex]
        || null
      );
      const requiredLabelFeatureIds = forcedLabelFeatureIds(selectedMutationOperations, {
        features: [...(candidateCatalogAdmission.renderedFeaturesByResult[nextSelectedResultIndex]?.values() || [])],
        diagramOptions: canonical.renderRequest.diagramOptions
      });
      const optionalLabelFeatureIds = new Set(labelOperationIds([
        ...(selectedMutationOperations?.labelText || []),
        ...(selectedMutationOperations?.labelVisibility || [])
      ]).filter((renderedId) => !requiredLabelFeatureIds.includes(renderedId)));
      const candidatePreviewReadiness = previewRuntime.registerReadinessExpectation({
        result: selectedCandidateResult,
        resultIndex: nextSelectedResultIndex,
        artifactIdentity: generationResponse.artifactIdentity,
        generationToken: String(generationToken),
        catalogState: candidateCatalogAdmission,
        phase: 'generate',
        bindingOptions: {
          isIncrementalEdit: false,
          requiredLabelFeatureIds,
          optionalLabelFeatureIds: Object.freeze([...optionalLabelFeatureIds])
        },
        isCurrent: () => (
          latestGenerationToken === generationToken
          && !generationCancelRequested.value
          && Number(selectedResultIndex.value) === nextSelectedResultIndex
        )
      });
      activatedGeneratedArtifactCandidate = generatedArtifactCandidate;
      const previousOrthogroups = orthogroups.value;
      generatedArtifactTransactionOwner.activate(generatedArtifactCandidate, {
        selectedResultIndex: nextSelectedResultIndex
      });
      if (resultGenerationKey) resultGenerationKey.value += 1;
      selectedFeatureRecordIdx.value = 0;
      selectedOrthogroupAlignmentFeature.value = migratedSelectedAlignmentFeature;
      const candidateGroupIds = candidateGroups
        .map((group) => String(group?.id || '').trim())
        .filter(Boolean);
      rekeyCommittedOrthogroupOverrides(drawing, previousOrthogroups, candidateGroups);
      if (
        !workingSelectedOrthogroupId
        || !candidateGroupIds.includes(String(workingSelectedOrthogroupId || '').trim())
      ) {
        workingSelectedOrthogroupId = candidateGroupIds[0] || '';
      }
      selectedOrthogroupId.value = workingSelectedOrthogroupId;
      featureExtractionPending.value = false;
      featureExtractionError.value = null;
      Object.keys(drawing.featureColorOverrides).forEach((key) => delete drawing.featureColorOverrides[key]);
      Object.assign(
        drawing.featureColorOverrides,
        cloneJsonValue(candidateCommit.featureColorOverrides, {})
      );
      // Feature strokes stay in the draft when this Result does not draw their
      // feature (the other mode's, or one this request leaves out); the next
      // Generate that draws the feature draws them again (R2, OIPC-C06, OV-84).
      setFeatureEditorStatus({
        status: candidateExtractedFeatures.length ? 'summary-ready' : 'idle',
        generationId: featureExtractionRequestId,
        error: null,
        summaryCount: candidateExtractedFeatures.length,
        detailsCacheSize: 0
      });
      resultPanelTab.value = 'preview';
      if (typeof resetPreviewViewport === 'function') {
        resetPreviewViewport({ resetZoom: true });
      } else {
        zoom.value = 1.0;
      }
      recordSessionLifecycleEvent('preview-result-commit-end');
      acceptedCandidateReadyReceipt = await candidatePreviewReadiness.promise;
      if (generationToken !== latestGenerationToken || generationCancelRequested.value) {
        if (generationAbortSignal?.aborted || generationCancelRequested.value) {
          return finishCanceledManualRun();
        }
        await restoreCommittedArtifact();
        return { status: 'stale' };
      }
      await waitForPostBindFrame();
      recordStructuralMetric('previewPostBindFrameCount', 1, { phase: 'generate' });
      recordSessionLifecycleEvent('preview.post-bind-frame-completed', {
        resultIndex: nextSelectedResultIndex,
        rootGeneration: acceptedCandidateReadyReceipt.rootGeneration
      });
      if (generationToken !== latestGenerationToken || generationCancelRequested.value) {
        if (generationAbortSignal?.aborted || generationCancelRequested.value) {
          return finishCanceledManualRun();
        }
        await restoreCommittedArtifact();
        return { status: 'stale' };
      }
      // Owner decision Q3 = A (design Q4 3.4, 6.3): a Generate that replaced a
      // source removes only the edits whose feature that source no longer has,
      // and the edits of records the request dropped; Python names them.
      featureEditRemovalCount.value = sourceReplaced ? pruneUnmatchedFeatureOverrides({
        featureOverrides: drawing.featureOverrides,
        featurePlacementOverrides: drawing.featurePlacementOverrides,
        featureStrokeOverrides: drawing.featureStrokeOverrides,
        notices: canonicalExecution.featureIdentityNotices,
        replacedRecordKeys: replacedFeatures.filter((feature) => feature.scope === canonical.renderRequest.mode)
          .map((feature) => String(feature.record_key)),
        previousRecords: previousRequestRecords,
        currentRecords: canonical.renderRequest.records || [],
        biologicalFeatures: candidateBiologicalFeatures
      }) : 0;
      // OV-120: a source replacement also retires the renames of the Legend
      // rows the new Result does not draw (as the OV-84 strokes).
      if (sourceReplaced) drawing.dormantLegendEntries.value = [];
      if (typeof setGeneratedArtifactIdentity === 'function') {
        setGeneratedArtifactIdentity(generationResponse.artifactIdentity, {
          results: candidateCommit.results
        });
      }
      if (typeof adoptCanonicalRenderArtifacts === 'function') {
        // The committed request keeps the computed typed analysis resource, so
        // target-only renders and Similarity alignment use it without LOSATP.
        adoptCanonicalRenderArtifacts(canonical, { adoptOwnedRequest: true });
      }
      if (!useCommittedComparison && drawing.importedComparisonIntent) {
        Object.assign(
          drawing.importedComparisonIntent,
          createImportedComparisonIntentState(),
          { disposition: IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE }
        );
      }
      colorCandidate.notifyChanges();
      return {
        status: 'ok',
        generatedArtifactCandidate: activatedGeneratedArtifactCandidate
      };
    } catch (e) {
      if (!legacyPromotionCommitted) {
        legacyPromotionTransaction.forEach(({ cacheMap, cacheKey }) => {
          cacheMap.delete(cacheKey);
        });
      }
      if (isDiagramGenerationCanceled(e)) {
        return finishCanceledManualRun();
      }
      if (generationToken !== latestGenerationToken) {
        await restoreCommittedArtifact();
        return { status: 'stale' };
      }
      return await failOperation(e, { handle: committedArtifactHandle, stage: failureStage,
        restore: restoreCommittedArtifact, isCurrent: isCurrentAlert, isCurrentOperation });
    } finally {
      if (activeLosatAbortController === generationAbortController) {
        activeLosatAbortController = null;
      }
      if (generationToken === latestGenerationToken || canceledAttemptOwnsPresentation) {
        generationCancelRequested.value = false;
      }
    }
  };

  /**
   * @param {Record<string, any> | null} [comparisonPlanSnapshot]
   * @param {Record<string, any> | null} [generatedArtifactHandle]
   * @param {Record<string, any> | null} [comparisonExecution]
   * @param {{ prepareGenerate?: (() => Promise<Record<string, any>>) | null,
   *   afterGenerate?: ((outcome: Record<string, any> | null) => any) | null }} [options]
   */
  const runAnalysis = async (
    comparisonPlanSnapshot = null,
    generatedArtifactHandle = null,
    comparisonExecution = null,
    { prepareGenerate = null, afterGenerate = null } = {}
  ) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    /** @type {Record<string, any> | null} */
    let outcome = null;
    const operationId = ++latestOperationId;
    const isCurrentOperation = () => operationId === latestOperationId;
    const previousAlert = errorLog.value;
    const isCurrentAlert = () => isCurrentOperation()
      && (errorLog.value === previousAlert || errorLog.value === null);
    /** @type {Record<string, any> | null} */
    let beforeHandle = null;
    const initialResults = results.value;
    /** @type {string | null} */
    let historyRecovery = null;
    processing.value = true;
    processingStatus.value = 'Preparing input files...';
    generationCancelRequested.value = false;
    recordSessionLifecycleEvent('generate.processing-published');
    try {
      const decorationContinuity = captureDecorationContinuity(getCommittedCanonicalSession?.(), projectCompositionRecordIdentity, drawing);
      await nextTick();
      await waitForAfterPaint();
      if (!isCurrentOperation()) return { status: 'stale' };
      recordSessionLifecycleEvent('generate.paint-opportunity-completed');
      // Cancel is honored before rendering, including during preparation.
      const cancelBeforeRender = () => {
        if (!generationCancelRequested.value) return null;
        processingStatus.value = 'Canceled.';
        outcome = { status: 'canceled' };
        failedGeneratePreservedResult.value = results.value.length > 0;
        return outcome;
      };
      if (cancelBeforeRender()) return outcome;
      await settleComparisonRecordLabels();
      if (!isCurrentOperation()) return { status: 'stale' };
      if (cancelBeforeRender()) return outcome;
      if (prepareGenerate) {
        const prepared = await prepareGenerate();
        // A preparation that finishes after supersession or Cancel is discarded.
        if (!isCurrentOperation()) return { status: 'stale' };
        if (cancelBeforeRender()) return outcome;
        if (prepared.status !== 'ready') return prepared;
        comparisonPlanSnapshot = prepared.comparisonPlanSnapshot;
        comparisonExecution = prepared.comparisonExecution;
      }
      const execute = (handle) => {
        beforeHandle = handle || generatedArtifactHandle;
        return runAnalysisInternal(drawing, {
          decorationContinuity, comparisonPlanSnapshot, generatedArtifactHandle: beforeHandle,
          comparisonExecution, isCurrentOperation, isCurrentAlert
        });
      };
      outcome = typeof runGeneratedArtifactReplacement === 'function'
        ? await runGeneratedArtifactReplacement(
            'Generate diagram',
            execute,
            {
              shouldCommit: (result) => result?.status === 'ok',
              // A Generate that replaced a source removes the feature placements
              // it no longer resolves (Q3 = A). Placements are draft intent, not
              // part of the Result, so Undo and Redo restore them, and the
              // removal count, with this step.
              captureIntentCheckpoint: () => ({
                placements: cloneJsonData(drawing.featurePlacementOverrides) || {},
                removed: Number(featureEditRemovalCount.value) || 0
              }),
              restoreIntentCheckpoint: ({ placements, removed }) => {
                restorePlacements(drawing.featurePlacementOverrides, placements);
                featureEditRemovalCount.value = removed;
              },
              onCheckpointCapture: onGeneratedArtifactCheckpointCapture,
              restoreAppliedArtifact: async (beforeHandle) => {
                historyRecovery = 'restore-failed';
                const activeReceipt = previewRuntime.getActiveRuntime?.()?.readyReceipt || null;
                if (activeReceipt) {
                  previewRuntime.invalidateReadyReceipt(
                    activeReceipt,
                    'History finalization failed.'
                  );
                }
                recordSessionLifecycleEvent('artifact.rollback-started', {
                  phase: 'history-finalization'
                });
                await previewRuntime.restorePreviousSelectedResult({
                  handle: beforeHandle,
                  phase: 'history-finalization-rollback',
                  restore: () => generatedArtifactTransactionOwner.restore(beforeHandle)
                });
                recordSessionLifecycleEvent('artifact.rollback-completed', {
                  phase: 'history-finalization'
                });
                historyRecovery = beforeHandle.ownerSet.results.length ? 'restored' : 'no-result';
              }
            }
          )
        : await execute(generatedArtifactHandle || await captureGeneratedArtifactHandle());
      if (outcome?.status === 'ok' && outcome.generatedArtifactCandidate) {
        generatedArtifactTransactionOwner.finalize();
        completedLosatSearch = null;
        // The committed Results draw every edit, so the last live edit failure
        // no longer applies (OV-36). A failed Generate keeps the Result it
        // could not replace, and the note with it.
        clearLabelBuildNotices({ rerender: true });
        recordSessionLifecycleEvent('generate.completed');
      }
      if (outcome && Object.prototype.hasOwnProperty.call(outcome, 'generatedArtifactCandidate')) {
        const { generatedArtifactCandidate, ...publicOutcome } = outcome;
        outcome = publicOutcome;
      }
      if (outcome?.status === 'canceled' && isCurrentOperation()) {
        failedGeneratePreservedResult.value = results.value.length > 0;
      } else if (outcome?.status === 'ok') {
        failedGeneratePreservedResult.value = false;
        if (generationFailureRecovery) generationFailureRecovery.value = null;
      }
      await afterGenerate?.(outcome);
      return outcome;
    } catch (cause) {
      outcome = await failOperation(cause, { handle: beforeHandle, isCurrent: isCurrentAlert, isCurrentOperation,
        recovery: cause?.artifactRestoreFailed ? 'restore-failed'
          : historyRecovery || (!beforeHandle && results.value === initialResults && initialResults.length ? 'preserved' : null) });
      return outcome;
    } finally {
      if (isCurrentOperation()) {
        if (outcome?.status !== 'canceled') processingStatus.value = '';
        generationCancelRequested.value = false;
        processing.value = false;
        recordSessionLifecycleEvent('generate.processing-cleared', {
          status: outcome?.status || 'error'
        });
        replayDeferredLabelReflow();
      }
    }
  };

  /**
   * @param {DrawingState} drawing
   * @param {{ canonical: Record<string, any>, decorationContinuity?: any, generatedArtifactHandle?: Record<string, any> | null,
   *   commitIntent?: (() => any) | null, alignmentResetBefore?: any, alignmentResetReceipt?: any, operation?: string,
   *   isCurrentOperation?: () => boolean, isCurrentAlert?: () => boolean }} options
   */
  const runCommittedCanonicalCandidateInternal = async (drawing, {
    canonical,
    decorationContinuity = null,
    generatedArtifactHandle = null,
    commitIntent = null,
    alignmentResetBefore = null,
    alignmentResetReceipt = undefined,
    operation = 'generate',
    isCurrentOperation = () => true,
    isCurrentAlert = () => true
  }) => {
    const generationToken = ++latestGenerationToken;
    const generationAbortController = typeof AbortController === 'function'
      ? new AbortController()
      : null;
    const generationAbortSignal = generationAbortController?.signal || null;
    activeLosatAbortController = generationAbortController;
    const committedArtifactHandle = generatedArtifactHandle
      || await captureGeneratedArtifactHandle();
    /** @type {ReturnType<typeof generatedArtifactTransactionOwner.build> | null} */
    let activatedCandidate = null;
    /** @type {ReadyReceipt | null} */
    let acceptedReadyReceipt = null;
    const restoreCommittedArtifact = async () => {
      if (!activatedCandidate) return false;
      if (captureGeneratedArtifactOwnerSet().results !== activatedCandidate.ownerSet.results) return false;
      if (acceptedReadyReceipt) {
        previewRuntime.invalidateReadyReceipt(
          acceptedReadyReceipt,
          'The target-only candidate entered rollback.'
        );
      } else {
        previewRuntime.invalidateReadinessExpectation(
          String(generationToken),
          new Error('The target-only candidate entered rollback.')
        );
      }
      recordSessionLifecycleEvent('artifact.rollback-started');
      await previewRuntime.restorePreviousSelectedResult({
        handle: committedArtifactHandle,
        restore: () => generatedArtifactTransactionOwner.restore(committedArtifactHandle)
      });
      activatedCandidate = null;
      acceptedReadyReceipt = null;
      recordSessionLifecycleEvent('artifact.rollback-completed');
      return true;
    };
    const finishCanceled = async () => {
      if (!isCurrentOperation()) return { status: 'stale' };
      await restoreCommittedArtifact();
      if (!isCurrentOperation()) return { status: 'stale' };
      if (isCurrentAlert()) errorLog.value = null;
      processingStatus.value = 'Canceled.';
      return { status: 'canceled' };
    };
    try {
      const timingEntries = [];
      const execution = await executeCanonicalRenderCandidate({
        canonical,
        decorationContinuity,
        mode: canonical.renderRequest.mode,
        kind: 'target-record-transform',
        shouldAdmit: () => generationToken === latestGenerationToken
          && !generationCancelRequested.value,
        onProgress: ({ stage }) => {
          const message = {
            'preparing-runtime': 'Preparing diagram runtime (first use)...',
            'preparing-resources': 'Preparing diagram input resources...',
            rendering: 'Rendering diagram...',
            finalizing: 'Finalizing diagram results...'
          }[stage];
          if (message && generationToken === latestGenerationToken) {
            processingStatus.value = message;
          }
        },
        prepareCommit: prepareCandidateCommit,
        prepareCommitInput: {
          sourceReplaced: false,
          featureColorOverrides: drawing.featureColorOverrides,
          featureStrokeOverrides: drawing.featureStrokeOverrides,
          featureOverrides: drawing.featureOverrides,
          legendEntries: drawing.legendEntries.value,
          deletedLegendEntries: drawing.deletedLegendEntries.value,
          dormantLegendEntries: drawing.dormantLegendEntries.value,
          originalLegendOrder: originalLegendOrder.value,
          addedLegendCaptions: drawing.addedLegendCaptions.value,
          unrequestedDepthCaptions: unrequestedDepthCaptions(drawing, canonical),
          legendColorOverrides: drawing.legendColorOverrides,
          legendStrokeOverrides: drawing.legendStrokeOverrides,
          manualSpecificRules: drawing.manualSpecificRules
        },
        timingEntries
      });
      if (execution.status === 'superseded'
        || generationToken !== latestGenerationToken) {
        if (generationAbortSignal?.aborted || generationCancelRequested.value) {
          return finishCanceled();
        }
        return { status: 'stale' };
      }
      if (execution.status === 'engine-error') {
        return await failOperation(execution.engineError, { handle: committedArtifactHandle, operation,
          restore: restoreCommittedArtifact, isCurrent: isCurrentAlert, isCurrentOperation });
      }
      processingStatus.value = 'Preparing preview...';
      const candidateCommit = execution.commit;
      const candidateCatalogAdmission = execution.catalogAdmission;
      const candidateCatalog = execution.catalog;
      const candidateGroups = Array.isArray(candidateCommit.featureState.orthogroups)
        ? candidateCommit.featureState.orthogroups
        : [];
      const receipt = alignmentResetBefore
        ? await buildSimilarityAlignmentResetReceipt({before:alignmentResetBefore, after:canonical})
        : alignmentResetReceipt === undefined ? state.similarityAlignmentResetReceipt?.value
          : alignmentResetReceipt;
      await validateSimilarityAlignmentResetReceipt(receipt, canonical);
      const currentOwnerSet = captureGeneratedArtifactOwnerSet();
      const candidateOwnerSet = {
        ...currentOwnerSet,
        similarityAlignmentPlan: canonical.renderRequest.layout?.similarityAlignment ?? null,
        linearRecordTranslations: canonical.renderRequest.layout?.recordTranslations || [],
        similarityAlignmentResetReceipt: receipt ?? null,
        results: candidateCommit.results,
        featureCatalog: candidateCatalog,
        extractedFeatures: candidateCommit.featureState.extractedFeatures,
        biologicalFeatures: candidateCommit.featureState.biologicalFeatures,
        featureRecordIds: candidateCommit.featureState.featureRecordIds,
        orthogroups: candidateGroups,
        featureOrthogroupIndex: candidateCommit.featureState.featureOrthogroupIndex,
        collinearGroups: Array.isArray(candidateCommit.featureState.collinearGroups)
          ? candidateCommit.featureState.collinearGroups
          : [],
        trackSlotResolvedGeometry:
          execution.generationMetadata.trackSlotGeometry || null,
        annotationWarnings: execution.annotationWarnings,
        featureIdentityNotices: execution.featureIdentityNotices,
        comparisonWarnings: execution.comparisonWarnings,
        matchSequenceOwner: matchSequenceRegistry?.buildTrustedOwner?.(
          candidateCommit.featureState.sequenceSources
        ) || currentOwnerSet.matchSequenceOwner,
        editableLabels: []
      };
      const preservedRuntimeState = captureGeneratedArtifactRuntimeState();
      const generatedCandidate = generatedArtifactTransactionOwner.build(
        candidateOwnerSet,
        { runtimeState: {
          files: preservedRuntimeState.latestCliHelperFiles,
          archiveName: preservedRuntimeState.latestCliHelperArchiveName,
          losatTelemetry: preservedRuntimeState.losatTelemetry
        } }
      );
      const nextSelectedResultIndex = Math.max(
        0,
        Math.min(
          selectedResultIndex.value,
          Math.max(0, candidateCommit.results.length - 1)
        )
      );
      const selectedCandidateResult = candidateCommit.results[nextSelectedResultIndex];
      if (!selectedCandidateResult) {
        throw new Error('The generated artifact has no selected preview Result.');
      }
      const candidatePreviewReadiness = previewRuntime.registerReadinessExpectation({
        result: selectedCandidateResult,
        resultIndex: nextSelectedResultIndex,
        artifactIdentity: execution.generationResponse.artifactIdentity,
        generationToken: String(generationToken),
        catalogState: candidateCatalogAdmission,
        phase: 'target-record-transform',
        bindingOptions: { isIncrementalEdit: false },
        isCurrent: () => generationToken === latestGenerationToken
          && !generationCancelRequested.value
          && Number(selectedResultIndex.value) === nextSelectedResultIndex
      });
      activatedCandidate = generatedCandidate;
      generatedArtifactTransactionOwner.activate(generatedCandidate, {
        selectedResultIndex: nextSelectedResultIndex
      });
      if (resultGenerationKey) resultGenerationKey.value += 1;
      featureExtractionPending.value = false;
      featureExtractionError.value = null;
      setFeatureEditorStatus({
        status: candidateOwnerSet.extractedFeatures.length ? 'summary-ready' : 'idle',
        generationId: featureExtractionRequestId,
        error: null,
        summaryCount: candidateOwnerSet.extractedFeatures.length,
        detailsCacheSize: 0
      });
      resultPanelTab.value = 'preview';
      acceptedReadyReceipt = await candidatePreviewReadiness.promise;
      if (generationToken !== latestGenerationToken || generationCancelRequested.value) {
        if (generationAbortSignal?.aborted || generationCancelRequested.value) {
          return finishCanceled();
        }
        await restoreCommittedArtifact();
        return { status: 'stale' };
      }
      await waitForPostBindFrame();
      if (generationToken !== latestGenerationToken || generationCancelRequested.value) {
        if (generationAbortSignal?.aborted || generationCancelRequested.value) {
          return finishCanceled();
        }
        await restoreCommittedArtifact();
        return { status: 'stale' };
      }
      if (typeof setGeneratedArtifactIdentity === 'function') {
        setGeneratedArtifactIdentity(execution.generationResponse.artifactIdentity, {
          results: candidateCommit.results
        });
      }
      if (typeof adoptCanonicalRenderArtifacts === 'function') {
        adoptCanonicalRenderArtifacts(canonical, { adoptOwnedRequest: true });
      }
      if (typeof commitIntent === 'function') await commitIntent();
      if (isCurrentAlert()) errorLog.value = null;
      failedGeneratePreservedResult.value = false;
      if (generationFailureRecovery) generationFailureRecovery.value = null;
      logPostGbdrawTimings(timingEntries);
      return { status: 'ok', generatedArtifactCandidate: activatedCandidate };
    } catch (error) {
      if (isDiagramGenerationCanceled(error)) return finishCanceled();
      if (generationToken !== latestGenerationToken) {
        await restoreCommittedArtifact();
        return { status: 'stale' };
      }
      return await failOperation(error, { handle: committedArtifactHandle, operation,
        restore: restoreCommittedArtifact, isCurrent: isCurrentAlert, isCurrentOperation });
    } finally {
      if (activeLosatAbortController === generationAbortController) {
        activeLosatAbortController = null;
      }
      if (generationToken === latestGenerationToken) {
        generationCancelRequested.value = false;
      }
    }
  };

  /**
   * @param {{ canonical: Record<string, any>, label?: string,
   *   captureIntentCheckpoint?: ((...args: any[]) => any) | null, restoreIntentCheckpoint?: ((...args: any[]) => any) | null,
   *   commitIntent?: (() => any) | null, alignmentResetBefore?: any, alignmentResetReceipt?: any, operation?: string }} options
   */
  const runCommittedCanonicalCandidate = async ({
    canonical,
    label = 'Rotate record to feature',
    captureIntentCheckpoint = null,
    restoreIntentCheckpoint = null,
    commitIntent = null,
    alignmentResetBefore = null,
    alignmentResetReceipt = undefined,
    operation = 'generate'
  }) => {
    const drawing = committedDrawing(canonical);
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    /** @type {Record<string, any> | null} */
    let outcome = null;
    const operationId = ++latestOperationId;
    const isCurrentOperation = () => operationId === latestOperationId;
    const previousAlert = errorLog.value;
    const isCurrentAlert = () => isCurrentOperation()
      && (errorLog.value === previousAlert || errorLog.value === null);
    /** @type {Record<string, any> | null} */
    let beforeHandle = null;
    const initialResults = results.value;
    /** @type {string | null} */
    let historyRecovery = null;
    processing.value = true;
    processingStatus.value = 'Preparing target record...';
    generationCancelRequested.value = false;
    try {
      const decorationContinuity = captureDecorationContinuity(getCommittedCanonicalSession?.(), projectCompositionRecordIdentity, drawing);
      await nextTick();
      await waitForAfterPaint();
      if (!isCurrentOperation()) return { status: 'stale' };
      const execute = (handle) => {
        beforeHandle = handle;
        return runCommittedCanonicalCandidateInternal(drawing, {
          canonical, decorationContinuity, generatedArtifactHandle: beforeHandle,
          commitIntent, alignmentResetBefore, alignmentResetReceipt, operation, isCurrentOperation, isCurrentAlert
        });
      };
      outcome = typeof runGeneratedArtifactReplacement === 'function'
        ? await runGeneratedArtifactReplacement(label, execute, {
            shouldCommit: (result) => result?.status === 'ok',
            captureIntentCheckpoint,
            restoreIntentCheckpoint,
            onCheckpointCapture: onGeneratedArtifactCheckpointCapture,
            restoreAppliedArtifact: async (beforeHandle) => {
              historyRecovery = 'restore-failed';
              await previewRuntime.restorePreviousSelectedResult({
                handle: beforeHandle,
                phase: 'target-history-finalization-rollback',
                restore: () => generatedArtifactTransactionOwner.restore(beforeHandle)
              });
              historyRecovery = beforeHandle.ownerSet.results.length ? 'restored' : 'no-result';
            }
          })
        : await execute(await captureGeneratedArtifactHandle());
      if (outcome?.status === 'ok' && outcome.generatedArtifactCandidate) {
        generatedArtifactTransactionOwner.finalize();
      }
      if (outcome && Object.prototype.hasOwnProperty.call(outcome, 'generatedArtifactCandidate')) {
        const { generatedArtifactCandidate, ...publicOutcome } = outcome;
        outcome = publicOutcome;
      }
      return outcome;
    } catch (cause) {
      outcome = await failOperation(cause, { handle: beforeHandle, operation, isCurrent: isCurrentAlert, isCurrentOperation,
        recovery: cause?.artifactRestoreFailed ? 'restore-failed'
          : historyRecovery || (!beforeHandle && results.value === initialResults && initialResults.length ? 'preserved' : null) });
      return outcome;
    } finally {
      if (isCurrentOperation()) {
        if (outcome?.status !== 'canceled') processingStatus.value = '';
        generationCancelRequested.value = false;
        processing.value = false;
        replayDeferredLabelReflow();
      }
    }
  };

  const cancelRunAnalysis = () => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const canceledGenerationToken = latestGenerationToken;
    latestGenerationToken += 1;
    generationCancelRequested.value = true;
    previewRuntime.invalidateReadinessExpectation(
      String(canceledGenerationToken),
      new DiagramGenerationCanceledError()
    );
    if (activeLosatAbortController && !activeLosatAbortController.signal.aborted) {
      activeLosatAbortController.abort(new DiagramGenerationCanceledError());
    }
    if (processing.value) {
      processingStatus.value = 'Canceling generation...';
    }
    return cancelDiagramGeneration();
  };

  // A reflow keeps the Result it draws, so the preview binder reports a forced
  // label that the Result does not draw instead of failing the binding, with
  // the diagnostic Generate raises for it (R6, OV-06). The rerender also draws
  // the Legend again, so the Legend it draws is the generated inventory, as
  // after Generate (OV-42, OV-43): the next edit compares against it.
  const expectReflowBindings = (commit, resultIndex, isCurrentReflow, diagramOptions) => {
    const result = commit.results[resultIndex];
    const featureIds = forcedLabelFeatureIds(commit.mutationPlan?.operationsByResult?.[resultIndex], {
      features: [...(commit.featureState?.renderedFeaturesByResult?.[resultIndex]?.values()
        || extractedFeatures.value || [])],
      diagramOptions
    });
    if (!result) return;
    const readinessId = `label-reflow:${latestGenerationToken}`;
    const resultIdentity = previewRuntime.getResultIdentity(result);
    previewRuntime.registerReadinessExpectation({
      result,
      resultIndex,
      artifactIdentity: readinessId,
      generationToken: readinessId,
      catalogState: featureCatalog?.value || null,
      phase: 'label-reflow',
      bindingOptions: {
        replaceGeneratedLegend: true,
        ...(featureIds.length === 0 ? {} : {
          reportedLabelBinding: Object.freeze({
            featureIds,
            report: (error) => {
              if (!isCurrentReflow()) return;
              labelReflowLastError.value = liveEditFailure(formatError(error, 'generate', 'render'));
            }
          })
        })
      },
      isCurrent: () => (
        previewRuntime.getResultIdentity(results.value[resultIndex]) === resultIdentity
        && Number(selectedResultIndex.value) === resultIndex
      )
    });
  };

  // A label reflow re-renders the committed Session with the current editor
  // tables (R1(c), N-16): it never reads the settings draft.
  /** @param {DrawingState} drawing */
  const runLabelReflowCandidate = async (drawing, { requestId, decorationContinuity }) => {
    if (mode.value === 'circular' && shouldDeferCircularPreviewUpdates.value) {
      return { status: 'skipped' };
    }
    const committed = getCommittedCanonicalSession?.();
    if (!committed) {
      labelReflowLastError.value = liveEditFailure(formatError(diagnosticError('LIVE_EDIT_REQUIRES_GENERATE')));
      return { status: 'skipped' };
    }
    // The reflow draws the committed Session's Results again, so they keep their names (OV-136).
    const committedResultNames = results.value.map((/** @type {{ name: string }} */ result) => result.name);
    const generationToken = ++latestGenerationToken;
    /** @type {Awaited<ReturnType<typeof prepareAndAdmitCandidate>>} */
    let colorCandidate = null;
    clearLabelBuildNotices({ rerender: true });
    skipCaptureBaseConfig.value = true;
    const isCurrent = () => generationToken === latestGenerationToken && requestId === pendingReflowRequestId;
    try {
      colorCandidate = await prepareAndAdmitCandidate(drawing, isCurrent);
      if (!colorCandidate) return { status: 'stale' };
      const candidateRules = colorCandidate.rules;
      const canonical = projectCommittedEditorIntent({
        committed,
        promotion: {
          featureCatalog: featureCatalog?.value ?? null,
          legacyOrthogroupState: { groups: cloneJsonData(orthogroups.value || []) }
        },
        state,
        drawing: {
          ...drawing,
          selectedPalette: appliedPaletteName,
          currentColors: appliedPaletteColors,
          manualSpecificRules: candidateRules
        }
      });
      const timingEntries = [];
      const execution = await executeCanonicalRenderCandidate({
        canonical,
        decorationContinuity,
        mode: canonical.renderRequest.mode,
        kind: 'reflow',
        shouldAdmit: colorCandidate.shouldAdmit,
        prepareCommit: prepareReflowResultCommit,
        prepareCommitInput: {
          featureColorOverrides: colorCandidate.featureColorOverrides,
          featureStrokeOverrides: drawing.featureStrokeOverrides,
          featureOverrides: drawing.featureOverrides,
          legendEntries: drawing.legendEntries.value,
          deletedLegendEntries: drawing.deletedLegendEntries.value,
          dormantLegendEntries: drawing.dormantLegendEntries.value,
          originalLegendOrder: originalLegendOrder.value,
          addedLegendCaptions: drawing.addedLegendCaptions.value,
          unrequestedDepthCaptions: unrequestedDepthCaptions(drawing, canonical),
          legendColorOverrides: drawing.legendColorOverrides,
          legendStrokeOverrides: drawing.legendStrokeOverrides,
          manualSpecificRules: candidateRules
        },
        timingEntries,
        resultNames: committedResultNames
      });
      console.info(`gbdraw ${canonical.renderRequest.mode} typed request render: ${formatDuration(execution.elapsedMs)}.`);
      if (execution.status === 'superseded' || !isCurrent()) return { status: 'stale' };
      if (execution.status === 'engine-error') {
        logPostGbdrawTimings(timingEntries);
        const error = formatError(execution.engineError, 'generate', 'render');
        labelReflowLastError.value = liveEditFailure(error);
        return { status: 'error', error };
      }
      const previousSelectedResultIndex = selectedResultIndex.value;
      const nextSelectedResultIndex = Math.max(
        0, Math.min(previousSelectedResultIndex, execution.commit.results.length - 1)
      );
      skipCaptureBaseConfig.value = true;
      // The rerender's catalog describes the Results it draws, so they replace
      // the previous pair together: a feature it draws again or no more keeps
      // its popup, label binding, Features list row, and History projection
      // (R-5). It renders the committed request, so the records and the
      // similarity groups the Generate committed stay.
      featureCatalog.value = execution.catalog;
      extractedFeatures.value = execution.commit.featureState.extractedFeatures;
      biologicalFeatures.value = execution.commit.featureState.biologicalFeatures;
      results.value = execution.commit.results;
      // The preview owner selects the Result (R13). It registers readiness only
      // when fewer Results leave the selection out of range; the reflow's own
      // expectation, registered next, replaces it, so the binder reports a
      // forced label the Result does not draw (OV-06).
      if (execution.commit.results.length > 0) {
        previewRuntime.selectResult(nextSelectedResultIndex);
      }
      expectReflowBindings(
        execution.commit, nextSelectedResultIndex, isCurrent, canonical.renderRequest.diagramOptions
      );
      logPostGbdrawTimings(timingEntries);
      colorCandidate.notifyChanges();
      return { status: 'ok' };
    } catch (e) {
      if (isDiagramGenerationCanceled(e)) {
        labelReflowLastError.value = null;
        return { status: 'canceled' };
      }
      const error = formatError(e, 'generate', 'render');
      labelReflowLastError.value = liveEditFailure(error);
      return { status: 'error', error };
    }
  };

  const runLabelReflow = async () => {
    const drawing = committedDrawing(getCommittedCanonicalSession?.());
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    pendingReflowRequestId += 1;
    if (activeReflowRequestId !== 0) return;
    if (processing.value) {
      reflowDeferredByProcessing = true;
      return;
    }

    labelReflowProcessing.value = true;
    try {
      // Generate takes priority: an iteration never starts while it runs.
      while (!processing.value && activeReflowRequestId < pendingReflowRequestId) {
        activeReflowRequestId = pendingReflowRequestId;
        let decorationContinuity;
        try {
          decorationContinuity = captureDecorationContinuity(getCommittedCanonicalSession?.(), projectCompositionRecordIdentity, drawing);
        } catch (error) {
          labelReflowLastError.value = liveEditFailure(formatError(error));
          return;
        }
        await runLabelReflowCandidate(drawing, {
          decorationContinuity,
          requestId: activeReflowRequestId
        });
      }
      if (processing.value && activeReflowRequestId < pendingReflowRequestId) {
        reflowDeferredByProcessing = true;
      }
    } finally {
      activeReflowRequestId = 0;
      labelReflowProcessing.value = false;
    }
  };

  // Called where a run clears `processing`: the requests it held back run once.
  const replayDeferredLabelReflow = () => {
    if (!reflowDeferredByProcessing || processing.value) return;
    reflowDeferredByProcessing = false;
    void runLabelReflow();
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
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    completedLosatSearch = null;
    losatCache.value = new Map();
    losatDerivedCache.value = new Map();
    proteinIdentityManifest.value = emptyProteinIdentityManifest();
    legacyProteinRawCandidates.value = { schema: 1, entries: [] };
    legacyProteinDerivedEvidence.value = { schema: 1, entries: [] };
    losatCacheInfo.value = [];
  };

  return {
    runAnalysis,
    runCommittedCanonicalCandidate,
    projectCommittedRecordTransform,
    projectCommittedSimilarityAlignment,
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
