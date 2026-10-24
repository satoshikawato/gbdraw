// @ts-check
/** @import { DrawingState, LinearSeq, WebState } from '../state.js' */
/** @import { InputFileRef } from '../services/history-files.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { SpecificColorRule } from '../services/specific-color-rules.js' */
/** @import { RulePreparation } from './rule-matching.js' */
/** @import { PreviousResultRestoreOptions, ReadinessExpectation, ReadinessExpectationOptions, ReadyReceipt } from './preview-runtime.js' */
/** @import { AskChoice, ShowNotice } from './app-messages.js' */
/** @import { ParsedSequenceRecord } from '../services/record-discovery.js' */
/** @import { AlignmentResetCanonical, AlignmentResetPlannedCanonical } from '../services/session-active-config-contract.js' */
/** @import { LosatSourceBatch } from '../services/linear-sources.js' */
/** @import { CommittedCanonicalSession } from '../services/session-request.js' */
/** @import { ConvertedResource } from '../services/main-session-comparison-frame.js' */
/** @import { RedrawnArtifact } from '../services/history-snapshot.js' */
/** @import { GenerationPayload } from '../services/diagram-generation.js' */
import { validateAnnotationWarnings } from '../services/session-feature-metadata.js';
import { msg, t, tm } from '../utils/i18n.js';
import { validateComparisonWarnings } from '../services/comparison-warnings.js';
import { draftDrawsDepth } from '../services/track-legend-rows.js';
import { normalizeResultLabelFits } from '../services/label-fit.js';
import { keepResultNames } from '../services/result-normalization.js';
import { regionWithinRecord } from '../services/record-options.js';
import {
  nameFeaturePlacementFailure,
  requestFeatureOverrides,
  validateFeatureIdentityNotices
} from '../services/feature-placement.js';
import { rekeyOrthogroupOverrides } from '../services/orthogroup-feature-metadata.js';
import { resolveLinearRegionBounds } from '../services/feature-metadata-extraction.js';
import { buildSimilarityAlignmentResetReceipt, validateSimilarityAlignmentResetReceipt } from '../services/session-active-config-contract.js';
import { prepareLosatRuntime, runLosatPairsParallel } from '../services/losat.js';
import {
  losatRecordGencode, losatSearchSequences, prepareLosatSourceBatches, splitLosatSourceResult
} from '../services/linear-sources.js';
import { offRecordRequestKeys } from '../services/record-draw-selection.js';
import {
  cancelDiagramGeneration,
  DIAGRAM_HELPER_OPERATIONS,
  DiagramGenerationCanceledError,
  isDiagramGenerationCanceled,
  runDiagramGeneration,
  resourceIdOwnedBy,
  runDiagramHelperOperation,
  selectHeldRenderState
} from '../services/diagram-generation.js';
import {
  projectCompositionRecordIdentity,
  buildCanonicalRenderRequest,
  bindCanonicalTypedResource,
  isCanonicalTypedResourceBound,
  committedFeatureVisibilityMatches,
  linearRecordSourceResources,
  projectCommittedEditorIntent,
  projectCommittedRecordTransform,
  projectLiveSettingWrites,
  projectCommittedSimilarityAlignment,
  readCanonicalResourceRecordCount,
  requestLabelProjection,
  requireLinearComparisonPlanSnapshot
} from '../services/session-request.js';
import { labelDrawingBlocker } from './feature-editor/label-actions.js';
import { applyCircularSuppressControlsToSlots } from './circular-track-slots.js';
import {
  clampCircularTrackAxisIndex,
  placeCircularTrackStack
} from '../services/circular-track-slot-model.js';
import {
  normalizeFileList,
  orderedConservationSources
} from '../services/conservation-series.js';
import {
  applyLinearTrackOrderPlacements,
  clampLinearTrackAxisIndex,
  normalizeLinearTrackSlots,
  resolveLinearTrackAxisIndex
} from '../services/linear-track-slot-model.js';
import {
  activeDepthTrackIndices,
  depthFileSlotsFromValue,
  depthTrackCoverageCount,
  depthTrackMatrixWidth,
  isRecordMajorDepthFileMatrix,
  normalizeRecordMajorDepthFileRows,
  depthSeriesLegendCaptions,
  representativeDepthFiles,
  shapeDepthTracksForGenerate,
  syncDepthSlotLabels
} from '../services/depth-track-state.js';
import {
  CustomTrackPlanValidationError,
  customTrackPlanIssues,
  validateCustomTrackPlan
} from '../services/track-slot-validation.js';
import { buildRunInfo, buildSourceRecipe, summarizeLosatRuntimes } from './run-info.js';
import {
  buildLosatJobSpecs,
  canonicalGridDisplayPairs,
  numberLosatDisplayPairs,
  resolveLinearComparisonPlan
} from '../services/linear-comparisons.js';
import {
  pruneUnmatchedFeatureOverrides,
  serializeFeatureVisibilityRules
} from '../services/feature-visibility.js';
import {
  normalizeDefinitionLineStyleState,
  refusedDefinitionLineFontSizePath
} from '../services/definition-line-style-state.js';
import { requireLinearLabelVisibilityMode } from '../services/linear-label-visibility.js';
import { createZipBlob } from '../utils/zip.js';
import { classifyOptionalPositiveNumber } from '../utils/optional-positive-number.js';
import { cloneJsonData, cloneJsonValue } from '../services/json-clone.js';
import { base64ToBytes, bytesToBase64, bytesToText, sha256Hex, textToBytes } from '../services/byte-utils.js';
import { readConvertedResources, withConvertedResources } from '../services/main-session-comparison-frame.js';
import { askBeforeLargeDownload, downloadBlob, downloadTextFile } from '../services/text-download.js';
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
  paralogLinksFor,
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
  discoveryErrorIsFinal,
  gffFastaReaderSelector,
  parseSequenceRecordsFast,
  readAllSequenceRecords,
  readComparisonSequence,
  readComparisonSequenceDigest,
  readOncePerFile,
  readParsedSequenceRecords
} from '../services/record-discovery.js';
import {
  classifyRawLosatCacheEntry,
  emptyProteinIdentityManifest,
  getCurrentRawLosatCacheEntry,
  isCurrentRawLosatCacheEntry,
  isDerivedAnalysisResult,
  isLosatDerivedCacheEntry,
  mergeProteinIdentityManifests,
  normalizeLosatArgs,
  sameLosatArgs,
  transitionLegacyProteinCandidate,
  losatEdgeFilename,
  validateDerivedProteinReferences,
  webLosatRuntimeRecord,
  writtenLegacyProteinArtifacts
} from '../services/losat-cache.js';
import {
  LEGACY_LOSAT_DERIVED_CACHE_SCHEMA, LEGACY_PROTEIN_CANDIDATE_SCHEMA,
  LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA,
  LOSAT_DERIVED_CACHE_SCHEMA,
  NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
  ORTHOGROUP_INFERENCE_VERSION,
  PROTEIN_LOSAT_CACHE_SCHEMA
} from '../contracts/schema.generated.js';
import { comparisonFiltersForMode, resolveComparisonThresholds } from '../mode-profiles.js';
import { diagnosticError, liveEditFailure, normalizeUserFacingError } from '../utils/error-normalization.js';
import { cloneFileBytesForTransfer } from '../services/file-content-cache.js';
import {
  admitFeatureCatalog
} from '../services/feature-catalog.js';
import {
  prepareCandidateRenderCommit,
  prepareReflowResultCommit
} from './candidate-render.js';
import { editorPaintState } from './result-paint-record.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from '../services/runtime-test-hooks.js';
import { GROUP_NAMES_SLICE } from '../services/drawing/group-names.js';
import { FEATURE_EDITS_SLICE } from '../services/drawing/feature-edits.js';
import { replaceDormantLegendRows, replaceFeatureFills } from '../services/drawing/paint.js';
import {
  IMPORTED_COMPARISON_DISPOSITIONS,
  createImportedComparisonIntentState,
  inheritCommittedComparisonIntent
} from '../services/imported-comparison-intent.js';
/** @import { LabelDiagramOptions, LabelFeature } from './feature-editor/label-actions.js' */
/** @import { SvgMutationOperations } from '../services/svg-result-ingestion.js' */
/** @import { FileLike } from '../services/session-resource-backing.js' */
/** @import { FeatureCatalog } from '../services/feature-catalog.js' */
/** @import { LosatJob, LosatJobResult, LosatRunOptions } from '../services/losat.js' */
/** @import { OrderedConservationSource } from '../services/conservation-series.js' */
/** @import { LosatJobSpec } from '../services/linear-comparisons.js' */
/** @import { FeatureCatalogRow } from '../services/feature-catalog.js' */
/** @import { LegacyProteinCandidateEnvelope, LosatDerivedCacheEntry, RawLosatCacheEntry, RawLosatCacheMatch } from '../services/losat-cache.js' */
/** @import { RecordKeyedFeature } from '../services/feature-identity.js' */
/** @import { RecipeGeneratedFile, RecipeHelperFile } from './run-info.js' */
/** @import { LinearTrackSlot } from '../services/linear-track-slot-model.js' */
/** @import { DepthSlotRow, RawDepthTrackConfig } from '../services/depth-track-state.js' */
/** @import { UserFacingError } from '../utils/error-normalization.js' */
/** @import { HistoryArtifactHandle, HistoryStepOptions } from '../services/history.js' */
/** @import { FeaturePlacementDraft } from '../services/feature-placement.js' */
/** @import { EditorPaintState } from './result-paint-record.js' */
/** @import { OrthogroupGroup } from '../services/orthogroup-feature-metadata.js' */
/** @import { CircularConservationLike } from '../services/conservation-series.js' */
/** @import { FeatureFill, LegendRow } from '../services/drawing/paint.js' */
/** @import { LinearComparisonResolvedEdge } from '../services/linear-comparisons.js' */
/** @import { RequestCanonicalComparison } from '../services/session-request.js' */
/** @import { CanonicalFeatureOverride, CanonicalGeneratedProteinComparison, CanonicalLosslessCollinearityParameters, CanonicalRecordSource, CanonicalRenderRequest } from '../contracts/request.generated.js' */

const DEFAULT_LINEAR_BLAST_FILTERS = Object.freeze(
  comparisonFiltersForMode('linear')
);

const DIAGRAM_PROGRESS_MESSAGES = Object.freeze({
  'preparing-runtime': msg('Preparing diagram runtime (first use)...'),
  'preparing-resources': msg('Preparing diagram input resources...'),
  rendering: msg('Rendering diagram...'),
  finalizing: msg('Finalizing diagram results...')
});
// The status line of a diagram Worker stage, shown in the processing overlay
// while it runs. The first-use runtime preparation names its real step:
// "step 2 of 5, loading the package installer" (the Worker marks the step
// names with `msg`).
/** @param {{ stage?: string, step?: number, total?: number, label?: string }} progress */
const diagramProgressMessage = ({ stage = '', step, total, label } = {}) => {
  const message = DIAGRAM_PROGRESS_MESSAGES[/** @type {keyof typeof DIAGRAM_PROGRESS_MESSAGES} */ (stage)] || '';
  return stage === 'preparing-runtime' && step && total && label
    ? t('Preparing diagram runtime (first use): step {step} of {total}, {label}...', { step, total, label: tm(label) })
    : tm(message);
};

export const afterPaint = ({ requestFrame = globalThis.requestAnimationFrame } = {}) => {
  if (typeof requestFrame !== 'function') return Promise.resolve();
  return new Promise((resolve) => {
    requestFrame(() => requestFrame(resolve));
  });
};

// One task boundary, no frame (PD-OI-101). Chromium dispatches input that
// arrived during the work before it first, so a Cancel clicked then is seen
// before the run goes on; the status set just before may not paint.
/** @returns {Promise<void>} */
export const afterTask = () => new Promise((resolve) => { setTimeout(resolve, 0); });

export const afterFrame = ({ requestFrame = globalThis.requestAnimationFrame } = {}) => {
  if (typeof requestFrame !== 'function') return Promise.resolve();
  return new Promise((resolve) => requestFrame(resolve));
};

// The SHA-256 of a text's UTF-8 bytes, from the one SHA-256 owner (OV-1566). Generate needs Web Crypto anyway:
// the Worker digests every reply (`buildGeneratedArtifactTransportIdentity`), so no weaker fallback key exists.
/** @param {string} text */
const hashText = (text) => sha256Hex(textToBytes(text));

// Label visibility On draws the label whatever Show Labels and the label
// filters select, so it must bind, unless the request cannot draw it: a feature
// drawn as underlay, or Label Rendering = Embedded Only. That On takes effect
// when the label can be drawn (Owner Q1, Q2). The operations reach only the
// features the request's catalog lists as rendered, which Python draws
// (should_render_feature), so a hidden feature has none. Label text alone and
// Off follow those settings, which may leave the feature unlabeled.
/** @param {readonly ({ renderedId?: unknown } | null | undefined)[] | null | undefined} operations */
const labelOperationIds = (operations) => (operations || [])
  .map((operation) => String(operation?.renderedId || '').trim())
  .filter(Boolean);
/**
 * @param {Pick<SvgMutationOperations, 'labelVisibility'> | null | undefined} operations
 * @param {{ features?: readonly (LabelFeature | null | undefined)[] | null, diagramOptions: LabelDiagramOptions | null | undefined }} context
 */
const forcedLabelFeatureIds = (operations, { features, diagramOptions }) => {
  const featuresById = new Map((features || []).map((feature) => [String(feature?.svg_id || '').trim(), feature]));
  return Object.freeze([...new Set(labelOperationIds(
    (operations?.labelVisibility || []).filter((operation) => operation?.mode === 'on')
  ))].filter((featureId) => !labelDrawingBlocker(featuresById.get(featureId), diagramOptions)));
};

const getNow = () => (typeof globalThis.performance?.now === 'function' ? performance.now() : Date.now());
/** @param {number} ms */
const formatDuration = (ms) => `${(ms / 1000).toFixed(2)}s`;
/**
 * One record's FASTA text as the LOSAT search reads it: the browser or Worker extraction, and for a protein
 * search the Worker's protein map and identity fields.
 * @typedef {object} LosatSeqEntry
 * @property {string} fasta
 * @property {string} recordId
 * @property {number} canonicalLength
 * @property {Record<string, unknown>} [proteinMap]
 * @property {number} [proteinCount]
 * @property {string} [proteinSetHash]
 * @property {string} [recordAnalysisId]
 * @property {string} [recordInstanceKey]
 * @property {string} [runtimeBindingHash]
 * @property {string} [displayBindingHash]
 * @property {unknown} [identityManifest]
 * @property {string} [hash] The SHA-256 of `fasta`, which a protein extraction carries.
 */
/**
 * Values kept per owner object (an admitted, frozen manifest) and the settings of the call that built them.
 * @template T
 * @typedef {{ byOwner: WeakMap<object, Map<string, T>>, limit: number }} FileReadCache
 */
// What Generate reads from an input File is kept per File in `services/record-discovery.js` (`readOncePerFile`,
// OV-1587): each sequence file parsed once by format (OV-1556), the LOSAT FASTA of a record and region and,
// once a search needs it, its SHA-256, the blastp protein extraction, and the Circular similarity ring reference
// with its SHA-256 (OV-1545, OV-1565). Each held value is frozen: its readers never write to it (OV-1594).
/** The blastp protein extractions kept per source File (or GFF3 and FASTA pair). */
const PROTEIN_EXTRACTIONS_PER_FILE = 16;
/** The LOSAT FASTA texts, and their digests, kept per source File: one per record and region a card names. */
const LOSAT_FASTA_READS_PER_FILE = 64;
/**
 * The checked BUILD_PROTEIN_LOSAT_CACHE_KEYS replies of an admitted merged manifest, by the rest of the request
 * (its pairs and derived records), which with the manifest is every input of the reply (OV-1553).
 * @type {FileReadCache<Readonly<{ keys: readonly string[], derivedKey: string }>>}
 */
const proteinCacheKeysCache = { byOwner: new WeakMap(), limit: 8 };
const LOSAT_DERIVED_CACHE_LIMIT = 16;
export const LOSAT_EXPORT_CONFIRM_THRESHOLD_BYTES = 50 * 1024 * 1024;

/** @param {unknown} hydratedResults */
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

/** @param {unknown} totalBytes */
export const losatExportBytesExceedLimit = (
  totalBytes,
  limitBytes = LOSAT_EXPORT_CONFIRM_THRESHOLD_BYTES
) => Number(totalBytes) > Number(limitBytes);

/**
 * @typedef {object} LosatCacheMetadata
 * @property {string} [identityKind] `protein` for a BLASTP search, otherwise nucleotide.
 * @property {string} [flow]
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
 * @property {string} [searchContext] The source batch's search context.
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
  /** @type {{ cacheSchema: number, program: string, outfmt: string, args: string[], queryCanonicalHash?: string, subjectCanonicalHash?: string, flow?: string, searchContext?: string }} */
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

/** @import { LosatCacheInfoEntry, MultiRecordPosition } from '../services/config.js' */
/** @import { LiveSettingWrite } from '../services/session-request.js' */

/**
 * A LOSAT search Generate queues: the fields the pool reads (`LosatJob`) and what Generate keeps on it to
 * hash and name the result.
 * @typedef {LosatJob & {
 *   queryCanonicalHash?: string, subjectCanonicalHash?: string, outfmt?: string, extraArgs?: readonly unknown[]
 * }} LosatQueuedJob
 */

/**
 * A Linear comparison search: the job of one resolved pair, with the batch that runs it.
 * @typedef {LosatQueuedJob & {
 *   cacheKey: string, ordinal: number, edgeKey: string, queryUid: string, subjectUid: string, queryIndex: number,
 *   subjectIndex: number, extraArgs: readonly unknown[], cacheMetadata: LosatCacheMetadata,
 *   batch: LosatSourceBatch<LosatJobSpec>
 * }} LosatLinearJob
 */

/**
 * A LOSAT search of one source-file pair: the first job of its batch, with the record pairs it covers.
 * @typedef {LosatLinearJob & {
 *   queryRecordIndexes: readonly number[], subjectRecordIndexes: readonly number[],
 *   recordPairs: number[][], members: LosatLinearJob[]
 * }} LosatSourceJob
 */

// The Worker conversion of Generate's ring LOSAT rows (`circular_ring_tables_in_search_frame_json`).
const RING_TABLE_CONVERSION = /** @type {const} */ ({ stage: 'helper', operation: 'convertRingTablesToSearchFrame' });

/**
 * The ring tables of a request in the search frame of its drawn records, by ring table resource ID. Python
 * converts each ring's whole-file LOSAT rows once per raw search and frame: the conversion is kept on the reference
 * File (or GFF3 and FASTA pair), whose bytes never change, under the ring's raw cache key (which names immutable
 * search output) and the records' selector, crop and orientation. Each ring keeps its last frame only, and a failed
 * conversion is not kept, nor one of a ring without current search output (`cacheKey` null).
 * @param {readonly [FileLike] | readonly [FileLike, FileLike]} referenceFiles The Files the request's records read.
 * @param {{ renderRequest: unknown, records: readonly unknown[], resources: GenerationPayload['resources'] }} request
 * @param {readonly { resourceId: string, cacheKey: string | null }[]} rings
 * @returns {Promise<Map<string, ConvertedResource>>}
 */
const convertRingTablesOncePerFrame = async (referenceFiles, { renderRequest, records, resources }, rings) => {
  const frameKey = JSON.stringify(records);
  /** @type {string[]} */
  const missing = [];
  /** @type {(resourceIds: string[]) => void} */
  let startConversion = () => {};
  const conversion = new Promise((resolve) => { startConversion = resolve; }).then(
    /** @param {string[]} resourceIds */
    async (resourceIds) => (resourceIds.length === 0 ? new Map() : readConvertedResources(
      await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.CONVERT_RING_TABLES_TO_SEARCH_FRAME, {
        resourceIds, canonicalRequest: renderRequest, resources
      }),
      RING_TABLE_CONVERSION
    ))
  );
  const tables = rings.map(({ resourceId, cacheKey }, index) => {
    const convert = async () => {
      missing.push(resourceId);
      const table = (await conversion).get(resourceId);
      if (!table) throw diagnosticError('RESULT_INVALID', {}, RING_TABLE_CONVERSION);
      return Object.freeze(table);
    };
    return cacheKey === null
      ? convert()
      : readOncePerFile(referenceFiles, `ring-search-frame-${index}:${cacheKey}\n${frameKey}`, convert, 1);
  });
  startConversion(missing);
  const converted = await Promise.all(tables);
  return new Map(rings.map(({ resourceId }, index) => [resourceId, converted[index]]));
};

/**
 * One circular conservation comparison, resolved to the BLAST text the request carries.
 * @typedef {{ name: string, text: string, sourceIndex: number, label: string, color: string }} ResolvedCircularConservation
 */

/**
 * A match in one of Generate's cache maps, whose keys are text (`RawLosatCacheMatch` of services/losat-cache.js
 * leaves the key open: its reader takes any `Map`).
 * @typedef {RawLosatCacheMatch & { key: string }} TextKeyedLosatCacheMatch
 */

/**
 * @param {Map<string, RawLosatCacheEntry> | null | undefined} cacheMap
 * @param {string} cacheKey
 * @param {LosatCacheMetadata} metadata
 * @param {Record<string, any> | null} [manifest]
 * @returns {TextKeyedLosatCacheMatch | null}
 */
const getRawLosatCacheEntry = (cacheMap, cacheKey, metadata, manifest = null) => {
  if (!cacheMap) return null;
  const direct = getCurrentRawLosatCacheEntry(cacheMap, cacheKey, metadata, manifest);
  // `getCurrentRawLosatCacheEntry` returns the `cacheKey` it was given (its `unknown` parameter widens the type).
  if (direct) return /** @type {TextKeyedLosatCacheMatch} */ (direct);
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

/**
 * @param {Map<string, RawLosatCacheEntry> | null | undefined} cacheMap
 * @param {string} cacheKey
 * @param {TextKeyedLosatCacheMatch | null | undefined} found
 * @param {LosatCacheMetadata} metadata
 */
const promoteRawLosatCacheEntry = (cacheMap, cacheKey, found, metadata) => {
  if (!cacheMap || !found?.entry) return found?.entry || null;
  if (classifyRawLosatCacheEntry(found.entry) !== 'nucleotide-current') {
    cacheMap.set(cacheKey, found.entry);
    return found.entry;
  }
  /** @type {RawLosatCacheEntry} */
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

/** @param {Map<string, unknown> | null | undefined} cacheMap */
const pruneLosatDerivedCache = (cacheMap) => {
  if (!cacheMap || typeof cacheMap.delete !== 'function') return;
  while (cacheMap.size > LOSAT_DERIVED_CACHE_LIMIT) {
    const oldestKey = cacheMap.keys().next().value;
    if (oldestKey === undefined) break;
    cacheMap.delete(oldestKey);
  }
};

/**
 * @param {Map<string, unknown> | null | undefined} cacheMap The derived cache; an entry is read as unchecked.
 * @param {string} key
 * @param {unknown} manifest
 */
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

/**
 * Whether a converted protein payload holds what a request needs from it: the typed result bound to its Python
 * bytes and each displayed pair's canonical comparison table. A payload read from a Session or written before
 * the table was canonical lacks one of them, so it is a cache miss and the analysis runs again.
 * @param {unknown} mode
 * @param {Record<string, unknown> | null | undefined} payload
 */
export const hasRequiredCanonicalAnalysisResource = (mode, payload) => {
  const normalizedMode = String(mode || '').trim().toLowerCase();
  if (!['orthogroup', 'collinear'].includes(normalizedMode)) return true;
  const resource = normalizedMode === 'collinear'
    ? payload?.collinearityResult
    : payload?.orthogroupResult;
  return Boolean(
    isDerivedAnalysisResult(normalizedMode, resource)
    && isCanonicalTypedResourceBound(resource)
    && Array.isArray(payload?.pairs)
    && payload.pairs.every((pair) => typeof pair?.canonicalTsv === 'string')
  );
};

/** @param {unknown} key */
const isHexDigest = (key) => typeof key === 'string' && /^[0-9a-f]{64}$/.test(key);

/**
 * The keys of the `buildProteinLosatCacheKeys` reply: one raw key per pair, and the key of the derived payload
 * they convert to when the request asked for it (`''` otherwise).
 * @typedef {{ keys: string[], derivedKey: string }} ProteinLosatCacheKeys
 */
// The reply of the `buildProteinLosatCacheKeys` helper is unchecked data. Generate uses each key as a raw LOSAT
// cache key (and the derived key as the derived cache key), so the reply holds SHA-256 hex digests, as text
// (`RegExp.test` alone would coerce).
/**
 * @param {unknown} reply
 * @param {number} pairCount
 * @param {{ derived?: boolean }} [options] Whether the request asked for the derived key.
 * @returns {ProteinLosatCacheKeys}
 */
export const requireProteinLosatCacheKeys = (reply, pairCount, { derived = false } = {}) => {
  /** @type {{ error?: unknown, keys?: unknown, derived?: { key?: unknown } | null }} */
  const result = reply !== null && typeof reply === 'object' ? reply : {};
  const keys = Array.isArray(result.keys) ? /** @type {unknown[]} */ (result.keys) : null;
  const derivedKey = derived ? result.derived?.key : '';
  if (
    result.error || !keys
    || keys.length !== pairCount
    || keys.some((key) => !isHexDigest(key))
    || (derived && !isHexDigest(derivedKey))
  ) {
    throw result.error || new Error('Protein cache key generation failed.');
  }
  // Every key passed the text and pattern check above.
  return { keys: /** @type {string[]} */ (keys), derivedKey: /** @type {string} */ (derivedKey) };
};

/** @param {unknown} candidateLimit */
export const resolveProteinBlastpCandidateLimit = (candidateLimit) => (
  requireCurrentProteinBlastpCandidateLimit(candidateLimit)
);

/** @param {unknown} comparisons */
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

/**
 * @param {unknown} left
 * @param {unknown} right
 */
const sameNumber = (left, right) => Number(left) === Number(right);

/**
 * @typedef {object} ResolvedProteinReuseActive The settings of this Generate that the committed Session must match.
 * @property {string} featureVisibility
 * @property {CanonicalFeatureOverride[]} featureOverrides The draft's Feature visibility rows of the compared records.
 * @property {string} mode
 * @property {string | number} [bitscore]
 * @property {string | number} [evalue]
 * @property {string | number} [identity]
 * @property {string | number} [alignmentLength]
 * @property {number | null} candidateLimit
 * @property {number} maxHits
 * @property {string} orthogroupMembershipMode
 * @property {number | null} memberMaxHits
 * @property {number} minAnchors
 * @property {number} maxUnitGap
 * @property {number} maxDiagonalDrift
 * @property {number} maxConflicts
 * @property {number} maxParalogLinks
 * @property {string} colorMode
 * @property {string} unitMode
 * @property {string} anchorMode
 * @property {string} mergeOrientation
 * @property {boolean} inferOrthogroups
 * @property {string} searchScope
 */

/**
 * @param {object} options
 * @param {RequestCanonicalComparison[]} options.canonicalComparisons The draft's `drawing.files.linearCanonicalComparisons`.
 *   The caller's value is untyped (`files` is a Vue reactive object), so this is not checked there; a hand-edited
 *   Session can still put a `null` entry in the list (OV-750), which the `?.` reads below tolerate.
 * @param {CommittedCanonicalSession | null} options.committedSession The committed Session: its request is of the
 *   current schema (Load reads a stored one to schema 9; Generate and the rerender commit the one they built).
 * @param {readonly LinearSeq[]} options.sequences
 * @param {ResolvedProteinReuseActive} options.active
 */
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
  // The test of `kind` is the type.
  const committedMarker = /** @type {CanonicalGeneratedProteinComparison | undefined} */ (committed.find((comparison) => (
    comparison?.kind === 'generatedProteinComparison' && comparison.mode === 'none'
  )));
  // A committed marker implies a committed request; the test only states it.
  if (!committedRequest || !persistedMarker || !committedMarker || !active) return false;
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
    const selector = record.region?.selector || record.selector;
    return record.recordKey !== seq.uid
      || String((selector && 'value' in selector ? selector.value : '') || '') !== String(seq.region_record_id || '').trim()
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
    String(settings.orthogroupMembershipMode || ORTHOGROUP_INFERENCE_VERSION)
      !== active.orthogroupMembershipMode
    || (settings.orthogroupMemberMaxHits ?? null) !== active.memberMaxHits
  ) {
    return false;
  }
  // Similarity groups and Collinear both draw with the paralog links per group (OV-1128 B).
  if (!sameNumber(settings.collinearMaxParalogLinksPerOrthogroup ?? 2, active.maxParalogLinks)) return false;
  if (active.mode === 'orthogroup') return true;

  /** @type {Partial<CanonicalLosslessCollinearityParameters>} */
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
  );
};

/** @param {Record<string, unknown>} payload */
const stripRuntimeCacheStats = (payload) => {
  // Canonical comparison data is read-only across rendering and cache reuse.
  // Only the transient statistics owner differs; do not duplicate the full result.
  const { cache: _cache, ...stored } = payload;
  return stored;
};

const LEGACY_PROTEIN_REFERENCE_RE = /p_[A-Za-z0-9._%+-]+?_\d+_\d+_(?:-1|0|1)_[0-9a-f]{12}(?:_[2-9][0-9]*)?/g;

/** @param {unknown[]} values */
const collectLegacyProteinReferences = (...values) => {
  /** @type {Set<string>} */
  const references = new Set();
  const visited = new WeakSet();
  /** @param {unknown} value */
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

// Rewrites strings and object keys and keeps the value's shape, so the result has the argument's type.
/**
 * @template T
 * @param {T} value
 * @param {Record<string, unknown> | null | undefined} idMap
 * @returns {T}
 */
const rewriteMappedProteinReferences = (value, idMap) => {
  if (!idMap || typeof idMap !== 'object' || Array.isArray(idMap)) return cloneJsonData(value);
  if (typeof value === 'string') {
    if (Object.prototype.hasOwnProperty.call(idMap, value)) return /** @type {T} */ (String(idMap[value]));
    return /** @type {T} */ (value.replace(
      LEGACY_PROTEIN_REFERENCE_RE,
      (reference) => (
        Object.prototype.hasOwnProperty.call(idMap, reference)
          ? String(idMap[reference])
          : reference
      )
    ));
  }
  if (Array.isArray(value)) {
    return /** @type {T} */ (value.map((item) => rewriteMappedProteinReferences(item, idMap)));
  }
  if (!value || typeof value !== 'object') return value;
  /** @type {Record<string, unknown>} */
  const rewritten = {};
  Object.entries(value).forEach(([key, item]) => {
    const nextKey = Object.prototype.hasOwnProperty.call(idMap, key) ? String(idMap[key]) : key;
    if (Object.prototype.hasOwnProperty.call(rewritten, nextKey)) {
      throw new Error(`Protein reference migration produced duplicate key '${nextKey}'.`);
    }
    rewritten[nextKey] = rewriteMappedProteinReferences(item, idMap);
  });
  return /** @type {T} */ (rewritten);
};

/**
 * @param {Map<string, LosatDerivedCacheEntry> | null | undefined} cacheMap
 * @param {string} key
 * @param {{ mode: unknown, payload: Record<string, unknown> | null | undefined, manifest: unknown }} entry
 */
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

/** @param {unknown} name */
const makeSafeFilename = (name) => {
  const cleaned = String(name || '').replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '');
  return cleaned || 'losat';
};
// The deferred replay builder outlives its Generate through the CLI helper
// files. Build it outside that scope so it retains only the published artifact,
// not the pre-Generate rollback handle and, through it, every earlier artifact.
/**
 * @param {{
 *   version: unknown, startedAtIso: unknown, renderRequest: unknown, resources: unknown,
 *   publishedArtifact: { results: unknown, featureCatalog: unknown }
 * }} replay Only `JSON.stringify` reads the first four.
 */
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
/** @param {unknown} value */
const normalizeRecordSelectorText = (value) => {
  const normalized = String(value ?? '').trim();
  if (!normalized || ['none', 'null', 'jsnull', 'undefined', 'jsundefined', '-'].includes(normalized.toLowerCase())) {
    return '';
  }
  return normalized;
};
/** @param {unknown} value */
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
/** @param {unknown} sequence */
const reverseComplementSequence = (sequence) => {
  /** @type {Record<string, string>} */
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
/** @param {string} sequence */
const wrapFastaSequence = (sequence) => {
  const lines = [];
  for (let i = 0; i < sequence.length; i += 60) lines.push(sequence.slice(i, i + 60));
  return lines.join('\n');
};
/**
 * A record's sequence as the LOSAT search reads it.
 * @typedef {{ id: string, sequence: string }} LosatSequence
 */
/** @param {LosatSequence} record */
const buildFastaText = (record) => `>${record.id}\n${wrapFastaSequence(record.sequence)}\n`;
/** @param {unknown} fasta */
const getFastaSequenceLength = (fasta) =>
  String(fasta || '')
    .split(/\r?\n/)
    .filter((line) => line && !line.startsWith('>'))
    .join('')
    .replace(/\s+/g, '').length;
/**
 * @param {readonly ParsedSequenceRecord[]} records
 * @param {unknown} selectorRaw
 */
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
/**
 * @param {ParsedSequenceRecord} record
 * @param {unknown} regionSpec
 * @returns {LosatSequence}
 */
const applyLosatSequenceTransforms = (record, regionSpec) => {
  if (record.sequence === null) throw new Error(`Record ${record.id} has no sequence for the browser reader.`);
  const region = parseRegionText(regionSpec);
  let sequence = record.sequence;
  if (region) {
    const drawn = regionWithinRecord(region.start, region.end, sequence.length);
    if (!drawn) {
      // The region starts after the record ends, as crop_genbank.py reports it (GX-09).
      throw diagnosticError('REGION_INVALID', { field: 'region', reason: 'RECORD_BOUNDS' });
    }
    sequence = sequence.slice(drawn.start - 1, drawn.drawnEnd);
    if (region.reverse) sequence = reverseComplementSequence(sequence);
  }
  return { id: record.id, sequence };
};
/**
 * @template T
 * @param {FileReadCache<T>} cache
 * @param {object} owner
 * @param {string} key
 * @returns {T | null}
 */
const getCachedFileRead = (cache, owner, key) => cache.byOwner.get(owner)?.get(key) ?? null;
/**
 * @template T
 * @param {FileReadCache<T>} cache
 * @param {object} owner
 * @param {string} key
 * @param {T} value
 */
const setCachedFileRead = (cache, owner, key, value) => {
  let byKey = cache.byOwner.get(owner);
  if (!byKey) {
    byKey = new Map();
    cache.byOwner.set(owner, byKey);
  }
  const oldestKey = byKey.keys().next().value;
  if (byKey.size >= cache.limit && oldestKey !== undefined) byKey.delete(oldestKey);
  byKey.set(key, value);
};
/**
 * @typedef {{ label: string, ms: number, details?: string }} TimingEntry
 */
/**
 * @template T
 * @param {TimingEntry[]} entries
 * @param {string} label
 * @param {() => T} fn
 * @returns {T}
 */
const measureTiming = (entries, label, fn) => {
  const startedAt = getNow();
  const result = fn();
  entries.push({ label, ms: getNow() - startedAt });
  return result;
};
/** @param {readonly TimingEntry[] | null | undefined} entries */
const logPostGbdrawTimings = (entries) => {
  if (!entries || entries.length === 0) return;
  console.groupCollapsed('post-gbdraw timing');
  entries.forEach(({ label, ms, details }) => {
    console.info(`${label}: ${formatDuration(ms)}${details ? ` (${details})` : ''}`);
  });
  console.groupEnd();
};
/**
 * The source of a Linear record as a request names it, with its resources, when the request draws these files
 * (a GenBank file, or a GFF3 file and its FASTA file); otherwise null.
 * @param {Readonly<Record<string, unknown>> | null | undefined} resources The committed Session's resources, or
 *   those the next Generate's request names (`linearRecordSourceResources`).
 * @param {FileLike} sourceFile
 * @param {FileLike | null} fastaFile
 * @returns {{ recordSource: CanonicalRecordSource, resources: NonNullable<GenerationPayload['resources']> } | null}
 */
const requestRecordSource = (resources, sourceFile, fastaFile) => {
  if (!resources) return null;
  /** @type {CanonicalRecordSource | null} */
  let recordSource = null;
  if (!fastaFile) {
    const resourceId = resourceIdOwnedBy(resources, sourceFile, 'genbank');
    if (resourceId) recordSource = { kind: 'genbank', resourceId };
  } else {
    const gffResourceId = resourceIdOwnedBy(resources, sourceFile, 'gff3');
    const fastaResourceId = resourceIdOwnedBy(resources, fastaFile, 'fasta');
    if (gffResourceId && fastaResourceId) recordSource = { kind: 'gffFasta', gffResourceId, fastaResourceId };
  }
  if (!recordSource) return null;
  const ids = recordSource.kind === 'genbank'
    ? [recordSource.resourceId] : [recordSource.gffResourceId, recordSource.fastaResourceId];
  // The transport validates each descriptor when it stages it.
  return {
    recordSource,
    resources: Object.fromEntries(ids.map((id) => [id, /** @type {NonNullable<GenerationPayload['resources']>[string]} */ (resources[id])]))
  };
};
/**
 * @param {ParsedSequenceRecord} record
 * @param {unknown} regionSpec
 * @returns {LosatSeqEntry}
 */
const losatSeqEntryFromRecord = (record, regionSpec) => {
  const transformed = applyLosatSequenceTransforms(record, regionSpec);
  return {
    fasta: buildFastaText(transformed),
    recordId: transformed.id,
    canonicalLength: transformed.sequence.length
  };
};
/**
 * The LOSAT FASTA of one record of a sequence text; Generate reads the same records through
 * `readParsedSequenceRecords`, and the Python twin is `losat_search_frame_fasta`.
 * @param {{ text: string, fmt: string, regionSpec: unknown, recordSelector: unknown }} options
 */
export const extractLosatFastaFast = async ({ text, fmt, regionSpec, recordSelector }) => (
  losatSeqEntryFromRecord(selectParsedRecord(parseSequenceRecordsFast(text, fmt), recordSelector), regionSpec)
);
/**
 * The LOSAT FASTA of every record of the Circular reference: the subject database of the conservation rings, the
 * records the CLI loader reads (OV-1599). A record without a sequence is left out.
 * @param {FileLike} file
 * @param {string} fmt
 */
const readAllLosatFasta = async (file, fmt) => {
  /** @type {LosatSequence[]} */
  const records = (await readAllSequenceRecords(file, fmt))
    .flatMap(({ id, sequence }) => (sequence === null ? [] : [{ id, sequence }]));
  if (!records.length) throw new Error('No records found for circular conservation reference.');
  return records.map((record) => buildFastaText(record)).join('');
};
/**
 * @param {InputFileRef | InputFileRef[] | null | undefined} sourceFiles
 * @param {CircularConservationLike | null | undefined} circularConservation
 */
const buildConservationSeries = (sourceFiles, circularConservation) => {
  return orderedConservationSources(sourceFiles, circularConservation).map((entry) => ({
    label: entry.label,
    color: entry.color
  }));
};
/** @param {unknown} value */
const normalizeBlastpMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['pairwise', 'orthogroup', 'collinear'].includes(normalized) ? normalized : 'orthogroup';
};
/** @param {unknown} value */
const normalizeCircularConservationLosatProgram = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return normalized === 'tblastx' ? 'tblastx' : 'blastn';
};
/**
 * @param {unknown} value
 * @param {number} [fallback]
 */
const normalizePositiveInteger = (value, fallback = 1) => {
  const parsed = Number(value);
  return Number.isInteger(parsed) && parsed > 0 ? parsed : fallback;
};
/**
 * @param {unknown} value
 * @param {{ maxRow?: number }} [options]
 */
const normalizeMultiRecordPositions = (value, { maxRow = Number.POSITIVE_INFINITY } = {}) => {
  if (!Array.isArray(value)) return [];
  /** @type {MultiRecordPosition[]} */
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
/** @param {MultiRecordPosition[]} positions */
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
/**
 * The default Record Order of a Circular source's records: a near-square grid.
 * @param {unknown} selectors
 * @returns {MultiRecordPosition[]}
 */
export const buildDefaultMultiRecordPositions = (selectors) => {
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
/**
 * @param {unknown} records
 * @param {unknown} currentPositions
 */
const mergeCircularRecordPositions = (records, currentPositions) => {
  const availableSelectors = Array.isArray(records)
    ? records.map((entry) => String(entry?.selector || '').trim()).filter(Boolean)
    : [];
  if (availableSelectors.length === 0) return [];
  const availableSet = new Set(availableSelectors);
  const defaultPositions = buildDefaultMultiRecordPositions(availableSelectors);
  const defaultRowBySelector = new Map(defaultPositions.map((entry) => [entry.selector, entry.row]));
  const normalizedCurrent = normalizeMultiRecordPositions(currentPositions, { maxRow: availableSelectors.length });
  /** @type {MultiRecordPosition[]} */
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
 * A file of the CLI reproducibility bundle: its content, or the builder that makes it when the bundle is written.
 * @typedef {RecipeGeneratedFile & { retainedBytes?: number }} GeneratedCliFile
 */

/**
 * @typedef {object} CanonicalRenderCandidateExecutionOptions
 * @property {Record<string, any>} canonical The candidate's request and resources.
 * @property {string} mode
 * @property {string} [kind]
 * @property {((progress: any) => void) | null} [onProgress]
 * @property {() => boolean} [shouldAdmit]
 * @property {((payload: Record<string, any>, options?: { onProgress?: ((progress: any) => void) | null, held?: ReturnType<typeof selectHeldRenderState> }) => Promise<any>) | null} [generationExecutor] Test seam.
 * @property {typeof admitFeatureCatalog} [catalogAdmission]
 * @property {typeof prepareCandidateRenderCommit} [prepareCommit]
 * @property {Record<string, any>} [prepareCommitInput]
 * @property {((canonical: Record<string, any>, catalogState: any) => any) | null} [decorationContinuity]
 * @property {any[]} [timingEntries]
 * @property {readonly string[] | null} [resultNames] The names of the Results a rerender draws again,
 *   which they keep (OV-136); Generate's Results take the names the engine gives them.
 * @property {FeatureCatalog | null} [heldFeatureCatalog] The catalog the page holds for the drawing it
 *   renders; when this render reuses its build, the reply says so and the page keeps it.
 * @property {readonly unknown[] | null} [heldResults] The drawing's committed Results; Python
 *   omits the children texts whose sanitized segments they keep.
 * @property {DrawingId | null} [drawingId] The drawing it renders: the Worker keeps the prepared
 *   inputs of each drawing's last request (A4b).
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
  resultNames = null,
  heldFeatureCatalog = null,
  heldResults = null,
  drawingId = null
}) => {
  if (!canonical?.renderRequest || !canonical?.resources) {
    throw new Error('Canonical candidate execution requires a request and resources.');
  }
  recordStructuralMetric('canonicalCandidateExecutionCount', 1, { kind });
  const startedAt = getNow();
  const held = selectHeldRenderState(heldFeatureCatalog, heldResults, resultNames);
  const generationResponse = generationExecutor
    ? await generationExecutor({
        request: canonical.renderRequest,
        resources: canonical.resources
      }, { onProgress, held })
    : await runDiagramGeneration({
        request: canonical.renderRequest,
        resources: canonical.resources
      }, { onProgress: onProgress, held, drawingId });
  const results = generationResponse.results;
  if (!shouldAdmit()) {
    return { status: /** @type {const} */ ('superseded'), generationResponse, elapsedMs: getNow() - startedAt };
  }
  if (results?.error) {
    return {
      status: /** @type {const} */ ('engine-error'),
      generationResponse,
      engineError: results.error,
      // Overlapping labels: Python's report, for the dialog (LP-D43).
      labelFit: generationResponse.labelFit ?? null,
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
  // Crossing label leaders drawn anyway: each Result's report (LP-D46).
  const labelFitReports = normalizeResultLabelFits(metadata.labelFitReports, results);
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
    labelFitReports,
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
 * A built generated artifact before it replaces the shown one: the owners' values and the runtime state to publish.
 * @typedef {{ readonly ownerSet: Readonly<Parameters<RunAnalysisOptions['installGeneratedArtifactOwnerSet']>[0]>,
 *   readonly runtimeState: { files?: unknown, archiveName?: unknown, losatTelemetry?: unknown } | null }} GeneratedArtifactCandidate
 */

/**
 * @typedef {object} RunAnalysisOptions
 * @property {WebState} state App state (state.js).
 * @property {RulePreparation | null} [rulePreparation] The rule owner's preparation, which Generate runs before it draws.
 * @property {(feature: RecordKeyedFeature) => boolean} isCurrentFeature Whether the feature belongs to the displayed Result.
 * @property {(comparisonPlanSnapshot: Record<string, any> | null, linearRecordCatalog: any, drawing: DrawingState) => any} serializeCanonicalFiles
 *   The Session owner's serialization of the active render files, chosen by the run's drawing.
 * @property {number} canonicalSessionVersion `SESSION_VERSION` of the Session owner.
 * @property {(canonical: Record<string, any>, options: { adoptOwnedRequest?: boolean, drawing: DrawingState }) => void} adoptCanonicalRenderArtifacts
 *   The Session owner's adoption of a candidate's request, resources, and files.
 * @property {(drawing: DrawingState, enabled: boolean) => void} setCircularConservationEnabled
 *   The ring owner's switch of a Circular drawing's comparison rings, with their managed stack rows.
 * @property {(() => (CommittedCanonicalSession | null)) | null} [getCommittedCanonicalSession] The committed canonical Session.
 * @property {() => number} [readCommittedCanonicalRevision] The Session owner's count of committed Session replacements.
 * @property {(writes: readonly LiveSettingWrite[]) => void} [commitLiveSettingValues]
 *   The Session owner's commit of the live setting values an automatic rerender drew (R1 (c)).
 * @property {(session: Record<string, any> | null | undefined, projectIdentity: (...args: any[]) => any, drawing: Record<string, any>) => any} [captureDecorationContinuity]
 *   The Legend layout owner's capture of the decoration positions a rerender of the drawing keeps.
 * @property {(drawing: DrawingState) => HistoryArtifactHandle | Promise<HistoryArtifactHandle>} captureGeneratedArtifactHandle History's capture of the current generated artifact, with the edits of `drawing`.
 * @property {() => DrawingState} displayedArtifactDrawing The drawing the displayed Result belongs to.
 * @property {(drawing: DrawingState) => ({ results: any[] } & Record<string, any>)} captureGeneratedArtifactOwnerSet
 *   The artifact owners' values, with the color rules and queued palette of `drawing`.
 * @property {(ownerSet: Record<string, any>, options: { drawing: DrawingState, selectedResultIndex: number, installResults: (results: any[]) => void }) => void} installGeneratedArtifactOwnerSet
 *   Installs an owner set as the displayed artifact, its color rules and queued palette into `drawing`.
 * @property {(drawingId: DrawingState['id'], redrawn: RedrawnArtifact, options: { installResults: (results: unknown[]) => void }) => void} installRedrawnArtifact
 *   The artifact owner's install of a label rerender into its Result's drawing (OV-111).
 * @property {(handle: HistoryArtifactHandle) => unknown} restoreGeneratedArtifactHandle
 * @property {((artifactIdentity: any, options: { results: any[] }) => void) | null} [setGeneratedArtifactIdentity]
 * @property {((label: string, execute: (before: HistoryArtifactHandle | null) => any, options?: HistoryStepOptions<Record<string, any>>) => Promise<Record<string, any>>) | null} [runGeneratedArtifactReplacement]
 *   History's undoable replacement of the generated artifact.
 * @property {RunAnalysisPreviewRuntime | null} [previewRuntime]
 * @property {() => Promise<void>} [nextTick] Vue `nextTick`.
 * @property {() => Promise<void>} [waitForAfterPaint]
 * @property {() => Promise<void>} [yieldToTask] The task yield after the later Generate statuses.
 * @property {() => Promise<void>} [waitForPostBindFrame]
 * @property {((capture: { phase?: string, diagnostics?: any }) => void) | null} [onGeneratedArtifactCheckpointCapture]
 * @property {(options?: { pan?: any, resetZoom?: boolean }) => void} resetPreviewViewport The preview owner's viewport reset.
 * @property {(() => ({ code: string, context?: any } | null | undefined)) | null} [validateAnnotationTargets]
 * @property {(() => string[]) | null} [omittedAnnotationRecordKeys] The record-catalog keys of the drawing's OFF records.
 * @property {((discovered: { primaryFile: InputFileRef | null, pairedFile: InputFileRef | null, recordCount: number }) => void) | null}
 *   [onCircularRecordsDiscovered] A discovery of the Circular input settled (success or error).
 * @property {(() => Promise<Record<string, any>>) | null} [prepareLinearRecordCatalog]
 * @property {{ value: Record<string, any>[] } | null} [recordDisplayRows]
 *   The draft record display rows the request reads (`recordDisplayControls.allRows`, R13).
 * @property {() => Promise<void>} [settlePendingInputs] The input reads Generate waits for before it reads the drawing:
 *   ring row names (D12) and auxiliary file imports still applying (OV-1193).
 * @property {((drawing: DrawingState) => void) | null} [assertDrawingInputs]
 * @property {(() => string) | null} [readDraftSignature] The History signature of the committed draft
 *   intent (UJ-10: a Cancel compares the canceled draft with the draft of the shown Result).
 *   services/config.js owns the active-mode input check shared with Save.
 * @property {() => number} [readHistoryStepSerial] History's `intentStepSerial`, read as a label rerender
 *   is requested: the Legend rows it draws join only the History steps begun after it (D-40).
 * @property {(write: () => void, request: { requestedAt: number }) => void} [adoptDrawnWrite] History's
 *   `adoptDrawnWrite`: a label rerender replaces the shown Results through it (OV-873).
 * @property {() => void} closeLabelTextScopeDialog The label owner's port (app/feature-editor/label-actions.js).
 * @property {(options?: { rerender?: boolean }) => void} clearLabelBuildNotices The label owner's port.
 * @property {typeof prepareCandidateRenderCommit} [prepareCandidateCommit] Test seam.
 * @property {AskChoice} [askChoice] The root's choice dialog: a large download asks first.
 * @property {ShowNotice} [showNotice] The root's notice banner.
 */

/** @param {RunAnalysisOptions} options */
export const createRunAnalysis = ({
  state,
  rulePreparation = null,
  isCurrentFeature,
  serializeCanonicalFiles,
  canonicalSessionVersion,
  adoptCanonicalRenderArtifacts,
  setCircularConservationEnabled,
  getCommittedCanonicalSession = null,
  readCommittedCanonicalRevision = () => 0,
  commitLiveSettingValues = () => {},
  captureDecorationContinuity = () => null,
  captureGeneratedArtifactHandle,
  displayedArtifactDrawing,
  captureGeneratedArtifactOwnerSet,
  installGeneratedArtifactOwnerSet,
  installRedrawnArtifact,
  restoreGeneratedArtifactHandle,
  setGeneratedArtifactIdentity = null,
  runGeneratedArtifactReplacement = null,
  previewRuntime = null,
  nextTick = () => Promise.resolve(),
  waitForAfterPaint = afterPaint,
  yieldToTask = afterTask,
  waitForPostBindFrame = afterFrame,
  onGeneratedArtifactCheckpointCapture = null,
  resetPreviewViewport,
  validateAnnotationTargets = null,
  omittedAnnotationRecordKeys = null,
  onCircularRecordsDiscovered = null,
  prepareLinearRecordCatalog = null,
  // The draft record display rows the request reads (`recordDisplayControls.allRows`, R13).
  recordDisplayRows = null,
  // D12, OV-1193: the input reads still applying settle before Generate reads the drawing.
  settlePendingInputs = async () => {},
  // services/config.js owns the active-mode input check shared with Save.
  assertDrawingInputs = null,
  readDraftSignature = null,
  readHistoryStepSerial = () => 0,
  adoptDrawnWrite = (write) => write(),
  // R13: the label owner's ports (app/feature-editor/label-actions.js).
  closeLabelTextScopeDialog,
  clearLabelBuildNotices,
  prepareCandidateCommit = prepareCandidateRenderCommit,
  askChoice = async () => 'cancel',
  showNotice = () => {}
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
    appliedPaletteName,
    appliedPaletteColors,
    mode,
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
    inputFiles,
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
    || typeof installRedrawnArtifact !== 'function'
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
  /**
   * Generate, a committed rerender and a label reflow render the shown
   * drawing: they name its feature catalog, so a reused build keeps it, and
   * the children texts its Results keep sanitized, so Python omits them.
   * @param {DrawingState} drawing The drawing the run renders, given by its caller.
   * @param {CanonicalRenderCandidateExecutionOptions} options
   */
  const executeShownDrawingCandidate = (drawing, options) => {
    const shown = captureGeneratedArtifactOwnerSet(drawing);
    return executeCanonicalRenderCandidate({
      ...options,
      heldFeatureCatalog: shown.featureCatalog ?? null,
      heldResults: shown.results ?? null,
      drawingId: drawing.id
    });
  };
  /**
   * @param {LosatJob[]} jobs
   * @param {LosatRunOptions} [options]
   * @returns {Promise<LosatJobResult[]>}
   */
  const executeLosatJobs = (jobs, options) => {
    const override = globalThis.__GBDRAW_LOSAT_EXECUTOR__;
    return (typeof override === 'function' ? override : runLosatPairsParallel)(jobs, options);
  };
  // Show Depth off leaves the Depth sources out of the request, so Python cannot name
  // the Depth rows a Legend style may still address; the draft names them (OV-81).
  /**
   * @param {DrawingState} drawing
   * @param {Pick<CanonicalRenderRequest, 'mode' | 'diagramOptions'>} request The candidate's request: the one
   *   `buildCanonicalRenderRequest` built, or the committed request that `projectCommittedEditorIntent` or the record
   *   rotation and alignment projections (`runCommittedCanonicalCandidate`) return. Only the first caller is checked:
   *   the other two pass a value typed `any` (the `canonical` option, and the return of `projectCommittedEditorIntent`),
   *   and for a loaded schema-9 request no reader checks `diagramOptions.depthTracks` yet (only its truth is read).
   */
  const unrequestedDepthCaptions = (drawing, request) => (
    request.diagramOptions?.depthTracks
      ? []
      : depthSeriesLegendCaptions({
          depthTracks: drawing.adv.depth_tracks,
          slots: request.mode === 'linear' ? drawing.adv.linear_track_slots : drawing.adv.circular_track_slots
        })
  );
  let pendingReflowRequestId = 0;
  let activeReflowRequestId = 0;
  // The latest request that makes the running iteration stale: an editor
  // request, or one for a live setting value that iteration applies. A
  // request for other live values queues the next iteration only: the running
  // one commits only its own values (OV-1075).
  let reflowSupersededAt = 0;
  /** @type {readonly string[]} The live setting request paths the running iteration applies. */
  let activeLiveSettings = [];
  // `readHistoryStepSerial()` at the latest request (D-40).
  let pendingReflowRequestedAt = 0;
  // A rerender requested while Generate (or a committed-candidate run) is in
  // progress waits for it: the rerender would take the newer generation token
  // and silently supersede the run (OV-48). It runs once after the run settles.
  // OV-1076: so do the requests of an action that holds the rerender while it
  // writes (`holdAutomaticRerender`): one action, one rerender.
  let reflowDeferred = false;
  let rerenderHolds = 0;
  const reflowWaits = () => processing.value || rerenderHolds > 0;
  // R1 (c), PD-OI-102: the request paths of live setting values each drawing's
  // actions committed and its Result does not show yet (R15: per drawing).
  // An automatic rerender applies them; it runs for a request of these alone
  // only while one of them still differs from the committed request. A
  // Generate keeps only the paths it did not draw, and a Load none (6.7).
  /** @type {Map<DrawingState, Set<string>>} */
  const pendingLiveSettings = new Map();
  /** @param {DrawingState} drawing */
  const pendingLiveSettingsOf = (drawing) => {
    const pending = pendingLiveSettings.get(drawing) || new Set();
    pendingLiveSettings.set(drawing, pending);
    return pending;
  };
  // Whether a request since the last iteration began asked for more than the
  // live setting values (an editor edit's rerender).
  let editorRerenderRequested = false;
  // The label rerender loop while it runs. A run that sets `processing` waits
  // for it before it renders: its iteration may already be in the diagram
  // Worker, which takes one request at a time (OV-349), and no further
  // iteration starts.
  /** @type {Promise<void> | null} */
  let activeReflowRun = null;
  let featureExtractionRequestId = 0;
  let latestGenerationToken = 0;
  let latestOperationId = 0;
  // Each drawing's Circular discovery: its read generation and the read still running.
  /** @type {Map<DrawingId, number>} */
  const circularRecordRefreshGenerations = new Map();
  /** @type {Map<DrawingId, { fingerprint: unknown[], promise: Promise<any> | null }>} */
  const activeCircularRecordRefreshes = new Map();
  /** @type {AbortController | null} */
  let activeLosatAbortController = null;
  /** @type {readonly Record<string, any>[]} */
  let latestCliHelperFiles = Object.freeze([]);
  let latestCliHelperArchiveName = 'out-cli-files.zip';
  let latestCliHelperRetainedBytes = 0;
  /** @param {unknown} files */
  const cloneCliHelperFiles = (files) => (
    (Array.isArray(files) ? files : []).map((file) => ({ ...file }))
  );
  /** @param {{ files?: unknown, archiveName?: unknown, losatTelemetry?: unknown }} runtimeState */
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
     * @param {{ runtimeState?: GeneratedArtifactCandidate['runtimeState'] }} [options]
     */
    build(ownerSet, { runtimeState = null } = {}) {
      recordSessionLifecycleEvent('artifact.candidate-completed');
      recordStructuralMetric('generatedArtifactCandidateBuildCount', 1);
      return Object.freeze({
        ownerSet: Object.freeze(ownerSet),
        runtimeState
      });
    },
    /**
     * @param {GeneratedArtifactCandidate} candidate
     * @param {{ drawing: DrawingState, selectedResultIndex?: number }} options `drawing`: the run's drawing.
     */
    activate(candidate, { drawing, selectedResultIndex = 0 }) {
      recordSessionLifecycleEvent('artifact.activation-started', { selectedResultIndex });
      installGeneratedArtifactOwnerSet(candidate.ownerSet, {
        drawing,
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
    /** @param {HistoryArtifactHandle} handle */
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

  /**
   * @param {DrawingState} drawing
   * @param {{ helperFiles?: readonly RecipeHelperFile[] } | null | undefined} runInfo
   * @param {Map<string, GeneratedCliFile>} generatedCliFileMap
   * @param {unknown} archiveBaseName
   */
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

  const downloadCliHelperFiles = async () => {
    if (!latestCliHelperFiles.length) {
      showNotice({ text: t('The latest run has no reproducibility files to download.') });
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
    if (totalChars > 50 * 1024 * 1024 && !await askBeforeLargeDownload(askChoice, {
      title: t('Download large reproducibility files?'),
      description: t('The ZIP holds about {size} MB of input and helper files.', { size: (totalChars / (1024 * 1024)).toFixed(1) })
    })) return;
    downloadBlob(createZipBlob(materializedFiles), latestCliHelperArchiveName);
  };

  /**
   * @param {unknown} cause A caught value or a diagnostic.
   * @param {string} [operation]
   * @param {string} [stage]
   */
  const formatError = (cause, operation = 'generate', stage = 'request-validation') =>
    // The argument is truthy (`cause || {...}`), so the model is not null.
    /** @type {UserFacingError} */ (normalizeUserFacingError(cause || { code: 'UNKNOWN' }, { operation, stage }));
  // UJ-08: the shown Generate failure about the source inputs, and its mode.
  /** @type {{ error: any, mode: string } | null} */
  let sourceInputFailure = null;
  /** @type {ReadonlySet<string>} */
  const SOURCE_INPUT_FAILURES = new Set(['INPUT_REQUIRED', 'FASTA_REQUIRED', 'INPUT_UNREADABLE', 'NO_RECORDS']);
  /**
   * Record discovery read every source input of `drawing`: a shown Generate
   * failure about those inputs no longer holds, so it ends here, whether the
   * inputs now give records or the inspection line says why they do not
   * (OV-1257).
   * @param {DrawingState} drawing
   */
  const releaseSourceInputFailure = (drawing) => {
    if (!sourceInputFailure || sourceInputFailure.mode !== drawing.mode
      || errorLog.value !== sourceInputFailure.error) return;
    try {
      assertDrawingInputs?.(drawing);
    } catch (_) {
      return;
    }
    sourceInputFailure = null;
    errorLog.value = null;
  };
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
        const currentResults = captureGeneratedArtifactOwnerSet(displayedArtifactDrawing()).results;
        const unchanged = currentResults.length === previousResults.length
          && currentResults.every((result, index) => result === previousResults[index]);
        recovery = !handle ? (results.value.length ? 'restore-failed' : 'no-result')
          : previousResults.length === 0 ? 'no-result'
          : restored ? 'restored' : unchanged ? 'preserved' : 'restore-failed';
      } catch (_) { recovery = 'restore-failed'; }
    }
    if (!isCurrent()) return { status: 'stale' };
    errorLog.value = error;
    sourceInputFailure = operation === 'generate' && SOURCE_INPUT_FAILURES.has(error?.code)
      ? { error: errorLog.value, mode: mode.value } : null;
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
   * @param {Map<string, RawLosatCacheEntry> | null | undefined} cacheMap
   * @param {string} cacheKey
   * @param {LosatCacheMetadata} metadata
   * @param {Record<string, any> | null} [manifest]
   */
  const getReusableLosatCacheEntry = (cacheMap, cacheKey, metadata, manifest = null) => {
    if (completedLosatSearch?.owner !== losatCache.value) completedLosatSearch = null;
    return getRawLosatCacheEntry(cacheMap, cacheKey, metadata, manifest)
      || getRawLosatCacheEntry(completedLosatSearch?.entries, cacheKey, metadata, manifest);
  };
  /**
   * @param {Map<string, RawLosatCacheEntry>} cacheMap
   * @param {{ cacheKey: string }[]} pairs
   */
  const retainCompletedLosatSearch = (cacheMap, pairs) => {
    completedLosatSearch = {
      owner: losatCache.value,
      entries: new Map(pairs.map(({ cacheKey }) => [cacheKey, cacheMap.get(cacheKey)]))
    };
  };

  /** @param {AbortSignal | null | undefined} signal */
  const getGenerationCancelReason = (signal) =>
    signal?.reason instanceof Error ? signal.reason : new DiagramGenerationCanceledError();

  /**
   * @template T
   * @param {Promise<T> | null} promise
   * @param {AbortSignal | null | undefined} signal
   * @returns {Promise<T | null> | null}
   */
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
  /**
   * @param {DrawingState} drawing
   * @param {readonly OrthogroupGroup[] | null | undefined} previousGroups
   * @param {readonly OrthogroupGroup[] | null | undefined} candidateGroups
   */
  const rekeyCommittedOrthogroupOverrides = (drawing, previousGroups, candidateGroups) => {
    const next = rekeyOrthogroupOverrides({
      previousGroups,
      candidateGroups,
      names: drawing.orthogroupNameOverrides,
      descriptions: drawing.orthogroupDescriptionOverrides,
      dormant: drawing.orthogroupDormantOverrides
    });
    GROUP_NAMES_SLICE.install(drawing.slices.groupNames, {
      orthogroupNameOverrides: next.names,
      orthogroupDescriptionOverrides: next.descriptions,
      orthogroupDormantOverrides: next.dormant
    });
  };

  /** @param {{ status?: string, generationId?: number, error?: unknown, summaryCount?: number, detailsCacheSize?: number }} [updates] */
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

  /**
   * @param {LinearSeq | undefined} seq
   * @param {string} [fallback]
   */
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

  /**
   * @param {DrawingState} drawing
   * @param {string} leftLabel
   * @param {string} rightLabel
   */
  const buildLosatFilename = (drawing, leftLabel, rightLabel) => (
    losatEdgeFilename(leftLabel, rightLabel, buildLosatSuffix(drawing))
  );

  /**
   * @param {DrawingState} drawing
   * @param {unknown} edgeKey
   */
  const getResolvedLinearEdge = (drawing, edgeKey) => {
    const normalizedKey = String(edgeKey || '').trim();
    if (!normalizedKey) return null;
    const edges = drawing.linearComparisonResolution?.value?.edges;
    return (Array.isArray(edges) ? edges : []).find((edge) => edge.edgeKey === normalizedKey) || null;
  };

  /** @param {string | number} edgeKey An edge key, or the entry's index. */
  const getLosatCacheInfoEntry = (edgeKey) => {
    /** @type {LosatCacheInfoEntry[]} */
    const entries = Array.isArray(losatCacheInfo.value) ? losatCacheInfo.value : [];
    // `Number.isInteger` is true only for a number.
    if (Number.isInteger(edgeKey)) return entries[/** @type {number} */ (edgeKey)] || null;
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
    const drawnSeqs = drawing.drawnLinearSeqs.value;
    const leftLabel = getSeqLabel(
      drawnSeqs[queryIndex],
      queryEntry?.recordId || `seq_${Number.isInteger(queryIndex) ? queryIndex + 1 : 1}`
    );
    const rightLabel = getSeqLabel(
      drawnSeqs[subjectIndex],
      subjectEntry?.recordId || `seq_${Number.isInteger(subjectIndex) ? subjectIndex + 1 : 2}`
    );
    return buildLosatFilename(drawing, leftLabel, rightLabel);
  };
  // LOSAT pairs are the Linear comparisons: the three pair handlers act for
  // the drawing their caller (the panel binding) gives, and do nothing for a
  // Circular one (D-15-89).
  /**
   * @param {DrawingState} drawing
   * @param {string} edgeKey
   * @param {{ recordId?: string } | null} [queryEntry]
   * @param {{ recordId?: string } | null} [subjectEntry]
   */
  const getLosatPairDefaultName = (drawing, edgeKey, queryEntry = null, subjectEntry = null) => {
    if (drawing.mode !== 'linear') return '';
    return losatPairDefaultName(drawing, edgeKey, queryEntry, subjectEntry);
  };

  /**
   * @param {unknown} name
   * @param {unknown} fallback
   */
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


  /**
   * @param {unknown} cacheKey
   * @param {RawLosatCacheEntry} cached
   * @param {unknown} [identityManifest] A copy of the current manifest; Download all copies it once for every table.
   */
  const hydrateLosatDownloadText = async (cacheKey, cached, identityManifest = cloneJsonData(proteinIdentityManifest.value)) => {
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
        identityManifest
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

  /**
   * @param {DrawingState} drawing
   * @param {string} edgeKey
   * @param {unknown} customName
   */
  const downloadLosatPair = async (drawing, edgeKey, customName) => {
    if (drawing.mode !== 'linear') return;
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
      losatCacheInfo.value = losatCacheInfo.value.map((/** @type {LosatCacheInfoEntry} */ candidate) => (
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

  /**
   * @param {DrawingState} drawing
   * @param {string} edgeKey
   * @param {unknown} customName
   */
  const setLosatPairFilename = (drawing, edgeKey, customName) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (drawing.mode !== 'linear') return;
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
    losatCacheInfo.value = losatCacheInfo.value.map((/** @type {LosatCacheInfoEntry} */ candidate) => (
      candidate === entry ? { ...candidate, filename } : candidate
    ));
  };

  /** @param {DrawingState} drawing */
  const circularDiscoveryTargetsCurrentInput = (drawing) => circularDiscoveryForInput(drawing).current;

  /** @param {DrawingState} drawing */
  const circularDiscoveryMatchesCurrentInput = (drawing) => (
    drawing.circularRecordDiscovery.status === 'ready' &&
    circularDiscoveryTargetsCurrentInput(drawing)
  );

  /** @param {DrawingState} drawing */
  const validateDepthInputPresence = (drawing) => {
    const linear = mode.value === 'linear';
    const slots = linear ? drawing.adv.linear_track_slots : drawing.adv.circular_track_slots;
    if (!draftDrawsDepth(linear ? 'linear' : 'circular', drawing)) return '';
    // The draft draws Depth; with the custom stack on, through an enabled Depth slot.
    const customStack = Boolean(linear ? drawing.adv.linear_track_slots_enabled : drawing.adv.circular_track_slots_enabled);
    let rows;
    if (linear) {
      rows = drawing.drawnLinearSeqs.value.map((seq) => depthFileSlotsFromValue(seq.depth));
    } else {
      const discoveredCount = (
        circularDiscoveryMatchesCurrentInput(drawing) &&
        Array.isArray(drawing.circularRecordList.value)
      )
        ? drawing.circularRecordList.value.length
        : 0;
      const recordCount = discoveredCount > 0
        ? discoveredCount
        : (
            isRecordMajorDepthFileMatrix(drawing.files.c_depth)
              ? Math.max(1, drawing.files.c_depth.length)
              : 1
          );
      if (isRecordMajorDepthFileMatrix(drawing.files.c_depth) && drawing.files.c_depth.length !== recordCount) {
        return `Circular Depth matrix has ${drawing.files.c_depth.length} record rows; expected ${recordCount}.`;
      }
      rows = normalizeRecordMajorDepthFileRows(drawing.files.c_depth, recordCount);
    }
    if (customStack) {
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

  /**
   * Reads the Circular records of `drawing`'s inputs into its Multi-Record
   * Canvas order. It writes only `drawing`, shown or not: a read of another
   * drawing neither supersedes it nor lands in it.
   * @param {DrawingState} drawing
   * @param {{ suppress?: boolean, automatic?: boolean, reuseFinalError?: boolean }} [options]
   */
  const runCircularRecordRefresh = async (drawing, { suppress = false, automatic = false, reuseFinalError = false } = {}) => {
    if (drawing.mode !== 'circular') return;
    const currentGeneration = () => circularRecordRefreshGenerations.get(drawing.id) ?? 0;
    if (automatic && circularDiscoveryTargetsCurrentInput(drawing)
      && ['deferred', 'ready', 'error'].includes(drawing.circularRecordDiscovery.status)) return;
    // Generate reports a final error of these exact files without reading them
    // again; Retry source inspection still reads.
    if (reuseFinalError && circularDiscoveryTargetsCurrentInput(drawing)
      && drawing.circularRecordDiscovery.status === 'error' && discoveryErrorIsFinal(drawing.circularRecordDiscovery.error)) return;
    const refreshGeneration = currentGeneration() + 1;
    circularRecordRefreshGenerations.set(drawing.id, refreshGeneration);
    if (suppress || recordDiscoverySuppressed()) return;
    if (!Array.isArray(drawing.adv.multi_record_positions)) {
      drawing.adv.multi_record_positions = [];
    }
    const inputType = drawing.inputType.value;
    const primaryFile = inputType === 'gff' ? drawing.files.c_gff : drawing.files.c_gb;
    const pairedFile = inputType === 'gff' ? drawing.files.c_fasta : null;
    const hasCompleteInput = Boolean(primaryFile && (inputType !== 'gff' || pairedFile));
    const preserveCanonicalRecordKeys = circularDiscoveryTargetsCurrentInput(drawing);
    const preservedRecordKeys = new Map();
    if (preserveCanonicalRecordKeys) {
      [
        ...(Array.isArray(drawing.circularRecordDiscovery.canonicalRecordIdentities)
          ? drawing.circularRecordDiscovery.canonicalRecordIdentities
          : []),
        ...drawing.circularRecordList.value
      ]
        .forEach((/** @type {Readonly<Record<string, unknown>>} */ record) => {
          const recordKey = String(record?.recordKey || '').trim();
          if (!recordKey) return;
          [record?.selector, record?.record_id, record?.recordId]
            .map((value) => String(value || '').trim())
            .filter(Boolean)
            .forEach((value) => preservedRecordKeys.set(value, recordKey));
        });
    }
    Object.assign(drawing.circularRecordDiscovery, {
      status: hasCompleteInput ? 'loading' : 'idle',
      error: '',
      inputType,
      primaryFile: primaryFile || null,
      pairedFile: pairedFile || null,
      canonicalRecordIdentities: preserveCanonicalRecordKeys
        ? drawing.circularRecordDiscovery.canonicalRecordIdentities
        : []
    });
    drawing.circularRecordList.value = [];
    if (!hasCompleteInput) {
      drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length);
      return;
    }
    if (
      refreshGeneration !== currentGeneration() ||
      recordDiscoverySuppressed() ||
      drawing.inputType.value !== inputType ||
      (inputType === 'gff' ? drawing.files.c_gff : drawing.files.c_gb) !== primaryFile ||
      (inputType === 'gff' ? drawing.files.c_fasta : null) !== pairedFile
    ) return;

    try {
      const records = inputType === 'gff'
        ? await discoverGffFastaRecords({ gffFile: inputFiles.fileOf(primaryFile), fastaFile: inputFiles.fileOf(pairedFile) })
        : await discoverSequenceRecords({ file: inputFiles.fileOf(primaryFile), format: 'genbank' });
      if (
        refreshGeneration !== currentGeneration() ||
        recordDiscoverySuppressed() ||
        drawing.inputType.value !== inputType ||
        (inputType === 'gff' ? drawing.files.c_gff : drawing.files.c_gb) !== primaryFile ||
        (inputType === 'gff' ? drawing.files.c_fasta : null) !== pairedFile
      ) return;
      const busy = state.sessionOperationAvailability?.();
      if (busy) {
        drawing.circularRecordDiscovery.status = 'deferred';
        drawing.circularRecordDiscovery.primaryFile = null;
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
      drawing.circularRecordList.value = nextRecords;
      drawing.circularRecordDiscovery.status = 'ready';
      // Identity belongs to the source, including while its mode is inactive.
      drawing.circularRecordDiscovery.canonicalRecordIdentities = nextRecords
        .filter((record) => record.recordKey)
        .map(({ selector, record_id, recordKey }) => ({ selector, record_id, recordKey }));
      const nextPositions = mergeCircularRecordPositions(nextRecords, drawing.adv.multi_record_positions);
      drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length, ...nextPositions);
      if (nextRecords.length > 0) releaseSourceInputFailure(drawing);
      onCircularRecordsDiscovered?.({ primaryFile, pairedFile, recordCount: nextRecords.length });
    } catch (error) {
      if (
        refreshGeneration !== currentGeneration() ||
        recordDiscoverySuppressed() ||
        drawing.inputType.value !== inputType ||
        (inputType === 'gff' ? drawing.files.c_gff : drawing.files.c_gb) !== primaryFile ||
        (inputType === 'gff' ? drawing.files.c_fasta : null) !== pairedFile
      ) return;
      const busy = state.sessionOperationAvailability?.();
      if (busy) {
        drawing.circularRecordDiscovery.status = 'deferred';
        drawing.circularRecordDiscovery.primaryFile = null;
        return busy;
      }

      drawing.circularRecordList.value = [];
      drawing.circularRecordDiscovery.status = 'error';
      drawing.circularRecordDiscovery.error = formatError(error, inputType === 'gff' ? 'listGffFastaRecords' : 'listSequenceRecords', 'helper');
      drawing.adv.multi_record_positions.splice(0, drawing.adv.multi_record_positions.length);
      releaseSourceInputFailure(drawing);
      onCircularRecordsDiscovered?.({ primaryFile, pairedFile, recordCount: 0 });
    }
  };

  /**
   * The owner of the Circular record refresh: it reads the drawing its caller
   * gives (Generate: the drawing it pinned; a transition that changed the
   * drawing's inputs: that drawing, R10) and shares a read still running for
   * the same drawing and inputs. The record editor gives the shown drawing
   * (app-setup's binding, D-15-89).
   * @param {{ suppress?: boolean, automatic?: boolean, reuseFinalError?: boolean } | undefined} options
   * @param {DrawingState} drawing
   */
  const refreshCircularRecordOrder = (options, drawing) => {
    options = options || {};
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const inputType = drawing.inputType.value;
    const fingerprint = [
      Boolean(options.suppress || recordDiscoverySuppressed()),
      drawing.id,
      inputType,
      inputType === 'gff' ? drawing.files.c_gff : drawing.files.c_gb,
      inputType === 'gff' ? drawing.files.c_fasta : null
    ];
    const inflightRefresh = activeCircularRecordRefreshes.get(drawing.id);
    if (
      inflightRefresh &&
      fingerprint.length === inflightRefresh.fingerprint.length &&
      fingerprint.every((value, index) => (
        Object.is(value, inflightRefresh.fingerprint[index])
      ))
    ) {
      // An entry is activated after its promise is assigned.
      return /** @type {NonNullable<typeof inflightRefresh.promise>} */ (inflightRefresh.promise);
    }
    const entry = { fingerprint, promise: /** @type {Promise<any> | null} */ (null) };
    entry.promise = runCircularRecordRefresh(drawing, options).finally(() => {
      if (activeCircularRecordRefreshes.get(drawing.id) === entry) activeCircularRecordRefreshes.delete(drawing.id);
    });
    activeCircularRecordRefreshes.set(drawing.id, entry);
    return entry.promise;
  };

  // R13: Generate and the label reflow prepare the specific-color rules they
  // render here, the Generate compiler's one trigger of the rule preparation.
  // `isCurrent` is the run's own staleness: the candidate is null when the run
  // went stale during the preparation, and is admitted (`shouldAdmit`) while
  // the run and the rule inputs it was prepared from stay current.
  /**
   * @param {DrawingState} drawing
   * @param {() => boolean} isCurrent
   * @param {Record<string, any>} [options]
   */
  const prepareAndAdmitCandidate = async (drawing, isCurrent, options) => {
    if (!rulePreparation) return { rules: drawing.manualSpecificRules, shouldAdmit: isCurrent };
    const candidate = await rulePreparation.prepareCandidate(drawing.manualSpecificRules, options);
    if (!candidate || !isCurrent()) return null;
    // Without `asGiven`, the preparation wrote each rule as `normalizeSpecificRule` does.
    return { rules: /** @type {readonly SpecificColorRule[]} */ (candidate.rules), shouldAdmit: () => rulePreparation.isCurrent(candidate.snapshot) && isCurrent() };
  };

  // A render of the committed Session draws the drawing's live editor paint (PD-OI-102): the
  // label rerender and a committed candidate (a record rotation, an alignment Apply, OV-1083)
  // render the committed Session with the editor intent and the rules prepared for it, so their
  // Result shows the live edits as a Generate of the same settings would. `drawnPaint` is the
  // editor state it draws from; an edit made later is shown once its Result is bound (OV-346).
  // Null when the run went stale during the rule preparation.
  /**
   * @param {DrawingState} drawing
   * @param {CommittedCanonicalSession & { config?: Record<string, unknown> | null }} committed The committed Session.
   * @param {() => boolean} isCurrent
   * @param {LiveSettingWrite[]} [liveSettingWrites] Draft values of live setting rows the render applies.
   */
  const prepareCommittedPaintCandidate = async (drawing, committed, isCurrent, liveSettingWrites = []) => {
    const drawnPaint = editorPaintState(/** @type {import('./result-paint-record.js').PaintStateSource} */ (state), drawing);
    const colorCandidate = await prepareAndAdmitCandidate(drawing, isCurrent);
    if (!colorCandidate) return null;
    const canonical = projectCommittedEditorIntent({
      committed,
      state,
      drawing: {
        ...drawing,
        selectedPalette: appliedPaletteName,
        currentColors: appliedPaletteColors,
        manualSpecificRules: colorCandidate.rules
      },
      liveSettingWrites
    });
    return { canonical, drawnPaint, rules: colorCandidate.rules, shouldAdmit: colorCandidate.shouldAdmit };
  };

  /**
   * @param {DrawingState} drawing
   * @param {{ decorationContinuity?: any, comparisonPlanSnapshot?: Record<string, any> | null,
   *   generatedArtifactHandle?: HistoryArtifactHandle | null, comparisonExecution?: Record<string, any> | null,
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
    // The records this Generate draws; every Linear index below is a position in it.
    const linearSeqs = drawing.drawnLinearSeqs.value;
    const useCommittedComparison = comparisonExecution?.mode === 'inherit';
    const forceEmptyComparison = useCommittedComparison || comparisonExecution?.mode === 'clear';
    const activeComparisonPlanSnapshot = mode.value === 'linear'
      ? (
          forceEmptyComparison
            ? resolveLinearComparisonPlan({
                plan: { mode: 'none', defaultSource: 'losat', edges: [] },
                sequences: drawing.linearSeqs,
                recordsOff: drawing.recordsOff,
                layout: drawing.linearRecordLayoutEnabled.value ? drawing.linearRecordRows : [],
                losatProgram: drawing.losatProgram.value,
                blastpMode: normalizeBlastpMode(drawing.losat.blastp?.mode)
              })
            : comparisonPlanSnapshot || resolveLinearComparisonPlan({
                plan: drawing.linearComparisonPlan,
                sequences: drawing.linearSeqs,
                recordsOff: drawing.recordsOff,
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
    const committedArtifactHandle = generatedArtifactHandle || await captureGeneratedArtifactHandle(drawing);
    /** @type {ReturnType<typeof generatedArtifactTransactionOwner.build> | null} */
    let activatedGeneratedArtifactCandidate = null;
    /** @type {ReadyReceipt | null} */
    let acceptedCandidateReadyReceipt = null;
    let workingProteinIdentityManifest = committedArtifactHandle?.ownerSet
      ?.proteinIdentityManifest ?? proteinIdentityManifest.value;
    /** @type {LegacyProteinCandidateEnvelope | null} */
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
      || drawing.legendEntries.value.some((/** @type {LegendRow} */ entry) => entry.originalCaption && /** @type {string} */ (entry.originalCaption) !== entry.caption)
      || drawing.dormantLegendEntries.value.length > 0;
    // The committed records whose source this Generate replaces or drops. The
    // metadata of a Session without a feature catalog has no record keys.
    const replacedFeatures = hasSourceBoundEditorIntent
      ? [...new Map((workingBiologicalFeatures || []).filter((/** @type {FeatureCatalogRow} */ feature) => feature.record_key)
        .map((/** @type {FeatureCatalogRow} */ feature) => [feature.record_key, feature])).values()]
        .filter((/** @type {FeatureCatalogRow} */ feature) => !isCurrentFeature(feature))
      : [];
    const sourceReplaced = replacedFeatures.length > 0;
    const previousCommittedRequest = typeof getCommittedCanonicalSession === 'function'
      ? getCommittedCanonicalSession()?.renderRequest : null;
    const previousRequestRecords = previousCommittedRequest && previousCommittedRequest.mode === mode.value
      ? previousCommittedRequest.records || [] : [];
    let workingLosatCacheInfo = committedArtifactHandle?.ownerSet?.losatCacheInfo
      ?? losatCacheInfo.value;
    let workingSelectedOrthogroupId = selectedOrthogroupId.value;
    let workingSelectedOrthogroupAlignmentFeature = selectedOrthogroupAlignmentFeature.value;
    const restoreCommittedArtifact = async () => {
      if (!committedArtifactHandle || !activatedGeneratedArtifactCandidate) return false;
      if (captureGeneratedArtifactOwnerSet(drawing).results !== activatedGeneratedArtifactCandidate.ownerSet.results) return false;
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
      processingStatus.value = t('Canceled.');
      generationCancelRequested.value = false;
      return { status: 'canceled' };
    };
    /** @param {unknown} message */
    const setProcessingStatus = (message) => {
      if (generationToken === latestGenerationToken) {
        processingStatus.value = String(message || '');
      }
    };
    /** @param {Parameters<typeof diagramProgressMessage>[0]} progress */
    const onDiagramProgress = (progress) => {
      const message = diagramProgressMessage(progress);
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
      processingStatus.value = t('Reading input records...');
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
          processingStatus.value = t('Preparing input files...');
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
      const inputType = drawing.inputType.value;
      const primaryFile = inputType === 'gff' ? drawing.files.c_gff : drawing.files.c_gb;
      const pairedFile = inputType === 'gff' ? drawing.files.c_fasta : null;
      const hasCompleteInput = Boolean(primaryFile && (inputType !== 'gff' || pairedFile));
      if (
        hasCompleteInput &&
        !circularDiscoveryMatchesCurrentInput(drawing)
      ) {
        processingStatus.value = t('Reading input records...');
        try {
          await refreshCircularRecordOrder({ reuseFinalError: true }, drawing);
        } finally {
          if (generationToken === latestGenerationToken) {
            processingStatus.value = t('Preparing input files...');
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
        if (!circularDiscoveryMatchesCurrentInput(drawing)) {
          const message = drawing.circularRecordDiscovery.error || diagnosticError('INPUT_UNREADABLE');
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
    /** @type {{ cacheMap: Map<string, RawLosatCacheEntry>, cacheKey: string, candidateIndex: number, proteinIdMap: Record<string, unknown> }[]} */
    const legacyPromotionTransaction = [];
    let legacyPromotionCommitted = false;
    /** @type {((candidateOwnerSet: Record<string, any>) => { ownerSet: Record<string, any>, selectedOrthogroupAlignmentFeature: any }) | null} */
    let commitProteinMigration = null;
    /** @type {{ cacheInfo: any[], cacheMap: Map<any, any>, derivedCacheMap: Map<any, any> | null } | null} */
    let pendingLosatCacheCommit = null;

    featureExtractionRequestId += 1;
    processingStatus.value = t('Preparing input files...');
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
        requireLinearComparisonPlanSnapshot(activeComparisonPlanSnapshot);
      }
      if (typeof validateAnnotationTargets === 'function' && drawing.annotationSets.length > 0) {
        const annotationError = validateAnnotationTargets();
        if (annotationError) throw diagnosticError(annotationError.code, annotationError.context);
      }
      /** @type {({ file: string, displayFile: string } | null)[]} */
      let regionSpecs = [];
      /** @type {string[]} */
      let recordSelectors = [];
      const resolvedComparisons = [];
      /** @type {ResolvedCircularConservation[]} */
      let resolvedCircularConservation = [];
      /**
       * Each ring LOSAT search, in ring order: its raw cache key (null without current search output) and the
       * Source recipe's generated BLAST table.
       * @type {{ cacheKey: string | null, path: string, name: string, slot: string }[]}
       */
      const circularRingSearches = [];
      /** @type {Map<string, GeneratedCliFile>} */
      const generatedCliFileMap = new Map();
      const textEncoder = new TextEncoder();
      const textDecoder = new TextDecoder();
      /**
       * @param {unknown} path
       * @param {string} [fallback]
       */
      const getPayloadName = (path, fallback = 'input') => {
        const name = String(path || '').split('/').filter(Boolean).pop();
        return name || fallback;
      };
      /**
       * @param {unknown} path
       * @param {unknown} data
       * @param {{ name?: unknown, slot?: unknown }} [options]
       */
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
      /**
       * @param {unknown} path
       * @param {NonNullable<GeneratedCliFile['buildData']>} buildData
       * @param {{ name?: unknown, slot?: unknown, retainedBytes?: unknown }} [options]
       */
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
      const normalizedOutputPrefix = String(drawing.form.prefix || '').trim();

      // The annotation, color, whitelist, priority, label and feature visibility tables reach Run Info from the
      // request (the Source recipe writes its own copies), so Generate stages no copy of them (OV-1573, OV-1589).
      // The request's label table holds label rules only; per-feature label
      // edits travel as featureOverrides rows (design Q4). One projection per
      // Generate serves the request (CW-02).
      generatedLabelProjection = requestLabelProjection(state, requestDrawing);
      // Per-feature Feature visibility decides LOSATP's proteins as the rules do;
      // the helper resolves one record's rows as Generate does (R4).
      /** @param {string[]} recordKeys */
      const proteinVisibilityRows = (recordKeys) => requestFeatureOverrides(
        drawing.featureOverrides, recordKeys.map((recordKey) => ({ recordKey }))
      ).filter((row) => row.featureVisibility !== null)
        .map((row) => ({ ...row, labelVisibility: null, labelText: null }));
      const featureVisibilityTsv = serializeFeatureVisibilityRules(drawing.featureVisibilityRules?.value || []);
      // The rules table the protein helper reads, and part of its cache key; '' when there are no rules.
      const featureVisibilityCacheKey = featureVisibilityTsv.trim() ? featureVisibilityTsv : '';
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
      /** @param {unknown} value */
      const normalizeDepthTrackValue = (value) => {
        if (value === null || value === undefined || value === '') return null;
        const numeric = Number(value);
        return Number.isFinite(numeric) && numeric > 0 ? String(numeric) : null;
      };
      /** @param {(InputFileRef | null)[][]} rows */
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
      // Generate gives each series the drawing has its settings through the
      // series owner, before a row reads them; it adds no series (OV-1510).
      /** @param {number} seriesCount */
      const shapeDepthSeriesForGenerate = (seriesCount) => {
        if (!Array.isArray(drawing.adv.depth_tracks)) drawing.adv.depth_tracks = [];
        shapeDepthTracksForGenerate(drawing.adv.depth_tracks, seriesCount, drawing.adv);
      };
      /**
       * @param {DepthSlotRow[]} slots
       * @param {number} seriesCount
       */
      const syncDepthSlotLegendLabelsFromTrackConfigs = (slots, seriesCount) => {
        syncDepthSlotLabels({ slots, depthTracks: drawing.adv.depth_tracks, activeCount: seriesCount });
      };
      // A Depth row draws its series' stored height, read without writing: a
      // series without settings, or one the drawing does not have, leaves the
      // row its own height (TK-03).
      /** @param {LinearTrackSlot | null | undefined} slot */
      const syncLinearDepthSlotHeightFromTrackConfig = (slot) => {
        if (!slot || String(slot.renderer || '') !== 'depth') return;
        const rawTrackIndex = Number(slot.params?.track_index);
        const trackIndex = Number.isInteger(rawTrackIndex) && rawTrackIndex >= 0 ? rawTrackIndex : 0;
        /** @type {unknown} */
        const config = Array.isArray(drawing.adv.depth_tracks) ? drawing.adv.depth_tracks[trackIndex] : null;
        if (!config || typeof config !== 'object') return;
        const configuredHeight = normalizeDepthTrackValue(/** @type {RawDepthTrackConfig} */ (config).height);
        if (configuredHeight) slot.height = configuredHeight;
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
        // The rows' sides follow one Axis index, and Generate stores that index (OV-1263).
        const placedCircularTrackStack = useCircularTrackSlots
          ? placeCircularTrackStack(
              drawing.adv.circular_track_slots,
              drawing.adv.nt,
              drawing.form.track_type,
              circularTrackAxisIndex
            )
          : null;
        const circularTrackSlots = placedCircularTrackStack
          ? applyCircularSuppressControlsToSlots(placedCircularTrackStack.slots, drawing.form, drawing.adv.nt)
          : [];
        if (placedCircularTrackStack) {
          drawing.adv.circular_track_slots.splice(0, drawing.adv.circular_track_slots.length, ...circularTrackSlots);
          drawing.adv.circular_track_slots_axis_index = placedCircularTrackStack.axisIndex;
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
            drawing.circularRecordList.value,
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
          circularDiscoveryMatchesCurrentInput(drawing) &&
          Array.isArray(drawing.circularRecordList.value)
        )
          ? drawing.circularRecordList.value.length
          : 0;
        const circularDepthRecordCount = discoveredCircularRecordCount > 0
          ? discoveredCircularRecordCount
          : (
              isRecordMajorDepthFileMatrix(drawing.files.c_depth)
                ? Math.max(1, drawing.files.c_depth.length)
                : 1
            );
        const circularDepthRows = normalizeRecordMajorDepthFileRows(
          drawing.files.c_depth,
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
          const circularDepthSeriesCount = representativeDepthFiles(circularDepthRows).length;
          shapeDepthSeriesForGenerate(circularDepthSeriesCount);
          syncDepthSlotLegendLabelsFromTrackConfigs(circularTrackSlots, circularDepthSeriesCount);
        }
        const hasCircularDepthFile = circularDepthEntries.length > 0;
        if (draftDrawsDepth('circular', drawing)) {
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

        const sourceMode = String(drawing.circularConservation.source || '').trim().toLowerCase() === 'upload'
          ? 'upload'
          : 'losat';
        const circularConservationSourceFiles = sourceMode === 'upload'
          ? normalizeFileList(drawing.files.c_conservation_blasts)
          : normalizeFileList(drawing.files.c_conservation_fastas);
        const shouldDrawCircularPairwiseComparisons = circularConservationSourceFiles.length > 0;
        setCircularConservationEnabled(drawing, shouldDrawCircularPairwiseComparisons);

        if (shouldDrawCircularPairwiseComparisons) {
          setProcessingStatus(t('Preparing sequence comparisons for the circular diagram...'));
          // Thresholds are evaluated before conservation work; the request
          // projects the same resolution, and the draft keeps what was typed.
          resolveComparisonThresholds(drawing.adv, 'circular');

          /** @param {OrderedConservationSource<InputFileRef>[]} comparisonEntries */
          const runCircularLosatConservation = async (comparisonEntries) => {
            const circularLosatProgram = normalizeCircularConservationLosatProgram(
              drawing.circularConservation.losat_program
            );
            drawing.circularConservation.losat_program = circularLosatProgram;
            const subjectGencode = normalizePositiveInteger(drawing.circularConservation.subject_gencode, 1);
            /** @param {unknown} comparisonGencode */
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
            // Generate checked above that the Circular input is complete.
            const subjectFile = /** @type {FileLike} */ (inputFiles.fileOf(drawing.inputType.value === 'gb' ? drawing.files.c_gb : drawing.files.c_fasta));
            const subjectFmt = drawing.inputType.value === 'gb' ? 'genbank' : 'fasta';
            const { fasta: subjectFasta, hash: subjectHash } = await readOncePerFile(
              [subjectFile], `ring-subject:${subjectFmt}`, async () => {
                const fasta = await readAllLosatFasta(subjectFile, subjectFmt);
                return Object.freeze({ fasta, hash: await hashText(fasta) });
              }
            );
            const subjectSequenceKey = `circular-subject:${subjectHash}`;
            const sequenceEntriesByKey = new Map([[subjectSequenceKey, subjectFasta]]);
            const cacheMap = new Map(losatCache.value || []);
            /** @type {LosatCacheInfoEntry[]} */
            const cacheInfo = [];
            const losatPairs = [];
            /** @type {LosatQueuedJob[]} */
            const losatJobs = [];
            const pendingJobKeys = new Set();
            const executionMode = getLosatExecutionMode();

            for (let index = 0; index < comparisonEntries.length; index += 1) {
              throwIfGenerationCanceled();
              const comparisonEntry = comparisonEntries[index];
              // A ring source names its bound file.
              const fileObj = /** @type {FileLike} */ (inputFiles.fileOf(comparisonEntry.file));
              // One sequence reader (D12): the Python helper reads FASTA, GenBank
              // or DDBJ and returns the query FASTA that the CLI ring hashes.
              const { fasta: queryFasta } = await readComparisonSequence({ file: fileObj });
              const queryHash = await readComparisonSequenceDigest({ file: fileObj });
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
              setProcessingStatus(t('Preparing comparison search runtime...'));
              const runtime = await prepareLosatRuntime({ includeThreaded: executionMode !== 'serial' }).catch((error) => {
                console.warn('LOSAT runtime warmup failed; execution will report the error.', formatError(error));
                return null;
              });
              if (runtime?.threaded && losatThreadingStatus) {
                const { wasmModule: _wasmModule, ...threadedStatus } = runtime.threaded;
                losatThreadingStatus.value = threadedStatus;
              }
              setProcessingStatus(t('Running {program} comparisons for the circular diagram: {completed}/{total} jobs complete', { program: circularLosatSuffix.toUpperCase(), completed: 0, total: losatJobs.length }));
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
                  setProcessingStatus(t('Running {program} comparisons for the circular diagram: {completed}/{total} jobs complete', { program: circularLosatSuffix.toUpperCase(), completed, total }));
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
              setProcessingStatus(t('Using cached LOSAT comparison results for the circular diagram...'));
            }
            throwIfGenerationCanceled();
            retainCompletedLosatSearch(cacheMap, losatPairs);

            const resolved = [];
            for (const pair of losatPairs) {
              const cached = cacheMap.get(pair.cacheKey);
              const current = isCurrentRawLosatCacheEntry(cached);
              const blastText = current ? cached.text : '';
              const blastPath = `/conservation_blast_${pair.sourceIndex}.txt`;
              // Its name is the Source recipe's name for this ring's generated BLAST table, recorded once the
              // request holds the table in the search frame.
              circularRingSearches.push({
                cacheKey: current ? pair.cacheKey : null,
                path: blastPath,
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
        if (useLinearTrackSlots) {
          linearTrackSlots = normalizeLinearTrackSlots(
            drawing.adv.linear_track_slots,
            drawing.adv.nt,
            drawing.form.linear_track_layout
          );
          const linearDepthSeriesCount = depthTrackMatrixWidth(
            linearSeqs.map((seq) => depthFileSlotsFromValue(seq.depth))
          );
          shapeDepthSeriesForGenerate(linearDepthSeriesCount);
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
          syncDepthSlotLegendLabelsFromTrackConfigs(linearTrackSlots, linearDepthSeriesCount);
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
          setProcessingStatus(t('Preparing comparisons...'));
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
          : ORTHOGROUP_INFERENCE_VERSION;
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
        const collinearMaxParalogLinksPerOrthogroup = paralogLinksFor(
          useProteinBlastp ? blastpMode : '',
          drawing.losat.blastp?.collinearMaxParalogLinksPerOrthogroup
        );
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
            canonicalComparisons: drawing.files.linearCanonicalComparisons,
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
        const refusedLineFontSize = refusedDefinitionLineFontSizePath(drawing.adv.linear_definition_line_styles);
        if (refusedLineFontSize) throw diagnosticError('INPUT_INVALID', { configPath: refusedLineFontSize, reason: 'POSITIVE' });
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

        /** @type {{ reverse: boolean }[]} */
        const viewTransformSpecs = [];
        /** @param {LinearSeq} seq */
        const buildRegionSpec = (seq) => {
          const recordIdRaw = seq.region_record_id ? String(seq.region_record_id).trim() : '';
          const wantsReverse = Boolean(seq.region_reverse);
          // An error names the record card, OFF cards counted.
          const bounds = resolveLinearRegionBounds(seq, drawing.linearSeqs.indexOf(seq));

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

        // LOSAT reads the drawn records, then the OFF records of the same
        // files, which stay in the searched databases (D-07); an index below
        // `losatSearch.drawnCount` is a position in the drawn list.
        const losatSearch = losatSearchSequences(linearSeqs, useLosat ? drawing.linearSeqs : linearSeqs);
        const losatSeqs = /** @type {typeof linearSeqs} */ (losatSearch.sequences);
        /**
         * Each record's LOSAT FASTA in this Generate, and the read of its SHA-256 (OV-1594).
         * @type {Map<number, { entry: LosatSeqEntry, readDigest: () => Promise<string> }>}
         */
        const fastaCache = new Map();
        const sequenceEntriesByKey = new Map();
        if (!reuseResolvedProteinArtifacts) {
          if (useOrthogroupBlastp) {
            workingOrthogroups = [];
          } else {
            workingOrthogroups = [];
            workingSelectedOrthogroupId = '';
            workingSelectedOrthogroupAlignmentFeature = '';
          }
        }
        /** @type {LosatCacheInfoEntry[]} */
        let cacheInfo = [];
        const cacheMap = new Map(losatCache.value || []);
        const derivedCacheMap = useProteinBlastp
          ? new Map(losatDerivedCache.value || [])
          : null;
        const losatTiming = useLosat
          ? {
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

        /** @type {string[]} */
        let proteinRecordInstanceKeys = [];
        /** @type {Readonly<Record<string, unknown>> | null | undefined} */
        let nextResources;
        // The record sources this Generate's request will name, read once (OV-1595).
        const nextRecordResources = () => {
          // The descriptors' payload owners are the cards' files (the store's), not their refs.
          if (nextResources === undefined) {
            nextResources = linearRecordSourceResources(linearSeqs.map((seq) => ({
              ...seq,
              gb: inputFiles.fileOf(seq.gb),
              gff: inputFiles.fileOf(seq.gff),
              fasta: inputFiles.fileOf(seq.fasta)
            })), drawing.inputType.value);
          }
          return nextResources;
        };

        const buildProteinRecordInstanceKeys = async () => {
          const used = new Set();
          return losatSeqs.map((sequence, index) => {
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

        /**
         * @param {number} idx
         * @returns {Promise<LosatSeqEntry>}
         */
        const getSeqEntry = async (idx) => {
          throwIfGenerationCanceled();
          const held = fastaCache.get(idx);
          if (held) return held.entry;
          const startedAt = getNow();
          const fmt = drawing.inputType.value === 'gb'
            ? 'genbank'
            : (useProteinBlastp ? 'gff' : 'fasta');
          const regionSpec = regionSpecs[idx]?.file || null;
          // The FASTA readers below count every entry of a GFF3 input's FASTA file (OV-1681).
          const recordSelector = !useProteinBlastp && drawing.inputType.value === 'gff'
            ? await gffFastaReaderSelector({
                gffFile: inputFiles.fileOf(losatSeqs[idx]?.gff), fastaFile: inputFiles.fileOf(losatSeqs[idx]?.fasta),
                selector: recordSelectors[idx] ?? ''
              })
            : (recordSelectors[idx] ?? '');
          const reverseFlag = /** @type {string} */ ('0');
          // The bytes of the card's files (the card names them by ref); a
          // searched card holds its source file (Generate checked each card).
          const sourceFile = /** @type {FileLike} */ (inputFiles.fileOf(drawing.inputType.value === 'gb'
            ? losatSeqs[idx]?.gb
            : (useProteinBlastp ? losatSeqs[idx]?.gff : losatSeqs[idx]?.fasta)));
          const pairedFastaFile = drawing.inputType.value === 'gff' && useProteinBlastp
            ? inputFiles.fileOf(losatSeqs[idx]?.fasta)
            : null;
          const recordInstanceKey = useProteinBlastp
            ? (proteinRecordInstanceKeys[idx] || `r_${idx + 1}`)
            : '';
          const recordVisibilityRows = useProteinBlastp ? proteinVisibilityRows([recordInstanceKey]) : [];
          // What is read from the source File (or GFF3 and FASTA pair), kept per File (OV-1587): the protein
          // extraction, or the LOSAT FASTA of this record and region.
          const persistentCacheKey = useProteinBlastp
            ? JSON.stringify({
                inputFormat: drawing.inputType.value,
                regionSpec,
                recordSelector,
                recordInstanceKey,
                recordIndex: idx,
                featureVisibility: featureVisibilityCacheKey,
                featureOverrides: recordVisibilityRows,
                proteinMapSchema: 4
              })
            : JSON.stringify({ fmt, regionSpec, recordSelector, reverseFlag });
          const readFiles = /** @type {[FileLike] | [FileLike, FileLike]} */ (
            pairedFastaFile ? [sourceFile, pairedFastaFile] : [sourceFile]
          );
          let extracted = false;
          /** @type {LosatSeqEntry} */
          let entry;
          if (useProteinBlastp) {
            entry = await readOncePerFile(readFiles, `protein-fasta:${persistentCacheKey}`, async () => {
              extracted = true;
              // A source is staged under the resource ids this Generate's request names, so it crosses once and
              // is parsed once for every card and the render (OV-1595); when the request repeats the last one, the
              // Worker already holds its bytes and the records the last Generate parsed (OV-1555). Only when the
              // request cannot be named (a card lacks a file) are the committed drawing's ids used. Another file
              // (of a record turned off only) is sent.
              const committed = requestRecordSource(
                nextRecordResources() ?? (
                  typeof getCommittedCanonicalSession === 'function' ? getCommittedCanonicalSession()?.resources : null
                ),
                sourceFile, pairedFastaFile ?? null
              );
              /** @type {{ role: string, bytes: ArrayBuffer }[]} */
              const helperFiles = [];
              if (!committed) {
                helperFiles.push({ role: 'source', bytes: await cloneFileBytesForTransfer(sourceFile) });
                if (pairedFastaFile) {
                  helperFiles.push({ role: 'fasta', bytes: await cloneFileBytesForTransfer(pairedFastaFile) });
                }
              }
              if (featureVisibilityCacheKey) {
                helperFiles.push({
                  role: 'visibility',
                  bytes: textEncoder.encode(featureVisibilityCacheKey).buffer
                });
              }
              const response = await runDiagramHelperOperation(
                DIAGRAM_HELPER_OPERATIONS.EXTRACT_CDS_PROTEIN_FASTA,
                {
                  files: helperFiles,
                  ...committed,
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
              return Object.freeze({
                fasta: res.fasta,
                recordId: res.record_id || `seq_${idx + 1}`,
                canonicalLength: Number(res.record_length || 0) || getFastaSequenceLength(res.fasta),
                proteinMap: res.protein_map || {},
                proteinCount: res.protein_count || 0,
                proteinSetHash: String(res.protein_set_hash || ''),
                recordAnalysisId: String(res.record_analysis_id || ''),
                recordInstanceKey: String(res.record_instance_key || recordInstanceKey),
                runtimeBindingHash: String(res.runtime_binding_hash || ''),
                displayBindingHash: String(res.display_binding_hash || ''),
                identityManifest: res.identity_manifest || null,
                hash: await hashText(res.fasta || '')
              });
            }, PROTEIN_EXTRACTIONS_PER_FILE);
            if (losatTiming) {
              if (extracted) losatTiming.fastaWorkerFallbacks += 1;
              else losatTiming.proteinExtractionCacheHits += 1;
            }
          } else {
            entry = await readOncePerFile(readFiles, `losat-fasta:${persistentCacheKey}`, async () => {
              extracted = true;
              /** @param {unknown} [reason] Why the browser reader did not serve the file (none: it has no reader for it). */
              const extractWithWorker = async (reason) => {
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
                if (losatTiming) {
                  losatTiming.fastaWorkerFallbacks += 1;
                  console.warn('LOSAT browser FASTA extraction fell back to the diagram Worker.', reason ? formatError(reason) : '');
                }
                return Object.freeze({
                  fasta: res.fasta,
                  recordId: res.record_id || `seq_${idx + 1}`,
                  canonicalLength: Number(res.record_length || 0) || getFastaSequenceLength(res.fasta)
                });
              };
              try {
                // One parse per file (OV-1556); each card slices its record and region from it.
                const { records, parsed } = await readParsedSequenceRecords(sourceFile, fmt);
                if (records) {
                  if (parsed && losatTiming) losatTiming.fastaJsExtractions += 1;
                  return Object.freeze(losatSeqEntryFromRecord(selectParsedRecord(records, recordSelector), regionSpec));
                }
              } catch (fastError) {
                return await extractWithWorker(fastError);
              }
              return await extractWithWorker();
            }, LOSAT_FASTA_READS_PER_FILE);
            if (!extracted && losatTiming) losatTiming.fastaCacheHits += 1;
          }
          const fasta = entry.fasta;
          const proteinHash = entry.hash;
          fastaCache.set(idx, {
            entry,
            // A protein extraction carries its digest. A nucleotide FASTA is digested on first need, once per File
            // and key beside its entry: the records of a searched database that no pair names are never digested.
            readDigest: proteinHash ? async () => proteinHash : async () => {
              let digested = false;
              const hash = await readOncePerFile(readFiles, `losat-fasta-sha256:${persistentCacheKey}`, async () => {
                digested = true;
                const digestStartedAt = getNow();
                const digest = await hashText(fasta);
                if (losatTiming) losatTiming.cacheHashMs += getNow() - digestStartedAt;
                return digest;
              }, LOSAT_FASTA_READS_PER_FILE);
              if (!digested && losatTiming) losatTiming.cacheHashHits += 1;
              return hash;
            }
          });
          if (losatTiming) {
            losatTiming.fastaExtractionMs += getNow() - startedAt;
            losatTiming.totalFastaChars += entry.fasta.length;
          }
          return entry;
        };

        /**
         * The SHA-256 of a record's LOSAT FASTA.
         * @param {number} idx
         * @returns {Promise<string>}
         */
        const getSeqHash = async (idx) => {
          await getSeqEntry(idx);
          return /** @type {{ readDigest: () => Promise<string> }} */ (fastaCache.get(idx)).readDigest();
        };
        /**
         * The key a record's FASTA is passed to the search under: its kind and SHA-256.
         * @param {string} hash
         */
        const losatSequenceKey = (hash) => `${useProteinBlastp ? 'protein' : 'seq'}:${hash}`;

        /** @param {number} idx */
        const getViewTransform = async (idx) => {
          const entry = await getSeqEntry(idx);
          const length = Number(entry.canonicalLength || 0) || getFastaSequenceLength(entry.fasta);
          return {
            length,
            reverse: Boolean(viewTransformSpecs[idx]?.reverse)
          };
        };

        /**
         * @param {string[]} argsKey
         * @param {number} queryIdx
         * @param {number} subjectIdx
         * @returns {Promise<LosatCacheMetadata>}
         */
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

        /**
         * @param {{ cacheKey: string, metadata: LosatCacheMetadata, queryEntry: LosatSeqEntry, subjectEntry: LosatSeqEntry }} candidate
         */
        const tryPromoteLegacyProteinEntry = async ({
          cacheKey,
          metadata,
          queryEntry,
          subjectEntry
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
            migratedFromSchema: LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA
          };
          delete promoted.key;
          delete promoted.filename;
          delete promoted.display;
          cacheMap.set(cacheKey, promoted);
          const verified = getCurrentRawLosatCacheEntry(
            cacheMap,
            cacheKey,
            metadata,
            workingProteinIdentityManifest
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
          // `getCurrentRawLosatCacheEntry` returns the `cacheKey` it was given (its `unknown` parameter widens the type).
          return /** @type {TextKeyedLosatCacheMatch} */ (verified);
        };

        /**
         * @param {LosatJobSpec} spec
         * @param {LosatSeqEntry} queryEntry
         * @param {LosatSeqEntry} subjectEntry
         */
        const buildCacheFilename = (spec, queryEntry, subjectEntry) => {
          const edge = comparisonResolution.edges.find(
            (/** @type {LinearComparisonResolvedEdge} */ candidate) => candidate.edgeKey === spec.edgeKey
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

        /**
         * @param {string[]} arr
         * @param {string} flag
         * @param {unknown} value
         */
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

        /**
         * @param {number} queryIdx
         * @param {number} subjectIdx
         */
        const buildLosatArgs = (queryIdx, subjectIdx) => {
          /** @type {string[]} */
          const args = [];
          if (drawing.losatProgram.value === 'blastn') {
            pushArg(args, '--task', drawing.losat.blastn.task);
          } else if (drawing.losatProgram.value === 'tblastx') {
            pushArg(args, '--query-gencode', losatRecordGencode(losatSeqs[queryIdx]));
            pushArg(args, '--db-gencode', losatRecordGencode(losatSeqs[subjectIdx]));
          } else {
            if (!useOrthogroupBlastp && !useCollinearBlastp) {
              pushArg(args, '--max-hsps-per-subject', 1);
            }
            pushArg(args, '--max-target-seqs', getBlastpCandidateLimit());
          }
          return args;
        };

        regionSpecs = losatSeqs.map((seq) => buildRegionSpec(seq));
        if (useLosat && useProteinBlastp) {
          proteinRecordInstanceKeys = await buildProteinRecordInstanceKeys();
          const proteinRecordIndexes = (useOrthogroupBlastp || useCollinearBlastp)
            ? linearSeqs.map((_, index) => index)
            : Array.from(new Set(comparisonResolution.edges
                .filter((/** @type {LinearComparisonResolvedEdge} */ edge) => edge.source === 'losat')
                .flatMap((/** @type {LinearComparisonResolvedEdge} */ edge) => [edge.queryIndex, edge.subjectIndex])))
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
          setProcessingStatus(t('Preparing LOSAT jobs...'));
          const losatPairs = [];
          /** @type {LosatLinearJob[]} */
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
            ...losatSearch,
            specs: jobSpecs,
            getEntry: getSeqEntry,
            getDigest: getSeqHash,
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
            // `prepareLosatSourceBatches` maps every spec of `jobSpecs` to its batch.
            const batch = /** @type {LosatSourceBatch<typeof spec>} */ (sourcePlan.bySpec.get(spec));
            if (batch.searchContext) cacheMetadata.searchContext = batch.searchContext;
            preparedJobs.push({ spec, losatArgs, cacheMetadata, batch });
          }
          // CLI Sessions can specify a complete two-dimensional row layout.
          // For collinearity, every pair across adjacent rows is displayed;
          // the Web comparison-plan edges alone do not represent those pairs.
          const canonicalGridPair = useCollinearBlastp && drawing.linearRecordLayoutEnabled.value
            ? canonicalGridDisplayPairs(linearSeqs, drawing.linearRecordRows)
            : null;
          const hasCanonicalGridRows = canonicalGridPair !== null;
          const pairNumbers = numberLosatDisplayPairs(
            preparedJobs.map(({ spec }) => spec), comparisonResolution.edges, canonicalGridPair
          );
          const useDerivedProteinPayloadCache = useOrthogroupBlastp || useCollinearBlastp;
          // The converter settings. The derived key and the converter read these same values.
          const derivedSettings = {
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
            explicitDisplayPairs: hasCanonicalGridRows
          };
          // The records the converter reads: every record for similarity groups and Collinear.
          const converterRecordIndexes = useDerivedProteinPayloadCache
            ? linearSeqs.map((_, index) => index)
            : Array.from(new Set(preparedJobs.flatMap(({ spec }) => [spec.queryIndex, spec.subjectIndex])))
              .sort((left, right) => left - right);
          /** @type {readonly string[] | null} */
          let proteinCacheKeys = null;
          // The key of the derived payload these raw results convert to (similarity groups and Collinear).
          let derivedCacheKey = '';
          if (useProteinBlastp && preparedJobs.length > 0) {
            throwIfGenerationCanceled();
            const derivedRecords = [];
            if (useDerivedProteinPayloadCache) {
              for (const recordIndex of converterRecordIndexes) {
                const entry = await getSeqEntry(recordIndex);
                derivedRecords.push({
                  recordIndex,
                  recordInstanceKey: entry.recordInstanceKey,
                  viewTransform: await getViewTransform(recordIndex)
                });
              }
            }
            const keyRequest = {
              pairs: preparedJobs.map(({ cacheMetadata }) => ({
                queryRecordInstanceKey: cacheMetadata.queryRecordInstanceKey,
                subjectRecordInstanceKey: cacheMetadata.subjectRecordInstanceKey,
                expectedOptions: {
                  program: cacheMetadata.program,
                  outfmt: cacheMetadata.outfmt,
                  args: normalizeLosatArgs(cacheMetadata.args),
                  ...(cacheMetadata.searchContext ? { searchContext: cacheMetadata.searchContext } : {})
                }
              })),
              ...(useDerivedProteinPayloadCache
                ? { derived: { settings: derivedSettings, records: derivedRecords, pairs: pairNumbers } }
                : {})
            };
            const keyRequestText = JSON.stringify(keyRequest);
            let builtKeys = getCachedFileRead(proteinCacheKeysCache, workingProteinIdentityManifest, keyRequestText);
            if (!builtKeys) {
              // The merged manifest is frozen, so the helper call reads it without a copy.
              const response = await runDiagramHelperOperation(
                DIAGRAM_HELPER_OPERATIONS.BUILD_PROTEIN_LOSAT_CACHE_KEYS,
                { identityManifest: workingProteinIdentityManifest, ...keyRequest }
              );
              throwIfGenerationCanceled();
              const { keys, derivedKey } = requireProteinLosatCacheKeys(
                response.result, preparedJobs.length, { derived: useDerivedProteinPayloadCache }
              );
              builtKeys = Object.freeze({ keys: Object.freeze(keys), derivedKey });
              setCachedFileRead(proteinCacheKeysCache, workingProteinIdentityManifest, keyRequestText, builtKeys);
            }
            ({ keys: proteinCacheKeys, derivedKey: derivedCacheKey } = builtKeys);
          }

          // The merged manifest is admitted (validated once and frozen), so the raw-entry checks below read its
          // identity from `losat-cache.js` without a per-run index (OV-1553).
          for (const [jobIndex, { spec, losatArgs, cacheMetadata, batch }] of preparedJobs.entries()) {
            throwIfGenerationCanceled();
            const queryEntry = await getSeqEntry(spec.queryIndex);
            throwIfGenerationCanceled();
            const subjectEntry = await getSeqEntry(spec.subjectIndex);
            throwIfGenerationCanceled();
            const cacheKey = useProteinBlastp
              // The loop runs only for a non-empty `preparedJobs`, and then `proteinCacheKeys` was set above for a protein run.
              ? /** @type {readonly string[]} */ (proteinCacheKeys)[jobIndex]
              : await hashText(JSON.stringify(buildLosatCachePayload(cacheMetadata)));
            throwIfGenerationCanceled();
            const queryCanonicalHash = await getSeqHash(spec.queryIndex);
            throwIfGenerationCanceled();
            const subjectCanonicalHash = await getSeqHash(spec.subjectIndex);
            throwIfGenerationCanceled();
            const querySequenceKey = losatSequenceKey(queryCanonicalHash);
            const subjectSequenceKey = losatSequenceKey(subjectCanonicalHash);
            sequenceEntriesByKey.set(querySequenceKey, queryEntry.fasta);
            sequenceEntriesByKey.set(subjectSequenceKey, subjectEntry.fasta);
            let cached = getReusableLosatCacheEntry(
              cacheMap,
              cacheKey,
              cacheMetadata,
              workingProteinIdentityManifest
            );
            if (!cached && useProteinBlastp && !cacheMetadata.searchContext) {
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
              (/** @type {LinearComparisonResolvedEdge} */ edge) => edge.edgeKey === spec.edgeKey
            );
            const { displayPair: isResolvedDisplayPair, pairIndex } = pairNumbers[jobIndex];
            const pair = {
              pairIndex,
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
                pairIndex,
                ordinal: spec.ordinal,
                edgeKey: spec.edgeKey,
                queryUid: spec.queryUid,
                subjectUid: spec.subjectUid,
                queryIndex: spec.queryIndex,
                subjectIndex: spec.subjectIndex,
                cacheKey,
                program: drawing.losatProgram.value,
                querySequenceKey,
                subjectSequenceKey,
                queryCanonicalHash,
                subjectCanonicalHash,
                outfmt: drawing.losat.outfmt || '6',
                extraArgs: losatArgs,
                cacheMetadata,
                batch
              });
            }
          }
          /** @type {LosatSourceJob[]} */
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
            setProcessingStatus(t('Preparing comparison search runtime...'));
            const runtimeWaitStartedAt = getNow();
            await waitForCancelablePromise(losatRuntimeWarmup, generationAbortSignal);
            throwIfGenerationCanceled();
            losatTiming.runtimeWaitMs += getNow() - runtimeWaitStartedAt;
            setProcessingStatus(t('Running LOSAT: {completed}/{total} source jobs complete', { completed: 0, total: sourceJobs.length }));
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
                setProcessingStatus(t('Running LOSAT: {completed}/{total} source jobs complete', { completed, total }));
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
              /** @type {Partial<LosatCacheMetadata>} */
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
            setProcessingStatus(t('Using cached LOSAT results...'));
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
            for (const i of converterRecordIndexes) {
              const entry = await getSeqEntry(i);
              recordPayloads.push({
                recordIndex: i,
                recordId: entry.recordId || `seq_${i + 1}`,
                proteinMap: entry.proteinMap || {},
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
                ? t('Converting LOSAT protein links to collinear ribbons...')
                : t('Converting LOSAT protein hits...')
            );
            /** @type {Record<string, any> | null} */
            let convertedPayload = null;
            if (useDerivedProteinPayloadCache) {
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
                  ...derivedSettings
                }
              );
              // The converter returns a JSON object (convert_losatp_blastp_pairs_to_genomic_payload in python-helpers.js).
              convertedPayload = /** @type {Record<string, any>} */ (response.result);
              const { identityKey } = convertedPayload;
              delete convertedPayload.identityKey;
              // The payload is stored under the key it was looked up by, so the converter must name the same identity.
              if (useDerivedProteinPayloadCache && !convertedPayload.error && identityKey !== derivedCacheKey) {
                throw diagnosticError('RESULT_INVALID', { reason: 'DERIVED_KEY_MISMATCH' }, {
                  stage: 'result-admission', operation: 'convertLosatpPairsToGenomicPayload'
                });
              }
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
              recordGeneratedCliFile(blastPath, converted.tsv || '', {
                name: blastName,
                slot: blastSlot
              });
              if (!useCollinearBlastp) {
                if (typeof converted.canonicalTsv !== 'string') {
                  throw diagnosticError('RESULT_INVALID', { reason: 'MISSING_COMPARISON_TABLE' }, { stage: 'result-admission' });
                }
                resolvedComparisons.push({
                  kind: 'precomputedProteinComparison',
                  edgeKey: sourcePair?.edgeKey || '',
                  ordinal: Number(sourcePair?.ordinal),
                  queryRecordIndex: Number(sourcePair?.queryIndex),
                  subjectRecordIndex: Number(sourcePair?.subjectIndex),
                  canonicalTsv: converted.canonicalTsv
                });
              }
            }
          } else {
            setProcessingStatus(t('Preparing nucleotide comparison results...'));
            for (const pair of losatPairs) {
              throwIfGenerationCanceled();
              const cached = cacheMap.get(pair.cacheKey);
              // Raw search-frame rows; the Python planner projects orientation (PD-OI-073).
              const blastText = isCurrentRawLosatCacheEntry(cached) ? String(cached.text || '') : '';
              const blastPath = `/blast_${pair.pairIndex}.txt`;
              const blastName = pair.filename || getPayloadName(blastPath);
              const blastSlot = `generatedFiles.losat_blasts[${pair.pairIndex}]`;
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
              const written = writtenLegacyProteinArtifacts({
                candidates: candidateOwnerSet.legacyProteinRawCandidates,
                evidence: candidateOwnerSet.legacyProteinDerivedEvidence,
                promoted: legacyPromotionTransaction.map(({ candidateIndex }) => candidateIndex),
                hasCurrentDerived: false
              });
              /** @type {Record<string, unknown>} */
              const promotedProteinIdMap = {};
              legacyPromotionTransaction.forEach(({ proteinIdMap }) => {
                Object.entries(proteinIdMap || {}).forEach(([oldId, runtimeHandle]) => {
                  const previous = promotedProteinIdMap[oldId];
                  if (previous && previous !== runtimeHandle) {
                    throw new Error(`Protein reference '${oldId}' migrated ambiguously.`);
                  }
                  promotedProteinIdMap[oldId] = runtimeHandle;
                });
              });
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
                  legacyProteinRawCandidates: written.candidates,
                  legacyProteinDerivedEvidence: written.evidence,
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
        if (draftDrawsDepth('linear', drawing)) {
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

      throwIfGenerationCanceled();
      setProcessingStatus(t('Preparing render inputs and session...'));
      await nextTick();
      await yieldToTask();
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
        resolvedCircularConservation: canonicalCircularConservation,
        omittedAnnotationRecordKeys: requestDrawing.annotationSets.length > 0
          ? (omittedAnnotationRecordKeys?.() || [])
          : []
      });
      if (useCommittedComparison) {
        if (typeof getCommittedCanonicalSession !== 'function') {
          throw new Error('The preserved comparison owner is unavailable.');
        }
        // Load read the committed request to schema 9 (`readStoredRenderRequest`).
        inheritCommittedComparisonIntent({
          candidate: canonical,
          committed: getCommittedCanonicalSession()
        });
      }
      // A ring LOSAT search reports positions on the whole reference file (D-25); Python puts its rows in the
      // search frame of the records this request draws, as every ring table is (PD-OI-073), before the request
      // renders or is saved. The raw LOSAT cache keeps the search output.
      if (circularRingSearches.length > 0) {
        const ringTables = canonical.renderRequest.mode === 'circular'
          ? canonical.renderRequest.diagramOptions.conservationBlastFiles : null;
        if (!Array.isArray(ringTables) || ringTables.length !== circularRingSearches.length) {
          throw diagnosticError('RESULT_INVALID', {}, RING_TABLE_CONVERSION);
        }
        const converted = await convertRingTablesOncePerFrame(
          /** @type {[FileLike] | [FileLike, FileLike]} */ (drawing.inputType.value === 'gb'
            ? [inputFiles.fileOf(drawing.files.c_gb)]
            : [inputFiles.fileOf(drawing.files.c_gff), inputFiles.fileOf(drawing.files.c_fasta)]),
          {
            renderRequest: canonical.renderRequest,
            // The frame of each drawn record: its source, selector, crop and orientation.
            records: canonical.renderRequest.records.map(({ source, selector, region, presentation }) => (
              [source, selector, region, presentation.reverseComplement]
            )),
            resources: canonical.resources
          },
          ringTables.map(({ resourceId }, index) => ({ resourceId, cacheKey: circularRingSearches[index].cacheKey }))
        );
        throwIfGenerationCanceled();
        canonical.resources = await withConvertedResources(canonical.resources, converted);
        ringTables.forEach(({ resourceId }, index) => {
          const { path, cacheKey: _cacheKey, ...file } = circularRingSearches[index];
          recordGeneratedCliFile(path, bytesToText(base64ToBytes(/** @type {ConvertedResource} */ (converted.get(resourceId)).data)), file);
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
      // Run Info names the files of its exact replay, the canonical Session (`buildRunInfo` keeps only those);
      // the Source recipe carries its own file names (OV-1591).
      const runInfoFileMap = new Map([[canonicalReplayPath, {
        path: canonicalReplayPath,
        name: canonicalReplayName,
        slot: 'generatedFiles.canonical_render_session',
        kind: 'generated'
      }]]);
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
          readResourceRecordCount: (/** @type {string} */ resourceId, /** @type {'genbank' | 'fasta'} */ kind) => (
            readCanonicalResourceRecordCount(canonical.resources, resourceId, kind)
          )
        }));
        throwIfGenerationCanceled();
      }
      /** @type {TimingEntry[]} */
      const postGbdrawTimingEntries = [];
      failureStage = 'render';
      const canonicalExecution = await executeShownDrawingCandidate(drawing, {
        canonical,
        decorationContinuity,
        mode: mode.value,
        kind: 'generate',
        shouldAdmit: colorCandidate.shouldAdmit,
        onProgress: onDiagramProgress,
        prepareCommit: prepareCandidateCommit,
        prepareCommitInput: {
          sourceReplaced,
          featureColorOverrides: drawing.featureColorOverrides,
          featureStrokeOverrides: drawing.featureStrokeOverrides,
          featureOverrides: drawing.featureOverrides,
          legendEntries: drawing.legendEntries.value,
          deletedLegendEntries: drawing.deletedLegendEntries.value,
          dormantLegendEntries: drawing.dormantLegendEntries.value,
          originalLegendOrder: originalLegendOrder.value,
          addedLegendCaptions: drawing.addedLegendCaptions.value,
          unrequestedDepthCaptions: unrequestedDepthCaptions(drawing, canonical.renderRequest),
          legendColorOverrides: drawing.legendColorOverrides,
          legendStrokeOverrides: drawing.legendStrokeOverrides,
          manualSpecificRules: candidateRules
        },
        timingEntries: postGbdrawTimingEntries
      });
      console.info(`gbdraw ${mode.value} typed request render: ${formatDuration(canonicalExecution.elapsedMs)}.`);
      setProcessingStatus(t('Preparing preview...'));
      await nextTick();
      await yieldToTask();
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
        const failure = await failOperation(nameFeaturePlacementFailure(canonicalExecution.engineError,
          canonical.renderRequest, extractedFeatures.value), { handle: committedArtifactHandle,
          restore: restoreCommittedArtifact, isCurrent: isCurrentAlert, isCurrentOperation });
        // The composition root asks how to make the labels fit (LP-D43).
        return failure.status === 'error' && canonicalExecution.labelFit
          ? { ...failure, labelFit: canonicalExecution.labelFit } : failure;
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

      /** @type {(ReturnType<typeof buildRunInfo> & { losatTelemetry?: unknown }) | null} */
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
      const currentOwnerSet = captureGeneratedArtifactOwnerSet(drawing);
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
        labelFitReports: canonicalExecution.labelFitReports,
        specificRules: candidateRules,
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
        generatedMode: drawing.mode,
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
        drawing,
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
      // The compiler keeps the edits that hold a color (`normalizedFeatureColorOverrides`).
      replaceFeatureFills(drawing.slices.paint, /** @type {Record<string, FeatureFill>} */ (cloneJsonValue(candidateCommit.featureColorOverrides, {})));
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
        featureEdits: drawing.slices.featureEdits,
        paint: drawing.slices.paint,
        notices: canonicalExecution.featureIdentityNotices,
        replacedRecordKeys: replacedFeatures.filter((feature) => feature.scope === canonical.renderRequest.mode)
          .map((feature) => String(feature.record_key)),
        previousRecords: previousRequestRecords,
        currentRecords: canonical.renderRequest.records || [],
        retainedRecordKeys: offRecordRequestKeys(canonical.renderRequest.mode, drawing.recordsOff),
        biologicalFeatures: candidateBiologicalFeatures
      }) : 0;
      // OV-120: a source replacement also retires the renames of the Legend
      // rows the new Result does not draw (as the OV-84 strokes).
      if (sourceReplaced) replaceDormantLegendRows(drawing.slices.paint, []);
      if (typeof setGeneratedArtifactIdentity === 'function') {
        setGeneratedArtifactIdentity(generationResponse.artifactIdentity, {
          results: candidateCommit.results
        });
      }
      if (typeof adoptCanonicalRenderArtifacts === 'function') {
        // The committed request keeps the computed typed analysis resource, so
        // target-only renders and Similarity alignment use it without LOSATP.
        adoptCanonicalRenderArtifacts(canonical, { adoptOwnedRequest: true, drawing });
      }
      settlePendingLiveSettings(drawing);
      if (!useCommittedComparison && drawing.importedComparisonIntent) {
        Object.assign(
          drawing.importedComparisonIntent,
          createImportedComparisonIntentState(),
          { disposition: IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE }
        );
      }
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

  // UJ-10: the Results a successful Generate showed and the draft it read.
  /** @type {{ results: any, signature: string } | null} */
  let shownGenerateDraft = null;
  /**
   * @param {Record<string, any> | null} [comparisonPlanSnapshot]
   * @param {HistoryArtifactHandle | null} [generatedArtifactHandle]
   * @param {Record<string, any> | null} [comparisonExecution]
   * @param {{ prepareGenerate?: (() => Promise<Record<string, any>>) | null,
   *   afterGenerate?: ((outcome: Record<string, any> | null, drawing: DrawingState) => any) | null }} [options]
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
    /** @type {HistoryArtifactHandle | null} */
    let beforeHandle = null;
    let initialResults = results.value;
    /** @type {string | null} */
    let historyRecovery = null;
    // UJ-10: whether the draft this run reads is the one of the shown Result.
    // It is read when the run settles its draft or when Cancel is observed,
    // whichever comes first, so an edit committed before either counts.
    const readRunDraftIsShown = () => {
      const signature = readDraftSignature?.() || '';
      return Boolean(signature) && shownGenerateDraft?.signature === signature
        && shownGenerateDraft.results === results.value;
    };
    /** @param {boolean} draftIsShown */
    const preserveCanceledResult = (draftIsShown) => {
      failedGeneratePreservedResult.value = results.value.length > 0 && (draftIsShown ? 'current' : true);
    };
    let runDraftIsShown = false;
    processing.value = true;
    processingStatus.value = t('Preparing input files...');
    generationCancelRequested.value = false;
    recordSessionLifecycleEvent('generate.processing-published');
    try {
      const decorationContinuity = captureDecorationContinuity(getCommittedCanonicalSession?.(), projectCompositionRecordIdentity, drawing);
      await nextTick();
      await waitForAfterPaint();
      // A rerender this run waited for may commit a Result; a failure keeps that one.
      if (activeReflowRun) {
        await activeReflowRun;
        initialResults = results.value;
      }
      if (!isCurrentOperation()) return { status: 'stale' };
      recordSessionLifecycleEvent('generate.paint-opportunity-completed');
      // Cancel is honored before rendering, including during preparation.
      const cancelBeforeRender = () => {
        if (!generationCancelRequested.value) return null;
        processingStatus.value = t('Canceled.');
        outcome = { status: 'canceled' };
        preserveCanceledResult(readRunDraftIsShown());
        return outcome;
      };
      if (cancelBeforeRender()) return outcome;
      // UI-07: the active mode's input check runs before anything here starts
      // the diagram Worker, so a missing input is named at once.
      try {
        assertDrawingInputs?.(drawing);
      } catch (cause) {
        outcome = await failOperation(cause, { handle: null, stage: 'request-validation',
          isCurrent: isCurrentAlert, isCurrentOperation, recovery: results.value.length ? 'preserved' : 'no-result' });
        return outcome;
      }
      await settlePendingInputs();
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
      /** @param {HistoryArtifactHandle | null} handle */
      const execute = (handle) => {
        beforeHandle = handle || generatedArtifactHandle;
        // History has settled the draft this run reads.
        runDraftIsShown = readRunDraftIsShown();
        return runAnalysisInternal(drawing, {
          decorationContinuity, comparisonPlanSnapshot, generatedArtifactHandle: beforeHandle,
          comparisonExecution, isCurrentOperation, isCurrentAlert
        });
      };
      outcome = typeof runGeneratedArtifactReplacement === 'function'
        ? await runGeneratedArtifactReplacement(
            msg('Generate diagram'),
            execute,
            {
              shouldCommit: (result) => result?.status === 'ok',
              // A Generate that replaced a source removes the feature placements
              // it no longer resolves (Q3 = A). Placements are draft intent, not
              // part of the Result, so Undo and Redo restore them, and the
              // removal count, with this step.
              captureIntentCheckpoint: () => ({
                placements: FEATURE_EDITS_SLICE.capture(drawing.slices.featureEdits).featurePlacementOverrides,
                removed: Number(featureEditRemovalCount.value) || 0
              }),
              // The checkpoint is what `captureIntentCheckpoint` above returned.
              restoreIntentCheckpoint: (/** @type {{ placements: FeaturePlacementDraft, removed: number }} */ { placements, removed }) => {
                FEATURE_EDITS_SLICE.install(drawing.slices.featureEdits, { featurePlacementOverrides: placements });
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
        : await execute(generatedArtifactHandle || await captureGeneratedArtifactHandle(drawing));
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
        preserveCanceledResult(runDraftIsShown);
      } else if (outcome?.status === 'ok') {
        failedGeneratePreservedResult.value = false;
        if (generationFailureRecovery) generationFailureRecovery.value = null;
      }
      await afterGenerate?.(outcome, drawing);
      if (outcome?.status === 'ok') shownGenerateDraft = { results: results.value, signature: readDraftSignature?.() || '' };
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
   * @param {{ canonical: Record<string, any>, decorationContinuity?: any, generatedArtifactHandle?: HistoryArtifactHandle | null,
   *   commitIntent?: (() => any) | null, alignmentResetBefore?: any, alignmentResetReceipt?: any, operation?: string,
   *   isCurrentOperation?: () => boolean, isCurrentAlert?: () => boolean }} options
   */
  const runCommittedCanonicalCandidateInternal = async (drawing, {
    canonical: committedCandidate,
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
      || await captureGeneratedArtifactHandle(drawing);
    /** @type {ReturnType<typeof generatedArtifactTransactionOwner.build> | null} */
    let activatedCandidate = null;
    /** @type {ReadyReceipt | null} */
    let acceptedReadyReceipt = null;
    const restoreCommittedArtifact = async () => {
      if (!activatedCandidate) return false;
      if (captureGeneratedArtifactOwnerSet(drawing).results !== activatedCandidate.ownerSet.results) return false;
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
      processingStatus.value = t('Canceled.');
      return { status: 'canceled' };
    };
    try {
      // The candidate draws the live editor paint, as the label rerender does (OV-1083).
      const paint = await prepareCommittedPaintCandidate(drawing, /** @type {CommittedCanonicalSession} */ (committedCandidate),
        () => generationToken === latestGenerationToken && !generationCancelRequested.value);
      if (!paint) return generationCancelRequested.value ? finishCanceled() : { status: 'stale' };
      const { canonical, drawnPaint } = paint;
      /** @type {TimingEntry[]} */
      const timingEntries = [];
      const execution = await executeShownDrawingCandidate(drawing, {
        canonical,
        decorationContinuity,
        mode: canonical.renderRequest.mode,
        kind: 'target-record-transform',
        shouldAdmit: () => paint.shouldAdmit() && !generationCancelRequested.value,
        onProgress: (progress) => {
          const message = diagramProgressMessage(progress);
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
          unrequestedDepthCaptions: unrequestedDepthCaptions(drawing, canonical.renderRequest),
          legendColorOverrides: drawing.legendColorOverrides,
          legendStrokeOverrides: drawing.legendStrokeOverrides,
          manualSpecificRules: paint.rules
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
        const failure = await failOperation(execution.engineError, { handle: committedArtifactHandle, operation,
          restore: restoreCommittedArtifact, isCurrent: isCurrentAlert, isCurrentOperation });
        return failure.status === 'error' && execution.labelFit
          ? { ...failure, labelFit: execution.labelFit } : failure;
      }
      processingStatus.value = t('Preparing preview...');
      const candidateCommit = execution.commit;
      const candidateCatalogAdmission = execution.catalogAdmission;
      const candidateCatalog = execution.catalog;
      const candidateGroups = Array.isArray(candidateCommit.featureState.orthogroups)
        ? candidateCommit.featureState.orthogroups
        : [];
      // Only Align sets `alignmentResetBefore`, and it projects its plan onto `canonical`, so `canonical` holds the plan and the translations.
      const receipt = alignmentResetBefore
        ? await buildSimilarityAlignmentResetReceipt({before:alignmentResetBefore, after:/** @type {AlignmentResetPlannedCanonical} */ (canonical)})
        : alignmentResetReceipt === undefined ? state.similarityAlignmentResetReceipt?.value
          : alignmentResetReceipt;
      await validateSimilarityAlignmentResetReceipt(receipt, /** @type {AlignmentResetCanonical} */ (canonical));
      const currentOwnerSet = captureGeneratedArtifactOwnerSet(drawing);
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
        labelFitReports: execution.labelFitReports,
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
        bindingOptions: { isIncrementalEdit: false, drawnPaint },
        isCurrent: () => generationToken === latestGenerationToken
          && !generationCancelRequested.value
          && Number(selectedResultIndex.value) === nextSelectedResultIndex
      });
      activatedCandidate = generatedCandidate;
      generatedArtifactTransactionOwner.activate(generatedCandidate, {
        drawing,
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
        adoptCanonicalRenderArtifacts(canonical, { adoptOwnedRequest: true, drawing });
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
   *   captureIntentCheckpoint?: HistoryStepOptions['captureIntentCheckpoint'], restoreIntentCheckpoint?: HistoryStepOptions['restoreIntentCheckpoint'],
   *   commitIntent?: (() => any) | null, alignmentResetBefore?: any, alignmentResetReceipt?: any, operation?: string }} options
   */
  const runCommittedCanonicalCandidate = async ({
    canonical,
    label = msg('Rotate record to feature'),
    captureIntentCheckpoint = null,
    restoreIntentCheckpoint = null,
    commitIntent = null,
    alignmentResetBefore = null,
    alignmentResetReceipt = undefined,
    operation = 'generate'
  }) => {
    // A target-record transform and a label reflow draw the committed Results
    // again, with the settings and edits of the drawing they belong to.
    const drawing = displayedArtifactDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    /** @type {Record<string, any> | null} */
    let outcome = null;
    const operationId = ++latestOperationId;
    const isCurrentOperation = () => operationId === latestOperationId;
    const previousAlert = errorLog.value;
    const isCurrentAlert = () => isCurrentOperation()
      && (errorLog.value === previousAlert || errorLog.value === null);
    /** @type {HistoryArtifactHandle | null} */
    let beforeHandle = null;
    let initialResults = results.value;
    /** @type {string | null} */
    let historyRecovery = null;
    processing.value = true;
    processingStatus.value = t('Preparing target record...');
    generationCancelRequested.value = false;
    try {
      const decorationContinuity = captureDecorationContinuity(getCommittedCanonicalSession?.(), projectCompositionRecordIdentity, drawing);
      await nextTick();
      await waitForAfterPaint();
      // A rerender this run waited for may commit a Result; a failure keeps that one.
      if (activeReflowRun) {
        await activeReflowRun;
        initialResults = results.value;
      }
      if (!isCurrentOperation()) return { status: 'stale' };
      /** @param {HistoryArtifactHandle | null} handle */
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
        : await execute(await captureGeneratedArtifactHandle(drawing));
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
      processingStatus.value = t('Canceling generation...');
    }
    return cancelDiagramGeneration();
  };

  // A reflow keeps the Result it draws, so the preview binder reports a forced
  // label that the Result does not draw instead of failing the binding, with
  // the diagnostic Generate raises for it (R6, OV-06). The rerender also draws
  // the Legend again, so the Legend it draws is the generated inventory, as
  // after Generate (OV-42, OV-43): the next edit compares against it.
  /**
   * @param {Extract<Awaited<ReturnType<typeof executeCanonicalRenderCandidate>>, { status: 'ok' }>['commit']} commit The rerender's commit.
   * @param {number} resultIndex
   * @param {() => boolean} isCurrentReflow
   * @param {LabelDiagramOptions | null | undefined} diagramOptions The `diagramOptions` of the request the rerender drew.
   * @param {EditorPaintState} drawnPaint
   * @param {number} historyRequestedAt The History step serial its request was made at.
   */
  const expectReflowBindings = (commit, resultIndex, isCurrentReflow, diagramOptions, drawnPaint, historyRequestedAt) => {
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
        drawnPaint,
        historyRequestedAt,
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

  // The automatic rerender re-renders the committed Session with the current
  // editor tables and the draft values of the live setting rows its requests
  // named (R1(c), N-16); no other draft setting reaches the Result.
  /**
   * @param {DrawingState} drawing
   * @param {{ requestId: number, historyRequestedAt: number, decorationContinuity: CanonicalRenderCandidateExecutionOptions['decorationContinuity'], settings: readonly string[], editor: boolean }} options
   *   `settings`: the live setting request paths it applies; `editor`: an editor edit asked for it too.
   * @returns {Promise<{ status: string, rerun?: boolean, error?: UserFacingError | null }>}
   */
  const runLabelReflowCandidate = async (drawing, { requestId, historyRequestedAt, decorationContinuity, settings, editor }) => {
    const committed = getCommittedCanonicalSession?.();
    if (!committed) {
      labelReflowLastError.value = liveEditFailure(formatError(diagnosticError('LIVE_EDIT_REQUIRES_GENERATE')));
      return { status: 'skipped' };
    }
    // It commits only onto the committed Session it drew from (an Undo of a
    // Generate step may replace it meanwhile): otherwise it is stale and runs
    // again from the current one.
    const committedRevision = readCommittedCanonicalRevision();
    /** @type {LiveSettingWrite[]} */
    let liveSettingWrites;
    try {
      liveSettingWrites = projectLiveSettingWrites({ committed, drawing, requestPaths: settings });
    } catch (e) {
      const error = formatError(e, 'generate', 'render');
      labelReflowLastError.value = liveEditFailure(error);
      return { status: 'error', error };
    }
    // Values equal to the committed ones (a Generate applied them) draw nothing.
    if (!editor && liveSettingWrites.length === 0) return { status: 'skipped' };
    // The reflow draws the committed Session's Results again, so they keep their names (OV-136).
    const committedResultNames = results.value.map((/** @type {{ name: string }} */ result) => result.name);
    const generationToken = ++latestGenerationToken;
    clearLabelBuildNotices({ rerender: true });
    skipCaptureBaseConfig.value = true;
    const isRequested = () => generationToken === latestGenerationToken && requestId >= reflowSupersededAt;
    const isCurrent = () => isRequested() && readCommittedCanonicalRevision() === committedRevision;
    // A candidate whose rule inputs changed while it was prepared or drawn (an
    // edit the Result shows live) is not admitted. No newer request replaces
    // it, so the rerender it was asked for (a Legend relabel) runs again (OV-351).
    const notAdmitted = () => ({ status: 'stale', rerun: isRequested() });
    try {
      const paint = await prepareCommittedPaintCandidate(drawing, committed, isCurrent, liveSettingWrites);
      if (!paint) return notAdmitted();
      const { canonical, drawnPaint, rules: candidateRules } = paint;
      /** @type {TimingEntry[]} */
      const timingEntries = [];
      const execution = await executeShownDrawingCandidate(drawing, {
        canonical,
        decorationContinuity,
        mode: canonical.renderRequest.mode,
        kind: 'reflow',
        shouldAdmit: paint.shouldAdmit,
        prepareCommit: prepareReflowResultCommit,
        prepareCommitInput: {
          featureColorOverrides: drawing.featureColorOverrides,
          featureStrokeOverrides: drawing.featureStrokeOverrides,
          featureOverrides: drawing.featureOverrides,
          legendEntries: drawing.legendEntries.value,
          deletedLegendEntries: drawing.deletedLegendEntries.value,
          dormantLegendEntries: drawing.dormantLegendEntries.value,
          originalLegendOrder: originalLegendOrder.value,
          addedLegendCaptions: drawing.addedLegendCaptions.value,
          unrequestedDepthCaptions: unrequestedDepthCaptions(drawing, canonical.renderRequest),
          legendColorOverrides: drawing.legendColorOverrides,
          legendStrokeOverrides: drawing.legendStrokeOverrides,
          manualSpecificRules: candidateRules
        },
        timingEntries,
        resultNames: committedResultNames
      });
      console.info(`gbdraw ${canonical.renderRequest.mode} typed request render: ${formatDuration(execution.elapsedMs)}.`);
      if (execution.status === 'superseded' || !isCurrent()) return notAdmitted();
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
      // OV-873: the new Results have new identities, which key the Legend and
      // title moves History records; the rerender draws the same moves, so a
      // step begun after its request takes the change as drawn (D-40).
      // OV-111: they install into the artifact of the drawing whose Result they draw.
      adoptDrawnWrite(() => {
        installRedrawnArtifact(drawing.id, {
          featureCatalog: execution.catalog,
          extractedFeatures: execution.commit.featureState.extractedFeatures,
          biologicalFeatures: execution.commit.featureState.biologicalFeatures,
          results: execution.commit.results,
          // The notice follows the Results the rerender drew: it goes when they fit (LP-D48).
          labelFitReports: execution.labelFitReports
        }, {
          installResults(nextResults) {
            results.value = nextResults;
          }
        });
        // The preview owner selects the Result (R13). It registers readiness only
        // when fewer Results leave the selection out of range; the reflow's own
        // expectation, registered next, replaces it, so the binder reports a
        // forced label the Result does not draw (OV-06).
        if (execution.commit.results.length > 0) {
          previewRuntime.selectResult(nextSelectedResultIndex);
        }
      }, { requestedAt: historyRequestedAt });
      // The committed Session names what the new Results draw (PD-OI-102); the
      // binder reports for these Results while that Session stays committed.
      commitLiveSettingValues(liveSettingWrites);
      const drawnRevision = readCommittedCanonicalRevision();
      const isShownReflow = () => isRequested() && readCommittedCanonicalRevision() === drawnRevision;
      if (canonical.renderRequest.mode === 'circular') {
        state.generatedCircularPlotTitlePosition.value = normalizeCircularPlotTitlePosition(
          canonical.renderRequest.diagramOptions.output.plotTitlePosition
        );
      }
      expectReflowBindings(
        execution.commit, nextSelectedResultIndex, isShownReflow, canonical.renderRequest.diagramOptions, drawnPaint,
        historyRequestedAt
      );
      logPostGbdrawTimings(timingEntries);
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

  /**
   * @param {{ settings?: readonly string[] }} [options] `settings`: the live
   *   setting request paths of a request to apply live values alone; without
   *   it, an editor edit's rerender.
   */
  const runLabelReflow = async ({ settings } = {}) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!settings) editorRerenderRequested = true;
    pendingReflowRequestId += 1;
    if (!settings || settings.some((path) => activeLiveSettings.includes(path))) reflowSupersededAt = pendingReflowRequestId;
    pendingReflowRequestedAt = readHistoryStepSerial();
    if (activeReflowRequestId !== 0) return;
    if (reflowWaits()) {
      reflowDeferred = true;
      return;
    }

    labelReflowProcessing.value = true;
    activeReflowRun = (async () => {
      // Generate takes priority: an iteration never starts while it runs.
      while (!reflowWaits() && activeReflowRequestId < pendingReflowRequestId) {
        activeReflowRequestId = pendingReflowRequestId;
        const historyRequestedAt = pendingReflowRequestedAt;
        // Each iteration draws the drawing of the committed Session it starts from.
        const drawing = displayedArtifactDrawing();
        let decorationContinuity;
        try {
          decorationContinuity = captureDecorationContinuity(getCommittedCanonicalSession?.(), projectCompositionRecordIdentity, drawing);
        } catch (error) {
          labelReflowLastError.value = liveEditFailure(formatError(error));
          return;
        }
        const editor = editorRerenderRequested;
        editorRerenderRequested = false;
        const pending = pendingLiveSettingsOf(drawing);
        const settings = [...pending];
        activeLiveSettings = settings;
        const outcome = await runLabelReflowCandidate(drawing, {
          decorationContinuity,
          requestId: activeReflowRequestId,
          historyRequestedAt,
          settings,
          editor
        });
        // A stale iteration leaves its values and its editor request to the
        // next one; an applied, failed or empty one is done with them (a
        // failed value waits for its next change or Generate, PD-OI-102).
        if (outcome.status === 'stale' || outcome.status === 'canceled') editorRerenderRequested ||= editor;
        else settings.forEach((path) => pending.delete(path));
        if (outcome.rerun && !processing.value && !state.sessionOperationAvailability?.()) pendingReflowRequestId += 1;
      }
      if (reflowWaits() && activeReflowRequestId < pendingReflowRequestId) {
        reflowDeferred = true;
      }
    })().finally(() => {
      activeReflowRequestId = 0;
      activeLiveSettings = [];
      labelReflowProcessing.value = false;
      activeReflowRun = null;
    });
    await activeReflowRun;
  };

  // A Generate drew the drawing's draft: the paths whose values now equal the
  // committed request are done; a value changed after its request was built
  // stays for the rerender that follows it.
  /** @param {DrawingState} drawing */
  const settlePendingLiveSettings = (drawing) => {
    const pending = pendingLiveSettings.get(drawing);
    const committed = getCommittedCanonicalSession?.();
    if (!pending || !committed) return;
    for (const path of [...pending]) {
      try {
        if (projectLiveSettingWrites({ committed, drawing, requestPaths: [path] }).length === 0) pending.delete(path);
      } catch {
        // A value Generate would refuse stays for the rerender, which reports it.
      }
    }
  };
  // A Load replaces both drawings' committed Sessions and applies nothing (6.7).
  const discardLiveSettingRequests = () => { pendingLiveSettings.clear(); };

  // R1 (c), PD-OI-102: an action committed the values of these live setting
  // request paths in the drawing of `mode` (the shown one). Its Result shows
  // them after the automatic rerender applies them; a mode without a Result
  // applies nothing until its Generate (R15).
  /** @param {readonly string[]} requestPaths @param {'circular' | 'linear'} mode */
  const requestLiveSettingApply = (requestPaths, mode) => {
    const committed = getCommittedCanonicalSession?.();
    if (requestPaths.length === 0 || !committed || committed.renderRequest?.mode !== mode) return false;
    const pending = pendingLiveSettingsOf(displayedArtifactDrawing());
    requestPaths.forEach((path) => pending.add(path));
    labelReflowLastError.value = null;
    void runLabelReflow({ settings: requestPaths });
    return true;
  };

  // Called where a run clears `processing`, and when a hold ends: the requests
  // held back run once. An editor request held back kept its flag; live values
  // the Generate committed are equal now, so a replay with nothing left is
  // skipped.
  const replayDeferredLabelReflow = () => {
    if (!reflowDeferred || reflowWaits()) return;
    reflowDeferred = false;
    void runLabelReflow({ settings: [] });
  };
  // OV-1076: Reset Settings asks for the rerender while it writes (its label
  // projection), again for the rules it removed, and its History checkpoint
  // asks to apply the live values it changed. A rerender started midway reads
  // a Legend the action has not finished writing and runs again. The action
  // holds the requests and one rerender draws what it wrote.
  /**
   * @template T
   * @param {() => T} run
   * @returns {Promise<Awaited<T>>}
   */
  const holdAutomaticRerender = async (run) => {
    rerenderHolds += 1;
    try {
      return await run();
    } finally {
      rerenderHolds -= 1;
      replayDeferredLabelReflow();
    }
  };

  const downloadLosatCache = async () => {
    if (!losatCacheInfo.value || losatCacheInfo.value.length === 0) return;
    const cacheMap = losatCache.value;
    if (!cacheMap || cacheMap.size === 0) return;
    const identityManifest = cloneJsonData(proteinIdentityManifest.value);
    const hydratedEntries = (await Promise.all(losatCacheInfo.value.map(async (/** @type {LosatCacheInfoEntry} */ entry, /** @type {number} */ idx) => {
      const cached = cacheMap.get(entry.key);
      if (!isCurrentRawLosatCacheEntry(cached)) return null;
      return {
        entry,
        cached,
        idx,
        hydrated: await hydrateLosatDownloadText(entry.key, cached, identityManifest)
      };
    }))).filter(Boolean);
    const totalBytes = totalHydratedLosatExportBytes(
      hydratedEntries.map((item) => item.hydrated)
    );

    if (losatExportBytesExceedLimit(totalBytes) && !await askBeforeLargeDownload(askChoice, {
      title: t('Download large LOSAT tables?'),
      description: t('The raw LOSAT TSV files total about {size} MB.', { size: (totalBytes / (1024 * 1024)).toFixed(1) })
    })) return;

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
    legacyProteinRawCandidates.value = { schema: LEGACY_PROTEIN_CANDIDATE_SCHEMA, entries: [] };
    legacyProteinDerivedEvidence.value = { schema: LEGACY_LOSAT_DERIVED_CACHE_SCHEMA, entries: [] };
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
    holdAutomaticRerender,
    requestLiveSettingApply,
    discardLiveSettingRequests,
    refreshCircularRecordOrder,
    releaseSourceInputFailure,
    downloadCliHelperFiles,
    downloadLosatCache,
    downloadLosatPair,
    setLosatPairFilename,
    clearLosatCache,
    getLosatPairDefaultName
  };
};
