// @ts-check
import { writeCanonicalRecordReverseComplement } from './record-display-model.js';
import {
  canonicalFeatureOverrides,
  canonicalFeaturePlacements,
  featureDraftMap,
  featureIdentityKeyOf,
  requestFeatureOverrides,
  requestFeaturePlacements
} from './feature-placement.js';
import { resultRenderedFeatures } from './feature-catalog.js';
export { canonicalFeaturePlacements } from './feature-placement.js';
import { buildDefaultColorOverrideTsv, normalizePaletteColors } from '../utils/color-utils.js';
import {
  EDITOR_TABLE_OPTIONS,
  editorTableRef,
  parseColorTable,
  parsePriorityRules,
  parseSpecificRules,
  parseWhitelistRules,
  serializeLabelWhitelistRules,
  serializeQualifierPriorityRules,
  serializeSpecificRules
} from './file-imports.js';
import {
  buildBulkLabelProjection,
  parseLabelOverrideTsv,
  serializeLabelOverrideRows
} from './label-override-table.js';
import {
  parseFeatureVisibilityRules,
  serializeFeatureVisibilityRules
} from './feature-visibility.js';
import {
  circularRendererName,
  applyCircularGeometryShortcuts,
  buildCircularTrackSlotPayload,
  CIRCULAR_TRACK_RENDERERS,
  createDefaultCircularTrackSlots,
  hasCircularGeometryShortcuts,
  inferLegacyAxisIndexFromFeature,
  parseCircularTrackSlotSpecs
} from './circular-track-slot-model.js';
import { countGenBankRecords } from './genbank-header.js';
import { downloadSafeName } from '../utils/download-names.js';
import { projectCircularMeasureDraft } from './circular-track-measure.js';
import {
  buildLinearTrackSlotPayload,
  LINEAR_TRACK_RENDERERS,
  LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  migrateLinearTrackSlotsToCurrentSchema,
  parseLinearTrackSlotSpecs
} from './linear-track-slot-model.js';
import {
  activeDepthTrackIndices,
  isRecordMajorDepthFileMatrix,
  normalizeRecordMajorDepthFileRows,
  parseDepthTrackIndexIdentity
} from './depth-track-state.js';
import {
  buildDisambiguatedRecordEntries,
  isSingleCircularJourney,
  resolveCircularRequestRecordSet,
  resolveDisambiguatedRecordSelection
} from './record-options.js';
import { circularRecordRequestKey, circularSingleRecordKey } from './record-draw-selection.js';
import {
  orderedConservationSources,
  orderedOptionalConservationFiles
} from './conservation-series.js';
import {
  resolveLinearRecordEffectiveDefinition,
  resolveLinearRecordEffectiveSubtitle
} from './linear-sources.js';
import {
  linearRecordLayoutHasSharedRow,
  resolveEffectiveLinearRecordRows
} from './linear-record-layout.js';
import { resolveLinearLabelVisibility } from './linear-label-visibility.js';
import {
  assertValidCustomTrackPlan,
  validateCustomTrackPlan,
  validateTrackSlotBindingInvariants
} from './track-slot-validation.js';
/** @import { CustomTrackPlan, CustomTrackPlanInput } from './track-slot-validation.js' */
import { annotationOptionsPayload, draftAnnotationSetsOfRequest } from './annotation-state.js';
import { classifyOptionalNumber, classifyOptionalPositiveNumber, projectOptionalNumber } from '../utils/optional-positive-number.js';
import { diagnosticError, normalizeCaughtError } from '../utils/error-normalization.js';
import { loadRepair } from '../utils/load-repairs.js';
import { repairByDescription } from '../utils/contract-repair.js';
import { isPlainObject } from '../utils/contract-shape.js';
import { CANONICAL_REQUEST_DESCRIPTION } from '../contracts/request-shape.generated.js';
import { DEFAULT_LINEAR_RECORD_GAP_PX, draftLinearRecordGap } from './config-field-readers.js';
import { ownEntry } from '../utils/own-entry.js';
import { toRaw } from '../utils/vue-raw.js';
import { materializeLegacySimilarityAlignment } from './legacy-similarity-alignment.js';
import { liftLegacyRenderRequest } from './legacy-render-request.js';
import {
  CANONICAL_REQUEST_SCHEMA,
  CONFIG_UNSET_NUMBER_PATHS,
  CONFIG_UNSET_NUMBER_SPELLINGS,
  SIMILARITY_ALIGNMENT_PLAN_SCHEMA,
  SUPPORTED_CANONICAL_REQUEST_SCHEMAS
} from '../contracts/schema.generated.js';
import {
  arrowHeadLengthRatioForState,
  defaultFeatureRendering,
  normalizeArrowHeadLengthRatio,
  normalizeArrowShaftWidthRatio,
  normalizeFeatureRenderingMap
} from '../utils/feature-rendering.js';
import {
  comparisonFiltersForMode,
  effectiveLinearAxisColor,
  isModeProfileDefault,
  linearAxisRulerApplies,
  MODE_DEFAULT_FEATURE_TYPES,
  modeProfile,
  resolveComparisonThresholds,
  trackDefaultsForMode
} from '../mode-profiles.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';
import {
  createDefaultLayoutPreferences,
  legendPositionOfMode,
  updateActiveLayoutPreference
} from './layout-preferences.js';
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
  requireCurrentCollinearInferOrthogroups,
  requireCurrentCollinearSearchScope,
  requireCurrentCollinearUnitMode,
  requireCurrentLinearLabelPlacement,
  requireCurrentLinearTrackLayout,
  requireCurrentOrthogroupMemberMaxHits,
  requireCurrentOrthogroupMembershipMode,
  requireCurrentProteinBlastpCandidateLimit,
  requireCurrentProteinBlastpMaxHits,
  requireCurrentProteinBlastpMode,
  requireCurrentWebStateFieldNames
} from './current-option-values.js';
import {
  readCollinearAnchorMode,
  readCollinearColorMode,
  readCollinearSearchScope,
  readCollinearUnitMode,
  readOrthogroupMembershipMode
} from './losat-normalization.js';
import {
  canonicalComparisonResourceKind,
  isResourceBackedCanonicalComparison,
  mapResourceBackedCanonicalComparison
} from './canonical-comparisons.js';
import {
  base64ToBytes,
  bytesToBase64,
  bytesToText,
  getSessionResourceSource,
  readFileText,
  textToBytes
} from './file-content-cache.js';
import {
  adoptedSessionResourceDescriptor,
  adoptCurrentSessionResources,
  createCombinedSessionResourceFileView,
  createSessionResourceFileView,
  validateWebFileBindings
} from './session-resource-backing.js';
import { normalizeLinearComparisonPlan } from './linear-comparisons.js';
import { createDefaultLosatExecution } from './session-active-config-contract.js';
import { LINEAR_SEQ_READ_FIELDS } from './config-field-readers.js';
import {
  getResourcePayloadOwner,
  setResourcePayloadOwner
} from './resource-payload-owner.js';
import { sha256Hex } from './byte-utils.js';
import { cloneJsonData } from './json-clone.js';
import { isCanonicalResourceReferenceField } from './canonical-resource-references.js';
import { recordStructuralMetric } from './runtime-test-hooks.js';
import { bindRecordDisplayDrafts, recordDisplayKey, requestedRecordTransform } from './record-display-model.js';

/** @import { LoadRepair, Repaired } from '../utils/load-repairs.js' */
/** @import { ContractRepair, RepairPolicy } from '../utils/contract-repair.js' */
/** @import { ContractShape } from '../utils/contract-shape.js' */
/** @import { SessionResourceSource } from './session-resources.js' */
/** @import { OrderedConservationSource } from './conservation-series.js' */
/** @import { LinearComparisonResolution, LinearComparisonResolvedEdge } from './linear-comparisons.js' */
/** @import { LinearTrackSlot } from './linear-track-slot-model.js' */
/**
 * @import { FileLike, RawSessionResourceDescriptor, SessionResourceFileMetadata, SessionResourceTable, WebFileBindingRow,
 *   WebFileBindings } from './session-resource-backing.js'
 */
/** @import { CircularRequestRecordSetOptions } from './record-options.js' */
/** @import { CircularTrackSlot } from './circular-track-slot-model.js' */
/** @import { RecordDisplayDraft, RecordDisplayDraftIdentity, RecordDisplayRow } from './record-display-model.js' */
/** @import { ResourceBackedComparison } from './canonical-comparisons.js' */
/** @import { DepthFileSlots } from './depth-track-state.js' */
/** @import { CatalogStateView } from './feature-catalog.js' */
/** @import { ModeProfileName } from '../mode-profiles.js' */
/** @import { ExportRef } from './svg-serialization.js' */
/** @import { JsonValue } from './json-clone.js' */
/** @import { UnvalidatedFields } from './session-authority.js' */
/** @import { TableRepair } from './file-imports.js' */

// The drawing a request reads its settings and edits from: a `DrawingState`
// (state.js), Generate's drawing with the rules it admits, or a Gallery
// Session's (`buildCanonicalRequestState`). Project inputs, the mode and the
// generated artifacts still come from `state`.
/** @typedef {Record<string, any>} RequestDrawing */

/**
 * A file slot's value as `addFile` reads it: a `File` or a Session resource
 * view, or the descriptor of an encoded Session resource. Nothing has checked
 * which, so every field is read as a descriptor's unvalidated field.
 * @typedef {RawSessionResourceDescriptor} ResourceEntry
 */
/**
 * What the draft reads of a request record's selector, crop, and orientation
 * (a record of a Session is unvalidated until the request is promoted).
 * @typedef {object} RequestRecordView
 * @property {{ kind?: unknown, value?: unknown, index?: unknown } | null} [selector]
 * @property {{ selector?: { kind?: unknown, value?: unknown, index?: unknown } | null, start?: unknown, end?: unknown, reverseComplement?: unknown } | null} [region]
 * @property {{ reverseComplement?: unknown } | null} [presentation]
 */
/**
 * A record of the request before its display options are bound.
 * @typedef {Omit<CanonicalRecord, 'display'> & { display?: CanonicalRecord['display'] }} RequestDraftRecord
 */
/**
 * A comparison ring a Generate resolved (LOSAT), as the request reads it.
 * @typedef {object} ResolvedConservationEntry
 * @property {string} [name]
 * @property {string} [text]
 * @property {ResourceEntry | null} [fasta]
 * @property {string} [label]
 * @property {string} [color]
 * @property {string} [sourceKey]
 * @property {string} [seriesKey]
 */
/**
 * A file slot of a Linear record, as the request reads it: a `LinearSeq`
 * (state.js, which this state-free service cannot import, R13) whose files
 * are serialized, or the row a Session's request state is rebuilt with.
 * @typedef {object} RequestLinearSeq
 * @property {string} uid
 * @property {ResourceEntry | null} [gb]
 * @property {ResourceEntry | null} [gff]
 * @property {ResourceEntry | null} [fasta]
 * @property {DepthFileSlots} [depth]
 * @property {number | ''} [losat_gencode]
 * @property {string} [definition]
 * @property {string} [record_subtitle]
 * @property {string} [file_definition]
 * @property {string} [file_subtitle]
 * @property {string} [inferred_definition]
 * @property {string} [region_record_id]
 * @property {unknown} [region_start] A draft, or the value a Session stores: the projection validates it (R7).
 * @property {unknown} [region_end]
 * @property {boolean} [region_reverse]
 * @property {CanonicalRecord['cardinality']} [cardinality]
 */
/**
 * The input files and records a request reads (the `filesData` of the draft).
 * @typedef {object} RequestFilesData
 * @property {RequestLinearSeq[]} [linearSeqs]
 * @property {ResourceEntry | null} [c_gb]
 * @property {ResourceEntry | null} [c_gff]
 * @property {ResourceEntry | null} [c_fasta]
 * @property {DepthFileSlots | DepthFileSlots[]} [c_depth]
 * @property {FileLike[]} [c_conservation_blasts]
 * @property {string | null} [c_conservation_blasts_source]
 * @property {FileLike[]} [c_conservation_fastas]
 * @property {(ResourceEntry | null)[]} [c_conservation_sequence_sources]
 * @property {RequestCircularFileRecord[]} [circularRecords]
 * @property {{ id?: unknown, file?: ResourceEntry | null }[]} [linearComparisons]
 * @property {RequestCanonicalComparison[]} [linearCanonicalComparisons]
 */
/**
 * A Circular record the draft keeps with its canonical request pieces.
 * @typedef {object} RequestCircularFileRecord
 * @property {string} [sourceKind]
 * @property {ResourceEntry | null} [gb]
 * @property {ResourceEntry | null} [gff]
 * @property {ResourceEntry | null} [fasta]
 * @property {string} [recordKey]
 * @property {'exactly_one' | 'first' | 'all'} [cardinality]
 * @property {CanonicalRecord['selector']} [selector]
 * @property {CanonicalRecord['region']} [region]
 * @property {CanonicalRecord['presentation']} [presentation]
 * @property {CanonicalRecord['display']} [display]
 */
/**
 * A comparison a Session committed with its resource, as the draft holds it.
 * @typedef {ResourceBackedComparison & {
 *   file?: ResourceEntry | null, edgeKey?: unknown, mode?: string, settings?: ProteinSettingsBaseline
 * }} RequestCanonicalComparison
 */
/**
 * A comparison entry that names a resource, as `buildComparisons` writes it
 * before the request is validated: it copies `encoding` and `valueKind` of a
 * committed comparison as strings, so it is not yet one of the generated
 * comparison types (`contracts/request.generated.js`).
 * @typedef {object} RequestComparisonRef
 * @property {string} kind
 * @property {string} resourceId
 * @property {string} [encoding]
 * @property {string} [valueKind]
 * @property {number} [queryRecordIndex]
 * @property {number} [subjectRecordIndex]
 */
/**
 * What the label projection of a request reads of `state`: the displayed
 * Result's catalog features (`resultRenderedFeatures`) and its labels.
 * @typedef {CatalogStateView & {
 *   editableLabels?: ExportRef<{ featureId?: unknown, sourceText?: unknown }[]>
 * }} RequestLabelState
 */
/**
 * What the generated tables of a request read of `state` (state.js), or of
 * the state a Gallery Session's request is rebuilt from.
 * @typedef {RequestLabelState & {
 *   paletteDefinitions: ExportRef<Record<string, Record<string, string>>>,
 *   canonicalPublicationFiles?: Partial<Record<'d_color' | 't_color' | 'whitelist' | 'qualifier_priority', ResourceEntry | null>>
 * }} RequestTableState
 */
/**
 * What the request reads of `state` (state.js), or of the state a Gallery
 * Session's request is rebuilt from (`buildCanonicalRequestState`).
 * @typedef {RequestTableState & {
 *   mode: ExportRef<ModeProfileName>, cInputType: ExportRef<string>, lInputType: ExportRef<string>,
 *   circularRecordList: ExportRef<CircularRequestRecordSetOptions['records']>,
 *   losatExecution?: { threadsPerJob?: unknown }
 * }} RequestState
 */
/**
 * The settings of a committed generated protein comparison, as
 * `generatedProteinSettings` reads them for its baseline.
 * @typedef {{
 *   alignOrthogroupFeature?: unknown, losatpBin?: string, ncbiBlastpBin?: string | null,
 *   collinearityParams?: { kind?: string, parameters?: { mergeOrientation?: string } }
 * } & Record<string, unknown>} ProteinSettingsBaseline
 */
/**
 * @typedef {object} ResourceBuilder
 * @property {Record<string, RawSessionResourceDescriptor>} resources The descriptors as written; no reader has checked them.
 * @property {Record<string, string>} resourceOriginalNames
 * @property {(resourceId: string, kind: string | null, entry: ResourceEntry | null | undefined) => string} addFile
 * @property {(resourceId: string, kind: string, name: string, text: unknown) => string} addText
 * @property {(resourceId: string, kind: string | null, name: string, value: WeakKey) => string} addJson
 */

export { CANONICAL_REQUEST_SCHEMA };
// The canonical render request of the current schema (the writer format only;
// readers take unvalidated data, R14) is generated from the Python codec.
/**
 * @import {
 *   CanonicalAlignmentAnchorIdentity, CanonicalAlignmentRecordDecision, CanonicalCircularDiagramOptions,
 *   CanonicalComparison, CanonicalFileRef, CanonicalLinearRenderRequest, CanonicalRecord, CanonicalRenderRequest,
 *   CanonicalSimilarityAlignmentPlan, CanonicalTableRef
 * } from '../contracts/request.generated.js'
 */
/**
 * The table fields of a request's `diagramOptions` that the generated-table
 * writers fill; Circular and Linear requests declare them alike.
 * @typedef {Pick<CanonicalCircularDiagramOptions, 'colors' | 'featureVisibilityTableFile' | 'labelWhitelistTable'
 *   | 'labelWhitelistFile' | 'qualifierPriorityTable' | 'qualifierPriorityFile' | 'labelOverrideFile'>} RequestTableOptions
 */
/**
 * What the request writer returns and the session stores beside the request:
 * the encoded resources by id and the Web-only file metadata.
 * @typedef {object} CanonicalRenderEnvelope
 * @property {CanonicalRenderRequest} renderRequest
 * @property {Record<string, any>} resources
 * @property {Record<string, any>} webFiles
 */
const SUPPORTED_CANONICAL_REQUEST_SCHEMA_SET = new Set(SUPPORTED_CANONICAL_REQUEST_SCHEMAS);

/**
 * A layout object of a request its reader read: its keys are those of its
 * description (`readStoredRenderRequest`, or the writer), so no key is
 * checked here.
 * @param {unknown} value
 * @param {string} path
 * @returns {Record<string, unknown>}
 */
const requireCanonicalObject = (value, path) => {
  if (!isPlainObject(value)) throw new Error(`${path} must be an object.`);
  return value;
};

/**
 * @param {unknown} value
 * @param {string} path
 * @returns {string}
 */
const requireCanonicalText = (value, path) => {
  if (typeof value !== 'string' || !value.trim() || value.includes('\0')) {
    throw new Error(`${path} must be non-empty text without NUL.`);
  }
  return value.trim();
};

/**
 * @param {unknown} value
 * @param {string} path
 * @returns {CanonicalAlignmentAnchorIdentity}
 */
const canonicalAlignmentAnchor = (value, path) => {
  const anchor = requireCanonicalObject(value, path);
  if (anchor.sourceFeatureIndex !== null && (
    // A safe integer is a number.
    !Number.isSafeInteger(anchor.sourceFeatureIndex) || /** @type {number} */ (anchor.sourceFeatureIndex) < 0
  )) throw new Error(`${path}.sourceFeatureIndex must be a non-negative integer or null.`);
  if (anchor.stableFeatureSvgId !== null) {
    requireCanonicalText(anchor.stableFeatureSvgId, `${path}.stableFeatureSvgId`);
  }
  return {
    recordKey: requireCanonicalText(anchor.recordKey, `${path}.recordKey`),
    biologicalFeatureId: requireCanonicalText(
      anchor.biologicalFeatureId,
      `${path}.biologicalFeatureId`
    ),
    // Checked above: null or a non-negative safe integer, and null or text.
    sourceFeatureIndex: /** @type {number | null} */ (anchor.sourceFeatureIndex),
    stableFeatureSvgId: anchor.stableFeatureSvgId === null
      ? null
      : /** @type {string} */ (anchor.stableFeatureSvgId).trim()
  };
};

/**
 * @param {unknown} value
 * @param {Iterable<unknown>} recordKeys The displayed record keys (a legacy promotion passes them unchecked).
 * @param {string} path
 * @returns {CanonicalSimilarityAlignmentPlan | null}
 */
const canonicalSimilarityAlignment = (value, recordKeys, path) => {
  if (value === null) return null;
  const plan = requireCanonicalObject(value, path);
  if (plan.schema !== SIMILARITY_ALIGNMENT_PLAN_SCHEMA) throw new Error(`${path}.schema must be ${SIMILARITY_ALIGNMENT_PLAN_SCHEMA}.`);
  if (!Array.isArray(plan.records) || plan.records.length === 0) {
    throw new Error(`${path}.records must be a non-empty array.`);
  }
  const reference = canonicalAlignmentAnchor(plan.reference, `${path}.reference`);
  const records = /** @type {unknown[]} */ (plan.records).map((raw, index) => {
    const decisionPath = `${path}.records[${index}]`;
    const decision = requireCanonicalObject(raw, decisionPath);
    const recordKey = requireCanonicalText(decision.recordKey, `${decisionPath}.recordKey`);
    const anchor = decision.anchor === null
      ? null
      : canonicalAlignmentAnchor(decision.anchor, `${decisionPath}.anchor`);
    if (anchor && anchor.recordKey !== recordKey) {
      throw new Error(`${decisionPath}.anchor belongs to another record.`);
    }
    /** @type {Set<unknown>} */
    const alignedRationales = new Set([
      'user_selected', 'only_usable_candidate', 'unique_direct_rbh'
    ]);
    /** @type {Set<unknown>} */
    const skippedRationales = new Set([
      'skipped_by_user', 'skipped_no_candidate', 'skipped_unmappable'
    ]);
    const valid = decision.status === 'reference'
      ? anchor !== null && decision.rationale === 'reference'
      : decision.status === 'aligned'
        ? anchor !== null && alignedRationales.has(decision.rationale)
        : decision.status === 'skipped'
          ? anchor === null && skippedRationales.has(decision.rationale)
          : false;
    if (!valid) throw new Error(`${decisionPath} contains an invalid plan combination.`);
    // `valid` holds only for a status with one of its rationales.
    return {
      recordKey,
      status: /** @type {CanonicalAlignmentRecordDecision['status']} */ (decision.status),
      rationale: /** @type {CanonicalAlignmentRecordDecision['rationale']} */ (decision.rationale),
      anchor
    };
  });
  const decisionKeys = records.map((decision) => decision.recordKey);
  if (new Set(decisionKeys).size !== decisionKeys.length) {
    throw new Error(`${path}.records contains duplicate record keys.`);
  }
  const references = records.filter((decision) => decision.status === 'reference');
  if (references.length !== 1 || references[0].recordKey !== reference.recordKey ||
      JSON.stringify(references[0].anchor) !== JSON.stringify(reference)) {
    throw new Error(`${path}.reference must match exactly one reference decision.`);
  }
  const expected = new Set(recordKeys);
  /** @type {Set<unknown>} */
  const actual = new Set(decisionKeys);
  if (expected.size !== actual.size || [...expected].some((key) => !actual.has(key))) {
    throw new Error(`${path}.records does not cover the displayed record keys.`);
  }
  return {
    schema: SIMILARITY_ALIGNMENT_PLAN_SCHEMA,
    groupId: requireCanonicalText(plan.groupId, `${path}.groupId`),
    reference,
    records
  };
};

/**
 * @param {unknown} value
 * @param {Iterable<string>} recordKeys
 * @param {string} path
 * @param {{ requireCoverage?: boolean }} [options]
 */
const canonicalRecordTranslations = (value, recordKeys, path, { requireCoverage = false } = {}) => {
  if (!Array.isArray(value)) throw new Error(`${path} must be an array.`);
  const translations = value.map((raw, index) => {
    const itemPath = `${path}[${index}]`;
    const item = requireCanonicalObject(raw, itemPath);
    if (!Number.isFinite(item.x) || !Number.isFinite(item.y)) {
      throw new Error(`${itemPath} x and y must be finite numbers.`);
    }
    return {
      recordKey: requireCanonicalText(item.recordKey, `${itemPath}.recordKey`),
      x: Number(item.x),
      y: Number(item.y)
    };
  });
  const keys = translations.map((item) => item.recordKey);
  if (new Set(keys).size !== keys.length) {
    throw new Error(`${path} contains duplicate record keys.`);
  }
  if (requireCoverage || translations.length > 0) {
    const expected = new Set(recordKeys);
    const actual = new Set(keys);
    if (expected.size !== actual.size || [...expected].some((key) => !actual.has(key))) {
      throw new Error(`${path} does not cover the displayed record keys.`);
    }
  }
  return translations;
};

// The Linear Accession and Length visibility a committed request resolved,
// by the `adv` key that selects it; null when the request is not Linear.
/**
 * @param {Record<string, any> | null} renderRequest
 * @returns {Record<string, boolean | null>}
 */
export const linearDefinitionVisibilityOf = (renderRequest) => {
  const overrides = renderRequest?.mode === 'linear'
    ? renderRequest.diagramOptions?.configOverrides || {}
    : {};
  /** @param {string} path @returns {boolean | null} */
  const read = (path) => (typeof overrides[path] === 'boolean' ? overrides[path] : null);
  return {
    linear_accession_visibility: read(CONFIG_OVERRIDE_PATHS.linearDefinitionShowAccession),
    linear_length_visibility: read(CONFIG_OVERRIDE_PATHS.linearDefinitionShowLength)
  };
};

// E1: the similarity alignment plan and record translations of a committed
// Linear request, admitted with the rules the request projection applies;
// Session Load takes them from the Linear Result set wherever it sits.
/** @param {Record<string, any>} renderRequest */
export const canonicalLinearRecordLayout = (renderRequest) => {
  const layout = renderRequest?.layout || {};
  if (renderRequest?.mode !== 'linear' || Object.keys(layout).length === 0) {
    return { similarityAlignment: null, recordTranslations: [] };
  }
  const recordKeys = (Array.isArray(renderRequest.records) ? renderRequest.records : [])
    .map((/** @type {Record<string, any>} */ record, /** @type {number} */ index) => requireCanonicalText(
      record?.recordKey, `renderRequest.records[${index}].recordKey`
    ));
  const similarityAlignment = canonicalSimilarityAlignment(
    layout.similarityAlignment ?? null, recordKeys, 'renderRequest.layout.similarityAlignment'
  );
  return {
    similarityAlignment,
    recordTranslations: canonicalRecordTranslations(
      layout.recordTranslations ?? [], recordKeys, 'renderRequest.layout.recordTranslations',
      { requireCoverage: similarityAlignment !== null }
    )
  };
};

/** @param {UnvalidatedFields<'isCircular' | 'startCoordinate'>} raw A record's `display`, checked here. */
const canonicalRecordDisplay = (raw) => {
  if (!raw || Object.keys(raw).sort().join(',') !== 'isCircular,startCoordinate'
    || (raw.isCircular !== null && typeof raw.isCircular !== 'boolean')
    // A safe integer is a number.
    || (raw.startCoordinate !== null && (!Number.isSafeInteger(raw.startCoordinate) || /** @type {number} */ (raw.startCoordinate) < 1))
    || (raw.isCircular === false && raw.startCoordinate !== null)) {
    throw new Error('Record display requires nullable isCircular and a positive integer or null startCoordinate.');
  }
  return { ...raw };
};

// Fresh canonical requests use only these typed GbdrawConfig leaf paths.
// The snake_case aliases derived below exist solely in the persisted-session
// reader; they are never emitted by buildConfigOverrides().
const CONFIG_OVERRIDE_PATHS = Object.freeze({
  arrowHeadLengthRatio: 'objects.features.arrow_geometry.head_length_ratio',
  arrowShaftWidthRatio: 'objects.features.arrow_geometry.shaft_width_ratio',
  blockStrokeColor: 'objects.features.block_stroke_color',
  circularAxisStrokeColor: 'objects.axis.circular.stroke_color',
  linearAxisStrokeColor: 'objects.axis.linear.stroke_color',
  lineStrokeColor: 'objects.features.line_stroke_color',
  circularDefinitionFontSize: 'objects.definition.circular.font_size', circularDefinitionInterval: 'objects.definition.circular.interval',
  plotTitleFontSize: 'objects.definition.circular.plot_title_font_size',
  showGc: 'canvas.show_gc',
  showSkew: 'canvas.show_skew',
  showDepth: 'canvas.show_depth',
  strandedness: 'canvas.strandedness',
  resolveOverlaps: 'canvas.resolve_overlaps',
  featureOverlapToleranceBp: 'canvas.feature_overlap_tolerance_bp',
  trackType: 'canvas.circular.track_type',
  alignCenter: 'canvas.linear.align_center',
  keepDefinitionLeftAligned: 'canvas.linear.keep_definition_left_aligned',
  linearTrackLayout: 'canvas.linear.track_layout',
  linearTrackAxisGap: 'canvas.linear.track_axis_gap',
  linearRulerOnAxis: 'canvas.linear.ruler_on_axis',
  comparisonHeight: 'canvas.linear.comparison_height',
  gcHeight: 'canvas.linear.default_gc_height',
  depthHeight: 'canvas.linear.depth_height',
  normalizeLength: 'canvas.linear.normalize_length',
  labelRendering: 'labels.rendering',
  circularLabelSpacing: 'labels.spacing.circular',
  circularLabelPlacement: 'labels.circular.placement',
  linearLabelSpacing: 'labels.spacing.linear',
  labelPlacement: 'labels.linear.placement',
  labelRotation: 'labels.linear.rotation',
  labelBlacklist: 'labels.filtering.blacklist_keywords',
  linearDefinitionShowReplicon: 'objects.definition.linear.show_replicon',
  linearDefinitionShowAccession: 'objects.definition.linear.show_accession',
  linearDefinitionShowLength: 'objects.definition.linear.show_length',
  gcContentMode: 'objects.gc_content.mode',
  gcContentMinPercent: 'objects.gc_content.min_percent',
  gcContentMaxPercent: 'objects.gc_content.max_percent',
  gcContentShowAxis: 'objects.gc_content.show_axis',
  gcContentShowTicks: 'objects.gc_content.show_ticks',
  gcContentLargeTickInterval: 'objects.gc_content.large_tick_interval',
  gcContentSmallTickInterval: 'objects.gc_content.small_tick_interval',
  gcContentTickFontSize: 'objects.gc_content.tick_font_size',
  depthColor: 'objects.depth.fill_color',
  depthMin: 'objects.depth.min_depth',
  depthMax: 'objects.depth.max_depth',
  depthNormalize: 'objects.depth.normalize',
  depthShowAxis: 'objects.depth.show_axis',
  depthShowTicks: 'objects.depth.show_ticks',
  depthLargeTickInterval: 'objects.depth.large_tick_interval',
  depthSmallTickInterval: 'objects.depth.small_tick_interval',
  depthTickFontSize: 'objects.depth.tick_font_size',
  depthShareAxis: 'objects.depth.share_axis',
  showScale: 'objects.scale.show',
  scaleStyle: 'objects.scale.style',
  scaleStrokeColor: 'objects.scale.stroke_color',
  scaleLabelColor: 'objects.scale.label_color',
  scaleStrokeWidth: 'objects.scale.stroke_width',
  scaleInterval: 'objects.scale.interval',
  tickLabelFontSize: 'objects.ticks.tick_labels.font_size',
  outerLabelXRadiusOffset: 'labels.unified_adjustment.outer_labels.x_radius_offset',
  outerLabelYRadiusOffset: 'labels.unified_adjustment.outer_labels.y_radius_offset',
  innerLabelXRadiusOffset: 'labels.unified_adjustment.inner_labels.x_radius_offset',
  innerLabelYRadiusOffset: 'labels.unified_adjustment.inner_labels.y_radius_offset',
  pairwiseMatchStyle: 'objects.blast_match.style',
  fontFamily: 'objects.text.font_family'
});

const SHARED_LENGTH_CONFIG_OVERRIDE_PATHS = Object.freeze({
  blockStrokeWidth: 'objects.features.block_stroke_width',
  circularAxisStrokeWidth: 'objects.axis.circular.stroke_width',
  linearAxisStrokeWidth: 'objects.axis.linear.stroke_width',
  lineStrokeWidth: 'objects.features.line_stroke_width',
  linearDefinitionFontSize: 'objects.definition.linear.font_size',
  defaultCdsHeight: 'canvas.linear.default_cds_height',
  legendBoxSize: 'objects.legends.color_rect_size',
  legendFontSize: 'objects.legends.font_size',
  scaleFontSize: 'objects.scale.font_size',
  rulerLabelFontSize: 'objects.scale.ruler_label_font_size'
});

const MODE_LABEL_SCOPE_PATHS = Object.freeze({
  circular: 'labels.circular.scope',
  linear: 'labels.linear.scope'
});

/** @type {Readonly<Record<string, string>>} */
const CIRCULAR_LABEL_SCOPES = Object.freeze({ none: 'none', out: 'outer', both: 'both' });
// The Labels choice a committed Circular label scope shows; any other scope shows None.
/** @type {Readonly<Record<string, string>>} */
const CIRCULAR_LABEL_MODES = Object.freeze({ outer: 'out', both: 'both' });
/** The Labels choice a Circular label scope shows. @param {unknown} scope */
const circularLabelsMode = (scope) => (
  Object.hasOwn(CIRCULAR_LABEL_MODES, String(scope)) ? CIRCULAR_LABEL_MODES[String(scope)] : 'none'
);

/**
 * @param {{ labels_mode?: unknown, show_labels_linear?: string }} form
 * @param {boolean} circular
 */
const labelScopeOverride = (form, circular) => (circular
  ? (Object.hasOwn(CIRCULAR_LABEL_SCOPES, String(form.labels_mode))
      ? CIRCULAR_LABEL_SCOPES[String(form.labels_mode)]
      : 'none')
  : form.show_labels_linear);

/** @param {RequestDrawing} drawing */
const labelBlacklistOverride = (drawing) => (drawing.filterMode.value === 'Blacklist'
  ? String(drawing.manualBlacklist.value || '').split(/[,\n]/)
    .map((keyword) => keyword.trim()).filter(Boolean)
  : []);

const LINEAR_DEFINITION_STYLE_PATHS = Object.freeze(
  Object.fromEntries(
    ['name', 'subtitle', 'replicon', 'accession', 'length'].map((kind) => [
      kind,
      `objects.definition.linear.line_styles.${kind}`
    ])
  )
);
const LINEAR_DEFINITION_STYLE_FIELDS = Object.freeze(['font_size', 'font_weight', 'fill']);
const LINEAR_DEFINITION_WEIGHT_PATHS = new Set(
  Object.values(LINEAR_DEFINITION_STYLE_PATHS).map((prefix) => `${prefix}.font_weight`)
);

const CIRCULAR_ONLY_GUI_CONFIG_OVERRIDE_PATHS = new Set([
  CONFIG_OVERRIDE_PATHS.circularAxisStrokeColor,
  CONFIG_OVERRIDE_PATHS.circularDefinitionFontSize,
  CONFIG_OVERRIDE_PATHS.circularDefinitionInterval,
  CONFIG_OVERRIDE_PATHS.plotTitleFontSize,
  CONFIG_OVERRIDE_PATHS.circularLabelSpacing,
  CONFIG_OVERRIDE_PATHS.circularLabelPlacement,
  CONFIG_OVERRIDE_PATHS.trackType,
  CONFIG_OVERRIDE_PATHS.tickLabelFontSize,
  CONFIG_OVERRIDE_PATHS.outerLabelXRadiusOffset,
  CONFIG_OVERRIDE_PATHS.outerLabelYRadiusOffset,
  CONFIG_OVERRIDE_PATHS.innerLabelXRadiusOffset,
  CONFIG_OVERRIDE_PATHS.innerLabelYRadiusOffset
]);
const LINEAR_ONLY_GUI_CONFIG_OVERRIDE_PATHS = new Set([
  CONFIG_OVERRIDE_PATHS.linearAxisStrokeColor,
  CONFIG_OVERRIDE_PATHS.linearDefinitionShowReplicon,
  CONFIG_OVERRIDE_PATHS.linearDefinitionShowAccession,
  CONFIG_OVERRIDE_PATHS.linearDefinitionShowLength,
  CONFIG_OVERRIDE_PATHS.linearLabelSpacing,
  CONFIG_OVERRIDE_PATHS.labelPlacement,
  CONFIG_OVERRIDE_PATHS.labelRotation,
  CONFIG_OVERRIDE_PATHS.alignCenter,
  CONFIG_OVERRIDE_PATHS.keepDefinitionLeftAligned,
  CONFIG_OVERRIDE_PATHS.linearTrackLayout,
  CONFIG_OVERRIDE_PATHS.linearTrackAxisGap,
  CONFIG_OVERRIDE_PATHS.linearRulerOnAxis,
  CONFIG_OVERRIDE_PATHS.comparisonHeight,
  CONFIG_OVERRIDE_PATHS.pairwiseMatchStyle,
  CONFIG_OVERRIDE_PATHS.gcHeight,
  CONFIG_OVERRIDE_PATHS.depthHeight,
  CONFIG_OVERRIDE_PATHS.scaleStyle,
  CONFIG_OVERRIDE_PATHS.scaleStrokeColor,
  CONFIG_OVERRIDE_PATHS.scaleLabelColor,
  CONFIG_OVERRIDE_PATHS.scaleStrokeWidth,
  CONFIG_OVERRIDE_PATHS.normalizeLength
]);

/** @param {ModeProfileName} mode */
export const managedConfigOverridePathsForMode = (mode) => {
  if (!['circular', 'linear'].includes(mode)) {
    throw new Error(`Unsupported config override mode: ${String(mode)}.`);
  }
  /** @type {Set<string>} */
  const excluded = mode === 'circular'
    ? LINEAR_ONLY_GUI_CONFIG_OVERRIDE_PATHS
    : CIRCULAR_ONLY_GUI_CONFIG_OVERRIDE_PATHS;
  /** @type {Set<string>} */
  const paths = new Set(
    Object.values(CONFIG_OVERRIDE_PATHS).filter((path) => !excluded.has(path))
  );
  paths.add(MODE_LABEL_SCOPE_PATHS[mode]);
  for (const path of [
    SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.blockStrokeWidth,
    SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.lineStrokeWidth,
    SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.legendBoxSize,
    SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.legendFontSize,
    ...(mode === 'circular'
      ? [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.circularAxisStrokeWidth, 'labels.font_size']
      : [
          SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.linearAxisStrokeWidth,
          SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.linearDefinitionFontSize,
          SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.defaultCdsHeight,
          SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.scaleFontSize,
          SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.rulerLabelFontSize,
          'labels.font_size.linear'
        ])
  ]) {
    paths.add(`${path}.short`);
    paths.add(`${path}.long`);
  }
  if (mode === 'linear') {
    Object.values(LINEAR_DEFINITION_STYLE_PATHS).forEach((prefix) => {
      LINEAR_DEFINITION_STYLE_FIELDS.forEach((field) => paths.add(`${prefix}.${field}`));
    });
  }
  return Object.freeze([...paths].sort());
};

// The draft (`drawing.adv`) field of a flat config key whose name differs
// from it; every other projected flat key is its own draft field.
/** @type {Readonly<Record<string, string>>} */
const DRAFT_FIELD_OF_FLAT_CONFIG_KEY = Object.freeze({
  outer_label_x_radius_offset: 'outer_label_x_offset',
  outer_label_y_radius_offset: 'outer_label_y_offset',
  inner_label_x_radius_offset: 'inner_label_x_offset',
  inner_label_y_radius_offset: 'inner_label_y_offset'
});

/**
 * The draft fields that the request mapping projects to `values` (config
 * paths and their values), as a label fit row writes them: the reverse of the
 * request's config overrides. The label scope goes to the form's label mode;
 * the rest to the drawing's `adv`. A path the mode does not manage, or one
 * whose draft field the drawing does not hold, is an `UNKNOWN_CONFIG_PATH`
 * error.
 * @param {Readonly<Record<string, number | string | boolean>>} values
 * @param {ModeProfileName} mode
 * @param {Readonly<Record<string, unknown>>} adv The drawing's draft, whose fields the result names.
 * @returns {{ form: Record<string, string>, adv: Record<string, number | string | boolean> }}
 */
export const draftFieldsForConfigValues = (values, mode, adv) => {
  const managed = new Set(managedConfigOverridePathsForMode(mode));
  const unknown = () => diagnosticError('INPUT_INVALID', { field: 'configOverrides', reason: 'UNKNOWN_CONFIG_PATH' });
  if (Object.keys(values).some((path) => !managed.has(path))) throw unknown();
  // A copy of the draft fields, not the drawing's own objects.
  /** @type {Record<string, string>} */
  const formFields = {};
  /** @type {Record<string, number | string | boolean>} */
  const advFields = {};
  for (const [key, value] of Object.entries(projectCanonicalConfigOverrides(values, mode))) {
    if (key === 'label_scope') {
      if (typeof value !== 'string') throw unknown();
      formFields[mode === 'circular' ? 'labels_mode' : 'show_labels_linear'] = mode === 'circular'
        ? circularLabelsMode(value) : value;
      continue;
    }
    const field = DRAFT_FIELD_OF_FLAT_CONFIG_KEY[key] || key;
    if (!Object.hasOwn(adv, field) || !['number', 'string', 'boolean'].includes(typeof value)) throw unknown();
    advFields[field] = value;
  }
  if (Object.keys(formFields).length + Object.keys(advFields).length === 0) throw unknown();
  return { form: formFields, adv: advFields };
};

// The font leaf that a Session kept as a GUI-unmanaged preserved override
// (PD-OI-009) before the GUI had a Font setting. A reader moves it into
// `adv.font_family`, so a preserved override is never a path the GUI manages.
const PROMOTED_FONT_FAMILY_PATH = CONFIG_OVERRIDE_PATHS.fontFamily;

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isPlainRecord = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

/**
 * A stored draft (`config` of a Session 27-44, or of a Session 46 slice) with
 * the preserved font override moved into the Font setting. A font the draft
 * already holds wins; the preserved leaf is dropped either way.
 * @template T
 * @param {T} config
 * @returns {T}
 */
export const adoptPromotedPreservedOverrides = (config) => {
  if (!isPlainRecord(config) || !isPlainRecord(config.unmanagedConfigOverrides)
    || !Object.hasOwn(config.unmanagedConfigOverrides, PROMOTED_FONT_FAMILY_PATH)) return config;
  const { [PROMOTED_FONT_FAMILY_PATH]: fontFamily, ...preserved } = config.unmanagedConfigOverrides;
  const adv = isPlainRecord(config.adv) ? { ...config.adv } : {};
  adv.font_family ??= fontFamily;
  return /** @type {T} */ ({ ...config, adv, unmanagedConfigOverrides: preserved });
};

/** @param {string} semanticName */
const legacyFlatConfigKey = (semanticName) => (
  semanticName
    .replace(/([A-Z]+)([A-Z][a-z])/g, '$1_$2')
    .replace(/([a-z0-9])([A-Z])/g, '$1_$2')
    .toLowerCase()
);

// FL-10: the typed Output Prefix goes to Python as typed (Python strips it as
// this does); Python rejects a value that is not one portable file name
// (INPUT_INVALID, field Output Prefix). A name derived from a record ID is
// made safe here, by the browser's download-name rule (utils/download-names.js).
/** @param {unknown} value */
const explicitOutputPrefix = (value) => {
  const raw = String(value ?? '').trim();
  return raw || null;
};

/**
 * @param {{ record_id?: unknown, recordId?: unknown } | null | undefined} record
 * @param {number} index
 */
const circularRecordId = (record, index) => (
  String(record?.record_id ?? record?.recordId ?? '').trim() || `Record_${index + 1}`
);

/**
 * @param {({ record_id?: unknown, recordId?: unknown } | null | undefined)[]} records
 * @param {string | null} explicitPrefix
 * @returns {string[]}
 */
const resolveCircularBatchPrefixes = (records, explicitPrefix) => {
  if (explicitPrefix !== null) {
    if (records.length === 1) return [explicitPrefix];
    return records.map((_, index) => `${explicitPrefix}_${index + 1}`);
  }
  /** @type {string[]} */
  const prefixes = [];
  const used = new Set();
  records.forEach((record, index) => {
    const base = downloadSafeName(circularRecordId(record, index), 'out');
    let candidate = base;
    let suffix = 2;
    while (used.has(candidate)) {
      candidate = `${base}_${suffix}`;
      suffix += 1;
    }
    used.add(candidate);
    prefixes.push(candidate);
  });
  return prefixes;
};

/** @param {string | null} prefix */
const renderOutputPayload = (prefix) => ({
  prefix,
  formats: ['svg'],
  overwrite: false,
  interactiveMetadataPolicy: 'auto'
});

/** @param {Record<string, any>} thresholds */
const projectComparisonThresholds = ({ evalue, bitscore, identity, alignmentLength }) => ({
  evalue: Number(evalue), bitscore, identity, alignmentLength
});
// An integer setting that the current request does not use for rendering
// (a mode-inactive protein parameter or a LOSAT execution choice) keeps its
// documented fallback; projected render values use projectOptionalNumber.
/**
 * @template F
 * @param {unknown} value
 * @param {F} fallback
 * @param {number} minimum
 * @returns {number | F}
 */
const integerSettingOr = (value, fallback, minimum) => {
  const { value: numeric } = classifyOptionalNumber(value);
  return typeof numeric === 'number' && Number.isInteger(numeric) && numeric >= minimum ? numeric : fallback;
};

// A record label or subtitle is saved resolved. Sessions written since file
// defaults were introduced also record the per-record value, which says directly
// whether the record inherited the file default or overrode it with the same
// text. Older sessions carry no such value, so fall back to comparing them.
/** @param {{ savedOverride?: unknown, fileDefault?: string | null, resolved: string }} input */
const resolveSavedRecordOverride = ({ savedOverride, fileDefault, resolved }) => {
  if (savedOverride !== undefined && savedOverride !== null) return String(savedOverride);
  if (!fileDefault) return resolved;
  return resolved === fileDefault ? '' : resolved;
};

// A Depth series that no record gives a file (TK-06). The one producer of
// this failure for the request builder and the Session projections.
/** @param {number} seriesIndex */
const depthSeriesWithoutSource = (seriesIndex) => (
  diagnosticError('DEPTH_INVALID', { seriesIndex, reason: 'DEPTH_SERIES_SOURCE' })
);

/**
 * @param {unknown[]} depthRows
 * @param {number} logicalTrackCount
 */
const validateProjectedDepthSources = (depthRows, logicalTrackCount) => {
  for (let trackIndex = 0; trackIndex < logicalTrackCount; trackIndex += 1) {
    const hasSource = depthRows.some((row) => (
      Array.isArray(row) && Boolean(row[trackIndex]?.resourceId)
    ));
    if (!hasSource) throw depthSeriesWithoutSource(trackIndex);
  }
};

/**
 * @param {string} resourceId
 * @param {unknown} name
 */
const normalizeResourceName = (resourceId, name) => {
  // split() always returns at least one element, so pop() is a string.
  const basename = /** @type {string} */ (String(name || 'resource.dat').replace(/\\/g, '/').split('/').pop());
  const safe = basename.replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^[._]+|[._]+$/g, '');
  const prefix = `${resourceId}-`;
  let leaf = safe || 'resource.dat';
  while (leaf.startsWith(prefix)) {
    leaf = leaf.slice(prefix.length);
  }
  return `${prefix}${leaf || 'resource.dat'}`;
};

/** @param {unknown} name */
const normalizeOriginalResourceName = (name) => {
  // split() always returns at least one element, so pop() is a string.
  const basename = /** @type {string} */ (String(name || '')
    .replace(/\\/g, '/')
    .split('/')
    .pop())
    .replace(/[\u0000-\u001f\u007f]/g, '')
    .trim();
  if (!basename || basename === '.' || basename === '..') return '';
  return basename.slice(0, 1024);
};

/**
 * The text or JSON value of each descriptor the request writer generated, keyed by that descriptor.
 * @type {WeakMap<object, string | { bindings: unknown[] }>}
 */
const generatedResourceValues = new WeakMap();
const canonicalTypedResourceBackings = new WeakMap();
// The one lookup of a typed result's Python bytes. A derived-cache hit hands the result back through the
// reactive cache as a Vue proxy, so the binding is keyed and read by the object behind it (OV-1040).
/** @param {unknown} value */
const canonicalTypedResourceBacking = (value) => (
  value && typeof value === 'object' ? canonicalTypedResourceBackings.get(toRaw(value)) : undefined
);
/**
 * Whether the analysis helper bound this typed result to its Python bytes. A result read from a Session or
 * cloned is not bound, so the analysis computes it again instead of encoding it in JS.
 * @param {unknown} value
 */
export const isCanonicalTypedResourceBound = (value) => canonicalTypedResourceBacking(value) !== undefined;

// The analysis helper supplies the bytes already emitted by the Python typed codec.
/**
 * @param {{ kind?: unknown, schema?: unknown } | null | undefined} value The decoded result, checked here.
 * @param {RawSessionResourceDescriptor | null | undefined} descriptor Checked here.
 */
export const bindCanonicalTypedResource = (value, descriptor) => {
  const expectedKind = value?.kind === "result" ? "collinearity-result" : "orthogroup-result";
  if (!value || value.schema !== 3 || !(/** @type {readonly unknown[]} */ (["result", "orthogroupResult"])).includes(value.kind)
    || descriptor?.kind !== expectedKind || descriptor.encoding !== "base64"
    || descriptor.type !== "application/json" || typeof descriptor.data !== "string"
    // `isSafeInteger` holds, so `size` is a number.
    || !Number.isSafeInteger(descriptor.size) || /** @type {number} */ (descriptor.size) < 0) {
    throw new Error("Analysis canonical resource backing is invalid.");
  }
  canonicalTypedResourceBackings.set(toRaw(value), descriptor);
};
/**
 * @param {{ encode?: boolean }} [options]
 * @returns {ResourceBuilder}
 */
const createResourceBuilder = ({ encode = true } = {}) => {
  // A resource id can be one a Session names (`publicationFileRef`), so the
  // table is read by own key (OV-564).
  /** @type {Record<string, RawSessionResourceDescriptor>} */
  const resources = {};
  /** @type {Record<string, string>} */
  const resourceOriginalNames = {};
  const fileResourceIds = new Map();

  /** @type {ResourceBuilder['addFile']} */
  const addFile = (resourceId, kind, entry) => {
    if (!entry || typeof entry !== 'object' || Array.isArray(entry)) {
      throw new Error(`Canonical resource ${resourceId} is missing.`);
    }
    if (Object.hasOwn(resources, resourceId)) return resourceId;
    const adoptedSource = /** @type {SessionResourceSource | null} */ (
      getSessionResourceSource(entry)
    );
    const encodedEntry = adoptedSource?.descriptor || entry;
    const effectiveKind = String(encodedEntry.kind || kind);
    const kindResourceIds = fileResourceIds.get(effectiveKind) || new WeakMap();
    const owner = getResourcePayloadOwner(entry);
    const identity = /** @type {SessionResourceSource | null} */ (
      getSessionResourceSource(owner)
    )?.descriptor || owner;
    const existingResourceId = kindResourceIds.get(identity);
    if (existingResourceId) return existingResourceId;
    /** @type {RawSessionResourceDescriptor} */
    const descriptor = {
      kind: effectiveKind,
      name: normalizeResourceName(resourceId, entry.name),
      type: String(entry.type || 'application/octet-stream'),
      size: Number(entry.size) || 0,
      lastModified: Number(entry.lastModified) || 0,
      encoding: encodedEntry.encoding || 'base64',
      data: encodedEntry.data
    };
    if (typeof encodedEntry.checksum === 'string' && encodedEntry.checksum.trim()) {
      descriptor.checksum = encodedEntry.checksum;
    }
    resources[resourceId] = setResourcePayloadOwner(
      descriptor,
      getResourcePayloadOwner(entry)
    );
    kindResourceIds.set(identity, resourceId);
    fileResourceIds.set(effectiveKind, kindResourceIds);
    const originalName = normalizeOriginalResourceName(entry.name);
    if (originalName) resourceOriginalNames[resourceId] = originalName;
    return resourceId;
  };

  /** @type {ResourceBuilder['addText']} */
  const addText = (resourceId, kind, name, text) => {
    if (Object.hasOwn(resources, resourceId)) return resourceId;
    const normalized = String(text || '');
    if (!encode) {
      resources[resourceId] = { kind };
      generatedResourceValues.set(resources[resourceId], normalized);
      return resourceId;
    }
    const bytes = textToBytes(normalized);
    resources[resourceId] = {
      kind,
      name: normalizeResourceName(resourceId, name),
      type: 'text/tab-separated-values',
      size: bytes.byteLength,
      lastModified: 0,
      encoding: 'base64',
      data: bytesToBase64(bytes)
    };
    generatedResourceValues.set(resources[resourceId], normalized);
    return resourceId;
  };

  /** @type {ResourceBuilder['addJson']} */
  const addJson = (resourceId, kind, name, value) => {
    if (Object.hasOwn(resources, resourceId)) return resourceId;
    if (!encode) {
      resources[resourceId] = { kind };
      generatedResourceValues.set(resources[resourceId], { bindings: [value] });
      return resourceId;
    }
    const backing = canonicalTypedResourceBacking(value);
    if (!backing) throw diagnosticError('RESOURCE_INVALID', { field: resourceId, reason: 'UNBOUND_CANONICAL_BYTES' }, { stage: 'resource-staging' });
    if (backing.kind !== kind) throw diagnosticError('RESOURCE_INVALID', { field: resourceId, reason: 'KIND_MISMATCH' }, { stage: 'resource-staging' });
    resources[resourceId] = { ...backing, name: normalizeResourceName(resourceId, name) };
    generatedResourceValues.set(resources[resourceId], { bindings: [value] });
    return resourceId;
  };

  return {
    resources,
    resourceOriginalNames,
    addFile,
    addText,
    addJson
  };
};

/**
 * @param {string} resourceId
 * @returns {CanonicalFileRef}
 */
const fileRef = (resourceId) => ({ resourceId, representation: 'file' });
// The table inputs whose committed file a rebuilt request keeps (`publicationFileRef`).
const REQUEST_TABLE_FILE_SLOTS = Object.freeze(['d_color', 't_color', 'whitelist', 'qualifier_priority']);
/**
 * @param {ResourceBuilder} resources
 * @param {Record<string, ResourceEntry | null | undefined>} files
 * @param {string} key
 * @param {string} fallbackId
 * @returns {CanonicalFileRef | CanonicalTableRef | null}
 */
const publicationFileRef = (resources, files, key, fallbackId) => {
  if (!files[key]) return null;
  const source = /** @type {SessionResourceSource | null} */ (getSessionResourceSource(files[key]));
  return {
    resourceId: resources.addFile(source?.resourceId || fallbackId, fallbackId, files[key]),
    representation: source?.descriptor?.kind === 'canonical-tsv' ? 'canonicalTsv' : 'file'
  };
};

/**
 * @param {unknown} rawValue
 * @returns {CanonicalRecord['selector']}
 */
const selectorPayload = (rawValue) => {
  const raw = String(rawValue || '').trim();
  if (!raw) return null;
  const indexMatch = raw.match(/^#(\d+)$/);
  if (indexMatch) {
    const index = Number(indexMatch[1]) - 1;
    return Number.isInteger(index) && index >= 0 ? { kind: 'recordIndex', index } : null;
  }
  return { kind: 'recordId', value: raw };
};

/** @param {{ label?: unknown, subtitle?: unknown, reverseComplement?: unknown, gridRow?: number | null, gridColumn?: number | null }} [options] */
const presentationPayload = ({
  label = null,
  subtitle = null,
  reverseComplement = false,
  gridRow = null,
  gridColumn = null
} = {}) => ({
  label: String(label || '').trim() || null,
  subtitle: String(subtitle || '').trim() || null,
  reverseComplement: Boolean(reverseComplement),
  gridRow,
  gridColumn
});

/**
 * @param {{ circular_record_label?: unknown, circular_record_subtitle?: unknown, circular_reverse?: unknown }} form
 * @param {{ hasRegion?: boolean }} [options]
 */
const circularPresentationPayload = (form, { hasRegion = false } = {}) => (
  presentationPayload({
    label: form.circular_record_label,
    subtitle: form.circular_record_subtitle,
    reverseComplement: !hasRegion && form.circular_reverse
  })
);

/**
 * @param {{ circular_region_start?: unknown, circular_region_end?: unknown, circular_reverse?: unknown }} form
 * @param {{ recordLength?: unknown, record_length?: unknown, value?: unknown } | null} record
 * @returns {CanonicalRecord['region']}
 */
const circularRegionPayload = (form, record) => {
  const rawStart = form.circular_region_start;
  const rawEnd = form.circular_region_end;
  const hasStart = rawStart !== null && rawStart !== undefined && String(rawStart).trim() !== '';
  const hasEnd = rawEnd !== null && rawEnd !== undefined && String(rawEnd).trim() !== '';
  if (!hasStart && !hasEnd) return null;
  if (hasStart !== hasEnd) {
    throw new Error('Circular region requires both Start and End coordinates.');
  }
  const start = Number(rawStart);
  const end = Number(rawEnd);
  if (!Number.isInteger(start) || !Number.isInteger(end) || start < 1 || end < 1) {
    throw new Error('Circular region Start and End must be positive integers.');
  }
  if (start > end) {
    throw new Error('Circular region Start must not exceed End. Use Reverse complement to change display orientation.');
  }
  const recordLength = Number(record?.recordLength ?? record?.record_length);
  if (Number.isInteger(recordLength) && recordLength > 0 && end > recordLength) {
    throw new Error(`Circular region End (${end}) exceeds the selected record length (${recordLength}).`);
  }
  return {
    selector: selectorPayload(record?.value),
    start,
    end,
    reverseComplement: Boolean(form.circular_reverse)
  };
};

/** @param {{ region_start?: unknown, region_end?: unknown, region_record_id?: unknown, region_reverse?: unknown } | null | undefined} seq */
const linearRegionPayload = (seq) => {
  // projectOptionalNumber yields null or a finite number: the only arm without
  // a number-or-null value is 'invalid', which it throws on.
  const start = /** @type {number | null} */ (projectOptionalNumber(seq?.region_start, { field: 'region' }));
  const end = /** @type {number | null} */ (projectOptionalNumber(seq?.region_end, { field: 'region' }));
  if (start === null && end === null) return null;
  if (start === null || end === null) {
    throw diagnosticError('REGION_INVALID', { field: 'region', reason: 'BOTH_ENDPOINTS' });
  }
  return {
    selector: selectorPayload(seq?.region_record_id),
    start: Math.min(start, end),
    end: Math.max(start, end),
    reverseComplement: Boolean(seq?.region_reverse) || start > end
  };
};

// The Linear card fields that hold a request record's selection, crop, and
// orientation.
/** @param {RequestRecordView | null | undefined} record */
const linearSeqRecordFields = (record) => ({
  region_record_id: canonicalRecordSelector(record) || '',
  region_start: record?.region?.start ?? null,
  region_end: record?.region?.end ?? null,
  region_reverse: Boolean(record?.region?.reverseComplement || record?.presentation?.reverseComplement)
});

/**
 * @param {{ state: RequestState, drawing: RequestDrawing, filesData: RequestFilesData, resources: ResourceBuilder }} input
 */
const buildRecords = ({ state, drawing, filesData, resources }) => {
  if (state.mode.value === 'linear') {
    const resolvedRows = /** @type {Record<string, any>[]} */ (resolveEffectiveLinearRecordRows(
      filesData.linearSeqs,
      drawing.linearRecordRows,
      { enabled: Boolean(drawing.linearRecordLayoutEnabled?.value) }
    ));
    const canonicalCardinalityByUid = new Map(
      (drawing.linearRecordRows || []).map(
        (/** @type {{ uid?: unknown, canonicalCardinality?: string }} */ entry) => [entry.uid, entry.canonicalCardinality]
      )
    );
    /** @type {RequestDraftRecord[]} */
    const records = (filesData.linearSeqs || []).map((seq, index) => {
      /** @type {CanonicalRecord['source']} */
      const source = state.lInputType.value === 'gff'
        ? {
            kind: 'gffFasta',
            gffResourceId: resources.addFile(`record-${index + 1}-gff3`, 'gff3', seq.gff),
            fastaResourceId: resources.addFile(`record-${index + 1}-fasta`, 'fasta', seq.fasta)
          }
        : {
            kind: 'genbank',
            resourceId: resources.addFile(`record-${index + 1}-genbank`, 'genbank', seq.gb)
          };
      const region = linearRegionPayload(seq);
      const selector = region ? null : selectorPayload(seq.region_record_id);
      return {
        recordKey: String(seq.uid || `record-${index + 1}`),
        cardinality: canonicalCardinalityByUid.get(seq.uid)
          || seq.cardinality || (selector || region ? 'exactly_one' : 'all'),
        source,
        selector,
        region,
        presentation: {
          ...presentationPayload({
            label: resolveLinearRecordEffectiveDefinition(seq),
            subtitle: resolveLinearRecordEffectiveSubtitle(seq),
            gridRow: drawing.linearRecordLayoutEnabled?.value
              ? (resolvedRows[index]?.row ?? null)
              : null,
            gridColumn: drawing.linearRecordLayoutEnabled?.value
              ? (resolvedRows[index]?.canonicalColumn ?? null)
              : null
          }),
          reverseComplement: region ? false : Boolean(seq.region_reverse)
        }
      };
    });
    return { records, circularSourceIndexes: null, circularSourceCount: null, circularPresentedCount: null };
  }

  if (Array.isArray(filesData.circularRecords) && filesData.circularRecords.length > 0) {
    /** @type {RequestDraftRecord[]} */
    const records = filesData.circularRecords.map((record, index) => {
      /** @type {CanonicalRecord['source']} */
      const source = record.sourceKind === 'gffFasta'
        ? { kind: 'gffFasta',
            gffResourceId: resources.addFile(`record-${index + 1}-gff3`, 'gff3', record.gff),
            fastaResourceId: resources.addFile(`record-${index + 1}-fasta`, 'fasta', record.fasta) }
        : { kind: 'genbank', resourceId: resources.addFile(
            `record-${index + 1}-genbank`, 'genbank', record.gb) };
      return {
        recordKey: String(record.recordKey || `record-${index + 1}`),
        cardinality: record.cardinality || 'exactly_one',
        source,
        selector: record.selector || null,
        region: record.region || null,
        presentation: publicationClone(record.presentation) || presentationPayload(),
        display: publicationClone(record.display) || { isCircular: null, startCoordinate: null }
      };
    });
    const singleJourney = (
      records.length === 1 &&
      !drawing.form.multi_record_canvas &&
      drawing.adv.circular_grouping_intent !== 'batch'
    );
    if (singleJourney) {
      const record = records[0];
      const savedSelector = canonicalRecordSelector(record);
      const requestedSelector = String(
        drawing.form.circular_record_selector || savedSelector || ''
      ).trim();
      const knownRecords = (Array.isArray(state.circularRecordList.value)
        ? state.circularRecordList.value
        : []).map((entry) => ({
          ...entry,
          recordId: entry?.record_id ?? entry?.recordId
        }));
      const selection = resolveDisambiguatedRecordSelection(
        knownRecords,
        requestedSelector
      );
      const selectedKnownRecord = selection.record || (
        !requestedSelector && selection.entries.length === 1
          ? selection.entries[0]
          : null
      );
      if (knownRecords.length > 0 && !selectedKnownRecord) {
        throw diagnosticError('RECORD_SELECTION', {
          reason: !requestedSelector ? 'SELECT_ONE' : selection.status === 'ambiguous' ? 'AMBIGUOUS' : 'NO_MATCH'
        });
      }
      const selected = selectedKnownRecord || {
        value: requestedSelector,
        selector: requestedSelector,
        recordId: requestedSelector
      };
      const region = circularRegionPayload(drawing.form, selected);
      records[0] = {
        ...record,
        recordKey: record.recordKey || circularSingleRecordKey(selected),
        selector: region ? null : selectorPayload(selected.value),
        region,
        presentation: circularPresentationPayload(drawing.form, {
          hasRegion: Boolean(region)
        })
      };
    }
    return {
      records,
      circularSourceIndexes: records.map((_, index) => index),
      circularSourceCount: records.length,
      circularPresentedCount: records.length
    };
  }
  /** @type {CanonicalRecord['source']} */
  const source = state.cInputType.value === 'gff'
    ? {
        kind: 'gffFasta',
        gffResourceId: resources.addFile('record-1-gff3', 'gff3', filesData.c_gff),
        fastaResourceId: resources.addFile('record-1-fasta', 'fasta', filesData.c_fasta)
      }
    : {
        kind: 'genbank',
        resourceId: resources.addFile('record-1-genbank', 'genbank', filesData.c_gb)
      };
  const recordSet = resolveCircularRequestRecordSet(/** @type {Record<string, any>} */ ({
    records: state.circularRecordList.value,
    selector: drawing.form.circular_record_selector,
    multiRecordCanvas: drawing.form.multi_record_canvas,
    groupingIntent: drawing.adv.circular_grouping_intent,
    recordsOff: drawing.recordsOff
  }));
  if (recordSet.selectionFailure) throw diagnosticError('RECORD_SELECTION', { reason: recordSet.selectionFailure });
  const { recordSelectors } = recordSet;
  const selectedRecords = recordSet.records.length > 0 ? recordSet.records : [null];
  // The offered records decide the journey, so one ON record of several
  // stays a batch of one or a grid of one.
  const singleJourney = isSingleCircularJourney(recordSet);
  /** @type {RequestDraftRecord[]} */
  const records = selectedRecords.map((record, index) => {
    const region = singleJourney
      ? circularRegionPayload(drawing.form, record)
      : null;
    // A record keeps its key while others go OFF: `record-<file position>`.
    return {
      recordKey: singleJourney && record
        ? circularSingleRecordKey(record)
        : circularRecordRequestKey(Number.isInteger(record?.sourceIndex) ? record.sourceIndex : index),
      cardinality: 'exactly_one',
      source,
      selector: region
        ? null
        : selectorPayload(record?.value ?? record?.selector),
      region,
      presentation: singleJourney
        ? circularPresentationPayload(drawing.form, { hasRegion: Boolean(region) })
        : presentationPayload()
    };
  });
  const circularSourceIndexes = selectedRecords.map((record, index) => (
    record && Number.isInteger(record.sourceIndex) ? record.sourceIndex : index
  ));
  return {
    records,
    circularSourceIndexes,
    circularSourceCount: recordSelectors.length || records.length,
    circularPresentedCount: recordSet.presentedCount || records.length
  };
};

// The label rendering a Generate draws (OV-159): Auto while the mode shows no
// labels or a Linear label sits above its feature. The draft keeps the
// choice, so the control shows it again where it applies (R7).
/**
 * @param {Record<string, any>} form
 * @param {Record<string, any>} adv
 * @param {boolean} circular
 * @param {string | null} linearLabelPlacement
 */
const effectiveLabelRendering = (form, adv, circular, linearLabelPlacement) => {
  const labelsMode = circular
    ? String((typeof form.labels_mode === 'string'
      ? form.labels_mode
      : (form.allow_inner_labels ? 'both' : (form.show_labels ? 'out' : 'none'))) || 'none').trim().toLowerCase()
    : String(form.show_labels_linear || 'none');
  if (labelsMode === 'none' || linearLabelPlacement === 'above_feature') return 'auto';
  const normalized = String(adv.label_rendering || '').trim().toLowerCase();
  return ['embedded_only', 'external_only'].includes(normalized) ? normalized : 'auto';
};

/**
 * @param {RequestState} state
 * @param {RequestDrawing} drawing
 * @param {{ depthRequested?: boolean, hasComparisonIntent?: boolean, linearHasSharedRow?: boolean }} [options]
 */
const buildConfigOverrides = (
  state,
  drawing,
  {
    depthRequested = Boolean(drawing.form.show_depth),
    hasComparisonIntent = false,
    linearHasSharedRow = false
  } = {}
) => {
  const { form, adv } = drawing;
  const circular = state.mode.value === 'circular';
  // R7: a config leaf is projected literally; Python validates its domain.
  /** @param {string} configPath @param {unknown} value */
  const leaf = (configPath, value) => projectOptionalNumber(value, { configPath });
  /** @param {string} path @param {unknown} value */
  const sharedLeaf = (path, value) => leaf(`${path}.short`, value);
  const linearLabelPlacement = circular
    ? null
    : requireCurrentLinearLabelPlacement(adv.label_placement);
  const linearTrackLayout = circular
    ? null
    : requireCurrentLinearTrackLayout(form.linear_track_layout);
  const comparisonHeight = !circular && hasComparisonIntent
    ? classifyOptionalPositiveNumber(adv.comparison_height)
    : null;
  if (comparisonHeight?.status === 'invalid') {
    throw diagnosticError('INPUT_INVALID', { field: 'match_height', reason: 'POSITIVE_OR_AUTO' });
  }
  // An axis color still at the Linear default follows the ruler.
  const linearAxisManaged = isModeProfileDefault('linear', 'axis_stroke_color', adv.axis_stroke_color);
  const linearAxisStrokeColor = circular
    ? null
    : effectiveLinearAxisColor({
        axisColor: adv.axis_stroke_color,
        rulerOnAxis: Boolean(form.linear_ruler_on_axis) && linearAxisRulerApplies(form),
        managed: linearAxisManaged
      });
  /** @type {Record<string, unknown>} */
  const overrides = {
    [CONFIG_OVERRIDE_PATHS.arrowHeadLengthRatio]:
      normalizeArrowHeadLengthRatio(adv.arrow_head_length_ratio),
    [CONFIG_OVERRIDE_PATHS.arrowShaftWidthRatio]:
      normalizeArrowShaftWidthRatio(adv.arrow_shaft_width_ratio),
    [CONFIG_OVERRIDE_PATHS.blockStrokeColor]: adv.block_stroke_color || null,
    [CONFIG_OVERRIDE_PATHS.lineStrokeColor]: adv.line_stroke_color || null,
    [CONFIG_OVERRIDE_PATHS.fontFamily]: adv.font_family || null,
    [CONFIG_OVERRIDE_PATHS.labelRendering]: effectiveLabelRendering(form, adv, circular, linearLabelPlacement),
    [CONFIG_OVERRIDE_PATHS.showGc]: circular ? !form.suppress_gc : Boolean(form.show_gc),
    [CONFIG_OVERRIDE_PATHS.showSkew]: circular
      ? !form.suppress_skew
      : Boolean(form.show_skew),
    [CONFIG_OVERRIDE_PATHS.showDepth]: Boolean(depthRequested),
    [MODE_LABEL_SCOPE_PATHS[state.mode.value]]: labelScopeOverride(form, circular),
    [CONFIG_OVERRIDE_PATHS.strandedness]: Boolean(form.separate_strands),
    [CONFIG_OVERRIDE_PATHS.resolveOverlaps]: Boolean(adv.resolve_overlaps),
    [CONFIG_OVERRIDE_PATHS.featureOverlapToleranceBp]: leaf(CONFIG_OVERRIDE_PATHS.featureOverlapToleranceBp, adv.feature_overlap_tolerance_bp),
    [CONFIG_OVERRIDE_PATHS.gcContentMode]: adv.gc_content_mode || 'deviation',
    [CONFIG_OVERRIDE_PATHS.gcContentMinPercent]: leaf(CONFIG_OVERRIDE_PATHS.gcContentMinPercent, adv.gc_content_min_percent),
    [CONFIG_OVERRIDE_PATHS.gcContentMaxPercent]: leaf(CONFIG_OVERRIDE_PATHS.gcContentMaxPercent, adv.gc_content_max_percent),
    [CONFIG_OVERRIDE_PATHS.gcContentShowAxis]: Boolean(adv.gc_content_show_axis),
    [CONFIG_OVERRIDE_PATHS.gcContentShowTicks]: Boolean(adv.gc_content_show_ticks),
    [CONFIG_OVERRIDE_PATHS.gcContentLargeTickInterval]:
      leaf(CONFIG_OVERRIDE_PATHS.gcContentLargeTickInterval, adv.gc_content_tick_interval),
    [CONFIG_OVERRIDE_PATHS.gcContentSmallTickInterval]:
      leaf(CONFIG_OVERRIDE_PATHS.gcContentSmallTickInterval, adv.gc_content_small_tick_interval),
    [CONFIG_OVERRIDE_PATHS.gcContentTickFontSize]:
      leaf(CONFIG_OVERRIDE_PATHS.gcContentTickFontSize, adv.gc_content_tick_font_size),
    [CONFIG_OVERRIDE_PATHS.depthColor]: adv.depth_color || null,
    [CONFIG_OVERRIDE_PATHS.depthMin]: leaf(CONFIG_OVERRIDE_PATHS.depthMin, adv.depth_min),
    [CONFIG_OVERRIDE_PATHS.depthMax]: leaf(CONFIG_OVERRIDE_PATHS.depthMax, adv.depth_max),
    [CONFIG_OVERRIDE_PATHS.depthNormalize]: Boolean(adv.depth_normalize),
    [CONFIG_OVERRIDE_PATHS.depthShowAxis]: Boolean(adv.depth_show_axis),
    [CONFIG_OVERRIDE_PATHS.depthShowTicks]: Boolean(adv.depth_show_ticks),
    [CONFIG_OVERRIDE_PATHS.depthLargeTickInterval]:
      leaf(CONFIG_OVERRIDE_PATHS.depthLargeTickInterval, adv.depth_large_tick_interval),
    [CONFIG_OVERRIDE_PATHS.depthSmallTickInterval]:
      leaf(CONFIG_OVERRIDE_PATHS.depthSmallTickInterval, adv.depth_small_tick_interval),
    [CONFIG_OVERRIDE_PATHS.depthTickFontSize]: leaf(CONFIG_OVERRIDE_PATHS.depthTickFontSize, adv.depth_tick_font_size),
    [CONFIG_OVERRIDE_PATHS.depthShareAxis]: Boolean(adv.depth_share_axis),
    [CONFIG_OVERRIDE_PATHS.showScale]: form.show_scale !== false,
    [CONFIG_OVERRIDE_PATHS.scaleInterval]: leaf(CONFIG_OVERRIDE_PATHS.scaleInterval, adv.scale_interval),
    [CONFIG_OVERRIDE_PATHS.labelBlacklist]: labelBlacklistOverride(drawing),
    ...(circular
      ? {
          [CONFIG_OVERRIDE_PATHS.circularAxisStrokeColor]:
            adv.axis_stroke_color || null,
          [CONFIG_OVERRIDE_PATHS.circularDefinitionFontSize]:
            leaf(CONFIG_OVERRIDE_PATHS.circularDefinitionFontSize, adv.def_font_size),
          [CONFIG_OVERRIDE_PATHS.circularDefinitionInterval]: leaf(CONFIG_OVERRIDE_PATHS.circularDefinitionInterval, adv.circular_definition_interval),
          [CONFIG_OVERRIDE_PATHS.plotTitleFontSize]:
            leaf(CONFIG_OVERRIDE_PATHS.plotTitleFontSize, adv.plot_title_font_size),
          [CONFIG_OVERRIDE_PATHS.circularLabelSpacing]:
            leaf(CONFIG_OVERRIDE_PATHS.circularLabelSpacing, adv.circular_label_spacing),
          [CONFIG_OVERRIDE_PATHS.circularLabelPlacement]:
            adv.circular_label_placement || 'horizontal',
          [CONFIG_OVERRIDE_PATHS.trackType]: form.track_type,
          [CONFIG_OVERRIDE_PATHS.tickLabelFontSize]:
            leaf(CONFIG_OVERRIDE_PATHS.tickLabelFontSize, adv.tick_label_font_size),
          [CONFIG_OVERRIDE_PATHS.outerLabelXRadiusOffset]:
            leaf(CONFIG_OVERRIDE_PATHS.outerLabelXRadiusOffset, adv.outer_label_x_offset),
          [CONFIG_OVERRIDE_PATHS.outerLabelYRadiusOffset]:
            leaf(CONFIG_OVERRIDE_PATHS.outerLabelYRadiusOffset, adv.outer_label_y_offset),
          [CONFIG_OVERRIDE_PATHS.innerLabelXRadiusOffset]:
            leaf(CONFIG_OVERRIDE_PATHS.innerLabelXRadiusOffset, adv.inner_label_x_offset),
          [CONFIG_OVERRIDE_PATHS.innerLabelYRadiusOffset]:
            leaf(CONFIG_OVERRIDE_PATHS.innerLabelYRadiusOffset, adv.inner_label_y_offset)
        }
      : {
          [CONFIG_OVERRIDE_PATHS.linearAxisStrokeColor]: linearAxisStrokeColor,
          [CONFIG_OVERRIDE_PATHS.linearDefinitionShowReplicon]:
            Boolean(adv.linear_show_replicon),
          [CONFIG_OVERRIDE_PATHS.linearDefinitionShowAccession]:
            resolveLinearLabelVisibility(adv.linear_accession_visibility, {
              hasSharedRow: linearHasSharedRow
            }),
          [CONFIG_OVERRIDE_PATHS.linearDefinitionShowLength]:
            resolveLinearLabelVisibility(adv.linear_length_visibility, {
              hasSharedRow: linearHasSharedRow
            }),
          [CONFIG_OVERRIDE_PATHS.linearLabelSpacing]:
            leaf(CONFIG_OVERRIDE_PATHS.linearLabelSpacing, adv.linear_label_spacing),
          [CONFIG_OVERRIDE_PATHS.labelPlacement]: linearLabelPlacement,
          [CONFIG_OVERRIDE_PATHS.labelRotation]: leaf(CONFIG_OVERRIDE_PATHS.labelRotation, adv.label_rotation),
          [CONFIG_OVERRIDE_PATHS.alignCenter]: Boolean(form.align_center),
          [CONFIG_OVERRIDE_PATHS.keepDefinitionLeftAligned]:
            Boolean(form.keep_definition_left_aligned),
          [CONFIG_OVERRIDE_PATHS.linearTrackLayout]: linearTrackLayout,
          [CONFIG_OVERRIDE_PATHS.linearTrackAxisGap]: leaf(CONFIG_OVERRIDE_PATHS.linearTrackAxisGap, adv.track_axis_gap),
          [CONFIG_OVERRIDE_PATHS.linearRulerOnAxis]: Boolean(form.linear_ruler_on_axis),
          ...(hasComparisonIntent
            ? {
                [CONFIG_OVERRIDE_PATHS.comparisonHeight]:
                  // hasComparisonIntent && !circular holds here, so comparisonHeight is set.
                  /** @type {NonNullable<typeof comparisonHeight>} */ (comparisonHeight).status === 'auto'
                    ? null
                    : /** @type {NonNullable<typeof comparisonHeight>} */ (comparisonHeight).value,
                [CONFIG_OVERRIDE_PATHS.pairwiseMatchStyle]: adv.pairwise_match_style
              }
            : {}),
          [CONFIG_OVERRIDE_PATHS.gcHeight]: leaf(CONFIG_OVERRIDE_PATHS.gcHeight, adv.gc_height),
          [CONFIG_OVERRIDE_PATHS.depthHeight]: leaf(CONFIG_OVERRIDE_PATHS.depthHeight, adv.depth_height),
          [CONFIG_OVERRIDE_PATHS.scaleStyle]: form.scale_style,
          [CONFIG_OVERRIDE_PATHS.scaleStrokeColor]: adv.scale_stroke_color || null,
          [CONFIG_OVERRIDE_PATHS.scaleLabelColor]: adv.ruler_label_color || null,
          [CONFIG_OVERRIDE_PATHS.scaleStrokeWidth]:
            leaf(CONFIG_OVERRIDE_PATHS.scaleStrokeWidth, adv.scale_stroke_width),
          [CONFIG_OVERRIDE_PATHS.normalizeLength]: Boolean(form.normalize_length)
        })
  };
  const sharedLengthValues = {
    [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.blockStrokeWidth]:
      sharedLeaf(SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.blockStrokeWidth, adv.block_stroke_width),
    [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.lineStrokeWidth]:
      sharedLeaf(SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.lineStrokeWidth, adv.line_stroke_width),
    [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.legendBoxSize]:
      sharedLeaf(SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.legendBoxSize, adv.legend_box_size),
    [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.legendFontSize]:
      sharedLeaf(SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.legendFontSize, adv.legend_font_size),
    ...(circular
      ? {
          [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.circularAxisStrokeWidth]:
            sharedLeaf(SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.circularAxisStrokeWidth, adv.axis_stroke_width),
          'labels.font_size': sharedLeaf('labels.font_size', adv.label_font_size)
        }
      : {
          [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.linearAxisStrokeWidth]:
            sharedLeaf(SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.linearAxisStrokeWidth, adv.axis_stroke_width),
          [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.linearDefinitionFontSize]:
            sharedLeaf(SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.linearDefinitionFontSize, adv.def_font_size),
          [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.defaultCdsHeight]:
            sharedLeaf(SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.defaultCdsHeight, adv.feature_height),
          [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.scaleFontSize]:
            sharedLeaf(SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.scaleFontSize, adv.scale_font_size),
          [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.rulerLabelFontSize]:
            sharedLeaf(SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.rulerLabelFontSize, adv.ruler_label_font_size),
          'labels.font_size.linear': sharedLeaf('labels.font_size.linear', adv.label_font_size)
        })
  };
  for (const [path, value] of Object.entries(sharedLengthValues)) {
    if (value === null || value === undefined) continue;
    overrides[`${path}.short`] = value;
    overrides[`${path}.long`] = value;
  }
  if (!circular) {
    for (const [kind, prefix] of Object.entries(LINEAR_DEFINITION_STYLE_PATHS)) {
      const style = adv.linear_definition_line_styles?.[kind];
      if (!style || typeof style !== 'object' || Array.isArray(style)) continue;
      for (const field of LINEAR_DEFINITION_STYLE_FIELDS) {
        const value = style[field];
        if (value !== null && value !== undefined) {
          overrides[`${prefix}.${field}`] = value;
        }
      }
    }
  }
  const managedOverrides = Object.fromEntries(
    Object.entries(overrides).filter(([, value]) => value !== null && value !== undefined)
  );
  const preservedOverrides = drawing.unmanagedConfigOverrides;
  const plainPreservedOverrides = (
    preservedOverrides && typeof preservedOverrides === 'object' && !Array.isArray(preservedOverrides)
  )
    ? cloneJsonData(preservedOverrides)
    : {};
  return {
    ...plainPreservedOverrides,
    ...managedOverrides
  };
};

/**
 * Whether request `comparisons` are only the protein settings that the CLI
 * writes without -b or a protein mode: disabled pipelines (mode `none`) with
 * no pairs, which draw no comparison.
 * @param {any[]} comparisons
 */
const drawsNoComparison = (comparisons) => comparisons.length > 0
  && comparisons.every((comparison) => comparison?.kind === 'generatedProteinComparison'
    && comparison.mode === 'none' && !comparison.pairs?.length);

/** @type {Readonly<Record<string, boolean>>} */
const ALL_GENERATED_TABLES = Object.freeze({
  colors: true, visibility: true, whitelist: true, priority: true, labelOverrides: true
});

// The label edits a request carries: the bulk label edits projected onto the
// displayed Result's labels and the recorded label sources, as identity row
// text, and the `* * label` rows of the rest (CW-02: built once per request).
/**
 * @param {RequestLabelState} state
 * @param {RequestDrawing} drawing
 */
export const requestLabelProjection = (state, drawing) => {
  const bulk = drawing.labelTextBulkOverrides || {};
  if (Object.keys(bulk).length === 0) return { bulkLabelText: {}, rows: [] };
  const displayed = resultRenderedFeatures(state);
  const labelTargets = (state.editableLabels?.value || []).map((entry) => ({
    identityKey: featureIdentityKeyOf(displayed?.get(String(entry?.featureId || '').trim())),
    sourceText: entry?.sourceText
  })).filter((target) => target.identityKey);
  return buildBulkLabelProjection(bulk, { labelTargets, featureOverrides: drawing.featureOverrides });
};

/**
 * @param {RequestTableState} state
 * @param {RequestDrawing} drawing
 * @param {ResourceBuilder} resources
 * @param {RequestTableOptions} diagramOptions
 * @param {Readonly<Record<string, boolean>>} [tables]
 * @param {ReturnType<typeof requestLabelProjection>} [labelProjection]
 */
const addGeneratedTableResources = (
  state, drawing, resources, diagramOptions, tables = ALL_GENERATED_TABLES,
  labelProjection = requestLabelProjection(state, drawing)
) => {
  recordStructuralMetric('generatedTableBuildCount');
  if (tables.colors) addColorTableResources(state, drawing, resources, diagramOptions);
  if (tables.visibility) {
    const visibility = serializeFeatureVisibilityRules(drawing.featureVisibilityRules?.value || []);
    if (visibility.trim()) {
      diagramOptions.featureVisibilityTableFile = fileRef(resources.addText(
        'feature-visibility-table-file', 'feature-visibility-table-file', 'feature-visibility.tsv', visibility
      ));
    }
  }
  if (tables.whitelist) addLabelWhitelistResource(state, drawing, resources, diagramOptions);
  if (tables.priority) addQualifierPriorityResource(state, drawing, resources, diagramOptions);
  if (tables.labelOverrides) addLabelOverrideResource(state, drawing, resources, diagramOptions, labelProjection);
};

/**
 * @param {RequestTableState} state
 * @param {RequestDrawing} drawing
 * @param {ResourceBuilder} resources
 * @param {RequestTableOptions} diagramOptions
 */
const addColorTableResources = (state, drawing, resources, diagramOptions) => {
  const paletteName = String(drawing.selectedPalette.value || 'default');
  const paletteColors = state.canonicalPublicationFiles && !state.canonicalPublicationFiles.d_color ? {}
    : normalizePaletteColors(state.paletteDefinitions.value?.[paletteName]
      || state.paletteDefinitions.value?.default || {});
  const defaultColors = buildDefaultColorOverrideTsv({
    colors: drawing.currentColors.value,
    paletteColors
  });
  const specificColors = serializeSpecificRules(drawing.manualSpecificRules);
  const publicationFiles = state.canonicalPublicationFiles || {};
  const defaultColorsFile = publicationFileRef(
    resources, publicationFiles, 'd_color', 'colors-default-colors-file'
  ) || (defaultColors.trim() ? fileRef(resources.addText(
        'colors-default-colors-file', 'colors-default-colors-file',
        'default-colors.tsv', `${defaultColors}\n`
      )) : null);
  const colorTableFile = publicationFileRef(
    resources, publicationFiles, 't_color', 'colors-color-table-file'
  ) || (specificColors.trim() ? fileRef(resources.addText(
        'colors-color-table-file', 'colors-color-table-file',
        'specific-colors.tsv', specificColors
      )) : null);
  diagramOptions.colors = {
    colorTable: colorTableFile?.representation === 'canonicalTsv' ? colorTableFile : null,
    colorTableFile: colorTableFile?.representation === 'canonicalTsv' ? null : colorTableFile,
    defaultColors: defaultColorsFile?.representation === 'canonicalTsv' ? defaultColorsFile : null,
    defaultColorsPalette: paletteName,
    defaultColorsFile:
      defaultColorsFile?.representation === 'canonicalTsv' ? null : defaultColorsFile
  };
};

/**
 * @param {RequestTableState} state
 * @param {RequestDrawing} drawing
 * @param {ResourceBuilder} resources
 * @param {RequestTableOptions} diagramOptions
 */
const addLabelWhitelistResource = (state, drawing, resources, diagramOptions) => {
  const publicationFiles = state.canonicalPublicationFiles || {};
  const preservedWhitelist = publicationFileRef(
    resources, publicationFiles, 'whitelist', 'label-whitelist-file'
  );
  if (preservedWhitelist) {
    if (preservedWhitelist.representation === 'canonicalTsv') diagramOptions.labelWhitelistTable = preservedWhitelist;
    else diagramOptions.labelWhitelistFile = preservedWhitelist;
  } else if (drawing.filterMode.value === 'Whitelist' && drawing.manualWhitelist.length > 0) {
    const whitelist = serializeLabelWhitelistRules(drawing.manualWhitelist);
    if (whitelist) {
      diagramOptions.labelWhitelistFile = fileRef(resources.addText(
        'label-whitelist-file', 'label-whitelist-file', 'label-whitelist.tsv', whitelist
      ));
    }
  }
};

/**
 * @param {RequestTableState} state
 * @param {RequestDrawing} drawing
 * @param {ResourceBuilder} resources
 * @param {RequestTableOptions} diagramOptions
 */
const addQualifierPriorityResource = (state, drawing, resources, diagramOptions) => {
  const publicationFiles = state.canonicalPublicationFiles || {};
  const priority = serializeQualifierPriorityRules(drawing.manualPriorityRules);
  const priorityRef = publicationFileRef(
    resources, publicationFiles, 'qualifier_priority', 'qualifier-priority-file'
  );
  if (priorityRef) {
    if (priorityRef.representation === 'canonicalTsv') diagramOptions.qualifierPriorityTable = priorityRef;
    else diagramOptions.qualifierPriorityFile = priorityRef;
  } else if (priority) {
    diagramOptions.qualifierPriorityFile = fileRef(resources.addText(
      'qualifier-priority-file', 'qualifier-priority-file', 'qualifier-priority.tsv', priority
    ));
  }
};

// The label table carries rules only: a saved table, else the bulk label edits
// no displayed label resolves. Per-feature label edits are `featureOverrides`.
/**
 * @param {RequestLabelState} state
 * @param {RequestDrawing} drawing
 * @param {{ rows: string[] }} [labelProjection]
 */
export const requestLabelTableTsv = (state, drawing, labelProjection = requestLabelProjection(state, drawing)) => (
  serializeLabelOverrideRows(drawing.canonicalLabelOverrideRows?.value)
    || (labelProjection.rows.length ? `${labelProjection.rows.join('\n')}\n` : '')
);

/**
 * @param {RequestTableState} state
 * @param {RequestDrawing} drawing
 * @param {ResourceBuilder} resources
 * @param {RequestTableOptions} diagramOptions
 * @param {{ rows: string[] }} labelProjection
 */
const addLabelOverrideResource = (state, drawing, resources, diagramOptions, labelProjection) => {
  const labelOverrideTsv = requestLabelTableTsv(state, drawing, labelProjection);
  if (labelOverrideTsv) {
    diagramOptions.labelOverrideFile = fileRef(resources.addText(
      'label-override-file', 'label-override-file', 'label-overrides.tsv', labelOverrideTsv
    ));
  }
};

/** @param {{ state: RequestState, filesData: RequestFilesData, recordCount: number }} input */
const depthSourceRowsForRequest = ({ state, filesData, recordCount }) => (
  state.mode.value === 'linear'
    ? (filesData.linearSeqs || []).map((seq) => (
        Array.isArray(seq.depth) ? seq.depth : (seq.depth ? [seq.depth] : [])
      ))
    : normalizeRecordMajorDepthFileRows(filesData.c_depth, recordCount)
);

/** @param {{ state: RequestState, filesData: RequestFilesData, recordCount: number }} input */
const logicalDepthTrackCountForRequest = ({ state, filesData, recordCount }) => {
  const rows = depthSourceRowsForRequest({ state, filesData, recordCount });
  return rows.reduce((maximum, row) => (
    Math.max(maximum, Array.isArray(row) ? row.length : 0)
  ), 0);
};

/**
 * @param {{ state: RequestState, drawing: RequestDrawing, filesData: RequestFilesData,
 *   resources: ResourceBuilder, diagramOptions: Record<string, any>, recordCount: number }} input
 */
const buildDepthResources = ({ state, drawing, filesData, resources, diagramOptions, recordCount }) => {
  if (
    state.mode.value === 'circular' &&
    isRecordMajorDepthFileMatrix(filesData.c_depth) &&
    filesData.c_depth.length !== recordCount
  ) {
    throw new Error(
      `Circular Depth matrix has ${filesData.c_depth.length} record rows; expected ${recordCount}.`
    );
  }
  const rows = depthSourceRowsForRequest({ state, filesData, recordCount });
  if (!rows.some((row) => row.some(Boolean))) return;
  if (rows.length !== recordCount || rows.some((row) => !Array.isArray(row))) {
    throw new Error(
      `${state.mode.value === 'circular' ? 'Circular' : 'Linear'} Depth sources must contain one row per record (${recordCount}).`
    );
  }
  const tracks = Array.isArray(drawing.adv.depth_tracks) ? drawing.adv.depth_tracks : [];
  const logicalTrackCount = Math.max(
    tracks.length,
    ...rows.map((row) => row.length)
  );
  diagramOptions.depthTracks = Array.from({ length: logicalTrackCount }, (_, trackIndex) => {
    const sources = rows.map((row) => row[trackIndex] || null);
    if (!sources.some(Boolean)) throw depthSeriesWithoutSource(trackIndex);
    const sharedSource = sources[0] && sources.every((source) => source === sources[0]);
    const sourceName = `depth-tracks-${trackIndex + 1}-source`;
    const source = sharedSource
      ? fileRef(resources.addFile(sourceName, 'depth-track-file', sources[0]))
      : sources.map((entry, recordIndex) => (
          entry
            ? fileRef(resources.addFile(
                `${sourceName}-record-${recordIndex + 1}`,
                'depth-track-file',
                entry
              ))
            : null
        ));
    const track = tracks[trackIndex] || {};
    return {
      source,
      label: String(track.label || (logicalTrackCount === 1 ? 'Depth' : `Depth ${trackIndex + 1}`)),
      color: String(track.color || drawing.adv.depth_color || '#4A90E2'),
      height: state.mode.value === 'linear'
        ? projectOptionalNumber(track.height, { field: 'height', seriesIndex: trackIndex })
        : null,
      largeTickInterval: projectOptionalNumber(track.large_tick_interval, { field: 'large_tick_interval', seriesIndex: trackIndex }),
      smallTickInterval: projectOptionalNumber(track.small_tick_interval, { field: 'small_tick_interval', seriesIndex: trackIndex }),
      tickFontSize: projectOptionalNumber(track.tick_font_size, { field: 'tick_font_size', seriesIndex: trackIndex })
    };
  });
};

/** @param {RequestDrawing} drawing */
const circularGeometryShortcutsForState = (drawing) => ({
  featureWidth: drawing.adv.feature_width_circular,
  depthWidth: drawing.adv.depth_width_circular,
  gcContentWidth: drawing.adv.gc_content_width_circular,
  gcContentRadius: drawing.adv.gc_content_radius_circular,
  gcSkewWidth: drawing.adv.gc_skew_width_circular,
  gcSkewRadius: drawing.adv.gc_skew_radius_circular
});

/** @param {RequestDrawing} drawing */
const annotationSetIdsForState = (drawing) => (
  (Array.isArray(drawing.annotationSets) ? drawing.annotationSets : [])
    .map((set) => String(set?.id || '').trim())
    .filter(Boolean)
);

/** @param {RequestDrawing} drawing */
const automaticCircularAnnotationSlots = (drawing) => {
  /** @type {CircularTrackSlot[]} */
  const slots = [];
  (Array.isArray(drawing.annotationSets) ? drawing.annotationSets : [])
    .forEach((set, index) => {
      const setId = String(set?.id || '').trim();
      if (!setId) return;
      const marks = Array.from(new Set(
        (Array.isArray(set.annotations) ? set.annotations : [])
          .map((/** @type {{ mark?: unknown } | null} */ annotation) => String(annotation?.mark || '').trim().toLowerCase())
          .filter(Boolean)
      ));
      const laneMarks = marks.filter((mark) => mark !== 'highlight');
      if (laneMarks.length > 0) {
        slots.push({
          id: `annotations_${index + 1}`,
          renderer: 'annotations',
          enabled: true,
          width: null,
          radius: null,
          inner_gap_px: null,
          outer_gap_px: null,
          side: 'outside',
          z: 0,
          params: { set_id: setId, marks: laneMarks }
        });
      }
      if (marks.includes('highlight')) {
        slots.push({
          id: `annotations_${index + 1}${laneMarks.length > 0 ? '_highlight' : ''}`,
          renderer: 'annotations',
          enabled: true,
          width: null,
          radius: null,
          inner_gap_px: null,
          outer_gap_px: null,
          side: 'overlay',
          z: -1,
          params: {
            set_id: setId,
            marks: ['highlight'],
            anchor_slot: 'features',
            layer: 'underlay',
            cover_anchor: true,
            padding_px: 0
          }
        });
      }
    });
  return slots;
};

/**
 * @param {{ reference?: unknown, ring_width?: unknown, ring_gap?: unknown }} conservation
 * @param {{ label: string, color: string }[]} entries
 * @param {string} referenceDefault
 */
const conservationDiagramOptions = (conservation, entries, referenceDefault) => ({
  conservationReference: String(conservation.reference || referenceDefault),
  conservationLabels: entries.map((entry) => entry.label),
  conservationColors: entries.map((entry) => entry.color),
  conservationRingWidth: projectOptionalNumber(conservation.ring_width, { field: 'conservation_ring_width' }),
  conservationRingGap: projectOptionalNumber(conservation.ring_gap, { field: 'conservation_ring_gap' })
});

/**
 * @param {{ drawing: RequestDrawing, filesData: RequestFilesData,
 *   resolvedCircularConservation: ResolvedConservationEntry[] | undefined }} input
 */
const conservationSeriesForValidation = ({
  drawing,
  filesData,
  resolvedCircularConservation
}) => {
  if (drawing.circularConservation?.enabled !== true) return [];
  const conservationSource = String(
    drawing.circularConservation?.source || ''
  ).trim().toLowerCase();
  const blastsAreDerived = filesData.c_conservation_blasts_source === 'losat-cache';
  const sourceFiles = (
    conservationSource === 'upload' || blastsAreDerived
  )
    ? filesData.c_conservation_blasts
    : filesData.c_conservation_fastas;
  /** @type {OrderedConservationSource[]} */
  let ordered = [];
  try {
    ordered = orderedConservationSources(
      Array.isArray(sourceFiles) ? sourceFiles : [],
      drawing.circularConservation || {}
    );
  } catch {
    ordered = [];
  }
  const resolved = (Array.isArray(resolvedCircularConservation)
    ? resolvedCircularConservation
    : []
  ).map((entry, index) => ({
    ...entry,
    sourceKey: String(entry?.sourceKey || entry?.seriesKey || `resolved-${index + 1}`),
    sourceIndex: index,
    orderIndex: index
  }));
  return [...ordered, ...resolved];
};

/**
 * The custom stack's plan; throws on a row issue. With `omitUnsourcedDepth`
 * (the request a Save builds for a mode not generated yet, BA-OD114), an
 * enabled Depth row whose series has no source in any record (the
 * validator's `depth_source_missing`) is left out, as it draws nothing; the
 * draft keeps the row. Generate stops on that row issue before it builds a
 * request (PD-OI-083).
 * @param {CustomTrackPlanInput} input
 * @param {{ omitUnsourcedDepth: boolean, depthRows: DepthFileSlots[] }} depth
 * @returns {CustomTrackPlan}
 */
const assertedCustomTrackPlan = (input, { omitUnsourcedDepth, depthRows }) => {
  if (!omitUnsourcedDepth) return assertValidCustomTrackPlan(validateCustomTrackPlan(input));
  const sourced = { ...input, depthSourcedTrackIndexes: activeDepthTrackIndices(depthRows) };
  const plan = validateCustomTrackPlan(sourced);
  const unsourced = new Set([...plan.rowIssues]
    .filter(([, issues]) => issues.some(({ code }) => code === 'depth_source_missing'))
    .map(([index]) => index));
  if (unsourced.size === 0) return assertValidCustomTrackPlan(plan);
  const slots = Array.isArray(input.slots) ? input.slots : [];
  return assertValidCustomTrackPlan(validateCustomTrackPlan({
    ...sourced,
    slots: slots.map((slot, index) => (unsourced.has(index) ? { ...slot, enabled: false } : slot))
  }));
};

/**
 * @param {{ state: RequestState, drawing: RequestDrawing, filesData: RequestFilesData,
 *   recordCount: number, resolvedCircularConservation: ResolvedConservationEntry[] | undefined,
 *   omitUnsourcedDepth: boolean }} input
 */
const buildTrackPlan = ({
  state,
  drawing,
  filesData,
  recordCount,
  resolvedCircularConservation,
  omitUnsourcedDepth
}) => {
  const circular = state.mode.value === 'circular';
  const depthTrackCount = logicalDepthTrackCountForRequest({
    state,
    filesData,
    recordCount
  });
  const depth = { omitUnsourcedDepth, depthRows: depthSourceRowsForRequest({ state, filesData, recordCount }) };
  const annotationSetIds = annotationSetIdsForState(drawing);

  if (circular && drawing.adv.circular_track_slots_enabled) {
    const validation = assertedCustomTrackPlan({
      mode: 'circular',
      slots: drawing.adv.circular_track_slots,
      axisIndex: drawing.adv.circular_track_slots_axis_index,
      trackType: drawing.form.track_type,
      depthTrackCount,
      annotationSetIds,
      conservationSeries: conservationSeriesForValidation({
        drawing,
        filesData,
        resolvedCircularConservation
      })
    }, depth);
    const depthRequested = validation.enabledSlots.some(
      (slot) => slot.renderer === 'depth'
    );
    return {
      depthRequested,
      tracks: {
        circularTrackSlots: validation.enabledSlots.map((slot) => (
          buildCircularTrackSlotPayload(slot, drawing.adv.nt, drawing.form.track_type)
        )),
        circularTrackAxisIndex: validation.emittedAxisIndex,
        linearTrackSlots: null,
        linearTrackAxisIndex: null,
        centerReservedRadius: projectOptionalNumber(drawing.adv.center_reserved_radius, { field: 'center_reserved_radius' })
      }
    };
  }

  if (!circular && drawing.adv.linear_track_slots_enabled) {
    const validation = assertedCustomTrackPlan({
      mode: 'linear',
      slots: drawing.adv.linear_track_slots,
      axisIndex: drawing.adv.linear_track_slots_axis_index,
      trackType: drawing.form.linear_track_layout,
      depthTrackCount,
      annotationSetIds,
      conservationSeries: []
    }, depth);
    const depthRequested = validation.enabledSlots.some(
      (slot) => slot.renderer === 'depth'
    );
    return {
      depthRequested,
      tracks: {
        circularTrackSlots: null,
        circularTrackAxisIndex: null,
        linearTrackSlots: validation.enabledSlots.map(buildLinearTrackSlotPayload),
        linearTrackAxisIndex: validation.emittedAxisIndex,
        centerReservedRadius: null
      }
    };
  }

  const depthRequested = Boolean(drawing.form.show_depth);
  if (circular) {
    const shortcuts = circularGeometryShortcutsForState(drawing);
    if (hasCircularGeometryShortcuts(shortcuts)) {
      const implicitSlots = applyCircularGeometryShortcuts(
        createDefaultCircularTrackSlots({
          nt: drawing.adv.nt,
          showDepth: depthRequested,
          depthTrackCount: Math.max(1, depthTrackCount),
          showGc: !drawing.form.suppress_gc,
          showSkew: !drawing.form.suppress_skew,
          showTicks: drawing.form.show_scale !== false,
          preset: drawing.form.track_type
        }),
        shortcuts
      );
      const slots = [
        ...automaticCircularAnnotationSlots(drawing),
        ...implicitSlots
      ];
      const axisIndex = inferLegacyAxisIndexFromFeature(
        slots,
        drawing.form.track_type
      );
      const validation = assertValidCustomTrackPlan(validateCustomTrackPlan({
        mode: 'circular',
        slots,
        axisIndex,
        trackType: drawing.form.track_type,
        depthTrackCount,
        annotationSetIds,
        conservationSeries: []
      }));
      return {
        depthRequested,
        tracks: {
          circularTrackSlots: validation.enabledSlots.map((slot) => (
            buildCircularTrackSlotPayload(slot, drawing.adv.nt, drawing.form.track_type)
          )),
          circularTrackAxisIndex: validation.emittedAxisIndex,
          linearTrackSlots: null,
          linearTrackAxisIndex: null,
          centerReservedRadius: projectOptionalNumber(drawing.adv.center_reserved_radius, { field: 'center_reserved_radius' })
        }
      };
    }
  }

  return {
    depthRequested,
    tracks: {
      circularTrackSlots: null,
      circularTrackAxisIndex: null,
      linearTrackSlots: null,
      linearTrackAxisIndex: null,
      centerReservedRadius: circular ? projectOptionalNumber(drawing.adv.center_reserved_radius, { field: 'center_reserved_radius' }) : null
    }
  };
};

// `losatExecution` is the app-level LOSAT execution setting (`state.losatExecution`).
/**
 * @param {RequestDrawing} drawing
 * @param {{ threadsPerJob?: unknown } | undefined} losatExecution
 * @param {ProteinSettingsBaseline} [baseline]
 */
const generatedProteinSettings = (drawing, losatExecution, baseline = {}) => {
  const { alignOrthogroupFeature: _legacyAlignment, ...currentBaseline } = baseline;
  const blastp = drawing.losat.blastp || {};
  const blastpMode = requireCurrentProteinBlastpMode(blastp.mode);
  /** @param {unknown} value @param {number} fallback */
  const positiveInteger = (value, fallback) => integerSettingOr(value, fallback, 1);
  /** @param {unknown} value @param {number} fallback */
  const nonNegativeInteger = (value, fallback) => integerSettingOr(value, fallback, 0);
  const collinearityUnitMode = blastpMode === 'collinear'
    ? requireCurrentCollinearUnitMode(blastp.collinearUnitMode)
    : readCollinearUnitMode(blastp.collinearUnitMode).value;
  const collinearityColorMode = blastpMode === 'collinear'
    ? requireCurrentCollinearColorMode(blastp.collinearColorMode)
    : readCollinearColorMode(blastp.collinearColorMode).value;
  const baselineCollinearity = baseline.collinearityParams &&
    typeof baseline.collinearityParams === 'object' &&
    !Array.isArray(baseline.collinearityParams)
    ? baseline.collinearityParams
    : {};
  const baselineParameters = baselineCollinearity.parameters &&
    typeof baselineCollinearity.parameters === 'object' &&
    !Array.isArray(baselineCollinearity.parameters)
    ? baselineCollinearity.parameters
    : {};
  return {
    ...currentBaseline,
    collinearityParams: {
      ...baselineCollinearity,
      kind: baselineCollinearity.kind || 'lossless',
      parameters: {
        ...baselineParameters,
        minAnchors: blastpMode === 'collinear'
          ? requireCurrentCollinearMinAnchors(blastp.collinearMinAnchors)
          : positiveInteger(blastp.collinearMinAnchors, 1),
        maxUnitGap: blastpMode === 'collinear'
          ? requireCurrentCollinearMaxUnitGap(blastp.collinearMaxUnitGap)
          : nonNegativeInteger(blastp.collinearMaxUnitGap, 0),
        maxDiagonalDrift: blastpMode === 'collinear'
          ? requireCurrentCollinearMaxDiagonalDrift(blastp.collinearMaxDiagonalDrift)
          : nonNegativeInteger(blastp.collinearMaxDiagonalDrift, 0),
        maxConflicts: blastpMode === 'collinear'
          ? requireCurrentCollinearMaxConflicts(blastp.collinearMaxConflictsInMergeGap)
          : nonNegativeInteger(blastp.collinearMaxConflictsInMergeGap, 1),
        mergeOrientation: blastpMode === 'collinear'
          ? requireCurrentCollinearMergeOrientation(blastp.collinearMergeOrientation)
          : (baselineParameters.mergeOrientation || 'either')
      }
    },
    collinearityUnitMode,
    collinearInferOrthogroups: requireCurrentCollinearInferOrthogroups(blastp.collinearInferOrthogroups),
    collinearityAnchorMode: blastpMode === 'collinear'
      ? requireCurrentCollinearAnchorMode(blastp.collinearAnchorMode)
      : readCollinearAnchorMode(blastp.collinearAnchorMode).value,
    collinearitySearchScope: blastpMode === 'collinear'
      ? requireCurrentCollinearSearchScope(blastp.collinearSearchScope)
      : readCollinearSearchScope(blastp.collinearSearchScope).value,
    collinearityColorMode,
    losatpBin: baseline.losatpBin || 'losat',
    ncbiBlastpBin: baseline.ncbiBlastpBin ?? null,
    losatpThreads: integerSettingOr(losatExecution?.threadsPerJob, null, 1),
    proteinBlastpMaxHits: blastpMode === 'pairwise'
      ? requireCurrentProteinBlastpMaxHits(blastp.maxHits)
      : positiveInteger(blastp.maxHits, 5),
    proteinBlastpCandidateLimit: requireCurrentProteinBlastpCandidateLimit(
      blastp.candidateLimit
    ),
    orthogroupMembershipMode: ['orthogroup', 'collinear'].includes(blastpMode)
      ? requireCurrentOrthogroupMembershipMode(blastp.orthogroupMembershipMode)
      : readOrthogroupMembershipMode(blastp.orthogroupMembershipMode).value,
    orthogroupMemberMaxHits: ['orthogroup', 'collinear'].includes(blastpMode)
      ? requireCurrentOrthogroupMemberMaxHits(blastp.orthogroupMemberMaxHits)
      : integerSettingOr(blastp.orthogroupMemberMaxHits, null, 1),
    collinearMaxParalogLinksPerOrthogroup:
      blastpMode === 'collinear'
        ? requireCurrentCollinearMaxParalogLinks(
            blastp.collinearMaxParalogLinksPerOrthogroup
          )
        : positiveInteger(blastp.collinearMaxParalogLinksPerOrthogroup, 2)
  };
};

// The correction of a resolver issue the comparison panel reports (CI-06).
/** @type {Readonly<Record<string, string>>} */
const PLAN_ISSUE_REASONS = Object.freeze({
  'missing-upload': 'BLAST_TSV_REQUIRED', 'selected-losat-requires-pairwise': 'SELECTED_PAIRS_NEED_PAIRWISE'
});

/**
 * @template {{ edges?: unknown, errors?: unknown, error?: unknown }} S
 * @param {S | null | undefined} snapshot The Generate run's plan (`resolveLinearComparisonPlan`, or the
 *   snapshot its caller passed), checked here.
 * @returns {S}
 */
export const requireLinearComparisonPlanSnapshot = (snapshot) => {
  if (!snapshot || !Array.isArray(snapshot.edges)) {
    throw new Error('A resolved Linear comparison plan is required.');
  }
  const issue = Array.isArray(snapshot.errors) ? snapshot.errors[0] : null;
  if (issue || String(snapshot.error || '').trim()) {
    // A LOSATP-mode issue names the first Upload pair, when the plan has one.
    const edge = snapshot.edges.find((/** @type {LinearComparisonResolvedEdge} */ item) => item?.id && item.id === issue?.edgeId)
      || (issue?.code === 'selected-losat-requires-pairwise'
        ? snapshot.edges.find((/** @type {LinearComparisonResolvedEdge} */ item) => item?.source === 'upload') : undefined);
    throw diagnosticError('COMPARISON_INPUT', {
      reason: PLAN_ISSUE_REASONS[issue?.code] || 'PAIR_TOPOLOGY',
      queryRecordIndex: edge?.queryIndex, subjectRecordIndex: edge?.subjectIndex
    });
  }
  return snapshot;
};

/** @param {LinearComparisonResolution | null | undefined} snapshot */
const orderedComparisonPlanEdges = (snapshot) => (
  (Array.isArray(snapshot?.edges) ? snapshot.edges : [])
    .slice()
    .sort((left, right) => Number(left?.ordinal) - Number(right?.ordinal))
);

/**
 * @param {unknown} queryRecordIndex
 * @param {unknown} subjectRecordIndex
 */
const comparisonEndpointKey = (queryRecordIndex, subjectRecordIndex) => (
  `${Number(queryRecordIndex)}->${Number(subjectRecordIndex)}`
);

/**
 * The resource id of an edge's protein comparison table. The analysis that writes the table and every later
 * request that reads it back from the committed Session use this id, so the Worker keeps the table's cache token
 * and the render stages keyed by it are reused (OV-1049).
 * @param {LinearComparisonResolvedEdge} edge
 */
const proteinComparisonResourceId = (edge) => `comparison-canonical-protein-${edge.ordinal + 1}`;
// The resource ids of the typed analysis results, for the same reason: the analysis and the committed file name
// each result alike.
const ORTHOGROUP_RESULT_RESOURCE_ID = 'comparison-canonical-orthogroups-1';
const COLLINEARITY_RESULT_RESOURCE_ID = 'comparison-canonical-collinearity-1';

/**
 * @param {{
 *   comparison: { kind?: unknown, canonicalTsv?: unknown, text?: unknown },
 *   edge: LinearComparisonResolvedEdge, resources: ResourceBuilder
 * }} input
 */
const addResolvedComparisonResource = ({
  comparison,
  edge,
  resources
}) => {
  if (comparison.kind === 'precomputedProteinComparison') {
    // The table is the Python canonical TSV text of the analysis helper, written unchanged.
    const resourceId = proteinComparisonResourceId(edge);
    return {
      kind: 'precomputedProteinComparison',
      resourceId: resources.addText(resourceId, 'canonical-tsv', `${resourceId}.tsv`, comparison.canonicalTsv),
      encoding: 'canonicalTsv',
      queryRecordIndex: edge.queryIndex,
      subjectRecordIndex: edge.subjectIndex
    };
  }
  if (comparison.kind !== 'nucleotideBlast') {
    throw new Error(`Unsupported resolved comparison kind for '${edge.edgeKey}'.`);
  }
  const resourceId = `comparison-resolved-${edge.ordinal + 1}`;
  return {
    kind: 'nucleotideBlast',
    resourceId: resources.addText(
      resourceId,
      'nucleotide-blast',
      `${resourceId}.tsv`,
      String(comparison.text || '')
    ),
    queryRecordIndex: edge.queryIndex,
    subjectRecordIndex: edge.subjectIndex
  };
};

/** @param {{ comparison: RequestCanonicalComparison, edge: LinearComparisonResolvedEdge, resources: ResourceBuilder }} input */
const addPersistedProteinComparisonResource = ({
  comparison,
  edge,
  resources
}) => ({
  kind: 'precomputedProteinComparison',
  resourceId: resources.addFile(
    proteinComparisonResourceId(edge),
    canonicalComparisonResourceKind(comparison),
    comparison.file
  ),
  encoding: String(comparison.encoding || 'canonicalTsv'),
  queryRecordIndex: edge.queryIndex,
  subjectRecordIndex: edge.subjectIndex
});

/** @param {{ comparison: RequestCanonicalComparison, index: number, resources: ResourceBuilder }} input */
const addCanonicalInputProteinComparisonResource = ({
  comparison,
  index,
  resources
}) => ({
  kind: 'precomputedProteinComparison',
  resourceId: resources.addFile(
    `comparison-canonical-input-protein-${index + 1}`,
    canonicalComparisonResourceKind(comparison),
    comparison.file
  ),
  encoding: String(comparison.encoding || 'canonicalTsv'),
  queryRecordIndex: Number(comparison.queryRecordIndex),
  subjectRecordIndex: Number(comparison.subjectRecordIndex)
});

/**
 * @param {{ state: RequestState, drawing: RequestDrawing, filesData: RequestFilesData,
 *   resources: ResourceBuilder, comparisonPlanSnapshot: LinearComparisonResolution | null, resolvedComparisons?: Record<string, any>[] }} input
 */
const buildComparisons = ({
  state,
  drawing,
  filesData,
  resources,
  comparisonPlanSnapshot,
  resolvedComparisons = []
}) => {
  if (state.mode.value !== 'linear') return [];
  const snapshot = requireLinearComparisonPlanSnapshot(comparisonPlanSnapshot);
  const edges = orderedComparisonPlanEdges(snapshot);
  // The plan is the permission boundary. Do not inspect dormant uploads,
  // committed artifacts, or generated-protein metadata for an empty plan.
  if (snapshot.mode === 'none' || edges.length === 0) return [];

  const comparisons = [];
  /** @type {RequestComparisonRef[]} */
  const metadataComparisons = [];
  const uploadFilesByEdgeId = new Map(
    (Array.isArray(filesData.linearComparisons) ? filesData.linearComparisons : [])
      .filter((binding) => binding?.file)
      .map((binding) => [String(binding.id || ''), binding.file])
  );
  const resolvedByEdgeKey = new Map();
  (Array.isArray(resolvedComparisons) ? resolvedComparisons : []).forEach((comparison) => {
    const edgeKey = String(comparison?.edgeKey || '').trim();
    if (!edgeKey) return;
    if (resolvedByEdgeKey.has(edgeKey)) {
      throw new Error(`Multiple resolved comparison results were produced for '${edgeKey}'.`);
    }
    resolvedByEdgeKey.set(edgeKey, comparison);
  });
  const resolvedAnalysisArtifacts = new Map();
  (Array.isArray(resolvedComparisons) ? resolvedComparisons : []).forEach((comparison) => {
    if (!['orthogroupResult', 'collinearityResult'].includes(comparison?.kind)) return;
    if (resolvedAnalysisArtifacts.has(comparison.kind)) {
      throw new Error(`Multiple resolved ${comparison.kind} artifacts were produced.`);
    }
    const expectedValueKind = comparison.kind === 'collinearityResult'
      ? 'result'
      : 'orthogroupResult';
    const typedResource = comparison.typedResource;
    if (
      !typedResource
      || typeof typedResource !== 'object'
      || Array.isArray(typedResource)
      || typedResource.schema !== 3
      || typedResource.kind !== expectedValueKind
      || !Object.prototype.hasOwnProperty.call(typedResource, 'value')
    ) {
      throw new Error(`Resolved ${comparison.kind} metadata is not a canonical typed resource.`);
    }
    resolvedAnalysisArtifacts.set(comparison.kind, typedResource);
  });

  const persistedCanonicalComparisons = Array.isArray(filesData.linearCanonicalComparisons)
    ? filesData.linearCanonicalComparisons
    : [];
  const persistedGeneratedComparison = persistedCanonicalComparisons.find(
    (comparison) => comparison?.kind === 'generatedProteinComparison'
  );
  const activeProteinPipeline = (
    snapshot.hasLosatIntent === true && drawing.losatProgram?.value === 'blastp'
  );
  const persistedCanonicalInputs = persistedCanonicalComparisons.filter(
    (comparison) => comparison?.canonicalInput === true
  );
  const persistedMetadata = persistedCanonicalComparisons.filter((comparison) => (
    (
      (
        comparison?.kind === 'orthogroupResult'
        && String(drawing.losat?.blastp?.mode || '').trim().toLowerCase() === 'orthogroup'
      ) || (
        comparison?.kind === 'collinearityResult'
        && String(drawing.losat?.blastp?.mode || '').trim().toLowerCase() === 'collinear'
      )
    ) && activeProteinPipeline && comparison.canonicalInput !== true
  ));
  let hasResolvedProteinAnalysis = false;
  [
    ...(resolvedAnalysisArtifacts.size === 0 ? persistedCanonicalInputs : []),
    ...(resolvedAnalysisArtifacts.size === 0 ? persistedMetadata : [])
  ].forEach((comparison, index) => {
    if (!comparison.file) return;
    if (comparison.kind === 'orthogroupResult') {
      metadataComparisons.push({
        kind: 'orthogroupResult',
        resourceId: resources.addFile(
          ORTHOGROUP_RESULT_RESOURCE_ID,
          canonicalComparisonResourceKind(comparison),
          comparison.file
        ),
        encoding: String(comparison.encoding || 'canonicalJson')
      });
      if (comparison.canonicalInput !== true) hasResolvedProteinAnalysis = true;
      return;
    }
    if (comparison.kind === 'collinearityResult') {
      metadataComparisons.push({
        kind: 'collinearityResult',
        resourceId: resources.addFile(
          COLLINEARITY_RESULT_RESOURCE_ID,
          canonicalComparisonResourceKind(comparison),
          comparison.file
        ),
        encoding: String(comparison.encoding || 'canonicalJson'),
        valueKind: String(comparison.valueKind || 'result')
      });
      if (comparison.canonicalInput !== true) hasResolvedProteinAnalysis = true;
      return;
    }
    if (comparison.kind !== 'precomputedProteinComparison') return;
    const queryRecordIndex = Number(comparison.queryRecordIndex);
    const subjectRecordIndex = Number(comparison.subjectRecordIndex);
    if (
      !Number.isInteger(queryRecordIndex) ||
      !Number.isInteger(subjectRecordIndex) ||
      !filesData.linearSeqs?.[queryRecordIndex] ||
      !filesData.linearSeqs?.[subjectRecordIndex]
    ) return;
    metadataComparisons.push(addCanonicalInputProteinComparisonResource({
      comparison,
      index,
      resources
    }));
  });
  resolvedAnalysisArtifacts.forEach((typedResource, kind) => {
    const collinearity = kind === 'collinearityResult';
    const resourceId = collinearity
      ? COLLINEARITY_RESULT_RESOURCE_ID
      : ORTHOGROUP_RESULT_RESOURCE_ID;
    metadataComparisons.push({
      kind,
      resourceId: resources.addJson(
        resourceId,
        canonicalComparisonResourceKind({ kind }),
        `${resourceId}.json`,
        typedResource
      ),
      encoding: 'canonicalJson',
      ...(collinearity ? { valueKind: 'result' } : {})
    });
    hasResolvedProteinAnalysis = true;
  });

  const persistedProteinByEdgeKey = new Map();
  const persistedProteinByEndpoints = new Map();
  if (activeProteinPipeline && resolvedAnalysisArtifacts.size === 0) {
    persistedCanonicalComparisons
      .filter((comparison) => (
        comparison?.kind === 'precomputedProteinComparison' &&
        comparison.canonicalInput !== true &&
        comparison.file
      ))
      .forEach((comparison) => {
        const edgeKey = String(comparison.edgeKey || '').trim();
        if (edgeKey && !persistedProteinByEdgeKey.has(edgeKey)) {
          persistedProteinByEdgeKey.set(edgeKey, comparison);
        }
        const endpointKey = comparisonEndpointKey(
          comparison.queryRecordIndex,
          comparison.subjectRecordIndex
        );
        if (!persistedProteinByEndpoints.has(endpointKey)) {
          persistedProteinByEndpoints.set(endpointKey, comparison);
        }
      });
  }

  let hasPrecomputedProteinComparisons = false;
  for (const edge of edges) {
    if (
      !Number.isInteger(edge?.ordinal) ||
      !Number.isInteger(edge?.queryIndex) ||
      !Number.isInteger(edge?.subjectIndex) ||
      !String(edge?.edgeKey || '').trim()
    ) {
      throw new Error('The resolved Linear comparison plan contains an invalid edge.');
    }
    if (edge.source === 'upload') {
      const file = uploadFilesByEdgeId.get(String(edge.id || ''));
      if (!file) {
        throw diagnosticError('COMPARISON_INPUT', {
          reason: 'BLAST_TSV_REQUIRED', queryRecordIndex: edge.queryIndex, subjectRecordIndex: edge.subjectIndex
        });
      }
      const resourceId = `comparison-nucleotide-${edge.ordinal + 1}`;
      comparisons.push({
        kind: 'nucleotideBlast',
        resourceId: resources.addFile(resourceId, 'nucleotide-blast', file),
        queryRecordIndex: edge.queryIndex,
        subjectRecordIndex: edge.subjectIndex
      });
      continue;
    }
    if (edge.source !== 'losat') {
      throw new Error(`Unsupported comparison source for '${edge.edgeKey}'.`);
    }
    const resolved = resolvedByEdgeKey.get(edge.edgeKey);
    if (resolved) {
      const descriptor = addResolvedComparisonResource({
        comparison: resolved,
        edge,
        resources
      });
      comparisons.push(descriptor);
      if (descriptor.kind === 'precomputedProteinComparison') {
        hasPrecomputedProteinComparisons = true;
      }
      continue;
    }
    if (activeProteinPipeline) {
      const persisted = persistedProteinByEdgeKey.get(edge.edgeKey) ||
        persistedProteinByEndpoints.get(comparisonEndpointKey(edge.queryIndex, edge.subjectIndex));
      if (persisted) {
        comparisons.push(addPersistedProteinComparisonResource({
          comparison: persisted,
          edge,
          resources
        }));
        hasPrecomputedProteinComparisons = true;
      }
    }
  }

  const activeLosatEdges = edges.filter((edge) => edge.source === 'losat');
  const selectedPairwiseLosat = (
    snapshot.mode === 'selected' &&
    activeProteinPipeline &&
    String(drawing.losat?.blastp?.mode || '').trim().toLowerCase() === 'pairwise'
  );
  const shouldEmitResolvedProteinMarker = (
    activeProteinPipeline &&
    activeLosatEdges.length > 0 &&
    (hasPrecomputedProteinComparisons || hasResolvedProteinAnalysis)
  );
  const shouldGenerateSelectedProteinPairs = (
    selectedPairwiseLosat &&
    activeLosatEdges.length > 0 &&
    !hasPrecomputedProteinComparisons
  );
  const shouldGenerateAdjacentProteinPipeline = (
    snapshot.mode === 'adjacent' &&
    activeProteinPipeline &&
    activeLosatEdges.length > 0 &&
    !hasPrecomputedProteinComparisons
  );
  comparisons.push(...metadataComparisons);
  if (shouldEmitResolvedProteinMarker || shouldGenerateSelectedProteinPairs || shouldGenerateAdjacentProteinPipeline) {
    const mode = shouldEmitResolvedProteinMarker
      ? 'none'
      : selectedPairwiseLosat
        ? 'pairwise'
        : String(drawing.losat?.blastp?.mode || persistedGeneratedComparison?.mode || 'orthogroup');
    comparisons.push({
      kind: 'generatedProteinComparison',
      mode,
      pairs: mode === 'pairwise'
        ? activeLosatEdges.map((edge) => ({
            queryRecordIndex: edge.queryIndex,
            subjectRecordIndex: edge.subjectIndex
          }))
        : [],
      settings: generatedProteinSettings(
        drawing,
        state.losatExecution,
        persistedGeneratedComparison?.settings || {}
      )
    });
  }
  return comparisons;
};

/**
 * @param {RequestState & { similarityAlignmentPlan?: ExportRef<unknown>, linearRecordTranslations?: ExportRef<unknown> }} state
 * @param {RequestDrawing} drawing
 * @param {{ recordKey: string, presentation: { gridRow?: number | null } }[]} [records]
 * @param {number[] | null} [circularSourceIndexes] The file position of each drawn Circular record.
 * @returns {CanonicalRenderRequest['layout']}
 */
const buildLayout = (state, drawing, records = [], circularSourceIndexes = null) => {
  if (state.mode.value === 'linear') {
    const recordKeys = records.map((record) => requireCanonicalText(
      record.recordKey,
      'renderRequest.records[].recordKey'
    ));
    const plan = canonicalSimilarityAlignment(
      state.similarityAlignmentPlan?.value ?? null,
      recordKeys,
      'renderRequest.layout.similarityAlignment'
    );
    const translations = canonicalRecordTranslations(
      state.linearRecordTranslations?.value || [],
      recordKeys,
      'renderRequest.layout.recordTranslations',
      { requireCoverage: plan !== null }
    );
    if (!drawing.linearRecordLayoutEnabled?.value &&
        !records.some((record) => record.presentation.gridRow != null) &&
        plan === null && translations.length === 0) return {};
    return {
      recordGapPx: Math.max(0, Number(draftLinearRecordGap(drawing.linearRecordGap?.value)) || 0),
      multiRecordPositions: null,
      recordTranslations: translations,
      similarityAlignment: plan
    };
  }
  if (!drawing.form.multi_record_canvas) return {};
  // The draft keys each record's row by its file selector `#N`; the request
  // names the drawn records by their position among them, so an OFF record
  // has no token and returns to its row when it is ON again.
  /** @param {string} selector */
  const drawnPosition = (selector) => {
    const match = /^#([1-9][0-9]*)$/.exec(selector);
    if (!match || !Array.isArray(circularSourceIndexes)) return selector;
    const position = circularSourceIndexes.indexOf(Number(match[1]) - 1);
    return position < 0 ? '' : `#${position + 1}`;
  };
  const positions = Array.isArray(drawing.adv.multi_record_positions)
    ? drawing.adv.multi_record_positions
        .map((/** @type {{ selector?: unknown, row?: unknown } | null} */ entry) => {
          const selector = drawnPosition(String(entry?.selector || '').trim());
          const row = Number(entry?.row);
          return selector && Number.isInteger(row) && row > 0 ? `${selector}@${row}` : null;
        })
        .filter(Boolean)
    : [];
  return {
    multiRecordSizeMode: requireCurrentCircularMultiRecordSizeMode(
      drawing.adv.multi_record_size_mode
    ),
    // Blank is the documented default (Python's); any number is sent literally.
    multiRecordMinRadiusRatio: projectOptionalNumber(drawing.adv.multi_record_min_radius_ratio, { field: 'multi_record_min_radius_ratio' }) ?? 0.55,
    multiRecordColumnGapRatio: projectOptionalNumber(drawing.adv.multi_record_column_gap_ratio, { field: 'multi_record_column_gap_ratio' }) ?? 0.10,
    multiRecordRowGapRatio: projectOptionalNumber(drawing.adv.multi_record_row_gap_ratio, { field: 'multi_record_row_gap_ratio' }) ?? 0.05,
    multiRecordPositions: positions.length > 0 ? positions : null
  };
};

// A draft without a per-feature edit field (a Session written by the CLI, the
// Python API, or the Gallery publication) takes the request projection's
// value for it. Gallery publication and Session Load both read edits this way
// (OV-216, OV-220).
const PROJECTED_FEATURE_EDIT_FIELDS = ['featureOverrides', 'featureVisibilityManualRules', 'labelOverrideRows'];
/**
 * @param {unknown} draft A stored feature-edit draft, unvalidated; a draft of no object is read as none.
 * @param {{ semanticFeatureState?: Record<string, any> } | null | undefined} projection
 * @returns {Record<string, any>}
 */
export const featureEditsOverProjection = (draft, projection) => {
  const saved = draft && typeof draft === 'object' ? draft : {};
  const projected = projection?.semanticFeatureState || {};
  return {
    ...saved,
    ...Object.fromEntries(PROJECTED_FEATURE_EDIT_FIELDS.filter((field) => !Object.hasOwn(saved, field)
      && projected[field] !== undefined).map((field) => [field, projected[field]]))
  };
};

/** @param {unknown} value */
const publicationClone = (value) => value === undefined ? undefined : JSON.parse(JSON.stringify(value));
/**
 * @template T
 * @param {T} value
 * @returns {ExportRef<T>}
 */
const publicationRef = (value) => ({ value });
/** @param {RequestRecordView | null | undefined} record */
const canonicalRecordSelector = (record) => {
  const selector = record?.region?.selector || record?.selector;
  if (selector?.kind === 'recordId') return String(selector.value || '');
  if (selector?.kind === 'recordIndex') return `#${Number(selector.index) + 1}`;
  return null;
};
/**
 * The flat draft a Gallery Session's request is rebuilt from, as
 * `buildCanonicalRequestState` reads it: unvalidated, each leaf passed on or
 * converted where it is read.
 * @typedef {{
 *   colors?: unknown, colorsAreOverrides?: unknown, palette?: unknown,
 *   linearRecordLayout?: UnvalidatedFields<'enabled' | 'recordGap' | 'rows'>,
 *   circularConservation?: unknown, losat?: UnvalidatedFields<string>, filterMode?: unknown,
 *   blacklistText?: unknown, losatProgram?: unknown, form?: unknown, adv?: unknown, rules?: unknown,
 *   whitelist?: unknown, qualifierPriorityRules?: unknown, linearComparisonPlan?: unknown,
 *   annotationSets?: unknown, recordDisplayDrafts?: unknown, featurePlacementOverrides?: unknown
 * }} RequestStateDraftView
 */
// A Gallery Session's request inputs: the project inputs, mode and artifacts
// as `state`, the settings and edits as the drawing the request reads, and the
// per-feature edits of that drawing as the draft stores them (`features`).
/**
 * @param {{
 *   session: { features?: unknown, renderRequest?: { records?: RequestRecordView[] } | null,
 *     ui?: UnvalidatedFields<'cInputType' | 'lInputType'> } | null | undefined,
 *   projection: ReturnType<typeof projectCanonicalSessionRequest>,
 *   config: RequestStateDraftView
 * }} input The Session (version 31-33 or 39, its request read), the projection of that request, and the draft
 *   its migration rebuilt.
 * @returns {{ state: Record<string, any>, drawing: RequestDrawing, features: Record<string, unknown> }}
 */
export const buildCanonicalRequestState = ({ session, projection, config }) => {
  // The request keeps the tables it was written with, not the files bound as
  // their inputs (OV-367).
  const canonicalPublicationFiles = { ...projection.requestTableFiles };
  /** @param {unknown} colors */
  const pythonColors = (colors) => Object.fromEntries(Object.entries(colors || {}).filter(
    ([key]) => !key.startsWith('collinear_block_')).sort(([left], [right]) => left.localeCompare(right)));
  const activeColors = pythonColors(config?.colors);
  const colorOverridesChanged = config?.colorsAreOverrides === true
    && JSON.stringify(activeColors) !== JSON.stringify(pythonColors(projection.config?.colors));
  if (colorOverridesChanged) delete canonicalPublicationFiles.d_color;
  const features = featureEditsOverProjection(session?.features, projection), layout = config?.linearRecordLayout || {};
  /** @type {{ similarityAlignment?: unknown, recordTranslations?: unknown }} */
  const canonicalLayout = projection.config?.linearRecordLayout || {};
  const palette = String(config?.palette || 'default');
  const records = session?.renderRequest?.records || [];
  const circularRecordList = records.length === 1 && !records[0]?.selector
    && !records[0]?.region?.selector ? []
    : records.map((record) => ({ selector: canonicalRecordSelector(record) }));
  const conservation = publicationClone(config?.circularConservation || {
    enabled: false, reference: 'auto', labels: '', series: []
  });
  if (!String(conservation.labels || '').trim() && conservation.series?.length)
    // A row is unvalidated: one without a label, null included, gives an empty label (OV-963).
    conservation.labels = conservation.series.map((/** @type {UnvalidatedFields<'label'>} */ row) => String(row?.label || '').trim()).join(',');
  const stateRefs = {
    mode: projection.mode, cInputType: session?.ui?.cInputType || projection.inputType,
    lInputType: session?.ui?.lInputType || projection.inputType,
    circularRecordList,
    paletteDefinitions: { [palette]: {} }, editableLabels: [],
    similarityAlignmentPlan: canonicalLayout.similarityAlignment || null,
    linearRecordTranslations: publicationClone(canonicalLayout.recordTranslations || [])
  };
  const drawingRefs = {
    currentColors: colorOverridesChanged ? activeColors : config.colors || {},
    selectedPalette: palette, featureVisibilityRules: publicationClone(features.featureVisibilityManualRules || []),
    filterMode: config.filterMode || 'None', manualBlacklist: String(config.blacklistText || ''),
    canonicalLabelOverrideRows: publicationClone(features.labelOverrideRows || []),
    losatProgram: config.losatProgram || 'blastn',
    linearRecordLayoutEnabled: Boolean(layout.enabled),
    linearRecordGap: layout.recordGap ?? DEFAULT_LINEAR_RECORD_GAP_PX
  };
  /** @param {Record<string, any>} values */
  const refsOf = (values) => Object.fromEntries(Object.entries(values).map(([key, value]) => [key, publicationRef(value)]));
  // How LOSAT runs is app-level (`state.losatExecution`): the request's thread
  // count, which the projection puts in `config.losat`, over the defaults.
  const savedLosat = config?.losat && typeof config.losat === 'object' ? config.losat : {};
  const losatExecution = Object.fromEntries(Object.entries(createDefaultLosatExecution())
    .map(([field, value]) => [field, Object.hasOwn(savedLosat, field) ? savedLosat[field] : value]));
  return { features, state: { ...refsOf(stateRefs), canonicalPublicationFiles, losatExecution }, drawing: {
    ...refsOf(drawingRefs),
    form: config.form || {}, adv: config.adv || {},
    manualSpecificRules: publicationClone(config.rules || []), manualWhitelist: publicationClone(config.whitelist || []), manualPriorityRules: publicationClone(config.qualifierPriorityRules || []),
    featureOverrides: publicationClone(features.featureOverrides || {}),
    labelTextBulkOverrides: publicationClone(features.labelTextBulkOverrides || {}),
    circularConservation: conservation, losat: publicationClone(config.losat || { blastp: {} }),
    linearRecordRows: publicationClone(layout.rows || []), linearComparisonPlan: normalizeLinearComparisonPlan(config.linearComparisonPlan),
    annotationSets: publicationClone(config.annotationSets || []),
    recordDisplayDrafts: publicationClone(config.recordDisplayDrafts || []),
    featurePlacementOverrides: publicationClone(config.featurePlacementOverrides || {})
  } };
};
// `recordDisplayRows` are the record display rows of the draft sources
// (`recordDisplayControls.allRows`); the caller reads them from the owner that
// the composition root wires, never from `state` (R13).
/**
 * @param {Record<string, any>} input
 * @returns {CanonicalRenderEnvelope}
 */
const projectCanonicalRenderInput = ({
  state,
  drawing,
  // Generate supplies the label projection it already built for this operation.
  generatedLabelProjection = null,
  filesData,
  recordDisplayRows = [],
  comparisonPlanSnapshot = null,
  resolvedComparisons = [],
  resolvedCircularConservation = [],
  // The record-catalog keys of OFF records: annotations bound to them wait in the draft.
  omittedAnnotationRecordKeys = [],
  // A Save of a mode not generated yet leaves out a Depth row without a source (BA-OD114).
  omitUnsourcedDepth = false,
  resources = createResourceBuilder()
}) => {
  recordStructuralMetric('canonicalRequestProjectionCount');
  requireCurrentWebStateFieldNames(drawing);
  requireCurrentCircularMultiRecordSizeMode(drawing.adv.multi_record_size_mode);
  requireCurrentLinearTrackLayout(drawing.form.linear_track_layout);
  requireCurrentLinearLabelPlacement(drawing.adv.label_placement);
  if (
    state.mode.value === 'linear' &&
    Boolean(drawing.form.normalize_length) &&
    Boolean(drawing.linearRecordLayoutEnabled?.value) &&
    linearRecordLayoutHasSharedRow(filesData?.linearSeqs, drawing.linearRecordRows)
  ) {
    throw new Error(
      'Normalize Record Lengths cannot be used when multiple records share the same Linear row. ' +
      'Turn Normalize off or assign each record to a separate row.'
    );
  }
  const hasLinearComparisonIntent = (
    state.mode.value === 'linear' &&
    comparisonPlanSnapshot?.hasComparisonIntent === true
  );
  const comparisonOptionsRequested = (
    state.mode.value === 'circular' || hasLinearComparisonIntent
  );
  /** @type {Record<string, any>} */
  const webFiles = {};
  const recordPlan = buildRecords({ state, drawing, filesData, resources });
  const drafts = bindRecordDisplayDrafts(drawing.recordDisplayDrafts || [], recordDisplayRows);
  /** @type {number[]} */
  const sourceInputIndexes = [];
  const records = recordPlan.records.flatMap((record, index) => {
    const sourceUid = state.mode.value === 'linear'
      ? String(filesData.linearSeqs?.[index]?.uid || record.recordKey) : 'circular';
    const sourceRows = recordDisplayRows.filter((/** @type {RecordDisplayRow} */ row) => row.scope === state.mode.value && row.sourceUid === sourceUid);
    const selector = record.region?.selector || record.selector;
    const selectedRows = sourceRows.filter((/** @type {RecordDisplayRow} */ row) => !selector
      || (selector.kind === 'recordIndex' ? row.selector === `#${selector.index + 1}` : row.recordId === selector.value));
    /** @param {RecordDisplayRow} row */
    const transformFor = (row) => requestedRecordTransform(row,
      drafts.find((/** @type {RecordDisplayDraftIdentity} */ draft) => recordDisplayKey(draft) === row.key) || {}, { ...row, cropped: Boolean(record.region) });
    if (record.cardinality === 'all' && selectedRows.length > 1) {
      const transforms = selectedRows.map(transformFor);
      if (new Set(transforms.map((/** @type {unknown} */ transform) => JSON.stringify(transform))).size > 1) {
        sourceInputIndexes.push(...selectedRows.map(() => index));
        return selectedRows.map((/** @type {RecordDisplayRow} */ row, /** @type {number} */ rowIndex) => ({ ...record,
          recordKey: `${record.recordKey}:${Number(row.selector.slice(1))}`,
          cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: Number(row.selector.slice(1)) - 1 },
          presentation: { ...record.presentation,
            reverseComplement: record.region ? false : transforms[rowIndex].reverseComplement,
            gridRow: record.presentation.gridRow ?? index + 1 },
          display: transforms[rowIndex].display }));
      }
    }
    const selected = selectedRows[0];
    const savedDraft = drafts.find((/** @type {RecordDisplayDraftIdentity & Partial<RecordDisplayDraft>} */ draft) => draft.sourceUid === sourceUid
      && (selector?.kind === 'recordIndex' ? draft.selector === `#${selector.index + 1}`
        : selector?.kind === 'recordId' ? draft.recordId === selector.value : true));
    const transform = selected ? transformFor(selected) : {
      display: record.display || { isCircular: savedDraft?.topologyOverride ?? null,
        startCoordinate: savedDraft?.startCoordinate ?? null },
      reverseComplement: record.region
        ? Boolean(record.region.reverseComplement)
        : savedDraft?.reverseComplementOverride ?? Boolean(record.presentation?.reverseComplement)
    };
    sourceInputIndexes.push(index);
    return [{ ...record,
      presentation: { ...record.presentation,
        reverseComplement: record.region ? false : transform.reverseComplement },
      display: transform.display }];
  });
  if (state.mode.value === 'linear' && records.length !== recordPlan.records.length) {
    records.forEach((record, index) => { record.presentation.gridRow ??= sourceInputIndexes[index] + 1; });
  }
  if (records.length === 0) throw new Error('A canonical request requires at least one record.');
  const selectedCircularFilesData = (
    state.mode.value === 'circular' &&
    isRecordMajorDepthFileMatrix(filesData.c_depth) &&
    Array.isArray(recordPlan.circularSourceIndexes) &&
    // `circularSourceCount` is set with `circularSourceIndexes` (both are null for Linear).
    recordPlan.circularSourceIndexes.length < /** @type {number} */ (recordPlan.circularSourceCount) &&
    filesData.c_depth.length === recordPlan.circularSourceCount
  )
    ? {
        ...filesData,
        c_depth: recordPlan.circularSourceIndexes.map((sourceIndex) => (
          filesData.c_depth[sourceIndex] || []
        ))
      }
    : filesData;
  const trackPlan = buildTrackPlan({
    state,
    drawing,
    filesData: state.mode.value === 'linear'
      ? { ...filesData, linearSeqs: sourceInputIndexes.map((index) => filesData.linearSeqs[index]) }
      : selectedCircularFilesData,
    recordCount: records.length,
    resolvedCircularConservation,
    omitUnsourcedDepth
  });
  const circularGroupingIntent = ['single', 'batch'].includes(
    drawing.adv.circular_grouping_intent
  )
    ? drawing.adv.circular_grouping_intent
    : null;
  const grouping = state.mode.value === 'linear'
    ? 'single'
    : (
        drawing.form.multi_record_canvas
          ? 'grid'
          : (
              (recordPlan.circularPresentedCount ?? records.length) === 1
                ? (
                    circularGroupingIntent ||
                    WEB_UX_PROFILE.circular.singleRecordGrouping
                  )
                : WEB_UX_PROFILE.circular.multiRecordGrouping
            )
      );
  const explicitPrefix = explicitOutputPrefix(drawing.form.prefix);
  if (state.mode.value === 'circular') {
    webFiles.circularOutputPrefixExplicit = explicitPrefix !== null;
  }
  const knownCircularRecords = Array.isArray(state.circularRecordList.value)
    ? state.circularRecordList.value
    : [];
  const knownCircularEntries = buildDisambiguatedRecordEntries(
    knownCircularRecords.map((/** @type {{ record_id?: unknown, recordId?: unknown } | null} */ record) => ({
      ...record,
      recordId: record?.record_id ?? record?.recordId
    }))
  );
  const circularOutputRecords = records.map((record, index) => {
    const selector = canonicalRecordSelector(record);
    const resolved = resolveDisambiguatedRecordSelection(
      knownCircularEntries,
      selector
    );
    return {
      record_id: resolved.record?.recordId || selector || `Record_${index + 1}`
    };
  });
  const defaultCircularPrefix = downloadSafeName(circularRecordId(circularOutputRecords[0], 0), 'out');
  const output = grouping === 'batch'
    ? resolveCircularBatchPrefixes(circularOutputRecords, explicitPrefix)
        .map(renderOutputPayload)
    : renderOutputPayload(
        explicitPrefix ?? (state.mode.value === 'circular' ? defaultCircularPrefix : 'out')
      );
  const labelProjection = generatedLabelProjection || requestLabelProjection(state, drawing);
  // The option keys are not checked yet; the two keys every request has are.
  /** @type {Record<string, any> & RequestTableOptions & Pick<CanonicalRenderRequest['diagramOptions'], 'featurePlacements' | 'featureOverrides'>} */
  const diagramOptions = {
    featurePlacements: requestFeaturePlacements(drawing.featurePlacementOverrides, state.mode.value, records),
    // Per-feature edits by source identity (design Q4); the rows of other
    // records stay in the drawing's draft (R2).
    featureOverrides: requestFeatureOverrides(drawing.featureOverrides, records, {
      bulkLabelText: labelProjection.bulkLabelText
    }),
    configOverrides: buildConfigOverrides(state, drawing, {
      depthRequested: trackPlan.depthRequested,
      hasComparisonIntent: hasLinearComparisonIntent,
      linearHasSharedRow: state.mode.value === 'linear' && linearRecordLayoutHasSharedRow(
        filesData.linearSeqs,
        drawing.linearRecordRows,
        { enabled: Boolean(drawing.linearRecordLayoutEnabled?.value) }
      )
    }),
    tracks: trackPlan.tracks,
    output: {
      legend: legendPositionOfMode(drawing.form.legend, state.mode.value)
        ?? WEB_UX_PROFILE[state.mode.value === 'linear' ? 'linear' : 'circular'].legend,
      plotTitlePosition: String(drawing.adv.plot_title_position || (state.mode.value === 'linear' ? 'bottom' : 'none'))
    },
    selectedFeaturesSet: Array.from(drawing.adv.features || []).map((value) => String(value)),
    featureShapes: {
      repeat_region: defaultFeatureRendering('repeat_region'),
      ...normalizeFeatureRenderingMap(drawing.adv.feature_shapes || {})
    },
    dinucleotide: String(drawing.adv.nt || 'GC').toUpperCase(),
    window: projectOptionalNumber(drawing.adv.window_size, { field: 'window' }),
    step: projectOptionalNumber(drawing.adv.step_size, { field: 'step' }),
    depthWindow: projectOptionalNumber(drawing.adv.depth_window_size, { field: 'depth_window' }),
    depthStep: projectOptionalNumber(drawing.adv.depth_step_size, { field: 'depth_step' }),
    plotTitle: String(drawing.form.plot_title || '').trim() || null,
    plotTitleFontSize: projectOptionalNumber(drawing.adv.plot_title_font_size, { field: 'plot_title_font_size' }),
    ...(comparisonOptionsRequested
      ? projectComparisonThresholds(resolveComparisonThresholds(drawing.adv, state.mode.value))
      : {})
  };
  if (Array.isArray(drawing.annotationSets) && drawing.annotationSets.length > 0) {
    diagramOptions.annotations = annotationOptionsPayload(drawing.annotationSets, records, omittedAnnotationRecordKeys);
  }
  if (state.mode.value === 'circular') {
    diagramOptions.keepFullDefinitionWithPlotTitle = Boolean(drawing.adv.keep_full_definition_with_plot_title);
    diagramOptions.species = String(drawing.form.species || '').trim() || null;
    diagramOptions.strain = String(drawing.form.strain || '').trim() || null;
    const conservationSource = String(
      drawing.circularConservation.source || ''
    ).trim().toLowerCase();
    const conservationBlastsAreDerived = (
      filesData.c_conservation_blasts_source === 'losat-cache'
    );
    const conservation = (
      conservationSource === 'upload' || conservationBlastsAreDerived
    ) && Array.isArray(filesData.c_conservation_blasts)
      ? filesData.c_conservation_blasts
      : [];
    if (conservation.length > 0) {
      const conservationEntries = orderedConservationSources(
        conservation,
        drawing.circularConservation
      );
      diagramOptions.conservationBlastFiles = conservationEntries.map((entry, index) => fileRef(
        resources.addFile(
          `conservation-blast-files-${index + 1}`,
          'conservation-blast-file',
          entry.file
        )
      ));
      const comparisonSources = Array.isArray(filesData.c_conservation_sequence_sources)
        ? filesData.c_conservation_sequence_sources
        : [];
      const orderedComparisonSources = conservationEntries.map(
        (entry) => comparisonSources[entry.sourceIndex] || null
      );
      if (orderedComparisonSources.some(Boolean)) {
        diagramOptions.conservationFastaFiles = orderedComparisonSources.map((entry, index) => (
          entry
            ? fileRef(resources.addFile(
              `conservation-fasta-files-${index + 1}`,
              'conservation-fasta-file',
              entry
            ))
            : null
        ));
      }
      Object.assign(diagramOptions, conservationDiagramOptions(drawing.circularConservation, conservationEntries, 'auto'));
      if (conservationBlastsAreDerived) {
        webFiles.conservationBlastSource = 'losat-cache';
      }
    }
    if (conservationSource === 'losat') {
      const comparisonFastas = orderedOptionalConservationFiles(
        filesData.c_conservation_fastas,
        drawing.circularConservation
      );
      if (comparisonFastas.length > 0) {
        webFiles.conservationLosatFastaSources = comparisonFastas.map(
          (entry, index) => (
            entry
              ? resources.addFile(
                  `conservation-losat-fasta-files-${index + 1}`,
                  'conservation-fasta-file',
                  entry
                )
              : null
          )
        );
      }
    }
    if (
      Array.isArray(resolvedCircularConservation) &&
      resolvedCircularConservation.length > 0
    ) {
      diagramOptions.conservationBlastFiles = resolvedCircularConservation.map(
        (entry, index) => fileRef(resources.addText(
          `conservation-resolved-blast-${index + 1}`,
          'conservation-blast-file',
          String(entry?.name || `comparison-${index + 1}.tsv`),
          String(entry?.text || '')
        ))
      );
      const comparisonFastas = resolvedCircularConservation.map(
        (entry, index) => (
          entry?.fasta
            ? fileRef(resources.addFile(
                `conservation-resolved-fasta-${index + 1}`,
                'conservation-fasta-file',
                entry.fasta
              ))
            : null
        )
      );
      if (comparisonFastas.some(Boolean)) {
        diagramOptions.conservationFastaFiles = comparisonFastas;
      }
      Object.assign(diagramOptions, conservationDiagramOptions(drawing.circularConservation,
        resolvedCircularConservation.map((entry, index) => ({
          label: String(entry?.label || `Comparison ${index + 1}`),
          color: String(entry?.color || '#D9EAF7')
        })), 'subject'));
      webFiles.conservationBlastSource = 'losat-cache';
    }
  } else if (hasLinearComparisonIntent) {
    diagramOptions.pairwiseMatchStyle = String(drawing.adv.pairwise_match_style || 'ribbon');
  }

  addGeneratedTableResources(state, drawing, resources, diagramOptions, ALL_GENERATED_TABLES, labelProjection);
  if (trackPlan.depthRequested) {
    buildDepthResources({
      state,
      drawing,
      filesData: state.mode.value === 'linear'
      ? { ...filesData, linearSeqs: sourceInputIndexes.map((index) => filesData.linearSeqs[index]) }
      : selectedCircularFilesData,
      resources,
      diagramOptions,
      recordCount: records.length
    });
  }

  [
    ['colors-default-colors-file', filesData.d_color],
    ['colors-color-table-file', filesData.t_color],
    ['label-whitelist-file', filesData.whitelist],
    ['qualifier-priority-file', filesData.qualifier_priority]
  ].forEach(([resourceId, entry]) => {
    if (!resources.resources[resourceId]) return;
    const originalName = normalizeOriginalResourceName(entry?.name);
    if (originalName) resources.resourceOriginalNames[resourceId] = originalName;
  });

  if (Object.keys(resources.resourceOriginalNames).length > 0) {
    webFiles.resourceOriginalNames = { ...resources.resourceOriginalNames };
  }
  if (state.mode.value === 'circular' && state.cInputType.value === 'gb') {
    const circularInputOriginalName = normalizeOriginalResourceName(filesData.c_gb?.name);
    if (circularInputOriginalName) webFiles.circularInputOriginalName = circularInputOriginalName;
  }
  if (state.mode.value === 'linear') {
    webFiles.linearRecordMetadata = sourceInputIndexes.map((sourceIndex, index) => {
      /** @type {{ recordKey: string, losatGencode: number, fileDefinition?: string, fileSubtitle?: string, recordDefinition?: string, recordSubtitle?: string }} */
      const entry = {
        recordKey: String(records[index]?.recordKey || filesData.linearSeqs[sourceIndex]?.uid || `record-${index + 1}`),
        losatGencode: integerSettingOr(filesData.linearSeqs[sourceIndex]?.losat_gencode, 1, 1)
      };
      const fileDefinition = String(filesData.linearSeqs[sourceIndex]?.file_definition || '').trim();
      const fileSubtitle = String(filesData.linearSeqs[sourceIndex]?.file_subtitle || '').trim();
      if (fileDefinition) entry.fileDefinition = fileDefinition;
      if (fileSubtitle) entry.fileSubtitle = fileSubtitle;
      // presentation.label and presentation.subtitle carry the resolved value.
      // Record the per-record value too, so an override that happens to equal
      // the file default is restored as an override rather than as inheritance.
      const recordDefinition = String(filesData.linearSeqs[sourceIndex]?.definition || '');
      const recordSubtitle = String(filesData.linearSeqs[sourceIndex]?.record_subtitle || '');
      if (fileDefinition) entry.recordDefinition = recordDefinition;
      if (fileSubtitle) entry.recordSubtitle = recordSubtitle;
      return entry;
    });
  }

  return /** @satisfies {CanonicalRenderEnvelope} */ ({
    // The mode decides the options, the layout, the comparisons and the output
    // together; tsc cannot see that the builders above agree on it.
    renderRequest: /** @type {CanonicalRenderRequest} */ ({
      schema: CANONICAL_REQUEST_SCHEMA,
      mode: state.mode.value,
      grouping,
      records,
      diagramOptions,
      layout: buildLayout(state, drawing, records, recordPlan.circularSourceIndexes),
      comparisons: buildComparisons({
      state,
      drawing,
      filesData,
      resources,
      comparisonPlanSnapshot,
      resolvedComparisons
    }),
      output
    }),
    resources: resources.resources,
    webFiles
  });
};

// The request writer and the inexpensive comparison use the same projection.
// The latter retains bindings and generated table text, never resource bytes.
/**
 * @param {unknown} resource A committed resource entry (unvalidated).
 * @returns {{ text: string, bindings?: undefined } | { bindings: unknown[] | null } | null}
 */
const generationResourceIdentity = (resource) => {
  if (!resource) return null;
  // Every key is a descriptor object the writer made, so a value that is not an object has no entry.
  const generated = typeof resource === 'object' ? generatedResourceValues.get(resource) : undefined;
  if (generated !== undefined) {
    return typeof generated === 'string' ? { text: generated } : generated;
  }
  // A view: the reads below test each member they use.
  const owner = /** @type {UnvalidatedFields<'arrayBuffer' | 'data' | 'encoding'>} */ (getResourcePayloadOwner(resource));
  const backing = /** @type {SessionResourceSource | null} */ (getSessionResourceSource(owner));
  if (backing?.descriptors) {
    return { bindings: backing.descriptors.map((entry) => entry.descriptor) };
  }
  if (backing?.descriptor) return { bindings: [backing.descriptor] };
  if (typeof owner?.arrayBuffer === 'function'
    || (typeof owner?.data === 'string' && owner?.encoding)) {
    return { bindings: [owner] };
  }
  return { bindings: null };
};

// Successful render/import already validates these backings. Compare immutable
// owners or their encoded payloads, without another genome read or digest.
/**
 * A record source by field: a resource reference field holds the parts of its
 * resource (null when it has none), any other field its value.
 * @typedef {{ payload: unknown, size: unknown, encoding: unknown }} CompositionSourcePart
 * @typedef {Array<[string, unknown] | [string, CompositionSourcePart[] | null]>} CompositionSourceIdentity
 */
/**
 * @param {Readonly<Record<string, unknown>> | null | undefined} source A record's `source`.
 * @param {Readonly<Record<string, unknown>> | undefined} resources
 */
const compositionSourceIdentity = (source, resources) => {
  if (!source || typeof source !== 'object') return null;
  const fields = Object.keys(source).sort();
  /** @type {CompositionSourceIdentity} */
  const identity = fields.map((key) => {
    if (!isCanonicalResourceReferenceField(key)) return [key, source[key]];
    // A reference field of a request its reader read holds a resource id.
    const descriptor = resources?.[/** @type {string} */ (source[key])];
    if (!descriptor) return [key, null];
    const parts = generationResourceIdentity(descriptor)?.bindings;
    if (!Array.isArray(parts) || !parts.length) return [key, null];
    return [key, parts.map((entry) => {
      const part = /** @type {UnvalidatedFields<'data' | 'size' | 'encoding'>} */ (entry);
      return {
        payload: typeof part?.data === 'string' ? part.data : part,
        size: part?.size,
        encoding: part?.encoding || 'file'
      };
    })];
  });
  return identity;
};

/**
 * @param {CommittedCanonicalSession | null | undefined} canonical A committed Session, or a candidate's.
 * @param {unknown} keys The record keys of one Result, unchecked; checked here.
 */
export const projectCompositionRecordIdentity = (canonical, keys) => {
  const request = canonical?.renderRequest;
  if (!request || !Array.isArray(keys) || !keys.length || new Set(keys).size !== keys.length) return null;
  const records = [];
  for (const key of [...keys].sort()) {
    // The renderer's validated catalog expands an ALL source as
    // <recordKey>:<one-based biological source selector>. This is a source
    // record binding, never a Result index or display order.
    const matches = (request.records || []).flatMap(record => {
      if (record.recordKey === key) {
        // With no region/selector, successful EXACTLY_ONE selects the sole
        // source record; FIRST selects #1, and unexpanded ALL also proves one.
        const selector = record.selector ?? (!record.region && ['exactly_one', 'first', 'all'].includes(record.cardinality)
          ? { kind: /** @type {const} */ ('recordIndex'), index: 0 } : null);
        return [{ record, selector }];
      }
      const suffix = key.startsWith(`${record.recordKey}:`) ? key.slice(record.recordKey.length + 1) : '';
      if (record.cardinality !== 'all' || !/^[1-9]\d*$/.test(suffix)) return [];
      return [{ record, selector: { kind: /** @type {const} */ ('recordIndex'), index: Number(suffix) - 1 } }];
    });
    if (matches.length !== 1) return null;
    const { record, selector } = matches[0];
    const source = compositionSourceIdentity(record.source, canonical.resources);
    if (!source || source.some(([field, value]) => isCanonicalResourceReferenceField(field)
      // A reference field holds its parts (`compositionSourceIdentity`).
      && (!value || /** @type {CompositionSourcePart[]} */ (value).some(part => !part.payload || !Number.isSafeInteger(part.size))))) return null;
    records.push({ key, source, selector: publicationClone(selector),
      region: publicationClone(record.region ?? null) });
  }
  return { mode: request.mode, grouping: request.grouping, records };
};

/** @type {typeof projectCanonicalRenderInput} */
export const buildCanonicalRenderRequest = (args) => projectCanonicalRenderInput(args);

/**
 * @param {{ source?: { kind?: unknown, resourceId?: unknown, gffResourceId?: unknown, fastaResourceId?: unknown } | null } | null | undefined} record
 * @param {string} field
 */
const recordSourceResourceId = (record, field) => {
  const source = record?.source || {};
  if (field === 'gb' && source.kind === 'genbank') return String(source.resourceId || '');
  if (field === 'gff' && source.kind === 'gffFasta') return String(source.gffResourceId || '');
  if (field === 'fasta' && source.kind === 'gffFasta') return String(source.fastaResourceId || '');
  return '';
};

/** @param {{ resourceId?: unknown } | null | undefined} ref */
const referencedResourceId = (ref) => String(ref?.resourceId || '').trim();

/**
 * @param {Record<string, string>} target
 * @param {unknown} resourceId
 * @param {unknown} name
 */
const addResourceOriginalNameHint = (target, resourceId, name) => {
  const id = String(resourceId || '').trim();
  const originalName = normalizeOriginalResourceName(name);
  if (!id || !originalName || Object.prototype.hasOwnProperty.call(target, id)) return;
  target[id] = originalName;
};

/**
 * The projection's unvalidated Session parts, read through views.
 * @param {{
 *   renderRequest: UnvalidatedFields<'records' | 'mode'> & UnvalidatedFields<'diagramOptions', UnvalidatedFields<'colors',
 *     UnvalidatedFields<'defaultColorsFile' | 'defaultColors' | 'colorTableFile' | 'colorTable', UnvalidatedFields<'resourceId'>>>>,
 *   legacyFiles: unknown,
 *   fileBindings: unknown
 * }} input
 */
const legacyResourceOriginalNames = ({ renderRequest, legacyFiles, fileBindings }) => {
  /** @type {Record<string, string>} */
  const hints = {};
  /** @type {readonly UnvalidatedFields<'source', UnvalidatedFields<'kind' | 'resourceId' | 'gffResourceId' | 'fastaResourceId'>>[]} */
  const records = Array.isArray(renderRequest?.records) ? renderRequest.records : [];
  const options = renderRequest?.diagramOptions || {};
  // The cast: an object other than an array, read through views of its properties.
  const files = legacyFiles && typeof legacyFiles === 'object' && !Array.isArray(legacyFiles)
    ? /** @type {Readonly<Record<string, UnvalidatedFields<'name'>>>} */ (legacyFiles)
    : {};
  const namedOptionResources = {
    d_color: referencedResourceId(options.colors?.defaultColorsFile || options.colors?.defaultColors),
    t_color: referencedResourceId(options.colors?.colorTableFile || options.colors?.colorTable),
    whitelist: referencedResourceId(editorTableRef(options, 'whitelist')),
    qualifier_priority: referencedResourceId(editorTableRef(options, 'priority'))
  };
  Object.entries(namedOptionResources).forEach(([slot, resourceId]) => {
    addResourceOriginalNameHint(hints, resourceId, files?.[slot]?.name);
  });

  if (renderRequest?.mode === 'linear') {
    /** @type {readonly UnvalidatedFields<'gb' | 'gff' | 'fasta', UnvalidatedFields<'name'>>[]} */
    const sequences = Array.isArray(files.linearSeqs) ? files.linearSeqs : [];
    records.forEach((record, index) => {
      const sequence = sequences[index] || {};
      /** @type {const} */ (['gb', 'gff', 'fasta']).forEach((field) => {
        addResourceOriginalNameHint(
          hints,
          recordSourceResourceId(record, field),
          sequence?.[field]?.name
        );
      });
    });
  } else {
    const record = records[0];
    ['gb', 'gff', 'fasta'].forEach((field) => {
      addResourceOriginalNameHint(
        hints,
        recordSourceResourceId(record, field),
        files?.[`c_${field}`]?.name
      );
    });
  }

  (Array.isArray(fileBindings) ? fileBindings : []).forEach((/** @type {UnvalidatedFields<'slot' | 'name'>} */ binding) => {
    const slot = String(binding?.slot || '');
    const normalizedSlot = slot.replace(/^files\./, '');
    if (Object.prototype.hasOwnProperty.call(namedOptionResources, normalizedSlot)) {
      addResourceOriginalNameHint(
        hints,
        // An own key of the table, checked above.
        namedOptionResources[/** @type {keyof typeof namedOptionResources} */ (normalizedSlot)],
        binding?.name
      );
      return;
    }
    const linearMatch = slot.match(/^(?:files\.)?linearSeqs\[(\d+)\]\.(gb|gff|fasta)$/);
    if (linearMatch) {
      const record = records[Number(linearMatch[1])];
      addResourceOriginalNameHint(
        hints,
        recordSourceResourceId(record, linearMatch[2]),
        binding?.name
      );
      return;
    }
    const circularMatch = slot.match(/^(?:files\.)?c_(gb|gff|fasta)$/);
    if (circularMatch) {
      addResourceOriginalNameHint(
        hints,
        recordSourceResourceId(records[0], circularMatch[1]),
        binding?.name
      );
    }
  });

  return hints;
};

/**
 * @param {StoredResources} resources
 * @param {Readonly<Record<string, unknown>>} originalNameHints
 */
const resourcesWithOriginalNames = (resources, originalNameHints) => {
  const hints = originalNameHints && typeof originalNameHints === 'object' && !Array.isArray(originalNameHints)
    ? originalNameHints
    : {};
  return Object.fromEntries(Object.entries(resources || {}).map(([resourceId, entry]) => {
    if (!entry || typeof entry !== 'object' || Array.isArray(entry)) return [resourceId, entry];
    const storedName = normalizeOriginalResourceName(entry.name);
    const prefix = `${resourceId}-`;
    let inferredName = storedName;
    while (inferredName.startsWith(prefix) && inferredName.length > prefix.length) {
      inferredName = inferredName.slice(prefix.length);
    }
    const originalName = normalizeOriginalResourceName(hints[resourceId]) || inferredName;
    return [resourceId, originalName && originalName !== entry.name ? { ...entry, name: originalName } : entry];
  }));
};

/**
 * The resources of a Session, as its file stores them: no reader has checked a descriptor.
 * @typedef {Readonly<Record<string, RawSessionResourceDescriptor>> | null | undefined} StoredResources
 */
// A resource id is Session text, so only an own key of the table names a
// resource (OV-613).
/**
 * @param {StoredResources} resources
 * @param {unknown} resourceId A request's resource id, or one a stored Web inventory names (unvalidated).
 */
const resourceAsLegacyFile = (resources, resourceId) => {
  // `String` is the property key the lookup used before: an id is text or a JSON value.
  const entry = ownEntry(resources, String(resourceId));
  if (!entry || typeof entry !== 'object') throw new Error(`Missing canonical resource: ${resourceId}`);
  const { kind: _kind, ...file } = entry;
  return file;
};

/**
 * @typedef {(resourceId: unknown, metadata?: SessionResourceFileMetadata) => unknown} ResourceFileResolver
 *   Builds the file a Session binding stands for from its resource id and the
 *   binding's recorded name, type, and modification time.
 */

/**
 * @param {StoredResources} resources
 * @param {UnvalidatedFields<'resourceId' | 'name' | 'type' | 'lastModified'>} binding
 * @param {ResourceFileResolver | null} [resolveResourceFile]
 */
const webBindingAsLegacyFile = (resources, binding, resolveResourceFile = null, schema = 1) => {
  if (binding === null || binding === undefined) return null;
  const resourceId = String(binding.resourceId || '').trim();
  const metadata = schema === 2
    ? { name: binding.name, type: binding.type, lastModified: binding.lastModified }
    : {
        name: normalizeOriginalResourceName(binding.name),
        type: String(binding.type || ''),
        lastModified: Number(binding.lastModified) || 0
      };
  if (resolveResourceFile) return resolveResourceFile(resourceId, metadata);
  const file = resourceAsLegacyFile(resources, resourceId);
  return { ...file, ...metadata, name: schema === 2 ? metadata.name : (metadata.name || file.name) };
};

/**
 * @param {StoredResources} resources
 * @param {unknown} value A binding, or a list of bindings at any depth.
 * @param {ResourceFileResolver | null} [resolveResourceFile]
 * @returns {unknown}
 */
const webBindingValueAsLegacyFile = (resources, value, resolveResourceFile = null, schema = 1) => (
  Array.isArray(value)
    ? /** @type {unknown[]} */ (value).map((item) => webBindingValueAsLegacyFile(resources, item, resolveResourceFile, schema))
    : webBindingAsLegacyFile(
      resources, /** @type {UnvalidatedFields<'resourceId' | 'name' | 'type' | 'lastModified'>} */ (value),
      resolveResourceFile, schema
    )
);

/**
 * @param {Record<string, any>} files
 * @param {{ readonly bindings?: WebFileBindings } | null | undefined} webMetadata `webFiles`, its bindings checked
 * @param {Record<string, any> | null | undefined} resources
 * @param {{
 *   resolveResourceFile?: ResourceFileResolver | null,
 *   sessionResourceTable?: SessionResourceTable | null,
 *   adoptCanonicalPayloads?: boolean
 * }} [options]
 */
const applyWebFileBindings = (
  files,
  webMetadata,
  resources,
  { resolveResourceFile = null, sessionResourceTable = null, adoptCanonicalPayloads = false } = {}
) => {
  const bindings = webMetadata?.bindings;
  if (bindings === undefined) return files;
  /** @type {ResourceFileResolver | null} */
  const resolveBinding = bindings.schema === 2 && sessionResourceTable
    ? (id, metadata) => createSessionResourceFileView(sessionResourceTable, id, metadata)
    : resolveResourceFile;
  const restored = { ...files };
  [
    'c_gb',
    'c_gff',
    'c_fasta',
    'c_depth',
    'c_conservation_blasts',
    'c_conservation_fastas',
    'c_conservation_sequence_sources',
    'd_color',
    't_color',
    'blacklist',
    'whitelist',
    'qualifier_priority'
  ].forEach((field) => {
    if (!Object.prototype.hasOwnProperty.call(bindings, field)) return;
    if (field === 'c_gb' && bindings.c_gb?.kind === 'composite') {
      const table = sessionResourceTable || adoptCurrentSessionResources(resources);
      restored.c_gb = createCombinedSessionResourceFileView(table, bindings.c_gb.components, bindings.c_gb);
      return;
    }
    restored[field] = webBindingValueAsLegacyFile(
      resources,
      bindings[field],
      resolveBinding,
      bindings.schema
    );
  });
  restored.c_conservation_blasts_source =
    bindings.c_conservation_blasts_source === 'losat-cache' ? 'losat-cache' : null;

  if (Array.isArray(bindings.linearSeqs)) {
    restored.linearSeqs = bindings.linearSeqs.map((/** @type {WebFileBindingRow} */ sequence, /** @type {number} */ index) => ({
      uid: String(sequence?.uid || `canonical-seq-${index + 1}`),
      gb: webBindingValueAsLegacyFile(resources, sequence?.gb, resolveBinding, bindings.schema),
      gff: webBindingValueAsLegacyFile(resources, sequence?.gff, resolveBinding, bindings.schema),
      fasta: webBindingValueAsLegacyFile(resources, sequence?.fasta, resolveBinding, bindings.schema),
      depth: webBindingValueAsLegacyFile(resources, sequence?.depth, resolveBinding, bindings.schema),
      blast: webBindingValueAsLegacyFile(resources, sequence?.blast, resolveBinding, bindings.schema),
      losat_filename: String(sequence?.losat_filename || ''),
      // As stored: the Files read them once, with a row for what they cannot restore (OV-568).
      ...Object.fromEntries(LINEAR_SEQ_READ_FIELDS.map((field) => [field, sequence?.[field]]))
    }));
  }
  if (Array.isArray(bindings.linearComparisons)) {
    restored.linearComparisons = bindings.linearComparisons.map((/** @type {WebFileBindingRow} */ comparison, /** @type {number} */ index) => ({
      ...(adoptCanonicalPayloads ? comparison : cloneCanonicalJsonValue(comparison)),
      id: String(comparison?.id || `linear-comparison-restored-${index + 1}`),
      queryUid: String(comparison?.queryUid || ''),
      subjectUid: String(comparison?.subjectUid || ''),
      source: String(comparison?.source || 'upload'),
      file: webBindingValueAsLegacyFile(resources, comparison?.file, resolveBinding, bindings.schema)
    }));
  }
  if (Array.isArray(bindings.linearCanonicalComparisons)) {
    restored.linearCanonicalComparisons = bindings.linearCanonicalComparisons.map((/** @type {WebFileBindingRow} */ comparison) => ({
      ...(adoptCanonicalPayloads ? comparison : cloneCanonicalJsonValue(comparison)),
      file: webBindingValueAsLegacyFile(resources, comparison?.file, resolveBinding, bindings.schema)
    }));
  }
  return restored;
};

// Called only after Session authority has admitted the explicit source-free
// document. No canonical request or committed owner is created by this projection.
// A settings-only Session's draft: the shown mode's slice (Session 46), or
// the flat draft of an older Session.
/**
 * What `projectSettingsOnlySession` reads of the admitted document: the shown
 * mode Load read (K19), the drafts as stored, and the checked Web inventory.
 * @typedef {{
 *   ui: { mode: 'circular' | 'linear', lInputType?: unknown, cInputType?: unknown },
 *   modes?: UnvalidatedFields<'circular' | 'linear', UnvalidatedFields<'config'>>, config?: unknown,
 *   webFiles?: { readonly bindings?: WebFileBindings } | null, resources?: Readonly<Record<string, unknown>> | null
 * }} SettingsOnlyProjectionInput
 */
/**
 * @param {SettingsOnlyProjectionInput} data
 * @param {SessionResourceTable | null} [sessionResourceTable]
 */
export const projectSettingsOnlySession = (data, sessionResourceTable) => ({
  mode: data.ui.mode,
  inputType: data.ui.mode === 'linear' ? data.ui.lInputType : data.ui.cInputType,
  config: cloneCanonicalJsonValue(data.modes && typeof data.modes === 'object'
    ? data.modes[data.ui.mode]?.config
    : data.config),
  files: applyWebFileBindings({}, data.webFiles, data.resources, {
    sessionResourceTable, adoptCanonicalPayloads: true
  }),
  semanticFeatureState: {}
});

/**
 * A request projection as Load reads it: a settings-only Session's has no
 * request, so it has none of the request's tables, layout preferences, legacy
 * table rows, or protein pipeline.
 * @typedef {ReturnType<typeof projectCanonicalSessionRequest>
 *   | (ReturnType<typeof projectSettingsOnlySession> & { requestTableFiles?: undefined, layoutPreferences?: undefined,
 *     legacyTableRepairs?: undefined, pipelineState?: undefined })} SessionRequestProjection
 */

/** @param {unknown} value */
const cloneCanonicalJsonValue = (value) => (
  value === undefined ? undefined : JSON.parse(JSON.stringify(value))
);

/**
 * A record of the request `projectCanonicalSessionRequest` projects. Its
 * reader (`readStoredRenderRequest`) or the writer made the request, so a
 * record is an object; the projection reads each field as unvalidated and
 * checks or converts it where it reads it.
 * @typedef {RequestRecordView & {
 *   recordKey?: unknown, cardinality?: unknown,
 *   source?: UnvalidatedFields<'kind' | 'resourceId' | 'gffResourceId' | 'fastaResourceId'>,
 *   display?: UnvalidatedFields<'isCircular' | 'startCoordinate'>,
 *   presentation?: UnvalidatedFields<'reverseComplement' | 'label' | 'subtitle' | 'gridRow' | 'gridColumn'>
 * }} ProjectionRecordView
 */
// Web controls own row membership and record order, not a distinct numeric
// column value. Preserve the engine's row/column render order in the projected
// record sequence, then retire the unsupported numeric column in the draft.
/**
 * @template {{ presentation?: UnvalidatedFields<'gridRow' | 'gridColumn'> }} R
 * @param {readonly R[]} [records]
 */
export const normalizeWebGridColumnOrdering = (records = []) => {
  const source = Array.isArray(records) ? records : [];
  const entries = source.map((record, sourceIndex) => {
    const rawRow = record?.presentation?.gridRow;
    const rawColumn = record?.presentation?.gridColumn;
    const row = Number(rawRow);
    const column = Number(rawColumn);
    return {
      record,
      sourceIndex,
      row: rawRow !== null && rawRow !== undefined && Number.isInteger(row) ? row : sourceIndex + 1,
      column: rawColumn !== null && rawColumn !== undefined && Number.isInteger(column)
        ? column
        : sourceIndex + 1
    };
  });
  const hasColumns = source.some((record) => (
    record?.presentation?.gridColumn !== null
    && record?.presentation?.gridColumn !== undefined
  ));
  const ordered = hasColumns
    ? entries.slice().sort((left, right) => (
        left.row - right.row
        || left.column - right.column
        || left.sourceIndex - right.sourceIndex
      ))
    : entries;
  const projectedIndexBySourceIndex = new Map(
    ordered.map((entry, projectedIndex) => [entry.sourceIndex, projectedIndex])
  );
  return {
    records: ordered.map(({ record }) => ({
      ...record,
      presentation: {
        ...(record?.presentation || {}),
        gridColumn: null
      }
    })),
    sourceIndexByProjectedIndex: ordered.map((entry) => entry.sourceIndex),
    projectedIndexBySourceIndex
  };
};

/** @param {{ collinearInferOrthogroups?: unknown }} settings @param {string} mode */
const resolvePipelineCollinearInference = (settings, mode) =>
  settings.collinearInferOrthogroups ?? ['collinear', 'none'].includes(mode);

/**
 * @param {UnvalidatedFields<'kind' | 'mode' | 'settings'>} comparison A request's comparison (unvalidated).
 * @param {{ adoptCanonicalPayloads?: boolean }} [options]
 */
const projectGeneratedProteinPipeline = (
  comparison,
  { adoptCanonicalPayloads = false } = {}
) => {
  if (
    !comparison ||
    comparison.kind !== 'generatedProteinComparison' ||
    !comparison.settings ||
    typeof comparison.settings !== 'object' ||
    Array.isArray(comparison.settings)
  ) return null;
  // The cast: an object other than an array (checked above), read through views.
  const settings = /** @type {Readonly<Record<string, unknown>> & {
    readonly collinearInferOrthogroups?: unknown,
    readonly collinearityParams?: UnvalidatedFields<'parameters', UnvalidatedFields<
      'minAnchors' | 'maxUnitGap' | 'maxDiagonalDrift' | 'maxConflicts' | 'mergeOrientation'>>
  }} */ (comparison.settings);
  if (Object.hasOwn(settings, 'alignOrthogroupFeature')) {
    throw new Error('Current canonical protein settings cannot contain legacy alignment state.');
  }
  const parameters = settings.collinearityParams?.parameters || {};
  const mode = String(comparison.mode || 'orthogroup');
  return {
    generatedProteinComparison: adoptCanonicalPayloads
      ? comparison
      : cloneCanonicalJsonValue(comparison),
    config: {
      blastSource: 'losat',
      losatProgram: 'blastp',
      losat: {
        threadsPerJob: settings.losatpThreads ?? 'auto',
        blastp: {
          // "none" means that the canonical request is reusing resolved
          // artifacts. It does not identify the Web generation mode that
          // should be selected for a future rerun, so leave the saved UI
          // setting authoritative in that case.
          ...(mode === 'none' ? {} : { mode }),
          maxHits: settings.proteinBlastpMaxHits,
          candidateLimit: settings.proteinBlastpCandidateLimit,
          orthogroupMembershipMode: settings.orthogroupMembershipMode,
          orthogroupMemberMaxHits: settings.orthogroupMemberMaxHits,
          collinearMinAnchors: parameters.minAnchors,
          collinearMaxUnitGap: parameters.maxUnitGap,
          collinearMaxDiagonalDrift: parameters.maxDiagonalDrift,
          collinearMaxConflictsInMergeGap: parameters.maxConflicts,
          collinearMaxParalogLinksPerOrthogroup:
            settings.collinearMaxParalogLinksPerOrthogroup,
          collinearColorMode: settings.collinearityColorMode,
          collinearUnitMode: settings.collinearityUnitMode,
          collinearAnchorMode: settings.collinearityAnchorMode,
          collinearMergeOrientation: parameters.mergeOrientation,
          collinearInferOrthogroups: resolvePipelineCollinearInference(settings, mode),
          collinearSearchScope: settings.collinearitySearchScope
        }
      }
    }
  };
};

const LEGACY_DEPTH_OPTION_FIELDS = Object.freeze([
  'depthTable',
  'depthFile',
  'depthTables',
  'depthFiles',
  'depthTrackTables',
  'depthTrackFiles',
  'depthTrackLabels',
  'depthTrackColors',
  'depthTrackHeights',
  'depthTrackLargeTickIntervals',
  'depthTrackSmallTickIntervals',
  'depthTrackTickFontSizes'
]);

/** @param {unknown} value */
const hasOptionValue = (value) => value !== null && value !== undefined;

/**
 * @param {unknown} ref
 * @param {StoredResources} resources
 * @param {string} fieldName
 * @param {ResourceFileResolver | null} [resolveResourceFile]
 */
const canonicalDepthResourceFile = (
  ref,
  resources,
  fieldName,
  resolveResourceFile = null
) => {
  if (!isPlainObject(ref)) {
    throw new Error(`${fieldName} must be a canonical resource reference.`);
  }
  if (ref.representation !== 'file' && ref.representation !== 'canonicalTsv') {
    throw new Error(`${fieldName} has unsupported representation '${ref.representation}'.`);
  }
  if (typeof ref.resourceId !== 'string' || !ref.resourceId.trim()) {
    throw new Error(`${fieldName}.resourceId must be a non-empty string.`);
  }
  const resourceId = ref.resourceId.trim();
  return resolveResourceFile
    ? resolveResourceFile(resourceId)
    : resourceAsLegacyFile(resources, resourceId);
};

/**
 * @param {unknown} value
 * @param {string} fieldName
 */
const canonicalDepthNumber = (value, fieldName) => {
  if (value === null) return null;
  if (typeof value !== 'number' || !Number.isFinite(value) || value <= 0) {
    throw new Error(`${fieldName} must be null or a positive finite number.`);
  }
  return value;
};

/**
 * @param {unknown} value
 * @param {string} fieldName
 */
const canonicalDepthText = (value, fieldName) => {
  if (value === null) return null;
  if (typeof value !== 'string' || !value.trim()) {
    throw new Error(`${fieldName} must be null or a non-empty string.`);
  }
  return value.trim();
};

/**
 * @param {{
 *   options: Record<string, any>,
 *   records: unknown[],
 *   resources: Record<string, any> | null | undefined,
 *   mode: string,
 *   resolveResourceFile?: ResourceFileResolver | null
 * }} params
 */
const projectCanonicalDepthTracks = ({
  options,
  records,
  resources,
  mode,
  resolveResourceFile = null
}) => {
  const canonicalPresent = hasOptionValue(options.depthTracks);
  const legacyFields = LEGACY_DEPTH_OPTION_FIELDS.filter(
    (fieldName) => hasOptionValue(options[fieldName])
  );
  if (canonicalPresent && legacyFields.length > 0) {
    throw new Error(
      `diagramOptions.depthTracks cannot be combined with legacy depth fields: ${legacyFields.join(', ')}.`
    );
  }
  if (!canonicalPresent) return null;
  if (!Array.isArray(options.depthTracks) || options.depthTracks.length === 0) {
    throw new Error('diagramOptions.depthTracks must be a non-empty array.');
  }

  /** @type {unknown[][]} */
  const sourceRows = Array.from({ length: records.length }, () => []);
  /** @type {unknown[][]} */
  const fileRows = Array.from({ length: records.length }, () => []);
  const tracks = options.depthTracks.map((track, trackIndex) => {
    const fieldName = `diagramOptions.depthTracks[${trackIndex}]`;
    if (!track || typeof track !== 'object' || Array.isArray(track)) {
      throw new Error(`${fieldName} must be an object.`);
    }
    /** @type {unknown[]} */
    let sourceRefs;
    if (Array.isArray(track.source)) {
      if (track.source.length !== records.length) {
        throw new Error(
          `${fieldName}.source must contain one source per displayed record (${records.length}).`
        );
      }
      sourceRefs = track.source;
    } else {
      sourceRefs = Array.from({ length: records.length }, () => track.source);
    }
    if (!sourceRefs.some((ref) => ref !== null && ref !== undefined)) {
      throw depthSeriesWithoutSource(trackIndex);
    }
    sourceRefs.forEach((ref, recordIndex) => {
      sourceRows[recordIndex][trackIndex] = ref ?? null;
      fileRows[recordIndex][trackIndex] = ref === null || ref === undefined
        ? null
        : canonicalDepthResourceFile(
          ref,
          resources,
          `${fieldName}.source${Array.isArray(track.source) ? `[${recordIndex}]` : ''}`,
          resolveResourceFile
        );
    });
    if (mode === 'circular' && track.height !== null) {
      throw new Error(`${fieldName}.height must be null for Circular requests.`);
    }
    const label = canonicalDepthText(track.label, `${fieldName}.label`);
    const color = canonicalDepthText(track.color, `${fieldName}.color`);
    return {
      label: label ?? (options.depthTracks.length === 1 ? 'Depth' : `Depth ${trackIndex + 1}`),
      color: color ?? '#4A90E2',
      height: mode === 'linear'
        ? canonicalDepthNumber(track.height, `${fieldName}.height`)
        : null,
      large_tick_interval: canonicalDepthNumber(
        track.largeTickInterval,
        `${fieldName}.largeTickInterval`
      ),
      small_tick_interval: canonicalDepthNumber(
        track.smallTickInterval,
        `${fieldName}.smallTickInterval`
      ),
      tick_font_size: canonicalDepthNumber(
        track.tickFontSize,
        `${fieldName}.tickFontSize`
      )
    };
  });
  return { sourceRows, fileRows, tracks };
};

/**
 * @param {Readonly<Record<string, unknown>> | null | undefined} resources A Session's resources, or a committed
 *   Session's: each entry is checked here.
 * @param {string} resourceId
 */
export const decodeCanonicalResourceText = (resources, resourceId) => {
  const ownValue = ownEntry(resources, resourceId);
  if (!ownValue || typeof ownValue !== 'object' || Array.isArray(ownValue)) {
    throw new Error(`Missing canonical resource: ${resourceId}`);
  }
  // An object other than an array (checked above), read through a view.
  const entry = /** @type {UnvalidatedFields<'encoding' | 'data'> & object} */ (ownValue);
  if (entry.encoding && entry.encoding !== 'base64') {
    throw new Error(`Unsupported canonical resource encoding: ${entry.encoding}`);
  }
  if (typeof entry.data !== 'string') {
    throw new Error(`Canonical resource ${resourceId} has no text payload.`);
  }
  let bytes;
  try {
    bytes = base64ToBytes(entry.data);
  } catch (error) {
    throw new Error(`Canonical resource ${resourceId} contains invalid base64 data.`, { cause: error });
  }
  return bytesToText(bytes, { fatal: true });
};

/** @type {WeakMap<WeakKey, { genbank: number, fasta: number }>} */
const sourceRecordCounts = new WeakMap();

/**
 * @param {Readonly<Record<string, unknown>> | null | undefined} resources A committed Session's resources.
 * @param {string} resourceId
 * @param {'genbank' | 'fasta'} kind
 */
export const readCanonicalResourceRecordCount = async (resources, resourceId, kind) => {
  const descriptor = resources?.[resourceId];
  const owner = getResourcePayloadOwner(descriptor);
  // A set owner other than the descriptor is a value `setResourcePayloadOwner` stored: a WeakKey.
  const fileOwner = owner && owner !== descriptor ? /** @type {WeakKey} */ (owner) : null;
  let counts = fileOwner ? sourceRecordCounts.get(fileOwner) : undefined;
  if (!counts) {
    const text = fileOwner
      ? await readFileText(fileOwner)
      : decodeCanonicalResourceText(resources, resourceId);
    counts = {
      genbank: countGenBankRecords(text),
      fasta: (text.match(/^>/gm) || []).length
    };
    // Immutable File/view identity outlives transferred bytes. Retain only counts;
    // a replacement source gets its own entry, and ownerless descriptors stay uncached.
    if (fileOwner) sourceRecordCounts.set(fileOwner, counts);
  }
  return counts[kind];
};

/**
 * @param {Readonly<Record<string, unknown>> | null | undefined} resources
 * @param {{ resourceId?: unknown } | null | undefined} ref
 */
const resourceTextFromRef = (resources, ref) => (
  // A property key is the value's string form, so `String` reads the same entry.
  ref?.resourceId ? decodeCanonicalResourceText(resources, String(ref.resourceId)) : null
);

// A table the Session stores is named in its failure: the user corrects the
// Session file, not a table in this page.
/**
 * @template T
 * @param {string} table
 * @param {() => T} read
 * @returns {T}
 */
const readSessionTable = (table, read) => {
  try {
    return read();
  } catch (error) {
    const model = normalizeCaughtError(error);
    if (['UNKNOWN', 'VALIDATION_UNCLASSIFIED'].includes(model.code)) throw error;
    throw diagnosticError(model.code, { ...model.context, sessionTable: table }, { stage: model.stage });
  }
};

// Whether the committed Session drew with the given Feature visibility rules
// (TSV), compared after normalization. Reusing resolved protein comparisons is
// valid only when it holds (CO-02); an unreadable committed table declines.
// Feature visibility decides which proteins LOSATP compares: the rules and the
// per-feature `featureVisibility` of the request's identity rows (design Q4 3.3).
/** @param {unknown} rows */
const featureVisibilityRowsKey = (rows) => JSON.stringify((Array.isArray(rows) ? rows : [])
  .filter((row) => row?.featureVisibility)
  .map((row) => [row.recordKey, row.biologicalFeatureId, row.featureVisibility]));
/**
 * @param {CommittedCanonicalSession | null | undefined} committedSession
 * @param {unknown} activeRulesTsv
 * @param {unknown} [activeOverrideRows]
 */
export const committedFeatureVisibilityMatches = (committedSession, activeRulesTsv, activeOverrideRows = []) => {
  /** @param {unknown} text */
  const normalize = (text) => serializeFeatureVisibilityRules(parseFeatureVisibilityRules(String(text ?? '')).rules);
  try {
    return normalize(resourceTextFromRef(
      committedSession?.resources,
      editorTableRef(committedSession?.renderRequest?.diagramOptions, 'visibility')
    )) === normalize(activeRulesTsv)
      && featureVisibilityRowsKey(committedSession?.renderRequest?.diagramOptions?.featureOverrides)
        === featureVisibilityRowsKey(activeOverrideRows);
  } catch {
    return false;
  }
};

/**
 * @param {unknown} config
 * @param {string} path
 */
const nestedConfigValue = (config, path) => {
  let current = config;
  for (const key of path.split('.')) {
    if (!current || typeof current !== 'object' || Array.isArray(current)) return undefined;
    // The cast: an object other than an array; each property is an unvalidated value.
    current = /** @type {Readonly<Record<string, unknown>>} */ (current)[key];
  }
  return current;
};

/**
 * @param {Record<string, unknown>} overrides
 * @param {string} path
 */
const sharedLengthOverrideValue = (overrides, path) => {
  // The Web UI has one control for values that Python stores per genome length.
  // Project them only when both variants encode the same explicit setting.
  const shortValue = overrides[`${path}.short`];
  const longValue = overrides[`${path}.long`];
  if (shortValue === undefined || longValue === undefined) return undefined;
  if (Object.is(shortValue, longValue)) return longValue;
  const shortNumber = shortValue === null || String(shortValue).trim() === ''
    ? Number.NaN
    : Number(shortValue);
  const longNumber = longValue === null || String(longValue).trim() === ''
    ? Number.NaN
    : Number(longValue);
  return Number.isFinite(shortNumber) && shortNumber === longNumber ? longValue : undefined;
};

/** @param {{ mode: ModeProfileName, showLabels?: unknown, allowInnerLabels?: unknown }} input */
const legacyLabelScope = ({ mode, showLabels, allowInnerLabels }) => {
  /** @type {string | undefined} */
  let scope;
  if (showLabels !== undefined) {
    if (typeof showLabels === 'boolean') {
      scope = showLabels ? (mode === 'circular' ? 'outer' : 'all') : 'none';
    } else {
      const normalized = String(showLabels).trim().toLowerCase();
      /** @type {Readonly<Record<string, string | undefined>>} */
      const aliases = {
        true: 'all', yes: 'all', on: 'all',
        false: 'none', no: 'none', off: 'none'
      };
      const policy = aliases[normalized] || normalized;
      if (!['all', 'first', 'orthogroup_top', 'none'].includes(policy)) {
        throw new Error(`Unsupported persisted label policy: ${showLabels}`);
      }
      if (mode === 'circular' && ['first', 'orthogroup_top'].includes(policy)) {
        throw new Error(`Circular labels cannot use Linear-only policy '${policy}'.`);
      }
      scope = mode === 'circular' ? (policy === 'all' ? 'outer' : 'none') : policy;
    }
  }
  if (mode === 'circular' && allowInnerLabels === true) {
    if (scope === undefined || scope === 'outer') scope = 'both';
  }
  return scope;
};

/** @param {unknown} config @param {ModeProfileName} mode */
const projectFullConfigOverrides = (config, mode) => {
  if (!isPlainObject(config)) return {};
  /** @type {Set<string>} */
  const paths = new Set(Object.values(CONFIG_OVERRIDE_PATHS));
  const labelScopePath = MODE_LABEL_SCOPE_PATHS[mode];
  paths.add(labelScopePath);
  Object.values(SHARED_LENGTH_CONFIG_OVERRIDE_PATHS).forEach((path) => {
    paths.add(`${path}.short`);
    paths.add(`${path}.long`);
  });
  const labelFontPath = mode === 'linear' ? 'labels.font_size.linear' : 'labels.font_size';
  paths.add(`${labelFontPath}.short`);
  paths.add(`${labelFontPath}.long`);
  Object.values(LINEAR_DEFINITION_STYLE_PATHS).forEach((prefix) => {
    LINEAR_DEFINITION_STYLE_FIELDS.forEach((field) => paths.add(`${prefix}.${field}`));
  });
  /** @type {Record<string, unknown>} */
  const projected = {};
  paths.forEach((path) => {
    const value = nestedConfigValue(config, path);
    if (value !== undefined) projected[path] = value;
  });
  // A number Python's config reads as unset ("none", "auto") is the Web's null (OV-791).
  CONFIG_UNSET_NUMBER_PATHS.forEach((path) => {
    const value = projected[path];
    if (typeof value === 'string' && CONFIG_UNSET_NUMBER_SPELLINGS.includes(value.trim().toLowerCase())) {
      projected[path] = null;
    }
  });
  if (!Object.prototype.hasOwnProperty.call(projected, labelScopePath)) {
    const modeShowLabels = nestedConfigValue(config, `canvas.${mode}.show_labels`);
    const sharedShowLabels = nestedConfigValue(config, 'canvas.show_labels');
    const scope = legacyLabelScope({
      mode,
      showLabels: modeShowLabels === undefined ? sharedShowLabels : modeShowLabels,
      allowInnerLabels: mode === 'circular'
        ? nestedConfigValue(config, 'canvas.circular.allow_inner_labels')
        : undefined
    });
    if (scope !== undefined) projected[labelScopePath] = scope;
  }
  return projected;
};

/**
 * @param {Record<string, unknown>} overrides
 * @param {ModeProfileName} mode
 */
const projectCanonicalConfigOverrides = (overrides, mode) => {
  /** @type {Record<string, any>} */
  const projected = {};
  for (const [semanticName, path] of Object.entries(CONFIG_OVERRIDE_PATHS)) {
    if (Object.prototype.hasOwnProperty.call(overrides, path)) {
      projected[legacyFlatConfigKey(semanticName)] = overrides[path];
    }
  }
  const labelScopePath = MODE_LABEL_SCOPE_PATHS[mode];
  if (Object.prototype.hasOwnProperty.call(overrides, labelScopePath)) {
    projected.label_scope = overrides[labelScopePath];
  }
  for (const [semanticName, path] of Object.entries(SHARED_LENGTH_CONFIG_OVERRIDE_PATHS)) {
    const value = sharedLengthOverrideValue(overrides, path);
    if (value !== undefined) projected[legacyFlatConfigKey(semanticName)] = value;
  }
  const labelFontPath = mode === 'linear' ? 'labels.font_size.linear' : 'labels.font_size';
  const labelFontSize = sharedLengthOverrideValue(overrides, labelFontPath);
  if (labelFontSize !== undefined) projected.label_font_size = labelFontSize;

  /** @type {Record<string, Record<string, unknown>>} */
  const lineStyles = {};
  for (const [kind, prefix] of Object.entries(LINEAR_DEFINITION_STYLE_PATHS)) {
    /** @type {Record<string, unknown>} */
    const style = {};
    for (const field of LINEAR_DEFINITION_STYLE_FIELDS) {
      const path = `${prefix}.${field}`;
      if (Object.prototype.hasOwnProperty.call(overrides, path)) {
        style[field] = overrides[path];
      }
    }
    if (Object.keys(style).length > 0) lineStyles[kind] = style;
  }
  if (Object.keys(lineStyles).length > 0) {
    projected.linear_definition_line_styles = lineStyles;
  }
  return projected;
};

/** @param {Record<string, unknown>} overrides */
const projectLegacyFlatConfigOverrides = (overrides) => Object.fromEntries(
  Object.entries(overrides).filter(([key]) => (
    !key.includes('.') && !['show_labels', 'allow_inner_labels'].includes(key)
  ))
);

/**
 * @param {Record<string, unknown>} overrides
 * @param {ModeProfileName} mode
 */
const projectExplicitConfigOverrides = (overrides, mode) => {
  const legacyLabelPaths = new Set([
    'canvas.show_labels',
    'canvas.circular.show_labels',
    'canvas.linear.show_labels',
    'canvas.circular.allow_inner_labels'
  ]);
  const projected = Object.fromEntries(
    Object.entries(overrides).filter(([key]) => (
      key.includes('.') && !legacyLabelPaths.has(key)
    ))
  );
  const labelScopePath = MODE_LABEL_SCOPE_PATHS[mode];
  if (!Object.prototype.hasOwnProperty.call(projected, labelScopePath)) {
    const modeShowLabels = overrides[`canvas.${mode}.show_labels`];
    const sharedShowLabels = overrides['canvas.show_labels'];
    const scope = legacyLabelScope({
      mode,
      showLabels: modeShowLabels ?? sharedShowLabels ?? overrides.show_labels,
      allowInnerLabels: mode === 'circular'
        ? (
            overrides['canvas.circular.allow_inner_labels']
            ?? overrides.allow_inner_labels
          )
        : undefined
    });
    if (scope !== undefined) projected[labelScopePath] = scope;
  }
  return projected;
};

/**
 * @param {NonNullable<UnvalidatedFields<'conservationLabels' | 'conservationColors' | 'conservationReference'
 *   | 'conservationRingWidth' | 'conservationRingGap'>>} options The request's `diagramOptions` (unvalidated), or `{}`.
 * @param {{ c_conservation_blasts?: FileLike[] }} files
 */
const projectCircularConservationConfig = (options, files) => {
  const sourceFiles = Array.isArray(files.c_conservation_blasts)
    ? files.c_conservation_blasts
    : [];
  if (sourceFiles.length === 0) return undefined;
  const labels = Array.isArray(options.conservationLabels)
    ? options.conservationLabels.map((/** @type {unknown} */ value) => String(value || '').trim())
    : [];
  const colors = Array.isArray(options.conservationColors)
    ? options.conservationColors.map((/** @type {unknown} */ value) => String(value || '').trim())
    : [];
  const series = sourceFiles.map((file, index) => {
    const fileName = String(file?.name || `comparison-${index + 1}.tsv`);
    // A Session without stored labels predates D-03 and drew the full file name,
    // as the Python decoder reads it; every Web Session stores its labels.
    const defaultLabel = Array.isArray(options.conservationLabels)
      ? fileName.replace(/\.[^.]+$/, '').trim() || `Comparison ${index + 1}`
      : fileName;
    return {
      fileName,
      sourceIndex: index,
      label: labels[index] || defaultLabel,
      color: colors[index] || '',
      losat_gencode: 1
    };
  });
  return {
    enabled: true,
    source: 'upload',
    losat_program: 'blastn',
    subject_gencode: 1,
    reference: String(options.conservationReference || 'auto'),
    labels: series.map((entry) => entry.label).join(','),
    series,
    ring_width: projectOptionalNumber(options.conservationRingWidth, { field: 'conservation_ring_width' }),
    ring_gap: projectOptionalNumber(options.conservationRingGap, { field: 'conservation_ring_gap' })
  };
};

/** @param {unknown} slot */
const projectCurrentCanonicalCircularSlot = (slot) => {
  if (
    !isPlainObject(slot) ||
    ['spacing', 'inner_gap_px', 'outer_gap_px', 'strict'].some((field) => (
      Object.prototype.hasOwnProperty.call(slot, field)
    ))
  ) throw new Error('Current canonical circular track slot uses an obsolete shape.');
  for (const field of ['innerGapPx', 'outerGapPx']) {
    const value = slot[field];
    if (value !== null && value !== undefined && (typeof value !== 'number' || !Number.isFinite(value) || value < 0)) {
      throw new Error(`Canonical Circular track ${field} must be a nonnegative finite number or null.`);
    }
  }
  return {
    id: String(slot.id || ''),
    renderer: circularRendererName(slot.renderer),
    enabled: slot.enabled !== false,
    width: projectCircularMeasureDraft(slot.width) ?? null, radius: projectCircularMeasureDraft(slot.radius) ?? null,
    inner_gap_px: slot.innerGapPx == null ? null : String(slot.innerGapPx),
    outer_gap_px: slot.outerGapPx == null ? null : String(slot.outerGapPx),
    side: slot.side ?? null, z: Number(slot.z) || 0,
    params: cloneCanonicalJsonValue(slot.params || {})
  };
};

/**
 * @param {Record<string, any> | null | undefined} resources
 * @param {Array<Record<string, any> | null | undefined>} records
 * @param {string} [originalName]
 * @param {{
 *   resolveResourceFile?: ResourceFileResolver | null,
 *   sessionResourceTable?: SessionResourceTable | null
 * }} [options]
 */
const combineCircularGenbankResources = (
  resources,
  records,
  originalName = '',
  { resolveResourceFile = null, sessionResourceTable = null } = {}
) => {
  /** @type {string[]} */
  const resourceIds = [];
  const seen = new Set();
  records.forEach((record) => {
    const source = record?.source || {};
    const resourceId = source.kind === 'genbank' ? String(source.resourceId || '').trim() : '';
    if (!resourceId || seen.has(resourceId)) return;
    seen.add(resourceId);
    resourceIds.push(resourceId);
  });
  if (resourceIds.length === 0) return null;

  if (!sessionResourceTable && resourceIds.length === 1) {
    return resolveResourceFile
      ? resolveResourceFile(resourceIds[0])
      : resourceAsLegacyFile(resources, resourceIds[0]);
  }
  return createCombinedSessionResourceFileView(
    sessionResourceTable || adoptCurrentSessionResources(resources),
    resourceIds.map(resourceId => ({ resourceId })),
    {
      name: normalizeOriginalResourceName(originalName) || 'canonical-circular-records.gb',
      type: 'text/plain'
    }
  );
};

/** @param {Readonly<Record<string, unknown>>} renderRequest @param {readonly unknown[]} records */
const canonicalGrouping = (renderRequest, records) => {
  const allowed = renderRequest.mode === 'circular'
    ? ['single', 'grid', 'batch']
    : ['single'];
  if (typeof renderRequest.grouping !== 'string' || !allowed.includes(renderRequest.grouping)) {
    throw new Error(`Unsupported canonical ${renderRequest.mode} grouping.`);
  }
  if (
    renderRequest.mode === 'circular' &&
    renderRequest.grouping === 'single' &&
    (
      records.length !== 1 ||
      Object.keys(renderRequest.layout || {}).length > 0
    )
  ) {
    throw new Error('A single Circular canonical request requires one record and no grid layout.');
  }
  if (
    renderRequest.grouping === 'batch' &&
    Object.keys(renderRequest.layout || {}).length > 0
  ) {
    throw new Error('A Circular batch canonical request cannot define a grid layout.');
  }
  return renderRequest.grouping;
};

/**
 * @param {{ output?: unknown }} renderRequest
 * @param {string} grouping
 * @param {number} recordCount
 */
const canonicalOutputPrefixes = (renderRequest, grouping, recordCount) => {
  const outputs = grouping === 'batch' ? renderRequest.output : [renderRequest.output];
  if (
    !Array.isArray(outputs) ||
    outputs.length !== (grouping === 'batch' ? recordCount : 1) ||
    outputs.some((/** @type {unknown} */ output) => !output || typeof output !== 'object' || Array.isArray(output))
  ) {
    throw new Error(
      grouping === 'batch'
        ? 'Canonical Circular batch output must contain one object per record.'
        : 'Canonical renderRequest output must be an object.'
    );
  }
  // Each output passed the object check above.
  const prefixes = /** @type {{ prefix?: unknown }[]} */ (outputs).map((output) => String(output.prefix || '').trim());
  if (prefixes.some((prefix) => !prefix)) {
    throw new Error('Canonical renderRequest output prefixes must be non-empty.');
  }
  return prefixes;
};

/** @param {({ selector?: { kind?: unknown, value?: unknown } | null } | null | undefined)[]} records */
const inferredImplicitBatchPrefixes = (records) => {
  const recordIds = records.map((record) => (
    record?.selector?.kind === 'recordId'
      ? String(record.selector.value || '').trim()
      : ''
  ));
  if (recordIds.some((recordId) => !recordId)) return null;
  return resolveCircularBatchPrefixes(
    recordIds.map((recordId) => ({ record_id: recordId })),
    null
  );
};

/**
 * @param {{ mode?: unknown }} renderRequest
 * @param {string} grouping
 * @param {({ selector?: { kind?: unknown, value?: unknown } | null } | null | undefined)[]} records
 * @param {string[]} prefixes
 * @param {boolean | null | undefined} circularPrefixExplicit
 */
const projectedOutputPrefix = (
  renderRequest,
  grouping,
  records,
  prefixes,
  circularPrefixExplicit
) => {
  if (renderRequest.mode === 'circular' && circularPrefixExplicit === false) {
    return '';
  }
  if (grouping !== 'batch') return prefixes[0] || 'out';
  const implicitPrefixes = inferredImplicitBatchPrefixes(records);
  if (
    circularPrefixExplicit !== true &&
    implicitPrefixes &&
    implicitPrefixes.every((prefix, index) => prefix === prefixes[index])
  ) {
    return '';
  }
  if (prefixes.length === 1) return prefixes[0];
  const firstMatch = prefixes[0]?.match(/^(.*)_1$/);
  const base = firstMatch?.[1] || '';
  return base && prefixes.every((prefix, index) => prefix === `${base}_${index + 1}`)
    ? base
    : '';
};

/**
 * A reader of unvalidated Session data (R14): the input is not typed.
 * @param {Record<string, any>} input
 */
export const projectCanonicalSessionRequest = ({
  renderRequest,
  resources: canonicalResources,
  webFiles = {},
  legacyFiles = null,
  storedConfig = null,
  initializeCliInputs = false,
  fileBindings = [],
  linearTrackSlotSchemaVersion = LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  repairLegacyTableRows = false,
  sessionResourceTable = null,
  deferResourceContent = false,
  adoptCanonicalPayloads = false
}) => {
  // Every request reaches the projection at the current schema: a stored one
  // through its reader (`readStoredRenderRequest`), a new one from the writer.
  if (!renderRequest || renderRequest.schema !== CANONICAL_REQUEST_SCHEMA) {
    throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'SESSION_FORMAT' });
  }
  if (!['circular', 'linear'].includes(renderRequest.mode)) {
    throw new Error('Unsupported canonical renderRequest mode.');
  }
  /** @type {ProjectionRecordView[]} */
  const sourceRecords = Array.isArray(renderRequest.records) ? renderRequest.records : [];
  sourceRecords.forEach((record) => canonicalRecordDisplay(record.display));
  if (!Array.isArray(renderRequest.diagramOptions?.featurePlacements)) {
    throw new Error('Canonical featurePlacements must be an array.');
  }
  canonicalFeaturePlacements(renderRequest.diagramOptions.featurePlacements, renderRequest.mode);
  if (!Array.isArray(renderRequest.diagramOptions?.featureOverrides)) {
    throw diagnosticError('INPUT_INVALID', { field: 'featureOverrides', reason: 'SESSION_FORMAT' });
  }
  const projectedFeatureOverrides = canonicalFeatureOverrides(renderRequest.diagramOptions.featureOverrides);
  const normalizedRecordOrdering = normalizeWebGridColumnOrdering(sourceRecords);
  const records = normalizedRecordOrdering.records;
  /**
   * @template T
   * @param {readonly T[]} values Indexed by the request's record order.
   * @returns {T[]}
   */
  const reorderRecordIndexedValues = (values) => (
    normalizedRecordOrdering.sourceIndexByProjectedIndex
      .map((sourceIndex) => values?.[sourceIndex])
  );
  if (records.length === 0) throw new Error('Canonical renderRequest records are required.');
  /** @type {ReturnType<typeof canonicalSimilarityAlignment> | null} */
  let similarityAlignment = null;
  /** @type {ReturnType<typeof canonicalRecordTranslations>} */
  let recordTranslations = [];
  if (renderRequest.mode === 'linear') {
    const layout = renderRequest.layout || {};
    if (Object.keys(layout).length > 0) {
      if (!Number.isFinite(layout.recordGapPx) || layout.recordGapPx < 0) {
        throw new Error('renderRequest.layout.recordGapPx must be a finite non-negative number.');
      }
      if (layout.multiRecordPositions !== null && (
        !Array.isArray(layout.multiRecordPositions) ||
        layout.multiRecordPositions.some((/** @type {unknown} */ token) => typeof token !== 'string' || !token.trim())
      )) throw new Error(
        'renderRequest.layout.multiRecordPositions must be null or an array of non-empty text.'
      );
      if (layout.multiRecordPositions?.length && sourceRecords.some(
        (record) => record.presentation?.gridRow != null
      )) throw new Error(
        'Linear row placement must use record presentation or layout positions, not both.'
      );
      ({ similarityAlignment, recordTranslations } = canonicalLinearRecordLayout(renderRequest));
    }
  }
  const grouping = canonicalGrouping(renderRequest, records);
  const sourceOutputPrefixes = canonicalOutputPrefixes(
    renderRequest,
    grouping,
    records.length
  );
  const outputPrefixes = grouping === 'batch'
    ? reorderRecordIndexedValues(sourceOutputPrefixes)
    : sourceOutputPrefixes;
  let webMetadata = webFiles && typeof webFiles === 'object' && !Array.isArray(webFiles)
    ? webFiles
    : {};
  // The checked bindings; the CLI initialization below rebuilds them field by field.
  /** @type {Readonly<Record<string, unknown>> | undefined} */
  let explicitBindings = validateWebFileBindings(webMetadata, canonicalResources);
  // A CLI sidecar has no saved Web draft. Its empty writer slots are initial
  // values, not a user's cleared inputs. Keep real original CLI bindings (in
  // particular GFF + FASTA), and initialize absent slots from the typed request.
  // A saved Web draft, including explicit null/[] inputs, always wins unchanged.
  if (initializeCliInputs && storedConfig == null && explicitBindings) {
    const originalInputFields = renderRequest.mode === 'circular'
      ? ['c_gb', 'c_gff', 'c_fasta']
      : ['linearSeqs'];
    // A const, so that the narrowing of `explicitBindings` reaches the closure.
    const bindings = explicitBindings;
    const hasOriginalInputs = originalInputFields.some((field) => {
      const value = bindings[field];
      return Array.isArray(value) ? value.length > 0 : value != null;
    });
    explicitBindings = Object.fromEntries(Object.entries(explicitBindings).filter(
      ([field, value]) => (hasOriginalInputs && originalInputFields.includes(field))
        || (value != null && (!Array.isArray(value) || value.length > 0))
    ));
    // A CLI binding uid (`cli-seq-N`) is only an initial value. The committed
    // request owns record identity, so each Linear file takes the recordKey
    // that Inherit matches against. A multi-record file, which Python expands
    // as `<fileKey>:<n>`, becomes one row per record. The request also owns
    // each record's selection, crop, and orientation: the CLI binds the input
    // file and writes them in the request, so a file of which the CLI drew
    // some records stays one File that draws only those records. A record
    // without a selector (a CLI Session that stored a copy of each drawn
    // record) takes the `#n` of its expanded recordKey.
    /** @type {Map<string, string[]>} */
    const recordKeysByFile = new Map();
    /** @type {Map<string, ProjectionRecordView>} */
    const recordsByKey = new Map();
    sourceRecords.forEach((record) => {
      const key = String(record.recordKey || '');
      const fileKey = key.replace(/:[1-9]\d*$/, '');
      recordKeysByFile.set(fileKey, [...(recordKeysByFile.get(fileKey) || []), key]);
      recordsByKey.set(key, record);
    });
    const fileRecordKeys = [...recordKeysByFile.keys()];
    if (renderRequest.mode === 'linear' && Array.isArray(explicitBindings.linearSeqs)
      && fileRecordKeys.every(Boolean)
      && fileRecordKeys.length === explicitBindings.linearSeqs.length) {
      explicitBindings = {
        ...explicitBindings,
        linearSeqs: explicitBindings.linearSeqs.flatMap((sequence, index) => (
          // A key of the map: `fileRecordKeys` lists its keys.
          /** @type {string[]} */ (recordKeysByFile.get(fileRecordKeys[index])).map((uid) => {
            const fields = linearSeqRecordFields(recordsByKey.get(uid));
            return {
              ...sequence,
              uid,
              ...fields,
              ...(fields.region_record_id || uid === fileRecordKeys[index] ? {} : {
                region_record_id: `#${uid.slice(fileRecordKeys[index].length + 1)}`
              })
            };
          })
        ))
      };
    }
    webMetadata = { ...webMetadata, bindings: explicitBindings };
  }
  const storedResourceOriginalNames = webMetadata.resourceOriginalNames;
  const originalNameHints = {
    ...legacyResourceOriginalNames({ renderRequest, legacyFiles, fileBindings }),
    ...(storedResourceOriginalNames && typeof storedResourceOriginalNames === 'object' &&
      !Array.isArray(storedResourceOriginalNames) ? storedResourceOriginalNames : {})
  };
  const resources = sessionResourceTable
    ? canonicalResources
    : resourcesWithOriginalNames(canonicalResources, originalNameHints);
  /** @type {ResourceFileResolver | null} */
  const resolveResourceFile = sessionResourceTable
    ? (resourceId, metadata = {}) => {
        const descriptor = adoptedSessionResourceDescriptor(
          sessionResourceTable,
          resourceId
        );
        const storedName = normalizeOriginalResourceName(descriptor.name);
        const prefix = `${resourceId}-`;
        let inferredName = storedName;
        while (inferredName.startsWith(prefix) && inferredName.length > prefix.length) {
          inferredName = inferredName.slice(prefix.length);
        }
        const displayName = normalizeOriginalResourceName(metadata.name)
          // `String` is the property key the read used before: an id is text or a JSON value.
          || normalizeOriginalResourceName(originalNameHints[String(resourceId)])
          || inferredName;
        return createSessionResourceFileView(
          sessionResourceTable,
          resourceId,
          {
            ...metadata,
            ...(displayName ? { name: displayName } : {})
          }
        );
      }
    : null;
  const legacyCircularInputBinding = (Array.isArray(fileBindings) ? fileBindings : [])
    .find((binding) => /^(?:files\.)?c_gb$/.test(String(binding?.slot || '')));
  const circularInputOriginalName = normalizeOriginalResourceName(
    webMetadata.circularInputOriginalName ||
    legacyFiles?.c_gb?.name ||
    legacyCircularInputBinding?.name ||
    (records.length === 1 ? resources?.[records[0]?.source?.resourceId]?.name : '')
  );
  /** @type {readonly unknown[]} */
  const savedLinearRecordMetadata = Array.isArray(webMetadata.linearRecordMetadata)
    ? webMetadata.linearRecordMetadata
    : [];
  const savedLinearRecordMetadataByKey = new Map(
    savedLinearRecordMetadata
      .map((entry) => /** @type {const} */ ([
        String(/** @type {UnvalidatedFields<'recordKey'>} */ (entry)?.recordKey || ''), entry
      ]))
      .filter(([recordKey]) => recordKey)
  );
  const legacyLinearSequences = Array.isArray(legacyFiles?.linearSeqs)
    ? legacyFiles.linearSeqs
    : [];
  const projectedProteinPipeline = projectGeneratedProteinPipeline(
    (renderRequest.comparisons || []).find(
      (/** @type {UnvalidatedFields<'kind'>} */ comparison) => comparison?.kind === 'generatedProteinComparison'
    ),
    { adoptCanonicalPayloads }
  );
  const comparisonsContainGeneratedProteinPipeline = (
    renderRequest.comparisons || []
  ).some((/** @type {UnvalidatedFields<'kind'>} */ comparison) => comparison?.kind === 'generatedProteinComparison');
  // A CLI sidecar has no Web comparison draft, so its plan is the Web default
  // (No comparison) unless the projection states one. Its editable protein
  // pipeline is the adjacent LOSATP comparison the CLI drew. Without -b or a
  // protein mode the CLI still writes its protein settings as a disabled
  // pipeline (mode `none`, no pairs): that draft selects no LOSATP program. A
  // read-only comparison (-b) reaches this projection as an empty list.
  const committedComparisons = renderRequest.comparisons || [];
  const cliLinearSidecar = initializeCliInputs && storedConfig == null
    && renderRequest.mode === 'linear';
  const cliDraftWithoutComparison = cliLinearSidecar && drawsNoComparison(committedComparisons);
  const cliProteinPlan = cliLinearSidecar && projectedProteinPipeline && !cliDraftWithoutComparison
    ? { linearComparisonPlan: normalizeLinearComparisonPlan({ mode: 'adjacent' }) }
    : {};
  const {
    blastSource: _cliBlastSource, losatProgram: _cliLosatProgram, ...cliProteinSettings
  } = projectedProteinPipeline?.config || {};
  // The legacy file bag that the Web form state reads; its fields depend on the mode.
  /** @type {Record<string, any>} */
  const files = { linearSeqs: [] };
  if (renderRequest.mode === 'circular') {
    files.circularRecords = records.map((record) => {
      const source = record.source || {};
      const resourceFile = (/** @type {string} */ id) => resolveResourceFile
        ? resolveResourceFile(id)
        : resourceAsLegacyFile(resources, id);
      return {
        recordKey: String(record.recordKey || ''),
        cardinality: record.cardinality || 'exactly_one',
        sourceKind: source.kind,
        gb: source.kind === 'genbank' ? resourceFile(source.resourceId) : null,
        gff: source.kind === 'gffFasta' ? resourceFile(source.gffResourceId) : null,
        fasta: source.kind === 'gffFasta' ? resourceFile(source.fastaResourceId) : null,
        selector: cloneCanonicalJsonValue(record.selector),
        region: cloneCanonicalJsonValue(record.region),
        presentation: cloneCanonicalJsonValue(record.presentation),
        display: cloneCanonicalJsonValue(record.display || { isCircular: null, startCoordinate: null })
      };
    });
    const source = records[0]?.source || {};
    if (source.kind === 'genbank' && !Object.hasOwn(explicitBindings || {}, 'c_gb')) {
      files.c_gb = combineCircularGenbankResources(
        resources,
        records,
        circularInputOriginalName,
        { resolveResourceFile, sessionResourceTable }
      );
    }
    if (source.kind === 'gffFasta') {
      files.c_gff = resolveResourceFile
        ? resolveResourceFile(source.gffResourceId)
        : resourceAsLegacyFile(resources, source.gffResourceId);
      files.c_fasta = resolveResourceFile
        ? resolveResourceFile(source.fastaResourceId)
        : resourceAsLegacyFile(resources, source.fastaResourceId);
    }
  } else {
    files.linearSeqs = records.map((record, index) => {
      const source = record.source || {};
      const sourceIndex = normalizedRecordOrdering.sourceIndexByProjectedIndex[index];
      const savedMetadata = savedLinearRecordMetadataByKey.get(String(record.recordKey || '')) ||
        savedLinearRecordMetadata[sourceIndex] || legacyLinearSequences[sourceIndex] || {};
      const fileDefinition = String(savedMetadata.fileDefinition ?? savedMetadata.file_definition ?? '');
      const fileSubtitle = String(savedMetadata.fileSubtitle ?? savedMetadata.file_subtitle ?? '');
      const recordLabel = String(record.presentation?.label || '');
      const recordSubtitle = String(record.presentation?.subtitle || '');
      return {
        uid: String(record.recordKey || `canonical-seq-${index + 1}`),
        gb: source.kind === 'genbank'
          ? (resolveResourceFile
              ? resolveResourceFile(source.resourceId)
              : resourceAsLegacyFile(resources, source.resourceId))
          : null,
        gff: source.kind === 'gffFasta'
          ? (resolveResourceFile
              ? resolveResourceFile(source.gffResourceId)
              : resourceAsLegacyFile(resources, source.gffResourceId))
          : null,
        fasta: source.kind === 'gffFasta'
          ? (resolveResourceFile
              ? resolveResourceFile(source.fastaResourceId)
              : resourceAsLegacyFile(resources, source.fastaResourceId))
          : null,
        depth: null,
        blast: null,
        losat_gencode: integerSettingOr(savedMetadata.losatGencode ?? savedMetadata.losat_gencode, 1, 1),
        losat_filename: String(
          savedMetadata.losatFilename ?? savedMetadata.losat_filename ?? ''
        ),
        definition: resolveSavedRecordOverride({
          savedOverride: savedMetadata.recordDefinition,
          fileDefault: fileDefinition,
          resolved: recordLabel
        }),
        record_subtitle: resolveSavedRecordOverride({
          savedOverride: savedMetadata.recordSubtitle,
          fileDefault: fileSubtitle,
          resolved: recordSubtitle
        }),
        file_definition: fileDefinition,
        file_subtitle: fileSubtitle,
        ...linearSeqRecordFields(record)
      };
    });
    files.linearComparisons = [];
    files.linearCanonicalComparisons = [];
    (renderRequest.comparisons || [])
      .filter((/** @type {UnvalidatedFields<'kind'>} */ comparison) => comparison?.kind === 'nucleotideBlast')
      .forEach((
        // The filter kept only comparisons with a kind: none is null.
        /** @type {NonNullable<UnvalidatedFields<'queryRecordIndex' | 'subjectRecordIndex' | 'resourceId'>>} */ comparison,
        /** @type {number} */ index
      ) => {
        const queryIndex = normalizedRecordOrdering.projectedIndexBySourceIndex
          .get(Number(comparison.queryRecordIndex));
        const subjectIndex = normalizedRecordOrdering.projectedIndexBySourceIndex
          .get(Number(comparison.subjectRecordIndex));
        const file = resolveResourceFile
          ? resolveResourceFile(comparison.resourceId)
          : resourceAsLegacyFile(resources, comparison.resourceId);
        // An unmapped source index is undefined; it never named a sequence here.
        if (queryIndex === undefined || subjectIndex === undefined) return;
        if (files.linearSeqs[queryIndex] && subjectIndex === queryIndex + 1) {
          files.linearSeqs[queryIndex].blast = file;
        }
        if (!files.linearSeqs[queryIndex] || !files.linearSeqs[subjectIndex]) return;
        files.linearComparisons.push({
          id: `linear-comparison-canonical-${index + 1}`,
          queryUid: files.linearSeqs[queryIndex].uid,
          subjectUid: files.linearSeqs[subjectIndex].uid,
          source: 'upload',
          file
        });
      });
    (renderRequest.comparisons || []).forEach((/** @type {CanonicalComparison} */ comparison) => {
      if (isResourceBackedCanonicalComparison(comparison)) {
        const sourceQueryRecordIndex = Number(comparison.queryRecordIndex);
        const sourceSubjectRecordIndex = Number(comparison.subjectRecordIndex);
        const queryRecordIndex = normalizedRecordOrdering.projectedIndexBySourceIndex
          .get(sourceQueryRecordIndex);
        const subjectRecordIndex = normalizedRecordOrdering.projectedIndexBySourceIndex
          .get(sourceSubjectRecordIndex);
        if (
          comparison.kind === 'precomputedProteinComparison' &&
          (
            queryRecordIndex === undefined ||
            subjectRecordIndex === undefined ||
            !Number.isInteger(queryRecordIndex) ||
            !Number.isInteger(subjectRecordIndex) ||
            !files.linearSeqs[queryRecordIndex] ||
            !files.linearSeqs[subjectRecordIndex]
          )
        ) return;
        files.linearCanonicalComparisons.push(
          {
            ...mapResourceBackedCanonicalComparison(
              comparison,
              () => resolveResourceFile
                ? resolveResourceFile(comparison.resourceId)
                : resourceAsLegacyFile(resources, comparison.resourceId)
            ),
            ...(comparison.kind === 'precomputedProteinComparison'
              ? { queryRecordIndex, subjectRecordIndex }
              : {}),
            // This is in-memory projection provenance, not a canonical-schema
            // field. Direct CLI/Python comparison options have no Web pipeline
            // marker and must survive a projection/rebuild unchanged.
            ...(
              !comparisonsContainGeneratedProteinPipeline
                ? { canonicalInput: true }
                : {}
            )
          }
        );
        return;
      }
      if (comparison?.kind === 'generatedProteinComparison') {
        const projectedComparison = adoptCanonicalPayloads
          ? { ...comparison }
          : cloneCanonicalJsonValue(comparison);
        projectedComparison.pairs = (Array.isArray(comparison.pairs) ? comparison.pairs : [])
          .map((/** @type {UnvalidatedFields<'queryRecordIndex' | 'subjectRecordIndex'>} */ pair) => ({
            queryRecordIndex: normalizedRecordOrdering.projectedIndexBySourceIndex
              .get(Number(pair?.queryRecordIndex)),
            subjectRecordIndex: normalizedRecordOrdering.projectedIndexBySourceIndex
              .get(Number(pair?.subjectRecordIndex))
          }));
        files.linearCanonicalComparisons.push(
          projectedComparison
        );
        (Array.isArray(comparison.pairs) ? comparison.pairs : [])
          .forEach((
            /** @type {UnvalidatedFields<'queryRecordIndex' | 'subjectRecordIndex'>} */ pair,
            /** @type {number} */ index
          ) => {
            const queryIndex = normalizedRecordOrdering.projectedIndexBySourceIndex
              .get(Number(pair?.queryRecordIndex));
            const subjectIndex = normalizedRecordOrdering.projectedIndexBySourceIndex
              .get(Number(pair?.subjectRecordIndex));
            if (queryIndex === undefined || subjectIndex === undefined) return;
            if (!files.linearSeqs[queryIndex] || !files.linearSeqs[subjectIndex]) return;
            files.linearComparisons.push({
              id: `linear-comparison-canonical-losat-${index + 1}`,
              queryUid: files.linearSeqs[queryIndex].uid,
              subjectUid: files.linearSeqs[subjectIndex].uid,
              source: 'losat',
              file: null
            });
          });
      }
    });
  }

  const options = renderRequest.diagramOptions || {};
  const canonicalDepth = projectCanonicalDepthTracks({
    options,
    records,
    resources,
    mode: renderRequest.mode,
    resolveResourceFile
  });
  const projectedCircularSizeMode = renderRequest.mode === 'circular'
    ? requireCurrentCircularMultiRecordSizeMode(renderRequest.layout?.multiRecordSizeMode)
    : 'auto';
  const sourceDepthRows = canonicalDepth?.sourceRows || (
    Array.isArray(options.depthTrackFiles) ? options.depthTrackFiles : []
  );
  const depthRows = sourceDepthRows.length > 0
    ? reorderRecordIndexedValues(sourceDepthRows)
    : sourceDepthRows;
  const depthFileRows = canonicalDepth?.fileRows
    ? reorderRecordIndexedValues(canonicalDepth.fileRows)
    : null;
  if (depthRows.length > 0) {
    if (depthRows.length !== records.length || depthRows.some((/** @type {unknown} */ row) => !Array.isArray(row))) {
      throw new Error(
        `Canonical ${renderRequest.mode === 'circular' ? 'Circular' : 'Linear'} Depth matrix must contain one row per record (${records.length}).`
      );
    }
  }
  if (renderRequest.mode === 'circular' && depthRows.length > 0) {
    files.c_depth = normalizeRecordMajorDepthFileRows(
      // Each row is an array: the check above throws otherwise.
      depthFileRows || depthRows.map((/** @type {readonly UnvalidatedFields<'resourceId'>[]} */ row) => row.map((ref) => (
        ref?.resourceId
          ? (resolveResourceFile
              ? resolveResourceFile(ref.resourceId)
              : resourceAsLegacyFile(resources, ref.resourceId))
          : null
      ))),
      records.length
    );
  }
  if (renderRequest.mode === 'linear') {
    // Each row is an array: the check above throws otherwise.
    depthRows.forEach((
      /** @type {readonly UnvalidatedFields<'resourceId'>[]} */ row,
      /** @type {number} */ index
    ) => {
      if (!files.linearSeqs[index] || !Array.isArray(row)) return;
      const depth = depthFileRows?.[index] || row
        .map((ref) => ref?.resourceId
          ? (resolveResourceFile
              ? resolveResourceFile(ref.resourceId)
              : resourceAsLegacyFile(resources, ref.resourceId))
          : null);
      files.linearSeqs[index].depth = depth.length > 1 ? depth : (depth[0] || null);
    });
  }
  // Sessions 31–39 (repairLegacyTableRows) read these three tables as the
  // current writer writes them; the Load dialog lists the rows read so.
  /** @type {LoadRepair[]} */
  const legacyTableRepairs = [];
  /**
   * @template {{ repairs?: readonly TableRepair[] | null }} P
   * @param {string} table
   * @param {{ resourceId: unknown }} ref Its `resourceId` is set: each caller checks it.
   * @param {(text: string, options: { legacyRows: boolean }) => P} parse
   * @returns {P}
   */
  const readLegacyTable = (table, ref, parse) => {
    // `resourceTextFromRef` read a set id so: a property key is the id's string form.
    const parsed = readSessionTable(table, () => parse(
      decodeCanonicalResourceText(resources, String(ref.resourceId)), { legacyRows: repairLegacyTableRows }
    ));
    legacyTableRepairs.push(...(parsed.repairs || []).map(({ row, repair }) => loadRepair({
      code: `table-row-${repair}`, path: table, effect: repair === 'joined' ? 'changed' : 'dropped', line: row
    })));
    return parsed;
  };
  const defaultColorsRef = options.colors?.defaultColorsFile || options.colors?.defaultColors;
  let projectedDefaultColors = deferResourceContent
    && storedConfig?.colors
    && typeof storedConfig.colors === 'object'
    && !Array.isArray(storedConfig.colors)
    ? storedConfig.colors
    : {};
  if (defaultColorsRef?.resourceId) {
    files.d_color = resolveResourceFile
      ? resolveResourceFile(defaultColorsRef.resourceId)
      : resourceAsLegacyFile(resources, defaultColorsRef.resourceId);
    if (!deferResourceContent) {
      projectedDefaultColors = readLegacyTable('default-colors', defaultColorsRef, parseColorTable).colors;
    }
  }
  const colorTableRef = options.colors?.colorTableFile || options.colors?.colorTable;
  let projectedSpecificRules = deferResourceContent && Array.isArray(storedConfig?.rules)
    ? storedConfig.rules
    : [];
  if (colorTableRef?.resourceId) {
    files.t_color = resolveResourceFile
      ? resolveResourceFile(colorTableRef.resourceId)
      : resourceAsLegacyFile(resources, colorTableRef.resourceId);
    if (!deferResourceContent) {
      projectedSpecificRules = readSessionTable('specific-colors', () => parseSpecificRules(
        resourceTextFromRef(resources, colorTableRef)
      )).rules.map(({ fromFile: _fromFile, ...rule }) => rule);
    }
  }
  let projectedWhitelist = deferResourceContent && Array.isArray(storedConfig?.whitelist)
    ? storedConfig.whitelist
    : [];
  const whitelistRef = editorTableRef(options, 'whitelist');
  if (whitelistRef) {
    files.whitelist = resolveResourceFile
      ? resolveResourceFile(whitelistRef.resourceId)
      : resourceAsLegacyFile(resources, whitelistRef.resourceId);
    if (!deferResourceContent) {
      projectedWhitelist = readLegacyTable('label-whitelist', whitelistRef, parseWhitelistRules).rules;
    }
  }
  const qualifierPriorityRef = editorTableRef(options, 'priority');
  let projectedPriorityRules = deferResourceContent
    && Array.isArray(storedConfig?.qualifierPriorityRules)
    ? storedConfig.qualifierPriorityRules
    : [];
  if (qualifierPriorityRef) {
    files.qualifier_priority = resolveResourceFile
      ? resolveResourceFile(qualifierPriorityRef.resourceId)
      : resourceAsLegacyFile(resources, qualifierPriorityRef.resourceId);
    if (!deferResourceContent) {
      projectedPriorityRules = readLegacyTable('qualifier-priority', qualifierPriorityRef, parsePriorityRules).rules;
    }
  }
  const visibilityRef = editorTableRef(options, 'visibility');
  const projectedFeatureVisibilityRules = !deferResourceContent && visibilityRef
    ? readSessionTable('feature-visibility', () => parseFeatureVisibilityRules(
        resourceTextFromRef(resources, visibilityRef)
      )).rules
    : [];
  const labelOverrideRef = editorTableRef(options, 'labelOverrides');
  const projectedLabelOverrideRows = !deferResourceContent && labelOverrideRef
    ? readSessionTable('label-overrides', () => parseLabelOverrideTsv(
      resourceTextFromRef(resources, labelOverrideRef)
    )).map((row) => ({
        recordId: row.recordId,
        featureType: row.featureType,
        qualifier: row.qualifier,
        valueRegex: row.valueRegex,
        labelText: row.labelText
      }))
    : [];
  if (renderRequest.mode === 'circular' && Array.isArray(options.conservationBlastFiles)) {
    files.c_conservation_blasts = options.conservationBlastFiles
      .map((/** @type {UnvalidatedFields<'resourceId'>} */ ref) => ref?.resourceId
        ? (resolveResourceFile
            ? resolveResourceFile(ref.resourceId)
            : resourceAsLegacyFile(resources, ref.resourceId))
        : null)
      .filter(Boolean);
    const storedConservationSource = String(
      storedConfig?.circularConservation?.source || ''
    ).trim().toLowerCase();
    if (
      webMetadata.conservationBlastSource === 'losat-cache' ||
      storedConservationSource === 'losat'
    ) {
      files.c_conservation_blasts_source = 'losat-cache';
    }
  }
  if (renderRequest.mode === 'circular' && Array.isArray(options.conservationFastaFiles)) {
    files.c_conservation_sequence_sources = options.conservationFastaFiles.map((/** @type {UnvalidatedFields<'resourceId'>} */ ref) => (
      ref?.resourceId
        ? (resolveResourceFile
            ? resolveResourceFile(ref.resourceId)
            : resourceAsLegacyFile(resources, ref.resourceId))
        : null
    ));
  } else if (
    renderRequest.mode === 'circular'
    && Array.isArray(webMetadata.conservationSequenceSources)
  ) {
    files.c_conservation_sequence_sources = webMetadata.conservationSequenceSources.map(
      (/** @type {unknown} */ resourceId) => (resourceId
        ? (resolveResourceFile
            ? resolveResourceFile(resourceId)
            : resourceAsLegacyFile(resources, resourceId))
        : null)
    );
  }
  if (
    renderRequest.mode === 'circular'
    && Array.isArray(webMetadata.conservationLosatFastaSources)
  ) {
    files.c_conservation_fastas = webMetadata.conservationLosatFastaSources
      .map((/** @type {unknown} */ resourceId) => (
        resourceId
          ? (resolveResourceFile
              ? resolveResourceFile(resourceId)
              : resourceAsLegacyFile(resources, resourceId))
          : null
      ));
  }
  // The tables the request reads. A Web or CLI binding names the file a table
  // was imported from; the writer stored the table derived from it (OV-367).
  const requestTableFiles = Object.fromEntries(REQUEST_TABLE_FILE_SLOTS.map((slot) => [slot, files[slot] || null]));
  Object.assign(files, applyWebFileBindings(
    files,
    webMetadata,
    resources,
    { resolveResourceFile, sessionResourceTable, adoptCanonicalPayloads }
  ));
  const explicitOverrides = Object.fromEntries(
    Object.entries(options.configOverrides || {}).filter(
      ([, value]) => value !== null && value !== undefined
    )
  );
  const currentTrackDefaults = trackDefaultsForMode(renderRequest.mode);
  const sparseConfigDefaults = {
    [CONFIG_OVERRIDE_PATHS.showGc]: currentTrackDefaults.gc,
    [CONFIG_OVERRIDE_PATHS.showSkew]: currentTrackDefaults.skew,
    ...(renderRequest.mode === 'linear'
      ? {
          [CONFIG_OVERRIDE_PATHS.linearAxisStrokeColor]:
            modeProfile('linear').linearAxisColor
        }
      : {})
  };
  const canonicalOverrides = {
    ...(
      options.config === null || options.config === undefined
        ? sparseConfigDefaults
        : {}
    ),
    ...projectFullConfigOverrides(options.config, renderRequest.mode),
    ...projectExplicitConfigOverrides(explicitOverrides, renderRequest.mode)
  };
  const overrides = {
    ...projectCanonicalConfigOverrides(canonicalOverrides, renderRequest.mode),
    ...projectLegacyFlatConfigOverrides(explicitOverrides)
  };
  const projectedLinearLabelPlacement = renderRequest.mode === 'linear'
    ? requireCurrentLinearLabelPlacement(overrides.label_placement)
    : 'auto';
  const projectedLinearTrackLayout = renderRequest.mode === 'linear'
    ? requireCurrentLinearTrackLayout(overrides.linear_track_layout)
    : 'middle';
  const comparisonHeight = classifyOptionalPositiveNumber(overrides.comparison_height);
  if (renderRequest.mode === 'linear' && comparisonHeight.status === 'invalid') {
    throw diagnosticError('INPUT_INVALID', { field: 'match_height', reason: 'POSITIVE_OR_AUTO' });
  }
  const tracks = options.tracks || {};
  const projectedCircularTrackSlots = renderRequest.mode === 'circular'
    ? (Array.isArray(tracks.circularTrackSlots)
        ? tracks.circularTrackSlots.map((/** @type {unknown} */ slot) => (
          slot && typeof slot === 'object' && !Array.isArray(slot)
              ? projectCurrentCanonicalCircularSlot(slot)
              : parseCircularTrackSlotSpecs(
                  [slot],
                  options.dinucleotide || 'GC',
                  overrides.track_type || 'tuckin'
                )[0]
          ))
        : [])
    : [];
  const projectedLinearTrackSlots = renderRequest.mode === 'linear'
    // An array of rows migrates to an array of rows.
    ? /** @type {LinearTrackSlot[]} */ (migrateLinearTrackSlotsToCurrentSchema(
        parseLinearTrackSlotSpecs(tracks.linearTrackSlots, options.dinucleotide || 'GC'),
        linearTrackSlotSchemaVersion
      ))
    : [];
  const depthMetadataFields = [
    options.depthTrackLabels,
    options.depthTrackColors,
    options.depthTrackHeights,
    options.depthTrackLargeTickIntervals,
    options.depthTrackSmallTickIntervals,
    options.depthTrackTickFontSizes
  ];
  const referencedDepthTrackWidth = [
    ...projectedCircularTrackSlots,
    ...projectedLinearTrackSlots
  ].reduce((width, slot) => {
    if (slot?.renderer !== 'depth') return width;
    const trackIndex = parseDepthTrackIndexIdentity(
      slot?.params?.track_index ?? 0,
      `Depth slot '${slot?.id || ''}' track_index`
    );
    return Math.max(width, trackIndex + 1);
  }, 0);
  const projectedDepthTrackCount = canonicalDepth
    ? canonicalDepth.tracks.length
    : Math.max(
        0,
        ...depthRows.map((/** @type {unknown} */ row) => Array.isArray(row) ? row.length : 0),
        ...depthMetadataFields.map((values) => Array.isArray(values) ? values.length : 0),
        referencedDepthTrackWidth
      );
  validateTrackSlotBindingInvariants(projectedCircularTrackSlots, {
    modeLabel: 'Circular',
    layoutKind: 'circular',
    supportedRenderers: CIRCULAR_TRACK_RENDERERS,
    supportedSides: ['inside', 'outside', 'overlay'],
    anchorlessRenderers: ['ticks', 'spacer'],
    depthTrackCount: projectedDepthTrackCount
  });
  validateTrackSlotBindingInvariants(projectedLinearTrackSlots, {
    modeLabel: 'Linear',
    layoutKind: 'linear',
    supportedRenderers: LINEAR_TRACK_RENDERERS,
    supportedSides: ['above', 'below', 'overlay'],
    anchorlessRenderers: ['spacer'],
    depthTrackCount: projectedDepthTrackCount
  });
  validateProjectedDepthSources(depthRows, projectedDepthTrackCount);
  if (
    !canonicalDepth &&
    renderRequest.mode === 'circular' &&
    Array.isArray(options.depthTrackHeights) &&
    options.depthTrackHeights.some((/** @type {unknown} */ height) => height !== null && height !== undefined)
  ) {
    throw new Error('diagramOptions.depthTrackHeights must contain only null values for Circular requests.');
  }
  const projectedDepthTracks = canonicalDepth?.tracks || Array.from(
    { length: projectedDepthTrackCount },
    (_, index) => ({
      label: String(options.depthTrackLabels?.[index] ?? (index === 0 ? 'Depth' : `Depth ${index + 1}`)),
      color: String(options.depthTrackColors?.[index] || (index === 0 ? overrides.depth_color : '') || '#4A90E2'),
      height: projectOptionalNumber(options.depthTrackHeights?.[index], { field: 'height' }),
      large_tick_interval: projectOptionalNumber(options.depthTrackLargeTickIntervals?.[index], { field: 'large_tick_interval' }),
      small_tick_interval: projectOptionalNumber(options.depthTrackSmallTickIntervals?.[index], { field: 'small_tick_interval' }),
      tick_font_size: projectOptionalNumber(options.depthTrackTickFontSizes?.[index], { field: 'tick_font_size' })
    })
  );
  const circularPresentationRecord = (
    renderRequest.mode === 'circular' && grouping === 'single' && records.length === 1
  ) ? records[0] : null;
  const circularPresentationRegion = circularPresentationRecord?.region || null;
  const form = {
    prefix: projectedOutputPrefix(
      renderRequest,
      grouping,
      records,
      outputPrefixes,
      webMetadata.circularOutputPrefixExplicit
    ),
    plot_title: options.plotTitle || '',
    // A Linear request has no Circular grouping; keep the fresh default.
    multi_record_canvas: renderRequest.mode === 'circular'
      ? grouping === 'grid'
      : WEB_UX_PROFILE.circular.gridByDefault,
    circular_record_selector: circularPresentationRecord
      ? (canonicalRecordSelector(circularPresentationRecord) || '')
      : '',
    circular_region_start: circularPresentationRegion?.start ?? null,
    circular_region_end: circularPresentationRegion?.end ?? null,
    circular_reverse: Boolean(
      circularPresentationRegion?.reverseComplement ||
      circularPresentationRecord?.presentation?.reverseComplement
    ),
    circular_record_label: circularPresentationRecord?.presentation?.label || '',
    circular_record_subtitle: circularPresentationRecord?.presentation?.subtitle || '',
    suppress_gc: renderRequest.mode === 'circular' ? overrides.show_gc === false : false,
    suppress_skew: renderRequest.mode === 'circular' ? overrides.show_skew === false : false,
    show_gc: renderRequest.mode === 'linear' ? Boolean(overrides.show_gc) : false,
    show_skew: renderRequest.mode === 'linear' ? Boolean(overrides.show_skew) : false,
    show_depth: Boolean(overrides.show_depth),
    separate_strands: Boolean(overrides.strandedness),
    labels_mode: renderRequest.mode === 'circular' ? circularLabelsMode(overrides.label_scope) : 'none',
    show_labels_linear: renderRequest.mode === 'linear'
      ? (overrides.label_scope || 'none')
      : 'none',
    track_type: overrides.track_type || 'tuckin',
    linear_track_layout: projectedLinearTrackLayout,
    show_scale: overrides.show_scale !== false,
    scale_style: overrides.scale_style || 'bar',
    align_center: Boolean(overrides.align_center),
    keep_definition_left_aligned: Boolean(overrides.keep_definition_left_aligned),
    linear_ruler_on_axis: Boolean(overrides.linear_ruler_on_axis),
    normalize_length: Boolean(overrides.normalize_length),
    species: renderRequest.mode === 'circular' ? (options.species || '') : '',
    strain: renderRequest.mode === 'circular' ? (options.strain || '') : ''
  };
  const projectedFeatureShapes = normalizeFeatureRenderingMap(options.featureShapes || {});
  if (
    !Object.prototype.hasOwnProperty.call(projectedFeatureShapes, 'repeat_region') &&
    Array.isArray(options.selectedFeaturesSet) &&
    options.selectedFeaturesSet.includes('repeat_region')
  ) {
    projectedFeatureShapes.repeat_region = defaultFeatureRendering('repeat_region');
  }
  const sparseFeatureTypes = MODE_DEFAULT_FEATURE_TYPES;
  const currentComparisonDefaults = comparisonFiltersForMode(renderRequest.mode);
  const sparseComparisonDefaults = {
    ...currentComparisonDefaults,
    alignmentLength: currentComparisonDefaults.alignment_length
  };
  const adv = {
    features: options.selectedFeaturesSet ?? [...sparseFeatureTypes],
    feature_shapes: projectedFeatureShapes,
    arrow_head_length_ratio: arrowHeadLengthRatioForState(
      overrides.arrow_head_length_ratio
    ),
    arrow_shaft_width_ratio: normalizeArrowShaftWidthRatio(
      overrides.arrow_shaft_width_ratio
    ),
    font_family: overrides.font_family ?? null,
    nt: options.dinucleotide || 'GC',
    window_size: options.window ?? null,
    step_size: options.step ?? null,
    label_rendering: overrides.label_rendering || 'auto',
    circular_label_placement: renderRequest.mode === 'circular'
      ? (overrides.circular_label_placement || 'horizontal')
      : 'horizontal',
    label_placement: projectedLinearLabelPlacement,
    circular_label_spacing: renderRequest.mode === 'circular'
      ? (overrides.circular_label_spacing ?? null)
      : null,
    linear_label_spacing: renderRequest.mode === 'linear'
      ? (overrides.linear_label_spacing ?? null)
      : null,
    plot_title_font_size: options.plotTitleFontSize ?? overrides.plot_title_font_size ?? null,
    def_font_size: renderRequest.mode === 'circular'
      ? (overrides.circular_definition_font_size ?? null)
      : (overrides.linear_definition_font_size ?? null),
    circular_definition_interval: renderRequest.mode === 'circular' ? (overrides.circular_definition_interval ?? null) : null,
    label_font_size: overrides.label_font_size ?? null,
    label_rotation: renderRequest.mode === 'linear' ? (overrides.label_rotation ?? null) : null,
    block_stroke_width: overrides.block_stroke_width ?? null,
    block_stroke_color: overrides.block_stroke_color ?? null,
    line_stroke_width: overrides.line_stroke_width ?? null,
    line_stroke_color: overrides.line_stroke_color ?? null,
    axis_stroke_width: renderRequest.mode === 'circular'
      ? (overrides.circular_axis_stroke_width ?? null)
      : (overrides.linear_axis_stroke_width ?? null),
    axis_stroke_color: renderRequest.mode === 'circular'
      ? (overrides.circular_axis_stroke_color ?? null)
      : (overrides.linear_axis_stroke_color ?? null),
    legend_box_size: overrides.legend_box_size ?? null,
    legend_font_size: overrides.legend_font_size ?? null,
    circular_grouping_intent: renderRequest.mode === 'circular'
      ? grouping
      : 'auto',
    multi_record_size_mode: projectedCircularSizeMode,
    multi_record_min_radius_ratio: renderRequest.layout?.multiRecordMinRadiusRatio ?? 0.55,
    multi_record_column_gap_ratio: renderRequest.layout?.multiRecordColumnGapRatio ?? 0.10,
    multi_record_row_gap_ratio: renderRequest.layout?.multiRecordRowGapRatio ?? 0.05,
    center_reserved_radius: tracks.centerReservedRadius ?? null,
    resolve_overlaps: Boolean(overrides.resolve_overlaps),
    feature_overlap_tolerance_bp: overrides.feature_overlap_tolerance_bp ?? 0,
    comparison_height: renderRequest.mode === 'linear' && comparisonHeight.status === 'valid'
      ? comparisonHeight.value
      : null,
    feature_height: overrides.default_cds_height ?? null,
    gc_height: overrides.gc_height ?? null,
    track_axis_gap: overrides.linear_track_axis_gap ?? null,
    linear_definition_line_styles: overrides.linear_definition_line_styles || {},
    linear_show_replicon: Boolean(overrides.linear_definition_show_replicon),
    // A Circular request has no Linear display values; keep the fresh Auto.
    linear_accession_visibility: renderRequest.mode !== 'linear' ? 'auto'
      : overrides.linear_definition_show_accession !== false ? 'show' : 'hide',
    linear_length_visibility: renderRequest.mode !== 'linear' ? 'auto'
      : overrides.linear_definition_show_length !== false ? 'show' : 'hide',
    keep_full_definition_with_plot_title: Boolean(options.keepFullDefinitionWithPlotTitle),
    gc_content_mode: overrides.gc_content_mode || 'deviation',
    gc_content_min_percent: overrides.gc_content_min_percent ?? 0,
    gc_content_max_percent: overrides.gc_content_max_percent ?? 100,
    gc_content_show_axis: overrides.gc_content_show_axis !== false,
    gc_content_show_ticks: overrides.gc_content_show_ticks !== false,
    gc_content_tick_interval: overrides.gc_content_large_tick_interval ?? null,
    gc_content_small_tick_interval: overrides.gc_content_small_tick_interval ?? null,
    gc_content_tick_font_size: overrides.gc_content_tick_font_size ?? null,
    depth_color: overrides.depth_color || '#4A90E2',
    depth_height: overrides.depth_height ?? null,
    depth_min: overrides.depth_min ?? null,
    depth_max: overrides.depth_max ?? null,
    depth_normalize: Boolean(overrides.depth_normalize),
    depth_show_axis: overrides.depth_show_axis !== false,
    depth_show_ticks: overrides.depth_show_ticks !== false,
    depth_large_tick_interval: overrides.depth_large_tick_interval ?? null,
    depth_small_tick_interval: overrides.depth_small_tick_interval ?? null,
    depth_tick_font_size: overrides.depth_tick_font_size ?? null,
    depth_share_axis: Boolean(overrides.depth_share_axis),
    depth_window_size: options.depthWindow ?? null,
    depth_step_size: options.depthStep ?? null,
    depth_tracks: projectedDepthTracks,
    min_bitscore: options.bitscore ?? sparseComparisonDefaults.bitscore,
    evalue: options.evalue === null || options.evalue === undefined
      ? String(sparseComparisonDefaults.evalue)
      : String(options.evalue),
    identity: options.identity ?? sparseComparisonDefaults.identity,
    alignment_length:
      options.alignmentLength ?? sparseComparisonDefaults.alignmentLength,
    scale_stroke_color: overrides.scale_stroke_color ?? null,
    ruler_label_color: overrides.scale_label_color ?? null,
    scale_stroke_width: overrides.scale_stroke_width ?? null,
    scale_font_size: overrides.scale_font_size ?? null,
    ruler_label_font_size: overrides.ruler_label_font_size ?? null,
    scale_interval: overrides.scale_interval ?? null,
    tick_label_font_size: overrides.tick_label_font_size ?? null,
    ...Object.fromEntries(Object.entries(DRAFT_FIELD_OF_FLAT_CONFIG_KEY).map(([key, field]) => [
      field, renderRequest.mode === 'circular' ? (overrides[key] ?? null) : null
    ])),
    pairwise_match_style: overrides.pairwise_match_style || options.pairwiseMatchStyle || 'ribbon',
    circular_track_slots_enabled: renderRequest.mode === 'circular' && Array.isArray(tracks.circularTrackSlots),
    circular_track_slots_schema_version: 4,
    circular_track_slots_axis_index: tracks.circularTrackAxisIndex ?? null,
    circular_track_slots: projectedCircularTrackSlots,
    linear_track_slots_enabled: renderRequest.mode === 'linear' && Array.isArray(tracks.linearTrackSlots),
    linear_track_slots_schema_version: LINEAR_TRACK_SLOT_SCHEMA_VERSION,
    linear_track_slots_axis_index: renderRequest.mode === 'linear'
      ? (tracks.linearTrackAxisIndex ?? null)
      : null,
    linear_track_slots: projectedLinearTrackSlots,
    multi_record_positions: (renderRequest.layout?.multiRecordPositions || []).map((/** @type {unknown} */ token) => {
      const split = String(token).lastIndexOf('@');
      return { selector: String(token).slice(0, split), row: Number(String(token).slice(split + 1)) };
    })
  };
  const presentationOwnsLinearRows = records.some(
    (record) => record.presentation?.gridRow != null
  );
  const linearLayoutRows = presentationOwnsLinearRows
    ? records.map((record, index) => {
        const row = Number(record.presentation?.gridRow) || index + 1;
        const sourceColumn = sourceRecords[
          normalizedRecordOrdering.sourceIndexByProjectedIndex[index]
        ]?.presentation?.gridColumn;
        return {
          uid: files.linearSeqs[index]?.uid || '', row,
          ...(initializeCliInputs && storedConfig == null
            && sourceRecords[normalizedRecordOrdering.sourceIndexByProjectedIndex[index]]?.cardinality === 'exactly_one'
            ? { canonicalCardinality: 'exactly_one' } : {}),
          // An integer is a number.
          ...(initializeCliInputs && storedConfig == null
            && Number.isInteger(sourceColumn) && /** @type {number} */ (sourceColumn) > 0
            ? { canonicalRow: row, canonicalColumn: sourceColumn }
            : {})
        };
      })
    : (renderRequest.layout?.multiRecordPositions || []).map((/** @type {unknown} */ token, /** @type {number} */ index) => {
        const split = String(token).lastIndexOf('@');
        return {
          uid: files.linearSeqs[index]?.uid || '',
          row: Number(String(token).slice(split + 1)) || index + 1
        };
      });
  const linearLayout = renderRequest.mode === 'linear'
    ? {
        enabled: presentationOwnsLinearRows ||
          (renderRequest.layout?.multiRecordPositions || []).length > 0,
        recordGap: renderRequest.layout?.recordGapPx ?? DEFAULT_LINEAR_RECORD_GAP_PX,
        rows: linearLayoutRows,
        recordTranslations,
        similarityAlignment
      }
    : undefined;
  const projectedBlacklistText = Array.isArray(overrides.label_blacklist)
    ? overrides.label_blacklist.join(', ')
    : String(overrides.label_blacklist || '');
  let inputType = records[0]?.source?.kind === 'gffFasta' ? 'gff' : 'gb';
  // CLI rendering can normalize GFF records to GenBank resources. Without a
  // saved Web draft, initialize the input selector from the original bindings.
  if (storedConfig == null && explicitBindings) {
    const source = renderRequest.mode === 'circular'
      ? { gb: files.c_gb, gff: files.c_gff, fasta: files.c_fasta }
      : files.linearSeqs[0];
    if (!source?.gb && source?.gff && source?.fasta) inputType = 'gff';
  }
  // Legend and plot-title positions are layout preferences, not form fields:
  // the committed request sets the slot of its own mode and grouping.
  const layoutPreferences = createDefaultLayoutPreferences();
  updateActiveLayoutPreference(
    layoutPreferences,
    renderRequest.mode,
    renderRequest.mode === 'circular' && grouping === 'grid',
    {
      legend: options.output?.legend || 'right',
      plotTitlePosition: options.output?.plotTitlePosition
        || (renderRequest.mode === 'linear' ? 'bottom' : 'none')
    }
  );
  return {
    mode: renderRequest.mode,
    inputType,
    files,
    requestTableFiles,
    layoutPreferences,
    config: {
      form,
      adv,
      ...(cliDraftWithoutComparison ? cliProteinSettings : projectedProteinPipeline?.config || {}),
      ...cliProteinPlan,
      colors: projectedDefaultColors,
      colorsAreOverrides: true,
      palette: options.colors?.defaultColorsPalette || 'default',
      rules: projectedSpecificRules,
      qualifierPriorityRules: projectedPriorityRules,
      filterMode: projectedWhitelist.length > 0
        ? 'Whitelist'
        : (projectedBlacklistText.trim() ? 'Blacklist' : 'None'),
      whitelist: projectedWhitelist,
      blacklistText: projectedBlacklistText,
      linearRecordLayout: linearLayout,
      annotationSets: draftAnnotationSetsOfRequest(options.annotations?.sets),
      recordDisplayDrafts: records.flatMap((record, index) => (record.display?.isCircular != null || record.display?.startCoordinate != null) ? [{
        sourceUid: renderRequest.mode === 'linear' ? String(files.linearSeqs[index]?.uid || record.recordKey) : 'circular',
        selector: record.selector?.kind === 'recordIndex' ? `#${record.selector.index + 1}` : '#1',
        recordId: record.selector?.kind === 'recordId' ? record.selector.value : '',
        topologyOverride: record.display?.isCircular ?? null,
        startCoordinate: record.display?.startCoordinate ?? null,
        reverseComplementOverride: null,
        anchorIntent: null
      }] : []),
      featurePlacementOverrides: featureDraftMap(canonicalFeaturePlacements(
        options.featurePlacements || [], renderRequest.mode
      )),
      circularConservation: renderRequest.mode === 'circular'
        ? projectCircularConservationConfig(options, files)
        : undefined
    },
    semanticFeatureState: {
      featureVisibilityManualRules: projectedFeatureVisibilityRules,
      // A request written without the Web draft (CLI, Python) carries its
      // per-feature edits only as these rows.
      featureOverrides: featureDraftMap(
        projectedFeatureOverrides.map((row) => ({ ...row, labelSourceText: null }))
      ),
      labelOverrideRows: projectedLabelOverrideRows
    },
    pipelineState: projectedProteinPipeline
      ? {
          generatedProteinComparison:
            projectedProteinPipeline.generatedProteinComparison
        }
      : null,
    legacyTableRepairs
  };
};
// Output fields that name files, not the figure; a batch request has one output per record.
const PUBLICATION_OUTPUT_ONLY_FIELDS = new Set(['prefix', 'formats', 'overwrite', 'artifactFilename']);
const PUBLICATION_COMPARISON_FILTER_FIELDS = new Set(['evalue', 'bitscore', 'identity', 'alignmentLength']);
// A request key outside this table reads no default: a JSON value is never an inherited member.
/** @type {Readonly<Record<string, JsonValue>>} */
const PUBLICATION_OPTION_DEFAULTS = { dinucleotide: 'GC', keepFullDefinitionWithPlotTitle: false, conservationReference: 'auto', 'objects.features.arrow_geometry.head_length_ratio': 'auto', 'objects.features.arrow_geometry.shaft_width_ratio': 1, 'objects.scale.show': true };
/**
 * A resource of a request under comparison: the descriptor a Session stores or
 * the writer built (no reader has checked it), or a backing that reads its bytes.
 * @typedef {RawSessionResourceDescriptor & Pick<SessionResourceSource, 'readBytes'>} PublicationResource
 */
/** @typedef {Readonly<Record<string, PublicationResource>> | null | undefined} PublicationResources */
/**
 * What one comparison keeps while it walks a request.
 * @typedef {object} PublicationContext
 * @property {{ path: string, kind: string, decodedPayloadSha256: string }[]} bindings
 * @property {boolean} normalize
 * @property {Map<unknown, Promise<string>>} resourceIdentities Digests by a resource's base64 text.
 * @property {boolean} ignoreComparisonFilters
 */
/**
 * @param {PublicationResource | undefined} resource
 * @param {string} id
 * @returns {Promise<Uint8Array>}
 */
const publicationBytes = async (resource, id) => {
  if (!resource || typeof resource !== 'object' || Array.isArray(resource)) throw new Error(`Canonical request resource '${id}' is missing.`);
  if (typeof resource.readBytes === 'function') return resource.readBytes();
  if (resource.data instanceof Uint8Array) return resource.data;
  if (resource.encoding === 'base64') { const binary = atob(String(resource.data || '')), bytes = new Uint8Array(binary.length); for (let index = 0; index < binary.length; index += 1) bytes[index] = binary.charCodeAt(index); return bytes; }
  if (typeof resource.data === 'string') return textToBytes(resource.data);
  throw new Error(`Canonical request resource '${id}' has no decodable payload.`);
};
const TSV_NUMBER = /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
// The tables a CLI Session stores with a header; the files the Web writes have
// none (OV-266, OV-367). Each compares by kind and rows
// without the header; default colors compare as a set, the others keep their
// rule order.
/** @type {Readonly<Record<string, readonly [kind: string, header: string]>>} */
const PUBLICATION_TABLES = Object.freeze({
  '$.diagramOptions.colors.defaultColors': ['default-colors', 'feature_type\tcolor'],
  '$.diagramOptions.colors.colorTable': ['color-table', 'feature_type\tqualifier_key\tvalue\tcolor\tcaption'],
  '$.diagramOptions.featureVisibilityTable': ['feature-visibility-table', 'record_id\tfeature_type\tqualifier\tvalue\taction'],
  '$.diagramOptions.labelWhitelistTable': ['label-whitelist-table', 'feature_type\tqualifier\tkeyword'],
  '$.diagramOptions.qualifierPriorityTable': ['qualifier-priority-table', 'feature_type\tpriorities'],
  '$.diagramOptions.labelOverrideTable': ['label-override-table', 'record_id\tfeature_type\tqualifier\tvalue\tlabel_text']
});
/** @param {string} path A path of the request from `$`, so never an inherited member's name. */
const publicationTable = (path) => PUBLICATION_TABLES[path.replace(/\.resourceId$/, '')] || null;
/**
 * @param {PublicationResource | undefined} resource
 * @param {string} id
 * @param {boolean} normalize
 * @param {string} path
 */
const normalizedPublicationBytes = async (resource, id, normalize, path) => {
  const bytes = await publicationBytes(resource, id);
  const table = publicationTable(path);
  if (table) {
    const [kind, header] = table, sorted = kind === 'default-colors';
    const rows = bytesToText(bytes).split(/\r?\n/).map((row) => (sorted ? row.trim() : row))
      .filter((row) => row.trim() && row !== header);
    return textToBytes(`${(sorted ? rows.sort() : rows).join('\n')}\n`);
  }
  if (!normalize) return bytes;
  // `publicationBytes` above refused a missing resource.
  if (!path.startsWith('$.comparisons') || /** @type {PublicationResource} */ (resource).kind !== 'canonical-tsv') return bytes;
  const rows = bytesToText(bytes).trimEnd().split(/\r?\n/).map((row, index) => index === 0
    ? row : row.split('\t').map((cell) => {
      const number = Number(cell);
      if (TSV_NUMBER.test(cell) && Number.isFinite(number)) return String(Number(number.toPrecision(15)));
      // The CLI replay writes booleans through pandas as True/False (OV-413).
      return cell === 'True' || cell === 'False' ? cell.toLowerCase() : cell;
    }).join('\t'));
  return textToBytes(`${rows.join('\n')}\n`);
};
/**
 * @param {PublicationResources} resources
 * @param {JsonValue} id
 * @param {boolean} normalize
 * @param {string} path
 * @param {PublicationContext['resourceIdentities']} cache
 */
const publicationResourceIdentity = async (resources, id, normalize, path, cache) => {
  const resourceId = String(id || '').trim();
  if (!resourceId) throw new Error('Canonical request contains an empty resourceId.');
  const table = publicationTable(path);
  // Own key (OV-614); a resource the table lacks is refused by `publicationBytes` below.
  const resource = ownEntry(resources, resourceId), key = resource?.encoding === 'base64' && !table && (!normalize || !path.startsWith('$.comparisons') || resource.kind !== 'canonical-tsv') ? resource.data : null;
  const kind = table ? table[0] : String(resource?.kind || '');
  let digest = key ? cache.get(key) : null;
  if (!digest) { digest = normalizedPublicationBytes(resource, resourceId, normalize, path).then((bytes) => sha256Hex(bytes)); if (key) cache.set(key, digest); }
  return { kind, decodedPayloadSha256: await digest };
};
/**
 * @param {JsonValue} value
 * @param {PublicationResources} resources
 * @param {PublicationContext} context
 * @param {string} [path]
 * @returns {Promise<JsonValue>}
 */
const canonicalizePublicationValue = async (value, resources, context, path = '$') => {
  if (Array.isArray(value)) return Promise.all(value.map((entry, index) =>
    canonicalizePublicationValue(entry, resources, context, `${path}[${index}]`)));
  if (!value || typeof value !== 'object') return value;
  /** @type {{ [key: string]: JsonValue }} */
  const output = {};
  const resourceReference = Object.hasOwn(value, 'resourceId');
  for (const key of Object.keys(value).sort()) {
    if (value[key] === null || Object.is(value[key], PUBLICATION_OPTION_DEFAULTS[key]) || key === 'plotTitleFontSize'
      || (path === '$.diagramOptions.configOverrides' && key === 'labels.filtering.blacklist_keywords'
        && Array.isArray(value[key]) && value[key].length === 0)
      || (path === '$.diagramOptions.configOverrides' && key === 'canvas.feature_overlap_tolerance_bp' && value[key] === 0)
      // A request may state a default the other writer omits (OV-398): the
      // definition-line weight `normal`, which a Web draft holds as null, and
      // a feature rendering, which a CLI request states only when changed.
      || (path === '$.diagramOptions.configOverrides' && LINEAR_DEFINITION_WEIGHT_PATHS.has(key) && value[key] === 'normal')
      || (path === '$.diagramOptions.featureShapes' && value[key] === defaultFeatureRendering(key))
      || (path === '$.diagramOptions' && key === 'featureShapes' && value[key] && typeof value[key] === 'object'
        && Object.entries(value[key]).every(([type, rendering]) => rendering === defaultFeatureRendering(type)))
      || (path === '$.diagramOptions.config.canvas' && key === 'feature_overlap_tolerance_bp' && value[key] === 0)
      || (context.ignoreComparisonFilters && path === '$.diagramOptions' && PUBLICATION_COMPARISON_FILTER_FIELDS.has(key))
      || (/^\$\.output(?:\[\d+\])?$/.test(path) && PUBLICATION_OUTPUT_ONLY_FIELDS.has(key))) continue;
    if (resourceReference && ['encoding', 'representation'].includes(key)) continue;
    const childPath = `${path}.${key}`;
    if (key === 'resourceId') {
      const identity = await publicationResourceIdentity(resources, value[key], context.normalize, childPath, context.resourceIdentities);
      context.bindings.push({ path: childPath, ...identity });
      output[key] = identity;
    } else output[key] = await canonicalizePublicationValue(value[key], resources, context, childPath);
  }
  return output;
};
/**
 * @typedef {{ path: string, expected: JsonValue | undefined, actual: JsonValue | undefined }} PublicationDifference
 */
/** @typedef {{ [key: string]: JsonValue }} JsonEntries An object, or an array read through its own keys. */
/**
 * @param {JsonValue} expected
 * @param {JsonValue} actual
 * @param {string} [path]
 * @returns {PublicationDifference | null}
 */
const firstPublicationDiff = (expected, actual, path = '$') => {
  if (Object.is(expected, actual)) return null;
  if (!expected || !actual || typeof expected !== 'object' || typeof actual !== 'object'
      || Array.isArray(expected) !== Array.isArray(actual)) return { path, expected, actual };
  for (const key of new Set([...Object.keys(expected), ...Object.keys(actual)])) {
    const childPath = Array.isArray(expected) ? `${path}[${key}]` : `${path}.${key}`;
    if (!Object.hasOwn(expected, key) || !Object.hasOwn(actual, key)) return {
      path: childPath, expected: /** @type {JsonEntries} */ (expected)[key], actual: /** @type {JsonEntries} */ (actual)[key] };
    const difference = firstPublicationDiff(
      /** @type {JsonEntries} */ (expected)[key], /** @type {JsonEntries} */ (actual)[key], childPath
    );
    if (difference) return difference;
  }
  return null;
};
/**
 * What a promotion of schemas 5 to 7 reads a legacy similarity alignment with.
 * `unreadAlignment`, when given, receives an alignment that cannot be read,
 * which is then left unset; without it the promotion throws.
 * @typedef {object} RequestPromotion
 * @property {Parameters<typeof materializeLegacySimilarityAlignment>[0]['featureCatalog']} [featureCatalog]
 * @property {Parameters<typeof materializeLegacySimilarityAlignment>[0]['legacyOrthogroupState']} [legacyOrthogroupState]
 * @property {(error: unknown) => void} [unreadAlignment]
 */
/**
 * @param {Record<string, any>} request A canonical request of schema 5 to 9.
 * @param {RequestPromotion} [promotion]
 * @returns {CanonicalRenderRequest}
 */
export const promoteCanonicalRenderRequestToCurrent = (request, promotion = {}) => {
  const promoted = cloneCanonicalJsonValue(request);
  if (promoted.schema === CANONICAL_REQUEST_SCHEMA) return promoted;
  if (![5, 6, 7, 8].includes(promoted.schema)) {
    throw new Error('Only canonical renderRequest schemas 5, 6, 7, and 8 can be promoted to schema 9.');
  }
  if (promoted.schema < 8) promoteRenderRequestToSchema8(promoted, promotion);
  // Schema 9 adds per-feature edits by source identity (design Q4).
  promoted.schema = CANONICAL_REQUEST_SCHEMA;
  promoted.diagramOptions = { ...promoted.diagramOptions, featureOverrides: [] };
  return promoted;
};

/**
 * A stored request of schema 5, 6 or 7 as the promotion reads and edits it (a JSON clone).
 * `readStoredRenderRequest` checked the root only: a plain object with an integer schema,
 * a `records` array, and `comparisons` absent or an array. Nothing below the root is checked;
 * the repair follows the promotion.
 * @typedef {{
 *   schema: number,
 *   mode?: unknown,
 *   records: unknown[],
 *   comparisons?: unknown[],
 *   layout?: unknown,
 *   diagramOptions?: unknown
 * }} StoredRequestBeforeSchema8
 */
/**
 * @param {StoredRequestBeforeSchema8} promoted Edited in place.
 * @param {RequestPromotion} promotion
 */
const promoteRenderRequestToSchema8 = (promoted, promotion) => {
  const { featureCatalog = null, legacyOrthogroupState = null, unreadAlignment } = promotion;
  const sourceSchema = promoted.schema;
  const storedLayout = /** @type {UnvalidatedFields<'multiRecordPositions' | 'recordGapPx'>} */ (promoted.layout);
  const linearRows = promoted.mode === 'linear'
    ? (Array.isArray(storedLayout?.multiRecordPositions) ? storedLayout.multiRecordPositions : []).map((/** @type {unknown} */ token) => {
        const split = String(token).lastIndexOf('@');
        return Number(String(token).slice(split + 1)) || null;
      })
    : [];
  promoted.schema = 8;
  // A malformed record is left to the repair that follows the promotion.
  (promoted.records || []).forEach((value, index) => {
    if (!value || typeof value !== 'object') return;
    // An object (checked above): its properties are read as unvalidated values and written.
    const record = /** @type {{ display?: unknown, cardinality?: unknown, selector?: unknown, region?: unknown,
      presentation?: unknown }} */ (value);
    if (sourceSchema < 7) record.display = { isCircular: null, startCoordinate: null };
    if (sourceSchema === 5) record.cardinality = promoted.mode === 'linear' &&
      !record.selector && !record.region
      ? 'all'
      : 'exactly_one';
    if (linearRows[index] && record.presentation && typeof record.presentation === 'object') {
      // An object (checked above).
      /** @type {{ gridRow?: unknown }} */ (record.presentation).gridRow = linearRows[index];
    }
  });
  if (sourceSchema < 7) {
    promoted.diagramOptions = {
      .../** @type {UnvalidatedFields<string>} */ (promoted.diagramOptions), featurePlacements: []
    };
  }
  if (promoted.mode === 'linear') {
    // A record that is not an object is the repair's, as above.
    const records = (promoted.records || []).filter((/** @type {unknown} */ record) => isPlainObject(record));
    // A view, and the settings are edited below; a value that is no object has no own `alignOrthogroupFeature`.
    const generated = /** @type {{ settings?: { alignOrthogroupFeature?: unknown } | null } | null | undefined} */ (
      (promoted.comparisons || []).find(
        (comparison) => /** @type {UnvalidatedFields<'kind'>} */ (comparison)?.kind === 'generatedProteinComparison'
      )
    );
    const rawLegacy = generated?.settings?.alignOrthogroupFeature;
    const readAlignment = () => canonicalSimilarityAlignment(materializeLegacySimilarityAlignment({
      target: rawLegacy,
      records,
      featureCatalog,
      legacyOrthogroupState
    }), records.map((record) => record.recordKey),
    'renderRequest.layout.similarityAlignment');
    /** @type {ReturnType<typeof canonicalSimilarityAlignment> | null} */
    let similarityAlignment = null;
    if (rawLegacy !== null && rawLegacy !== undefined) {
      try {
        similarityAlignment = readAlignment();
      } catch (error) {
        if (!unreadAlignment) throw error;
        unreadAlignment(error);
      }
    }
    if (generated?.settings) delete generated.settings.alignOrthogroupFeature;
    const hadLayout = Object.keys(promoted.layout || {}).length > 0;
    if (hadLayout || similarityAlignment !== null) {
      const recordGapPx = Number(storedLayout?.recordGapPx ?? DEFAULT_LINEAR_RECORD_GAP_PX);
      promoted.layout = {
        recordGapPx: Number.isFinite(recordGapPx) && recordGapPx >= 0 ? recordGapPx : DEFAULT_LINEAR_RECORD_GAP_PX,
        multiRecordPositions: null,
        recordTranslations: records.map((record) => ({
          recordKey: String(record.recordKey || ''),
          x: 0,
          y: 0
        })),
        similarityAlignment
      };
    } else {
      promoted.layout = {};
    }
  }
  return promoted;
};

/**
 * What a stored request is read with: the catalog and orthogroup state of its
 * Result set (a legacy alignment of schemas 5 to 7 names their features), the
 * resource ids of the Session, and where the request is in the file.
 * @typedef {object} StoredRequestContext
 * @property {unknown} featureCatalog
 * @property {unknown} legacyOrthogroupState
 * @property {ReadonlySet<string>} resourceIds
 * @property {string} path
 * @property {boolean} [holdComparisons] Load holds the request's comparisons
 *   (a saved comparison kept read-only, or one that needs a decision): their
 *   unresolved record references are what the Comparisons notice after Load
 *   shows (receipt 〔今までの知らせ〕), so they are not dropped. A reference to a
 *   resource the Session lacks is dropped as any other is: no Save can write
 *   it (OV-899).
 */
/**
 * The committed Session: the request that drew the shown Result, with its
 * resources and Web file bindings. Generate and the automatic rerender commit
 * the schema 9 request they built; Load commits the request the Session
 * stored, read to schema 9 (`readStoredRenderRequest`); History and the Load
 * rollback put a captured one back.
 * @typedef {object} CommittedCanonicalSession
 * @property {CanonicalRenderRequest} renderRequest
 * @property {Readonly<Record<string, unknown>>} resources
 * @property {Readonly<Record<string, unknown>>} webFiles
 */
/**
 * A stored request that cannot be read at all (K22): its Session loads as
 * settings only, or its Result set is left out.
 * @typedef {{ kind: 'request-dropped', path: string, reason: string }} StoredRequestDropped
 */
/** @typedef {'unknown-resource' | 'unknown-record' | 'unknown-slot'} StoredReferenceReason */

const STORED_REQUEST_ROOTS = /** @type {const} */ ({ circular: 'CanonicalCircularRenderRequest', linear: 'CanonicalLinearRenderRequest' });

/**
 * Per mode, from the description: the places a value of the other diagram
 * type is read as this one's (a placeholder that holds only null or no item,
 * and an option of the other type that this one holds at its default), the
 * unions picked by their `kind`, and the places that have a default (a repair
 * puts it there in place of a value it cannot read, so each is a part).
 * @param {string} rootType
 * @param {string} otherOptionsType
 */
const storedRequestPlaces = (rootType, otherOptionsType) => {
  const { types } = CANONICAL_REQUEST_DESCRIPTION;
  /** @type {Map<string, unknown>} */
  const otherType = new Map();
  /** @type {Set<string>} */
  const tagged = new Set();
  /** @type {Set<string>} */
  const defaulted = new Set();
  const otherOptions = /** @type {import('../utils/contract-shape.js').ObjectShape} */ (types[otherOptionsType]);
  /** @param {ContractShape} shape @returns {ContractShape} */
  const resolve = (shape) => (shape.kind === 'ref' ? resolve(types[shape.name]) : shape);
  /**
   * @param {ContractShape} shape @param {string} pattern @param {number} depth
   * @param {boolean} [option] the shape is one option of the place's union
   */
  const walk = (shape, pattern, depth, option = false) => {
    const resolved = resolve(shape);
    if (depth > 40) return;
    if (resolved.kind === 'null' && pattern && !option) otherType.set(pattern, null);
    // A list only the other type fills (a Circular request's `comparisons`).
    if (resolved.kind === 'tuple' && resolved.items.length === 0 && pattern && !option) otherType.set(pattern, []);
    if (resolved.kind === 'literal' && resolved.values.length === 1 && /^diagramOptions\.[^.]+$/.test(pattern)) {
      const other = otherOptions.fields[pattern.slice('diagramOptions.'.length)];
      if (other && JSON.stringify(resolve(other)) !== JSON.stringify(resolved)) otherType.set(pattern, resolved.values[0]);
    }
    if (resolved.kind === 'union') {
      const options = resolved.options.map(resolve);
      if (options.every((option) => option.kind === 'object' && Object.hasOwn(option.fields, 'kind'))) tagged.add(pattern);
      options.forEach((item) => walk(item, pattern, depth + 1, true));
    } else if (resolved.kind === 'array') walk(resolved.items, `${pattern}[*]`, depth + 1);
    else if (resolved.kind === 'tuple') resolved.items.forEach((item) => walk(item, `${pattern}[*]`, depth + 1));
    else if (resolved.kind === 'record') walk(resolved.values, `${pattern}.*`, depth + 1);
    else if (resolved.kind === 'object') {
      Object.entries(resolved.fields).forEach(([key, field]) => {
        const place = pattern ? `${pattern}.${key}` : key;
        if (resolved.defaults && Object.hasOwn(resolved.defaults, key)) defaulted.add(place);
        walk(field, place, depth + 1);
      });
    }
  };
  walk(types[rootType], '', 0);
  return { otherType, tagged, defaulted, otherOptionKeys: new Set(Object.keys(otherOptions.fields)) };
};
const STORED_REQUEST_PLACES = {
  circular: storedRequestPlaces(STORED_REQUEST_ROOTS.circular, 'CanonicalLinearDiagramOptions'),
  linear: storedRequestPlaces(STORED_REQUEST_ROOTS.linear, 'CanonicalCircularDiagramOptions')
};
/**
 * The optional places without a default that a repair drops as one part
 * (rule 6); a place with a default is a part too. Each can go alone.
 */
const STORED_REQUEST_PARTS = new Set(['diagramOptions.tracks', 'diagramOptions.annotations']);
/** The conservation ring columns, read by index together (B16). */
const CONSERVATION_RING_COLUMNS = ['conservationBlastFiles', 'conservationLabels', 'conservationColors', 'conservationFastaFiles']
  .map((key) => `diagramOptions.${key}`);
/** The comparisons Load can hold: Linear comparisons and the conservation rings. */
const HELD_COMPARISON_PLACE = /^(comparisons\[\*\]|diagramOptions\.conservation[A-Za-z]+\[\*\])/;
const RECORD_INDEX_PLACE = /^comparisons\[\*\](\.pairs\[\*\])?\.(query|subject)RecordIndex$/;
const RECORD_KEY_PLACE = /^layout\.(recordTranslations\[\*\]|similarityAlignment\.(reference|records\[\*\](\.anchor)?))\.recordKey$/;

/**
 * The request policy of the repair leaf: what the description cannot say
 * about a stored request of `mode`.
 * @param {Readonly<Record<string, unknown>>} request The promoted request (schema 9 shape, unchecked).
 * @param {'circular' | 'linear'} mode
 * @param {ReadonlySet<string>} resourceIds
 * @param {boolean} holdComparisons The record references of the comparisons are kept.
 * @returns {RepairPolicy<StoredReferenceReason>}
 */
const storedRequestPolicy = (request, mode, resourceIds, holdComparisons) => {
  const places = STORED_REQUEST_PLACES[mode];
  const records = Array.isArray(request.records) ? request.records : [];
  const recordKeys = new Set(records.map((/** @type {unknown} */ record) => (isPlainObject(record) ? record.recordKey : undefined)));
  const options = isPlainObject(request.diagramOptions) ? request.diagramOptions : {};
  const tracks = isPlainObject(options.tracks) ? options.tracks : {};
  const slotKey = mode === 'linear' ? 'linearTrackSlots' : 'circularTrackSlots';
  const slots = Array.isArray(tracks[slotKey]) ? tracks[slotKey] : null;
  // An axis index names a slot by position, so no slot is dropped alone then.
  const axisIndex = tracks[mode === 'linear' ? 'linearTrackAxisIndex' : 'circularTrackAxisIndex'];
  return {
    positionalArrays: [
      // An alignment plan decides every record: one decision lost loses the plan.
      'records', 'layout.multiRecordPositions', 'diagramOptions.depthTracks[*].source', 'layout.similarityAlignment.records',
      ...(axisIndex === null || axisIndex === undefined ? [] : [`diagramOptions.tracks.${slotKey}`])
    ],
    partBoundary: (pattern) => STORED_REQUEST_PARTS.has(pattern) || places.defaulted.has(pattern),
    parallelArrays: mode === 'circular' ? [CONSERVATION_RING_COLUMNS] : [],
    lossless: (pattern) => (places.otherType.has(pattern) ? { value: places.otherType.get(pattern) } : null),
    discriminator: (pattern) => (pattern === '' ? 'mode' : places.tagged.has(pattern) ? 'kind' : null),
    reference: (pattern, value) => {
      if (/(^|\.)(resourceId|gffResourceId|fastaResourceId)$/.test(pattern)) {
        return typeof value === 'string' && !resourceIds.has(value) ? 'unknown-resource' : null;
      }
      if (holdComparisons && HELD_COMPARISON_PLACE.test(pattern)) return null;
      if (RECORD_INDEX_PLACE.test(pattern)) {
        return Number.isInteger(value) && Number(value) >= 0 && Number(value) < records.length ? null : 'unknown-record';
      }
      if (RECORD_KEY_PLACE.test(pattern)) return recordKeys.has(value) ? null : 'unknown-record';
      if (pattern === 'diagramOptions.tracks.linearTrackAxisIndex' && value !== null) {
        return slots && Number.isInteger(value) && Number(value) >= 0 && Number(value) < slots.length ? null : 'unknown-slot';
      }
      return null;
    }
  };
};

/**
 * The Load rows of the leaf's rows. A setting of the other diagram type had
 * no effect on the figure and gets none (〔出さないもの〕).
 * @param {readonly ContractRepair<StoredReferenceReason>[]} rows
 * @param {ReadonlySet<string>} otherOptionKeys
 * @param {readonly string[]} droppedSlots The slots of dropped rings (B16), by their stored place.
 * @returns {LoadRepair[]}
 */
const storedRequestRepairs = (rows, otherOptionKeys, droppedSlots) => {
  /** @type {Map<string, number>} */
  const droppedEntries = new Map();
  /** @type {LoadRepair[]} */
  const repairs = [];
  rows.forEach((row) => {
    const reference = row.reason === 'unknown-resource' || row.reason === 'unknown-record' || row.reason === 'unknown-slot';
    if (row.rule === 2 && /^diagramOptions\.[^.]+$/.test(row.pattern)
      && otherOptionKeys.has(row.pattern.slice('diagramOptions.'.length))) return;
    if (reference) {
      repairs.push(loadRepair({ code: 'request-reference-dropped', path: row.path, effect: row.effect }));
    } else if (row.rule === 6) {
      repairs.push(loadRepair({ code: 'request-part-dropped', path: row.path, effect: row.effect }));
    } else if (row.rule === 2 || row.rule === 5) {
      // One row per list or object: an element's path without its index, an
      // unknown key's path without the key.
      const list = row.path.replace(row.rule === 5 ? /\[\d+\]$/ : /\.[^.[\]]+$/, '');
      droppedEntries.set(list, (droppedEntries.get(list) ?? 0) + 1);
    } else {
      const used = row.used;
      repairs.push(loadRepair({
        code: 'setting-replaced', path: row.path, effect: row.effect,
        ...(typeof used === 'string' || typeof used === 'number' || typeof used === 'boolean' ? { used: String(used) } : {})
      }));
    }
  });
  droppedSlots.forEach((slot) => {
    const list = slot.replace(/\[\d+\]$/, '');
    droppedEntries.set(list, (droppedEntries.get(list) ?? 0) + 1);
  });
  droppedEntries.forEach((count, path) => {
    repairs.push(loadRepair({ code: 'setting-entries-dropped', path, effect: 'dropped', count }));
  });
  return repairs;
};

/** The element places of the ring columns: a rule-5 row there drops a column. */
const RING_COLUMN_ELEMENTS = new Set(CONSERVATION_RING_COLUMNS.map((pattern) => `${pattern}[*]`));
/**
 * A slot parameter as Python's `int()` reads it, or null.
 * @param {unknown} value
 */
const slotIndexParam = (value) => (
  typeof value === 'number' && Number.isFinite(value) ? Math.trunc(value)
    : typeof value === 'string' && /^\s*[+-]?\d+\s*$/.test(value) ? Number(value) : null);

/**
 * Keeps each managed conservation slot on its ring when the repair dropped
 * ring columns (B16). Python binds a `sequence_conservation` slot to the ring
 * whose index is its `source_index` (`gbdraw/api/diagram.py`), so the slot of
 * a dropped ring is left out and a later slot's `source_index` and
 * `track_index` go down by the rings dropped before it, in their stored type;
 * the axis boundary goes down by the slots left out before it. Copies what it
 * changes: the repaired request may share the stored objects.
 * @param {Record<string, unknown>} request The repaired Circular request.
 * @param {readonly ContractRepair<StoredReferenceReason>[]} rows
 * @param {string} path
 * @returns {{ request: Record<string, unknown>, droppedSlots: string[] }}
 */
const keepConservationSlotsOnTheirRings = (request, rows, path) => {
  const droppedRings = rows.filter((row) => row.rule === 5 && RING_COLUMN_ELEMENTS.has(row.pattern))
    .map((row) => Number(/\[(\d+)\]$/.exec(row.path)?.[1]));
  const options = request.diagramOptions;
  const tracks = isPlainObject(options) ? options.tracks : null;
  const slots = isPlainObject(tracks) ? tracks.circularTrackSlots : null;
  if (!droppedRings.length || !isPlainObject(options) || !isPlainObject(tracks) || !Array.isArray(slots)) {
    return { request, droppedSlots: [] };
  }
  /** @param {number} index */
  const droppedBefore = (index) => new Set(droppedRings.filter((ring) => ring < index)).size;
  /** @param {unknown} stored @param {number} lower */
  const lowered = (stored, lower) => {
    const index = slotIndexParam(stored);
    return index === null ? stored : typeof stored === 'number' ? index - lower : String(index - lower);
  };
  /** @type {string[]} */
  const droppedSlots = [];
  /** @type {unknown[]} */
  const kept = [];
  let axis = tracks.circularTrackAxisIndex;
  slots.forEach((/** @type {unknown} */ slot, index) => {
    const ringSlot = isPlainObject(slot) && slot.renderer === 'sequence_conservation' ? slot : null;
    const params = ringSlot && isPlainObject(ringSlot.params) ? ringSlot.params : null;
    const ring = params ? slotIndexParam(params.source_index) : null;
    if (!ringSlot || !params || ring === null || !droppedBefore(ring + 1)) {
      kept.push(slot);
    } else if (droppedRings.includes(ring)) {
      droppedSlots.push(`${path}.diagramOptions.tracks.circularTrackSlots[${index}]`);
      if (typeof axis === 'number' && index < axis) axis -= 1;
    } else {
      const lower = droppedBefore(ring);
      kept.push({
        ...ringSlot,
        params: {
          ...params,
          source_index: lowered(params.source_index, lower),
          ...(Object.hasOwn(params, 'track_index') ? { track_index: lowered(params.track_index, lower) } : {})
        }
      });
    }
  });
  if (kept.length === slots.length && kept.every((slot, index) => slot === slots[index])) return { request, droppedSlots };
  return {
    request: {
      ...request,
      diagramOptions: {
        ...options,
        tracks: { ...tracks, circularTrackSlots: kept, ...(Object.hasOwn(tracks, 'circularTrackAxisIndex') ? { circularTrackAxisIndex: axis } : {}) }
      }
    },
    droppedSlots
  };
};

/**
 * The one reader of a stored render request: the only place in `js/` where an
 * unchecked value becomes a `CanonicalRenderRequest`. (1) An integer schema
 * outside the supported set is refused as a Session-format error; (2) schemas
 * 1 and 2 are lifted to schema 5 (`liftLegacyRenderRequest`, PR-D), with a row
 * for what schema 5 cannot carry; (3) schemas 5 to 8 are promoted once, with
 * the Load context; (4) the request is repaired against its description. What
 * cannot be read at its root, a schema that is no integer included, is a
 * dropped request (K22).
 * @param {unknown} value
 * @param {StoredRequestContext} context
 * @returns {Repaired<CanonicalRenderRequest> | StoredRequestDropped}
 */
export const readStoredRenderRequest = (value, context) => {
  const { path } = context;
  if (!isPlainObject(value)) return { kind: 'request-dropped', path, reason: 'wrong-kind' };
  // No writer wrote a schema that is no integer: the request is unreadable
  // (K22, OV-863), not another format.
  if (!Number.isInteger(value.schema)) return { kind: 'request-dropped', path: `${path}.schema`, reason: 'wrong-kind' };
  if (!SUPPORTED_CANONICAL_REQUEST_SCHEMA_SET.has(/** @type {number} */ (value.schema))) {
    throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'SESSION_FORMAT' });
  }
  const lifted = liftLegacyRenderRequest(value, path);
  /** @type {LoadRepair[]} */
  const repairs = [...lifted.repairs];
  const stored = lifted.value;
  // The promotion walks the records and comparisons; lists of another kind
  // leave nothing to read, as the repair would find at the root; nor does a
  // schema 1 or 2 request the lift cannot read (no options object or mode).
  if (stored.schema === 1 || stored.schema === 2
    || !Array.isArray(stored.records) || (Object.hasOwn(stored, 'comparisons') && !Array.isArray(stored.comparisons))) {
    return { kind: 'request-dropped', path, reason: 'wrong-kind' };
  }
  const promoted = stored.schema === CANONICAL_REQUEST_SCHEMA ? stored : promoteCanonicalRenderRequestToCurrent(stored, {
    // The casts: the legacy alignment reader checks the catalog and the
    // orthogroup state where it reads them (materializeLegacySimilarityAlignment).
    featureCatalog: /** @type {RequestPromotion['featureCatalog']} */ (context.featureCatalog),
    legacyOrthogroupState: /** @type {RequestPromotion['legacyOrthogroupState']} */ (context.legacyOrthogroupState),
    unreadAlignment: () => repairs.push(loadRepair({
      code: 'request-part-dropped', path: `${path}.layout.similarityAlignment`, effect: 'dropped'
    }))
  });
  const mode = promoted.mode === 'linear' ? 'linear' : promoted.mode === 'circular' ? 'circular' : null;
  if (mode === null) return { kind: 'request-dropped', path: `${path}.mode`, reason: 'not-in-set' };
  const repaired = repairByDescription(
    promoted, { type: STORED_REQUEST_ROOTS[mode], path }, CANONICAL_REQUEST_DESCRIPTION,
    storedRequestPolicy(promoted, mode, context.resourceIds, context.holdComparisons === true)
  );
  const root = repaired.rows.find((row) => row.pattern === '');
  if (repaired.value === undefined || root) return { kind: 'request-dropped', path, reason: String(root?.reason ?? 'wrong-kind') };
  // The Linear layout's cross-record rules (translations cover the drawn
  // records, the alignment plan names them) are its owner's, which Load and
  // Generate run; a layout it refuses leaves no request to draw (C-OD20).
  try {
    canonicalLinearRecordLayout(/** @type {Record<string, unknown>} */ (repaired.value));
  } catch (error) {
    // The owner refuses a layout with a plain `Error` and no code; anything
    // else (a TypeError, a typed error) is a defect and is not read as a refusal.
    if (!(error instanceof Error) || Object.getPrototypeOf(error) !== Error.prototype || Object.hasOwn(error, 'code')) throw error;
    return { kind: 'request-dropped', path: `${path}.layout`, reason: 'unknown-record' };
  }
  const rebound = keepConservationSlotsOnTheirRings(
    /** @type {Record<string, unknown>} */ (repaired.value), repaired.rows, path
  );
  // The leaf returns the stored object itself when nothing changed (RK-3).
  // The cast: repairByDescription checked it against the request description.
  return {
    value: /** @type {CanonicalRenderRequest} */ (/** @type {unknown} */ (rebound.request)),
    repairs: [
      ...repairs,
      ...storedRequestRepairs(repaired.rows, STORED_REQUEST_PLACES[mode].otherOptionKeys, rebound.droppedSlots)
    ]
  };
};

/**
 * A stored request as Gallery publication reads it: through the one reader,
 * refusing any repair (a published Session is written by the current writer,
 * so a row is a defect of its source) and any request the reader drops.
 * Returns a copy the caller may change.
 * @param {unknown} request
 * @param {Omit<StoredRequestContext, 'path' | 'holdComparisons'>} context
 * @returns {CanonicalRenderRequest}
 */
export const readPublishedRenderRequest = (request, context) => {
  const read = readStoredRenderRequest(request, { ...context, path: 'renderRequest' });
  if ('kind' in read || read.repairs.length > 0) {
    throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
  }
  return /** @type {CanonicalRenderRequest} */ (cloneCanonicalJsonValue(read.value));
};

/**
 * @param {unknown} left
 * @param {unknown} right
 */
const sameCanonicalValue = (left, right) => (
  JSON.stringify(left) === JSON.stringify(right)
);

/** @param {unknown} selector */
const materializedRecordSelector = (selector) => {
  const match = String(selector || '').match(/^#([1-9]\d*)$/);
  if (!match) {
    throw new Error('Target record materialization requires an exact record selector.');
  }
  return { kind: 'recordIndex', index: Number(match[1]) - 1 };
};

/**
 * A member of a record collection that a feature rotation targets: the
 * committed row of one source record (`targetForFeature`, record display),
 * which `materializeCanonicalRecordCollection` checks.
 * @typedef {object} RecordTransformTargetMember
 * @property {unknown} [canonicalRecordKey]
 * @property {unknown} [recordKey]
 * @property {unknown} [selector]
 * @property {unknown} [recordLength]
 * @property {unknown} [committedDisplay]
 * @property {unknown} [committedReverseComplement]
 */
/**
 * The record a feature rotation targets, as the record display's
 * `targetForFeature` port gives it. Its fields are unchecked:
 * `validateRecordTransformTarget` checks them against the committed request.
 * @typedef {object} RecordTransformTarget
 * @property {unknown} [scope]
 * @property {unknown} [recordKey]
 * @property {unknown} [canonicalRecordKey]
 * @property {unknown} [source]
 * @property {unknown} [effectiveCircular]
 * @property {unknown} [cropped]
 * @property {unknown} [recordLength]
 * @property {readonly (RecordTransformTargetMember | null | undefined)[]} [members]
 */
/**
 * The transform a feature rotation resolved, before it is checked.
 * @typedef {{ recordLength: unknown, startCoordinate: number | null, reverseComplement: boolean | null }} RecordTransformInput
 */
/**
 * @param {RecordTransformTarget} target
 * @param {RecordTransformInput} transform
 * @param {unknown} mode The mode of the committed request, compared as it is.
 */
const validateRecordTransformTarget = (target, transform, mode) => {
  if (!target || target.scope !== mode || typeof target.recordKey !== 'string'
    || !target.recordKey || typeof target.canonicalRecordKey !== 'string'
    || !target.canonicalRecordKey || !target.source || typeof target.source !== 'object') {
    throw new Error('Record transform target is stale or incomplete.');
  }
  if (target.cropped) {
    throw new Error('A cropped record cannot be rotated from a feature.');
  }
  if (target.effectiveCircular !== true) {
    throw new Error('Feature-based record rotation requires a circular record.');
  }
  const recordLength = Number(target.recordLength);
  if (!Number.isSafeInteger(recordLength) || recordLength < 1
    || Number(transform?.recordLength) !== recordLength) {
    throw new Error('The target record length changed after the feature popup opened.');
  }
  const startCoordinate = transform?.startCoordinate;
  if (startCoordinate === null || !Number.isSafeInteger(startCoordinate)
    || startCoordinate < 1 || startCoordinate > recordLength
    || typeof transform?.reverseComplement !== 'boolean') {
    throw new Error('Resolved record transform is invalid for the target record.');
  }
};

/**
 * @param {CanonicalRecord} record
 * @param {RecordTransformTarget} target
 * @param {number} recordIndex
 */
const materializeCanonicalRecordCollection = (record, target, recordIndex) => {
  const members = Array.isArray(target.members) ? target.members : [];
  const identities = new Set();
  const selectors = new Set();
  const materialized = members.map((member) => {
    const selector = materializedRecordSelector(member?.selector);
    const recordKey = String(member?.recordKey || '');
    if (!recordKey || identities.has(recordKey) || selectors.has(selector.index)
      || member?.canonicalRecordKey !== target.canonicalRecordKey
      || !Number.isSafeInteger(Number(member?.recordLength))
      || Number(member.recordLength) < 1
      || !member?.committedDisplay || typeof member.committedDisplay !== 'object'
      || typeof member?.committedReverseComplement !== 'boolean') {
      throw new Error('Committed record collection cannot be materialized safely.');
    }
    identities.add(recordKey);
    selectors.add(selector.index);
    return {
      ...cloneCanonicalJsonValue(record),
      recordKey,
      cardinality: 'exactly_one',
      selector,
      display: cloneCanonicalJsonValue(member.committedDisplay),
      presentation: {
        ...(cloneCanonicalJsonValue(record.presentation) || {}),
        reverseComplement: member.committedReverseComplement,
        gridRow: record.presentation?.gridRow ?? recordIndex + 1
      }
    };
  });
  if (materialized.length < 1
    || materialized.filter((entry) => entry.recordKey === target.recordKey).length !== 1) {
    throw new Error('Target record collection does not resolve to exactly one record.');
  }
  return materialized;
};

/**
 * @param {unknown} comparisons
 * @param {number} recordIndex
 * @param {number} expansion
 */
const shiftCanonicalComparisonIndexes = (comparisons, recordIndex, expansion) => (
  /** @type {unknown[]} */ (Array.isArray(comparisons) ? comparisons : []).map((comparison) => {
    const shifted = cloneCanonicalJsonValue(comparison);
    for (const field of ['queryRecordIndex', 'subjectRecordIndex']) {
      const index = Number(shifted?.[field]);
      if (Number.isInteger(index) && index > recordIndex) shifted[field] = index + expansion;
    }
    if (Array.isArray(shifted?.pairs)) {
      shifted.pairs = shifted.pairs.map((/** @type {Record<string, unknown>} */ pair) => {
        const next = { ...pair };
        for (const field of ['queryRecordIndex', 'subjectRecordIndex']) {
          const index = Number(next[field]);
          if (Number.isInteger(index) && index > recordIndex) next[field] = index + expansion;
        }
        return next;
      });
    }
    return shifted;
  })
);

/**
 * The committed request (Load, Generate and the rerender commit schema 9 only)
 * of a committed Session that holds its resources.
 * @param {CommittedCanonicalSession | null | undefined} committed
 * @returns {CanonicalRenderRequest}
 */
const requireCurrentCommittedRequest = (committed) => {
  const request = committed?.renderRequest;
  if (!request || !['circular', 'linear'].includes(request.mode)
    || !committed?.resources || typeof committed.resources !== 'object') {
    throw new Error('A current committed canonical Session is required.');
  }
  return request;
};

/**
 * Clone the last committed canonical Session and overlay one record transform.
 * No live form state participates in this projection.
 * @param {{ committed: CommittedCanonicalSession | null, target: RecordTransformTarget, transform: RecordTransformInput }} input
 */
export const projectCommittedRecordTransform = ({ committed, target, transform }) => {
  const request = requireCurrentCommittedRequest(committed);
  validateRecordTransformTarget(target, transform, request.mode);
  const matchingIndexes = request.records
    .map((record, index) => ({ record, index }))
    .filter(({ record }) => record?.recordKey === target.canonicalRecordKey);
  if (matchingIndexes.length !== 1) {
    throw new Error('Target record does not resolve uniquely in the committed request.');
  }
  const { record: sourceRecord, index: sourceIndex } = matchingIndexes[0];
  if (!sameCanonicalValue(sourceRecord.source, target.source)) {
    throw new Error('Target record resource binding is stale.');
  }
  if (sourceRecord.region || target.cropped) {
    throw new Error('A cropped record cannot be rotated from a feature.');
  }

  const candidate = cloneCanonicalJsonValue(committed);
  let records = candidate.renderRequest.records;
  let targetIndex = sourceIndex;
  let materialized = false;
  if (sourceRecord.cardinality === 'all') {
    const replacements = materializeCanonicalRecordCollection(
      sourceRecord,
      target,
      sourceIndex
    );
    records.splice(sourceIndex, 1, ...replacements);
    targetIndex = sourceIndex + replacements.findIndex(
      (entry) => entry.recordKey === target.recordKey
    );
    const expansion = replacements.length - 1;
    candidate.renderRequest.comparisons = shiftCanonicalComparisonIndexes(
      candidate.renderRequest.comparisons,
      sourceIndex,
      expansion
    );
    const depthFiles = candidate.renderRequest.diagramOptions?.depthTrackFiles;
    if (Array.isArray(depthFiles) && depthFiles.length === request.records.length) {
      depthFiles.splice(
        sourceIndex,
        1,
        ...replacements.map(() => cloneCanonicalJsonValue(depthFiles[sourceIndex]))
      );
    }
    const metadata = candidate.webFiles?.linearRecordMetadata;
    if (Array.isArray(metadata) && metadata.length === request.records.length) {
      const sourceMetadata = metadata[sourceIndex] || {};
      metadata.splice(sourceIndex, 1, ...replacements.map((entry) => ({
        ...cloneCanonicalJsonValue(sourceMetadata),
        recordKey: entry.recordKey
      })));
    }
    materialized = true;
  } else if (sourceRecord.cardinality !== 'exactly_one'
    || target.recordKey !== target.canonicalRecordKey) {
    throw new Error('Target record identity is stale.');
  }

  const candidateTarget = records[targetIndex];
  candidateTarget.display = {
    ...(candidateTarget.display || { isCircular: null }),
    startCoordinate: transform.startCoordinate
  };
  candidateTarget.presentation = {
    ...(candidateTarget.presentation || {}),
    reverseComplement: transform.reverseComplement
  };

  projectCanonicalSessionRequest({
    renderRequest: candidate.renderRequest,
    resources: candidate.resources,
    webFiles: candidate.webFiles || {},
    storedConfig: candidate.config || null,
    deferResourceContent: true
  });
  return {
    canonical: candidate,
    receipt: Object.freeze({
      recordKey: target.recordKey,
      canonicalRecordKey: target.canonicalRecordKey,
      recordIndex: targetIndex,
      materialized,
      startCoordinate: transform.startCoordinate,
      reverseComplement: transform.reverseComplement
    })
  };
};

const COMMITTED_EDITOR_TABLE_OPTIONS = Object.freeze({
  colors: ['colors'],
  ...EDITOR_TABLE_OPTIONS
});

/**
 * @param {unknown} value
 * @param {Set<string>} [ids]
 * @returns {Set<string>}
 */
const referencedResourceIds = (value, ids = new Set()) => {
  if (Array.isArray(value)) value.forEach((entry) => referencedResourceIds(entry, ids));
  else if (value && typeof value === 'object') {
    Object.entries(value).forEach(([key, entry]) => {
      if (key === 'resourceId' && typeof entry === 'string' && entry) ids.add(entry);
      else referencedResourceIds(entry, ids);
    });
  }
  return ids;
};

/**
 * Clone the last committed canonical Session and replace only the tables of
 * live editor intent: feature colors (rules and the applied palette), feature
 * visibility, and label overrides. A label reflow renders this request, so no
 * other draft setting reaches the Result before Generate (R1(c), N-16).
 * @param {{ committed: any, state: Record<string, any>, drawing: RequestDrawing }} input
 */
export const projectCommittedEditorIntent = ({
  committed,
  state,
  drawing
}) => {
  const request = committed?.renderRequest;
  requireCurrentCommittedRequest(committed);
  const candidate = {
    ...committed,
    renderRequest: cloneCanonicalJsonValue(request),
    resources: { ...committed.resources }
  };
  const options = candidate.renderRequest.diagramOptions || {};
  candidate.renderRequest.diagramOptions = options;
  /** @type {Readonly<Record<string, boolean>>} */
  const tables = {
    colors: true, visibility: true, whitelist: false, priority: false, labelOverrides: true
  };
  const retired = new Set();
  Object.entries(COMMITTED_EDITOR_TABLE_OPTIONS).forEach(([table, keys]) => {
    if (!tables[table]) return;
    keys.forEach((key) => {
      referencedResourceIds(options[key], retired);
      delete options[key];
    });
  });
  const stillReferenced = referencedResourceIds(candidate.renderRequest);
  retired.forEach((resourceId) => {
    if (!stillReferenced.has(resourceId)) delete candidate.resources[resourceId];
  });
  const resources = createResourceBuilder();
  const labelProjection = requestLabelProjection(state, drawing);
  // Per-feature edits are live editor intent too (design Q4).
  options.featureOverrides = requestFeatureOverrides(
    drawing.featureOverrides, candidate.renderRequest.records,
    { bulkLabelText: labelProjection.bulkLabelText }
  );
  // `state` is the app state (state.js), which holds the refs the tables read.
  addGeneratedTableResources(
    /** @type {RequestTableState} */ (state), drawing, resources, options, tables, labelProjection
  );
  Object.assign(candidate.resources, resources.resources);
  projectCanonicalSessionRequest({
    renderRequest: candidate.renderRequest,
    resources: candidate.resources,
    webFiles: candidate.webFiles || {},
    storedConfig: candidate.config || null,
    deferResourceContent: true
  });
  return candidate;
};

/**
 * Project only alignment-owned fields of the last committed canonical artifact.
 * @param {{
 *   committed: CommittedCanonicalSession | null, plan: unknown, translations: unknown,
 *   orientations: readonly { recordKey: string, reverseComplement: boolean }[]
 * }} input The plan and translations are cloned as they are.
 */
export const projectCommittedSimilarityAlignment = ({ committed, plan, translations, orientations }) => {
  const committedRequest = committed?.renderRequest;
  if (!committedRequest || committedRequest.mode !== 'linear' || !committed?.resources) {
    throw new Error('Alignment requires a committed canonical Linear artifact.');
  }
  /** @type {CanonicalLinearRenderRequest} */
  const request = cloneCanonicalJsonValue(committedRequest);
  const canonical = { ...committed, renderRequest: request };
  const keys = request.records.map(({ recordKey }) => recordKey);
  const byKey = new Map((orientations || []).map((entry) => [entry.recordKey, entry.reverseComplement]));
  if (byKey.size !== keys.length || orientations.length !== keys.length
    || keys.some(key => typeof byKey.get(key) !== 'boolean')) {
    throw new Error('Alignment orientation coverage changed.');
  }
  // The test above gave every record key a boolean.
  request.records.forEach(record => writeCanonicalRecordReverseComplement(record, /** @type {boolean} */ (byKey.get(record.recordKey))));
  request.layout = { recordGapPx: DEFAULT_LINEAR_RECORD_GAP_PX, multiRecordPositions: null, ...request.layout,
    similarityAlignment: cloneCanonicalJsonValue(plan ?? null),
    recordTranslations: cloneCanonicalJsonValue(translations) };
  projectCanonicalSessionRequest({ ...canonical, deferResourceContent: true });
  return canonical;
};

/**
 * @param {{ [key: string]: JsonValue }} request A canonical request of schema 5 to 9, as stored or as written.
 * @param {PublicationResources} resources
 */
const normalizePublicationRequestAliases = (request, resources) => {
  // A request of schema 5 to 8 is read by the one reader (no context: a legacy
  // alignment target refuses); any other is compared as written, with its aliases.
  const normalized = /** @type {Record<string, any>} */ (!/** @type {readonly unknown[]} */ ([5, 6, 7, 8]).includes(request?.schema)
    ? cloneCanonicalJsonValue(request)
    : readPublishedRenderRequest(request, {
        featureCatalog: null, legacyOrthogroupState: null, resourceIds: new Set(Object.keys(resources || {}))
      }));
  for (const comparison of normalized.comparisons || []) {
    if (comparison.kind === 'generatedProteinComparison' && comparison.settings) {
      comparison.settings.collinearInferOrthogroups = resolvePipelineCollinearInference(
        comparison.settings, comparison.mode
      );
    }
  }
  if (drawsNoComparison(normalized.comparisons || [])) normalized.comparisons = [];
  const colors = normalized.diagramOptions?.colors;
  if (colors) {
    colors.defaultColors = colors.defaultColors || colors.defaultColorsFile || null;
    colors.colorTable = colors.colorTable || colors.colorTableFile || null;
    delete colors.defaultColorsFile;
    delete colors.colorTableFile;
  }
  const options = normalized.diagramOptions;
  if (options) Object.entries(EDITOR_TABLE_OPTIONS).forEach(([table, [key, fileKey]]) => {
    options[key] = editorTableRef(options, /** @type {keyof typeof EDITOR_TABLE_OPTIONS} */ (table));
    delete options[fileKey];
  });
  return normalized;
};
/**
 * @param {JsonValue} request
 * @param {PublicationResources} resources
 * @param {boolean} normalize
 * @param {PublicationContext['resourceIdentities']} resourceIdentities
 */
const publicationRequestIdentity = async (request, resources, normalize, resourceIdentities) => {
  if (!request || typeof request !== 'object' || Array.isArray(request)) throw new Error('Canonical request equivalence requires a renderRequest object.');
  const normalizedRequest = normalizePublicationRequestAliases(request, resources);
  const conservation = normalizedRequest.diagramOptions?.conservationBlastFiles;
  /** @type {PublicationContext} */
  const context = { bindings: [], normalize, resourceIdentities,
    ignoreComparisonFilters: !(normalizedRequest.comparisons?.length || (Array.isArray(conservation) ? conservation.length : conservation)) };
  const canonical = await canonicalizePublicationValue(normalizedRequest, resources, context);
  return { canonical, digest: await sha256Hex(textToBytes(JSON.stringify(canonical))),
    resourceBindings: context.bindings };
};
/**
 * The two requests of a Gallery publication check, each with its resource table.
 * @typedef {object} PublicationComparisonInput
 * @property {JsonValue} expectedRequest
 * @property {PublicationResources} expectedResources
 * @property {JsonValue} actualRequest
 * @property {PublicationResources} actualResources
 * @property {boolean} [normalizeReplayGeneratedResources]
 */
/** @param {PublicationComparisonInput} input */
export const compareCanonicalRenderRequests = async (input) => {
  const normalize = input.normalizeReplayGeneratedResources === true;
  const resourceIdentities = new Map(), [expected, actual] = await Promise.all([publicationRequestIdentity(input.expectedRequest, input.expectedResources, normalize, resourceIdentities), publicationRequestIdentity(input.actualRequest, input.actualResources, normalize, resourceIdentities)]);
  const difference = firstPublicationDiff(expected.canonical, actual.canonical);
  return { equivalent: expected.digest === actual.digest, expected, actual, differences: difference ? [difference] : [] };
};
/** @param {PublicationComparisonInput} input */
export const assertCanonicalRenderRequestsEquivalent = async (input) => {
  const comparison = await compareCanonicalRenderRequests(input);
  if (comparison.equivalent) return comparison;
  const error = /** @type {Error & { comparison?: typeof comparison }} */ (new Error(`Gallery publication request differs at ${comparison.differences[0]?.path || '$'} (committed ${comparison.expected.digest}, rebuilt ${comparison.actual.digest}).`));
  error.comparison = comparison;
  throw error;
};
