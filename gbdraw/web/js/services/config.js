// @ts-check
import { diagnosticError, normalizeCaughtError, normalizeUserFacingError } from '../utils/error-normalization.js';
import { state, sessionOperationAvailability, normalizeLinearSeqList, collapseEmptyLinearSeqList } from '../state.js';
import { normalizePaletteColors, resolveColorToHex } from '../utils/color-utils.js';
import {
  captureRightDrawerState,
  resetRightDrawerState,
  restoreRightDrawerState
} from './right-drawer-state.js';
import { resetLayoutState, resetSettings as resetSettingsState } from './reset.js';
import { serializeCleanSvg } from './svg-serialization.js';
import { cloneJsonData, cloneJsonValue } from './json-clone.js';
import {
  applyCircularTrackOrderPlacements,
  clampCircularTrackAxisIndex,
  inferLegacyAxisIndexFromFeature,
  migrateLegacyCircularTrackSlot,
  migrateLegacyCircularTrackSlotSpec,
  parseCircularTrackSlotSpec,
  normalizeCircularTrackSlots
} from '../app/circular-track-slots.js';
import {
  applyLinearTrackOrderPlacements,
  clampLinearTrackAxisIndex,
  LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  migrateLinearTrackSlotsToCurrentSchema,
  normalizeLinearTrackSlots,
  resolveLinearTrackAxisIndex
} from '../app/linear-track-slots.js';
import {
  depthFileSlotsFromValue,
  depthTrackMatrixWidth,
  depthTrackSessionWidth,
  dropInvalidManagedDepthSlots,
  padDepthFileSlots,
  reconcileDepthTracksToFiles,
  representativeDepthFiles,
  syncDepthSlotLabels
} from './depth-track-state.js';
import {
  decodeDepthText,
  isEncodedDepthFileEntry
} from './depth-file-codec.js';
import {
  normalizeCollinearAnchorMode,
  normalizeCollinearSearchScope,
  normalizeOrthogroupMembershipMode
} from './losat-normalization.js';
import { normalizeDefinitionLineStyleState } from './definition-line-style-state.js';
import {
  migrateLegacyLinearLabelVisibility,
  requireLinearLabelVisibilityMode
} from './linear-label-visibility.js';
import { migrateLegacyOrthogroupMembers } from './legacy-similarity-alignment.js';
import {
  normalizeCircularPlotTitlePosition,
  normalizeLayoutPreferences,
  replaceLayoutPreferences,
  resolveActiveLayoutPreference,
  restoredLayoutPreferences
} from './layout-preferences.js';
import {
  serializeFeatureVisibilityRules,
  normalizeFeatureVisibilityRule,
  splitLegacyVisibilityRules
} from './feature-visibility.js';
import { canonicalFeatureOverrides, featureDraftMap } from './feature-placement.js';
import {
  ANNOTATION_TARGET_MIGRATION_NOTICE,
  FEATURE_EDIT_MIGRATION_WARNING,
  FEATURE_VISIBILITY_NARROWED_NOTICE,
  RENDERED_ID_FEATURE_EDIT_FIELDS,
  hasRenderedIdFeatureEdits,
  migrateSessionAnnotationTargets,
  migrateSessionFeatureEdits,
  migrateSessionFeaturePlacements
} from './feature-edit-migration.js';
import {
  buildSessionFeatureRecoveryPlan,
  extractSessionSourceFeatures
} from '../app/session-feature-metadata.js';
import {
  analyzeCatalogSequenceSourceCoverage,
  buildRestoredMatchSequenceSources,
  resolveCircularComparisonSequenceAvailability
} from './match-sequences.js';
import {
  CANONICAL_REQUEST_SCHEMA,
  buildCanonicalRenderRequest,
  canonicalLinearRecordLayout,
  legacyTableRowsNotice,
  managedConfigOverridePathsForMode,
  promoteCanonicalRenderRequestToCurrent,
  projectCanonicalSessionRequest,
  projectSettingsOnlySession,
  readCanonicalResourceRecordCount
} from './session-request.js';
import {
  createDefaultLinearComparisonPlan,
  linearComparisonEdgeKey,
  normalizeLinearComparisonPlan,
  reconcileLinearComparisonPlan,
  resolveLinearComparisonPlan
} from './linear-comparisons.js';
import { buildSessionResources as assembleSessionResources } from './session-resources.js';
import { ARTIFACT_SLOT_KEYS, createArtifactSlot } from './artifact-slot.js';
import {
  base64ToBytes,
  bytesToBase64,
  getSessionResourceSource,
  readFileBytes
} from './file-content-cache.js';
import {
  adoptCurrentSessionResources,
  isSessionResourceFileView
} from './session-resource-backing.js';
import { normalizeLogicalResults } from './result-normalization.js';
import {
  admitCurrentSessionResults,
  admitLegacyImportedResults,
  createCurrentSessionResultSource,
  createEmptySvgMutationPlan,
  createLegacyImportResultSource,
  isCommittedSvgResult
} from './svg-result-ingestion.js';
import {
  admitFeatureCatalog,
  featureStateFromCatalog,
  isAdoptedFeatureCatalog,
  validateFeatureCatalog,
  validateFeatureCatalogForImport
} from './feature-catalog.js';
import { migrateLegacyRecordDisplayDrafts } from '../app/record-display-options.js';
import {
  buildOrthogroupFeatureIndex,
  enrichFeaturesWithOrthogroups,
  normalizeOrthogroupDormantOverrides
} from './orthogroup-feature-metadata.js';
import {
  isResourceBackedCanonicalComparison,
  mapResourceBackedCanonicalComparison
} from './canonical-comparisons.js';
import {
  migrateLegacyLinearComparisonDraft,
  promoteGallerySessionToCurrent
} from './gallery-session-migration.js';
import {
  compressSessionData,
  confirmLargeSessionBlob
} from './session-file.js';
import { importSessionFile } from './session-import-client.js';
import { convertMainSessionComparisonFrames } from './main-session-comparison-frame.js';
import { downloadBlob } from './text-download.js';
import { normalizeAnnotationSets } from './annotation-state.js';
import { applySpecificRuleProvenance } from './specific-color-rules.js';
import { normalizeLegacyLegendEntryGroups } from './svg-result-normalization.js';
import {
  LOSAT_DERIVED_CACHE_SCHEMA,
  buildValidatedProteinIdentityIndex,
  classifyRawLosatCacheEntry,
  createLegacyProteinCandidateEnvelope,
  emptyProteinIdentityManifest,
  isCurrentRawLosatCacheEntry,
  isLosatDerivedCacheEntry,
  normalizeLegacyProteinCandidateEnvelope,
  releaseValidatedProteinIdentityIndex,
  serializableLegacyProteinCandidateEnvelope,
  validateDerivedProteinReferences,
  validateProteinRawEntryReferences,
  validateProteinIdentityManifest
} from '../app/losat-cache.js';
import {
  arrowHeadLengthRatioForState,
  defaultFeatureRendering,
  normalizeArrowShaftWidthRatio,
  normalizeFeatureRenderingMap
} from '../utils/feature-rendering.js';
import {
  adoptCurrentSessionDocument,
  adoptRuntimeCanonicalSession,
  isAdoptedCanonicalSession,
  hasBiologicalSessionInputs,
  isSettingsOnlySessionDocument,
  MODE_SCOPED_SESSION_VERSION,
  projectArtifactState,
  projectDocumentMetadata,
  projectWebOnlyEditorMetadata,
  TYPED_DRAFT_SESSION_VERSION,
  validateSessionAuthorityInventory
} from './session-authority.js';
import { assertSafeObjectKeysForImport } from './safe-object-keys.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from './runtime-test-hooks.js';
import { setResourcePayloadOwner } from './resource-payload-owner.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';
import {
  migratePersistedCircularMultiRecordSizeMode,
  migratePersistedLinearLabelPlacement,
  migratePersistedLinearTrackLayout,
  migratePersistedWebStateFieldNames,
  normalizeCurrentPairwiseMatchStyle,
  requireCurrentCircularMultiRecordSizeMode,
  requireCurrentLinearLabelPlacement,
  requireCurrentLinearTrackLayout,
  requireCurrentWebStateFieldNames
} from './current-option-values.js';
import {
  validateSimilarityAlignmentResetReceipt,
  CIRCULAR_TRACK_SLOT_SCHEMA_VERSION,
  CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS,
  createDefaultLosatpHitLimits,
  LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION,
  reconcileImportedLinearTypographyLink,
  validateCurrentWriterActiveConfig,
  validateImportedCircularTrackSlots,
  validateImportedLinearTrackSlots
} from './session-active-config-contract.js';
import {
  IMPORTED_COMPARISON_ACTIONS,
  IMPORTED_COMPARISON_DISPOSITIONS,
  classifyImportedComparisonIntent,
  createImportedComparisonIntentState,
  restoreImportedComparisonIntent,
  serializeImportedComparisonResolution
} from './imported-comparison-intent.js';

import { MODE_SCOPED_SETTINGS } from '../mode-scoped-settings.generated.js';
import {
  SLICE_MODES,
  sessionDepthSourceWidths,
  splitDraftIntoModes,
  unscopedDraftRows,
  validateModeSliceFields
} from './mode-scoped-migration.js';

const { nextTick } = window.Vue;

// How LOSAT runs: one app-level setting for both drawings (Session 46
// `ui.losatExecution`; a Session 44 or older draft keeps it in `config.losat`).
const LOSAT_EXECUTION_FIELDS = MODE_SCOPED_SETTINGS.losatExecutionFields;
/**
 * Installs saved LOSAT execution settings; a field the value lacks takes its
 * normalized default, as a saved draft always read.
 * @param {Record<string, any>} source
 */
export const applyLosatExecutionData = (source) => {
  const target = state.losatExecution;
  const rawParallelWorkers = String(source.parallelWorkers ?? '').trim().toLowerCase();
  const parsedParallelWorkers = Number(rawParallelWorkers);
  target.parallelWorkers = Number.isInteger(parsedParallelWorkers) && parsedParallelWorkers >= 1
    ? rawParallelWorkers
    : undefined;
  const rawExecutionMode = String(source.executionMode ?? '').trim().toLowerCase();
  target.executionMode = ['auto', 'serial', 'threaded'].includes(rawExecutionMode)
    ? rawExecutionMode
    : 'auto';
  const rawThreadsPerJob = String(source.threadsPerJob ?? 'auto').trim().toLowerCase();
  const parsedThreadsPerJob = Number(rawThreadsPerJob);
  target.threadsPerJob = rawThreadsPerJob === 'auto' ||
    (Number.isInteger(parsedThreadsPerJob) && parsedThreadsPerJob >= 1)
    ? rawThreadsPerJob
    : 'auto';
  const rawTotalThreadBudget = String(source.totalThreadBudget ?? 'safe').trim().toLowerCase();
  const parsedTotalThreadBudget = Number(rawTotalThreadBudget);
  target.totalThreadBudget = ['safe', 'auto', 'available'].includes(rawTotalThreadBudget) ||
    (Number.isInteger(parsedTotalThreadBudget) && parsedTotalThreadBudget >= 1)
    ? (rawTotalThreadBudget === 'auto' ? 'safe' : rawTotalThreadBudget)
    : 'safe';
};
/** @returns {Record<string, any>} The settings as Save writes them (`ui.losatExecution`). */
export const buildLosatExecutionData = () => cloneJsonData(Object.fromEntries(
  LOSAT_EXECUTION_FIELDS.map((/** @type {string} */ field) => [field, state.losatExecution[field]])
));

/** @import { ActiveWebConfig } from './session-active-config-contract.js' */
/** @import { CanonicalRenderEnvelope, CanonicalRenderRequest } from './session-request.js' */
/** @import { FeatureCatalog } from './feature-catalog.js' */
/** @import { SessionResourceSource } from './session-resources.js' */
/** @import { ArtifactSlot } from './artifact-slot.js' */
/** @import { FeatureOverrideDraft } from './feature-placement.js' */
/** @import { DrawingState } from '../state.js' */

/**
 * The editor state of a Session's committed Result: the Legend's original
 * order and colors, the original stroke, the reset receipt, and the admitted
 * feature catalog. The Legend and stroke edits are each mode's
 * (`SessionModeSlice.editorState`, R7).
 * @typedef {object} SessionEditorState
 * @property {Record<string, any>} legend
 * @property {Record<string, any>} originalSvgStroke
 * @property {Record<string, any> | null} alignmentResetReceipt
 * @property {FeatureCatalog | null} featureCatalog
 */

/**
 * The feature edits of a Session. `featureOverrides` carries the identity-keyed
 * rows that JavaScript builds itself; the other fields are open.
 * @typedef {object} SessionFeatureState
 * @property {number} selectedFeatureRecordIdx
 * @property {Record<string, any>} featureColorOverrides
 * @property {Record<string, any>[]} featureVisibilityManualRules
 * @property {FeatureOverrideDraft} featureOverrides
 * @property {Record<string, any>[]} labelOverrideRows
 * @property {Record<string, any>} labelTextBulkOverrides
 */

/**
 * The other diagram mode's Result set (E1): written only when both
 * modes keep a Result. Its fields mirror the top-level fields of one committed
 * set; its request's resources are in the top-level `resources` table.
 * @typedef {object} SessionOtherModeResult
 * @property {CanonicalRenderRequest} renderRequest Its mode differs from the top-level request's.
 * @property {Record<string, any>[]} results At least one.
 * @property {{ featureCatalog: FeatureCatalog, alignmentResetReceipt: Record<string, any> | null,
 *   legend: { originalOrder: string[], originalColors: Record<string, string> },
 *   originalSvgStroke: { color: string | null, width: number | null } }} editorState
 * @property {Record<string, any>} ui The selected Result and the generated layout and palette.
 * @property {Record<string, any>} runMetadata
 * @property {Record<string, any>} [cliInvocation]
 */

/**
 * A gbdraw Session in the current writer format (`SESSION_VERSION`; the
 * writer format only, readers take unvalidated data, R14). The render fields
 * stay in `CanonicalRenderRequest` and `ActiveWebConfig`; the top-level fields
 * equal `CURRENT_SESSION_TOP_LEVEL_FIELDS` in `gbdraw/session_io.py` (a pytest
 * checks this list). `title`, `runMetadata`, `legacyArtifacts`,
 * `otherModeResult`, and `cliInvocation` are written only when they exist.
 * @typedef {object} GbdrawSession
 * @property {string} format `gbdraw-session`.
 * @property {number} version
 * @property {string} createdAt
 * @property {string} [title]
 * @property {CanonicalRenderRequest} renderRequest
 * @property {CanonicalRenderEnvelope['resources']} resources
 * @property {CanonicalRenderEnvelope['webFiles']} webFiles
 * @property {Record<string, any>} ui App-level UI: the selected Result, the shown layout, and `losatExecution`.
 * @property {Record<string, any>[]} results
 * @property {SessionEditorState} editorState The committed Result's artifacts (catalog, receipt, original colors).
 * @property {Record<'circular' | 'linear', SessionModeSlice>} modes Each diagram mode's drawing (PR-1, Session 46).
 * @property {Record<string, any>} orthogroupState
 * @property {Record<string, any>} losatCache
 * @property {Record<string, any>} losatDerivedCache
 * @property {Record<string, any>} proteinIdentityManifest
 * @property {Record<string, any>} [legacyArtifacts]
 * @property {Record<string, any>} [runMetadata]
 * @property {SessionOtherModeResult} [otherModeResult]
 * @property {Record<string, any>} [cliInvocation]
 * @property {Record<string, any>} [cliOptions] CLI options kept for both modes (written by the CLI).
 */

/**
 * One diagram mode's drawing in a Session 46 (`modes.<mode>`): the mode's
 * settings, edits, and per-mode UI.
 * @typedef {object} SessionModeSlice
 * @property {ActiveWebConfig} config
 * @property {Omit<SessionFeatureState, 'selectedFeatureRecordIdx'>} features
 * @property {{ legend: Record<string, any>, featureStrokes: Record<string, any> }} editorState
 * @property {Record<string, any>} ui
 */

/**
 * What an older Session's Result needs from the Session: the committed legend
 * and title sides with the saved user offsets, and the extracted features with
 * the two saved stroke override maps.
 * @typedef {object} LegacyResultSvgData
 * @property {{ legendSide: string, titleSide: string, userDeltas: Record<string, number[] | null> }} composition
 * @property {{
 *   features: Record<string, any>[],
 *   legendStrokeOverrides: Record<string, any>,
 *   featureStrokeOverrides: Record<string, any>
 * }} strokes
 */

/**
 * The composition root's transform of an older Session's Result (R13 port):
 * it gives a Result without composition metadata the legacy composition and
 * projects the saved strokes. Returns whether the SVG changed.
 * @typedef {(svg: Element, data: LegacyResultSvgData) => boolean} LegacyResultSvgTransform
 */
export const SESSION_VERSION = 46;
const CURRENT_AUTHORITY_SESSION_MIN_VERSION = 40;
const LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION = 32;
const SUPPORTED_SESSION_VERSIONS = new Set([
  27, 28, 29, 30, 31, 32, 33, 39, 40, 41, 42, 44, SESSION_VERSION
]);
const CURRENT_ARTIFACT_SESSION_MIN_VERSION = 39;
const LOSAT_DERIVED_CACHE_LIMIT = 16;
// D-25 (PD-OI-079): a Result without current feature metadata (a legacy
// Session) is saved only after one Generate; the error offers that Generate.
// E1: `diagramMode` names the mode whose Result it is when that mode is not
// shown; the error's Generate action then runs in that mode.
/** @param {Record<string, any>} [context] */
const sessionSaveRequiresGenerate = (context = {}) => diagnosticError(
  'SESSION_SAVE_REQUIRES_GENERATE', context, /** @type {Record<string, any>} */ ({
    operation: 'session-save', stage: 'result-admission'
  })
);
// A Save wrapper keeps a recognized cause; otherwise it reports the bounded fallback.
const recognizedCauseOr = (error, fallback) => {
  const model = normalizeUserFacingError(error);
  return model && !['UNKNOWN', 'VALIDATION_UNCLASSIFIED'].includes(model.code)
    ? diagnosticError(model.code, model.context, { stage: model.stage }) : fallback;
};
const isPlainObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

const cloneColors = (colors) => ({ ...(colors || {}) });

const hasColorEntries = (colors) =>
  Boolean(colors) && typeof colors === 'object' && !Array.isArray(colors) && Object.keys(colors).length > 0;

const normalizeColorMap = (colors) => {
  const normalized = {};
  if (!colors || typeof colors !== 'object' || Array.isArray(colors)) return normalized;
  Object.entries(colors).forEach(([key, value]) => {
    normalized[key] = resolveColorToHex(String(value || '').trim());
  });
  return normalized;
};

const paletteColorsFromDefinitions = (paletteName) => {
  const name = String(paletteName || '').trim();
  if (!name) return null;
  const definitions = state.paletteDefinitions?.value || {};
  const colors = definitions[name];
  return hasColorEntries(colors) ? normalizePaletteColors(cloneColors(colors)) : null;
};

const cloneStringMap = (source) => {
  const cloned = {};
  if (!source || typeof source !== 'object' || Array.isArray(source)) return cloned;
  Object.entries(source).forEach(([key, value]) => {
    const normalizedKey = String(key || '').trim();
    if (!normalizedKey) return;
    cloned[normalizedKey] = String(value ?? '');
  });
  return cloned;
};

const normalizeFeatureVisibilityRulesForSession = (rules) => (
  Array.isArray(rules) ? rules.map((rule) => normalizeFeatureVisibilityRule(rule)) : []
);

// A Session before manual rules had their own field stored the editor's
// per-feature rows among the rules (`featureVisibilityRules`); its rows keyed
// by rendered ID join the other rendered-ID edits for migration.
const splitLegacyFeatureVisibilityRules = (features = {}) => {
  if (Array.isArray(features.featureVisibilityManualRules) || !Array.isArray(features.featureVisibilityRules)) {
    return features;
  }
  const { manualRules, overrides } = splitLegacyVisibilityRules(features.featureVisibilityRules);
  return {
    ...features,
    featureVisibilityManualRules: manualRules,
    featureVisibilityOverrides: { ...overrides, ...(features.featureVisibilityOverrides || {}) }
  };
};

// The identity-keyed per-feature edit draft; a Session or History value is
// checked as the request rows are (services/feature-placement.js).
const featureOverridesForState = (value) => featureDraftMap(
  canonicalFeatureOverrides(isPlainObject(value) ? value : {})
);

const sanitizeExtractedFeatureForSession = (feature) => {
  if (!feature || typeof feature !== 'object' || Array.isArray(feature)) return feature;
  const {
    nucleotide_sequence: _nucleotideSequence,
    amino_acid_sequence: _aminoAcidSequence,
    nucleotideSequence: _nucleotideSequenceAlias,
    aminoAcidSequence: _aminoAcidSequenceAlias,
    ...rest
  } = feature;
  return rest;
};

const sanitizeExtractedFeaturesForSession = (features) => {
  if (!Array.isArray(features)) return [];
  return features.map((feature) => sanitizeExtractedFeatureForSession(feature));
};

const replaceStringMap = (target, source) => {
  Object.keys(target).forEach((key) => delete target[key]);
  Object.entries(cloneStringMap(source)).forEach(([key, value]) => {
    target[key] = value;
  });
};

const cloneQualifierPriorityRules = (rules) => {
  const cloned = [];
  if (!Array.isArray(rules)) return cloned;

  rules.forEach((rule) => {
    if (!rule || typeof rule !== 'object' || Array.isArray(rule)) return;
    const feat = String(rule.feat ?? '').trim();
    const order = String(rule.order ?? '').trim();
    if (!feat || !order) return;

    const existingIndex = cloned.findIndex((entry) => entry.feat === feat);
    if (existingIndex >= 0) {
      cloned[existingIndex].order = order;
    } else {
      cloned.push({ feat, order });
    }
  });

  return cloned;
};

/** @param {DrawingState} drawing */
const replaceQualifierPriorityRules = (drawing, rules) => {
  drawing.manualPriorityRules.splice(
    0,
    drawing.manualPriorityRules.length,
    ...cloneQualifierPriorityRules(rules)
  );
};

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

    // 5. Null restores an explicit value to Auto, as well as accepting one from Auto.
    if (typeof targetValue === typeof sourceValue || targetValue === null || sourceValue === null) {
      target[key] = sourceValue;
    }
  });
};

const parseMultiRecordPositionToken = (value) => {
  const raw = String(value ?? '').trim();
  const separatorIndex = raw.lastIndexOf('@');
  if (separatorIndex <= 0 || separatorIndex === raw.length - 1) return null;

  const selector = raw.slice(0, separatorIndex).trim();
  const row = Number(raw.slice(separatorIndex + 1).trim());
  if (!selector || !Number.isInteger(row) || row <= 0) return null;
  return { selector, row };
};

const multiRecordPositionsFromCliInvocation = (cliInvocation) => {
  const args = Array.isArray(cliInvocation?.args) ? cliInvocation.args : [];
  const positions = [];
  const seenSelectors = new Set();

  args.forEach((arg, index) => {
    if (arg !== '--multi_record_position') return;
    const position = parseMultiRecordPositionToken(args[index + 1]);
    if (!position || seenSelectors.has(position.selector)) return;
    seenSelectors.add(position.selector);
    positions.push(position);
  });

  return positions;
};

const hydrateMissingMultiRecordPositionsFromCliInvocation = (config, cliInvocation) => {
  if (!isPlainObject(config) || !isPlainObject(config.form) || config.form.multi_record_canvas !== true) {
    return;
  }

  const adv = isPlainObject(config.adv) ? config.adv : {};
  if (Array.isArray(adv.multi_record_positions)) return;

  const positions = multiRecordPositionsFromCliInvocation(cliInvocation);
  if (positions.length === 0) return;

  adv.multi_record_positions = positions;
  config.adv = adv;
};

/**
 * @param {Record<string, any> | null | undefined} invocation
 * @returns {boolean}
 */
const isCliInvocationSessionExportable = (invocation) => {
  if (!invocation || typeof invocation !== 'object') return false;
  if (invocation.sessionExportable === false) return false;
  const bindings = Array.isArray(invocation.fileBindings) ? invocation.fileBindings : [];
  return bindings.every((binding) => String(binding?.slot || '').startsWith('files.'));
};

const makeSafeFilename = (name) => {
  const cleaned = String(name || '')
    .replace(/[^\w.-]+/g, '_')
    .replace(/^_+|_+$/g, '');
  return cleaned || 'gbdraw_session';
};

const buildSessionFilename = (title) => {
  const base = String(title || '').trim();
  if (!base) return 'gbdraw_session.json.gz';
  const safe = makeSafeFilename(base);
  return `${safe}.gbdraw-session.json.gz`;
};

const normalizeLegendPosition = (value, fallback = 'left') => {
  const normalized = String(value || '').trim().toLowerCase();
  return normalized || fallback;
};

const normalizeLabelRendering = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['auto', 'embedded_only', 'external_only'].includes(normalized) ? normalized : 'auto';
};

const withCurrentLinearLabelVisibility = (configData) => {
  if (!isPlainObject(configData) || !isPlainObject(configData.adv)) return configData;
  return {
    ...configData,
    adv: migrateLegacyLinearLabelVisibility(configData.adv)
  };
};

const normalizePositiveNumberOrNull = (value) => {
  if (
    value === null ||
    value === undefined ||
    value === '' ||
    String(value).trim().toLowerCase() === 'auto'
  ) {
    return null;
  }
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric > 0 ? numeric : null;
};

/**
 * @param {Record<string, any>} [configData] The stored active configuration.
 * @returns {Record<string, any>}
 */
const migrateImportedCircularTrackSlots = (configData = {}) => {
  const adv = configData && typeof configData === 'object' ? configData.adv : null;
  if (!adv || typeof adv !== 'object' || Array.isArray(adv)) return configData;
  if (!Object.prototype.hasOwnProperty.call(adv, 'circular_track_slots')) return configData;
  if (!Array.isArray(adv.circular_track_slots)) {
    throw new Error('Custom Track Slots must be an array.');
  }
  const sourceSchemaVersion = adv.circular_track_slots_schema_version;
  if (sourceSchemaVersion === CIRCULAR_TRACK_SLOT_SCHEMA_VERSION) return configData;
  if (sourceSchemaVersion !== LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION) {
    throw new Error(
      `Custom Track Slots use an obsolete schema. Recreate the slots with schema version ${CIRCULAR_TRACK_SLOT_SCHEMA_VERSION}.`
    );
  }

  const defaultNt = adv.nt || 'GC';
  const preset = configData.form?.track_type || 'tuckin';
  return {
    ...configData,
    adv: {
      ...adv,
      circular_track_slots_schema_version: CIRCULAR_TRACK_SLOT_SCHEMA_VERSION,
      circular_track_slots: adv.circular_track_slots.map((slot, index) => (
        typeof slot === 'string'
          ? parseCircularTrackSlotSpec(
              migrateLegacyCircularTrackSlotSpec(slot),
              index,
              defaultNt,
              preset
            )
          : migrateLegacyCircularTrackSlot(slot)
      ))
    }
  };
};

const migratePersistedWebOptionValues = (configData = {}) => {
  if (!configData || typeof configData !== 'object' || Array.isArray(configData)) {
    return configData;
  }
  const migratedNames = migratePersistedWebStateFieldNames(configData);
  const form = isPlainObject(migratedNames.form) ? { ...migratedNames.form } : migratedNames.form;
  const adv = isPlainObject(migratedNames.adv) ? { ...migratedNames.adv } : migratedNames.adv;
  if (isPlainObject(form) && Object.prototype.hasOwnProperty.call(form, 'linear_track_layout')) {
    form.linear_track_layout = migratePersistedLinearTrackLayout(form.linear_track_layout);
  }
  if (isPlainObject(adv) && Object.prototype.hasOwnProperty.call(adv, 'label_placement')) {
    adv.label_placement = migratePersistedLinearLabelPlacement(adv.label_placement);
  }
  if (isPlainObject(adv) && Object.prototype.hasOwnProperty.call(adv, 'multi_record_size_mode')) {
    adv.multi_record_size_mode = migratePersistedCircularMultiRecordSizeMode(
      adv.multi_record_size_mode
    );
  }
  return withCurrentLinearLabelVisibility({
    ...migratedNames,
    ...(form === undefined ? {} : { form }),
    ...(adv === undefined ? {} : { adv })
  });
};

/**
 * @param {Record<string, any>} [configData] The stored active configuration.
 * @param {number | null} [sourceSessionVersion]
 * @returns {Record<string, any>}
 */
const migrateImportedLinearTrackSlots = (configData = {}, sourceSessionVersion = null) => {
  const adv = configData && typeof configData === 'object' ? configData.adv : null;
  if (!adv || typeof adv !== 'object' || Array.isArray(adv)) return configData;
  if (!Object.prototype.hasOwnProperty.call(adv, 'linear_track_slots')) return configData;
  if (!Array.isArray(adv.linear_track_slots)) {
    throw new Error('Custom Track Slots must be an array.');
  }

  const hasStoredSchemaVersion = Object.prototype.hasOwnProperty.call(
    adv,
    'linear_track_slots_schema_version'
  );
  const storedSchemaVersion = adv.linear_track_slots_schema_version;
  if (
    hasStoredSchemaVersion &&
    (
      !Number.isInteger(storedSchemaVersion) ||
      ![LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION, LINEAR_TRACK_SLOT_SCHEMA_VERSION].includes(storedSchemaVersion)
    )
  ) {
    throw new Error(
      `Custom Track Slots use an obsolete schema. Recreate the slots with schema version ${LINEAR_TRACK_SLOT_SCHEMA_VERSION}.`
    );
  }
  const sessionUsesLegacySemantics = (
    Number.isInteger(sourceSessionVersion) &&
    // Number.isInteger above is false for null, so the version is a number here.
    /** @type {number} */ (sourceSessionVersion) <= LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION
  );
  const sourceSchemaVersion = sessionUsesLegacySemantics
    ? LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION
    : storedSchemaVersion;
  if (![LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION, LINEAR_TRACK_SLOT_SCHEMA_VERSION].includes(sourceSchemaVersion)) {
    throw new Error(
      `Custom Track Slots use an obsolete schema. Recreate the slots with schema version ${LINEAR_TRACK_SLOT_SCHEMA_VERSION}.`
    );
  }

  return {
    ...configData,
    adv: {
      ...adv,
      linear_track_slots_schema_version: LINEAR_TRACK_SLOT_SCHEMA_VERSION,
      linear_track_slots: migrateLinearTrackSlotsToCurrentSchema(
        adv.linear_track_slots,
        sourceSchemaVersion
      )
    }
  };
};

const hasStoredLayoutValue = (value) => typeof value === 'string' && value.trim() !== '';

const normalizePositiveInteger = (value, fallback) => {
  const numeric = Number(value);
  return Number.isInteger(numeric) && numeric > 0 ? numeric : fallback;
};

const normalizeBlastpMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['pairwise', 'orthogroup', 'collinear'].includes(normalized) ? normalized : 'orthogroup';
};

const normalizeCollinearColorMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
  if (normalized === 'identity') return 'average_identity';
  return ['average_identity', 'orientation', 'orientation_identity'].includes(normalized) ? normalized : 'orientation';
};

const normalizeCircularConservationSource = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return normalized === 'upload' ? 'upload' : 'losat';
};

const normalizeCircularConservationLosatProgram = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return normalized === 'tblastx' ? 'tblastx' : 'blastn';
};

const normalizeCircularConservationReference = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['auto', 'query', 'subject'].includes(normalized) ? normalized : 'auto';
};

const normalizeHexColor = (value, fallback = '#4e79a7') => {
  const resolved = resolveColorToHex(String(value || fallback).trim());
  const color = String(resolved || fallback).trim();
  const shortMatch = color.match(/^#([0-9a-fA-F]{3})$/);
  if (shortMatch) {
    return `#${shortMatch[1].split('').map((char) => char + char).join('').toLowerCase()}`;
  }
  return /^#[0-9a-fA-F]{6}$/.test(color) ? color.toLowerCase() : fallback;
};

const normalizeOptionalHexColor = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const resolved = resolveColorToHex(String(value).trim());
  const color = String(resolved || '').trim();
  const shortMatch = color.match(/^#([0-9a-fA-F]{3})$/);
  if (shortMatch) {
    return `#${shortMatch[1].split('').map((char) => char + char).join('').toLowerCase()}`;
  }
  return /^#[0-9a-fA-F]{6}$/.test(color) ? color.toLowerCase() : null;
};

const normalizeSessionLegendColor = (value) => {
  const color = String(value || '').trim();
  if (!color) return null;
  if (/^#(?:[0-9a-f]{3}|[0-9a-f]{4}|[0-9a-f]{6}|[0-9a-f]{8})$/i.test(color)) {
    return color.toLowerCase();
  }
  if (/^[a-z]+$/i.test(color)) return color.toLowerCase();
  if (/^rgba?\(\s*[-+.\d%]+(?:\s*[,/]\s*|\s+)[-+.\d%]+(?:\s*[,/]\s*|\s+)[-+.\d%]+(?:\s*[,/]\s*[-+.\d%]+)?\s*\)$/i.test(color)) {
    return color;
  }
  if (/^hsla?\(\s*[-+.\d]+(?:deg|grad|rad|turn)?(?:\s*[,/]\s*|\s+)[-+.\d%]+(?:\s*[,/]\s*|\s+)[-+.\d%]+(?:\s*[,/]\s*[-+.\d%]+)?\s*\)$/i.test(color)) {
    return color;
  }
  return null;
};

const normalizeSessionLegendEntries = (entries) => {
  if (!Array.isArray(entries)) return [];
  const normalized = [];
  const captions = new Set();
  entries.forEach((entry) => {
    if (!isPlainObject(entry)) return;
    const caption = String(entry.caption || '').trim();
    const color = normalizeSessionLegendColor(entry.color);
    if (!caption || !color || captions.has(caption)) return;
    captions.add(caption);
    normalized.push({
      ...cloneJsonData(entry),
      caption,
      color,
      showStroke: Boolean(entry.showStroke),
      featureIds: Array.isArray(entry.featureIds)
        ? entry.featureIds.map((id) => String(id || '').trim()).filter(Boolean)
        : []
    });
  });
  return normalized;
};

const normalizeStrokeWidth = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric >= 0 ? numeric : null;
};

const cloneJsonArray = (value) => {
  if (!Array.isArray(value)) return [];
  return cloneJsonValue(value, []);
};

const cloneJsonObject = (value) => {
  if (!isPlainObject(value)) return {};
  return cloneJsonValue(value, {});
};

const normalizeLegendColorOverrides = (source) => {
  const normalized = {};
  if (!isPlainObject(source)) return normalized;
  Object.entries(source).forEach(([key, value]) => {
    const caption = String(key || '').trim();
    const color = normalizeOptionalHexColor(value);
    if (!caption || !color) return;
    normalized[caption] = color;
  });
  return normalized;
};

const normalizeStrokeOverride = (source, { requireOverride = false } = {}) => {
  if (!isPlainObject(source)) return null;
  const normalized = {};
  const strokeColor = normalizeOptionalHexColor(source.strokeColor);
  const strokeWidth = normalizeStrokeWidth(source.strokeWidth);
  const originalStrokeColor = normalizeOptionalHexColor(source.originalStrokeColor);
  const originalStrokeWidth = normalizeStrokeWidth(source.originalStrokeWidth);

  if (strokeColor !== null) normalized.strokeColor = strokeColor;
  if (strokeWidth !== null) normalized.strokeWidth = strokeWidth;
  if (Object.prototype.hasOwnProperty.call(source, 'originalStrokeColor')) {
    normalized.originalStrokeColor = originalStrokeColor;
  }
  if (Object.prototype.hasOwnProperty.call(source, 'originalStrokeWidth')) {
    normalized.originalStrokeWidth = originalStrokeWidth;
  }

  if (
    requireOverride &&
    !Object.prototype.hasOwnProperty.call(normalized, 'strokeColor') &&
    !Object.prototype.hasOwnProperty.call(normalized, 'strokeWidth')
  ) {
    return null;
  }
  return Object.keys(normalized).length > 0 ? normalized : null;
};

const normalizeStrokeOverrideMap = (source, { requireOverride = false } = {}) => {
  const normalized = {};
  if (!isPlainObject(source)) return normalized;
  Object.entries(source).forEach(([key, value]) => {
    const normalizedKey = String(key || '').trim();
    if (!normalizedKey) return;
    const override = normalizeStrokeOverride(value, { requireOverride });
    if (!override) return;
    normalized[normalizedKey] = override;
  });
  return normalized;
};

const normalizeStringArray = (source) => {
  if (!Array.isArray(source)) return [];
  return source
    .map((value) => String(value || '').trim())
    .filter(Boolean);
};

const normalizeCircularConservationSeries = (series) => {
  if (!Array.isArray(series)) return [];
  return series
    .filter((entry) => entry && typeof entry === 'object')
    .map((entry, index) => ({
      sourceKey: String(entry.sourceKey || ''),
      fileName: String(entry.fileName || ''),
      sourceIndex: Number.isInteger(Number(entry.sourceIndex)) ? Number(entry.sourceIndex) : index,
      label: String(entry.label ?? entry.name ?? ''),
      color: normalizeHexColor(entry.color, '#4e79a7'),
      losat_gencode: normalizePositiveInteger(entry.losat_gencode, 1)
    }));
};

const DEPTH_TRACK_FALLBACK_COLORS = [
  '#4A90E2',
  '#E45756',
  '#2CA02C',
  '#F28E2B',
  '#9467BD',
  '#8C564B',
  '#17BECF',
  '#7F7F7F'
];

const normalizeDepthTrackConfig = (entry, index, legacyAdv = {}) => {
  const source = entry && typeof entry === 'object' && !Array.isArray(entry) ? entry : {};
  const hasHeight = Object.prototype.hasOwnProperty.call(source, 'height');
  const fallbackColor =
    index === 0
      ? String(legacyAdv.depth_color || DEPTH_TRACK_FALLBACK_COLORS[0])
      : DEPTH_TRACK_FALLBACK_COLORS[index % DEPTH_TRACK_FALLBACK_COLORS.length];
  return {
    label: String(source.label ?? (index === 0 ? 'Depth' : `Depth ${index + 1}`)),
    color: resolveColorToHex(String(source.color || fallbackColor)),
    height: normalizePositiveNumberOrNull(hasHeight ? source.height : legacyAdv.depth_height),
    large_tick_interval: normalizePositiveNumberOrNull(
      source.large_tick_interval ?? (index === 0 ? legacyAdv.depth_large_tick_interval : null)
    ),
    small_tick_interval: normalizePositiveNumberOrNull(
      source.small_tick_interval ?? (index === 0 ? legacyAdv.depth_small_tick_interval : null)
    ),
    tick_font_size: normalizePositiveNumberOrNull(
      source.tick_font_size ?? (index === 0 ? legacyAdv.depth_tick_font_size : null)
    )
  };
};

const normalizeDepthTracks = (tracks, legacyAdv = {}) => {
  const rawTracks = Array.isArray(tracks) ? tracks : [];
  const normalized = rawTracks.map((entry, index) => normalizeDepthTrackConfig(entry, index, legacyAdv));
  if (normalized.length === 0) {
    normalized.push(normalizeDepthTrackConfig(null, 0, legacyAdv));
  }
  return normalized;
};

/** @type {string | null} */
let lastSessionFilename = null;
/** @type {Record<string, any> | null} */
let preservedCliOptions = null;
/** @type {Record<string, any> | null} */
let committedCanonicalSession = null;
/** @type {ReturnType<typeof adoptCurrentSessionResources> | null} */
let activeSessionResourceTable = null;
/** @type {Record<string, any> | null} */
let adoptedProteinIdentityManifest = null;
// Current-session preflight receipts pair each adopted cache value with the
// exact manifest that already validated its protein references and raw text.
const adoptedLosatCacheValues = new WeakMap();

const rawReactiveValue = (value) => (
  globalThis.window?.Vue?.toRaw?.(value) ?? value
);

const cloneCanonicalSession = (canonical) => {
  if (
    !canonical
    || !isPlainObject(canonical.renderRequest)
    || !isPlainObject(canonical.resources)
  ) return null;
  return {
    renderRequest: cloneJsonData(canonical.renderRequest),
    resources: cloneJsonData(canonical.resources),
    webFiles: isPlainObject(canonical.webFiles)
      ? cloneJsonData(canonical.webFiles)
      : {}
  };
};

const normalizedArrowGeometryState = (adv = {}) => ({
  arrow_head_length_ratio: arrowHeadLengthRatioForState(
    adv?.arrow_head_length_ratio
  ),
  arrow_shaft_width_ratio: normalizeArrowShaftWidthRatio(
    adv?.arrow_shaft_width_ratio
  )
});

const persistedArrowGeometryValue = (value) => {
  if (typeof value !== 'string' || value.trim() === '') return value;
  const numeric = Number(value.trim());
  return Number.isFinite(numeric) ? numeric : value;
};

const normalizedPersistedArrowGeometryState = (adv = {}) =>
  normalizedArrowGeometryState({
    arrow_head_length_ratio: persistedArrowGeometryValue(
      adv?.arrow_head_length_ratio
    ),
    arrow_shaft_width_ratio: persistedArrowGeometryValue(
      adv?.arrow_shaft_width_ratio
    )
  });

const replaceLinearComparisonPlan = (target, source) => {
  const normalized = normalizeLinearComparisonPlan(source);
  target.mode = normalized.mode;
  target.defaultSource = normalized.defaultSource;
  if (!Array.isArray(target.edges)) target.edges = [];
  target.edges.splice(0, target.edges.length, ...normalized.edges);
};

const serializeLinearComparisonPlan = (plan) => {
  const normalized = normalizeLinearComparisonPlan(plan);
  return {
    mode: normalized.mode,
    defaultSource: normalized.defaultSource,
    edges: normalized.edges.map((edge) => ({
      id: edge.id,
      queryUid: edge.queryUid,
      subjectUid: edge.subjectUid,
      included: edge.included,
      fileActive: edge.fileActive,
      losatFilenameActive: edge.losatFilenameActive,
      source: edge.source,
      losatFilename: edge.losatFilename
    }))
  };
};

// A drawing's settings as its Session 46 slice `config` holds them: the
// registry rows of `config` (mode-scoped-settings.generated.js) and nothing
// else. App-level settings (LOSAT execution, Instant Preview, the rich popup)
// and the Session's CLI provenance are saved outside the slices.
/** @param {DrawingState} drawing */
export const buildConfigData = (drawing) => ({
  form: drawing.form,
  adv: {
    ...drawing.adv,
    feature_shapes: {
      repeat_region: defaultFeatureRendering('repeat_region'),
      ...normalizeFeatureRenderingMap(drawing.adv.feature_shapes || {})
    }
  },
  losat: cloneJsonData(drawing.losat || {}),
  colors: drawing.currentColors.value,
  palette: drawing.selectedPalette.value,
  rules: drawing.manualSpecificRules,
  qualifierPriorityRules: cloneQualifierPriorityRules(drawing.manualPriorityRules),
  filterMode: drawing.filterMode.value,
  whitelist: drawing.manualWhitelist,
  blacklistText: drawing.manualBlacklist.value,
  losatProgram: drawing.losatProgram.value,
  circularConservation: drawing.circularConservation,
  annotationSets: normalizeAnnotationSets(drawing.annotationSets),
  recordDisplayDrafts: cloneJsonData(drawing.recordDisplayDrafts),
  featurePlacementOverrides: cloneJsonData(drawing.featurePlacementOverrides),
  linearRecordLayout: {
    enabled: Boolean(drawing.linearRecordLayoutEnabled.value),
    recordGap: Number(drawing.linearRecordGap.value) || 0,
    rows: (drawing.linearRecordRows || []).map((entry) => ({
      uid: String(entry?.uid || ''),
      row: Number(entry?.row) || 1,
      ...(entry?.canonicalCardinality === 'exactly_one'
        ? { canonicalCardinality: 'exactly_one' } : {}),
      ...(Number(entry?.canonicalRow) === Number(entry?.row)
        && Number.isInteger(entry?.canonicalColumn) && entry.canonicalColumn > 0
        ? { canonicalRow: entry.canonicalRow, canonicalColumn: entry.canonicalColumn }
        : {})
    }))
  },
  linearComparisonPlan: serializeLinearComparisonPlan(drawing.linearComparisonPlan),
  importedComparisonResolution: serializeImportedComparisonResolution(
    drawing.importedComparisonIntent
  ),
  unmanagedConfigOverrides: cloneJsonData(drawing.unmanagedConfigOverrides || {}),
  webEdits: {
    orthogroupNameOverrides: cloneStringMap(drawing.orthogroupNameOverrides),
    orthogroupDescriptionOverrides: cloneStringMap(drawing.orthogroupDescriptionOverrides),
    orthogroupDormantOverrides: normalizeOrthogroupDormantOverrides(drawing.orthogroupDormantOverrides)
  }
});

const defaultEditorStateData = () => ({
  legend: {
    entries: [],
    deletedEntries: [],
    dormantEntries: [],
    originalOrder: [],
    originalColors: {},
    colorOverrides: {},
    strokeOverrides: {},
    addedCaptions: []
  },
  featureStrokes: {
    overrides: {}
  },
  originalSvgStroke: {
    color: null,
    width: null
  },
  featureCatalog: null
});

// State, History, and Session rollback hold the Generate-owned catalog by
// reference; only an admitted catalog may enter state. Vue proxies obscure the
// identity that admission records.
const admittedFeatureCatalog = (catalog) => {
  const rawCatalog = rawReactiveValue(catalog);
  return isAdoptedFeatureCatalog(rawCatalog) ? rawCatalog : null;
};

/** @param {DrawingState} drawing */
export const buildEditorStateData = (drawing) => ({
  legend: {
    entries: cloneJsonArray(drawing.legendEntries.value),
    deletedEntries: cloneJsonArray(drawing.deletedLegendEntries.value),
    dormantEntries: cloneJsonArray(drawing.dormantLegendEntries.value),
    originalOrder: cloneJsonArray(state.originalLegendOrder.value),
    originalColors: cloneStringMap(state.originalLegendColors.value),
    colorOverrides: cloneJsonObject(drawing.legendColorOverrides),
    strokeOverrides: cloneJsonObject(drawing.legendStrokeOverrides),
    addedCaptions: Array.from(drawing.addedLegendCaptions.value || [])
      .map((caption) => String(caption || '').trim())
      .filter(Boolean)
  },
  featureStrokes: {
    overrides: cloneJsonObject(drawing.featureStrokeOverrides)
  },
  originalSvgStroke: {
    color: state.originalSvgStroke.value?.color ?? null,
    width: state.originalSvgStroke.value?.width ?? null
  },
  alignmentResetReceipt: cloneJsonValue(state.similarityAlignmentResetReceipt?.value, null),
  featureCatalog: admittedFeatureCatalog(state.featureCatalog?.value)
});

/**
 * @param {Record<string, any>} [editorState]
 * @param {{ featureCatalog?: unknown }} [options]
 */
const normalizeEditorStateData = (editorState = {}, { featureCatalog = undefined } = {}) => {
  const defaults = defaultEditorStateData();
  const source = isPlainObject(editorState) ? editorState : {};
  const legend = isPlainObject(source.legend) ? source.legend : {};
  const featureStrokes = isPlainObject(source.featureStrokes) ? source.featureStrokes : {};
  const originalSvgStroke = isPlainObject(source.originalSvgStroke) ? source.originalSvgStroke : {};
  // A Session 46 slice lists the renamed rows its Result does not draw
  // (OV-120) after the shown rows, marked `dormant`.
  const entries = normalizeSessionLegendEntries(legend.entries);
  const dormantEntries = [
    ...normalizeSessionLegendEntries(legend.dormantEntries),
    ...entries.filter((entry) => entry.dormant === true)
  ].map(({ dormant: _dormant, ...entry }) => entry);

  return {
    legend: {
      entries: entries.filter((entry) => entry.dormant !== true),
      deletedEntries: normalizeSessionLegendEntries(legend.deletedEntries),
      dormantEntries,
      originalOrder: normalizeStringArray(legend.originalOrder),
      originalColors: normalizeLegendColorOverrides(legend.originalColors),
      colorOverrides: normalizeLegendColorOverrides(legend.colorOverrides),
      strokeOverrides: normalizeStrokeOverrideMap(legend.strokeOverrides, { requireOverride: true }),
      addedCaptions: normalizeStringArray(legend.addedCaptions)
    },
    featureStrokes: {
      overrides: normalizeStrokeOverrideMap(featureStrokes.overrides, { requireOverride: true })
    },
    originalSvgStroke: {
      color: Object.prototype.hasOwnProperty.call(originalSvgStroke, 'color')
        ? normalizeOptionalHexColor(originalSvgStroke.color)
        : defaults.originalSvgStroke.color,
      width: Object.prototype.hasOwnProperty.call(originalSvgStroke, 'width')
        ? normalizeStrokeWidth(originalSvgStroke.width)
        : defaults.originalSvgStroke.width
    },
    alignmentResetReceipt: cloneJsonValue(source.alignmentResetReceipt, null),
    featureCatalog: admittedFeatureCatalog(
      featureCatalog !== undefined ? featureCatalog : source.featureCatalog
    )
  };
};

const replacePlainObject = (target, source) => {
  Object.keys(target).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = value;
  });
};

/** @type {((request: Record<string, any>) => any) | null} */
let unmanagedConfigOverrideValidator = null;

export const setUnmanagedConfigOverrideValidator = (validator) => {
  unmanagedConfigOverrideValidator = typeof validator === 'function'
    ? validator
    : null;
};

const validateUnmanagedConfigOverrides = async ({
  mode,
  config = null,
  configOverrides = {},
  requireUnmanagedOnly = false
}) => {
  const managedPaths = managedConfigOverridePathsForMode(mode);
  const overrides = isPlainObject(configOverrides)
    ? cloneJsonData(configOverrides)
    : configOverrides;
  const overridePaths = isPlainObject(overrides) ? Object.keys(overrides) : [];
  if (
    config === null
    && overridePaths.length === 0
  ) {
    return {};
  }
  if (
    config === null
    && !requireUnmanagedOnly
    && overridePaths.every((path) => managedPaths.includes(path))
  ) {
    return {};
  }
  if (!unmanagedConfigOverrideValidator) {
    throw new Error('Configuration validation service is unavailable.');
  }
  const result = await unmanagedConfigOverrideValidator({
    mode,
    config,
    configOverrides: overrides,
    managedPaths,
    requireUnmanagedOnly
  });
  if (result?.result?.error && typeof result.result.error === 'object') {
    throw Object.assign(new Error('Configuration validation failed.'), result.result.error);
  }
  if (typeof result?.result?.error === 'string' && result.result.error.trim()) {
    throw new Error(result.result.error.trim());
  }
  if (!isPlainObject(result?.result?.overrides)) {
    throw new Error('Configuration validation returned an invalid preserved-settings result.');
  }
  return cloneJsonData(result.result.overrides);
};

/** @param {DrawingState} drawing */
export const applyEditorStateData = (
  drawing,
  editorState = {},
  { normalized: alreadyNormalized = false } = {}
) => {
  const normalized = alreadyNormalized
    ? editorState
    : normalizeEditorStateData(editorState);

  applyEditorArtifactData(normalized);
  applyDrawingEditorData(drawing, normalized);
};

// The shown Result's editor artifacts: its alignment Reset receipt, generated
// Legend inventory and colors, stroke defaults, and feature catalog.
/** @param {Record<string, any>} normalized The result of `normalizeEditorStateData`. */
const applyEditorArtifactData = (normalized) => {
  if (state.similarityAlignmentResetReceipt) {
    state.similarityAlignmentResetReceipt.value = normalized.alignmentResetReceipt ?? null;
  }
  state.originalLegendOrder.value = normalized.legend.originalOrder;
  state.originalLegendColors.value = normalized.legend.originalColors;
  state.originalSvgStroke.value = normalized.originalSvgStroke;
  if (state.featureCatalog) {
    state.featureCatalog.value = admittedFeatureCatalog(normalized.featureCatalog);
  }
};

// The Legend and stroke edits of a drawing, already normalized.
/**
 * @param {DrawingState} drawing
 * @param {Record<string, any>} normalized The result of `normalizeEditorStateData`.
 */
const applyDrawingEditorData = (drawing, normalized) => {
  drawing.legendEntries.value = normalized.legend.entries;
  drawing.deletedLegendEntries.value = normalized.legend.deletedEntries;
  drawing.dormantLegendEntries.value = normalized.legend.dormantEntries;
  replacePlainObject(drawing.legendColorOverrides, normalized.legend.colorOverrides);
  replacePlainObject(drawing.legendStrokeOverrides, normalized.legend.strokeOverrides);
  drawing.addedLegendCaptions.value = new Set(normalized.legend.addedCaptions);
  replacePlainObject(drawing.featureStrokeOverrides, normalized.featureStrokes.overrides);
};

// ---- Session 46: one drawing's slice (`modes.<m>`, PD-OI-086) ----

// A drawing's complete Session 46 slice: every registry row, as Save writes
// it. `selectedFeatureRecordIdx` is the Features-list record of the drawing's
// shown Result (0 for a mode that is not shown).
/**
 * @param {DrawingState} drawing
 * @param {'circular' | 'linear'} mode
 * @param {{ selectedFeatureRecordIdx?: number }} [options]
 */
export const buildModeSliceData = (drawing, mode, { selectedFeatureRecordIdx = 0 } = {}) => {
  const { legend, featureStrokes } = buildEditorStateData(drawing);
  return {
    config: cloneJsonData(buildConfigData(drawing)),
    features: {
      featureOverrides: featureOverridesForState(drawing.featureOverrides),
      featureColorOverrides: cloneJsonData(drawing.featureColorOverrides),
      featureVisibilityManualRules: normalizeFeatureVisibilityRulesForSession(drawing.featureVisibilityManualRules),
      labelOverrideRows: cloneJsonData(drawing.canonicalLabelOverrideRows.value),
      labelTextBulkOverrides: cloneJsonData(drawing.labelTextBulkOverrides)
    },
    editorState: {
      legend: {
        entries: [...legend.entries, ...legend.dormantEntries.map((entry) => ({ ...entry, dormant: true }))],
        deletedEntries: legend.deletedEntries,
        colorOverrides: legend.colorOverrides,
        strokeOverrides: legend.strokeOverrides,
        addedCaptions: legend.addedCaptions
      },
      featureStrokes
    },
    ui: {
      layoutPreferences: cloneJsonData(drawing.layoutPreferences[mode]),
      canvasPadding: { ...drawing.canvasPadding },
      pendingPaletteName: drawing.pendingPaletteName.value,
      pendingPaletteColors: cloneColors(drawing.pendingPaletteColors.value),
      linearTypographyLinked: Boolean(drawing.linearTypographyLinked.value),
      selectedFeatureRecordIdx
    }
  };
};

const CONFIG_SLICE_ROWS = MODE_SCOPED_SETTINGS.rows.filter((/** @type {{ domain: string }} */ row) => (
  row.domain === 'config' || row.domain.startsWith('config.')
));
// A drawing's settings on Load (plan 4.1): the projection of its mode's
// committed request, then each registry row its slice holds. A slice may
// omit any row; the drawing's defaults fill what neither holds.
/**
 * @param {Record<string, any> | null | undefined} projectedConfig
 * @param {Record<string, any> | null | undefined} sliceConfig
 * @returns {Record<string, any>}
 */
const overlayModeSliceConfig = (projectedConfig, sliceConfig) => {
  /** @type {Record<string, any>} */
  const config = isPlainObject(projectedConfig) ? cloneJsonData(projectedConfig) : {};
  const slice = isPlainObject(sliceConfig) ? /** @type {Record<string, any>} */ (sliceConfig) : {};
  CONFIG_SLICE_ROWS.forEach((/** @type {{ domain: string, path: string }} */ row) => {
    const container = row.domain === 'config' ? slice : slice[row.domain.slice('config.'.length)];
    if (!isPlainObject(container) || !Object.hasOwn(container, row.path)) return;
    if (row.domain === 'config') {
      config[row.path] = cloneJsonData(container[row.path]);
    } else {
      const domain = row.domain.slice('config.'.length);
      config[domain] = { ...(isPlainObject(config[domain]) ? config[domain] : {}), [row.path]: cloneJsonData(container[row.path]) };
    }
  });
  if (Object.hasOwn(slice, 'colors')) delete config.colorsAreOverrides;
  if (isPlainObject(config.adv)) delete config.adv.losatProgram;
  return config;
};

// Installs one drawing's slice over the projection of its mode's committed
// request (`projectedConfig`, or none); the drawing holds its mode's defaults
// before. `ui` keys apply as Session Load reads them: the pending palette
// only while Instant Preview is off, and the Linear typography link only
// while the two Linear sizes are equal.
/**
 * @param {DrawingState} drawing
 * @param {'circular' | 'linear'} mode
 * @param {Record<string, any>} slice
 * @param {{ projectedConfig?: Record<string, any> | null, resolveTrackPlacements?: boolean, applyCanvasPadding?: boolean }} [options]
 *   `applyCanvasPadding: false` leaves the padding to the caller (Load pads the shown Result after it mounts).
 */
const applyModeSliceData = (drawing, mode, slice, {
  projectedConfig = null, resolveTrackPlacements = true, applyCanvasPadding = true
} = {}) => {
  const config = overlayModeSliceConfig(projectedConfig, slice.config);
  if (Object.keys(config).length) applyConfigData(drawing, config, { resolveTrackPlacements });
  applyDrawingFeatureData(drawing, isPlainObject(slice.features) ? slice.features : {});
  applyDrawingEditorData(drawing, normalizeEditorStateData(isPlainObject(slice.editorState) ? slice.editorState : {}));
  const ui = isPlainObject(slice.ui) ? slice.ui : {};
  if (isPlainObject(ui.layoutPreferences)) {
    replaceLayoutPreferences(drawing.layoutPreferences, {
      ...cloneJsonData(drawing.layoutPreferences), [mode]: cloneJsonData(ui.layoutPreferences)
    });
  }
  if (applyCanvasPadding && isPlainObject(ui.canvasPadding)) {
    ['top', 'right', 'bottom', 'left'].forEach((side) => {
      drawing.canvasPadding[side] = Number(ui.canvasPadding[side]) || 0;
    });
  }
  restorePendingPaletteFromSession(drawing, ui);
  reconcileImportedLinearTypographyLink({ adv: drawing.adv, linked: drawing.linearTypographyLinked, ui });
};

const SESSION_FORMAT_ERROR = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'SESSION_FORMAT' });
const SESSION_FIELDS_ERROR = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
const validateSessionVersion = version => {
  if (!Number.isInteger(version)) throw SESSION_FORMAT_ERROR();
  if (version > SESSION_VERSION || !SUPPORTED_SESSION_VERSIONS.has(version)) {
    throw SESSION_FIELDS_ERROR();
  }
};

const normalizeSessionData = (data) => {
  if (!isPlainObject(data) || data.format !== 'gbdraw-session') throw SESSION_FORMAT_ERROR();
  const version = data.version;
  validateSessionVersion(version);
  if (version >= CURRENT_AUTHORITY_SESSION_MIN_VERSION && Object.prototype.hasOwnProperty.call(data, 'files')) {
    throw new Error(
      `Session version ${version} cannot contain legacy files; use resources and webFiles.`
    );
  }
  if (version >= 31) {
    if (!isPlainObject(data.renderRequest)) {
      throw new Error(`Session version ${version} requires a canonical renderRequest object.`);
    }
    if (!isPlainObject(data.resources)) {
      throw new Error(`Session version ${version} requires a canonical resources object.`);
    }
  }

  return {
    ...data,
    editorState: normalizeEditorStateData(data.editorState)
  };
};

const migrateLegacyFeatureRenderingConfig = (configData, legacy) => {
  if (!legacy || !isPlainObject(configData)) return configData;
  const adv = isPlainObject(configData.adv) ? configData.adv : null;
  if (!adv) return configData;
  const features = Array.isArray(adv.features) ? adv.features : null;
  if (features && !features.includes('repeat_region')) return configData;
  const featureShapes = isPlainObject(adv.feature_shapes) ? adv.feature_shapes : {};
  if (Object.prototype.hasOwnProperty.call(featureShapes, 'repeat_region')) return configData;
  return {
    ...configData,
    adv: {
      ...adv,
      feature_shapes: { ...featureShapes, repeat_region: 'rectangle' }
    }
  };
};

const sessionArtifactEntries = (data, field) => {
  const container = data[field];
  if (container === undefined || container === null) return [];
  if (!isPlainObject(container)) throw SESSION_FIELDS_ERROR();
  const entries = Object.prototype.hasOwnProperty.call(container, 'entries')
    ? container.entries
    : [];
  if (!Array.isArray(entries)) throw SESSION_FIELDS_ERROR();
  return entries;
};

// A raw key names one cache entry; a repeated or missing key is a malformed
// Session (R6: classified, without document values).
const rejectInvalidLosatCacheKeys = (entries, { requireKey = false } = {}) => {
  const seen = new Set();
  for (const entry of entries) {
    const key = isPlainObject(entry) && typeof entry.key === 'string'
      ? entry.key
      : '';
    if (!key) {
      if (requireKey) throw SESSION_FIELDS_ERROR();
      continue;
    }
    if (seen.has(key)) throw SESSION_FIELDS_ERROR();
    seen.add(key);
  }
};

function* sessionLosatArtifactSteps(data, sourceSessionVersion) {
  if (sourceSessionVersion < CURRENT_ARTIFACT_SESSION_MIN_VERSION) return;
  const rawEntries = sessionArtifactEntries(data, 'losatCache');
  const derivedEntries = sessionArtifactEntries(data, 'losatDerivedCache');
  const manifest = data.proteinIdentityManifest;
  rejectInvalidLosatCacheKeys(rawEntries, { requireKey: true });
  rejectInvalidLosatCacheKeys(derivedEntries);

  if (sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION) {
    recordStructuralMetric('currentSessionPreflightProteinManifestValidationCount');
  }
  const identityIndex = buildValidatedProteinIdentityIndex(manifest);
  if (!identityIndex) throw SESSION_FIELDS_ERROR();
  let invalidDerivedEntry = false;
  try {
    for (const entry of rawEntries) {
      const classification = classifyRawLosatCacheEntry(entry);
      if (!['protein-current', 'nucleotide-current'].includes(classification)) {
        throw SESSION_FIELDS_ERROR();
      }
      yield;
      if (classification !== 'protein-current') continue;
      if (sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION) {
        recordStructuralMetric('currentSessionPreflightProteinRawTextValidationCount');
      }
      if (
        !validateProteinRawEntryReferences(entry, manifest, { identityIndex })
      ) {
        throw SESSION_FIELDS_ERROR();
      }
    }
    invalidDerivedEntry = derivedEntries.some(
      (entry) => !validateDerivedProteinReferences(entry, manifest, { identityIndex })
    );
  } finally {
    releaseValidatedProteinIdentityIndex(identityIndex);
  }
  if (invalidDerivedEntry) throw SESSION_FIELDS_ERROR();
};

export const validateSessionLosatArtifacts = (data, sourceSessionVersion) => {
  for (const _step of sessionLosatArtifactSteps(data, sourceSessionVersion)) { /* exhaust validation */ }
};

const validateSessionLosatArtifactsForImport = async (data, sourceSessionVersion) => {
  let deadline = performance.now() + 16;
  for (const _step of sessionLosatArtifactSteps(data, sourceSessionVersion)) {
    if (performance.now() >= deadline) {
      await new Promise((resolve) => setTimeout(resolve, 0));
      deadline = performance.now() + 16;
    }
  }
};

export const buildSessionLegacyArtifacts = ({
  legacyRawCandidates,
  legacyDerivedEvidence
}) => {
  const legacyArtifacts = {};
  if (legacyRawCandidates?.entries?.length) {
    legacyArtifacts.proteinRawCandidates = legacyRawCandidates;
  }
  if (legacyDerivedEvidence?.entries?.length) {
    legacyArtifacts.proteinDerivedEvidence = legacyDerivedEvidence;
  }
  return Object.keys(legacyArtifacts).length > 0 ? legacyArtifacts : null;
};

// The Web writers of Sessions 27–33 saved every schema-4 Circular slot row with
// `spacing: null`. With Custom Track Slots off that null is lossless and is
// dropped; any other obsolete slot field fails validateImportedCircularTrackSlots.
const withoutLegacyNullCircularSlotSpacing = (configData) => {
  const adv = configData?.adv;
  if (!isPlainObject(adv) || adv.circular_track_slots_enabled
    || adv.circular_track_slots_schema_version !== CIRCULAR_TRACK_SLOT_SCHEMA_VERSION
    || !Array.isArray(adv.circular_track_slots)) return configData;
  return {
    ...configData,
    adv: {
      ...adv,
      circular_track_slots: adv.circular_track_slots.map((slot) => {
        if (!isPlainObject(slot) || slot.spacing !== null) return slot;
        const { spacing: _spacing, ...current } = slot;
        return current;
      })
    }
  };
};

/**
 * @param {Record<string, any>} data An unvalidated Session 27-39.
 * @param {number} sourceSessionVersion
 * @returns {Record<string, any>} Its flat draft in the current option names;
 *   Load splits it into mode slices.
 */
const migrateSessionDataToCurrent = (data, sourceSessionVersion) => {
  const readsLegacyOptionValues = sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION;
  const migratedOptions = readsLegacyOptionValues
    ? migratePersistedWebOptionValues(data.config)
    : data.config;
  const circularSlotConfig = migrateImportedCircularTrackSlots(migratedOptions);
  return ({
    ...data,
    version: SESSION_VERSION,
    config: migrateLegacyFeatureRenderingConfig(
      migrateImportedLinearTrackSlots(
        sourceSessionVersion <= 33
          ? withoutLegacyNullCircularSlotSpacing(circularSlotConfig)
          : circularSlotConfig,
        sourceSessionVersion
      ),
      sourceSessionVersion <= 33
    )
  });
};

const LEGACY_CONFIG_KEYS = new Set([
  'form',
  'adv',
  'losat',
  'cliOptions',
  'colors',
  'palette',
  'rules',
  'qualifierPriorityRules',
  'filterMode',
  'whitelist',
  'blacklistText',
  'blastSource',
  'losatProgram',
  'circularConservation'
]);

const isLegacyConfigPayload = (data) =>
  isPlainObject(data) &&
  !Object.prototype.hasOwnProperty.call(data, 'format') &&
  Object.keys(data).some((key) => LEGACY_CONFIG_KEYS.has(key));

/** @param {DrawingState} drawing */
const applyLegacyConfigPayload = (drawing, data) => {
  const circularSlotSchema = data?.adv?.circular_track_slots_schema_version;
  const migratedOptions = circularSlotSchema === CIRCULAR_TRACK_SLOT_SCHEMA_VERSION
    ? data
    : migratePersistedWebOptionValues(data);
  const migrated = migrateLegacyFeatureRenderingConfig(
    migrateImportedLinearTrackSlots(
      migrateImportedCircularTrackSlots(migratedOptions)
    ),
    true
  );
  validateImportedCircularTrackSlots(migrated);
  validateImportedLinearTrackSlots(migrated);
  state.suppressCircularMultiRecordDefaults.value = shouldSuppressCircularMultiRecordDefaults(drawing, migrated.form);
  applyConfigData(drawing, migrated);
  restorePaletteStateAfterConfigImport(drawing);
};

/** @param {DrawingState} drawing */
const shouldSuppressCircularMultiRecordDefaults = (drawing, incomingForm) => {
  if (state.mode.value !== 'circular') return false;
  if (!incomingForm || typeof incomingForm !== 'object' || Array.isArray(incomingForm)) return false;
  if (!Object.prototype.hasOwnProperty.call(incomingForm, 'multi_record_canvas')) return false;
  return drawing.form.multi_record_canvas === false && incomingForm.multi_record_canvas === true;
};

const overlayCanonicalObject = (stored, canonical) => {
  const merged = isPlainObject(stored) ? cloneJsonData(stored) : {};
  if (!isPlainObject(canonical)) return merged;
  Object.entries(canonical).forEach(([key, value]) => {
    if (['__proto__', 'constructor', 'prototype'].includes(key)) return;
    merged[key] = isPlainObject(value)
      ? overlayCanonicalObject(merged[key], value)
      : cloneJsonData(value);
  });
  return merged;
};

const restoreStoredNonCanonicalConfig = (
  projectedConfig,
  storedConfig,
  { hasCanonicalProteinPipeline = false } = {}
) => {
  const restored = cloneJsonData(projectedConfig);
  if (!isPlainObject(storedConfig)) return restored;
  // Keep canonical drawing values authoritative; only supplement state that the
  // canonical request does not currently represent.
  [
    'losat',
    'blastSource',
    'losatProgram',
    'cliOptions',
    'paletteInstantPreviewEnabled',
    'modeProfiles',
    'linearRecordLayout',
    'linearComparisonPlan',
    'webEdits'
  ].forEach((key) => {
    if (hasCanonicalProteinPipeline && key === 'losat') {
      if (!isPlainObject(storedConfig.losat)) return;
      const canonicalLosat = isPlainObject(restored.losat)
        ? cloneJsonData(restored.losat)
        : {};
      restored.losat = overlayCanonicalObject(
        storedConfig.losat,
        canonicalLosat
      );
      return;
    }
    if (
      hasCanonicalProteinPipeline &&
      ['blastSource', 'losatProgram'].includes(key)
    ) return;
    if (Object.prototype.hasOwnProperty.call(storedConfig, key)) {
      restored[key] = cloneJsonData(storedConfig[key]);
    }
  });
  const storedAdv = storedConfig.adv;
  if (isPlainObject(storedAdv) && isPlainObject(restored.adv)) {
    [
      'rich_feature_popup',
      'feature_width_circular',
      'depth_width_circular',
      'gc_content_width_circular',
      'gc_content_radius_circular',
      'gc_skew_width_circular',
      'gc_skew_radius_circular'
    ].forEach((key) => {
      if (Object.prototype.hasOwnProperty.call(storedAdv, key)) {
        restored.adv[key] = cloneJsonData(storedAdv[key]);
      }
    });
  }
  const storedConservation = storedConfig?.circularConservation;
  if (!isPlainObject(storedConservation)) {
    return restored;
  }
  if (!isPlainObject(restored?.circularConservation)) {
    restored.circularConservation = cloneJsonData(storedConservation);
    return restored;
  }

  ['enabled', 'source', 'losat_program', 'subject_gencode'].forEach((key) => {
    if (Object.prototype.hasOwnProperty.call(storedConservation, key)) {
      restored.circularConservation[key] = cloneJsonData(storedConservation[key]);
    }
  });
  if (
    Array.isArray(restored.circularConservation.series) &&
    Array.isArray(storedConservation.series)
  ) {
    restored.circularConservation.series = restored.circularConservation.series.map((entry, index) => {
      const storedEntry = storedConservation.series[index];
      if (!isPlainObject(entry) || !isPlainObject(storedEntry)) return entry;
      if (!Object.prototype.hasOwnProperty.call(storedEntry, 'losat_gencode')) return entry;
      return { ...entry, losat_gencode: cloneJsonData(storedEntry.losat_gencode) };
    });
  }
  return restored;
};

// Version-40 config owns the active controls and the next Generate. The
// canonical projection owns the committed artifact and only fills fields the
// current writer deliberately omitted.
export const restoreCurrentWriterActiveConfig = ({
  mode,
  projectedConfig,
  storedConfig
}) => {
  if (!isPlainObject(projectedConfig)) {
    throw new Error('Current session is missing its canonical configuration projection.');
  }
  validateCurrentWriterActiveConfig({ mode, storedConfig, scopedDrafts: true });
  const restored = cloneJsonData(projectedConfig);
  const restoredDomains = [];
  CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS.forEach((domain) => {
    if (!Object.prototype.hasOwnProperty.call(storedConfig, domain)) return;
    restored[domain] = cloneJsonData(storedConfig[domain]);
    restoredDomains.push(domain);
  });
  if (Object.prototype.hasOwnProperty.call(storedConfig, 'colors')) {
    delete restored.colorsAreOverrides;
  }
  if (isPlainObject(restored.adv)) delete restored.adv.losatProgram;

  const projectedRulerLabelFontSize = projectedConfig?.adv?.ruler_label_font_size;
  if (
    projectedRulerLabelFontSize !== undefined
    && isPlainObject(restored.adv)
    && !Object.prototype.hasOwnProperty.call(restored.adv, 'ruler_label_font_size')
  ) {
    restored.adv.ruler_label_font_size = cloneJsonData(projectedRulerLabelFontSize);
  }

  // Current sessions written before the anchor control became editable omitted
  // this value. Preserve their committed projection only when active config has
  // no value of its own; new sessions store the explicit editor value above.
  const projectedAnchorMode = projectedConfig?.losat?.blastp?.collinearAnchorMode;
  if (
    projectedAnchorMode !== undefined &&
    isPlainObject(restored.losat) &&
    isPlainObject(restored.losat.blastp) &&
    !Object.prototype.hasOwnProperty.call(restored.losat.blastp, 'collinearAnchorMode')
  ) {
    restored.losat.blastp.collinearAnchorMode = cloneJsonData(projectedAnchorMode);
  }

  recordStructuralMetric('currentWriterActiveConfigRestoreCount', 1, {
    domains: restoredDomains
  });
  recordStructuralMetric('activeConfigCanonicalOverwriteCount', 0, {
    domains: restoredDomains
  });
  return restored;
};

const validateCurrentWriterFeatureCatalog = async (data, { adopt = false } = {}) => {
  const results = normalizeLogicalResults(
    (Array.isArray(data.results) ? data.results : []).map((result, index) => ({
      name: result?.name || `Result ${index + 1}`,
      content: result?.content || ''
    }))
  );
  const catalog = data.editorState?.featureCatalog ?? null;
  if (catalog === null) return null;
  return validateFeatureCatalogForImport(catalog, results, {
    adopt,
    mode: data.renderRequest?.mode || ''
  });
};

const preflightSessionImport = async (sessionData) => {
  const sourceSessionVersion = sessionData?.version;
  validateSessionVersion(sourceSessionVersion);
  const rawData = await convertMainSessionComparisonFrames(sessionData);
  const currentSession = sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION;
  /** @type {ReturnType<typeof adoptCurrentSessionDocument> | null} */
  let adoptedSession = null;
  /** @type {ReturnType<typeof adoptCurrentSessionResources> | null} */
  let currentResourceTable = null;
  /** @type {Awaited<ReturnType<typeof validateCurrentWriterFeatureCatalog>> | null} */
  let validatedFeatureCatalog = null;
  /** @type {Awaited<ReturnType<typeof validateCurrentWriterFeatureCatalog>> | null} */
  let otherModeCatalog = null;
  /** @type {ReturnType<typeof normalizeEditorStateData> | null} */
  let normalizedEditorState = null;
  let normalizedData;

  if (currentSession) {
    if (!isPlainObject(rawData) || rawData.format !== 'gbdraw-session') throw SESSION_FORMAT_ERROR();
    recordSessionLifecycleEvent('session-authority-validation-start');
    adoptedSession = adoptCurrentSessionDocument(rawData, sourceSessionVersion);
    recordSessionLifecycleEvent('session-authority-validation-end');
    recordSessionLifecycleEvent('feature-catalog-validation-start');
    validatedFeatureCatalog = await validateCurrentWriterFeatureCatalog(rawData, {
      adopt: true
    });
    // E1: the other mode's set is validated with its own mode.
    if (isPlainObject(rawData.otherModeResult)) {
      otherModeCatalog = await validateCurrentWriterFeatureCatalog(rawData.otherModeResult, { adopt: true });
    }
    recordSessionLifecycleEvent('feature-catalog-validation-end');
    recordSessionLifecycleEvent('editor-state-normalization-start');
    normalizedEditorState = normalizeEditorStateData(rawData.editorState, {
      featureCatalog: validatedFeatureCatalog
    });
    recordSessionLifecycleEvent('editor-state-normalization-end');
    recordSessionLifecycleEvent('resource-table-adoption-start');
    currentResourceTable = adoptCurrentSessionResources(rawData.resources);
    recordSessionLifecycleEvent('resource-table-adoption-end');
    normalizedData = rawData;
  } else {
    validateSessionAuthorityInventory(rawData, sourceSessionVersion);
    normalizedData = normalizeSessionData(rawData);
    migrateImportedLinearTrackSlots(normalizedData.config, sourceSessionVersion);
  }

  const promotedData = !currentSession && (
    sourceSessionVersion >= 31 &&
    Number(normalizedData.renderRequest?.schema) === 2
  )
    ? promoteGallerySessionToCurrent(normalizedData)
    : normalizedData;
  if (currentSession) recordSessionLifecycleEvent('losat-artifact-validation-start');
  await validateSessionLosatArtifactsForImport(promotedData, sourceSessionVersion);
  if (currentSession) recordSessionLifecycleEvent('losat-artifact-validation-end');
  const data = currentSession
    ? promotedData
    : migrateSessionDataToCurrent(promotedData, sourceSessionVersion);
  const comparisonClassification = classifyImportedComparisonIntent(/** @type {Record<string, any>} */ ({
    renderRequest: data.renderRequest,
    resources: data.resources
  }));
  const projectionRenderRequest = (
    data.renderRequest?.mode === 'linear'
    && comparisonClassification.disposition
      !== IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE
  )
    ? { ...data.renderRequest, comparisons: [] }
    : data.renderRequest;
  let currentStoredConfig = sourceSessionVersion < TYPED_DRAFT_SESSION_VERSION
    ? withCurrentLinearLabelVisibility(data.config)
    : data.config;
  if (sourceSessionVersion < TYPED_DRAFT_SESSION_VERSION && isPlainObject(currentStoredConfig)
    && Object.prototype.hasOwnProperty.call(currentStoredConfig, 'recordDisplayDrafts')) {
    currentStoredConfig = {
      ...currentStoredConfig,
      recordDisplayDrafts: migrateLegacyRecordDisplayDrafts(
        currentStoredConfig.recordDisplayDrafts
      )
    };
  }
  // A Session 41-44 placement row reached both modes; the migration names
  // its mode (Main rows both), and the Session 46 split moves it there (R2).
  if (sourceSessionVersion < SESSION_VERSION && isPlainObject(currentStoredConfig)
    && Object.prototype.hasOwnProperty.call(currentStoredConfig, 'featurePlacementOverrides')) {
    currentStoredConfig = {
      ...currentStoredConfig,
      featurePlacementOverrides: migrateSessionFeaturePlacements(currentStoredConfig.featurePlacementOverrides)
    };
  }
  // Session 46 keeps each mode's draft in its slice (PD-OI-086): the committed
  // request's projection reads the slice of the request's mode as a Session
  // 40-44 read its stored draft, and a Session without slices is a CLI or
  // Python writer's.
  const modeScopedSession = sourceSessionVersion >= MODE_SCOPED_SESSION_VERSION;
  const committedSliceConfig = modeScopedSession
    ? data.modes?.[data.renderRequest?.mode === 'linear' ? 'linear' : 'circular']?.config
    : undefined;
  const hasStoredDraft = modeScopedSession
    ? isPlainObject(committedSliceConfig)
    : Object.prototype.hasOwnProperty.call(data, 'config');
  const runtimeStoredConfig = currentSession && hasStoredDraft
    ? migrateImportedLinearTrackSlots(
        migrateImportedCircularTrackSlots(modeScopedSession ? committedSliceConfig : currentStoredConfig),
        sourceSessionVersion
      )
    : data.config;
  if (currentSession) recordSessionLifecycleEvent('canonical-request-projection-start');
  const settingsOnly = isSettingsOnlySessionDocument(data);
  const canonicalProjection = settingsOnly
    ? projectSettingsOnlySession(data, currentResourceTable)
    : sourceSessionVersion >= 31
    ? projectCanonicalSessionRequest({
        renderRequest: projectionRenderRequest,
        resources: data.resources,
        webFiles: data.webFiles,
        legacyFiles: data.files,
        storedConfig: runtimeStoredConfig,
        initializeCliInputs: !hasStoredDraft
          && data.cliInvocation?.generatedBy === 'gbdraw',
        fileBindings: data.cliInvocation?.fileBindings,
        linearTrackSlotSchemaVersion: sourceSessionVersion <= LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION
          ? LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION
          : LINEAR_TRACK_SLOT_SCHEMA_VERSION,
        repairInvalidComparisonHeight: sourceSessionVersion >= 31 && sourceSessionVersion <= 33,
        repairLegacyTableRows: sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION,
        sessionResourceTable: currentResourceTable,
        deferResourceContent: currentSession,
        adoptCanonicalPayloads: currentSession
      })
    : null;
  if (currentSession) recordSessionLifecycleEvent('canonical-request-projection-end');
  if (canonicalProjection?.pipelineState?.legacySimilarityAlignment) {
    const promotedRequest = promoteCanonicalRenderRequestToCurrent(data.renderRequest, {
      featureCatalog: data.editorState?.featureCatalog || null,
      legacyOrthogroupState: data.orthogroupState || null
    });
    canonicalProjection.config.linearRecordLayout = {
      ...canonicalProjection.config.linearRecordLayout,
      recordTranslations: promotedRequest.layout.recordTranslations,
      similarityAlignment: promotedRequest.layout.similarityAlignment
    };
  }
  let restoredConfig = canonicalProjection
    ? sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION
      ? cloneJsonData(canonicalProjection.config)
      : {
          ...restoreStoredNonCanonicalConfig(
            canonicalProjection.config,
            runtimeStoredConfig,
            { hasCanonicalProteinPipeline: Boolean(canonicalProjection.pipelineState) }
          ),
          rules: applySpecificRuleProvenance(
            canonicalProjection.config.rules,
            runtimeStoredConfig?.rules
          )
        }
    : data.config;
  if (sourceSessionVersion < TYPED_DRAFT_SESSION_VERSION && restoredConfig) {
    restoredConfig = withCurrentLinearLabelVisibility(restoredConfig);
  }
  if (sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION && restoredConfig) {
    const sourceStoredConfig = isPlainObject(normalizedData.config)
      ? normalizedData.config
      : {};
    const forceWebDraft = isPlainObject(sourceStoredConfig.linearRecordLayout)
      || !isPlainObject(sourceStoredConfig.cliOptions);
    const comparisonMigrationConfig = cloneJsonData(restoredConfig);
    if (
      !Object.prototype.hasOwnProperty.call(comparisonMigrationConfig, 'blastSource')
      && normalizedData.ui?.blastSource
    ) {
      comparisonMigrationConfig.blastSource = normalizedData.ui.blastSource;
    }
    const migratedComparisonDraft = migrateLegacyLinearComparisonDraft({
      config: comparisonMigrationConfig,
      filesData: canonicalProjection?.files || data.files || {},
      forceWebDraft
    });
    restoredConfig = migratedComparisonDraft.config;
    if (canonicalProjection) {
      canonicalProjection.files = migratedComparisonDraft.filesData;
    } else {
      data.files = migratedComparisonDraft.filesData;
    }
  }
  if (canonicalProjection && !modeScopedSession && sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION
    && hasStoredDraft) {
    recordSessionLifecycleEvent('current-draft-validation-start');
    restoredConfig = restoreCurrentWriterActiveConfig({
      mode: canonicalProjection.mode,
      projectedConfig: canonicalProjection.config,
      storedConfig: runtimeStoredConfig
    });
    recordSessionLifecycleEvent('current-draft-validation-end');
  }
  // Session 46: each slice is a current-writer draft of its own mode; the
  // committed mode's draft is its slice over its request's projection.
  if (modeScopedSession) {
    recordSessionLifecycleEvent('current-draft-validation-start');
    SLICE_MODES.forEach((mode) => {
      const slice = data.modes?.[mode];
      if (slice === undefined) return;
      validateModeSliceFields(slice);
      if (isPlainObject(slice.config)) {
        validateCurrentWriterActiveConfig({
          mode,
          storedConfig: migrateImportedLinearTrackSlots(migrateImportedCircularTrackSlots(
            { form: {}, adv: {}, ...slice.config }
          ), sourceSessionVersion)
        });
      }
    });
    restoredConfig = canonicalProjection && !settingsOnly
      ? overlayModeSliceConfig(canonicalProjection.config, runtimeStoredConfig)
      : null;
    if (restoredConfig) {
      // As for Sessions 40-44: the stored draft owns the active controls and the
      // projection fills only what it omits.
      const restoredDomains = CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS
        .filter((domain) => isPlainObject(runtimeStoredConfig) && Object.hasOwn(runtimeStoredConfig, domain));
      recordStructuralMetric('currentWriterActiveConfigRestoreCount', 1, { domains: restoredDomains });
      recordStructuralMetric('activeConfigCanonicalOverwriteCount', 0, { domains: restoredDomains });
    }
    recordSessionLifecycleEvent('current-draft-validation-end');
  }
  if (!canonicalProjection && restoredConfig) {
    hydrateMissingMultiRecordPositionsFromCliInvocation(restoredConfig, data.cliInvocation);
  }
  // Saved omission means historical OFF, independently of fresh/reset defaults.
  restoredConfig = {
    ...restoredConfig,
    form: { keep_definition_left_aligned: false, ...restoredConfig?.form }
  };
  const hasCurrentStoredUnmanagedOverrides = currentSession && hasStoredDraft
    && isPlainObject(runtimeStoredConfig)
    && Object.prototype.hasOwnProperty.call(
      runtimeStoredConfig,
      'unmanagedConfigOverrides'
    );
  const unmanagedConfigValidation = canonicalProjection
    ? hasCurrentStoredUnmanagedOverrides
      ? {
          mode: canonicalProjection.mode,
          configOverrides: runtimeStoredConfig.unmanagedConfigOverrides,
          requireUnmanagedOnly: true
        }
      : {
          mode: canonicalProjection.mode,
          config: projectionRenderRequest?.diagramOptions?.config ?? null,
          configOverrides:
            projectionRenderRequest?.diagramOptions?.configOverrides ?? {}
        }
    : null;
  if (restoredConfig) {
    const activeDepthTrackCount = Array.isArray(restoredConfig?.adv?.depth_tracks)
      ? restoredConfig.adv.depth_tracks.length
      : null;
    validateImportedCircularTrackSlots(restoredConfig, {
      depthTrackCount: canonicalProjection?.mode === 'circular' &&
        restoredConfig.adv?.circular_track_slots_enabled
        ? activeDepthTrackCount
        : null
    });
    validateImportedLinearTrackSlots(restoredConfig, {
      depthTrackCount: canonicalProjection?.mode === 'linear' &&
        restoredConfig.adv?.linear_track_slots_enabled
        ? activeDepthTrackCount
        : null
    });
  }
  if (currentSession) recordSessionLifecycleEvent('artifact-projection-start');
  const projectionResult = canonicalProjection
    ? (() => {
        const artifactState = projectArtifactState(data);
        if (currentSession) artifactState.editorState = normalizedEditorState;
        artifactState.legacySimilarityAlignment =
          canonicalProjection.pipelineState?.legacySimilarityAlignment || null;
        return {
          documentMetadata: projectDocumentMetadata(data),
          renderState: {
            mode: canonicalProjection.mode,
            inputType: canonicalProjection.inputType,
            config: restoredConfig,
            semanticFeatureState: canonicalProjection.semanticFeatureState
          },
          editorMetadata: projectWebOnlyEditorMetadata(data),
          artifactState,
          restoredFiles: canonicalProjection.files,
          validatedFeatureCatalog
        };
      })()
    : null;
  if (currentSession) recordSessionLifecycleEvent('artifact-projection-end');
  return {
    data,
    sourceSessionVersion,
    canonicalProjection,
    restoredConfig,
    projectionResult,
    adoptedCanonicalSession: adoptedSession?.canonical || null,
    currentResourceTable,
    otherModeCatalog,
    comparisonClassification,
    unmanagedConfigValidation
  };
};

// E1: the match-sequence sources one saved Result set (the top-level set or
// `otherModeResult`) needs beyond its catalog's: when the catalog does not
// cover the set's request, or a Session older than 40 has none, the sources
// are read again from the restored files in the set's own mode.
/**
 * @param {{ mode: 'circular' | 'linear', catalog: FeatureCatalog | null, renderRequest: Record<string, any> | undefined,
 *   olderSession: boolean, settingsOnly: boolean, cInputType: string, lInputType: string,
 *   files: Record<string, any>, linearSeqs: Record<string, any>[], circularConservation: Record<string, any> }} set
 * @returns {Promise<{ restored: Record<string, any>[], error: unknown }>}
 */
const restoreLoadedSetSequenceSources = async ({
  mode, catalog, renderRequest, olderSession, settingsOnly, cInputType, lInputType, files, linearSeqs,
  circularConservation
}) => {
  const coverage = catalog
    ? analyzeCatalogSequenceSourceCoverage({
        mode,
        catalogFeatureState: catalog,
        renderRequest,
        comparisonSourceAvailability: mode === 'circular'
          ? resolveCircularComparisonSequenceAvailability({ files, circularConservation })
          : undefined
      })
    : null;
  if (settingsOnly || (!olderSession && coverage?.complete)) return { restored: [], error: null };
  recordStructuralMetric('sourceRecoveryCount');
  try {
    return {
      restored: await buildRestoredMatchSequenceSources({
        mode, cInputType, lInputType, files, linearSeqs, circularConservation
      }),
      error: null
    };
  } catch (sequenceError) {
    console.warn('Session match sequence preparation failed.', normalizeUserFacingError(sequenceError));
    return { restored: [], error: sequenceError };
  }
};

// E1: the Results of one saved Result set of a current Session (the top-level
// set or `otherModeResult`), admitted with that set's catalog in its mode.
/**
 * @param {Record<string, any>[]} logicalResults
 * @param {{ featureCatalog: FeatureCatalog, mode: 'circular' | 'linear', selectedFeatureTypes: readonly string[] | null | undefined }} set
 */
const admitLoadedSetResults = (logicalResults, { featureCatalog, mode, selectedFeatureTypes }) => (
  admitCurrentSessionResults(
    createCurrentSessionResultSource(
      logicalResults,
      admitFeatureCatalog(featureCatalog, logicalResults, { adopt: true, mode })
    ),
    { mutationPlan: createEmptySvgMutationPlan(logicalResults.length), selectedFeatureTypes }
  )
);

// E1: one saved Result set of a current Session (the top-level set or
// `otherModeResult`) as the artifact slot of its mode, built off-line during
// Load: its Results admitted with its catalog (or the Results Load already
// admitted), its request adopted as that mode's committed Session over the
// shared resource table, and its match-sequence sources.
/**
 * @param {{ renderRequest: Record<string, any>, results: Record<string, any>[], editorState?: unknown,
 *   ui?: unknown, runMetadata?: unknown }} set
 * @param {{ featureCatalog: FeatureCatalog, admittedResults?: Record<string, any>[] | null,
 *   restoredSequenceSources: Record<string, any>[], resources: Record<string, any>, resourceTable: any,
 *   webFiles: any, retainedBytes: number }} options
 */
const buildLoadedArtifactSlot = (set, {
  featureCatalog, admittedResults = null, restoredSequenceSources, resources, resourceTable, webFiles, retainedBytes
}) => {
  const mode = set.renderRequest.mode === 'linear' ? 'linear' : 'circular';
  const results = admittedResults || admitLoadedSetResults(normalizeLogicalResults(set.results.map(
    (/** @type {Record<string, any>} */ result, /** @type {number} */ index) => ({
      name: result?.name || `Result ${index + 1}`,
      content: result?.content || ''
    })
  )), { featureCatalog, mode, selectedFeatureTypes: set.renderRequest.diagramOptions?.selectedFeaturesSet });
  const features = featureStateFromCatalog(featureCatalog, { mode });
  const editorState = normalizeEditorStateData(
    isPlainObject(set.editorState) ? /** @type {Record<string, any>} */ (set.editorState) : {},
    { featureCatalog }
  );
  const ui = isPlainObject(set.ui) ? /** @type {Record<string, any>} */ (set.ui) : {};
  const runMetadata = isPlainObject(set.runMetadata) ? /** @type {Record<string, any>} */ (set.runMetadata) : {};
  return createArtifactSlot({
    mode,
    values: {
      results,
      selectedResultIndex: Number.isInteger(ui.selectedResultIndex) && ui.selectedResultIndex >= 0
        ? Math.min(ui.selectedResultIndex, Math.max(0, results.length - 1)) : 0,
      featureCatalog,
      extractedFeatures: features.extractedFeatures || [],
      biologicalFeatures: features.biologicalFeatures || [],
      featureRecordIds: features.featureRecordIds || [],
      orthogroups: features.orthogroups || [],
      featureOrthogroupIndex: features.featureOrthogroupIndex || new Map(),
      collinearGroups: features.collinearGroups || [],
      trackSlotResolvedGeometry: cloneJsonData(runMetadata.trackSlotGeometry ?? null),
      annotationWarnings: cloneJsonData(runMetadata.annotationWarnings || []),
      featureIdentityNotices: cloneJsonData(runMetadata.featureIdentityNotices || []),
      comparisonWarnings: cloneJsonData(runMetadata.comparisonWarnings || []),
      generatedLegendPosition: normalizeLegendPosition(
        ui.generatedLegendPosition, mode === 'linear' ? 'bottom' : 'left'
      ),
      generatedMode: mode,
      generatedMultiRecordCanvas: Boolean(ui.generatedMultiRecordCanvas),
      generatedCircularPlotTitlePosition: hasStoredLayoutValue(ui.generatedCircularPlotTitlePosition)
        ? normalizeCircularPlotTitlePosition(ui.generatedCircularPlotTitlePosition)
        : normalizeCircularPlotTitlePosition(ui.circularPlotTitlePosition),
      appliedPaletteName: String(ui.appliedPaletteName || 'default'),
      appliedPaletteColors: cloneColors(ui.appliedPaletteColors),
      similarityAlignmentResetReceipt: editorState.alignmentResetReceipt ?? null,
      originalLegendColors: editorState.legend.originalColors,
      originalSvgStroke: editorState.originalSvgStroke
    },
    legendInventory: editorState.legend.originalOrder,
    matchSequenceOwner: state.matchSequenceRegistry?.buildSourceOwner?.([
      ...(features.sequenceSources || []), ...restoredSequenceSources
    ]) || null,
    runtimeState: { canonical: {
      committedCanonicalSession: adoptRuntimeCanonicalSession({
        renderRequest: set.renderRequest, resources, webFiles: isPlainObject(webFiles) ? webFiles : {}
      }),
      activeSessionResourceTable: resourceTable
    } },
    retainedBytes
  });
};

// A saved layout owner (current or legacy ui fields) wins; without one, a
// canonical Session takes the layout projected from its committed request
// (`restoredLayoutPreferences`).
/** @param {DrawingState} drawing */
const restoreLayoutPreferences = (drawing, ui = {}, { projected = null } = {}) => {
  replaceLayoutPreferences(drawing.layoutPreferences, restoredLayoutPreferences(ui, {
    mode: state.mode.value,
    multiRecord: Boolean(drawing.form.multi_record_canvas),
    projected,
    active: { legend: drawing.form.legend, plotTitlePosition: drawing.adv.plot_title_position }
  }));
};

// Captured or stored settings (History Undo and Redo, the failed Session Load
// rollback, a settings-only Session) pass `resolveTrackPlacements: false`: the
// track stacks are installed as given, so unset slot sides, lane directions,
// and axis indexes stay unset (R11).
/** @param {DrawingState} drawing */
export const applyConfigData = (drawing, data, { resolveTrackPlacements = true } = {}) => {
  requireCurrentWebStateFieldNames(data);
  if (isPlainObject(data.form) && Object.prototype.hasOwnProperty.call(data.form, 'linear_track_layout')) {
    requireCurrentLinearTrackLayout(data.form.linear_track_layout);
  }
  if (isPlainObject(data.adv) && Object.prototype.hasOwnProperty.call(data.adv, 'label_placement')) {
    requireCurrentLinearLabelPlacement(data.adv.label_placement);
  }
  if (isPlainObject(data.adv) && Object.prototype.hasOwnProperty.call(data.adv, 'multi_record_size_mode')) {
    requireCurrentCircularMultiRecordSizeMode(data.adv.multi_record_size_mode);
  }
  if (data.form) safeDeepMerge(drawing.form, data.form);
  if (data.adv) {
    safeDeepMerge(drawing.adv, data.adv);
    ['scale_font_size', 'ruler_label_font_size'].forEach((field) => {
      if (Object.prototype.hasOwnProperty.call(data.adv, field)) {
        drawing.adv[field] = cloneJsonData(data.adv[field]);
      }
    });
  }
  replacePlainObject(
    drawing.unmanagedConfigOverrides,
    isPlainObject(data.unmanagedConfigOverrides)
      ? cloneJsonData(data.unmanagedConfigOverrides)
      : {}
  );
  drawing.recordDisplayDrafts.splice(0, drawing.recordDisplayDrafts.length, ...cloneJsonData(data.recordDisplayDrafts || []));
  replacePlainObject(drawing.featurePlacementOverrides, cloneJsonData(data.featurePlacementOverrides || {}));
  drawing.annotationSets.splice(
    0,
    drawing.annotationSets.length,
    ...normalizeAnnotationSets(data.annotationSets)
  );
  const linearLayout = data.linearRecordLayout && typeof data.linearRecordLayout === 'object'
    ? data.linearRecordLayout
    : null;
  // Omission takes the fresh default in every Session version.
  drawing.linearRecordLayoutEnabled.value = typeof linearLayout?.enabled === 'boolean'
    ? linearLayout.enabled
    : WEB_UX_PROFILE.linear.arrangeInRowsByDefault;
  const linearRecordGap = Number(linearLayout?.recordGap);
  drawing.linearRecordGap.value = Number.isFinite(linearRecordGap) && linearRecordGap >= 0
    ? linearRecordGap
    : 24;
  drawing.linearRecordRows.splice(
    0,
    drawing.linearRecordRows.length,
    ...(Array.isArray(linearLayout?.rows) ? linearLayout.rows : [])
      .map((entry) => ({
        uid: String(entry?.uid || ''), row: Number(entry?.row),
        ...(entry?.canonicalCardinality === 'exactly_one'
          ? { canonicalCardinality: 'exactly_one' } : {}),
        ...(Number(entry?.canonicalRow) === Number(entry?.row)
          && Number.isInteger(entry?.canonicalColumn) && entry.canonicalColumn > 0
          ? { canonicalRow: entry.canonicalRow, canonicalColumn: entry.canonicalColumn }
          : {})
      }))
      .filter((entry) => entry.uid && Number.isInteger(entry.row) && entry.row > 0)
  );
  // The record translations and alignment plan are the Linear mode's
  // artifact; only the Linear drawing's settings carry them.
  if (state.linearRecordTranslations && drawing === state.drawings.linear) {
    state.linearRecordTranslations.value = cloneJsonData(
      Array.isArray(linearLayout?.recordTranslations)
        ? linearLayout.recordTranslations
        : []
    );
  }
  if (state.similarityAlignmentPlan && drawing === state.drawings.linear) {
    state.similarityAlignmentPlan.value = linearLayout?.similarityAlignment
      ? cloneJsonData(linearLayout.similarityAlignment)
      : null;
  }
  replaceLinearComparisonPlan(
    drawing.linearComparisonPlan,
    data.linearComparisonPlan || createDefaultLinearComparisonPlan()
  );
  drawing.importedComparisonIntent.action = Object.values(IMPORTED_COMPARISON_ACTIONS)
    .includes(data.importedComparisonResolution?.action)
    ? data.importedComparisonResolution.action
    : null;
  drawing.adv.label_placement = requireCurrentLinearLabelPlacement(
    drawing.adv.label_placement
  );
  drawing.adv.label_rendering = normalizeLabelRendering(drawing.adv.label_rendering);
  drawing.adv.circular_label_placement =
    String(drawing.adv.circular_label_placement || '').trim().toLowerCase() === 'radial'
      ? 'radial'
      : 'horizontal';
  drawing.adv.circular_label_spacing = normalizePositiveNumberOrNull(drawing.adv.circular_label_spacing);
  drawing.adv.linear_label_spacing = normalizePositiveNumberOrNull(drawing.adv.linear_label_spacing);
  const rawTrackAxisGap = drawing.adv.track_axis_gap;
  if (
    rawTrackAxisGap === null ||
    rawTrackAxisGap === undefined ||
    rawTrackAxisGap === '' ||
    String(rawTrackAxisGap).trim().toLowerCase() === 'auto'
  ) {
    drawing.adv.track_axis_gap = null;
  } else {
    const numericTrackAxisGap = Number(rawTrackAxisGap);
    drawing.adv.track_axis_gap = Number.isFinite(numericTrackAxisGap) && numericTrackAxisGap >= 0
      ? numericTrackAxisGap
      : null;
  }
  drawing.form.linear_track_layout = requireCurrentLinearTrackLayout(
    drawing.form.linear_track_layout
  );
  drawing.form.plot_title = String(drawing.form.plot_title || '');
  // `form.legend` and `adv.plot_title_position` are accessors over
  // `layoutPreferences` (state.js) whose setters normalize a merged value.
  // Writing the resolved value back would only pin an unset Circular
  // multi-record preference, so neither is rewritten here.
  drawing.adv.feature_shapes = normalizeFeatureRenderingMap(drawing.adv.feature_shapes);
  Object.assign(drawing.adv, normalizedPersistedArrowGeometryState(data.adv));
  drawing.adv.multi_record_size_mode = requireCurrentCircularMultiRecordSizeMode(
    drawing.adv.multi_record_size_mode
  );
  const numericMinRadiusRatio = Number(drawing.adv.multi_record_min_radius_ratio);
  drawing.adv.multi_record_min_radius_ratio =
    Number.isFinite(numericMinRadiusRatio) && numericMinRadiusRatio > 0 && numericMinRadiusRatio <= 1
      ? numericMinRadiusRatio
      : 0.55;
  const numericColumnGapRatio = Number(drawing.adv.multi_record_column_gap_ratio);
  drawing.adv.multi_record_column_gap_ratio =
    Number.isFinite(numericColumnGapRatio) && numericColumnGapRatio >= 0
      ? numericColumnGapRatio
      : 0.10;
  const numericRowGapRatio = Number(drawing.adv.multi_record_row_gap_ratio);
  drawing.adv.multi_record_row_gap_ratio =
    Number.isFinite(numericRowGapRatio) && numericRowGapRatio >= 0
      ? numericRowGapRatio
      : 0.05;
  const rawMultiRecordPositions = Array.isArray(drawing.adv.multi_record_positions)
    ? drawing.adv.multi_record_positions
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
  drawing.adv.multi_record_positions = dedupedMultiRecordPositions
    .map((entry, index) => ({ ...entry, __index: index }))
    .sort((left, right) => {
      if (left.row !== right.row) return left.row - right.row;
      return left.__index - right.__index;
    })
    .map(({ __index, ...entry }) => entry);
  const rawPlotTitleFontSize = drawing.adv.plot_title_font_size;
  if (
    rawPlotTitleFontSize === null ||
    rawPlotTitleFontSize === undefined ||
    rawPlotTitleFontSize === ''
  ) {
    drawing.adv.plot_title_font_size = null;
  } else {
    const numericPlotTitleFontSize = Number(rawPlotTitleFontSize);
    drawing.adv.plot_title_font_size =
      Number.isFinite(numericPlotTitleFontSize) && numericPlotTitleFontSize > 0
        ? numericPlotTitleFontSize
        : null;
  }
  drawing.adv.keep_full_definition_with_plot_title =
    drawing.adv.keep_full_definition_with_plot_title === true;
  drawing.adv.depth_color = resolveColorToHex(String(drawing.adv.depth_color || '#4A90E2'));
  drawing.adv.depth_normalize = drawing.adv.depth_normalize === true;
  drawing.adv.depth_show_axis = drawing.adv.depth_show_axis !== false;
  drawing.adv.depth_show_ticks = drawing.adv.depth_show_ticks !== false;
  drawing.adv.depth_share_axis = drawing.adv.depth_share_axis === true;
  drawing.adv.depth_height = normalizePositiveNumberOrNull(drawing.adv.depth_height);
  drawing.adv.depth_width_circular = normalizePositiveNumberOrNull(drawing.adv.depth_width_circular);
  drawing.adv.circular_track_slots_schema_version = CIRCULAR_TRACK_SLOT_SCHEMA_VERSION;
  drawing.adv.circular_track_slots_enabled = drawing.adv.circular_track_slots_enabled === true;
  if (resolveTrackPlacements) {
    {
      const normalizedSlots = normalizeCircularTrackSlots(
        drawing.adv.circular_track_slots,
        drawing.adv.nt,
        drawing.form.track_type
      );
      const importedAxis = clampCircularTrackAxisIndex(
        drawing.adv.circular_track_slots_axis_index,
        normalizedSlots.length
      );
      drawing.adv.circular_track_slots_axis_index = importedAxis === null
        ? inferLegacyAxisIndexFromFeature(normalizedSlots, drawing.form.track_type)
        : importedAxis;
    }
    drawing.adv.circular_track_slots.splice(
      0,
      drawing.adv.circular_track_slots.length,
      ...applyCircularTrackOrderPlacements(
        drawing.adv.circular_track_slots,
        drawing.adv.nt,
        drawing.form.track_type,
        drawing.adv.circular_track_slots_axis_index
      )
    );
    drawing.adv.linear_track_slots_schema_version = LINEAR_TRACK_SLOT_SCHEMA_VERSION;
    drawing.adv.linear_track_slots_enabled = drawing.adv.linear_track_slots_enabled === true;
    {
      const normalizedLinearSlots = normalizeLinearTrackSlots(
        drawing.adv.linear_track_slots,
        drawing.adv.nt,
        drawing.form.linear_track_layout
      );
      drawing.adv.linear_track_slots_axis_index = clampLinearTrackAxisIndex(
        drawing.adv.linear_track_slots_axis_index,
        normalizedLinearSlots.length
      );
      drawing.adv.linear_track_slots_axis_index = resolveLinearTrackAxisIndex(
        normalizedLinearSlots,
        drawing.adv.linear_track_slots_axis_index
      );
      drawing.adv.linear_track_slots.splice(
        0,
        drawing.adv.linear_track_slots.length,
        ...applyLinearTrackOrderPlacements(
          normalizedLinearSlots,
          drawing.adv.linear_track_slots_axis_index,
          drawing.adv.nt,
          drawing.form.linear_track_layout
        )
      );
    }
  } else {
    // The merge installs the caller's slot objects; copy them so a later edit
    // never writes into a History entry or another caller-owned snapshot.
    ['circular_track_slots', 'linear_track_slots'].forEach((key) => {
      drawing.adv[key].splice(0, drawing.adv[key].length, ...cloneJsonData(drawing.adv[key]));
    });
  }
  drawing.adv.depth_window_size = normalizePositiveNumberOrNull(drawing.adv.depth_window_size);
  drawing.adv.depth_step_size = normalizePositiveNumberOrNull(drawing.adv.depth_step_size);
  drawing.adv.depth_large_tick_interval = normalizePositiveNumberOrNull(
    drawing.adv.depth_large_tick_interval
  );
  drawing.adv.depth_small_tick_interval = normalizePositiveNumberOrNull(drawing.adv.depth_small_tick_interval);
  drawing.adv.depth_tick_font_size = normalizePositiveNumberOrNull(drawing.adv.depth_tick_font_size);
  drawing.adv.depth_tracks.splice(
    0,
    drawing.adv.depth_tracks.length,
    ...normalizeDepthTracks(drawing.adv.depth_tracks, drawing.adv)
  );
  drawing.adv.gc_content_mode = String(drawing.adv.gc_content_mode || '').trim().toLowerCase() === 'percent'
    ? 'percent'
    : 'deviation';
  drawing.adv.gc_content_show_axis = drawing.adv.gc_content_show_axis !== false;
  drawing.adv.gc_content_show_ticks = drawing.adv.gc_content_show_ticks !== false;
  drawing.adv.gc_content_tick_interval = normalizePositiveNumberOrNull(drawing.adv.gc_content_tick_interval);
  drawing.adv.gc_content_small_tick_interval = normalizePositiveNumberOrNull(drawing.adv.gc_content_small_tick_interval);
  drawing.adv.gc_content_tick_font_size = normalizePositiveNumberOrNull(drawing.adv.gc_content_tick_font_size);
  const normalizeNonNegativeNumberOrNull = (value) => {
    if (
      value === null ||
      value === undefined ||
      value === '' ||
      String(value).trim().toLowerCase() === 'auto'
    ) {
      return null;
    }
    const numeric = Number(value);
    return Number.isFinite(numeric) && numeric >= 0 ? numeric : null;
  };
  drawing.adv.center_reserved_radius = normalizeNonNegativeNumberOrNull(drawing.adv.center_reserved_radius);
  drawing.adv.depth_min = normalizeNonNegativeNumberOrNull(drawing.adv.depth_min);
  drawing.adv.depth_max = normalizeNonNegativeNumberOrNull(drawing.adv.depth_max);
  if (
    drawing.adv.depth_min !== null &&
    drawing.adv.depth_max !== null &&
    drawing.adv.depth_min > drawing.adv.depth_max
  ) {
    drawing.adv.depth_max = null;
  }
  const normalizeFiniteNumberOrFallback = (value, fallback) => {
    if (value === null || value === undefined || value === '') return fallback;
    const numeric = Number(value);
    return Number.isFinite(numeric) ? numeric : fallback;
  };
  drawing.adv.gc_content_min_percent = normalizeFiniteNumberOrFallback(drawing.adv.gc_content_min_percent, 0);
  drawing.adv.gc_content_max_percent = normalizeFiniteNumberOrFallback(drawing.adv.gc_content_max_percent, 100);
  if (drawing.adv.gc_content_min_percent > drawing.adv.gc_content_max_percent) {
    drawing.adv.gc_content_max_percent = drawing.adv.gc_content_min_percent;
  }
  drawing.adv.linear_show_replicon = drawing.adv.linear_show_replicon === true;
  drawing.adv.linear_accession_visibility = requireLinearLabelVisibilityMode(
    drawing.adv.linear_accession_visibility,
    'Linear Accession visibility'
  );
  drawing.adv.linear_length_visibility = requireLinearLabelVisibilityMode(
    drawing.adv.linear_length_visibility,
    'Linear Length / Coordinates visibility'
  );
  drawing.adv.linear_definition_line_styles = normalizeDefinitionLineStyleState(
    drawing.adv.linear_definition_line_styles
  );
  drawing.adv.pairwise_match_style = normalizeCurrentPairwiseMatchStyle(
    drawing.adv.pairwise_match_style,
    'ribbon'
  );
  if (data.losat) {
    // The execution settings of a Session 44 or older draft or of a request
    // projection are app-level (`state.losatExecution`); the drawing keeps the
    // search settings.
    const executionFields = LOSAT_EXECUTION_FIELDS.filter((field) => Object.hasOwn(data.losat, field));
    if (executionFields.length) {
      applyLosatExecutionData(Object.fromEntries(executionFields.map((field) => [field, data.losat[field]])));
    }
    safeDeepMerge(drawing.losat, data.losat);
    drawing.losat.blastp.mode = normalizeBlastpMode(drawing.losat.blastp?.mode);
    drawing.losat.blastp.collinearInferOrthogroups = data.losat.blastp?.collinearInferOrthogroups ?? (drawing.losat.blastp.mode === 'collinear');
    drawing.losat.blastp.hitLimitsByMode = {
      ...createDefaultLosatpHitLimits(), ...cloneJsonData(data.losat.blastp?.hitLimitsByMode || {})
    };
    drawing.losat.blastp.maxHits = normalizePositiveInteger(drawing.losat.blastp?.maxHits, 5);
    drawing.losat.blastp.candidateLimit = normalizePositiveInteger(
      drawing.losat.blastp?.candidateLimit,
      null
    );
    if (
      drawing.losat.blastp.orthogroupMemberMaxHits === undefined &&
      drawing.losat.blastp.orthogroupMaxHits !== null &&
      drawing.losat.blastp.orthogroupMaxHits !== undefined
    ) {
      drawing.losat.blastp.orthogroupMemberMaxHits = drawing.losat.blastp.orthogroupMaxHits;
    }
    drawing.losat.blastp.orthogroupMembershipMode = normalizeOrthogroupMembershipMode(drawing.losat.blastp?.orthogroupMembershipMode);
    drawing.losat.blastp.orthogroupMemberMaxHits = normalizePositiveInteger(drawing.losat.blastp?.orthogroupMemberMaxHits, null);
    drawing.losat.blastp.collinearMinAnchors = normalizePositiveInteger(drawing.losat.blastp?.collinearMinAnchors, 1);
    {
      const maxGap = Number(drawing.losat.blastp?.collinearMaxUnitGap);
      drawing.losat.blastp.collinearMaxUnitGap = Number.isInteger(maxGap) && maxGap >= 0 ? maxGap : 0;
      const diagonalDrift = Number(drawing.losat.blastp?.collinearMaxDiagonalDrift);
      drawing.losat.blastp.collinearMaxDiagonalDrift = Number.isInteger(diagonalDrift) && diagonalDrift >= 0 ? diagonalDrift : 0;
      const mergeConflicts = Number(drawing.losat.blastp?.collinearMaxConflictsInMergeGap);
      drawing.losat.blastp.collinearMaxConflictsInMergeGap = Number.isInteger(mergeConflicts) && mergeConflicts >= 0 ? mergeConflicts : 1;
      const paralogLinks = Number(drawing.losat.blastp?.collinearMaxParalogLinksPerOrthogroup);
      drawing.losat.blastp.collinearMaxParalogLinksPerOrthogroup = Number.isInteger(paralogLinks) && paralogLinks > 0 ? paralogLinks : 2;
      drawing.losat.blastp.collinearColorMode = normalizeCollinearColorMode(drawing.losat.blastp?.collinearColorMode);
      const unitMode = String(drawing.losat.blastp?.collinearUnitMode || '').trim().toLowerCase();
      drawing.losat.blastp.collinearUnitMode = ['auto', 'cds', 'locus'].includes(unitMode) ? unitMode : 'auto';
      drawing.losat.blastp.collinearAnchorMode = normalizeCollinearAnchorMode(drawing.losat.blastp?.collinearAnchorMode);
      const mergeOrientation = String(
        drawing.losat.blastp?.collinearMergeOrientation || ''
      ).trim().toLowerCase();
      drawing.losat.blastp.collinearMergeOrientation = [
        'strand',
        'order',
        'either'
      ].includes(mergeOrientation) ? mergeOrientation : 'either';
      drawing.losat.blastp.collinearSearchScope = normalizeCollinearSearchScope(drawing.losat.blastp?.collinearSearchScope);
    }
    delete drawing.losat.blastp.collinearBlockMergeGap;
    delete drawing.losat.blastp.collinearSingletonMergeGap;
    delete drawing.losat.blastp.orthogroupHitPolicy;
    delete drawing.losat.blastp.orthogroupMaxHits;
  }
  if (typeof data.paletteInstantPreviewEnabled === 'boolean') {
    state.paletteInstantPreviewEnabled.value = data.paletteInstantPreviewEnabled;
  }
  const importedPalette = String(data.palette || '').trim();
  if (importedPalette) drawing.selectedPalette.value = importedPalette;
  if (hasColorEntries(data.colors)) {
    if (data.colorsAreOverrides) {
      const paletteColors = paletteColorsFromDefinitions(drawing.selectedPalette.value) || {};
      drawing.currentColors.value = normalizePaletteColors({
        ...paletteColors,
        ...normalizeColorMap(data.colors)
      });
    } else {
      drawing.currentColors.value = normalizePaletteColors(normalizeColorMap(data.colors));
    }
  } else {
    const paletteColors = paletteColorsFromDefinitions(drawing.selectedPalette.value);
    if (paletteColors) drawing.currentColors.value = paletteColors;
  }

  if (data.rules && Array.isArray(data.rules)) {
    drawing.manualSpecificRules.length = 0;
    data.rules.forEach((r) => {
      drawing.manualSpecificRules.push({
        feat: String(r.feat || ''),
        qual: String(r.qual || ''),
        val: String(r.val || ''),
        color: resolveColorToHex(String(r.color || '#000000')),
        cap: String(r.cap || ''),
        fromFile: !!r.fromFile
      });
    });
    drawing.fileLegendCaptions.value = new Set(
      drawing.manualSpecificRules
        .filter((rule) => rule.fromFile && rule.cap)
        .map((rule) => rule.cap)
    );
  }
  if (Object.prototype.hasOwnProperty.call(data, 'qualifierPriorityRules')) {
    replaceQualifierPriorityRules(drawing, data.qualifierPriorityRules);
  } else if (Object.prototype.hasOwnProperty.call(data, 'priorityRules')) {
    replaceQualifierPriorityRules(drawing, data.priorityRules);
  }
  if (data.filterMode) drawing.filterMode.value = data.filterMode;
  if (data.whitelist && Array.isArray(data.whitelist)) {
    drawing.manualWhitelist.length = 0;
    data.whitelist.forEach((w) => {
      drawing.manualWhitelist.push({
        feat: String(w.feat || ''),
        qual: String(w.qual || ''),
        key: String(w.key || '')
      });
    });
  }
  if (data.blacklistText !== undefined) drawing.manualBlacklist.value = String(data.blacklistText || '');
  if (data.losatProgram) {
    const program = String(data.losatProgram);
    drawing.losatProgram.value = ['blastn', 'tblastx', 'blastp'].includes(program) ? program : 'blastn';
  }
  if (data.circularConservation) {
    safeDeepMerge(drawing.circularConservation, data.circularConservation);
  }
  drawing.circularConservation.enabled = drawing.circularConservation.enabled === true;
  drawing.circularConservation.source = normalizeCircularConservationSource(drawing.circularConservation.source);
  drawing.circularConservation.losat_program = normalizeCircularConservationLosatProgram(
    drawing.circularConservation.losat_program
  );
  drawing.circularConservation.subject_gencode = normalizePositiveInteger(drawing.circularConservation.subject_gencode, 1);
  drawing.circularConservation.reference = normalizeCircularConservationReference(drawing.circularConservation.reference);
  drawing.circularConservation.labels = String(drawing.circularConservation.labels || '');
  drawing.circularConservation.series.splice(
    0,
    drawing.circularConservation.series.length,
    ...normalizeCircularConservationSeries(drawing.circularConservation.series)
  );
  drawing.circularConservation.ring_width = normalizePositiveNumberOrNull(drawing.circularConservation.ring_width);
  drawing.circularConservation.ring_gap = normalizePositiveNumberOrNull(drawing.circularConservation.ring_gap);
  const webEdits = data.webEdits && typeof data.webEdits === 'object' ? data.webEdits : {};
  if (Object.prototype.hasOwnProperty.call(webEdits, 'orthogroupNameOverrides')) {
    replaceStringMap(drawing.orthogroupNameOverrides, webEdits.orthogroupNameOverrides);
  }
  if (Object.prototype.hasOwnProperty.call(webEdits, 'orthogroupDescriptionOverrides')) {
    replaceStringMap(drawing.orthogroupDescriptionOverrides, webEdits.orthogroupDescriptionOverrides);
  }
  // Absent in older Sessions: no dormant names (D-21).
  clearObject(drawing.orthogroupDormantOverrides);
  Object.assign(drawing.orthogroupDormantOverrides, normalizeOrthogroupDormantOverrides(webEdits.orthogroupDormantOverrides));
};

/** @param {DrawingState} drawing */
const restorePaletteStateAfterConfigImport = (drawing) => {
  const draftPaletteName = String(drawing.selectedPalette.value || state.appliedPaletteName.value || 'default');
  const draftColors = normalizePaletteColors(cloneColors(drawing.currentColors.value));
  const hasPreviewResults = Array.isArray(state.results.value) && state.results.value.length > 0;

  if (
    !hasPreviewResults ||
    state.paletteInstantPreviewEnabled.value ||
    draftPaletteName === String(state.appliedPaletteName.value || '')
  ) {
    state.appliedPaletteName.value = draftPaletteName;
    state.appliedPaletteColors.value = draftColors;
    drawing.pendingPaletteName.value = '';
    drawing.pendingPaletteColors.value = {};
    return;
  }

  drawing.pendingPaletteName.value = draftPaletteName;
  drawing.pendingPaletteColors.value = draftColors;
};

// The palette the shown Result was drawn with (an artifact value), or else
// the drawing's palette.
/** @param {DrawingState} drawing */
const restoreAppliedPaletteFromSession = (drawing, ui = {}) => {
  const draftPaletteName = String(drawing.selectedPalette.value || state.appliedPaletteName.value || 'default');
  const draftColors = normalizePaletteColors(cloneColors(drawing.currentColors.value));
  const savedAppliedPaletteName = String(ui.appliedPaletteName || draftPaletteName || 'default');
  const savedAppliedPaletteColors =
    ui.appliedPaletteColors && typeof ui.appliedPaletteColors === 'object'
      ? Object.fromEntries(
          Object.entries(ui.appliedPaletteColors).map(([key, value]) => [
            key,
            resolveColorToHex(String(value || '').trim())
          ])
        )
      : draftColors;
  state.appliedPaletteName.value = savedAppliedPaletteName;
  state.appliedPaletteColors.value = normalizePaletteColors(cloneColors(savedAppliedPaletteColors));
};

// A drawing's pending palette applies only while Instant Preview is off.
/** @param {DrawingState} drawing */
const restorePendingPaletteFromSession = (drawing, ui = {}) => {
  const draftColors = normalizePaletteColors(cloneColors(drawing.currentColors.value));
  const savedPendingPaletteName = String(ui.pendingPaletteName || '').trim();
  const savedPendingPaletteColors =
    ui.pendingPaletteColors && typeof ui.pendingPaletteColors === 'object'
      ? Object.fromEntries(
          Object.entries(ui.pendingPaletteColors).map(([key, value]) => [
            key,
            resolveColorToHex(String(value || '').trim())
          ])
        )
      : draftColors;
  if (!state.paletteInstantPreviewEnabled.value && savedPendingPaletteName) {
    drawing.pendingPaletteName.value = savedPendingPaletteName;
    drawing.pendingPaletteColors.value = normalizePaletteColors(cloneColors(savedPendingPaletteColors));
  } else {
    drawing.pendingPaletteName.value = '';
    drawing.pendingPaletteColors.value = {};
  }
};

/** @param {DrawingState} drawing */
const restorePaletteStateFromSession = (drawing, ui = {}) => {
  restoreAppliedPaletteFromSession(drawing, ui);
  restorePendingPaletteFromSession(drawing, ui);
};

const serializedFileDescriptors = new WeakMap();
const cacheSerializedFileDescriptor = (file, descriptor) => {
  setResourcePayloadOwner(descriptor, file);
  serializedFileDescriptors.set(file, descriptor);
  return descriptor;
};

const serializeFile = async (file) => {
  if (!file) return null;
  const source = /** @type {SessionResourceSource | null} */ (getSessionResourceSource(file));
  if (source?.descriptor) {
    const visibleName = String(file.name || '').trim();
    const descriptor = visibleName && visibleName !== source.descriptor.name
      ? { ...source.descriptor, name: visibleName }
      : source.descriptor;
    setResourcePayloadOwner(descriptor, file);
    return descriptor;
  }
  const cached = serializedFileDescriptors.get(file);
  if (cached) return cached;
  const bytes = await readFileBytes(file);
  const shared = serializedFileDescriptors.get(file);
  if (shared) return shared;
  recordStructuralMetric('resourceReencodeCount', 1, {
    resourceName: String(file.name || 'file')
  });
  recordStructuralMetric('base64EncodeCount', 1, {
    resourceName: String(file.name || 'file')
  });
  recordStructuralMetric('encodedByteCount', bytes.byteLength, {
    resourceName: String(file.name || 'file')
  });
  const descriptor = {
    name: file.name || 'file',
    type: file.type || '',
    size: bytes.byteLength,
    lastModified: file.lastModified ?? Date.now(),
    encoding: 'base64',
    data: bytesToBase64(bytes)
  };
  return cacheSerializedFileDescriptor(file, descriptor);
};

const serializeFileValue = async (value) => (
  Array.isArray(value)
    ? Promise.all(value.map((item) => serializeFileValue(item)))
    : serializeFile(value)
);

const deserializeFile = (entry) => {
  if (Array.isArray(entry)) {
    return entry.map((item) => deserializeFile(item));
  }
  if (isSessionResourceFileView(entry)) return entry;
  if (
    !entry ||
    !Object.prototype.hasOwnProperty.call(entry, 'data') ||
    entry.data === null ||
    entry.data === undefined
  ) return null;
  if (isEncodedDepthFileEntry(entry)) {
    const text = decodeDepthText(entry.data);
    recordStructuralMetric('fileConstructionCount');
    return new File([text], entry.name || 'depth.tsv', {
      type: entry.type || 'text/tab-separated-values',
      lastModified: entry.lastModified ?? Date.now()
    });
  }
  if (typeof entry.data !== 'string') return null;
  const bytes = base64ToBytes(entry.data);
  recordStructuralMetric('base64DecodeCount');
  recordStructuralMetric('decodedByteCount', bytes.byteLength);
  recordStructuralMetric('fileConstructionCount');
  const file = new File([bytes], entry.name || 'file', {
    type: entry.type || 'application/octet-stream',
    lastModified: entry.lastModified ?? Date.now()
  });
  cacheSerializedFileDescriptor(file, entry);
  return file;
};

/** @type {{ clearActiveRuntime?: () => void } | null} */
let activePreviewRuntime = null;

export const setPreviewRuntime = (runtime) => {
  activePreviewRuntime = runtime || null;
};

export const serializeResults = () => {
  if (activePreviewRuntime) {
    return normalizeLogicalResults(state.results.value.map((res, idx) => ({
      name: res.name || `Result ${idx + 1}`,
      content: res.content
    })));
  }

  const currentSvg = (() => {
    if (!state.svgContainer.value) return null;
    const svg = state.svgContainer.value.querySelector('svg');
    if (!svg) return null;
    return serializeCleanSvg(svg);
  })();

  return normalizeLogicalResults(state.results.value.map((res, idx) => ({
    name: res.name || `Result ${idx + 1}`,
    content: idx === state.selectedResultIndex.value && currentSvg ? currentSvg : res.content
  })));
};

const LOSAT_CACHE_INFO_STRING_FIELDS = ['edgeKey', 'queryUid', 'subjectUid'];
const LOSAT_CACHE_INFO_INTEGER_FIELDS = ['ordinal', 'queryIndex', 'subjectIndex'];
const losatCacheInfoIdentity = (entry) => {
  const identity = {};
  LOSAT_CACHE_INFO_STRING_FIELDS.forEach((field) => {
    if (typeof entry?.[field] === 'string' && entry[field].trim()) {
      identity[field] = entry[field];
    }
  });
  LOSAT_CACHE_INFO_INTEGER_FIELDS.forEach((field) => {
    if (Number.isInteger(entry?.[field])) identity[field] = entry[field];
  });
  return identity;
};

/** @param {DrawingState} drawing */
const restoredLosatCacheInfoIdentity = (drawing, entry) => {
  const identity = losatCacheInfoIdentity(entry);
  if (state.mode.value !== 'linear') return identity;

  const queryInstanceUid = String(entry?.queryRecordInstanceKey || '').trim();
  const subjectInstanceUid = String(entry?.subjectRecordInstanceKey || '').trim();
  if (
    Boolean(queryInstanceUid) !== Boolean(subjectInstanceUid)
    || (queryInstanceUid && identity.queryUid && queryInstanceUid !== identity.queryUid)
    || (subjectInstanceUid && identity.subjectUid && subjectInstanceUid !== identity.subjectUid)
  ) return {};

  const queryUid = queryInstanceUid || identity.queryUid || '';
  const subjectUid = subjectInstanceUid || identity.subjectUid || '';
  if (!queryUid || !subjectUid || queryUid === subjectUid) return {};

  const indexByUid = new Map(
    state.linearSeqs.map((sequence, index) => [String(sequence?.uid || ''), index])
  );
  const queryIndex = indexByUid.get(queryUid);
  const subjectIndex = indexByUid.get(subjectUid);
  if (!Number.isInteger(queryIndex) || !Number.isInteger(subjectIndex)) return {};

  const edgeKey = linearComparisonEdgeKey(queryUid, subjectUid);
  if (identity.edgeKey && identity.edgeKey !== edgeKey) return {};
  const resolved = drawing.linearComparisonResolution.value.edges.find(
    (edge) => edge.edgeKey === edgeKey
  );
  return {
    ...identity,
    edgeKey,
    queryUid,
    subjectUid,
    queryIndex,
    subjectIndex,
    ...(Number.isInteger(resolved?.ordinal) ? { ordinal: resolved.ordinal } : {})
  };
};

const serializeLosatCache = () => {
  const cacheMap = state.losatCache?.value;
  if (!cacheMap || cacheMap.size === 0) {
    return { entries: [], validatedManifest: null, manifestValidated: false };
  }
  const manifest = state.proteinIdentityManifest.value;
  const rawManifest = rawReactiveValue(manifest);
  const info = Array.isArray(state.losatCacheInfo.value) ? state.losatCacheInfo.value : [];
  const entries = [];
  const seen = new Set();
  const prevalidatedProteinEntries = new WeakSet();

  const buildEntry = (key, cached, infoEntry = {}) => {
    const rawCached = rawReactiveValue(cached);
    let serialized = rawCached;
    if (!adoptedLosatCacheValues.has(rawCached)) {
      const { text, ...metadata } = rawCached;
      serialized = {
        ...cloneJsonData(metadata),
        text: String(text ?? '')
      };
    }
    const entry = {
      ...serialized,
      key: String(key),
      filename: String(infoEntry.filename || ''),
      display: Boolean(infoEntry.display),
      ...losatCacheInfoIdentity(infoEntry)
    };
    if (
      classifyRawLosatCacheEntry(entry) === 'protein-current'
      && adoptedLosatCacheValues.get(rawCached) === rawManifest
    ) {
      prevalidatedProteinEntries.add(entry);
    }
    return entry;
  };

  // One entry per raw key: rows that share a key (two rings of one sequence)
  // keep the first row's filename, as the CLI raw cache does.
  info.forEach((entry, idx) => {
    if (!entry || !entry.key || seen.has(entry.key)) return;
    const cached = cacheMap.get(entry.key);
    if (!isCurrentRawLosatCacheEntry(cached)) return;
    entries.push(buildEntry(entry.key, cached, {
      ...entry,
      filename: entry.filename || `losat_pair_${idx + 1}.tsv`,
      display: entry.display !== false
    }));
    seen.add(entry.key);
  });

  cacheMap.forEach((value, key) => {
    if (seen.has(key)) return;
    if (!isCurrentRawLosatCacheEntry(value)) return;
    entries.push(buildEntry(key, value));
  });

  // A replaced source can leave earlier bindings in the live cache. Persist
  // only protein evidence that the Session's current manifest can resolve.
  /** @type {ReturnType<typeof buildValidatedProteinIdentityIndex>} */
  let identityIndex = null;
  let reusedValidation = false;
  try {
    return {
      entries: entries.filter((entry) => {
        if (classifyRawLosatCacheEntry(entry) !== 'protein-current') return true;
        if (prevalidatedProteinEntries.has(entry)) {
          reusedValidation = true;
          recordStructuralMetric('sessionSaveProteinRawTextValidationReuseCount');
          return true;
        }
        if (!identityIndex) {
          identityIndex = buildValidatedProteinIdentityIndex(manifest);
        }
        if (!identityIndex) {
          throw new Error('Save Session requires a valid protein identity manifest.');
        }
        recordStructuralMetric('sessionSaveProteinRawTextValidationCount');
        return validateProteinRawEntryReferences(entry, manifest, { identityIndex });
      }),
      validatedManifest: identityIndex || reusedValidation ? manifest : null,
      manifestValidated: Boolean(identityIndex || reusedValidation)
    };
  } finally {
    if (identityIndex) releaseValidatedProteinIdentityIndex(identityIndex);
  }
};

/** @param {DrawingState} drawing */
const applyLosatCache = (
  drawing,
  entries,
  legacyEnvelope = null,
  { adoptCurrent = false, validatedManifest = null } = {}
) => {
  const map = new Map();
  const info = [];
  const legacyEntries = [];

  if (Array.isArray(entries)) {
    entries.forEach((entry, idx) => {
      const classification = classifyRawLosatCacheEntry(entry);
      if (classification === 'protein-legacy') {
        legacyEntries.push(entry);
        return;
      }
      if (
        !entry ||
        !isCurrentRawLosatCacheEntry(entry) ||
        !entry.key
      ) {
        return;
      }
      const excludedFields = new Set([
        'key',
        'filename',
        'display',
        ...LOSAT_CACHE_INFO_STRING_FIELDS,
        ...LOSAT_CACHE_INFO_INTEGER_FIELDS
      ]);
      const restored = adoptCurrent
        ? Object.fromEntries(
            Object.entries(entry).filter(([field]) => !excludedFields.has(field))
          )
        : cloneJsonData(entry);
      if (!adoptCurrent) {
        excludedFields.forEach((field) => delete restored[field]);
      } else {
        adoptedLosatCacheValues.set(restored, rawReactiveValue(validatedManifest));
      }
      map.set(entry.key, restored);
      if (entry.display === false) return;
      info.push({
        key: entry.key,
        filename: entry.filename || `losat_pair_${idx + 1}.tsv`,
        display: true,
        ...restoredLosatCacheInfoIdentity(drawing, entry)
      });
    });
  }

  state.losatCache.value = map;
  state.losatCacheInfo.value = info;
  const restoredEnvelope = normalizeLegacyProteinCandidateEnvelope(legacyEnvelope);
  const importedEnvelope = createLegacyProteinCandidateEnvelope(legacyEntries);
  state.legacyProteinRawCandidates.value = {
    schema: 1,
    entries: [...restoredEnvelope.entries, ...importedEnvelope.entries]
  };
};

const pruneLosatDerivedCache = (map) => {
  if (!map || typeof map.delete !== 'function') return;
  while (map.size > LOSAT_DERIVED_CACHE_LIMIT) {
    const oldestKey = map.keys().next().value;
    if (oldestKey === undefined) break;
    map.delete(oldestKey);
  }
};

const normalizeLegacyDerivedEvidence = (value, fallbackEntries = []) => ({
  schema: 1,
  entries: [
    ...(
      isPlainObject(value) && value.schema === 1 && Array.isArray(value.entries)
        ? value.entries
        : []
    ),
    ...fallbackEntries
  ]
    .filter((entry) => isLosatDerivedCacheEntry(entry) && entry.schema === 1)
    .map((entry) => cloneJsonData(entry))
});

const applyLosatDerivedCache = (
  entries,
  legacyEvidence = null,
  { adoptCurrent = false } = {}
) => {
  const map = new Map();
  const legacyEntries = [];

  if (Array.isArray(entries)) {
    entries.forEach((entry) => {
      if (!isLosatDerivedCacheEntry(entry)) return;
      if (entry.schema === 1) {
        legacyEntries.push(entry);
        return;
      }
      map.set(entry.key, {
        schema: LOSAT_DERIVED_CACHE_SCHEMA,
        kind: 'derived-losatp-payload',
        idEncoding: 'runtime-handle-v1',
        key: entry.key,
        mode: String(entry.mode || ''),
        payload: adoptCurrent ? entry.payload : cloneJsonData(entry.payload)
      });
    });
  }

  pruneLosatDerivedCache(map);
  state.losatDerivedCache.value = map;
  state.legacyProteinDerivedEvidence.value = normalizeLegacyDerivedEvidence(
    legacyEvidence,
    legacyEntries
  );
};

const applyProteinIdentityManifest = (manifest, { adoptCurrent = false } = {}) => {
  if (!validateProteinIdentityManifest(manifest)) {
    state.proteinIdentityManifest.value = emptyProteinIdentityManifest();
    adoptedProteinIdentityManifest = null;
    return;
  }
  state.proteinIdentityManifest.value = adoptCurrent
    ? manifest
    : cloneJsonData(manifest);
  adoptedProteinIdentityManifest = adoptCurrent ? manifest : null;
};

/**
 * @param {DrawingState} drawing
 * @param {Record<string, any>} [orthogroupState]
 * @param {{ legacyRecords?: any, catalogFeatureState?: Record<string, any> | null }} [options]
 */
export const applyOrthogroupStateData = (
  drawing, orthogroupState = {}, { legacyRecords = null, catalogFeatureState = null } = {}
) => {
  const storedGroups = Array.isArray(orthogroupState.groups) ? orthogroupState.groups : [];
  const groups = legacyRecords
    ? migrateLegacyOrthogroupMembers(storedGroups, legacyRecords)
    : storedGroups;
  const groupIds = groups
    .map((group) => String(group?.id || '').trim())
    .filter(Boolean);
  const groupIdSet = new Set(groupIds);
  const index = catalogFeatureState?.featureOrthogroupIndex || buildOrthogroupFeatureIndex(groups);

  state.orthogroups.value = groups;
  state.featureOrthogroupIndex.value = index;
  // Current catalog admission has already projected these exact groups onto
  // its features. Reuse that owner rather than deriving the same metadata twice.
  if (!catalogFeatureState) {
    state.extractedFeatures.value = enrichFeaturesWithOrthogroups(
      state.extractedFeatures.value,
      index
    );
    if (state.biologicalFeatures) {
      state.biologicalFeatures.value = enrichFeaturesWithOrthogroups(
        state.biologicalFeatures.value,
        index
      );
    }
  }
  const selectedId = String(orthogroupState.selectedOrthogroupId || '').trim();
  state.selectedOrthogroupId.value = selectedId && groupIdSet.has(selectedId) ? selectedId : (groupIds[0] || '');
  state.selectedOrthogroupAlignmentFeature.value = String(orthogroupState.selectedOrthogroupAlignmentFeature || '').trim();

  replaceStringMap(drawing.orthogroupNameOverrides, orthogroupState.orthogroupNameOverrides);
  replaceStringMap(drawing.orthogroupDescriptionOverrides, orthogroupState.orthogroupDescriptionOverrides);
  Object.keys(drawing.orthogroupNameOverrides).forEach((id) => {
    if (!groupIdSet.has(id)) delete drawing.orthogroupNameOverrides[id];
  });
  Object.keys(drawing.orthogroupDescriptionOverrides).forEach((id) => {
    if (!groupIdSet.has(id)) delete drawing.orthogroupDescriptionOverrides[id];
  });
  if (Object.hasOwn(orthogroupState, 'orthogroupDormantOverrides')) {
    clearObject(drawing.orthogroupDormantOverrides);
    Object.assign(drawing.orthogroupDormantOverrides,
      normalizeOrthogroupDormantOverrides(orthogroupState.orthogroupDormantOverrides));
  }
};

/** @param {string} mode @param {DrawingState} drawing */
const customDepthRequested = (mode, drawing) => {
  const adv = drawing?.adv || {};
  const form = drawing?.form || {};
  const customEnabled = mode === 'linear'
    ? Boolean(adv.linear_track_slots_enabled)
    : Boolean(adv.circular_track_slots_enabled);
  if (!customEnabled) return Boolean(form.show_depth);
  const slots = mode === 'linear'
    ? adv.linear_track_slots
    : adv.circular_track_slots;
  return (Array.isArray(slots) ? slots : []).some(
    (slot) => slot?.enabled !== false && slot?.renderer === 'depth'
  );
};

/**
 * The one check that the active mode has its biological inputs; Generate and
 * Save from the draft both use it, so both explain the same missing input.
 */
export const assertActiveModeInputs = (mode = state.mode.value, sourceState = state) => {
  const files = sourceState.files || {};
  if (mode === 'circular') {
    const gff = sourceState.cInputType?.value === 'gff';
    if (!(gff ? files.c_gff : files.c_gb)) throw diagnosticError('INPUT_REQUIRED');
    if (gff && !files.c_fasta) throw diagnosticError('FASTA_REQUIRED');
    return;
  }
  const gff = sourceState.lInputType?.value === 'gff';
  (Array.isArray(sourceState.linearSeqs) ? sourceState.linearSeqs : []).forEach((seq, index) => {
    if (!(gff ? seq?.gff : seq?.gb)) throw diagnosticError('INPUT_REQUIRED', { inputOrdinal: index + 1 });
    if (gff && !seq?.fasta) throw diagnosticError('FASTA_REQUIRED', { inputOrdinal: index + 1 });
  });
};

export const materializeLinearRecordFiles = (
  sequences,
  catalog,
  _options = {}
) => {
  const sourceSequences = Array.isArray(sequences) ? sequences : [];
  if (catalog == null) return sourceSequences;
  if (catalog?.mode !== 'linear' || catalog?.status !== 'ready') {
    const issue = Array.isArray(catalog?.issues) ? catalog.issues[0] : null;
    throw issue ? diagnosticError(issue.code, issue.context)
      : diagnosticError('RECORD_SELECTION', { reason: 'DISCOVERY_PENDING' });
  }
  const records = Array.isArray(catalog.records) ? catalog.records : [];
  if (records.length === 0) throw diagnosticError('NO_RECORDS');
  const recordCountBySource = new Map();
  records.forEach((record) => {
    const sourceIndex = Number(record?.sourceIndex);
    recordCountBySource.set(sourceIndex, (recordCountBySource.get(sourceIndex) || 0) + 1);
  });
  sourceSequences.forEach((source, sourceIndex) => {
    const count = recordCountBySource.get(sourceIndex) || 0;
    if (count === 0) throw diagnosticError('NO_RECORDS', { inputOrdinal: sourceIndex + 1 });
    if (count <= 1) return;
    const hasRegion = [source.region_start, source.region_end].some(
      (value) => value !== null && value !== undefined && value !== ''
    );
    if (hasRegion) {
      throw diagnosticError('REGION_INVALID', { inputOrdinal: sourceIndex + 1, reason: 'SELECT_RECORD_FOR_REGION' });
    }
  });
  return sourceSequences;
};

/**
 * @param {string} mode
 * @param {Record<string, any>} sourceState The project inputs (`state`).
 * @param {DrawingState} drawing The drawing whose settings choose the inputs.
 * @param {Record<string, any> | null} [comparisonPlanOrOptions]
 */
export const serializeActiveRenderFiles = async (
  mode,
  sourceState,
  drawing,
  comparisonPlanOrOptions = null
) => {
  if (!['circular', 'linear'].includes(mode)) {
    throw new Error(`Unsupported render mode: ${String(mode)}.`);
  }
  const sourceFiles = sourceState.files || {};
  const normalizedLinearSeqs = mode === 'linear'
    ? normalizeLinearSeqList(sourceState.linearSeqs)
    : [];
  const depthRequested = customDepthRequested(mode, drawing);
  const serializedLinearSeqs = await Promise.all(
    normalizedLinearSeqs.map(async (seq) => ({
      uid: seq.uid,
      gb: await serializeFile(seq.gb),
      gff: await serializeFile(seq.gff),
      fasta: await serializeFile(seq.fasta),
      depth: depthRequested ? await serializeFileValue(seq.depth) : null,
      losat_gencode: seq.losat_gencode ?? 1,
      definition: seq.definition ?? '',
      record_subtitle: seq.record_subtitle ?? '',
      file_definition: seq.file_definition ?? '',
      file_subtitle: seq.file_subtitle ?? '',
      inferred_definition: seq.inferred_definition ?? '',
      region_record_id: seq.region_record_id ?? '',
      region_start: seq.region_start ?? null,
      region_end: seq.region_end ?? null,
      region_reverse: !!seq.region_reverse
    }))
  );
  const optionBag = comparisonPlanOrOptions
    && typeof comparisonPlanOrOptions === 'object'
    && (
      Object.prototype.hasOwnProperty.call(comparisonPlanOrOptions, 'comparisonPlan') ||
      Object.prototype.hasOwnProperty.call(comparisonPlanOrOptions, 'linearRecordCatalog')
    )
    ? comparisonPlanOrOptions
    : null;
  const suppliedComparisonPlan = optionBag
    ? optionBag.comparisonPlan
    : comparisonPlanOrOptions;
  const linearSeqs = materializeLinearRecordFiles(
    serializedLinearSeqs,
    optionBag?.linearRecordCatalog ?? null,
    { layoutEnabled: Boolean(drawing.linearRecordLayoutEnabled?.value) }
  );
  const resolvedComparisonPlan = mode === 'linear'
    ? suppliedComparisonPlan || resolveLinearComparisonPlan({
        plan: drawing.linearComparisonPlan,
        sequences: normalizedLinearSeqs,
        layout: drawing.linearRecordLayoutEnabled?.value
          ? drawing.linearRecordRows
          : [],
        losatProgram: drawing.losatProgram?.value,
        blastpMode: drawing.losat?.blastp?.mode
      })
    : null;
  const linearComparisons = resolvedComparisonPlan ? await Promise.all(
    resolvedComparisonPlan.edges
      .filter((edge) => edge.source === 'upload' && edge.fileActive && edge.file)
      .map(async (edge) => ({
      id: String(edge.id || ''),
      file: await serializeFile(edge.file)
    }))
  ) : [];
  const linearCanonicalComparisons = mode === 'linear'
    && resolvedComparisonPlan?.mode !== 'none'
    && resolvedComparisonPlan?.edges?.length > 0
    ? await Promise.all(
    (Array.isArray(sourceFiles.linearCanonicalComparisons)
      ? sourceFiles.linearCanonicalComparisons
      : []
    ).map(async (comparison) => (
      isResourceBackedCanonicalComparison(comparison)
        ? {
            ...mapResourceBackedCanonicalComparison(comparison),
            file: await serializeFile(comparison.file)
          }
        : cloneJsonData(comparison)
    ))
    )
    : [];
  const conservationEnabled = mode === 'circular'
    && Boolean(drawing.circularConservation?.enabled);
  const conservationSource = String(drawing.circularConservation?.source || 'upload');
  const includeConservationBlasts = conservationEnabled && (
    conservationSource === 'upload'
    || sourceFiles.c_conservation_blasts_source === 'losat-cache'
  );
  const includeConservationFastas =
    conservationEnabled && conservationSource === 'losat';

  return {
    c_gb: mode === 'circular' ? await serializeFile(sourceFiles.c_gb) : null,
    c_gff: mode === 'circular' ? await serializeFile(sourceFiles.c_gff) : null,
    c_fasta: mode === 'circular' ? await serializeFile(sourceFiles.c_fasta) : null,
    c_depth: mode === 'circular' && depthRequested
      ? await serializeFileValue(sourceFiles.c_depth)
      : null,
    c_conservation_blasts: includeConservationBlasts
      ? await serializeFileValue(sourceFiles.c_conservation_blasts)
      : [],
    c_conservation_blasts_source: includeConservationBlasts
      && sourceFiles.c_conservation_blasts_source === 'losat-cache'
      ? 'losat-cache'
      : null,
    c_conservation_fastas: includeConservationFastas
      ? await serializeFileValue(sourceFiles.c_conservation_fastas)
      : [],
    c_conservation_sequence_sources: includeConservationBlasts
      ? await serializeFileValue(sourceFiles.c_conservation_sequence_sources)
      : [],
    d_color: null,
    t_color: null,
    blacklist: null,
    whitelist: null,
    qualifier_priority: null,
    linearSeqs,
    linearComparisons,
    linearCanonicalComparisons
  };
};

const deserializeCanonicalComparisons = (
  comparisons,
  { adoptCanonicalPayloads = false } = {}
) => (
  Array.isArray(comparisons)
    ? comparisons.map((comparison) => (
        isResourceBackedCanonicalComparison(comparison)
          ? mapResourceBackedCanonicalComparison(comparison, deserializeFile)
          : adoptCanonicalPayloads ? comparison : cloneJsonData(comparison)
      ))
    : []
);

export const adoptCanonicalRenderArtifacts = (
  canonical,
  { adoptOwnedRequest = false } = {}
) => {
  const ownedCanonical = adoptOwnedRequest
    ? adoptRuntimeCanonicalSession(canonical)
    : canonical;
  const preserveAdoptedResources = isAdoptedCanonicalSession(ownedCanonical);
  const sessionResourceTable = preserveAdoptedResources
    ? adoptCurrentSessionResources(ownedCanonical?.resources)
    : null;
  const projection = projectCanonicalSessionRequest({
    ...ownedCanonical,
    sessionResourceTable,
    deferResourceContent: preserveAdoptedResources,
    adoptCanonicalPayloads: preserveAdoptedResources
  });
  const nextCommittedCanonicalSession = preserveAdoptedResources
    ? adoptRuntimeCanonicalSession(ownedCanonical)
    : cloneCanonicalSession(ownedCanonical);
  /** @type {ReturnType<typeof deserializeCanonicalComparisons> | null} */
  let nextLinearComparisons = null;
  // The Circular conservation settings belong to the Circular drawing.
  const drawing = state.drawings.circular;
  /** @type {{
   *   blasts: any[], fastas: any[], sequenceSources: any[],
   *   projectedConservation: Record<string, any> | undefined, series: Record<string, any>[]
   * } | null} */
  let nextCircularState = null;
  if (projection.mode === 'linear') {
    nextLinearComparisons = deserializeCanonicalComparisons(
      projection.files.linearCanonicalComparisons,
      { adoptCanonicalPayloads: preserveAdoptedResources }
    );
  } else if (projection.files.c_conservation_blasts_source === 'losat-cache') {
    const projectedConservation = projection.config.circularConservation;
    const currentSeries = Array.isArray(drawing.circularConservation.series)
      ? drawing.circularConservation.series.map((entry) => cloneJsonData(entry))
      : [];
    const nextSeries = Array.isArray(projectedConservation?.series)
      ? projectedConservation.series.map((entry, index) => ({
          ...entry,
          ...(currentSeries[index] || {}),
          sourceIndex: index,
          fileName: entry.fileName,
          label: entry.label,
          color: entry.color
        }))
      : [];
    nextCircularState = {
      blasts: (projection.files.c_conservation_blasts || [])
        .map((entry) => deserializeFile(entry))
        .filter(Boolean),
      fastas: (projection.files.c_conservation_fastas || [])
        .map((entry) => deserializeFile(entry)),
      sequenceSources: (projection.files.c_conservation_sequence_sources || [])
        .map((entry) => deserializeFile(entry)),
      projectedConservation,
      series: nextSeries
    };
  }

  // Everything that can validate, project, or deserialize finishes before the
  // committed request and its comparison-backed draft resources are replaced.
  activeSessionResourceTable = sessionResourceTable;
  committedCanonicalSession = nextCommittedCanonicalSession;
  if (nextLinearComparisons) {
    state.files.linearCanonicalComparisons = nextLinearComparisons;
  }
  if (nextCircularState) {
    state.files.c_conservation_blasts = nextCircularState.blasts;
    state.files.c_conservation_blasts_source = 'losat-cache';
    state.files.c_conservation_fastas = nextCircularState.fastas;
    state.files.c_conservation_sequence_sources = nextCircularState.sequenceSources;
    if (nextCircularState.projectedConservation) {
      drawing.circularConservation.enabled = true;
      drawing.circularConservation.source = 'losat';
      drawing.circularConservation.reference = nextCircularState.projectedConservation.reference;
      drawing.circularConservation.labels = nextCircularState.series
        .map((entry) => entry.label).join(',');
      drawing.circularConservation.ring_width = nextCircularState.projectedConservation.ring_width;
      drawing.circularConservation.ring_gap = nextCircularState.projectedConservation.ring_gap;
      drawing.circularConservation.series.splice(
        0,
        drawing.circularConservation.series.length,
        ...nextCircularState.series
      );
    }
  }
};

export const getCommittedCanonicalRenderRequest = () => (
  committedCanonicalSession?.renderRequest || null
);

export const getCommittedCanonicalSession = () => committedCanonicalSession;

// The record count of a committed source, as the Source recipe reads it before
// it names records by #index (Export Feature Edits TSV).
export const readCommittedResourceRecordCount = (resourceId, kind) => (
  readCanonicalResourceRecordCount(committedCanonicalSession?.resources, resourceId, kind)
);

// `restore(null)` installs an artifact without a committed Session (E1: an
// empty mode slot).
export const canonicalRenderArtifactOwner = Object.freeze({
  capture: () => Object.freeze({ committedCanonicalSession, activeSessionResourceTable }),
  /** @param {{ committedCanonicalSession: any, activeSessionResourceTable: any } | null} snapshot */
  restore: (snapshot) => {
    committedCanonicalSession = snapshot?.committedCanonicalSession ?? null;
    activeSessionResourceTable = snapshot?.activeSessionResourceTable ?? null;
  }
});

/**
 * `targetDrawing` holds the Linear record layout and comparison plan that the
 * Files reconcile: a drawing, or the Load candidate.
 * @param {Record<string, any> | null} filesData
 * @param {{
 *   adoptCanonicalPayloads?: boolean, resolveRecordInputs?: boolean, targetState?: Record<string, any>,
 *   targetDrawing: Pick<DrawingState, 'linearRecordRows' | 'linearComparisonPlan'>
 * }} options
 */
const applyFiles = (filesData, {
  adoptCanonicalPayloads = false, resolveRecordInputs = true, targetState = state, targetDrawing
}) => {
  targetState.matchSequenceRegistry?.reset?.();
  targetState.circularRecordList.value = [];
  Object.assign(targetState.circularRecordDiscovery, {
    status: 'idle',
    error: '',
    inputType: '',
    primaryFile: null,
    pairedFile: null,
    canonicalRecordIdentities: []
  });
  targetState.files.c_gb = null;
  targetState.files.c_gff = null;
  targetState.files.c_fasta = null;
  targetState.files.c_depth = null;
  targetState.files.c_conservation_blasts = [];
  targetState.files.c_conservation_blasts_source = null;
  targetState.files.c_conservation_fastas = [];
  targetState.files.c_conservation_sequence_sources = [];
  targetState.files.linearCanonicalComparisons = [];
  targetState.files.d_color = null;
  targetState.files.t_color = null;
  targetState.files.blacklist = null;
  targetState.files.whitelist = null;
  targetState.files.qualifier_priority = null;
  targetState.linearReorderNotice.value = '';

  if (!filesData) {
    targetState.linearSeqs.splice(0, targetState.linearSeqs.length, ...normalizeLinearSeqList([]));
    targetDrawing.linearRecordRows.splice(0);
    replaceLinearComparisonPlan(
      targetDrawing.linearComparisonPlan,
      reconcileLinearComparisonPlan(targetDrawing.linearComparisonPlan, targetState.linearSeqs)
    );
    return { collapsedLinearSeqs: false };
  }

  targetState.files.c_gb = deserializeFile(filesData.c_gb);
  targetState.files.c_gff = deserializeFile(filesData.c_gff);
  targetState.files.c_fasta = deserializeFile(filesData.c_fasta);
  targetState.files.c_depth = deserializeFile(filesData.c_depth);
  targetState.files.c_conservation_blasts = Array.isArray(filesData.c_conservation_blasts)
    ? filesData.c_conservation_blasts.map((entry) => deserializeFile(entry)).filter(Boolean)
    : [];
  targetState.files.c_conservation_blasts_source = filesData.c_conservation_blasts_source === 'losat-cache'
    ? 'losat-cache'
    : null;
  targetState.files.c_conservation_fastas = Array.isArray(filesData.c_conservation_fastas)
    ? filesData.c_conservation_fastas.map((entry) => deserializeFile(entry))
    : [];
  targetState.files.c_conservation_sequence_sources = Array.isArray(filesData.c_conservation_sequence_sources)
    ? filesData.c_conservation_sequence_sources.map((entry) => deserializeFile(entry))
    : [];
  targetState.files.linearCanonicalComparisons = deserializeCanonicalComparisons(
    filesData.linearCanonicalComparisons,
    { adoptCanonicalPayloads }
  );
  targetState.files.d_color = deserializeFile(filesData.d_color);
  targetState.files.t_color = deserializeFile(filesData.t_color);
  targetState.files.blacklist = deserializeFile(filesData.blacklist);
  targetState.files.whitelist = deserializeFile(filesData.whitelist);
  targetState.files.qualifier_priority = deserializeFile(filesData.qualifier_priority);

  const canonicalCircularRecords = Array.isArray(filesData.circularRecords)
    ? filesData.circularRecords
    : [];
  if (canonicalCircularRecords.length > 0) {
    targetState.circularRecordDiscovery.canonicalRecordIdentities = canonicalCircularRecords
      .map((record, index) => {
        const selector = record?.region?.selector || record?.selector;
        const selectorValue = selector?.kind === 'recordId'
          ? String(selector.value || '')
          : selector?.kind === 'recordIndex'
            ? `#${Number(selector.index) + 1}`
            : `#${index + 1}`;
        return {
          selector: selectorValue,
          record_id: selector?.kind === 'recordId' ? selectorValue : '',
          recordKey: String(record?.recordKey || '')
        };
      });
    Object.assign(targetState.circularRecordDiscovery, {
      status: (targetState.cInputType.value === 'gff'
        ? targetState.files.c_gff && targetState.files.c_fasta
        : targetState.files.c_gb) ? 'deferred' : 'idle',
      error: '',
      inputType: targetState.cInputType.value,
      primaryFile: targetState.cInputType.value === 'gff'
        ? targetState.files.c_gff
        : targetState.files.c_gb,
      pairedFile: targetState.cInputType.value === 'gff'
        ? targetState.files.c_fasta
        : null
    });
  }

  if (Array.isArray(filesData.linearSeqs)) {
    const loadedLinearSeqs = filesData.linearSeqs.map((seq) => ({
      uid: seq.uid,
      gb: deserializeFile(seq.gb),
      gff: deserializeFile(seq.gff),
      fasta: deserializeFile(seq.fasta),
      depth: deserializeFile(seq.depth),
      losat_gencode: seq.losat_gencode ?? 1,
      definition: seq.definition ?? '',
      record_subtitle: seq.record_subtitle ?? '',
      file_definition: seq.file_definition ?? '',
      file_subtitle: seq.file_subtitle ?? '',
      inferred_definition: seq.inferred_definition ?? '',
      region_record_id: seq.region_record_id ?? '',
      region_start: seq.region_start ?? null,
      region_end: seq.region_end ?? null,
      region_reverse: !!seq.region_reverse
    }));
    const normalized = normalizeLinearSeqList(loadedLinearSeqs);
    const collapsed = resolveRecordInputs ? collapseEmptyLinearSeqList(loadedLinearSeqs) : loadedLinearSeqs;
    const collapsedLinearSeqs = collapsed.length !== normalized.length;
    targetState.linearSeqs.splice(0, targetState.linearSeqs.length, ...collapsed);
    if (resolveRecordInputs) {
      const rowByUid = new Map(targetDrawing.linearRecordRows.map((entry) => [String(entry?.uid || ''), entry]));
      targetDrawing.linearRecordRows.splice(
        0,
        targetDrawing.linearRecordRows.length,
        ...targetState.linearSeqs.map((seq, index) => {
          const saved = rowByUid.get(seq.uid);
          const row = Number.isInteger(Number(saved?.row)) && Number(saved.row) > 0
            ? Number(saved.row) : index + 1;
          return {
            uid: seq.uid, row,
            ...(saved?.canonicalCardinality === 'exactly_one'
              ? { canonicalCardinality: 'exactly_one' } : {}),
            ...(Number(saved?.canonicalRow) === row
              && Number.isInteger(saved?.canonicalColumn) && saved.canonicalColumn > 0
              ? { canonicalRow: row, canonicalColumn: saved.canonicalColumn }
              : {})
          };
        })
      );
    }
    const comparisonFiles = new Map(
      /** @type {[string, any][]} */ ((Array.isArray(filesData.linearComparisons) ? filesData.linearComparisons : [])
        .map((comparison) => [
          String(comparison?.id || ''),
          deserializeFile(comparison?.file)
        ]))
        .filter(([id]) => id)
    );
    const planWithFiles = normalizeLinearComparisonPlan(targetDrawing.linearComparisonPlan);
    planWithFiles.edges.forEach((edge) => {
      edge.file = comparisonFiles.get(edge.id) || null;
    });
    replaceLinearComparisonPlan(
      targetDrawing.linearComparisonPlan,
      reconcileLinearComparisonPlan(planWithFiles, targetState.linearSeqs)
    );
    return { collapsedLinearSeqs };
  }

  targetState.linearSeqs.splice(0, targetState.linearSeqs.length, ...normalizeLinearSeqList([]));
  replaceLinearComparisonPlan(
    targetDrawing.linearComparisonPlan,
    reconcileLinearComparisonPlan(targetDrawing.linearComparisonPlan, targetState.linearSeqs)
  );
  return { collapsedLinearSeqs: false };
};

// A loaded drawing's Depth series fit its own mode's Depth sources (plan 4.2
// `depth`): each drawing runs this once with its mode.
/**
 * @param {DrawingState} drawing
 * @param {'circular' | 'linear'} mode
 */
const reconcileDepthTrackStateAfterSessionFiles = (drawing, mode) => {
  const circularDepthFiles = representativeDepthFiles(state.files.c_depth);
  const circularDepthCount = circularDepthFiles.some(Boolean) ? circularDepthFiles.length : 0;
  const linearRows = state.linearSeqs.map((seq) => depthFileSlotsFromValue(seq.depth));
  const linearDepthCount = mode === 'linear'
    ? depthTrackSessionWidth({
        rows: linearRows,
        depthTracks: drawing.adv.depth_tracks,
        slots: drawing.adv.linear_track_slots
      })
    : depthTrackMatrixWidth(linearRows);
  if (mode === 'linear' && linearDepthCount > 0) {
    state.linearSeqs.forEach((seq) => {
      seq.depth = padDepthFileSlots(seq.depth, linearDepthCount);
    });
  }

  const defaults = {
    depthColor: drawing.adv.depth_color,
    depthHeight: drawing.adv.depth_height,
    largeTickInterval: drawing.adv.depth_large_tick_interval,
    smallTickInterval: drawing.adv.depth_small_tick_interval,
    tickFontSize: drawing.adv.depth_tick_font_size
  };
  let normalizedTracks;
  if (mode === 'linear') {
    normalizedTracks = normalizeDepthTracks(drawing.adv.depth_tracks, drawing.adv);
    while (normalizedTracks.length < Math.max(1, linearDepthCount)) {
      normalizedTracks.push(normalizeDepthTrackConfig(null, normalizedTracks.length, drawing.adv));
    }
  } else {
    normalizedTracks = reconcileDepthTracksToFiles(/** @type {Record<string, any>} */ ({
      files: circularDepthFiles,
      depthTracks: drawing.adv.depth_tracks,
      targetCount: Math.max(1, circularDepthCount),
      defaults
    }));
  }
  drawing.adv.depth_tracks.splice(0, drawing.adv.depth_tracks.length, ...normalizedTracks);

  drawing.adv.circular_track_slots.splice(
    0,
    drawing.adv.circular_track_slots.length,
    ...dropInvalidManagedDepthSlots(/** @type {Record<string, any>} */ ({
      slots: drawing.adv.circular_track_slots,
      activeCount: circularDepthCount
    }))
  );
  syncDepthSlotLabels(/** @type {Record<string, any>} */ ({
    slots: drawing.adv.circular_track_slots,
    depthTracks: drawing.adv.depth_tracks,
    activeCount: circularDepthCount
  }));
  drawing.adv.circular_track_slots.splice(
    0,
    drawing.adv.circular_track_slots.length,
    ...applyCircularTrackOrderPlacements(
      drawing.adv.circular_track_slots,
      drawing.adv.nt,
      drawing.form.track_type,
      drawing.adv.circular_track_slots_axis_index
    )
  );

  drawing.adv.linear_track_slots.splice(
    0,
    drawing.adv.linear_track_slots.length,
    ...dropInvalidManagedDepthSlots(/** @type {Record<string, any>} */ ({
      slots: drawing.adv.linear_track_slots,
      activeCount: linearDepthCount
    }))
  );
  syncDepthSlotLabels(/** @type {Record<string, any>} */ ({
    slots: drawing.adv.linear_track_slots,
    depthTracks: drawing.adv.depth_tracks,
    activeCount: linearDepthCount
  }));
  drawing.adv.linear_track_slots.splice(
    0,
    drawing.adv.linear_track_slots.length,
    ...applyLinearTrackOrderPlacements(
      drawing.adv.linear_track_slots,
      drawing.adv.linear_track_slots_axis_index,
      drawing.adv.nt,
      drawing.form.linear_track_layout
    )
  );
};

/** @param {DrawingState} drawing */
const cloneLiveFileState = (drawing) => ({
  files: {
    ...state.files,
    c_conservation_blasts: Array.isArray(state.files.c_conservation_blasts)
      ? [...state.files.c_conservation_blasts]
      : [],
    c_conservation_fastas: Array.isArray(state.files.c_conservation_fastas)
      ? [...state.files.c_conservation_fastas]
      : [],
    c_conservation_sequence_sources: Array.isArray(state.files.c_conservation_sequence_sources)
      ? [...state.files.c_conservation_sequence_sources]
      : [],
    linearCanonicalComparisons: (
      Array.isArray(state.files.linearCanonicalComparisons)
        ? state.files.linearCanonicalComparisons
        : []
    ).map((comparison) => (
      isResourceBackedCanonicalComparison(comparison)
        ? mapResourceBackedCanonicalComparison(comparison)
        : cloneJsonData(comparison)
    ))
  },
  circularRecordList: cloneJsonData(state.circularRecordList.value),
  circularRecordDiscovery: { ...state.circularRecordDiscovery },
  linearSeqs: state.linearSeqs.map((seq) => ({
    ...seq,
    depth: Array.isArray(seq.depth) ? [...seq.depth] : seq.depth
  })),
  linearRecordRows: drawing.linearRecordRows.map((entry) => ({ ...entry })),
  linearComparisonPlan: {
    mode: drawing.linearComparisonPlan.mode,
    defaultSource: drawing.linearComparisonPlan.defaultSource,
    edges: drawing.linearComparisonPlan.edges.map((edge) => ({ ...edge }))
  }
});

/** @param {DrawingState} drawing */
const restoreLiveFileState = (drawing, snapshot) => {
  state.matchSequenceRegistry?.reset?.();
  Object.keys(state.files).forEach((key) => {
    state.files[key] = snapshot.files[key] ?? null;
  });
  state.circularRecordList.value = cloneJsonData(snapshot.circularRecordList);
  Object.assign(state.circularRecordDiscovery, snapshot.circularRecordDiscovery);
  state.linearSeqs.splice(0, state.linearSeqs.length, ...snapshot.linearSeqs);
  drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...snapshot.linearRecordRows);
  replaceLinearComparisonPlan(drawing.linearComparisonPlan, snapshot.linearComparisonPlan);
};

const captureSessionImportTransientState = () => ({
  semanticFileWatchersSuppressed: Boolean(
    state.semanticFileWatchersSuppressed.value
  ),
  skipCaptureBaseConfig: Boolean(state.skipCaptureBaseConfig.value),
  suppressCircularMultiRecordDefaults: Boolean(
    state.suppressCircularMultiRecordDefaults.value
  ),
  linearReorderNotice: state.linearReorderNotice.value,
  rightDrawer: captureRightDrawerState(state),
  showCanvasControls: Boolean(state.showCanvasControls.value),
  isPanning: Boolean(state.isPanning.value),
  panStart: cloneJsonData(state.panStart),
  selectedAnnotation: state.selectedAnnotation.value,
  selectedSpecificPreset: state.selectedSpecificPreset.value,
  specificRulePresetLoading: Boolean(state.specificRulePresetLoading.value),
  newSpecRule: cloneJsonData(state.newSpecRule),
  newPriorityRule: cloneJsonData(state.newPriorityRule),
  newColorFeat: state.newColorFeat.value,
  newColorVal: state.newColorVal.value,
  newFeatureToAdd: state.newFeatureToAdd.value,
  newLegendCaption: state.newLegendCaption.value,
  newLegendColor: state.newLegendColor.value,
  featureSearch: state.featureSearch.value,
  labelSearch: state.labelSearch.value,
  selectedFeatureIds: Array.from(state.selectedFeatureIds.value || []),
  selectedFeatureAnchorId: state.selectedFeatureAnchorId.value,
  featureSelectionStatus: state.featureSelectionStatus.value,
  featureSelectionSuppressNextClick: Boolean(
    state.featureSelectionSuppressNextClick.value
  ),
  featureSelectionDrag: cloneJsonData(state.featureSelectionDrag),
  labelReflowLastError: state.labelReflowLastError.value,
  clickedFeature: state.clickedFeature.value,
  clickedPairwiseMatch: state.clickedPairwiseMatch.value,
  clickedLabel: state.clickedLabel.value,
  featureExtractionPending: Boolean(state.featureExtractionPending.value),
  featureExtractionError: state.featureExtractionError.value,
  featureEditorStatus: cloneJsonData(state.featureEditorStatus),
  matchSequenceSources: cloneJsonData(state.matchSequenceRegistry?.values?.() || []),
  diagramElements: [...state.diagramElements.value],
  diagramElementIds: [...state.diagramElementIds.value],
  diagramElementOriginalTransforms: new Map(
    [...state.diagramElementOriginalTransforms.value].map(([element, transform]) => [
      element,
      cloneJsonData(transform)
    ])
  ),
  legendDragging: Boolean(state.legendDragging.value),
  legendDragStart: cloneJsonData(state.legendDragStart),
  legendOriginalTransform: cloneJsonData(state.legendOriginalTransform.value),
  legendInitialTransform: cloneJsonData(state.legendInitialTransform.value),
  diagramDragging: Boolean(state.diagramDragging.value),
  diagramDragStart: cloneJsonData(state.diagramDragStart),
  lengthBarElement: state.lengthBarElement.value,
  lengthBarOriginalTransform: cloneJsonData(state.lengthBarOriginalTransform.value),
  plotTitleElement: state.plotTitleElement.value,
  plotTitleDragging: Boolean(state.plotTitleDragging.value),
  plotTitleDragStart: cloneJsonData(state.plotTitleDragStart),
  plotTitleAutoTransform: cloneJsonData(state.plotTitleAutoTransform.value)
});

const restoreSessionImportTransientState = (snapshot) => {
  state.suppressCircularMultiRecordDefaults.value =
    snapshot.suppressCircularMultiRecordDefaults;
  state.linearReorderNotice.value = snapshot.linearReorderNotice;
  restoreRightDrawerState(state, snapshot.rightDrawer);
  state.showCanvasControls.value = snapshot.showCanvasControls;
  state.isPanning.value = snapshot.isPanning;
  Object.assign(state.panStart, cloneJsonData(snapshot.panStart));
  state.selectedAnnotation.value = snapshot.selectedAnnotation;
  state.selectedSpecificPreset.value = snapshot.selectedSpecificPreset;
  state.specificRulePresetLoading.value = snapshot.specificRulePresetLoading;
  replacePlainObject(state.newSpecRule, cloneJsonData(snapshot.newSpecRule));
  replacePlainObject(state.newPriorityRule, cloneJsonData(snapshot.newPriorityRule));
  state.newColorFeat.value = snapshot.newColorFeat;
  state.newColorVal.value = snapshot.newColorVal;
  state.newFeatureToAdd.value = snapshot.newFeatureToAdd;
  state.newLegendCaption.value = snapshot.newLegendCaption;
  state.newLegendColor.value = snapshot.newLegendColor;
  state.featureSearch.value = snapshot.featureSearch;
  state.labelSearch.value = snapshot.labelSearch;
  state.selectedFeatureIds.value = new Set(snapshot.selectedFeatureIds);
  state.selectedFeatureAnchorId.value = snapshot.selectedFeatureAnchorId;
  state.featureSelectionStatus.value = snapshot.featureSelectionStatus;
  state.featureSelectionSuppressNextClick.value =
    snapshot.featureSelectionSuppressNextClick;
  Object.assign(
    state.featureSelectionDrag,
    cloneJsonData(snapshot.featureSelectionDrag)
  );
  state.labelReflowLastError.value = snapshot.labelReflowLastError;
  state.clickedFeature.value = snapshot.clickedFeature;
  state.clickedPairwiseMatch.value = snapshot.clickedPairwiseMatch;
  state.clickedLabel.value = snapshot.clickedLabel;
  state.featureExtractionPending.value = snapshot.featureExtractionPending;
  state.featureExtractionError.value = snapshot.featureExtractionError;
  Object.assign(state.featureEditorStatus, cloneJsonData(snapshot.featureEditorStatus));
  state.matchSequenceRegistry?.reset?.(cloneJsonData(snapshot.matchSequenceSources));
  state.diagramElements.value = [...snapshot.diagramElements];
  state.diagramElementIds.value = [...snapshot.diagramElementIds];
  state.diagramElementOriginalTransforms.value = new Map(
    snapshot.diagramElementOriginalTransforms
  );
  state.legendDragging.value = snapshot.legendDragging;
  Object.assign(state.legendDragStart, cloneJsonData(snapshot.legendDragStart));
  state.legendOriginalTransform.value = cloneJsonData(snapshot.legendOriginalTransform);
  state.legendInitialTransform.value = cloneJsonData(snapshot.legendInitialTransform);
  state.diagramDragging.value = snapshot.diagramDragging;
  Object.assign(state.diagramDragStart, cloneJsonData(snapshot.diagramDragStart));
  state.lengthBarElement.value = snapshot.lengthBarElement;
  state.lengthBarOriginalTransform.value = cloneJsonData(
    snapshot.lengthBarOriginalTransform
  );
  state.plotTitleElement.value = snapshot.plotTitleElement;
  state.plotTitleDragging.value = snapshot.plotTitleDragging;
  Object.assign(state.plotTitleDragStart, cloneJsonData(snapshot.plotTitleDragStart));
  state.plotTitleAutoTransform.value = cloneJsonData(snapshot.plotTitleAutoTransform);
  state.semanticFileWatchersSuppressed.value =
    snapshot.semanticFileWatchersSuppressed;
  state.skipCaptureBaseConfig.value = snapshot.skipCaptureBaseConfig;
};

// Session Load's rollback holds both drawings; the shown one is restored
// last, so the shown mode's artifact values are its own.
/** @param {DrawingState} drawing */
const captureDrawingImportSnapshot = (drawing) => ({
  config: cloneJsonData(buildConfigData(drawing)),
  ui: cloneJsonData(buildUiStateData(drawing)),
  features: buildFeatureStateData(drawing),
  editorState: buildEditorStateData(drawing),
  orthogroupState: buildOrthogroupStateData(drawing),
  importedComparisonIntent: cloneJsonData(drawing.importedComparisonIntent),
  fileLegendCaptions: Array.from(drawing.fileLegendCaptions.value || [])
});

const captureSessionImportSnapshot = () => ({
  drawings: {
    circular: captureDrawingImportSnapshot(state.drawings.circular),
    linear: captureDrawingImportSnapshot(state.drawings.linear)
  },
  ui: cloneJsonData(buildUiStateData(state.drawings[state.mode.value])),
  files: cloneLiveFileState(state.drawings.linear),
  results: state.results.value,
  collinearGroups: state.collinearGroups.value,
  runState: buildRunStateData(),
  losatCache: state.losatCache.value,
  losatDerivedCache: state.losatDerivedCache.value,
  proteinIdentityManifest: rawReactiveValue(state.proteinIdentityManifest.value)
    === adoptedProteinIdentityManifest
    ? rawReactiveValue(state.proteinIdentityManifest.value)
    : cloneJsonData(state.proteinIdentityManifest.value),
  adoptedProteinIdentityManifest,
  legacyProteinRawCandidates: cloneJsonData(state.legacyProteinRawCandidates.value),
  legacyProteinDerivedEvidence: cloneJsonData(state.legacyProteinDerivedEvidence.value),
  losatCacheInfo: cloneJsonData(state.losatCacheInfo.value),
  committedCanonicalSession: isAdoptedCanonicalSession(committedCanonicalSession)
    ? committedCanonicalSession
    : cloneCanonicalSession(committedCanonicalSession),
  activeSessionResourceTable,
  errorLog: state.errorLog.value,
  resultPanelTab: state.resultPanelTab.value,
  transients: captureSessionImportTransientState()
});

/** @param {ReturnType<typeof captureSessionImportSnapshot>} snapshot */
const restoreSessionImportSnapshot = async (snapshot) => {
  const shownMode = snapshot.ui.mode === 'linear' ? 'linear' : 'circular';
  /** @type {Array<'circular' | 'linear'>} */
  const restoreOrder = shownMode === 'linear' ? ['circular', 'linear'] : ['linear', 'circular'];
  state.sessionImportRollbackInProgress.value = true;
  try {
    state.semanticFileWatchersSuppressed.value = true;
    resetSessionBaseline();
    state.mode.value = shownMode;
    restoreOrder.forEach((drawingMode) => {
      const drawing = state.drawings[drawingMode];
      const saved = snapshot.drawings[drawingMode];
      applyConfigData(drawing, saved.config, { resolveTrackPlacements: false });
      applyUiStateData(drawing, saved.ui);
    });
    applyUiStateData(state.drawings[shownMode], snapshot.ui);
    restoreLiveFileState(state.drawings.linear, snapshot.files);
    state.losatCache.value = snapshot.losatCache;
    state.losatDerivedCache.value = snapshot.losatDerivedCache;
    adoptedProteinIdentityManifest = snapshot.adoptedProteinIdentityManifest;
    state.proteinIdentityManifest.value = snapshot.proteinIdentityManifest === adoptedProteinIdentityManifest
      ? snapshot.proteinIdentityManifest
      : cloneJsonData(snapshot.proteinIdentityManifest);
    state.legacyProteinRawCandidates.value = cloneJsonData(snapshot.legacyProteinRawCandidates);
    state.legacyProteinDerivedEvidence.value = cloneJsonData(snapshot.legacyProteinDerivedEvidence);
    state.losatCacheInfo.value = cloneJsonData(snapshot.losatCacheInfo);
    committedCanonicalSession = isAdoptedCanonicalSession(snapshot.committedCanonicalSession)
      ? snapshot.committedCanonicalSession
      : cloneCanonicalSession(snapshot.committedCanonicalSession);
    activeSessionResourceTable = snapshot.activeSessionResourceTable;
    state.skipCaptureBaseConfig.value = true;
    applyResultsData(snapshot.results, snapshot.ui);
    restoreOrder.forEach((drawingMode) => {
      const drawing = state.drawings[drawingMode];
      const saved = snapshot.drawings[drawingMode];
      Object.assign(
        drawing.importedComparisonIntent,
        createImportedComparisonIntentState(),
        cloneJsonData(saved.importedComparisonIntent)
      );
      applyFeatureStateData(drawing, saved.features);
      applyOrthogroupStateData(drawing, saved.orthogroupState);
      applyEditorStateData(drawing, saved.editorState, { normalized: true });
      drawing.fileLegendCaptions.value = new Set(saved.fileLegendCaptions);
    });
    state.collinearGroups.value = snapshot.collinearGroups;
    applyRunStateData(snapshot.runState);
    state.errorLog.value = snapshot.errorLog;
    state.resultPanelTab.value = snapshot.resultPanelTab;
    await nextTick();
    recordSessionLifecycleEvent('session-rollback-source-restored');
    restoreSessionImportTransientState(snapshot.transients);
    recordSessionLifecycleEvent('session-rollback-transients-reconciled');
    await nextTick();
  } finally {
    state.sessionImportRollbackInProgress.value = false;
  }
};

const clearObject = (target) => {
  Object.keys(target).forEach((key) => {
    delete target[key];
  });
};

// Both drawings return to their defaults before a Session installs its own.
const resetSessionBaseline = () => {
  activePreviewRuntime?.clearActiveRuntime?.();
  preservedCliOptions = null;
  committedCanonicalSession = null;
  activeSessionResourceTable = null;
  adoptedProteinIdentityManifest = null;
  state.sessionResourceDiscoveryDeferred.value = false;
  resetSettingsState(state);
  resetLayoutState(state);
  resetRightDrawerState(state);
  state.mode.value = 'circular';
  state.cInputType.value = 'gb';
  state.lInputType.value = 'gb';
  state.sessionTitle.value = '';
  state.errorLog.value = null;
  state.results.value = [];
  state.failedGeneratePreservedResult.value = false;
  state.selectedResultIndex.value = 0;
  state.resultPanelTab.value = 'preview';
  state.lastRunInfo.value = null;
  state.trackSlotResolvedGeometry.value = null;
  state.annotationWarnings.value = [];
  state.featureIdentityNotices.value = [];
  state.featureEditRemovalCount.value = 0;
  state.comparisonWarnings.value = [];
  applyFiles(null, { targetDrawing: state.drawings.linear });
  state.losatCache.value = new Map();
  state.losatDerivedCache.value = new Map();
  state.proteinIdentityManifest.value = emptyProteinIdentityManifest();
  state.legacyProteinRawCandidates.value = { schema: 1, entries: [] };
  state.legacyProteinDerivedEvidence.value = { schema: 1, entries: [] };
  state.losatCacheInfo.value = [];
  state.orthogroups.value = [];
  state.collinearGroups.value = [];
  state.featureOrthogroupIndex.value = new Map();
  state.selectedOrthogroupId.value = '';
  state.selectedOrthogroupAlignmentFeature.value = '';
  state.extractedFeatures.value = [];
  if (state.biologicalFeatures) state.biologicalFeatures.value = [];
  state.featureRecordIds.value = [];
  state.selectedFeatureRecordIdx.value = 0;
  Object.values(state.drawings).forEach((/** @type {DrawingState} */ drawing) => {
    Object.assign(drawing.importedComparisonIntent, createImportedComparisonIntentState());
    clearObject(drawing.orthogroupNameOverrides);
    clearObject(drawing.orthogroupDescriptionOverrides);
    clearObject(drawing.orthogroupDormantOverrides);
    clearObject(drawing.featureColorOverrides);
    drawing.featureVisibilityManualRules.splice(0);
    clearObject(drawing.featureOverrides);
    clearObject(drawing.featureStrokeOverrides);
    drawing.canonicalLabelOverrideRows.value = [];
    clearObject(drawing.labelTextBulkOverrides);
    drawing.legendEntries.value = [];
    drawing.dormantLegendEntries.value = [];
  });
  state.generatedMode.value = 'circular';
  state.generatedLegendPosition.value = 'left';
  state.generatedMultiRecordCanvas.value = false;
  state.generatedCircularPlotTitlePosition.value = 'none';
};

/** @param {DrawingState} drawing */
export const buildUiStateData = (drawing, { includePreviewNavigation = true } = {}) => {
  const ui = {
    title: String(state.sessionTitle.value || ''),
    mode: state.mode.value,
    canvasPadding: { ...drawing.canvasPadding },
    selectedResultIndex: state.selectedResultIndex.value,
    generatedLegendPosition: state.generatedLegendPosition.value,
    generatedMode: state.generatedMode.value,
    generatedMultiRecordCanvas: Boolean(state.generatedMultiRecordCanvas.value),
    generatedCircularPlotTitlePosition: normalizeCircularPlotTitlePosition(
      state.generatedCircularPlotTitlePosition.value
    ),
    layoutPreferences: cloneJsonData(drawing.layoutPreferences),
    featurePanelTab: state.featurePanelTab.value,
    cInputType: state.cInputType.value,
    lInputType: state.lInputType.value,
    losatProgram: drawing.losatProgram.value,
    downloadDpi: state.downloadDpi.value,
    autoLabelReflow: Boolean(state.autoLabelReflowEnabled.value),
    linearTypographyLinked: Boolean(drawing.linearTypographyLinked.value),
    paletteInstantPreviewEnabled: Boolean(state.paletteInstantPreviewEnabled.value),
    // App-level settings (History and the Load rollback restore them with the UI).
    losatExecution: buildLosatExecutionData(),
    richFeaturePopup: Boolean(state.richFeaturePopup.value),
    appliedPaletteName: state.appliedPaletteName.value,
    appliedPaletteColors: cloneColors(state.appliedPaletteColors.value),
    pendingPaletteName: drawing.pendingPaletteName.value,
    pendingPaletteColors: cloneColors(drawing.pendingPaletteColors.value),
    legendCurrentOffset: { ...state.legendCurrentOffset },
    diagramOffset: { ...state.diagramOffset },
    lengthBarUserOffset: { ...state.lengthBarUserOffset },
    plotTitleUserOffset: { ...state.plotTitleUserOffset }
  };

  if (includePreviewNavigation) {
    ui.zoom = state.zoom.value;
    ui.canvasPan = { x: state.canvasPan.x, y: state.canvasPan.y };
  }

  return ui;
};

// `mode` is not restored here: History restores run the mode transition
// first, and Session rollback sets it with the restored document.
/** @param {DrawingState} drawing */
export const applyUiStateData = (drawing, ui = {}, { restorePreviewNavigation = true } = {}) => {
  if (typeof ui.title === 'string') state.sessionTitle.value = ui.title;
  if (ui.cInputType) state.cInputType.value = ui.cInputType;
  if (ui.lInputType) state.lInputType.value = ui.lInputType;
  if (ui.losatProgram) {
    const program = String(ui.losatProgram);
    drawing.losatProgram.value = ['blastn', 'tblastx', 'blastp'].includes(program) ? program : 'blastn';
  }
  if (ui.downloadDpi) state.downloadDpi.value = ui.downloadDpi;
  state.autoLabelReflowEnabled.value = Boolean(ui.autoLabelReflow);
  reconcileImportedLinearTypographyLink({
    adv: drawing.adv,
    linked: drawing.linearTypographyLinked,
    ui
  });
  state.paletteInstantPreviewEnabled.value = Boolean(ui.paletteInstantPreviewEnabled);
  if (isPlainObject(ui.losatExecution)) applyLosatExecutionData(ui.losatExecution);
  if (typeof ui.richFeaturePopup === 'boolean') state.richFeaturePopup.value = ui.richFeaturePopup;
  if (ui.featurePanelTab === 'labels' || ui.featurePanelTab === 'colors') {
    state.featurePanelTab.value = ui.featurePanelTab;
  }

  if (ui.generatedMode) state.generatedMode.value = ui.generatedMode === 'linear' ? 'linear' : 'circular';
  if (ui.generatedLegendPosition) {
    state.generatedLegendPosition.value = normalizeLegendPosition(
      ui.generatedLegendPosition,
      state.generatedMode.value === 'linear' ? 'bottom' : 'left'
    );
  }
  if (Object.prototype.hasOwnProperty.call(ui, 'generatedMultiRecordCanvas')) {
    state.generatedMultiRecordCanvas.value = Boolean(ui.generatedMultiRecordCanvas);
  }
  if (ui.generatedCircularPlotTitlePosition || ui.circularPlotTitlePosition) {
    state.generatedCircularPlotTitlePosition.value = hasStoredLayoutValue(ui.generatedCircularPlotTitlePosition)
      ? normalizeCircularPlotTitlePosition(ui.generatedCircularPlotTitlePosition)
      : normalizeCircularPlotTitlePosition(ui.circularPlotTitlePosition);
  }

  restorePaletteStateFromSession(drawing, ui);
  restoreLayoutPreferences(drawing, ui);

  if (ui.legendCurrentOffset) {
    state.legendCurrentOffset.x = Number(ui.legendCurrentOffset.x) || 0;
    state.legendCurrentOffset.y = Number(ui.legendCurrentOffset.y) || 0;
  }
  if (ui.diagramOffset) {
    state.diagramOffset.x = Number(ui.diagramOffset.x) || 0;
    state.diagramOffset.y = Number(ui.diagramOffset.y) || 0;
  }
  if (ui.lengthBarUserOffset) {
    state.lengthBarUserOffset.x = Number(ui.lengthBarUserOffset.x) || 0;
    state.lengthBarUserOffset.y = Number(ui.lengthBarUserOffset.y) || 0;
  }
  if (ui.plotTitleUserOffset) {
    state.plotTitleUserOffset.x = Number(ui.plotTitleUserOffset.x) || 0;
    state.plotTitleUserOffset.y = Number(ui.plotTitleUserOffset.y) || 0;
  }

  if (ui.canvasPadding) {
    drawing.canvasPadding.top = Number(ui.canvasPadding.top) || 0;
    drawing.canvasPadding.right = Number(ui.canvasPadding.right) || 0;
    drawing.canvasPadding.bottom = Number(ui.canvasPadding.bottom) || 0;
    drawing.canvasPadding.left = Number(ui.canvasPadding.left) || 0;
  }
  if (restorePreviewNavigation) {
    if (ui.canvasPan) {
      state.canvasPan.x = Number(ui.canvasPan.x) || 0;
      state.canvasPan.y = Number(ui.canvasPan.y) || 0;
    }
    if (typeof ui.zoom === 'number') state.zoom.value = ui.zoom;
  }

};

export const applyResultsData = (resultsData = [], ui = {}) => {
  state.failedGeneratePreservedResult.value = false;
  if (Array.isArray(resultsData)) {
    const logicalResults = resultsData.every(isCommittedSvgResult)
      ? resultsData
      : normalizeLogicalResults(resultsData.map((res, idx) => (
          isCommittedSvgResult(res)
            ? res
            : { name: res?.name || `Result ${idx + 1}`, content: res?.content || '' }
        )));
    const committedCount = logicalResults.filter(isCommittedSvgResult).length;
    if (committedCount !== 0 && committedCount !== logicalResults.length) {
      throw new Error('Committed and imported SVG Results cannot be mixed.');
    }
    // Results restored without their committed metadata (a History
    // checkpoint) carry the feature types of the committed request (R13).
    state.results.value = committedCount === logicalResults.length
      ? logicalResults
      : admitLegacyImportedResults(createLegacyImportResultSource(logicalResults), {
          selectedFeatureTypes: committedCanonicalSession?.renderRequest?.diagramOptions?.selectedFeaturesSet
        });
  } else {
    state.results.value = [];
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
};

/** @param {DrawingState} drawing */
export const buildFeatureStateData = (drawing) => ({
  extractedFeatures: sanitizeExtractedFeaturesForSession(state.extractedFeatures.value),
  biologicalFeatures: sanitizeExtractedFeaturesForSession(state.biologicalFeatures?.value),
  featureRecordIds: cloneJsonData(state.featureRecordIds.value),
  selectedFeatureRecordIdx: state.selectedFeatureRecordIdx.value,
  featureColorOverrides: cloneJsonData(drawing.featureColorOverrides),
  featureVisibilityManualRules: normalizeFeatureVisibilityRulesForSession(drawing.featureVisibilityManualRules),
  featureOverrides: featureOverridesForState(drawing.featureOverrides),
  labelOverrideRows: cloneJsonData(drawing.canonicalLabelOverrideRows.value),
  labelTextBulkOverrides: cloneJsonData(drawing.labelTextBulkOverrides)
});

// The catalog features of the shown Result (artifact values).
const applyFeatureArtifactData = (features = {}) => {
  state.extractedFeatures.value = Array.isArray(features.extractedFeatures)
    ? features.extractedFeatures
    : [];
  if (state.biologicalFeatures) {
    state.biologicalFeatures.value = Array.isArray(features.biologicalFeatures)
      ? features.biologicalFeatures
      : [];
  }
  state.featureRecordIds.value = Array.isArray(features.featureRecordIds)
    ? features.featureRecordIds
    : [];
};

/** @param {DrawingState} drawing */
export const applyFeatureStateData = (drawing, features = {}) => {
  applyFeatureArtifactData(features);
  state.selectedFeatureRecordIdx.value = Number.isInteger(features.selectedFeatureRecordIdx)
    ? features.selectedFeatureRecordIdx
    : 0;
  applyDrawingFeatureData(drawing, features);
};

// A drawing's per-feature edits as saved (a Session 46 slice's `features`, a
// flat draft's, or a History checkpoint's).
/** @param {DrawingState} drawing */
const applyDrawingFeatureData = (drawing, features = {}) => {
  replacePlainObject(drawing.featureColorOverrides, cloneJsonObject(features.featureColorOverrides));
  drawing.featureVisibilityManualRules.splice(
    0,
    drawing.featureVisibilityManualRules.length,
    ...normalizeFeatureVisibilityRulesForSession(features.featureVisibilityManualRules || [])
  );
  replacePlainObject(drawing.featureOverrides, featureOverridesForState(features.featureOverrides));
  drawing.canonicalLabelOverrideRows.value = Array.isArray(features.labelOverrideRows)
    ? cloneJsonData(features.labelOverrideRows)
    : [];
  replacePlainObject(drawing.labelTextBulkOverrides, cloneStringMap(features.labelTextBulkOverrides));
};

/** @param {DrawingState} drawing */
export const buildOrthogroupStateData = (drawing) => ({
  groups: Array.isArray(state.orthogroups.value) ? cloneJsonData(state.orthogroups.value) : [],
  selectedOrthogroupId: String(state.selectedOrthogroupId.value || ''),
  orthogroupNameOverrides: cloneStringMap(drawing.orthogroupNameOverrides),
  orthogroupDescriptionOverrides: cloneStringMap(drawing.orthogroupDescriptionOverrides),
  orthogroupDormantOverrides: normalizeOrthogroupDormantOverrides(drawing.orthogroupDormantOverrides)
});

export const buildRunStateData = () => ({
  lastRunInfo: cloneJsonData(state.lastRunInfo.value),
  pairwiseMatchFactors: cloneJsonObject(state.pairwiseMatchFactors.value)
});

export const applyRunStateData = (runState = {}) => {
  state.lastRunInfo.value = runState.lastRunInfo ? cloneJsonData(runState.lastRunInfo) : null;
  state.pairwiseMatchFactors.value = cloneJsonObject(runState.pairwiseMatchFactors);
};

const setFeatureEditorStatusData = (updates = {}) => {
  if (!state.featureEditorStatus || typeof state.featureEditorStatus !== 'object') return;
  Object.assign(state.featureEditorStatus, {
    status: updates.status ?? state.featureEditorStatus.status,
    generationId: updates.generationId ?? state.featureEditorStatus.generationId,
    error: updates.error === undefined ? state.featureEditorStatus.error : updates.error,
    summaryCount: updates.summaryCount ?? state.featureEditorStatus.summaryCount,
    detailsCacheSize: updates.detailsCacheSize ?? state.featureEditorStatus.detailsCacheSize
  });
};

const synchronizeRestoredFeatureSummaryStatus = ({ generationId = 'session-load' } = {}) => {
  const summaryCount = Array.isArray(state.extractedFeatures.value)
    ? state.extractedFeatures.value.length
    : 0;
  if (summaryCount === 0) return false;
  state.featureExtractionPending.value = false;
  state.featureExtractionError.value = null;
  setFeatureEditorStatusData({
    status: 'summary-ready',
    generationId,
    error: null,
    summaryCount
  });
  return true;
};

/** @param {DrawingState} drawing */
const applySessionFeatureRecoveryPlan = (drawing, plan, { generationId = 'session-feature-recovery' } = {}) => {
  state.featureExtractionPending.value = false;

  if (plan?.status === 'recovered' || plan?.status === 'aligned') {
    if (plan.recoveredFeatureState) applyFeatureStateData(drawing, plan.recoveredFeatureState);
    if (plan.migratedEditorState) applyEditorStateData(drawing, plan.migratedEditorState);
    state.featureExtractionError.value = null;
    setFeatureEditorStatusData({
      status: 'summary-ready',
      generationId,
      error: plan.warning || null,
      summaryCount: Array.isArray(plan.recoveredFeatureState?.extractedFeatures)
        ? plan.recoveredFeatureState.extractedFeatures.length
        : state.extractedFeatures.value.length
    });
    return;
  }

  if (plan?.status === 'unrecoverable' || plan?.status === 'failed') {
    const warning = plan.warning || 'Feature metadata recovery failed. The SVG preview remains available.';
    state.featureExtractionError.value = { summary: warning, details: [] };
    setFeatureEditorStatusData({
      status: 'failed',
      generationId,
      error: warning,
      summaryCount: 0
    });
    return;
  }

  synchronizeRestoredFeatureSummaryStatus({ generationId });
};

// E1: one mode's generated artifact as Save reads it: the installed artifact
// of the displayed mode, or the other mode's slot, which the composition root
// keeps while that mode is not shown.
/**
 * @typedef {Readonly<Record<string, any>> & {
 *   mode: 'circular' | 'linear',
 *   committedCanonicalSession: any,
 *   activeSessionResourceTable: any
 * }} SessionArtifactView
 */
/** @returns {SessionArtifactView} */
const displayedArtifactView = () => Object.freeze({
  ...Object.fromEntries(ARTIFACT_SLOT_KEYS.map((key) => [key, state[key].value])),
  originalLegendOrder: state.originalLegendOrder.value,
  mode: state.generatedMode.value === 'linear' ? 'linear' : 'circular',
  committedCanonicalSession,
  activeSessionResourceTable
});
/**
 * @param {Readonly<ArtifactSlot>} slot
 * @returns {SessionArtifactView}
 */
const stashedArtifactView = (slot) => Object.freeze({
  ...slot.values,
  originalLegendOrder: slot.legendInventory,
  mode: slot.mode,
  committedCanonicalSession: slot.runtimeState?.canonical?.committedCanonicalSession ?? null,
  activeSessionResourceTable: slot.runtimeState?.canonical?.activeSessionResourceTable ?? null
});
/** @param {SessionArtifactView | null} view */
const artifactHasResult = (view) => Array.isArray(view?.results) && view.results.length > 0;
// Save writes every Result (Owner, 2026-10-07). The top level holds the shown
// mode's set when it has a Result, otherwise the other mode's set;
// `otherModeResult` holds the remaining set only when it has a Result.
/** @param {Readonly<ArtifactSlot> | null} otherSlot */
const chooseSessionArtifacts = (otherSlot) => {
  const displayed = displayedArtifactView();
  const other = otherSlot && otherSlot.mode !== displayed.mode ? stashedArtifactView(otherSlot) : null;
  if (!other || !artifactHasResult(other)) return { artifact: displayed, otherArtifact: null };
  if (!artifactHasResult(displayed)) return { artifact: other, otherArtifact: null };
  return { artifact: displayed, otherArtifact: other };
};
// Live edits are committed into a Result's content at once, so its content is
// what Save writes, also for the displayed Result.
/** @param {SessionArtifactView} view */
const serializeArtifactResults = (view) => normalizeLogicalResults(view.results.map(
  (/** @type {Record<string, any>} */ res, /** @type {number} */ idx) => ({
    name: res.name || `Result ${idx + 1}`,
    content: res.content
  })
));
// The artifact's own catalog, generated Legend inventory and colors, stroke
// defaults and alignment Reset receipt; the Legend and stroke edits are its
// mode's drawing's (Session 46 `modes.<m>.editorState`).
/** @param {SessionArtifactView} view */
const buildArtifactEditorState = (view) => {
  return {
    legend: {
      originalOrder: cloneJsonArray(view.originalLegendOrder),
      originalColors: cloneStringMap(view.originalLegendColors)
    },
    originalSvgStroke: {
      color: view.originalSvgStroke?.color ?? null,
      width: view.originalSvgStroke?.width ?? null
    },
    alignmentResetReceipt: cloneJsonValue(view.similarityAlignmentResetReceipt, null),
    featureCatalog: admittedFeatureCatalog(view.featureCatalog)
  };
};
// D-25 (PD-OI-079): Results need their admitted catalog; a Result without one
// (a legacy Session) is saved after a Generate in its mode. The error names that
// mode when it is not the one shown, and its Generate action runs there.
/**
 * @param {Record<string, any>} editorState
 * @param {Record<string, any>[]} logicalResults
 * @param {'circular' | 'linear'} mode
 */
const admitSavedFeatureCatalog = (editorState, logicalResults, mode) => {
  if (logicalResults.length === 0) {
    editorState.featureCatalog = null;
    return;
  }
  const requiresGenerate = () => sessionSaveRequiresGenerate(
    mode === state.mode.value ? {} : { diagramMode: mode }
  );
  if (!editorState.featureCatalog) throw requiresGenerate();
  try {
    editorState.featureCatalog = validateFeatureCatalog(
      editorState.featureCatalog,
      logicalResults,
      { adopt: true, mode }
    );
  } catch (error) {
    console.warn('Session feature catalog validation failed.', normalizeUserFacingError(error));
    throw requiresGenerate();
  }
};
/**
 * @param {Record<string, any> | null} committed
 * @param {SessionArtifactView} view
 * @param {FeatureCatalog | null} featureCatalog
 */
const promoteSavedCanonicalSession = (committed, view, featureCatalog) => {
  if (!committed || committed.renderRequest.schema >= CANONICAL_REQUEST_SCHEMA) return committed;
  const promoted = {
    ...committed,
    renderRequest: promoteCanonicalRenderRequestToCurrent(committed.renderRequest, {
      featureCatalog,
      legacyOrthogroupState: { groups: cloneJsonData(view.orthogroups || []) }
    })
  };
  return isAdoptedCanonicalSession(committed) ? adoptRuntimeCanonicalSession(promoted) : promoted;
};
/** @param {SessionArtifactView} view */
const artifactRunMetadata = (view) => ({
  ...(view.trackSlotResolvedGeometry
    ? { trackSlotGeometry: cloneJsonData(view.trackSlotResolvedGeometry) } : {}),
  annotationWarnings: cloneJsonData(view.annotationWarnings),
  ...(view.featureIdentityNotices?.length
    ? { featureIdentityNotices: cloneJsonData(view.featureIdentityNotices) } : {}),
  ...(view.comparisonWarnings?.length
    ? { comparisonWarnings: cloneJsonData(view.comparisonWarnings) } : {})
});
/** @param {SessionArtifactView} view */
const artifactCliInvocation = (view) => {
  const invocation = view.lastRunInfo?.invocation;
  return isCliInvocationSessionExportable(invocation) ? cloneJsonData(invocation) : undefined;
};
// The generated layout and palette of an artifact's Results (`ui` fields).
/** @param {SessionArtifactView} view */
const artifactUiState = (view) => ({
  selectedResultIndex: view.selectedResultIndex,
  generatedLegendPosition: view.generatedLegendPosition,
  generatedMultiRecordCanvas: Boolean(view.generatedMultiRecordCanvas),
  generatedCircularPlotTitlePosition: normalizeCircularPlotTitlePosition(
    view.generatedCircularPlotTitlePosition
  ),
  appliedPaletteName: view.appliedPaletteName,
  appliedPaletteColors: cloneColors(view.appliedPaletteColors)
});
// The other mode's set: its Results, admitted catalog and committed request
// are written as one unit beside the top-level set (PD-OI-045).
/** @param {SessionArtifactView} view */
const prepareOtherModeArtifact = (view) => {
  const results = serializeArtifactResults(view);
  const editorState = buildArtifactEditorState(view);
  admitSavedFeatureCatalog(editorState, results, view.mode);
  if (!view.committedCanonicalSession) throw sessionSaveRequiresGenerate({ diagramMode: view.mode });
  const committed = promoteSavedCanonicalSession(
    isAdoptedCanonicalSession(view.committedCanonicalSession)
      ? view.committedCanonicalSession
      : cloneCanonicalSession(view.committedCanonicalSession),
    view,
    editorState.featureCatalog
  );
  return { view, results, editorState, committed };
};
/**
 * @param {ReturnType<typeof prepareOtherModeArtifact>} other
 * @param {CanonicalRenderRequest} renderRequest Its resource references name the Session's table.
 * @returns {SessionOtherModeResult}
 */
const buildOtherModeResult = (other, renderRequest) => {
  const cliInvocation = artifactCliInvocation(other.view);
  return {
    renderRequest,
    results: other.results,
    editorState: {
      featureCatalog: other.editorState.featureCatalog,
      alignmentResetReceipt: other.editorState.alignmentResetReceipt ?? null,
      legend: {
        originalOrder: other.editorState.legend.originalOrder,
        originalColors: other.editorState.legend.originalColors
      },
      originalSvgStroke: other.editorState.originalSvgStroke
    },
    ui: artifactUiState(other.view),
    runMetadata: artifactRunMetadata(other.view),
    ...(cliInvocation ? { cliInvocation } : {})
  };
};

/**
 * @typedef {{
 *   drawing: DrawingState,
 *   linearRecordCatalog?: any,
 *   recordDisplayRows?: any,
 *   modes: Record<'circular' | 'linear', SessionModeSlice>,
 *   savedUi: Record<string, any>,
 *   isCurrent: () => boolean,
 *   artifact: SessionArtifactView,
 *   otherArtifact: SessionArtifactView | null
 * }} ExportSessionDocumentOptions
 */

/**
 * @param {string | null | undefined} titleOverride
 * @param {ExportSessionDocumentOptions} options
 */
const exportSessionDocument = async (
  titleOverride = null,
  {
    drawing, linearRecordCatalog = null, recordDisplayRows = null, modes, savedUi, isCurrent, artifact, otherArtifact
  }
) => {
  const resolvedTitle =
    typeof titleOverride === 'string'
      ? titleOverride.trim()
      : typeof state.sessionTitle?.value === 'string'
        ? state.sessionTitle.value.trim()
        : '';
  const sessionFilename = buildSessionFilename(resolvedTitle);
  if (lastSessionFilename && lastSessionFilename === sessionFilename) {
    const proceed = confirm(`Download "${sessionFilename}" again? Your browser may overwrite or rename the file.`);
    if (!proceed) {
      recordSessionLifecycleEvent('session-save-download-canceled', {
        reason: 'repeat-download'
      });
      return { status: 'canceled' };
    }
  }

  recordSessionLifecycleEvent('session-save-projection-start');
  const logicalResults = serializeArtifactResults(artifact);
  const editorState = buildArtifactEditorState(artifact);
  admitSavedFeatureCatalog(editorState, logicalResults, artifact.mode);
  const other = otherArtifact ? prepareOtherModeArtifact(otherArtifact) : null;

  const {
    entries: losatEntries,
    validatedManifest,
    manifestValidated
  } = serializeLosatCache();
  const exportableCliInvocation = artifactCliInvocation(artifact);
  // Each slice's unmanaged overrides are checked against its own mode (PD-OI-009, OV-106).
  for (const mode of SLICE_MODES) {
    const config = modes[mode].config;
    Object.assign(config.adv, normalizedArrowGeometryState(config.adv));
    config.unmanagedConfigOverrides = await validateUnmanagedConfigOverrides({
      mode,
      configOverrides: config.unmanagedConfigOverrides,
      requireUnmanagedOnly: true
    });
  }
  let committed = isAdoptedCanonicalSession(artifact.committedCanonicalSession)
    ? artifact.committedCanonicalSession
    : cloneCanonicalSession(artifact.committedCanonicalSession);
  const settingsOnly = !artifact.committedCanonicalSession && logicalResults.length === 0
    && !hasBiologicalSessionInputs({ ...state.files, linearSeqs: state.linearSeqs });
  if (committed) {
    try {
      const adoptedCommitted = isAdoptedCanonicalSession(committed);
      const projected = projectCanonicalSessionRequest({
        ...committed,
        sessionResourceTable: adoptedCommitted ? artifact.activeSessionResourceTable : null,
        deferResourceContent: adoptedCommitted,
        adoptCanonicalPayloads: adoptedCommitted
      });
      validateCurrentWriterActiveConfig({
        mode: projected.mode,
        storedConfig: modes[projected.mode === 'linear' ? 'linear' : 'circular'].config
      });
    } catch (error) {
      console.warn('Session active configuration validation failed.', normalizeUserFacingError(error));
      throw recognizedCauseOr(error, diagnosticError('INPUT_INVALID', { field: 'config', reason: 'FIELDS' }));
    }
  }
  if (!committed && !settingsOnly) {
    const comparisonPlanSnapshot = state.mode.value === 'linear'
      ? resolveLinearComparisonPlan({
          plan: drawing.linearComparisonPlan,
          sequences: normalizeLinearSeqList(state.linearSeqs),
          layout: drawing.linearRecordLayoutEnabled?.value
            ? drawing.linearRecordRows
            : [],
          losatProgram: drawing.losatProgram?.value,
          blastpMode: drawing.losat?.blastp?.mode
        })
      : null;
    const activeFiles = await serializeActiveRenderFiles(
      state.mode.value,
      state,
      drawing,
      {
        comparisonPlan: comparisonPlanSnapshot,
        linearRecordCatalog
      }
    );
    committed = buildCanonicalRenderRequest({
      state,
      drawing,
      filesData: activeFiles,
      recordDisplayRows: recordDisplayRows?.value || [],
      comparisonPlanSnapshot
    });
  }
  committed = promoteSavedCanonicalSession(committed, artifact, editorState.featureCatalog);
  const canonical = await assembleSessionResources(state, committed, state.drawings.linear, other?.committed ?? null);
  await validateSimilarityAlignmentResetReceipt(editorState.alignmentResetReceipt, canonical);
  if (other) {
    await validateSimilarityAlignmentResetReceipt(other.editorState.alignmentResetReceipt, {
      renderRequest: canonical.otherRenderRequest, resources: canonical.resources
    });
  }
  const legacyRawCandidates = serializableLegacyProteinCandidateEnvelope(
    state.legacyProteinRawCandidates.value
  );
  const legacyDerivedEvidence = normalizeLegacyDerivedEvidence(
    state.legacyProteinDerivedEvidence.value
  );
  if (
    (!manifestValidated || validatedManifest !== state.proteinIdentityManifest.value)
    && !validateProteinIdentityManifest(state.proteinIdentityManifest.value)
  ) {
    throw new Error('Save Session requires a valid protein identity manifest.');
  }
  /** @type {GbdrawSession} */
  const sessionData = /** @satisfies {GbdrawSession} */ ({
    format: 'gbdraw-session',
    version: SESSION_VERSION,
    createdAt: new Date().toISOString(),
    title: resolvedTitle || undefined,
    ui: savedUi,
    renderRequest: canonical.renderRequest,
    resources: canonical.resources,
    webFiles: canonical.webFiles,
    results: logicalResults,
    ...(!settingsOnly ? { runMetadata: artifactRunMetadata(artifact) } : {}),
    ...(other ? { otherModeResult: buildOtherModeResult(other, canonical.otherRenderRequest) } : {}),
    editorState,
    // The group names are the Linear drawing's (`modes.linear.config.webEdits`).
    orthogroupState: {
      selectedOrthogroupId: String(state.selectedOrthogroupId.value || ''),
      selectedOrthogroupAlignmentFeature: String(state.selectedOrthogroupAlignmentFeature.value || '')
    },
    modes,
    ...(preservedCliOptions ? { cliOptions: cloneJsonData(preservedCliOptions) } : {}),
    losatCache: {
      entries: losatEntries
    },
    losatDerivedCache: {
      entries: []
    },
    proteinIdentityManifest: rawReactiveValue(state.proteinIdentityManifest.value)
      === adoptedProteinIdentityManifest
      ? rawReactiveValue(state.proteinIdentityManifest.value)
      : cloneJsonData(state.proteinIdentityManifest.value),
    cliInvocation: exportableCliInvocation
  });

  const legacyArtifacts = buildSessionLegacyArtifacts({
    legacyRawCandidates,
    legacyDerivedEvidence
  });
  if (legacyArtifacts) {
    sessionData.legacyArtifacts = legacyArtifacts;
  }
  try {
    validateSessionAuthorityInventory(sessionData, SESSION_VERSION);
  } catch (error) {
    console.error('Session writer validation failed.', normalizeUserFacingError(error));
    throw recognizedCauseOr(error, diagnosticError('INPUT_INVALID', { field: 'schema' }));
  }

  recordSessionLifecycleEvent('session-save-projection-end');
  recordSessionLifecycleEvent('session-save-compression-start');
  const compressed = await compressSessionData(sessionData);
  recordSessionLifecycleEvent('session-save-compression-end', {
    compressedSize: compressed.size
  });
  if (!isCurrent()) return { status: 'canceled' };
  if (!confirmLargeSessionBlob(compressed)) {
    recordSessionLifecycleEvent('session-save-download-canceled', {
      reason: 'large-download',
      compressedSize: compressed.size
    });
    return { status: 'canceled', compressedSize: compressed.size };
  }
  if (!isCurrent()) return { status: 'canceled' };
  downloadBlob(compressed, sessionFilename);
  recordSessionLifecycleEvent('session-save-download-handoff-completed', {
    compressedSize: compressed.size
  });
  lastSessionFilename = sessionFilename;
  return { status: 'saved', blob: compressed, filename: sessionFilename };
};

/**
 * @param {Record<string, any>} [options] The options of `importSession` with
 * `isCurrent` and `signal`; `transformLegacyResultSvg` is the composition
 * root's `LegacyResultSvgTransform` port.
 */
const importSessionDocument = async (e, options = {}) => {
  const file = e.target.files[0];
  if (!file) return { status: 'skipped' };
  recordSessionLifecycleEvent('sessionSelection');

  const semanticFileWatchersSuppressedBeforeImport = Boolean(
    state.semanticFileWatchersSuppressed.value
  );
  const rollbackStateExtension = options?.rollbackState;
  /** @type {ReturnType<typeof captureSessionImportSnapshot> | null} */
  let rollbackSnapshot = null;
  let rollbackExtensionSnapshot;
  let commitStarted = false;
  const previousAlert = state.errorLog.value;

  try {
    const candidate = await importSessionFile(file, { signal: options.signal });
    if (!options.isCurrent()) return { status: 'canceled' };
    recordSessionLifecycleEvent('session-import-codec-completed', candidate.timings);
    let data = candidate.data;
    await assertSafeObjectKeysForImport(data, 'Session');
    if (isLegacyConfigPayload(data)) {
      // A legacy configuration names no mode: it is the shown mode's settings.
      applyLegacyConfigPayload(state.drawings[state.mode.value === 'linear' ? 'linear' : 'circular'], data);
      alert('Legacy configuration loaded. Save as a session to use the current format.');
      return { status: 'legacy' };
    }

    recordSessionLifecycleEvent('current-session-preflight-start');
    const preflight = await preflightSessionImport(data);
    await validateSimilarityAlignmentResetReceipt(
      data.editorState?.alignmentResetReceipt,
      { renderRequest: data.renderRequest, resources: data.resources }
    );
    if (isPlainObject(data.otherModeResult)) {
      await validateSimilarityAlignmentResetReceipt(
        data.otherModeResult.editorState?.alignmentResetReceipt,
        { renderRequest: data.otherModeResult.renderRequest, resources: data.resources }
      );
    }
    recordSessionLifecycleEvent('current-session-preflight-end');
    data = preflight.data;
    const {
      sourceSessionVersion,
      canonicalProjection,
      restoredConfig,
      projectionResult,
      adoptedCanonicalSession,
      currentResourceTable,
      otherModeCatalog,
      comparisonClassification,
      unmanagedConfigValidation
    } = preflight;
    if (restoredConfig && unmanagedConfigValidation) {
      restoredConfig.unmanagedConfigOverrides = await validateUnmanagedConfigOverrides(
        unmanagedConfigValidation
      );
    }
    const canonicalSession = projectionResult !== null;
    const settingsOnly = isSettingsOnlySessionDocument(data);
    const currentSchemaSession = sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION;
    const committedMode = projectionResult?.renderState.mode;
    const savedCurrentWriterMode = sourceSessionVersion >= TYPED_DRAFT_SESSION_VERSION
      ? data.ui?.mode
      : null;
    const ui = canonicalSession
      ? {
          ...projectionResult.editorMetadata.ui,
          ...projectionResult.artifactState.ui,
          mode: ['circular', 'linear'].includes(savedCurrentWriterMode)
            ? savedCurrentWriterMode
            : committedMode
        }
      : (data.ui || {});
    const candidateMode = ui.mode === 'linear' ? 'linear' : 'circular';
    // PD-OI-086: each mode's drawing gets its own draft. A Session 46 keeps it
    // in its slice; the committed mode's slice lies over its request's
    // projection (`restoredConfig`, plan 4.1). A Session 27-44 has one flat
    // draft, split below once its edits are migrated.
    const modeScopedSession = sourceSessionVersion >= MODE_SCOPED_SESSION_VERSION;
    const committedDraftMode = committedMode === 'linear' || (!committedMode && candidateMode === 'linear')
      ? 'linear' : 'circular';
    // A Session 46 committed mode with a Result: its slice over its request's
    // projection (`restoredConfig`); any other slice over its mode's defaults.
    const committedSliceProjected = modeScopedSession && canonicalSession && !settingsOnly && Boolean(restoredConfig);
    /** @param {'circular' | 'linear'} mode @returns {Record<string, any> | null} */
    const draftConfigOf = (mode) => (modeScopedSession
      ? (mode === committedDraftMode && committedSliceProjected ? restoredConfig : overlayModeSliceConfig(null, data.modes?.[mode]?.config))
      : restoredConfig);
    const linearDraftConfig = draftConfigOf('linear');
    const circularDraftConfig = draftConfigOf('circular');
    const candidateInputType = (candidateMode === 'linear'
      ? ui.lInputType : ui.cInputType) || canonicalProjection?.inputType || 'gb';
    const candidateFiles = {
      files: {}, linearSeqs: [], circularRecordList: { value: [] },
      circularRecordDiscovery: {}, linearReorderNotice: { value: '' },
      cInputType: { value: candidateMode === 'circular' ? candidateInputType : (ui.cInputType || 'gb') },
      linearRecordRows: cloneJsonData(linearDraftConfig?.linearRecordLayout?.rows || []),
      linearComparisonPlan: normalizeLinearComparisonPlan(linearDraftConfig?.linearComparisonPlan)
    };
    recordSessionLifecycleEvent('session-candidate-files-start');
    const { collapsedLinearSeqs } = applyFiles(
      canonicalSession ? projectionResult.restoredFiles : data.files,
      {
        adoptCanonicalPayloads: currentSchemaSession, resolveRecordInputs: !settingsOnly,
        targetState: candidateFiles, targetDrawing: candidateFiles
      }
    );
    recordSessionLifecycleEvent('session-candidate-files-end');
    const importedResults = canonicalSession
      ? projectionResult.artifactState.results
      : data.results;
    const logicalImportedResults = normalizeLogicalResults(
      (Array.isArray(importedResults) ? importedResults : []).map((result, index) => ({
        name: result?.name || `Result ${index + 1}`,
        content: result?.content || ''
      }))
    );
    const storedEditorState = canonicalSession
      ? projectionResult.artifactState.editorState
      : data.editorState;
    /** @type {ReturnType<typeof featureStateFromCatalog> | null} */
    let currentCatalogFeatureState = null;
    let validatedSessionCatalog = currentSchemaSession
      ? projectionResult?.validatedFeatureCatalog || null
      : null;
    let restoredEditorState = storedEditorState;
    if (currentSchemaSession && validatedSessionCatalog) {
      currentCatalogFeatureState = featureStateFromCatalog(
        validatedSessionCatalog,
        { mode: committedMode }
      );
    }

    const artifactFeatureState = canonicalSession
      ? currentSchemaSession
        ? { ...projectionResult.artifactState.features }
        : cloneJsonData(projectionResult.artifactState.features)
      : {};
    if (currentSchemaSession) {
      [
        'extractedFeatures',
        'biologicalFeatures',
        'featureSelectorSafetyScope',
        'featureRecordIds'
      ].forEach((field) => delete artifactFeatureState[field]);
    }
    let features = splitLegacyFeatureVisibilityRules(canonicalSession
      ? {
          ...projectionResult.renderState.semanticFeatureState,
          ...(currentCatalogFeatureState || {}),
          ...artifactFeatureState
        }
      : (data.features || {}));
    const catalogSequenceSources = currentSchemaSession
      ? (currentCatalogFeatureState?.sequenceSources || [])
      : [];
    recordSessionLifecycleEvent('session-candidate-sequences-start');
    // E1: each Result set reads its sources in its own mode.
    const setSequenceSourceOptions = {
      olderSession: sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION,
      settingsOnly,
      cInputType: candidateFiles.cInputType.value,
      lInputType: candidateMode === 'linear' ? candidateInputType : (ui.lInputType || 'gb'),
      files: candidateFiles.files,
      linearSeqs: candidateFiles.linearSeqs,
      circularConservation: circularDraftConfig?.circularConservation || {}
    };
    const { restored: restoredFileSequenceSources, error: currentRecoveryError } = await restoreLoadedSetSequenceSources({
      ...setSequenceSourceOptions,
      mode: committedMode === 'linear' || (!committedMode && candidateMode === 'linear') ? 'linear' : 'circular',
      catalog: currentSchemaSession ? validatedSessionCatalog : null,
      renderRequest: data.renderRequest
    });

    recordSessionLifecycleEvent('session-candidate-sequences-end');
    const restoredFeatureState = currentCatalogFeatureState || features || {};
    // A saved Result is laid out at its committed legend and title sides.
    const committedLayout = canonicalProjection?.layoutPreferences
      ? resolveActiveLayoutPreference(
          canonicalProjection.layoutPreferences,
          canonicalProjection.mode,
          Boolean(canonicalProjection.config?.form?.multi_record_canvas)
        )
      : null;
    // The composition and stroke owners transform an older Session's Result
    // through the composition root's port (R13); this module passes data only.
    /** @type {LegacyResultSvgTransform} */
    const transformLegacyResultSvg = options.transformLegacyResultSvg;
    const transformRestoredSessionSvg = (svg) => {
      const legendGroupsChanged = normalizeLegacyLegendEntryGroups(svg);
      const legacyResultChanged = transformLegacyResultSvg(svg, {
        composition: {
          legendSide: committedLayout?.legend || restoredConfig?.form?.legend || 'none',
          titleSide: committedLayout?.plotTitlePosition
            || restoredConfig?.adv?.plot_title_position || 'none',
          userDeltas: {
            primary: ui.diagramOffset ? [ui.diagramOffset.x, ui.diagramOffset.y] : null,
            legend: ui.legendCurrentOffset
              ? [ui.legendCurrentOffset.x, ui.legendCurrentOffset.y]
              : null,
            lengthBar: ui.lengthBarUserOffset
              ? [ui.lengthBarUserOffset.x, ui.lengthBarUserOffset.y]
              : null,
            title: ui.plotTitleUserOffset
              ? [ui.plotTitleUserOffset.x, ui.plotTitleUserOffset.y]
              : null
          }
        },
        strokes: {
          features: restoredFeatureState.extractedFeatures || [],
          legendStrokeOverrides: restoredEditorState?.legend?.strokeOverrides || {},
          featureStrokeOverrides: restoredEditorState?.featureStrokes?.overrides || {}
        }
      });
      return legendGroupsChanged || legacyResultChanged;
    };

    recordSessionLifecycleEvent('svg-admission-start');
    // The loaded Results carry the feature types of the request that becomes
    // the committed request below (R13). An older Session adopts no catalog,
    // so its Features list lists what each Result renders and reads none.
    const selectedFeatureTypes = currentSchemaSession
      ? adoptedCanonicalSession?.renderRequest?.diagramOptions?.selectedFeaturesSet
      : null;
    const committedImportedResults = currentSchemaSession && validatedSessionCatalog
      ? admitLoadedSetResults(logicalImportedResults, {
          featureCatalog: validatedSessionCatalog, mode: committedMode, selectedFeatureTypes
        })
      : admitLegacyImportedResults(
          createLegacyImportResultSource(logicalImportedResults),
          { transformSvg: transformRestoredSessionSvg, selectedFeatureTypes }
        );
    recordSessionLifecycleEvent('svg-admission-end');
    // E1: a Session with a Result set of each mode builds both sets' artifact
    // slots with one function. It opens on its shown mode when that mode has a
    // Result, otherwise on the mode that has one; the opening slot is installed
    // last, before the preview mounts, and the other slot waits in its mode.
    const otherModeCatalogAdmitted = currentSchemaSession ? otherModeCatalog : null;
    const otherModeResult = otherModeCatalogAdmitted ? data.otherModeResult : null;
    const otherSetMode = otherModeResult?.renderRequest?.mode === 'linear' ? 'linear' : 'circular';
    const resultModes = new Set([
      ...(committedImportedResults.length > 0 && committedMode ? [committedMode] : []),
      ...(otherModeResult ? [otherSetMode] : [])
    ]);
    const displayMode = resultModes.size === 0 || resultModes.has(candidateMode)
      ? candidateMode
      : [...resultModes][0];
    // Each set's History byte estimate: its share of the decompressed file,
    // split by the size of its Results (a serialization of the other set only
    // to size it costs about 20 ms for a Vibrio-size set).
    /** @param {unknown} setResults */
    const resultCharacters = (setResults) => (Array.isArray(setResults) ? setResults : [])
      .reduce((sum, result) => sum + String(result?.content || '').length, 0);
    const otherResultCharacters = otherModeResult ? resultCharacters(otherModeResult.results) : 0;
    const otherSetCharacters = otherModeResult ? Math.round(candidate.characters * otherResultCharacters
      / Math.max(1, otherResultCharacters + resultCharacters(data.results))) : 0;
    const topSetCharacters = Math.max(0, candidate.characters - otherSetCharacters);
    const slotOptions = { resources: data.resources, resourceTable: currentResourceTable, webFiles: data.webFiles };
    const otherArtifactSlot = otherModeResult && otherModeCatalogAdmitted
      ? buildLoadedArtifactSlot(otherModeResult, {
          ...slotOptions,
          featureCatalog: otherModeCatalogAdmitted,
          restoredSequenceSources: (await restoreLoadedSetSequenceSources({
            ...setSequenceSourceOptions,
            mode: otherSetMode,
            catalog: otherModeCatalogAdmitted,
            renderRequest: otherModeResult.renderRequest
          })).restored,
          retainedBytes: otherSetCharacters * 2
        })
      : null;
    const opensOtherSet = Boolean(otherArtifactSlot) && displayMode === otherSetMode;
    // The Linear set's committed request owns the Linear draft's alignment plan
    // and record translations, wherever the set sits; both are admitted as the
    // request projection admits the top-level set's.
    const linearSetLayout = otherModeResult && otherSetMode === 'linear'
      ? canonicalLinearRecordLayout(otherModeResult.renderRequest) : null;

    /** @type {Record<string, any> | null} */
    let legacyFeatureRecoveryPlan = null;
    const legacyFeatureSnapshot = {
      mode: candidateMode, cInputType: candidateFiles.cInputType.value,
      lInputType: candidateInputType, files: candidateFiles.files,
      linearSeqs: candidateFiles.linearSeqs, results: committedImportedResults,
      selectedResultIndex: ui.selectedResultIndex || 0,
      featureState: features, editorState: restoredEditorState
    };
    if (sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION) {
      try {
        legacyFeatureRecoveryPlan = await buildSessionFeatureRecoveryPlan(/** @type {any} */ ({
          snapshot: legacyFeatureSnapshot,
          featureVisibilityTsv: serializeFeatureVisibilityRules(
            features.featureVisibilityManualRules || features.featureVisibilityRules || []
          )
        }));
      } catch (error) {
        legacyFeatureRecoveryPlan = { status: 'failed', warning: normalizeCaughtError(error).summary };
      }
    }
    // Session 46 keys per-feature edits by source identity: an older Session's
    // saved rendered-ID edits are mapped through its saved catalog, or without
    // one through its sources read again with its crops and orientations (the
    // saved metadata when they cannot be read), once (design Q4 4.3).
    let droppedFeatureEditCount = 0;
    let narrowedFeatureVisibilityCount = 0;
    let migratedAnnotationTargetCount = 0;
    if (sourceSessionVersion < SESSION_VERSION) {
      const recovered = legacyFeatureRecoveryPlan?.recoveredFeatureState;
      const sourceFeatures = !validatedSessionCatalog && hasRenderedIdFeatureEdits(features)
        ? await extractSessionSourceFeatures(/** @type {any} */ ({ snapshot: legacyFeatureSnapshot }))
        : null;
      const migration = migrateSessionFeatureEdits({
        features,
        mode: data.renderRequest?.mode || candidateMode,
        catalog: validatedSessionCatalog,
        legacy: validatedSessionCatalog ? null : {
          records: data.renderRequest?.records,
          features: [
            sourceFeatures?.extractedFeatures, recovered?.biologicalFeatures, recovered?.extractedFeatures,
            features.biologicalFeatures, features.extractedFeatures
          ].find((list) => Array.isArray(list) && list.length > 0) || [],
          biologicalFeatures: sourceFeatures?.biologicalFeatures || []
        }
      });
      features = migration.features;
      droppedFeatureEditCount = migration.droppedCount;
      narrowedFeatureVisibilityCount = migration.narrowedVisibilityCount;
      // R-7: an annotation's `hash=` target moves to its source feature only
      // where the saved catalog makes the figure certain.
      if (restoredConfig) {
        const annotationMigration = migrateSessionAnnotationTargets({
          annotationSets: restoredConfig.annotationSets,
          mode: data.renderRequest?.mode || candidateMode,
          catalog: validatedSessionCatalog,
          records: data.renderRequest?.records
        });
        migratedAnnotationTargetCount = annotationMigration.migratedCount;
        if (migratedAnnotationTargetCount > 0) restoredConfig.annotationSets = annotationMigration.annotationSets;
      }
      if (recovered) {
        const recoveredWithoutRenderedIdEdits = { ...recovered };
        RENDERED_ID_FEATURE_EDIT_FIELDS.forEach((field) => delete recoveredWithoutRenderedIdEdits[field]);
        // `recovered` came from legacyFeatureRecoveryPlan.recoveredFeatureState, so the plan is set.
        // Its edit rows are the committed mode's drawing's (R1-3).
        /** @type {Record<string, any>} */ (legacyFeatureRecoveryPlan).recoveredFeatureState = {
          ...recoveredWithoutRenderedIdEdits,
          featureOverrides: unscopedDraftRows(features.featureOverrides, committedDraftMode),
          labelOverrideRows: features.labelOverrideRows
        };
      }
    }
    // Each drawing's draft (PD-OI-086, plan 4.1). A Session 46 slice lies over
    // its mode's committed projection. The flat draft of a Session 27-44, as
    // Load has read it so far, is split by the registry (plan 4.2); a Session
    // written without a draft (CLI, Python) splits only what its migrations
    // wrote, and its committed mode keeps its request's projection.
    const hasFlatDraft = !modeScopedSession && isPlainObject(data.config);
    /** @type {Record<'circular' | 'linear', Record<string, any>>} */
    let modeSlices = { circular: {}, linear: {} };
    /** @type {Record<string, any>} */
    let appUi = modeScopedSession && isPlainObject(data.ui) ? data.ui : {};
    let sessionCliOptions = modeScopedSession ? data.cliOptions : undefined;
    if (modeScopedSession) {
      modeSlices = {
        circular: isPlainObject(data.modes?.circular) ? data.modes.circular : {},
        linear: isPlainObject(data.modes?.linear) ? data.modes.linear : {}
      };
    } else {
      const flatConfig = hasFlatDraft
        ? restoredConfig
        : migratedAnnotationTargetCount > 0 && restoredConfig ? { annotationSets: restoredConfig.annotationSets } : {};
      const legend = isPlainObject(restoredEditorState?.legend) ? restoredEditorState.legend : {};
      // The split writes what the Session saved and its migrations wrote
      // (plan 4.2), not the empty defaults Load reads in their place.
      const savedDraft = (/** @type {unknown} */ value) => isPlainObject(value);
      /**
       * @param {Record<string, any> | null | undefined} value
       * @param {unknown} saved
       * @param {string[]} keys
       */
      const withoutUnsavedEmpty = (value, saved, keys) => {
        if (!value || !isPlainObject(value)) return value;
        const draft = /** @type {Record<string, any>} */ (value);
        const unsaved = keys.filter((key) => {
          if (!Object.hasOwn(draft, key) || (isPlainObject(saved) && Object.hasOwn(/** @type {Record<string, any>} */ (saved), key))) return false;
          const field = draft[key];
          return isPlainObject(field) ? Object.keys(field).length === 0
            : Array.isArray(field) ? field.length === 0 : field === 0;
        });
        return Object.fromEntries(Object.entries(draft).filter(([key]) => !unsaved.includes(key)));
      };
      const split = splitDraftIntoModes({
        config: withoutUnsavedEmpty(flatConfig, data.config, ['unmanagedConfigOverrides']) || {},
        // With a catalog, the feature-edit migration writes `featureOverrides`.
        ...(savedDraft(data.features) ? { features: withoutUnsavedEmpty(features, data.features, [
          'featureColorOverrides', 'featureVisibilityManualRules', 'labelOverrideRows', 'labelTextBulkOverrides',
          'selectedFeatureRecordIdx', ...(validatedSessionCatalog ? [] : ['featureOverrides'])
        ]) } : {}),
        editorState: {
          ...(savedDraft(data.editorState?.legend) ? { legend } : {}),
          ...(savedDraft(data.editorState?.featureStrokes)
            ? { featureStrokes: restoredEditorState?.featureStrokes || { overrides: {} } } : {})
        },
        // A saved layout goes to its slots (registry `layout`); legacy layout
        // fields are restored below, as before.
        ui,
        renderRequest: data.renderRequest
      }, {
        committedMode: committedDraftMode,
        depthSources: sessionDepthSourceWidths({
          c_depth: candidateFiles.files.c_depth, linearSeqs: candidateFiles.linearSeqs
        }),
        paletteColors: paletteColorsFromDefinitions(String(flatConfig?.palette || '').trim() || 'default')
      });
      // The shared migration vectors read the split as Load made it
      // (tests/web/mode-split-vectors.test.mjs).
      recordSessionLifecycleEvent('session-draft-split', { split });
      modeSlices = split.modes || modeSlices;
      appUi = split.ui || {};
      sessionCliOptions = split.cliOptions;
    }
    // The draft each drawing installs: its slice over its committed projection.
    /** @type {Record<'circular' | 'linear', Record<string, any>>} */
    const modeConfigs = {
      circular: {}, linear: {}
    };
    for (const mode of SLICE_MODES) {
      const committed = mode === committedDraftMode;
      modeConfigs[mode] = modeScopedSession
        ? (committed && committedSliceProjected ? restoredConfig : overlayModeSliceConfig(null, modeSlices[mode].config))
        : overlayModeSliceConfig(committed && !hasFlatDraft ? restoredConfig : null, modeSlices[mode].config);
      // A saved unmanaged override is checked against its own drawing's mode (PD-OI-009, OV-106).
      const savedOverrides = modeSlices[mode].config?.unmanagedConfigOverrides;
      if (isPlainObject(savedOverrides) && Object.keys(savedOverrides).length > 0) {
        modeConfigs[mode].unmanagedConfigOverrides = await validateUnmanagedConfigOverrides({
          mode, configOverrides: savedOverrides, requireUnmanagedOnly: true
        });
      }
    }
    recordSessionLifecycleEvent('session-candidate-prepared');
    // The admitted candidate is still private. Let the browser handle input
    // after SVG sanitation, before the one atomic live-state transaction.
    await new Promise((resolve) => setTimeout(resolve, 0));
    if (!options.isCurrent()) return { status: 'canceled' };
    rollbackSnapshot = captureSessionImportSnapshot();
    if (typeof rollbackStateExtension?.capture === 'function') {
      rollbackExtensionSnapshot = rollbackStateExtension.capture();
    }
    commitStarted = true;
    state.semanticFileWatchersSuppressed.value = true;
    resetSessionBaseline();
    state.sessionResourceDiscoveryDeferred.value = currentSchemaSession;
    if (currentSchemaSession) {
      committedCanonicalSession = adoptedCanonicalSession;
      activeSessionResourceTable = currentResourceTable;
    }
    state.sessionTitle.value = canonicalSession
      ? projectionResult.documentMetadata.title
      : (typeof data.title === 'string' ? data.title : '');
    if (ui.mode) state.mode.value = ui.mode;
    if (canonicalProjection) {
      if (canonicalProjection.mode === 'circular') {
        state.cInputType.value = canonicalProjection.inputType;
      } else {
        state.lInputType.value = canonicalProjection.inputType;
      }
    }
    if (ui.cInputType) state.cInputType.value = ui.cInputType;
    if (ui.lInputType) state.lInputType.value = ui.lInputType;
    if (ui.downloadDpi) state.downloadDpi.value = ui.downloadDpi;
    // Suppressed mode watchers observe only the complete admitted document.
    state.autoLabelReflowEnabled.value = Boolean(ui.autoLabelReflow);
    // A Session 27-44 draft's boolean Instant Preview wins over its `ui` (plan 4.2).
    state.paletteInstantPreviewEnabled.value = Boolean(
      typeof appUi.paletteInstantPreviewEnabled === 'boolean'
        ? appUi.paletteInstantPreviewEnabled : ui.paletteInstantPreviewEnabled
    );
    if (ui.featurePanelTab === 'labels' || ui.featurePanelTab === 'colors') {
      state.featurePanelTab.value = ui.featurePanelTab;
    } else {
      state.featurePanelTab.value = 'colors';
    }
    state.generatedMode.value = (committedMode || ui.mode) === 'linear' ? 'linear' : 'circular';
    if (ui.generatedLegendPosition) {
      state.generatedLegendPosition.value = normalizeLegendPosition(
        ui.generatedLegendPosition,
        (committedMode || ui.mode) === 'linear' ? 'bottom' : 'left'
      );
    }
    state.generatedMultiRecordCanvas.value = Boolean(ui.generatedMultiRecordCanvas);
    state.generatedCircularPlotTitlePosition.value = hasStoredLayoutValue(ui.generatedCircularPlotTitlePosition)
      ? normalizeCircularPlotTitlePosition(ui.generatedCircularPlotTitlePosition)
      : normalizeCircularPlotTitlePosition(ui.circularPlotTitlePosition);

    state.suppressCircularMultiRecordDefaults.value = shouldSuppressCircularMultiRecordDefaults(
      state.drawings.circular,
      modeConfigs.circular.form
    );
    // Each drawing over its mode's defaults (resetSessionBaseline); the shown
    // drawing's Canvas padding follows once its Result is mounted. A Session 46
    // slice without its mode's request projection is stored settings: its track
    // stacks install as saved, unset axis indexes included (R11).
    for (const mode of SLICE_MODES) {
      const storedSlice = modeScopedSession && !(mode === committedDraftMode && committedSliceProjected);
      applyModeSliceData(state.drawings[mode], mode, { ...modeSlices[mode], config: modeConfigs[mode] }, {
        resolveTrackPlacements: !settingsOnly && !storedSlice,
        applyCanvasPadding: mode !== displayMode
      });
    }
    // App-level settings: LOSAT execution (after the drafts, whose request
    // projection may name a thread count), the rich popup, CLI provenance.
    if (isPlainObject(appUi.losatExecution)) applyLosatExecutionData(appUi.losatExecution);
    state.richFeaturePopup.value = appUi.richFeaturePopup !== false;
    preservedCliOptions = isPlainObject(sessionCliOptions) ? cloneJsonData(sessionCliOptions) : null;
    state.mode.value = displayMode;
    const shownDrawing = state.drawings[displayMode];
    const canonicalLinearLayout = linearSetLayout || canonicalProjection?.config?.linearRecordLayout;
    if (canonicalLinearLayout && state.linearRecordTranslations) {
      state.linearRecordTranslations.value = cloneJsonData(
        canonicalLinearLayout.recordTranslations || []
      );
    }
    if (canonicalLinearLayout && state.similarityAlignmentPlan) {
      state.similarityAlignmentPlan.value = canonicalLinearLayout.similarityAlignment
        ? cloneJsonData(canonicalLinearLayout.similarityAlignment)
        : null;
    }
    if (state.legacySimilarityAlignment) {
      state.legacySimilarityAlignment.value = canonicalSession
        ? cloneJsonData(
            projectionResult.artifactState.legacySimilarityAlignment || null
          )
        : null;
    }
    restoreAppliedPaletteFromSession(shownDrawing, ui);
    // Layout preferences: a Session 46 slice holds its mode's slot (the
    // committed mode's request sets it when the slice has none); a Session
    // 27-44 keeps one set, read once, which both drawings take.
    const projectedLayout = canonicalSession ? canonicalProjection?.layoutPreferences : null;
    if (modeScopedSession) {
      if (projectedLayout && !isPlainObject(modeSlices[committedDraftMode].ui?.layoutPreferences)) {
        const drawingOfCommit = state.drawings[committedDraftMode];
        replaceLayoutPreferences(drawingOfCommit.layoutPreferences, {
          ...cloneJsonData(drawingOfCommit.layoutPreferences),
          [committedDraftMode]: normalizeLayoutPreferences(projectedLayout)[committedDraftMode]
        });
      }
    } else {
      SLICE_MODES.forEach((mode) => restoreLayoutPreferences(state.drawings[mode], ui, { projected: projectedLayout }));
    }

    restoreLiveFileState(state.drawings.linear, {
      files: candidateFiles.files,
      linearSeqs: candidateFiles.linearSeqs,
      circularRecordList: candidateFiles.circularRecordList.value,
      circularRecordDiscovery: candidateFiles.circularRecordDiscovery,
      linearRecordRows: candidateFiles.linearRecordRows,
      linearComparisonPlan: candidateFiles.linearComparisonPlan
    });
    // The comparisons are the Linear drawing's.
    restoreImportedComparisonIntent(
      state.drawings.linear.importedComparisonIntent,
      comparisonClassification,
      linearDraftConfig?.importedComparisonResolution
    );
    if (!settingsOnly) {
      SLICE_MODES.forEach((mode) => reconcileDepthTrackStateAfterSessionFiles(state.drawings[mode], mode));
    }
    if (canonicalSession) {
      applyLosatCache(
        state.drawings.linear,
        projectionResult.artifactState.losatCache?.entries,
        projectionResult.artifactState.legacyArtifacts?.proteinRawCandidates,
        {
          adoptCurrent: currentSchemaSession,
          validatedManifest: projectionResult.artifactState.proteinIdentityManifest
        }
      );
      applyLosatDerivedCache(
        projectionResult.artifactState.losatDerivedCache?.entries,
        projectionResult.artifactState.legacyArtifacts?.proteinDerivedEvidence,
        { adoptCurrent: currentSchemaSession }
      );
      applyProteinIdentityManifest(
        projectionResult.artifactState.proteinIdentityManifest,
        { adoptCurrent: currentSchemaSession }
      );
    } else {
      applyLosatCache(
        state.drawings.linear,
        data.losatCache?.entries,
        data.legacyArtifacts?.proteinRawCandidates
      );
      applyLosatDerivedCache(
        data.losatDerivedCache?.entries,
        data.legacyArtifacts?.proteinDerivedEvidence
      );
      applyProteinIdentityManifest(data.proteinIdentityManifest);
    }
    if (collapsedLinearSeqs) {
      state.losatCacheInfo.value = [];
    }

    state.skipCaptureBaseConfig.value = false;

    // The drawings hold their edits; the Result's catalog features are artifacts.
    applyFeatureArtifactData(features);
    state.selectedFeatureRecordIdx.value = Number.isInteger(modeSlices[displayMode].ui?.selectedFeatureRecordIdx)
      ? modeSlices[displayMode].ui.selectedFeatureRecordIdx : 0;
    if (currentSchemaSession && currentCatalogFeatureState) {
      state.collinearGroups.value = currentCatalogFeatureState.collinearGroups;
      synchronizeRestoredFeatureSummaryStatus({ generationId: 'session-load' });
    }
    state.matchSequenceRegistry?.reset?.([
      ...catalogSequenceSources,
      ...restoredFileSequenceSources
    ]);

    // Group names are the Linear drawing's (Session 46 `webEdits`).
    applyOrthogroupStateData(
      state.drawings.linear,
      canonicalSession
        ? {
            ...projectionResult.artifactState.orthogroupState,
            ...(modeScopedSession ? {
              orthogroupNameOverrides: modeConfigs.linear.webEdits?.orthogroupNameOverrides || {},
              orthogroupDescriptionOverrides: modeConfigs.linear.webEdits?.orthogroupDescriptionOverrides || {}
            } : {}),
            ...(currentCatalogFeatureState
              ? { groups: currentCatalogFeatureState.orthogroups }
              : {})
          }
        : data.orthogroupState && typeof data.orthogroupState === 'object'
          ? data.orthogroupState
        : {
            groups: Array.isArray(data.orthogroups) ? data.orthogroups : [],
            selectedOrthogroupId: features.selectedOrthogroupId,
            selectedOrthogroupAlignmentFeature: features.selectedOrthogroupAlignmentFeature,
            orthogroupNameOverrides:
              features.orthogroupNameOverrides ||
              data.config?.webEdits?.orthogroupNameOverrides ||
              {},
            orthogroupDescriptionOverrides:
              features.orthogroupDescriptionOverrides ||
              data.config?.webEdits?.orthogroupDescriptionOverrides ||
              {}
        },
      {
        legacyRecords: sourceSessionVersion <= 39 ? data.renderRequest?.records : null,
        catalogFeatureState: currentCatalogFeatureState
      }
    );
    applyEditorArtifactData(currentSchemaSession ? restoredEditorState : normalizeEditorStateData(restoredEditorState));
    if (legacyFeatureRecoveryPlan) {
      applySessionFeatureRecoveryPlan(state.drawings[committedDraftMode], legacyFeatureRecoveryPlan, { generationId: 'session-load' });
    }

    let desiredResultIndex = (
      Number.isInteger(ui.selectedResultIndex) && ui.selectedResultIndex >= 0
    )
      ? Math.min(ui.selectedResultIndex, Math.max(0, committedImportedResults.length - 1))
      : 0;
    if (!options.isCurrent()) throw new Error('Session loading was canceled.');
    state.skipCaptureBaseConfig.value = true;
    recordSessionLifecycleEvent('session-candidate-adopted');
    recordSessionLifecycleEvent('preview-mount-start');
    applyResultsData(committedImportedResults, ui);
    state.annotationWarnings.value = cloneJsonData(
      projectionResult?.artifactState?.runMetadata?.annotationWarnings || []
    );
    state.featureIdentityNotices.value = cloneJsonData(
      projectionResult?.artifactState?.runMetadata?.featureIdentityNotices || []
    );
    state.comparisonWarnings.value = cloneJsonData(
      projectionResult?.artifactState?.runMetadata?.comparisonWarnings || []
    );
    state.trackSlotResolvedGeometry.value = cloneJsonData(
      projectionResult?.artifactState?.runMetadata?.trackSlotGeometry ?? null
    );
    {
      // The top-level set, built like the other set, waits in its mode when the
      // Session opens on the other set's mode. The stash starts again with the
      // loaded Session.
      const topArtifactSlot = opensOtherSet && validatedSessionCatalog
        ? buildLoadedArtifactSlot({
            renderRequest: data.renderRequest,
            results: data.results,
            editorState: data.editorState,
            ui,
            runMetadata: projectionResult?.artifactState?.runMetadata
          }, {
            ...slotOptions,
            featureCatalog: validatedSessionCatalog,
            admittedResults: committedImportedResults,
            restoredSequenceSources: restoredFileSequenceSources,
            retainedBytes: topSetCharacters * 2
          })
        : null;
      /** @type {{ installLoadedArtifactSlots?: (slots: Record<string, any>) => void }} */ (options)
        .installLoadedArtifactSlots?.(opensOtherSet
        ? { opening: otherArtifactSlot, stashed: topArtifactSlot }
        : { opening: null, stashed: otherArtifactSlot });
      if (opensOtherSet) desiredResultIndex = Number(state.selectedResultIndex.value) || 0;
    }
    const previewReadiness = typeof options?.beforePreviewMount === 'function'
      ? options.beforePreviewMount({
          results: state.results.value,
          resultIndex: desiredResultIndex,
          data,
          ui
        })
      : null;
    recordSessionLifecycleEvent('firstCommittedPreview', {
      resultCount: state.results.value.length
    });
    await nextTick();
    recordSessionLifecycleEvent('preview-mount-end');
    if (previewReadiness?.promise) await previewReadiness.promise;
    else if (previewReadiness?.then) await previewReadiness;

    if (typeof options?.afterLoad === 'function') {
      await options.afterLoad({ data, ui });
    }

    const shownPadding = modeSlices[displayMode].ui?.canvasPadding;
    if (isPlainObject(shownPadding)) {
      shownDrawing.canvasPadding.top = shownPadding.top || 0;
      shownDrawing.canvasPadding.right = shownPadding.right || 0;
      shownDrawing.canvasPadding.bottom = shownPadding.bottom || 0;
      shownDrawing.canvasPadding.left = shownPadding.left || 0;
    }
    if (ui.canvasPan) {
      state.canvasPan.x = ui.canvasPan.x || 0;
      state.canvasPan.y = ui.canvasPan.y || 0;
    }
    if (typeof ui.zoom === 'number') {
      state.zoom.value = ui.zoom;
    }

    state.semanticFileWatchersSuppressed.value =
      semanticFileWatchersSuppressedBeforeImport;
    await nextTick();
    if (!currentSchemaSession) {
      committedCanonicalSession = cloneCanonicalSession(data);
      activeSessionResourceTable = null;
    }
    recordSessionLifecycleEvent('interactiveReady', {
      status: 'success',
      degradedRecovery: Boolean(currentRecoveryError)
    });
    if (!options.isCurrent()) throw new Error('Session loading was canceled.');
    await options.afterImport?.({
      status: 'ok',
      decompressedCharacters: opensOtherSet ? otherSetCharacters : topSetCharacters,
      isCurrent: options.isCurrent
    });
    if (!options.isCurrent()) throw new Error('Session loading was canceled.');
    alert(['Session loaded successfully!',
      droppedFeatureEditCount > 0 ? FEATURE_EDIT_MIGRATION_WARNING(droppedFeatureEditCount) : '',
      narrowedFeatureVisibilityCount > 0 ? FEATURE_VISIBILITY_NARROWED_NOTICE(narrowedFeatureVisibilityCount) : '',
      migratedAnnotationTargetCount > 0 ? ANNOTATION_TARGET_MIGRATION_NOTICE(migratedAnnotationTargetCount) : '',
      legacyTableRowsNotice(canonicalProjection?.legacyTableRepairs)
    ].filter(Boolean).join(' '));
    return {
      status: 'ok',
      data,
      decompressedCharacters: candidate.characters,
      degradedRecovery: Boolean(currentRecoveryError),
      comparisonDisposition: state.drawings.linear.importedComparisonIntent.disposition
    };
  } catch (err) {
    if (err?.name === 'AbortError' && !commitStarted) return { status: 'canceled' };
    const error = normalizeCaughtError(err, { stage: 'request-validation' });
    const currentAlert = state.errorLog.value;
    const canNotify = currentAlert === previousAlert || currentAlert === null;
    if (commitStarted && rollbackSnapshot) {
      try {
        await restoreSessionImportSnapshot(rollbackSnapshot);
        if (typeof rollbackStateExtension?.restore === 'function') {
          await rollbackStateExtension.restore(rollbackExtensionSnapshot);
        }
      } catch (rollbackError) {
        state.generationFailureRecovery.value = 'restore-failed';
      }
    }
    const restoredPreviousAlert = JSON.stringify(normalizeUserFacingError(state.errorLog.value))
      === JSON.stringify(normalizeUserFacingError(previousAlert));
    if (!canNotify) {
      if (restoredPreviousAlert) state.errorLog.value = currentAlert;
      return { status: 'stale' };
    }
    if (state.errorLog.value !== previousAlert && state.errorLog.value !== null && !restoredPreviousAlert) return { status: 'stale' };
    state.errorLog.value = error;
    recordSessionLifecycleEvent('interactiveReady', {
      status: 'error',
      error: error.code
    });
    return { status: 'error', error };
  } finally {
    state.semanticFileWatchersSuppressed.value =
      semanticFileWatchersSuppressedBeforeImport;
    e.target.value = '';
  }
};


/** @type {{ canceled: boolean, promise: Promise<any> | null } | null} */
let sessionSaveInFlight = null;
/** @type {{ canceled: boolean, controller: AbortController } | null} */
let activeSessionImport = null;

export const disposeSessionOperations = () => {
  if (sessionSaveInFlight) sessionSaveInFlight.canceled = true;
  if (activeSessionImport) {
    activeSessionImport.canceled = true;
    activeSessionImport.controller.abort();
  }
  sessionSaveInFlight = null;
  activeSessionImport = null;
  state.sessionSavePending.value = false;
  state.sessionImportPending.value = false;
};

// `options.availability` is the Save and Load availability the composition
// root composes from `sessionOperationAvailability` and the edits still
// applying (R13); `options.recordDisplayRows` are the draft record display
// rows a Save without a committed request projects;
// `options.readOtherModeArtifact` returns the other mode's artifact slot that
// the root keeps (E1), and `options.beforeExport` learns whether Save projects
// the draft request.
export const exportSession = (titleOverride = null, options = {}) => {
  if (sessionSaveInFlight) {
    recordSessionLifecycleEvent('session-save-joined');
    return sessionSaveInFlight.promise;
  }
  const { availability = sessionOperationAvailability } = options;
  const busy = availability('save');
  if (busy) return Promise.resolve(busy);
  /** @type {{ canceled: boolean, promise: Promise<any> | null }} */
  const operation = { canceled: false, promise: null };
  const previousAlert = state.errorLog.value;
  const isCurrent = () => sessionSaveInFlight === operation && !operation.canceled;
  operation.promise = Promise.resolve().then(async () => {
    const busy = availability('save');
    if (busy) return busy;
    if (!isCurrent()) return { status: 'canceled' };
    const title = options.resolveTitle ? options.resolveTitle() : titleOverride;
    if (title === null && options.resolveTitle) return;
    state.sessionSavePending.value = true;
    recordSessionLifecycleEvent('session-save-pending-published');
    // Only mutable draft configuration/navigation is copied. Adopted biological
    // payloads, resources, catalogs, caches and Results retain their existing owner.
    // Session 46 keeps each drawing as its mode's slice (PD-OI-086); the
    // Features-list record is the shown Result's. `drawing` is the shown mode's.
    const drawing = state.drawings[state.mode.value];
    // Each drawing is checked as it is, before the JSON copy, so a value JSON
    // cannot write (Infinity, NaN) fails Save instead of being written as null.
    SLICE_MODES.forEach((mode) => validateCurrentWriterActiveConfig({
      mode, storedConfig: buildConfigData(state.drawings[mode])
    }));
    /** @param {'circular' | 'linear'} mode */
    const sliceOf = (mode) => buildModeSliceData(state.drawings[mode], mode, {
      selectedFeatureRecordIdx: mode === state.mode.value ? Number(state.selectedFeatureRecordIdx.value) || 0 : 0
    });
    const modes = { circular: sliceOf('circular'), linear: sliceOf('linear') };
    SLICE_MODES.forEach((mode) => validateModeSliceFields(modes[mode]));
    const { artifact, otherArtifact } = chooseSessionArtifacts(
      /** @type {{ readOtherModeArtifact?: () => Readonly<ArtifactSlot> | null }} */ (options).readOtherModeArtifact?.() ?? null
    );
    if (!artifact.committedCanonicalSession
      && hasBiologicalSessionInputs({ ...state.files, linearSeqs: state.linearSeqs })) {
      assertActiveModeInputs();
    }
    const artifactUi = artifactUiState(artifact);
    const savedUi = {
      mode: state.mode.value,
      zoom: state.zoom.value,
      canvasPan: { x: state.canvasPan.x, y: state.canvasPan.y },
      selectedResultIndex: artifactUi.selectedResultIndex,
      generatedLegendPosition: artifactUi.generatedLegendPosition,
      generatedMultiRecordCanvas: artifactUi.generatedMultiRecordCanvas,
      generatedCircularPlotTitlePosition: artifactUi.generatedCircularPlotTitlePosition,
      featurePanelTab: state.featurePanelTab.value,
      cInputType: state.cInputType.value,
      lInputType: state.lInputType.value,
      downloadDpi: state.downloadDpi.value,
      autoLabelReflow: Boolean(state.autoLabelReflowEnabled.value),
      paletteInstantPreviewEnabled: Boolean(state.paletteInstantPreviewEnabled.value),
      appliedPaletteName: artifactUi.appliedPaletteName,
      appliedPaletteColors: artifactUi.appliedPaletteColors,
      // App-level settings shared by both drawings.
      losatExecution: buildLosatExecutionData(),
      richFeaturePopup: Boolean(state.richFeaturePopup.value)
    };
    const prepared = await options.beforeExport?.({ draftRequest: !artifact.committedCanonicalSession });
    if (!isCurrent()) return { status: 'canceled' };
    return exportSessionDocument(title, {
      ...options, ...prepared, drawing, modes, savedUi, isCurrent, artifact, otherArtifact
    });
  }).catch((error) => {
    recordSessionLifecycleEvent('session-save-error');
    if (!isCurrent() || (state.errorLog.value !== previousAlert && state.errorLog.value !== null)) {
      return { status: 'stale' };
    }
    if (typeof options.onError !== 'function') throw error;
    const model = normalizeUserFacingError(error, { operation: 'session-save' });
    options.onError(model);
    return { status: 'error', error: model };
  }).finally(() => {
    if (sessionSaveInFlight === operation) {
      state.sessionSavePending.value = false;
      sessionSaveInFlight = null;
    }
    recordSessionLifecycleEvent('session-save-pending-cleared');
  });
  sessionSaveInFlight = operation;
  return operation.promise;
};

export const importSession = async (event, options = {}) => {
  const input = event?.target;
  const file = input?.files?.[0];
  if (!file) return { status: 'skipped' };
  const { availability = sessionOperationAvailability } = options;
  const busy = availability('load');
  if (busy) {
    input.value = '';
    return busy;
  }
  const operation = { canceled: false, controller: new AbortController() };
  activeSessionImport = operation;
  const isCurrent = () => activeSessionImport === operation && !operation.canceled;
  state.sessionImportPending.value = true;
  recordSessionLifecycleEvent('session-import-pending-published');
  try {
    await options.beforeImport?.();
    if (!isCurrent()) return { status: 'canceled' };
    const result = await importSessionDocument({ target: { files: [file], value: input.value } }, {
      ...options, isCurrent, signal: operation.controller.signal
    });
    return result;
  } finally {
    if (activeSessionImport === operation) {
      activeSessionImport = null;
      state.sessionImportPending.value = false;
    }
    input.value = '';
    recordSessionLifecycleEvent('session-import-pending-cleared');
    // Discovery watchers see the cleared pending flag while an adopted
    // current-schema Session still defers resource reads; clear it afterwards.
    await nextTick();
    if (!activeSessionImport) state.sessionResourceDiscoveryDeferred.value = false;
  }
};
