import { diagnosticError, normalizeUserFacingError } from './error-normalization.js';
import { state, sessionOperationAvailability, normalizeLinearSeqList, collapseEmptyLinearSeqList } from '../state.js';
import { resolveColorToHex } from '../app/color-utils.js';
import {
  captureRightDrawerState,
  resetRightDrawerState,
  restoreRightDrawerState
} from '../app/right-drawer.js';
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
} from '../app/depth-track-state.js';
import {
  decodeDepthText,
  isEncodedDepthFileEntry
} from './depth-file-codec.js';
import {
  normalizeCollinearAnchorMode,
  normalizeCollinearSearchScope,
  normalizeOrthogroupMembershipMode
} from '../app/losat-normalization.js';
import { normalizeDefinitionLineStyleState } from '../app/definition-line-style-state.js';
import {
  migrateLegacyLinearLabelVisibility,
  requireLinearLabelVisibilityMode
} from '../app/linear-label-visibility.js';
import { isCliInvocationSessionExportable } from '../app/run-info.js';
import { migrateLegacyOrthogroupMembers } from './legacy-similarity-alignment.js';
import { normalizeCircularPlotTitlePosition } from '../app/plot-title-position.js';
import {
  migrateLegacyLayoutPreferences,
  normalizeLayoutPreferences,
  replaceLayoutPreferences,
  resolveActiveLayoutPreference
} from '../app/layout-preferences.js';
import { reconcileImportedLinearTypographyLink } from '../app/linear-typography.js';
import {
  serializeFeatureVisibilityRules,
  normalizeFeatureVisibilityRule,
  normalizeVisibilityMode,
  splitLegacyVisibilityRules
} from '../app/feature-visibility.js';
import {
  buildSessionFeatureRecoveryPlan,
  classifyFeatureMetadataState,
  hasUsableBiologicalFeatureCatalog
} from '../app/session-feature-metadata.js';
import {
  analyzeCatalogSequenceSourceCoverage,
  buildRestoredMatchSequenceSources,
  resolveCircularComparisonSequenceAvailability
} from '../app/match-sequences.js';
import {
  CANONICAL_REQUEST_SCHEMA,
  buildCanonicalRenderRequest,
  managedConfigOverridePathsForMode,
  promoteCanonicalRenderRequestToCurrent,
  projectCanonicalSessionRequest,
  projectSettingsOnlySession
} from './session-request.js';
import {
  createDefaultLinearComparisonPlan,
  createLinearComparisonEdge,
  linearComparisonEdgeKey,
  normalizeLinearComparisonPlan,
  reconcileLinearComparisonPlan,
  resolveLinearComparisonPlan
} from '../app/linear-comparisons.js';
import { buildSessionResources as assembleSessionResources } from './session-resources.js';
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
import { normalizeAnnotationSets } from '../app/annotations/state.js';
import { applySpecificRuleProvenance } from '../app/specific-color-rules.js';
import { applyStrokeOverridesToSvg } from '../app/legend/stroke-actions.js';
import { normalizeLegacyLegendEntryGroups } from './svg-result-normalization.js';
import {
  COMPOSITION_METADATA_ATTRIBUTE,
  COMPOSITION_SCHEMA_ATTRIBUTE,
  normalizeLegacyComposition
} from '../app/legend-layout/composition-actions.js';
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
  projectArtifactState,
  projectDocumentMetadata,
  projectWebOnlyEditorMetadata,
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
} from '../app/current-option-values.js';
import {
  validateSimilarityAlignmentResetReceipt,
  CIRCULAR_TRACK_SLOT_SCHEMA_VERSION,
  CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS,
  createDefaultLosatpHitLimits,
  LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION,
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

const { nextTick } = window.Vue;

export const SESSION_VERSION = 44;
const CURRENT_AUTHORITY_SESSION_MIN_VERSION = 40;
const LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION = 32;
const SUPPORTED_SESSION_VERSIONS = new Set([
  27, 28, 29, 30, 31, 32, 33, 39, 40, 41, 42, SESSION_VERSION
]);
const CURRENT_ARTIFACT_SESSION_MIN_VERSION = 39;
const LOSAT_DERIVED_CACHE_LIMIT = 16;
// D-25 (PD-OI-079): a Result without current feature metadata (a legacy
// Session) is saved only after one Generate; the error offers that Generate.
const sessionSaveRequiresGenerate = () => diagnosticError(
  'SESSION_SAVE_REQUIRES_GENERATE', {}, { operation: 'session-save', stage: 'result-admission' }
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
  return hasColorEntries(colors) ? state.normalizePaletteColors(cloneColors(colors)) : null;
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

const normalizeFeatureVisibilityOverridesForSession = (overrides) => {
  const normalized = {};
  if (!overrides || typeof overrides !== 'object' || Array.isArray(overrides)) return normalized;
  Object.entries(overrides).forEach(([featureIdRaw, modeRaw]) => {
    const featureId = String(featureIdRaw || '').trim();
    const mode = normalizeVisibilityMode(modeRaw);
    if (!featureId || mode === 'default') return;
    normalized[featureId] = mode;
  });
  return normalized;
};

const splitFeatureVisibilityStateForSession = (features = {}) => {
  if (Array.isArray(features.featureVisibilityManualRules)) {
    return {
      manualRules: normalizeFeatureVisibilityRulesForSession(features.featureVisibilityManualRules),
      overrides: normalizeFeatureVisibilityOverridesForSession(features.featureVisibilityOverrides)
    };
  }
  if (Array.isArray(features.featureVisibilityRules)) {
    return splitLegacyVisibilityRules(features.featureVisibilityRules);
  }
  return {
    manualRules: [],
    overrides: normalizeFeatureVisibilityOverridesForSession(features.featureVisibilityOverrides)
  };
};

const replaceFeatureVisibilityState = (features = {}) => {
  const { manualRules, overrides } = splitFeatureVisibilityStateForSession(features);
  state.featureVisibilityManualRules.splice(
    0,
    state.featureVisibilityManualRules.length,
    ...normalizeFeatureVisibilityRulesForSession(manualRules)
  );
  replacePlainObject(state.featureVisibilityOverrides, overrides);
};

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

const replaceQualifierPriorityRules = (rules) => {
  state.manualPriorityRules.splice(
    0,
    state.manualPriorityRules.length,
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

const withHistoricalPairwiseMatchStyleFallback = (configData, mode = null) => {
  if (!isPlainObject(configData) || !isPlainObject(configData.adv)) return configData;
  const adv = Object.prototype.hasOwnProperty.call(configData.adv, 'pairwise_match_style')
    ? configData.adv
    : { ...configData.adv, pairwise_match_style: 'ribbon' };
  const profiles = configData.modeProfiles;
  const activeMode = ['circular', 'linear'].includes(mode)
    ? mode
    : profiles?.activeMode;
  if (
    !isPlainObject(profiles)
    || !isPlainObject(profiles.profiles)
    || !['circular', 'linear'].includes(activeMode)
    || !isPlainObject(profiles.profiles[activeMode])
  ) {
    return adv === configData.adv ? configData : { ...configData, adv };
  }
  const activeProfile = profiles.profiles[activeMode];
  const values = isPlainObject(activeProfile.values) ? activeProfile.values : {};
  if (Object.prototype.hasOwnProperty.call(values, 'pairwise_match_style')) {
    return adv === configData.adv ? configData : { ...configData, adv };
  }
  const managed = isPlainObject(activeProfile.managed) ? activeProfile.managed : {};
  return {
    ...configData,
    adv,
    modeProfiles: {
      ...profiles,
      profiles: {
        ...profiles.profiles,
        [activeMode]: {
          ...activeProfile,
          values: { ...values, pairwise_match_style: adv.pairwise_match_style },
          managed: { ...managed, pairwise_match_style: false }
        }
      }
    }
  };
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
  return withCurrentLinearLabelVisibility(withHistoricalPairwiseMatchStyleFallback({
    ...migratedNames,
    ...(form === undefined ? {} : { form }),
    ...(adv === undefined ? {} : { adv })
  }));
};

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
    sourceSessionVersion <= LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION
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

let lastSessionFilename = null;
let preservedCliOptions = null;
let committedCanonicalSession = null;
let activeSessionResourceTable = null;
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

export const buildConfigData = () => ({
  form: state.form,
  adv: {
    ...state.adv,
    feature_shapes: {
      repeat_region: defaultFeatureRendering('repeat_region'),
      ...normalizeFeatureRenderingMap(state.adv.feature_shapes || {})
    }
  },
  losat: cloneJsonData(state.losat || {}),
  cliOptions: preservedCliOptions ? cloneJsonData(preservedCliOptions) : undefined,
  colors: state.currentColors.value,
  palette: state.selectedPalette.value,
  paletteInstantPreviewEnabled: Boolean(state.paletteInstantPreviewEnabled.value),
  rules: state.manualSpecificRules,
  qualifierPriorityRules: cloneQualifierPriorityRules(state.manualPriorityRules),
  filterMode: state.filterMode.value,
  whitelist: state.manualWhitelist,
  blacklistText: state.manualBlacklist.value,
  losatProgram: state.losatProgram.value,
  circularConservation: state.circularConservation,
  annotationSets: normalizeAnnotationSets(state.annotationSets),
  recordDisplayDrafts: cloneJsonData(state.recordDisplayDrafts),
  featurePlacementOverrides: cloneJsonData(state.featurePlacementOverrides),
  modeProfiles: state.modeProfileStateManager?.exportState?.(),
  linearRecordLayout: {
    enabled: Boolean(state.linearRecordLayoutEnabled.value),
    recordGap: Number(state.linearRecordGap.value) || 0,
    rows: (state.linearRecordRows || []).map((entry) => ({
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
  linearComparisonPlan: serializeLinearComparisonPlan(state.linearComparisonPlan),
  importedComparisonResolution: serializeImportedComparisonResolution(
    state.importedComparisonIntent
  ),
  unmanagedConfigOverrides: cloneJsonData(state.unmanagedConfigOverrides || {}),
  webEdits: {
    orthogroupNameOverrides: cloneStringMap(state.orthogroupNameOverrides),
    orthogroupDescriptionOverrides: cloneStringMap(state.orthogroupDescriptionOverrides),
    orthogroupDormantOverrides: normalizeOrthogroupDormantOverrides(state.orthogroupDormantOverrides)
  }
});

const defaultEditorStateData = () => ({
  legend: {
    entries: [],
    deletedEntries: [],
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

export const buildEditorStateData = () => ({
  legend: {
    entries: cloneJsonArray(state.legendEntries.value),
    deletedEntries: cloneJsonArray(state.deletedLegendEntries.value),
    originalOrder: cloneJsonArray(state.originalLegendOrder.value),
    originalColors: cloneStringMap(state.originalLegendColors.value),
    colorOverrides: cloneJsonObject(state.legendColorOverrides),
    strokeOverrides: cloneJsonObject(state.legendStrokeOverrides),
    addedCaptions: Array.from(state.addedLegendCaptions.value || [])
      .map((caption) => String(caption || '').trim())
      .filter(Boolean)
  },
  featureStrokes: {
    overrides: cloneJsonObject(state.featureStrokeOverrides)
  },
  originalSvgStroke: {
    color: state.originalSvgStroke.value?.color ?? null,
    width: state.originalSvgStroke.value?.width ?? null
  },
  alignmentResetReceipt: cloneJsonValue(state.similarityAlignmentResetReceipt?.value, null),
  featureCatalog: admittedFeatureCatalog(state.featureCatalog?.value)
});

const normalizeEditorStateData = (editorState = {}, { featureCatalog = undefined } = {}) => {
  const defaults = defaultEditorStateData();
  const source = isPlainObject(editorState) ? editorState : {};
  const legend = isPlainObject(source.legend) ? source.legend : {};
  const featureStrokes = isPlainObject(source.featureStrokes) ? source.featureStrokes : {};
  const originalSvgStroke = isPlainObject(source.originalSvgStroke) ? source.originalSvgStroke : {};

  return {
    legend: {
      entries: normalizeSessionLegendEntries(legend.entries),
      deletedEntries: normalizeSessionLegendEntries(legend.deletedEntries),
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

export const applyEditorStateData = (
  editorState = {},
  { normalized: alreadyNormalized = false } = {}
) => {
  const normalized = alreadyNormalized
    ? editorState
    : normalizeEditorStateData(editorState);

  if (state.similarityAlignmentResetReceipt) {
    state.similarityAlignmentResetReceipt.value = normalized.alignmentResetReceipt ?? null;
  }
  state.legendEntries.value = normalized.legend.entries;
  state.deletedLegendEntries.value = normalized.legend.deletedEntries;
  state.originalLegendOrder.value = normalized.legend.originalOrder;
  state.originalLegendColors.value = normalized.legend.originalColors;
  replacePlainObject(state.legendColorOverrides, normalized.legend.colorOverrides);
  replacePlainObject(state.legendStrokeOverrides, normalized.legend.strokeOverrides);
  state.addedLegendCaptions.value = new Set(normalized.legend.addedCaptions);
  replacePlainObject(state.featureStrokeOverrides, normalized.featureStrokes.overrides);
  state.originalSvgStroke.value = normalized.originalSvgStroke;
  if (state.featureCatalog) {
    state.featureCatalog.value = admittedFeatureCatalog(normalized.featureCatalog);
  }
};

const SESSION_FORMAT_ERROR = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'SESSION_FORMAT' });
const validateSessionVersion = version => {
  if (!Number.isInteger(version)) throw SESSION_FORMAT_ERROR();
  if (version > SESSION_VERSION || !SUPPORTED_SESSION_VERSIONS.has(version)) {
    throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
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
  if (!isPlainObject(container)) {
    throw new Error(`Session ${field} must be an object when present.`);
  }
  const entries = Object.prototype.hasOwnProperty.call(container, 'entries')
    ? container.entries
    : [];
  if (!Array.isArray(entries)) {
    throw new Error(`Session ${field}.entries must be an array.`);
  }
  return entries;
};

const rejectInvalidLosatCacheKeys = (entries, owner, { requireKey = false } = {}) => {
  const seen = new Set();
  for (const [index, entry] of entries.entries()) {
    const key = isPlainObject(entry) && typeof entry.key === 'string'
      ? entry.key
      : '';
    if (!key) {
      if (requireKey) {
        throw new Error(
          `LOSAT cache entry at losatCache.entries[${index}] requires a key.`
        );
      }
      continue;
    }
    if (seen.has(key)) {
      throw new Error(`Duplicate ${owner} cache key: ${JSON.stringify(key)}.`);
    }
    seen.add(key);
  }
};

function* sessionLosatArtifactSteps(data, sourceSessionVersion) {
  if (sourceSessionVersion < CURRENT_ARTIFACT_SESSION_MIN_VERSION) return;
  const rawEntries = sessionArtifactEntries(data, 'losatCache');
  const derivedEntries = sessionArtifactEntries(data, 'losatDerivedCache');
  const manifest = data.proteinIdentityManifest;
  rejectInvalidLosatCacheKeys(rawEntries, 'LOSAT', { requireKey: true });
  rejectInvalidLosatCacheKeys(derivedEntries, 'derived LOSATP');

  if (sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION) {
    recordStructuralMetric('currentSessionPreflightProteinManifestValidationCount');
  }
  const identityIndex = buildValidatedProteinIdentityIndex(manifest);
  if (!identityIndex) {
    throw new Error(
      `Session version ${sourceSessionVersion} requires a valid schema-2 protein manifest.`
    );
  }
  let invalidDerivedEntry = false;
  try {
    for (const entry of rawEntries) {
      const classification = classifyRawLosatCacheEntry(entry);
      if (!['protein-current', 'nucleotide-current'].includes(classification)) {
        throw new Error(
          `Session version ${sourceSessionVersion} contains a non-current raw LOSAT entry.`
        );
      }
      yield;
      if (classification !== 'protein-current') continue;
      if (sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION) {
        recordStructuralMetric('currentSessionPreflightProteinRawTextValidationCount');
      }
      if (
        !validateProteinRawEntryReferences(entry, manifest, { identityIndex })
      ) {
        throw new Error(
          `Session version ${sourceSessionVersion} contains an unresolved protein raw entry.`
        );
      }
    }
    invalidDerivedEntry = derivedEntries.some(
      (entry) => !validateDerivedProteinReferences(entry, manifest, { identityIndex })
    );
  } finally {
    releaseValidatedProteinIdentityIndex(identityIndex);
  }
  if (invalidDerivedEntry) {
    throw new Error(
      `Session version ${sourceSessionVersion} contains an invalid derived LOSATP entry.`
    );
  }
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

const migrateSessionDataToCurrent = (data, sourceSessionVersion) => {
  const readsLegacyOptionValues = sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION;
  const migratedOptions = readsLegacyOptionValues
    ? migratePersistedWebOptionValues(data.config)
    : data.config;
  return {
    ...data,
    version: SESSION_VERSION,
    config: migrateLegacyFeatureRenderingConfig(
      migrateImportedLinearTrackSlots(
        migrateImportedCircularTrackSlots(migratedOptions),
        sourceSessionVersion
      ),
      sourceSessionVersion <= 33
    )
  };
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

const applyLegacyConfigPayload = (data) => {
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
  state.suppressCircularMultiRecordDefaults.value = shouldSuppressCircularMultiRecordDefaults(migrated.form);
  applyConfigData(migrated);
  restorePaletteStateAfterConfigImport();
};

const shouldSuppressCircularMultiRecordDefaults = (incomingForm) => {
  if (state.mode.value !== 'circular') return false;
  if (!incomingForm || typeof incomingForm !== 'object' || Array.isArray(incomingForm)) return false;
  if (!Object.prototype.hasOwnProperty.call(incomingForm, 'multi_record_canvas')) return false;
  return state.form.multi_record_canvas === false && incomingForm.multi_record_canvas === true;
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
  validateCurrentWriterActiveConfig({ mode, storedConfig });
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
  let adoptedSession = null;
  let currentResourceTable = null;
  let validatedFeatureCatalog = null;
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
  const comparisonClassification = classifyImportedComparisonIntent({
    renderRequest: data.renderRequest,
    resources: data.resources
  });
  const projectionRenderRequest = (
    data.renderRequest?.mode === 'linear'
    && comparisonClassification.disposition
      !== IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE
  )
    ? { ...data.renderRequest, comparisons: [] }
    : data.renderRequest;
  let currentStoredConfig = sourceSessionVersion < SESSION_VERSION
    ? withCurrentLinearLabelVisibility(data.config)
    : data.config;
  if (sourceSessionVersion < SESSION_VERSION && isPlainObject(currentStoredConfig)
    && Object.prototype.hasOwnProperty.call(currentStoredConfig, 'recordDisplayDrafts')) {
    currentStoredConfig = {
      ...currentStoredConfig,
      recordDisplayDrafts: migrateLegacyRecordDisplayDrafts(
        currentStoredConfig.recordDisplayDrafts
      )
    };
  }
  const runtimeStoredConfig = currentSession && Object.prototype.hasOwnProperty.call(data, 'config')
    ? migrateImportedLinearTrackSlots(
        migrateImportedCircularTrackSlots(withHistoricalPairwiseMatchStyleFallback(
          currentStoredConfig,
          data.ui?.mode || data.renderRequest?.mode
        )),
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
        initializeCliInputs: !Object.hasOwn(data, 'config')
          && data.cliInvocation?.generatedBy === 'gbdraw',
        fileBindings: data.cliInvocation?.fileBindings,
        linearTrackSlotSchemaVersion: sourceSessionVersion <= LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION
          ? LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION
          : LINEAR_TRACK_SLOT_SCHEMA_VERSION,
        repairInvalidComparisonHeight: sourceSessionVersion >= 31 && sourceSessionVersion <= 33,
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
  if (sourceSessionVersion < SESSION_VERSION && restoredConfig) {
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
  if (canonicalProjection && sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION
    && Object.prototype.hasOwnProperty.call(data, 'config')) {
    recordSessionLifecycleEvent('current-draft-validation-start');
    restoredConfig = restoreCurrentWriterActiveConfig({
      mode: canonicalProjection.mode,
      projectedConfig: canonicalProjection.config,
      storedConfig: runtimeStoredConfig
    });
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
  const hasCurrentStoredUnmanagedOverrides = currentSession
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
    comparisonClassification,
    unmanagedConfigValidation
  };
};

const LEGACY_LAYOUT_PREFERENCE_FIELDS = Object.freeze([
  'legend',
  'circularLegendPosition',
  'linearLegendPosition',
  'circularPlotTitlePosition',
  'linearPlotTitlePosition',
  'circularSingleRecordLegendPosition',
  'circularSingleRecordPlotTitlePosition',
  'circularMultiRecordLegendPosition',
  'circularMultiRecordPlotTitlePosition'
]);

// Partial objects remain authoritative for session compatibility; normalization
// supplies the current defaults for omitted branches.
const hasStoredLayoutPreferences = (ui) => (
  isPlainObject(ui?.layoutPreferences) ||
  LEGACY_LAYOUT_PREFERENCE_FIELDS.some((field) => hasStoredLayoutValue(ui?.[field]))
);

// A saved layout owner (current or legacy ui fields) wins. Without one, a
// canonical Session takes the layout preferences projected from its committed
// request. Legacy fields migrate with the committed values (canonical) or the
// active values (other payloads) as their fallback.
const restoreLayoutPreferences = (ui = {}, { projected = null } = {}) => {
  if (projected && !hasStoredLayoutPreferences(ui)) {
    replaceLayoutPreferences(state.layoutPreferences, normalizeLayoutPreferences(projected));
    return;
  }
  const active = projected
    ? resolveActiveLayoutPreference(projected, state.mode.value, Boolean(state.form.multi_record_canvas))
    : { legend: state.form.legend, plotTitlePosition: state.adv.plot_title_position };
  const migrationUi = (
    !isPlainObject(ui.layoutPreferences) &&
    state.mode.value === 'linear' &&
    !hasStoredLayoutValue(ui.linearLegendPosition) &&
    hasStoredLayoutValue(ui.legend)
  )
    ? { ...ui, linearLegendPosition: ui.legend }
    : ui;
  replaceLayoutPreferences(
    state.layoutPreferences,
    migrateLegacyLayoutPreferences(migrationUi, {
      mode: state.mode.value,
      multiRecord: Boolean(state.form.multi_record_canvas),
      activeLegend: active.legend,
      activePlotTitlePosition: active.plotTitlePosition
    })
  );
};

// Captured or stored settings (History Undo and Redo, the failed Session Load
// rollback, a settings-only Session) pass `resolveTrackPlacements: false`: the
// track stacks are installed as given, so unset slot sides, lane directions,
// and axis indexes stay unset (R11).
export const applyConfigData = (data, { resolveTrackPlacements = true } = {}) => {
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
  if (data.form) safeDeepMerge(state.form, data.form);
  if (data.adv) {
    safeDeepMerge(state.adv, data.adv);
    ['scale_font_size', 'ruler_label_font_size'].forEach((field) => {
      if (Object.prototype.hasOwnProperty.call(data.adv, field)) {
        state.adv[field] = cloneJsonData(data.adv[field]);
      }
    });
  }
  replacePlainObject(
    state.unmanagedConfigOverrides,
    isPlainObject(data.unmanagedConfigOverrides)
      ? cloneJsonData(data.unmanagedConfigOverrides)
      : {}
  );
  state.recordDisplayDrafts.splice(0, state.recordDisplayDrafts.length, ...cloneJsonData(data.recordDisplayDrafts || []));
  replacePlainObject(state.featurePlacementOverrides, cloneJsonData(data.featurePlacementOverrides || {}));
  state.annotationSets.splice(
    0,
    state.annotationSets.length,
    ...normalizeAnnotationSets(data.annotationSets)
  );
  const linearLayout = data.linearRecordLayout && typeof data.linearRecordLayout === 'object'
    ? data.linearRecordLayout
    : null;
  // Omission takes the fresh default in every Session version.
  state.linearRecordLayoutEnabled.value = typeof linearLayout?.enabled === 'boolean'
    ? linearLayout.enabled
    : WEB_UX_PROFILE.linear.arrangeInRowsByDefault;
  const linearRecordGap = Number(linearLayout?.recordGap);
  state.linearRecordGap.value = Number.isFinite(linearRecordGap) && linearRecordGap >= 0
    ? linearRecordGap
    : 24;
  state.linearRecordRows.splice(
    0,
    state.linearRecordRows.length,
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
  if (state.linearRecordTranslations) {
    state.linearRecordTranslations.value = cloneJsonData(
      Array.isArray(linearLayout?.recordTranslations)
        ? linearLayout.recordTranslations
        : []
    );
  }
  if (state.similarityAlignmentPlan) {
    state.similarityAlignmentPlan.value = linearLayout?.similarityAlignment
      ? cloneJsonData(linearLayout.similarityAlignment)
      : null;
  }
  replaceLinearComparisonPlan(
    state.linearComparisonPlan,
    data.linearComparisonPlan || createDefaultLinearComparisonPlan()
  );
  state.importedComparisonIntent.action = Object.values(IMPORTED_COMPARISON_ACTIONS)
    .includes(data.importedComparisonResolution?.action)
    ? data.importedComparisonResolution.action
    : null;
  state.adv.rich_feature_popup = data?.adv?.rich_feature_popup !== false;
  state.adv.label_placement = requireCurrentLinearLabelPlacement(
    state.adv.label_placement
  );
  state.adv.label_rendering = normalizeLabelRendering(state.adv.label_rendering);
  state.adv.circular_label_placement =
    String(state.adv.circular_label_placement || '').trim().toLowerCase() === 'radial'
      ? 'radial'
      : 'horizontal';
  if (state.adv.label_placement === 'above_feature') {
    state.adv.label_rendering = 'auto';
  }
  state.adv.circular_label_spacing = normalizePositiveNumberOrNull(state.adv.circular_label_spacing);
  state.adv.linear_label_spacing = normalizePositiveNumberOrNull(state.adv.linear_label_spacing);
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
  state.form.linear_track_layout = requireCurrentLinearTrackLayout(
    state.form.linear_track_layout
  );
  state.form.plot_title = String(state.form.plot_title || '');
  // `form.legend` and `adv.plot_title_position` are accessors over
  // `layoutPreferences` (state.js) whose setters normalize a merged value.
  // Writing the resolved value back would only pin an unset Circular
  // multi-record preference, so neither is rewritten here.
  state.adv.feature_shapes = normalizeFeatureRenderingMap(state.adv.feature_shapes);
  Object.assign(state.adv, normalizedPersistedArrowGeometryState(data.adv));
  state.adv.multi_record_size_mode = requireCurrentCircularMultiRecordSizeMode(
    state.adv.multi_record_size_mode
  );
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
  state.adv.depth_color = resolveColorToHex(String(state.adv.depth_color || '#4A90E2'));
  state.adv.depth_normalize = state.adv.depth_normalize === true;
  state.adv.depth_show_axis = state.adv.depth_show_axis !== false;
  state.adv.depth_show_ticks = state.adv.depth_show_ticks !== false;
  state.adv.depth_share_axis = state.adv.depth_share_axis === true;
  state.adv.depth_height = normalizePositiveNumberOrNull(state.adv.depth_height);
  state.adv.depth_width_circular = normalizePositiveNumberOrNull(state.adv.depth_width_circular);
  state.adv.circular_track_slots_schema_version = CIRCULAR_TRACK_SLOT_SCHEMA_VERSION;
  state.adv.circular_track_slots_enabled = state.adv.circular_track_slots_enabled === true;
  if (resolveTrackPlacements) {
    {
      const normalizedSlots = normalizeCircularTrackSlots(
        state.adv.circular_track_slots,
        state.adv.nt,
        state.form.track_type
      );
      const importedAxis = clampCircularTrackAxisIndex(
        state.adv.circular_track_slots_axis_index,
        normalizedSlots.length
      );
      state.adv.circular_track_slots_axis_index = importedAxis === null
        ? inferLegacyAxisIndexFromFeature(normalizedSlots, state.form.track_type)
        : importedAxis;
    }
    state.adv.circular_track_slots.splice(
      0,
      state.adv.circular_track_slots.length,
      ...applyCircularTrackOrderPlacements(
        state.adv.circular_track_slots,
        state.adv.nt,
        state.form.track_type,
        state.adv.circular_track_slots_axis_index
      )
    );
    state.adv.linear_track_slots_schema_version = LINEAR_TRACK_SLOT_SCHEMA_VERSION;
    state.adv.linear_track_slots_enabled = state.adv.linear_track_slots_enabled === true;
    {
      const normalizedLinearSlots = normalizeLinearTrackSlots(
        state.adv.linear_track_slots,
        state.adv.nt,
        state.form.linear_track_layout
      );
      state.adv.linear_track_slots_axis_index = clampLinearTrackAxisIndex(
        state.adv.linear_track_slots_axis_index,
        normalizedLinearSlots.length
      );
      state.adv.linear_track_slots_axis_index = resolveLinearTrackAxisIndex(
        normalizedLinearSlots,
        state.adv.linear_track_slots_axis_index
      );
      state.adv.linear_track_slots.splice(
        0,
        state.adv.linear_track_slots.length,
        ...applyLinearTrackOrderPlacements(
          normalizedLinearSlots,
          state.adv.linear_track_slots_axis_index,
          state.adv.nt,
          state.form.linear_track_layout
        )
      );
    }
  } else {
    // The merge installs the caller's slot objects; copy them so a later edit
    // never writes into a History entry or another caller-owned snapshot.
    ['circular_track_slots', 'linear_track_slots'].forEach((key) => {
      state.adv[key].splice(0, state.adv[key].length, ...cloneJsonData(state.adv[key]));
    });
  }
  state.adv.depth_window_size = normalizePositiveNumberOrNull(state.adv.depth_window_size);
  state.adv.depth_step_size = normalizePositiveNumberOrNull(state.adv.depth_step_size);
  state.adv.depth_large_tick_interval = normalizePositiveNumberOrNull(
    state.adv.depth_large_tick_interval
  );
  state.adv.depth_small_tick_interval = normalizePositiveNumberOrNull(state.adv.depth_small_tick_interval);
  state.adv.depth_tick_font_size = normalizePositiveNumberOrNull(state.adv.depth_tick_font_size);
  state.adv.depth_tracks.splice(
    0,
    state.adv.depth_tracks.length,
    ...normalizeDepthTracks(state.adv.depth_tracks, state.adv)
  );
  state.adv.gc_content_mode = String(state.adv.gc_content_mode || '').trim().toLowerCase() === 'percent'
    ? 'percent'
    : 'deviation';
  state.adv.gc_content_show_axis = state.adv.gc_content_show_axis !== false;
  state.adv.gc_content_show_ticks = state.adv.gc_content_show_ticks !== false;
  state.adv.gc_content_tick_interval = normalizePositiveNumberOrNull(state.adv.gc_content_tick_interval);
  state.adv.gc_content_small_tick_interval = normalizePositiveNumberOrNull(state.adv.gc_content_small_tick_interval);
  state.adv.gc_content_tick_font_size = normalizePositiveNumberOrNull(state.adv.gc_content_tick_font_size);
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
  state.adv.center_reserved_radius = normalizeNonNegativeNumberOrNull(state.adv.center_reserved_radius);
  state.adv.depth_min = normalizeNonNegativeNumberOrNull(state.adv.depth_min);
  state.adv.depth_max = normalizeNonNegativeNumberOrNull(state.adv.depth_max);
  if (
    state.adv.depth_min !== null &&
    state.adv.depth_max !== null &&
    state.adv.depth_min > state.adv.depth_max
  ) {
    state.adv.depth_max = null;
  }
  const normalizeFiniteNumberOrFallback = (value, fallback) => {
    if (value === null || value === undefined || value === '') return fallback;
    const numeric = Number(value);
    return Number.isFinite(numeric) ? numeric : fallback;
  };
  state.adv.gc_content_min_percent = normalizeFiniteNumberOrFallback(state.adv.gc_content_min_percent, 0);
  state.adv.gc_content_max_percent = normalizeFiniteNumberOrFallback(state.adv.gc_content_max_percent, 100);
  if (state.adv.gc_content_min_percent > state.adv.gc_content_max_percent) {
    state.adv.gc_content_max_percent = state.adv.gc_content_min_percent;
  }
  state.adv.linear_show_replicon = state.adv.linear_show_replicon === true;
  state.adv.linear_accession_visibility = requireLinearLabelVisibilityMode(
    state.adv.linear_accession_visibility,
    'Linear Accession visibility'
  );
  state.adv.linear_length_visibility = requireLinearLabelVisibilityMode(
    state.adv.linear_length_visibility,
    'Linear Length / Coordinates visibility'
  );
  state.adv.linear_definition_line_styles = normalizeDefinitionLineStyleState(
    state.adv.linear_definition_line_styles
  );
  state.adv.pairwise_match_style = normalizeCurrentPairwiseMatchStyle(
    state.adv.pairwise_match_style,
    'ribbon'
  );
  if (data.losat) {
    safeDeepMerge(state.losat, data.losat);
    const rawParallelWorkers = String(data.losat.parallelWorkers ?? '').trim().toLowerCase();
    const parsedParallelWorkers = Number(rawParallelWorkers);
    state.losat.parallelWorkers = Number.isInteger(parsedParallelWorkers) && parsedParallelWorkers >= 1
      ? rawParallelWorkers
      : undefined;
    const rawExecutionMode = String(data.losat.executionMode ?? '').trim().toLowerCase();
    state.losat.executionMode = ['auto', 'serial', 'threaded'].includes(rawExecutionMode)
      ? rawExecutionMode
      : 'auto';
    const rawThreadsPerJob = String(data.losat.threadsPerJob ?? 'auto').trim().toLowerCase();
    const parsedThreadsPerJob = Number(rawThreadsPerJob);
    state.losat.threadsPerJob = rawThreadsPerJob === 'auto' ||
      (Number.isInteger(parsedThreadsPerJob) && parsedThreadsPerJob >= 1)
      ? rawThreadsPerJob
      : 'auto';
    const rawTotalThreadBudget = String(data.losat.totalThreadBudget ?? 'safe').trim().toLowerCase();
    const parsedTotalThreadBudget = Number(rawTotalThreadBudget);
    state.losat.totalThreadBudget = ['safe', 'auto', 'available'].includes(rawTotalThreadBudget) ||
      (Number.isInteger(parsedTotalThreadBudget) && parsedTotalThreadBudget >= 1)
      ? (rawTotalThreadBudget === 'auto' ? 'safe' : rawTotalThreadBudget)
      : 'safe';
    state.losat.blastp.mode = normalizeBlastpMode(state.losat.blastp?.mode);
    state.losat.blastp.collinearInferOrthogroups = data.losat.blastp?.collinearInferOrthogroups ?? (state.losat.blastp.mode === 'collinear');
    state.losat.blastp.hitLimitsByMode = {
      ...createDefaultLosatpHitLimits(), ...cloneJsonData(data.losat.blastp?.hitLimitsByMode || {})
    };
    state.losat.blastp.maxHits = normalizePositiveInteger(state.losat.blastp?.maxHits, 5);
    state.losat.blastp.candidateLimit = normalizePositiveInteger(
      state.losat.blastp?.candidateLimit,
      null
    );
    if (
      state.losat.blastp.orthogroupMemberMaxHits === undefined &&
      state.losat.blastp.orthogroupMaxHits !== null &&
      state.losat.blastp.orthogroupMaxHits !== undefined
    ) {
      state.losat.blastp.orthogroupMemberMaxHits = state.losat.blastp.orthogroupMaxHits;
    }
    state.losat.blastp.orthogroupMembershipMode = normalizeOrthogroupMembershipMode(state.losat.blastp?.orthogroupMembershipMode);
    state.losat.blastp.orthogroupMemberMaxHits = normalizePositiveInteger(state.losat.blastp?.orthogroupMemberMaxHits, null);
    state.losat.blastp.collinearMinAnchors = normalizePositiveInteger(state.losat.blastp?.collinearMinAnchors, 1);
    {
      const maxGap = Number(state.losat.blastp?.collinearMaxUnitGap);
      state.losat.blastp.collinearMaxUnitGap = Number.isInteger(maxGap) && maxGap >= 0 ? maxGap : 0;
      const diagonalDrift = Number(state.losat.blastp?.collinearMaxDiagonalDrift);
      state.losat.blastp.collinearMaxDiagonalDrift = Number.isInteger(diagonalDrift) && diagonalDrift >= 0 ? diagonalDrift : 0;
      const mergeConflicts = Number(state.losat.blastp?.collinearMaxConflictsInMergeGap);
      state.losat.blastp.collinearMaxConflictsInMergeGap = Number.isInteger(mergeConflicts) && mergeConflicts >= 0 ? mergeConflicts : 1;
      const paralogLinks = Number(state.losat.blastp?.collinearMaxParalogLinksPerOrthogroup);
      state.losat.blastp.collinearMaxParalogLinksPerOrthogroup = Number.isInteger(paralogLinks) && paralogLinks > 0 ? paralogLinks : 2;
      state.losat.blastp.collinearColorMode = normalizeCollinearColorMode(state.losat.blastp?.collinearColorMode);
      const unitMode = String(state.losat.blastp?.collinearUnitMode || '').trim().toLowerCase();
      state.losat.blastp.collinearUnitMode = ['auto', 'cds', 'locus'].includes(unitMode) ? unitMode : 'auto';
      state.losat.blastp.collinearAnchorMode = normalizeCollinearAnchorMode(state.losat.blastp?.collinearAnchorMode);
      const mergeOrientation = String(
        state.losat.blastp?.collinearMergeOrientation || ''
      ).trim().toLowerCase();
      state.losat.blastp.collinearMergeOrientation = [
        'strand',
        'order',
        'either'
      ].includes(mergeOrientation) ? mergeOrientation : 'either';
      state.losat.blastp.collinearSearchScope = normalizeCollinearSearchScope(state.losat.blastp?.collinearSearchScope);
    }
    delete state.losat.blastp.collinearBlockMergeGap;
    delete state.losat.blastp.collinearSingletonMergeGap;
    delete state.losat.blastp.orthogroupHitPolicy;
    delete state.losat.blastp.orthogroupMaxHits;
  }
  if (typeof data.paletteInstantPreviewEnabled === 'boolean') {
    state.paletteInstantPreviewEnabled.value = data.paletteInstantPreviewEnabled;
  }
  const importedPalette = String(data.palette || '').trim();
  if (importedPalette) state.selectedPalette.value = importedPalette;
  if (hasColorEntries(data.colors)) {
    if (data.colorsAreOverrides) {
      const paletteColors = paletteColorsFromDefinitions(state.selectedPalette.value) || {};
      state.currentColors.value = state.normalizePaletteColors({
        ...paletteColors,
        ...normalizeColorMap(data.colors)
      });
    } else {
      state.currentColors.value = state.normalizePaletteColors(normalizeColorMap(data.colors));
    }
  } else {
    const paletteColors = paletteColorsFromDefinitions(state.selectedPalette.value);
    if (paletteColors) state.currentColors.value = paletteColors;
  }

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
    state.fileLegendCaptions.value = new Set(
      state.manualSpecificRules
        .filter((rule) => rule.fromFile && rule.cap)
        .map((rule) => rule.cap)
    );
  }
  if (Object.prototype.hasOwnProperty.call(data, 'qualifierPriorityRules')) {
    replaceQualifierPriorityRules(data.qualifierPriorityRules);
  } else if (Object.prototype.hasOwnProperty.call(data, 'priorityRules')) {
    replaceQualifierPriorityRules(data.priorityRules);
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
  if (data.losatProgram) {
    const program = String(data.losatProgram);
    state.losatProgram.value = ['blastn', 'tblastx', 'blastp'].includes(program) ? program : 'blastn';
  }
  if (data.circularConservation) {
    safeDeepMerge(state.circularConservation, data.circularConservation);
  }
  state.circularConservation.enabled = state.circularConservation.enabled === true;
  state.circularConservation.source = normalizeCircularConservationSource(state.circularConservation.source);
  state.circularConservation.losat_program = normalizeCircularConservationLosatProgram(
    state.circularConservation.losat_program
  );
  state.circularConservation.subject_gencode = normalizePositiveInteger(state.circularConservation.subject_gencode, 1);
  state.circularConservation.reference = normalizeCircularConservationReference(state.circularConservation.reference);
  state.circularConservation.labels = String(state.circularConservation.labels || '');
  state.circularConservation.series.splice(
    0,
    state.circularConservation.series.length,
    ...normalizeCircularConservationSeries(state.circularConservation.series)
  );
  state.circularConservation.ring_width = normalizePositiveNumberOrNull(state.circularConservation.ring_width);
  state.circularConservation.ring_gap = normalizePositiveNumberOrNull(state.circularConservation.ring_gap);
  preservedCliOptions = isPlainObject(data.cliOptions) ? cloneJsonData(data.cliOptions) : null;
  const webEdits = data.webEdits && typeof data.webEdits === 'object' ? data.webEdits : {};
  if (Object.prototype.hasOwnProperty.call(webEdits, 'orthogroupNameOverrides')) {
    replaceStringMap(state.orthogroupNameOverrides, webEdits.orthogroupNameOverrides);
  }
  if (Object.prototype.hasOwnProperty.call(webEdits, 'orthogroupDescriptionOverrides')) {
    replaceStringMap(state.orthogroupDescriptionOverrides, webEdits.orthogroupDescriptionOverrides);
  }
  // Absent in older Sessions: no dormant names (D-21).
  clearObject(state.orthogroupDormantOverrides);
  Object.assign(state.orthogroupDormantOverrides, normalizeOrthogroupDormantOverrides(webEdits.orthogroupDormantOverrides));
  state.modeProfileStateManager?.importState?.(
    data.modeProfiles ?? null,
    state.mode.value,
    state.adv
  );
};

const restorePaletteStateAfterConfigImport = () => {
  const draftPaletteName = String(state.selectedPalette.value || state.appliedPaletteName.value || 'default');
  const draftColors = state.normalizePaletteColors(cloneColors(state.currentColors.value));
  const hasPreviewResults = Array.isArray(state.results.value) && state.results.value.length > 0;

  if (
    !hasPreviewResults ||
    state.paletteInstantPreviewEnabled.value ||
    draftPaletteName === String(state.appliedPaletteName.value || '')
  ) {
    state.appliedPaletteName.value = draftPaletteName;
    state.appliedPaletteColors.value = draftColors;
    state.pendingPaletteName.value = '';
    state.pendingPaletteColors.value = {};
    return;
  }

  state.pendingPaletteName.value = draftPaletteName;
  state.pendingPaletteColors.value = draftColors;
};

const restorePaletteStateFromSession = (ui = {}) => {
  const draftPaletteName = String(state.selectedPalette.value || state.appliedPaletteName.value || 'default');
  const draftColors = state.normalizePaletteColors(cloneColors(state.currentColors.value));
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

  state.appliedPaletteName.value = savedAppliedPaletteName;
  state.appliedPaletteColors.value = state.normalizePaletteColors(cloneColors(savedAppliedPaletteColors));

  if (!state.paletteInstantPreviewEnabled.value && savedPendingPaletteName) {
    state.pendingPaletteName.value = savedPendingPaletteName;
    state.pendingPaletteColors.value = state.normalizePaletteColors(cloneColors(savedPendingPaletteColors));
  } else {
    state.pendingPaletteName.value = '';
    state.pendingPaletteColors.value = {};
  }
};

const serializedFileDescriptors = new WeakMap();
const cacheSerializedFileDescriptor = (file, descriptor) => {
  setResourcePayloadOwner(descriptor, file);
  serializedFileDescriptors.set(file, descriptor);
  return descriptor;
};

const serializeFile = async (file) => {
  if (!file) return null;
  const source = getSessionResourceSource(file);
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

let activePreviewRuntime = null;

export const setPreviewRuntime = (runtime) => {
  activePreviewRuntime = runtime || null;
};

export const serializeResults = () => {
  if (activePreviewRuntime?.flushActiveResult) {
    activePreviewRuntime.flushActiveResult();
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

const restoredLosatCacheInfoIdentity = (entry) => {
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
  const resolved = state.linearComparisonResolution.value.edges.find(
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

  info.forEach((entry, idx) => {
    if (!entry || !entry.key) return;
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

const applyLosatCache = (
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
        ...restoredLosatCacheInfoIdentity(entry)
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

export const applyOrthogroupStateData = (
  orthogroupState = {}, { legacyRecords = null, catalogFeatureState = null } = {}
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

  replaceStringMap(state.orthogroupNameOverrides, orthogroupState.orthogroupNameOverrides);
  replaceStringMap(state.orthogroupDescriptionOverrides, orthogroupState.orthogroupDescriptionOverrides);
  Object.keys(state.orthogroupNameOverrides).forEach((id) => {
    if (!groupIdSet.has(id)) delete state.orthogroupNameOverrides[id];
  });
  Object.keys(state.orthogroupDescriptionOverrides).forEach((id) => {
    if (!groupIdSet.has(id)) delete state.orthogroupDescriptionOverrides[id];
  });
  if (Object.hasOwn(orthogroupState, 'orthogroupDormantOverrides')) {
    clearObject(state.orthogroupDormantOverrides);
    Object.assign(state.orthogroupDormantOverrides,
      normalizeOrthogroupDormantOverrides(orthogroupState.orthogroupDormantOverrides));
  }
};

const customDepthRequested = (mode, sourceState) => {
  const adv = sourceState?.adv || {};
  const form = sourceState?.form || {};
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

export const serializeActiveRenderFiles = async (
  mode = state.mode.value,
  sourceState = state,
  comparisonPlanOrOptions = null
) => {
  if (!['circular', 'linear'].includes(mode)) {
    throw new Error(`Unsupported render mode: ${String(mode)}.`);
  }
  const sourceFiles = sourceState.files || {};
  const normalizedLinearSeqs = mode === 'linear'
    ? normalizeLinearSeqList(sourceState.linearSeqs)
    : [];
  const depthRequested = customDepthRequested(mode, sourceState);
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
    { layoutEnabled: Boolean(sourceState.linearRecordLayoutEnabled?.value) }
  );
  const resolvedComparisonPlan = mode === 'linear'
    ? suppliedComparisonPlan || resolveLinearComparisonPlan({
        plan: sourceState.linearComparisonPlan,
        sequences: normalizedLinearSeqs,
        layout: sourceState.linearRecordLayoutEnabled?.value
          ? sourceState.linearRecordRows
          : [],
        losatProgram: sourceState.losatProgram?.value,
        blastpMode: sourceState.losat?.blastp?.mode
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
    && Boolean(sourceState.circularConservation?.enabled);
  const conservationSource = String(sourceState.circularConservation?.source || 'upload');
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

export const buildSessionResources = (sourceState, committedRequest) => (
  assembleSessionResources(sourceState, committedRequest)
);

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
  let nextLinearComparisons = null;
  let nextCircularState = null;
  if (projection.mode === 'linear') {
    nextLinearComparisons = deserializeCanonicalComparisons(
      projection.files.linearCanonicalComparisons,
      { adoptCanonicalPayloads: preserveAdoptedResources }
    );
  } else if (projection.files.c_conservation_blasts_source === 'losat-cache') {
    const projectedConservation = projection.config.circularConservation;
    const currentSeries = Array.isArray(state.circularConservation.series)
      ? state.circularConservation.series.map((entry) => cloneJsonData(entry))
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
      state.circularConservation.enabled = true;
      state.circularConservation.source = 'losat';
      state.circularConservation.reference = nextCircularState.projectedConservation.reference;
      state.circularConservation.labels = nextCircularState.series
        .map((entry) => entry.label).join(',');
      state.circularConservation.ring_width = nextCircularState.projectedConservation.ring_width;
      state.circularConservation.ring_gap = nextCircularState.projectedConservation.ring_gap;
      state.circularConservation.series.splice(
        0,
        state.circularConservation.series.length,
        ...nextCircularState.series
      );
    }
  }
};

export const getCommittedCanonicalRenderRequest = () => (
  committedCanonicalSession?.renderRequest || null
);

export const getCommittedCanonicalSession = () => committedCanonicalSession;

export const canonicalRenderArtifactOwner = Object.freeze({
  capture: () => Object.freeze({ committedCanonicalSession, activeSessionResourceTable }),
  restore: (snapshot) => {
    committedCanonicalSession = snapshot.committedCanonicalSession;
    activeSessionResourceTable = snapshot.activeSessionResourceTable;
  }
});

const applyFiles = (filesData, { adoptCanonicalPayloads = false, resolveRecordInputs = true, targetState = state } = {}) => {
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
    targetState.linearRecordRows.splice(0);
    replaceLinearComparisonPlan(
      targetState.linearComparisonPlan,
      reconcileLinearComparisonPlan(targetState.linearComparisonPlan, targetState.linearSeqs)
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
      const rowByUid = new Map(targetState.linearRecordRows.map((entry) => [String(entry?.uid || ''), entry]));
      targetState.linearRecordRows.splice(
        0,
        targetState.linearRecordRows.length,
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
      (Array.isArray(filesData.linearComparisons) ? filesData.linearComparisons : [])
        .map((comparison) => [
          String(comparison?.id || ''),
          deserializeFile(comparison?.file)
        ])
        .filter(([id]) => id)
    );
    const planWithFiles = normalizeLinearComparisonPlan(targetState.linearComparisonPlan);
    planWithFiles.edges.forEach((edge) => {
      edge.file = comparisonFiles.get(edge.id) || null;
    });
    replaceLinearComparisonPlan(
      targetState.linearComparisonPlan,
      reconcileLinearComparisonPlan(planWithFiles, targetState.linearSeqs)
    );
    return { collapsedLinearSeqs };
  }

  targetState.linearSeqs.splice(0, targetState.linearSeqs.length, ...normalizeLinearSeqList([]));
  replaceLinearComparisonPlan(
    targetState.linearComparisonPlan,
    reconcileLinearComparisonPlan(targetState.linearComparisonPlan, targetState.linearSeqs)
  );
  return { collapsedLinearSeqs: false };
};

const reconcileDepthTrackStateAfterSessionFiles = () => {
  const circularDepthFiles = representativeDepthFiles(state.files.c_depth);
  const circularDepthCount = circularDepthFiles.some(Boolean) ? circularDepthFiles.length : 0;
  const linearRows = state.linearSeqs.map((seq) => depthFileSlotsFromValue(seq.depth));
  const linearDepthCount = state.mode.value === 'linear'
    ? depthTrackSessionWidth({
        rows: linearRows,
        depthTracks: state.adv.depth_tracks,
        slots: state.adv.linear_track_slots
      })
    : depthTrackMatrixWidth(linearRows);
  if (state.mode.value === 'linear' && linearDepthCount > 0) {
    state.linearSeqs.forEach((seq) => {
      seq.depth = padDepthFileSlots(seq.depth, linearDepthCount);
    });
  }

  const defaults = {
    depthColor: state.adv.depth_color,
    depthHeight: state.adv.depth_height,
    largeTickInterval: state.adv.depth_large_tick_interval,
    smallTickInterval: state.adv.depth_small_tick_interval,
    tickFontSize: state.adv.depth_tick_font_size
  };
  let normalizedTracks;
  if (state.mode.value === 'linear') {
    normalizedTracks = normalizeDepthTracks(state.adv.depth_tracks, state.adv);
    while (normalizedTracks.length < Math.max(1, linearDepthCount)) {
      normalizedTracks.push(normalizeDepthTrackConfig(null, normalizedTracks.length, state.adv));
    }
  } else {
    normalizedTracks = reconcileDepthTracksToFiles({
      files: circularDepthFiles,
      depthTracks: state.adv.depth_tracks,
      targetCount: Math.max(1, circularDepthCount),
      defaults
    });
  }
  state.adv.depth_tracks.splice(0, state.adv.depth_tracks.length, ...normalizedTracks);

  state.adv.circular_track_slots.splice(
    0,
    state.adv.circular_track_slots.length,
    ...dropInvalidManagedDepthSlots({
      slots: state.adv.circular_track_slots,
      activeCount: circularDepthCount
    })
  );
  syncDepthSlotLabels({
    slots: state.adv.circular_track_slots,
    depthTracks: state.adv.depth_tracks,
    activeCount: circularDepthCount
  });
  state.adv.circular_track_slots.splice(
    0,
    state.adv.circular_track_slots.length,
    ...applyCircularTrackOrderPlacements(
      state.adv.circular_track_slots,
      state.adv.nt,
      state.form.track_type,
      state.adv.circular_track_slots_axis_index
    )
  );

  state.adv.linear_track_slots.splice(
    0,
    state.adv.linear_track_slots.length,
    ...dropInvalidManagedDepthSlots({
      slots: state.adv.linear_track_slots,
      activeCount: linearDepthCount
    })
  );
  syncDepthSlotLabels({
    slots: state.adv.linear_track_slots,
    depthTracks: state.adv.depth_tracks,
    activeCount: linearDepthCount
  });
  state.adv.linear_track_slots.splice(
    0,
    state.adv.linear_track_slots.length,
    ...applyLinearTrackOrderPlacements(
      state.adv.linear_track_slots,
      state.adv.linear_track_slots_axis_index,
      state.adv.nt,
      state.form.linear_track_layout
    )
  );
};

const cloneLiveFileState = () => ({
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
  linearRecordRows: state.linearRecordRows.map((entry) => ({ ...entry })),
  linearComparisonPlan: {
    mode: state.linearComparisonPlan.mode,
    defaultSource: state.linearComparisonPlan.defaultSource,
    edges: state.linearComparisonPlan.edges.map((edge) => ({ ...edge }))
  }
});

const restoreLiveFileState = (snapshot) => {
  state.matchSequenceRegistry?.reset?.();
  Object.keys(state.files).forEach((key) => {
    state.files[key] = snapshot.files[key] ?? null;
  });
  state.circularRecordList.value = cloneJsonData(snapshot.circularRecordList);
  Object.assign(state.circularRecordDiscovery, snapshot.circularRecordDiscovery);
  state.linearSeqs.splice(0, state.linearSeqs.length, ...snapshot.linearSeqs);
  state.linearRecordRows.splice(0, state.linearRecordRows.length, ...snapshot.linearRecordRows);
  replaceLinearComparisonPlan(state.linearComparisonPlan, snapshot.linearComparisonPlan);
};

const captureSessionImportTransientState = () => ({
  semanticFileWatchersSuppressed: Boolean(
    state.semanticFileWatchersSuppressed.value
  ),
  skipCaptureBaseConfig: Boolean(state.skipCaptureBaseConfig.value),
  skipPositionReapply: Boolean(state.skipPositionReapply.value),
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
  fileLegendCaptions: Array.from(state.fileLegendCaptions.value || []),
  featureSearch: state.featureSearch.value,
  labelSearch: state.labelSearch.value,
  featureVisibilitySelectorCache: cloneJsonData(state.featureVisibilitySelectorCache),
  selectedFeatureIds: Array.from(state.selectedFeatureIds.value || []),
  selectedFeatureAnchorId: state.selectedFeatureAnchorId.value,
  featureSelectionStatus: state.featureSelectionStatus.value,
  featureSelectionSuppressNextClick: Boolean(
    state.featureSelectionSuppressNextClick.value
  ),
  featureSelectionDrag: cloneJsonData(state.featureSelectionDrag),
  labelReflowLastError: state.labelReflowLastError.value,
  labelOverrideBuildWarning: state.labelOverrideBuildWarning.value,
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
  state.fileLegendCaptions.value = new Set(snapshot.fileLegendCaptions);
  state.featureSearch.value = snapshot.featureSearch;
  state.labelSearch.value = snapshot.labelSearch;
  if (typeof state.replaceFeatureVisibilitySelectorCacheOwner === 'function') {
    state.replaceFeatureVisibilitySelectorCacheOwner(
      cloneJsonData(snapshot.featureVisibilitySelectorCache)
    );
  } else {
    replacePlainObject(
      state.featureVisibilitySelectorCache,
      cloneJsonData(snapshot.featureVisibilitySelectorCache)
    );
  }
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
  state.labelOverrideBuildWarning.value = snapshot.labelOverrideBuildWarning;
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
  state.skipPositionReapply.value = snapshot.skipPositionReapply;
};

const captureSessionImportSnapshot = () => ({
  config: cloneJsonData(buildConfigData()),
  ui: cloneJsonData(buildUiStateData()),
  files: cloneLiveFileState(),
  results: state.results.value,
  features: buildFeatureStateData(),
  editorState: buildEditorStateData(),
  orthogroupState: buildOrthogroupStateData(),
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
  importedComparisonIntent: cloneJsonData(state.importedComparisonIntent),
  errorLog: state.errorLog.value,
  resultPanelTab: state.resultPanelTab.value,
  transients: captureSessionImportTransientState()
});

const restoreSessionImportSnapshot = async (snapshot) => {
  state.sessionImportRollbackInProgress.value = true;
  try {
    state.semanticFileWatchersSuppressed.value = true;
    resetSessionBaseline();
    state.mode.value = snapshot.ui.mode === 'linear' ? 'linear' : 'circular';
    applyConfigData(snapshot.config, { resolveTrackPlacements: false });
    applyUiStateData(snapshot.ui);
    restoreLiveFileState(snapshot.files);
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
    Object.assign(
      state.importedComparisonIntent,
      createImportedComparisonIntentState(),
      cloneJsonData(snapshot.importedComparisonIntent)
    );
    state.skipCaptureBaseConfig.value = true;
    state.skipPositionReapply.value = true;
    applyResultsData(snapshot.results, snapshot.ui);
    applyFeatureStateData(snapshot.features);
    applyOrthogroupStateData(snapshot.orthogroupState);
    state.collinearGroups.value = snapshot.collinearGroups;
    applyEditorStateData(snapshot.editorState, { normalized: true });
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

const resetSessionBaseline = () => {
  activePreviewRuntime?.clearActiveRuntime?.();
  preservedCliOptions = null;
  committedCanonicalSession = null;
  activeSessionResourceTable = null;
  Object.assign(
    state.importedComparisonIntent,
    createImportedComparisonIntentState()
  );
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
  state.comparisonWarnings.value = [];
  applyFiles(null);
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
  clearObject(state.orthogroupNameOverrides);
  clearObject(state.orthogroupDescriptionOverrides);
  clearObject(state.orthogroupDormantOverrides);
  state.extractedFeatures.value = [];
  if (state.biologicalFeatures) state.biologicalFeatures.value = [];
  state.featureSelectorSafetyScope.value = [];
  state.featureRecordIds.value = [];
  state.selectedFeatureRecordIdx.value = 0;
  clearObject(state.featureColorOverrides);
  state.featureVisibilityManualRules.splice(0);
  clearObject(state.featureVisibilityOverrides);
  if (typeof state.replaceFeatureVisibilitySelectorCacheOwner === 'function') {
    state.replaceFeatureVisibilitySelectorCacheOwner({});
  } else {
    clearObject(state.featureVisibilitySelectorCache);
  }
  clearObject(state.featureStrokeOverrides);
  clearObject(state.labelTextFeatureOverrides);
  state.canonicalLabelOverrideRows.value = [];
  clearObject(state.labelTextBulkOverrides);
  clearObject(state.labelTextFeatureOverrideSources);
  clearObject(state.labelVisibilityOverrides);
  state.labelOverrideBuildWarning.value = '';
  state.generatedMode.value = 'circular';
  state.generatedLegendPosition.value = 'left';
  state.generatedMultiRecordCanvas.value = false;
  state.generatedCircularPlotTitlePosition.value = 'none';
};

export const buildUiStateData = ({ includePreviewNavigation = true } = {}) => {
  const ui = {
    title: String(state.sessionTitle.value || ''),
    mode: state.mode.value,
    canvasPadding: { ...state.canvasPadding },
    selectedResultIndex: state.selectedResultIndex.value,
    generatedLegendPosition: state.generatedLegendPosition.value,
    generatedMode: state.generatedMode.value,
    generatedMultiRecordCanvas: Boolean(state.generatedMultiRecordCanvas.value),
    generatedCircularPlotTitlePosition: normalizeCircularPlotTitlePosition(
      state.generatedCircularPlotTitlePosition.value
    ),
    layoutPreferences: cloneJsonData(state.layoutPreferences),
    featurePanelTab: state.featurePanelTab.value,
    cInputType: state.cInputType.value,
    lInputType: state.lInputType.value,
    losatProgram: state.losatProgram.value,
    downloadDpi: state.downloadDpi.value,
    autoLabelReflow: Boolean(state.autoLabelReflowEnabled.value),
    linearTypographyLinked: Boolean(state.linearTypographyLinked.value),
    paletteInstantPreviewEnabled: Boolean(state.paletteInstantPreviewEnabled.value),
    appliedPaletteName: state.appliedPaletteName.value,
    appliedPaletteColors: cloneColors(state.appliedPaletteColors.value),
    pendingPaletteName: state.pendingPaletteName.value,
    pendingPaletteColors: cloneColors(state.pendingPaletteColors.value),
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

export const applyUiStateData = (ui = {}, { restorePreviewNavigation = true } = {}) => {
  if (typeof ui.title === 'string') state.sessionTitle.value = ui.title;
  if (ui.mode) state.mode.value = ui.mode === 'linear' ? 'linear' : 'circular';
  if (ui.cInputType) state.cInputType.value = ui.cInputType;
  if (ui.lInputType) state.lInputType.value = ui.lInputType;
  if (ui.losatProgram) {
    const program = String(ui.losatProgram);
    state.losatProgram.value = ['blastn', 'tblastx', 'blastp'].includes(program) ? program : 'blastn';
  }
  if (ui.downloadDpi) state.downloadDpi.value = ui.downloadDpi;
  state.autoLabelReflowEnabled.value = Boolean(ui.autoLabelReflow);
  reconcileImportedLinearTypographyLink({
    adv: state.adv,
    linked: state.linearTypographyLinked,
    ui
  });
  state.paletteInstantPreviewEnabled.value = Boolean(ui.paletteInstantPreviewEnabled);
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

  restorePaletteStateFromSession(ui);
  restoreLayoutPreferences(ui);

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
    state.canvasPadding.top = Number(ui.canvasPadding.top) || 0;
    state.canvasPadding.right = Number(ui.canvasPadding.right) || 0;
    state.canvasPadding.bottom = Number(ui.canvasPadding.bottom) || 0;
    state.canvasPadding.left = Number(ui.canvasPadding.left) || 0;
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
    state.results.value = committedCount === logicalResults.length
      ? logicalResults
      : admitLegacyImportedResults(createLegacyImportResultSource(logicalResults));
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

export const buildFeatureStateData = () => ({
  extractedFeatures: sanitizeExtractedFeaturesForSession(state.extractedFeatures.value),
  biologicalFeatures: sanitizeExtractedFeaturesForSession(state.biologicalFeatures?.value),
  featureSelectorSafetyScope: cloneJsonData(state.featureSelectorSafetyScope.value),
  featureRecordIds: cloneJsonData(state.featureRecordIds.value),
  selectedFeatureRecordIdx: state.selectedFeatureRecordIdx.value,
  featureColorOverrides: cloneJsonData(state.featureColorOverrides),
  featureVisibilityManualRules: normalizeFeatureVisibilityRulesForSession(state.featureVisibilityManualRules),
  featureVisibilityOverrides: normalizeFeatureVisibilityOverridesForSession(state.featureVisibilityOverrides),
  labelTextFeatureOverrides: cloneJsonData(state.labelTextFeatureOverrides),
  labelOverrideRows: cloneJsonData(state.canonicalLabelOverrideRows.value),
  labelTextBulkOverrides: cloneJsonData(state.labelTextBulkOverrides),
  labelTextFeatureOverrideSources: cloneJsonData(state.labelTextFeatureOverrideSources),
  labelVisibilityOverrides: cloneJsonData(state.labelVisibilityOverrides)
});

export const applyFeatureStateData = (features = {}) => {
  state.extractedFeatures.value = Array.isArray(features.extractedFeatures)
    ? features.extractedFeatures
    : [];
  if (state.biologicalFeatures) {
    state.biologicalFeatures.value = Array.isArray(features.biologicalFeatures)
      ? features.biologicalFeatures
      : [];
  }
  state.featureSelectorSafetyScope.value = Array.isArray(features.featureSelectorSafetyScope)
    ? features.featureSelectorSafetyScope
    : [];
  state.featureRecordIds.value = Array.isArray(features.featureRecordIds)
    ? features.featureRecordIds
    : [];
  state.selectedFeatureRecordIdx.value = Number.isInteger(features.selectedFeatureRecordIdx)
    ? features.selectedFeatureRecordIdx
    : 0;
  replacePlainObject(state.featureColorOverrides, cloneJsonObject(features.featureColorOverrides));
  replaceFeatureVisibilityState(features);
  replacePlainObject(state.labelTextFeatureOverrides, cloneStringMap(features.labelTextFeatureOverrides));
  state.canonicalLabelOverrideRows.value = Array.isArray(features.labelOverrideRows)
    ? cloneJsonData(features.labelOverrideRows)
    : [];
  replacePlainObject(state.labelTextBulkOverrides, cloneStringMap(features.labelTextBulkOverrides));
  replacePlainObject(state.labelTextFeatureOverrideSources, cloneStringMap(features.labelTextFeatureOverrideSources));
  replacePlainObject(state.labelVisibilityOverrides, cloneJsonObject(features.labelVisibilityOverrides));
};

export const buildOrthogroupStateData = () => ({
  groups: Array.isArray(state.orthogroups.value) ? cloneJsonData(state.orthogroups.value) : [],
  selectedOrthogroupId: String(state.selectedOrthogroupId.value || ''),
  orthogroupNameOverrides: cloneStringMap(state.orthogroupNameOverrides),
  orthogroupDescriptionOverrides: cloneStringMap(state.orthogroupDescriptionOverrides),
  orthogroupDormantOverrides: normalizeOrthogroupDormantOverrides(state.orthogroupDormantOverrides)
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

const applySessionFeatureRecoveryPlan = (plan, { generationId = 'session-feature-recovery' } = {}) => {
  state.featureExtractionPending.value = false;

  if (plan?.status === 'recovered' || plan?.status === 'aligned') {
    if (plan.recoveredFeatureState) applyFeatureStateData(plan.recoveredFeatureState);
    if (plan.migratedEditorState) applyEditorStateData(plan.migratedEditorState);
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

const exportSessionDocument = async (
  titleOverride = null,
  { linearRecordCatalog = null, storedConfig, savedUi, isCurrent } = {}
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
  const logicalResults = serializeResults();
  const editorState = buildEditorStateData();
  if (logicalResults.length > 0) {
    if (!editorState.featureCatalog) throw sessionSaveRequiresGenerate();
    try {
      editorState.featureCatalog = validateFeatureCatalog(
        editorState.featureCatalog,
        logicalResults,
        { adopt: true, mode: state.mode.value }
      );
    } catch (error) {
      console.warn('Session feature catalog validation failed.', normalizeUserFacingError(error));
      throw sessionSaveRequiresGenerate();
    }
  } else {
    editorState.featureCatalog = null;
  }

  const {
    entries: losatEntries,
    validatedManifest,
    manifestValidated
  } = serializeLosatCache();
  const lastRunInvocation = state.lastRunInfo.value?.invocation;
  const exportableCliInvocation = isCliInvocationSessionExportable(lastRunInvocation)
    ? cloneJsonData(lastRunInvocation)
    : undefined;
  Object.assign(
    storedConfig.adv,
    normalizedArrowGeometryState(storedConfig.adv)
  );
  storedConfig.unmanagedConfigOverrides = await validateUnmanagedConfigOverrides({
    mode: state.mode.value,
    configOverrides: storedConfig.unmanagedConfigOverrides,
    requireUnmanagedOnly: true
  });
  let committed = isAdoptedCanonicalSession(committedCanonicalSession)
    ? committedCanonicalSession
    : cloneCanonicalSession(committedCanonicalSession);
  const settingsOnly = !committedCanonicalSession && logicalResults.length === 0
    && !hasBiologicalSessionInputs({ ...state.files, linearSeqs: state.linearSeqs });
  if (settingsOnly) validateCurrentWriterActiveConfig({ mode: state.mode.value, storedConfig });
  if (committed) {
    try {
      const adoptedCommitted = isAdoptedCanonicalSession(committed);
      const projected = projectCanonicalSessionRequest({
        ...committed,
        sessionResourceTable: adoptedCommitted ? activeSessionResourceTable : null,
        deferResourceContent: adoptedCommitted,
        adoptCanonicalPayloads: adoptedCommitted
      });
      validateCurrentWriterActiveConfig({
        mode: projected.mode,
        storedConfig
      });
    } catch (error) {
      console.warn('Session active configuration validation failed.', normalizeUserFacingError(error));
      throw recognizedCauseOr(error, diagnosticError('INPUT_INVALID', { field: 'config', reason: 'FIELDS' }));
    }
  }
  if (!committed && !settingsOnly) {
    const comparisonPlanSnapshot = state.mode.value === 'linear'
      ? resolveLinearComparisonPlan({
          plan: state.linearComparisonPlan,
          sequences: normalizeLinearSeqList(state.linearSeqs),
          layout: state.linearRecordLayoutEnabled?.value
            ? state.linearRecordRows
            : [],
          losatProgram: state.losatProgram?.value,
          blastpMode: state.losat?.blastp?.mode
        })
      : null;
    const activeFiles = await serializeActiveRenderFiles(
      state.mode.value,
      state,
      {
        comparisonPlan: comparisonPlanSnapshot,
        linearRecordCatalog
      }
    );
    committed = buildCanonicalRenderRequest({
      state,
      filesData: activeFiles,
      comparisonPlanSnapshot
    });
  }
  if (committed && committed.renderRequest.schema < CANONICAL_REQUEST_SCHEMA) {
    const promoted = {
      ...committed,
      renderRequest: promoteCanonicalRenderRequestToCurrent(
        committed.renderRequest,
        {
          featureCatalog: editorState.featureCatalog,
          legacyOrthogroupState: {
            groups: cloneJsonData(state.orthogroups.value || [])
          }
        }
      )
    };
    committed = isAdoptedCanonicalSession(committed)
      ? adoptRuntimeCanonicalSession(promoted)
      : promoted;
  }
  const canonical = await assembleSessionResources(state, committed);
  await validateSimilarityAlignmentResetReceipt(editorState.alignmentResetReceipt, canonical);
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
  const sessionData = {
    format: 'gbdraw-session',
    version: SESSION_VERSION,
    createdAt: new Date().toISOString(),
    title: resolvedTitle || undefined,
    config: storedConfig,
    ui: savedUi,
    renderRequest: canonical.renderRequest,
    resources: canonical.resources,
    webFiles: canonical.webFiles,
    results: logicalResults,
    ...(!settingsOnly ? { runMetadata: {
      ...(state.trackSlotResolvedGeometry.value
        ? { trackSlotGeometry: cloneJsonData(state.trackSlotResolvedGeometry.value) } : {}),
      annotationWarnings: cloneJsonData(state.annotationWarnings.value),
      ...(state.comparisonWarnings.value.length
        ? { comparisonWarnings: cloneJsonData(state.comparisonWarnings.value) } : {})
    } } : {}),
    features: {
      selectedFeatureRecordIdx: state.selectedFeatureRecordIdx.value,
      featureColorOverrides: cloneJsonData(state.featureColorOverrides),
      featureVisibilityManualRules: normalizeFeatureVisibilityRulesForSession(state.featureVisibilityManualRules),
      featureVisibilityOverrides: normalizeFeatureVisibilityOverridesForSession(state.featureVisibilityOverrides),
      labelTextFeatureOverrides: cloneJsonData(state.labelTextFeatureOverrides),
      labelOverrideRows: cloneJsonData(state.canonicalLabelOverrideRows.value),
      labelTextBulkOverrides: cloneJsonData(state.labelTextBulkOverrides),
      labelTextFeatureOverrideSources: cloneJsonData(state.labelTextFeatureOverrideSources),
      labelVisibilityOverrides: cloneJsonData(state.labelVisibilityOverrides)
    },
    editorState,
    orthogroupState: {
      selectedOrthogroupId: String(state.selectedOrthogroupId.value || ''),
      orthogroupNameOverrides: cloneStringMap(state.orthogroupNameOverrides),
      orthogroupDescriptionOverrides: cloneStringMap(state.orthogroupDescriptionOverrides)
    },
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
  };

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

const importSessionDocument = async (e, options = {}) => {
  const file = e.target.files[0];
  if (!file) return { status: 'skipped' };
  recordSessionLifecycleEvent('sessionSelection');

  const semanticFileWatchersSuppressedBeforeImport = Boolean(
    state.semanticFileWatchersSuppressed.value
  );
  const rollbackStateExtension = options?.rollbackState;
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
      applyLegacyConfigPayload(data);
      alert('Legacy configuration loaded. Save as a session to use the current format.');
      return { status: 'legacy' };
    }

    recordSessionLifecycleEvent('current-session-preflight-start');
    const preflight = await preflightSessionImport(data);
    await validateSimilarityAlignmentResetReceipt(
      data.editorState?.alignmentResetReceipt,
      { renderRequest: data.renderRequest, resources: data.resources }
    );
    recordSessionLifecycleEvent('current-session-preflight-end');
    data = preflight.data;
    const {
      sourceSessionVersion,
      canonicalProjection,
      restoredConfig,
      projectionResult,
      adoptedCanonicalSession,
      currentResourceTable,
      comparisonClassification,
      unmanagedConfigValidation
    } = preflight;
    if (restoredConfig && unmanagedConfigValidation) {
      restoredConfig.unmanagedConfigOverrides = await validateUnmanagedConfigOverrides(
        unmanagedConfigValidation
      );
    }
    const canonicalSession = Boolean(projectionResult);
    const settingsOnly = isSettingsOnlySessionDocument(data);
    const currentSchemaSession = sourceSessionVersion >= CURRENT_AUTHORITY_SESSION_MIN_VERSION;
    const committedMode = projectionResult?.renderState.mode;
    const savedCurrentWriterMode = sourceSessionVersion === SESSION_VERSION
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
    const candidateInputType = (candidateMode === 'linear'
      ? ui.lInputType : ui.cInputType) || canonicalProjection?.inputType || 'gb';
    const candidateFiles = {
      files: {}, linearSeqs: [], circularRecordList: { value: [] },
      circularRecordDiscovery: {}, linearReorderNotice: { value: '' },
      cInputType: { value: candidateMode === 'circular' ? candidateInputType : (ui.cInputType || 'gb') },
      linearRecordRows: cloneJsonData(restoredConfig?.linearRecordLayout?.rows || []),
      linearComparisonPlan: normalizeLinearComparisonPlan(restoredConfig?.linearComparisonPlan)
    };
    recordSessionLifecycleEvent('session-candidate-files-start');
    const { collapsedLinearSeqs } = applyFiles(
      canonicalSession ? projectionResult.restoredFiles : data.files,
      { adoptCanonicalPayloads: currentSchemaSession, resolveRecordInputs: !settingsOnly, targetState: candidateFiles }
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
    const features = canonicalSession
      ? {
          ...projectionResult.renderState.semanticFeatureState,
          ...(currentCatalogFeatureState || {}),
          ...artifactFeatureState
        }
      : (data.features || {});
    const catalogSequenceSources = currentSchemaSession
      ? (currentCatalogFeatureState?.sequenceSources || [])
      : [];
    recordSessionLifecycleEvent('session-candidate-sequences-start');
    const comparisonSourceAvailability = committedMode === 'circular'
      ? resolveCircularComparisonSequenceAvailability({
          files: candidateFiles.files,
          circularConservation: restoredConfig?.circularConservation || {}
        })
      : undefined;
    const catalogSequenceSourceCoverage = (
      currentSchemaSession
      && validatedSessionCatalog
    )
      ? analyzeCatalogSequenceSourceCoverage({
          mode: committedMode,
          catalogFeatureState: validatedSessionCatalog,
          renderRequest: data.renderRequest,
          comparisonSourceAvailability
        })
      : null;
    const missingCatalogSequenceSources = sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION
      || !catalogSequenceSourceCoverage?.complete;
    let restoredFileSequenceSources = [];
    let currentRecoveryError = null;
    if (missingCatalogSequenceSources && !settingsOnly) {
      recordStructuralMetric('sourceRecoveryCount');
      try {
        restoredFileSequenceSources = await buildRestoredMatchSequenceSources({
          mode: candidateMode,
          cInputType: candidateFiles.cInputType.value,
          lInputType: candidateInputType,
          files: candidateFiles.files,
          linearSeqs: candidateFiles.linearSeqs,
          circularConservation: restoredConfig?.circularConservation || {}
        });
      } catch (sequenceError) {
        currentRecoveryError = sequenceError;
        console.warn('Session match sequence preparation failed.', normalizeUserFacingError(sequenceError));
      }
    }

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
    const transformRestoredSessionSvg = (svg, { applyStrokes = true } = {}) => {
      const legendGroupsChanged = normalizeLegacyLegendEntryGroups(svg);
      let compositionChanged = false;
      if (
        svg.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) === null
        && svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE) === null
      ) {
        normalizeLegacyComposition(svg, {
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
        });
        compositionChanged = true;
      }
      const strokeCount = applyStrokes
        ? applyStrokeOverridesToSvg({
            svg,
            features: restoredFeatureState.extractedFeatures || [],
            legendStrokeOverrides: restoredEditorState?.legend?.strokeOverrides || {},
            featureStrokeOverrides: restoredEditorState?.featureStrokes?.overrides || {}
          })
        : 0;
      return legendGroupsChanged || compositionChanged || strokeCount > 0;
    };

    recordSessionLifecycleEvent('svg-admission-start');
    const committedImportedResults = currentSchemaSession && validatedSessionCatalog
      ? (() => {
          const catalogAdmission = admitFeatureCatalog(
            validatedSessionCatalog,
            logicalImportedResults,
            { adopt: true, mode: committedMode }
          );
          return admitCurrentSessionResults(
            createCurrentSessionResultSource(logicalImportedResults, catalogAdmission),
            { mutationPlan: createEmptySvgMutationPlan(logicalImportedResults.length) }
          );
        })()
      : admitLegacyImportedResults(
          createLegacyImportResultSource(logicalImportedResults),
          { transformSvg: transformRestoredSessionSvg }
        );
    recordSessionLifecycleEvent('svg-admission-end');

    let legacyFeatureRecoveryPlan = null;
    if (sourceSessionVersion < CURRENT_AUTHORITY_SESSION_MIN_VERSION) {
      try {
        legacyFeatureRecoveryPlan = await buildSessionFeatureRecoveryPlan({
          snapshot: {
            mode: candidateMode, cInputType: candidateFiles.cInputType.value,
            lInputType: candidateInputType, files: candidateFiles.files,
            linearSeqs: candidateFiles.linearSeqs, results: committedImportedResults,
            selectedResultIndex: ui.selectedResultIndex || 0,
            featureState: features, editorState: restoredEditorState
          },
          featureVisibilityTsv: serializeFeatureVisibilityRules(
            features.featureVisibilityManualRules || features.featureVisibilityRules || []
          )
        });
      } catch (error) {
        legacyFeatureRecoveryPlan = { status: 'failed', warning: normalizeUserFacingError(error).summary };
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
    state.paletteInstantPreviewEnabled.value = Boolean(ui.paletteInstantPreviewEnabled);
    state.labelOverrideBuildWarning.value = '';
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

    if (restoredConfig) {
      state.suppressCircularMultiRecordDefaults.value = shouldSuppressCircularMultiRecordDefaults(
        restoredConfig.form
      );
      applyConfigData(restoredConfig, { resolveTrackPlacements: !settingsOnly });
    }
    const canonicalLinearLayout = canonicalProjection?.config?.linearRecordLayout;
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
    reconcileImportedLinearTypographyLink({
      adv: state.adv,
      linked: state.linearTypographyLinked,
      ui
    });
    restorePaletteStateFromSession(ui);
    restoreLayoutPreferences(ui, {
      projected: canonicalSession ? canonicalProjection?.layoutPreferences : null
    });

    restoreLiveFileState({
      files: candidateFiles.files,
      linearSeqs: candidateFiles.linearSeqs,
      circularRecordList: candidateFiles.circularRecordList.value,
      circularRecordDiscovery: candidateFiles.circularRecordDiscovery,
      linearRecordRows: candidateFiles.linearRecordRows,
      linearComparisonPlan: candidateFiles.linearComparisonPlan
    });
    restoreImportedComparisonIntent(
      state.importedComparisonIntent,
      comparisonClassification,
      restoredConfig?.importedComparisonResolution
    );
    if (!settingsOnly) reconcileDepthTrackStateAfterSessionFiles();
    if (canonicalSession) {
      applyLosatCache(
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
    state.skipPositionReapply.value = false;

    applyFeatureStateData(features);
    if (currentSchemaSession && currentCatalogFeatureState) {
      state.collinearGroups.value = currentCatalogFeatureState.collinearGroups;
      synchronizeRestoredFeatureSummaryStatus({ generationId: 'session-load' });
    }
    state.matchSequenceRegistry?.reset?.([
      ...catalogSequenceSources,
      ...restoredFileSequenceSources
    ]);

    applyOrthogroupStateData(
      canonicalSession
        ? {
            ...projectionResult.artifactState.orthogroupState,
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
    applyEditorStateData(restoredEditorState, { normalized: currentSchemaSession });
    if (legacyFeatureRecoveryPlan) {
      applySessionFeatureRecoveryPlan(legacyFeatureRecoveryPlan, { generationId: 'session-load' });
    }

    const desiredResultIndex = (
      Number.isInteger(ui.selectedResultIndex) && ui.selectedResultIndex >= 0
    )
      ? Math.min(ui.selectedResultIndex, Math.max(0, committedImportedResults.length - 1))
      : 0;
    if (!options.isCurrent()) throw new Error('Session loading was canceled.');
    state.skipCaptureBaseConfig.value = true;
    state.skipPositionReapply.value = true;
    recordSessionLifecycleEvent('session-candidate-adopted');
    recordSessionLifecycleEvent('preview-mount-start');
    applyResultsData(committedImportedResults, ui);
    state.annotationWarnings.value = cloneJsonData(
      projectionResult?.artifactState?.runMetadata?.annotationWarnings || []
    );
    state.comparisonWarnings.value = cloneJsonData(
      projectionResult?.artifactState?.runMetadata?.comparisonWarnings || []
    );
    state.trackSlotResolvedGeometry.value = cloneJsonData(
      projectionResult?.artifactState?.runMetadata?.trackSlotGeometry ?? null
    );
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
    await options.afterImport?.({ status: 'ok', decompressedCharacters: candidate.characters, isCurrent: options.isCurrent });
    if (!options.isCurrent()) throw new Error('Session loading was canceled.');
    alert('Session loaded successfully!');
    return {
      status: 'ok',
      data,
      decompressedCharacters: candidate.characters,
      degradedRecovery: Boolean(currentRecoveryError),
      comparisonDisposition: state.importedComparisonIntent.disposition
    };
  } catch (err) {
    if (err?.name === 'AbortError' && !commitStarted) return { status: 'canceled' };
    const error = normalizeUserFacingError(err, { stage: 'request-validation' });
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


let sessionSaveInFlight = null;
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

export const exportSession = (titleOverride = null, options = {}) => {
  if (sessionSaveInFlight) {
    recordSessionLifecycleEvent('session-save-joined');
    return sessionSaveInFlight.promise;
  }
  const busy = sessionOperationAvailability('save');
  if (busy) return Promise.resolve(busy);
  const operation = { canceled: false, promise: null };
  const previousAlert = state.errorLog.value;
  const isCurrent = () => sessionSaveInFlight === operation && !operation.canceled;
  operation.promise = Promise.resolve().then(async () => {
    const busy = sessionOperationAvailability('save');
    if (busy) return busy;
    if (!isCurrent()) return { status: 'canceled' };
    const title = options.resolveTitle ? options.resolveTitle() : titleOverride;
    if (title === null && options.resolveTitle) return;
    state.sessionSavePending.value = true;
    recordSessionLifecycleEvent('session-save-pending-published');
    // Only mutable draft configuration/navigation is copied. Adopted biological
    // payloads, resources, catalogs, caches and Results retain their existing owner.
    const activeConfig = buildConfigData();
    validateCurrentWriterActiveConfig({ mode: state.mode.value, storedConfig: activeConfig });
    if (!committedCanonicalSession
      && hasBiologicalSessionInputs({ ...state.files, linearSeqs: state.linearSeqs })) {
      assertActiveModeInputs();
    }
    const storedConfig = cloneJsonData(activeConfig);
    const savedUi = {
      mode: state.mode.value,
      zoom: state.zoom.value,
      canvasPan: { x: state.canvasPan.x, y: state.canvasPan.y },
      canvasPadding: { ...state.canvasPadding },
      selectedResultIndex: state.selectedResultIndex.value,
      generatedLegendPosition: state.generatedLegendPosition.value,
      generatedMultiRecordCanvas: Boolean(state.generatedMultiRecordCanvas.value),
      generatedCircularPlotTitlePosition: normalizeCircularPlotTitlePosition(
        state.generatedCircularPlotTitlePosition.value
      ),
      layoutPreferences: cloneJsonData(state.layoutPreferences),
      featurePanelTab: state.featurePanelTab.value,
      cInputType: state.cInputType.value,
      lInputType: state.lInputType.value,
      downloadDpi: state.downloadDpi.value,
      autoLabelReflow: Boolean(state.autoLabelReflowEnabled.value),
      linearTypographyLinked: Boolean(state.linearTypographyLinked.value),
      paletteInstantPreviewEnabled: Boolean(state.paletteInstantPreviewEnabled.value),
      appliedPaletteName: state.appliedPaletteName.value,
      appliedPaletteColors: cloneColors(state.appliedPaletteColors.value),
      pendingPaletteName: state.pendingPaletteName.value,
      pendingPaletteColors: cloneColors(state.pendingPaletteColors.value)
    };
    const prepared = await options.beforeExport?.();
    if (!isCurrent()) return { status: 'canceled' };
    return exportSessionDocument(title, {
      ...options, ...prepared, storedConfig, savedUi, isCurrent
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
  const busy = sessionOperationAvailability('load');
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
