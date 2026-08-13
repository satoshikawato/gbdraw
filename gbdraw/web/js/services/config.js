import { state, normalizeLinearSeqList, collapseEmptyLinearSeqList } from '../state.js';
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
  buildCircularTrackSlotPayload,
  CIRCULAR_TRACK_RENDERERS,
  circularTrackAxisIndexForEnabledSlots,
  clampCircularTrackAxisIndex,
  inferLegacyAxisIndexFromFeature,
  migrateLegacyCircularTrackSlot,
  migrateLegacyCircularTrackSlotSpec,
  parseCircularTrackSlotSpec,
  normalizeCircularTrackSlots
} from '../app/circular-track-slots.js';
import {
  applyLinearTrackOrderPlacements,
  buildLinearTrackSlotPayload,
  clampLinearTrackAxisIndex,
  LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  LINEAR_TRACK_RENDERERS,
  LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  migrateLinearTrackSlotsToCurrentSchema,
  normalizeLinearTrackSlots,
  linearTrackAxisIndexForEnabledSlots,
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
import { validateTrackSlotBindingInvariants } from '../app/track-slot-validation.js';
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
import { isCliInvocationSessionExportable } from '../app/run-info.js';
import {
  normalizeCircularPlotTitlePosition,
  normalizeLinearPlotTitlePosition
} from '../app/plot-title-position.js';
import {
  migrateLegacyLayoutPreferences,
  replaceLayoutPreferences
} from '../app/layout-preferences.js';
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
  buildCanonicalRenderRequest,
  projectCanonicalSessionRequest
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
import { base64ToBytes, bytesToBase64, readFileBytes } from './file-content-cache.js';
import { normalizeLogicalResults } from './result-normalization.js';
import {
  ingestSvgResults,
  isCommittedSvgResult
} from './svg-result-ingestion.js';
import {
  featureStateFromCatalog,
  validateFeatureCatalog
} from './feature-catalog.js';
import {
  buildOrthogroupFeatureIndex,
  enrichFeaturesWithOrthogroups
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
  confirmLargeSessionBlob,
  readSessionText
} from './session-file.js';
import { downloadBlob } from './text-download.js';
import { normalizeAnnotationSets } from '../app/annotations/state.js';
import { applySpecificRuleProvenance } from '../app/specific-color-rules.js';
import { prepareCandidateRenderCommit } from '../app/candidate-render.js';
import { applyStrokeOverridesToSvg } from '../app/legend/stroke-actions.js';
import { normalizeLegacyLegendEntryGroups } from './svg-result-normalization.js';
import {
  COMPOSITION_METADATA_ATTRIBUTE,
  COMPOSITION_SCHEMA_ATTRIBUTE,
  normalizeLegacyComposition
} from '../app/legend-layout/composition-actions.js';
import {
  LOSAT_DERIVED_CACHE_SCHEMA,
  classifyRawLosatCacheEntry,
  createLegacyProteinCandidateEnvelope,
  emptyProteinIdentityManifest,
  isCurrentRawLosatCacheEntry,
  isLosatDerivedCacheEntry,
  normalizeLegacyProteinCandidateEnvelope,
  proteinRuntimeIdSets,
  rawProteinTextMatchesBindings,
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
  projectArtifactState,
  projectDocumentMetadata,
  projectWebOnlyEditorMetadata,
  validateSessionAuthorityInventory
} from './session-authority.js';
import {
  migratePersistedCircularMultiRecordSizeMode,
  migratePersistedLinearLabelPlacement,
  migratePersistedLinearTrackLayout,
  migratePersistedWebStateFieldNames,
  requireCurrentCircularMultiRecordSizeMode,
  requireCurrentLinearLabelPlacement,
  requireCurrentLinearTrackLayout,
  requireCurrentWebStateFieldNames
} from '../app/current-option-values.js';

const { nextTick } = window.Vue;

export const SESSION_VERSION = 40;
const LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION = 32;
const SUPPORTED_SESSION_VERSIONS = new Set([
  27, 28, 29, 30, 31, 32, 33, 39, SESSION_VERSION
]);
const CURRENT_ARTIFACT_SESSION_MIN_VERSION = 39;
const LOSAT_DERIVED_CACHE_LIMIT = 16;
const SESSION_FEATURE_CATALOG_SAVE_ERROR =
  'Generate again before using Save Session. The current results are missing compatible feature metadata.';
const SESSION_ACTIVE_DRAFT_SAVE_ERROR =
  'Generate again before using Save Session. The active Custom Track draft has changed since the committed render.';
const CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 4;
const LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 3;
const OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS = [
  'gapAfter',
  'gap_after',
  'innerRadius',
  'inner_radius',
  'outerRadius',
  'outer_radius',
  'placement',
  'spacing',
  'strict',
  'compress',
  'reserve'
];
const OBSOLETE_CIRCULAR_TRACK_SLOT_PARAM_KEYS = [
  'side',
  'radius',
  'width',
  'spacing',
  'inner_gap_px',
  'outer_gap_px',
  'strict',
  'compress',
  'reserve'
];

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

    // 5. Update value only if types match or initial value is null
    if (typeof targetValue === typeof sourceValue || targetValue === null) {
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

const findObsoleteCircularTrackSlotShape = (slots) => {
  if (!Array.isArray(slots)) return null;

  for (let index = 0; index < slots.length; index += 1) {
    const slot = slots[index];
    if (!slot || typeof slot !== 'object' || Array.isArray(slot)) continue;

    for (const key of OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS) {
      if (Object.prototype.hasOwnProperty.call(slot, key)) {
        return `circular_track_slots[${index}].${key}`;
      }
    }

    const params = slot.params;
    if (!params || typeof params !== 'object' || Array.isArray(params)) continue;
    for (const key of OBSOLETE_CIRCULAR_TRACK_SLOT_PARAM_KEYS) {
      if (Object.prototype.hasOwnProperty.call(params, key)) {
        return `circular_track_slots[${index}].params.${key}`;
      }
    }
  }

  return null;
};

const validateImportedCircularTrackSlots = (configData = {}, { depthTrackCount = null } = {}) => {
  const adv = configData && typeof configData === 'object' ? configData.adv : null;
  if (!adv || typeof adv !== 'object' || Array.isArray(adv)) return;
  if (!Object.prototype.hasOwnProperty.call(adv, 'circular_track_slots')) return;

  if (adv.circular_track_slots_schema_version !== CIRCULAR_TRACK_SLOT_SCHEMA_VERSION) {
    throw new Error(
      `Custom Track Slots use an obsolete schema. Recreate the slots with schema version ${CIRCULAR_TRACK_SLOT_SCHEMA_VERSION}.`
    );
  }

  const obsoletePath = findObsoleteCircularTrackSlotShape(adv.circular_track_slots);
  if (obsoletePath) {
    throw new Error(
      `Custom Track Slots use obsolete field '${obsoletePath}'. Use slot-level radius, width, inner_gap_px, outer_gap_px, side, and z fields.`
    );
  }
  validateTrackSlotBindingInvariants(adv.circular_track_slots, {
    modeLabel: 'Circular',
    layoutKind: 'circular',
    supportedRenderers: CIRCULAR_TRACK_RENDERERS,
    supportedSides: ['inside', 'outside', 'overlay'],
    anchorlessRenderers: ['ticks', 'spacer'],
    depthTrackCount
  });
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
  return {
    ...migratedNames,
    ...(form === undefined ? {} : { form }),
    ...(adv === undefined ? {} : { adv })
  };
};

const validateImportedLinearTrackSlots = (configData = {}, { depthTrackCount = null } = {}) => {
  const adv = configData && typeof configData === 'object' ? configData.adv : null;
  if (!adv || typeof adv !== 'object' || Array.isArray(adv)) return;
  if (!Object.prototype.hasOwnProperty.call(adv, 'linear_track_slots')) return;

  if (adv.linear_track_slots_schema_version !== LINEAR_TRACK_SLOT_SCHEMA_VERSION) {
    throw new Error(
      `Custom Track Slots use an obsolete schema. Recreate the slots with schema version ${LINEAR_TRACK_SLOT_SCHEMA_VERSION}.`
    );
  }
  if (!Array.isArray(adv.linear_track_slots)) {
    throw new Error('Custom Track Slots must be an array.');
  }
  validateTrackSlotBindingInvariants(adv.linear_track_slots, {
    modeLabel: 'Linear',
    layoutKind: 'linear',
    supportedRenderers: LINEAR_TRACK_RENDERERS,
    supportedSides: ['above', 'below', 'overlay'],
    anchorlessRenderers: ['spacer'],
    depthTrackCount
  });
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

const normalizePairwiseMatchStyle = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['ribbon', 'curve'].includes(normalized) ? normalized : 'ribbon';
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

const cloneLosatForConfig = () => {
  const cloned = cloneJsonData(state.losat || {});
  if (cloned.blastp && typeof cloned.blastp === 'object' && !Array.isArray(cloned.blastp)) {
    delete cloned.blastp.collinearAnchorMode;
  }
  return cloned;
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
  losat: cloneLosatForConfig(),
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
  modeProfiles: state.modeProfileStateManager?.exportState?.(),
  linearRecordLayout: {
    enabled: Boolean(state.linearRecordLayoutEnabled.value),
    recordGap: Number(state.linearRecordGap.value) || 0,
    rows: (state.linearRecordRows || []).map((entry) => ({
      uid: String(entry?.uid || ''),
      row: Number(entry?.row) || 1
    }))
  },
  linearComparisonPlan: serializeLinearComparisonPlan(state.linearComparisonPlan),
  webEdits: {
    orthogroupNameOverrides: cloneStringMap(state.orthogroupNameOverrides),
    orthogroupDescriptionOverrides: cloneStringMap(state.orthogroupDescriptionOverrides)
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
  featureCatalog: cloneJsonValue(state.featureCatalog?.value, null)
});

const normalizeEditorStateData = (editorState = {}) => {
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
    featureCatalog: isPlainObject(source.featureCatalog)
      ? cloneJsonData(source.featureCatalog)
      : null
  };
};

const replacePlainObject = (target, source) => {
  Object.keys(target).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = value;
  });
};

export const applyEditorStateData = (editorState = {}, { trusted = false } = {}) => {
  const normalized = trusted
    ? cloneJsonData(editorState)
    : normalizeEditorStateData(editorState);

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
    state.featureCatalog.value = cloneJsonValue(normalized.featureCatalog, null);
  }
};

const normalizeSessionData = (data) => {
  if (!isPlainObject(data) || data.format !== 'gbdraw-session') {
    throw new Error('Invalid session file.');
  }
  const version = data.version;
  if (!Number.isInteger(version)) {
    throw new Error('Session version is required and must be an integer.');
  }
  if (version > SESSION_VERSION) {
    throw new Error(`Session version ${version} is newer than this gbdraw supports (${SESSION_VERSION}).`);
  }
  if (!SUPPORTED_SESSION_VERSIONS.has(version)) {
    throw new Error(`Unsupported session version: ${version}.`);
  }
  if (version === SESSION_VERSION && Object.prototype.hasOwnProperty.call(data, 'files')) {
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

export const validateSessionLosatArtifacts = (data, sourceSessionVersion) => {
  if (sourceSessionVersion < CURRENT_ARTIFACT_SESSION_MIN_VERSION) return;
  const rawEntries = sessionArtifactEntries(data, 'losatCache');
  const derivedEntries = sessionArtifactEntries(data, 'losatDerivedCache');
  const manifest = data.proteinIdentityManifest;
  rejectInvalidLosatCacheKeys(rawEntries, 'LOSAT', { requireKey: true });
  rejectInvalidLosatCacheKeys(derivedEntries, 'derived LOSATP');

  if (!validateProteinIdentityManifest(manifest)) {
    throw new Error(
      `Session version ${sourceSessionVersion} requires a valid schema-2 protein manifest.`
    );
  }
  for (const entry of rawEntries) {
    const classification = classifyRawLosatCacheEntry(entry);
    if (!['protein-current', 'nucleotide-current'].includes(classification)) {
      throw new Error(
        `Session version ${sourceSessionVersion} contains a non-current raw LOSAT entry.`
      );
    }
    if (classification !== 'protein-current') continue;
    const ids = proteinRuntimeIdSets(
      manifest,
      entry.queryRecordInstanceKey,
      entry.subjectRecordInstanceKey
    );
    if (
      !validateProteinRawEntryReferences(entry, manifest) ||
      !ids ||
      !rawProteinTextMatchesBindings(entry.text, ids.query, ids.subject)
    ) {
      throw new Error(
        `Session version ${sourceSessionVersion} contains an unresolved protein raw entry.`
      );
    }
  }
  if (
    derivedEntries.some(
      (entry) => !validateDerivedProteinReferences(entry, manifest)
    )
  ) {
    throw new Error(
      `Session version ${sourceSessionVersion} contains an invalid derived LOSATP entry.`
    );
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
  const readsLegacyOptionValues = sourceSessionVersion < SESSION_VERSION;
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

const CURRENT_WRITER_DRAFT_ADV_FIELDS = Object.freeze([
  'circular_track_slots_enabled',
  'circular_track_slots',
  'circular_track_slots_axis_index',
  'circular_track_slots_schema_version',
  'linear_track_slots_enabled',
  'linear_track_slots',
  'linear_track_slots_axis_index',
  'linear_track_slots_schema_version',
  'feature_width_circular',
  'depth_width_circular',
  'gc_content_width_circular',
  'gc_content_radius_circular',
  'gc_skew_width_circular',
  'gc_skew_radius_circular',
  'comparison_height',
  'pairwise_match_style',
  'min_bitscore',
  'evalue',
  'identity',
  'alignment_length'
]);

export const overlayCurrentWriterDraftConfig = (projectedConfig, storedConfig) => {
  const restored = cloneJsonData(projectedConfig);
  if (!isPlainObject(storedConfig)) return restored;
  if (isPlainObject(storedConfig.adv)) {
    if (!isPlainObject(restored.adv)) restored.adv = {};
    CURRENT_WRITER_DRAFT_ADV_FIELDS.forEach((field) => {
      if (Object.prototype.hasOwnProperty.call(storedConfig.adv, field)) {
        restored.adv[field] = cloneJsonData(storedConfig.adv[field]);
      }
    });
  }
  if (isPlainObject(storedConfig.linearComparisonPlan)) {
    if (isPlainObject(storedConfig.losat)) {
      restored.losat = cloneJsonData(storedConfig.losat);
    }
    if (Object.prototype.hasOwnProperty.call(storedConfig, 'losatProgram')) {
      restored.losatProgram = cloneJsonData(storedConfig.losatProgram);
    }
  }
  return restored;
};

const stableJsonValue = (value) => {
  if (Array.isArray(value)) return value.map((item) => stableJsonValue(item));
  if (!isPlainObject(value)) return value;
  return Object.fromEntries(
    Object.keys(value).sort().map((key) => [key, stableJsonValue(value[key])])
  );
};

const trackPayloadSignature = (payload) => JSON.stringify(stableJsonValue(payload));

export const validateCurrentWriterActiveDraft = ({
  mode,
  projectedConfig,
  storedConfig
}) => {
  const storedAdv = storedConfig?.adv;
  const projectedAdv = projectedConfig?.adv;
  if (!isPlainObject(storedAdv) || !isPlainObject(projectedAdv)) {
    throw new Error('Current session is missing its Web track draft state.');
  }

  const circular = mode === 'circular';
  const enabledField = circular
    ? 'circular_track_slots_enabled'
    : 'linear_track_slots_enabled';
  if (!Boolean(storedAdv[enabledField])) return;

  const slotsField = circular ? 'circular_track_slots' : 'linear_track_slots';
  const axisField = circular
    ? 'circular_track_slots_axis_index'
    : 'linear_track_slots_axis_index';
  const storedDraft = Array.isArray(storedAdv[slotsField]) ? storedAdv[slotsField] : [];
  const storedEnabled = storedDraft.filter((slot) => slot?.enabled !== false);
  const projectedEnabled = Array.isArray(projectedAdv[slotsField])
    ? projectedAdv[slotsField]
    : [];
  const storedPayloads = circular
    ? storedEnabled.map((slot) => buildCircularTrackSlotPayload(
        slot,
        storedAdv.nt || 'GC',
        storedConfig?.form?.track_type || 'tuckin'
      ))
    : storedEnabled.map((slot) => buildLinearTrackSlotPayload(slot));
  const projectedPayloads = circular
    ? projectedEnabled.map((slot) => buildCircularTrackSlotPayload(
        slot,
        projectedAdv.nt || 'GC',
        projectedConfig?.form?.track_type || 'tuckin'
      ))
    : projectedEnabled.map((slot) => buildLinearTrackSlotPayload(slot));
  const storedAxis = circular
    ? circularTrackAxisIndexForEnabledSlots(storedDraft, storedAdv[axisField])
    : linearTrackAxisIndexForEnabledSlots(storedDraft, storedAdv[axisField]);
  const projectedAxis = projectedAdv[axisField];
  if (
    trackPayloadSignature(storedPayloads) !== trackPayloadSignature(projectedPayloads)
    || Number(storedAxis) !== Number(projectedAxis)
  ) {
    throw new Error(
      `Current ${mode} Web track draft does not match the committed render request.`
    );
  }
};

const validateCurrentWriterFeatureCatalog = (data) => {
  const results = normalizeLogicalResults(
    (Array.isArray(data.results) ? data.results : []).map((result, index) => ({
      name: result?.name || `Result ${index + 1}`,
      content: result?.content || ''
    }))
  );
  const catalog = data.editorState?.featureCatalog ?? null;
  if (catalog === null) return;
  data.editorState.featureCatalog = validateFeatureCatalog(catalog, results);
};

const preflightSessionImport = (rawData) => {
  const sourceSessionVersion = rawData.version;
  validateSessionAuthorityInventory(rawData, sourceSessionVersion);
  const normalizedData = normalizeSessionData(rawData);
  if (sourceSessionVersion === SESSION_VERSION) {
    validateCurrentWriterFeatureCatalog(normalizedData);
  }
  migrateImportedLinearTrackSlots(normalizedData.config, sourceSessionVersion);
  const promotedData = (
    sourceSessionVersion >= 31 &&
    Number(normalizedData.renderRequest?.schema) === 2
    )
    ? promoteGallerySessionToCurrent(normalizedData)
    : normalizedData;
  validateSessionLosatArtifacts(promotedData, sourceSessionVersion);
  const data = migrateSessionDataToCurrent(promotedData, sourceSessionVersion);
  const canonicalProjection = sourceSessionVersion >= 31
    ? projectCanonicalSessionRequest({
        renderRequest: data.renderRequest,
        resources: data.resources,
        webFiles: data.webFiles,
        legacyFiles: data.files,
        storedConfig: data.config,
        fileBindings: data.cliInvocation?.fileBindings,
        linearTrackSlotSchemaVersion: sourceSessionVersion <= LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION
          ? LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION
          : LINEAR_TRACK_SLOT_SCHEMA_VERSION,
        repairInvalidComparisonHeight: sourceSessionVersion >= 31 && sourceSessionVersion <= 33
      })
    : null;
  let restoredConfig = canonicalProjection
    ? {
        ...restoreStoredNonCanonicalConfig(
          canonicalProjection.config,
          data.config,
          { hasCanonicalProteinPipeline: Boolean(canonicalProjection.pipelineState) }
        ),
        rules: applySpecificRuleProvenance(
          canonicalProjection.config.rules,
          data.config?.rules
        )
      }
    : data.config;
  if (sourceSessionVersion < SESSION_VERSION && restoredConfig) {
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
  if (canonicalProjection && sourceSessionVersion === SESSION_VERSION) {
    validateCurrentWriterActiveDraft({
      mode: canonicalProjection.mode,
      projectedConfig: canonicalProjection.config,
      storedConfig: data.config
    });
    restoredConfig = overlayCurrentWriterDraftConfig(restoredConfig, data.config);
  }
  if (!canonicalProjection && restoredConfig) {
    hydrateMissingMultiRecordPositionsFromCliInvocation(restoredConfig, data.cliInvocation);
  }
  if (restoredConfig) {
    const projectedDepthTrackCount = Array.isArray(canonicalProjection?.config?.adv?.depth_tracks)
      ? canonicalProjection.config.adv.depth_tracks.length
      : null;
    validateImportedCircularTrackSlots(restoredConfig, {
      depthTrackCount: canonicalProjection?.mode === 'circular' &&
        restoredConfig.adv?.circular_track_slots_enabled
        ? projectedDepthTrackCount
        : null
    });
    validateImportedLinearTrackSlots(restoredConfig, {
      depthTrackCount: canonicalProjection?.mode === 'linear' &&
        restoredConfig.adv?.linear_track_slots_enabled
        ? projectedDepthTrackCount
        : null
    });
  }
  const projectionResult = canonicalProjection
    ? (() => {
        const artifactState = projectArtifactState(data);
        if (canonicalProjection.pipelineState) {
          artifactState.orthogroupState = {
            ...(artifactState.orthogroupState || {}),
            selectedOrthogroupAlignmentFeature:
              canonicalProjection.pipelineState.selectedOrthogroupAlignmentFeature
          };
        }
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
          restoredFiles: canonicalProjection.files
        };
      })()
    : null;
  return {
    data,
    sourceSessionVersion,
    canonicalProjection,
    restoredConfig,
    projectionResult
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

const restoreLayoutPreferences = (ui = {}, { preserveActive = false } = {}) => {
  const activeBeforeRestore = {
    legend: state.form.legend,
    plotTitlePosition: state.adv.plot_title_position
  };
  const migrationUi = (
    !isPlainObject(ui.layoutPreferences) &&
    state.mode.value === 'linear' &&
    !hasStoredLayoutValue(ui.linearLegendPosition) &&
    hasStoredLayoutValue(ui.legend)
  )
    ? { ...ui, linearLegendPosition: ui.legend }
    : ui;
  const migrated = migrateLegacyLayoutPreferences(migrationUi, {
    mode: state.mode.value,
    multiRecord: Boolean(state.form.multi_record_canvas),
    activeLegend: activeBeforeRestore.legend,
    activePlotTitlePosition: activeBeforeRestore.plotTitlePosition
  });
  // Current-session render requests describe the last generated artifact, while
  // stored layout preferences own the editor's active semantic position. Use the
  // projected request only when neither the current nor legacy owner is present.
  if (preserveActive && !hasStoredLayoutPreferences(ui)) {
    if (state.mode.value === 'linear') {
      migrated.linear = {
        legend: normalizeLegendPosition(activeBeforeRestore.legend, 'bottom'),
        plotTitlePosition: normalizeLinearPlotTitlePosition(
          activeBeforeRestore.plotTitlePosition
        )
      };
    } else {
      const key = state.form.multi_record_canvas ? 'multi' : 'single';
      migrated.circular[key] = {
        legend: normalizeLegendPosition(activeBeforeRestore.legend, 'left'),
        plotTitlePosition: normalizeCircularPlotTitlePosition(
          activeBeforeRestore.plotTitlePosition
        )
      };
    }
  }
  replaceLayoutPreferences(
    state.layoutPreferences,
    migrated
  );
};

export const applyConfigData = (data) => {
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
  if (data.adv) safeDeepMerge(state.adv, data.adv);
  state.annotationSets.splice(
    0,
    state.annotationSets.length,
    ...normalizeAnnotationSets(data.annotationSets)
  );
  const linearLayout = data.linearRecordLayout && typeof data.linearRecordLayout === 'object'
    ? data.linearRecordLayout
    : null;
  state.linearRecordLayoutEnabled.value = Boolean(linearLayout?.enabled);
  const linearRecordGap = Number(linearLayout?.recordGap);
  state.linearRecordGap.value = Number.isFinite(linearRecordGap) && linearRecordGap >= 0
    ? linearRecordGap
    : 24;
  state.linearRecordRows.splice(
    0,
    state.linearRecordRows.length,
    ...(Array.isArray(linearLayout?.rows) ? linearLayout.rows : [])
      .map((entry) => ({ uid: String(entry?.uid || ''), row: Number(entry?.row) }))
      .filter((entry) => entry.uid && Number.isInteger(entry.row) && entry.row > 0)
  );
  replaceLinearComparisonPlan(
    state.linearComparisonPlan,
    data.linearComparisonPlan || createDefaultLinearComparisonPlan()
  );
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
  state.form.legend = normalizeLegendPosition(state.form.legend, state.mode.value === 'linear' ? 'bottom' : 'left');
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
  if (state.mode.value === 'linear') {
    state.adv.plot_title_position = normalizeLinearPlotTitlePosition(state.adv.plot_title_position);
  } else {
    state.adv.plot_title_position = normalizeCircularPlotTitlePosition(state.adv.plot_title_position);
  }
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
  state.adv.linear_show_accession = state.adv.linear_show_accession !== false;
  state.adv.linear_show_length = state.adv.linear_show_length !== false;
  state.adv.linear_definition_line_styles = normalizeDefinitionLineStyleState(
    state.adv.linear_definition_line_styles
  );
  state.adv.pairwise_match_style = normalizePairwiseMatchStyle(state.adv.pairwise_match_style);
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
    state.losat.blastp.maxHits = normalizePositiveInteger(state.losat.blastp?.maxHits, 5);
    state.losat.blastp.candidateLimit = null;
    if (
      (state.losat.blastp.orthogroupMemberMaxHits === null ||
        state.losat.blastp.orthogroupMemberMaxHits === undefined) &&
      state.losat.blastp.orthogroupMaxHits !== null &&
      state.losat.blastp.orthogroupMaxHits !== undefined
    ) {
      state.losat.blastp.orthogroupMemberMaxHits = state.losat.blastp.orthogroupMaxHits;
    }
    state.losat.blastp.orthogroupMembershipMode = normalizeOrthogroupMembershipMode(state.losat.blastp?.orthogroupMembershipMode);
    state.losat.blastp.orthogroupMemberMaxHits = normalizePositiveInteger(state.losat.blastp?.orthogroupMemberMaxHits, 5);
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

const serializeFile = async (file) => {
  if (!file) return null;
  const bytes = await readFileBytes(file);
  return {
    name: file.name || 'file',
    type: file.type || '',
    size: bytes.byteLength,
    lastModified: file.lastModified ?? Date.now(),
    encoding: 'base64',
    data: bytesToBase64(bytes)
  };
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
  if (
    !entry ||
    !Object.prototype.hasOwnProperty.call(entry, 'data') ||
    entry.data === null ||
    entry.data === undefined
  ) return null;
  if (isEncodedDepthFileEntry(entry)) {
    const text = decodeDepthText(entry.data);
    return new File([text], entry.name || 'depth.tsv', {
      type: entry.type || 'text/tab-separated-values',
      lastModified: entry.lastModified ?? Date.now()
    });
  }
  if (typeof entry.data !== 'string') return null;
  const bytes = base64ToBytes(entry.data);
  return new File([bytes], entry.name || 'file', {
    type: entry.type || 'application/octet-stream',
    lastModified: entry.lastModified ?? Date.now()
  });
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
  if (!cacheMap || cacheMap.size === 0) return [];
  const info = Array.isArray(state.losatCacheInfo.value) ? state.losatCacheInfo.value : [];
  const entries = [];
  const seen = new Set();

  const buildEntry = (key, cached, infoEntry = {}) => ({
    ...cloneJsonData(cached),
    key: String(key),
    filename: String(infoEntry.filename || ''),
    display: Boolean(infoEntry.display),
    ...losatCacheInfoIdentity(infoEntry)
  });

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

  return entries;
};

const applyLosatCache = (entries, legacyEnvelope = null) => {
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
      const restored = cloneJsonData(entry);
      delete restored.key;
      delete restored.filename;
      delete restored.display;
      [...LOSAT_CACHE_INFO_STRING_FIELDS, ...LOSAT_CACHE_INFO_INTEGER_FIELDS]
        .forEach((field) => delete restored[field]);
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

const applyLosatDerivedCache = (entries, legacyEvidence = null) => {
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
        payload: cloneJsonData(entry.payload)
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

const applyProteinIdentityManifest = (manifest) => {
  state.proteinIdentityManifest.value = validateProteinIdentityManifest(manifest)
    ? cloneJsonData(manifest)
    : emptyProteinIdentityManifest();
};

export const applyOrthogroupStateData = (orthogroupState = {}) => {
  const groups = Array.isArray(orthogroupState.groups) ? orthogroupState.groups : [];
  const groupIds = groups
    .map((group) => String(group?.id || '').trim())
    .filter(Boolean);
  const groupIdSet = new Set(groupIds);
  const index = buildOrthogroupFeatureIndex(groups);

  state.orthogroups.value = groups;
  state.featureOrthogroupIndex.value = index;
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

export const materializeLinearRecordFiles = (
  sequences,
  catalog,
  { layoutEnabled = false } = {}
) => {
  const sourceSequences = Array.isArray(sequences) ? sequences : [];
  if (catalog == null) return sourceSequences;
  if (catalog?.mode !== 'linear' || catalog?.status !== 'ready') {
    const issue = Array.isArray(catalog?.issues) ? catalog.issues[0] : '';
    throw new Error(issue || 'Linear record discovery is not ready.');
  }
  const records = Array.isArray(catalog.records) ? catalog.records : [];
  if (records.length === 0) {
    throw new Error('Linear record discovery did not find any records.');
  }
  const recordCountBySource = new Map();
  records.forEach((record) => {
    const sourceIndex = Number(record?.sourceIndex);
    recordCountBySource.set(sourceIndex, (recordCountBySource.get(sourceIndex) || 0) + 1);
  });
  sourceSequences.forEach((source, sourceIndex) => {
    const count = recordCountBySource.get(sourceIndex) || 0;
    if (count === 0) {
      throw new Error(`Sequence #${sourceIndex + 1}: no records were found.`);
    }
    if (count <= 1) return;
    if (layoutEnabled) {
      throw new Error(
        'Arrange in rows requires one selected record per input file. ' +
        'Turn it off or choose a Record for each multi-record file.'
      );
    }
    const hasRegion = [source.region_start, source.region_end].some(
      (value) => value !== null && value !== undefined && value !== ''
    );
    if (hasRegion) {
      throw new Error(
        `Sequence #${sourceIndex + 1}: choose a Record before setting a region on a multi-record file.`
      );
    }
  });
  const emittedBySource = new Map();
  return records.map((record) => {
    const sourceIndex = Number(record?.sourceIndex);
    const localIndex = Number(record?.localIndex);
    if (
      !Number.isInteger(sourceIndex) || sourceIndex < 0 || !sourceSequences[sourceIndex] ||
      !Number.isInteger(localIndex) || localIndex < 0
    ) {
      throw new Error('Linear record discovery returned an invalid source mapping.');
    }
    const source = sourceSequences[sourceIndex];
    const sourceUid = String(source.uid || `linear-${sourceIndex + 1}`);
    const expanded = recordCountBySource.get(sourceIndex) > 1;
    const occurrence = emittedBySource.get(sourceIndex) || 0;
    emittedBySource.set(sourceIndex, occurrence + 1);
    return {
      ...source,
      uid: expanded ? `${sourceUid}::record-${localIndex + 1}` : sourceUid,
      depth: Array.isArray(source.depth) ? [...source.depth] : source.depth,
      definition: expanded && occurrence > 0 ? '' : source.definition,
      record_subtitle: expanded && occurrence > 0 ? '' : source.record_subtitle,
      region_record_id: `#${localIndex + 1}`
    };
  });
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

const deserializeCanonicalComparisons = (comparisons) => (
  Array.isArray(comparisons)
    ? comparisons.map((comparison) => (
        isResourceBackedCanonicalComparison(comparison)
          ? mapResourceBackedCanonicalComparison(comparison, deserializeFile)
          : cloneJsonData(comparison)
      ))
    : []
);

export const adoptCanonicalRenderArtifacts = (canonical) => {
  const projection = projectCanonicalSessionRequest(canonical);
  const nextCommittedCanonicalSession = cloneCanonicalSession(canonical);
  if (projection.mode === 'linear') {
    const nextComparisons = deserializeCanonicalComparisons(
      projection.files.linearCanonicalComparisons
    );
    state.files.linearCanonicalComparisons = nextComparisons;
    committedCanonicalSession = nextCommittedCanonicalSession;
    return;
  }
  if (projection.files.c_conservation_blasts_source !== 'losat-cache') {
    committedCanonicalSession = nextCommittedCanonicalSession;
    return;
  }

  const nextBlastFiles = (projection.files.c_conservation_blasts || [])
    .map((entry) => deserializeFile(entry))
    .filter(Boolean);
  const nextFastaFiles = (projection.files.c_conservation_fastas || [])
    .map((entry) => deserializeFile(entry));
  const nextSequenceSources = (projection.files.c_conservation_sequence_sources || [])
    .map((entry) => deserializeFile(entry));
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

  state.files.c_conservation_blasts = nextBlastFiles;
  state.files.c_conservation_blasts_source = 'losat-cache';
  state.files.c_conservation_fastas = nextFastaFiles;
  state.files.c_conservation_sequence_sources = nextSequenceSources;
  if (projectedConservation) {
    state.circularConservation.enabled = true;
    state.circularConservation.source = 'losat';
    state.circularConservation.reference = projectedConservation.reference;
    state.circularConservation.labels = nextSeries.map((entry) => entry.label).join(',');
    state.circularConservation.ring_width = projectedConservation.ring_width;
    state.circularConservation.ring_gap = projectedConservation.ring_gap;
    state.circularConservation.series.splice(
      0,
      state.circularConservation.series.length,
      ...nextSeries
    );
  }
  committedCanonicalSession = nextCommittedCanonicalSession;
};

const applyFiles = (filesData) => {
  state.matchSequenceRegistry?.reset?.();
  state.files.c_gb = null;
  state.files.c_gff = null;
  state.files.c_fasta = null;
  state.files.c_depth = null;
  state.files.c_conservation_blasts = [];
  state.files.c_conservation_blasts_source = null;
  state.files.c_conservation_fastas = [];
  state.files.c_conservation_sequence_sources = [];
  state.files.linearCanonicalComparisons = [];
  state.files.d_color = null;
  state.files.t_color = null;
  state.files.blacklist = null;
  state.files.whitelist = null;
  state.files.qualifier_priority = null;
  state.linearReorderNotice.value = '';

  if (!filesData) {
    state.linearSeqs.splice(0, state.linearSeqs.length, ...normalizeLinearSeqList([]));
    state.linearRecordRows.splice(0);
    replaceLinearComparisonPlan(
      state.linearComparisonPlan,
      reconcileLinearComparisonPlan(state.linearComparisonPlan, state.linearSeqs)
    );
    return { collapsedLinearSeqs: false };
  }

  state.files.c_gb = deserializeFile(filesData.c_gb);
  state.files.c_gff = deserializeFile(filesData.c_gff);
  state.files.c_fasta = deserializeFile(filesData.c_fasta);
  state.files.c_depth = deserializeFile(filesData.c_depth);
  state.files.c_conservation_blasts = Array.isArray(filesData.c_conservation_blasts)
    ? filesData.c_conservation_blasts.map((entry) => deserializeFile(entry)).filter(Boolean)
    : [];
  state.files.c_conservation_blasts_source = filesData.c_conservation_blasts_source === 'losat-cache'
    ? 'losat-cache'
    : null;
  state.files.c_conservation_fastas = Array.isArray(filesData.c_conservation_fastas)
    ? filesData.c_conservation_fastas.map((entry) => deserializeFile(entry))
    : [];
  state.files.c_conservation_sequence_sources = Array.isArray(filesData.c_conservation_sequence_sources)
    ? filesData.c_conservation_sequence_sources.map((entry) => deserializeFile(entry))
    : [];
  state.files.linearCanonicalComparisons = deserializeCanonicalComparisons(
    filesData.linearCanonicalComparisons
  );
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
      depth: deserializeFile(seq.depth),
      losat_gencode: seq.losat_gencode ?? 1,
      definition: seq.definition ?? '',
      record_subtitle: seq.record_subtitle ?? '',
      region_record_id: seq.region_record_id ?? '',
      region_start: seq.region_start ?? null,
      region_end: seq.region_end ?? null,
      region_reverse: !!seq.region_reverse
    }));
    const normalized = normalizeLinearSeqList(loadedLinearSeqs);
    const collapsed = collapseEmptyLinearSeqList(loadedLinearSeqs);
    const collapsedLinearSeqs = collapsed.length !== normalized.length;
    state.linearSeqs.splice(0, state.linearSeqs.length, ...collapsed);
    const rowByUid = new Map(state.linearRecordRows.map((entry) => [String(entry?.uid || ''), entry]));
    state.linearRecordRows.splice(
      0,
      state.linearRecordRows.length,
      ...state.linearSeqs.map((seq, index) => ({
        uid: seq.uid,
        row: Number.isInteger(Number(rowByUid.get(seq.uid)?.row)) && Number(rowByUid.get(seq.uid)?.row) > 0
          ? Number(rowByUid.get(seq.uid).row)
          : index + 1
      }))
    );
    const comparisonFiles = new Map(
      (Array.isArray(filesData.linearComparisons) ? filesData.linearComparisons : [])
        .map((comparison) => [
          String(comparison?.id || ''),
          deserializeFile(comparison?.file)
        ])
        .filter(([id]) => id)
    );
    const planWithFiles = normalizeLinearComparisonPlan(state.linearComparisonPlan);
    planWithFiles.edges.forEach((edge) => {
      edge.file = comparisonFiles.get(edge.id) || null;
    });
    replaceLinearComparisonPlan(
      state.linearComparisonPlan,
      reconcileLinearComparisonPlan(planWithFiles, state.linearSeqs)
    );
    return { collapsedLinearSeqs };
  }

  state.linearSeqs.splice(0, state.linearSeqs.length, ...normalizeLinearSeqList([]));
  replaceLinearComparisonPlan(
    state.linearComparisonPlan,
    reconcileLinearComparisonPlan(state.linearComparisonPlan, state.linearSeqs)
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
  labelLayoutDirtyReason: state.labelLayoutDirtyReason.value,
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
  replacePlainObject(
    state.featureVisibilitySelectorCache,
    cloneJsonData(snapshot.featureVisibilitySelectorCache)
  );
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
  state.labelLayoutDirtyReason.value = snapshot.labelLayoutDirtyReason;
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
  results: serializeResults(),
  features: buildFeatureStateData(),
  editorState: buildEditorStateData(),
  orthogroupState: buildOrthogroupStateData(),
  runState: buildRunStateData(),
  losatCache: new Map(state.losatCache.value),
  losatDerivedCache: new Map(state.losatDerivedCache.value),
  proteinIdentityManifest: cloneJsonData(state.proteinIdentityManifest.value),
  legacyProteinRawCandidates: cloneJsonData(state.legacyProteinRawCandidates.value),
  legacyProteinDerivedEvidence: cloneJsonData(state.legacyProteinDerivedEvidence.value),
  losatCacheInfo: cloneJsonData(state.losatCacheInfo.value),
  committedCanonicalSession: cloneCanonicalSession(committedCanonicalSession),
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
    await nextTick();
    applyConfigData(snapshot.config);
    applyUiStateData(snapshot.ui);
    restoreLiveFileState(snapshot.files);
    state.losatCache.value = new Map(snapshot.losatCache);
    state.losatDerivedCache.value = new Map(snapshot.losatDerivedCache);
    state.proteinIdentityManifest.value = cloneJsonData(snapshot.proteinIdentityManifest);
    state.legacyProteinRawCandidates.value = cloneJsonData(snapshot.legacyProteinRawCandidates);
    state.legacyProteinDerivedEvidence.value = cloneJsonData(snapshot.legacyProteinDerivedEvidence);
    state.losatCacheInfo.value = cloneJsonData(snapshot.losatCacheInfo);
    committedCanonicalSession = cloneCanonicalSession(snapshot.committedCanonicalSession);
    state.skipCaptureBaseConfig.value = true;
    state.skipPositionReapply.value = true;
    applyResultsData(snapshot.results, snapshot.ui);
    applyFeatureStateData(snapshot.features);
    applyOrthogroupStateData(snapshot.orthogroupState);
    applyEditorStateData(snapshot.editorState);
    applyRunStateData(snapshot.runState);
    state.errorLog.value = snapshot.errorLog;
    state.resultPanelTab.value = snapshot.resultPanelTab;
    await nextTick();
    restoreSessionImportTransientState(snapshot.transients);
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
  resetSettingsState(state);
  resetLayoutState(state);
  resetRightDrawerState(state);
  state.mode.value = 'circular';
  state.cInputType.value = 'gb';
  state.lInputType.value = 'gb';
  state.sessionTitle.value = '';
  state.errorLog.value = null;
  state.results.value = [];
  state.selectedResultIndex.value = 0;
  state.resultPanelTab.value = 'preview';
  state.lastRunInfo.value = null;
  applyFiles(null);
  state.losatCache.value = new Map();
  state.losatDerivedCache.value = new Map();
  state.proteinIdentityManifest.value = emptyProteinIdentityManifest();
  state.legacyProteinRawCandidates.value = { schema: 1, entries: [] };
  state.legacyProteinDerivedEvidence.value = { schema: 1, entries: [] };
  state.losatCacheInfo.value = [];
  state.orthogroups.value = [];
  state.featureOrthogroupIndex.value = new Map();
  state.selectedOrthogroupId.value = '';
  state.selectedOrthogroupAlignmentFeature.value = '';
  clearObject(state.orthogroupNameOverrides);
  clearObject(state.orthogroupDescriptionOverrides);
  state.extractedFeatures.value = [];
  if (state.biologicalFeatures) state.biologicalFeatures.value = [];
  state.featureSelectorSafetyScope.value = [];
  state.featureRecordIds.value = [];
  state.selectedFeatureRecordIdx.value = 0;
  clearObject(state.featureColorOverrides);
  state.featureVisibilityManualRules.splice(0);
  clearObject(state.featureVisibilityOverrides);
  clearObject(state.featureVisibilitySelectorCache);
  clearObject(state.featureStrokeOverrides);
  clearObject(state.labelTextFeatureOverrides);
  state.canonicalLabelOverrideRows.value = [];
  clearObject(state.labelTextBulkOverrides);
  clearObject(state.labelTextFeatureOverrideSources);
  clearObject(state.labelVisibilityOverrides);
  state.labelOverrideContextKey.value = '';
  state.labelOverrideBuildWarning.value = '';
  state.labelLayoutDirtyReason.value = '';
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
  if (Array.isArray(resultsData)) {
    const logicalResults = normalizeLogicalResults(resultsData.map((res, idx) => (
      isCommittedSvgResult(res)
        ? { ...res, name: res?.name || `Result ${idx + 1}` }
        : {
            name: res?.name || `Result ${idx + 1}`,
            content: res?.content || ''
          }
    )));
    state.results.value = ingestSvgResults(logicalResults);
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
  labelVisibilityOverrides: cloneJsonData(state.labelVisibilityOverrides),
  labelOverrideContextKey: String(state.labelOverrideContextKey.value || '')
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
  state.labelOverrideContextKey.value = String(features.labelOverrideContextKey || '');
};

export const buildOrthogroupStateData = () => ({
  groups: Array.isArray(state.orthogroups.value) ? cloneJsonData(state.orthogroups.value) : [],
  selectedOrthogroupId: String(state.selectedOrthogroupId.value || ''),
  selectedOrthogroupAlignmentFeature: String(state.selectedOrthogroupAlignmentFeature.value || ''),
  orthogroupNameOverrides: cloneStringMap(state.orthogroupNameOverrides),
  orthogroupDescriptionOverrides: cloneStringMap(state.orthogroupDescriptionOverrides)
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

const buildSessionFeatureRecoverySnapshot = () => ({
  mode: state.mode.value,
  cInputType: state.cInputType.value,
  lInputType: state.lInputType.value,
  files: state.files,
  linearSeqs: state.linearSeqs,
  results: state.results.value,
  selectedResultIndex: state.selectedResultIndex.value,
  featureState: buildFeatureStateData(),
  editorState: buildEditorStateData(),
  orthogroupIndex: state.featureOrthogroupIndex.value
});

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

const recoverSessionFeatureMetadataIfNeeded = async ({ generationId = 'session-feature-recovery' } = {}) => {
  const validation = classifyFeatureMetadataState({
    results: state.results.value,
    selectedResultIndex: state.selectedResultIndex.value,
    extractedFeatures: state.extractedFeatures.value
  });

  if (validation.state === 'missing' || validation.state === 'alignable' || validation.state === 'stale') {
    state.featureExtractionPending.value = true;
    state.featureExtractionError.value = null;
    setFeatureEditorStatusData({
      status: 'pending-summary',
      generationId,
      error: null,
      summaryCount: state.extractedFeatures.value.length
    });
  }

  const featureVisibilityTsv = serializeFeatureVisibilityRules(state.featureVisibilityRules?.value || []);
  let plan;
  try {
    plan = await buildSessionFeatureRecoveryPlan({
      snapshot: buildSessionFeatureRecoverySnapshot(),
      featureVisibilityTsv
    });
  } catch (error) {
    console.warn('Session feature metadata recovery failed.', error);
    plan = {
      status: 'failed',
      reason: 'recovery-plan-failed',
      validation,
      warning: 'Feature metadata recovery failed. The SVG preview and pairwise popups remain available.',
      errors: [error]
    };
  }
  applySessionFeatureRecoveryPlan(plan, { generationId });
  return plan;
};

export const exportSession = async (
  titleOverride = null,
  { linearRecordCatalog = null } = {}
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
    if (!proceed) return { status: 'canceled' };
  }

  const logicalResults = serializeResults();
  const editorState = buildEditorStateData();
  if (logicalResults.length > 0) {
    if (!editorState.featureCatalog) {
      throw new Error(SESSION_FEATURE_CATALOG_SAVE_ERROR);
    }
    try {
      editorState.featureCatalog = validateFeatureCatalog(
        editorState.featureCatalog,
        logicalResults
      );
    } catch (error) {
      console.warn('Session feature catalog validation failed.', error);
      throw new Error(SESSION_FEATURE_CATALOG_SAVE_ERROR);
    }
  } else {
    editorState.featureCatalog = null;
  }

  const losatEntries = serializeLosatCache();
  const lastRunInvocation = state.lastRunInfo.value?.invocation;
  const exportableCliInvocation = isCliInvocationSessionExportable(lastRunInvocation)
    ? cloneJsonData(lastRunInvocation)
    : undefined;
  const storedConfig = buildConfigData();
  Object.assign(
    storedConfig.adv,
    normalizedArrowGeometryState(storedConfig.adv)
  );
  let committed = cloneCanonicalSession(committedCanonicalSession);
  if (committed) {
    try {
      const projected = projectCanonicalSessionRequest(committed);
      validateCurrentWriterActiveDraft({
        mode: projected.mode,
        projectedConfig: projected.config,
        storedConfig
      });
    } catch (error) {
      console.warn('Session active draft validation failed.', error);
      throw new Error(SESSION_ACTIVE_DRAFT_SAVE_ERROR);
    }
  }
  if (!committed) {
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
  const canonical = await assembleSessionResources(state, committed);
  const legacyRawCandidates = serializableLegacyProteinCandidateEnvelope(
    state.legacyProteinRawCandidates.value
  );
  const legacyDerivedEvidence = normalizeLegacyDerivedEvidence(
    state.legacyProteinDerivedEvidence.value
  );
  const sessionData = {
    format: 'gbdraw-session',
    version: SESSION_VERSION,
    createdAt: new Date().toISOString(),
    title: resolvedTitle || undefined,
    config: storedConfig,
    ui: {
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
      paletteInstantPreviewEnabled: Boolean(state.paletteInstantPreviewEnabled.value),
      appliedPaletteName: state.appliedPaletteName.value,
      appliedPaletteColors: cloneColors(state.appliedPaletteColors.value),
      pendingPaletteName: state.pendingPaletteName.value,
      pendingPaletteColors: cloneColors(state.pendingPaletteColors.value)
    },
    renderRequest: canonical.renderRequest,
    resources: canonical.resources,
    webFiles: canonical.webFiles,
    results: logicalResults,
    features: {
      selectedFeatureRecordIdx: state.selectedFeatureRecordIdx.value,
      featureColorOverrides: cloneJsonData(state.featureColorOverrides),
      featureVisibilityManualRules: normalizeFeatureVisibilityRulesForSession(state.featureVisibilityManualRules),
      featureVisibilityOverrides: normalizeFeatureVisibilityOverridesForSession(state.featureVisibilityOverrides),
      labelTextFeatureOverrides: cloneJsonData(state.labelTextFeatureOverrides),
      labelOverrideRows: cloneJsonData(state.canonicalLabelOverrideRows.value),
      labelTextBulkOverrides: cloneJsonData(state.labelTextBulkOverrides),
      labelTextFeatureOverrideSources: cloneJsonData(state.labelTextFeatureOverrideSources),
      labelVisibilityOverrides: cloneJsonData(state.labelVisibilityOverrides),
      labelOverrideContextKey: String(state.labelOverrideContextKey.value || '')
    },
    editorState,
    orthogroupState: {
      selectedOrthogroupId: String(state.selectedOrthogroupId.value || ''),
      selectedOrthogroupAlignmentFeature: String(state.selectedOrthogroupAlignmentFeature.value || ''),
      orthogroupNameOverrides: cloneStringMap(state.orthogroupNameOverrides),
      orthogroupDescriptionOverrides: cloneStringMap(state.orthogroupDescriptionOverrides)
    },
    losatCache: {
      entries: losatEntries
    },
    losatDerivedCache: {
      entries: []
    },
    proteinIdentityManifest: validateProteinIdentityManifest(state.proteinIdentityManifest.value)
      ? cloneJsonData(state.proteinIdentityManifest.value)
      : emptyProteinIdentityManifest(),
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
    console.error('Session writer validation failed.', error);
    throw new Error('Save Session could not validate the session data.');
  }

  const compressed = await compressSessionData(sessionData);
  if (!confirmLargeSessionBlob(compressed)) {
    return { status: 'canceled', compressedSize: compressed.size };
  }
  downloadBlob(compressed, sessionFilename);
  lastSessionFilename = sessionFilename;
  return { status: 'saved', blob: compressed, filename: sessionFilename };
};

export const importSession = async (e, options = {}) => {
  const file = e.target.files[0];
  if (!file) return { status: 'skipped' };

  const semanticFileWatchersSuppressedBeforeImport = Boolean(
    state.semanticFileWatchersSuppressed.value
  );
  const rollbackStateExtension = options?.rollbackState;
  let rollbackSnapshot = null;
  let rollbackExtensionSnapshot;
  let commitStarted = false;

  try {
    const text = await readSessionText(file);
    let data = JSON.parse(text, (key, value) => {
      if (key === '__proto__' || key === 'constructor' || key === 'prototype') {
        return undefined;
      }
      return value;
    });
    if (isLegacyConfigPayload(data)) {
      applyLegacyConfigPayload(data);
      alert('Legacy configuration loaded. Save as a session to use the current format.');
      return { status: 'legacy' };
    }

    const preflight = preflightSessionImport(data);
    data = preflight.data;
    const {
      sourceSessionVersion,
      canonicalProjection,
      restoredConfig,
      projectionResult
    } = preflight;
    const canonicalSession = Boolean(projectionResult);
    rollbackSnapshot = captureSessionImportSnapshot();
    if (typeof rollbackStateExtension?.capture === 'function') {
      rollbackExtensionSnapshot = rollbackStateExtension.capture();
    }
    commitStarted = true;
    state.semanticFileWatchersSuppressed.value = true;
    resetSessionBaseline();
    const ui = canonicalSession
      ? {
          ...projectionResult.editorMetadata.ui,
          ...projectionResult.artifactState.ui,
          mode: projectionResult.renderState.mode
        }
      : (data.ui || {});
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
    } else {
      if (ui.cInputType) state.cInputType.value = ui.cInputType;
      if (ui.lInputType) state.lInputType.value = ui.lInputType;
    }
    if (ui.downloadDpi) state.downloadDpi.value = ui.downloadDpi;
    // The mode watcher clears feature/editor state. Let that reset finish before
    // restoring session-owned metadata such as extractedFeatures.
    await nextTick();
    state.autoLabelReflowEnabled.value = Boolean(ui.autoLabelReflow);
    state.paletteInstantPreviewEnabled.value = Boolean(ui.paletteInstantPreviewEnabled);
    state.labelOverrideBuildWarning.value = '';
    state.labelLayoutDirtyReason.value = '';
    if (ui.featurePanelTab === 'labels' || ui.featurePanelTab === 'colors') {
      state.featurePanelTab.value = ui.featurePanelTab;
    } else {
      state.featurePanelTab.value = 'colors';
    }
    state.generatedMode.value = ui.mode === 'linear' ? 'linear' : 'circular';
    if (ui.generatedLegendPosition) {
      state.generatedLegendPosition.value = normalizeLegendPosition(
        ui.generatedLegendPosition,
        ui.mode === 'linear' ? 'bottom' : 'left'
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
      applyConfigData(restoredConfig);
    }
    restorePaletteStateFromSession(ui);
    restoreLayoutPreferences(ui, { preserveActive: Boolean(canonicalSession) });

    const { collapsedLinearSeqs } = applyFiles(
      canonicalSession ? projectionResult.restoredFiles : data.files
    );
    reconcileDepthTrackStateAfterSessionFiles();
    if (canonicalSession) {
      applyLosatCache(
        projectionResult.artifactState.losatCache?.entries,
        projectionResult.artifactState.legacyArtifacts?.proteinRawCandidates
      );
      applyLosatDerivedCache(
        projectionResult.artifactState.losatDerivedCache?.entries,
        projectionResult.artifactState.legacyArtifacts?.proteinDerivedEvidence
      );
      applyProteinIdentityManifest(projectionResult.artifactState.proteinIdentityManifest);
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
    let validatedSessionCatalog = null;
    let restoredEditorState = storedEditorState;
    if (sourceSessionVersion === SESSION_VERSION) {
      if (storedEditorState?.featureCatalog) {
        try {
          validatedSessionCatalog = validateFeatureCatalog(
            storedEditorState.featureCatalog,
            logicalImportedResults
          );
          currentCatalogFeatureState = featureStateFromCatalog(
            validatedSessionCatalog,
            { mode: state.mode.value }
          );
        } catch (catalogError) {
          console.warn('Session feature catalog is unavailable or incompatible.', catalogError);
        }
      }
      restoredEditorState = {
        ...(isPlainObject(storedEditorState) ? storedEditorState : {}),
        featureCatalog: validatedSessionCatalog
      };
    }

    const artifactFeatureState = canonicalSession
      ? cloneJsonData(projectionResult.artifactState.features)
      : {};
    if (sourceSessionVersion === SESSION_VERSION) {
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
    applyFeatureStateData(features);
    if (sourceSessionVersion === SESSION_VERSION && currentCatalogFeatureState) {
      synchronizeRestoredFeatureSummaryStatus({ generationId: 'session-load' });
    }
    const catalogSequenceSources = sourceSessionVersion === SESSION_VERSION
      ? (currentCatalogFeatureState?.sequenceSources || [])
      : [];
    const comparisonSourceAvailability = state.mode.value === 'circular'
      ? resolveCircularComparisonSequenceAvailability({
          files: state.files,
          circularConservation: state.circularConservation
        })
      : undefined;
    const catalogSequenceSourceCoverage = (
      sourceSessionVersion === SESSION_VERSION
      && validatedSessionCatalog
    )
      ? analyzeCatalogSequenceSourceCoverage({
          mode: state.mode.value,
          catalogFeatureState: validatedSessionCatalog,
          renderRequest: data.renderRequest,
          comparisonSourceAvailability
        })
      : null;
    const missingCatalogSequenceSources = sourceSessionVersion !== SESSION_VERSION
      || !catalogSequenceSourceCoverage?.complete;
    let restoredFileSequenceSources = [];
    if (missingCatalogSequenceSources) {
      try {
        restoredFileSequenceSources = await buildRestoredMatchSequenceSources({
          mode: state.mode.value,
          cInputType: state.cInputType.value,
          lInputType: state.lInputType.value,
          files: state.files,
          linearSeqs: state.linearSeqs,
          circularConservation: state.circularConservation
        });
      } catch (sequenceError) {
        console.warn('Session loaded, but match sequence recovery failed.', sequenceError);
      }
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
        }
    );
    applyEditorStateData(restoredEditorState);

    const restoredFeatureState = currentCatalogFeatureState || features || {};
    const transformRestoredSessionSvg = (svg, { applyStrokes = true } = {}) => {
      const legendGroupsChanged = normalizeLegacyLegendEntryGroups(svg);
      let compositionChanged = false;
      if (
        svg.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) === null
        && svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE) === null
      ) {
        normalizeLegacyComposition(svg, {
          legendSide: state.form.legend || 'none',
          titleSide: state.adv.plot_title_position || 'none',
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

    const committedImportedResults = validatedSessionCatalog
      ? prepareCandidateRenderCommit({
          results: logicalImportedResults,
          catalog: validatedSessionCatalog,
          mode: state.mode.value,
          featureColorOverrides: state.featureColorOverrides,
          featureStrokeOverrides: state.featureStrokeOverrides,
          legendStrokeOverrides: state.legendStrokeOverrides,
          manualSpecificRules: state.manualSpecificRules,
          legacyFeatures: features?.extractedFeatures || [],
          preparedFeatureState: currentCatalogFeatureState,
          transformSvg: (svg) => transformRestoredSessionSvg(svg, { applyStrokes: false })
        }).results
      : ingestSvgResults(logicalImportedResults, {
          transformSvg: transformRestoredSessionSvg
        });

    await nextTick();
    state.skipCaptureBaseConfig.value = true;
    state.skipPositionReapply.value = true;
    applyResultsData(committedImportedResults, ui);
    await nextTick();

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

    if (sourceSessionVersion !== SESSION_VERSION) {
      try {
        await recoverSessionFeatureMetadataIfNeeded({ generationId: 'session-load' });
      } catch (recoveryError) {
        console.warn('Session loaded, but feature metadata recovery failed.', recoveryError);
      }
    }

    state.semanticFileWatchersSuppressed.value =
      semanticFileWatchersSuppressedBeforeImport;
    await nextTick();
    committedCanonicalSession = cloneCanonicalSession(data);
    alert('Session loaded successfully!');
    return { status: 'ok', data };
  } catch (err) {
    console.error(err);
    if (commitStarted && rollbackSnapshot) {
      try {
        await restoreSessionImportSnapshot(rollbackSnapshot);
        if (typeof rollbackStateExtension?.restore === 'function') {
          await rollbackStateExtension.restore(rollbackExtensionSnapshot);
        }
      } catch (rollbackError) {
        console.error('Failed to roll back the interrupted session import.', rollbackError);
      }
    }
    const message = err?.message || 'Invalid JSON structure.';
    alert(`Failed to load session: ${message}`);
    return { status: 'error', error: err };
  } finally {
    state.semanticFileWatchersSuppressed.value =
      semanticFileWatchersSuppressedBeforeImport;
    e.target.value = '';
  }
};
