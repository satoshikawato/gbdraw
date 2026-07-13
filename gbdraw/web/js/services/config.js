import { state, normalizeLinearSeqList, collapseEmptyLinearSeqList } from '../state.js';
import { resolveColorToHex } from '../app/color-utils.js';
import { resetLayoutState, resetSettings as resetSettingsState } from './reset.js';
import { serializeCleanSvg } from './svg-serialization.js';
import { cloneJsonData } from './history-snapshot.js';
import {
  applyCircularTrackOrderPlacements,
  clampCircularTrackAxisIndex,
  inferLegacyAxisIndexFromFeature,
  normalizeCircularTrackSlots
} from '../app/circular-track-slots.js';
import {
  applyLinearTrackOrderPlacements,
  clampLinearTrackAxisIndex,
  normalizeLinearTrackSlots,
  resolveLinearTrackAxisIndex
} from '../app/linear-track-slots.js';
import {
  depthFileSlotsFromValue,
  dropInvalidManagedDepthSlots,
  reconcileDepthTracksToFiles,
  syncDepthSlotLabels,
  uploadedDepthFileCount
} from '../app/depth-track-state.js';
import {
  DEPTH_FILE_ENCODING,
  decodeDepthText,
  encodeDepthText,
  isEncodedDepthFileEntry
} from './depth-file-codec.js';
import {
  normalizeCollinearAnchorMode,
  normalizeCollinearSearchScope,
  normalizeGroupMetadataScope,
  normalizeOrthogroupMembershipMode
} from '../app/losat-normalization.js';
import { normalizeDefinitionLineStyleState } from '../app/cli-args.js';
import { isCliInvocationSessionExportable } from '../app/run-info.js';
import {
  normalizeCircularPlotTitlePosition,
  normalizeLinearPlotTitlePosition
} from '../app/plot-title-position.js';
import {
  serializeFeatureVisibilityRules,
  normalizeFeatureVisibilityRule,
  normalizeVisibilityMode,
  splitLegacyVisibilityRules
} from '../app/feature-visibility.js';
import {
  buildSessionFeatureRecoveryPlan,
  classifyFeatureMetadataState
} from '../app/session-feature-metadata.js';

const { nextTick } = window.Vue;

const SESSION_VERSION = 30;
const SUPPORTED_SESSION_VERSIONS = new Set([27, 28, 29, SESSION_VERSION]);
const LOSAT_CACHE_SCHEMA = 2;
const LOSAT_DERIVED_CACHE_SCHEMA = 1;
const LOSAT_DERIVED_CACHE_LIMIT = 16;
const CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 4;
const LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 3;
const LINEAR_TRACK_SLOT_SCHEMA_VERSION = 1;
const OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS = [
  'gapAfter',
  'gap_after',
  'innerRadius',
  'inner_radius',
  'outerRadius',
  'outer_radius',
  'placement',
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

const isRawLosatCacheEntry = (entry) =>
  Boolean(entry) &&
  entry.schema === LOSAT_CACHE_SCHEMA &&
  entry.kind === 'raw-losat' &&
  typeof entry.text === 'string';

const isPlainObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

const isLosatDerivedCacheEntry = (entry) =>
  Boolean(entry) &&
  entry.schema === LOSAT_DERIVED_CACHE_SCHEMA &&
  entry.kind === 'derived-losatp-payload' &&
  typeof entry.key === 'string' &&
  isPlainObject(entry.payload);

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

const downloadJson = (data, filename, { pretty = true } = {}) => {
  const blob = new Blob([pretty ? JSON.stringify(data, null, 2) : JSON.stringify(data)], {
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

const validateImportedCircularTrackSlots = (configData = {}) => {
  const adv = configData && typeof configData === 'object' ? configData.adv : null;
  if (!adv || typeof adv !== 'object' || Array.isArray(adv)) return;
  if (!Object.prototype.hasOwnProperty.call(adv, 'circular_track_slots')) return;

  if (
    adv.circular_track_slots_schema_version !== CIRCULAR_TRACK_SLOT_SCHEMA_VERSION &&
    adv.circular_track_slots_schema_version !== LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION
  ) {
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
};

const validateImportedLinearTrackSlots = (configData = {}) => {
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
};

const hasStoredLayoutValue = (value) => typeof value === 'string' && value.trim() !== '';

const normalizeFeatureShape = (value) => (String(value || '').trim().toLowerCase() === 'arrow' ? 'arrow' : 'rectangle');

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

const normalizeStrokeWidth = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric >= 0 ? numeric : null;
};

const cloneJsonArray = (value) => {
  if (!Array.isArray(value)) return [];
  try {
    return JSON.parse(JSON.stringify(value));
  } catch (_err) {
    return [];
  }
};

const cloneJsonObject = (value) => {
  if (!isPlainObject(value)) return {};
  try {
    return JSON.parse(JSON.stringify(value));
  } catch (_err) {
    return {};
  }
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
      source.large_tick_interval ?? source.tick_interval ?? (index === 0 ? legacyAdv.depth_tick_interval : null)
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

const cloneLosatForConfig = () => {
  const cloned = JSON.parse(JSON.stringify(state.losat || {}));
  if (cloned.blastp && typeof cloned.blastp === 'object' && !Array.isArray(cloned.blastp)) {
    delete cloned.blastp.collinearAnchorMode;
  }
  return cloned;
};

export const buildConfigData = () => ({
  form: state.form,
  adv: state.adv,
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
  blastSource: state.blastSource.value,
  losatProgram: state.losatProgram.value,
  circularConservation: state.circularConservation,
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
  }
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
  }
});

const normalizeEditorStateData = (editorState = {}) => {
  const defaults = defaultEditorStateData();
  const source = isPlainObject(editorState) ? editorState : {};
  const legend = isPlainObject(source.legend) ? source.legend : {};
  const featureStrokes = isPlainObject(source.featureStrokes) ? source.featureStrokes : {};
  const originalSvgStroke = isPlainObject(source.originalSvgStroke) ? source.originalSvgStroke : {};

  return {
    legend: {
      entries: cloneJsonArray(legend.entries),
      deletedEntries: cloneJsonArray(legend.deletedEntries),
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
    }
  };
};

const replacePlainObject = (target, source) => {
  Object.keys(target).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = value;
  });
};

export const applyEditorStateData = (editorState = {}) => {
  const normalized = normalizeEditorStateData(editorState);

  state.legendEntries.value = normalized.legend.entries;
  state.deletedLegendEntries.value = normalized.legend.deletedEntries;
  state.originalLegendOrder.value = normalized.legend.originalOrder;
  state.originalLegendColors.value = normalized.legend.originalColors;
  replacePlainObject(state.legendColorOverrides, normalized.legend.colorOverrides);
  replacePlainObject(state.legendStrokeOverrides, normalized.legend.strokeOverrides);
  state.addedLegendCaptions.value = new Set(normalized.legend.addedCaptions);
  replacePlainObject(state.featureStrokeOverrides, normalized.featureStrokes.overrides);
  state.originalSvgStroke.value = normalized.originalSvgStroke;
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

  return {
    ...data,
    version: SESSION_VERSION,
    editorState: normalizeEditorStateData(data.editorState)
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
  state.suppressCircularMultiRecordDefaults.value = shouldSuppressCircularMultiRecordDefaults(data.form);
  validateImportedCircularTrackSlots(data);
  validateImportedLinearTrackSlots(data);
  applyConfigData(data);
  restorePaletteStateAfterConfigImport();
};

const shouldSuppressCircularMultiRecordDefaults = (incomingForm) => {
  if (state.mode.value !== 'circular') return false;
  if (!incomingForm || typeof incomingForm !== 'object' || Array.isArray(incomingForm)) return false;
  if (!Object.prototype.hasOwnProperty.call(incomingForm, 'multi_record_canvas')) return false;
  return state.form.multi_record_canvas === false && incomingForm.multi_record_canvas === true;
};

const syncActiveCircularLayoutCache = () => {
  const normalizedLegend = normalizeLegendPosition(state.form.legend, 'left');
  const normalizedPlotTitlePosition = normalizeCircularPlotTitlePosition(state.adv.plot_title_position);

  state.circularLegendPosition.value = normalizedLegend;
  state.circularPlotTitlePosition.value = normalizedPlotTitlePosition;

  if (state.form.multi_record_canvas) {
    state.circularMultiRecordLegendPosition.value = normalizedLegend;
    state.circularMultiRecordPlotTitlePosition.value = normalizedPlotTitlePosition;
  } else {
    state.circularSingleRecordLegendPosition.value = normalizedLegend;
    state.circularSingleRecordPlotTitlePosition.value = normalizedPlotTitlePosition;
  }

  return {
    legend: normalizedLegend,
    plotTitlePosition: normalizedPlotTitlePosition
  };
};

const getStoredCircularLayout = (useMultiRecord) => {
  const singleLegend = normalizeLegendPosition(state.circularSingleRecordLegendPosition.value, 'left');
  const singlePlotTitlePosition = normalizeCircularPlotTitlePosition(state.circularSingleRecordPlotTitlePosition.value);

  if (useMultiRecord) {
    return {
      legend: hasStoredLayoutValue(state.circularMultiRecordLegendPosition.value)
        ? normalizeLegendPosition(state.circularMultiRecordLegendPosition.value, singleLegend)
        : singleLegend,
      plotTitlePosition: hasStoredLayoutValue(state.circularMultiRecordPlotTitlePosition.value)
        ? normalizeCircularPlotTitlePosition(state.circularMultiRecordPlotTitlePosition.value)
        : singlePlotTitlePosition
    };
  }

  return {
    legend: singleLegend,
    plotTitlePosition: singlePlotTitlePosition
  };
};

const resolveActiveCircularLayout = () => getStoredCircularLayout(Boolean(state.form.multi_record_canvas));

const syncActiveModePlotTitlePosition = () => {
  if (state.mode.value === 'linear') {
    state.form.legend = normalizeLegendPosition(state.form.legend, 'bottom');
    state.linearLegendPosition.value = state.form.legend;
    state.linearPlotTitlePosition.value = normalizeLinearPlotTitlePosition(state.adv.plot_title_position);
    state.adv.plot_title_position = state.linearPlotTitlePosition.value;
    return;
  }

  state.form.legend = normalizeLegendPosition(state.form.legend, 'left');
  state.adv.plot_title_position = normalizeCircularPlotTitlePosition(state.adv.plot_title_position);
  syncActiveCircularLayoutCache();
};

const restoreSessionPlotTitlePositions = (ui = {}) => {
  const activeMode = state.mode.value === 'linear' ? 'linear' : 'circular';
  const activePosition = normalizeLinearPlotTitlePosition(state.adv.plot_title_position);
  const hasLinearPlotTitlePosition =
    typeof ui.linearPlotTitlePosition === 'string' && ui.linearPlotTitlePosition.trim() !== '';

  state.linearPlotTitlePosition.value = hasLinearPlotTitlePosition
    ? normalizeLinearPlotTitlePosition(ui.linearPlotTitlePosition)
    : activeMode === 'linear'
      ? activePosition
      : 'bottom';
};

const restoreSessionCircularLayoutCaches = (ui = {}) => {
  const activeLegend = normalizeLegendPosition(state.form.legend, 'left');
  const activePlotTitlePosition = normalizeCircularPlotTitlePosition(state.adv.plot_title_position);
  const legacyLegend = hasStoredLayoutValue(ui.circularLegendPosition)
    ? normalizeLegendPosition(ui.circularLegendPosition, activeLegend)
    : hasStoredLayoutValue(ui.legend)
      ? normalizeLegendPosition(ui.legend, activeLegend)
    : activeLegend;
  const legacyPlotTitlePosition = hasStoredLayoutValue(ui.circularPlotTitlePosition)
    ? normalizeCircularPlotTitlePosition(ui.circularPlotTitlePosition)
    : activePlotTitlePosition;

  state.circularSingleRecordLegendPosition.value = hasStoredLayoutValue(ui.circularSingleRecordLegendPosition)
    ? normalizeLegendPosition(ui.circularSingleRecordLegendPosition, legacyLegend)
    : (state.form.multi_record_canvas ? legacyLegend : activeLegend);
  state.circularSingleRecordPlotTitlePosition.value = hasStoredLayoutValue(ui.circularSingleRecordPlotTitlePosition)
    ? normalizeCircularPlotTitlePosition(ui.circularSingleRecordPlotTitlePosition)
    : (state.form.multi_record_canvas ? legacyPlotTitlePosition : activePlotTitlePosition);
  state.circularMultiRecordLegendPosition.value = hasStoredLayoutValue(ui.circularMultiRecordLegendPosition)
    ? normalizeLegendPosition(ui.circularMultiRecordLegendPosition, legacyLegend)
    : (state.form.multi_record_canvas ? activeLegend : legacyLegend);
  state.circularMultiRecordPlotTitlePosition.value = hasStoredLayoutValue(ui.circularMultiRecordPlotTitlePosition)
    ? normalizeCircularPlotTitlePosition(ui.circularMultiRecordPlotTitlePosition)
    : (state.form.multi_record_canvas ? activePlotTitlePosition : legacyPlotTitlePosition);

  const activeLayout = resolveActiveCircularLayout();
  state.circularLegendPosition.value = activeLayout.legend;
  state.circularPlotTitlePosition.value = activeLayout.plotTitlePosition;
};

export const applyConfigData = (data) => {
  if (data.form) safeDeepMerge(state.form, data.form);
  if (data.adv) safeDeepMerge(state.adv, data.adv);
  state.adv.rich_feature_popup = data?.adv?.rich_feature_popup !== false;
  if (state.adv.label_placement === 'on_feature') {
    state.adv.label_placement = 'above_feature';
  }
  state.adv.label_rendering = normalizeLabelRendering(state.adv.label_rendering);
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
  if (state.form.linear_track_layout === 'spreadout') {
    state.form.linear_track_layout = 'above';
  } else if (state.form.linear_track_layout === 'tuckin') {
    state.form.linear_track_layout = 'below';
  } else if (!['above', 'middle', 'below'].includes(state.form.linear_track_layout)) {
    state.form.linear_track_layout = 'middle';
  }
  state.form.plot_title = String(state.form.plot_title || '');
  state.form.legend = normalizeLegendPosition(state.form.legend, state.mode.value === 'linear' ? 'bottom' : 'left');
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
  state.adv.depth_tick_interval = normalizePositiveNumberOrNull(state.adv.depth_tick_interval);
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
      const maxGap = Number(state.losat.blastp?.collinearMaxGeneGap);
      state.losat.blastp.collinearMaxGeneGap = Number.isInteger(maxGap) && maxGap >= 0 ? maxGap : 0;
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
    state.currentColors.value = state.normalizePaletteColors(normalizeColorMap(data.colors));
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
  if (data.blastSource) state.blastSource.value = String(data.blastSource);
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

const serializeFileArray = async (files) => {
  const items = Array.isArray(files) ? files.filter(Boolean) : [];
  return Promise.all(items.map((file) => serializeFile(file)));
};

const serializeDepthFile = async (file) => {
  if (Array.isArray(file)) {
    return Promise.all(file.map((item) => serializeDepthFile(item)));
  }
  if (!file) return null;
  try {
    const text = await file.text();
    const encodedDepth = encodeDepthText(text);
    if (encodedDepth) {
      return {
        name: file.name || 'depth.tsv',
        type: file.type || 'text/tab-separated-values',
        size: file.size || text.length,
        lastModified: file.lastModified || Date.now(),
        encoding: DEPTH_FILE_ENCODING,
        data: encodedDepth
      };
    }
  } catch (err) {
    console.warn('Failed to encode depth file for session storage; falling back to base64.', err);
  }
  return serializeFile(file);
};

const fileSizeOf = (fileOrFiles) => (
  Array.isArray(fileOrFiles)
    ? fileOrFiles.reduce((sum, file) => sum + (file?.size || 0), 0)
    : (fileOrFiles?.size || 0)
);

const deserializeFile = (entry) => {
  if (Array.isArray(entry)) {
    return entry.map((item) => deserializeFile(item));
  }
  if (!entry || !entry.data) return null;
  if (isEncodedDepthFileEntry(entry)) {
    const text = decodeDepthText(entry.data);
    return new File([text], entry.name || 'depth.tsv', {
      type: entry.type || 'text/tab-separated-values',
      lastModified: entry.lastModified || Date.now()
    });
  }
  if (typeof entry.data !== 'string') return null;
  const bytes = base64ToUint8(entry.data);
  return new File([bytes], entry.name || 'file', {
    type: entry.type || 'application/octet-stream',
    lastModified: entry.lastModified || Date.now()
  });
};

let activePreviewRuntime = null;

export const setPreviewRuntime = (runtime) => {
  activePreviewRuntime = runtime || null;
};

export const serializeResults = () => {
  if (activePreviewRuntime?.flushActiveResult) {
    activePreviewRuntime.flushActiveResult();
    return state.results.value.map((res, idx) => ({
      name: res.name || `Result ${idx + 1}`,
      content: res.content
    }));
  }

  const currentSvg = (() => {
    if (!state.svgContainer.value) return null;
    const svg = state.svgContainer.value.querySelector('svg');
    if (!svg) return null;
    return serializeCleanSvg(svg);
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

  const buildEntry = (key, cached, { filename = '', display = false } = {}) => {
    const entry = {
      schema: LOSAT_CACHE_SCHEMA,
      kind: 'raw-losat',
      key,
      filename,
      display,
      text: cached.text,
      program: cached.program || '',
      queryCanonicalHash: cached.queryCanonicalHash || '',
      subjectCanonicalHash: cached.subjectCanonicalHash || ''
    };
    if (cached.flow) entry.flow = cached.flow;
    if (cached.outfmt) entry.outfmt = cached.outfmt;
    if (Array.isArray(cached.args)) entry.args = cached.args.map((arg) => String(arg));
    return entry;
  };

  info.forEach((entry, idx) => {
    if (!entry || !entry.key) return;
    const cached = cacheMap.get(entry.key);
    if (!isRawLosatCacheEntry(cached)) return;
    entries.push(buildEntry(entry.key, cached, {
      filename: entry.filename || `losat_pair_${idx + 1}.tsv`,
      display: entry.display !== false
    }));
    seen.add(entry.key);
  });

  cacheMap.forEach((value, key) => {
    if (seen.has(key)) return;
    if (!isRawLosatCacheEntry(value)) return;
    entries.push(buildEntry(key, value));
  });

  return entries;
};

const applyLosatCache = (entries) => {
  const map = new Map();
  const info = [];

  if (Array.isArray(entries)) {
    entries.forEach((entry, idx) => {
      if (
        !entry ||
        entry.schema !== LOSAT_CACHE_SCHEMA ||
        entry.kind !== 'raw-losat' ||
        !entry.key ||
        typeof entry.text !== 'string'
      ) {
        return;
      }
      const restored = {
        schema: LOSAT_CACHE_SCHEMA,
        kind: 'raw-losat',
        text: entry.text,
        program: entry.program || '',
        queryCanonicalHash: entry.queryCanonicalHash || '',
        subjectCanonicalHash: entry.subjectCanonicalHash || ''
      };
      if (entry.flow) restored.flow = entry.flow;
      if (entry.outfmt) restored.outfmt = entry.outfmt;
      if (Array.isArray(entry.args)) restored.args = entry.args.map((arg) => String(arg));
      map.set(entry.key, restored);
      if (entry.display === false) return;
      info.push({
        key: entry.key,
        filename: entry.filename || `losat_pair_${idx + 1}.tsv`,
        display: true
      });
    });
  }

  state.losatCache.value = map;
  state.losatCacheInfo.value = info;
};

const pruneLosatDerivedCache = (map) => {
  if (!map || typeof map.delete !== 'function') return;
  while (map.size > LOSAT_DERIVED_CACHE_LIMIT) {
    const oldestKey = map.keys().next().value;
    if (oldestKey === undefined) break;
    map.delete(oldestKey);
  }
};

const serializeLosatDerivedCache = () => {
  const cacheMap = state.losatDerivedCache?.value;
  if (!cacheMap || cacheMap.size === 0) return [];
  const entries = [];

  cacheMap.forEach((value, key) => {
    const entry = {
      schema: LOSAT_DERIVED_CACHE_SCHEMA,
      kind: 'derived-losatp-payload',
      key: String(key || value?.key || ''),
      mode: String(value?.mode || ''),
      payload: cloneJsonData(value?.payload)
    };
    if (!isLosatDerivedCacheEntry(entry)) return;
    entries.push(entry);
  });

  return entries.slice(-LOSAT_DERIVED_CACHE_LIMIT);
};

const applyLosatDerivedCache = (entries) => {
  const map = new Map();

  if (Array.isArray(entries)) {
    entries.forEach((entry) => {
      if (!isLosatDerivedCacheEntry(entry)) return;
      map.set(entry.key, {
        schema: LOSAT_DERIVED_CACHE_SCHEMA,
        kind: 'derived-losatp-payload',
        key: entry.key,
        mode: String(entry.mode || ''),
        payload: cloneJsonData(entry.payload)
      });
    });
  }

  pruneLosatDerivedCache(map);
  state.losatDerivedCache.value = map;
};

const buildOrthogroupIndexKey = (recordIndex, svgId) => `${Number(recordIndex)}:${String(svgId || '').trim()}`;

const enrichExtractedFeaturesWithOrthogroups = (index) => {
  if (!(index instanceof Map) || !Array.isArray(state.extractedFeatures.value) || state.extractedFeatures.value.length === 0) {
    return;
  }
  state.extractedFeatures.value = state.extractedFeatures.value.map((feature) => {
    const ids = [
      feature?.svg_id,
      feature?.stable_svg_id,
      feature?.stableFeatureSvgId,
      feature?.stable_feature_id
    ].map((value) => String(value || '').trim()).filter(Boolean);
    const uniqueIds = Array.from(new Set(ids));
    if (uniqueIds.length === 0) return feature;

    const recordIndexes = [
      feature?.fileIdx,
      feature?.recordIndex,
      feature?.record_index,
      feature?.record_idx
    ].map((value) => Number(value)).filter((value) => Number.isInteger(value));
    let entry = null;
    for (const recordIndex of recordIndexes) {
      entry = uniqueIds
        .map((id) => index.get(buildOrthogroupIndexKey(recordIndex, id)))
        .find(Boolean);
      if (entry) break;
    }
    entry = entry || uniqueIds.map((id) => index.get(id)).find(Boolean);
    if (!entry) return feature;
    return {
      ...feature,
      proteinId: entry.proteinId,
      sourceProteinId: entry.sourceProteinId,
      orthogroupId: entry.orthogroupId,
      orthogroupMemberCount: entry.orthogroupMemberCount,
      orthogroupRecordCoverage: entry.orthogroupRecordCoverage,
      orthogroupRepresentative: entry.orthogroupRepresentative,
      orthogroupScope: entry.orthogroupScope,
      orthogroupSourceRecordIndex: entry.orthogroupSourceRecordIndex,
      orthogroupMember: entry.orthogroupMember
    };
  });
};

export const applyOrthogroupStateData = (orthogroupState = {}) => {
  const groups = Array.isArray(orthogroupState.groups) ? orthogroupState.groups : [];
  const groupIds = groups
    .map((group) => String(group?.id || '').trim())
    .filter(Boolean);
  const groupIdSet = new Set(groupIds);
  const index = new Map();

  groups.forEach((group) => {
    const orthogroupId = String(group?.id || '').trim();
    const members = Array.isArray(group?.members) ? group.members : [];
    const memberCount = Number(group?.member_count || members.length || 0);
    const recordCoverage = Number(group?.record_coverage_count || new Set(
      members.map((member) => Number(member?.recordIndex)).filter((recordIndex) => Number.isInteger(recordIndex))
    ).size || 0);
    const orthogroupScope = normalizeGroupMetadataScope(group?.scope);
    const sourceRecordIndex = Number(group?.source_record_index);
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
        orthogroupScope,
        orthogroupSourceRecordIndex: Number.isInteger(sourceRecordIndex) ? sourceRecordIndex : null,
        orthogroupMember: member
      };
      index.set(buildOrthogroupIndexKey(recordIndex, featureSvgId), entry);
      if (!index.has(featureSvgId)) index.set(featureSvgId, entry);
    });
  });

  state.orthogroups.value = groups;
  state.featureOrthogroupIndex.value = index;
  enrichExtractedFeaturesWithOrthogroups(index);
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

export const serializeFiles = async () => {
  const normalizedLinearSeqs = normalizeLinearSeqList(state.linearSeqs);
  const linearSeqs = await Promise.all(
    normalizedLinearSeqs.map(async (seq) => ({
      uid: seq.uid,
      gb: await serializeFile(seq.gb),
      gff: await serializeFile(seq.gff),
      fasta: await serializeFile(seq.fasta),
      depth: await serializeDepthFile(seq.depth),
      blast: await serializeFile(seq.blast),
      losat_gencode: seq.losat_gencode ?? 1,
      losat_filename: seq.losat_filename ?? '',
      definition: seq.definition ?? '',
      record_subtitle: seq.record_subtitle ?? '',
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
    c_depth: await serializeDepthFile(state.files.c_depth),
    c_conservation_blasts: await serializeFileArray(state.files.c_conservation_blasts),
    c_conservation_fastas: await serializeFileArray(state.files.c_conservation_fastas),
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
  state.files.c_depth = null;
  state.files.c_conservation_blasts = [];
  state.files.c_conservation_fastas = [];
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
  state.files.c_depth = deserializeFile(filesData.c_depth);
  state.files.c_conservation_blasts = Array.isArray(filesData.c_conservation_blasts)
    ? filesData.c_conservation_blasts.map((entry) => deserializeFile(entry)).filter(Boolean)
    : [];
  state.files.c_conservation_fastas = Array.isArray(filesData.c_conservation_fastas)
    ? filesData.c_conservation_fastas.map((entry) => deserializeFile(entry)).filter(Boolean)
    : [];
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
      blast: deserializeFile(seq.blast),
      losat_gencode: seq.losat_gencode ?? 1,
      losat_filename: seq.losat_filename ?? '',
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
    return { collapsedLinearSeqs };
  }

  state.linearSeqs.splice(0, state.linearSeqs.length, ...normalizeLinearSeqList([]));
  return { collapsedLinearSeqs: false };
};

const representativeLinearDepthFiles = () => {
  const rows = state.linearSeqs.map((seq) => depthFileSlotsFromValue(seq.depth));
  const maxDepthTracks = Math.max(...rows.map((row) => row.length), 0);
  return Array.from({ length: maxDepthTracks }, (_, trackIndex) => (
    rows.map((row) => row[trackIndex]).find(Boolean) || null
  ));
};

const reconcileDepthTrackStateAfterSessionFiles = () => {
  const circularDepthCount = uploadedDepthFileCount(state.files.c_depth);
  const linearDepthCount = state.linearSeqs.reduce(
    (maxCount, seq) => Math.max(maxCount, uploadedDepthFileCount(seq.depth)),
    0
  );
  const targetDepthTrackCount = Math.max(1, circularDepthCount, linearDepthCount);
  const depthFilesForTracks = circularDepthCount > 0
    ? depthFileSlotsFromValue(state.files.c_depth).filter(Boolean)
    : representativeLinearDepthFiles();
  const normalizedTracks = reconcileDepthTracksToFiles({
    files: depthFilesForTracks,
    depthTracks: state.adv.depth_tracks,
    targetCount: targetDepthTrackCount,
    defaults: {
      depthColor: state.adv.depth_color,
      depthHeight: state.adv.depth_height,
      largeTickInterval: state.adv.depth_tick_interval,
      smallTickInterval: state.adv.depth_small_tick_interval,
      tickFontSize: state.adv.depth_tick_font_size
    }
  });
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

const clearObject = (target) => {
  Object.keys(target).forEach((key) => {
    delete target[key];
  });
};

const resetSessionBaseline = () => {
  activePreviewRuntime?.clearActiveRuntime?.();
  preservedCliOptions = null;
  resetSettingsState(state);
  resetLayoutState(state);
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
  state.losatCacheInfo.value = [];
  state.orthogroups.value = [];
  state.featureOrthogroupIndex.value = new Map();
  state.selectedOrthogroupId.value = '';
  state.selectedOrthogroupAlignmentFeature.value = '';
  clearObject(state.orthogroupNameOverrides);
  clearObject(state.orthogroupDescriptionOverrides);
  state.extractedFeatures.value = [];
  state.featureSelectorSafetyScope.value = [];
  state.featureRecordIds.value = [];
  state.selectedFeatureRecordIdx.value = 0;
  clearObject(state.featureColorOverrides);
  state.featureVisibilityManualRules.splice(0);
  clearObject(state.featureVisibilityOverrides);
  clearObject(state.featureVisibilitySelectorCache);
  clearObject(state.featureStrokeOverrides);
  clearObject(state.labelTextFeatureOverrides);
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
  state.showFeaturePanel.value = false;
  state.showLegendPanel.value = false;
};

export const buildUiStateData = ({ includePreviewNavigation = true } = {}) => {
  const currentLegend = state.form.legend;
  const isLinear = state.mode.value === 'linear';
  const currentPlotTitlePosition = state.adv.plot_title_position;
  const activeCircularLegend = isLinear
    ? normalizeLegendPosition(state.circularLegendPosition.value, 'left')
    : normalizeLegendPosition(currentLegend, 'left');
  const activeCircularPlotTitlePosition = isLinear
    ? normalizeCircularPlotTitlePosition(state.circularPlotTitlePosition.value)
    : normalizeCircularPlotTitlePosition(currentPlotTitlePosition);
  const savedCircularSingleRecordLegendPosition =
    !isLinear && !state.form.multi_record_canvas
      ? activeCircularLegend
      : normalizeLegendPosition(state.circularSingleRecordLegendPosition.value, activeCircularLegend);
  const savedCircularSingleRecordPlotTitlePosition =
    !isLinear && !state.form.multi_record_canvas
      ? activeCircularPlotTitlePosition
      : normalizeCircularPlotTitlePosition(state.circularSingleRecordPlotTitlePosition.value);
  const savedCircularMultiRecordLegendPosition =
    !isLinear && state.form.multi_record_canvas
      ? activeCircularLegend
      : hasStoredLayoutValue(state.circularMultiRecordLegendPosition.value)
        ? normalizeLegendPosition(state.circularMultiRecordLegendPosition.value, savedCircularSingleRecordLegendPosition)
        : savedCircularSingleRecordLegendPosition;
  const savedCircularMultiRecordPlotTitlePosition =
    !isLinear && state.form.multi_record_canvas
      ? activeCircularPlotTitlePosition
      : hasStoredLayoutValue(state.circularMultiRecordPlotTitlePosition.value)
        ? normalizeCircularPlotTitlePosition(state.circularMultiRecordPlotTitlePosition.value)
        : savedCircularSingleRecordPlotTitlePosition;

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
    legend: currentLegend,
    circularLegendPosition: activeCircularLegend,
    linearLegendPosition: isLinear
      ? normalizeLegendPosition(currentLegend, 'bottom')
      : normalizeLegendPosition(state.linearLegendPosition.value, 'bottom'),
    circularPlotTitlePosition: activeCircularPlotTitlePosition,
    linearPlotTitlePosition: isLinear
      ? normalizeLinearPlotTitlePosition(currentPlotTitlePosition)
      : normalizeLinearPlotTitlePosition(state.linearPlotTitlePosition.value),
    circularSingleRecordLegendPosition: savedCircularSingleRecordLegendPosition,
    circularSingleRecordPlotTitlePosition: savedCircularSingleRecordPlotTitlePosition,
    circularMultiRecordLegendPosition: savedCircularMultiRecordLegendPosition,
    circularMultiRecordPlotTitlePosition: savedCircularMultiRecordPlotTitlePosition,
    featurePanelTab: state.featurePanelTab.value,
    cInputType: state.cInputType.value,
    lInputType: state.lInputType.value,
    blastSource: state.blastSource.value,
    losatProgram: state.losatProgram.value,
    downloadDpi: state.downloadDpi.value,
    autoLabelReflow: Boolean(state.autoLabelReflowEnabled.value),
    paletteInstantPreviewEnabled: Boolean(state.paletteInstantPreviewEnabled.value),
    appliedPaletteName: state.appliedPaletteName.value,
    appliedPaletteColors: cloneColors(state.appliedPaletteColors.value),
    pendingPaletteName: state.pendingPaletteName.value,
    pendingPaletteColors: cloneColors(state.pendingPaletteColors.value),
    circularBaseConfig: cloneJsonData(state.circularBaseConfig.value),
    linearBaseConfig: {
      ...cloneJsonData(state.linearBaseConfig.value),
      diagramBaseTransforms: []
    },
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
  if (ui.blastSource) state.blastSource.value = String(ui.blastSource);
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
  restoreSessionCircularLayoutCaches(ui);
  restoreSessionPlotTitlePositions(ui);

  if (ui.circularBaseConfig && typeof ui.circularBaseConfig === 'object') {
    state.circularBaseConfig.value = cloneJsonData(ui.circularBaseConfig);
  }
  if (ui.linearBaseConfig && typeof ui.linearBaseConfig === 'object') {
    state.linearBaseConfig.value = {
      ...cloneJsonData(ui.linearBaseConfig),
      diagramBaseTransforms: new Map()
    };
  }
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

  if (state.mode.value === 'linear') {
    const nextLinearLegend = hasStoredLayoutValue(ui.linearLegendPosition)
      ? normalizeLegendPosition(ui.linearLegendPosition, 'bottom')
      : hasStoredLayoutValue(ui.legend)
        ? normalizeLegendPosition(ui.legend, 'bottom')
        : normalizeLegendPosition(state.form.legend, 'bottom');
    state.form.legend = nextLinearLegend;
    state.linearLegendPosition.value = nextLinearLegend;
    state.adv.plot_title_position = state.linearPlotTitlePosition.value;
  } else {
    const activeCircularLayout = resolveActiveCircularLayout();
    state.form.legend = activeCircularLayout.legend;
    state.adv.plot_title_position = activeCircularLayout.plotTitlePosition;
    state.circularLegendPosition.value = activeCircularLayout.legend;
    state.circularPlotTitlePosition.value = activeCircularLayout.plotTitlePosition;
  }
};

export const applyResultsData = (resultsData = [], ui = {}) => {
  if (Array.isArray(resultsData)) {
    state.results.value = resultsData.map((res, idx) => ({
      name: res.name || `Result ${idx + 1}`,
      content: res.content || ''
    }));
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
  featureSelectorSafetyScope: cloneJsonData(state.featureSelectorSafetyScope.value),
  featureRecordIds: cloneJsonData(state.featureRecordIds.value),
  selectedFeatureRecordIdx: state.selectedFeatureRecordIdx.value,
  featureColorOverrides: cloneJsonData(state.featureColorOverrides),
  featureVisibilityManualRules: normalizeFeatureVisibilityRulesForSession(state.featureVisibilityManualRules),
  featureVisibilityOverrides: normalizeFeatureVisibilityOverridesForSession(state.featureVisibilityOverrides),
  labelTextFeatureOverrides: cloneJsonData(state.labelTextFeatureOverrides),
  labelTextBulkOverrides: cloneJsonData(state.labelTextBulkOverrides),
  labelTextFeatureOverrideSources: cloneJsonData(state.labelTextFeatureOverrideSources),
  labelVisibilityOverrides: cloneJsonData(state.labelVisibilityOverrides),
  labelOverrideContextKey: String(state.labelOverrideContextKey.value || '')
});

export const applyFeatureStateData = (features = {}) => {
  state.extractedFeatures.value = Array.isArray(features.extractedFeatures)
    ? features.extractedFeatures
    : [];
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

  if (state.extractedFeatures.value.length > 0) {
    state.featureExtractionError.value = null;
    setFeatureEditorStatusData({
      status: 'summary-ready',
      generationId,
      error: null,
      summaryCount: state.extractedFeatures.value.length
    });
  }
};

const recoverSessionFeatureMetadataIfNeeded = async ({ generationId = 'session-feature-recovery' } = {}) => {
  const validation = classifyFeatureMetadataState({
    results: state.results.value,
    selectedResultIndex: state.selectedResultIndex.value,
    extractedFeatures: state.extractedFeatures.value
  });

  if (validation.state === 'ready' || validation.state === 'not-needed') {
    applySessionFeatureRecoveryPlan({ status: 'ready', validation }, { generationId });
    return { status: 'ready', validation };
  }

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

const guardSessionFeatureMetadataForExport = async () => {
  try {
    return await recoverSessionFeatureMetadataIfNeeded({ generationId: 'session-export' });
  } catch (error) {
    console.warn('Session feature metadata export guard failed.', error);
    state.featureExtractionPending.value = false;
    const warning = 'Feature metadata recovery failed. The SVG preview and pairwise popups remain available.';
    state.featureExtractionError.value = { summary: warning, details: [] };
    setFeatureEditorStatusData({
      status: 'failed',
      generationId: 'session-export',
      error: warning,
      summaryCount: 0
    });
    return { status: 'failed', error };
  }
};

export const exportSession = async (titleOverride = null) => {
  const losatEntries = serializeLosatCache();
  const losatBytes = losatEntries.reduce((sum, entry) => sum + (entry.text ? entry.text.length : 0), 0);
  const currentLegend = state.form.legend;
  const isLinear = state.mode.value === 'linear';
  const currentPlotTitlePosition = state.adv.plot_title_position;
  const activeCircularLegend = isLinear
    ? normalizeLegendPosition(state.circularLegendPosition.value, 'left')
    : normalizeLegendPosition(currentLegend, 'left');
  const activeCircularPlotTitlePosition = isLinear
    ? normalizeCircularPlotTitlePosition(state.circularPlotTitlePosition.value)
    : normalizeCircularPlotTitlePosition(currentPlotTitlePosition);
  const savedCircularLegend = activeCircularLegend;
  const savedLinearLegend = isLinear
    ? normalizeLegendPosition(currentLegend, 'bottom')
    : normalizeLegendPosition(state.linearLegendPosition.value, 'bottom');
  const savedCircularPlotTitlePosition = activeCircularPlotTitlePosition;
  const savedLinearPlotTitlePosition = isLinear
    ? normalizeLinearPlotTitlePosition(currentPlotTitlePosition)
    : normalizeLinearPlotTitlePosition(state.linearPlotTitlePosition.value);
  const savedCircularSingleRecordLegendPosition =
    !isLinear && !state.form.multi_record_canvas
      ? activeCircularLegend
      : normalizeLegendPosition(state.circularSingleRecordLegendPosition.value, activeCircularLegend);
  const savedCircularSingleRecordPlotTitlePosition =
    !isLinear && !state.form.multi_record_canvas
      ? activeCircularPlotTitlePosition
      : normalizeCircularPlotTitlePosition(state.circularSingleRecordPlotTitlePosition.value);
  const savedCircularMultiRecordLegendPosition =
    !isLinear && state.form.multi_record_canvas
      ? activeCircularLegend
      : hasStoredLayoutValue(state.circularMultiRecordLegendPosition.value)
        ? normalizeLegendPosition(state.circularMultiRecordLegendPosition.value, savedCircularSingleRecordLegendPosition)
        : savedCircularSingleRecordLegendPosition;
  const savedCircularMultiRecordPlotTitlePosition =
    !isLinear && state.form.multi_record_canvas
      ? activeCircularPlotTitlePosition
      : hasStoredLayoutValue(state.circularMultiRecordPlotTitlePosition.value)
        ? normalizeCircularPlotTitlePosition(state.circularMultiRecordPlotTitlePosition.value)
        : savedCircularSingleRecordPlotTitlePosition;
  const resolvedTitle =
    typeof titleOverride === 'string'
      ? titleOverride.trim()
      : typeof state.sessionTitle?.value === 'string'
        ? state.sessionTitle.value.trim()
        : '';
  const sessionFilename = buildSessionFilename(resolvedTitle);
  await guardSessionFeatureMetadataForExport();
  const totalBytes =
    (state.files.c_gb?.size || 0) +
    (state.files.c_gff?.size || 0) +
    (state.files.c_fasta?.size || 0) +
    fileSizeOf(state.files.c_depth) +
    (Array.isArray(state.files.c_conservation_blasts)
      ? state.files.c_conservation_blasts.reduce((sum, file) => sum + (file?.size || 0), 0)
      : 0) +
    (Array.isArray(state.files.c_conservation_fastas)
      ? state.files.c_conservation_fastas.reduce((sum, file) => sum + (file?.size || 0), 0)
      : 0) +
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
        fileSizeOf(seq.depth) +
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

  const lastRunInvocation = state.lastRunInfo.value?.invocation;
  const exportableCliInvocation = isCliInvocationSessionExportable(lastRunInvocation)
    ? JSON.parse(JSON.stringify(lastRunInvocation))
    : undefined;
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
      generatedMultiRecordCanvas: Boolean(state.generatedMultiRecordCanvas.value),
      generatedCircularPlotTitlePosition: normalizeCircularPlotTitlePosition(
        state.generatedCircularPlotTitlePosition.value
      ),
      legend: currentLegend,
      circularLegendPosition: savedCircularLegend,
      linearLegendPosition: savedLinearLegend,
      circularPlotTitlePosition: savedCircularPlotTitlePosition,
      linearPlotTitlePosition: savedLinearPlotTitlePosition,
      circularSingleRecordLegendPosition: savedCircularSingleRecordLegendPosition,
      circularSingleRecordPlotTitlePosition: savedCircularSingleRecordPlotTitlePosition,
      circularMultiRecordLegendPosition: savedCircularMultiRecordLegendPosition,
      circularMultiRecordPlotTitlePosition: savedCircularMultiRecordPlotTitlePosition,
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
    files: await serializeFiles(),
    results: serializeResults(),
    features: {
      extractedFeatures: sanitizeExtractedFeaturesForSession(state.extractedFeatures.value),
      featureSelectorSafetyScope: cloneJsonData(state.featureSelectorSafetyScope.value),
      featureRecordIds: state.featureRecordIds.value,
      selectedFeatureRecordIdx: state.selectedFeatureRecordIdx.value,
      featureColorOverrides: JSON.parse(JSON.stringify(state.featureColorOverrides)),
      featureVisibilityManualRules: normalizeFeatureVisibilityRulesForSession(state.featureVisibilityManualRules),
      featureVisibilityOverrides: normalizeFeatureVisibilityOverridesForSession(state.featureVisibilityOverrides),
      labelTextFeatureOverrides: JSON.parse(JSON.stringify(state.labelTextFeatureOverrides)),
      labelTextBulkOverrides: JSON.parse(JSON.stringify(state.labelTextBulkOverrides)),
      labelTextFeatureOverrideSources: JSON.parse(JSON.stringify(state.labelTextFeatureOverrideSources)),
      labelVisibilityOverrides: JSON.parse(JSON.stringify(state.labelVisibilityOverrides)),
      labelOverrideContextKey: String(state.labelOverrideContextKey.value || '')
    },
    editorState: buildEditorStateData(),
    orthogroupState: {
      groups: Array.isArray(state.orthogroups.value) ? JSON.parse(JSON.stringify(state.orthogroups.value)) : [],
      selectedOrthogroupId: String(state.selectedOrthogroupId.value || ''),
      selectedOrthogroupAlignmentFeature: String(state.selectedOrthogroupAlignmentFeature.value || ''),
      orthogroupNameOverrides: cloneStringMap(state.orthogroupNameOverrides),
      orthogroupDescriptionOverrides: cloneStringMap(state.orthogroupDescriptionOverrides)
    },
    losatCache: {
      entries: losatEntries
    },
    losatDerivedCache: {
      entries: serializeLosatDerivedCache()
    },
    cliInvocation: exportableCliInvocation
  };

  if (lastSessionFilename && lastSessionFilename === sessionFilename) {
    const proceed = confirm(`Download "${sessionFilename}" again? Your browser may overwrite or rename the file.`);
    if (!proceed) return;
  }
  lastSessionFilename = sessionFilename;
  downloadJson(sessionData, sessionFilename, { pretty: false });
};

export const importSession = async (e, options = {}) => {
  const file = e.target.files[0];
  if (!file) return { status: 'skipped' };

  if (file.size > 200 * 1024 * 1024) {
    alert('Session file is too large.');
    return { status: 'error' };
  }

  try {
    const text = await file.text();
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

    data = normalizeSessionData(data);
    resetSessionBaseline();

    const ui = data.ui || {};
    state.sessionTitle.value = typeof data.title === 'string' ? data.title : '';
    if (ui.mode) state.mode.value = ui.mode;
    if (ui.cInputType) state.cInputType.value = ui.cInputType;
    if (ui.lInputType) state.lInputType.value = ui.lInputType;
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
    if (ui.linearLegendPosition) {
      state.linearLegendPosition.value = normalizeLegendPosition(ui.linearLegendPosition, 'bottom');
    }
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

    if (data.config) {
      hydrateMissingMultiRecordPositionsFromCliInvocation(data.config, data.cliInvocation);
      state.suppressCircularMultiRecordDefaults.value = shouldSuppressCircularMultiRecordDefaults(data.config.form);
      validateImportedCircularTrackSlots(data.config);
      validateImportedLinearTrackSlots(data.config);
      applyConfigData(data.config);
    }
    restorePaletteStateFromSession(ui);
    restoreSessionCircularLayoutCaches(ui);
    restoreSessionPlotTitlePositions(ui);

    const { collapsedLinearSeqs } = applyFiles(data.files);
    reconcileDepthTrackStateAfterSessionFiles();
    applyLosatCache(data.losatCache?.entries);
    applyLosatDerivedCache(data.losatDerivedCache?.entries);
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
    applyFeatureStateData(features);

    applyOrthogroupStateData(
      data.orthogroupState && typeof data.orthogroupState === 'object'
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
    applyEditorStateData(data.editorState);

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

    if (state.mode.value === 'linear') {
      const nextLinearLegend = hasStoredLayoutValue(ui.linearLegendPosition)
        ? normalizeLegendPosition(ui.linearLegendPosition, 'bottom')
        : hasStoredLayoutValue(ui.legend)
          ? normalizeLegendPosition(ui.legend, 'bottom')
          : normalizeLegendPosition(state.form.legend, 'bottom');
      state.form.legend = nextLinearLegend;
      state.linearLegendPosition.value = nextLinearLegend;
      state.adv.plot_title_position = state.linearPlotTitlePosition.value;
    } else {
      const activeCircularLayout = resolveActiveCircularLayout();
      state.form.legend = activeCircularLayout.legend;
      state.adv.plot_title_position = activeCircularLayout.plotTitlePosition;
      state.circularLegendPosition.value = activeCircularLayout.legend;
      state.circularPlotTitlePosition.value = activeCircularLayout.plotTitlePosition;
    }

    await recoverSessionFeatureMetadataIfNeeded({ generationId: 'session-load' });

    if (typeof options?.afterLoad === 'function') {
      try {
        await options.afterLoad({ data, ui });
      } catch (callbackError) {
        console.warn('Session loaded, but post-load refresh failed.', callbackError);
      }
    }

    alert('Session loaded successfully!');
    return { status: 'ok', data };
  } catch (err) {
    console.error(err);
    state.suppressCircularMultiRecordDefaults.value = false;
    const message = err?.message || 'Invalid JSON structure.';
    alert(`Failed to load session: ${message}`);
    return { status: 'error', error: err };
  } finally {
    e.target.value = '';
  }
};
