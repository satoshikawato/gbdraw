import { writeCanonicalRecordReverseComplement } from '../app/record-display-options.js';
import { canonicalFeaturePlacements } from './feature-placement.js';
export { canonicalFeaturePlacements } from './feature-placement.js';
import { buildDefaultColorOverrideTsv, normalizePaletteColors } from '../app/color-utils.js';
import {
  parseColorTable,
  parsePriorityRules,
  parseSpecificRules,
  parseWhitelistRules,
  serializeSpecificRules
} from '../app/file-imports.js';
import {
  buildLabelOverrideTsv,
  parseLabelOverrideTsv,
  serializeLabelOverrideRows
} from '../app/feature-editor/label-override-table.js';
import {
  parseFeatureVisibilityRules,
  serializeFeatureVisibilityRules
} from '../app/feature-visibility.js';
import {
  applyCircularGeometryShortcuts,
  buildCircularTrackSlotPayload,
  CIRCULAR_TRACK_RENDERERS,
  createDefaultCircularTrackSlots,
  hasCircularGeometryShortcuts,
  inferLegacyAxisIndexFromFeature,
  migrateLegacyCircularTrackSlot,
  migrateLegacyCircularTrackSlotSpec,
  normalizeCircularTrackSlot,
  parseCircularTrackSlotSpecs
} from '../app/circular-track-slots.js';
import { projectCircularMeasureDraft } from '../app/circular-track-slots/measure-editor.js';
import {
  buildLinearTrackSlotPayload,
  LINEAR_TRACK_RENDERERS,
  LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  migrateLinearTrackSlotsToCurrentSchema,
  parseLinearTrackSlotSpecs
} from '../app/linear-track-slots.js';
import {
  isRecordMajorDepthFileMatrix,
  normalizeRecordMajorDepthFileRows,
  parseDepthTrackIndexIdentity
} from '../app/depth-track-state.js';
import {
  buildDisambiguatedRecordEntries,
  resolveCircularRequestRecordSet,
  resolveDisambiguatedRecordSelection
} from '../app/record-options.js';
import {
  orderedConservationSources,
  orderedOptionalConservationFiles
} from '../app/conservation-series.js';
import {
  resolveLinearRecordEffectiveDefinition,
  resolveLinearRecordEffectiveSubtitle
} from '../app/linear-sources.js';
import {
  linearRecordLayoutHasSharedRow,
  resolveEffectiveLinearRecordRows
} from '../app/linear-record-layout.js';
import { resolveLinearLabelVisibility } from '../app/linear-label-visibility.js';
import {
  assertValidCustomTrackPlan,
  parseOptionalPixel,
  validateCustomTrackPlan,
  validateTrackSlotBindingInvariants
} from '../app/track-slot-validation.js';
import { annotationOptionsPayload, normalizeAnnotationSets } from '../app/annotations/state.js';
import { classifyOptionalNumber, classifyOptionalPositiveNumber, projectOptionalNumber } from '../utils/optional-positive-number.js';
import { diagnosticError } from './error-normalization.js';
import { materializeLegacySimilarityAlignment } from './legacy-similarity-alignment.js';
import {
  arrowHeadLengthRatioForState,
  defaultFeatureRendering,
  normalizeArrowHeadLengthRatio,
  normalizeArrowShaftWidthRatio,
  normalizeFeatureRenderingMap,
  visibleFeatureUnderlaysForState
} from '../utils/feature-rendering.js';
import {
  comparisonFiltersForMode,
  effectiveLinearAxisColor,
  MODE_DEFAULT_FEATURE_TYPES,
  modeProfile,
  resolveComparisonThresholds,
  trackDefaultsForMode
} from '../mode-profiles.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';
import {
  createDefaultLayoutPreferences,
  updateActiveLayoutPreference
} from '../app/layout-preferences.js';
import {
  migratePersistedCircularMultiRecordSizeMode,
  migratePersistedLinearLabelPlacement,
  migratePersistedLinearTrackLayout,
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
} from '../app/current-option-values.js';
import {
  normalizeCollinearAnchorMode,
  normalizeCollinearSearchScope,
  normalizeOrthogroupMembershipMode
} from '../app/losat-normalization.js';
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
  textToBase64,
  textToBytes
} from './file-content-cache.js';
import {
  adoptedSessionResourceDescriptor,
  adoptCurrentSessionResources,
  createCombinedSessionResourceFileView,
  createSessionResourceFileView,
  validateWebFileBindings
} from './session-resource-backing.js';
import { normalizeLinearComparisonPlan } from '../app/linear-comparisons.js';
import {
  getResourcePayloadOwner,
  setResourcePayloadOwner
} from './resource-payload-owner.js';
import { sha256Hex } from './byte-utils.js';
import { cloneJsonData } from './json-clone.js';
import { isCanonicalResourceReferenceField } from './canonical-resource-references.js';
import { recordStructuralMetric } from './runtime-test-hooks.js';
import { recordDisplayKey, requestedRecordTransform } from '../app/record-display-options.js';

export const CANONICAL_REQUEST_SCHEMA = 8;
const SUPPORTED_CANONICAL_REQUEST_SCHEMAS = new Set([
  1, 2, 5, 6, 7, CANONICAL_REQUEST_SCHEMA
]);

const requireExactCanonicalKeys = (value, keys, path) => {
  if (!value || typeof value !== 'object' || Array.isArray(value)) {
    throw new Error(`${path} must be an object.`);
  }
  const expected = [...keys].sort();
  const actual = Object.keys(value).sort();
  if (actual.length !== expected.length || actual.some((key, index) => key !== expected[index])) {
    throw new Error(`${path} contains missing or unknown fields.`);
  }
  return value;
};

const requireCanonicalText = (value, path) => {
  if (typeof value !== 'string' || !value.trim() || value.includes('\0')) {
    throw new Error(`${path} must be non-empty text without NUL.`);
  }
  return value.trim();
};

const canonicalAlignmentAnchor = (value, path) => {
  const anchor = requireExactCanonicalKeys(value, [
    'recordKey',
    'biologicalFeatureId',
    'sourceFeatureIndex',
    'stableFeatureSvgId'
  ], path);
  if (anchor.sourceFeatureIndex !== null && (
    !Number.isSafeInteger(anchor.sourceFeatureIndex) || anchor.sourceFeatureIndex < 0
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
    sourceFeatureIndex: anchor.sourceFeatureIndex,
    stableFeatureSvgId: anchor.stableFeatureSvgId === null
      ? null
      : anchor.stableFeatureSvgId.trim()
  };
};

const canonicalSimilarityAlignment = (value, recordKeys, path) => {
  if (value === null) return null;
  const plan = requireExactCanonicalKeys(value, [
    'schema', 'groupId', 'reference', 'records'
  ], path);
  if (plan.schema !== 2) throw new Error(`${path}.schema must be 2.`);
  if (!Array.isArray(plan.records) || plan.records.length === 0) {
    throw new Error(`${path}.records must be a non-empty array.`);
  }
  const reference = canonicalAlignmentAnchor(plan.reference, `${path}.reference`);
  const records = plan.records.map((raw, index) => {
    const decisionPath = `${path}.records[${index}]`;
    const decision = requireExactCanonicalKeys(raw, [
      'recordKey', 'status', 'rationale', 'anchor'
    ], decisionPath);
    const recordKey = requireCanonicalText(decision.recordKey, `${decisionPath}.recordKey`);
    const anchor = decision.anchor === null
      ? null
      : canonicalAlignmentAnchor(decision.anchor, `${decisionPath}.anchor`);
    if (anchor && anchor.recordKey !== recordKey) {
      throw new Error(`${decisionPath}.anchor belongs to another record.`);
    }
    const alignedRationales = new Set([
      'user_selected', 'only_usable_candidate', 'unique_direct_rbh'
    ]);
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
    return {
      recordKey,
      status: decision.status,
      rationale: decision.rationale,
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
  const actual = new Set(decisionKeys);
  if (expected.size !== actual.size || [...expected].some((key) => !actual.has(key))) {
    throw new Error(`${path}.records does not cover the displayed record keys.`);
  }
  return {
    schema: 2,
    groupId: requireCanonicalText(plan.groupId, `${path}.groupId`),
    reference,
    records
  };
};

const canonicalRecordTranslations = (value, recordKeys, path, { requireCoverage = false } = {}) => {
  if (!Array.isArray(value)) throw new Error(`${path} must be an array.`);
  const translations = value.map((raw, index) => {
    const itemPath = `${path}[${index}]`;
    const item = requireExactCanonicalKeys(raw, ['recordKey', 'x', 'y'], itemPath);
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

const canonicalRecordDisplay = (raw) => {
  if (!raw || Object.keys(raw).sort().join(',') !== 'isCircular,startCoordinate'
    || (raw.isCircular !== null && typeof raw.isCircular !== 'boolean')
    || (raw.startCoordinate !== null && (!Number.isSafeInteger(raw.startCoordinate) || raw.startCoordinate < 1))
    || (raw.isCircular === false && raw.startCoordinate !== null)) {
    throw new Error('Record display requires nullable isCircular and a positive integer or null startCoordinate.');
  }
  return { ...raw };
};

// Canonical schemas 1-2 omitted values that matched the former shared API
// defaults. Keep those values stable when reading sparse persisted requests.
const HISTORICAL_COMPARISON_DEFAULTS = Object.freeze({
  bitscore: 50,
  evalue: 1e-5,
  identity: 70,
  alignmentLength: 0
});
const HISTORICAL_FEATURE_TYPES = Object.freeze([
  'CDS',
  'rRNA',
  'tRNA',
  'tmRNA',
  'ncRNA',
  'misc_RNA',
  'repeat_region'
]);
const HISTORICAL_CONFIG_OVERRIDES = Object.freeze({
  circular: Object.freeze({
    'canvas.show_gc': true,
    'canvas.show_skew': true
  }),
  linear: Object.freeze({
    'canvas.show_gc': true,
    'canvas.show_skew': true,
    'objects.axis.linear.stroke_color': 'gray'
  })
});

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
  pairwiseMatchStyle: 'objects.blast_match.style'
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

const LINEAR_DEFINITION_STYLE_PATHS = Object.freeze(
  Object.fromEntries(
    ['name', 'subtitle', 'replicon', 'accession', 'length'].map((kind) => [
      kind,
      `objects.definition.linear.line_styles.${kind}`
    ])
  )
);
const LINEAR_DEFINITION_STYLE_FIELDS = Object.freeze(['font_size', 'font_weight', 'fill']);

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

export const managedConfigOverridePathsForMode = (mode) => {
  if (!['circular', 'linear'].includes(mode)) {
    throw new Error(`Unsupported config override mode: ${String(mode)}.`);
  }
  const excluded = mode === 'circular'
    ? LINEAR_ONLY_GUI_CONFIG_OVERRIDE_PATHS
    : CIRCULAR_ONLY_GUI_CONFIG_OVERRIDE_PATHS;
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

const legacyFlatConfigKey = (semanticName) => (
  semanticName
    .replace(/([A-Z]+)([A-Z][a-z])/g, '$1_$2')
    .replace(/([a-z0-9])([A-Z])/g, '$1_$2')
    .toLowerCase()
);

const safePrefix = (value, fallback = 'out') => {
  const normalized = String(value || '').trim().replace(/[\\/]+/g, '_');
  return normalized && normalized !== '.' && normalized !== '..' ? normalized : fallback;
};

const explicitOutputPrefix = (value) => {
  const raw = String(value || '').trim();
  return raw ? safePrefix(raw) : null;
};

const circularRecordId = (record, index) => (
  String(record?.record_id ?? record?.recordId ?? '').trim() || `Record_${index + 1}`
);

const resolveCircularBatchPrefixes = (records, explicitPrefix) => {
  if (explicitPrefix !== null) {
    if (records.length === 1) return [explicitPrefix];
    return records.map((_, index) => `${explicitPrefix}_${index + 1}`);
  }
  const prefixes = [];
  const used = new Set();
  records.forEach((record, index) => {
    const base = safePrefix(circularRecordId(record, index));
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

const renderOutputPayload = (prefix) => ({
  prefix,
  formats: ['svg'],
  overwrite: false,
  interactiveMetadataPolicy: 'auto'
});

const projectComparisonThresholds = ({ evalue, bitscore, identity, alignmentLength }) => ({
  evalue: Number(evalue), bitscore, identity, alignmentLength
});
// An integer setting that the current request does not use for rendering
// (a mode-inactive protein parameter or a LOSAT execution choice) keeps its
// documented fallback; projected render values use projectOptionalNumber.
const integerSettingOr = (value, fallback, minimum) => {
  const { value: numeric } = classifyOptionalNumber(value);
  return Number.isInteger(numeric) && numeric >= minimum ? numeric : fallback;
};

// A record label or subtitle is saved resolved. Sessions written since file
// defaults were introduced also record the per-record value, which says directly
// whether the record inherited the file default or overrode it with the same
// text. Older sessions carry no such value, so fall back to comparing them.
const resolveSavedRecordOverride = ({ savedOverride, fileDefault, resolved }) => {
  if (savedOverride !== undefined && savedOverride !== null) return String(savedOverride);
  if (!fileDefault) return resolved;
  return resolved === fileDefault ? '' : resolved;
};

const validateProjectedDepthSources = (depthRows, logicalTrackCount) => {
  for (let trackIndex = 0; trackIndex < logicalTrackCount; trackIndex += 1) {
    const hasSource = depthRows.some((row) => (
      Array.isArray(row) && Boolean(row[trackIndex]?.resourceId)
    ));
    if (!hasSource) {
      throw new Error(
        `Depth series #${trackIndex + 1} (logical track index ${trackIndex}) has no source in any record.`
      );
    }
  }
};

const normalizeResourceName = (resourceId, name) => {
  const basename = String(name || 'resource.dat').replace(/\\/g, '/').split('/').pop();
  const safe = basename.replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^[._]+|[._]+$/g, '');
  const prefix = `${resourceId}-`;
  let leaf = safe || 'resource.dat';
  while (leaf.startsWith(prefix)) {
    leaf = leaf.slice(prefix.length);
  }
  return `${prefix}${leaf || 'resource.dat'}`;
};

const normalizeOriginalResourceName = (name) => {
  const basename = String(name || '')
    .replace(/\\/g, '/')
    .split('/')
    .pop()
    .replace(/[\u0000-\u001f\u007f]/g, '')
    .trim();
  if (!basename || basename === '.' || basename === '..') return '';
  return basename.slice(0, 1024);
};

const generatedResourceValues = new WeakMap();
const canonicalTypedResourceBackings = new WeakMap();

// The analysis helper supplies the bytes already emitted by the Python typed codec.
export const bindCanonicalTypedResource = (value, descriptor) => {
  const expectedKind = value?.kind === "result" ? "collinearity-result" : "orthogroup-result";
  if (!value || value.schema !== 3 || !["result", "orthogroupResult"].includes(value.kind)
    || descriptor?.kind !== expectedKind || descriptor.encoding !== "base64"
    || descriptor.type !== "application/json" || typeof descriptor.data !== "string"
    || !Number.isSafeInteger(descriptor.size) || descriptor.size < 0) {
    throw new Error("Analysis canonical resource backing is invalid.");
  }
  canonicalTypedResourceBackings.set(value, descriptor);
};
const createResourceBuilder = ({ encode = true } = {}) => {
  const resources = {};
  const resourceOriginalNames = {};
  const fileResourceIds = new Map();

  const addFile = (resourceId, kind, entry) => {
    if (!entry || typeof entry !== 'object' || Array.isArray(entry)) {
      throw new Error(`Canonical resource ${resourceId} is missing.`);
    }
    if (resources[resourceId]) return resourceId;
    const adoptedSource = getSessionResourceSource(entry);
    const encodedEntry = adoptedSource?.descriptor || entry;
    const effectiveKind = String(encodedEntry.kind || kind);
    const kindResourceIds = fileResourceIds.get(effectiveKind) || new WeakMap();
    const owner = getResourcePayloadOwner(entry);
    const identity = getSessionResourceSource(owner)?.descriptor || owner;
    const existingResourceId = kindResourceIds.get(identity);
    if (existingResourceId) return existingResourceId;
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

  const addText = (resourceId, kind, name, text) => {
    if (resources[resourceId]) return resourceId;
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

  const addJson = (resourceId, kind, name, value) => {
    if (resources[resourceId]) return resourceId;
    if (!encode) {
      resources[resourceId] = { kind };
      generatedResourceValues.set(resources[resourceId], { bindings: [value] });
      return resourceId;
    }
    const backing = canonicalTypedResourceBackings.get(value);
    if (backing) {
      if (backing.kind !== kind) throw new Error("Analysis canonical resource kind does not match.");
      resources[resourceId] = { ...backing, name: normalizeResourceName(resourceId, name) };
      generatedResourceValues.set(resources[resourceId], { bindings: [value] });
      return resourceId;
    }
    const normalized = JSON.stringify(value);
    const bytes = textToBytes(normalized);
    resources[resourceId] = {
      kind,
      name: normalizeResourceName(resourceId, name),
      type: 'application/json',
      size: bytes.byteLength,
      lastModified: 0,
      encoding: 'base64',
      data: bytesToBase64(bytes)
    };
    generatedResourceValues.set(resources[resourceId], { bindings: [value] });
    return resourceId;
  };

  const addCanonicalTable = (resourceId, rows, columns = null) => {
    if (resources[resourceId]) return resourceId;
    const normalizedRows = Array.isArray(rows)
      ? rows.filter((row) => row && typeof row === 'object' && !Array.isArray(row))
      : [];
    const preferredColumns = [
      'query',
      'subject',
      'identity',
      'alignment_length',
      'mismatches',
      'gap_opens',
      'qstart',
      'qend',
      'sstart',
      'send',
      'evalue',
      'bitscore'
    ];
    const requestedColumns = Array.isArray(columns)
      ? columns.map((column) => String(column || '').trim()).filter(Boolean)
      : [];
    const discoveredColumns = new Set(
      normalizedRows.flatMap((row) => Object.keys(row))
    );
    const orderedColumns = [
      ...requestedColumns,
      ...preferredColumns.filter((column) => discoveredColumns.has(column)),
      ...Array.from(discoveredColumns)
        .filter((column) => !requestedColumns.includes(column) && !preferredColumns.includes(column))
        .sort()
    ];
    const uniqueColumns = Array.from(new Set(orderedColumns));
    if (uniqueColumns.length === 0) {
      uniqueColumns.push(...preferredColumns);
    }
    const escapeCell = (value) => {
      if (value === null || value === undefined) return '';
      const text = String(value);
      return /[\t\n\r"]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
    };
    const lines = [
      uniqueColumns.map(escapeCell).join('\t'),
      ...normalizedRows.map((row) => (
        uniqueColumns.map((column) => escapeCell(row[column])).join('\t')
      ))
    ];
    return addText(
      resourceId,
      'canonical-tsv',
      `${resourceId}.tsv`,
      `${lines.join('\n')}\n`
    );
  };

  return {
    resources,
    resourceOriginalNames,
    addFile,
    addText,
    addJson,
    addCanonicalTable
  };
};

const fileRef = (resourceId) => ({ resourceId, representation: 'file' });
const publicationFileRef = (resources, files, key, fallbackId) => {
  if (!files[key]) return null;
  const source = getSessionResourceSource(files[key]);
  return {
    resourceId: resources.addFile(source?.resourceId || fallbackId, fallbackId, files[key]),
    representation: source?.descriptor?.kind === 'canonical-tsv' ? 'canonicalTsv' : 'file'
  };
};

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

const circularPresentationPayload = (form, { hasRegion = false } = {}) => (
  presentationPayload({
    label: form.circular_record_label,
    subtitle: form.circular_record_subtitle,
    reverseComplement: !hasRegion && form.circular_reverse
  })
);

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
    selector: selectorPayload(record.value),
    start,
    end,
    reverseComplement: Boolean(form.circular_reverse)
  };
};

const circularRecordKey = (record) => {
  const preserved = String(record?.recordKey || '').trim();
  if (preserved) return preserved;
  const recordId = safePrefix(record?.recordId, 'record');
  const selector = safePrefix(record?.selector, '1');
  return `circular-${recordId}-${selector}`;
};

const linearRegionPayload = (seq) => {
  const start = projectOptionalNumber(seq?.region_start, { field: 'region' });
  const end = projectOptionalNumber(seq?.region_end, { field: 'region' });
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

const buildRecords = ({ state, filesData, resources }) => {
  if (state.mode.value === 'linear') {
    const resolvedRows = resolveEffectiveLinearRecordRows(
      filesData.linearSeqs,
      state.linearRecordRows,
      { enabled: Boolean(state.linearRecordLayoutEnabled?.value) }
    );
    const canonicalCardinalityByUid = new Map(
      (state.linearRecordRows || []).map((entry) => [entry.uid, entry.canonicalCardinality])
    );
    const records = (filesData.linearSeqs || []).map((seq, index) => {
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
            gridRow: state.linearRecordLayoutEnabled?.value
              ? (resolvedRows[index]?.row ?? null)
              : null,
            gridColumn: state.linearRecordLayoutEnabled?.value
              ? (resolvedRows[index]?.canonicalColumn ?? null)
              : null
          }),
          reverseComplement: region ? false : Boolean(seq.region_reverse)
        }
      };
    });
    return { records, circularSourceIndexes: null, circularSourceCount: null };
  }

  if (Array.isArray(filesData.circularRecords) && filesData.circularRecords.length > 0) {
    const records = filesData.circularRecords.map((record, index) => {
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
      !state.form.multi_record_canvas &&
      state.adv.circular_grouping_intent !== 'batch'
    );
    if (singleJourney) {
      const record = records[0];
      const savedSelector = canonicalRecordSelector(record);
      const requestedSelector = String(
        state.form.circular_record_selector || savedSelector || ''
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
      const region = circularRegionPayload(state.form, selected);
      records[0] = {
        ...record,
        recordKey: record.recordKey || circularRecordKey(selected),
        selector: region ? null : selectorPayload(selected.value),
        region,
        presentation: circularPresentationPayload(state.form, {
          hasRegion: Boolean(region)
        })
      };
    }
    return {
      records,
      circularSourceIndexes: records.map((_, index) => index),
      circularSourceCount: records.length
    };
  }
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
  const recordSet = resolveCircularRequestRecordSet({
    records: state.circularRecordList.value,
    selector: state.form.circular_record_selector,
    multiRecordCanvas: state.form.multi_record_canvas,
    groupingIntent: state.adv.circular_grouping_intent
  });
  if (recordSet.selectionFailure) throw diagnosticError('RECORD_SELECTION', { reason: recordSet.selectionFailure });
  const { recordSelectors } = recordSet;
  const selectedRecords = recordSet.records.length > 0 ? recordSet.records : [null];
  const singleJourney = (
    selectedRecords.length === 1 &&
    recordSet.singlePresentation
  );
  const records = selectedRecords.map((record, index) => {
    const region = singleJourney
      ? circularRegionPayload(state.form, record)
      : null;
    return {
      recordKey: singleJourney && record
        ? circularRecordKey(record)
        : `record-${index + 1}`,
      cardinality: 'exactly_one',
      source,
      selector: region
        ? null
        : selectorPayload(record?.value ?? record?.selector),
      region,
      presentation: singleJourney
        ? circularPresentationPayload(state.form, { hasRegion: Boolean(region) })
        : presentationPayload()
    };
  });
  const circularSourceIndexes = selectedRecords.map((record, index) => (
    record && Number.isInteger(record.sourceIndex) ? record.sourceIndex : index
  ));
  return {
    records,
    circularSourceIndexes,
    circularSourceCount: recordSelectors.length || records.length
  };
};

const buildConfigOverrides = (
  state,
  {
    depthRequested = Boolean(state.form.show_depth),
    hasComparisonIntent = false,
    linearHasSharedRow = false
  } = {}
) => {
  const { form, adv } = state;
  const circular = state.mode.value === 'circular';
  // R7: a config leaf is projected literally; Python validates its domain.
  const leaf = (configPath, value) => projectOptionalNumber(value, { configPath });
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
  const linearAxisManaged = state.modeProfileStateManager?.isManaged?.(
    adv,
    'axis_stroke_color'
  ) === true;
  const linearAxisStrokeColor = circular
    ? null
    : effectiveLinearAxisColor({
        axisColor: adv.axis_stroke_color,
        rulerOnAxis:
          form.show_scale !== false && Boolean(form.linear_ruler_on_axis),
        managed: linearAxisManaged
      });
  const overrides = {
    [CONFIG_OVERRIDE_PATHS.arrowHeadLengthRatio]:
      normalizeArrowHeadLengthRatio(adv.arrow_head_length_ratio),
    [CONFIG_OVERRIDE_PATHS.arrowShaftWidthRatio]:
      normalizeArrowShaftWidthRatio(adv.arrow_shaft_width_ratio),
    [CONFIG_OVERRIDE_PATHS.blockStrokeColor]: adv.block_stroke_color || null,
    [CONFIG_OVERRIDE_PATHS.lineStrokeColor]: adv.line_stroke_color || null,
    [CONFIG_OVERRIDE_PATHS.labelRendering]: adv.label_rendering || 'auto',
    [CONFIG_OVERRIDE_PATHS.showGc]: circular ? !form.suppress_gc : Boolean(form.show_gc),
    [CONFIG_OVERRIDE_PATHS.showSkew]: circular
      ? !form.suppress_skew
      : Boolean(form.show_skew),
    [CONFIG_OVERRIDE_PATHS.showDepth]: Boolean(depthRequested),
    [MODE_LABEL_SCOPE_PATHS[state.mode.value]]: circular
      ? ({ none: 'none', out: 'outer', both: 'both' }[form.labels_mode] || 'none')
      : form.show_labels_linear,
    [CONFIG_OVERRIDE_PATHS.strandedness]: Boolean(form.separate_strands),
    [CONFIG_OVERRIDE_PATHS.resolveOverlaps]: Boolean(adv.resolve_overlaps),
    [CONFIG_OVERRIDE_PATHS.featureOverlapToleranceBp]: adv.feature_overlap_tolerance_bp ?? 0,
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
    [CONFIG_OVERRIDE_PATHS.labelBlacklist]: state.filterMode.value === 'Blacklist'
      ? String(state.manualBlacklist.value || '').split(/[,\n]/)
        .map((keyword) => keyword.trim()).filter(Boolean)
      : [],
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
                  comparisonHeight.status === 'auto' ? null : comparisonHeight.value,
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
  const preservedOverrides = state.unmanagedConfigOverrides;
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

const addGeneratedTableResources = (state, resources, diagramOptions) => {
  recordStructuralMetric('generatedTableBuildCount');
  const paletteName = String(state.selectedPalette.value || 'default');
  const paletteColors = state.canonicalPublicationFiles && !state.canonicalPublicationFiles.d_color ? {}
    : state.normalizePaletteColors(state.paletteDefinitions.value?.[paletteName]
      || state.paletteDefinitions.value?.default || {});
  const defaultColors = buildDefaultColorOverrideTsv({
    colors: state.currentColors.value,
    paletteColors
  });
  const specificColors = serializeSpecificRules(state.manualSpecificRules);
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

  const visibility = serializeFeatureVisibilityRules(state.featureVisibilityRules?.value || []);
  if (visibility.trim()) {
    diagramOptions.featureVisibilityTableFile = fileRef(resources.addText(
      'feature-visibility-table-file', 'feature-visibility-table-file', 'feature-visibility.tsv', visibility
    ));
  }
  const preservedWhitelist = publicationFileRef(
    resources, publicationFiles, 'whitelist', 'label-whitelist-file'
  );
  if (preservedWhitelist) {
    diagramOptions.labelWhitelistFile = preservedWhitelist;
  } else if (state.filterMode.value === 'Whitelist' && state.manualWhitelist.length > 0) {
    const whitelist = state.manualWhitelist
      .filter((rule) => rule?.feat && rule?.qual)
      .map((rule) => `${rule.feat}\t${rule.qual}\t${rule.key || ''}`)
      .join('\n');
    if (whitelist) {
      diagramOptions.labelWhitelistFile = fileRef(resources.addText(
        'label-whitelist-file', 'label-whitelist-file', 'label-whitelist.tsv', `${whitelist}\n`
      ));
    }
  }
  const priority = state.manualPriorityRules
    .filter((rule) => rule?.feat && rule?.order)
    .map((rule) => `${rule.feat}\t${rule.order}`)
    .join('\n');
  const priorityRef = publicationFileRef(
    resources, publicationFiles, 'qualifier_priority', 'qualifier-priority-file'
  );
  if (priorityRef) {
    diagramOptions[
      priorityRef.representation === 'canonicalTsv'
        ? 'qualifierPriorityTable' : 'qualifierPriorityFile'
    ] = priorityRef;
  } else if (priority) {
    diagramOptions.qualifierPriorityFile = fileRef(resources.addText(
      'qualifier-priority-file', 'qualifier-priority-file', 'qualifier-priority.tsv', `${priority}\n`
    ));
  }
  // Generate supplies the table it already built for this operation (CW-02).
  const labelOverrideTsv = serializeLabelOverrideRows(state.canonicalLabelOverrideRows?.value)
    || (typeof state.generatedLabelOverrideTsv === 'string' ? state.generatedLabelOverrideTsv
      : buildLabelOverrideTsv(
        state.labelTextFeatureOverrides,
        state.labelTextBulkOverrides,
        {
          editableLabels: state.editableLabels?.value || [],
          extractedFeatures: state.extractedFeatures.value,
          featureOverrideSources: state.labelTextFeatureOverrideSources,
          visibilityOverrides: state.labelVisibilityOverrides
        }
      ).tsv);
  if (labelOverrideTsv) {
    diagramOptions.labelOverrideFile = fileRef(resources.addText(
      'label-override-file', 'label-override-file', 'label-overrides.tsv', labelOverrideTsv
    ));
  }
};

const depthSourceRowsForRequest = ({ state, filesData, recordCount }) => (
  state.mode.value === 'linear'
    ? (filesData.linearSeqs || []).map((seq) => (
        Array.isArray(seq.depth) ? seq.depth : (seq.depth ? [seq.depth] : [])
      ))
    : normalizeRecordMajorDepthFileRows(filesData.c_depth, recordCount)
);

const logicalDepthTrackCountForRequest = ({ state, filesData, recordCount }) => {
  const rows = depthSourceRowsForRequest({ state, filesData, recordCount });
  return rows.reduce((maximum, row) => (
    Math.max(maximum, Array.isArray(row) ? row.length : 0)
  ), 0);
};

const buildDepthResources = ({ state, filesData, resources, diagramOptions, recordCount }) => {
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
  const tracks = Array.isArray(state.adv.depth_tracks) ? state.adv.depth_tracks : [];
  const logicalTrackCount = Math.max(
    tracks.length,
    ...rows.map((row) => row.length)
  );
  diagramOptions.depthTracks = Array.from({ length: logicalTrackCount }, (_, trackIndex) => {
    const sources = rows.map((row) => row[trackIndex] || null);
    if (!sources.some(Boolean)) {
      throw new Error(
        `Depth series #${trackIndex + 1} (logical track index ${trackIndex}) has no source in any record.`
      );
    }
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
      color: String(track.color || state.adv.depth_color || '#4A90E2'),
      height: state.mode.value === 'linear'
        ? projectOptionalNumber(track.height, { field: 'height', seriesIndex: trackIndex })
        : null,
      largeTickInterval: projectOptionalNumber(track.large_tick_interval, { field: 'large_tick_interval', seriesIndex: trackIndex }),
      smallTickInterval: projectOptionalNumber(track.small_tick_interval, { field: 'small_tick_interval', seriesIndex: trackIndex }),
      tickFontSize: projectOptionalNumber(track.tick_font_size, { field: 'tick_font_size', seriesIndex: trackIndex })
    };
  });
};

const circularGeometryShortcutsForState = (state) => ({
  featureWidth: state.adv.feature_width_circular,
  depthWidth: state.adv.depth_width_circular,
  gcContentWidth: state.adv.gc_content_width_circular,
  gcContentRadius: state.adv.gc_content_radius_circular,
  gcSkewWidth: state.adv.gc_skew_width_circular,
  gcSkewRadius: state.adv.gc_skew_radius_circular
});

const annotationSetIdsForState = (state) => (
  (Array.isArray(state.annotationSets) ? state.annotationSets : [])
    .map((set) => String(set?.id || '').trim())
    .filter(Boolean)
);

const automaticCircularAnnotationSlots = (state) => {
  const slots = [];
  (Array.isArray(state.annotationSets) ? state.annotationSets : [])
    .forEach((set, index) => {
      const setId = String(set?.id || '').trim();
      if (!setId) return;
      const marks = Array.from(new Set(
        (Array.isArray(set.annotations) ? set.annotations : [])
          .map((annotation) => String(annotation?.mark || '').trim().toLowerCase())
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

const conservationDiagramOptions = (conservation, entries, referenceDefault) => ({
  conservationReference: String(conservation.reference || referenceDefault),
  conservationLabels: entries.map((entry) => entry.label),
  conservationColors: entries.map((entry) => entry.color),
  conservationRingWidth: projectOptionalNumber(conservation.ring_width, { field: 'conservation_ring_width' }),
  conservationRingGap: projectOptionalNumber(conservation.ring_gap, { field: 'conservation_ring_gap' })
});

const conservationSeriesForValidation = ({
  state,
  filesData,
  resolvedCircularConservation
}) => {
  if (state.circularConservation?.enabled !== true) return [];
  const conservationSource = String(
    state.circularConservation?.source || ''
  ).trim().toLowerCase();
  const blastsAreDerived = filesData.c_conservation_blasts_source === 'losat-cache';
  const sourceFiles = (
    conservationSource === 'upload' || blastsAreDerived
  )
    ? filesData.c_conservation_blasts
    : filesData.c_conservation_fastas;
  let ordered = [];
  try {
    ordered = orderedConservationSources(
      Array.isArray(sourceFiles) ? sourceFiles : [],
      state.circularConservation || {}
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

const buildTrackPlan = ({
  state,
  filesData,
  recordCount,
  resolvedCircularConservation
}) => {
  const circular = state.mode.value === 'circular';
  const depthTrackCount = logicalDepthTrackCountForRequest({
    state,
    filesData,
    recordCount
  });
  const annotationSetIds = annotationSetIdsForState(state);
  const visibleFeatureUnderlays = visibleFeatureUnderlaysForState(state);

  if (circular && state.adv.circular_track_slots_enabled) {
    const validation = assertValidCustomTrackPlan(validateCustomTrackPlan({
      mode: 'circular',
      slots: state.adv.circular_track_slots,
      axisIndex: state.adv.circular_track_slots_axis_index,
      trackType: state.form.track_type,
      depthTrackCount,
      annotationSetIds,
      visibleFeatureUnderlays,
      conservationSeries: conservationSeriesForValidation({
        state,
        filesData,
        resolvedCircularConservation
      })
    }));
    const depthRequested = validation.enabledSlots.some(
      (slot) => slot.renderer === 'depth'
    );
    return {
      depthRequested,
      tracks: {
        circularTrackSlots: validation.enabledSlots.map((slot) => (
          buildCircularTrackSlotPayload(slot, state.adv.nt, state.form.track_type)
        )),
        circularTrackAxisIndex: validation.emittedAxisIndex,
        linearTrackSlots: null,
        linearTrackAxisIndex: null,
        centerReservedRadius: projectOptionalNumber(state.adv.center_reserved_radius, { field: 'center_reserved_radius' })
      }
    };
  }

  if (!circular && state.adv.linear_track_slots_enabled) {
    const validation = assertValidCustomTrackPlan(validateCustomTrackPlan({
      mode: 'linear',
      slots: state.adv.linear_track_slots,
      axisIndex: state.adv.linear_track_slots_axis_index,
      trackType: state.form.linear_track_layout,
      depthTrackCount,
      annotationSetIds,
      visibleFeatureUnderlays,
      conservationSeries: []
    }));
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

  const depthRequested = Boolean(state.form.show_depth);
  if (circular) {
    const shortcuts = circularGeometryShortcutsForState(state);
    if (hasCircularGeometryShortcuts(shortcuts)) {
      const implicitSlots = applyCircularGeometryShortcuts(
        createDefaultCircularTrackSlots({
          nt: state.adv.nt,
          showDepth: depthRequested,
          depthTrackCount: Math.max(1, depthTrackCount),
          showGc: !state.form.suppress_gc,
          showSkew: !state.form.suppress_skew,
          showTicks: state.form.show_scale !== false,
          preset: state.form.track_type
        }),
        shortcuts
      );
      const slots = [
        ...automaticCircularAnnotationSlots(state),
        ...implicitSlots
      ];
      const axisIndex = inferLegacyAxisIndexFromFeature(
        slots,
        state.form.track_type
      );
      const validation = assertValidCustomTrackPlan(validateCustomTrackPlan({
        mode: 'circular',
        slots,
        axisIndex,
        trackType: state.form.track_type,
        depthTrackCount,
        annotationSetIds,
        visibleFeatureUnderlays,
        conservationSeries: []
      }));
      return {
        depthRequested,
        tracks: {
          circularTrackSlots: validation.enabledSlots.map((slot) => (
            buildCircularTrackSlotPayload(slot, state.adv.nt, state.form.track_type)
          )),
          circularTrackAxisIndex: validation.emittedAxisIndex,
          linearTrackSlots: null,
          linearTrackAxisIndex: null,
          centerReservedRadius: projectOptionalNumber(state.adv.center_reserved_radius, { field: 'center_reserved_radius' })
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
      centerReservedRadius: circular ? projectOptionalNumber(state.adv.center_reserved_radius, { field: 'center_reserved_radius' }) : null
    }
  };
};

const generatedProteinSettings = (state, baseline = {}) => {
  const { alignOrthogroupFeature: _legacyAlignment, ...currentBaseline } = baseline;
  const blastp = state.losat.blastp || {};
  const blastpMode = requireCurrentProteinBlastpMode(blastp.mode);
  const positiveInteger = (value, fallback) => integerSettingOr(value, fallback, 1);
  const nonNegativeInteger = (value, fallback) => integerSettingOr(value, fallback, 0);
  const rawCollinearityUnitMode = String(blastp.collinearUnitMode || '').trim().toLowerCase();
  const collinearityUnitMode = blastpMode === 'collinear'
    ? requireCurrentCollinearUnitMode(rawCollinearityUnitMode)
    : (['auto', 'cds', 'locus'].includes(rawCollinearityUnitMode)
        ? rawCollinearityUnitMode
        : 'auto');
  const rawCollinearityColorMode = String(
    blastp.collinearColorMode || ''
  ).trim().toLowerCase().replace(/-/g, '_');
  const resolvedCollinearityColorMode = rawCollinearityColorMode === 'identity'
    ? 'average_identity'
    : rawCollinearityColorMode;
  const collinearityColorMode = blastpMode === 'collinear'
    ? requireCurrentCollinearColorMode(resolvedCollinearityColorMode)
    : ([
        'average_identity',
        'orientation',
        'orientation_identity'
      ].includes(resolvedCollinearityColorMode)
        ? resolvedCollinearityColorMode
        : 'orientation');
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
      : normalizeCollinearAnchorMode(blastp.collinearAnchorMode),
    collinearitySearchScope: blastpMode === 'collinear'
      ? requireCurrentCollinearSearchScope(blastp.collinearSearchScope)
      : normalizeCollinearSearchScope(blastp.collinearSearchScope),
    collinearityColorMode,
    losatpBin: baseline.losatpBin || 'losat',
    ncbiBlastpBin: baseline.ncbiBlastpBin ?? null,
    losatpThreads: integerSettingOr(state.losat.threadsPerJob, null, 1),
    proteinBlastpMaxHits: blastpMode === 'pairwise'
      ? requireCurrentProteinBlastpMaxHits(blastp.maxHits)
      : positiveInteger(blastp.maxHits, 5),
    proteinBlastpCandidateLimit: requireCurrentProteinBlastpCandidateLimit(
      blastp.candidateLimit
    ),
    orthogroupMembershipMode: ['orthogroup', 'collinear'].includes(blastpMode)
      ? requireCurrentOrthogroupMembershipMode(blastp.orthogroupMembershipMode)
      : normalizeOrthogroupMembershipMode(blastp.orthogroupMembershipMode),
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

const comparisonPlanErrorMessage = (snapshot) => {
  const direct = String(snapshot?.error || '').trim();
  if (direct) return direct;
  const issue = Array.isArray(snapshot?.errors) ? snapshot.errors[0] : null;
  return String(issue?.message || issue || '').trim();
};

const requireLinearComparisonPlanSnapshot = (snapshot) => {
  if (!snapshot || !Array.isArray(snapshot.edges)) {
    throw new Error('A resolved Linear comparison plan is required.');
  }
  const error = comparisonPlanErrorMessage(snapshot);
  if (error) throw new Error(error);
  return snapshot;
};

const orderedComparisonPlanEdges = (snapshot) => (
  (Array.isArray(snapshot?.edges) ? snapshot.edges : [])
    .slice()
    .sort((left, right) => Number(left?.ordinal) - Number(right?.ordinal))
);

const comparisonEndpointKey = (queryRecordIndex, subjectRecordIndex) => (
  `${Number(queryRecordIndex)}->${Number(subjectRecordIndex)}`
);

const addResolvedComparisonResource = ({
  comparison,
  edge,
  resources
}) => {
  const resourceId = `comparison-resolved-${edge.ordinal + 1}`;
  if (comparison.kind === 'precomputedProteinComparison') {
    return {
      kind: 'precomputedProteinComparison',
      resourceId: resources.addCanonicalTable(
        resourceId,
        comparison.rows,
        comparison.columns
      ),
      encoding: 'canonicalTsv',
      queryRecordIndex: edge.queryIndex,
      subjectRecordIndex: edge.subjectIndex
    };
  }
  if (comparison.kind !== 'nucleotideBlast') {
    throw new Error(`Unsupported resolved comparison kind for '${edge.edgeKey}'.`);
  }
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

const addPersistedProteinComparisonResource = ({
  comparison,
  edge,
  resources
}) => ({
  kind: 'precomputedProteinComparison',
  resourceId: resources.addFile(
    `comparison-canonical-protein-${edge.ordinal + 1}`,
    canonicalComparisonResourceKind(comparison),
    comparison.file
  ),
  encoding: String(comparison.encoding || 'canonicalTsv'),
  queryRecordIndex: edge.queryIndex,
  subjectRecordIndex: edge.subjectIndex
});

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

const buildComparisons = ({
  state,
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
      || ![1, 2, 3].includes(typedResource.schema)
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
    snapshot.hasLosatIntent === true && state.losatProgram?.value === 'blastp'
  );
  const persistedCanonicalInputs = persistedCanonicalComparisons.filter(
    (comparison) => comparison?.canonicalInput === true
  );
  const persistedMetadata = persistedCanonicalComparisons.filter((comparison) => (
    (
      (
        comparison?.kind === 'orthogroupResult'
        && String(state.losat?.blastp?.mode || '').trim().toLowerCase() === 'orthogroup'
      ) || (
        comparison?.kind === 'collinearityResult'
        && String(state.losat?.blastp?.mode || '').trim().toLowerCase() === 'collinear'
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
          'comparison-canonical-orthogroups-1',
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
          'comparison-canonical-collinearity-1',
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
      ? 'comparison-canonical-collinearity-1'
      : 'comparison-canonical-orthogroups-1';
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
        throw new Error(`The uploaded comparison '${edge.edgeKey}' has no active BLAST TSV file.`);
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
    String(state.losat?.blastp?.mode || '').trim().toLowerCase() === 'pairwise'
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
        : String(state.losat?.blastp?.mode || persistedGeneratedComparison?.mode || 'orthogroup');
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
        state,
        persistedGeneratedComparison?.settings || {}
      )
    });
  }
  return comparisons;
};

const buildLayout = (state, filesData, records = []) => {
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
    if (!state.linearRecordLayoutEnabled?.value &&
        !records.some((record) => record.presentation.gridRow != null) &&
        plan === null && translations.length === 0) return {};
    return {
      recordGapPx: Math.max(0, Number(state.linearRecordGap?.value) || 0),
      multiRecordPositions: null,
      recordTranslations: translations,
      similarityAlignment: plan
    };
  }
  if (!state.form.multi_record_canvas) return {};
  const positions = Array.isArray(state.adv.multi_record_positions)
    ? state.adv.multi_record_positions
        .map((entry) => {
          const selector = String(entry?.selector || '').trim();
          const row = Number(entry?.row);
          return selector && Number.isInteger(row) && row > 0 ? `${selector}@${row}` : null;
        })
        .filter(Boolean)
    : [];
  return {
    multiRecordSizeMode: requireCurrentCircularMultiRecordSizeMode(
      state.adv.multi_record_size_mode
    ),
    // Blank is the documented default (Python's); any number is sent literally.
    multiRecordMinRadiusRatio: projectOptionalNumber(state.adv.multi_record_min_radius_ratio, { field: 'multi_record_min_radius_ratio' }) ?? 0.55,
    multiRecordColumnGapRatio: projectOptionalNumber(state.adv.multi_record_column_gap_ratio, { field: 'multi_record_column_gap_ratio' }) ?? 0.10,
    multiRecordRowGapRatio: projectOptionalNumber(state.adv.multi_record_row_gap_ratio, { field: 'multi_record_row_gap_ratio' }) ?? 0.05,
    multiRecordPositions: positions.length > 0 ? positions : null
  };
};

const publicationClone = (value) => value === undefined ? undefined : JSON.parse(JSON.stringify(value));
const publicationRef = (value) => ({ value });
const canonicalRecordSelector = (record) => {
  const selector = record?.region?.selector || record?.selector;
  if (selector?.kind === 'recordId') return String(selector.value || '');
  if (selector?.kind === 'recordIndex') return `#${Number(selector.index) + 1}`;
  return null;
};
export const buildCanonicalRequestState = ({ session, projection, config,
  filesData = projection.files }) => {
  const canonicalPublicationFiles = { ...filesData };
  const pythonColors = (colors) => Object.fromEntries(Object.entries(colors || {}).filter(
    ([key]) => !key.startsWith('collinear_block_')).sort(([left], [right]) => left.localeCompare(right)));
  const activeColors = pythonColors(config?.colors);
  const colorOverridesChanged = config?.colorsAreOverrides === true
    && JSON.stringify(activeColors) !== JSON.stringify(pythonColors(projection.config?.colors));
  if (colorOverridesChanged) delete canonicalPublicationFiles.d_color;
  const features = session?.features || {}, layout = config?.linearRecordLayout || {};
  const canonicalLayout = projection.config?.linearRecordLayout || {};
  const legacyAlignment = projection.pipelineState?.legacySimilarityAlignment;
  const materializedLegacyPlan = legacyAlignment
      ? canonicalSimilarityAlignment(materializeLegacySimilarityAlignment({
        target: legacyAlignment.target,
        records: session?.renderRequest?.records || [],
        featureCatalog: session?.editorState?.featureCatalog || null,
        legacyOrthogroupState: session?.orthogroupState || null
      }), (session?.renderRequest?.records || []).map((record) => record.recordKey),
      'renderRequest.layout.similarityAlignment')
    : null;
  const palette = String(config?.palette || 'default');
  const records = session?.renderRequest?.records || [];
  const circularRecordList = records.length === 1 && !records[0]?.selector
    && !records[0]?.region?.selector ? []
    : records.map((record) => ({ selector: canonicalRecordSelector(record) }));
  const conservation = publicationClone(config?.circularConservation || {
    enabled: false, reference: 'auto', labels: '', series: []
  });
  if (!String(conservation.labels || '').trim() && conservation.series?.length)
    conservation.labels = conservation.series.map(({ label }) => String(label || '').trim()).join(',');
  const refs = {
    mode: projection.mode, cInputType: session?.ui?.cInputType || projection.inputType,
    lInputType: session?.ui?.lInputType || projection.inputType,
    circularRecordList,
    paletteDefinitions: { [palette]: {} }, currentColors: colorOverridesChanged ? activeColors : config.colors || {},
    selectedPalette: palette, featureVisibilityRules: publicationClone(features.featureVisibilityManualRules
      || projection.semanticFeatureState?.featureVisibilityManualRules || []),
    filterMode: config.filterMode || 'None', manualBlacklist: String(config.blacklistText || ''),
    canonicalLabelOverrideRows: publicationClone(features.labelOverrideRows || []), editableLabels: [],
    extractedFeatures: features.extractedFeatures || [],
    losatProgram: config.losatProgram || 'blastn',
    linearRecordLayoutEnabled: Boolean(layout.enabled),
    linearRecordGap: layout.recordGap ?? 24,
    similarityAlignmentPlan: materializedLegacyPlan ||
      canonicalLayout.similarityAlignment || null,
    linearRecordTranslations: materializedLegacyPlan
      ? (session?.renderRequest?.records || []).map((record) => ({
          recordKey: String(record.recordKey || ''), x: 0, y: 0
        }))
      : publicationClone(canonicalLayout.recordTranslations || [])
  };
  return {
    ...Object.fromEntries(Object.entries(refs).map(([key, value]) => [key, publicationRef(value)])),
    form: config.form || {}, adv: config.adv || {}, normalizePaletteColors,
    manualSpecificRules: publicationClone(config.rules || []), manualWhitelist: publicationClone(config.whitelist || []), manualPriorityRules: publicationClone(config.qualifierPriorityRules || []),
    labelTextFeatureOverrides: publicationClone(features.labelTextFeatureOverrides || {}), labelTextBulkOverrides: publicationClone(features.labelTextBulkOverrides || {}),
    labelTextFeatureOverrideSources: publicationClone(features.labelTextFeatureOverrideSources || {}), labelVisibilityOverrides: publicationClone(features.labelVisibilityOverrides || {}),
    circularConservation: conservation, losat: publicationClone(config.losat || { blastp: {} }),
    linearRecordRows: publicationClone(layout.rows || []), linearComparisonPlan: normalizeLinearComparisonPlan(config.linearComparisonPlan || { mode: 'none', defaultSource: 'losat', edges: [] }),
    annotationSets: publicationClone(config.annotationSets || []),
    recordDisplayDrafts: publicationClone(config.recordDisplayDrafts || []),
    featurePlacementOverrides: publicationClone(config.featurePlacementOverrides || {}), canonicalPublicationFiles
  };
};
const projectCanonicalRenderInput = ({
  state,
  filesData,
  comparisonPlanSnapshot = null,
  resolvedComparisons = [],
  resolvedCircularConservation = [],
  resources = createResourceBuilder()
}) => {
  recordStructuralMetric('canonicalRequestProjectionCount');
  requireCurrentWebStateFieldNames(state);
  requireCurrentCircularMultiRecordSizeMode(state.adv.multi_record_size_mode);
  requireCurrentLinearTrackLayout(state.form.linear_track_layout);
  requireCurrentLinearLabelPlacement(state.adv.label_placement);
  if (
    state.mode.value === 'linear' &&
    Boolean(state.form.normalize_length) &&
    Boolean(state.linearRecordLayoutEnabled?.value) &&
    linearRecordLayoutHasSharedRow(filesData?.linearSeqs, state.linearRecordRows)
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
  const webFiles = {};
  const recordPlan = buildRecords({ state, filesData, resources });
  const drafts = state.recordDisplayDrafts || [];
  const displayRows = state.recordDisplayRows?.value || [];
  const sourceInputIndexes = [];
  const records = recordPlan.records.flatMap((record, index) => {
    const sourceUid = state.mode.value === 'linear'
      ? String(filesData.linearSeqs?.[index]?.uid || record.recordKey) : 'circular';
    const sourceRows = displayRows.filter((row) => row.scope === state.mode.value && row.sourceUid === sourceUid);
    const selector = record.region?.selector || record.selector;
    const selectedRows = sourceRows.filter((row) => !selector
      || (selector.kind === 'recordIndex' ? row.selector === `#${selector.index + 1}` : row.recordId === selector.value));
    const transformFor = (row) => requestedRecordTransform(row,
      drafts.find((draft) => recordDisplayKey(draft) === row.key) || {}, { ...row, cropped: Boolean(record.region) });
    if (record.cardinality === 'all' && selectedRows.length > 1) {
      const transforms = selectedRows.map(transformFor);
      if (new Set(transforms.map((transform) => JSON.stringify(transform))).size > 1) {
        sourceInputIndexes.push(...selectedRows.map(() => index));
        return selectedRows.map((row, rowIndex) => ({ ...record,
          recordKey: `${record.recordKey}:${Number(row.selector.slice(1))}`,
          cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: Number(row.selector.slice(1)) - 1 },
          presentation: { ...record.presentation,
            reverseComplement: record.region ? false : transforms[rowIndex].reverseComplement,
            gridRow: record.presentation.gridRow ?? index + 1 },
          display: transforms[rowIndex].display }));
      }
    }
    const selected = selectedRows[0];
    const savedDraft = drafts.find((draft) => draft.scope === state.mode.value && draft.sourceUid === sourceUid
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
    recordPlan.circularSourceIndexes.length < recordPlan.circularSourceCount &&
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
    filesData: state.mode.value === 'linear'
      ? { ...filesData, linearSeqs: sourceInputIndexes.map((index) => filesData.linearSeqs[index]) }
      : selectedCircularFilesData,
    recordCount: records.length,
    resolvedCircularConservation
  });
  const circularGroupingIntent = ['single', 'batch'].includes(
    state.adv.circular_grouping_intent
  )
    ? state.adv.circular_grouping_intent
    : null;
  const grouping = state.mode.value === 'linear'
    ? 'single'
    : (
        state.form.multi_record_canvas
          ? 'grid'
          : (
              records.length === 1
                ? (
                    circularGroupingIntent ||
                    WEB_UX_PROFILE.circular.singleRecordGrouping
                  )
                : WEB_UX_PROFILE.circular.multiRecordGrouping
            )
      );
  const explicitPrefix = explicitOutputPrefix(state.form.prefix);
  if (state.mode.value === 'circular') {
    webFiles.circularOutputPrefixExplicit = explicitPrefix !== null;
  }
  const knownCircularRecords = Array.isArray(state.circularRecordList.value)
    ? state.circularRecordList.value
    : [];
  const knownCircularEntries = buildDisambiguatedRecordEntries(
    knownCircularRecords.map((record) => ({
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
  const defaultCircularPrefix = safePrefix(circularRecordId(circularOutputRecords[0], 0));
  const output = grouping === 'batch'
    ? resolveCircularBatchPrefixes(circularOutputRecords, explicitPrefix)
        .map(renderOutputPayload)
    : renderOutputPayload(
        explicitPrefix ?? (state.mode.value === 'circular' ? defaultCircularPrefix : 'out')
      );
  const diagramOptions = {
    featurePlacements: canonicalFeaturePlacements(state.featurePlacementOverrides || {}, state.mode.value),
    configOverrides: buildConfigOverrides(state, {
      depthRequested: trackPlan.depthRequested,
      hasComparisonIntent: hasLinearComparisonIntent,
      linearHasSharedRow: state.mode.value === 'linear' && linearRecordLayoutHasSharedRow(
        filesData.linearSeqs,
        state.linearRecordRows,
        { enabled: Boolean(state.linearRecordLayoutEnabled?.value) }
      )
    }),
    tracks: trackPlan.tracks,
    output: {
      legend: String(state.form.legend || 'right'),
      plotTitlePosition: String(state.adv.plot_title_position || (state.mode.value === 'linear' ? 'bottom' : 'none'))
    },
    selectedFeaturesSet: Array.from(state.adv.features || []).map((value) => String(value)),
    featureShapes: {
      repeat_region: defaultFeatureRendering('repeat_region'),
      ...normalizeFeatureRenderingMap(state.adv.feature_shapes || {})
    },
    dinucleotide: String(state.adv.nt || 'GC').toUpperCase(),
    window: projectOptionalNumber(state.adv.window_size, { field: 'window' }),
    step: projectOptionalNumber(state.adv.step_size, { field: 'step' }),
    depthWindow: projectOptionalNumber(state.adv.depth_window_size, { field: 'depth_window' }),
    depthStep: projectOptionalNumber(state.adv.depth_step_size, { field: 'depth_step' }),
    plotTitle: String(state.form.plot_title || '').trim() || null,
    plotTitleFontSize: projectOptionalNumber(state.adv.plot_title_font_size, { field: 'plot_title_font_size' }),
    ...(comparisonOptionsRequested
      ? projectComparisonThresholds(resolveComparisonThresholds(state.adv, state.mode.value))
      : {})
  };
  if (Array.isArray(state.annotationSets) && state.annotationSets.length > 0) {
    diagramOptions.annotations = annotationOptionsPayload(state.annotationSets);
  }
  if (state.mode.value === 'circular') {
    diagramOptions.keepFullDefinitionWithPlotTitle = Boolean(state.adv.keep_full_definition_with_plot_title);
    diagramOptions.species = String(state.form.species || '').trim() || null;
    diagramOptions.strain = String(state.form.strain || '').trim() || null;
    const conservationSource = String(
      state.circularConservation.source || ''
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
        state.circularConservation
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
      Object.assign(diagramOptions, conservationDiagramOptions(state.circularConservation, conservationEntries, 'auto'));
      if (conservationBlastsAreDerived) {
        webFiles.conservationBlastSource = 'losat-cache';
      }
    }
    if (conservationSource === 'losat') {
      const comparisonFastas = orderedOptionalConservationFiles(
        filesData.c_conservation_fastas,
        state.circularConservation
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
      Object.assign(diagramOptions, conservationDiagramOptions(state.circularConservation,
        resolvedCircularConservation.map((entry, index) => ({
          label: String(entry?.label || `Comparison ${index + 1}`),
          color: String(entry?.color || '#D9EAF7')
        })), 'subject'));
      webFiles.conservationBlastSource = 'losat-cache';
    }
  } else if (hasLinearComparisonIntent) {
    diagramOptions.pairwiseMatchStyle = String(state.adv.pairwise_match_style || 'ribbon');
  }

  addGeneratedTableResources(state, resources, diagramOptions);
  if (trackPlan.depthRequested) {
    buildDepthResources({
      state,
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

  return {
    renderRequest: {
      schema: CANONICAL_REQUEST_SCHEMA,
      mode: state.mode.value,
      grouping,
      records,
      diagramOptions,
      layout: buildLayout(state, filesData, records),
      comparisons: buildComparisons({
      state,
      filesData,
      resources,
      comparisonPlanSnapshot,
      resolvedComparisons
    }),
      output
    },
    resources: resources.resources,
    webFiles
  };
};

// The request writer and the inexpensive comparison use the same projection.
// The latter retains bindings and generated table text, never resource bytes.
const generationResourceIdentity = (resource) => {
  if (!resource) return null;
  if (generatedResourceValues.has(resource)) {
    const value = generatedResourceValues.get(resource);
    return typeof value === 'string' ? { text: value } : value;
  }
  const owner = getResourcePayloadOwner(resource);
  const backing = getSessionResourceSource(owner);
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
const compositionSourceIdentity = (source, resources) => {
  if (!source || typeof source !== 'object') return null;
  const fields = Object.keys(source).sort();
  const identity = fields.map((key) => {
    if (!isCanonicalResourceReferenceField(key)) return [key, source[key]];
    const descriptor = resources?.[source[key]];
    if (!descriptor) return [key, null];
    const parts = generationResourceIdentity(descriptor)?.bindings;
    if (!Array.isArray(parts) || !parts.length) return [key, null];
    return [key, parts.map(part => ({
      payload: typeof part?.data === 'string' ? part.data : part,
      size: part?.size,
      encoding: part?.encoding || 'file'
    }))];
  });
  return identity;
};

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
          ? { kind: 'recordIndex', index: 0 } : null);
        return [{ record, selector }];
      }
      const suffix = key.startsWith(`${record.recordKey}:`) ? key.slice(record.recordKey.length + 1) : '';
      if (record.cardinality !== 'all' || !/^[1-9]\d*$/.test(suffix)) return [];
      return [{ record, selector: { kind: 'recordIndex', index: Number(suffix) - 1 } }];
    });
    if (matches.length !== 1) return null;
    const { record, selector } = matches[0];
    const source = compositionSourceIdentity(record.source, canonical.resources);
    if (!source || source.some(([field, value]) => isCanonicalResourceReferenceField(field)
      && (!value || value.some(part => !part.payload || !Number.isSafeInteger(part.size))))) return null;
    records.push({ key, source, selector: publicationClone(selector),
      region: publicationClone(record.region ?? null) });
  }
  return { mode: request.mode, grouping: request.grouping, records };
};

export const buildCanonicalRenderRequest = (args) => projectCanonicalRenderInput(args);

const recordSourceResourceId = (record, field) => {
  const source = record?.source || {};
  if (field === 'gb' && source.kind === 'genbank') return String(source.resourceId || '');
  if (field === 'gff' && source.kind === 'gffFasta') return String(source.gffResourceId || '');
  if (field === 'fasta' && source.kind === 'gffFasta') return String(source.fastaResourceId || '');
  return '';
};

const referencedResourceId = (ref) => String(ref?.resourceId || '').trim();

const addResourceOriginalNameHint = (target, resourceId, name) => {
  const id = String(resourceId || '').trim();
  const originalName = normalizeOriginalResourceName(name);
  if (!id || !originalName || Object.prototype.hasOwnProperty.call(target, id)) return;
  target[id] = originalName;
};

const legacyResourceOriginalNames = ({ renderRequest, legacyFiles, fileBindings }) => {
  const hints = {};
  const records = Array.isArray(renderRequest?.records) ? renderRequest.records : [];
  const options = renderRequest?.diagramOptions || {};
  const files = legacyFiles && typeof legacyFiles === 'object' && !Array.isArray(legacyFiles)
    ? legacyFiles
    : {};
  const namedOptionResources = {
    d_color: referencedResourceId(options.colors?.defaultColorsFile || options.colors?.defaultColors),
    t_color: referencedResourceId(options.colors?.colorTableFile || options.colors?.colorTable),
    whitelist: referencedResourceId(options.labelWhitelistFile),
    qualifier_priority: referencedResourceId(
      options.qualifierPriorityFile || options.qualifierPriorityTable
    )
  };
  Object.entries(namedOptionResources).forEach(([slot, resourceId]) => {
    addResourceOriginalNameHint(hints, resourceId, files?.[slot]?.name);
  });

  if (renderRequest?.mode === 'linear') {
    const sequences = Array.isArray(files.linearSeqs) ? files.linearSeqs : [];
    records.forEach((record, index) => {
      const sequence = sequences[index] || {};
      ['gb', 'gff', 'fasta'].forEach((field) => {
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

  (Array.isArray(fileBindings) ? fileBindings : []).forEach((binding) => {
    const slot = String(binding?.slot || '');
    const normalizedSlot = slot.replace(/^files\./, '');
    if (Object.prototype.hasOwnProperty.call(namedOptionResources, normalizedSlot)) {
      addResourceOriginalNameHint(
        hints,
        namedOptionResources[normalizedSlot],
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

const resourceAsLegacyFile = (resources, resourceId) => {
  const entry = resources?.[resourceId];
  if (!entry || typeof entry !== 'object') throw new Error(`Missing canonical resource: ${resourceId}`);
  const { kind: _kind, ...file } = entry;
  return file;
};

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

const webBindingValueAsLegacyFile = (resources, value, resolveResourceFile = null, schema = 1) => (
  Array.isArray(value)
    ? value.map((item) => webBindingValueAsLegacyFile(resources, item, resolveResourceFile, schema))
    : webBindingAsLegacyFile(resources, value, resolveResourceFile, schema)
);

const applyWebFileBindings = (
  files,
  webMetadata,
  resources,
  { resolveResourceFile = null, sessionResourceTable = null, adoptCanonicalPayloads = false } = {}
) => {
  const bindings = webMetadata?.bindings;
  if (bindings === undefined) return files;
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
    restored.linearSeqs = bindings.linearSeqs.map((sequence, index) => ({
      uid: String(sequence?.uid || `canonical-seq-${index + 1}`),
      gb: webBindingValueAsLegacyFile(resources, sequence?.gb, resolveBinding, bindings.schema),
      gff: webBindingValueAsLegacyFile(resources, sequence?.gff, resolveBinding, bindings.schema),
      fasta: webBindingValueAsLegacyFile(resources, sequence?.fasta, resolveBinding, bindings.schema),
      depth: webBindingValueAsLegacyFile(resources, sequence?.depth, resolveBinding, bindings.schema),
      blast: webBindingValueAsLegacyFile(resources, sequence?.blast, resolveBinding, bindings.schema),
      losat_gencode: integerSettingOr(sequence?.losat_gencode, 1, 1),
      losat_filename: String(sequence?.losat_filename || ''),
      definition: String(sequence?.definition || ''),
      record_subtitle: String(sequence?.record_subtitle || ''),
      file_definition: String(sequence?.file_definition || ''),
      file_subtitle: String(sequence?.file_subtitle || ''),
      inferred_definition: String(sequence?.inferred_definition || ''),
      region_record_id: String(sequence?.region_record_id || ''),
      region_start: sequence?.region_start ?? null,
      region_end: sequence?.region_end ?? null,
      region_reverse: Boolean(sequence?.region_reverse)
    }));
  }
  if (Array.isArray(bindings.linearComparisons)) {
    restored.linearComparisons = bindings.linearComparisons.map((comparison, index) => ({
      ...(adoptCanonicalPayloads ? comparison : cloneCanonicalJsonValue(comparison)),
      id: String(comparison?.id || `linear-comparison-restored-${index + 1}`),
      queryUid: String(comparison?.queryUid || ''),
      subjectUid: String(comparison?.subjectUid || ''),
      source: String(comparison?.source || 'upload'),
      file: webBindingValueAsLegacyFile(resources, comparison?.file, resolveBinding, bindings.schema)
    }));
  }
  if (Array.isArray(bindings.linearCanonicalComparisons)) {
    restored.linearCanonicalComparisons = bindings.linearCanonicalComparisons.map((comparison) => ({
      ...(adoptCanonicalPayloads ? comparison : cloneCanonicalJsonValue(comparison)),
      file: webBindingValueAsLegacyFile(resources, comparison?.file, resolveBinding, bindings.schema)
    }));
  }
  return restored;
};

// Called only after Session authority has admitted the explicit source-free
// document. No canonical request or committed owner is created by this projection.
export const projectSettingsOnlySession = (data, sessionResourceTable) => ({
  mode: data.ui.mode,
  inputType: data.ui.mode === 'linear' ? data.ui.lInputType : data.ui.cInputType,
  config: cloneCanonicalJsonValue(data.config),
  files: applyWebFileBindings({}, data.webFiles, data.resources, {
    sessionResourceTable, adoptCanonicalPayloads: true
  }),
  semanticFeatureState: {}
});

const cloneCanonicalJsonValue = (value) => (
  value === undefined ? undefined : JSON.parse(JSON.stringify(value))
);

// Web controls own row membership and record order, not a distinct numeric
// column value. Preserve the engine's row/column render order in the projected
// record sequence, then retire the unsupported numeric column in the draft.
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

const resolvePipelineCollinearInference = (settings, mode) =>
  settings.collinearInferOrthogroups ?? ['collinear', 'none'].includes(mode);

const projectGeneratedProteinPipeline = (
  comparison,
  { adoptCanonicalPayloads = false, requestSchema = CANONICAL_REQUEST_SCHEMA } = {}
) => {
  if (
    !comparison ||
    comparison.kind !== 'generatedProteinComparison' ||
    !comparison.settings ||
    typeof comparison.settings !== 'object' ||
    Array.isArray(comparison.settings)
  ) return null;
  const settings = comparison.settings;
  if (requestSchema >= 8 && Object.hasOwn(settings, 'alignOrthogroupFeature')) {
    throw new Error('Current canonical protein settings cannot contain legacy alignment state.');
  }
  let legacySimilarityAlignment = null;
  if (requestSchema <= 7 && settings.alignOrthogroupFeature !== null &&
      settings.alignOrthogroupFeature !== undefined) {
    legacySimilarityAlignment = {
      target: requireCanonicalText(
        settings.alignOrthogroupFeature,
        'renderRequest.comparisons[].settings.alignOrthogroupFeature'
      ),
      sourceSchema: requestSchema
    };
  }
  const parameters = settings.collinearityParams?.parameters || {};
  const mode = String(comparison.mode || 'orthogroup');
  return {
    generatedProteinComparison: adoptCanonicalPayloads
      ? comparison
      : cloneCanonicalJsonValue(comparison),
    legacySimilarityAlignment,
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

const hasOptionValue = (value) => value !== null && value !== undefined;

const requireExactCanonicalFields = (value, required, fieldName) => {
  const keys = Object.keys(value);
  const missing = required.filter(
    (key) => !Object.prototype.hasOwnProperty.call(value, key)
  );
  if (missing.length > 0) {
    throw new Error(`${fieldName} is missing required field(s): ${missing.join(', ')}.`);
  }
  const unknown = keys.filter((key) => !required.includes(key));
  if (unknown.length > 0) {
    throw new Error(`${fieldName} contains unknown field(s): ${unknown.join(', ')}.`);
  }
};

const canonicalDepthResourceFile = (
  ref,
  resources,
  fieldName,
  resolveResourceFile = null
) => {
  if (!ref || typeof ref !== 'object' || Array.isArray(ref)) {
    throw new Error(`${fieldName} must be a canonical resource reference.`);
  }
  requireExactCanonicalFields(ref, ['resourceId', 'representation'], fieldName);
  if (!['file', 'canonicalTsv'].includes(ref.representation)) {
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

const canonicalDepthNumber = (value, fieldName) => {
  if (value === null) return null;
  if (typeof value !== 'number' || !Number.isFinite(value) || value <= 0) {
    throw new Error(`${fieldName} must be null or a positive finite number.`);
  }
  return value;
};

const canonicalDepthText = (value, fieldName) => {
  if (value === null) return null;
  if (typeof value !== 'string' || !value.trim()) {
    throw new Error(`${fieldName} must be null or a non-empty string.`);
  }
  return value.trim();
};

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

  const sourceRows = Array.from({ length: records.length }, () => []);
  const fileRows = Array.from({ length: records.length }, () => []);
  const tracks = options.depthTracks.map((track, trackIndex) => {
    const fieldName = `diagramOptions.depthTracks[${trackIndex}]`;
    if (!track || typeof track !== 'object' || Array.isArray(track)) {
      throw new Error(`${fieldName} must be an object.`);
    }
    requireExactCanonicalFields(
      track,
      [
        'source',
        'label',
        'color',
        'height',
        'largeTickInterval',
        'smallTickInterval',
        'tickFontSize'
      ],
      fieldName
    );
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
      throw new Error(
        `Depth series #${trackIndex + 1} (logical track index ${trackIndex}) has no source in any record.`
      );
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

export const decodeCanonicalResourceText = (resources, resourceId) => {
  const entry = resources?.[resourceId];
  if (!entry || typeof entry !== 'object' || Array.isArray(entry)) {
    throw new Error(`Missing canonical resource: ${resourceId}`);
  }
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

const sourceRecordCounts = new WeakMap();

export const readCanonicalResourceRecordCount = async (resources, resourceId, kind) => {
  const descriptor = resources?.[resourceId];
  const owner = getResourcePayloadOwner(descriptor);
  const hasFileOwner = owner && owner !== descriptor;
  let counts = hasFileOwner && sourceRecordCounts.get(owner);
  if (!counts) {
    const text = hasFileOwner
      ? await readFileText(owner)
      : decodeCanonicalResourceText(resources, resourceId);
    counts = {
      genbank: (text.match(/^LOCUS\s+/gm) || []).length,
      fasta: (text.match(/^>/gm) || []).length
    };
    // Immutable File/view identity outlives transferred bytes. Retain only counts;
    // a replacement source gets its own entry, and ownerless descriptors stay uncached.
    if (hasFileOwner) sourceRecordCounts.set(owner, counts);
  }
  return counts[kind];
};

const resourceTextFromRef = (resources, ref) => (
  ref?.resourceId ? decodeCanonicalResourceText(resources, ref.resourceId) : null
);

const nestedConfigValue = (config, path) => {
  let current = config;
  for (const key of path.split('.')) {
    if (!current || typeof current !== 'object' || Array.isArray(current)) return undefined;
    current = current[key];
  }
  return current;
};

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

const legacyLabelScope = ({ mode, showLabels, allowInnerLabels }) => {
  let scope;
  if (showLabels !== undefined) {
    if (typeof showLabels === 'boolean') {
      scope = showLabels ? (mode === 'circular' ? 'outer' : 'all') : 'none';
    } else {
      const normalized = String(showLabels).trim().toLowerCase();
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

const projectFullConfigOverrides = (config, mode) => {
  if (!config || typeof config !== 'object' || Array.isArray(config)) return {};
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
  const projected = {};
  paths.forEach((path) => {
    const value = nestedConfigValue(config, path);
    if (value !== undefined) projected[path] = value;
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

const projectCanonicalConfigOverrides = (overrides, mode) => {
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

  const lineStyles = {};
  for (const [kind, prefix] of Object.entries(LINEAR_DEFINITION_STYLE_PATHS)) {
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

const projectLegacyFlatConfigOverrides = (overrides) => Object.fromEntries(
  Object.entries(overrides).filter(([key]) => (
    !key.includes('.') && !['show_labels', 'allow_inner_labels'].includes(key)
  ))
);

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

const projectCircularConservationConfig = (options, files) => {
  const sourceFiles = Array.isArray(files.c_conservation_blasts)
    ? files.c_conservation_blasts
    : [];
  if (sourceFiles.length === 0) return undefined;
  const labels = Array.isArray(options.conservationLabels)
    ? options.conservationLabels.map((value) => String(value || '').trim())
    : [];
  const colors = Array.isArray(options.conservationColors)
    ? options.conservationColors.map((value) => String(value || '').trim())
    : [];
  const series = sourceFiles.map((file, index) => {
    const fileName = String(file?.name || `comparison-${index + 1}.tsv`);
    const defaultLabel = fileName.replace(/\.[^.]+$/, '').trim() || `Comparison ${index + 1}`;
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

const projectCanonicalCircularPixel = (measure) => {
  // Historical structured gaps/spacing are physical pixels, including zero.
  if (measure && typeof measure === 'object' && !Array.isArray(measure)) {
    if (String(measure.unit || '').trim().toLowerCase() !== 'px') {
      throw new Error('Circular gap/spacing must use pixels.');
    }
    measure = measure.value;
  }
  const value = parseOptionalPixel(measure, 'Circular gap/spacing', { allowZero: true });
  return value === null ? null : String(value);
};

const projectCanonicalCircularSlot = (slot) => ({
  ...slot,
  width: projectCircularMeasureDraft(slot?.width),
  radius: projectCircularMeasureDraft(slot?.radius),
  inner_gap_px: projectCanonicalCircularPixel(
    slot?.innerGapPx ?? slot?.inner_gap_px
  ),
  outer_gap_px: projectCanonicalCircularPixel(
    slot?.outerGapPx ?? slot?.outer_gap_px
  )
});

const projectCurrentCanonicalCircularSlot = (slot) => {
  if (
    !slot || typeof slot !== 'object' || Array.isArray(slot) ||
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
  const projected = projectCanonicalCircularSlot(slot);
  return {
    id: String(projected.id || ''),
    renderer: String(projected.renderer || ''),
    enabled: projected.enabled !== false,
    width: projected.width ?? null, radius: projected.radius ?? null,
    inner_gap_px: slot.innerGapPx == null ? null : String(slot.innerGapPx),
    outer_gap_px: slot.outerGapPx == null ? null : String(slot.outerGapPx),
    side: projected.side ?? null, z: Number(projected.z) || 0,
    params: cloneCanonicalJsonValue(projected.params || {})
  };
};
const projectLegacyCanonicalCircularSlot = (slot) => {
  const projected = projectCanonicalCircularSlot(slot);
  if (Object.prototype.hasOwnProperty.call(slot, 'spacing')) {
    projected.spacing = projectCanonicalCircularPixel(slot.spacing);
  }
  return migrateLegacyCircularTrackSlot(projected);
};

const combineCircularGenbankResources = (
  resources,
  records,
  originalName = '',
  { resolveResourceFile = null, sessionResourceTable = null } = {}
) => {
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

const validateCanonicalAssemblyOutput = (value, schema) => {
  if (value === undefined || value === null) return;
  if (!value || typeof value !== 'object' || Array.isArray(value)) {
    throw new Error('Canonical diagramOptions.output must be an object.');
  }
  const required = schema <= 2
    ? ['outputPrefix', 'legend', 'plotTitlePosition']
    : ['legend', 'plotTitlePosition'];
  const keys = Object.keys(value);
  const missing = required.filter((key) => !Object.prototype.hasOwnProperty.call(value, key));
  if (missing.length > 0) {
    throw new Error(
      `Missing required canonical diagramOptions.output field(s): ${missing.join(', ')}.`
    );
  }
  const unknown = keys.filter((key) => !required.includes(key));
  if (unknown.length > 0) {
    throw new Error(
      `Unknown canonical diagramOptions.output field(s): ${unknown.join(', ')}.`
    );
  }
};

const canonicalGrouping = (renderRequest, records) => {
  if (renderRequest.schema <= 2) {
    if (renderRequest.mode === 'linear') return 'single';
    return Object.keys(renderRequest.layout || {}).length > 0 || records.length > 1
      ? 'grid'
      : 'single';
  }
  const allowed = renderRequest.mode === 'circular'
    ? ['single', 'grid', 'batch']
    : ['single'];
  if (!allowed.includes(renderRequest.grouping)) {
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

const canonicalOutputPrefixes = (renderRequest, grouping, recordCount) => {
  const outputs = grouping === 'batch' ? renderRequest.output : [renderRequest.output];
  if (
    !Array.isArray(outputs) ||
    outputs.length !== (grouping === 'batch' ? recordCount : 1) ||
    outputs.some((output) => !output || typeof output !== 'object' || Array.isArray(output))
  ) {
    throw new Error(
      grouping === 'batch'
        ? 'Canonical Circular batch output must contain one object per record.'
        : 'Canonical renderRequest output must be an object.'
    );
  }
  const prefixes = outputs.map((output) => String(output.prefix || '').trim());
  if (prefixes.some((prefix) => !prefix)) {
    throw new Error('Canonical renderRequest output prefixes must be non-empty.');
  }
  return prefixes;
};

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

export const projectCanonicalSessionRequest = ({
  renderRequest,
  resources: canonicalResources,
  webFiles = {},
  legacyFiles = null,
  storedConfig = null,
  initializeCliInputs = false,
  fileBindings = [],
  linearTrackSlotSchemaVersion = LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  repairInvalidComparisonHeight = false,
  sessionResourceTable = null,
  deferResourceContent = false,
  adoptCanonicalPayloads = false
}) => {
  if (!renderRequest || !SUPPORTED_CANONICAL_REQUEST_SCHEMAS.has(renderRequest.schema)) {
    throw new Error('Unsupported canonical renderRequest schema.');
  }
  if (!['circular', 'linear'].includes(renderRequest.mode)) {
    throw new Error('Unsupported canonical renderRequest mode.');
  }
  const sourceRecords = Array.isArray(renderRequest.records) ? renderRequest.records : [];
  if (renderRequest.schema >= 7) {
    sourceRecords.forEach((record) => canonicalRecordDisplay(record.display));
    if (!Array.isArray(renderRequest.diagramOptions?.featurePlacements)) {
      throw new Error('Canonical featurePlacements must be an array.');
    }
    canonicalFeaturePlacements(renderRequest.diagramOptions.featurePlacements, renderRequest.mode);
  } else if (sourceRecords.some((record) => Object.hasOwn(record, 'display'))
    || Object.hasOwn(renderRequest.diagramOptions || {}, 'featurePlacements')) {
    throw new Error('Record display and feature placements require canonical schema 7.');
  }
  const normalizedRecordOrdering = normalizeWebGridColumnOrdering(sourceRecords);
  const records = normalizedRecordOrdering.records;
  const reorderRecordIndexedValues = (values) => (
    normalizedRecordOrdering.sourceIndexByProjectedIndex
      .map((sourceIndex) => values?.[sourceIndex])
  );
  if (records.length === 0) throw new Error('Canonical renderRequest records are required.');
  let similarityAlignment = null;
  let recordTranslations = [];
  if (renderRequest.mode === 'linear' && renderRequest.schema >= 8) {
    const layout = renderRequest.layout || {};
    if (Object.keys(layout).length > 0) {
      requireExactCanonicalKeys(layout, [
        'recordGapPx', 'multiRecordPositions', 'recordTranslations', 'similarityAlignment'
      ], 'renderRequest.layout');
      if (!Number.isFinite(layout.recordGapPx) || layout.recordGapPx < 0) {
        throw new Error('renderRequest.layout.recordGapPx must be a finite non-negative number.');
      }
      if (layout.multiRecordPositions !== null && (
        !Array.isArray(layout.multiRecordPositions) ||
        layout.multiRecordPositions.some((token) => typeof token !== 'string' || !token.trim())
      )) throw new Error(
        'renderRequest.layout.multiRecordPositions must be null or an array of non-empty text.'
      );
      if (layout.multiRecordPositions?.length && sourceRecords.some(
        (record) => record.presentation?.gridRow != null
      )) throw new Error(
        'Linear row placement must use record presentation or layout positions, not both.'
      );
      const recordKeys = sourceRecords.map((record, index) => requireCanonicalText(
        record.recordKey,
        `renderRequest.records[${index}].recordKey`
      ));
      similarityAlignment = canonicalSimilarityAlignment(
        layout.similarityAlignment,
        recordKeys,
        'renderRequest.layout.similarityAlignment'
      );
      recordTranslations = canonicalRecordTranslations(
        layout.recordTranslations,
        recordKeys,
        'renderRequest.layout.recordTranslations',
        { requireCoverage: similarityAlignment !== null }
      );
    }
  } else if (renderRequest.mode === 'linear' && (
    Object.hasOwn(renderRequest.layout || {}, 'recordTranslations') ||
    Object.hasOwn(renderRequest.layout || {}, 'similarityAlignment')
  )) {
    throw new Error('Typed similarity alignment requires canonical schema 8.');
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
  let explicitBindings = validateWebFileBindings(webMetadata, canonicalResources);
  // A CLI sidecar has no saved Web draft. Its empty writer slots are initial
  // values, not a user's cleared inputs. Keep real original CLI bindings (in
  // particular GFF + FASTA), and initialize absent slots from the typed request.
  // A saved Web draft, including explicit null/[] inputs, always wins unchanged.
  if (initializeCliInputs && storedConfig == null && explicitBindings) {
    const originalInputFields = renderRequest.mode === 'circular'
      ? ['c_gb', 'c_gff', 'c_fasta']
      : ['linearSeqs'];
    const hasOriginalInputs = originalInputFields.some((field) => (
      Array.isArray(explicitBindings[field])
        ? explicitBindings[field].length > 0
        : explicitBindings[field] != null
    ));
    explicitBindings = Object.fromEntries(Object.entries(explicitBindings).filter(
      ([field, value]) => (hasOriginalInputs && originalInputFields.includes(field))
        || (value != null && (!Array.isArray(value) || value.length > 0))
    ));
    // A CLI binding uid (`cli-seq-N`) is only an initial value. The committed
    // request owns record identity, so each Linear file takes the file-level
    // recordKey that Inherit matches against.
    const fileRecordKeys = [...new Set(sourceRecords.map((record) => (
      String(record.recordKey || '').replace(/:[1-9]\d*$/, '')
    )))];
    if (renderRequest.mode === 'linear' && Array.isArray(explicitBindings.linearSeqs)
      && fileRecordKeys.every(Boolean)
      && fileRecordKeys.length === explicitBindings.linearSeqs.length) {
      explicitBindings = {
        ...explicitBindings,
        linearSeqs: explicitBindings.linearSeqs.map((sequence, index) => ({
          ...sequence,
          uid: fileRecordKeys[index]
        }))
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
          || normalizeOriginalResourceName(originalNameHints[resourceId])
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
  const savedLinearRecordMetadata = Array.isArray(webMetadata.linearRecordMetadata)
    ? webMetadata.linearRecordMetadata
    : [];
  const savedLinearRecordMetadataByKey = new Map(
    savedLinearRecordMetadata
      .map((entry) => [String(entry?.recordKey || ''), entry])
      .filter(([recordKey]) => recordKey)
  );
  const legacyLinearSequences = Array.isArray(legacyFiles?.linearSeqs)
    ? legacyFiles.linearSeqs
    : [];
  const projectedProteinPipeline = projectGeneratedProteinPipeline(
    (renderRequest.comparisons || []).find(
      (comparison) => comparison?.kind === 'generatedProteinComparison'
    ),
    { adoptCanonicalPayloads, requestSchema: renderRequest.schema }
  );
  const comparisonsContainGeneratedProteinPipeline = (
    renderRequest.comparisons || []
  ).some((comparison) => comparison?.kind === 'generatedProteinComparison');
  const files = { linearSeqs: [] };
  if (renderRequest.mode === 'circular') {
    files.circularRecords = records.map((record) => {
      const source = record.source || {};
      const resourceFile = (id) => resolveResourceFile
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
      const region = record.region || null;
      const selector = region?.selector || record.selector;
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
        region_record_id: selector?.kind === 'recordId' ? selector.value : (selector?.kind === 'recordIndex' ? `#${selector.index + 1}` : ''),
        region_start: region?.start ?? null,
        region_end: region?.end ?? null,
        region_reverse: Boolean(region?.reverseComplement || record.presentation?.reverseComplement)
      };
    });
    files.linearComparisons = [];
    files.linearCanonicalComparisons = [];
    (renderRequest.comparisons || [])
      .filter((comparison) => comparison?.kind === 'nucleotideBlast')
      .forEach((comparison, index) => {
        const sourceQueryIndex = Number.isInteger(Number(comparison.queryRecordIndex))
          ? Number(comparison.queryRecordIndex)
          : index;
        const sourceSubjectIndex = Number.isInteger(Number(comparison.subjectRecordIndex))
          ? Number(comparison.subjectRecordIndex)
          : index + 1;
        const queryIndex = normalizedRecordOrdering.projectedIndexBySourceIndex
          .get(sourceQueryIndex);
        const subjectIndex = normalizedRecordOrdering.projectedIndexBySourceIndex
          .get(sourceSubjectIndex);
        const file = resolveResourceFile
          ? resolveResourceFile(comparison.resourceId)
          : resourceAsLegacyFile(resources, comparison.resourceId);
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
    (renderRequest.comparisons || []).forEach((comparison) => {
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
          .map((pair) => ({
            queryRecordIndex: normalizedRecordOrdering.projectedIndexBySourceIndex
              .get(Number(pair?.queryRecordIndex)),
            subjectRecordIndex: normalizedRecordOrdering.projectedIndexBySourceIndex
              .get(Number(pair?.subjectRecordIndex))
          }));
        files.linearCanonicalComparisons.push(
          projectedComparison
        );
        (Array.isArray(comparison.pairs) ? comparison.pairs : [])
          .forEach((pair, index) => {
            const queryIndex = normalizedRecordOrdering.projectedIndexBySourceIndex
              .get(Number(pair?.queryRecordIndex));
            const subjectIndex = normalizedRecordOrdering.projectedIndexBySourceIndex
              .get(Number(pair?.subjectRecordIndex));
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
  validateCanonicalAssemblyOutput(options.output, renderRequest.schema);
  const canonicalDepth = projectCanonicalDepthTracks({
    options,
    records,
    resources,
    mode: renderRequest.mode,
    resolveResourceFile
  });
  const projectedCircularSizeMode = renderRequest.mode === 'circular'
    ? (
        renderRequest.schema <= 2
          ? migratePersistedCircularMultiRecordSizeMode(
              renderRequest.layout?.multiRecordSizeMode
            )
          : requireCurrentCircularMultiRecordSizeMode(
              renderRequest.layout?.multiRecordSizeMode
            )
      )
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
    if (depthRows.length !== records.length || depthRows.some((row) => !Array.isArray(row))) {
      throw new Error(
        `Canonical ${renderRequest.mode === 'circular' ? 'Circular' : 'Linear'} Depth matrix must contain one row per record (${records.length}).`
      );
    }
  }
  if (renderRequest.mode === 'circular' && depthRows.length > 0) {
    files.c_depth = normalizeRecordMajorDepthFileRows(
      depthFileRows || depthRows.map((row) => row.map((ref) => (
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
    depthRows.forEach((row, index) => {
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
      projectedDefaultColors = parseColorTable(
        resourceTextFromRef(resources, defaultColorsRef)
      ).colors;
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
      projectedSpecificRules = parseSpecificRules(
        resourceTextFromRef(resources, colorTableRef)
      ).rules.map(({ fromFile: _fromFile, ...rule }) => rule);
    }
  }
  let projectedWhitelist = deferResourceContent && Array.isArray(storedConfig?.whitelist)
    ? storedConfig.whitelist
    : [];
  if (options.labelWhitelistFile?.resourceId) {
    files.whitelist = resolveResourceFile
      ? resolveResourceFile(options.labelWhitelistFile.resourceId)
      : resourceAsLegacyFile(resources, options.labelWhitelistFile.resourceId);
    if (!deferResourceContent) {
      projectedWhitelist = parseWhitelistRules(
        resourceTextFromRef(resources, options.labelWhitelistFile)
      ).rules;
    }
  }
  const qualifierPriorityRef = options.qualifierPriorityFile || options.qualifierPriorityTable;
  let projectedPriorityRules = deferResourceContent
    && Array.isArray(storedConfig?.qualifierPriorityRules)
    ? storedConfig.qualifierPriorityRules
    : [];
  if (qualifierPriorityRef?.resourceId) {
    files.qualifier_priority = resolveResourceFile
      ? resolveResourceFile(qualifierPriorityRef.resourceId)
      : resourceAsLegacyFile(resources, qualifierPriorityRef.resourceId);
    if (!deferResourceContent) {
      projectedPriorityRules = parsePriorityRules(
        resourceTextFromRef(resources, qualifierPriorityRef)
      ).rules;
    }
  }
  const projectedFeatureVisibilityRules = !deferResourceContent
    && options.featureVisibilityTableFile?.resourceId
    ? parseFeatureVisibilityRules(
        resourceTextFromRef(resources, options.featureVisibilityTableFile)
      ).rules
    : [];
  const projectedLabelOverrideRows = !deferResourceContent
    && options.labelOverrideFile?.resourceId
    ? parseLabelOverrideTsv(resourceTextFromRef(resources, options.labelOverrideFile)).map((row) => ({
        recordId: row.recordId,
        featureType: row.featureType,
        qualifier: row.qualifier,
        valueRegex: row.valueRegex,
        labelText: row.labelText
      }))
    : [];
  if (renderRequest.mode === 'circular' && Array.isArray(options.conservationBlastFiles)) {
    files.c_conservation_blasts = options.conservationBlastFiles
      .map((ref) => ref?.resourceId
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
    files.c_conservation_sequence_sources = options.conservationFastaFiles.map((ref) => (
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
      (resourceId) => (resourceId
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
      .map((resourceId) => (
        resourceId
          ? (resolveResourceFile
              ? resolveResourceFile(resourceId)
              : resourceAsLegacyFile(resources, resourceId))
          : null
      ));
  }
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
  const legacySparseDefaults = renderRequest.schema <= 2;
  const currentTrackDefaults = trackDefaultsForMode(renderRequest.mode);
  const sparseConfigDefaults = legacySparseDefaults
    ? HISTORICAL_CONFIG_OVERRIDES[renderRequest.mode]
    : {
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
    ? (
        renderRequest.schema <= 2
          ? migratePersistedLinearLabelPlacement(overrides.label_placement)
          : requireCurrentLinearLabelPlacement(overrides.label_placement)
      )
    : 'auto';
  const projectedLinearTrackLayout = renderRequest.mode === 'linear'
    ? (
        renderRequest.schema <= 2
          ? migratePersistedLinearTrackLayout(overrides.linear_track_layout)
          : requireCurrentLinearTrackLayout(overrides.linear_track_layout)
      )
    : 'middle';
  const comparisonHeight = classifyOptionalPositiveNumber(overrides.comparison_height);
  if (renderRequest.mode === 'linear' && comparisonHeight.status === 'invalid') {
    if (!repairInvalidComparisonHeight) {
      throw diagnosticError('INPUT_INVALID', { field: 'match_height', reason: 'POSITIVE_OR_AUTO' });
    }
  }
  const tracks = options.tracks || {};
  const projectedCircularTrackSlots = renderRequest.mode === 'circular'
    ? (Array.isArray(tracks.circularTrackSlots)
        ? tracks.circularTrackSlots.map((slot, index) => (
          slot && typeof slot === 'object' && !Array.isArray(slot)
              ? (
                  renderRequest.schema <= 2
                    ? normalizeCircularTrackSlot(
                        projectLegacyCanonicalCircularSlot(slot),
                        index,
                        options.dinucleotide || 'GC',
                        overrides.track_type || 'tuckin'
                      )
                    : projectCurrentCanonicalCircularSlot(slot)
                )
              : parseCircularTrackSlotSpecs(
                  [
                    renderRequest.schema <= 2
                      ? migrateLegacyCircularTrackSlotSpec(slot)
                      : slot
                  ],
                  options.dinucleotide || 'GC',
                  overrides.track_type || 'tuckin'
                )[0]
          ))
        : [])
    : [];
  const projectedLinearTrackSlots = renderRequest.mode === 'linear'
    ? migrateLinearTrackSlotsToCurrentSchema(
        parseLinearTrackSlotSpecs(tracks.linearTrackSlots),
        linearTrackSlotSchemaVersion
      )
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
        ...depthRows.map((row) => Array.isArray(row) ? row.length : 0),
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
    options.depthTrackHeights.some((height) => height !== null && height !== undefined)
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
    labels_mode: renderRequest.mode === 'circular'
      ? ({ outer: 'out', both: 'both' }[overrides.label_scope] || 'none')
      : 'none',
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
    (
      renderRequest.schema <= 2 ||
      (
        Array.isArray(options.selectedFeaturesSet) &&
        options.selectedFeaturesSet.includes('repeat_region')
      )
    )
  ) {
    projectedFeatureShapes.repeat_region = renderRequest.schema <= 2
      ? 'rectangle'
      : defaultFeatureRendering('repeat_region');
  }
  const sparseFeatureTypes = legacySparseDefaults
    ? HISTORICAL_FEATURE_TYPES
    : MODE_DEFAULT_FEATURE_TYPES;
  const currentComparisonDefaults = comparisonFiltersForMode(renderRequest.mode);
  const sparseComparisonDefaults = legacySparseDefaults
    ? HISTORICAL_COMPARISON_DEFAULTS
    : {
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
    outer_label_x_offset: renderRequest.mode === 'circular'
      ? (overrides.outer_label_x_radius_offset ?? null)
      : null,
    outer_label_y_offset: renderRequest.mode === 'circular'
      ? (overrides.outer_label_y_radius_offset ?? null)
      : null,
    inner_label_x_offset: renderRequest.mode === 'circular'
      ? (overrides.inner_label_x_radius_offset ?? null)
      : null,
    inner_label_y_offset: renderRequest.mode === 'circular'
      ? (overrides.inner_label_y_radius_offset ?? null)
      : null,
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
    multi_record_positions: (renderRequest.layout?.multiRecordPositions || []).map((token) => {
      const split = String(token).lastIndexOf('@');
      return { selector: String(token).slice(0, split), row: Number(String(token).slice(split + 1)) };
    })
  };
  const presentationOwnsLinearRows = renderRequest.schema >= 6 && records.some(
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
          ...(initializeCliInputs && storedConfig == null
            && Number.isInteger(sourceColumn) && sourceColumn > 0
            ? { canonicalRow: row, canonicalColumn: sourceColumn }
            : {})
        };
      })
    : (renderRequest.layout?.multiRecordPositions || []).map((token, index) => {
        const split = String(token).lastIndexOf('@');
        return {
          uid: files.linearSeqs[index]?.uid || '',
          row: Number(String(token).slice(split + 1)) || index + 1
        };
      });
  const linearLayout = renderRequest.mode === 'linear' && renderRequest.schema >= 2
    ? {
        enabled: renderRequest.schema >= 6
          ? presentationOwnsLinearRows ||
            (renderRequest.layout?.multiRecordPositions || []).length > 0
          : Object.keys(renderRequest.layout || {}).length > 0,
        recordGap: renderRequest.layout?.recordGapPx ?? 24,
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
    layoutPreferences,
    config: {
      form,
      adv,
      ...(projectedProteinPipeline?.config || {}),
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
      annotationSets: normalizeAnnotationSets(options.annotations?.sets),
      recordDisplayDrafts: records.flatMap((record, index) => (record.display?.isCircular != null || record.display?.startCoordinate != null) ? [{
        scope: renderRequest.mode,
        sourceUid: renderRequest.mode === 'linear' ? String(files.linearSeqs[index]?.uid || record.recordKey) : 'circular',
        selector: record.selector?.kind === 'recordIndex' ? `#${record.selector.index + 1}` : '#1',
        recordId: record.selector?.kind === 'recordId' ? record.selector.value : '',
        topologyOverride: record.display?.isCircular ?? null,
        startCoordinate: record.display?.startCoordinate ?? null,
        reverseComplementOverride: null,
        anchorIntent: null
      }] : []),
      featurePlacementOverrides: Object.fromEntries(canonicalFeaturePlacements(
        options.featurePlacements || [], renderRequest.mode
      ).map((row) => [JSON.stringify([row.recordKey, row.biologicalFeatureId]), row])),
      circularConservation: renderRequest.mode === 'circular'
        ? projectCircularConservationConfig(options, files)
        : undefined
    },
    semanticFeatureState: {
      featureVisibilityManualRules: projectedFeatureVisibilityRules,
      featureVisibilityOverrides: {},
      labelOverrideRows: projectedLabelOverrideRows
    },
    pipelineState: projectedProteinPipeline
      ? {
          generatedProteinComparison:
            projectedProteinPipeline.generatedProteinComparison,
          legacySimilarityAlignment:
            projectedProteinPipeline.legacySimilarityAlignment
        }
      : null
  };
};
const PUBLICATION_OUTPUT_ONLY_FIELDS = new Set(['prefix', 'formats', 'overwrite', 'artifactFilename']);
const PUBLICATION_COMPARISON_FILTER_FIELDS = new Set(['evalue', 'bitscore', 'identity', 'alignmentLength']);
const PUBLICATION_OPTION_DEFAULTS = { dinucleotide: 'GC', keepFullDefinitionWithPlotTitle: false, conservationReference: 'auto', 'objects.features.arrow_geometry.head_length_ratio': 'auto', 'objects.features.arrow_geometry.shaft_width_ratio': 1, 'objects.scale.show': true };
const publicationBytes = async (resource, id) => {
  if (!resource || typeof resource !== 'object' || Array.isArray(resource)) throw new Error(`Canonical request resource '${id}' is missing.`);
  if (typeof resource.readBytes === 'function') return resource.readBytes();
  if (resource.data instanceof Uint8Array) return resource.data;
  if (resource.encoding === 'base64') { const binary = atob(String(resource.data || '')), bytes = new Uint8Array(binary.length); for (let index = 0; index < binary.length; index += 1) bytes[index] = binary.charCodeAt(index); return bytes; }
  if (typeof resource.data === 'string') return textToBytes(resource.data);
  throw new Error(`Canonical request resource '${id}' has no decodable payload.`);
};
const TSV_NUMBER = /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
const normalizedPublicationBytes = async (resource, id, normalize, path) => {
  const bytes = await publicationBytes(resource, id);
  if (path.includes('.diagramOptions.colors.defaultColors')) {
    const rows = bytesToText(bytes).split(/\r?\n/).map((row) => row.trim()).filter(
      (row) => row && row !== 'feature_type\tcolor').sort();
    return textToBytes(`${rows.join('\n')}\n`);
  }
  if (!normalize) return bytes;
  if (!path.startsWith('$.comparisons') || resource.kind !== 'canonical-tsv') return bytes;
  const rows = bytesToText(bytes).trimEnd().split(/\r?\n/).map((row, index) => index === 0
    ? row : row.split('\t').map((cell) => {
      const number = Number(cell);
      return TSV_NUMBER.test(cell) && Number.isFinite(number)
        ? String(Number(number.toPrecision(15))) : cell;
    }).join('\t'));
  return textToBytes(`${rows.join('\n')}\n`);
};
const publicationResourceIdentity = async (resources, id, normalize, path, cache) => {
  const resourceId = String(id || '').trim();
  if (!resourceId) throw new Error('Canonical request contains an empty resourceId.');
  const resource = resources?.[resourceId], key = resource?.encoding === 'base64' && !path.includes('.diagramOptions.colors.defaultColors') && (!normalize || !path.startsWith('$.comparisons') || resource.kind !== 'canonical-tsv') ? resource.data : null;
  const kind = path.includes('.diagramOptions.colors.defaultColors') ? 'default-colors'
    : (path.includes('.diagramOptions.colors.colorTable') ? 'color-table' : String(resource.kind || ''));
  let digest = key ? cache.get(key) : null;
  if (!digest) { digest = normalizedPublicationBytes(resource, resourceId, normalize, path).then((bytes) => sha256Hex(bytes)); if (key) cache.set(key, digest); }
  return { kind, decodedPayloadSha256: await digest };
};
const canonicalizePublicationValue = async (value, resources, context, path = '$') => {
  if (Array.isArray(value)) return Promise.all(value.map((entry, index) =>
    canonicalizePublicationValue(entry, resources, context, `${path}[${index}]`)));
  if (!value || typeof value !== 'object') return value;
  const output = {};
  const resourceReference = Object.hasOwn(value, 'resourceId');
  for (const key of Object.keys(value).sort()) {
    if (value[key] === null || Object.is(value[key], PUBLICATION_OPTION_DEFAULTS[key]) || key === 'plotTitleFontSize'
      || (path === '$.diagramOptions.configOverrides' && key === 'labels.filtering.blacklist_keywords'
        && Array.isArray(value[key]) && value[key].length === 0)
      || (path === '$.diagramOptions.configOverrides' && key === 'canvas.feature_overlap_tolerance_bp' && value[key] === 0)
      || (path === '$.diagramOptions.config.canvas' && key === 'feature_overlap_tolerance_bp' && value[key] === 0)
      || (context.ignoreComparisonFilters && path === '$.diagramOptions' && PUBLICATION_COMPARISON_FILTER_FIELDS.has(key))
      || (path === '$.output' && PUBLICATION_OUTPUT_ONLY_FIELDS.has(key))) continue;
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
const firstPublicationDiff = (expected, actual, path = '$') => {
  if (Object.is(expected, actual)) return null;
  if (!expected || !actual || typeof expected !== 'object' || typeof actual !== 'object'
      || Array.isArray(expected) !== Array.isArray(actual)) return { path, expected, actual };
  for (const key of new Set([...Object.keys(expected), ...Object.keys(actual)])) {
    const childPath = Array.isArray(expected) ? `${path}[${key}]` : `${path}.${key}`;
    if (!Object.hasOwn(expected, key) || !Object.hasOwn(actual, key)) return {
      path: childPath, expected: expected[key], actual: actual[key] };
    const difference = firstPublicationDiff(expected[key], actual[key], childPath);
    if (difference) return difference;
  }
  return null;
};
export const promoteCanonicalRenderRequestToCurrent = (
  request,
  { featureCatalog = null, legacyOrthogroupState = null } = {}
) => {
  const promoted = cloneCanonicalJsonValue(request);
  if (promoted.schema === CANONICAL_REQUEST_SCHEMA) return promoted;
  if (![5, 6, 7].includes(promoted.schema)) {
    throw new Error('Only canonical renderRequest schemas 5, 6, and 7 can be promoted to schema 8.');
  }
  const sourceSchema = promoted.schema;
  const linearRows = promoted.mode === 'linear'
    ? (promoted.layout?.multiRecordPositions || []).map((token) => {
        const split = String(token).lastIndexOf('@');
        return Number(String(token).slice(split + 1)) || null;
      })
    : [];
  promoted.schema = CANONICAL_REQUEST_SCHEMA;
  (promoted.records || []).forEach((record, index) => {
    if (sourceSchema < 7) record.display = { isCircular: null, startCoordinate: null };
    if (sourceSchema === 5) record.cardinality = promoted.mode === 'linear' &&
      !record.selector && !record.region
      ? 'all'
      : 'exactly_one';
    if (linearRows[index]) record.presentation.gridRow = linearRows[index];
  });
  if (sourceSchema < 7) {
    promoted.diagramOptions = { ...promoted.diagramOptions, featurePlacements: [] };
  }
  if (promoted.mode === 'linear') {
    const generated = (promoted.comparisons || []).find(
      (comparison) => comparison?.kind === 'generatedProteinComparison'
    );
    const rawLegacy = generated?.settings?.alignOrthogroupFeature;
    const similarityAlignment = rawLegacy === null || rawLegacy === undefined
      ? null
      : canonicalSimilarityAlignment(materializeLegacySimilarityAlignment({
          target: rawLegacy,
          records: promoted.records || [],
          featureCatalog,
          legacyOrthogroupState
        }), (promoted.records || []).map((record) => record.recordKey),
        'renderRequest.layout.similarityAlignment');
    if (generated?.settings) delete generated.settings.alignOrthogroupFeature;
    const hadLayout = Object.keys(promoted.layout || {}).length > 0;
    if (hadLayout || similarityAlignment !== null) {
      const recordGapPx = Number(promoted.layout?.recordGapPx ?? 24);
      promoted.layout = {
        recordGapPx: Number.isFinite(recordGapPx) && recordGapPx >= 0 ? recordGapPx : 24,
        multiRecordPositions: null,
        recordTranslations: (promoted.records || []).map((record) => ({
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

const sameCanonicalValue = (left, right) => (
  JSON.stringify(left) === JSON.stringify(right)
);

const materializedRecordSelector = (selector) => {
  const match = String(selector || '').match(/^#([1-9]\d*)$/);
  if (!match) {
    throw new Error('Target record materialization requires an exact record selector.');
  }
  return { kind: 'recordIndex', index: Number(match[1]) - 1 };
};

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
  if (!Number.isSafeInteger(transform?.startCoordinate)
    || transform.startCoordinate < 1 || transform.startCoordinate > recordLength
    || typeof transform?.reverseComplement !== 'boolean') {
    throw new Error('Resolved record transform is invalid for the target record.');
  }
};

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

const shiftCanonicalComparisonIndexes = (comparisons, recordIndex, expansion) => (
  (Array.isArray(comparisons) ? comparisons : []).map((comparison) => {
    const shifted = cloneCanonicalJsonValue(comparison);
    for (const field of ['queryRecordIndex', 'subjectRecordIndex']) {
      const index = Number(shifted?.[field]);
      if (Number.isInteger(index) && index > recordIndex) shifted[field] = index + expansion;
    }
    if (Array.isArray(shifted?.pairs)) {
      shifted.pairs = shifted.pairs.map((pair) => {
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
 * Clone the last committed canonical Session and overlay one record transform.
 * No live form state participates in this projection.
 */
export const projectCommittedRecordTransform = ({ committed, target, transform }) => {
  const request = committed?.renderRequest;
  if (!request || request.schema !== CANONICAL_REQUEST_SCHEMA
    || !['circular', 'linear'].includes(request.mode)
    || !committed?.resources || typeof committed.resources !== 'object') {
    throw new Error('A current committed canonical Session is required.');
  }
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

/** Project only alignment-owned fields of the last committed canonical artifact. */
export const projectCommittedSimilarityAlignment = ({ committed, plan, translations, orientations }) => {
  if (committed?.renderRequest?.schema !== CANONICAL_REQUEST_SCHEMA
    || committed.renderRequest.mode !== 'linear' || !committed.resources) {
    throw new Error('Alignment requires a committed canonical Linear artifact.');
  }
  const canonical = { ...committed, renderRequest: cloneCanonicalJsonValue(committed.renderRequest) };
  const request = canonical.renderRequest;
  const keys = request.records.map(({ recordKey }) => recordKey);
  const byKey = new Map((orientations || []).map((entry) => [entry.recordKey, entry.reverseComplement]));
  if (byKey.size !== keys.length || orientations.length !== keys.length
    || keys.some(key => typeof byKey.get(key) !== 'boolean')) {
    throw new Error('Alignment orientation coverage changed.');
  }
  request.records.forEach(record => writeCanonicalRecordReverseComplement(record, byKey.get(record.recordKey)));
  request.layout = { recordGapPx: 24, multiRecordPositions: null, ...request.layout,
    similarityAlignment: cloneCanonicalJsonValue(plan ?? null),
    recordTranslations: cloneCanonicalJsonValue(translations) };
  projectCanonicalSessionRequest({ ...canonical, deferResourceContent: true });
  return canonical;
};

const normalizePublicationRequestAliases = (request) => {
  const normalized = [5, 6, 7].includes(request?.schema)
    ? promoteCanonicalRenderRequestToCurrent(request)
    : cloneCanonicalJsonValue(request);
  for (const comparison of normalized.comparisons || []) {
    if (comparison.kind === 'generatedProteinComparison' && comparison.settings) {
      comparison.settings.collinearInferOrthogroups = resolvePipelineCollinearInference(
        comparison.settings, comparison.mode
      );
    }
  }
  const colors = normalized.diagramOptions?.colors;
  if (colors) {
    colors.defaultColors = colors.defaultColors || colors.defaultColorsFile || null;
    colors.colorTable = colors.colorTable || colors.colorTableFile || null;
    delete colors.defaultColorsFile;
    delete colors.colorTableFile;
  }
  return normalized;
};
const publicationRequestIdentity = async (request, resources, normalize, resourceIdentities) => {
  if (!request || typeof request !== 'object' || Array.isArray(request)) throw new Error('Canonical request equivalence requires a renderRequest object.');
  const normalizedRequest = normalizePublicationRequestAliases(request);
  const conservation = normalizedRequest.diagramOptions?.conservationBlastFiles;
  const context = { bindings: [], normalize, resourceIdentities,
    ignoreComparisonFilters: !(normalizedRequest.comparisons?.length || (Array.isArray(conservation) ? conservation.length : conservation)) };
  const canonical = await canonicalizePublicationValue(normalizedRequest, resources, context);
  return { canonical, digest: await sha256Hex(textToBytes(JSON.stringify(canonical))),
    resourceBindings: context.bindings };
};
export const compareCanonicalRenderRequests = async (input) => {
  const normalize = input.normalizeReplayGeneratedResources === true;
  const resourceIdentities = new Map(), [expected, actual] = await Promise.all([publicationRequestIdentity(input.expectedRequest, input.expectedResources, normalize, resourceIdentities), publicationRequestIdentity(input.actualRequest, input.actualResources, normalize, resourceIdentities)]);
  const difference = firstPublicationDiff(expected.canonical, actual.canonical);
  return { equivalent: expected.digest === actual.digest, expected, actual, differences: difference ? [difference] : [] };
};
export const assertCanonicalRenderRequestsEquivalent = async (input) => {
  const comparison = await compareCanonicalRenderRequests(input);
  if (comparison.equivalent) return comparison;
  const error = new Error(`Gallery publication request differs at ${comparison.differences[0]?.path || '$'} (committed ${comparison.expected.digest}, rebuilt ${comparison.actual.digest}).`);
  error.comparison = comparison;
  throw error;
};
