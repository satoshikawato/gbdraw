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
  resolveDisambiguatedRecordSelection
} from '../app/record-options.js';
import {
  orderedConservationSources,
  orderedOptionalConservationFiles
} from '../app/conservation-series.js';
import {
  assertValidCustomTrackPlan,
  validateCustomTrackPlan,
  validateTrackSlotBindingInvariants
} from '../app/track-slot-validation.js';
import { annotationOptionsPayload, normalizeAnnotationSets } from '../app/annotations/state.js';
import { classifyOptionalPositiveNumber } from '../utils/optional-positive-number.js';
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
  trackDefaultsForMode
} from '../mode-profiles.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';
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
  textToBase64,
  textToBytes
} from './file-content-cache.js';
import {
  adoptedSessionResourceDescriptor,
  createCombinedSessionResourceFileView,
  createSessionResourceFileView
} from './session-resource-backing.js';
import { normalizeLinearComparisonPlan } from '../app/linear-comparisons.js';
import {
  getResourcePayloadOwner,
  setResourcePayloadOwner
} from './resource-payload-owner.js';
import { sha256Hex } from './byte-utils.js';

export const CANONICAL_REQUEST_SCHEMA = 6;
const SUPPORTED_CANONICAL_REQUEST_SCHEMAS = new Set([
  1, 2, 5, CANONICAL_REQUEST_SCHEMA
]);

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

const optionalNumber = (value) => {
  if (value === null || value === undefined || String(value).trim() === '') return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) ? numeric : null;
};

const optionalPositiveInteger = (value) => {
  const numeric = optionalNumber(value);
  return Number.isInteger(numeric) && numeric > 0 ? numeric : null;
};

const canonicalOptionalPositiveNumber = (value, fieldName) => {
  if (value === null || value === undefined || String(value).trim() === '') return null;
  const numeric = Number(value);
  if (!Number.isFinite(numeric) || numeric <= 0) {
    throw new Error(`${fieldName} must be null or a positive finite number.`);
  }
  return numeric;
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

const createResourceBuilder = () => {
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
    const existingResourceId = kindResourceIds.get(entry);
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
    kindResourceIds.set(entry, resourceId);
    fileResourceIds.set(effectiveKind, kindResourceIds);
    const originalName = normalizeOriginalResourceName(entry.name);
    if (originalName) resourceOriginalNames[resourceId] = originalName;
    return resourceId;
  };

  const addText = (resourceId, kind, name, text) => {
    if (resources[resourceId]) return resourceId;
    const normalized = String(text || '');
    const bytes = textToBytes(normalized);
    resources[resourceId] = {
      kind,
      name: normalizeResourceName(resourceId, name),
      type: 'text/tab-separated-values',
      size: bytes.byteLength,
      lastModified: 0,
      encoding: 'base64',
      data: textToBase64(normalized)
    };
    return resourceId;
  };

  const addJson = (resourceId, kind, name, value) => {
    if (resources[resourceId]) return resourceId;
    const normalized = JSON.stringify(value);
    const bytes = textToBytes(normalized);
    resources[resourceId] = {
      kind,
      name: normalizeResourceName(resourceId, name),
      type: 'application/json',
      size: bytes.byteLength,
      lastModified: 0,
      encoding: 'base64',
      data: textToBase64(normalized)
    };
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
  gridRow = null
} = {}) => ({
  label: String(label || '').trim() || null,
  subtitle: String(subtitle || '').trim() || null,
  reverseComplement: Boolean(reverseComplement),
  gridRow,
  gridColumn: null
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
  const recordId = safePrefix(record?.recordId, 'record');
  const selector = safePrefix(record?.selector, '1');
  return `circular-${recordId}-${selector}`;
};

const linearRegionPayload = (seq) => {
  const start = optionalPositiveInteger(seq?.region_start);
  const end = optionalPositiveInteger(seq?.region_end);
  if (start === null && end === null) return null;
  if (start === null || end === null) {
    throw new Error('Canonical linear regions require both start and end coordinates.');
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
    const resolvedRows = state.linearRecordLayoutEnabled?.value
      ? resolvedLinearRecordRows(filesData.linearSeqs, state.linearRecordRows)
      : [];
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
        cardinality: seq.cardinality || (selector || region ? 'exactly_one' : 'all'),
        source,
        selector,
        region,
        presentation: {
          ...presentationPayload({
            label: seq.definition,
            subtitle: seq.record_subtitle,
            gridRow: resolvedRows[index] ?? null
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
        presentation: publicationClone(record.presentation) || presentationPayload()
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
        const reason = selection.status === 'ambiguous' ? 'is ambiguous' : 'was not found';
        const label = requestedSelector || '(automatic)';
        throw new Error(`Circular record selector '${label}' ${reason} in the current input.`);
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
  const knownRecords = Array.isArray(state.circularRecordList.value)
    ? state.circularRecordList.value
    : [];
  const recordSelectors = buildDisambiguatedRecordEntries(
    knownRecords.map((record) => ({
      ...record,
      recordId: record?.record_id ?? record?.recordId
    }))
  );
  const requestedSelector = String(state.form.circular_record_selector || '').trim();
  const selection = resolveDisambiguatedRecordSelection(
    recordSelectors,
    requestedSelector
  );
  const singlePresentationRequested = (
    !state.form.multi_record_canvas &&
    state.adv.circular_grouping_intent !== 'batch'
  );
  if (
    singlePresentationRequested &&
    requestedSelector &&
    selection.status !== 'resolved'
  ) {
    const reason = selection.status === 'ambiguous' ? 'is ambiguous' : 'was not found';
    throw new Error(`Circular record selector '${requestedSelector}' ${reason} in the current input.`);
  }
  const selectedRecords = singlePresentationRequested && selection.record
    ? [selection.record]
    : (recordSelectors.length > 0 ? recordSelectors : [null]);
  const singleJourney = (
    selectedRecords.length === 1 &&
    singlePresentationRequested
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
    hasComparisonIntent = false
  } = {}
) => {
  const { form, adv } = state;
  const circular = state.mode.value === 'circular';
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
    throw new Error('Pairwise Match Height must be Auto or a positive finite number.');
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
    [CONFIG_OVERRIDE_PATHS.gcContentMode]: adv.gc_content_mode || 'deviation',
    [CONFIG_OVERRIDE_PATHS.gcContentMinPercent]: optionalNumber(adv.gc_content_min_percent),
    [CONFIG_OVERRIDE_PATHS.gcContentMaxPercent]: optionalNumber(adv.gc_content_max_percent),
    [CONFIG_OVERRIDE_PATHS.gcContentShowAxis]: Boolean(adv.gc_content_show_axis),
    [CONFIG_OVERRIDE_PATHS.gcContentShowTicks]: Boolean(adv.gc_content_show_ticks),
    [CONFIG_OVERRIDE_PATHS.gcContentLargeTickInterval]:
      optionalNumber(adv.gc_content_tick_interval),
    [CONFIG_OVERRIDE_PATHS.gcContentSmallTickInterval]:
      optionalNumber(adv.gc_content_small_tick_interval),
    [CONFIG_OVERRIDE_PATHS.gcContentTickFontSize]:
      optionalNumber(adv.gc_content_tick_font_size),
    [CONFIG_OVERRIDE_PATHS.depthColor]: adv.depth_color || null,
    [CONFIG_OVERRIDE_PATHS.depthMin]: optionalNumber(adv.depth_min),
    [CONFIG_OVERRIDE_PATHS.depthMax]: optionalNumber(adv.depth_max),
    [CONFIG_OVERRIDE_PATHS.depthNormalize]: Boolean(adv.depth_normalize),
    [CONFIG_OVERRIDE_PATHS.depthShowAxis]: Boolean(adv.depth_show_axis),
    [CONFIG_OVERRIDE_PATHS.depthShowTicks]: Boolean(adv.depth_show_ticks),
    [CONFIG_OVERRIDE_PATHS.depthLargeTickInterval]:
      optionalNumber(adv.depth_large_tick_interval),
    [CONFIG_OVERRIDE_PATHS.depthSmallTickInterval]:
      optionalNumber(adv.depth_small_tick_interval),
    [CONFIG_OVERRIDE_PATHS.depthTickFontSize]: optionalNumber(adv.depth_tick_font_size),
    [CONFIG_OVERRIDE_PATHS.depthShareAxis]: Boolean(adv.depth_share_axis),
    [CONFIG_OVERRIDE_PATHS.showScale]: form.show_scale !== false,
    [CONFIG_OVERRIDE_PATHS.scaleInterval]: optionalNumber(adv.scale_interval),
    [CONFIG_OVERRIDE_PATHS.labelBlacklist]: state.filterMode.value === 'Blacklist'
      ? String(state.manualBlacklist.value || '').split(/[,\n]/)
        .map((keyword) => keyword.trim()).filter(Boolean)
      : [],
    ...(circular
      ? {
          [CONFIG_OVERRIDE_PATHS.circularAxisStrokeColor]:
            adv.axis_stroke_color || null,
          [CONFIG_OVERRIDE_PATHS.circularDefinitionFontSize]:
            optionalNumber(adv.def_font_size),
          [CONFIG_OVERRIDE_PATHS.circularDefinitionInterval]: optionalNumber(adv.circular_definition_interval),
          [CONFIG_OVERRIDE_PATHS.plotTitleFontSize]:
            optionalNumber(adv.plot_title_font_size),
          [CONFIG_OVERRIDE_PATHS.circularLabelSpacing]:
            optionalNumber(adv.circular_label_spacing),
          [CONFIG_OVERRIDE_PATHS.circularLabelPlacement]:
            adv.circular_label_placement || 'horizontal',
          [CONFIG_OVERRIDE_PATHS.trackType]: form.track_type,
          [CONFIG_OVERRIDE_PATHS.tickLabelFontSize]:
            optionalNumber(adv.tick_label_font_size),
          [CONFIG_OVERRIDE_PATHS.outerLabelXRadiusOffset]:
            optionalNumber(adv.outer_label_x_offset),
          [CONFIG_OVERRIDE_PATHS.outerLabelYRadiusOffset]:
            optionalNumber(adv.outer_label_y_offset),
          [CONFIG_OVERRIDE_PATHS.innerLabelXRadiusOffset]:
            optionalNumber(adv.inner_label_x_offset),
          [CONFIG_OVERRIDE_PATHS.innerLabelYRadiusOffset]:
            optionalNumber(adv.inner_label_y_offset)
        }
      : {
          [CONFIG_OVERRIDE_PATHS.linearAxisStrokeColor]: linearAxisStrokeColor,
          [CONFIG_OVERRIDE_PATHS.linearDefinitionShowReplicon]:
            Boolean(adv.linear_show_replicon),
          [CONFIG_OVERRIDE_PATHS.linearDefinitionShowAccession]:
            Boolean(adv.linear_show_accession),
          [CONFIG_OVERRIDE_PATHS.linearDefinitionShowLength]:
            Boolean(adv.linear_show_length),
          [CONFIG_OVERRIDE_PATHS.linearLabelSpacing]:
            optionalNumber(adv.linear_label_spacing),
          [CONFIG_OVERRIDE_PATHS.labelPlacement]: linearLabelPlacement,
          [CONFIG_OVERRIDE_PATHS.labelRotation]: optionalNumber(adv.label_rotation),
          [CONFIG_OVERRIDE_PATHS.alignCenter]: Boolean(form.align_center),
          [CONFIG_OVERRIDE_PATHS.keepDefinitionLeftAligned]:
            Boolean(form.keep_definition_left_aligned),
          [CONFIG_OVERRIDE_PATHS.linearTrackLayout]: linearTrackLayout,
          [CONFIG_OVERRIDE_PATHS.linearTrackAxisGap]: optionalNumber(adv.track_axis_gap),
          [CONFIG_OVERRIDE_PATHS.linearRulerOnAxis]: Boolean(form.linear_ruler_on_axis),
          ...(hasComparisonIntent
            ? {
                [CONFIG_OVERRIDE_PATHS.comparisonHeight]:
                  comparisonHeight.status === 'auto' ? null : comparisonHeight.value,
                [CONFIG_OVERRIDE_PATHS.pairwiseMatchStyle]: adv.pairwise_match_style
              }
            : {}),
          [CONFIG_OVERRIDE_PATHS.gcHeight]: optionalNumber(adv.gc_height),
          [CONFIG_OVERRIDE_PATHS.depthHeight]: optionalNumber(adv.depth_height),
          [CONFIG_OVERRIDE_PATHS.scaleStyle]: form.scale_style,
          [CONFIG_OVERRIDE_PATHS.scaleStrokeColor]: adv.scale_stroke_color || null,
          [CONFIG_OVERRIDE_PATHS.scaleLabelColor]: adv.ruler_label_color || null,
          [CONFIG_OVERRIDE_PATHS.scaleStrokeWidth]:
            optionalNumber(adv.scale_stroke_width),
          [CONFIG_OVERRIDE_PATHS.normalizeLength]: Boolean(form.normalize_length)
        })
  };
  const sharedLengthValues = {
    [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.blockStrokeWidth]:
      optionalNumber(adv.block_stroke_width),
    [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.lineStrokeWidth]:
      optionalNumber(adv.line_stroke_width),
    [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.legendBoxSize]:
      optionalNumber(adv.legend_box_size),
    [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.legendFontSize]:
      optionalNumber(adv.legend_font_size),
    ...(circular
      ? {
          [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.circularAxisStrokeWidth]:
            optionalNumber(adv.axis_stroke_width),
          'labels.font_size': optionalNumber(adv.label_font_size)
        }
      : {
          [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.linearAxisStrokeWidth]:
            optionalNumber(adv.axis_stroke_width),
          [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.linearDefinitionFontSize]:
            optionalNumber(adv.def_font_size),
          [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.defaultCdsHeight]:
            optionalNumber(adv.feature_height),
          [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.scaleFontSize]:
            optionalNumber(adv.scale_font_size),
          [SHARED_LENGTH_CONFIG_OVERRIDE_PATHS.rulerLabelFontSize]:
            optionalNumber(adv.ruler_label_font_size),
          'labels.font_size.linear': optionalNumber(adv.label_font_size)
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
  return {
    ...(
      preservedOverrides && typeof preservedOverrides === 'object' && !Array.isArray(preservedOverrides)
        ? preservedOverrides
        : {}
    ),
    ...managedOverrides
  };
};

const addGeneratedTableResources = (state, resources, diagramOptions) => {
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
  const labelOverride = buildLabelOverrideTsv(
    state.labelTextFeatureOverrides,
    state.labelTextBulkOverrides,
    {
      editableLabels: state.editableLabels?.value || [],
      extractedFeatures: state.extractedFeatures.value,
      featureOverrideSources: state.labelTextFeatureOverrideSources,
      visibilityOverrides: state.labelVisibilityOverrides
    }
  );
  const labelOverrideTsv = serializeLabelOverrideRows(state.canonicalLabelOverrideRows?.value) || labelOverride.tsv;
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
        ? canonicalOptionalPositiveNumber(
            track.height,
            `depthTracks[${trackIndex}].height`
          )
        : null,
      largeTickInterval: canonicalOptionalPositiveNumber(
        track.large_tick_interval,
        `depthTracks[${trackIndex}].largeTickInterval`
      ),
      smallTickInterval: canonicalOptionalPositiveNumber(
        track.small_tick_interval,
        `depthTracks[${trackIndex}].smallTickInterval`
      ),
      tickFontSize: canonicalOptionalPositiveNumber(
        track.tick_font_size,
        `depthTracks[${trackIndex}].tickFontSize`
      )
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
        centerReservedRadius: optionalNumber(state.adv.center_reserved_radius)
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
          centerReservedRadius: optionalNumber(state.adv.center_reserved_radius)
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
      centerReservedRadius: circular ? optionalNumber(state.adv.center_reserved_radius) : null
    }
  };
};

const generatedProteinSettings = (state, baseline = {}) => {
  const blastp = state.losat.blastp || {};
  const blastpMode = requireCurrentProteinBlastpMode(blastp.mode);
  const positiveInteger = (value, fallback) => optionalPositiveInteger(value) ?? fallback;
  const nonNegativeInteger = (value, fallback) => {
    const numeric = optionalNumber(value);
    return Number.isInteger(numeric) && numeric >= 0 ? numeric : fallback;
  };
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
    ...baseline,
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
    collinearityAnchorMode: blastpMode === 'collinear'
      ? requireCurrentCollinearAnchorMode(blastp.collinearAnchorMode)
      : normalizeCollinearAnchorMode(blastp.collinearAnchorMode),
    collinearitySearchScope: blastpMode === 'collinear'
      ? requireCurrentCollinearSearchScope(blastp.collinearSearchScope)
      : normalizeCollinearSearchScope(blastp.collinearSearchScope),
    collinearityColorMode,
    losatpBin: baseline.losatpBin || 'losat',
    ncbiBlastpBin: baseline.ncbiBlastpBin ?? null,
    losatpThreads: optionalPositiveInteger(state.losat.threadsPerJob),
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
      : positiveInteger(blastp.orthogroupMemberMaxHits, 5),
    collinearMaxParalogLinksPerOrthogroup:
      blastpMode === 'collinear'
        ? requireCurrentCollinearMaxParalogLinks(
            blastp.collinearMaxParalogLinksPerOrthogroup
          )
        : positiveInteger(blastp.collinearMaxParalogLinksPerOrthogroup, 2),
    alignOrthogroupFeature:
      String(state.selectedOrthogroupAlignmentFeature.value || '').trim() || null
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
      || ![1, 2].includes(typedResource.schema)
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

const resolvedLinearRecordRows = (sequences, layoutRows) => {
  const rowsByUid = new Map(
    (Array.isArray(layoutRows) ? layoutRows : [])
      .map((entry) => [String(entry?.uid || ''), Number(entry?.row)])
  );
  return (Array.isArray(sequences) ? sequences : []).map((sequence, index) => {
    const row = rowsByUid.get(String(sequence?.uid || ''));
    return Number.isInteger(row) && row > 0 ? row : index + 1;
  });
};

export const linearRecordLayoutHasSharedRow = (sequences, layoutRows) => {
  const seenRows = new Set();
  return resolvedLinearRecordRows(sequences, layoutRows).some((row) => {
    if (seenRows.has(row)) return true;
    seenRows.add(row);
    return false;
  });
};

const buildLayout = (state, filesData) => {
  if (state.mode.value === 'linear') {
    if (!state.linearRecordLayoutEnabled?.value) return {};
    return {
      recordGapPx: Math.max(0, Number(state.linearRecordGap?.value) || 0)
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
    multiRecordMinRadiusRatio: Number(state.adv.multi_record_min_radius_ratio) || 0.55,
    multiRecordColumnGapRatio: Number(state.adv.multi_record_column_gap_ratio) || 0,
    multiRecordRowGapRatio: Number(state.adv.multi_record_row_gap_ratio) || 0,
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
    selectedOrthogroupAlignmentFeature: session?.orthogroupState?.selectedOrthogroupAlignmentFeature || '',
    linearRecordLayoutEnabled: Boolean(layout.enabled), linearRecordGap: layout.recordGap ?? 24
  };
  return {
    ...Object.fromEntries(Object.entries(refs).map(([key, value]) => [key, publicationRef(value)])),
    form: config.form || {}, adv: config.adv || {}, normalizePaletteColors,
    manualSpecificRules: publicationClone(config.rules || []), manualWhitelist: publicationClone(config.whitelist || []), manualPriorityRules: publicationClone(config.qualifierPriorityRules || []),
    labelTextFeatureOverrides: publicationClone(features.labelTextFeatureOverrides || {}), labelTextBulkOverrides: publicationClone(features.labelTextBulkOverrides || {}),
    labelTextFeatureOverrideSources: publicationClone(features.labelTextFeatureOverrideSources || {}), labelVisibilityOverrides: publicationClone(features.labelVisibilityOverrides || {}),
    circularConservation: conservation, losat: publicationClone(config.losat || { blastp: {} }),
    linearRecordRows: publicationClone(layout.rows || []), linearComparisonPlan: normalizeLinearComparisonPlan(config.linearComparisonPlan || { mode: 'none', defaultSource: 'losat', edges: [] }),
    annotationSets: publicationClone(config.annotationSets || []), canonicalPublicationFiles
  };
};
export const buildCanonicalRenderRequest = ({
  state,
  filesData,
  comparisonPlanSnapshot = null,
  resolvedComparisons = [],
  resolvedCircularConservation = []
}) => {
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
  const resources = createResourceBuilder();
  const webFiles = {};
  const recordPlan = buildRecords({ state, filesData, resources });
  const records = recordPlan.records;
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
    filesData: selectedCircularFilesData,
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
    configOverrides: buildConfigOverrides(state, {
      depthRequested: trackPlan.depthRequested,
      hasComparisonIntent: hasLinearComparisonIntent
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
    window: optionalPositiveInteger(state.adv.window_size),
    step: optionalPositiveInteger(state.adv.step_size),
    depthWindow: optionalPositiveInteger(state.adv.depth_window_size),
    depthStep: optionalPositiveInteger(state.adv.depth_step_size),
    plotTitle: String(state.form.plot_title || '').trim() || null,
    plotTitleFontSize: optionalNumber(state.adv.plot_title_font_size),
    ...(comparisonOptionsRequested
      ? {
          evalue: Number(state.adv.evalue),
          bitscore: Number(state.adv.min_bitscore),
          identity: Number(state.adv.identity),
          alignmentLength: Number(state.adv.alignment_length) || 0
        }
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
      diagramOptions.conservationReference = String(state.circularConservation.reference || 'auto');
      diagramOptions.conservationLabels = conservationEntries.map((entry) => entry.label);
      diagramOptions.conservationColors = conservationEntries.map((entry) => entry.color);
      diagramOptions.conservationRingWidth = optionalNumber(state.circularConservation.ring_width);
      diagramOptions.conservationRingGap = optionalNumber(state.circularConservation.ring_gap);
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
      diagramOptions.conservationReference = String(
        state.circularConservation.reference || 'subject'
      );
      diagramOptions.conservationLabels = resolvedCircularConservation.map(
        (entry, index) => String(entry?.label || `Comparison ${index + 1}`)
      );
      diagramOptions.conservationColors = resolvedCircularConservation.map(
        (entry) => String(entry?.color || '#D9EAF7')
      );
      diagramOptions.conservationRingWidth = optionalNumber(
        state.circularConservation.ring_width
      );
      diagramOptions.conservationRingGap = optionalNumber(
        state.circularConservation.ring_gap
      );
      webFiles.conservationBlastSource = 'losat-cache';
    }
  } else if (hasLinearComparisonIntent) {
    diagramOptions.pairwiseMatchStyle = String(state.adv.pairwise_match_style || 'ribbon');
  }

  addGeneratedTableResources(state, resources, diagramOptions);
  if (trackPlan.depthRequested) {
    buildDepthResources({
      state,
      filesData: selectedCircularFilesData,
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
    webFiles.linearRecordMetadata = (filesData.linearSeqs || []).map((sequence, index) => ({
      recordKey: String(records[index]?.recordKey || sequence?.uid || `record-${index + 1}`),
      losatGencode: optionalPositiveInteger(sequence?.losat_gencode) || 1
    }));
  }

  return {
    renderRequest: {
      schema: CANONICAL_REQUEST_SCHEMA,
      mode: state.mode.value,
      grouping,
      records,
      diagramOptions,
      layout: buildLayout(state, filesData),
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

const webBindingAsLegacyFile = (resources, binding, resolveResourceFile = null) => {
  if (binding === null || binding === undefined) return null;
  if (!binding || typeof binding !== 'object' || Array.isArray(binding)) {
    throw new Error('Web file bindings must be objects or null.');
  }
  const resourceId = String(binding.resourceId || '').trim();
  if (!resourceId) throw new Error('A Web file binding requires a resourceId.');
  const metadata = {
    name: normalizeOriginalResourceName(binding.name),
    type: String(binding.type || ''),
    lastModified: Number(binding.lastModified) || 0
  };
  if (resolveResourceFile) return resolveResourceFile(resourceId, metadata);
  const file = resourceAsLegacyFile(resources, resourceId);
  return { ...file, ...metadata, name: metadata.name || file.name };
};

const webBindingValueAsLegacyFile = (resources, value, resolveResourceFile = null) => (
  Array.isArray(value)
    ? value.map((item) => webBindingValueAsLegacyFile(resources, item, resolveResourceFile))
    : webBindingAsLegacyFile(resources, value, resolveResourceFile)
);

const applyWebFileBindings = (
  files,
  webMetadata,
  resources,
  { resolveResourceFile = null, adoptCanonicalPayloads = false } = {}
) => {
  const bindings = webMetadata?.bindings;
  if (bindings === undefined || bindings === null) return files;
  if (!bindings || typeof bindings !== 'object' || Array.isArray(bindings)) {
    throw new Error('Session webFiles.bindings must be an object.');
  }
  if (bindings.schema !== 1) {
    throw new Error('Unsupported Web file binding schema.');
  }

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
    restored[field] = webBindingValueAsLegacyFile(
      resources,
      bindings[field],
      resolveResourceFile
    );
  });
  restored.c_conservation_blasts_source =
    bindings.c_conservation_blasts_source === 'losat-cache' ? 'losat-cache' : null;

  if (Array.isArray(bindings.linearSeqs)) {
    restored.linearSeqs = bindings.linearSeqs.map((sequence, index) => ({
      uid: String(sequence?.uid || `canonical-seq-${index + 1}`),
      gb: webBindingValueAsLegacyFile(resources, sequence?.gb, resolveResourceFile),
      gff: webBindingValueAsLegacyFile(resources, sequence?.gff, resolveResourceFile),
      fasta: webBindingValueAsLegacyFile(resources, sequence?.fasta, resolveResourceFile),
      depth: webBindingValueAsLegacyFile(resources, sequence?.depth, resolveResourceFile),
      blast: webBindingValueAsLegacyFile(resources, sequence?.blast, resolveResourceFile),
      losat_gencode: optionalPositiveInteger(sequence?.losat_gencode) || 1,
      losat_filename: String(sequence?.losat_filename || ''),
      definition: String(sequence?.definition || ''),
      record_subtitle: String(sequence?.record_subtitle || ''),
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
      file: webBindingValueAsLegacyFile(resources, comparison?.file, resolveResourceFile)
    }));
  }
  if (Array.isArray(bindings.linearCanonicalComparisons)) {
    restored.linearCanonicalComparisons = bindings.linearCanonicalComparisons.map((comparison) => ({
      ...(adoptCanonicalPayloads ? comparison : cloneCanonicalJsonValue(comparison)),
      file: webBindingValueAsLegacyFile(resources, comparison?.file, resolveResourceFile)
    }));
  }
  return restored;
};

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
  const settings = comparison.settings;
  const parameters = settings.collinearityParams?.parameters || {};
  const mode = String(comparison.mode || 'orthogroup');
  return {
    generatedProteinComparison: adoptCanonicalPayloads
      ? comparison
      : cloneCanonicalJsonValue(comparison),
    selectedOrthogroupAlignmentFeature:
      String(settings.alignOrthogroupFeature || '').trim(),
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
    ring_width: optionalNumber(options.conservationRingWidth),
    ring_gap: optionalNumber(options.conservationRingGap)
  };
};

const projectCanonicalCircularMeasure = (measure) => {
  if (measure === null || measure === undefined) return null;
  if (!measure || typeof measure !== 'object' || Array.isArray(measure)) return measure;
  const value = Number(measure.value);
  if (!Number.isFinite(value)) return measure;
  const unit = String(measure.unit || '').trim().toLowerCase();
  if (!unit || unit === 'factor') return String(value);
  return `${value}${unit}`;
};

const projectCanonicalCircularSlot = (slot) => ({
  ...slot,
  width: projectCanonicalCircularMeasure(slot?.width),
  radius: projectCanonicalCircularMeasure(slot?.radius),
  inner_gap_px: projectCanonicalCircularMeasure(
    slot?.innerGapPx ?? slot?.inner_gap_px
  )?.replace?.(/px$/i, ''),
  outer_gap_px: projectCanonicalCircularMeasure(
    slot?.outerGapPx ?? slot?.outer_gap_px
  )?.replace?.(/px$/i, '')
});

const projectCurrentCanonicalCircularSlot = (slot) => {
  if (
    !slot || typeof slot !== 'object' || Array.isArray(slot) ||
    ['spacing', 'inner_gap_px', 'outer_gap_px', 'strict'].some((field) => (
      Object.prototype.hasOwnProperty.call(slot, field)
    ))
  ) throw new Error('Current canonical circular track slot uses an obsolete shape.');
  const projected = projectCanonicalCircularSlot(slot);
  return {
    id: String(projected.id || ''),
    renderer: String(projected.renderer || ''),
    enabled: projected.enabled !== false,
    width: projected.width ?? null, radius: projected.radius ?? null,
    inner_gap_px: projected.inner_gap_px ?? null,
    outer_gap_px: projected.outer_gap_px ?? null,
    side: projected.side ?? null, z: Number(projected.z) || 0,
    params: cloneCanonicalJsonValue(projected.params || {})
  };
};
const projectLegacyCanonicalCircularSlot = (slot) => {
  const projected = projectCanonicalCircularSlot(slot);
  if (Object.prototype.hasOwnProperty.call(slot, 'spacing')) {
    projected.spacing = projectCanonicalCircularMeasure(slot.spacing);
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

  if (sessionResourceTable) {
    return createCombinedSessionResourceFileView(
      sessionResourceTable,
      resourceIds,
      {
        name: normalizeOriginalResourceName(originalName)
          || 'canonical-circular-records.gb',
        type: 'text/plain'
      }
    );
  }

  const files = resourceIds.map((resourceId) => (
    resolveResourceFile
      ? resolveResourceFile(resourceId)
      : resourceAsLegacyFile(resources, resourceId)
  ));
  if (files.length === 1) return files[0];
  const chunks = files.map((file) => {
    if (file.encoding && file.encoding !== 'base64') {
      throw new Error(`Unsupported canonical resource encoding: ${file.encoding}`);
    }
    const decoded = base64ToBytes(file.data);
    if (decoded[decoded.length - 1] === 0x0A) return decoded;
    const terminated = new Uint8Array(decoded.length + 1);
    terminated.set(decoded);
    terminated[decoded.length] = 0x0A;
    return terminated;
  });
  const bytes = new Uint8Array(chunks.reduce((size, chunk) => size + chunk.length, 0));
  let offset = 0;
  chunks.forEach((chunk) => {
    bytes.set(chunk, offset);
    offset += chunk.length;
  });
  return {
    name: normalizeOriginalResourceName(originalName) || 'canonical-circular-records.gb',
    type: 'text/plain',
    size: bytes.length,
    lastModified: Math.max(0, ...files.map((file) => Number(file.lastModified) || 0)),
    encoding: 'base64',
    data: bytesToBase64(bytes)
  };
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
  const normalizedRecordOrdering = normalizeWebGridColumnOrdering(sourceRecords);
  const records = normalizedRecordOrdering.records;
  const reorderRecordIndexedValues = (values) => (
    normalizedRecordOrdering.sourceIndexByProjectedIndex
      .map((sourceIndex) => values?.[sourceIndex])
  );
  if (records.length === 0) throw new Error('Canonical renderRequest records are required.');
  const grouping = canonicalGrouping(renderRequest, records);
  const sourceOutputPrefixes = canonicalOutputPrefixes(
    renderRequest,
    grouping,
    records.length
  );
  const outputPrefixes = grouping === 'batch'
    ? reorderRecordIndexedValues(sourceOutputPrefixes)
    : sourceOutputPrefixes;
  const webMetadata = webFiles && typeof webFiles === 'object' && !Array.isArray(webFiles)
    ? webFiles
    : {};
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
    { adoptCanonicalPayloads }
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
        presentation: cloneCanonicalJsonValue(record.presentation)
      };
    });
    const source = records[0]?.source || {};
    if (source.kind === 'genbank') {
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
        losat_gencode: optionalPositiveInteger(
          savedMetadata.losatGencode ?? savedMetadata.losat_gencode
        ) || 1,
        losat_filename: String(
          savedMetadata.losatFilename ?? savedMetadata.losat_filename ?? ''
        ),
        definition: record.presentation?.label || '',
        record_subtitle: record.presentation?.subtitle || '',
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
    { resolveResourceFile, adoptCanonicalPayloads }
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
      throw new Error('Pairwise Match Height must be Auto or a positive finite number.');
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
      height: optionalNumber(options.depthTrackHeights?.[index]),
      large_tick_interval: optionalNumber(options.depthTrackLargeTickIntervals?.[index]),
      small_tick_interval: optionalNumber(options.depthTrackSmallTickIntervals?.[index]),
      tick_font_size: optionalNumber(options.depthTrackTickFontSizes?.[index])
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
    legend: options.output?.legend || 'right',
    multi_record_canvas: renderRequest.mode === 'circular' && grouping === 'grid',
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
    plot_title_position: options.output?.plotTitlePosition || (renderRequest.mode === 'linear' ? 'bottom' : 'none'),
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
    comparison_height: renderRequest.mode === 'linear' && comparisonHeight.status === 'valid'
      ? comparisonHeight.value
      : null,
    feature_height: overrides.default_cds_height ?? null,
    gc_height: overrides.gc_height ?? null,
    track_axis_gap: overrides.linear_track_axis_gap ?? null,
    linear_definition_line_styles: overrides.linear_definition_line_styles || {},
    linear_show_replicon: Boolean(overrides.linear_definition_show_replicon),
    linear_show_accession: overrides.linear_definition_show_accession !== false,
    linear_show_length: overrides.linear_definition_show_length !== false,
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
  const linearLayoutRows = renderRequest.schema >= 6
    ? records.map((record, index) => ({
        uid: files.linearSeqs[index]?.uid || '',
        row: Number(record.presentation?.gridRow) || index + 1
      }))
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
          ? records.some((record) => record.presentation?.gridRow != null)
          : Object.keys(renderRequest.layout || {}).length > 0,
        recordGap: renderRequest.layout?.recordGapPx ?? 24,
        rows: linearLayoutRows
      }
    : undefined;
  const projectedBlacklistText = Array.isArray(overrides.label_blacklist)
    ? overrides.label_blacklist.join(', ')
    : String(overrides.label_blacklist || '');
  return {
    mode: renderRequest.mode,
    inputType: records[0]?.source?.kind === 'gffFasta' ? 'gff' : 'gb',
    files,
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
          selectedOrthogroupAlignmentFeature:
            projectedProteinPipeline.selectedOrthogroupAlignmentFeature
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
export const promoteCanonicalRenderRequestToCurrent = (request) => {
  const promoted = cloneCanonicalJsonValue(request);
  if (promoted.schema === CANONICAL_REQUEST_SCHEMA) return promoted;
  if (promoted.schema !== 5) {
    throw new Error('Only canonical renderRequest schema 5 can be promoted to schema 6.');
  }
  const linearRows = promoted.mode === 'linear'
    ? (promoted.layout?.multiRecordPositions || []).map((token) => {
        const split = String(token).lastIndexOf('@');
        return Number(String(token).slice(split + 1)) || null;
      })
    : [];
  promoted.schema = CANONICAL_REQUEST_SCHEMA;
  (promoted.records || []).forEach((record, index) => {
    record.cardinality = promoted.mode === 'linear' &&
      !record.selector && !record.region
      ? 'all'
      : 'exactly_one';
    if (linearRows[index]) record.presentation.gridRow = linearRows[index];
  });
  if (promoted.mode === 'linear' && promoted.layout) {
    delete promoted.layout.multiRecordPositions;
  }
  return promoted;
};
const normalizePublicationRequestAliases = (request) => {
  const normalized = request?.schema === 5
    ? promoteCanonicalRenderRequestToCurrent(request)
    : cloneCanonicalJsonValue(request);
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
