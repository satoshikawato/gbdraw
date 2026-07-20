import { buildDefaultColorOverrideTsv } from '../app/color-utils.js';
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
  buildCircularTrackSlotSpec,
  CIRCULAR_TRACK_RENDERERS,
  normalizeCircularTrackSlot,
  parseCircularTrackSlotSpecs
} from '../app/circular-track-slots.js';
import {
  buildLinearTrackSlotSpec,
  LINEAR_TRACK_RENDERERS,
  LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  linearTrackAxisIndexForEnabledSlots,
  migrateLinearTrackSlotsToCurrentSchema,
  parseLinearTrackSlotSpecs
} from '../app/linear-track-slots.js';
import {
  isRecordMajorDepthFileMatrix,
  normalizeRecordMajorDepthFileRows,
  parseDepthTrackIndexIdentity
} from '../app/depth-track-state.js';
import { validateTrackSlotBindingInvariants } from '../app/track-slot-validation.js';
import { annotationOptionsPayload, normalizeAnnotationSets } from '../app/annotations/state.js';
import { classifyOptionalPositiveNumber } from '../utils/optional-positive-number.js';
import {
  defaultFeatureRendering,
  normalizeFeatureRenderingMap
} from '../utils/feature-rendering.js';

export const CANONICAL_REQUEST_SCHEMA = 3;

const safePrefix = (value) => {
  const normalized = String(value || '').trim().replace(/[\\/]+/g, '_');
  return normalized && normalized !== '.' && normalized !== '..' ? normalized : 'out';
};

const optionalNumber = (value) => {
  if (value === null || value === undefined || String(value).trim() === '') return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) ? numeric : null;
};

const optionalPositiveInteger = (value) => {
  const numeric = optionalNumber(value);
  return Number.isInteger(numeric) && numeric > 0 ? numeric : null;
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

const textToBase64 = (text) => {
  const bytes = new TextEncoder().encode(String(text));
  let binary = '';
  const chunkSize = 0x8000;
  for (let index = 0; index < bytes.length; index += chunkSize) {
    binary += String.fromCharCode(...bytes.subarray(index, index + chunkSize));
  }
  return btoa(binary);
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

  const addFile = (resourceId, kind, entry) => {
    if (!entry || typeof entry !== 'object' || Array.isArray(entry)) {
      throw new Error(`Canonical resource ${resourceId} is missing.`);
    }
    if (resources[resourceId]) return resourceId;
    resources[resourceId] = {
      kind,
      name: normalizeResourceName(resourceId, entry.name),
      type: String(entry.type || 'application/octet-stream'),
      size: Number(entry.size) || 0,
      lastModified: Number(entry.lastModified) || 0,
      encoding: entry.encoding || 'base64',
      data: entry.data
    };
    const originalName = normalizeOriginalResourceName(entry.name);
    if (originalName) resourceOriginalNames[resourceId] = originalName;
    return resourceId;
  };

  const addText = (resourceId, kind, name, text) => {
    if (resources[resourceId]) return resourceId;
    const normalized = String(text || '');
    const bytes = new TextEncoder().encode(normalized);
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

  return { resources, resourceOriginalNames, addFile, addText };
};

const fileRef = (resourceId) => ({ resourceId, representation: 'file' });

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

const presentationPayload = ({ label = null, subtitle = null } = {}) => ({
  label: String(label || '').trim() || null,
  subtitle: String(subtitle || '').trim() || null,
  reverseComplement: false,
  gridRow: null,
  gridColumn: null
});

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
    return (filesData.linearSeqs || []).map((seq, index) => {
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
      return {
        recordKey: String(seq.uid || `record-${index + 1}`),
        source,
        selector: region ? null : selectorPayload(seq.region_record_id),
        region,
        presentation: {
          ...presentationPayload({ label: seq.definition, subtitle: seq.record_subtitle }),
          reverseComplement: region ? false : Boolean(seq.region_reverse)
        }
      };
    });
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
  if (!state.form.multi_record_canvas && knownRecords.length > 1) {
    throw new Error(
      'A non-canvas Circular run with multiple outputs cannot be represented by canonical request schema 1.'
    );
  }
  const selectedRecords = state.form.multi_record_canvas && knownRecords.length > 0
    ? knownRecords
    : [null];
  return selectedRecords.map((record, index) => ({
    recordKey: `record-${index + 1}`,
    source,
    selector: selectorPayload(record?.selector),
    region: null,
    presentation: presentationPayload()
  }));
};

const buildConfigOverrides = (state) => {
  const { form, adv } = state;
  const circular = state.mode.value === 'circular';
  const comparisonHeight = classifyOptionalPositiveNumber(adv.comparison_height);
  if (!circular && comparisonHeight.status === 'invalid') {
    throw new Error('Pairwise Match Height must be Auto or a positive finite number.');
  }
  return {
    block_stroke_width: optionalNumber(adv.block_stroke_width),
    block_stroke_color: adv.block_stroke_color || null,
    circular_axis_stroke_color: circular ? (adv.axis_stroke_color || null) : null,
    circular_axis_stroke_width: circular ? optionalNumber(adv.axis_stroke_width) : null,
    linear_axis_stroke_color: circular ? null : (adv.axis_stroke_color || null),
    linear_axis_stroke_width: circular ? null : optionalNumber(adv.axis_stroke_width),
    line_stroke_color: adv.line_stroke_color || null,
    line_stroke_width: optionalNumber(adv.line_stroke_width),
    circular_definition_font_size: circular ? optionalNumber(adv.def_font_size) : null,
    linear_definition_font_size: circular ? null : optionalNumber(adv.def_font_size),
    linear_definition_line_styles: circular ? null : adv.linear_definition_line_styles,
    linear_definition_show_replicon: circular ? null : Boolean(adv.linear_show_replicon),
    linear_definition_show_accession: circular ? null : Boolean(adv.linear_show_accession),
    linear_definition_show_length: circular ? null : Boolean(adv.linear_show_length),
    plot_title_font_size: optionalNumber(adv.plot_title_font_size),
    label_font_size: optionalNumber(adv.label_font_size),
    circular_label_spacing: circular ? optionalNumber(adv.circular_label_spacing) : null,
    linear_label_spacing: circular ? null : optionalNumber(adv.linear_label_spacing),
    label_rendering: adv.label_rendering || 'auto',
    circular_label_placement: circular ? (adv.circular_label_placement || 'horizontal') : null,
    label_placement: circular ? null : (adv.label_placement || 'auto'),
    label_rotation: circular ? null : optionalNumber(adv.label_rotation),
    show_gc: circular ? !form.suppress_gc : Boolean(form.show_gc),
    show_skew: circular ? !form.suppress_skew : Boolean(form.show_skew),
    show_depth: Boolean(form.show_depth),
    gc_content_mode: adv.gc_content_mode || 'deviation',
    gc_content_min_percent: optionalNumber(adv.gc_content_min_percent),
    gc_content_max_percent: optionalNumber(adv.gc_content_max_percent),
    gc_content_show_axis: Boolean(adv.gc_content_show_axis),
    gc_content_show_ticks: Boolean(adv.gc_content_show_ticks),
    gc_content_large_tick_interval: optionalNumber(adv.gc_content_tick_interval),
    gc_content_small_tick_interval: optionalNumber(adv.gc_content_small_tick_interval),
    gc_content_tick_font_size: optionalNumber(adv.gc_content_tick_font_size),
    depth_color: adv.depth_color || null,
    depth_height: optionalNumber(adv.depth_height),
    depth_min: optionalNumber(adv.depth_min),
    depth_max: optionalNumber(adv.depth_max),
    depth_normalize: Boolean(adv.depth_normalize),
    depth_show_axis: Boolean(adv.depth_show_axis),
    depth_show_ticks: Boolean(adv.depth_show_ticks),
    depth_large_tick_interval: optionalNumber(adv.depth_tick_interval),
    depth_small_tick_interval: optionalNumber(adv.depth_small_tick_interval),
    depth_tick_font_size: optionalNumber(adv.depth_tick_font_size),
    depth_share_axis: Boolean(adv.depth_share_axis),
    show_labels: circular ? form.labels_mode !== 'none' : form.show_labels_linear,
    align_center: circular ? null : Boolean(form.align_center),
    keep_definition_left_aligned: circular ? null : Boolean(form.keep_definition_left_aligned),
    linear_track_layout: circular ? null : (form.linear_track_layout || 'middle'),
    linear_track_axis_gap: circular ? null : optionalNumber(adv.track_axis_gap),
    linear_ruler_on_axis: circular ? null : Boolean(form.linear_ruler_on_axis),
    track_type: circular ? form.track_type : null,
    strandedness: Boolean(form.separate_strands),
    resolve_overlaps: Boolean(adv.resolve_overlaps),
    allow_inner_labels: circular ? form.labels_mode === 'both' : null,
    label_blacklist: state.filterMode.value === 'Blacklist' ? state.manualBlacklist.value : null,
    comparison_height: circular || comparisonHeight.status === 'auto' ? null : comparisonHeight.value,
    default_cds_height: circular ? null : optionalNumber(adv.feature_height),
    gc_height: circular ? null : optionalNumber(adv.gc_height),
    scale_style: circular ? null : form.scale_style,
    scale_stroke_color: circular ? null : (adv.scale_stroke_color || null),
    scale_label_color: circular ? null : (adv.ruler_label_color || null),
    scale_stroke_width: circular ? null : optionalNumber(adv.scale_stroke_width),
    scale_font_size: circular ? null : optionalNumber(adv.scale_font_size),
    ruler_label_font_size: circular ? null : optionalNumber(adv.scale_font_size),
    scale_interval: optionalNumber(adv.scale_interval),
    tick_label_font_size: circular ? optionalNumber(adv.tick_label_font_size) : null,
    outer_label_x_radius_offset: circular ? optionalNumber(adv.outer_label_x_offset) : null,
    outer_label_y_radius_offset: circular ? optionalNumber(adv.outer_label_y_offset) : null,
    inner_label_x_radius_offset: circular ? optionalNumber(adv.inner_label_x_offset) : null,
    inner_label_y_radius_offset: circular ? optionalNumber(adv.inner_label_y_offset) : null,
    pairwise_match_style: circular ? null : adv.pairwise_match_style,
    legend_box_size: optionalNumber(adv.legend_box_size),
    legend_font_size: optionalNumber(adv.legend_font_size),
    normalize_length: circular ? null : Boolean(form.normalize_length)
  };
};

const addGeneratedTableResources = (state, resources, diagramOptions) => {
  const paletteName = String(state.selectedPalette.value || 'default');
  const paletteColors = state.normalizePaletteColors(
    state.paletteDefinitions.value?.[paletteName] || state.paletteDefinitions.value?.default || {}
  );
  const defaultColors = buildDefaultColorOverrideTsv({
    colors: state.currentColors.value,
    paletteColors
  });
  const specificColors = serializeSpecificRules(state.manualSpecificRules);
  let defaultColorsFile = null;
  let colorTableFile = null;
  if (defaultColors.trim()) {
    defaultColorsFile = fileRef(resources.addText(
      'colors-default-colors-file', 'colors-default-colors-file', 'default-colors.tsv', `${defaultColors}\n`
    ));
  }
  if (specificColors.trim()) {
    colorTableFile = fileRef(resources.addText(
      'colors-color-table-file', 'colors-color-table-file', 'specific-colors.tsv', specificColors
    ));
  }
  diagramOptions.colors = {
    colorTable: null,
    colorTableFile,
    defaultColors: null,
    defaultColorsPalette: paletteName,
    defaultColorsFile
  };

  const visibility = serializeFeatureVisibilityRules(state.featureVisibilityRules?.value || []);
  if (visibility.trim()) {
    diagramOptions.featureVisibilityTableFile = fileRef(resources.addText(
      'feature-visibility-table-file', 'feature-visibility-table-file', 'feature-visibility.tsv', visibility
    ));
  }
  if (state.filterMode.value === 'Whitelist' && state.manualWhitelist.length > 0) {
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
  if (priority) {
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

const buildDepthResources = ({ state, filesData, resources, diagramOptions, recordCount }) => {
  const rows = state.mode.value === 'linear'
    ? (filesData.linearSeqs || []).map((seq) => Array.isArray(seq.depth) ? seq.depth : (seq.depth ? [seq.depth] : []))
    : normalizeRecordMajorDepthFileRows(filesData.c_depth, recordCount);
  if (rows.every((row) => row.length === 0)) return;
  if (
    state.mode.value === 'circular' &&
    isRecordMajorDepthFileMatrix(filesData.c_depth) &&
    filesData.c_depth.length !== recordCount
  ) {
    throw new Error(
      `Circular Depth matrix has ${filesData.c_depth.length} record rows; expected ${recordCount}.`
    );
  }
  diagramOptions.depthTrackFiles = rows.map((row, rowIndex) => row.map((entry, columnIndex) => {
    if (!entry) return null;
    const id = `depth-track-files-${rowIndex + 1}-${columnIndex + 1}`;
    return fileRef(resources.addFile(id, 'depth-track-file', entry));
  }));
  const tracks = Array.isArray(state.adv.depth_tracks) ? state.adv.depth_tracks : [];
  diagramOptions.depthTrackLabels = tracks.map((track, index) => String(track?.label || `Depth ${index + 1}`));
  diagramOptions.depthTrackColors = tracks.map((track) => String(track?.color || state.adv.depth_color));
  if (state.mode.value === 'linear') {
    diagramOptions.depthTrackHeights = tracks.map((track) => optionalNumber(track?.height));
  }
  diagramOptions.depthTrackLargeTickIntervals = tracks.map((track) => optionalNumber(track?.large_tick_interval));
  diagramOptions.depthTrackSmallTickIntervals = tracks.map((track) => optionalNumber(track?.small_tick_interval));
  diagramOptions.depthTrackTickFontSizes = tracks.map((track) => optionalNumber(track?.tick_font_size));
};

const buildTracks = (state) => {
  const circular = state.mode.value === 'circular';
  const storedLinearAxisIndex = optionalNumber(state.adv.linear_track_slots_axis_index);
  const emittedLinearAxisIndex = (
    !circular &&
    state.adv.linear_track_slots_enabled &&
    storedLinearAxisIndex !== null
  )
    ? linearTrackAxisIndexForEnabledSlots(
        state.adv.linear_track_slots,
        storedLinearAxisIndex
      )
    : null;
  return {
    circularTrackSlots: circular && state.adv.circular_track_slots_enabled
      ? state.adv.circular_track_slots
          .filter((slot) => slot?.enabled !== false)
          .map((slot) => buildCircularTrackSlotSpec(slot, state.adv.nt, state.form.track_type))
      : null,
    circularTrackAxisIndex: circular && state.adv.circular_track_slots_enabled
      ? optionalNumber(state.adv.circular_track_slots_axis_index)
      : null,
    linearTrackSlots: !circular && state.adv.linear_track_slots_enabled
      ? state.adv.linear_track_slots
          .filter((slot) => slot?.enabled !== false)
          .map((slot) => buildLinearTrackSlotSpec(slot))
      : null,
    linearTrackAxisIndex: circular ? null : emittedLinearAxisIndex,
    centerReservedRadius: circular ? optionalNumber(state.adv.center_reserved_radius) : null
  };
};

const buildComparisons = ({ state, filesData, resources }) => {
  if (state.mode.value !== 'linear') return [];
  const comparisons = [];
  const indexByUid = new Map((filesData.linearSeqs || []).map((seq, index) => [String(seq.uid || ''), index]));
  const explicitComparisons = Array.isArray(filesData.linearComparisons)
    ? filesData.linearComparisons
    : [];
  explicitComparisons.forEach((comparison, index) => {
    if (!comparison?.file) return;
    const queryRecordIndex = indexByUid.get(String(comparison.queryUid || ''));
    const subjectRecordIndex = indexByUid.get(String(comparison.subjectUid || ''));
    if (!Number.isInteger(queryRecordIndex) || !Number.isInteger(subjectRecordIndex)) return;
    const id = `comparison-nucleotide-${index + 1}`;
    comparisons.push({
      kind: 'nucleotideBlast',
      resourceId: resources.addFile(id, 'nucleotide-blast', comparison.file),
      queryRecordIndex,
      subjectRecordIndex
    });
  });
  if (explicitComparisons.length === 0) {
    (filesData.linearSeqs || []).forEach((seq, index) => {
      if (!seq.blast || index + 1 >= filesData.linearSeqs.length) return;
      const id = `comparison-nucleotide-${index + 1}`;
      comparisons.push({
        kind: 'nucleotideBlast',
        resourceId: resources.addFile(id, 'nucleotide-blast', seq.blast),
        queryRecordIndex: index,
        subjectRecordIndex: index + 1
      });
    });
  }
  const generatedPairs = explicitComparisons
    .filter((comparison) => comparison?.source === 'losat')
    .map((comparison) => ({
      queryRecordIndex: indexByUid.get(String(comparison.queryUid || '')),
      subjectRecordIndex: indexByUid.get(String(comparison.subjectUid || ''))
    }))
    .filter((pair) => Number.isInteger(pair.queryRecordIndex) && Number.isInteger(pair.subjectRecordIndex));
  if (
    state.losatProgram.value === 'blastp' &&
    (generatedPairs.length > 0 || state.blastSource.value === 'losat')
  ) {
    const blastp = state.losat.blastp || {};
    comparisons.push({
      kind: 'generatedProteinComparison',
      mode: String(blastp.mode || 'orthogroup'),
      pairs: generatedPairs,
      settings: {
        collinearityParams: {
          kind: 'lossless',
          parameters: {
            minAnchors: Number(blastp.collinearMinAnchors) || 1,
            maxUnitGap: Number(blastp.collinearMaxGeneGap) || 0,
            maxDiagonalDrift: Number(blastp.collinearMaxDiagonalDrift) || 0,
            maxConflicts: Number(blastp.collinearMaxConflictsInMergeGap) || 0,
            mergeOrientation: 'either'
          }
        },
        collinearityUnitMode: String(blastp.collinearUnitMode || 'auto'),
        collinearityAnchorMode: 'rbh',
        collinearitySearchScope: String(blastp.collinearSearchScope || 'adjacent'),
        collinearityColorMode: String(blastp.collinearColorMode || 'orientation'),
        losatpBin: 'losat',
        ncbiBlastpBin: null,
        losatpThreads: optionalPositiveInteger(state.losat.threadsPerJob),
        proteinBlastpMaxHits: Number(blastp.maxHits) || 5,
        proteinBlastpCandidateLimit: optionalPositiveInteger(blastp.candidateLimit),
        orthogroupMembershipMode: String(blastp.orthogroupMembershipMode || 'anchor_core_v1'),
        orthogroupMemberMaxHits: Number(blastp.orthogroupMemberMaxHits) || 5,
        collinearMaxParalogLinksPerOrthogroup: Number(blastp.collinearMaxParalogLinksPerOrthogroup) || 2,
        alignOrthogroupFeature: String(state.selectedOrthogroupAlignmentFeature.value || '').trim() || null
      }
    });
  }
  return comparisons;
};

const buildLayout = (state, filesData) => {
  if (state.mode.value === 'linear') {
    if (!state.linearRecordLayoutEnabled?.value) return {};
    const rows = new Map(
      (state.linearRecordRows || []).map((entry) => [String(entry?.uid || ''), Number(entry?.row)])
    );
    return {
      recordGapPx: Math.max(0, Number(state.linearRecordGap?.value) || 0),
      multiRecordPositions: (filesData.linearSeqs || []).map((sequence, index) => {
        const row = rows.get(String(sequence?.uid || ''));
        return `#${index + 1}@${Number.isInteger(row) && row > 0 ? row : index + 1}`;
      })
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
    multiRecordSizeMode: String(state.adv.multi_record_size_mode || 'auto'),
    multiRecordMinRadiusRatio: Number(state.adv.multi_record_min_radius_ratio) || 0.55,
    multiRecordColumnGapRatio: Number(state.adv.multi_record_column_gap_ratio) || 0,
    multiRecordRowGapRatio: Number(state.adv.multi_record_row_gap_ratio) || 0,
    multiRecordPositions: positions.length > 0 ? positions : null
  };
};

export const buildCanonicalSessionRequest = ({ state, filesData }) => {
  const resources = createResourceBuilder();
  const webFiles = {};
  const records = buildRecords({ state, filesData, resources });
  if (records.length === 0) throw new Error('A canonical request requires at least one record.');

  const diagramOptions = {
    configOverrides: buildConfigOverrides(state),
    tracks: buildTracks(state),
    output: {
      outputPrefix: safePrefix(state.form.prefix),
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
    evalue: Number(state.adv.evalue),
    bitscore: Number(state.adv.min_bitscore),
    identity: Number(state.adv.identity),
    alignmentLength: Number(state.adv.alignment_length) || 0
  };
  if (Array.isArray(state.annotationSets) && state.annotationSets.length > 0) {
    diagramOptions.annotations = annotationOptionsPayload(state.annotationSets);
  }
  if (state.mode.value === 'circular') {
    diagramOptions.keepFullDefinitionWithPlotTitle = Boolean(state.adv.keep_full_definition_with_plot_title);
    diagramOptions.species = String(state.form.species || '').trim() || null;
    diagramOptions.strain = String(state.form.strain || '').trim() || null;
    const conservation = Array.isArray(filesData.c_conservation_blasts)
      ? filesData.c_conservation_blasts
      : [];
    if (conservation.length > 0) {
      diagramOptions.conservationBlastFiles = conservation.map((entry, index) => fileRef(
        resources.addFile(`conservation-blast-files-${index + 1}`, 'conservation-blast-file', entry)
      ));
      const comparisonSources = Array.isArray(filesData.c_conservation_sequence_sources)
        ? filesData.c_conservation_sequence_sources
        : [];
      if (comparisonSources.some(Boolean)) {
        webFiles.conservationSequenceSources = comparisonSources.map((entry, index) => (
          entry
            ? resources.addFile(
                `conservation-sequence-sources-${index + 1}`,
                'conservation-sequence-source',
                entry
              )
            : null
        ));
      }
      diagramOptions.conservationReference = String(state.circularConservation.reference || 'auto');
      diagramOptions.conservationLabels = String(state.circularConservation.labels || '')
        .split(',').map((value) => value.trim()).filter(Boolean);
      diagramOptions.conservationColors = (state.circularConservation.series || [])
        .map((entry) => String(entry?.color || '').trim()).filter(Boolean);
      diagramOptions.conservationRingWidth = optionalNumber(state.circularConservation.ring_width);
      diagramOptions.conservationRingGap = optionalNumber(state.circularConservation.ring_gap);
    }
  } else {
    diagramOptions.pairwiseMatchStyle = String(state.adv.pairwise_match_style || 'ribbon');
  }

  addGeneratedTableResources(state, resources, diagramOptions);
  buildDepthResources({ state, filesData, resources, diagramOptions, recordCount: records.length });

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
      losatGencode: optionalPositiveInteger(sequence?.losat_gencode) || 1,
      losatFilename: String(sequence?.losat_filename || '')
    }));
  }

  return {
    renderRequest: {
      schema: CANONICAL_REQUEST_SCHEMA,
      mode: state.mode.value,
      records,
      diagramOptions,
      layout: buildLayout(state, filesData),
      comparisons: buildComparisons({ state, filesData, resources }),
      output: {
        prefix: safePrefix(state.form.prefix),
        formats: ['interactive_svg'],
        overwrite: true,
        interactiveMetadataPolicy: 'auto'
      }
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
  let binary;
  try {
    binary = atob(entry.data);
  } catch (error) {
    throw new Error(`Canonical resource ${resourceId} contains invalid base64 data.`, { cause: error });
  }
  const bytes = Uint8Array.from(binary, (character) => character.charCodeAt(0));
  return new TextDecoder('utf-8', { fatal: true }).decode(bytes);
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

const sharedLengthConfigValue = (config, path) => {
  const value = nestedConfigValue(config, path);
  if (!value || typeof value !== 'object' || Array.isArray(value)) return value;
  // The Web UI has one control for values that Python stores per genome length.
  // Project them only when both variants encode the same explicit setting.
  const shortValue = value.short;
  const longValue = value.long;
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

const projectFullConfigOverrides = (config, mode) => {
  if (!config || typeof config !== 'object' || Array.isArray(config)) return {};
  const paths = {
    block_stroke_color: 'objects.features.block_stroke_color',
    circular_axis_stroke_color: 'objects.axis.circular.stroke_color',
    linear_axis_stroke_color: 'objects.axis.linear.stroke_color',
    line_stroke_color: 'objects.features.line_stroke_color',
    circular_definition_font_size: 'objects.definition.circular.font_size',
    plot_title_font_size: 'objects.definition.circular.plot_title_font_size',
    show_gc: 'canvas.show_gc',
    show_skew: 'canvas.show_skew',
    show_depth: 'canvas.show_depth',
    show_labels: 'canvas.show_labels',
    strandedness: 'canvas.strandedness',
    resolve_overlaps: 'canvas.resolve_overlaps',
    track_type: 'canvas.circular.track_type',
    allow_inner_labels: 'canvas.circular.allow_inner_labels',
    align_center: 'canvas.linear.align_center',
    keep_definition_left_aligned: 'canvas.linear.keep_definition_left_aligned',
    linear_track_layout: 'canvas.linear.track_layout',
    linear_track_axis_gap: 'canvas.linear.track_axis_gap',
    linear_ruler_on_axis: 'canvas.linear.ruler_on_axis',
    comparison_height: 'canvas.linear.comparison_height',
    gc_height: 'canvas.linear.default_gc_height',
    depth_height: 'canvas.linear.depth_height',
    normalize_length: 'canvas.linear.normalize_length',
    label_rendering: 'labels.rendering',
    circular_label_spacing: 'labels.spacing.circular',
    circular_label_placement: 'labels.circular.placement',
    linear_label_spacing: 'labels.spacing.linear',
    label_placement: 'labels.linear.placement',
    label_rotation: 'labels.linear.rotation',
    label_blacklist: 'labels.filtering.blacklist_keywords',
    linear_definition_line_styles: 'objects.definition.linear.line_styles',
    linear_definition_show_replicon: 'objects.definition.linear.show_replicon',
    linear_definition_show_accession: 'objects.definition.linear.show_accession',
    linear_definition_show_length: 'objects.definition.linear.show_length',
    gc_content_mode: 'objects.gc_content.mode',
    gc_content_min_percent: 'objects.gc_content.min_percent',
    gc_content_max_percent: 'objects.gc_content.max_percent',
    gc_content_show_axis: 'objects.gc_content.show_axis',
    gc_content_show_ticks: 'objects.gc_content.show_ticks',
    gc_content_large_tick_interval: 'objects.gc_content.large_tick_interval',
    gc_content_small_tick_interval: 'objects.gc_content.small_tick_interval',
    gc_content_tick_font_size: 'objects.gc_content.tick_font_size',
    depth_color: 'objects.depth.fill_color',
    depth_min: 'objects.depth.min_depth',
    depth_max: 'objects.depth.max_depth',
    depth_normalize: 'objects.depth.normalize',
    depth_show_axis: 'objects.depth.show_axis',
    depth_show_ticks: 'objects.depth.show_ticks',
    depth_large_tick_interval: 'objects.depth.large_tick_interval',
    depth_small_tick_interval: 'objects.depth.small_tick_interval',
    depth_tick_font_size: 'objects.depth.tick_font_size',
    depth_share_axis: 'objects.depth.share_axis',
    scale_style: 'objects.scale.style',
    scale_stroke_color: 'objects.scale.stroke_color',
    scale_label_color: 'objects.scale.label_color',
    scale_stroke_width: 'objects.scale.stroke_width',
    scale_interval: 'objects.scale.interval',
    tick_label_font_size: 'objects.ticks.tick_labels.font_size',
    outer_label_x_radius_offset: 'labels.unified_adjustment.outer_labels.x_radius_offset',
    outer_label_y_radius_offset: 'labels.unified_adjustment.outer_labels.y_radius_offset',
    inner_label_x_radius_offset: 'labels.unified_adjustment.inner_labels.x_radius_offset',
    inner_label_y_radius_offset: 'labels.unified_adjustment.inner_labels.y_radius_offset',
    pairwise_match_style: 'objects.blast_match.style'
  };
  const projected = Object.fromEntries(
    Object.entries(paths)
      .map(([name, path]) => [name, nestedConfigValue(config, path)])
      .filter(([, value]) => value !== undefined)
  );
  const sharedLengthPaths = {
    block_stroke_width: 'objects.features.block_stroke_width',
    circular_axis_stroke_width: 'objects.axis.circular.stroke_width',
    linear_axis_stroke_width: 'objects.axis.linear.stroke_width',
    line_stroke_width: 'objects.features.line_stroke_width',
    linear_definition_font_size: 'objects.definition.linear.font_size',
    default_cds_height: 'canvas.linear.default_cds_height',
    legend_box_size: 'objects.legends.color_rect_size',
    legend_font_size: 'objects.legends.font_size',
    scale_font_size: 'objects.scale.font_size',
    ruler_label_font_size: 'objects.scale.ruler_label_font_size',
    label_font_size: mode === 'linear' ? 'labels.font_size.linear' : 'labels.font_size'
  };
  Object.entries(sharedLengthPaths).forEach(([name, path]) => {
    const value = sharedLengthConfigValue(config, path);
    if (value !== undefined) projected[name] = value;
  });
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
  spacing: projectCanonicalCircularMeasure(slot?.spacing),
  inner_gap_px: projectCanonicalCircularMeasure(slot?.innerGapPx ?? slot?.inner_gap_px),
  outer_gap_px: projectCanonicalCircularMeasure(slot?.outerGapPx ?? slot?.outer_gap_px)
});

const combineCircularGenbankResources = (resources, records, originalName = '') => {
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

  const files = resourceIds.map((resourceId) => resourceAsLegacyFile(resources, resourceId));
  if (files.length === 1) return files[0];
  const binary = files
    .map((file) => {
      if (file.encoding && file.encoding !== 'base64') {
        throw new Error(`Unsupported canonical resource encoding: ${file.encoding}`);
      }
      const decoded = atob(String(file.data || ''));
      return decoded.endsWith('\n') ? decoded : `${decoded}\n`;
    })
    .join('');
  return {
    name: normalizeOriginalResourceName(originalName) || 'canonical-circular-records.gb',
    type: 'text/plain',
    size: binary.length,
    lastModified: Math.max(0, ...files.map((file) => Number(file.lastModified) || 0)),
    encoding: 'base64',
    data: btoa(binary)
  };
};

export const projectCanonicalSessionRequest = ({
  renderRequest,
  resources: canonicalResources,
  webFiles = {},
  legacyFiles = null,
  fileBindings = [],
  linearTrackSlotSchemaVersion = LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  repairInvalidComparisonHeight = false
}) => {
  if (!renderRequest || ![1, 2, CANONICAL_REQUEST_SCHEMA].includes(renderRequest.schema)) {
    throw new Error('Unsupported canonical renderRequest schema.');
  }
  if (!['circular', 'linear'].includes(renderRequest.mode)) {
    throw new Error('Unsupported canonical renderRequest mode.');
  }
  const records = Array.isArray(renderRequest.records) ? renderRequest.records : [];
  if (records.length === 0) throw new Error('Canonical renderRequest records are required.');
  const webMetadata = webFiles && typeof webFiles === 'object' && !Array.isArray(webFiles)
    ? webFiles
    : {};
  const storedResourceOriginalNames = webMetadata.resourceOriginalNames;
  const resources = resourcesWithOriginalNames(canonicalResources, {
    ...legacyResourceOriginalNames({ renderRequest, legacyFiles, fileBindings }),
    ...(storedResourceOriginalNames && typeof storedResourceOriginalNames === 'object' &&
      !Array.isArray(storedResourceOriginalNames) ? storedResourceOriginalNames : {})
  });
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
  const files = { linearSeqs: [] };
  if (renderRequest.mode === 'circular') {
    const source = records[0]?.source || {};
    if (source.kind === 'genbank') {
      files.c_gb = combineCircularGenbankResources(resources, records, circularInputOriginalName);
    }
    if (source.kind === 'gffFasta') {
      files.c_gff = resourceAsLegacyFile(resources, source.gffResourceId);
      files.c_fasta = resourceAsLegacyFile(resources, source.fastaResourceId);
    }
  } else {
    files.linearSeqs = records.map((record, index) => {
      const source = record.source || {};
      const region = record.region || null;
      const selector = region?.selector || record.selector;
      const savedMetadata = savedLinearRecordMetadataByKey.get(String(record.recordKey || '')) ||
        savedLinearRecordMetadata[index] || legacyLinearSequences[index] || {};
      return {
        uid: String(record.recordKey || `canonical-seq-${index + 1}`),
        gb: source.kind === 'genbank' ? resourceAsLegacyFile(resources, source.resourceId) : null,
        gff: source.kind === 'gffFasta' ? resourceAsLegacyFile(resources, source.gffResourceId) : null,
        fasta: source.kind === 'gffFasta' ? resourceAsLegacyFile(resources, source.fastaResourceId) : null,
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
    (renderRequest.comparisons || [])
      .filter((comparison) => comparison?.kind === 'nucleotideBlast')
      .forEach((comparison, index) => {
        const queryIndex = Number.isInteger(Number(comparison.queryRecordIndex))
          ? Number(comparison.queryRecordIndex)
          : index;
        const subjectIndex = Number.isInteger(Number(comparison.subjectRecordIndex))
          ? Number(comparison.subjectRecordIndex)
          : index + 1;
        const file = resourceAsLegacyFile(resources, comparison.resourceId);
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
    (renderRequest.comparisons || [])
      .filter((comparison) => comparison?.kind === 'generatedProteinComparison')
      .flatMap((comparison) => Array.isArray(comparison.pairs) ? comparison.pairs : [])
      .forEach((pair, index) => {
        const queryIndex = Number(pair?.queryRecordIndex);
        const subjectIndex = Number(pair?.subjectRecordIndex);
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

  const options = renderRequest.diagramOptions || {};
  const depthRows = Array.isArray(options.depthTrackFiles) ? options.depthTrackFiles : [];
  if (depthRows.length > 0) {
    if (depthRows.length !== records.length || depthRows.some((row) => !Array.isArray(row))) {
      throw new Error(
        `Canonical ${renderRequest.mode === 'circular' ? 'Circular' : 'Linear'} Depth matrix must contain one row per record (${records.length}).`
      );
    }
  }
  if (renderRequest.mode === 'circular' && depthRows.length > 0) {
    files.c_depth = normalizeRecordMajorDepthFileRows(
      depthRows.map((row) => row.map((ref) => (
        ref?.resourceId ? resourceAsLegacyFile(resources, ref.resourceId) : null
      ))),
      records.length
    );
  }
  if (renderRequest.mode === 'linear') {
    depthRows.forEach((row, index) => {
      if (!files.linearSeqs[index] || !Array.isArray(row)) return;
      const depth = row
        .map((ref) => ref?.resourceId ? resourceAsLegacyFile(resources, ref.resourceId) : null);
      files.linearSeqs[index].depth = depth.length > 1 ? depth : (depth[0] || null);
    });
  }
  const defaultColorsRef = options.colors?.defaultColorsFile || options.colors?.defaultColors;
  let projectedDefaultColors = {};
  if (defaultColorsRef?.resourceId) {
    files.d_color = resourceAsLegacyFile(resources, defaultColorsRef.resourceId);
    projectedDefaultColors = parseColorTable(
      resourceTextFromRef(resources, defaultColorsRef)
    ).colors;
  }
  const colorTableRef = options.colors?.colorTableFile || options.colors?.colorTable;
  let projectedSpecificRules = [];
  if (colorTableRef?.resourceId) {
    files.t_color = resourceAsLegacyFile(resources, colorTableRef.resourceId);
    projectedSpecificRules = parseSpecificRules(
      resourceTextFromRef(resources, colorTableRef)
    ).rules.map(({ fromFile: _fromFile, ...rule }) => rule);
  }
  let projectedWhitelist = [];
  if (options.labelWhitelistFile?.resourceId) {
    files.whitelist = resourceAsLegacyFile(resources, options.labelWhitelistFile.resourceId);
    projectedWhitelist = parseWhitelistRules(
      resourceTextFromRef(resources, options.labelWhitelistFile)
    ).rules;
  }
  const qualifierPriorityRef = options.qualifierPriorityFile || options.qualifierPriorityTable;
  let projectedPriorityRules = [];
  if (qualifierPriorityRef?.resourceId) {
    files.qualifier_priority = resourceAsLegacyFile(resources, qualifierPriorityRef.resourceId);
    projectedPriorityRules = parsePriorityRules(
      resourceTextFromRef(resources, qualifierPriorityRef)
    ).rules;
  }
  const projectedFeatureVisibilityRules = options.featureVisibilityTableFile?.resourceId
    ? parseFeatureVisibilityRules(
        resourceTextFromRef(resources, options.featureVisibilityTableFile)
      ).rules
    : [];
  const projectedLabelOverrideRows = options.labelOverrideFile?.resourceId
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
      .map((ref) => ref?.resourceId ? resourceAsLegacyFile(resources, ref.resourceId) : null)
      .filter(Boolean);
  }
  if (renderRequest.mode === 'circular' && Array.isArray(webMetadata.conservationSequenceSources)) {
    files.c_conservation_sequence_sources = webMetadata.conservationSequenceSources.map((resourceId) => (
      resourceId ? resourceAsLegacyFile(resources, resourceId) : null
    ));
  }
  const overrides = {
    ...projectFullConfigOverrides(options.config, renderRequest.mode),
    ...(options.configOverrides || {})
  };
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
              ? normalizeCircularTrackSlot(
                  projectCanonicalCircularSlot(slot),
                  index,
                  options.dinucleotide || 'GC',
                  overrides.track_type || 'tuckin'
                )
              : parseCircularTrackSlotSpecs(
                  [slot],
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
  const projectedDepthTrackCount = Math.max(
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
  const projectedDepthTracks = Array.from({ length: projectedDepthTrackCount }, (_, index) => ({
    label: String(options.depthTrackLabels?.[index] ?? (index === 0 ? 'Depth' : `Depth ${index + 1}`)),
    color: String(options.depthTrackColors?.[index] || (index === 0 ? overrides.depth_color : '') || '#4A90E2'),
    height: optionalNumber(options.depthTrackHeights?.[index]),
    large_tick_interval: optionalNumber(options.depthTrackLargeTickIntervals?.[index]),
    small_tick_interval: optionalNumber(options.depthTrackSmallTickIntervals?.[index]),
    tick_font_size: optionalNumber(options.depthTrackTickFontSizes?.[index])
  }));
  const form = {
    prefix: renderRequest.output?.prefix || 'out',
    plot_title: options.plotTitle || '',
    legend: options.output?.legend || 'right',
    multi_record_canvas: renderRequest.mode === 'circular' && Object.keys(renderRequest.layout || {}).length > 0,
    suppress_gc: renderRequest.mode === 'circular' ? overrides.show_gc === false : false,
    suppress_skew: renderRequest.mode === 'circular' ? overrides.show_skew === false : false,
    show_gc: renderRequest.mode === 'linear' ? Boolean(overrides.show_gc) : false,
    show_skew: renderRequest.mode === 'linear' ? Boolean(overrides.show_skew) : false,
    show_depth: Boolean(overrides.show_depth),
    separate_strands: Boolean(overrides.strandedness),
    labels_mode: renderRequest.mode === 'circular' ? (overrides.allow_inner_labels ? 'both' : (overrides.show_labels ? 'out' : 'none')) : 'none',
    show_labels_linear: renderRequest.mode === 'linear' ? (overrides.show_labels || 'none') : 'none',
    track_type: overrides.track_type || 'tuckin',
    linear_track_layout: overrides.linear_track_layout || 'middle',
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
  const adv = {
    features: options.selectedFeaturesSet || [],
    feature_shapes: projectedFeatureShapes,
    nt: options.dinucleotide || 'GC',
    window_size: options.window ?? null,
    step_size: options.step ?? null,
    label_rendering: overrides.label_rendering || 'auto',
    circular_label_placement: renderRequest.mode === 'circular'
      ? (overrides.circular_label_placement || 'horizontal')
      : 'horizontal',
    label_placement: renderRequest.mode === 'linear' ? (overrides.label_placement || 'auto') : 'auto',
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
    multi_record_size_mode: renderRequest.layout?.multiRecordSizeMode || 'auto',
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
    depth_tick_interval: overrides.depth_large_tick_interval ?? null,
    depth_small_tick_interval: overrides.depth_small_tick_interval ?? null,
    depth_tick_font_size: overrides.depth_tick_font_size ?? null,
    depth_share_axis: Boolean(overrides.depth_share_axis),
    depth_window_size: options.depthWindow ?? null,
    depth_step_size: options.depthStep ?? null,
    depth_tracks: projectedDepthTracks,
    min_bitscore: options.bitscore ?? 50,
    evalue: options.evalue === null || options.evalue === undefined
      ? '1e-2'
      : String(options.evalue),
    identity: options.identity ?? 0,
    alignment_length: options.alignmentLength ?? 0,
    scale_stroke_color: overrides.scale_stroke_color ?? null,
    ruler_label_color: overrides.scale_label_color ?? null,
    scale_stroke_width: overrides.scale_stroke_width ?? null,
    scale_font_size: form.scale_style === 'ruler'
      ? (overrides.ruler_label_font_size ?? overrides.scale_font_size ?? null)
      : (overrides.scale_font_size ?? overrides.ruler_label_font_size ?? null),
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
  const linearLayout = renderRequest.mode === 'linear' && renderRequest.schema >= 2
    ? {
        enabled: Object.keys(renderRequest.layout || {}).length > 0,
        recordGap: renderRequest.layout?.recordGapPx ?? 24,
        rows: (renderRequest.layout?.multiRecordPositions || []).map((token, index) => {
          const split = String(token).lastIndexOf('@');
          return {
            uid: files.linearSeqs[index]?.uid || '',
            row: Number(String(token).slice(split + 1)) || index + 1
          };
        }),
        comparisons: (files.linearComparisons || []).map(({ file: _file, ...comparison }) => comparison)
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
    }
  };
};
