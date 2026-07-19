import { buildDefaultColorOverrideTsv } from '../app/color-utils.js';
import { serializeSpecificRules } from '../app/file-imports.js';
import { buildLabelOverrideTsv } from '../app/feature-editor/label-override-table.js';
import { serializeFeatureVisibilityRules } from '../app/feature-visibility.js';
import {
  buildCircularTrackSlotSpec,
  parseCircularTrackSlotSpecs
} from '../app/circular-track-slots.js';
import { buildLinearTrackSlotSpec } from '../app/linear-track-slots.js';
import { annotationOptionsPayload, normalizeAnnotationSets } from '../app/annotations/state.js';

export const CANONICAL_REQUEST_SCHEMA = 2;

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
  return `${resourceId}-${safe || 'resource.dat'}`;
};

const createResourceBuilder = () => {
  const resources = {};

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

  return { resources, addFile, addText };
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
    comparison_height: circular ? null : optionalNumber(adv.comparison_height),
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
  if (labelOverride.tsv) {
    diagramOptions.labelOverrideFile = fileRef(resources.addText(
      'label-override-file', 'label-override-file', 'label-overrides.tsv', labelOverride.tsv
    ));
  }
};

const buildDepthResources = ({ state, filesData, resources, diagramOptions }) => {
  const rows = state.mode.value === 'linear'
    ? (filesData.linearSeqs || []).map((seq) => Array.isArray(seq.depth) ? seq.depth : (seq.depth ? [seq.depth] : []))
    : [Array.isArray(filesData.c_depth) ? filesData.c_depth : (filesData.c_depth ? [filesData.c_depth] : [])];
  if (rows.every((row) => row.length === 0)) return;
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
  return {
    circularTrackSlots: circular && state.adv.circular_track_slots_enabled
      ? state.adv.circular_track_slots
          .filter((slot) => slot?.enabled !== false)
          .map((slot) => buildCircularTrackSlotSpec(slot, state.adv.nt, state.form.track_type))
      : null,
    circularTrackAxisIndex: circular ? optionalNumber(state.adv.circular_track_slots_axis_index) : null,
    linearTrackSlots: !circular && state.adv.linear_track_slots_enabled
      ? state.adv.linear_track_slots
          .filter((slot) => slot?.enabled !== false)
          .map((slot) => buildLinearTrackSlotSpec(slot))
      : null,
    linearTrackAxisIndex: circular ? null : optionalNumber(state.adv.linear_track_slots_axis_index),
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
    featureShapes: { ...(state.adv.feature_shapes || {}) },
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
  buildDepthResources({ state, filesData, resources, diagramOptions });

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

const resourceAsLegacyFile = (resources, resourceId) => {
  const entry = resources?.[resourceId];
  if (!entry || typeof entry !== 'object') throw new Error(`Missing canonical resource: ${resourceId}`);
  const { kind: _kind, ...file } = entry;
  return file;
};

const combineCircularGenbankResources = (resources, records) => {
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
    name: 'canonical-circular-records.gb',
    type: 'text/plain',
    size: binary.length,
    lastModified: Math.max(0, ...files.map((file) => Number(file.lastModified) || 0)),
    encoding: 'base64',
    data: btoa(binary)
  };
};

export const projectCanonicalSessionRequest = ({ renderRequest, resources, webFiles = {} }) => {
  if (!renderRequest || ![1, CANONICAL_REQUEST_SCHEMA].includes(renderRequest.schema)) {
    throw new Error('Unsupported canonical renderRequest schema.');
  }
  if (!['circular', 'linear'].includes(renderRequest.mode)) {
    throw new Error('Unsupported canonical renderRequest mode.');
  }
  const records = Array.isArray(renderRequest.records) ? renderRequest.records : [];
  if (records.length === 0) throw new Error('Canonical renderRequest records are required.');
  const files = { linearSeqs: [] };
  if (renderRequest.mode === 'circular') {
    const source = records[0]?.source || {};
    if (source.kind === 'genbank') files.c_gb = combineCircularGenbankResources(resources, records);
    if (source.kind === 'gffFasta') {
      files.c_gff = resourceAsLegacyFile(resources, source.gffResourceId);
      files.c_fasta = resourceAsLegacyFile(resources, source.fastaResourceId);
    }
  } else {
    files.linearSeqs = records.map((record, index) => {
      const source = record.source || {};
      const region = record.region || null;
      const selector = region?.selector || record.selector;
      return {
        uid: String(record.recordKey || `canonical-seq-${index + 1}`),
        gb: source.kind === 'genbank' ? resourceAsLegacyFile(resources, source.resourceId) : null,
        gff: source.kind === 'gffFasta' ? resourceAsLegacyFile(resources, source.gffResourceId) : null,
        fasta: source.kind === 'gffFasta' ? resourceAsLegacyFile(resources, source.fastaResourceId) : null,
        depth: null,
        blast: null,
        losat_gencode: 1,
        losat_filename: '',
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
  if (renderRequest.mode === 'circular' && Array.isArray(depthRows[0])) {
    const depth = depthRows[0]
      .map((ref) => ref?.resourceId ? resourceAsLegacyFile(resources, ref.resourceId) : null)
      .filter(Boolean);
    files.c_depth = depth.length > 1 ? depth : (depth[0] || null);
  }
  if (renderRequest.mode === 'linear') {
    depthRows.forEach((row, index) => {
      if (!files.linearSeqs[index] || !Array.isArray(row)) return;
      const depth = row
        .map((ref) => ref?.resourceId ? resourceAsLegacyFile(resources, ref.resourceId) : null)
        .filter(Boolean);
      files.linearSeqs[index].depth = depth.length > 1 ? depth : (depth[0] || null);
    });
  }
  const defaultColorsRef = options.colors?.defaultColorsFile || options.colors?.defaultColors;
  if (defaultColorsRef?.resourceId) {
    files.d_color = resourceAsLegacyFile(resources, defaultColorsRef.resourceId);
  }
  const colorTableRef = options.colors?.colorTableFile || options.colors?.colorTable;
  if (colorTableRef?.resourceId) {
    files.t_color = resourceAsLegacyFile(resources, colorTableRef.resourceId);
  }
  if (options.labelWhitelistFile?.resourceId) {
    files.whitelist = resourceAsLegacyFile(resources, options.labelWhitelistFile.resourceId);
  }
  const qualifierPriorityRef = options.qualifierPriorityFile || options.qualifierPriorityTable;
  if (qualifierPriorityRef?.resourceId) {
    files.qualifier_priority = resourceAsLegacyFile(resources, qualifierPriorityRef.resourceId);
  }
  if (renderRequest.mode === 'circular' && Array.isArray(options.conservationBlastFiles)) {
    files.c_conservation_blasts = options.conservationBlastFiles
      .map((ref) => ref?.resourceId ? resourceAsLegacyFile(resources, ref.resourceId) : null)
      .filter(Boolean);
  }
  if (renderRequest.mode === 'circular' && Array.isArray(webFiles.conservationSequenceSources)) {
    files.c_conservation_sequence_sources = webFiles.conservationSequenceSources.map((resourceId) => (
      resourceId ? resourceAsLegacyFile(resources, resourceId) : null
    ));
  }
  const overrides = options.configOverrides || {};
  const tracks = options.tracks || {};
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
    scale_style: overrides.scale_style || 'bar'
  };
  const adv = {
    features: options.selectedFeaturesSet || [],
    feature_shapes: options.featureShapes || {},
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
    plot_title_font_size: options.plotTitleFontSize ?? null,
    multi_record_size_mode: renderRequest.layout?.multiRecordSizeMode || 'auto',
    multi_record_min_radius_ratio: renderRequest.layout?.multiRecordMinRadiusRatio ?? 0.55,
    multi_record_column_gap_ratio: renderRequest.layout?.multiRecordColumnGapRatio ?? 0.10,
    multi_record_row_gap_ratio: renderRequest.layout?.multiRecordRowGapRatio ?? 0.05,
    center_reserved_radius: tracks.centerReservedRadius ?? null,
    circular_track_slots_enabled: renderRequest.mode === 'circular' && Array.isArray(tracks.circularTrackSlots),
    circular_track_slots_schema_version: 4,
    circular_track_slots_axis_index: tracks.circularTrackAxisIndex ?? null,
    circular_track_slots: renderRequest.mode === 'circular'
      ? parseCircularTrackSlotSpecs(tracks.circularTrackSlots, options.dinucleotide || 'GC', overrides.track_type || 'tuckin')
      : [],
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
  return {
    mode: renderRequest.mode,
    inputType: records[0]?.source?.kind === 'gffFasta' ? 'gff' : 'gb',
    files,
    config: {
      form,
      adv,
      linearRecordLayout: linearLayout,
      annotationSets: normalizeAnnotationSets(options.annotations?.sets)
    }
  };
};
