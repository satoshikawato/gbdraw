import {
  COMPARISON_COLOR_KEYS,
  normalizePaletteColors,
  normalizePaletteDefinitions
} from './app/color-utils.js';
import { createDefaultLinearDefinitionLineStyles } from './app/cli-args.js';
import { createDefaultCircularTrackSlots } from './app/circular-track-slots.js';
import { createDefaultLinearTrackSlots } from './app/linear-track-slots.js';
import { collectSpecificColorQualifierSuggestions } from './app/feature-selector.js';
import { deriveFeatureVisibilityRulesForBoundary } from './app/feature-visibility.js';
import { normalizeCircularPlotTitlePosition } from './app/plot-title-position.js';
import { createSequenceSourceRegistry } from './app/match-sequences.js';
const { ref, reactive, computed } = window.Vue;
const DOMPurify = window.DOMPurify;
const getNow = () => (globalThis.performance?.now ? performance.now() : Date.now());
const formatTimingMs = (ms) => `${ms.toFixed(1)}ms`;

// System State
const pyodideReady = ref(false);
const diagramGenerationWorkerReady = ref(false);
const diagramGenerationWorkerStatus = ref('Preparing diagram engine...');
const diagramGenerationWorkerError = ref(null);
const processing = ref(false);
const processingStatus = ref('');
const generationCancelRequested = ref(false);
const loadingStatus = ref('Initializing...');
const errorLog = ref(null);
const sessionTitle = ref('');

const results = ref([]);
const selectedResultIndex = ref(0);
const resultPanelTab = ref('preview');
const lastRunInfo = ref(null);
const trackSlotResolvedGeometry = ref(null);
// Store original pairwise match factors for re-interpolation
const pairwiseMatchFactors = ref({}); // { pathId: factor }
// Analysis-scoped materialized nucleotide sources used by match span popups.
const matchSequenceRegistry = createSequenceSourceRegistry();
const svgContent = computed(() => {
  if (results.value.length > 0) {
    const rawSvg = results.value[selectedResultIndex.value].content;

    // Sanitize the SVG output from svgwrite to ensure safety and prevent DOMXSS
    const sanitizeStartedAt = getNow();
    const sanitizedSvg = DOMPurify.sanitize(rawSvg, {
      USE_PROFILES: { svg: true },
      // Only allow main tags used by svgwrite
      ADD_TAGS: [
        'use',
        'g',
        'defs',
        'linearGradient',
        'radialGradient',
        'stop',
        'path',
        'rect',
        'circle',
        'line',
        'polyline',
        'polygon',
        'text',
        'tspan'
      ],
      // Allowed attributes list (event handlers are automatically excluded, but can be explicitly forbidden as well)
      ADD_ATTR: [
        'xlink:href',
        'href',
        'id',
        'class',
        'data-definition-line-kind',
        'data-gbdraw-feature-id',
        'data-legend-key',
        'data-label-key',
        'data-label-feature-id',
        'data-label-source-text',
        'data-label-editable',
        'data-collinearity-block-id',
        'data-collinearity-block-kind',
        'data-collinearity-orientation',
        'data-collinearity-block-evalue',
        'data-collinearity-color-mode',
        'data-group-kind',
        'data-group-scope',
        'data-collinear-group-scope',
        'data-orthogroup-id',
        'data-query-protein-id',
        'data-subject-protein-id',
        'data-query-feature-svg-id',
        'data-subject-feature-svg-id',
        'data-query-unit-id',
        'data-subject-unit-id',
        'data-gbdraw-match-id',
        'data-gbdraw-pairwise-match-id',
        'data-match-kind',
        'data-query-record-index',
        'data-subject-record-index',
        'data-query-record-id',
        'data-subject-record-id',
        'data-qstart',
        'data-qend',
        'data-sstart',
        'data-send',
        'data-alignment-length',
        'data-mismatches',
        'data-gap-opens',
        'data-collinearity-block-score',
        'data-collinearity-anchor-index',
        'data-collinearity-anchor-count',
        'data-query-locus-id',
        'data-subject-locus-id',
        'data-query-display-name',
        'data-subject-display-name',
        'data-pairwise-match-style',
        'data-identity-factor',
        'data-source-index',
        'data-track-index',
        'data-track-label',
        'data-track-color',
        'data-reference-side',
        'data-identity',
        'data-query',
        'data-subject',
        'data-evalue',
        'data-bitscore',
        'data-orientation',
        'data-reference-record-id',
        'data-gbdraw-annotation-id',
        'data-gbdraw-annotation-set-id',
        'data-gbdraw-annotation-track-id',
        'data-gbdraw-record-id',
        'data-gbdraw-record-index',
        'data-gbdraw-annotation-mark',
        'data-gbdraw-annotation-label',
        'data-annotation-set-id',
        'data-annotation-track-id',
        'data-annotation-record-id',
        'data-annotation-record-index',
        'data-annotation-mark',
        'data-annotation-label',
        'fill',
        'fill-opacity',
        'stroke',
        'stroke-width',
        'stroke-opacity',
        'stroke-dasharray',
        'stroke-linecap',
        'stroke-linejoin',
        'd',
        'x',
        'y',
        'x1',
        'y1',
        'x2',
        'y2',
        'cx',
        'cy',
        'r',
        'rx',
        'ry',
        'width',
        'height',
        'transform',
        'viewBox',
        'preserveAspectRatio',
        'font-family',
        'font-size',
        'font-weight',
        'text-anchor',
        'dominant-baseline',
        'writing-mode',
        'letter-spacing',
        'display'
      ],
      // Forbid <style> tags (to prevent breaking parent page CSS). Also forbid script-related tags.
      FORBID_TAGS: [
        'style',
        'script',
        'foreignObject',
        'iframe',
        'embed',
        'object',
        'animate',
        'set',
        'animateTransform',
        'image'
      ],
      FORBID_ATTR: ['name', 'onload', 'onclick', 'onmouseover', 'onfocus', 'onerror']
    });
    console.info(
      `post-gbdraw timing: DOMPurify.sanitize SVG: ${formatTimingMs(getNow() - sanitizeStartedAt)} (${String(rawSvg || '').length} chars)`
    );
    return sanitizedSvg;
  }
  return null;
});

const zoom = ref(1.0);
const layoutRepositionMode = ref(false);

// Pan (drag) functionality
const isPanning = ref(false);
const panStart = reactive({ x: 0, y: 0, panX: 0, panY: 0 });
const canvasPan = reactive({ x: 0, y: 0 });
const canvasContainerRef = ref(null);

// App State
const mode = ref('circular');
const circularLegendPosition = ref('left'); // Separate legend position for circular mode
const linearLegendPosition = ref('bottom'); // Separate legend position for linear mode
const circularPlotTitlePosition = ref('none'); // Separate plot title position for circular mode
const linearPlotTitlePosition = ref('bottom'); // Separate plot title position for linear mode
const circularSingleRecordLegendPosition = ref('left');
const circularSingleRecordPlotTitlePosition = ref('none');
const circularMultiRecordLegendPosition = ref(null);
const circularMultiRecordPlotTitlePosition = ref(null);
const suppressCircularMultiRecordDefaults = ref(false);
const cInputType = ref('gb');
const lInputType = ref('gb');
const blastSource = ref('losat'); // 'upload' | 'losat'
const losatProgram = ref('blastn'); // 'blastn' | 'tblastx' | 'blastp'
const files = reactive({
  c_gb: null,
  c_gff: null,
  c_fasta: null,
  c_depth: null,
  c_conservation_blasts: [],
  c_conservation_fastas: [],
  c_conservation_sequence_sources: [],
  d_color: null,
  t_color: null,
  blacklist: null,
  whitelist: null,
  qualifier_priority: null
});

let linearSeqUidCounter = 0;

const generateLinearSeqUid = () => {
  linearSeqUidCounter += 1;
  if (globalThis.crypto?.randomUUID) {
    return `linear-seq-${globalThis.crypto.randomUUID()}`;
  }
  return `linear-seq-${Date.now()}-${linearSeqUidCounter}`;
};

const normalizeLinearSeqNumber = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) ? numeric : null;
};

export const createLinearSeq = (overrides = {}) => {
  const source = overrides && typeof overrides === 'object' && !Array.isArray(overrides) ? overrides : {};
  const rawUid = String(source.uid ?? '').trim();
  const rawLosatGencode = Number(source.losat_gencode);
  return {
    uid: rawUid || generateLinearSeqUid(),
    gb: source.gb ?? null,
    gff: source.gff ?? null,
    fasta: source.fasta ?? null,
    depth: source.depth ?? null,
    blast: source.blast ?? null,
    losat_gencode: Number.isFinite(rawLosatGencode) && rawLosatGencode > 0 ? rawLosatGencode : 1,
    losat_filename: String(source.losat_filename ?? ''),
    definition: String(source.definition ?? ''),
    record_subtitle: String(source.record_subtitle ?? ''),
    region_record_id: String(source.region_record_id ?? ''),
    region_start: normalizeLinearSeqNumber(source.region_start),
    region_end: normalizeLinearSeqNumber(source.region_end),
    region_reverse: Boolean(source.region_reverse)
  };
};

export const clearLinearSeqGapData = (seq) => ({
  ...createLinearSeq(seq),
  blast: null,
  losat_filename: ''
});

const getLinearPairKey = (leftUid, rightUid) => `${String(leftUid || '').trim()}->${String(rightUid || '').trim()}`;

export const normalizeLinearSeqList = (items) => {
  const baseItems = Array.isArray(items) && items.length > 0 ? items : [null];
  const seenUids = new Set();
  const normalized = baseItems.map((item) => {
    const next = createLinearSeq(item);
    if (!next.uid || seenUids.has(next.uid)) {
      next.uid = generateLinearSeqUid();
    }
    seenUids.add(next.uid);
    return next;
  });
  const lastIndex = normalized.length - 1;
  return normalized.map((seq, index) => (index === lastIndex ? clearLinearSeqGapData(seq) : seq));
};

export const reconcileLinearSeqPairData = (previousItems, nextItems) => {
  const previous = normalizeLinearSeqList(previousItems);
  const next = normalizeLinearSeqList(nextItems);
  const previousPairs = new Map();

  for (let index = 0; index < previous.length - 1; index += 1) {
    const left = previous[index];
    const right = previous[index + 1];
    previousPairs.set(getLinearPairKey(left.uid, right.uid), {
      blast: left.blast ?? null,
      losatFilename: String(left.losat_filename ?? '')
    });
  }

  for (let index = 0; index < next.length; index += 1) {
    next[index] = clearLinearSeqGapData(next[index]);
  }

  const restoredPairKeys = new Set();
  for (let index = 0; index < next.length - 1; index += 1) {
    const left = next[index];
    const right = next[index + 1];
    const pairKey = getLinearPairKey(left.uid, right.uid);
    const previousPair = previousPairs.get(pairKey);
    if (!previousPair) continue;
    restoredPairKeys.add(pairKey);
    left.blast = previousPair.blast ?? null;
    left.losat_filename = String(previousPair.losatFilename ?? '');
  }

  let clearedBlastSlots = 0;
  let clearedLosatNames = 0;
  previousPairs.forEach((previousPair, pairKey) => {
    if (restoredPairKeys.has(pairKey)) return;
    if (previousPair.blast) clearedBlastSlots += 1;
    if (String(previousPair.losatFilename || '').trim()) clearedLosatNames += 1;
  });

  return {
    linearSeqs: next,
    clearedBlastSlots,
    clearedLosatNames
  };
};

const hasLinearSeqPrimaryInput = (seq) => Boolean(seq?.gb || seq?.gff || seq?.fasta);

export const collapseEmptyLinearSeqList = (items) => {
  const previous = normalizeLinearSeqList(items);
  if (previous.length <= 1) return previous;
  const collapsed = previous.filter((seq) => hasLinearSeqPrimaryInput(seq));
  if (collapsed.length === previous.length) return previous;
  return reconcileLinearSeqPairData(previous, collapsed).linearSeqs;
};

const linearSeqs = reactive(normalizeLinearSeqList([]));
const linearRecordLayoutEnabled = ref(false);
const linearRecordGap = ref(24);
const linearRecordRows = reactive([]);
const linearComparisons = reactive([]);
const annotationSets = reactive([]);
const selectedAnnotation = ref(null);

const defaultDirectionalFeatureTypes = ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'misc_RNA'];

export const createDefaultFeatureShapes = () =>
  Object.fromEntries(defaultDirectionalFeatureTypes.map((featureType) => [featureType, 'arrow']));

export const createDefaultForm = () => ({
  prefix: '',
  species: '',
  strain: '',
  plot_title: '',
  track_type: 'tuckin',
  linear_track_layout: 'middle',
  legend: 'left',
  scale_style: 'bar',
  linear_ruler_on_axis: false,
  labels_mode: 'none',
  show_labels_linear: 'none',
  multi_record_canvas: false,
  separate_strands: true,
  suppress_gc: false,
  suppress_skew: false,
  align_center: false,
  keep_definition_left_aligned: false,
  show_gc: false,
  show_skew: false,
  show_depth: false,
  normalize_length: false
});

export const createDefaultAdv = () => ({
  rich_feature_popup: true,
  features: ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'repeat_region'],
  feature_shapes: createDefaultFeatureShapes(),
  window_size: null,
  step_size: null,
  nt: 'GC',
  def_font_size: null,
  label_font_size: null,
  circular_label_spacing: null,
  linear_label_spacing: null,
  label_rendering: 'auto',
  circular_label_placement: 'horizontal',
  label_placement: 'auto',
  label_rotation: null,

  // Styles
  block_stroke_width: null,
  block_stroke_color: null,
  line_stroke_width: null,
  line_stroke_color: null,
  axis_stroke_width: null,
  axis_stroke_color: null,

  // Legend
  legend_box_size: null,
  legend_font_size: null,

  // Shared / overlap handling
  resolve_overlaps: false,

  // Linear Specific
  feature_height: null,
  track_axis_gap: null,
  linear_show_replicon: false,
  linear_show_accession: true,
  linear_show_length: true,
  linear_definition_line_styles: createDefaultLinearDefinitionLineStyles(),
  gc_height: null,
  depth_height: null,
  depth_color: '#4A90E2',
  depth_tracks: [],
  depth_window_size: null,
  depth_step_size: null,
  depth_share_axis: false,
  depth_min: null,
  depth_max: null,
  depth_normalize: false,
  depth_show_axis: true,
  depth_show_ticks: true,
  depth_tick_interval: null,
  depth_small_tick_interval: null,
  depth_tick_font_size: null,
  linear_track_slots_enabled: false,
  linear_track_slots_schema_version: 2,
  linear_track_slots_axis_index: null,
  linear_track_slots: createDefaultLinearTrackSlots(),
  gc_content_mode: 'deviation',
  gc_content_min_percent: 0,
  gc_content_max_percent: 100,
  gc_content_show_axis: true,
  gc_content_show_ticks: true,
  gc_content_tick_interval: 20,
  gc_content_small_tick_interval: null,
  gc_content_tick_font_size: null,
  comparison_height: null,
  pairwise_match_style: 'ribbon',
  min_bitscore: 50,
  evalue: '1e-2',
  identity: 0,
  alignment_length: 0,
  scale_interval: null,
  scale_font_size: null,
  scale_stroke_width: null,
  scale_stroke_color: null,
  ruler_label_color: null,

  // Circular Specific
  multi_record_size_mode: 'auto',
  multi_record_min_radius_ratio: 0.55,
  multi_record_column_gap_ratio: 0.10,
  multi_record_row_gap_ratio: 0.05,
  multi_record_positions: [],
  tick_label_font_size: null,
  plot_title_position: 'none',
  plot_title_font_size: null,
  keep_full_definition_with_plot_title: false,
  center_reserved_radius: null,
  feature_width_circular: null,
  depth_width_circular: null,
  gc_content_width_circular: null,
  gc_content_radius_circular: null,
  gc_skew_width_circular: null,
  gc_skew_radius_circular: null,
  circular_track_slots_enabled: false,
  circular_track_slots_schema_version: 4,
  circular_track_slots_axis_index: null,
  circular_track_slots: createDefaultCircularTrackSlots(),
  outer_label_x_offset: null,
  outer_label_y_offset: null,
  inner_label_x_offset: null,
  inner_label_y_offset: null
});

export const createDefaultLosat = () => ({
  outfmt: '6',
  parallelWorkers: undefined,
  executionMode: 'auto',
  totalThreadBudget: 'safe',
  threadsPerJob: 'auto',
  blastn: {
    task: 'megablast'
  },
  blastp: {
    mode: 'orthogroup',
    maxHits: 5,
    candidateLimit: null,
    orthogroupMembershipMode: 'anchor_core_v1',
    orthogroupMemberMaxHits: 5,
    collinearMinAnchors: 1,
    collinearMaxGeneGap: 0,
    collinearMaxDiagonalDrift: 0,
    collinearMaxConflictsInMergeGap: 1,
    collinearMaxParalogLinksPerOrthogroup: 2,
    collinearColorMode: 'orientation',
    collinearUnitMode: 'auto',
    collinearAnchorMode: 'rbh',
    collinearSearchScope: 'adjacent'
  }
});

export const createDefaultCircularConservation = () => ({
  enabled: false,
  source: 'losat',
  losat_program: 'blastn',
  subject_gencode: 1,
  reference: 'auto',
  labels: '',
  series: [],
  ring_width: null,
  ring_gap: null
});

export const createDefaultSpecificRule = () => ({
  feat: 'CDS',
  qual: 'product',
  val: '',
  color: '#ff0000',
  cap: ''
});

export const createDefaultPriorityRule = () => ({
  feat: 'CDS',
  order: 'product,gene,locus_tag'
});

export const createDefaultPaletteDraftState = () => ({
  selectedPalette: 'default',
  currentColors: {},
  paletteInstantPreviewEnabled: false,
  appliedPaletteName: 'default',
  appliedPaletteColors: {},
  pendingPaletteName: '',
  pendingPaletteColors: {}
});

export const createDefaultLabelFilterState = () => ({
  filterMode: 'None',
  manualBlacklist: 'hypothetical, uncharacterized, putative, unknown',
  manualWhitelist: []
});

export const createDefaultEditorDraftState = () => ({
  selectedSpecificPreset: '',
  specificRulePresetLoading: false,
  downloadDpi: 300,
  featurePanelTab: 'colors',
  newColorFeat: 'gene',
  newColorVal: '#d3d3d3',
  newFeatureToAdd: 'mobile_element',
  newLegendCaption: '',
  newLegendColor: '#808080'
});

const circularConservation = reactive(createDefaultCircularConservation());
const defaultEditorDraftState = createDefaultEditorDraftState();

// Configuration Forms
const form = reactive(createDefaultForm());

// Extended Advanced Config
const adv = reactive(createDefaultAdv());

const losat = reactive(createDefaultLosat());

const losatCacheInfo = ref([]);
const losatThreadingStatus = ref({
  state: 'unknown',
  message: ''
});
const losatCache = ref(new Map());
const losatDerivedCache = ref(new Map());
const orthogroups = ref([]);
const featureOrthogroupIndex = ref(new Map());
const selectedOrthogroupAlignmentFeature = ref('');
const orthogroupNameOverrides = reactive({});
const orthogroupDescriptionOverrides = reactive({});
const selectedOrthogroupId = ref('');
const orthogroupSearch = ref('');
const orthogroupSortMode = ref('id');
const showRightDrawer = ref(false);
const rightDrawerTab = ref('features'); // 'legend' | 'features' | 'orthogroups'
const linearReorderNotice = ref('');
const circularRecordList = ref([]); // [{ selector: '#1', record_id: 'NC_xxx' }]
const circularRecordDiscovery = reactive({
  status: 'idle',
  error: '',
  inputType: '',
  primaryFile: null,
  pairedFile: null
});

// Color & Filter State
const paletteDefinitions = ref({});
const paletteNames = ref(['default']);
const defaultPaletteDraftState = createDefaultPaletteDraftState();
const selectedPalette = ref(defaultPaletteDraftState.selectedPalette);
const currentColors = ref(defaultPaletteDraftState.currentColors);
const paletteInstantPreviewEnabled = ref(defaultPaletteDraftState.paletteInstantPreviewEnabled);
const appliedPaletteName = ref(defaultPaletteDraftState.appliedPaletteName);
const appliedPaletteColors = ref(defaultPaletteDraftState.appliedPaletteColors);
const pendingPaletteName = ref(defaultPaletteDraftState.pendingPaletteName);
const pendingPaletteColors = ref(defaultPaletteDraftState.pendingPaletteColors);
const hasPendingPaletteDraft = computed(
  () => !paletteInstantPreviewEnabled.value && String(pendingPaletteName.value || '').trim() !== ''
);
const defaultLabelFilterState = createDefaultLabelFilterState();
const filterMode = ref(defaultLabelFilterState.filterMode);
const manualBlacklist = ref(defaultLabelFilterState.manualBlacklist);
const manualWhitelist = reactive(defaultLabelFilterState.manualWhitelist);
const manualSpecificRules = reactive([]);
const newSpecRule = reactive(createDefaultSpecificRule());
const specificRulePresets = [
  {
    id: 'pharokka',
    label: 'Pharokka (function/phrog/vfdb)',
    path: 'presets/pharokka_color_table.txt'
  },
  {
    id: 'bakta',
    label: 'Bakta (COG via note)',
    path: 'presets/bakta_color_table.txt'
  }
];
const selectedSpecificPreset = ref(defaultEditorDraftState.selectedSpecificPreset);
const specificRulePresetLoading = ref(defaultEditorDraftState.specificRulePresetLoading);
const downloadDpi = ref(defaultEditorDraftState.downloadDpi);

// Feature Color Editor state
const extractedFeatures = ref([]); // Features from last generation
const specificRuleQualifierSuggestions = computed(() =>
  collectSpecificColorQualifierSuggestions(extractedFeatures.value, manualSpecificRules)
);
const featureSelectorSafetyScope = ref([]); // Python selector scope before feature visibility filtering
const featuresBySvgId = computed(() => {
  const indexed = new Map();
  const features = Array.isArray(extractedFeatures.value) ? extractedFeatures.value : [];
  for (const feat of features) {
    const svgId = String(feat?.svg_id || '').trim();
    if (!svgId || indexed.has(svgId)) continue;
    indexed.set(svgId, feat);
  }
  return indexed;
});
const selectedFeatureIds = ref(new Set());
const selectedFeatureAnchorId = ref('');
const featureSelectionStatus = ref('');
const featureSelectionSuppressNextClick = ref(false);
const featureSelectionDrag = reactive({
  active: false,
  committed: false,
  startX: 0,
  startY: 0,
  currentX: 0,
  currentY: 0,
  additive: false
});
const selectedFeatureCount = computed(() => selectedFeatureIds.value.size);
const selectedFeatures = computed(() => {
  const ids = Array.from(selectedFeatureIds.value || [])
    .map((id) => String(id || '').trim())
    .filter(Boolean);
  if (ids.length === 0) return [];

  const bySvgId = featuresBySvgId.value instanceof Map ? featuresBySvgId.value : new Map();
  const fallback = new Map();
  (Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []).forEach((feature) => {
    const svgId = String(feature?.svg_id || '').trim();
    if (svgId && !fallback.has(svgId)) fallback.set(svgId, feature);
  });

  return ids
    .map((id) => bySvgId.get(id) || fallback.get(id) || null)
    .filter(Boolean);
});
const hasFeatureSelection = computed(() => selectedFeatureCount.value > 0);
const featureEditorStatus = reactive({
  status: 'idle',
  generationId: 0,
  error: null,
  summaryCount: 0,
  detailsCacheSize: 0
});
const featureExtractionPending = ref(false);
const featureExtractionError = ref(null);
const featureRecordIds = ref([]); // Record IDs for multi-record files
const selectedFeatureRecordIdx = ref(0); // Currently selected record index
const resultGenerationKey = ref(0);
const showFeaturePanel = ref(false);
const featurePanelTab = ref(defaultEditorDraftState.featurePanelTab); // 'colors' | 'labels'
const featureSearchInput = ref('');
const featureSearch = ref('');
const previewFeatureSearchInput = ref('');
const previewFeatureSearchQuery = ref('');
const previewFeatureSearchField = ref('all');
const previewFeatureSearchQualifierKey = ref('');
const previewFeatureSearchUseRegex = ref(false);
const previewFeatureSearchMatches = ref([]);
const previewFeatureSearchMatchDetails = ref({});
const previewFeatureSearchActiveIndex = ref(-1);
const previewFeatureSearchError = ref('');
const previewFeatureSearchRenderedCount = ref(0);
const featureColorOverrides = reactive({}); // {featureKey: color}
const featureVisibilityManualRules = reactive([]);
const featureVisibilityOverrides = reactive({}); // {svg_id: 'on' | 'off' | 'exclude_matching'}
const featureVisibilitySelectorCache = reactive({});
const featureVisibilityRules = computed(() => deriveFeatureVisibilityRulesForBoundary(
  featureVisibilityManualRules,
  featureVisibilityOverrides,
  featureVisibilitySelectorCache
));
const featureStrokeOverrides = reactive({}); // {featureKey: { strokeColor, strokeWidth, originalStrokeColor, originalStrokeWidth }}
const labelSearch = ref('');
const editableLabels = ref([]); // [{key, text, sourceText, featureId, draftText}]
const labelTextFeatureOverrides = reactive({}); // { featureId: text }
const labelTextBulkOverrides = reactive({}); // { sourceText: text }
const labelTextFeatureOverrideSources = reactive({}); // { featureId: sourceText }
const labelVisibilityOverrides = reactive({}); // { featureId: 'on' | 'off' }
const labelOverrideContextKey = ref('');
const labelOverrideBuildWarning = ref('');
const autoLabelReflowEnabled = ref(false);
const labelReflowProcessing = ref(false);
const labelReflowRequestSeq = ref(0);
const labelReflowRequestReason = ref('');
const labelReflowForceRequestSeq = ref(0);
const labelReflowForceRequestReason = ref('');
const labelReflowLastError = ref(null);
const labelLayoutDirtyReason = ref('');

// SVG Feature Click state
const svgContainer = ref(null);
const clickedFeature = ref(null); // {id, svg_id, label, location, color, feat}
const clickedFeaturePos = reactive({ x: 0, y: 0 });
const clickedPairwiseMatch = ref(null); // {title, subtitle, sections}
const clickedPairwiseMatchPos = reactive({ x: 0, y: 0 });
const pairwiseMatchPopupRef = ref(null);
const pairwiseMatchPopupDrag = reactive({ active: false, offsetX: 0, offsetY: 0 });
const pairwiseMatchPopupSize = reactive({ width: 0, height: 0 });
const pairwiseMatchPopupResize = reactive({
  active: false,
  startX: 0,
  startY: 0,
  startWidth: 0,
  startHeight: 0
});
const featurePopupRef = ref(null);
const featurePopupDrag = reactive({ active: false, offsetX: 0, offsetY: 0 });
const featurePopupSize = reactive({ width: 0, height: 0 });
const featurePopupResize = reactive({
  active: false,
  startX: 0,
  startY: 0,
  startWidth: 0,
  startHeight: 0
});
const clickedLabel = ref(null); // { key, text, sourceText, featureId }
const clickedLabelPos = reactive({ x: 0, y: 0 });

// Color Change Scope Dialog state
const colorScopeDialog = reactive({
  show: false,
  feat: null,
  color: null,
  matchingRule: null, // Existing regex rule from -t table
  ruleMatchCount: 0, // Number of features matching the rule
  legendName: null,
  siblingCount: 0, // Number of other features with same caption
  displayLabel: null, // Current rendered label text in SVG (edited label)
  displayLabelSiblingCount: 0, // Number of other features sharing display label
  annotationLabel: null, // Feature's source annotation label (product/gene/locus_tag)
  annotationLabelSiblingCount: 0, // Number of other features sharing source annotation label
  // Backward-compatible alias fields for old template/method references
  individualLabel: null,
  individualLabelSiblingCount: 0,
  existingCaptionRule: null, // Existing hash rule for same caption (already colored)
  existingCaptionColor: null, // Color of existing caption rule
  resolve: null // Promise resolver
});

// Reset Color Dialog state
const resetColorDialog = reactive({
  show: false,
  caption: '',
  defaultColor: '',
  siblingCount: 0
});

const legendRenameDialog = reactive({
  show: false,
  mode: 'scope', // 'scope' | 'target'
  oldCaption: '',
  newCaption: '',
  targetCaption: '',
  targetColor: '',
  currentColor: '',
  siblingCount: 0,
  pendingRequest: null
});

// Label text scope dialog state
const labelTextScopeDialog = reactive({
  show: false,
  labelKey: '',
  newText: '',
  sourceText: '',
  featureId: '',
  matchingCount: 0
});

const featureVisibilityScopeDialog = reactive({
  show: false,
  feat: null,
  mode: 'default',
  previousMode: 'default',
  scopes: []
});

const globalLabelModeDialog = reactive({
  show: false,
  featureId: '',
  featureType: '',
  resolve: null
});

// Sidebar resize state
const sidebarWidth = ref(320); // Initial width in pixels
const isResizing = ref(false);

// Legend Editor state
const showLegendPanel = ref(false);
const legendEntries = ref([]); // [{caption, originalCaption, color, yPos, showStroke, featureIds}]
const deletedLegendEntries = ref([]); // Track deleted entries for restoration
const originalLegendOrder = ref([]); // Store original order from generation
const originalLegendColors = ref({}); // Store original colors: { caption: color }
const newLegendCaption = ref(defaultEditorDraftState.newLegendCaption);
const newLegendColor = ref(defaultEditorDraftState.newLegendColor);

// Legend stroke overrides: { caption: { strokeColor, strokeWidth, originalStrokeColor, originalStrokeWidth } }
const legendStrokeOverrides = reactive({});

// Legend color overrides: { caption: color } - tracks custom colors set via Legend Editor
const legendColorOverrides = reactive({});

// Original stroke values from SVG generation (gbdraw's auto-determined defaults)
const originalSvgStroke = ref({ color: null, width: null });

// Legend Drag state
const legendDragging = ref(false);
const legendDragStart = reactive({ x: 0, y: 0 });
const legendOriginalTransform = ref({ x: 0, y: 0 });
const legendInitialTransform = ref({ x: 0, y: 0 }); // Store SVG's original legend position
const legendCurrentOffset = reactive({ x: 0, y: 0 });

// Main Diagram Drag state (tick, labels, axis, definition, records as one group)
const diagramDragging = ref(false);
const diagramDragStart = reactive({ x: 0, y: 0 });
const diagramOffset = reactive({ x: 0, y: 0 }); // Cumulative drag offset
const diagramElementIds = ref([]); // IDs of elements that move together
const diagramElementOriginalTransforms = ref(new Map()); // Store original transforms for each element
const diagramElements = ref([]);
const lengthBarElement = ref(null);
const lengthBarOriginalTransform = ref({ x: 0, y: 0 });
const lengthBarUserOffset = reactive({ x: 0, y: 0 });
const plotTitleElement = ref(null);
const plotTitleDragging = ref(false);
const plotTitleDragStart = reactive({ x: 0, y: 0 });
const plotTitleAutoTransform = ref({ x: 0, y: 0 });
const plotTitleUserOffset = reactive({ x: 0, y: 0 });

// Canvas size state
const canvasPadding = reactive({ top: 0, right: 0, bottom: 0, left: 0 });
const showCanvasControls = ref(false);

// Track legend position at generation time (for repositioning without regeneration)
const generatedLegendPosition = ref('left');
const generatedMode = ref('circular');
const generatedMultiRecordCanvas = ref(false);
const generatedCircularPlotTitlePosition = ref('none');
const normalizeCircularLegendPosition = (value) => String(value || '').trim().toLowerCase() || 'left';
const shouldDeferCircularPreviewUpdates = computed(
  () =>
    generatedMode.value === 'circular' &&
    mode.value === 'circular' &&
    (
      Boolean(form.multi_record_canvas) !== Boolean(generatedMultiRecordCanvas.value) ||
      normalizeCircularLegendPosition(form.legend) !== normalizeCircularLegendPosition(generatedLegendPosition.value) ||
      normalizeCircularPlotTitlePosition(adv.plot_title_position) !==
        normalizeCircularPlotTitlePosition(generatedCircularPlotTitlePosition.value)
    )
);

// Flag to skip captureBaseConfig when editing SVG (repositioning legend, adding legend entries, etc.)
// This prevents base config from being overwritten during incremental edits
const skipCaptureBaseConfig = ref(false);

// Flag to skip position reapply after repositionForLegendChange is called
// This prevents infinite loop when repositionForLegendChange triggers watch(svgContent)
const skipPositionReapply = ref(false);

// Flag to skip extractLegendEntries in watch(svgContent) when setFeatureColor is handling it
// This prevents race condition where watcher overwrites correct legend state
const skipExtractOnSvgChange = ref(false);

// Base configuration for circular mode (stored when SVG is generated)
// Used to calculate absolute positions without accumulation
const circularBaseConfig = ref({
  viewBoxWidth: 0, // Base viewBox width (without legend-based adjustments)
  viewBoxHeight: 0, // Base viewBox height
  diagramCenterX: 0, // Original diagram center X coordinate
  diagramCenterY: 0, // Original diagram center Y coordinate
  legendWidth: 0, // Legend width
  legendHeight: 0, // Legend height
  generatedPosition: 'right' // The position when SVG was generated
});

// Base configuration for linear mode (stored when SVG is generated)
const linearBaseConfig = ref({
  verticalViewBox: { x: 0, y: 0, w: 0, h: 0 }, // ViewBox for vertical legend (left/right)
  horizontalViewBox: { x: 0, y: 0, w: 0, h: 0 }, // ViewBox for horizontal legend (top/bottom)
  diagramBaseTransforms: new Map(), // Base transforms for diagram elements
  horizontalLegendWidth: 0,
  horizontalLegendHeight: 0,
  verticalLegendWidth: 0,
  verticalLegendHeight: 0,
  generatedPosition: 'right' // The position when SVG was generated
});

// Store base transforms for diagram elements (separate from current offsets)
const diagramElementBaseTransforms = ref(new Map());

const featureKeys = [
  'assembly_gap',
  'C_region',
  'CDS',
  'centromere',
  'D-loop',
  'D_segment',
  'exon',
  'gap',
  'gene',
  'intron',
  'J_segment',
  'mat_peptide',
  'misc_binding',
  'misc_difference',
  'misc_feature',
  'misc_RNA',
  'misc_structure',
  'mobile_element',
  'modified_base',
  'mRNA',
  'ncRNA',
  'operon',
  'oriT',
  'precursor_RNA',
  'primer_bind',
  'propeptide',
  'protein_bind',
  'regulatory',
  'repeat_region',
  'rep_origin',
  'rRNA',
  'sig_peptide',
  'stem_loop',
  'telomere',
  'tmRNA',
  'transcript',
  'transit_peptide',
  'tRNA',
  'unsure',
  'V_region',
  'V_segment',
  'variation',
  "3'UTR",
  "5'UTR"
];

const defaultColorKeys = [...featureKeys, 'default', 'skew_high', 'skew_low', 'gc_content', ...COMPARISON_COLOR_KEYS];

const newColorFeat = ref(defaultEditorDraftState.newColorFeat);
const newColorVal = ref(defaultEditorDraftState.newColorVal);

const manualPriorityRules = reactive([]);
const newPriorityRule = reactive(createDefaultPriorityRule());

const newFeatureToAdd = ref(defaultEditorDraftState.newFeatureToAdd);

const addedLegendCaptions = ref(new Set());
const fileLegendCaptions = ref(new Set());

const filteredFeatures = computed(() => {
  let features = [...extractedFeatures.value];

  // Filter by selected record (if multiple records exist)
  if (featureRecordIds.value.length > 1) {
    const selectedIdx = selectedFeatureRecordIdx.value;
    if (mode.value === 'circular') {
      // For circular: filter by record_idx within the file
      features = features.filter((f) => f.record_idx === selectedIdx);
    } else {
      // For linear: filter by the combined file+record label
      const selectedLabel = featureRecordIds.value[selectedIdx];
      features = features.filter((f) => f.displayRecordId === selectedLabel);
    }
  }

  // Then apply search filter
  if (featureSearch.value) {
    const q = featureSearch.value.toLowerCase();
    features = features.filter(
      (f) =>
        (f.product || '').toLowerCase().includes(q) ||
        (f.gene || '').toLowerCase().includes(q) ||
        (f.locus_tag || '').toLowerCase().includes(q) ||
        (f.note || '').toLowerCase().includes(q) ||
        f.type.toLowerCase().includes(q)
    );
  }
  return features;
});

const FEATURE_ROW_HEIGHT_PX = 64;
const FEATURE_LIST_OVERSCAN = 8;
const featureListScrollTop = ref(0);
const featureListViewportHeight = ref(520);
const isFeatureDrawerMounted = computed(
  () => Boolean(svgContent.value && showRightDrawer.value && rightDrawerTab.value === 'features')
);
const featureListStartIndex = computed(() => {
  if (!isFeatureDrawerMounted.value) return 0;
  return Math.max(0, Math.floor(featureListScrollTop.value / FEATURE_ROW_HEIGHT_PX) - FEATURE_LIST_OVERSCAN);
});
const featureListEndIndex = computed(() => {
  if (!isFeatureDrawerMounted.value) return 0;
  const visibleCount = Math.ceil(featureListViewportHeight.value / FEATURE_ROW_HEIGHT_PX) + (FEATURE_LIST_OVERSCAN * 2);
  return Math.min(filteredFeatures.value.length, featureListStartIndex.value + visibleCount);
});
const visibleFeatureRows = computed(() => {
  if (!isFeatureDrawerMounted.value) return [];
  return filteredFeatures.value.slice(featureListStartIndex.value, featureListEndIndex.value);
});
const featureListTopSpacerPx = computed(() => featureListStartIndex.value * FEATURE_ROW_HEIGHT_PX);
const featureListBottomSpacerPx = computed(
  () => Math.max(0, (filteredFeatures.value.length - featureListEndIndex.value) * FEATURE_ROW_HEIGHT_PX)
);
const featureEditorStatusText = computed(() => {
  const count = Number(featureEditorStatus.summaryCount || extractedFeatures.value.length || 0);
  if (featureEditorStatus.status === 'pending-summary' || featureExtractionPending.value) {
    return count > 0 ? `${count.toLocaleString()} features updating...` : 'Preparing features...';
  }
  if (featureEditorStatus.status === 'failed' || featureExtractionError.value) {
    return 'Feature editor unavailable';
  }
  if (count > 0) return `${count.toLocaleString()} features ready`;
  return 'No features ready';
});

const filteredEditableLabels = computed(() => {
  const query = String(labelSearch.value || '').trim().toLowerCase();
  if (!query) return editableLabels.value;

  return editableLabels.value.filter((entry) => {
    return (
      String(entry.text || '').toLowerCase().includes(query) ||
      String(entry.sourceText || '').toLowerCase().includes(query) ||
      String(entry.featureId || '').toLowerCase().includes(query)
    );
  });
});

export const state = {
  pyodideReady,
  diagramGenerationWorkerReady,
  diagramGenerationWorkerStatus,
  diagramGenerationWorkerError,
  processing,
  processingStatus,
  generationCancelRequested,
  loadingStatus,
  errorLog,
  sessionTitle,
  results,
  selectedResultIndex,
  resultPanelTab,
  lastRunInfo,
  trackSlotResolvedGeometry,
  pairwiseMatchFactors,
  matchSequenceRegistry,
  svgContent,
  zoom,
  layoutRepositionMode,
  isPanning,
  panStart,
  canvasPan,
  canvasContainerRef,
  mode,
  circularLegendPosition,
  linearLegendPosition,
  circularPlotTitlePosition,
  linearPlotTitlePosition,
  circularSingleRecordLegendPosition,
  circularSingleRecordPlotTitlePosition,
  circularMultiRecordLegendPosition,
  circularMultiRecordPlotTitlePosition,
  suppressCircularMultiRecordDefaults,
  cInputType,
  lInputType,
  blastSource,
  losatProgram,
  files,
  circularConservation,
  annotationSets,
  selectedAnnotation,
  linearSeqs,
  linearRecordLayoutEnabled,
  linearRecordGap,
  linearRecordRows,
  linearComparisons,
  form,
  adv,
  losat,
  losatCacheInfo,
  losatThreadingStatus,
  losatCache,
  losatDerivedCache,
  orthogroups,
  featureOrthogroupIndex,
  selectedOrthogroupAlignmentFeature,
  orthogroupNameOverrides,
  orthogroupDescriptionOverrides,
  selectedOrthogroupId,
  orthogroupSearch,
  orthogroupSortMode,
  showRightDrawer,
  rightDrawerTab,
  linearReorderNotice,
  circularRecordList,
  circularRecordDiscovery,
  paletteDefinitions,
  paletteNames,
  selectedPalette,
  currentColors,
  paletteInstantPreviewEnabled,
  appliedPaletteName,
  appliedPaletteColors,
  pendingPaletteName,
  pendingPaletteColors,
  hasPendingPaletteDraft,
  filterMode,
  manualBlacklist,
  manualWhitelist,
  manualSpecificRules,
  newSpecRule,
  specificRulePresets,
  specificRuleQualifierSuggestions,
  selectedSpecificPreset,
  specificRulePresetLoading,
  downloadDpi,
  extractedFeatures,
  featureSelectorSafetyScope,
  featuresBySvgId,
  selectedFeatureIds,
  selectedFeatureAnchorId,
  featureSelectionStatus,
  featureSelectionSuppressNextClick,
  featureSelectionDrag,
  selectedFeatureCount,
  selectedFeatures,
  hasFeatureSelection,
  featureEditorStatus,
  featureEditorStatusText,
  featureExtractionPending,
  featureExtractionError,
  featureRecordIds,
  selectedFeatureRecordIdx,
  resultGenerationKey,
  showFeaturePanel,
  featurePanelTab,
  featureSearchInput,
  featureSearch,
  previewFeatureSearchInput,
  previewFeatureSearchQuery,
  previewFeatureSearchField,
  previewFeatureSearchQualifierKey,
  previewFeatureSearchUseRegex,
  previewFeatureSearchMatches,
  previewFeatureSearchMatchDetails,
  previewFeatureSearchActiveIndex,
  previewFeatureSearchError,
  previewFeatureSearchRenderedCount,
  featureListScrollTop,
  featureListViewportHeight,
  isFeatureDrawerMounted,
  visibleFeatureRows,
  featureListTopSpacerPx,
  featureListBottomSpacerPx,
  featureColorOverrides,
  featureVisibilityManualRules,
  featureVisibilityRules,
  featureVisibilityOverrides,
  featureVisibilitySelectorCache,
  featureStrokeOverrides,
  labelSearch,
  editableLabels,
  labelTextFeatureOverrides,
  labelTextBulkOverrides,
  labelTextFeatureOverrideSources,
  labelVisibilityOverrides,
  labelOverrideContextKey,
  labelOverrideBuildWarning,
  autoLabelReflowEnabled,
  labelReflowProcessing,
  labelReflowRequestSeq,
  labelReflowRequestReason,
  labelReflowForceRequestSeq,
  labelReflowForceRequestReason,
  labelReflowLastError,
  labelLayoutDirtyReason,
  svgContainer,
  clickedFeature,
  clickedFeaturePos,
  clickedPairwiseMatch,
  clickedPairwiseMatchPos,
  pairwiseMatchPopupRef,
  pairwiseMatchPopupDrag,
  pairwiseMatchPopupSize,
  pairwiseMatchPopupResize,
  featurePopupRef,
  featurePopupDrag,
  featurePopupSize,
  featurePopupResize,
  clickedLabel,
  clickedLabelPos,
  colorScopeDialog,
  resetColorDialog,
  legendRenameDialog,
  labelTextScopeDialog,
  featureVisibilityScopeDialog,
  globalLabelModeDialog,
  sidebarWidth,
  isResizing,
  showLegendPanel,
  legendEntries,
  deletedLegendEntries,
  originalLegendOrder,
  originalLegendColors,
  newLegendCaption,
  newLegendColor,
  legendStrokeOverrides,
  legendColorOverrides,
  originalSvgStroke,
  legendDragging,
  legendDragStart,
  legendOriginalTransform,
  legendInitialTransform,
  legendCurrentOffset,
  diagramDragging,
  diagramDragStart,
  diagramOffset,
  diagramElementIds,
  diagramElementOriginalTransforms,
  diagramElements,
  lengthBarElement,
  lengthBarOriginalTransform,
  lengthBarUserOffset,
  plotTitleElement,
  plotTitleDragging,
  plotTitleDragStart,
  plotTitleAutoTransform,
  plotTitleUserOffset,
  canvasPadding,
  showCanvasControls,
  generatedLegendPosition,
  generatedMode,
  generatedMultiRecordCanvas,
  generatedCircularPlotTitlePosition,
  shouldDeferCircularPreviewUpdates,
  skipCaptureBaseConfig,
  skipPositionReapply,
  skipExtractOnSvgChange,
  circularBaseConfig,
  linearBaseConfig,
  diagramElementBaseTransforms,
  normalizePaletteColors,
  normalizePaletteDefinitions,
  featureKeys,
  defaultColorKeys,
  newColorFeat,
  newColorVal,
  manualPriorityRules,
  newPriorityRule,
  newFeatureToAdd,
  addedLegendCaptions,
  fileLegendCaptions,
  filteredFeatures,
  filteredEditableLabels
};
