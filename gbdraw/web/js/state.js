const { ref, reactive, computed } = window.Vue;
const DOMPurify = window.DOMPurify;

// System State
const pyodideReady = ref(false);
const processing = ref(false);
const loadingStatus = ref('Initializing...');
const errorLog = ref(null);
const sessionTitle = ref('');

const results = ref([]);
const selectedResultIndex = ref(0);
// Store original pairwise match factors for re-interpolation
const pairwiseMatchFactors = ref({}); // { pathId: factor }
const svgContent = computed(() => {
  if (results.value.length > 0) {
    const rawSvg = results.value[selectedResultIndex.value].content;

    // Sanitize the SVG output from svgwrite to ensure safety and prevent DOMXSS
    return DOMPurify.sanitize(rawSvg, {
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
        'data-legend-key',
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
  }
  return null;
});

const zoom = ref(1.0);

// Pan (drag) functionality
const isPanning = ref(false);
const panStart = reactive({ x: 0, y: 0, panX: 0, panY: 0 });
const canvasPan = reactive({ x: 0, y: 0 });
const canvasContainerRef = ref(null);

// App State
const mode = ref('circular');
const circularLegendPosition = ref('left'); // Separate legend position for circular mode
const linearLegendPosition = ref('bottom'); // Separate legend position for linear mode
const cInputType = ref('gb');
const lInputType = ref('gb');
const blastSource = ref('upload'); // 'upload' | 'losat'
const losatProgram = ref('blastn'); // 'blastn' | 'tblastx'
const files = reactive({
  c_gb: null,
  c_gff: null,
  c_fasta: null,
  d_color: null,
  t_color: null,
  blacklist: null,
  whitelist: null,
  qualifier_priority: null
});
const linearSeqs = reactive([
  {
    gb: null,
    gff: null,
    fasta: null,
    blast: null,
    losat_gencode: 1,
    losat_filename: '',
    definition: '',
    region_record_id: '',
    region_start: null,
    region_end: null,
    region_reverse: false
  }
]);

const defaultDirectionalFeatureTypes = ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'misc_RNA'];
const defaultFeatureShapes = Object.fromEntries(defaultDirectionalFeatureTypes.map((featureType) => [featureType, 'arrow']));

// Configuration Forms
const form = reactive({
  prefix: '',
  species: '',
  strain: '',
  track_type: 'tuckin',
  linear_track_layout: 'middle',
  legend: 'left',
  scale_style: 'bar',
  linear_ruler_on_axis: false,
  show_labels: false,
  show_labels_linear: 'none',
  separate_strands: true,
  allow_inner_labels: false,
  suppress_gc: false,
  suppress_skew: false,
  align_center: false,
  show_gc: false,
  show_skew: false,
  normalize_length: false
});

// Extended Advanced Config
const adv = reactive({
  features: ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'repeat_region'],
  feature_shapes: { ...defaultFeatureShapes },
  window_size: null,
  step_size: null,
  nt: 'GC',
  def_font_size: null,
  label_font_size: null,
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
  gc_height: null,
  comparison_height: null,
  min_bitscore: 50,
  evalue: '1e-2',
  identity: 0,
  scale_interval: null,
  scale_font_size: null,
  scale_stroke_width: null,
  scale_stroke_color: null,
  ruler_label_color: null,

  // Circular Specific
  feature_width_circular: null,
  gc_content_width_circular: null,
  gc_content_radius_circular: null,
  gc_skew_width_circular: null,
  gc_skew_radius_circular: null,
  outer_label_x_offset: null,
  outer_label_y_offset: null,
  inner_label_x_offset: null,
  inner_label_y_offset: null
});

const losat = reactive({
  outfmt: '6',
  blastn: {
    task: 'megablast'
  }
});

const losatCacheInfo = ref([]);
const losatCache = ref(new Map());

// Color & Filter State
const paletteNames = ref(['default']);
const selectedPalette = ref('default');
const currentColors = ref({});
const filterMode = ref('None');
const manualBlacklist = ref('hypothetical, uncharacterized, putative, unknown');
const manualWhitelist = reactive([]);
const manualSpecificRules = reactive([]);
const newSpecRule = reactive({
  feat: 'CDS',
  qual: 'product',
  val: '',
  color: '#ff0000',
  cap: ''
});
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
const selectedSpecificPreset = ref('');
const specificRulePresetLoading = ref(false);
const downloadDpi = ref(300);

// Feature Color Editor state
const extractedFeatures = ref([]); // Features from last generation
const featureRecordIds = ref([]); // Record IDs for multi-record files
const selectedFeatureRecordIdx = ref(0); // Currently selected record index
const showFeaturePanel = ref(false);
const featureSearch = ref('');
const featureColorOverrides = reactive({}); // {featureKey: color}

// SVG Feature Click state
const svgContainer = ref(null);
const clickedFeature = ref(null); // {id, svg_id, label, location, color, feat}
const clickedFeaturePos = reactive({ x: 0, y: 0 });
const featurePopupRef = ref(null);
const featurePopupDrag = reactive({ active: false, offsetX: 0, offsetY: 0 });

// Color Change Scope Dialog state
const colorScopeDialog = reactive({
  show: false,
  feat: null,
  color: null,
  matchingRule: null, // Existing regex rule from -t table
  ruleMatchCount: 0, // Number of features matching the rule
  legendName: null,
  siblingCount: 0, // Number of other features with same caption
  // For individual label option when feature belongs to multi-feature rule
  individualLabel: null, // Feature's own label (product/gene/locus_tag)
  individualLabelSiblingCount: 0, // Number of other features with same individual label
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

// Sidebar resize state
const sidebarWidth = ref(320); // Initial width in pixels
const isResizing = ref(false);

// Legend Editor state
const showLegendPanel = ref(false);
const legendEntries = ref([]); // [{caption, originalCaption, color, yPos, showStroke, featureIds}]
const deletedLegendEntries = ref([]); // Track deleted entries for restoration
const originalLegendOrder = ref([]); // Store original order from generation
const originalLegendColors = ref({}); // Store original colors: { caption: color }
const newLegendCaption = ref('');
const newLegendColor = ref('#808080');

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

// Canvas size state
const canvasPadding = reactive({ top: 0, right: 0, bottom: 0, left: 0 });
const showCanvasControls = ref(false);

// Track legend position at generation time (for repositioning without regeneration)
const generatedLegendPosition = ref('left');

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
  'transit_peptide',
  'tRNA',
  'unsure',
  'V_region',
  'V_segment',
  'variation',
  "3'UTR",
  "5'UTR"
];

const newColorFeat = ref('gene');
const newColorVal = ref('#d3d3d3');

const manualPriorityRules = reactive([]);
const newPriorityRule = reactive({ feat: 'CDS', order: 'product,gene,locus_tag' });

const newFeatureToAdd = ref('mobile_element');

const addedLegendCaptions = ref(new Set());
const fileLegendCaptions = ref(new Set());

const filteredFeatures = computed(() => {
  // First filter by feature type (only show features that are actually drawn)
  let features = extractedFeatures.value.filter((f) => adv.features.includes(f.type));

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

export const state = {
  pyodideReady,
  processing,
  loadingStatus,
  errorLog,
  sessionTitle,
  results,
  selectedResultIndex,
  pairwiseMatchFactors,
  svgContent,
  zoom,
  isPanning,
  panStart,
  canvasPan,
  canvasContainerRef,
  mode,
  circularLegendPosition,
  linearLegendPosition,
  cInputType,
  lInputType,
  blastSource,
  losatProgram,
  files,
  linearSeqs,
  form,
  adv,
  losat,
  losatCacheInfo,
  losatCache,
  paletteNames,
  selectedPalette,
  currentColors,
  filterMode,
  manualBlacklist,
  manualWhitelist,
  manualSpecificRules,
  newSpecRule,
  specificRulePresets,
  selectedSpecificPreset,
  specificRulePresetLoading,
  downloadDpi,
  extractedFeatures,
  featureRecordIds,
  selectedFeatureRecordIdx,
  showFeaturePanel,
  featureSearch,
  featureColorOverrides,
  svgContainer,
  clickedFeature,
  clickedFeaturePos,
  featurePopupRef,
  featurePopupDrag,
  colorScopeDialog,
  resetColorDialog,
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
  canvasPadding,
  showCanvasControls,
  generatedLegendPosition,
  skipCaptureBaseConfig,
  skipPositionReapply,
  skipExtractOnSvgChange,
  circularBaseConfig,
  linearBaseConfig,
  diagramElementBaseTransforms,
  featureKeys,
  newColorFeat,
  newColorVal,
  manualPriorityRules,
  newPriorityRule,
  newFeatureToAdd,
  addedLegendCaptions,
  fileLegendCaptions,
  filteredFeatures
};
