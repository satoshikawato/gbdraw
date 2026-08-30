import {
  COMPARISON_COLOR_KEYS,
  normalizePaletteColors,
  normalizePaletteDefinitions
} from './app/color-utils.js';
import { collectSpecificColorQualifierSuggestions } from './app/feature-selector.js';
import { deriveFeatureVisibilityRulesForBoundary } from './app/feature-visibility.js';
import { normalizeCircularPlotTitlePosition } from './app/plot-title-position.js';
import {
  createDefaultLayoutPreferences,
  resolveActiveLayoutPreference,
  updateActiveLayoutPreference
} from './app/layout-preferences.js';
import { createSequenceSourceRegistry } from './app/match-sequences.js';
import {
  createDefaultLinearComparisonPlan,
  resolveLinearComparisonPlan
} from './app/linear-comparisons.js';
import { createModeProfileStateManager } from './mode-profiles.js';
import { createDefaultFeatureRenderings } from './utils/feature-rendering.js';
import { getCommittedSvgContent } from './services/svg-result-ingestion.js';
import { createImportedComparisonIntentState } from './services/imported-comparison-intent.js';
import {
  createDefaultAdv,
  createDefaultCircularConservation,
  createDefaultForm,
  createDefaultLosat
} from './services/session-active-config-contract.js';
export { createDefaultAdv, createDefaultCircularConservation, createDefaultForm, createDefaultLosat };
const { ref, reactive, computed } = window.Vue;

// System State
const processing = ref(false);
const processingStatus = ref('');
const generationCancelRequested = ref(false);
const errorLog = ref(null);
const sessionTitle = ref('');
const semanticFileWatchersSuppressed = ref(false);
const sessionResourceDiscoveryDeferred = ref(false);
const sessionImportRollbackInProgress = ref(false);
const importedComparisonIntent = reactive(createImportedComparisonIntentState());

const results = ref([]);
const selectedResultIndex = ref(0);
const failedGeneratePreservedResult = ref(false);
const resultPanelTab = ref('preview');
const lastRunInfo = ref(null);
const trackSlotResolvedGeometry = ref(null);
// Store original pairwise match factors for re-interpolation
const pairwiseMatchFactors = ref({}); // { pathId: factor }
// Analysis-scoped materialized nucleotide sources used by match span popups.
const matchSequenceRegistry = createSequenceSourceRegistry();
const svgContent = computed(() => {
  if (results.value.length === 0) return null;
  return getCommittedSvgContent(results.value[selectedResultIndex.value]);
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
const layoutPreferences = reactive(createDefaultLayoutPreferences());
const suppressCircularMultiRecordDefaults = ref(false);
const cInputType = ref('gb');
const lInputType = ref('gb');
const losatProgram = ref('blastn'); // 'blastn' | 'tblastx' | 'blastp'
const files = reactive({
  c_gb: null,
  c_gff: null,
  c_fasta: null,
  c_depth: null,
  c_conservation_blasts: [],
  c_conservation_blasts_source: null,
  c_conservation_fastas: [],
  c_conservation_sequence_sources: [],
  linearCanonicalComparisons: [],
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
    losat_gencode: Number.isFinite(rawLosatGencode) && rawLosatGencode > 0 ? rawLosatGencode : 1,
    definition: String(source.definition ?? ''),
    record_subtitle: String(source.record_subtitle ?? ''),
    region_record_id: String(source.region_record_id ?? ''),
    region_start: normalizeLinearSeqNumber(source.region_start),
    region_end: normalizeLinearSeqNumber(source.region_end),
    region_reverse: Boolean(source.region_reverse)
  };
};

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
  return normalized;
};

const hasLinearSeqPrimaryInput = (seq) => Boolean(seq?.gb || seq?.gff || seq?.fasta);

export const collapseEmptyLinearSeqList = (items) => {
  const previous = normalizeLinearSeqList(items);
  if (previous.length <= 1) return previous;
  const collapsed = previous.filter((seq) => hasLinearSeqPrimaryInput(seq));
  if (collapsed.length === previous.length) return previous;
  return normalizeLinearSeqList(collapsed);
};

const linearSeqs = reactive(normalizeLinearSeqList([]));
const linearRecordLayoutEnabled = ref(false);
const linearRecordGap = ref(24);
const linearRecordRows = reactive([]);
const linearComparisonPlan = reactive(createDefaultLinearComparisonPlan());
const annotationSets = reactive([]);
const selectedAnnotation = ref(null);

export const createDefaultFeatureShapes = () => createDefaultFeatureRenderings();

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
const adv = reactive(createDefaultAdv(mode.value));
const modeProfileStateManager = createModeProfileStateManager(mode.value, adv);
const activeLayoutPreferences = computed(() => resolveActiveLayoutPreference(
  layoutPreferences,
  mode.value,
  form.multi_record_canvas
));
Object.defineProperty(form, 'legend', {
  enumerable: false,
  get: () => activeLayoutPreferences.value.legend,
  set: (value) => {
    updateActiveLayoutPreference(
      layoutPreferences,
      mode.value,
      form.multi_record_canvas,
      { legend: value }
    );
  }
});
Object.defineProperty(adv, 'plot_title_position', {
  enumerable: false,
  get: () => activeLayoutPreferences.value.plotTitlePosition,
  set: (value) => {
    updateActiveLayoutPreference(
      layoutPreferences,
      mode.value,
      form.multi_record_canvas,
      { plotTitlePosition: value }
    );
  }
});

const losat = reactive(createDefaultLosat());

const linearComparisonResolution = computed(() => resolveLinearComparisonPlan({
  plan: linearComparisonPlan,
  sequences: linearSeqs,
  layout: linearRecordLayoutEnabled.value ? linearRecordRows : [],
  losatProgram: losatProgram.value,
  blastpMode: losat.blastp?.mode
}));
const hasLinearComparisonIntent = computed(() => linearComparisonResolution.value.hasComparisonIntent);
const hasActiveLinearLosatIntent = computed(() => linearComparisonResolution.value.hasLosatIntent);
const hasActiveLinearUploadIntent = computed(() => linearComparisonResolution.value.hasUploadIntent);

const losatCacheInfo = ref([]);
const losatThreadingStatus = ref({
  state: 'unknown',
  message: ''
});
const losatCache = ref(new Map());
const losatDerivedCache = ref(new Map());
const proteinIdentityManifest = ref({
  schema: 2,
  proteinSets: {},
  recordAnalyses: {},
  recordInstances: {}
});
const legacyProteinRawCandidates = ref({ schema: 1, entries: [] });
const legacyProteinDerivedEvidence = ref({ schema: 1, entries: [] });
const orthogroups = ref([]);
const collinearGroups = ref([]);
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
const biologicalFeatures = ref([]); // Complete source catalog, including non-rendered features
const featureCatalog = ref(null); // Validated schema-3 metadata for committed Results
const specificRuleQualifierSuggestions = computed(() =>
  collectSpecificColorQualifierSuggestions(extractedFeatures.value, manualSpecificRules)
);
const featureSelectorSafetyScope = ref([]); // Python selector scope before feature visibility filtering
const renderedFeatureSvgId = (feature) => {
  return String(
    feature?.rendered_svg_id ||
    feature?.renderedSvgId ||
    feature?.rendered_feature_svg_id ||
    feature?.renderedFeatureSvgId ||
    feature?.svg_id ||
    ''
  ).trim();
};
const featuresBySvgId = computed(() => {
  const indexed = new Map();
  const features = Array.isArray(extractedFeatures.value) ? extractedFeatures.value : [];
  for (const feat of features) {
    const svgId = renderedFeatureSvgId(feat);
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
    const svgId = renderedFeatureSvgId(feature);
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
const canonicalLabelOverrideRows = ref([]);
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

// Feature style scope dialog state
const featureStyleScopeDialog = reactive({
  show: false,
  kind: 'fill', // 'fill' | 'stroke'
  feat: null,
  color: null,
  strokeColor: null,
  strokeWidth: null,
  matchingRule: null, // Existing regex rule from -t table
  ruleMatchCount: 0, // Number of features matching the rule
  legendName: null,
  siblingCount: 0, // Number of other features with same caption
  displayLabel: null, // Current rendered label text in SVG (edited label)
  displayLabelSiblingCount: 0, // Number of other features sharing display label
  annotationLabel: null, // Feature's source annotation label (product/gene/locus_tag)
  annotationLabelSiblingCount: 0, // Number of other features sharing source annotation label
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

// Trusted Artifact History restores already-validated generated owners directly.
// Watchers still mount the SVG and handlers, but must not rebuild metadata from it.
const trustedArtifactRestoreInProgress = ref(false);

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
  processing,
  processingStatus,
  generationCancelRequested,
  errorLog,
  sessionTitle,
  semanticFileWatchersSuppressed,
  sessionResourceDiscoveryDeferred,
  sessionImportRollbackInProgress,
  importedComparisonIntent,
  results,
  selectedResultIndex,
  failedGeneratePreservedResult,
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
  layoutPreferences,
  activeLayoutPreferences,
  suppressCircularMultiRecordDefaults,
  cInputType,
  lInputType,
  losatProgram,
  files,
  circularConservation,
  annotationSets,
  selectedAnnotation,
  linearSeqs,
  linearRecordLayoutEnabled,
  linearRecordGap,
  linearRecordRows,
  linearComparisonPlan,
  linearComparisonResolution,
  hasLinearComparisonIntent,
  hasActiveLinearLosatIntent,
  hasActiveLinearUploadIntent,
  form,
  adv,
  modeProfileStateManager,
  losat,
  losatCacheInfo,
  losatThreadingStatus,
  losatCache,
  losatDerivedCache,
  proteinIdentityManifest,
  legacyProteinRawCandidates,
  legacyProteinDerivedEvidence,
  orthogroups,
  collinearGroups,
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
  biologicalFeatures,
  featureCatalog,
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
  canonicalLabelOverrideRows,
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
  featureStyleScopeDialog,
  resetColorDialog,
  legendRenameDialog,
  labelTextScopeDialog,
  featureVisibilityScopeDialog,
  globalLabelModeDialog,
  sidebarWidth,
  isResizing,
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
  trustedArtifactRestoreInProgress,
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
