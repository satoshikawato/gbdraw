// @ts-check
/** @import { WebFontEntry } from './services/font-store.js' */
import { COMPARISON_COLOR_KEYS } from './utils/color-utils.js';
import { collectFeatureQualifierNames, specificColorQualifierSuggestions } from './services/feature-selector.js';
import { DEFAULT_LINEAR_RECORD_GAP_PX } from './services/config-field-readers.js';
import {
  featureDrawnContext,
  listFeatureRows,
  normalizeFeatureVisibilityRule
} from './services/feature-visibility.js';
import { resultCatalogFeatures } from './services/feature-catalog.js';
import {
  createDefaultLayoutPreferences,
  resolveActiveLayoutPreference,
  updateActiveLayoutPreference
} from './services/layout-preferences.js';
import { createSequenceSourceRegistry } from './services/match-sequences.js';
import {
  createDefaultLinearComparisonPlan,
  resolveLinearComparisonPlan
} from './services/linear-comparisons.js';
import { WEB_UX_PROFILE } from './web-ux-profile.js';
import {
  LEGACY_LOSAT_DERIVED_CACHE_SCHEMA,
  LEGACY_PROTEIN_CANDIDATE_SCHEMA,
  PROTEIN_IDENTITY_MANIFEST_SCHEMA
} from './contracts/schema.generated.js';
import { createDefaultFeatureRenderings } from './utils/feature-rendering.js';
import {
  getCommittedSvgContent,
  getCommittedSvgResultMetadata,
  getCommittedSvgResultRuntimeIdentity
} from './services/svg-result-ingestion.js';
import { createImportedComparisonIntentState } from './services/imported-comparison-intent.js';
import { drawnLinearSequences } from './services/record-draw-selection.js';
import { DEFAULT_SETTINGS_PANE_DISPLAY } from './services/settings-pane-state.js';
import { toRaw } from './utils/vue-raw.js';
import { computed as typedComputed } from './utils/reactivity.js';
import {
  DEFAULT_PNG_DPI,
  createDefaultAdv,
  createDefaultCircularConservation,
  createDefaultForm,
  createDefaultLosat,
  createDefaultLosatExecution
} from './services/session-active-config-contract.js';
export { createDefaultAdv, createDefaultCircularConservation, createDefaultForm, createDefaultLosat };
const { ref, reactive, computed } = window.Vue;

// System State
const processing = ref(false);
const processingStatus = ref('');
const generationCancelRequested = ref(false);
const errorLog = ref(null);
const sessionTitle = ref('');
const sessionSavePending = ref(false);
const sessionImportPending = ref(false);
const semanticFileWatchersSuppressed = ref(false);
const sessionResourceDiscoveryDeferred = ref(false);
const sessionImportRollbackInProgress = ref(false);

const results = ref([]);
const selectedResultIndex = ref(0);
const failedGeneratePreservedResult = ref(false);
const generationFailureRecovery = ref(null);
const resultPanelTab = ref('preview');
const lastRunInfo = ref(null);
const trackSlotResolvedGeometry = ref(null);
const annotationWarnings = ref([]);
// Python's report of the identity-addressed edits a Result does not draw
// (`featureIdentityNotices`, design Q4 3.4); artifact-owned like annotationWarnings.
const featureIdentityNotices = ref([]);
// The edits the last source-replacing Generate removed (Q3 = A); a notice only.
const featureEditRemovalCount = ref(0);
const comparisonWarnings = ref([]);
/** @import { ResultLabelFit } from './services/label-fit.js' */
// Python's label fit report of each drawn Result with crossing label leaders
// (LP-D46): the notice bar of the shown Result reads it. Artifact-owned like
// comparisonWarnings; not saved in a Session.
const labelFitReports = ref(/** @type {readonly ResultLabelFit[]} */ ([]));
// Store original pairwise match factors for re-interpolation
const pairwiseMatchFactors = ref({}); // { pathId: factor }
// Analysis-scoped materialized nucleotide sources used by match span popups.
const matchSequenceRegistry = createSequenceSourceRegistry();
const svgContent = computed(() => {
  if (results.value.length === 0) return null;
  return getCommittedSvgContent(results.value[selectedResultIndex.value]);
});

// Distinguish a new artifact even when its initial SVG bytes repeat a prior
// artifact whose mounted DOM has since been edited. Live edits keep this ID.
const svgResultIdentity = computed(() => (
  getCommittedSvgResultRuntimeIdentity(results.value[selectedResultIndex.value])
));

const zoom = ref(1.0);
const layoutRepositionMode = ref(false);

// Pan (drag) functionality
const isPanning = ref(false);
const panStart = reactive({ x: 0, y: 0, panX: 0, panY: 0 });
const canvasPan = reactive({ x: 0, y: 0 });
const canvasContainerRef = ref(null);

// App State
const mode = ref('circular');
const suppressCircularMultiRecordDefaults = ref(false);
const cInputType = ref('gb');
const lInputType = ref('gb');
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

/** @param {unknown} value */
const normalizeLinearSeqNumber = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) ? numeric : null;
};

/**
 * One linear record slot as `createLinearSeq` builds it. The file fields hold
 * browser `File` objects or null.
 * @typedef {object} LinearSeq
 * @property {string} uid
 * @property {any} gb
 * @property {any} gff
 * @property {any} fasta
 * @property {any} depth
 * @property {number | ''} losat_gencode `v-model.number` writes an empty string when the field is cleared
 * @property {string} definition
 * @property {string} record_subtitle
 * @property {string} file_definition
 * @property {string} file_subtitle
 * @property {string} inferred_definition
 * @property {string} region_record_id
 * @property {unknown} region_start a draft: `createLinearSeq` and the crop control write a number or null,
 *   a Session load stores the saved value as it is; the request projection validates it (R7)
 * @property {unknown} region_end a draft, as `region_start`
 * @property {boolean} region_reverse
 */

/**
 * @param {Partial<LinearSeq> | Record<string, any> | null} [overrides] a slot,
 *   or unvalidated data (a session, an import) that is normalized here
 * @returns {LinearSeq}
 */
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
    file_definition: String(source.file_definition ?? ''),
    file_subtitle: String(source.file_subtitle ?? ''),
    inferred_definition: String(source.inferred_definition ?? ''),
    region_record_id: String(source.region_record_id ?? ''),
    region_start: normalizeLinearSeqNumber(source.region_start),
    region_end: normalizeLinearSeqNumber(source.region_end),
    region_reverse: Boolean(source.region_reverse)
  };
};

/**
 * @param {readonly (Partial<LinearSeq> | Record<string, unknown> | null)[] | null | undefined} items
 *   slots, or unvalidated data (a session, an import), as `createLinearSeq` takes them
 * @returns {LinearSeq[]}
 */
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

/** @param {LinearSeq | null | undefined} seq */
const hasLinearSeqPrimaryInput = (seq) => Boolean(seq?.gb || seq?.gff || seq?.fasta);

/**
 * @param {readonly (Partial<LinearSeq> | Record<string, unknown> | null)[] | null | undefined} items
 * @returns {LinearSeq[]}
 */
export const collapseEmptyLinearSeqList = (items) => {
  const previous = normalizeLinearSeqList(items);
  if (previous.length <= 1) return previous;
  const collapsed = previous.filter((seq) => hasLinearSeqPrimaryInput(seq));
  if (collapsed.length === previous.length) return previous;
  return normalizeLinearSeqList(collapsed);
};

const linearSeqs = /** @type {LinearSeq[]} */ (reactive(normalizeLinearSeqList([])));
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
  downloadDpi: DEFAULT_PNG_DPI,
  featurePanelTab: 'colors',
  newColorFeat: 'gene',
  newColorVal: '#d3d3d3',
  newFeatureToAdd: 'mobile_element'
});

const defaultEditorDraftState = createDefaultEditorDraftState();

// App-level settings shared by both drawings: the LOSAT execution and thread
// settings (Session `ui.losatExecution`) and the rich feature popup
// (`ui.richFeaturePopup`), PD-OI-086.
const losatExecution = reactive(createDefaultLosatExecution());
const richFeaturePopup = ref(true);

const losatCacheInfo = ref([]);
const losatThreadingStatus = ref({
  state: 'unknown',
  message: ''
});
const losatCache = ref(new Map());
const losatDerivedCache = ref(new Map());
const proteinIdentityManifest = ref({
  schema: PROTEIN_IDENTITY_MANIFEST_SCHEMA,
  proteinSets: {},
  recordAnalyses: {},
  recordInstances: {}
});
const legacyProteinRawCandidates = ref({ schema: LEGACY_PROTEIN_CANDIDATE_SCHEMA, entries: [] });
const legacyProteinDerivedEvidence = ref({ schema: LEGACY_LOSAT_DERIVED_CACHE_SCHEMA, entries: [] });
const orthogroups = ref([]);
const collinearGroups = ref([]);
const featureOrthogroupIndex = ref(new Map());
const selectedOrthogroupAlignmentFeature = ref('');
const similarityAlignmentPlan = ref(null);
const similarityAlignmentResetReceipt = ref(null);
const linearRecordTranslations = ref([]);
const selectedOrthogroupId = ref('');
const orthogroupSearch = ref('');
const orthogroupSortMode = ref('id');
const showRightDrawer = ref(false);
const rightDrawerTab = ref('features'); // 'legend' | 'features' | 'orthogroups'
const linearReorderNotice = ref('');
// The Linear File cards folded to their header line, by source uid: app-level
// view state. Its transitions are in services/linear-sources.js.
/** @import { FoldedLinearSources } from './services/linear-sources.js' */
/** @type {FoldedLinearSources} */
const linearFoldedSources = reactive(new Set());
const circularRecordList = ref([]); // [{ selector: '#1', record_id: 'NC_xxx' }]
/** @import { CircularRecordDiscovery } from './services/record-discovery.js' */
/** @type {CircularRecordDiscovery} */
const circularRecordDiscovery = reactive({
  status: 'idle',
  error: '',
  inputType: '',
  primaryFile: null,
  pairedFile: null,
  canonicalRecordIdentities: []
});

// The fonts registered in this page (custom fonts D-01, services/font-store.js):
// frozen entries with their bytes, never saved (D-02) and not in History.
/** @type {WebFontEntry[]} */
const loadedFontsInitial = [];
const loadedFonts = ref(loadedFontsInitial);

// Color & Filter State
const paletteDefinitions = ref({});
const paletteNames = ref(['default']);
const defaultPaletteDraftState = createDefaultPaletteDraftState();
const paletteInstantPreviewEnabled = ref(defaultPaletteDraftState.paletteInstantPreviewEnabled);
const appliedPaletteName = ref(defaultPaletteDraftState.appliedPaletteName);
const appliedPaletteColors = ref(defaultPaletteDraftState.appliedPaletteColors);
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
/** @import { FeatureCatalog } from './services/feature-catalog.js' */
// Null, or the admitted catalog of the committed Results (`admittedFeatureCatalog`
// in services/config.js; R11), held by reference.
const featureCatalog = /** @type {{ value: FeatureCatalog | null }} */ (ref(null));
// Whether the rule preparation (app/rule-matching.js) is matching rules; it
// records the matches, which live on the feature objects, before this falls back.
const ruleMatchingPending = ref(false);
// OV-924: the qualifier names of the drawn features, read from the raw rows
// once per feature list; a rule edit or a mode switch only merges the rules.
const extractedFeatureQualifierNames = computed(() => collectFeatureQualifierNames(toRaw(extractedFeatures.value)));
const specificRuleQualifierSuggestions = computed(() =>
  specificColorQualifierSuggestions(extractedFeatureQualifierNames.value, activeDrawing().manualSpecificRules)
);
/**
 * What state.js reads of a listed feature (the Features list rows and the
 * extracted features, whose owner has no typedef).
 * @typedef {{
 *   id?: unknown, svg_id?: unknown, record_idx?: unknown, displayRecordId?: unknown,
 *   rendered_svg_id?: unknown, renderedSvgId?: unknown, rendered_feature_svg_id?: unknown,
 *   renderedFeatureSvgId?: unknown
 * }} FeatureListRow
 */
/** @param {FeatureListRow | null | undefined} feature */
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
  (Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []).forEach((/** @type {FeatureListRow} */ feature) => {
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
// Whether the Preview search bar is shown; transient view state (search-bar hide D-04).
const previewFeatureSearchVisible = ref(true);
const labelSearch = ref('');
const editableLabels = ref([]); // [{key, text, sourceText, featureId, draftText}]
const autoLabelReflowEnabled = ref(false);
const labelReflowProcessing = ref(false);
const labelReflowRequestSeq = ref(0);
const labelReflowForceRequestSeq = ref(0);
const labelReflowLastError = ref(null);

// SVG Feature Click state
const svgContainer = ref(null);
const clickedFeature = ref(null); // {id, svg_id, label, location, color, feat}
// A popup's position. A popup opened at a clicked point keeps that `anchor`
// until it is placed with its rendered size (OV-371); `maxHeight` (0: none)
// keeps a popup placed above the point from growing over it.
/** @import { PopupAnchor } from './utils/popup-bounds.js' */
/** @typedef {{ x: number, y: number, maxHeight: number, anchor: PopupAnchor | null }} PopupPosition */
const clickedFeaturePos = /** @type {PopupPosition} */ (reactive({ x: 0, y: 0, maxHeight: 0, anchor: null }));
const clickedPairwiseMatch = ref(null); // {title, subtitle, sections}
const clickedPairwiseMatchPos = /** @type {PopupPosition} */ (reactive({ x: 0, y: 0, maxHeight: 0, anchor: null }));
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
  // OV-648: whether each label scope selects other features than the Legend
  // item scope (the annotation label also than the display label scope).
  displayLabelOffered: false,
  annotationLabelOffered: false,
  existingCaptionRule: null, // Existing hash rule for same caption (already colored)
  existingCaptionColor: null, // Color of existing caption rule
  defaultColorType: null, // D-15: the feature type whose default color "Apply to all" sets
  replacedDefaultColor: null, // D-15: the user default color of that type it replaces
  resolve: null // Promise resolver
});

/**
 * The palette dialog (D-15): a palette switch or Default colors Reset while
 * user default colors exist. Display values only (`createResultsManager`).
 * @typedef {object} PaletteColorsDialog
 * @property {boolean} show
 * @property {'switch' | 'reset'} kind
 * @property {string} fromPalette
 * @property {string} toPalette
 * @property {string} droppedPalette The queued palette a switch back to the applied one drops, or ''.
 * @property {number} count The number of user default colors.
 * @property {string} keysText Their keys, for the dialog body.
 */
/** @type {PaletteColorsDialog} */
const paletteColorsDialog = reactive({
  show: false,
  kind: 'switch',
  fromPalette: '',
  toPalette: '',
  droppedPalette: '',
  count: 0,
  keysText: ''
});

// Reset Color Dialog state
const resetColorDialog = reactive({
  show: false,
  caption: '',
  siblingCount: 0
});

const legendRenameDialog = reactive({
  show: false,
  mode: 'scope', // 'scope' | 'target'
  // The key of the renamed row and of the target row ('' for a color rule's
  // name that no row shows); the dialog shows each row's drawn name.
  sourceKey: /** @type {PythonLegendKey | ''} */ (''),
  newCaption: '',
  targetKey: /** @type {PythonLegendKey | ''} */ (''),
  targetCaption: '', // the name a Merge adopts
  targetColor: '',
  currentColor: '',
  siblingCount: 0,
  mergeAvailable: true, // false unless both rows draw features of one same type
  // R15-3 (OV-285): the Python key of the deleted row the new caption names;
  // '' when the target is listed. Its Merge restores that row first.
  deletedTargetKey: /** @type {PythonLegendKey | ''} */ (''),
  // PD-OI-061 revision 4: the new caption is the name of a color rule that no
  // Legend row shows.
  targetRuleOnly: false,
  // PD-OI-061 scenario revision 5: the name the renamed row is drawn with
  // after "Keep both", and the target row's drawn name after it ('' when the
  // target is not drawn or its drawn name stays).
  keepBothName: '',
  targetNameAfter: '',
  // D-28: "Keep both" onto a deleted row in the renamed row's color restores
  // that row first, as Merge does: rows of one name and one swatch are one row.
  keepBothRestores: false,
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

// A label text edit on a feature that has no label in the displayed Result.
// `reason` is the label absence reason (label-actions.js) the dialog words.
const hiddenLabelTextDialog = reactive({
  show: false,
  featureId: '',
  reason: ''
});

// Label visibility On that the diagram cannot draw (Owner Q1, Q2): why, and
// the feature type the dialog names. The label actions own the pending edit.
const labelOnDialog = reactive({
  show: false,
  reason: '',
  featureType: ''
});

const featureVisibilityScopeDialog = reactive({
  show: false,
  feat: null,
  mode: 'default',
  previousMode: 'default',
  scopes: []
});

// Sidebar resize state
const sidebarWidth = ref(320); // Initial width in pixels
const isResizing = ref(false);
// The display of the settings pane and the entry docked beside its icon rail.
// services/settings-pane-state.js holds every transition; the browser, not
// the Session, remembers the display (app/settings-pane.js).
/** @import { SettingsPaneDisplay, SettingsPopup } from './services/settings-pane-state.js' */
/** @import { SettingsSectionKey } from './services/settings-sections.js' */
/** @type {{ value: SettingsPaneDisplay }} */
const settingsPaneDisplay = ref(DEFAULT_SETTINGS_PANE_DISPLAY);
/** @type {{ value: SettingsSectionKey | null }} */
const settingsPaneDocked = ref(null);
// The settings section a click on its drawn element moved to a popup beside the Result.
/** @type {{ value: SettingsPopup | null }} */
const settingsPopup = ref(null);

// Legend Editor state
const originalLegendOrder = ref([]); // Store original order from generation
const originalLegendColors = ref({}); // Store original colors: { caption: color }

// The Legend editor rows whose stroke options are shown, by caption: app-level
// view state, which History and the Session do not hold (OV-157).
const legendStrokeOptionsOpen = reactive(new Set());

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
const showCanvasControls = ref(false);

// Track legend position at generation time (for repositioning without regeneration)
const generatedLegendPosition = ref('left');
const generatedMode = ref('circular');
const generatedMultiRecordCanvas = ref(false);
const generatedCircularPlotTitlePosition = ref('none');

// Flag to skip captureBaseConfig when editing SVG (repositioning legend, adding legend entries, etc.)
// This prevents base config from being overwritten during incremental edits
const skipCaptureBaseConfig = ref(false);

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

const newPriorityRule = reactive(createDefaultPriorityRule());

const newFeatureToAdd = ref(defaultEditorDraftState.newFeatureToAdd);

// The committed metadata of the displayed Result (svg-result-ingestion.js).
const displayedResultMetadata = () => getCommittedSvgResultMetadata(toRaw(results.value[selectedResultIndex.value]));

// The Features drawer and Search features list the displayed Result's catalog
// features (`listFeatureRows`, R-5): a hidden feature stays listed so it can be
// shown again. The selected feature types are those of the request that drew
// the Result, which its committed metadata carries (R13). A Result
// without a catalog lists the features it renders. A batch Result shows one
// record, so the record picker appears only when the displayed Result shows
// several records (FE-03).
// The list reads of the Results only their names (the catalog admission matches
// them) and, of the displayed Result's committed metadata, its rendered
// identities and selected feature types. A write of a Result's content (an
// editor edit, R1) keeps them, and so does a Generate that reuses the catalog
// (its new metadata holds the same identities), so neither lists the features
// again (CW-06 row M3). The rule matches the drawn state reads are recorded
// while `ruleMatchingPending` holds, so the list follows its fall.
/** @import { CommittedSvgResultMetadata, SvgResultEntry } from './services/svg-result-ingestion.js' */
/** @typedef {{ list: readonly SvgResultEntry[], names: unknown[], metadata: CommittedSvgResultMetadata | null }} ListedResults */
/**
 * Whether two committed metadata give the list the same inputs: the same
 * rendered identities and the same selected feature types (none only equals none).
 * @param {CommittedSvgResultMetadata | null} previous
 * @param {CommittedSvgResultMetadata | null} metadata
 */
const sameListedMetadata = (previous, metadata) => {
  if (previous === metadata) return true;
  if (!previous || !metadata || previous.renderedFeatureIdentities !== metadata.renderedFeatureIdentities) return false;
  const before = previous.selectedFeatureTypes ?? null;
  const after = metadata.selectedFeatureTypes ?? null;
  return before === after || (before !== null && after !== null && before.length === after.length
    && before.every((type, index) => type === after[index]));
};
const listedResults = computed((/** @type {ListedResults | undefined} */ previous) => {
  /** @type {readonly SvgResultEntry[]} */
  const list = results.value;
  const names = list.map((result) => result?.name);
  const metadata = displayedResultMetadata();
  return previous !== undefined && sameListedMetadata(previous.metadata, metadata) && previous.names.length === names.length
    && previous.names.every((name, index) => name === names[index])
    ? previous : { list, names, metadata };
});
const featureList = computed(() => {
  void ruleMatchingPending.value;
  const { list, names, metadata } = listedResults.value;
  const catalogFeatures = resultCatalogFeatures(state, selectedResultIndex.value, list);
  if (catalogFeatures) {
    return listFeatureRows(catalogFeatures, featureDrawnContext(activeDrawing(), {
      diagramOptions: { selectedFeaturesSet: metadata?.selectedFeatureTypes }
    }));
  }
  const features = extractedFeatures.value;
  const renderedIds = names.length < 2 ? null : metadata?.renderedFeatureIdentities?.renderedIds;
  return {
    rows: renderedIds instanceof Set ? features.filter((/** @type {FeatureListRow} */ feature) => renderedIds.has(feature.svg_id)) : features,
    drawn: null,
    rendered: null
  };
});
// A listed feature's drawn state and the rendered ID it is drawn with ('' when
// the Result does not draw it).
/** @param {VisibilityFeature | null | undefined} feature */
const featureListState = (feature) => {
  const { drawn, rendered } = featureList.value;
  const raw = toRaw(feature);
  return {
    drawn: drawn?.get(raw) ?? true,
    renderedId: !rendered || rendered.has(raw) ? String(feature?.svg_id || '') : ''
  };
};
const featureRecordPickerVisible = computed(() => (
  featureRecordIds.value.length > 1
  && new Set(featureList.value.rows.map((/** @type {FeatureListRow} */ feature) => (
    mode.value === 'circular' ? feature.record_idx : feature.displayRecordId
  ))).size > 1
));

const filteredFeatures = computed(() => {
  let features = [...featureList.value.rows];

  // Filter by the selected record when the displayed Result shows several.
  if (featureRecordPickerVisible.value) {
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

// The height of a compact Features list row (`h-[36px]` in index.html), and
// the extra height of the one row expanded to edit (`h-[44px]`).
const FEATURE_ROW_HEIGHT_PX = 36;
const FEATURE_ROW_EXPANDED_EXTRA_PX = 44;
const FEATURE_LIST_OVERSCAN = 8;
const featureListScrollTop = ref(0);
const featureListViewportHeight = ref(520);
const featureListExpandedId = ref(null);
const isFeatureDrawerMounted = computed(
  () => Boolean(svgContent.value && showRightDrawer.value && rightDrawerTab.value === 'features')
);
const featureListExpandedIndex = computed(() => (
  featureListExpandedId.value === null
    ? -1
    : filteredFeatures.value.findIndex((/** @type {FeatureListRow} */ feature) => feature.id === featureListExpandedId.value)
));
/** @param {number} index */
const featureListRowTop = (index) => (
  index * FEATURE_ROW_HEIGHT_PX
  + (featureListExpandedIndex.value >= 0 && index > featureListExpandedIndex.value ? FEATURE_ROW_EXPANDED_EXTRA_PX : 0)
);
const featureListStartIndex = computed(() => {
  if (!isFeatureDrawerMounted.value) return 0;
  const top = featureListScrollTop.value;
  const expanded = featureListExpandedIndex.value;
  let index = Math.floor(top / FEATURE_ROW_HEIGHT_PX);
  if (expanded >= 0 && index > expanded) {
    index = Math.max(expanded, Math.floor((top - FEATURE_ROW_EXPANDED_EXTRA_PX) / FEATURE_ROW_HEIGHT_PX));
  }
  return Math.max(0, index - FEATURE_LIST_OVERSCAN);
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
const featureListTopSpacerPx = computed(() => featureListRowTop(featureListStartIndex.value));
const featureListBottomSpacerPx = computed(
  () => Math.max(0, featureListRowTop(filteredFeatures.value.length) - featureListRowTop(featureListEndIndex.value))
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

  return editableLabels.value.filter((/** @type {{ text?: unknown, sourceText?: unknown, featureId?: unknown }} */ entry) => {
    return (
      String(entry.text || '').toLowerCase().includes(query) ||
      String(entry.sourceText || '').toLowerCase().includes(query) ||
      String(entry.featureId || '').toLowerCase().includes(query)
    );
  });
});

/** @typedef {{ status: 'busy', reason: string }} SessionBusy What `sessionOperationAvailability` returns while a Session operation blocks edits. */

// Pending is owned by the existing Session lifecycle; availability stores no lock.
// 'history' words why History rejects Undo/Redo; History decides when (D-28).
// `preparationBusyReason` words an edit still applying, whose owners only the
// composition root holds: it passes it for Save and Load (R13). History words
// its own pending edit the same way, so it passes none.
/** @returns {SessionBusy | null} */
export const sessionOperationAvailability = (
  operation = 'mutation',
  preparationBusyReason = /** @type {(() => string | null | undefined) | null} */ (null)
) => {
  let reason = '';
  if (sessionImportPending.value) reason = 'Loading session. Retry after loading finishes.';
  else if (sessionSavePending.value) reason = 'Saving session. Retry after saving finishes.';
  else if (operation === 'save' || operation === 'load' || operation === 'history') {
    if (processing.value) reason = 'Generating diagram. Retry after generation finishes.';
    else if (labelReflowProcessing.value) reason = 'Updating diagram. Retry after the update finishes.';
    else reason = preparationBusyReason?.() || '';
  }
  return reason ? { status: 'busy', reason } : null;
};

/** @import { ModeProfileName } from './mode-profiles.js' */
/** @typedef {ModeProfileName} DiagramMode The diagram modes; `mode-profiles.js` owns the names. */

/** @import { PythonLegendKey, RenderedFeatureId } from './services/legend-svg.js' */
/** @import { FeatureVisibilityRuleDraft, VisibilityFeature } from './services/feature-visibility.js' */
/** @typedef {string & { readonly __brand: 'LegendCaption' }} LegendCaption The caption a Legend row shows, which a rename writes. */
/**
 * The color a Legend row records, never its identity (OV-288): one the user
 * set (a Legend color, an editor row's color, the color of a rule's row), or,
 * for a row the palette colors, the swatch the displayed Result showed when
 * its Legend list was last extracted (a Generate, a Result display, Settings
 * Reset). That can be a draft palette color, not Python's swatch, and a
 * palette, Default colors or Reset edit does not refresh it by itself. A row
 * shows `displayedLegendRowColors` (app/candidate-render.js; OV-282, D-26).
 * @typedef {string & { readonly __brand: 'RecordedLegendColor' }} RecordedLegendColor
 */
/**
 * A row of a drawing's Legend list. Its identity is Python's key
 * (`originalCaption`; an editor row's own caption); the caption and the
 * recorded color are display values their writers own (rename, Legend color),
 * never a row's identity or the features it reaches (OV-288).
 * @typedef {object} LegendEntry
 * @property {LegendCaption} caption
 * @property {PythonLegendKey} originalCaption
 * @property {RecordedLegendColor} color
 * @property {number} xPos
 * @property {number} yPos
 * @property {RenderedFeatureId[]} featureIds The features the row lists (older Sessions); empty: Python's row decides.
 */

// The drawing of a diagram mode: its settings, its editor edits, and the
// values derived from them, under their former `state` names and kinds (a ref
// stays a ref). Only the drawing holds them: a service reads them from the
// drawing it is given, and an owner resolves `state.activeDrawing()` or the
// drawing of the mode its action is about (PD-OI-086). Each mode has its own
// drawing: an edit in one mode never reaches the other, and a mode switch
// writes neither (R10). A drawing holds the complete default `form` and `adv`
// of its mode, so every reader stays valid; its `form.legend` and
// `adv.plot_title_position` read and write its own mode's layout slot.
/** @param {DiagramMode} drawingMode */
const createDrawingState = (drawingMode) => {
  const form = reactive(createDefaultForm());
  const adv = reactive(createDefaultAdv(drawingMode));
  const layoutPreferences = reactive(createDefaultLayoutPreferences());
  const activeLayoutPreferences = typedComputed(() => resolveActiveLayoutPreference(
    layoutPreferences,
    drawingMode,
    form.multi_record_canvas
  ));
  Object.defineProperty(form, 'legend', {
    enumerable: false,
    get: () => activeLayoutPreferences.value.legend,
    set: (value) => {
      updateActiveLayoutPreference(layoutPreferences, drawingMode, form.multi_record_canvas, { legend: value });
    }
  });
  Object.defineProperty(adv, 'plot_title_position', {
    enumerable: false,
    get: () => activeLayoutPreferences.value.plotTitlePosition,
    set: (value) => {
      updateActiveLayoutPreference(layoutPreferences, drawingMode, form.multi_record_canvas, { plotTitlePosition: value });
    }
  });
  const losat = reactive(createDefaultLosat());
  const losatProgram = ref('blastn'); // 'blastn' | 'tblastx' | 'blastp'
  const linearComparisonPlan = reactive(createDefaultLinearComparisonPlan());
  const linearRecordLayoutEnabled = ref(WEB_UX_PROFILE.linear.arrangeInRowsByDefault);
  const linearRecordRows = reactive([]);
  // The keys of the records this drawing leaves out (record-selection D-02):
  // Linear card uids or Circular source selectors `#N`
  // (services/record-draw-selection.js).
  /** @type {string[]} */
  const recordsOff = reactive([]);
  const linearComparisonResolution = typedComputed(() => resolveLinearComparisonPlan({
    plan: linearComparisonPlan,
    sequences: linearSeqs,
    recordsOff,
    layout: linearRecordLayoutEnabled.value ? linearRecordRows : [],
    losatProgram: losatProgram.value,
    blastpMode: losat.blastp?.mode
  }));
  /** @type {FeatureVisibilityRuleDraft[]} */
  const featureVisibilityManualRules = reactive([]);
  const pendingPaletteName = ref(defaultPaletteDraftState.pendingPaletteName);
  const labelFilters = createDefaultLabelFilterState();
  return Object.freeze({
    // Settings
    form,
    adv,
    losat,
    losatProgram,
    circularConservation: reactive(createDefaultCircularConservation()),
    linearComparisonPlan,
    linearRecordLayoutEnabled,
    linearRecordGap: ref(DEFAULT_LINEAR_RECORD_GAP_PX),
    linearRecordRows,
    recordDisplayDrafts: reactive([]),
    recordsOff,
    unmanagedConfigOverrides: reactive({}),
    importedComparisonIntent: reactive(createImportedComparisonIntentState()),
    layoutPreferences,
    linearTypographyLinked: ref(true),
    // Editor edits
    // Per-feature edits keyed by original-source identity (design Q4, R2):
    // {JSON.stringify([recordKey, biologicalFeatureId]): {recordKey,
    // biologicalFeatureId, featureVisibility, labelVisibility, labelText,
    // labelSourceText}} (services/feature-placement.js owns the row shape).
    featureOverrides: reactive({}),
    featurePlacementOverrides: reactive({}),
    featureColorOverrides: reactive({}), // {featureKey: color}
    // {featureKey: { strokeColor, strokeWidth, originalStrokeColor, originalStrokeWidth }}
    featureStrokeOverrides: reactive({}),
    featureVisibilityManualRules,
    labelTextBulkOverrides: reactive({}), // { sourceText: text }
    canonicalLabelOverrideRows: ref([]),
    legendEntries: ref([]), // LegendEntry[]; typed with its writers (OV-288 STATUS)
    deletedLegendEntries: ref([]),
    // OV-120: renamed rows the last Generate did not draw (GC off, Show Depth
    // off); their renames apply again when a Generate draws the row.
    dormantLegendEntries: ref([]),
    legendColorOverrides: reactive({}), // { caption: color }
    // { caption: { strokeColor, strokeWidth, originalStrokeColor, originalStrokeWidth } }
    legendStrokeOverrides: reactive({}),
    addedLegendCaptions: ref(new Set()),
    fileLegendCaptions: ref(new Set()),
    manualSpecificRules: reactive([]),
    manualPriorityRules: reactive([]),
    filterMode: ref(labelFilters.filterMode),
    manualBlacklist: ref(labelFilters.manualBlacklist),
    manualWhitelist: reactive(labelFilters.manualWhitelist),
    selectedPalette: ref(defaultPaletteDraftState.selectedPalette),
    currentColors: ref({ ...defaultPaletteDraftState.currentColors }),
    pendingPaletteName,
    pendingPaletteColors: ref({ ...defaultPaletteDraftState.pendingPaletteColors }),
    annotationSets: reactive([]),
    orthogroupNameOverrides: reactive({}),
    orthogroupDescriptionOverrides: reactive({}),
    // User names of groups whose members no longer form one group, keyed by member set (D-21).
    orthogroupDormantOverrides: reactive({}),
    canvasPadding: reactive({ top: 0, right: 0, bottom: 0, left: 0 }),
    // Derived from the members above
    activeLayoutPreferences,
    // The Linear record cards a Generate draws, in file order (read-only; the
    // cards are project inputs that both drawings read).
    // Frozen: a reader never edits the drawn list (the cards are `linearSeqs`).
    drawnLinearSeqs: typedComputed(() => Object.freeze(drawnLinearSequences(linearSeqs, recordsOff))),
    linearComparisonResolution,
    hasLinearComparisonIntent: typedComputed(() => linearComparisonResolution.value.hasComparisonIntent),
    hasActiveLinearLosatIntent: typedComputed(() => linearComparisonResolution.value.hasLosatIntent),
    hasActiveLinearUploadIntent: typedComputed(() => linearComparisonResolution.value.hasUploadIntent),
    // The visibility table carries rules only; per-feature edits are featureOverrides.
    featureVisibilityRules: typedComputed(() => featureVisibilityManualRules
      .map((rule) => normalizeFeatureVisibilityRule(rule))),
    hasPendingPaletteDraft: typedComputed(
      () => !paletteInstantPreviewEnabled.value && String(pendingPaletteName.value || '').trim() !== ''
    )
  });
};
const circularDrawing = createDrawingState('circular');
/** @typedef {typeof circularDrawing} DrawingState */
/** @type {Readonly<Record<DiagramMode, DrawingState>>} */
const drawings = Object.freeze({ circular: circularDrawing, linear: createDrawingState('linear') });
/** @returns {DrawingState} */
const activeDrawing = () => drawings[mode.value === 'linear' ? 'linear' : 'circular'];

export const state = Object.freeze({
  processing,
  processingStatus,
  sessionSavePending,
  sessionImportPending,
  sessionOperationAvailability,
  generationCancelRequested,
  errorLog,
  sessionTitle,
  semanticFileWatchersSuppressed,
  sessionResourceDiscoveryDeferred,
  sessionImportRollbackInProgress,
  results,
  selectedResultIndex,
  displayedResultMetadata,
  failedGeneratePreservedResult,
  generationFailureRecovery,
  resultPanelTab,
  lastRunInfo,
  trackSlotResolvedGeometry,
  annotationWarnings,
  featureIdentityNotices,
  featureEditRemovalCount,
  comparisonWarnings,
  labelFitReports,
  pairwiseMatchFactors,
  matchSequenceRegistry,
  svgContent,
  svgResultIdentity,
  zoom,
  layoutRepositionMode,
  isPanning,
  panStart,
  canvasPan,
  canvasContainerRef,
  mode,
  drawings,
  activeDrawing,
  suppressCircularMultiRecordDefaults,
  cInputType,
  lInputType,
  files,
  selectedAnnotation,
  linearSeqs,
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
  similarityAlignmentPlan,
  similarityAlignmentResetReceipt,
  linearRecordTranslations,
  selectedOrthogroupId,
  orthogroupSearch,
  orthogroupSortMode,
  showRightDrawer,
  rightDrawerTab,
  linearReorderNotice,
  linearFoldedSources,
  circularRecordList,
  circularRecordDiscovery,
  paletteDefinitions,
  paletteNames,
  paletteInstantPreviewEnabled,
  loadedFonts,
  losatExecution,
  richFeaturePopup,
  appliedPaletteName,
  appliedPaletteColors,
  newSpecRule,
  specificRulePresets,
  specificRuleQualifierSuggestions,
  selectedSpecificPreset,
  specificRulePresetLoading,
  downloadDpi,
  extractedFeatures,
  biologicalFeatures,
  featureCatalog,
  ruleMatchingPending,
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
  previewFeatureSearchVisible,
  featureListScrollTop,
  featureListViewportHeight,
  featureListExpandedId,
  isFeatureDrawerMounted,
  visibleFeatureRows,
  featureRecordPickerVisible,
  featureListTopSpacerPx,
  featureListBottomSpacerPx,
  labelSearch,
  editableLabels,
  autoLabelReflowEnabled,
  labelReflowProcessing,
  labelReflowRequestSeq,
  labelReflowForceRequestSeq,
  labelReflowLastError,
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
  paletteColorsDialog,
  resetColorDialog,
  legendRenameDialog,
  labelTextScopeDialog,
  featureVisibilityScopeDialog,
  hiddenLabelTextDialog,
  labelOnDialog,
  sidebarWidth,
  isResizing,
  settingsPaneDisplay,
  settingsPaneDocked,
  settingsPopup,
  originalLegendOrder,
  originalLegendColors,
  legendStrokeOptionsOpen,
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
  showCanvasControls,
  generatedLegendPosition,
  generatedMode,
  generatedMultiRecordCanvas,
  generatedCircularPlotTitlePosition,
  skipCaptureBaseConfig,
  skipExtractOnSvgChange,
  trustedArtifactRestoreInProgress,
  featureKeys,
  defaultColorKeys,
  newColorFeat,
  newColorVal,
  newPriorityRule,
  newFeatureToAdd,
  featureList,
  featureListState,
  filteredFeatures,
  filteredEditableLabels
});

/**
 * The Web state, for a module that receives it as a parameter. A member whose
 * declaration above has no type is `any` (Vue's `ref`, `reactive` and
 * `computed` are untyped globals).
 * @typedef {typeof state} WebState
 */
