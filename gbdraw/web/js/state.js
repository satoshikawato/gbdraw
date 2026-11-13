// @ts-check
/** @import { WebFontEntry } from './services/font-store.js' */
/** @import { FileImports, ParsedSpecificRule } from './services/file-imports.js' */
/** @typedef {import('./services/drawing/inputs.js').LinearSeq} LinearSeq */
import { msg, numberText, t, tn } from './utils/i18n.js';
import { COMPARISON_COLOR_KEYS } from './utils/color-utils.js';
import { collectFeatureQualifierNames, specificColorQualifierSuggestions } from './services/feature-selector.js';
import { BUILTIN_COLOR_TABLES } from './builtin-color-tables.generated.js';
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
import { createInputFileStore } from './services/history-files.js';
import { keepFileCopy, releaseFileCopy, residentBytesOf, useInputCopyStore } from './services/file-content-cache.js';
import { createBrowserInputCopyStore } from './services/input-copies.js';
import { INPUTS_SLICE } from './services/drawing/inputs.js';
import { sessionResourceSource } from './services/session-resource-backing.js';
import {
  createDefaultLinearComparisonPlan,
  resolveLinearComparisonPlan
} from './services/linear-comparisons.js';
import { WEB_UX_PROFILE } from './web-ux-profile.js';
import {
  DEFAULT_PNG_DPI,
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
import {
  computed as typedComputed, reactive as typedReactive, ref as typedRef, shallowRef as typedShallowRef
} from './utils/reactivity.js';
import { CANVAS_SLICE } from './services/drawing/canvas.js';
import { GROUP_NAMES_SLICE } from './services/drawing/group-names.js';
import { ANNOTATION_SETS_SLICE } from './services/drawing/annotation-sets.js';
import { VISIBILITY_RULES_SLICE } from './services/drawing/visibility-rules.js';
import { COLOR_RULES_SLICE } from './services/drawing/color-rules.js';
import { PALETTE_SLICE } from './services/drawing/palette.js';
import { LABEL_TABLES_SLICE, labelTextEditsData } from './services/drawing/label-tables.js';
import { FEATURE_EDITS_SLICE } from './services/drawing/feature-edits.js';
import { PAINT_SLICE } from './services/drawing/paint.js';
import { SETTINGS_SLICE } from './services/drawing/settings.js';
import { allocateDrawingId, defaultDrawingName } from './services/drawing-view.js';
import { FEATURE_EDITOR_STATUS_IDLE } from './services/artifact-slot.js';
/** @import { FeatureEditorStatus } from './services/artifact-slot.js' */
/** @import { LosatRuntimeStatus } from './services/losat.js' */
import { diagnosticError } from './utils/error-normalization.js';
import {
  createDefaultAdv,
  createDefaultCircularConservation,
  createDefaultForm,
  createDefaultLosat,
  createDefaultLosatExecution
} from './services/session-active-config-contract.js';
export { createDefaultAdv, createDefaultCircularConservation, createDefaultForm, createDefaultLosat };
const { ref, reactive, computed } = window.Vue;

// System State
const processing = typedRef(false);
const processingStatus = ref('');
// The labels that may not fit (CAPACITY-WARN), as Python sends their values, shown under the status line until the
// run ends; the next run that shows the overlay starts without them.
const labelCapacityNotices = /** @type {{ value: readonly LabelCapacityNotice[] }} */ (ref([]));
const generationCancelRequested = ref(false);
/** @import { UserFacingError } from './utils/error-normalization.js' */
// The shown alert: the one diagnostic model the app shows (`normalizeUserFacingError`), or none.
const errorLog = typedRef(/** @type {UserFacingError | null} */ (null));
const sessionTitle = ref('');
const sessionSavePending = ref(false);
const sessionImportPending = typedRef(false);
const semanticFileWatchersSuppressed = typedRef(false);
const sessionImportRollbackInProgress = typedRef(false);

const results = typedRef(/** @type {SvgResultEntry[]} */ ([]));
const selectedResultIndex = typedRef(0);
const failedGeneratePreservedResult = ref(false);
const generationFailureRecovery = ref(null);
const resultPanelTab = typedRef('preview');
const lastRunInfo = ref(null);
const trackSlotResolvedGeometry = ref(null);
const annotationWarnings = ref([]);
// Python's report of the identity-addressed edits a Result does not draw
// (`featureIdentityNotices`, design Q4 3.4); artifact-owned like annotationWarnings.
const featureIdentityNotices = ref([]);
// The edits the last source-replacing Generate removed (Q3 = A); a notice only.
const featureEditRemovalCount = ref(0);
const comparisonWarnings = ref([]);
/** @import { LabelCapacityNotice, ResultLabelFit } from './services/label-fit.js' */
// Python's label fit report of each drawn Result with crossing label leaders
// (LP-D46): the notice bar of the shown Result reads it. Artifact-owned like
// comparisonWarnings; not saved in a Session.
const labelFitReports = ref(/** @type {readonly ResultLabelFit[]} */ ([]));
/** @import { DrawnInputSummary } from './services/input-summary.js' */
// Python's input summary of each drawn Result: the records and comparison rows
// it read and drew, and its warnings, with the Feature visibility it drew.
// Artifact-owned like labelFitReports, so History keeps it with the Result;
// not saved in a Session.
const inputSummaries = ref(/** @type {readonly DrawnInputSummary[]} */ ([]));
// Store original pairwise match factors for re-interpolation
const pairwiseMatchFactors = typedRef(/** @type {Record<string, number>} */ ({})); // { pathId: factor }
// Analysis-scoped materialized nucleotide sources used by match span popups.
const matchSequenceRegistry = createSequenceSourceRegistry();
const svgContent = typedComputed(() => {
  if (results.value.length === 0) return null;
  return getCommittedSvgContent(results.value[selectedResultIndex.value]);
});

// Distinguish a new artifact even when its initial SVG bytes repeat a prior
// artifact whose mounted DOM has since been edited. Live edits keep this ID.
const svgResultIdentity = typedComputed(() => (
  getCommittedSvgResultRuntimeIdentity(results.value[selectedResultIndex.value])
));

// The Preview's view: its zoom and pan, whether a pan is under way, and whether the figure is
// hidden or shown without its transform transition. The pan and zoom owner (app/ui.js) and the
// Load, Reset, and feature search write them; `v-preview-transform` shows them on the surface.
const zoom = typedRef(1.0);
const layoutRepositionMode = typedRef(false);

// Pan (drag) functionality: a pan writes its live position to `canvasPan` (OV-2722).
const isPanning = typedRef(false);
const panStart = typedReactive({ x: 0, y: 0, panX: 0, panY: 0 });
const canvasPan = typedReactive({ x: 0, y: 0 });
// The figure a file Session Load mounted stays hidden until the Load shows it fitted (OV-2702).
const previewHiddenUntilShown = typedRef(false);
// While it holds, the surface shows the view without the 0.2 s transition: the Fit measures the
// current transform, not a frame of its transition, and the Fit shown at once shows its view at once.
const previewTransitionHeld = typedRef(false);
const canvasContainerRef = typedRef(/** @type {HTMLElement | null} */ (null));

// App State
// The drawing a new Session starts with: one Circular drawing (PD-OI-090
// 〔drawing bar と切り替え〕); the Linear offer tab makes a Linear one.
const FIRST_DRAWING_ID = allocateDrawingId('circular', []);
// The shown drawing, by its ID: `installShownDrawing` (services/config.js) is
// its one writer. The shown mode is the shown drawing's (`drawing.mode`), so
// nothing writes `mode`.
/** @type {import('./utils/reactivity.js').Ref<DrawingId>} */
const activeDrawingId = typedRef(FIRST_DRAWING_ID);
const mode = typedComputed(() => activeDrawing().mode);

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

// The app's palette state; a drawing's palette is its palette slice's
// (services/drawing/palette.js).
export const createDefaultPaletteDraftState = () => ({
  paletteInstantPreviewEnabled: false,
  appliedPaletteName: 'default',
  appliedPaletteColors: /** @type {Record<string, string>} */ ({})
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

// App-level settings shared by every drawing: the LOSAT execution and thread
// settings (Session `ui.losatExecution`) and the rich feature popup
// (`ui.richFeaturePopup`), PD-OI-090.
const losatExecution = typedReactive(createDefaultLosatExecution());
const richFeaturePopup = ref(true);

const losatCacheInfo = ref([]);
const losatThreadingStatus = typedRef(/** @type {LosatRuntimeStatus} */ ({
  state: 'unknown',
  message: ''
}));
// The temporary copies of the chosen inputs (PD-OI-121): creating the store takes
// this tab's folder and lock and deletes the folders of closed tabs.
useInputCopyStore(createBrowserInputCopyStore());
// The project's input file store: History steps and the input bindings name
// its files by `InputFileRef` (services/history-files.js). The views of one
// Session resource share its descriptor, so their refs share a `sourceId`.
const inputFiles = createInputFileStore({
  sourceOf: (file) => /** @type {{ descriptor?: object } | null} */ (sessionResourceSource(file))?.descriptor,
  // A held File is copied when it is chosen, and its copy is deleted when no step or binding names it (PD-OI-121).
  onHold: keepFileCopy,
  onRelease: releaseFileCopy,
  // The History budget counts the bytes the page keeps in memory, not those of a File on disk (OV-2750).
  residentBytesOf
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
/** @import { OrthogroupGroup } from './services/orthogroup-feature-metadata.js' */
const orthogroups = typedRef(/** @type {OrthogroupGroup[]} */ ([]));
const collinearGroups = ref([]);
const featureOrthogroupIndex = ref(new Map());
const selectedOrthogroupAlignmentFeature = typedRef('');
/** @import { CanonicalSimilarityAlignmentPlan } from './contracts/request.generated.js' */
/** @import { StoredRecordTranslation } from './services/linear-record-layout.js' */
/** @import { AlignmentResetReceipt } from './services/session-active-config-contract.js' */
const similarityAlignmentPlan = typedRef(/** @type {CanonicalSimilarityAlignmentPlan | null} */ (null));
const similarityAlignmentResetReceipt = typedRef(/** @type {AlignmentResetReceipt | null} */ (null));
const linearRecordTranslations = typedRef(/** @type {StoredRecordTranslation[]} */ ([]));
const selectedOrthogroupId = typedRef('');
const orthogroupSearch = typedRef('');
const orthogroupSortMode = typedRef('id');
const showRightDrawer = typedRef(false);
/** @import { RightDrawerTab } from './services/right-drawer-state.js' */
const rightDrawerTab = typedRef(/** @type {RightDrawerTab} */ ('features'));
const linearReorderNotice = typedRef('');
// The Linear File cards folded to their header line, by source uid: app-level
// view state. Its transitions are in services/linear-sources.js.
/** @import { FoldedLinearSources } from './services/linear-sources.js' */
/** @type {FoldedLinearSources} */
const linearFoldedSources = reactive(new Set());

// The fonts registered in this page (custom fonts D-01, services/font-store.js):
// frozen entries with their bytes, never saved (D-02) and not in History.
const loadedFonts = typedRef(/** @type {WebFontEntry[]} */ ([]));

// Color & Filter State
const paletteDefinitions = typedRef(/** @type {Record<string, Record<string, string>>} */ ({}));
const paletteNames = typedRef(['default']);
const defaultPaletteDraftState = createDefaultPaletteDraftState();
const paletteInstantPreviewEnabled = ref(defaultPaletteDraftState.paletteInstantPreviewEnabled);
const appliedPaletteName = typedRef(defaultPaletteDraftState.appliedPaletteName);
const appliedPaletteColors = typedRef(defaultPaletteDraftState.appliedPaletteColors);
const newSpecRule = typedReactive(createDefaultSpecificRule());
/**
 * A built-in color table offered as a preset of the specific rules.
 * @typedef {{ id: string, label: string, rules: readonly Readonly<ParsedSpecificRule>[] }} SpecificRulePreset
 */
// The built-in color tables (gbdraw/data/color_tables/, read by Python), by name
// in Python's order; the CLI names the same table `-t builtin:<id>`.
/** @type {readonly SpecificRulePreset[]} */
const specificRulePresets = Object.freeze(Object.entries(BUILTIN_COLOR_TABLES).map(([id, table]) => Object.freeze({
  id,
  label: table.title,
  rules: Object.freeze(table.rows.map(([feat, qual, val, color, cap]) => Object.freeze({ feat, qual, val, color, cap })))
})));
const selectedSpecificPreset = typedRef(defaultEditorDraftState.selectedSpecificPreset);
const specificRulePresetLoading = typedRef(defaultEditorDraftState.specificRulePresetLoading);
const downloadDpi = ref(defaultEditorDraftState.downloadDpi);
// A choice of the PNG download, as the DPI beside it, kept for the page only:
// off, the PNG keeps the drawing's transparent background (as the CLI's PNG
// does); on, it is painted white. Neither a Session nor History holds it.
const downloadWhiteBackground = ref(false);

// Feature Color Editor state
/** @import { ExtractedFeature } from './services/feature-metadata-extraction.js' */
const extractedFeatures = typedRef(/** @type {ExtractedFeature[]} */ ([])); // Features from last generation
const biologicalFeatures = typedRef(/** @type {ExtractedFeature[]} */ ([])); // Complete source catalog, including non-rendered features
/** @import { FeatureCatalog } from './services/feature-catalog.js' */
// Null, or the admitted catalog of the committed Results (`admittedFeatureCatalog`
// in services/config.js; R11), held by reference.
const featureCatalog = typedRef(/** @type {FeatureCatalog | null} */ (null));
// Whether the rule preparation (app/rule-matching.js) is matching rules; it
// records the matches, which live on the feature objects, before this falls back.
const ruleMatchingPending = typedRef(false);
// OV-924: the qualifier names of the drawn features, read from the raw rows
// once per feature list; a rule edit or a mode switch only merges the rules.
const extractedFeatureQualifierNames = computed(() => collectFeatureQualifierNames(toRaw(extractedFeatures.value)));
const specificRuleQualifierSuggestions = computed(() =>
  specificColorQualifierSuggestions(extractedFeatureQualifierNames.value, activeDrawing().manualSpecificRules)
);
/** @param {ExtractedFeature | null | undefined} feature */
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
const featuresBySvgId = typedComputed(() => {
  /** @type {Map<string, ExtractedFeature>} */
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
const selectedFeatures = typedComputed(() => {
  const ids = Array.from(selectedFeatureIds.value || [])
    .map((id) => String(id || '').trim())
    .filter(Boolean);
  if (ids.length === 0) return [];

  const bySvgId = featuresBySvgId.value;
  /** @type {Map<string, ExtractedFeature>} */
  const fallback = new Map();
  (Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []).forEach((feature) => {
    const svgId = renderedFeatureSvgId(feature);
    if (svgId && !fallback.has(svgId)) fallback.set(svgId, feature);
  });

  return ids
    .map((id) => bySvgId.get(id) || fallback.get(id) || null)
    .filter((feature) => feature !== null);
});
const hasFeatureSelection = computed(() => selectedFeatureCount.value > 0);
// The shown drawing's feature-editor status: a field of its artifact slot, so a drawing switch swaps it
// with the drawing's `extractedFeatures` (OV-1764). Never saved in a Session.
const featureEditorStatus = ref(FEATURE_EDITOR_STATUS_IDLE);
/**
 * The one writer of the feature-editor status (Generate and Load).
 * @param {FeatureEditorStatus['status']} status
 * @param {string | null} [error] The failure, or the warning of a recovered Load (the Editor header's title).
 */
const setFeatureEditorStatus = (status, error = null) => {
  featureEditorStatus.value = Object.freeze({ status, error });
};
const featureExtractionPending = ref(false);
const featureRecordIds = ref([]); // Record IDs for multi-record files
const selectedFeatureRecordIdx = typedRef(0); // Currently selected record index
const resultGenerationKey = typedRef(0);
const featurePanelTab = typedRef(defaultEditorDraftState.featurePanelTab); // 'colors' | 'labels'
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
const labelSearch = typedRef('');
/** @import { EditableLabelEntry } from './services/label-override-table.js' */
// The editable labels of the displayed Result (the label actions read them from its SVG).
const editableLabels = typedRef(/** @type {EditableLabelEntry[]} */ ([]));
const autoLabelReflowEnabled = typedRef(false);
const labelReflowProcessing = typedRef(false);
const labelReflowRequestSeq = typedRef(0);
const labelReflowForceRequestSeq = typedRef(0);
// The last label rerender's failure (`liveEditFailure`), or none.
const labelReflowLastError = typedRef(
  /** @type {ReturnType<typeof import('./utils/error-normalization.js').liveEditFailure>} */ (null)
);
// The GC and skew hides an Undo took back, per drawing id: the automatic rerender
// (run-analysis.js) keeps them; its track projection and the note read them.
/** @import { UndoneTrackHide } from './services/setting-application.js' */
const undoneTrackHides = ref(/** @type {Map<string, ReadonlySet<UndoneTrackHide>>} */ (new Map()));

// SVG Feature Click state
const svgContainer = typedRef(/** @type {HTMLElement | null} */ (null));
/** @import { FeaturePopup } from './services/feature-popup.js' */
// The open feature popup (the SVG owner writes it; its fields are edited while it is open).
const clickedFeature = typedRef(/** @type {FeaturePopup | null} */ (null));
// A popup's position. A popup opened at a clicked point keeps that `anchor`
// until it is placed with its rendered size (OV-371); `maxHeight` (0: none)
// keeps a popup placed above the point from growing over it.
/** @import { PopupAnchor } from './utils/popup-bounds.js' */
/** @typedef {{ x: number, y: number, maxHeight: number, anchor: PopupAnchor | null }} PopupPosition */
const clickedFeaturePos = /** @type {PopupPosition} */ (reactive({ x: 0, y: 0, maxHeight: 0, anchor: null }));
/** @import { MatchFramePopup } from './services/match-frame-popup.js' */
// The open match-frame popup: a match (services/pairwise-match-popup.js) or a preview
// item. A deep ref: a Legend row's popup sets its `matchesReady` in place.
const clickedPairwiseMatch = typedRef(/** @type {MatchFramePopup | null} */ (null));
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

// Feature style scope dialog state (app/feature-editor/color-actions.js writes it).
/** @import { ColorFeature } from './services/feature-utils.js' */
/** @import { ColorRuleLegendGroup, SpecificColorRule } from './services/specific-color-rules.js' */
const featureStyleScopeDialog = typedReactive({
  show: false,
  kind: /** @type {'fill' | 'stroke'} */ ('fill'),
  feat: /** @type {ColorFeature | null} */ (null),
  color: /** @type {string | null} */ (null),
  strokeColor: /** @type {string | null} */ (null),
  strokeWidth: /** @type {number | null} */ (null),
  matchingRule: /** @type {Partial<SpecificColorRule> | null} */ (null), // Existing regex rule from -t table
  ruleMatchCount: 0, // Number of features matching the rule
  legendName: /** @type {string | null} */ (null),
  siblingCount: 0, // Number of other features with same caption
  displayLabel: /** @type {string | null} */ (null), // Current rendered label text in SVG (edited label)
  displayLabelSiblingCount: 0, // Number of other features sharing display label
  annotationLabel: /** @type {string | null} */ (null), // Feature's source annotation label (product/gene/locus_tag)
  annotationLabelSiblingCount: 0, // Number of other features sharing source annotation label
  // OV-648: whether each label scope selects other features than the Legend
  // item scope (the annotation label also than the display label scope).
  displayLabelOffered: false,
  annotationLabelOffered: false,
  existingCaptionColor: /** @type {string | null} */ (null), // Color of existing caption rule
  // D-15: the feature type whose default color "Apply to all" sets
  defaultColorType: /** @type {string | null} */ (null),
  // D-15: the user default color of that type it replaces
  replacedDefaultColor: /** @type {string | null} */ (null)
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
const resetColorDialog = typedReactive({
  show: false,
  caption: '',
  siblingCount: 0
});

/**
 * A rename of a Legend row from the Legend editor or the feature popup.
 * `oldCaption` is the caption the source shows; `sourceKey` is its row's
 * Python key (the Legend editor's `legendEntryKey`, the popup's effective
 * caption), by which the rules the row draws are read. A renamed row's shown
 * caption is not its key (OV-294 residual, review M1 of 2b96350f).
 * @typedef {object} LegendRenameRequest
 * @property {'popup' | 'legend'} source Where the rename was typed.
 * @property {ColorFeature} [feat] The popup's feature.
 * @property {string} oldCaption
 * @property {PythonLegendKey} sourceKey
 * @property {string} newCaption
 * @property {ColorRuleLegendGroup} [joinGroup] The group picked from the popup's list, which the rename
 *   joins: `newCaption` is its name, and its color comes with it (PD-OI-061 revision 4).
 * @property {string} [currentColor]
 * @property {number} [siblingCount]
 * @property {ColorFeature[]} [features]
 * @property {'single' | 'group' | 'manual' | null} [sourceScope]
 * @property {'merge' | 'keep-both' | null} [targetResolution]
 * @property {string} [finalCaption]
 * @property {string} [finalColor]
 * @property {PythonLegendKey} [keepBothAfter] Keep both of features split off their row: the key of the
 *   target row the new row is placed right after (D-29).
 * @property {PythonLegendKey} [joins] The key of the listed row, not the renamed one, whose name the rename
 *   takes (a join or a Merge): that row keeps its place, so the new rules are not put in the renamed row's.
 */
const legendRenameDialog = typedReactive({
  show: false,
  mode: /** @type {'scope' | 'target'} */ ('scope'),
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
  // Why the dialog asks: its one sentence (`legendNameConflictReason`, OV-649).
  reason: '',
  // PD-OI-061 scenario revision 5: the name the renamed row is drawn with
  // after "Keep both", and the target row's drawn name after it ('' when the
  // target is not drawn or its drawn name stays).
  keepBothName: '',
  targetNameAfter: '',
  // D-28: "Keep both" onto a deleted row in the renamed row's color restores
  // that row first, as Merge does: rows of one name and one swatch are one row.
  keepBothRestores: false,
  pendingRequest: /** @type {LegendRenameRequest | null} */ (null)
});

// Label text scope dialog state
const labelTextScopeDialog = typedReactive({
  show: false,
  labelKey: '',
  newText: '',
  sourceText: '',
  featureId: '',
  matchingCount: 0
});

// A label text edit on a feature that has no label in the displayed Result.
// `reason` is the label absence reason (label-actions.js) the dialog words.
const hiddenLabelTextDialog = typedReactive({
  show: false,
  featureId: '',
  reason: ''
});

// Label visibility On that the diagram cannot draw (Owner Q1, Q2): why, and
// the feature type the dialog names. The label actions own the pending edit.
const labelOnDialog = typedReactive({
  show: false,
  reason: '',
  featureType: ''
});

/**
 * The Change Scope dialog of a feature visibility edit
 * (`app/feature-editor/visibility-actions.js`). Display values only.
 * @typedef {object} FeatureVisibilityScopeDialog
 * @property {boolean} show
 * @property {VisibilityTarget | null} feat The feature the edit is for.
 * @property {string} mode The visibility mode the edit sets.
 * @property {string} previousMode The feature's mode before it, which Cancel shows again.
 * @property {VisibilityScope[]} scopes The scopes the dialog offers.
 */
const featureVisibilityScopeDialog = typedReactive(/** @type {FeatureVisibilityScopeDialog} */ ({
  show: false,
  feat: null,
  mode: 'default',
  previousMode: 'default',
  scopes: []
}));

// Sidebar resize state
const sidebarWidth = typedRef(320); // Initial width in pixels
const isResizing = typedRef(false);
// The display of the settings pane and the entry docked beside its icon rail.
// services/settings-pane-state.js holds every transition; the browser, not
// the Session, remembers the display (app/settings-pane.js).
/** @import { SettingsPaneDisplay, SettingsPopup } from './services/settings-pane-state.js' */
/** @import { SettingsSectionKey } from './services/settings-sections.js' */
const settingsPaneDisplay = typedRef(/** @type {SettingsPaneDisplay} */ (DEFAULT_SETTINGS_PANE_DISPLAY));
const settingsPaneDocked = typedRef(/** @type {SettingsSectionKey | null} */ (null));
// The settings section a click on its drawn element moved to a popup beside the Result.
const settingsPopup = typedRef(/** @type {SettingsPopup | null} */ (null));

// Legend Editor state
// The shown Result's generated Legend inventory: Python's keys of its rows, in the order it drew them.
const originalLegendOrder = typedRef(/** @type {PythonLegendKey[]} */ ([]));
// The shown Result's generated Legend colors, by Python's key of each row.
const originalLegendColors = typedRef(/** @type {Record<string, string>} */ ({}));

// The Legend editor rows whose stroke options are shown, by the row's key: app-level
// view state, which History and the Session do not hold (OV-157).
const legendStrokeOptionsOpen = typedReactive(/** @type {Set<PythonLegendKey>} */ (new Set()));

// Original stroke values from SVG generation (gbdraw's auto-determined defaults)
const originalSvgStroke = ref({ color: null, width: null });

// Legend Drag state
const legendDragging = typedRef(false);
const legendDragStart = typedReactive({ x: 0, y: 0 });
const legendOriginalTransform = typedRef({ x: 0, y: 0 });
const legendInitialTransform = typedRef({ x: 0, y: 0 }); // Store SVG's original legend position
const legendCurrentOffset = typedReactive({ x: 0, y: 0 });

// Main Diagram Drag state (tick, labels, axis, definition, records as one group)
const diagramDragging = typedRef(false);
const diagramDragStart = typedReactive({ x: 0, y: 0 });
const diagramOffset = typedReactive({ x: 0, y: 0 }); // Cumulative drag offset
const diagramElementIds = typedRef(/** @type {string[]} */ ([])); // IDs of elements that move together
// The automatic translation of each bound primary group, by group.
const diagramElementOriginalTransforms = typedRef(
  /** @type {Map<SVGGraphicsElement, { x: number, y: number }>} */ (new Map())
);
const diagramElements = typedRef(/** @type {SVGGraphicsElement[]} */ ([]));
const lengthBarElement = typedRef(/** @type {SVGGraphicsElement | null} */ (null));
const lengthBarOriginalTransform = typedRef({ x: 0, y: 0 });
const lengthBarUserOffset = typedReactive({ x: 0, y: 0 });
const plotTitleElement = typedRef(/** @type {SVGGraphicsElement | null} */ (null));
const plotTitleDragging = typedRef(false);
const plotTitleDragStart = typedReactive({ x: 0, y: 0 });
const plotTitleAutoTransform = typedRef({ x: 0, y: 0 });
const plotTitleUserOffset = typedReactive({ x: 0, y: 0 });

// Canvas size state
const showCanvasControls = ref(false);

// Track legend position at generation time (for repositioning without regeneration)
const generatedLegendPosition = typedRef('left');
// The drawing the displayed Result belongs to, by ID (S6): the artifact slot
// installs write it (`installArtifactSlot`, the History owner set), and Load.
/** @type {import('./utils/reactivity.js').Ref<DrawingId>} */
const generatedDrawingId = typedRef(FIRST_DRAWING_ID);
const generatedMultiRecordCanvas = ref(false);
const generatedCircularPlotTitlePosition = ref('none');

// Flag to skip captureBaseConfig when editing SVG (repositioning legend, adding legend entries, etc.)
// This prevents base config from being overwritten during incremental edits
const skipCaptureBaseConfig = typedRef(false);

// Flag to skip extractLegendEntries in watch(svgContent) when setFeatureColor is handling it
// This prevents race condition where watcher overwrites correct legend state
const skipExtractOnSvgChange = typedRef(false);

// Trusted Artifact History restores already-validated generated owners directly.
// Watchers still mount the SVG and handlers, but must not rebuild metadata from it.
const trustedArtifactRestoreInProgress = typedRef(false);

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

const newColorFeat = typedRef(defaultEditorDraftState.newColorFeat);
const newColorVal = typedRef(defaultEditorDraftState.newColorVal);

const newPriorityRule = typedReactive(createDefaultPriorityRule());

const newFeatureToAdd = typedRef(defaultEditorDraftState.newFeatureToAdd);

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
  const drawing = activeDrawing();
  const catalogFeatures = resultCatalogFeatures(state, drawing.mode, selectedResultIndex.value, list);
  if (catalogFeatures) {
    return listFeatureRows(catalogFeatures, featureDrawnContext(drawing, {
      diagramOptions: { selectedFeaturesSet: metadata?.selectedFeatureTypes }
    }));
  }
  const features = extractedFeatures.value;
  const renderedIds = names.length < 2 ? null : metadata?.renderedFeatureIdentities?.renderedIds;
  return {
    rows: renderedIds instanceof Set ? features.filter((feature) => renderedIds.has(feature.svg_id)) : features,
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
  && new Set(featureList.value.rows.map((/** @type {ExtractedFeature} */ feature) => (
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
const isFeatureDrawerMounted = typedComputed(
  () => Boolean(svgContent.value && showRightDrawer.value && rightDrawerTab.value === 'features')
);
const featureListExpandedIndex = computed(() => (
  featureListExpandedId.value === null
    ? -1
    : filteredFeatures.value.findIndex((/** @type {ExtractedFeature} */ feature) => feature.id === featureListExpandedId.value)
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
// The count is the shown drawing's: its `extractedFeatures` switch with it (the artifact slot, OV-1762).
const featureEditorStatusText = computed(() => {
  const count = extractedFeatures.value.length;
  if (featureExtractionPending.value) {
    return count > 0 ? tn(count, '{count} feature updating...', '{count} features updating...', { count: numberText(count) }) : t('Preparing features...');
  }
  if (featureEditorStatus.value.status === 'failed') {
    return t('Feature editor unavailable');
  }
  if (count > 0) return tn(count, '{count} feature ready', '{count} features ready', { count: numberText(count) });
  return t('No features ready');
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

/**
 * What `sessionOperationAvailability` returns while a Session operation blocks
 * edits. `reason` is English marked with `msg`; show it with `tm`.
 * @typedef {{ status: 'busy', reason: string }} SessionBusy
 */

// Pending is owned by the existing Session lifecycle; availability stores no lock.
// 'history' words why History rejects Undo/Redo; History decides when (D-28).
// `preparationBusyReason` words an edit still applying, whose owners only the
// composition root holds: it passes it for Save and Load (R13). History words
// its own pending edit the same way, so it passes none. 'not-during-run' also
// refuses while a run is processing (a Generate, or a Rotate, Align, or Align at
// feature run): Reset Settings derives the ring series again, from which the
// running Generate built its ring searches, and those runs replace the Result
// that Generate adopts (OV-2446).
/** @returns {SessionBusy | null} */
export const sessionOperationAvailability = (
  operation = 'mutation',
  preparationBusyReason = /** @type {(() => string | null | undefined) | null} */ (null)
) => {
  let reason = '';
  const generating = processing.value ? msg('Generating diagram. Retry after generation finishes.') : '';
  if (sessionImportPending.value) reason = msg('Loading session. Retry after loading finishes.');
  else if (sessionSavePending.value) reason = msg('Saving session. Retry after saving finishes.');
  else if (operation === 'not-during-run') reason = generating;
  else if (operation === 'save' || operation === 'load' || operation === 'history') {
    if (generating) reason = generating;
    else if (labelReflowProcessing.value) reason = msg('Updating diagram. Retry after the update finishes.');
    else reason = preparationBusyReason?.() || '';
  }
  return reason ? { status: 'busy', reason } : null;
};

/** @import { ModeProfileName } from './mode-profiles.js' */
/** @typedef {ModeProfileName} DiagramMode The diagram modes; `mode-profiles.js` owns the names. */
/** @import { DrawingId } from './services/drawing-view.js' */

/** @import { PythonLegendKey } from './services/legend-svg.js' */
/** @import { VisibilityFeature, VisibilityScope, VisibilityTarget } from './services/feature-visibility.js' */

// The drawing of a diagram mode: its settings, its editor edits, and the
// values derived from them, under their former `state` names and kinds (a ref
// stays a ref). Only the drawing holds them: a service reads them from the
// drawing it is given, and an owner resolves `state.activeDrawing()` or the
// drawing its action is about (PD-OI-090). Each drawing holds its own: an
// edit in one drawing never reaches another, and a drawing switch writes
// none (R10). A drawing holds the complete default `form` and `adv`
// of its mode, so every reader stays valid; its `form.legend` and
// `adv.plot_title_position` read and write its own mode's layout slot.
// A drawing's settings as a new drawing and Reset Settings hold them, fresh
// on each call (services/drawing/settings.js keeps those Reset does not return).
/** @param {DiagramMode} drawingMode */
const settingsDefaults = (drawingMode) => () => ({
  form: createDefaultForm(),
  adv: createDefaultAdv(drawingMode),
  losat: createDefaultLosat(),
  losatProgram: 'blastn', // 'blastn' | 'tblastx' | 'blastp'
  circularConservation: createDefaultCircularConservation(),
  linearComparisonPlan: createDefaultLinearComparisonPlan(),
  linearRecordLayoutEnabled: /** @type {boolean} */ (WEB_UX_PROFILE.linear.arrangeInRowsByDefault),
  linearRecordGap: DEFAULT_LINEAR_RECORD_GAP_PX,
  linearRecordRows: [],
  recordDisplayDrafts: [],
  recordsOff: [],
  unmanagedConfigOverrides: {},
  importedComparisonIntent: createImportedComparisonIntentState(),
  layoutPreferences: createDefaultLayoutPreferences(),
  linearTypographyLinked: true,
  // What the Default colors and Blacklist files set, for their Remove (OV-1280).
  fileImports: /** @type {FileImports} */ ({})
});

/**
 * @param {DrawingId} id
 * @param {DiagramMode} drawingMode
 */
export const createDrawingState = (id, drawingMode) => {
  const defaults = settingsDefaults(drawingMode);
  const initial = defaults();
  const form = reactive(initial.form);
  const adv = reactive(initial.adv);
  const layoutPreferences = reactive(initial.layoutPreferences);
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
  const losat = typedReactive(initial.losat);
  const losatProgram = typedRef(initial.losatProgram);
  const linearComparisonPlan = reactive(initial.linearComparisonPlan);
  const linearRecordLayoutEnabled = typedRef(initial.linearRecordLayoutEnabled);
  const linearRecordRows = reactive(initial.linearRecordRows);
  // The keys of the records this drawing leaves out (record-selection D-02):
  // Linear card uids or Circular source selectors `#N`
  // (services/record-draw-selection.js).
  /** @type {string[]} */
  const recordsOff = reactive(initial.recordsOff);
  // The drawing's inputs (DESIGN §2): its input fields, Linear record rows,
  // input type, and Circular record list and discovery. Its cards, Load, and
  // History write the fields and rows through the slice's writers.
  const inputs = INPUTS_SLICE.create();
  const { linearSeqs } = inputs.view;
  const linearComparisonResolution = typedComputed(() => resolveLinearComparisonPlan({
    plan: linearComparisonPlan,
    sequences: linearSeqs,
    recordsOff,
    layout: linearRecordLayoutEnabled.value ? linearRecordRows : [],
    losatProgram: losatProgram.value,
    blastpMode: losat.blastp?.mode
  }));
  // The settings: the controls write them; Load, History, and Reset install
  // them whole through their slice.
  const settings = SETTINGS_SLICE.create({
    form,
    adv,
    losat,
    losatProgram,
    circularConservation: reactive(initial.circularConservation),
    linearComparisonPlan,
    linearRecordLayoutEnabled,
    linearRecordGap: ref(initial.linearRecordGap),
    linearRecordRows,
    recordDisplayDrafts: reactive(initial.recordDisplayDrafts),
    recordsOff,
    unmanagedConfigOverrides: reactive(initial.unmanagedConfigOverrides),
    importedComparisonIntent: reactive(initial.importedComparisonIntent),
    layoutPreferences,
    linearTypographyLinked: ref(initial.linearTypographyLinked),
    fileImports: reactive(initial.fileImports)
  }, defaults);
  const visibilityRules = VISIBILITY_RULES_SLICE.create();
  const colorRules = COLOR_RULES_SLICE.create();
  const { featureVisibilityManualRules } = visibilityRules.view;
  const palette = PALETTE_SLICE.create();
  const { pendingPaletteName } = palette.view;
  const labelTables = LABEL_TABLES_SLICE.create();
  // The editor edits a slice owns (services/drawing/): the drawing shows its
  // view's members, and only the slice writes them, through its handle.
  const canvas = CANVAS_SLICE.create();
  const annotationSets = ANNOTATION_SETS_SLICE.create();
  const groupNames = GROUP_NAMES_SLICE.create();
  const featureEdits = FEATURE_EDITS_SLICE.create();
  const paint = PAINT_SLICE.create();
  return Object.freeze({
    // Identity, constant: the drawing's key in `state.drawings` and its mode.
    id,
    mode: drawingMode,
    // Settings
    ...settings.view,
    // Inputs
    ...inputs.view,
    // Editor edits, slice by slice (services/drawing/): the visibility and
    // color rules, the palette, the label tables, the annotation sets, the
    // Similarity group names and dormant names keyed by member set (D-21), the
    // canvas padding, the per-feature edits and Feature placements keyed by
    // original-source identity (design Q4, R2), and the paint intent (R2):
    // feature and Legend colors and strokes, and the Legend rows, with the
    // renamed rows the last Generate did not draw (OV-120).
    ...visibilityRules.view,
    ...colorRules.view,
    ...palette.view,
    ...labelTables.view,
    ...annotationSets.view,
    ...groupNames.view,
    ...canvas.view,
    ...featureEdits.view,
    ...paint.view,
    // Derived from the members above
    activeLayoutPreferences,
    // The Linear record cards a Generate draws, in file order (read-only).
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
    ),
    // The handles a slice's writers take.
    slices: Object.freeze({
      visibilityRules: visibilityRules.handle, colorRules: colorRules.handle, palette: palette.handle, labelTables: labelTables.handle, annotationSets: annotationSets.handle, groupNames: groupNames.handle, canvas: canvas.handle,
      featureEdits: featureEdits.handle, paint: paint.handle, settings: settings.handle, inputs: inputs.handle })
  });
};
/**
 * The edits of a drawing that its settings and inputs do not give (R17-6;
 * the items PD-OI-090 〔作成〕 does not copy): its feature edits and Feature
 * placements, its paint (feature fills and strokes, Legend edits), its
 * annotation sets, its label bulk edits, and its group names.
 * @param {DrawingState} drawing
 * @returns {readonly unknown[]}
 */
const editsOf = (drawing) => [
  FEATURE_EDITS_SLICE.capture(drawing.slices.featureEdits),
  PAINT_SLICE.capture(drawing.slices.paint),
  ANNOTATION_SETS_SLICE.capture(drawing.slices.annotationSets),
  labelTextEditsData(drawing.slices.labelTables),
  GROUP_NAMES_SLICE.capture(drawing.slices.groupNames)
];
/** @param {unknown} value The edits as JSON text (a Set or Map as its entries). */
const editsText = (value) => JSON.stringify(value, (_key, item) => (
  item instanceof Set || item instanceof Map ? [...item] : item
));
/**
 * Whether a drawing has edits: those a new drawing of its mode does not have
 * (R17-6: Delete asks first for a drawing with a Result or edits).
 * @param {DrawingState} drawing
 */
export const drawingHasEdits = (drawing) => (
  editsText(editsOf(drawing)) !== editsText(editsOf(createDrawingState(drawing.id, drawing.mode)))
);
const firstDrawing = createDrawingState(FIRST_DRAWING_ID, 'circular');
/** @typedef {typeof firstDrawing} DrawingState */
/**
 * A drawing tab: the drawing's ID and its name, in tab order.
 * @typedef {Readonly<{ id: DrawingId, name: string }>} DrawingTab
 */
/**
 * Every drawing by ID: a read-only map that each store write below replaces, in
 * a shallow ref (a reactive Map would unwrap the members' refs), so a reader
 * of `drawingOf` or `activeDrawing()` follows a write that replaces a drawing
 * under its ID (a Load commits new drawing states; OV-1778). `drawingTabs`
 * holds the order and names, and only the store writers below change either.
 * @type {import('./utils/reactivity.js').Ref<ReadonlyMap<DrawingId, DrawingState>>}
 */
const drawings = typedShallowRef(new Map([[firstDrawing.id, firstDrawing]]));
/** @param {readonly DrawingTab[]} tabs @returns {readonly DrawingTab[]} */
const freezeTabs = (tabs) => Object.freeze(tabs.map(({ id, name }) => Object.freeze({ id, name: String(name) })));
/**
 * The tab order and names: a frozen list, replaced on each write, so a reader
 * of `drawingList()` reacts to it.
 * @type {import('./utils/reactivity.js').Ref<readonly DrawingTab[]>}
 */
const drawingTabs = typedRef(freezeTabs([...drawings.value.values()].map((d) => ({ id: d.id, name: defaultDrawingName(d.mode, d.id) }))));
/** @returns {readonly DrawingState[]} Every drawing, in tab order. */
const drawingList = () => drawingTabs.value.map((tab) => drawingOf(tab.id));
/**
 * The one lookup invariant: a drawing ID names a drawing.
 * @param {DrawingId} id
 * @returns {DrawingState}
 */
const drawingOf = (id) => {
  const drawing = drawings.value.get(id);
  if (!drawing) throw new Error(`Unknown drawing: ${id}`);
  return drawing;
};
/**
 * A refused store write.
 * @param {string} operation @param {string} reason @param {Record<string, unknown>} context
 */
const drawingStoreError = (operation, reason, context) => diagnosticError(
  'INPUT_INVALID', { field: 'drawing', reason, ...context }, { stage: 'drawing-store', operation }
);
/** @param {string} operation @param {DrawingId} id */
const requireDrawingTab = (operation, id) => {
  const index = drawingTabs.value.findIndex((tab) => tab.id === id);
  if (index < 0) throw drawingStoreError(operation, 'UNKNOWN_DRAWING', { drawingId: id });
  return index;
};
/**
 * A drawing's tab (its ID and name): the one reader of a drawing's name for
 * Save and Run info. An unknown ID is a refused read.
 * @param {DrawingId} id
 * @returns {DrawingTab}
 */
const drawingTabOf = (id) => drawingTabs.value[requireDrawingTab('drawingTabOf', id)];
/** @param {number} index @param {number} length */
const clampTabIndex = (index, length) => Math.max(0, Math.min(length, Math.trunc(Number(index) || 0)));
/**
 * Load and a new Session: the drawings and their tabs replace the store. The
 * caller installs the shown drawing next (`installShownDrawing`).
 * @param {readonly DrawingState[]} list @param {readonly DrawingTab[]} tabs
 */
const replaceDrawings = (list, tabs) => {
  /** @type {Map<DrawingId, DrawingState>} */
  const byId = new Map();
  for (const drawing of list) {
    if (byId.has(drawing.id)) throw drawingStoreError('replaceDrawings', 'DUPLICATE_DRAWING', { drawingId: drawing.id });
    byId.set(drawing.id, drawing);
  }
  const tabIds = new Set();
  for (const tab of tabs) {
    if (tabIds.has(tab.id)) throw drawingStoreError('replaceDrawings', 'DUPLICATE_DRAWING', { drawingId: tab.id });
    if (!byId.has(tab.id)) throw drawingStoreError('replaceDrawings', 'UNKNOWN_DRAWING', { drawingId: tab.id });
    tabIds.add(tab.id);
  }
  const untabbed = [...byId.keys()].find((id) => !tabIds.has(id));
  if (untabbed !== undefined || byId.size === 0) {
    throw drawingStoreError('replaceDrawings', 'UNKNOWN_DRAWING', { drawingId: untabbed ?? '' });
  }
  drawings.value = byId;
  drawingTabs.value = freezeTabs(tabs);
};
/** @param {DrawingState} drawing @param {string} name @param {number} index */
const addDrawing = (drawing, name, index) => {
  if (drawings.value.has(drawing.id)) throw drawingStoreError('addDrawing', 'DUPLICATE_DRAWING', { drawingId: drawing.id });
  drawings.value = new Map([...drawings.value, [drawing.id, drawing]]);
  const tabs = [...drawingTabs.value];
  tabs.splice(clampTabIndex(index, tabs.length), 0, { id: drawing.id, name });
  drawingTabs.value = freezeTabs(tabs);
};
/**
 * Removes a drawing that is not shown and not the last one (the caller
 * switches first through `installShownDrawing`), with its artifact slot,
 * which it returns for the step that holds the drawing (R16-3).
 * @param {DrawingId} id
 * @param {Pick<import('./services/artifact-slot.js').ArtifactSlotStore, 'takeSlot'>} slots
 * @returns {Readonly<import('./services/artifact-slot.js').ArtifactSlot> | null}
 */
const removeDrawing = (id, slots) => {
  const index = requireDrawingTab('removeDrawing', id);
  if (drawingTabs.value.length <= 1) throw drawingStoreError('removeDrawing', 'LAST_DRAWING', { drawingId: id });
  if (id === activeDrawingId.value) throw drawingStoreError('removeDrawing', 'SHOWN_DRAWING', { drawingId: id });
  const tabs = [...drawingTabs.value];
  tabs.splice(index, 1);
  drawings.value = new Map([...drawings.value].filter(([key]) => key !== id));
  const slot = slots.takeSlot(id);
  drawingTabs.value = freezeTabs(tabs);
  return slot;
};
/** @param {DrawingId} id @param {number} index */
const moveDrawing = (id, index) => {
  const from = requireDrawingTab('moveDrawing', id);
  const tabs = [...drawingTabs.value];
  const [tab] = tabs.splice(from, 1);
  tabs.splice(clampTabIndex(index, tabs.length), 0, tab);
  drawingTabs.value = freezeTabs(tabs);
};
/** @param {DrawingId} id @param {string} name */
const renameDrawing = (id, name) => {
  const index = requireDrawingTab('renameDrawing', id);
  drawingTabs.value = freezeTabs(drawingTabs.value.map((tab, i) => (i === index ? { id, name } : tab)));
};
/** @returns {DrawingState} The shown drawing. */
const activeDrawing = () => drawingOf(activeDrawingId.value);

export const state = Object.freeze({
  processing,
  processingStatus,
  labelCapacityNotices,
  sessionSavePending,
  sessionImportPending,
  sessionOperationAvailability,
  generationCancelRequested,
  errorLog,
  sessionTitle,
  semanticFileWatchersSuppressed,
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
  inputSummaries,
  pairwiseMatchFactors,
  matchSequenceRegistry,
  svgContent,
  svgResultIdentity,
  zoom,
  layoutRepositionMode,
  isPanning,
  panStart,
  canvasPan,
  previewHiddenUntilShown,
  previewTransitionHeld,
  canvasContainerRef,
  activeDrawingId,
  mode,
  drawingTabs,
  drawingTabOf,
  drawingList,
  drawingOf,
  replaceDrawings,
  addDrawing,
  removeDrawing,
  moveDrawing,
  renameDrawing,
  activeDrawing,
  selectedAnnotation,
  losatCacheInfo,
  losatThreadingStatus,
  inputFiles,
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
  downloadWhiteBackground,
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
  setFeatureEditorStatus,
  featureEditorStatusText,
  featureExtractionPending,
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
  undoneTrackHides,
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
  generatedDrawingId,
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
