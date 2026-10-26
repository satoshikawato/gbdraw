// The drawing context of the Web state (gbdraw/web/js/state.js): every `state`
// key and drawing member in exactly one class. A drawing holds one diagram
// mode's settings, its editor edits, and the values derived from them; only
// the drawing holds them (tests/web/drawing-context.test.mjs). Every `state`
// key is a generated artifact, a project input or cache, an app setting or
// catalog, transient, or the drawing store. Each diagram mode has its own
// drawing (PR-1, Session 46 `modes`).
import { drawingIdOfMode } from '../../../gbdraw/web/js/services/drawing-view.js';
import { CANVAS_SLICE } from '../../../gbdraw/web/js/services/drawing/canvas.js';
import { GROUP_NAMES_SLICE } from '../../../gbdraw/web/js/services/drawing/group-names.js';
import { ANNOTATION_SETS_SLICE } from '../../../gbdraw/web/js/services/drawing/annotation-sets.js';
import { VISIBILITY_RULES_SLICE } from '../../../gbdraw/web/js/services/drawing/visibility-rules.js';
import { COLOR_RULES_SLICE } from '../../../gbdraw/web/js/services/drawing/color-rules.js';
import { PALETTE_SLICE } from '../../../gbdraw/web/js/services/drawing/palette.js';
import { LABEL_TABLES_SLICE } from '../../../gbdraw/web/js/services/drawing/label-tables.js';
import { FEATURE_EDITS_SLICE } from '../../../gbdraw/web/js/services/drawing/feature-edits.js';
import { PAINT_SLICE } from '../../../gbdraw/web/js/services/drawing/paint.js';
import { SETTINGS_SLICE } from '../../../gbdraw/web/js/services/drawing/settings.js';
import { INPUTS_SLICE } from '../../../gbdraw/web/js/services/drawing/inputs.js';

// The slices of a drawing (gbdraw/web/js/services/drawing/), by the name of
// their handle in `drawing.slices`: the editor slices, which alone write their
// members, and the settings, which their controls write (plan-r2 §13 item 3).
export const DRAWING_SLICES = Object.freeze({
  visibilityRules: VISIBILITY_RULES_SLICE, colorRules: COLOR_RULES_SLICE, palette: PALETTE_SLICE, labelTables: LABEL_TABLES_SLICE,
  annotationSets: ANNOTATION_SETS_SLICE,
  groupNames: GROUP_NAMES_SLICE, canvas: CANVAS_SLICE, featureEdits: FEATURE_EDITS_SLICE, paint: PAINT_SLICE,
  settings: SETTINGS_SLICE, inputs: INPUTS_SLICE
});

// The calls of a slice module that only read or build a value: each slice's
// `capture`, under the name its module exports the slice by, and the exports
// listed. Given a handle, each changes neither drawing (the M6 table in
// tests/web/per-mode-drawings.test.mjs, which also fails on an export it does
// not list), so a call of one is no write of its caller (R2,
// tests/web/feature-visibility-rule-writers.test.mjs). Any other call that
// receives a handle writes.
// A registered design-rule allowlist (tools/web-design-rule-guards.json, R2,
// kind `set`): a new entry is an authority-only change.
export const SLICE_READERS = Object.freeze([
  'CANVAS_SLICE.capture', 'GROUP_NAMES_SLICE.capture', 'ANNOTATION_SETS_SLICE.capture', 'VISIBILITY_RULES_SLICE.capture',
  'COLOR_RULES_SLICE.capture', 'PALETTE_SLICE.capture', 'LABEL_TABLES_SLICE.capture', 'FEATURE_EDITS_SLICE.capture',
  'PAINT_SLICE.capture', 'SETTINGS_SLICE.capture', 'INPUTS_SLICE.capture',
  'legendEntryKey', 'fillRecordOf', 'paintOf', 'captureLegendPaint', 'captureFeatureFills', 'featureEditsOf',
  'labelTextEditsData', 'queuedPaletteData', 'createDefaultLabelFilters', 'createDefaultPaletteData'
]);

// The members of a drawing, as its slices name them: the members its controls
// write (the draft keys: the settings, then the inputs) and the editor slices'
// members, which only their slice writes, in the order the drawing spreads them.
const CONTROL_WRITTEN_SLICES = Object.freeze(['settings', 'inputs']);
export const DRAWING_DRAFT_KEYS = Object.freeze(CONTROL_WRITTEN_SLICES.flatMap((name) => DRAWING_SLICES[name].members));
export const DRAWING_EDITOR_KEYS = Object.freeze(Object.entries(DRAWING_SLICES)
  .filter(([name]) => !CONTROL_WRITTEN_SLICES.includes(name)).flatMap(([, slice]) => slice.members));

// Computed from the members above (and project inputs); each drawing derives its own.
export const DRAWING_DERIVED_KEYS = Object.freeze([
  'activeLayoutPreferences', 'drawnLinearSeqs', 'linearComparisonResolution', 'hasLinearComparisonIntent',
  'hasActiveLinearLosatIntent', 'hasActiveLinearUploadIntent', 'featureVisibilityRules', 'hasPendingPaletteDraft'
]);

// The handles of the drawing's slices, which a slice's writers take.
export const DRAWING_HANDLE_KEYS = Object.freeze(['slices']);

// A drawing's identity, constant: its key in `state.drawings` and its mode.
// They hold no settings or edits, so the drawing keys leave them out.
export const DRAWING_IDENTITY_KEYS = Object.freeze(['id', 'mode']);

export const DRAWING_KEYS = Object.freeze([
  ...DRAWING_DRAFT_KEYS, ...DRAWING_EDITOR_KEYS, ...DRAWING_DERIVED_KEYS, ...DRAWING_HANDLE_KEYS
]);

// What a Generate produces and shows: the Results and their catalog, run facts,
// generated Legend inventory and colors, and the user's offsets on the Result.
export const ARTIFACT_KEYS = Object.freeze([
  'results', 'selectedResultIndex', 'lastRunInfo', 'trackSlotResolvedGeometry', 'annotationWarnings',
  'featureIdentityNotices', 'featureEditRemovalCount', 'comparisonWarnings', 'labelFitReports', 'pairwiseMatchFactors',
  'orthogroups', 'collinearGroups', 'featureOrthogroupIndex', 'similarityAlignmentPlan',
  'similarityAlignmentResetReceipt', 'linearRecordTranslations',
  'appliedPaletteName', 'appliedPaletteColors', 'extractedFeatures', 'biologicalFeatures', 'featureCatalog',
  'featureRecordIds', 'editableLabels', 'originalLegendOrder', 'originalLegendColors', 'originalSvgStroke',
  'legendCurrentOffset', 'diagramOffset', 'lengthBarUserOffset', 'plotTitleUserOffset',
  'generatedLegendPosition', 'generatedDrawingId', 'generatedMultiRecordCanvas', 'generatedCircularPlotTitlePosition'
]);

// Inputs and caches that every drawing shares.
export const PROJECT_KEYS = Object.freeze([
  'sessionTitle', 'inputFiles', 'losatCache', 'losatDerivedCache', 'losatCacheInfo', 'proteinIdentityManifest',
  'legacyProteinRawCandidates', 'legacyProteinDerivedEvidence'
]);

// App settings and catalogs outside any drawing, and the shown drawing: its ID
// (`activeDrawingId`, which `installShownDrawing` writes) and its mode (`mode`,
// read from the drawing). LOSAT execution and the rich feature popup are app
// settings (Session `ui.losatExecution`, `ui.richFeaturePopup`).
export const UI_KEYS = Object.freeze([
  'activeDrawingId', 'mode', 'downloadDpi', 'autoLabelReflowEnabled', 'paletteInstantPreviewEnabled', 'paletteDefinitions',
  'paletteNames', 'specificRulePresets', 'featureKeys', 'defaultColorKeys', 'losatExecution', 'richFeaturePopup'
]);

// Never saved: progress and operation flags, input drafts, selection, search,
// dialogs, popups, drag and view state, and views derived from the others.
export const TRANSIENT_KEYS = Object.freeze([
  'loadedFonts', 'processing', 'processingStatus', 'sessionSavePending', 'sessionImportPending', 'sessionOperationAvailability',
  'generationCancelRequested', 'errorLog', 'semanticFileWatchersSuppressed',
  'sessionImportRollbackInProgress', 'failedGeneratePreservedResult', 'generationFailureRecovery',
  'resultPanelTab', 'matchSequenceRegistry', 'svgContent', 'svgResultIdentity', 'zoom', 'layoutRepositionMode',
  'isPanning', 'panStart', 'canvasPan', 'canvasContainerRef',
  'selectedAnnotation', 'losatThreadingStatus', 'selectedOrthogroupAlignmentFeature', 'selectedOrthogroupId',
  'orthogroupSearch', 'orthogroupSortMode', 'showRightDrawer', 'rightDrawerTab', 'linearReorderNotice', 'linearFoldedSources',
  'newSpecRule', 'specificRuleQualifierSuggestions', 'selectedSpecificPreset', 'specificRulePresetLoading',
  'featuresBySvgId', 'selectedFeatureIds', 'selectedFeatureAnchorId', 'featureSelectionStatus',
  'featureSelectionSuppressNextClick', 'featureSelectionDrag', 'selectedFeatureCount', 'selectedFeatures',
  'hasFeatureSelection', 'featureEditorStatus', 'setFeatureEditorStatus', 'featureEditorStatusText', 'featureExtractionPending',
  'selectedFeatureRecordIdx', 'resultGenerationKey', 'featurePanelTab',
  'featureSearchInput', 'featureSearch', 'previewFeatureSearchInput', 'previewFeatureSearchQuery',
  'previewFeatureSearchField', 'previewFeatureSearchQualifierKey', 'previewFeatureSearchUseRegex',
  'previewFeatureSearchMatches', 'previewFeatureSearchMatchDetails', 'previewFeatureSearchActiveIndex',
  'previewFeatureSearchError', 'previewFeatureSearchRenderedCount', 'previewFeatureSearchVisible', 'featureListScrollTop',
  'featureListViewportHeight', 'featureListExpandedId', 'isFeatureDrawerMounted', 'visibleFeatureRows', 'featureRecordPickerVisible',
  'featureListTopSpacerPx', 'featureListBottomSpacerPx', 'labelSearch', 'labelReflowProcessing',
  'labelReflowRequestSeq', 'labelReflowForceRequestSeq', 'labelReflowLastError', 'svgContainer',
  'clickedFeature', 'clickedFeaturePos', 'clickedPairwiseMatch', 'clickedPairwiseMatchPos',
  'pairwiseMatchPopupRef', 'pairwiseMatchPopupDrag', 'pairwiseMatchPopupSize', 'pairwiseMatchPopupResize',
  'featurePopupRef', 'featurePopupDrag', 'featurePopupSize', 'featurePopupResize', 'clickedLabel',
  'clickedLabelPos', 'featureStyleScopeDialog', 'paletteColorsDialog', 'resetColorDialog', 'legendRenameDialog',
  'labelTextScopeDialog', 'featureVisibilityScopeDialog', 'hiddenLabelTextDialog', 'labelOnDialog',
  'sidebarWidth', 'isResizing', 'settingsPaneDisplay', 'settingsPaneDocked', 'settingsPopup', 'legendStrokeOptionsOpen', 'legendDragging', 'legendDragStart',
  'legendOriginalTransform', 'legendInitialTransform', 'diagramDragging', 'diagramDragStart',
  'diagramElementIds', 'diagramElementOriginalTransforms', 'diagramElements', 'lengthBarElement',
  'lengthBarOriginalTransform', 'plotTitleElement', 'plotTitleDragging', 'plotTitleDragStart',
  'plotTitleAutoTransform', 'showCanvasControls', 'skipCaptureBaseConfig',
  'skipExtractOnSvgChange', 'trustedArtifactRestoreInProgress', 'newColorFeat', 'newColorVal',
  'newPriorityRule', 'newFeatureToAdd', 'displayedResultMetadata', 'featureList', 'featureListState', 'filteredFeatures',
  'filteredEditableLabels', 'ruleMatchingPending'
]);

// The drawing store itself: the drawings by ID, their tab order and names
// (`drawingTabs`), and the store writers.
export const STORE_KEYS = Object.freeze([
  'drawingTabs', 'drawingTabOf', 'drawingList', 'drawingOf', 'drawingOfMode', 'activeDrawing',
  'replaceDrawings', 'addDrawing', 'removeDrawing', 'moveDrawing', 'renameDrawing'
]);

// A `state`-shaped node fixture's drawing: its own drawing members, by reference.
export const fixtureDrawing = (fixture) => Object.freeze(withSliceHandles(Object.fromEntries(
  DRAWING_KEYS.filter((key) => Object.hasOwn(fixture, key)).map((key) => [key, fixture[key]])
)));

// A fixture drawing's slice handles (`drawing.slices`), over the members the
// fixture holds: a writer writes the fixture's own objects. A handle is made
// when it is first read, so a fixture without Vue makes none it does not use.
// Not enumerable, so a test that clones or serializes the fixture sees only
// its members; a drawing that has handles keeps them.
const withSliceHandles = (drawing) => {
  if (!drawing || typeof drawing !== 'object' || 'slices' in drawing || Object.isFrozen(drawing)) return drawing;
  // Every drawing holds its inputs (its record rows: one empty row by
  // default); a fixture that names none gets the default, not enumerable.
  const inputs = INPUTS_SLICE.create().view;
  INPUTS_SLICE.members.filter((key) => !(key in drawing)).forEach((key) => {
    Object.defineProperty(drawing, key, { value: inputs[key], configurable: true, writable: true });
  });
  const handles = {};
  const slices = {};
  for (const [name, slice] of Object.entries(DRAWING_SLICES)) {
    Object.defineProperty(slices, name, {
      enumerable: true,
      get: () => {
        if (!handles[name]) {
          const members = Object.fromEntries(slice.members.filter((key) => key in drawing).map((key) => [key, drawing[key]]));
          handles[name] = slice.create(members).handle;
        }
        return handles[name];
      }
    });
  }
  return Object.defineProperty(drawing, 'slices', { value: Object.freeze(slices), configurable: true });
};

// The drawing of a diagram mode in a `state`-shaped value (the real state or a
// fixture from `withDrawings` / `withModeDrawings`): the drawing whose ID the
// mode names.
export const modeDrawing = (state, mode) => state.drawingOf(drawingIdOfMode(mode));

// The drawing store as state.js builds it: the drawings by ID in creation
// order, read through `drawingList`, `drawingOf`, and `drawingOfMode`. A fake
// `state` spreads it.
export const drawingStore = (circular, linear) => {
  const drawings = new Map([[drawingIdOfMode('circular'), circular], [drawingIdOfMode('linear'), linear]]);
  const drawingOf = (id) => {
    if (!drawings.has(id)) throw new Error(`Unknown drawing ${id}`);
    return drawings.get(id);
  };
  return {
    drawingList: () => [...drawings.values()],
    drawingOf,
    drawingOfMode: (mode) => drawingOf(drawingIdOfMode(mode))
  };
};

// The store's members as non-enumerable properties of a fixture.
const storeProperties = (circular, linear) => Object.fromEntries(Object.entries(drawingStore(circular, linear))
  .map(([key, value]) => [key, { value, configurable: true }]));

// The drawing the displayed Result belongs to (state.generatedDrawingId), the
// Circular drawing unless the fixture names one, as in state.js.
const generatedDrawingProperty = (fixture) => (Object.hasOwn(fixture, 'generatedDrawingId') ? {} : {
  generatedDrawingId: { value: { value: drawingIdOfMode('circular') }, configurable: true, writable: true }
});

// The fixture with the drawing store that services and owners read
// (`state.drawingOf`, `state.drawingList()`, `state.activeDrawing()`): the
// fixture is the drawing of both modes, so a member the test replaces later is
// read. A test of per-mode behavior uses `withModeDrawings`. The store is not
// enumerable, so a test that clones or serializes the fixture sees only its
// members.
export const withDrawings = (fixture) => Object.defineProperties(withSliceHandles(fixture), {
  ...storeProperties(fixture, fixture),
  ...(Object.hasOwn(fixture, 'id') ? {} : { id: { value: drawingIdOfMode('circular'), configurable: true } }),
  ...generatedDrawingProperty(fixture),
  activeDrawing: { value: () => fixture, configurable: true }
});

// The fixture with one drawing per mode, as in state.js: `activeDrawing()`
// follows `fixture.mode.value`.
export const withModeDrawings = (fixture, { circular, linear }) => Object.defineProperties(fixture, {
  ...storeProperties(withSliceHandles(circular), withSliceHandles(linear)),
  ...generatedDrawingProperty(fixture),
  activeDrawing: { value: () => (fixture.mode?.value === 'linear' ? linear : circular), configurable: true }
});
