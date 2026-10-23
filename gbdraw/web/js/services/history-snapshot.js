// @ts-check
/** @import { InputFileStore } from './history-files.js' */
/** @import { ArtifactSlot } from './artifact-slot.js' */
/** @import { ResourceBackedComparison } from './canonical-comparisons.js' */
/** @import { AlignmentResetCanonical, AlignmentResetReceipt } from './session-active-config-contract.js' */
/** @import { DrawingId, DrawingIdentity, DrawingView } from './drawing-view.js' */
import { boundInputFileIds } from './drawing/inputs.js';
import { isInputFileRef } from './history-files.js';
import { FILE_IMPORT_CARDS } from './file-imports.js';
/** @import { InputBindings } from './drawing/inputs.js' */
import { validateSimilarityAlignmentResetReceipt } from './session-active-config-contract.js';
import { inputTypeKeyOf } from './drawing-view.js';
import { normalizeFeatureVisibilityRule } from './feature-visibility.js';
import { serializeCleanSvg } from './svg-serialization.js';
import { cloneJsonData } from './json-clone.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { normalizeLayoutPreferences } from './layout-preferences.js';
import {
  isResourceBackedCanonicalComparison,
  mapResourceBackedCanonicalComparison
} from './canonical-comparisons.js';
import { recordStructuralMetric } from './runtime-test-hooks.js';
import {
  ARTIFACT_SLOT_EMPTY, ARTIFACT_SLOT_KEYS, createArtifactSlot, createArtifactSlotStore
} from './artifact-slot.js';
import { toRaw } from '../utils/vue-raw.js';
import { CANVAS_SLICE } from './drawing/canvas.js';
/** @import { CanvasPadding, CanvasSlice } from './drawing/canvas.js' */
import { GROUP_NAMES_SLICE } from './drawing/group-names.js';
import { VISIBILITY_RULES_SLICE } from './drawing/visibility-rules.js';
import { COLOR_RULES_SLICE } from './drawing/color-rules.js';
import { PALETTE_SLICE, queuedPaletteData, restoreArtifactQueuedPalette } from './drawing/palette.js';
import { LABEL_TABLES_SLICE, labelTextEditsData } from './drawing/label-tables.js';
import { FEATURE_EDITS_SLICE } from './drawing/feature-edits.js';
import { PAINT_SLICE, captureFeatureFills, captureLegendPaint } from './drawing/paint.js';
import { SETTINGS_SLICE } from './drawing/settings.js';
/** @import { FeatureFill, KeptLegendEntry, LegendEntry } from './drawing/paint.js' */
/** @import { StrokeOverride } from './style-overrides.js' */
/** @import { StoredLabelOverrideRow } from './label-override-table.js' */
/** @import { FeatureOverrideDraft } from './feature-placement.js' */
/** @import { PaletteColors, PaletteSlice } from './drawing/palette.js' */
/** @import { GroupNamesData } from './drawing/group-names.js' */
/** @import { HistoryArtifactHandle } from './history.js' */

export { cloneJsonData };

// Where an intent holds the composition offsets, keyed by Result identity (B17).
const COMPOSITION_INTENT_PATH = Object.freeze(['ui', 'compositionUserDeltas']);
// Where an intent lists the Legend rows (`editorIntentOf`). A row's `xPos`
// and `yPos` are where the drawn figure lays it out; a redraw moves them,
// which is no edit of the list (OV-1655).
const LEGEND_ROWS_INTENT_PATH = Object.freeze(['editorState', 'legend', 'entries']);
/**
 * The part of an intent value that an edit changes (History's `editedPartOf`):
 * a Legend row list without its rows' layout; any other value whole.
 * @param {readonly string[]} path
 * @param {unknown} value
 * @returns {unknown}
 */
const editedIntentPartOf = (path, value) => (
  Array.isArray(value) && path.length === LEGEND_ROWS_INTENT_PATH.length
    && path.every((key, index) => key === LEGEND_ROWS_INTENT_PATH[index])
    ? value.map((row) => {
        if (!row || typeof row !== 'object') return row;
        const { xPos: _x, yPos: _y, ...edited } = /** @type {Record<string, unknown>} */ (row);
        return edited;
      })
    : value
);

/**
 * The per-feature edits a History intent, checkpoint, or generated artifact holds.
 * @typedef {{
 *   selectedFeatureRecordIdx?: unknown,
 *   featureColorOverrides?: unknown,
 *   featureVisibilityManualRules?: unknown,
 *   featureOverrides?: unknown,
 *   labelOverrideRows?: unknown,
 *   labelTextBulkOverrides?: unknown
 * }} FeatureEditsData
 */
/**
 * The Legend entries and stroke overrides the editor intent holds.
 * @typedef {{
 *   legend?: {
 *     entries?: unknown,
 *     deletedEntries?: unknown,
 *     dormantEntries?: unknown,
 *     colorOverrides?: unknown,
 *     strokeOverrides?: unknown,
 *     addedCaptions?: unknown
 *   } | null,
 *   featureStrokes?: { overrides?: unknown } | null
 * }} EditorIntentData
 */
/** @typedef {{ selectedOrthogroupId?: unknown, selectedOrthogroupAlignmentFeature?: unknown }} OrthogroupIntentData */
/** @typedef {{ orthogroupNameOverrides?: unknown, orthogroupDescriptionOverrides?: unknown, orthogroupDormantOverrides?: unknown }} GroupNameIntentData */
/** @typedef {{ mode?: unknown, defaultSource?: unknown, edges?: unknown }} LinearComparisonPlanData */
/** @typedef {{ x?: unknown, y?: unknown }} OffsetData */
/** @typedef {HistorySnapshotServiceOptions['state']} HistorySnapshotState The app state (state.js has no typedef for it yet). */
/** @typedef {{ name?: unknown, content?: unknown }} ResultData One generated Result, as the fallback capture and restore read it. */
/** @typedef {{ receipt?: AlignmentResetReceipt | null, plan?: unknown, recordTranslations?: unknown }} AlignmentStateData */
/**
 * The `ui` values of an intent, checkpoint, or generated artifact: the fallback fields, the shown
 * mode (a checkpoint's and a handle's only), and the keys the owners of `ui` add.
 * @typedef {FallbackUiData & Record<string, unknown> & { mode?: string, selectedFeatureRecordIdx?: unknown }} IntentUiData
 */
/**
 * A History intent: one drawing's settings and edits, the file bindings of its
 * slots, and the shared selections its steps record (PD-OI-090); no
 * Session-level value (R15-23).
 * @typedef {object} HistoryIntentData
 * @property {HistorySnapshotData} [config] The Settings capture.
 * @property {IntentUiData} [ui]
 * @property {FeatureEditsData} [features]
 * @property {EditorIntentData} [editorState]
 * @property {GroupNameIntentData} [groupNames]
 * @property {FilesData} [files]
 * @property {AlignmentStateData} [alignmentState] A Linear drawing's.
 * @property {OrthogroupIntentData} [orthogroupState]
 */
/**
 * The steps of an intent that `applyHistoryIntent` is told changed.
 * @typedef {{ changes?: { path?: unknown }[], [key: string]: unknown }} HistoryIntentContext
 */
/**
 * An artifact checkpoint: the shown drawing's intent and the generated Result domains.
 * @typedef {object} ArtifactCheckpointData
 * @property {HistorySnapshotData} [config]
 * @property {FilesData} [files]
 * @property {AlignmentStateData} [alignmentState]
 * @property {IntentUiData} [ui]
 * @property {HistorySnapshotData[]} [results]
 * @property {FeatureEditsData & { extractedFeatures?: unknown, featureRecordIds?: unknown }} [features]
 * @property {EditorIntentData} [editorState]
 * @property {OrthogroupIntentData & GroupNameIntentData & { groups?: unknown }} [orthogroupState]
 * @property {{ lastRunInfo?: unknown, pairwiseMatchFactors?: unknown }} [runState]
 */
/**
 * The presentation values a generated artifact handle keeps with its intent.
 * @typedef {{
 *   resultGenerationKey?: unknown, resultPanelTab?: unknown, errorLog?: unknown, labelTextScopeDialog?: unknown,
 *   featureEditorStatus?: unknown, featureExtractionPending?: unknown, featureExtractionError?: unknown
 * }} PresentationData
 */
/**
 * The mutable intent of a generated artifact handle.
 * @typedef {object} GeneratedArtifactIntent
 * @property {IntentUiData} ui
 * @property {FeatureEditsData} features
 * @property {EditorIntentData & { originalSvgStroke?: unknown }} editorState
 * @property {OrthogroupIntentData & GroupNameIntentData} orthogroupState
 * @property {unknown} alignmentState
 * @property {PresentationData} presentation
 */
/**
 * One Linear record row of an intent's files: the File descriptors and the metadata the row keeps.
 * @typedef {{
 *   uid?: unknown, gb?: unknown, gff?: unknown, fasta?: unknown, depth?: unknown,
 *   losat_gencode?: unknown, definition?: unknown, record_subtitle?: unknown, file_definition?: unknown,
 *   file_subtitle?: unknown, inferred_definition?: unknown, region_record_id?: unknown,
 *   region_start?: unknown, region_end?: unknown, region_reverse?: unknown
 * }} FilesDataSequence
 */
/**
 * The files of an intent or checkpoint: File descriptors (`fileStore.describeValue`) by input.
 * @typedef {{
 *   c_gb?: unknown, c_gff?: unknown, c_fasta?: unknown, c_depth?: unknown,
 *   c_conservation_blasts?: unknown, c_conservation_blasts_source?: unknown,
 *   c_conservation_fastas?: unknown, c_conservation_sequence_sources?: unknown,
 *   d_color?: unknown, t_color?: unknown, blacklist?: unknown, whitelist?: unknown, qualifier_priority?: unknown,
 *   linearSeqs?: FilesDataSequence[],
 *   linearComparisons?: { id?: unknown, file?: unknown }[],
 *   linearCanonicalComparisons?: ResourceBackedComparison[]
 * }} FilesData
 */
/**
 * The `ui` values that belong to a drawing (Session 46 `modes.<m>.ui`, and its LOSAT program).
 * @typedef {{
 *   canvasPadding?: Record<string, unknown> | null,
 *   layoutPreferences?: unknown,
 *   linearTypographyLinked?: unknown,
 *   pendingPaletteName?: unknown,
 *   pendingPaletteColors?: unknown,
 *   losatProgram?: unknown
 * }} DrawingUiData
 */
/**
 * The UI state the fallback restore reads when the service has no UI owner.
 * @typedef {DrawingUiData & {
 *   title?: unknown,
 *   cInputType?: unknown,
 *   lInputType?: unknown,
 *   downloadDpi?: unknown,
 *   generatedLegendPosition?: unknown,
 *   generatedMode?: unknown,
 *   generatedMultiRecordCanvas?: unknown,
 *   generatedCircularPlotTitlePosition?: unknown,
 *   autoLabelReflow?: unknown,
 *   paletteInstantPreviewEnabled?: unknown,
 *   appliedPaletteName?: unknown,
 *   appliedPaletteColors?: unknown,
 *   legendCurrentOffset?: OffsetData | null,
 *   diagramOffset?: OffsetData | null,
 *   lengthBarUserOffset?: OffsetData | null,
 *   plotTitleUserOffset?: OffsetData | null,
 *   selectedResultIndex?: number
 * }} FallbackUiData
 * @typedef {{
 *   title?: unknown, downloadDpi?: unknown, autoLabelReflow?: unknown, paletteInstantPreviewEnabled?: unknown,
 *   losatExecution?: unknown, featurePanelTab?: unknown
 * }} SessionUiData
 *   The Session-level `ui` values (`SESSION_UI_KEYS`), which the History router's project stream records.
 */
/**
 * The draft inputs (annotation, preset, and rule forms) a History intent holds.
 * @typedef {{
 *   selectedAnnotation?: unknown,
 *   selectedSpecificPreset?: unknown,
 *   newSpecRule?: unknown,
 *   newPriorityRule?: unknown,
 *   newColorFeat?: unknown,
 *   newColorVal?: unknown,
 *   newFeatureToAdd?: unknown
 * }} DraftIntentData
 */

/**
 * @param {unknown} value
 * @returns {Record<string, unknown>}
 */
const clonePlainObject = (value) => {
  if (!value || typeof value !== 'object' || Array.isArray(value)) return {};
  // A non-null, non-array object that survives the JSON clone is a plain JSON object.
  return /** @type {Record<string, unknown>} */ (cloneJsonData(value));
};

/**
 * @param {Record<string, unknown> | null | undefined} target
 * @param {Record<string, unknown> | null | undefined} source
 */
const replacePlainObject = (target, source) => {
  if (!target || typeof target !== 'object') return;
  Object.keys(target).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = value;
  });
};

/** @param {FeatureEditsData} [features] */
const buildFeatureIntentData = (features = {}) => ({
  featureColorOverrides: clonePlainObject(features.featureColorOverrides),
  featureVisibilityManualRules: cloneFeatureVisibilityRules(features.featureVisibilityManualRules),
  featureOverrides: clonePlainObject(features.featureOverrides),
  labelOverrideRows: cloneLabelOverrideRows(features.labelOverrideRows),
  labelTextBulkOverrides: clonePlainObject(features.labelTextBulkOverrides)
});

// A drawing's Legend paint as a History step holds it (`editorState`).
/** @param {ReturnType<typeof captureLegendPaint>} paint A capture (a copy). */
const editorIntentOf = (paint) => ({
  legend: {
    entries: paint.legendEntries,
    deletedEntries: paint.deletedLegendEntries,
    dormantEntries: paint.dormantLegendEntries,
    colorOverrides: paint.legendColorOverrides,
    strokeOverrides: paint.legendStrokeOverrides,
    addedCaptions: paint.addedLegendCaptions
  },
  featureStrokes: {
    overrides: paint.featureStrokeOverrides
  }
});

/** @param {OrthogroupIntentData} [orthogroupState] */
const buildOrthogroupIntentData = (orthogroupState = {}) => ({
  selectedOrthogroupId: String(orthogroupState.selectedOrthogroupId || ''),
  selectedOrthogroupAlignmentFeature: String(orthogroupState.selectedOrthogroupAlignmentFeature || '')
});

// A drawing's group names (its Session 46 `config.webEdits`).
/** @param {DrawingView} drawing */
const buildGroupNameIntentData = (drawing) => GROUP_NAMES_SLICE.capture(drawing.slices.groupNames);

// The `ui` values of the Session and the app, which no drawing owns
// (tests/web/helpers/drawing-state.mjs PROJECT_KEYS, UI_KEYS, TRANSIENT_KEYS):
// the router's project stream records them with the drafts (R15-23), and a
// drawing's intent and checkpoint hold none of them.
const SESSION_UI_KEYS = Object.freeze([
  'title', 'downloadDpi', 'autoLabelReflow', 'paletteInstantPreviewEnabled', 'losatExecution', 'featurePanelTab'
]);
/** @param {Record<string, unknown>} ui */
const withoutSessionUi = (ui) => {
  const drawingUi = { ...ui };
  SESSION_UI_KEYS.forEach((key) => { delete drawingUi[key]; });
  return drawingUi;
};

// Neither the selected Result nor the offsets read from its mounted root is
// an edit, so neither is intent; composition offsets are recorded per Result
// in `compositionUserDeltas` (B17).
// A drawing switch is no step (PD-OI-090), so the intent holds no `mode`; a
// drawing's intent holds the input type of its own slots.
/**
 * @param {Record<string, unknown>} ui
 * @param {DrawingIdentity} drawing
 */
const buildUiIntentData = (ui, drawing) => {
  const intent = withoutSessionUi(clonePlainObject(ui));
  delete intent.mode;
  delete intent[drawing.mode === 'linear' ? 'cInputType' : 'lInputType'];
  delete intent.generatedLegendPosition;
  delete intent.generatedMode;
  delete intent.generatedMultiRecordCanvas;
  delete intent.generatedCircularPlotTitlePosition;
  delete intent.selectedResultIndex;
  delete intent.legendCurrentOffset;
  delete intent.diagramOffset;
  delete intent.plotTitleUserOffset;
  return intent;
};

// A generated artifact's features hold the Features-list record it showed;
// an intent's per-mode features do not (the intent's `ui` holds it).
/**
 * @param {Record<string, any>} state
 * @param {DrawingView} drawing
 * @param {FeatureEditsData} [features]
 */
const applyFeatureIntentData = (state, drawing, features = {}) => {
  if (Object.hasOwn(features, 'selectedFeatureRecordIdx')) {
    setRef(
      state.selectedFeatureRecordIdx,
      Number.isInteger(features.selectedFeatureRecordIdx) ? features.selectedFeatureRecordIdx : 0
    );
  }
  PAINT_SLICE.install(drawing.slices.paint, { featureColorOverrides: capturedFills(features) });
  replaceFeatureEditState(drawing, features);
};

// A captured step's feature color edits (History installs its copy as is, R11;
// the cast: a captured copy of the drawing's fills is fills).
/** @param {FeatureEditsData} features */
const capturedFills = (features) => /** @type {Record<string, FeatureFill>} */ (clonePlainObject(features.featureColorOverrides));

/** @param {unknown} rows @returns {KeptLegendEntry[]} A captured list of Legend rows (the cast), else none. */
const capturedRows = (rows) => (Array.isArray(rows) ? /** @type {KeptLegendEntry[]} */ (rows) : []);

/**
 * @param {DrawingView} drawing
 * @param {EditorIntentData} [editorState]
 */
const applyEditorIntentData = (drawing, editorState = {}) => {
  const legend = editorState.legend || {};
  // A captured copy of the drawing's paint is its paint (the casts).
  PAINT_SLICE.install(drawing.slices.paint, {
    legendEntries: /** @type {LegendEntry[]} */ (capturedRows(legend.entries)),
    deletedLegendEntries: capturedRows(legend.deletedEntries),
    dormantLegendEntries: capturedRows(legend.dormantEntries),
    legendColorOverrides: /** @type {Record<string, string>} */ (clonePlainObject(legend.colorOverrides)),
    legendStrokeOverrides: /** @type {Record<string, StrokeOverride>} */ (clonePlainObject(legend.strokeOverrides)),
    addedLegendCaptions: Array.isArray(legend.addedCaptions) ? legend.addedCaptions.map((entry) => String(entry || '')) : [],
    featureStrokeOverrides: /** @type {Record<string, StrokeOverride>} */ (clonePlainObject(editorState?.featureStrokes?.overrides))
  });
};

/**
 * @param {Record<string, any>} state
 * @param {OrthogroupIntentData} [orthogroupState]
 */
const applyOrthogroupIntentData = (state, orthogroupState = {}) => {
  setRef(state.selectedOrthogroupId, String(orthogroupState.selectedOrthogroupId || ''));
  setRef(
    state.selectedOrthogroupAlignmentFeature,
    String(orthogroupState.selectedOrthogroupAlignmentFeature || '')
  );
};

/**
 * @param {DrawingView} drawing
 * @param {GroupNameIntentData} [names]
 */
const applyGroupNameIntentData = (drawing, names = {}) => {
  // An intent, checkpoint, or generated artifact holds what the slice captured; a missing map is empty.
  GROUP_NAMES_SLICE.install(drawing.slices.groupNames, /** @type {GroupNamesData} */ ({
    orthogroupNameOverrides: clonePlainObject(names?.orthogroupNameOverrides),
    orthogroupDescriptionOverrides: clonePlainObject(names?.orthogroupDescriptionOverrides),
    orthogroupDormantOverrides: clonePlainObject(names?.orthogroupDormantOverrides)
  }));
};

/** @param {LinearComparisonPlanData | null} [plan] */
const cloneLinearComparisonPlanMetadata = (plan = {}) => ({
  mode: String(plan?.mode || 'none'),
  defaultSource: String(plan?.defaultSource || 'losat'),
  edges: (Array.isArray(plan?.edges) ? plan.edges : []).map((edge) => ({
    id: String(edge?.id || ''),
    queryUid: String(edge?.queryUid || ''),
    subjectUid: String(edge?.subjectUid || ''),
    included: edge?.included === true,
    fileActive: edge?.fileActive === true,
    losatFilenameActive: edge?.losatFilenameActive === true,
    source: String(edge?.source || 'upload'),
    losatFilename: String(edge?.losatFilename || '')
  }))
});

// A History step restores a comparison plan without its files (a later
// files restore or the retained files give them back).
/**
 * @param {DrawingView} drawing
 * @param {LinearComparisonPlanData} [plan]
 */
const installComparisonPlan = (drawing, plan = {}) => {
  const cloned = cloneLinearComparisonPlanMetadata(plan);
  SETTINGS_SLICE.install(drawing.slices.settings, {
    linearComparisonPlan: { ...cloned, edges: cloned.edges.map((edge) => ({ ...edge, file: null })) }
  });
};

/** @param {unknown} rows The saved label table (`canonicalLabelOverrideRows`). */
const cloneLabelOverrideRows = (rows) => (Array.isArray(rows) ? cloneJsonData(rows) : []);

/** @param {unknown} rules */
const cloneFeatureVisibilityRules = (rules) => (
  Array.isArray(rules) ? rules.map((rule) => normalizeFeatureVisibilityRule(rule)) : []
);

// Manual visibility rules, identity-keyed per-feature edits, and bulk label
// edits with the saved label table they replace (OV-1102): the feature edit
// intent History restores as captured (R11).
/**
 * @param {DrawingView} drawing
 * @param {FeatureEditsData} [features]
 */
const replaceFeatureEditState = (drawing, features = {}) => {
  VISIBILITY_RULES_SLICE.install(drawing.slices.visibilityRules, {
    featureVisibilityManualRules: cloneFeatureVisibilityRules(features.featureVisibilityManualRules)
  });
  // A History step holds the draft the slice captured (the cast).
  FEATURE_EDITS_SLICE.install(drawing.slices.featureEdits, {
    featureOverrides: /** @type {FeatureOverrideDraft} */ (clonePlainObject(features.featureOverrides))
  });
  LABEL_TABLES_SLICE.install(drawing.slices.labelTables, {
    // A History step holds the rows the slice captured (the cast).
    canonicalLabelOverrideRows: /** @type {StoredLabelOverrideRow[]} */ (cloneLabelOverrideRows(features.labelOverrideRows)),
    labelTextBulkOverrides: /** @type {Record<string, string>} */ (clonePlainObject(features.labelTextBulkOverrides))
  });
};

/**
 * @param {unknown} target A Vue ref, or any other value (a non-ref is left as it is).
 * @param {unknown} value
 */
const setRef = (target, value) => {
  if (target && typeof target === 'object' && 'value' in target) {
    target.value = value;
  }
};

/**
 * @template T, F
 * @param {{ value: T } | null | undefined} target A Vue ref (the fallback is returned for a non-ref).
 * @param {F} fallback The value for a target that is not a ref.
 * @returns {T | F}
 */
const getRef = (target, fallback) => (
  target && typeof target === 'object' && 'value' in target
    ? target.value
    : fallback
);

/**
 * @param {unknown} target A Vue ref, or any other value (a non-ref is left as it is).
 * @param {unknown} value
 */
const setGeneratedArtifactRef = (target, value) => {
  if (target && typeof target === 'object' && 'value' in target) target.value = value;
};

/**
 * @template T, F
 * @param {{ value: T } | null | undefined} target A Vue ref (the fallback is returned for a non-ref).
 * @param {F} fallback The value for a target that is not a ref.
 * @returns {T | F}
 */
const getGeneratedArtifactRef = (target, fallback) => (
  target && typeof target === 'object' && 'value' in target ? target.value : fallback
);

/** @param {Record<string, any>} state */
const buildDraftIntentData = (state) => ({
  selectedAnnotation: cloneJsonData(getRef(state.selectedAnnotation, null)),
  selectedSpecificPreset: String(getRef(state.selectedSpecificPreset, '') || ''),
  newSpecRule: clonePlainObject(state.newSpecRule),
  newPriorityRule: clonePlainObject(state.newPriorityRule),
  newColorFeat: String(getRef(state.newColorFeat, '') || ''),
  newColorVal: String(getRef(state.newColorVal, '') || ''),
  newFeatureToAdd: String(getRef(state.newFeatureToAdd, '') || '')
});

/**
 * @param {Record<string, any>} state
 * @param {DraftIntentData} [drafts]
 */
const applyDraftIntentData = (state, drafts = {}) => {
  setRef(state.selectedAnnotation, cloneJsonData(drafts.selectedAnnotation) || null);
  setRef(state.selectedSpecificPreset, String(drafts.selectedSpecificPreset || ''));
  replacePlainObject(state.newSpecRule, clonePlainObject(drafts.newSpecRule));
  replacePlainObject(state.newPriorityRule, clonePlainObject(drafts.newPriorityRule));
  setRef(state.newColorFeat, String(drafts.newColorFeat || ''));
  setRef(state.newColorVal, String(drafts.newColorVal || ''));
  setRef(state.newFeatureToAdd, String(drafts.newFeatureToAdd || ''));
};

// A History step holds the padding the canvas slice captured (the cast).
/**
 * @param {{ readonly slices: { readonly canvas: CanvasSlice } }} drawing
 * @param {DrawingUiData} ui
 */
const installCanvasPadding = (drawing, ui) => {
  if (!ui.canvasPadding || typeof ui.canvasPadding !== 'object') return;
  CANVAS_SLICE.install(drawing.slices.canvas, { canvasPadding: /** @type {CanvasPadding} */ (clonePlainObject(ui.canvasPadding)) });
};
// A History step holds the queued palette the palette slice captured (the cast).
/**
 * @param {{ readonly slices: { readonly palette: PaletteSlice } }} drawing
 * @param {{ pendingPaletteName?: unknown, pendingPaletteColors?: unknown }} ui
 */
const installQueuedPalette = (drawing, { pendingPaletteName, pendingPaletteColors }) => PALETTE_SLICE.install(drawing.slices.palette, {
  pendingPaletteName: pendingPaletteName === undefined ? undefined : String(pendingPaletteName || ''),
  pendingPaletteColors: pendingPaletteColors && typeof pendingPaletteColors === 'object'
    ? /** @type {PaletteColors} */ (clonePlainObject(pendingPaletteColors)) : undefined
});
const nextFrame = () => /** @type {Promise<void>} */ (new Promise((resolve) => {
  if (typeof window !== 'undefined' && typeof window.requestAnimationFrame === 'function') {
    window.requestAnimationFrame(() => resolve());
  } else {
    setTimeout(resolve, 0);
  }
}));

/** @param {HistorySnapshotState} state */
const closeTransientState = (state) => {
  setRef(state.clickedFeature, null);
  setRef(state.clickedPairwiseMatch, null);
  setRef(state.clickedLabel, null);
  if (state.featureStyleScopeDialog) state.featureStyleScopeDialog.show = false;
  if (state.resetColorDialog) state.resetColorDialog.show = false;
  if (state.legendRenameDialog) state.legendRenameDialog.show = false;
  if (state.labelTextScopeDialog) state.labelTextScopeDialog.show = false;
  if (state.featureVisibilityScopeDialog) state.featureVisibilityScopeDialog.show = false;
  if (state.hiddenLabelTextDialog) state.hiddenLabelTextDialog.show = false;
  if (state.featurePopupDrag) state.featurePopupDrag.active = false;
  if (state.featurePopupResize) state.featurePopupResize.active = false;
  if (state.pairwiseMatchPopupDrag) state.pairwiseMatchPopupDrag.active = false;
  if (state.pairwiseMatchPopupResize) state.pairwiseMatchPopupResize.active = false;
};

/**
 * @param {Record<string, any>} state
 * @param {DrawingView} drawing
 */
const buildFallbackUiStateData = (state, drawing) => ({
  title: getRef(state.sessionTitle, ''),
  mode: getRef(state.mode, 'circular'),
  cInputType: getRef(state.cInputType, 'gb'),
  lInputType: getRef(state.lInputType, 'gb'),
  ...SETTINGS_SLICE.capture(drawing.slices.settings, ['losatProgram', 'layoutPreferences']),
  selectedResultIndex: getRef(state.selectedResultIndex, 0),
  downloadDpi: getRef(state.downloadDpi, 300),
  canvasPadding: CANVAS_SLICE.capture(drawing.slices.canvas).canvasPadding,
  generatedLegendPosition: getRef(state.generatedLegendPosition, 'left'),
  generatedMode: getRef(state.generatedMode, 'circular'),
  generatedMultiRecordCanvas: Boolean(getRef(state.generatedMultiRecordCanvas, false)),
  generatedCircularPlotTitlePosition: getRef(state.generatedCircularPlotTitlePosition, 'none'),
  autoLabelReflow: Boolean(getRef(state.autoLabelReflowEnabled, false)),
  paletteInstantPreviewEnabled: Boolean(getRef(state.paletteInstantPreviewEnabled, false)),
  appliedPaletteName: getRef(state.appliedPaletteName, 'default'),
  appliedPaletteColors: clonePlainObject(getRef(state.appliedPaletteColors, {})),
  ...queuedPaletteData(drawing.slices.palette),
  legendCurrentOffset: { ...(state.legendCurrentOffset || {}) },
  diagramOffset: { ...(state.diagramOffset || {}) },
  lengthBarUserOffset: { ...(state.lengthBarUserOffset || {}) },
  plotTitleUserOffset: { ...(state.plotTitleUserOffset || {}) }
});

// The fallback's Session-level values (R15-23).
/**
 * @param {HistorySnapshotState} state
 * @param {SessionUiData} ui
 */
const applyFallbackSessionUiData = (state, ui) => {
  if (typeof ui.title === 'string') setRef(state.sessionTitle, ui.title);
  if (ui.downloadDpi) setRef(state.downloadDpi, ui.downloadDpi);
  if (Object.hasOwn(ui, 'autoLabelReflow')) setRef(state.autoLabelReflowEnabled, Boolean(ui.autoLabelReflow));
  if (Object.hasOwn(ui, 'paletteInstantPreviewEnabled')) {
    setRef(state.paletteInstantPreviewEnabled, Boolean(ui.paletteInstantPreviewEnabled));
  }
};

// `mode` is not restored: a drawing switch is no History step (PD-OI-090).
/**
 * @param {HistorySnapshotState} state
 * @param {DrawingView} drawing
 * @param {FallbackUiData} [ui]
 */
const applyFallbackUiStateData = (state, drawing, ui = {}) => {
  applyFallbackSessionUiData(state, ui);
  const inputType = inputTypeKeyOf(drawing);
  if (ui[inputType]) setRef(state[inputType], ui[inputType]);
  SETTINGS_SLICE.install(drawing.slices.settings, {
    losatProgram: ui.losatProgram || undefined, layoutPreferences: normalizeLayoutPreferences(ui.layoutPreferences)
  });
  if (ui.generatedLegendPosition) setRef(state.generatedLegendPosition, ui.generatedLegendPosition);
  if (ui.generatedMode) setRef(state.generatedMode, ui.generatedMode);
  if (Object.prototype.hasOwnProperty.call(ui, 'generatedMultiRecordCanvas')) {
    setRef(state.generatedMultiRecordCanvas, Boolean(ui.generatedMultiRecordCanvas));
  }
  if (ui.generatedCircularPlotTitlePosition) {
    setRef(state.generatedCircularPlotTitlePosition, ui.generatedCircularPlotTitlePosition);
  }
  installCanvasPadding(drawing, ui);
  if (ui.appliedPaletteName !== undefined) setRef(state.appliedPaletteName, String(ui.appliedPaletteName || 'default'));
  if (ui.appliedPaletteColors) setRef(state.appliedPaletteColors, clonePlainObject(ui.appliedPaletteColors));
  installQueuedPalette(drawing, ui);
  if (state.legendCurrentOffset && ui.legendCurrentOffset) {
    state.legendCurrentOffset.x = Number(ui.legendCurrentOffset.x) || 0;
    state.legendCurrentOffset.y = Number(ui.legendCurrentOffset.y) || 0;
  }
  if (state.diagramOffset && ui.diagramOffset) {
    state.diagramOffset.x = Number(ui.diagramOffset.x) || 0;
    state.diagramOffset.y = Number(ui.diagramOffset.y) || 0;
  }
  if (state.lengthBarUserOffset && ui.lengthBarUserOffset) {
    state.lengthBarUserOffset.x = Number(ui.lengthBarUserOffset.x) || 0;
    state.lengthBarUserOffset.y = Number(ui.lengthBarUserOffset.y) || 0;
  }
  if (state.plotTitleUserOffset && ui.plotTitleUserOffset) {
    state.plotTitleUserOffset.x = Number(ui.plotTitleUserOffset.x) || 0;
    state.plotTitleUserOffset.y = Number(ui.plotTitleUserOffset.y) || 0;
  }
  if (Number.isInteger(ui.selectedResultIndex)) {
    const count = Array.isArray(getRef(state.results, [])) ? getRef(state.results, []).length : 0;
    // The guard above passes only integers; `Number.isInteger` does not narrow the type.
    setRef(state.selectedResultIndex, count > 0 ? Math.max(0, Math.min(/** @type {number} */ (ui.selectedResultIndex), count - 1)) : 0);
  }
};

/**
 * @param {Record<string, any>} state
 * @param {DrawingView} drawing
 */
const buildFallbackFeatureStateData = (state, drawing) => ({
  extractedFeatures: cloneJsonData(getRef(state.extractedFeatures, [])) || [],
  featureRecordIds: cloneJsonData(getRef(state.featureRecordIds, [])) || [],
  selectedFeatureRecordIdx: getRef(state.selectedFeatureRecordIdx, 0),
  featureColorOverrides: captureFeatureFills(drawing.slices.paint),
  featureVisibilityManualRules: cloneFeatureVisibilityRules(VISIBILITY_RULES_SLICE.capture(drawing.slices.visibilityRules).featureVisibilityManualRules),
  featureOverrides: FEATURE_EDITS_SLICE.capture(drawing.slices.featureEdits).featureOverrides,
  ...labelTextEditsData(drawing.slices.labelTables)
});

/**
 * @param {Record<string, any>} state
 * @param {DrawingView} drawing
 * @param {FeatureEditsData & { extractedFeatures?: unknown, featureRecordIds?: unknown }} [features]
 */
const applyFallbackFeatureStateData = (state, drawing, features = {}) => {
  setRef(state.extractedFeatures, cloneJsonData(features.extractedFeatures) || []);
  setRef(state.featureRecordIds, cloneJsonData(features.featureRecordIds) || []);
  setRef(
    state.selectedFeatureRecordIdx,
    Number.isInteger(features.selectedFeatureRecordIdx) ? features.selectedFeatureRecordIdx : 0
  );
  PAINT_SLICE.install(drawing.slices.paint, { featureColorOverrides: capturedFills(features) });
  replaceFeatureEditState(drawing, features);
};

/**
 * @param {Record<string, any>} state
 * @param {DrawingView} drawing
 */
const buildFallbackOrthogroupStateData = (state, drawing) => ({
  groups: cloneJsonData(getRef(state.orthogroups, [])) || [],
  selectedOrthogroupId: getRef(state.selectedOrthogroupId, ''),
  selectedOrthogroupAlignmentFeature: getRef(state.selectedOrthogroupAlignmentFeature, ''),
  ...buildGroupNameIntentData(drawing)
});

/**
 * @param {Record<string, any>} state
 * @param {DrawingView} drawing
 * @param {OrthogroupIntentData & GroupNameIntentData & { groups?: unknown }} [data]
 */
const applyFallbackOrthogroupStateData = (state, drawing, data = {}) => {
  setRef(state.orthogroups, cloneJsonData(data.groups) || []);
  setRef(state.selectedOrthogroupId, String(data.selectedOrthogroupId || ''));
  setRef(state.selectedOrthogroupAlignmentFeature, String(data.selectedOrthogroupAlignmentFeature || ''));
  applyGroupNameIntentData(drawing, data);
};

/** @param {HistorySnapshotState} state */
const buildFallbackResultsData = (state) => {
  const currentSvg = (() => {
    const svg = state.svgContainer?.value?.querySelector?.('svg');
    if (!svg || typeof XMLSerializer === 'undefined') return null;
    return serializeCleanSvg(svg);
  })();
  const selected = getRef(state.selectedResultIndex, 0);
  return (getRef(state.results, []) || []).map((/** @type {ResultData | null | undefined} */ result, /** @type {number} */ index) => ({
    name: result?.name || `Result ${index + 1}`,
    content: index === selected && currentSvg ? currentSvg : String(result?.content || '')
  }));
};

/**
 * @param {HistorySnapshotState} state
 * @param {unknown} [results]
 */
const applyFallbackResultsData = (state, results = []) => {
  setRef(
    state.results,
    Array.isArray(results)
      ? results.map((result, index) => ({
          name: result?.name || `Result ${index + 1}`,
          content: String(result?.content || '')
        }))
      : []
  );
};

// A drawing's intent holds the file bindings of its slots by reference
// (PD-OI-090): Circular's `c_*` inputs, including the BLAST rows and comparison
// sequences of a LOSAT-cache replay (its ring rows and managed track slots name
// those rows, so Undo restores them together, B23), or Linear's records and
// comparison files. The five table files stay shared until each drawing binds
// its own; every drawing's intent holds them.
/**
 * @param {Record<string, any>} state
 * @param {DrawingView} drawing
 * @param {HistorySnapshotFileStore} fileStore
 * @returns {FilesData}
 */
const buildIntentFilesData = (state, drawing, fileStore) => {
  const shared = Object.fromEntries(FILE_IMPORT_CARDS.map((key) => [key, fileStore.describeValue(state.files?.[key])]));
  if (drawing.mode !== 'linear') {
    return {
      c_gb: fileStore.describeValue(state.files?.c_gb),
      c_gff: fileStore.describeValue(state.files?.c_gff),
      c_fasta: fileStore.describeValue(state.files?.c_fasta),
      c_depth: fileStore.describeValue(state.files?.c_depth),
      c_conservation_blasts: fileStore.describeValue(state.files?.c_conservation_blasts || []),
      c_conservation_blasts_source: state.files?.c_conservation_blasts_source === 'losat-cache'
        ? 'losat-cache'
        : null,
      c_conservation_fastas: fileStore.describeValue(state.files?.c_conservation_fastas || []),
      c_conservation_sequence_sources: fileStore.describeValue(state.files?.c_conservation_sequence_sources || []),
      ...shared
    };
  }
  return {
    ...shared,
    linearSeqs: Array.from(state.linearSeqs || []).map((seq) => ({
      uid: seq.uid,
      gb: fileStore.describeValue(seq.gb),
      gff: fileStore.describeValue(seq.gff),
      fasta: fileStore.describeValue(seq.fasta),
      depth: fileStore.describeValue(seq.depth),
      losat_gencode: seq.losat_gencode ?? 1,
      definition: seq.definition ?? '',
      record_subtitle: seq.record_subtitle ?? '',
      file_definition: seq.file_definition ?? '',
      file_subtitle: seq.file_subtitle ?? '',
      inferred_definition: seq.inferred_definition ?? '',
      region_record_id: seq.region_record_id ?? '',
      region_start: seq.region_start ?? null,
      region_end: seq.region_end ?? null,
      region_reverse: Boolean(seq.region_reverse)
    })),
    linearComparisons: Array.from(drawing.linearComparisonPlan?.edges || []).map((edge) => ({
      id: String(edge?.id || ''),
      file: fileStore.describeValue(edge?.file)
    }))
  };
};

// A Linear drawing's artifact checkpoint adds its generated comparisons.
/**
 * @param {Record<string, any>} state
 * @param {DrawingView} drawing
 * @param {HistorySnapshotFileStore} fileStore
 * @returns {FilesData}
 */
const buildFilesData = (state, drawing, fileStore) => {
  const files = buildIntentFilesData(state, drawing, fileStore);
  if (drawing.mode !== 'linear') return files;
  return {
    ...files,
    linearCanonicalComparisons: (
      Array.isArray(state.files?.linearCanonicalComparisons)
        ? state.files.linearCanonicalComparisons
        : []
    ).map((/** @type {ResourceBackedComparison} */ comparison) => (
      isResourceBackedCanonicalComparison(comparison)
        ? {
            ...mapResourceBackedCanonicalComparison(comparison),
            file: fileStore.describeValue(comparison.file)
          }
        : cloneJsonData(comparison)
    ))
  };
};

/**
 * The input bindings of the app state.
 * @param {HistorySnapshotState} state
 * @returns {InputBindings}
 */
const inputBindingsOf = (state) => ({ files: state.files, linearSeqs: state.linearSeqs, drawingList: state.drawingList });

/**
 * @param {Record<string, any>} state
 * @param {DrawingView} drawing
 * @param {FilesData} filesData
 * @param {HistorySnapshotFileStore} fileStore
 * @param {HistorySnapshotServiceOptions['normalizeLinearSeqList']} [normalizeLinearSeqList]
 */
const applyFilesData = (state, drawing, filesData, fileStore, normalizeLinearSeqList = null) => {
  if (!state.files) return;
  state.matchSequenceRegistry?.reset?.();
  /** @param {unknown} value */
  const restore = (value) => fileStore.restoreValue(value);
  // A table card holds its file's ref: the step's ref, kept bound by History.
  FILE_IMPORT_CARDS.forEach((key) => {
    const value = filesData?.[/** @type {keyof FilesData} */ (key)];
    state.files[key] = isInputFileRef(value) ? fileStore.registerFile(value) : null;
  });
  if (drawing.mode !== 'linear') {
    state.files.c_gb = restore(filesData?.c_gb);
    state.files.c_gff = restore(filesData?.c_gff);
    state.files.c_fasta = restore(filesData?.c_fasta);
    state.files.c_depth = restore(filesData?.c_depth);
    state.files.c_conservation_blasts = Array.isArray(filesData?.c_conservation_blasts)
      ? restore(filesData.c_conservation_blasts).filter(Boolean)
      : [];
    state.files.c_conservation_blasts_source = filesData?.c_conservation_blasts_source === 'losat-cache'
      ? 'losat-cache'
      : null;
    state.files.c_conservation_fastas = Array.isArray(filesData?.c_conservation_fastas)
      ? restore(filesData.c_conservation_fastas)
      : [];
    state.files.c_conservation_sequence_sources = Array.isArray(filesData?.c_conservation_sequence_sources)
      ? restore(filesData.c_conservation_sequence_sources)
      : [];
    return;
  }
  if (Object.prototype.hasOwnProperty.call(filesData || {}, 'linearCanonicalComparisons')) {
    state.files.linearCanonicalComparisons = Array.isArray(filesData?.linearCanonicalComparisons)
      ? filesData.linearCanonicalComparisons.map((comparison) => (
          isResourceBackedCanonicalComparison(comparison)
            ? mapResourceBackedCanonicalComparison(comparison, restore)
            : cloneJsonData(comparison)
        ))
      : [];
  }

  if (!state.linearSeqs || typeof state.linearSeqs.splice !== 'function') return;
  const rows = Array.isArray(filesData?.linearSeqs)
    ? filesData.linearSeqs.map((seq) => ({
        uid: seq.uid,
        gb: restore(seq.gb),
        gff: restore(seq.gff),
        fasta: restore(seq.fasta),
        depth: restore(seq.depth),
        losat_gencode: seq.losat_gencode ?? 1,
        definition: seq.definition ?? '',
        record_subtitle: seq.record_subtitle ?? '',
        file_definition: seq.file_definition ?? '',
        file_subtitle: seq.file_subtitle ?? '',
        inferred_definition: seq.inferred_definition ?? '',
        region_record_id: seq.region_record_id ?? '',
        region_start: seq.region_start ?? null,
        region_end: seq.region_end ?? null,
        region_reverse: Boolean(seq.region_reverse)
      }))
    : [];
  const normalized = typeof normalizeLinearSeqList === 'function'
    ? normalizeLinearSeqList(rows)
    : rows;
  state.linearSeqs.splice(0, state.linearSeqs.length, ...normalized);
  if (drawing.linearComparisonPlan && Array.isArray(drawing.linearComparisonPlan.edges)) {
    const comparisonFiles = new Map(
      (Array.isArray(filesData?.linearComparisons) ? filesData.linearComparisons : [])
        .map((comparison) => [String(comparison?.id || ''), restore(comparison?.file)])
    );
    drawing.linearComparisonPlan.edges.forEach((edge) => {
      edge.file = comparisonFiles.get(String(edge?.id || '')) || null;
    });
  }
};

/**
 * @typedef {object} HistorySnapshotFileStore
 *   The functions of the project's input file store (`state.inputFiles`)
 *   that the snapshot service calls (F-05).
 * @property {(value: any) => any} describeValue
 *   A file, or an array of files, as the refs an intent holds (`null` for a non-file).
 * @property {InputFileStore['registerFile']} registerFile
 *   Stores a file by reference and returns its one ref; `null` for a non-file.
 * @property {(value: any) => any} restoreValue
 *   The stored file, or an array of them, that the refs name (`null` when none).
 *
 * @typedef {object} HistorySnapshotRuntimeOwner
 *   The generated-artifact runtime owner that the composition root registers (`app/app-setup.js`).
 * @property {() => Record<string, any>} capture The runtime state, with its canonical session.
 * @property {(snapshot: any, options?: any) => unknown} restore Installs a captured runtime state.
 *
 * @typedef {object} HistorySnapshotTransportIdentity
 *   The frozen identity of the generated artifact that a transport carries.
 * @property {1} schema
 * @property {'SHA-256'} algorithm
 * @property {string} fingerprint
 * @property {Readonly<Record<string, any>>} ownerReferences
 *
 * @typedef {Record<string, any>} HistorySnapshotData
 *   A domain payload that the snapshot service holds by value (config, ui, features, editor state,
 *   orthogroup state, run state): the owner that builds it declares its shape.
 */

/**
 * The ports take the drawing type `D` the composition root's drawings have
 * (state.js `DrawingState`); the service reads them as a `DrawingView`.
 * @template {DrawingView} [D=DrawingView]
 * @typedef {object} HistorySnapshotServiceOptions
 * @property {Record<string, any> & { drawingOf: (id: DrawingId) => D, drawingList: () => readonly D[] }} state App state (state.js; not yet typed).
 * @property {HistorySnapshotFileStore} fileStore
 *   The file store that holds the files an intent and a checkpoint name.
 * @property {() => Promise<unknown>} [nextTick] Vue `nextTick`, awaited between restore steps.
 * @property {((rows: HistorySnapshotData[]) => HistorySnapshotData[]) | null} [normalizeLinearSeqList]
 *   The Linear record rows' normalization after a files restore.
 * @property {((drawing: D) => HistorySnapshotData) | null} [buildConfigData] The Settings capture.
 * @property {((drawing: D, config: HistorySnapshotData, options: { resolveTrackPlacements: boolean, applyAlignment?: boolean }) => unknown) | null} [applyConfigData]
 *   The Settings restore.
 * @property {((drawing: D, options: { includePreviewNavigation: boolean }) => HistorySnapshotData) | null} [buildUiStateData]
 *   The UI capture; `includePreviewNavigation` adds pan, zoom, and tab state.
 * @property {((drawing: D, ui: HistorySnapshotData, options?: { restorePreviewNavigation: boolean, artifact?: boolean }) => unknown) | null} [applyUiStateData]
 *   The UI restore.
 * @property {((ui: SessionUiData) => unknown) | null} [applySessionUiData]
 *   The restore of the Session-level `ui` values, which writes each one present.
 * @property {((drawing: D) => HistorySnapshotData) | null} [buildFeatureStateData] The feature-override capture.
 * @property {((drawing: D, features: HistorySnapshotData) => unknown) | null} [applyFeatureStateData]
 *   The feature-override restore.
 * @property {((drawing: D) => HistorySnapshotData) | null} [buildEditorStateData] The legend and stroke editor capture.
 * @property {((drawing: D, editorState: HistorySnapshotData, options: { normalized: boolean }) => unknown) | null} [applyEditorStateData]
 *   The editor restore; `normalized` marks a state that is installed as captured.
 * @property {((drawing: D) => HistorySnapshotData) | null} [buildOrthogroupStateData] The orthogroup capture.
 * @property {((drawing: D, orthogroupState: HistorySnapshotData) => unknown) | null} [applyOrthogroupStateData]
 *   The orthogroup restore.
 * @property {(() => HistorySnapshotData[]) | null} [serializeResults] The generated Results' capture.
 * @property {((results: HistorySnapshotData[], ui: HistorySnapshotData) => unknown) | null} [applyResultsData]
 *   The generated Results' restore, with the checkpoint's UI state.
 * @property {(() => HistorySnapshotData) | null} [buildRunStateData] The last-run capture.
 * @property {((runState: HistorySnapshotData) => unknown) | null} [applyRunStateData] The last-run restore.
 * @property {(() => void) | null} [clearSimilarityGroupCopyFeedback]
 *   Clears the Similarity groups "Copied!" feedback when a Generate's artifact replaces the shown
 *   Results (D-37): a reused feature catalog keeps the groups array, so its watcher would not fire.
 */

/**
 * @typedef {Readonly<{ featureCatalog: unknown, extractedFeatures: readonly unknown[],
 *   biologicalFeatures: readonly unknown[], results: unknown[], labelFitReports: readonly unknown[] }>} RedrawnArtifact
 *   What a label rerender draws again of its Result's drawing (OV-111).
 */

/**
 * @template {DrawingView} D
 * @param {HistorySnapshotServiceOptions<D>} options
 */
export const createHistorySnapshotService = ({
  state,
  fileStore,
  nextTick = async () => {},
  normalizeLinearSeqList = null,
  buildConfigData = null,
  applyConfigData = null,
  buildUiStateData = null,
  applyUiStateData = null,
  applySessionUiData = null,
  buildFeatureStateData = null,
  applyFeatureStateData = null,
  buildEditorStateData = null,
  applyEditorStateData = null,
  buildOrthogroupStateData = null,
  applyOrthogroupStateData = null,
  serializeResults = null,
  applyResultsData = null,
  buildRunStateData = null,
  applyRunStateData = null,
  clearSimilarityGroupCopyFeedback = null
}) => {
  if (!state || !fileStore) {
    throw new Error('createHistorySnapshotService requires state and fileStore.');
  }

  // A generated artifact holds the settings and edits of the drawing of its
  // mode (`generatedMode`, Circular when unset): the one reading of that mode.
  /**
   * @param {unknown} generatedMode
   * @returns {D}
   */
  const drawingOfGeneratedMode = (generatedMode) => state.drawingOfMode(generatedMode === 'linear' ? 'linear' : 'circular');
  // The drawing the displayed artifact belongs to: the drawing of its mode
  // (`generatedMode`). It differs from the shown drawing after a Load of a
  // Session with no Result whose request mode is not its `ui.mode`. Every
  // capture of the displayed artifact names its drawing so.
  const displayedArtifactDrawing = () => drawingOfGeneratedMode(getGeneratedArtifactRef(state.generatedMode, null));

  /** @type {((intent: HistorySnapshotData, context: Record<string, any>) => unknown) | null} */
  let afterApplyHistoryIntent = null;
  /** @type {HistorySnapshotRuntimeOwner | null} */
  let generatedArtifactRuntimeOwner = null;
  /** @type {HistorySnapshotTransportIdentity | null} */
  let currentGeneratedArtifactIdentity = null;
  let currentGeneratedArtifactRetainedBytes = 0;
  let generatedArtifactRestoreDepth = 0;
  let restoreSemanticSuppressionBaseline = false;
  let restoreTrustedStateBaseline = false;

  /** @param {typeof afterApplyHistoryIntent} callback */
  const setAfterApplyHistoryIntent = (callback) => {
    afterApplyHistoryIntent = typeof callback === 'function' ? callback : null;
  };

  // Intent captures of owners created after this service (R13): the
  // composition root registers each once its owner exists. An unregistered
  // capture adds nothing to the intent, as for a service without that owner.
  /** @type {Record<string, ((drawingId: DrawingId) => unknown) | null>} */
  const captures = { composition: null };
  /**
   * @param {string} name The capture's name; an unknown name throws.
   * @param {((drawingId: DrawingId) => unknown) | null} capture The capture of one drawing's part.
   */
  const registerCapture = (name, capture) => {
    if (!Object.hasOwn(captures, name)) throw new Error(`Unknown History intent capture: ${name}`);
    captures[name] = typeof capture === 'function' ? capture : null;
  };

  /** @param {HistorySnapshotRuntimeOwner | null} owner */
  const setGeneratedArtifactRuntimeOwner = (owner) => {
    generatedArtifactRuntimeOwner = (
      owner
      && typeof owner.capture === 'function'
      && typeof owner.restore === 'function'
    ) ? owner : null;
  };

  const captureLinearRecordOrientations = () => (
    (Array.isArray(state.linearSeqs) ? state.linearSeqs : []).map((sequence) => ({
      recordKey: String(sequence?.uid || ''),
      reverseComplement: Boolean(sequence?.region_reverse)
    })).filter(({ recordKey }) => recordKey)
  );

  /** @param {unknown} orientations */
  const installLinearRecordOrientations = (orientations) => {
    const byRecord = new Map(
      (Array.isArray(orientations) ? orientations : []).map((entry) => [
        String(entry?.recordKey || ''),
        Boolean(entry?.reverseComplement)
      ])
    );
    (Array.isArray(state.linearSeqs) ? state.linearSeqs : []).forEach((sequence) => {
      const recordKey = String(sequence?.uid || '');
      if (byRecord.has(recordKey)) sequence.region_reverse = byRecord.get(recordKey);
    });
  };

  /** @typedef {ReturnType<typeof captureGeneratedArtifactOwnerSet>} GeneratedArtifactOwnerSet The owners' values, by reference, key by key. */
  const captureGeneratedArtifactOwnerSet = () => {
    const drawing = displayedArtifactDrawing();
    const results = toRaw(getGeneratedArtifactRef(state.results, []));
    const colorRules = COLOR_RULES_SLICE.capture(drawing.slices.colorRules);
    return Object.freeze({
      results,
      featureCatalog: toRaw(getGeneratedArtifactRef(state.featureCatalog, null)),
      extractedFeatures: toRaw(getGeneratedArtifactRef(state.extractedFeatures, null)),
      biologicalFeatures: toRaw(getGeneratedArtifactRef(state.biologicalFeatures, null)),
      featureRecordIds: toRaw(getGeneratedArtifactRef(state.featureRecordIds, null)),
      orthogroups: toRaw(getGeneratedArtifactRef(state.orthogroups, null)),
      featureOrthogroupIndex: toRaw(
        getGeneratedArtifactRef(state.featureOrthogroupIndex, null)
      ),
      collinearGroups: toRaw(getGeneratedArtifactRef(state.collinearGroups, null)),
      similarityAlignmentResetReceipt: toRaw(
        getGeneratedArtifactRef(state.similarityAlignmentResetReceipt, null)
      ),
      similarityAlignmentPlan: toRaw(
        getGeneratedArtifactRef(state.similarityAlignmentPlan, null)
      ),
      linearRecordTranslations: toRaw(
        getGeneratedArtifactRef(state.linearRecordTranslations, null)
      ),
      linearRecordOrientations: toRaw(captureLinearRecordOrientations()),
      annotationWarnings: toRaw(getGeneratedArtifactRef(state.annotationWarnings, null)),
      featureIdentityNotices: toRaw(
        getGeneratedArtifactRef(state.featureIdentityNotices, null)
      ),
      comparisonWarnings: toRaw(getGeneratedArtifactRef(state.comparisonWarnings, null)),
      labelFitReports: toRaw(getGeneratedArtifactRef(state.labelFitReports, null)),
      specificRules: colorRules.manualSpecificRules,
      trackSlotResolvedGeometry: toRaw(
        getGeneratedArtifactRef(state.trackSlotResolvedGeometry, null)
      ),
      proteinIdentityManifest: toRaw(
        getGeneratedArtifactRef(state.proteinIdentityManifest, null)
      ),
      legacyProteinRawCandidates: toRaw(
        getGeneratedArtifactRef(state.legacyProteinRawCandidates, null)
      ),
      legacyProteinDerivedEvidence: toRaw(
        getGeneratedArtifactRef(state.legacyProteinDerivedEvidence, null)
      ),
      losatCache: toRaw(getGeneratedArtifactRef(state.losatCache, null)),
      losatDerivedCache: toRaw(
        getGeneratedArtifactRef(state.losatDerivedCache, null)
      ),
      losatCacheInfo: toRaw(getGeneratedArtifactRef(state.losatCacheInfo, null)),
      matchSequenceOwner: state.matchSequenceRegistry?.captureTrustedOwner?.() || null,
      lastRunInfo: toRaw(getGeneratedArtifactRef(state.lastRunInfo, null)),
      pairwiseMatchFactors: toRaw(
        getGeneratedArtifactRef(state.pairwiseMatchFactors, null)
      ),
      editableLabels: toRaw(getGeneratedArtifactRef(state.editableLabels, null)),
      generatedLegendPosition: getGeneratedArtifactRef(state.generatedLegendPosition, null),
      generatedMode: getGeneratedArtifactRef(state.generatedMode, null),
      generatedMultiRecordCanvas: getGeneratedArtifactRef(
        state.generatedMultiRecordCanvas,
        false
      ),
      generatedCircularPlotTitlePosition: getGeneratedArtifactRef(
        state.generatedCircularPlotTitlePosition,
        null
      ),
      appliedPaletteName: getGeneratedArtifactRef(state.appliedPaletteName, ''),
      appliedPaletteColors: toRaw(
        getGeneratedArtifactRef(state.appliedPaletteColors, null)
      ),
      // By reference, as every artifact value (`restoreArtifactQueuedPalette`).
      pendingPaletteName: drawing.pendingPaletteName.value,
      pendingPaletteColors: toRaw(drawing.pendingPaletteColors.value)
    });
  };

  const installGeneratedArtifactOwnerSet = (
    /** @type {Partial<GeneratedArtifactOwnerSet> | null | undefined} */
    ownerSet,
    /** @type {{ selectedResultIndex?: number, installResults?: ((results: any[]) => unknown) | null }} */
    { selectedResultIndex = 0, installResults = null } = {}
  ) => {
    if (!ownerSet || typeof ownerSet !== 'object') {
      throw new Error('A generated artifact owner set is required.');
    }
    const drawing = drawingOfGeneratedMode(ownerSet.generatedMode);
    const results = Array.isArray(ownerSet.results) ? ownerSet.results : [];
    const replacesShownResults = toRaw(results)
      !== toRaw(getGeneratedArtifactRef(state.results, []));
    if (typeof installResults === 'function') installResults(results);
    else setGeneratedArtifactRef(state.results, results);
    setGeneratedArtifactRef(
      state.selectedResultIndex,
      results.length > 0
        ? Math.max(0, Math.min(Number(selectedResultIndex) || 0, results.length - 1))
        : 0
    );
    setGeneratedArtifactRef(state.featureCatalog, ownerSet.featureCatalog ?? null);
    setGeneratedArtifactRef(state.extractedFeatures, ownerSet.extractedFeatures || []);
    setGeneratedArtifactRef(state.biologicalFeatures, ownerSet.biologicalFeatures || []);
    setGeneratedArtifactRef(state.featureRecordIds, ownerSet.featureRecordIds || []);
    setGeneratedArtifactRef(state.orthogroups, ownerSet.orthogroups || []);
    setGeneratedArtifactRef(
      state.featureOrthogroupIndex,
      ownerSet.featureOrthogroupIndex || new Map()
    );
    setGeneratedArtifactRef(state.collinearGroups, ownerSet.collinearGroups || []);
    setGeneratedArtifactRef(
      state.similarityAlignmentResetReceipt,
      ownerSet.similarityAlignmentResetReceipt ?? null
    );
    setGeneratedArtifactRef(
      state.similarityAlignmentPlan,
      ownerSet.similarityAlignmentPlan ?? null
    );
    setGeneratedArtifactRef(
      state.linearRecordTranslations,
      ownerSet.linearRecordTranslations || []
    );
    installLinearRecordOrientations(ownerSet.linearRecordOrientations);
    setGeneratedArtifactRef(state.annotationWarnings, ownerSet.annotationWarnings || []);
    setGeneratedArtifactRef(state.featureIdentityNotices, ownerSet.featureIdentityNotices || []);
    setGeneratedArtifactRef(state.comparisonWarnings, ownerSet.comparisonWarnings || []);
    setGeneratedArtifactRef(state.labelFitReports, ownerSet.labelFitReports || []);
    COLOR_RULES_SLICE.install(drawing.slices.colorRules, { manualSpecificRules: ownerSet.specificRules });
    setGeneratedArtifactRef(
      state.trackSlotResolvedGeometry,
      ownerSet.trackSlotResolvedGeometry ?? null
    );
    setGeneratedArtifactRef(state.proteinIdentityManifest, ownerSet.proteinIdentityManifest);
    setGeneratedArtifactRef(
      state.legacyProteinRawCandidates,
      ownerSet.legacyProteinRawCandidates
    );
    setGeneratedArtifactRef(
      state.legacyProteinDerivedEvidence,
      ownerSet.legacyProteinDerivedEvidence
    );
    setGeneratedArtifactRef(state.losatCache, ownerSet.losatCache || new Map());
    setGeneratedArtifactRef(state.losatDerivedCache, ownerSet.losatDerivedCache || new Map());
    setGeneratedArtifactRef(state.losatCacheInfo, ownerSet.losatCacheInfo || []);
    if (ownerSet.matchSequenceOwner) {
      state.matchSequenceRegistry?.replaceTrustedOwner?.(ownerSet.matchSequenceOwner);
    }
    setGeneratedArtifactRef(state.lastRunInfo, ownerSet.lastRunInfo ?? null);
    setGeneratedArtifactRef(state.pairwiseMatchFactors, ownerSet.pairwiseMatchFactors || {});
    setGeneratedArtifactRef(state.editableLabels, ownerSet.editableLabels || []);
    setGeneratedArtifactRef(state.generatedLegendPosition, ownerSet.generatedLegendPosition);
    setGeneratedArtifactRef(state.generatedMode, ownerSet.generatedMode);
    setGeneratedArtifactRef(
      state.generatedMultiRecordCanvas,
      Boolean(ownerSet.generatedMultiRecordCanvas)
    );
    setGeneratedArtifactRef(
      state.generatedCircularPlotTitlePosition,
      ownerSet.generatedCircularPlotTitlePosition
    );
    setGeneratedArtifactRef(state.appliedPaletteName, ownerSet.appliedPaletteName || '');
    setGeneratedArtifactRef(state.appliedPaletteColors, ownerSet.appliedPaletteColors || {});
    restoreArtifactQueuedPalette(drawing.slices.palette, ownerSet.pendingPaletteName || '', ownerSet.pendingPaletteColors || {});
    if (replacesShownResults) clearSimilarityGroupCopyFeedback?.();
  };

  // OV-111: a label rerender's Results and the catalog that describes them
  // install into the artifact of `drawingId`, the drawing whose Result they
  // draw again. A rerender draws the displayed artifact (a switch waits for
  // it, OIPC-C07), so the artifact a slot keeps is refused. The caller commits
  // the Results (R1).
  /**
   * @param {DrawingId} drawingId
   * @param {RedrawnArtifact} redrawn
   * @param {{ installResults: (results: unknown[]) => void }} options
   */
  const installRedrawnArtifact = (drawingId, redrawn, { installResults }) => {
    const displayed = displayedArtifactDrawing();
    if (displayed.id !== drawingId) {
      throw diagnosticError('RESULT_INVALID', { reason: 'DRAWING_MISMATCH', drawingId, artifactDrawingId: displayed.id }, { stage: 'result-admission', operation: 'installRedrawnArtifact' });
    }
    setGeneratedArtifactRef(state.featureCatalog, redrawn.featureCatalog);
    setGeneratedArtifactRef(state.extractedFeatures, redrawn.extractedFeatures);
    setGeneratedArtifactRef(state.biologicalFeatures, redrawn.biologicalFeatures);
    installResults(redrawn.results);
    setGeneratedArtifactRef(state.labelFitReports, redrawn.labelFitReports);
  };

  // E1: the displayed artifact, by reference, key by key, as the artifact of
  // the drawing it belongs to (`displayedArtifactDrawing`).
  const captureArtifactSlot = () => createArtifactSlot({
    drawing: displayedArtifactDrawing(),
    values: Object.fromEntries(ARTIFACT_SLOT_KEYS.map((key) => [
      key, toRaw(getGeneratedArtifactRef(state[key], null))
    ])),
    legendInventory: getGeneratedArtifactRef(state.originalLegendOrder, []) || [],
    matchSequenceOwner: state.matchSequenceRegistry?.captureTrustedOwner?.() || null,
    runtimeState: generatedArtifactRuntimeOwner?.capture?.() || null,
    transportIdentity: currentGeneratedArtifactIdentity,
    retainedBytes: currentGeneratedArtifactRetainedBytes
  });

  // Installs `drawing`'s slot, or its empty slot (no Result, no committed
  // Session, the defaults of the drawing's mode), as the displayed artifact;
  // a null value of the slot installs the empty slot's. Nothing outside the
  // slot changes, so the draft, the editor intent and the caches stay shared.
  /**
   * @param {Readonly<ArtifactSlot> | null} slot
   * @param {{ drawing: DrawingIdentity, ui?: { canvasPan?: { x: number, y: number } } }} options
   *   `ui.canvasPan` is the pan the runtime owner's viewport reset keeps.
   */
  const installArtifactSlot = (slot, { drawing, ui = {} }) => {
    if (slot && slot.drawingId !== drawing.id) {
      throw diagnosticError('RESULT_INVALID', { reason: 'DRAWING_MISMATCH', drawingId: drawing.id, artifactDrawingId: slot.drawingId }, { stage: 'result-admission', operation: 'installArtifactSlot' });
    }
    ARTIFACT_SLOT_KEYS.forEach((key) => {
      setGeneratedArtifactRef(state[key], slot?.values[key]
        ?? ARTIFACT_SLOT_EMPTY[key](drawing.mode, getGeneratedArtifactRef(state[key], null)));
    });
    if (slot?.matchSequenceOwner) state.matchSequenceRegistry?.replaceTrustedOwner?.(slot.matchSequenceOwner);
    else state.matchSequenceRegistry?.reset?.();
    generatedArtifactRuntimeOwner?.restore?.(slot?.runtimeState ?? null, { ui });
    currentGeneratedArtifactIdentity = /** @type {HistorySnapshotTransportIdentity | null} */ (
      slot?.transportIdentity ?? null
    );
    currentGeneratedArtifactRetainedBytes = slot?.retainedBytes ?? 0;
  };

  /**
   * @param {any} identity
   * @param {{ results?: any[] | null }} [options]
   */
  const setGeneratedArtifactIdentity = (identity, { results = null } = {}) => {
    const fingerprint = String(identity?.fingerprint || '').toLowerCase();
    if (
      identity?.schema !== 1
      || identity?.algorithm !== 'SHA-256'
      || !/^[0-9a-f]{64}$/.test(fingerprint)
    ) {
      currentGeneratedArtifactIdentity = null;
      currentGeneratedArtifactRetainedBytes = 0;
      return false;
    }
    const committedResults = Array.isArray(results)
      ? results
      : getGeneratedArtifactRef(state.results, []);
    const ownerReferences = captureGeneratedArtifactOwnerSet();
    currentGeneratedArtifactIdentity = Object.freeze({
      schema: 1,
      algorithm: 'SHA-256',
      fingerprint,
      ownerReferences: Object.freeze({
        ...ownerReferences,
        results: toRaw(committedResults)
      })
    });
    currentGeneratedArtifactRetainedBytes = Math.max(
      0,
      (Number(identity.retainedBytes) || 0) + (Number(identity.resultBytes) || 0)
    );
    return true;
  };

  const clearGeneratedArtifactIdentity = ({ retainedBytes = 0 } = {}) => {
    currentGeneratedArtifactIdentity = null;
    currentGeneratedArtifactRetainedBytes = Math.max(0, Number(retainedBytes) || 0);
  };

  /**
   * @param {Readonly<Record<string, unknown>>} ownerSet
   * @param {HistorySnapshotTransportIdentity | null} identity
   */
  const currentArtifactMatchesIdentity = (ownerSet, identity) => {
    const expected = identity?.ownerReferences;
    if (!expected || expected.results !== ownerSet.results) return false;
    return [
      'featureCatalog',
      'extractedFeatures',
      'biologicalFeatures',
      'featureRecordIds',
      'orthogroups',
      'featureOrthogroupIndex',
      'collinearGroups',
      'similarityAlignmentResetReceipt',
      'similarityAlignmentPlan',
      'linearRecordTranslations',
      'trackSlotResolvedGeometry',
      'annotationWarnings',
      'featureIdentityNotices',
      'comparisonWarnings',
      'labelFitReports',
      'proteinIdentityManifest',
      'legacyProteinRawCandidates',
      'legacyProteinDerivedEvidence',
      'losatCache',
      'losatDerivedCache',
      'losatCacheInfo',
      'matchSequenceOwner',
      'lastRunInfo',
      'pairwiseMatchFactors',
      'editableLabels',
      'generatedLegendPosition',
      'generatedMode',
      'generatedMultiRecordCanvas',
      'generatedCircularPlotTitlePosition'
    ].every((key) => expected[key] === ownerSet[key]);
  };

  /** @param {GeneratedArtifactIntent} mutableIntent */
  const compactGeneratedArtifactSignature = (mutableIntent) => {
    const {
      zoom: _zoom,
      canvasPan: _canvasPan,
      generatedLegendPosition: _generatedLegendPosition,
      generatedMode: _generatedMode,
      generatedMultiRecordCanvas: _generatedMultiRecordCanvas,
      generatedCircularPlotTitlePosition: _generatedCircularPlotTitlePosition,
      ...ui
    } = mutableIntent.ui || {};
    return JSON.stringify({
      ui,
      featureEdits: {
        selectedFeatureRecordIdx: mutableIntent.features.selectedFeatureRecordIdx,
        featureColorOverrides: mutableIntent.features.featureColorOverrides,
        featureVisibilityManualRules: mutableIntent.features.featureVisibilityManualRules,
        featureOverrides: mutableIntent.features.featureOverrides,
        labelOverrideRows: mutableIntent.features.labelOverrideRows,
        labelTextBulkOverrides: mutableIntent.features.labelTextBulkOverrides
      },
      editor: {
        legend: mutableIntent.editorState.legend,
        featureStrokes: mutableIntent.editorState.featureStrokes,
        originalSvgStroke: mutableIntent.editorState.originalSvgStroke
      },
      orthogroupEdits: {
        selectedOrthogroupId: mutableIntent.orthogroupState.selectedOrthogroupId,
        selectedOrthogroupAlignmentFeature:
          mutableIntent.orthogroupState.selectedOrthogroupAlignmentFeature,
        orthogroupNameOverrides: mutableIntent.orthogroupState.orthogroupNameOverrides,
        orthogroupDescriptionOverrides:
          mutableIntent.orthogroupState.orthogroupDescriptionOverrides,
        orthogroupDormantOverrides: mutableIntent.orthogroupState.orthogroupDormantOverrides
      },
      alignmentState: mutableIntent.alignmentState,
      presentation: {
        resultPanelTab: mutableIntent.presentation.resultPanelTab
      }
    });
  };

  // The handle holds the edits of the drawing whose step or Generate captures
  // it (the shown one), and its restore installs them in that drawing.
  /** @param {D} drawing */
  const captureGeneratedArtifactHandle = (drawing) => {
    recordStructuralMetric('generatedArtifactFullSignatureCount', 0);
    recordStructuralMetric('generatedArtifactHeavyTraversalCount', 0);
    recordStructuralMetric('generatedArtifactMutableIntentSnapshotCount', 1);
    const ownerSet = captureGeneratedArtifactOwnerSet();
    const ui = typeof buildUiStateData === 'function'
      ? buildUiStateData(drawing, { includePreviewNavigation: true })
      : buildFallbackUiStateData(state, drawing);
    const features = {
      selectedFeatureRecordIdx: getGeneratedArtifactRef(state.selectedFeatureRecordIdx, 0),
      featureColorOverrides: captureFeatureFills(drawing.slices.paint),
      featureVisibilityManualRules: cloneFeatureVisibilityRules(
        VISIBILITY_RULES_SLICE.capture(drawing.slices.visibilityRules).featureVisibilityManualRules
      ),
      featureOverrides: FEATURE_EDITS_SLICE.capture(drawing.slices.featureEdits).featureOverrides,
      ...labelTextEditsData(drawing.slices.labelTables)
    };
    const paint = captureLegendPaint(drawing.slices.paint);
    const editorState = {
      legend: {
        entries: paint.legendEntries,
        deletedEntries: paint.deletedLegendEntries,
        dormantEntries: paint.dormantLegendEntries,
        originalOrder: cloneJsonData(
          getGeneratedArtifactRef(state.originalLegendOrder, [])
        ) || [],
        originalColors: clonePlainObject(
          getGeneratedArtifactRef(state.originalLegendColors, {})
        ),
        colorOverrides: paint.legendColorOverrides,
        strokeOverrides: paint.legendStrokeOverrides,
        addedCaptions: paint.addedLegendCaptions
      },
      featureStrokes: { overrides: paint.featureStrokeOverrides },
      originalSvgStroke: cloneJsonData(
        getGeneratedArtifactRef(state.originalSvgStroke, null)
      ) || {
        color: null,
        width: null
      }
    };
    const orthogroupState = {
      selectedOrthogroupId: String(
        getGeneratedArtifactRef(state.selectedOrthogroupId, '') || ''
      ),
      selectedOrthogroupAlignmentFeature: String(
        getGeneratedArtifactRef(state.selectedOrthogroupAlignmentFeature, '') || ''
      ),
      ...buildGroupNameIntentData(drawing)
    };
    const mutableIntent = Object.freeze({
      ui,
      features,
      editorState,
      orthogroupState,
      alignmentState: Object.freeze({
        receipt: cloneJsonData(getGeneratedArtifactRef(state.similarityAlignmentResetReceipt, null)),
        plan: cloneJsonData(getGeneratedArtifactRef(state.similarityAlignmentPlan, null)),
        recordTranslations: cloneJsonData(
          getGeneratedArtifactRef(state.linearRecordTranslations, [])
        ) || [],
        recordOrientations: captureLinearRecordOrientations()
      }),
      presentation: Object.freeze({
        resultGenerationKey: getGeneratedArtifactRef(state.resultGenerationKey, 0),
        resultPanelTab: getGeneratedArtifactRef(state.resultPanelTab, 'preview'),
        errorLog: cloneJsonData(getGeneratedArtifactRef(state.errorLog, null)),
        labelTextScopeDialog: cloneJsonData(state.labelTextScopeDialog || {}) || {},
        featureEditorStatus: cloneJsonData(state.featureEditorStatus || {}) || {},
        featureExtractionPending: Boolean(
          getGeneratedArtifactRef(state.featureExtractionPending, false)
        ),
        featureExtractionError: getGeneratedArtifactRef(state.featureExtractionError, null)
      })
    });
    const runtimeState = generatedArtifactRuntimeOwner?.capture?.() || null;
    const compactSignature = compactGeneratedArtifactSignature(mutableIntent);
    const transportIdentity = currentArtifactMatchesIdentity(
      ownerSet,
      currentGeneratedArtifactIdentity
    ) ? currentGeneratedArtifactIdentity : null;
    return Object.freeze({
      kind: 'GeneratedArtifactHandle',
      drawingId: drawing.id,
      ownerSet,
      mutableIntent,
      runtimeState,
      transportIdentity,
      baseRetainedBytes: currentGeneratedArtifactRetainedBytes,
      fileIds: Object.freeze([...boundInputFileIds(inputBindingsOf(state), fileStore)]),
      identity: Object.freeze({
        fingerprint: transportIdentity?.fingerprint || '',
        compactSignature
      }),
      retainedBytes: currentGeneratedArtifactRetainedBytes
        + Math.max(0, Number(runtimeState?.retainedBytes) || 0)
        + compactSignature.length * 2
    });
  };

  /**
   * The handle History stored and passes back: History's own type names the fields History reads; these are
   * the fields `captureGeneratedArtifactHandle` wrote beside them, as a restore reads them.
   * @typedef {Pick<HistoryArtifactHandle, 'retainedBytes' | 'fileIds' | 'identity'> & {
   *   kind?: unknown,
   *   drawingId?: DrawingId,
   *   ownerSet?: GeneratedArtifactOwnerSet | null,
   *   mutableIntent?: Partial<GeneratedArtifactIntent>,
   *   runtimeState?: { canonical?: { committedCanonicalSession?: AlignmentResetCanonical | null } | null } | null,
   *   transportIdentity?: HistorySnapshotTransportIdentity | null,
   *   baseRetainedBytes?: unknown
   * }} StoredGeneratedArtifactHandle
   */
  const restoreGeneratedArtifactHandle = async (
    /** @type {StoredGeneratedArtifactHandle | null | undefined} */
    handle,
    { clearFailedGeneratePresentation = false } = {}
  ) => {
    if (!handle || handle.kind !== 'GeneratedArtifactHandle' || !handle.drawingId) return false;
    const { drawingId } = handle;
    await validateSimilarityAlignmentResetReceipt(
      handle.ownerSet?.similarityAlignmentResetReceipt,
      handle.runtimeState?.canonical?.committedCanonicalSession
    );
    if (clearFailedGeneratePresentation && state.failedGeneratePreservedResult) {
      state.failedGeneratePreservedResult.value = false;
    }
    if (generatedArtifactRestoreDepth === 0) {
      restoreSemanticSuppressionBaseline = Boolean(
        getGeneratedArtifactRef(state.semanticFileWatchersSuppressed, false)
      );
      restoreTrustedStateBaseline = Boolean(
        getGeneratedArtifactRef(state.trustedArtifactRestoreInProgress, false)
      );
    }
    generatedArtifactRestoreDepth += 1;
    setGeneratedArtifactRef(state.semanticFileWatchersSuppressed, true);
    setGeneratedArtifactRef(state.trustedArtifactRestoreInProgress, true);
    try {
      closeTransientState(state);
      const mutableIntent = handle.mutableIntent || {};
      const ui = mutableIntent.ui || {};
      const drawing = state.drawingOf(drawingId);
      const inputType = inputTypeKeyOf(drawing);
      if (ui[inputType]) setGeneratedArtifactRef(state[inputType], ui[inputType]);
      await nextTick();

      if (state.skipCaptureBaseConfig) state.skipCaptureBaseConfig.value = true;
      installGeneratedArtifactOwnerSet(handle.ownerSet, {
        selectedResultIndex: ui.selectedResultIndex
      });
      applyFeatureIntentData(state, drawing, mutableIntent.features || {});
      setGeneratedArtifactRef(
        state.selectedOrthogroupId,
        String(mutableIntent.orthogroupState?.selectedOrthogroupId || '')
      );
      setGeneratedArtifactRef(
        state.selectedOrthogroupAlignmentFeature,
        String(mutableIntent.orthogroupState?.selectedOrthogroupAlignmentFeature || '')
      );
      applyGroupNameIntentData(drawing, mutableIntent.orthogroupState);

      const trustedEditorState = {
        ...cloneJsonData({
          legend: mutableIntent.editorState?.legend || {},
          featureStrokes: mutableIntent.editorState?.featureStrokes || {},
          originalSvgStroke: mutableIntent.editorState?.originalSvgStroke || {}
        }),
        alignmentResetReceipt: handle.ownerSet?.similarityAlignmentResetReceipt ?? null,
        featureCatalog: handle.ownerSet?.featureCatalog || null
      };
      if (typeof applyEditorStateData === 'function') {
        applyEditorStateData(drawing, trustedEditorState, { normalized: true });
      } else {
        applyEditorIntentData(drawing, trustedEditorState);
        setGeneratedArtifactRef(state.featureCatalog, trustedEditorState.featureCatalog);
      }
      const presentation = mutableIntent.presentation || {};
      setGeneratedArtifactRef(state.resultGenerationKey, presentation.resultGenerationKey ?? 0);
      setGeneratedArtifactRef(state.resultPanelTab, presentation.resultPanelTab || 'preview');
      setGeneratedArtifactRef(state.errorLog, cloneJsonData(presentation.errorLog));
      if (state.labelTextScopeDialog) {
        Object.assign(
          state.labelTextScopeDialog,
          cloneJsonData(presentation.labelTextScopeDialog) || {}
        );
      }
      if (state.featureEditorStatus) {
        Object.assign(
          state.featureEditorStatus,
          cloneJsonData(presentation.featureEditorStatus) || {}
        );
      }
      setGeneratedArtifactRef(
        state.featureExtractionPending,
        Boolean(presentation.featureExtractionPending)
      );
      setGeneratedArtifactRef(
        state.featureExtractionError,
        presentation.featureExtractionError ?? null
      );

      // A Generate step restores no Session-level value (R15-23).
      const boundedUi = withoutSessionUi(ui);
      [
        'generatedLegendPosition',
        'generatedMode',
        'generatedMultiRecordCanvas',
        'generatedCircularPlotTitlePosition',
        'appliedPaletteName',
        'appliedPaletteColors',
        'pendingPaletteName',
        'pendingPaletteColors'
      ].forEach((key) => delete boundedUi[key]);
      if (typeof applyUiStateData === 'function') {
        applyUiStateData(drawing, boundedUi);
      } else {
        applyFallbackUiStateData(state, drawing, boundedUi);
      }
      await generatedArtifactRuntimeOwner?.restore?.(handle.runtimeState, { ui });
      currentGeneratedArtifactIdentity = handle.transportIdentity || null;
      currentGeneratedArtifactRetainedBytes = Number(handle.baseRetainedBytes) || 0;
      await nextTick();
      // The SVG watcher schedules its trusted mount work on a nested nextTick.
      await nextTick();
      return true;
    } finally {
      generatedArtifactRestoreDepth = Math.max(0, generatedArtifactRestoreDepth - 1);
      if (generatedArtifactRestoreDepth === 0) {
        setGeneratedArtifactRef(
          state.trustedArtifactRestoreInProgress,
          restoreTrustedStateBaseline
        );
        setGeneratedArtifactRef(
          state.semanticFileWatchersSuppressed,
          restoreSemanticSuppressionBaseline
        );
      }
    }
  };

  /**
   * @param {{ identity?: { fingerprint?: string, compactSignature?: string } } | null | undefined} before
   * @param {{ identity?: { fingerprint?: string, compactSignature?: string } } | null | undefined} after
   */
  const compareGeneratedArtifactHandles = (before, after) => Boolean(
    before?.identity?.fingerprint
    && before.identity.fingerprint === after?.identity?.fingerprint
    && before.identity.compactSignature === after?.identity?.compactSignature
  );

  // The alignment plan, its Reset receipt, and record translations of
  // drawing `drawingId`; intents and checkpoints capture and restore them
  // together. The receipt is an artifact slot value, the drawing's own
  // (`readSlotValues`, OV-1195).
  /** @param {DrawingId} drawingId */
  const captureAlignmentState = (drawingId) => ({
    receipt: readSlotValues(drawingId, ['similarityAlignmentResetReceipt']).similarityAlignmentResetReceipt,
    plan: getGeneratedArtifactRef(state.similarityAlignmentPlan, null),
    recordTranslations: getGeneratedArtifactRef(state.linearRecordTranslations, [])
  });
  /** @param {AlignmentStateData | null | undefined} alignmentState @param {DrawingId} drawingId */
  const installAlignmentState = (alignmentState, drawingId) => {
    writeArtifactValues(drawingId, { similarityAlignmentResetReceipt: cloneJsonData(alignmentState?.receipt) || null });
    setGeneratedArtifactRef(state.similarityAlignmentPlan, cloneJsonData(alignmentState?.plan) || null);
    setGeneratedArtifactRef(
      state.linearRecordTranslations,
      cloneJsonData(alignmentState?.recordTranslations) || []
    );
  };

  /** @param {D} drawing */
  const buildGeneratedArtifactSnapshot = (drawing) => {
    const ui = typeof buildUiStateData === 'function'
      ? buildUiStateData(drawing, { includePreviewNavigation: false })
      : buildFallbackUiStateData(state, drawing);
    const features = typeof buildFeatureStateData === 'function'
      ? buildFeatureStateData(drawing)
      : buildFallbackFeatureStateData(state, drawing);
    const editorState = typeof buildEditorStateData === 'function'
      ? buildEditorStateData(drawing)
      : {};
    const orthogroupState = typeof buildOrthogroupStateData === 'function'
      ? buildOrthogroupStateData(drawing)
      : buildFallbackOrthogroupStateData(state, drawing);
    const results = typeof serializeResults === 'function'
      ? serializeResults()
      : buildFallbackResultsData(state);
    const runState = typeof buildRunStateData === 'function'
      ? buildRunStateData()
      : {
          lastRunInfo: cloneJsonData(getRef(state.lastRunInfo, null)),
          pairwiseMatchFactors: clonePlainObject(getRef(state.pairwiseMatchFactors, {}))
        };

    return { ui, results, features, editorState, orthogroupState, runState };
  };

  // SE-01/N-20 (R11): a checkpoint never copies or signs the Generate-owned
  // feature catalog. It keeps the admitted catalog by reference, and History
  // retains each checkpoint object as captured.
  const checkpointFeatureCatalogs = new WeakMap();

  /**
   * @param {D} drawing
   * @param {ArtifactCheckpointData} snapshot
   */
  const applyArtifactDomains = (drawing, snapshot) => {
    const ui = snapshot?.ui || {};
    if (typeof applyResultsData === 'function') {
      applyResultsData(snapshot?.results || [], ui);
    } else {
      applyFallbackResultsData(state, snapshot?.results || []);
      applyFallbackUiStateData(state, drawing, { selectedResultIndex: ui.selectedResultIndex });
    }

    if (typeof applyFeatureStateData === 'function') {
      applyFeatureStateData(drawing, snapshot?.features || {});
    } else {
      applyFallbackFeatureStateData(state, drawing, snapshot?.features || {});
    }

    if (typeof applyOrthogroupStateData === 'function') {
      applyOrthogroupStateData(drawing, snapshot?.orthogroupState || {});
    } else {
      applyFallbackOrthogroupStateData(state, drawing, snapshot?.orthogroupState || {});
    }

    if (typeof applyEditorStateData === 'function') {
      // The checkpoint holds the editor state as captured: install a copy as
      // is, so a value such as the Result's named stroke color stays unchanged.
      applyEditorStateData(drawing, {
        ...cloneJsonData(snapshot?.editorState || {}),
        featureCatalog: checkpointFeatureCatalogs.get(snapshot) ?? null
      }, { normalized: true });
    }

    if (typeof applyRunStateData === 'function') {
      applyRunStateData(snapshot?.runState || {});
    } else {
      setRef(state.lastRunInfo, cloneJsonData(snapshot?.runState?.lastRunInfo) || null);
      setRef(state.pairwiseMatchFactors, clonePlainObject(snapshot?.runState?.pairwiseMatchFactors));
    }
  };

  // PD-OI-090: each drawing has its own History, so an intent and a checkpoint
  // hold one drawing: its settings and edits, the bindings of its slots
  // (`buildIntentFilesData`), and the shared selections; the Session-level
  // values are the router's project stream's (R15-23). Undo and Redo walk the
  // shown drawing's History and the project stream.
  /** @param {D} drawing */
  const buildDrawingConfig = (drawing) => (typeof buildConfigData === 'function'
    ? buildConfigData(drawing)
    : {
        form: drawing.form,
        adv: drawing.adv,
        linearComparisonPlan: cloneLinearComparisonPlanMetadata(drawing.linearComparisonPlan)
      });

  /** @param {DrawingId} drawingId */
  const isShown = (drawingId) => drawingId === state.activeDrawingId?.value;
  // E1: the artifacts of the drawings not shown (the slot store). History
  // writes a hidden drawing's slot values (`writeArtifactValues`); the mode
  // transition and Session Load (with its rollback) use its other named
  // operations through the service's `artifactSlots`.
  const artifactSlots = createArtifactSlotStore();
  // OV-1168, OV-1195: the artifact slot values in a drawing's intent (the
  // applied palette in `ui`, the alignment Reset receipt) are its own
  // artifact's: the displayed one's while it is shown, else the slot that
  // keeps it. A drawing with no slot has the empty slot's (null: it takes
  // the displayed palette when it is shown).
  /** @param {DrawingId} drawingId @param {readonly (keyof typeof ARTIFACT_SLOT_EMPTY)[]} keys */
  const readSlotValues = (drawingId, keys) => Object.fromEntries(keys.map((key) => [key, artifactSlots.readOfDrawing(
    drawingId, state.activeDrawingId?.value,
    () => getGeneratedArtifactRef(state[key], null), (slot) => slot.values[key]
  ) ?? null]));
  /**
   * Writes slot values of drawing `drawingId`: into the displayed artifact
   * while it is shown, else into its slot.
   * @param {DrawingId} drawingId
   * @param {Readonly<Partial<Record<keyof typeof ARTIFACT_SLOT_EMPTY, unknown>>>} values
   */
  const writeArtifactValues = (drawingId, values) => {
    if (isShown(drawingId)) {
      Object.entries(values).forEach(([key, value]) => setGeneratedArtifactRef(state[key], value));
      return;
    }
    artifactSlots.writeValues(state.drawingOf(drawingId), values);
  };
  const UI_SLOT_KEYS = /** @type {const} */ (['appliedPaletteName', 'appliedPaletteColors']);

  // Synchronous: History captures it without an await (OV-922).
  /** @param {DrawingId} drawingId */
  const buildHistoryIntent = (drawingId) => {
    const drawing = state.drawingOf(drawingId);
    const ui = typeof buildUiStateData === 'function'
      ? buildUiStateData(drawing, { includePreviewNavigation: false })
      : buildFallbackUiStateData(state, drawing);
    const uiIntent = buildUiIntentData({ ...ui, ...readSlotValues(drawingId, UI_SLOT_KEYS) }, drawing);
    uiIntent.selectedFeatureRecordIdx = getRef(state.selectedFeatureRecordIdx, 0);
    // B17, D-15-50: the composition offsets of this drawing's Results.
    const compositionDeltas = captures.composition ? cloneJsonData(captures.composition(drawingId)) : null;
    if (compositionDeltas) uiIntent.compositionUserDeltas = compositionDeltas;

    return cloneJsonData({
      config: buildDrawingConfig(drawing),
      ui: uiIntent,
      features: buildFeatureIntentData({
        featureColorOverrides: captureFeatureFills(drawing.slices.paint),
        featureVisibilityManualRules: VISIBILITY_RULES_SLICE.capture(drawing.slices.visibilityRules).featureVisibilityManualRules,
        featureOverrides: FEATURE_EDITS_SLICE.capture(drawing.slices.featureEdits).featureOverrides,
        ...labelTextEditsData(drawing.slices.labelTables)
      }),
      editorState: editorIntentOf(captureLegendPaint(drawing.slices.paint)),
      groupNames: buildGroupNameIntentData(drawing),
      files: buildIntentFilesData(state, drawing, fileStore),
      ...(drawing.mode === 'linear' ? { alignmentState: captureAlignmentState(drawingId) } : {}),
      orthogroupState: buildOrthogroupIntentData({
        selectedOrthogroupId: getRef(state.selectedOrthogroupId, ''),
        selectedOrthogroupAlignmentFeature: getRef(state.selectedOrthogroupAlignmentFeature, '')
      })
    });
  };

  const INTENT_DOMAINS = Object.freeze([
    'config', 'ui', 'features', 'editorState', 'groupNames', 'files', 'alignmentState', 'orthogroupState'
  ]);
  // The domains a step changed: the first key of each change path.
  /** @param {HistoryIntentContext} context */
  const intentDomains = (context) => {
    const changes = Array.isArray(context.changes) && context.changes.length > 0 ? context.changes : null;
    if (!changes) return new Set(INTENT_DOMAINS);
    /** @type {Set<string>} */
    const domains = new Set();
    changes.forEach((change) => {
      const path = Array.isArray(change?.path) ? change.path : [];
      if (path[0]) domains.add(String(path[0]));
    });
    return domains;
  };

  /**
   * @param {DrawingId} drawingId
   * @param {HistoryIntentData | null | undefined} intent
   * @param {HistoryIntentContext} [context]
   */
  const applyHistoryIntent = async (drawingId, intent, context = {}) => {
    if (!intent || typeof intent !== 'object') return;
    const drawing = state.drawingOf(drawingId);
    closeTransientState(state);

    const domains = intentDomains(context);
    // The alignment receipt is read against the committed Session of the
    // shown drawing's artifact before anything is applied.
    if (domains.has('alignmentState') && isShown(drawingId)) {
      const canonical = generatedArtifactRuntimeOwner?.capture?.()?.canonical?.committedCanonicalSession;
      await validateSimilarityAlignmentResetReceipt(intent.alignmentState?.receipt, canonical && {
        ...canonical,
        renderRequest: { ...canonical.renderRequest, layout: {
          ...canonical.renderRequest.layout, similarityAlignment: intent.alignmentState?.plan
        } }
      });
    }
    const retainedComparisonFiles = domains.has('config') && !domains.has('files')
      ? new Map(
          (Array.isArray(drawing.linearComparisonPlan?.edges)
            ? drawing.linearComparisonPlan.edges
            : []
          )
            .map((edge) => /** @type {[string, unknown]} */ ([String(edge?.id || ''), edge?.file ?? null]))
            .filter(([edgeId]) => edgeId)
        )
      : null;
    const suppressRef = state.semanticFileWatchersSuppressed;
    const previousSuppressed = getGeneratedArtifactRef(suppressRef, false);
    if (domains.has('files')) setGeneratedArtifactRef(suppressRef, true);
    try {
      // A drawing that is not shown takes its slot values in its slot.
      const { appliedPaletteName, appliedPaletteColors, ...drawingUi } = intent.ui || {};
      const shown = isShown(drawingId);
      const ui = shown ? intent.ui || {} : drawingUi;
      if (domains.has('ui')) {
        const inputType = inputTypeKeyOf(drawing);
        if (ui[inputType]) setRef(state[inputType], ui[inputType]);
        await nextTick();
      }

      if (domains.has('config')) {
        const config = intent.config;
        if (typeof applyConfigData === 'function' && config) {
          // The alignmentState domain below owns the plan and translations (OV-383).
          applyConfigData(drawing, config, { resolveTrackPlacements: false, applyAlignment: false });
        } else if (config?.linearComparisonPlan) {
          installComparisonPlan(drawing, config.linearComparisonPlan);
        }
      }
      if (retainedComparisonFiles) {
        (drawing.linearComparisonPlan?.edges || []).forEach((edge) => {
          const edgeId = String(edge?.id || '');
          if (edgeId && retainedComparisonFiles.has(edgeId)) {
            edge.file = retainedComparisonFiles.get(edgeId);
          }
        });
      }

      if (domains.has('ui')) {
        if (typeof applyUiStateData === 'function') {
          // A null applied palette (a drawing that had no slot) keeps the displayed one, shown or not (L1).
          applyUiStateData(drawing, ui, {
            restorePreviewNavigation: false, artifact: shown && intent.ui?.appliedPaletteName != null
          });
        } else {
          applyFallbackUiStateData(state, drawing, ui);
        }
        if (!shown && intent.ui && Object.hasOwn(intent.ui, 'appliedPaletteName')) {
          writeArtifactValues(drawingId, { appliedPaletteName, appliedPaletteColors });
        }
        if (Object.hasOwn(ui, 'selectedFeatureRecordIdx')) {
          setRef(state.selectedFeatureRecordIdx, Number.isInteger(ui.selectedFeatureRecordIdx)
            ? ui.selectedFeatureRecordIdx : 0);
        }
      }

      if (domains.has('files')) {
        applyFilesData(state, drawing, intent.files || {}, fileStore, normalizeLinearSeqList);
      }
      if (domains.has('alignmentState')) installAlignmentState(intent.alignmentState, drawingId);
      if (domains.has('features')) applyFeatureIntentData(state, drawing, intent.features || {});
      if (domains.has('editorState')) applyEditorIntentData(drawing, intent.editorState || {});
      if (domains.has('groupNames')) applyGroupNameIntentData(drawing, intent.groupNames || {});
      if (domains.has('orthogroupState')) applyOrthogroupIntentData(state, intent.orthogroupState || {});
      await nextTick();
      // The mounted Result is the shown drawing's: a drawing that is not shown
      // takes the restored intent when it is shown again (D-15-51).
      if (afterApplyHistoryIntent && shown) {
        await afterApplyHistoryIntent(intent, { ...context, domains });
      }
    } finally {
      if (domains.has('files')) setGeneratedArtifactRef(suppressRef, previousSuppressed);
    }
  };

  // A checkpoint holds the shown drawing's artifact with the drawing.
  /** @param {DrawingId} drawingId */
  const buildArtifactCheckpoint = (drawingId) => {
    const drawing = state.drawingOf(drawingId);
    const generated = buildGeneratedArtifactSnapshot(drawing);
    const { featureCatalog = null, ...editorState } = generated.editorState || {};
    const checkpoint = cloneJsonData({
      config: buildDrawingConfig(drawing),
      files: buildFilesData(state, drawing, fileStore),
      ...(drawing.mode === 'linear' ? { alignmentState: captureAlignmentState(drawingId) } : {}),
      ...generated,
      ...(generated.ui ? { ui: withoutSessionUi(generated.ui) } : {}),
      editorState
    });
    checkpointFeatureCatalogs.set(checkpoint, featureCatalog);
    return checkpoint;
  };

  /**
   * @param {DrawingId} drawingId
   * @param {ArtifactCheckpointData | null | undefined} snapshot
   */
  const applyArtifactCheckpoint = async (drawingId, snapshot) => {
    if (!snapshot || typeof snapshot !== 'object') return;
    const drawing = state.drawingOf(drawingId);
    const suppressRef = state.semanticFileWatchersSuppressed;
    const previousSuppressed = getGeneratedArtifactRef(suppressRef, false);
    setGeneratedArtifactRef(suppressRef, true);
    try {
      closeTransientState(state);

      // A3 review M1: the checkpoint restores what its drawing owns and its
      // own input type (it holds no Session-level value, R15-23).
      const ui = { ...(snapshot.ui || {}) };
      delete ui[inputTypeKeyOf(drawing) === 'lInputType' ? 'cInputType' : 'lInputType'];
      const inputType = inputTypeKeyOf(drawing);
      if (ui[inputType]) setRef(state[inputType], ui[inputType]);
      await nextTick();

      if (typeof applyConfigData === 'function' && snapshot.config) {
        applyConfigData(drawing, snapshot.config, { resolveTrackPlacements: false });
      } else if (snapshot.config?.linearComparisonPlan) {
        installComparisonPlan(drawing, snapshot.config.linearComparisonPlan);
      }

      if (typeof applyUiStateData === 'function') {
        applyUiStateData(drawing, ui, { restorePreviewNavigation: false });
      } else {
        applyFallbackUiStateData(state, drawing, ui);
      }
      await nextTick();

      applyFilesData(state, drawing, snapshot.files || {}, fileStore, normalizeLinearSeqList);
      if (Object.hasOwn(snapshot, 'alignmentState')) installAlignmentState(snapshot.alignmentState, drawingId);

      if (state.skipCaptureBaseConfig) state.skipCaptureBaseConfig.value = true;
      if (state.skipExtractOnSvgChange) state.skipExtractOnSvgChange.value = false;

      applyArtifactDomains(drawing, snapshot);

      await nextTick();
      await nextFrame();
      if (typeof applyUiStateData === 'function') {
        applyUiStateData(drawing, ui, { restorePreviewNavigation: false });
      } else {
        applyFallbackUiStateData(state, drawing, ui);
      }
    } finally {
      setGeneratedArtifactRef(suppressRef, previousSuppressed);
    }
  };

  // The History ports of one drawing (`createHistoryRouter`'s managers).
  /** @param {DrawingId} drawingId */
  const drawingHistory = (drawingId) => ({
    buildIntent: () => buildHistoryIntent(drawingId),
    /** @param {HistoryIntentData | null | undefined} intent @param {HistoryIntentContext} [context] */
    applyIntent: (intent, context) => applyHistoryIntent(drawingId, intent, context),
    buildCheckpoint: () => buildArtifactCheckpoint(drawingId),
    /** @param {ArtifactCheckpointData | null | undefined} snapshot */
    applyCheckpoint: (snapshot) => applyArtifactCheckpoint(drawingId, snapshot),
    // OV-110: the intent's composition offsets are keyed by Result identity.
    identityMapPath: COMPOSITION_INTENT_PATH,
    editedPartOf: editedIntentPartOf
  });

  // R15-23: the History ports of the Session-level values (the router's
  // project stream): the `ui` values no drawing owns and the drafts. They do
  // not depend on the drawing that `buildUiStateData` reads them with.
  const sessionHistory = () => ({
    buildIntent: () => {
      const drawing = state.drawingList()[0];
      const ui = /** @type {Record<string, unknown>} */ (typeof buildUiStateData === 'function'
        ? buildUiStateData(drawing, { includePreviewNavigation: false })
        : buildFallbackUiStateData(state, drawing));
      return cloneJsonData({
        ui: Object.fromEntries(SESSION_UI_KEYS.filter((key) => Object.hasOwn(ui, key)).map((key) => [key, ui[key]])),
        drafts: buildDraftIntentData(state)
      });
    },
    /** @param {{ ui?: SessionUiData, drafts?: DraftIntentData }} intent */
    applyIntent: async (intent) => {
      closeTransientState(state);
      if (typeof applySessionUiData === 'function') applySessionUiData(intent.ui || {});
      else applyFallbackSessionUiData(state, intent.ui || {});
      applyDraftIntentData(state, intent.drafts || {});
      await nextTick();
    }
  });

  /** @param {unknown} snapshot */
  const snapshotSignature = (snapshot) => JSON.stringify(snapshot);

  return {
    captureGeneratedArtifactHandle,
    captureGeneratedArtifactOwnerSet,
    clearGeneratedArtifactIdentity,
    collectCurrentFileIds: () => boundInputFileIds(inputBindingsOf(state), fileStore),
    compareGeneratedArtifactHandles,
    drawingHistory,
    sessionHistory,
    artifactSlots,
    captureArtifactSlot,
    installArtifactSlot,
    writeArtifactValues,
    installGeneratedArtifactOwnerSet,
    installRedrawnArtifact,
    displayedArtifactDrawing,
    registerCapture,
    restoreGeneratedArtifactHandle,
    setGeneratedArtifactIdentity,
    setGeneratedArtifactRuntimeOwner,
    setAfterApplyHistoryIntent,
    snapshotSignature
  };
};
