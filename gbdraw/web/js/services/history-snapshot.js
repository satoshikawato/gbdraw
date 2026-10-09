// @ts-check
/** @import { ArtifactSlot } from './artifact-slot.js' */
import { validateSimilarityAlignmentResetReceipt } from './session-active-config-contract.js';
import { normalizeFeatureVisibilityRule } from './feature-visibility.js';
import { serializeCleanSvg } from './svg-serialization.js';
import { cloneJsonData } from './json-clone.js';
import { replaceLayoutPreferences } from './layout-preferences.js';
import {
  isResourceBackedCanonicalComparison,
  mapResourceBackedCanonicalComparison
} from './canonical-comparisons.js';
import { recordStructuralMetric } from './runtime-test-hooks.js';
import { ARTIFACT_SLOT_EMPTY, ARTIFACT_SLOT_KEYS, createArtifactSlot } from './artifact-slot.js';

export { cloneJsonData };

const clonePlainObject = (value) => {
  if (!value || typeof value !== 'object' || Array.isArray(value)) return {};
  return cloneJsonData(value);
};

const replacePlainObject = (target, source) => {
  if (!target || typeof target !== 'object') return;
  Object.keys(target).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = value;
  });
};

const replaceRefArray = (target, source) => {
  if (!target || typeof target !== 'object' || !('value' in target)) return;
  target.value = Array.isArray(source) ? cloneJsonData(source) : [];
};

const buildFeatureIntentData = (features = {}) => ({
  featureColorOverrides: clonePlainObject(features.featureColorOverrides),
  featureVisibilityManualRules: cloneFeatureVisibilityRules(features.featureVisibilityManualRules),
  featureOverrides: clonePlainObject(features.featureOverrides),
  labelTextBulkOverrides: clonePlainObject(features.labelTextBulkOverrides)
});

const buildEditorIntentData = (editorState = {}) => ({
  legend: {
    entries: cloneJsonData(editorState?.legend?.entries) || [],
    ...(editorState?.legend?.entryOwners ? { entryOwners: cloneJsonData(editorState.legend.entryOwners) } : {}),
    deletedEntries: cloneJsonData(editorState?.legend?.deletedEntries) || [],
    dormantEntries: cloneJsonData(editorState?.legend?.dormantEntries) || [],
    colorOverrides: clonePlainObject(editorState?.legend?.colorOverrides),
    strokeOverrides: clonePlainObject(editorState?.legend?.strokeOverrides),
    addedCaptions: cloneJsonData(editorState?.legend?.addedCaptions) || []
  },
  featureStrokes: {
    overrides: clonePlainObject(editorState?.featureStrokes?.overrides)
  }
});

const buildOrthogroupIntentData = (orthogroupState = {}) => ({
  selectedOrthogroupId: String(orthogroupState.selectedOrthogroupId || ''),
  selectedOrthogroupAlignmentFeature: String(orthogroupState.selectedOrthogroupAlignmentFeature || '')
});

// A drawing's group names (its Session 46 `config.webEdits`).
/** @param {HistorySnapshotDrawing} drawing */
const buildGroupNameIntentData = (drawing) => ({
  orthogroupNameOverrides: clonePlainObject(drawing.orthogroupNameOverrides),
  orthogroupDescriptionOverrides: clonePlainObject(drawing.orthogroupDescriptionOverrides),
  orthogroupDormantOverrides: clonePlainObject(drawing.orthogroupDormantOverrides)
});

// Neither the selected Result nor the offsets read from its mounted root is
// an edit, so neither is intent; composition offsets are recorded per Result
// in `compositionUserDeltas` (B17).
const buildUiIntentData = (ui = {}) => {
  const intent = clonePlainObject(ui);
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
 * @param {HistorySnapshotDrawing} drawing
 */
const applyFeatureIntentData = (state, drawing, features = {}) => {
  if (Object.hasOwn(features, 'selectedFeatureRecordIdx')) {
    setRef(
      state.selectedFeatureRecordIdx,
      Number.isInteger(features.selectedFeatureRecordIdx) ? features.selectedFeatureRecordIdx : 0
    );
  }
  replacePlainObject(drawing.featureColorOverrides, clonePlainObject(features.featureColorOverrides));
  replaceFeatureEditState(drawing, features);
};

/** @param {HistorySnapshotDrawing} drawing */
const applyEditorIntentData = (drawing, editorState = {}) => {
  const legend = editorState.legend || {};
  replaceRefArray(drawing.legendEntries, legend.entries);
  replaceRefArray(drawing.deletedLegendEntries, legend.deletedEntries);
  replaceRefArray(drawing.dormantLegendEntries, legend.dormantEntries);
  replacePlainObject(drawing.legendColorOverrides, clonePlainObject(legend.colorOverrides));
  replacePlainObject(drawing.legendStrokeOverrides, clonePlainObject(legend.strokeOverrides));
  if (drawing.addedLegendCaptions?.value !== undefined) {
    drawing.addedLegendCaptions.value = new Set(
      Array.isArray(legend.addedCaptions) ? legend.addedCaptions.map((entry) => String(entry || '')) : []
    );
  }
  replacePlainObject(
    drawing.featureStrokeOverrides,
    clonePlainObject(editorState?.featureStrokes?.overrides)
  );
};

/** @param {Record<string, any>} state */
const applyOrthogroupIntentData = (state, orthogroupState = {}) => {
  setRef(state.selectedOrthogroupId, String(orthogroupState.selectedOrthogroupId || ''));
  setRef(
    state.selectedOrthogroupAlignmentFeature,
    String(orthogroupState.selectedOrthogroupAlignmentFeature || '')
  );
};

/** @param {HistorySnapshotDrawing} drawing */
const applyGroupNameIntentData = (drawing, names = {}) => {
  replacePlainObject(drawing.orthogroupNameOverrides, clonePlainObject(names.orthogroupNameOverrides));
  replacePlainObject(drawing.orthogroupDescriptionOverrides, clonePlainObject(names.orthogroupDescriptionOverrides));
  replacePlainObject(drawing.orthogroupDormantOverrides, clonePlainObject(names.orthogroupDormantOverrides));
};

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

const replaceLinearComparisonPlan = (target, plan = {}) => {
  if (!target || typeof target !== 'object') return;
  const cloned = cloneLinearComparisonPlanMetadata(plan);
  target.mode = cloned.mode;
  target.defaultSource = cloned.defaultSource;
  if (!Array.isArray(target.edges)) target.edges = [];
  target.edges.splice(0, target.edges.length, ...cloned.edges.map((edge) => ({
    ...edge,
    file: null
  })));
};

const cloneFeatureVisibilityRules = (rules) => (
  Array.isArray(rules) ? rules.map((rule) => normalizeFeatureVisibilityRule(rule)) : []
);

// Manual visibility rules, identity-keyed per-feature edits, and bulk label
// edits: the feature edit intent History restores as captured (R11).
/** @param {HistorySnapshotDrawing} drawing */
const replaceFeatureEditState = (drawing, features = {}) => {
  if (Array.isArray(drawing.featureVisibilityManualRules)) {
    drawing.featureVisibilityManualRules.splice(
      0,
      drawing.featureVisibilityManualRules.length,
      ...cloneFeatureVisibilityRules(features.featureVisibilityManualRules)
    );
  }
  replacePlainObject(drawing.featureOverrides, clonePlainObject(features.featureOverrides));
  replacePlainObject(drawing.labelTextBulkOverrides, clonePlainObject(features.labelTextBulkOverrides));
};

const setRef = (target, value) => {
  if (target && typeof target === 'object' && 'value' in target) {
    target.value = value;
  }
};

/**
 * @param target A Vue ref, or any other value (the fallback is returned for a non-ref).
 * @param {unknown} [fallback] The value for a target that is not a ref.
 */
const getRef = (target, fallback = null) => (
  target && typeof target === 'object' && 'value' in target
    ? target.value
    : fallback
);

const setGeneratedArtifactRef = (target, value) => {
  if (target && typeof target === 'object' && 'value' in target) target.value = value;
};

/**
 * @param target A Vue ref, or any other value (the fallback is returned for a non-ref).
 * @param {unknown} [fallback] The value for a target that is not a ref.
 */
const getGeneratedArtifactRef = (target, fallback = null) => (
  target && typeof target === 'object' && 'value' in target ? target.value : fallback
);

const artifactOwnedValue = (value) => (
  typeof globalThis.Vue?.toRaw === 'function' ? globalThis.Vue.toRaw(value) : value
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

/** @param {HistorySnapshotDrawing} drawing */
const buildFileLegendCaptionData = (drawing) => Array.from(getRef(drawing.fileLegendCaptions, new Set()) || [])
  .map((caption) => String(caption || '').trim())
  .filter(Boolean);

/** @param {HistorySnapshotDrawing} drawing */
const applyFileLegendCaptionData = (drawing, captions) => {
  if (drawing.fileLegendCaptions?.value !== undefined) {
    drawing.fileLegendCaptions.value = new Set(
      Array.isArray(captions) ? captions.map((caption) => String(caption || '').trim()).filter(Boolean) : []
    );
  }
};

/** @param {Record<string, any>} state */
const applyDraftIntentData = (state, drafts = {}) => {
  setRef(state.selectedAnnotation, cloneJsonData(drafts.selectedAnnotation) || null);
  setRef(state.selectedSpecificPreset, String(drafts.selectedSpecificPreset || ''));
  replacePlainObject(state.newSpecRule, clonePlainObject(drafts.newSpecRule));
  replacePlainObject(state.newPriorityRule, clonePlainObject(drafts.newPriorityRule));
  setRef(state.newColorFeat, String(drafts.newColorFeat || ''));
  setRef(state.newColorVal, String(drafts.newColorVal || ''));
  setRef(state.newFeatureToAdd, String(drafts.newFeatureToAdd || ''));
};

// The `ui` values that belong to a drawing (Session 46 `modes.<m>.ui`, and
// its LOSAT program); a History intent and checkpoint keep them per mode.
const DRAWING_UI_KEYS = Object.freeze([
  'canvasPadding', 'layoutPreferences', 'linearTypographyLinked', 'pendingPaletteName', 'pendingPaletteColors',
  'losatProgram'
]);
/** @param {Record<string, any>} ui */
const drawingUiPart = (ui = {}) => Object.fromEntries(
  DRAWING_UI_KEYS.filter((key) => Object.hasOwn(ui, key)).map((key) => [key, ui[key]])
);
/** @param {Record<string, any>} ui */
const appUiPart = (ui = {}) => Object.fromEntries(
  Object.entries(ui).filter(([key]) => !DRAWING_UI_KEYS.includes(key))
);
// Installs a drawing's `ui` values as captured. The shown drawing's values go
// through the UI restore with the app-level values instead.
/** @param {HistorySnapshotDrawing} drawing */
const applyDrawingUiIntentData = (drawing, ui = {}) => {
  if (ui.layoutPreferences) replaceLayoutPreferences(drawing.layoutPreferences, ui.layoutPreferences);
  if (drawing.canvasPadding && ui.canvasPadding) {
    ['top', 'right', 'bottom', 'left'].forEach((side) => {
      drawing.canvasPadding[side] = Number(ui.canvasPadding[side]) || 0;
    });
  }
  if (typeof ui.linearTypographyLinked === 'boolean') setRef(drawing.linearTypographyLinked, ui.linearTypographyLinked);
  if (ui.pendingPaletteName !== undefined) setRef(drawing.pendingPaletteName, String(ui.pendingPaletteName || ''));
  if (ui.pendingPaletteColors) setRef(drawing.pendingPaletteColors, clonePlainObject(ui.pendingPaletteColors));
  if (ui.losatProgram) setRef(drawing.losatProgram, ui.losatProgram);
};

// A History intent, checkpoint, or generated artifact holds the settings and
// edits of one diagram mode: the drawing of that mode.
/**
 * @param {Record<string, any>} state
 * @param {unknown} mode
 * @returns {HistorySnapshotDrawing}
 */
const drawingOfMode = (state, mode) => state.drawings[mode === 'linear' ? 'linear' : 'circular'];

const nextFrame = () => /** @type {Promise<void>} */ (new Promise((resolve) => {
  if (typeof window !== 'undefined' && typeof window.requestAnimationFrame === 'function') {
    window.requestAnimationFrame(() => resolve());
  } else {
    setTimeout(resolve, 0);
  }
}));

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
 * @param {HistorySnapshotDrawing} drawing
 */
const buildFallbackUiStateData = (state, drawing) => ({
  title: getRef(state.sessionTitle, ''),
  mode: getRef(state.mode, 'circular'),
  cInputType: getRef(state.cInputType, 'gb'),
  lInputType: getRef(state.lInputType, 'gb'),
  losatProgram: getRef(drawing.losatProgram, 'blastn'),
  selectedResultIndex: getRef(state.selectedResultIndex, 0),
  downloadDpi: getRef(state.downloadDpi, 300),
  canvasPadding: { ...(drawing.canvasPadding || {}) },
  generatedLegendPosition: getRef(state.generatedLegendPosition, 'left'),
  generatedMode: getRef(state.generatedMode, 'circular'),
  generatedMultiRecordCanvas: Boolean(getRef(state.generatedMultiRecordCanvas, false)),
  generatedCircularPlotTitlePosition: getRef(state.generatedCircularPlotTitlePosition, 'none'),
  layoutPreferences: clonePlainObject(drawing.layoutPreferences),
  autoLabelReflow: Boolean(getRef(state.autoLabelReflowEnabled, false)),
  paletteInstantPreviewEnabled: Boolean(getRef(state.paletteInstantPreviewEnabled, false)),
  appliedPaletteName: getRef(state.appliedPaletteName, 'default'),
  appliedPaletteColors: clonePlainObject(getRef(state.appliedPaletteColors, {})),
  pendingPaletteName: getRef(drawing.pendingPaletteName, ''),
  pendingPaletteColors: clonePlainObject(getRef(drawing.pendingPaletteColors, {})),
  legendCurrentOffset: { ...(state.legendCurrentOffset || {}) },
  diagramOffset: { ...(state.diagramOffset || {}) },
  lengthBarUserOffset: { ...(state.lengthBarUserOffset || {}) },
  plotTitleUserOffset: { ...(state.plotTitleUserOffset || {}) }
});

// `mode` is restored by the mode transition (`restoreMode`), not here.
/**
 * @param {Record<string, any>} state
 * @param {HistorySnapshotDrawing} drawing
 */
const applyFallbackUiStateData = (state, drawing, ui = {}) => {
  if (typeof ui.title === 'string') setRef(state.sessionTitle, ui.title);
  if (ui.cInputType) setRef(state.cInputType, ui.cInputType);
  if (ui.lInputType) setRef(state.lInputType, ui.lInputType);
  if (ui.losatProgram) setRef(drawing.losatProgram, ui.losatProgram);
  if (ui.downloadDpi) setRef(state.downloadDpi, ui.downloadDpi);
  if (ui.generatedLegendPosition) setRef(state.generatedLegendPosition, ui.generatedLegendPosition);
  if (ui.generatedMode) setRef(state.generatedMode, ui.generatedMode);
  if (Object.prototype.hasOwnProperty.call(ui, 'generatedMultiRecordCanvas')) {
    setRef(state.generatedMultiRecordCanvas, Boolean(ui.generatedMultiRecordCanvas));
  }
  if (ui.generatedCircularPlotTitlePosition) {
    setRef(state.generatedCircularPlotTitlePosition, ui.generatedCircularPlotTitlePosition);
  }
  replaceLayoutPreferences(drawing.layoutPreferences, ui.layoutPreferences);
  if (drawing.canvasPadding && ui.canvasPadding) {
    drawing.canvasPadding.top = Number(ui.canvasPadding.top) || 0;
    drawing.canvasPadding.right = Number(ui.canvasPadding.right) || 0;
    drawing.canvasPadding.bottom = Number(ui.canvasPadding.bottom) || 0;
    drawing.canvasPadding.left = Number(ui.canvasPadding.left) || 0;
  }
  setRef(state.autoLabelReflowEnabled, Boolean(ui.autoLabelReflow));
  setRef(state.paletteInstantPreviewEnabled, Boolean(ui.paletteInstantPreviewEnabled));
  if (ui.appliedPaletteName !== undefined) setRef(state.appliedPaletteName, String(ui.appliedPaletteName || 'default'));
  if (ui.appliedPaletteColors) setRef(state.appliedPaletteColors, clonePlainObject(ui.appliedPaletteColors));
  if (ui.pendingPaletteName !== undefined) setRef(drawing.pendingPaletteName, String(ui.pendingPaletteName || ''));
  if (ui.pendingPaletteColors) setRef(drawing.pendingPaletteColors, clonePlainObject(ui.pendingPaletteColors));
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
    setRef(state.selectedResultIndex, count > 0 ? Math.max(0, Math.min(ui.selectedResultIndex, count - 1)) : 0);
  }
};

/**
 * @param {Record<string, any>} state
 * @param {HistorySnapshotDrawing} drawing
 */
const buildFallbackFeatureStateData = (state, drawing) => ({
  extractedFeatures: cloneJsonData(getRef(state.extractedFeatures, [])) || [],
  featureRecordIds: cloneJsonData(getRef(state.featureRecordIds, [])) || [],
  selectedFeatureRecordIdx: getRef(state.selectedFeatureRecordIdx, 0),
  featureColorOverrides: clonePlainObject(drawing.featureColorOverrides),
  featureVisibilityManualRules: cloneFeatureVisibilityRules(drawing.featureVisibilityManualRules),
  featureOverrides: clonePlainObject(drawing.featureOverrides),
  labelTextBulkOverrides: clonePlainObject(drawing.labelTextBulkOverrides)
});

/**
 * @param {Record<string, any>} state
 * @param {HistorySnapshotDrawing} drawing
 */
const applyFallbackFeatureStateData = (state, drawing, features = {}) => {
  setRef(state.extractedFeatures, cloneJsonData(features.extractedFeatures) || []);
  setRef(state.featureRecordIds, cloneJsonData(features.featureRecordIds) || []);
  setRef(
    state.selectedFeatureRecordIdx,
    Number.isInteger(features.selectedFeatureRecordIdx) ? features.selectedFeatureRecordIdx : 0
  );
  replacePlainObject(drawing.featureColorOverrides, clonePlainObject(features.featureColorOverrides));
  replaceFeatureEditState(drawing, features);
};

/**
 * @param {Record<string, any>} state
 * @param {HistorySnapshotDrawing} drawing
 */
const buildFallbackOrthogroupStateData = (state, drawing) => ({
  groups: cloneJsonData(getRef(state.orthogroups, [])) || [],
  selectedOrthogroupId: getRef(state.selectedOrthogroupId, ''),
  selectedOrthogroupAlignmentFeature: getRef(state.selectedOrthogroupAlignmentFeature, ''),
  orthogroupNameOverrides: clonePlainObject(drawing.orthogroupNameOverrides),
  orthogroupDescriptionOverrides: clonePlainObject(drawing.orthogroupDescriptionOverrides),
  orthogroupDormantOverrides: clonePlainObject(drawing.orthogroupDormantOverrides)
});

/**
 * @param {Record<string, any>} state
 * @param {HistorySnapshotDrawing} drawing
 */
const applyFallbackOrthogroupStateData = (state, drawing, data = {}) => {
  setRef(state.orthogroups, cloneJsonData(data.groups) || []);
  setRef(state.selectedOrthogroupId, String(data.selectedOrthogroupId || ''));
  setRef(state.selectedOrthogroupAlignmentFeature, String(data.selectedOrthogroupAlignmentFeature || ''));
  replacePlainObject(drawing.orthogroupNameOverrides, clonePlainObject(data.orthogroupNameOverrides));
  replacePlainObject(drawing.orthogroupDescriptionOverrides, clonePlainObject(data.orthogroupDescriptionOverrides));
  replacePlainObject(drawing.orthogroupDormantOverrides, clonePlainObject(data.orthogroupDormantOverrides));
};

const buildFallbackResultsData = (state) => {
  const currentSvg = (() => {
    const svg = state.svgContainer?.value?.querySelector?.('svg');
    if (!svg || typeof XMLSerializer === 'undefined') return null;
    return serializeCleanSvg(svg);
  })();
  const selected = getRef(state.selectedResultIndex, 0);
  return (getRef(state.results, []) || []).map((result, index) => ({
    name: result?.name || `Result ${index + 1}`,
    content: index === selected && currentSvg ? currentSvg : String(result?.content || '')
  }));
};

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

// The History intent holds every file binding by reference, including the
// BLAST rows and comparison sequences of a LOSAT-cache replay: its ring rows and
// managed track slots name those rows, so Undo restores them together (B23).
/**
 * @param {Record<string, any>} state
 * @param {HistorySnapshotDrawing} drawing
 */
const buildIntentFilesData = (state, drawing, fileStore) => ({
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
  d_color: fileStore.describeValue(state.files?.d_color),
  t_color: fileStore.describeValue(state.files?.t_color),
  blacklist: fileStore.describeValue(state.files?.blacklist),
  whitelist: fileStore.describeValue(state.files?.whitelist),
  qualifier_priority: fileStore.describeValue(state.files?.qualifier_priority),
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
});

// An artifact checkpoint adds the generated Linear comparisons.
/**
 * @param {Record<string, any>} state
 * @param {HistorySnapshotDrawing} drawing
 */
const buildFilesData = (state, drawing, fileStore) => ({
  ...buildIntentFilesData(state, drawing, fileStore),
  linearCanonicalComparisons: (
    Array.isArray(state.files?.linearCanonicalComparisons)
      ? state.files.linearCanonicalComparisons
      : []
  ).map((comparison) => (
    isResourceBackedCanonicalComparison(comparison)
      ? {
          ...mapResourceBackedCanonicalComparison(comparison),
          file: fileStore.describeValue(comparison.file)
        }
      : cloneJsonData(comparison)
  ))
});

const collectCurrentFileIds = (state, fileStore) => {
  const fileIds = new Set();
  const register = (value) => {
    if (Array.isArray(value)) {
      value.forEach(register);
      return;
    }
    const fileId = fileStore.registerFile(value);
    if (fileId) fileIds.add(fileId);
  };
  const files = state.files || {};
  [
    files.c_gb,
    files.c_gff,
    files.c_fasta,
    files.c_depth,
    files.c_conservation_blasts,
    files.c_conservation_fastas,
    files.c_conservation_sequence_sources,
    files.d_color,
    files.t_color,
    files.blacklist,
    files.whitelist,
    files.qualifier_priority
  ].forEach(register);
  (Array.isArray(files.linearCanonicalComparisons)
    ? files.linearCanonicalComparisons
    : []
  ).forEach((comparison) => register(comparison?.file));
  Array.from(state.linearSeqs || []).forEach((sequence) => {
    register(sequence?.gb);
    register(sequence?.gff);
    register(sequence?.fasta);
    register(sequence?.depth);
  });
  // Every drawing's comparison plan holds its uploaded comparison files.
  new Set(Object.values(state.drawings)).forEach((drawing) => {
    Array.from(drawing.linearComparisonPlan?.edges || []).forEach((comparison) => {
      register(comparison?.file);
    });
  });
  return fileIds;
};

/**
 * @param {Record<string, any>} state
 * @param {HistorySnapshotDrawing} drawing
 * @param {HistorySnapshotFileStore} fileStore
 * @param {HistorySnapshotServiceOptions['normalizeLinearSeqList']} [normalizeLinearSeqList]
 */
const applyFilesData = (state, drawing, filesData, fileStore, normalizeLinearSeqList = null) => {
  if (!state.files) return;
  state.matchSequenceRegistry?.reset?.();
  const restore = (value) => fileStore.restoreValue(value);
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
  if (Object.prototype.hasOwnProperty.call(filesData || {}, 'linearCanonicalComparisons')) {
    state.files.linearCanonicalComparisons = Array.isArray(filesData?.linearCanonicalComparisons)
      ? filesData.linearCanonicalComparisons.map((comparison) => (
          isResourceBackedCanonicalComparison(comparison)
            ? mapResourceBackedCanonicalComparison(comparison, restore)
            : cloneJsonData(comparison)
        ))
      : [];
  }
  state.files.d_color = restore(filesData?.d_color);
  state.files.t_color = restore(filesData?.t_color);
  state.files.blacklist = restore(filesData?.blacklist);
  state.files.whitelist = restore(filesData?.whitelist);
  state.files.qualifier_priority = restore(filesData?.qualifier_priority);

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
 *   The three functions of the History file store that the snapshot service
 *   calls (F-05). `services/history.js` calls two other functions of the same
 *   store, so the receivers share the object and no function.
 * @property {(value: any) => any} describeValue
 *   A file, or an array of files, as the descriptors an intent holds (`null` for a non-file).
 * @property {(file: any) => string | null} registerFile
 *   Stores a file by reference and returns its id; `null` for a non-file.
 * @property {(value: any) => any} restoreValue
 *   The stored file, or an array of them, that the descriptors name (`null` when none).
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
 * @typedef {any} HistorySnapshotDrawing
 *   The drawing (`DrawingState` of state.js, a higher layer) whose settings and edits
 *   a History intent, checkpoint, or generated artifact holds.
 *
 * @typedef {Record<string, any>} HistorySnapshotData
 *   A domain payload that the snapshot service holds by value (config, ui, features, editor state,
 *   orthogroup state, run state): the owner that builds it declares its shape.
 *
 * @typedef {object} HistorySnapshotServiceOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {HistorySnapshotFileStore} fileStore
 *   The file store that holds the files an intent and a checkpoint name.
 * @property {() => Promise<unknown>} [nextTick] Vue `nextTick`, awaited between restore steps.
 * @property {((rows: HistorySnapshotData[]) => HistorySnapshotData[]) | null} [normalizeLinearSeqList]
 *   The Linear record rows' normalization after a files restore.
 * @property {((drawing: HistorySnapshotDrawing) => HistorySnapshotData) | null} [buildConfigData] The Settings capture.
 * @property {((drawing: HistorySnapshotDrawing, config: HistorySnapshotData, options: { resolveTrackPlacements: boolean }) => unknown) | null} [applyConfigData]
 *   The Settings restore.
 * @property {((drawing: HistorySnapshotDrawing, options: { includePreviewNavigation: boolean }) => HistorySnapshotData) | null} [buildUiStateData]
 *   The UI capture; `includePreviewNavigation` adds pan, zoom, and tab state.
 * @property {((drawing: HistorySnapshotDrawing, ui: HistorySnapshotData, options?: { restorePreviewNavigation: boolean }) => unknown) | null} [applyUiStateData]
 *   The UI restore.
 * @property {((drawing: HistorySnapshotDrawing) => HistorySnapshotData) | null} [buildFeatureStateData] The feature-override capture.
 * @property {((drawing: HistorySnapshotDrawing, features: HistorySnapshotData) => unknown) | null} [applyFeatureStateData]
 *   The feature-override restore.
 * @property {((drawing: HistorySnapshotDrawing) => HistorySnapshotData) | null} [buildEditorStateData] The legend and stroke editor capture.
 * @property {((drawing: HistorySnapshotDrawing, editorState: HistorySnapshotData, options: { normalized: boolean }) => unknown) | null} [applyEditorStateData]
 *   The editor restore; `normalized` marks a state that is installed as captured.
 * @property {((drawing: HistorySnapshotDrawing) => HistorySnapshotData) | null} [buildOrthogroupStateData] The orthogroup capture.
 * @property {((drawing: HistorySnapshotDrawing, orthogroupState: HistorySnapshotData) => unknown) | null} [applyOrthogroupStateData]
 *   The orthogroup restore.
 * @property {(() => HistorySnapshotData[]) | null} [serializeResults] The generated Results' capture.
 * @property {((results: HistorySnapshotData[], ui: HistorySnapshotData) => unknown) | null} [applyResultsData]
 *   The generated Results' restore, with the checkpoint's UI state.
 * @property {(() => HistorySnapshotData) | null} [buildRunStateData] The last-run capture.
 * @property {((runState: HistorySnapshotData) => unknown) | null} [applyRunStateData] The last-run restore.
 */

/** @param {HistorySnapshotServiceOptions} options */
export const createHistorySnapshotService = ({
  state,
  fileStore,
  nextTick = async () => {},
  normalizeLinearSeqList = null,
  buildConfigData = null,
  applyConfigData = null,
  buildUiStateData = null,
  applyUiStateData = null,
  buildFeatureStateData = null,
  applyFeatureStateData = null,
  buildEditorStateData = null,
  applyEditorStateData = null,
  buildOrthogroupStateData = null,
  applyOrthogroupStateData = null,
  serializeResults = null,
  applyResultsData = null,
  buildRunStateData = null,
  applyRunStateData = null
}) => {
  if (!state || !fileStore) {
    throw new Error('createHistorySnapshotService requires state and fileStore.');
  }

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

  const setAfterApplyHistoryIntent = (callback) => {
    afterApplyHistoryIntent = typeof callback === 'function' ? callback : null;
  };

  // Intent captures of owners created after this service (R13): the
  // composition root registers each once its owner exists. An unregistered
  // capture adds nothing to the intent, as for a service without that owner.
  /** @type {{ legend: (() => unknown) | null, composition: (() => unknown) | null }} */
  const captures = { legend: null, composition: null };
  const registerCapture = (name, capture) => {
    if (!Object.hasOwn(captures, name)) throw new Error(`Unknown History intent capture: ${name}`);
    captures[name] = typeof capture === 'function' ? capture : null;
  };

  // E1: the composition root's mode transition and the committed Session of
  // each mode's artifact, registered once they exist (R13). A History restore
  // of a mode switch runs the transition, so Undo and Redo swap the modes'
  // artifact slots. `restoreMode` is the one mode writer of the History
  // restores; it writes `mode` itself only for a service without a
  // composition root (unit tests).
  /** @type {((mode: 'circular' | 'linear') => unknown) | null} */
  let modeTransition = null;
  /** @type {((mode: 'circular' | 'linear') => Record<string, any> | null) | null} */
  let modeCommittedSession = null;
  /** @param {((mode: 'circular' | 'linear') => unknown) | null} transition */
  const registerModeTransition = (transition) => {
    modeTransition = typeof transition === 'function' ? transition : null;
  };
  /** @param {((mode: 'circular' | 'linear') => Record<string, any> | null) | null} reader */
  const registerModeCommittedSession = (reader) => {
    modeCommittedSession = typeof reader === 'function' ? reader : null;
  };
  /** @param {string} mode */
  const restoreMode = async (mode) => {
    const next = mode === 'linear' ? 'linear' : 'circular';
    if (getGeneratedArtifactRef(state.mode, 'circular') === next) return;
    if (modeTransition) modeTransition(next);
    else setGeneratedArtifactRef(state.mode, next);
    await nextTick();
  };
  // An artifact handle or checkpoint restores the displayed mode's artifact.
  // History is last in, first out, so it finds the mode it was captured in; a
  // mismatch is counted, and the transition shows that mode first.
  /** @param {string} mode */
  const restoreArtifactMode = async (mode) => {
    if (getGeneratedArtifactRef(state.mode, 'circular') === (mode === 'linear' ? 'linear' : 'circular')) return;
    recordStructuralMetric('artifactRestoreModeMismatchCount', 1);
    await restoreMode(mode);
  };

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

  const captureGeneratedArtifactOwnerSet = () => {
    const drawing = drawingOfMode(state, getGeneratedArtifactRef(state.generatedMode, null));
    const results = artifactOwnedValue(getGeneratedArtifactRef(state.results, []));
    return Object.freeze({
      results,
      featureCatalog: artifactOwnedValue(getGeneratedArtifactRef(state.featureCatalog, null)),
      extractedFeatures: artifactOwnedValue(getGeneratedArtifactRef(state.extractedFeatures, null)),
      biologicalFeatures: artifactOwnedValue(getGeneratedArtifactRef(state.biologicalFeatures, null)),
      featureRecordIds: artifactOwnedValue(getGeneratedArtifactRef(state.featureRecordIds, null)),
      orthogroups: artifactOwnedValue(getGeneratedArtifactRef(state.orthogroups, null)),
      featureOrthogroupIndex: artifactOwnedValue(
        getGeneratedArtifactRef(state.featureOrthogroupIndex, null)
      ),
      collinearGroups: artifactOwnedValue(getGeneratedArtifactRef(state.collinearGroups, null)),
      similarityAlignmentResetReceipt: artifactOwnedValue(
        getGeneratedArtifactRef(state.similarityAlignmentResetReceipt, null)
      ),
      similarityAlignmentPlan: artifactOwnedValue(
        getGeneratedArtifactRef(state.similarityAlignmentPlan, null)
      ),
      linearRecordTranslations: artifactOwnedValue(
        getGeneratedArtifactRef(state.linearRecordTranslations, null)
      ),
      linearRecordOrientations: artifactOwnedValue(captureLinearRecordOrientations()),
      annotationWarnings: artifactOwnedValue(getGeneratedArtifactRef(state.annotationWarnings, null)),
      featureIdentityNotices: artifactOwnedValue(
        getGeneratedArtifactRef(state.featureIdentityNotices, null)
      ),
      comparisonWarnings: artifactOwnedValue(getGeneratedArtifactRef(state.comparisonWarnings, null)),
      specificRules: (drawing.manualSpecificRules || []).map(rule => ({ ...rule })),
      fileLegendCaptions: new Set(drawing.fileLegendCaptions?.value || []),
      trackSlotResolvedGeometry: artifactOwnedValue(
        getGeneratedArtifactRef(state.trackSlotResolvedGeometry, null)
      ),
      proteinIdentityManifest: artifactOwnedValue(
        getGeneratedArtifactRef(state.proteinIdentityManifest, null)
      ),
      legacyProteinRawCandidates: artifactOwnedValue(
        getGeneratedArtifactRef(state.legacyProteinRawCandidates, null)
      ),
      legacyProteinDerivedEvidence: artifactOwnedValue(
        getGeneratedArtifactRef(state.legacyProteinDerivedEvidence, null)
      ),
      losatCache: artifactOwnedValue(getGeneratedArtifactRef(state.losatCache, null)),
      losatDerivedCache: artifactOwnedValue(
        getGeneratedArtifactRef(state.losatDerivedCache, null)
      ),
      losatCacheInfo: artifactOwnedValue(getGeneratedArtifactRef(state.losatCacheInfo, null)),
      matchSequenceOwner: state.matchSequenceRegistry?.captureTrustedOwner?.() || null,
      lastRunInfo: artifactOwnedValue(getGeneratedArtifactRef(state.lastRunInfo, null)),
      pairwiseMatchFactors: artifactOwnedValue(
        getGeneratedArtifactRef(state.pairwiseMatchFactors, null)
      ),
      editableLabels: artifactOwnedValue(getGeneratedArtifactRef(state.editableLabels, null)),
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
      appliedPaletteColors: artifactOwnedValue(
        getGeneratedArtifactRef(state.appliedPaletteColors, null)
      ),
      pendingPaletteName: getGeneratedArtifactRef(drawing.pendingPaletteName, ''),
      pendingPaletteColors: artifactOwnedValue(
        getGeneratedArtifactRef(drawing.pendingPaletteColors, null)
      )
    });
  };

  const installGeneratedArtifactOwnerSet = (
    ownerSet,
    /** @type {{ selectedResultIndex?: number, installResults?: ((results: any[]) => unknown) | null }} */
    { selectedResultIndex = 0, installResults = null } = {}
  ) => {
    if (!ownerSet || typeof ownerSet !== 'object') {
      throw new Error('A generated artifact owner set is required.');
    }
    const drawing = drawingOfMode(state, ownerSet.generatedMode);
    const results = Array.isArray(ownerSet.results) ? ownerSet.results : [];
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
    if (drawing.manualSpecificRules && ownerSet.specificRules) {
      drawing.manualSpecificRules.splice(0, drawing.manualSpecificRules.length, ...ownerSet.specificRules.map(rule => ({ ...rule })));
    }
    if (drawing.fileLegendCaptions && ownerSet.fileLegendCaptions) drawing.fileLegendCaptions.value = new Set(ownerSet.fileLegendCaptions);
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
    setGeneratedArtifactRef(drawing.pendingPaletteName, ownerSet.pendingPaletteName || '');
    setGeneratedArtifactRef(drawing.pendingPaletteColors, ownerSet.pendingPaletteColors || {});
  };

  // E1: the displayed mode's generated artifact, by reference, key by key.
  // The state's `generatedMode` names the slot's mode.
  const captureArtifactSlot = () => createArtifactSlot({
    mode: getGeneratedArtifactRef(state.generatedMode, null),
    values: Object.fromEntries(ARTIFACT_SLOT_KEYS.map((key) => [
      key, artifactOwnedValue(getGeneratedArtifactRef(state[key], null))
    ])),
    legendInventory: getGeneratedArtifactRef(state.originalLegendOrder, []) || [],
    matchSequenceOwner: state.matchSequenceRegistry?.captureTrustedOwner?.() || null,
    runtimeState: generatedArtifactRuntimeOwner?.capture?.() || null,
    transportIdentity: currentGeneratedArtifactIdentity,
    retainedBytes: currentGeneratedArtifactRetainedBytes
  });

  // Installs one mode's slot, or an empty slot of `mode` (no Result, no
  // committed Session), as the displayed artifact. Nothing outside the slot
  // changes, so the draft, the editor intent and the caches stay shared.
  /**
   * @param {Readonly<ArtifactSlot> | null} slot
   * @param {{ mode: 'circular' | 'linear', ui?: { canvasPan?: { x: number, y: number } } }} options
   *   `ui.canvasPan` is the pan the runtime owner's viewport reset keeps.
   */
  const installArtifactSlot = (slot, { mode, ui = {} }) => {
    const slotMode = mode === 'linear' ? 'linear' : 'circular';
    if (slot && slot.mode !== slotMode) {
      throw new Error(`A ${slot.mode} artifact cannot be installed for ${slotMode}.`);
    }
    ARTIFACT_SLOT_KEYS.forEach((key) => {
      setGeneratedArtifactRef(state[key], slot
        ? slot.values[key]
        : ARTIFACT_SLOT_EMPTY[key](slotMode, getGeneratedArtifactRef(state[key], null)));
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
        results: artifactOwnedValue(committedResults)
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

  const captureGeneratedArtifactHandle = () => {
    recordStructuralMetric('generatedArtifactFullSignatureCount', 0);
    recordStructuralMetric('generatedArtifactHeavyTraversalCount', 0);
    recordStructuralMetric('generatedArtifactMutableIntentSnapshotCount', 1);
    const ownerSet = captureGeneratedArtifactOwnerSet();
    const drawing = drawingOfMode(state, getRef(state.mode, 'circular'));
    const ui = typeof buildUiStateData === 'function'
      ? buildUiStateData(drawing, { includePreviewNavigation: true })
      : buildFallbackUiStateData(state, drawing);
    const features = {
      selectedFeatureRecordIdx: getGeneratedArtifactRef(state.selectedFeatureRecordIdx, 0),
      featureColorOverrides: clonePlainObject(drawing.featureColorOverrides),
      featureVisibilityManualRules: cloneFeatureVisibilityRules(
        drawing.featureVisibilityManualRules
      ),
      featureOverrides: clonePlainObject(drawing.featureOverrides),
      labelOverrideRows: cloneJsonData(
        getGeneratedArtifactRef(drawing.canonicalLabelOverrideRows, [])
      ) || [],
      labelTextBulkOverrides: clonePlainObject(drawing.labelTextBulkOverrides)
    };
    const editorState = {
      legend: {
        entries: cloneJsonData(getGeneratedArtifactRef(drawing.legendEntries, [])) || [],
        deletedEntries: cloneJsonData(
          getGeneratedArtifactRef(drawing.deletedLegendEntries, [])
        ) || [],
        dormantEntries: cloneJsonData(
          getGeneratedArtifactRef(drawing.dormantLegendEntries, [])
        ) || [],
        originalOrder: cloneJsonData(
          getGeneratedArtifactRef(state.originalLegendOrder, [])
        ) || [],
        originalColors: clonePlainObject(
          getGeneratedArtifactRef(state.originalLegendColors, {})
        ),
        colorOverrides: clonePlainObject(drawing.legendColorOverrides),
        strokeOverrides: clonePlainObject(drawing.legendStrokeOverrides),
        addedCaptions: Array.from(
          getGeneratedArtifactRef(drawing.addedLegendCaptions, new Set()) || []
        )
      },
      featureStrokes: { overrides: clonePlainObject(drawing.featureStrokeOverrides) },
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
      orthogroupNameOverrides: clonePlainObject(drawing.orthogroupNameOverrides),
      orthogroupDescriptionOverrides: clonePlainObject(
        drawing.orthogroupDescriptionOverrides
      ),
      orthogroupDormantOverrides: clonePlainObject(drawing.orthogroupDormantOverrides)
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
      ownerSet,
      mutableIntent,
      runtimeState,
      transportIdentity,
      baseRetainedBytes: currentGeneratedArtifactRetainedBytes,
      fileIds: Object.freeze([...collectCurrentFileIds(state, fileStore)]),
      identity: Object.freeze({
        fingerprint: transportIdentity?.fingerprint || '',
        compactSignature
      }),
      retainedBytes: currentGeneratedArtifactRetainedBytes
        + Math.max(0, Number(runtimeState?.retainedBytes) || 0)
        + compactSignature.length * 2
    });
  };

  const restoreGeneratedArtifactHandle = async (
    handle,
    { clearFailedGeneratePresentation = false } = {}
  ) => {
    if (!handle || handle.kind !== 'GeneratedArtifactHandle') return false;
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
      if (ui.mode) await restoreArtifactMode(ui.mode);
      const drawing = drawingOfMode(state, getRef(state.mode, 'circular'));
      if (ui.cInputType) setGeneratedArtifactRef(state.cInputType, ui.cInputType);
      if (ui.lInputType) setGeneratedArtifactRef(state.lInputType, ui.lInputType);
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
      replacePlainObject(
        drawing.orthogroupNameOverrides,
        clonePlainObject(mutableIntent.orthogroupState?.orthogroupNameOverrides)
      );
      replacePlainObject(
        drawing.orthogroupDescriptionOverrides,
        clonePlainObject(mutableIntent.orthogroupState?.orthogroupDescriptionOverrides)
      );
      replacePlainObject(
        drawing.orthogroupDormantOverrides,
        clonePlainObject(mutableIntent.orthogroupState?.orthogroupDormantOverrides)
      );

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

      const boundedUi = { ...ui };
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

  const compareGeneratedArtifactHandles = (before, after) => Boolean(
    before?.identity?.fingerprint
    && before.identity.fingerprint === after?.identity?.fingerprint
    && before.identity.compactSignature === after?.identity?.compactSignature
  );

  // The alignment plan, its Reset receipt, and record translations; intents
  // and checkpoints capture and restore them together.
  const captureAlignmentState = () => ({
    receipt: getGeneratedArtifactRef(state.similarityAlignmentResetReceipt, null),
    plan: getGeneratedArtifactRef(state.similarityAlignmentPlan, null),
    recordTranslations: getGeneratedArtifactRef(state.linearRecordTranslations, [])
  });
  const installAlignmentState = (alignmentState) => {
    setGeneratedArtifactRef(
      state.similarityAlignmentResetReceipt,
      cloneJsonData(alignmentState?.receipt) || null
    );
    setGeneratedArtifactRef(state.similarityAlignmentPlan, cloneJsonData(alignmentState?.plan) || null);
    setGeneratedArtifactRef(
      state.linearRecordTranslations,
      cloneJsonData(alignmentState?.recordTranslations) || []
    );
  };

  /** @param {HistorySnapshotDrawing} drawing */
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

  /** @param {HistorySnapshotDrawing} drawing */
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

  // A History intent holds both drawings (PD-OI-086): an edit changes the
  // shown mode's drawing only, while Reset Settings and Session Load change
  // both, so Undo restores each drawing as captured. The Legend entry owners
  // belong to the mounted Result, so only the shown drawing records them.
  /**
   * @param {'circular' | 'linear'} drawingMode
   * @param {boolean} shown
   */
  const buildModeIntent = (drawingMode, shown) => {
    const drawing = drawingOfMode(state, drawingMode);
    const config = typeof buildConfigData === 'function'
      ? buildConfigData(drawing)
      : {
          form: drawing.form,
          adv: drawing.adv,
          linearComparisonPlan: cloneLinearComparisonPlanMetadata(drawing.linearComparisonPlan)
        };
    const ui = typeof buildUiStateData === 'function'
      ? buildUiStateData(drawing, { includePreviewNavigation: false })
      : buildFallbackUiStateData(state, drawing);
    return {
      config,
      ui: drawingUiPart(ui),
      features: buildFeatureIntentData({
        featureColorOverrides: drawing.featureColorOverrides,
        featureVisibilityManualRules: drawing.featureVisibilityManualRules,
        featureOverrides: drawing.featureOverrides,
        labelTextBulkOverrides: drawing.labelTextBulkOverrides
      }),
      editorState: buildEditorIntentData({
        legend: {
          entries: getRef(drawing.legendEntries, []),
          ...(shown && captures.legend ? { entryOwners: captures.legend() } : {}),
          deletedEntries: getRef(drawing.deletedLegendEntries, []),
          dormantEntries: getRef(drawing.dormantLegendEntries, []),
          colorOverrides: drawing.legendColorOverrides,
          strokeOverrides: drawing.legendStrokeOverrides,
          addedCaptions: Array.from(getRef(drawing.addedLegendCaptions, new Set()) || [])
        },
        featureStrokes: { overrides: drawing.featureStrokeOverrides }
      }),
      groupNames: buildGroupNameIntentData(drawing),
      fileLegendCaptions: buildFileLegendCaptionData(drawing)
    };
  };

  const buildHistoryIntent = async () => {
    const shownMode = getRef(state.mode, 'circular') === 'linear' ? 'linear' : 'circular';
    const drawing = drawingOfMode(state, shownMode);
    const ui = typeof buildUiStateData === 'function'
      ? buildUiStateData(drawing, { includePreviewNavigation: false })
      : buildFallbackUiStateData(state, drawing);
    const uiIntent = buildUiIntentData(appUiPart(ui));
    uiIntent.selectedFeatureRecordIdx = getRef(state.selectedFeatureRecordIdx, 0);
    const compositionDeltas = captures.composition ? cloneJsonData(captures.composition()) : null;
    if (compositionDeltas) uiIntent.compositionUserDeltas = compositionDeltas;

    return cloneJsonData({
      modes: {
        circular: buildModeIntent('circular', shownMode === 'circular'),
        linear: buildModeIntent('linear', shownMode === 'linear')
      },
      files: buildIntentFilesData(state, drawingOfMode(state, 'linear'), fileStore),
      alignmentState: captureAlignmentState(),
      ui: uiIntent,
      drafts: buildDraftIntentData(state),
      orthogroupState: buildOrthogroupIntentData({
        selectedOrthogroupId: getRef(state.selectedOrthogroupId, ''),
        selectedOrthogroupAlignmentFeature: getRef(state.selectedOrthogroupAlignmentFeature, '')
      })
    });
  };

  const INTENT_DOMAINS = Object.freeze(['files', 'alignmentState', 'ui', 'drafts', 'orthogroupState']);
  const MODE_INTENT_DOMAINS = Object.freeze(['config', 'ui', 'features', 'editorState', 'groupNames', 'fileLegendCaptions']);
  // The domains a step changed: a top-level domain, or `modes.<m>.<domain>`.
  /** @param {Record<string, any>} context */
  const intentDomains = (context) => {
    const changes = Array.isArray(context.changes) && context.changes.length > 0 ? context.changes : null;
    if (!changes) {
      return new Set([...INTENT_DOMAINS, ...['circular', 'linear'].flatMap((drawingMode) => (
        MODE_INTENT_DOMAINS.map((domain) => `modes.${drawingMode}.${domain}`)
      ))]);
    }
    /** @type {Set<string>} */
    const domains = new Set();
    changes.forEach((change) => {
      const path = Array.isArray(change?.path) ? change.path : [];
      if (path[0] !== 'modes') {
        if (path[0]) domains.add(path[0]);
        return;
      }
      const drawingModes = path.length > 1 ? [path[1]] : ['circular', 'linear'];
      drawingModes.forEach((drawingMode) => {
        (path.length > 2 ? [path[2]] : MODE_INTENT_DOMAINS)
          .forEach((domain) => domains.add(`modes.${drawingMode}.${domain}`));
      });
    });
    return domains;
  };
  // The shown mode's part of a step as the composition root reads it: its
  // domains and change paths without the `modes.<mode>` prefix.
  /**
   * @param {Record<string, any>} intent
   * @param {Set<string>} domains
   * @param {any[] | undefined} changes
   * @param {'circular' | 'linear'} shownMode
   */
  const shownModeStep = (intent, domains, changes, shownMode) => {
    const prefix = `modes.${shownMode}.`;
    const shownDomains = new Set([...domains].flatMap((domain) => (
      domain.startsWith('modes.')
        ? (domain.startsWith(prefix) ? [domain.slice(prefix.length)] : [])
        : [domain]
    )));
    const shownChanges = Array.isArray(changes)
      ? changes.flatMap((change) => {
          const path = Array.isArray(change?.path) ? change.path : [];
          if (path[0] !== 'modes') return [change];
          return path[1] === shownMode ? [{ ...change, path: path.slice(2) }] : [];
        })
      : changes;
    const shownIntent = { ...intent, ...(intent.modes?.[shownMode] || {}), ui: intent.ui };
    return { shownIntent, shownDomains, shownChanges };
  };

  const applyHistoryIntent = async (intent, context = {}) => {
    if (!intent || typeof intent !== 'object') return;
    closeTransientState(state);

    const domains = intentDomains(context);
    // The alignment receipt is read against the committed Session of the
    // step's mode before anything is applied; a mode switch then shows that
    // mode with its own artifact.
    const currentMode = getGeneratedArtifactRef(state.mode, 'circular') === 'linear' ? 'linear' : 'circular';
    const stepMode = domains.has('ui') && intent.ui?.mode
      ? (intent.ui.mode === 'linear' ? 'linear' : 'circular')
      : currentMode;
    if (domains.has('alignmentState')) {
      const canonical = stepMode === currentMode || !modeCommittedSession
        ? generatedArtifactRuntimeOwner?.capture?.()?.canonical?.committedCanonicalSession
        : modeCommittedSession(stepMode);
      await validateSimilarityAlignmentResetReceipt(intent.alignmentState?.receipt, canonical && {
        ...canonical,
        renderRequest: { ...canonical.renderRequest, layout: {
          ...canonical.renderRequest.layout, similarityAlignment: intent.alignmentState?.plan
        } }
      });
    }
    await restoreMode(stepMode);
    const shownMode = getGeneratedArtifactRef(state.mode, 'circular') === 'linear' ? 'linear' : 'circular';
    const drawingModes = /** @type {const} */ (['circular', 'linear']);
    /** @param {'circular' | 'linear'} drawingMode @param {string} domain */
    const changed = (drawingMode, domain) => domains.has(`modes.${drawingMode}.${domain}`);
    const linearDrawing = drawingOfMode(state, 'linear');
    const retainedComparisonFiles = changed('linear', 'config') && !domains.has('files')
      ? new Map(
          (Array.isArray(linearDrawing.linearComparisonPlan?.edges)
            ? linearDrawing.linearComparisonPlan.edges
            : []
          )
            .map((edge) => [String(edge?.id || ''), edge?.file ?? null])
            .filter(([edgeId]) => edgeId)
        )
      : null;
    const suppressRef = state.semanticFileWatchersSuppressed;
    const previousSuppressed = getGeneratedArtifactRef(suppressRef, false);
    if (domains.has('files')) setGeneratedArtifactRef(suppressRef, true);
    try {
      if (domains.has('ui')) {
        const ui = intent.ui || {};
        if (ui.cInputType) setRef(state.cInputType, ui.cInputType);
        if (ui.lInputType) setRef(state.lInputType, ui.lInputType);
        await nextTick();
      }

      drawingModes.forEach((drawingMode) => {
        if (!changed(drawingMode, 'config')) return;
        const drawing = drawingOfMode(state, drawingMode);
        const config = intent.modes?.[drawingMode]?.config;
        if (typeof applyConfigData === 'function' && config) {
          applyConfigData(drawing, config, { resolveTrackPlacements: false });
        } else if (config?.linearComparisonPlan) {
          replaceLinearComparisonPlan(drawing.linearComparisonPlan, config.linearComparisonPlan);
        }
      });
      if (retainedComparisonFiles) {
        (linearDrawing.linearComparisonPlan?.edges || []).forEach((edge) => {
          const edgeId = String(edge?.id || '');
          if (edgeId && retainedComparisonFiles.has(edgeId)) {
            edge.file = retainedComparisonFiles.get(edgeId);
          }
        });
      }

      // The shown drawing's `ui` values are restored with the app-level ones.
      if (domains.has('ui') || changed(shownMode, 'ui')) {
        const shownUi = { ...(intent.ui || {}), ...(intent.modes?.[shownMode]?.ui || {}) };
        if (typeof applyUiStateData === 'function') {
          applyUiStateData(drawingOfMode(state, shownMode), shownUi, { restorePreviewNavigation: false });
        } else {
          applyFallbackUiStateData(state, drawingOfMode(state, shownMode), shownUi);
        }
      }
      if (domains.has('ui') && Object.hasOwn(intent.ui || {}, 'selectedFeatureRecordIdx')) {
        setRef(state.selectedFeatureRecordIdx, Number.isInteger(intent.ui.selectedFeatureRecordIdx)
          ? intent.ui.selectedFeatureRecordIdx : 0);
      }
      drawingModes.forEach((drawingMode) => {
        if (drawingMode !== shownMode && changed(drawingMode, 'ui')) {
          applyDrawingUiIntentData(drawingOfMode(state, drawingMode), intent.modes?.[drawingMode]?.ui || {});
        }
      });

      if (domains.has('files')) {
        applyFilesData(state, linearDrawing, intent.files || {}, fileStore, normalizeLinearSeqList);
      }
      if (domains.has('alignmentState')) installAlignmentState(intent.alignmentState);
      if (domains.has('drafts')) applyDraftIntentData(state, intent.drafts || {});
      drawingModes.forEach((drawingMode) => {
        const drawing = drawingOfMode(state, drawingMode);
        const modeIntent = intent.modes?.[drawingMode] || {};
        if (changed(drawingMode, 'fileLegendCaptions')) applyFileLegendCaptionData(drawing, modeIntent.fileLegendCaptions);
        if (changed(drawingMode, 'features')) applyFeatureIntentData(state, drawing, modeIntent.features || {});
        if (changed(drawingMode, 'editorState')) applyEditorIntentData(drawing, modeIntent.editorState || {});
        if (changed(drawingMode, 'groupNames')) applyGroupNameIntentData(drawing, modeIntent.groupNames || {});
      });
      if (domains.has('orthogroupState')) applyOrthogroupIntentData(state, intent.orthogroupState || {});
      await nextTick();
      if (afterApplyHistoryIntent) {
        const { shownIntent, shownDomains, shownChanges } = shownModeStep(intent, domains, context.changes, shownMode);
        await afterApplyHistoryIntent(shownIntent, { ...context, changes: shownChanges, domains: shownDomains });
      }
    } finally {
      if (domains.has('files')) setGeneratedArtifactRef(suppressRef, previousSuppressed);
    }
  };

  // A checkpoint holds the shown mode's artifact with its drawing, and the
  // other mode's drawing as an intent does (Reset Settings resets both).
  const buildArtifactCheckpoint = () => {
    const shownMode = getRef(state.mode, 'circular') === 'linear' ? 'linear' : 'circular';
    const otherMode = shownMode === 'linear' ? 'circular' : 'linear';
    const drawing = drawingOfMode(state, shownMode);
    const config = typeof buildConfigData === 'function'
      ? buildConfigData(drawing)
      : {
          form: drawing.form,
          adv: drawing.adv,
          linearComparisonPlan: cloneLinearComparisonPlanMetadata(drawing.linearComparisonPlan)
        };

    const generated = buildGeneratedArtifactSnapshot(drawing);
    const { featureCatalog = null, ...editorState } = generated.editorState || {};
    const checkpoint = cloneJsonData({
      config,
      files: buildFilesData(state, drawingOfMode(state, 'linear'), fileStore),
      alignmentState: captureAlignmentState(),
      drafts: buildDraftIntentData(state),
      fileLegendCaptions: buildFileLegendCaptionData(drawing),
      otherMode: { mode: otherMode, drawing: buildModeIntent(otherMode, false) },
      ...generated,
      editorState
    });
    checkpointFeatureCatalogs.set(checkpoint, featureCatalog);
    return checkpoint;
  };

  const applyArtifactCheckpoint = async (snapshot) => {
    if (!snapshot || typeof snapshot !== 'object') return;
    const suppressRef = state.semanticFileWatchersSuppressed;
    const previousSuppressed = getGeneratedArtifactRef(suppressRef, false);
    setGeneratedArtifactRef(suppressRef, true);
    try {
      closeTransientState(state);

      const ui = snapshot.ui || {};
      if (ui.mode) await restoreArtifactMode(ui.mode);
      const drawing = drawingOfMode(state, getRef(state.mode, 'circular'));
      if (ui.cInputType) setRef(state.cInputType, ui.cInputType);
      if (ui.lInputType) setRef(state.lInputType, ui.lInputType);
      await nextTick();

      if (typeof applyConfigData === 'function' && snapshot.config) {
        applyConfigData(drawing, snapshot.config, { resolveTrackPlacements: false });
      } else if (snapshot.config?.linearComparisonPlan) {
        replaceLinearComparisonPlan(
          drawing.linearComparisonPlan,
          snapshot.config.linearComparisonPlan
        );
      }
      applyDraftIntentData(state, snapshot.drafts || {});
      applyFileLegendCaptionData(drawing, snapshot.fileLegendCaptions);
      const other = snapshot.otherMode;
      if (other?.drawing && (other.mode === 'circular' || other.mode === 'linear')) {
        const otherDrawing = drawingOfMode(state, other.mode);
        if (typeof applyConfigData === 'function' && other.drawing.config) {
          applyConfigData(otherDrawing, other.drawing.config, { resolveTrackPlacements: false });
        }
        applyDrawingUiIntentData(otherDrawing, other.drawing.ui || {});
        applyFeatureIntentData(state, otherDrawing, other.drawing.features || {});
        applyEditorIntentData(otherDrawing, other.drawing.editorState || {});
        applyGroupNameIntentData(otherDrawing, other.drawing.groupNames || {});
        applyFileLegendCaptionData(otherDrawing, other.drawing.fileLegendCaptions);
      }

      if (typeof applyUiStateData === 'function') {
        applyUiStateData(drawing, ui, { restorePreviewNavigation: false });
      } else {
        applyFallbackUiStateData(state, drawing, ui);
      }
      await nextTick();

      applyFilesData(state, drawingOfMode(state, 'linear'), snapshot.files || {}, fileStore, normalizeLinearSeqList);
      installAlignmentState(snapshot.alignmentState);

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

  const snapshotSignature = (snapshot) => JSON.stringify(snapshot);

  return {
    applyArtifactCheckpoint,
    applyHistoryIntent,
    buildArtifactCheckpoint,
    captureGeneratedArtifactHandle,
    captureGeneratedArtifactOwnerSet,
    clearGeneratedArtifactIdentity,
    collectCurrentFileIds: () => collectCurrentFileIds(state, fileStore),
    compareGeneratedArtifactHandles,
    buildHistoryIntent,
    captureArtifactSlot,
    installArtifactSlot,
    installGeneratedArtifactOwnerSet,
    registerCapture,
    registerModeCommittedSession,
    registerModeTransition,
    restoreGeneratedArtifactHandle,
    setGeneratedArtifactIdentity,
    setGeneratedArtifactRuntimeOwner,
    setAfterApplyHistoryIntent,
    snapshotSignature
  };
};
