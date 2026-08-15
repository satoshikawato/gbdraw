import {
  normalizeFeatureVisibilityRule,
  normalizeVisibilityMode,
  splitLegacyVisibilityRules
} from '../app/feature-visibility.js';
import { serializeCleanSvg } from './svg-serialization.js';
import { cloneJsonData } from './json-clone.js';
import { replaceLayoutPreferences } from '../app/layout-preferences.js';
import {
  isResourceBackedCanonicalComparison,
  mapResourceBackedCanonicalComparison
} from './canonical-comparisons.js';

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
  if (!target || typeof target !== 'object' || !Object.prototype.hasOwnProperty.call(target, 'value')) return;
  target.value = Array.isArray(source) ? cloneJsonData(source) : [];
};

const buildFeatureIntentData = (features = {}) => ({
  selectedFeatureRecordIdx: Number.isInteger(features.selectedFeatureRecordIdx)
    ? features.selectedFeatureRecordIdx
    : 0,
  featureColorOverrides: clonePlainObject(features.featureColorOverrides),
  featureVisibilityManualRules: cloneFeatureVisibilityRules(features.featureVisibilityManualRules),
  featureVisibilityOverrides: cloneFeatureVisibilityOverrides(features.featureVisibilityOverrides),
  labelTextFeatureOverrides: clonePlainObject(features.labelTextFeatureOverrides),
  labelTextBulkOverrides: clonePlainObject(features.labelTextBulkOverrides),
  labelTextFeatureOverrideSources: clonePlainObject(features.labelTextFeatureOverrideSources),
  labelVisibilityOverrides: clonePlainObject(features.labelVisibilityOverrides),
  labelOverrideContextKey: String(features.labelOverrideContextKey || '')
});

const buildEditorIntentData = (editorState = {}) => ({
  legend: {
    entries: cloneJsonData(editorState?.legend?.entries) || [],
    deletedEntries: cloneJsonData(editorState?.legend?.deletedEntries) || [],
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
  selectedOrthogroupAlignmentFeature: String(orthogroupState.selectedOrthogroupAlignmentFeature || ''),
  orthogroupNameOverrides: clonePlainObject(orthogroupState.orthogroupNameOverrides),
  orthogroupDescriptionOverrides: clonePlainObject(orthogroupState.orthogroupDescriptionOverrides)
});

const buildUiIntentData = (ui = {}) => {
  const intent = clonePlainObject(ui);
  delete intent.generatedLegendPosition;
  delete intent.generatedMode;
  delete intent.generatedMultiRecordCanvas;
  delete intent.generatedCircularPlotTitlePosition;
  return intent;
};

const applyFeatureIntentData = (state, features = {}) => {
  setRef(
    state.selectedFeatureRecordIdx,
    Number.isInteger(features.selectedFeatureRecordIdx) ? features.selectedFeatureRecordIdx : 0
  );
  replacePlainObject(state.featureColorOverrides, clonePlainObject(features.featureColorOverrides));
  replaceFeatureVisibilityState(state, features);
  replacePlainObject(state.labelTextFeatureOverrides, clonePlainObject(features.labelTextFeatureOverrides));
  replacePlainObject(state.labelTextBulkOverrides, clonePlainObject(features.labelTextBulkOverrides));
  replacePlainObject(
    state.labelTextFeatureOverrideSources,
    clonePlainObject(features.labelTextFeatureOverrideSources)
  );
  replacePlainObject(state.labelVisibilityOverrides, clonePlainObject(features.labelVisibilityOverrides));
  setRef(state.labelOverrideContextKey, String(features.labelOverrideContextKey || ''));
};

const applyEditorIntentData = (state, editorState = {}) => {
  const legend = editorState.legend || {};
  replaceRefArray(state.legendEntries, legend.entries);
  replaceRefArray(state.deletedLegendEntries, legend.deletedEntries);
  replacePlainObject(state.legendColorOverrides, clonePlainObject(legend.colorOverrides));
  replacePlainObject(state.legendStrokeOverrides, clonePlainObject(legend.strokeOverrides));
  if (state.addedLegendCaptions?.value !== undefined) {
    state.addedLegendCaptions.value = new Set(
      Array.isArray(legend.addedCaptions) ? legend.addedCaptions.map((entry) => String(entry || '')) : []
    );
  }
  replacePlainObject(
    state.featureStrokeOverrides,
    clonePlainObject(editorState?.featureStrokes?.overrides)
  );
};

const applyOrthogroupIntentData = (state, orthogroupState = {}) => {
  setRef(state.selectedOrthogroupId, String(orthogroupState.selectedOrthogroupId || ''));
  setRef(
    state.selectedOrthogroupAlignmentFeature,
    String(orthogroupState.selectedOrthogroupAlignmentFeature || '')
  );
  replacePlainObject(state.orthogroupNameOverrides, clonePlainObject(orthogroupState.orthogroupNameOverrides));
  replacePlainObject(
    state.orthogroupDescriptionOverrides,
    clonePlainObject(orthogroupState.orthogroupDescriptionOverrides)
  );
};

const cloneLinearComparisonPlanMetadata = (plan = {}) => ({
  mode: String(plan?.mode || 'adjacent'),
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

const cloneFeatureVisibilityOverrides = (overrides) => {
  const normalized = {};
  if (!overrides || typeof overrides !== 'object' || Array.isArray(overrides)) return normalized;
  Object.entries(overrides).forEach(([featureIdRaw, modeRaw]) => {
    const featureId = String(featureIdRaw || '').trim();
    const mode = normalizeVisibilityMode(modeRaw);
    if (!featureId || mode === 'default') return;
    normalized[featureId] = mode;
  });
  return normalized;
};

const splitFeatureVisibilityState = (features = {}) => {
  if (Array.isArray(features.featureVisibilityManualRules)) {
    return {
      manualRules: cloneFeatureVisibilityRules(features.featureVisibilityManualRules),
      overrides: cloneFeatureVisibilityOverrides(features.featureVisibilityOverrides)
    };
  }
  if (Array.isArray(features.featureVisibilityRules)) {
    return splitLegacyVisibilityRules(features.featureVisibilityRules);
  }
  return {
    manualRules: [],
    overrides: cloneFeatureVisibilityOverrides(features.featureVisibilityOverrides)
  };
};

const replaceFeatureVisibilityState = (state, features = {}) => {
  const { manualRules, overrides } = splitFeatureVisibilityState(features);
  if (Array.isArray(state.featureVisibilityManualRules)) {
    state.featureVisibilityManualRules.splice(
      0,
      state.featureVisibilityManualRules.length,
      ...cloneFeatureVisibilityRules(manualRules)
    );
  }
  replacePlainObject(state.featureVisibilityOverrides, overrides);
};

const setRef = (target, value) => {
  if (target && typeof target === 'object' && Object.prototype.hasOwnProperty.call(target, 'value')) {
    target.value = value;
  }
};

const getRef = (target, fallback = null) => (
  target && typeof target === 'object' && Object.prototype.hasOwnProperty.call(target, 'value')
    ? target.value
    : fallback
);

const setGeneratedArtifactRef = (target, value) => {
  if (target && typeof target === 'object' && 'value' in target) target.value = value;
};

const getGeneratedArtifactRef = (target, fallback = null) => (
  target && typeof target === 'object' && 'value' in target ? target.value : fallback
);

const artifactOwnedValue = (value) => (
  typeof globalThis.Vue?.toRaw === 'function' ? globalThis.Vue.toRaw(value) : value
);

const buildDraftIntentData = (state) => ({
  selectedAnnotation: cloneJsonData(getRef(state.selectedAnnotation, null)),
  selectedSpecificPreset: String(getRef(state.selectedSpecificPreset, '') || ''),
  newSpecRule: clonePlainObject(state.newSpecRule),
  newPriorityRule: clonePlainObject(state.newPriorityRule),
  newColorFeat: String(getRef(state.newColorFeat, '') || ''),
  newColorVal: String(getRef(state.newColorVal, '') || ''),
  newFeatureToAdd: String(getRef(state.newFeatureToAdd, '') || ''),
  newLegendCaption: String(getRef(state.newLegendCaption, '') || ''),
  newLegendColor: String(getRef(state.newLegendColor, '') || ''),
  fileLegendCaptions: Array.from(getRef(state.fileLegendCaptions, new Set()) || [])
    .map((caption) => String(caption || '').trim())
    .filter(Boolean)
});

const applyDraftIntentData = (state, drafts = {}) => {
  setRef(state.selectedAnnotation, cloneJsonData(drafts.selectedAnnotation) || null);
  setRef(state.selectedSpecificPreset, String(drafts.selectedSpecificPreset || ''));
  replacePlainObject(state.newSpecRule, clonePlainObject(drafts.newSpecRule));
  replacePlainObject(state.newPriorityRule, clonePlainObject(drafts.newPriorityRule));
  setRef(state.newColorFeat, String(drafts.newColorFeat || ''));
  setRef(state.newColorVal, String(drafts.newColorVal || ''));
  setRef(state.newFeatureToAdd, String(drafts.newFeatureToAdd || ''));
  setRef(state.newLegendCaption, String(drafts.newLegendCaption || ''));
  setRef(state.newLegendColor, String(drafts.newLegendColor || ''));
  if (state.fileLegendCaptions?.value !== undefined) {
    state.fileLegendCaptions.value = new Set(
      Array.isArray(drafts.fileLegendCaptions)
        ? drafts.fileLegendCaptions.map((caption) => String(caption || '').trim()).filter(Boolean)
        : []
    );
  }
};

const nextFrame = () => new Promise((resolve) => {
  if (typeof window !== 'undefined' && typeof window.requestAnimationFrame === 'function') {
    window.requestAnimationFrame(() => resolve());
  } else {
    setTimeout(resolve, 0);
  }
});

const closeTransientState = (state) => {
  setRef(state.clickedFeature, null);
  setRef(state.clickedPairwiseMatch, null);
  setRef(state.clickedLabel, null);
  if (state.featureStyleScopeDialog) state.featureStyleScopeDialog.show = false;
  if (state.resetColorDialog) state.resetColorDialog.show = false;
  if (state.legendRenameDialog) state.legendRenameDialog.show = false;
  if (state.labelTextScopeDialog) state.labelTextScopeDialog.show = false;
  if (state.featureVisibilityScopeDialog) state.featureVisibilityScopeDialog.show = false;
  if (state.globalLabelModeDialog) state.globalLabelModeDialog.show = false;
  if (state.featurePopupDrag) state.featurePopupDrag.active = false;
  if (state.featurePopupResize) state.featurePopupResize.active = false;
  if (state.pairwiseMatchPopupDrag) state.pairwiseMatchPopupDrag.active = false;
  if (state.pairwiseMatchPopupResize) state.pairwiseMatchPopupResize.active = false;
};

const buildFallbackUiStateData = (state) => ({
  title: getRef(state.sessionTitle, ''),
  mode: getRef(state.mode, 'circular'),
  cInputType: getRef(state.cInputType, 'gb'),
  lInputType: getRef(state.lInputType, 'gb'),
  losatProgram: getRef(state.losatProgram, 'blastn'),
  selectedResultIndex: getRef(state.selectedResultIndex, 0),
  downloadDpi: getRef(state.downloadDpi, 300),
  canvasPadding: { ...(state.canvasPadding || {}) },
  generatedLegendPosition: getRef(state.generatedLegendPosition, 'left'),
  generatedMode: getRef(state.generatedMode, 'circular'),
  generatedMultiRecordCanvas: Boolean(getRef(state.generatedMultiRecordCanvas, false)),
  generatedCircularPlotTitlePosition: getRef(state.generatedCircularPlotTitlePosition, 'none'),
  layoutPreferences: clonePlainObject(state.layoutPreferences),
  autoLabelReflow: Boolean(getRef(state.autoLabelReflowEnabled, false)),
  paletteInstantPreviewEnabled: Boolean(getRef(state.paletteInstantPreviewEnabled, false)),
  appliedPaletteName: getRef(state.appliedPaletteName, 'default'),
  appliedPaletteColors: clonePlainObject(getRef(state.appliedPaletteColors, {})),
  pendingPaletteName: getRef(state.pendingPaletteName, ''),
  pendingPaletteColors: clonePlainObject(getRef(state.pendingPaletteColors, {})),
  legendCurrentOffset: { ...(state.legendCurrentOffset || {}) },
  diagramOffset: { ...(state.diagramOffset || {}) },
  lengthBarUserOffset: { ...(state.lengthBarUserOffset || {}) },
  plotTitleUserOffset: { ...(state.plotTitleUserOffset || {}) }
});

const applyFallbackUiStateData = (state, ui = {}) => {
  if (typeof ui.title === 'string') setRef(state.sessionTitle, ui.title);
  if (ui.mode) setRef(state.mode, ui.mode === 'linear' ? 'linear' : 'circular');
  if (ui.cInputType) setRef(state.cInputType, ui.cInputType);
  if (ui.lInputType) setRef(state.lInputType, ui.lInputType);
  if (ui.losatProgram) setRef(state.losatProgram, ui.losatProgram);
  if (ui.downloadDpi) setRef(state.downloadDpi, ui.downloadDpi);
  if (ui.generatedLegendPosition) setRef(state.generatedLegendPosition, ui.generatedLegendPosition);
  if (ui.generatedMode) setRef(state.generatedMode, ui.generatedMode);
  if (Object.prototype.hasOwnProperty.call(ui, 'generatedMultiRecordCanvas')) {
    setRef(state.generatedMultiRecordCanvas, Boolean(ui.generatedMultiRecordCanvas));
  }
  if (ui.generatedCircularPlotTitlePosition) {
    setRef(state.generatedCircularPlotTitlePosition, ui.generatedCircularPlotTitlePosition);
  }
  replaceLayoutPreferences(state.layoutPreferences, ui.layoutPreferences);
  if (state.canvasPadding && ui.canvasPadding) {
    state.canvasPadding.top = Number(ui.canvasPadding.top) || 0;
    state.canvasPadding.right = Number(ui.canvasPadding.right) || 0;
    state.canvasPadding.bottom = Number(ui.canvasPadding.bottom) || 0;
    state.canvasPadding.left = Number(ui.canvasPadding.left) || 0;
  }
  setRef(state.autoLabelReflowEnabled, Boolean(ui.autoLabelReflow));
  setRef(state.paletteInstantPreviewEnabled, Boolean(ui.paletteInstantPreviewEnabled));
  if (ui.appliedPaletteName !== undefined) setRef(state.appliedPaletteName, String(ui.appliedPaletteName || 'default'));
  if (ui.appliedPaletteColors) setRef(state.appliedPaletteColors, clonePlainObject(ui.appliedPaletteColors));
  if (ui.pendingPaletteName !== undefined) setRef(state.pendingPaletteName, String(ui.pendingPaletteName || ''));
  if (ui.pendingPaletteColors) setRef(state.pendingPaletteColors, clonePlainObject(ui.pendingPaletteColors));
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

const buildFallbackFeatureStateData = (state) => ({
  extractedFeatures: cloneJsonData(getRef(state.extractedFeatures, [])) || [],
  featureSelectorSafetyScope: cloneJsonData(getRef(state.featureSelectorSafetyScope, [])) || [],
  featureRecordIds: cloneJsonData(getRef(state.featureRecordIds, [])) || [],
  selectedFeatureRecordIdx: getRef(state.selectedFeatureRecordIdx, 0),
  featureColorOverrides: clonePlainObject(state.featureColorOverrides),
  featureVisibilityManualRules: cloneFeatureVisibilityRules(state.featureVisibilityManualRules),
  featureVisibilityOverrides: cloneFeatureVisibilityOverrides(state.featureVisibilityOverrides),
  labelTextFeatureOverrides: clonePlainObject(state.labelTextFeatureOverrides),
  labelTextBulkOverrides: clonePlainObject(state.labelTextBulkOverrides),
  labelTextFeatureOverrideSources: clonePlainObject(state.labelTextFeatureOverrideSources),
  labelVisibilityOverrides: clonePlainObject(state.labelVisibilityOverrides),
  labelOverrideContextKey: getRef(state.labelOverrideContextKey, '')
});

const applyFallbackFeatureStateData = (state, features = {}) => {
  setRef(state.extractedFeatures, cloneJsonData(features.extractedFeatures) || []);
  setRef(state.featureSelectorSafetyScope, cloneJsonData(features.featureSelectorSafetyScope) || []);
  setRef(state.featureRecordIds, cloneJsonData(features.featureRecordIds) || []);
  setRef(
    state.selectedFeatureRecordIdx,
    Number.isInteger(features.selectedFeatureRecordIdx) ? features.selectedFeatureRecordIdx : 0
  );
  replacePlainObject(state.featureColorOverrides, clonePlainObject(features.featureColorOverrides));
  replaceFeatureVisibilityState(state, features);
  replacePlainObject(state.labelTextFeatureOverrides, clonePlainObject(features.labelTextFeatureOverrides));
  replacePlainObject(state.labelTextBulkOverrides, clonePlainObject(features.labelTextBulkOverrides));
  replacePlainObject(
    state.labelTextFeatureOverrideSources,
    clonePlainObject(features.labelTextFeatureOverrideSources)
  );
  replacePlainObject(state.labelVisibilityOverrides, clonePlainObject(features.labelVisibilityOverrides));
  setRef(state.labelOverrideContextKey, String(features.labelOverrideContextKey || ''));
};

const buildFallbackOrthogroupStateData = (state) => ({
  groups: cloneJsonData(getRef(state.orthogroups, [])) || [],
  selectedOrthogroupId: getRef(state.selectedOrthogroupId, ''),
  selectedOrthogroupAlignmentFeature: getRef(state.selectedOrthogroupAlignmentFeature, ''),
  orthogroupNameOverrides: clonePlainObject(state.orthogroupNameOverrides),
  orthogroupDescriptionOverrides: clonePlainObject(state.orthogroupDescriptionOverrides)
});

const applyFallbackOrthogroupStateData = (state, data = {}) => {
  setRef(state.orthogroups, cloneJsonData(data.groups) || []);
  setRef(state.selectedOrthogroupId, String(data.selectedOrthogroupId || ''));
  setRef(state.selectedOrthogroupAlignmentFeature, String(data.selectedOrthogroupAlignmentFeature || ''));
  replacePlainObject(state.orthogroupNameOverrides, clonePlainObject(data.orthogroupNameOverrides));
  replacePlainObject(state.orthogroupDescriptionOverrides, clonePlainObject(data.orthogroupDescriptionOverrides));
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

const buildFilesData = (state, fileStore) => ({
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
  )),
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
    region_record_id: seq.region_record_id ?? '',
    region_start: seq.region_start ?? null,
    region_end: seq.region_end ?? null,
    region_reverse: Boolean(seq.region_reverse)
  })),
  linearComparisons: Array.from(state.linearComparisonPlan?.edges || []).map((edge) => ({
    id: String(edge?.id || ''),
    file: fileStore.describeValue(edge?.file)
  }))
});

const buildIntentFilesData = (state, fileStore) => {
  const files = {
    c_gb: fileStore.describeValue(state.files?.c_gb),
    c_gff: fileStore.describeValue(state.files?.c_gff),
    c_fasta: fileStore.describeValue(state.files?.c_fasta),
    c_depth: fileStore.describeValue(state.files?.c_depth),
    c_conservation_blasts_source: state.files?.c_conservation_blasts_source === 'losat-cache'
      ? 'losat-cache'
      : null,
    c_conservation_fastas: fileStore.describeValue(state.files?.c_conservation_fastas || []),
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
      region_record_id: seq.region_record_id ?? '',
      region_start: seq.region_start ?? null,
      region_end: seq.region_end ?? null,
      region_reverse: Boolean(seq.region_reverse)
    })),
    linearComparisons: Array.from(state.linearComparisonPlan?.edges || []).map((edge) => ({
      id: String(edge?.id || ''),
      file: fileStore.describeValue(edge?.file)
    }))
  };
  if (state.files?.c_conservation_blasts_source !== 'losat-cache') {
    files.c_conservation_blasts = fileStore.describeValue(state.files?.c_conservation_blasts || []);
  }
  return files;
};

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
  Array.from(state.linearComparisonPlan?.edges || []).forEach((comparison) => {
    register(comparison?.file);
  });
  return fileIds;
};

const applyFilesData = (state, filesData, fileStore, normalizeLinearSeqList = null) => {
  if (!state.files) return;
  state.matchSequenceRegistry?.reset?.();
  const restore = (value) => fileStore.restoreValue(value);
  state.files.c_gb = restore(filesData?.c_gb);
  state.files.c_gff = restore(filesData?.c_gff);
  state.files.c_fasta = restore(filesData?.c_fasta);
  state.files.c_depth = restore(filesData?.c_depth);
  if (Object.prototype.hasOwnProperty.call(filesData || {}, 'c_conservation_blasts')) {
    state.files.c_conservation_blasts = Array.isArray(filesData?.c_conservation_blasts)
      ? restore(filesData.c_conservation_blasts).filter(Boolean)
      : [];
  }
  if (Object.prototype.hasOwnProperty.call(filesData || {}, 'c_conservation_blasts_source')) {
    state.files.c_conservation_blasts_source = filesData?.c_conservation_blasts_source === 'losat-cache'
      ? 'losat-cache'
      : null;
  }
  state.files.c_conservation_fastas = Array.isArray(filesData?.c_conservation_fastas)
    ? restore(filesData.c_conservation_fastas)
    : [];
  if (Object.prototype.hasOwnProperty.call(filesData || {}, 'c_conservation_sequence_sources')) {
    state.files.c_conservation_sequence_sources = Array.isArray(filesData?.c_conservation_sequence_sources)
      ? restore(filesData.c_conservation_sequence_sources)
      : [];
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
  if (state.linearComparisonPlan && Array.isArray(state.linearComparisonPlan.edges)) {
    const comparisonFiles = new Map(
      (Array.isArray(filesData?.linearComparisons) ? filesData.linearComparisons : [])
        .map((comparison) => [String(comparison?.id || ''), restore(comparison?.file)])
    );
    state.linearComparisonPlan.edges.forEach((edge) => {
      edge.file = comparisonFiles.get(String(edge?.id || '')) || null;
    });
  }
};

export const createHistorySnapshotService = ({
  state,
  fileStore,
  nextTick = async () => {},
  normalizeLinearSeqList = null,
  buildConfigData = null,
  applyConfigData = null,
  buildUiStateData = null,
  applyUiStateData = null,
  buildCompositionIntent = null,
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

  let afterApplyHistoryIntent = null;
  let generatedArtifactRuntimeOwner = null;
  let currentGeneratedArtifactIdentity = null;
  let currentGeneratedArtifactRetainedBytes = 0;
  let generatedArtifactRestoreDepth = 0;
  let restoreSemanticSuppressionBaseline = false;
  let restoreTrustedStateBaseline = false;

  const setAfterApplyHistoryIntent = (callback) => {
    afterApplyHistoryIntent = typeof callback === 'function' ? callback : null;
  };

  const setGeneratedArtifactRuntimeOwner = (owner) => {
    generatedArtifactRuntimeOwner = (
      owner
      && typeof owner.capture === 'function'
      && typeof owner.restore === 'function'
    ) ? owner : null;
  };

  const artifactOwnerReferences = (results) => Object.freeze({
    results: Object.freeze(
      (Array.isArray(results) ? results : []).map(artifactOwnedValue)
    ),
    featureCatalog: artifactOwnedValue(getGeneratedArtifactRef(state.featureCatalog, null)),
    extractedFeatures: artifactOwnedValue(getGeneratedArtifactRef(state.extractedFeatures, null)),
    biologicalFeatures: artifactOwnedValue(getGeneratedArtifactRef(state.biologicalFeatures, null)),
    featureSelectorSafetyScope: artifactOwnedValue(
      getGeneratedArtifactRef(state.featureSelectorSafetyScope, null)
    ),
    orthogroups: artifactOwnedValue(getGeneratedArtifactRef(state.orthogroups, null)),
    featureOrthogroupIndex: artifactOwnedValue(
      getGeneratedArtifactRef(state.featureOrthogroupIndex, null)
    ),
    collinearGroups: artifactOwnedValue(getGeneratedArtifactRef(state.collinearGroups, null)),
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
    )
  });

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
    currentGeneratedArtifactIdentity = Object.freeze({
      schema: 1,
      algorithm: 'SHA-256',
      fingerprint,
      ownerReferences: artifactOwnerReferences(committedResults)
    });
    currentGeneratedArtifactRetainedBytes = Math.max(
      0,
      Number(identity.retainedBytes) || 0
    );
    return true;
  };

  const clearGeneratedArtifactIdentity = ({ retainedBytes = 0 } = {}) => {
    currentGeneratedArtifactIdentity = null;
    currentGeneratedArtifactRetainedBytes = Math.max(0, Number(retainedBytes) || 0);
  };

  const currentArtifactMatchesIdentity = (results, identity) => {
    const expected = identity?.ownerReferences;
    if (!expected || expected.results.length !== results.length) return false;
    if (!expected.results.every((result, index) => (
      result === artifactOwnedValue(results[index])
    ))) return false;
    return (
      expected.featureCatalog === artifactOwnedValue(
        getGeneratedArtifactRef(state.featureCatalog, null)
      )
      && expected.extractedFeatures === artifactOwnedValue(
        getGeneratedArtifactRef(state.extractedFeatures, null)
      )
      && expected.biologicalFeatures === artifactOwnedValue(
        getGeneratedArtifactRef(state.biologicalFeatures, null)
      )
      && expected.featureSelectorSafetyScope
        === artifactOwnedValue(getGeneratedArtifactRef(state.featureSelectorSafetyScope, null))
      && expected.orthogroups === artifactOwnedValue(
        getGeneratedArtifactRef(state.orthogroups, null)
      )
      && expected.featureOrthogroupIndex
        === artifactOwnedValue(getGeneratedArtifactRef(state.featureOrthogroupIndex, null))
      && expected.collinearGroups === artifactOwnedValue(
        getGeneratedArtifactRef(state.collinearGroups, null)
      )
      && expected.trackSlotResolvedGeometry
        === artifactOwnedValue(getGeneratedArtifactRef(state.trackSlotResolvedGeometry, null))
      && expected.proteinIdentityManifest
        === artifactOwnedValue(getGeneratedArtifactRef(state.proteinIdentityManifest, null))
      && expected.legacyProteinRawCandidates
        === artifactOwnedValue(getGeneratedArtifactRef(state.legacyProteinRawCandidates, null))
      && expected.legacyProteinDerivedEvidence
        === artifactOwnedValue(getGeneratedArtifactRef(state.legacyProteinDerivedEvidence, null))
    );
  };

  const copyMapEntries = (map) => {
    const ownedMap = artifactOwnedValue(map);
    return Object.freeze(
      Array.from(ownedMap instanceof Map ? ownedMap : new Map(), ([key, value]) => (
        Object.freeze([key, artifactOwnedValue(value)])
      ))
    );
  };

  const copyObjectEntries = (value) => Object.freeze(
    Object.entries(value && typeof value === 'object' ? value : {}).map(([key, entry]) => (
      Object.freeze([key, artifactOwnedValue(entry)])
    ))
  );

  const estimateCacheEntries = (entries) => entries.reduce((total, pair) => {
    const entry = pair?.[1];
    const rawTextBytes = typeof entry?.text === 'string' ? entry.text.length * 2 : 0;
    const declaredBytes = Number(entry?.retainedBytes || entry?.byteLength || entry?.size) || 0;
    const complexFloor = entry && typeof entry === 'object' ? 4096 : 128;
    return total + Math.max(rawTextBytes, declaredBytes, complexFloor);
  }, 0);

  const compactGeneratedArtifactSignature = (handle) => JSON.stringify({
    ui: handle.ui,
    featureEdits: {
      selectedFeatureRecordIdx: handle.features.selectedFeatureRecordIdx,
      featureRecordIds: handle.features.featureRecordIds,
      featureColorOverrides: handle.features.featureColorOverrides,
      featureVisibilityManualRules: handle.features.featureVisibilityManualRules,
      featureVisibilityOverrides: handle.features.featureVisibilityOverrides,
      labelTextFeatureOverrides: handle.features.labelTextFeatureOverrides,
      labelOverrideRows: handle.features.labelOverrideRows,
      labelTextBulkOverrides: handle.features.labelTextBulkOverrides,
      labelTextFeatureOverrideSources: handle.features.labelTextFeatureOverrideSources,
      labelVisibilityOverrides: handle.features.labelVisibilityOverrides,
      labelOverrideContextKey: handle.features.labelOverrideContextKey
    },
    editor: {
      legend: handle.editorState.legend,
      featureStrokes: handle.editorState.featureStrokes,
      originalSvgStroke: handle.editorState.originalSvgStroke
    },
    orthogroupEdits: {
      selectedOrthogroupId: handle.orthogroupState.selectedOrthogroupId,
      selectedOrthogroupAlignmentFeature:
        handle.orthogroupState.selectedOrthogroupAlignmentFeature,
      orthogroupNameOverrides: handle.orthogroupState.orthogroupNameOverrides,
      orthogroupDescriptionOverrides:
        handle.orthogroupState.orthogroupDescriptionOverrides
    },
    presentation: {
      resultPanelTab: handle.resultPanelTab
    },
    cacheEntries: {
      raw: handle.losatCache.map(([key, value]) => [
        String(key),
        String(value?.kind || ''),
        String(value?.program || ''),
        String(value?.flow || ''),
        String(value?.queryCanonicalHash || ''),
        String(value?.subjectCanonicalHash || ''),
        String(value?.text || '').length
      ]),
      derived: handle.losatDerivedCache.map(([key, value]) => [
        String(key),
        String(value?.kind || ''),
        String(value?.mode || '')
      ]),
      display: handle.losatCacheInfo
    },
    matchSequences: handle.matchSequenceSources.map((source) => [
      String(source?.key || ''),
      String(source?.recordId || ''),
      String(source?.origin || ''),
      Number(source?.recordIndex ?? -1),
      Number(source?.sourceIndex ?? -1),
      String(source?.sequence || '').length
    ])
  });

  const captureGeneratedArtifactHandle = () => {
    // Vue exposes these owners through reactive proxies. History retains their
    // raw committed values only after the production mutation audit established
    // replacement/copy-on-write ownership for every heavy domain below.
    const results = Object.freeze(
      [...(getGeneratedArtifactRef(state.results, []) || [])].map(artifactOwnedValue)
    );
    const ui = typeof buildUiStateData === 'function'
      ? buildUiStateData({ includePreviewNavigation: true })
      : buildFallbackUiStateData(state);
    const features = {
      extractedFeatures: artifactOwnedValue(
        getGeneratedArtifactRef(state.extractedFeatures, [])
      ),
      biologicalFeatures: artifactOwnedValue(
        getGeneratedArtifactRef(state.biologicalFeatures, [])
      ),
      featureSelectorSafetyScope: artifactOwnedValue(
        getGeneratedArtifactRef(state.featureSelectorSafetyScope, [])
      ),
      featureRecordIds: Object.freeze([
        ...(getGeneratedArtifactRef(state.featureRecordIds, []) || [])
      ]),
      featureVisibilitySelectorCache: copyObjectEntries(
        state.featureVisibilitySelectorCache
      ),
      selectedFeatureRecordIdx: getGeneratedArtifactRef(state.selectedFeatureRecordIdx, 0),
      featureColorOverrides: clonePlainObject(state.featureColorOverrides),
      featureVisibilityManualRules: cloneFeatureVisibilityRules(
        state.featureVisibilityManualRules
      ),
      featureVisibilityOverrides: cloneFeatureVisibilityOverrides(
        state.featureVisibilityOverrides
      ),
      labelTextFeatureOverrides: clonePlainObject(state.labelTextFeatureOverrides),
      labelOverrideRows: cloneJsonData(
        getGeneratedArtifactRef(state.canonicalLabelOverrideRows, [])
      ) || [],
      labelTextBulkOverrides: clonePlainObject(state.labelTextBulkOverrides),
      labelTextFeatureOverrideSources: clonePlainObject(
        state.labelTextFeatureOverrideSources
      ),
      labelVisibilityOverrides: clonePlainObject(state.labelVisibilityOverrides),
      labelOverrideContextKey: String(
        getGeneratedArtifactRef(state.labelOverrideContextKey, '') || ''
      )
    };
    const editorState = {
      legend: {
        entries: cloneJsonData(getGeneratedArtifactRef(state.legendEntries, [])) || [],
        deletedEntries: cloneJsonData(
          getGeneratedArtifactRef(state.deletedLegendEntries, [])
        ) || [],
        originalOrder: cloneJsonData(
          getGeneratedArtifactRef(state.originalLegendOrder, [])
        ) || [],
        originalColors: clonePlainObject(
          getGeneratedArtifactRef(state.originalLegendColors, {})
        ),
        colorOverrides: clonePlainObject(state.legendColorOverrides),
        strokeOverrides: clonePlainObject(state.legendStrokeOverrides),
        addedCaptions: Array.from(
          getGeneratedArtifactRef(state.addedLegendCaptions, new Set()) || []
        )
      },
      featureStrokes: { overrides: clonePlainObject(state.featureStrokeOverrides) },
      originalSvgStroke: cloneJsonData(
        getGeneratedArtifactRef(state.originalSvgStroke, null)
      ) || {
        color: null,
        width: null
      },
      featureCatalog: artifactOwnedValue(
        getGeneratedArtifactRef(state.featureCatalog, null)
      )
    };
    const orthogroupState = {
      groups: artifactOwnedValue(getGeneratedArtifactRef(state.orthogroups, [])),
      featureIndex: artifactOwnedValue(
        getGeneratedArtifactRef(state.featureOrthogroupIndex, new Map())
      ),
      selectedOrthogroupId: String(
        getGeneratedArtifactRef(state.selectedOrthogroupId, '') || ''
      ),
      selectedOrthogroupAlignmentFeature: String(
        getGeneratedArtifactRef(state.selectedOrthogroupAlignmentFeature, '') || ''
      ),
      orthogroupNameOverrides: clonePlainObject(state.orthogroupNameOverrides),
      orthogroupDescriptionOverrides: clonePlainObject(
        state.orthogroupDescriptionOverrides
      )
    };
    const runState = typeof buildRunStateData === 'function'
      ? buildRunStateData()
      : {
          lastRunInfo: cloneJsonData(getGeneratedArtifactRef(state.lastRunInfo, null)),
          pairwiseMatchFactors: clonePlainObject(
            getGeneratedArtifactRef(state.pairwiseMatchFactors, {})
          )
        };
    const losatCache = copyMapEntries(
      getGeneratedArtifactRef(state.losatCache, new Map())
    );
    const losatDerivedCache = copyMapEntries(
      getGeneratedArtifactRef(state.losatDerivedCache, new Map())
    );
    const handle = {
      kind: 'GeneratedArtifactHandle',
      ui,
      results,
      features,
      editorState,
      orthogroupState,
      runState,
      trackSlotResolvedGeometry: artifactOwnedValue(
        getGeneratedArtifactRef(state.trackSlotResolvedGeometry, null)
      ),
      resultGenerationKey: getGeneratedArtifactRef(state.resultGenerationKey, 0),
      resultPanelTab: getGeneratedArtifactRef(state.resultPanelTab, 'preview'),
      errorLog: cloneJsonData(getGeneratedArtifactRef(state.errorLog, null)),
      editableLabels: cloneJsonData(
        getGeneratedArtifactRef(state.editableLabels, [])
      ) || [],
      labelTextScopeDialog: cloneJsonData(state.labelTextScopeDialog || {}) || {},
      featureEditorStatus: cloneJsonData(state.featureEditorStatus || {}) || {},
      featureExtractionPending: Boolean(
        getGeneratedArtifactRef(state.featureExtractionPending, false)
      ),
      featureExtractionError: getGeneratedArtifactRef(state.featureExtractionError, null),
      labelOverrideBuildWarning: String(
        getGeneratedArtifactRef(state.labelOverrideBuildWarning, '') || ''
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
      losatCache,
      losatDerivedCache,
      losatCacheInfo: cloneJsonData(
        getGeneratedArtifactRef(state.losatCacheInfo, [])
      ) || [],
      collinearGroups: artifactOwnedValue(
        getGeneratedArtifactRef(state.collinearGroups, [])
      ),
      matchSequenceSources: Object.freeze([
        ...(state.matchSequenceRegistry?.values?.() || [])
      ]),
      runtimeState: generatedArtifactRuntimeOwner?.capture?.() || null,
      transportIdentity: null,
      baseRetainedBytes: currentGeneratedArtifactRetainedBytes,
      fileIds: Object.freeze([...collectCurrentFileIds(state, fileStore)])
    };
    const compactSignature = compactGeneratedArtifactSignature(handle);
    const transportIdentity = currentArtifactMatchesIdentity(
      results,
      currentGeneratedArtifactIdentity
    ) ? currentGeneratedArtifactIdentity : null;
    handle.transportIdentity = transportIdentity;
    const resultUtf16Bytes = results.reduce(
      (total, result) => total + String(result?.content || '').length * 2,
      0
    );
    const derivedCollectionBytes = (
      (Array.isArray(features.extractedFeatures) ? features.extractedFeatures.length : 0) * 512
      + (Array.isArray(features.biologicalFeatures) ? features.biologicalFeatures.length : 0) * 512
      + features.featureVisibilitySelectorCache.length * 512
      + (Array.isArray(orthogroupState.groups) ? orthogroupState.groups.length : 0) * 1024
      + (Array.isArray(handle.collinearGroups) ? handle.collinearGroups.length : 0) * 1024
      + handle.matchSequenceSources.reduce(
          (total, source) => total + String(source?.sequence || '').length * 2 + 256,
          0
        )
    );
    handle.identity = Object.freeze({
      fingerprint: transportIdentity?.fingerprint || '',
      compactSignature
    });
    handle.retainedBytes = currentGeneratedArtifactRetainedBytes
      + resultUtf16Bytes + derivedCollectionBytes + estimateCacheEntries(losatCache)
      + estimateCacheEntries(losatDerivedCache)
      + Math.max(0, Number(handle.runtimeState?.retainedBytes) || 0)
      + compactSignature.length * 2;
    return Object.freeze(handle);
  };

  const restoreGeneratedArtifactHandle = async (handle) => {
    if (!handle || handle.kind !== 'GeneratedArtifactHandle') return false;
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
      const ui = handle.ui || {};
      if (ui.mode) {
        setGeneratedArtifactRef(state.mode, ui.mode === 'linear' ? 'linear' : 'circular');
      }
      if (ui.cInputType) setGeneratedArtifactRef(state.cInputType, ui.cInputType);
      if (ui.lInputType) setGeneratedArtifactRef(state.lInputType, ui.lInputType);
      await nextTick();

      if (state.skipCaptureBaseConfig) state.skipCaptureBaseConfig.value = true;
      if (state.skipPositionReapply) state.skipPositionReapply.value = true;
      setGeneratedArtifactRef(state.results, [...handle.results]);
      setGeneratedArtifactRef(
        state.selectedResultIndex,
        handle.results.length > 0
          ? Math.max(0, Math.min(Number(ui.selectedResultIndex) || 0, handle.results.length - 1))
          : 0
      );

      if (typeof applyFeatureStateData === 'function') {
        applyFeatureStateData(handle.features || {});
      } else {
        applyFallbackFeatureStateData(state, handle.features || {});
      }
      replacePlainObject(
        state.featureVisibilitySelectorCache,
        Object.fromEntries(handle.features?.featureVisibilitySelectorCache || [])
      );
      setGeneratedArtifactRef(state.orthogroups, handle.orthogroupState?.groups || []);
      setGeneratedArtifactRef(
        state.featureOrthogroupIndex,
        handle.orthogroupState?.featureIndex || new Map()
      );
      setGeneratedArtifactRef(
        state.selectedOrthogroupId,
        String(handle.orthogroupState?.selectedOrthogroupId || '')
      );
      setGeneratedArtifactRef(
        state.selectedOrthogroupAlignmentFeature,
        String(handle.orthogroupState?.selectedOrthogroupAlignmentFeature || '')
      );
      replacePlainObject(
        state.orthogroupNameOverrides,
        clonePlainObject(handle.orthogroupState?.orthogroupNameOverrides)
      );
      replacePlainObject(
        state.orthogroupDescriptionOverrides,
        clonePlainObject(handle.orthogroupState?.orthogroupDescriptionOverrides)
      );

      const trustedEditorState = {
        ...cloneJsonData({
          legend: handle.editorState?.legend || {},
          featureStrokes: handle.editorState?.featureStrokes || {},
          originalSvgStroke: handle.editorState?.originalSvgStroke || {}
        }),
        featureCatalog: handle.editorState?.featureCatalog || null
      };
      if (typeof applyEditorStateData === 'function') {
        applyEditorStateData(trustedEditorState, {
          normalized: true,
          adoptCatalog: true
        });
      } else {
        applyEditorIntentData(state, trustedEditorState);
        setGeneratedArtifactRef(state.featureCatalog, trustedEditorState.featureCatalog);
      }
      if (typeof applyRunStateData === 'function') {
        applyRunStateData(handle.runState || {});
      } else {
        setGeneratedArtifactRef(
          state.lastRunInfo,
          cloneJsonData(handle.runState?.lastRunInfo) || null
        );
        setGeneratedArtifactRef(
          state.pairwiseMatchFactors,
          clonePlainObject(handle.runState?.pairwiseMatchFactors)
        );
      }

      setGeneratedArtifactRef(
        state.trackSlotResolvedGeometry,
        handle.trackSlotResolvedGeometry ?? null
      );
      setGeneratedArtifactRef(state.resultGenerationKey, handle.resultGenerationKey ?? 0);
      setGeneratedArtifactRef(state.resultPanelTab, handle.resultPanelTab || 'preview');
      setGeneratedArtifactRef(state.errorLog, cloneJsonData(handle.errorLog));
      setGeneratedArtifactRef(
        state.editableLabels,
        cloneJsonData(handle.editableLabels) || []
      );
      if (state.labelTextScopeDialog) {
        Object.assign(
          state.labelTextScopeDialog,
          cloneJsonData(handle.labelTextScopeDialog) || {}
        );
      }
      if (state.featureEditorStatus) {
        Object.assign(
          state.featureEditorStatus,
          cloneJsonData(handle.featureEditorStatus) || {}
        );
      }
      setGeneratedArtifactRef(
        state.featureExtractionPending,
        Boolean(handle.featureExtractionPending)
      );
      setGeneratedArtifactRef(
        state.featureExtractionError,
        handle.featureExtractionError ?? null
      );
      setGeneratedArtifactRef(
        state.labelOverrideBuildWarning,
        handle.labelOverrideBuildWarning || ''
      );
      setGeneratedArtifactRef(state.proteinIdentityManifest, handle.proteinIdentityManifest);
      setGeneratedArtifactRef(
        state.legacyProteinRawCandidates,
        handle.legacyProteinRawCandidates
      );
      setGeneratedArtifactRef(
        state.legacyProteinDerivedEvidence,
        handle.legacyProteinDerivedEvidence
      );
      setGeneratedArtifactRef(state.losatCache, new Map(handle.losatCache || []));
      setGeneratedArtifactRef(
        state.losatDerivedCache,
        new Map(handle.losatDerivedCache || [])
      );
      setGeneratedArtifactRef(
        state.losatCacheInfo,
        cloneJsonData(handle.losatCacheInfo) || []
      );
      setGeneratedArtifactRef(state.collinearGroups, handle.collinearGroups || []);
      if (typeof state.matchSequenceRegistry?.resetTrusted === 'function') {
        state.matchSequenceRegistry.resetTrusted(handle.matchSequenceSources || []);
      } else {
        state.matchSequenceRegistry?.reset?.(handle.matchSequenceSources || []);
      }

      if (typeof applyUiStateData === 'function') {
        applyUiStateData(ui);
      } else {
        applyFallbackUiStateData(state, ui);
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

  const buildGeneratedArtifactSnapshot = ({ includePreviewNavigation = true } = {}) => {
    const ui = typeof buildUiStateData === 'function'
      ? buildUiStateData({ includePreviewNavigation })
      : buildFallbackUiStateData(state);
    const features = typeof buildFeatureStateData === 'function'
      ? buildFeatureStateData()
      : buildFallbackFeatureStateData(state);
    const editorState = typeof buildEditorStateData === 'function'
      ? buildEditorStateData()
      : {};
    const orthogroupState = typeof buildOrthogroupStateData === 'function'
      ? buildOrthogroupStateData()
      : buildFallbackOrthogroupStateData(state);
    const results = typeof serializeResults === 'function'
      ? serializeResults()
      : buildFallbackResultsData(state);
    const runState = typeof buildRunStateData === 'function'
      ? buildRunStateData()
      : {
          lastRunInfo: cloneJsonData(getRef(state.lastRunInfo, null)),
          pairwiseMatchFactors: clonePlainObject(getRef(state.pairwiseMatchFactors, {}))
        };

    return cloneJsonData({
      ui,
      results,
      features,
      editorState,
      orthogroupState,
      runState
    });
  };

  const applyArtifactDomains = (
    snapshot,
    { trusted = false, restoreResults = true } = {}
  ) => {
    const ui = snapshot?.ui || {};
    if (restoreResults) {
      if (typeof applyResultsData === 'function') {
        applyResultsData(snapshot?.results || [], ui);
      } else {
        applyFallbackResultsData(state, snapshot?.results || []);
        applyFallbackUiStateData(state, { selectedResultIndex: ui.selectedResultIndex });
      }
    }

    if (typeof applyFeatureStateData === 'function') {
      applyFeatureStateData(snapshot?.features || {});
    } else {
      applyFallbackFeatureStateData(state, snapshot?.features || {});
    }

    if (typeof applyOrthogroupStateData === 'function') {
      applyOrthogroupStateData(snapshot?.orthogroupState || {});
    } else {
      applyFallbackOrthogroupStateData(state, snapshot?.orthogroupState || {});
    }

    if (typeof applyEditorStateData === 'function') {
      applyEditorStateData(snapshot?.editorState || {}, { trusted });
    }

    if (typeof applyRunStateData === 'function') {
      applyRunStateData(snapshot?.runState || {});
    } else {
      setRef(state.lastRunInfo, cloneJsonData(snapshot?.runState?.lastRunInfo) || null);
      setRef(state.pairwiseMatchFactors, clonePlainObject(snapshot?.runState?.pairwiseMatchFactors));
    }
  };

  const applyGeneratedArtifactSnapshot = (snapshot, { restoreResults = true } = {}) => {
    if (!snapshot || typeof snapshot !== 'object') return;
    if (state.skipCaptureBaseConfig) state.skipCaptureBaseConfig.value = true;
    if (state.skipPositionReapply) state.skipPositionReapply.value = true;
    applyArtifactDomains(snapshot, { trusted: true, restoreResults });
    if (typeof applyUiStateData === 'function') {
      applyUiStateData(snapshot.ui || {});
    } else {
      applyFallbackUiStateData(state, snapshot.ui || {});
    }
  };

  const buildHistoryIntent = async () => {
    const config = typeof buildConfigData === 'function'
      ? buildConfigData()
      : {
          form: state.form,
          adv: state.adv,
          linearComparisonPlan: cloneLinearComparisonPlanMetadata(state.linearComparisonPlan)
        };
    const ui = typeof buildUiStateData === 'function'
      ? buildUiStateData({ includePreviewNavigation: false })
      : buildFallbackUiStateData(state);
    const features = {
      selectedFeatureRecordIdx: getRef(state.selectedFeatureRecordIdx, 0),
      featureColorOverrides: state.featureColorOverrides,
      featureVisibilityManualRules: state.featureVisibilityManualRules,
      featureVisibilityOverrides: state.featureVisibilityOverrides,
      labelTextFeatureOverrides: state.labelTextFeatureOverrides,
      labelTextBulkOverrides: state.labelTextBulkOverrides,
      labelTextFeatureOverrideSources: state.labelTextFeatureOverrideSources,
      labelVisibilityOverrides: state.labelVisibilityOverrides,
      labelOverrideContextKey: getRef(state.labelOverrideContextKey, '')
    };
    const editorState = {
      legend: {
        entries: getRef(state.legendEntries, []),
        deletedEntries: getRef(state.deletedLegendEntries, []),
        colorOverrides: state.legendColorOverrides,
        strokeOverrides: state.legendStrokeOverrides,
        addedCaptions: Array.from(getRef(state.addedLegendCaptions, new Set()) || [])
      },
      featureStrokes: { overrides: state.featureStrokeOverrides }
    };
    const orthogroupState = {
      selectedOrthogroupId: getRef(state.selectedOrthogroupId, ''),
      selectedOrthogroupAlignmentFeature: getRef(state.selectedOrthogroupAlignmentFeature, ''),
      orthogroupNameOverrides: state.orthogroupNameOverrides,
      orthogroupDescriptionOverrides: state.orthogroupDescriptionOverrides
    };

    const uiIntent = buildUiIntentData(ui);
    const compositionDeltas = typeof buildCompositionIntent === 'function'
      ? cloneJsonData(buildCompositionIntent())
      : null;
    if (compositionDeltas) uiIntent.compositionUserDeltas = compositionDeltas;

    return cloneJsonData({
      config,
      files: buildIntentFilesData(state, fileStore),
      ui: uiIntent,
      drafts: buildDraftIntentData(state),
      features: buildFeatureIntentData(features),
      editorState: buildEditorIntentData(editorState),
      orthogroupState: buildOrthogroupIntentData(orthogroupState)
    });
  };

  const applyHistoryIntent = async (intent, context = {}) => {
    if (!intent || typeof intent !== 'object') return;
    closeTransientState(state);

    const domains = new Set(
      Array.isArray(context.changes) && context.changes.length > 0
        ? context.changes.map((change) => change?.path?.[0]).filter(Boolean)
        : ['config', 'files', 'ui', 'drafts', 'features', 'editorState', 'orthogroupState']
    );
    const retainedComparisonFiles = domains.has('config') && !domains.has('files')
      ? new Map(
          (Array.isArray(state.linearComparisonPlan?.edges)
            ? state.linearComparisonPlan.edges
            : []
          )
            .map((edge) => [String(edge?.id || ''), edge?.file ?? null])
            .filter(([edgeId]) => edgeId)
        )
      : null;
    const suppressRef = state.semanticFileWatchersSuppressed;
    const previousSuppressed = getRef(suppressRef, false);
    setRef(suppressRef, true);
    try {
      if (domains.has('ui')) {
        const ui = intent.ui || {};
        if (ui.mode) setRef(state.mode, ui.mode === 'linear' ? 'linear' : 'circular');
        if (ui.cInputType) setRef(state.cInputType, ui.cInputType);
        if (ui.lInputType) setRef(state.lInputType, ui.lInputType);
        await nextTick();
      }

      if (domains.has('config')) {
        if (typeof applyConfigData === 'function' && intent.config) {
          applyConfigData(intent.config);
        } else if (intent.config?.linearComparisonPlan) {
          replaceLinearComparisonPlan(state.linearComparisonPlan, intent.config.linearComparisonPlan);
        }
        if (retainedComparisonFiles) {
          (state.linearComparisonPlan?.edges || []).forEach((edge) => {
            const edgeId = String(edge?.id || '');
            if (edgeId && retainedComparisonFiles.has(edgeId)) {
              edge.file = retainedComparisonFiles.get(edgeId);
            }
          });
        }
      }

      if (domains.has('ui')) {
        if (typeof applyUiStateData === 'function') {
          applyUiStateData(intent.ui || {}, { restorePreviewNavigation: false });
        } else {
          applyFallbackUiStateData(state, intent.ui || {});
        }
        if (Number.isInteger(intent.ui?.selectedResultIndex)) {
          const resultCount = Array.isArray(getRef(state.results, []))
            ? getRef(state.results, []).length
            : 0;
          setRef(
            state.selectedResultIndex,
            resultCount > 0
              ? Math.max(0, Math.min(intent.ui.selectedResultIndex, resultCount - 1))
              : 0
          );
        }
      }

      if (domains.has('files')) {
        applyFilesData(state, intent.files || {}, fileStore, normalizeLinearSeqList);
      }
      if (domains.has('drafts')) applyDraftIntentData(state, intent.drafts || {});
      if (domains.has('features')) applyFeatureIntentData(state, intent.features || {});
      if (domains.has('editorState')) applyEditorIntentData(state, intent.editorState || {});
      if (domains.has('orthogroupState')) {
        applyOrthogroupIntentData(state, intent.orthogroupState || {});
      }
      await nextTick();
      if (afterApplyHistoryIntent) await afterApplyHistoryIntent(intent, { ...context, domains });
    } finally {
      setRef(suppressRef, previousSuppressed);
    }
  };

  const buildArtifactCheckpoint = () => {
    const config = typeof buildConfigData === 'function'
      ? buildConfigData()
      : {
          form: state.form,
          adv: state.adv,
          linearComparisonPlan: cloneLinearComparisonPlanMetadata(state.linearComparisonPlan)
        };

    return cloneJsonData({
      config,
      files: buildFilesData(state, fileStore),
      drafts: buildDraftIntentData(state),
      ...buildGeneratedArtifactSnapshot({ includePreviewNavigation: false })
    });
  };

  const applyArtifactCheckpoint = async (snapshot) => {
    if (!snapshot || typeof snapshot !== 'object') return;
    closeTransientState(state);

    const ui = snapshot.ui || {};
    if (ui.mode) setRef(state.mode, ui.mode === 'linear' ? 'linear' : 'circular');
    if (ui.cInputType) setRef(state.cInputType, ui.cInputType);
    if (ui.lInputType) setRef(state.lInputType, ui.lInputType);
    // The mode watcher clears generated metadata. Let that reset finish before
    // restoring snapshot-owned feature, label, and orthogroup state.
    await nextTick();

    if (typeof applyConfigData === 'function' && snapshot.config) {
      applyConfigData(snapshot.config);
    } else if (snapshot.config?.linearComparisonPlan) {
      replaceLinearComparisonPlan(
        state.linearComparisonPlan,
        snapshot.config.linearComparisonPlan
      );
    }
    applyDraftIntentData(state, snapshot.drafts || {});

    if (typeof applyUiStateData === 'function') {
      applyUiStateData(ui, { restorePreviewNavigation: false });
    } else {
      applyFallbackUiStateData(state, ui);
    }
    await nextTick();

    applyFilesData(state, snapshot.files || {}, fileStore, normalizeLinearSeqList);

    if (state.skipCaptureBaseConfig) state.skipCaptureBaseConfig.value = true;
    if (state.skipPositionReapply) state.skipPositionReapply.value = true;
    if (state.skipExtractOnSvgChange) state.skipExtractOnSvgChange.value = false;

    applyArtifactDomains(snapshot);

    await nextTick();
    await nextFrame();
    if (typeof applyUiStateData === 'function') {
      applyUiStateData(ui, { restorePreviewNavigation: false });
    } else {
      applyFallbackUiStateData(state, ui);
    }
  };

  const snapshotSignature = (snapshot) => JSON.stringify(snapshot);

  return {
    applyArtifactCheckpoint,
    applyGeneratedArtifactSnapshot,
    applyHistoryIntent,
    buildArtifactCheckpoint,
    captureGeneratedArtifactHandle,
    clearGeneratedArtifactIdentity,
    collectCurrentFileIds: () => collectCurrentFileIds(state, fileStore),
    compareGeneratedArtifactHandles,
    buildGeneratedArtifactSnapshot,
    buildHistoryIntent,
    restoreGeneratedArtifactHandle,
    setGeneratedArtifactIdentity,
    setGeneratedArtifactRuntimeOwner,
    setAfterApplyHistoryIntent,
    snapshotSignature
  };
};
