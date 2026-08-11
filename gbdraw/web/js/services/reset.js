import {
  createDefaultAdv,
  createDefaultCircularConservation,
  createDefaultEditorDraftState,
  createDefaultForm,
  createDefaultLabelFilterState,
  createDefaultLosat,
  createDefaultPaletteDraftState,
  createDefaultPriorityRule,
  createDefaultSpecificRule
} from '../state.js';
import { createDefaultLayoutPreferences } from '../app/layout-preferences.js';

const clonePlain = (value) => {
  if (Array.isArray(value)) return value.map((entry) => clonePlain(entry));
  if (!value || typeof value !== 'object') return value;
  return Object.fromEntries(
    Object.entries(value).map(([key, entry]) => [key, clonePlain(entry)])
  );
};

const replaceReactiveObject = (target, source) => {
  Object.keys(target).forEach((key) => {
    delete target[key];
  });
  Object.entries(source).forEach(([key, value]) => {
    target[key] = clonePlain(value);
  });
};

const replaceReactiveArray = (target, source = []) => {
  target.splice(0, target.length, ...source.map((entry) => clonePlain(entry)));
};

const clearReactiveObject = (target) => {
  Object.keys(target).forEach((key) => {
    delete target[key];
  });
};

const defaultPaletteColors = (state) => {
  const definitions = state.paletteDefinitions?.value || {};
  const colors = definitions.default && typeof definitions.default === 'object'
    ? definitions.default
    : {};
  return state.normalizePaletteColors(clonePlain(colors));
};

const resetPaletteState = (state) => {
  const defaults = createDefaultPaletteDraftState();
  const colors = defaultPaletteColors(state);

  state.selectedPalette.value = defaults.selectedPalette;
  state.paletteInstantPreviewEnabled.value = defaults.paletteInstantPreviewEnabled;
  state.currentColors.value = colors;
  state.appliedPaletteName.value = defaults.appliedPaletteName;
  state.appliedPaletteColors.value = { ...colors };
  state.pendingPaletteName.value = defaults.pendingPaletteName;
  state.pendingPaletteColors.value = defaults.pendingPaletteColors;
};

const resetLayoutPreferenceState = (state) => {
  const defaults = createDefaultLayoutPreferences();
  Object.assign(state.layoutPreferences.circular.single, defaults.circular.single);
  Object.assign(state.layoutPreferences.circular.multi, defaults.circular.multi);
  Object.assign(state.layoutPreferences.linear, defaults.linear);
  state.suppressCircularMultiRecordDefaults.value = false;
};

const resetRuleDraftState = (state) => {
  const labelDefaults = createDefaultLabelFilterState();
  const editorDefaults = createDefaultEditorDraftState();
  state.filterMode.value = labelDefaults.filterMode;
  state.manualBlacklist.value = labelDefaults.manualBlacklist;
  replaceReactiveArray(state.manualWhitelist, labelDefaults.manualWhitelist);
  replaceReactiveArray(state.manualSpecificRules);
  replaceReactiveObject(state.newSpecRule, createDefaultSpecificRule());
  state.selectedSpecificPreset.value = editorDefaults.selectedSpecificPreset;
  state.specificRulePresetLoading.value = editorDefaults.specificRulePresetLoading;
  replaceReactiveArray(state.manualPriorityRules);
  replaceReactiveObject(state.newPriorityRule, createDefaultPriorityRule());
  state.newColorFeat.value = editorDefaults.newColorFeat;
  state.newColorVal.value = editorDefaults.newColorVal;
  state.newFeatureToAdd.value = editorDefaults.newFeatureToAdd;
};

const resetEditorDraftState = (state) => {
  const editorDefaults = createDefaultEditorDraftState();
  clearReactiveObject(state.featureColorOverrides);
  replaceReactiveArray(state.featureVisibilityManualRules);
  clearReactiveObject(state.featureVisibilityOverrides);
  clearReactiveObject(state.featureVisibilitySelectorCache);
  clearReactiveObject(state.featureStrokeOverrides);
  clearReactiveObject(state.legendColorOverrides);
  clearReactiveObject(state.legendStrokeOverrides);
  state.deletedLegendEntries.value = [];
  state.newLegendCaption.value = editorDefaults.newLegendCaption;
  state.newLegendColor.value = editorDefaults.newLegendColor;
  state.addedLegendCaptions.value = new Set();
  state.fileLegendCaptions.value = new Set();

  clearReactiveObject(state.labelTextFeatureOverrides);
  clearReactiveObject(state.labelTextBulkOverrides);
  clearReactiveObject(state.labelTextFeatureOverrideSources);
  clearReactiveObject(state.labelVisibilityOverrides);
  state.labelOverrideContextKey.value = '';
  state.labelOverrideBuildWarning.value = '';
  state.autoLabelReflowEnabled.value = false;
  state.labelReflowLastError.value = null;
  state.labelLayoutDirtyReason.value = '';

  state.featureSearch.value = '';
  state.labelSearch.value = '';
  state.featurePanelTab.value = editorDefaults.featurePanelTab;
  state.downloadDpi.value = editorDefaults.downloadDpi;
  state.clickedFeature.value = null;
  state.clickedLabel.value = null;

  state.selectedOrthogroupAlignmentFeature.value = '';
  clearReactiveObject(state.orthogroupNameOverrides);
  clearReactiveObject(state.orthogroupDescriptionOverrides);
};

const resetLinearComparisonPlan = (state) => {
  const plan = state.linearComparisonPlan;
  if (!plan || typeof plan !== 'object') return;
  const retainedEdges = (Array.isArray(plan.edges) ? plan.edges : [])
    .filter((edge) => Boolean(edge?.file) || String(edge?.losatFilename || '').trim())
    .map((edge) => ({
      ...edge,
      included: false,
      fileActive: false,
      losatFilenameActive: false
    }));
  plan.mode = 'none';
  plan.defaultSource = 'losat';
  if (!Array.isArray(plan.edges)) plan.edges = [];
  plan.edges.splice(0, plan.edges.length, ...retainedEdges);
};

export const resetSettings = (state) => {
  replaceReactiveObject(state.form, createDefaultForm());
  replaceReactiveObject(state.adv, createDefaultAdv(state.mode.value));
  state.modeProfileStateManager?.reset?.(state.mode.value, state.adv);
  replaceReactiveObject(state.losat, createDefaultLosat());
  replaceReactiveObject(state.circularConservation, createDefaultCircularConservation());
  replaceReactiveArray(state.annotationSets);
  state.selectedAnnotation.value = null;
  resetLinearComparisonPlan(state);
  state.losatProgram.value = 'blastn';

  resetLayoutPreferenceState(state);
  resetPaletteState(state);
  resetRuleDraftState(state);
  resetEditorDraftState(state);
};

export const resetLayoutState = (state) => {
  state.zoom.value = 1.0;
  state.isPanning.value = false;
  state.panStart.x = 0;
  state.panStart.y = 0;
  state.panStart.panX = 0;
  state.panStart.panY = 0;
  state.canvasPan.x = 0;
  state.canvasPan.y = 0;
  state.showCanvasControls.value = false;
};
