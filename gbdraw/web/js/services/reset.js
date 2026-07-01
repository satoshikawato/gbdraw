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
  state.circularLegendPosition.value = 'left';
  state.linearLegendPosition.value = 'bottom';
  state.circularPlotTitlePosition.value = 'none';
  state.linearPlotTitlePosition.value = 'bottom';
  state.circularSingleRecordLegendPosition.value = 'left';
  state.circularSingleRecordPlotTitlePosition.value = 'none';
  state.circularMultiRecordLegendPosition.value = null;
  state.circularMultiRecordPlotTitlePosition.value = null;
  state.suppressCircularMultiRecordDefaults.value = false;

  if (state.mode.value === 'linear') {
    state.form.legend = 'bottom';
    state.adv.plot_title_position = 'bottom';
  } else {
    state.form.legend = 'left';
    state.adv.plot_title_position = 'none';
  }
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
  replaceReactiveArray(state.featureVisibilityRules);
  clearReactiveObject(state.featureVisibilityOverrides);
  clearReactiveObject(state.featureStrokeOverrides);
  clearReactiveObject(state.legendColorOverrides);
  clearReactiveObject(state.legendStrokeOverrides);
  state.deletedLegendEntries.value = [];
  state.newLegendCaption.value = editorDefaults.newLegendCaption;
  state.newLegendColor.value = editorDefaults.newLegendColor;
  state.addedLegendCaptions.value = new Set();

  clearReactiveObject(state.labelTextFeatureOverrides);
  clearReactiveObject(state.labelTextBulkOverrides);
  clearReactiveObject(state.labelTextFeatureOverrideSources);
  clearReactiveObject(state.labelVisibilityOverrides);
  state.labelOverrideContextKey.value = '';
  state.labelOverrideBuildWarning.value = '';
  state.autoLabelReflowEnabled.value = false;
  state.labelReflowLastError.value = null;

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

export const resetSettings = (state) => {
  replaceReactiveObject(state.form, createDefaultForm());
  replaceReactiveObject(state.adv, createDefaultAdv());
  replaceReactiveObject(state.losat, createDefaultLosat());
  replaceReactiveObject(state.circularConservation, createDefaultCircularConservation());
  state.blastSource.value = 'losat';
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
