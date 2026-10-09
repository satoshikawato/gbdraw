// @ts-check
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
/** @import { DrawingState } from '../state.js' */
import { createDefaultLayoutPreferences } from './layout-preferences.js';
import { createDefaultLosatExecution } from './session-active-config-contract.js';
import { normalizePaletteColors } from '../utils/color-utils.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';

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
  return normalizePaletteColors(clonePlain(colors));
};

/**
 * @param {Record<string, any>} state
 * @param {DrawingState} drawing
 */
const resetPaletteState = (state, drawing) => {
  const defaults = createDefaultPaletteDraftState();
  const colors = defaultPaletteColors(state);

  drawing.selectedPalette.value = defaults.selectedPalette;
  state.paletteInstantPreviewEnabled.value = defaults.paletteInstantPreviewEnabled;
  drawing.currentColors.value = colors;
  state.appliedPaletteName.value = defaults.appliedPaletteName;
  state.appliedPaletteColors.value = { ...colors };
  drawing.pendingPaletteName.value = defaults.pendingPaletteName;
  drawing.pendingPaletteColors.value = defaults.pendingPaletteColors;
};

/**
 * @param {Record<string, any>} state
 * @param {DrawingState} drawing
 */
const resetLayoutPreferenceState = (state, drawing) => {
  const defaults = createDefaultLayoutPreferences();
  Object.assign(drawing.layoutPreferences.circular.single, defaults.circular.single);
  Object.assign(drawing.layoutPreferences.circular.multi, defaults.circular.multi);
  Object.assign(drawing.layoutPreferences.linear, defaults.linear);
  state.suppressCircularMultiRecordDefaults.value = false;
};

/**
 * @param {Record<string, any>} state
 * @param {DrawingState} drawing
 */
const resetRuleDraftState = (state, drawing) => {
  const labelDefaults = createDefaultLabelFilterState();
  const editorDefaults = createDefaultEditorDraftState();
  drawing.filterMode.value = labelDefaults.filterMode;
  drawing.manualBlacklist.value = labelDefaults.manualBlacklist;
  replaceReactiveArray(drawing.manualWhitelist, labelDefaults.manualWhitelist);
  replaceReactiveArray(drawing.manualSpecificRules);
  replaceReactiveObject(state.newSpecRule, createDefaultSpecificRule());
  state.selectedSpecificPreset.value = editorDefaults.selectedSpecificPreset;
  state.specificRulePresetLoading.value = editorDefaults.specificRulePresetLoading;
  replaceReactiveArray(drawing.manualPriorityRules);
  replaceReactiveObject(state.newPriorityRule, createDefaultPriorityRule());
  state.newColorFeat.value = editorDefaults.newColorFeat;
  state.newColorVal.value = editorDefaults.newColorVal;
  state.newFeatureToAdd.value = editorDefaults.newFeatureToAdd;
};

/**
 * @param {Record<string, any>} state
 * @param {DrawingState} drawing
 */
const resetEditorDraftState = (state, drawing) => {
  const editorDefaults = createDefaultEditorDraftState();
  clearReactiveObject(drawing.featureColorOverrides);
  replaceReactiveArray(drawing.featureVisibilityManualRules);
  clearReactiveObject(drawing.featureOverrides);
  clearReactiveObject(drawing.featureStrokeOverrides);
  clearReactiveObject(drawing.legendColorOverrides);
  clearReactiveObject(drawing.legendStrokeOverrides);
  drawing.deletedLegendEntries.value = [];
  drawing.dormantLegendEntries.value = [];
  drawing.addedLegendCaptions.value = new Set();
  drawing.fileLegendCaptions.value = new Set();

  clearReactiveObject(drawing.labelTextBulkOverrides);
  state.autoLabelReflowEnabled.value = false;
  state.labelReflowLastError.value = null;

  state.featureSearch.value = '';
  state.labelSearch.value = '';
  state.featurePanelTab.value = editorDefaults.featurePanelTab;
  state.downloadDpi.value = editorDefaults.downloadDpi;
  state.clickedFeature.value = null;
  state.clickedLabel.value = null;

  state.selectedOrthogroupAlignmentFeature.value = '';
  clearReactiveObject(drawing.orthogroupNameOverrides);
  clearReactiveObject(drawing.orthogroupDescriptionOverrides);
};

/** @param {DrawingState} drawing */
const resetLinearComparisonPlan = (drawing) => {
  const plan = drawing.linearComparisonPlan;
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

// Reset Settings returns both drawings to their own mode's defaults
// (PD-OI-070), and the app-level LOSAT execution and popup settings.
/** @param {Record<string, any>} state */
export const resetSettings = (state) => {
  Object.assign(state.losatExecution, createDefaultLosatExecution());
  state.richFeaturePopup.value = true;
  /** @type {[string, DrawingState][]} */ (Object.entries(state.drawings)).forEach(([mode, drawing]) => {
    replaceReactiveObject(drawing.form, createDefaultForm());
    drawing.linearRecordLayoutEnabled.value = WEB_UX_PROFILE.linear.arrangeInRowsByDefault;
    replaceReactiveObject(drawing.adv, createDefaultAdv(mode));
    drawing.linearTypographyLinked.value = true;
    replaceReactiveObject(drawing.losat, createDefaultLosat());
    replaceReactiveObject(drawing.circularConservation, createDefaultCircularConservation());
    clearReactiveObject(drawing.unmanagedConfigOverrides);
    replaceReactiveArray(drawing.annotationSets);
    replaceReactiveArray(drawing.recordDisplayDrafts);
    Object.keys(drawing.featurePlacementOverrides).forEach((key) => delete drawing.featurePlacementOverrides[key]);
    state.selectedAnnotation.value = null;
    resetLinearComparisonPlan(drawing);
    drawing.losatProgram.value = 'blastn';

    resetLayoutPreferenceState(state, drawing);
    resetPaletteState(state, drawing);
    resetRuleDraftState(state, drawing);
    resetEditorDraftState(state, drawing);
  });
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
