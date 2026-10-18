// @ts-check
import {
  createDefaultEditorDraftState,
  createDefaultPaletteDraftState,
  createDefaultPriorityRule,
  createDefaultSpecificRule
} from '../state.js';
/** @import { DrawingState, state as AppState } from '../state.js' */
import { GROUP_NAMES_SLICE } from './drawing/group-names.js';
import { ANNOTATION_SETS_SLICE } from './drawing/annotation-sets.js';
import { VISIBILITY_RULES_SLICE } from './drawing/visibility-rules.js';
import { COLOR_RULES_SLICE } from './drawing/color-rules.js';
import { PALETTE_SLICE, createDefaultPaletteData } from './drawing/palette.js';
import { clearLabelTextEdits, resetLabelFilters } from './drawing/label-tables.js';
import { FEATURE_EDITS_SLICE } from './drawing/feature-edits.js';
import { resetPaintEdits } from './drawing/paint.js';
import { SETTINGS_SLICE } from './drawing/settings.js';
import { createDefaultLosatExecution } from './session-active-config-contract.js';
import { normalizePaletteColors } from '../utils/color-utils.js';

// Every caller passes plain data (defaults, palette colors, JSON-like rows), so
// the copy has the shape of the value.
/**
 * @template T
 * @param {T} value
 * @returns {T}
 */
const clonePlain = (value) => {
  if (Array.isArray(value)) return /** @type {T} */ (value.map((entry) => clonePlain(entry)));
  if (!value || typeof value !== 'object') return value;
  return /** @type {T} */ (Object.fromEntries(
    Object.entries(value).map(([key, entry]) => [key, clonePlain(entry)])
  ));
};

/**
 * @param {Record<string, unknown>} target
 * @param {Record<string, unknown>} source
 */
const replaceReactiveObject = (target, source) => {
  Object.keys(target).forEach((key) => {
    delete target[key];
  });
  Object.entries(source).forEach(([key, value]) => {
    target[key] = clonePlain(value);
  });
};

/** @param {Pick<typeof AppState, 'paletteDefinitions'>} state */
const defaultPaletteColors = (state) => {
  const definitions = state.paletteDefinitions?.value || {};
  const colors = definitions.default && typeof definitions.default === 'object'
    ? definitions.default
    : {};
  return normalizePaletteColors(clonePlain(colors));
};

/**
 * @param {typeof AppState} state
 * @param {DrawingState} drawing
 */
const resetPaletteState = (state, drawing) => {
  const defaults = createDefaultPaletteDraftState();
  const colors = defaultPaletteColors(state);

  // The default palette, with its colors.
  PALETTE_SLICE.install(drawing.slices.palette, { ...createDefaultPaletteData(), currentColors: colors });
  state.paletteInstantPreviewEnabled.value = defaults.paletteInstantPreviewEnabled;
  state.appliedPaletteName.value = defaults.appliedPaletteName;
  state.appliedPaletteColors.value = { ...colors };
};


/**
 * @param {typeof AppState} state
 * @param {DrawingState} drawing
 */
const resetRuleDraftState = (state, drawing) => {
  const editorDefaults = createDefaultEditorDraftState();
  resetLabelFilters(drawing.slices.labelTables);
  COLOR_RULES_SLICE.reset(drawing.slices.colorRules);
  replaceReactiveObject(state.newSpecRule, createDefaultSpecificRule());
  state.selectedSpecificPreset.value = editorDefaults.selectedSpecificPreset;
  state.specificRulePresetLoading.value = editorDefaults.specificRulePresetLoading;
  replaceReactiveObject(state.newPriorityRule, createDefaultPriorityRule());
  state.newColorFeat.value = editorDefaults.newColorFeat;
  state.newColorVal.value = editorDefaults.newColorVal;
  state.newFeatureToAdd.value = editorDefaults.newFeatureToAdd;
};

/**
 * @param {typeof AppState} state
 * @param {DrawingState} drawing
 */
const resetEditorDraftState = (state, drawing) => {
  const editorDefaults = createDefaultEditorDraftState();
  VISIBILITY_RULES_SLICE.reset(drawing.slices.visibilityRules);
  // R15-14: every renamed row, a feature row's or a track row's, returns to
  // the caption Python draws, its key.
  resetPaintEdits(drawing.slices.paint);

  // The saved label table of a loaded Session goes with the bulk label edits,
  // also when the bulk map is already empty (R15-21).
  clearLabelTextEdits(drawing.slices.labelTables);
  state.autoLabelReflowEnabled.value = false;
  state.labelReflowLastError.value = null;

  state.featureSearch.value = '';
  state.labelSearch.value = '';
  state.featurePanelTab.value = editorDefaults.featurePanelTab;
  state.downloadDpi.value = editorDefaults.downloadDpi;
  state.clickedFeature.value = null;
  state.clickedLabel.value = null;

  state.selectedOrthogroupAlignmentFeature.value = '';
  // Reset clears the dormant group names too (R15-21).
  GROUP_NAMES_SLICE.reset(drawing.slices.groupNames);
};

// Reset Settings returns both drawings to their own mode's defaults
// (PD-OI-070), and the app-level LOSAT execution and popup settings.
/** @param {typeof AppState} state */
export const resetSettings = (state) => {
  Object.assign(state.losatExecution, createDefaultLosatExecution());
  state.richFeaturePopup.value = true;
  /** @type {DrawingState[]} */ (Object.values(state.drawings)).forEach((drawing) => {
    SETTINGS_SLICE.reset(drawing.slices.settings);
    ANNOTATION_SETS_SLICE.reset(drawing.slices.annotationSets);
    FEATURE_EDITS_SLICE.reset(drawing.slices.featureEdits);
    state.selectedAnnotation.value = null;
    state.suppressCircularMultiRecordDefaults.value = false;
    resetPaletteState(state, drawing);
    resetRuleDraftState(state, drawing);
    resetEditorDraftState(state, drawing);
  });
};

/** @param {Pick<typeof AppState, 'zoom' | 'isPanning' | 'panStart' | 'canvasPan' | 'showCanvasControls'>} state */
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
