// @ts-check
import {
  createDefaultEditorDraftState,
  createDefaultPaletteDraftState,
  createDefaultPriorityRule,
  createDefaultSpecificRule
} from '../state.js';
/** @import { DrawingState, state as AppState } from '../state.js' */
/**
 * The applied palette a drawing's Result takes when the drawing is reset; the
 * caller writes it into that drawing's artifact (the shown one, or its slot).
 * @typedef {{ appliedPaletteName: string, appliedPaletteColors: Record<string, string> }} DrawingResetArtifact
 */
import { GROUP_NAMES_SLICE } from './drawing/group-names.js';
import { ANNOTATION_SETS_SLICE } from './drawing/annotation-sets.js';
import { VISIBILITY_RULES_SLICE } from './drawing/visibility-rules.js';
import { COLOR_RULES_SLICE } from './drawing/color-rules.js';
import { PALETTE_SLICE, createDefaultPaletteData } from './drawing/palette.js';
import { LABEL_TABLES_SLICE } from './drawing/label-tables.js';
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
 * @returns {DrawingResetArtifact}
 */
const resetPaletteState = (state, drawing) => {
  const colors = defaultPaletteColors(state);

  // The default palette, with its colors.
  PALETTE_SLICE.install(drawing.slices.palette, { ...createDefaultPaletteData(), currentColors: colors });
  return { appliedPaletteName: createDefaultPaletteDraftState().appliedPaletteName, appliedPaletteColors: { ...colors } };
};

// The rule drafts of the editor forms (one form for every drawing).
/** @param {typeof AppState} state */
const resetRuleDrafts = (state) => {
  const editorDefaults = createDefaultEditorDraftState();
  replaceReactiveObject(state.newSpecRule, createDefaultSpecificRule());
  state.selectedSpecificPreset.value = editorDefaults.selectedSpecificPreset;
  state.specificRulePresetLoading.value = editorDefaults.specificRulePresetLoading;
  replaceReactiveObject(state.newPriorityRule, createDefaultPriorityRule());
  state.newColorFeat.value = editorDefaults.newColorFeat;
  state.newColorVal.value = editorDefaults.newColorVal;
  state.newFeatureToAdd.value = editorDefaults.newFeatureToAdd;
};

/** @param {DrawingState} drawing */
const resetEditorDraftState = (drawing) => {
  VISIBILITY_RULES_SLICE.reset(drawing.slices.visibilityRules);
  // R15-14: every renamed row, a feature row's or a track row's, returns to
  // the caption Python draws, its key.
  resetPaintEdits(drawing.slices.paint);
};

// The editor's page-level state: searches, the open panel tab, and the
// clicked feature, label, and alignment feature.
/** @param {typeof AppState} state */
const resetEditorUiState = (state) => {
  const editorDefaults = createDefaultEditorDraftState();
  state.selectedAnnotation.value = null;
  state.labelReflowLastError.value = null;

  state.featureSearch.value = '';
  state.labelSearch.value = '';
  state.featurePanelTab.value = editorDefaults.featurePanelTab;
  state.clickedFeature.value = null;
  state.clickedLabel.value = null;

  state.selectedOrthogroupAlignmentFeature.value = '';
};

/**
 * The uploaded tables whose rows become settings (Default colors, Specific
 * color rules, Blacklist, Whitelist, Priority File). Reset Settings, which a
 * Session replace also runs, clears those settings, so the cards go with them
 * (OV-1267); a card that stayed would name a table its rules no longer hold.
 * @param {Record<'d_color' | 't_color' | 'blacklist' | 'whitelist' | 'qualifier_priority', unknown>} files
 */
const clearSettingTableFiles = (files) => {
  files.d_color = null;
  files.t_color = null;
  files.blacklist = null;
  files.whitelist = null;
  files.qualifier_priority = null;
};

/**
 * The project part of Reset all drawings: the one-per-Session values (LOSAT
 * execution, Rich Feature Popup, Auto Reflow, PNG DPI, palette Instant
 * Preview) return to their defaults, as 0.14.0 Reset Settings did, and so do
 * the page-level editor drafts and the setting-table cards. Reset this drawing
 * keeps them (R17-5); it writes no drawing.
 * @param {typeof AppState} state
 */
export const resetProject = (state) => {
  clearSettingTableFiles(state.cardFiles);
  Object.assign(state.losatExecution, createDefaultLosatExecution());
  state.richFeaturePopup.value = true;
  state.autoLabelReflowEnabled.value = false;
  state.downloadDpi.value = createDefaultEditorDraftState().downloadDpi;
  state.paletteInstantPreviewEnabled.value = createDefaultPaletteDraftState().paletteInstantPreviewEnabled;
  resetRuleDrafts(state);
  resetEditorUiState(state);
};

/**
 * Reset Settings of one drawing: its settings and edits return to a new
 * drawing's of its mode (PD-OI-070); no other drawing and no one-per-Session
 * value changes. Reset all drawings runs `resetProject` and this for every
 * drawing; the caller runs the mode's follow-ups and writes the returned
 * applied palette into the drawing's artifact.
 * @param {typeof AppState} state
 * @param {DrawingState} drawing
 * @returns {DrawingResetArtifact}
 */
export const resetSettings = (state, drawing) => {
  SETTINGS_SLICE.reset(drawing.slices.settings);
  ANNOTATION_SETS_SLICE.reset(drawing.slices.annotationSets);
  FEATURE_EDITS_SLICE.reset(drawing.slices.featureEdits);
  // The label tables (the bulk label edits, a loaded Session's saved label
  // table, and the filters) and the Similarity group names, the dormant
  // ones included, return to a new drawing's (R15-21).
  LABEL_TABLES_SLICE.reset(drawing.slices.labelTables);
  GROUP_NAMES_SLICE.reset(drawing.slices.groupNames);
  COLOR_RULES_SLICE.reset(drawing.slices.colorRules);
  resetEditorDraftState(drawing);
  return resetPaletteState(state, drawing);
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
