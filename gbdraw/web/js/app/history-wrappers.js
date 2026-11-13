// @ts-check
// The editor's actions as History steps: each wrapper runs one action as one
// undoable step under its label (`undoableAction`), as a paint-intent edit
// that shows what changed (`editEditorIntent`, the root's), or as a dialog
// choice that is one step and records none on Cancel (`scopeChoiceWithHistory`).
// The wrappers here are those of the feature visibility and color edits and
// their scope dialogs, the Legend rename and reorder and the palette and
// Reset color dialogs, and the clicked feature's stroke, the selected
// features' bulk edits and the label text edits. The composition root passes
// History's step runners, the dialog choice's History wrapper, and the editor
// owners' actions as ports (R11: each step keeps its label, `joinsOpenStep`
// and the calls inside it).
import { countUnresolvedFeatureEdits, removeUnresolvedFeatureEdits } from '../services/feature-visibility.js';
import { legendEntryKey } from '../services/drawing/paint.js';
import { msg } from '../utils/i18n.js';
/** @import { WebState } from '../state.js' */
/** @import { KeptLegendRow } from '../services/drawing/paint.js' */
/** @import { UnresolvedNotice, VisibilityTarget } from '../services/feature-visibility.js' */
/** @import { LegendRowName } from '../services/legend-row-names.js' */
/** @import { ColorRuleLegendGroup } from '../services/specific-color-rules.js' */
/** @import { PythonLegendKey } from '../services/legend-svg.js' */
/** @import { PopupOrigin } from '../utils/popup-bounds.js' */
/** @import { ColorFeature } from './feature-editor/color-actions.js' */
/** @import { LabelTableInput } from './feature-editor/label-table.js' */
/** @import { EditableLabelEntry } from '../services/label-override-table.js' */
/** @import { VisibilityScopeChoice } from './feature-editor/visibility-actions.js' */

const { computed, reactive, ref } = window.Vue;

/**
 * History's undoable step: `fn` runs as one step under `label`; `joinsOpenStep` joins the step already open.
 * @typedef {(label: string, fn: () => unknown, options?: { joinsOpenStep?: boolean }) => Promise<unknown>} RunUndoable
 */
/**
 * An action as one undoable step under `label` (`undoableAction`).
 * @typedef {<A extends unknown[]>(label: string, fn: (...args: A) => unknown, options?: { joinsOpenStep?: boolean }) => (...args: A) => Promise<unknown>} UndoableAction
 */
/**
 * A live edit of the paint intent as one History step that shows the domains it changed (the root's `editEditorIntent`).
 * @typedef {<A extends unknown[]>(label: string, mutate: (...args: A) => unknown, options?: { domains?: readonly string[] }) => (...args: A) => Promise<unknown>} EditEditorIntent
 */
/**
 * The action that writes the paint intent and, when it changed, shows `domains` (the root's `showingEditorIntent`).
 * @typedef {<A extends unknown[]>(mutate: (...args: A) => unknown, domains: readonly string[]) => (...args: A) => Promise<unknown>} ShowingEditorIntent
 */
/** @typedef {(choice: string, ...rest: unknown[]) => unknown} ChoiceHandler A dialog's handler of the chosen option. */
/**
 * The dialog choice's History wrapper: the choice is one step, Cancel records none (`createDialogChoice`).
 * @typedef {(label: () => string, handler: ChoiceHandler, cancel: () => unknown) => ChoiceHandler} WithHistory
 */
/** @typedef {(label: () => string, handler: ChoiceHandler, cancel?: () => void) => ChoiceHandler} ScopeChoiceWithHistory */

/**
 * @typedef {object} FeatureEditHistoryPorts
 * @property {WebState} state
 * @property {RunUndoable} runUndoable
 * @property {WithHistory} withHistory
 * @property {EditEditorIntent} editEditorIntent
 * @property {ShowingEditorIntent} showingEditorIntent
 * @property {readonly string[]} LEGEND_FILL_DOMAINS The paint domains a Legend fill edit shows.
 * @property {readonly string[]} STROKE_DOMAINS The paint domains a stroke edit shows.
 * @property {() => unknown} addFeatureVisibilityRule
 * @property {(index: number) => unknown} moveFeatureVisibilityRuleDown
 * @property {(index: number) => unknown} moveFeatureVisibilityRuleUp
 * @property {(index: number) => unknown} removeFeatureVisibilityRule
 * @property {(index: number, field: string, value: unknown) => unknown} setFeatureVisibilityRuleField
 * @property {(feat: VisibilityTarget, modeRaw: string, options?: { triggerReflow?: boolean, scope?: VisibilityScopeChoice }) => unknown} setFeatureVisibility
 * @property {(modeRaw: unknown) => unknown} updateClickedFeatureVisibility
 * @property {(scopeId: string) => unknown} handleFeatureVisibilityScopeChoice
 * @property {(feat: ColorFeature, color: string, requestedLegendName?: string | null) => unknown} requestFeatureColorChange
 * @property {(feature: ColorFeature, value: string | null, customCaption?: LegendRowName | null) => unknown} setFeatureColorValue
 * @property {(color: string) => unknown} updateClickedFeatureColor
 * @property {(index: number, color: string) => unknown} updateLegendEntryColor
 * @property {(index: number) => boolean} legendRowHasRules
 * @property {(key: PythonLegendKey, type: string, color: string) => unknown} setLegendRowDefaultColor Sets the default
 *   color of the type whose default color draws the Legend row of `key`, which releases its own color.
 * @property {(key: string) => unknown} preparedLegendRowDefaultColorType
 *   The type whose default color draws a Legend row, null when none does, undefined when the matches went stale.
 * @property {(idx: number, value: string | null) => unknown} setLegendEntryStrokeColorValue
 * @property {() => void} cancelFeatureStyleScope
 * @property {(choice: string) => unknown} handleColorScopeChoice
 * @property {(choice: string) => unknown} handleFeatureStyleScopeChoice
 */

/**
 * The feature visibility and color edits and their scope dialogs as History
 * steps. Registers the `featureIdentityNoticeSummary` and
 * `featureStyleScopeLabelOptions` computeds; the root calls it where they stood.
 * @param {FeatureEditHistoryPorts} ports
 */
export const createFeatureEditHistory = ({
  state, runUndoable, withHistory, editEditorIntent, showingEditorIntent, LEGEND_FILL_DOMAINS, STROKE_DOMAINS,
  addFeatureVisibilityRule, moveFeatureVisibilityRuleDown, moveFeatureVisibilityRuleUp, removeFeatureVisibilityRule,
  setFeatureVisibilityRuleField, setFeatureVisibility, updateClickedFeatureVisibility, handleFeatureVisibilityScopeChoice,
  requestFeatureColorChange, setFeatureColorValue, updateClickedFeatureColor, updateLegendEntryColor, legendRowHasRules,
  setLegendRowDefaultColor, preparedLegendRowDefaultColorType, setLegendEntryStrokeColorValue, cancelFeatureStyleScope,
  handleColorScopeChoice, handleFeatureStyleScopeChoice
}) => {
  const { featureIdentityNotices, featureEditRemovalCount, featureStyleScopeDialog } = state;
  /**
   * @template {unknown[]} A
   * @template R
   * @param {string} label
   * @param {(...args: A) => R} fn
   * @param {{ joinsOpenStep?: boolean }} [options] History's options for the step.
   */
  const undoableAction = (label, fn, options) => (/** @type {A} */ ...args) => runUndoable(label, () => fn(...args), options);
  // Python's report of the per-feature edits the Results do not draw (design
  // Q4 3.4): dormant edits outside the crop or display, edits whose feature
  // the source does not have, and the edits a source-replacing Generate removed.
  const featureIdentityNoticeSummary = computed(() => {
    const drawing = state.activeDrawing();
    const notices = Array.isArray(featureIdentityNotices.value) ? featureIdentityNotices.value : [];
    const summary = {
      dormant: notices.filter((/** @type {UnresolvedNotice} */ notice) => notice.status !== 'unresolved').length,
      unmatched: countUnresolvedFeatureEdits({ featureEdits: drawing.slices.featureEdits, notices }),
      removed: Number(featureEditRemovalCount.value) || 0
    };
    return summary.dormant || summary.unmatched || summary.removed ? summary : null;
  });
  // An explicit removal of the edits Python reported unresolved (R2).
  const removeUnmatchedFeatureEdits = undoableAction(msg('Remove unmatched feature edits'), () => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability();
    if (busy) return busy;
    return removeUnresolvedFeatureEdits({
      featureEdits: drawing.slices.featureEdits, notices: featureIdentityNotices.value
    }) > 0;
  });
  const addFeatureVisibilityRuleWithHistory = undoableAction(msg('Add feature visibility rule'), addFeatureVisibilityRule);
  const moveFeatureVisibilityRuleDownWithHistory = undoableAction(msg('Move feature visibility rule'), moveFeatureVisibilityRuleDown);
  const moveFeatureVisibilityRuleUpWithHistory = undoableAction(msg('Move feature visibility rule'), moveFeatureVisibilityRuleUp);
  const removeFeatureVisibilityRuleWithHistory = undoableAction(msg('Remove feature visibility rule'), removeFeatureVisibilityRule);
  const setFeatureVisibilityRuleFieldWithHistory = undoableAction(
    msg('Edit feature visibility rule'),
    setFeatureVisibilityRuleField
  );
  const setFeatureVisibilityWithHistory = undoableAction(msg('Change feature visibility'), setFeatureVisibility);
  const updateClickedFeatureVisibilityWithHistory = undoableAction(
    msg('Change feature visibility'),
    updateClickedFeatureVisibility
  );
  const handleFeatureVisibilityScopeChoiceWithHistory = undoableAction(
    msg('Change feature visibility'),
    handleFeatureVisibilityScopeChoice
  );
  const requestFeatureColorChangeWithHistory = undoableAction(msg('Change feature color'), requestFeatureColorChange);
  const setFeatureColorValueWithHistory = undoableAction(msg('Change feature color'), setFeatureColorValue);
  const updateClickedFeatureColorWithHistory = undoableAction(msg('Change feature color'), updateClickedFeatureColor);
  const editLegendEntryColor = editEditorIntent(msg('Change legend color'), updateLegendEntryColor, { domains: LEGEND_FILL_DOMAINS });
  /** @param {number} index @param {string} color */
  const updateLegendEntryColorShown = (index, color) => (
    legendRowHasRules(index) ? updateLegendEntryColor(index, color) : editLegendEntryColor(index, color)
  );
  const setLegendRowDefaultColorWithHistory = undoableAction(msg('Change feature color'), setLegendRowDefaultColor);
  // OV-1281: the Legend item popup's Color on a row a type's default color draws sets that color (OIC-029).
  // The template hands the commit to its `v-committed-value` field, which shows the row's color once it has
  // settled: a refused edit leaves the view as it was.
  /** @param {number} index @param {string} color */
  const updateLegendEntryColorInPopup = async (index, color) => {
    const entry = state.activeDrawing().legendEntries.value[index];
    const type = entry ? await preparedLegendRowDefaultColorType(legendEntryKey(entry)) : null;
    if (entry && typeof type === 'string') return setLegendRowDefaultColorWithHistory(legendEntryKey(entry), type, color);
    // null: no default color draws the row (the swatch only); undefined: the matches went stale, so do nothing.
    return type === null ? updateLegendEntryColorShown(index, color) : type;
  };
  const setLegendEntryStrokeColorValueWithHistory = editEditorIntent(
    msg('Change legend stroke color'),
    setLegendEntryStrokeColorValue
  );
  // OV-161, D-12: a dialog choice is one History step, Cancel records none, and
  // neither starts while History's step is open (`createDialogChoice`).
  /**
   * @param {() => string} label
   * @param {(choice: string, ...rest: unknown[]) => unknown} handler
   * @param {() => void} [cancel]
   */
  const scopeChoiceWithHistory = (label, handler, cancel = cancelFeatureStyleScope) => (
    withHistory(label, handler, cancel)
  );
  const handleColorScopeChoiceWithHistory = scopeChoiceWithHistory(() => msg('Change feature color'), handleColorScopeChoice);
  const handleFeatureColorScopeChoiceWithHistory = scopeChoiceWithHistory(() => msg('Change feature color'), handleFeatureStyleScopeChoice);
  const handleFeatureStrokeScopeChoiceWithHistory = scopeChoiceWithHistory(
    () => msg('Change feature stroke'),
    showingEditorIntent(handleFeatureStyleScopeChoice, STROKE_DOMAINS)
  );
  /** @param {string} choice @param {unknown[]} rest */
  const handleFeatureStyleScopeChoiceWithHistory = (choice, ...rest) => (
    featureStyleScopeDialog.kind === 'stroke'
      ? handleFeatureStrokeScopeChoiceWithHistory(choice, ...rest)
      : handleFeatureColorScopeChoiceWithHistory(choice, ...rest)
  );
  // The label scopes the Color/Stroke Change Scope dialog offers; a
  // label that repeats the Legend name, or the other label, is not offered.
  // OV-648: the color owner decides which label scopes select other features.
  const featureStyleScopeLabelOptions = computed(() => ({
    displayLabel: featureStyleScopeDialog.displayLabelOffered === true,
    annotationLabel: featureStyleScopeDialog.annotationLabelOffered === true
  }));
  return {
    undoableAction, featureIdentityNoticeSummary, removeUnmatchedFeatureEdits, addFeatureVisibilityRuleWithHistory,
    moveFeatureVisibilityRuleDownWithHistory, moveFeatureVisibilityRuleUpWithHistory, removeFeatureVisibilityRuleWithHistory,
    setFeatureVisibilityRuleFieldWithHistory, setFeatureVisibilityWithHistory, updateClickedFeatureVisibilityWithHistory,
    handleFeatureVisibilityScopeChoiceWithHistory, requestFeatureColorChangeWithHistory, setFeatureColorValueWithHistory,
    updateClickedFeatureColorWithHistory, updateLegendEntryColorShown, updateLegendEntryColorInPopup,
    setLegendEntryStrokeColorValueWithHistory, scopeChoiceWithHistory, handleColorScopeChoiceWithHistory,
    handleFeatureStyleScopeChoiceWithHistory, featureStyleScopeLabelOptions
  };
};

/**
 * @typedef {object} LegendEditHistoryPorts
 * @property {WebState} state
 * @property {UndoableAction} undoableAction
 * @property {ScopeChoiceWithHistory} scopeChoiceWithHistory
 * @property {(label: string, fn: () => unknown) => unknown} runUndoableCheckpoint
 *   History's step whose rule commits join it, so one Undo returns them all.
 * @property {() => unknown} handleLegendNameCommit
 * @property {(group: ColorRuleLegendGroup | null | undefined) => unknown} joinLegendGroup
 * @property {(from: number, to: number) => unknown} moveLegendEntry
 * @property {(choice: string) => unknown} handleLegendRenameChoice
 * @property {(restore: () => unknown) => unknown} restoreLegendRows Runs a Restore of deleted Legend rows.
 * @property {(indexes?: number[] | null) => unknown} restoreDeletedLegendEntries
 * @property {() => void} cancelLegendRename
 * @property {(choice: string) => unknown} handlePaletteColorsChoice
 * @property {() => void} cancelPaletteColorsDialog
 * @property {(choice: string) => unknown} handleResetColorChoice
 * @property {() => void} cancelResetColor
 */

/**
 * The Legend rename and reorder, and the palette and Reset color dialog
 * choices, as History steps. Registers nothing.
 * @param {LegendEditHistoryPorts} ports
 */
export const createLegendEditHistory = ({
  state, undoableAction, scopeChoiceWithHistory, runUndoableCheckpoint, handleLegendNameCommit, joinLegendGroup,
  moveLegendEntry, handleLegendRenameChoice, restoreLegendRows, restoreDeletedLegendEntries, cancelLegendRename,
  handlePaletteColorsChoice, cancelPaletteColorsDialog, handleResetColorChoice, cancelResetColor
}) => {
  const { legendRenameDialog } = state;
  const handleLegendNameCommitWithHistory = undoableAction(msg('Rename legend item'), handleLegendNameCommit);
  const joinLegendGroupWithHistory = undoableAction(msg('Rename legend item'), joinLegendGroup);
  // The Legend tab's grip reorders a row by drag or by the Up and Down arrow
  // keys; neither is a click, so each move records its own History step. The
  // list keys its rows by identity and keeps the moved row's focus
  // (`v-keeps-list-focus`, OV-3070), so the grip the keys moved stays focused.
  const moveLegendEntryWithHistory = undoableAction(msg('Move legend item'), moveLegendEntry);
  const legendReorder = reactive({ from: /** @type {number | null} */ (null), over: /** @type {number | null} */ (null) });
  /** @param {number} from @param {number} to */
  const reorderLegendEntry = (from, to) => {
    const count = state.activeDrawing().legendEntries.value.length;
    if (from === to || to < 0 || to >= count) return;
    moveLegendEntryWithHistory(from, to);
  };
  /** @param {KeyboardEvent} event @param {number} index */
  const handleLegendReorderKey = (event, index) => {
    const step = event.key === 'ArrowUp' ? -1 : (event.key === 'ArrowDown' ? 1 : 0);
    if (!step) return;
    event.preventDefault();
    reorderLegendEntry(index, index + step);
  };
  /** @param {DragEvent} event @param {number} index */
  const startLegendReorder = (event, index) => {
    legendReorder.from = index;
    legendReorder.over = index;
    const row = event.target instanceof Element ? event.target.closest('[data-legend-drawer-row]') : null;
    if (event.dataTransfer) {
      event.dataTransfer.effectAllowed = 'move';
      event.dataTransfer.setData('text/plain', String(index));
      if (row) event.dataTransfer.setDragImage(row, 12, 12);
    }
  };
  /** @param {DragEvent} event @param {number} index */
  const overLegendReorder = (event, index) => {
    if (legendReorder.from === null) return;
    event.preventDefault();
    legendReorder.over = index;
  };
  /** @param {DragEvent} event @param {number} index */
  const dropLegendReorder = (event, index) => {
    const from = legendReorder.from;
    if (from === null) return;
    event.preventDefault();
    legendReorder.from = null;
    legendReorder.over = null;
    reorderLegendEntry(from, index);
  };
  const endLegendReorder = () => {
    legendReorder.from = null;
    legendReorder.over = null;
  };
  // R15-3 (OV-285): the Merge of a rename onto a deleted row's caption is the
  // Restore of that row (by its Python key), then the merge, in one History
  // checkpoint: the merge's rule commit joins it, so one Undo returns both.
  // Keep both onto a deleted row in the renamed row's color restores it the
  // same way: the two are one row (D-28).
  /** @param {string} choice */
  const chooseLegendRename = (choice) => {
    const key = legendRenameDialog.deletedTargetKey;
    const restores = choice === 'merge' || (choice === 'keep-both' && legendRenameDialog.keepBothRestores);
    if (!restores || !key) return handleLegendRenameChoice(choice);
    return runUndoableCheckpoint(msg('Rename legend item'), async () => {
      const deleted = state.activeDrawing().deletedLegendEntries.value;
      const index = deleted.findIndex((/** @type {KeptLegendRow} */ entry) => (entry.originalCaption || entry.caption) === key);
      if (index >= 0) restoreLegendRows(() => restoreDeletedLegendEntries([index]));
      return handleLegendRenameChoice(choice);
    });
  };
  const handleLegendRenameChoiceWithHistory = scopeChoiceWithHistory(
    () => msg('Rename legend item'),
    chooseLegendRename,
    cancelLegendRename
  );
  // D-15: the palette dialog's choice is one History step too.
  const handlePaletteColorsChoiceWithHistory = scopeChoiceWithHistory(
    () => msg('Change setting'),
    handlePaletteColorsChoice,
    cancelPaletteColorsDialog
  );
  const handleResetColorChoiceWithHistory = scopeChoiceWithHistory(
    () => msg('Reset feature color'),
    handleResetColorChoice,
    cancelResetColor
  );
  return {
    handleLegendNameCommitWithHistory, joinLegendGroupWithHistory, legendReorder, handleLegendReorderKey,
    startLegendReorder, overLegendReorder, dropLegendReorder, endLegendReorder, handleLegendRenameChoiceWithHistory,
    handlePaletteColorsChoiceWithHistory, handleResetColorChoiceWithHistory
  };
};

/**
 * A History command: what it applies and reverts.
 * @typedef {{ label: string, apply: () => Promise<boolean>, revert: () => Promise<boolean>, estimateBytes: () => number }} SelectedFeaturesVisibilityCommand
 */

/**
 * @typedef {object} SelectedFeatureHistoryPorts
 * @property {WebState} state
 * @property {UndoableAction} undoableAction
 * @property {EditEditorIntent} editEditorIntent
 * @property {RunUndoable} runUndoable
 * @property {(label: string, buildCommand: () => SelectedFeaturesVisibilityCommand | null) => Promise<unknown>} runUndoableCommand
 *   History's step of a command the builder returns; no command records no step.
 * @property {() => unknown} resetClickedFeatureFillColor
 * @property {(strokeColor: string | null | undefined, strokeWidth: unknown) => unknown} updateClickedFeatureStroke
 * @property {(value: string | null) => unknown} setClickedFeatureStrokeColorValue
 * @property {(value: unknown) => unknown} setClickedFeatureStrokeWidthValue
 * @property {() => unknown} resetClickedFeatureStroke
 * @property {(feature: Record<string, unknown>, color: string, customCaption?: LegendRowName | null) => unknown} setFeatureColor
 * @property {(features: unknown, color: string, caption: unknown) => unknown} applyColorToSelectedFeatures
 * @property {(svg?: SVGSVGElement | null) => void} syncFeatureSelectionClasses
 * @property {(features: VisibilityTarget[], modeRaw: unknown) => SelectedFeaturesVisibilityCommand | null} buildSelectedFeaturesVisibilityCommand
 * @property {(options?: { keepAnchorId?: unknown, clearStatus?: boolean, syncDom?: boolean }) => void} clearFeatureSelection
 * @property {(features: unknown, strokeColor: string, strokeWidth: unknown) => unknown} applyStrokeToSelectedFeatures
 * @property {(feat: Record<string, unknown>, eventLike?: PopupOrigin | null) => unknown} openFeatureEditorForFeature
 * @property {() => unknown} updateClickedFeatureLabelText
 * @property {(choice: string) => unknown} handleLabelTextScopeChoice
 * @property {(choice: string) => unknown} handleHiddenLabelTextChoice
 * @property {(featureId: unknown, nextTextRaw: unknown) => unknown} requestLabelTextChangeByFeatureId
 * @property {(labelKey: string | null | undefined, nextTextRaw: unknown) => unknown} requestLabelTextChangeByKey
 * @property {(featureId: unknown) => EditableLabelEntry | null} getEditableLabelByFeatureId
 * @property {() => unknown} resetAllLabelTextOverrides
 * @property {(event: { target: LabelTableInput | null, sourceInput?: LabelTableInput }) => unknown} loadLabelOverrideTable
 */

/**
 * The clicked feature's fill reset and stroke edits, the selected features'
 * bulk edits, and the label text edits as History steps. Registers nothing.
 * @param {SelectedFeatureHistoryPorts} ports
 */
export const createSelectedFeatureHistory = ({
  state, undoableAction, editEditorIntent, runUndoable, runUndoableCommand, resetClickedFeatureFillColor,
  updateClickedFeatureStroke, setClickedFeatureStrokeColorValue, setClickedFeatureStrokeWidthValue,
  resetClickedFeatureStroke, setFeatureColor, applyColorToSelectedFeatures, syncFeatureSelectionClasses,
  buildSelectedFeaturesVisibilityCommand, clearFeatureSelection, applyStrokeToSelectedFeatures, openFeatureEditorForFeature,
  updateClickedFeatureLabelText, handleLabelTextScopeChoice, handleHiddenLabelTextChoice, requestLabelTextChangeByFeatureId,
  requestLabelTextChangeByKey, getEditableLabelByFeatureId, resetAllLabelTextOverrides, loadLabelOverrideTable
}) => {
  const { selectedFeatures } = state;
  const resetClickedFeatureFillColorWithHistory = undoableAction(msg('Reset feature color'), resetClickedFeatureFillColor);
  const updateClickedFeatureStrokeWithHistory = editEditorIntent(msg('Change feature stroke'), updateClickedFeatureStroke);
  const setClickedFeatureStrokeColorValueWithHistory = editEditorIntent(
    msg('Change feature stroke'),
    setClickedFeatureStrokeColorValue
  );
  const setClickedFeatureStrokeWidthValueWithHistory = editEditorIntent(
    msg('Change feature stroke'),
    setClickedFeatureStrokeWidthValue
  );
  const resetClickedFeatureStrokeWithHistory = editEditorIntent(msg('Reset feature stroke'), resetClickedFeatureStroke);
  const setFeatureColorWithHistory = undoableAction(msg('Change feature color'), setFeatureColor);
  const selectedFeatureBulkColor = ref('#2563eb');
  const selectedFeatureBulkCaption = ref('Selected features');
  const selectedFeatureBulkVisibility = ref('off');
  const selectedFeatureBulkStrokeColor = ref('#1f2937');
  const selectedFeatureBulkStrokeWidth = ref(1.5);
  const applySelectedFeatureColor = () => runUndoable(msg('Change selected feature color'), async () => {
    const changed = await applyColorToSelectedFeatures(
      selectedFeatures.value,
      selectedFeatureBulkColor.value,
      selectedFeatureBulkCaption.value
    );
    if (changed) syncFeatureSelectionClasses();
    return changed;
  });
  const applySelectedFeatureVisibility = async () => {
    const changed = await runUndoableCommand(msg('Change selected feature visibility'), () =>
      buildSelectedFeaturesVisibilityCommand(selectedFeatures.value, selectedFeatureBulkVisibility.value)
    );
    if (changed) clearFeatureSelection({ clearStatus: true });
    return changed;
  };
  const editSelectedFeatureStroke = editEditorIntent(msg('Change selected feature stroke'), () => applyStrokeToSelectedFeatures(
    selectedFeatures.value,
    selectedFeatureBulkStrokeColor.value,
    selectedFeatureBulkStrokeWidth.value
  ));
  const applySelectedFeatureStroke = async () => {
    const changed = await editSelectedFeatureStroke();
    if (changed === true) syncFeatureSelectionClasses();
    return changed;
  };
  const openFirstSelectedFeature = (event = null) => {
    const first = selectedFeatures.value[0] || null;
    if (!first) return null;
    return openFeatureEditorForFeature(first, event);
  };
  const updateClickedFeatureLabelTextWithHistory = undoableAction(msg('Change label text'), updateClickedFeatureLabelText);
  const handleLabelTextScopeChoiceWithHistory = undoableAction(msg('Change label text'), handleLabelTextScopeChoice);
  // The Label Not Shown choice answers the label edit whose step is open: one step, under that step's name.
  // A registered test reads the label from this line (feature-label-visual-unit.test.mjs).
  const handleHiddenLabelTextChoiceWithHistory = undoableAction(msg('Change label visibility'), handleHiddenLabelTextChoice, { joinsOpenStep: true });
  const requestLabelTextChangeByFeatureIdWithHistory = undoableAction(
    msg('Change label text'),
    requestLabelTextChangeByFeatureId
  );
  const requestLabelTextChangeByKeyWithHistory = undoableAction(msg('Change label text'), requestLabelTextChangeByKey);
  // The Features list's label input commits what was typed; its draft then follows the label.
  /** @param {string} featureId @param {string} text */
  const commitListLabelDraft = (featureId, text) => {
    const committed = requestLabelTextChangeByFeatureIdWithHistory(featureId, text);
    const entry = getEditableLabelByFeatureId(featureId);
    if (entry) entry.draftText = entry.text;
    return committed;
  };
  const resetAllLabelTextOverridesWithHistory = undoableAction(msg('Reset label edits'), resetAllLabelTextOverrides);
  const loadLabelOverrideTableWithHistory = undoableAction(msg('Load label edits'), loadLabelOverrideTable);
  return {
    resetClickedFeatureFillColorWithHistory, updateClickedFeatureStrokeWithHistory,
    setClickedFeatureStrokeColorValueWithHistory, setClickedFeatureStrokeWidthValueWithHistory,
    resetClickedFeatureStrokeWithHistory, setFeatureColorWithHistory, selectedFeatureBulkColor, selectedFeatureBulkCaption,
    selectedFeatureBulkVisibility, selectedFeatureBulkStrokeColor, selectedFeatureBulkStrokeWidth, applySelectedFeatureColor,
    applySelectedFeatureVisibility, applySelectedFeatureStroke, openFirstSelectedFeature,
    updateClickedFeatureLabelTextWithHistory, handleLabelTextScopeChoiceWithHistory, handleHiddenLabelTextChoiceWithHistory,
    requestLabelTextChangeByFeatureIdWithHistory, requestLabelTextChangeByKeyWithHistory, commitListLabelDraft,
    resetAllLabelTextOverridesWithHistory, loadLabelOverrideTableWithHistory
  };
};
