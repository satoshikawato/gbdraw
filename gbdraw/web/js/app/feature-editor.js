// @ts-check
/** @import { RulePreparation } from './rule-matching.js' */
/** @import { PreparedFileLegend, RuleActionsPorts } from './feature-editor/rule-actions.js' */
/** @import { VisibilityActionsPorts } from './feature-editor/visibility-actions.js' */
/** @import { FeatureSelectionPort, PrepareDrawnFeatureMatchesPort, PreviewTransformInteractionPort } from './feature-editor/svg-actions.js' */
/** @import { ReadUserDefaultColorPort, SetDefaultColorPort } from './feature-editor/color-actions.js' */
import { createFeatureColorActions } from './feature-editor/color-actions.js';
import { createFeatureLabelActions } from './feature-editor/label-actions.js';
import { createFeatureRuleActions } from './feature-editor/rule-actions.js';
import { createFeatureSvgActions } from './feature-editor/svg-actions.js';
import { createFeatureVisibilityActions } from './feature-editor/visibility-actions.js';
import { createFeaturePlacementActions } from './feature-editor/placement-actions.js';
import { createFeatureEditTableActions } from './feature-editor/feature-edit-table.js';

/**
 * The reactions registered in `editorPorts` once their owner exists (the label
 * owner's, and the visibility owner's preparation); the owners created before
 * it call them through the bag only.
 * @typedef {RuleActionsPorts & VisibilityActionsPorts & {
 *   syncLabelEditor: (options?: Record<string, any>) => void,
 *   prepareDrawnFeatureMatches: PrepareDrawnFeatureMatchesPort
 * }} FeatureEditorPorts
 */

/**
 * @typedef {object} FeatureEditorOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {RulePreparation} rulePreparation
 * @property {(label: string, fn: () => any, options?: Record<string, any>) => any} runUndoable History's undoable step.
 * @property {(label: string, fn: () => any, options?: Record<string, any>) => any} runUndoableCheckpoint
 *   History's undoable step that stores a checkpoint of the Result.
 * @property {(close: () => unknown) => void} [closeAfterDialogChoice]
 *   Closes a choice dialog once the History step of its choice ends (D-12).
 * @property {() => ({ mode?: string, diagramOptions?: Record<string, any> } | null)} getCommittedRequest
 *   The committed canonical request (Python owns the option fields, R7).
 * @property {() => Record<string, any> | null} [getCommittedSession] The committed canonical Session.
 * @property {((resourceId: string, kind: string) => any) | null} [readResourceRecordCount]
 *   Counts the records of a committed resource.
 * @property {((payload: Record<string, any>) => Promise<{ result?: any }>) | null} [readFeatureOverrideTable]
 *   The diagram helper that reads a Feature Edits TSV (R7).
 * @property {(feature: Record<string, any>) => boolean} isCurrentFeature Whether the feature belongs to the displayed Result.
 * @property {(callback?: () => void) => Promise<void>} nextTick Vue `nextTick`
 * @property {(intents: Record<string, any>[], options?: { previousFileIntents?: Record<string, any>[], isCurrent?: () => boolean }) => Promise<PreparedFileLegend | false>} prepareFileLegendEntries
 *   The Legend owner's preparation of the rows the rules draw.
 * @property {(options?: { replaceGeneratedInventory?: boolean }) => any} extractLegendEntries
 *   The Legend owner's reading of the mounted Legend rows.
 * @property {() => void} onLegendGeometryChanged The Legend owner's reaction to a change of Legend geometry.
 * @property {FeatureSelectionPort | null} [featureSelection]
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 * @property {((changes: { featureId: string, mode: string }[], options?: { reason?: string }) => boolean) | null} [applyFeatureVisibilityChanges]
 *   The preview owner's projection of feature visibility changes.
 * @property {(index: number) => any} selectResult The preview owner's selection of a Result.
 * @property {() => boolean} [isPatternEditAvailable] False while a Session import is pending.
 * @property {PreviewTransformInteractionPort | null} [previewTransformInteraction]
 * @property {(options?: { recolor?: Record<string, any>, prepareRules?: boolean }) => boolean | Promise<boolean>} projectPaletteAndRules
 *   The root's projection of the palette and the rules (R3).
 * @property {() => any} projectFeatureEdits
 *   The root's projection of loaded feature edits onto the displayed Result (R3).
 * @property {(payload: Record<string, any>, options?: Record<string, any>) => Promise<Record<string, any>>} evaluateLabelRules
 *   The root's Python evaluation of Label TSV rows against the displayed labels (R7). * @property {ReadUserDefaultColorPort} readUserDefaultColor
 * @property {SetDefaultColorPort} setDefaultColor
 */

/** @param {FeatureEditorOptions} options */
export const createFeatureEditor = ({
  state,
  rulePreparation,
  // R13: History's undoable steps and the preview owner's commit, visibility
  // projection, and Result selection arrive as ports.
  runUndoable,
  runUndoableCheckpoint,
  closeAfterDialogChoice,
  getCommittedRequest,
  getCommittedSession = () => null,
  readResourceRecordCount = null,
  readFeatureOverrideTable = null,
  isCurrentFeature,
  nextTick,
  prepareFileLegendEntries,
  extractLegendEntries,
  onLegendGeometryChanged,
  featureSelection = null,
  commitActiveResultEdit = null,
  applyFeatureVisibilityChanges = null,
  selectResult,
  isPatternEditAvailable = () => true,
  previewTransformInteraction = null,
  projectPaletteAndRules,
  projectFeatureEdits,
  evaluateLabelRules,
  readUserDefaultColor,
  setDefaultColor
}) => {
  const { ref, computed, watch, reactive } = window.Vue;
  // R13: the label owner's reactions, registered once it exists; the owners
  // created before it call them through this object only.
  const editorPorts = /** @type {FeatureEditorPorts} */ ({});
  const ruleActions = createFeatureRuleActions({
    state, prepareFileLegendEntries, rulePreparation, runUndoable, runUndoableCheckpoint, projectPaletteAndRules,
    ports: editorPorts, getCommittedRequest, ref, computed, watch, isPatternEditAvailable
  });
  const featureSvgActions = createFeatureSvgActions({
    state,
    getFeatureColor: ruleActions.getFeatureColor,
    getEffectiveLegendCaption: ruleActions.getEffectiveLegendCaption,
    // R13: the visibility owner prepares the matches the popup reads; it is
    // created after this owner, so the port registers in `editorPorts`.
    prepareDrawnFeatureMatches: (...args) => editorPorts.prepareDrawnFeatureMatches(...args),
    onFeaturePopupOpened: (...args) => editorPorts.syncLabelEditor(...args),
    featureSelection,
    applyFeatureVisibilityChanges,
    previewTransformInteraction
  });
  const colorActions = createFeatureColorActions({
    state,
    extractLegendEntries,
    onLegendGeometryChanged,
    ruleActions,
    getFeatureElements: featureSvgActions.getFeatureElements,
    getFeatureFillElements: featureSvgActions.getFeatureFillElements,
    commitActiveResultEdit,
    closeAfterDialogChoice,
    readUserDefaultColor,
    setDefaultColor
  });
  // The visibility owner comes before the label owner, which receives its
  // transition as a port: Show feature and label (Owner Q2) sets Feature
  // visibility through it.
  const visibilityActions = createFeatureVisibilityActions({
    state,
    applyVisibilityPreviewChanges: featureSvgActions.applyVisibilityPreviewChanges,
    ports: editorPorts,
    selectResult,
    rulePreparation,
    getCommittedRequest
  });
  const labelActions = createFeatureLabelActions({
    state, commitActiveResultEdit, rulePreparation, ref, computed, watch, nextTick, getCommittedRequest,
    setFeatureVisibility: visibilityActions.setFeatureVisibility,
    prepareDrawnFeatureMatches: visibilityActions.prepareDrawnFeatureMatches,
    evaluateLabelRules
  });
  editorPorts.prepareDrawnFeatureMatches = visibilityActions.prepareDrawnFeatureMatches;
  editorPorts.applyFeatureVisibilityToLabels = labelActions.applyFeatureVisibilityToLabels;
  editorPorts.requestAutomaticRerender = labelActions.requestAutomaticRerender;
  editorPorts.syncLabelEditor = labelActions.syncLabelEditor;
  // A loaded table shows on the displayed Result as a History apply does (R3):
  // the root's projection (`projectFeatureEdits`) projects it.
  const featureEditTableActions = createFeatureEditTableActions({
    state, ref, computed, getCommittedSession, readResourceRecordCount, readFeatureOverrideTable,
    projectFeatureEdits
  });
  /** @param {Record<string, any>} feat
   *  @param {{ clientX: number, clientY: number } | null} [eventLike] */
  const openFeatureEditorForFeature = (feat, eventLike = null) => {
    return featureSvgActions.openFeatureEditorForFeature(feat, eventLike);
  };

  return {
    specificRulePattern: ruleActions.specificRulePattern,
    specificRulePatternDraft: ruleActions.specificRulePatternDraft,
    specificRulePatternFieldId: ruleActions.specificRulePatternFieldId,
    editSpecificRulePattern: ruleActions.editSpecificRulePattern,
    retrySpecificRulePattern: ruleActions.retrySpecificRulePattern,
    revertSpecificRulePattern: ruleActions.revertSpecificRulePattern,
    suspendSpecificRulePatternDrafts: ruleActions.suspendSpecificRulePatternDrafts,
    clearSpecificRulePatternDrafts: ruleActions.clearSpecificRulePatternDrafts,
    captureSpecificRulePatternDrafts: ruleActions.captureSpecificRulePatternDrafts,
    restoreSpecificRulePatternDrafts: ruleActions.restoreSpecificRulePatternDrafts,
    placementActions: createFeaturePlacementActions({ state, runUndoable, getCommittedRequest, isCurrentFeature, reactive, nextTick }),
    canRetrySpecificRuleFailure: ruleActions.canRetrySpecificRuleFailure,
    canEditSpecificRuleFailure: ruleActions.canEditSpecificRuleFailure,
    retrySpecificRuleFailure: ruleActions.retrySpecificRuleFailure,
    editSpecificRuleFailure: ruleActions.editSpecificRuleFailure,
    addCustomColor: ruleActions.addCustomColor,
    addPriorityRule: ruleActions.addPriorityRule,
    setLabelFilterMode: ruleActions.setLabelFilterMode,
    addWhitelistRule: ruleActions.addWhitelistRule,
    removeWhitelistRule: ruleActions.removeWhitelistRule,
    removePriorityRule: ruleActions.removePriorityRule,
    addFeature: ruleActions.addFeature,
    removeFeature: ruleActions.removeFeature,
    getFeatureShape: ruleActions.getFeatureShape,
    setFeatureShape: ruleActions.setFeatureShape,
    addSpecificRule: ruleActions.addSpecificRule,
    commitSpecificRules: ruleActions.commitSpecificRules,
    followRestoredSpecificRules: ruleActions.followRestoredRules,
    applySpecificRulePreset: ruleActions.applySpecificRulePreset,
    clearAllSpecificRules: ruleActions.clearAllSpecificRules,
    downloadSpecificRulesTsv: ruleActions.downloadSpecificRulesTsv,
    moveSpecificRuleDown: ruleActions.moveSpecificRuleDown,
    moveSpecificRuleUp: ruleActions.moveSpecificRuleUp,
    removeSpecificRule: ruleActions.removeSpecificRule,
    setSpecificRuleField: ruleActions.setSpecificRuleField,
    getFeatureColor: ruleActions.getFeatureColor,
    getFeatureColorValue: ruleActions.getFeatureColorValue,
    canEditFeatureColor: ruleActions.canEditFeatureColor,
    addFeatureVisibilityRule: visibilityActions.addFeatureVisibilityRule,
    downloadFeatureVisibilityRulesTsv: visibilityActions.downloadFeatureVisibilityRulesTsv,
    featureVisibilityQualifierSuggestions: visibilityActions.featureVisibilityQualifierSuggestions,
    featureVisibilityRuleDetail: visibilityActions.featureVisibilityRuleDetail,
    projectFeatureVisibility: visibilityActions.projectFeatureVisibility,
    handleFeatureVisibilityScopeChoice: visibilityActions.handleFeatureVisibilityScopeChoice,
    moveFeatureVisibilityRuleDown: visibilityActions.moveFeatureVisibilityRuleDown,
    moveFeatureVisibilityRuleUp: visibilityActions.moveFeatureVisibilityRuleUp,
    removeFeatureVisibilityRule: visibilityActions.removeFeatureVisibilityRule,
    setFeatureVisibility: visibilityActions.setFeatureVisibility,
    setSelectedFeaturesVisibility: visibilityActions.setSelectedFeaturesVisibility,
    buildSelectedFeaturesVisibilityCommand: visibilityActions.buildSelectedFeaturesVisibilityCommand,
    setFeatureVisibilityRuleField: visibilityActions.setFeatureVisibilityRuleField,
    updateClickedFeatureVisibility: visibilityActions.updateClickedFeatureVisibility,
    requestFeatureColorChange: colorActions.requestFeatureColorChange,
    setFeatureColorValue: colorActions.setFeatureColorValue,
    updateClickedFeatureColor: colorActions.updateClickedFeatureColor,
    cancelFeatureStyleScope: colorActions.cancelFeatureStyleScope,
    cancelLegendRename: colorActions.cancelLegendRename,
    cancelResetColor: colorActions.cancelResetColor,
    handleColorScopeChoice: colorActions.handleColorScopeChoice,
    handleFeatureStyleScopeChoice: colorActions.handleFeatureStyleScopeChoice,
    handleLegendNameCommit: colorActions.handleLegendNameCommit,
    handleLegendRenameChoice: colorActions.handleLegendRenameChoice,
    selectLegendNameOption: colorActions.selectLegendNameOption,
    renameLegendEntry: colorActions.renameLegendEntry,
    handleResetColorChoice: colorActions.handleResetColorChoice,
    resetClickedFeatureFillColor: colorActions.resetClickedFeatureFillColor,
    getFeatureStrokeColorValue: colorActions.getFeatureStrokeColorValue,
    setClickedFeatureStrokeColorValue: colorActions.setClickedFeatureStrokeColorValue,
    setClickedFeatureStrokeWidthValue: colorActions.setClickedFeatureStrokeWidthValue,
    updateClickedFeatureStroke: colorActions.updateClickedFeatureStroke,
    resetClickedFeatureStroke: colorActions.resetClickedFeatureStroke,
    applyColorToSelectedFeatures: colorActions.applyColorToSelectedFeatures,
    applyStrokeToSelectedFeatures: colorActions.applyStrokeToSelectedFeatures,
    setFeatureColor: colorActions.setFeatureColor,
    attachSvgFeatureHandlers: featureSvgActions.attachSvgFeatureHandlers,
    preparePairwiseInteractionAffordances:
      featureSvgActions.preparePairwiseInteractionAffordances,
    previewAlignmentCandidate: featureSvgActions.previewAlignmentCandidate,
    clearAlignmentCandidatePreview: featureSvgActions.clearAlignmentCandidatePreview,
    showAlignmentOverlay: featureSvgActions.showAlignmentOverlay,
    clearAlignmentOverlay: featureSvgActions.clearAlignmentOverlay,
    dispose: featureSvgActions.dispose,
    openFeatureEditorForFeature,
    refreshFeatureOverrides: ruleActions.refreshFeatureOverrides,
    getEditableLabelByFeatureId: labelActions.getEditableLabelByFeatureId,
    syncLabelEditor: labelActions.syncLabelEditor,
    // R13: label owner ports the root hands Generate and the watchers.
    closeLabelTextScopeDialog: labelActions.closeLabelTextScopeDialog,
    clearLabelBuildNotices: labelActions.clearLabelBuildNotices,
    ...featureEditTableActions,
    downloadLabelOverrideTable: labelActions.downloadLabelOverrideTable,
    loadLabelOverrideTable: labelActions.loadLabelOverrideTable,
    canRetryLabelImportFailure: labelActions.canRetryLabelImportFailure,
    retryLabelImportFailure: labelActions.retryLabelImportFailure,
    editLabelImportFailure: labelActions.editLabelImportFailure,
    updateClickedFeatureLabelText: labelActions.updateClickedFeatureLabelText,
    handleLabelTextScopeChoice: labelActions.handleLabelTextScopeChoice,
    handleHiddenLabelTextChoice: labelActions.handleHiddenLabelTextChoice,
    handleLabelOnChoice: labelActions.handleLabelOnChoice,
    clickedFeatureLabelHint: labelActions.clickedFeatureLabelHint,
    hiddenLabelTextMessage: labelActions.hiddenLabelTextMessage,
    requestLabelTextChangeByFeatureId: labelActions.requestLabelTextChangeByFeatureId,
    requestLabelTextChangeByKey: labelActions.requestLabelTextChangeByKey,
    reconcileLabelOverrides: labelActions.reconcileLabelOverrides,
    resetAllLabelTextOverrides: labelActions.resetAllLabelTextOverrides
  };
};
