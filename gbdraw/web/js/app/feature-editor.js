import { createFeatureColorActions } from './feature-editor/color-actions.js';
import { createFeatureLabelActions } from './feature-editor/label-actions.js';
import { createFeatureRuleActions } from './feature-editor/rule-actions.js';
import { createFeatureSvgActions } from './feature-editor/svg-actions.js';
import { createFeatureVisibilityActions } from './feature-editor/visibility-actions.js';
import { createFeaturePlacementActions } from './feature-editor/placement-actions.js';
import { createFeatureEditTableActions } from './feature-editor/feature-edit-table.js';

export const createFeatureEditor = ({
  state,
  rulePreparation,
  // R13: History's undoable steps and the preview owner's commit, visibility
  // projection, and Result selection arrive as ports.
  runUndoable,
  runUndoableCheckpoint,
  getCommittedRequest,
  getCommittedSession = () => null,
  readResourceRecordCount = null,
  readFeatureOverrideTable = null,
  isCurrentFeature,
  nextTick,
  prepareFileLegendEntries,
  compactLegendEntries,
  extractLegendEntries,
  onLegendGeometryChanged,
  featureSelection = null,
  commitActiveResultEdit = null,
  applyFeatureVisibilityChanges = null,
  selectResult,
  isPatternEditAvailable = () => true,
  previewTransformInteraction = null,
  projectPaletteAndRules,
  projectFeatureEdits
}) => {
  const { ref, computed, watch, reactive } = window.Vue;
  // R13: the label owner's reactions, registered once it exists; the owners
  // created before it call them through this object only.
  const editorPorts = {};
  const ruleActions = createFeatureRuleActions({
    state, nextTick, prepareFileLegendEntries, rulePreparation, runUndoable, runUndoableCheckpoint, projectPaletteAndRules,
    ports: editorPorts, getCommittedRequest, ref, computed, isPatternEditAvailable
  });
  const featureSvgActions = createFeatureSvgActions({
    state,
    getFeatureColor: ruleActions.getFeatureColor,
    getEffectiveLegendCaption: ruleActions.getEffectiveLegendCaption,
    // R13: the rule preparation's runs reach the owners that read its matches
    // as ports; they never hold the preparation.
    runWithDrawnMatches: rulePreparation.runDrawn,
    onFeaturePopupOpened: (...args) => editorPorts.syncLabelEditor(...args),
    featureSelection,
    applyFeatureVisibilityChanges,
    previewTransformInteraction
  });
  const colorActions = createFeatureColorActions({
    state,
    runWithRuleMatches: rulePreparation.run,
    nextTick,
    compactLegendEntries,
    extractLegendEntries,
    onLegendGeometryChanged,
    ruleActions,
    getFeatureElements: featureSvgActions.getFeatureElements,
    getFeatureFillElements: featureSvgActions.getFeatureFillElements,
    commitActiveResultEdit
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
    setFeatureVisibility: visibilityActions.setFeatureVisibility
  });
  editorPorts.applyFeatureVisibilityToLabels = labelActions.applyFeatureVisibilityToLabels;
  editorPorts.requestAutomaticRerender = labelActions.requestAutomaticRerender;
  editorPorts.syncLabelEditor = labelActions.syncLabelEditor;
  // A loaded table shows on the displayed Result as a History apply does (R3):
  // the root's projection (`projectFeatureEdits`) projects it.
  const featureEditTableActions = createFeatureEditTableActions({
    state, ref, computed, getCommittedSession, readResourceRecordCount, readFeatureOverrideTable,
    projectFeatureEdits
  });
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
