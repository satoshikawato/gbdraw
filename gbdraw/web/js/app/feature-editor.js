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
  history,
  getCommittedRequest,
  getCommittedSession = () => null,
  readResourceRecordCount = null,
  readFeatureOverrideTable = null,
  isCurrentFeature,
  nextTick,
  legendActions,
  svgActions,
  featureSelection = null,
  previewRuntime = null,
  isPatternEditAvailable = () => true,
  previewTransformInteraction = null
}) => {
  const { ref, computed, watch, reactive } = window.Vue;
  const ruleActions = createFeatureRuleActions({ state, nextTick, legendActions, rulePreparation, history, svgActions, ref, computed, isPatternEditAvailable });
  // Show feature and label (Owner Q2) sets Feature visibility through its owner.
  const labelActions = createFeatureLabelActions({
    state, previewRuntime, rulePreparation, ref, computed, watch, nextTick, getCommittedRequest,
    setFeatureVisibility: (...args) => visibilityActions.setFeatureVisibility(...args)
  });
  const featureSvgActions = createFeatureSvgActions({
    state,
    getFeatureColor: ruleActions.getFeatureColor,
    getEffectiveLegendCaption: ruleActions.getEffectiveLegendCaption,
    rulePreparation,
    onFeaturePopupOpened: labelActions.syncLabelEditor,
    featureSelection,
    previewRuntime,
    previewTransformInteraction
  });
  const colorActions = createFeatureColorActions({
    state,
    rulePreparation,
    nextTick,
    legendActions,
    svgActions,
    ruleActions,
    featureSvgActions,
    previewRuntime
  });
  const visibilityActions = createFeatureVisibilityActions({
    state,
    featureSvgActions,
    labelActions,
    previewRuntime
  });
  // A loaded table shows on the displayed Result as a History apply does (R3).
  const featureEditTableActions = createFeatureEditTableActions({
    state, ref, computed, getCommittedSession, readResourceRecordCount, readFeatureOverrideTable,
    projectFeatureEdits: () => {
      visibilityActions.reconcileFeatureVisibility();
      labelActions.reconcileLabelOverrides();
      labelActions.applyFeatureVisibilityToLabels();
    }
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
    placementActions: createFeaturePlacementActions({ state, history, getCommittedRequest, isCurrentFeature, reactive, nextTick }),
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
    getFeatureVisibility: visibilityActions.getFeatureVisibility,
    reconcileFeatureVisibility: visibilityActions.reconcileFeatureVisibility,
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
    requestLabelTextChangeByFeatureId: labelActions.requestLabelTextChangeByFeatureId,
    requestLabelTextChangeByKey: labelActions.requestLabelTextChangeByKey,
    reconcileLabelOverrides: labelActions.reconcileLabelOverrides,
    resetAllLabelTextOverrides: labelActions.resetAllLabelTextOverrides
  };
};
