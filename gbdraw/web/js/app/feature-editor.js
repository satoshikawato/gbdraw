import { createFeatureLabelActions } from './feature-editor/label-actions.js';
import { createFeatureRuleActions } from './feature-editor/rule-actions.js';
import { createFeatureSvgActions } from './feature-editor/svg-actions.js';
import { createFeatureVisibilityActions } from './feature-editor/visibility-actions.js';
import { createFeatureFillActions } from './feature-editor/fill-actions.js';
import { createFeatureStrokeActions } from './feature-editor/stroke-actions.js';
import { createFeatureBulkStyleActions } from './feature-editor/bulk-style-actions.js';
import { createFeatureLabelStyleActions } from './feature-editor/label-style-actions.js';

export const createFeatureEditor = ({
  state,
  nextTick,
  legendActions,
  svgActions,
  featureSelection = null,
  previewRuntime = null,
  history = null
}) => {
  let bulkStyleOwner = null;
  const bulkStyleActions = {
    requestFeatureBulkStyleChange: (...args) => {
      if (!bulkStyleOwner) throw new Error('Bulk Feature style actions are unavailable.');
      return bulkStyleOwner.requestFeatureBulkStyleChange(...args);
    }
  };
  const ruleActions = createFeatureRuleActions({
    state,
    nextTick,
    legendActions,
    bulkStyleActions
  });
  const labelActions = createFeatureLabelActions({ state, previewRuntime });
  const featureSvgActions = createFeatureSvgActions({
    state,
    getFeatureColor: ruleActions.getFeatureColor,
    getEffectiveLegendCaption: ruleActions.getEffectiveLegendCaption,
    onFeaturePopupOpened: labelActions.syncLabelEditor,
    featureSelection,
    previewRuntime
  });
  const fillActions = createFeatureFillActions({
    state,
    history,
    previewRuntime,
    featureSvgActions,
    legendActions
  });
  bulkStyleOwner = createFeatureBulkStyleActions({
    state,
    history,
    previewRuntime,
    featureSvgActions,
    legendActions
  });
  const strokeActions = createFeatureStrokeActions({
    state,
    history,
    previewRuntime,
    featureSvgActions,
    legendActions,
    pendingFeatureStrokePlan: state.pendingFeatureStrokePlan,
    featureStrokePlanStatus: state.featureStrokePlanStatus,
    featureStrokePlanProgress: state.featureStrokePlanProgress,
    getEffectiveLegendCaption: ruleActions.getEffectiveLegendCaption
  });
  const labelStyleActions = createFeatureLabelStyleActions({
    state,
    history,
    previewRuntime,
    featureSvgActions,
    legendActions,
    pendingFeatureLabelPlan: state.pendingFeatureLabelPlan,
    featureLabelPlanStatus: state.featureLabelPlanStatus,
    featureLabelPlanProgress: state.featureLabelPlanProgress
  });
  const requestedFillValue = (feat) => {
    const color = ruleActions.getFeatureFillViewModel(feat)?.effectiveColor || '#cccccc';
    return color === 'none' ? { kind: 'none' } : { kind: 'color', color };
  };
  const handleLegendNameCommit = () => {
    const clicked = state.clickedFeature?.value;
    if (!clicked?.feat) return false;
    return fillActions.requestFeatureFillChange(
      clicked.feat,
      requestedFillValue(clicked.feat),
      clicked.legendName,
      { source: 'popup', closePopupOnDialog: true }
    );
  };
  const selectLegendNameOption = (caption) => {
    if (!state.clickedFeature?.value) return false;
    state.clickedFeature.value.legendName = String(caption || '').trim();
    return handleLegendNameCommit();
  };
  const visibilityActions = createFeatureVisibilityActions({
    state,
    featureSvgActions,
    previewRuntime
  });
  const openFeatureEditorForFeature = (feat, eventLike = null) => {
    return featureSvgActions.openFeatureEditorForFeature(feat, eventLike);
  };

  return {
    addPriorityRule: ruleActions.addPriorityRule,
    addFeature: ruleActions.addFeature,
    removeFeature: ruleActions.removeFeature,
    getFeatureShape: ruleActions.getFeatureShape,
    setFeatureShape: ruleActions.setFeatureShape,
    addSpecificRule: ruleActions.addSpecificRule,
    applySpecificRulePreset: ruleActions.applySpecificRulePreset,
    clearAllSpecificRules: ruleActions.clearAllSpecificRules,
    downloadSpecificRulesTsv: ruleActions.downloadSpecificRulesTsv,
    canMoveSpecificRule: ruleActions.canMoveSpecificRule,
    isOpaqueSpecificRule: ruleActions.isOpaqueSpecificRule,
    moveSpecificRuleDown: ruleActions.moveSpecificRuleDown,
    moveSpecificRuleUp: ruleActions.moveSpecificRuleUp,
    removeSpecificRule: ruleActions.removeSpecificRule,
    setSpecificRuleField: ruleActions.setSpecificRuleField,
    requestFeatureBulkStyleChange: bulkStyleActions.requestFeatureBulkStyleChange,
    getFeatureColor: ruleActions.getFeatureColor,
    getFeatureColorValue: ruleActions.getFeatureColorValue,
    getFeatureFillViewModel: ruleActions.getFeatureFillViewModel,
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
    requestFeatureFillChange: fillActions.requestFeatureFillChange,
    resolveFeatureFillScope: fillActions.resolveFeatureFillScope,
    cancelFeatureFillScope: fillActions.cancelFeatureFillScope,
    requestSelectedFeatureFillChange: fillActions.requestSelectedFeatureFillChange,
    handleLegendNameCommit,
    selectLegendNameOption,
    getFeatureStrokeViewModel: strokeActions.getFeatureStrokeViewModel,
    requestFeatureStrokeChange: strokeActions.requestFeatureStrokeChange,
    resolveFeatureStrokeScope: strokeActions.resolveFeatureStrokeScope,
    cancelFeatureStrokeScope: strokeActions.cancelFeatureStrokeScope,
    requestSelectedFeatureStrokeChange: strokeActions.requestSelectedFeatureStrokeChange,
    attachSvgFeatureHandlers: featureSvgActions.attachSvgFeatureHandlers,
    openFeatureEditorForFeature,
    getEditableLabelByFeatureId: labelActions.getEditableLabelByFeatureId,
    getFeatureLabelViewModel: labelStyleActions.getFeatureLabelViewModel,
    syncLabelEditor: labelActions.syncLabelEditor,
    downloadLabelOverrideTable: labelActions.downloadLabelOverrideTable,
    loadLabelOverrideTable: labelActions.loadLabelOverrideTable,
    requestFeatureLabelTextChange: labelStyleActions.requestFeatureLabelTextChange,
    requestSelectedFeatureLabelTextChange: labelStyleActions.requestSelectedFeatureLabelTextChange,
    resolveFeatureLabelScope: labelStyleActions.resolveFeatureLabelScope,
    cancelFeatureLabelScope: labelStyleActions.cancelFeatureLabelScope,
    handleGlobalLabelModeChoice: labelActions.handleGlobalLabelModeChoice,
    updateClickedFeatureLabelVisibility: labelActions.updateClickedFeatureLabelVisibility,
    reconcileLabelOverrides: labelActions.reconcileLabelOverrides,
    resetAllLabelTextOverrides: labelActions.resetAllLabelTextOverrides
  };
};
