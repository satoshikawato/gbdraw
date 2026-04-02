import { createFeatureColorActions } from './feature-editor/color-actions.js';
import { createFeatureLabelActions } from './feature-editor/label-actions.js';
import { createFeatureRuleActions } from './feature-editor/rule-actions.js';
import { createFeatureSvgActions } from './feature-editor/svg-actions.js';
import { createFeatureVisibilityActions } from './feature-editor/visibility-actions.js';

export const createFeatureEditor = ({ state, nextTick, legendActions, svgActions }) => {
  const ruleActions = createFeatureRuleActions({ state, nextTick, legendActions });
  const labelActions = createFeatureLabelActions({ state });
  const featureSvgActions = createFeatureSvgActions({
    state,
    getFeatureColor: ruleActions.getFeatureColor,
    getEffectiveLegendCaption: ruleActions.getEffectiveLegendCaption,
    onFeaturePopupOpened: labelActions.syncLabelEditor
  });
  const colorActions = createFeatureColorActions({
    state,
    nextTick,
    legendActions,
    svgActions,
    ruleActions,
    featureSvgActions
  });
  const visibilityActions = createFeatureVisibilityActions({
    state,
    featureSvgActions
  });
  const openFeatureEditorForFeature = (feat, eventLike = null) => {
    return featureSvgActions.openFeatureEditorForFeature(feat, eventLike);
  };

  return {
    addCustomColor: ruleActions.addCustomColor,
    addPriorityRule: ruleActions.addPriorityRule,
    addFeature: ruleActions.addFeature,
    removeFeature: ruleActions.removeFeature,
    getFeatureShape: ruleActions.getFeatureShape,
    setFeatureShape: ruleActions.setFeatureShape,
    addSpecificRule: ruleActions.addSpecificRule,
    applySpecificRulePreset: ruleActions.applySpecificRulePreset,
    clearAllSpecificRules: ruleActions.clearAllSpecificRules,
    getFeatureColor: ruleActions.getFeatureColor,
    canEditFeatureColor: ruleActions.canEditFeatureColor,
    getFeatureVisibility: visibilityActions.getFeatureVisibility,
    setFeatureVisibility: visibilityActions.setFeatureVisibility,
    updateClickedFeatureVisibility: visibilityActions.updateClickedFeatureVisibility,
    requestFeatureColorChange: colorActions.requestFeatureColorChange,
    updateClickedFeatureColor: colorActions.updateClickedFeatureColor,
    handleColorScopeChoice: colorActions.handleColorScopeChoice,
    handleLegendNameCommit: colorActions.handleLegendNameCommit,
    handleLegendRenameChoice: colorActions.handleLegendRenameChoice,
    selectLegendNameOption: colorActions.selectLegendNameOption,
    renameLegendEntry: colorActions.renameLegendEntry,
    handleResetColorChoice: colorActions.handleResetColorChoice,
    resetClickedFeatureFillColor: colorActions.resetClickedFeatureFillColor,
    updateClickedFeatureStroke: colorActions.updateClickedFeatureStroke,
    resetClickedFeatureStroke: colorActions.resetClickedFeatureStroke,
    applyStrokeToAllSiblings: colorActions.applyStrokeToAllSiblings,
    setFeatureColor: colorActions.setFeatureColor,
    attachSvgFeatureHandlers: featureSvgActions.attachSvgFeatureHandlers,
    openFeatureEditorForFeature,
    refreshFeatureOverrides: ruleActions.refreshFeatureOverrides,
    getEditableLabelByFeatureId: labelActions.getEditableLabelByFeatureId,
    syncLabelEditor: labelActions.syncLabelEditor,
    downloadLabelOverrideTable: labelActions.downloadLabelOverrideTable,
    loadLabelOverrideTable: labelActions.loadLabelOverrideTable,
    updateClickedFeatureLabelText: labelActions.updateClickedFeatureLabelText,
    handleLabelTextScopeChoice: labelActions.handleLabelTextScopeChoice,
    handleGlobalLabelModeChoice: labelActions.handleGlobalLabelModeChoice,
    requestLabelTextChangeByFeatureId: labelActions.requestLabelTextChangeByFeatureId,
    requestLabelTextChangeByKey: labelActions.requestLabelTextChangeByKey,
    resetAllLabelTextOverrides: labelActions.resetAllLabelTextOverrides
  };
};
