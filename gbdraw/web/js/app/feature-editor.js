import { createFeatureColorActions } from './feature-editor/color-actions.js';
import { createFeatureRuleActions } from './feature-editor/rule-actions.js';
import { createFeatureSvgActions } from './feature-editor/svg-actions.js';

export const createFeatureEditor = ({ state, nextTick, legendActions, svgActions }) => {
  const ruleActions = createFeatureRuleActions({ state, nextTick, legendActions });
  const featureSvgActions = createFeatureSvgActions({ state, getFeatureColor: ruleActions.getFeatureColor });
  const colorActions = createFeatureColorActions({
    state,
    nextTick,
    legendActions,
    svgActions,
    ruleActions,
    featureSvgActions
  });

  return {
    addCustomColor: ruleActions.addCustomColor,
    addPriorityRule: ruleActions.addPriorityRule,
    addFeature: ruleActions.addFeature,
    addSpecificRule: ruleActions.addSpecificRule,
    applySpecificRulePreset: ruleActions.applySpecificRulePreset,
    clearAllSpecificRules: ruleActions.clearAllSpecificRules,
    getFeatureColor: ruleActions.getFeatureColor,
    canEditFeatureColor: ruleActions.canEditFeatureColor,
    updateClickedFeatureColor: colorActions.updateClickedFeatureColor,
    handleColorScopeChoice: colorActions.handleColorScopeChoice,
    handleResetColorChoice: colorActions.handleResetColorChoice,
    resetClickedFeatureFillColor: colorActions.resetClickedFeatureFillColor,
    updateClickedFeatureStroke: colorActions.updateClickedFeatureStroke,
    resetClickedFeatureStroke: colorActions.resetClickedFeatureStroke,
    applyStrokeToAllSiblings: colorActions.applyStrokeToAllSiblings,
    setFeatureColor: colorActions.setFeatureColor,
    attachSvgFeatureHandlers: featureSvgActions.attachSvgFeatureHandlers,
    refreshFeatureOverrides: ruleActions.refreshFeatureOverrides
  };
};
