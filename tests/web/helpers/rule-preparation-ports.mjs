// The members of a rule preparation (app/rule-matching.js) that the rule
// owner takes as ports (createFeatureRuleActions), picked by name: a test
// that builds a real preparation passes these as feature-editor.js does.
/** @param {import('../../../gbdraw/web/js/app/rule-matching.js').RulePreparation} preparation */
export const ruleOwnerPorts = (preparation) => ({
  isRuleInputCurrent: preparation.isCurrent,
  snapshotRuleInputs: preparation.snapshot,
  isRulePrepared: preparation.isPrepared,
  prepareRuleCandidate: preparation.prepareCandidate,
  synthesizeGroupRules: preparation.synthesizeGroup
});
