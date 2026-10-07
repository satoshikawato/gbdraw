// @ts-check
// Ephemeral Python results belong to feature objects, never a session or a SVG.
// An absent result is pending, not a non-match; a declined one (`matches:
// null`) is settled but unknown until the next Generate. A result is a fact of
// one feature and one rule's content, so a later edit cannot make it stale.
// A reactive copy of a feature reads the results of its raw catalog object.
const matchesByFeature = new WeakMap();
const rawFeature = (feature) => globalThis.window?.Vue?.toRaw?.(feature) ?? feature;
export const cacheOf = (feature) => matchesByFeature.get(rawFeature(feature));
export const cacheFor = (feature) => {
  const raw = rawFeature(feature);
  if (!matchesByFeature.has(raw)) matchesByFeature.set(raw, new Map());
  return matchesByFeature.get(raw);
};
export const ruleKey = (rule) => JSON.stringify([rule.feat, rule.qual, rule.val]);
export const ruleMatchesFeature = (feature, rule) => {
  if (!rule || (rule.feat !== '*' && rule.feat !== feature?.type)) return false;
  return cacheOf(feature)?.get(ruleKey(rule))?.matches ?? null;
};
// `firstMatchingRule` for a reader that must not take an unknown match for a
// miss: undefined while a match of a rule of the feature's type is pending or
// declined.
export const firstMatchingRuleIfKnown = (feature, rules) => {
  const cache = cacheOf(feature);
  /** @type {Record<string, any> | null} */
  let winner = null;
  let priority = Infinity;
  for (const rule of rules) {
    if (!rule || (rule.feat !== '*' && rule.feat !== feature?.type)) continue;
    const result = cache?.get(ruleKey(rule));
    if (!result || result.matches === null) return undefined;
    if (result.matches && result.priority < priority) {
      winner = rule;
      priority = result.priority;
    }
  }
  return winner;
};
// Feature visibility rules (gbdraw/features/visibility.py) are matched by the
// same Python helper on the same feature payloads; `resolveFeatureDrawn` in
// services/feature-visibility.js reads the results (R4).
export const visibilityRuleKey = (rule) => JSON.stringify([
  'visibility', rule.recordId, rule.featureType, rule.qualifier, rule.value
]);
export const visibilityRuleMatchesFeature = (feature, rule) => (
  cacheOf(feature)?.get(visibilityRuleKey(rule))?.matches ?? null
);
