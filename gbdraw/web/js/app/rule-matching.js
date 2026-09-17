import { normalizeFeatureSelectorMetadata } from './feature-selector.js';

// Ephemeral Python results belong to feature objects, never a session or a SVG.
// An absent result is pending, not a non-match.
const matchesByFeature = new WeakMap();
const ruleKey = (rule) => JSON.stringify([rule.feat, rule.qual, rule.val]);
export const ruleMatchesFeature = (feature, rule) => {
  if (!rule || (rule.feat !== '*' && rule.feat !== feature?.type)) return false;
  return matchesByFeature.get(feature)?.get(ruleKey(rule))?.matches ?? null;
};
export const firstMatchingRule = (feature, rules) => {
  let winner = null;
  let priority = Infinity;
  for (const rule of rules) {
    const result = matchesByFeature.get(feature)?.get(ruleKey(rule));
    if (result?.matches && result.priority < priority) {
      winner = rule;
      priority = result.priority;
    }
  }
  return winner;
};
export const ruleMatchesReady = (features, rules) => features.every((feature) =>
  rules.every((rule) => ruleMatchesFeature(feature, rule) !== null)
);
export const ruleFeaturePayload = (feature, label = '') => {
  const metadata = normalizeFeatureSelectorMetadata(feature);
  return {
    type: metadata.featureType,
    record: metadata.record,
    qualifiers: Object.fromEntries(Object.entries(feature.selector?.qualifiers || feature.qualifiers || metadata.qualifiers)
      .map(([key, values]) => [key, (Array.isArray(values) ? values : [values]).filter(value => value != null).map(String)])),
    selector: {
      hash: metadata.stableFeatureId,
      location: metadata.location,
      record_location: metadata.recordLocation || `${metadata.record}:${metadata.position}`
    },
    label
  };
};

export const createRulePreparation = ({ state, evaluate, pending = { value: false } }) => {
  let validated = new Set();
  let pendingCount = 0;
  const features = () => [...new Set([
    ...(state.extractedFeatures.value || []), ...(state.biologicalFeatures?.value || [])
  ])];
  const snapshot = () => ({
    catalog: state.extractedFeatures.value,
    biology: state.biologicalFeatures?.value,
    result: state.svgResultIdentity?.value,
    mode: state.mode?.value,
    rules: JSON.stringify(state.manualSpecificRules)
  });
  const isCurrent = (before) => {
    const after = snapshot();
    return Object.keys(before).every((key) => before[key] === after[key]);
  };
  const prepare = (rules = state.manualSpecificRules) => {
    const targets = features();
    const draft = [...new Map(rules.map((rule) => [ruleKey(rule), { feat: rule.feat, qual: rule.qual, val: rule.val }])).values()];
    if (!draft.length) return true;
    // Empty catalogs still require syntax validation at input boundaries.
    if (draft.every((rule) => validated.has(ruleKey(rule))) && ruleMatchesReady(targets, draft)) return true;
    const before = snapshot();
    pending.value = ++pendingCount > 0;
    return evaluate({ features: targets.map((feature) => ruleFeaturePayload(feature)), rules: draft, kind: 'color' })
      .then((result) => {
        if (!isCurrent(before)) return false;
        validated = new Set(draft.map(ruleKey));
        targets.forEach((feature, index) => {
          const cache = new Map();
          matchesByFeature.set(feature, cache);
          const matches = new Set(result.matches[index]);
          draft.forEach((rule, ruleIndex) => cache.set(ruleKey(rule), { matches: matches.has(ruleIndex), priority: result.priorities[index][ruleIndex] }));
        });
        return true;
      }).finally(() => { pending.value = --pendingCount > 0; });
  };
  const run = (rules, commit) => {
    const prepared = prepare(rules);
    if (prepared === true) return commit();
    return Promise.resolve(prepared).then((current) => current ? commit() : undefined);
  };
  return { prepare, run, evaluate, snapshot, isCurrent };
};
