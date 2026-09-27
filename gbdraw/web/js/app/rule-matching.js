import { normalizeSpecificRule, buildLegendIntents } from './specific-color-rules.js';
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

export const createRulePreparation = ({ state, evaluate, pending = { value: false }, notify = () => {} }) => {
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
    inputKinds: JSON.stringify([state.cInputType?.value, state.lInputType?.value]),
    rules: JSON.stringify(state.manualSpecificRules),
    file: state.files?.t_color,
    linearFiles: [...(state.linearSeqs || [])].flatMap(sequence => [sequence.gb, sequence.gff, sequence.fasta]),
    resultNames: JSON.stringify((state.results?.value || []).map(result => result.name)),
    selectedResult: state.selectedResultIndex?.value,
    legend: JSON.stringify(state.legendEntries?.value || []),
    legendColors: JSON.stringify(state.legendColorOverrides || {}),
    legendStrokes: JSON.stringify(state.legendStrokeOverrides || {}),
    featureColors: JSON.stringify(state.featureColorOverrides || {}),
    featureVisibility: JSON.stringify(state.featureVisibilityOverrides || {}),
    ...Object.fromEntries(Object.entries(state.files || {}).map(([key, value]) => [`file:${key}`, value]))
  });
  const isCurrent = (before) => {
    const after = snapshot();
    return Object.keys(before).every((key) => key === 'linearFiles'
      ? before[key].length === after[key].length && before[key].every((file, index) => file === after[key][index])
      : before[key] === after[key]);
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
  const prepareCandidate = async (rules = state.manualSpecificRules) => {
    const before = snapshot();
    const source = rules.map(rule => normalizeSpecificRule(rule));
    const response = await evaluate({ features: [], rules: source, kind: 'color-captions' });
    if (!isCurrent(before)) return null;
    const normalized = response.rules;
    if (!await prepare(normalized) || !isCurrent(before)) return null;
    const rendered = (state.extractedFeatures.value || []).filter(feature =>
      state.featureVisibilityOverrides?.[feature.svg_id] !== 'off');
    const used = new Set(rendered.map(feature => firstMatchingRule(feature, normalized)).filter(Boolean));
    const intents = buildLegendIntents(normalized.filter(rule => used.has(rule))).intents;
    const changes = normalized.flatMap((rule, index) => rule.cap !== source[index].cap
      ? [{ index, before: source[index].cap, after: rule.cap }] : []);
    // Rebind existing rule-derived overrides by their source row, never by suffix parsing.
    const featureColorOverrides = Object.fromEntries(Object.entries(state.featureColorOverrides || {}).map(([key, override]) => {
      const index = source.findIndex(rule => rule.cap === override?.caption
        && rule.color === String(override?.color || '').toLowerCase());
      return [key, index < 0 ? override : { ...override, caption: normalized[index].cap }];
    }));
    return { rules: normalized, intents, changes, featureColorOverrides, snapshot: before };
  };
  const notifyChanges = (candidate) => {
    if (candidate?.changes.length) notify(`Updated ${candidate.changes.length} specific-color caption(s) to distinguish their colors.`);
  };
  return { prepare, prepareCandidate, notifyChanges, run, evaluate, snapshot, isCurrent };
};
