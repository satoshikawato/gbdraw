// @ts-check
// Ephemeral Python results belong to feature objects, never a session or a SVG.
// An absent result is pending, not a non-match; a declined one is settled but
// unknown until the next Generate. A result is a fact of one feature and one
// rule's content, so a later edit cannot make it stale. A reactive copy of a
// feature reads the results of its raw catalog object.
//
// A feature holds only its color matches (rule key -> priority) and the keys
// it declines; a key it was evaluated for and does not hold is a non-match.
// The features of one evaluation share one set of evaluated keys (OV-193), so
// no entry is written per feature and rule.
/**
 * @typedef {object} FeatureRuleResults
 * @property {Map<string, number>} matches Matched color rule key -> priority.
 * @property {Set<string>} declined Color rule keys settled as unknown.
 * @property {Set<string>} evaluated Color rule keys with a result (shared).
 * @property {Map<string, { matches: boolean | null, declined?: boolean }>} visibility
 */
/** @type {WeakMap<object, FeatureRuleResults>} */
const resultsByFeature = new WeakMap();
/** @type {Set<string>} */
const NOTHING_EVALUATED = new Set();
const raw = (value) => globalThis.window?.Vue?.toRaw?.(value) ?? value;
const resultsOf = (feature) => resultsByFeature.get(raw(feature));
const resultsFor = (feature) => {
  const key = raw(feature);
  let results = resultsByFeature.get(key);
  if (!results) {
    results = { matches: new Map(), declined: new Set(), evaluated: NOTHING_EVALUATED, visibility: new Map() };
    resultsByFeature.set(key, results);
  }
  return results;
};

// A rule's key is built once per content: a rule edited in place
// (`Object.assign`) gets a new key on its next read.
/**
 * @param {WeakMap<object, { fields: unknown[], key: string }>} memo
 * @param {(rule: Record<string, any>) => unknown[]} fieldsOf
 */
const memoizedKey = (memo, fieldsOf) => (rule) => {
  const fields = fieldsOf(rule);
  const source = raw(rule);
  const known = memo.get(source);
  if (known && known.fields.every((field, index) => field === fields[index])) return known.key;
  const key = JSON.stringify(fields);
  memo.set(source, { fields, key });
  return key;
};
export const ruleKey = memoizedKey(new WeakMap(), (rule) => [rule.feat, rule.qual, rule.val]);
const applies = (rule, feature) => Boolean(rule) && (rule.feat === '*' || rule.feat === feature?.type);
export const ruleMatchesFeature = (feature, rule) => {
  if (!applies(rule, feature)) return false;
  const results = resultsOf(feature);
  const key = ruleKey(rule);
  if (!results || results.declined.has(key)) return null;
  if (results.matches.has(key)) return true;
  return results.evaluated.has(key) ? false : null;
};

/**
 * Records one color rule evaluation. `keys` are the distinct keys evaluated;
 * `resultOf(index)` gives, for `features[index]`, the indexes of the keys it
 * matches with their priorities, and the indexes of the keys it declines.
 * @param {Record<string, any>[]} features
 * @param {string[]} keys
 * @param {(index: number) => { matched: number[], priorities: number[], declined: number[] }} resultOf
 */
export const recordRuleMatches = (features, keys, resultOf) => {
  const evaluated = new Set(keys);
  /** @type {Map<Set<string>, Set<string>>} */
  const merged = new Map();
  features.forEach((feature, index) => {
    const results = resultsFor(feature);
    for (const key of results.matches.keys()) if (evaluated.has(key)) results.matches.delete(key);
    for (const key of results.declined) if (evaluated.has(key)) results.declined.delete(key);
    const { matched, priorities, declined } = resultOf(index);
    declined.forEach((keyIndex) => results.declined.add(keys[keyIndex]));
    matched.forEach((keyIndex, position) => {
      if (!results.declined.has(keys[keyIndex])) results.matches.set(keys[keyIndex], priorities[position]);
    });
    let union = merged.get(results.evaluated);
    if (!union) {
      union = new Set([...results.evaluated, ...evaluated]);
      merged.set(results.evaluated, union);
    }
    results.evaluated = union;
  });
};
// The keys of `keys` without a result for some feature of `features`.
/** @param {Record<string, any>[]} features @param {string[]} keys */
export const ruleKeysPending = (features, keys) => {
  const sets = [...new Set(features.map((feature) => resultsOf(feature)?.evaluated || NOTHING_EVALUATED))];
  return keys.filter((key) => sets.some((set) => !set.has(key)));
};

/**
 * The color rule matches of one rule list, read for many features (one pass):
 * the rule keys are built once, and each feature is read in its own matches.
 * Ties of priority go to the earliest rule.
 * @param {Record<string, any>[]} rules
 */
export const ruleMatcher = (rules) => {
  /** @type {Map<string, { rule: Record<string, any>, index: number }>} */
  const firstByKey = new Map();
  /** @type {Map<string, Set<string>>} */
  const keysByFeat = new Map();
  (rules || []).forEach((rule, index) => {
    if (!rule) return;
    const key = ruleKey(rule);
    if (!firstByKey.has(key)) firstByKey.set(key, { rule, index });
    if (!keysByFeat.has(rule.feat)) keysByFeat.set(rule.feat, new Set());
    keysByFeat.get(rule.feat)?.add(key);
  });
  /** @type {Map<Set<string>, Map<string, boolean>>} */
  const readiness = new Map();
  const ready = (feature) => {
    const evaluated = resultsOf(feature)?.evaluated || NOTHING_EVALUATED;
    const type = feature?.type;
    let byType = readiness.get(evaluated);
    if (!byType) readiness.set(evaluated, byType = new Map());
    if (!byType.has(type)) {
      const keys = [...(keysByFeat.get('*') || []), ...(keysByFeat.get(type) || [])];
      byType.set(type, keys.every((key) => evaluated.has(key)));
    }
    return byType.get(type);
  };
  const declinesApplying = (feature) => [...(resultsOf(feature)?.declined || [])]
    .some((key) => applies(firstByKey.get(key)?.rule, feature));
  const first = (feature) => {
    /** @type {{ rule: Record<string, any>, index: number } | null} */
    let winner = null;
    let priority = Infinity;
    for (const [key, rank] of resultsOf(feature)?.matches || []) {
      const entry = firstByKey.get(key);
      if (entry && (rank < priority || (rank === priority && winner && entry.index < winner.index))) {
        winner = entry;
        priority = rank;
      }
    }
    return winner?.rule ?? null;
  };
  return {
    // Every rule of the feature's type has a result (a match, a non-match, or declined).
    ready: (features) => features.every(ready),
    // Some rule of the list is declined for the feature: it keeps what Generate drew.
    declined: (feature) => [...(resultsOf(feature)?.declined || [])].some((key) => firstByKey.has(key)),
    // The rule that colors the feature, or null.
    first,
    // `first`, for a reader that must not take an unknown match for a miss:
    // undefined while a match of a rule of the feature's type is pending or declined.
    firstIfKnown: (feature) => (ready(feature) && !declinesApplying(feature) ? first(feature) : undefined),
    // Whether some rule matches the feature: true, false, or null while unknown.
    matchesAny: (feature) => {
      for (const key of resultsOf(feature)?.matches.keys() || []) {
        if (applies(firstByKey.get(key)?.rule, feature)) return true;
      }
      return ready(feature) && !declinesApplying(feature) ? false : null;
    }
  };
};
/** @typedef {ReturnType<typeof ruleMatcher>} RuleMatcher */

// Feature visibility rules (gbdraw/features/visibility.py) are matched by the
// same Python helper on the same feature payloads; `resolveFeatureDrawn` in
// services/feature-visibility.js reads the results (R4).
export const visibilityRuleKey = memoizedKey(new WeakMap(), (rule) => [
  'visibility', rule.recordId, rule.featureType, rule.qualifier, rule.value
]);
/** @param {Record<string, any>} feature @param {Record<string, any>} rule */
export const visibilityRuleKnown = (feature, rule) => Boolean(resultsOf(feature)?.visibility.has(visibilityRuleKey(rule)));
export const visibilityRuleMatchesFeature = (feature, rule) => (
  resultsOf(feature)?.visibility.get(visibilityRuleKey(rule))?.matches ?? null
);
/**
 * @param {Record<string, any>} feature
 * @param {Record<string, any>} rule
 * @param {{ matches: boolean | null, declined?: boolean }} result
 */
export const recordVisibilityMatch = (feature, rule, result) => {
  resultsFor(feature).visibility.set(visibilityRuleKey(rule), result);
};
// The keys of the color rules a feature matches.
/** @param {Record<string, any>} feature @returns {Iterable<string>} */
export const matchedRuleKeys = (feature) => resultsOf(feature)?.matches.keys() || [];
