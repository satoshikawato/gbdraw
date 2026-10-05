import {
  normalizeSpecificRule, buildLegendIntents, createRuleLegendCaptions, rendererLegendRows
} from './specific-color-rules.js';
import { normalizeFeatureSelectorMetadata } from './feature-selector.js';
import { getFeatureColorRuleHash } from './feature-utils.js';
import { featureOverrideValue } from '../services/feature-placement.js';

// Ephemeral Python results belong to feature objects, never a session or a SVG.
// An absent result is pending, not a non-match; a declined one (`matches:
// null`) is settled but unknown until the next Generate. A result is a fact of
// one feature and one rule's content, so a later edit cannot make it stale.
// A reactive copy of a feature reads the results of its raw catalog object.
const matchesByFeature = new WeakMap();
const rawFeature = (feature) => globalThis.window?.Vue?.toRaw?.(feature) ?? feature;
const cacheOf = (feature) => matchesByFeature.get(rawFeature(feature));
const cacheFor = (feature) => {
  const raw = rawFeature(feature);
  if (!matchesByFeature.has(raw)) matchesByFeature.set(raw, new Map());
  return matchesByFeature.get(raw);
};
const ruleKey = (rule) => JSON.stringify([rule.feat, rule.qual, rule.val]);
export const ruleMatchesFeature = (feature, rule) => {
  if (!rule || (rule.feat !== '*' && rule.feat !== feature?.type)) return false;
  return cacheOf(feature)?.get(ruleKey(rule))?.matches ?? null;
};
// A feature of a catalog before schema 5 has no drawn selector values. Where
// its rendered ID carries its source hash, its record was drawn with the
// source's coordinates and the source values are the drawn ones; otherwise (a
// cropped, reverse-complemented, or rotated record) a `location` or
// `record_location` rule is not matched live (R4: the fast path declines what
// it cannot read exactly).
export const DRAWN_SELECTOR_QUALIFIERS = new Set(['location', 'record_location']);
export const drawnSelectorUnknown = (feature) => !feature?.drawnSelector
  && Object.prototype.hasOwnProperty.call(feature || {}, 'drawnSelector')
  && getFeatureColorRuleHash(feature) !== String(feature?.selector?.hash || '');
const declinesLiveMatch = (feature, rule) => drawnSelectorUnknown(feature)
  && DRAWN_SELECTOR_QUALIFIERS.has(String(rule?.qual || '').toLowerCase());
export const ruleMatchDeclined = (feature, rules) => rules
  .some((rule) => cacheOf(feature)?.get(ruleKey(rule))?.declined === true);
export const firstMatchingRule = (feature, rules) => {
  let winner = null;
  let priority = Infinity;
  for (const rule of rules) {
    const result = cacheOf(feature)?.get(ruleKey(rule));
    if (result?.matches && result.priority < priority) {
      winner = rule;
      priority = result.priority;
    }
  }
  return winner;
};
export const ruleMatchesReady = (features, rules) => features.every((feature) =>
  rules.every((rule) => ruleMatchesFeature(feature, rule) !== null || ruleMatchDeclined(feature, [rule]))
);
// Feature visibility rules (gbdraw/features/visibility.py) are matched by the
// same Python helper on the same feature payloads; `resolveFeatureDrawn` in
// app/feature-visibility.js reads the results (R4).
const visibilityRuleKey = (rule) => JSON.stringify([
  'visibility', rule.recordId, rule.featureType, rule.qualifier, rule.value
]);
export const visibilityRuleMatchesFeature = (feature, rule) => (
  cacheOf(feature)?.get(visibilityRuleKey(rule))?.matches ?? null
);
// A catalog feature that Python did not render carries no drawn values (feature
// catalog 5 has them on rendered features only), so a `hash`, `location`, or
// `record_location` rule is not matched live for it.
const DRAWN_VALUE_QUALIFIERS = new Set(['hash', ...DRAWN_SELECTOR_QUALIFIERS]);
const declinesVisibilityMatch = (feature, rule) => {
  const qualifier = String(rule?.qualifier || '').toLowerCase();
  return Object.prototype.hasOwnProperty.call(feature || {}, 'drawnSelector')
    ? drawnSelectorUnknown(feature) && DRAWN_SELECTOR_QUALIFIERS.has(qualifier)
    : DRAWN_VALUE_QUALIFIERS.has(qualifier);
};
// Generate matches `hash`, `location`, and `record_location` rules against the
// drawn feature (D-14, PD-OI-069): a cropped or reverse-complemented record
// draws other coordinates than its source. The catalog gives those drawn
// values (`drawnSelector`, feature catalog 5, OV-02); a feature of an older
// catalog sends the hash its rendered ID carries, and its source values only
// where they are the drawn ones.
export const ruleFeaturePayload = (feature, label = '') => {
  const metadata = normalizeFeatureSelectorMetadata(feature);
  const qualifiers = Object.fromEntries(Object.entries(feature.selector?.qualifiers || feature.qualifiers || metadata.qualifiers)
    .map(([key, values]) => [key, (Array.isArray(values) ? values : [values]).filter(value => value != null).map(String)]));
  const drawn = feature?.drawnSelector;
  const selector = drawn
    ? { hash: drawn.hash, location: drawn.location, record_location: drawn.recordLocation }
    : drawnSelectorUnknown(feature)
      ? { hash: getFeatureColorRuleHash(feature) || null, location: null, record_location: null }
      : {
          hash: getFeatureColorRuleHash(feature) || metadata.stableFeatureId,
          location: metadata.location,
          record_location: metadata.recordLocation || `${metadata.record}:${metadata.position}`
        };
  return { type: metadata.featureType, record: metadata.record, qualifiers, selector, label };
};

export const createRulePreparation = ({
  state, evaluate, pending = { value: false }, notify = () => {}, visibilityRules = () => []
}) => {
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
    featureVisibility: JSON.stringify(Object.values(state.featureOverrides || {})
      .map((row) => [row.scope, row.recordKey, row.biologicalFeatureId, row.featureVisibility])),
    // Physical source, palette, and selector inputs retain their identity while
    // request-owned comparison artifacts are published independently.
    inputFiles: [
      state.files?.c_gb, state.files?.c_gff, state.files?.c_fasta, state.files?.c_depth,
      state.files?.d_color, state.files?.blacklist, state.files?.whitelist, state.files?.qualifier_priority,
      state.files?.c_conservation_blasts, state.files?.c_conservation_blasts_source,
      state.files?.c_conservation_fastas, state.files?.c_conservation_sequence_sources
    ].flatMap(input => Array.isArray(input) ? [input, ...input] : [input])
  });
  const isCurrent = (before) => {
    const after = snapshot();
    return Object.keys(before).every((key) => key === 'linearFiles' || key === 'inputFiles'
      ? before[key].length === after[key].length && before[key].every((file, index) => file === after[key][index])
      : before[key] === after[key]);
  };
  const matchesPrepared = (targets, draft) => !draft.length || (
    draft.every((rule) => validated.has(ruleKey(rule))) && ruleMatchesReady(targets, draft)
  );
  const isPrepared = (rules = state.manualSpecificRules) => {
    const draft = [...new Map(rules.map((rule) => [ruleKey(rule), rule])).values()];
    return matchesPrepared(features(), draft);
  };
  const prepare = (rules = state.manualSpecificRules, options = {}) => {
    const targets = features();
    const draft = [...new Map(rules.map((rule) => [ruleKey(rule), { feat: rule.feat, qual: rule.qual, val: rule.val }])).values()];
    // Empty catalogs still require syntax validation at input boundaries.
    if (matchesPrepared(targets, draft)) return true;
    const before = snapshot();
    pending.value = ++pendingCount > 0;
    return evaluate({ features: targets.map((feature) => ruleFeaturePayload(feature)), rules: draft, kind: 'color' }, options)
      .then((result) => {
        if (!isCurrent(before) || state.sessionOperationAvailability?.()) return false;
        validated = new Set(draft.map(ruleKey));
        targets.forEach((feature, index) => {
          const cache = cacheFor(feature);
          const matches = new Set(result.matches[index]);
          draft.forEach((rule, ruleIndex) => cache.set(ruleKey(rule), declinesLiveMatch(feature, rule)
            ? { matches: null, declined: true, priority: Infinity }
            : { matches: matches.has(ruleIndex), priority: result.priorities[index][ruleIndex] }));
        });
        return true;
      }).finally(() => { pending.value = --pendingCount > 0; });
  };
  // The visibility rule matches of every catalog feature. A rule Generate
  // rejects fails the evaluation and leaves its matches unknown.
  const prepareVisibility = (rules = visibilityRules(), options = {}) => {
    const draft = [...new Map(rules.map((rule) => [visibilityRuleKey(rule), {
      recordId: rule.recordId, featureType: rule.featureType, qualifier: rule.qualifier,
      value: rule.value, action: rule.action
    }])).values()];
    if (draft.length === 0) return true;
    const targets = features().filter((feature) => draft
      .some((rule) => !cacheOf(feature)?.has(visibilityRuleKey(rule))));
    if (targets.length === 0) return true;
    pending.value = ++pendingCount > 0;
    return Promise.resolve()
      .then(() => evaluate({
        features: targets.map((feature) => ruleFeaturePayload(feature)), rules: draft, kind: 'visibility'
      }, options))
      .then((result) => {
        targets.forEach((feature, index) => {
          const cache = cacheFor(feature);
          const matches = new Set(result.matches[index]);
          draft.forEach((rule, ruleIndex) => cache.set(visibilityRuleKey(rule), declinesVisibilityMatch(feature, rule)
            ? { matches: null, declined: true }
            : { matches: matches.has(ruleIndex) }));
        });
        return true;
      }, () => false)
      .finally(() => { pending.value = --pendingCount > 0; });
  };
  // `drawn` also prepares the visibility rule matches `resolveFeatureDrawn` reads.
  const run = (rules, commit, { drawn = false } = {}) => {
    if (state.sessionOperationAvailability?.()) return state.sessionOperationAvailability();
    const prepared = prepare(rules);
    const visibility = drawn ? prepareVisibility() : true;
    if (prepared === true && visibility === true) return commit();
    return Promise.all([prepared, visibility]).then(([current]) =>
      state.sessionOperationAvailability?.() || (current ? commit() : undefined));
  };
  // Everything `resolveFeatureDrawn` reads: the visibility rule matches and,
  // for a feature of a type the request does not select, the color rule
  // matches. Never rejects; what stays unknown is resolved as unknown.
  const prepareDrawn = (options = {}) => {
    const colors = prepare(state.manualSpecificRules || [], options);
    const visibility = prepareVisibility(undefined, options);
    if (colors === true && visibility === true) return true;
    return Promise.all([Promise.resolve(colors).catch(() => false), visibility]).then(() => true);
  };
  // `retiredLegendIntents` are rows this commit replaces; they are no renderer
  // rows for the N-06 caption allocation.
  const prepareCandidate = async (rules = state.manualSpecificRules, options = {}, { retiredLegendIntents = [] } = {}) => {
    const before = snapshot();
    const source = rules.map(rule => normalizeSpecificRule(rule));
    const response = await evaluate({ features: [], rules: source, kind: 'color-captions' }, options);
    if (!isCurrent(before)) return null;
    const normalized = response.rules;
    if (!await prepare(normalized, options) || !isCurrent(before)) return null;
    const rendered = (state.extractedFeatures.value || []).filter(feature =>
      featureOverrideValue(state.featureOverrides, feature, 'featureVisibility') !== 'off');
    const used = new Set(rendered.map(feature => firstMatchingRule(feature, normalized)).filter(Boolean));
    const current = state.manualSpecificRules || [];
    const rendererRows = rendererLegendRows({
      legendEntries: state.legendEntries?.value,
      originalLegendOrder: state.originalLegendOrder?.value,
      rules: [...current, ...normalized,
        ...retiredLegendIntents.map(intent => ({ cap: intent?.caption, color: intent?.color }))]
    });
    const intents = buildLegendIntents(normalized.filter(rule => used.has(rule)), rendererRows).intents;
    // The rows the current rules draw, so the commit retires the row Generate drew.
    const currentCaption = createRuleLegendCaptions(current, rendererRows);
    const previousIntents = current.filter(rule => rule.cap)
      .map(rule => ({ caption: currentCaption(rule), color: rule.color }));
    const changes = normalized.flatMap((rule, index) => rule.cap !== source[index].cap
      ? [{ index, before: source[index].cap, after: rule.cap }] : []);
    // Rebind existing rule-derived overrides by their source row, never by suffix parsing.
    const featureColorOverrides = Object.fromEntries(Object.entries(state.featureColorOverrides || {}).map(([key, override]) => {
      const index = source.findIndex(rule => rule.cap === override?.caption
        && rule.color === String(override?.color || '').toLowerCase());
      return [key, index < 0 ? override : { ...override, caption: normalized[index].cap }];
    }));
    return { rules: normalized, intents, previousIntents, changes, featureColorOverrides, snapshot: before };
  };
  const notifyChanges = (candidate) => {
    if (candidate?.changes.length) notify(`Updated ${candidate.changes.length} specific-color caption(s) to distinguish their colors.`);
  };
  return {
    prepare, prepareVisibility, prepareDrawn, isPrepared, prepareCandidate, notifyChanges, run, evaluate, snapshot,
    isCurrent
  };
};
