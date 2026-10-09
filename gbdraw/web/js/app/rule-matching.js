// @ts-check
import { normalizeSpecificRule } from '../services/specific-color-rules.js';
import { normalizeFeatureSelectorMetadata } from '../services/feature-selector.js';
import { getFeatureColorRuleHash } from '../services/feature-utils.js';
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import {
  recordRuleMatches, recordVisibilityMatch, ruleKey, ruleKeysPending, ruleMatcher, visibilityRuleKnown
} from '../services/rule-matchers.js';

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
// The indexes of the `rules` a feature declines (`drawnSelectorUnknown`).
const liveMatchDeclines = (rules) => {
  const drawnRules = rules.flatMap((rule, index) => (
    DRAWN_SELECTOR_QUALIFIERS.has(String(rule?.qual || '').toLowerCase()) ? [index] : []
  ));
  return (feature) => (drawnSelectorUnknown(feature) ? drawnRules : []);
};
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
// The payload of a feature without a label is built once per feature object.
/** @type {WeakMap<object, Record<string, any>>} */
const payloads = new WeakMap();
export const ruleFeaturePayload = (feature, label = '') => {
  const source = globalThis.window?.Vue?.toRaw?.(feature) ?? feature;
  const known = label === '' ? payloads.get(source) : undefined;
  if (known) return known;
  const payload = buildRuleFeaturePayload(feature, label);
  if (label === '') payloads.set(source, payload);
  return payload;
};
const buildRuleFeaturePayload = (feature, label) => {
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

// The feature color overrides a rule drew follow the rule's normalized caption
// (`prepareCandidate`): rebound by their source row, never by suffix parsing.
// A direct override keeps its caption.
export const rebindRuleColorOverrides = (overrides, rules, normalized) => {
  const source = rules.map((rule) => normalizeSpecificRule(rule));
  return Object.fromEntries(Object.entries(overrides || {}).map(([key, override]) => {
    const index = source.findIndex((rule) => rule.cap === override?.caption
      && rule.color === String(override?.color || '').toLowerCase());
    return [key, index < 0 ? override : { ...override, caption: normalized[index].cap }];
  }));
};

// A run (`runWhenPrepared`) that fails, in its preparation or its commit,
// shows its error as the caller's `operation`, unless another error replaced
// the one shown when it started. Every owner that runs through the rule
// preparation reports a failure this way.
export const reportRuleRunFailure = (state, operation, run) => {
  const previousAlert = state.errorLog?.value;
  const result = run();
  return result?.catch ? result.catch((error) => {
    if (state.errorLog && state.errorLog.value === previousAlert) {
      state.errorLog.value = normalizeUserFacingError(error, { operation, stage: 'helper' });
    }
  }) : result;
};

// Runs `commit` once the matches it reads are prepared (`preparations`, the
// color rule matches first): at once when they are, after them when the color
// matches are still current, and not while a session operation runs. A failed
// color preparation rejects (`reportRuleRunFailure`).
/**
 * @param {Record<string, any>} state App state (state.js; not yet typed).
 * @param {() => (boolean | Promise<any>)[]} preparations
 * @param {() => any} commit
 */
export const runWhenPrepared = (state, preparations, commit) => {
  if (state.sessionOperationAvailability?.()) return state.sessionOperationAvailability();
  const prepared = preparations();
  if (prepared.every((value) => value === true)) return commit();
  return Promise.all(prepared).then(([current]) =>
    state.sessionOperationAvailability?.() || (current ? commit() : undefined));
};

/**
 * Python's answer to a rule evaluation (R7): per feature, the indexes of the
 * rules it matches and, for color rules, the priority of each match (aligned
 * with the indexes); or the normalized rows of `color-captions`.
 * @typedef {Record<string, any>} RuleEvaluationResult
 */

/**
 * @typedef {object} RulePreparationOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(payload: Record<string, any>, options?: Record<string, any>) => Promise<RuleEvaluationResult>} evaluate
 *   The diagram helper that matches rules against feature payloads (R7).
 * @property {{ value: boolean }} [pending] Receives whether a preparation is running.
 * @property {(notice: string) => void} [notify] Shows the notice of a caption Python changed.
 * @property {() => Record<string, any>[]} [visibilityRules] The Feature visibility rule rows, in table order.
 */

/**
 * The candidate rules a commit or a run admits (`prepareCandidate`).
 * @typedef {object} RuleCandidate
 * @property {Record<string, any>[]} rules Python's normalized rule rows.
 * @property {{ index: number, before: string, after: string }[]} changes The captions it changed.
 * @property {Record<string, any>} snapshot The inputs the candidate was prepared from.
 */

/**
 * The rule preparation: every owner that reads rule matches runs through it
 * (R13), as the whole-object `rulePreparation` or as one function of it.
 * @typedef {object} RulePreparation
 * @property {(rules?: Record<string, any>[], options?: Record<string, any>) => boolean | Promise<boolean>} prepare
 *   True at once when the matches of `rules` are prepared, else a promise of whether they are now.
 * @property {(rules?: Record<string, any>[]) => void} retain
 *   Keeps the rules a History restore replaces, so the next preparation matches them with the restored ones.
 * @property {(options?: { strict?: boolean }) => boolean | Promise<boolean | { error: any }>} prepareDrawn
 *   Prepares what `resolveFeatureDrawn` reads; resolves to `{ error }` when Generate rejects the visibility rule table.
 *   `strict` rejects when the color preparation fails and resolves to false when it is stale, as a run does.
 * @property {(rules?: Record<string, any>[]) => boolean} isPrepared
 * @property {(rules?: Record<string, any>[], options?: Record<string, any>) => Promise<RuleCandidate | null>} prepareCandidate
 *   `options.captions` false prepares the rules as given, without Python's caption normalization.
 * @property {(candidate: RuleCandidate | null) => void} notifyChanges
 * @property {() => Record<string, any>} snapshot The inputs the matches depend on.
 * @property {(before: Record<string, any>) => boolean} isCurrent Whether `before` is still the current inputs.
 * @property {{ value: boolean }} pending
 */

/**
 * @param {RulePreparationOptions} options
 * @returns {RulePreparation}
 */
export const createRulePreparation = ({
  state, evaluate, pending = { value: false }, notify = () => {}, visibilityRules = () => []
}) => {
  // The rule keys Python accepted: a key's syntax is a fact of its content.
  /** @type {Set<string>} */
  const validated = new Set();
  let pendingCount = 0;
  const features = () => [...new Set([
    ...(state.extractedFeatures.value || []), ...(state.biologicalFeatures?.value || [])
  ])];
  // The drawing inputs the matches depend on, as JSON. Each is a computed of
  // the active drawing: it is stringified again only after one of the reactive
  // values it read changed, not on each currency check (a rule commit runs
  // about ten), so the N-sized feature color overrides are stringified at most
  // once per change.
  const computed = globalThis.window?.Vue?.computed ?? ((getter) => ({ get value() { return getter(); } }));
  /** @type {Record<string, (drawing: Record<string, any>) => any>} */
  const drawingInputReaders = {
    rules: (drawing) => drawing.manualSpecificRules,
    legend: (drawing) => drawing.legendEntries?.value || [],
    legendColors: (drawing) => drawing.legendColorOverrides || {},
    legendStrokes: (drawing) => drawing.legendStrokeOverrides || {},
    featureColors: (drawing) => drawing.featureColorOverrides || {},
    featureVisibility: (drawing) => Object.values(drawing.featureOverrides || {})
      .map((row) => [row.recordKey, row.biologicalFeatureId, row.featureVisibility])
  };
  const drawingInputs = Object.entries(drawingInputReaders).map(([key, read]) => (
    /** @type {[string, { value: string }]} */ ([key, computed(() => JSON.stringify(read(state.activeDrawing())))])
  ));
  const snapshot = () => ({
    catalog: state.extractedFeatures.value,
    biology: state.biologicalFeatures?.value,
    result: state.svgResultIdentity?.value,
    mode: state.mode?.value,
    inputKinds: JSON.stringify([state.cInputType?.value, state.lInputType?.value]),
    file: state.files?.t_color,
    linearFiles: [...(state.linearSeqs || [])].flatMap(sequence => [sequence.gb, sequence.gff, sequence.fasta]),
    resultNames: JSON.stringify((state.results?.value || []).map(result => result.name)),
    selectedResult: state.selectedResultIndex?.value,
    ...Object.fromEntries(drawingInputs.map(([key, input]) => [key, input.value])),
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
    draft.every((rule) => validated.has(ruleKey(rule))) && ruleMatcher(draft).ready(targets)
  );
  const isPrepared = (rules = state.activeDrawing().manualSpecificRules) => {
    const draft = [...new Map(rules.map((rule) => [ruleKey(rule), rule])).values()];
    return matchesPrepared(features(), draft);
  };
  // Rules a restore replaces: the next preparation matches them with the
  // restored ones, in one evaluation, so the Legend change of the restore is
  // read from known matches (`retain`).
  let retained = [];
  const retain = (rules = []) => { retained = rules; };
  const prepare = (rules = state.activeDrawing().manualSpecificRules, options = {}) => {
    const targets = features();
    const draft = [...new Map([...rules, ...retained].map((rule) => [ruleKey(rule), { feat: rule.feat, qual: rule.qual, val: rule.val }])).values()];
    // Empty catalogs still require syntax validation at input boundaries.
    if (matchesPrepared(targets, draft)) return true;
    // Only the rules without a result for every target, or not yet validated,
    // go to Python; the others' results are already known.
    const pendingKeys = new Set(ruleKeysPending(targets, draft.map(ruleKey)));
    const sent = draft.filter((rule) => pendingKeys.has(ruleKey(rule)) || !validated.has(ruleKey(rule)));
    const before = snapshot();
    pending.value = ++pendingCount > 0;
    return evaluate({ features: targets.map((feature) => ruleFeaturePayload(feature)), rules: sent, kind: 'color' }, options)
      .then((result) => {
        if (!isCurrent(before) || state.sessionOperationAvailability?.()) return false;
        const keys = sent.map(ruleKey);
        keys.forEach((key) => validated.add(key));
        const declines = liveMatchDeclines(sent);
        recordRuleMatches(targets, keys, (index) => ({
          matched: result.matches[index], priorities: result.priorities[index], declined: declines(targets[index])
        }));
        return true;
      }).finally(() => { pending.value = --pendingCount > 0; });
  };
  // The visibility rule matches of every catalog feature. The rules go to
  // Python in table order, so a table Generate rejects (an invalid regex)
  // fails with Generate's error and row: its matches stay unknown, and the
  // preparation resolves to `{ error }`.
  /** @returns {boolean | Promise<boolean | { error: any }>} */
  const prepareVisibility = () => {
    const draft = visibilityRules().map((rule) => ({
      recordId: rule.recordId, featureType: rule.featureType, qualifier: rule.qualifier,
      value: rule.value, action: rule.action
    }));
    if (draft.length === 0) return true;
    const targets = features().filter((feature) => draft.some((rule) => !visibilityRuleKnown(feature, rule)));
    if (targets.length === 0) return true;
    pending.value = ++pendingCount > 0;
    return Promise.resolve()
      .then(() => evaluate({
        features: targets.map((feature) => ruleFeaturePayload(feature)), rules: draft, kind: 'visibility'
      }))
      .then((result) => {
        targets.forEach((feature, index) => {
          const matches = new Set(result.matches[index]);
          draft.forEach((rule, ruleIndex) => recordVisibilityMatch(feature, rule, declinesVisibilityMatch(feature, rule)
            ? { matches: null, declined: true }
            : { matches: matches.has(ruleIndex) }));
        });
        return true;
      }, (error) => ({ error }))
      .finally(() => { pending.value = --pendingCount > 0; });
  };
  // Everything `resolveFeatureDrawn` reads: the visibility rule matches and,
  // for a feature of a type the request does not select, the color rule
  // matches. Never rejects; what stays unknown is resolved as unknown. Resolves
  // to `{ error }` when Generate would reject the visibility rule table. A
  // `strict` preparation is a run's: a failed color preparation rejects and a
  // stale one resolves to false (`runWhenPrepared`).
  const prepareDrawn = ({ strict = false } = {}) => {
    const drawing = state.activeDrawing();
    const colors = prepare(drawing.manualSpecificRules || []);
    const visibility = prepareVisibility();
    if (colors === true && visibility === true) return true;
    return Promise.all([strict ? colors : Promise.resolve(colors).catch(() => false), visibility])
      .then(([current, outcome]) => strict && !current ? false : /** @type {any} */ (outcome)?.error ? outcome : true);
  };
  // The rules a commit or a run admits: Python's normalized captions, their
  // matches prepared, and the captions it changed (`notifyChanges`). `options`
  // are the caller's evaluation options (a run's progress observer). A run
  // that only reads the matches of rules it keeps as they are asks for
  // `captions: false`: the rules are prepared as given and no caption changes.
  const prepareCandidate = async (rules = state.activeDrawing().manualSpecificRules, { captions = true, ...options } = {}) => {
    const before = snapshot();
    if (!captions) {
      return await prepare(rules, options) && isCurrent(before)
        ? { rules, changes: [], snapshot: before } : null;
    }
    const source = rules.map(rule => normalizeSpecificRule(rule));
    const response = await evaluate({ features: [], rules: source, kind: 'color-captions' }, options);
    if (!isCurrent(before)) return null;
    const normalized = response.rules;
    if (!await prepare(normalized, options) || !isCurrent(before)) return null;
    const changes = normalized.flatMap((rule, index) => rule.cap !== source[index].cap
      ? [{ index, before: source[index].cap, after: rule.cap }] : []);
    return { rules: normalized, changes, snapshot: before };
  };
  const notifyChanges = (candidate) => {
    if (candidate?.changes.length) notify(`Updated ${candidate.changes.length} specific-color caption(s) to distinguish their colors.`);
  };
  return {
    prepare, retain, prepareDrawn, isPrepared, prepareCandidate, notifyChanges, snapshot,
    isCurrent, pending
  };
};
