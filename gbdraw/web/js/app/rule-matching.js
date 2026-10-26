// @ts-check
import { normalizeSpecificRule } from '../services/specific-color-rules.js';
import { normalizeFeatureSelectorMetadata } from '../services/feature-selector.js';
import { getFeatureColorRuleHash, getFeatureInstance, getRenderedFeatureHash } from '../services/feature-utils.js';
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { toRaw } from '../utils/vue-raw.js';
import {
  recordRuleMatches, recordVisibilityMatch, ruleKey, ruleKeysPending, ruleMatcher, visibilityRuleKnown
} from '../services/rule-matchers.js';
/** @import { FeatureLike } from '../services/feature-utils.js' */
/** @import { ColorRule, RuleFeature, VisibilityRule } from '../services/rule-matchers.js' */
/** @import { SpecificColorRule } from '../services/specific-color-rules.js' */
/** @import { RuleFeaturePayload, SynthesizeColorRulesRequest, SynthesizedColorRules } from '../contracts/rule-synthesis.generated.js' */

// A feature of a catalog before schema 5 has no drawn selector values. Where
// its rendered ID carries its source hash, its record was drawn with the
// source's coordinates and the source values are the drawn ones; otherwise (a
// cropped, reverse-complemented, or rotated record) a `location` or
// `record_location` rule is not matched live (R4: the fast path declines what
// it cannot read exactly).
export const DRAWN_SELECTOR_QUALIFIERS = new Set(['location', 'record_location']);
/**
 * The values a feature of feature catalog 5 was drawn with (`null` once unknown).
 * @typedef {{ hash: string | null, location: string | null, recordLocation: string | null }} DrawnSelector
 */
/**
 * A feature the rule helpers read: its catalog row, with the drawn selector
 * when it has one (owned by `services/feature-catalog.js`).
 * @typedef {FeatureLike & { drawnSelector?: DrawnSelector | null }} RuleMatchedFeature
 */
/** @param {RuleMatchedFeature | null | undefined} feature */
export const drawnSelectorUnknown = (feature) => !feature?.drawnSelector
  && Object.prototype.hasOwnProperty.call(feature || {}, 'drawnSelector')
  && getRenderedFeatureHash(feature) !== String(feature?.selector?.hash || '');
// The indexes of the `rules` a feature declines (`drawnSelectorUnknown`).
/** @param {readonly ColorRule[]} rules */
const liveMatchDeclines = (rules) => {
  const drawnRules = rules.flatMap((rule, index) => (
    DRAWN_SELECTOR_QUALIFIERS.has(String(rule?.qual || '').toLowerCase()) ? [index] : []
  ));
  /** @param {RuleMatchedFeature} feature */
  return (feature) => (drawnSelectorUnknown(feature) ? drawnRules : []);
};
// A catalog feature that Python did not render carries no drawn values (feature
// catalog 5 has them on rendered features only), so a `location` or
// `record_location` rule is not matched live for it.
/**
 * @param {RuleMatchedFeature} feature
 * @param {VisibilityRule} rule
 */
const declinesVisibilityMatch = (feature, rule) => {
  const qualifier = String(rule?.qualifier || '').toLowerCase();
  return Object.prototype.hasOwnProperty.call(feature || {}, 'drawnSelector')
    ? drawnSelectorUnknown(feature) && DRAWN_SELECTOR_QUALIFIERS.has(qualifier)
    : DRAWN_SELECTOR_QUALIFIERS.has(qualifier);
};
// Generate matches `location` and `record_location` rules against the drawn
// feature (D-14): a cropped or reverse-complemented record draws other
// coordinates than its source. The catalog gives those drawn values
// (`drawnSelector`, feature catalog 5, OV-02); a feature of an older catalog
// sends its source values only where they are the drawn ones. A `hash` rule
// matches the feature's source-record hash on every record (OV-401), and its
// feature instance (`instance`), which names one feature of one record.
// The payload of a feature without a label is built once per feature object;
// its type is Python's (gbdraw/web_support/rule_matching.py, R7).
/** @type {WeakMap<object, RuleFeaturePayload>} */
const payloads = new WeakMap();
/**
 * @param {RuleMatchedFeature} feature
 * @param {string} [label]
 * @returns {RuleFeaturePayload}
 */
export const ruleFeaturePayload = (feature, label = '') => {
  const source = toRaw(feature);
  const known = label === '' ? payloads.get(source) : undefined;
  if (known) return known;
  const payload = buildRuleFeaturePayload(feature, label);
  if (label === '') payloads.set(source, payload);
  return payload;
};
/**
 * @param {RuleMatchedFeature} feature
 * @param {string} label
 * @returns {RuleFeaturePayload}
 */
const buildRuleFeaturePayload = (feature, label) => {
  const metadata = normalizeFeatureSelectorMetadata(feature);
  const qualifiers = Object.fromEntries(Object.entries(feature.selector?.qualifiers || feature.qualifiers || metadata.qualifiers)
    .map(([key, values]) => [key, (Array.isArray(values) ? values : [values]).filter(value => value != null).map(String)]));
  const drawn = feature?.drawnSelector;
  const instance = getFeatureInstance(feature) || null;
  const selector = drawn
    ? { hash: getFeatureColorRuleHash(feature) || null, instance, location: drawn.location, record_location: drawn.recordLocation }
    : drawnSelectorUnknown(feature)
      ? { hash: getFeatureColorRuleHash(feature) || null, instance, location: null, record_location: null }
      : {
          hash: getFeatureColorRuleHash(feature) || metadata.stableFeatureId,
          instance,
          location: metadata.location,
          record_location: metadata.recordLocation || `${metadata.record}:${metadata.position}`
        };
  return { type: metadata.featureType, record: metadata.record, qualifiers, selector, label };
};

// A run (`runWhenPrepared`) that fails, in its preparation or its commit,
// shows its error as the caller's `operation`, unless another error replaced
// the one shown when it started. Every owner that runs through the rule
// preparation reports a failure this way.
/**
 * @param {{ errorLog?: { value: unknown } }} state The error log ref (state.js).
 * @param {string} operation
 * @param {() => unknown} run
 * @returns {unknown} What `run` returned; a promise it returned, with the failure reported.
 */
export const reportRuleRunFailure = (state, operation, run) => {
  const previousAlert = state.errorLog?.value;
  // Every value but null and undefined reads `catch` (a promise's method, or undefined), and `?.` covers those two.
  const result = /** @type {{ catch?: Promise<unknown>['catch'] } | null | undefined} */ (run());
  return result?.catch ? result.catch((error) => {
    if (state.errorLog && state.errorLog.value === previousAlert) {
      state.errorLog.value = normalizeUserFacingError(error, { operation, stage: 'helper' });
    }
  }) : result;
};

// Runs `commit` once the matches it reads are prepared (`preparations`, the
// color rule matches first): at once when they are, after them when the color
// matches are still current (else `stale` answers), and not while a session
// operation runs. A failed color preparation rejects (`reportRuleRunFailure`).
/**
 * @param {Record<string, any>} state App state (state.js; not yet typed).
 * @param {() => (boolean | Promise<any>)[]} preparations
 * @param {() => any} commit
 * @param {() => unknown} [stale]
 */
export const runWhenPrepared = (state, preparations, commit, stale = () => undefined) => {
  if (state.sessionOperationAvailability?.()) return state.sessionOperationAvailability();
  const prepared = preparations();
  if (prepared.every((value) => value === true)) return commit();
  return Promise.all(prepared).then(([current]) =>
    state.sessionOperationAvailability?.() || (current ? commit() : stale()));
};

/**
 * Python's answer to a rule evaluation (R7): per feature, the indexes of the
 * rules it matches and, for color rules, the priority of each match (aligned
 * with the indexes).
 * @typedef {Record<string, any>} RuleEvaluationResult
 */

/**
 * @typedef {object} RulePreparationOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(payload: Record<string, any>, options?: Record<string, any>) => Promise<RuleEvaluationResult>} evaluate
 *   The diagram helper that matches rules against feature payloads (R7).
 * @property {(request: SynthesizeColorRulesRequest) => Promise<SynthesizedColorRules>} synthesize
 *   The diagram helper that writes the rules of a group color edit (OV-244, Owner R-07).
 * @property {{ value: boolean }} [pending] Receives whether a preparation is running.
 * @property {{ value: boolean }} [blocking] Receives whether an edit's preparation is running: one
 *   that is not a `view` preparation (OV-377). Save and Load wait for it.
 * @property {() => Record<string, any>[]} [visibilityRules] The Feature visibility rule rows, in table order.
 */

/**
 * The inputs the matches depend on (`snapshot`): the file lists keep their identity;
 * the other members compare with `===`.
 * @typedef {{ linearFiles: unknown[], inputFiles: unknown[] } & Record<string, unknown>} RuleInputsSnapshot
 */

/**
 * The candidate rules a commit or a run admits (`prepareCandidate`).
 * @typedef {object} RuleCandidate
 * @property {readonly Partial<SpecificColorRule>[]} rules The rule rows, as `normalizeSpecificRule` writes them
 *   (as given with `options.asGiven`); a caption stays as written.
 * @property {RuleInputsSnapshot} snapshot The inputs the candidate was prepared from.
 */

/**
 * The rule preparation: every owner that reads rule matches runs through it
 * (R13), as the whole-object `rulePreparation` or as one function of it.
 * @typedef {object} RulePreparation
 * @property {(rules?: readonly Partial<SpecificColorRule>[], options?: Record<string, any>) => boolean | Promise<boolean>} prepare
 *   True at once when the matches of `rules` are prepared, else a promise of whether they are now.
 *   `options.blocked` says whether a Session operation discards the matches (default: any Session operation);
 *   the other options are the evaluation's.
 * @property {(rules?: readonly Partial<SpecificColorRule>[]) => void} retain
 *   Keeps the rules a History restore replaces, so the next preparation matches them with the restored ones.
 * @property {(options?: { strict?: boolean, view?: boolean }) => boolean | Promise<boolean | { error: any }>} prepareDrawn
 *   Prepares what `resolveFeatureDrawn` reads; resolves to `{ error }` when Generate rejects the visibility rule table.
 *   `strict` rejects when the color preparation fails and resolves to false when it is stale, as a run does.
 *   `view` marks a preparation that only reads the matches (a popup opening): it does not raise `blocking`.
 * @property {(rules?: readonly Partial<SpecificColorRule>[]) => boolean} isPrepared
 * @property {(rules?: readonly Partial<SpecificColorRule>[], options?: Record<string, any>) => Promise<RuleCandidate | null>} prepareCandidate
 *   `options.asGiven` prepares the rules as they are, without writing them as `normalizeSpecificRule` does.
 * @property {(group: GroupColorEdit, rules: readonly Partial<SpecificColorRule>[]) => Promise<SynthesizedColorRules>} synthesizeGroup
 *   Python's rules for coloring `group.features` as one group over the rule table `rules` and the whole catalog.
 * @property {() => RuleInputsSnapshot} snapshot The inputs the matches depend on.
 * @property {(before: RuleInputsSnapshot) => boolean} isCurrent Whether `before` is still the current inputs.
 * @property {{ value: boolean }} pending
 */

/**
 * One group color edit: the selected features, found in the catalog by their
 * rendered id (`svg_id`; a feature without one is not drawn and takes no rule),
 * the label texts each shows (aligned with `features`), and the color and
 * caption they take; `scope` "features" colors each feature by itself.
 * @typedef {object} GroupColorEdit
 * @property {SynthesizeColorRulesRequest['scope']} scope
 * @property {readonly (RuleMatchedFeature & { svg_id?: string | null })[]} features
 * @property {readonly string[][]} labels
 * @property {string} color
 * @property {string} caption
 */

/**
 * @param {RulePreparationOptions} options
 * @returns {RulePreparation}
 */
export const createRulePreparation = ({
  state, evaluate, synthesize, pending = { value: false }, blocking = { value: false }, visibilityRules = () => []
}) => {
  // The rule keys Python accepted: a key's syntax is a fact of its content.
  /** @type {Set<string>} */
  const validated = new Set();
  let pendingCount = 0;
  let blockingCount = 0;
  // Raises `pending` (and `blocking`, unless `view`) until the returned release runs.
  /** @param {boolean} view */
  const track = (view) => {
    pending.value = ++pendingCount > 0;
    if (!view) blocking.value = ++blockingCount > 0;
    return () => {
      pending.value = --pendingCount > 0;
      if (!view) blocking.value = --blockingCount > 0;
    };
  };
  // The catalog and biological features, each once, built again only when
  // either array is replaced (CW-06 row M1). Both arrays are replaced
  // wholesale, never changed in place after adoption, so their raw identities
  // key the union; callers read it and do not change it.
  /** @type {{ catalog: unknown, biology: unknown, union: RuleFeature[] } | null} */
  let union = null;
  /** @returns {RuleFeature[]} */
  const features = () => {
    const catalog = toRaw(state.extractedFeatures.value);
    const biology = toRaw(state.biologicalFeatures?.value);
    if (!union || union.catalog !== catalog || union.biology !== biology) {
      union = { catalog, biology, union: [...new Set([
        ...(state.extractedFeatures.value || []), ...(state.biologicalFeatures?.value || [])
      ])] };
    }
    return union.union;
  };
  // The drawing inputs the matches depend on, as JSON. Each is a computed of
  // the active drawing: it is stringified again only after one of the reactive
  // values it read changed, not on each currency check (a rule commit runs
  // about ten), so the N-sized feature color overrides are stringified at most
  // once per change.
  const computed = globalThis.window?.Vue?.computed ?? ((/** @type {() => string} */ getter) => ({ get value() { return getter(); } }));
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
  /** @returns {RuleInputsSnapshot} */
  const snapshot = () => ({
    catalog: state.extractedFeatures.value,
    biology: state.biologicalFeatures?.value,
    result: state.svgResultIdentity?.value,
    mode: state.mode?.value,
    inputKinds: JSON.stringify([state.cInputType?.value, state.lInputType?.value]),
    file: state.files?.t_color,
    linearFiles: [...(state.linearSeqs || [])].flatMap(sequence => [sequence.gb, sequence.gff, sequence.fasta]),
    resultNames: JSON.stringify((state.results?.value || []).map((/** @type {{ name?: unknown }} */ result) => result.name)),
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
  /** @type {RulePreparation['isCurrent']} */
  const isCurrent = (before) => {
    const after = snapshot();
    return Object.keys(before).every((key) => key === 'linearFiles' || key === 'inputFiles'
      ? before[key].length === after[key].length && before[key].every((file, index) => file === after[key][index])
      : before[key] === after[key]);
  };
  /**
   * @param {RuleFeature[]} targets
   * @param {ColorRule[]} draft
   */
  const matchesPrepared = (targets, draft) => !draft.length || (
    draft.every((rule) => validated.has(ruleKey(rule))) && ruleMatcher(draft).ready(targets)
  );
  /** @type {RulePreparation['isPrepared']} */
  const isPrepared = (rules = state.activeDrawing().manualSpecificRules) => {
    const draft = [...new Map(rules.map((rule) => [ruleKey(rule), rule])).values()];
    return matchesPrepared(features(), draft);
  };
  // Rules a restore replaces: the next preparation matches them with the
  // restored ones, in one evaluation, so the Legend change of the restore is
  // read from known matches (`retain`).
  /** @type {readonly Partial<SpecificColorRule>[]} */
  let retained = [];
  /** @type {RulePreparation['retain']} */
  const retain = (rules = []) => { retained = rules; };
  /** @type {RulePreparation['prepare']} */
  const prepare = (rules = state.activeDrawing().manualSpecificRules, {
    blocked = () => Boolean(state.sessionOperationAvailability?.()), view = false, ...options
  } = /** @type {{ blocked?: () => boolean, view?: boolean } & Record<string, unknown>} */ ({})) => {
    const targets = features();
    const draft = [...new Map([...rules, ...retained].map((rule) => [ruleKey(rule), { feat: rule.feat, qual: rule.qual, val: rule.val }])).values()];
    // Empty catalogs still require syntax validation at input boundaries.
    if (matchesPrepared(targets, draft)) return true;
    // Only the rules without a result for every target, or not yet validated,
    // go to Python; the others' results are already known.
    const pendingKeys = new Set(ruleKeysPending(targets, draft.map(ruleKey)));
    const sent = draft.filter((rule) => pendingKeys.has(ruleKey(rule)) || !validated.has(ruleKey(rule)));
    const before = snapshot();
    const release = track(view);
    return evaluate({ features: targets.map((feature) => ruleFeaturePayload(feature)), rules: sent, kind: 'color' }, options)
      .then((result) => {
        if (!isCurrent(before) || blocked()) return false;
        const keys = sent.map(ruleKey);
        keys.forEach((key) => validated.add(key));
        const declines = liveMatchDeclines(sent);
        recordRuleMatches(targets, keys, (index) => ({
          matched: result.matches[index], priorities: result.priorities[index], declined: declines(targets[index])
        }));
        return true;
      }).finally(release);
  };
  // The visibility rule matches of every catalog feature. The rules go to
  // Python in table order, so a table Generate rejects (an invalid regex)
  // fails with Generate's error and row: its matches stay unknown, and the
  // preparation resolves to `{ error }`.
  /** @returns {boolean | Promise<boolean | { error: any }>} */
  const prepareVisibility = (view = false) => {
    const draft = visibilityRules().map((rule) => ({
      recordId: rule.recordId, featureType: rule.featureType, qualifier: rule.qualifier,
      value: rule.value, action: rule.action
    }));
    if (draft.length === 0) return true;
    const targets = features().filter((feature) => draft.some((rule) => !visibilityRuleKnown(feature, rule)));
    if (targets.length === 0) return true;
    const release = track(view);
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
      .finally(release);
  };
  // Everything `resolveFeatureDrawn` reads: the visibility rule matches and,
  // for a feature of a type the request does not select, the color rule
  // matches. Never rejects; what stays unknown is resolved as unknown. Resolves
  // to `{ error }` when Generate would reject the visibility rule table. A
  // `strict` preparation is a run's: a failed color preparation rejects and a
  // stale one resolves to false (`runWhenPrepared`).
  /** @type {RulePreparation['prepareDrawn']} */
  const prepareDrawn = ({ strict = false, view = false } = {}) => {
    const drawing = state.activeDrawing();
    const colors = prepare(drawing.manualSpecificRules || [], { view });
    const visibility = prepareVisibility(view);
    if (colors === true && visibility === true) return true;
    return Promise.all([strict ? colors : Promise.resolve(colors).catch(() => false), visibility])
      .then(([current, outcome]) => strict && !current ? false : /** @type {any} */ (outcome)?.error ? outcome : true);
  };
  // The rules a commit or a run admits, with their matches prepared: the rows
  // as `normalizeSpecificRule` writes them, each caption as written (gbdraw
  // never rewrites one; the Legend numbers the names it draws). `options` are
  // the caller's evaluation options (a run's progress observer). A run that
  // only reads the matches of rules it keeps as they are asks for `asGiven`.
  /** @type {RulePreparation['prepareCandidate']} */
  const prepareCandidate = async (rules = state.activeDrawing().manualSpecificRules, { asGiven = false, ...options } = {}) => {
    const before = snapshot();
    const candidate = asGiven ? rules : rules.map(rule => normalizeSpecificRule(rule));
    return await prepare(candidate, options) && isCurrent(before) ? { rules: candidate, snapshot: before } : null;
  };
  // Every catalog feature goes to Python, which refuses a condition that
  // matches an unselected one.
  /** @type {RulePreparation['synthesizeGroup']} */
  const synthesizeGroup = ({ scope, features: selected, labels, color, caption }, rules) => {
    const targets = features();
    /** @type {Map<string, string[]>} */
    const labelsById = new Map();
    selected.forEach((feature, index) => { if (feature.svg_id) labelsById.set(feature.svg_id, [...labels[index]]); });
    const picked = targets.flatMap((feature, index) => {
      const id = /** @type {{ svg_id?: unknown }} */ (feature).svg_id;
      return typeof id === 'string' && labelsById.has(id) ? [index] : [];
    });
    return synthesize({
      scope,
      features: targets.map((feature) => ruleFeaturePayload(feature)),
      selected: picked,
      labels: picked.map((index) => labelsById.get(/** @type {{ svg_id: string }} */ (targets[index]).svg_id) || []),
      rules: rules.map((rule) => ({ feat: String(rule.feat ?? ''), qual: String(rule.qual ?? ''), val: String(rule.val ?? '') })),
      color,
      caption
    });
  };
  return {
    prepare, retain, prepareDrawn, isPrepared, prepareCandidate, synthesizeGroup, snapshot,
    isCurrent, pending
  };
};
