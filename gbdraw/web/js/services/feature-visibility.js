// @ts-check
import { exactRegexValue } from './feature-selector.js';
import { ruleMatcher, visibilityRuleMatchesFeature } from './rule-matchers.js';
import { resultCatalogFeatures, stableFeatureOverrideKey as stableKeyOf } from './feature-catalog.js';
import {
  FEATURE_OVERRIDE_EDIT_FIELDS,
  featureIdentityKey,
  featureIdentityKeyOf,
  featureOverrideValue,
  parseFeatureIdentityKey,
  recordKeyBelongsToRequest,
  rowBelongsToRequest,
  updateFeatureOverride
} from './feature-placement.js';
import { normalizeTsvCell as normalizeCell } from '../utils/tsv-cell.js';
import { cloneJsonData } from './json-clone.js';
import { pythonLegendRows } from './legend-svg.js';
import { featureEditsOf, setFeatureOverrideRow, setFeaturePlacementRow } from './drawing/feature-edits.js';
import { paintOf, setFeatureStroke } from './drawing/paint.js';
import { replaceAnnotationTarget } from './drawing/annotation-sets.js';
import { featureOverrideKey } from './feature-override-identity.js';
import { createRuleLegendCaptions } from './specific-color-rules.js';
/** @import { PaintSlice } from './drawing/paint.js' */
/** @import { FeatureEditsSlice } from './drawing/feature-edits.js' */
/** @import { AnnotationRow, AnnotationSetRow, AnnotationSetsSlice } from './drawing/annotation-sets.js' */
/** @import { AnnotationTarget } from './annotation-state.js' */
/** @import { FeatureOverrideDraft, FeatureRequestRecord } from './feature-placement.js' */
/** @import { RuleMatcher } from './rule-matchers.js' */
/** @import { LegendRowContext, SpecificColorRule } from './specific-color-rules.js' */
/** @import { PythonLegendKey } from './legend-svg.js' */
/** @import { SvgLegendRowFacts } from './svg-result-ingestion.js' */
export { escapeRegexLiteral, exactRegexValue } from './feature-selector.js';

const REQUIRED_COLUMNS = ['record_id', 'feature_type', 'qualifier', 'value', 'action'];
const COMMON_QUALIFIERS = ['product', 'gene', 'protein_id', 'locus_tag', 'hash', 'location', 'record_location'];
const SHOW_ACTIONS = new Set(['show', 'on']);
const HIDE_ACTIONS = new Set(['hide', 'off', 'false', '0']);
const EXCLUDE_MATCHING_ACTIONS = new Set(['exclude_matching', 'suppress']);

let generatedRuleId = 0;

/** @param {unknown} value */
const normalizeSource = (value) => {
  const normalized = normalizeCell(value).toLowerCase();
  return ['manual', 'editor', 'file'].includes(normalized) ? normalized : 'manual';
};

export const featureVisibilityQualifierSuggestions = COMMON_QUALIFIERS;

/** @param {unknown} value */
export const normalizeVisibilityMode = (value) => {
  const normalized = normalizeCell(value).toLowerCase();
  if (normalized === 'suppress') return 'exclude_matching';
  return ['on', 'off', 'exclude_matching'].includes(normalized) ? normalized : 'default';
};

/** @param {unknown} value */
export const normalizeFeatureVisibilityAction = (value) => {
  const normalized = normalizeCell(value).toLowerCase();
  if (SHOW_ACTIONS.has(normalized)) return 'show';
  if (HIDE_ACTIONS.has(normalized)) return 'off';
  if (EXCLUDE_MATCHING_ACTIONS.has(normalized)) return 'exclude_matching';
  return '';
};

/** @param {unknown} value */
export const featureVisibilityActionToMode = (value) => {
  const action = normalizeFeatureVisibilityAction(value);
  if (action === 'show') return 'on';
  if (action === 'off') return 'off';
  if (action === 'exclude_matching') return 'exclude_matching';
  return 'default';
};

/** @param {unknown} value */
export const featureVisibilityModeToAction = (value) => {
  const normalized = normalizeCell(value).toLowerCase();
  if (normalized === 'on') return 'show';
  if (normalized === 'off') return 'off';
  if (normalized === 'exclude_matching' || normalized === 'suppress') return 'exclude_matching';
  if (normalized === 'default') return '';
  return normalizeFeatureVisibilityAction(normalized);
};

const nextRuleId = (prefix = 'feature-visibility-rule') => {
  generatedRuleId += 1;
  return `${prefix}-${generatedRuleId}`;
};

/** @param {string[]} fields */
const isHeaderRow = (fields) => {
  const normalized = fields.map((field) => normalizeCell(field).toLowerCase());
  if (normalized.join('\t') === REQUIRED_COLUMNS.join('\t')) return true;
  return (
    normalized[0] === 'record_id' ||
    normalized[0] === 'record'
  ) &&
    normalized[1] === 'feature_type' &&
    (normalized[2] === 'qualifier' || normalized[2] === 'qualifier_key') &&
    (normalized[3] === 'value' || normalized[3] === 'qualifier_value_regex') &&
    normalized[4] === 'action';
};

/**
 * A normalized Feature visibility rule row (`normalizeFeatureVisibilityRule`).
 * @typedef {object} FeatureVisibilityRule
 * @property {string} id
 * @property {string} source
 * @property {string} featureId
 * @property {string} label
 * @property {string} recordId
 * @property {string} featureType
 * @property {string} qualifier
 * @property {string} value
 * @property {string} action
 */

/**
 * A rule row as the draft or a Session holds it: the row of
 * `FeatureVisibilityRule` under its camel-case or snake-case names, unvalidated.
 * @typedef {object} FeatureVisibilityRuleDraft
 * @property {unknown} [id]
 * @property {unknown} [source]
 * @property {unknown} [featureId]
 * @property {unknown} [feature_id]
 * @property {unknown} [label]
 * @property {unknown} [recordId]
 * @property {unknown} [record_id]
 * @property {unknown} [featureType]
 * @property {unknown} [feature_type]
 * @property {unknown} [qualifier]
 * @property {unknown} [value]
 * @property {unknown} [action]
 */

/**
 * The fields an exact-qualifier rule is built or removed by.
 * @typedef {object} FeatureVisibilityRuleInput
 * @property {string} [featureType]
 * @property {string} [qualifier]
 * @property {string} [value]
 * @property {unknown} [action]
 * @property {string} [label]
 */

/**
 * A feature as the resolver reads it: a catalog row, an extracted feature, or a
 * change that names one (JSON keyed by name).
 * @typedef {{
 *   type?: unknown, record_key?: unknown, recordKey?: unknown,
 *   biological_feature_id?: unknown, biologicalFeatureId?: unknown,
 *   svg_id?: unknown, svgId?: unknown, featureId?: unknown, feature_id?: unknown, id?: unknown
 * }} VisibilityFeature
 */

/**
 * One feature's visibility change: its identity, the modes before and after
 * (`mode` replaces `after` when a History restore sets it).
 * @typedef {VisibilityFeature & { before?: string, after?: unknown, mode?: unknown }} FeatureVisibilityChange
 */

/**
 * What Python reports for an edit of a request; an `unresolved` one names the
 * kinds of edit whose feature the source does not have.
 * @typedef {{ status?: unknown, recordKey: string, biologicalFeatureId: string, kinds?: readonly string[] }} UnresolvedNotice
 */

/** @returns {FeatureVisibilityRule} */
export const createDefaultFeatureVisibilityRule = () => ({
  id: nextRuleId(),
  source: 'manual',
  featureId: '',
  label: '',
  recordId: '*',
  featureType: '*',
  qualifier: 'product',
  value: '',
  action: 'off'
});

/**
 * @param {FeatureVisibilityRuleDraft} [raw]
 * @returns {FeatureVisibilityRule}
 */
export const normalizeFeatureVisibilityRule = (raw = {}) => {
  const source = normalizeSource(raw.source);
  const action = normalizeFeatureVisibilityAction(raw.action) || 'off';
  return {
    id: normalizeCell(raw.id) || nextRuleId(source === 'file' ? 'feature-visibility-file-rule' : 'feature-visibility-rule'),
    source,
    featureId: normalizeCell(raw.featureId ?? raw.feature_id),
    label: normalizeCell(raw.label),
    recordId: normalizeCell(raw.recordId ?? raw.record_id) || '*',
    featureType: normalizeCell(raw.featureType ?? raw.feature_type) || '*',
    qualifier: normalizeCell(raw.qualifier) || 'product',
    value: normalizeCell(raw.value),
    action
  };
};

/** @param {FeatureVisibilityRuleDraft} rule */
const isSerializableRule = (rule) => {
  const normalized = normalizeFeatureVisibilityRule(rule);
  return Boolean(
    normalized.recordId &&
    normalized.featureType &&
    normalized.qualifier &&
    normalized.value &&
    normalizeFeatureVisibilityAction(normalized.action)
  );
};

// The rules a request carries, in table order (Generate reads no other row).
/** @param {readonly FeatureVisibilityRuleDraft[] | null | undefined} rules */
export const requestFeatureVisibilityRules = (rules) => (Array.isArray(rules) ? rules : [])
  .map((rule) => normalizeFeatureVisibilityRule(rule))
  .filter(isSerializableRule);

/** @param {readonly FeatureVisibilityRuleDraft[] | null | undefined} rules */
export const serializeFeatureVisibilityRules = (rules) => {
  const rows = requestFeatureVisibilityRules(rules)
    .map((rule) => [
      rule.recordId,
      rule.featureType,
      rule.qualifier,
      rule.value,
      rule.action
    ].map(normalizeCell).join('\t'));

  return rows.length > 0 ? `${rows.join('\n')}\n` : '';
};

/** @param {unknown} text */
export const parseFeatureVisibilityRules = (text) => {
  /** @type {FeatureVisibilityRule[]} */
  const rules = [];
  const lines = String(text ?? '').split(/\r?\n/);

  lines.forEach((rawLine, index) => {
    const lineNo = index + 1;
    const trimmed = rawLine.trim();
    if (!trimmed || trimmed.startsWith('#')) return;

    const fields = rawLine.replace(/\r$/, '').split('\t');
    if (fields.length > REQUIRED_COLUMNS.length) {
      throw new Error(`Malformed feature visibility row at line ${lineNo}: expected ${REQUIRED_COLUMNS.length} columns.`);
    }
    if (fields.length < REQUIRED_COLUMNS.length) {
      throw new Error(`Missing feature visibility columns at line ${lineNo}.`);
    }
    if (isHeaderRow(fields)) return;

    const [recordId, featureType, qualifier, value, actionRaw] = fields.map(normalizeCell);
    const missingIndex = [recordId, featureType, qualifier, value, actionRaw].findIndex((field) => !field);
    if (missingIndex >= 0) {
      throw new Error(`Missing ${REQUIRED_COLUMNS[missingIndex]} in feature visibility row at line ${lineNo}.`);
    }
    const action = normalizeFeatureVisibilityAction(actionRaw);
    if (!action) {
      throw new Error(`Invalid feature visibility action at line ${lineNo}: ${actionRaw}`);
    }

    rules.push(normalizeFeatureVisibilityRule({
      id: `feature-visibility-file-${lineNo}`,
      source: 'file',
      recordId,
      featureType,
      qualifier,
      value,
      action
    }));
  });

  return { rules, count: rules.length };
};

/** @param {VisibilityFeature | null} [feat] */
const getFeatureId = (feat = {}) => normalizeCell(feat?.svg_id ?? feat?.svgId ?? feat?.featureId ?? feat?.feature_id ?? feat?.id);

// Feature visibility edits live in the identity-keyed featureOverrides draft
// (design Q4); these read and write one feature's `featureVisibility`.
/**
 * @param {FeatureOverrideDraft | null | undefined} featureOverrides
 * @param {VisibilityFeature | null | undefined} feature
 */
export const getFeatureVisibilityOverride = (featureOverrides, feature) => normalizeVisibilityMode(
  featureOverrideValue(featureOverrides, feature, 'featureVisibility')
);

/**
 * @param {FeatureEditsSlice} slice `drawing.slices.featureEdits`.
 * @param {VisibilityFeature | null | undefined} feature
 * @param {unknown} modeRaw
 */
export const setFeatureVisibilityOverride = (slice, feature, modeRaw) => {
  const previous = getFeatureVisibilityOverride(featureEditsOf(slice).featureOverrides, feature);
  const mode = normalizeVisibilityMode(modeRaw);
  updateFeatureOverride(slice, feature, { featureVisibility: mode === 'default' ? null : mode });
  return previous;
};

// One change per feature identity: {scope, recordKey, biologicalFeatureId,
// featureId (the rendered ID it is drawn with), before, after}.
/**
 * @param {readonly VisibilityFeature[] | null | undefined} features
 * @param {unknown} modeRaw
 * @param {FeatureOverrideDraft} [featureOverrides]
 */
export const buildFeatureVisibilityChanges = (features, modeRaw, featureOverrides = {}) => {
  const mode = normalizeVisibilityMode(modeRaw);
  /** @type {Set<string>} */
  const seen = new Set();
  /** @type {FeatureVisibilityChange[]} */
  const changes = [];
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const key = featureIdentityKeyOf(feature);
    if (!key || seen.has(key)) return;
    seen.add(key);
    const before = getFeatureVisibilityOverride(featureOverrides, feature);
    if (before === mode) return;
    changes.push({ ...parseFeatureIdentityKey(key), featureId: getFeatureId(feature), before, after: mode });
  });
  return changes;
};

/**
 * @param {FeatureEditsSlice} slice `drawing.slices.featureEdits`.
 * @param {readonly FeatureVisibilityChange[] | null | undefined} changes
 */
export const applyFeatureVisibilityOverrideChanges = (slice, changes) => {
  if (!Array.isArray(changes)) return 0;
  let applied = 0;
  changes.forEach((change) => {
    if (!featureIdentityKeyOf(change)) return;
    const mode = Object.prototype.hasOwnProperty.call(change || {}, 'mode') ? change.mode : change?.after;
    setFeatureVisibilityOverride(slice, change, mode);
    applied += 1;
  });
  return applied;
};

const KIND_FIELDS = Object.freeze({
  feature_visibility: 'featureVisibility', label_visibility: 'labelVisibility', label_text: 'labelText'
});

// Python's `unresolved` notices of a request by identity key: the edit kinds
// whose feature the source does not have (design Q4 3.4).
/**
 * @param {readonly UnresolvedNotice[] | null | undefined} notices
 * @param {(recordKey: string) => boolean} [recordKeyFilter]
 */
const unresolvedNoticeKinds = (notices, recordKeyFilter = () => true) => {
  /** @type {Map<string, Set<string>>} */
  const unresolved = new Map();
  (Array.isArray(notices) ? notices : []).forEach((notice) => {
    if (notice?.status !== 'unresolved' || !recordKeyFilter(notice.recordKey)) return;
    const key = featureIdentityKey(notice.recordKey, notice.biologicalFeatureId);
    if (key) unresolved.set(key, new Set([...(unresolved.get(key) || []), ...(notice.kinds || [])]));
  });
  return unresolved;
};

// Removes edits from the identity-keyed drafts of a drawing's feature-edits
// slice: every edit of a dropped record, the unresolved edit kinds of an
// identity, and (`sourceGone`) the label source text kept for a bulk edit.
// Returns the count of removed edits.
/**
 * @typedef {object} FeatureEditRemoval
 * @property {FeatureEditsSlice} featureEdits `drawing.slices.featureEdits`.
 * @property {Map<string, Set<string>>} [unresolved]
 * @property {(row: { readonly recordKey: string }) => boolean} [dropped]
 * @property {(key: string, row: { readonly recordKey: string }) => boolean} [sourceGone]
 */

/** @param {FeatureEditRemoval} options */
const removeFeatureEdits = ({
  featureEdits,
  unresolved = new Map(),
  dropped = () => false,
  sourceGone = () => false
}) => {
  const { featureOverrides, featurePlacementOverrides } = featureEditsOf(featureEdits);
  let removed = 0;
  Object.entries(featurePlacementOverrides).forEach(([key, row]) => {
    if (!dropped(row) && !unresolved.get(key)?.has('placement')) return;
    setFeaturePlacementRow(featureEdits, key, null);
    removed += 1;
  });
  Object.entries(featureOverrides).forEach(([key, row]) => {
    const edits = FEATURE_OVERRIDE_EDIT_FIELDS.filter((field) => row[field] !== null && row[field] !== undefined);
    if (dropped(row)) {
      removed += edits.length;
      setFeatureOverrideRow(featureEdits, key, null);
      return;
    }
    const kinds = unresolved.get(key);
    const gone = kinds ? edits.filter((field) => Object.entries(KIND_FIELDS)
      .some(([kind, name]) => name === field && kinds.has(kind))) : [];
    const source = sourceGone(key, row);
    if (gone.length === 0 && !source) return;
    removed += gone.length;
    const patch = Object.fromEntries(gone.map((field) => [field, null]));
    if (source) patch.labelSourceText = null;
    updateFeatureOverride(featureEdits, row, patch);
  });
  return removed;
};

// Record key moves (DESIGN 2-3, OV-1949): a Generate moves the per-feature
// rows of each record whose key changes (`recordKeyMoves` in
// services/record-draw-selection.js) to its new key. The moves apply together
// from one snapshot. Where rows meet under one key, the row of the earlier
// move wins field by field over a later one (an alias), and a row that stays
// at its key comes last; a field the winner leaves null keeps the other value.

/**
 * @param {unknown} winner
 * @param {unknown} other
 */
const mergeRowFields = (winner, other) => {
  if (!winner || typeof winner !== 'object' || !other || typeof other !== 'object') return winner ?? other;
  const merged = { .../** @type {Record<string, unknown>} */ (other) };
  Object.entries(/** @type {Record<string, unknown>} */ (winner)).forEach(([field, value]) => {
    if (value !== null && value !== undefined) merged[field] = value;
  });
  return merged;
};

/**
 * The rows of `rows` under their moved keys.
 * @template T
 * @param {Readonly<Record<string, T>>} rows
 * @param {ReadonlyMap<string, string>} moves
 * @param {(key: string) => { recordKey: string, biologicalFeatureId: string } | null} identityOf
 * @param {(recordKey: string, biologicalFeatureId: string, row: T) => [string, T]} rekey
 * @returns {Record<string, T>}
 */
const movedRows = (rows, moves, identityOf, rekey) => {
  if (moves.size === 0) return { ...rows };
  const priorities = new Map([...moves.keys()].map((key, index) => [key, index]));
  /** @type {Map<string, { priority: number, row: T }[]>} */
  const byKey = new Map();
  Object.entries(rows).forEach(([key, row]) => {
    const identity = identityOf(key);
    const to = identity ? moves.get(identity.recordKey) : undefined;
    const [nextKey, nextRow] = identity && to !== undefined && to !== identity.recordKey
      ? rekey(to, identity.biologicalFeatureId, row) : [key, row];
    const priority = identity && to !== undefined ? Number(priorities.get(identity.recordKey)) : Infinity;
    byKey.set(nextKey, [...(byKey.get(nextKey) || []), { priority, row: nextRow }]);
  });
  return Object.fromEntries([...byKey].map(([key, entries]) => [key, /** @type {T} */ (entries
    .sort((left, right) => left.priority - right.priority)
    .map(({ row }) => row)
    .reduce((winner, other) => /** @type {T} */ (mergeRowFields(winner, other))))]));
};

/** @param {string} key A draft key (`featureIdentityKey`). */
const draftIdentityOf = (key) => parseFeatureIdentityKey(key);
/** @param {string} key A stroke or fill key (`featureOverrideKey`): `<recordKey>\0<biologicalFeatureId>`. */
const strokeIdentityOf = (key) => {
  const at = key.indexOf('\0');
  return at > 0 ? { recordKey: key.slice(0, at), biologicalFeatureId: key.slice(at + 1) } : null;
};

/**
 * @template {{ recordKey: string }} T
 * @param {string} recordKey
 * @param {string} biologicalFeatureId
 * @param {T} row
 * @returns {[string, T]}
 */
const rekeyDraftRow = (recordKey, biologicalFeatureId, row) => [
  featureIdentityKey(recordKey, biologicalFeatureId), { ...row, recordKey }
];
/**
 * @template T
 * @param {string} recordKey
 * @param {string} biologicalFeatureId
 * @param {T} value
 * @returns {[string, T]}
 */
const rekeyStroke = (recordKey, biologicalFeatureId, value) => [
  featureOverrideKey({ recordKey, biologicalFeatureId }), value
];

/**
 * The identity-keyed values of `values` (`bulkLabelText` of a label
 * projection) under their moved keys.
 * @template T
 * @param {Readonly<Record<string, T>>} values
 * @param {ReadonlyMap<string, string>} moves
 */
export const movedIdentityValues = (values, moves) => movedRows(values, moves, draftIdentityOf,
  (recordKey, biologicalFeatureId, value) => [featureIdentityKey(recordKey, biologicalFeatureId), value]);

/**
 * The per-feature edits a request of this Generate reads: the drawing's rows
 * with the record key moves applied, without writing the drawing (the
 * drawing takes the moves when the Generate succeeds,
 * `pruneUnmatchedFeatureOverrides`). Feature fills are a cache the Generate
 * rebuilds from these.
 * @template {{ readonly recordKey: string }} O
 * @template {{ readonly recordKey: string }} P
 * @template S, F
 * @param {{ featureOverrides: Readonly<Record<string, O>>, featurePlacementOverrides: Readonly<Record<string, P>>,
 *   featureStrokeOverrides: Readonly<Record<string, S>>, featureColorOverrides: Readonly<Record<string, F>>,
 *   annotationSets: readonly AnnotationSetRow[] }} edits
 * @param {ReadonlyMap<string, string>} moves
 */
export const movedRecordEdits = (edits, moves) => ({
  featureOverrides: movedRows(edits.featureOverrides, moves, draftIdentityOf, rekeyDraftRow),
  featurePlacementOverrides: movedRows(edits.featurePlacementOverrides, moves, draftIdentityOf, rekeyDraftRow),
  // Copies: the request reads them, the drawing keeps its own.
  featureStrokeOverrides: movedRows(cloneJsonData({ ...edits.featureStrokeOverrides }), moves, strokeIdentityOf, rekeyStroke),
  featureColorOverrides: movedRows(cloneJsonData({ ...edits.featureColorOverrides }), moves, strokeIdentityOf, rekeyStroke),
  annotationSets: moves.size === 0 ? edits.annotationSets : edits.annotationSets.map((set) => ({
    ...set,
    annotations: set.annotations.map((item) => {
      const target = movedAnnotationTarget(item.target, moves);
      return target ? { ...item, target } : item;
    })
  }))
});

/**
 * A selected-feature annotation target under its record's moved key, or null
 * when it does not move.
 * @param {AnnotationRow['target']} target
 * @param {ReadonlyMap<string, string>} moves
 * @returns {AnnotationTarget | null}
 */
const movedAnnotationTarget = (target, moves) => {
  const to = target.kind === 'featureIdentity' && typeof target.recordKey === 'string' ? moves.get(target.recordKey) : undefined;
  // A copy of a readonly target with one field replaced: the same target type.
  return to === undefined ? null : /** @type {AnnotationTarget} */ ({ ...target, recordKey: to });
};

/**
 * Writes `after` over `before` key by key through `write` (null removes a key).
 * @template T
 * @param {Readonly<Record<string, T>>} before
 * @param {Readonly<Record<string, T>>} after
 * @param {(key: string, row: T | null) => void} write
 */
const replaceRows = (before, after, write) => {
  Object.keys(before).filter((key) => !Object.hasOwn(after, key)).forEach((key) => write(key, null));
  Object.entries(after).filter(([key, row]) => before[key] !== row).forEach(([key, row]) => write(key, row));
};

// Owner decision Q3 = A (design Q4 3.4, 6.3). A successful Generate that
// replaced a source removes the edits Python reported `unresolved` for a
// replaced record (Feature placement rows included), and every edit of a
// record that the previous request of this mode had and this request has not.
// A label source text kept for a bulk edit goes with its feature, and so does
// a per-feature stroke (OV-84, Owner-delegated 2026-10-07): a stroke of a
// feature the replaced source still has stays. Edits of features outside the
// crop or display stay dormant; the drafts are those of the request's drawing,
// so the other mode's edits are not read (R2).
// Then every successful Generate moves the rows of the records whose key
// changed (`moves`, `recordKeyMoves`; DESIGN A5): the removal reads the keys
// before the move, so a row whose record is drawn under another key is not
// gone, and a record that no request draws keeps its rows under its identity.
// Returns the count of removed edits.
/**
 * @param {{
 *   featureEdits: FeatureEditsSlice,
 *   paint?: PaintSlice,
 *   annotationSets?: AnnotationSetsSlice,
 *   annotationSetRows?: readonly AnnotationSetRow[],
 *   notices?: UnresolvedNotice[],
 *   replacedRecordKeys?: string[],
 *   previousRecords?: FeatureRequestRecord[],
 *   currentRecords?: FeatureRequestRecord[],
 *   retainedRecords?: FeatureRequestRecord[],
 *   biologicalFeatures?: Record<string, any>[],
 *   moves?: ReadonlyMap<string, string>
 * }} options `retainedRecords`: the records of OFF records, which the request
 *   leaves out without dropping their rows (record selection D-05), with their
 *   identities (`offRecordRequestRecords`).
 */
export const pruneUnmatchedFeatureOverrides = ({
  featureEdits,
  paint,
  annotationSets,
  annotationSetRows = [],
  notices = [],
  replacedRecordKeys = [],
  previousRecords = [],
  currentRecords = [],
  retainedRecords = [],
  biologicalFeatures = [],
  moves = new Map()
}) => {
  /** @param {string} recordKey */
  const movedKey = (recordKey) => moves.get(recordKey) ?? recordKey;
  /** @param {string} recordKey */
  const isRetained = (recordKey) => recordKeyBelongsToRequest(recordKey, retainedRecords)
    || recordKeyBelongsToRequest(movedKey(recordKey), retainedRecords);
  const replaced = new Set(replacedRecordKeys.filter((recordKey) => !isRetained(recordKey)));
  // The keys before the move of each key the new catalog and Python's notices name.
  /** @type {Map<string, string[]>} */
  const sourcesOf = new Map();
  moves.forEach((to, from) => sourcesOf.set(to, [...(sourcesOf.get(to) || []), from]));
  /** @param {string} recordKey */
  const keysBeforeMove = (recordKey) => [recordKey, ...(sourcesOf.get(recordKey) || [])];
  const presentFeatures = Array.isArray(biologicalFeatures) ? biologicalFeatures : [];
  const presentIdentities = presentFeatures.flatMap((feature) => {
    const recordKey = String(feature?.record_key ?? feature?.recordKey ?? '');
    const biologicalFeatureId = String(feature?.biological_feature_id ?? feature?.biologicalFeatureId ?? '');
    return keysBeforeMove(recordKey).map((key) => ({ recordKey: key, biologicalFeatureId }));
  });
  const present = new Set(presentIdentities.map(featureIdentityKeyOf).filter(Boolean));
  const presentStrokeKeys = new Set(presentIdentities.map(stableKeyOf).filter(Boolean));
  let removedStrokes = 0;
  Object.keys(paint ? paintOf(paint).featureStrokeOverrides : {}).forEach((key) => {
    const recordKey = key.split('\0')[0];
    if (!paint || !key.includes('\0') || !replaced.has(recordKey) || presentStrokeKeys.has(key)) return;
    setFeatureStroke(paint, key, null);
    removedStrokes += 1;
  });
  const removed = removedStrokes + removeFeatureEdits({
    featureEdits,
    unresolved: unresolvedNoticeKinds(
      (Array.isArray(notices) ? notices : []).flatMap((notice) => keysBeforeMove(String(notice?.recordKey ?? ''))
        .map((recordKey) => ({ ...notice, recordKey }))),
      (recordKey) => replaced.has(recordKey)
    ),
    dropped: (row) => rowBelongsToRequest(row, previousRecords)
      && !recordKeyBelongsToRequest(movedKey(row.recordKey), currentRecords)
      && !isRetained(row.recordKey),
    sourceGone: (key, row) => replaced.has(row?.recordKey) && !present.has(key)
  });
  if (moves.size === 0) return removed;
  const { featureOverrides, featurePlacementOverrides } = featureEditsOf(featureEdits);
  replaceRows(featureOverrides, movedRows(featureOverrides, moves, draftIdentityOf, rekeyDraftRow),
    (key, row) => setFeatureOverrideRow(featureEdits, key, row));
  replaceRows(featurePlacementOverrides, movedRows(featurePlacementOverrides, moves, draftIdentityOf, rekeyDraftRow),
    (key, row) => setFeaturePlacementRow(featureEdits, key, row));
  if (paint) {
    const strokes = paintOf(paint).featureStrokeOverrides;
    replaceRows(strokes, movedRows(strokes, moves, strokeIdentityOf, rekeyStroke),
      (key, stroke) => setFeatureStroke(paint, key, stroke));
  }
  annotationSetRows.forEach((set) => set.annotations.forEach((item) => {
    const target = movedAnnotationTarget(item.target, moves);
    if (annotationSets && target) replaceAnnotationTarget(annotationSets, set, item, target);
  }));
  return removed;
};

// The unresolved edits the drafts still hold after a Generate that replaced no
// source; "Remove N unmatched feature edits" removes them (an explicit Reset,
// R2). The drafts are those of the drawing of the request whose notices these are.
/** @typedef {{ featureEdits: FeatureEditsSlice, notices?: readonly UnresolvedNotice[] }} UnresolvedEditDrafts */

/** @param {UnresolvedEditDrafts} drafts */
export const countUnresolvedFeatureEdits = ({ featureEdits, notices = [] }) => {
  const { featureOverrides, featurePlacementOverrides } = featureEditsOf(featureEdits);
  let count = 0;
  unresolvedNoticeKinds(notices).forEach((kinds, key) => {
    if (kinds.has('placement') && featurePlacementOverrides[key]) count += 1;
    const row = featureOverrides[key];
    Object.entries(KIND_FIELDS).forEach(([kind, field]) => {
      if (kinds.has(kind) && row?.[field] !== null && row?.[field] !== undefined) count += 1;
    });
  });
  return count;
};

/** @param {UnresolvedEditDrafts} drafts */
export const removeUnresolvedFeatureEdits = ({ featureEdits, notices = [] }) => (
  removeFeatureEdits({ featureEdits, unresolved: unresolvedNoticeKinds(notices) })
);

/** @param {readonly FeatureVisibilityRuleDraft[] | null | undefined} rules */
export const splitLegacyVisibilityRules = (rules) => {
  /** @type {Record<string, string>} */
  const overrides = {};
  /** @type {FeatureVisibilityRule[]} */
  const manualRules = [];
  /** @type {string[]} */
  const warnings = [];
  (Array.isArray(rules) ? rules : []).forEach((rule) => {
    const normalized = normalizeFeatureVisibilityRule(rule);
    const mode = featureVisibilityActionToMode(normalized.action);
    if (normalized.source === 'editor' && normalized.featureId) {
      if (mode === 'default') {
        warnings.push(`Feature visibility rule ${normalized.id} has no supported action and was kept as manual.`);
        manualRules.push(normalized);
        return;
      }
      overrides[normalized.featureId] = mode;
      return;
    }
    manualRules.push(normalized);
  });
  return { overrides, manualRules, warnings };
};

/**
 * The inputs `resolveFeatureDrawn` reads.
 * @typedef {object} FeatureDrawnContext
 * @property {FeatureOverrideDraft} featureOverrides
 * @property {FeatureVisibilityRule[]} rules
 * @property {Set<string> | null} selectedTypes `null` while the request's types are unknown.
 * @property {readonly Partial<SpecificColorRule>[]} colorRules
 */

// The inputs `resolveFeatureDrawn` reads. The feature types come from the
// request's diagram options (a label rerender keeps the last Generate's); the
// edits and rules are the current ones, which a label rerender and Generate
// both carry.
/**
 * @param {Record<string, any>} drawing The drawing whose edits and rules decide.
 * @param {{ diagramOptions?: { selectedFeaturesSet?: unknown } | null, featureOverrides?: FeatureOverrideDraft }} [inputs]
 * @returns {FeatureDrawnContext}
 */
export const featureDrawnContext = (drawing, {
  diagramOptions = /** @type {{ selectedFeaturesSet?: unknown } | null} */ (null),
  featureOverrides = drawing?.featureOverrides
} = {}) => ({
  featureOverrides: featureOverrides || {},
  rules: requestFeatureVisibilityRules(drawing?.featureVisibilityManualRules),
  selectedTypes: Array.isArray(diagramOptions?.selectedFeaturesSet)
    ? new Set(diagramOptions.selectedFeaturesSet.map(String))
    : null,
  colorRules: Array.isArray(drawing?.manualSpecificRules) ? drawing.manualSpecificRules : []
});
// The color rule matches of a context's rules, built once per context: a
// context is built for one pass over the features (OV-193).
/** @type {WeakMap<object, RuleMatcher>} */
const colorMatchersByContext = new WeakMap();
/** @param {FeatureDrawnContext} context */
const colorRuleMatcher = (context) => {
  let matcher = colorMatchersByContext.get(context);
  if (!matcher) colorMatchersByContext.set(context, matcher = ruleMatcher(context.colorRules));
  return matcher;
};

// "Is this feature drawn?" as Generate answers it
// (gbdraw/features/visibility.py::should_render_feature): the feature's
// identity-keyed override decides when On or Off, and Leave out of protein
// search skips the rules; else the first matching visibility rule decides when Show
// or Off; else the selected feature types, where a feature of another type is
// drawn when a specific color rule matches it. Python matches the rules
// (app/rule-matching.js, R4); the shared vectors in
// tests/fixtures/feature_drawn_cases.json hold both sides to one answer.
// Returns true, false, or null while a match it needs is unknown.
/**
 * @param {VisibilityFeature | null | undefined} feature
 * @param {FeatureDrawnContext} context
 * @returns {boolean | null}
 */
export const resolveFeatureDrawn = (feature, context) => {
  const { featureOverrides, rules, selectedTypes } = context;
  if (!feature) return null;
  const override = getFeatureVisibilityOverride(featureOverrides, feature);
  if (override === 'on' || override === 'off') return override === 'on';
  if (override !== 'exclude_matching') {
    for (const rule of rules) {
      const matched = visibilityRuleMatchesFeature(feature, rule);
      if (matched === null) return null;
      if (!matched) continue;
      if (rule.action === 'show' || rule.action === 'off') return rule.action === 'show';
      break;
    }
  }
  if (!selectedTypes) return null;
  if (selectedTypes.size === 0 || selectedTypes.has(String(feature?.type ?? ''))) return true;
  return colorRuleMatcher(context).matchesAny(feature);
};

/**
 * What the Legend readers take of the displayed Result's committed metadata (svg-result-ingestion.js
 * `CommittedSvgResultMetadata`).
 * @typedef {{ renderedFeatureIdentities?: { renderedIds?: unknown } | null, legendRowFacts?: SvgLegendRowFacts | null }}
 *   DisplayedResultLegendMetadata
 */

// The features the displayed Result draws, as its Legend rows are read
// (`LegendRowContext.features`, N-06): those Python drew in it at its last
// render (the rendered IDs of its committed metadata; every read feature for a
// Result without them) that are still drawn (`resolveFeatureDrawn`; an unknown
// answer keeps a feature). Read when iterated.
/**
 * @param {{ extractedFeatures?: { value?: unknown },
 *   displayedResultMetadata?: () => DisplayedResultLegendMetadata | null | undefined }} state
 * @param {Parameters<typeof featureDrawnContext>[0]} drawing
 * @returns {Iterable<{ svg_id?: unknown, type?: unknown, fill_color?: unknown }>}
 */
export const displayedDrawnFeatures = (state, drawing) => ({
  * [Symbol.iterator]() {
    const features = Array.isArray(state.extractedFeatures?.value) ? state.extractedFeatures.value : [];
    const renderedIds = state.displayedResultMetadata?.()?.renderedFeatureIdentities?.renderedIds;
    const context = featureDrawnContext(drawing);
    for (const feature of features) {
      if (renderedIds instanceof Set && !renderedIds.has(feature?.svg_id)) continue;
      if (resolveFeatureDrawn(feature, context) !== false) yield feature;
    }
  }
});

/**
 * What `displayedLegendRowContext` reads of the app state.
 * @typedef {Parameters<typeof displayedDrawnFeatures>[0] & {
 *   svgContainer?: { value?: { querySelector?: (selector: 'svg') => SVGSVGElement | null } | null },
 *   originalLegendOrder?: { value?: readonly string[] }
 * }} DisplayedLegendState
 */

/** @type {ReadonlySet<PythonLegendKey>} */
const NO_LEGEND_KEYS = new Set();

/**
 * The Legend state of the displayed Result, and Python's keys of it that no row draws (its Legend row
 * facts' `suppressed`, which hold the shared keys `trackLegendRows` reads; none for a Result without facts).
 * @typedef {LegendRowContext & { sharedLegendKeys: ReadonlySet<PythonLegendKey> }} DisplayedLegendRowContext
 */

// The Legend state of the displayed Result that every caller of the N-06
// allocation (specific-color-rules.js) and of the track rows (track-legend-rows.js) reads: the compile,
// the Legend panel, the feature popup, the Legend editor and the rule commit.
/**
 * @param {DisplayedLegendState} state
 * @param {Parameters<typeof displayedDrawnFeatures>[1] & { manualSpecificRules: readonly Partial<SpecificColorRule>[] }} drawing
 * @returns {DisplayedLegendRowContext}
 */
export const displayedLegendRowContext = (state, drawing) => ({
  rules: drawing.manualSpecificRules,
  pythonRows: pythonLegendRows(state.svgContainer?.value?.querySelector?.('svg')),
  features: displayedDrawnFeatures(state, drawing),
  originalLegendOrder: state.originalLegendOrder?.value || [],
  sharedLegendKeys: /** @type {ReadonlySet<PythonLegendKey> | undefined} */ (state.displayedResultMetadata?.()?.legendRowFacts?.suppressed)
    ?? NO_LEGEND_KEYS
});

// Whether a Result draws a catalog feature (R-5): the resolver's answer, else
// (unknown) whether the Result's catalog (`resultCatalogFeatures`) renders it,
// which Python decided when it last rendered the Result; null without one.
/**
 * @param {VisibilityFeature | null | undefined} feature
 * @param {FeatureDrawnContext} context
 * @param {LegendCatalogFeatures | null | undefined} catalogFeatures
 * @returns {boolean | null}
 */
export const featureDrawnInResult = (feature, context, catalogFeatures) => (
  resolveFeatureDrawn(feature, context)
  ?? (catalogFeatures ? catalogFeatures.renderedByIdentity.has(stableKeyOf(feature)) : null)
);

// The Features list of one Result (Owner decision 2026-10-05, R-5): from the
// Result's catalog features, each biological feature that is of a selected
// type, has its own Feature visibility, or is drawn, with whether it is drawn.
// A feature the Result draws is listed as its rendered feature.
/**
 * @param {LegendCatalogFeatures} catalogFeatures
 * @param {FeatureDrawnContext} context
 */
export const listFeatureRows = (catalogFeatures, context) => {
  /** @type {VisibilityFeature[]} */
  const rows = [];
  /** @type {Map<VisibilityFeature, boolean | null>} */
  const drawn = new Map();
  /** @type {Set<VisibilityFeature>} */
  const rendered = new Set();
  const { selectedTypes, featureOverrides } = context;
  // The catalog's rows are untyped (`any`) at their provider; this reads the fields of `VisibilityFeature`.
  catalogFeatures.biological.forEach((/** @type {VisibilityFeature} */ feature) => {
    const shown = catalogFeatures.renderedByIdentity.get(stableKeyOf(feature));
    const row = shown || feature;
    const isDrawn = featureDrawnInResult(row, context, catalogFeatures);
    if (!isDrawn
      && selectedTypes?.size && !selectedTypes.has(String(feature.type ?? ''))
      && getFeatureVisibilityOverride(featureOverrides, feature) === 'default') return;
    rows.push(row);
    drawn.set(row, isDrawn);
    if (shown) rendered.add(row);
  });
  return { rows, drawn, rendered };
};

// The features a Result draws, in catalog order, each with the color rule
// that wins it when its type has a captioned rule (`undefined` for a type
// without one: Python keeps its default row whatever its rules match, so no
// match is read). Null while a match it needs is unknown. A Result is read
// from its catalog features: `asRendered` reads the features Python drew at
// the last render, else `featureDrawnInResult` answers.
/** @typedef {NonNullable<ReturnType<typeof resultCatalogFeatures>>} LegendCatalogFeatures */
/** @typedef {{ feat?: unknown, cap?: unknown, color?: unknown }} LegendRuleWinner The color rule that wins a feature (a row of `colorRules`). */
/**
 * @param {LegendCatalogFeatures} catalogFeatures
 * @param {FeatureDrawnContext} context
 * @param {{ asRendered: boolean }} options
 */
const drawnLegendFeatures = (catalogFeatures, context, { asRendered }) => {
  const rules = context.colorRules;
  const ruleMatches = colorRuleMatcher(context);
  const captionedTypes = new Set(rules.filter((rule) => legendRuleCaption(rule)).map((rule) => rule.feat));
  const drawn = [];
  for (const feature of catalogFeatures.biological) {
    const shown = catalogFeatures.renderedByIdentity.get(stableKeyOf(feature));
    const row = shown || feature;
    if (!(asRendered ? shown : featureDrawnInResult(row, context, catalogFeatures))) continue;
    const type = String(feature.type ?? '');
    if (!captionedTypes.has(type) && !captionedTypes.has('*')) {
      drawn.push({ feature, type, winner: undefined });
      continue;
    }
    const winner = ruleMatches.firstIfKnown(row);
    if (winner === undefined) return null;
    drawn.push({ feature, type, winner });
  }
  return drawn;
};
/** @param {{ cap?: unknown } | null | undefined} rule */
const legendRuleCaption = (rule) => String(rule?.cap ?? '').trim();

// What Python derives the feature rows of a Result's Legend from
// (gbdraw/legend/table.py::prepare_legend_table): per record, the types of
// the drawn features in first-drawn order; for a type with a captioned rule,
// the rules whose caption a drawn feature uses (in table order) and whether a
// drawn feature keeps the default color (the "other <type>s" row); and every
// rule caption, which Python reserves. A type without a captioned rule keeps
// its default row whatever the rules match, so no match is read for it, and
// neither is one for a type that its rules, all of the caption `<type>` and one
// color, recolor as a whole: that row keeps its caption and the Legend follows
// the color live. This is the input of the rows, not the rows: Python stays
// their only derivation (Owner decision 2026-10-06, option A), and an edit
// that changes a Result's source asks for the automatic rerender (OV-42).
// A rule is used by its caption and its color, as Python draws a row for each
// (caption, color) pair a drawn feature uses: hiding the only feature of one
// color of a caption written in several colors removes that row (OV-1236).
// A source is null while a rule match it needs is unknown; null equals no
// source (`sameLegendSources`).
/**
 * @param {LegendCatalogFeatures} catalogFeatures
 * @param {FeatureDrawnContext} context
 * @param {{ asRendered: boolean }} options
 */
const legendSourceOf = (catalogFeatures, context, { asRendered }) => {
  const caption = legendRuleCaption;
  /** @param {{ color?: unknown } | null | undefined} rule */
  const colorOf = (rule) => String(rule?.color ?? '').trim().toLowerCase();
  /** @param {{ cap?: unknown, color?: unknown } | null | undefined} rule */
  const usage = (rule) => JSON.stringify([caption(rule), colorOf(rule)]);
  const rules = context.colorRules;
  const drawn = drawnLegendFeatures(catalogFeatures, context, { asRendered });
  if (!drawn) return null;
  /** @type {Map<string, string[]>} */
  const typesByRecord = new Map();
  /** @type {Map<string, (LegendRuleWinner | null)[]>} */
  const winnersByType = new Map();
  for (const { feature, type, winner } of drawn) {
    const record = String(feature.recordKey ?? '');
    const recordTypes = typesByRecord.get(record) || [];
    typesByRecord.set(record, recordTypes);
    if (!recordTypes.includes(type)) recordTypes.push(type);
    const winners = winnersByType.get(type) || [];
    winnersByType.set(type, winners);
    if (winner !== undefined) winners.push(winner);
  }
  /** @type {Set<string>} */
  const used = new Set();
  /** @type {Set<string>} */
  const defaultTypes = new Set();
  /** @type {Set<string>} */
  const plainTypes = new Set();
  /** @type {Set<string | undefined>} */
  const recolored = new Set();
  for (const [type, winners] of winnersByType) {
    if (winners.length === 0) {
      plainTypes.add(type);
    } else if (winners.every((winner) => winner && caption(winner) === type && winner.feat === type)
      && new Set(winners.map(colorOf)).size === 1) {
      plainTypes.add(type);
      recolored.add(type);
    } else {
      winners.forEach((winner) => {
        if (!winner) defaultTypes.add(type);
        else if (caption(winner)) used.add(usage(winner));
      });
    }
  }
  /** @param {Partial<SpecificColorRule>} rule */
  const counted = (rule) => !(recolored.has(rule.feat) && caption(rule) === rule.feat);
  return JSON.stringify([
    [...typesByRecord],
    rules.filter((rule) => caption(rule) && counted(rule) && used.has(usage(rule)))
      .map((rule) => [String(rule.feat ?? ''), caption(rule), colorOf(rule)]),
    [...defaultTypes].sort(),
    [...plainTypes].sort(),
    [...new Set(rules.filter(counted).map(caption).filter(Boolean))].sort()
  ]);
};

// The Legend source of every committed Result, in Result order; `mode`: the
// mode of the drawing that owns the Results (`resultCatalogFeatures`).
/**
 * @param {{ results?: { value?: unknown } } | null | undefined} state
 * @param {string} mode
 * @param {FeatureDrawnContext} context
 * @param {{ asRendered?: boolean }} [options]
 */
export const resultLegendSources = (state, mode, context, { asRendered = false } = {}) => {
  const results = Array.isArray(state?.results?.value) ? state.results.value : [];
  return results.map((_, index) => {
    const catalogFeatures = resultCatalogFeatures(state, mode, index);
    return catalogFeatures ? legendSourceOf(catalogFeatures, context, { asRendered }) : '';
  });
};

/** @param {readonly (string | null)[]} left @param {readonly (string | null)[]} right */
export const sameLegendSources = (left, right) => left.length === right.length
  && left.every((source, index) => source !== null && source === right[index]);

// U3a 1d (R14-8): the Legend row of each feature a Result draws, as Python
// keys it (prepare_legend_table): the key the whole rule table allocates to
// the row of the captioned rule that wins it (`<caption> [<fill>]` for a
// caption written in several colors, OV-1235); else the "other <type>s" row of
// a type whose drawn features use a captioned rule, or no row when an
// uncaptioned rule colors it; else its type's row. This tells which rows a
// rule change regroups; Python still derives the rows and their order. Null
// while a match is unknown, when a wildcard rule has a caption, or when a used
// caption names a row Python generates (it suffixes the caption by color,
// PD-OI-042). Without `byColor`, a rule row is keyed by its caption (a History
// restore, D-15-38: a restored table whose caption is written in several
// colors shows its rows live, accepted residual OV-357).
/**
 * @param {LegendCatalogFeatures} catalogFeatures @param {FeatureDrawnContext} context
 * @param {{ byColor: boolean }} options
 */
const legendRowKeysOf = (catalogFeatures, context, { byColor }) => {
  const drawn = drawnLegendFeatures(catalogFeatures, context, { asRendered: false });
  if (!drawn || context.colorRules.some((rule) => rule?.feat === '*' && legendRuleCaption(rule))) return null;
  const ruleTypes = new Set(drawn.filter(({ winner }) => legendRuleCaption(winner)).map(({ type }) => type));
  const ruleRowKey = byColor ? createRuleLegendCaptions(context.colorRules) : legendRuleCaption;
  /** @type {Map<string, string | null>} */
  const keys = new Map();
  const generated = new Set();
  for (const { feature, type, winner } of drawn) {
    const caption = legendRuleCaption(winner);
    /** @type {string | null} */
    let key = caption ? ruleRowKey(winner) : null;
    if (!caption && !ruleTypes.has(type)) key = type;
    else if (!caption && !winner) key = type === 'CDS' ? 'other proteins' : `other ${type}s`;
    if (key && !caption) generated.add(key);
    keys.set(stableKeyOf(feature), key);
  }
  return drawn.some(({ winner }) => generated.has(legendRuleCaption(winner))) ? null : keys;
};

// The row keys of every committed Result, in Result order (an empty map for
// a Result without a catalog); `byColor`: as `legendRowKeysOf`; `mode`: as
// `resultLegendSources`.
/**
 * @param {{ results?: { value?: unknown } }} state @param {string} mode @param {FeatureDrawnContext} context
 * @param {{ byColor?: boolean }} [options]
 */
export const resultLegendRowKeys = (state, mode, context, { byColor = true } = {}) => {
  const results = Array.isArray(state?.results?.value) ? state.results.value : [];
  return results.map((_, index) => {
    const catalogFeatures = resultCatalogFeatures(state, mode, index);
    return catalogFeatures ? legendRowKeysOf(catalogFeatures, context, { byColor }) : new Map();
  });
};

// The place of each rule row in Python's Legend: the index of the first rule
// drawn in the row of the key, by the key the rule table allocates to its
// rule rows (`legendRowKeysOf`; prepare_legend_table walks the rules in table
// order and skips a row it has drawn); -1 for a row no rule names.
// `byColor`: as `legendRowKeysOf`.
/**
 * @param {readonly Partial<SpecificColorRule>[]} rules @param {{ byColor?: boolean }} [options]
 * @returns {(key: string) => number}
 */
export const legendRuleOrder = (rules, { byColor = true } = {}) => {
  const ruleRowKey = byColor ? createRuleLegendCaptions(rules) : legendRuleCaption;
  return (key) => rules.findIndex((rule) => Boolean(legendRuleCaption(rule)) && ruleRowKey(rule) === key);
};

// Whether the displayed Result shows the Legend rows a change of the row keys
// regroups, so Python need not draw them (R14-8). Every drawn feature stays
// drawn and each row keeps its features together: a row relabeled to a new
// key (`shows` the new key and not the old), merged into a row it joins
// (`shows` the joined row and not the merged one), or unchanged (a recolor).
// A merge keeps the joined row in its place, so a row merged in from a rule
// before the joined row's first rule (`firstRule`, the rules before the
// change) moves it in Python's Legend; an "other <type>s" row has no rule.
// Another batch Result shows only its committed rows, so its keys must not
// change. `shows(key)`: whether the displayed Result shows a row of the key
// once the change is shown; `placed(key, nextKey)`: whether the row of a new
// key takes the place of the row it relabels, as Python keeps it (OV-158).
// `restored(key)` (a History restore): whether the restored Result shows
// Python's own row of the key, so the rows a restore splits are exact when it
// shows Python's rows of every key after the change and no row of the keys
// that leave.
/**
 * @param {(Map<string, string | null> | null)[]} before
 * @param {(Map<string, string | null> | null)[]} after
 * @param {{ displayed: number, shows: (key: string) => boolean, firstRule: (key: string) => number,
 *   placed?: (key: string, nextKey: string) => boolean, restored?: ((key: string) => boolean) | null }} options
 */
export const legendRowsShowable = (before, after, { displayed, shows, firstRule, placed = () => true, restored = null }) => (
  before.length === after.length
  && before.every((keys, index) => {
    const next = after[index];
    if (!keys || !next || keys.size !== next.size) return false;
    if (![...keys.keys()].every((feature) => next.has(feature))) return false;
    const nextKeys = new Set([...next.values()].filter((key) => key !== null));
    const exact = Boolean(restored) && index === displayed && [...nextKeys].every((key) => restored?.(key))
      && [...keys.values()].every((key) => key === null || nextKeys.has(key) || !shows(key));
    /** @type {Map<string, string>} */
    const moved = new Map();
    /** @type {Map<string, Set<string>>} */
    const joined = new Map();
    for (const [feature, key] of keys) {
      const nextKey = next.get(feature) ?? null;
      if (key === null || nextKey === null) {
        if (key !== nextKey) return exact;
        continue;
      }
      if (moved.has(key) ? moved.get(key) !== nextKey : (index !== displayed && key !== nextKey)) return exact;
      moved.set(key, nextKey);
      if (!joined.has(nextKey)) joined.set(nextKey, new Set());
      joined.get(nextKey)?.add(key);
    }
    return [...joined].every(([nextKey, keysBefore]) => {
      if (keysBefore.size === 1 && keysBefore.has(nextKey)) return true;
      const [first] = keysBefore;
      const place = firstRule(nextKey);
      const target = keysBefore.size === 1
        ? !moved.has(nextKey) && placed(first, nextKey)
        : moved.get(nextKey) === nextKey && [...keysBefore].every((key) => key === nextKey || place < 0 || firstRule(key) < 0 || firstRule(key) > place);
      return target && shows(nextKey) && [...keysBefore].every((key) => key === nextKey || !shows(key));
    });
  })
);

/** @param {FeatureVisibilityRuleDraft} rule */
const isEditorFeatureRule = (rule) => {
  const normalized = normalizeFeatureVisibilityRule(rule);
  return normalized.source === 'editor' && Boolean(normalized.featureId);
};

/** @param {FeatureVisibilityRuleDraft} rule */
const isEditorExactQualifierRule = (rule) => {
  const normalized = normalizeFeatureVisibilityRule(rule);
  const qualifier = normalized.qualifier.toLowerCase();
  return normalized.source === 'editor' &&
    normalized.recordId === '*' &&
    normalized.featureType &&
    normalized.featureType !== '*' &&
    ['product', 'protein_id'].includes(qualifier) &&
    normalized.value.startsWith('^') &&
    normalized.value.endsWith('$');
};

/** @param {FeatureVisibilityRuleDraft[] | null | undefined} rules */
const reorderEditorVisibilityRules = (rules) => {
  if (!Array.isArray(rules)) return;
  /** @type {FeatureVisibilityRuleDraft[]} */
  const featureRules = [];
  /** @type {FeatureVisibilityRuleDraft[]} */
  const qualifierRules = [];
  /** @type {FeatureVisibilityRuleDraft[]} */
  const otherRules = [];
  rules.forEach((rule) => {
    if (isEditorFeatureRule(rule)) {
      featureRules.push(rule);
    } else if (isEditorExactQualifierRule(rule)) {
      qualifierRules.push(rule);
    } else {
      otherRules.push(rule);
    }
  });
  rules.splice(0, rules.length, ...featureRules, ...qualifierRules, ...otherRules);
};

/**
 * @param {FeatureVisibilityRuleInput} [input]
 * @returns {FeatureVisibilityRule | null}
 */
export const buildExactQualifierFeatureVisibilityRule = ({
  featureType,
  qualifier,
  value,
  action: actionRaw,
  label = ''
} = {}) => {
  const normalizedFeatureType = normalizeCell(featureType);
  const normalizedQualifier = normalizeCell(qualifier).toLowerCase();
  const normalizedValue = normalizeCell(value);
  const action = featureVisibilityModeToAction(actionRaw);
  if (!normalizedFeatureType || !normalizedQualifier || !normalizedValue || !action) return null;
  return normalizeFeatureVisibilityRule({
    source: 'editor',
    featureId: '',
    label: normalizeCell(label) || `${normalizedQualifier}: ${normalizedValue}`,
    recordId: '*',
    featureType: normalizedFeatureType,
    qualifier: normalizedQualifier,
    value: exactRegexValue(normalizedValue),
    action
  });
};

/**
 * @param {FeatureVisibilityRuleDraft[] | null | undefined} rules The table, edited in place.
 * @param {FeatureVisibilityRuleInput | undefined} ruleInput
 * @param {unknown} actionRaw
 */
export const upsertEditorQualifierFeatureVisibilityRule = (rules, ruleInput, actionRaw) => {
  if (!Array.isArray(rules)) return null;
  const nextRule = buildExactQualifierFeatureVisibilityRule({ ...ruleInput, action: actionRaw });
  if (!nextRule) {
    removeEditorQualifierFeatureVisibilityRule(rules, ruleInput);
    return null;
  }
  const existingIndex = rules.findIndex((rule) => {
    const normalized = normalizeFeatureVisibilityRule(rule);
    return normalized.source === 'editor' &&
      normalized.recordId === nextRule.recordId &&
      normalized.featureType === nextRule.featureType &&
      normalized.qualifier.toLowerCase() === nextRule.qualifier.toLowerCase() &&
      normalized.value === nextRule.value;
  });

  if (existingIndex >= 0) {
    nextRule.id = normalizeFeatureVisibilityRule(rules[existingIndex]).id;
    rules.splice(existingIndex, 1, nextRule);
  } else {
    rules.unshift(nextRule);
  }
  reorderEditorVisibilityRules(rules);
  return nextRule;
};

/**
 * @param {FeatureVisibilityRuleDraft[] | null | undefined} rules The table, edited in place.
 * @param {FeatureVisibilityRuleInput} [ruleInput]
 */
export const removeEditorQualifierFeatureVisibilityRule = (rules, ruleInput = {}) => {
  if (!Array.isArray(rules)) return 0;
  const featureType = normalizeCell(ruleInput.featureType);
  const qualifier = normalizeCell(ruleInput.qualifier).toLowerCase();
  const value = normalizeCell(ruleInput.value);
  if (!featureType || !qualifier || !value) return 0;
  const expectedValue = exactRegexValue(value);
  let removed = 0;
  for (let index = rules.length - 1; index >= 0; index -= 1) {
    const rule = normalizeFeatureVisibilityRule(rules[index]);
    if (rule.source !== 'editor') continue;
    if (rule.recordId !== '*') continue;
    if (rule.featureType !== featureType) continue;
    if (rule.qualifier.toLowerCase() !== qualifier) continue;
    if (rule.value !== expectedValue) continue;
    rules.splice(index, 1);
    removed += 1;
  }
  return removed;
};
