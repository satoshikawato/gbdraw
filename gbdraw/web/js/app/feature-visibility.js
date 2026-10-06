// @ts-check
import { exactRegexValue } from '../services/feature-selector.js';
import { firstMatchingRuleIfKnown, ruleMatchesFeature, visibilityRuleMatchesFeature } from './rule-matching.js';
import { resultCatalogFeatures, stableFeatureOverrideKey as stableKeyOf } from '../services/feature-catalog.js';
import {
  featureIdentityKey,
  featureIdentityKeyOf,
  featureOverrideValue,
  parseFeatureIdentityKey,
  rowBelongsToRequest,
  updateFeatureOverride
} from '../services/feature-placement.js';
import { normalizeTsvCell as normalizeCell } from '../utils/tsv-cell.js';
export { escapeRegexLiteral, exactRegexValue } from '../services/feature-selector.js';

const REQUIRED_COLUMNS = ['record_id', 'feature_type', 'qualifier', 'value', 'action'];
const COMMON_QUALIFIERS = ['product', 'gene', 'protein_id', 'locus_tag', 'hash', 'location', 'record_location'];
const SHOW_ACTIONS = new Set(['show', 'on']);
const HIDE_ACTIONS = new Set(['hide', 'off', 'false', '0']);
const EXCLUDE_MATCHING_ACTIONS = new Set(['exclude_matching', 'suppress']);

let generatedRuleId = 0;

const normalizeSource = (value) => {
  const normalized = normalizeCell(value).toLowerCase();
  return ['manual', 'editor', 'file'].includes(normalized) ? normalized : 'manual';
};

export const featureVisibilityQualifierSuggestions = COMMON_QUALIFIERS;

export const normalizeVisibilityMode = (value) => {
  const normalized = normalizeCell(value).toLowerCase();
  if (normalized === 'suppress') return 'exclude_matching';
  return ['on', 'off', 'exclude_matching'].includes(normalized) ? normalized : 'default';
};

export const normalizeFeatureVisibilityAction = (value) => {
  const normalized = normalizeCell(value).toLowerCase();
  if (SHOW_ACTIONS.has(normalized)) return 'show';
  if (HIDE_ACTIONS.has(normalized)) return 'off';
  if (EXCLUDE_MATCHING_ACTIONS.has(normalized)) return 'exclude_matching';
  return '';
};

export const featureVisibilityActionToMode = (value) => {
  const action = normalizeFeatureVisibilityAction(value);
  if (action === 'show') return 'on';
  if (action === 'off') return 'off';
  if (action === 'exclude_matching') return 'exclude_matching';
  return 'default';
};

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
 * The fields an exact-qualifier rule is built or removed by.
 * @typedef {object} FeatureVisibilityRuleInput
 * @property {string} [featureType]
 * @property {string} [qualifier]
 * @property {string} [value]
 * @property {string} [action]
 * @property {string} [label]
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
 * @param {Record<string, any>} [raw]
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
export const requestFeatureVisibilityRules = (rules) => (Array.isArray(rules) ? rules : [])
  .map((rule) => normalizeFeatureVisibilityRule(rule))
  .filter(isSerializableRule);

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

export const parseFeatureVisibilityRules = (text) => {
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

const getFeatureId = (feat = {}) => normalizeCell(feat?.svg_id ?? feat?.svgId ?? feat?.featureId ?? feat?.feature_id ?? feat?.id);

// Feature visibility edits live in the identity-keyed featureOverrides draft
// (design Q4); these read and write one feature's `featureVisibility`.
export const getFeatureVisibilityOverride = (featureOverrides, feature) => normalizeVisibilityMode(
  featureOverrideValue(featureOverrides, feature, 'featureVisibility')
);

export const setFeatureVisibilityOverride = (featureOverrides, feature, modeRaw) => {
  const previous = getFeatureVisibilityOverride(featureOverrides, feature);
  const mode = normalizeVisibilityMode(modeRaw);
  updateFeatureOverride(featureOverrides, feature, { featureVisibility: mode === 'default' ? null : mode });
  return previous;
};

// One change per feature identity: {scope, recordKey, biologicalFeatureId,
// featureId (the rendered ID it is drawn with), before, after}.
export const buildFeatureVisibilityChanges = (features, modeRaw, featureOverrides = {}) => {
  const mode = normalizeVisibilityMode(modeRaw);
  const seen = new Set();
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

export const applyFeatureVisibilityOverrideChanges = (featureOverrides, changes) => {
  if (!featureOverrides || !Array.isArray(changes)) return 0;
  let applied = 0;
  changes.forEach((change) => {
    if (!featureIdentityKeyOf(change)) return;
    const mode = Object.prototype.hasOwnProperty.call(change || {}, 'mode') ? change.mode : change?.after;
    setFeatureVisibilityOverride(featureOverrides, change, mode);
    applied += 1;
  });
  return applied;
};

const EDIT_FIELDS = ['featureVisibility', 'labelVisibility', 'labelText'];
const KIND_FIELDS = Object.freeze({
  feature_visibility: 'featureVisibility', label_visibility: 'labelVisibility', label_text: 'labelText'
});

// Python's `unresolved` notices of a request of mode `scope` by identity key:
// the edit kinds whose feature the source does not have (design Q4 3.4).
/**
 * @param {any} notices
 * @param {string} scope
 * @param {(recordKey: string) => boolean} [recordKeyFilter]
 */
const unresolvedNoticeKinds = (notices, scope, recordKeyFilter = () => true) => {
  const unresolved = new Map();
  (Array.isArray(notices) ? notices : []).forEach((notice) => {
    if (notice?.status !== 'unresolved' || !recordKeyFilter(notice.recordKey)) return;
    const key = featureIdentityKey(scope, notice.recordKey, notice.biologicalFeatureId);
    if (key) unresolved.set(key, new Set([...(unresolved.get(key) || []), ...(notice.kinds || [])]));
  });
  return unresolved;
};

// Removes edits from the identity-keyed drafts: every edit of a dropped
// record, the unresolved edit kinds of an identity, and (`sourceGone`) the
// label source text kept for a bulk edit. Returns the count of removed edits.
/**
 * @typedef {object} FeatureEditRemoval
 * @property {Record<string, any>} [featureOverrides]
 * @property {Record<string, any>} [featurePlacementOverrides]
 * @property {Map<string, Set<string>>} [unresolved]
 * @property {(row: any) => boolean} [dropped]
 * @property {(key: string, row: any) => boolean} [sourceGone]
 */

/** @param {FeatureEditRemoval} options */
const removeFeatureEdits = ({
  featureOverrides = {},
  featurePlacementOverrides = {},
  unresolved = new Map(),
  dropped = () => false,
  sourceGone = () => false
}) => {
  let removed = 0;
  Object.entries(featurePlacementOverrides || {}).forEach(([key, row]) => {
    if (!dropped(row) && !unresolved.get(key)?.has('placement')) return;
    delete featurePlacementOverrides[key];
    removed += 1;
  });
  Object.entries(featureOverrides || {}).forEach(([key, row]) => {
    const edits = EDIT_FIELDS.filter((field) => row?.[field] !== null && row?.[field] !== undefined);
    if (dropped(row)) {
      removed += edits.length;
      delete featureOverrides[key];
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
    updateFeatureOverride(featureOverrides, row, patch);
  });
  return removed;
};

// Owner decision Q3 = A (design Q4 3.4, 6.3). A successful Generate that
// replaced a source removes the edits Python reported `unresolved` for a
// replaced record (Feature placement rows included), and every edit of a
// record that the previous request of this mode had and this request has not.
// A label source text kept for a bulk edit goes with its feature. Edits of
// features outside the crop or display stay dormant, and the other mode's
// edits wait for their mode (R2). Returns the count of removed edits.
export const pruneUnmatchedFeatureOverrides = ({
  featureOverrides = {},
  featurePlacementOverrides = {},
  notices = [],
  scope = '',
  replacedRecordKeys = [],
  previousRecords = [],
  currentRecords = [],
  biologicalFeatures = []
} = {}) => {
  const replaced = new Set(replacedRecordKeys);
  const present = new Set((Array.isArray(biologicalFeatures) ? biologicalFeatures : [])
    .map(featureIdentityKeyOf).filter(Boolean));
  return removeFeatureEdits({
    featureOverrides,
    featurePlacementOverrides,
    unresolved: unresolvedNoticeKinds(notices, scope, (recordKey) => replaced.has(recordKey)),
    dropped: (row) => rowBelongsToRequest(row, scope, previousRecords)
      && !rowBelongsToRequest(row, scope, currentRecords),
    sourceGone: (key, row) => row?.scope === scope && replaced.has(row?.recordKey) && !present.has(key)
  });
};

// The unresolved edits the drafts still hold after a Generate that replaced no
// source; "Remove N unmatched feature edits" removes them (an explicit Reset,
// R2). `scope` is the mode of the request whose notices these are.
export const countUnresolvedFeatureEdits = ({
  featureOverrides = {}, featurePlacementOverrides = {}, notices = [], scope = ''
} = {}) => {
  let count = 0;
  unresolvedNoticeKinds(notices, scope).forEach((kinds, key) => {
    if (kinds.has('placement') && featurePlacementOverrides?.[key]) count += 1;
    const row = featureOverrides?.[key];
    Object.entries(KIND_FIELDS).forEach(([kind, field]) => {
      if (kinds.has(kind) && row?.[field] !== null && row?.[field] !== undefined) count += 1;
    });
  });
  return count;
};

export const removeUnresolvedFeatureEdits = ({
  featureOverrides = {}, featurePlacementOverrides = {}, notices = [], scope = ''
} = {}) => (
  removeFeatureEdits({ featureOverrides, featurePlacementOverrides, unresolved: unresolvedNoticeKinds(notices, scope) })
);

export const splitLegacyVisibilityRules = (rules) => {
  const overrides = {};
  const manualRules = [];
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

// The inputs `resolveFeatureDrawn` reads. The feature types come from the
// request's diagram options (a label rerender keeps the last Generate's); the
// edits and rules are the current ones, which a label rerender and Generate
// both carry.
export const featureDrawnContext = (state, {
  diagramOptions = null,
  featureOverrides = state?.featureOverrides
} = {}) => ({
  featureOverrides: featureOverrides || {},
  rules: requestFeatureVisibilityRules(state?.featureVisibilityManualRules),
  selectedTypes: Array.isArray(diagramOptions?.selectedFeaturesSet)
    ? new Set(diagramOptions.selectedFeaturesSet.map(String))
    : null,
  colorRules: Array.isArray(state?.manualSpecificRules) ? state.manualSpecificRules : []
});

// "Is this feature drawn?" as Generate answers it
// (gbdraw/features/visibility.py::should_render_feature): the feature's
// identity-keyed override decides when On or Off, and Exclude from matching
// skips the rules; else the first matching visibility rule decides when Show
// or Off; else the selected feature types, where a feature of another type is
// drawn when a specific color rule matches it. Python matches the rules
// (app/rule-matching.js, R4); the shared vectors in
// tests/fixtures/feature_drawn_cases.json hold both sides to one answer.
// Returns true, false, or null while a match it needs is unknown.
export const resolveFeatureDrawn = (feature, { featureOverrides, rules, selectedTypes, colorRules }) => {
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
  const colorMatches = colorRules.map((rule) => ruleMatchesFeature(feature, rule));
  if (colorMatches.includes(true)) return true;
  return colorMatches.includes(null) ? null : false;
};

// Whether a Result draws a catalog feature (R-5): the resolver's answer, else
// (unknown) whether the Result's catalog (`resultCatalogFeatures`) renders it,
// which Python decided when it last rendered the Result; null without one.
export const featureDrawnInResult = (feature, context, catalogFeatures) => (
  resolveFeatureDrawn(feature, context)
  ?? (catalogFeatures ? catalogFeatures.renderedByIdentity.has(stableKeyOf(feature)) : null)
);

// The Features list of one Result (Owner decision 2026-10-05, R-5): from the
// Result's catalog features, each biological feature that is of a selected
// type, has its own Feature visibility, or is drawn, with whether it is drawn.
// A feature the Result draws is listed as its rendered feature.
export const listFeatureRows = (catalogFeatures, context) => {
  const rows = [];
  const drawn = new Map();
  const rendered = new Set();
  const { selectedTypes, featureOverrides } = context;
  catalogFeatures.biological.forEach((feature) => {
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
// that changes a Result's source asks for the automatic rerender (OV-42,
// OV-43). The color of a rule is in it only for a batch: with one Result the
// rows keep their captions and the Legend follows the color live
// (`ruleLegendCaption` gives the caption of a rule whose caption names another
// row), while another batch Result shows the row it was drawn with until
// Python draws it again (OV-44). A Result is read from its catalog features:
// `asRendered` reads the features Python drew at the last render, else
// `featureDrawnInResult` answers. A source is null while a rule match it needs
// is unknown; null equals no source (`sameLegendSources`).
const legendSourceOf = (catalogFeatures, context, { asRendered, withColors }) => {
  const caption = (rule) => String(rule?.cap ?? '').trim();
  const colorOf = (rule) => String(rule?.color ?? '').trim().toLowerCase();
  const color = (rule) => (withColors ? colorOf(rule) : '');
  const usage = (rule) => JSON.stringify([caption(rule), color(rule)]);
  const rules = context.colorRules;
  const typesByRecord = new Map();
  const winnersByType = new Map();
  for (const feature of catalogFeatures.biological) {
    const shown = catalogFeatures.renderedByIdentity.get(stableKeyOf(feature));
    const row = shown || feature;
    if (!(asRendered ? shown : featureDrawnInResult(row, context, catalogFeatures))) continue;
    const type = String(feature.type ?? '');
    const record = String(feature.recordKey ?? '');
    if (!typesByRecord.has(record)) typesByRecord.set(record, []);
    if (!typesByRecord.get(record).includes(type)) typesByRecord.get(record).push(type);
    const winners = winnersByType.get(type) || [];
    winnersByType.set(type, winners);
    // A type without a captioned rule needs no match.
    if (!rules.some((rule) => caption(rule) && (rule.feat === type || rule.feat === '*'))) continue;
    const winner = firstMatchingRuleIfKnown(row, rules);
    if (winner === undefined) return null;
    winners.push(winner);
  }
  const used = new Set();
  const defaultTypes = new Set();
  const plainTypes = new Set();
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
  const counted = (rule) => !(recolored.has(rule.feat) && caption(rule) === rule.feat);
  return JSON.stringify([
    [...typesByRecord],
    rules.filter((rule) => caption(rule) && counted(rule) && used.has(usage(rule)))
      .map((rule) => [String(rule.feat ?? ''), caption(rule), color(rule)]),
    [...defaultTypes].sort(),
    [...plainTypes].sort(),
    [...new Set(rules.filter(counted).map(caption).filter(Boolean))].sort()
  ]);
};

// The Legend source of every committed Result, in Result order.
export const resultLegendSources = (state, context, { asRendered = false } = {}) => {
  const results = Array.isArray(state?.results?.value) ? state.results.value : [];
  return results.map((_, index) => {
    const catalogFeatures = resultCatalogFeatures(state, index);
    return catalogFeatures
      ? legendSourceOf(catalogFeatures, context, { asRendered, withColors: results.length > 1 })
      : '';
  });
};

export const sameLegendSources = (left, right) => left.length === right.length
  && left.every((source, index) => source !== null && source === right[index]);

const isEditorFeatureRule = (rule) => {
  const normalized = normalizeFeatureVisibilityRule(rule);
  return normalized.source === 'editor' && Boolean(normalized.featureId);
};

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

const reorderEditorVisibilityRules = (rules) => {
  if (!Array.isArray(rules)) return;
  const featureRules = [];
  const qualifierRules = [];
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
 * @param {any} rules
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
