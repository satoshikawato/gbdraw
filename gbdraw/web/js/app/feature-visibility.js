import { exactRegexValue } from './feature-selector.js';
import {
  featureIdentityKey,
  featureIdentityKeyOf,
  featureOverrideValue,
  recordKeyBelongsToRequest,
  updateFeatureOverride
} from '../services/feature-placement.js';
export { escapeRegexLiteral, exactRegexValue } from './feature-selector.js';

const REQUIRED_COLUMNS = ['record_id', 'feature_type', 'qualifier', 'value', 'action'];
const COMMON_QUALIFIERS = ['product', 'gene', 'protein_id', 'locus_tag', 'hash', 'location', 'record_location'];
const SHOW_ACTIONS = new Set(['show', 'on']);
const HIDE_ACTIONS = new Set(['hide', 'off', 'false', '0']);
const EXCLUDE_MATCHING_ACTIONS = new Set(['exclude_matching', 'suppress']);

let generatedRuleId = 0;

const normalizeCell = (value) => String(value ?? '').replace(/[\t\r\n]+/g, ' ').trim();
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

export const serializeFeatureVisibilityRules = (rules) => {
  const rows = (Array.isArray(rules) ? rules : [])
    .map((rule) => normalizeFeatureVisibilityRule(rule))
    .filter(isSerializableRule)
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

const getFeatureType = (feat = {}) => normalizeCell(feat?.type ?? feat?.featureType ?? feat?.feature_type) || '*';

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

// One change per feature identity: {recordKey, biologicalFeatureId, featureId
// (the rendered ID it is drawn with), before, after}.
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
    const [recordKey, biologicalFeatureId] = JSON.parse(key);
    changes.push({ recordKey, biologicalFeatureId, featureId: getFeatureId(feature), before, after: mode });
  });
  return changes;
};

export const applyFeatureVisibilityOverrideChanges = (featureOverrides, changes) => {
  if (!featureOverrides || !Array.isArray(changes)) return 0;
  let applied = 0;
  changes.forEach((change) => {
    if (!featureIdentityKey(change?.recordKey, change?.biologicalFeatureId)) return;
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

// Python's `unresolved` notices by identity key: the edit kinds whose feature
// the source does not have (design Q4 3.4).
const unresolvedNoticeKinds = (notices, recordKeyFilter = () => true) => {
  const unresolved = new Map();
  (Array.isArray(notices) ? notices : []).forEach((notice) => {
    if (notice?.status !== 'unresolved' || !recordKeyFilter(notice.recordKey)) return;
    const key = featureIdentityKey(notice.recordKey, notice.biologicalFeatureId);
    if (key) unresolved.set(key, new Set([...(unresolved.get(key) || []), ...(notice.kinds || [])]));
  });
  return unresolved;
};

// Removes edits from the identity-keyed drafts: every edit of a dropped
// record, the unresolved edit kinds of an identity, and (`sourceGone`) the
// label source text kept for a bulk edit. Returns the count of removed edits.
const removeFeatureEdits = ({
  featureOverrides = {},
  featurePlacementOverrides = {},
  unresolved = new Map(),
  dropped = () => false,
  sourceGone = () => false
}) => {
  let removed = 0;
  Object.entries(featurePlacementOverrides || {}).forEach(([key, row]) => {
    if (!dropped(row?.recordKey) && !unresolved.get(key)?.has('placement')) return;
    delete featurePlacementOverrides[key];
    removed += 1;
  });
  Object.entries(featureOverrides || {}).forEach(([key, row]) => {
    const edits = EDIT_FIELDS.filter((field) => row?.[field] !== null && row?.[field] !== undefined);
    if (dropped(row?.recordKey)) {
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
// features outside the crop or display stay dormant. Returns the count of
// removed edits.
export const pruneUnmatchedFeatureOverrides = ({
  featureOverrides = {},
  featurePlacementOverrides = {},
  notices = [],
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
    unresolved: unresolvedNoticeKinds(notices, (recordKey) => replaced.has(recordKey)),
    dropped: (recordKey) => recordKeyBelongsToRequest(recordKey, previousRecords)
      && !recordKeyBelongsToRequest(recordKey, currentRecords),
    sourceGone: (key, row) => replaced.has(row?.recordKey) && !present.has(key)
  });
};

// The unresolved edits the drafts still hold after a Generate that replaced no
// source; "Remove N unmatched feature edits" removes them (an explicit Reset, R2).
export const countUnresolvedFeatureEdits = ({ featureOverrides = {}, featurePlacementOverrides = {}, notices = [] } = {}) => {
  let count = 0;
  unresolvedNoticeKinds(notices).forEach((kinds, key) => {
    if (kinds.has('placement') && featurePlacementOverrides?.[key]) count += 1;
    const row = featureOverrides?.[key];
    Object.entries(KIND_FIELDS).forEach(([kind, field]) => {
      if (kinds.has(kind) && row?.[field] !== null && row?.[field] !== undefined) count += 1;
    });
  });
  return count;
};

export const removeUnresolvedFeatureEdits = ({ featureOverrides = {}, featurePlacementOverrides = {}, notices = [] } = {}) => (
  removeFeatureEdits({ featureOverrides, featurePlacementOverrides, unresolved: unresolvedNoticeKinds(notices) })
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

// Feature visibility for the live preview: the feature's identity-keyed
// override, then the first matching manual rule. With the feature, the
// editor's exact-qualifier rules match as Generate does, so an action and a
// later reconcile agree.
export const resolveEffectiveFeatureVisibility = (feature, featureOverrides = {}, manualRules = []) => {
  const override = getFeatureVisibilityOverride(featureOverrides, feature);
  if (override !== 'default') return override;
  const featureId = getFeatureId(feature || {});
  if (!featureId) return 'default';
  for (const rule of Array.isArray(manualRules) ? manualRules : []) {
    const normalized = normalizeFeatureVisibilityRule(rule);
    if (isEditorExactQualifierRule(normalized)) {
      if (featureMatchesExactQualifier(feature, normalized)) {
        return featureVisibilityActionToMode(normalized.action);
      }
      continue;
    }
    if (normalized.qualifier.toLowerCase() !== 'hash') continue;
    if (normalized.recordId !== '*' || normalized.featureType !== '*') continue;
    try {
      if (new RegExp(normalized.value).test(featureId)) {
        return featureVisibilityActionToMode(normalized.action);
      }
    } catch (_err) {
      return 'default';
    }
  }
  return 'on';
};

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

const qualifierValues = (feat, qualifier) => {
  const key = normalizeCell(qualifier).toLowerCase();
  const qualifiers = feat?.qualifiers && typeof feat.qualifiers === 'object' ? feat.qualifiers : {};
  const entry = Object.entries(qualifiers).find(([name]) => String(name).toLowerCase() === key);
  const raw = entry ? entry[1] : feat?.[key];
  return (Array.isArray(raw) ? raw : [raw]).map((value) => normalizeCell(value)).filter(Boolean);
};

// Python matches qualifier rules with a case-insensitive regex search over
// every value (gbdraw/features/visibility.py). The editor's exact-qualifier
// scope and its rule use this one matcher.
export const featureMatchesExactQualifier = (feat, { featureType, qualifier, value } = {}) => {
  if (getFeatureType(feat) !== normalizeCell(featureType)) return false;
  let pattern;
  try {
    pattern = new RegExp(normalizeCell(value), 'i');
  } catch (_err) {
    return false;
  }
  return qualifierValues(feat, qualifier).some((candidate) => pattern.test(candidate));
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
