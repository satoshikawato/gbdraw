import {
  buildSelectorSafetyUniquenessIndex,
  exactRegexValue,
  selectFeatureSelector
} from './feature-selector.js';
export { escapeRegexLiteral, exactRegexValue } from './feature-selector.js';

const REQUIRED_COLUMNS = ['record_id', 'feature_type', 'qualifier', 'value', 'action'];
const COMMON_QUALIFIERS = ['product', 'gene', 'protein_id', 'locus_tag', 'hash', 'location', 'record_location'];
const EDITOR_FEATURE_SELECTOR_PRIORITY = ['protein_id', 'locus_tag', 'gene_id', 'old_locus_tag'];
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

const firstText = (...values) => {
  for (const value of values) {
    if (Array.isArray(value)) {
      const found = firstText(...value);
      if (found) return found;
      continue;
    }
    const normalized = normalizeCell(value);
    if (normalized) return normalized;
  }
  return '';
};

const featureLabel = (feat = {}) => {
  const qualifiers = feat.qualifiers && typeof feat.qualifiers === 'object' ? feat.qualifiers : {};
  return firstText(
    feat.label,
    feat.product,
    feat.gene,
    feat.locus_tag,
    qualifiers.product,
    qualifiers.gene,
    qualifiers.locus_tag,
    feat.type && (feat.start || feat.end) ? `${feat.type} ${feat.start ?? '?'}..${feat.end ?? '?'}` : feat.type
  );
};

const getFeatureId = (feat = {}) => normalizeCell(feat?.svg_id ?? feat?.svgId ?? feat?.featureId ?? feat?.feature_id ?? feat?.id);

const getFeatureStableHashValue = (feat = {}) => normalizeCell(
  feat?.selector?.hash ??
  feat?.stable_svg_id ??
  feat?.stableSvgId ??
  feat?.stableFeatureSvgId ??
  feat?.stable_feature_id ??
  feat?.stableFeatureId ??
  feat?.feature_hash ??
  feat?.hash ??
  getFeatureId(feat)
);

const getFeatureRecordId = (feat = {}) => normalizeCell(feat?.record_id ?? feat?.recordId ?? feat?.record) || '*';

const getFeatureType = (feat = {}) => normalizeCell(feat?.type ?? feat?.featureType ?? feat?.feature_type) || '*';

const getSelectorCacheEntry = (selectorCache, featureIdRaw) => {
  const featureId = normalizeCell(featureIdRaw);
  if (!featureId) return null;
  if (selectorCache instanceof Map) return selectorCache.get(featureId) || null;
  if (!selectorCache || typeof selectorCache !== 'object' || Array.isArray(selectorCache)) return null;
  return selectorCache[featureId] || null;
};

const hasOverrideKey = (overrides, featureId) =>
  Boolean(overrides && typeof overrides === 'object' && Object.prototype.hasOwnProperty.call(overrides, featureId));

export const getFeatureVisibilityOverride = (overrides, featureIdRaw) => {
  const featureId = normalizeCell(featureIdRaw);
  if (!featureId || !hasOverrideKey(overrides, featureId)) return 'default';
  return normalizeVisibilityMode(overrides[featureId]);
};

export const setFeatureVisibilityOverride = (overrides, featureIdRaw, modeRaw) => {
  if (!overrides || typeof overrides !== 'object') return 'default';
  const featureId = normalizeCell(featureIdRaw);
  if (!featureId) return 'default';
  const previous = getFeatureVisibilityOverride(overrides, featureId);
  const mode = normalizeVisibilityMode(modeRaw);
  if (mode === 'default') {
    delete overrides[featureId];
  } else {
    overrides[featureId] = mode;
  }
  return previous;
};

export const applyFeatureVisibilityOverrideChanges = (overrides, changes) => {
  if (!overrides || typeof overrides !== 'object' || !Array.isArray(changes)) return 0;
  let applied = 0;
  changes.forEach((change) => {
    const featureId = normalizeCell(change?.featureId ?? change?.svgId ?? change?.id);
    if (!featureId) return;
    const mode = Object.prototype.hasOwnProperty.call(change || {}, 'mode')
      ? change.mode
      : change?.after;
    setFeatureVisibilityOverride(overrides, featureId, mode);
    applied += 1;
  });
  return applied;
};

export const buildFeatureVisibilityChanges = (features, modeRaw, currentOverrides = {}) => {
  const mode = normalizeVisibilityMode(modeRaw);
  const seen = new Set();
  const changes = [];
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const featureId = getFeatureId(feature);
    if (!featureId || seen.has(featureId)) return;
    seen.add(featureId);
    const before = getFeatureVisibilityOverride(currentOverrides, featureId);
    if (before === mode) return;
    changes.push({ featureId, before, after: mode });
  });
  return changes;
};

export const buildFeatureVisibilitySelectorCache = (features, selectorSafetyScope) => {
  const cache = {};
  const selectorUniquenessIndex = buildSelectorSafetyUniquenessIndex(selectorSafetyScope);
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const featureId = getFeatureId(feature);
    if (!featureId || cache[featureId]) return;
    const selector = selectFeatureSelector(
      feature,
      selectorUniquenessIndex,
      {
        priority: EDITOR_FEATURE_SELECTOR_PRIORITY,
        requireSelector: true,
        requireSafetyScope: true
      }
    );
    const value = normalizeCell(selector?.value || featureId);
    if (!value) return;
    cache[featureId] = {
      recordId: getFeatureRecordId(feature),
      featureType: getFeatureType(feature),
      qualifier: normalizeCell(selector?.qualifier || 'hash').toLowerCase() || 'hash',
      value,
      label: featureLabel(feature)
    };
  });
  return cache;
};

export const preserveFeatureVisibilitySelectorCacheForOverrides = (
  nextCache = {},
  previousCache = {},
  overrides = {}
) => {
  const merged = { ...(nextCache && typeof nextCache === 'object' && !Array.isArray(nextCache) ? nextCache : {}) };
  if (!overrides || typeof overrides !== 'object' || Array.isArray(overrides)) return merged;

  Object.entries(overrides).forEach(([featureIdRaw, modeRaw]) => {
    if (normalizeVisibilityMode(modeRaw) === 'default') return;
    const featureId = normalizeCell(featureIdRaw);
    if (!featureId || merged[featureId]) return;
    const previous = getSelectorCacheEntry(previousCache, featureId);
    if (previous) merged[featureId] = previous;
  });

  return merged;
};

const buildRuleFromSelectorCacheEntry = (featureIdRaw, modeRaw, selectorCache) => {
  const featureId = normalizeCell(featureIdRaw);
  const action = featureVisibilityModeToAction(modeRaw);
  if (!featureId || !action) return null;
  const cached = getSelectorCacheEntry(selectorCache, featureId) || {};
  const qualifier = normalizeCell(cached.qualifier).toLowerCase() || 'hash';
  const value = normalizeCell(cached.value) || featureId;
  return normalizeFeatureVisibilityRule({
    source: 'editor',
    featureId,
    label: normalizeCell(cached.label) || featureId,
    recordId: normalizeCell(cached.recordId) || '*',
    featureType: normalizeCell(cached.featureType) || '*',
    qualifier,
    value: exactRegexValue(value),
    action
  });
};

export const featureVisibilityOverridesToRules = (overrides, selectorCache = {}) => {
  if (!overrides || typeof overrides !== 'object' || Array.isArray(overrides)) return [];
  return Object.entries(overrides)
    .map(([featureId, mode]) => buildRuleFromSelectorCacheEntry(featureId, mode, selectorCache))
    .filter(Boolean);
};

export const deriveFeatureVisibilityRulesForBoundary = (
  manualRules = [],
  overrides = {},
  selectorCache = {}
) => [
  ...featureVisibilityOverridesToRules(overrides, selectorCache),
  ...(Array.isArray(manualRules) ? manualRules.map((rule) => normalizeFeatureVisibilityRule(rule)) : [])
];

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

const getCacheValue = (cache, featureId) => {
  if (cache instanceof Map) return cache.get(featureId);
  if (!cache || typeof cache !== 'object' || Array.isArray(cache)) return undefined;
  return cache[featureId];
};

export const resolveEffectiveFeatureVisibility = (
  featureIdRaw,
  overrides = {},
  baseVisibilityCache = null,
  manualRules = []
) => {
  const featureId = normalizeCell(featureIdRaw);
  if (!featureId) return 'default';
  const override = getFeatureVisibilityOverride(overrides, featureId);
  if (override !== 'default') return override;
  const cached = normalizeVisibilityMode(getCacheValue(baseVisibilityCache, featureId));
  if (cached !== 'default') return cached;
  for (const rule of Array.isArray(manualRules) ? manualRules : []) {
    const normalized = normalizeFeatureVisibilityRule(rule);
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

const isEditorRuleForFeatureId = (rule, featureId) => {
  const normalized = normalizeFeatureVisibilityRule(rule);
  const expectedFeatureId = normalizeCell(featureId || normalized.featureId);
  if (!expectedFeatureId) return false;
  return normalized.source === 'editor' && normalized.featureId === expectedFeatureId;
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

export const buildExactHashFeatureVisibilityRule = (feat, actionRaw) => {
  const featureId = getFeatureId(feat);
  const hashValue = getFeatureStableHashValue(feat);
  const action = featureVisibilityModeToAction(actionRaw);
  if (!featureId || !hashValue || !action) return null;
  return normalizeFeatureVisibilityRule({
    source: 'editor',
    featureId,
    label: featureLabel(feat),
    recordId: '*',
    featureType: '*',
    qualifier: 'hash',
    value: exactRegexValue(hashValue),
    action
  });
};

const resolveSelectorUniquenessIndex = (selectorContext = {}) => {
  if (selectorContext?.selectorUniquenessIndex instanceof Map) {
    return selectorContext.selectorUniquenessIndex;
  }
  return buildSelectorSafetyUniquenessIndex(selectorContext?.selectorSafetyScope);
};

export const buildEditorFeatureVisibilityRule = (feat, selectorContext = {}, mode) => {
  const featureId = getFeatureId(feat);
  const action = featureVisibilityModeToAction(mode);
  if (!featureId || !action) return null;

  const selector = selectFeatureSelector(
    feat,
    resolveSelectorUniquenessIndex(selectorContext),
    {
      priority: EDITOR_FEATURE_SELECTOR_PRIORITY,
      requireSelector: true,
      requireSafetyScope: true
    }
  );
  const selectedValue = normalizeCell(selector?.value || featureId);
  if (!selectedValue) return null;

  return normalizeFeatureVisibilityRule({
    source: 'editor',
    featureId,
    label: featureLabel(feat),
    recordId: getFeatureRecordId(feat),
    featureType: getFeatureType(feat),
    qualifier: normalizeCell(selector?.qualifier || 'hash').toLowerCase() || 'hash',
    value: exactRegexValue(selectedValue),
    action
  });
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

export const upsertEditorFeatureVisibilityRule = (rules, feat, actionRaw, selectorContext = {}) => {
  if (!Array.isArray(rules)) return null;
  const featureId = getFeatureId(feat);
  if (!featureId) return null;
  const action = featureVisibilityModeToAction(actionRaw);
  if (!action) {
    removeEditorFeatureVisibilityRule(rules, featureId);
    return null;
  }

  const nextRule = buildEditorFeatureVisibilityRule(feat, selectorContext, action);
  if (!nextRule) return null;
  const existingIndex = rules.findIndex((rule) => normalizeFeatureVisibilityRule(rule).source === 'editor' &&
    normalizeFeatureVisibilityRule(rule).featureId === featureId);

  if (existingIndex >= 0) {
    nextRule.id = normalizeFeatureVisibilityRule(rules[existingIndex]).id;
    rules.splice(existingIndex, 1, nextRule);
  } else {
    rules.unshift(nextRule);
  }
  reorderEditorVisibilityRules(rules);
  return nextRule;
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

export const removeEditorFeatureVisibilityRule = (rules, featureIdRaw) => {
  if (!Array.isArray(rules)) return 0;
  const featureId = normalizeCell(featureIdRaw);
  if (!featureId) return 0;
  let removed = 0;
  for (let index = rules.length - 1; index >= 0; index -= 1) {
    const rule = normalizeFeatureVisibilityRule(rules[index]);
    if (rule.source !== 'editor' || rule.featureId !== featureId) continue;
    rules.splice(index, 1);
    removed += 1;
  }
  return removed;
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

export const getEditorFeatureVisibilityMode = (rules, featureIdRaw) => {
  const featureId = normalizeCell(featureIdRaw);
  if (!featureId) return 'default';
  for (const rule of Array.isArray(rules) ? rules : []) {
    if (!isEditorRuleForFeatureId(rule, featureId)) continue;
    return featureVisibilityActionToMode(rule.action);
  }
  return 'default';
};

export const buildFeatureVisibilityOverrideCache = (rules) => {
  const cache = {};
  for (const rule of Array.isArray(rules) ? rules : []) {
    const normalized = normalizeFeatureVisibilityRule(rule);
    if (!isEditorRuleForFeatureId(normalized, normalized.featureId)) continue;
    const mode = featureVisibilityActionToMode(normalized.action);
    if (mode === 'default') continue;
    cache[normalized.featureId] = mode;
  }
  return cache;
};

export const featureVisibilityRulesFromOverrideCache = (overrides) => {
  if (!overrides || typeof overrides !== 'object' || Array.isArray(overrides)) return [];
  const rules = [];
  Object.entries(overrides).forEach(([featureIdRaw, modeRaw]) => {
    const featureId = normalizeCell(featureIdRaw);
    const action = featureVisibilityModeToAction(modeRaw);
    if (!featureId || !action) return;
    rules.push(normalizeFeatureVisibilityRule({
      source: 'editor',
      featureId,
      label: featureId,
      recordId: '*',
      featureType: '*',
      qualifier: 'hash',
      value: exactRegexValue(featureId),
      action
    }));
  });
  return rules;
};
