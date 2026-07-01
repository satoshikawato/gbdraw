const REQUIRED_COLUMNS = ['record_id', 'feature_type', 'qualifier', 'value', 'action'];
const COMMON_QUALIFIERS = ['product', 'gene', 'locus_tag', 'hash', 'location', 'record_location'];
const SHOW_ACTIONS = new Set(['show', 'on', 'display', 'include', 'true', '1']);
const HIDE_ACTIONS = new Set(['hide', 'off', 'false', '0']);
const SUPPRESS_ACTIONS = new Set(['suppress', 'exclude']);

let generatedRuleId = 0;

const normalizeCell = (value) => String(value ?? '').replace(/[\t\r\n]+/g, ' ').trim();
const normalizeSource = (value) => {
  const normalized = normalizeCell(value).toLowerCase();
  return ['manual', 'editor', 'file'].includes(normalized) ? normalized : 'manual';
};

export const featureVisibilityQualifierSuggestions = COMMON_QUALIFIERS;

export const escapeRegexLiteral = (value) => String(value ?? '').replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

export const normalizeFeatureVisibilityAction = (value) => {
  const normalized = normalizeCell(value).toLowerCase();
  if (SHOW_ACTIONS.has(normalized)) return 'show';
  if (HIDE_ACTIONS.has(normalized)) return 'hide';
  if (SUPPRESS_ACTIONS.has(normalized)) return 'suppress';
  return '';
};

export const featureVisibilityActionToMode = (value) => {
  const action = normalizeFeatureVisibilityAction(value);
  if (action === 'show') return 'on';
  if (action === 'hide') return 'off';
  if (action === 'suppress') return 'suppress';
  return 'default';
};

export const featureVisibilityModeToAction = (value) => {
  const normalized = normalizeCell(value).toLowerCase();
  if (normalized === 'on') return 'show';
  if (normalized === 'off') return 'hide';
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
  action: 'hide'
});

export const normalizeFeatureVisibilityRule = (raw = {}) => {
  const source = normalizeSource(raw.source);
  const action = normalizeFeatureVisibilityAction(raw.action) || 'hide';
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

const isExactEditorHashRule = (rule, featureId) => {
  const normalized = normalizeFeatureVisibilityRule(rule);
  const expectedFeatureId = normalizeCell(featureId || normalized.featureId);
  if (!expectedFeatureId) return false;
  return normalized.source === 'editor' &&
    normalized.featureId === expectedFeatureId &&
    normalized.recordId === '*' &&
    normalized.featureType === '*' &&
    normalized.qualifier.toLowerCase() === 'hash' &&
    normalized.value === `^${escapeRegexLiteral(expectedFeatureId)}$`;
};

export const upsertEditorFeatureVisibilityRule = (rules, feat, actionRaw) => {
  if (!Array.isArray(rules)) return null;
  const featureId = normalizeCell(feat?.svg_id ?? feat?.svgId ?? feat?.id);
  if (!featureId) return null;
  const action = featureVisibilityModeToAction(actionRaw);
  if (!action) {
    removeEditorFeatureVisibilityRule(rules, featureId);
    return null;
  }

  const nextRule = normalizeFeatureVisibilityRule({
    source: 'editor',
    featureId,
    label: featureLabel(feat),
    recordId: '*',
    featureType: '*',
    qualifier: 'hash',
    value: `^${escapeRegexLiteral(featureId)}$`,
    action
  });
  const existingIndex = rules.findIndex((rule) => normalizeFeatureVisibilityRule(rule).source === 'editor' &&
    normalizeFeatureVisibilityRule(rule).featureId === featureId);

  if (existingIndex >= 0) {
    nextRule.id = normalizeFeatureVisibilityRule(rules[existingIndex]).id;
    rules.splice(existingIndex, 1, nextRule);
  } else {
    rules.unshift(nextRule);
  }
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

export const getEditorFeatureVisibilityMode = (rules, featureIdRaw) => {
  const featureId = normalizeCell(featureIdRaw);
  if (!featureId) return 'default';
  for (const rule of Array.isArray(rules) ? rules : []) {
    if (!isExactEditorHashRule(rule, featureId)) continue;
    return featureVisibilityActionToMode(rule.action);
  }
  return 'default';
};

export const buildFeatureVisibilityOverrideCache = (rules) => {
  const cache = {};
  for (const rule of Array.isArray(rules) ? rules : []) {
    const normalized = normalizeFeatureVisibilityRule(rule);
    if (!isExactEditorHashRule(normalized, normalized.featureId)) continue;
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
      value: `^${escapeRegexLiteral(featureId)}$`,
      action
    }));
  });
  return rules;
};
