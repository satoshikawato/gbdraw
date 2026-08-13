import { FEATURE_INSTANCE_HASH_QUALIFIER } from '../../services/feature-instance-identity.js';
import {
  resolveOrderedSpecificRule,
  specificRuleMatchesFeature
} from '../feature-utils.js';

export {
  FEATURE_INSTANCE_HASH_QUALIFIER,
  resolveOrderedSpecificRule,
  specificRuleMatchesFeature
};

const RULE_COLOR_PATTERN = /^#(?:[0-9a-f]{3}|[0-9a-f]{4}|[0-9a-f]{6}|[0-9a-f]{8})$/i;
const SIX_DIGIT_HEX_PATTERN = /^#[0-9a-f]{6}$/i;

const text = (value) => String(value ?? '').trim();

export const isRuleColor = (value) => {
  const normalized = text(value);
  return normalized.toLowerCase() === 'none' || RULE_COLOR_PATTERN.test(normalized);
};

export const normalizeRuleColor = (value) => {
  const normalized = text(value).toLowerCase();
  return isRuleColor(normalized) ? normalized : null;
};

export const normalizePickerColor = (value) => {
  const normalized = text(value).toLowerCase();
  return SIX_DIGIT_HEX_PATTERN.test(normalized) ? normalized : null;
};

export const normalizeFeatureFillValue = (value) => {
  const kind = value && typeof value === 'object'
    ? text(value.kind).toLowerCase()
    : '';
  if (value === null || kind === 'inherit') return { kind: 'inherit' };
  if (kind === 'none' || text(value).toLowerCase() === 'none') return { kind: 'none' };
  const rawColor = kind === 'color' ? value.color : value;
  const color = normalizeRuleColor(rawColor);
  return color && color !== 'none' ? { kind: 'color', color } : null;
};

const ruleColor = (rule) => normalizeRuleColor(rule?.color);

const isExactLocalRule = (rule) => (
  text(rule?.qual) === FEATURE_INSTANCE_HASH_QUALIFIER
  && text(rule?.match).toLowerCase() !== 'regex'
);

const ruleOriginLabel = (rule) => {
  const identity = text(rule?.cap) || text(rule?.feat) || text(rule?.qual);
  return identity ? `Inherited from rule: ${identity}` : 'Inherited from specific rule';
};

const explicitFillValue = ({ explicitValue, explicitOrigin, localRule, feature }) => {
  if (explicitValue !== undefined && explicitOrigin === 'local') {
    const normalized = normalizeFeatureFillValue(explicitValue);
    return normalized?.kind === 'inherit'
      ? null
      : (normalized?.kind === 'none' ? 'none' : normalized?.color || null);
  }
  const featureValue = feature?.explicitFillValue ?? feature?.explicit_fill_value;
  if (featureValue !== undefined) {
    const normalized = normalizeFeatureFillValue(featureValue);
    return normalized?.kind === 'inherit'
      ? null
      : (normalized?.kind === 'none' ? 'none' : normalized?.color || null);
  }
  return ruleColor(localRule);
};

export const resolveFeatureFillViewModel = ({
  feature = null,
  explicitValue = undefined,
  explicitOrigin = 'local',
  localRule = null,
  effectiveRule = undefined,
  specificRules = [],
  paletteColors = {},
  svgDefaultColor = '#cccccc',
  allowNone = true
} = {}) => {
  const resolved = effectiveRule === undefined
    ? resolveOrderedSpecificRule(feature, specificRules)
    : null;
  const resolvedRule = effectiveRule === undefined ? resolved?.rule || null : effectiveRule;
  const reservedLocalRule = isExactLocalRule(resolvedRule) ? resolvedRule : null;
  const effectiveLocalRule = isExactLocalRule(localRule) ? localRule : reservedLocalRule;
  const explicit = explicitFillValue({
    explicitValue,
    explicitOrigin: text(explicitOrigin).toLowerCase(),
    localRule: effectiveLocalRule,
    feature
  });
  const inheritedRuleColor = explicit === null ? ruleColor(resolvedRule) : null;
  const paletteColor = normalizeRuleColor(paletteColors?.[feature?.type]);
  const defaultColor = normalizeRuleColor(svgDefaultColor) || '#cccccc';

  let effectiveColor;
  let origin;
  let originLabel;
  if (explicit !== null) {
    effectiveColor = explicit;
    origin = 'local';
    originLabel = 'Feature override';
  } else if (inheritedRuleColor) {
    effectiveColor = inheritedRuleColor;
    origin = 'specific-rule';
    originLabel = ruleOriginLabel(resolvedRule);
  } else if (paletteColor) {
    effectiveColor = paletteColor;
    origin = 'palette';
    originLabel = 'Inherited from palette';
  } else {
    effectiveColor = defaultColor;
    origin = 'svg-default';
    originLabel = 'Inherited from SVG default';
  }

  return Object.freeze({
    effectiveColor,
    explicitValue: explicit,
    origin,
    originLabel,
    canReset: explicit !== null,
    allowNone: Boolean(allowNone)
  });
};

export const getFeatureFillViewModel = resolveFeatureFillViewModel;
