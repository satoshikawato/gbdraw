// @ts-check
import { FEATURE_INSTANCE_SEPARATOR } from '../contracts/schema.generated.js';
/** @import { FeatureLike } from './feature-utils.js' */
const FEATURE_ID_KEY_EMPTY = '';

export const SPECIFIC_COLOR_QUALIFIER_PRESETS = Object.freeze([
  'product',
  'gene',
  'locus_tag',
  'protein_id',
  'note',
  'function',
  'phrog',
  'vfdb_short_name',
  'AMR_Gene_Family',
  'color'
]);

/** @param {unknown} value */
export const escapeRegexLiteral = (value) => String(value ?? '').replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
/** @param {unknown} value */
export const exactRegexValue = (value) => `^${escapeRegexLiteral(value)}$`;
/**
 * The text an exact pattern matches (`exactRegexValue` of it), or null for any
 * other pattern.
 * @param {string} pattern
 * @returns {string | null}
 */
export const exactRegexLiteral = (pattern) => {
  const literal = /^\^([\s\S]*)\$$/.exec(pattern)?.[1].replace(/\\([\s\S])/g, '$1') ?? null;
  return literal !== null && exactRegexValue(literal) === pattern ? literal : null;
};
/**
 * The one literal text of a pattern: an exact pattern's (`exactRegexLiteral`),
 * or a text with no regex syntax as written; null for any other pattern.
 * @param {string} pattern
 * @returns {string | null}
 */
export const literalPatternText = (pattern) => exactRegexLiteral(pattern) ?? (escapeRegexLiteral(pattern) === pattern ? pattern : null);

/**
 * A feature instance a `hash` rule names: `<biological feature ID>` and
 * `<record key>` joined by `FEATURE_INSTANCE_SEPARATOR` (Python's
 * `feature_instance`, gbdraw/features/ids.py).
 * @typedef {{ featureId: string, recordKey: string }} NamedFeatureInstance
 */
/**
 * The feature instance a `hash` rule names: its value is one literal text (as
 * written or as an exact pattern) with the separator, and Python's matcher
 * matches it to the whole instance only (gbdraw/features/selector_values.py).
 * Null for any other rule, a `hash` rule naming a stable hash included.
 * @param {{ qual?: unknown, val?: unknown } | null | undefined} rule
 * @returns {NamedFeatureInstance | null}
 */
export const namedFeatureInstance = (rule) => {
  if (String(rule?.qual ?? '').trim().toLowerCase() !== 'hash') return null;
  const text = literalPatternText(String(rule?.val ?? '').trim());
  const at = text === null ? -1 : text.indexOf(FEATURE_INSTANCE_SEPARATOR);
  return text !== null && at > 0 && at < text.length - 1
    ? { featureId: text.slice(0, at), recordKey: text.slice(at + FEATURE_INSTANCE_SEPARATOR.length) }
    : null;
};

/** @param {unknown} value */
export const normalizeSelectorText = (value) => String(value ?? '').trim();
/** @param {unknown} value */
export const normalizeFeatureIdKey = (value) => normalizeSelectorText(value).toLowerCase();

/**
 * The qualifier names of the features, own and selector qualifiers.
 * @param {(FeatureLike | null | undefined)[] | null} [features]
 * @returns {Set<string>}
 */
export const collectFeatureQualifierNames = (features = []) => {
  /** @type {Set<string>} */
  const names = new Set();
  /** @param {unknown} qualifiers */
  const addQualifierMap = (qualifiers) => {
    if (!qualifiers || typeof qualifiers !== 'object' || Array.isArray(qualifiers)) return;
    Object.keys(qualifiers).forEach((key) => {
      const qualifier = normalizeSelectorText(key);
      if (qualifier) names.add(qualifier);
    });
  };
  if (Array.isArray(features)) {
    features.forEach((feature) => {
      addQualifierMap(feature?.qualifiers);
      addQualifierMap(feature?.selector?.qualifiers);
    });
  }
  return names;
};

/**
 * The presets, then the other qualifier names of the features and the rules, sorted.
 * @param {Iterable<string>} featureQualifierNames from `collectFeatureQualifierNames`
 * @param {readonly ({ qual?: unknown } | null | undefined)[] | null} [rules]
 * @returns {string[]}
 */
export const specificColorQualifierSuggestions = (featureQualifierNames, rules = []) => {
  /** @type {Set<string>} */
  const suggestions = new Set([...SPECIFIC_COLOR_QUALIFIER_PRESETS, ...featureQualifierNames]);
  if (Array.isArray(rules)) {
    rules.forEach((rule) => {
      const qualifier = normalizeSelectorText(rule?.qual);
      if (qualifier) suggestions.add(qualifier);
    });
  }

  const presets = new Set(SPECIFIC_COLOR_QUALIFIER_PRESETS);
  const custom = [...suggestions]
    .filter((qualifier) => !presets.has(qualifier))
    .sort((left, right) => left.localeCompare(right));
  return [...SPECIFIC_COLOR_QUALIFIER_PRESETS, ...custom];
};

/**
 * @param {...unknown} values
 * @returns {string}
 */
const firstText = (...values) => {
  for (const value of values) {
    if (Array.isArray(value)) {
      const found = firstText(...value);
      if (found) return found;
      continue;
    }
    const normalized = normalizeSelectorText(value);
    if (normalized) return normalized;
  }
  return '';
};

/** @param {unknown} value */
const normalizeStrandToken = (value) => {
  const token = normalizeSelectorText(value).toLowerCase();
  if (token === '+' || token === 'positive' || token === 'plus' || token === 'forward' || token === '1') return '+';
  if (token === '-' || token === 'negative' || token === 'minus' || token === 'reverse' || token === '-1') return '-';
  return 'undefined';
};

/**
 * @param {unknown} qualifiers
 * @returns {Record<string, string[]>}
 */
export const normalizeQualifierMap = (qualifiers) => {
  /** @type {Record<string, string[]>} */
  const normalized = {};
  if (!qualifiers || typeof qualifiers !== 'object') return normalized;
  Object.entries(qualifiers).forEach(([keyRaw, valuesRaw]) => {
    const key = normalizeSelectorText(keyRaw).toLowerCase();
    if (!key) return;
    const values = Array.isArray(valuesRaw) ? valuesRaw : [valuesRaw];
    normalized[key] = values
      .filter((value) => value !== null && value !== undefined)
      .map((value) => String(value));
  });
  return normalized;
};

/** @param {Record<string, unknown> | null | undefined} obj */
const hasOwnValues = (obj) => Object.keys(obj || {}).length > 0;

/**
 * @param {unknown} recordId
 * @param {unknown} recordLocation
 */
const recordLocationPosition = (recordId, recordLocation) => {
  const record = normalizeSelectorText(recordId);
  const value = normalizeSelectorText(recordLocation);
  if (!record || !value) return '';
  const prefix = `${record}:`;
  return value.startsWith(prefix) ? value.slice(prefix.length) : '';
};

/**
 * @param {FeatureLike | null | undefined} feature
 * @param {{ requireSelector?: boolean, preferSelector?: boolean }} [options]
 */
export const normalizeFeatureSelectorMetadata = (feature, options = {}) => {
  const requireSelector = options.requireSelector === true;
  const preferSelector = options.preferSelector !== false;
  const selector = feature?.selector && typeof feature.selector === 'object' ? feature.selector : null;
  const hasSelectorMetadata = Boolean(selector);
  const selectorQualifiers = normalizeQualifierMap(selector?.qualifiers);
  const fallbackQualifiers = requireSelector ? {} : normalizeQualifierMap(feature?.qualifiers);
  const qualifiers = preferSelector || requireSelector || hasOwnValues(selectorQualifiers)
    ? selectorQualifiers
    : fallbackQualifiers;

  const recordId = firstText(
    feature?.record_id,
    feature?.recordId,
    feature?.record,
    feature?.record_id_text
  );
  const featureType = firstText(
    feature?.type,
    feature?.featureType,
    feature?.feature_type
  );
  const featureId = firstText(
    feature?.svg_id,
    feature?.svgId,
    selector?.hash,
    feature?.featureId,
    feature?.feature_id,
    feature?.id
  );
  const stableFeatureId = firstText(
    selector?.hash,
    feature?.stable_svg_id,
    feature?.stableSvgId,
    feature?.stableFeatureSvgId,
    feature?.stable_feature_id,
    feature?.stableFeatureId,
    feature?.feature_hash,
    feature?.hash,
    featureId
  );
  const location = firstText(
    selector?.location,
    feature?.location,
    feature?.start !== undefined && feature?.end !== undefined
      ? `${feature.start}..${feature.end}`
      : ''
  );
  const recordLocation = firstText(selector?.record_location, feature?.recordLocation);
  const position = firstText(
    feature?.position,
    recordLocationPosition(recordId, recordLocation),
    location ? `${location}:${normalizeStrandToken(feature?.strand)}` : ''
  );

  return {
    featureId,
    featureIdKey: normalizeFeatureIdKey(featureId) || FEATURE_ID_KEY_EMPTY,
    stableFeatureId,
    stableFeatureIdKey: normalizeFeatureIdKey(stableFeatureId) || FEATURE_ID_KEY_EMPTY,
    recordId,
    record: recordId,
    featureType,
    location,
    position,
    recordLocation,
    qualifiers,
    hasSelectorMetadata
  };
};
