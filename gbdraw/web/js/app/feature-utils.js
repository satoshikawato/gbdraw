import {
  FEATURE_INSTANCE_HASH_QUALIFIER as INSTANCE_HASH_QUALIFIER,
  FEATURE_INSTANCE_HASH_PATTERN as INSTANCE_HASH_PATTERN
} from '../services/feature-instance-identity.js';
import {
  featureSemanticSelectorMatches,
  isFeatureSemanticScopeRule,
  parseFeatureSemanticSelector
} from './feature-editor/semantic-fill-selectors.js';

export const normalizeStringArray = (value) => {
  if (Array.isArray(value)) {
    return value
      .filter((item) => item !== null && item !== undefined)
      .map((item) => String(item));
  }
  if (value === null || value === undefined || value === '') return [];
  return [String(value)];
};

export const getFeatureCaption = (feature) => {
  const caption = firstNonInternalProteinDisplayText(
    feature?.label,
    feature?.display_label,
    feature?.displayLabel,
    feature?.product,
    feature?.gene,
    feature?.locus_tag,
    feature?.note
  );
  return caption || `${feature?.type} at ${feature?.start}..${feature?.end}`;
};

const firstFeatureText = (...values) => {
  for (const value of values) {
    if (Array.isArray(value)) {
      const nested = firstFeatureText(...value);
      if (nested) return nested;
      continue;
    }
    const text = String(value === null || value === undefined ? '' : value).trim();
    if (text) return text;
  }
  return '';
};

const LINEAR_RENDERED_RECORD_SUFFIX = /_record_\d+$/i;

export const getFeatureHashCandidates = (feature) => {
  if (!feature || typeof feature !== 'object') return [];

  const renderedId = firstFeatureText(
    feature.rendered_svg_id,
    feature.renderedSvgId,
    feature.rendered_feature_svg_id,
    feature.renderedFeatureSvgId,
    feature.svg_id,
    feature.svgId
  );
  const generationId = renderedId.replace(LINEAR_RENDERED_RECORD_SUFFIX, '');

  return [...new Set([generationId, renderedId].filter(Boolean))];
};

export const getFeatureGenerationHash = (feature) => getFeatureHashCandidates(feature)[0] || '';

const directFeatureValue = (feature, ...keys) => {
  if (!feature || typeof feature !== 'object') return '';
  for (const key of keys) {
    const value = firstFeatureText(feature[key]);
    if (value) return value;
  }
  return '';
};

const RUNTIME_PROTEIN_HANDLE_RE = /^h_[a-z2-7]{26}$/;
const FEATURE_ANALYSIS_ID_RE = /^f_[0-9a-f]{64}$/;
// Keep unsupported historical shapes from leaking through display-only fallbacks.
const UNSUPPORTED_HISTORICAL_PROTEIN_ID_RE =
  /@[^|]+\|.+~f_[0-9a-f]{64}$/;
const GENERATED_PROTEIN_ID_RE =
  /^(?:gbd_r\d+_(?:cds\d+|unit\d+)|p_.+_\d+_\d+_-?\d+_[0-9a-f]{12}(?:_\d+)?)$/i;

export const isInternalProteinDisplayId = (value) => {
  const text = firstFeatureText(value);
  return Boolean(
    text &&
    (
      RUNTIME_PROTEIN_HANDLE_RE.test(text) ||
      FEATURE_ANALYSIS_ID_RE.test(text) ||
      text.startsWith('p_r_') ||
      UNSUPPORTED_HISTORICAL_PROTEIN_ID_RE.test(text) ||
      GENERATED_PROTEIN_ID_RE.test(text)
    )
  );
};

const firstNonInternalProteinDisplayText = (...values) => {
  for (const value of values) {
    if (Array.isArray(value)) {
      const nested = firstNonInternalProteinDisplayText(...value);
      if (nested) return nested;
      continue;
    }
    const text = firstFeatureText(value);
    if (text && !isInternalProteinDisplayId(text)) return text;
  }
  return '';
};

const directFeatureDisplayValue = (feature, ...keys) => {
  if (!feature || typeof feature !== 'object') return '';
  for (const key of keys) {
    const value = firstNonInternalProteinDisplayText(feature[key]);
    if (value) return value;
  }
  return '';
};

const getFeatureQualifierDisplayValue = (feature, key) => {
  const normalizedKey = String(key || '').trim().toLowerCase();
  if (!feature || !normalizedKey) return '';
  const qualifiers = feature.qualifiers && typeof feature.qualifiers === 'object' && !Array.isArray(feature.qualifiers)
    ? feature.qualifiers
    : {};
  if (Object.prototype.hasOwnProperty.call(qualifiers, normalizedKey)) {
    return firstNonInternalProteinDisplayText(qualifiers[normalizedKey]);
  }
  const matchingKey = Object.keys(qualifiers).find((candidate) => candidate.toLowerCase() === normalizedKey);
  return matchingKey
    ? firstNonInternalProteinDisplayText(qualifiers[matchingKey])
    : '';
};

export const resolveDisplayProteinId = (feature, member = null, fallback = '') =>
  firstNonInternalProteinDisplayText(
    directFeatureDisplayValue(feature, 'displayProteinId', 'display_protein_id'),
    directFeatureDisplayValue(member, 'displayProteinId', 'display_protein_id'),
    directFeatureDisplayValue(feature, 'sourceProteinId', 'source_protein_id'),
    directFeatureDisplayValue(member, 'sourceProteinId', 'source_protein_id'),
    getFeatureQualifierDisplayValue(feature, 'protein_id'),
    directFeatureDisplayValue(feature, 'locusTag', 'locus_tag'),
    getFeatureQualifierDisplayValue(feature, 'locus_tag'),
    directFeatureDisplayValue(member, 'locusTag', 'locus_tag'),
    directFeatureDisplayValue(feature, 'geneId', 'gene_id'),
    getFeatureQualifierDisplayValue(feature, 'gene_id'),
    directFeatureDisplayValue(member, 'geneId', 'gene_id'),
    directFeatureDisplayValue(feature, 'oldLocusTag', 'old_locus_tag'),
    getFeatureQualifierDisplayValue(feature, 'old_locus_tag'),
    directFeatureDisplayValue(member, 'oldLocusTag', 'old_locus_tag'),
    directFeatureDisplayValue(feature, 'ID'),
    getFeatureQualifierDisplayValue(feature, 'ID'),
    directFeatureDisplayValue(feature, 'Name', 'name'),
    getFeatureQualifierDisplayValue(feature, 'Name'),
    directFeatureDisplayValue(feature, 'Parent', 'parent'),
    getFeatureQualifierDisplayValue(feature, 'Parent'),
    directFeatureDisplayValue(feature, 'gene'),
    getFeatureQualifierDisplayValue(feature, 'gene'),
    directFeatureDisplayValue(member, 'gene'),
    directFeatureDisplayValue(member, 'label'),
    directFeatureDisplayValue(feature, 'proteinId', 'protein_id'),
    directFeatureDisplayValue(member, 'proteinId', 'protein_id'),
    fallback
  );

export const resolveInternalProteinId = (feature, member = null, fallback = '') => firstFeatureText(
  directFeatureValue(feature, 'proteinId', 'protein_id'),
  directFeatureValue(member, 'proteinId', 'protein_id'),
  fallback
);

const getFeatureQualifierValue = (feat, qual) => {
  if (!qual) return null;
  const key = qual.toLowerCase();
  const qualifiers = feat?.qualifiers && typeof feat.qualifiers === 'object'
    && !Array.isArray(feat.qualifiers)
    ? feat.qualifiers
    : {};
  const sourceKey = Object.keys(qualifiers).find(
    (candidate) => candidate.toLowerCase() === key
  );
  if (sourceKey) {
    return qualifiers[sourceKey];
  }
  if (key === 'product') return feat.product || null;
  if (key === 'gene') return feat.gene || null;
  if (key === 'locus_tag') return feat.locus_tag || null;
  if (key === 'note') return feat.note || null;
  return null;
};

const matchRuleValue = (value, ruleVal, strictEquals = false) => {
  if (value === null || value === undefined || ruleVal === null || ruleVal === undefined) {
    return false;
  }
  const values = Array.isArray(value) ? value : [value];
  const needle = String(ruleVal);
  let regex = null;
  if (!strictEquals) {
    try {
      regex = new RegExp(needle, 'i');
    } catch {
      return false;
    }
  }
  for (const item of values) {
    if (item === null || item === undefined) continue;
    const text = String(item);
    if (strictEquals ? text === needle : regex.test(text)) return true;
  }
  return false;
};

const normalizeFeatureStrand = (value) => {
  if (value === null || value === undefined) return 'undefined';
  const normalized = String(value).trim().toLowerCase();
  if (['positive', 'plus', '+', 'forward', '1'].includes(normalized)) return '+';
  if (['negative', 'minus', '-', 'reverse', '-1'].includes(normalized)) return '-';
  return 'undefined';
};

const getFeatureLocationValue = (feat) => firstFeatureText(
  feat?.selector?.location,
  feat?.location,
  feat?.start !== undefined && feat?.end !== undefined ? `${feat.start}..${feat.end}` : ''
);

const getFeatureRecordLocationValue = (feat) => {
  const explicit = firstFeatureText(
    feat?.selector?.record_location,
    feat?.selector?.recordLocation,
    feat?.record_location,
    feat?.recordLocation
  );
  if (explicit) return explicit;
  const recordId = firstFeatureText(feat?.record_id, feat?.recordId, feat?.record);
  const location = getFeatureLocationValue(feat);
  return recordId && location
    ? `${recordId}:${location}:${normalizeFeatureStrand(feat?.strand)}`
    : '';
};

const getRuleHashCandidates = (feature) => [...new Set([
  feature?.selector?.hash,
  feature?.hash,
  feature?.stableFeatureId,
  feature?.stable_feature_id,
  ...getFeatureHashCandidates(feature)
].map((value) => firstFeatureText(value)).filter(Boolean))];

const isInstanceLiteralRule = (rule) => (
  String(rule?.qual || '').trim() === INSTANCE_HASH_QUALIFIER
  && String(rule?.match || '').trim().toLowerCase() !== 'regex'
);

const ruleKind = (rule) => {
  if (isInstanceLiteralRule(rule)) return 'instance-literal';
  if (isFeatureSemanticScopeRule(rule)) return 'semantic-literal';
  const qualifier = String(rule?.qual || '').trim().toLowerCase();
  if (qualifier === 'hash') return 'hash';
  if (qualifier === 'record_location') return 'record-location';
  if (qualifier === 'location') return 'location';
  return 'source-qualifier';
};

const ruleTypeMatchesFeature = (feat, rule) => {
  const featureType = String(feat?.type || '').trim();
  const ruleType = String(rule?.feat || '').trim();
  return Boolean(featureType && (ruleType === featureType || ruleType === '*'));
};

/** Match one normalized specific-color rule without applying precedence. */
export const specificRuleMatchesFeature = (
  feat,
  rule,
  { baseLegendCaption = undefined } = {}
) => {
  if (!feat || !rule || !ruleTypeMatchesFeature(feat, rule)) return false;
  const rawQual = String(rule.qual || '').trim();
  const qualKey = rawQual.toLowerCase();

  if (isInstanceLiteralRule(rule)) {
    const instanceHash = firstFeatureText(feat.instanceHash, feat.instance_hash);
    const literal = String(rule.val || '').trim();
    return INSTANCE_HASH_PATTERN.test(literal) && instanceHash === literal;
  }

  if (isFeatureSemanticScopeRule(rule)) {
    const resolvedBaseLegendCaption = baseLegendCaption === undefined
      ? firstFeatureText(
        feat.baseLegendCaption,
        feat.base_legend_caption,
        feat.type
      )
      : String(baseLegendCaption ?? '');
    return featureSemanticSelectorMatches(
      feat,
      parseFeatureSemanticSelector(rule.val),
      { baseLegendCaption: resolvedBaseLegendCaption }
    );
  }

  if (qualKey === 'hash') {
    return getRuleHashCandidates(feat).some((candidate) => matchRuleValue(candidate, rule.val));
  }

  if (qualKey === 'record_location') {
    return matchRuleValue(getFeatureRecordLocationValue(feat), rule.val);
  }

  if (qualKey === 'location') {
    return matchRuleValue(getFeatureLocationValue(feat), rule.val);
  }

  const value = getFeatureQualifierValue(feat, qualKey);
  if (!value) return false;
  return matchRuleValue(value, rule.val);
};

const sourceQualifierOrder = (feature) => {
  const qualifiers = feature?.qualifiers && typeof feature.qualifiers === 'object'
    && !Array.isArray(feature.qualifiers)
    ? feature.qualifiers
    : {};
  const ordered = [];
  const seen = new Set();
  Object.keys(qualifiers).forEach((key) => {
    const normalized = String(key).trim().toLowerCase();
    if (!normalized || seen.has(normalized)) return;
    seen.add(normalized);
    ordered.push(normalized);
  });
  for (const key of ['product', 'gene', 'locus_tag', 'note']) {
    if (seen.has(key) || !firstFeatureText(feature?.[key])) continue;
    seen.add(key);
    ordered.push(key);
  }
  return ordered;
};

const firstMatchingRuleOfKind = (feature, rules, kind, matchOptions = undefined) => {
  for (let index = 0; index < rules.length; index += 1) {
    const rule = rules[index];
    if (
      ruleKind(rule) !== kind
      || !specificRuleMatchesFeature(feature, rule, matchOptions)
    ) {
      continue;
    }
    return { rule, ruleIndex: index };
  }
  return null;
};

const firstMatchingLegacyRule = (feature, rules) => {
  for (const kind of ['hash', 'record-location']) {
    const matched = firstMatchingRuleOfKind(feature, rules, kind);
    if (matched) return matched;
  }

  for (const qualifier of sourceQualifierOrder(feature)) {
    for (let index = 0; index < rules.length; index += 1) {
      const rule = rules[index];
      if (
        ruleKind(rule) !== 'source-qualifier'
        || String(rule?.qual || '').trim().toLowerCase() !== qualifier
        || !specificRuleMatchesFeature(feature, rule)
      ) {
        continue;
      }
      return { rule, ruleIndex: index };
    }
  }

  return firstMatchingRuleOfKind(feature, rules, 'location');
};

/**
 * Resolve renderer-visible rule precedence: exact Feature type, then wildcard;
 * inside either set, instance literal, hash, record-location, source qualifier
 * order, then location. Rule-array order is stable within one selector class.
 */
export const resolveOrderedSpecificRule = (feature, rules = []) => {
  const source = Array.isArray(rules) ? rules : [];
  const featureType = String(feature?.type || '').trim();
  if (!featureType) return null;

  const typedRuleSets = [featureType, '*'].map((candidateType) => ({
    candidateType,
    rules: source.filter(
      (rule) => String(rule?.feat || '').trim() === candidateType
    )
  }));
  const underlyingLegacyRule = typedRuleSets
    .map(({ rules: typedRules }) => firstMatchingLegacyRule(feature, typedRules))
    .find(Boolean);
  const baseLegendCaption = firstFeatureText(
    underlyingLegacyRule?.rule?.cap,
    featureType
  );

  for (const { rules: typedRules } of typedRuleSets) {
    for (const kind of [
      'instance-literal',
      'semantic-literal',
      'hash',
      'record-location'
    ]) {
      const matched = firstMatchingRuleOfKind(
        feature,
        typedRules,
        kind,
        kind === 'semantic-literal' ? { baseLegendCaption } : undefined
      );
      if (matched) {
        return {
          rule: matched.rule,
          ruleIndex: source.indexOf(matched.rule)
        };
      }
    }

    for (const qualifier of sourceQualifierOrder(feature)) {
      for (const rule of typedRules) {
        if (
          ruleKind(rule) !== 'source-qualifier'
          || String(rule?.qual || '').trim().toLowerCase() !== qualifier
          || !specificRuleMatchesFeature(feature, rule)
        ) {
          continue;
        }
        return { rule, ruleIndex: source.indexOf(rule) };
      }
    }

    const locationRule = firstMatchingRuleOfKind(feature, typedRules, 'location');
    if (locationRule) {
      return {
        rule: locationRule.rule,
        ruleIndex: source.indexOf(locationRule.rule)
      };
    }
  }
  return null;
};

// Compatibility name for callers that only need one-rule matching.
export const ruleMatchesFeature = specificRuleMatchesFeature;
