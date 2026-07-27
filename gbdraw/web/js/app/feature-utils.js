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

export const firstFeatureText = (...values) => {
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

export const getFeatureQualifierFirstValue = (feature, key) => {
  const normalizedKey = String(key || '').trim().toLowerCase();
  if (!feature || !normalizedKey) return '';
  const qualifiers = feature.qualifiers && typeof feature.qualifiers === 'object' && !Array.isArray(feature.qualifiers)
    ? feature.qualifiers
    : {};
  if (Object.prototype.hasOwnProperty.call(qualifiers, normalizedKey)) {
    return firstFeatureText(qualifiers[normalizedKey]);
  }
  const matchingKey = Object.keys(qualifiers).find((candidate) => candidate.toLowerCase() === normalizedKey);
  return matchingKey ? firstFeatureText(qualifiers[matchingKey]) : '';
};

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

export const getFeatureQualifierValue = (feat, qual) => {
  if (!qual) return null;
  const key = qual.toLowerCase();
  if (feat.qualifiers && Object.prototype.hasOwnProperty.call(feat.qualifiers, key)) {
    return feat.qualifiers[key];
  }
  if (key === 'product') return feat.product || null;
  if (key === 'gene') return feat.gene || null;
  if (key === 'locus_tag') return feat.locus_tag || null;
  if (key === 'note') return feat.note || null;
  return null;
};

export const matchRuleValue = (value, ruleVal, strictEquals = false) => {
  if (!value || !ruleVal) return false;
  const values = Array.isArray(value) ? value : [value];
  const needle = String(ruleVal);
  for (const item of values) {
    if (item === null || item === undefined) continue;
    const text = String(item);
    try {
      const regex = new RegExp(needle, 'i');
      if (regex.test(text)) return true;
    } catch {
      if (strictEquals) {
        if (text === needle) return true;
      } else if (text.toLowerCase().includes(needle.toLowerCase())) {
        return true;
      }
    }
  }
  return false;
};

export const ruleMatchesFeature = (feat, rule) => {
  if (!rule || rule.feat !== feat.type) return false;
  const qualKey = (rule.qual || '').toLowerCase();

  if (qualKey === 'hash') {
    return getFeatureHashCandidates(feat).some((candidate) => matchRuleValue(candidate, rule.val, true));
  }

  if (qualKey === 'record_location') {
    const recordId = firstFeatureText(feat.record_id, feat.recordId, feat.record);
    const location = firstFeatureText(
      feat.selector?.location,
      feat.location,
      feat.start !== undefined && feat.end !== undefined ? `${feat.start}..${feat.end}` : ''
    );
    const strandToken = String(feat.strand ?? '').trim().toLowerCase();
    const strand = ['positive', 'plus', 'forward', '1'].includes(strandToken)
      ? '+'
      : (['negative', 'minus', 'reverse', '-1'].includes(strandToken) ? '-' : strandToken);
    const value = firstFeatureText(
      feat.selector?.record_location,
      feat.record_location,
      feat.recordLocation,
      recordId && location && strand ? `${recordId}:${location}:${strand}` : ''
    );
    return matchRuleValue(value, rule.val);
  }

  if (qualKey === 'location') {
    const location = `${feat.start}..${feat.end}`;
    return matchRuleValue(location, rule.val);
  }

  const value = getFeatureQualifierValue(feat, qualKey);
  if (!value) return false;
  const strictEquals = qualKey === 'locus_tag';
  return matchRuleValue(value, rule.val, strictEquals);
};
