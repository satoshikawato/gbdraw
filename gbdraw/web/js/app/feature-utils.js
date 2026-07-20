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
  return (
    feature.label ||
    feature.display_label ||
    feature.displayLabel ||
    feature.product ||
    feature.gene ||
    feature.locus_tag ||
    feature.note ||
    `${feature.type} at ${feature.start}..${feature.end}`
  );
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

export const resolveDisplayProteinId = (feature, member = null, fallback = '') => firstFeatureText(
  directFeatureValue(feature, 'sourceProteinId', 'source_protein_id'),
  directFeatureValue(member, 'sourceProteinId', 'source_protein_id'),
  getFeatureQualifierFirstValue(feature, 'protein_id'),
  directFeatureValue(feature, 'locusTag', 'locus_tag'),
  getFeatureQualifierFirstValue(feature, 'locus_tag'),
  directFeatureValue(member, 'locusTag', 'locus_tag'),
  directFeatureValue(feature, 'geneId', 'gene_id'),
  getFeatureQualifierFirstValue(feature, 'gene_id'),
  directFeatureValue(member, 'geneId', 'gene_id'),
  directFeatureValue(feature, 'oldLocusTag', 'old_locus_tag'),
  getFeatureQualifierFirstValue(feature, 'old_locus_tag'),
  directFeatureValue(member, 'oldLocusTag', 'old_locus_tag'),
  directFeatureValue(feature, 'ID'),
  getFeatureQualifierFirstValue(feature, 'ID'),
  directFeatureValue(feature, 'Name', 'name'),
  getFeatureQualifierFirstValue(feature, 'Name'),
  directFeatureValue(feature, 'Parent', 'parent'),
  getFeatureQualifierFirstValue(feature, 'Parent'),
  directFeatureValue(feature, 'gene'),
  getFeatureQualifierFirstValue(feature, 'gene'),
  directFeatureValue(member, 'gene'),
  directFeatureValue(feature, 'proteinId', 'protein_id'),
  directFeatureValue(member, 'proteinId', 'protein_id'),
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
