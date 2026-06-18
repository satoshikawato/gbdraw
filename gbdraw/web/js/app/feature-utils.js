export const getFeatureCaption = (feature) => {
  return (
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
    if (!feat.svg_id) return false;
    return matchRuleValue(feat.svg_id, rule.val, true);
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
