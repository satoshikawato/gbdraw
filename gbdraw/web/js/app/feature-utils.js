export const getFeatureCaption = (feature) => {
  return (
    feature.product ||
    feature.gene ||
    feature.locus_tag ||
    feature.note ||
    `${feature.type} at ${feature.start}..${feature.end}`
  );
};

export const getFeatureQualifierValue = (feat, qual) => {
  if (!qual) return null;
  const key = qual.toLowerCase();
  if (key === 'product') return feat.product || null;
  if (key === 'gene') return feat.gene || null;
  if (key === 'locus_tag') return feat.locus_tag || null;
  if (key === 'note') return feat.note || null;
  return feat.qualifiers && Object.prototype.hasOwnProperty.call(feat.qualifiers, key)
    ? feat.qualifiers[key]
    : null;
};

export const matchRuleValue = (value, ruleVal, strictEquals = false) => {
  if (!value || !ruleVal) return false;
  try {
    const regex = new RegExp(ruleVal, 'i');
    return regex.test(value);
  } catch {
    if (strictEquals) return value === ruleVal;
    return value.toLowerCase().includes(String(ruleVal).toLowerCase());
  }
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
