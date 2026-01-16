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
