export const FEATURE_SEMANTIC_SCOPE_QUALIFIER = '__gbdraw_semantic_scope__';
export const FEATURE_SEMANTIC_SCOPE_PREFIX = 'fs1:';
export const FEATURE_SEMANTIC_SCOPE_KINDS = Object.freeze([
  'feature-type',
  'base-legend-caption',
  'rendered-label',
  'source-annotation-label',
  'similarity-group'
]);

const text = (value) => String(value ?? '').trim();

const firstText = (...values) => {
  for (const value of values) {
    if (Array.isArray(value)) {
      const nested = firstText(...value);
      if (nested) return nested;
      continue;
    }
    const normalized = text(value);
    if (normalized) return normalized;
  }
  return '';
};

export const encodeFeatureSemanticSelector = (kind, value) => {
  const normalizedKind = text(kind);
  const normalizedValue = text(value);
  if (!FEATURE_SEMANTIC_SCOPE_KINDS.includes(normalizedKind)) {
    throw new Error(`Unsupported Feature semantic selector kind: ${normalizedKind}`);
  }
  if (!normalizedValue) throw new Error('Feature semantic selector values must be non-empty.');
  return `${FEATURE_SEMANTIC_SCOPE_PREFIX}${normalizedKind}:${encodeURIComponent(normalizedValue)}`;
};

export const parseFeatureSemanticSelector = (value) => {
  const literal = String(value ?? '');
  if (!literal.startsWith(FEATURE_SEMANTIC_SCOPE_PREFIX)) return null;
  const remainder = literal.slice(FEATURE_SEMANTIC_SCOPE_PREFIX.length);
  const separator = remainder.indexOf(':');
  if (separator <= 0) return null;
  const kind = remainder.slice(0, separator);
  let decoded;
  try {
    decoded = decodeURIComponent(remainder.slice(separator + 1));
  } catch {
    return null;
  }
  if (!FEATURE_SEMANTIC_SCOPE_KINDS.includes(kind) || !text(decoded)) return null;
  try {
    if (encodeFeatureSemanticSelector(kind, decoded) !== literal) return null;
  } catch {
    return null;
  }
  return Object.freeze({ kind, value: decoded });
};

export const isFeatureSemanticSelector = (value) => (
  parseFeatureSemanticSelector(value) !== null
);

export const getFeatureSourceAnnotation = (feature) => {
  const explicitLabel = firstText(
    feature?.sourceAnnotationLabel,
    feature?.source_annotation_label
  );
  const explicitQualifier = firstText(
    feature?.sourceAnnotationQualifier,
    feature?.source_annotation_qualifier
  );
  if (explicitLabel) {
    return Object.freeze({ label: explicitLabel, qualifier: explicitQualifier || 'label' });
  }
  const qualifiers = feature?.qualifiers && typeof feature.qualifiers === 'object'
    && !Array.isArray(feature.qualifiers)
    ? feature.qualifiers
    : {};
  for (const qualifier of ['product', 'gene', 'locus_tag', 'note']) {
    const key = Object.keys(qualifiers).find(
      (candidate) => candidate.toLowerCase() === qualifier
    );
    const label = firstText(feature?.[qualifier], key ? qualifiers[key] : '');
    if (label) return Object.freeze({ label, qualifier });
  }
  return Object.freeze({ label: '', qualifier: '' });
};

export const getFeatureRenderedLabel = (feature) => firstText(
  feature?.renderedLabel,
  feature?.rendered_label,
  feature?.displayLabel,
  feature?.display_label,
  feature?.labelText,
  feature?.label_text,
  feature?.label
);

export const getFeatureSimilarityGroupIds = (feature) => [...new Set([
  feature?.orthogroupId,
  feature?.orthogroup_id,
  feature?.orthogroupIds,
  feature?.orthogroup_ids,
  feature?.similarityGroupId,
  feature?.similarity_group_id
].flatMap((value) => (
  Array.isArray(value) ? value : String(value ?? '').split(';')
)).map(text).filter(Boolean))].sort();

export const featureSemanticSelectorMatches = (
  feature,
  rawSelector,
  { baseLegendCaption = '' } = {}
) => {
  const selector = rawSelector && typeof rawSelector === 'object'
    ? rawSelector
    : parseFeatureSemanticSelector(rawSelector);
  if (!selector) return false;
  if (selector.kind === 'feature-type') return text(feature?.type) === selector.value;
  if (selector.kind === 'source-annotation-label') {
    return getFeatureSourceAnnotation(feature).label === selector.value;
  }
  if (selector.kind === 'base-legend-caption') {
    return text(baseLegendCaption).toLowerCase() === selector.value.toLowerCase();
  }
  if (selector.kind === 'rendered-label') {
    return getFeatureRenderedLabel(feature) === selector.value;
  }
  if (selector.kind === 'similarity-group') {
    return getFeatureSimilarityGroupIds(feature).includes(selector.value);
  }
  return false;
};

export const isFeatureSemanticScopeRule = (rule) => (
  text(rule?.qual) === FEATURE_SEMANTIC_SCOPE_QUALIFIER
  && text(rule?.match).toLowerCase() !== 'regex'
  && isFeatureSemanticSelector(rule?.val)
);
