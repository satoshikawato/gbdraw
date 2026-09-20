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
