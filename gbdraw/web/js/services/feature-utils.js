// @ts-check
import { classifyOptionalNumber } from '../utils/optional-positive-number.js';

/**
 * A feature or Orthogroup member row: Python and JSON fields read by name.
 * @typedef {{
 *   [field: string]: unknown,
 *   selector?: FeatureSelector | null,
 *   qualifiers?: unknown,
 *   location_parts?: unknown
 * }} FeatureLike
 */
/**
 * The selector a Python feature row carries beside its fields.
 * @typedef {{
 *   [field: string]: unknown,
 *   hash?: unknown,
 *   location?: unknown,
 *   record_location?: unknown,
 *   qualifiers?: unknown
 * }} FeatureSelector
 */

/**
 * @param {unknown} value
 * @returns {string[]}
 */
export const normalizeStringArray = (value) => {
  if (Array.isArray(value)) {
    return value
      .filter((item) => item !== null && item !== undefined)
      .map((item) => String(item));
  }
  if (value === null || value === undefined || value === '') return [];
  return [String(value)];
};

/** @param {FeatureLike | null | undefined} feature */
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
  return caption || `${feature?.type} at ${formatFeatureRange(feature)}`;
};

// Python owns `location_parts[].display` (1-based INSDC). A feature without
// parts is one part spanning its own interval.
/**
 * @param {FeatureLike | null | undefined} feature
 * @returns {(FeatureLike | null | undefined)[]}
 */
const featureLocationParts = (feature) => (
  Array.isArray(feature?.location_parts) && feature.location_parts.length
    ? feature.location_parts
    : [feature]
);

const NON_ASCII = /[^\x00-\x7F]/;

/**
 * The one reader of a start or end coordinate (PT-OD14, PT-OD23): a number or
 * ASCII number text, truncated toward zero. Anything else is missing (null):
 * a missing, null or blank value, a boolean, other text (other digits,
 * underscores, hex, any non-ASCII character), a non-finite number, and a
 * truncation beyond the safe integers. A reader must not take a missing bound
 * as 0 (OV-968, OV-1010). The twin of `_location_part_number` in
 * gbdraw/render/interactive_svg.py (tests/fixtures/location-parts-vectors.json);
 * the Interactive SVG runtime cannot import it and keeps one copy
 * (`locationPartBound` in standalone-interactivity-assets.js).
 * @param {unknown} value
 * @returns {number | null}
 */
export const locationPartNumber = (value) => {
  // Any non-ASCII text is missing, spaces included: Python and JavaScript strip different Unicode spaces.
  if (typeof value === 'string' && NON_ASCII.test(value)) return null;
  const number = classifyOptionalNumber(value).value ?? null;
  if (number === null) return null;
  // `|| 0` turns the -0 of a negative fraction into 0, as Python's int() does.
  const truncated = Math.trunc(number) || 0;
  return Number.isSafeInteger(truncated) ? truncated : null;
};

/**
 * `locationPartNumber` of a part's field, NaN when missing.
 * @param {unknown} part A location part, a feature read as its one part, or a
 *   Similarity Group member row.
 * @param {'start' | 'end'} field
 * @returns {number}
 */
export const locationPartBound = (part, field) => (
  locationPartNumber(part && typeof part === 'object' ? /** @type {Record<string, unknown>} */ (part)[field] : undefined)
  ?? NaN
);

// The 1-based coordinates of every part, without the strand. Also the location
// of the fallback caption of an unnamed feature, here and in the rule editor.
/** @param {FeatureLike | null | undefined} feature */
export const formatFeatureRange = (feature) => featureLocationParts(feature).map((part) => {
  const display = String(part?.display || '').trim();
  if (display) return display;
  const start = locationPartBound(part, 'start');
  const end = locationPartBound(part, 'end');
  return Number.isFinite(start) && Number.isFinite(end) ? `${start + 1}..${end}` : '';
}).filter(Boolean).join(', ');

/** @param {FeatureLike | null | undefined} feature */
export const formatFeatureLocation = (feature) => {
  const range = formatFeatureRange(feature);
  const strand = String(feature?.strand || '').trim();
  return range && strand ? `${range} (${strand})` : range;
};

/** @param {FeatureLike | null | undefined} feature */
export const formatFeatureLength = (feature) => {
  let length = 0;
  for (const part of featureLocationParts(feature)) {
    const start = locationPartBound(part, 'start');
    const end = locationPartBound(part, 'end');
    if (!Number.isFinite(start) || !Number.isFinite(end) || end < start) return '';
    length += end - start;
  }
  return `${Math.round(length).toLocaleString()} bp`;
};

/**
 * @param {...unknown} values
 * @returns {string}
 */
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
// Python `instance_svg_id` appends `__instance_<instance>_<digest>` to a
// duplicated record's or feature's rendered id.
const RENDERED_INSTANCE_SUFFIX = /__instance_.*$/;

/** @param {FeatureLike} feature */
const renderedFeatureId = (feature) => firstFeatureText(
  feature.rendered_svg_id,
  feature.renderedSvgId,
  feature.rendered_feature_svg_id,
  feature.renderedFeatureSvgId,
  feature.svg_id,
  feature.svgId
);

/** @param {FeatureLike | null | undefined} feature */
export const getFeatureHashCandidates = (feature) => {
  if (!feature || typeof feature !== 'object') return [];

  const renderedId = renderedFeatureId(feature);
  const generationId = renderedId.replace(LINEAR_RENDERED_RECORD_SUFFIX, '');

  return [...new Set([generationId, renderedId].filter(Boolean))];
};

// The hash in a feature's rendered ID: no record or instance suffix. A feature
// catalog before the one feature-hash builder (OV-401) drew a cropped or
// reverse-complemented record with its drawn location's hash here.
/** @param {FeatureLike | null | undefined} feature */
export const getRenderedFeatureHash = (feature) => (
  !feature || typeof feature !== 'object'
    ? ''
    : renderedFeatureId(feature)
      .replace(RENDERED_INSTANCE_SUFFIX, '')
      .replace(LINEAR_RENDERED_RECORD_SUFFIX, '')
);

// The hash Python's `hash` rule qualifier matches (`get_feature_hash`): the
// feature's source-record hash, which duplicates of one record share
// (PD-OI-069), whatever crop or orientation its record is drawn with.
/** @param {FeatureLike | null | undefined} feature */
export const getFeatureColorRuleHash = (feature) => (
  !feature || typeof feature !== 'object'
    ? ''
    : firstFeatureText(feature.selector?.hash) || getRenderedFeatureHash(feature)
);

/**
 * @param {FeatureLike | null | undefined} feature
 * @param {...string} keys
 * @returns {string}
 */
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

/** @param {unknown} value */
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

/**
 * @param {...unknown} values
 * @returns {string}
 */
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

/**
 * @param {FeatureLike | null | undefined} feature
 * @param {...string} keys
 * @returns {string}
 */
const directFeatureDisplayValue = (feature, ...keys) => {
  if (!feature || typeof feature !== 'object') return '';
  for (const key of keys) {
    const value = firstNonInternalProteinDisplayText(feature[key]);
    if (value) return value;
  }
  return '';
};

/**
 * @param {FeatureLike | null | undefined} feature
 * @param {unknown} key
 * @returns {string}
 */
const getFeatureQualifierDisplayValue = (feature, key) => {
  const normalizedKey = String(key || '').trim().toLowerCase();
  if (!feature || !normalizedKey) return '';
  /** @type {{ [name: string]: unknown }} */
  const qualifiers = feature.qualifiers && typeof feature.qualifiers === 'object' && !Array.isArray(feature.qualifiers)
    ? /** @type {{ [name: string]: unknown }} */ (feature.qualifiers) // an object of qualifier values, checked just above
    : {};
  if (Object.prototype.hasOwnProperty.call(qualifiers, normalizedKey)) {
    return firstNonInternalProteinDisplayText(qualifiers[normalizedKey]);
  }
  const matchingKey = Object.keys(qualifiers).find((candidate) => candidate.toLowerCase() === normalizedKey);
  return matchingKey
    ? firstNonInternalProteinDisplayText(qualifiers[matchingKey])
    : '';
};

/**
 * @param {FeatureLike | null | undefined} feature
 * @param {FeatureLike | null} [member]
 * @param {string} [fallback]
 * @returns {string}
 */
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

// A feature's own protein ID, for a row labeled Protein ID: unlike
// resolveDisplayProteinId, never its locus tag, gene, or other identifier.
// The Interactive SVG mirrors it (featureProteinId).
/**
 * @param {Record<string, any> | null | undefined} feature
 * @param {Record<string, any> | null} [member]
 * @returns {string}
 */
export const resolveFeatureProteinId = (feature, member = null) =>
  firstNonInternalProteinDisplayText(
    directFeatureDisplayValue(feature, 'displayProteinId', 'display_protein_id'),
    directFeatureDisplayValue(member, 'displayProteinId', 'display_protein_id'),
    directFeatureDisplayValue(feature, 'sourceProteinId', 'source_protein_id'),
    directFeatureDisplayValue(member, 'sourceProteinId', 'source_protein_id'),
    getFeatureQualifierDisplayValue(feature, 'protein_id'),
    directFeatureDisplayValue(feature, 'proteinId', 'protein_id'),
    directFeatureDisplayValue(member, 'proteinId', 'protein_id')
  );

/**
 * @param {FeatureLike | null | undefined} feature
 * @param {FeatureLike | null} [member]
 * @param {string} [fallback]
 * @returns {string}
 */
export const resolveInternalProteinId = (feature, member = null, fallback = '') => firstFeatureText(
  directFeatureValue(feature, 'proteinId', 'protein_id'),
  directFeatureValue(member, 'proteinId', 'protein_id'),
  fallback
);
