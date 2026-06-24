const CLI_DEFAULT_FEATURES = Object.freeze([
  'CDS',
  'rRNA',
  'tRNA',
  'tmRNA',
  'ncRNA',
  'misc_RNA',
  'repeat_region'
]);

const DEFAULT_DIRECTIONAL_FEATURE_TYPES = new Set([
  'CDS',
  'rRNA',
  'tRNA',
  'tmRNA',
  'ncRNA',
  'misc_RNA'
]);

export const DEFAULT_LINEAR_BLAST_FILTERS = Object.freeze({
  bitscore: 50,
  evalue: '1e-2',
  identity: 0,
  alignment_length: 0
});

export const DEFAULT_CIRCULAR_CONSERVATION_BLAST_FILTERS = Object.freeze({
  bitscore: 50,
  evalue: '1e-5',
  identity: 70,
  alignment_length: 0
});

const normalizeFeatureList = (features) => {
  const values = Array.isArray(features) ? features : String(features || '').split(',');
  return values
    .map((featureType) => String(featureType || '').trim())
    .filter(Boolean);
};

const normalizeFeatureShape = (value) => (
  String(value || '').trim().toLowerCase() === 'arrow' ? 'arrow' : 'rectangle'
);

const defaultFeatureShape = (featureType) => (
  DEFAULT_DIRECTIONAL_FEATURE_TYPES.has(String(featureType || '').trim()) ? 'arrow' : 'rectangle'
);

const numericEquivalent = (left, right) => {
  const leftNumber = Number(left);
  const rightNumber = Number(right);
  if (Number.isFinite(leftNumber) && Number.isFinite(rightNumber)) {
    const tolerance = Number.EPSILON * Math.max(1, Math.abs(leftNumber), Math.abs(rightNumber)) * 8;
    return Math.abs(leftNumber - rightNumber) <= tolerance;
  }
  return String(left ?? '').trim() === String(right ?? '').trim();
};

const trimTrailingDefaults = (values, isDefault) => {
  const items = Array.isArray(values) ? [...values] : [];
  let length = items.length;
  while (length > 0 && isDefault(items[length - 1])) length -= 1;
  return items.slice(0, length);
};

const buildTrimmedRepeatedOptionArgs = (option, values, isDefault, formatValue = (value) => value) => (
  trimTrailingDefaults(values, isDefault).flatMap((value) => [option, formatValue(value)])
);

export const isCliDefaultFeatureList = (features) => {
  const normalized = normalizeFeatureList(features);
  return normalized.length === CLI_DEFAULT_FEATURES.length &&
    normalized.every((featureType, index) => featureType === CLI_DEFAULT_FEATURES[index]);
};

export const buildFeatureShapeAssignments = (features, featureShapes) => (
  normalizeFeatureList(features)
    .map((featureType) => {
      const shape = normalizeFeatureShape(featureShapes?.[featureType]);
      return shape === defaultFeatureShape(featureType) ? null : `${featureType}=${shape}`;
    })
    .filter(Boolean)
);

export const buildBlastFilterArgs = (filters, defaults = DEFAULT_LINEAR_BLAST_FILTERS) => {
  const args = [];
  const normalizedFilters = filters || {};
  if (!numericEquivalent(normalizedFilters.bitscore, defaults.bitscore)) {
    args.push('--bitscore', normalizedFilters.bitscore);
  }
  if (!numericEquivalent(normalizedFilters.evalue, defaults.evalue)) {
    args.push('--evalue', normalizedFilters.evalue);
  }
  if (!numericEquivalent(normalizedFilters.identity, defaults.identity)) {
    args.push('--identity', normalizedFilters.identity);
  }
  if (!numericEquivalent(normalizedFilters.alignment_length, defaults.alignment_length)) {
    args.push('--alignment_length', normalizedFilters.alignment_length);
  }
  return args;
};

export const buildRecordSelectorArgs = (selectors) => buildTrimmedRepeatedOptionArgs(
  '--record_id',
  selectors,
  (selector) => String(selector || '').trim() === '',
  (selector) => String(selector || '').trim()
);

export const buildReverseComplementArgs = (flags) => buildTrimmedRepeatedOptionArgs(
  '--reverse_complement',
  flags,
  (flag) => flag !== true,
  (flag) => (flag ? '1' : '0')
);
