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

export const DEFINITION_LINE_STYLE_KINDS = Object.freeze(['name', 'replicon', 'accession', 'length']);

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

const normalizeDefinitionLineFontSize = (value) => {
  if (value === null || value === undefined) return null;
  const raw = String(value).trim();
  if (!raw || ['auto', 'none', 'null', 'default'].includes(raw.toLowerCase())) return null;
  const numeric = Number(raw);
  return Number.isFinite(numeric) && numeric > 0 ? numeric : null;
};

const normalizeDefinitionLineFontWeight = (value) => {
  if (value === null || value === undefined) return null;
  const normalized = String(value).trim().toLowerCase();
  if (!normalized || ['auto', 'none', 'null', 'default', 'normal'].includes(normalized)) return null;
  if (['normal', 'bold', 'lighter', 'bolder'].includes(normalized)) return normalized;
  return /^(100|200|300|400|500|600|700|800|900)$/.test(normalized) ? normalized : null;
};

const normalizeDefinitionLineFill = (value) => {
  if (value === null || value === undefined) return null;
  const normalized = String(value).trim();
  return normalized || null;
};

export const createDefaultDefinitionLineStyle = () => ({
  font_size: null,
  font_weight: null,
  fill: null
});

export const createDefaultLinearDefinitionLineStyles = () => Object.fromEntries(
  DEFINITION_LINE_STYLE_KINDS.map((kind) => [kind, createDefaultDefinitionLineStyle()])
);

export const normalizeDefinitionLineStyleState = (source) => {
  const normalized = createDefaultLinearDefinitionLineStyles();
  if (!source || typeof source !== 'object' || Array.isArray(source)) return normalized;
  DEFINITION_LINE_STYLE_KINDS.forEach((kind) => {
    const entry = source[kind];
    if (!entry || typeof entry !== 'object' || Array.isArray(entry)) return;
    normalized[kind] = {
      font_size: normalizeDefinitionLineFontSize(entry.font_size),
      font_weight: normalizeDefinitionLineFontWeight(entry.font_weight),
      fill: normalizeDefinitionLineFill(entry.fill)
    };
  });
  return normalized;
};

export const buildDefinitionLineStyleAssignments = (source) => {
  const styles = normalizeDefinitionLineStyleState(source);
  const assignments = [];
  DEFINITION_LINE_STYLE_KINDS.forEach((kind) => {
    const style = styles[kind];
    const parts = [];
    if (style.font_size !== null) parts.push(`size=${style.font_size}`);
    if (style.font_weight !== null) parts.push(`weight=${style.font_weight}`);
    if (style.fill !== null) parts.push(`color=${style.fill}`);
    if (parts.length > 0) assignments.push(`${kind}:${parts.join(',')}`);
  });
  return assignments;
};

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
