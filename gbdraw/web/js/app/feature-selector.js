const FEATURE_ID_KEY_EMPTY = '';

export const SPECIFIC_COLOR_QUALIFIER_PRESETS = Object.freeze([
  'product',
  'gene',
  'locus_tag',
  'protein_id',
  'note',
  'function',
  'phrog',
  'vfdb_short_name',
  'AMR_Gene_Family',
  'color'
]);

export const escapeRegexLiteral = (value) => String(value ?? '').replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
export const exactRegexValue = (value) => `^${escapeRegexLiteral(value)}$`;

export const normalizeSelectorText = (value) => String(value ?? '').trim();
export const normalizeFeatureIdKey = (value) => normalizeSelectorText(value).toLowerCase();

export const collectSpecificColorQualifierSuggestions = (features = [], rules = []) => {
  const suggestions = new Set(SPECIFIC_COLOR_QUALIFIER_PRESETS);
  const addQualifier = (value) => {
    const qualifier = normalizeSelectorText(value);
    if (qualifier) suggestions.add(qualifier);
  };
  const addQualifierMap = (qualifiers) => {
    if (!qualifiers || typeof qualifiers !== 'object' || Array.isArray(qualifiers)) return;
    Object.keys(qualifiers).forEach(addQualifier);
  };

  if (Array.isArray(features)) {
    features.forEach((feature) => {
      addQualifierMap(feature?.qualifiers);
      addQualifierMap(feature?.selector?.qualifiers);
    });
  }
  if (Array.isArray(rules)) {
    rules.forEach((rule) => addQualifier(rule?.qual));
  }

  const presets = new Set(SPECIFIC_COLOR_QUALIFIER_PRESETS);
  const custom = [...suggestions]
    .filter((qualifier) => !presets.has(qualifier))
    .sort((left, right) => left.localeCompare(right));
  return [...SPECIFIC_COLOR_QUALIFIER_PRESETS, ...custom];
};

const firstText = (...values) => {
  for (const value of values) {
    if (Array.isArray(value)) {
      const found = firstText(...value);
      if (found) return found;
      continue;
    }
    const normalized = normalizeSelectorText(value);
    if (normalized) return normalized;
  }
  return '';
};

const normalizeStrandToken = (value) => {
  const token = normalizeSelectorText(value).toLowerCase();
  if (token === '+' || token === 'positive' || token === 'plus' || token === 'forward' || token === '1') return '+';
  if (token === '-' || token === 'negative' || token === 'minus' || token === 'reverse' || token === '-1') return '-';
  return 'undefined';
};

export const normalizeQualifierMap = (qualifiers) => {
  const normalized = {};
  if (!qualifiers || typeof qualifiers !== 'object') return normalized;
  Object.entries(qualifiers).forEach(([keyRaw, valuesRaw]) => {
    const key = normalizeSelectorText(keyRaw).toLowerCase();
    if (!key) return;
    const values = Array.isArray(valuesRaw) ? valuesRaw : [valuesRaw];
    normalized[key] = values
      .filter((value) => value !== null && value !== undefined)
      .map((value) => String(value));
  });
  return normalized;
};

const hasOwnValues = (obj) => Object.keys(obj || {}).length > 0;

const recordLocationPosition = (recordId, recordLocation) => {
  const record = normalizeSelectorText(recordId);
  const value = normalizeSelectorText(recordLocation);
  if (!record || !value) return '';
  const prefix = `${record}:`;
  return value.startsWith(prefix) ? value.slice(prefix.length) : '';
};

export const normalizeFeatureSelectorMetadata = (feature, options = {}) => {
  const requireSelector = options.requireSelector === true;
  const preferSelector = options.preferSelector !== false;
  const selector = feature?.selector && typeof feature.selector === 'object' ? feature.selector : null;
  const hasSelectorMetadata = Boolean(selector);
  const selectorQualifiers = normalizeQualifierMap(selector?.qualifiers);
  const fallbackQualifiers = requireSelector ? {} : normalizeQualifierMap(feature?.qualifiers);
  const qualifiers = preferSelector || requireSelector || hasOwnValues(selectorQualifiers)
    ? selectorQualifiers
    : fallbackQualifiers;

  const recordId = firstText(
    feature?.record_id,
    feature?.recordId,
    feature?.record,
    feature?.record_id_text
  );
  const featureType = firstText(
    feature?.type,
    feature?.featureType,
    feature?.feature_type
  );
  const featureId = firstText(
    feature?.svg_id,
    feature?.svgId,
    selector?.hash,
    feature?.featureId,
    feature?.feature_id,
    feature?.id
  );
  const stableFeatureId = firstText(
    selector?.hash,
    feature?.stable_svg_id,
    feature?.stableSvgId,
    feature?.stableFeatureSvgId,
    feature?.stable_feature_id,
    feature?.stableFeatureId,
    feature?.feature_hash,
    feature?.hash,
    featureId
  );
  const location = firstText(
    selector?.location,
    feature?.location,
    feature?.start !== undefined && feature?.end !== undefined
      ? `${feature.start}..${feature.end}`
      : ''
  );
  const recordLocation = firstText(selector?.record_location, feature?.recordLocation);
  const position = firstText(
    feature?.position,
    recordLocationPosition(recordId, recordLocation),
    location ? `${location}:${normalizeStrandToken(feature?.strand)}` : ''
  );

  return {
    featureId,
    featureIdKey: normalizeFeatureIdKey(featureId) || FEATURE_ID_KEY_EMPTY,
    stableFeatureId,
    stableFeatureIdKey: normalizeFeatureIdKey(stableFeatureId) || FEATURE_ID_KEY_EMPTY,
    recordId,
    record: recordId,
    featureType,
    location,
    position,
    recordLocation,
    qualifiers,
    hasSelectorMetadata
  };
};

export const buildFeatureMetadataMap = (features, options = {}) => {
  const metadataByFeatureId = new Map();
  if (!Array.isArray(features)) return metadataByFeatureId;

  features.forEach((feature) => {
    const metadata = normalizeFeatureSelectorMetadata(feature, options);
    const key = normalizeFeatureIdKey(metadata.featureId);
    if (!key || metadataByFeatureId.has(key)) return;
    metadataByFeatureId.set(key, metadata);
  });

  return metadataByFeatureId;
};

export const makeSelectorUniquenessKey = (recordId, featureType, qualifier, value) =>
  `${normalizeSelectorText(recordId)}\u0000${normalizeSelectorText(featureType)}\u0000` +
  `${normalizeSelectorText(qualifier).toLowerCase()}\u0000${normalizeSelectorText(value)}`;

const markSelectorSafetyScopeAvailability = (counts, available) => {
  Object.defineProperty(counts, 'selectorSafetyScopeAvailable', {
    value: Boolean(available),
    enumerable: false,
    configurable: true
  });
  return counts;
};

const incrementSelectorCount = (counts, recordId, featureType, qualifier, value) => {
  const record = normalizeSelectorText(recordId);
  const type = normalizeSelectorText(featureType);
  const normalizedQualifier = normalizeSelectorText(qualifier).toLowerCase();
  const normalizedValue = normalizeSelectorText(value);
  if (!record || !type || !normalizedQualifier || !normalizedValue) return;
  const key = makeSelectorUniquenessKey(record, type, normalizedQualifier, normalizedValue);
  counts.set(key, (counts.get(key) || 0) + 1);
};

const addMetadataToSelectorCounts = (counts, metadata) => {
  const recordId = normalizeSelectorText(metadata?.recordId || metadata?.record);
  const featureType = normalizeSelectorText(metadata?.featureType);
  if (!recordId || !featureType) return;

  const seenForFeature = new Set();
  Object.entries(metadata?.qualifiers || {}).forEach(([qualifierRaw, valuesRaw]) => {
    const qualifier = normalizeSelectorText(qualifierRaw).toLowerCase();
    if (!qualifier) return;
    const values = Array.isArray(valuesRaw) ? valuesRaw : [valuesRaw];
    values.forEach((valueRaw) => {
      const value = normalizeSelectorText(valueRaw);
      if (!value) return;
      const key = makeSelectorUniquenessKey(recordId, featureType, qualifier, value);
      if (seenForFeature.has(key)) return;
      seenForFeature.add(key);
      incrementSelectorCount(counts, recordId, featureType, qualifier, value);
    });
  });

  const recordLocation = normalizeSelectorText(metadata?.recordLocation || getRecordLocationFromMetadata(metadata));
  if (recordLocation) {
    const key = makeSelectorUniquenessKey(recordId, featureType, 'record_location', recordLocation);
    if (!seenForFeature.has(key)) {
      seenForFeature.add(key);
      incrementSelectorCount(counts, recordId, featureType, 'record_location', recordLocation);
    }
  }
};

export const getRecordLocationFromMetadata = (metadata) => {
  const existing = normalizeSelectorText(metadata?.recordLocation);
  if (existing) return existing;
  const record = normalizeSelectorText(metadata?.recordId || metadata?.record);
  const position = normalizeSelectorText(metadata?.position);
  if (!record || !position) return '';
  return `${record}:${position}`;
};

export const buildFeatureSelectorUniquenessIndexFromMetadata = (metadataByFeatureId) => {
  const counts = markSelectorSafetyScopeAvailability(new Map(), false);
  if (!(metadataByFeatureId instanceof Map)) return counts;
  metadataByFeatureId.forEach((metadata) => addMetadataToSelectorCounts(counts, metadata));
  return counts;
};

export const buildFeatureSelectorUniquenessIndex = (features, options = {}) => {
  const metadataByFeatureId = buildFeatureMetadataMap(features, options);
  return buildFeatureSelectorUniquenessIndexFromMetadata(metadataByFeatureId);
};

export const buildSelectorSafetyUniquenessIndex = (selectorSafetyScope) => {
  const counts = markSelectorSafetyScopeAvailability(new Map(), Array.isArray(selectorSafetyScope));
  if (!Array.isArray(selectorSafetyScope)) return counts;

  selectorSafetyScope.forEach((entry) => {
    const metadata = normalizeFeatureSelectorMetadata(
      {
        ...entry,
        type: entry?.feature_type ?? entry?.featureType ?? entry?.type,
        record_id: entry?.record_id ?? entry?.recordId ?? entry?.record,
        selector: entry?.selector
      },
      { requireSelector: true }
    );
    addMetadataToSelectorCounts(counts, metadata);
  });
  return counts;
};

export const hasSelectorSafetyScope = (uniquenessIndex) =>
  uniquenessIndex instanceof Map && uniquenessIndex.selectorSafetyScopeAvailable === true;

export const selectFeatureSelector = (featureMeta, uniquenessIndex, options = {}) => {
  const priority = Array.isArray(options.priority) ? options.priority : [];
  const requireSelector = options.requireSelector === true;
  const requireSafetyScope = options.requireSafetyScope === true;
  const allowRecordLocation = options.allowRecordLocation !== false;
  const metadata = normalizeFeatureSelectorMetadata(featureMeta, {
    requireSelector,
    preferSelector: options.preferSelector !== false
  });
  const fallbackHashValue = metadata.stableFeatureId || metadata.featureId;
  const fallback = {
    qualifier: 'hash',
    value: fallbackHashValue,
    isFallbackHash: true
  };

  if (!metadata.featureId && !fallbackHashValue) return fallback;
  if (requireSelector && !metadata.hasSelectorMetadata) return fallback;
  if (requireSafetyScope && !hasSelectorSafetyScope(uniquenessIndex)) return fallback;

  const recordId = normalizeSelectorText(metadata.recordId || metadata.record);
  const featureType = normalizeSelectorText(metadata.featureType);
  if (!recordId || !featureType || !(uniquenessIndex instanceof Map)) return fallback;

  for (const qualifierRaw of priority) {
    const qualifier = normalizeSelectorText(qualifierRaw).toLowerCase();
    if (!qualifier) continue;
    const values = Array.isArray(metadata.qualifiers?.[qualifier]) ? metadata.qualifiers[qualifier] : [];
    for (const valueRaw of values) {
      const value = normalizeSelectorText(valueRaw);
      if (!value) continue;
      const key = makeSelectorUniquenessKey(recordId, featureType, qualifier, value);
      if ((uniquenessIndex.get(key) || 0) === 1) {
        return {
          qualifier,
          value,
          isFallbackHash: false
        };
      }
    }
  }

  const recordLocation = getRecordLocationFromMetadata(metadata);
  if (allowRecordLocation && recordLocation) {
    const key = makeSelectorUniquenessKey(recordId, featureType, 'record_location', recordLocation);
    if ((uniquenessIndex.get(key) || 0) === 1) {
      return {
        qualifier: 'record_location',
        value: recordLocation,
        isFallbackHash: false
      };
    }
  }

  return fallback;
};
