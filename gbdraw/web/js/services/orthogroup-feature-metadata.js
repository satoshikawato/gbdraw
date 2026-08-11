import { normalizeGroupMetadataScope } from '../app/losat-normalization.js';
import {
  RECORD_INDEX_KEYS,
  RENDERED_FEATURE_ID_KEYS,
  SOURCE_FEATURE_INDEX_KEYS,
  nonnegativeIntegerAliasStatus,
  textAliasStatus,
  uniqueOrthogroupEntries
} from './feature-identity.js';

const normalizeText = (value) => String(value || '').trim();

const normalizeRecordIndex = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const text = String(value).trim();
  if (!/^\d+$/.test(text) || typeof value === 'boolean') return null;
  const normalized = Number(text);
  return Number.isSafeInteger(normalized) ? normalized : null;
};

const orthogroupFeatureIndexKey = (recordIndex, featureSvgId) => (
  `record-id:${Number(recordIndex)}:${normalizeText(featureSvgId)}`
);

const orthogroupRenderedFeatureIndexKey = (recordIndex, featureSvgId) => (
  `record-rendered-id:${Number(recordIndex)}:${normalizeText(featureSvgId)}`
);

const orthogroupFeatureIndexSourceKey = (recordIndex, featureIndex) => (
  `record-feature-index:${Number(recordIndex)}:${Number(featureIndex)}`
);

const orthogroupFeatureCanonicalKey = (recordKey, biologicalFeatureId) => (
  `canonical:${normalizeText(recordKey)}\u0000${normalizeText(biologicalFeatureId)}`
);

const consistentTextAlias = (payload, keys) => {
  const status = textAliasStatus(payload, keys);
  return {
    valid: status.valid,
    value: status.value
  };
};

const consistentIntegerAlias = (payload, keys, fallback = null) => {
  const fallbackKey = '__gbdrawFallbackIdentity';
  const source = fallback !== null && fallback !== undefined && fallback !== ''
    ? { ...(payload || {}), [fallbackKey]: fallback }
    : payload;
  const status = nonnegativeIntegerAliasStatus(source, [
    ...keys,
    ...(source === payload ? [] : [fallbackKey])
  ]);
  return {
    valid: status.valid,
    value: status.value
  };
};

const sourceIdentity = (payload) => {
  const stable = consistentTextAlias(payload, [
    'stableFeatureSvgId', 'stable_feature_svg_id', 'stableFeatureId',
    'stable_feature_id', 'stableSvgId', 'stable_svg_id'
  ]);
  const legacy = consistentTextAlias(payload, ['featureSvgId', 'feature_svg_id']);
  return {
    valid: stable.valid && legacy.valid && (
      !stable.value || !legacy.value || stable.value === legacy.value
    ),
    value: stable.value || legacy.value
  };
};

const renderedIdentity = (payload) => consistentTextAlias(payload, [
  ...RENDERED_FEATURE_ID_KEYS
]);

const memberIdentity = (payload) => {
  const source = sourceIdentity(payload);
  const rendered = renderedIdentity(payload);
  return { valid: source.valid && rendered.valid, source, rendered };
};

const featureIdentity = (payload) => {
  const source = sourceIdentity(payload);
  const rendered = renderedIdentity(payload);
  const svg = consistentTextAlias(payload, ['svgId', 'svg_id']);
  if (!source.valid || !rendered.valid || !svg.valid) {
    return { valid: false, source, rendered };
  }
  if (svg.value && rendered.value && svg.value !== rendered.value) {
    return { valid: false, source, rendered };
  }
  if (svg.value && !rendered.value && source.value && svg.value !== source.value) {
    rendered.value = svg.value;
  } else if (svg.value && !rendered.value && !source.value) {
    source.value = svg.value;
  }
  return { valid: true, source, rendered };
};

const canonicalIdentity = (payload) => {
  const record = consistentTextAlias(payload, ['recordKey', 'record_key']);
  const feature = consistentTextAlias(
    payload,
    ['biologicalFeatureId', 'biological_feature_id']
  );
  return {
    valid: record.valid && feature.valid && Boolean(record.value) === Boolean(feature.value),
    value: record.value && feature.value ? [record.value, feature.value] : null
  };
};

const sourceFeatureIndex = (payload) => consistentIntegerAlias(
  payload,
  SOURCE_FEATURE_INDEX_KEYS
);

const recordIdentity = (payload, fallback = null) => consistentIntegerAlias(
  payload,
  RECORD_INDEX_KEYS,
  fallback
);

const ORTHOGROUP_FEATURE_FIELDS = Object.freeze([
  'orthogroupId',
  'orthogroup_id',
  'orthogroupMemberCount',
  'orthogroup_member_count',
  'orthogroupRecordCoverage',
  'orthogroup_record_coverage',
  'orthogroupRepresentative',
  'orthogroup_representative',
  'orthogroupScope',
  'orthogroup_scope',
  'orthogroupSourceRecordIndex',
  'orthogroup_source_record_index',
  'orthogroupMember',
  'orthogroup_member'
]);

const clearDerivedOrthogroupMetadata = (feature) => {
  if (!feature || !ORTHOGROUP_FEATURE_FIELDS.some((field) => (
    Object.prototype.hasOwnProperty.call(feature, field)
  ))) return feature;
  const cleared = { ...feature };
  ORTHOGROUP_FEATURE_FIELDS.forEach((field) => delete cleared[field]);
  return cleared;
};

export const buildOrthogroupFeatureIndex = (orthogroups) => {
  const index = new Map();
  const owners = new Map();
  const ambiguousKeys = new Set();
  const recordsWithRenderedIds = new Set();
  const addUniqueEntry = (key, owner, entry) => {
    if (!key || ambiguousKeys.has(key)) return;
    const existingOwner = owners.get(key);
    if (existingOwner && existingOwner !== owner) {
      index.delete(key);
      ambiguousKeys.add(key);
      return;
    }
    owners.set(key, owner);
    index.set(key, entry);
  };

  uniqueOrthogroupEntries(orthogroups).forEach(({ group, id: orthogroupId }) => {
    const members = Array.isArray(group?.members) ? group.members : [];
    const memberCount = Number(group?.member_count || members.length || 0);
    const recordCoverage = Number(group?.record_coverage_count || new Set(
      members
        .map((member) => normalizeRecordIndex(member?.recordIndex ?? member?.record_index))
        .filter((recordIndex) => recordIndex !== null)
    ).size || 0);
    const orthogroupScope = normalizeGroupMetadataScope(group?.scope);
    const sourceRecordIndex = normalizeRecordIndex(group?.source_record_index);

    members.forEach((member) => {
      const record = recordIdentity(member);
      const featureIndex = sourceFeatureIndex(member);
      const identity = memberIdentity(member);
      const canonical = canonicalIdentity(member);
      const hasRecordIdentity = Boolean(
        identity.source.value || identity.rendered.value
      ) || featureIndex.value !== null;
      if (
        !record.valid
        || !featureIndex.valid
        || !identity.valid
        || !canonical.valid
        || (hasRecordIdentity && record.value === null)
        || (!hasRecordIdentity && canonical.value === null)
      ) return;
      const entry = {
        orthogroupId,
        orthogroupMemberCount: memberCount,
        orthogroupRecordCoverage: recordCoverage,
        proteinId: normalizeText(member?.proteinId),
        sourceProteinId: normalizeText(member?.sourceProteinId),
        orthogroupRepresentative: Boolean(member?.representative),
        orthogroupScope,
        orthogroupSourceRecordIndex: sourceRecordIndex,
        orthogroupMember: member
      };
      const owner = {};
      if (identity.source.value) {
        addUniqueEntry(
          orthogroupFeatureIndexKey(record.value, identity.source.value),
          owner,
          entry
        );
      }
      if (identity.rendered.value) {
        addUniqueEntry(
          orthogroupRenderedFeatureIndexKey(record.value, identity.rendered.value),
          owner,
          entry
        );
        recordsWithRenderedIds.add(record.value);
      }
      if (featureIndex.value !== null) {
        addUniqueEntry(
          orthogroupFeatureIndexSourceKey(record.value, featureIndex.value),
          owner,
          entry
        );
      }
      if (canonical.value) {
        addUniqueEntry(
          orthogroupFeatureCanonicalKey(...canonical.value),
          owner,
          entry
        );
      }
    });
  });

  index.recordsWithRenderedIds = recordsWithRenderedIds;
  return index;
};

const resolveOrthogroupFeatureMetadata = (
  index,
  feature,
  fallbackRecordIndex = null
) => {
  if (!(index instanceof Map)) return null;
  const record = recordIdentity(feature, fallbackRecordIndex);
  const featureIndex = sourceFeatureIndex(feature);
  const identity = featureIdentity(feature);
  const canonical = canonicalIdentity(feature);
  const hasRecordIdentity = Boolean(
    identity.source.value || identity.rendered.value
  ) || featureIndex.value !== null;
  if (
    !record.valid
    || !featureIndex.valid
    || !identity.valid
    || !canonical.valid
    || (hasRecordIdentity && record.value === null)
    || (!hasRecordIdentity && canonical.value === null)
  ) return null;
  const keys = identity.source.value
    ? [orthogroupFeatureIndexKey(record.value, identity.source.value)]
    : [];
  if (
    identity.rendered.value
    && index.recordsWithRenderedIds?.has(record.value)
  ) {
    keys.push(
      orthogroupRenderedFeatureIndexKey(record.value, identity.rendered.value)
    );
  }
  if (featureIndex.value !== null) {
    keys.push(orthogroupFeatureIndexSourceKey(record.value, featureIndex.value));
  }
  if (canonical.value) {
    keys.push(orthogroupFeatureCanonicalKey(...canonical.value));
  }
  const matches = new Set();
  for (const key of keys) {
    const entry = index.get(key);
    if (!entry) return null;
    matches.add(entry);
  }
  return keys.length > 0 && matches.size === 1
    ? matches.values().next().value
    : null;
};

export const enrichFeatureWithOrthogroup = (
  index,
  feature,
  fallbackRecordIndex = null
) => {
  const entry = resolveOrthogroupFeatureMetadata(index, feature, fallbackRecordIndex);
  const baseFeature = clearDerivedOrthogroupMetadata(feature);
  if (!entry) return baseFeature;
  return {
    ...baseFeature,
    proteinId: entry.proteinId,
    sourceProteinId: entry.sourceProteinId,
    orthogroupId: entry.orthogroupId,
    orthogroupMemberCount: entry.orthogroupMemberCount,
    orthogroupRecordCoverage: entry.orthogroupRecordCoverage,
    orthogroupRepresentative: entry.orthogroupRepresentative,
    orthogroupScope: entry.orthogroupScope,
    orthogroupSourceRecordIndex: entry.orthogroupSourceRecordIndex,
    orthogroupMember: entry.orthogroupMember
  };
};

export const enrichFeaturesWithOrthogroups = (features, index) => (
  (Array.isArray(features) ? features : []).map((feature) => (
    enrichFeatureWithOrthogroup(index, feature)
  ))
);
