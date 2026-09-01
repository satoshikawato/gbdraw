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

const applyDerivedOrthogroupMetadata = (descriptor, entry) => {
  const feature = descriptor?.feature;
  if (!feature || typeof feature !== 'object') return feature;
  ORTHOGROUP_FEATURE_FIELDS.forEach((field) => delete feature[field]);
  ['proteinId', 'sourceProteinId'].forEach((field) => {
    const original = descriptor.originalProteinFields[field];
    if (original.supplied) {
      feature[field] = original.value;
    } else {
      delete feature[field];
    }
  });
  if (!entry) return feature;
  Object.assign(feature, {
    proteinId: entry.proteinId,
    sourceProteinId: entry.sourceProteinId,
    orthogroupId: entry.orthogroupId,
    orthogroupMemberCount: entry.orthogroupMemberCount,
    orthogroupRecordCoverage: entry.orthogroupRecordCoverage,
    orthogroupRepresentative: entry.orthogroupRepresentative,
    orthogroupScope: entry.orthogroupScope,
    orthogroupSourceRecordIndex: entry.orthogroupSourceRecordIndex,
    orthogroupMember: entry.orthogroupMember
  });
  return feature;
};

const indexCandidateForMember = ({
  group,
  member,
  orthogroupId,
  memberCount,
  recordCoverage
}) => {
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
  ) return null;
  const entry = {
    orthogroupId,
    orthogroupMemberCount: memberCount,
    orthogroupRecordCoverage: recordCoverage,
    proteinId: normalizeText(member?.proteinId),
    sourceProteinId: normalizeText(member?.sourceProteinId),
    orthogroupRepresentative: Boolean(member?.representative),
    orthogroupScope: normalizeGroupMetadataScope(group?.scope),
    orthogroupSourceRecordIndex: normalizeRecordIndex(group?.source_record_index),
    orthogroupMember: member
  };
  const keys = [];
  if (identity.source.value) {
    keys.push(orthogroupFeatureIndexKey(record.value, identity.source.value));
  }
  if (identity.rendered.value) {
    keys.push(orthogroupRenderedFeatureIndexKey(record.value, identity.rendered.value));
  }
  if (featureIndex.value !== null) {
    keys.push(orthogroupFeatureIndexSourceKey(record.value, featureIndex.value));
  }
  if (canonical.value) {
    keys.push(orthogroupFeatureCanonicalKey(...canonical.value));
  }
  return {
    entry,
    keys,
    renderedRecordIndex: identity.rendered.value ? record.value : null
  };
};

/**
 * Incrementally project orthogroup metadata while the catalog owner is already
 * visiting groups and Features. This keeps the canonical admission path from
 * rebuilding an index and remapping the complete Feature population afterward.
 */
export const createOrthogroupFeatureProjection = () => {
  const index = new Map();
  const recordsWithRenderedIds = new Set();
  const candidatesByKey = new Map();
  const featureDescriptorsByKey = new Map();
  const featureDescriptorsByRecord = new Map();
  const groupCandidatesById = new Map();
  const duplicateGroupIds = new Set();
  const renderedRecordCounts = new Map();
  index.recordsWithRenderedIds = recordsWithRenderedIds;

  const refreshIndexKey = (key) => {
    const candidates = candidatesByKey.get(key);
    if (candidates?.size === 1) {
      index.set(key, candidates.values().next().value.entry);
    } else {
      index.delete(key);
    }
  };

  const refreshFeatureDescriptor = (descriptor) => {
    applyDerivedOrthogroupMetadata(
      descriptor,
      resolveOrthogroupFeatureMetadata(
        index,
        descriptor.feature,
        descriptor.fallbackRecordIndex
      )
    );
  };

  const refreshFeatureKeys = (keys) => {
    const descriptors = new Set();
    keys.forEach((key) => {
      featureDescriptorsByKey.get(key)?.forEach((descriptor) => descriptors.add(descriptor));
    });
    descriptors.forEach(refreshFeatureDescriptor);
  };

  const updateRenderedRecordCount = (recordIndex, delta, refreshDescriptors) => {
    if (recordIndex === null) return;
    const previous = renderedRecordCounts.get(recordIndex) || 0;
    const next = previous + delta;
    if (next > 0) {
      renderedRecordCounts.set(recordIndex, next);
      recordsWithRenderedIds.add(recordIndex);
    } else {
      renderedRecordCounts.delete(recordIndex);
      recordsWithRenderedIds.delete(recordIndex);
    }
    if (Boolean(previous) !== Boolean(next)) {
      featureDescriptorsByRecord.get(recordIndex)?.forEach((descriptor) => (
        refreshDescriptors.add(descriptor)
      ));
    }
  };

  const addCandidate = (candidate, refreshKeys, refreshDescriptors) => {
    candidate.keys.forEach((key) => {
      if (!candidatesByKey.has(key)) candidatesByKey.set(key, new Map());
      candidatesByKey.get(key).set(candidate.owner, candidate);
      refreshIndexKey(key);
      refreshKeys.add(key);
    });
    updateRenderedRecordCount(candidate.renderedRecordIndex, 1, refreshDescriptors);
  };

  const removeCandidate = (candidate, refreshKeys, refreshDescriptors) => {
    candidate.keys.forEach((key) => {
      const candidates = candidatesByKey.get(key);
      candidates?.delete(candidate.owner);
      if (candidates?.size === 0) candidatesByKey.delete(key);
      refreshIndexKey(key);
      refreshKeys.add(key);
    });
    updateRenderedRecordCount(candidate.renderedRecordIndex, -1, refreshDescriptors);
  };

  const addGroup = (group) => {
    const [{ id: orthogroupId } = {}] = uniqueOrthogroupEntries([group]);
    if (!orthogroupId || duplicateGroupIds.has(orthogroupId)) return;
    const refreshKeys = new Set();
    const refreshDescriptors = new Set();
    const previousCandidates = groupCandidatesById.get(orthogroupId);
    if (previousCandidates) {
      previousCandidates.forEach((candidate) => (
        removeCandidate(candidate, refreshKeys, refreshDescriptors)
      ));
      groupCandidatesById.delete(orthogroupId);
      duplicateGroupIds.add(orthogroupId);
      refreshFeatureKeys(refreshKeys);
      refreshDescriptors.forEach(refreshFeatureDescriptor);
      return;
    }

    const members = Array.isArray(group?.members) ? group.members : [];
    const memberCount = Number(group?.member_count || members.length || 0);
    const recordIndexes = new Set();
    const candidates = [];
    members.forEach((member) => {
      const recordIndex = normalizeRecordIndex(member?.recordIndex ?? member?.record_index);
      if (recordIndex !== null) recordIndexes.add(recordIndex);
      candidates.push(member);
    });
    const recordCoverage = Number(
      group?.record_coverage_count || recordIndexes.size || 0
    );
    const projectedCandidates = candidates
      .map((member) => indexCandidateForMember({
        group,
        member,
        orthogroupId,
        memberCount,
        recordCoverage
      }))
      .filter(Boolean)
      .map((candidate) => ({
        ...candidate,
        owner: {},
        recordIndex: normalizeRecordIndex(
          candidate.entry.orthogroupMember?.recordIndex
          ?? candidate.entry.orthogroupMember?.record_index
        )
    }));
    groupCandidatesById.set(orthogroupId, projectedCandidates);
    projectedCandidates.forEach((candidate) => (
      addCandidate(candidate, refreshKeys, refreshDescriptors)
    ));
    refreshFeatureKeys(refreshKeys);
    refreshDescriptors.forEach(refreshFeatureDescriptor);
  };

  const registerFeature = (feature, fallbackRecordIndex = null) => {
    const record = recordIdentity(feature, fallbackRecordIndex);
    const sourceIndex = sourceFeatureIndex(feature);
    const identity = featureIdentity(feature);
    const canonical = canonicalIdentity(feature);
    const keys = [];
    const descriptor = {
      feature,
      fallbackRecordIndex,
      originalProteinFields: Object.fromEntries(
        ['proteinId', 'sourceProteinId'].map((field) => [field, {
          supplied: Object.prototype.hasOwnProperty.call(feature, field),
          value: feature[field]
        }])
      )
    };
    if (record.value !== null) {
      if (!featureDescriptorsByRecord.has(record.value)) {
        featureDescriptorsByRecord.set(record.value, new Set());
      }
      featureDescriptorsByRecord.get(record.value).add(descriptor);
    }
    if (record.value !== null && identity.source.value) {
      keys.push(orthogroupFeatureIndexKey(record.value, identity.source.value));
    }
    if (record.value !== null && identity.rendered.value) {
      keys.push(orthogroupRenderedFeatureIndexKey(record.value, identity.rendered.value));
    }
    if (record.value !== null && sourceIndex.value !== null) {
      keys.push(orthogroupFeatureIndexSourceKey(record.value, sourceIndex.value));
    }
    if (canonical.value) keys.push(orthogroupFeatureCanonicalKey(...canonical.value));
    keys.forEach((key) => {
      if (!featureDescriptorsByKey.has(key)) featureDescriptorsByKey.set(key, new Set());
      featureDescriptorsByKey.get(key).add(descriptor);
    });
    refreshFeatureDescriptor(descriptor);
    return feature;
  };

  return Object.freeze({ index, addGroup, registerFeature });
};

export const buildOrthogroupFeatureIndex = (orthogroups) => {
  const projection = createOrthogroupFeatureProjection();
  (Array.isArray(orthogroups) ? orthogroups : []).forEach(projection.addGroup);
  return projection.index;
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
