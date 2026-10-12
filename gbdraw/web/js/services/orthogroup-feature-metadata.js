// @ts-check
import { normalizeGroupMetadataScope } from './losat-normalization.js';
import {
  RECORD_INDEX_KEYS,
  RENDERED_FEATURE_ID_KEYS,
  SOURCE_FEATURE_INDEX_KEYS,
  nonnegativeIntegerAliasStatus,
  textAliasStatus,
  uniqueOrthogroupEntries
} from './feature-identity.js';

/** @import { AliasSource } from './feature-identity.js' */

/** @typedef {AliasSource} OrthogroupFeature A Feature record; the alias fields are read, the orthogroup fields written. */
/** @typedef {AliasSource} OrthogroupMember A member record of a saved orthogroup. */
/**
 * @typedef {object} OrthogroupGroup
 * @property {unknown} [id]
 * @property {unknown} [members]
 * @property {unknown} [member_count]
 * @property {unknown} [record_coverage_count]
 * @property {unknown} [scope]
 * @property {unknown} [source_record_index]
 */
/**
 * The orthogroup fields one Feature takes from its member.
 * @typedef {object} OrthogroupMetadataEntry
 * @property {string} orthogroupId
 * @property {number} orthogroupMemberCount
 * @property {number} orthogroupRecordCoverage
 * @property {string} proteinId
 * @property {string} sourceProteinId
 * @property {boolean} orthogroupRepresentative
 * @property {string} orthogroupScope
 * @property {number | null} orthogroupSourceRecordIndex
 * @property {OrthogroupMember} orthogroupMember
 */
/**
 * @typedef {object} IndexCandidate
 * @property {OrthogroupMetadataEntry} entry
 * @property {string[]} keys The index keys the member answers to.
 * @property {number | null} renderedRecordIndex The record whose rendered IDs the member names, else null.
 */
/** @typedef {IndexCandidate & { owner: Record<string, never>, recordIndex: number | null }} ProjectedCandidate `owner` is the identity of the candidate in a key's candidate map. */
/** @typedef {{ supplied: boolean, value: unknown }} OriginalProteinField */
/**
 * A registered Feature and what restores it before a refresh.
 * @typedef {object} FeatureDescriptor
 * @property {OrthogroupFeature} feature
 * @property {number | null} fallbackRecordIndex
 * @property {Record<string, OriginalProteinField>} originalProteinFields
 */
/** @typedef {{ name?: string, description?: string }} OrthogroupOverrideText */

/** @param {unknown} value */
const normalizeText = (value) => String(value || '').trim();

/**
 * @param {unknown} value
 * @returns {number | null}
 */
const normalizeRecordIndex = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const text = String(value).trim();
  if (!/^\d+$/.test(text) || typeof value === 'boolean') return null;
  const normalized = Number(text);
  return Number.isSafeInteger(normalized) ? normalized : null;
};

/**
 * @param {number | null} recordIndex
 * @param {string} featureSvgId
 */
const orthogroupFeatureIndexKey = (recordIndex, featureSvgId) => (
  `record-id:${Number(recordIndex)}:${normalizeText(featureSvgId)}`
);

/**
 * @param {number | null} recordIndex
 * @param {string} featureSvgId
 */
const orthogroupRenderedFeatureIndexKey = (recordIndex, featureSvgId) => (
  `record-rendered-id:${Number(recordIndex)}:${normalizeText(featureSvgId)}`
);

/**
 * @param {number | null} recordIndex
 * @param {number} featureIndex
 */
const orthogroupFeatureIndexSourceKey = (recordIndex, featureIndex) => (
  `record-feature-index:${Number(recordIndex)}:${Number(featureIndex)}`
);

/**
 * @param {string} recordKey
 * @param {string} biologicalFeatureId
 */
const orthogroupFeatureCanonicalKey = (recordKey, biologicalFeatureId) => (
  `canonical:${normalizeText(recordKey)}\u0000${normalizeText(biologicalFeatureId)}`
);

/**
 * @param {AliasSource | null | undefined} payload
 * @param {readonly string[]} keys
 */
const consistentTextAlias = (payload, keys) => {
  const status = textAliasStatus(payload, keys);
  return {
    valid: status.valid,
    value: status.value
  };
};

/**
 * @param {AliasSource | null | undefined} payload
 * @param {readonly string[]} keys
 * @param {number | string | null} [fallback]
 */
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

/** @param {AliasSource | null | undefined} payload */
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

/** @param {AliasSource | null | undefined} payload */
const renderedIdentity = (payload) => consistentTextAlias(payload, [
  ...RENDERED_FEATURE_ID_KEYS
]);

/** @param {AliasSource | null | undefined} payload */
const memberIdentity = (payload) => {
  const source = sourceIdentity(payload);
  const rendered = renderedIdentity(payload);
  return { valid: source.valid && rendered.valid, source, rendered };
};

/** @param {AliasSource | null | undefined} payload */
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

/**
 * @param {AliasSource | null | undefined} payload
 * @returns {{ valid: boolean, value: [string, string] | null }}
 */
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

/** @param {AliasSource | null | undefined} payload */
const sourceFeatureIndex = (payload) => consistentIntegerAlias(
  payload,
  SOURCE_FEATURE_INDEX_KEYS
);

/**
 * @param {AliasSource | null | undefined} payload
 * @param {number | string | null} [fallback]
 */
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

/** @param {OrthogroupFeature | null | undefined} feature */
const clearDerivedOrthogroupMetadata = (feature) => {
  if (!feature || !ORTHOGROUP_FEATURE_FIELDS.some((field) => (
    Object.prototype.hasOwnProperty.call(feature, field)
  ))) return feature;
  const cleared = { ...feature };
  ORTHOGROUP_FEATURE_FIELDS.forEach((field) => delete cleared[field]);
  return cleared;
};

/**
 * @param {FeatureDescriptor} descriptor
 * @param {OrthogroupMetadataEntry | null} entry
 */
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

/**
 * @param {{
 *   group: OrthogroupGroup,
 *   member: OrthogroupMember,
 *   orthogroupId: string,
 *   memberCount: number,
 *   recordCoverage: number
 * }} input
 * @returns {IndexCandidate | null}
 */
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
  /** @type {string[]} */
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

/** @typedef {Map<string, OrthogroupMetadataEntry> & { recordsWithRenderedIds?: Set<number> }} OrthogroupFeatureIndex */

/**
 * Incrementally project orthogroup metadata while the catalog owner is already
 * visiting groups and Features. This keeps the canonical admission path from
 * rebuilding an index and remapping the complete Feature population afterward.
 */
export const createOrthogroupFeatureProjection = () => {
  /** @type {OrthogroupFeatureIndex} */
  const index = new Map();
  /** @type {Set<number>} */
  const recordsWithRenderedIds = new Set();
  /** @type {Map<string, Map<Record<string, never>, ProjectedCandidate>>} */
  const candidatesByKey = new Map();
  /** @type {Map<string, Set<FeatureDescriptor>>} */
  const featureDescriptorsByKey = new Map();
  /** @type {Map<number, Set<FeatureDescriptor>>} */
  const featureDescriptorsByRecord = new Map();
  /** @type {Map<string, ProjectedCandidate[]>} */
  const groupCandidatesById = new Map();
  /** @type {Set<string>} */
  const duplicateGroupIds = new Set();
  /** @type {Map<number, number>} */
  const renderedRecordCounts = new Map();
  index.recordsWithRenderedIds = recordsWithRenderedIds;

  /** @param {string} key */
  const refreshIndexKey = (key) => {
    const candidates = candidatesByKey.get(key);
    if (candidates?.size === 1) {
      // The one candidate of a key whose map holds one.
      index.set(key, /** @type {ProjectedCandidate} */ (candidates.values().next().value).entry);
    } else {
      index.delete(key);
    }
  };

  /** @param {FeatureDescriptor} descriptor */
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

  /** @param {Set<string>} keys */
  const refreshFeatureKeys = (keys) => {
    /** @type {Set<FeatureDescriptor>} */
    const descriptors = new Set();
    keys.forEach((key) => {
      featureDescriptorsByKey.get(key)?.forEach((descriptor) => descriptors.add(descriptor));
    });
    descriptors.forEach(refreshFeatureDescriptor);
  };

  /**
   * @param {number | null} recordIndex
   * @param {number} delta
   * @param {Set<FeatureDescriptor>} refreshDescriptors
   */
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

  /**
   * @param {ProjectedCandidate} candidate
   * @param {Set<string>} refreshKeys
   * @param {Set<FeatureDescriptor>} refreshDescriptors
   */
  const addCandidate = (candidate, refreshKeys, refreshDescriptors) => {
    candidate.keys.forEach((key) => {
      if (!candidatesByKey.has(key)) candidatesByKey.set(key, new Map());
      /** @type {Map<Record<string, never>, ProjectedCandidate>} */ (candidatesByKey.get(key)).set(candidate.owner, candidate); // set above
      refreshIndexKey(key);
      refreshKeys.add(key);
    });
    updateRenderedRecordCount(candidate.renderedRecordIndex, 1, refreshDescriptors);
  };

  /**
   * @param {ProjectedCandidate} candidate
   * @param {Set<string>} refreshKeys
   * @param {Set<FeatureDescriptor>} refreshDescriptors
   */
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

  /** @param {OrthogroupGroup} group */
  const addGroup = (group) => {
    const [{ id: orthogroupId } = {}] = uniqueOrthogroupEntries([group]);
    if (!orthogroupId || duplicateGroupIds.has(orthogroupId)) return;
    /** @type {Set<string>} */
    const refreshKeys = new Set();
    /** @type {Set<FeatureDescriptor>} */
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

    /** @type {OrthogroupMember[]} */
    const members = Array.isArray(group?.members) ? group.members : [];
    const memberCount = Number(group?.member_count || members.length || 0);
    /** @type {Set<number>} */
    const recordIndexes = new Set();
    /** @type {OrthogroupMember[]} */
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
      .filter((candidate) => candidate !== null)
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

  /**
   * @param {OrthogroupFeature} feature
   * @param {number | null} [fallbackRecordIndex]
   */
  const registerFeature = (feature, fallbackRecordIndex = null) => {
    const record = recordIdentity(feature, fallbackRecordIndex);
    const sourceIndex = sourceFeatureIndex(feature);
    const identity = featureIdentity(feature);
    const canonical = canonicalIdentity(feature);
    /** @type {string[]} */
    const keys = [];
    /** @type {FeatureDescriptor} */
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
      /** @type {Set<FeatureDescriptor>} */ (featureDescriptorsByRecord.get(record.value)).add(descriptor); // set above
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
      /** @type {Set<FeatureDescriptor>} */ (featureDescriptorsByKey.get(key)).add(descriptor); // set above
    });
    refreshFeatureDescriptor(descriptor);
    return feature;
  };

  return Object.freeze({ index, addGroup, registerFeature });
};

/** @param {readonly OrthogroupGroup[] | null | undefined} orthogroups */
export const buildOrthogroupFeatureIndex = (orthogroups) => {
  const projection = createOrthogroupFeatureProjection();
  (Array.isArray(orthogroups) ? orthogroups : []).forEach(projection.addGroup);
  return projection.index;
};

/**
 * @param {OrthogroupFeatureIndex} index
 * @param {OrthogroupFeature} feature
 * @param {number | null} [fallbackRecordIndex]
 * @returns {OrthogroupMetadataEntry | null}
 */
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
    && index.recordsWithRenderedIds?.has(/** @type {number} */ (record.value)) // a rendered ID names a record above
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
  /** @type {Set<OrthogroupMetadataEntry>} */
  const matches = new Set();
  for (const key of keys) {
    const entry = index.get(key);
    if (!entry) return null;
    matches.add(entry);
  }
  return keys.length > 0 && matches.size === 1
    ? /** @type {OrthogroupMetadataEntry} */ (matches.values().next().value) // one match
    : null;
};

/**
 * @param {OrthogroupFeatureIndex} index
 * @param {OrthogroupFeature} feature
 * @param {number | null} [fallbackRecordIndex]
 */
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

/**
 * @param {readonly OrthogroupFeature[] | null | undefined} features
 * @param {OrthogroupFeatureIndex} index
 */
export const enrichFeaturesWithOrthogroups = (features, index) => (
  (Array.isArray(features) ? features : []).map((feature) => (
    enrichFeatureWithOrthogroup(index, feature)
  ))
);

/** @param {OrthogroupMember | null | undefined} member */
const memberSetEntry = (member) => normalizeText(member?.proteinId)
  || [member?.recordId, member?.featureSvgId, member?.start, member?.end].map(normalizeText).join('|');

// The member set that a user name belongs to (D-21, PD-OI-075).
/** @param {OrthogroupGroup | null | undefined} group */
export const orthogroupMemberSetKey = (group) => JSON.stringify(
  [...new Set((Array.isArray(group?.members) ? group.members : []).map(memberSetEntry).filter(Boolean))].sort()
);

/**
 * @param {unknown} value
 * @returns {Record<string, OrthogroupOverrideText>}
 */
export const normalizeOrthogroupDormantOverrides = (value) => {
  /** @type {Record<string, OrthogroupOverrideText>} */
  const result = {};
  if (!value || typeof value !== 'object' || Array.isArray(value)) return result;
  // A saved override table maps a member-set key to a `{ name, description }` record.
  Object.entries(/** @type {Record<string, { name?: unknown, description?: unknown } | null | undefined>} */ (value)).forEach(([key, entry]) => {
    const name = normalizeText(entry?.name);
    const description = normalizeText(entry?.description);
    if (!key || key === '[]' || (!name && !description)) return;
    result[key] = { ...(name ? { name } : {}), ...(description ? { description } : {}) };
  });
  return result;
};

// D-21 (PD-OI-075): a group name and description follow the exact member set
// across a regrouping. Names whose members no longer form one group stay
// dormant and return when that group forms again; nothing moves by ID alone.
/**
 * @param {{
 *   previousGroups?: readonly OrthogroupGroup[] | null,
 *   candidateGroups?: readonly OrthogroupGroup[] | null,
 *   names?: Record<string, unknown>,
 *   descriptions?: Record<string, unknown>,
 *   dormant?: unknown
 * }} input
 */
export const rekeyOrthogroupOverrides = ({
  previousGroups, candidateGroups, names = {}, descriptions = {}, dormant = {}
}) => {
  const pool = new Map(Object.entries(normalizeOrthogroupDormantOverrides(dormant)));
  (Array.isArray(previousGroups) ? previousGroups : []).forEach((group) => {
    const id = normalizeText(group?.id);
    const name = normalizeText(names[id]);
    const description = normalizeText(descriptions[id]);
    const key = orthogroupMemberSetKey(group);
    if (!id || key === '[]' || (!name && !description)) return;
    pool.set(key, { ...(name ? { name } : {}), ...(description ? { description } : {}) });
  });
  /** @type {{ names: Record<string, string>, descriptions: Record<string, string> }} */
  const next = { names: {}, descriptions: {} };
  (Array.isArray(candidateGroups) ? candidateGroups : []).forEach((group) => {
    const id = normalizeText(group?.id);
    const key = orthogroupMemberSetKey(group);
    const entry = id && pool.get(key);
    if (!entry) return;
    if (entry.name) next.names[id] = entry.name;
    if (entry.description) next.descriptions[id] = entry.description;
    pool.delete(key);
  });
  return { ...next, dormant: Object.fromEntries(pool) };
};
