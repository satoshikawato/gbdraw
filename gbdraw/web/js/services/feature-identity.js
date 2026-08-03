const normalizeText = (value) => String(value ?? '').trim();

export const RECORD_INDEX_KEYS = Object.freeze([
  'recordIndex', 'record_index', 'recordIdx', 'record_idx', 'fileIdx', 'file_idx'
]);
export const SOURCE_FEATURE_INDEX_KEYS = Object.freeze([
  'featureIndex', 'feature_index', 'sourceFeatureIndex', 'source_feature_index'
]);
export const STABLE_FEATURE_ID_KEYS = Object.freeze([
  'stableFeatureSvgId', 'stable_feature_svg_id', 'stableFeatureId',
  'stable_feature_id', 'stableSvgId', 'stable_svg_id', 'featureSvgId',
  'feature_svg_id'
]);
export const RENDERED_FEATURE_ID_KEYS = Object.freeze([
  'renderedFeatureSvgId', 'rendered_feature_svg_id', 'renderedSvgId',
  'rendered_svg_id'
]);
export const RECORD_KEY_KEYS = Object.freeze(['recordKey', 'record_key']);
export const BIOLOGICAL_FEATURE_ID_KEYS = Object.freeze([
  'biologicalFeatureId', 'biological_feature_id'
]);
export const ORTHOGROUP_ID_KEYS = Object.freeze([
  'id', 'orthogroupId', 'orthogroup_id'
]);

export const textAliasStatus = (source, keys) => {
  let supplied = false;
  let value = '';
  for (const key of keys) {
    const normalized = normalizeText(source?.[key]);
    if (!normalized) continue;
    supplied = true;
    if (value && value !== normalized) return { valid: false, supplied, value: '' };
    value = normalized;
  }
  return { valid: true, supplied, value };
};

export const nonnegativeIntegerAliasStatus = (source, keys) => {
  let supplied = false;
  let value = null;
  for (const key of keys) {
    const raw = source?.[key];
    if (raw === null || raw === undefined || raw === '') continue;
    supplied = true;
    const normalized = typeof raw === 'number'
      ? raw
      : (typeof raw === 'string' && /^\d+$/.test(raw.trim()) ? Number(raw.trim()) : Number.NaN);
    if (!Number.isSafeInteger(normalized) || normalized < 0) {
      return { valid: false, supplied, value: null };
    }
    if (value !== null && value !== normalized) {
      return { valid: false, supplied, value: null };
    }
    value = normalized;
  }
  return { valid: true, supplied, value };
};

export const orthogroupIdStatus = (source) => (
  textAliasStatus(source, ORTHOGROUP_ID_KEYS)
);

export const uniqueOrthogroupEntries = (orthogroups) => {
  const counts = new Map();
  const entries = [];
  (Array.isArray(orthogroups) ? orthogroups : []).forEach((group) => {
    const status = orthogroupIdStatus(group);
    if (!status.valid || !status.supplied) return;
    counts.set(status.value, (counts.get(status.value) || 0) + 1);
    entries.push({ group, id: status.value });
  });
  return entries.filter(({ id }) => counts.get(id) === 1);
};

export const featureIdentity = (source, { allowLegacySvgStable = false } = {}) => {
  const recordIndex = nonnegativeIntegerAliasStatus(source, RECORD_INDEX_KEYS);
  const sourceIndex = nonnegativeIntegerAliasStatus(source, SOURCE_FEATURE_INDEX_KEYS);
  const recordKey = textAliasStatus(source, RECORD_KEY_KEYS);
  const biologicalId = textAliasStatus(source, BIOLOGICAL_FEATURE_ID_KEYS);
  const renderedId = textAliasStatus(source, RENDERED_FEATURE_ID_KEYS);
  let stableId = textAliasStatus(source, STABLE_FEATURE_ID_KEYS);
  if (allowLegacySvgStable && !renderedId.supplied) {
    const legacySvgId = textAliasStatus(source, ['svg_id', 'svgId']);
    if (
      !legacySvgId.valid ||
      (stableId.supplied && legacySvgId.supplied && stableId.value !== legacySvgId.value)
    ) {
      stableId = { valid: false, supplied: true, value: '' };
    } else if (!stableId.supplied) {
      stableId = legacySvgId;
    }
  }
  const valid = [recordIndex, sourceIndex, recordKey, biologicalId, renderedId, stableId]
    .every((status) => status.valid)
    && recordKey.supplied === biologicalId.supplied;
  const scoped = recordIndex.supplied || recordKey.supplied;
  const identified = stableId.supplied || sourceIndex.supplied || biologicalId.supplied;
  return {
    valid,
    usable: valid && scoped && identified,
    recordIndex,
    sourceIndex,
    recordKey,
    biologicalId,
    renderedId,
    stableId
  };
};

export const renderedFeatureIdentity = (feature) => {
  const identity = featureIdentity(feature);
  const renderedId = textAliasStatus(feature, [
    ...RENDERED_FEATURE_ID_KEYS,
    'svg_id',
    'svgId'
  ]);
  return {
    ...identity,
    valid: identity.valid && renderedId.valid,
    usable: identity.usable && renderedId.valid && renderedId.supplied,
    renderedId
  };
};

export const identityMatches = (reference, candidate, { includeRendered = false } = {}) => {
  if (!reference?.usable || !candidate?.usable) return false;
  const fields = ['recordIndex', 'sourceIndex', 'recordKey', 'biologicalId', 'stableId'];
  if (includeRendered) fields.push('renderedId');
  return fields.every((field) => (
    !reference[field].supplied || (
      candidate[field].supplied && candidate[field].value === reference[field].value
    )
  ));
};

export const resolveUniqueOrthogroupMemberForFeature = (feature, members = []) => {
  const featureIdentityValue = renderedFeatureIdentity(feature);
  if (!featureIdentityValue.usable) return null;
  const matches = (Array.isArray(members) ? members : []).filter((member) => {
    const memberIdentityValue = featureIdentity(member);
    return (
      identityMatches(featureIdentityValue, memberIdentityValue, { includeRendered: true }) &&
      identityMatches(memberIdentityValue, featureIdentityValue, { includeRendered: true })
    );
  });
  return matches.length === 1 ? matches[0] : null;
};
