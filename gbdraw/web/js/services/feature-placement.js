// Requested placement wire validation shared by the codec and editable drafts.
export const canonicalFeaturePlacements = (overrides, mode) => {
  const rows = Array.isArray(overrides) ? overrides : Object.values(overrides);
  const identities = new Set();
  const sides = mode === 'circular' ? ['outward', 'inward'] : ['above', 'below'];
  return rows.map((row) => {
    if (!row || Object.keys(row).sort().join(',') !== 'biologicalFeatureId,placement,recordKey'
      || [row.recordKey, row.biologicalFeatureId].some((id) => typeof id !== 'string' || !id.trim() || id.includes('\0'))) {
      throw new Error('Feature placement requires an exact record and biological feature identity.');
    }
    const key = JSON.stringify([row.recordKey, row.biologicalFeatureId]);
    if (!Array.isArray(overrides) && overrides[key] !== row) {
      throw new Error('Feature placement draft keys must encode their exact identity as a JSON pair.');
    }
    if (identities.has(key)) throw new Error('Duplicate feature placement identity.');
    identities.add(key);
    const target = row.placement;
    if (!target || (target.kind === 'main'
      ? Object.keys(target).join(',') !== 'kind'
      : target.kind !== 'lane' || Object.keys(target).sort().join(',') !== 'kind,level,side'
        || !sides.includes(target.side) || target.level !== 1)) {
      throw new Error('Feature placement must be Main or a mode-compatible lane 1; Auto removes the override.');
    }
    return { recordKey: row.recordKey, biologicalFeatureId: row.biologicalFeatureId, placement: { ...target } };
  }).sort((a, b) => compareSourceIdentity(a.recordKey, b.recordKey)
    || compareSourceIdentity(a.biologicalFeatureId, b.biologicalFeatureId));
};

const compareSourceIdentity = (left, right) => {
  const a = Array.from(left, (value) => value.codePointAt(0));
  const b = Array.from(right, (value) => value.codePointAt(0));
  for (let index = 0; index < Math.min(a.length, b.length); index += 1) {
    if (a[index] !== b[index]) return a[index] - b[index];
  }
  return a.length - b.length;
};

