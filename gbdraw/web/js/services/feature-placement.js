import { getFeatureCaption } from '../app/feature-utils.js';
import { diagnosticError } from './error-normalization.js';

const SIDES = Object.freeze({ circular: ['outward', 'inward'], linear: ['above', 'below'] });
// A malformed row can come only from a Session file; the editor writes exact rows.
const invalidRow = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });

// Requested placement wire validation shared by the codec and editable drafts.
// Without a mode it validates the draft, which keeps each mode's rows (R2).
export const canonicalFeaturePlacements = (overrides, mode = null) => {
  const rows = Array.isArray(overrides) ? overrides : Object.values(overrides);
  const identities = new Set();
  const sides = mode ? SIDES[mode] || [] : [...SIDES.circular, ...SIDES.linear];
  return rows.map((row) => {
    if (!row || Object.keys(row).sort().join(',') !== 'biologicalFeatureId,placement,recordKey'
      || [row.recordKey, row.biologicalFeatureId].some((id) => typeof id !== 'string' || !id.trim() || id.includes('\0'))) {
      throw invalidRow();
    }
    const key = JSON.stringify([row.recordKey, row.biologicalFeatureId]);
    if (!Array.isArray(overrides) && overrides[key] !== row) throw invalidRow();
    if (identities.has(key)) throw invalidRow();
    identities.add(key);
    const target = row.placement;
    if (!target || (target.kind === 'main'
      ? Object.keys(target).join(',') !== 'kind'
      : target.kind !== 'lane' || Object.keys(target).sort().join(',') !== 'kind,level,side'
        || !sides.includes(target.side) || target.level !== 1)) {
      throw invalidRow();
    }
    return { recordKey: row.recordKey, biologicalFeatureId: row.biologicalFeatureId, placement: { ...target } };
  }).sort((a, b) => compareSourceIdentity(a.recordKey, b.recordKey)
    || compareSourceIdentity(a.biologicalFeatureId, b.biologicalFeatureId));
};

// A draft row applies to a mode unless it holds the other mode's lane.
export const placementAppliesToMode = (row, mode) => row?.placement?.kind !== 'lane'
  || Boolean(SIDES[mode]?.includes(row.placement.side));

// Record keys are mode-specific, so a request carries only the draft rows of its
// own records (an ALL record also owns its <recordKey>:<n> expansions); the other
// mode's rows stay in the draft for that mode (OV-08, R2).
export const requestFeaturePlacements = (overrides, mode, records = []) => {
  const belongs = (recordKey) => typeof recordKey === 'string' && records.some((record) => (
    record.recordKey === recordKey || (record.cardinality === 'all'
      && recordKey.startsWith(`${record.recordKey}:`) && /^[1-9]\d*$/.test(recordKey.slice(record.recordKey.length + 1)))));
  return canonicalFeaturePlacements(Object.fromEntries(Object.entries(overrides || {})
    .filter(([, row]) => belongs(row?.recordKey) && placementAppliesToMode(row, mode))), mode);
};

// Python reports an unsupported lane by its row in the request's canonical
// featurePlacements (R6); the Web adds the caption of that row's feature.
export const nameFeaturePlacementFailure = (error, request, features = []) => {
  const index = error?.context?.placementIndex;
  if (error?.code !== 'FEATURE_PLACEMENT' || !Number.isSafeInteger(index)) return error;
  let row;
  try { row = canonicalFeaturePlacements(request?.diagramOptions?.featurePlacements || [], request?.mode)[index]; }
  catch { return error; }
  const feature = row && features.find((item) => item?.record_key === row.recordKey
    && item?.biological_feature_id === row.biologicalFeatureId);
  return feature ? { ...error, context: { ...error.context, featureCaption: getFeatureCaption(feature) } } : error;
};

const compareSourceIdentity = (left, right) => {
  const a = Array.from(left, (value) => value.codePointAt(0));
  const b = Array.from(right, (value) => value.codePointAt(0));
  for (let index = 0; index < Math.min(a.length, b.length); index += 1) {
    if (a[index] !== b[index]) return a[index] - b[index];
  }
  return a.length - b.length;
};
