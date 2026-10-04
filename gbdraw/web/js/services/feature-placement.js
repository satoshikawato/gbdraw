import { getFeatureCaption } from '../app/feature-utils.js';
import { diagnosticError } from './error-normalization.js';

// Drafts keyed by original-source feature identity: Feature placement rows and
// per-feature edits (design Q4). Each draft map is
// {JSON.stringify([recordKey, biologicalFeatureId]): row}; this module owns the
// key, the row checks, and the projection onto a render request (R2).

const SIDES = Object.freeze({ circular: ['outward', 'inward'], linear: ['above', 'below'] });
// A malformed row can come only from a Session file; the editor writes exact rows.
const invalidRow = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });

const validIdentityText = (value) => typeof value === 'string' && Boolean(value.trim()) && !value.includes('\0');

export const featureIdentityKey = (recordKey, biologicalFeatureId) => (
  validIdentityText(recordKey) && validIdentityText(biologicalFeatureId)
    ? JSON.stringify([recordKey, biologicalFeatureId])
    : ''
);

// The identity key of a catalog feature (every rendered and biological
// feature of an admitted catalog carries both fields).
export const featureIdentityKeyOf = (feature) => featureIdentityKey(
  feature?.record_key ?? feature?.recordKey,
  feature?.biological_feature_id ?? feature?.biologicalFeatureId
);

const compareSourceIdentity = (left, right) => {
  const a = Array.from(left, (value) => value.codePointAt(0));
  const b = Array.from(right, (value) => value.codePointAt(0));
  for (let index = 0; index < Math.min(a.length, b.length); index += 1) {
    if (a[index] !== b[index]) return a[index] - b[index];
  }
  return a.length - b.length;
};

// Validates a draft map (or a request array) of identity rows with one row
// check, and returns the rows in code-point order of their identity.
const canonicalIdentityRows = (overrides, projectRow) => {
  const draft = !Array.isArray(overrides);
  const rows = draft ? Object.entries(overrides || {}) : overrides.map((row) => [null, row]);
  const identities = new Set();
  return rows.map(([key, row]) => {
    if (!row || typeof row !== 'object'
      || !validIdentityText(row.recordKey) || !validIdentityText(row.biologicalFeatureId)) {
      throw invalidRow();
    }
    const identity = featureIdentityKey(row.recordKey, row.biologicalFeatureId);
    if ((draft && key !== identity) || identities.has(identity)) throw invalidRow();
    identities.add(identity);
    return projectRow(row);
  }).sort((a, b) => compareSourceIdentity(a.recordKey, b.recordKey)
    || compareSourceIdentity(a.biologicalFeatureId, b.biologicalFeatureId));
};

// Record keys are mode-specific, so a request carries only the draft rows of its
// own records (an ALL record also owns its <recordKey>:<n> expansions); the other
// mode's rows stay in the draft for that mode (OV-08, R2).
export const recordKeyBelongsToRequest = (recordKey, records = []) => (
  typeof recordKey === 'string' && records.some((record) => (
    record.recordKey === recordKey || (record.cardinality === 'all'
      && recordKey.startsWith(`${record.recordKey}:`) && /^[1-9]\d*$/.test(recordKey.slice(record.recordKey.length + 1)))))
);

// Requested placement wire validation shared by the codec and editable drafts.
// Without a mode it validates the draft, which keeps each mode's rows (R2).
export const canonicalFeaturePlacements = (overrides, mode = null) => {
  const sides = mode ? SIDES[mode] || [] : [...SIDES.circular, ...SIDES.linear];
  return canonicalIdentityRows(overrides, (row) => {
    if (Object.keys(row).sort().join(',') !== 'biologicalFeatureId,placement,recordKey') throw invalidRow();
    const target = row.placement;
    if (!target || (target.kind === 'main'
      ? Object.keys(target).join(',') !== 'kind'
      : target.kind !== 'lane' || Object.keys(target).sort().join(',') !== 'kind,level,side'
        || !sides.includes(target.side) || target.level !== 1)) {
      throw invalidRow();
    }
    return { recordKey: row.recordKey, biologicalFeatureId: row.biologicalFeatureId, placement: { ...target } };
  });
};

// A draft row applies to a mode unless it holds the other mode's lane.
export const placementAppliesToMode = (row, mode) => row?.placement?.kind !== 'lane'
  || Boolean(SIDES[mode]?.includes(row.placement.side));

export const requestFeaturePlacements = (overrides, mode, records = []) => canonicalFeaturePlacements(
  Object.fromEntries(Object.entries(overrides || {})
    .filter(([, row]) => recordKeyBelongsToRequest(row?.recordKey, records) && placementAppliesToMode(row, mode))),
  mode
);

// Per-feature edits (request schema 9 `featureOverrides`). A null field keeps
// the rule-based result; the draft row adds the Web-only `labelSourceText`,
// the label's text before any edit, which a bulk label edit and the label
// projection read (B6). A draft row may hold only that source text.
export const FEATURE_OVERRIDE_EDIT_FIELDS = Object.freeze(['featureVisibility', 'labelVisibility', 'labelText']);
const REQUEST_OVERRIDE_FIELDS = 'biologicalFeatureId,featureVisibility,labelText,labelVisibility,recordKey';
const DRAFT_OVERRIDE_FIELDS = 'biologicalFeatureId,featureVisibility,labelSourceText,labelText,labelVisibility,recordKey';
const FEATURE_VISIBILITY_VALUES = new Set(['on', 'off', 'exclude_matching']);
const LABEL_VISIBILITY_VALUES = new Set(['on', 'off']);

// One line of label text, as the request accepts it; a blank text is null.
export const normalizeFeatureOverrideLabelText = (value) => {
  if (typeof value !== 'string') return null;
  const text = value.replace(/[\t\r\n\0]+/g, ' ');
  return text.trim() ? text : null;
};

const hasEdit = (row) => FEATURE_OVERRIDE_EDIT_FIELDS.some((field) => row[field] !== null);

const checkedOverrideRow = (row, { draft }) => {
  if (Object.keys(row).sort().join(',') !== (draft ? DRAFT_OVERRIDE_FIELDS : REQUEST_OVERRIDE_FIELDS)
    || !(row.featureVisibility === null || FEATURE_VISIBILITY_VALUES.has(row.featureVisibility))
    || !(row.labelVisibility === null || LABEL_VISIBILITY_VALUES.has(row.labelVisibility))
    || !(row.labelText === null || normalizeFeatureOverrideLabelText(row.labelText) === row.labelText)
    || (draft && !(row.labelSourceText === null
      || (typeof row.labelSourceText === 'string' && row.labelSourceText && !row.labelSourceText.includes('\0'))))
    || (!hasEdit(row) && !(draft && row.labelSourceText !== null))) {
    throw invalidRow();
  }
  return {
    recordKey: row.recordKey,
    biologicalFeatureId: row.biologicalFeatureId,
    featureVisibility: row.featureVisibility,
    labelVisibility: row.labelVisibility,
    labelText: row.labelText,
    ...(draft ? { labelSourceText: row.labelSourceText } : {})
  };
};

// Validates request `featureOverrides` rows (an array) or the editor draft map.
export const canonicalFeatureOverrides = (overrides) => canonicalIdentityRows(
  overrides, (row) => checkedOverrideRow(row, { draft: !Array.isArray(overrides) })
);

// The request rows of the current request's records: the draft's edits, and
// the label text that a bulk label edit gives a feature without its own
// (`bulkLabelText`, {identityKey: text}). Rows of other records stay in the
// draft (R2).
export const requestFeatureOverrides = (overrides, records = [], { bulkLabelText = {} } = {}) => {
  const rows = new Map();
  Object.entries(overrides || {}).forEach(([key, row]) => {
    if (recordKeyBelongsToRequest(row?.recordKey, records)) rows.set(key, row);
  });
  Object.entries(bulkLabelText || {}).forEach(([key, text]) => {
    const row = rows.get(key);
    const labelText = normalizeFeatureOverrideLabelText(text);
    if (row?.labelText || !labelText) return;
    let identity;
    try { identity = JSON.parse(key); } catch { return; }
    if (!Array.isArray(identity) || featureIdentityKey(identity[0], identity[1]) !== key
      || !recordKeyBelongsToRequest(identity[0], records)) return;
    rows.set(key, { ...(row || emptyOverrideRow(identity[0], identity[1])), labelText });
  });
  return canonicalFeatureOverrides(Object.fromEntries([...rows].filter(([, row]) => hasEdit(row))))
    .map(({ labelSourceText: _source, ...row }) => row);
};

const emptyOverrideRow = (recordKey, biologicalFeatureId) => ({
  recordKey,
  biologicalFeatureId,
  featureVisibility: null,
  labelVisibility: null,
  labelText: null,
  labelSourceText: null
});

// The draft value of one field for a feature, or null.
export const featureOverrideValue = (overrides, feature, field) => {
  const key = featureIdentityKeyOf(feature);
  return key ? overrides?.[key]?.[field] ?? null : null;
};

// Sets fields of one feature's draft row and deletes a row that no longer
// holds an edit or a label source text. Returns whether the draft changed.
export const updateFeatureOverride = (overrides, feature, patch) => {
  const recordKey = feature?.record_key ?? feature?.recordKey;
  const biologicalFeatureId = feature?.biological_feature_id ?? feature?.biologicalFeatureId;
  const key = featureIdentityKey(recordKey, biologicalFeatureId);
  if (!key || !overrides) return false;
  const next = { ...(overrides[key] || emptyOverrideRow(recordKey, biologicalFeatureId)) };
  Object.entries(patch || {}).forEach(([field, value]) => {
    if (!(field in next) || field === 'recordKey' || field === 'biologicalFeatureId') return;
    next[field] = value === undefined ? null : value;
  });
  const keep = hasEdit(next) || next.labelSourceText !== null;
  if (JSON.stringify(overrides[key] || null) === JSON.stringify(keep ? next : null)) return false;
  if (keep) overrides[key] = next;
  else delete overrides[key];
  return true;
};

const FEATURE_IDENTITY_NOTICE_STATUSES = new Set(['crop_excluded', 'absent', 'unresolved']);
const FEATURE_IDENTITY_NOTICE_KINDS = new Set(['placement', 'feature_visibility', 'label_visibility', 'label_text']);
// Python's featureIdentityNotices metadata, as Generate admitted it (design Q4 3.4).
export const validateFeatureIdentityNotices = (notices, results) => {
  if (notices === undefined) return;
  if (!Array.isArray(notices)) throw new Error('runMetadata.featureIdentityNotices must be an array.');
  const resultCount = Array.isArray(results) ? results.length : 0;
  notices.forEach((notice) => {
    if (!notice || typeof notice !== 'object' || Array.isArray(notice)
      || Object.keys(notice).sort().join(',') !== 'biologicalFeatureId,kinds,recordKey,resultIndex,status'
      || typeof notice.recordKey !== 'string' || !notice.recordKey
      || typeof notice.biologicalFeatureId !== 'string' || !notice.biologicalFeatureId
      || !FEATURE_IDENTITY_NOTICE_STATUSES.has(notice.status)
      || !Array.isArray(notice.kinds) || notice.kinds.length === 0
      || notice.kinds.some((kind) => !FEATURE_IDENTITY_NOTICE_KINDS.has(kind))
      || !Number.isInteger(notice.resultIndex) || notice.resultIndex < 0 || notice.resultIndex >= resultCount) {
      throw new Error('runMetadata.featureIdentityNotices contains an invalid notice.');
    }
  });
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
