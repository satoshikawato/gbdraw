// @ts-check
import { getFeatureCaption } from './feature-utils.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { cloneJsonData } from './json-clone.js';

// Drafts keyed by original-source feature identity: Feature placement rows and
// per-feature edits (design Q4). A request chooses its record keys, so both
// modes can use the same key for the same feature (`record-1` in a Gallery or
// Python Session); each draft row therefore names its mode (`scope`, as record
// display drafts do). Each draft map is
// {JSON.stringify([scope, recordKey, biologicalFeatureId]): row}; this module
// owns the key, the row checks, and the projection onto a render request (R2).

/**
 * @typedef {'circular' | 'linear'} FeatureMode
 */

/**
 * The identity a draft key names.
 * @typedef {object} FeatureIdentity
 * @property {string} scope The mode of the request that chose the record key.
 * @property {string} recordKey
 * @property {string} biologicalFeatureId
 */

/**
 * Where a feature is drawn: the main track, or one lane beside it.
 * @typedef {{ kind: 'main' } | { kind: 'lane', side: string, level: number }} FeaturePlacementTarget
 */

/**
 * A request row of `diagramOptions.featurePlacements`.
 * @typedef {object} FeaturePlacementRow
 * @property {string} recordKey
 * @property {string} biologicalFeatureId
 * @property {FeaturePlacementTarget} placement
 */

/**
 * A draft row: a request row that names its mode.
 * @typedef {FeaturePlacementRow & { scope: string }} FeaturePlacementDraftRow
 */

/**
 * A request row of `featureOverrides`. A null field keeps the rule-based result.
 * @typedef {object} FeatureOverrideRow
 * @property {string} recordKey
 * @property {string} biologicalFeatureId
 * @property {string | null} featureVisibility `on`, `off`, `exclude_matching`, or null.
 * @property {string | null} labelVisibility `on`, `off`, or null.
 * @property {string | null} labelText
 */

/**
 * A draft row: a request row that names its mode and keeps the Web-only
 * `labelSourceText`, the label's text before any edit.
 * @typedef {FeatureOverrideRow & { scope: string, labelSourceText: string | null }} FeatureOverrideDraftRow
 */

/** @typedef {Record<string, FeaturePlacementDraftRow>} FeaturePlacementDraft */
/** @typedef {Record<string, FeatureOverrideDraftRow>} FeatureOverrideDraft */

/**
 * The records of a request, as a draft row's record key is matched against them.
 * @typedef {{ recordKey: string, cardinality?: string }} FeatureRequestRecord
 */

const SIDES = Object.freeze({ circular: ['outward', 'inward'], linear: ['above', 'below'] });
// A malformed row can come only from a Session file; the editor writes exact rows.
const invalidRow = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });

const validIdentityText = (value) => typeof value === 'string' && Boolean(value.trim()) && !value.includes('\0');

/**
 * @param {string} scope
 * @param {string} recordKey
 * @param {string} biologicalFeatureId
 * @returns {string} The draft key, or '' for an invalid identity.
 */
export const featureIdentityKey = (scope, recordKey, biologicalFeatureId) => (
  Object.hasOwn(SIDES, scope) && validIdentityText(recordKey) && validIdentityText(biologicalFeatureId)
    ? JSON.stringify([scope, recordKey, biologicalFeatureId])
    : ''
);

// The identity key of a draft row or a catalog feature (an admitted catalog
// gives each rendered and biological feature its Result's mode as `scope`).
/**
 * @param {Record<string, any> | null | undefined} feature
 * @returns {string}
 */
export const featureIdentityKeyOf = (feature) => featureIdentityKey(
  feature?.scope,
  feature?.record_key ?? feature?.recordKey,
  feature?.biological_feature_id ?? feature?.biologicalFeatureId
);

// The identity a draft key names ({scope, recordKey, biologicalFeatureId}), or null.
/**
 * @param {string} key
 * @returns {FeatureIdentity | null}
 */
export const parseFeatureIdentityKey = (key) => {
  let parts;
  try { parts = JSON.parse(key); } catch { return null; }
  if (!Array.isArray(parts) || featureIdentityKey(.../** @type {[string, string, string]} */ (parts)) !== key) return null;
  const [scope, recordKey, biologicalFeatureId] = parts;
  return { scope, recordKey, biologicalFeatureId };
};

// A draft map of rows that name their mode.
/**
 * @param {Array<Record<string, any>>} rows
 * @returns {Record<string, any>}
 */
export const featureDraftMap = (rows) => Object.fromEntries(rows.map((row) => [featureIdentityKeyOf(row), row]));

// The draft rows of a request's rows, which are of the request's mode.
/**
 * @param {Array<Record<string, any>>} rows Request rows.
 * @param {string} mode
 * @returns {Record<string, any>}
 */
export const draftRowsOfRequest = (rows, mode) => featureDraftMap(rows.map((row) => ({ scope: mode, ...row })));

const compareSourceIdentity = (left, right) => {
  const a = Array.from(left, (value) => value.codePointAt(0));
  const b = Array.from(right, (value) => value.codePointAt(0));
  for (let index = 0; index < Math.min(a.length, b.length); index += 1) {
    if (a[index] !== b[index]) return a[index] - b[index];
  }
  return a.length - b.length;
};

// Validates a draft map (or a request array) of identity rows with one row
// check, and returns the rows in code-point order of their identity. A draft
// row names its mode; a request row is of the request's mode.
const canonicalIdentityRows = (overrides, projectRow) => {
  const draft = !Array.isArray(overrides);
  const rows = draft ? Object.entries(overrides || {}) : overrides.map((row) => [null, row]);
  const identities = new Set();
  return rows.map(([key, row]) => {
    if (!row || typeof row !== 'object'
      || !validIdentityText(row.recordKey) || !validIdentityText(row.biologicalFeatureId)) {
      throw invalidRow();
    }
    const identity = draft ? featureIdentityKeyOf(row) : JSON.stringify([row.recordKey, row.biologicalFeatureId]);
    if (!identity || (draft && key !== identity) || identities.has(identity)) throw invalidRow();
    identities.add(identity);
    return projectRow(row, draft);
  }).sort((a, b) => compareSourceIdentity(a.scope || '', b.scope || '')
    || compareSourceIdentity(a.recordKey, b.recordKey)
    || compareSourceIdentity(a.biologicalFeatureId, b.biologicalFeatureId));
};

// Record keys are chosen per request, so a request carries only the draft rows
// of its own mode and records (an ALL record also owns its <recordKey>:<n>
// expansions); the other mode's rows stay in the draft for that mode (OV-08, R2).
/**
 * @param {unknown} recordKey
 * @param {FeatureRequestRecord[]} [records]
 * @returns {boolean}
 */
export const recordKeyBelongsToRequest = (recordKey, records = []) => (
  typeof recordKey === 'string' && records.some((record) => (
    record.recordKey === recordKey || (record.cardinality === 'all'
      && recordKey.startsWith(`${record.recordKey}:`) && /^[1-9]\d*$/.test(recordKey.slice(record.recordKey.length + 1)))))
);

/**
 * @param {Record<string, any> | null | undefined} row
 * @param {string} mode
 * @param {FeatureRequestRecord[]} [records]
 * @returns {boolean}
 */
export const rowBelongsToRequest = (row, mode, records = []) => (
  row?.scope === mode && recordKeyBelongsToRequest(row?.recordKey, records)
);

// Placement validation shared by the codec and editable drafts: a request row
// takes the sides of the request's mode, a draft row those of its own mode.
/**
 * @param {FeaturePlacementDraft | FeaturePlacementRow[] | null | undefined} overrides
 *   An editor draft map, or request rows.
 * @param {string | null} [mode] The request's mode, for request rows.
 * @returns {Array<FeaturePlacementRow | FeaturePlacementDraftRow>}
 */
export const canonicalFeaturePlacements = (overrides, mode = null) => canonicalIdentityRows(overrides, (row, draft) => {
  if (Object.keys(row).sort().join(',')
    !== `biologicalFeatureId,placement,recordKey${draft ? ',scope' : ''}`) throw invalidRow();
  const sides = SIDES[draft ? row.scope : mode] || [];
  const target = row.placement;
  if (!target || (target.kind === 'main'
    ? Object.keys(target).join(',') !== 'kind'
    : target.kind !== 'lane' || Object.keys(target).sort().join(',') !== 'kind,level,side'
      || !sides.includes(target.side) || target.level !== 1)) {
    throw invalidRow();
  }
  return {
    ...(draft ? { scope: row.scope } : {}),
    recordKey: row.recordKey,
    biologicalFeatureId: row.biologicalFeatureId,
    placement: { ...target }
  };
});

/**
 * @template {Record<string, any>} T
 * @param {T} row
 * @returns {Omit<T, 'scope'>}
 */
const requestRow = ({ scope: _scope, ...row }) => row;

/**
 * @param {FeaturePlacementDraft | null | undefined} overrides
 * @param {string} mode
 * @param {FeatureRequestRecord[]} [records]
 * @returns {FeaturePlacementRow[]}
 */
export const requestFeaturePlacements = (overrides, mode, records = []) => canonicalFeaturePlacements(
  Object.fromEntries(Object.entries(overrides || {}).filter(([, row]) => rowBelongsToRequest(row, mode, records)))
).map(requestRow);

// Replaces the Feature placement draft with copies of a checkpoint's draft rows:
// Undo, Redo, and the rollback of a Generate restore the placements it removed
// (Q3 = A).
/**
 * @param {FeaturePlacementDraft} overrides The draft, replaced in place.
 * @param {FeaturePlacementDraft | null | undefined} placements
 * @returns {void}
 */
export const restorePlacements = (overrides, placements) => {
  Object.keys(overrides).forEach((key) => delete overrides[key]);
  Object.assign(overrides, cloneJsonData(placements) || {});
};

// Per-feature edits (request schema 9 `featureOverrides`). A null field keeps
// the rule-based result; the draft row adds the Web-only `labelSourceText`,
// the label's text before any edit, which a bulk label edit and the label
// projection read (B6). A draft row may hold only that source text.
export const FEATURE_OVERRIDE_EDIT_FIELDS = Object.freeze(['featureVisibility', 'labelVisibility', 'labelText']);
const REQUEST_OVERRIDE_FIELDS = 'biologicalFeatureId,featureVisibility,labelText,labelVisibility,recordKey';
const DRAFT_OVERRIDE_FIELDS = 'biologicalFeatureId,featureVisibility,labelSourceText,labelText,labelVisibility,recordKey,scope';
const FEATURE_VISIBILITY_VALUES = new Set(['on', 'off', 'exclude_matching']);
const LABEL_VISIBILITY_VALUES = new Set(['on', 'off']);

// One line of label text, as the request accepts it; a blank text is null.
/**
 * @param {unknown} value
 * @returns {string | null}
 */
export const normalizeFeatureOverrideLabelText = (value) => {
  if (typeof value !== 'string') return null;
  const text = value.replace(/[\t\r\n\0]+/g, ' ');
  return text.trim() ? text : null;
};

const hasEdit = (row) => FEATURE_OVERRIDE_EDIT_FIELDS.some((field) => row[field] !== null);

const checkedOverrideRow = (row, draft) => {
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
    ...(draft ? { scope: row.scope } : {}),
    recordKey: row.recordKey,
    biologicalFeatureId: row.biologicalFeatureId,
    featureVisibility: row.featureVisibility,
    labelVisibility: row.labelVisibility,
    labelText: row.labelText,
    ...(draft ? { labelSourceText: row.labelSourceText } : {})
  };
};

// Validates request `featureOverrides` rows (an array) or the editor draft map.
/**
 * @param {FeatureOverrideDraft | FeatureOverrideRow[] | null | undefined} overrides
 *   An editor draft map, or request rows.
 * @returns {Array<FeatureOverrideRow | FeatureOverrideDraftRow>}
 */
export const canonicalFeatureOverrides = (overrides) => canonicalIdentityRows(overrides, checkedOverrideRow);

// The request rows of the current request's mode and records: the draft's
// edits, and the label text that a bulk label edit gives a feature without its
// own (`bulkLabelText`, {identityKey: text}). Rows of other records and of the
// other mode stay in the draft (R2).
/**
 * @param {FeatureOverrideDraft | null | undefined} overrides
 * @param {string} mode
 * @param {FeatureRequestRecord[]} [records]
 * @param {{ bulkLabelText?: Record<string, string> }} [options]
 * @returns {FeatureOverrideRow[]}
 */
export const requestFeatureOverrides = (overrides, mode, records = [], { bulkLabelText = {} } = {}) => {
  const rows = new Map();
  Object.entries(overrides || {}).forEach(([key, row]) => {
    if (rowBelongsToRequest(row, mode, records)) rows.set(key, row);
  });
  Object.entries(bulkLabelText || {}).forEach(([key, text]) => {
    const row = rows.get(key);
    const labelText = normalizeFeatureOverrideLabelText(text);
    if (row?.labelText || !labelText) return;
    const identity = parseFeatureIdentityKey(key);
    if (!rowBelongsToRequest(identity, mode, records)) return;
    rows.set(key, { ...(row || emptyOverrideRow(identity)), labelText });
  });
  return /** @type {FeatureOverrideDraftRow[]} */ (
    canonicalFeatureOverrides(Object.fromEntries([...rows].filter(([, row]) => hasEdit(row))))
  ).map(({ labelSourceText: _source, ...row }) => requestRow(row));
};

const emptyOverrideRow = ({ scope, recordKey, biologicalFeatureId }) => ({
  scope,
  recordKey,
  biologicalFeatureId,
  featureVisibility: null,
  labelVisibility: null,
  labelText: null,
  labelSourceText: null
});

// The draft value of one field for a feature, or null.
/**
 * @param {FeatureOverrideDraft | null | undefined} overrides
 * @param {Record<string, any> | null | undefined} feature
 * @param {string} field
 * @returns {any}
 */
export const featureOverrideValue = (overrides, feature, field) => {
  const key = featureIdentityKeyOf(feature);
  return key ? overrides?.[key]?.[field] ?? null : null;
};

// Sets fields of one feature's (or draft row's) draft row and deletes a row
// that no longer holds an edit or a label source text. Returns whether the
// draft changed.
/**
 * @param {FeatureOverrideDraft | null | undefined} overrides The draft, edited in place.
 * @param {Record<string, any> | null | undefined} feature
 * @param {Record<string, any> | null | undefined} patch
 * @returns {boolean}
 */
export const updateFeatureOverride = (overrides, feature, patch) => {
  const key = featureIdentityKeyOf(feature);
  if (!key || !overrides) return false;
  const next = { ...(overrides[key] || emptyOverrideRow(parseFeatureIdentityKey(key))) };
  Object.entries(patch || {}).forEach(([field, value]) => {
    if (!(field in next) || ['scope', 'recordKey', 'biologicalFeatureId'].includes(field)) return;
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
/**
 * @param {unknown} notices
 * @param {unknown} results
 * @returns {void}
 */
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
/**
 * @param {Record<string, any>} error
 * @param {Record<string, any> | null | undefined} request
 * @param {Array<Record<string, any>>} [features]
 * @returns {Record<string, any>}
 */
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
