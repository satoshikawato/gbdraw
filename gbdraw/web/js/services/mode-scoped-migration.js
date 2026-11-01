// @ts-check
// Session 46 keeps each diagram mode's settings and edits in its own slice
// (`modes.circular`, `modes.linear`, PD-OI-086). This state-free module owns
// the registry's slice shape (mode-scoped-settings.generated.js, written from
// gbdraw/web_support/mode_scoped_settings.py). Python alone splits a Session
// 27-44 draft into the two slices (`split_draft_into_modes` in
// gbdraw/session_io.py). The Session 46 Legend rekey is session-drawings-upgrade.js's.
import { MODE_SCOPED_SETTINGS } from '../mode-scoped-settings.generated.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { featureIdentityKey } from './feature-placement.js';

/** @typedef {'circular' | 'linear'} SliceMode */
/**
 * One registry row: a value that each mode's slice holds on its own.
 * @typedef {{ domain: string, path: string, modes: string, key: string | null, migrate: string }} ModeScopedRow
 */

/** @type {readonly SliceMode[]} */
export const SLICE_MODES = Object.freeze(['circular', 'linear']);
/** @type {readonly ModeScopedRow[]} */
export const MODE_SCOPED_ROWS = MODE_SCOPED_SETTINGS.rows;
/**
 * The draft fields (dotted, under `config`) of the other mode's `own`
 * registry rows: Python's split gives such a row to its mode's drawing alone
 * (`_split_setting`), so `mode`'s drawing reads a flat draft, and takes its
 * request's projection, without them.
 * @param {SliceMode} mode
 * @returns {readonly string[]}
 */
const otherModeOwnFields = (mode) => Object.freeze(MODE_SCOPED_ROWS
  .filter((row) => row.migrate === 'own' && row.modes !== mode && (row.domain === 'config' || row.domain.startsWith('config.')))
  .map((row) => [row.domain.slice('config.'.length), row.path].filter(Boolean).join('.')));
/** @type {Readonly<Record<SliceMode, readonly string[]>>} */
export const OTHER_MODE_OWN_FIELDS = Object.freeze({ circular: otherModeOwnFields('circular'), linear: otherModeOwnFields('linear') });
/** @type {readonly string[]} */
export const LOSAT_EXECUTION_FIELDS = MODE_SCOPED_SETTINGS.losatExecutionFields;

/** @param {unknown} value */
const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);
const fieldsError = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });

// Each slice container ('' is the slice itself) and the keys it may hold.
/** @type {Readonly<Record<string, ReadonlySet<string>>>} */
export const MODE_SLICE_CONTAINERS = (() => {
  /** @type {Record<string, Set<string>>} */
  const containers = { '': new Set() };
  MODE_SCOPED_ROWS.forEach((row) => {
    const parts = row.domain.split('.');
    parts.forEach((part, depth) => {
      const parent = parts.slice(0, depth).join('.');
      (containers[parent] ||= new Set()).add(part);
      containers[parts.slice(0, depth + 1).join('.')] ||= new Set();
    });
    containers[row.domain].add(row.path);
  });
  return Object.freeze(containers);
})();

// A slice holds only registry fields: the writer and both validators reject
// any other key (plan 4.2).
/**
 * The places of a slice that hold no registry field: a container that is no
 * object (`field` null), or a key it may not hold.
 * @param {unknown} value
 * @param {string} [container]
 * @returns {{ container: string, field: string | null }[]}
 */
export const modeSliceFieldViolations = (value, container = '') => {
  if (!isObject(value)) return [{ container, field: null }];
  const allowed = MODE_SLICE_CONTAINERS[container];
  return Object.entries(/** @type {Record<string, unknown>} */ (value)).flatMap(([field, child]) => {
    if (!allowed.has(field)) return [{ container, field }];
    const childContainer = container ? `${container}.${field}` : field;
    return Object.hasOwn(MODE_SLICE_CONTAINERS, childContainer) ? modeSliceFieldViolations(child, childContainer) : [];
  });
};
/**
 * The registry fields of `value` read as the slice `container` (its nested
 * containers the same way; a value of another kind is kept for its reader):
 * what a drawing's draft takes of its Session 50 entry's `editorState` and
 * `ui`, which also hold its Result set's fields.
 * @param {unknown} value
 * @param {string} container
 * @returns {unknown}
 */
export const modeSliceFieldsOf = (value, container) => {
  if (!isObject(value)) return value;
  const allowed = MODE_SLICE_CONTAINERS[container];
  return Object.fromEntries(Object.entries(/** @type {Record<string, unknown>} */ (value))
    .filter(([field]) => allowed.has(field))
    .map(([field, child]) => {
      const childContainer = container ? `${container}.${field}` : field;
      return [field, Object.hasOwn(MODE_SLICE_CONTAINERS, childContainer) ? modeSliceFieldsOf(child, childContainer) : child];
    }));
};
/** Save: a slice holds only registry fields. @param {unknown} value */
export const validateModeSliceFields = (value) => {
  if (modeSliceFieldViolations(value).length) throw fieldsError();
};

// The scoped draft rows of `mode` (a Session 41-44 migration's rows), keyed by
// identity pair and without `scope`; a row that names no mode is kept as it is.
/**
 * @param {unknown} rows
 * @param {SliceMode} mode
 */
export const unscopedDraftRows = (rows, mode) => {
  if (!isObject(rows)) return rows;
  /** @type {Array<[string, any]>} */
  const kept = [];
  Object.entries(/** @type {Record<string, any>} */ (rows)).forEach(([key, row]) => {
    if (isObject(row) && Object.hasOwn(row, 'scope')) {
      if (row.scope !== mode) return;
      const { scope: _scope, ...fields } = row;
      kept.push([featureIdentityKey(fields.recordKey, fields.biologicalFeatureId) || key, fields]);
    } else {
      kept.push([key, row]);
    }
  });
  return Object.fromEntries(kept);
};
