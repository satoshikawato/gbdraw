// @ts-check
// Python's readings of the stored config leaves a form cannot hold as stored (a wrap length typed with a unit such
// as `"5kb"`, a figure width written `"180 MM"`). Load asks Python, through the Worker's configuration check, for
// every such leaf before its readers run (services/config.js `readStoredLeaves`); each reader of a leaf then looks
// its stored value up here (services/linear-wrap-form.js `wrapFormValue`, services/figure-size-form.js
// `figureLengthFormValue`). A value Python refuses has no reading, and its reader drops it with a Load repair row.
import { storedRecordSourceResourceIds } from './input-types.js';

/**
 * Python's canonical value of each stored leaf the form cannot hold as stored, by `storedLeafKey` of the leaf's
 * config path and its stored value.
 * @typedef {ReadonlyMap<string, unknown>} StoredLeafReadings
 */
/** @type {StoredLeafReadings} */
export const NO_STORED_READINGS = new Map();

/**
 * @param {string} path The leaf's config path (`canvas.figure_width`).
 * @param {unknown} raw The value as stored.
 */
export const storedLeafKey = (path, raw) => JSON.stringify([path, raw]);

/**
 * A record of a stored render request as Load holds it before the request is read: any field may be missing or of
 * another type.
 * @typedef {{ readonly selector?: Readonly<Record<string, unknown>> | null,
 *   readonly source?: Readonly<Record<string, unknown>> | null }} UnreadStoredRecord
 */
/**
 * A stored request record whose selector is a record ID, with the key of its reading: a record text (a Linear card,
 * the Circular Record select, a Circular record range) reads an ID such as `#2` or `#+2` as a record number, so Load
 * reads the record's place in its file (`#k`, the text the Record select gives such a record;
 * `readStoredRecordIdPlaces`, services/record-discovery.js) and the text takes it (OV-3311). The one builder of the
 * key, for its writer and its reader. Null for any other record, and for a source of no file.
 * @param {UnreadStoredRecord | null | undefined} record
 * @returns {{ recordId: string, resourceIds: string[], key: string } | null}
 */
export const storedRecordIdOf = (record) => {
  const recordId = record?.selector?.kind === 'recordId' ? record.selector.value : null;
  const resourceIds = storedRecordSourceResourceIds(record?.source);
  if (typeof recordId !== 'string' || resourceIds.length === 0) return null;
  return { recordId, resourceIds, key: storedLeafKey('records.selector', [resourceIds, recordId]) };
};
