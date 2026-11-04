// @ts-check
// Python's readings of the stored config leaves a form cannot hold as stored (a wrap length typed with a unit such
// as `"5kb"`, a figure width written `"180 MM"`). Load asks Python, through the Worker's configuration check, for
// every such leaf before its readers run (services/config.js `readStoredLeaves`); each reader of a leaf then looks
// its stored value up here (services/linear-wrap-form.js `wrapFormValue`, services/figure-size-form.js
// `figureLengthFormValue`). A value Python refuses has no reading, and its reader drops it with a Load repair row.

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
