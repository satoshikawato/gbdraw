// @ts-check

// Readers of one stored value: the value the app uses, and whether reading it
// lost anything (BA-D25). Load lists a lost reading before it applies a file;
// the apply and a History restore read with the same reader. An unset value
// (absent, null, empty text) takes its default without a loss, and so do
// case, outer spaces and a listed alias.

/**
 * @template T
 * @typedef {{ value: T, lost: boolean, dropped?: number }} FieldReading
 *   `dropped`: how many entries of a list the reading dropped.
 */
/**
 * @template T
 * @typedef {(value: unknown) => FieldReading<T>} FieldReader
 */

/** @param {unknown} value */
export const isUnsetValue = (value) => value === undefined || value === null || value === '';

/**
 * One of `values`, read without case and outer spaces (and with `-` read as
 * `_` when `foldHyphens`); an alias reads as its value; anything else is
 * `fallback`. A falsy value is unset.
 * @template {string} V
 * @template F
 * @param {readonly V[]} values
 * @param {F} fallback
 * @param {{ aliases?: Readonly<Record<string, V>>, foldHyphens?: boolean }} [options]
 * @returns {FieldReader<V | F>}
 */
export const textChoice = (values, fallback, { aliases = {}, foldHyphens = false } = {}) => (value) => {
  const folded = String(value || '').trim().toLowerCase();
  const normalized = foldHyphens ? folded.replace(/-/g, '_') : folded;
  const resolved = Object.hasOwn(aliases, normalized) ? aliases[normalized] : normalized;
  const known = values.find((candidate) => candidate === resolved);
  return known === undefined ? { value: fallback, lost: normalized !== '' } : { value: known, lost: false };
};

/**
 * A number `accepts` takes, else null; unset and `auto` read as null.
 * @param {(numeric: number) => boolean} accepts
 * @returns {FieldReader<number | null>}
 */
const optionalNumber = (accepts) => (value) => {
  if (isUnsetValue(value) || String(value).trim().toLowerCase() === 'auto') return { value: null, lost: false };
  const numeric = Number(value);
  return Number.isFinite(numeric) && accepts(numeric) ? { value: numeric, lost: false } : { value: null, lost: true };
};
export const positiveNumberOrNull = optionalNumber((numeric) => numeric > 0);
export const nonNegativeNumberOrNull = optionalNumber((numeric) => numeric >= 0);

/** @param {number} numeric */
export const isPositiveInteger = (numeric) => Number.isInteger(numeric) && numeric > 0;
/** @param {number} numeric */
export const isNonNegativeInteger = (numeric) => Number.isInteger(numeric) && numeric >= 0;

/**
 * A finite number `accepts` takes, else `fallback`. An unset value is
 * `fallback`, or with `coerceUnset` the number `Number()` makes of it (0 for
 * null and empty text), as the reader always took it.
 * @template F
 * @param {(numeric: number) => boolean} accepts
 * @param {F} fallback
 * @param {{ coerceUnset?: boolean }} [options]
 * @returns {FieldReader<number | F>}
 */
export const numberOr = (accepts, fallback, { coerceUnset = false } = {}) => (value) => {
  const unset = isUnsetValue(value);
  if (unset && !coerceUnset) return { value: fallback, lost: false };
  const numeric = Number(value);
  return Number.isFinite(numeric) && accepts(numeric) ? { value: numeric, lost: false } : { value: fallback, lost: !unset };
};

/** @param {unknown} value */
const notBoolean = (value) => !isUnsetValue(value) && typeof value !== 'boolean';
/** True only for `true`. @type {FieldReader<boolean>} */
export const onlyTrue = (value) => ({ value: value === true, lost: notBoolean(value) });
/** False only for `false`. @type {FieldReader<boolean>} */
export const trueUnlessFalse = (value) => ({ value: value !== false, lost: notBoolean(value) });
/** The value's truth. @type {FieldReader<boolean>} */
export const truthy = (value) => ({ value: Boolean(value), lost: notBoolean(value) });

/**
 * One of `values` exactly, else `fallback`; a falsy value is unset.
 * @template {string} V
 * @template F
 * @param {readonly V[]} values
 * @param {F} fallback
 * @returns {FieldReader<V | F>}
 */
export const exactChoice = (values, fallback) => (value) => {
  const known = values.find((candidate) => candidate === value);
  return known === undefined ? { value: fallback, lost: Boolean(value) } : { value: known, lost: false };
};
