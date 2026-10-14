// @ts-check

// Readers of one stored value: the value the app uses, and whether reading it
// lost anything (BA-D25). Load lists a lost reading before it applies a file;
// the apply and a History restore read with the same reader. An unset value
// (absent, null, empty text) takes its default without a loss, and so do
// case, outer spaces and a listed alias.

/**
 * A value inside a read list or object that its reader replaced; `at` is its
 * place under the field (`[0].height`, `.name.font_size`).
 * @typedef {{ at: string, stored: unknown, used: unknown }} FieldLeaf
 */
/**
 * @template T
 * @typedef {{ value: T, lost: boolean, dropped?: number, leaves?: readonly FieldLeaf[], kept?: boolean }} FieldReading
 *   `dropped`: how many entries of a list the reading dropped; `leaves`: the
 *   values inside it the reading replaced (a reading that lost something and
 *   names neither lost the field as a whole); `kept`: the stored value is not
 *   used, and the field keeps the value it holds.
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
  // A boolean reads as 1 or 0, as `numberOr` reads it, but is no number (rule 2).
  return Number.isFinite(numeric) && accepts(numeric) ? { value: numeric, lost: typeof value === 'boolean' } : { value: null, lost: true };
};
export const positiveNumberOrNull = optionalNumber((numeric) => numeric > 0);
export const nonNegativeNumberOrNull = optionalNumber((numeric) => numeric >= 0);

/** @param {number} numeric */
export const isPositiveInteger = (numeric) => Number.isInteger(numeric) && numeric > 0;
/** @param {number} numeric */
export const isNonNegativeInteger = (numeric) => Number.isInteger(numeric) && numeric >= 0;

/**
 * A finite number `accepts` takes, else `fallback`. An unset value (a blank
 * is null, R7) is `fallback` without a loss, never the 0 `Number()` makes of it.
 * @template F
 * @param {(numeric: number) => boolean} accepts
 * @param {F} fallback
 * @returns {FieldReader<number | F>}
 */
export const numberOr = (accepts, fallback) => (value) => {
  if (isUnsetValue(value)) return { value: fallback, lost: false };
  const numeric = Number(value);
  // A boolean reads as 1 or 0, as it always did, but is no number (rule 2).
  return Number.isFinite(numeric) && accepts(numeric)
    ? { value: numeric, lost: typeof value === 'boolean' } : { value: fallback, lost: true };
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

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isObjectValue = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

/** @param {unknown} value Its text without case and outer spaces. */
const foldedText = (value) => String(value).trim().toLowerCase();

/**
 * Whether the value `used` for a stored value keeps its meaning: the same
 * value, or its own text (`"5"` as 5, case and outer spaces aside). A boolean
 * read as text (`true` as "true") or a list or object read as a value is
 * another kind (rule 2).
 * @param {unknown} stored
 * @param {unknown} used
 */
export const keepsMeaning = (stored, used) => stored === used || (
  used !== undefined && used !== null && typeof used !== 'object'
  && stored !== null && typeof stored !== 'object' && typeof stored !== 'boolean'
  && foldedText(stored) === foldedText(used));

/**
 * The values of `stored` that its reading `read` does not keep: a value read
 * as other text, a list or object read as another kind, and a key or entry
 * the reading dropped. An unset stored value takes its default without a
 * loss, and a value read as its own text (`"5"` as 5, case and outer spaces
 * aside) keeps its meaning.
 * @param {unknown} stored
 * @param {unknown} read
 * @param {string} [at] Its place under the field (`[0].mark`).
 * @returns {FieldLeaf[]}
 */
export const changedLeaves = (stored, read, at = '') => {
  if (isUnsetValue(stored)) return [];
  const list = Array.isArray(stored);
  if (list || isObjectValue(stored)) {
    if (list ? !Array.isArray(read) : !isObjectValue(read)) return [{ at, stored, used: read }];
    const source = /** @type {Record<string, unknown>} */ (stored);
    const target = /** @type {Record<string, unknown>} */ (read);
    return Object.keys(source).flatMap((key) => {
      const keyAt = list ? `${at}[${key}]` : `${at}.${key}`;
      return Object.hasOwn(target, key) ? changedLeaves(source[key], target[key], keyAt)
        : changedLeaves(source[key], undefined, keyAt);
    });
  }
  return keepsMeaning(stored, read) ? [] : [{ at, stored, used: read }];
};

/**
 * Records in `losses` a stored value a reader replaced with `used`, unless the
 * stored value was unset or `used` keeps its meaning (OV-758: the track-slot
 * readers report what they replace this way).
 * @param {FieldLeaf[] | null | undefined} losses
 * @param {string} at
 * @param {unknown} stored
 * @param {unknown} used
 */
export const noteReplacedLeaf = (losses, at, stored, used) => {
  if (losses && !isUnsetValue(stored) && !keepsMeaning(stored, used)) losses.push({ at, stored, used });
};
