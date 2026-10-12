// @ts-check
/**
 * @template T
 * @param {T} value
 * @returns {NonNullable<T>} The clone; a nullish `value` is returned as it is, so callers pass a value they know is set.
 */
export const cloneJsonData = (value) => {
  // A nullish value passes through; callers type their input as set (see @returns).
  if (value === null || value === undefined) return /** @type {NonNullable<T>} */ (value);
  return JSON.parse(JSON.stringify(value));
};

/**
 * @template T, F
 * @param {T | undefined} value
 * @param {F} fallback
 * @returns {T | F}
 */
export const cloneJsonValue = (value, fallback) => {
  if (value === undefined) return fallback;
  try {
    return cloneJsonData(value);
  } catch (_err) {
    return fallback;
  }
};
