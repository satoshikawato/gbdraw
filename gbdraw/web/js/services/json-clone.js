// @ts-check
/**
 * A value as `JSON.parse` returns it.
 * @typedef {null | boolean | number | string | JsonValue[] | { [key: string]: JsonValue }} JsonValue
 */
/**
 * @template T
 * @param {T} value
 * @returns {T} The clone; a nullish `value` is returned as it is.
 */
export const cloneJsonData = (value) => {
  if (value === null || value === undefined) return value;
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
