// @ts-check
// The entry of a table whose keys are data (file text, a Session id): only an
// own key names an entry, so `constructor` or `__proto__` reads none.
/**
 * @template T
 * @param {Readonly<Record<string, T>> | null | undefined} table
 * @param {string} key
 * @returns {T | undefined} The entry, or `undefined` when the table has no such own key.
 */
export const ownEntry = (table, key) => {
  if (!table) return undefined;
  // Read first. On a Vue reactive table `table[key]` is the dependency (also for
  // a key that is not there yet); `Object.hasOwn` alone is none, so a computed
  // that found no entry would not run again when the entry is first written.
  const value = table[key];
  return Object.hasOwn(table, key) ? value : undefined;
};
