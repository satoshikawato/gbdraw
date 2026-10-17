// @ts-check
// The one reader of Vue's `toRaw`. The vendored Vue is a global of the page
// (index.html); a Worker or a unit test without a page may not set it.
/**
 * The object behind a Vue proxy, or the value itself outside Vue.
 * @template T
 * @param {T} value
 * @returns {T}
 */
export const toRaw = (value) => globalThis.window?.Vue?.toRaw?.(value) ?? value;
