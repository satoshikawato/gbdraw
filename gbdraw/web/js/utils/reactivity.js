// @ts-check
// Typed reactivity over the page's vendored Vue (index.html loads it as the
// global `window.Vue`, which web-globals.d.ts leaves untyped). The state
// module and the drawing slices (services/drawing/) build their members here,
// so a member's type is written once, where it is created. Vue is read at each
// call: a unit test sets its own `window.Vue` before it creates a member.

/**
 * A writable Vue ref.
 * @template T
 * @typedef {{ value: T }} Ref
 */
/**
 * A value whose members, and their members, are read-only to the type
 * checker. Functions stay callable; a set or a map is its read-only view.
 * Only the type is read-only: no Vue `readonly` proxy is made, because the
 * production build drops a write to one without a warning.
 * @template T
 * @typedef {T extends string | number | boolean | bigint | symbol | null | undefined ? T
 *   : T extends (...args: never[]) => unknown ? T
 *   : T extends ReadonlySet<infer V> ? ReadonlySet<DeepReadonly<V>>
 *   : T extends ReadonlyMap<infer K, infer V> ? ReadonlyMap<K, DeepReadonly<V>>
 *   : { readonly [K in keyof T]: DeepReadonly<T[K]> }} DeepReadonly
 */
/**
 * A ref whose value its readers may not change.
 * @template T
 * @typedef {{ readonly value: DeepReadonly<T> }} ReadonlyRef
 */

/**
 * @template T
 * @param {T} value
 * @returns {Ref<T>}
 */
export const ref = (value) => window.Vue.ref(value);

/**
 * A reactive object over `value` (Vue `reactive`).
 * @template {object} T
 * @param {T} value
 * @returns {T}
 */
export const reactive = (value) => window.Vue.reactive(value);

/**
 * @template T
 * @param {() => T} getter
 * @returns {{ readonly value: T }}
 */
export const computed = (getter) => window.Vue.computed(getter);
