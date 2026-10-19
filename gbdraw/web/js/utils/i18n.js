// @ts-check
// The Web UI translation runtime. English is the source and the key: the first
// argument of `t()` is the English text the UI shows, so the English wording
// has one owner, the place that writes it. There is no English catalog. A
// catalog (`./i18n-catalogs/<locale>.js`, a frozen `Record<string, string>`)
// maps that English text to the shown text of one locale; a missing entry
// shows the English. Templates reach `t`, `tn` and `tx` through
// `app.config.globalProperties` (app.js).
//
// Never translated: the figure and every export, Session JSON, CLI / Source
// recipe / Exact replay text, diagnostics text, console output, identifiers
// (History step labels, `:key` values, `data-*`, ids), file and biological
// names. Keys are literal arguments, so tests/web/i18n-catalogs.test.mjs can
// list them.

/**
 * The BCP 47 tag of a UI locale: `en`, or the name of a catalog module.
 * @typedef {string} UiLocale
 */
/**
 * A locale's translations, keyed by the English text (`t`, and the `other`
 * form of `tn`) or by `context + '\u0004' + text` (`tx`).
 * @typedef {Readonly<Record<string, string>>} TranslationCatalog
 */
/**
 * Values of the `{name}` placeholders of a text.
 * @typedef {Readonly<Record<string, string | number>>} TranslationParams
 */
/**
 * Loads the catalog module of a locale.
 * @typedef {(locale: UiLocale) => Promise<{ default: TranslationCatalog }>} CatalogLoader
 */

export const SOURCE_LOCALE = 'en';
// Separates the context of `tx` from its text in a catalog key.
export const CONTEXT_SEPARATOR = '\u0004';
const PLACEHOLDER = /\{([A-Za-z_][A-Za-z0-9_]*)\}/g;
const LOCALE_TAG = /^[a-z]{2,3}(?:-[A-Za-z0-9]{2,8})*$/;

// The vendored Vue is a global of the page; a unit test without a page gets
// a plain holder, which is not reactive.
const vueRef = globalThis.window?.Vue?.ref;
/** @type {{ value: UiLocale }} */
const locale = vueRef ? vueRef(SOURCE_LOCALE) : { value: SOURCE_LOCALE };
/** @type {TranslationCatalog} */
let catalog = Object.freeze({});

/**
 * The shown UI locale. Templates that call `t()` read it, so they render
 * again when it changes; only `setUiLocale` writes it.
 * @type {{ readonly value: UiLocale }}
 */
export const uiLocale = locale;

/**
 * Replaces each `{name}` whose value is given; any other placeholder stays.
 * @param {string} text
 * @param {TranslationParams | undefined} params
 */
const substitute = (text, params) => (params
  ? text.replace(PLACEHOLDER, (/** @type {string} */ whole, /** @type {string} */ name) => (Object.hasOwn(params, name) ? String(params[name]) : whole))
  : text);

/**
 * The shown text of a catalog key, or `fallback` when the shown locale has none.
 * @param {string} key
 * @param {string} fallback
 */
const lookup = (key, fallback) => {
  // Read the locale first: a render that calls `t()` depends on it.
  if (locale.value === SOURCE_LOCALE) return fallback;
  return Object.hasOwn(catalog, key) ? catalog[key] : fallback;
};

/**
 * The shown text of an English UI text.
 * @param {string} text
 * @param {TranslationParams} [params]
 * @returns {string}
 */
export const t = (text, params) => substitute(lookup(text, text), params);

/**
 * The shown text of a count: `one` when `count` is 1 in English, else
 * `other`. The catalog key is `other`; `{count}` is the count unless `params`
 * gives it.
 * @param {number} count
 * @param {string} one
 * @param {string} other
 * @param {TranslationParams} [params]
 * @returns {string}
 */
export const tn = (count, one, other, params) => (
  substitute(lookup(other, count === 1 ? one : other), { count, ...params })
);

/**
 * The shown text of an English text that needs a different translation in
 * `context` (two meanings of one English word).
 * @param {string} context
 * @param {string} text
 * @param {TranslationParams} [params]
 * @returns {string}
 */
export const tx = (context, text, params) => substitute(lookup(context + CONTEXT_SEPARATOR + text, text), params);

/** @type {CatalogLoader} */
const importCatalog = (tag) => import(`./i18n-catalogs/${tag}.js`);

/**
 * Shows the UI in `next`: loads its catalog, then switches. `en` needs no
 * catalog. Nothing in the UI calls it yet (PD-OI-106 decides the menu).
 * @param {UiLocale} next
 * @param {CatalogLoader} [load] the catalog source; tests pass an in-memory one
 * @returns {Promise<void>}
 */
export const setUiLocale = async (next, load = importCatalog) => {
  if (!LOCALE_TAG.test(next)) throw new Error(`Unknown UI locale: ${next}`);
  catalog = next === SOURCE_LOCALE ? Object.freeze({}) : Object.freeze({ ...(await load(next)).default });
  locale.value = next;
};
