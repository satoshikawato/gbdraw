// @ts-check
// The Web UI translation runtime. English is the source and the key: the first
// argument of `t()` is the English text the UI shows, so the English wording
// has one owner, the place that writes it. There is no English catalog. A
// catalog (`./i18n-catalogs/<locale>.js`, a frozen `Record<string, string>`)
// maps that English text to the shown text of one locale; a missing entry
// shows the English. Templates reach `t`, `tn`, `tx` and `tm` through
// `app.config.globalProperties` (app.js). A text kept until it is shown (a
// table evaluated at import, a value in state, a value code compares) stays
// English, marked with `msg('…')`, and is translated where it is shown with
// `tm(value)`.
//
// Never translated: the figure and every export, Session JSON, CLI / Source
// recipe / Exact replay text, diagnostics text, console output, identifiers
// (History step labels, `:key` values, `data-*`, ids), file and biological
// names. Keys are literal arguments of `t`, `tn`, `tx` and `msg`, so
// tests/web/i18n-catalogs.test.mjs can list them.

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
 * A piece of a text shown with elements inside it (`I18nText` in
 * components.js): plain text, or the slot that a `{name}` placeholder names.
 * @typedef {{ text: string } | { slot: string }} TextPart
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

/**
 * Marks an English UI text that is kept until it is shown: a value of a table
 * evaluated at import, a value kept in state, or a value that code compares
 * (a History step label, a section title). Returns `text` unchanged, so the
 * value stays English and comparisons keep working; the catalog check reads
 * its literal argument as a key. Translate the value where it is shown with
 * `tm(value)`.
 * @param {string} text
 * @returns {string}
 */
export const msg = (text) => text;

/**
 * The shown text of a value that `msg()` marked, read where it is shown. A
 * value that no `msg()` or `t()` call makes a key (a file name, a warning
 * Python words) shows as written.
 * @param {string} text
 * @param {TranslationParams} [params]
 * @returns {string}
 */
export const tm = (text, params) => substitute(lookup(text, text), params);

/**
 * Splits a shown text at the placeholders that `slots` names. Another
 * placeholder takes its value from `params`, or stays as written. The values
 * are filled after the split, so a value never becomes a slot.
 * @param {string} text the shown text of a key (`t('…')` without values)
 * @param {readonly string[]} slots
 * @param {TranslationParams} [params]
 * @returns {TextPart[]}
 */
export const textParts = (text, slots, params) => {
  /** @type {TextPart[]} */
  const parts = [];
  let pending = '';
  let at = 0;
  for (const match of text.matchAll(PLACEHOLDER)) {
    const [whole, name] = match;
    pending += text.slice(at, match.index);
    at = (match.index ?? 0) + whole.length;
    if (slots.includes(name)) {
      if (pending) parts.push({ text: pending });
      parts.push({ slot: name });
      pending = '';
    } else {
      pending += params && Object.hasOwn(params, name) ? String(params[name]) : whole;
    }
  }
  pending += text.slice(at);
  if (pending) parts.push({ text: pending });
  return parts;
};

/** @type {CatalogLoader} */
const importCatalog = (tag) => import(`./i18n-catalogs/${tag}.js`);

// The latest `setUiLocale` call; an earlier call whose catalog arrives later is dropped.
let latestRequest = 0;

/**
 * Shows the UI in `next`: loads its catalog, then switches. `en` needs no
 * catalog. When calls overlap, only the latest one applies. Nothing in the UI
 * calls it yet (PD-OI-106 decides the menu).
 * @param {UiLocale} next
 * @param {CatalogLoader} [load] the catalog source; tests pass an in-memory one
 * @returns {Promise<boolean>} whether this call switched the locale
 */
export const setUiLocale = async (next, load = importCatalog) => {
  if (!LOCALE_TAG.test(next)) throw new Error(`Unknown UI locale: ${next}`);
  latestRequest += 1;
  const request = latestRequest;
  const loaded = next === SOURCE_LOCALE ? {} : (await load(next)).default;
  if (request !== latestRequest) return false;
  catalog = Object.freeze({ ...loaded });
  locale.value = next;
  return true;
};
