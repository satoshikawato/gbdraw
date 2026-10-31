// @ts-check
// The Web UI translation runtime. English is the source and the key: the first
// argument of `t()` is the English text the UI shows, so the English wording
// has one owner, the place that writes it. There is no English catalog. A
// catalog (`./i18n-catalogs/<locale>.js`, a frozen `TranslationCatalog`)
// maps that English text to the shown text of one locale; a missing entry
// shows the English. A form of a `tn()` pair may instead map to its text per
// plural category (`PluralForms`). Templates reach `t`, `tn`, `tx` and `tm` through
// `app.config.globalProperties` (app.js). A text kept until it is shown (a
// table evaluated at import, a value in state, a value code compares) stays
// English, marked with `msg('…')`, and is translated where it is shown with
// `tm(value)`.
//
// Never translated: the figure and every export, Session JSON, CLI / Source
// recipe / Exact replay text, diagnostics text, console output, identifiers
// (History step labels, `:key` values, `data-*`, ids), file and biological
// names. Keys are literal arguments of `t`, `tx`, `msg` and both forms of
// `tn`, so tests/web/i18n-catalogs.test.mjs can list them.
//
// The catalogs are generated, never edited by hand; their wording follows
// docs/internal/I18N_STYLE.md and docs/internal/I18N_GLOSSARY.tsv.

/**
 * The BCP 47 tag of a UI locale: `en`, or the name of a catalog module. It is
 * also the page's `lang`.
 * @typedef {'en' | 'id' | 'es' | 'fr' | 'pt' | 'vi' | 'th' | 'ja' | 'zh-Hans' | 'zh-Hant' | 'ko'} UiLocale
 */
/**
 * A language the Web UI can be shown in (PD-OI-106).
 * @typedef {object} UiLocaleOption
 * @property {UiLocale} tag
 * @property {string} name The language's name in itself (an endonym); never translated.
 * @property {string} sentenceSeparator What separates two sentences shown one after the other.
 */
/**
 * The translation of a form of a `tn` pair in a language whose plural rule
 * has more categories than English (`few` and `many` in Polish, Russian and
 * Ukrainian): the text per CLDR plural category (Intl.PluralRules); a
 * category without a text takes `other`. Only the value of a key that only
 * `tn` calls produce may be an object; every other value is one string.
 * @typedef {Readonly<{ other: string, zero?: string, one?: string, two?: string, few?: string, many?: string }>} PluralForms
 */
/**
 * A locale's translations, keyed by the English text (`t`, and each form
 * of `tn`) or by `context + '\u0004' + text` (`tx`).
 * @typedef {Readonly<Record<string, string | PluralForms>>} TranslationCatalog
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

/** @type {UiLocale} */
export const SOURCE_LOCALE = 'en';
/**
 * The languages of the header Language menu, in its order: English, the
 * other Latin-script languages by name, Thai, then East Asian (PD-OI-106).
 * Japanese and Chinese put no space between sentences.
 * @type {readonly Readonly<UiLocaleOption>[]}
 */
export const UI_LOCALES = Object.freeze([
  Object.freeze({ tag: SOURCE_LOCALE, name: 'English', sentenceSeparator: ' ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('id'), name: 'Bahasa Indonesia', sentenceSeparator: ' ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('es'), name: 'Español', sentenceSeparator: ' ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('fr'), name: 'Français', sentenceSeparator: ' ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('pt'), name: 'Português', sentenceSeparator: ' ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('vi'), name: 'Tiếng Việt', sentenceSeparator: ' ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('th'), name: 'ไทย', sentenceSeparator: ' ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('ja'), name: '日本語', sentenceSeparator: '' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('zh-Hans'), name: '简体中文', sentenceSeparator: '' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('zh-Hant'), name: '繁體中文', sentenceSeparator: '' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('ko'), name: '한국어', sentenceSeparator: ' ' })
]);
// Separates the context of `tx` from its text in a catalog key.
export const CONTEXT_SEPARATOR = '\u0004';
/**
 * The CLDR plural categories, the keys a `PluralForms` may have, in CLDR
 * order; tools/i18n checks and prints catalogs by it.
 * @type {readonly Intl.LDMLPluralRule[]}
 */
export const PLURAL_CATEGORIES = Object.freeze(/** @type {Intl.LDMLPluralRule[]} */ (['zero', 'one', 'two', 'few', 'many', 'other']));
const PLACEHOLDER = /\{([A-Za-z_][A-Za-z0-9_]*)\}/g;

// The vendored Vue is a global of the page; a unit test without a page gets
// a plain holder, which is not reactive.
const vueRef = globalThis.window?.Vue?.ref;
/** @type {{ value: UiLocale }} */
const locale = vueRef ? vueRef(SOURCE_LOCALE) : { value: SOURCE_LOCALE };

/**
 * The UI locale a stored or chosen value names, or null for any other value.
 * @param {unknown} value
 * @returns {UiLocale | null}
 */
export const uiLocaleOf = (value) => UI_LOCALES.find((option) => option.tag === value)?.tag ?? null;

/** The shown locale's entry of UI_LOCALES; reading it makes a render depend on the locale. */
const shownLocale = () => UI_LOCALES.find((option) => option.tag === locale.value) ?? UI_LOCALES[0];
/** @type {TranslationCatalog} */
let catalog = Object.freeze({});
// The plural rules of English, for the English a missing translation shows,
// and of the shown locale; `setUiLocale` replaces the second.
const ENGLISH_PLURALS = new Intl.PluralRules(SOURCE_LOCALE);
let shownPlurals = ENGLISH_PLURALS;

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
  if (locale.value === SOURCE_LOCALE || !Object.hasOwn(catalog, key)) return fallback;
  const value = catalog[key];
  // Plural forms are a catalog error on a key that is not only a `tn` form (the catalog test refuses them).
  return typeof value === 'string' ? value : value.other;
};

/**
 * The shown text of an English UI text.
 * @param {string} text
 * @param {TranslationParams} [params]
 * @returns {string}
 */
export const t = (text, params) => substitute(lookup(text, text), params);

/**
 * The English form of `count` under `rules`: `one` when they select the
 * category `one`, else `other` (every other category takes `other`).
 * @param {Intl.PluralRules} rules
 * @param {number} count
 * @param {string} one
 * @param {string} other
 */
const pluralForm = (rules, count, one, other) => (rules.select(count) === 'one' ? one : other);

/**
 * The translation of `count` in `entries` under `rules`, or null when
 * `entries` lack the key. The category the rules select picks the key, the
 * English text of a form (`one` for the category `one`, else `other`), and,
 * when the value is `PluralForms`, the text of that category, else its
 * `other`. A string value is the text of every category.
 * @param {Intl.PluralRules} rules
 * @param {TranslationCatalog} entries
 * @param {number} count
 * @param {string} one
 * @param {string} other
 * @returns {string | null}
 */
export const pluralTranslation = (rules, entries, count, one, other) => {
  const category = rules.select(count);
  const key = category === 'one' ? one : other;
  if (!Object.hasOwn(entries, key)) return null;
  const value = entries[key];
  return typeof value === 'string' ? value : value[category] ?? value.other;
};

/**
 * The shown text of a count. The shown locale's plural rules
 * (Intl.PluralRules) choose the form, and that form's English text is the
 * catalog key: `fr` and `pt` take `one` for 0 and 1, `en` and `es` for 1
 * only, `ja`, `zh-Hans`, `zh-Hant` and `ko` never. A value with plural
 * forms gives the text of the category (`pluralTranslation`). A missing
 * translation shows the form English chooses. `{count}` is the count unless
 * `params` gives it.
 * @param {number} count
 * @param {string} one
 * @param {string} other
 * @param {TranslationParams} [params]
 * @returns {string}
 */
export const tn = (count, one, other, params) => {
  const english = pluralForm(ENGLISH_PLURALS, count, one, other);
  // Read the locale first: a render that calls `tn()` depends on it.
  const shown = locale.value === SOURCE_LOCALE ? english : pluralTranslation(shownPlurals, catalog, count, one, other) ?? english;
  return substitute(shown, { count, ...params });
};

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
 * The shown text, in `context`, of an English text derived from a value set
 * that another owner generates, so that no literal of it is in the source:
 * the `tx()` counterpart of `tm()`. The module that derives the text also
 * lists its keys, which the catalog check reads.
 * @param {string} context
 * @param {string} text
 * @returns {string}
 */
export const tmx = (context, text) => lookup(context + CONTEXT_SEPARATOR + text, text);

/**
 * What separates two sentences shown one after the other in the shown locale:
 * a space in English, Spanish, French, Portuguese and Korean, nothing in
 * Japanese and Chinese (after 。). Templates
 * put it between two keys: `{{ t('A.') }}{{ sentenceSeparator() }}{{ t('B.') }}`.
 * @returns {string}
 */
export const sentenceSeparator = () => shownLocale().sentenceSeparator;

/**
 * Sentences shown one after the other, each already a shown text.
 * @param {...string} sentences
 * @returns {string}
 */
export const joinSentences = (...sentences) => sentences.join(sentenceSeparator());

/**
 * Shown texts as one list: "A", "A and B", "A, B, and C" in English, the
 * locale's own list pattern (Intl.ListFormat) in the other locales.
 * @param {readonly string[]} items
 * @returns {string}
 */
export const listText = (items) => {
  const tag = shownLocale().tag;
  if (tag !== SOURCE_LOCALE) return new Intl.ListFormat(tag, { type: 'conjunction' }).format(items);
  if (items.length < 3) return items.join(' and ');
  return `${items.slice(0, -1).join(', ')}, and ${items.at(-1)}`;
};

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
 * Shows the UI in `next`: loads its catalog, then switches, and sets the
 * page's `lang` (the one writer of `<html lang>`). `en` needs no catalog. When
 * calls overlap, only the latest one applies; a catalog that fails to load
 * rejects and changes nothing. The header Language menu and the stored
 * language at boot call it (app/ui-language.js).
 * @param {string} next
 * @param {CatalogLoader} [load] the catalog source; tests pass an in-memory one
 * @returns {Promise<boolean>} whether this call switched the locale
 */
export const setUiLocale = async (next, load = importCatalog) => {
  const tag = uiLocaleOf(next);
  if (!tag) throw new Error(`Unknown UI locale: ${next}`);
  latestRequest += 1;
  const request = latestRequest;
  const loaded = tag === SOURCE_LOCALE ? {} : (await load(tag)).default;
  if (request !== latestRequest) return false;
  catalog = Object.freeze({ ...loaded });
  shownPlurals = tag === SOURCE_LOCALE ? ENGLISH_PLURALS : new Intl.PluralRules(tag);
  locale.value = tag;
  // A unit test without a page has no document.
  globalThis.document?.documentElement.setAttribute('lang', tag);
  return true;
};
