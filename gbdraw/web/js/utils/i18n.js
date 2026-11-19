// @ts-check
// The Web UI translation runtime. English is the source and the key: the first
// argument of `t()` is the English text the UI shows, so the English wording
// has one owner, the place that writes it. There is no English catalog. A
// catalog (`./i18n-catalogs/<locale>.js`, a frozen `TranslationCatalog`)
// maps that English text to the shown text of one locale; a missing entry
// shows the English. The `other` form of a `tn()` pair may instead map to its
// text per plural category (`PluralForms`). Templates reach `t`, `tn`, `tx` and `tm` through
// `app.config.globalProperties` (app.js). A text kept until it is shown (a
// table evaluated at import, a value in state, a value code compares) stays
// English, marked with `msg('…')`, and is translated where it is shown with
// `tm(value)`. A text that stays shown while the UI language may change (a
// notice) is kept as a `Wording`, a function that words it, and is worded
// where it is shown.
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
 * @typedef {'en' | 'id' | 'de' | 'es' | 'fr' | 'it' | 'pl' | 'pt' | 'vi' | 'tr' | 'ru' | 'uk' | 'ar' | 'fa' | 'th' | 'ja' | 'zh-Hans' | 'zh-Hant' | 'ko'} UiLocale
 */
/**
 * A language the Web UI can be shown in (PD-OI-106).
 * @typedef {object} UiLocaleOption
 * @property {UiLocale} tag
 * @property {string} name The language's name in itself (an endonym); never translated.
 * @property {string} sentenceSeparator What separates two sentences shown one after the other.
 * @property {UiTextStyle} textStyle How the page sets its small UI text (the page's `data-ui-text-style`).
 * @property {TextDirection} dir The direction the page reads in this language: the page's `dir`.
 * @property {string} listSeparator What separates shown items listed without a conjunction (`itemListText`).
 */
/**
 * How the page sets the UI text of a locale; index.html reads it from the
 * page's `data-ui-text-style`, and the figure is never affected.
 * `cased`: as designed, with upper-casing and tracking.
 * `no-caps`: no upper-casing or tracking (no letter case, or stacked diacritics).
 * `tall-marks`: no-caps, and no text under 12px and tight line heights widened
 * to 1.5, for the marks above and below the letters (Thai).
 * @typedef {'cased' | 'no-caps' | 'tall-marks'} UiTextStyle
 */
/**
 * The direction a language is read in (the HTML `dir` attribute).
 * @typedef {'ltr' | 'rtl'} TextDirection
 */
/**
 * The translation of the `other` form of a `tn` pair in a language whose
 * plural rule has more categories than English (`few` and `many` in Polish,
 * Russian and Ukrainian): the text per CLDR plural category
 * (Intl.PluralRules) for every count except exactly 1; a category without a
 * text takes `other`. The `one` form's value is then the text for exactly 1
 * only, and `one` here is the text of the other counts the rule calls `one`
 * (21, 31, … in Russian and Ukrainian). Only the value of a key that is only
 * ever the `other` form of `tn` calls may be an object; every other value is
 * one string.
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
 * How to word a text that stays shown while the UI language may change (a
 * notice and its lines): a function that returns the shown text by calling
 * `t`, `tn`, `tx`, `tm` or `joinSentences` with literal keys and the values
 * it holds (file and font names, counts, Python's words). The place that shows
 * the text calls it, so a UI language change rewords the shown text.
 * @typedef {() => string} Wording
 */
/**
 * Loads the catalog module of a locale.
 * @typedef {(locale: UiLocale) => Promise<{ default: TranslationCatalog }>} CatalogLoader
 */

/** @type {UiLocale} */
export const SOURCE_LOCALE = 'en';
/**
 * The languages of the header Language menu, in its order: English, the
 * Latin-script languages by name, the Cyrillic ones by name, Thai, then East
 * Asian (PD-OI-106). Japanese and Chinese put no space between sentences;
 * Persian and Arabic separate listed items with the Arabic comma ، (I18N_STYLE:
 * Persian and Arabic punctuation).
 * @type {readonly Readonly<UiLocaleOption>[]}
 */
export const UI_LOCALES = Object.freeze([
  Object.freeze({ tag: SOURCE_LOCALE, name: 'English', sentenceSeparator: ' ', textStyle: 'cased', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('id'), name: 'Bahasa Indonesia', sentenceSeparator: ' ', textStyle: 'cased', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('de'), name: 'Deutsch', sentenceSeparator: ' ', textStyle: 'cased', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('es'), name: 'Español', sentenceSeparator: ' ', textStyle: 'cased', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('fr'), name: 'Français', sentenceSeparator: ' ', textStyle: 'cased', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('it'), name: 'Italiano', sentenceSeparator: ' ', textStyle: 'cased', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('pl'), name: 'Polski', sentenceSeparator: ' ', textStyle: 'cased', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('pt'), name: 'Português', sentenceSeparator: ' ', textStyle: 'cased', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('vi'), name: 'Tiếng Việt', sentenceSeparator: ' ', textStyle: 'no-caps', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('tr'), name: 'Türkçe', sentenceSeparator: ' ', textStyle: 'cased', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('ru'), name: 'Русский', sentenceSeparator: ' ', textStyle: 'cased', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('uk'), name: 'Українська', sentenceSeparator: ' ', textStyle: 'cased', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('ar'), name: 'العربية', sentenceSeparator: ' ', textStyle: 'tall-marks', dir: 'rtl', listSeparator: '، ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('fa'), name: 'فارسی', sentenceSeparator: ' ', textStyle: 'tall-marks', dir: 'rtl', listSeparator: '، ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('th'), name: 'ไทย', sentenceSeparator: ' ', textStyle: 'tall-marks', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('ja'), name: '日本語', sentenceSeparator: '', textStyle: 'no-caps', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('zh-Hans'), name: '简体中文', sentenceSeparator: '', textStyle: 'no-caps', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('zh-Hant'), name: '繁體中文', sentenceSeparator: '', textStyle: 'no-caps', dir: 'ltr', listSeparator: ', ' }),
  Object.freeze({ tag: /** @type {UiLocale} */ ('ko'), name: '한국어', sentenceSeparator: ' ', textStyle: 'no-caps', dir: 'ltr', listSeparator: ', ' })
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

// First Strong Isolate, Left-to-Right Isolate, Right-to-Left Isolate and Pop Directional Isolate (Unicode bidi).
const FSI = '\u2068';
const LRI = '\u2066';
const RLI = '\u2067';
const PDI = '\u2069';
// The direction marks Intl.ListFormat writes into a right-to-left list (fa); isolates take their place.
const DIRECTION_MARKS = /[\u200E\u200F]/g;
// A conjunction و that a list pattern writes against the next item (Arabic: `A وB`).
const GLUED_AND = /و$/u;
// A shown text that begins with an Arabic letter, after the isolates and marks in front of it.
const ARABIC_LETTER_START = /^[\u061C\u200E\u200F\u2066-\u2069]*(?=\p{L})\p{Script=Arabic}/u;
// A left-to-right run of a text: placeholders, Latin letters, digits and the number signs # % ‰ + - − /,
// with one . , or : between two of them (`#{n}`, `{value}%`, `{done}/{total}`, `{record}:{start}`), a
// range or a size between two of them (`{start}–{end}`, `{from}..{to}`, `{width} × {height}`), and the
// Latin words written beside them after one space (`{nt} skew`, `{n} bp`, `the record's`). Without
// isolates the bidi algorithm reads such a run with a left-to-right value as one run, so the run is
// isolated whole: an isolate of the value alone would show `5#` for `#5`, `skew GC` for `GC skew` and
// `9–1` for `1–9`. Two values with only a space between them (`{record} {selector}`) are two runs, so a
// right-to-left sentence keeps the order its translator wrote; a value whose unit the code writes is one
// placeholder (`{value}` = `12.5 px`, services/track-slot-display.js).
const LTR_CHAR = String.raw`[\p{Script=Latin}0-9#%‰+\-\u2212/]`;
const LTR_UNIT = String.raw`(?:\{[A-Za-z_][A-Za-z0-9_]*\}|${LTR_CHAR})`;
const LTR_WORD = String.raw`\p{Script=Latin}(?:[.,:']?${LTR_CHAR})*`;
const LTR_NEXT = String.raw`(?:(?:\.\.|[.,:\u2013]|[ \u00A0]?\u00D7[ \u00A0]?)?${LTR_UNIT})`;
const LTR_RUN = new RegExp(String.raw`(?:${LTR_WORD}[ \u00A0])*${LTR_UNIT}${LTR_NEXT}*(?:[ \u00A0]${LTR_WORD})*`, 'gu');
const ONE_PLACEHOLDER = /^\{[A-Za-z_][A-Za-z0-9_]*\}$/;

/**
 * A value as a text shows it: a number as `numberText` writes it, so every
 * number of a sentence reads alike (`1,234 of 4,321`) whichever call gives
 * it; a text value keeps its characters.
 * @param {string | number} value
 * @returns {string}
 */
const valueText = (value) => (typeof value === 'number' ? numberText(value) : value);

/**
 * Puts each value that `params` gives in place of its `{name}` (`valueText`);
 * any other placeholder stays. No isolates: the English text of an error
 * (utils/error-normalization.js), a stored text and every left-to-right
 * language.
 * @param {string} text
 * @param {TranslationParams | undefined} params
 * @returns {string}
 */
export const fillText = (text, params) => (params
  ? text.replace(PLACEHOLDER, (/** @type {string} */ whole, /** @type {string} */ name) => (Object.hasOwn(params, name) ? valueText(params[name]) : whole))
  : text);

/**
 * @param {string} run
 * @param {TranslationParams} params
 */
const holdsValue = (run, params) => [...run.matchAll(PLACEHOLDER)].some(([, name]) => Object.hasOwn(params, name));

/**
 * A left-to-right run of a right-to-left text with its values: a run that is
 * one value is that value's isolate (FSI ... PDI), so a name keeps the
 * direction of its own letters; a longer run is one left-to-right isolate
 * (LRI ... PDI) with each value isolated inside it, so `#2`, `1/3`, `1–9`,
 * `5 px` and `GC skew` read left to right whatever letters a value has.
 * @param {string} run
 * @param {TranslationParams} params
 */
const isolatedRun = (run, params) => {
  if (!holdsValue(run, params)) return run;
  if (ONE_PLACEHOLDER.test(run)) return `${FSI}${fillText(run, params)}${PDI}`;
  return `${LRI}${run.replace(PLACEHOLDER, (/** @type {string} */ whole, /** @type {string} */ name) => (Object.hasOwn(params, name) ? `${FSI}${valueText(params[name])}${PDI}` : whole))}${PDI}`;
};

/**
 * The shown text with its values (`fillText`). In a right-to-left language
 * each left-to-right run that holds a given value is isolated
 * (`isolatedRun`), so a file name, a number or a term keeps its own order
 * inside the sentence; a left-to-right language shows the values as they
 * are, so its texts keep their bytes. A text that is stored (a Session
 * value) is filled without isolates: `fillText(t(key), params)`.
 * @param {string} text
 * @param {TranslationParams | undefined} params
 */
const substitute = (text, params) => (params && shownLocale().dir === 'rtl'
  ? text.replace(LTR_RUN, (run) => isolatedRun(run, params))
  : fillText(text, params));

/**
 * The text of a catalog value that is read as one string. Plural forms are a
 * catalog error on any key but the `other` form of a `tn` pair (the catalog
 * test refuses them); they read as their `other` text.
 * @param {string | PluralForms} value
 * @returns {string}
 */
const textOf = (value) => (typeof value === 'string' ? value : value.other);

/**
 * The shown text of a catalog key, or `fallback` when the shown locale has none.
 * @param {string} key
 * @param {string} fallback
 */
const lookup = (key, fallback) => {
  // Read the locale first: a render that calls `t()` depends on it.
  if (locale.value === SOURCE_LOCALE || !Object.hasOwn(catalog, key)) return fallback;
  return textOf(catalog[key]);
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
 * `entries` lack the key it needs. The keys are the English texts of the two
 * forms. When the `other` form's value is `PluralForms`, the `one` form's
 * value is the text for exactly 1 (ICU `=1`), and every other count takes the
 * object's text of the category the rules select, else its `other`: Russian
 * selects `one` for 21 too, where a `one` text without `{count}` would be
 * wrong. When it is a string, the category picks the key (`one` for the
 * category `one`, else `other`).
 * @param {Intl.PluralRules} rules
 * @param {TranslationCatalog} entries
 * @param {number} count
 * @param {string} one
 * @param {string} other
 * @returns {string | null}
 */
export const pluralTranslation = (rules, entries, count, one, other) => {
  const forms = Object.hasOwn(entries, other) ? entries[other] : undefined;
  if (typeof forms === 'object' && count !== 1) return forms[rules.select(count)] ?? forms.other;
  const key = typeof forms === 'object' ? one : pluralForm(rules, count, one, other);
  return Object.hasOwn(entries, key) ? textOf(entries[key]) : null;
};

// The placeholder a slot keeps in place of a count (`I18nText` fills it with an element).
const COUNT_SLOT = '{count}';

/**
 * The shown text of a count. The shown locale's plural rules
 * (Intl.PluralRules) choose the form, and that form's English text is the
 * catalog key: `fr` and `pt` take `one` for 0 and 1, `en` and `es` for 1
 * only, `ja`, `zh-Hans`, `zh-Hant` and `ko` never. A value with plural
 * forms gives the text of the category (`pluralTranslation`). A missing
 * translation shows the form English chooses. `{count}` is always the count,
 * as `numberText` writes it (`1,234`); a `{count}` that `params` gives is
 * not shown, except the kept placeholder `{ count: '{count}' }` of a slot.
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
  return substitute(shown, { ...params, count: params?.count === COUNT_SLOT ? COUNT_SLOT : count });
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

// The full stops of the UI languages' scripts: . (Latin, Cyrillic, Korean)
// and 。 (Japanese, Chinese). Thai ends a sentence with none.
const FULL_STOP_END = /[.。]$/u;

/**
 * A shown sentence without its final full stop, for a template that ends the
 * sentence with its own mark (`t('{pair}: {text}.', { text: sentenceBody(x) })`):
 * the locale's mark then shows once, never 。。. Only a full stop goes;
 * ?, ! and … keep their meaning and stay.
 * @param {string} text a shown text
 * @returns {string}
 */
export const sentenceBody = (text) => text.replace(FULL_STOP_END, '');

/**
 * Sentences shown one after the other, each already a shown text.
 * @param {...string} sentences
 * @returns {string}
 */
export const joinSentences = (...sentences) => sentences.join(sentenceSeparator());

/**
 * The arrow between two values that means "then" or "to" where a template or
 * a module writes it outside a catalog text (`#1 → #2`, current → restored):
 * it points along the reading order, → in a left-to-right language and ← in
 * a right-to-left one, as the catalogs write that arrow. An arrow that names
 * a side of the figure is not this one.
 * @returns {string}
 */
export const thenArrow = () => (shownLocale().dir === 'rtl' ? '←' : '→');

/**
 * "a, then b" written outside a catalog text (`#1 → #2`, current →
 * restored, a path `5,000 → 9,000 → 1 → 200`): `a → b` in a left-to-right
 * language; in a right-to-left one `a ← b` as one right-to-left isolate
 * (RLI ... PDI) with each step isolated (FSI ... PDI), so it reads right to
 * left also inside a sentence's isolate or a list item, and each step keeps
 * its own order. A key never writes this arrow between two placeholders. Use
 * `thenArrow()` only where the arrow stands between two elements rendered
 * apart.
 * @param {string} first
 * @param {...string} then the steps after `first`, at least one
 * @returns {string}
 */
export const thenText = (first, ...then) => (shownLocale().dir === 'rtl'
  ? `${RLI}${[first, ...then].map((step) => `${FSI}${step}${PDI}`).join(` ${thenArrow()} `)}${PDI}`
  : [first, ...then].join(` ${thenArrow()} `));

/**
 * A part of a list: an item, or the text the list pattern writes between items.
 * @typedef {{ type: 'element' | 'literal', value: string }} ListPart
 */

/**
 * A list in a right-to-left language: one right-to-left isolate (RLI ... PDI)
 * with each item an isolate of its own (FSI ... PDI), so the list reads right
 * to left also inside a sentence's isolate and an item keeps its own order.
 * The direction marks a pattern may write go. A conjunction و that the pattern
 * writes against the next item (Arabic) stays there before an item that
 * begins with an Arabic letter and is set apart with a space before any other
 * (a digit, a Latin name): I18N_STYLE 21.4 (`وملف`, `و GenBank`, `و 5`).
 * No parts show as the empty text.
 * @param {readonly ListPart[]} parts
 * @returns {string}
 */
const rightToLeftList = (parts) => (parts.length ? `${RLI}${parts.map(({ type, value }, index) => {
  if (type === 'element') return `${FSI}${value}${PDI}`;
  const literal = value.replace(DIRECTION_MARKS, '');
  const next = parts[index + 1];
  return next && GLUED_AND.test(literal) && !ARABIC_LETTER_START.test(next.value) ? `${literal} ` : literal;
}).join('')}${PDI}` : '');

/**
 * Shown texts as one list: "A", "A and B", "A, B, and C" in English, the
 * locale's own list pattern (Intl.ListFormat) in the other locales, read
 * right to left in a right-to-left language (`rightToLeftList`).
 * @param {readonly string[]} items
 * @returns {string}
 */
export const listText = (items) => {
  const { tag, dir } = shownLocale();
  if (dir === 'rtl') return rightToLeftList(new Intl.ListFormat(tag, { type: 'conjunction' }).formatToParts(items));
  if (tag !== SOURCE_LOCALE) return new Intl.ListFormat(tag, { type: 'conjunction' }).format(items);
  return englishListText(items);
};

/**
 * A list that ends in "and" as English writes it ("A", "A and B", "A, B, and
 * C"): `listText` in English, and the English text of an error in every UI
 * language (utils/error-normalization.js).
 * @param {readonly string[]} items
 * @returns {string}
 */
export const englishListText = (items) => (items.length < 3
  ? items.join(' and ')
  : `${items.slice(0, -1).join(', ')}, and ${items.at(-1)}`);

/**
 * Shown items one after the other with no conjunction (record ids, gene
 * names, file names, ranges, translated names), separated by the shown
 * language's `listSeparator`: `A, B, C` in English and every left-to-right
 * language (the bytes of `items.join(', ')`), `A، B، C` in Persian and
 * Arabic, read right to left there (`rightToLeftList`). No items show as the
 * empty text. A list that ends in "and" is `listText`.
 * @param {readonly string[]} items each already a shown text
 * @returns {string}
 */
export const itemListText = (items) => {
  const { dir, listSeparator } = shownLocale();
  if (dir !== 'rtl') return items.join(listSeparator);
  return rightToLeftList(items.flatMap((value, index) => {
    /** @type {ListPart} */
    const element = { type: 'element', value };
    return index ? [{ type: 'literal', value: listSeparator }, element] : [element];
  }));
};

// Every number of the UI and of an exported interactive file reads alike in every
// language: Western digits, the English decimal point and grouping by thousands.
const NUMBER_TEXT = new Intl.NumberFormat('en-US', { maximumFractionDigits: 20 });
/** @type {Map<number, Intl.NumberFormat>} */
const fixedNumberTexts = new Map();

/**
 * A number as the UI writes it in every language (PD-OI-106: "Numbers keep the
 * English decimal point and digit grouping in every language"): `4,641,652`,
 * `0.25`, never the browser's or the shown locale's digits and grouping. Every
 * fraction digit is kept, or exactly `fractionDigits` of them. The one owner of
 * number text; the exported interactive file's viewer keeps a copy
 * (`numberText` in services/standalone-interactivity-assets.js), pinned to this
 * one by tests/web/feature-search-scope.test.mjs.
 * @param {number} value
 * @param {number} [fractionDigits]
 * @returns {string}
 */
export const numberText = (value, fractionDigits) => {
  if (fractionDigits === undefined) return NUMBER_TEXT.format(value);
  let format = fixedNumberTexts.get(fractionDigits);
  if (!format) {
    format = new Intl.NumberFormat('en-US', { minimumFractionDigits: fractionDigits, maximumFractionDigits: fractionDigits });
    fixedNumberTexts.set(fractionDigits, format);
  }
  return format.format(value);
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
  let at = 0;
  /** @param {number} end */
  const pushText = (end) => {
    const shown = substitute(text.slice(at, end), params);
    if (shown) parts.push({ text: shown });
  };
  for (const match of text.matchAll(PLACEHOLDER)) {
    const [whole, name] = match;
    if (!slots.includes(name)) continue;
    pushText(match.index ?? 0);
    parts.push({ slot: name });
    at = (match.index ?? 0) + whole.length;
  }
  pushText(text.length);
  return parts;
};

/** @type {CatalogLoader} */
const importCatalog = (tag) => import(`./i18n-catalogs/${tag}.js`);

// The latest `setUiLocale` call; an earlier call whose catalog arrives later is dropped.
let latestRequest = 0;

/**
 * Shows the UI in `next`: loads its catalog, then switches, and sets the
 * page's `lang`, `dir` and `data-ui-text-style` (the one writer of all three). `en` needs no catalog. When
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
  const root = globalThis.document?.documentElement;
  root?.setAttribute('lang', tag);
  root?.setAttribute('dir', shownLocale().dir);
  root?.setAttribute('data-ui-text-style', shownLocale().textStyle);
  return true;
};
