// @ts-check
// The locale catalogs (gbdraw/web/js/utils/i18n-catalogs/<locale>.js) as data,
// and the glossary (docs/internal/I18N_GLOSSARY.tsv). A catalog file is one
// canonical text: its header comment, then `export default Object.freeze({`
// with one `"key": value` row per key (keys in code-unit order; a value is a
// JSON string, or the plural forms of the `other` form of a `tn()` pair as
// one JSON object in CLDR category order), then `});`. Reading parses that text; writing prints it, so a
// catalog read and written unchanged keeps every byte.
import { readFileSync } from 'node:fs';
import { join } from 'node:path';

import { PLURAL_CATEGORIES, SOURCE_LOCALE, UI_LOCALES } from '../../gbdraw/web/js/utils/i18n.js';
import { REPOSITORY_ROOT, placeholdersOf } from './source-keys.mjs';

/** @import { PluralForms, TranslationCatalog } from '../../gbdraw/web/js/utils/i18n.js' */
/** @typedef {import('./source-keys.mjs').SourceKey} SourceKey */
/**
 * What the catalog check reads of a source key.
 * @typedef {Pick<SourceKey, 'key' | 'english' | 'plural'>} CatalogKey
 */

export const CATALOG_ROOT = join(REPOSITORY_ROOT, 'gbdraw', 'web', 'js', 'utils', 'i18n-catalogs');
export const GLOSSARY_PATH = join(REPOSITORY_ROOT, 'docs', 'internal', 'I18N_GLOSSARY.tsv');
/**
 * The locales with a catalog: the UI languages except English, in menu order.
 * @type {readonly string[]}
 */
export const CATALOG_LOCALES = Object.freeze(UI_LOCALES.map(({ tag }) => String(tag)).filter((tag) => tag !== SOURCE_LOCALE));

const BODY_START = 'export default Object.freeze({\n';
const BODY_END = '\n});\n';

/**
 * A catalog file as data.
 * @typedef {object} CatalogFile
 * @property {string} header the text before `export default`
 * @property {Record<string, string | PluralForms>} entries
 */

/**
 * Whether a value is an object (not an array, not null): the shape of plural forms.
 * @param {unknown} value
 * @returns {value is Readonly<Record<string, unknown>>}
 */
export const isRecord = (value) => typeof value === 'object' && value !== null && !Array.isArray(value);

/**
 * The plural forms of a value in CLDR category order, as [category, text] pairs.
 * @param {PluralForms} forms
 * @returns {[Intl.LDMLPluralRule, string][]}
 */
export const pluralFormEntries = (forms) => PLURAL_CATEGORIES.flatMap((category) => {
  const text = forms[category];
  return text === undefined ? [] : [[category, text]];
});

/**
 * The canonical text of a catalog value: a JSON string, or the plural forms
 * as one JSON object in CLDR category order.
 * @param {string | PluralForms} value
 * @returns {string}
 */
export const formatCatalogValue = (value) => (typeof value === 'string' ? JSON.stringify(value)
  : `{${pluralFormEntries(value).map(([category, text]) => `${JSON.stringify(category)}: ${JSON.stringify(text)}`).join(', ')}}`);

/**
 * @param {string} text the catalog file
 * @returns {CatalogFile}
 */
export const parseCatalog = (text) => {
  const start = text.indexOf(BODY_START);
  if (start < 0 || !text.endsWith(BODY_END)) throw new Error('not a catalog in the canonical format');
  /** @type {Record<string, unknown>} */
  const entries = JSON.parse(`{${text.slice(start + BODY_START.length, text.length - BODY_END.length)}}`);
  for (const [key, value] of Object.entries(entries)) {
    if (typeof value !== 'string' && !(isRecord(value) && Object.values(value).every((form) => typeof form === 'string'))) {
      throw new Error(`${JSON.stringify(key)}: the value is not a string or an object of strings`);
    }
  }
  // Each value is a string or an object of strings; catalogProblems checks the categories.
  return { header: text.slice(0, start), entries: /** @type {Record<string, string | PluralForms>} */ (entries) };
};

/**
 * @param {CatalogFile} catalog
 * @returns {string} the canonical text
 */
export const formatCatalog = ({ header, entries }) => {
  const rows = Object.keys(entries).sort().map((key) => `  ${JSON.stringify(key)}: ${formatCatalogValue(entries[key])}`);
  return `${header}${BODY_START}${rows.join(',\n')}${BODY_END}`;
};

/** @param {string} locale @param {string} [root] */
export const catalogPath = (locale, root = CATALOG_ROOT) => {
  if (!CATALOG_LOCALES.includes(locale)) throw new Error(`no catalog for locale ${JSON.stringify(locale)}`);
  return join(root, `${locale}.js`);
};

/** @param {string} locale @param {string} [root] @returns {CatalogFile} */
export const readCatalog = (locale, root = CATALOG_ROOT) => parseCatalog(readFileSync(catalogPath(locale, root), 'utf8'));

/**
 * The problem of one text against its English text, or null.
 * @param {string} where the key, quoted, and the plural category of a form
 * @param {string} english
 * @param {unknown} value
 */
const textProblem = (where, english, value) => {
  if (typeof value !== 'string') return `${where}: the value is not a string`;
  if (!value.trim()) return `${where}: the translation is empty`;
  const expected = placeholdersOf(english);
  const actual = placeholdersOf(value);
  return expected.join() === actual.join() ? null
    : `${where}: placeholders {${actual.join('}, {')}} differ from the English {${expected.join('}, {')}}`;
};

/**
 * The problems of one translation against its source key: a value that is
 * not a string, empty, or with placeholders that differ from the English.
 * The value of the `other` form of a `tn()` pair (`plural.other`) may
 * instead be plural forms (`PluralForms`): an object of CLDR categories with
 * `other`, each form checked like a string. The `one` form's value is the
 * text for exactly 1 and stays a string.
 * @param {CatalogKey} source
 * @param {unknown} value
 * @returns {string[]}
 */
export const translationProblems = ({ key, english, plural }, value) => {
  const where = JSON.stringify(key);
  if (!isRecord(value)) {
    const problem = textProblem(where, english, value);
    return problem ? [problem] : [];
  }
  if (!plural || plural.other !== key) return [`${where}: plural forms are allowed only for the other form of a tn() pair`];
  /** @type {readonly string[]} */
  const categories = PLURAL_CATEGORIES;
  const unknown = Object.keys(value).filter((category) => !categories.includes(category)).sort();
  if (unknown.length) return [`${where}: ${unknown.join(', ')} are not plural categories`];
  if (!Object.hasOwn(value, 'other')) return [`${where}: the plural forms lack "other"`];
  return Object.entries(value).flatMap(([category, text]) => {
    const problem = textProblem(`${where} ${category}`, english, text);
    return problem ? [problem] : [];
  });
};

/**
 * The problems of one catalog against the keys of the sources
 * (`translationProblems`), on a key the sources have. A key the sources no
 * longer have is not a problem here (it is stale until the translation batch
 * drops it; `catalogStatus` lists it and `sync.mjs status --require-complete`
 * fails on it).
 * @param {ReadonlyMap<string, CatalogKey>} keys the keys of the sources
 * @param {Readonly<Record<string, unknown>>} catalog
 * @returns {string[]}
 */
export const catalogProblems = (keys, catalog) => Object.entries(catalog).flatMap(([key, value]) => {
  const source = keys.get(key);
  return source ? translationProblems(source, value) : [];
});

/**
 * The keys of the sources a catalog translates, lacks, and holds without a source.
 * @param {Iterable<string>} keys the keys of the sources
 * @param {TranslationCatalog} catalog
 * @returns {{ translated: string[], missing: string[], stale: string[] }}
 */
export const catalogStatus = (keys, catalog) => {
  const wanted = new Set(keys);
  const all = [...wanted];
  return {
    translated: all.filter((key) => Object.hasOwn(catalog, key)),
    missing: all.filter((key) => !Object.hasOwn(catalog, key)),
    stale: Object.keys(catalog).filter((key) => !wanted.has(key))
  };
};

/**
 * One row of the glossary: the English term, a cell per catalog locale, and the note.
 * @typedef {Record<string, string> & { en: string, note: string }} GlossaryRow
 */

/**
 * @param {string} text the TSV: a header row `en`, the locales, `note`
 * @returns {GlossaryRow[]}
 */
export const parseGlossary = (text) => {
  const [header, ...rows] = text.split('\n').filter((line) => line !== '');
  const columns = header.split('\t');
  if (columns[0] !== 'en' || columns.at(-1) !== 'note') throw new Error(`glossary header ${header}`);
  return rows.map((row) => {
    const cells = row.split('\t');
    return /** @type {GlossaryRow} */ (Object.fromEntries(columns.map((column, at) => [column, cells[at] ?? ''])));
  });
};

/** @returns {GlossaryRow[]} */
export const readGlossary = () => parseGlossary(readFileSync(GLOSSARY_PATH, 'utf8'));

const UNITS = new Set(['bp', 'kb', 'Mb', 'pt', 'px', '×R']);

/**
 * The pattern that finds a glossary term in an English text: the term without
 * a trailing parenthesis, case-insensitive, whole words, an optional plural s.
 * @param {string} en
 * @returns {RegExp | null}
 */
const termPattern = (en) => {
  const base = en.replace(/\s*\([^)]*\)\s*$/, '').trim();
  if (base.length < 2 || UNITS.has(base)) return null;
  return new RegExp(`(?<![A-Za-z])${base.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}s?(?![A-Za-z])`, 'i');
};

/**
 * The glossary rows whose English term occurs in a text, longest term first,
 * at most six; a row kept as written in every locale (a name) is no hint.
 * @param {readonly GlossaryRow[]} glossary
 * @returns {(english: string) => number[]} the row indexes for a text
 */
export const glossaryMatcher = (glossary) => {
  const patterns = glossary.flatMap((row, index) => {
    if (CATALOG_LOCALES.every((locale) => !row[locale] || row[locale] === row.en)) return [];
    const pattern = termPattern(row.en);
    return pattern ? [{ index, pattern, length: row.en.length }] : [];
  });
  return (english) => patterns.filter(({ pattern }) => pattern.test(english))
    .sort((a, b) => b.length - a.length).slice(0, 6).map(({ index }) => index);
};
