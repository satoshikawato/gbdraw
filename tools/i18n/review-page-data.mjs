#!/usr/bin/env node
// @ts-check
// Writes the data files of the translation review page (branch
// translation-review, served by GitHub Pages) from this checkout's sources:
//
//   node tools/i18n/review-page-data.mjs --out <translation-review>/data [--shots <shots.json>]
//
// locales.json (the languages with a catalog, {tag, label}, in menu order),
// keys.json (each key with its id, context, first call site and glossary
// hints; no line numbers, so the data changes only when text changes),
// draft-<locale>.json (the catalog by key id, each value one string:
// `reviewDraftText`), notes-<locale>.json
// (empty: drafting notes have no owner in the repository), glossary.json
// (I18N_GLOSSARY.tsv), style.md (I18N_STYLE.md), and shots.json: the
// screenshot map, read from --shots or from the old <out>/shots.json, with
// the ids of keys the sources no longer have left out. The output folder is
// replaced whole. A key id is the first 12 hex digits of sha1(key); the page
// stores review rows under it (db/current/rev-<locale>/c<first digit>.json).
import { createHash } from 'node:crypto';
import { existsSync, mkdirSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

import { UI_LOCALES } from '../../gbdraw/web/js/utils/i18n.js';
import { CATALOG_LOCALES, glossaryMatcher, pluralFormEntries, readCatalog, readGlossary } from './catalogs.mjs';
import { REPOSITORY_ROOT, sourceKeys } from './source-keys.mjs';

/** @import { PluralForms } from '../../gbdraw/web/js/utils/i18n.js' */
/** @typedef {import('./source-keys.mjs').KeyRole} KeyRole */
/**
 * One string of the review page (data/keys.json).
 * @typedef {object} ReviewKey
 * @property {string} id sha1(key)[:12]
 * @property {string} en the English text
 * @property {string} sense the tx() context, or ''
 * @property {string} key the catalog key
 * @property {string} area the part of the app of the first call site
 * @property {string} where the file of the first call site
 * @property {KeyRole} kind the role of the first call site
 * @property {number} uses the number of call sites
 * @property {string} ctx the source line of the first call site
 * @property {number[]} g glossary.json indexes of the terms in the English
 */
/**
 * The screenshot map: screenshots, and per key id the rectangles [shot index, x, y, w, h].
 * @typedef {{ scale: number, shots: { file: string, w: number, h: number, state: string }[], keys: Record<string, number[][]> }} ShotMap
 */

/** @param {string} key */
export const reviewKeyId = (key) => createHash('sha1').update(key, 'utf8').digest('hex').slice(0, 12);

/**
 * The text the review page shows for a catalog value: the string, or the
 * plural forms joined in CLDR category order (`few: … | many: … | other: …`).
 * @param {string | PluralForms} value
 * @returns {string}
 */
export const reviewDraftText = (value) => (typeof value === 'string' ? value
  : pluralFormEntries(value).map(([category, text]) => `${category}: ${text}`).join(' | '));

/** @param {unknown} value */
const json = (value) => JSON.stringify(value, null, 0);

/**
 * The files of data/, by name.
 * @param {{ shots: ShotMap | null }} input
 * @returns {Promise<Map<string, string>>}
 */
export const reviewPageData = async ({ shots }) => {
  const { keys, problems } = await sourceKeys();
  if (problems.length) throw new Error(`the sources have calls whose key cannot be read:\n${problems.join('\n')}`);
  const glossary = readGlossary();
  const termsOf = glossaryMatcher(glossary);
  /** @type {ReviewKey[]} */
  const rows = [...keys.values()].map((entry) => {
    const [site] = entry.sites;
    return {
      id: reviewKeyId(entry.key), en: entry.english, sense: entry.context ?? '', key: entry.key, area: site.area,
      where: site.file, kind: site.role, uses: entry.sites.length, ctx: site.snippet, g: termsOf(entry.english)
    };
  });
  if (new Set(rows.map(({ id }) => id)).size !== rows.length) throw new Error('two keys share an id');
  /** @type {Map<string, string>} */
  const files = new Map([
    ['locales.json', json(CATALOG_LOCALES.map((tag) => ({ tag, label: UI_LOCALES.find((option) => option.tag === tag)?.name ?? tag })))],
    ['keys.json', json(rows)]
  ]);
  for (const locale of CATALOG_LOCALES) {
    const { entries } = readCatalog(locale);
    files.set(`draft-${locale}.json`, json(Object.fromEntries(rows.filter(({ key }) => Object.hasOwn(entries, key)).map(({ id, key }) => [id, reviewDraftText(entries[key])]))));
    files.set(`notes-${locale}.json`, '{}');
  }
  files.set('glossary.json', json(glossary.map((row) => Object.fromEntries(['en', ...CATALOG_LOCALES, 'note'].map((column) => [column, row[column] ?? ''])))));
  files.set('style.md', readFileSync(join(REPOSITORY_ROOT, 'docs', 'internal', 'I18N_STYLE.md'), 'utf8'));
  if (shots) {
    const ids = new Set(rows.map(({ id }) => id));
    files.set('shots.json', json({ ...shots, keys: Object.fromEntries(Object.entries(shots.keys).filter(([id]) => ids.has(id))) }));
  }
  return files;
};

/** @param {string[]} argv */
const main = async (argv) => {
  const option = (/** @type {string} */ name) => {
    const at = argv.indexOf(name);
    return at >= 0 ? argv[at + 1] : undefined;
  };
  const out = option('--out');
  if (!out) {
    process.stderr.write('usage: review-page-data.mjs --out <data folder> [--shots <shots.json>]\n');
    return 2;
  }
  const shotsPath = option('--shots') ?? join(out, 'shots.json');
  /** @type {ShotMap | null} */
  const shots = existsSync(shotsPath) ? JSON.parse(readFileSync(shotsPath, 'utf8')) : null;
  const files = await reviewPageData({ shots });
  rmSync(out, { recursive: true, force: true });
  mkdirSync(out, { recursive: true });
  for (const [name, text] of files) writeFileSync(join(out, name), text);
  process.stdout.write(`${out}: ${files.size} files, ${JSON.parse(String(files.get('keys.json'))).length} keys\n`);
  return 0;
};

if (import.meta.url === pathToFileURL(process.argv[1] ?? '').href) process.exitCode = await main(process.argv.slice(2));
