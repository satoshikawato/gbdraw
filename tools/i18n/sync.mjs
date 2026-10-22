#!/usr/bin/env node
// @ts-check
// Keeps the Web UI catalogs in step with the keys of the sources
// (.claude/skills/translate-web-ui/SKILL.md runs it before each release).
//
//   node tools/i18n/sync.mjs status [--locale <x>] [--require-complete]
//   node tools/i18n/sync.mjs export --locale <x> [--all] > work.jsonl
//   node tools/i18n/sync.mjs apply --locale <x> <work.jsonl>
//
// status: per locale the translated, missing and stale keys; with
//   --require-complete it exits 1 when a catalog lacks a key or holds one the
//   sources no longer have (the release tier runs it).
// export: one JSON line per missing key (every key with --all) for a drafter:
//   the English, its placeholders, the tn() pair, the tx() context, the call
//   sites with their roles, the glossary rows of its terms with this locale's
//   cell, the current translation, and, for a missing key whose English is a
//   small edit of a stale key, that key and its translation as `previous`.
// apply: reads lines of {"key", "translation"}, refuses a line that is not an
//   object with those two strings, a key the sources do not have, and an empty
//   translation or one whose placeholders differ, then writes the
//   catalog in its canonical format without its stale keys. Nothing is written
//   when a line is refused; a catalog that does not change is not rewritten.
import { readFileSync, writeFileSync } from 'node:fs';
import { pathToFileURL } from 'node:url';

import {
  CATALOG_LOCALES, CATALOG_ROOT, catalogPath, catalogStatus, formatCatalog, glossaryMatcher, readCatalog, readGlossary, translationProblem
} from './catalogs.mjs';
import { placeholdersOf, sourceKeys } from './source-keys.mjs';

/** @typedef {import('./source-keys.mjs').SourceKey} SourceKey */
/** @typedef {import('./source-keys.mjs').KeyRole} KeyRole */

/**
 * One line of an export file.
 * @typedef {object} ExportRecord
 * @property {string} key the catalog key
 * @property {string} english
 * @property {string | null} context the tx() context
 * @property {string[]} placeholders
 * @property {{ one: string, other: string } | null} plural the tn() pair this key is a form of
 * @property {{ file: string, line: number, role: KeyRole }[]} sites
 * @property {{ en: string, term: string, note: string }[]} glossary
 * @property {string | null} translation the catalog's translation, or null
 * @property {{ key: string, translation: string, status: 'needs-check' } | null} previous
 */
/**
 * One line of an apply file; other fields of an export line are ignored.
 * @typedef {{ key: string, translation: string }} ApplyRecord
 */
/**
 * @typedef {object} SyncOptions
 * @property {string} [catalogRoot] the folder of the catalogs (tests pass a copy)
 * @property {(line: string) => void} [print]
 */

const USAGE = 'usage: sync.mjs status [--locale <x>] [--require-complete] | export --locale <x> [--all] | apply --locale <x> <file.jsonl>';

/** Levenshtein distance, for the `previous` hint. @param {string} a @param {string} b */
const distance = (a, b) => {
  let row = Array.from({ length: b.length + 1 }, (_, at) => at);
  for (let i = 1; i <= a.length; i += 1) {
    const next = [i];
    for (let j = 1; j <= b.length; j += 1) next.push(Math.min(row[j] + 1, next[j - 1] + 1, row[j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1)));
    row = next;
  }
  return row[b.length];
};

/**
 * The stale key whose text a missing key's text is a small edit of (at most
 * a third of the longer text changed), or null.
 * @param {string} key
 * @param {readonly string[]} stale
 */
export const previousKeyOf = (key, stale) => {
  let best = null;
  let bestScore = 1 / 3;
  for (const candidate of stale) {
    const longer = Math.max(key.length, candidate.length);
    if (Math.abs(key.length - candidate.length) > longer * bestScore) continue;
    const score = distance(key, candidate) / longer;
    if (score <= bestScore) [best, bestScore] = [candidate, score];
  }
  return best;
};

/**
 * @param {unknown} value one parsed line of an apply file
 * @returns {ApplyRecord | null} null unless a plain object with a string key and a string translation
 */
const applyRecordOf = (value) => {
  if (typeof value !== 'object' || value === null || Array.isArray(value)) return null;
  const { key, translation } = /** @type {Record<string, unknown>} */ (value);
  return typeof key === 'string' && typeof translation === 'string' ? { key, translation } : null;
};

/**
 * @param {string[]} argv
 * @returns {{ command: string, locale: string | null, flags: Set<string>, files: string[] }}
 */
const parseArgs = (argv) => {
  const [command = '', ...rest] = argv;
  let locale = null;
  const flags = new Set();
  const files = [];
  for (let at = 0; at < rest.length; at += 1) {
    if (rest[at] === '--locale') locale = rest[++at] ?? '';
    else if (rest[at].startsWith('--')) flags.add(rest[at]);
    else files.push(rest[at]);
  }
  return { command, locale, flags, files };
};

/**
 * Runs one command; returns the exit status.
 * @param {string[]} argv
 * @param {SyncOptions} [options]
 * @returns {Promise<number>}
 */
export const runSync = async (argv, { catalogRoot = CATALOG_ROOT, print = (line) => process.stdout.write(`${line}\n`) } = {}) => {
  const { command, locale, flags, files } = parseArgs(argv);
  /** @type {Record<string, string[]>} */
  const commands = { status: ['--require-complete'], export: ['--all'], apply: [] };
  const known = Object.hasOwn(commands, command) ? commands[command] : null;
  if (!known || [...flags].some((flag) => !known.includes(flag)) || (command !== 'status' && !locale)
      || (command === 'apply') !== (files.length === 1) || (command !== 'apply' && files.length)) {
    print(USAGE);
    return 2;
  }
  if (locale !== null && !CATALOG_LOCALES.includes(locale)) {
    print(`unknown locale ${JSON.stringify(locale)}; catalogs: ${CATALOG_LOCALES.join(', ')}`);
    return 2;
  }
  const { keys, problems } = await sourceKeys();
  if (command === 'status') {
    for (const problem of problems) print(`source: ${problem}`);
    let incomplete = problems.length > 0;
    for (const tag of locale ? [locale] : CATALOG_LOCALES) {
      const { translated, missing, stale } = catalogStatus(keys.keys(), readCatalog(tag, catalogRoot).entries);
      print(`${tag}: ${translated.length} of ${keys.size} translated, ${missing.length} missing, ${stale.length} stale`);
      for (const key of missing) print(`  missing ${JSON.stringify(key)}`);
      for (const key of stale) print(`  stale ${JSON.stringify(key)}`);
      incomplete ||= missing.length > 0 || stale.length > 0;
    }
    if (incomplete && flags.has('--require-complete')) {
      print('Every catalog needs every key of the sources and no other key: translate them with .claude/skills/translate-web-ui/SKILL.md.');
      return 1;
    }
    return 0;
  }
  const tag = /** @type {string} */ (locale);
  const catalog = readCatalog(tag, catalogRoot);
  if (command === 'export') {
    const glossary = readGlossary();
    const termsOf = glossaryMatcher(glossary);
    const { stale } = catalogStatus(keys.keys(), catalog.entries);
    for (const entry of keys.values()) {
      const translation = Object.hasOwn(catalog.entries, entry.key) ? catalog.entries[entry.key] : null;
      if (translation !== null && !flags.has('--all')) continue;
      const previous = translation === null ? previousKeyOf(entry.key, stale) : null;
      /** @type {ExportRecord} */
      const record = {
        key: entry.key,
        english: entry.english,
        context: entry.context,
        placeholders: placeholdersOf(entry.english),
        plural: entry.plural,
        sites: entry.sites.map(({ file, line, role }) => ({ file, line, role })),
        glossary: termsOf(entry.english).map((index) => ({ en: glossary[index].en, term: glossary[index][tag] ?? '', note: glossary[index].note })),
        translation,
        previous: previous === null ? null : { key: previous, translation: catalog.entries[previous], status: 'needs-check' }
      };
      print(JSON.stringify(record));
    }
    return 0;
  }
  // apply
  /** @type {Record<string, string>} */
  const applied = {};
  const refused = [];
  const lines = readFileSync(files[0], 'utf8').split('\n');
  for (let at = 0; at < lines.length; at += 1) {
    if (!lines[at].trim()) continue;
    const where = `${files[0]}:${at + 1}`;
    /** @type {unknown} */
    let parsed;
    try {
      parsed = JSON.parse(lines[at]);
    } catch {
      refused.push(`${where}: not a JSON line`);
      continue;
    }
    const record = applyRecordOf(parsed);
    if (!record) {
      refused.push(`${where}: not an object with a string "key" and a string "translation"`);
      continue;
    }
    const entry = keys.get(record.key);
    if (!entry) {
      refused.push(`${where}: ${JSON.stringify(record.key)} is not a key of the sources`);
      continue;
    }
    if (Object.hasOwn(applied, entry.key)) {
      refused.push(`${where}: ${JSON.stringify(entry.key)} is given twice`);
      continue;
    }
    const problem = translationProblem(entry.key, entry.english, record.translation);
    if (problem) refused.push(`${where}: ${problem}`);
    else applied[entry.key] = record.translation;
  }
  if (refused.length) {
    for (const line of refused) print(`refused ${line}`);
    print(`${tag}: nothing written; fix the ${refused.length} refused line(s)`);
    return 1;
  }
  const { stale } = catalogStatus(keys.keys(), catalog.entries);
  const entries = { ...catalog.entries, ...applied };
  for (const key of stale) delete entries[key];
  const path = catalogPath(tag, catalogRoot);
  const before = readFileSync(path, 'utf8');
  const after = formatCatalog({ header: catalog.header, entries });
  const changed = Object.keys(applied).filter((key) => catalog.entries[key] !== applied[key]).length;
  if (after !== before) writeFileSync(path, after);
  print(`${tag}: ${Object.keys(applied).length} lines applied (${changed} changed), ${stale.length} stale keys dropped`);
  return 0;
};

if (import.meta.url === pathToFileURL(process.argv[1] ?? '').href) process.exitCode = await runSync(process.argv.slice(2));
