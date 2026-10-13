// @ts-check

import { SESSION_TABLE_LABELS } from './error-normalization.js';

/**
 * What Load did with a value it could not restore as stored.
 * @typedef {'dropped' | 'defaulted' | 'changed'} LoadRepairEffect
 */
/**
 * The rule of the producer that could not restore a value: one member per rule.
 * @typedef {'feature-edit-unmatched'
 *   | 'feature-visibility-narrowed'
 *   | 'hash-selector-unmapped'
 *   | 'table-row-joined'
 *   | 'table-row-dropped'
 *   | 'table-row-invalid'
 *   | 'legend-position-defaulted'
 *   | 'legend-side-defaulted'
 *   | 'legend-color-dropped'
 *   | 'legend-entry-duplicate'
 *   | 'legend-edit-overwritten'
 *   | 'legend-edit-overwritten-session46'
 *   | 'legend-edit-rekey-dropped'
 *   | 'record-display-drafts-dropped'
 * } LoadRepairCode
 */
/**
 * One thing Load could not restore as stored. The producer owns the meaning
 * (a code and bounded fields); this module owns the wording. `stored` is shown
 * in the dialog only: it never reaches the error log, a diagnostic, or the console.
 * @typedef {object} LoadRepair
 * @property {LoadRepairCode} code
 * @property {string} path Where in the file: dotted, with [i]; at most 200 characters.
 * @property {LoadRepairEffect} effect
 * @property {string} [stored] A bounded excerpt (at most 60 characters) of the stored value; never free text.
 * @property {string} [used] The value used instead.
 * @property {number} [count] When one record stands for several values.
 * @property {string} [label] The name the app shows for the place (a Legend caption), bounded.
 * @property {number} [line] The line of a table.
 */
/**
 * A value read from a file, with what its reader could not restore.
 * @template T
 * @typedef {{ value: T, repairs: readonly LoadRepair[] }} Repaired
 */
/**
 * The confirmation Load shows before it applies a file with repairs.
 * @typedef {object} LoadRepairDialog
 * @property {string} title
 * @property {string} description
 * @property {string[]} items
 * @property {string} note
 */
/**
 * The Load port of the caller that loads: it shows the dialog and resolves to
 * true for [Load], false for [Cancel].
 * @typedef {(dialog: LoadRepairDialog) => Promise<boolean>} ConfirmLoadRepairs
 */

const PATH_LIMIT = 200;
const EXCERPT_LIMIT = 60;
const LINE_LIMIT = 12;

/**
 * Codes whose values had no effect on the figure: Load drops them without a
 * dialog row (Owner answer Q1 of the Load repair receipt). A Session 46 drew
 * the later of two Legend edits of one row.
 * @type {ReadonlySet<LoadRepairCode>}
 */
const UNLISTED_CODES = new Set(/** @type {LoadRepairCode[]} */ (['legend-edit-overwritten-session46']));

/** @param {string} value @param {number} limit */
const bounded = (value, limit) => (value.length > limit ? `${value.slice(0, limit - 1)}…` : value);

/**
 * A repair record with its text fields bounded.
 * @param {LoadRepair} repair
 * @returns {LoadRepair}
 */
export const loadRepair = (repair) => ({
  ...repair,
  path: bounded(repair.path, PATH_LIMIT),
  ...(repair.stored === undefined ? {} : { stored: bounded(repair.stored, EXCERPT_LIMIT) }),
  ...(repair.used === undefined ? {} : { used: bounded(repair.used, EXCERPT_LIMIT) }),
  ...(repair.label === undefined ? {} : { label: bounded(repair.label, EXCERPT_LIMIT) })
});

/**
 * @param {string} singular
 * @param {string} plural
 * @param {readonly unknown[]} items
 */
const cappedList = (singular, plural, items) => (items.length === 1 ? `${singular} ${items[0]}`
  : `${plural} ${items.slice(0, 5).join(', ')}${items.length > 5 ? `, and ${items.length - 5} more` : ''}`);
/** @param {number} count */
const wasWere = (count) => (count === 1 ? 'was' : 'were');
/** @param {readonly LoadRepair[]} repairs */
const total = (repairs) => repairs.reduce((sum, repair) => sum + (repair.count ?? 1), 0);
/** @type {Readonly<Record<string, string>>} */
const LEGEND_LIST_MARKS = { entries: '', dormantEntries: ' (not drawn)', deletedEntries: ' (deleted)' };
/** @param {LoadRepair} repair The list is the last part of the path. */
const legendListMark = (repair) => LEGEND_LIST_MARKS[repair.path.split('.').pop() || ''] ?? '';

/** @type {Readonly<Record<string, string>>} */
const LEGEND_RECORD_KINDS = { colorOverrides: 'color edit', strokeOverrides: 'stroke edit' };
/** @param {LoadRepair} repair The collection is the last part of the path. */
const legendRecordKind = (repair) => LEGEND_RECORD_KINDS[repair.path.split('.').pop() || ''] ?? 'entry';
/** @param {readonly LoadRepair[]} repairs */
const overwrittenLegendEdits = (repairs) => repairs.map((repair) => (
  `Legend: the ${legendRecordKind(repair)} stored under "${repair.label ?? ''}" (${repair.path}) sets the row "${repair.used ?? ''}",`
  + ` which a later edit also sets; the later edit is used and this one was dropped.`));

/**
 * One line per stored value (the same value at one path, once).
 * @param {readonly LoadRepair[]} repairs
 * @param {(repair: LoadRepair) => string} describe
 */
const valueLines = (repairs, describe) => [...new Map(repairs.map((repair) => [
  JSON.stringify([repair.path, repair.stored, repair.used]), repair
])).values()].map(describe);

/**
 * The lines of one table's records of one code.
 * @param {readonly LoadRepair[]} repairs
 * @param {(table: string, lines: string, count: number) => string} describe
 */
const tableLines = (repairs, describe) => [...new Set(repairs.map(({ path }) => path))].map((table) => {
  const lines = repairs.filter((repair) => repair.path === table).map(({ line }) => line);
  return describe(SESSION_TABLE_LABELS[table] ?? table, cappedList('line', 'lines', lines), lines.length);
});

/**
 * The wording of each code, in the order the dialog lists them, by area:
 * Request settings, Legend, Feature edits, Tables, Files, Results, Other.
 * Every code has an entry (the type requires it).
 * @type {Readonly<Record<LoadRepairCode, (repairs: readonly LoadRepair[]) => string[]>>}
 */
const WORDING = {
  'legend-position-defaulted': (repairs) => valueLines(repairs, ({ path, stored, used }) => (
    `Legend: the position "${stored}" (${path}) is not one of its mode's positions; `
    + `${used === undefined ? 'it is not used' : `${used} is used`}.`)),
  'legend-side-defaulted': (repairs) => valueLines(repairs, ({ path, stored, used }) => (
    `Legend of a saved Result: the side "${stored}" (${path}) is not one of its mode's positions; `
    + `the Result is read with ${used}.`)),
  'legend-color-dropped': (repairs) => {
    const named = [...new Set(repairs.map((repair) => `"${repair.label ?? ''}"${legendListMark(repair)}`))];
    return [`Legend: ${cappedList('entry', 'entries', named)} had a color the app does not accept and ${wasWere(named.length)} dropped.`];
  },
  'legend-entry-duplicate': (repairs) => {
    const named = [...new Set(repairs.map((repair) => `"${repair.label ?? ''}"${legendListMark(repair)}`))];
    return [`Legend: ${cappedList('entry', 'entries', named)} ${named.length === 1 ? 'repeats an earlier entry of its list and was'
      : 'repeat an earlier entry of their list and were'} dropped.`];
  },
  'legend-edit-overwritten': overwrittenLegendEdits,
  'legend-edit-overwritten-session46': overwrittenLegendEdits,
  'legend-edit-rekey-dropped': (repairs) => repairs.map((repair) => {
    const kind = legendRecordKind(repair);
    return `Legend: the ${kind} stored under "${repair.label ?? ''}" (${repair.path}) moves to the row "${repair.used ?? ''}",`
      + ` which has ${kind === 'entry' ? 'an entry' : 'an edit'} of its own; that ${kind === 'entry' ? 'entry' : 'edit'} is used and this one was dropped.`;
  }),
  'feature-edit-unmatched': (repairs) => [
    `${total(repairs)} feature edit(s) from an older Session could not be matched to a feature of its saved diagram and were dropped.`
  ],
  'feature-visibility-narrowed': (repairs) => [
    `${total(repairs)} Feature visibility edit(s) from an older Session hid every feature with the same hash, such as each copy`
    + ' of a duplicated record. Each now applies only to the feature that was edited, so the next Generate draws the others.'
  ],
  'hash-selector-unmapped': (repairs) => [
    `${total(repairs)} hash= value(s) from an older Session could not be renamed to the hash of a feature in its source record:`
    + ' the saved diagram drew a cropped or reverse-complemented record, and the value names features of several'
    + ' source hashes there, or the sources were not read again. hash= now names a feature by its hash in the source'
    + ' record, so each may now match another feature or none.'
  ],
  'table-row-joined': (repairs) => tableLines(repairs, (table, lines) => (
    `${table}: ${lines} had extra cells, joined into the last column with one space.`)),
  'table-row-dropped': (repairs) => tableLines(repairs, (table, lines, count) => (
    `${table}: ${lines} lacked a required column and ${wasWere(count)} dropped.`)),
  'table-row-invalid': (repairs) => tableLines(repairs, (table, lines, count) => (
    `${table}: ${lines} had a color the table does not accept and ${wasWere(count)} dropped.`)),
  'record-display-drafts-dropped': () => [
    'Record display (recordDisplayDrafts): the stored record display settings could not be read and were dropped.'
  ]
};

/**
 * The records the dialog lists.
 * @param {readonly LoadRepair[]} repairs
 */
export const listedLoadRepairs = (repairs) => repairs.filter((repair) => !UNLISTED_CODES.has(repair.code));

/**
 * The dialog for a Load's repairs, or null when it lists none.
 * @param {readonly LoadRepair[]} repairs
 * @returns {LoadRepairDialog | null}
 */
export const describeLoadRepairs = (repairs) => {
  const listed = listedLoadRepairs(repairs);
  if (!listed.length) return null;
  const lines = /** @type {[LoadRepairCode, (repairs: readonly LoadRepair[]) => string[]][]} */ (Object.entries(WORDING))
    .flatMap(([code, describe]) => {
      const ofCode = listed.filter((repair) => repair.code === code);
      return ofCode.length ? describe(ofCode) : [];
    });
  const items = lines.length > LINE_LIMIT
    ? [...lines.slice(0, LINE_LIMIT - 1), `and ${lines.length - (LINE_LIMIT - 1)} more`]
    : lines;
  return {
    title: 'Load this file?',
    description: 'Some parts of this file cannot be fully restored:',
    items,
    note: 'Save Session afterwards writes the current format. The file itself is not changed.'
  };
};
