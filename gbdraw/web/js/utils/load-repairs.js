// @ts-check

import { SESSION_TABLE_LABELS } from './error-normalization.js';

/**
 * What Load did with a value it could not restore as stored. `kept`: the
 * value is not used and the drawing keeps its current value (a configuration
 * file is installed over the drawing).
 * @typedef {'dropped' | 'defaulted' | 'changed' | 'kept'} LoadRepairEffect
 */
/**
 * The rule of the producer that could not restore a value: one member per rule.
 * @typedef {'request-dropped'
 *   | 'session-part-dropped'
 *   | 'session-field-dropped'
 *   | 'other-result-set-dropped'
 *   | 'request-part-dropped'
 *   | 'request-reference-dropped'
 *   | 'setting-replaced'
 *   | 'setting-entries-dropped'
 *   | 'file-entry-dropped'
 *   | 'comparison-series-dropped'
 *   | 'pending-palette-dropped'
 *   | 'depth-row-disabled'
 *   | 'orthogroup-edit-dropped'
 *   | 'cache-entry-dropped'
 *   | 'feature-edit-unmatched'
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
 *   | 'feature-recovery-failed'
 *   | 'match-sequences-unrestored'
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
 * dialog row (Owner answer Q1 of the Load repair receipt, BA-OD88). A Session
 * 46 drew the later of two Legend edits of one row; the next comparison
 * computes a LOSAT cache entry again; a pending palette waits while Instant
 * Preview, which the file turns on, applies a palette at once; a similarity
 * group name names no group the Result has. (Load keeps a Depth series
 * without a Depth file, as the app does, BA-OD94.)
 * @type {ReadonlySet<LoadRepairCode>}
 */
const UNLISTED_CODES = new Set(/** @type {LoadRepairCode[]} */ ([
  'legend-edit-overwritten-session46', 'cache-entry-dropped', 'pending-palette-dropped', 'orthogroup-edit-dropped'
]));

/**
 * The content parts of a Session the reader may drop or replace whole
 * (`readSessionDocument`), how the dialog names each, and the verb its name takes.
 * @type {Readonly<Record<string, readonly [name: string, verb: 'is' | 'are']>>}
 */
const SESSION_PART_NAMES = Object.freeze({
  results: ['the saved Results', 'are'],
  editorState: ['the feature list and the Legend and stroke edits of the saved diagram', 'are'],
  resources: ['the stored input files', 'are'],
  otherModeResult: ['the saved diagram of the other mode', 'is'],
  modes: ['the settings of both drawings', 'are'],
  'modes.circular': ['the settings of the Circular drawing', 'are'],
  'modes.linear': ['the settings of the Linear drawing', 'are']
});
/** The paths of `SESSION_PART_NAMES`. @type {readonly string[]} */
export const SESSION_PART_PATHS = Object.freeze(Object.keys(SESSION_PART_NAMES));

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
  // A stored request (services/session-request.js `readStoredRenderRequest`).
  'request-dropped': (repairs) => repairs.map(({ path }) => (
    `Saved diagram: the settings of the saved diagram (${path}) could not be read. The file loads as settings only;`
    + ' the saved Results, the GenBank, GFF and FASTA inputs, the feature list and the caches are left out.')),
  // The Session document (services/session-authority.js `readSessionDocument`).
  'session-part-dropped': (repairs) => repairs.map(({ path }) => {
    const [name, verb] = SESSION_PART_NAMES[path] ?? [path, 'is'];
    return `Session file: ${name} (${path}) could not be read and ${verb} left out.`;
  }),
  'session-field-dropped': (repairs) => [
    `Session file: ${cappedList('the field', 'the fields', repairs.map(({ path, label }) => (label === undefined ? path : `${path} (${label})`)))}`
    + ` could not be read and ${wasWere(repairs.length)} left out.`
  ],
  'other-result-set-dropped': (repairs) => repairs.map(({ path, label }) => (
    `Saved ${label ?? 'other'} diagram: its settings (${path}) could not be read; that diagram and its Results are left out.`)),
  'request-part-dropped': (repairs) => repairs.map(({ path, label, effect }) => (
    `Settings: ${label ?? 'a part'} (${path}) could not be read${effect === 'defaulted'
      ? '; its default is used.' : ' and was left out, with the settings inside it.'}`)),
  'request-reference-dropped': (repairs) => repairs.map(({ path, label }) => (
    `Settings: ${path} names a file or record this Session does not hold; ${label ?? 'it'} is left out.`)),
  // Free text is named by its place only (a changed one, such as an ID made
  // unique, names neither value); a value replaced by a list or an object
  // names no value used.
  'setting-replaced': (repairs) => repairs.map(({ path, stored, used, effect }) => (
    `Settings: ${stored === undefined ? `the value of ${path}` : `the value "${stored}" (${path})`} cannot be used; `
    + `${used !== undefined ? `${used} is used` : {
      defaulted: 'its default is used', changed: 'it is changed', dropped: 'it is left unset', kept: 'the drawing keeps its current value'
    }[effect]}.`)),
  'setting-entries-dropped': (repairs) => repairs.map(({ path, count = 1 }) => (
    `Settings: ${count} ${count === 1 ? 'entry' : 'entries'} of ${path} could not be read and ${wasWere(count)} dropped.`)),
  'file-entry-dropped': (repairs) => repairs.map(({ path, count = 1 }) => (
    `Input files: ${count} ${count === 1 ? 'entry' : 'entries'} of ${path} could not be read and ${wasWere(count)} dropped.`)),
  'comparison-series-dropped': (repairs) => repairs.map(({ path, label }) => (
    `Comparison series: the series ${label === undefined ? '' : `"${label}" `}(${path}) has no BLAST table that could be read;`
    + ' it is left out with its sequence files.')),
  'pending-palette-dropped': (repairs) => repairs.map(({ path, label }) => (
    `Palette: the palette "${label ?? ''}" (${path}) waits to be applied, but Instant Preview is on, which applies a palette at once;`
    + ' it was dropped.')),
  'depth-row-disabled': (repairs) => repairs.map(({ path, label }) => (
    `Track stack: the Depth row "${label ?? ''}" (${path}) names a Depth series its drawing does not have;`
    + ' the row is turned off and no longer names that series.')),
  'legend-position-defaulted': (repairs) => repairs.map(({ path, stored, used }) => (
    `Legend: the position "${stored}" (${path}) is not one of its mode's positions; `
    + `${used === undefined ? 'it is not used' : `${used} is used`}.`)),
  'legend-side-defaulted': (repairs) => repairs.map(({ path, stored, used }) => (
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
  'orthogroup-edit-dropped': (repairs) => repairs.map(({ path, count = 1 }) => (
    `Similarity groups: ${count} ${path.endsWith('DescriptionOverrides') ? 'description' : 'name'} edit(s) (${path})`
    + ' name a group the saved Result does not have and were dropped.')),
  'cache-entry-dropped': (repairs) => repairs.map(({ path, count = 1 }) => (
    `LOSAT cache: ${count} saved value(s) (${path}) could not be read and were dropped.`)),
  'feature-recovery-failed': () => [
    'Features: the feature information of the saved diagram could not be restored from this file;'
    + ' the SVG preview stays, and the next Generate restores it.'
  ],
  'match-sequences-unrestored': () => [
    'Comparisons: the sequences of the saved comparison matches could not be read again from the input files;'
    + ' the SVG preview stays, and the next Generate reads them.'
  ]
};

/**
 * A Load's records, one per place in the file and reading: a Session 27-44
 * keeps one draft that both drawings read, and each drawing's reader records
 * what it cannot restore at that one place (OV-765). The first of equal
 * records is kept, in order.
 * @param {readonly LoadRepair[]} repairs
 * @returns {LoadRepair[]}
 */
export const distinctLoadRepairs = (repairs) => [...new Map(repairs.map((repair) => [
  JSON.stringify([repair.code, repair.path, repair.effect, repair.stored, repair.used, repair.count, repair.label, repair.line]), repair
])).values()];

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
