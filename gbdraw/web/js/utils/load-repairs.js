// @ts-check

import { SESSION_TABLE_LABELS } from './error-normalization.js';
import { t, tm, tn } from './i18n.js';

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
 * (`readSessionDocument`), and the sentence the dialog says of each.
 * @type {Readonly<Record<string, (path: string) => string>>}
 */
const SESSION_PART_WORDING = Object.freeze({
  results: (path) => t('Session file: the saved Results ({path}) could not be read and are left out.', { path }),
  editorState: (path) => t('Session file: the feature list and the Legend and stroke edits of the saved diagram ({path}) could not be read and are left out.', { path }),
  resources: (path) => t('Session file: the stored input files ({path}) could not be read and are left out.', { path }),
  otherModeResult: (path) => t('Session file: the saved diagram of the other mode ({path}) could not be read and is left out.', { path }),
  modes: (path) => t('Session file: the settings of both drawings ({path}) could not be read and are left out.', { path }),
  'modes.circular': (path) => t('Session file: the settings of the Circular drawing ({path}) could not be read and are left out.', { path }),
  'modes.linear': (path) => t('Session file: the settings of the Linear drawing ({path}) could not be read and are left out.', { path })
});
/** The paths of `SESSION_PART_WORDING`. @type {readonly string[]} */
export const SESSION_PART_PATHS = Object.freeze(Object.keys(SESSION_PART_WORDING));

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
 * The first five items, then how many more.
 * @param {readonly unknown[]} items
 */
const cappedList = (items) => {
  const head = items.slice(0, 5).join(', ');
  return items.length > 5 ? t('{list}, and {count} more', { list: head, count: items.length - 5 }) : head;
};
/** @param {readonly LoadRepair[]} repairs */
const total = (repairs) => repairs.reduce((sum, repair) => sum + (repair.count ?? 1), 0);
/**
 * A Legend list entry as the dialog names it; the list is the last part of the path.
 * @param {LoadRepair} repair
 */
const legendListEntry = (repair) => {
  const label = repair.label ?? '';
  const list = repair.path.split('.').pop();
  if (list === 'dormantEntries') return t('"{label}" (not drawn)', { label });
  if (list === 'deletedEntries') return t('"{label}" (deleted)', { label });
  return `"${label}"`;
};

/** @param {LoadRepair} repair The collection is the last part of the path. */
const legendRecordKind = (repair) => {
  const collection = repair.path.split('.').pop();
  return collection === 'colorOverrides' || collection === 'strokeOverrides' ? collection : 'entry';
};
/** @param {readonly LoadRepair[]} repairs */
const overwrittenLegendEdits = (repairs) => repairs.map((repair) => {
  const params = { label: repair.label ?? '', path: repair.path, used: repair.used ?? '' };
  const kind = legendRecordKind(repair);
  if (kind === 'colorOverrides') return t('Legend: the color edit stored under "{label}" ({path}) sets the row "{used}", which a later edit also sets; the later edit is used and this one was dropped.', params);
  if (kind === 'strokeOverrides') return t('Legend: the stroke edit stored under "{label}" ({path}) sets the row "{used}", which a later edit also sets; the later edit is used and this one was dropped.', params);
  return t('Legend: the entry stored under "{label}" ({path}) sets the row "{used}", which a later edit also sets; the later edit is used and this one was dropped.', params);
});

/**
 * The lines of one table's records of one code.
 * @param {readonly LoadRepair[]} repairs
 * @param {(table: string, lines: string, count: number) => string} describe
 */
const tableLines = (repairs, describe) => [...new Set(repairs.map(({ path }) => path))].map((table) => {
  const lines = repairs.filter((repair) => repair.path === table).map(({ line }) => line);
  // A table's label is a key of the error wording owner; an unknown table is shown by its name.
  return describe(Object.hasOwn(SESSION_TABLE_LABELS, table) ? tm(SESSION_TABLE_LABELS[table]) : table, cappedList(lines), lines.length);
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
    t('Saved diagram: the settings of the saved diagram ({path}) could not be read. The file loads as settings only; the saved Results, the GenBank, GFF and FASTA inputs, the feature list and the caches are left out.', { path }))),
  // The Session document (services/session-authority.js `readSessionDocument`).
  'session-part-dropped': (repairs) => repairs.map(({ path }) => (
    SESSION_PART_WORDING[path]?.(path)
    ?? t('Session file: {name} ({path}) could not be read and is left out.', { name: path, path }))),
  'session-field-dropped': (repairs) => {
    const items = cappedList(repairs.map(({ path, label }) => (label === undefined ? path : `${path} (${label})`)));
    return [tn(repairs.length,
      'Session file: the field {items} could not be read and was left out.',
      'Session file: the fields {items} could not be read and were left out.', { items })];
  },
  'other-result-set-dropped': (repairs) => repairs.map(({ path, label }) => (
    t('Saved {label} diagram: its settings ({path}) could not be read; that diagram and its Results are left out.', { label: label ?? t('other'), path }))),
  'request-part-dropped': (repairs) => repairs.map(({ path, label, effect }) => (
    effect === 'defaulted'
      ? t('Settings: {label} ({path}) could not be read; its default is used.', { label: label ?? t('a part'), path })
      : t('Settings: {label} ({path}) could not be read and was left out, with the settings inside it.', { label: label ?? t('a part'), path }))),
  'request-reference-dropped': (repairs) => repairs.map(({ path, label }) => (
    t('Settings: {path} names a file or record this Session does not hold; {label} is left out.', { path, label: label ?? t('it') }))),
  // Free text is named by its place only (a changed one, such as an ID made
  // unique, names neither value); a value replaced by a list or an object
  // names no value used.
  'setting-replaced': (repairs) => repairs.map(({ path, stored, used, effect }) => {
    const params = { path, stored: stored ?? '', used: used ?? '' };
    if (stored === undefined) {
      if (used !== undefined) return t('Settings: the value of {path} cannot be used; {used} is used.', params);
      if (effect === 'defaulted') return t('Settings: the value of {path} cannot be used; its default is used.', params);
      if (effect === 'changed') return t('Settings: the value of {path} cannot be used; it is changed.', params);
      if (effect === 'dropped') return t('Settings: the value of {path} cannot be used; it is left unset.', params);
      return t('Settings: the value of {path} cannot be used; the drawing keeps its current value.', params);
    }
    if (used !== undefined) return t('Settings: the value "{stored}" ({path}) cannot be used; {used} is used.', params);
    if (effect === 'defaulted') return t('Settings: the value "{stored}" ({path}) cannot be used; its default is used.', params);
    if (effect === 'changed') return t('Settings: the value "{stored}" ({path}) cannot be used; it is changed.', params);
    if (effect === 'dropped') return t('Settings: the value "{stored}" ({path}) cannot be used; it is left unset.', params);
    return t('Settings: the value "{stored}" ({path}) cannot be used; the drawing keeps its current value.', params);
  }),
  'setting-entries-dropped': (repairs) => repairs.map(({ path, count = 1 }) => tn(count,
    'Settings: {count} entry of {path} could not be read and was dropped.',
    'Settings: {count} entries of {path} could not be read and were dropped.', { path })),
  'file-entry-dropped': (repairs) => repairs.map(({ path, count = 1 }) => tn(count,
    'Input files: {count} entry of {path} could not be read and was dropped.',
    'Input files: {count} entries of {path} could not be read and were dropped.', { path })),
  'comparison-series-dropped': (repairs) => repairs.map(({ path, label }) => (label === undefined
    ? t('Comparison series: the series ({path}) has no BLAST table that could be read; it is left out with its sequence files.', { path })
    : t('Comparison series: the series "{label}" ({path}) has no BLAST table that could be read; it is left out with its sequence files.', { label, path }))),
  'pending-palette-dropped': (repairs) => repairs.map(({ path, label }) => (
    t('Palette: the palette "{label}" ({path}) waits to be applied, but Instant Preview is on, which applies a palette at once; it was dropped.', { label: label ?? '', path }))),
  'depth-row-disabled': (repairs) => repairs.map(({ path, label }) => (
    t('Track stack: the Depth row "{label}" ({path}) names a Depth series its drawing does not have; the row is turned off and no longer names that series.', { label: label ?? '', path }))),
  'legend-position-defaulted': (repairs) => repairs.map(({ path, stored, used }) => (used === undefined
    ? t('Legend: the position "{stored}" ({path}) is not one of its mode\'s positions; it is not used.', { stored: stored ?? '', path })
    : t('Legend: the position "{stored}" ({path}) is not one of its mode\'s positions; {used} is used.', { stored: stored ?? '', path, used }))),
  'legend-side-defaulted': (repairs) => repairs.map(({ path, stored, used }) => (
    t('Legend of a saved Result: the side "{stored}" ({path}) is not one of its mode\'s positions; the Result is read with {used}.', { stored: stored ?? '', path, used: used ?? '' }))),
  'legend-color-dropped': (repairs) => {
    const named = [...new Set(repairs.map(legendListEntry))];
    return [tn(named.length,
      'Legend: entry {items} had a color the app does not accept and was dropped.',
      'Legend: entries {items} had a color the app does not accept and were dropped.', { items: cappedList(named) })];
  },
  'legend-entry-duplicate': (repairs) => {
    const named = [...new Set(repairs.map(legendListEntry))];
    return [tn(named.length,
      'Legend: entry {items} repeats an earlier entry of its list and was dropped.',
      'Legend: entries {items} repeat an earlier entry of their list and were dropped.', { items: cappedList(named) })];
  },
  'legend-edit-overwritten': overwrittenLegendEdits,
  'legend-edit-overwritten-session46': overwrittenLegendEdits,
  'legend-edit-rekey-dropped': (repairs) => repairs.map((repair) => {
    const params = { label: repair.label ?? '', path: repair.path, used: repair.used ?? '' };
    const kind = legendRecordKind(repair);
    if (kind === 'colorOverrides') return t('Legend: the color edit stored under "{label}" ({path}) moves to the row "{used}", which has an edit of its own; that edit is used and this one was dropped.', params);
    if (kind === 'strokeOverrides') return t('Legend: the stroke edit stored under "{label}" ({path}) moves to the row "{used}", which has an edit of its own; that edit is used and this one was dropped.', params);
    return t('Legend: the entry stored under "{label}" ({path}) moves to the row "{used}", which has an entry of its own; that entry is used and this one was dropped.', params);
  }),
  'feature-edit-unmatched': (repairs) => [
    t('{count} feature edit(s) from an older Session could not be matched to a feature of its saved diagram and were dropped.', { count: total(repairs) })
  ],
  'feature-visibility-narrowed': (repairs) => [
    t('{count} Feature visibility edit(s) from an older Session hid every feature with the same hash, such as each copy of a duplicated record. Each now applies only to the feature that was edited, so the next Generate draws the others.', { count: total(repairs) })
  ],
  'hash-selector-unmapped': (repairs) => [
    t('{count} hash= value(s) from an older Session could not be renamed to the hash of a feature in its source record: the saved diagram drew a cropped or reverse-complemented record, and the value names features of several source hashes there, or the sources were not read again. hash= now names a feature by its hash in the source record, so each may now match another feature or none.', { count: total(repairs) })
  ],
  'table-row-joined': (repairs) => tableLines(repairs, (table, lines, count) => tn(count,
    '{table}: line {lines} had extra cells, joined into the last column with one space.',
    '{table}: lines {lines} had extra cells, joined into the last column with one space.', { table, lines })),
  'table-row-dropped': (repairs) => tableLines(repairs, (table, lines, count) => tn(count,
    '{table}: line {lines} lacked a required column and was dropped.',
    '{table}: lines {lines} lacked a required column and were dropped.', { table, lines })),
  'table-row-invalid': (repairs) => tableLines(repairs, (table, lines, count) => tn(count,
    '{table}: line {lines} had a color the table does not accept and was dropped.',
    '{table}: lines {lines} had a color the table does not accept and were dropped.', { table, lines })),
  'orthogroup-edit-dropped': (repairs) => repairs.map(({ path, count = 1 }) => (path.endsWith('DescriptionOverrides')
    ? t('Similarity groups: {count} description edit(s) ({path}) name a group the saved Result does not have and were dropped.', { count, path })
    : t('Similarity groups: {count} name edit(s) ({path}) name a group the saved Result does not have and were dropped.', { count, path }))),
  'cache-entry-dropped': (repairs) => repairs.map(({ path, count = 1 }) => (
    t('LOSAT cache: {count} saved value(s) ({path}) could not be read and were dropped.', { count, path }))),
  'feature-recovery-failed': () => [
    t('Features: the feature information of the saved diagram could not be restored from this file; the SVG preview stays, and the next Generate restores it.')
  ],
  'match-sequences-unrestored': () => [
    t('Comparisons: the sequences of the saved comparison matches could not be read again from the input files; the SVG preview stays, and the next Generate reads them.')
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
    ? [...lines.slice(0, LINE_LIMIT - 1), t('and {count} more', { count: lines.length - (LINE_LIMIT - 1) })]
    : lines;
  return {
    title: t('Load this file?'),
    description: t('Some parts of this file cannot be fully restored:'),
    items,
    note: t('Save Session afterwards writes the current format. The file itself is not changed.')
  };
};
