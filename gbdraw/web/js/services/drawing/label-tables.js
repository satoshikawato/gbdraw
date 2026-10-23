// @ts-check
// The label-tables slice of a drawing: its bulk label edits
// (`labelTextBulkOverrides`, {source text: text}) with the saved label table
// they replace (`canonicalLabelOverrideRows`, Session 46
// `modes.<m>.features`), and its label filters: the qualifier priority rules
// (`manualPriorityRules`), the filter mode, the blacklist words and the
// whitelist rules (`config`). The drawing shows them read-only to the type
// checker; only this module writes them, through `drawing.slices.labelTables`:
// the label editor and filter controls' writers below, the auxiliary table
// imports, and Load, a History restore, and Reset (`install`, `reset`).
import { reactive, ref } from '../../utils/reactivity.js';
import { cloneJsonData } from '../json-clone.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';
/** @import { ReadonlyRef, Ref } from '../../utils/reactivity.js' */
/** @import { SliceHandle } from './slice-handle.js' */
/** @import { StoredLabelOverrideRow } from '../label-override-table.js' */
/** @import { QualifierPriorityRule, StoredWhitelistRule } from '../config-field-readers.js' */

/** @typedef {Record<string, string>} BulkLabelEdits A label text by the source text it replaces. */
/** @typedef {string} LabelFilterMode 'None', 'Blacklist' or 'Whitelist' (`readFilterMode` reads it). */
/**
 * @typedef {{
 *   labelTextBulkOverrides: BulkLabelEdits, canonicalLabelOverrideRows: Ref<StoredLabelOverrideRow[]>,
 *   manualPriorityRules: QualifierPriorityRule[], filterMode: Ref<LabelFilterMode>, manualBlacklist: Ref<string>,
 *   manualWhitelist: StoredWhitelistRule[]
 * }} LabelTablesMembers
 */
/**
 * What `capture` returns and `install` writes.
 * @typedef {{
 *   labelTextBulkOverrides: BulkLabelEdits, canonicalLabelOverrideRows: StoredLabelOverrideRow[],
 *   manualPriorityRules: QualifierPriorityRule[], filterMode: LabelFilterMode, manualBlacklist: string,
 *   manualWhitelist: StoredWhitelistRule[]
 * }} LabelTablesData
 */
/**
 * @typedef {{
 *   readonly labelTextBulkOverrides: Readonly<BulkLabelEdits>,
 *   readonly canonicalLabelOverrideRows: ReadonlyRef<readonly Readonly<StoredLabelOverrideRow>[]>,
 *   readonly manualPriorityRules: readonly Readonly<QualifierPriorityRule>[],
 *   readonly filterMode: ReadonlyRef<LabelFilterMode>, readonly manualBlacklist: ReadonlyRef<string>,
 *   readonly manualWhitelist: readonly Readonly<StoredWhitelistRule>[]
 * }} LabelTablesView
 */
/**
 * What `install` copies: `LabelTablesData`, read-only.
 * @typedef {{
 *   labelTextBulkOverrides: Readonly<BulkLabelEdits>, canonicalLabelOverrideRows: readonly Readonly<StoredLabelOverrideRow>[],
 *   manualPriorityRules: readonly Readonly<QualifierPriorityRule>[], filterMode: LabelFilterMode, manualBlacklist: string,
 *   manualWhitelist: readonly Readonly<StoredWhitelistRule>[]
 * }} LabelTablesInput
 */
/** @typedef {SliceHandle<'labelTables'>} LabelTablesSlice The handle a writer takes. */

/** The label filters of a new drawing. */
export const createDefaultLabelFilters = () => ({
  manualPriorityRules: /** @type {QualifierPriorityRule[]} */ ([]),
  filterMode: 'None',
  manualBlacklist: 'hypothetical, uncharacterized, putative, unknown',
  manualWhitelist: /** @type {StoredWhitelistRule[]} */ ([])
});

// A label-tables handle exists only through `create`, which gives it these members.
/** @param {LabelTablesSlice} slice */
const membersOf = (slice) => /** @type {LabelTablesMembers} */ (sliceMembers(slice, 'labelTables'));

/**
 * @template T
 * @param {T[]} list
 * @param {readonly T[]} items
 */
const replaceList = (list, items) => { list.splice(0, list.length, ...items); };

/**
 * @param {Record<string, string>} target
 * @param {Readonly<Record<string, string>>} entries
 */
const replaceEntries = (target, entries) => {
  const copy = { ...entries };
  Object.keys(target).forEach((key) => { delete target[key]; });
  Object.assign(target, copy);
};

/**
 * A new slice, or one over a drawing's existing members (a test fixture).
 * @param {Partial<LabelTablesMembers>} [members]
 * @returns {{ handle: LabelTablesSlice, view: LabelTablesView }}
 */
const create = (members = {}) => {
  const defaults = createDefaultLabelFilters();
  /** @type {LabelTablesMembers} */
  const data = {
    labelTextBulkOverrides: members.labelTextBulkOverrides ?? reactive(/** @type {BulkLabelEdits} */ ({})),
    canonicalLabelOverrideRows: members.canonicalLabelOverrideRows ?? ref(/** @type {StoredLabelOverrideRow[]} */ ([])),
    manualPriorityRules: members.manualPriorityRules ?? reactive(defaults.manualPriorityRules),
    filterMode: members.filterMode ?? ref(defaults.filterMode),
    manualBlacklist: members.manualBlacklist ?? ref(defaults.manualBlacklist),
    manualWhitelist: members.manualWhitelist ?? reactive(defaults.manualWhitelist)
  };
  return { handle: createSliceHandle('labelTables', data), view: Object.freeze({ ...data }) };
};

/**
 * @param {LabelTablesSlice} slice
 * @returns {LabelTablesData}
 */
const capture = (slice) => {
  const members = membersOf(slice);
  return {
    labelTextBulkOverrides: { ...members.labelTextBulkOverrides },
    canonicalLabelOverrideRows: cloneJsonData([...members.canonicalLabelOverrideRows.value]),
    manualPriorityRules: cloneJsonData([...members.manualPriorityRules]),
    filterMode: members.filterMode.value,
    manualBlacklist: members.manualBlacklist.value,
    manualWhitelist: cloneJsonData([...members.manualWhitelist])
  };
};

// Each member `data` holds becomes a copy of it; the others are kept.
/**
 * @param {LabelTablesSlice} slice
 * @param {Partial<LabelTablesInput>} data
 */
const install = (slice, data) => {
  const members = membersOf(slice);
  // Copies of the rows `data` holds (`cloneJsonData` gives writable rows).
  if (data.labelTextBulkOverrides) replaceEntries(members.labelTextBulkOverrides, data.labelTextBulkOverrides);
  if (data.canonicalLabelOverrideRows) {
    members.canonicalLabelOverrideRows.value = /** @type {StoredLabelOverrideRow[]} */ (cloneJsonData([...data.canonicalLabelOverrideRows]));
  }
  if (data.manualPriorityRules) {
    replaceList(members.manualPriorityRules, /** @type {QualifierPriorityRule[]} */ (cloneJsonData([...data.manualPriorityRules])));
  }
  if (data.filterMode !== undefined) members.filterMode.value = data.filterMode;
  if (data.manualBlacklist !== undefined) members.manualBlacklist.value = String(data.manualBlacklist);
  if (data.manualWhitelist) {
    replaceList(members.manualWhitelist, /** @type {StoredWhitelistRule[]} */ (cloneJsonData([...data.manualWhitelist])));
  }
};

/**
 * The saved label table and the bulk label edits as Save and a History step
 * hold them (`features.labelOverrideRows`, `features.labelTextBulkOverrides`).
 * @param {LabelTablesSlice} slice
 */
export const labelTextEditsData = (slice) => {
  const { canonicalLabelOverrideRows, labelTextBulkOverrides } = capture(slice);
  return { labelOverrideRows: canonicalLabelOverrideRows, labelTextBulkOverrides };
};

/** @param {LabelTablesSlice} slice */
const reset = (slice) => install(slice, { labelTextBulkOverrides: {}, canonicalLabelOverrideRows: [], ...createDefaultLabelFilters() });

/**
 * The one way an edit replaces a drawing's bulk label edits. A change also
 * clears the saved label table, which the editor's rows replace, in the
 * action's own History step, so Undo returns both (OV-1102). Session Load and
 * History restore install both members as saved.
 * @param {LabelTablesSlice} slice
 * @param {Readonly<BulkLabelEdits>} edits
 * @returns {boolean} Whether the bulk label edits changed.
 */
export const replaceBulkLabelEdits = (slice, edits) => {
  const members = membersOf(slice);
  const current = members.labelTextBulkOverrides;
  const sourceTexts = Object.keys(edits);
  if (sourceTexts.length === Object.keys(current).length
    && sourceTexts.every((sourceText) => Object.hasOwn(current, sourceText) && current[sourceText] === edits[sourceText])) {
    return false;
  }
  replaceEntries(current, edits);
  members.canonicalLabelOverrideRows.value = [];
  return true;
};

/**
 * Clears the bulk label edits and the saved label table (a new Session).
 * @param {LabelTablesSlice} slice
 */
export const clearLabelTextEdits = (slice) => install(slice, { labelTextBulkOverrides: {}, canonicalLabelOverrideRows: [] });

/**
 * @param {LabelTablesSlice} slice
 * @param {LabelFilterMode} mode
 */
export const setLabelFilterMode = (slice, mode) => { membersOf(slice).filterMode.value = mode; };

/**
 * @param {LabelTablesSlice} slice
 * @param {string} text Comma-separated words.
 */
export const setBlacklistText = (slice, text) => { membersOf(slice).manualBlacklist.value = text; };

/**
 * @param {LabelTablesSlice} slice
 * @param {readonly StoredWhitelistRule[]} rules New rows (copied).
 */
export const addWhitelistRules = (slice, rules) => {
  membersOf(slice).manualWhitelist.push(...rules.map((rule) => ({ ...rule })));
};

/**
 * @param {LabelTablesSlice} slice
 * @param {number} index
 */
export const removeWhitelistRule = (slice, index) => {
  if (index >= 0) membersOf(slice).manualWhitelist.splice(index, 1);
};

/**
 * One field of a whitelist row, which makes the row the user's: it loses its
 * file card (OV-1280). A row of no list of this slice changes nothing.
 * @param {LabelTablesSlice} slice
 * @param {Readonly<StoredWhitelistRule>} rule
 * @param {'feat' | 'qual' | 'key'} field
 * @param {string} value
 */
export const setWhitelistRuleField = (slice, rule, field, value) => {
  const row = membersOf(slice).manualWhitelist.find((candidate) => candidate === rule);
  if (!row) return;
  row[field] = value;
  // The mark of the Whitelist file card (services/file-imports.js).
  delete row.importedFrom;
};

/**
 * Sets the qualifier order of each rule's feature type by hand: a type with
 * a row takes a new row of the user's (it loses its file card, OV-1280),
 * another type is added.
 * @param {LabelTablesSlice} slice
 * @param {readonly Readonly<QualifierPriorityRule>[]} rules
 */
export const setPriorityRules = (slice, rules) => {
  const list = membersOf(slice).manualPriorityRules;
  rules.forEach(({ feat, order }) => {
    const index = list.findIndex((candidate) => candidate.feat === feat);
    if (index >= 0) list.splice(index, 1, { feat, order });
    else list.push({ feat, order });
  });
};

// OV-1280 (D-19): a Whitelist or Qualifier priority file's import and its
// Remove on the card: the rows services/file-imports.js marks and takes back.
/**
 * @param {LabelTablesSlice} slice
 * @param {readonly StoredWhitelistRule[]} rules The whole list (copied).
 */
export const replaceWhitelistRules = (slice, rules) => {
  replaceList(membersOf(slice).manualWhitelist, /** @type {StoredWhitelistRule[]} */ (cloneJsonData([...rules])));
};

/**
 * @param {LabelTablesSlice} slice
 * @param {readonly Readonly<QualifierPriorityRule>[]} rules The whole list (copied).
 */
export const replacePriorityRules = (slice, rules) => {
  replaceList(membersOf(slice).manualPriorityRules, /** @type {QualifierPriorityRule[]} */ (cloneJsonData([...rules])));
};

/**
 * @param {LabelTablesSlice} slice
 * @param {number} index
 */
export const removePriorityRule = (slice, index) => {
  if (index >= 0) membersOf(slice).manualPriorityRules.splice(index, 1);
};

export const LABEL_TABLES_SLICE = Object.freeze({
  members: Object.freeze([
    'labelTextBulkOverrides', 'canonicalLabelOverrideRows', 'manualPriorityRules', 'filterMode', 'manualBlacklist',
    'manualWhitelist'
  ]),
  create, capture, install, reset
});
