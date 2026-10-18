// @ts-check
// The color-rules slice of a drawing: its specific color rules
// (`manualSpecificRules`, Session 46 `modes.<m>.config.rules`) and the Legend
// captions of the rules a rules file brought (`fileLegendCaptions`, never
// saved). The drawing shows both read-only to the type checker; only this
// module writes them, through `drawing.slices.colorRules`: a rule table write
// (`replaceColorRules`: the rule commit and Load), and a History restore and
// Reset (`install`, `reset`). A caller passes rules it normalized.
import { reactive, ref } from '../../utils/reactivity.js';
import { cloneJsonData } from '../json-clone.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';
/** @import { ReadonlyRef, Ref } from '../../utils/reactivity.js' */
/** @import { SliceHandle } from './slice-handle.js' */
/** @import { SpecificColorRule } from '../specific-color-rules.js' */

/** @typedef {{ manualSpecificRules: SpecificColorRule[], fileLegendCaptions: Ref<Set<string>> }} ColorRulesMembers */
/** @typedef {{ fileLegendCaptions: string[], manualSpecificRules: SpecificColorRule[] }} ColorRulesData What `capture` returns. */
/**
 * @typedef {{
 *   readonly manualSpecificRules: readonly Readonly<SpecificColorRule>[],
 *   readonly fileLegendCaptions: ReadonlyRef<Set<string>>
 * }} ColorRulesView
 */
/** @typedef {SliceHandle<'colorRules'>} ColorRulesSlice The handle a writer takes. */

// A color-rules handle exists only through `create`, which gives it these members.
/** @param {ColorRulesSlice} slice */
const membersOf = (slice) => /** @type {ColorRulesMembers} */ (sliceMembers(slice, 'colorRules'));

/**
 * The Legend captions of the rules a rules file brought.
 * @param {readonly Readonly<SpecificColorRule>[]} rules
 * @returns {Set<string>}
 */
export const fileLegendCaptionsOf = (rules) => new Set(rules.filter((rule) => rule.fromFile && rule.cap).map((rule) => rule.cap));

/**
 * A new slice, or one over a drawing's existing members (a test fixture).
 * @param {Partial<ColorRulesMembers>} [members]
 * @returns {{ handle: ColorRulesSlice, view: ColorRulesView }}
 */
const create = (members = {}) => {
  /** @type {ColorRulesMembers} */
  const data = {
    fileLegendCaptions: members.fileLegendCaptions ?? ref(new Set()),
    manualSpecificRules: members.manualSpecificRules ?? reactive(/** @type {SpecificColorRule[]} */ ([]))
  };
  return { handle: createSliceHandle('colorRules', data), view: Object.freeze({ ...data }) };
};

/**
 * @param {ColorRulesSlice} slice
 * @returns {ColorRulesData}
 */
const capture = (slice) => {
  const { manualSpecificRules, fileLegendCaptions } = membersOf(slice);
  return { fileLegendCaptions: [...fileLegendCaptions.value], manualSpecificRules: cloneJsonData([...manualSpecificRules]) };
};

/**
 * The rule table becomes `rules`, and the file captions the ones of its file
 * rules. A rule with a row in `sourceRows` (at its index) is written into that
 * row, so the row keeps its identity (a pattern draft and an input hold it).
 * @param {ColorRulesSlice} slice
 * @param {readonly SpecificColorRule[]} rules
 * @param {readonly (Readonly<SpecificColorRule> | null)[]} [sourceRows] Rows of this table, as the view shows them.
 */
export const replaceColorRules = (slice, rules, sourceRows = []) => {
  const members = membersOf(slice);
  members.manualSpecificRules.splice(0, members.manualSpecificRules.length, ...rules.map((rule, index) => {
    // A source row is a row of this slice's table, which it writes.
    const row = /** @type {SpecificColorRule | null | undefined} */ (sourceRows[index]);
    if (!row) return rule;
    Object.assign(row, rule);
    if (!Object.hasOwn(rule, 'fromFile')) delete row.fromFile;
    return row;
  }));
  members.fileLegendCaptions.value = fileLegendCaptionsOf(rules);
};

// Each member `data` holds becomes a copy of it; the others are kept (a
// History step restores the file captions on their own).
/**
 * @param {ColorRulesSlice} slice
 * @param {{ manualSpecificRules?: readonly Readonly<SpecificColorRule>[], fileLegendCaptions?: readonly string[] }} data
 */
const install = (slice, { manualSpecificRules, fileLegendCaptions }) => {
  const members = membersOf(slice);
  if (manualSpecificRules) {
    members.manualSpecificRules.splice(0, members.manualSpecificRules.length, ...cloneJsonData([...manualSpecificRules]));
  }
  if (fileLegendCaptions) members.fileLegendCaptions.value = new Set(fileLegendCaptions);
};

/** @param {ColorRulesSlice} slice */
const reset = (slice) => install(slice, { manualSpecificRules: [], fileLegendCaptions: [] });

export const COLOR_RULES_SLICE = Object.freeze({
  members: Object.freeze(['fileLegendCaptions', 'manualSpecificRules']), create, capture, install, reset
});
