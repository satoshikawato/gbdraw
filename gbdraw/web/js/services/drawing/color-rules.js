// @ts-check
// The color-rules slice of a drawing: its specific color rules
// (`manualSpecificRules`, Session 46 `modes.<m>.config.rules`). The drawing
// shows them read-only to the type checker; only this module writes them,
// through `drawing.slices.colorRules`: a rule table write
// (`replaceColorRules`: the rule commit and Load), and a History restore and
// Reset (`install`, `reset`). A caller passes rules it normalized.
import { reactive } from '../../utils/reactivity.js';
import { cloneJsonData } from '../json-clone.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';
/** @import { SliceHandle } from './slice-handle.js' */
/** @import { SpecificColorRule } from '../specific-color-rules.js' */

/** @typedef {{ manualSpecificRules: SpecificColorRule[] }} ColorRulesData What `capture` returns. */
/** @typedef {{ readonly manualSpecificRules: readonly Readonly<SpecificColorRule>[] }} ColorRulesView */
/** @typedef {SliceHandle<'colorRules'>} ColorRulesSlice The handle a writer takes. */

// A color-rules handle exists only through `create`, which gives it these members.
/** @param {ColorRulesSlice} slice */
const membersOf = (slice) => /** @type {ColorRulesData} */ (sliceMembers(slice, 'colorRules'));

/**
 * A new slice, or one over a drawing's existing members (a test fixture).
 * @param {Partial<ColorRulesData>} [members]
 * @returns {{ handle: ColorRulesSlice, view: ColorRulesView }}
 */
const create = (members = {}) => {
  /** @type {ColorRulesData} */
  const data = { manualSpecificRules: members.manualSpecificRules ?? reactive(/** @type {SpecificColorRule[]} */ ([])) };
  return { handle: createSliceHandle('colorRules', data), view: Object.freeze({ ...data }) };
};

/**
 * @param {ColorRulesSlice} slice
 * @returns {ColorRulesData}
 */
const capture = (slice) => ({ manualSpecificRules: cloneJsonData([...membersOf(slice).manualSpecificRules]) });

/**
 * The rule table becomes `rules`. A rule with a row in `sourceRows` (at its
 * index) is written into that row, so the row keeps its identity (a pattern
 * draft and an input hold it).
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
};

// The rule table becomes a copy of the rules `data` holds; without them it is kept.
/**
 * @param {ColorRulesSlice} slice
 * @param {{ manualSpecificRules?: readonly Readonly<SpecificColorRule>[] }} data
 */
const install = (slice, { manualSpecificRules }) => {
  const rules = membersOf(slice).manualSpecificRules;
  if (manualSpecificRules) rules.splice(0, rules.length, ...cloneJsonData([...manualSpecificRules]));
};

/** @param {ColorRulesSlice} slice */
const reset = (slice) => install(slice, { manualSpecificRules: [] });

export const COLOR_RULES_SLICE = Object.freeze({ members: Object.freeze(['manualSpecificRules']), create, capture, install, reset });
