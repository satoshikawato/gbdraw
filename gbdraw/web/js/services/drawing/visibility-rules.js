// @ts-check
// The visibility-rules slice of a drawing: its Feature Visibility rules
// (`featureVisibilityManualRules`, Session 46 `modes.<m>.features`). The
// drawing shows the rule list read-only to the type checker; only this module
// writes it, through `drawing.slices.visibilityRules`: the rule edit
// transition (`editFeatureVisibilityRules`, R10) and Load, a History restore,
// and Reset (`install`, `reset`). A caller passes rules it normalized.
import { reactive } from '../../utils/reactivity.js';
import { cloneJsonData } from '../json-clone.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';
/** @import { SliceHandle } from './slice-handle.js' */
/** @import { FeatureVisibilityRuleDraft } from '../feature-visibility.js' */

/** @typedef {{ featureVisibilityManualRules: FeatureVisibilityRuleDraft[] }} VisibilityRulesData What `capture` returns. */
/** @typedef {{ readonly featureVisibilityManualRules: readonly Readonly<FeatureVisibilityRuleDraft>[] }} VisibilityRulesView */
/** @typedef {SliceHandle<'visibilityRules'>} VisibilityRulesSlice The handle a writer takes. */

// A visibility-rules handle exists only through `create`, which gives it these members.
/** @param {VisibilityRulesSlice} slice */
const membersOf = (slice) => /** @type {VisibilityRulesData} */ (sliceMembers(slice));

/**
 * A new slice, or one over a drawing's existing members (a test fixture).
 * @param {Partial<VisibilityRulesData>} [members]
 * @returns {{ handle: VisibilityRulesSlice, view: VisibilityRulesView }}
 */
const create = (members = {}) => {
  const data = {
    featureVisibilityManualRules: members.featureVisibilityManualRules
      ?? reactive(/** @type {FeatureVisibilityRuleDraft[]} */ ([]))
  };
  return {
    handle: createSliceHandle('visibilityRules', data),
    view: Object.freeze({ featureVisibilityManualRules: data.featureVisibilityManualRules })
  };
};

/**
 * @param {VisibilityRulesSlice} slice
 * @returns {VisibilityRulesData}
 */
const capture = (slice) => ({ featureVisibilityManualRules: cloneJsonData([...membersOf(slice).featureVisibilityManualRules]) });

/**
 * The rule list becomes `rules` (their objects).
 * @param {VisibilityRulesSlice} slice
 * @param {readonly FeatureVisibilityRuleDraft[]} rules
 */
export const replaceVisibilityRules = (slice, rules) => {
  const list = membersOf(slice).featureVisibilityManualRules;
  list.splice(0, list.length, ...rules);
};

/**
 * @param {VisibilityRulesSlice} slice
 * @param {{ featureVisibilityManualRules: readonly FeatureVisibilityRuleDraft[] }} data
 */
const install = (slice, { featureVisibilityManualRules }) => replaceVisibilityRules(slice, featureVisibilityManualRules);

/** @param {VisibilityRulesSlice} slice */
const reset = (slice) => replaceVisibilityRules(slice, []);

export const VISIBILITY_RULES_SLICE = Object.freeze({
  members: Object.freeze(['featureVisibilityManualRules']), create, capture, install, reset
});
