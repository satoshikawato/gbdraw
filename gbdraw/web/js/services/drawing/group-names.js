// @ts-check
// The group-names slice of a drawing: the user's names and descriptions of
// Similarity groups, keyed by group id, and the dormant names of groups whose
// members no longer form one group, keyed by member set (D-21; Session 46
// `modes.<m>.config.webEdits`). A drawing shows them read-only to the type
// checker; only this module writes them, through `drawing.slices.groupNames`:
// the rename writers, and `install`, `reset`, and `capture` for Load, a
// History restore, a Generate's regrouping, and Reset.
import { reactive } from '../../utils/reactivity.js';
import { cloneJsonData } from '../json-clone.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';
/** @import { SliceHandle } from './slice-handle.js' */
/** @import { OrthogroupOverrideText } from '../orthogroup-feature-metadata.js' */

/**
 * @typedef {object} GroupNamesData What `capture` returns and `install` writes.
 * @property {Record<string, string>} orthogroupNameOverrides User names by group id.
 * @property {Record<string, string>} orthogroupDescriptionOverrides User descriptions by group id.
 * @property {Record<string, OrthogroupOverrideText>} orthogroupDormantOverrides Dormant names by member-set key.
 */
/**
 * @typedef {{
 *   readonly orthogroupNameOverrides: Readonly<Record<string, string>>,
 *   readonly orthogroupDescriptionOverrides: Readonly<Record<string, string>>,
 *   readonly orthogroupDormantOverrides: Readonly<Record<string, Readonly<OrthogroupOverrideText>>>
 * }} GroupNamesView
 */
/** @typedef {SliceHandle<'groupNames'>} GroupNamesSlice The handle a writer takes. */

const MEMBERS = /** @type {const} */ (['orthogroupNameOverrides', 'orthogroupDescriptionOverrides', 'orthogroupDormantOverrides']);

// A group-names handle exists only through `create`, which gives it these members.
/** @param {GroupNamesSlice} slice */
const membersOf = (slice) => /** @type {GroupNamesData} */ (sliceMembers(slice, 'groupNames'));

// The copy is made first, so a member installed into itself keeps its entries.
/**
 * @template T
 * @param {Record<string, T>} target
 * @param {Readonly<Record<string, T>>} source
 */
const replaceEntries = (target, source) => {
  const entries = cloneJsonData({ ...source });
  Object.keys(target).forEach((key) => delete target[key]);
  Object.assign(target, entries);
};

/**
 * A new slice, or one over a drawing's existing members (a test fixture).
 * @param {Partial<GroupNamesData>} [members]
 * @returns {{ handle: GroupNamesSlice, view: GroupNamesView }}
 */
const create = (members = {}) => {
  /** @type {GroupNamesData} */
  const data = {
    orthogroupNameOverrides: members.orthogroupNameOverrides ?? reactive({}),
    orthogroupDescriptionOverrides: members.orthogroupDescriptionOverrides ?? reactive({}),
    orthogroupDormantOverrides: members.orthogroupDormantOverrides ?? reactive({})
  };
  return { handle: createSliceHandle('groupNames', data), view: Object.freeze({ ...data }) };
};

/**
 * @param {GroupNamesSlice} slice
 * @returns {GroupNamesData}
 */
const capture = (slice) => cloneJsonData({ ...membersOf(slice) });

// Replaces each member `data` holds, with a copy; a member it does not hold
// keeps its entries (an older Session holds no dormant names).
/**
 * @param {GroupNamesSlice} slice
 * @param {Partial<GroupNamesData>} data
 */
const install = (slice, data) => {
  const members = membersOf(slice);
  MEMBERS.forEach((member) => {
    const source = data[member];
    if (source) replaceEntries(/** @type {Record<string, unknown>} */ (members[member]), source);
  });
};

/** @param {GroupNamesSlice} slice */
const reset = (slice) => install(slice, { orthogroupNameOverrides: {}, orthogroupDescriptionOverrides: {}, orthogroupDormantOverrides: {} });

/**
 * A group's user name (an empty name removes it).
 * @param {GroupNamesSlice} slice
 * @param {string} groupId
 * @param {string} name
 */
export const setGroupName = (slice, groupId, name) => {
  const names = membersOf(slice).orthogroupNameOverrides;
  if (name) names[groupId] = name;
  else delete names[groupId];
};

/**
 * A group's user description (an empty description removes it).
 * @param {GroupNamesSlice} slice
 * @param {string} groupId
 * @param {string} description
 */
export const setGroupDescription = (slice, groupId, description) => {
  const descriptions = membersOf(slice).orthogroupDescriptionOverrides;
  if (description) descriptions[groupId] = description;
  else delete descriptions[groupId];
};

/** Removes the dormant names (Clear dormant names). @param {GroupNamesSlice} slice */
export const clearDormantGroupNames = (slice) => install(slice, { orthogroupDormantOverrides: {} });

/**
 * Removes a group's user name and description.
 * @param {GroupNamesSlice} slice
 * @param {string} groupId
 */
export const clearGroupRename = (slice, groupId) => {
  setGroupName(slice, groupId, '');
  setGroupDescription(slice, groupId, '');
};

export const GROUP_NAMES_SLICE = Object.freeze({ members: Object.freeze([...MEMBERS]), create, capture, install, reset });
