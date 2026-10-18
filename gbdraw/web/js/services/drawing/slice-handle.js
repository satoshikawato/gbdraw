// @ts-check
// The handle of a drawing slice (services/drawing/): the value its writers
// take (`drawing.slices.<name>`). To the type checker a handle is opaque
// (`{ readonly __slice: name }`); the slice that made it reads its writable
// members from it. The members are kept under a registered symbol, so a copy
// of a slice module (a unit test's copy of the tree) reads the handles of the
// other copy.
const MEMBERS = Symbol.for('gbdraw.drawing-slice.members');

/**
 * @template {string} N
 * @typedef {{ readonly __slice: N }} SliceHandle
 */

/**
 * @template {string} N
 * @template {object} M
 * @param {N} name
 * @param {M} members The slice's writable members.
 * @returns {SliceHandle<N>}
 */
export const createSliceHandle = (name, members) => Object.freeze(
  Object.defineProperty({ __slice: name }, MEMBERS, { value: members })
);

/**
 * The writable members of a handle, for the slice that made it.
 * @param {SliceHandle<string>} handle
 * @returns {unknown}
 */
export const sliceMembers = (handle) => Reflect.get(handle, MEMBERS);
