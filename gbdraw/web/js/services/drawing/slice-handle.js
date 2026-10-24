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
 * The writable members of a handle, for the slice that made it: a handle of
 * another slice has none, so a write through it fails at once instead of
 * writing the other slice's members.
 * @template {string} N
 * @param {SliceHandle<N>} handle
 * @param {N} name The slice's name, as it made the handle.
 * @returns {unknown}
 */
export const sliceMembers = (handle, name) => (handle.__slice === name ? Reflect.get(handle, MEMBERS) : undefined);
