// @ts-check
import { getSessionResourceSource } from './file-content-cache.js';

/**
 * What owns a resource's payload bytes: the File the user chose, or another
 * object (the descriptor itself when no owner was set).
 * @typedef {WeakKey} ResourcePayloadOwner
 */

/** @type {WeakMap<object, ResourcePayloadOwner>} */
const payloadOwners = new WeakMap();

/** @param {unknown} value */
const objectKey = (value) => (
  value && (typeof value === 'object' || typeof value === 'function') ? value : null
);

/**
 * @template T
 * @param {T} descriptor
 * @param {unknown} owner
 * @returns {T}
 */
export const setResourcePayloadOwner = (descriptor, owner) => {
  const descriptorKey = objectKey(descriptor);
  const ownerKey = objectKey(owner);
  if (descriptorKey && ownerKey) payloadOwners.set(descriptorKey, ownerKey);
  return descriptor;
};

/**
 * @template T
 * @param {T} descriptor
 * @returns {T | ResourcePayloadOwner | null}
 */
export const getResourcePayloadOwner = (descriptor) => {
  const key = objectKey(descriptor);
  return (key && payloadOwners.get(key)) || key;
};

/**
 * The object whose bytes a File, a Session resource view or a descriptor
 * stands for: every view of one Session resource stands for its source
 * descriptor, anything else for its payload owner.
 * @param {object} value
 * @returns {object}
 */
export const getResourcePayloadIdentity = (value) => {
  const owner = /** @type {object} */ (getResourcePayloadOwner(value));
  return /** @type {{ descriptor?: object } | null} */ (getSessionResourceSource(owner))?.descriptor || owner;
};
