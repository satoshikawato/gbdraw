// @ts-check
import { diagnosticError } from '../utils/error-normalization.js';

const UNSAFE_OBJECT_KEYS = Object.freeze(['__proto__', 'constructor', 'prototype']);

/**
 * Whether a key of stored data names an object's prototype or constructor
 * rather than a field: a merge skips it, and Load refuses a Session that holds it.
 * @param {string} key
 */
export const isUnsafeObjectKey = (key) => UNSAFE_OBJECT_KEYS.includes(key);

/**
 * Walks every own key; an unsafe key is the Session-schema diagnostic that
 * Load and Save give for any other malformed Session structure (OV-952).
 * Session Load runs this yielding walk in the import Worker (`workers/session-import-worker.js`); Load's authority
 * inventory (`validateSessionAuthorityInventory(..., { load: true })`) still runs one synchronous walk on the main thread.
 * @param {unknown} value
 */
export const assertSafeObjectKeys = (value) => {
  if (!value || typeof value !== 'object') return;
  /** @type {unknown[]} */
  const pending = [value];
  const seen = new WeakSet();
  while (pending.length) {
    const current = pending.pop();
    if (!current || typeof current !== 'object' || seen.has(current)) continue;
    seen.add(current);
    for (const key of Object.keys(current)) {
      if (isUnsafeObjectKey(key)) throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
      // `key` comes from `Object.keys(current)`.
      pending.push(/** @type {Record<string, unknown>} */ (current)[key]);
    }
  }
};
