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

const QUOTE = 0x22;
const BACKSLASH = 0x5C;
const OPEN_OBJECT = 0x7B;
const CLOSE_OBJECT = 0x7D;
const OPEN_ARRAY = 0x5B;
const CLOSE_ARRAY = 0x5D;
const COMMA = 0x2C;

/**
 * Refuses JSON text whose object holds one key twice, as Python's Session
 * reader refuses it (`json` with a duplicate-key hook): the "not valid JSON"
 * refusal of Load (OV-2694). `text` is JSON that `JSON.parse` read, so one
 * linear pass over it finds each object's keys: a string is a key when it
 * follows `{` or a comma inside an object. Keys compare by their decoded
 * text, so `"a"` and `"a"` are one key, as Python reads them.
 * @param {string} text
 */
export const assertUniqueJsonKeys = (text) => {
  /** @type {(Set<string> | null)[]} The keys of each open object; null for an open array. */
  const open = [];
  let expectsKey = false;
  for (let at = 0; at < text.length; at += 1) {
    const code = text.charCodeAt(at);
    if (code === QUOTE) {
      let end = text.indexOf('"', at + 1);
      for (; ; end = text.indexOf('"', end + 1)) {
        let slashes = 0;
        for (let back = end - 1; text.charCodeAt(back) === BACKSLASH; back -= 1) slashes += 1;
        if (slashes % 2 === 0) break;
      }
      if (expectsKey) {
        const raw = text.slice(at + 1, end);
        // A key with no backslash is its own text.
        const key = raw.includes('\\') ? /** @type {string} */ (JSON.parse(text.slice(at, end + 1))) : raw;
        const keys = /** @type {Set<string>} */ (open[open.length - 1]);
        if (keys.has(key)) {
          throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'JSON_FORMAT' }, { stage: 'parse', operation: 'session-load' });
        }
        keys.add(key);
        expectsKey = false;
      }
      at = end;
    } else if (code === OPEN_OBJECT) {
      open.push(new Set());
      expectsKey = true;
    } else if (code === OPEN_ARRAY) {
      open.push(null);
      expectsKey = false;
    } else if (code === CLOSE_OBJECT || code === CLOSE_ARRAY) {
      open.pop();
      expectsKey = false;
    } else if (code === COMMA) {
      expectsKey = open[open.length - 1] instanceof Set;
    }
  }
};
