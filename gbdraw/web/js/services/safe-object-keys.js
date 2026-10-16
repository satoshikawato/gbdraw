// @ts-check
import { diagnosticError } from '../utils/error-normalization.js';

const UNSAFE_OBJECT_KEYS = new Set(['__proto__', 'constructor', 'prototype']);

/**
 * Walks every own key; an unsafe key is the Session-schema diagnostic that
 * Load and Save give for any other malformed Session structure (OV-952).
 * @param {unknown} value
 */
function* safeObjectKeySteps(value) {
  if (!value || typeof value !== 'object') return;
  /** @type {unknown[]} */
  const pending = [value];
  const seen = new WeakSet();
  while (pending.length) {
    const current = pending.pop();
    if (!current || typeof current !== 'object' || seen.has(current)) continue;
    seen.add(current);
    for (const key of Object.keys(current)) {
      if (UNSAFE_OBJECT_KEYS.has(key)) throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
      // `key` comes from `Object.keys(current)`.
      pending.push(/** @type {Record<string, unknown>} */ (current)[key]);
      yield;
    }
  }
};


/** @param {unknown} value */
export const assertSafeObjectKeys = (value) => {
  for (const _step of safeObjectKeySteps(value)) { /* exhaust validation */ }
};

/** @param {unknown} value */
export const assertSafeObjectKeysForImport = async (value) => {
  let deadline = performance.now() + 16;
  for (const _step of safeObjectKeySteps(value)) {
    if (performance.now() >= deadline) {
      await new Promise((resolve) => setTimeout(resolve, 0));
      deadline = performance.now() + 16;
    }
  }
};
