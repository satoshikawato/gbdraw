// @ts-check
const UNSAFE_OBJECT_KEYS = new Set(['__proto__', 'constructor', 'prototype']);

function* safeObjectKeySteps(value, path) {
  if (!value || typeof value !== 'object') return;
  const pending = [value];
  const seen = new WeakSet();
  while (pending.length) {
    const current = pending.pop();
    if (!current || typeof current !== 'object' || seen.has(current)) continue;
    seen.add(current);
    for (const key of Object.keys(current)) {
      if (UNSAFE_OBJECT_KEYS.has(key)) {
        throw new Error(`${path} contains unsafe key ${key}.`);
      }
      pending.push(current[key]);
      yield;
    }
  }
};


export const assertSafeObjectKeys = (value, path = 'value') => {
  for (const _step of safeObjectKeySteps(value, path)) { /* exhaust validation */ }
};

export const assertSafeObjectKeysForImport = async (value, path = 'value') => {
  let deadline = performance.now() + 16;
  for (const _step of safeObjectKeySteps(value, path)) {
    if (performance.now() >= deadline) {
      await new Promise((resolve) => setTimeout(resolve, 0));
      deadline = performance.now() + 16;
    }
  }
};
