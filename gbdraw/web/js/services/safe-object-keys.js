const UNSAFE_OBJECT_KEYS = new Set(['__proto__', 'constructor', 'prototype']);

export const assertSafeObjectKeys = (value, path = 'value') => {
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
    }
  }
};
