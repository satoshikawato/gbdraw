// @ts-check
const testHooks = () => globalThis.__GBDRAW_TEST_HOOKS__;

export const runtimeTestHooksEnabled = () => Boolean(testHooks());

/**
 * A test's stand-in for the cap on a Session's expanded JSON, so a browser
 * case needs no Session of 512 MiB.
 * @returns {number | undefined}
 */
export const testExpandedSessionByteCap = () => {
  const cap = testHooks()?.maxExpandedSessionBytes;
  return Number.isSafeInteger(cap) && cap > 0 ? cap : undefined;
};

/**
 * @param {unknown} name
 * @param {unknown} [value]
 * @param {unknown} [detail]
 */
export const recordStructuralMetric = (name, value = 1, detail = {}) => {
  const callback = testHooks()?.onStructuralMetric;
  if (typeof callback !== 'function') return;
  callback({
    name: String(name || ''),
    value: Number(value) || 0,
    ...(detail && typeof detail === 'object' ? detail : {})
  });
};

/**
 * @param {unknown} name
 * @param {unknown} [detail]
 */
export const recordSessionLifecycleEvent = (name, detail = {}) => {
  const callback = testHooks()?.onSessionLifecycleEvent;
  if (typeof callback !== 'function') return;
  callback({
    name: String(name || ''),
    timestamp: globalThis.performance?.now?.() ?? Date.now(),
    ...(detail && typeof detail === 'object' ? detail : {})
  });
};
