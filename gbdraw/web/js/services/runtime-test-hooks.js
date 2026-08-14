const testHooks = () => globalThis.__GBDRAW_TEST_HOOKS__;

export const runtimeTestHooksEnabled = () => Boolean(testHooks());

export const recordStructuralMetric = (name, value = 1, detail = {}) => {
  const callback = testHooks()?.onStructuralMetric;
  if (typeof callback !== 'function') return;
  callback({
    name: String(name || ''),
    value: Number(value) || 0,
    ...(detail && typeof detail === 'object' ? detail : {})
  });
};

export const recordSessionLifecycleEvent = (name, detail = {}) => {
  const callback = testHooks()?.onSessionLifecycleEvent;
  if (typeof callback !== 'function') return;
  callback({
    name: String(name || ''),
    timestamp: globalThis.performance?.now?.() ?? Date.now(),
    ...(detail && typeof detail === 'object' ? detail : {})
  });
};
