const FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id';
const FEATURE_PART_SUFFIX_RE = /__part\d+$/;
const FEATURE_SELECTOR = [
  `path[${FEATURE_ID_ATTRIBUTE}]`,
  `polygon[${FEATURE_ID_ATTRIBUTE}]`,
  `rect[${FEATURE_ID_ATTRIBUTE}]`,
  'path[id^="f"]',
  'polygon[id^="f"]',
  'rect[id^="f"]'
].join(', ');

const normalizeFeatureIdentity = (value) =>
  String(value || '').trim().replace(FEATURE_PART_SUFFIX_RE, '');

const getFeatureIdentity = (element) =>
  normalizeFeatureIdentity(
    element?.getAttribute?.(FEATURE_ID_ATTRIBUTE) ||
    element?.getAttribute?.('id') ||
    element?.id ||
    ''
  );

const normalizeVisibilityMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  if (normalized === 'suppress') return 'exclude_matching';
  return ['on', 'off', 'exclude_matching'].includes(normalized) ? normalized : 'default';
};

const normalizeChanges = (changes) => {
  if (!Array.isArray(changes)) return [];
  const byFeatureId = new Map();
  changes.forEach((change) => {
    const featureId = normalizeFeatureIdentity(change?.featureId ?? change?.svgId ?? change?.id);
    if (!featureId) return;
    byFeatureId.set(featureId, { ...change, featureId });
  });
  return Array.from(byFeatureId.values());
};

const makeRuntime = (resultIndex, svg) => ({
  resultIndex,
  svg,
  dirty: false,
  dirtyReasons: new Set(),
  indexes: {
    features: null,
    legend: null,
    pairwiseMatches: null,
    orthogroupComparisons: null
  }
});

export const createPreviewRuntime = ({ state, serializeSvg }) => {
  if (!state) throw new Error('createPreviewRuntime requires state.');
  if (typeof serializeSvg !== 'function') throw new Error('createPreviewRuntime requires serializeSvg.');

  let activeRuntime = null;

  const getMountedSvg = () => state.svgContainer?.value?.querySelector?.('svg') || null;

  const mountResultSvg = (resultIndex = state.selectedResultIndex?.value || 0, svg = getMountedSvg()) => {
    if (!svg) {
      activeRuntime = null;
      return null;
    }
    activeRuntime = makeRuntime(Number(resultIndex) || 0, svg);
    return activeRuntime;
  };

  const clearActiveRuntime = () => {
    activeRuntime = null;
  };

  const getActiveRuntime = () => activeRuntime;

  const ensureRuntimeForCurrentSvg = () => {
    const svg = getMountedSvg();
    if (!svg) return null;
    const resultIndex = Number(state.selectedResultIndex?.value || 0);
    if (!activeRuntime || activeRuntime.svg !== svg || activeRuntime.resultIndex !== resultIndex) {
      activeRuntime = makeRuntime(resultIndex, svg);
    }
    return activeRuntime;
  };

  const invalidatePreviewIndexes = (reason = 'unknown', keys = null) => {
    const runtime = activeRuntime;
    if (!runtime) return;
    const targetKeys = Array.isArray(keys) && keys.length > 0
      ? keys
      : Object.keys(runtime.indexes);
    targetKeys.forEach((key) => {
      if (Object.prototype.hasOwnProperty.call(runtime.indexes, key)) {
        runtime.indexes[key] = null;
      }
    });
    runtime.lastInvalidationReason = String(reason || 'unknown');
  };

  const markActiveResultDirty = (reason = 'preview-edit') => {
    const runtime = activeRuntime || ensureRuntimeForCurrentSvg();
    if (!runtime?.svg) return false;
    runtime.dirty = true;
    runtime.dirtyReasons.add(String(reason || 'preview-edit'));
    return true;
  };

  const flushActiveResult = ({ force = false, markIncremental = true } = {}) => {
    const runtime = activeRuntime || (force ? ensureRuntimeForCurrentSvg() : null);
    if (!runtime?.svg) return false;
    if (!force && !runtime.dirty) return false;

    const resultIndex = Number(runtime.resultIndex);
    if (!Number.isInteger(resultIndex) || resultIndex < 0 || resultIndex >= state.results.value.length) {
      runtime.dirty = false;
      runtime.dirtyReasons.clear();
      return false;
    }

    const content = serializeSvg(runtime.svg);
    if (state.results.value[resultIndex]?.content === content) {
      runtime.dirty = false;
      runtime.dirtyReasons.clear();
      return false;
    }
    if (markIncremental && state.skipCaptureBaseConfig) state.skipCaptureBaseConfig.value = true;
    state.results.value[resultIndex] = {
      ...state.results.value[resultIndex],
      content
    };
    runtime.dirty = false;
    runtime.dirtyReasons.clear();
    return true;
  };

  const selectResult = (index) => {
    const count = Array.isArray(state.results.value) ? state.results.value.length : 0;
    const numeric = Number(index);
    const nextIndex = Number.isInteger(numeric) ? Math.max(0, Math.min(numeric, Math.max(0, count - 1))) : 0;
    if (state.selectedResultIndex.value === nextIndex) {
      flushActiveResult();
      return false;
    }
    flushActiveResult({ markIncremental: false });
    state.selectedResultIndex.value = nextIndex;
    clearActiveRuntime();
    return true;
  };

  const buildFeatureIndex = (runtime) => {
    const indexed = new Map();
    if (!runtime?.svg) return indexed;
    Array.from(runtime.svg.querySelectorAll?.(FEATURE_SELECTOR) || []).forEach((element) => {
      const featureId = getFeatureIdentity(element);
      if (!featureId) return;
      if (!indexed.has(featureId)) indexed.set(featureId, []);
      indexed.get(featureId).push(element);
    });
    runtime.indexes.features = indexed;
    return indexed;
  };

  const getFeatureElements = (featureId) => {
    const normalizedId = normalizeFeatureIdentity(featureId);
    const runtime = activeRuntime || ensureRuntimeForCurrentSvg();
    if (!runtime?.svg || !normalizedId) return [];

    const featureIndex = runtime.indexes.features || buildFeatureIndex(runtime);
    const indexed = featureIndex.get(normalizedId);
    if (indexed?.length) return indexed;

    const byId = runtime.svg.getElementById?.(normalizedId);
    return byId ? [byId] : [];
  };

  const applyFeatureFillChanges = (changes, { reason = 'feature-fill' } = {}) => {
    const normalized = normalizeChanges(changes);
    if (normalized.length === 0) return false;

    let updated = 0;
    normalized.forEach((change) => {
      const color = String(change?.color || '').trim();
      if (!color) return;
      getFeatureElements(change.featureId).forEach((element) => {
        element.setAttribute('fill', color);
        updated += 1;
      });
    });

    if (updated > 0) markActiveResultDirty(reason);
    return updated > 0;
  };

  const applyFeatureVisibilityChanges = (changes, { reason = 'feature-visibility' } = {}) => {
    const normalized = normalizeChanges(changes);
    if (normalized.length === 0) return false;

    let updated = 0;
    normalized.forEach((change) => {
      const mode = normalizeVisibilityMode(change?.mode);
      getFeatureElements(change.featureId).forEach((element) => {
        if (mode === 'off') {
          element.setAttribute('display', 'none');
        } else {
          element.removeAttribute('display');
        }
        updated += 1;
      });
    });

    if (updated > 0) markActiveResultDirty(reason);
    return updated > 0;
  };

  const applyFeatureStrokeChanges = (changes, { reason = 'feature-stroke' } = {}) => {
    const normalized = normalizeChanges(changes);
    if (normalized.length === 0) return false;

    let updated = 0;
    normalized.forEach((change) => {
      const strokeColor = String(change?.strokeColor || '').trim();
      const strokeWidth = change?.strokeWidth;
      const hasStrokeWidth = strokeWidth !== null && strokeWidth !== undefined && strokeWidth !== '';
      getFeatureElements(change.featureId).forEach((element) => {
        if (strokeColor) element.setAttribute('stroke', strokeColor);
        if (hasStrokeWidth) element.setAttribute('stroke-width', Number(strokeWidth));
        updated += 1;
      });
    });

    if (updated > 0) markActiveResultDirty(reason);
    return updated > 0;
  };

  const applyLegendChanges = (_changes, { reason = 'legend' } = {}) => markActiveResultDirty(reason);

  return {
    applyFeatureFillChanges,
    applyFeatureStrokeChanges,
    applyFeatureVisibilityChanges,
    applyLegendChanges,
    clearActiveRuntime,
    flushActiveResult,
    getActiveRuntime,
    getFeatureElements,
    invalidatePreviewIndexes,
    markActiveResultDirty,
    mountResultSvg,
    selectResult
  };
};
