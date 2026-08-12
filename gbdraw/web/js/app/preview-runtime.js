import {
  filterFeatureFillTargets,
  getFeatureElementIndex,
  getFeatureIdentity,
  normalizeFeatureIdentity
} from './feature-dom.js';
import {
  markCommittedSvgResultMounted,
  markCommittedSvgResultUnmounted
} from '../services/svg-result-ingestion.js';

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

  const releaseActiveResult = () => {
    if (!activeRuntime) return;
    markCommittedSvgResultUnmounted(state.results.value[activeRuntime.resultIndex]);
  };

  const mountResultSvg = (resultIndex = state.selectedResultIndex?.value || 0, svg = getMountedSvg()) => {
    if (!svg) {
      releaseActiveResult();
      activeRuntime = null;
      return null;
    }
    const normalizedIndex = Number(resultIndex) || 0;
    if (
      activeRuntime?.svg === svg
      && activeRuntime.resultIndex === normalizedIndex
    ) {
      return activeRuntime;
    }
    releaseActiveResult();
    activeRuntime = makeRuntime(normalizedIndex, svg);
    markCommittedSvgResultMounted(state.results.value[normalizedIndex]);
    return activeRuntime;
  };

  const clearActiveRuntime = () => {
    releaseActiveResult();
    activeRuntime = null;
  };

  const getActiveRuntime = () => activeRuntime;

  const ensureRuntimeForCurrentSvg = () => {
    const svg = getMountedSvg();
    if (!svg) return null;
    const resultIndex = Number(state.selectedResultIndex?.value || 0);
    if (!activeRuntime || activeRuntime.svg !== svg || activeRuntime.resultIndex !== resultIndex) {
      return mountResultSvg(resultIndex, svg);
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
    releaseActiveResult();
    activeRuntime = null;
    state.selectedResultIndex.value = nextIndex;
    return true;
  };

  const buildFeatureIndex = (runtime) => {
    const indexed = runtime?.svg ? getFeatureElementIndex(runtime.svg) : new Map();
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
      filterFeatureFillTargets(getFeatureElements(change.featureId)).forEach((element) => {
        if (element.getAttribute?.('fill') === color) return;
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
          if (element.getAttribute?.('display') === 'none') return;
          element.setAttribute('display', 'none');
        } else {
          if (element.getAttribute?.('display') === null) return;
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
        let changed = false;
        if (strokeColor && element.getAttribute?.('stroke') !== strokeColor) {
          element.setAttribute('stroke', strokeColor);
          changed = true;
        }
        const normalizedWidth = hasStrokeWidth ? String(Number(strokeWidth)) : '';
        if (hasStrokeWidth && element.getAttribute?.('stroke-width') !== normalizedWidth) {
          element.setAttribute('stroke-width', normalizedWidth);
          changed = true;
        }
        if (changed) updated += 1;
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
