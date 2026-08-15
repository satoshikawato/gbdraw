import {
  filterFeatureFillTargets,
  getFeatureElementIndex,
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

/**
 * @typedef {Object} DomEditCommitResult
 * @property {boolean} changed Whether the mutation changed the mounted SVG.
 * @property {boolean} flushed Whether the serialized Result artifact changed.
 * @property {number|null} resultIndex Index of the mounted Result, when available.
 * @property {string} reason Stable diagnostic reason for the edit.
 */

/**
 * @typedef {Object} DomEditTarget
 * @property {object} runtime Mounted preview runtime.
 * @property {object} svg Mounted SVG element.
 * @property {object} result Mounted Result artifact.
 * @property {number} resultIndex Mounted Result index.
 */

export const createPreviewRuntime = ({ state, serializeSvg }) => {
  if (!state) throw new Error('createPreviewRuntime requires state.');
  if (typeof serializeSvg !== 'function') throw new Error('createPreviewRuntime requires serializeSvg.');

  let activeRuntime = null;
  let activeDomEditBatch = null;

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

  const resolveActiveCommitTarget = () => {
    const runtime = activeRuntime || ensureRuntimeForCurrentSvg();
    if (!runtime?.svg) return null;
    const resultIndex = Number(runtime.resultIndex);
    const results = state.results.value;
    if (
      !Array.isArray(results)
      || !Number.isInteger(resultIndex)
      || resultIndex < 0
      || resultIndex >= results.length
      || !results[resultIndex]
    ) {
      return null;
    }
    return {
      runtime,
      svg: runtime.svg,
      result: results[resultIndex],
      resultIndex
    };
  };

  const normalizeInvalidationKeys = (keys) => (
    Array.isArray(keys)
      ? Array.from(new Set(keys.filter((key) => typeof key === 'string' && key)))
      : []
  );

  const invalidateRuntimeIndexes = (runtime, reason, keys) => {
    normalizeInvalidationKeys(keys).forEach((key) => {
      if (Object.prototype.hasOwnProperty.call(runtime.indexes, key)) {
        runtime.indexes[key] = null;
      }
    });
    runtime.lastInvalidationReason = String(reason || 'unknown');
  };

  const recordDomChange = (target, reason, invalidateIndexes) => {
    target.runtime.dirty = true;
    target.runtime.dirtyReasons.add(String(reason || 'preview-edit'));
    invalidateRuntimeIndexes(target.runtime, reason, invalidateIndexes);
  };

  const flushCommitTarget = (target, { markIncremental = true } = {}) => {
    if (!target?.runtime?.dirty) return false;

    // Serialization may throw. Keep dirty state intact so the caller can retry or
    // report the failure without silently treating an uncommitted edit as saved.
    const content = serializeSvg(target.svg);
    if (target.result.content === content) {
      target.runtime.dirty = false;
      target.runtime.dirtyReasons.clear();
      return false;
    }
    if (markIncremental && state.skipCaptureBaseConfig) state.skipCaptureBaseConfig.value = true;
    state.results.value[target.resultIndex] = {
      ...target.result,
      content
    };
    target.runtime.dirty = false;
    target.runtime.dirtyReasons.clear();
    return true;
  };

  const mutationChanged = (value) => {
    if (value && typeof value === 'object' && Object.hasOwn(value, 'changed')) {
      return Boolean(value.changed);
    }
    return value !== false && value !== 0 && value !== null && value !== undefined;
  };

  /**
   * Sole production owner for applying a mounted-SVG edit to its selected Result.
   * It never invokes a Worker, Python helper, Generate, or History checkpoint.
   *
   * Apply one synchronous mutation to the mounted SVG and commit its Result.
   * The callback must not return a Promise. A false, zero, null, or undefined
   * result is a no-op and causes no serialization or Result replacement.
   * A mutation error commits nothing. A serialization error leaves the runtime
   * dirty so the unchanged Result artifact cannot be mistaken for the edit.
   *
   * @param {Object} options
   * @param {string} [options.reason]
   * @param {(target: DomEditTarget) => (boolean|number|Object)} options.mutate
   * @param {string[]} [options.invalidateIndexes]
   * @returns {DomEditCommitResult}
   */
  const commitDomEdit = ({
    reason = 'preview-edit',
    mutate,
    invalidateIndexes = []
  } = {}) => {
    if (typeof mutate !== 'function') {
      throw new TypeError('commitDomEdit requires a synchronous mutate callback.');
    }
    const target = activeDomEditBatch?.target || resolveActiveCommitTarget();
    if (!target) {
      return { changed: false, flushed: false, resultIndex: null, reason };
    }

    const outcome = mutate(target);
    if (outcome && typeof outcome.then === 'function') {
      throw new TypeError('commitDomEdit mutate callback must be synchronous.');
    }
    const changed = mutationChanged(outcome);
    if (!changed) {
      return { changed: false, flushed: false, resultIndex: target.resultIndex, reason };
    }

    recordDomChange(target, reason, invalidateIndexes);
    if (activeDomEditBatch) {
      activeDomEditBatch.changed = true;
      return { changed: true, flushed: false, resultIndex: target.resultIndex, reason };
    }
    return {
      changed: true,
      flushed: flushCommitTarget(target),
      resultIndex: target.resultIndex,
      reason
    };
  };

  /**
   * Batch nested synchronous DOM commits across one asynchronous user action.
   * Nested commitDomEdit calls share one resolved target and serialize at most
   * once after the action completes successfully.
   * If the action fails, no partial Result replacement occurs and dirty state is
   * retained for explicit recovery.
   *
   * @template T
   * @param {Object} options
   * @param {string} [options.reason]
   * @param {() => (T|Promise<T>)} options.action
   * @returns {Promise<T>}
   */
  const runDomEdit = async ({ reason = 'preview-edit', action } = {}) => {
    if (typeof action !== 'function') throw new TypeError('runDomEdit requires an action callback.');
    if (activeDomEditBatch) return action();

    const target = resolveActiveCommitTarget();
    if (!target) return action();
    activeDomEditBatch = { target, changed: false, reason };
    try {
      const result = await action();
      if (activeDomEditBatch.changed) flushCommitTarget(target);
      return result;
    } finally {
      activeDomEditBatch = null;
    }
  };

  const invalidatePreviewIndexes = (reason = 'unknown', keys = null) => {
    const runtime = activeRuntime;
    if (!runtime) return;
    const targetKeys = Array.isArray(keys) && keys.length > 0
      ? keys
      : Object.keys(runtime.indexes);
    invalidateRuntimeIndexes(runtime, reason, targetKeys);
  };

  const flushActiveResult = ({ force = false, markIncremental = true } = {}) => {
    const target = resolveActiveCommitTarget();
    if (!target) return false;
    if (!force && !target.runtime.dirty) return false;
    if (force && !target.runtime.dirty) target.runtime.dirty = true;
    return flushCommitTarget(target, { markIncremental });
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
    return commitDomEdit({
      reason,
      invalidateIndexes: ['features'],
      mutate: () => {
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
        return updated;
      }
    }).changed;
  };

  const applyFeatureVisibilityChanges = (changes, { reason = 'feature-visibility' } = {}) => {
    const normalized = normalizeChanges(changes);
    if (normalized.length === 0) return false;
    return commitDomEdit({
      reason,
      invalidateIndexes: ['features'],
      mutate: () => {
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
        return updated;
      }
    }).changed;
  };

  const applyFeatureStrokeChanges = (changes, { reason = 'feature-stroke' } = {}) => {
    const normalized = normalizeChanges(changes);
    if (normalized.length === 0) return false;
    return commitDomEdit({
      reason,
      invalidateIndexes: ['features'],
      mutate: () => {
        let updated = 0;
        normalized.forEach((change) => {
          const hasStrokeColor = Object.hasOwn(change, 'strokeColor');
          const strokeColor = change.strokeColor === null || change.strokeColor === undefined || change.strokeColor === ''
            ? null
            : String(change.strokeColor || '').trim();
          const strokeWidth = change?.strokeWidth;
          const hasStrokeWidth = Object.hasOwn(change, 'strokeWidth');
          const normalizedWidth = strokeWidth === null || strokeWidth === undefined || strokeWidth === ''
            ? null
            : String(Number(strokeWidth));
          getFeatureElements(change.featureId).forEach((element) => {
            let changed = false;
            if (hasStrokeColor && element.getAttribute?.('stroke') !== strokeColor) {
              if (strokeColor === null) element.removeAttribute('stroke');
              else element.setAttribute('stroke', strokeColor);
              changed = true;
            }
            if (hasStrokeWidth && element.getAttribute?.('stroke-width') !== normalizedWidth) {
              if (normalizedWidth === null) element.removeAttribute('stroke-width');
              else element.setAttribute('stroke-width', normalizedWidth);
              changed = true;
            }
            if (changed) updated += 1;
          });
        });
        return updated;
      }
    }).changed;
  };

  return {
    applyFeatureFillChanges,
    applyFeatureStrokeChanges,
    applyFeatureVisibilityChanges,
    clearActiveRuntime,
    commitDomEdit,
    flushActiveResult,
    getActiveRuntime,
    getFeatureElements,
    invalidatePreviewIndexes,
    mountResultSvg,
    runDomEdit,
    selectResult
  };
};
