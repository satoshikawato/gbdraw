import {
  filterFeatureFillTargets,
  getFeatureElementIndex,
  normalizeFeatureIdentity
} from './feature-dom.js';
import {
  markCommittedSvgResultMounted,
  markCommittedSvgResultUnmounted
} from '../services/svg-result-ingestion.js';
import { recordStructuralMetric } from '../services/runtime-test-hooks.js';
import { createDomMutationJournal } from './dom-mutation-journal.js';

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

const makeRuntime = (resultIndex, svg, resultOwner, mountRevision) => ({
  resultIndex,
  svg,
  resultOwner,
  mountRevision,
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
 * @property {object[]} resultsOwner Mounted Result collection identity.
 * @property {number} resultIndex Mounted Result index.
 * @property {number} mountRevision Mounted SVG ownership revision.
 * @property {string} generationKey Generation identity at target resolution.
 * @property {ReturnType<typeof createDomMutationJournal>} mutation Compact mutation journal.
 */

export const createPreviewRuntime = ({ state, serializeSvg }) => {
  if (!state) throw new Error('createPreviewRuntime requires state.');
  if (typeof serializeSvg !== 'function') throw new Error('createPreviewRuntime requires serializeSvg.');

  let activeRuntime = null;
  let activeDomEditBatch = null;
  let mountRevision = 0;

  const getMountedSvg = () => state.svgContainer?.value?.querySelector?.('svg') || null;

  const releaseActiveResult = () => {
    if (!activeRuntime) return;
    markCommittedSvgResultUnmounted(activeRuntime.resultOwner);
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
      && activeRuntime.resultOwner === state.results.value[normalizedIndex]
    ) {
      return activeRuntime;
    }
    releaseActiveResult();
    mountRevision += 1;
    const resultOwner = state.results.value[normalizedIndex];
    activeRuntime = makeRuntime(normalizedIndex, svg, resultOwner, mountRevision);
    markCommittedSvgResultMounted(resultOwner);
    return activeRuntime;
  };

  const clearActiveRuntime = () => {
    releaseActiveResult();
    mountRevision += 1;
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
    if (activeRuntime.resultOwner !== state.results.value[resultIndex]) return null;
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
      resultsOwner: results,
      resultIndex,
      mountRevision: runtime.mountRevision,
      generationKey: String(state.resultGenerationKey?.value ?? '')
    };
  };

  const isCurrentCommitTarget = (target) => Boolean(
    target
    && activeRuntime === target.runtime
    && target.runtime?.mountRevision === target.mountRevision
    && getMountedSvg() === target.svg
    && Number(state.selectedResultIndex?.value || 0) === target.resultIndex
    && state.results.value === target.resultsOwner
    && target.resultsOwner[target.resultIndex] === target.result
    && String(state.resultGenerationKey?.value ?? '') === target.generationKey
  );

  const rejectStaleTarget = (target, reason = 'preview-edit') => {
    recordStructuralMetric('staleDomEditTargetRejectionCount', 1, {
      reason: String(reason || 'preview-edit'),
      resultIndex: Number(target?.resultIndex ?? -1)
    });
    if (target?.runtime) {
      target.runtime.dirty = false;
      target.runtime.dirtyReasons.clear();
    }
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
    if (!isCurrentCommitTarget(target)) {
      rejectStaleTarget(target, target.runtime.dirtyReasons.values().next().value);
      return false;
    }

    // Serialization may throw. Keep dirty state intact so the caller can retry or
    // report the failure without silently treating an uncommitted edit as saved.
    const content = serializeSvg(target.svg);
    if (!isCurrentCommitTarget(target)) {
      rejectStaleTarget(target, target.runtime.dirtyReasons.values().next().value);
      return false;
    }
    if (target.result.content === content) {
      target.runtime.dirty = false;
      target.runtime.dirtyReasons.clear();
      return false;
    }
    if (markIncremental && state.skipCaptureBaseConfig) state.skipCaptureBaseConfig.value = true;
    const nextResult = {
      ...target.result,
      content
    };
    state.results.value[target.resultIndex] = nextResult;
    recordStructuralMetric('domEditResultReplacementCount', 1, {
      resultIndex: target.resultIndex,
      resultOwnerBefore: target.result,
      resultOwnerAfter: nextResult,
      resultsOwner: target.resultsOwner,
      svgOwner: target.svg
    });
    target.result = nextResult;
    target.runtime.resultOwner = nextResult;
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
   * Mutations made through target.mutation are rolled back if the callback
   * throws. Callers that write directly must validate every target and value
   * before their first write. A serialization error leaves the runtime dirty so
   * the unchanged Result artifact cannot be mistaken for the edit and can be
   * retried with flushActiveResult().
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
    if (!isCurrentCommitTarget(target)) {
      rejectStaleTarget(target, reason);
      return { changed: false, flushed: false, resultIndex: target.resultIndex, reason };
    }
    const mutation = createDomMutationJournal();
    let outcome;
    try {
      outcome = mutate({ ...target, mutation });
    } catch (error) {
      mutation.rollback();
      throw error;
    }
    if (outcome && typeof outcome.then === 'function') {
      mutation.rollback();
      throw new TypeError('commitDomEdit mutate callback must be synchronous.');
    }
    const changed = mutationChanged(outcome) || mutation.changed;
    if (!changed) {
      mutation.commit();
      return { changed: false, flushed: false, resultIndex: target.resultIndex, reason };
    }
    if (!isCurrentCommitTarget(target)) {
      mutation.rollback();
      rejectStaleTarget(target, reason);
      return { changed: false, flushed: false, resultIndex: target.resultIndex, reason };
    }

    recordDomChange(target, reason, invalidateIndexes);
    if (activeDomEditBatch) {
      activeDomEditBatch.changed = true;
      activeDomEditBatch.mutations.push(mutation);
      return { changed: true, flushed: false, resultIndex: target.resultIndex, reason };
    }
    try {
      const result = {
        changed: true,
        flushed: flushCommitTarget(target),
        resultIndex: target.resultIndex,
        reason
      };
      recordStructuralMetric('domEditCommitCount', 1, {
        reason: String(reason || 'preview-edit'),
        resultIndex: target.resultIndex,
        resultOwner: target.result,
        resultsOwner: target.resultsOwner,
        svgOwner: target.svg
      });
      return result;
    } finally {
      mutation.commit();
    }
  };

  /**
   * Hold an explicit synchronous-event lease for pointer drags. The lease owns
   * one mounted SVG/Result pair, revalidates it before every mutation and before
   * settlement, and serializes once when the interaction ends.
   */
  const beginDomEditLease = ({
    reason = 'preview-edit',
    invalidateIndexes = []
  } = {}) => {
    const target = resolveActiveCommitTarget();
    if (!target) return null;
    const mutation = createDomMutationJournal();
    let changed = false;
    let closed = false;

    const closeStale = () => {
      mutation.rollback();
      rejectStaleTarget(target, reason);
      closed = true;
      return false;
    };

    return {
      get current() {
        return !closed && isCurrentCommitTarget(target);
      },
      get target() {
        return target;
      },
      mutate(callback) {
        if (closed) return false;
        if (!isCurrentCommitTarget(target)) return closeStale();
        if (typeof callback !== 'function') throw new TypeError('DOM edit lease requires a mutate callback.');
        let outcome;
        try {
          outcome = callback({ ...target, mutation });
        } catch (error) {
          mutation.rollback();
          closed = true;
          throw error;
        }
        if (outcome && typeof outcome.then === 'function') {
          mutation.rollback();
          closed = true;
          throw new TypeError('DOM edit lease mutate callback must be synchronous.');
        }
        changed = mutationChanged(outcome) || mutation.changed || changed;
        if (!isCurrentCommitTarget(target)) return closeStale();
        return mutationChanged(outcome) || mutation.changed;
      },
      commit() {
        if (closed) return false;
        if (!isCurrentCommitTarget(target)) return closeStale();
        closed = true;
        if (!changed) {
          mutation.commit();
          return false;
        }
        recordDomChange(target, reason, invalidateIndexes);
        try {
          const flushed = flushCommitTarget(target);
          recordStructuralMetric('domEditLeaseCommitCount', 1, {
            reason: String(reason || 'preview-edit'),
            resultIndex: target.resultIndex,
            resultOwner: target.result,
            resultsOwner: target.resultsOwner,
            svgOwner: target.svg
          });
          return flushed;
        } finally {
          mutation.commit();
        }
      },
      cancel() {
        if (closed) return false;
        closed = true;
        const hadChanges = changed || mutation.changed;
        mutation.rollback();
        return hadChanges;
      }
    };
  };

  /**
   * Batch only the synchronous prefix of a user action. If the action returns a
   * Promise, direct edits are committed before the first asynchronous settlement;
   * later continuation edits resolve a new current target. Mounted ownership is
   * never retained across await.
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
    const runtimeState = {
      dirty: target.runtime.dirty,
      dirtyReasons: new Set(target.runtime.dirtyReasons),
      indexes: { ...target.runtime.indexes },
      lastInvalidationReason: target.runtime.lastInvalidationReason
    };
    activeDomEditBatch = { target, changed: false, reason, mutations: [] };
    let result;
    try {
      result = action();
    } catch (error) {
      for (let index = activeDomEditBatch.mutations.length - 1; index >= 0; index -= 1) {
        activeDomEditBatch.mutations[index].rollback();
      }
      target.runtime.dirty = runtimeState.dirty;
      target.runtime.dirtyReasons = runtimeState.dirtyReasons;
      target.runtime.indexes = runtimeState.indexes;
      target.runtime.lastInvalidationReason = runtimeState.lastInvalidationReason;
      activeDomEditBatch = null;
      throw error;
    }
    const completedBatch = activeDomEditBatch;
    activeDomEditBatch = null;
    try {
      if (completedBatch.changed) flushCommitTarget(target);
    } finally {
      completedBatch.mutations.forEach((mutation) => mutation.commit());
    }
    return await result;
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
      mutate: ({ mutation }) => {
        let updated = 0;
        normalized.forEach((change) => {
          const color = String(change?.color || '').trim();
          if (!color) return;
          filterFeatureFillTargets(getFeatureElements(change.featureId)).forEach((element) => {
            if (element.getAttribute?.('fill') === color) return;
            mutation.setAttribute(element, 'fill', color);
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
      mutate: ({ mutation }) => {
        let updated = 0;
        normalized.forEach((change) => {
          const mode = normalizeVisibilityMode(change?.mode);
          getFeatureElements(change.featureId).forEach((element) => {
            if (mode === 'off') {
              if (element.getAttribute?.('display') === 'none') return;
              mutation.setAttribute(element, 'display', 'none');
            } else {
              if (element.getAttribute?.('display') === null) return;
              mutation.removeAttribute(element, 'display');
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
      mutate: ({ mutation }) => {
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
              if (strokeColor === null) mutation.removeAttribute(element, 'stroke');
              else mutation.setAttribute(element, 'stroke', strokeColor);
              changed = true;
            }
            if (hasStrokeWidth && element.getAttribute?.('stroke-width') !== normalizedWidth) {
              if (normalizedWidth === null) mutation.removeAttribute(element, 'stroke-width');
              else mutation.setAttribute(element, 'stroke-width', normalizedWidth);
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
    beginDomEditLease,
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
