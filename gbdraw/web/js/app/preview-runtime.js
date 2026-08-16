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
import { applyFeatureVisibility } from './feature-visibility-dom.js';

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

const newDocumentAdmissionIdentity = () => Object.freeze({});

const makeRuntime = (
  resultIndex,
  svg,
  resultOwner,
  mountedResultContent,
  mountRevision,
  documentAdmissionIdentity = newDocumentAdmissionIdentity()
) => ({
  resultIndex,
  svg,
  resultOwner,
  mountedResultContent,
  mountRevision,
  documentAdmissionIdentity,
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
 * @property {*} resultContent Immutable Result content captured with the target.
 * @property {object[]} resultsOwner Mounted Result collection identity.
 * @property {number} resultIndex Mounted Result index.
 * @property {number} mountRevision Mounted SVG ownership revision.
 * @property {object} documentAdmissionIdentity Mounted document admission identity.
 * @property {string} generationKey Generation identity at target resolution.
 * @property {ReturnType<typeof createDomMutationJournal>} mutation Compact mutation journal.
 */

export const createPreviewRuntime = ({ state, serializeSvg }) => {
  if (!state) throw new Error('createPreviewRuntime requires state.');
  if (typeof serializeSvg !== 'function') throw new Error('createPreviewRuntime requires serializeSvg.');

  let activeRuntime = null;
  let activeDomEditBatch = null;
  let activeDomEditLease = null;
  let mountRevision = 0;

  const liveEditBusyError = () => {
    const error = new Error('Another live edit transaction is still in progress.');
    error.name = 'LiveEditBusyError';
    return error;
  };

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
      activeRuntime.documentAdmissionIdentity = newDocumentAdmissionIdentity();
      activeRuntime.mountedResultContent = activeRuntime.resultOwner?.content;
      return activeRuntime;
    }
    releaseActiveResult();
    mountRevision += 1;
    const resultOwner = state.results.value[normalizedIndex];
    activeRuntime = makeRuntime(
      normalizedIndex,
      svg,
      resultOwner,
      resultOwner?.content,
      mountRevision
    );
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
    if (
      activeRuntime.resultOwner !== state.results.value[resultIndex]
      || activeRuntime.mountedResultContent !== state.results.value[resultIndex]?.content
    ) return null;
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
      || runtime.resultOwner !== results[resultIndex]
      || runtime.mountedResultContent !== results[resultIndex].content
    ) {
      return null;
    }
    return {
      runtime,
      svg: runtime.svg,
      result: results[resultIndex],
      resultContent: results[resultIndex].content,
      resultsOwner: results,
      resultIndex,
      mountRevision: runtime.mountRevision,
      documentAdmissionIdentity: runtime.documentAdmissionIdentity,
      generationKey: String(state.resultGenerationKey?.value ?? '')
    };
  };

  const isCurrentCommitTarget = (target) => Boolean(
    target
    && activeRuntime === target.runtime
    && target.runtime?.resultOwner === target.result
    && target.runtime?.mountedResultContent === target.resultContent
    && target.runtime?.mountRevision === target.mountRevision
    && target.runtime?.documentAdmissionIdentity === target.documentAdmissionIdentity
    && getMountedSvg() === target.svg
    && Number(state.selectedResultIndex?.value || 0) === target.resultIndex
    && state.results.value === target.resultsOwner
    && target.resultsOwner[target.resultIndex] === target.result
    && target.result.content === target.resultContent
    && String(state.resultGenerationKey?.value ?? '') === target.generationKey
  );

  const captureDomEditToken = () => {
    const target = resolveActiveCommitTarget();
    return isCurrentCommitTarget(target) ? Object.freeze({ ...target }) : null;
  };

  const isDomEditTokenCurrent = (targetToken) => isCurrentCommitTarget(targetToken);

  const recordTargetRejection = (target, reason = 'preview-edit') => {
    const staleContent = Boolean(
      target?.result
      && Object.prototype.hasOwnProperty.call(target, 'resultContent')
      && target.result.content !== target.resultContent
    );
    recordStructuralMetric('staleDomEditTargetRejectionCount', 1, {
      reason: String(reason || 'preview-edit'),
      resultIndex: Number(target?.resultIndex ?? -1),
      staleContent
    });
    if (staleContent) {
      recordStructuralMetric('liveEditStaleContentTokenCount', 1, {
        reason: String(reason || 'preview-edit'),
        resultIndex: Number(target?.resultIndex ?? -1)
      });
    }
  };

  const rejectStaleTarget = (target, reason = 'preview-edit') => {
    recordTargetRejection(target, reason);
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

  const captureRuntimeState = (runtime) => ({
    dirty: runtime.dirty,
    dirtyReasons: new Set(runtime.dirtyReasons),
    indexes: { ...runtime.indexes },
    lastInvalidationReason: runtime.lastInvalidationReason,
    resultOwner: runtime.resultOwner,
    mountedResultContent: runtime.mountedResultContent
  });

  const restoreRuntimeState = (runtime, snapshot) => {
    runtime.dirty = snapshot.dirty;
    runtime.dirtyReasons = new Set(snapshot.dirtyReasons);
    runtime.indexes = { ...snapshot.indexes };
    runtime.lastInvalidationReason = snapshot.lastInvalidationReason;
    runtime.resultOwner = snapshot.resultOwner;
    runtime.mountedResultContent = snapshot.mountedResultContent;
  };

  const captureCommitState = (target) => ({
    resultOwner: target.resultsOwner[target.resultIndex],
    runtime: captureRuntimeState(target.runtime),
    skipCaptureBaseConfig: state.skipCaptureBaseConfig?.value
  });

  const captureCanonicalStateTargets = (mutation, targets = []) => {
    (Array.isArray(targets) ? targets : []).forEach((descriptor) => {
      const target = descriptor?.target;
      if (target === null || target === undefined) return;
      if (Object.prototype.hasOwnProperty.call(descriptor, 'key')) {
        mutation.captureProperty(target, descriptor.key, { deep: Boolean(descriptor.deep) });
      } else {
        mutation.captureState(target, { deep: descriptor?.deep !== false });
      }
    });
  };

  const captureCanonicalState = (targets = []) => {
    const mutation = createDomMutationJournal();
    captureCanonicalStateTargets(mutation, targets);
    return mutation;
  };

  const attachRollbackErrors = (primaryError, errors) => {
    if (!Array.isArray(errors) || errors.length === 0) return primaryError;
    let error = primaryError instanceof Error
      ? primaryError
      : new Error(String(primaryError), { cause: primaryError });
    const secondaryErrors = Object.freeze([
      ...(Array.isArray(error.rollbackErrors) ? error.rollbackErrors : []),
      ...errors
    ]);
    try {
      Object.defineProperty(error, 'rollbackErrors', {
        configurable: true,
        enumerable: false,
        value: secondaryErrors
      });
    } catch (_attachError) {
      const wrapped = new Error(error.message, { cause: error });
      wrapped.name = error.name;
      Object.defineProperty(wrapped, 'rollbackErrors', {
        enumerable: false,
        value: secondaryErrors
      });
      error = wrapped;
    }
    return error;
  };

  // A Result collection setter can install the replacement and then throw. Keep
  // that otherwise-unobservable replacement associated with the primary error
  // so every outer transaction owner can retry the Result restore after it has
  // rolled back the independent DOM, canonical, and runtime state.
  const failedResultReplacements = new WeakMap();
  const failedResultReplacement = (error, fallback = null) => (
    (error && typeof error === 'object' && failedResultReplacements.get(error)) || fallback
  );

  const recordMetricSafely = (name, value, detail, errors) => {
    try {
      recordStructuralMetric(name, value, detail);
    } catch (error) {
      errors.push(error instanceof Error ? error : new Error(String(error)));
    }
  };

  const rollbackJournal = (mutation, { reason, resultIndex } = {}) => {
    if (!mutation) return [];
    let errors;
    try {
      errors = [...mutation.rollback()];
    } catch (error) {
      errors = [error];
    }
    const detail = {
      reason: String(reason || 'preview-edit'),
      resultIndex: Number(resultIndex ?? -1)
    };
    const rollbackFailureCount = errors.length;
    if (mutation.domChangeCount > 0) {
      recordMetricSafely('liveEditDomRollbackCount', 1, detail, errors);
    }
    if (mutation.stateChangeCount > 0 || mutation.capturedStateOwnerCount > 0) {
      recordMetricSafely('liveEditCanonicalRollbackCount', 1, detail, errors);
      if (rollbackFailureCount > 0) {
        recordMetricSafely('liveEditCanonicalRollbackFailureCount', rollbackFailureCount, {
          ...detail,
          errorCount: rollbackFailureCount
        }, errors);
      }
    }
    return errors;
  };

  const rollbackJournals = (mutations, canonicalMutation, detail) => {
    const errors = [];
    for (let index = mutations.length - 1; index >= 0; index -= 1) {
      errors.push(...rollbackJournal(mutations[index], detail));
    }
    errors.push(...rollbackJournal(canonicalMutation, detail));
    return errors;
  };

  const restoreCommitState = (target, snapshot, replacement = null, reason = 'preview-edit') => {
    const errors = [];
    let resultRestored = false;
    try {
      if (
        replacement
        && target.resultsOwner[target.resultIndex] === replacement
        && state.results.value === target.resultsOwner
      ) {
        target.resultsOwner[target.resultIndex] = snapshot.resultOwner;
        resultRestored = true;
      }
    } catch (error) {
      errors.push(error instanceof Error ? error : new Error(String(error)));
    }
    const runtimeRestorations = [
      () => { target.runtime.dirty = snapshot.runtime.dirty; },
      () => { target.runtime.dirtyReasons = new Set(snapshot.runtime.dirtyReasons); },
      () => { target.runtime.indexes = { ...snapshot.runtime.indexes }; },
      () => { target.runtime.lastInvalidationReason = snapshot.runtime.lastInvalidationReason; },
      () => { target.runtime.resultOwner = snapshot.runtime.resultOwner; },
      () => { target.runtime.mountedResultContent = snapshot.runtime.mountedResultContent; }
    ];
    runtimeRestorations.forEach((restore) => {
      try {
        restore();
      } catch (error) {
        errors.push(error instanceof Error ? error : new Error(String(error)));
      }
    });
    try {
      if (state.skipCaptureBaseConfig) {
        state.skipCaptureBaseConfig.value = snapshot.skipCaptureBaseConfig;
      }
    } catch (error) {
      errors.push(error instanceof Error ? error : new Error(String(error)));
    }
    recordMetricSafely('liveEditResultRestoreCount', 1, {
      reason: String(reason || 'preview-edit'),
      resultIndex: target.resultIndex,
      resultRestored
    }, errors);
    return errors;
  };

  const restoreCommitStateSafely = (target, snapshot, replacement, reason, errors) => {
    errors.push(...restoreCommitState(target, snapshot, replacement, reason));
  };

  const staleRollbackError = (reason) => {
    const error = new Error(`Live edit target became stale during ${String(reason || 'preview-edit')}.`);
    error.name = 'StaleDomEditError';
    return error;
  };

  const flushCommitTarget = (target, { markIncremental = true } = {}) => {
    if (!target?.runtime?.dirty) return { flushed: false, replacement: null, stale: false };
    if (!isCurrentCommitTarget(target)) {
      rejectStaleTarget(target, target.runtime.dirtyReasons.values().next().value);
      return { flushed: false, replacement: null, stale: true };
    }

    // Serialization may throw. Keep dirty state intact so the caller can retry or
    // report the failure without silently treating an uncommitted edit as saved.
    const content = serializeSvg(target.svg);
    if (!isCurrentCommitTarget(target)) {
      rejectStaleTarget(target, target.runtime.dirtyReasons.values().next().value);
      return { flushed: false, replacement: null, stale: true };
    }
    if (target.result.content === content) {
      target.runtime.dirty = false;
      target.runtime.dirtyReasons.clear();
      return { flushed: false, replacement: null, stale: false };
    }
    if (markIncremental && state.skipCaptureBaseConfig) state.skipCaptureBaseConfig.value = true;
    const nextResult = {
      ...target.result,
      content
    };
    if (!isCurrentCommitTarget(target)) {
      rejectStaleTarget(target, target.runtime.dirtyReasons.values().next().value);
      return { flushed: false, replacement: null, stale: true };
    }
    try {
      state.results.value[target.resultIndex] = nextResult;
      // Vue deep refs expose a reactive proxy rather than the raw object that was
      // assigned. The installed owner is the identity every later token and
      // rollback comparison must use.
      const installedReplacement = target.resultsOwner[target.resultIndex];
      recordStructuralMetric('domEditResultReplacementCount', 1, {
        resultIndex: target.resultIndex,
        resultOwnerBefore: target.result,
        resultOwnerAfter: installedReplacement,
        resultsOwner: target.resultsOwner,
        svgOwner: target.svg
      });
      target.runtime.resultOwner = installedReplacement;
      target.runtime.mountedResultContent = installedReplacement?.content;
      target.runtime.dirty = false;
      target.runtime.dirtyReasons.clear();
      return { flushed: true, replacement: installedReplacement, stale: false };
    } catch (error) {
      const rollbackErrors = [];
      const installedReplacement = state.results.value === target.resultsOwner
        && target.resultsOwner[target.resultIndex] !== target.result
        ? target.resultsOwner[target.resultIndex]
        : null;
      if (installedReplacement) {
        try {
          target.resultsOwner[target.resultIndex] = target.result;
        } catch (rollbackError) {
          rollbackErrors.push(
            rollbackError instanceof Error ? rollbackError : new Error(String(rollbackError))
          );
        }
      }
      const primaryError = attachRollbackErrors(error, rollbackErrors);
      if (installedReplacement) failedResultReplacements.set(primaryError, installedReplacement);
      throw primaryError;
    }
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
   * Mutations made through target.mutation are rolled back if mutation,
   * serialization, Result replacement, or metric recording fails. The rollback
   * restores the DOM, Result, runtime dirty state, and invalidated indexes.
   *
   * @param {Object} options
   * @param {string} [options.reason]
   * @param {(target: DomEditTarget) => (boolean|number|Object)} options.mutate
   * @param {string[]} [options.invalidateIndexes]
   * @param {boolean} [options.journalChangesAffectResult] Whether journal writes
   * imply a serialized Result change. Editor-only annotations set this false.
   * @param {DomEditTarget|null} [options.targetToken]
   * @returns {DomEditCommitResult}
   */
  const commitDomEdit = ({
    reason = 'preview-edit',
    mutate,
    invalidateIndexes = [],
    journalChangesAffectResult = true,
    targetToken = null
  } = {}) => {
    if (typeof mutate !== 'function') {
      throw new TypeError('commitDomEdit requires a synchronous mutate callback.');
    }
    if (activeDomEditLease) throw liveEditBusyError();
    const target = targetToken || activeDomEditBatch?.target || resolveActiveCommitTarget();
    if (!target) {
      return { changed: false, flushed: false, resultIndex: null, reason };
    }
    if (!isCurrentCommitTarget(target)) {
      rejectStaleTarget(target, reason);
      if (activeDomEditBatch?.target === target) activeDomEditBatch.stale = true;
      return { changed: false, flushed: false, resultIndex: target.resultIndex, reason };
    }
    if (activeDomEditBatch && target !== activeDomEditBatch.target) {
      recordTargetRejection(target, reason);
      return { changed: false, flushed: false, resultIndex: target.resultIndex, reason };
    }
    const commitState = activeDomEditBatch ? null : captureCommitState(target);
    const mutation = createDomMutationJournal();
    let outcome;
    try {
      outcome = mutate({ ...target, mutation });
    } catch (error) {
      const rollbackErrors = rollbackJournal(mutation, { reason, resultIndex: target.resultIndex });
      if (commitState) {
        restoreCommitStateSafely(target, commitState, null, reason, rollbackErrors);
      }
      throw attachRollbackErrors(error, rollbackErrors);
    }
    if (outcome && typeof outcome.then === 'function') {
      const error = new TypeError('commitDomEdit mutate callback must be synchronous.');
      const rollbackErrors = rollbackJournal(mutation, { reason, resultIndex: target.resultIndex });
      if (commitState) {
        restoreCommitStateSafely(target, commitState, null, reason, rollbackErrors);
      }
      throw attachRollbackErrors(error, rollbackErrors);
    }
    const changed = mutationChanged(outcome) || (journalChangesAffectResult && mutation.changed);
    if (!isCurrentCommitTarget(target)) {
      const rollbackErrors = rollbackJournal(mutation, { reason, resultIndex: target.resultIndex });
      if (activeDomEditBatch) activeDomEditBatch.stale = true;
      else restoreCommitStateSafely(target, commitState, null, reason, rollbackErrors);
      try {
        rejectStaleTarget(target, reason);
      } catch (error) {
        rollbackErrors.push(error instanceof Error ? error : new Error(String(error)));
      }
      if (rollbackErrors.length > 0) {
        throw attachRollbackErrors(staleRollbackError(reason), rollbackErrors);
      }
      return { changed: false, flushed: false, resultIndex: target.resultIndex, reason };
    }
    if (!changed) {
      if (activeDomEditBatch && mutation.changed) activeDomEditBatch.mutations.push(mutation);
      else mutation.commit();
      return { changed: false, flushed: false, resultIndex: target.resultIndex, reason };
    }

    recordDomChange(target, reason, invalidateIndexes);
    if (activeDomEditBatch) {
      activeDomEditBatch.changed = true;
      activeDomEditBatch.mutations.push(mutation);
      return { changed: true, flushed: false, resultIndex: target.resultIndex, reason };
    }
    let flushStatus = { flushed: false, replacement: null, stale: false };
    let rollbackComplete = false;
    try {
      flushStatus = flushCommitTarget(target);
      if (flushStatus.stale) {
        const rollbackErrors = rollbackJournal(mutation, { reason, resultIndex: target.resultIndex });
        restoreCommitStateSafely(target, commitState, null, reason, rollbackErrors);
        rollbackComplete = true;
        if (rollbackErrors.length > 0) {
          throw attachRollbackErrors(staleRollbackError(reason), rollbackErrors);
        }
        return { changed: false, flushed: false, resultIndex: target.resultIndex, reason };
      }
      const result = {
        changed: true,
        flushed: flushStatus.flushed,
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
      mutation.commit();
      return result;
    } catch (error) {
      if (rollbackComplete) throw error;
      const rollbackErrors = rollbackJournal(mutation, { reason, resultIndex: target.resultIndex });
      restoreCommitStateSafely(
        target,
        commitState,
        failedResultReplacement(error, flushStatus.replacement),
        reason,
        rollbackErrors
      );
      throw attachRollbackErrors(error, rollbackErrors);
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
    if (activeDomEditLease) throw liveEditBusyError();
    const target = resolveActiveCommitTarget();
    if (!target) return null;
    const commitState = captureCommitState(target);
    const mutation = createDomMutationJournal();
    let changed = false;
    let closed = false;
    let terminalError = null;
    let lease = null;
    const release = () => {
      if (activeDomEditLease === lease) activeDomEditLease = null;
    };

    const closeStale = () => {
      if (terminalError) throw terminalError;
      const rollbackErrors = rollbackJournal(mutation, { reason, resultIndex: target.resultIndex });
      restoreCommitStateSafely(target, commitState, null, reason, rollbackErrors);
      try {
        rejectStaleTarget(target, reason);
      } catch (error) {
        rollbackErrors.push(error instanceof Error ? error : new Error(String(error)));
      }
      closed = true;
      release();
      terminalError = attachRollbackErrors(staleRollbackError(reason), rollbackErrors);
      throw terminalError;
    };

    lease = {
      get current() {
        return !closed && isCurrentCommitTarget(target);
      },
      get target() {
        return target;
      },
      mutate(callback, { journalChangesAffectResult = true } = {}) {
        if (terminalError) throw terminalError;
        if (closed) return false;
        if (!isCurrentCommitTarget(target)) return closeStale();
        if (typeof callback !== 'function') throw new TypeError('DOM edit lease requires a mutate callback.');
        const domChangeCountBefore = mutation.domChangeCount;
        const stateChangeCountBefore = mutation.stateChangeCount;
        let outcome;
        try {
          outcome = callback({ ...target, mutation });
        } catch (error) {
          const rollbackErrors = rollbackJournal(mutation, { reason, resultIndex: target.resultIndex });
          restoreCommitStateSafely(target, commitState, null, reason, rollbackErrors);
          closed = true;
          release();
          throw attachRollbackErrors(error, rollbackErrors);
        }
        if (outcome && typeof outcome.then === 'function') {
          const error = new TypeError('DOM edit lease mutate callback must be synchronous.');
          const rollbackErrors = rollbackJournal(mutation, { reason, resultIndex: target.resultIndex });
          restoreCommitStateSafely(target, commitState, null, reason, rollbackErrors);
          closed = true;
          release();
          throw attachRollbackErrors(error, rollbackErrors);
        }
        const outcomeChanged = mutationChanged(outcome);
        const journalChanged = (
          mutation.domChangeCount > domChangeCountBefore
          || mutation.stateChangeCount > stateChangeCountBefore
        );
        changed = outcomeChanged || (journalChangesAffectResult && journalChanged) || changed;
        if (!isCurrentCommitTarget(target)) return closeStale();
        return outcomeChanged || journalChanged;
      },
      commit() {
        if (terminalError) throw terminalError;
        if (closed) return false;
        if (!isCurrentCommitTarget(target)) return closeStale();
        closed = true;
        try {
          if (!changed) {
            mutation.commit();
            return false;
          }
          recordDomChange(target, reason, invalidateIndexes);
          let flushStatus = { flushed: false, replacement: null, stale: false };
          let rollbackComplete = false;
          try {
            flushStatus = flushCommitTarget(target);
            if (flushStatus.stale) {
              const rollbackErrors = rollbackJournal(mutation, { reason, resultIndex: target.resultIndex });
              restoreCommitStateSafely(target, commitState, null, reason, rollbackErrors);
              rollbackComplete = true;
              throw attachRollbackErrors(staleRollbackError(reason), rollbackErrors);
            }
            recordStructuralMetric('domEditLeaseCommitCount', 1, {
              reason: String(reason || 'preview-edit'),
              resultIndex: target.resultIndex,
              resultOwner: target.result,
              resultsOwner: target.resultsOwner,
              svgOwner: target.svg
            });
            mutation.commit();
            return flushStatus.flushed;
          } catch (error) {
            if (rollbackComplete) throw error;
            const rollbackErrors = rollbackJournal(mutation, { reason, resultIndex: target.resultIndex });
            restoreCommitStateSafely(
              target,
              commitState,
              failedResultReplacement(error, flushStatus.replacement),
              reason,
              rollbackErrors
            );
            throw attachRollbackErrors(error, rollbackErrors);
          }
        } finally {
          release();
        }
      },
      cancel() {
        if (closed) return false;
        closed = true;
        try {
          const hadChanges = changed || mutation.changed;
          const rollbackErrors = rollbackJournal(mutation, { reason, resultIndex: target.resultIndex });
          if (rollbackErrors.length > 0) {
            throw attachRollbackErrors(
              new Error(`Live edit lease cancellation failed during ${String(reason || 'preview-edit')}.`),
              rollbackErrors
            );
          }
          return hadChanges;
        } finally {
          release();
        }
      }
    };
    activeDomEditLease = lease;
    return lease;
  };

  const settleDomEditBatch = ({ completedBatch, target, commitState, reason }) => {
    let flushStatus = { flushed: false, replacement: null, stale: false };
    let rollbackComplete = false;
    try {
      if (!completedBatch.stale && !isCurrentCommitTarget(target)) {
        completedBatch.stale = true;
        rejectStaleTarget(target, reason);
      }
      if (completedBatch.stale) {
        const rollbackErrors = rollbackJournals(
          completedBatch.mutations,
          completedBatch.canonicalMutation,
          { reason, resultIndex: target.resultIndex }
        );
        restoreCommitStateSafely(target, commitState, null, reason, rollbackErrors);
        rollbackComplete = true;
        if (rollbackErrors.length > 0) {
          throw attachRollbackErrors(staleRollbackError(reason), rollbackErrors);
        }
        return false;
      }
      if (completedBatch.changed) {
        flushStatus = flushCommitTarget(target);
        if (flushStatus.stale) {
          const rollbackErrors = rollbackJournals(
            completedBatch.mutations,
            completedBatch.canonicalMutation,
            { reason, resultIndex: target.resultIndex }
          );
          restoreCommitStateSafely(target, commitState, null, reason, rollbackErrors);
          rollbackComplete = true;
          if (rollbackErrors.length > 0) {
            throw attachRollbackErrors(staleRollbackError(reason), rollbackErrors);
          }
          return false;
        }
      }
      completedBatch.mutations.forEach((mutation) => mutation.commit());
      completedBatch.canonicalMutation.commit();
      return true;
    } catch (error) {
      if (rollbackComplete) throw error;
      const rollbackErrors = rollbackJournals(
        completedBatch.mutations,
        completedBatch.canonicalMutation,
        { reason, resultIndex: target.resultIndex }
      );
      restoreCommitStateSafely(
        target,
        commitState,
        failedResultReplacement(error, flushStatus.replacement),
        reason,
        rollbackErrors
      );
      throw attachRollbackErrors(error, rollbackErrors);
    }
  };

  const startDomEditBatch = ({ reason, canonicalState }) => {
    const target = resolveActiveCommitTarget();
    if (!target) return null;
    const batch = {
      target,
      changed: false,
      stale: false,
      reason,
      mutations: [],
      canonicalMutation: captureCanonicalState(canonicalState)
    };
    return { batch, commitState: captureCommitState(target), target };
  };

  const rollbackActionFailure = ({ batch, commitState, target, reason }, error) => {
    const rollbackErrors = rollbackJournals(
      batch.mutations,
      batch.canonicalMutation,
      { reason, resultIndex: target.resultIndex }
    );
    restoreCommitStateSafely(target, commitState, null, reason, rollbackErrors);
    throw attachRollbackErrors(error, rollbackErrors);
  };

  /**
   * Batch only the synchronous prefix of a user action. If the action returns a
   * Promise, direct edits are committed before the first asynchronous settlement;
   * later continuation edits resolve a new current target. Mounted ownership is
   * never retained across await. The optional canonicalState descriptors capture
   * only small editor-owned arrays, maps, sets, refs, and plain objects.
   *
   * @template T
   * @param {Object} options
   * @param {string} [options.reason]
   * @param {() => (T|Promise<T>)} options.action
   * @param {Array<object>} [options.canonicalState]
   * @returns {Promise<T>}
   */
  const runDomEdit = async ({
    reason = 'preview-edit',
    action,
    canonicalState = []
  } = {}) => {
    if (typeof action !== 'function') throw new TypeError('runDomEdit requires an action callback.');
    if (activeDomEditBatch) {
      captureCanonicalStateTargets(activeDomEditBatch.canonicalMutation, canonicalState);
      return action();
    }

    const context = startDomEditBatch({ reason, canonicalState });
    if (!context) {
      if (getMountedSvg()) throw staleRollbackError(reason);
      return action();
    }
    activeDomEditBatch = context.batch;
    let result;
    try {
      result = action();
    } catch (error) {
      activeDomEditBatch = null;
      return rollbackActionFailure(context, error);
    }
    const completedBatch = activeDomEditBatch;
    activeDomEditBatch = null;
    settleDomEditBatch({ ...context, completedBatch, reason });
    return await result;
  };

  /**
   * Keep one explicit mounted-document lease across an asynchronous user action.
   * Callers must route every DOM mutation through the supplied lease; canonical
   * editor state is captured at action entry and the Result is serialized once,
   * after the final continuation has settled.
   */
  const runDomEditTransaction = async ({
    reason = 'preview-edit',
    action,
    canonicalState = [],
    invalidateIndexes = []
  } = {}) => {
    if (typeof action !== 'function') {
      throw new TypeError('runDomEditTransaction requires an action callback.');
    }
    const lease = beginDomEditLease({ reason, invalidateIndexes });
    if (!lease) {
      if (getMountedSvg()) throw staleRollbackError(reason);
      return action({ lease: null, target: null });
    }
    lease.mutate(({ mutation }) => {
      captureCanonicalStateTargets(mutation, canonicalState);
      return false;
    });
    try {
      const result = await action({ lease, target: lease.target });
      lease.commit();
      return result;
    } catch (error) {
      const rollbackErrors = [];
      try {
        lease.cancel();
      } catch (rollbackError) {
        if (Array.isArray(rollbackError?.rollbackErrors)) {
          rollbackErrors.push(...rollbackError.rollbackErrors);
        } else {
          rollbackErrors.push(
            rollbackError instanceof Error ? rollbackError : new Error(String(rollbackError))
          );
        }
      }
      throw attachRollbackErrors(error, rollbackErrors);
    }
  };

  const runDomEditSync = ({
    reason = 'preview-edit',
    action,
    canonicalState = []
  } = {}) => {
    if (typeof action !== 'function') throw new TypeError('runDomEditSync requires an action callback.');
    if (activeDomEditBatch) {
      captureCanonicalStateTargets(activeDomEditBatch.canonicalMutation, canonicalState);
      return action();
    }

    const context = startDomEditBatch({ reason, canonicalState });
    if (!context) {
      if (getMountedSvg()) throw staleRollbackError(reason);
      return action();
    }
    activeDomEditBatch = context.batch;
    let result;
    try {
      result = action();
      if (result && typeof result.then === 'function') {
        throw new TypeError('runDomEditSync action callback must be synchronous.');
      }
    } catch (error) {
      activeDomEditBatch = null;
      return rollbackActionFailure(context, error);
    }
    const completedBatch = activeDomEditBatch;
    activeDomEditBatch = null;
    settleDomEditBatch({ ...context, completedBatch, reason });
    return result;
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
    const commitState = captureCommitState(target);
    if (force && !target.runtime.dirty) target.runtime.dirty = true;
    try {
      return flushCommitTarget(target, { markIncremental }).flushed;
    } catch (error) {
      const rollbackErrors = restoreCommitState(
        target,
        commitState,
        failedResultReplacement(error),
        'flush-active-result'
      );
      throw attachRollbackErrors(error, rollbackErrors);
    }
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

  const getFeatureElementsForRuntime = (runtime, featureId) => {
    const normalizedId = normalizeFeatureIdentity(featureId);
    if (!runtime?.svg || !normalizedId) return [];

    const featureIndex = runtime.indexes.features || buildFeatureIndex(runtime);
    const indexed = featureIndex.get(normalizedId);
    if (indexed?.length) return indexed;

    const byId = runtime.svg.getElementById?.(normalizedId);
    return byId ? [byId] : [];
  };

  const getFeatureElements = (featureId) => (
    getFeatureElementsForRuntime(activeRuntime || ensureRuntimeForCurrentSvg(), featureId)
  );

  const applyFeatureFillChanges = (changes, { reason = 'feature-fill', lease = null } = {}) => {
    const normalized = normalizeChanges(changes);
    if (normalized.length === 0) return false;
    const mutate = ({ mutation }) => {
        let updated = 0;
        normalized.forEach((change) => {
          const color = String(change?.color || '').trim();
          if (!color) return;
          const elements = lease
            ? getFeatureElementsForRuntime(lease.target.runtime, change.featureId)
            : getFeatureElements(change.featureId);
          filterFeatureFillTargets(elements).forEach((element) => {
            if (element.getAttribute?.('fill') === color) return;
            mutation.setAttribute(element, 'fill', color);
            updated += 1;
          });
        });
        return updated;
      };
    if (lease) return lease.mutate(mutate);
    return commitDomEdit({ reason, invalidateIndexes: ['features'], mutate }).changed;
  };

  const applyFeatureVisibilityChanges = (changes, { reason = 'feature-visibility', lease = null } = {}) => {
    const normalized = normalizeChanges(changes);
    if (normalized.length === 0) return false;
    const mutate = ({ mutation }) => {
        let updated = 0;
        normalized.forEach((change) => {
          const elements = lease
            ? getFeatureElementsForRuntime(lease.target.runtime, change.featureId)
            : getFeatureElements(change.featureId);
          elements.forEach((element) => {
            if (applyFeatureVisibility(element, change?.mode, {
              markPreview: true,
              mutation,
              reason
            })) updated += 1;
          });
        });
        return updated;
      };
    if (lease) return lease.mutate(mutate);
    return commitDomEdit({ reason, invalidateIndexes: ['features'], mutate }).changed;
  };

  const applyFeatureStrokeChanges = (changes, { reason = 'feature-stroke', lease = null } = {}) => {
    const normalized = normalizeChanges(changes);
    if (normalized.length === 0) return false;
    const mutate = ({ mutation }) => {
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
          const elements = lease
            ? getFeatureElementsForRuntime(lease.target.runtime, change.featureId)
            : getFeatureElements(change.featureId);
          elements.forEach((element) => {
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
      };
    if (lease) return lease.mutate(mutate);
    return commitDomEdit({ reason, invalidateIndexes: ['features'], mutate }).changed;
  };

  return {
    applyFeatureFillChanges,
    applyFeatureStrokeChanges,
    applyFeatureVisibilityChanges,
    beginDomEditLease,
    captureDomEditToken,
    clearActiveRuntime,
    commitDomEdit,
    flushActiveResult,
    getActiveRuntime,
    getFeatureElements,
    invalidatePreviewIndexes,
    isDomEditTokenCurrent,
    mountResultSvg,
    runDomEdit,
    runDomEditTransaction,
    runDomEditSync,
    selectResult
  };
};
