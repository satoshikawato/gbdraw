import {
  filterFeatureFillTargets,
  getFeatureElementIndex,
  getFeatureIdentity,
  normalizeFeatureIdentity
} from './feature-dom.js';
import {
  getCommittedSvgResultMetadata,
  getCommittedSvgResultRuntimeIdentity,
  markCommittedSvgResultMounted,
  markCommittedSvgResultUnmounted
} from '../services/svg-result-ingestion.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from '../services/runtime-test-hooks.js';

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

const REQUIRED_BINDING_FLAGS = Object.freeze([
  'rootAdopted',
  'legendReady',
  'compositionReady',
  'dragReady',
  'interactionsReady',
  'labelEditorReady',
  'strokeCanvasReady',
  'selectionReady'
]);

const readinessError = (message, code = 'PREVIEW_READINESS_REJECTED') => {
  const error = new Error(message);
  error.code = code;
  return error;
};

const makeRuntime = ({
  resultIndex,
  svg,
  result,
  resultIdentity,
  artifactIdentity,
  generationToken,
  rootGeneration,
  bindSequence
}) => ({
  resultIndex,
  svg,
  result,
  resultIdentity,
  artifactIdentity,
  generationToken,
  rootGeneration,
  bindSequence,
  readyReceipt: null,
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
  let activeExpectation = null;
  let pendingBind = null;
  let bindingSteps = Object.freeze({});
  let nextExpectationIdentity = 1;
  let nextFallbackResultIdentity = 1;
  let nextRootGeneration = 1;
  let nextBindSequence = 1;
  let nextRestoreToken = 1;
  let lastObservedMount = null;
  const fallbackResultIdentities = new WeakMap();
  const invalidatedReceipts = new WeakSet();

  const getMountedSvg = () => state.svgContainer?.value?.querySelector?.('svg') || null;

  const resultRuntimeIdentity = (result) => {
    const committedIdentity = getCommittedSvgResultRuntimeIdentity(result);
    if (committedIdentity !== null) return `result:${committedIdentity}`;
    if (!result || typeof result !== 'object') return '';
    if (!fallbackResultIdentities.has(result)) {
      fallbackResultIdentities.set(result, `runtime-result:${nextFallbackResultIdentity++}`);
    }
    return fallbackResultIdentities.get(result);
  };

  const normalizeArtifactIdentity = (value, resultIdentity = '') => {
    const fingerprint = String(value?.fingerprint || '').trim().toLowerCase();
    if (fingerprint) return fingerprint;
    const normalized = String(value || '').trim();
    return normalized || `preview-artifact:${resultIdentity || 'unknown'}`;
  };

  const rejectExpectation = (expectation, reason) => {
    if (!expectation || expectation.settled) return false;
    expectation.settled = true;
    expectation.reject(
      reason instanceof Error
        ? reason
        : readinessError(String(reason || 'Preview readiness was invalidated.'))
    );
    if (activeExpectation === expectation) activeExpectation = null;
    return true;
  };

  const invalidateReadinessExpectation = (target = null, reason = null) => {
    const expectation = activeExpectation;
    if (!expectation) return false;
    if (
      target
      && target !== expectation
      && target !== expectation.generationToken
      && target?.generationToken !== expectation.generationToken
    ) return false;
    const rejected = rejectExpectation(
      expectation,
      reason || readinessError('Preview readiness was superseded.', 'PREVIEW_READINESS_SUPERSEDED')
    );
    if (rejected) {
      recordStructuralMetric('previewReadyReceiptRejectedCount', 1, {
        phase: expectation.phase,
        rootGeneration: expectation.rootGeneration || 0
      });
      recordSessionLifecycleEvent('preview.ready-receipt-rejected', {
        phase: expectation.phase,
        resultIndex: expectation.resultIndex
      });
    }
    return rejected;
  };

  const registerReadinessExpectation = ({
    result,
    resultIndex = 0,
    artifactIdentity = '',
    generationToken = '',
    catalogState = null,
    phase = 'preview',
    bindingOptions = {},
    isCurrent = () => true
  } = {}) => {
    if (!result || typeof result !== 'object') {
      throw new Error('Preview readiness requires a selected Result.');
    }
    if (typeof isCurrent !== 'function') {
      throw new Error('Preview readiness requires a current-operation predicate.');
    }
    invalidateReadinessExpectation(
      null,
      readinessError('Preview readiness was replaced.', 'PREVIEW_READINESS_REPLACED')
    );
    const normalizedIndex = Number(resultIndex);
    const resultIdentity = resultRuntimeIdentity(result);
    const normalizedGenerationToken = String(generationToken || '').trim()
      || `preview:${nextExpectationIdentity}`;
    let resolvePromise;
    let rejectPromise;
    const promise = new Promise((resolve, reject) => {
      resolvePromise = resolve;
      rejectPromise = reject;
    });
    // Readiness may be owned by a UI transition whose caller does not await it.
    // Keep the rejection observed without changing what awaiting callers receive.
    void promise.catch(() => {});
    const expectation = {
      identity: nextExpectationIdentity++,
      result,
      resultIndex: Number.isInteger(normalizedIndex) ? normalizedIndex : 0,
      resultIdentity,
      artifactIdentity: normalizeArtifactIdentity(artifactIdentity, resultIdentity),
      generationToken: normalizedGenerationToken,
      catalogState,
      phase: String(phase || 'preview'),
      bindingOptions: Object.freeze({ ...(bindingOptions || {}) }),
      bindSequence: nextBindSequence++,
      rootGeneration: 0,
      isCurrent,
      settled: false,
      promise,
      resolve: resolvePromise,
      reject: rejectPromise
    };
    activeExpectation = expectation;
    recordStructuralMetric('generatedArtifactReadinessExpectationCount', 1, {
      phase: expectation.phase
    });
    recordSessionLifecycleEvent('preview.readiness-expectation-registered', {
      phase: expectation.phase,
      resultIndex: expectation.resultIndex,
      bindSequence: expectation.bindSequence
    });
    return expectation;
  };

  const releaseActiveResult = () => {
    if (!activeRuntime) return;
    bindingSteps.disposeMountedResult?.(activeRuntime);
    markCommittedSvgResultUnmounted(activeRuntime.result);
  };

  const mountResultSvg = (resultIndex = state.selectedResultIndex?.value || 0, svg = getMountedSvg()) => {
    if (!svg) {
      releaseActiveResult();
      activeRuntime = null;
      return null;
    }
    const normalizedIndex = Number(resultIndex) || 0;
    const result = state.results.value[normalizedIndex] || null;
    const resultIdentity = resultRuntimeIdentity(result);
    if (
      activeRuntime?.svg === svg
      && activeRuntime.resultIndex === normalizedIndex
    ) {
      return activeRuntime;
    }
    releaseActiveResult();
    activeRuntime = makeRuntime({
      resultIndex: normalizedIndex,
      svg,
      result,
      resultIdentity,
      artifactIdentity: activeExpectation?.artifactIdentity || `mounted:${resultIdentity}`,
      generationToken: activeExpectation?.generationToken || '',
      rootGeneration: nextRootGeneration++,
      bindSequence: activeExpectation?.bindSequence || nextBindSequence++
    });
    markCommittedSvgResultMounted(result);
    return activeRuntime;
  };

  const clearActiveRuntime = () => {
    invalidateReadinessExpectation(
      null,
      readinessError('The mounted preview was cleared.', 'PREVIEW_RUNTIME_CLEARED')
    );
    releaseActiveResult();
    activeRuntime = null;
    pendingBind = null;
    lastObservedMount = null;
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

  const configureMountedResultBinder = (steps = {}) => {
    if (!steps || typeof steps !== 'object' || Array.isArray(steps)) {
      throw new Error('PreviewRuntime binder steps must be an object.');
    }
    bindingSteps = Object.freeze({ ...steps });
  };

  const createMountedResultContext = ({
    root,
    result = state.results.value[state.selectedResultIndex?.value || 0],
    resultIndex = state.selectedResultIndex?.value || 0,
    catalogState = state.featureCatalog?.value || null,
    bindingOptions = {}
  } = {}) => {
    if (!root) throw new Error('Mounted preview observation requires an SVG root.');
    const normalizedIndex = Number(resultIndex) || 0;
    const resultIdentity = resultRuntimeIdentity(result);
    if (
      activeRuntime?.svg === root
      && activeRuntime.resultIdentity === resultIdentity
      && activeRuntime.resultIndex === normalizedIndex
      && activeRuntime.readyReceipt
    ) {
      return Object.freeze({
        root,
        result,
        sourceClass: getCommittedSvgResultMetadata(result)?.sourceClass || '',
        resultIndex: normalizedIndex,
        resultIdentity,
        artifactIdentity: activeRuntime.artifactIdentity,
        generationToken: activeRuntime.generationToken,
        catalogState,
        bindSequence: activeRuntime.bindSequence,
        rootGeneration: activeRuntime.rootGeneration,
        expectationIdentity: 0,
        phase: activeRuntime.readyReceipt.phase,
        bindingOptions: Object.freeze({ ...(bindingOptions || {}) })
      });
    }
    let expectation = activeExpectation;
    if (
      expectation
      && (
        expectation.resultIndex !== normalizedIndex
        || expectation.resultIdentity !== resultIdentity
      )
    ) {
      invalidateReadinessExpectation(
        expectation,
        readinessError('The mounted Result does not match the readiness expectation.')
      );
      expectation = null;
    }
    if (!expectation) {
      expectation = registerReadinessExpectation({
        result,
        resultIndex: normalizedIndex,
        artifactIdentity: `passive:${resultIdentity}`,
        generationToken: `passive:${nextExpectationIdentity}`,
        catalogState,
        phase: 'passive-preview',
        bindingOptions,
        isCurrent: () => (
          resultRuntimeIdentity(state.results.value[normalizedIndex]) === resultIdentity
          && Number(state.selectedResultIndex?.value || 0) === normalizedIndex
        )
      });
    }
    const rootGeneration = activeRuntime?.svg === root
      ? activeRuntime.rootGeneration
      : nextRootGeneration++;
    expectation.rootGeneration = rootGeneration;
    const context = Object.freeze({
      root,
      result,
      sourceClass: getCommittedSvgResultMetadata(result)?.sourceClass || '',
      resultIndex: normalizedIndex,
      resultIdentity,
      artifactIdentity: expectation.artifactIdentity,
      generationToken: expectation.generationToken,
      catalogState: expectation.catalogState || catalogState,
      bindSequence: expectation.bindSequence,
      rootGeneration,
      expectationIdentity: expectation.identity,
      phase: expectation.phase,
      bindingOptions: Object.freeze({
        ...expectation.bindingOptions,
        ...(bindingOptions || {})
      })
    });
    if (
      !lastObservedMount
      || lastObservedMount.root !== root
      || lastObservedMount.resultIdentity !== resultIdentity
      || lastObservedMount.bindSequence !== context.bindSequence
    ) {
      lastObservedMount = {
        root,
        resultIdentity,
        bindSequence: context.bindSequence
      };
      recordStructuralMetric('previewMaterializationObservedCount', 1, {
        phase: context.phase,
        rootGeneration
      });
      recordSessionLifecycleEvent('preview.mount-observed', {
        phase: context.phase,
        resultIndex: normalizedIndex,
        rootGeneration,
        bindSequence: context.bindSequence
      });
    }
    return context;
  };

  const rejectReadyReceipt = (receipt, message) => {
    recordStructuralMetric('previewReadyReceiptRejectedCount', 1, {
      phase: receipt?.phase || activeExpectation?.phase || 'preview',
      rootGeneration: Number(receipt?.rootGeneration) || 0
    });
    recordSessionLifecycleEvent('preview.ready-receipt-rejected', {
      phase: receipt?.phase || activeExpectation?.phase || 'preview',
      resultIndex: Number(receipt?.resultIndex) || 0
    });
    return { accepted: false, error: readinessError(message) };
  };

  const acceptReadyReceipt = (receipt) => {
    const expectation = activeExpectation;
    if (!expectation || !receipt || invalidatedReceipts.has(receipt)) {
      return rejectReadyReceipt(receipt, 'No matching preview readiness expectation is active.');
    }
    const flags = receipt.requiredBindingFlags || {};
    const complete = REQUIRED_BINDING_FLAGS.every((flag) => flags[flag] === true);
    const matches = (
      receipt.artifactIdentity === expectation.artifactIdentity
      && receipt.resultIdentity === expectation.resultIdentity
      && receipt.resultIndex === expectation.resultIndex
      && receipt.generationToken === expectation.generationToken
      && receipt.bindSequence === expectation.bindSequence
      && receipt.rootGeneration === expectation.rootGeneration
    );
    let current = false;
    try {
      current = Boolean(expectation.isCurrent());
    } catch (_error) {
      current = false;
    }
    if (!matches || !complete || !current || expectation.settled) {
      const rejection = rejectReadyReceipt(
        receipt,
        !complete
          ? 'Preview binding did not complete every required substep.'
          : 'The preview readiness receipt is stale or mismatched.'
      );
      rejectExpectation(expectation, rejection.error);
      return rejection;
    }
    expectation.settled = true;
    expectation.resolve(receipt);
    activeExpectation = null;
    recordStructuralMetric('previewReadyReceiptAcceptedCount', 1, {
      phase: expectation.phase,
      rootGeneration: receipt.rootGeneration
    });
    recordSessionLifecycleEvent('preview.ready-receipt-accepted', {
      phase: expectation.phase,
      resultIndex: receipt.resultIndex,
      rootGeneration: receipt.rootGeneration,
      bindSequence: receipt.bindSequence
    });
    return { accepted: true, receipt };
  };

  const assertBindContextCurrent = (context) => {
    const expectation = activeExpectation;
    if (
      !expectation
      || expectation.identity !== context.expectationIdentity
      || expectation.generationToken !== context.generationToken
      || expectation.bindSequence !== context.bindSequence
    ) {
      throw readinessError('Preview binding was superseded.', 'PREVIEW_BIND_SUPERSEDED');
    }
    if (expectation.settled || !expectation.isCurrent()) {
      throw readinessError('Preview binding is stale or canceled.', 'PREVIEW_BIND_STALE');
    }
    if (
      resultRuntimeIdentity(state.results.value[context.resultIndex]) !== context.resultIdentity
      || Number(state.selectedResultIndex?.value || 0) !== context.resultIndex
      || getMountedSvg() !== context.root
      || context.root?.isConnected === false
    ) {
      throw readinessError('The mounted preview root is not current.', 'PREVIEW_ROOT_MISMATCH');
    }
  };

  const bindMountedResult = (context = {}) => {
    if (
      pendingBind
      && pendingBind.root === context.root
      && pendingBind.bindSequence === context.bindSequence
    ) {
      recordStructuralMetric('previewDuplicateBindRejectedCount', 1, {
        phase: context.phase,
        rootGeneration: context.rootGeneration
      });
      return pendingBind.promise;
    }
    if (
      activeRuntime?.svg === context.root
      && activeRuntime.bindSequence === context.bindSequence
      && activeRuntime.readyReceipt
    ) {
      recordStructuralMetric('previewDuplicateBindRejectedCount', 1, {
        phase: context.phase,
        rootGeneration: context.rootGeneration
      });
      return Promise.resolve(activeRuntime.readyReceipt);
    }

    const promise = (async () => {
      assertBindContextCurrent(context);
      const previousRuntime = activeRuntime;
      if (
        previousRuntime
        && (
          previousRuntime.svg !== context.root
          || previousRuntime.resultIdentity !== context.resultIdentity
        )
      ) {
        releaseActiveResult();
        activeRuntime = null;
      }
      if (!activeRuntime) {
        activeRuntime = makeRuntime({
          resultIndex: context.resultIndex,
          svg: context.root,
          result: context.result,
          resultIdentity: context.resultIdentity,
          artifactIdentity: context.artifactIdentity,
          generationToken: context.generationToken,
          rootGeneration: context.rootGeneration,
          bindSequence: context.bindSequence
        });
        markCommittedSvgResultMounted(context.result);
        recordStructuralMetric('previewMountAdoptionCount', 1, {
          phase: context.phase,
          rootGeneration: context.rootGeneration
        });
      }
      recordStructuralMetric('previewBinderInvocationCount', 1, {
        phase: context.phase,
        rootGeneration: context.rootGeneration
      });
      recordStructuralMetric('previewPollingIterationCount', 0, {
        phase: context.phase,
        rootGeneration: context.rootGeneration
      });
      recordSessionLifecycleEvent('preview.bind-started', {
        phase: context.phase,
        resultIndex: context.resultIndex,
        rootGeneration: context.rootGeneration,
        bindSequence: context.bindSequence
      });

      const flags = {};
      const completeStep = async (flag, stepName) => {
        assertBindContextCurrent(context);
        if (typeof bindingSteps[stepName] === 'function') {
          await bindingSteps[stepName](context);
        }
        flags[flag] = true;
      };
      flags.rootAdopted = true;
      await completeStep('legendReady', 'adoptLegend');
      await completeStep('compositionReady', 'bindComposition');
      await completeStep('dragReady', 'setupDragAffordances');
      await completeStep('interactionsReady', 'installDelegatedInteractions');
      await completeStep('labelEditorReady', 'synchronizeLabelEditor');
      await completeStep('strokeCanvasReady', 'initializeStrokeAndCanvas');
      await completeStep('selectionReady', 'reconcileSelection');
      assertBindContextCurrent(context);

      const receipt = Object.freeze({
        artifactIdentity: context.artifactIdentity,
        resultIdentity: context.resultIdentity,
        resultIndex: context.resultIndex,
        generationToken: context.generationToken,
        rootIdentity: context.rootGeneration,
        rootGeneration: context.rootGeneration,
        bindSequence: context.bindSequence,
        requiredBindingFlags: Object.freeze({ ...flags }),
        readyTimestamp: globalThis.performance?.now?.() ?? Date.now(),
        phase: context.phase
      });
      recordSessionLifecycleEvent('preview.bind-completed', {
        phase: context.phase,
        resultIndex: context.resultIndex,
        rootGeneration: context.rootGeneration,
        bindSequence: context.bindSequence
      });
      recordStructuralMetric('previewReadyReceiptEmittedCount', 1, {
        phase: context.phase,
        rootGeneration: context.rootGeneration
      });
      recordSessionLifecycleEvent('preview.ready-receipt-emitted', {
        phase: context.phase,
        resultIndex: context.resultIndex,
        rootGeneration: context.rootGeneration,
        bindSequence: context.bindSequence
      });
      const acceptance = acceptReadyReceipt(receipt);
      if (!acceptance.accepted) throw acceptance.error;
      activeRuntime.readyReceipt = receipt;
      if (typeof bindingSteps.afterReady === 'function') {
        queueMicrotask(() => {
          Promise.resolve(bindingSteps.afterReady(context)).catch((error) => {
            console.error('Post-ready preview work failed.', error);
          });
        });
      }
      return receipt;
    })().catch((error) => {
      const expectation = activeExpectation;
      if (
        expectation
        && expectation.identity === context.expectationIdentity
        && rejectExpectation(expectation, error)
      ) {
        recordStructuralMetric('previewReadyReceiptRejectedCount', 1, {
          phase: context.phase,
          rootGeneration: context.rootGeneration
        });
        recordSessionLifecycleEvent('preview.ready-receipt-rejected', {
          phase: context.phase,
          resultIndex: context.resultIndex
        });
      }
      throw error;
    });
    pendingBind = {
      root: context.root,
      bindSequence: context.bindSequence,
      promise
    };
    void promise.finally(() => {
      if (pendingBind?.promise === promise) pendingBind = null;
    }).catch(() => {});
    return promise;
  };

  const invalidateReadyReceipt = (receipt, reason = 'Preview readiness entered rollback.') => {
    if (!receipt || typeof receipt !== 'object') return false;
    invalidatedReceipts.add(receipt);
    if (activeRuntime?.readyReceipt === receipt) activeRuntime.readyReceipt = null;
    recordStructuralMetric('previewReadyReceiptRejectedCount', 1, {
      phase: receipt.phase || 'rollback-restoration',
      rootGeneration: receipt.rootGeneration || 0
    });
    recordSessionLifecycleEvent('preview.ready-receipt-rejected', {
      phase: receipt.phase || 'rollback-restoration',
      resultIndex: receipt.resultIndex,
      reason: String(reason || 'rollback')
    });
    return true;
  };

  const restorePreviousSelectedResult = async ({
    handle,
    restore,
    phase = 'rollback-restoration'
  } = {}) => {
    if (typeof restore !== 'function') {
      throw new Error('Preview restoration requires the generated-artifact restore owner.');
    }
    invalidateReadinessExpectation(
      null,
      readinessError('Candidate preview readiness entered rollback.', 'PREVIEW_ROLLBACK')
    );
    const beforeResults = Array.isArray(handle?.ownerSet?.results)
      ? handle.ownerSet.results
      : [];
    const requestedIndex = Number(handle?.mutableIntent?.ui?.selectedResultIndex);
    const resultIndex = Number.isInteger(requestedIndex)
      ? Math.max(0, Math.min(requestedIndex, Math.max(0, beforeResults.length - 1)))
      : 0;
    const result = beforeResults[resultIndex] || null;
    const expectation = result
      ? registerReadinessExpectation({
          result,
          resultIndex,
          artifactIdentity: handle?.identity?.fingerprint || `restored:${resultRuntimeIdentity(result)}`,
          generationToken: `restore:${nextRestoreToken++}`,
          catalogState: handle?.ownerSet?.featureCatalog || null,
          phase,
          bindingOptions: { trustedRestore: true, isIncrementalEdit: true },
          isCurrent: () => (
            resultRuntimeIdentity(state.results.value[resultIndex]) === resultRuntimeIdentity(result)
            && Number(state.selectedResultIndex?.value || 0) === resultIndex
          )
        })
      : null;
    recordSessionLifecycleEvent('preview.restore-bind-started', {
      phase,
      resultIndex
    });
    await restore();
    if (!expectation) {
      clearActiveRuntime();
      recordSessionLifecycleEvent('preview.restore-bind-completed', { phase, resultIndex });
      return null;
    }
    const receipt = await expectation.promise;
    recordSessionLifecycleEvent('preview.restore-bind-completed', {
      phase,
      resultIndex,
      rootGeneration: receipt.rootGeneration
    });
    recordSessionLifecycleEvent('preview.restore-ready-receipt-accepted', {
      phase,
      resultIndex,
      rootGeneration: receipt.rootGeneration
    });
    return receipt;
  };

  const isActiveResultReady = () => {
    if (!activeRuntime?.readyReceipt) return false;
    const resultIndex = Number(state.selectedResultIndex?.value || 0);
    return (
      activeRuntime.resultIndex === resultIndex
      && activeRuntime.resultIdentity === resultRuntimeIdentity(state.results.value[resultIndex])
      && activeRuntime.svg === getMountedSvg()
      && activeRuntime.svg?.isConnected !== false
    );
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
    const nextResults = [...state.results.value];
    nextResults[resultIndex] = {
      ...state.results.value[resultIndex],
      content
    };
    state.results.value = nextResults;
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
    const nextResult = state.results.value[nextIndex];
    const expectation = nextResult
      ? registerReadinessExpectation({
          result: nextResult,
          resultIndex: nextIndex,
          artifactIdentity: activeRuntime?.artifactIdentity || `selection:${resultRuntimeIdentity(nextResult)}`,
          generationToken: `result-selection:${nextBindSequence}`,
          catalogState: state.featureCatalog?.value || null,
          phase: 'result-selection',
          isCurrent: () => (
            state.results.value[nextIndex] === nextResult
            && Number(state.selectedResultIndex.value) === nextIndex
          )
        })
      : null;
    if (expectation) void expectation.promise.catch(() => {});
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

    if (updated === 0) return false;
    markActiveResultDirty(reason);
    flushActiveResult();
    return true;
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
    acceptReadyReceipt,
    applyFeatureFillChanges,
    applyFeatureStrokeChanges,
    applyFeatureVisibilityChanges,
    applyLegendChanges,
    bindMountedResult,
    clearActiveRuntime,
    configureMountedResultBinder,
    createMountedResultContext,
    flushActiveResult,
    getActiveRuntime,
    getFeatureElements,
    invalidateReadinessExpectation,
    invalidateReadyReceipt,
    invalidatePreviewIndexes,
    isActiveResultReady,
    markActiveResultDirty,
    mountResultSvg,
    registerReadinessExpectation,
    restorePreviousSelectedResult,
    selectResult
  };
};
