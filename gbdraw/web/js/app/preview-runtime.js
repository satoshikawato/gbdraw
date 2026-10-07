// @ts-check
import { normalizeUserFacingError } from '../services/error-normalization.js';
import {
  getFeatureElementIndex,
  normalizeFeatureIdentity
} from './feature-dom.js';
import {
  applyEditorOperationsToMountedSvg,
  getCommittedSvgResultMetadata,
  getCommittedSvgResultRuntimeIdentity,
  markCommittedSvgResultMounted,
  markCommittedSvgResultUnmounted
} from '../services/svg-result-ingestion.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from '../services/runtime-test-hooks.js';

/**
 * What the Label editor and the mounted-preview binder steps read from
 * `bindingOptions`. A caller names only the options its phase needs.
 * @typedef {object} PreviewBindingOptions
 * @property {boolean} [isIncrementalEdit] The Result is edited in place, not newly generated.
 * @property {boolean} [skipLegendExtraction]
 * @property {boolean} [trustedRestore] A restored artifact: its editor state is already current.
 * @property {boolean} [replaceGeneratedLegend]
 * @property {readonly string[]} [requiredLabelFeatureIds]
 * @property {readonly string[]} [optionalLabelFeatureIds]
 * @property {{ featureIds: readonly string[], report: (error: unknown) => void }} [reportedLabelBinding]
 */

/**
 * The refs the runtime reads and writes (state.js owns them).
 * @typedef {object} PreviewRuntimeState
 * @property {{ value: any }} [svgContainer] The element that holds the mounted SVG.
 * @property {{ value: Record<string, any>[] }} results
 * @property {{ value: number }} selectedResultIndex
 * @property {{ value: any }} [featureCatalog]
 * @property {{ value: boolean }} [skipCaptureBaseConfig]
 * @property {() => any} [sessionOperationAvailability]
 */

/**
 * @typedef {object} PreviewRuntimeOptions
 * @property {PreviewRuntimeState} state
 * @property {(svg: Element) => string} serializeSvg The clean serializer of services/svg-serialization.js.
 */

/**
 * The preview readiness a caller registers before it changes the mounted
 * Result, and the receipt the binder must return for it.
 * @typedef {object} ReadinessExpectationOptions
 * @property {Record<string, any>} result
 * @property {number} [resultIndex]
 * @property {any} [artifactIdentity] A fingerprint object or string.
 * @property {string} [generationToken]
 * @property {any} [catalogState] The feature catalog admission of the Result.
 * @property {string} [phase]
 * @property {PreviewBindingOptions} [bindingOptions]
 * @property {() => boolean} [isCurrent] Whether the operation that registered it is still current.
 */

/**
 * @typedef {object} ReadinessExpectation
 * @property {number} identity
 * @property {Record<string, any>} result
 * @property {number} resultIndex
 * @property {string} resultIdentity
 * @property {string} artifactIdentity
 * @property {string} generationToken
 * @property {any} catalogState
 * @property {string} phase
 * @property {Readonly<PreviewBindingOptions>} bindingOptions
 * @property {number} bindSequence
 * @property {number} rootGeneration
 * @property {() => boolean} isCurrent
 * @property {boolean} settled
 * @property {Promise<ReadyReceipt>} promise
 * @property {(receipt: ReadyReceipt) => void} resolve
 * @property {(reason: Error) => void} reject
 */

/**
 * @typedef {object} ReadyReceipt
 * @property {string} artifactIdentity
 * @property {string} resultIdentity
 * @property {number} resultIndex
 * @property {string} generationToken
 * @property {number} rootIdentity
 * @property {number} rootGeneration
 * @property {number} bindSequence
 * @property {Readonly<Record<string, boolean>>} requiredBindingFlags
 * @property {number} readyTimestamp
 * @property {string} phase
 */

/**
 * One mounted Result and the lazily built indexes of its SVG.
 * @typedef {object} PreviewResultRuntime
 * @property {number} resultIndex
 * @property {SVGSVGElement} svg
 * @property {Record<string, any>} result
 * @property {string} resultIdentity
 * @property {string} artifactIdentity
 * @property {string} generationToken
 * @property {number} rootGeneration
 * @property {number} bindSequence
 * @property {ReadyReceipt | null} readyReceipt
 * @property {boolean} dirty
 * @property {Set<string>} dirtyReasons
 * @property {{ features: Map<string, Element[]> | null, legend: any, pairwiseMatches: any, orthogroupComparisons: any }} indexes
 * @property {string} [lastInvalidationReason]
 */

/**
 * The mount a caller observes. A missing `result` and `resultIndex` mean the
 * selected Result.
 * @typedef {object} MountedResultContextOptions
 * @property {SVGSVGElement} root
 * @property {Record<string, any>} [result]
 * @property {number} [resultIndex]
 * @property {any} [catalogState]
 * @property {PreviewBindingOptions} [bindingOptions]
 */

/**
 * What a binder step receives for one observed mount.
 * @typedef {object} MountedResultContext
 * @property {SVGSVGElement} root
 * @property {Record<string, any>} result
 * @property {string} sourceClass
 * @property {number} resultIndex
 * @property {string} resultIdentity
 * @property {string} artifactIdentity
 * @property {string} generationToken
 * @property {any} catalogState
 * @property {number} bindSequence
 * @property {number} rootGeneration
 * @property {number} expectationIdentity `0` for a mount no expectation owns.
 * @property {string} phase
 * @property {Readonly<PreviewBindingOptions>} bindingOptions
 */

/**
 * The steps of the mounted-preview binder, in the order they run. The root
 * wires them with `configureMountedResultBinder`; each is optional.
 * @typedef {object} MountedResultBinderSteps
 * @property {(context: MountedResultContext) => unknown} [adoptLegend]
 * @property {(context: MountedResultContext) => unknown} [bindComposition]
 * @property {(context: MountedResultContext) => unknown} [setupDragAffordances]
 * @property {(context: MountedResultContext) => unknown} [installDelegatedInteractions]
 * @property {(context: MountedResultContext) => unknown} [synchronizeLabelEditor]
 * @property {(context: MountedResultContext) => unknown} [initializeStrokeAndCanvas]
 * @property {(context: MountedResultContext) => unknown} [reconcileSelection]
 * @property {(context: MountedResultContext) => unknown} [afterReady] Runs after the receipt is accepted.
 * @property {(runtime: PreviewResultRuntime) => unknown} [disposeMountedResult]
 */

/**
 * @typedef {object} PreviousResultRestoreOptions
 * @property {Record<string, any>} [handle] The generated-artifact handle that holds the owner set to restore.
 * @property {() => unknown} restore The artifact owner's restore.
 * @property {string} [phase]
 */

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
  const error = /** @type {Error & { code?: string }} */ (new Error(message));
  error.code = code;
  return error;
};

/**
 * @param {Pick<PreviewResultRuntime,
 *   'resultIndex' | 'svg' | 'result' | 'resultIdentity' | 'artifactIdentity'
 *   | 'generationToken' | 'rootGeneration' | 'bindSequence'>} identity
 * @returns {PreviewResultRuntime}
 */
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

/** @param {PreviewRuntimeOptions} options */
export const createPreviewRuntime = ({ state, serializeSvg }) => {
  if (!state) throw new Error('createPreviewRuntime requires state.');
  if (typeof serializeSvg !== 'function') throw new Error('createPreviewRuntime requires serializeSvg.');

  /** @type {PreviewResultRuntime | null} */
  let activeRuntime = null;
  /** @type {ReadinessExpectation | null} */
  let activeExpectation = null;
  /** @type {{ root: SVGSVGElement, bindSequence: number, promise: Promise<ReadyReceipt> } | null} */
  let pendingBind = null;
  /** @type {Readonly<MountedResultBinderSteps>} */
  let bindingSteps = Object.freeze({});
  let nextExpectationIdentity = 1;
  let nextFallbackResultIdentity = 1;
  let nextRootGeneration = 1;
  let nextBindSequence = 1;
  let nextRestoreToken = 1;
  /** @type {{ root: SVGSVGElement, resultIdentity: string, bindSequence: number } | null} */
  let lastObservedMount = null;
  const fallbackResultIdentities = new WeakMap();
  const invalidatedReceipts = new WeakSet();

  /** @returns {SVGSVGElement | null} */
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

  /**
   * @param {ReadinessExpectation | ReadyReceipt | string | null} [target] An expectation, a receipt, or a generation token.
   * @param {Error | string | null} [reason]
   */
  const invalidateReadinessExpectation = (target = null, reason = null) => {
    const expectation = activeExpectation;
    if (!expectation) return false;
    if (
      target
      && target !== expectation
      && target !== expectation.generationToken
      // A string target reads `generationToken` as undefined, which is the token mismatch below.
      && /** @type {{ generationToken?: string }} */ (target)?.generationToken !== expectation.generationToken
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

  /**
   * @param {ReadinessExpectationOptions} [options]
   * @returns {ReadinessExpectation}
   */
  const registerReadinessExpectation = ({
    result,
    resultIndex = 0,
    artifactIdentity = '',
    generationToken = '',
    catalogState = null,
    phase = 'preview',
    bindingOptions = {},
    isCurrent = () => true
  } = /** @type {ReadinessExpectationOptions} */ ({})) => {
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

  /** @param {MountedResultBinderSteps} [steps] */
  const configureMountedResultBinder = (steps = {}) => {
    if (!steps || typeof steps !== 'object' || Array.isArray(steps)) {
      throw new Error('PreviewRuntime binder steps must be an object.');
    }
    bindingSteps = Object.freeze({ ...steps });
  };

  /**
   * @param {MountedResultContextOptions} [options]
   * @returns {MountedResultContext}
   */
  const createMountedResultContext = ({
    root,
    result = state.results.value[state.selectedResultIndex?.value || 0],
    resultIndex = state.selectedResultIndex?.value || 0,
    catalogState = state.featureCatalog?.value || null,
    bindingOptions = {}
  } = /** @type {MountedResultContextOptions} */ ({})) => {
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

  /**
   * @param {MountedResultContext} context
   * @returns {Promise<ReadyReceipt>}
   */
  const bindMountedResult = (context = /** @type {MountedResultContext} */ ({})) => {
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
          Promise.resolve(bindingSteps.afterReady?.(context)).catch((error) => {
            console.error('Post-ready preview work failed.', normalizeUserFacingError(error));
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

  /**
   * @param {ReadyReceipt | null | undefined} receipt
   * @param {string} [reason]
   */
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

  /**
   * @param {PreviousResultRestoreOptions} options
   * @returns {Promise<ReadyReceipt | null>}
   */
  const restorePreviousSelectedResult = async ({
    handle,
    restore,
    phase = 'rollback-restoration'
  } = /** @type {PreviousResultRestoreOptions} */ ({})) => {
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
    const rootBeforeRestore = getMountedSvg();
    await restore();
    if (!expectation) {
      clearActiveRuntime();
      recordSessionLifecycleEvent('preview.restore-bind-completed', { phase, resultIndex });
      return null;
    }
    // R10: a failure before any candidate replaced the Result restores the
    // Result that is still mounted. The mount watcher observes no change, so
    // bind that root here, as the watcher binds a remounted one.
    const root = getMountedSvg();
    if (
      root
      && root === rootBeforeRestore
      && activeExpectation === expectation
      && activeRuntime?.svg === root
      && activeRuntime.resultIdentity === expectation.resultIdentity
      && expectation.isCurrent()
    ) {
      if (activeRuntime.readyReceipt) invalidateReadyReceipt(activeRuntime.readyReceipt);
      if (state.skipCaptureBaseConfig) state.skipCaptureBaseConfig.value = false;
      void bindMountedResult(createMountedResultContext({
        root,
        result,
        resultIndex,
        catalogState: expectation.catalogState
      })).catch(() => {});
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

  /**
   * @param {string} [reason]
   * @param {string[] | null} [keys]
   */
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

  // The one write of a Result's committed content; the Result keeps its
  // committed identity.
  const writeResultContent = (resultIndex, content) => {
    const nextResults = [...state.results.value];
    nextResults[resultIndex] = {
      ...state.results.value[resultIndex],
      content
    };
    state.results.value = nextResults;
  };

  // Only commitActiveResultEdit marks the runtime dirty, and it flushes at once.
  const flushActiveResult = () => {
    const runtime = activeRuntime;
    if (!runtime?.svg) return false;
    if (!runtime.dirty) return false;

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
    if (state.skipCaptureBaseConfig) state.skipCaptureBaseConfig.value = true;
    writeResultContent(resultIndex, content);
    runtime.dirty = false;
    runtime.dirtyReasons.clear();
    return true;
  };

  const selectResult = (index) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    const count = Array.isArray(state.results.value) ? state.results.value.length : 0;
    const numeric = Number(index);
    const nextIndex = Number.isInteger(numeric) ? Math.max(0, Math.min(numeric, Math.max(0, count - 1))) : 0;
    if (state.selectedResultIndex.value === nextIndex) return false;
    const nextResult = state.results.value[nextIndex];
    const expectation = nextResult
      ? registerReadinessExpectation({
          result: nextResult,
          resultIndex: nextIndex,
          artifactIdentity: activeRuntime?.artifactIdentity || `selection:${resultRuntimeIdentity(nextResult)}`,
          generationToken: `result-selection:${nextBindSequence}`,
          catalogState: state.featureCatalog?.value || null,
          phase: 'result-selection',
          // Editor intent projected while binding replaces the Result content,
          // never its committed identity (D-07).
          isCurrent: () => (
            resultRuntimeIdentity(state.results.value[nextIndex]) === resultRuntimeIdentity(nextResult)
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
    commitActiveResultEdit(reason);
    return true;
  };

  // R1: the one commit for an editor's edit of the displayed Result's SVG.
  // Serializes the mounted root into its Result at once, so no edit waits for
  // a Result switch; unchanged content is not written.
  const commitActiveResultEdit = (reason) => {
    const runtime = activeRuntime || ensureRuntimeForCurrentSvg();
    if (!runtime?.svg) return false;
    runtime.dirty = true;
    runtime.dirtyReasons.add(String(reason || 'preview-edit'));
    return flushActiveResult();
  };

  // B17 (R1, R11): an edit of one Result's SVG. The displayed Result is
  // edited in place and committed like an editor edit. Another Result's
  // committed content is parsed, edited, and written once, so History restores
  // the Result a step was made on while a different Result is displayed.
  const commitResultEdit = (resultIndex, edit, reason = 'result-edit') => {
    const index = Number(resultIndex);
    const result = state.results.value[index];
    if (!result || typeof edit !== 'function') return false;
    const runtime = activeRuntime || ensureRuntimeForCurrentSvg();
    if (runtime?.svg && runtime.resultIndex === index) {
      return edit(runtime.svg, { mounted: true }) ? commitActiveResultEdit(reason) : false;
    }
    const Parser = globalThis.DOMParser;
    if (typeof Parser !== 'function' || typeof result.content !== 'string') return false;
    const svg = new Parser().parseFromString(result.content, 'image/svg+xml').documentElement;
    if (String(svg?.localName || '').toLowerCase() !== 'svg' || !edit(svg, { mounted: false })) return false;
    const content = serializeSvg(svg);
    if (content === result.content) return false;
    writeResultContent(index, content);
    return true;
  };

  // D-07: show the canonical editor operations on the displayed Result with
  // the executor that Generate admission uses, then persist the Result once.
  /**
   * @param {Record<string, any> | null | undefined} operations
   * @param {{ afterApply?: ((svg: SVGSVGElement) => void) | null }} [options]
   */
  const applyEditorOperations = (operations, { afterApply = null } = {}) => {
    const runtime = activeRuntime || ensureRuntimeForCurrentSvg();
    if (!runtime?.svg) return false;
    if (operations) {
      applyEditorOperationsToMountedSvg(runtime.svg, operations, { resultIndex: runtime.resultIndex });
    }
    afterApply?.(runtime.svg);
    invalidatePreviewIndexes('editor-intent-display');
    return commitActiveResultEdit('editor-intent-display');
  };

  return {
    acceptReadyReceipt,
    applyEditorOperations,
    applyFeatureVisibilityChanges,
    bindMountedResult,
    clearActiveRuntime,
    commitActiveResultEdit,
    commitResultEdit,
    configureMountedResultBinder,
    createMountedResultContext,
    getActiveRuntime,
    getFeatureElements,
    getResultIdentity: resultRuntimeIdentity,
    invalidateReadinessExpectation,
    invalidateReadyReceipt,
    invalidatePreviewIndexes,
    isActiveResultReady,
    mountResultSvg,
    registerReadinessExpectation,
    restorePreviousSelectedResult,
    selectResult
  };
};
