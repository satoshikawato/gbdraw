// Shared browser operations for Node and Python acceptance.
// Assert after the same render/rollback completion boundary in either runner.

export async function settleAppRender() {
  await window.Vue.nextTick();
  await new Promise((resolveFrame) => requestAnimationFrame(resolveFrame));
  await window.Vue.nextTick();
}

export async function cancelDuringRender() {
  const app = window.__GBDRAW_APP__;
  const { state } = await import('/gbdraw/web/js/state.js');
  const {
    CANONICAL_REQUEST_SCHEMA
  } = await import('/gbdraw/web/js/services/session-request.js');
  const before = {
    proteinIdentityManifest: state.proteinIdentityManifest.value,
    legacyProteinRawCandidates: state.legacyProteinRawCandidates.value,
    legacyProteinDerivedEvidence: state.legacyProteinDerivedEvidence.value,
    losatCache: Array.from(state.losatCache.value.entries()),
    losatDerivedCache: Array.from(state.losatDerivedCache.value.entries()),
    losatCacheInfo: state.losatCacheInfo.value
  };
  const originalWorkerPostMessage = Worker.prototype.postMessage;
  let releaseRunRequest;
  const runRequestIssued = new Promise((resolve) => {
    releaseRunRequest = resolve;
  });
  let heldCanonicalRun = null;
  let cancelInvoked = false;
  Worker.prototype.postMessage = function holdFirstCanonicalRun(message, ...args) {
    if (
      !heldCanonicalRun &&
      message?.type === 'run' &&
      message?.payload?.request?.schema === CANONICAL_REQUEST_SCHEMA
    ) {
      heldCanonicalRun = { worker: this, message, args };
      releaseRunRequest();
      return;
    }
    return originalWorkerPostMessage.call(this, message, ...args);
  };
  try {
    const runPromise = app.runAnalysis();
    await runRequestIssued;
    app.cancelGeneration();
    cancelInvoked = true;
    const result = await runPromise;
    const sameMapEntries = (entries, current) => (
      entries.length === current.size &&
      entries.every(([key, value]) => current.get(key) === value)
    );
    const authorityDomains = {
      proteinIdentityManifestSame:
        state.proteinIdentityManifest.value === before.proteinIdentityManifest,
      legacyProteinRawCandidatesSame:
        state.legacyProteinRawCandidates.value === before.legacyProteinRawCandidates,
      legacyProteinDerivedEvidenceSame:
        state.legacyProteinDerivedEvidence.value === before.legacyProteinDerivedEvidence,
      losatCacheValuesSame: sameMapEntries(before.losatCache, state.losatCache.value),
      losatDerivedCacheValuesSame:
        sameMapEntries(before.losatDerivedCache, state.losatDerivedCache.value),
      losatCacheInfoSame: state.losatCacheInfo.value === before.losatCacheInfo
    };
    return {
      result,
      runRequestIssued: Boolean(heldCanonicalRun),
      cancelInvoked,
      errorSummary: String(app.errorLog?.summary || ''),
      executorCalls: Number(window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ || 0),
      ...authorityDomains,
      authorityRestored: Object.values(authorityDomains).every(Boolean)
    };
  } finally {
    Worker.prototype.postMessage = originalWorkerPostMessage;
  }
}

export async function failRendererAfterMigration() {
  const app = window.__GBDRAW_APP__;
  const { state } = await import('/gbdraw/web/js/state.js');
  const {
    CANONICAL_REQUEST_SCHEMA
  } = await import('/gbdraw/web/js/services/session-request.js');
  const before = {
    proteinIdentityManifest: state.proteinIdentityManifest.value,
    legacyProteinRawCandidates: state.legacyProteinRawCandidates.value,
    legacyProteinDerivedEvidence: state.legacyProteinDerivedEvidence.value,
    losatCache: Array.from(state.losatCache.value.entries()),
    losatDerivedCache: Array.from(state.losatDerivedCache.value.entries()),
    losatCacheInfo: state.losatCacheInfo.value
  };
  const originalWorkerPostMessage = Worker.prototype.postMessage;
  let rendererFailureInjected = false;
  Worker.prototype.postMessage = function (...args) {
    const message = args[0];
    if (
      !rendererFailureInjected &&
      message?.type === 'run' &&
      message?.payload?.request?.schema === CANONICAL_REQUEST_SCHEMA &&
      Array.isArray(message?.payload?.resourceManifest) &&
      Array.isArray(message?.payload?.stagedResources) &&
      !Object.hasOwn(message.payload, 'resources')
    ) {
      rendererFailureInjected = true;
      message.payload.request = null;
    }
    return originalWorkerPostMessage.apply(this, args);
  };
  try {
    const result = await app.runAnalysis();
    const sameMapEntries = (entries, current) => (
      entries.length === current.size &&
      entries.every(([key, value]) => current.get(key) === value)
    );
    const authorityDomains = {
      proteinIdentityManifestSame:
        state.proteinIdentityManifest.value === before.proteinIdentityManifest,
      legacyProteinRawCandidatesSame:
        state.legacyProteinRawCandidates.value === before.legacyProteinRawCandidates,
      legacyProteinDerivedEvidenceSame:
        state.legacyProteinDerivedEvidence.value === before.legacyProteinDerivedEvidence,
      losatCacheValuesSame: sameMapEntries(before.losatCache, state.losatCache.value),
      losatDerivedCacheValuesSame:
        sameMapEntries(before.losatDerivedCache, state.losatDerivedCache.value),
      losatCacheInfoSame: state.losatCacheInfo.value === before.losatCacheInfo
    };
    return {
      result,
      errorSummary: String(app.errorLog?.summary || ''),
      executorCalls: Number(window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ || 0),
      rendererFailureInjected,
      ...authorityDomains,
      authorityRestored: Object.values(authorityDomains).every(Boolean)
    };
  } finally {
    Worker.prototype.postMessage = originalWorkerPostMessage;
  }
}
