import assert from 'node:assert/strict';
import test from 'node:test';
import contract from './helpers/session-regeneration-contract.cjs';

const { evaluateSessionRegenerationContract } = contract;

const successEvents = () => [
  'generate.processing-published',
  'generate.paint-opportunity-completed',
  'history.before-capture-started',
  'history.before-capture-completed',
  'worker-post-start',
  'python-wrapper-start',
  'python-wrapper-end',
  'catalog.admission-started',
  'catalog.admission-completed',
  'svg.admission-started',
  'svg.admission-completed',
  'artifact.candidate-completed',
  'preview.readiness-expectation-registered',
  'artifact.activation-started',
  'artifact.activation-completed',
  'preview.mount-observed',
  'preview.bind-started',
  'preview.bind-completed',
  'preview.ready-receipt-emitted',
  'preview.ready-receipt-accepted',
  'preview.post-bind-frame-completed',
  'history.finalization-started',
  'history.finalization-completed',
  'artifact.finalization-started',
  'artifact.finalization-completed',
  'generate.completed',
  'generate.processing-cleared'
];

const successMetrics = (profile = 'first-generate') => ({
  generatedArtifactHeavyTraversalCount: 0,
  generatedArtifactFullSignatureCount: 0,
  generatedArtifactMutableIntentSnapshotCount: 2,
  rendererExecutionCount: 1,
  featureCatalogAdmissionCount: 1,
  featureCatalogSecondaryTraversalCount: 0,
  svgSanitizationCount: 1,
  applicationSvgParseCount: 0,
  svgMutationIndexBuildCount: 0,
  admissionFeatureDomFullScanCount: 0,
  admissionLegendDomFullScanCount: 0,
  svgIdentityScanCount: 0,
  svgSerializationCount: 0,
  currentLegacyNormalizationCount: 0,
  legacyOverrideMigrationCount: 0,
  manualRuleFeatureMatchCount: 0,
  generatedArtifactCandidateBuildCount: 1,
  generatedArtifactReadinessExpectationCount: 1,
  generatedArtifactActivationCount: 1,
  generatedArtifactRollbackCount: 0,
  generatedArtifactFinalizeCount: 1,
  historyReplacementCount: 1,
  previewMaterializationObservedCount: 1,
  previewMountAdoptionCount: 1,
  previewBinderInvocationCount: 1,
  previewDuplicateBindRejectedCount: 0,
  previewReadyReceiptEmittedCount: 1,
  previewReadyReceiptAcceptedCount: 1,
  previewReadyReceiptRejectedCount: 0,
  previewPostBindFrameCount: 1,
  preReadyFeatureDomFullScanCount: 0,
  preReadyLegendDomFullScanCount: 1,
  preReadyFeatureSearchIndexBuildCount: 0,
  perFeatureListenerRegistrationCount: 0,
  previewPollingIterationCount: 0,
  workerConstructionCount: profile === 'first-generate' ? 1 : 0,
  workerInitializationCount: profile === 'first-generate' ? 1 : 0,
  resourceTransferredByteCount: 0
});

const evaluateSuccess = (profile = 'first-generate', overrides = {}, events = successEvents()) => (
  evaluateSessionRegenerationContract({
    profile,
    metrics: { ...successMetrics(profile), ...overrides },
    events,
    resultCount: 1,
    state: { terminalOutcome: 'success', interactiveProbePassed: true }
  })
);

const expectFailure = ({ code, profile = 'first-generate', overrides = {}, events, state }) => {
  const failures = evaluateSessionRegenerationContract({
    profile,
    metrics: { ...successMetrics(profile), ...overrides },
    events: events || successEvents(),
    resultCount: 1,
    state: state || { terminalOutcome: 'success', interactiveProbePassed: true }
  });
  assert.equal(failures.some((failure) => failure.startsWith(`${code} `)), true, failures.join('\n'));
};

test('accepts saved-preview, first Generate, and second Generate matrices', () => {
  assert.deepEqual(evaluateSessionRegenerationContract({
    profile: 'saved-preview',
    metrics: {
      workerConstructionCount: 0,
      workerInitializationCount: 0,
      previewMountAdoptionCount: 1,
      previewBinderInvocationCount: 1,
      previewReadyReceiptAcceptedCount: 1
    },
    state: {
      selectedResultCount: 1,
      currentCatalog: true,
      activeDraftCommittedDistinct: true,
      interactiveProbePassed: true
    }
  }), []);
  assert.deepEqual(evaluateSuccess('first-generate'), []);
  assert.deepEqual(evaluateSuccess('second-generate'), []);
});

test('accepts failure-before-activation, rollback-restoration, and stale matrices', () => {
  assert.deepEqual(evaluateSessionRegenerationContract({
    profile: 'failure-before-activation',
    metrics: {
      generatedArtifactActivationCount: 0,
      generatedArtifactRollbackCount: 0,
      generatedArtifactFinalizeCount: 0,
      historyReplacementCount: 0,
      candidatePromotionPublishCount: 0,
      previewReadyReceiptAcceptedCount: 0
    },
    events: ['generate.failed', 'generate.processing-cleared'],
    state: { lastSuccessfulOwnerUnchanged: true }
  }), []);
  assert.deepEqual(evaluateSessionRegenerationContract({
    profile: 'failure-after-activation',
    metrics: {
      generatedArtifactActivationCount: 1,
      generatedArtifactRollbackCount: 1,
      generatedArtifactFinalizeCount: 0,
      historyReplacementCount: 0,
      candidatePromotionPublishCount: 0,
      previewReadyReceiptAcceptedCount: 0,
      beforeArtifactRestorationCount: 1,
      restoredPreviewBinderInvocationCount: 1,
      restoredPreviewReadyReceiptAcceptedCount: 1
    },
    events: [
      'artifact.activation-completed',
      'generate.failed',
      'artifact.rollback-started',
      'artifact.rollback-completed',
      'preview.restore-bind-started',
      'preview.restore-bind-completed',
      'preview.restore-ready-receipt-accepted',
      'generate.processing-cleared'
    ],
    state: {
      hadPriorReadyPreview: true,
      restoredPreviewReady: true,
      lastSuccessfulOwnerUnchanged: true,
      lastSuccessfulArtifactRestored: true
    }
  }), []);
  assert.deepEqual(evaluateSessionRegenerationContract({
    profile: 'stale-completion',
    metrics: {
      generatedArtifactActivationCount: 0,
      generatedArtifactRollbackCount: 0,
      generatedArtifactFinalizeCount: 0,
      historyReplacementCount: 0,
      candidatePromotionPublishCount: 0,
      previewReadyReceiptAcceptedCount: 0
    },
    events: ['generate.stale-discarded', 'generate.processing-cleared'],
    state: {
      lastSuccessfulOwnerUnchanged: true,
      newerArtifactUnchanged: true,
      candidateResourcesReleased: true
    }
  }), []);
});

test('rejects every prohibited structural work condition', () => {
  [
    ['REGEN_HEAVY_HISTORY_TRAVERSAL', { generatedArtifactHeavyTraversalCount: 1 }],
    ['REGEN_FULL_ARTIFACT_SIGNATURE', { generatedArtifactFullSignatureCount: 1 }],
    ['REGEN_CATALOG_ADMISSION', { featureCatalogAdmissionCount: 2 }],
    ['REGEN_EMPTY_APPLICATION_PARSE', { applicationSvgParseCount: 1 }],
    ['REGEN_CURRENT_LEGACY_NORMALIZATION', { currentLegacyNormalizationCount: 1 }],
    ['REGEN_LEGACY_OVERRIDE_MIGRATION', { legacyOverrideMigrationCount: 1 }],
    ['REGEN_MANUAL_RULE_FEATURE_PASS', { manualRuleFeatureMatchCount: 1 }],
    ['PREVIEW_BINDER', { previewBinderInvocationCount: 2 }],
    ['PREVIEW_EAGER_SEARCH_INDEX', { preReadyFeatureSearchIndexBuildCount: 1 }],
    ['PREVIEW_PER_FEATURE_LISTENER', { perFeatureListenerRegistrationCount: 1 }],
    ['REGEN_ACTIVATION', { generatedArtifactActivationCount: 2 }]
  ].forEach(([code, overrides]) => expectFailure({ code, overrides }));
});

test('rejects success and processing-clear lifecycle inversions', () => {
  const successBeforeReady = successEvents();
  successBeforeReady.splice(successBeforeReady.indexOf('generate.completed'), 1);
  successBeforeReady.splice(successBeforeReady.indexOf('preview.ready-receipt-accepted'), 0, 'generate.completed');
  expectFailure({ code: 'SUCCESS_BEFORE_READY', events: successBeforeReady });

  const clearBeforeFrame = successEvents();
  clearBeforeFrame.splice(clearBeforeFrame.indexOf('generate.processing-cleared'), 1);
  clearBeforeFrame.splice(clearBeforeFrame.indexOf('preview.post-bind-frame-completed'), 0, 'generate.processing-cleared');
  expectFailure({ code: 'PROCESSING_CLEAR_BEFORE_POST_BIND', events: clearBeforeFrame });
});

test('rejects failure commits and finalization after rollback', () => {
  const failureMetrics = {
    generatedArtifactActivationCount: 0,
    generatedArtifactRollbackCount: 0,
    generatedArtifactFinalizeCount: 1,
    historyReplacementCount: 1,
    candidatePromotionPublishCount: 1,
    previewReadyReceiptAcceptedCount: 0
  };
  const failures = evaluateSessionRegenerationContract({
    profile: 'failure-before-activation',
    metrics: failureMetrics,
    events: ['generate.failed'],
    state: { terminalOutcome: 'failure' }
  });
  assert.equal(failures.some((failure) => failure.startsWith('FAILURE_FINALIZATION_COMMIT ')), true);
  assert.equal(failures.some((failure) => failure.startsWith('FAILURE_HISTORY_COMMIT ')), true);
  assert.equal(failures.some((failure) => failure.startsWith('FAILURE_PROMOTION_COMMIT ')), true);

  expectFailure({
    code: 'FINALIZATION_AFTER_ROLLBACK',
    events: [...successEvents(), 'artifact.rollback-started', 'artifact.finalization-completed']
  });
});

test('rejects second-Generate Worker, renderer, and resource reuse regressions', () => {
  expectFailure({ code: 'SECOND_WORKER_CONSTRUCTION', profile: 'second-generate', overrides: { workerConstructionCount: 1 } });
  expectFailure({ code: 'REGEN_RENDERER_EXECUTION', profile: 'second-generate', overrides: { rendererExecutionCount: 0 } });
  expectFailure({ code: 'SECOND_RESOURCE_RETRANSMISSION', profile: 'second-generate', overrides: { resourceTransferredByteCount: 1 } });
});

test('rejects stale receipt acceptance', () => {
  expectFailure({
    code: 'STALE_RECEIPT_ACCEPTANCE',
    events: successEvents().map((event) => event === 'preview.ready-receipt-accepted'
      ? { name: event, stale: true }
      : event)
  });
});

test('rejects post-activation failure without restored ready preview', () => {
  const failures = evaluateSessionRegenerationContract({
    profile: 'failure-after-activation',
    metrics: {
      generatedArtifactActivationCount: 1,
      generatedArtifactRollbackCount: 1,
      generatedArtifactFinalizeCount: 0,
      historyReplacementCount: 0,
      candidatePromotionPublishCount: 0,
      previewReadyReceiptAcceptedCount: 0,
      beforeArtifactRestorationCount: 1,
      restoredPreviewBinderInvocationCount: 1,
      restoredPreviewReadyReceiptAcceptedCount: 0
    },
    events: [
      'artifact.activation-completed',
      'artifact.rollback-started',
      'artifact.rollback-completed',
      'preview.restore-bind-started',
      'preview.restore-bind-completed',
      'generate.processing-cleared'
    ],
    state: { hadPriorReadyPreview: true, restoredPreviewReady: false }
  });
  assert.equal(failures.some((failure) => failure.startsWith('ROLLBACK_READY_PREVIEW_MISSING ')), true);
});

test('returns deterministic failures for identical inputs', () => {
  const input = {
    profile: 'second-generate',
    metrics: { ...successMetrics('second-generate'), rendererExecutionCount: 0, workerConstructionCount: 1 },
    events: successEvents(),
    resultCount: 1,
    state: { terminalOutcome: 'success', interactiveProbePassed: true }
  };
  assert.deepEqual(
    evaluateSessionRegenerationContract(input),
    evaluateSessionRegenerationContract(structuredClone(input))
  );
});
