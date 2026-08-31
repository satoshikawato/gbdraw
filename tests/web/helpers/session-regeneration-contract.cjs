'use strict';

const SUCCESS_PROFILES = new Set(['first-generate', 'second-generate']);
const FAILURE_PROFILES = new Set([
  'failure-before-activation',
  'failure-after-activation',
  'stale-completion'
]);

const number = (value) => Number.isFinite(Number(value)) ? Number(value) : 0;

const eventName = (event) => typeof event === 'string' ? event : String(event?.name || '');

const eventDetail = (event) => (
  event && typeof event === 'object' ? event : {}
);

const eventIndexes = (events) => {
  const indexes = new Map();
  events.forEach((event, index) => {
    const name = eventName(event);
    if (name && !indexes.has(name)) indexes.set(name, index);
  });
  return indexes;
};

const evaluateSessionRegenerationContract = ({
  profile,
  metrics = {},
  events = [],
  resultCount = 1,
  state = {}
} = {}) => {
  const failures = [];
  const indexes = eventIndexes(events);
  const add = (code, contract, message) => {
    failures.push(`${code} [${contract}] ${message}`);
  };
  const exact = (key, expected, code, contract) => {
    const actual = number(metrics[key]);
    if (actual !== expected) add(code, contract, `${key} must be ${expected}; received ${actual}.`);
  };
  const atMost = (key, limit, code, contract) => {
    const actual = number(metrics[key]);
    if (actual > limit) add(code, contract, `${key} must be <= ${limit}; received ${actual}.`);
  };
  const before = (earlier, later, code) => {
    const earlierIndex = indexes.get(earlier);
    const laterIndex = indexes.get(later);
    if (earlierIndex === undefined || laterIndex === undefined || earlierIndex >= laterIndex) {
      add(code, 'LIFECYCLE', `${earlier} must occur before ${later}.`);
    }
  };

  if (profile === 'saved-preview') {
    exact('workerConstructionCount', 0, 'SAVED_WORKER_CONSTRUCTION', 'SWC-WEB-REGEN-001');
    exact('workerInitializationCount', 0, 'SAVED_WORKER_INITIALIZATION', 'SWC-WEB-REGEN-001');
    exact('previewMountAdoptionCount', 1, 'SAVED_PREVIEW_ADOPTION', 'SWC-WEB-PREVIEW-001');
    exact('previewBinderInvocationCount', 1, 'SAVED_PREVIEW_BINDER', 'SWC-WEB-PREVIEW-001');
    exact('previewReadyReceiptAcceptedCount', 1, 'SAVED_PREVIEW_READY', 'SWC-WEB-PREVIEW-001');
    if (number(state.selectedResultCount) !== 1) {
      add('SAVED_SELECTED_RESULT', 'SWC-WEB-PREVIEW-001', 'Saved Session must select exactly one Result.');
    }
    if (state.currentCatalog !== true) {
      add('SAVED_CURRENT_CATALOG', 'SWC-WEB-REGEN-001', 'Saved Result must retain an aligned schema-3 catalog.');
    }
    if (state.activeDraftCommittedDistinct !== true) {
      add('SAVED_DRAFT_COMMITTED_DISTINCTION', 'SWC-WEB-REGEN-001', 'Active draft and committed artifact distinction must remain represented.');
    }
    if (state.interactiveProbePassed !== true) {
      add('SAVED_INTERACTION', 'SWC-WEB-PREVIEW-001', 'Saved preview must be interactive before Generate.');
    }
    return failures;
  }

  if (SUCCESS_PROFILES.has(profile)) {
    exact('generatedArtifactHeavyTraversalCount', 0, 'REGEN_HEAVY_HISTORY_TRAVERSAL', 'SWC-WEB-REGEN-001');
    exact('generatedArtifactFullSignatureCount', 0, 'REGEN_FULL_ARTIFACT_SIGNATURE', 'SWC-WEB-REGEN-001');
    exact('generatedArtifactMutableIntentSnapshotCount', 2, 'REGEN_MUTABLE_INTENT_SNAPSHOTS', 'SWC-WEB-REGEN-001');
    exact('rendererExecutionCount', 1, 'REGEN_RENDERER_EXECUTION', 'SWC-WEB-REGEN-001');
    exact('featureCatalogAdmissionCount', 1, 'REGEN_CATALOG_ADMISSION', 'SWC-WEB-REGEN-001');
    exact('featureCatalogSecondaryTraversalCount', 0, 'REGEN_CATALOG_SECONDARY_TRAVERSAL', 'SWC-WEB-REGEN-001');
    exact('svgSanitizationCount', number(resultCount), 'REGEN_SVG_SANITIZATION', 'SWC-WEB-REGEN-001');

    [
      ['applicationSvgParseCount', 'REGEN_EMPTY_APPLICATION_PARSE'],
      ['svgMutationIndexBuildCount', 'REGEN_EMPTY_MUTATION_INDEX'],
      ['admissionFeatureDomFullScanCount', 'REGEN_EMPTY_FEATURE_SCAN'],
      ['admissionLegendDomFullScanCount', 'REGEN_EMPTY_LEGEND_SCAN'],
      ['svgIdentityScanCount', 'REGEN_EMPTY_IDENTITY_SCAN'],
      ['svgSerializationCount', 'REGEN_EMPTY_SERIALIZATION'],
      ['currentLegacyNormalizationCount', 'REGEN_CURRENT_LEGACY_NORMALIZATION'],
      ['legacyOverrideMigrationCount', 'REGEN_LEGACY_OVERRIDE_MIGRATION'],
      ['manualRuleFeatureMatchCount', 'REGEN_MANUAL_RULE_FEATURE_PASS']
    ].forEach(([key, code]) => exact(key, 0, code, 'SWC-WEB-REGEN-001'));

    exact('generatedArtifactCandidateBuildCount', 1, 'REGEN_CANDIDATE_BUILD', 'SWC-WEB-REGEN-001');
    exact('generatedArtifactReadinessExpectationCount', 1, 'REGEN_READINESS_EXPECTATION', 'SWC-WEB-REGEN-001');
    exact('generatedArtifactActivationCount', 1, 'REGEN_ACTIVATION', 'SWC-WEB-REGEN-001');
    exact('generatedArtifactRollbackCount', 0, 'REGEN_UNEXPECTED_ROLLBACK', 'SWC-WEB-REGEN-001');
    exact('generatedArtifactFinalizeCount', 1, 'REGEN_FINALIZATION', 'SWC-WEB-REGEN-001');
    atMost('historyReplacementCount', 1, 'REGEN_HISTORY_REPLACEMENT', 'SWC-WEB-REGEN-001');

    exact('previewMaterializationObservedCount', 1, 'PREVIEW_MATERIALIZATION', 'SWC-WEB-PREVIEW-001');
    exact('previewMountAdoptionCount', 1, 'PREVIEW_MOUNT_ADOPTION', 'SWC-WEB-PREVIEW-001');
    exact('previewBinderInvocationCount', 1, 'PREVIEW_BINDER', 'SWC-WEB-PREVIEW-001');
    exact('previewDuplicateBindRejectedCount', 0, 'PREVIEW_DUPLICATE_BINDER', 'SWC-WEB-PREVIEW-001');
    exact('previewReadyReceiptEmittedCount', 1, 'PREVIEW_READY_EMITTED', 'SWC-WEB-PREVIEW-001');
    exact('previewReadyReceiptAcceptedCount', 1, 'PREVIEW_READY_ACCEPTED', 'SWC-WEB-PREVIEW-001');
    exact('previewReadyReceiptRejectedCount', 0, 'PREVIEW_READY_REJECTED', 'SWC-WEB-PREVIEW-001');
    exact('previewPostBindFrameCount', 1, 'PREVIEW_POST_BIND_FRAME', 'SWC-WEB-PREVIEW-001');
    exact('preReadyFeatureDomFullScanCount', 0, 'PREVIEW_EAGER_FEATURE_SCAN', 'SWC-WEB-PREVIEW-001');
    atMost('preReadyLegendDomFullScanCount', 1, 'PREVIEW_LEGEND_SCAN', 'SWC-WEB-PREVIEW-001');
    exact('preReadyFeatureSearchIndexBuildCount', 0, 'PREVIEW_EAGER_SEARCH_INDEX', 'SWC-WEB-PREVIEW-001');
    exact('perFeatureListenerRegistrationCount', 0, 'PREVIEW_PER_FEATURE_LISTENER', 'SWC-WEB-PREVIEW-001');
    exact('previewPollingIterationCount', 0, 'PREVIEW_POLLING_AUTHORITY', 'SWC-WEB-PREVIEW-001');

    if (profile === 'first-generate') {
      exact('workerConstructionCount', 1, 'FIRST_WORKER_CONSTRUCTION', 'SWC-WEB-REGEN-001');
      exact('workerInitializationCount', 1, 'FIRST_WORKER_INITIALIZATION', 'SWC-WEB-REGEN-001');
    } else {
      exact('workerConstructionCount', 0, 'SECOND_WORKER_CONSTRUCTION', 'SWC-WEB-REGEN-001');
      exact('workerInitializationCount', 0, 'SECOND_WORKER_INITIALIZATION', 'SWC-WEB-REGEN-001');
      exact('resourceTransferredByteCount', 0, 'SECOND_RESOURCE_RETRANSMISSION', 'SWC-WEB-REGEN-001');
    }

    [
      ['generate.processing-published', 'generate.paint-opportunity-completed', 'ORDER_PROCESSING_BEFORE_PAINT'],
      ['generate.paint-opportunity-completed', 'history.before-capture-started', 'ORDER_PAINT_BEFORE_CAPTURE'],
      ['history.before-capture-started', 'worker-post-start', 'ORDER_CAPTURE_BEFORE_WORKER'],
      ['worker-post-start', 'python-wrapper-start', 'ORDER_WORKER_BEFORE_RENDERER'],
      ['python-wrapper-start', 'python-wrapper-end', 'ORDER_RENDERER_COMPLETION'],
      ['python-wrapper-end', 'catalog.admission-started', 'ORDER_RENDERER_BEFORE_ADMISSION'],
      ['catalog.admission-started', 'svg.admission-started', 'ORDER_CATALOG_BEFORE_SVG'],
      ['svg.admission-completed', 'artifact.candidate-completed', 'ORDER_ADMISSION_BEFORE_CANDIDATE'],
      ['artifact.candidate-completed', 'preview.readiness-expectation-registered', 'ORDER_CANDIDATE_BEFORE_EXPECTATION'],
      ['preview.readiness-expectation-registered', 'artifact.activation-started', 'ORDER_EXPECTATION_BEFORE_ACTIVATION'],
      ['artifact.activation-completed', 'preview.mount-observed', 'ORDER_ACTIVATION_BEFORE_MOUNT'],
      ['preview.mount-observed', 'preview.bind-started', 'ORDER_MOUNT_BEFORE_BIND'],
      ['preview.bind-completed', 'preview.ready-receipt-emitted', 'ORDER_BIND_BEFORE_READY_EMIT'],
      ['preview.ready-receipt-emitted', 'preview.ready-receipt-accepted', 'ORDER_READY_EMIT_BEFORE_ACCEPT'],
      ['preview.ready-receipt-accepted', 'preview.post-bind-frame-completed', 'ORDER_READY_BEFORE_POST_BIND_FRAME'],
      ['preview.ready-receipt-accepted', 'generate.completed', 'SUCCESS_BEFORE_READY'],
      ['preview.post-bind-frame-completed', 'generate.processing-cleared', 'PROCESSING_CLEAR_BEFORE_POST_BIND'],
      ['preview.post-bind-frame-completed', 'history.finalization-started', 'ORDER_POST_BIND_BEFORE_HISTORY_FINALIZE'],
      ['history.finalization-completed', 'artifact.finalization-completed', 'ORDER_HISTORY_BEFORE_ARTIFACT_FINALIZE'],
      ['artifact.finalization-completed', 'generate.completed', 'ORDER_FINALIZE_BEFORE_SUCCESS'],
      ['generate.completed', 'generate.processing-cleared', 'ORDER_SUCCESS_BEFORE_PROCESSING_CLEAR']
    ].forEach(([earlier, later, code]) => before(earlier, later, code));

    if (state.terminalOutcome !== 'success') {
      add('SUCCESS_TERMINAL_OUTCOME', 'LIFECYCLE', 'Successful Generate profile requires a success terminal outcome.');
    }
    if (state.interactiveProbePassed !== true) {
      add('SUCCESS_INTERACTION', 'SWC-WEB-PREVIEW-001', 'Fresh preview must be interactive after READY.');
    }
  }

  if (FAILURE_PROFILES.has(profile)) {
    const afterActivation = profile === 'failure-after-activation';
    const stale = profile === 'stale-completion';
    const terminalNames = stale
      ? ['generate.stale-discarded']
      : ['generate.failed', 'generate.canceled'];
    const terminalIndex = events.findIndex((event) => terminalNames.includes(eventName(event)));
    const processingClearIndex = indexes.get('generate.processing-cleared');
    exact('generatedArtifactActivationCount', afterActivation ? 1 : 0, 'FAILURE_ACTIVATION_MATRIX', 'SWC-WEB-REGEN-001');
    exact('generatedArtifactRollbackCount', afterActivation ? 1 : 0, 'FAILURE_ROLLBACK_MATRIX', 'SWC-WEB-REGEN-001');
    exact('generatedArtifactFinalizeCount', 0, 'FAILURE_FINALIZATION_COMMIT', 'SWC-WEB-REGEN-001');
    exact('historyReplacementCount', 0, 'FAILURE_HISTORY_COMMIT', 'SWC-WEB-REGEN-001');
    exact('candidatePromotionPublishCount', 0, 'FAILURE_PROMOTION_COMMIT', 'SWC-WEB-REGEN-001');
    exact('previewReadyReceiptAcceptedCount', 0, stale ? 'STALE_RECEIPT_ACCEPTANCE' : 'FAILURE_RECEIPT_ACCEPTANCE', 'SWC-WEB-PREVIEW-001');
    if (terminalIndex < 0) {
      add('FAILURE_TERMINAL_EVENT_MISSING', 'LIFECYCLE', `Expected one of: ${terminalNames.join(', ')}.`);
    } else if (processingClearIndex !== undefined && terminalIndex >= processingClearIndex) {
      add('FAILURE_TERMINAL_AFTER_PROCESSING_CLEAR', 'LIFECYCLE', 'Failure, cancellation, or stale detection must precede processing clear.');
    }
    if (state.lastSuccessfulOwnerUnchanged !== true) {
      add('FAILURE_LAST_SUCCESS_CHANGED', 'SWC-WEB-REGEN-001', 'A failed, canceled, or stale candidate must not replace the last successful owner.');
    }

    if (afterActivation) {
      exact('beforeArtifactRestorationCount', 1, 'ROLLBACK_ARTIFACT_RESTORATION', 'SWC-WEB-REGEN-001');
      exact('restoredPreviewBinderInvocationCount', 1, 'ROLLBACK_RESTORED_BINDER', 'SWC-WEB-PREVIEW-001');
      exact('restoredPreviewReadyReceiptAcceptedCount', state.hadPriorReadyPreview === false ? 0 : 1, 'ROLLBACK_RESTORED_READY', 'SWC-WEB-PREVIEW-001');
      if (state.hadPriorReadyPreview !== false && state.restoredPreviewReady !== true) {
        add('ROLLBACK_READY_PREVIEW_MISSING', 'SWC-WEB-PREVIEW-001', 'Post-activation failure must restore a ready preview before processing clears.');
      }
      if (state.lastSuccessfulArtifactRestored !== true) {
        add('ROLLBACK_ARTIFACT_IDENTITY_NOT_RESTORED', 'SWC-WEB-REGEN-001', 'Rollback must restore the last successful semantic artifact identity.');
      }
      before('artifact.activation-completed', 'artifact.rollback-started', 'ORDER_ACTIVATION_BEFORE_ROLLBACK');
      if (
        terminalIndex < 0
        || indexes.get('artifact.activation-completed') >= terminalIndex
        || terminalIndex >= indexes.get('artifact.rollback-started')
      ) {
        add('ORDER_FAILURE_DETECTION_BEFORE_ROLLBACK', 'LIFECYCLE', 'Post-activation failure detection must occur before rollback starts.');
      }
      before('artifact.rollback-started', 'artifact.rollback-completed', 'ORDER_ROLLBACK_COMPLETION');
      before('artifact.rollback-completed', 'preview.restore-bind-started', 'ORDER_ROLLBACK_BEFORE_RESTORE_BIND');
      before('preview.restore-bind-completed', 'preview.restore-ready-receipt-accepted', 'ORDER_RESTORE_BIND_BEFORE_READY');
      before('preview.restore-ready-receipt-accepted', 'generate.processing-cleared', 'ORDER_RESTORED_READY_BEFORE_PROCESSING_CLEAR');
    }
    if (stale && (state.newerArtifactUnchanged !== true || state.candidateResourcesReleased !== true)) {
      add('STALE_COMPLETION_STATE', 'SWC-WEB-REGEN-001', 'Stale completion must leave the newer artifact unchanged and release candidate-local resources.');
    }
  }

  const rollbackIndex = indexes.get('artifact.rollback-started');
  const finalizationAfterRollback = rollbackIndex !== undefined && events.some((event, index) => (
    eventName(event) === 'artifact.finalization-completed' && index > rollbackIndex
  ));
  if (finalizationAfterRollback) {
    add('FINALIZATION_AFTER_ROLLBACK', 'SWC-WEB-REGEN-001', 'Candidate finalization must not occur after rollback starts.');
  }

  const staleAccepted = events.some((event) => (
    eventName(event) === 'preview.ready-receipt-accepted'
    && (eventDetail(event).stale === true
      || eventDetail(event).mismatched === true
      || eventDetail(event).incomplete === true)
  ));
  if (staleAccepted) {
    add('STALE_RECEIPT_ACCEPTANCE', 'SWC-WEB-PREVIEW-001', 'A stale, mismatched, or incomplete receipt must not be accepted.');
  }

  return failures;
};

module.exports = {
  evaluateSessionRegenerationContract
};
