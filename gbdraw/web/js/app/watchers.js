// @ts-check
/** @import { DrawingState } from '../state.js' */
/** @import { MountedResultContext, MountedResultContextOptions } from './preview-runtime.js' */
/** @import { RulePreparation } from './rule-matching.js' */
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import {
  parseBlacklistWords,
  parseColorTable,
  parsePriorityRules,
  parseWhitelistRules
} from '../services/file-imports.js';
import {
  prepareSpecificColorImport
} from '../services/specific-color-rules.js';
import {
  normalizeCircularPlotTitlePosition,
  resolveCircularLayoutPreference
} from '../services/layout-preferences.js';
import { readFileText } from '../services/file-content-cache.js';
import { isCommittedSvgResultMounted } from '../services/svg-result-ingestion.js';

export const runRecordDiscoveryWatcher = async ({
  rollbackInProgress,
  semanticWatchersSuppressed,
  sessionResourceDiscoveryDeferred,
  refresh
}) => {
  const suppress = Boolean(
    rollbackInProgress?.value
    || semanticWatchersSuppressed?.value
    || sessionResourceDiscoveryDeferred?.value
  );
  await refresh({ suppress });
  return !suppress;
};

/**
 * The members of the preview owner's runtime (app/preview-runtime.js) that the
 * mounted-Result watcher reads.
 * @typedef {object} WatchersPreviewRuntime
 * @property {() => void} clearActiveRuntime
 * @property {(options: MountedResultContextOptions) => MountedResultContext} createMountedResultContext
 * @property {(context: MountedResultContext) => Promise<any>} bindMountedResult
 */

/**
 * @typedef {object} WatchersOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {RulePreparation} rulePreparation
 * @property {(value?: any) => { value: any }} ref Vue `ref`.
 * @property {(getter: () => any) => { readonly value: any }} computed Vue `computed`.
 * @property {(source: any, callback: (...args: any[]) => any, options?: Record<string, any>) => any} watch Vue `watch`.
 * @property {(callback?: () => void) => Promise<void>} nextTick Vue `nextTick`.
 * @property {(callback: () => any) => void} onMounted Vue `onMounted`.
 * @property {{ refreshLegendDragAffordances: () => void }} legendActions
 *   The Legend owner's reaction to a rebuilt Legend.
 * @property {{
 *   syncLabelEditor: (options?: Record<string, any>) => void,
 *   commitSpecificRules: (rules: Record<string, any>[], label?: string, options?: Record<string, any>) => Promise<any>
 * }} featureActions The feature editor's label reaction and its commit of the specific color rules.
 * @property {{ applyCanvasPadding: () => void, refreshDiagramDragAffordances: () => void }} legendLayout
 *   The Legend layout owner's canvas padding and drag affordances.
 * @property {{ applyPaletteDraftToPreview: () => void, syncPaletteDraftState: () => void }} resultsManager
 *   The Results owner's palette draft projections.
 * @property {(() => Promise<any>) | null} [runLabelReflow] Generate's label reflow.
 * @property {((options?: { suppress?: boolean, automatic?: boolean }) => Promise<any>) | null} [refreshCircularRecordOrder]
 * @property {((options?: { suppress?: boolean }) => Promise<any>) | null} [refreshLinearRecordSelectors]
 * @property {(options?: { pan?: any, resetZoom?: boolean }) => void} [resetPreviewViewport] The preview owner's viewport reset.
 * @property {() => void} resetRightDrawer The right drawer owner's reset.
 * @property {() => void} closeLabelTextScopeDialog The label owner's port (app/feature-editor/label-actions.js).
 * @property {(options?: { rerender?: boolean }) => void} clearLabelBuildNotices The label owner's port.
 * @property {WatchersPreviewRuntime} previewRuntime The preview owner's mounted-Result binding.
 * @property {(() => Promise<any>) | null} [preparePaletteDefinitions] The palette loader's load of the browser palette definitions.
 */

/** @param {WatchersOptions} options */
export const setupWatchers = ({
  state,
  rulePreparation,
  ref, computed, watch,
  nextTick,
  onMounted,
  legendActions,
  featureActions,
  legendLayout,
  resultsManager,
  runLabelReflow,
  refreshCircularRecordOrder,
  refreshLinearRecordSelectors,
  resetPreviewViewport,
  resetRightDrawer,
  // R13: the label owner's ports (app/feature-editor/label-actions.js).
  closeLabelTextScopeDialog,
  clearLabelBuildNotices,
  previewRuntime,
  preparePaletteDefinitions = null
}) => {
  const {
    layoutRepositionMode,
    results,
    svgContent,
    selectedResultIndex,
    shouldDeferCircularPreviewUpdates,
    mode,
    cInputType,
    lInputType,
    skipCaptureBaseConfig,
    skipExtractOnSvgChange,
    trustedArtifactRestoreInProgress,
    svgContainer,
    suppressCircularMultiRecordDefaults,
    selectedFeatureRecordIdx,
    featurePanelTab,
    labelSearch,
    selectedOrthogroupAlignmentFeature,
    selectedOrthogroupId,
    orthogroupSearch,
    isFeatureDrawerMounted,
    clickedPairwiseMatch,
    clickedLabel,
    hiddenLabelTextDialog,
    files,
    paletteInstantPreviewEnabled,
    semanticFileWatchersSuppressed,
    sessionResourceDiscoveryDeferred,
    sessionImportRollbackInProgress,
    linearSeqs,
    linearReorderNotice,
    autoLabelReflowEnabled,
    labelReflowRequestSeq,
    labelReflowForceRequestSeq,
    errorLog
  } = state;

  const { refreshLegendDragAffordances } = legendActions;

  const { syncLabelEditor } = featureActions;
  const {
    applyCanvasPadding,
    refreshDiagramDragAffordances
  } = legendLayout;
  const {
    applyPaletteDraftToPreview,
    syncPaletteDraftState
  } = resultsManager;

  const normalizeLegendPosition = (value, fallback = 'left') => {
    const normalized = String(value || '').trim().toLowerCase();
    return normalized || fallback;
  };

  const hasStoredLayoutValue = (value) => typeof value === 'string' && value.trim() !== '';

  /** @param {DrawingState} drawing */
  const hasStoredCircularMultiRecordLayout = (drawing) =>
    hasStoredLayoutValue(drawing.layoutPreferences.circular.multi.legend) ||
    hasStoredLayoutValue(drawing.layoutPreferences.circular.multi.plotTitlePosition);

  /** @param {DrawingState} drawing */
  const applyCircularMultiRecordSmartDefaults = (drawing) => {
    const singleLayout = resolveCircularLayoutPreference(drawing.layoutPreferences, false);
    drawing.layoutPreferences.circular.multi.legend =
      singleLayout.legend === 'left' ? 'bottom' : singleLayout.legend;
    drawing.layoutPreferences.circular.multi.plotTitlePosition =
      singleLayout.plotTitlePosition === 'none' ? 'bottom' : singleLayout.plotTitlePosition;
  };

  watch(
    () => state.activeDrawing().currentColors.value,
    () => {
      syncPaletteDraftState();
    },
    { deep: true }
  );

  watch(
    () => paletteInstantPreviewEnabled.value,
    (enabled) => {
      const drawing = state.activeDrawing();
      if (!enabled) return;
      if (String(drawing.pendingPaletteName.value || '').trim() === '') return;
      applyPaletteDraftToPreview();
    }
  );

  watch(
    () => state.activeDrawing().canvasPadding,
    () => {
      if (semanticFileWatchersSuppressed.value || state.sessionOperationAvailability?.()) return;
      applyCanvasPadding();
    },
    { deep: true }
  );

  watch(
    () => layoutRepositionMode.value,
    () => {
      nextTick(() => {
        refreshLegendDragAffordances();
        refreshDiagramDragAffordances();
      });
    }
  );

  watch(
    () => state.activeDrawing().form.multi_record_canvas,
    (enabled, previousEnabled) => {
      const drawing = state.activeDrawing();
      if (mode.value !== 'circular') return;
      if (enabled === previousEnabled) return;

      if (enabled && !hasStoredCircularMultiRecordLayout(drawing)) {
        if (suppressCircularMultiRecordDefaults.value) {
          drawing.layoutPreferences.circular.multi.legend = normalizeLegendPosition(
            drawing.form.legend,
            'left'
          );
          drawing.layoutPreferences.circular.multi.plotTitlePosition =
            normalizeCircularPlotTitlePosition(drawing.adv.plot_title_position);
        } else {
          applyCircularMultiRecordSmartDefaults(drawing);
        }
      }

      if (suppressCircularMultiRecordDefaults.value) {
        suppressCircularMultiRecordDefaults.value = false;
      }

    }
  );

  // Batch Result and template-ref changes after the replacement root is mounted.
  watch([svgContent, svgContainer, () => results.value[selectedResultIndex.value]], () => {
    const isIncrementalEdit = Boolean(skipCaptureBaseConfig.value);
    skipCaptureBaseConfig.value = false;

    nextTick(async () => {
      const root = svgContainer.value?.querySelector('svg') || null;
      if (!root) {
        previewRuntime.clearActiveRuntime();
        return;
      }
      const resultIndex = Number(selectedResultIndex.value) || 0;
      const result = results.value[resultIndex] || null;
      try {
        const context = previewRuntime.createMountedResultContext({
          root,
          result,
          resultIndex,
          catalogState: state.featureCatalog?.value || null,
          bindingOptions: {
            isIncrementalEdit,
            skipLegendExtraction: Boolean(skipExtractOnSvgChange.value),
            trustedRestore: Boolean(trustedArtifactRestoreInProgress.value)
          }
        });
        await previewRuntime.bindMountedResult(context);
      } catch (error) {
        if ([
          'PREVIEW_BIND_STALE',
          'PREVIEW_BIND_SUPERSEDED',
          'PREVIEW_ROOT_MISMATCH'
        ].includes(error?.code)) return;
        errorLog.value = normalizeUserFacingError(error, { operation: 'generate', stage: 'result-admission' });
      }
    });
  }, { flush: 'post' });

  // Persisting the current live DOM changes Result text but deliberately leaves the
  // mounted root in place. Consume the old remount-only flags at that boundary.
  watch(
    () => results.value[selectedResultIndex.value]?.content,
    () => {
      const result = results.value[selectedResultIndex.value];
      if (!isCommittedSvgResultMounted(result)) return;
      skipCaptureBaseConfig.value = false;
    }
  );

  // The saved label table is the rules a loaded Session sent; a bulk label
  // edit replaces it with the rows the editor builds. Per-feature label edits
  // are identity rows, which apply before the table (design Q4).
  watch(
    () => state.activeDrawing().labelTextBulkOverrides,
    () => {
      const drawing = state.activeDrawing();
      if (semanticFileWatchersSuppressed.value) return;
      drawing.canonicalLabelOverrideRows.value = [];
    },
    { deep: true }
  );

  watch(
    () => labelReflowRequestSeq.value,
    async (nextSeq, prevSeq) => {
      if (nextSeq === prevSeq) return;
      if (!autoLabelReflowEnabled.value) return;
      if (mode.value === 'circular' && shouldDeferCircularPreviewUpdates.value) return;
      if (typeof runLabelReflow !== 'function') return;
      await runLabelReflow();
    }
  );

  watch(
    () => labelReflowForceRequestSeq.value,
    async (nextSeq, prevSeq) => {
      if (nextSeq === prevSeq) return;
      if (mode.value === 'circular' && shouldDeferCircularPreviewUpdates.value) return;
      if (typeof runLabelReflow !== 'function') return;
      await runLabelReflow();
    }
  );

  // E1: the transient UI of the departing mode. The mode transition in the
  // composition root calls it after it installed the arriving mode's artifact
  // (R10); durable label and visibility intent stays (R2).
  const resetModeTransientUi = () => {
    if (typeof resetPreviewViewport === 'function') resetPreviewViewport();
    selectedFeatureRecordIdx.value = 0;
    selectedOrthogroupAlignmentFeature.value = '';
    selectedOrthogroupId.value = '';
    orthogroupSearch.value = '';
    clearLabelBuildNotices();
    labelSearch.value = '';
    featurePanelTab.value = 'colors';
    clickedPairwiseMatch.value = null;
    clickedLabel.value = null;
    closeLabelTextScopeDialog();
    hiddenLabelTextDialog.show = false;
    hiddenLabelTextDialog.featureId = '';
    hiddenLabelTextDialog.reason = '';
    resetRightDrawer();
    linearReorderNotice.value = '';
  };

  const auxiliaryImportFailure = ref(null);
  const canRetryAuxiliaryImportFailure = computed(() => Boolean(auxiliaryImportFailure.value
    && errorLog.value === auxiliaryImportFailure.value.error
    && files[auxiliaryImportFailure.value.key] === auxiliaryImportFailure.value.selection
    && (!auxiliaryImportFailure.value.snapshot || rulePreparation.isCurrent(auxiliaryImportFailure.value.snapshot))));
  const retryAuxiliaryImportFailure = () => {
    const busy = state.sessionOperationAvailability?.();
    return busy || (canRetryAuxiliaryImportFailure.value ? auxiliaryImportFailure.value.retry() : false);
  };
  const pendingFileImports = window.Vue.reactive(new Map());
  /** @type {Promise<any>} */
  let fileImportApplications = Promise.resolve();
  const restoredFileSelections = new Map();
  const applyFileImport = (key, apply, file, previousFile) => {
    const restored = restoredFileSelections.has(key) && restoredFileSelections.get(key) === file;
    restoredFileSelections.delete(key);
    if (restored) return;
    if (semanticFileWatchersSuppressed.value || !file) return;
    const ruleContext = key === 't_color' ? rulePreparation.snapshot() : null;
    const isCurrent = () => files[key] === file && !semanticFileWatchersSuppressed.value
      && !state.sessionOperationAvailability?.()
      && (!ruleContext || rulePreparation.isCurrent(ruleContext));
    const previousError = errorLog.value;
    const retainFailure = () => {
      auxiliaryImportFailure.value = { error: errorLog.value, key, selection: files[key], snapshot: ruleContext ? rulePreparation.snapshot() : null,
        retry: async () => {
          if (files[key] === file) return applyFileImport(key, apply, file, previousFile);
          files[key] = file;
          await nextTick();
          return pendingFileImports.get(file);
        } };
    };
    const pending = readFileText(file).then((text) => {
      // Reads may finish out of order; serialize only their live application.
      const application = fileImportApplications.then(async () => {
        if (!isCurrent()) return false;
        if (await apply(text, isCurrent) === false && isCurrent()) {
          restoredFileSelections.set(key, previousFile);
          files[key] = previousFile;
          if (errorLog.value !== previousError) retainFailure();
          return false;
        } else if (files[key] === file && errorLog.value === auxiliaryImportFailure.value?.error) {
          errorLog.value = null;
          auxiliaryImportFailure.value = null;
        }
        return true;
      });
      fileImportApplications = application.catch(() => {});
      return application;
    }).catch((error) => {
      if (isCurrent() && errorLog.value === previousError) {
        errorLog.value = normalizeUserFacingError(error, { operation: key === 't_color' ? 'evaluateRules' : 'unknown', stage: 'resource-staging' });
        restoredFileSelections.set(key, previousFile);
        files[key] = previousFile;
        retainFailure();
      }
      return false;
    });
    pendingFileImports.set(file, pending);
    void pending.finally(() => {
      if (pendingFileImports.get(file) === pending) pendingFileImports.delete(file);
    });
    return pending;
  };
  const watchFileImport = (key, apply) => watch(() => files[key], (file, previousFile) => applyFileImport(key, apply, file, previousFile));
  const waitForAuxiliaryFileImport = (file) => pendingFileImports.get(file);

  watchFileImport('d_color', (text) => {
    const drawing = state.activeDrawing();
    try {
      const { colors, count } = parseColorTable(text);
      Object.entries(colors).forEach(([key, color]) => {
        drawing.currentColors.value[key] = color;
      });
      console.log(`Loaded ${count} colors from file.`);
    } catch (e) {
      errorLog.value = normalizeUserFacingError(e, { stage: 'request-validation' });
      return false;
    }
  });

  watchFileImport('t_color', async (text, isCurrent) => {
    const drawing = state.activeDrawing();
    try {
      const prepared = prepareSpecificColorImport(text, drawing.manualSpecificRules);
      if (!await featureActions.commitSpecificRules(prepared.nextRules, 'Import specific color rules', { isCurrent })) return false;
      console.log(`Loaded ${prepared.importedCount} rules from file.`);
    } catch (e) {
      if (isCurrent()) errorLog.value = normalizeUserFacingError(e, { operation: 'evaluateRules', stage: 'rule-validation' });
      return false;
    }
  });

  watchFileImport('qualifier_priority', (text) => {
    const drawing = state.activeDrawing();
    try {
      const { rules, count } = parsePriorityRules(text);
      rules.forEach((rule) => {
        const idx = drawing.manualPriorityRules.findIndex((r) => r.feat === rule.feat);
        if (idx >= 0) {
          drawing.manualPriorityRules[idx].order = rule.order;
        } else {
          drawing.manualPriorityRules.push({ feat: rule.feat, order: rule.order });
        }
      });
      console.log(`Loaded ${count} priority rules.`);
    } catch (e) {
      errorLog.value = normalizeUserFacingError(e, { stage: 'request-validation' });
      return false;
    }
  });

  watchFileImport('whitelist', (text) => {
    const drawing = state.activeDrawing();
    try {
      const { rules, count } = parseWhitelistRules(text);
      rules.forEach((rule) => drawing.manualWhitelist.push(rule));
      console.log(`Loaded ${count} whitelist rules.`);
    } catch (e) {
      errorLog.value = normalizeUserFacingError(e, { stage: 'request-validation' });
      return false;
    }
  });

  watchFileImport('blacklist', (text) => {
    const drawing = state.activeDrawing();
    try {
      const { words, count } = parseBlacklistWords(text);
      if (words.length > 0) {
        const existing = drawing.manualBlacklist.value ? drawing.manualBlacklist.value.trim() : '';
        const separator = existing && !existing.endsWith(',') ? ', ' : '';
        drawing.manualBlacklist.value = existing + separator + words.join(', ');
        console.log(`Loaded ${count} blacklist words.`);
      }
    } catch (e) {
      errorLog.value = normalizeUserFacingError(e, { stage: 'request-validation' });
      return false;
    }
  });

  watch(
    () => isFeatureDrawerMounted.value,
    (visible) => {
      if (!visible) return;
      nextTick(() => {
        syncLabelEditor();
      });
    }
  );

  watch(
    () => [
      semanticFileWatchersSuppressed.value,
      state.sessionSavePending?.value,
      state.sessionImportPending?.value,
      mode.value,
      cInputType.value,
      files.c_gb,
      files.c_gff,
      files.c_fasta
    ],
    async () => {
      if (typeof refreshCircularRecordOrder !== 'function') return;
      await runRecordDiscoveryWatcher({
        rollbackInProgress: sessionImportRollbackInProgress,
        semanticWatchersSuppressed: semanticFileWatchersSuppressed,
        sessionResourceDiscoveryDeferred,
        refresh: ({ suppress }) => refreshCircularRecordOrder({ suppress, automatic: true })
      });
    }
  );
  watch(
    () => [
      semanticFileWatchersSuppressed.value,
      state.sessionSavePending?.value,
      state.sessionImportPending?.value,
      mode.value,
      lInputType.value,
      ...linearSeqs.flatMap((seq) => [
        seq.uid,
        lInputType.value === 'gff' ? seq.gff : seq.gb,
        lInputType.value === 'gff' ? seq.fasta : null
      ])
    ],
    async () => {
      if (typeof refreshLinearRecordSelectors !== 'function') return;
      await runRecordDiscoveryWatcher({
        rollbackInProgress: sessionImportRollbackInProgress,
        semanticWatchersSuppressed: semanticFileWatchersSuppressed,
        sessionResourceDiscoveryDeferred,
        refresh: refreshLinearRecordSelectors
      });
    },
    { immediate: true }
  );

  onMounted(async () => {
    await nextTick();
    if (typeof preparePaletteDefinitions !== 'function') return;
    try {
      await preparePaletteDefinitions();
    } catch (error) {
      console.warn('Could not load browser palette definitions.', normalizeUserFacingError(error, { stage: 'initialization' }));
    }
  });
  return { waitForAuxiliaryFileImport, auxiliaryFileImportPending: () => pendingFileImports.size > 0, canRetryAuxiliaryImportFailure, retryAuxiliaryImportFailure, resetModeTransientUi };
};
