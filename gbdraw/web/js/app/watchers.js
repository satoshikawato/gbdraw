// @ts-check
/** @import { DrawingState, LinearSeq, SessionBusy } from '../state.js' */
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
  legendPositionOrDefault,
  normalizeCircularPlotTitlePosition,
  resolveCircularLayoutPreference
} from '../services/layout-preferences.js';
import { readFileText } from '../services/file-content-cache.js';
import { isCommittedSvgResultMounted } from '../services/svg-result-ingestion.js';
/** @import { FileLike } from '../services/session-resource-backing.js' */

/**
 * @typedef {object} RecordDiscoveryWatcherOptions
 * @property {{ readonly value: unknown } | null | undefined} [rollbackInProgress] A Session import is rolling back.
 * @property {{ readonly value: unknown } | null | undefined} [semanticWatchersSuppressed]
 * @property {{ readonly value: unknown } | null | undefined} [sessionResourceDiscoveryDeferred]
 * @property {(options: { suppress: boolean }) => unknown} refresh The record discovery; `suppress` tells it to leave the derived state alone.
 */

/** @param {RecordDiscoveryWatcherOptions} options */
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
 *   importSpecificRules: (text: string, options?: { drawing?: DrawingState, isCurrent?: () => boolean }) => Promise<number | false>
 * }} featureActions The feature editor's label reaction and the rule owner's import of a rules TSV.
 * @property {{ applyCanvasPadding: () => void, refreshDiagramDragAffordances: () => void }} legendLayout
 *   The Legend layout owner's canvas padding and drag affordances.
 * @property {{ applyPaletteDraftToPreview: () => void, syncPaletteDraftState: () => void }} resultsManager
 *   The Results owner's palette draft projections.
 * @property {(() => Promise<any>) | null} [runLabelReflow] Generate's label reflow.
 * @property {((options?: { suppress?: boolean, automatic?: boolean }) => Promise<any> | SessionBusy) | null} [refreshCircularRecordOrder]
 * @property {((options?: { suppress?: boolean }) => Promise<any>) | null} [refreshLinearRecordSelectors]
 * @property {(options?: { pan?: any, resetZoom?: boolean }) => void} [resetPreviewViewport] The preview owner's viewport reset.
 * @property {() => void} resetRightDrawer The right drawer owner's reset.
 * @property {() => void} closeLabelTextScopeDialog The label owner's port (app/feature-editor/label-actions.js).
 * @property {() => void} closeHiddenLabelTextDialog The label owner's port: closes the Label Not Shown dialog as Cancel.
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
  closeHiddenLabelTextDialog,
  clearLabelBuildNotices,
  previewRuntime,
  preparePaletteDefinitions = null
}) => {
  const {
    layoutRepositionMode,
    results,
    svgContent,
    selectedResultIndex,
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

  /** @param {unknown} value */
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

  // Each drawing's watchers watch that drawing, so a mode switch, which
  // changes neither drawing, fires none of them (R10); an edit of the drawing
  // that is not shown (a Load, a History restore) touches no shown Result.
  /** @type {DrawingState[]} */
  const drawingsOfModes = [state.drawings.circular, state.drawings.linear];
  drawingsOfModes.forEach((drawing) => watch(
    () => drawing.currentColors.value,
    () => {
      if (drawing === state.activeDrawing()) syncPaletteDraftState();
    },
    { deep: true }
  ));

  watch(
    () => paletteInstantPreviewEnabled.value,
    (enabled) => {
      const drawing = state.activeDrawing();
      if (!enabled) return;
      if (String(drawing.pendingPaletteName.value || '').trim() === '') return;
      applyPaletteDraftToPreview();
    }
  );

  drawingsOfModes.forEach((drawing) => watch(
    () => drawing.canvasPadding,
    () => {
      if (drawing !== state.activeDrawing()) return;
      if (semanticFileWatchersSuppressed.value || state.sessionOperationAvailability?.()) return;
      applyCanvasPadding();
    },
    { deep: true }
  ));

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
    () => state.drawings.circular.form.multi_record_canvas,
    (enabled, previousEnabled) => {
      const drawing = state.drawings.circular;
      if (mode.value !== 'circular') return;
      if (enabled === previousEnabled) return;

      if (enabled && !hasStoredCircularMultiRecordLayout(drawing)) {
        if (suppressCircularMultiRecordDefaults.value) {
          drawing.layoutPreferences.circular.multi.legend =
            legendPositionOrDefault(drawing.form.legend, 'circular');
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

  watch(
    () => labelReflowRequestSeq.value,
    async (nextSeq, prevSeq) => {
      if (nextSeq === prevSeq) return;
      if (!autoLabelReflowEnabled.value) return;
      if (typeof runLabelReflow !== 'function') return;
      await runLabelReflow();
    }
  );

  watch(
    () => labelReflowForceRequestSeq.value,
    async (nextSeq, prevSeq) => {
      if (nextSeq === prevSeq) return;
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
    closeHiddenLabelTextDialog();
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
  /** @type {Map<FileLike, Promise<boolean>>} */
  const pendingFileImports = window.Vue.reactive(new Map());
  /** @type {Promise<any>} */
  let fileImportApplications = Promise.resolve();
  /** @type {Map<string, FileLike | null | undefined>} */
  const restoredFileSelections = new Map();
  /**
   * Applies the text of an auxiliary file to its drawing; false rejects the file.
   * @callback FileImportApply
   * @param {string} text
   * @param {() => boolean} isCurrent
   * @returns {unknown}
   */
  /**
   * @param {string} key The `files` member the watcher follows.
   * @param {FileImportApply} apply
   * @param {FileLike | null | undefined} file
   * @param {FileLike | null | undefined} previousFile
   * @returns {Promise<boolean> | undefined}
   */
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
  /**
   * @param {string} key
   * @param {FileImportApply} apply
   */
  const watchFileImport = (key, apply) => watch(() => files[key], (file, previousFile) => applyFileImport(key, apply, file, previousFile));
  /** @param {FileLike} file */
  const waitForAuxiliaryFileImport = (file) => pendingFileImports.get(file);
  // Generate reads the drawing once every import still applying has settled,
  // including one that starts meanwhile (OV-1193).
  /** @returns {Promise<void>} */
  const settleAuxiliaryFileImports = async () => {
    while (pendingFileImports.size > 0) await Promise.allSettled([...pendingFileImports.values()]);
  };

  watchFileImport('d_color', (text) => {
    const drawing = state.activeDrawing();
    try {
      const { colors, count } = parseColorTable(text);
      drawing.currentColors.value = { ...drawing.currentColors.value, ...colors };
      console.log(`Loaded ${count} colors from file.`);
    } catch (e) {
      errorLog.value = normalizeUserFacingError(e, { stage: 'request-validation' });
      return false;
    }
  });

  watchFileImport('t_color', async (text, isCurrent) => {
    const drawing = state.activeDrawing();
    try {
      const importedCount = await featureActions.importSpecificRules(text, { drawing, isCurrent });
      if (importedCount === false) return false;
      console.log(`Loaded ${importedCount} rules from file.`);
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
        const idx = drawing.manualPriorityRules.findIndex((/** @type {{ feat: string }} */ r) => r.feat === rule.feat);
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
      // OV-923: an unchanged Result keeps its label binding.
      nextTick(() => {
        syncLabelEditor({ reuseCurrentBinding: true });
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
      ...linearSeqs.flatMap((/** @type {LinearSeq} */ seq) => [
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
  return { waitForAuxiliaryFileImport, settleAuxiliaryFileImports, auxiliaryFileImportPending: () => pendingFileImports.size > 0, canRetryAuxiliaryImportFailure, retryAuxiliaryImportFailure, resetModeTransientUi };
};
