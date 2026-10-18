// @ts-check
/** @import { DrawingState, LinearSeq, SessionBusy } from '../state.js' */
/** @import { MountedResultContext, MountedResultContextOptions } from './preview-runtime.js' */
/** @import { RulePreparation } from './rule-matching.js' */
import { diagnosticError, normalizeUserFacingError } from '../utils/error-normalization.js';
import {
  parseBlacklistWords,
  parseColorTable,
  parsePriorityRules,
  parseWhitelistRules
} from '../services/file-imports.js';
import { readFileText } from '../services/file-content-cache.js';
import { isCommittedSvgResultMounted } from '../services/svg-result-ingestion.js';
import { addWhitelistRules, appendBlacklistWords, setPriorityRules } from '../services/drawing/label-tables.js';
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
 * The auxiliary file imports still applying (OV-1193). The file watchers
 * track them; Save and Load read whether one is pending, and Generate waits
 * for them. The composition root creates it before all of them.
 * @typedef {object} AuxiliaryFileImports
 * @property {(file: FileLike, pending: Promise<boolean>) => void} track Holds `file`'s import until it settles.
 * @property {(file: FileLike) => Promise<boolean> | undefined} wait The import of `file`, while it applies.
 * @property {() => boolean} isPending
 * @property {() => Promise<void>} settle Waits until no import is pending, including one that starts meanwhile.
 */

/**
 * @param {<T extends object>(target: T) => T} reactive Vue `reactive`, so the readers' computeds follow the imports.
 * @returns {AuxiliaryFileImports}
 */
export const createAuxiliaryFileImports = (reactive) => {
  /** @type {Map<FileLike, Promise<boolean>>} */
  const pending = reactive(new Map());
  return {
    track: (file, promise) => {
      pending.set(file, promise);
      void promise.finally(() => {
        if (pending.get(file) === promise) pending.delete(file);
      });
    },
    wait: (file) => pending.get(file),
    isPending: () => pending.size > 0,
    settle: async () => {
      while (pending.size > 0) await Promise.allSettled([...pending.values()]);
    }
  };
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
 * @property {AuxiliaryFileImports['track']} trackAuxiliaryFileImport Holds an import until it settles (`createAuxiliaryFileImports`).
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
 * @property {{ refreshDiagramDragAffordances: () => void }} legendLayout
 *   The Legend layout owner's drag affordances.
 * @property {{ editDefaultColors: (drawing: DrawingState, colors: Record<string, string>) => void }} resultsManager
 *   The Results owner's palette draft projection and default color edit.
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
  trackAuxiliaryFileImport,
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
  const { refreshDiagramDragAffordances } = legendLayout;
  const { editDefaultColors } = resultsManager;

  watch(
    () => layoutRepositionMode.value,
    () => {
      nextTick(() => {
        refreshLegendDragAffordances();
        refreshDiagramDragAffordances();
      });
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
  /** @type {Promise<any>} */
  let fileImportApplications = Promise.resolve();
  /**
   * Applies the text of an auxiliary file to the drawing it was chosen in;
   * false rejects the file.
   * @callback FileImportApply
   * @param {DrawingState} drawing
   * @param {string} text
   * @param {() => boolean} isCurrent
   * @returns {unknown}
   */
  /**
   * A parsed table, as its slice writer takes it; a parse error rejects the file.
   * @template T
   * @param {(text: string) => T} parse
   * @param {(drawing: DrawingState, table: T) => void} write
   * @returns {FileImportApply}
   */
  const tableImport = (parse, write) => (drawing, text) => {
    try {
      write(drawing, parse(text));
    } catch (e) {
      errorLog.value = normalizeUserFacingError(e, { stage: 'request-validation' });
      return false;
    }
  };
  /** @type {Record<string, FileImportApply>} */
  const AUXILIARY_TABLE_IMPORTS = {
    d_color: tableImport(parseColorTable, (drawing, { colors }) => editDefaultColors(drawing, colors)),
    t_color: async (drawing, text, isCurrent) => {
      try {
        const importedCount = await featureActions.importSpecificRules(text, { drawing, isCurrent });
        if (importedCount === false) return false;
      } catch (e) {
        if (isCurrent()) errorLog.value = normalizeUserFacingError(e, { operation: 'evaluateRules', stage: 'rule-validation' });
        return false;
      }
    },
    qualifier_priority: tableImport(parsePriorityRules, (drawing, { rules }) => setPriorityRules(drawing.slices.labelTables, rules)),
    whitelist: tableImport(parseWhitelistRules, (drawing, { rules }) => addWhitelistRules(drawing.slices.labelTables, rules)),
    blacklist: tableImport(parseBlacklistWords, (drawing, { words }) => {
      if (words.length > 0) appendBlacklistWords(drawing.slices.labelTables, words);
    })
  };
  /**
   * Reads `file` and applies it to the drawing shown now; a rejected file
   * gives the input its previous file back.
   * @param {string} key The `files` member of the auxiliary table.
   * @param {FileLike | null | undefined} file
   * @param {FileLike | null | undefined} previousFile
   * @returns {Promise<boolean> | undefined}
   */
  const applyFileImport = (key, file, previousFile) => {
    if (semanticFileWatchersSuppressed.value || !file) return;
    return startFileImport(key, file, previousFile);
  };
  /**
   * @param {string} key
   * @param {FileLike} file
   * @param {FileLike | null | undefined} previousFile
   * @returns {Promise<boolean>}
   */
  const startFileImport = (key, file, previousFile) => {
    const drawing = state.activeDrawing();
    const apply = AUXILIARY_TABLE_IMPORTS[key];
    // OV-1205, OV-1206: the import is current while its file stays selected.
    // The Specific rules owner (`runSpecificRuleWrite`) prepares and retries
    // the rules on the inputs current when it runs; one it gives up on comes
    // back as false with no error, and is reported here with Retry.
    const isCurrent = () => files[key] === file && !semanticFileWatchersSuppressed.value
      && !state.sessionOperationAvailability?.();
    const previousError = errorLog.value;
    // A Specific rules failure found on given inputs offers Retry while they
    // stay; one about inputs that kept changing offers it on any inputs, since
    // Retry is meant to run on the new ones.
    const retainFailure = (inputsBound = key === 't_color') => {
      auxiliaryImportFailure.value = { error: errorLog.value, key, selection: files[key], snapshot: inputsBound ? rulePreparation.snapshot() : null,
        retry: () => {
          files[key] = file;
          return startFileImport(key, file, previousFile);
        } };
    };
    const pending = readFileText(file).then((text) => {
      // Reads may finish out of order; serialize only their live application.
      const application = fileImportApplications.then(async () => {
        if (!isCurrent()) return false;
        if (await apply(drawing, text, isCurrent) === false && isCurrent()) {
          const gaveUp = key === 't_color' && errorLog.value === previousError;
          if (gaveUp) {
            errorLog.value = normalizeUserFacingError(diagnosticError('RULE_INPUTS_CHANGED', {},
              { operation: 'evaluateRules', stage: 'rule-validation' }));
          }
          files[key] = previousFile;
          if (errorLog.value !== previousError) retainFailure(key === 't_color' && !gaveUp);
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
        files[key] = previousFile;
        retainFailure();
      }
      return false;
    });
    trackAuxiliaryFileImport(file, pending);
    return pending;
  };
  // The file input of an auxiliary table (Default colors, Specific rules,
  // Priority, Whitelist, Blacklist): it stores the file and imports it into
  // the shown drawing. The `file-uploader` owns its History step and waits
  // for the import (`waitForAuxiliaryFileImport`, R11). A Load or a History
  // restore sets `files` without importing.
  /**
   * @param {string} key
   * @param {FileLike | null | undefined} file
   */
  const importAuxiliaryTable = (key, file) => {
    const previousFile = files[key];
    files[key] = file;
    return applyFileImport(key, file, previousFile);
  };

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
  return { importAuxiliaryTable, canRetryAuxiliaryImportFailure, retryAuxiliaryImportFailure, resetModeTransientUi };
};
