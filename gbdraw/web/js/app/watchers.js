import { normalizeUserFacingError } from '../services/error-normalization.js';
import {
  parseBlacklistWords,
  parseColorTable,
  parsePriorityRules,
  parseWhitelistRules
} from './file-imports.js';
import {
  prepareSpecificColorImport
} from './specific-color-rules.js';
import {
  normalizeCircularPlotTitlePosition
} from './plot-title-position.js';
import { resolveCircularLayoutPreference } from './layout-preferences.js';
import { readFileText } from '../services/file-content-cache.js';
import { isCommittedSvgResultMounted } from '../services/svg-result-ingestion.js';
import { featureStateFromCatalog } from '../services/feature-catalog.js';

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

export const setupWatchers = ({
  state,
  rulePreparation,
  ref, computed, watch,
  nextTick,
  onMounted,
  legendActions,
  svgActions,
  featureActions,
  legendLayout,
  resultsManager,
  runLabelReflow,
  refreshCircularRecordOrder,
  refreshLinearRecordSelectors,
  resetPreviewViewport,
  resetRightDrawer,
  previewRuntime = null,
  preparePaletteDefinitions = null
}) => {
  const {
    manualSpecificRules,
    extractedFeatures,
    biologicalFeatures,
    layoutRepositionMode,
    editableLabels,
    results,
    svgContent,
    selectedResultIndex,
    form,
    generatedMode,
    shouldDeferCircularPreviewUpdates,
    mode,
    cInputType,
    lInputType,
    canvasPadding,
    skipCaptureBaseConfig,
    skipExtractOnSvgChange,
    trustedArtifactRestoreInProgress,
    svgContainer,
    layoutPreferences,
    suppressCircularMultiRecordDefaults,
    featureRecordIds,
    selectedFeatureRecordIdx,
    featurePanelTab,
    labelSearch,
    orthogroups,
    collinearGroups,
    featureOrthogroupIndex,
    selectedOrthogroupAlignmentFeature,
    selectedOrthogroupId,
    orthogroupSearch,
    labelTextBulkOverrides,
    canonicalLabelOverrideRows,
    labelOverrideBuildWarning,
    isFeatureDrawerMounted,
    clickedFeature,
    clickedPairwiseMatch,
    clickedLabel,
    labelTextScopeDialog,
    hiddenLabelTextDialog,
    files,
    currentColors,
    paletteInstantPreviewEnabled,
    pendingPaletteName,
    semanticFileWatchersSuppressed,
    sessionResourceDiscoveryDeferred,
    sessionImportRollbackInProgress,
    manualPriorityRules,
    manualWhitelist,
    manualBlacklist,
    linearSeqs,
    linearReorderNotice,
    autoLabelReflowEnabled,
    labelReflowRequestSeq,
    labelReflowForceRequestSeq,
    errorLog
  } = state;

  const {
    addLegendEntry,
    extractLegendEntries,
    refreshLegendDragAffordances
  } = legendActions;

  const { applyPaletteToSvg } = svgActions;
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

  const hasStoredCircularMultiRecordLayout = () =>
    hasStoredLayoutValue(layoutPreferences.circular.multi.legend) ||
    hasStoredLayoutValue(layoutPreferences.circular.multi.plotTitlePosition);

  const applyCircularMultiRecordSmartDefaults = () => {
    const singleLayout = resolveCircularLayoutPreference(layoutPreferences, false);
    layoutPreferences.circular.multi.legend =
      singleLayout.legend === 'left' ? 'bottom' : singleLayout.legend;
    layoutPreferences.circular.multi.plotTitlePosition =
      singleLayout.plotTitlePosition === 'none' ? 'bottom' : singleLayout.plotTitlePosition;
  };

  watch(
    currentColors,
    () => {
      syncPaletteDraftState();
    },
    { deep: true }
  );

  watch(
    () => paletteInstantPreviewEnabled.value,
    (enabled) => {
      if (!enabled) return;
      if (String(pendingPaletteName.value || '').trim() === '') return;
      applyPaletteDraftToPreview();
    }
  );

  watch(
    canvasPadding,
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
    () => form.multi_record_canvas,
    (enabled, previousEnabled) => {
      if (mode.value !== 'circular') return;
      if (enabled === previousEnabled) return;

      if (enabled && !hasStoredCircularMultiRecordLayout()) {
        if (suppressCircularMultiRecordDefaults.value) {
          layoutPreferences.circular.multi.legend = normalizeLegendPosition(
            form.legend,
            'left'
          );
          layoutPreferences.circular.multi.plotTitlePosition =
            normalizeCircularPlotTitlePosition(state.adv.plot_title_position);
        } else {
          applyCircularMultiRecordSmartDefaults();
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
        previewRuntime?.clearActiveRuntime?.();
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
    labelTextBulkOverrides,
    () => {
      if (semanticFileWatchersSuppressed.value) return;
      canonicalLabelOverrideRows.value = [];
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

  watch(
    () => mode.value,
    () => {
      if (semanticFileWatchersSuppressed.value) return;

      // Vue replaces the mode-keyed container. Release its frozen live-edit
      // payload first so the new root materializes the selected Result content.
      previewRuntime?.clearActiveRuntime?.();

      if (typeof resetPreviewViewport === 'function') {
        resetPreviewViewport();
      }

      // The retained Result's catalog owns this projection across mode inactivity.
      const features = state.featureCatalog?.value && generatedMode.value === mode.value
        ? featureStateFromCatalog(window.Vue.toRaw(state.featureCatalog.value), { mode: mode.value })
        : {};
      extractedFeatures.value = features.extractedFeatures || [];
      if (biologicalFeatures) biologicalFeatures.value = features.biologicalFeatures || [];
      featureRecordIds.value = features.featureRecordIds || [];
      selectedFeatureRecordIdx.value = 0;
      // Mode inactivity clears projections, not durable label/visibility intent.
      editableLabels.value = [];
      orthogroups.value = features.orthogroups || [];
      collinearGroups.value = features.collinearGroups || [];
      featureOrthogroupIndex.value = features.featureOrthogroupIndex || new Map();
      selectedOrthogroupAlignmentFeature.value = '';
      selectedOrthogroupId.value = '';
      orthogroupSearch.value = '';
      labelOverrideBuildWarning.value = '';
      labelSearch.value = '';
      featurePanelTab.value = 'colors';
      clickedPairwiseMatch.value = null;
      clickedLabel.value = null;
      labelTextScopeDialog.show = false;
      labelTextScopeDialog.labelKey = '';
      labelTextScopeDialog.newText = '';
      labelTextScopeDialog.sourceText = '';
      labelTextScopeDialog.featureId = '';
      labelTextScopeDialog.matchingCount = 0;
      hiddenLabelTextDialog.show = false;
      hiddenLabelTextDialog.featureId = '';
      hiddenLabelTextDialog.reason = '';
      resetRightDrawer();
      linearReorderNotice.value = '';
    }
  );

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
    try {
      const { colors, count } = parseColorTable(text);
      Object.entries(colors).forEach(([key, color]) => {
        currentColors.value[key] = color;
      });
      console.log(`Loaded ${count} colors from file.`);
    } catch (e) {
      errorLog.value = normalizeUserFacingError(e, { stage: 'request-validation' });
      return false;
    }
  });

  watchFileImport('t_color', async (text, isCurrent) => {
    try {
      const prepared = prepareSpecificColorImport(text, manualSpecificRules);
      if (!await featureActions.commitSpecificRules(prepared.nextRules, 'Import specific color rules', { isCurrent })) return false;
      console.log(`Loaded ${prepared.importedCount} rules from file.`);
    } catch (e) {
      if (isCurrent()) errorLog.value = normalizeUserFacingError(e, { operation: 'evaluateRules', stage: 'rule-validation' });
      return false;
    }
  });

  watchFileImport('qualifier_priority', (text) => {
    try {
      const { rules, count } = parsePriorityRules(text);
      rules.forEach((rule) => {
        const idx = manualPriorityRules.findIndex((r) => r.feat === rule.feat);
        if (idx >= 0) {
          manualPriorityRules[idx].order = rule.order;
        } else {
          manualPriorityRules.push({ feat: rule.feat, order: rule.order });
        }
      });
      console.log(`Loaded ${count} priority rules.`);
    } catch (e) {
      errorLog.value = normalizeUserFacingError(e, { stage: 'request-validation' });
      return false;
    }
  });

  watchFileImport('whitelist', (text) => {
    try {
      const { rules, count } = parseWhitelistRules(text);
      rules.forEach((rule) => manualWhitelist.push(rule));
      console.log(`Loaded ${count} whitelist rules.`);
    } catch (e) {
      errorLog.value = normalizeUserFacingError(e, { stage: 'request-validation' });
      return false;
    }
  });

  watchFileImport('blacklist', (text) => {
    try {
      const { words, count } = parseBlacklistWords(text);
      if (words.length > 0) {
        const existing = manualBlacklist.value ? manualBlacklist.value.trim() : '';
        const separator = existing && !existing.endsWith(',') ? ', ' : '';
        manualBlacklist.value = existing + separator + words.join(', ');
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
  return { waitForAuxiliaryFileImport, auxiliaryFileImportPending: () => pendingFileImports.size > 0, canRetryAuxiliaryImportFailure, retryAuxiliaryImportFailure };
};
