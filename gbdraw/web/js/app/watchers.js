import {
  parseBlacklistWords,
  parseColorTable,
  parsePriorityRules,
  parseWhitelistRules
} from './file-imports.js';
import {
  buildLegendIntents,
  prepareSpecificColorImport
} from './specific-color-rules.js';
import {
  buildFeatureVisibilitySelectorCache,
  preserveFeatureVisibilitySelectorCacheForOverrides
} from './feature-visibility.js';
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
  watch,
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
    featureSelectorSafetyScope,
    addedLegendCaptions,
    layoutRepositionMode,
    editableLabels,
    results,
    svgContent,
    selectedResultIndex,
    form,
    generatedLegendPosition,
    generatedMode,
    shouldDeferCircularPreviewUpdates,
    mode,
    cInputType,
    lInputType,
    canvasPadding,
    skipCaptureBaseConfig,
    skipPositionReapply,
    skipExtractOnSvgChange,
    trustedArtifactRestoreInProgress,
    svgContainer,
    layoutPreferences,
    suppressCircularMultiRecordDefaults,
    featureRecordIds,
    selectedFeatureRecordIdx,
    featureVisibilityOverrides,
    replaceFeatureVisibilitySelectorCacheOwner,
    featurePanelTab,
    labelSearch,
    orthogroups,
    collinearGroups,
    featureOrthogroupIndex,
    selectedOrthogroupAlignmentFeature,
    selectedOrthogroupId,
    orthogroupSearch,
    labelTextBulkOverrides,
    labelTextFeatureOverrides,
    canonicalLabelOverrideRows,
    labelVisibilityOverrides,
    labelOverrideBuildWarning,
    isFeatureDrawerMounted,
    clickedFeature,
    clickedPairwiseMatch,
    clickedLabel,
    labelTextScopeDialog,
    globalLabelModeDialog,
    files,
    currentColors,
    paletteInstantPreviewEnabled,
    pendingPaletteName,
    fileLegendCaptions,
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
    labelReflowRequestReason,
    labelReflowForceRequestSeq,
    labelReflowForceRequestReason,
    labelLayoutDirtyReason,
    errorLog
  } = state;

  const {
    removeLegendEntry,
    addLegendEntry,
    extractLegendEntries,
    refreshLegendDragAffordances,
    syncFileLegendEntries
  } = legendActions;

  const { applyPaletteToSvg, applySpecificRulesToSvg } = svgActions;
  const { refreshFeatureOverrides, syncLabelEditor } = featureActions;
  const {
    applyCanvasPadding,
    repositionForLegendChange,
    refreshDiagramDragAffordances
  } = legendLayout;
  const {
    applyPaletteDraftToPreview,
    scheduleDefinitionUpdate,
    cancelDefinitionUpdate,
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

  const replacePlainObject = (target, source = {}) => {
    Object.keys(target || {}).forEach((key) => delete target[key]);
    Object.entries(source || {}).forEach(([key, value]) => {
      target[key] = value;
    });
  };

  const refreshFeatureVisibilitySelectorCache = () => {
    const nextCache = preserveFeatureVisibilitySelectorCacheForOverrides(
      buildFeatureVisibilitySelectorCache(extractedFeatures.value, featureSelectorSafetyScope.value),
      state.featureVisibilitySelectorCache,
      featureVisibilityOverrides
    );
    if (typeof replaceFeatureVisibilitySelectorCacheOwner === 'function') {
      replaceFeatureVisibilitySelectorCacheOwner(nextCache);
    } else {
      replacePlainObject(state.featureVisibilitySelectorCache, nextCache);
    }
  };

  const scheduleCircularDefinitionUpdate = () => {
    if (mode.value !== 'circular') return;
    if (generatedMode.value !== mode.value) return;
    if (shouldDeferCircularPreviewUpdates.value) {
      cancelDefinitionUpdate();
      return;
    }
    scheduleDefinitionUpdate();
  };

  watch(
    () => [...manualSpecificRules],
    async (newRules, oldRules) => {
      if (semanticFileWatchersSuppressed.value) return;
      if (extractedFeatures.value.length > 0) {
        refreshFeatureOverrides(extractedFeatures.value);
      }
      applyPaletteToSvg();
      applySpecificRulesToSvg();

      const currentCaptions = new Set(newRules.filter((r) => r.cap).map((r) => r.cap));
      const oldCaptions = new Set((oldRules || []).filter((r) => r.cap).map((r) => r.cap));

      const removedFromRules = [...oldCaptions].filter((cap) => !currentCaptions.has(cap));
      const removedFromTracked = [...addedLegendCaptions.value].filter((cap) => !currentCaptions.has(cap));

      const allRemovedCaptions = new Set([...removedFromRules, ...removedFromTracked]);

      for (const cap of allRemovedCaptions) {
        removeLegendEntry(cap);
        addedLegendCaptions.value.delete(cap);
      }
    },
    { deep: true }
  );

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
    () => form.legend,
    (newPos, oldPos) => {
      if (generatedMode.value !== mode.value) return;
      if (mode.value === 'circular' && shouldDeferCircularPreviewUpdates.value) return;
      if (
        svgContent.value &&
        oldPos !== undefined &&
        newPos !== oldPos &&
        newPos !== generatedLegendPosition.value
      ) {
        nextTick(() => {
          if (mode.value === 'circular' && shouldDeferCircularPreviewUpdates.value) return;
          repositionForLegendChange(newPos, generatedLegendPosition.value);
        });
      }
    }
  );

  watch(
    () => form.multi_record_canvas,
    (enabled, previousEnabled) => {
      cancelDefinitionUpdate();
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
    skipPositionReapply.value = false;

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
        errorLog.value = {
          summary: error?.message || 'The mounted preview could not be prepared.',
          details: []
        };
        console.error('Could not bind the mounted SVG Result.', error);
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
      skipPositionReapply.value = false;
    }
  );

  watch(extractedFeatures, () => {
    if (semanticFileWatchersSuppressed.value || trustedArtifactRestoreInProgress.value) return;
    refreshFeatureVisibilitySelectorCache();
  });

  watch(featureSelectorSafetyScope, () => {
    if (trustedArtifactRestoreInProgress.value) return;
    refreshFeatureVisibilitySelectorCache();
  }, { immediate: true });

  watch(
    [labelTextFeatureOverrides, labelTextBulkOverrides, labelVisibilityOverrides],
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
      await runLabelReflow(labelReflowRequestReason.value || 'label-edit');
    }
  );

  watch(
    () => labelReflowForceRequestSeq.value,
    async (nextSeq, prevSeq) => {
      if (nextSeq === prevSeq) return;
      if (mode.value === 'circular' && shouldDeferCircularPreviewUpdates.value) return;
      if (typeof runLabelReflow !== 'function') return;
      await runLabelReflow(labelReflowForceRequestReason.value || 'label-edit');
    }
  );

  watch(
    () => mode.value,
    () => {
      if (semanticFileWatchersSuppressed.value) return;
      cancelDefinitionUpdate();

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
      featureSelectorSafetyScope.value = features.featureSelectorSafetyScope || [];
      featureRecordIds.value = features.featureRecordIds || [];
      selectedFeatureRecordIdx.value = 0;
      // Mode inactivity clears projections, not durable label/visibility intent.
      // The editor's existing identity reconciliation handles source replacement.
      if (typeof replaceFeatureVisibilitySelectorCacheOwner === 'function') {
        replaceFeatureVisibilitySelectorCacheOwner({});
      } else {
        replacePlainObject(state.featureVisibilitySelectorCache, {});
      }
      editableLabels.value = [];
      orthogroups.value = features.orthogroups || [];
      collinearGroups.value = features.collinearGroups || [];
      featureOrthogroupIndex.value = features.featureOrthogroupIndex || new Map();
      selectedOrthogroupAlignmentFeature.value = '';
      selectedOrthogroupId.value = '';
      orthogroupSearch.value = '';
      labelOverrideBuildWarning.value = '';
      labelLayoutDirtyReason.value = '';
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
      globalLabelModeDialog.show = false;
      globalLabelModeDialog.featureId = '';
      globalLabelModeDialog.featureType = '';
      globalLabelModeDialog.resolve = null;
      resetRightDrawer();
      linearReorderNotice.value = '';
    }
  );

  const pendingFileImports = new WeakMap();
  let fileImportApplications = Promise.resolve();
  const restoredFileSelections = new Map();
  const watchFileImport = (key, apply) => watch(() => files[key], (file, previousFile) => {
    const restored = restoredFileSelections.has(key) && restoredFileSelections.get(key) === file;
    restoredFileSelections.delete(key);
    if (restored) return;
    if (semanticFileWatchersSuppressed.value || !file) return;
    const ruleContext = key === 't_color' ? rulePreparation.snapshot() : null;
    const isCurrent = () => files[key] === file && !semanticFileWatchersSuppressed.value
      && (!ruleContext || rulePreparation.isCurrent(ruleContext));
    const pending = readFileText(file).then((text) => {
      // Reads may finish out of order; serialize only their live application.
      const application = fileImportApplications.then(async () => {
        if (!isCurrent()) return;
        if (await apply(text, isCurrent) === false && isCurrent()) {
          restoredFileSelections.set(key, previousFile);
          files[key] = previousFile;
        }
      });
      fileImportApplications = application.catch(() => {});
      return application;
    }).catch((error) => {
      if (isCurrent()) alert(`Failed to read ${file.name || 'uploaded table'}: ${error.message}`);
    });
    pendingFileImports.set(file, pending);
  });
  const waitForAuxiliaryFileImport = (file) => pendingFileImports.get(file);

  watchFileImport('d_color', (text) => {
    try {
      const { colors, count } = parseColorTable(text);
      Object.entries(colors).forEach(([key, color]) => {
        currentColors.value[key] = color;
      });
      console.log(`Loaded ${count} colors from file.`);
    } catch (e) {
      console.error('Failed to load color file:', e);
      alert('Failed to load color file. Please check the TSV format.');
    }
  });

  watchFileImport('t_color', async (text, isCurrent) => {
    try {
      const prepared = prepareSpecificColorImport(text, manualSpecificRules);
      if (!await rulePreparation.prepare(prepared.nextRules) || !isCurrent()) return;
      const previousCaptions = Array.from(fileLegendCaptions.value);
      const previousFileIntents = buildLegendIntents(
        manualSpecificRules.filter((rule) => rule.fromFile),
        { conflictPolicy: 'last-wins' }
      ).intents;
      if (svgContent.value) {
        await nextTick();
        await syncFileLegendEntries(prepared.intents, { previousFileIntents });
      }
      if (!isCurrent()) return;

      manualSpecificRules.splice(0, manualSpecificRules.length, ...prepared.nextRules);
      previousCaptions.forEach((caption) => addedLegendCaptions.value.delete(caption));
      fileLegendCaptions.value = new Set(prepared.fileLegendCaptions);
      prepared.fileLegendCaptions.forEach((caption) => addedLegendCaptions.value.add(caption));
      extractLegendEntries();
      console.log(`Loaded ${prepared.importedCount} rules from file.`);
    } catch (e) {
      console.error('Failed to load rules file:', e);
      if (isCurrent()) alert(`Failed to load rules file. ${e?.message || 'Please check the TSV format.'}`);
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
      console.error('Failed to load priority file:', e);
      alert('Failed to load priority file.');
    }
  });

  watchFileImport('whitelist', (text) => {
    try {
      const { rules, count } = parseWhitelistRules(text);
      rules.forEach((rule) => manualWhitelist.push(rule));
      console.log(`Loaded ${count} whitelist rules.`);
    } catch (e) {
      console.error('Failed to load whitelist file:', e);
      alert('Failed to load whitelist file.');
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
      console.error('Failed to load blacklist file:', e);
      alert('Failed to load blacklist file.');
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

  watch(() => form.species, scheduleCircularDefinitionUpdate);
  watch(() => form.strain, scheduleCircularDefinitionUpdate);
  watch(() => form.plot_title, scheduleCircularDefinitionUpdate);
  watch(() => state.adv.def_font_size, scheduleCircularDefinitionUpdate);
  watch(
    () => state.adv.plot_title_position,
    () => {
      scheduleCircularDefinitionUpdate();
    }
  );
  watch(() => state.adv.plot_title_font_size, scheduleCircularDefinitionUpdate);
  watch(() => state.adv.keep_full_definition_with_plot_title, scheduleCircularDefinitionUpdate);
  watch(
    () => [
      semanticFileWatchersSuppressed.value,
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
        refresh: ({ suppress }) => refreshCircularRecordOrder({ suppress })
      });
    }
  );
  watch(
    () => [
      semanticFileWatchersSuppressed.value,
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
      console.warn('Could not load browser palette definitions.', error);
    }
  });
  return { waitForAuxiliaryFileImport };
};
