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
import {
  COMPOSITION_METADATA_ATTRIBUTE,
  COMPOSITION_SCHEMA_ATTRIBUTE
} from './legend-layout/composition-actions.js';
import { isCommittedSvgResultMounted } from '../services/svg-result-ingestion.js';

export const runRecordDiscoveryWatcher = async ({
  rollbackInProgress,
  semanticWatchersSuppressed,
  refresh
}) => {
  const suppress = Boolean(
    rollbackInProgress?.value || semanticWatchersSuppressed?.value
  );
  await refresh({ suppress });
  return !suppress;
};

export const setupWatchers = ({
  state,
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
    svgContainer,
    layoutPreferences,
    suppressCircularMultiRecordDefaults,
    featureRecordIds,
    selectedFeatureRecordIdx,
    featureVisibilityManualRules,
    featureVisibilityOverrides,
    featureVisibilitySelectorCache,
    featurePanelTab,
    labelSearch,
    orthogroups,
    collinearGroups,
    featureOrthogroupIndex,
    selectedOrthogroupAlignmentFeature,
    orthogroupNameOverrides,
    orthogroupDescriptionOverrides,
    selectedOrthogroupId,
    orthogroupSearch,
    labelOverrideContextKey,
    labelTextBulkOverrides,
    labelTextFeatureOverrides,
    canonicalLabelOverrideRows,
    labelTextFeatureOverrideSources,
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
    setupLegendDrag,
    refreshLegendDragAffordances,
    syncFileLegendEntries
  } = legendActions;

  const { applyPaletteToSvg, applySpecificRulesToSvg } = svgActions;
  const { attachSvgFeatureHandlers, refreshFeatureOverrides, syncLabelEditor } = featureActions;
  const {
    applyCanvasPadding,
    captureBaseConfig,
    captureOriginalStroke,
    repositionForLegendChange,
    refreshDiagramDragAffordances,
    setupDiagramDrag
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

  const hasLabelOverrides = () =>
    Object.keys(labelTextFeatureOverrides).length > 0 ||
    Object.keys(labelTextBulkOverrides).length > 0 ||
    Object.keys(labelVisibilityOverrides).length > 0;

  const shouldSyncLabelEditor = () =>
    isFeatureDrawerMounted.value ||
    Boolean(clickedFeature.value) ||
    labelTextScopeDialog.show ||
    globalLabelModeDialog.show ||
    hasLabelOverrides();

  const getNow = () => (globalThis.performance?.now ? performance.now() : Date.now());
  const formatDuration = (ms) => `${ms.toFixed(1)}ms`;
  const measureTiming = (entries, label, fn) => {
    const startedAt = getNow();
    const result = fn();
    entries.push({ label, ms: getNow() - startedAt });
    return result;
  };
  const logPostGbdrawTimings = (entries) => {
    if (!entries || entries.length === 0) return;
    console.groupCollapsed('post-gbdraw timing');
    entries.forEach(({ label, ms, details }) => {
      console.info(`${label}: ${formatDuration(ms)}${details ? ` (${details})` : ''}`);
    });
    console.groupEnd();
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
      featureVisibilitySelectorCache,
      featureVisibilityOverrides
    );
    replacePlainObject(
      featureVisibilitySelectorCache,
      nextCache
    );
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
      applyPaletteToSvg();
      applySpecificRulesToSvg();
      if (extractedFeatures.value.length > 0) {
        refreshFeatureOverrides(extractedFeatures.value);
      }

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

  watch(svgContent, () => {
    const isIncrementalEdit = Boolean(skipCaptureBaseConfig.value);
    skipCaptureBaseConfig.value = false;
    skipPositionReapply.value = false;

    nextTick(() => {
      const timingEntries = [];
      const svg = svgContainer.value?.querySelector('svg') || null;
      if (svg) {
        previewRuntime?.mountResultSvg?.(selectedResultIndex.value, svg);
      } else {
        previewRuntime?.clearActiveRuntime?.();
      }

      if (!skipExtractOnSvgChange.value) {
        measureTiming(timingEntries, 'watch(svgContent) extractLegendEntries', extractLegendEntries);
      }
      const shouldBindComposition = Boolean(
        svg && (
          !isIncrementalEdit ||
          svg.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) !== null ||
          svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE) !== null
        )
      );
      if (shouldBindComposition) {
        try {
          measureTiming(timingEntries, 'watch(svgContent) bind composition metadata', captureBaseConfig);
        } catch (error) {
          errorLog.value = {
            summary: error?.message || 'The SVG composition metadata is invalid.',
            details: []
          };
          console.error('Could not bind SVG composition metadata.', error);
          return;
        }
      }
      measureTiming(timingEntries, 'watch(svgContent) setupLegendDrag', setupLegendDrag);
      measureTiming(timingEntries, 'watch(svgContent) setupDiagramDrag', () => setupDiagramDrag(isIncrementalEdit));
      measureTiming(timingEntries, 'watch(svgContent) attachSvgFeatureHandlers', attachSvgFeatureHandlers);
      if (shouldSyncLabelEditor()) {
        measureTiming(timingEntries, 'watch(svgContent) syncLabelEditor', syncLabelEditor);
      }

      if (!isIncrementalEdit) {
        measureTiming(timingEntries, 'watch(svgContent) captureOriginalStroke', captureOriginalStroke);
        canvasPadding.top = 0;
        canvasPadding.right = 0;
        canvasPadding.bottom = 0;
        canvasPadding.left = 0;

      }

      logPostGbdrawTimings(timingEntries);
    });
  });

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
    if (semanticFileWatchersSuppressed.value) return;
    refreshFeatureVisibilitySelectorCache();
    if (!svgContent.value) return;
    nextTick(() => {
      if (semanticFileWatchersSuppressed.value) return;
      const timingEntries = [];
      measureTiming(timingEntries, 'watch(extractedFeatures) apply palette colors', applyPaletteToSvg);
      measureTiming(timingEntries, 'watch(extractedFeatures) apply specific rules', applySpecificRulesToSvg);
      measureTiming(timingEntries, 'watch(extractedFeatures) refresh delegated feature handlers', attachSvgFeatureHandlers);
      logPostGbdrawTimings(timingEntries);
    });
  });

  watch(featureSelectorSafetyScope, refreshFeatureVisibilitySelectorCache, { immediate: true });

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

      if (typeof resetPreviewViewport === 'function') {
        resetPreviewViewport();
      }

      extractedFeatures.value = [];
      if (biologicalFeatures) biologicalFeatures.value = [];
      featureSelectorSafetyScope.value = [];
      featureRecordIds.value = [];
      selectedFeatureRecordIdx.value = 0;
      featureVisibilityManualRules.splice(0);
      Object.keys(featureVisibilityOverrides).forEach((k) => delete featureVisibilityOverrides[k]);
      Object.keys(featureVisibilitySelectorCache).forEach((k) => delete featureVisibilitySelectorCache[k]);
      editableLabels.value = [];
      Object.keys(labelTextFeatureOverrides).forEach((k) => delete labelTextFeatureOverrides[k]);
      Object.keys(labelTextBulkOverrides).forEach((k) => delete labelTextBulkOverrides[k]);
      Object.keys(labelTextFeatureOverrideSources).forEach((k) => delete labelTextFeatureOverrideSources[k]);
      Object.keys(labelVisibilityOverrides).forEach((k) => delete labelVisibilityOverrides[k]);
      orthogroups.value = [];
      collinearGroups.value = [];
      featureOrthogroupIndex.value = new Map();
      selectedOrthogroupAlignmentFeature.value = '';
      selectedOrthogroupId.value = '';
      orthogroupSearch.value = '';
      Object.keys(orthogroupNameOverrides).forEach((k) => delete orthogroupNameOverrides[k]);
      Object.keys(orthogroupDescriptionOverrides).forEach((k) => delete orthogroupDescriptionOverrides[k]);
      labelOverrideContextKey.value = '';
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

  watch(() => files.d_color, async (newFile) => {
    if (semanticFileWatchersSuppressed.value) return;
    if (!newFile) return;
    try {
      const text = await readFileText(newFile);
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

  watch(() => files.t_color, async (newFile) => {
    if (semanticFileWatchersSuppressed.value) return;
    if (!newFile) return;
    try {
      const text = await readFileText(newFile);
      const prepared = prepareSpecificColorImport(text, manualSpecificRules);
      const previousCaptions = Array.from(fileLegendCaptions.value);
      const previousFileIntents = buildLegendIntents(
        manualSpecificRules.filter((rule) => rule.fromFile),
        { conflictPolicy: 'last-wins' }
      ).intents;
      if (svgContent.value) {
        await nextTick();
        await syncFileLegendEntries(prepared.intents, { previousFileIntents });
      }

      manualSpecificRules.splice(0, manualSpecificRules.length, ...prepared.nextRules);
      previousCaptions.forEach((caption) => addedLegendCaptions.value.delete(caption));
      fileLegendCaptions.value = new Set(prepared.fileLegendCaptions);
      prepared.fileLegendCaptions.forEach((caption) => addedLegendCaptions.value.add(caption));
      extractLegendEntries();
      console.log(`Loaded ${prepared.importedCount} rules from file.`);
    } catch (e) {
      console.error('Failed to load rules file:', e);
      alert(`Failed to load rules file. ${e?.message || 'Please check the TSV format.'}`);
    }
  });

  watch(() => files.qualifier_priority, async (newFile) => {
    if (semanticFileWatchersSuppressed.value) return;
    if (!newFile) return;
    try {
      const text = await readFileText(newFile);
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

  watch(() => files.whitelist, async (newFile) => {
    if (semanticFileWatchersSuppressed.value) return;
    if (!newFile) return;
    try {
      const text = await readFileText(newFile);
      const { rules, count } = parseWhitelistRules(text);
      rules.forEach((rule) => manualWhitelist.push(rule));
      console.log(`Loaded ${count} whitelist rules.`);
    } catch (e) {
      console.error('Failed to load whitelist file:', e);
      alert('Failed to load whitelist file.');
    }
  });

  watch(() => files.blacklist, async (newFile) => {
    if (semanticFileWatchersSuppressed.value) return;
    if (!newFile) return;
    try {
      const text = await readFileText(newFile);
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
        refresh: ({ suppress }) => refreshCircularRecordOrder({ suppress })
      });
    }
  );
  watch(
    () => [
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
};
