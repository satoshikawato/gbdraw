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
  normalizeCircularPlotTitlePosition,
  normalizeLinearPlotTitlePosition
} from './plot-title-position.js';

export const setupWatchers = ({
  state,
  watch,
  nextTick,
  onMounted,
  debugLog,
  pyodideManager,
  legendActions,
  svgActions,
  featureActions,
  legendLayout,
  resultsManager,
  runLabelReflow,
  refreshCircularRecordOrder,
  refreshLinearRecordSelectors,
  resetPreviewViewport,
  previewRuntime = null,
  prepareDiagramGenerationWorker
}) => {
  const {
    manualSpecificRules,
    extractedFeatures,
    featureSelectorSafetyScope,
    addedLegendCaptions,
    layoutRepositionMode,
    editableLabels,
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
    diagramElementBaseTransforms,
    svgContainer,
    diagramElements,
    linearBaseConfig,
    circularLegendPosition,
    linearLegendPosition,
    circularPlotTitlePosition,
    linearPlotTitlePosition,
    circularSingleRecordLegendPosition,
    circularSingleRecordPlotTitlePosition,
    circularMultiRecordLegendPosition,
    circularMultiRecordPlotTitlePosition,
    suppressCircularMultiRecordDefaults,
    featureRecordIds,
    selectedFeatureRecordIdx,
    featureColorOverrides,
    featureVisibilityManualRules,
    featureVisibilityOverrides,
    featureVisibilitySelectorCache,
    featureStrokeOverrides,
    featurePanelTab,
    labelSearch,
    orthogroups,
    featureOrthogroupIndex,
    selectedOrthogroupAlignmentFeature,
    orthogroupNameOverrides,
    orthogroupDescriptionOverrides,
    selectedOrthogroupId,
    orthogroupSearch,
    showRightDrawer,
    rightDrawerTab,
    labelOverrideContextKey,
    labelTextBulkOverrides,
    labelTextFeatureOverrides,
    canonicalLabelOverrideRows,
    labelTextFeatureOverrideSources,
    labelVisibilityOverrides,
    labelOverrideBuildWarning,
    showFeaturePanel,
    clickedFeature,
    clickedPairwiseMatch,
    clickedLabel,
    labelTextScopeDialog,
    globalLabelModeDialog,
    files,
    currentColors,
    paletteInstantPreviewEnabled,
    pendingPaletteName,
    pyodideReady,
    fileLegendCaptions,
    semanticFileWatchersSuppressed,
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
    labelLayoutDirtyReason
  } = state;

  const {
    removeLegendEntry,
    addLegendEntry,
    extractLegendEntries,
    setupLegendDrag,
    refreshLegendDragAffordances,
    reapplyStrokeOverrides,
    syncFileLegendEntries
  } = legendActions;

  const {
    applyPaletteToSvg,
    applySpecificRulesToSvg,
    ensureUniquePairwiseGradientIds,
    ensureUniqueSkewClipPathIds
  } = svgActions;
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

  const getStoredCircularLayout = (useMultiRecord) => {
    if (useMultiRecord) {
      const fallbackLegend = normalizeLegendPosition(circularSingleRecordLegendPosition.value, 'left');
      const fallbackPlotTitlePosition = normalizeCircularPlotTitlePosition(circularSingleRecordPlotTitlePosition.value);
      return {
        legend: hasStoredLayoutValue(circularMultiRecordLegendPosition.value)
          ? normalizeLegendPosition(circularMultiRecordLegendPosition.value, fallbackLegend)
          : fallbackLegend,
        plotTitlePosition: hasStoredLayoutValue(circularMultiRecordPlotTitlePosition.value)
          ? normalizeCircularPlotTitlePosition(circularMultiRecordPlotTitlePosition.value)
          : fallbackPlotTitlePosition
      };
    }

    return {
      legend: normalizeLegendPosition(circularSingleRecordLegendPosition.value, 'left'),
      plotTitlePosition: normalizeCircularPlotTitlePosition(circularSingleRecordPlotTitlePosition.value)
    };
  };

  const syncCurrentCircularLayoutCache = () => {
    const normalizedLegend = normalizeLegendPosition(form.legend, 'left');
    const normalizedPlotTitlePosition = normalizeCircularPlotTitlePosition(state.adv.plot_title_position);

    circularLegendPosition.value = normalizedLegend;
    circularPlotTitlePosition.value = normalizedPlotTitlePosition;

    if (form.multi_record_canvas) {
      circularMultiRecordLegendPosition.value = normalizedLegend;
      circularMultiRecordPlotTitlePosition.value = normalizedPlotTitlePosition;
    } else {
      circularSingleRecordLegendPosition.value = normalizedLegend;
      circularSingleRecordPlotTitlePosition.value = normalizedPlotTitlePosition;
    }

    return {
      legend: normalizedLegend,
      plotTitlePosition: normalizedPlotTitlePosition
    };
  };

  const restoreCircularLayoutCache = (useMultiRecord) => {
    const nextLayout = getStoredCircularLayout(useMultiRecord);
    circularLegendPosition.value = nextLayout.legend;
    circularPlotTitlePosition.value = nextLayout.plotTitlePosition;

    if (form.legend !== nextLayout.legend) {
      form.legend = nextLayout.legend;
    }
    if (state.adv.plot_title_position !== nextLayout.plotTitlePosition) {
      state.adv.plot_title_position = nextLayout.plotTitlePosition;
    }

    return nextLayout;
  };

  const hasStoredCircularMultiRecordLayout = () =>
    hasStoredLayoutValue(circularMultiRecordLegendPosition.value) ||
    hasStoredLayoutValue(circularMultiRecordPlotTitlePosition.value);

  const applyCircularMultiRecordSmartDefaults = () => {
    const singleLayout = getStoredCircularLayout(false);
    circularMultiRecordLegendPosition.value = singleLayout.legend === 'left' ? 'bottom' : singleLayout.legend;
    circularMultiRecordPlotTitlePosition.value =
      singleLayout.plotTitlePosition === 'none' ? 'bottom' : singleLayout.plotTitlePosition;
  };

  const hasLabelOverrides = () =>
    Object.keys(labelTextFeatureOverrides).length > 0 ||
    Object.keys(labelTextBulkOverrides).length > 0 ||
    Object.keys(labelVisibilityOverrides).length > 0;

  const shouldSyncLabelEditor = () =>
    showFeaturePanel.value ||
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
      if (mode.value === 'circular') {
        const normalizedLegend = normalizeLegendPosition(newPos, 'left');
        circularLegendPosition.value = normalizedLegend;
        if (form.multi_record_canvas) {
          circularMultiRecordLegendPosition.value = normalizedLegend;
        } else {
          circularSingleRecordLegendPosition.value = normalizedLegend;
        }
      } else if (mode.value === 'linear') {
        linearLegendPosition.value = normalizeLegendPosition(newPos, 'bottom');
      }

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
          syncCurrentCircularLayoutCache();
        } else {
          applyCircularMultiRecordSmartDefaults();
        }
      }

      if (suppressCircularMultiRecordDefaults.value) {
        suppressCircularMultiRecordDefaults.value = false;
      }

      restoreCircularLayoutCache(Boolean(enabled));
    }
  );

  watch(svgContent, () => {
    const isIncrementalEdit = Boolean(skipCaptureBaseConfig.value);
    const shouldSkipPositionReapply = Boolean(skipPositionReapply.value);
    skipCaptureBaseConfig.value = false;
    skipPositionReapply.value = false;

    let savedBaseTransformsById = null;
    if (isIncrementalEdit && diagramElementBaseTransforms.value.size > 0) {
      savedBaseTransformsById = new Map();
      for (const [el, transform] of diagramElementBaseTransforms.value) {
        const id = el.id || '';
        if (id) {
          if (!savedBaseTransformsById.has(id)) {
            savedBaseTransformsById.set(id, []);
          }
          savedBaseTransformsById.get(id).push(transform);
        }
      }
      debugLog('Saved', savedBaseTransformsById.size, 'base transform IDs before DOM update');
    }

    nextTick(() => {
      const timingEntries = [];
      const svg = svgContainer.value?.querySelector('svg') || null;
      if (svg) {
        previewRuntime?.mountResultSvg?.(selectedResultIndex.value, svg);
      } else {
        previewRuntime?.clearActiveRuntime?.();
      }

      if (svgContainer.value) {
        if (svg) {
          const tickEl = svg.getElementById('tick');
          if (tickEl) {
            debugLog(`After DOM update - tick transform: ${tickEl.getAttribute('transform')}`);
          }
        }
      }

      if (svg && !isIncrementalEdit) {
        const normalizedSvgChanged = measureTiming(timingEntries, 'watch(svgContent) normalize unique SVG ids', () => {
          const skewChanged = ensureUniqueSkewClipPathIds(svg);
          const gradientChanged = ensureUniquePairwiseGradientIds(svg);
          return Boolean(skewChanged || gradientChanged);
        });

        if (normalizedSvgChanged) {
          previewRuntime?.markActiveResultDirty?.('svg-normalization');
          timingEntries.push({
            label: 'watch(svgContent) persist normalized SVG',
            ms: 0,
            details: 'live DOM only'
          });
        }
      }

      if (!skipExtractOnSvgChange.value) {
        measureTiming(timingEntries, 'watch(svgContent) extractLegendEntries', extractLegendEntries);
      }
      measureTiming(timingEntries, 'watch(svgContent) setupLegendDrag', setupLegendDrag);
      measureTiming(timingEntries, 'watch(svgContent) setupDiagramDrag', () => setupDiagramDrag(isIncrementalEdit));
      measureTiming(timingEntries, 'watch(svgContent) attachSvgFeatureHandlers', attachSvgFeatureHandlers);
      if (shouldSyncLabelEditor()) {
        measureTiming(timingEntries, 'watch(svgContent) syncLabelEditor', syncLabelEditor);
      }

      if (!isIncrementalEdit) {
        measureTiming(timingEntries, 'watch(svgContent) captureBaseConfig', captureBaseConfig);
        measureTiming(timingEntries, 'watch(svgContent) captureOriginalStroke', captureOriginalStroke);
        debugLog('Full base config capture (fresh generation)');

        canvasPadding.top = 0;
        canvasPadding.right = 0;
        canvasPadding.bottom = 0;
        canvasPadding.left = 0;

        measureTiming(timingEntries, 'watch(svgContent) reapplyStrokeOverrides', reapplyStrokeOverrides);
      } else if (savedBaseTransformsById && savedBaseTransformsById.size > 0) {
        debugLog('Incremental edit - remapping base transforms');

        const newBaseTransforms = new Map();
        const idCounters = new Map();

        diagramElements.value.forEach((newEl) => {
          const id = newEl.id || '';
          const transforms = savedBaseTransformsById.get(id);
          if (transforms && transforms.length > 0) {
            const idx = idCounters.get(id) || 0;
            if (idx < transforms.length) {
              newBaseTransforms.set(newEl, transforms[idx]);
              idCounters.set(id, idx + 1);
            }
          }
        });

        diagramElementBaseTransforms.value = newBaseTransforms;
        debugLog('Remapped', newBaseTransforms.size, 'base transforms to new elements');

        if (!shouldSkipPositionReapply) {
          const isLinear = mode.value === 'linear';
          const currentLegendPos = form.legend;
          const currentDiagramPos = generatedLegendPosition.value;

          debugLog('Checking position after incremental edit:', {
            currentLegendPos,
            currentDiagramPos,
            isLinear,
            shouldSkipPositionReapply
          });

          if (isLinear && currentLegendPos && currentDiagramPos && currentLegendPos !== currentDiagramPos) {
            debugLog('Position differs from current - reapplying position shift');
            skipCaptureBaseConfig.value = true;
            const originalGeneratedPos = linearBaseConfig.value.generatedPosition;
            repositionForLegendChange(currentLegendPos, originalGeneratedPos);
            return;
          }
        } else {
          debugLog('Skipping position reapply (triggered by repositionForLegendChange)');
        }
      } else {
        debugLog('Incremental edit but no base transforms to remap');
      }

      logPostGbdrawTimings(timingEntries);
    });
  });

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
    (newMode, oldMode) => {
      if (semanticFileWatchersSuppressed.value) return;
      cancelDefinitionUpdate();

      if (oldMode === 'circular') {
        syncCurrentCircularLayoutCache();
      } else if (oldMode === 'linear') {
        linearLegendPosition.value = normalizeLegendPosition(form.legend, 'bottom');
        linearPlotTitlePosition.value = normalizeLinearPlotTitlePosition(state.adv.plot_title_position);
      }

      if (newMode === 'circular') {
        restoreCircularLayoutCache(Boolean(form.multi_record_canvas));
      } else if (newMode === 'linear') {
        form.legend = linearLegendPosition.value;
        state.adv.plot_title_position = linearPlotTitlePosition.value;
      }

      if (typeof resetPreviewViewport === 'function') {
        resetPreviewViewport();
      }

      extractedFeatures.value = [];
      featureSelectorSafetyScope.value = [];
      featureRecordIds.value = [];
      selectedFeatureRecordIdx.value = 0;
      Object.keys(featureColorOverrides).forEach((k) => delete featureColorOverrides[k]);
      featureVisibilityManualRules.splice(0);
      Object.keys(featureVisibilityOverrides).forEach((k) => delete featureVisibilityOverrides[k]);
      Object.keys(featureVisibilitySelectorCache).forEach((k) => delete featureVisibilitySelectorCache[k]);
      Object.keys(featureStrokeOverrides).forEach((k) => delete featureStrokeOverrides[k]);
      editableLabels.value = [];
      Object.keys(labelTextFeatureOverrides).forEach((k) => delete labelTextFeatureOverrides[k]);
      Object.keys(labelTextBulkOverrides).forEach((k) => delete labelTextBulkOverrides[k]);
      Object.keys(labelTextFeatureOverrideSources).forEach((k) => delete labelTextFeatureOverrideSources[k]);
      Object.keys(labelVisibilityOverrides).forEach((k) => delete labelVisibilityOverrides[k]);
      orthogroups.value = [];
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
      showFeaturePanel.value = false;
      showRightDrawer.value = false;
      rightDrawerTab.value = 'features';
      linearReorderNotice.value = '';
    }
  );

  watch(() => files.d_color, async (newFile) => {
    if (semanticFileWatchersSuppressed.value) return;
    if (!newFile) return;
    try {
      const text = await newFile.text();
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
      const text = await newFile.text();
      const prepared = prepareSpecificColorImport(text, manualSpecificRules);
      const previousCaptions = Array.from(fileLegendCaptions.value);
      const previousFileIntents = buildLegendIntents(
        manualSpecificRules.filter((rule) => rule.fromFile),
        { conflictPolicy: 'last-wins' }
      ).intents;
      if (pyodideReady.value && svgContent.value) {
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
      const text = await newFile.text();
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
      const text = await newFile.text();
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
      const text = await newFile.text();
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
    () => showFeaturePanel.value,
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
    (newPos) => {
      if (mode.value === 'circular') {
        const normalizedPlotTitlePosition = normalizeCircularPlotTitlePosition(newPos);
        circularPlotTitlePosition.value = normalizedPlotTitlePosition;
        if (form.multi_record_canvas) {
          circularMultiRecordPlotTitlePosition.value = normalizedPlotTitlePosition;
        } else {
          circularSingleRecordPlotTitlePosition.value = normalizedPlotTitlePosition;
        }
      } else if (mode.value === 'linear') {
        linearPlotTitlePosition.value = normalizeLinearPlotTitlePosition(newPos);
      }

      scheduleCircularDefinitionUpdate();
    }
  );
  watch(() => state.adv.plot_title_font_size, scheduleCircularDefinitionUpdate);
  watch(() => state.adv.keep_full_definition_with_plot_title, scheduleCircularDefinitionUpdate);
  watch(
    () => [mode.value, cInputType.value, files.c_gb, files.c_gff, files.c_fasta, pyodideReady.value],
    async () => {
      if (typeof refreshCircularRecordOrder !== 'function') return;
      await refreshCircularRecordOrder();
    }
  );
  watch(
    () => [
      mode.value,
      lInputType.value,
      pyodideReady.value,
      ...linearSeqs.flatMap((seq) => [
        seq.uid,
        lInputType.value === 'gff' ? seq.gff : seq.gb,
        lInputType.value === 'gff' ? seq.fasta : null
      ])
    ],
    async () => {
      if (typeof refreshLinearRecordSelectors !== 'function') return;
      await refreshLinearRecordSelectors();
    },
    { immediate: true }
  );

  onMounted(async () => {
    const diagramWorkerPromise = typeof prepareDiagramGenerationWorker === 'function'
      ? (async () => {
          await nextTick();
          await prepareDiagramGenerationWorker();
        })()
      : Promise.resolve();

    await pyodideManager.initPyodide();
    await diagramWorkerPromise;
  });
};
