import {
  parseBlacklistWords,
  parseColorTable,
  parsePriorityRules,
  parseSpecificRules,
  parseWhitelistRules
} from './file-imports.js';

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
  resultsManager
}) => {
  const {
    manualSpecificRules,
    extractedFeatures,
    addedLegendCaptions,
    svgContent,
    form,
    generatedLegendPosition,
    mode,
    canvasPadding,
    skipCaptureBaseConfig,
    skipPositionReapply,
    skipExtractOnSvgChange,
    diagramElementBaseTransforms,
    svgContainer,
    results,
    selectedResultIndex,
    diagramElements,
    linearBaseConfig,
    circularBaseConfig,
    circularLegendPosition,
    linearLegendPosition,
    featureRecordIds,
    selectedFeatureRecordIdx,
    featureColorOverrides,
    showFeaturePanel,
    files,
    currentColors,
    pyodideReady,
    fileLegendCaptions,
    manualPriorityRules,
    manualWhitelist,
    manualBlacklist,
    linearSeqs
  } = state;

  const {
    removeLegendEntry,
    addLegendEntry,
    extractLegendEntries,
    setupLegendDrag,
    reapplyStrokeOverrides
  } = legendActions;

  const { applySpecificRulesToSvg, ensureUniquePairwiseGradientIds } = svgActions;
  const { attachSvgFeatureHandlers, refreshFeatureOverrides } = featureActions;
  const {
    applyCanvasPadding,
    captureBaseConfig,
    captureOriginalStroke,
    repositionForLegendChange,
    setupDiagramDrag
  } = legendLayout;
  const { scheduleDefinitionUpdate } = resultsManager;

  watch(
    () => [...manualSpecificRules],
    async (newRules, oldRules) => {
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
    canvasPadding,
    () => {
      applyCanvasPadding();
    },
    { deep: true }
  );

  watch(
    () => form.legend,
    (newPos, oldPos) => {
      if (
        svgContent.value &&
        oldPos !== undefined &&
        newPos !== oldPos &&
        newPos !== generatedLegendPosition.value
      ) {
        nextTick(() => {
          repositionForLegendChange(newPos, generatedLegendPosition.value);
        });
      }
    }
  );

  watch(svgContent, () => {
    const isIncrementalEdit = skipCaptureBaseConfig.value;
    const shouldSkipPositionReapply = skipPositionReapply.value;
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
      if (svgContainer.value) {
        const svg = svgContainer.value.querySelector('svg');
        if (svg) {
          const tickEl = svg.getElementById('tick');
          if (tickEl) {
            console.log(`[DEBUG] After DOM update - tick transform: ${tickEl.getAttribute('transform')}`);
          }
        }
      }

      if (!skipExtractOnSvgChange.value) {
        extractLegendEntries();
      }
      setupLegendDrag();
      setupDiagramDrag(isIncrementalEdit);
      attachSvgFeatureHandlers();

      if (!isIncrementalEdit) {
        captureBaseConfig();
        captureOriginalStroke();
        debugLog('Full base config capture (fresh generation)');

        canvasPadding.top = 0;
        canvasPadding.right = 0;
        canvasPadding.bottom = 0;
        canvasPadding.left = 0;

        reapplyStrokeOverrides();

        if (svgContainer.value) {
          const svgEl = svgContainer.value.querySelector('svg');
          if (svgEl) {
            ensureUniquePairwiseGradientIds(svgEl);
            skipCaptureBaseConfig.value = true;
            const idx = selectedResultIndex.value;
            if (idx >= 0 && results.value.length > idx) {
              const serializer = new XMLSerializer();
              results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svgEl) };
            }
          }
        }
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

          if (currentLegendPos && currentDiagramPos && currentLegendPos !== currentDiagramPos) {
            debugLog('Position differs from current - reapplying position shift');
            skipCaptureBaseConfig.value = true;
            const originalGeneratedPos = isLinear
              ? linearBaseConfig.value.generatedPosition
              : circularBaseConfig.value.generatedPosition;
            repositionForLegendChange(currentLegendPos, originalGeneratedPos);
            return;
          }
        } else {
          debugLog('Skipping position reapply (triggered by repositionForLegendChange)');
        }
      } else {
        debugLog('Incremental edit but no base transforms to remap');
      }
    });
  });

  watch(
    () => mode.value,
    (newMode, oldMode) => {
      if (oldMode === 'circular') {
        circularLegendPosition.value = form.legend;
      } else if (oldMode === 'linear') {
        linearLegendPosition.value = form.legend;
      }

      if (newMode === 'circular') {
        form.legend = circularLegendPosition.value;
      } else if (newMode === 'linear') {
        form.legend = linearLegendPosition.value;
      }

      extractedFeatures.value = [];
      featureRecordIds.value = [];
      selectedFeatureRecordIdx.value = 0;
      Object.keys(featureColorOverrides).forEach((k) => delete featureColorOverrides[k]);
      showFeaturePanel.value = false;
    }
  );

  watch(() => files.d_color, async (newFile) => {
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
    if (!newFile) return;
    try {
      const previousFileRules = manualSpecificRules.filter((r) => r.fromFile);
      for (const rule of previousFileRules) {
        if (rule.cap) {
          await removeLegendEntry(rule.cap);
          fileLegendCaptions.value.delete(rule.cap);
        }
      }
      for (let i = manualSpecificRules.length - 1; i >= 0; i--) {
        if (manualSpecificRules[i].fromFile) {
          manualSpecificRules.splice(i, 1);
        }
      }

      const text = await newFile.text();
      const { rules, rulesWithCaptions, count } = parseSpecificRules(text);

      rules.forEach((rule) => manualSpecificRules.push(rule));
      console.log(`Loaded ${count} rules from file.`);

      if (rulesWithCaptions.length > 0 && pyodideReady.value) {
        await nextTick();
        for (const rule of rulesWithCaptions) {
          const actualCaption = await addLegendEntry(rule.cap, rule.color);
          if (actualCaption && typeof actualCaption === 'string') {
            addedLegendCaptions.value.add(actualCaption);
            fileLegendCaptions.value.add(actualCaption);
          }
        }
        extractLegendEntries();
        console.log(`Added ${rulesWithCaptions.length} legend entries from specific table.`);
      }
    } catch (e) {
      console.error('Failed to load rules file:', e);
      alert('Failed to load rules file. Please check the TSV format.');
    }
  });

  watch(() => files.qualifier_priority, async (newFile) => {
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

  watch(() => form.species, scheduleDefinitionUpdate);
  watch(() => form.strain, scheduleDefinitionUpdate);
  watch(() => state.adv.def_font_size, scheduleDefinitionUpdate);
  watch(() => linearSeqs.map((seq) => seq.definition), scheduleDefinitionUpdate);

  onMounted(async () => {
    await pyodideManager.initPyodide();
  });
};
