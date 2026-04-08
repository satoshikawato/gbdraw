import { resolveColorToHex } from '../color-utils.js';
import { getFeatureCaption } from '../feature-utils.js';

export const createFeatureColorActions = ({
  state,
  nextTick,
  legendActions,
  svgActions,
  ruleActions,
  featureSvgActions
}) => {
  const {
    pyodideReady,
    results,
    selectedResultIndex,
    appliedPaletteColors,
    manualSpecificRules,
    extractedFeatures,
    featureColorOverrides,
    svgContainer,
    clickedFeature,
    colorScopeDialog,
    resetColorDialog,
    legendRenameDialog,
    legendEntries,
    legendStrokeOverrides,
    legendColorOverrides,
    originalLegendOrder,
    originalLegendColors,
    originalSvgStroke,
    skipCaptureBaseConfig,
    skipExtractOnSvgChange,
    addedLegendCaptions
  } = state;

  const {
    addLegendEntry,
    removeLegendEntry,
    updateLegendEntryColorByCaption,
    compactLegendEntries,
    recenterCurrentLegendRoot,
    extractLegendEntries,
    getAllFeatureLegendGroups
  } = legendActions;
  const { applySpecificRulesToSvg } = svgActions;
  const {
    countFeaturesMatchingRule,
    findExistingColorForCaption,
    findFeaturesWithSameDisplayedLabel,
    findFeaturesWithSameIndividualLabel,
    findFeaturesWithSameLegendItem,
    findMatchingRegexRule,
    getDisplayedFeatureLabel,
    getEffectiveLegendCaption,
    getIndividualFeatureLabel,
    getFeatureQualifier
  } = ruleActions;
  const { applyInstantPreview } = featureSvgActions;
  const normalizeCaption = (value) => String(value || '').trim();
  const normalizeCaptionKey = (value) => normalizeCaption(value).toLowerCase();
  const normalizeColor = (value) => String(value || '').trim().toLowerCase();
  const captionsMatch = (left, right) => normalizeCaptionKey(left) === normalizeCaptionKey(right);
  const colorsMatch = (left, right) => normalizeColor(left) === normalizeColor(right);
  const SUFFIXED_CAPTION_PATTERN = /^(.*?)\s*\((\d+)\)$/;

  const findLegendEntryByCaption = (caption) => {
    const normalizedCaption = normalizeCaptionKey(caption);
    if (!normalizedCaption) return null;

    return (
      legendEntries.value.find(
        (entry) => normalizeCaptionKey(entry?.caption) === normalizedCaption
      ) || null
    );
  };

  const findExistingCaptionColor = (feat, caption) => {
    const existingCaption = findExistingColorForCaption(feat, caption);
    if (existingCaption?.color) {
      return {
        caption: existingCaption.rule?.cap || caption,
        color: existingCaption.color,
        rule: existingCaption.rule || null
      };
    }

    const legendEntry = findLegendEntryByCaption(caption);
    if (legendEntry?.color) {
      return {
        caption: legendEntry.caption,
        color: legendEntry.color,
        rule: null
      };
    }

    return null;
  };

  const clearColorScopeDialog = () => {
    colorScopeDialog.show = false;
    colorScopeDialog.feat = null;
    colorScopeDialog.color = null;
    colorScopeDialog.matchingRule = null;
    colorScopeDialog.ruleMatchCount = 0;
    colorScopeDialog.legendName = null;
    colorScopeDialog.siblingCount = 0;
    colorScopeDialog.displayLabel = null;
    colorScopeDialog.displayLabelSiblingCount = 0;
    colorScopeDialog.annotationLabel = null;
    colorScopeDialog.annotationLabelSiblingCount = 0;
    colorScopeDialog.individualLabel = null;
    colorScopeDialog.individualLabelSiblingCount = 0;
    colorScopeDialog.existingCaptionRule = null;
    colorScopeDialog.existingCaptionColor = null;
  };

  const getCurrentSvg = () => svgContainer.value?.querySelector('svg') || null;

  const persistCurrentSvg = (svg = getCurrentSvg()) => {
    if (!svg) return;
    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }
  };

  const findCaptionKey = (store, caption) => {
    if (!store) return null;
    return Object.keys(store).find((key) => captionsMatch(key, caption)) || null;
  };

  const moveCaptionStateKey = (store, oldCaption, newCaption) => {
    if (!store || !oldCaption || !newCaption || oldCaption === newCaption) return;
    const oldKey = findCaptionKey(store, oldCaption);
    if (!oldKey) return;
    const newKey = findCaptionKey(store, newCaption);
    if (!newKey || newKey === oldKey) {
      store[newCaption] = store[oldKey];
    }
    if (oldKey !== newCaption) {
      delete store[oldKey];
    }
  };

  const removeCaptionStateKey = (store, caption) => {
    if (!store || !caption) return;
    const matchingKey = findCaptionKey(store, caption);
    if (matchingKey) {
      delete store[matchingKey];
    }
  };

  const moveAddedLegendCaption = (oldCaption, newCaption) => {
    if (!oldCaption || !newCaption || oldCaption === newCaption) return;
    let matchedCaption = null;
    for (const caption of addedLegendCaptions.value) {
      if (captionsMatch(caption, oldCaption)) {
        matchedCaption = caption;
        break;
      }
    }
    if (!matchedCaption) return;
    addedLegendCaptions.value.delete(matchedCaption);
    addedLegendCaptions.value.add(newCaption);
  };

  const removeAddedLegendCaption = (caption) => {
    if (!caption) return;
    for (const existingCaption of addedLegendCaptions.value) {
      if (captionsMatch(existingCaption, caption)) {
        addedLegendCaptions.value.delete(existingCaption);
        break;
      }
    }
  };

  const syncOriginalLegendMetadataRename = (oldCaption, newCaption, color = null) => {
    if (!oldCaption || !newCaption || oldCaption === newCaption) return;

    const orderIdx = originalLegendOrder.value.findIndex((caption) => captionsMatch(caption, oldCaption));
    if (orderIdx >= 0) {
      originalLegendOrder.value.splice(orderIdx, 1, newCaption);
    }

    const oldColorKey = findCaptionKey(originalLegendColors.value, oldCaption);
    if (!oldColorKey) return;

    const newColorKey = findCaptionKey(originalLegendColors.value, newCaption);
    if (!newColorKey || newColorKey === oldColorKey) {
      originalLegendColors.value[newCaption] = color || originalLegendColors.value[oldColorKey];
    }
    if (oldColorKey !== newCaption) {
      delete originalLegendColors.value[oldColorKey];
    }
  };

  const removeOriginalLegendMetadata = (caption) => {
    if (!caption) return;
    originalLegendOrder.value = originalLegendOrder.value.filter((entryCaption) => !captionsMatch(entryCaption, caption));
    const matchingColorKey = findCaptionKey(originalLegendColors.value, caption);
    if (matchingColorKey) {
      delete originalLegendColors.value[matchingColorKey];
    }
  };

  const updateClickedFeatureLegendState = (feat, caption, color = null) => {
    if (!clickedFeature.value || !feat || clickedFeature.value.svg_id !== feat.svg_id) return;
    if (color) {
      clickedFeature.value.color = color;
    }
    clickedFeature.value.legendName = caption;
    clickedFeature.value.appliedLegendName = caption;
  };

  const clearLegendRenameDialog = ({ restoreInput = false } = {}) => {
    const pendingRequest = legendRenameDialog.pendingRequest;

    if (restoreInput) {
      if (
        pendingRequest?.source === 'popup' &&
        clickedFeature.value &&
        pendingRequest?.feat &&
        clickedFeature.value.svg_id === pendingRequest.feat.svg_id
      ) {
        const fallbackCaption =
          normalizeCaption(clickedFeature.value.appliedLegendName) ||
          normalizeCaption(pendingRequest.oldCaption) ||
          normalizeCaption(getEffectiveLegendCaption(pendingRequest.feat));
        clickedFeature.value.legendName = fallbackCaption;
      } else if (pendingRequest?.source === 'legend') {
        legendEntries.value = [...legendEntries.value];
      }
    }

    legendRenameDialog.show = false;
    legendRenameDialog.mode = 'scope';
    legendRenameDialog.oldCaption = '';
    legendRenameDialog.newCaption = '';
    legendRenameDialog.targetCaption = '';
    legendRenameDialog.targetColor = '';
    legendRenameDialog.currentColor = '';
    legendRenameDialog.siblingCount = 0;
    legendRenameDialog.pendingRequest = null;
  };

  const getCurrentFeatureFillColor = (feat) => {
    if (!feat) return '#cccccc';

    if (clickedFeature.value && clickedFeature.value.svg_id === feat.svg_id && clickedFeature.value.color) {
      return resolveColorToHex(clickedFeature.value.color) || clickedFeature.value.color;
    }

    const svg = getCurrentSvg();
    if (svg && feat.svg_id) {
      const element = svg.querySelector(`#${CSS.escape(feat.svg_id)}`);
      const fill = element?.getAttribute('fill');
      if (fill) {
        return resolveColorToHex(fill) || fill;
      }
    }

    const overrideColor = featureColorOverrides[feat.id]?.color;
    if (overrideColor) {
      return resolveColorToHex(overrideColor) || overrideColor;
    }

    const fallbackColor = appliedPaletteColors.value[feat.type] || '#cccccc';
    return resolveColorToHex(fallbackColor) || fallbackColor;
  };

  const getFeaturesForLegendCaption = (caption) => {
    const normalizedCaption = normalizeCaption(caption);
    if (!normalizedCaption) return [];
    return extractedFeatures.value.filter((feat) => captionsMatch(getEffectiveLegendCaption(feat), normalizedCaption));
  };

  const getUniqueLegendCaption = (caption, options = {}) => {
    const normalizedCaption = normalizeCaption(caption);
    if (!normalizedCaption) return '';

    const ignoreCaptionKeys = new Set(
      (Array.isArray(options.ignoreCaptions) ? options.ignoreCaptions : [])
        .map((value) => normalizeCaptionKey(value))
        .filter(Boolean)
    );

    const existingKeys = new Set();
    legendEntries.value.forEach((entry) => {
      const key = normalizeCaptionKey(entry?.caption);
      if (key) existingKeys.add(key);
    });
    manualSpecificRules.forEach((rule) => {
      const key = normalizeCaptionKey(rule?.cap);
      if (key) existingKeys.add(key);
    });
    Object.values(featureColorOverrides).forEach((override) => {
      const key = normalizeCaptionKey(override?.caption);
      if (key) existingKeys.add(key);
    });

    let finalCaption = normalizedCaption;
    const baseCaption = normalizeCaption(normalizedCaption.replace(/\s*\(\d+\)\s*$/, '')) || normalizedCaption;
    let counter = 1;

    while (existingKeys.has(normalizeCaptionKey(finalCaption)) && !ignoreCaptionKeys.has(normalizeCaptionKey(finalCaption))) {
      finalCaption = `${baseCaption} (${counter})`;
      counter += 1;
    }

    return finalCaption;
  };

  const upsertFeatureHashRule = (feat, color, caption) => {
    if (!feat?.svg_id) return;
    const existingIdx = manualSpecificRules.findIndex(
      (rule) => String(rule.qual || '').toLowerCase() === 'hash' && rule.val === feat.svg_id
    );

    if (existingIdx >= 0) {
      manualSpecificRules[existingIdx].feat = feat.type;
      manualSpecificRules[existingIdx].color = color;
      manualSpecificRules[existingIdx].cap = caption;
      return;
    }

    manualSpecificRules.push({
      feat: feat.type,
      qual: 'hash',
      val: feat.svg_id,
      color,
      cap: caption
    });
  };

  const syncFeatureLegendOverrides = (features, caption, color) => {
    for (const feature of features) {
      upsertFeatureHashRule(feature, color, caption);
      featureColorOverrides[feature.id] = { color, caption };
      updateClickedFeatureLegendState(feature, caption, color);
      applyInstantPreview(feature, color, caption);
    }
  };

  const refreshLegendEntryFeatureIds = (captions = []) => {
    const requestedCaptionKeys =
      Array.isArray(captions) && captions.length > 0
        ? new Set(captions.map((caption) => normalizeCaptionKey(caption)).filter(Boolean))
        : null;

    legendEntries.value.forEach((entry) => {
      if (!entry) return;
      if (requestedCaptionKeys && !requestedCaptionKeys.has(normalizeCaptionKey(entry.caption))) return;
      entry.featureIds = getFeaturesForLegendCaption(entry.caption).map((feature) => feature.svg_id);
    });
  };

  const renameLegendEntryInSvg = (oldCaption, newCaption, color = null) => {
    const svg = getCurrentSvg();
    if (!svg) return false;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;

    let updated = false;

    for (const targetGroup of targetGroups) {
      const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(oldCaption)}"]`);
      if (!entryGroup) continue;

      entryGroup.setAttribute('data-legend-key', newCaption);
      const textEl = entryGroup.querySelector('text');
      if (textEl) {
        textEl.textContent = newCaption;
      }

      if (color) {
        const paths = entryGroup.querySelectorAll('path');
        for (const path of paths) {
          const fill = path.getAttribute('fill');
          if (fill && fill !== 'none' && !fill.startsWith('url(')) {
            path.setAttribute('fill', color);
            break;
          }
        }
      }

      updated = true;
    }

    if (!updated) return false;

    moveCaptionStateKey(legendColorOverrides, oldCaption, newCaption);
    moveCaptionStateKey(legendStrokeOverrides, oldCaption, newCaption);
    moveAddedLegendCaption(oldCaption, newCaption);
    syncOriginalLegendMetadataRename(oldCaption, newCaption, color);

    const legendEntry = legendEntries.value.find((entry) => captionsMatch(entry?.caption, oldCaption));
    if (legendEntry) {
      legendEntry.caption = newCaption;
      legendEntry.originalCaption = newCaption;
      if (color) {
        legendEntry.color = color;
      }
    }

    compactLegendEntries(svg);
    recenterCurrentLegendRoot(svg);
    persistCurrentSvg(svg);
    return true;
  };

  const removeLegendCaptionArtifacts = (caption) => {
    removeCaptionStateKey(legendColorOverrides, caption);
    removeCaptionStateKey(legendStrokeOverrides, caption);
    removeAddedLegendCaption(caption);
    removeOriginalLegendMetadata(caption);
  };

  const ensureLegendEntry = async (caption, color) => {
    const normalizedCaption = normalizeCaption(caption);
    if (!normalizedCaption) return '';

    const existingEntry = findLegendEntryByCaption(normalizedCaption);
    if (existingEntry) {
      if (color && !colorsMatch(existingEntry.color, color)) {
        updateLegendEntryColorByCaption(existingEntry.caption, color);
      }
      return existingEntry.caption;
    }

    const addedCaption = await addLegendEntry(normalizedCaption, color);
    if (addedCaption && typeof addedCaption === 'string') {
      addedLegendCaptions.value.add(addedCaption);
      return addedCaption;
    }

    return normalizedCaption;
  };

  const applyLegendRenameRequest = async (request) => {
    const oldCaption = normalizeCaption(request.oldCaption);
    let finalCaption = normalizeCaption(request.finalCaption || request.newCaption);
    if (!finalCaption) return false;

    const features = Array.isArray(request.features) ? request.features.filter(Boolean) : [];
    const finalColor =
      resolveColorToHex(request.finalColor || request.currentColor) || request.finalColor || request.currentColor || '#cccccc';
    const selectedFeatureIds = new Set(features.map((feature) => feature.svg_id));
    const remainingOldFeatures = oldCaption
      ? getFeaturesForLegendCaption(oldCaption).filter((feature) => !selectedFeatureIds.has(feature.svg_id))
      : [];
    const existingTargetEntry = findLegendEntryByCaption(finalCaption);
    const hasDistinctTargetEntry = existingTargetEntry && !captionsMatch(existingTargetEntry.caption, oldCaption);
    const canRenameInPlace =
      Boolean(oldCaption) && !hasDistinctTargetEntry && remainingOldFeatures.length === 0 && finalCaption !== oldCaption;

    if (features.length > 0) {
      syncFeatureLegendOverrides(features, finalCaption, finalColor);
    }

    let actualCaption = finalCaption;
    if (canRenameInPlace) {
      renameLegendEntryInSvg(oldCaption, finalCaption, finalColor);
    } else if (!hasDistinctTargetEntry) {
      actualCaption = await ensureLegendEntry(finalCaption, finalColor);
      if (actualCaption && actualCaption !== finalCaption && features.length > 0) {
        finalCaption = actualCaption;
        syncFeatureLegendOverrides(features, actualCaption, finalColor);
      }
    } else {
      actualCaption = existingTargetEntry.caption;
      if (features.length > 0) {
        syncFeatureLegendOverrides(features, actualCaption, finalColor);
      }
    }

    applySpecificRulesToSvg();
    await nextTick();

    if (oldCaption && !captionsMatch(actualCaption, oldCaption)) {
      const currentOldUsers = getFeaturesForLegendCaption(oldCaption);
      if (currentOldUsers.length === 0 && !canRenameInPlace) {
        removeLegendEntry(oldCaption);
        removeLegendCaptionArtifacts(oldCaption);
      }
    }

    await nextTick();
    extractLegendEntries();
    refreshLegendEntryFeatureIds([oldCaption, actualCaption]);
    await reclaimOrphanedBaseCaptions();
    await nextTick();
    extractLegendEntries();
    refreshLegendEntryFeatureIds([oldCaption, actualCaption]);

    if (request.source === 'popup' && request.feat) {
      updateClickedFeatureLegendState(request.feat, actualCaption, finalColor);
    }

    legendEntries.value = [...legendEntries.value];
    return true;
  };

  const openLegendRenameScopeDialog = (request, siblingCount) => {
    legendRenameDialog.show = true;
    legendRenameDialog.mode = 'scope';
    legendRenameDialog.oldCaption = request.oldCaption;
    legendRenameDialog.newCaption = request.newCaption;
    legendRenameDialog.targetCaption = '';
    legendRenameDialog.targetColor = '';
    legendRenameDialog.currentColor = request.currentColor || '';
    legendRenameDialog.siblingCount = Math.max(0, siblingCount);
    legendRenameDialog.pendingRequest = request;
  };

  const openLegendRenameTargetDialog = (request, targetEntry) => {
    legendRenameDialog.show = true;
    legendRenameDialog.mode = 'target';
    legendRenameDialog.oldCaption = request.oldCaption;
    legendRenameDialog.newCaption = request.newCaption;
    legendRenameDialog.targetCaption = targetEntry.caption;
    legendRenameDialog.targetColor = targetEntry.color || '';
    legendRenameDialog.currentColor = request.currentColor || '';
    legendRenameDialog.siblingCount = request.siblingCount || 0;
    legendRenameDialog.pendingRequest = request;
  };

  const continueLegendRenameRequest = async (request) => {
    if (!request) return;

    const oldCaption = normalizeCaption(request.oldCaption);
    const newCaption = normalizeCaption(request.newCaption);
    if (!oldCaption || !newCaption || newCaption === oldCaption) {
      clearLegendRenameDialog({ restoreInput: true });
      return;
    }

    const currentColor =
      resolveColorToHex(request.currentColor) ||
      (request.feat ? getCurrentFeatureFillColor(request.feat) : resolveColorToHex(findLegendEntryByCaption(oldCaption)?.color)) ||
      '#cccccc';

    let features = Array.isArray(request.features) ? request.features.filter(Boolean) : [];
    if (!request.sourceScope) {
      const availableFeatures = request.source === 'popup' ? getFeaturesForLegendCaption(oldCaption) : features;
      const siblingCount = Math.max(0, availableFeatures.length - 1);
      if (request.source === 'popup' && siblingCount > 0) {
        openLegendRenameScopeDialog(
          {
            ...request,
            currentColor,
            features: availableFeatures,
            siblingCount
          },
          siblingCount
        );
        return;
      }
      request.sourceScope = request.source === 'popup' ? 'single' : availableFeatures.length > 0 ? 'group' : 'manual';
      features = availableFeatures;
    }

    if (request.sourceScope === 'single') {
      features = request.feat ? [request.feat] : [];
    } else if (request.sourceScope === 'group') {
      features = features.length > 0 ? features : getFeaturesForLegendCaption(oldCaption);
    } else {
      features = [];
    }

    const targetEntry = findLegendEntryByCaption(newCaption);
    const isDistinctTargetEntry = targetEntry && !captionsMatch(targetEntry.caption, oldCaption);

    if (isDistinctTargetEntry && !colorsMatch(targetEntry.color, currentColor)) {
      if (!request.targetResolution) {
        openLegendRenameTargetDialog(
          {
            ...request,
            currentColor,
            features
          },
          targetEntry
        );
        return;
      }

      if (request.targetResolution === 'merge') {
        await applyLegendRenameRequest({
          ...request,
          currentColor,
          features,
          finalCaption: targetEntry.caption,
          finalColor: targetEntry.color
        });
        clearLegendRenameDialog();
        return;
      }

      if (request.targetResolution === 'suffix') {
        await applyLegendRenameRequest({
          ...request,
          currentColor,
          features,
          finalCaption: getUniqueLegendCaption(newCaption, { ignoreCaptions: [oldCaption] }),
          finalColor: currentColor
        });
        clearLegendRenameDialog();
        return;
      }
    }

    const finalCaption =
      isDistinctTargetEntry && colorsMatch(targetEntry.color, currentColor) ? targetEntry.caption : newCaption;
    const finalColor =
      isDistinctTargetEntry && colorsMatch(targetEntry.color, currentColor) ? targetEntry.color : currentColor;

    await applyLegendRenameRequest({
      ...request,
      currentColor,
      features,
      finalCaption,
      finalColor
    });
    clearLegendRenameDialog();
  };

  const reclaimOrphanedBaseCaptions = async () => {
    if (extractedFeatures.value.length === 0) return false;

    const catalog = new Map();
    const rememberCaption = (caption) => {
      const normalized = normalizeCaption(caption);
      if (!normalized) return;
      const key = normalizeCaptionKey(normalized);
      if (!catalog.has(key)) {
        catalog.set(key, normalized);
      }
    };

    legendEntries.value.forEach((entry) => rememberCaption(entry?.caption));
    manualSpecificRules.forEach((rule) => rememberCaption(rule?.cap));
    Object.values(featureColorOverrides).forEach((override) => rememberCaption(override?.caption));

    const usageByCaptionKey = new Map();
    extractedFeatures.value.forEach((feat) => {
      const effectiveCaption = normalizeCaption(getEffectiveLegendCaption(feat));
      if (!effectiveCaption) return;
      rememberCaption(effectiveCaption);
      const key = normalizeCaptionKey(effectiveCaption);
      usageByCaptionKey.set(key, (usageByCaptionKey.get(key) || 0) + 1);
    });

    const suffixCandidatesByBase = new Map();
    for (const [captionKey, captionRaw] of catalog.entries()) {
      const match = captionRaw.match(SUFFIXED_CAPTION_PATTERN);
      if (!match) continue;

      const baseRaw = normalizeCaption(match[1]);
      if (!baseRaw) continue;

      const baseKey = normalizeCaptionKey(baseRaw);
      rememberCaption(baseRaw);
      const parsedIndex = Number.parseInt(match[2], 10);
      const index = Number.isFinite(parsedIndex) ? parsedIndex : Number.MAX_SAFE_INTEGER;

      if (!suffixCandidatesByBase.has(baseKey)) {
        suffixCandidatesByBase.set(baseKey, []);
      }
      suffixCandidatesByBase.get(baseKey).push({
        captionKey,
        captionRaw,
        baseRaw,
        index
      });
    }

    let changed = false;

    for (const [baseKey, candidates] of suffixCandidatesByBase.entries()) {
      const baseUsage = usageByCaptionKey.get(baseKey) || 0;
      if (baseUsage > 0) continue;

      const ordered = [...candidates].sort((a, b) => {
        if (a.index !== b.index) return a.index - b.index;
        return a.captionRaw.localeCompare(b.captionRaw);
      });

      const activeCandidate =
        ordered.find((candidate) => (usageByCaptionKey.get(candidate.captionKey) || 0) > 0) || ordered[0];
      if (!activeCandidate) continue;
      if (activeCandidate.captionKey === baseKey) continue;

      const baseCaption = catalog.get(baseKey) || activeCandidate.baseRaw;
      const suffixCaption = activeCandidate.captionRaw;

      let sourceColor = null;
      const suffixLegendEntry = findLegendEntryByCaption(suffixCaption);
      if (suffixLegendEntry?.color) {
        sourceColor = suffixLegendEntry.color;
      }

      manualSpecificRules.forEach((rule) => {
        if (!captionsMatch(rule.cap, suffixCaption)) return;
        if (!sourceColor && rule.color) sourceColor = rule.color;
        rule.cap = baseCaption;
        changed = true;
      });

      Object.values(featureColorOverrides).forEach((override) => {
        if (!override || !captionsMatch(override.caption, suffixCaption)) return;
        if (!sourceColor && override.color) sourceColor = override.color;
        override.caption = baseCaption;
        changed = true;
      });

      const baseLegendEntry = findLegendEntryByCaption(baseCaption);
      if (baseLegendEntry) {
        if (sourceColor) {
          updateLegendEntryColorByCaption(baseLegendEntry.caption, sourceColor);
          changed = true;
        }
      } else {
        const colorToUse = sourceColor || '#cccccc';
        const addedCaption = await addLegendEntry(baseCaption, colorToUse);
        if (addedCaption && typeof addedCaption === 'string') {
          addedLegendCaptions.value.add(addedCaption);
          const addedKey = normalizeCaptionKey(addedCaption);
          if (addedKey !== baseKey) {
            manualSpecificRules.forEach((rule) => {
              if (captionsMatch(rule.cap, baseCaption)) {
                rule.cap = addedCaption;
              }
            });
            Object.values(featureColorOverrides).forEach((override) => {
              if (override && captionsMatch(override.caption, baseCaption)) {
                override.caption = addedCaption;
              }
            });
          }
          changed = true;
        }
      }

      removeLegendEntry(suffixLegendEntry?.caption || suffixCaption);
      changed = true;
    }

    if (changed) {
      applySpecificRulesToSvg();
      await nextTick();
      extractLegendEntries();
    }

    return changed;
  };

  const applyColorToFeatureGroup = async (features, targetCaption, color) => {
    if (!Array.isArray(features) || features.length === 0) return;

    const normalizedTargetCaption = normalizeCaption(targetCaption);
    if (!normalizedTargetCaption) return;

    const existingLegendEntry = findLegendEntryByCaption(normalizedTargetCaption);
    let finalCaption = existingLegendEntry?.caption || normalizedTargetCaption;

    if (existingLegendEntry) {
      updateLegendEntryColorByCaption(existingLegendEntry.caption, color);
    } else {
      const addedCaption = await addLegendEntry(normalizedTargetCaption, color);
      if (addedCaption && typeof addedCaption === 'string') {
        finalCaption = addedCaption;
        addedLegendCaptions.value.add(addedCaption);
      }
    }

    for (const feature of features) {
      const existingIdx = manualSpecificRules.findIndex(
        (rule) => String(rule.qual || '').toLowerCase() === 'hash' && rule.val === feature.svg_id
      );
      if (existingIdx >= 0) {
        manualSpecificRules[existingIdx].color = color;
        manualSpecificRules[existingIdx].cap = finalCaption;
      } else {
        manualSpecificRules.push({
          feat: feature.type,
          qual: 'hash',
          val: feature.svg_id,
          color: color,
          cap: finalCaption
        });
      }
      featureColorOverrides[feature.id] = { color, caption: finalCaption };
      updateClickedFeatureLegendState(feature, finalCaption, color);
    }

    applySpecificRulesToSvg();
    await reclaimOrphanedBaseCaptions();
    extractLegendEntries();
  };

  const requestFeatureColorChange = async (feat, color, requestedLegendName = null, options = {}) => {
    if (!feat) return;
    const requestedCaption = normalizeCaption(requestedLegendName);
    const effectiveCaption = normalizeCaption(getEffectiveLegendCaption(feat));
    const fallbackCaption = normalizeCaption(getFeatureCaption(feat));
    const legendName = requestedCaption || effectiveCaption || fallbackCaption;
    if (!legendName) return;

    const matchingRule = findMatchingRegexRule(feat);
    const ruleMatchCount = matchingRule ? countFeaturesMatchingRule(matchingRule) : 0;

    const siblings = findFeaturesWithSameLegendItem(feat, legendName);
    const siblingCount = siblings.length;

    const existingCaption = findExistingCaptionColor(feat, legendName);
    const displayLabel = normalizeCaption(getDisplayedFeatureLabel(feat));
    const displayLabelSiblingCount = displayLabel ? findFeaturesWithSameDisplayedLabel(feat, displayLabel).length : 0;
    const annotationLabel = normalizeCaption(getIndividualFeatureLabel(feat));
    const annotationLabelSiblingCount = annotationLabel
      ? findFeaturesWithSameIndividualLabel(feat, annotationLabel).length
      : 0;

    const needsDialog =
      Boolean(matchingRule) || siblingCount > 0 || displayLabelSiblingCount > 0 || annotationLabelSiblingCount > 0;
    if (needsDialog) {
      colorScopeDialog.show = true;
      colorScopeDialog.feat = feat;
      colorScopeDialog.color = color;
      colorScopeDialog.matchingRule = matchingRule;
      colorScopeDialog.ruleMatchCount = ruleMatchCount;
      colorScopeDialog.legendName = legendName;
      colorScopeDialog.siblingCount = siblingCount;
      colorScopeDialog.displayLabel = displayLabel;
      colorScopeDialog.displayLabelSiblingCount = displayLabelSiblingCount;
      colorScopeDialog.annotationLabel = annotationLabel;
      colorScopeDialog.annotationLabelSiblingCount = annotationLabelSiblingCount;
      // Backward-compatible aliases
      colorScopeDialog.individualLabel = annotationLabel;
      colorScopeDialog.individualLabelSiblingCount = annotationLabelSiblingCount;
      colorScopeDialog.existingCaptionRule = existingCaption?.rule || null;
      colorScopeDialog.existingCaptionColor = existingCaption?.color || null;
      if (options.closePopupOnDialog) {
        clickedFeature.value = null;
      }
      return;
    }

    if (clickedFeature.value && clickedFeature.value.svg_id === feat.svg_id) {
      clickedFeature.value.color = color;
      if (requestedCaption) {
        clickedFeature.value.legendName = requestedCaption;
      }
    }
    await setFeatureColor(feat, color, legendName);
  };

  const updateClickedFeatureColor = async (color) => {
    if (!clickedFeature.value) return;
    const feat = clickedFeature.value.feat;
    if (!feat) return;
    const customName = normalizeCaption(clickedFeature.value.legendName);
    await requestFeatureColorChange(feat, color, customName, { closePopupOnDialog: true });
  };

  const handleLegendNameCommit = async () => {
    if (!clickedFeature.value) return;

    const feat = clickedFeature.value.feat;
    if (!feat) return;

    const requestedCaption = normalizeCaption(clickedFeature.value.legendName);
    const currentCaption =
      normalizeCaption(clickedFeature.value.appliedLegendName) || normalizeCaption(getEffectiveLegendCaption(feat));

    if (!requestedCaption) {
      clickedFeature.value.legendName = currentCaption;
      return;
    }

    if (!currentCaption || requestedCaption === currentCaption) {
      clickedFeature.value.legendName = currentCaption || requestedCaption;
      return;
    }

    await continueLegendRenameRequest({
      source: 'popup',
      feat,
      oldCaption: currentCaption,
      newCaption: requestedCaption,
      currentColor: getCurrentFeatureFillColor(feat),
      sourceScope: null
    });
  };

  const selectLegendNameOption = async (caption) => {
    if (!clickedFeature.value) return;
    const selectedCaption = String(caption || '').trim();
    if (!selectedCaption) return;
    clickedFeature.value.legendName = selectedCaption;
    await handleLegendNameCommit();
  };

  const handleLegendRenameChoice = async (choice) => {
    const pendingRequest = legendRenameDialog.pendingRequest;
    if (!pendingRequest || choice === 'cancel') {
      clearLegendRenameDialog({ restoreInput: true });
      return;
    }

    if (legendRenameDialog.mode === 'scope') {
      if (choice === 'single') {
        await continueLegendRenameRequest({
          ...pendingRequest,
          sourceScope: 'single',
          targetResolution: null
        });
        return;
      }

      if (choice === 'group') {
        await continueLegendRenameRequest({
          ...pendingRequest,
          sourceScope: 'group',
          targetResolution: null
        });
        return;
      }
    }

    if (legendRenameDialog.mode === 'target') {
      if (choice === 'merge') {
        await continueLegendRenameRequest({
          ...pendingRequest,
          targetResolution: 'merge'
        });
        return;
      }

      if (choice === 'suffix') {
        await continueLegendRenameRequest({
          ...pendingRequest,
          targetResolution: 'suffix'
        });
        return;
      }
    }

    clearLegendRenameDialog({ restoreInput: true });
  };

  const renameLegendEntry = async (idx, newCaption) => {
    const entry = legendEntries.value[idx];
    if (!entry) return;

    const requestedCaption = normalizeCaption(newCaption);
    if (!requestedCaption || requestedCaption === normalizeCaption(entry.caption)) {
      legendEntries.value = [...legendEntries.value];
      return;
    }

    await continueLegendRenameRequest({
      source: 'legend',
      oldCaption: normalizeCaption(entry.caption),
      newCaption: requestedCaption,
      currentColor: resolveColorToHex(entry.color) || entry.color || '#cccccc',
      features: getFeaturesForLegendCaption(entry.caption),
      sourceScope: getFeaturesForLegendCaption(entry.caption).length > 0 ? 'group' : 'manual'
    });
  };

  const handleColorScopeChoice = async (choice) => {
    const { feat, color, matchingRule, legendName, existingCaptionColor } = colorScopeDialog;
    if (choice === 'cancel' || !feat || !color) {
      clearColorScopeDialog();
      return;
    }

    if (choice === 'rule') {
      if (matchingRule) {
        matchingRule.color = color;
        if (matchingRule.cap) {
          updateLegendEntryColorByCaption(matchingRule.cap, color);
          extractLegendEntries();
        }
      }
      applySpecificRulesToSvg();
    } else if (choice === 'caption') {
      const targetLegendName = normalizeCaption(legendName) || normalizeCaption(getEffectiveLegendCaption(feat));
      if (!targetLegendName) {
        clearColorScopeDialog();
        return;
      }
      const siblings = findFeaturesWithSameLegendItem(feat, targetLegendName);
      const allFeatures = [feat, ...siblings];
      await applyColorToFeatureGroup(allFeatures, targetLegendName, color);
    } else if (choice === 'displayLabel') {
      const displayLabel =
        normalizeCaption(colorScopeDialog.displayLabel) || normalizeCaption(getDisplayedFeatureLabel(feat));
      if (!displayLabel) {
        clearColorScopeDialog();
        return;
      }
      const displaySiblings = findFeaturesWithSameDisplayedLabel(feat, displayLabel);
      const allFeatures = [feat, ...displaySiblings];
      await applyColorToFeatureGroup(allFeatures, displayLabel, color);
    } else if (choice === 'single') {
      let singleCaption = legendName;
      if (matchingRule && colorScopeDialog.ruleMatchCount > 1) {
        const ruleCaption = matchingRule.cap || matchingRule.val;
        if (legendName === ruleCaption) {
          singleCaption = getIndividualFeatureLabel(feat);
        }
      }
      await setFeatureColor(feat, color, singleCaption);
    } else if (choice === 'annotationLabel' || choice === 'individualLabel') {
      const annotationLabel =
        normalizeCaption(colorScopeDialog.annotationLabel) || normalizeCaption(getIndividualFeatureLabel(feat));
      if (!annotationLabel) {
        clearColorScopeDialog();
        return;
      }
      const annotationSiblings = findFeaturesWithSameIndividualLabel(feat, annotationLabel);
      const allFeatures = [feat, ...annotationSiblings];
      await applyColorToFeatureGroup(allFeatures, annotationLabel, color);
    } else if (choice === 'useExisting') {
      if (existingCaptionColor) {
        const targetLegendName = normalizeCaption(legendName) || normalizeCaption(getEffectiveLegendCaption(feat));
        const existingRule = manualSpecificRules.find(
          (r) => String(r.qual || '').toLowerCase() === 'hash' && r.val === feat.svg_id
        );
        if (existingRule) {
          existingRule.color = existingCaptionColor;
        } else {
          manualSpecificRules.push({
            feat: feat.type,
            qual: 'hash',
            val: feat.svg_id,
            color: existingCaptionColor,
            cap: targetLegendName
          });
        }
        if (clickedFeature.value && clickedFeature.value.svg_id === feat.svg_id) {
          clickedFeature.value.color = existingCaptionColor;
          clickedFeature.value.legendName = targetLegendName;
          clickedFeature.value.appliedLegendName = targetLegendName;
        }
        featureColorOverrides[feat.id] = { color: existingCaptionColor, caption: targetLegendName };
        applyInstantPreview(feat, existingCaptionColor, targetLegendName);
        applySpecificRulesToSvg();
      }
    }

    clearColorScopeDialog();
  };

  const updateClickedFeatureStroke = (strokeColor, strokeWidth) => {
    if (!clickedFeature.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const svgId = clickedFeature.value.svg_id;
    const elements = svg.querySelectorAll(`#${CSS.escape(svgId)}`);

    elements.forEach((el) => {
      if (strokeColor !== null) {
        el.setAttribute('stroke', strokeColor);
        clickedFeature.value.strokeColor = strokeColor;
      }
      if (strokeWidth !== null) {
        const widthVal = parseFloat(strokeWidth);
        if (!isNaN(widthVal)) {
          el.setAttribute('stroke-width', widthVal);
          clickedFeature.value.strokeWidth = widthVal;
        }
      }
    });

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }
  };

  const resetClickedFeatureStroke = () => {
    if (!clickedFeature.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const svgId = clickedFeature.value.svg_id;
    const elements = svg.querySelectorAll(`#${CSS.escape(svgId)}`);

    const originalColor = originalSvgStroke.value.color;
    const originalWidth = originalSvgStroke.value.width;

    elements.forEach((el) => {
      if (originalColor === null) {
        el.removeAttribute('stroke');
      } else {
        el.setAttribute('stroke', originalColor);
      }
      if (originalWidth === null) {
        el.removeAttribute('stroke-width');
      } else {
        el.setAttribute('stroke-width', originalWidth);
      }
    });

    clickedFeature.value.strokeColor = originalColor || '';
    clickedFeature.value.strokeWidth = originalWidth ?? '';

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }
  };

  const resetClickedFeatureFillColor = () => {
    if (!clickedFeature.value) return;
    if (!svgContainer.value) return;

    const feat = clickedFeature.value.feat;
    if (!feat) return;

    const defaultColor = appliedPaletteColors.value[feat.type];
    if (!defaultColor) {
      console.warn('No default color found for feature type:', feat.type);
      return;
    }

    const caption = getFeatureCaption(feat);

    const siblings = extractedFeatures.value.filter(
      (f) => getFeatureCaption(f) === caption && f.svg_id !== clickedFeature.value.svg_id
    );

    if (siblings.length > 0) {
      resetColorDialog.show = true;
      resetColorDialog.caption = caption;
      resetColorDialog.defaultColor = defaultColor;
      resetColorDialog.siblingCount = siblings.length;
    } else {
      doResetFillColor('this');
    }
  };

  const handleResetColorChoice = async (choice) => {
    resetColorDialog.show = false;
    await doResetFillColor(choice);
  };

  const doResetFillColor = async (choice) => {
    if (!clickedFeature.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const feat = clickedFeature.value.feat;
    if (!feat) return;

    const defaultColor = resetColorDialog.defaultColor || appliedPaletteColors.value[feat.type];
    const caption =
      resetColorDialog.caption ||
      feat.product ||
      feat.gene ||
      feat.locus_tag ||
      feat.note ||
      `${feat.type} at ${feat.start}..${feat.end}`;
    const svgId = clickedFeature.value.svg_id;

    if (choice === 'this' || choice === 'this_with_legend') {
      const elements = svg.querySelectorAll(`#${CSS.escape(svgId)}`);
      elements.forEach((el) => {
        el.setAttribute('fill', defaultColor);
      });

      clickedFeature.value.color = defaultColor;

      const ruleIdx = manualSpecificRules.findIndex((r) => r.qual === 'hash' && r.val === svgId);
      if (ruleIdx !== -1) {
        manualSpecificRules.splice(ruleIdx, 1);
      }

      if (choice === 'this_with_legend') {
        console.log(`Attempting to add legend entry: caption="${caption}", color="${defaultColor}"`);
        const addedCaption = await addLegendEntry(caption, defaultColor);
        console.log(`addLegendEntry returned: ${addedCaption}`);
        if (addedCaption) {
          extractLegendEntries();
          console.log(`Added legend entry: ${addedCaption} with color: ${defaultColor}`);
        } else {
          console.error(`Failed to add legend entry for caption="${caption}"`);
        }
      }

      console.log(`Reset fill color to default (${defaultColor}) for feature: ${svgId}`);
    } else if (choice === 'all') {
      const matchingFeatures = extractedFeatures.value.filter((f) => getFeatureCaption(f) === caption);

      for (const matchFeat of matchingFeatures) {
        const elements = svg.querySelectorAll(`#${CSS.escape(matchFeat.svg_id)}`);
        elements.forEach((el) => {
          el.setAttribute('fill', defaultColor);
        });
      }

      clickedFeature.value.color = defaultColor;

      for (let i = manualSpecificRules.length - 1; i >= 0; i--) {
        const rule = manualSpecificRules[i];
        if (rule.cap === caption) {
          manualSpecificRules.splice(i, 1);
        }
      }

      const legendIdx = legendEntries.value.findIndex((e) => e.caption === caption);
      if (legendIdx !== -1) {
        await removeLegendEntry(caption);
        extractLegendEntries();
      }

      console.log(
        `Reset fill color to default (${defaultColor}) for all ${matchingFeatures.length} features with caption: ${caption}`
      );
    }

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }

    clickedFeature.value = null;
  };

  const applyStrokeToAllSiblings = async () => {
    if (!clickedFeature.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const currentStrokeColor = clickedFeature.value.strokeColor;
    const currentStrokeWidth = clickedFeature.value.strokeWidth;
    const currentFillColor = clickedFeature.value.color;
    const clickedSvgId = clickedFeature.value.svg_id;

    const feat = clickedFeature.value.feat;
    if (!feat) {
      console.warn('No feature data found in clickedFeature');
      return;
    }
    const caption = getFeatureCaption(feat);

    let targetLegendEntry = legendEntries.value.find((e) => e.featureIds && e.featureIds.includes(clickedSvgId));

    if (!targetLegendEntry) {
      targetLegendEntry = legendEntries.value.find(
        (e) => e.color && e.color.toLowerCase() === currentFillColor?.toLowerCase()
      );
    }

    let siblingFeatureIds = [];
    if (targetLegendEntry && targetLegendEntry.featureIds && targetLegendEntry.featureIds.length > 0) {
      siblingFeatureIds = targetLegendEntry.featureIds;
      console.log(`Found legend entry "${targetLegendEntry.caption}" with ${siblingFeatureIds.length} features`);
    } else {
      const siblingFeatures = extractedFeatures.value.filter((f) => {
        if (getFeatureCaption(f) !== caption) return false;
        const el = svg.querySelector(`#${CSS.escape(f.svg_id)}`);
        if (!el) return false;
        const fillColor = el.getAttribute('fill');
        return fillColor && fillColor.toLowerCase() === currentFillColor?.toLowerCase();
      });
      siblingFeatureIds = siblingFeatures.map((f) => f.svg_id);
      console.log(
        `Fallback: Found ${siblingFeatureIds.length} features by caption="${caption}" and color="${currentFillColor}"`
      );
    }

    for (const svgId of siblingFeatureIds) {
      const elements = svg.querySelectorAll(`#${CSS.escape(svgId)}`);
      elements.forEach((el) => {
        if (currentStrokeColor) {
          el.setAttribute('stroke', currentStrokeColor);
        }
        if (currentStrokeWidth !== null && currentStrokeWidth !== '') {
          el.setAttribute('stroke-width', parseFloat(currentStrokeWidth));
        }
      });
    }

    console.log(`Applied stroke to ${siblingFeatureIds.length} features`);

    if (targetLegendEntry) {
      const legendGroups = getAllFeatureLegendGroups(svg);
      for (const targetGroup of legendGroups) {
        const entryGroup = targetGroup.querySelector(
          `g[data-legend-key="${CSS.escape(targetLegendEntry.caption)}"]`
        );
        if (entryGroup) {
          const paths = entryGroup.querySelectorAll('path');
          for (const path of paths) {
            const fill = path.getAttribute('fill');
            if (fill && fill !== 'none' && !fill.startsWith('url(')) {
              if (currentStrokeColor) {
                path.setAttribute('stroke', currentStrokeColor);
              }
              if (currentStrokeWidth !== null && currentStrokeWidth !== '') {
                path.setAttribute('stroke-width', parseFloat(currentStrokeWidth));
              }
              console.log(`Updated legend rect stroke for "${targetLegendEntry.caption}"`);
              break;
            }
          }
        }
      }
    }

    const overrideKey = targetLegendEntry ? targetLegendEntry.caption : caption;

    if (!targetLegendEntry) {
      const fillColor = clickedFeature.value.color || '#cccccc';
      const addedCaption = await addLegendEntry(caption, fillColor);
      if (addedCaption) {
        extractLegendEntries();
        targetLegendEntry = legendEntries.value.find((e) => e.caption === addedCaption);
        if (targetLegendEntry) {
          targetLegendEntry.showStroke = true;
          if (!targetLegendEntry.featureIds) targetLegendEntry.featureIds = [];
          if (!targetLegendEntry.featureIds.includes(clickedSvgId)) {
            targetLegendEntry.featureIds.push(clickedSvgId);
          }
        }
        console.log(`Created new legend entry for caption: ${addedCaption}`);
      } else {
        targetLegendEntry = {
          caption: caption,
          originalCaption: caption,
          color: fillColor,
          visible: true,
          showStroke: true,
          featureIds: [clickedSvgId]
        };
        legendEntries.value.push(targetLegendEntry);
      }
    }

    if (!legendStrokeOverrides[overrideKey]) {
      legendStrokeOverrides[overrideKey] = {};
    }
    if (currentStrokeColor) {
      legendStrokeOverrides[overrideKey].strokeColor = currentStrokeColor;
    }
    if (currentStrokeWidth !== null && currentStrokeWidth !== '') {
      legendStrokeOverrides[overrideKey].strokeWidth = parseFloat(currentStrokeWidth);
    }

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }

    console.log(
      `Applied stroke (color: ${currentStrokeColor}, width: ${currentStrokeWidth}) to ${siblingFeatureIds.length} features`
    );
  };

  const setFeatureColor = async (feat, color, customCaption = null) => {
    const qualInfo = getFeatureQualifier(feat);
    if (!qualInfo) {
      console.warn(
        `Cannot identify feature: ${feat.type} at ${feat.start}..${feat.end} (no locus_tag, gene, or product)`
      );
      return;
    }
    const { qual, val } = qualInfo;

    const featureKey = feat.id;

    const caption = normalizeCaption(
      customCaption || feat.product || feat.gene || feat.locus_tag || `${feat.type} at ${feat.start}..${feat.end}`
    );
    if (!caption) return;

    const oldOverride = featureColorOverrides[featureKey];
    const oldCaption = normalizeCaption(oldOverride?.caption);
    featureColorOverrides[featureKey] = { color, caption };

    applyInstantPreview(feat, color, caption);

    await nextTick();
    let actualCaption = caption;
    if (pyodideReady.value && caption) {
      if (oldCaption) {
        if (captionsMatch(oldCaption, caption)) {
          const hasNonHashRule = manualSpecificRules.some(
            (rule) => captionsMatch(rule.cap, caption) && String(rule.qual || '').toLowerCase() !== 'hash'
          );
          const hasOtherUses = extractedFeatures.value.some((f) => {
            if (f.svg_id === feat.svg_id) return false;
            return captionsMatch(getEffectiveLegendCaption(f), caption);
          });

          if (hasNonHashRule || hasOtherUses) {
            actualCaption = await addLegendEntry(caption, color);
            if (actualCaption && typeof actualCaption === 'string') {
              addedLegendCaptions.value.add(actualCaption);
            }
          } else {
            updateLegendEntryColorByCaption(oldCaption, color);
          }
        } else {
          removeLegendEntry(oldCaption);
          actualCaption = await addLegendEntry(caption, color);
          if (actualCaption && typeof actualCaption === 'string') {
            addedLegendCaptions.value.add(actualCaption);
          }
        }
      } else {
        actualCaption = await addLegendEntry(caption, color);
        if (actualCaption && typeof actualCaption === 'string') {
          addedLegendCaptions.value.add(actualCaption);
        }
      }

      if (actualCaption && typeof actualCaption === 'string' && actualCaption !== caption) {
        featureColorOverrides[featureKey] = { color, caption: actualCaption };
      }
    }

    skipExtractOnSvgChange.value = true;

    const finalCaption = actualCaption && typeof actualCaption === 'string' ? actualCaption : caption;
    const existingIdx = manualSpecificRules.findIndex(
      (r) => r.feat === feat.type && String(r.qual || '').toLowerCase() === String(qual || '').toLowerCase() && r.val === val
    );

    if (existingIdx >= 0) {
      manualSpecificRules[existingIdx].color = color;
      manualSpecificRules[existingIdx].cap = finalCaption;
    } else {
      manualSpecificRules.push({
        feat: feat.type,
        qual: qual,
        val: val,
        color: color,
        cap: finalCaption
      });
    }

    featureColorOverrides[featureKey] = { color, caption: finalCaption };
    updateClickedFeatureLegendState(feat, finalCaption, color);

    await reclaimOrphanedBaseCaptions();

    await nextTick();
    await nextTick();
    skipExtractOnSvgChange.value = false;
    extractLegendEntries();
  };

  return {
    applyStrokeToAllSiblings,
    handleColorScopeChoice,
    handleLegendNameCommit,
    handleLegendRenameChoice,
    renameLegendEntry,
    requestFeatureColorChange,
    selectLegendNameOption,
    handleResetColorChoice,
    resetClickedFeatureFillColor,
    resetClickedFeatureStroke,
    setFeatureColor,
    updateClickedFeatureColor,
    updateClickedFeatureStroke
  };
};
