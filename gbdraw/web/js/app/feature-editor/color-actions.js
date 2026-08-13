import { resolveColorToHex } from '../color-utils.js';
import { getFeatureCaption, getFeatureHashCandidates, ruleMatchesFeature } from '../feature-utils.js';
import { exactRegexValue } from '../feature-selector.js';
import {
  featureOverrideKey,
  getFeatureOverride
} from '../../services/feature-override-identity.js';

export const createFeatureColorActions = ({
  state,
  nextTick,
  legendActions,
  svgActions,
  ruleActions,
  featureSvgActions,
  previewRuntime = null
}) => {
  const {
    appliedPaletteColors,
    manualSpecificRules,
    extractedFeatures,
    biologicalFeatures,
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
    featureStrokeOverrides,
    skipExtractOnSvgChange,
    addedLegendCaptions
  } = state;

  const {
    addLegendEntry: addLegendEntryRaw,
    removeLegendEntry: removeLegendEntryRaw,
    updateLegendEntryColorByCaption: updateLegendEntryColorByCaptionRaw,
    compactLegendEntries,
    extractLegendEntries,
    getAllFeatureLegendGroups,
    onLegendGeometryChanged
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
    getFeatureQualifier,
    getLabelSpecificRule
  } = ruleActions;
  const { getFeatureElements, getFeatureFillElements } = featureSvgActions;
  const normalizeCaption = (value) => String(value || '').trim();
  const normalizeCaptionKey = (value) => normalizeCaption(value).toLowerCase();
  const normalizeColor = (value) => String(value || '').trim().toLowerCase();
  const captionsMatch = (left, right) => normalizeCaptionKey(left) === normalizeCaptionKey(right);
  const colorsMatch = (left, right) => normalizeColor(left) === normalizeColor(right);
  const isHashSpecificRule = (rule) => String(rule?.qual || '').toLowerCase() === 'hash';
  const SUFFIXED_CAPTION_PATTERN = /^(.*?)\s*\((\d+)\)$/;
  const hasOwn = (object, key) => Object.prototype.hasOwnProperty.call(object || {}, key);
  let colorActionDepth = 0;
  let colorActionChangedPreview = false;

  const markColorPreviewDirty = (reason = 'feature-color') => {
    const marked = previewRuntime?.markActiveResultDirty?.(reason) === true;
    colorActionChangedPreview = colorActionChangedPreview || marked;
    return marked;
  };

  const applyFeatureColorPreview = (feature, color) => {
    const featureId = String(
      feature?.rendered_svg_id
      || feature?.renderedSvgId
      || feature?.rendered_feature_svg_id
      || feature?.renderedFeatureSvgId
      || feature?.svg_id
      || ''
    ).trim();
    if (!featureId || !previewRuntime?.applyFeatureFillChanges) return false;
    const updated = previewRuntime.applyFeatureFillChanges(
      [{ featureId, color }],
      { reason: 'feature-color' }
    ) === true;
    colorActionChangedPreview = colorActionChangedPreview || updated;
    return updated;
  };

  const runColorAction = async (action) => {
    const isOuterAction = colorActionDepth === 0;
    const wasDirty = Boolean(previewRuntime?.getActiveRuntime?.()?.dirty);
    if (isOuterAction) colorActionChangedPreview = false;
    colorActionDepth += 1;
    let completed = false;
    try {
      const result = await action();
      completed = true;
      return result;
    } finally {
      colorActionDepth -= 1;
      if (isOuterAction) {
        const becameDirty = !wasDirty && Boolean(previewRuntime?.getActiveRuntime?.()?.dirty);
        if (completed && (colorActionChangedPreview || becameDirty) && previewRuntime?.flushActiveResult) {
          previewRuntime.flushActiveResult();
        }
        colorActionChangedPreview = false;
      }
    }
  };

  const colorAction = (action) => (...args) => runColorAction(() => action(...args));

  const hashRuleTargetsFeatureExactly = (rule, feature) => {
    if (!isHashSpecificRule(rule) || rule?.feat !== feature?.type) return false;
    const ruleValue = String(rule?.val || '').trim();
    const candidates = getFeatureHashCandidates(feature);
    const generationHash = candidates[0] || '';
    const renderedId = candidates[candidates.length - 1] || '';
    const isExact = (candidate) => (
      Boolean(candidate) && (ruleValue === candidate || ruleValue === exactRegexValue(candidate))
    );
    if (renderedId !== generationHash && isExact(renderedId)) return true;
    if (!isExact(generationHash)) return false;
    const collisionCount = extractedFeatures.value.filter(
      (candidate) => candidate?.type === feature?.type && getFeatureHashCandidates(candidate)[0] === generationHash
    ).length;
    return collisionCount <= 1;
  };

  const normalizeStrokeWidthValue = (value) => {
    if (value === null || value === undefined || value === '') return null;
    const numeric = Number(value);
    return Number.isFinite(numeric) && numeric >= 0 ? numeric : null;
  };

  const strokeColorAttributeMatches = (element, color) => {
    const current = element?.getAttribute?.('stroke');
    return color === null ? current === null : colorsMatch(current, color);
  };

  const strokeWidthAttributeMatches = (element, width) => {
    const current = element?.getAttribute?.('stroke-width');
    if (width === null) return current === null;
    const currentNumeric = Number(current);
    return Number.isFinite(currentNumeric) && currentNumeric === Number(width);
  };

  const findFeatureBySvgId = (svgId) => {
    const normalizedSvgId = String(svgId || '').trim();
    if (!normalizedSvgId) return null;
    return extractedFeatures.value.find((feature) => String(feature?.svg_id || '').trim() === normalizedSvgId) || null;
  };

  const featureStrokeKey = (featureLike, fallbackSvgId = '') => {
    const feature = featureLike?.feat || featureLike;
    return featureOverrideKey(feature) || String(fallbackSvgId || '').trim();
  };

  const recordFeatureStrokeOverride = (
    featureLike,
    { strokeColor = null, strokeWidth = null, originalStrokeColor = null, originalStrokeWidth = null } = {}
  ) => {
    const key = featureStrokeKey(featureLike, featureLike?.svg_id);
    if (!key) return;

    const existing = featureStrokeOverrides[key] || {};
    const next = { ...existing };
    if (!hasOwn(next, 'originalStrokeColor')) {
      next.originalStrokeColor = originalStrokeColor;
    }
    if (!hasOwn(next, 'originalStrokeWidth')) {
      next.originalStrokeWidth = normalizeStrokeWidthValue(originalStrokeWidth);
    }
    if (strokeColor !== null && strokeColor !== undefined && strokeColor !== '') {
      next.strokeColor = strokeColor;
    }
    const widthVal = normalizeStrokeWidthValue(strokeWidth);
    if (widthVal !== null) {
      next.strokeWidth = widthVal;
    }
    if (hasOwn(next, 'strokeColor') || hasOwn(next, 'strokeWidth')) {
      featureStrokeOverrides[key] = next;
    }
  };

  const clearFeatureStrokeOverride = (featureLike, fallbackSvgId = '') => {
    const key = featureStrokeKey(featureLike, fallbackSvgId);
    if (key) delete featureStrokeOverrides[key];
  };

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
    colorScopeDialog.existingCaptionRule = null;
    colorScopeDialog.existingCaptionColor = null;
  };

  const getCurrentSvg = () => svgContainer.value?.querySelector('svg') || null;

  const persistCurrentSvg = (svg = getCurrentSvg(), reason = 'feature-color') => {
    if (!svg) return false;
    return markColorPreviewDirty(reason);
  };

  const getLiveLegendColor = (caption) => {
    const svg = getCurrentSvg();
    const normalizedCaption = normalizeCaption(caption);
    if (!svg || !normalizedCaption) return null;
    const escapedCaption = globalThis.CSS?.escape
      ? globalThis.CSS.escape(normalizedCaption)
      : normalizedCaption.replace(/["\\]/g, '\\$&');
    for (const targetGroup of getAllFeatureLegendGroups(svg)) {
      const entryGroup = targetGroup.querySelector(`g[data-legend-key="${escapedCaption}"]`);
      if (!entryGroup) continue;
      const colorPath = Array.from(entryGroup.querySelectorAll('path')).find((path) => {
        const fill = path.getAttribute('fill');
        return fill && fill !== 'none' && !fill.startsWith('url(');
      });
      if (colorPath) return colorPath.getAttribute('fill');
    }
    return null;
  };

  const addLegendEntry = async (caption, color, options = {}) => {
    const beforeColor = getLiveLegendColor(caption);
    const addedCaption = await addLegendEntryRaw(caption, color, { ...options, commit: false });
    if (!addedCaption) return addedCaption;
    if (
      !beforeColor
      || !captionsMatch(addedCaption, caption)
      || !colorsMatch(beforeColor, color)
    ) {
      markColorPreviewDirty('feature-color-legend');
    }
    return addedCaption;
  };

  const updateLegendEntryColorByCaption = (caption, color) => {
    const existingEntry = findLegendEntryByCaption(caption);
    const resolvedCaption = existingEntry?.caption || caption;
    const overrideChanged = Boolean(
      existingEntry && !colorsMatch(legendColorOverrides[resolvedCaption], color)
    );
    if (overrideChanged) {
      legendColorOverrides[resolvedCaption] = color;
    }
    const updated = updateLegendEntryColorByCaptionRaw(resolvedCaption, color, { commit: false }) === true;
    if (updated) markColorPreviewDirty('feature-color-legend');
    return updated || overrideChanged;
  };

  const removeLegendEntry = (caption) => {
    const removed = removeLegendEntryRaw(caption, { commit: false }) === true;
    if (removed) markColorPreviewDirty('feature-color-legend');
    return removed;
  };

  const exactHashRulesForFeature = (feature) => manualSpecificRules.filter(
    (rule) => hashRuleTargetsFeatureExactly(rule, feature)
  );

  const liveFeatureColorMatches = (feature, color) => {
    const svg = getCurrentSvg();
    if (!svg) return true;
    const featureId = String(
      feature?.rendered_svg_id
      || feature?.renderedSvgId
      || feature?.rendered_feature_svg_id
      || feature?.renderedFeatureSvgId
      || feature?.svg_id
      || ''
    ).trim();
    if (!featureId) return false;
    const fillElements = getFeatureFillElements(svg, featureId);
    return fillElements.length > 0 && fillElements.every(
      (element) => colorsMatch(element.getAttribute('fill'), color)
    );
  };

  const featureColorAssignmentMatches = (
    feature,
    color,
    caption,
    { requireLegend = true } = {}
  ) => {
    const override = getFeatureOverride(featureColorOverrides, feature);
    if (!override || !colorsMatch(override.color, color) || !captionsMatch(override.caption, caption)) {
      return false;
    }
    const matchingRules = exactHashRulesForFeature(feature);
    if (!matchingRules.some(
      (rule) => colorsMatch(rule.color, color) && captionsMatch(rule.cap, caption)
    )) {
      return false;
    }
    if (!liveFeatureColorMatches(feature, color)) return false;
    if (!requireLegend) return true;
    const legendEntry = findLegendEntryByCaption(caption);
    return Boolean(legendEntry && colorsMatch(legendEntry.color, color));
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
      const element = getFeatureFillElements(svg, feat.svg_id)[0] || null;
      const fill = element?.getAttribute('fill');
      if (fill) {
        return resolveColorToHex(fill) || fill;
      }
    }

    const overrideColor = getFeatureOverride(featureColorOverrides, feat)?.color;
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
    const qualifier = getFeatureQualifier(feat);
    if (!qualifier?.val) return;
    const existingIdx = manualSpecificRules.findIndex(
      (rule) => hashRuleTargetsFeatureExactly(rule, feat)
    );

    if (existingIdx >= 0) {
      manualSpecificRules[existingIdx].feat = feat.type;
      manualSpecificRules[existingIdx].qual = 'hash';
      manualSpecificRules[existingIdx].val = qualifier.val;
      manualSpecificRules[existingIdx].color = color;
      manualSpecificRules[existingIdx].cap = caption;
      delete manualSpecificRules[existingIdx].fromFile;
      return;
    }

    const nextRule = {
      feat: feat.type,
      qual: 'hash',
      val: qualifier.val,
      color,
      cap: caption
    };
    const firstConflictingHashIdx = manualSpecificRules.findIndex(
      (rule) => rule?.feat === feat.type && isHashSpecificRule(rule) && ruleMatchesFeature(feat, rule)
    );
    if (firstConflictingHashIdx >= 0) {
      manualSpecificRules.splice(firstConflictingHashIdx, 0, nextRule);
    } else {
      manualSpecificRules.push(nextRule);
    }
  };

  const removeFeatureHashRules = (feature) => {
    for (let index = manualSpecificRules.length - 1; index >= 0; index -= 1) {
      if (hashRuleTargetsFeatureExactly(manualSpecificRules[index], feature)) {
        manualSpecificRules.splice(index, 1);
      }
    }
  };

  const upsertLabelSpecificRule = (rule, color, caption) => {
    if (!rule?.feat || !rule?.qual || !rule?.val) return false;
    const qualifier = String(rule.qual).toLowerCase();
    const existingIdx = manualSpecificRules.findIndex(
      (candidate) =>
        candidate.feat === rule.feat &&
        String(candidate.qual || '').toLowerCase() === qualifier &&
        candidate.val === rule.val
    );
    if (existingIdx >= 0) {
      const [existing] = manualSpecificRules.splice(existingIdx, 1);
      Object.assign(existing, { ...rule, color, cap: caption });
      delete existing.fromFile;
      const firstSameQualifierIdx = manualSpecificRules.findIndex(
        (candidate) =>
          candidate.feat === rule.feat && String(candidate.qual || '').toLowerCase() === qualifier
      );
      manualSpecificRules.splice(
        firstSameQualifierIdx >= 0 ? firstSameQualifierIdx : manualSpecificRules.length,
        0,
        existing
      );
      return true;
    }

    const firstSameQualifierIdx = manualSpecificRules.findIndex(
      (candidate) =>
        candidate.feat === rule.feat && String(candidate.qual || '').toLowerCase() === qualifier
    );
    const nextRule = { ...rule, color, cap: caption };
    if (firstSameQualifierIdx >= 0) {
      manualSpecificRules.splice(firstSameQualifierIdx, 0, nextRule);
    } else {
      manualSpecificRules.push(nextRule);
    }
    return true;
  };

  const featureIdentityKey = (feature) => String(
    feature?.id || `${feature?.type || ''}:${feature?.svg_id || ''}`
  );

  const getSafeLabelSpecificRule = (features, label) => {
    if (typeof getLabelSpecificRule !== 'function') return null;
    const candidates = features.map((feature) => getLabelSpecificRule(feature, label));
    if (candidates.some((rule) => !rule)) return null;

    const first = candidates[0];
    const sameSelector = candidates.every(
      (rule) =>
        rule.feat === first.feat &&
        String(rule.qual || '').toLowerCase() === String(first.qual || '').toLowerCase() &&
        rule.val === first.val
    );
    if (!sameSelector) return null;

    const selectedKeys = new Set(features.map(featureIdentityKey));
    const matchedFeatures = extractedFeatures.value.filter((feature) => ruleMatchesFeature(feature, first));
    if (
      matchedFeatures.length !== selectedKeys.size ||
      matchedFeatures.some((feature) => !selectedKeys.has(featureIdentityKey(feature)))
    ) {
      return null;
    }

    const safetyFeatures = Array.isArray(biologicalFeatures?.value) && biologicalFeatures.value.length > 0
      ? biologicalFeatures.value
      : extractedFeatures.value;
    if (safetyFeatures.filter((feature) => ruleMatchesFeature(feature, first)).length !== selectedKeys.size) {
      return null;
    }

    const hasPrecedenceConflict = manualSpecificRules.some((existing) => {
      const matchingSelected = features.filter((feature) => ruleMatchesFeature(feature, existing));
      if (matchingSelected.length === 0) return false;
      if (isHashSpecificRule(existing)) {
        return matchingSelected.some((feature) => !hashRuleTargetsFeatureExactly(existing, feature));
      }
      if (existing.feat !== first.feat) return false;
      const existingQualifier = String(existing.qual || '').toLowerCase();
      const candidateQualifier = String(first.qual || '').toLowerCase();
      return existingQualifier !== candidateQualifier;
    });
    return hasPrecedenceConflict ? null : first;
  };

  const syncFeatureLegendOverrides = (features, caption, color) => {
    for (const feature of features) {
      upsertFeatureHashRule(feature, color, caption);
      featureColorOverrides[featureOverrideKey(feature)] = { color, caption };
      updateClickedFeatureLegendState(feature, caption, color);
      applyFeatureColorPreview(feature, color);
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
    onLegendGeometryChanged();
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

  const applyColorToFeatureGroup = async (
    features,
    targetCaption,
    color,
    { preferLabelRules = false } = {}
  ) => {
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

    const labelRule = preferLabelRules
      ? getSafeLabelSpecificRule(features, normalizedTargetCaption)
      : null;
    if (labelRule) {
      for (const feature of features) {
        removeFeatureHashRules(feature);
      }
      upsertLabelSpecificRule(labelRule, color, finalCaption);
    } else {
      for (const feature of features) {
        upsertFeatureHashRule(feature, color, finalCaption);
      }
    }
    for (const feature of features) {
      featureColorOverrides[featureOverrideKey(feature)] = { color, caption: finalCaption };
      updateClickedFeatureLegendState(feature, finalCaption, color);
    }

    applySpecificRulesToSvg();
    await reclaimOrphanedBaseCaptions();
    extractLegendEntries();
  };

  const applyColorToLegendSpecificRules = async (targetCaption, color, captionFeatures = []) => {
    const normalizedTargetCaption = normalizeCaption(targetCaption);
    if (!normalizedTargetCaption) return false;
    const normalizedColor = String(color || '').trim();
    if (!normalizedColor) return false;

    const specificRules = manualSpecificRules.filter(
      (rule) => !isHashSpecificRule(rule) && captionsMatch(rule.cap, normalizedTargetCaption)
    );
    if (specificRules.length === 0) return false;

    const existingLegendEntry = findLegendEntryByCaption(normalizedTargetCaption);
    const finalCaption = existingLegendEntry?.caption || normalizeCaption(specificRules[0].cap) || normalizedTargetCaption;
    const coveredFeatures = extractedFeatures.value.filter((feature) =>
      specificRules.some((rule) => ruleMatchesFeature(feature, rule))
    );
    const affectedFeatures = new Map();
    [...captionFeatures, ...coveredFeatures].forEach((feature) => {
      const key = String(feature?.svg_id || feature?.id || '').trim();
      if (key) affectedFeatures.set(key, feature);
    });
    const captionHashRules = manualSpecificRules.filter(
      (rule) => isHashSpecificRule(rule) && captionsMatch(rule.cap, finalCaption)
    );
    const removesCoveredHashRule = captionHashRules.some((rule) => (
      coveredFeatures.some((feature) => hashRuleTargetsFeatureExactly(rule, feature))
    ));
    const rulesAlreadyMatch = specificRules.every(
      (rule) => colorsMatch(rule.color, normalizedColor) && captionsMatch(rule.cap, finalCaption)
    ) && captionHashRules.every(
      (rule) => colorsMatch(rule.color, normalizedColor) && captionsMatch(rule.cap, finalCaption)
    );
    const overridesAlreadyMatch = Array.from(affectedFeatures.values()).every((feature) => {
      if (!feature?.id) return true;
      const override = getFeatureOverride(featureColorOverrides, feature);
      return Boolean(
        override
        && colorsMatch(override.color, normalizedColor)
        && captionsMatch(override.caption, finalCaption)
      );
    });
    const liveColorsAlreadyMatch = Array.from(affectedFeatures.values()).every(
      (feature) => liveFeatureColorMatches(feature, normalizedColor)
    );
    const legendAlreadyMatches = !existingLegendEntry || colorsMatch(existingLegendEntry.color, normalizedColor);
    if (
      !removesCoveredHashRule
      && rulesAlreadyMatch
      && overridesAlreadyMatch
      && liveColorsAlreadyMatch
      && legendAlreadyMatches
    ) {
      return false;
    }

    for (const rule of specificRules) {
      rule.color = normalizedColor;
      rule.cap = finalCaption;
    }
    updateLegendEntryColorByCaption(finalCaption, normalizedColor);

    for (let i = manualSpecificRules.length - 1; i >= 0; i--) {
      const rule = manualSpecificRules[i];
      if (!isHashSpecificRule(rule) || !captionsMatch(rule.cap, finalCaption)) continue;
      if (coveredFeatures.some((feature) => hashRuleTargetsFeatureExactly(rule, feature))) {
        manualSpecificRules.splice(i, 1);
      } else {
        rule.color = normalizedColor;
        rule.cap = finalCaption;
      }
    }

    affectedFeatures.forEach((feature) => {
      if (!feature?.id) return;
      featureColorOverrides[featureOverrideKey(feature)] = { color: normalizedColor, caption: finalCaption };
      updateClickedFeatureLegendState(feature, finalCaption, normalizedColor);
    });

    applySpecificRulesToSvg();
    await nextTick();
    extractLegendEntries();
    refreshLegendEntryFeatureIds([finalCaption]);
    return true;
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
      if (!(await applyColorToLegendSpecificRules(targetLegendName, color, allFeatures))) {
        await applyColorToFeatureGroup(allFeatures, targetLegendName, color);
      }
    } else if (choice === 'displayLabel') {
      const displayLabel =
        normalizeCaption(colorScopeDialog.displayLabel) || normalizeCaption(getDisplayedFeatureLabel(feat));
      if (!displayLabel) {
        clearColorScopeDialog();
        return;
      }
      const displaySiblings = findFeaturesWithSameDisplayedLabel(feat, displayLabel);
      const allFeatures = [feat, ...displaySiblings];
      await applyColorToFeatureGroup(allFeatures, displayLabel, color, { preferLabelRules: true });
    } else if (choice === 'single') {
      let singleCaption = legendName;
      if (matchingRule && colorScopeDialog.ruleMatchCount > 1) {
        const ruleCaption = matchingRule.cap || matchingRule.val;
        if (legendName === ruleCaption) {
          singleCaption = getIndividualFeatureLabel(feat);
        }
      }
      await setFeatureColor(feat, color, singleCaption);
    } else if (choice === 'annotationLabel') {
      const annotationLabel =
        normalizeCaption(colorScopeDialog.annotationLabel) || normalizeCaption(getIndividualFeatureLabel(feat));
      if (!annotationLabel) {
        clearColorScopeDialog();
        return;
      }
      const annotationSiblings = findFeaturesWithSameIndividualLabel(feat, annotationLabel);
      const allFeatures = [feat, ...annotationSiblings];
      await applyColorToFeatureGroup(allFeatures, annotationLabel, color, { preferLabelRules: true });
    } else if (choice === 'useExisting') {
      if (existingCaptionColor) {
        const targetLegendName = normalizeCaption(legendName) || normalizeCaption(getEffectiveLegendCaption(feat));
        upsertFeatureHashRule(feat, existingCaptionColor, targetLegendName);
        if (clickedFeature.value && clickedFeature.value.svg_id === feat.svg_id) {
          clickedFeature.value.color = existingCaptionColor;
          clickedFeature.value.legendName = targetLegendName;
          clickedFeature.value.appliedLegendName = targetLegendName;
        }
        featureColorOverrides[featureOverrideKey(feat)] = {
          color: existingCaptionColor,
          caption: targetLegendName
        };
        applyFeatureColorPreview(feat, existingCaptionColor);
        applySpecificRulesToSvg();
      }
    }

    clearColorScopeDialog();
  };

  const updateClickedFeatureStroke = (strokeColor, strokeWidth) => {
    if (!clickedFeature.value) return false;
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const svgId = clickedFeature.value.svg_id;
    const elements = getFeatureElements(svg, svgId);
    if (elements.length === 0) return false;
    const normalizedStrokeColor = strokeColor === null || strokeColor === undefined
      ? null
      : String(strokeColor).trim();
    const normalizedStrokeWidth = normalizeStrokeWidthValue(strokeWidth);
    if (normalizedStrokeColor === null && normalizedStrokeWidth === null) return false;
    const firstElement = elements[0] || null;
    let changed = false;

    elements.forEach((element) => {
      if (normalizedStrokeColor !== null && !strokeColorAttributeMatches(element, normalizedStrokeColor)) {
        element.setAttribute('stroke', normalizedStrokeColor);
        changed = true;
      }
      if (normalizedStrokeWidth !== null && !strokeWidthAttributeMatches(element, normalizedStrokeWidth)) {
        element.setAttribute('stroke-width', normalizedStrokeWidth);
        changed = true;
      }
    });

    if (!changed) return false;
    recordFeatureStrokeOverride(clickedFeature.value.feat || clickedFeature.value, {
      strokeColor: normalizedStrokeColor,
      strokeWidth: normalizedStrokeWidth,
      originalStrokeColor: clickedFeature.value.originalStrokeColor ?? null,
      originalStrokeWidth: clickedFeature.value.originalStrokeWidth ?? firstElement?.getAttribute('stroke-width') ?? null
    });

    if (normalizedStrokeColor !== null) clickedFeature.value.strokeColor = normalizedStrokeColor;
    if (normalizedStrokeWidth !== null) clickedFeature.value.strokeWidth = normalizedStrokeWidth;

    persistCurrentSvg(svg, 'feature-stroke');
    return true;
  };

  const resetClickedFeatureStroke = () => {
    if (!clickedFeature.value) return false;
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const svgId = clickedFeature.value.svg_id;
    const elements = getFeatureElements(svg, svgId);

    const originalColor = originalSvgStroke.value.color;
    const originalWidth = normalizeStrokeWidthValue(originalSvgStroke.value.width);
    let changed = false;

    elements.forEach((element) => {
      if (!strokeColorAttributeMatches(element, originalColor)) {
        if (originalColor === null) element.removeAttribute('stroke');
        else element.setAttribute('stroke', originalColor);
        changed = true;
      }
      if (!strokeWidthAttributeMatches(element, originalWidth)) {
        if (originalWidth === null) element.removeAttribute('stroke-width');
        else element.setAttribute('stroke-width', originalWidth);
        changed = true;
      }
    });

    const feature = clickedFeature.value.feat || clickedFeature.value;
    const overrideKey = featureStrokeKey(feature, svgId);
    const hadOverride = Boolean(overrideKey && featureStrokeOverrides[overrideKey]);
    if (!changed && !hadOverride) return false;
    clickedFeature.value.strokeColor = originalColor || '';
    clickedFeature.value.strokeWidth = originalWidth ?? '';
    clearFeatureStrokeOverride(feature, svgId);

    if (changed) persistCurrentSvg(svg, 'feature-stroke');
    return true;
  };

  const getFeatureStrokeColorValue = (featureLike) => {
    const feature = featureLike?.feat || featureLike;
    const override = getFeatureOverride(featureStrokeOverrides, feature);
    return override && hasOwn(override, 'strokeColor')
      ? override.strokeColor
      : null;
  };

  const setClickedFeatureStrokeColorValue = (value) => {
    if (value !== null) {
      return updateClickedFeatureStroke(String(value || '').trim(), null);
    }
    if (!clickedFeature.value || !svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;
    const feature = clickedFeature.value.feat || clickedFeature.value;
    const key = featureStrokeKey(feature, clickedFeature.value.svg_id);
    const override = key ? featureStrokeOverrides[key] : null;
    const inheritedColor = override && hasOwn(override, 'originalStrokeColor')
      ? override.originalStrokeColor
      : (clickedFeature.value.originalStrokeColor ?? originalSvgStroke.value.color);
    const elements = getFeatureElements(svg, clickedFeature.value.svg_id);
    const domChanged = elements.some((element) => !strokeColorAttributeMatches(element, inheritedColor));
    const stateChanged = Boolean(override && hasOwn(override, 'strokeColor'));
    if (!domChanged && !stateChanged) return false;
    if (override) {
      delete override.strokeColor;
      if (!hasOwn(override, 'strokeWidth')) delete featureStrokeOverrides[key];
    }
    elements.forEach((element) => {
      if (strokeColorAttributeMatches(element, inheritedColor)) return;
      if (inheritedColor === null || inheritedColor === '') {
        element.removeAttribute('stroke');
      } else {
        element.setAttribute('stroke', inheritedColor);
      }
    });
    clickedFeature.value.strokeColor = inheritedColor || '';
    if (domChanged) persistCurrentSvg(svg, 'feature-stroke');
    return true;
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
      const elements = getFeatureFillElements(svg, svgId);
      elements.forEach((el) => {
        el.setAttribute('fill', defaultColor);
      });

      clickedFeature.value.color = defaultColor;

      removeFeatureHashRules(feat);

      let resetCaption = feat.type === 'CDS' ? 'other proteins' : `other ${feat.type}s`;
      if (choice === 'this_with_legend') {
        console.log(`Attempting to add legend entry: caption="${caption}", color="${defaultColor}"`);
        const addedCaption = await addLegendEntry(caption, defaultColor);
        console.log(`addLegendEntry returned: ${addedCaption}`);
        if (addedCaption) {
          resetCaption = addedCaption;
          extractLegendEntries();
          console.log(`Added legend entry: ${addedCaption} with color: ${defaultColor}`);
        } else {
          console.error(`Failed to add legend entry for caption="${caption}"`);
        }
      }

      const remainsCoveredBySpecificRule = manualSpecificRules.some(
        (rule) => ruleMatchesFeature(feat, rule)
      );
      if (remainsCoveredBySpecificRule || choice === 'this_with_legend') {
        upsertFeatureHashRule(feat, defaultColor, resetCaption);
        featureColorOverrides[featureOverrideKey(feat)] = {
          color: defaultColor,
          caption: resetCaption
        };
      } else {
        delete featureColorOverrides[featureOverrideKey(feat)];
      }
      applySpecificRulesToSvg();

      console.log(`Reset fill color to default (${defaultColor}) for feature: ${svgId}`);
    } else if (choice === 'all') {
      const matchingFeatures = extractedFeatures.value.filter((f) => getFeatureCaption(f) === caption);

      for (const matchFeat of matchingFeatures) {
        const elements = getFeatureFillElements(svg, matchFeat.svg_id);
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

    persistCurrentSvg(svg);

    clickedFeature.value = null;
  };

  const applyStrokeToAllSiblings = async () => {
    if (!clickedFeature.value) return false;
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const currentStrokeColor = String(clickedFeature.value.strokeColor || '').trim();
    const currentStrokeWidth = normalizeStrokeWidthValue(clickedFeature.value.strokeWidth);
    if (!currentStrokeColor && currentStrokeWidth === null) return false;
    const currentFillColor = clickedFeature.value.color;
    const clickedSvgId = clickedFeature.value.svg_id;

    const feat = clickedFeature.value.feat;
    if (!feat) {
      console.warn('No feature data found in clickedFeature');
      return false;
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
        const el = getFeatureFillElements(svg, f.svg_id)[0] || null;
        if (!el) return false;
        const fillColor = el.getAttribute('fill');
        return fillColor && fillColor.toLowerCase() === currentFillColor?.toLowerCase();
      });
      siblingFeatureIds = siblingFeatures.map((f) => f.svg_id);
      console.log(
        `Fallback: Found ${siblingFeatureIds.length} features by caption="${caption}" and color="${currentFillColor}"`
      );
    }

    let domChanged = false;
    let stateChanged = false;
    for (const svgId of siblingFeatureIds) {
      const elements = getFeatureElements(svg, svgId);
      const featureNeedsUpdate = elements.some((element) => (
        (currentStrokeColor && !strokeColorAttributeMatches(element, currentStrokeColor))
        || (currentStrokeWidth !== null && !strokeWidthAttributeMatches(element, currentStrokeWidth))
      ));
      if (!featureNeedsUpdate) continue;
      const firstElement = elements[0] || null;
      const siblingFeature = findFeatureBySvgId(svgId) || { id: svgId, svg_id: svgId };
      recordFeatureStrokeOverride(siblingFeature, {
        strokeColor: currentStrokeColor,
        strokeWidth: currentStrokeWidth,
        originalStrokeColor: firstElement?.getAttribute('stroke') ?? null,
        originalStrokeWidth: firstElement?.getAttribute('stroke-width') ?? null
      });
      elements.forEach((el) => {
        if (currentStrokeColor && !strokeColorAttributeMatches(el, currentStrokeColor)) {
          el.setAttribute('stroke', currentStrokeColor);
          domChanged = true;
        }
        if (currentStrokeWidth !== null && !strokeWidthAttributeMatches(el, currentStrokeWidth)) {
          el.setAttribute('stroke-width', currentStrokeWidth);
          domChanged = true;
        }
      });
      stateChanged = true;
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
                if (!strokeColorAttributeMatches(path, currentStrokeColor)) {
                  path.setAttribute('stroke', currentStrokeColor);
                  domChanged = true;
                }
              }
              if (currentStrokeWidth !== null && !strokeWidthAttributeMatches(path, currentStrokeWidth)) {
                path.setAttribute('stroke-width', currentStrokeWidth);
                domChanged = true;
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
        domChanged = true;
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
        stateChanged = true;
      }
    }

    if (!domChanged && !stateChanged) return false;

    if (!legendStrokeOverrides[overrideKey]) {
      legendStrokeOverrides[overrideKey] = {};
      stateChanged = true;
    }
    if (currentStrokeColor && !colorsMatch(legendStrokeOverrides[overrideKey].strokeColor, currentStrokeColor)) {
      legendStrokeOverrides[overrideKey].strokeColor = currentStrokeColor;
      stateChanged = true;
    }
    if (
      currentStrokeWidth !== null
      && normalizeStrokeWidthValue(legendStrokeOverrides[overrideKey].strokeWidth) !== currentStrokeWidth
    ) {
      legendStrokeOverrides[overrideKey].strokeWidth = currentStrokeWidth;
      stateChanged = true;
    }

    if (domChanged) persistCurrentSvg(svg, 'feature-stroke');

    console.log(
      `Applied stroke (color: ${currentStrokeColor}, width: ${currentStrokeWidth}) to ${siblingFeatureIds.length} features`
    );
    return domChanged || stateChanged;
  };

  const uniqueFeaturesBySvgId = (features) => {
    const seen = new Set();
    return (Array.isArray(features) ? features : []).filter((feature) => {
      const svgId = String(feature?.svg_id || '').trim();
      if (!svgId || seen.has(svgId)) return false;
      seen.add(svgId);
      return true;
    });
  };

  const applyColorToSelectedFeatures = async (features, color, caption) => {
    const targetFeatures = uniqueFeaturesBySvgId(features);
    const targetColor = resolveColorToHex(color) || String(color || '').trim();
    const targetCaption = normalizeCaption(caption);
    if (targetFeatures.length === 0 || !targetColor || !targetCaption) return false;
    await applyColorToFeatureGroup(targetFeatures, targetCaption, targetColor);
    return true;
  };

  const applyStrokeToSelectedFeatures = (features, strokeColor, strokeWidth) => {
    const targetFeatures = uniqueFeaturesBySvgId(features);
    if (targetFeatures.length === 0 || !svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const normalizedStrokeColor = String(strokeColor || '').trim();
    const normalizedStrokeWidth = normalizeStrokeWidthValue(strokeWidth);
    if (!normalizedStrokeColor && normalizedStrokeWidth === null) return false;

    let updatedCount = 0;
    targetFeatures.forEach((feature) => {
      const elements = getFeatureElements(svg, feature.svg_id);
      if (elements.length === 0) return;
      const needsUpdate = elements.some((element) => (
        (normalizedStrokeColor && !strokeColorAttributeMatches(element, normalizedStrokeColor))
        || (normalizedStrokeWidth !== null && !strokeWidthAttributeMatches(element, normalizedStrokeWidth))
      ));
      if (!needsUpdate) return;
      const firstElement = elements[0] || null;
      recordFeatureStrokeOverride(feature, {
        strokeColor: normalizedStrokeColor || null,
        strokeWidth: normalizedStrokeWidth,
        originalStrokeColor: firstElement?.getAttribute('stroke') ?? null,
        originalStrokeWidth: firstElement?.getAttribute('stroke-width') ?? null
      });
      elements.forEach((element) => {
        let changed = false;
        if (normalizedStrokeColor && !strokeColorAttributeMatches(element, normalizedStrokeColor)) {
          element.setAttribute('stroke', normalizedStrokeColor);
          changed = true;
        }
        if (normalizedStrokeWidth !== null && !strokeWidthAttributeMatches(element, normalizedStrokeWidth)) {
          element.setAttribute('stroke-width', normalizedStrokeWidth);
          changed = true;
        }
        if (changed) updatedCount += 1;
      });
      if (clickedFeature.value?.svg_id === feature.svg_id) {
        if (normalizedStrokeColor) clickedFeature.value.strokeColor = normalizedStrokeColor;
        if (normalizedStrokeWidth !== null) clickedFeature.value.strokeWidth = normalizedStrokeWidth;
      }
    });

    if (updatedCount > 0) {
      persistCurrentSvg(svg, 'feature-stroke');
    }
    return updatedCount > 0;
  };

  const setFeatureColor = async (feat, color, customCaption = null) => {
    if (!feat) return false;
    const qualInfo = getFeatureQualifier(feat);
    if (!qualInfo) {
      console.warn(
        `Cannot identify feature: ${feat.type} at ${feat.start}..${feat.end} (no locus_tag, gene, or product)`
      );
      return false;
    }
    const normalizedColor = String(color || '').trim();
    if (!normalizedColor) return false;
    const featureKey = featureOverrideKey(feat);

    const caption = normalizeCaption(
      customCaption || feat.product || feat.gene || feat.locus_tag || `${feat.type} at ${feat.start}..${feat.end}`
    );
    if (!caption) return false;
    if (featureColorAssignmentMatches(feat, normalizedColor, caption)) return false;

    const oldOverride = featureColorOverrides[featureKey];
    const oldCaption = normalizeCaption(oldOverride?.caption);
    featureColorOverrides[featureKey] = { color: normalizedColor, caption };

    applyFeatureColorPreview(feat, normalizedColor);

    await nextTick();
    let actualCaption = caption;
    if (caption) {
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
            actualCaption = await addLegendEntry(caption, normalizedColor);
            if (actualCaption && typeof actualCaption === 'string') {
              addedLegendCaptions.value.add(actualCaption);
            }
          } else {
            updateLegendEntryColorByCaption(oldCaption, normalizedColor);
          }
        } else {
          removeLegendEntry(oldCaption);
          actualCaption = await addLegendEntry(caption, normalizedColor);
          if (actualCaption && typeof actualCaption === 'string') {
            addedLegendCaptions.value.add(actualCaption);
          }
        }
      } else {
        actualCaption = await addLegendEntry(caption, normalizedColor);
        if (actualCaption && typeof actualCaption === 'string') {
          addedLegendCaptions.value.add(actualCaption);
        }
      }

      if (actualCaption && typeof actualCaption === 'string' && actualCaption !== caption) {
        featureColorOverrides[featureKey] = { color: normalizedColor, caption: actualCaption };
      }
    }

    skipExtractOnSvgChange.value = true;

    const finalCaption = actualCaption && typeof actualCaption === 'string' ? actualCaption : caption;
    upsertFeatureHashRule(feat, normalizedColor, finalCaption);

    featureColorOverrides[featureKey] = { color: normalizedColor, caption: finalCaption };
    updateClickedFeatureLegendState(feat, finalCaption, normalizedColor);

    await reclaimOrphanedBaseCaptions();

    await nextTick();
    await nextTick();
    skipExtractOnSvgChange.value = false;
    extractLegendEntries();
    return true;
  };

  const setFeatureColorValue = async (feat, value, customCaption = null) => {
    if (!feat) return false;
    const featureKey = featureOverrideKey(feat);
    if (!featureKey) return false;
    if (value === null) {
      if (!getFeatureOverride(featureColorOverrides, feat) && exactHashRulesForFeature(feat).length === 0) {
        return false;
      }
      delete featureColorOverrides[featureKey];
      removeFeatureHashRules(feat);
      const inheritedColor = appliedPaletteColors.value[feat.type] || '#cccccc';
      applyFeatureColorPreview(feat, inheritedColor);
      applySpecificRulesToSvg();
      if (clickedFeature.value?.feat === feat || clickedFeature.value?.svg_id === feat.svg_id) {
        clickedFeature.value.color = inheritedColor;
      }
      return true;
    }
    const normalizedValue = String(value || '').trim();
    if (normalizedValue.toLowerCase() === 'none') {
      const caption = normalizeCaption(
        customCaption
        || getEffectiveLegendCaption(feat)
        || getFeatureCaption(feat)
        || feat.type
      );
      if (featureColorAssignmentMatches(feat, 'none', caption, { requireLegend: false })) {
        return false;
      }
      featureColorOverrides[featureKey] = { color: 'none', caption };
      upsertFeatureHashRule(feat, 'none', caption);
      applyFeatureColorPreview(feat, 'none');
      applySpecificRulesToSvg();
      updateClickedFeatureLegendState(feat, caption, 'none');
      return true;
    }
    return setFeatureColor(feat, normalizedValue, customCaption);
  };

  return {
    applyStrokeToAllSiblings: colorAction(applyStrokeToAllSiblings),
    handleColorScopeChoice: colorAction(handleColorScopeChoice),
    handleLegendNameCommit: colorAction(handleLegendNameCommit),
    handleLegendRenameChoice: colorAction(handleLegendRenameChoice),
    renameLegendEntry: colorAction(renameLegendEntry),
    requestFeatureColorChange: colorAction(requestFeatureColorChange),
    selectLegendNameOption: colorAction(selectLegendNameOption),
    handleResetColorChoice: colorAction(handleResetColorChoice),
    applyColorToSelectedFeatures: colorAction(applyColorToSelectedFeatures),
    applyStrokeToSelectedFeatures: colorAction(applyStrokeToSelectedFeatures),
    resetClickedFeatureFillColor: colorAction(resetClickedFeatureFillColor),
    resetClickedFeatureStroke: colorAction(resetClickedFeatureStroke),
    getFeatureStrokeColorValue,
    setClickedFeatureStrokeColorValue: colorAction(setClickedFeatureStrokeColorValue),
    setFeatureColor: colorAction(setFeatureColor),
    setFeatureColorValue: colorAction(setFeatureColorValue),
    updateClickedFeatureColor: colorAction(updateClickedFeatureColor),
    updateClickedFeatureStroke: colorAction(updateClickedFeatureStroke)
  };
};
