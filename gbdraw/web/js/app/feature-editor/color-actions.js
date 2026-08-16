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
  if (
    typeof previewRuntime?.commitDomEdit !== 'function'
    || typeof previewRuntime?.runDomEdit !== 'function'
    || typeof previewRuntime?.applyFeatureFillChanges !== 'function'
    || typeof previewRuntime?.applyFeatureStrokeChanges !== 'function'
  ) {
    throw new Error('createFeatureColorActions requires the preview runtime edit protocol.');
  }
  const {
    appliedPaletteColors,
    manualSpecificRules,
    extractedFeatures,
    biologicalFeatures,
    featureColorOverrides,
    svgContainer,
    clickedFeature,
    featureStyleScopeDialog,
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
    return updated;
  };

  const runColorAction = (action) => previewRuntime.runDomEdit({
    reason: 'feature-style-action',
    action
  });

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

  const clearFeatureStyleScopeDialog = () => {
    featureStyleScopeDialog.show = false;
    featureStyleScopeDialog.kind = 'fill';
    featureStyleScopeDialog.feat = null;
    featureStyleScopeDialog.color = null;
    featureStyleScopeDialog.strokeColor = null;
    featureStyleScopeDialog.strokeWidth = null;
    featureStyleScopeDialog.matchingRule = null;
    featureStyleScopeDialog.ruleMatchCount = 0;
    featureStyleScopeDialog.legendName = null;
    featureStyleScopeDialog.siblingCount = 0;
    featureStyleScopeDialog.displayLabel = null;
    featureStyleScopeDialog.displayLabelSiblingCount = 0;
    featureStyleScopeDialog.annotationLabel = null;
    featureStyleScopeDialog.annotationLabelSiblingCount = 0;
    featureStyleScopeDialog.existingCaptionRule = null;
    featureStyleScopeDialog.existingCaptionColor = null;
  };

  const getFeatureStyleScope = (feat, requestedLegendName = null) => {
    const requestedCaption = normalizeCaption(requestedLegendName);
    const effectiveCaption = normalizeCaption(getEffectiveLegendCaption(feat));
    const fallbackCaption = normalizeCaption(getFeatureCaption(feat));
    const legendName = requestedCaption || effectiveCaption || fallbackCaption;
    if (!legendName) return null;

    const matchingRule = findMatchingRegexRule(feat);
    const ruleMatchCount = matchingRule ? countFeaturesMatchingRule(matchingRule) : 0;
    const siblingCount = findFeaturesWithSameLegendItem(feat, legendName).length;
    const displayLabel = normalizeCaption(getDisplayedFeatureLabel(feat));
    const displayLabelSiblingCount = displayLabel
      ? findFeaturesWithSameDisplayedLabel(feat, displayLabel).length
      : 0;
    const annotationLabel = normalizeCaption(getIndividualFeatureLabel(feat));
    const annotationLabelSiblingCount = annotationLabel
      ? findFeaturesWithSameIndividualLabel(feat, annotationLabel).length
      : 0;

    return {
      requestedCaption,
      legendName,
      matchingRule,
      ruleMatchCount,
      siblingCount,
      displayLabel,
      displayLabelSiblingCount,
      annotationLabel,
      annotationLabelSiblingCount,
      needsDialog: Boolean(
        matchingRule || siblingCount > 0 || displayLabelSiblingCount > 0 || annotationLabelSiblingCount > 0
      )
    };
  };

  const openFeatureStyleScopeDialog = ({
    kind,
    feat,
    scope,
    color = null,
    strokeColor = null,
    strokeWidth = null,
    existingCaption = null,
    closePopup = false
  }) => {
    featureStyleScopeDialog.show = true;
    featureStyleScopeDialog.kind = kind;
    featureStyleScopeDialog.feat = feat;
    featureStyleScopeDialog.color = color;
    featureStyleScopeDialog.strokeColor = strokeColor;
    featureStyleScopeDialog.strokeWidth = strokeWidth;
    featureStyleScopeDialog.matchingRule = scope.matchingRule;
    featureStyleScopeDialog.ruleMatchCount = scope.ruleMatchCount;
    featureStyleScopeDialog.legendName = scope.legendName;
    featureStyleScopeDialog.siblingCount = scope.siblingCount;
    featureStyleScopeDialog.displayLabel = scope.displayLabel;
    featureStyleScopeDialog.displayLabelSiblingCount = scope.displayLabelSiblingCount;
    featureStyleScopeDialog.annotationLabel = scope.annotationLabel;
    featureStyleScopeDialog.annotationLabelSiblingCount = scope.annotationLabelSiblingCount;
    featureStyleScopeDialog.existingCaptionRule = existingCaption?.rule || null;
    featureStyleScopeDialog.existingCaptionColor = existingCaption?.color || null;
    if (closePopup) clickedFeature.value = null;
  };

  const getCurrentSvg = () => svgContainer.value?.querySelector('svg') || null;

  const addLegendEntry = async (caption, color, options = {}) => {
    return addLegendEntryRaw(caption, color, options);
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
    const updated = updateLegendEntryColorByCaptionRaw(resolvedCaption, color) === true;
    return updated || overrideChanged;
  };

  const removeLegendEntry = (caption) => {
    return removeLegendEntryRaw(caption) === true;
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
    const committed = previewRuntime.commitDomEdit({
      reason: 'legend-rename',
      invalidateIndexes: ['legend'],
      mutate: ({ mutation }) => {
        let updated = false;
        for (const targetGroup of getAllFeatureLegendGroups(svg)) {
          const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(oldCaption)}"]`);
          if (!entryGroup) continue;

          mutation.setAttribute(entryGroup, 'data-legend-key', newCaption);
          const textEl = entryGroup.querySelector('text');
          if (textEl) mutation.setTextContent(textEl, newCaption);

          if (color) {
            const paths = entryGroup.querySelectorAll('path');
            for (const path of paths) {
              const fill = path.getAttribute('fill');
              if (fill && fill !== 'none' && !fill.startsWith('url(')) {
                mutation.setAttribute(path, 'fill', color);
                break;
              }
            }
          }
          updated = true;
        }
        if (updated) compactLegendEntries(svg);
        return updated;
      }
    });
    if (!committed.changed) return false;

    moveCaptionStateKey(legendColorOverrides, oldCaption, newCaption);
    moveCaptionStateKey(legendStrokeOverrides, oldCaption, newCaption);
    moveAddedLegendCaption(oldCaption, newCaption);
    syncOriginalLegendMetadataRename(oldCaption, newCaption, color);

    const legendEntry = legendEntries.value.find((entry) => captionsMatch(entry?.caption, oldCaption));
    if (legendEntry) {
      legendEntry.caption = newCaption;
      if (color) {
        legendEntry.color = color;
      }
    }

    onLegendGeometryChanged();
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

  const applyColorToLegendSpecificRules = (targetCaption, color, captionFeatures = []) => {
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
    extractLegendEntries();
    refreshLegendEntryFeatureIds([finalCaption]);
    return true;
  };

  const requestFeatureColorChange = async (feat, color, requestedLegendName = null, options = {}) => {
    if (!feat) return;
    const scope = getFeatureStyleScope(feat, requestedLegendName);
    if (!scope) return;

    if (scope.needsDialog) {
      openFeatureStyleScopeDialog({
        kind: 'fill',
        feat,
        scope,
        color,
        existingCaption: findExistingCaptionColor(feat, scope.legendName),
        closePopup: options.closePopupOnDialog
      });
      return;
    }

    if (clickedFeature.value && clickedFeature.value.svg_id === feat.svg_id) {
      clickedFeature.value.color = color;
      if (scope.requestedCaption) {
        clickedFeature.value.legendName = scope.requestedCaption;
      }
    }
    await setFeatureColor(feat, color, scope.legendName);
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
    const { feat, color, matchingRule, legendName, existingCaptionColor } = featureStyleScopeDialog;
    if (choice === 'cancel' || !feat || !color) {
      clearFeatureStyleScopeDialog();
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
        clearFeatureStyleScopeDialog();
        return;
      }
      const siblings = findFeaturesWithSameLegendItem(feat, targetLegendName);
      const allFeatures = [feat, ...siblings];
      if (!applyColorToLegendSpecificRules(targetLegendName, color, allFeatures)) {
        await applyColorToFeatureGroup(allFeatures, targetLegendName, color);
      }
    } else if (choice === 'displayLabel') {
      const displayLabel =
        normalizeCaption(featureStyleScopeDialog.displayLabel) || normalizeCaption(getDisplayedFeatureLabel(feat));
      if (!displayLabel) {
        clearFeatureStyleScopeDialog();
        return;
      }
      const displaySiblings = findFeaturesWithSameDisplayedLabel(feat, displayLabel);
      const allFeatures = [feat, ...displaySiblings];
      await applyColorToFeatureGroup(allFeatures, displayLabel, color, { preferLabelRules: true });
    } else if (choice === 'single') {
      let singleCaption = legendName;
      if (matchingRule && featureStyleScopeDialog.ruleMatchCount > 1) {
        const ruleCaption = matchingRule.cap || matchingRule.val;
        if (legendName === ruleCaption) {
          singleCaption = getIndividualFeatureLabel(feat);
        }
      }
      await setFeatureColor(feat, color, singleCaption);
    } else if (choice === 'annotationLabel') {
      const annotationLabel =
        normalizeCaption(featureStyleScopeDialog.annotationLabel) || normalizeCaption(getIndividualFeatureLabel(feat));
      if (!annotationLabel) {
        clearFeatureStyleScopeDialog();
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

    clearFeatureStyleScopeDialog();
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
    const changed = previewRuntime.applyFeatureStrokeChanges([{
      featureId: svgId,
      ...(normalizedStrokeColor !== null ? { strokeColor: normalizedStrokeColor } : {}),
      ...(normalizedStrokeWidth !== null ? { strokeWidth: normalizedStrokeWidth } : {})
    }], { reason: 'feature-stroke' });
    if (!changed) return false;
    recordFeatureStrokeOverride(clickedFeature.value.feat || clickedFeature.value, {
      strokeColor: normalizedStrokeColor,
      strokeWidth: normalizedStrokeWidth,
      originalStrokeColor: clickedFeature.value.originalStrokeColor ?? null,
      originalStrokeWidth: clickedFeature.value.originalStrokeWidth ?? firstElement?.getAttribute('stroke-width') ?? null
    });

    if (normalizedStrokeColor !== null) clickedFeature.value.strokeColor = normalizedStrokeColor;
    if (normalizedStrokeWidth !== null) clickedFeature.value.strokeWidth = normalizedStrokeWidth;
    return true;
  };

  const requestClickedFeatureStrokeChange = (strokeColor, strokeWidth) => {
    if (!clickedFeature.value) return false;
    const feat = clickedFeature.value.feat;
    if (!feat) return false;

    const normalizedStrokeColor = String(strokeColor || '').trim();
    const normalizedStrokeWidth = normalizeStrokeWidthValue(strokeWidth);
    if (!normalizedStrokeColor && normalizedStrokeWidth === null) return false;

    const scope = getFeatureStyleScope(feat, clickedFeature.value.legendName);
    if (!scope) return false;
    if (scope.needsDialog) {
      openFeatureStyleScopeDialog({
        kind: 'stroke',
        feat,
        scope,
        strokeColor: normalizedStrokeColor,
        strokeWidth: normalizedStrokeWidth,
        closePopup: true
      });
      return false;
    }

    return updateClickedFeatureStroke(normalizedStrokeColor, normalizedStrokeWidth);
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
    const changed = previewRuntime.applyFeatureStrokeChanges([{
      featureId: svgId,
      strokeColor: originalColor,
      strokeWidth: originalWidth
    }], { reason: 'feature-stroke-reset' });

    const feature = clickedFeature.value.feat || clickedFeature.value;
    const overrideKey = featureStrokeKey(feature, svgId);
    const hadOverride = Boolean(overrideKey && featureStrokeOverrides[overrideKey]);
    if (!changed && !hadOverride) return false;
    clickedFeature.value.strokeColor = originalColor || '';
    clickedFeature.value.strokeWidth = originalWidth ?? '';
    clearFeatureStrokeOverride(feature, svgId);

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
      if (!clickedFeature.value) return false;
      const feature = clickedFeature.value.feat || clickedFeature.value;
      const override = getFeatureOverride(featureStrokeOverrides, feature);
      if (
        !hasOwn(override, 'strokeColor')
        && colorsMatch(resolveColorToHex(value), resolveColorToHex(clickedFeature.value.strokeColor))
      ) {
        const normalizedValue = String(value || '').trim();
        if (updateClickedFeatureStroke(normalizedValue, null)) return true;
        const svg = getCurrentSvg();
        if (!svg) return false;
        const elements = getFeatureElements(svg, clickedFeature.value.svg_id);
        if (elements.length === 0) return false;
        recordFeatureStrokeOverride(feature, {
          strokeColor: normalizedValue,
          originalStrokeColor: clickedFeature.value.originalStrokeColor ?? elements[0]?.getAttribute('stroke') ?? null,
          originalStrokeWidth: clickedFeature.value.originalStrokeWidth ?? elements[0]?.getAttribute('stroke-width') ?? null
        });
        return true;
      }
      return requestClickedFeatureStrokeChange(value, clickedFeature.value.strokeWidth);
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
    if (domChanged) {
      previewRuntime.applyFeatureStrokeChanges([{
        featureId: clickedFeature.value.svg_id,
        strokeColor: inheritedColor || null
      }], { reason: 'feature-stroke-color-reset' });
    }
    clickedFeature.value.strokeColor = inheritedColor || '';
    return true;
  };

  const setClickedFeatureStrokeWidthValue = (value) => {
    if (!clickedFeature.value) return false;
    const normalizedStrokeWidth = normalizeStrokeWidthValue(value);
    const currentStrokeWidth = normalizeStrokeWidthValue(clickedFeature.value.strokeWidth);
    if (normalizedStrokeWidth === null || normalizedStrokeWidth === currentStrokeWidth) return false;
    return requestClickedFeatureStrokeChange(clickedFeature.value.strokeColor, normalizedStrokeWidth);
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
      applyFeatureColorPreview(feat, defaultColor);

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
      previewRuntime.applyFeatureFillChanges(
        matchingFeatures.map((matchFeat) => ({
          featureId: matchFeat.svg_id,
          color: defaultColor
        })),
        { reason: 'bulk-feature-fill-reset' }
      );

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

    clickedFeature.value = null;
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

    const changes = [];
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
      changes.push({
        featureId: feature.svg_id,
        ...(normalizedStrokeColor ? { strokeColor: normalizedStrokeColor } : {}),
        ...(normalizedStrokeWidth !== null ? { strokeWidth: normalizedStrokeWidth } : {})
      });
      if (clickedFeature.value?.svg_id === feature.svg_id) {
        if (normalizedStrokeColor) clickedFeature.value.strokeColor = normalizedStrokeColor;
        if (normalizedStrokeWidth !== null) clickedFeature.value.strokeWidth = normalizedStrokeWidth;
      }
    });
    return previewRuntime.applyFeatureStrokeChanges(changes, { reason: 'bulk-feature-stroke' });
  };

  const applyStrokeToLegendEntry = (caption, strokeColor, strokeWidth) => {
    const targetLegendEntry = findLegendEntryByCaption(caption);
    const svg = getCurrentSvg();
    if (!targetLegendEntry || !svg) return false;

    const normalizedStrokeColor = String(strokeColor || '').trim();
    const normalizedStrokeWidth = normalizeStrokeWidthValue(strokeWidth);
    if (!normalizedStrokeColor && normalizedStrokeWidth === null) return false;

    let originalSwatchStroke = null;
    const escapedCaption = globalThis.CSS?.escape
      ? globalThis.CSS.escape(targetLegendEntry.caption)
      : String(targetLegendEntry.caption).replace(/["\\]/g, '\\$&');
    const domChanged = previewRuntime.commitDomEdit({
      reason: 'legend-stroke',
      invalidateIndexes: ['legend'],
      mutate: ({ mutation }) => {
        let updated = false;
        for (const targetGroup of getAllFeatureLegendGroups(svg)) {
          const entryGroup = targetGroup.querySelector(
            `g[data-legend-key="${escapedCaption}"]`
          );
          if (!entryGroup) continue;
          const swatch = Array.from(entryGroup.querySelectorAll('path')).find((path) => {
            const fill = path.getAttribute('fill');
            return fill && fill !== 'none' && !fill.startsWith('url(');
          });
          if (!swatch) continue;
          if (!originalSwatchStroke) {
            originalSwatchStroke = {
              color: swatch.getAttribute('stroke'),
              width: normalizeStrokeWidthValue(swatch.getAttribute('stroke-width'))
            };
          }
          if (normalizedStrokeColor && !strokeColorAttributeMatches(swatch, normalizedStrokeColor)) {
            mutation.setAttribute(swatch, 'stroke', normalizedStrokeColor);
            updated = true;
          }
          if (normalizedStrokeWidth !== null && !strokeWidthAttributeMatches(swatch, normalizedStrokeWidth)) {
            mutation.setAttribute(swatch, 'stroke-width', normalizedStrokeWidth);
            updated = true;
          }
        }
        return updated;
      }
    }).changed;

    const overrideKey = targetLegendEntry.caption;
    const previousOverride = legendStrokeOverrides[overrideKey] || {};
    const nextOverride = { ...previousOverride };
    if (
      !hasOwn(previousOverride, 'strokeColor')
      && !hasOwn(previousOverride, 'strokeWidth')
    ) {
      nextOverride.originalStrokeColor = originalSwatchStroke?.color ?? null;
      nextOverride.originalStrokeWidth = originalSwatchStroke?.width ?? null;
    }
    if (normalizedStrokeColor) nextOverride.strokeColor = normalizedStrokeColor;
    if (normalizedStrokeWidth !== null) nextOverride.strokeWidth = normalizedStrokeWidth;
    const stateChanged =
      !colorsMatch(previousOverride.strokeColor, nextOverride.strokeColor)
      || normalizeStrokeWidthValue(previousOverride.strokeWidth)
        !== normalizeStrokeWidthValue(nextOverride.strokeWidth);
    if (stateChanged) legendStrokeOverrides[overrideKey] = nextOverride;
    return domChanged || stateChanged;
  };

  const handleStrokeScopeChoice = (choice) => {
    const {
      feat,
      strokeColor,
      strokeWidth,
      matchingRule,
      legendName,
      displayLabel,
      annotationLabel
    } = featureStyleScopeDialog;
    if (choice === 'cancel' || !feat || (!strokeColor && strokeWidth === null)) {
      clearFeatureStyleScopeDialog();
      return false;
    }

    let targetFeatures = [];
    let legendCaption = '';
    if (choice === 'rule' && matchingRule) {
      targetFeatures = extractedFeatures.value.filter((candidate) => (
        candidate?.type === matchingRule.feat && ruleMatchesFeature(candidate, matchingRule)
      ));
      legendCaption = normalizeCaption(matchingRule.cap);
    } else if (choice === 'caption') {
      targetFeatures = [feat, ...findFeaturesWithSameLegendItem(feat, legendName)];
      legendCaption = normalizeCaption(legendName);
    } else if (choice === 'displayLabel') {
      targetFeatures = [feat, ...findFeaturesWithSameDisplayedLabel(feat, displayLabel)];
    } else if (choice === 'annotationLabel') {
      targetFeatures = [feat, ...findFeaturesWithSameIndividualLabel(feat, annotationLabel)];
    } else if (choice === 'single') {
      targetFeatures = [feat];
    }

    const featureChanged = applyStrokeToSelectedFeatures(targetFeatures, strokeColor, strokeWidth);
    const legendChanged = legendCaption
      ? applyStrokeToLegendEntry(legendCaption, strokeColor, strokeWidth)
      : false;
    clearFeatureStyleScopeDialog();
    return featureChanged || legendChanged;
  };

  const handleFeatureStyleScopeChoice = (choice) => (
    featureStyleScopeDialog.kind === 'stroke'
      ? handleStrokeScopeChoice(choice)
      : handleColorScopeChoice(choice)
  );

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
    handleColorScopeChoice: colorAction(handleColorScopeChoice),
    handleFeatureStyleScopeChoice: colorAction(handleFeatureStyleScopeChoice),
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
    setClickedFeatureStrokeWidthValue: colorAction(setClickedFeatureStrokeWidthValue),
    setFeatureColor: colorAction(setFeatureColor),
    setFeatureColorValue: colorAction(setFeatureColorValue),
    updateClickedFeatureColor: colorAction(updateClickedFeatureColor),
    updateClickedFeatureStroke: colorAction(updateClickedFeatureStroke)
  };
};
