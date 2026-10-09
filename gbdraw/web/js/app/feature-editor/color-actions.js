// @ts-check
/** @import { DrawingState } from '../../state.js' */
import { reportRuleRunFailure } from '../rule-matching.js';
import { matchedRuleKeys, ruleKey, ruleMatcher, ruleMatchesFeature } from '../../services/rule-matchers.js';
import { resolveColorToHex } from '../../utils/color-utils.js';
import { getFeatureCaption, getFeatureColorRuleHash, getFeatureHashCandidates } from '../../services/feature-utils.js';
import { exactRegexValue } from '../../services/feature-selector.js';
import { getAllFeatureLegendGroups, mountedLegendRowFeatureIds, setsFeatureStroke } from '../../services/legend-svg.js';
import {
  featureOverrideKey,
  getFeatureOverride
} from '../../services/feature-override-identity.js';

/**
 * The rule owner's actions this owner reads (a whole-object port, R13).
 * @typedef {object} ColorActionsRuleActions
 * @property {(rules: Record<string, any>[], label?: string, options?: Record<string, any>) => Promise<any>} commitSpecificRules
 * @property {(rule: Record<string, any>) => number} countFeaturesMatchingRule
 * @property {(feature: Record<string, any>, caption: string) => { rule: Record<string, any>, color: string } | null} findExistingColorForCaption
 * @property {(feature: Record<string, any>, label?: string | null) => Record<string, any>[]} findFeaturesWithSameDisplayedLabel
 * @property {(feature: Record<string, any>, label?: string | null) => Record<string, any>[]} findFeaturesWithSameIndividualLabel
 * @property {(feature: Record<string, any>, caption?: string | null) => Record<string, any>[]} findFeaturesWithSameLegendItem
 * @property {(feature: Record<string, any>) => Record<string, any> | null} findMatchingRegexRule
 * @property {(feature: Record<string, any>) => string} getDisplayedFeatureLabel
 * @property {() => (feature: Record<string, any>) => string} effectiveLegendCaptions
 *   The Legend item of a feature, read for a pass over many features.
 * @property {(feature: Record<string, any>) => string} getIndividualFeatureLabel
 * @property {(feature: Record<string, any>) => { qual: string, val: string } | null} getFeatureQualifier
 * @property {(feature: Record<string, any>, label: string) => { feat: string, qual: string, val: string } | null} getLabelSpecificRule
 * @property {(caption: string) => Record<string, any>[]} getLegendRowRules
 * @property {(rules: Record<string, any>[], commit: () => any) => any} runWithRuleMatches
 *   Runs an action once the color rule matches of `rules` are prepared: the rule owner prepares the rules it builds.
 */

/**
 * @typedef {object} FeatureColorActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(options?: { replaceGeneratedInventory?: boolean }) => any} extractLegendEntries
 *   The Legend owner's reading of the mounted Legend rows.
 * @property {() => void} onLegendGeometryChanged The Legend owner's reaction to a change of Legend geometry.
 * @property {ColorActionsRuleActions} ruleActions The rule owner's lookups and commit of specific-color rules.
 * @property {(svg: Element, featureId: string) => Element[]} getFeatureElements The mounted elements of a feature.
 * @property {(svg: Element, featureId: string) => Element[]} getFeatureFillElements The mounted fill elements of a feature.
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 */

/** @param {FeatureColorActionsOptions} options */
export const createFeatureColorActions = ({
  state,
  extractLegendEntries,
  onLegendGeometryChanged,
  ruleActions,
  // R13: the mounted feature element lookups and the preview owner's commit
  // of an edit to the displayed Result.
  getFeatureElements,
  getFeatureFillElements,
  commitActiveResultEdit = null
}) => {
  const {
    appliedPaletteColors,
    extractedFeatures,
    biologicalFeatures,
    svgContainer,
    clickedFeature,
    featureStyleScopeDialog,
    resetColorDialog,
    legendRenameDialog,
    originalLegendOrder,
    originalLegendColors,
    originalSvgStroke
  } = state;

  const {
    countFeaturesMatchingRule,
    findExistingColorForCaption,
    findFeaturesWithSameDisplayedLabel,
    findFeaturesWithSameIndividualLabel,
    findFeaturesWithSameLegendItem,
    findMatchingRegexRule,
    getDisplayedFeatureLabel,
    effectiveLegendCaptions,
    getIndividualFeatureLabel,
    getFeatureQualifier,
    getLabelSpecificRule,
    getLegendRowRules,
    runWithRuleMatches
  } = ruleActions;
  const getEffectiveLegendCaption = (feature) => effectiveLegendCaptions()(feature);
  const normalizeCaption = (value) => String(value || '').trim();
  const normalizeCaptionKey = (value) => normalizeCaption(value).toLowerCase();
  const normalizeColor = (value) => String(value || '').trim().toLowerCase();
  const captionsMatch = (left, right) => normalizeCaptionKey(left) === normalizeCaptionKey(right);
  const colorsMatch = (left, right) => normalizeColor(left) === normalizeColor(right);
  const isHashSpecificRule = (rule) => String(rule?.qual || '').toLowerCase() === 'hash';
  const hasOwn = (object, key) => Object.prototype.hasOwnProperty.call(object || {}, key);
  // The DOM edits of one color action reach the Result in one commit, when the
  // last running color action settles.
  let colorActionDepth = 0;
  let pendingCommitReason = '';

  const runColorAction = async (action) => {
    colorActionDepth += 1;
    try {
      return await action();
    } finally {
      colorActionDepth -= 1;
      if (colorActionDepth === 0 && pendingCommitReason) {
        const reason = pendingCommitReason;
        pendingCommitReason = '';
        commitActiveResultEdit?.(reason);
      }
    }
  };

  // Runs `run` once the rules a color edit may add for the features in `args`,
  // the clicked feature and the scope dialog's feature (their hash and label
  // rules) are prepared with the saved ones.
  /**
   * @param {DrawingState} drawing
   * @param {any[]} args
   * @param {() => any} run
   */
  const withTargetRules = (drawing, args, run) => {
    const targets = new Set();
    const add = (feature) => { if (feature?.type && feature?.svg_id) targets.add(feature); };
    args.forEach((arg) => Array.isArray(arg) ? arg.forEach(add) : add(arg));
    add(clickedFeature.value?.feat);
    add(featureStyleScopeDialog.feat);
    const candidates = [...drawing.manualSpecificRules];
    targets.forEach((feature) => {
      const hash = getFeatureQualifier(feature);
      if (hash) candidates.push({ feat: feature.type, ...hash });
      for (const label of [getDisplayedFeatureLabel(feature), getIndividualFeatureLabel(feature)]) {
        const rule = getLabelSpecificRule(feature, label);
        if (rule) candidates.push(rule);
      }
    });
    return runWithRuleMatches(candidates, run);
  };

  // A color action also prepares the rules it may add (`withTargetRules`); a
  // stroke action reads only the saved rules (OV-198), and a color request
  // prepares its rules only when it commits, not to open the scope dialog,
  // which reads the saved rules only (OV-225).
  const colorAction = (action, { targetRules = true } = {}) => (...args) => {
    const drawing = state.activeDrawing();
    const run = () => runColorAction(() => action(drawing, ...args));
    const prepareTargets = () => (targetRules ? withTargetRules(drawing, args, run) : run());
    return reportRuleRunFailure(state, 'evaluateRules', () => runWithRuleMatches(drawing.manualSpecificRules, prepareTargets));
  };
  const strokeAction = (action) => colorAction(action, { targetRules: false });

  // The `hash` rule values that target a feature exactly.
  const exactHashValues = (feature) => {
    const candidates = getFeatureHashCandidates(feature);
    const generationHash = getFeatureColorRuleHash(feature);
    const renderedId = candidates[candidates.length - 1] || '';
    // Python matches only the stable hash, which duplicate records share.
    return [generationHash, ...(renderedId !== generationHash ? [renderedId] : [])]
      .filter(Boolean).flatMap((candidate) => [candidate, exactRegexValue(candidate)]);
  };
  const exactHashKey = (type, value) => JSON.stringify([type, value]);
  const ruleHashKey = (rule) => exactHashKey(rule.feat, String(rule?.val || '').trim());
  const hashRuleTargetsFeatureExactly = (rule, feature) => isHashSpecificRule(rule) && rule?.feat === feature?.type
    && exactHashValues(feature).includes(String(rule?.val || '').trim());
  // Whether a rule targets one of `features` exactly, read once per feature.
  const targetsSomeFeatureExactly = (features) => {
    const keys = new Set(features.flatMap((feature) => exactHashValues(feature)
      .map((value) => exactHashKey(feature.type, value))));
    return (rule) => isHashSpecificRule(rule) && keys.has(ruleHashKey(rule));
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

  /**
   * @param {DrawingState} drawing
   * @param {Record<string, any>} featureLike
   * @param {{ strokeColor?: string | null, strokeWidth?: number | null, originalStrokeColor?: string | null, originalStrokeWidth?: string | number | null }} [overrides]
   */
  const recordFeatureStrokeOverride = (
    drawing, featureLike,
    { strokeColor = null, strokeWidth = null, originalStrokeColor = null, originalStrokeWidth = null } = {}
  ) => {
    const key = featureStrokeKey(featureLike, featureLike?.svg_id);
    if (!key) return;

    const existing = drawing.featureStrokeOverrides[key] || {};
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
      drawing.featureStrokeOverrides[key] = next;
    }
  };

  /** @param {DrawingState} drawing */
  const clearFeatureStrokeOverride = (drawing, featureLike, fallbackSvgId = '') => {
    const key = featureStrokeKey(featureLike, fallbackSvgId);
    if (key) delete drawing.featureStrokeOverrides[key];
  };

  // The stroke a Legend row edit gives a feature without a stroke edit of its
  // own, which the feature shows once its own edit is removed, as Generate
  // draws it (`legendRowFeatureIds`, OV-123). Null when no row's stroke reaches it.
  /**
   * @param {DrawingState} drawing
   * @param {Element} svg
   * @param {Record<string, any>} feature
   * @param {string} svgId
   * @returns {Record<string, any> | null}
   */
  const legendRowStrokeOf = (drawing, svg, feature, svgId) => {
    const namedCaption = normalizeCaption(getFeatureOverride(drawing.featureColorOverrides, feature)?.caption);
    return Object.entries(drawing.legendStrokeOverrides).find(([caption]) => mountedLegendRowFeatureIds(
      svg, caption, drawing.legendEntries.value, { namedIds: caption === namedCaption ? [svgId] : [] }
    ).includes(svgId))?.[1] || null;
  };

  /** @param {DrawingState} drawing */
  const findLegendEntryByCaption = (drawing, caption) => {
    const normalizedCaption = normalizeCaptionKey(caption);
    if (!normalizedCaption) return null;

    return (
      drawing.legendEntries.value.find(
        (entry) => normalizeCaptionKey(entry?.caption) === normalizedCaption
      ) || null
    );
  };

  /** @param {DrawingState} drawing */
  const findExistingCaptionColor = (drawing, feat, caption) => {
    const existingCaption = findExistingColorForCaption(feat, caption);
    if (existingCaption?.color) {
      return {
        caption: existingCaption.rule?.cap || caption,
        color: existingCaption.color,
        rule: existingCaption.rule || null
      };
    }

    const legendEntry = findLegendEntryByCaption(drawing, caption);
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

  /**
   * @param {Record<string, any>} feat
   * @param {string | null} [requestedLegendName]
   */
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

  /**
   * @param {{
   *   kind: string,
   *   feat: Record<string, any>,
   *   scope: Record<string, any>,
   *   color?: string | null,
   *   strokeColor?: string | null,
   *   strokeWidth?: number | null,
   *   existingCaption?: { caption: string, color: string, rule: Record<string, any> | null } | null,
   *   closePopup?: boolean
   * }} options
   */
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

  const persistCurrentSvg = (svg = getCurrentSvg(), reason = 'feature-color') => {
    if (!svg) return;
    pendingCommitReason ||= reason;
  };

  /** @param {DrawingState} drawing */
  const exactHashRulesForFeature = (drawing, feature) => drawing.manualSpecificRules.filter(
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

  /** @param {DrawingState} drawing */
  const featureColorAssignmentMatches = (
    drawing, feature,
    color,
    caption,
    { requireLegend = true } = {}
  ) => {
    const override = getFeatureOverride(drawing.featureColorOverrides, feature);
    if (!override || !colorsMatch(override.color, color) || !captionsMatch(override.caption, caption)) {
      return false;
    }
    const matchingRules = exactHashRulesForFeature(drawing, feature);
    if (!matchingRules.some(
      (rule) => colorsMatch(rule.color, color) && captionsMatch(rule.cap, caption)
    )) {
      return false;
    }
    if (!liveFeatureColorMatches(feature, color)) return false;
    if (!requireLegend) return true;
    const legendEntry = findLegendEntryByCaption(drawing, caption);
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


  /** @param {DrawingState} drawing */
  const moveAddedLegendCaption = (drawing, oldCaption, newCaption) => {
    if (!oldCaption || !newCaption || oldCaption === newCaption) return;
    let matchedCaption = null;
    for (const caption of drawing.addedLegendCaptions.value) {
      if (captionsMatch(caption, oldCaption)) {
        matchedCaption = caption;
        break;
      }
    }
    if (!matchedCaption) return;
    drawing.addedLegendCaptions.value.delete(matchedCaption);
    drawing.addedLegendCaptions.value.add(newCaption);
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


  const updateClickedFeatureLegendState = (feat, caption, color = null) => {
    if (!clickedFeature.value || !feat || clickedFeature.value.svg_id !== feat.svg_id) return;
    if (color) {
      clickedFeature.value.color = color;
    }
    clickedFeature.value.legendName = caption;
    clickedFeature.value.appliedLegendName = caption;
  };

  /** @param {DrawingState} drawing */
  const clearLegendRenameDialog = (drawing, { restoreInput = false } = {}) => {
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
        drawing.legendEntries.value = [...drawing.legendEntries.value];
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
    legendRenameDialog.mergeAvailable = true;
    legendRenameDialog.pendingRequest = null;
  };

  /** @param {DrawingState} drawing */
  const getCurrentFeatureFillColor = (drawing, feat) => {
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

    const overrideColor = getFeatureOverride(drawing.featureColorOverrides, feat)?.color;
    if (overrideColor) {
      return resolveColorToHex(overrideColor) || overrideColor;
    }

    const fallbackColor = appliedPaletteColors.value[feat.type] || '#cccccc';
    return resolveColorToHex(fallbackColor) || fallbackColor;
  };

  const getFeaturesForLegendCaption = (caption) => {
    const normalizedCaption = normalizeCaption(caption);
    if (!normalizedCaption) return [];
    const captionOf = effectiveLegendCaptions();
    return extractedFeatures.value.filter((feat) => captionsMatch(captionOf(feat), normalizedCaption));
  };

  /** @param {DrawingState} drawing */
  const getUniqueLegendCaption = (drawing, caption, options = {}) => {
    const normalizedCaption = normalizeCaption(caption);
    if (!normalizedCaption) return '';

    const ignoreCaptionKeys = new Set(
      (Array.isArray(options.ignoreCaptions) ? options.ignoreCaptions : [])
        .map((value) => normalizeCaptionKey(value))
        .filter(Boolean)
    );

    const existingKeys = new Set();
    drawing.legendEntries.value.forEach((entry) => {
      const key = normalizeCaptionKey(entry?.caption);
      if (key) existingKeys.add(key);
    });
    drawing.manualSpecificRules.forEach((rule) => {
      const key = normalizeCaptionKey(rule?.cap);
      if (key) existingKeys.add(key);
    });
    Object.values(drawing.featureColorOverrides).forEach((override) => {
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

  /** @param {DrawingState} drawing */
  // Each feature's `hash` rule replaces the first rule that targets it
  // exactly, else goes before the first `hash` rule of its type that matches
  // it, else last. The rules are found by index, read once per feature (OV-193).
  const featureRuleCandidate = (drawing, features, color, caption, { preferLabelRules = false } = {}) => {
    let rules = drawing.manualSpecificRules.map(rule => ({ ...rule }));
    const labelRule = preferLabelRules ? getSafeLabelSpecificRule(drawing, features, caption) : null;
    if (labelRule) {
      const targeted = targetsSomeFeatureExactly(features);
      rules = rules.filter(rule => !(targeted(rule)
        || (rule.feat === labelRule.feat && rule.qual === labelRule.qual && rule.val === labelRule.val)));
      const first = rules.findIndex(rule => rule.feat === labelRule.feat && rule.qual === labelRule.qual);
      rules.splice(first < 0 ? rules.length : first, 0, { ...labelRule, color, cap: caption });
      return rules;
    }
    /** @type {Map<string, Record<string, any>[]>} */
    const byExactHash = new Map();
    /** @type {Map<string, Record<string, any>[]>} */
    const byKey = new Map();
    /** @param {Record<string, any>} rule @param {boolean} add */
    const index = (rule, add) => {
      /** @param {Map<string, Record<string, any>[]>} map @param {string} key */
      const update = (map, key) => {
        const list = map.get(key);
        if (add) {
          if (list) list.push(rule);
          else map.set(key, [rule]);
        } else if (list) {
          map.set(key, list.filter(other => other !== rule));
        }
      };
      update(byKey, ruleKey(rule));
      if (isHashSpecificRule(rule)) update(byExactHash, ruleHashKey(rule));
    };
    rules.forEach(rule => index(rule, true));
    // Rule positions, rebuilt after an insertion moved the rules after it.
    /** @type {{ at: Map<Record<string, any>, number> | null }} */
    const order = { at: null };
    const positionOf = (rule) => (order.at ||= new Map(rules.map((other, at) => [other, at]))).get(rule) ?? -1;
    /** @param {Record<string, any>[]} candidates @returns {Record<string, any> | null} */
    const firstOf = (candidates) => candidates.reduce((best, rule) => (
      best === null || positionOf(rule) < positionOf(best) ? rule : best), /** @type {Record<string, any> | null} */ (null));
    for (const feature of features) {
      const qualifier = getFeatureQualifier(feature);
      if (!qualifier) continue;
      const next = { feat: feature.type, ...qualifier, color, cap: caption };
      const existing = firstOf(exactHashValues(feature).flatMap(value => byExactHash.get(exactHashKey(feature.type, value)) || []));
      if (existing) {
        const at = positionOf(existing);
        rules[at] = next;
        order.at?.delete(existing);
        order.at?.set(next, at);
        index(existing, false);
      } else {
        const conflicting = firstOf([...matchedRuleKeys(feature)].flatMap(key => byKey.get(key) || [])
          .filter(rule => rule.feat === feature.type && isHashSpecificRule(rule) && ruleMatchesFeature(feature, rule)));
        if (conflicting) {
          rules.splice(positionOf(conflicting), 0, next);
          order.at = null;
        } else {
          rules.push(next);
          order.at?.set(next, rules.length - 1);
        }
      }
      index(next, true);
    }
    return rules;
  };

  const featureIdentityKey = (feature) => String(
    feature?.id || `${feature?.type || ''}:${feature?.svg_id || ''}`
  );

  /** @param {DrawingState} drawing */
  const getSafeLabelSpecificRule = (drawing, features, label) => {
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

    // The selected features by the keys of the rules they match, read once each.
    /** @type {Map<string, Record<string, any>[]>} */
    const selectedByMatchedKey = new Map();
    features.forEach((feature) => {
      for (const key of matchedRuleKeys(feature)) {
        const selected = selectedByMatchedKey.get(key);
        if (selected) selected.push(feature);
        else selectedByMatchedKey.set(key, [feature]);
      }
    });
    const hasPrecedenceConflict = drawing.manualSpecificRules.some((existing) => {
      const matchingSelected = (selectedByMatchedKey.get(ruleKey(existing)) || [])
        .filter((feature) => ruleMatchesFeature(feature, existing));
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

  /** @param {DrawingState} drawing */
  const renameLegendEntryInSvg = (drawing, oldCaption, newCaption, color = null) => {
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

    // The renamed row takes the caption: a style that an earlier row left under
    // that caption does not follow it, as Generate would otherwise apply it (OV-60).
    if (!findLegendEntryByCaption(drawing, newCaption)) {
      for (const store of [drawing.legendColorOverrides, drawing.legendStrokeOverrides]) {
        const staleKey = findCaptionKey(store, newCaption);
        if (staleKey) delete store[staleKey];
      }
    }
    moveCaptionStateKey(drawing.legendColorOverrides, oldCaption, newCaption);
    moveCaptionStateKey(drawing.legendStrokeOverrides, oldCaption, newCaption);
    moveAddedLegendCaption(drawing, oldCaption, newCaption);

    // A renderer-generated row keeps its generated caption as its identity, so
    // Generate replays the rename onto the regenerated row (PV-02). Rows the
    // editor added are identified by their current caption.
    const legendEntry = drawing.legendEntries.value.find((entry) => captionsMatch(entry?.caption, oldCaption));
    const generatedRow = Boolean(legendEntry) && originalLegendOrder.value.some(
      (caption) => captionsMatch(caption, legendEntry.originalCaption || legendEntry.caption)
    );
    if (!generatedRow) syncOriginalLegendMetadataRename(oldCaption, newCaption, color);
    if (legendEntry) {
      legendEntry.caption = newCaption;
      if (!generatedRow) legendEntry.originalCaption = newCaption;
      if (color) {
        legendEntry.color = color;
      }
    }

    // The layout owner lays the Legend out as Python would with the renamed
    // row (zero shift; OV-127) and docks it.
    onLegendGeometryChanged();
    persistCurrentSvg(svg);
    return true;
  };

  /** @param {DrawingState} drawing */
  const applyLegendRenameRequest = async (drawing, request) => {
    const oldCaption = normalizeCaption(request.oldCaption);
    const caption = normalizeCaption(request.finalCaption || request.newCaption);
    const color = resolveColorToHex(request.finalColor || request.currentColor) || '#cccccc';
    if (!caption) return false;
    const features = (request.features || []).filter(Boolean);
    const sourceRules = getLegendRowRules(oldCaption);
    if (sourceRules.length || features.length) {
      const rules = request.sourceScope === 'group' && sourceRules.length
        ? drawing.manualSpecificRules.map(rule => sourceRules.includes(rule) ? { ...rule, cap: caption, color } : { ...rule })
        : featureRuleCandidate(drawing, features, color, caption);
      const selected = new Set(features.map(feature => feature.svg_id));
      const retireOld = features.length > 0 && getFeaturesForLegendCaption(oldCaption)
        .every(feature => selected.has(feature.svg_id));
      const oldEntry = findLegendEntryByCaption(drawing, oldCaption);
      return ruleActions.commitSpecificRules(rules, 'Rename legend item', {
        drawing,
        previousLegendIntents: retireOld && oldEntry ? [{ caption: oldCaption, color: oldEntry.color }] : [],
        // OV-158 (Owner decision 2026-10-07): a renamed row keeps its place;
        // a row that joins another row keeps that row's place.
        legendPlacement: retireOld && oldEntry && !findLegendEntryByCaption(drawing, caption)
          ? { caption, at: oldEntry.caption } : null,
        afterCommit: () => {
          if (retireOld && !sourceRules.length) {
            const adoptedCaption = getEffectiveLegendCaption(features[0]);
            moveCaptionStateKey(drawing.legendColorOverrides, oldCaption, adoptedCaption);
            moveCaptionStateKey(drawing.legendStrokeOverrides, oldCaption, adoptedCaption);
            moveAddedLegendCaption(drawing, oldCaption, adoptedCaption);
            syncOriginalLegendMetadataRename(oldCaption, adoptedCaption, color);
          }
          const committedCaptionOf = effectiveLegendCaptions();
          for (const feature of features) updateClickedFeatureLegendState(feature, committedCaptionOf(feature), color);
        }
      });
    }
    // Unrelated manual legend rows retain their existing editor semantics.
    renameLegendEntryInSvg(drawing, oldCaption, caption, color);
    extractLegendEntries();
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

  const openLegendRenameTargetDialog = (request, targetEntry, mergeAvailable) => {
    legendRenameDialog.show = true;
    legendRenameDialog.mode = 'target';
    legendRenameDialog.oldCaption = request.oldCaption;
    legendRenameDialog.newCaption = request.newCaption;
    legendRenameDialog.targetCaption = targetEntry.caption;
    legendRenameDialog.targetColor = targetEntry.color || '';
    legendRenameDialog.currentColor = request.currentColor || '';
    legendRenameDialog.siblingCount = request.siblingCount || 0;
    legendRenameDialog.mergeAvailable = mergeAvailable;
    legendRenameDialog.pendingRequest = request;
  };

  /** @param {DrawingState} drawing */
  const continueLegendRenameRequest = async (drawing, request) => {
    if (!request) return;

    const oldCaption = normalizeCaption(request.oldCaption);
    const newCaption = normalizeCaption(request.newCaption);
    if (!oldCaption || !newCaption || newCaption === oldCaption) {
      clearLegendRenameDialog(drawing, { restoreInput: true });
      return;
    }

    const currentColor =
      resolveColorToHex(request.currentColor) ||
      (request.feat ? getCurrentFeatureFillColor(drawing, request.feat) : resolveColorToHex(findLegendEntryByCaption(drawing, oldCaption)?.color)) ||
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

    // D-06 (PD-OI-061): a rename onto another entry of a different color asks
    // Merge, Suffix, or Cancel, with or without features. A target owned by a
    // specific-color rule keeps the PD-OI-042 caption disambiguation instead.
    const targetEntry = findLegendEntryByCaption(drawing, newCaption);
    const isDistinctTargetEntry = targetEntry && !captionsMatch(targetEntry.caption, oldCaption);
    const ruleOwnedTarget = isDistinctTargetEntry && getLegendRowRules(targetEntry.caption).length > 0;
    const featureOrRuleRename = features.length > 0 || getLegendRowRules(oldCaption).length > 0;
    // OV-62 (PD-OI-061 amended): two rows merge only when each draws features of
    // one type and the type is the same. A row without features (GC content,
    // GC skew), rows of different types, and a row that spans several types
    // offer Suffix and Cancel only, even onto a rule-owned caption. The types
    // are those of the features that take each row's caption, as the live editor
    // knows them: a generated row such as `other proteins` has none live, so it
    // is not merged into.
    const featureTypes = (rowFeatures) => new Set(rowFeatures.map((feature) => feature?.type));
    const sourceTypes = featureTypes(features);
    const targetTypes = isDistinctTargetEntry ? featureTypes(getFeaturesForLegendCaption(targetEntry.caption)) : new Set();
    const mergeAllowed = isDistinctTargetEntry && sourceTypes.size === 1 && targetTypes.size === 1
      && [...sourceTypes][0] === [...targetTypes][0];

    if (featureOrRuleRename && (!isDistinctTargetEntry
      || (mergeAllowed && (ruleOwnedTarget || colorsMatch(targetEntry.color, currentColor))))) {
      await applyLegendRenameRequest(drawing, { ...request, currentColor, features,
        finalCaption: newCaption, finalColor: currentColor });
      clearLegendRenameDialog(drawing);
      return;
    }

    if (isDistinctTargetEntry && (!mergeAllowed || !colorsMatch(targetEntry.color, currentColor))) {
      if (!request.targetResolution) {
        openLegendRenameTargetDialog(
          {
            ...request,
            currentColor,
            features
          },
          targetEntry,
          mergeAllowed
        );
        return;
      }

      if (request.targetResolution === 'merge' && !mergeAllowed) {
        clearLegendRenameDialog(drawing, { restoreInput: true });
        return;
      }

      if (request.targetResolution === 'merge') {
        await applyLegendRenameRequest(drawing, {
          ...request,
          currentColor,
          features,
          finalCaption: targetEntry.caption,
          finalColor: targetEntry.color
        });
        clearLegendRenameDialog(drawing);
        return;
      }

      if (request.targetResolution === 'suffix') {
        await applyLegendRenameRequest(drawing, {
          ...request,
          currentColor,
          features,
          finalCaption: getUniqueLegendCaption(drawing, newCaption, { ignoreCaptions: [oldCaption] }),
          finalColor: currentColor
        });
        clearLegendRenameDialog(drawing);
        return;
      }
    }

    const finalCaption =
      isDistinctTargetEntry && colorsMatch(targetEntry.color, currentColor) ? targetEntry.caption : newCaption;
    const finalColor =
      isDistinctTargetEntry && colorsMatch(targetEntry.color, currentColor) ? targetEntry.color : currentColor;

    await applyLegendRenameRequest(drawing, {
      ...request,
      currentColor,
      features,
      finalCaption,
      finalColor
    });
    clearLegendRenameDialog(drawing);
  };

  /** @param {DrawingState} drawing */
  const applyColorToFeatureGroup = async (drawing, features, caption, color, options = {}) => {
    if (!features?.length || !normalizeCaption(caption)) return false;
    const existingEntry = findLegendEntryByCaption(drawing, caption);
    const contributors = getFeaturesForLegendCaption(caption);
    const captionOf = effectiveLegendCaptions();
    const selectedIds = new Set(features.map(feature => feature.svg_id));
    const replacesExistingGroup = existingEntry && contributors.length > 0
      && features.every(feature => captionsMatch(captionOf(feature), caption))
      && contributors.every(feature => selectedIds.has(feature.svg_id));
    const rules = featureRuleCandidate(drawing, features, color, normalizeCaption(caption), options);
    return ruleActions.commitSpecificRules(rules, 'Change feature color', {
      drawing,
      previousLegendIntents: replacesExistingGroup
        ? [{ caption: existingEntry.caption, color: existingEntry.color }] : [],
      afterCommit: intents => {
        if (replacesExistingGroup) {
          const intent = intents.find(entry => captionsMatch(entry.caption, caption)
            && colorsMatch(entry.color, color));
          if (intent) drawing.legendColorOverrides[intent.caption] = intent.color;
        }
        const committedCaptionOf = effectiveLegendCaptions();
        for (const feature of features) updateClickedFeatureLegendState(feature, committedCaptionOf(feature), color);
      }
    });
  };

  /** @param {DrawingState} drawing */
  const applyColorToLegendSpecificRules = async (drawing, caption, color) => {
    const rowRules = getLegendRowRules(caption);
    const specificRules = rowRules.filter(rule => !isHashSpecificRule(rule));
    if (!specificRules.length) return false;
    const specificMatches = ruleMatcher(specificRules);
    const coveredExactly = targetsSomeFeatureExactly(
      extractedFeatures.value.filter(feature => specificMatches.matchesAny(feature) === true)
    );
    const rowRuleSet = new Set(rowRules);
    const rules = drawing.manualSpecificRules.filter(rule => !(rowRuleSet.has(rule) && coveredExactly(rule)))
      .map(rule => rowRuleSet.has(rule) ? { ...rule, color } : { ...rule });
    return ruleActions.commitSpecificRules(rules, 'Change legend color', { drawing });
  };

  /**
   * @param {DrawingState} drawing
   * @param {Record<string, any>} feat
   * @param {string} color
   * @param {string | null} [requestedLegendName]
   * @param {{ closePopupOnDialog?: boolean }} [options]
   */
  const requestFeatureColorChange = async (drawing, feat, color, requestedLegendName = null, options = {}) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!feat) return;
    const scope = getFeatureStyleScope(feat, requestedLegendName);
    if (!scope) return;

    if (scope.needsDialog) {
      openFeatureStyleScopeDialog({
        kind: 'fill',
        feat,
        scope,
        color,
        existingCaption: findExistingCaptionColor(drawing, feat, scope.legendName),
        closePopup: options.closePopupOnDialog
      });
      return;
    }

    return withTargetRules(drawing, [feat], async () => {
      if (clickedFeature.value && clickedFeature.value.svg_id === feat.svg_id) {
        clickedFeature.value.color = color;
        if (scope.requestedCaption) {
          clickedFeature.value.legendName = scope.requestedCaption;
        }
      }
      await setFeatureColor(drawing, feat, color, scope.legendName);
    });
  };

  /** @param {DrawingState} drawing */
  const updateClickedFeatureColor = async (drawing, color) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!clickedFeature.value) return;
    const feat = clickedFeature.value.feat;
    if (!feat) return;
    const customName = normalizeCaption(clickedFeature.value.legendName);
    await requestFeatureColorChange(drawing, feat, color, customName, { closePopupOnDialog: true });
  };

  /** @param {DrawingState} drawing */
  const handleLegendNameCommit = async (drawing) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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

    await continueLegendRenameRequest(drawing, {
      source: 'popup',
      feat,
      oldCaption: currentCaption,
      newCaption: requestedCaption,
      currentColor: getCurrentFeatureFillColor(drawing, feat),
      sourceScope: null
    });
  };

  /** @param {DrawingState} drawing */
  const selectLegendNameOption = async (drawing, caption) => {
    if (!clickedFeature.value) return;
    const selectedCaption = String(caption || '').trim();
    if (!selectedCaption) return;
    clickedFeature.value.legendName = selectedCaption;
    await handleLegendNameCommit(drawing);
  };

  /** @param {DrawingState} drawing */
  const handleLegendRenameChoice = async (drawing, choice) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const pendingRequest = legendRenameDialog.pendingRequest;
    if (!pendingRequest || choice === 'cancel') {
      clearLegendRenameDialog(drawing, { restoreInput: true });
      return;
    }

    if (legendRenameDialog.mode === 'scope') {
      if (choice === 'single') {
        await continueLegendRenameRequest(drawing, {
          ...pendingRequest,
          sourceScope: 'single',
          targetResolution: null
        });
        return;
      }

      if (choice === 'group') {
        await continueLegendRenameRequest(drawing, {
          ...pendingRequest,
          sourceScope: 'group',
          targetResolution: null
        });
        return;
      }
    }

    if (legendRenameDialog.mode === 'target') {
      if (choice === 'merge') {
        await continueLegendRenameRequest(drawing, {
          ...pendingRequest,
          targetResolution: 'merge'
        });
        return;
      }

      if (choice === 'suffix') {
        await continueLegendRenameRequest(drawing, {
          ...pendingRequest,
          targetResolution: 'suffix'
        });
        return;
      }
    }

    clearLegendRenameDialog(drawing, { restoreInput: true });
  };

  /** @param {DrawingState} drawing */
  const renameLegendEntry = async (drawing, idx, newCaption) => {
    const entry = drawing.legendEntries.value[idx];
    if (!entry) return;

    const requestedCaption = normalizeCaption(newCaption);
    if (!requestedCaption || requestedCaption === normalizeCaption(entry.caption)) {
      drawing.legendEntries.value = [...drawing.legendEntries.value];
      return;
    }

    await continueLegendRenameRequest(drawing, {
      source: 'legend',
      oldCaption: normalizeCaption(entry.caption),
      newCaption: requestedCaption,
      currentColor: resolveColorToHex(entry.color) || entry.color || '#cccccc',
      features: getFeaturesForLegendCaption(entry.caption),
      sourceScope: getFeaturesForLegendCaption(entry.caption).length > 0 ? 'group' : 'manual'
    });
  };

  /** @param {DrawingState} drawing */
  const handleColorScopeChoice = async (drawing, choice) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const { feat, color, matchingRule, legendName, existingCaptionColor } = featureStyleScopeDialog;
    if (choice === 'cancel' || !feat || !color) {
      clearFeatureStyleScopeDialog();
      return;
    }

    if (choice === 'rule') {
      if (matchingRule) await ruleActions.commitSpecificRules(drawing.manualSpecificRules.map(rule => rule === matchingRule
        ? { ...rule, color } : { ...rule }), 'Change specific color rule', { drawing });
    } else if (choice === 'caption') {
      const targetLegendName = normalizeCaption(legendName) || normalizeCaption(getEffectiveLegendCaption(feat));
      if (!targetLegendName) {
        clearFeatureStyleScopeDialog();
        return;
      }
      const siblings = findFeaturesWithSameLegendItem(feat, targetLegendName);
      const allFeatures = [feat, ...siblings];
      if (!(await applyColorToLegendSpecificRules(drawing, targetLegendName, color))) {
        await applyColorToFeatureGroup(drawing, allFeatures, targetLegendName, color);
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
      await applyColorToFeatureGroup(drawing, allFeatures, displayLabel, color, { preferLabelRules: true });
    } else if (choice === 'single') {
      let singleCaption = legendName;
      if (featureStyleScopeDialog.siblingCount > 0 || (matchingRule && featureStyleScopeDialog.ruleMatchCount > 1)) {
        const ruleCaption = matchingRule ? (matchingRule.cap || matchingRule.val) : legendName;
        if (legendName === ruleCaption) {
          singleCaption = getIndividualFeatureLabel(feat);
        }
      }
      await setFeatureColor(drawing, feat, color, singleCaption);
    } else if (choice === 'annotationLabel') {
      const annotationLabel =
        normalizeCaption(featureStyleScopeDialog.annotationLabel) || normalizeCaption(getIndividualFeatureLabel(feat));
      if (!annotationLabel) {
        clearFeatureStyleScopeDialog();
        return;
      }
      const annotationSiblings = findFeaturesWithSameIndividualLabel(feat, annotationLabel);
      const allFeatures = [feat, ...annotationSiblings];
      await applyColorToFeatureGroup(drawing, allFeatures, annotationLabel, color, { preferLabelRules: true });
    } else if (choice === 'useExisting') {
      if (existingCaptionColor) {
        const targetLegendName = normalizeCaption(legendName) || normalizeCaption(getEffectiveLegendCaption(feat));
        await setFeatureColor(drawing, feat, existingCaptionColor, targetLegendName);
      }
    }

    clearFeatureStyleScopeDialog();
  };

  /** @param {DrawingState} drawing */
  const updateClickedFeatureStroke = (drawing, strokeColor, strokeWidth) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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
        element.setAttribute('stroke-width', /** @type {any} */ (normalizedStrokeWidth));
        changed = true;
      }
    });

    if (!changed) return false;
    recordFeatureStrokeOverride(drawing, clickedFeature.value.feat || clickedFeature.value, {
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

  /** @param {DrawingState} drawing */
  const requestClickedFeatureStrokeChange = (drawing, strokeColor, strokeWidth) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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

    return updateClickedFeatureStroke(drawing, normalizedStrokeColor, normalizedStrokeWidth);
  };

  /** @param {DrawingState} drawing */
  const resetClickedFeatureStroke = (drawing) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!clickedFeature.value) return false;
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const svgId = clickedFeature.value.svg_id;
    const elements = getFeatureElements(svg, svgId);

    // Without its own stroke edit, the feature shows its Legend row's stroke.
    const rowStroke = legendRowStrokeOf(drawing, svg, clickedFeature.value.feat || clickedFeature.value, svgId);
    const originalColor = String(rowStroke?.strokeColor || '').trim() || originalSvgStroke.value.color;
    const originalWidth = normalizeStrokeWidthValue(rowStroke?.strokeWidth)
      ?? normalizeStrokeWidthValue(originalSvgStroke.value.width);
    let changed = false;

    elements.forEach((element) => {
      if (!strokeColorAttributeMatches(element, originalColor)) {
        if (originalColor === null) element.removeAttribute('stroke');
        else element.setAttribute('stroke', originalColor);
        changed = true;
      }
      if (!strokeWidthAttributeMatches(element, originalWidth)) {
        if (originalWidth === null) element.removeAttribute('stroke-width');
        else element.setAttribute('stroke-width', /** @type {any} */ (originalWidth));
        changed = true;
      }
    });

    const feature = clickedFeature.value.feat || clickedFeature.value;
    const overrideKey = featureStrokeKey(feature, svgId);
    const hadOverride = Boolean(overrideKey && drawing.featureStrokeOverrides[overrideKey]);
    if (!changed && !hadOverride) return false;
    clickedFeature.value.strokeColor = originalColor || '';
    clickedFeature.value.strokeWidth = originalWidth ?? '';
    clearFeatureStrokeOverride(drawing, feature, svgId);

    if (changed) persistCurrentSvg(svg, 'feature-stroke');
    return true;
  };

  const getFeatureStrokeColorValue = (featureLike) => {
    const drawing = state.activeDrawing();
    const feature = featureLike?.feat || featureLike;
    const override = getFeatureOverride(drawing.featureStrokeOverrides, feature);
    return override && hasOwn(override, 'strokeColor')
      ? override.strokeColor
      : null;
  };

  /** @param {DrawingState} drawing */
  const setClickedFeatureStrokeColorValue = (drawing, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (value !== null) {
      if (!clickedFeature.value) return false;
      const feature = clickedFeature.value.feat || clickedFeature.value;
      const override = getFeatureOverride(drawing.featureStrokeOverrides, feature);
      if (
        !hasOwn(override, 'strokeColor')
        && colorsMatch(resolveColorToHex(value), resolveColorToHex(clickedFeature.value.strokeColor))
      ) {
        const normalizedValue = String(value || '').trim();
        if (updateClickedFeatureStroke(drawing, normalizedValue, null)) return true;
        const svg = getCurrentSvg();
        if (!svg) return false;
        const elements = getFeatureElements(svg, clickedFeature.value.svg_id);
        if (elements.length === 0) return false;
        recordFeatureStrokeOverride(drawing, feature, {
          strokeColor: normalizedValue,
          originalStrokeColor: clickedFeature.value.originalStrokeColor ?? elements[0]?.getAttribute('stroke') ?? null,
          originalStrokeWidth: clickedFeature.value.originalStrokeWidth ?? elements[0]?.getAttribute('stroke-width') ?? null
        });
        return true;
      }
      return requestClickedFeatureStrokeChange(drawing, value, clickedFeature.value.strokeWidth);
    }
    if (!clickedFeature.value || !svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;
    const feature = clickedFeature.value.feat || clickedFeature.value;
    const key = featureStrokeKey(feature, clickedFeature.value.svg_id);
    const override = key ? drawing.featureStrokeOverrides[key] : null;
    // Without a stroke edit of its own left, the feature shows its Legend row's stroke.
    const rowStroke = setsFeatureStroke({ strokeWidth: override?.strokeWidth })
      ? null
      : legendRowStrokeOf(drawing, svg, feature, clickedFeature.value.svg_id);
    const inheritedColor = String(rowStroke?.strokeColor || '').trim() || (override && hasOwn(override, 'originalStrokeColor')
      ? override.originalStrokeColor
      : (clickedFeature.value.originalStrokeColor ?? originalSvgStroke.value.color));
    const elements = getFeatureElements(svg, clickedFeature.value.svg_id);
    const domChanged = elements.some((element) => !strokeColorAttributeMatches(element, inheritedColor));
    const stateChanged = Boolean(override && hasOwn(override, 'strokeColor'));
    if (!domChanged && !stateChanged) return false;
    if (override) {
      delete override.strokeColor;
      if (!hasOwn(override, 'strokeWidth')) delete drawing.featureStrokeOverrides[key];
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

  /** @param {DrawingState} drawing */
  const setClickedFeatureStrokeWidthValue = (drawing, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!clickedFeature.value) return false;
    const normalizedStrokeWidth = normalizeStrokeWidthValue(value);
    const currentStrokeWidth = normalizeStrokeWidthValue(clickedFeature.value.strokeWidth);
    if (normalizedStrokeWidth === null || normalizedStrokeWidth === currentStrokeWidth) return false;
    return requestClickedFeatureStrokeChange(drawing, clickedFeature.value.strokeColor, normalizedStrokeWidth);
  };

  /** @param {DrawingState} drawing */
  const resetClickedFeatureFillColor = (drawing) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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
      resetColorDialog.siblingCount = siblings.length;
    } else {
      return doResetFillColor(drawing, 'this');
    }
  };

  /** @param {DrawingState} drawing */
  const handleResetColorChoice = async (drawing, choice) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    resetColorDialog.show = false;
    await doResetFillColor(drawing, choice);
  };

  /** @param {DrawingState} drawing */
  const doResetFillColor = async (drawing, choice) => {
    const feature = clickedFeature.value?.feat;
    if (!feature || choice === 'cancel') return false;
    const caption = getEffectiveLegendCaption(feature);
    // The reset color is the palette default of the feature being reset.
    const color = appliedPaletteColors.value[feature.type];
    if (choice === 'this_with_legend') return setFeatureColor(drawing, feature, color, caption);
    let rules = drawing.manualSpecificRules.filter(rule => choice === 'all'
      ? rule.cap !== caption : !hashRuleTargetsFeatureExactly(rule, feature));
    if (choice === 'this' && rules.some(rule => ruleMatchesFeature(feature, rule))) {
      const qualifier = getFeatureQualifier(feature);
      if (qualifier) rules.push({ feat: feature.type, ...qualifier, color,
        cap: feature.type === 'CDS' ? 'other proteins' : `other ${feature.type}s` });
    }
    const applied = await ruleActions.commitSpecificRules(rules, 'Reset feature color', { drawing });
    if (applied) clickedFeature.value = null;
    return applied;
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

  /** @param {DrawingState} drawing */
  const applyColorToSelectedFeatures = async (drawing, features, color, caption) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const targetFeatures = uniqueFeaturesBySvgId(features);
    const targetColor = resolveColorToHex(color) || String(color || '').trim();
    const targetCaption = normalizeCaption(caption);
    if (targetFeatures.length === 0 || !targetColor || !targetCaption) return false;
    await applyColorToFeatureGroup(drawing, targetFeatures, targetCaption, targetColor);
    return true;
  };

  /** @param {DrawingState} drawing */
  const applyStrokeToSelectedFeatures = (drawing, features, strokeColor, strokeWidth) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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
      recordFeatureStrokeOverride(drawing, feature, {
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
          element.setAttribute('stroke-width', /** @type {any} */ (normalizedStrokeWidth));
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

  /** @param {DrawingState} drawing */
  const applyStrokeToLegendEntry = (drawing, caption, strokeColor, strokeWidth) => {
    const targetLegendEntry = findLegendEntryByCaption(drawing, caption);
    const svg = getCurrentSvg();
    if (!targetLegendEntry || !svg) return false;

    const normalizedStrokeColor = String(strokeColor || '').trim();
    const normalizedStrokeWidth = normalizeStrokeWidthValue(strokeWidth);
    if (!normalizedStrokeColor && normalizedStrokeWidth === null) return false;

    let domChanged = false;
    /** @type {{ color: string | null, width: number | null } | null} */
    let originalSwatchStroke = null;
    const escapedCaption = globalThis.CSS?.escape
      ? globalThis.CSS.escape(targetLegendEntry.caption)
      : String(targetLegendEntry.caption).replace(/["\\]/g, '\\$&');
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
        swatch.setAttribute('stroke', normalizedStrokeColor);
        domChanged = true;
      }
      if (normalizedStrokeWidth !== null && !strokeWidthAttributeMatches(swatch, normalizedStrokeWidth)) {
        swatch.setAttribute('stroke-width', normalizedStrokeWidth);
        domChanged = true;
      }
    }

    const overrideKey = targetLegendEntry.caption;
    const previousOverride = drawing.legendStrokeOverrides[overrideKey] || {};
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
    if (stateChanged) drawing.legendStrokeOverrides[overrideKey] = nextOverride;
    if (domChanged) persistCurrentSvg(svg, 'feature-stroke');
    return domChanged || stateChanged;
  };

  /** @param {DrawingState} drawing */
  const handleStrokeScopeChoice = (drawing, choice) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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
        ruleMatchesFeature(candidate, matchingRule)
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

    const featureChanged = applyStrokeToSelectedFeatures(drawing, targetFeatures, strokeColor, strokeWidth);
    const legendChanged = legendCaption
      ? applyStrokeToLegendEntry(drawing, legendCaption, strokeColor, strokeWidth)
      : false;
    clearFeatureStyleScopeDialog();
    return featureChanged || legendChanged;
  };

  /**
   * @param {DrawingState} drawing
   * @param {Record<string, any>} feature
   * @param {string} color
   * @param {string | null} [customCaption]
   */
  const setFeatureColor = async (drawing, feature, color, customCaption = null) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!feature || !getFeatureQualifier(feature)) return false;
    const caption = normalizeCaption(customCaption || getIndividualFeatureLabel(feature));
    if (!caption || featureColorAssignmentMatches(drawing, feature, color, caption)) return false;
    return applyColorToFeatureGroup(drawing, [feature], caption, resolveColorToHex(color) || color);
  };

  /** @param {DrawingState} drawing */
  const setFeatureColorValue = async (drawing, feature, value, customCaption = null) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!feature) return false;
    if (value === null) {
      return ruleActions.commitSpecificRules(drawing.manualSpecificRules.filter(rule => !hashRuleTargetsFeatureExactly(rule, feature)), 'Reset feature color', { drawing });
    }
    if (String(value).trim().toLowerCase() === 'none') {
      return applyColorToFeatureGroup(drawing, [feature], normalizeCaption(customCaption || getEffectiveLegendCaption(feature) || feature.type), 'none');
    }
    return setFeatureColor(drawing, feature, value, customCaption);
  };

  // The dialog's kind decides at call time which rules the choice prepares.
  const strokeScopeChoice = strokeAction(handleStrokeScopeChoice);
  const colorScopeChoice = colorAction(handleColorScopeChoice);

  return {
    // OV-161: Cancel of a scope dialog only closes it: no rule matching, no
    // History step (app-setup.js answers Cancel outside `runUndoable`).
    cancelFeatureStyleScope: clearFeatureStyleScopeDialog,
    cancelLegendRename: () => clearLegendRenameDialog(state.activeDrawing(), { restoreInput: true }),
    handleColorScopeChoice: colorScopeChoice,
    handleFeatureStyleScopeChoice: (...args) => (featureStyleScopeDialog.kind === 'stroke' ? strokeScopeChoice : colorScopeChoice)(...args),
    handleLegendNameCommit: colorAction(handleLegendNameCommit),
    handleLegendRenameChoice: colorAction(handleLegendRenameChoice),
    renameLegendEntry: colorAction(renameLegendEntry),
    requestFeatureColorChange: colorAction(requestFeatureColorChange, { targetRules: false }),
    selectLegendNameOption: colorAction(selectLegendNameOption),
    handleResetColorChoice: colorAction(handleResetColorChoice),
    applyColorToSelectedFeatures: colorAction(applyColorToSelectedFeatures),
    applyStrokeToSelectedFeatures: strokeAction(applyStrokeToSelectedFeatures),
    resetClickedFeatureFillColor: colorAction(resetClickedFeatureFillColor),
    resetClickedFeatureStroke: strokeAction(resetClickedFeatureStroke),
    getFeatureStrokeColorValue,
    setClickedFeatureStrokeColorValue: strokeAction(setClickedFeatureStrokeColorValue),
    setClickedFeatureStrokeWidthValue: strokeAction(setClickedFeatureStrokeWidthValue),
    setFeatureColor: colorAction(setFeatureColor),
    setFeatureColorValue: colorAction(setFeatureColorValue),
    updateClickedFeatureColor: colorAction(updateClickedFeatureColor, { targetRules: false }),
    updateClickedFeatureStroke: strokeAction(updateClickedFeatureStroke)
  };
};
