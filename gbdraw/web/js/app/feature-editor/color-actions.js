// @ts-check
/** @import { DrawingState, LegendEntry } from '../../state.js' */
/** @import { LegendRowReach, PythonLegendKey, RenderedFeatureId } from '../../services/legend-svg.js' */
import { draftLegendRowColors, draftLegendRows, LIVE_EDIT_DOMAINS, namedLegendCaption } from '../candidate-render.js';
import { reportRuleRunFailure } from '../rule-matching.js';
import { matchedRuleKeys, ruleKey, ruleMatcher, ruleMatchesFeature } from '../../services/rule-matchers.js';
import { appliedFeatureColors, resolveColorToHex } from '../../utils/color-utils.js';
import { getFeatureCaption, getFeatureColorRuleHash, getFeatureHashCandidates } from '../../services/feature-utils.js';
import { exactRegexValue } from '../../services/feature-selector.js';
import {
  drawnLegendRowStroke,
  legendRowFeatureIds,
  pythonLegendRows,
  setsFeatureStroke
} from '../../services/legend-svg.js';
import { isAutoFeatureUnderlay } from '../../services/feature-dom.js';
import { pythonDrawnAttribute } from '../../services/result-paint-bases.js';
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
 * The palette owner's user default color of a key (D-15): its Default colors
 * value when it differs from the selected palette's color, else null.
 * @typedef {(drawing: DrawingState, key: string) => string | null} ReadUserDefaultColorPort
 */
/**
 * The palette owner's write of a key's default color, applied as an edit in
 * the Default colors list is (D-15).
 * @typedef {(drawing: DrawingState, key: string, color: string) => void} SetDefaultColorPort
 */

/**
 * @typedef {object} FeatureColorActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {((options: { domains: readonly string[] }) => unknown) | null} [showEditorIntent]
 *   The root's port of the editor intent onto the displayed Result (R1).
 * @property {ColorActionsRuleActions} ruleActions The rule owner's lookups and commit of specific-color rules.
 * @property {(svg: Element, featureId: string) => Element[]} getFeatureElements The mounted elements of a feature.
 * @property {(svg: Element, featureId: string) => Element[]} getFeatureFillElements The mounted fill elements of a feature.
 * @property {(close: () => unknown) => void} [closeAfterDialogChoice]
 *   Runs a dialog's close at once, or once the History step of the dialog's
 *   choice in flight ends (D-12, OIC-028).
 * @property {ReadUserDefaultColorPort} readUserDefaultColor
 * @property {SetDefaultColorPort} setDefaultColor
 */

/** @param {FeatureColorActionsOptions} options */
export const createFeatureColorActions = ({
  state,
  showEditorIntent = null,
  ruleActions,
  // R13: the mounted feature element lookups.
  getFeatureElements,
  getFeatureFillElements,
  closeAfterDialogChoice = (close) => { close(); },
  // R13, D-15: the palette owner's user default colors.
  readUserDefaultColor,
  setDefaultColor
}) => {
  const {
    extractedFeatures,
    biologicalFeatures,
    svgContainer,
    clickedFeature,
    featureStyleScopeDialog,
    resetColorDialog,
    legendRenameDialog,
    originalLegendOrder,
    originalLegendColors
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
  // The domains a Legend row rename shows: the row's structure, and the fill
  // and stroke an editor row takes under its new caption.
  const RENAMED_ROW_DOMAINS = [...LIVE_EDIT_DOMAINS.legendStructure, ...LIVE_EDIT_DOMAINS.legendFills, ...LIVE_EDIT_DOMAINS.strokes];

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
  // stroke action reads only the saved rules (OV-198). A color request, a
  // popup Legend rename, and a fill reset prepare them only when they commit,
  // not to open their dialog, which reads the saved rules only (OV-225).
  const colorAction = (action, { targetRules = true } = {}) => (...args) => {
    const drawing = state.activeDrawing();
    const run = async () => action(drawing, ...args);
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

  // The stroke Python drew on a feature's block (its first part that is not
  // an automatic underlay), which a Session 46 stroke edit keeps.
  /** @param {Element[]} elements */
  const drawnFeatureStroke = (elements) => {
    const block = elements.find((element) => !isAutoFeatureUnderlay(element)) || elements[0] || null;
    return {
      originalStrokeColor: pythonDrawnAttribute(block, 'stroke'),
      originalStrokeWidth: pythonDrawnAttribute(block, 'stroke-width')
    };
  };

  // The stroke a Legend row edit gives a feature without a stroke edit of its
  // own, which the feature shows once its own edit is removed, as the executor
  // draws it from Python's row (`legendRowFeatureIds`, OV-123, OV-288). The
  // stroked rows are the compile's (`draftLegendRows`): a deleted row strokes
  // nothing. Null when no row's stroke reaches it.
  /**
   * @param {DrawingState} drawing
   * @param {Element} svg
   * @param {Record<string, any>} feature
   * @param {string} svgId
   * @returns {{ strokeColor?: unknown, strokeWidth?: unknown } | null}
   */
  const legendRowStrokeOf = (drawing, svg, feature, svgId) => {
    const id = /** @type {RenderedFeatureId} */ (svgId);
    const drawnFills = [/** @type {[RenderedFeatureId, string]} */ ([id, getFeatureFillElements(svg, svgId)[0]?.getAttribute('fill') || ''])];
    const namedCaption = normalizeCaption(namedLegendCaption(getFeatureOverride(drawing.featureColorOverrides, feature), drawing.manualSpecificRules));
    const pythonRows = pythonLegendRows(svg);
    const originalOrder = originalLegendOrder.value || [];
    const rows = draftLegendRows({
      legendEntries: drawing.legendEntries.value, deletedLegendEntries: drawing.deletedLegendEntries.value,
      dormantLegendEntries: drawing.dormantLegendEntries.value, originalLegendOrder: originalOrder
    });
    const draft = draftLegendRowColors({
      rules: drawing.manualSpecificRules, pythonRows, features: state.extractedFeatures.value || [],
      originalLegendOrder: originalOrder, paletteColors: appliedFeatureColors(state)
    });
    return Object.entries(drawing.legendStrokeOverrides).find(([caption]) => {
      const row = rows.styledRow(caption);
      if (!row) return false;
      /** @type {LegendRowReach} */
      const reach = {
        listedIds: /** @type {RenderedFeatureId[]} */ (row.entry ? row.entry.featureIds : []),
        namedIds: caption === namedCaption ? [id] : [],
        ownStrokeIds: [],
        draftColor: draft.colorOf(row.targetCaption)
      };
      return legendRowFeatureIds(reach, pythonRows.get(row.targetCaption)?.color ?? null, drawnFills).length > 0;
    })?.[1] || null;
  };

  /** @param {LegendEntry[]} entries @param {string} caption @returns {LegendEntry | null} */
  const entryByCaption = (entries, caption) => {
    const normalizedCaption = normalizeCaptionKey(caption);
    if (!normalizedCaption) return null;
    return entries.find((entry) => normalizeCaptionKey(entry?.caption) === normalizedCaption) || null;
  };
  /** @param {DrawingState} drawing @param {string} caption */
  const findLegendEntryByCaption = (drawing, caption) => entryByCaption(drawing.legendEntries.value, caption);
  // R15-3 (OV-285): the deleted row a rename's caption names. The caption only
  // detects the conflict; the row is its Python key (an editor row's own
  // caption), which the composition root's Restore acts on.
  /** @param {DrawingState} drawing @param {string} caption */
  const findDeletedLegendEntryByCaption = (drawing, caption) => entryByCaption(drawing.deletedLegendEntries.value, caption);
  /** @param {LegendEntry} entry @returns {PythonLegendKey} */
  const deletedEntryKey = (entry) => /** @type {PythonLegendKey} */ (entry.originalCaption || entry.caption);

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

  // A dialog's choice closes it once the choice's History step ends.
  const clearFeatureStyleScopeDialog = () => closeAfterDialogChoice(() => {
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
    featureStyleScopeDialog.defaultColorType = null;
    featureStyleScopeDialog.replacedDefaultColor = null;
  });

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
    const siblings = findFeaturesWithSameLegendItem(feat, legendName);
    const siblingCount = siblings.length;
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
      siblings,
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
   *   defaultColorType?: string | null,
   *   replacedDefaultColor?: string | null,
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
    defaultColorType = null,
    replacedDefaultColor = null,
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
    featureStyleScopeDialog.defaultColorType = defaultColorType;
    featureStyleScopeDialog.replacedDefaultColor = replacedDefaultColor;
    if (closePopup) clickedFeature.value = null;
  };

  const getCurrentSvg = () => svgContainer.value?.querySelector('svg') || null;

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
  const clearLegendRenameDialog = (drawing, { restoreInput = false } = {}) => closeAfterDialogChoice(() => {
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
    legendRenameDialog.deletedTargetKey = '';
    legendRenameDialog.pendingRequest = null;
  });

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

    const fallbackColor = appliedFeatureColors(state)[feat.type] || '#cccccc';
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
    // A deleted row keeps its caption for Restore (R15-3).
    [...drawing.legendEntries.value, ...drawing.deletedLegendEntries.value].forEach((entry) => {
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

  // The rename of a Legend row without rules or features writes the intent
  // (U3a gap 3): the row's styles move to the new caption, a row of Python's
  // keeps its generated caption, so Generate and the port rename it
  // (`legendRenames`, PV-02), and an editor row is identified by its caption,
  // so the port removes the old row and adds the new one with its styles, in
  // its place in the order. The text control's open step holds the edit.
  /** @param {DrawingState} drawing @param {string} oldCaption @param {string} newCaption */
  const renameLegendRow = (drawing, oldCaption, newCaption) => {
    const legendEntry = drawing.legendEntries.value.find((entry) => captionsMatch(entry?.caption, oldCaption));
    if (!legendEntry) return false;
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
    const generatedRow = originalLegendOrder.value.some(
      (caption) => captionsMatch(caption, legendEntry.originalCaption || legendEntry.caption)
    );
    if (!generatedRow) {
      syncOriginalLegendMetadataRename(oldCaption, newCaption);
      legendEntry.originalCaption = newCaption;
    }
    legendEntry.caption = newCaption;
    showEditorIntent?.({ domains: RENAMED_ROW_DOMAINS });
    return true;
  };

  // A rename commits once the rules it may add are prepared; its dialogs read
  // only the saved rules (OV-225).
  /** @param {DrawingState} drawing */
  const applyLegendRenameRequest = (drawing, request) => withTargetRules(drawing, [], () => applyPreparedLegendRename(drawing, request));
  /** @param {DrawingState} drawing */
  const applyPreparedLegendRename = async (drawing, request) => {
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
    return renameLegendRow(drawing, oldCaption, caption);
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
    legendRenameDialog.deletedTargetKey = '';
    legendRenameDialog.pendingRequest = request;
  };

  /**
   * @param {{ oldCaption: string, newCaption: string, currentColor?: string, siblingCount?: number, [field: string]: unknown }} request
   * @param {LegendEntry} targetEntry
   * @param {boolean | null} mergeAvailable
   * @param {PythonLegendKey | ''} deletedTargetKey The key of a deleted target row, else ''.
   */
  const openLegendRenameTargetDialog = (request, targetEntry, mergeAvailable, deletedTargetKey) => {
    legendRenameDialog.show = true;
    legendRenameDialog.mode = 'target';
    legendRenameDialog.oldCaption = request.oldCaption;
    legendRenameDialog.newCaption = request.newCaption;
    legendRenameDialog.targetCaption = targetEntry.caption;
    legendRenameDialog.targetColor = targetEntry.color || '';
    legendRenameDialog.currentColor = request.currentColor || '';
    legendRenameDialog.siblingCount = request.siblingCount || 0;
    legendRenameDialog.mergeAvailable = mergeAvailable;
    legendRenameDialog.deletedTargetKey = deletedTargetKey;
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
    // R15-3 (OV-285): a rename onto a deleted row's caption always asks, as onto
    // a listed row; its Merge is a Restore of that row, then the merge.
    const listedTarget = findLegendEntryByCaption(drawing, newCaption);
    const deletedTarget = listedTarget ? null : findDeletedLegendEntryByCaption(drawing, newCaption);
    const targetEntry = listedTarget || deletedTarget;
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

    // A chosen Merge adopts the target's color, also once a Restore made a
    // rule-owned deleted target a listed one.
    if (featureOrRuleRename && request.targetResolution !== 'merge' && (!isDistinctTargetEntry
      || (!deletedTarget && mergeAllowed && (ruleOwnedTarget || colorsMatch(targetEntry.color, currentColor))))) {
      await applyLegendRenameRequest(drawing, { ...request, currentColor, features,
        finalCaption: newCaption, finalColor: currentColor });
      clearLegendRenameDialog(drawing);
      return;
    }

    if (isDistinctTargetEntry && (deletedTarget || !mergeAllowed || !colorsMatch(targetEntry.color, currentColor))) {
      if (!request.targetResolution) {
        openLegendRenameTargetDialog(
          {
            ...request,
            currentColor,
            features
          },
          targetEntry,
          mergeAllowed,
          deletedTarget ? deletedEntryKey(deletedTarget) : ''
        );
        return;
      }

      // A deleted row takes the merge only once the Restore returned it.
      if (request.targetResolution === 'merge' && (!mergeAllowed || deletedTarget)) {
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

  // D-15: a Legend row that names a feature type and that no Specific color
  // rule or per-feature edit draws is that type's palette row; "Apply to all"
  // on it sets the type's default color. Null for any other row.
  /**
   * @param {DrawingState} drawing
   * @param {string} caption
   * @param {ReadonlyArray<{ type?: unknown }>} features The row's features.
   * @returns {string | null}
   */
  const paletteRowType = (drawing, caption, features) => {
    const type = normalizeCaption(caption);
    if (!features.every((feature) => feature?.type === type) || getLegendRowRules(type).length > 0) return null;
    const matches = ruleMatcher(drawing.manualSpecificRules);
    return features.every((feature) => matches.matchesAny(feature) === false
      && !getFeatureOverride(drawing.featureColorOverrides, feature)) ? type : null;
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
      const defaultColorType = scope.siblingCount > 0
        ? paletteRowType(drawing, scope.legendName, [feat, ...scope.siblings]) : null;
      const replacedDefaultColor = defaultColorType ? readUserDefaultColor(drawing, defaultColorType) : null;
      openFeatureStyleScopeDialog({
        kind: 'fill',
        feat,
        scope,
        color,
        existingCaption: findExistingCaptionColor(drawing, feat, scope.legendName),
        // The dialog's line under "Apply to all" (D-15, Q2 B).
        defaultColorType,
        replacedDefaultColor,
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
    const { feat, color, matchingRule, legendName, existingCaptionColor, defaultColorType } = featureStyleScopeDialog;
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
      // The choice does what the dialog showed: its default-color line, or rules.
      if (defaultColorType) {
        // The row's swatch follows the default color again (svg-styles.js).
        delete drawing.legendColorOverrides[defaultColorType];
        setDefaultColor(drawing, defaultColorType, color);
        // OV-264: the Legend row takes the color in the same step, as a rule
        // commit's Legend change does; its color feeds the rule captions.
        drawing.legendEntries.value = drawing.legendEntries.value.map((entry) => (
          captionsMatch(entry.caption, defaultColorType) ? { ...entry, color } : entry));
      } else if (!(await applyColorToLegendSpecificRules(drawing, targetLegendName, color))) {
        const siblings = findFeaturesWithSameLegendItem(feat, targetLegendName);
        await applyColorToFeatureGroup(drawing, [feat, ...siblings], targetLegendName, color);
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

  // The stroke edits write the editor intent only. The composition root shows
  // them on the Result through the executor as one History step
  // (`editEditorIntent` in app/app-setup.js), which records Python's strokes.
  /** @param {DrawingState} drawing */
  const updateClickedFeatureStroke = (drawing, strokeColor, strokeWidth) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!clickedFeature.value) return false;
    const svg = getCurrentSvg();
    if (!svg) return false;

    const elements = getFeatureElements(svg, clickedFeature.value.svg_id);
    if (elements.length === 0) return false;
    const normalizedStrokeColor = strokeColor === null || strokeColor === undefined
      ? null
      : String(strokeColor).trim();
    const normalizedStrokeWidth = normalizeStrokeWidthValue(strokeWidth);
    if (normalizedStrokeColor === null && normalizedStrokeWidth === null) return false;
    const changed = elements.some((element) => (
      (normalizedStrokeColor !== null && !strokeColorAttributeMatches(element, normalizedStrokeColor))
      || (normalizedStrokeWidth !== null && !strokeWidthAttributeMatches(element, normalizedStrokeWidth))
    ));
    if (!changed) return false;
    recordFeatureStrokeOverride(drawing, clickedFeature.value.feat || clickedFeature.value, {
      strokeColor: normalizedStrokeColor,
      strokeWidth: normalizedStrokeWidth,
      ...drawnFeatureStroke(elements)
    });

    if (normalizedStrokeColor !== null) clickedFeature.value.strokeColor = normalizedStrokeColor;
    if (normalizedStrokeWidth !== null) clickedFeature.value.strokeWidth = normalizedStrokeWidth;
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
    const svg = getCurrentSvg();
    if (!svg) return false;

    const svgId = clickedFeature.value.svg_id;
    const feature = clickedFeature.value.feat || clickedFeature.value;
    const overrideKey = featureStrokeKey(feature, svgId);
    if (!overrideKey || !drawing.featureStrokeOverrides[overrideKey]) return false;
    // Without its own stroke edit, the feature shows its Legend row's stroke,
    // else the stroke Python drew.
    const rowStroke = legendRowStrokeOf(drawing, svg, feature, svgId);
    const drawn = drawnFeatureStroke(getFeatureElements(svg, svgId));
    clickedFeature.value.strokeColor = String(rowStroke?.strokeColor || '').trim() || drawn.originalStrokeColor || '';
    clickedFeature.value.strokeWidth = normalizeStrokeWidthValue(rowStroke?.strokeWidth)
      ?? normalizeStrokeWidthValue(drawn.originalStrokeWidth) ?? '';
    clearFeatureStrokeOverride(drawing, feature, svgId);
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
    if (!clickedFeature.value) return false;
    const svg = getCurrentSvg();
    if (!svg) return false;
    const feature = clickedFeature.value.feat || clickedFeature.value;
    const elements = getFeatureElements(svg, clickedFeature.value.svg_id);
    if (value !== null) {
      const override = getFeatureOverride(drawing.featureStrokeOverrides, feature);
      if (
        !hasOwn(override, 'strokeColor')
        && colorsMatch(resolveColorToHex(value), resolveColorToHex(clickedFeature.value.strokeColor))
      ) {
        const normalizedValue = String(value || '').trim();
        if (updateClickedFeatureStroke(drawing, normalizedValue, null)) return true;
        if (elements.length === 0) return false;
        recordFeatureStrokeOverride(drawing, feature, { strokeColor: normalizedValue, ...drawnFeatureStroke(elements) });
        return true;
      }
      return requestClickedFeatureStrokeChange(drawing, value, clickedFeature.value.strokeWidth);
    }
    const key = featureStrokeKey(feature, clickedFeature.value.svg_id);
    const override = key ? drawing.featureStrokeOverrides[key] : null;
    if (!override || !hasOwn(override, 'strokeColor')) return false;
    // Without a stroke edit of its own left, the feature shows its Legend row's stroke.
    const rowStroke = setsFeatureStroke({ strokeWidth: override.strokeWidth })
      ? null
      : legendRowStrokeOf(drawing, svg, feature, clickedFeature.value.svg_id);
    delete override.strokeColor;
    if (!hasOwn(override, 'strokeWidth')) delete drawing.featureStrokeOverrides[key];
    clickedFeature.value.strokeColor = String(rowStroke?.strokeColor || '').trim()
      || drawnFeatureStroke(elements).originalStrokeColor || '';
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

    const defaultColor = appliedFeatureColors(state)[feat.type];
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
      return;
    }
    return withTargetRules(drawing, [], () => doResetFillColor(drawing, 'this'));
  };

  const closeResetColorDialog = () => closeAfterDialogChoice(() => { resetColorDialog.show = false; });

  /** @param {DrawingState} drawing */
  const handleResetColorChoice = async (drawing, choice) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    await doResetFillColor(drawing, choice);
    closeResetColorDialog();
  };

  /** @param {DrawingState} drawing */
  const doResetFillColor = async (drawing, choice) => {
    const feature = clickedFeature.value?.feat;
    if (!feature || choice === 'cancel') return false;
    const caption = getEffectiveLegendCaption(feature);
    // The reset color is the palette default of the feature being reset.
    const color = appliedFeatureColors(state)[feature.type];
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
    const svg = getCurrentSvg();
    if (targetFeatures.length === 0 || !svg) return false;

    const normalizedStrokeColor = String(strokeColor || '').trim();
    const normalizedStrokeWidth = normalizeStrokeWidthValue(strokeWidth);
    if (!normalizedStrokeColor && normalizedStrokeWidth === null) return false;

    let updatedCount = 0;
    targetFeatures.forEach((feature) => {
      const elements = getFeatureElements(svg, feature.svg_id);
      const needsUpdate = elements.some((element) => (
        (normalizedStrokeColor && !strokeColorAttributeMatches(element, normalizedStrokeColor))
        || (normalizedStrokeWidth !== null && !strokeWidthAttributeMatches(element, normalizedStrokeWidth))
      ));
      if (!needsUpdate) return;
      recordFeatureStrokeOverride(drawing, feature, {
        strokeColor: normalizedStrokeColor || null,
        strokeWidth: normalizedStrokeWidth,
        ...drawnFeatureStroke(elements)
      });
      updatedCount += 1;
      if (clickedFeature.value?.svg_id === feature.svg_id) {
        if (normalizedStrokeColor) clickedFeature.value.strokeColor = normalizedStrokeColor;
        if (normalizedStrokeWidth !== null) clickedFeature.value.strokeWidth = normalizedStrokeWidth;
      }
    });
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

    const overrideKey = targetLegendEntry.caption;
    const previousOverride = drawing.legendStrokeOverrides[overrideKey] || {};
    const nextOverride = { ...previousOverride };
    if (
      !hasOwn(previousOverride, 'strokeColor')
      && !hasOwn(previousOverride, 'strokeWidth')
    ) {
      Object.assign(nextOverride, drawnLegendRowStroke(svg, overrideKey));
    }
    if (normalizedStrokeColor) nextOverride.strokeColor = normalizedStrokeColor;
    if (normalizedStrokeWidth !== null) nextOverride.strokeWidth = normalizedStrokeWidth;
    const stateChanged =
      !colorsMatch(previousOverride.strokeColor, nextOverride.strokeColor)
      || normalizeStrokeWidthValue(previousOverride.strokeWidth)
        !== normalizeStrokeWidthValue(nextOverride.strokeWidth);
    if (stateChanged) drawing.legendStrokeOverrides[overrideKey] = nextOverride;
    return stateChanged;
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
    cancelResetColor: closeResetColorDialog,
    handleColorScopeChoice: colorScopeChoice,
    handleFeatureStyleScopeChoice: (...args) => (featureStyleScopeDialog.kind === 'stroke' ? strokeScopeChoice : colorScopeChoice)(...args),
    handleLegendNameCommit: colorAction(handleLegendNameCommit, { targetRules: false }),
    handleLegendRenameChoice: colorAction(handleLegendRenameChoice),
    renameLegendEntry: colorAction(renameLegendEntry),
    requestFeatureColorChange: colorAction(requestFeatureColorChange, { targetRules: false }),
    selectLegendNameOption: colorAction(selectLegendNameOption, { targetRules: false }),
    handleResetColorChoice: colorAction(handleResetColorChoice),
    applyColorToSelectedFeatures: colorAction(applyColorToSelectedFeatures),
    applyStrokeToSelectedFeatures: strokeAction(applyStrokeToSelectedFeatures),
    resetClickedFeatureFillColor: colorAction(resetClickedFeatureFillColor, { targetRules: false }),
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
