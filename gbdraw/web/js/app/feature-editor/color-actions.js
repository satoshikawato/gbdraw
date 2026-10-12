// @ts-check
/** @import { DrawingState, LegendEntry } from '../../state.js' */
/** @import { LegendRowReach, PythonLegendKey, RenderedFeatureId } from '../../services/legend-svg.js' */
/** @import { ColorRuleLegendGroup, SpecificColorRule, SpecificRuleBuild } from '../../services/specific-color-rules.js' */
/** @import { DrawnTrack } from '../../services/track-legend-rows.js' */
import { displayedLegendRowColors, draftLegendRowColors, draftLegendRows, LIVE_EDIT_DOMAINS, namedLegendCaption } from '../candidate-render.js';
import { reportRuleRunFailure } from '../rule-matching.js';
import { matchedRuleKeys, ruleKey, ruleMatcher, ruleMatchesFeature } from '../../services/rule-matchers.js';
import { appliedFeatureColorOf, appliedFeatureColors, resolveColorToHex } from '../../utils/color-utils.js';
import { getFeatureCaption, getFeatureColorRuleHash, getFeatureHashCandidates } from '../../services/feature-utils.js';
import { exactRegexValue } from '../../services/feature-selector.js';
import {
  drawnLegendRowStroke,
  legendEntryKey,
  legendRowFeatureIds,
  setsFeatureStroke
} from '../../services/legend-svg.js';
import { displayedLegendRowContext } from '../../services/feature-visibility.js';
import { colorRuleLegendGroups, ruleLegendCaptions } from '../../services/specific-color-rules.js';
import { isAutoFeatureUnderlay } from '../../services/feature-dom.js';
import { pythonDrawnAttribute } from '../../services/result-paint-bases.js';
import {
  featureOverrideKey,
  getFeatureOverride
} from '../../services/feature-override-identity.js';

/**
 * The rule owner's actions this owner reads (a whole-object port, R13).
 * @typedef {object} ColorActionsRuleActions
 * @property {(build: SpecificRuleBuild, label?: string, options?: Record<string, any>) => Promise<any>} commitSpecificRules
 *   The rule owner's one write of the rule table: how the next table follows from the table at the write's turn.
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
 * @property {(feature: Record<string, any>, label: string) => { feat: string | undefined, qual: string, val: string } | null} getLabelSpecificRule
 * @property {(key: PythonLegendKey) => Record<string, any>[]} getLegendRowRules
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
 * A feature row of the Features list, as this owner reads it: its type and its
 * rendered id (null for the row-less placeholder of the Legend popup's count),
 * and the label fields Python writes as text, beside the Python fields that rule
 * matching reads by name.
 * @typedef {Readonly<Record<string, unknown>> & { type?: string, svg_id?: string | null, product?: string, gene?: string, locus_tag?: string }} ColorFeature
 */
/**
 * A feature, or the feature popup's state that holds it as `feat`.
 * @typedef {ColorFeature & { feat?: ColorFeature | null }} ColorFeatureOrPopup
 */
/**
 * A color rule as the editor reads it (a draft rule may lack fields).
 * @typedef {Partial<SpecificColorRule>} ColorRule
 */

/**
 * A rename of a Legend row from the Legend editor or the feature popup.
 * `oldCaption` is the caption the source shows; `sourceKey` is its row's
 * Python key (the Legend editor's `legendEntryKey`, the popup's effective
 * caption), by which the rules the row draws are read. A renamed row's shown
 * caption is not its key (OV-294 residual, review M1 of 2b96350f).
 * @typedef {object} LegendRenameRequest
 * @property {'popup' | 'legend'} source Where the rename was typed.
 * @property {ColorFeature} [feat] The popup's feature.
 * @property {string} oldCaption
 * @property {PythonLegendKey} sourceKey
 * @property {string} newCaption
 * @property {ColorRuleLegendGroup} [joinGroup] The group picked from the popup's list, which the rename
 *   joins: `newCaption` is its name, and its color comes with it (PD-OI-061 revision 4).
 * @property {string} [currentColor]
 * @property {number} [siblingCount]
 * @property {ColorFeature[]} [features]
 * @property {'single' | 'group' | 'manual' | null} [sourceScope]
 * @property {'merge' | 'suffix' | null} [targetResolution]
 * @property {string} [finalCaption]
 * @property {string} [finalColor]
 */
/**
 * What a rename collides with: a listed or deleted Legend row, or a color-rule
 * group that draws no row. `caption` and `color` are what a Merge adopts.
 * @typedef {object} RenameTarget
 * @property {string} caption
 * @property {string} color
 * @property {ReadonlySet<string | undefined>} types The feature types it draws.
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
 * @property {ReadDrawnTracksPort} readDrawnTracks
 */

/**
 * The tracks of the request that drew the displayed Result (`drawnTracks`), which give a track's
 * Legend row its color (OV-484).
 * @typedef {() => readonly DrawnTrack[]} ReadDrawnTracksPort
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
  setDefaultColor,
  readDrawnTracks
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
  // This owner's rule writes are part of its caller's History step (a popup,
  // Features panel, or dialog action): they join it and keep its name.
  /** @type {ColorActionsRuleActions['commitSpecificRules']} */
  const commitRuleWrite = (build, label, options) => ruleActions.commitSpecificRules(build, label, { ...options, joinsOpenStep: true });
  /** @param {ColorFeature} feature */
  const getEffectiveLegendCaption = (feature) => effectiveLegendCaptions()(feature);
  /** @param {unknown} value */
  const normalizeCaption = (value) => String(value || '').trim();
  /** @param {unknown} value */
  const normalizeCaptionKey = (value) => normalizeCaption(value).toLowerCase();
  /** @param {unknown} value */
  const normalizeColor = (value) => String(value || '').trim().toLowerCase();
  /** @param {unknown} left @param {unknown} right */
  const captionsMatch = (left, right) => normalizeCaptionKey(left) === normalizeCaptionKey(right);
  /** @param {unknown} left @param {unknown} right */
  const colorsMatch = (left, right) => normalizeColor(left) === normalizeColor(right);
  /** @param {ColorRule | null | undefined} rule */
  const isHashSpecificRule = (rule) => String(rule?.qual || '').toLowerCase() === 'hash';
  /** @param {Readonly<Record<string, unknown>> | null | undefined} object @param {string} key */
  const hasOwn = (object, key) => Object.prototype.hasOwnProperty.call(object || {}, key);
  // The domains a Legend row rename shows: the row's structure, and the fill
  // and stroke an editor row takes under its new caption.
  const RENAMED_ROW_DOMAINS = [...LIVE_EDIT_DOMAINS.legendStructure, ...LIVE_EDIT_DOMAINS.legendFills, ...LIVE_EDIT_DOMAINS.strokes];

  // Runs `run` once the rules a color edit may add for the features in `args`,
  // the clicked feature and the scope dialog's feature (their hash and label
  // rules) are prepared with the saved ones.
  /**
   * @param {DrawingState} drawing
   * @param {readonly unknown[]} args
   * @param {() => any} run
   */
  const withTargetRules = (drawing, args, run) => {
    /** @type {Set<ColorFeature>} */
    const targets = new Set();
    // An argument that is not a feature (a color, a caption) has no `type`: it reads as undefined and adds nothing.
    /** @param {ColorFeature | null | undefined} feature */
    const add = (feature) => { if (feature?.type && feature?.svg_id) targets.add(feature); };
    args.forEach((arg) => Array.isArray(arg) ? arg.forEach(add) : add(/** @type {ColorFeature | null | undefined} */ (arg)));
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
  /**
   * @template {unknown[]} Args
   * @param {(drawing: DrawingState, ...args: Args) => unknown} action
   * @param {{ targetRules?: boolean }} [options]
   * @returns {(...args: Args) => unknown}
   */
  const colorAction = (action, { targetRules = true } = {}) => (...args) => {
    const drawing = state.activeDrawing();
    const run = async () => action(drawing, ...args);
    const prepareTargets = () => (targetRules ? withTargetRules(drawing, args, run) : run());
    return reportRuleRunFailure(state, 'evaluateRules', () => runWithRuleMatches(drawing.manualSpecificRules, prepareTargets));
  };
  /**
   * @template {unknown[]} Args
   * @param {(drawing: DrawingState, ...args: Args) => unknown} action
   */
  const strokeAction = (action) => colorAction(action, { targetRules: false });

  // The `hash` rule values that target a feature exactly.
  /** @param {ColorFeature} feature */
  const exactHashValues = (feature) => {
    const candidates = getFeatureHashCandidates(feature);
    const generationHash = getFeatureColorRuleHash(feature);
    const renderedId = candidates[candidates.length - 1] || '';
    // Python matches only the stable hash, which duplicate records share.
    return [generationHash, ...(renderedId !== generationHash ? [renderedId] : [])]
      .filter(Boolean).flatMap((candidate) => [candidate, exactRegexValue(candidate)]);
  };
  /** @param {unknown} type @param {unknown} value */
  const exactHashKey = (type, value) => JSON.stringify([type, value]);
  /** @param {ColorRule} rule */
  const ruleHashKey = (rule) => exactHashKey(rule.feat, String(rule?.val || '').trim());
  /** @param {ColorRule} rule @param {ColorFeature} feature */
  const hashRuleTargetsFeatureExactly = (rule, feature) => isHashSpecificRule(rule) && rule?.feat === feature?.type
    && exactHashValues(feature).includes(String(rule?.val || '').trim());
  // Whether a rule targets one of `features` exactly, read once per feature.
  /**
   * @param {readonly ColorFeature[]} features
   * @returns {(rule: ColorRule) => boolean}
   */
  const targetsSomeFeatureExactly = (features) => {
    const keys = new Set(features.flatMap((feature) => exactHashValues(feature)
      .map((value) => exactHashKey(feature.type, value))));
    return (rule) => isHashSpecificRule(rule) && keys.has(ruleHashKey(rule));
  };

  /** @param {unknown} value */
  const normalizeStrokeWidthValue = (value) => {
    if (value === null || value === undefined || value === '') return null;
    const numeric = Number(value);
    return Number.isFinite(numeric) && numeric >= 0 ? numeric : null;
  };

  /** @param {Element | null | undefined} element @param {string | null} color */
  const strokeColorAttributeMatches = (element, color) => {
    const current = element?.getAttribute?.('stroke');
    return color === null ? current === null : colorsMatch(current, color);
  };

  /** @param {Element | null | undefined} element @param {number | null} width */
  const strokeWidthAttributeMatches = (element, width) => {
    const current = element?.getAttribute?.('stroke-width');
    if (width === null) return current === null;
    const currentNumeric = Number(current);
    return Number.isFinite(currentNumeric) && currentNumeric === Number(width);
  };

  /**
   * @param {ColorFeatureOrPopup | null | undefined} featureLike
   * @param {unknown} [fallbackSvgId]
   */
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

  /**
   * @param {DrawingState} drawing
   * @param {ColorFeatureOrPopup | null | undefined} featureLike
   * @param {unknown} [fallbackSvgId]
   */
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
    const context = displayedLegendRowContext(state, drawing);
    const { pythonRows } = context;
    const rows = draftLegendRows({
      legendEntries: drawing.legendEntries.value, deletedLegendEntries: drawing.deletedLegendEntries.value,
      dormantLegendEntries: drawing.dormantLegendEntries.value, originalLegendOrder: originalLegendOrder.value || []
    });
    const draft = draftLegendRowColors({ ...context, paletteColors: appliedFeatureColors(state), drawnTracks: readDrawnTracks() });
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
  // OV-282 (D-26): the color a listed row shows, the palette's for a row the
  // palette colors; never the color the row recorded at the last Generate.
  // Without a row, nothing is derived (review L4).
  /** @param {DrawingState} drawing @param {LegendEntry | null} entry */
  const listedRowColor = (drawing, entry) => (entry ? displayedLegendRowColors(state, drawing, readDrawnTracks()).listed(entry) : undefined);
  /** @param {DrawingState} drawing @param {ColorFeature} feat @param {string} caption */
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
    const legendColor = listedRowColor(drawing, legendEntry);
    if (legendEntry && legendColor) {
      return {
        caption: legendEntry.caption,
        color: legendColor,
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

  /** @param {DrawingState} drawing @param {ColorFeature} feature @returns {ColorRule[]} */
  const exactHashRulesForFeature = (drawing, feature) => drawing.manualSpecificRules.filter(
    (/** @type {ColorRule} */ rule) => hashRuleTargetsFeatureExactly(rule, feature)
  );

  /** @param {ColorFeature} feature @param {unknown} color */
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

  /** @param {DrawingState} drawing @param {ColorFeature} feature @param {string} color @param {string} caption */
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
    return Boolean(legendEntry && colorsMatch(listedRowColor(drawing, legendEntry), color));
  };

  /** @param {Readonly<Record<string, unknown>> | null | undefined} store @param {unknown} caption */
  const findCaptionKey = (store, caption) => {
    if (!store) return null;
    return Object.keys(store).find((key) => captionsMatch(key, caption)) || null;
  };

  /** @param {Record<string, unknown> | null | undefined} store @param {string} oldCaption @param {string} newCaption */
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


  /** @param {DrawingState} drawing @param {string} oldCaption @param {string} newCaption */
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


  /** @param {string} oldCaption @param {string} newCaption @param {string | null} [color] */
  const syncOriginalLegendMetadataRename = (oldCaption, newCaption, color = null) => {
    if (!oldCaption || !newCaption || oldCaption === newCaption) return;

    const orderIdx = originalLegendOrder.value.findIndex((/** @type {string} */ caption) => captionsMatch(caption, oldCaption));
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


  /** @param {ColorFeature | null | undefined} feat @param {string} caption @param {string | null} [color] */
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
    legendRenameDialog.suffixCaption = '';
    legendRenameDialog.pendingRequest = null;
  });

  /** @param {DrawingState} drawing @param {ColorFeature | null | undefined} feat */
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

    const fallbackColor = appliedFeatureColorOf(state, feat.type) || '#cccccc';
    return resolveColorToHex(fallbackColor) || fallbackColor;
  };

  /** @param {unknown} caption */
  const getFeaturesForLegendCaption = (caption) => {
    const normalizedCaption = normalizeCaption(caption);
    if (!normalizedCaption) return [];
    const captionOf = effectiveLegendCaptions();
    return /** @type {ColorFeature[]} */ (extractedFeatures.value).filter((feat) => captionsMatch(captionOf(feat), normalizedCaption));
  };

  /** @param {DrawingState} drawing @param {string} caption @param {{ ignoreCaptions?: unknown }} [options] */
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
    drawing.manualSpecificRules.forEach((/** @type {ColorRule} */ rule) => {
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

  // Each feature's `hash` rule replaces the first rule that targets it
  // exactly, else goes before the first `hash` rule of its type that matches
  // it, else last. The rules are found by index, read once per feature (OV-193).
  /**
   * @param {readonly SpecificColorRule[]} current The rule table at the write's turn.
   * @param {readonly ColorFeature[]} features
   * @param {string} color
   * @param {string} caption
   * @param {{ preferLabelRules?: boolean }} [options]
   */
  const featureRuleCandidate = (current, features, color, caption, { preferLabelRules = false } = {}) => {
    /** @type {Partial<SpecificColorRule>[]} */
    let rules = [...current];
    const labelRule = preferLabelRules ? getSafeLabelSpecificRule(current, features, caption) : null;
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
    /** @param {ColorRule} rule */
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

  /** @param {ColorFeature | null | undefined} feature */
  const featureIdentityKey = (feature) => String(
    feature?.id || `${feature?.type || ''}:${feature?.svg_id || ''}`
  );

  /** @param {readonly SpecificColorRule[]} current @param {readonly ColorFeature[]} features @param {string} label */
  const getSafeLabelSpecificRule = (current, features, label) => {
    if (typeof getLabelSpecificRule !== 'function') return null;
    const candidates = features.map((feature) => getLabelSpecificRule(feature, label));
    if (candidates.some((rule) => !rule)) return null;

    // The check above returned for a Result without a rule, so every candidate is one.
    const labelRules = /** @type {Array<NonNullable<typeof candidates[number]>>} */ (candidates);
    const first = labelRules[0];
    const sameSelector = labelRules.every(
      (rule) =>
        rule.feat === first.feat &&
        String(rule.qual || '').toLowerCase() === String(first.qual || '').toLowerCase() &&
        rule.val === first.val
    );
    if (!sameSelector) return null;

    const selectedKeys = new Set(features.map(featureIdentityKey));
    const matchedFeatures = /** @type {ColorFeature[]} */ (extractedFeatures.value).filter((feature) => ruleMatchesFeature(feature, first));
    if (
      matchedFeatures.length !== selectedKeys.size ||
      matchedFeatures.some((feature) => !selectedKeys.has(featureIdentityKey(feature)))
    ) {
      return null;
    }

    /** @type {ColorFeature[]} */
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
    const hasPrecedenceConflict = current.some((existing) => {
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
    const legendEntry = drawing.legendEntries.value.find((/** @type {LegendEntry} */ entry) => captionsMatch(entry?.caption, oldCaption));
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
      (/** @type {string} */ caption) => captionsMatch(caption, legendEntry.originalCaption || legendEntry.caption)
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
  /** @param {DrawingState} drawing @param {LegendRenameRequest} request */
  const applyLegendRenameRequest = (drawing, request) => withTargetRules(drawing, [], () => applyPreparedLegendRename(drawing, request));
  /** @param {DrawingState} drawing @param {LegendRenameRequest} request */
  const applyPreparedLegendRename = async (drawing, request) => {
    const oldCaption = normalizeCaption(request.oldCaption);
    const caption = normalizeCaption(request.finalCaption || request.newCaption);
    const color = resolveColorToHex(request.finalColor || request.currentColor) || '#cccccc';
    if (!caption) return false;
    const features = (request.features || []).filter(Boolean);
    const sourceRules = getLegendRowRules(request.sourceKey);
    if (sourceRules.length || features.length) {
      /** @type {SpecificRuleBuild} */
      const rules = request.sourceScope === 'group' && sourceRules.length
        ? (current, revise) => current.map(rule => (sourceRules.includes(rule) ? revise(rule, { cap: caption, color }) : rule))
        : (current) => featureRuleCandidate(current, features, color, caption);
      const selected = new Set(features.map(feature => feature.svg_id));
      const retireOld = features.length > 0 && getFeaturesForLegendCaption(oldCaption)
        .every(feature => selected.has(feature.svg_id));
      const oldEntry = findLegendEntryByCaption(drawing, oldCaption);
      return commitRuleWrite(rules, 'Rename legend item', {
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

  /** @param {LegendRenameRequest} request @param {number} siblingCount */
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
   * @param {LegendRenameRequest} request
   * @param {RenameTarget} target
   * @param {boolean} mergeAvailable
   * @param {PythonLegendKey | ''} deletedTargetKey The key of a deleted target row, else ''.
   * @param {string} suffixCaption The name "Keep both" gives the renamed row.
   */
  const openLegendRenameTargetDialog = (request, target, mergeAvailable, deletedTargetKey, suffixCaption) => {
    legendRenameDialog.show = true;
    legendRenameDialog.mode = 'target';
    legendRenameDialog.oldCaption = request.oldCaption;
    legendRenameDialog.newCaption = request.newCaption;
    legendRenameDialog.targetCaption = target.caption;
    legendRenameDialog.targetColor = target.color;
    legendRenameDialog.currentColor = request.currentColor || '';
    legendRenameDialog.siblingCount = request.siblingCount || 0;
    legendRenameDialog.mergeAvailable = mergeAvailable;
    legendRenameDialog.deletedTargetKey = deletedTargetKey;
    legendRenameDialog.suffixCaption = suffixCaption;
    legendRenameDialog.pendingRequest = request;
  };

  /** @param {DrawingState} drawing @param {LegendRenameRequest} request */
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
      (request.feat ? getCurrentFeatureFillColor(drawing, request.feat) : resolveColorToHex(listedRowColor(drawing, findLegendEntryByCaption(drawing, oldCaption)))) ||
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

    // D-06 (PD-OI-061): a rename onto another name of a different color asks
    // Merge, Keep both (suffix), or Cancel, with or without features.
    // Revision 4 (OV-580 to OV-582): the name of a color rule is such a name,
    // whether its rules draw a row or not, and a pick from the popup's list
    // (`joinGroup`) is its own Merge where Merge is allowed.
    // R15-3 (OV-285): a rename onto a deleted row's caption always asks, as onto
    // a listed row; its Merge is a Restore of that row, then the merge. The
    // caption only detects the conflict: the row is its Python key (an editor
    // row's own caption), which the composition root's Restore acts on.
    const group = request.joinGroup || null;
    // A picked group's row is the row the allocation draws for it (N-06), never the row of its name.
    const groupRowKey = group && ruleLegendCaptions(displayedLegendRowContext(state, drawing))({ cap: group.name, color: group.color });
    /** @param {LegendEntry[]} entries */
    const rowOf = (entries) => (group
      ? entries.find((entry) => legendEntryKey(entry) === groupRowKey) || null : entryByCaption(entries, newCaption));
    const listedTarget = rowOf(drawing.legendEntries.value);
    const deletedTarget = listedTarget ? null : rowOf(drawing.deletedLegendEntries.value);
    const targetRow = listedTarget || deletedTarget;
    const targetGroup = group || (targetRow ? null
      : colorRuleLegendGroups(drawing.manualSpecificRules).find((candidate) => captionsMatch(candidate.name, newCaption)) || null);
    /** @param {readonly ColorFeature[]} rowFeatures */
    const featureTypes = (rowFeatures) => new Set(rowFeatures.map((feature) => feature?.type));
    // The name and the color a Merge adopts, and the feature types the target
    // draws: a group's own (a picked one also when its row is listed), else
    // the row's caption, the color it shows (OV-282: the palette's for a row
    // the palette colors), and the types of the features that take its caption.
    /** @type {RenameTarget | null} */
    const target = targetRow && captionsMatch(targetRow.caption, oldCaption) ? null : targetRow ? {
      caption: group ? group.name : targetRow.caption,
      color: group ? group.color : String((deletedTarget
        ? displayedLegendRowColors(state, drawing, readDrawnTracks()).deleted(deletedTarget) : listedRowColor(drawing, listedTarget)) || ''),
      types: featureTypes(getFeaturesForLegendCaption(targetRow.caption))
    } : targetGroup && { caption: targetGroup.name, color: targetGroup.color, types: new Set(targetGroup.types) };
    const featureOrRuleRename = features.length > 0 || getLegendRowRules(request.sourceKey).length > 0;
    // OV-62 (PD-OI-061 amended): two rows merge only when each draws features of
    // one type and the type is the same. A row without features (GC content,
    // GC skew), rows of different types, and a row that spans several types
    // offer Keep both and Cancel only. The types are those of the features
    // that take each row's caption, as the live editor knows them: a generated
    // row such as `other proteins` has none live, so it is not merged into. A
    // color-rule group that draws no row has the types of its rules.
    const sourceTypes = featureTypes(features);
    const mergeAllowed = target !== null && sourceTypes.size === 1 && target.types.size === 1
      && [...sourceTypes][0] === [...target.types][0];
    const sameColor = target !== null && colorsMatch(target.color, currentColor);
    const resolution = request.targetResolution || (group && mergeAllowed && !deletedTarget ? 'merge' : null);

    // A chosen Merge adopts the target's color, also once a Restore made a
    // deleted target a listed one.
    if (featureOrRuleRename && resolution !== 'merge' && (!target || (!deletedTarget && mergeAllowed && sameColor))) {
      await applyLegendRenameRequest(drawing, { ...request, currentColor, features,
        finalCaption: newCaption, finalColor: currentColor });
      clearLegendRenameDialog(drawing);
      return;
    }

    if (target && (deletedTarget || !mergeAllowed || !sameColor)) {
      if (!resolution) {
        openLegendRenameTargetDialog(
          { ...request, currentColor, features },
          target,
          mergeAllowed,
          deletedTarget ? legendEntryKey(deletedTarget) : '',
          getUniqueLegendCaption(drawing, newCaption, { ignoreCaptions: [oldCaption] })
        );
        return;
      }

      // A deleted row takes the merge only once the Restore returned it.
      if (resolution === 'merge' && (!mergeAllowed || deletedTarget)) {
        clearLegendRenameDialog(drawing, { restoreInput: true });
        return;
      }

      await applyLegendRenameRequest(drawing, {
        ...request,
        currentColor,
        features,
        ...(resolution === 'merge' ? { finalCaption: target.caption, finalColor: target.color }
          : { finalCaption: getUniqueLegendCaption(drawing, newCaption, { ignoreCaptions: [oldCaption] }), finalColor: currentColor })
      });
      clearLegendRenameDialog(drawing);
      return;
    }

    await applyLegendRenameRequest(drawing, {
      ...request,
      currentColor,
      features,
      finalCaption: target && sameColor ? target.caption : newCaption,
      finalColor: target && sameColor ? target.color : currentColor
    });
    clearLegendRenameDialog(drawing);
  };

  /**
   * @param {DrawingState} drawing
   * @param {readonly ColorFeature[]} features
   * @param {string} caption
   * @param {string} color
   * @param {{ preferLabelRules?: boolean }} [options]
   */
  const applyColorToFeatureGroup = async (drawing, features, caption, color, options = {}) => {
    if (!features?.length || !normalizeCaption(caption)) return false;
    const existingEntry = findLegendEntryByCaption(drawing, caption);
    const contributors = getFeaturesForLegendCaption(caption);
    const captionOf = effectiveLegendCaptions();
    const selectedIds = new Set(features.map(feature => feature.svg_id));
    const replacesExistingGroup = existingEntry && contributors.length > 0
      && features.every(feature => captionsMatch(captionOf(feature), caption))
      && contributors.every(feature => selectedIds.has(feature.svg_id));
    return commitRuleWrite((current) => featureRuleCandidate(current, features, color, normalizeCaption(caption), options), 'Change feature color', {
      drawing,
      previousLegendIntents: replacesExistingGroup
        ? [{ caption: existingEntry.caption, color: existingEntry.color }] : [],
      afterCommit: (/** @type {Array<{ caption: string, color: string }>} */ intents) => {
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
    if (!features.every((feature) => feature?.type === type) || getLegendRowRules(/** @type {PythonLegendKey} */ (type)).length > 0) return null;
    const matches = ruleMatcher(drawing.manualSpecificRules);
    return features.every((feature) => matches.matchesAny(feature) === false
      && !getFeatureOverride(drawing.featureColorOverrides, feature)) ? type : null;
  };

  /** @param {DrawingState} drawing @param {PythonLegendKey} caption @param {string} color */
  const applyColorToLegendSpecificRules = async (drawing, caption, color) => {
    const rowRules = getLegendRowRules(caption);
    const specificRules = rowRules.filter(rule => !isHashSpecificRule(rule));
    if (!specificRules.length) return false;
    const specificMatches = ruleMatcher(specificRules);
    const coveredExactly = targetsSomeFeatureExactly(
      /** @type {ColorFeature[]} */ (extractedFeatures.value).filter(feature => specificMatches.matchesAny(feature) === true)
    );
    const rowRuleSet = new Set(rowRules);
    return commitRuleWrite((current, revise) => current.filter(rule => !(rowRuleSet.has(rule) && coveredExactly(rule)))
      .map(rule => (rowRuleSet.has(rule) ? revise(rule, { color }) : rule)), 'Change legend color', { drawing });
  };

  /**
   * @param {DrawingState} drawing
   * @param {ColorFeature} feat
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

  /** @param {DrawingState} drawing @param {string} color */
  const updateClickedFeatureColor = async (drawing, color) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!clickedFeature.value) return;
    const feat = clickedFeature.value.feat;
    if (!feat) return;
    const customName = normalizeCaption(clickedFeature.value.legendName);
    await requestFeatureColorChange(drawing, feat, color, customName, { closePopupOnDialog: true });
  };

  // The popup's Legend name, typed or picked from the list of color-rule groups.
  /** @param {DrawingState} drawing @param {ColorRuleLegendGroup} [joinGroup] The picked group. */
  const commitLegendName = async (drawing, joinGroup) => {
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
      sourceKey: /** @type {PythonLegendKey} */ (currentCaption),
      newCaption: requestedCaption,
      joinGroup,
      currentColor: getCurrentFeatureFillColor(drawing, feat),
      sourceScope: null
    });
  };

  // A pick from the list under the popup's Legend name joins the picked group:
  // the feature takes its name and its color (PD-OI-061 revision 4, OV-580).
  /** @param {DrawingState} drawing @param {ColorRuleLegendGroup | null | undefined} group */
  const joinLegendGroup = async (drawing, group) => {
    if (!clickedFeature.value || !group) return;
    clickedFeature.value.legendName = group.name;
    await commitLegendName(drawing, group);
  };

  /** @param {DrawingState} drawing @param {string} choice */
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

  /** @param {DrawingState} drawing @param {number} idx @param {string} newCaption */
  const renameLegendEntry = async (drawing, idx, newCaption) => {
    const entry = drawing.legendEntries.value[idx];
    if (!entry) return;

    const requestedCaption = normalizeCaption(newCaption);
    if (!requestedCaption || requestedCaption === normalizeCaption(entry.caption)) {
      drawing.legendEntries.value = [...drawing.legendEntries.value];
      return;
    }

    const shownColor = listedRowColor(drawing, entry);
    await continueLegendRenameRequest(drawing, {
      source: 'legend',
      oldCaption: normalizeCaption(entry.caption),
      sourceKey: legendEntryKey(entry),
      newCaption: requestedCaption,
      currentColor: resolveColorToHex(shownColor) || shownColor || '#cccccc',
      features: getFeaturesForLegendCaption(entry.caption),
      sourceScope: getFeaturesForLegendCaption(entry.caption).length > 0 ? 'group' : 'manual'
    });
  };

  /** @param {DrawingState} drawing @param {string} choice */
  const handleColorScopeChoice = async (drawing, choice) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const { feat, color, matchingRule, legendName, existingCaptionColor, defaultColorType } = featureStyleScopeDialog;
    if (choice === 'cancel' || !feat || !color) {
      clearFeatureStyleScopeDialog();
      return;
    }

    if (choice === 'rule') {
      if (matchingRule) await commitRuleWrite((current, revise) => current.map(rule => (rule === matchingRule
        ? revise(rule, { color }) : rule)), 'Change specific color rule', { drawing });
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
      } else if (!(await applyColorToLegendSpecificRules(drawing, /** @type {PythonLegendKey} */ (targetLegendName), color))) {
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
  /** @param {DrawingState} drawing @param {string | null | undefined} strokeColor @param {unknown} strokeWidth */
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

  /** @param {DrawingState} drawing @param {string | null | undefined} strokeColor @param {unknown} strokeWidth */
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

  /** @param {ColorFeatureOrPopup | null | undefined} featureLike */
  const getFeatureStrokeColorValue = (featureLike) => {
    const drawing = state.activeDrawing();
    const feature = featureLike?.feat || featureLike;
    const override = getFeatureOverride(drawing.featureStrokeOverrides, feature);
    return override && hasOwn(override, 'strokeColor')
      ? override.strokeColor
      : null;
  };

  /** @param {DrawingState} drawing @param {string | null} value */
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

  /** @param {DrawingState} drawing @param {unknown} value */
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

    const defaultColor = appliedFeatureColorOf(state, feat.type);
    if (!defaultColor) {
      console.warn('No default color found for feature type:', feat.type);
      return;
    }

    const caption = getFeatureCaption(feat);

    const siblings = /** @type {ColorFeature[]} */ (extractedFeatures.value).filter(
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

  /** @param {DrawingState} drawing @param {string} choice */
  const handleResetColorChoice = async (drawing, choice) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    await doResetFillColor(drawing, choice);
    closeResetColorDialog();
  };

  /** @param {DrawingState} drawing @param {string} choice */
  const doResetFillColor = async (drawing, choice) => {
    const feature = clickedFeature.value?.feat;
    if (!feature || choice === 'cancel') return false;
    const caption = getEffectiveLegendCaption(feature);
    // The reset color is the palette default of the feature being reset.
    const color = appliedFeatureColorOf(state, feature.type);
    // A type without a palette color has nothing to reset to.
    if (choice === 'this_with_legend') return color ? setFeatureColor(drawing, feature, color, caption) : false;
    /** @type {SpecificRuleBuild} */
    const reset = (current) => {
      /** @type {Partial<SpecificColorRule>[]} */
      const rules = current.filter(rule => (choice === 'all'
        ? rule.cap !== caption : !hashRuleTargetsFeatureExactly(rule, feature)));
      if (choice === 'this' && rules.some(rule => ruleMatchesFeature(feature, rule))) {
        const qualifier = getFeatureQualifier(feature);
        if (qualifier && color) rules.push({ feat: feature.type, ...qualifier, color,
          cap: feature.type === 'CDS' ? 'other proteins' : `other ${feature.type}s` });
      }
      return rules;
    };
    const applied = await commitRuleWrite(reset, 'Reset feature color', { drawing });
    if (applied) clickedFeature.value = null;
    return applied;
  };

  /** @param {unknown} features @returns {Array<ColorFeature & { svg_id: string }>} */
  const uniqueFeaturesBySvgId = (features) => {
    const seen = new Set();
    // The filter keeps only features with a non-empty svg_id.
    return /** @type {Array<ColorFeature & { svg_id: string }>} */ ((Array.isArray(features) ? features : []).filter((/** @type {ColorFeature | null | undefined} */ feature) => {
      const svgId = String(feature?.svg_id || '').trim();
      if (!svgId || seen.has(svgId)) return false;
      seen.add(svgId);
      return true;
    }));
  };

  /** @param {DrawingState} drawing @param {unknown} features @param {string} color @param {unknown} caption */
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

  /** @param {DrawingState} drawing @param {unknown} features @param {string} strokeColor @param {unknown} strokeWidth */
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

  /** @param {DrawingState} drawing @param {string} caption @param {string} strokeColor @param {unknown} strokeWidth */
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

  /** @param {DrawingState} drawing @param {string} choice */
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
      targetFeatures = /** @type {ColorFeature[]} */ (extractedFeatures.value).filter((candidate) => (
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

  /** @param {DrawingState} drawing @param {ColorFeature} feature @param {string | null} value @param {string | null} [customCaption] */
  const setFeatureColorValue = async (drawing, feature, value, customCaption = null) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!feature) return false;
    if (value === null) {
      return commitRuleWrite((current) => current.filter(rule => !hashRuleTargetsFeatureExactly(rule, feature)), 'Reset feature color', { drawing });
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
    /** @param {[choice: string]} args */
    handleFeatureStyleScopeChoice: (...args) => (featureStyleScopeDialog.kind === 'stroke' ? strokeScopeChoice : colorScopeChoice)(...args),
    handleLegendNameCommit: colorAction((drawing) => commitLegendName(drawing), { targetRules: false }),
    handleLegendRenameChoice: colorAction(handleLegendRenameChoice),
    renameLegendEntry: colorAction(renameLegendEntry),
    requestFeatureColorChange: colorAction(requestFeatureColorChange, { targetRules: false }),
    joinLegendGroup: colorAction(joinLegendGroup, { targetRules: false }),
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
