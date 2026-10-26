// @ts-check
/** @import { DrawingState } from '../../state.js' */
/** @import { KeptLegendRow, LegendCaption, LegendRow } from '../../services/drawing/paint.js' */
/** @import { LegendRowReach, PythonLegendKey, RenderedFeatureId } from '../../services/legend-svg.js' */
/** @import { LegendRowName } from '../../services/legend-row-names.js' */
/** @import { ColorRuleLegendGroup, SpecificColorRule, SpecificRuleBuild } from '../../services/specific-color-rules.js' */
/** @import { DrawnTrack } from '../../services/track-legend-rows.js' */
/** @typedef {import('../rule-matching.js').GroupColorEdit['scope']} RuleSynthesisScope */
import { displayedLegendRowColors, draftLegendRowColors, draftLegendRows, LIVE_EDIT_DOMAINS, namedLegendCaption } from '../candidate-render.js';
import { reportRuleRunFailure } from '../rule-matching.js';
import { ruleMatcher, ruleMatchesFeature } from '../../services/rule-matchers.js';
import { appliedFeatureColorOf, appliedFeatureColors, resolveColorToHex } from '../../utils/color-utils.js';
import { getFeatureCaption, getFeatureColorRuleHash, getFeatureHashCandidates } from '../../services/feature-utils.js';
import { exactRegexValue, literalPatternText, namedFeatureInstance } from '../../services/feature-selector.js';
import {
  drawnLegendRowStroke,
  legendRowFeatureIds,
  legendRowKeys,
  setsFeatureStroke
} from '../../services/legend-svg.js';
import {
  fillRecordOf,
  legendEntryKey,
  renameListedLegendRow,
  replaceLegendRows,
  setFeatureStroke,
  setLegendColor,
  setLegendStroke
} from '../../services/drawing/paint.js';
import { displayedLegendRowContext } from '../../services/feature-visibility.js';
import { colorRuleGroupTakes, colorRuleLegendGroups, groupLegendRowKeys } from '../../services/specific-color-rules.js';
import { isAutoFeatureUnderlay } from '../../services/feature-dom.js';
import { pythonDrawnAttribute } from '../../services/result-paint-bases.js';
import {
  featureOverrideKey,
  getFeatureOverride
} from '../../services/feature-override-identity.js';
import { msg } from '../../utils/i18n.js';

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
 * @property {() => (feature: Record<string, any>) => PythonLegendKey} effectiveLegendCaptions
 *   The rule owner's reader of a feature's Legend row, by the row's key (never a name to write).
 *   The Legend item of a feature, read for a pass over many features.
 * @property {(feature: Record<string, any>) => string} getIndividualFeatureLabel
 * @property {(feature: Record<string, any>) => { qual: string, val: string } | null} getFeatureQualifier
 * @property {(key: PythonLegendKey) => Record<string, any>[]} getLegendRowRules
 * @property {(features: readonly ColorFeature[], color: string, caption: string, scope: RuleSynthesisScope) => SpecificRuleBuild} groupColorRules
 *   The next rule table of a color edit: Python's rules for coloring exactly `features` (OV-244).
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
 * @property {'merge' | 'keep-both' | null} [targetResolution]
 * @property {string} [finalCaption]
 * @property {string} [finalColor]
 * @property {PythonLegendKey} [keepBothAfter] Keep both of features split off their row: the key of the
 *   target row the new row is placed right after (D-29).
 */
/**
 * What a rename collides with: a listed or deleted Legend row, or a color-rule
 * group that draws no row. `caption` and `color` are what a Merge adopts.
 * @typedef {object} RenameTarget
 * @property {PythonLegendKey | null} key The target row's key; null for a color rule's name that no row shows.
 * @property {string} caption The name a Merge adopts.
 * @property {string} color
 * @property {ReadonlySet<string | undefined>} types The feature types it draws.
 * @property {boolean} ruleOnly A color rule's name that no Legend row shows.
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
 * @property {LegendRowIsNamedPort} [legendRowIsNamed]
 * @property {LegendRowNamePort} [legendRowName]
 * @property {FindLegendRowNamedPort} [findLegendRowNamed]
 * @property {KeepBothNamesPort} [keepBothNamesOf]
 * @property {LegendEntryNamePort} [legendEntryName]
 * @property {WriteLegendRowColorPort} writeLegendRowColor The Legend owner's writer of a row's own color: a type's
 *   default color releases the row's own color through it.
 */

/**
 * The Legend color of the rows of some keys (the Legend owner's `writeLegendRowColor`, D-35).
 * @typedef {(drawing: DrawingState, keys: readonly PythonLegendKey[], color: string, options?: { inherit?: readonly PythonLegendKey[] }) => void} WriteLegendRowColorPort
 */
/**
 * Whether `typed` names the Legend row that shows the key `caption`: its key, its name, or the name
 * the Legend draws for it (the Legend owner's reading of the displayed Result). Typing one of them
 * into the row's name field is no rename.
 * @typedef {(caption: string, typed: string) => boolean} LegendRowIsNamedPort
 */
/**
 * The name a caption stands for (the Legend owner's `legendRowName`): the name of the row it is the
 * key of, else the caption itself. A caption read as a row is written as this name (OV-640).
 * @typedef {(caption: string) => LegendRowName} LegendRowNamePort
 */

/**
 * The name of a row of the drawing's Legend lists (the Legend owner's `legendEntryName`).
 * @typedef {(entry: KeptLegendRow) => LegendRowName} LegendEntryNamePort
 */
/**
 * The key of the Legend row a typed name names (the Legend owner's `findLegendRowNamed`, PD-OI-061
 * scenario revision 5): by the name the Legend draws for it, else by its name, the row in `color`
 * first; null when no listed or deleted row has the name.
 * @typedef {(typed: string, color: string) => PythonLegendKey | null} FindLegendRowNamedPort
 */
/**
 * The names "Keep both" gives (the Legend owner's `keepBothNamesOf`, PD-OI-061 scenario revision 5):
 * the name the renamed row is drawn with, and the target row's drawn name after the rename when the
 * target is drawn and that name changes, else null. `renamedKey` is the row renamed in place; null
 * for features split off their row, which are placed right after `targetKey` (D-29) in `color`.
 * `restores`: Keep both restores the deleted target `deletedKey` first, as the two rows have one
 * swatch and so are one row (D-28).
 * @typedef {(rename: { name: string, renamedKey: PythonLegendKey | null, targetKey: PythonLegendKey | null, color: string,
 *   deletedKey?: PythonLegendKey | null }) => { keepBothName: string, targetNameAfter: string | null, restores: boolean }} KeepBothNamesPort
 */
/**
 * The tracks of the request that drew the displayed Result (`drawnTracks`), which give a track's
 * Legend row its color (OV-484).
 * @typedef {() => readonly DrawnTrack[]} ReadDrawnTracksPort
 */

// Why the Legend Name Conflict dialog asks; a rule's name in the row's own color: its type (OV-649).
/** @param {{ name: string, ruleOnly: boolean, deleted: boolean, mergeAvailable: boolean, sameColor: boolean }} conflict */
const legendNameConflictReason = ({ name, ruleOnly, deleted, mergeAvailable, sameColor }) => {
  const named = `"${name}"`;
  if (ruleOnly) return `${named} is the name of a color rule ${sameColor ? 'for another feature type' : 'with another color'}.`;
  const row = deleted ? 'is a deleted Legend item' : 'already exists';
  if (!mergeAvailable) return `${named} ${row} and draws a different feature type, so the two cannot be merged.`;
  return deleted ? `${named} ${row}.` : `${named} ${row} with another color.`;
};

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
  readDrawnTracks,
  legendRowIsNamed = (caption, typed) => caption === typed,
  legendRowName = (caption) => /** @type {LegendRowName} */ (String(caption ?? '')),
  legendEntryName = (entry) => legendRowName(entry.caption),
  // Without a Legend owner (a test of this owner), a row is named by its
  // caption, the one in `color` first.
  findLegendRowNamed = (typed, color) => {
    const drawing = state.activeDrawing();
    const named = [...drawing.legendEntries.value, ...drawing.deletedLegendEntries.value]
      .filter((entry) => String(entry.caption).toLowerCase() === typed.toLowerCase());
    const wanted = resolveColorToHex(String(color || '')) || String(color || '').trim().toLowerCase();
    const found = named.find((entry) => (
      (resolveColorToHex(String(entry.color || '')) || String(entry.color || '').trim().toLowerCase()) === wanted)) || named[0];
    return found ? legendEntryKey(found) : null;
  },
  keepBothNamesOf = ({ name }) => ({ keepBothName: name, targetNameAfter: null, restores: false }),
  writeLegendRowColor
}) => {
  const {
    extractedFeatures,
    svgContainer,
    clickedFeature,
    featureStyleScopeDialog,
    resetColorDialog,
    legendRenameDialog,
    originalLegendOrder
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
    getLegendRowRules,
    groupColorRules,
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
    /** @type {ColorRule[]} */
    const candidates = [...drawing.manualSpecificRules];
    targets.forEach((feature) => {
      const hash = getFeatureQualifier(feature);
      if (hash) candidates.push({ feat: feature.type, ...hash });
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

  // The `hash` rule values that name a feature by its rendered ID, as rules
  // older than source hashes did; Python matches none of them.
  /** @param {ColorFeature} feature */
  const renderedIdValues = (feature) => {
    const candidates = getFeatureHashCandidates(feature);
    const renderedId = candidates[candidates.length - 1] || '';
    return renderedId && renderedId !== getFeatureColorRuleHash(feature) ? [renderedId, exactRegexValue(renderedId)] : [];
  };
  // A `hash` rule targets a feature exactly when its value is one literal
  // text (as written or as an exact pattern) that names the feature: its
  // feature instance, which Python composes and matches whole (PD-OI-069
  // revision 2; the rule's prepared match says whether it is this feature's),
  // or its stable hash, which older rules name it by and which the other
  // copies of the feature share.
  /** @param {ColorRule} rule @param {ColorFeature} feature */
  const hashRuleTargetsFeatureExactly = (rule, feature) => {
    if (!isHashSpecificRule(rule) || rule?.feat !== feature?.type) return false;
    const value = String(rule?.val || '').trim();
    if (namedFeatureInstance(rule)) return ruleMatchesFeature(feature, rule) === true;
    const hash = getFeatureColorRuleHash(feature);
    return (Boolean(hash) && literalPatternText(value) === hash) || renderedIdValues(feature).includes(value);
  };
  // Whether a rule targets one of `features` exactly.
  /**
   * @param {readonly ColorFeature[]} features
   * @returns {(rule: ColorRule) => boolean}
   */
  const targetsSomeFeatureExactly = (features) => (rule) => isHashSpecificRule(rule)
    && features.some((feature) => hashRuleTargetsFeatureExactly(rule, feature));

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
      setFeatureStroke(drawing.slices.paint, key, next);
    }
  };

  /**
   * @param {DrawingState} drawing
   * @param {ColorFeatureOrPopup | null | undefined} featureLike
   * @param {unknown} [fallbackSvgId]
   */
  const clearFeatureStrokeOverride = (drawing, featureLike, fallbackSvgId = '') => {
    const key = featureStrokeKey(featureLike, fallbackSvgId);
    if (key) setFeatureStroke(drawing.slices.paint, key, null);
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

  // The listed row of a key (a feature's Legend item is its row's key).
  /** @param {DrawingState} drawing @param {unknown} key @returns {LegendRow | null} */
  const findLegendRow = (drawing, key) => {
    const text = normalizeCaption(key);
    if (!text) return null;
    const entries = drawing.legendEntries.value;
    return entries.find((entry) => legendEntryKey(entry) === text) || null;
  };
  // The listed row a name names, by the key the Legend owner finds for it (its
  // drawn name, else its name, the row in `color` first; review L-3).
  /** @param {DrawingState} drawing @param {LegendRowName} name @param {string} color @returns {LegendRow | null} */
  const findLegendRowOfName = (drawing, name, color) => findLegendRow(drawing, findLegendRowNamed(String(name), color));
  // The row of a style scope's Legend name, which is the feature's row (a key)
  // or a name typed in the popup (`popupLegendCaption`): that row, else the
  // first row of that name.
  /** @param {DrawingState} drawing @param {string} legendName @returns {LegendRow | null} */
  const findScopeLegendRow = (drawing, legendName) => findLegendRow(drawing, legendName)
    || findLegendRowOfName(drawing, legendRowName(legendName), '');
  // OV-282 (D-26): the color a listed row shows, the palette's for a row the
  // palette colors; never the color the row recorded at the last Generate.
  // Without a row, nothing is derived (review L4).
  /** @param {DrawingState} drawing @param {LegendRow | null} entry */
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

    const legendEntry = findScopeLegendRow(drawing, caption);
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
    featureStyleScopeDialog.displayLabelOffered = false;
    featureStyleScopeDialog.annotationLabelOffered = false;
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
    const displayLabelSiblings = displayLabel ? findFeaturesWithSameDisplayedLabel(feat, displayLabel) : [];
    const displayLabelSiblingCount = displayLabelSiblings.length;
    const annotationLabel = normalizeCaption(getIndividualFeatureLabel(feat));
    const annotationLabelSiblings = annotationLabel ? findFeaturesWithSameIndividualLabel(feat, annotationLabel) : [];
    const annotationLabelSiblingCount = annotationLabelSiblings.length;
    // OV-648: a label scope is offered only when it selects other features
    // than the Legend item scope (and, for the annotation label, than the
    // display label scope); names are not compared, as a row key is no label
    // and a label that repeats a row's name may select other features.
    /** @param {readonly ColorFeature[]} features */
    const selection = (features) => [...new Set(features.map((feature) => feature?.svg_id))].sort().join('\u0000');
    const legendSelection = selection(siblings);
    const displayLabelOffered = displayLabelSiblingCount > 0 && selection(displayLabelSiblings) !== legendSelection;
    const annotationSelection = selection(annotationLabelSiblings);
    const annotationLabelOffered = annotationLabelSiblingCount > 0 && annotationSelection !== legendSelection
      && annotationSelection !== selection(displayLabelSiblings);

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
      displayLabelOffered,
      annotationLabelOffered,
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
   *   replacedDefaultColor?: string | null
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
    replacedDefaultColor = null
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
    featureStyleScopeDialog.displayLabelOffered = scope.displayLabelOffered;
    featureStyleScopeDialog.annotationLabelOffered = scope.annotationLabelOffered;
    featureStyleScopeDialog.existingCaptionRule = existingCaption?.rule || null;
    featureStyleScopeDialog.existingCaptionColor = existingCaption?.color || null;
    featureStyleScopeDialog.defaultColorType = defaultColorType;
    featureStyleScopeDialog.replacedDefaultColor = replacedDefaultColor;
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
    const override = fillRecordOf(getFeatureOverride(drawing.featureColorOverrides, feature));
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
    const legendEntry = findLegendRowOfName(drawing, /** @type {LegendRowName} */ (caption), color);
    return Boolean(legendEntry && colorsMatch(listedRowColor(drawing, legendEntry), color));
  };

  // The popup shows the fill an edit gave its feature; scope dialogs leave it open (OV-1278).
  /** @param {ColorFeature | null | undefined} feat @param {string} color */
  const showClickedFeatureColor = (feat, color) => {
    if (feat && clickedFeature.value?.svg_id === feat.svg_id) clickedFeature.value.color = color;
  };

  /** @param {ColorFeature | null | undefined} feat @param {string} caption @param {string | null} [color] */
  const updateClickedFeatureLegendState = (feat, caption, color = null) => {
    if (!clickedFeature.value || !feat || clickedFeature.value.svg_id !== feat.svg_id) return;
    if (color) showClickedFeatureColor(feat, color);
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
        replaceLegendRows(drawing.slices.paint, drawing.legendEntries.value);
      }
    }

    legendRenameDialog.show = false;
    legendRenameDialog.mode = 'scope';
    legendRenameDialog.sourceKey = '';
    legendRenameDialog.newCaption = '';
    legendRenameDialog.targetKey = '';
    legendRenameDialog.targetCaption = '';
    legendRenameDialog.targetColor = '';
    legendRenameDialog.currentColor = '';
    legendRenameDialog.siblingCount = 0;
    legendRenameDialog.mergeAvailable = true;
    legendRenameDialog.deletedTargetKey = '';
    legendRenameDialog.reason = '';
    legendRenameDialog.keepBothName = '';
    legendRenameDialog.targetNameAfter = '';
    legendRenameDialog.keepBothRestores = false;
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

    const overrideColor = fillRecordOf(getFeatureOverride(drawing.featureColorOverrides, feat))?.color;
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

  // Whether a rename takes every feature of the row of `sourceKey` (also a
  // row without features): the row is renamed whole, in its place; else the
  // features split off into a row of their own.
  /** @param {PythonLegendKey} sourceKey @param {readonly ColorFeature[]} features */
  const takesWholeRow = (sourceKey, features) => {
    const selected = new Set(features.map(feature => feature.svg_id));
    return getFeaturesForLegendCaption(sourceKey).every(feature => selected.has(feature.svg_id));
  };

  // The rename of a Legend row without rules or features writes the intent
  // (U3a gap 3): the row keeps its key, under which its styles stay (D-25),
  // and takes the new name; Generate and the port rename a row of Python's
  // (`legendRenames`, PV-02) and name an editor row (`legendAdds`), in its
  // place in the order. The text control's open step holds the edit.
  /** @param {DrawingState} drawing @param {PythonLegendKey} sourceKey @param {string} newCaption */
  const renameLegendRow = (drawing, sourceKey, newCaption) => {
    const legendEntry = drawing.legendEntries.value.find((entry) => legendEntryKey(entry) === sourceKey);
    if (!legendEntry) return false;
    // An entry without a key record (an older editor row) records its key
    // first. The rename's own step holds it. PAINT_INTENT_WRITERS does not
    // list this writer: before the paint slice it wrote through the row, and
    // the R2 guard sees neither that nor a handle held under another name.
    // Listing it is an authority change (open: A1 STATUS M3.10, R2 rename).
    const legendPaint = drawing.slices.paint;
    renameListedLegendRow(legendPaint, legendEntry, { caption: /** @type {LegendCaption} */ (newCaption), originalCaption: sourceKey });
    showEditorIntent?.({ domains: RENAMED_ROW_DOMAINS });
    return true;
  };

  // A rename commits once the rules it may add are prepared; its dialogs read
  // only the saved rules (OV-225).
  /** @param {DrawingState} drawing @param {LegendRenameRequest} request */
  const applyLegendRenameRequest = (drawing, request) => withTargetRules(drawing, [], () => applyPreparedLegendRename(drawing, request));
  /** @param {DrawingState} drawing @param {LegendRenameRequest} request */
  const applyPreparedLegendRename = async (drawing, request) => {
    const caption = normalizeCaption(request.finalCaption || request.newCaption);
    const color = resolveColorToHex(request.finalColor || request.currentColor) || '#cccccc';
    if (!caption) return false;
    const features = (request.features || []).filter(Boolean);
    const sourceRules = getLegendRowRules(request.sourceKey);
    if (sourceRules.length || features.length) {
      /** @type {SpecificRuleBuild} */
      const rules = request.sourceScope === 'group' && sourceRules.length
        ? (current, revise) => current.map(rule => (sourceRules.includes(rule) ? revise(rule, { cap: caption, color }) : rule))
        : groupColorRules(features, color, caption, 'group');
      const retireOld = features.length > 0 && takesWholeRow(request.sourceKey, features);
      const oldEntry = drawing.legendEntries.value.find((entry) => legendEntryKey(entry) === request.sourceKey);
      return commitRuleWrite(rules, msg('Rename legend item'), {
        drawing,
        previousLegendIntents: retireOld && oldEntry ? [{ caption: request.sourceKey, color: listedRowColor(drawing, oldEntry) }] : [],
        // OV-158 (Owner decision 2026-10-07): a renamed row keeps its place
        // and every record stored under its key; a row that joins another row
        // keeps that row's place.
        // Features that Keep both split off go right after the target row (D-29).
        legendPlacement: retireOld && oldEntry ? { name: caption, color, at: request.sourceKey, replaces: true }
          : request.keepBothAfter ? { name: caption, color, at: request.keepBothAfter } : null,
        afterCommit: () => {
          const committedCaptionOf = effectiveLegendCaptions();
          for (const feature of features) updateClickedFeatureLegendState(feature, committedCaptionOf(feature), color);
        }
      });
    }
    return renameLegendRow(drawing, request.sourceKey, caption);
  };

  /** @param {LegendRenameRequest} request @param {number} siblingCount */
  const openLegendRenameScopeDialog = (request, siblingCount) => {
    legendRenameDialog.show = true;
    legendRenameDialog.mode = 'scope';
    legendRenameDialog.sourceKey = request.sourceKey;
    legendRenameDialog.newCaption = request.newCaption;
    legendRenameDialog.targetKey = '';
    legendRenameDialog.targetCaption = '';
    legendRenameDialog.targetColor = '';
    legendRenameDialog.currentColor = request.currentColor || '';
    legendRenameDialog.siblingCount = Math.max(0, siblingCount);
    legendRenameDialog.deletedTargetKey = '';
    legendRenameDialog.reason = '';
    legendRenameDialog.pendingRequest = request;
  };

  /**
   * @param {LegendRenameRequest} request
   * @param {RenameTarget} target
   * @param {boolean} mergeAvailable
   * @param {PythonLegendKey | ''} deletedTargetKey The key of a deleted target row, else ''.
   * @param {{ keepBothName: string, targetNameAfter: string | null, keepBothRestores: boolean }} keepBoth
   *   What "Keep both" names the rows, and whether it restores the deleted target first.
   */
  const openLegendRenameTargetDialog = (request, target, mergeAvailable, deletedTargetKey, keepBoth) => {
    legendRenameDialog.show = true;
    legendRenameDialog.mode = 'target';
    legendRenameDialog.sourceKey = request.sourceKey;
    legendRenameDialog.newCaption = request.newCaption;
    legendRenameDialog.targetKey = target.key ?? '';
    legendRenameDialog.targetCaption = target.caption;
    legendRenameDialog.targetColor = target.color;
    legendRenameDialog.currentColor = request.currentColor || '';
    legendRenameDialog.siblingCount = request.siblingCount || 0;
    legendRenameDialog.mergeAvailable = mergeAvailable;
    legendRenameDialog.deletedTargetKey = deletedTargetKey;
    legendRenameDialog.reason = legendNameConflictReason({ name: request.newCaption, ruleOnly: target.ruleOnly,
      deleted: deletedTargetKey !== '', mergeAvailable, sameColor: colorsMatch(target.color, request.currentColor) });
    legendRenameDialog.keepBothName = keepBoth.keepBothName;
    legendRenameDialog.targetNameAfter = keepBoth.targetNameAfter ?? '';
    legendRenameDialog.keepBothRestores = keepBoth.keepBothRestores;
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
      (request.feat ? getCurrentFeatureFillColor(drawing, request.feat) : resolveColorToHex(listedRowColor(drawing, findLegendRow(drawing, request.sourceKey)))) ||
      '#cccccc';

    let features = Array.isArray(request.features) ? request.features.filter(Boolean) : [];
    if (!request.sourceScope) {
      const availableFeatures = request.source === 'popup' ? getFeaturesForLegendCaption(request.sourceKey) : features;
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
      features = features.length > 0 ? features : getFeaturesForLegendCaption(request.sourceKey);
    } else {
      features = [];
    }

    // D-06 (PD-OI-061): a rename onto another name of a different color asks
    // Merge, Keep both, or Cancel, with or without features.
    // Revision 4 (OV-580 to OV-582): the name of a color rule is such a name,
    // whether its rules draw a row or not, and a pick from the popup's list
    // (`joinGroup`) is its own Merge where Merge is allowed.
    // R15-3 (OV-285): a rename onto a deleted row's caption always asks, as onto
    // a listed row; its Merge is a Restore of that row, then the merge. The
    // caption only detects the conflict: the row is its Python key (an editor
    // row's own caption), which the composition root's Restore acts on.
    const group = request.joinGroup || null;
    // A picked group's row is the row the allocation draws for it (N-06), never the row of its name; a
    // typed name's row is the one the Legend owner finds by its drawn name, else by its name (PD-OI-061
    // scenario revision 5: a shown numbered name such as `CDS (2)` names that row).
    const targetKey = group ? groupLegendRowKeys(displayedLegendRowContext(state, drawing))(group) : findLegendRowNamed(newCaption, currentColor);
    /** @template {KeptLegendRow} R @param {readonly R[]} entries @returns {R | null} */
    const rowOf = (entries) => entries.find((entry) => legendEntryKey(entry) === targetKey) || null;
    const listedTarget = rowOf(drawing.legendEntries.value);
    const deletedTarget = listedTarget ? null : rowOf(drawing.deletedLegendEntries.value);
    const targetRow = listedTarget || deletedTarget;
    // A typed name that no row shows is a color-rule group's: of the groups
    // of that name, the one in the renamed row's color (one name and one
    // color share a Legend row), else the first.
    const namedGroups = group || targetRow ? []
      : colorRuleLegendGroups(drawing.manualSpecificRules).filter((candidate) => captionsMatch(candidate.name, newCaption));
    const targetGroup = group || namedGroups.find((candidate) => colorsMatch(candidate.color, currentColor)) || namedGroups[0] || null;
    /** @param {readonly ColorFeature[]} rowFeatures */
    const featureTypes = (rowFeatures) => new Set(rowFeatures.map((feature) => feature?.type));
    // The name and the color a Merge adopts, and the feature types the target
    // draws: a group's own (a picked one also when its row is listed), else
    // the row's name (never its key, OV-640), the color it shows (OV-282: the
    // palette's for a row the palette colors), and the types of the features
    // of its row. The renamed row itself is no target.
    /** @type {RenameTarget | null} */
    const target = targetRow && legendEntryKey(targetRow) === request.sourceKey ? null : targetRow ? {
      key: legendEntryKey(targetRow),
      caption: group ? group.name : legendEntryName(targetRow),
      color: group ? group.color : String((deletedTarget
        ? displayedLegendRowColors(state, drawing, readDrawnTracks()).deleted(deletedTarget) : listedRowColor(drawing, listedTarget)) || ''),
      types: featureTypes(getFeaturesForLegendCaption(legendEntryKey(targetRow))),
      ruleOnly: false
    } : targetGroup && { key: null, caption: targetGroup.name, color: targetGroup.color, types: new Set(targetGroup.types), ruleOnly: true };
    const featureOrRuleRename = features.length > 0 || getLegendRowRules(request.sourceKey).length > 0;
    // OV-62 (PD-OI-061 amended): two rows merge only when each draws features of
    // one type and the type is the same. A row without features (GC content,
    // GC skew), rows of different types, and a row that spans several types
    // offer Keep both and Cancel only. The types are those of the features
    // that take each row's caption, as the live editor knows them: a generated
    // row such as `other proteins` has none live, so it is not merged into. A
    // color-rule group that draws no row has the types of its rules. A picked
    // group takes each type one of its rules is for, whatever its row draws:
    // the popup lists it for that type, and every listed group can be joined.
    const sourceTypes = featureTypes(features);
    const sourceType = [...sourceTypes][0];
    const mergeAllowed = target !== null && sourceTypes.size === 1
      && (group ? colorRuleGroupTakes(group, String(sourceType)) : target.types.size === 1 && target.types.has(sourceType));
    const sameColor = target !== null && colorsMatch(target.color, currentColor);
    const resolution = request.targetResolution || (group && mergeAllowed && !deletedTarget ? 'merge' : null);

    // A chosen Merge adopts the target's color, also once a Restore made a
    // deleted target a listed one. A row that joins a target of its color
    // without asking takes the target's name, never a typed drawn name such
    // as `Shared (2)` (review M-3).
    if (featureOrRuleRename && resolution !== 'merge' && (!target || (!deletedTarget && mergeAllowed && sameColor))) {
      await applyLegendRenameRequest(drawing, { ...request, currentColor, features,
        finalCaption: target ? target.caption : newCaption, finalColor: currentColor });
      clearLegendRenameDialog(drawing);
      return;
    }

    // PD-OI-061 scenario revision 5: Keep both keeps the typed name, and
    // every row of the name is numbered in legend order. A row renamed whole
    // keeps its place; features split off their row go right after the target
    // row (D-29). Rows of one name and one swatch are one row (D-28): onto a
    // deleted target with the renamed row's swatch (the Legend owner's
    // `restores`), Keep both restores the target first (the composition root,
    // as for Merge), so the row shows.
    const splitsRow = !takesWholeRow(request.sourceKey, features);
    const keepBothAfter = splitsRow && listedTarget ? legendEntryKey(listedTarget) : undefined;
    if (target && (deletedTarget || !mergeAllowed || !sameColor)) {
      if (!resolution) {
        openLegendRenameTargetDialog(
          { ...request, currentColor, features },
          target,
          mergeAllowed,
          deletedTarget ? legendEntryKey(deletedTarget) : '',
          (({ restores, ...names }) => ({ ...names, keepBothRestores: restores }))(keepBothNamesOf({
            name: newCaption, renamedKey: splitsRow ? null : request.sourceKey,
            targetKey: listedTarget ? legendEntryKey(listedTarget) : null, color: currentColor,
            deletedKey: deletedTarget ? legendEntryKey(deletedTarget) : null
          }))
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
          : { finalCaption: newCaption, finalColor: currentColor, ...(keepBothAfter ? { keepBothAfter } : {}) })
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
   * @param {LegendRowName} caption The name the features' rules are captioned with (never a row key, OV-640).
   * @param {string} color
   * @param {RuleSynthesisScope} [scope] "group" (default), or "features" for an edit of each feature by itself.
   */
  const applyColorToFeatureGroup = async (drawing, features, caption, color, scope = 'group') => {
    if (!features?.length || !normalizeCaption(caption)) return false;
    const existingEntry = findLegendRowOfName(drawing, caption, color);
    const contributors = getFeaturesForLegendCaption(caption);
    const captionOf = effectiveLegendCaptions();
    const selectedIds = new Set(features.map(feature => feature.svg_id));
    const replacesExistingGroup = existingEntry && contributors.length > 0
      && features.every(feature => captionsMatch(captionOf(feature), caption))
      && contributors.every(feature => selectedIds.has(feature.svg_id));
    return commitRuleWrite(groupColorRules(features, color, normalizeCaption(caption), scope), msg('Change feature color'), {
      drawing,
      previousLegendIntents: replacesExistingGroup
        ? [{ caption: legendEntryKey(existingEntry), color: listedRowColor(drawing, existingEntry) }] : [],
      afterCommit: (/** @type {Array<{ caption: string, name?: string, color: string }>} */ intents) => {
        if (replacesExistingGroup) {
          // An intent's `caption` is its row's key; its `name` the rule's caption.
          const intent = intents.find(entry => captionsMatch(entry.name ?? entry.caption, caption)
            && colorsMatch(entry.color, color));
          if (intent) setLegendColor(drawing.slices.paint, /** @type {PythonLegendKey} */ (intent.caption), intent.color);
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

  // OIC-029: a palette row's color is its type's default color; the rows merged into it (D-35) follow. The one
  // route of both the scope choice and the Legend popup (`setTypeDefaultColor`): the row first releases its own
  // Legend color through the Legend owner, so it follows its type's default color again.
  /** @param {DrawingState} drawing @param {string} type @param {string} color */
  const applyDefaultColorToLegendRow = async (drawing, type, color) => {
    writeLegendRowColor(drawing, [], color, { inherit: [/** @type {PythonLegendKey} */ (type)] });
    const keys = legendRowKeys(getCurrentSvg(), /** @type {PythonLegendKey} */ (type));
    setDefaultColor(drawing, type, color);
    const merged = keys.filter((key) => key !== type);
    if (merged.length && !(await applyColorToLegendSpecificRules(drawing, merged, color))) writeLegendRowColor(drawing, merged, color);
  };

  // OV-1281: the type whose default color draws the row of `key`, else null (reads prepared rule matches).
  /** @param {DrawingState} drawing @param {string} key */
  const legendRowDefaultColorType = (drawing, key) => {
    const features = getFeaturesForLegendCaption(key);
    return features.length > 0 ? paletteRowType(drawing, key, features) : null;
  };

  // The rules of the rows of `keys` (a Legend row and the rows merged into it,
  // `legendRowKeys`) take `color`; a key without rules takes it as its Legend
  // color in the rule commit's step, so a merged row stays one row (D-35).
  /** @param {DrawingState} drawing @param {readonly PythonLegendKey[]} keys @param {string} color */
  const applyColorToLegendSpecificRules = async (drawing, keys, color) => {
    const rows = keys.map((key) => ({ key, rules: getLegendRowRules(key) }));
    const rowRules = rows.flatMap((row) => row.rules);
    const legendOnly = rows.filter((row) => row.rules.length === 0).map((row) => row.key);
    const specificRules = rowRules.filter(rule => !isHashSpecificRule(rule));
    if (!specificRules.length) return false;
    const specificMatches = ruleMatcher(specificRules);
    const coveredExactly = targetsSomeFeatureExactly(
      /** @type {ColorFeature[]} */ (extractedFeatures.value).filter(feature => specificMatches.matchesAny(feature) === true)
    );
    const rowRuleSet = new Set(rowRules);
    return commitRuleWrite((current, revise) => current.filter(rule => !(rowRuleSet.has(rule) && coveredExactly(rule)))
      .map(rule => (rowRuleSet.has(rule) ? revise(rule, { color }) : rule)), msg('Change legend color'),
    { drawing, afterCommit: () => { if (legendOnly.length) writeLegendRowColor(drawing, legendOnly, color); } });
  };

  /**
   * @param {DrawingState} drawing
   * @param {ColorFeature} feat
   * @param {string} color
   * @param {string | null} [requestedLegendName]
   */
  const requestFeatureColorChange = async (drawing, feat, color, requestedLegendName = null) => {
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
        replacedDefaultColor
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
      await setFeatureColor(drawing, feat, color, legendRowName(scope.legendName));
    });
  };

  /** @param {DrawingState} drawing @param {string} color */
  const updateClickedFeatureColor = async (drawing, color) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!clickedFeature.value) return;
    const feat = clickedFeature.value.feat;
    if (!feat) return;
    const customName = popupLegendCaption();
    await requestFeatureColorChange(drawing, feat, color, customName);
  };

  // The caption the popup's Legend name field holds: the key of the feature's
  // row while the field names that row (its key, its name, or its drawn name,
  // which the field shows), else the typed text.
  const popupLegendCaption = () => {
    const typed = normalizeCaption(clickedFeature.value?.legendName);
    const applied = normalizeCaption(clickedFeature.value?.appliedLegendName);
    return applied && legendRowIsNamed(applied, typed) ? applied : typed;
  };

  // The popup's Legend name, typed or picked from the list of color-rule groups.
  /** @param {DrawingState} drawing @param {ColorRuleLegendGroup} [joinGroup] The picked group. */
  const commitLegendName = async (drawing, joinGroup) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!clickedFeature.value) return;

    const feat = clickedFeature.value.feat;
    if (!feat) return;

    const requestedCaption = popupLegendCaption();
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

      if (choice === 'keep-both') {
        await continueLegendRenameRequest(drawing, {
          ...pendingRequest,
          targetResolution: 'keep-both'
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
    const key = legendEntryKey(entry);
    // The row's name or drawn name is no rename; its key names a renamed row back (OV-1216).
    if (!requestedCaption || requestedCaption === normalizeCaption(entry.caption)
      || (requestedCaption !== key && legendRowIsNamed(key, requestedCaption))) {
      replaceLegendRows(drawing.slices.paint, drawing.legendEntries.value);
      return;
    }

    const shownColor = listedRowColor(drawing, entry);
    await continueLegendRenameRequest(drawing, {
      source: 'legend',
      oldCaption: normalizeCaption(entry.caption),
      sourceKey: key,
      newCaption: requestedCaption,
      currentColor: resolveColorToHex(shownColor) || shownColor || '#cccccc',
      features: getFeaturesForLegendCaption(key),
      sourceScope: getFeaturesForLegendCaption(key).length > 0 ? 'group' : 'manual'
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
        ? revise(rule, { color }) : rule)), msg('Change specific color rule'), { drawing, afterCommit: () => showClickedFeatureColor(feat, color) });
    } else if (choice === 'caption') {
      const targetLegendName = normalizeCaption(legendName) || normalizeCaption(getEffectiveLegendCaption(feat));
      if (!targetLegendName) {
        clearFeatureStyleScopeDialog();
        return;
      }
      // The rows the Legend row stands for (D-35), so it stays one row.
      const keys = legendRowKeys(getCurrentSvg(), /** @type {PythonLegendKey} */ (targetLegendName));
      // The choice does what the dialog showed: its default-color line, or rules.
      if (defaultColorType) {
        await applyDefaultColorToLegendRow(drawing, defaultColorType, color);
        showClickedFeatureColor(feat, color);
      } else if (await applyColorToLegendSpecificRules(drawing, keys, color)) {
        showClickedFeatureColor(feat, color);
      } else {
        const siblings = keys.flatMap((key) => findFeaturesWithSameLegendItem(feat, key));
        await applyColorToFeatureGroup(drawing, [feat, ...siblings], legendRowName(targetLegendName), color);
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
      await applyColorToFeatureGroup(drawing, allFeatures, legendRowName(displayLabel), color);
    } else if (choice === 'single') {
      // The dialog holds the feature's row (a key) or a typed name: the name
      // is compared and written, so a feature that leaves its rule's row is
      // named by itself, never by the key of the row it leaves (OV-640).
      let singleCaption = legendRowName(legendName);
      if (featureStyleScopeDialog.siblingCount > 0 || (matchingRule && featureStyleScopeDialog.ruleMatchCount > 1)) {
        const ruleCaption = matchingRule ? (matchingRule.cap || matchingRule.val) : singleCaption;
        if (singleCaption === ruleCaption) {
          singleCaption = legendRowName(getIndividualFeatureLabel(feat));
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
      await applyColorToFeatureGroup(drawing, allFeatures, legendRowName(annotationLabel), color);
    } else if (choice === 'useExisting') {
      if (existingCaptionColor) {
        const targetLegendName = normalizeCaption(legendName) || normalizeCaption(getEffectiveLegendCaption(feat));
        await setFeatureColor(drawing, feat, existingCaptionColor, legendRowName(targetLegendName));
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

    const scope = getFeatureStyleScope(feat, popupLegendCaption());
    if (!scope) return false;
    if (scope.needsDialog) {
      openFeatureStyleScopeDialog({
        kind: 'stroke',
        feat,
        scope,
        strokeColor: normalizedStrokeColor,
        strokeWidth: normalizedStrokeWidth
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
    const { strokeColor: _removed, ...kept } = override;
    setFeatureStroke(drawing.slices.paint, key, hasOwn(kept, 'strokeWidth') ? kept : null);
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
    if (choice === 'this_with_legend') return color ? setFeatureColor(drawing, feature, color, legendRowName(caption)) : false;
    // "Reset all" resets the features of the row: the rules that draw it, by its key.
    const rowRules = new Set(getLegendRowRules(caption));
    // "This feature only" takes back the `hash` rules that name the feature, in
    // Python: a rule another copy of it shares keeps coloring that copy. A
    // feature that another rule still colors takes its palette color by itself.
    /** @type {SpecificRuleBuild} */
    const reset = choice === 'all'
      ? (current) => current.filter(rule => !rowRules.has(rule))
      : groupColorRules([feature], color || '', feature.type === 'CDS' ? 'other proteins' : `other ${feature.type}s`, 'reset');
    const applied = await commitRuleWrite(reset, msg('Reset feature color'), { drawing });
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
    await applyColorToFeatureGroup(drawing, targetFeatures, legendRowName(targetCaption), targetColor);
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
    const targetLegendEntry = findScopeLegendRow(drawing, caption);
    const svg = getCurrentSvg();
    if (!targetLegendEntry || !svg) return false;

    const normalizedStrokeColor = String(strokeColor || '').trim();
    const normalizedStrokeWidth = normalizeStrokeWidthValue(strokeWidth);
    if (!normalizedStrokeColor && normalizedStrokeWidth === null) return false;

    // Every row the Legend row stands for (D-35), so it stays one row.
    return legendRowKeys(svg, legendEntryKey(targetLegendEntry)).map((overrideKey) => {
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
      if (stateChanged) setLegendStroke(drawing.slices.paint, overrideKey, nextOverride);
      return stateChanged;
    }).some(Boolean);
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
      targetFeatures = [feat, ...legendRowKeys(getCurrentSvg(), legendName).flatMap((key) => findFeaturesWithSameLegendItem(feat, key))];
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
   * @param {LegendRowName | null} [customCaption] The name to caption the feature's rule with.
   */
  const setFeatureColor = async (drawing, feature, color, customCaption = null) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!feature || !getFeatureQualifier(feature)) return false;
    const caption = legendRowName(normalizeCaption(customCaption || getIndividualFeatureLabel(feature)));
    if (!caption || featureColorAssignmentMatches(drawing, feature, color, caption)) return false;
    return applyColorToFeatureGroup(drawing, [feature], caption, resolveColorToHex(color) || color, 'features');
  };

  /** @param {DrawingState} drawing @param {ColorFeature} feature @param {string | null} value @param {LegendRowName | null} [customCaption] */
  const setFeatureColorValue = async (drawing, feature, value, customCaption = null) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!feature) return false;
    if (value === null) {
      return commitRuleWrite(groupColorRules([feature], '', '', 'reset'), msg('Reset feature color'), { drawing });
    }
    if (String(value).trim().toLowerCase() === 'none') {
      return applyColorToFeatureGroup(drawing, [feature], legendRowName(normalizeCaption(customCaption || getEffectiveLegendCaption(feature) || feature.type)), 'none', 'features');
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
    // Its rules are Python's synthesis, which Python verified and the commit
    // prepares; a hash rule per selected feature would be matched for nothing.
    applyColorToSelectedFeatures: colorAction(applyColorToSelectedFeatures, { targetRules: false }),
    applyStrokeToSelectedFeatures: strokeAction(applyStrokeToSelectedFeatures),
    resetClickedFeatureFillColor: colorAction(resetClickedFeatureFillColor, { targetRules: false }),
    resetClickedFeatureStroke: strokeAction(resetClickedFeatureStroke),
    getFeatureStrokeColorValue,
    setClickedFeatureStrokeColorValue: strokeAction(setClickedFeatureStrokeColorValue),
    setClickedFeatureStrokeWidthValue: strokeAction(setClickedFeatureStrokeWidthValue),
    setFeatureColor: colorAction(setFeatureColor),
    setFeatureColorValue: colorAction(setFeatureColorValue),
    updateClickedFeatureColor: colorAction(updateClickedFeatureColor, { targetRules: false }),
    updateClickedFeatureStroke: strokeAction(updateClickedFeatureStroke),
    /** @param {string} key */
    legendRowDefaultColorType: (key) => legendRowDefaultColorType(state.activeDrawing(), key),
    preparedLegendRowDefaultColorType: colorAction(legendRowDefaultColorType, { targetRules: false }),
    setTypeDefaultColor: colorAction(applyDefaultColorToLegendRow, { targetRules: false })
  };
};
