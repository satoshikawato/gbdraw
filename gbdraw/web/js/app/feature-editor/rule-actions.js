// @ts-check
/** @import { DrawingState } from '../../state.js' */
import { createSpecificRulePatternDrafts } from './pattern-drafts.js';
import { normalizeUserFacingError } from '../../utils/error-normalization.js';
import { runWhenPrepared } from '../rule-matching.js';
import { ruleMatcher, ruleMatchesFeature } from '../../services/rule-matchers.js';
import { resolveColorToHex } from '../../utils/color-utils.js';
import { parseSpecificRules, serializeSpecificRules } from '../../services/file-imports.js';
import { formatFeatureRange, getFeatureColorRuleHash } from '../../services/feature-utils.js';
import {
  buildLegendIntents, createRuleLegendCaptions, legendRowRules, rendererLegendRows, ruleLegendCaptions
} from '../../services/specific-color-rules.js';
import { resolveFeatureLabelSelector } from '../../services/feature-selector.js';
import { downloadTextFile } from '../../services/text-download.js';
import {
  featureOverrideKey,
  getFeatureOverride,
  migrateLegacyFeatureOverrides
} from '../../services/feature-override-identity.js';
import {
  defaultFeatureRendering,
  normalizeFeatureRendering
} from '../../utils/feature-rendering.js';
import { featureOverrideValue } from '../../services/feature-placement.js';
import { featureDrawnContext, resultLegendSources, sameLegendSources } from '../../services/feature-visibility.js';

// R13: `projectPaletteAndRules` is the composition root's projection of the
// palette and the specific-color rules (R3); this owner calls it after a rule
// commit, whose candidate rules it has prepared. The composition root
// registers the label owner's `requestAutomaticRerender` in `ports` once that
// owner exists; this owner only calls it.
/**
 * The label owner's reaction the composition root registers once that owner exists.
 * @typedef {object} RuleActionsPorts
 * @property {() => boolean} requestAutomaticRerender
 *   Asks for the automatic rerender of an edit Python must draw again (OV-43).
 */

/**
 * A Legend row a rule commit draws in the place of the row `at` (OV-158).
 * @typedef {{ caption: string, at: string }} LegendPlacement
 */

/**
 * The Legend rows a rule commit draws, prepared and not yet applied.
 * @typedef {object} PreparedFileLegend
 * @property {{ add: any[], update: any[], remove: any[], unchanged: any[] }} diff
 * @property {() => boolean} isCurrent
 * @property {() => void} apply
 */

/**
 * @typedef {object} FeatureRuleActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(intents: Record<string, any>[], options?: { drawing?: DrawingState, previousFileIntents?: Record<string, any>[], isCurrent?: () => boolean, placement?: LegendPlacement | null }) => Promise<PreparedFileLegend | false>} prepareFileLegendEntries
 *   The Legend owner's preparation of the rows the rules draw.
 * @property {import('../rule-matching.js').RulePreparation} rulePreparation
 * @property {(label: string, fn: () => any, options?: Record<string, any>) => any} runUndoable History's undoable step.
 * @property {(label: string, fn: () => any, options?: Record<string, any>) => any} runUndoableCheckpoint
 *   History's undoable step that stores a checkpoint of the Result.
 * @property {(options?: { recolor?: Record<string, any>, prepareRules?: boolean }) => boolean | Promise<boolean>} projectPaletteAndRules
 *   The root's projection of the palette and the rules (R3).
 * @property {RuleActionsPorts} ports
 * @property {() => ({ diagramOptions?: Record<string, any> } | null)} [getCommittedRequest]
 *   The committed canonical request (Python owns the option fields, R7).
 * @property {(value?: any) => { value: any }} ref Vue `ref`
 * @property {<T>(getter: () => T) => { value: T }} computed Vue `computed`
 * @property {(source: { value: any }, callback: (value: any) => void) => () => void} watch
 *   Vue `watch`; a rule commit waits on it for an automatic rerender in flight.
 * @property {() => boolean} [isPatternEditAvailable] False while a Session import is pending.
 */

/** @param {FeatureRuleActionsOptions} options */
export const createFeatureRuleActions = ({ state, prepareFileLegendEntries, rulePreparation, runUndoable, runUndoableCheckpoint, projectPaletteAndRules, ports, getCommittedRequest = () => null, ref, computed, watch, isPatternEditAvailable = () => true }) => {
  const {
    appliedPaletteColors,
    newColorFeat,
    newColorVal,
    newSpecRule,
    specificRulePresets,
    selectedSpecificPreset,
    specificRulePresetLoading,
    newPriorityRule,
    newFeatureToAdd,
    extractedFeatures,
    editableLabels
  } = state;

  const normalizeCaption = (value) => String(value || '').trim();
  const normalizeCaptionKey = (value) => normalizeCaption(value).toLowerCase();
  const normalizeFeatureIdKey = (value) => String(value || '').trim().toLowerCase();
  const captionMatches = (value, target) => normalizeCaptionKey(value) === normalizeCaptionKey(target);
  const specificRuleFields = new Set(['feat', 'qual', 'val', 'color', 'cap']);

  let preparationRevision = 0;
  const patternDrafts = createSpecificRulePatternDrafts({
    rules: () => state.activeDrawing().manualSpecificRules, ref, invalidate: () => { preparationRevision += 1; }
  });
  const ruleFailure = ref(null);
  const canRetrySpecificRuleFailure = computed(() => Boolean(ruleFailure.value
    && state.errorLog?.value === ruleFailure.value.error
    && rulePreparation.isCurrent(ruleFailure.value.snapshot)
    && (!ruleFailure.value.presetId || selectedSpecificPreset.value === ruleFailure.value.presetId)));
  const canEditSpecificRuleFailure = computed(() => canRetrySpecificRuleFailure.value
    && Boolean(ruleFailure.value.input?.isConnected));
  const retrySpecificRuleFailure = () => {
    const busy = state.sessionOperationAvailability?.();
    return busy || (canRetrySpecificRuleFailure.value ? ruleFailure.value.retry() : false);
  };
  const editSpecificRuleFailure = () => {
    if (canEditSpecificRuleFailure.value) ruleFailure.value.input.focus();
  };
  // The legend rows the candidate rules draw on rendered features, and the rows
  // the commit retires: those the current rules draw (so it retires the row
  // Generate drew) and `retiredLegendIntents`, rows this commit replaces, which
  // are no renderer rows for the N-06 caption allocation. `removedRuleRows` are
  // the rows of current rules that no candidate rule names, drawn or not.
  /** @param {DrawingState} drawing */
  const candidateLegendIntents = (drawing, candidateRules, retiredLegendIntents) => {
    const rendered = (extractedFeatures.value || []).filter(feature =>
      featureOverrideValue(drawing.featureOverrides, feature, 'featureVisibility') !== 'off');
    const candidateMatches = ruleMatcher(candidateRules);
    const used = new Set(rendered.map(feature => candidateMatches.first(feature)).filter(Boolean));
    const rendererRows = rendererLegendRows({
      legendEntries: drawing.legendEntries?.value,
      originalLegendOrder: state.originalLegendOrder?.value,
      rules: [...drawing.manualSpecificRules, ...candidateRules,
        ...retiredLegendIntents.map(intent => ({ cap: intent?.caption, color: intent?.color }))]
    });
    const currentCaption = createRuleLegendCaptions(drawing.manualSpecificRules, rendererRows);
    const ruleRows = drawing.manualSpecificRules.filter(rule => rule.cap)
      .map(rule => ({ caption: currentCaption(rule), color: rule.color }));
    const candidateCaptions = new Set(buildLegendIntents(candidateRules, rendererRows).intents.map(intent => intent.caption));
    return {
      intents: buildLegendIntents(candidateRules.filter(rule => used.has(rule)), rendererRows).intents,
      previousIntents: [...ruleRows, ...retiredLegendIntents],
      removedRuleRows: ruleRows.filter(row => !candidateCaptions.has(row.caption))
    };
  };
  // OV-43 (Owner decision 2026-10-06, option A): whether a rule table change
  // changes a Result's Legend source, read with the current Feature
  // visibility before and after. Python redraws the Legend rows, their order,
  // and the "other <type>s" rows in the automatic rerender.
  /** @param {DrawingState} drawing */
  const legendSourceContext = (drawing, rules) => ({
    ...featureDrawnContext(drawing, { diagramOptions: getCommittedRequest()?.diagramOptions }),
    colorRules: rules
  });
  /** @param {DrawingState} drawing */
  const changesLegendSource = (drawing, before, after) => !sameLegendSources(
    resultLegendSources(state, legendSourceContext(drawing, before)),
    resultLegendSources(state, legendSourceContext(drawing, after))
  );
  // OV-152: rules that recolor a whole type keep its caption, so their Legend
  // source is that of the type's default row and the live Legend follows them
  // (`legendSourceOf`). A change that removes the last such rule of a type
  // makes Python draw the default row again, which the live removal of the
  // rule row does not; the automatic rerender draws it.
  /** @param {Record<string, any>[]} rules */
  const typeRowCaptions = (rules) => new Set(rules
    .filter(rule => String(rule?.cap ?? '').trim() && rule.cap === rule.feat)
    .map(rule => rule.cap));
  /**
   * @param {DrawingState} drawing
   * @param {Record<string, any>[]} before
   * @param {Record<string, any>[]} after
   * @param {(caption: string) => boolean} [shown] Whether the Legend already shows
   *   the default row of a type (a History restore returns the rows it held).
   */
  const redrawsLegendFor = (drawing, before, after, shown = () => false) => {
    const kept = typeRowCaptions(after);
    return changesLegendSource(drawing, before, after)
      || [...typeRowCaptions(before)].some(caption => !kept.has(caption) && !shown(caption));
  };
  // The Legend rows a History restore returned, shown or deleted in the editor.
  /** @param {DrawingState} drawing @param {string} caption */
  const legendHoldsRow = (drawing, caption) => [
    ...(drawing.legendEntries.value || []), ...(drawing.deletedLegendEntries.value || [])
  ].some((entry) => captionMatches(entry?.caption, caption) || captionMatches(entry?.originalCaption, caption));
  // Undo and Redo of a rule edit restore the rules; the composition root
  // passes the rules they replaced, and a changed Legend source asks for the
  // rerender, as the edit did. A removed whole-type rule asks for it only
  // when the restored Legend lacks the type's default row (OV-169): Undo of
  // an Apply to all returns the row it recolored, so it needs no Python run.
  const followRestoredRules = (previousRules) => {
    const drawing = state.activeDrawing();
    return (
      JSON.stringify(previousRules) !== JSON.stringify(drawing.manualSpecificRules)
      && redrawsLegendFor(drawing, previousRules, drawing.manualSpecificRules,
        (caption) => legendHoldsRow(drawing, caption))
      && ports.requestAutomaticRerender()
    );
  };
  // The automatic rerender replaces the Results the candidate was prepared
  // against, which makes the candidate stale (#857). An edit made while one runs
  // waits for it and, when it replaced the Results under an unchanged rule
  // table, prepares again instead of dropping the edit. The rerender is idle
  // before the binder has read the Result it wrote, so a preparation can also
  // go stale on that binding with the Results unchanged (OV-166): any change
  // of the inputs the candidate was prepared from counts, unless another
  // commit, the caller, or a Session operation superseded this one.
  /** @returns {Promise<void>} */
  const rerenderIdle = () => (state.labelReflowProcessing?.value
    ? new Promise((resolve) => {
      const stop = watch(state.labelReflowProcessing, (busy) => { if (!busy) { stop(); resolve(); } });
    })
    : Promise.resolve());
  // OV-61: the rule table is the latest explicit color of the rows its rules
  // draw. A Legend-only color an earlier popup edit left on such a row (the
  // popup stores the caption's color when it recolors a whole group) would
  // win at Generate over the recolored rule, so a commit that recolors the row
  // retires it; `afterCommit` may store a new one.
  // OV-152: a commit that removes a rule row (Clear All, a deleted or
  // recaptioned rule) retires that copy too; left behind, it would paint the
  // row Generate draws for the features' palette color. A copy has the removed
  // rule's caption and color. A color set in the Legend editor is no copy: the
  // editor recolors the rules of a row a rule draws (legend.js), so its colors
  // sit on rows no rule draws and stay.
  /**
   * @param {DrawingState} drawing
   * @param {{ caption: string, color: string }[]} intents
   * @param {{ caption: string, color: string }[]} removedRuleRows
   */
  const retireSupersededLegendColors = (drawing, intents, removedRuleRows) => {
    const overrides = drawing.legendColorOverrides;
    for (const { caption, color } of intents) {
      if (Object.hasOwn(overrides, caption) && resolveColorToHex(overrides[caption]) !== resolveColorToHex(color)) {
        delete overrides[caption];
      }
    }
    for (const { caption, color } of removedRuleRows) {
      if (Object.hasOwn(overrides, caption) && resolveColorToHex(overrides[caption]) === resolveColorToHex(color)) {
        delete overrides[caption];
      }
    }
  };
  /**
   * `drawing` is the drawing of the action (default: the shown mode's); an
   * action that awaits before it commits passes the drawing it started on.
   * @typedef {{ drawing?: DrawingState, isCurrent?: () => boolean, afterCommit?: (intents: Record<string, any>[]) => void, previousLegendIntents?: Record<string, any>[], sourceRows?: (Record<string, any> | null)[], legendPlacement?: LegendPlacement | null }} CommitSpecificRulesOptions
   */
  /**
   * @param {Record<string, any>[]} rules
   * @param {string} [label]
   * @param {CommitSpecificRulesOptions} [options]
   */
  const commitSpecificRules = async (rules, label = 'Change specific color rules', options = {}) => {
    const { drawing = state.activeDrawing(), ...commitOptions } = options;
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    await rerenderIdle();
    for (let attempt = 0; ; attempt += 1) {
      const results = state.results.value;
      const table = JSON.stringify(drawing.manualSpecificRules);
      const inputs = rulePreparation.snapshot();
      const revision = preparationRevision;
      const outcome = await commitOnce(drawing, rules, label, commitOptions);
      const stale = state.results.value !== results || !rulePreparation.isCurrent(inputs);
      if (outcome || attempt === 2 || !stale
        || JSON.stringify(drawing.manualSpecificRules) !== table
        || preparationRevision !== revision + 1
        || commitOptions.isCurrent?.() === false
        || state.sessionOperationAvailability?.()) return outcome;
      await rerenderIdle();
    }
  };
  // R13: this owner builds the candidate rules, so it is the one place that asks
  // the rule preparation for their matches: a commit prepares the candidate
  // with its captions, a color action prepares the rules it reads as they are.
  const prepareRules = (rules, options) => rulePreparation.prepareCandidate(rules, options);
  // Runs `commit` once the color rule matches of `rules` are prepared; a
  // prepared table answers at once, so the commit runs in the caller's tick.
  /**
   * @param {Record<string, any>[]} rules
   * @param {() => any} commit
   */
  const runWithRuleMatches = (rules, commit) => runWhenPrepared(
    state, () => [rulePreparation.isPrepared(rules) || prepareRules(rules, { captions: false })], commit
  );
  /**
   * @param {DrawingState} drawing
   * @param {Record<string, any>[]} rules
   * @param {string} label
   * @param {CommitSpecificRulesOptions} [options]
   */
  const commitOnce = async (drawing, rules, label, { isCurrent = () => true, afterCommit = () => {}, previousLegendIntents = [], sourceRows = rules.map(rule => drawing.manualSpecificRules.includes(rule) ? rule : null), legendPlacement = null } = {}) => {
    const revision = ++preparationRevision;
    const candidate = await prepareRules(rules);
    if (!candidate) return false;
    const current = () => revision === preparationRevision && !state.sessionOperationAvailability?.() && isCurrent()
      && rulePreparation.isCurrent(candidate.snapshot);
    if (!current()) return false;
    const { intents, previousIntents, removedRuleRows } = candidateLegendIntents(drawing, candidate.rules, previousLegendIntents);
    const previousCaptions = new Set(previousIntents.map(intent => intent.caption));
    const legend = await prepareFileLegendEntries(intents.filter(intent => !(drawing.deletedLegendEntries?.value || [])
      .some(entry => (entry.originalCaption || entry.caption) === intent.caption)), {
      drawing,
      previousFileIntents: previousIntents,
      isCurrent: current,
      placement: legendPlacement
    });
    if (!legend) return false;
    const redrawsLegend = redrawsLegendFor(drawing, [...drawing.manualSpecificRules], candidate.rules);
    let applied = false;
    // One History step: the rule transition first, then the legend rows it
    // draws (R13); a checkpoint when the legend gains or loses a row.
    const transact = legend.diff.add.length || legend.diff.remove.length
      ? runUndoableCheckpoint : runUndoable;
    await transact(label, () => {
      if (!current() || !legend.isCurrent()) return false;
      drawing.manualSpecificRules.splice(0, drawing.manualSpecificRules.length, ...candidate.rules.map((rule, index) => {
        const row = sourceRows[index];
        if (!row) return rule;
        Object.assign(row, rule);
        if (!Object.hasOwn(rule, 'fromFile')) delete row.fromFile;
        return row;
      }));
      patternDrafts.reconcile();
      drawing.fileLegendCaptions.value = new Set(candidate.rules.filter(rule => rule.fromFile && rule.cap).map(rule => rule.cap));
      drawing.addedLegendCaptions.value = new Set([
        ...[...drawing.addedLegendCaptions.value].filter(caption => !previousCaptions.has(caption)),
        ...intents.map(intent => intent.caption)
      ]);
      applyRulePreview();
      retireSupersededLegendColors(drawing, intents, removedRuleRows);
      afterCommit(intents);
      applied = true;
      legend.apply();
      return legend.diff;
    });
    if (applied) rulePreparation.notifyChanges(candidate);
    if (applied && redrawsLegend) ports.requestAutomaticRerender();
    return applied;
  };
  /**
   * @param {DrawingState} drawing
   * @param {Element | null} [input] The input a failed commit offers for editing.
   */
  const commitPrepared = async (drawing, rules, label, afterCommit = () => {}, input = null, sourceRows = rules.map(rule => drawing.manualSpecificRules.includes(rule) ? rule : null)) => {
    patternDrafts.suspend();
    const snapshot = rulePreparation.snapshot();
    const previousError = state.errorLog?.value;
    try {
      const applied = await commitSpecificRules(rules, label, { drawing, afterCommit, sourceRows });
      if (applied && state.errorLog?.value === ruleFailure.value?.error) state.errorLog.value = null;
      if (applied) ruleFailure.value = null;
      return applied;
    } catch (cause) {
      if (!rulePreparation.isCurrent(snapshot) || state.errorLog?.value !== previousError) return false;
      const error = normalizeUserFacingError(cause, { operation: 'evaluateRules', stage: 'resource-staging' });
      if (state.errorLog) state.errorLog.value = error;
      ruleFailure.value = { error, snapshot, input, retry: () => commitPrepared(drawing, rules, label, afterCommit, input, sourceRows) };
      return false;
    }
  };
  const applyRulePreview = () => {
    refreshFeatureOverrides(extractedFeatures.value);
    projectPaletteAndRules({ prepareRules: false });
  };

  const editSpecificRulePattern = (row, value) => {
    const busy = state.sessionOperationAvailability?.();
    return busy || patternDrafts.edit(row, value);
  };
  /** @param {HTMLInputElement | null} [input] */
  const revertSpecificRulePattern = (row, input = null) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    patternDrafts.revert(row);
    if (input?.isConnected) {
      input.value = row.val;
      input.focus();
    }
  };
  /** @param {DrawingState} drawing */
  const applySpecificRulePattern = async (drawing, row, value) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    if (!drawing.manualSpecificRules.includes(row) || !isPatternEditAvailable()
      || (state.generatedMode?.value && state.mode?.value !== state.generatedMode.value)) return false;
    const token = patternDrafts.begin(row, value);
    const snapshot = rulePreparation.snapshot();
    const attemptRevision = preparationRevision + 1;
    const current = () => preparationRevision === attemptRevision && patternDrafts.isCurrent(row, token) && rulePreparation.isCurrent(snapshot)
      && isPatternEditAvailable();
    const sourceRows = [...drawing.manualSpecificRules];
    const rules = sourceRows.map(rule => {
      if (rule !== row) return { ...rule };
      const next = { ...rule, val: String(value ?? '') };
      delete next.fromFile;
      return next;
    });
    try {
      const applied = await commitSpecificRules(rules, 'Edit specific color rule', { drawing, isCurrent: current, sourceRows });
      if (applied && patternDrafts.isCurrent(row, token)) patternDrafts.revert(row);
      return applied;
    } catch (cause) {
      if (!current() || cause?.canceled || cause?.name === 'AbortError' || ['canceled', 'stale', 'superseded'].includes(cause?.status)) return false;
      patternDrafts.get(row).error = normalizeUserFacingError(cause, { operation: 'evaluateRules', stage: 'resource-staging' });
      return false;
    } finally {
      if (patternDrafts.isCurrent(row, token)) patternDrafts.get(row).pending = false;
    }
  };
  const retrySpecificRulePattern = (row) => {
    const drawing = state.activeDrawing();
    const draft = patternDrafts.get(row);
    return draft && !draft.pending ? applySpecificRulePattern(drawing, row, draft.text) : false;
  };
  /** @param {HTMLInputElement | null} [input] */
  const setSpecificRuleField = (index, field, value, input = null) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!specificRuleFields.has(field)) return;
    const current = drawing.manualSpecificRules[index];
    if (!current) return;
    if (field === 'val') return applySpecificRulePattern(drawing, current, value);
    const nextValue = field === 'color' ? resolveColorToHex(String(value || '#000000')) : String(value ?? '');
    const nextRule = { ...current, [field]: nextValue };
    delete nextRule.fromFile;
    const sourceRows = [...drawing.manualSpecificRules];
    const rules = sourceRows.map(rule => rule === current ? nextRule : { ...rule });
    return commitPrepared(drawing, rules, 'Edit specific color rule', () => {}, input, sourceRows)?.finally(() => {
      if (input?.isConnected && drawing.manualSpecificRules.includes(current)) input.value = current[field] ?? '';
    });
  };

  /** @param {DrawingState} drawing */
  const moveSpecificRule = (drawing, index, offset) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const sourceRows = [...drawing.manualSpecificRules];
    const target = index + offset;
    if (target < 0 || target >= sourceRows.length) return;
    const [row] = sourceRows.splice(index, 1);
    sourceRows.splice(target, 0, row);
    return commitPrepared(drawing, sourceRows.map(rule => ({ ...rule })), 'Move specific color rule', () => {}, null, sourceRows);
  };

  const downloadSpecificRulesTsv = () => {
    const drawing = state.activeDrawing();
    const text = serializeSpecificRules(drawing.manualSpecificRules);
    if (!text.trim()) {
      alert('No specific rules to export.');
      return;
    }
    downloadTextFile('gbdraw_specific_rules.tsv', text);
  };

  const getIndividualFeatureLabel = (feat) => {
    return feat.product || feat.gene || feat.locus_tag || `${feat.type} at ${formatFeatureRange(feat)}`;
  };

  const getEditableLabelEntryForFeature = (feat) => {
    if (!feat || !Array.isArray(editableLabels.value)) return null;
    const featureIdKey = normalizeFeatureIdKey(feat.svg_id || feat.id);
    if (!featureIdKey) return null;
    return (
      editableLabels.value.find((entry) => normalizeFeatureIdKey(entry?.featureId) === featureIdKey) || null
    );
  };

  const getDisplayedFeatureLabel = (feat) => {
    const drawing = state.activeDrawing();
    if (!feat) return '';

    const editableEntry = getEditableLabelEntryForFeature(feat);
    const editableText = normalizeCaption(editableEntry?.text);
    if (editableText) return editableText;

    const normalizedOverride = normalizeCaption(featureOverrideValue(drawing.featureOverrides, feat, 'labelText'));
    if (normalizedOverride) return normalizedOverride;

    const sourceText = normalizeCaption(editableEntry?.sourceText);
    if (sourceText) {
      const normalizedBulk = normalizeCaption(drawing.labelTextBulkOverrides[sourceText]);
      if (normalizedBulk) return normalizedBulk;
    }

    return normalizeCaption(getIndividualFeatureLabel(feat));
  };

  /** @param {DrawingState} drawing */
  const legendRowContext = (drawing) => ({
    rules: drawing.manualSpecificRules,
    legendEntries: drawing.legendEntries?.value || [],
    originalLegendOrder: state.originalLegendOrder?.value || []
  });
  // The rules a legend row draws; editing the row edits them (N-06).
  const getLegendRowRules = (caption) => {
    const drawing = state.activeDrawing();
    return legendRowRules(caption, legendRowContext(drawing));
  };

  // Resolve the effective legend item label used by current SVG coloring
  // priority: a rule's feature belongs to the row Generate draws for it (N-06).
  // One reader serves a pass over many features (OV-193).
  /** @returns {(feat: Record<string, any> | null | undefined) => string} */
  const effectiveLegendCaptions = () => {
    const drawing = state.activeDrawing();
    const ruleMatches = ruleMatcher(drawing.manualSpecificRules);
    const legendCaption = ruleLegendCaptions(legendRowContext(drawing));
    return (feat) => {
      if (!feat) return '';
      const rule = ruleMatches.first(feat);
      if (rule && normalizeCaption(rule.cap)) return legendCaption(rule);

      const overrideCaption = normalizeCaption(
        getFeatureOverride(drawing.featureColorOverrides, feat)?.caption
      );
      if (overrideCaption) return overrideCaption;

      return normalizeCaption(feat.type) || normalizeCaption(getIndividualFeatureLabel(feat));
    };
  };
  const getEffectiveLegendCaption = (feat) => effectiveLegendCaptions()(feat);

  const addCustomColor = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!newColorFeat.value) return;
    drawing.currentColors.value = {
      ...drawing.currentColors.value,
      [newColorFeat.value]: newColorVal.value
    };
  };

  const setLabelFilterMode = (value) => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    drawing.filterMode.value = value;
  };
  const addWhitelistRule = () => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    drawing.manualWhitelist.push({ feat: 'CDS', qual: 'product', key: '' });
  };
  const removeWhitelistRule = (index) => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    drawing.manualWhitelist.splice(index, 1);
  };
  const removePriorityRule = (index) => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    drawing.manualPriorityRules.splice(index, 1);
  };

  const addPriorityRule = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!newPriorityRule.order) return;
    const idx = drawing.manualPriorityRules.findIndex((r) => r.feat === newPriorityRule.feat);
    if (idx >= 0) {
      drawing.manualPriorityRules[idx].order = newPriorityRule.order;
    } else {
      drawing.manualPriorityRules.push({ feat: newPriorityRule.feat, order: newPriorityRule.order });
    }
  };

  const addSpecificRule = async () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!newSpecRule.val) return;

    if (newSpecRule.val.length > 50) {
      if (!confirm('Regular expression is quite long (>50 chars). This might impact performance. Continue?')) {
        return;
      }
    }

    if (/\(.+[\+\*]\)[\+\*]/.test(newSpecRule.val) || /\(.*\)\+/.test(newSpecRule.val)) {
      if (
        !confirm(
          'This regular expression contains patterns that may freeze the browser (ReDoS risk). Are you sure you want to add it?'
        )
      ) {
        return;
      }
    }

    const rule = {
      feat: String(newSpecRule.feat || ''),
      qual: String(newSpecRule.qual || ''),
      val: String(newSpecRule.val),
      color: String(newSpecRule.color || '#000000'),
      cap: String(newSpecRule.cap || '')
    };
    return commitPrepared(drawing, [...drawing.manualSpecificRules, rule], 'Add specific color rule', () => {
      if (newSpecRule.val === rule.val) newSpecRule.val = '';
    }, document.querySelector?.('[data-new-specific-rule-pattern]'));
  };

  const applySpecificRulePreset = async () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (specificRulePresetLoading.value) return;
    const presetId = selectedSpecificPreset.value;
    if (!presetId) return;
    const preset = specificRulePresets.find((entry) => entry.id === presetId);
    if (!preset) {
      alert('Unknown preset selected.');
      return;
    }

    const presetContext = rulePreparation.snapshot();
    const previousError = state.errorLog?.value;
    specificRulePresetLoading.value = true;
    try {
      const response = await fetch(preset.path, { cache: 'no-store' });
      if (!response.ok) {
        throw { code: 'INPUT_UNREADABLE', operation: 'evaluateRules', stage: 'resource-staging' };
      }
      const text = await response.text();
      const { rules } = parseSpecificRules(text);
      if (selectedSpecificPreset.value !== presetId || !rulePreparation.isCurrent(presetContext)) return;
      return await commitSpecificRules(rules, 'Apply specific color preset', {
        drawing,
        isCurrent: () => selectedSpecificPreset.value === presetId,
        afterCommit: () => {
          if (presetId === 'bakta') {
            drawing.currentColors.value = { ...drawing.currentColors.value, CDS: '#cccccc' };
            drawing.adv.legend_box_size = 12;
            drawing.adv.legend_font_size = 12;
          }
        }
      });
    } catch (e) {
      if (selectedSpecificPreset.value !== presetId || !rulePreparation.isCurrent(presetContext)
        || state.errorLog?.value !== previousError) return;
      const error = normalizeUserFacingError(e, { operation: 'evaluateRules', stage: 'resource-staging' });
      if (state.errorLog) state.errorLog.value = error;
      ruleFailure.value = { error, snapshot: presetContext, presetId, retry: applySpecificRulePreset };
    } finally {
      specificRulePresetLoading.value = false;
    }
  };

  const addFeature = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (newFeatureToAdd.value && !drawing.adv.features.includes(newFeatureToAdd.value)) {
      drawing.adv.features.push(newFeatureToAdd.value);
      if (!drawing.adv.feature_shapes || typeof drawing.adv.feature_shapes !== 'object') {
        drawing.adv.feature_shapes = {};
      }
      if (!Object.prototype.hasOwnProperty.call(drawing.adv.feature_shapes, newFeatureToAdd.value)) {
        drawing.adv.feature_shapes[newFeatureToAdd.value] = defaultFeatureRendering(newFeatureToAdd.value);
      }
    }
  };

  const removeFeature = (featureType) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const idx = drawing.adv.features.indexOf(featureType);
    if (idx >= 0) {
      drawing.adv.features.splice(idx, 1);
    }
  };

  const getFeatureShape = (featureType) => {
    const drawing = state.activeDrawing();
    if (!drawing.adv.feature_shapes || typeof drawing.adv.feature_shapes !== 'object') {
      return defaultFeatureRendering(featureType);
    }
    return Object.prototype.hasOwnProperty.call(drawing.adv.feature_shapes, featureType)
      ? normalizeFeatureRendering(drawing.adv.feature_shapes[featureType])
      : defaultFeatureRendering(featureType);
  };

  const setFeatureShape = (featureType, shape) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!drawing.adv.feature_shapes || typeof drawing.adv.feature_shapes !== 'object') {
      drawing.adv.feature_shapes = {};
    }
    drawing.adv.feature_shapes[featureType] = normalizeFeatureRendering(shape);
  };

  const getFeatureColor = (feat) => {
    const drawing = state.activeDrawing();
    const override = getFeatureOverride(drawing.featureColorOverrides, feat);
    if (override) {
      return resolveColorToHex(override.color || override);
    }
    return resolveColorToHex(appliedPaletteColors.value[feat.type]) || '#cccccc';
  };

  const getFeatureColorValue = (feat) => {
    const drawing = state.activeDrawing();
    const override = getFeatureOverride(drawing.featureColorOverrides, feat);
    if (!override) return null;
    return Object.prototype.hasOwnProperty.call(override, 'color')
      ? override.color
      : override;
  };

  const canEditFeatureColor = () => true;

  // Python matches a single-feature color rule only by the stable generation
  // hash, so duplicate record instances that share one hash share the rule.
  const getFeatureQualifier = (feat) => {
    const colorRuleHash = getFeatureColorRuleHash(feat);
    return colorRuleHash ? { qual: 'hash', val: colorRuleHash } : null;
  };

  const getLabelSpecificRule = (feat, label) => {
    const drawing = state.activeDrawing();
    if (!feat) return null;
    const priorityRule = drawing.manualPriorityRules.find((rule) => rule.feat === feat.type);
    const priority = String(priorityRule?.order || '')
      .split(',')
      .map((qualifier) => qualifier.trim())
      .filter(Boolean);
    const selector = resolveFeatureLabelSelector(feat, label, { priority });
    if (!selector) return null;
    return {
      feat: feat.type,
      qual: selector.qualifier,
      val: selector.pattern
    };
  };

  const refreshFeatureOverrides = (features) => {
    const drawing = state.activeDrawing();
    const ruleMatches = ruleMatcher(drawing.manualSpecificRules);
    if (!features || features.length === 0 || !ruleMatches.ready(features)) return;
    migrateLegacyFeatureOverrides(drawing.featureColorOverrides, features);

    for (const feat of features) {
      const rule = ruleMatches.first(feat);
      const key = featureOverrideKey(feat);
      if (key) {
        if (rule) drawing.featureColorOverrides[key] = { color: rule.color, caption: rule.cap };
        else delete drawing.featureColorOverrides[key];
      }
    }
  };

  const findMatchingRegexRule = (feat) => {
    const drawing = state.activeDrawing();
    return ruleMatcher(drawing.manualSpecificRules.filter((rule) => rule.qual !== 'hash')).first(feat);
  };

  const countFeaturesMatchingRule = (rule) => {
    if (!rule || String(rule.qual || '').toLowerCase() === 'hash') return 0;

    let count = 0;
    for (const feat of extractedFeatures.value) {
      if (rule.feat !== '*' && feat.type !== rule.feat) continue;

      if (ruleMatchesFeature(feat, rule)) count++;
    }
    return count;
  };

  /** @param {string | null} [caption] */
  const findFeaturesWithSameLegendItem = (currentFeat, caption = null) => {
    const captionOf = effectiveLegendCaptions();
    const targetCaption = normalizeCaption(caption || captionOf(currentFeat));
    if (!targetCaption) return [];
    return extractedFeatures.value.filter((f) => {
      if (f.svg_id === currentFeat.svg_id) return false;
      return captionMatches(captionOf(f), targetCaption);
    });
  };

  /** @param {string | null} [label] */
  const findFeaturesWithSameIndividualLabel = (currentFeat, label = null) => {
    const targetLabel = normalizeCaption(label || getIndividualFeatureLabel(currentFeat));
    if (!targetLabel) return [];

    return extractedFeatures.value.filter((f) => {
      if (f.svg_id === currentFeat.svg_id) return false;
      return captionMatches(getIndividualFeatureLabel(f), targetLabel);
    });
  };

  /** @param {string | null} [label] */
  const findFeaturesWithSameDisplayedLabel = (currentFeat, label = null) => {
    const targetLabel = normalizeCaption(label || getDisplayedFeatureLabel(currentFeat));
    if (!targetLabel) return [];

    return extractedFeatures.value.filter((f) => {
      if (f.svg_id === currentFeat.svg_id) return false;
      return captionMatches(getDisplayedFeatureLabel(f), targetLabel);
    });
  };

  const findExistingColorForCaption = (_currentFeat, caption) => {
    const drawing = state.activeDrawing();
    const targetCaption = normalizeCaption(caption);
    if (!targetCaption) return null;

    for (const rule of drawing.manualSpecificRules) {
      if (captionMatches(rule.cap, targetCaption) && String(rule.qual || '').toLowerCase() === 'hash') {
        return { rule, color: rule.color };
      }
    }

    for (const rule of drawing.manualSpecificRules) {
      if (captionMatches(rule.cap, targetCaption) && String(rule.qual || '').toLowerCase() !== 'hash') {
        return { rule, color: rule.color };
      }
    }

    return null;
  };

  return {
    specificRulePattern: patternDrafts.text,
    specificRulePatternDraft: patternDrafts.get,
    specificRulePatternFieldId: patternDrafts.fieldId,
    editSpecificRulePattern, retrySpecificRulePattern, revertSpecificRulePattern,
    suspendSpecificRulePatternDrafts: patternDrafts.suspend,
    clearSpecificRulePatternDrafts: patternDrafts.clear,
    captureSpecificRulePatternDrafts: patternDrafts.capture,
    restoreSpecificRulePatternDrafts: patternDrafts.restore,
    canRetrySpecificRuleFailure, canEditSpecificRuleFailure, retrySpecificRuleFailure, editSpecificRuleFailure,
    addCustomColor,
    addFeature,
    removeFeature,
    getFeatureShape,
    setFeatureShape,
    addPriorityRule,
    setLabelFilterMode,
    addWhitelistRule,
    removeWhitelistRule,
    removePriorityRule,
    addSpecificRule,
    applySpecificRulePreset,
    canEditFeatureColor,
    clearAllSpecificRules: () => {
      const drawing = state.activeDrawing();
      return commitPrepared(drawing, [], 'Clear specific color rules');
    },
    commitSpecificRules,
    followRestoredRules,
    countFeaturesMatchingRule,
    downloadSpecificRulesTsv,
    findExistingColorForCaption,
    findFeaturesWithSameCaption: findFeaturesWithSameLegendItem,
    findFeaturesWithSameDisplayedLabel,
    findFeaturesWithSameIndividualLabel,
    findFeaturesWithSameLegendItem,
    findMatchingRegexRule,
    getFeatureColor,
    getFeatureColorValue,
    getDisplayedFeatureLabel,
    effectiveLegendCaptions,
    getEffectiveLegendCaption,
    getIndividualFeatureLabel,
    getFeatureQualifier,
    getLabelSpecificRule,
    getLegendRowRules,
    moveSpecificRuleDown: (index) => {
      const drawing = state.activeDrawing();
      return moveSpecificRule(drawing, index, 1);
    },
    moveSpecificRuleUp: (index) => {
      const drawing = state.activeDrawing();
      return moveSpecificRule(drawing, index, -1);
    },
    refreshFeatureOverrides,
    runWithRuleMatches,
    removeSpecificRule: (index) => {
      const drawing = state.activeDrawing();
      const sourceRows = drawing.manualSpecificRules.filter((_, i) => i !== index);
      return commitPrepared(drawing, sourceRows.map(rule => ({ ...rule })), 'Remove specific color rule', () => {}, null, sourceRows);
    },
    setSpecificRuleField
  };
};
