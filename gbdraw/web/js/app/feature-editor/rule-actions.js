// @ts-check
import { createSpecificRulePatternDrafts } from './pattern-drafts.js';
import { normalizeUserFacingError } from '../../services/error-normalization.js';
import { ruleMatchesFeature, firstMatchingRule, ruleMatchesReady, runWhenPrepared } from '../rule-matching.js';
import { resolveColorToHex } from '../color-utils.js';
import { parseSpecificRules, serializeSpecificRules } from '../file-imports.js';
import { getFeatureColorRuleHash } from '../feature-utils.js';
import {
  buildLegendIntents, createRuleLegendCaptions, legendRowRules, rendererLegendRows, ruleLegendCaption
} from '../specific-color-rules.js';
import { resolveFeatureLabelSelector } from '../feature-selector.js';
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
import { featureDrawnContext, resultLegendSources, sameLegendSources } from '../feature-visibility.js';

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
 * The Legend rows a rule commit draws, prepared and not yet applied.
 * @typedef {object} PreparedFileLegend
 * @property {{ add: any[], update: any[], remove: any[], unchanged: any[] }} diff
 * @property {() => boolean} isCurrent
 * @property {() => void} apply
 */

/**
 * @typedef {object} FeatureRuleActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(intents: Record<string, any>[], options?: { previousFileIntents?: Record<string, any>[], isCurrent?: () => boolean }) => Promise<PreparedFileLegend | false>} prepareFileLegendEntries
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
    currentColors,
    appliedPaletteColors,
    newColorFeat,
    newColorVal,
    manualSpecificRules,
    newSpecRule,
    specificRulePresets,
    selectedSpecificPreset,
    specificRulePresetLoading,
    manualPriorityRules,
    newPriorityRule,
    adv,
    newFeatureToAdd,
    extractedFeatures,
    featureColorOverrides,
    editableLabels,
    featureOverrides,
    labelTextBulkOverrides,
    addedLegendCaptions,
    fileLegendCaptions
  } = state;

  const normalizeCaption = (value) => String(value || '').trim();
  const normalizeCaptionKey = (value) => normalizeCaption(value).toLowerCase();
  const normalizeFeatureIdKey = (value) => String(value || '').trim().toLowerCase();
  const captionMatches = (value, target) => normalizeCaptionKey(value) === normalizeCaptionKey(target);
  const specificRuleFields = new Set(['feat', 'qual', 'val', 'color', 'cap']);

  let preparationRevision = 0;
  const patternDrafts = createSpecificRulePatternDrafts({
    rules: manualSpecificRules, ref, invalidate: () => { preparationRevision += 1; }
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
  // are no renderer rows for the N-06 caption allocation.
  const candidateLegendIntents = (candidateRules, retiredLegendIntents) => {
    const rendered = (extractedFeatures.value || []).filter(feature =>
      featureOverrideValue(featureOverrides, feature, 'featureVisibility') !== 'off');
    const used = new Set(rendered.map(feature => firstMatchingRule(feature, candidateRules)).filter(Boolean));
    const rendererRows = rendererLegendRows({
      legendEntries: state.legendEntries?.value,
      originalLegendOrder: state.originalLegendOrder?.value,
      rules: [...manualSpecificRules, ...candidateRules,
        ...retiredLegendIntents.map(intent => ({ cap: intent?.caption, color: intent?.color }))]
    });
    const currentCaption = createRuleLegendCaptions(manualSpecificRules, rendererRows);
    return {
      intents: buildLegendIntents(candidateRules.filter(rule => used.has(rule)), rendererRows).intents,
      previousIntents: [...manualSpecificRules.filter(rule => rule.cap)
        .map(rule => ({ caption: currentCaption(rule), color: rule.color })), ...retiredLegendIntents]
    };
  };
  // OV-43 (Owner decision 2026-10-06, option A): whether a rule table change
  // changes a Result's Legend source, read with the current Feature
  // visibility before and after. Python redraws the Legend rows, their order,
  // and the "other <type>s" rows in the automatic rerender.
  const legendSourceContext = (rules) => ({
    ...featureDrawnContext(state, { diagramOptions: getCommittedRequest()?.diagramOptions }),
    colorRules: rules
  });
  const changesLegendSource = (before, after) => !sameLegendSources(
    resultLegendSources(state, legendSourceContext(before)),
    resultLegendSources(state, legendSourceContext(after))
  );
  // Undo and Redo of a rule edit restore the rules; the composition root
  // passes the rules they replaced, and a changed Legend source asks for the
  // rerender, as the edit did.
  const followRestoredRules = (previousRules) => (
    JSON.stringify(previousRules) !== JSON.stringify(manualSpecificRules)
    && changesLegendSource(previousRules, manualSpecificRules)
    && ports.requestAutomaticRerender()
  );
  // The automatic rerender replaces the Results the candidate was prepared
  // against, which makes the candidate stale (#857). An edit made while one runs
  // waits for it and, when it replaced the Results under an unchanged rule
  // table, prepares again instead of dropping the edit.
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
  const retireSupersededLegendColors = (intents) => {
    const overrides = state.legendColorOverrides;
    for (const { caption, color } of intents) {
      if (Object.hasOwn(overrides, caption) && resolveColorToHex(overrides[caption]) !== resolveColorToHex(color)) {
        delete overrides[caption];
      }
    }
  };
  /**
   * @typedef {{ isCurrent?: () => boolean, afterCommit?: (intents: Record<string, any>[]) => void, previousLegendIntents?: Record<string, any>[], sourceRows?: (Record<string, any> | null)[] }} CommitSpecificRulesOptions
   */
  /**
   * @param {Record<string, any>[]} rules
   * @param {string} [label]
   * @param {CommitSpecificRulesOptions} [options]
   */
  const commitSpecificRules = async (rules, label = 'Change specific color rules', options = {}) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    await rerenderIdle();
    for (let attempt = 0; ; attempt += 1) {
      const results = state.results.value;
      const table = JSON.stringify(manualSpecificRules);
      const outcome = await commitOnce(rules, label, options);
      if (outcome || attempt === 2 || state.results.value === results
        || JSON.stringify(manualSpecificRules) !== table) return outcome;
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
   * @param {Record<string, any>[]} rules
   * @param {string} label
   * @param {CommitSpecificRulesOptions} [options]
   */
  const commitOnce = async (rules, label, { isCurrent = () => true, afterCommit = () => {}, previousLegendIntents = [], sourceRows = rules.map(rule => manualSpecificRules.includes(rule) ? rule : null) } = {}) => {
    const revision = ++preparationRevision;
    const candidate = await prepareRules(rules);
    if (!candidate) return false;
    const current = () => revision === preparationRevision && !state.sessionOperationAvailability?.() && isCurrent()
      && rulePreparation.isCurrent(candidate.snapshot);
    if (!current()) return false;
    const { intents, previousIntents } = candidateLegendIntents(candidate.rules, previousLegendIntents);
    const previousCaptions = new Set(previousIntents.map(intent => intent.caption));
    const legend = await prepareFileLegendEntries(intents.filter(intent => !(state.deletedLegendEntries?.value || [])
      .some(entry => (entry.originalCaption || entry.caption) === intent.caption)), {
      previousFileIntents: previousIntents,
      isCurrent: current
    });
    if (!legend) return false;
    const redrawsLegend = changesLegendSource([...manualSpecificRules], candidate.rules);
    let applied = false;
    // One History step: the rule transition first, then the legend rows it
    // draws (R13); a checkpoint when the legend gains or loses a row.
    const transact = legend.diff.add.length || legend.diff.remove.length
      ? runUndoableCheckpoint : runUndoable;
    await transact(label, () => {
      if (!current() || !legend.isCurrent()) return false;
      manualSpecificRules.splice(0, manualSpecificRules.length, ...candidate.rules.map((rule, index) => {
        const row = sourceRows[index];
        if (!row) return rule;
        Object.assign(row, rule);
        if (!Object.hasOwn(rule, 'fromFile')) delete row.fromFile;
        return row;
      }));
      patternDrafts.reconcile();
      fileLegendCaptions.value = new Set(candidate.rules.filter(rule => rule.fromFile && rule.cap).map(rule => rule.cap));
      addedLegendCaptions.value = new Set([
        ...[...addedLegendCaptions.value].filter(caption => !previousCaptions.has(caption)),
        ...intents.map(intent => intent.caption)
      ]);
      applyRulePreview();
      retireSupersededLegendColors(intents);
      afterCommit(intents);
      applied = true;
      legend.apply();
      return legend.diff;
    });
    if (applied) rulePreparation.notifyChanges(candidate);
    if (applied && redrawsLegend) ports.requestAutomaticRerender();
    return applied;
  };
  /** @param {Element | null} [input] The input a failed commit offers for editing. */
  const commitPrepared = async (rules, label, afterCommit = () => {}, input = null, sourceRows = rules.map(rule => manualSpecificRules.includes(rule) ? rule : null)) => {
    patternDrafts.suspend();
    const snapshot = rulePreparation.snapshot();
    const previousError = state.errorLog?.value;
    try {
      const applied = await commitSpecificRules(rules, label, { afterCommit, sourceRows });
      if (applied && state.errorLog?.value === ruleFailure.value?.error) state.errorLog.value = null;
      if (applied) ruleFailure.value = null;
      return applied;
    } catch (cause) {
      if (!rulePreparation.isCurrent(snapshot) || state.errorLog?.value !== previousError) return false;
      const error = normalizeUserFacingError(cause, { operation: 'evaluateRules', stage: 'resource-staging' });
      if (state.errorLog) state.errorLog.value = error;
      ruleFailure.value = { error, snapshot, input, retry: () => commitPrepared(rules, label, afterCommit, input, sourceRows) };
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
  const applySpecificRulePattern = async (row, value) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    if (!manualSpecificRules.includes(row) || !isPatternEditAvailable()
      || (state.generatedMode?.value && state.mode?.value !== state.generatedMode.value)) return false;
    const token = patternDrafts.begin(row, value);
    const snapshot = rulePreparation.snapshot();
    const attemptRevision = preparationRevision + 1;
    const current = () => preparationRevision === attemptRevision && patternDrafts.isCurrent(row, token) && rulePreparation.isCurrent(snapshot)
      && isPatternEditAvailable();
    const sourceRows = [...manualSpecificRules];
    const rules = sourceRows.map(rule => {
      if (rule !== row) return { ...rule };
      const next = { ...rule, val: String(value ?? '') };
      delete next.fromFile;
      return next;
    });
    try {
      const applied = await commitSpecificRules(rules, 'Edit specific color rule', { isCurrent: current, sourceRows });
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
    const draft = patternDrafts.get(row);
    return draft && !draft.pending ? applySpecificRulePattern(row, draft.text) : false;
  };
  /** @param {HTMLInputElement | null} [input] */
  const setSpecificRuleField = (index, field, value, input = null) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!specificRuleFields.has(field)) return;
    const current = manualSpecificRules[index];
    if (!current) return;
    if (field === 'val') return applySpecificRulePattern(current, value);
    const nextValue = field === 'color' ? resolveColorToHex(String(value || '#000000')) : String(value ?? '');
    const nextRule = { ...current, [field]: nextValue };
    delete nextRule.fromFile;
    const sourceRows = [...manualSpecificRules];
    const rules = sourceRows.map(rule => rule === current ? nextRule : { ...rule });
    return commitPrepared(rules, 'Edit specific color rule', () => {}, input, sourceRows)?.finally(() => {
      if (input?.isConnected && manualSpecificRules.includes(current)) input.value = current[field] ?? '';
    });
  };

  const moveSpecificRule = (index, offset) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const sourceRows = [...manualSpecificRules];
    const target = index + offset;
    if (target < 0 || target >= sourceRows.length) return;
    const [row] = sourceRows.splice(index, 1);
    sourceRows.splice(target, 0, row);
    return commitPrepared(sourceRows.map(rule => ({ ...rule })), 'Move specific color rule', () => {}, null, sourceRows);
  };

  const downloadSpecificRulesTsv = () => {
    const text = serializeSpecificRules(manualSpecificRules);
    if (!text.trim()) {
      alert('No specific rules to export.');
      return;
    }
    downloadTextFile('gbdraw_specific_rules.tsv', text);
  };

  const getIndividualFeatureLabel = (feat) => {
    return feat.product || feat.gene || feat.locus_tag || `${feat.type} at ${feat.start}..${feat.end}`;
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
    if (!feat) return '';

    const editableEntry = getEditableLabelEntryForFeature(feat);
    const editableText = normalizeCaption(editableEntry?.text);
    if (editableText) return editableText;

    const normalizedOverride = normalizeCaption(featureOverrideValue(featureOverrides, feat, 'labelText'));
    if (normalizedOverride) return normalizedOverride;

    const sourceText = normalizeCaption(editableEntry?.sourceText);
    if (sourceText) {
      const normalizedBulk = normalizeCaption(labelTextBulkOverrides[sourceText]);
      if (normalizedBulk) return normalizedBulk;
    }

    return normalizeCaption(getIndividualFeatureLabel(feat));
  };

  const legendRowContext = () => ({
    rules: manualSpecificRules,
    legendEntries: state.legendEntries?.value || [],
    originalLegendOrder: state.originalLegendOrder?.value || []
  });
  // The rules a legend row draws; editing the row edits them (N-06).
  const getLegendRowRules = (caption) => legendRowRules(caption, legendRowContext());

  // Resolve the effective legend item label used by current SVG coloring
  // priority: a rule's feature belongs to the row Generate draws for it (N-06).
  const getEffectiveLegendCaption = (feat) => {
    if (!feat) return '';

    const rule = firstMatchingRule(feat, manualSpecificRules);
    if (rule && normalizeCaption(rule.cap)) return ruleLegendCaption(rule, legendRowContext());

    const overrideCaption = normalizeCaption(
      getFeatureOverride(featureColorOverrides, feat)?.caption
    );
    if (overrideCaption) return overrideCaption;

    return normalizeCaption(feat.type) || normalizeCaption(getIndividualFeatureLabel(feat));
  };

  const addCustomColor = () => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!newColorFeat.value) return;
    currentColors.value = {
      ...currentColors.value,
      [newColorFeat.value]: newColorVal.value
    };
  };

  const setLabelFilterMode = (value) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    state.filterMode.value = value;
  };
  const addWhitelistRule = () => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    state.manualWhitelist.push({ feat: 'CDS', qual: 'product', key: '' });
  };
  const removeWhitelistRule = (index) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    state.manualWhitelist.splice(index, 1);
  };
  const removePriorityRule = (index) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    manualPriorityRules.splice(index, 1);
  };

  const addPriorityRule = () => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!newPriorityRule.order) return;
    const idx = manualPriorityRules.findIndex((r) => r.feat === newPriorityRule.feat);
    if (idx >= 0) {
      manualPriorityRules[idx].order = newPriorityRule.order;
    } else {
      manualPriorityRules.push({ feat: newPriorityRule.feat, order: newPriorityRule.order });
    }
  };

  const addSpecificRule = async () => {
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
    return commitPrepared([...manualSpecificRules, rule], 'Add specific color rule', () => {
      if (newSpecRule.val === rule.val) newSpecRule.val = '';
    }, document.querySelector?.('[data-new-specific-rule-pattern]'));
  };

  const applySpecificRulePreset = async () => {
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
        isCurrent: () => selectedSpecificPreset.value === presetId,
        afterCommit: () => {
          if (presetId === 'bakta') {
            currentColors.value = { ...currentColors.value, CDS: '#cccccc' };
            adv.legend_box_size = 12;
            adv.legend_font_size = 12;
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
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (newFeatureToAdd.value && !adv.features.includes(newFeatureToAdd.value)) {
      adv.features.push(newFeatureToAdd.value);
      if (!adv.feature_shapes || typeof adv.feature_shapes !== 'object') {
        adv.feature_shapes = {};
      }
      if (!Object.prototype.hasOwnProperty.call(adv.feature_shapes, newFeatureToAdd.value)) {
        adv.feature_shapes[newFeatureToAdd.value] = defaultFeatureRendering(newFeatureToAdd.value);
      }
    }
  };

  const removeFeature = (featureType) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const idx = adv.features.indexOf(featureType);
    if (idx >= 0) {
      adv.features.splice(idx, 1);
    }
  };

  const getFeatureShape = (featureType) => {
    if (!adv.feature_shapes || typeof adv.feature_shapes !== 'object') {
      return defaultFeatureRendering(featureType);
    }
    return Object.prototype.hasOwnProperty.call(adv.feature_shapes, featureType)
      ? normalizeFeatureRendering(adv.feature_shapes[featureType])
      : defaultFeatureRendering(featureType);
  };

  const setFeatureShape = (featureType, shape) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!adv.feature_shapes || typeof adv.feature_shapes !== 'object') {
      adv.feature_shapes = {};
    }
    adv.feature_shapes[featureType] = normalizeFeatureRendering(shape);
  };

  const getFeatureColor = (feat) => {
    const override = getFeatureOverride(featureColorOverrides, feat);
    if (override) {
      return resolveColorToHex(override.color || override);
    }
    return resolveColorToHex(appliedPaletteColors.value[feat.type]) || '#cccccc';
  };

  const getFeatureColorValue = (feat) => {
    const override = getFeatureOverride(featureColorOverrides, feat);
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
    if (!feat) return null;
    const priorityRule = manualPriorityRules.find((rule) => rule.feat === feat.type);
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
    if (!features || features.length === 0 || !ruleMatchesReady(features, manualSpecificRules)) return;
    migrateLegacyFeatureOverrides(featureColorOverrides, features);

    for (const feat of features) {
      const rule = firstMatchingRule(feat, manualSpecificRules);
      const key = featureOverrideKey(feat);
      if (key) {
        if (rule) featureColorOverrides[key] = { color: rule.color, caption: rule.cap };
        else delete featureColorOverrides[key];
      }
    }
  };

  const findMatchingRegexRule = (feat) => {
    return firstMatchingRule(feat, manualSpecificRules.filter((rule) => rule.qual !== 'hash'));
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
    const targetCaption = normalizeCaption(caption || getEffectiveLegendCaption(currentFeat));
    if (!targetCaption) return [];
    return extractedFeatures.value.filter((f) => {
      if (f.svg_id === currentFeat.svg_id) return false;
      return captionMatches(getEffectiveLegendCaption(f), targetCaption);
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
    const targetCaption = normalizeCaption(caption);
    if (!targetCaption) return null;

    for (const rule of manualSpecificRules) {
      if (captionMatches(rule.cap, targetCaption) && String(rule.qual || '').toLowerCase() === 'hash') {
        return { rule, color: rule.color };
      }
    }

    for (const rule of manualSpecificRules) {
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
    clearAllSpecificRules: () => commitPrepared([], 'Clear specific color rules'),
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
    getEffectiveLegendCaption,
    getIndividualFeatureLabel,
    getFeatureQualifier,
    getLabelSpecificRule,
    getLegendRowRules,
    moveSpecificRuleDown: (index) => moveSpecificRule(index, 1),
    moveSpecificRuleUp: (index) => moveSpecificRule(index, -1),
    refreshFeatureOverrides,
    runWithRuleMatches,
    removeSpecificRule: (index) => {
      const sourceRows = manualSpecificRules.filter((_, i) => i !== index);
      return commitPrepared(sourceRows.map(rule => ({ ...rule })), 'Remove specific color rule', () => {}, null, sourceRows);
    },
    setSpecificRuleField
  };
};
