import { createSpecificRulePatternDrafts } from './pattern-drafts.js';
import { normalizeUserFacingError } from '../../services/error-normalization.js';
import { ruleMatchesFeature, firstMatchingRule, ruleMatchesReady } from '../rule-matching.js';
import { resolveColorToHex } from '../color-utils.js';
import { parseSpecificRules, serializeSpecificRules } from '../file-imports.js';
import { getFeatureGenerationHash } from '../feature-utils.js';
import { legendRowRules, ruleLegendCaption } from '../specific-color-rules.js';
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

export const createFeatureRuleActions = ({ state, nextTick, legendActions, rulePreparation, history, svgActions, ref, computed, isPatternEditAvailable = () => true }) => {
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
    labelTextFeatureOverrides,
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
  const commitSpecificRules = async (rules, label = 'Change specific color rules', { isCurrent = () => true, afterCommit = () => {}, previousLegendIntents = [], sourceRows = rules.map(rule => manualSpecificRules.includes(rule) ? rule : null) } = {}) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    const revision = ++preparationRevision;
    const candidate = await rulePreparation.prepareCandidate(rules, {}, { retiredLegendIntents: previousLegendIntents });
    if (!candidate) return false;
    const current = () => revision === preparationRevision && !state.sessionOperationAvailability?.() && isCurrent()
      && rulePreparation.isCurrent(candidate.snapshot);
    if (!current()) return false;
    const previousIntents = [...candidate.previousIntents, ...previousLegendIntents];
    const previousCaptions = new Set(previousIntents.map(intent => intent.caption));
    let applied = false;
    await legendActions.syncFileLegendEntries(candidate.intents.filter(intent => !(state.deletedLegendEntries?.value || [])
      .some(entry => (entry.originalCaption || entry.caption) === intent.caption)), {
      previousFileIntents: previousIntents,
      isCurrent: current,
      transact: (diff, apply) => (diff.add.length || diff.remove.length
        ? history.runUndoableCheckpoint : history.runUndoable)(label, apply),
      commit: () => {
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
          ...candidate.intents.map(intent => intent.caption)
        ]);
        applyRulePreview();
        afterCommit(candidate);
        applied = true;
      }
    });
    if (applied) rulePreparation.notifyChanges(candidate);
    return applied;
  };
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
    svgActions.applyPaletteToSvg();
    svgActions.applySpecificRulesToSvg();
  };

  const editSpecificRulePattern = (row, value) => {
    const busy = state.sessionOperationAvailability?.();
    return busy || patternDrafts.edit(row, value);
  };
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

    const featureIdKey = normalizeFeatureIdKey(feat.svg_id || feat.id);
    if (featureIdKey) {
      for (const [overrideFeatureId, overrideText] of Object.entries(labelTextFeatureOverrides)) {
        if (normalizeFeatureIdKey(overrideFeatureId) !== featureIdKey) continue;
        const normalizedOverride = normalizeCaption(overrideText);
        if (normalizedOverride) return normalizedOverride;
        break;
      }
    }

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
    const generationHash = getFeatureGenerationHash(feat);
    return generationHash ? { qual: 'hash', val: generationHash } : null;
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

  const findFeaturesWithSameLegendItem = (currentFeat, caption = null) => {
    const targetCaption = normalizeCaption(caption || getEffectiveLegendCaption(currentFeat));
    if (!targetCaption) return [];
    return extractedFeatures.value.filter((f) => {
      if (f.svg_id === currentFeat.svg_id) return false;
      return captionMatches(getEffectiveLegendCaption(f), targetCaption);
    });
  };

  const findFeaturesWithSameIndividualLabel = (currentFeat, label = null) => {
    const targetLabel = normalizeCaption(label || getIndividualFeatureLabel(currentFeat));
    if (!targetLabel) return [];

    return extractedFeatures.value.filter((f) => {
      if (f.svg_id === currentFeat.svg_id) return false;
      return captionMatches(getIndividualFeatureLabel(f), targetLabel);
    });
  };

  const findFeaturesWithSameDisplayedLabel = (currentFeat, label = null) => {
    const targetLabel = normalizeCaption(label || getDisplayedFeatureLabel(currentFeat));
    if (!targetLabel) return [];

    return extractedFeatures.value.filter((f) => {
      if (f.svg_id === currentFeat.svg_id) return false;
      return captionMatches(getDisplayedFeatureLabel(f), targetLabel);
    });
  };

  const findExistingColorForCaption = (currentFeat, caption) => {
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
    removeSpecificRule: (index) => {
      const sourceRows = manualSpecificRules.filter((_, i) => i !== index);
      return commitPrepared(sourceRows.map(rule => ({ ...rule })), 'Remove specific color rule', () => {}, null, sourceRows);
    },
    setSpecificRuleField
  };
};
