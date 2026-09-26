import { ruleMatchesFeature, firstMatchingRule, ruleMatchesReady } from '../rule-matching.js';
import { resolveColorToHex } from '../color-utils.js';
import { parseSpecificRules, serializeSpecificRules } from '../file-imports.js';
import { getFeatureGenerationHash } from '../feature-utils.js';
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

export const createFeatureRuleActions = ({ state, nextTick, legendActions, rulePreparation, history, svgActions }) => {
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
  const commitSpecificRules = async (rules, label = 'Change specific color rules', { isCurrent = () => true, afterCommit = () => {}, previousLegendIntents = [] } = {}) => {
    const revision = ++preparationRevision;
    const candidate = await rulePreparation.prepareCandidate(rules);
    if (!candidate) return false;
    const current = () => revision === preparationRevision && isCurrent()
      && rulePreparation.isCurrent(candidate.snapshot);
    if (!current()) return false;
    const previousIntents = [...manualSpecificRules.filter(rule => rule.cap)
      .map(rule => ({ caption: rule.cap, color: rule.color })), ...previousLegendIntents];
    const previousCaptions = new Set(previousIntents.map(intent => intent.caption));
    let applied = false;
    await history.runUndoableCheckpoint(label, async () => {
      if (!current()) return;
      await legendActions.syncFileLegendEntries(candidate.intents.filter(intent => !(state.deletedLegendEntries?.value || [])
        .some(entry => (entry.originalCaption || entry.caption) === intent.caption)), {
        previousFileIntents: previousIntents,
        isCurrent: current,
        commit: () => {
          manualSpecificRules.splice(0, manualSpecificRules.length, ...candidate.rules);
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
    });
    if (applied) rulePreparation.notifyChanges(candidate);
    return applied;
  };
  const commitPrepared = (rules, label, afterCommit = () => {}) => commitSpecificRules(rules, label, { afterCommit })
    .catch(error => { alert(`Invalid rule: ${error.message}`); return false; });
  const applyRulePreview = () => {
    refreshFeatureOverrides(extractedFeatures.value);
    svgActions.applyPaletteToSvg();
    svgActions.applySpecificRulesToSvg();
  };

  const setSpecificRuleField = (index, field, value, input = null) => {
    if (!specificRuleFields.has(field)) return;
    const current = manualSpecificRules[index];
    if (!current) return;
    const nextValue = field === 'color' ? resolveColorToHex(String(value || '#000000')) : String(value ?? '');

    const nextRule = { ...current, [field]: nextValue };
    delete nextRule.fromFile;
    const rules = manualSpecificRules.map((rule, i) => i === index ? nextRule : { ...rule });
    return commitPrepared(rules, 'Edit specific color rule')?.finally(() => {
      if (input?.isConnected) input.value = manualSpecificRules[index]?.[field] ?? '';
    });
  };

  const moveSpecificRule = (index, offset) => {
    const rules = manualSpecificRules.map(rule => ({ ...rule }));
    const target = index + offset;
    if (target < 0 || target >= rules.length) return;
    const [rule] = rules.splice(index, 1);
    rules.splice(target, 0, rule);
    return commitPrepared(rules, 'Move specific color rule');
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

  // Resolve the effective legend item label used by current SVG coloring priority.
  const getEffectiveLegendCaption = (feat) => {
    if (!feat) return '';

    const rule = firstMatchingRule(feat, manualSpecificRules);
    if (rule && normalizeCaption(rule.cap)) return normalizeCaption(rule.cap);

    const overrideCaption = normalizeCaption(
      getFeatureOverride(featureColorOverrides, feat)?.caption
    );
    if (overrideCaption) return overrideCaption;

    return normalizeCaption(feat.type) || normalizeCaption(getIndividualFeatureLabel(feat));
  };

  const addCustomColor = () => {
    if (!newColorFeat.value) return;
    currentColors.value = {
      ...currentColors.value,
      [newColorFeat.value]: newColorVal.value
    };
  };

  const addPriorityRule = () => {
    if (!newPriorityRule.order) return;
    const idx = manualPriorityRules.findIndex((r) => r.feat === newPriorityRule.feat);
    if (idx >= 0) {
      manualPriorityRules[idx].order = newPriorityRule.order;
    } else {
      manualPriorityRules.push({ feat: newPriorityRule.feat, order: newPriorityRule.order });
    }
  };

  const addSpecificRule = async () => {
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
    });
  };

  const applySpecificRulePreset = async () => {
    if (specificRulePresetLoading.value) return;
    const presetId = selectedSpecificPreset.value;
    if (!presetId) return;
    const preset = specificRulePresets.find((entry) => entry.id === presetId);
    if (!preset) {
      alert('Unknown preset selected.');
      return;
    }

    const presetContext = rulePreparation.snapshot();
    specificRulePresetLoading.value = true;
    try {
      const response = await fetch(preset.path, { cache: 'no-store' });
      if (!response.ok) {
        throw new Error(`Preset fetch failed: ${response.status}`);
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
      console.error('Failed to load specific rule preset:', e);
      alert(`Invalid rule: ${e.message}`);
    } finally {
      specificRulePresetLoading.value = false;
    }
  };

  const addFeature = () => {
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

  const getFeatureQualifier = (feat) => {
    const generationHash = getFeatureGenerationHash(feat);
    if (!generationHash) return null;
    const collisionCount = extractedFeatures.value.filter(
      (candidate) => candidate?.type === feat?.type && getFeatureGenerationHash(candidate) === generationHash
    ).length;
    const renderedId = String(feat?.svg_id || '').trim();
    // Preserve the rendered instance when duplicate records share one generation hash.
    const value = collisionCount > 1 && renderedId ? renderedId : generationHash;
    return { qual: 'hash', val: value };
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
    addCustomColor,
    addFeature,
    removeFeature,
    getFeatureShape,
    setFeatureShape,
    addPriorityRule,
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
    moveSpecificRuleDown: (index) => moveSpecificRule(index, 1),
    moveSpecificRuleUp: (index) => moveSpecificRule(index, -1),
    refreshFeatureOverrides,
    removeSpecificRule: (index) => commitPrepared(manualSpecificRules.filter((_, i) => i !== index), 'Remove specific color rule'),
    setSpecificRuleField
  };
};
