import { resolveColorToHex } from '../color-utils.js';
import { parseSpecificRules, serializeSpecificRules } from '../file-imports.js';
import {
  getFeatureGenerationHash,
  resolveOrderedSpecificRule,
  ruleMatchesFeature
} from '../feature-utils.js';
import { resolveFeatureLabelSelector } from '../feature-selector.js';
import { downloadTextFile } from '../../services/text-download.js';
import {
  getFeatureOverride
} from '../../services/feature-override-identity.js';
import {
  defaultFeatureRendering,
  normalizeFeatureRendering
} from '../../utils/feature-rendering.js';
import { resolveFeatureFillViewModel } from './fill-view-model.js';
import { FEATURE_INSTANCE_HASH_QUALIFIER } from '../../services/feature-instance-identity.js';
import {
  FEATURE_SEMANTIC_SCOPE_QUALIFIER,
  isFeatureSemanticScopeRule
} from './semantic-fill-selectors.js';

export const createFeatureRuleActions = ({
  state,
  nextTick,
  legendActions,
  bulkStyleActions = null
}) => {
  const {
    appliedPaletteColors,
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
    fileLegendCaptions
  } = state;

  void nextTick;
  void legendActions;
  const normalizeCaption = (value) => String(value || '').trim();
  const normalizeCaptionKey = (value) => normalizeCaption(value).toLowerCase();
  const normalizeFeatureIdKey = (value) => String(value || '').trim().toLowerCase();
  const captionMatches = (value, target) => normalizeCaptionKey(value) === normalizeCaptionKey(target);
  const specificRuleFields = new Set(['feat', 'qual', 'val', 'color', 'cap']);
  const isExactFeatureRule = (rule) => (
    String(rule?.qual || '').trim() === FEATURE_INSTANCE_HASH_QUALIFIER
    && String(rule?.match || '').trim().toLowerCase() !== 'regex'
  );
  const isOpaqueSpecificRule = (rule) => (
    isExactFeatureRule(rule) || isFeatureSemanticScopeRule(rule)
  );
  const canMoveSpecificRule = (index, offset) => {
    const target = index + offset;
    return Number.isInteger(index)
      && Number.isInteger(offset)
      && target >= 0
      && target < manualSpecificRules.length
      && !isOpaqueSpecificRule(manualSpecificRules[index])
      && !isOpaqueSpecificRule(manualSpecificRules[target]);
  };

  const commitSpecificRules = async (rules, {
    writerKind,
    label,
    appliedPaletteColors: nextPaletteColors = undefined,
    presentationPatch = null
  } = {}) => {
    if (typeof bulkStyleActions?.requestFeatureBulkStyleChange !== 'function') {
      throw new Error('Bulk Feature style actions are unavailable.');
    }
    return bulkStyleActions.requestFeatureBulkStyleChange({
      writerKind,
      label,
      presentationPatch,
      replacement: {
        rules,
        ...(nextPaletteColors === undefined
          ? {}
          : { appliedPaletteColors: nextPaletteColors })
      }
    });
  };

  const syncFileCaptionTracking = () => {
    fileLegendCaptions.value = new Set(
      manualSpecificRules
        .filter((rule) => rule?.fromFile && normalizeCaption(rule?.cap))
        .map((rule) => normalizeCaption(rule.cap))
    );
  };

  const setSpecificRuleField = async (index, field, value) => {
    if (!specificRuleFields.has(field)) return;
    const current = manualSpecificRules[index];
    if (!current) return;
    if (isOpaqueSpecificRule(current)) return false;
    const nextValue = field === 'color' ? resolveColorToHex(String(value || '#000000')) : String(value ?? '');

    if (field === 'val') {
      try {
        new RegExp(nextValue);
      } catch (e) {
        alert('Invalid Regular Expression: ' + e.message);
        return;
      }
    }

    const nextRule = { ...current, [field]: nextValue };
    delete nextRule.fromFile;
    const nextRules = manualSpecificRules.map((rule, ruleIndex) => (
      ruleIndex === index ? nextRule : rule
    ));
    const committed = await commitSpecificRules(nextRules, {
      writerKind: 'rule-field',
      label: 'Edit specific color rule'
    });
    if (committed) syncFileCaptionTracking();
    return committed;
  };

  const moveSpecificRule = async (index, offset) => {
    const target = index + offset;
    if (!canMoveSpecificRule(index, offset)) return false;
    const nextRules = [...manualSpecificRules];
    const [rule] = nextRules.splice(index, 1);
    nextRules.splice(target, 0, rule);
    return commitSpecificRules(nextRules, {
      writerKind: 'rule-reorder',
      label: 'Reorder specific color rules'
    });
  };

  const removeSpecificRule = async (index) => {
    if (!manualSpecificRules[index] || isOpaqueSpecificRule(manualSpecificRules[index])) return false;
    const nextRules = manualSpecificRules.filter((_rule, ruleIndex) => ruleIndex !== index);
    const committed = await commitSpecificRules(nextRules, {
      writerKind: 'rule-remove',
      label: 'Remove specific color rule'
    });
    if (committed) syncFileCaptionTracking();
    return committed;
  };

  const downloadSpecificRulesTsv = () => {
    const text = serializeSpecificRules(manualSpecificRules);
    if (!text.trim()) {
      alert('No specific rules to export.');
      return;
    }
    if (
      manualSpecificRules.some(isOpaqueSpecificRule)
      && typeof confirm === 'function'
      && !confirm(
        'Managed Feature rows require gbdraw 0.14.0b0 or later. Exact rows only '
        + 'replay against the same canonical request/session record keys. Download this TSV?'
      )
    ) {
      return false;
    }
    downloadTextFile('gbdraw_specific_rules.tsv', text);
    return true;
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
    const resolvedCaption = normalizeCaption(
      resolveOrderedSpecificRule(feat, manualSpecificRules)?.rule?.cap
    );
    if (resolvedCaption) return resolvedCaption;

    return normalizeCaption(feat.type) || normalizeCaption(getIndividualFeatureLabel(feat));
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
    if ([
      FEATURE_INSTANCE_HASH_QUALIFIER,
      FEATURE_SEMANTIC_SCOPE_QUALIFIER
    ].includes(String(newSpecRule.qual || '').trim())) {
      alert('Managed Feature selectors cannot be added in the generic rule editor.');
      return false;
    }

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

    try {
      new RegExp(newSpecRule.val);
    } catch (e) {
      alert('Invalid Regular Expression: ' + e.message);
      return;
    }

    const rule = {
      feat: String(newSpecRule.feat || ''),
      qual: String(newSpecRule.qual || ''),
      val: String(newSpecRule.val),
      color: String(newSpecRule.color || '#000000'),
      cap: String(newSpecRule.cap || '')
    };
    const committed = await commitSpecificRules([...manualSpecificRules, rule], {
      writerKind: 'rule-add',
      label: 'Add specific color rule'
    });
    if (committed) newSpecRule.val = '';
    return committed;
  };

  const clearAllSpecificRules = async () => {
    const retainedRules = manualSpecificRules.filter(isOpaqueSpecificRule);
    if (retainedRules.length === manualSpecificRules.length) return false;
    const committed = await commitSpecificRules(retainedRules, {
      writerKind: 'rule-clear',
      label: 'Clear editable specific color rules'
    });
    if (committed) syncFileCaptionTracking();
    return committed;
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

    specificRulePresetLoading.value = true;
    try {
      const response = await fetch(preset.path, { cache: 'no-store' });
      if (!response.ok) {
        throw new Error(`Preset fetch failed: ${response.status}`);
      }
      const text = await response.text();
      const { rules } = parseSpecificRules(text);
      const retainedOpaqueRules = manualSpecificRules.filter(isOpaqueSpecificRule);
      const presetColors = presetId === 'bakta'
        ? { ...appliedPaletteColors.value, CDS: '#cccccc' }
        : appliedPaletteColors.value;
      const committed = await commitSpecificRules([...retainedOpaqueRules, ...rules], {
        writerKind: 'preset',
        label: 'Apply specific color preset',
        appliedPaletteColors: presetColors,
        presentationPatch: presetId === 'bakta'
          ? { legendBoxSize: 12, legendFontSize: 12 }
          : null
      });
      if (committed) {
        syncFileCaptionTracking();
      }
      return committed;
    } catch (e) {
      console.error('Failed to load specific rule preset:', e);
      alert('Failed to load preset. Please check the preset file and format.');
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
    const svgDefaultColor = feat?.fill_color
      || feat?.fillColor
      || getFeatureOverride(featureColorOverrides, feat)?.color
      || getFeatureOverride(featureColorOverrides, feat)
      || '#cccccc';
    return resolveFeatureFillViewModel({
      feature: feat,
      specificRules: manualSpecificRules,
      paletteColors: appliedPaletteColors.value,
      svgDefaultColor
    }).effectiveColor;
  };

  const getFeatureColorValue = (feat) => {
    const rule = resolveOrderedSpecificRule(feat, manualSpecificRules)?.rule;
    return isExactFeatureRule(rule) ? rule.color : null;
  };

  const canEditFeatureColor = () => true;

  const getFeatureFillViewModel = (feat) => resolveFeatureFillViewModel({
    feature: feat,
    specificRules: manualSpecificRules,
    paletteColors: appliedPaletteColors.value,
    svgDefaultColor: feat?.fill_color || feat?.fillColor || '#cccccc'
  });

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

  const findMatchingRegexRule = (feat) => {
    const rule = resolveOrderedSpecificRule(feat, manualSpecificRules)?.rule || null;
    return rule && !isOpaqueSpecificRule(rule) ? rule : null;
  };

  const countFeaturesMatchingRule = (rule) => {
    if (!rule || String(rule.qual || '').toLowerCase() === 'hash') return 0;

    let count = 0;
    for (const feat of extractedFeatures.value) {
      if (feat.type !== rule.feat) continue;

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
    addFeature,
    removeFeature,
    getFeatureShape,
    setFeatureShape,
    addPriorityRule,
    addSpecificRule,
    applySpecificRulePreset,
    canMoveSpecificRule,
    canEditFeatureColor,
    clearAllSpecificRules,
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
    getFeatureFillViewModel,
    getDisplayedFeatureLabel,
    getEffectiveLegendCaption,
    getIndividualFeatureLabel,
    getFeatureQualifier,
    getLabelSpecificRule,
    isOpaqueSpecificRule,
    moveSpecificRuleDown: (index) => moveSpecificRule(index, 1),
    moveSpecificRuleUp: (index) => moveSpecificRule(index, -1),
    removeSpecificRule,
    setSpecificRuleField
  };
};
