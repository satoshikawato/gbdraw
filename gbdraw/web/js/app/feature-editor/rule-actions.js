import { resolveColorToHex } from '../color-utils.js';
import { parseSpecificRules } from '../file-imports.js';
import { ruleMatchesFeature } from '../feature-utils.js';

export const createFeatureRuleActions = ({ state, nextTick, legendActions }) => {
  const {
    pyodideReady,
    currentColors,
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
    addedLegendCaptions,
    fileLegendCaptions
  } = state;

  const { addLegendEntry, removeLegendEntry, extractLegendEntries } = legendActions;
  const normalizeFeatureShape = (value) => (String(value || '').trim().toLowerCase() === 'arrow' ? 'arrow' : 'rectangle');

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
    manualSpecificRules.push(rule);

    if (rule.cap && pyodideReady.value) {
      await nextTick();
      const actualCaption = await addLegendEntry(rule.cap, rule.color);
      if (actualCaption && typeof actualCaption === 'string') {
        addedLegendCaptions.value.add(actualCaption);
      }
      extractLegendEntries();
    }

    newSpecRule.val = '';
  };

  const clearAllSpecificRules = async () => {
    const captionsToRemove = manualSpecificRules.filter((rule) => rule.cap).map((rule) => rule.cap);

    for (const cap of captionsToRemove) {
      await removeLegendEntry(cap);
    }

    manualSpecificRules.splice(0);
    extractLegendEntries();
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
      const { rules, rulesWithCaptions } = parseSpecificRules(text);

      const captionsToRemove = manualSpecificRules.filter((rule) => rule.cap).map((rule) => rule.cap);
      for (const cap of captionsToRemove) {
        await removeLegendEntry(cap);
        addedLegendCaptions.value.delete(cap);
        fileLegendCaptions.value.delete(cap);
      }

      manualSpecificRules.splice(0);
      rules.forEach((rule) => manualSpecificRules.push(rule));

      if (presetId === 'bakta') {
        currentColors.value = { ...currentColors.value, CDS: '#cccccc' };
        adv.legend_box_size = 12;
        adv.legend_font_size = 12;
      }

      if (rulesWithCaptions.length > 0 && pyodideReady.value) {
        await nextTick();
        for (const rule of rulesWithCaptions) {
          const actualCaption = await addLegendEntry(rule.cap, rule.color);
          if (actualCaption && typeof actualCaption === 'string') {
            addedLegendCaptions.value.add(actualCaption);
            fileLegendCaptions.value.add(actualCaption);
          }
        }
        extractLegendEntries();
      } else {
        extractLegendEntries();
      }
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
        adv.feature_shapes[newFeatureToAdd.value] = 'rectangle';
      }
    }
  };

  const removeFeature = (featureType) => {
    const idx = adv.features.indexOf(featureType);
    if (idx >= 0) {
      adv.features.splice(idx, 1);
    }
    if (adv.feature_shapes && typeof adv.feature_shapes === 'object') {
      delete adv.feature_shapes[featureType];
    }
  };

  const getFeatureShape = (featureType) => {
    if (!adv.feature_shapes || typeof adv.feature_shapes !== 'object') {
      return 'rectangle';
    }
    return normalizeFeatureShape(adv.feature_shapes[featureType]);
  };

  const setFeatureShape = (featureType, shape) => {
    if (!adv.feature_shapes || typeof adv.feature_shapes !== 'object') {
      adv.feature_shapes = {};
    }
    adv.feature_shapes[featureType] = normalizeFeatureShape(shape);
  };

  const getFeatureColor = (feat) => {
    const override = featureColorOverrides[feat.id];
    if (override) {
      return resolveColorToHex(override.color || override);
    }
    return resolveColorToHex(currentColors.value[feat.type]) || '#cccccc';
  };

  const canEditFeatureColor = () => true;

  const getFeatureQualifier = (feat) => {
    return { qual: 'hash', val: feat.svg_id };
  };

  const refreshFeatureOverrides = (features) => {
    Object.keys(featureColorOverrides).forEach((k) => delete featureColorOverrides[k]);
    if (!features || features.length === 0) return;

    for (const feat of features) {
      for (const rule of manualSpecificRules) {
        if (!ruleMatchesFeature(feat, rule)) continue;
        featureColorOverrides[feat.id] = { color: rule.color, caption: rule.cap };
        break;
      }
    }
  };

  const findMatchingRegexRule = (feat) => {
    for (const rule of manualSpecificRules) {
      if (rule.feat !== feat.type) continue;
      if (rule.qual === 'hash') continue;

      if (ruleMatchesFeature(feat, rule)) {
        return rule;
      }
    }
    return null;
  };

  const countFeaturesMatchingRule = (rule) => {
    if (!rule || rule.qual === 'hash') return 0;

    let count = 0;
    for (const feat of extractedFeatures.value) {
      if (feat.type !== rule.feat) continue;

      if (ruleMatchesFeature(feat, rule)) count++;
    }
    return count;
  };

  const findFeaturesWithSameCaption = (currentFeat, caption) => {
    if (!caption) return [];
    return extractedFeatures.value.filter((f) => {
      if (f.svg_id === currentFeat.svg_id) return false;
      const featCaption = f.product || f.gene || f.locus_tag || f.type;
      return featCaption === caption;
    });
  };

  const findExistingColorForCaption = (currentFeat, caption) => {
    if (!caption) return null;

    for (const rule of manualSpecificRules) {
      if (rule.cap === caption && rule.qual === 'hash') {
        return { rule, color: rule.color };
      }
    }

    for (const rule of manualSpecificRules) {
      if (rule.cap === caption && rule.qual !== 'hash') {
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
    clearAllSpecificRules,
    countFeaturesMatchingRule,
    findExistingColorForCaption,
    findFeaturesWithSameCaption,
    findMatchingRegexRule,
    getFeatureColor,
    getFeatureQualifier,
    refreshFeatureOverrides
  };
};
