import {
  buildFeatureVisibilityOverrideCache,
  createDefaultFeatureVisibilityRule,
  featureVisibilityQualifierSuggestions,
  getEditorFeatureVisibilityMode,
  normalizeFeatureVisibilityRule,
  removeEditorFeatureVisibilityRule,
  serializeFeatureVisibilityRules,
  upsertEditorFeatureVisibilityRule
} from '../feature-visibility.js';

export const createFeatureVisibilityActions = ({ state, featureSvgActions }) => {
  const {
    clickedFeature,
    featureVisibilityRules,
    featureVisibilityOverrides,
    labelReflowForceRequestSeq,
    labelReflowForceRequestReason
  } = state;

  const { applyVisibilityPreviewBySvgId } = featureSvgActions;
  const ruleFields = new Set(['recordId', 'featureType', 'qualifier', 'value', 'action']);

  const normalizeVisibilityMode = (value) => {
    const normalized = String(value || '').trim().toLowerCase();
    return ['on', 'off', 'suppress'].includes(normalized) ? normalized : 'default';
  };

  const queueFeatureVisibilityReflow = (reason = 'feature-visibility-apply') => {
    labelReflowForceRequestReason.value = String(reason || 'feature-visibility-apply');
    labelReflowForceRequestSeq.value += 1;
  };

  const rebuildFeatureVisibilityOverrideCache = () => {
    const cache = buildFeatureVisibilityOverrideCache(featureVisibilityRules);
    Object.keys(featureVisibilityOverrides).forEach((key) => delete featureVisibilityOverrides[key]);
    Object.assign(featureVisibilityOverrides, cache);
  };

  const setFeatureVisibility = (feat, modeRaw, options = {}) => {
    const svgId = String(feat?.svg_id || '').trim();
    if (!svgId) return false;

    const triggerReflow = options.triggerReflow !== false;
    const nextMode = normalizeVisibilityMode(modeRaw);
    const previousMode = getEditorFeatureVisibilityMode(featureVisibilityRules, svgId);

    if (nextMode === 'default') {
      removeEditorFeatureVisibilityRule(featureVisibilityRules, svgId);
    } else {
      upsertEditorFeatureVisibilityRule(featureVisibilityRules, feat, nextMode);
    }
    rebuildFeatureVisibilityOverrideCache();

    if (clickedFeature.value?.svg_id === svgId) {
      clickedFeature.value.featureVisibility = nextMode;
    }

    applyVisibilityPreviewBySvgId(svgId, nextMode);

    if (triggerReflow && previousMode !== nextMode) {
      queueFeatureVisibilityReflow();
    }

    return previousMode !== nextMode;
  };

  const updateClickedFeatureVisibility = (modeRaw) => {
    if (!clickedFeature.value?.feat) return false;
    return setFeatureVisibility(clickedFeature.value.feat, modeRaw, { triggerReflow: true });
  };

  const getFeatureVisibility = (feat) => {
    const svgId = String(feat?.svg_id || '').trim();
    if (!svgId) return 'default';
    return getEditorFeatureVisibilityMode(featureVisibilityRules, svgId);
  };

  const setFeatureVisibilityRuleField = (index, field, value) => {
    if (!ruleFields.has(field)) return;
    const current = featureVisibilityRules[index];
    if (!current) return;
    const nextRule = normalizeFeatureVisibilityRule({ ...current, [field]: value });
    featureVisibilityRules.splice(index, 1, nextRule);
    rebuildFeatureVisibilityOverrideCache();
  };

  const moveFeatureVisibilityRule = (index, offset) => {
    const target = index + offset;
    if (target < 0 || target >= featureVisibilityRules.length) return;
    const [rule] = featureVisibilityRules.splice(index, 1);
    featureVisibilityRules.splice(target, 0, rule);
    rebuildFeatureVisibilityOverrideCache();
  };

  const addFeatureVisibilityRule = () => {
    featureVisibilityRules.push(createDefaultFeatureVisibilityRule());
    rebuildFeatureVisibilityOverrideCache();
  };

  const removeFeatureVisibilityRule = (index) => {
    if (index < 0 || index >= featureVisibilityRules.length) return;
    featureVisibilityRules.splice(index, 1);
    rebuildFeatureVisibilityOverrideCache();
  };

  const downloadTextFile = (filename, text) => {
    const blob = new Blob([text], { type: 'text/tab-separated-values;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    link.addEventListener('click', (event) => event.stopPropagation(), { once: true });
    link.click();
    URL.revokeObjectURL(url);
  };

  const downloadFeatureVisibilityRulesTsv = () => {
    const text = serializeFeatureVisibilityRules(featureVisibilityRules);
    if (!text.trim()) {
      alert('No valid feature visibility rules to export.');
      return;
    }
    downloadTextFile('gbdraw_feature_visibility_table.tsv', text);
  };

  const featureVisibilityRuleDetail = (rule) => {
    const normalized = normalizeFeatureVisibilityRule(rule);
    if (normalized.source === 'editor' && normalized.featureId) {
      return normalized.label ? `${normalized.label} (${normalized.featureId})` : normalized.featureId;
    }
    if (normalized.qualifier.toLowerCase() === 'hash') return normalized.value;
    return '';
  };

  return {
    addFeatureVisibilityRule,
    downloadFeatureVisibilityRulesTsv,
    featureVisibilityQualifierSuggestions,
    featureVisibilityRuleDetail,
    getFeatureVisibility,
    moveFeatureVisibilityRuleDown: (index) => moveFeatureVisibilityRule(index, 1),
    moveFeatureVisibilityRuleUp: (index) => moveFeatureVisibilityRule(index, -1),
    removeFeatureVisibilityRule,
    setFeatureVisibility,
    setFeatureVisibilityRuleField,
    updateClickedFeatureVisibility
  };
};
