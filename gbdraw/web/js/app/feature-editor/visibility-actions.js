import {
  buildFeatureVisibilityOverrideCache,
  createDefaultFeatureVisibilityRule,
  featureVisibilityQualifierSuggestions,
  getEditorFeatureVisibilityMode,
  normalizeFeatureVisibilityRule,
  removeEditorQualifierFeatureVisibilityRule,
  removeEditorFeatureVisibilityRule,
  serializeFeatureVisibilityRules,
  upsertEditorQualifierFeatureVisibilityRule,
  upsertEditorFeatureVisibilityRule
} from '../feature-visibility.js';

export const createFeatureVisibilityActions = ({ state, featureSvgActions }) => {
  const {
    clickedFeature,
    extractedFeatures,
    featureSelectorSafetyScope,
    orthogroups,
    featureVisibilityRules,
    featureVisibilityOverrides,
    featureVisibilityScopeDialog,
    labelReflowForceRequestSeq,
    labelReflowForceRequestReason
  } = state;

  const { applyVisibilityPreviewBySvgId } = featureSvgActions;
  const ruleFields = new Set(['recordId', 'featureType', 'qualifier', 'value', 'action']);

  const normalizeVisibilityMode = (value) => {
    const normalized = String(value || '').trim().toLowerCase();
    if (normalized === 'suppress') return 'exclude_matching';
    return ['on', 'off', 'exclude_matching'].includes(normalized) ? normalized : 'default';
  };

  const normalizeText = (value) => String(value ?? '').trim();

  const firstText = (...values) => {
    for (const value of values) {
      if (Array.isArray(value)) {
        const found = firstText(...value);
        if (found) return found;
        continue;
      }
      const text = normalizeText(value);
      if (text) return text;
    }
    return '';
  };

  const getFeatureType = (feat) => normalizeText(feat?.type || feat?.featureType || feat?.feature_type);

  const getQualifierValue = (feat, key) => {
    const normalizedKey = normalizeText(key).toLowerCase();
    if (!normalizedKey) return '';
    const qualifiers = feat?.qualifiers && typeof feat.qualifiers === 'object' ? feat.qualifiers : {};
    return firstText(feat?.[normalizedKey], qualifiers[normalizedKey]);
  };

  const memberFeatureSvgId = (member) => firstText(member?.featureSvgId, member?.feature_svg_id);

  const memberRecordIndex = (member) => {
    const recordIndex = Number(member?.recordIndex ?? member?.record_index);
    return Number.isInteger(recordIndex) ? recordIndex : null;
  };

  const featureRecordIndex = (feat) => {
    const recordIndex = Number(feat?.record_idx ?? feat?.recordIndex ?? feat?.record_index ?? feat?.fileIdx);
    return Number.isInteger(recordIndex) ? recordIndex : null;
  };

  const buildFeatureLookup = () => {
    const bySvgId = new Map();
    const byRecordAndSvgId = new Map();
    (Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []).forEach((feat) => {
      const svgId = normalizeText(feat?.svg_id ?? feat?.svgId);
      if (!svgId) return;
      if (!bySvgId.has(svgId)) bySvgId.set(svgId, feat);
      const recordIndex = featureRecordIndex(feat);
      if (recordIndex !== null) byRecordAndSvgId.set(`${recordIndex}:${svgId}`, feat);
    });
    return { bySvgId, byRecordAndSvgId };
  };

  const getFeatureForMember = (member, lookup) => {
    const svgId = memberFeatureSvgId(member);
    if (!svgId) return null;
    const recordIndex = memberRecordIndex(member);
    if (recordIndex !== null) {
      const matched = lookup.byRecordAndSvgId.get(`${recordIndex}:${svgId}`);
      if (matched) return matched;
    }
    const direct = lookup.bySvgId.get(svgId);
    if (direct) return direct;
    return {
      id: svgId,
      svg_id: svgId,
      label: firstText(member?.product, member?.proteinId, member?.sourceProteinId, svgId),
      type: 'CDS',
      product: firstText(member?.product),
      protein_id: firstText(member?.proteinId, member?.sourceProteinId)
    };
  };

  const findOrthogroup = (orthogroupId) => {
    const id = normalizeText(orthogroupId);
    if (!id) return null;
    return (Array.isArray(orthogroups?.value) ? orthogroups.value : [])
      .find((group) => normalizeText(group?.id || group?.orthogroupId || group?.orthogroup_id) === id) || null;
  };

  const uniqueFeaturesBySvgId = (features) => {
    const seen = new Set();
    return features.filter((feat) => {
      const svgId = normalizeText(feat?.svg_id ?? feat?.svgId ?? feat?.id);
      if (!svgId || seen.has(svgId)) return false;
      seen.add(svgId);
      return true;
    });
  };

  const getOrthogroupMemberFeatures = (feat) => {
    const orthogroupId = normalizeText(feat?.orthogroupId || feat?.orthogroup_id);
    if (!orthogroupId) return [];
    const lookup = buildFeatureLookup();
    const group = findOrthogroup(orthogroupId);
    if (group) {
      return uniqueFeaturesBySvgId(
        (Array.isArray(group.members) ? group.members : [])
          .map((member) => getFeatureForMember(member, lookup))
          .filter(Boolean)
      );
    }
    return uniqueFeaturesBySvgId(
      (Array.isArray(extractedFeatures.value) ? extractedFeatures.value : [])
        .filter((candidate) => normalizeText(candidate?.orthogroupId || candidate?.orthogroup_id) === orthogroupId)
    );
  };

  const getMatchingQualifierFeatures = ({ featureType, qualifier, value }) => {
    const normalizedType = normalizeText(featureType);
    const normalizedQualifier = normalizeText(qualifier).toLowerCase();
    const normalizedValue = normalizeText(value).toLowerCase();
    if (!normalizedType || !normalizedQualifier || !normalizedValue) return [];
    return (Array.isArray(extractedFeatures.value) ? extractedFeatures.value : [])
      .filter((feat) => getFeatureType(feat) === normalizedType)
      .filter((feat) => getQualifierValue(feat, normalizedQualifier).toLowerCase() === normalizedValue);
  };

  const buildVisibilityScopes = (feat) => {
    const featureType = getFeatureType(feat);
    const scopes = [{
      id: 'feature',
      label: 'This feature',
      description: 'One rule for only this feature.'
    }];
    const orthogroupMembers = getOrthogroupMemberFeatures(feat);
    if (orthogroupMembers.length > 1) {
      scopes.push({
        id: 'orthogroup',
        label: `Current orthogroup members (${orthogroupMembers.length})`,
        description: 'One rule per current orthogroup member.',
        features: orthogroupMembers
      });
    }
    const product = getQualifierValue(feat, 'product');
    if (featureType && product) {
      scopes.push({
        id: 'product',
        label: `Exact product: ${product}`,
        description: 'Exact product qualifier rule.',
        featureType,
        qualifier: 'product',
        value: product
      });
    }
    const proteinId = firstText(
      feat?.proteinId,
      feat?.protein_id,
      feat?.sourceProteinId,
      feat?.source_protein_id,
      getQualifierValue(feat, 'protein_id')
    );
    if (featureType && proteinId) {
      scopes.push({
        id: 'protein_id',
        label: `Exact protein ID: ${proteinId}`,
        description: 'Exact protein_id qualifier rule.',
        featureType,
        qualifier: 'protein_id',
        value: proteinId
      });
    }
    return scopes;
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

  const collectAffectedFeatureIds = (scope, feat) => {
    if (scope?.id === 'orthogroup') {
      return uniqueFeaturesBySvgId(scope.features || [])
        .map((member) => normalizeText(member?.svg_id ?? member?.svgId ?? member?.id))
        .filter(Boolean);
    }
    if (scope?.id === 'product' || scope?.id === 'protein_id') {
      return uniqueFeaturesBySvgId(getMatchingQualifierFeatures(scope))
        .map((candidate) => normalizeText(candidate?.svg_id ?? candidate?.svgId ?? candidate?.id))
        .filter(Boolean);
    }
    const svgId = normalizeText(feat?.svg_id ?? feat?.svgId ?? feat?.id);
    return svgId ? [svgId] : [];
  };

  const applyVisibilityPreviewForScope = (scope, feat, mode) => {
    collectAffectedFeatureIds(scope, feat).forEach((svgId) => {
      const specificHashMode = (scope?.id === 'product' || scope?.id === 'protein_id')
        ? featureVisibilityOverrides[svgId]
        : '';
      applyVisibilityPreviewBySvgId(svgId, specificHashMode || mode);
    });
  };

  const applyFeatureVisibilityScope = (feat, modeRaw, scope) => {
    const nextMode = normalizeVisibilityMode(modeRaw);
    const selectedScope = scope || { id: 'feature' };
    const targetFeatures = selectedScope.id === 'orthogroup'
      ? uniqueFeaturesBySvgId(selectedScope.features || [])
      : [feat];

    if (selectedScope.id === 'product' || selectedScope.id === 'protein_id') {
      const ruleInput = {
        featureType: selectedScope.featureType,
        qualifier: selectedScope.qualifier,
        value: selectedScope.value,
        label: selectedScope.label
      };
      if (nextMode === 'default') {
        removeEditorQualifierFeatureVisibilityRule(featureVisibilityRules, ruleInput);
      } else {
        upsertEditorQualifierFeatureVisibilityRule(featureVisibilityRules, ruleInput, nextMode);
      }
    } else {
      const selectorContext = {
        selectorSafetyScope: featureSelectorSafetyScope?.value
      };
      targetFeatures.forEach((targetFeat) => {
        const svgId = normalizeText(targetFeat?.svg_id ?? targetFeat?.svgId ?? targetFeat?.id);
        if (!svgId) return;
        if (nextMode === 'default') {
          removeEditorFeatureVisibilityRule(featureVisibilityRules, svgId);
        } else {
          upsertEditorFeatureVisibilityRule(featureVisibilityRules, targetFeat, nextMode, selectorContext);
        }
      });
    }

    rebuildFeatureVisibilityOverrideCache();
    applyVisibilityPreviewForScope(selectedScope, feat, nextMode);

    const clickedSvgId = normalizeText(clickedFeature.value?.svg_id);
    if (clickedSvgId && collectAffectedFeatureIds(selectedScope, feat).includes(clickedSvgId)) {
      clickedFeature.value.featureVisibility = nextMode;
    }
  };

  const clearFeatureVisibilityScopeDialog = ({ restorePrevious = false } = {}) => {
    if (restorePrevious && clickedFeature.value) {
      clickedFeature.value.featureVisibility = featureVisibilityScopeDialog.previousMode || 'default';
    }
    featureVisibilityScopeDialog.show = false;
    featureVisibilityScopeDialog.feat = null;
    featureVisibilityScopeDialog.mode = 'default';
    featureVisibilityScopeDialog.previousMode = 'default';
    featureVisibilityScopeDialog.scopes = [];
  };

  const setFeatureVisibility = (feat, modeRaw, options = {}) => {
    const svgId = String(feat?.svg_id || '').trim();
    if (!svgId) return false;

    const triggerReflow = options.triggerReflow !== false;
    const scope = options.scope || { id: 'feature' };
    const nextMode = normalizeVisibilityMode(modeRaw);
    const previousMode = getEditorFeatureVisibilityMode(featureVisibilityRules, svgId);

    applyFeatureVisibilityScope(feat, nextMode, scope);

    if (triggerReflow && previousMode !== nextMode) {
      queueFeatureVisibilityReflow();
    }

    return previousMode !== nextMode;
  };

  const updateClickedFeatureVisibility = (modeRaw) => {
    if (!clickedFeature.value?.feat) return false;
    const feat = clickedFeature.value.feat;
    const scopes = buildVisibilityScopes(feat);
    const nextMode = normalizeVisibilityMode(modeRaw);
    const previousMode = getEditorFeatureVisibilityMode(featureVisibilityRules, normalizeText(feat.svg_id));
    if (scopes.length <= 1) {
      return setFeatureVisibility(feat, nextMode, { triggerReflow: true, scope: scopes[0] });
    }
    featureVisibilityScopeDialog.show = true;
    featureVisibilityScopeDialog.feat = feat;
    featureVisibilityScopeDialog.mode = nextMode;
    featureVisibilityScopeDialog.previousMode = previousMode;
    featureVisibilityScopeDialog.scopes = scopes;
    return false;
  };

  const handleFeatureVisibilityScopeChoice = (scopeId) => {
    if (scopeId === 'cancel' || !featureVisibilityScopeDialog.show) {
      clearFeatureVisibilityScopeDialog({ restorePrevious: true });
      return false;
    }
    const feat = featureVisibilityScopeDialog.feat;
    const nextMode = featureVisibilityScopeDialog.mode;
    const scope = (featureVisibilityScopeDialog.scopes || []).find((entry) => entry.id === scopeId);
    if (!feat || !scope) {
      clearFeatureVisibilityScopeDialog({ restorePrevious: true });
      return false;
    }
    const previousMode = featureVisibilityScopeDialog.previousMode;
    applyFeatureVisibilityScope(feat, nextMode, scope);
    clearFeatureVisibilityScopeDialog();
    if (previousMode !== nextMode) queueFeatureVisibilityReflow();
    return previousMode !== nextMode;
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
      const qualifier = normalized.qualifier.toLowerCase() === 'hash'
        ? 'hash fallback'
        : normalized.qualifier.toLowerCase();
      return normalized.label ? `${normalized.label} (${qualifier})` : `${normalized.featureId} (${qualifier})`;
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
    handleFeatureVisibilityScopeChoice,
    moveFeatureVisibilityRuleDown: (index) => moveFeatureVisibilityRule(index, 1),
    moveFeatureVisibilityRuleUp: (index) => moveFeatureVisibilityRule(index, -1),
    removeFeatureVisibilityRule,
    setFeatureVisibility,
    setFeatureVisibilityRuleField,
    updateClickedFeatureVisibility
  };
};
