import {
  applyFeatureVisibilityOverrideChanges,
  buildFeatureVisibilityChanges,
  createDefaultFeatureVisibilityRule,
  deriveFeatureVisibilityRulesForBoundary,
  featureVisibilityQualifierSuggestions,
  getFeatureVisibilityOverride,
  normalizeFeatureVisibilityRule,
  normalizeVisibilityMode,
  removeEditorQualifierFeatureVisibilityRule,
  resolveEffectiveFeatureVisibility,
  serializeFeatureVisibilityRules,
  setFeatureVisibilityOverride,
  upsertEditorQualifierFeatureVisibilityRule,
} from '../feature-visibility.js';
import { downloadTextFile } from '../../services/text-download.js';

export const createFeatureVisibilityActions = ({ state, featureSvgActions, previewRuntime = null }) => {
  const {
    clickedFeature,
    extractedFeatures,
    orthogroups,
    featureVisibilityManualRules,
    featureVisibilityRules,
    featureVisibilityOverrides,
    featureVisibilitySelectorCache,
    featureVisibilityScopeDialog,
    labelLayoutDirtyReason,
    resultGenerationKey,
    results,
    selectedResultIndex,
    svgContainer
  } = state;

  const { applyVisibilityPreviewBySvgId, applyVisibilityPreviewChanges } = featureSvgActions;
  const ruleFields = new Set(['recordId', 'featureType', 'qualifier', 'value', 'action']);

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

  const markFeatureVisibilityLabelLayoutDirty = (reason = 'feature-visibility') => {
    if (labelLayoutDirtyReason) labelLayoutDirtyReason.value = String(reason || 'feature-visibility');
  };

  const boundaryFeatureVisibilityRules = () => (
    Array.isArray(featureVisibilityRules?.value)
      ? featureVisibilityRules.value
      : deriveFeatureVisibilityRulesForBoundary(
          featureVisibilityManualRules,
          featureVisibilityOverrides,
          featureVisibilitySelectorCache
        )
  );

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
      const effectiveMode = specificHashMode ||
        (mode === 'default'
          ? resolveEffectiveFeatureVisibility(svgId, featureVisibilityOverrides, null, featureVisibilityManualRules)
          : mode);
      applyVisibilityPreviewBySvgId(svgId, effectiveMode);
    });
  };

  const updateClickedFeatureVisibilityFromRules = (featureIds) => {
    const clickedSvgId = normalizeText(clickedFeature.value?.svg_id);
    if (!clickedSvgId || !featureIds.includes(clickedSvgId)) return;
    clickedFeature.value.featureVisibility = getFeatureVisibilityOverride(featureVisibilityOverrides, clickedSvgId);
  };

  const cloneOverrides = () => ({ ...(featureVisibilityOverrides || {}) });

  const resolveEffectiveWithOverrides = (featureId, overrides) =>
    resolveEffectiveFeatureVisibility(featureId, overrides, null, featureVisibilityManualRules);

  const nextFrame = () => new Promise((resolve) => {
    if (typeof window !== 'undefined' && typeof window.requestAnimationFrame === 'function') {
      window.requestAnimationFrame(() => resolve());
    } else {
      setTimeout(resolve, 0);
    }
  });

  const ensureCommandTargetResult = async (resultIndex, generationKey) => {
    if (String(resultGenerationKey?.value ?? '') !== String(generationKey ?? '')) return false;
    const targetIndex = Number(resultIndex);
    if (!Number.isInteger(targetIndex) || targetIndex < 0) return false;
    const resultCount = Array.isArray(results?.value) ? results.value.length : 0;
    if (targetIndex >= resultCount) return false;
    if (Number(selectedResultIndex?.value || 0) !== targetIndex) {
      if (previewRuntime?.selectResult) {
        previewRuntime.selectResult(targetIndex);
      } else if (selectedResultIndex) {
        selectedResultIndex.value = targetIndex;
      }
      await nextFrame();
    }
    return Boolean(svgContainer?.value?.querySelector?.('svg'));
  };

  const buildSelectedFeaturesVisibilityCommand = (features, modeRaw) => {
    const targetFeatures = uniqueFeaturesBySvgId(Array.isArray(features) ? features : []);
    if (targetFeatures.length === 0) return null;

    const affectedFeatureIds = targetFeatures
      .map((feature) => normalizeText(feature?.svg_id ?? feature?.svgId ?? feature?.id))
      .filter(Boolean);
    if (affectedFeatureIds.length === 0) return null;

    const baseOverrides = cloneOverrides();
    const rawChanges = buildFeatureVisibilityChanges(targetFeatures, modeRaw, baseOverrides);
    if (rawChanges.length === 0) return null;
    const afterOverrides = cloneOverrides();
    applyFeatureVisibilityOverrideChanges(
      afterOverrides,
      rawChanges.map((change) => ({ featureId: change.featureId, mode: change.after }))
    );
    const changes = rawChanges.map((change) => ({
      ...change,
      beforeEffective: resolveEffectiveWithOverrides(change.featureId, baseOverrides),
      afterEffective: resolveEffectiveWithOverrides(change.featureId, afterOverrides)
    }));
    const commandResultIndex = Number(selectedResultIndex?.value || 0);
    const commandGenerationKey = String(resultGenerationKey?.value ?? '');

    const applyChangeSet = async (direction, reason) => {
      if (!(await ensureCommandTargetResult(commandResultIndex, commandGenerationKey))) return false;
      const useAfter = direction === 'apply';
      const overrideChanged = changes.some((change) => change.before !== change.after);
      applyFeatureVisibilityOverrideChanges(
        featureVisibilityOverrides,
        changes.map((change) => ({
          featureId: change.featureId,
          mode: useAfter ? change.after : change.before
        }))
      );
      const updated = applyVisibilityPreviewChanges(
        changes.map((change) => ({
          featureId: change.featureId,
          mode: useAfter ? change.afterEffective : change.beforeEffective
        })),
        { reason }
      );
      if (!updated && !overrideChanged) return false;
      if (updated && previewRuntime?.flushActiveResult) {
        previewRuntime.flushActiveResult({ force: true });
      }
      updateClickedFeatureVisibilityFromRules(affectedFeatureIds);
      markFeatureVisibilityLabelLayoutDirty(reason);
      return true;
    };

    return {
      label: 'Change selected feature visibility',
      resultIndex: commandResultIndex,
      resultGenerationKey: commandGenerationKey,
      changes,
      apply: () => applyChangeSet('apply', 'bulk-feature-visibility-apply'),
      revert: () => applyChangeSet('revert', 'bulk-feature-visibility-undo'),
      estimateBytes: () => JSON.stringify({ resultIndex: commandResultIndex, resultGenerationKey: commandGenerationKey, changes }).length * 2
    };
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
        removeEditorQualifierFeatureVisibilityRule(featureVisibilityManualRules, ruleInput);
      } else {
        upsertEditorQualifierFeatureVisibilityRule(featureVisibilityManualRules, ruleInput, nextMode);
      }
    } else {
      targetFeatures.forEach((targetFeat) => {
        const svgId = normalizeText(targetFeat?.svg_id ?? targetFeat?.svgId ?? targetFeat?.id);
        if (!svgId) return;
        setFeatureVisibilityOverride(featureVisibilityOverrides, svgId, nextMode);
      });
    }

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
    const previousMode = getFeatureVisibilityOverride(featureVisibilityOverrides, svgId);

    applyFeatureVisibilityScope(feat, nextMode, scope);

    if (triggerReflow && previousMode !== nextMode) {
      markFeatureVisibilityLabelLayoutDirty();
    }

    return previousMode !== nextMode;
  };

  const setSelectedFeaturesVisibility = async (features, modeRaw) => {
    const command = buildSelectedFeaturesVisibilityCommand(features, modeRaw);
    if (!command) return false;
    return command.apply();
  };

  const updateClickedFeatureVisibility = (modeRaw) => {
    if (!clickedFeature.value?.feat) return false;
    const feat = clickedFeature.value.feat;
    const scopes = buildVisibilityScopes(feat);
    const nextMode = normalizeVisibilityMode(modeRaw);
    const previousMode = getFeatureVisibilityOverride(featureVisibilityOverrides, normalizeText(feat.svg_id));
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
    if (previousMode !== nextMode) markFeatureVisibilityLabelLayoutDirty();
    return previousMode !== nextMode;
  };

  const getFeatureVisibility = (feat) => {
    const svgId = String(feat?.svg_id || '').trim();
    if (!svgId) return 'default';
    return getFeatureVisibilityOverride(featureVisibilityOverrides, svgId);
  };

  const setFeatureVisibilityRuleField = (index, field, value) => {
    if (!ruleFields.has(field)) return;
    const current = featureVisibilityManualRules[index];
    if (!current) return;
    const nextRule = normalizeFeatureVisibilityRule({ ...current, [field]: value });
    featureVisibilityManualRules.splice(index, 1, nextRule);
  };

  const moveFeatureVisibilityRule = (index, offset) => {
    const target = index + offset;
    if (target < 0 || target >= featureVisibilityManualRules.length) return;
    const [rule] = featureVisibilityManualRules.splice(index, 1);
    featureVisibilityManualRules.splice(target, 0, rule);
  };

  const addFeatureVisibilityRule = () => {
    featureVisibilityManualRules.push(createDefaultFeatureVisibilityRule());
  };

  const removeFeatureVisibilityRule = (index) => {
    if (index < 0 || index >= featureVisibilityManualRules.length) return;
    featureVisibilityManualRules.splice(index, 1);
  };

  const downloadFeatureVisibilityRulesTsv = () => {
    const text = serializeFeatureVisibilityRules(boundaryFeatureVisibilityRules());
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
    buildSelectedFeaturesVisibilityCommand,
    getFeatureVisibility,
    handleFeatureVisibilityScopeChoice,
    moveFeatureVisibilityRuleDown: (index) => moveFeatureVisibilityRule(index, 1),
    moveFeatureVisibilityRuleUp: (index) => moveFeatureVisibilityRule(index, -1),
    removeFeatureVisibilityRule,
    setFeatureVisibility,
    setSelectedFeaturesVisibility,
    setFeatureVisibilityRuleField,
    updateClickedFeatureVisibility
  };
};
