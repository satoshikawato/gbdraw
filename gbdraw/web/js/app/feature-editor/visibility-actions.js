import {
  applyFeatureVisibilityOverrideChanges,
  buildExactQualifierFeatureVisibilityRule,
  buildFeatureVisibilityChanges,
  createDefaultFeatureVisibilityRule,
  featureDrawnContext,
  featureVisibilityQualifierSuggestions,
  getFeatureVisibilityOverride,
  normalizeFeatureVisibilityRule,
  normalizeVisibilityMode,
  removeEditorQualifierFeatureVisibilityRule,
  resolveFeatureDrawn,
  serializeFeatureVisibilityRules,
  setFeatureVisibilityOverride,
  upsertEditorQualifierFeatureVisibilityRule,
} from '../feature-visibility.js';
import {
  isInternalProteinDisplayId,
  resolveDisplayProteinId
} from '../feature-utils.js';
import { downloadTextFile } from '../../services/text-download.js';
import { resolveUniqueOrthogroupMemberForFeature } from '../../services/feature-identity.js';
import { featureIdentityKeyOf } from '../../services/feature-placement.js';
import { resultRenderedFeatures } from '../../services/feature-catalog.js';

export const createFeatureVisibilityActions = ({
  state,
  featureSvgActions,
  labelActions = null,
  previewRuntime = null,
  rulePreparation = null,
  getCommittedRequest = () => null
}) => {
  const {
    clickedFeature,
    extractedFeatures,
    orthogroups,
    featureVisibilityManualRules,
    featureVisibilityRules,
    featureOverrides,
    featureVisibilityScopeDialog,
    resultGenerationKey,
    results,
    selectedResultIndex,
    svgContainer
  } = state;

  const { applyVisibilityPreviewChanges } = featureSvgActions;
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

  const orthogroupIdFor = (source) => {
    const ids = new Set([
      source?.id,
      source?.orthogroupId,
      source?.orthogroup_id
    ].map(normalizeText).filter(Boolean));
    return ids.size === 1 ? ids.values().next().value : '';
  };

  const findOrthogroup = (orthogroupId) => {
    const id = normalizeText(orthogroupId);
    if (!id) return null;
    const matches = (Array.isArray(orthogroups?.value) ? orthogroups.value : [])
      .filter((group) => orthogroupIdFor(group) === id);
    return matches.length === 1 ? matches[0] : null;
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
    const orthogroupId = orthogroupIdFor({
      orthogroupId: feat?.orthogroupId,
      orthogroup_id: feat?.orthogroup_id
    });
    if (!orthogroupId) return [];
    const group = findOrthogroup(orthogroupId);
    if (!group) return [];
    const members = Array.isArray(group.members) ? group.members : [];
    if (!resolveUniqueOrthogroupMemberForFeature(feat, members)) return [];

    const featuresByMember = new Map();
    (Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []).forEach((candidate) => {
      const candidateGroupId = orthogroupIdFor({
        orthogroupId: candidate?.orthogroupId,
        orthogroup_id: candidate?.orthogroup_id
      });
      if (candidateGroupId !== orthogroupId) return;
      const member = resolveUniqueOrthogroupMemberForFeature(candidate, members);
      if (!member) return;
      const matches = featuresByMember.get(member) || [];
      matches.push(candidate);
      featuresByMember.set(member, matches);
    });
    return members.flatMap((member) => {
      const matches = featuresByMember.get(member) || [];
      return matches.length === 1 ? matches : [];
    });
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
        label: `Current similarity-group members (${orthogroupMembers.length})`,
        description: 'One rule per current similarity-group member.',
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
    const proteinIdCandidate = firstText(
      feat?.sourceProteinId,
      feat?.source_protein_id,
      getQualifierValue(feat, 'protein_id'),
      feat?.proteinId,
      feat?.protein_id
    );
    const proteinId = isInternalProteinDisplayId(proteinIdCandidate) ? '' : proteinIdCandidate;
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

  // The label owner hides a hidden feature's label with it, as Generate does,
  // and queues the label reflow unless the caller declines it (F-3).
  const applyFeatureVisibilityToLabels = (options = {}) => (
    labelActions?.applyFeatureVisibilityToLabels?.(options) ?? false
  );

  // The visibility table holds rules only; per-feature edits are identity rows.
  const visibilityRuleRows = () => (
    Array.isArray(featureVisibilityRules?.value)
      ? featureVisibilityRules.value
      : featureVisibilityManualRules.map((rule) => normalizeFeatureVisibilityRule(rule))
  );

  const featureSvgId = (feature) => normalizeText(feature?.svg_id ?? feature?.svgId ?? feature?.id);
  const sameFeature = (left, right) => {
    const key = featureIdentityKeyOf(left);
    return Boolean(key) && key === featureIdentityKeyOf(right);
  };

  const isRuleScope = (scope) => scope?.id === 'product' || scope?.id === 'protein_id';

  // The live projection of what the label rerender and Generate draw (R3): a
  // feature shows the resolver's answer, and an unknown answer leaves it as
  // it is. The resolver reads Python's rule matches, so an action prepares
  // them first (`prepareDrawn`).
  const drawnChanges = (features) => {
    const context = featureDrawnContext(state, { diagramOptions: getCommittedRequest()?.diagramOptions });
    return features.flatMap((feature) => {
      const featureId = featureSvgId(feature);
      const drawn = featureId ? resolveFeatureDrawn(feature, context) : null;
      return drawn === null ? [] : [{ featureId, mode: drawn ? 'on' : 'off' }];
    });
  };
  const displayedFeatures = () => {
    const displayed = resultRenderedFeatures(state);
    return displayed ? [...displayed.values()]
      : uniqueFeaturesBySvgId(Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []);
  };
  const prepareDrawn = async (rules = []) => {
    if (rules.length) await rulePreparation?.prepareVisibility?.(rules);
    await rulePreparation?.prepareDrawn?.();
  };

  const updateClickedFeatureVisibilityFromRules = (features) => {
    const clicked = clickedFeature.value?.feat;
    if (!clicked || !features.some((feature) => sameFeature(feature, clicked))) return;
    clickedFeature.value.featureVisibility = getFeatureVisibilityOverride(featureOverrides, clicked);
  };

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
    const targetFeatures = uniqueFeaturesBySvgId(Array.isArray(features) ? features : [])
      .filter((feature) => featureIdentityKeyOf(feature));
    if (targetFeatures.length === 0) return null;

    const changes = buildFeatureVisibilityChanges(targetFeatures, modeRaw, featureOverrides);
    if (changes.length === 0) return null;
    const changedKeys = new Set(changes.map(featureIdentityKeyOf));
    const changedFeatures = targetFeatures.filter((feature) => changedKeys.has(featureIdentityKeyOf(feature)));
    const commandResultIndex = Number(selectedResultIndex?.value || 0);
    const commandGenerationKey = String(resultGenerationKey?.value ?? '');

    const applyChangeSet = async (direction, reason) => {
      if (!(await ensureCommandTargetResult(commandResultIndex, commandGenerationKey))) return false;
      await prepareDrawn();
      const useAfter = direction === 'apply';
      const overrideChanged = changes.some((change) => change.before !== change.after);
      applyFeatureVisibilityOverrideChanges(
        featureOverrides,
        changes.map((change) => ({
          ...change,
          mode: useAfter ? change.after : change.before
        }))
      );
      const updated = applyVisibilityPreviewChanges(drawnChanges(changedFeatures), { reason });
      if (!updated && !overrideChanged) return false;
      updateClickedFeatureVisibilityFromRules(targetFeatures);
      applyFeatureVisibilityToLabels();
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

    if (isRuleScope(selectedScope)) {
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
        setFeatureVisibilityOverride(featureOverrides, targetFeat, nextMode);
      });
    }

    // A rule can reach any feature, so it projects onto every displayed one.
    applyVisibilityPreviewChanges(drawnChanges(isRuleScope(selectedScope) ? displayedFeatures() : targetFeatures));

    const clicked = clickedFeature.value?.feat;
    if (clicked && targetFeatures.some((feature) => sameFeature(feature, clicked))) {
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
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!featureIdentityKeyOf(feat)) return false;

    const triggerReflow = options.triggerReflow !== false;
    const scope = options.scope || { id: 'feature' };
    const nextMode = normalizeVisibilityMode(modeRaw);
    const previousMode = getFeatureVisibilityOverride(featureOverrides, feat);

    applyFeatureVisibilityScope(feat, nextMode, scope);

    if (previousMode !== nextMode) {
      applyFeatureVisibilityToLabels({ reflow: triggerReflow });
    }

    return previousMode !== nextMode;
  };

  const setSelectedFeaturesVisibility = async (features, modeRaw) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const command = buildSelectedFeaturesVisibilityCommand(features, modeRaw);
    if (!command) return false;
    return command.apply();
  };

  const updateClickedFeatureVisibility = async (modeRaw) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!clickedFeature.value?.feat) return false;
    const feat = clickedFeature.value.feat;
    const scopes = buildVisibilityScopes(feat);
    const nextMode = normalizeVisibilityMode(modeRaw);
    const previousMode = getFeatureVisibilityOverride(featureOverrides, feat);
    if (scopes.length <= 1) {
      await prepareDrawn();
      return setFeatureVisibility(feat, nextMode, { triggerReflow: true, scope: scopes[0] });
    }
    featureVisibilityScopeDialog.show = true;
    featureVisibilityScopeDialog.feat = feat;
    featureVisibilityScopeDialog.mode = nextMode;
    featureVisibilityScopeDialog.previousMode = previousMode;
    featureVisibilityScopeDialog.scopes = scopes;
    return false;
  };

  const handleFeatureVisibilityScopeChoice = async (scopeId) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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
    clearFeatureVisibilityScopeDialog();
    const rule = isRuleScope(scope) ? buildExactQualifierFeatureVisibilityRule({ ...scope, action: nextMode }) : null;
    await prepareDrawn(rule ? [rule] : []);
    applyFeatureVisibilityScope(feat, nextMode, scope);
    if (previousMode !== nextMode) applyFeatureVisibilityToLabels();
    return previousMode !== nextMode;
  };

  const getFeatureVisibility = (feat) => getFeatureVisibilityOverride(featureOverrides, feat);

  const setFeatureVisibilityRuleField = (index, field, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    featureVisibilityManualRules.push(createDefaultFeatureVisibilityRule());
  };

  const removeFeatureVisibilityRule = (index) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (index < 0 || index >= featureVisibilityManualRules.length) return;
    featureVisibilityManualRules.splice(index, 1);
  };

  const downloadFeatureVisibilityRulesTsv = () => {
    const text = serializeFeatureVisibilityRules(visibilityRuleRows());
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

  // Reconcile with the same resolver as the visibility action, so Undo and
  // Redo of a scoped hide leave the preview as the action left it (FE-04).
  // The displayed Result's features carry the identities its edits name (R3).
  const reconcileFeatureVisibility = () => applyVisibilityPreviewChanges(drawnChanges(displayedFeatures()));

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
    reconcileFeatureVisibility,
    removeFeatureVisibilityRule,
    setFeatureVisibility,
    setSelectedFeaturesVisibility,
    setFeatureVisibilityRuleField,
    updateClickedFeatureVisibility
  };
};
