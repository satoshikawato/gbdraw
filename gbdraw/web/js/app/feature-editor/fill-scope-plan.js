import {
  FEATURE_INSTANCE_HASH_QUALIFIER,
  normalizeFeatureFillValue,
  normalizeRuleColor,
  resolveFeatureFillViewModel,
  resolveOrderedSpecificRule,
  specificRuleMatchesFeature
} from './fill-view-model.js';
import { isFeatureInstanceHash } from '../../services/feature-instance-identity.js';
import { catalogResultKey } from '../../services/feature-catalog.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
import {
  getFeatureRenderedLabel,
  getFeatureSimilarityGroupIds,
  getFeatureSourceAnnotation,
  isFeatureSemanticScopeRule
} from './semantic-fill-selectors.js';

export const FEATURE_FILL_SCOPE_ORDER = Object.freeze([
  'matching-rule',
  'feature-type',
  'legend-caption',
  'rendered-label',
  'source-annotation-label',
  'similarity-group',
  'selected-features',
  'single'
]);

export const FEATURE_FILL_RESULT_EXTENTS = Object.freeze([
  'current-result',
  'all-affected-results'
]);

const EXACT_SCOPES = new Set(['selected-features', 'single']);
const GROUP_SCOPES = new Set([
  'matching-rule',
  'feature-type',
  'legend-caption',
  'rendered-label',
  'source-annotation-label',
  'similarity-group'
]);
const VALID_SOURCES = new Set([
  'popup',
  'feature-list',
  'selection-toolbar',
  'legend-editor'
]);
const GENERATE_EXACT_DIAGNOSTIC = 'Generate to enable exact feature edits';

const text = (value) => String(value ?? '').trim();
const normalizeCaption = (value) => text(value);
const captionKey = (value) => normalizeCaption(value).toLowerCase();

export const normalizeFeatureFillSemanticScope = (value) => {
  const normalized = text(value);
  return FEATURE_FILL_SCOPE_ORDER.includes(normalized) ? normalized : null;
};

export const normalizeFeatureFillResultExtent = (value) => {
  const normalized = text(value);
  return FEATURE_FILL_RESULT_EXTENTS.includes(normalized) ? normalized : null;
};

const stableValue = (value) => {
  if (Array.isArray(value)) return value.map(stableValue);
  if (value && typeof value === 'object') {
    return Object.keys(value).sort().reduce((result, key) => {
      if (value[key] !== undefined) result[key] = stableValue(value[key]);
      return result;
    }, {});
  }
  return value;
};

export const stableFeatureFillStringify = (value) => {
  const normalized = stableValue(value);
  return JSON.stringify(normalized === undefined ? null : normalized);
};

const signatureHash = (value) => {
  const source = typeof value === 'string' ? value : stableFeatureFillStringify(value);
  let hash = 0x811c9dc5;
  for (let index = 0; index < source.length; index += 1) {
    hash ^= source.charCodeAt(index);
    hash = Math.imul(hash, 0x01000193);
  }
  return (hash >>> 0).toString(16).padStart(8, '0');
};

export const specificRuleIdentitySignature = (rule, ruleIndex = null) => (
  stableFeatureFillStringify({
    ...(Number.isInteger(ruleIndex) ? { ruleIndex } : {}),
    feat: text(rule?.feat),
    qual: text(rule?.qual),
    val: text(rule?.val),
    color: normalizeRuleColor(rule?.color) || text(rule?.color).toLowerCase(),
    cap: normalizeCaption(rule?.cap)
  })
);

const rawResultKey = (feature) => text(feature?.resultKey ?? feature?.result_key);
const rawRecordKey = (feature) => text(feature?.recordKey ?? feature?.record_key);
const rawBiologicalFeatureId = (feature) => text(
  feature?.biologicalFeatureId ?? feature?.biological_feature_id
);

export const stableFeatureTargetKey = (feature) => {
  const explicit = text(feature?.targetFeatureKey ?? feature?.target_feature_key);
  if (explicit) return explicit;
  const resultKey = rawResultKey(feature);
  const recordKey = rawRecordKey(feature);
  const featureId = rawBiologicalFeatureId(feature);
  return resultKey && recordKey && featureId
    ? stableFeatureFillStringify([resultKey, recordKey, featureId])
    : '';
};

export const targetSetSignature = (featureKeys) => stableFeatureFillStringify(
  [...new Set((Array.isArray(featureKeys) ? featureKeys : []).map(text).filter(Boolean))].sort()
);

const firstText = (...values) => {
  for (const value of values) {
    if (Array.isArray(value)) {
      const nested = firstText(...value);
      if (nested) return nested;
      continue;
    }
    const normalized = text(value);
    if (normalized) return normalized;
  }
  return '';
};

const rawFeatureAliases = (feature) => [...new Set([
  stableFeatureTargetKey(feature),
  feature?.featureKey,
  feature?.feature_key,
  feature?.id,
  feature?.stable_override_key,
  rawRecordKey(feature) && rawBiologicalFeatureId(feature)
    ? `${rawRecordKey(feature)}\u0000${rawBiologicalFeatureId(feature)}`
    : '',
  feature?.instanceHash
].map(text).filter(Boolean))];

const normalizedResultIndex = (feature) => {
  const value = Number(feature?.resultIndex ?? feature?.result_index);
  return Number.isSafeInteger(value) && value >= 0 ? value : null;
};

const normalizeFeature = (feature, rules, paletteColors) => {
  const featureKey = stableFeatureTargetKey(feature);
  const resultKey = rawResultKey(feature);
  if (!featureKey || !resultKey) return null;
  const resolvedRule = resolveOrderedSpecificRule(feature, rules);
  const viewModel = resolveFeatureFillViewModel({
    feature,
    specificRules: rules,
    paletteColors,
    explicitValue: feature?.explicitFillValue ?? feature?.explicit_fill_value
  });
  const annotation = getFeatureSourceAnnotation(feature);
  const effectiveCaption = firstText(
    feature?.effectiveLegendCaption,
    feature?.effective_legend_caption,
    feature?.legendCaption,
    feature?.legend_caption,
    resolvedRule?.rule?.cap,
    feature?.type
  );
  return {
    featureKey,
    aliases: rawFeatureAliases(feature),
    resultKey,
    resultName: text(feature?.resultName ?? feature?.result_name),
    resultIndex: normalizedResultIndex(feature),
    recordKey: rawRecordKey(feature),
    biologicalFeatureId: rawBiologicalFeatureId(feature),
    instanceHash: isFeatureInstanceHash(feature?.instanceHash)
      ? text(feature?.instanceHash)
      : '',
    type: text(feature?.type),
    similarityGroupIds: getFeatureSimilarityGroupIds(feature),
    legendCaption: effectiveCaption,
    renderedLabel: getFeatureRenderedLabel(feature),
    sourceAnnotationLabel: annotation.label,
    sourceAnnotationQualifier: annotation.qualifier,
    rendered: feature?.rendered !== false && feature?.hidden !== true && Boolean(firstText(
      feature?.renderedFeatureSvgId,
      feature?.rendered_feature_svg_id,
      feature?.renderedSvgId,
      feature?.rendered_svg_id,
      feature?.svgId,
      feature?.svg_id
    )),
    effectiveColor: viewModel.effectiveColor,
    explicitValue: viewModel.explicitValue,
    matchingRule: resolvedRule ? cloneJson(resolvedRule.rule) : null,
    matchingRuleIndex: resolvedRule?.ruleIndex ?? null,
    raw: cloneJson(feature)
  };
};

const mergeCatalogFeatures = (snapshot) => {
  const catalogFeatures = [];
  const similarityGroupsByTargetKey = new Map();
  if (Array.isArray(snapshot?.catalog?.items)) {
    snapshot.catalog.items.forEach((item, resultIndex) => {
      const resultKey = catalogResultKey(item);
      const resultName = text(item?.resultName);
      const renderedByPair = new Map(
        (Array.isArray(item?.features) ? item.features : []).map((feature) => [
          stableFeatureFillStringify([
            text(feature?.recordKey),
            text(feature?.biologicalFeatureId)
          ]),
          feature
        ])
      );
      (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : [])
        .forEach((feature) => {
          const pair = stableFeatureFillStringify([
            text(feature?.recordKey),
            text(feature?.biologicalFeatureId)
          ]);
          catalogFeatures.push({
            ...cloneJson(feature),
            ...cloneJson(renderedByPair.get(pair) || {}),
            resultKey,
            resultName,
            resultIndex
          });
        });
      (Array.isArray(item?.orthogroups) ? item.orthogroups : [])
        .filter((group) => (
          text(group?.presentationScope) !== 'adjacent_local'
          && text(group?.groupKind) !== 'collinear_gene_group'
        ))
        .forEach((group) => {
          const groupId = firstText(group?.id, group?.orthogroupId, group?.orthogroup_id);
          if (!groupId) return;
          (Array.isArray(group?.members) ? group.members : []).forEach((member) => {
            const targetKey = stableFeatureTargetKey({ ...member, resultKey });
            if (!targetKey) return;
            if (!similarityGroupsByTargetKey.has(targetKey)) {
              similarityGroupsByTargetKey.set(targetKey, new Set());
            }
            similarityGroupsByTargetKey.get(targetKey).add(groupId);
          });
        });
    });
  }

  const biological = Array.isArray(snapshot?.features)
    ? snapshot.features
    : (Array.isArray(snapshot?.biologicalFeatures)
      ? snapshot.biologicalFeatures
      : catalogFeatures);
  const rendered = Array.isArray(snapshot?.extractedFeatures)
    ? snapshot.extractedFeatures
    : [];
  const renderedByIdentity = new Map();
  rendered.forEach((feature) => {
    const key = stableFeatureTargetKey(feature);
    if (key && !renderedByIdentity.has(key)) renderedByIdentity.set(key, feature);
  });
  const sourceFeatures = biological.length > 0 ? biological : rendered;
  const labelsBySvgId = new Map();
  (Array.isArray(snapshot?.editableLabels) ? snapshot.editableLabels : []).forEach((entry) => {
    const svgId = text(entry?.featureId ?? entry?.feature_id);
    if (!svgId || labelsBySvgId.has(svgId)) return;
    labelsBySvgId.set(svgId, String(entry?.text ?? ''));
  });

  return sourceFeatures.map((feature) => {
    const key = stableFeatureTargetKey(feature);
    const merged = {
      ...cloneJson(feature),
      ...cloneJson(renderedByIdentity.get(key) || {})
    };
    const svgIds = [...new Set([
      merged?.renderedFeatureSvgId,
      merged?.rendered_feature_svg_id,
      merged?.renderedSvgId,
      merged?.rendered_svg_id,
      merged?.svgId,
      merged?.svg_id
    ].map(text).filter(Boolean))];
    const label = svgIds.map((svgId) => labelsBySvgId.get(svgId)).find(
      (value) => value !== undefined
    );
    if (!getFeatureRenderedLabel(merged) && label !== undefined) merged.renderedLabel = label;
    const catalogGroups = similarityGroupsByTargetKey.get(key);
    if (catalogGroups?.size) {
      merged.orthogroupIds = [...new Set([
        ...getFeatureSimilarityGroupIds(merged),
        ...catalogGroups
      ])].sort();
    }
    return merged;
  });
};

const normalizeSnapshot = (snapshot) => {
  const rules = cloneJson(
    Array.isArray(snapshot?.specificRules)
      ? snapshot.specificRules
      : (Array.isArray(snapshot?.manualSpecificRules) ? snapshot.manualSpecificRules : [])
  );
  const paletteColors = cloneJson(snapshot?.paletteColors ?? snapshot?.appliedPaletteColors ?? {});
  const byKey = new Map();
  mergeCatalogFeatures(snapshot).forEach((feature) => {
    const normalized = normalizeFeature(feature, rules, paletteColors);
    if (!normalized || byKey.has(normalized.featureKey)) return;
    byKey.set(normalized.featureKey, normalized);
  });
  const features = [...byKey.values()].sort((left, right) => (
    left.resultKey.localeCompare(right.resultKey)
    || left.featureKey.localeCompare(right.featureKey)
  ));
  return {
    features,
    rules,
    paletteColors,
    resultGenerationKey: text(snapshot?.resultGenerationKey),
    documentEpoch: text(snapshot?.documentEpoch),
    styleFingerprint: text(snapshot?.styleFingerprint),
    exactScopeAvailable: snapshot?.exactScopeAvailable !== false,
    selectedFeatureKeys: Array.isArray(snapshot?.selectedFeatureKeys)
      ? snapshot.selectedFeatureKeys.map(text).filter(Boolean)
      : []
  };
};

export const normalizeFeatureFillIntent = (intent = {}) => ({
  targetFeatureKey: text(intent?.targetFeatureKey) || text(intent?.targetFeatureKeys?.[0]),
  targetFeatureKeys: Array.isArray(intent?.targetFeatureKeys)
    ? intent.targetFeatureKeys.map(text).filter(Boolean)
    : [],
  value: normalizeFeatureFillValue(intent?.value),
  requestedCaption: normalizeCaption(intent?.requestedCaption),
  source: text(intent?.source)
    ? (VALID_SOURCES.has(text(intent.source)) ? text(intent.source) : null)
    : 'popup',
  originResultIndex: intent?.originResultIndex !== null
    && intent?.originResultIndex !== undefined
    && Number.isSafeInteger(Number(intent.originResultIndex))
    ? Number(intent.originResultIndex)
    : null,
  originResultKey: text(intent?.originResultKey),
  resultGenerationKey: text(intent?.resultGenerationKey),
  documentEpoch: text(intent?.documentEpoch),
  styleFingerprint: text(intent?.styleFingerprint),
  requestedSemanticScope: normalizeFeatureFillSemanticScope(
    intent?.semanticScope ?? intent?.requestedSemanticScope
  ),
  requestedResultExtent: normalizeFeatureFillResultExtent(intent?.resultExtent),
  token: text(intent?.token)
});

const findFeature = (features, wantedKey) => features.find(
  (feature) => feature.featureKey === wantedKey || feature.aliases.includes(wantedKey)
) || null;

const staleDiagnostic = (snapshot, intent) => {
  for (const [field, label] of [
    ['resultGenerationKey', 'Result generation changed'],
    ['documentEpoch', 'Document changed'],
    ['styleFingerprint', 'Feature style changed']
  ]) {
    if (snapshot[field] && intent[field] && snapshot[field] !== intent[field]) return label;
  }
  return '';
};

const exactRulesForFeatures = (features) => ({
  kind: 'exact-feature-rules',
  selectors: features.map((feature) => ({
    featureKey: feature.featureKey,
    feat: feature.type,
    qual: FEATURE_INSTANCE_HASH_QUALIFIER,
    val: feature.instanceHash
  }))
});

const sharedSemanticSelectorLiteral = (targets) => {
  const literals = new Set();
  for (const feature of targets) {
    if (!isFeatureSemanticScopeRule(feature.matchingRule)) return '';
    literals.add(text(feature.matchingRule.val));
  }
  return literals.size === 1 ? [...literals][0] : '';
};

const durableRuleIntent = (scope, target, targets, semanticValue = '') => {
  if (EXACT_SCOPES.has(scope)) return exactRulesForFeatures(targets);
  if (scope === 'matching-rule') {
    return {
      kind: 'update-specific-rule',
      ruleIndex: target.matchingRuleIndex,
      ruleIdentity: specificRuleIdentitySignature(
        target.matchingRule,
        target.matchingRuleIndex
      )
    };
  }
  if (scope === 'feature-type') {
    return { kind: 'feature-type', featureType: target.type };
  }
  if (scope === 'legend-caption') {
    const selector = sharedSemanticSelectorLiteral(targets);
    if (selector) return { kind: 'semantic-selector', selector };
    return { kind: 'legend-caption', caption: target.legendCaption };
  }
  if (scope === 'rendered-label') {
    return { kind: 'rendered-label', label: target.renderedLabel };
  }
  if (scope === 'similarity-group') {
    return {
      kind: 'similarity-group',
      groupId: text(semanticValue) || target.similarityGroupIds[0]
    };
  }
  return {
    kind: 'source-annotation-label',
    qualifier: target.sourceAnnotationQualifier,
    label: target.sourceAnnotationLabel
  };
};

const targetsByResult = (targets) => {
  const groups = new Map();
  targets.forEach((feature) => {
    if (!groups.has(feature.resultKey)) {
      groups.set(feature.resultKey, {
        resultKey: feature.resultKey,
        resultName: feature.resultName,
        resultIndex: feature.resultIndex,
        featureKeys: []
      });
    }
    groups.get(feature.resultKey).featureKeys.push(feature.featureKey);
  });
  return [...groups.values()]
    .map((group) => ({
      resultKey: group.resultKey,
      resultName: group.resultName,
      resultIndex: group.resultIndex,
      featureKeys: [...new Set(group.featureKeys)].sort()
    }))
    .sort((left, right) => (
      (left.resultIndex ?? Number.MAX_SAFE_INTEGER) - (right.resultIndex ?? Number.MAX_SAFE_INTEGER)
      || left.resultKey.localeCompare(right.resultKey)
    ));
};

const candidateLabel = (scope, targetCount, resultCount, target, semanticValue = '') => {
  if (scope === 'single') {
    return `Update only this feature in ${target.resultName || 'this output'}`;
  }
  const scopeLabel = {
    'matching-rule': 'Update matching rule',
    'feature-type': `Update all ${target.type} features`,
    'legend-caption': `Update legend caption “${target.legendCaption}”`,
    'rendered-label': `Update rendered label “${target.renderedLabel}”`,
    'source-annotation-label': `Update source label “${target.sourceAnnotationLabel}”`,
    'similarity-group': `Update similarity group “${text(semanticValue) || target.similarityGroupIds[0]}”`,
    'selected-features': 'Update selected features'
  }[scope] || 'Update features';
  return `${scopeLabel} — ${targetCount} features in ${resultCount} outputs`;
};

const buildCandidate = (scope, target, rawTargets, semanticValue = '') => {
  const targets = [...rawTargets].sort((left, right) => (
    left.resultKey.localeCompare(right.resultKey)
    || left.featureKey.localeCompare(right.featureKey)
  ));
  const groupedTargets = targetsByResult(targets);
  const resultExtent = GROUP_SCOPES.has(scope) || groupedTargets.length > 1
    ? 'all-affected-results'
    : 'current-result';
  const signature = targetSetSignature(targets.map((feature) => feature.featureKey));
  const durableValue = text(semanticValue);
  return {
    id: `${scope}:${signatureHash([durableValue, signature])}`,
    semanticScope: scope,
    resultExtent,
    label: candidateLabel(
      scope,
      targets.length,
      groupedTargets.length,
      target,
      durableValue
    ),
    targetCount: targets.length,
    affectedResultCount: groupedTargets.length,
    targetsByResult: groupedTargets,
    targetFeatureKeys: targets.map((feature) => feature.featureKey),
    targetSetSignature: signature,
    durableRuleIntent: durableRuleIntent(scope, target, targets, durableValue)
  };
};

const featuresMatchingRule = (features, target) => {
  if (!target.matchingRule) return [];
  return features.filter((feature) => specificRuleMatchesFeature(feature.raw, target.matchingRule));
};

const selectedFeatures = (snapshot, intent) => {
  const selectedKeys = intent.targetFeatureKeys.length > 0
    ? intent.targetFeatureKeys
    : snapshot.selectedFeatureKeys;
  const selected = [];
  selectedKeys.forEach((key) => {
    const feature = findFeature(snapshot.features, key);
    if (feature && !selected.some((candidate) => candidate.featureKey === feature.featureKey)) {
      selected.push(feature);
    }
  });
  return selected;
};

const captionConflict = (snapshot, intent, target) => {
  const caption = intent.requestedCaption;
  const requestedColor = intent.value?.kind === 'color'
    ? intent.value.color
    : (intent.value?.kind === 'none' ? 'none' : target.effectiveColor);
  if (!caption || !requestedColor) return null;
  const colors = new Set();
  const captionFeatureKeys = [];
  snapshot.features.forEach((feature) => {
    if (captionKey(feature.legendCaption) === captionKey(caption)) {
      colors.add(feature.effectiveColor);
      captionFeatureKeys.push(feature.featureKey);
    }
  });
  snapshot.rules.forEach((rule) => {
    if (captionKey(rule?.cap) === captionKey(caption)) {
      const color = normalizeRuleColor(rule?.color);
      if (color) colors.add(color);
    }
  });
  const existingColor = [...colors].find((color) => color !== requestedColor);
  if (!existingColor) return null;
  const usedCaptions = new Set(snapshot.features.map((feature) => captionKey(feature.legendCaption)));
  let suffix = 2;
  let separateCaption = `${caption} (${suffix})`;
  while (usedCaptions.has(captionKey(separateCaption))) {
    suffix += 1;
    separateCaption = `${caption} (${suffix})`;
  }
  return {
    kind: 'caption-color',
    caption,
    requestedColor,
    existingColor,
    captionFeatureKeys: [...new Set(captionFeatureKeys)].sort(),
    choices: [{
      id: 'reuse-existing',
      label: `Use existing ${existingColor}`,
      value: { kind: existingColor === 'none' ? 'none' : 'color', ...(existingColor === 'none' ? {} : { color: existingColor }) },
      caption
    }, {
      id: 'separate-caption',
      label: `Create “${separateCaption}”`,
      value: cloneJson(intent.value),
      caption: separateCaption
    }]
  };
};

const invalidPlan = (intent, diagnostics) => ({
  token: intent.token || `fill-plan:${signatureHash({ intent, diagnostics })}`,
  status: 'invalid',
  intent,
  candidates: [],
  diagnostics
});

export const planFeatureFillChange = (rawSnapshot = {}, rawIntent = {}) => {
  const snapshot = normalizeSnapshot(rawSnapshot);
  let intent = normalizeFeatureFillIntent(rawIntent);
  const rawRequestedScope = text(rawIntent?.semanticScope ?? rawIntent?.requestedSemanticScope);
  const rawRequestedExtent = text(rawIntent?.resultExtent);
  if (text(rawIntent?.source) && !intent.source) {
    return invalidPlan(intent, [{
      code: 'unsupported-source',
      message: `Feature fill source ${text(rawIntent.source)} is unavailable`
    }]);
  }
  if (rawRequestedScope && !intent.requestedSemanticScope) {
    return invalidPlan(intent, [{
      code: 'unsupported-scope',
      message: `Scope ${rawRequestedScope} is unavailable`
    }]);
  }
  if (rawRequestedExtent && !intent.requestedResultExtent) {
    return invalidPlan(intent, [{
      code: 'unsupported-result-extent',
      message: `Result extent ${rawRequestedExtent} is unavailable`
    }]);
  }
  const stale = staleDiagnostic(snapshot, intent);
  if (stale) return invalidPlan(intent, [{ code: 'stale-intent', message: stale }]);
  if (!intent.value) {
    return invalidPlan(intent, [{ code: 'invalid-value', message: 'Invalid Feature fill value' }]);
  }
  const target = findFeature(snapshot.features, intent.targetFeatureKey);
  if (!target) {
    return invalidPlan(intent, [{ code: 'missing-target', message: 'Feature target is unavailable' }]);
  }
  if (!intent.requestedCaption && target.legendCaption) {
    intent = { ...intent, requestedCaption: target.legendCaption };
  }
  if (
    intent.originResultKey && intent.originResultKey !== target.resultKey
    || intent.originResultIndex !== null && intent.originResultIndex !== target.resultIndex
  ) {
    return invalidPlan(intent, [{ code: 'stale-result', message: 'Origin Result changed' }]);
  }
  if (
    intent.requestedResultExtent === 'current-result'
    && GROUP_SCOPES.has(intent.requestedSemanticScope)
  ) {
    return invalidPlan(intent, [{
      code: 'unsupported-current-result-snapshot',
      message: 'Broad current-Result snapshots are not supported'
    }]);
  }

  const candidates = [];
  const addGroup = (scope, targets, semanticValue = '') => {
    if (targets.length > 1) {
      candidates.push(buildCandidate(scope, target, targets, semanticValue));
    }
  };
  addGroup('matching-rule', featuresMatchingRule(snapshot.features, target));
  if (target.type) {
    addGroup('feature-type', snapshot.features.filter(
      (feature) => feature.type === target.type
    ));
  }
  if (target.legendCaption) {
    addGroup('legend-caption', snapshot.features.filter(
      (feature) => captionKey(feature.legendCaption) === captionKey(target.legendCaption)
    ));
  }
  if (target.renderedLabel) {
    addGroup('rendered-label', snapshot.features.filter(
      (feature) => feature.renderedLabel === target.renderedLabel
    ));
  }
  if (target.sourceAnnotationLabel) {
    addGroup('source-annotation-label', snapshot.features.filter((feature) => (
      feature.sourceAnnotationLabel === target.sourceAnnotationLabel
    )));
  }
  target.similarityGroupIds.forEach((groupId) => {
    addGroup('similarity-group', snapshot.features.filter(
      (feature) => feature.similarityGroupIds.includes(groupId)
    ), groupId);
  });

  const selected = selectedFeatures(snapshot, intent);
  const selectedHasExactIdentity = selected.length > 0
    && selected.every((feature) => feature.instanceHash);
  if (selected.length > 1 && snapshot.exactScopeAvailable && selectedHasExactIdentity) {
    candidates.push(buildCandidate('selected-features', target, selected));
  }
  const exactAvailable = snapshot.exactScopeAvailable && Boolean(target.instanceHash);
  if (exactAvailable) candidates.push(buildCandidate('single', target, [target]));

  candidates.sort((left, right) => (
    FEATURE_FILL_SCOPE_ORDER.indexOf(left.semanticScope)
    - FEATURE_FILL_SCOPE_ORDER.indexOf(right.semanticScope)
  ));
  const diagnostics = [];
  if (!exactAvailable || (selected.length > 1 && !selectedHasExactIdentity)) {
    diagnostics.push({ code: 'exact-scope-unavailable', message: GENERATE_EXACT_DIAGNOSTIC });
  }
  if (candidates.length === 0) return invalidPlan(intent, diagnostics.length ? diagnostics : [{
    code: 'no-targets', message: 'No durable Feature fill target is available'
  }]);

  const groupCandidates = candidates.filter((candidate) => GROUP_SCOPES.has(candidate.semanticScope));
  const requestedCandidate = intent.requestedSemanticScope
    ? candidates.find((candidate) => candidate.semanticScope === intent.requestedSemanticScope)
    : null;
  if (intent.requestedSemanticScope && !requestedCandidate) {
    return invalidPlan(intent, [...diagnostics, {
      code: 'unsupported-scope',
      message: `Scope ${intent.requestedSemanticScope} is unavailable`
    }]);
  }
  if (
    requestedCandidate
    && intent.requestedResultExtent
    && intent.requestedResultExtent !== requestedCandidate.resultExtent
  ) {
    return invalidPlan(intent, [...diagnostics, {
      code: 'unsupported-result-extent',
      message: `Scope ${intent.requestedSemanticScope} requires ${requestedCandidate.resultExtent}`
    }]);
  }

  const token = intent.token || `fill-plan:${signatureHash({
    intent,
    target: target.featureKey,
    candidates: candidates.map((candidate) => candidate.id),
    resultGenerationKey: snapshot.resultGenerationKey,
    documentEpoch: snapshot.documentEpoch,
    styleFingerprint: snapshot.styleFingerprint
  })}`;
  const conflict = captionConflict(snapshot, intent, target);
  const selectedCandidate = requestedCandidate
    || (intent.source === 'selection-toolbar'
      ? candidates.find((candidate) => candidate.semanticScope === 'selected-features')
      : null)
    || (groupCandidates.length === 0 ? candidates.find((candidate) => candidate.semanticScope === 'single') : null);
  const targetNoop = selectedCandidate?.semanticScope === 'single' && (
    intent.value.kind === 'inherit' && target.explicitValue === null
  ) || (
    intent.value.kind === 'none' && target.explicitValue === 'none'
  ) || (
    intent.value.kind === 'color'
    && target.explicitValue === intent.value.color
    && (!intent.requestedCaption || captionKey(intent.requestedCaption) === captionKey(target.legendCaption))
  );

  return {
    token,
    status: targetNoop ? 'noop' : (conflict ? 'conflict' : (selectedCandidate ? 'ready' : 'needs-scope')),
    intent: { ...intent, targetFeatureKey: target.featureKey },
    candidates,
    diagnostics,
    ...(conflict ? { conflict } : {}),
    ...(selectedCandidate ? { selectedCandidateId: selectedCandidate.id } : {})
  };
};

const resolveConflict = (plan, conflictChoiceId) => {
  if (!plan.conflict) return { intent: cloneJson(plan.intent), conflictChoice: null };
  const choice = plan.conflict.choices.find((candidate) => candidate.id === conflictChoiceId);
  if (!choice) return null;
  return {
    intent: {
      ...cloneJson(plan.intent),
      value: cloneJson(choice.value),
      requestedCaption: choice.caption
    },
    conflictChoice: choice.id
  };
};

export const resolveFeatureFillPlan = (
  rawPlan,
  selectedCandidateId = rawPlan?.selectedCandidateId,
  { conflictChoiceId = null } = {}
) => {
  const plan = cloneJson(rawPlan);
  if (!plan || !['ready', 'needs-scope', 'conflict'].includes(plan.status)) return plan;
  const candidate = plan.candidates.find((entry) => entry.id === selectedCandidateId);
  if (!candidate) {
    return {
      ...plan,
      status: 'invalid',
      diagnostics: [...plan.diagnostics, {
        code: 'missing-scope',
        message: 'The selected Feature fill scope is unavailable'
      }]
    };
  }
  const candidateTargets = new Set(candidate.targetFeatureKeys || []);
  const requiresCaptionConflictChoice = Boolean(
    plan.conflict
    && (plan.conflict.captionFeatureKeys || []).some((key) => !candidateTargets.has(key))
  );
  const conflictResolution = requiresCaptionConflictChoice
    ? resolveConflict(plan, conflictChoiceId)
    : { intent: cloneJson(plan.intent), conflictChoice: null };
  if (!conflictResolution) {
    return { ...plan, selectedCandidateId: candidate.id };
  }
  const affectedResultKeys = candidate.targetsByResult.map((entry) => entry.resultKey);
  return {
    ...plan,
    status: 'ready',
    intent: conflictResolution.intent,
    selectedCandidateId: candidate.id,
    semanticScope: candidate.semanticScope,
    resultExtent: candidate.resultExtent,
    affectedResultKeys,
    targetsByResult: cloneJson(candidate.targetsByResult),
    semanticBefore: { styleFingerprint: text(plan.intent.styleFingerprint) },
    semanticAfter: {
      value: cloneJson(conflictResolution.intent.value),
      requestedCaption: conflictResolution.intent.requestedCaption,
      durableRuleIntent: cloneJson(candidate.durableRuleIntent)
    },
    legendProjectionByResult: candidate.targetsByResult.map((entry) => ({
      resultKey: entry.resultKey,
      featureKeys: cloneJson(entry.featureKeys)
    })),
    ...(conflictResolution.conflictChoice
      ? { conflictChoiceId: conflictResolution.conflictChoice }
      : {})
  };
};
