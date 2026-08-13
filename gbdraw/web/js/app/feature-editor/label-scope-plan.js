import { biologicalFeatureKey, catalogResultKey } from '../../services/feature-catalog.js';
import { isFeatureInstanceHash } from '../../services/feature-instance-identity.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
import {
  stableFeatureFillStringify,
  stableFeatureTargetKey,
  targetSetSignature
} from './fill-scope-plan.js';
import {
  normalizeFeatureLabelText,
  resolveFeatureLabelViewModel
} from './label-view-model.js';

export const FEATURE_LABEL_SCOPE_ORDER = Object.freeze([
  'rendered-label',
  'source-annotation-label',
  'selected-features',
  'single'
]);

const VALID_SOURCES = new Set(['popup', 'feature-list', 'selection-toolbar']);
const text = (value) => String(value ?? '').trim();

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

const deepFreeze = (value) => {
  if (!value || typeof value !== 'object' || Object.isFrozen(value)) return value;
  Object.values(value).forEach(deepFreeze);
  return Object.freeze(value);
};

const immutable = (value) => deepFreeze(cloneJson(value));

const signatureHash = (value) => {
  const source = typeof value === 'string' ? value : JSON.stringify(stableValue(value));
  let hash = 0x811c9dc5;
  for (let index = 0; index < source.length; index += 1) {
    hash ^= source.charCodeAt(index);
    hash = Math.imul(hash, 0x01000193);
  }
  return (hash >>> 0).toString(16).padStart(8, '0');
};

export const featureLabelStateFingerprint = ({
  labelTextFeatureOverrides = {},
  labelTextBulkOverrides = {},
  labelTextFeatureOverrideSources = {}
} = {}) => `ls1_${signatureHash(stableValue({
  labelTextFeatureOverrides,
  labelTextBulkOverrides,
  labelTextFeatureOverrideSources
}))}`;

export const stableFeatureLabelStringify = (value) => JSON.stringify(stableValue(value));

export const normalizeFeatureLabelSemanticScope = (value) => {
  const normalized = text(value);
  return FEATURE_LABEL_SCOPE_ORDER.includes(normalized) ? normalized : null;
};

const resultIndex = (feature, fallback = null) => {
  const normalized = Number(feature?.resultIndex ?? feature?.result_index ?? fallback);
  return Number.isSafeInteger(normalized) && normalized >= 0 ? normalized : null;
};

const presentationIndex = (snapshot) => {
  const bySvgId = new Map();
  const source = snapshot?.presentationLabelsBySvgId;
  Object.entries(source && typeof source === 'object' ? source : {}).forEach(([svgId, entry]) => {
    if (!text(svgId)) return;
    bySvgId.set(text(svgId), {
      text: String(entry?.text ?? ''),
      sourceText: String(entry?.sourceText ?? entry?.source_text ?? ''),
      labelKey: text(entry?.labelKey ?? entry?.label_key)
    });
  });
  (Array.isArray(snapshot?.editableLabels) ? snapshot.editableLabels : []).forEach((entry) => {
    const svgId = text(entry?.featureId ?? entry?.feature_id);
    if (!svgId) return;
    bySvgId.set(svgId, {
      text: String(entry?.text ?? ''),
      sourceText: String(entry?.sourceText ?? entry?.source_text ?? ''),
      labelKey: text(entry?.key ?? entry?.labelKey)
    });
  });
  return bySvgId;
};

const overlayIndex = (snapshot) => {
  const byPair = new Map();
  const features = [
    ...(Array.isArray(snapshot?.features) ? snapshot.features : []),
    ...(Array.isArray(snapshot?.biologicalFeatures) ? snapshot.biologicalFeatures : []),
    ...(Array.isArray(snapshot?.extractedFeatures) ? snapshot.extractedFeatures : [])
  ];
  features.forEach((feature) => {
    const resultKey = text(feature?.resultKey ?? feature?.result_key);
    const pairKey = biologicalFeatureKey(
      feature?.recordKey ?? feature?.record_key,
      feature?.biologicalFeatureId ?? feature?.biological_feature_id
    );
    if (!pairKey) return;
    const key = stableFeatureFillStringify([resultKey, pairKey]);
    byPair.set(key, { ...(byPair.get(key) || {}), ...cloneJson(feature) });
  });
  return byPair;
};

const catalogFeatures = (snapshot) => {
  if (!Array.isArray(snapshot?.catalog?.items)) {
    return (Array.isArray(snapshot?.features) ? snapshot.features : []).map((feature) => ({
      feature,
      rendered: [feature],
      resultKey: text(feature?.resultKey ?? feature?.result_key),
      resultName: text(feature?.resultName ?? feature?.result_name),
      resultIndex: resultIndex(feature)
    }));
  }
  const overlays = overlayIndex(snapshot);
  return snapshot.catalog.items.flatMap((item, itemIndex) => {
    const resultKey = catalogResultKey(item);
    const renderedByPair = new Map();
    (Array.isArray(item?.features) ? item.features : []).forEach((rendered) => {
      const pairKey = biologicalFeatureKey(rendered?.recordKey, rendered?.biologicalFeatureId);
      if (!renderedByPair.has(pairKey)) renderedByPair.set(pairKey, []);
      renderedByPair.get(pairKey).push(rendered);
    });
    return (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).map((feature) => {
      const pairKey = biologicalFeatureKey(feature?.recordKey, feature?.biologicalFeatureId);
      const overlay = overlays.get(stableFeatureFillStringify([resultKey, pairKey])) || {};
      return {
        feature: { ...cloneJson(feature), ...cloneJson(overlay) },
        rendered: renderedByPair.get(pairKey) || [],
        resultKey,
        resultName: text(item?.resultName),
        resultIndex: itemIndex
      };
    });
  });
};

const normalizeFeature = ({
  feature,
  rendered,
  resultKey,
  resultName,
  resultIndex: fallbackResultIndex,
  presentation,
  featureOverrides,
  bulkOverrides,
  featureOverrideSources
}) => {
  const recordKey = text(feature?.recordKey ?? feature?.record_key);
  const biologicalFeatureId = text(
    feature?.biologicalFeatureId ?? feature?.biological_feature_id
  );
  const renderedSvgIds = [...new Set([
    ...(Array.isArray(rendered)
      ? rendered.map((entry) => text(entry?.svgId ?? entry?.svg_id))
      : []),
    text(feature?.renderedFeatureSvgId ?? feature?.rendered_feature_svg_id),
    text(feature?.svgId ?? feature?.svg_id)
  ].filter(Boolean))];
  const featureKey = stableFeatureTargetKey({
    ...feature,
    resultKey,
    recordKey,
    biologicalFeatureId
  });
  if (!resultKey || !recordKey || !biologicalFeatureId || !featureKey || renderedSvgIds.length === 0) {
    return null;
  }
  const presented = renderedSvgIds.map((svgId) => presentation.get(svgId)).find(Boolean) || null;
  const viewModel = resolveFeatureLabelViewModel({
    feature,
    renderedSvgIds,
    featureOverrides,
    bulkOverrides,
    featureOverrideSources,
    presentationText: presented ? presented.text : null,
    presentationSourceText: presented?.sourceText || ''
  });
  return {
    featureKey,
    aliases: [...new Set([
      featureKey,
      biologicalFeatureKey(recordKey, biologicalFeatureId),
      text(feature?.instanceHash),
      text(feature?.targetFeatureKey ?? feature?.target_feature_key),
      ...renderedSvgIds
    ].filter(Boolean))],
    resultKey,
    resultName,
    resultIndex: resultIndex(feature, fallbackResultIndex),
    recordKey,
    biologicalFeatureId,
    instanceHash: isFeatureInstanceHash(feature?.instanceHash)
      ? text(feature.instanceHash)
      : '',
    renderedSvgIds,
    labelKey: presented?.labelKey || '',
    renderedLabel: viewModel.effectiveText,
    sourceAnnotationLabel: viewModel.sourceText,
    sourceAnnotationQualifier: viewModel.sourceQualifier,
    viewModel,
    raw: cloneJson(feature)
  };
};

const normalizeSnapshot = (snapshot) => {
  const featureOverrides = cloneJson(snapshot?.labelTextFeatureOverrides || {});
  const bulkOverrides = cloneJson(snapshot?.labelTextBulkOverrides || {});
  const featureOverrideSources = cloneJson(snapshot?.labelTextFeatureOverrideSources || {});
  const presentation = presentationIndex(snapshot);
  const byKey = new Map();
  catalogFeatures(snapshot).forEach((entry) => {
    const feature = normalizeFeature({
      ...entry,
      presentation,
      featureOverrides,
      bulkOverrides,
      featureOverrideSources
    });
    if (feature && !byKey.has(feature.featureKey)) byKey.set(feature.featureKey, feature);
  });
  return {
    features: [...byKey.values()].sort((left, right) => (
      left.resultKey.localeCompare(right.resultKey)
      || left.featureKey.localeCompare(right.featureKey)
    )),
    featureOverrides,
    bulkOverrides,
    featureOverrideSources,
    resultGenerationKey: Number(snapshot?.resultGenerationKey ?? 0),
    documentEpoch: Number(snapshot?.documentEpoch ?? 0),
    styleRevision: Number(snapshot?.styleRevision ?? snapshot?.semanticStyleRevision ?? 0),
    styleFingerprint: text(snapshot?.styleFingerprint),
    exactScopeAvailable: snapshot?.exactScopeAvailable !== false,
    selectedFeatureKeys: Array.isArray(snapshot?.selectedFeatureKeys)
      ? snapshot.selectedFeatureKeys.map(text).filter(Boolean)
      : []
  };
};

export const normalizeFeatureLabelIntent = (intent = {}) => ({
  targetFeatureKey: text(intent?.targetFeatureKey) || text(intent?.targetFeatureKeys?.[0]),
  targetFeatureKeys: Array.isArray(intent?.targetFeatureKeys)
    ? intent.targetFeatureKeys.map(text).filter(Boolean)
    : [],
  newText: normalizeFeatureLabelText(intent?.newText),
  source: text(intent?.source)
    ? (VALID_SOURCES.has(text(intent.source)) ? text(intent.source) : null)
    : 'popup',
  originResultKey: text(intent?.originResultKey),
  originResultIndex: resultIndex(intent, intent?.originResultIndex),
  resultGenerationKey: Number(intent?.resultGenerationKey ?? 0),
  documentEpoch: Number(intent?.documentEpoch ?? 0),
  styleRevision: Number(intent?.styleRevision ?? intent?.semanticStyleRevision ?? 0),
  styleFingerprint: text(intent?.styleFingerprint),
  requestedSemanticScope: normalizeFeatureLabelSemanticScope(
    intent?.semanticScope ?? intent?.requestedSemanticScope
  ),
  token: text(intent?.token)
});

const findFeature = (features, wanted) => features.find(
  (feature) => feature.featureKey === wanted || feature.aliases.includes(wanted)
) || null;

const selectedFeatures = (snapshot, intent) => {
  const keys = intent.targetFeatureKeys.length > 0
    ? intent.targetFeatureKeys
    : snapshot.selectedFeatureKeys;
  const selected = [];
  keys.forEach((key) => {
    const feature = findFeature(snapshot.features, key);
    if (feature && !selected.some((entry) => entry.featureKey === feature.featureKey)) {
      selected.push(feature);
    }
  });
  return selected;
};

const targetsByResult = (features) => {
  const groups = new Map();
  features.forEach((feature) => {
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
  return [...groups.values()].map((group) => ({
    ...group,
    featureKeys: [...new Set(group.featureKeys)].sort()
  })).sort((left, right) => (
    (left.resultIndex ?? Number.MAX_SAFE_INTEGER)
      - (right.resultIndex ?? Number.MAX_SAFE_INTEGER)
    || left.resultKey.localeCompare(right.resultKey)
  ));
};

const candidateLabel = (scope, target, count, resultCount) => {
  if (scope === 'single') {
    return `Update only this label in ${target.resultName || 'this output'}`;
  }
  const prefix = {
    'rendered-label': `Update rendered label “${target.renderedLabel}”`,
    'source-annotation-label': `Update source label “${target.sourceAnnotationLabel}”`,
    'selected-features': 'Update selected feature labels'
  }[scope];
  return `${prefix} — ${count} features in ${resultCount} outputs`;
};

const selector = (feature) => ({
  featureKey: feature.featureKey,
  instanceHash: feature.instanceHash,
  renderedSvgIds: feature.renderedSvgIds,
  sourceText: feature.sourceAnnotationLabel,
  fromText: feature.renderedLabel
});

const buildCandidate = (scope, target, rawTargets) => {
  const targets = [...rawTargets].sort((left, right) => (
    left.resultKey.localeCompare(right.resultKey)
    || left.featureKey.localeCompare(right.featureKey)
  ));
  const grouped = targetsByResult(targets);
  const signature = targetSetSignature(targets.map((feature) => feature.featureKey));
  return {
    id: `${scope}:${signatureHash({ signature, value: (
      scope === 'rendered-label' ? target.renderedLabel : target.sourceAnnotationLabel
    ) })}`,
    semanticScope: scope,
    resultExtent: scope === 'single' ? 'current-result' : 'all-affected-results',
    label: candidateLabel(scope, target, targets.length, grouped.length),
    targetCount: targets.length,
    affectedResultCount: grouped.length,
    targetsByResult: grouped,
    targetFeatureKeys: targets.map((feature) => feature.featureKey),
    targetSetSignature: signature,
    durableLabelIntent: scope === 'source-annotation-label'
      ? {
          kind: 'source-label-group',
          sourceText: target.sourceAnnotationLabel,
          qualifier: target.sourceAnnotationQualifier,
          selectors: targets.map(selector)
        }
      : {
          kind: 'exact-feature-labels',
          matchKind: scope,
          matchText: scope === 'rendered-label' ? target.renderedLabel : '',
          selectors: targets.map(selector)
        }
  };
};

const invalidPlan = (intent, diagnostics) => immutable({
  token: intent.token || `label-plan:${signatureHash({ intent, diagnostics })}`,
  status: 'invalid',
  intent,
  candidates: [],
  diagnostics
});

export const planFeatureLabelChange = (rawSnapshot = {}, rawIntent = {}) => {
  const snapshot = normalizeSnapshot(rawSnapshot);
  const intent = normalizeFeatureLabelIntent(rawIntent);
  const rawScope = text(rawIntent?.semanticScope ?? rawIntent?.requestedSemanticScope);
  if (text(rawIntent?.source) && !intent.source) {
    return invalidPlan(intent, [{ code: 'unsupported-source', message: 'Feature label source is unavailable' }]);
  }
  if (rawScope && !intent.requestedSemanticScope) {
    return invalidPlan(intent, [{ code: 'unsupported-scope', message: `Label scope ${rawScope} is unavailable` }]);
  }
  if (intent.newText === null) {
    return invalidPlan(intent, [{ code: 'invalid-value', message: 'Feature label text must be a string' }]);
  }
  for (const [field, code, message] of [
    ['resultGenerationKey', 'stale-generation', 'Result generation changed'],
    ['documentEpoch', 'stale-document', 'Document changed'],
    ['styleRevision', 'stale-revision', 'Feature style revision changed'],
    ['styleFingerprint', 'stale-style', 'Feature style changed']
  ]) {
    if (snapshot[field] !== intent[field] && (snapshot[field] || intent[field])) {
      return invalidPlan(intent, [{ code, message }]);
    }
  }
  const target = findFeature(snapshot.features, intent.targetFeatureKey);
  if (!target) {
    return invalidPlan(intent, [{ code: 'missing-target', message: 'Feature label target is unavailable' }]);
  }
  if (
    (intent.originResultKey && intent.originResultKey !== target.resultKey)
    || (intent.originResultIndex !== null && intent.originResultIndex !== target.resultIndex)
  ) {
    return invalidPlan(intent, [{ code: 'stale-result', message: 'Origin Result changed' }]);
  }

  const candidates = [];
  if (target.renderedLabel) {
    const targets = snapshot.features.filter(
      (feature) => feature.renderedLabel === target.renderedLabel
    );
    if (targets.length > 1) candidates.push(buildCandidate('rendered-label', target, targets));
  }
  if (target.sourceAnnotationLabel) {
    const targets = snapshot.features.filter(
      (feature) => feature.sourceAnnotationLabel === target.sourceAnnotationLabel
    );
    if (targets.length > 1) {
      candidates.push(buildCandidate('source-annotation-label', target, targets));
    }
  }
  const selected = selectedFeatures(snapshot, intent);
  const selectedExact = selected.length > 1
    && selected.every((feature) => Boolean(feature.instanceHash));
  if (snapshot.exactScopeAvailable && selectedExact) {
    candidates.push(buildCandidate('selected-features', target, selected));
  }
  const exactAvailable = snapshot.exactScopeAvailable && Boolean(target.instanceHash);
  if (exactAvailable) candidates.push(buildCandidate('single', target, [target]));
  candidates.sort((left, right) => (
    FEATURE_LABEL_SCOPE_ORDER.indexOf(left.semanticScope)
      - FEATURE_LABEL_SCOPE_ORDER.indexOf(right.semanticScope)
  ));

  const diagnostics = [];
  if (!exactAvailable || (selected.length > 1 && !selectedExact)) {
    diagnostics.push({
      code: 'exact-scope-unavailable',
      message: 'Generate to enable exact feature edits'
    });
  }
  if (candidates.length === 0) {
    return invalidPlan(intent, diagnostics.length ? diagnostics : [{
      code: 'no-targets', message: 'No durable Feature label target is available'
    }]);
  }
  const requested = intent.requestedSemanticScope
    ? candidates.find((candidate) => candidate.semanticScope === intent.requestedSemanticScope)
    : null;
  if (intent.requestedSemanticScope && !requested) {
    return invalidPlan(intent, [...diagnostics, {
      code: 'unsupported-scope', message: `Label scope ${intent.requestedSemanticScope} is unavailable`
    }]);
  }
  const selectedCandidate = requested
    || (intent.source === 'selection-toolbar'
      ? candidates.find((candidate) => candidate.semanticScope === 'selected-features')
      : null)
    || (candidates.every((candidate) => candidate.semanticScope === 'single')
      ? candidates[0]
      : null);
  const token = intent.token || `label-plan:${signatureHash({
    intent,
    target: target.featureKey,
    candidates: candidates.map((candidate) => candidate.id)
  })}`;
  return immutable({
    token,
    status: selectedCandidate ? 'ready' : 'needs-scope',
    intent: { ...intent, targetFeatureKey: target.featureKey },
    candidates,
    diagnostics,
    labelStateFingerprint: featureLabelStateFingerprint({
      labelTextFeatureOverrides: snapshot.featureOverrides,
      labelTextBulkOverrides: snapshot.bulkOverrides,
      labelTextFeatureOverrideSources: snapshot.featureOverrideSources
    }),
    ...(selectedCandidate ? { selectedCandidateId: selectedCandidate.id } : {})
  });
};

export const resolveFeatureLabelPlan = (
  rawPlan,
  selectedCandidateId = rawPlan?.selectedCandidateId
) => {
  const plan = cloneJson(rawPlan);
  if (!plan || !['ready', 'needs-scope'].includes(plan.status)) return plan;
  const candidate = plan.candidates.find((entry) => entry.id === selectedCandidateId);
  if (!candidate) {
    return immutable({
      ...plan,
      status: 'invalid',
      diagnostics: [...(plan.diagnostics || []), {
        code: 'missing-scope', message: 'The selected Feature label scope is unavailable'
      }]
    });
  }
  return immutable({
    ...plan,
    status: 'ready',
    selectedCandidateId: candidate.id,
    semanticScope: candidate.semanticScope,
    resultExtent: candidate.resultExtent,
    affectedResultKeys: candidate.targetsByResult.map((entry) => entry.resultKey),
    targetsByResult: candidate.targetsByResult,
    semanticBefore: {
      styleFingerprint: plan.intent.styleFingerprint,
      styleRevision: plan.intent.styleRevision,
      labelStateFingerprint: plan.labelStateFingerprint
    },
    semanticAfter: {
      newText: plan.intent.newText,
      durableLabelIntent: candidate.durableLabelIntent
    }
  });
};

export { stableFeatureTargetKey };
