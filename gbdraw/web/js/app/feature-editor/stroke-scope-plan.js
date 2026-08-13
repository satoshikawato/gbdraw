import {
  biologicalFeatureKey,
  catalogResultKey
} from '../../services/feature-catalog.js';
import { isFeatureInstanceHash } from '../../services/feature-instance-identity.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
import {
  normalizeFeatureStrokeValue,
  resolveFeatureStrokeViewModel
} from './stroke-view-model.js';
import {
  stableFeatureFillStringify,
  stableFeatureTargetKey,
  targetSetSignature
} from './fill-scope-plan.js';

export const FEATURE_STROKE_SCOPE_ORDER = Object.freeze([
  'legend-caption',
  'selected-features',
  'single'
]);

const EXACT_SCOPES = new Set(['selected-features', 'single']);
const VALID_SOURCES = new Set(['popup', 'feature-list', 'selection-toolbar', 'legend-editor']);
const EXACT_SCOPE_DIAGNOSTIC = 'Generate to enable exact feature edits';
const text = (value) => String(value ?? '').trim();
const captionKey = (value) => text(value).toLowerCase();

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

export const featureStrokeStateFingerprint = ({
  featureStrokeOverrides = {},
  legendStrokeOverrides = {},
  originalSvgStroke = {}
} = {}) => `ss1_${signatureHash(stableValue({
  featureStrokeOverrides,
  legendStrokeOverrides,
  originalSvgStroke
}))}`;

export const stableFeatureStrokeStringify = (value) => JSON.stringify(stableValue(value));

export const normalizeFeatureStrokeSemanticScope = (value) => {
  const normalized = text(value);
  return FEATURE_STROKE_SCOPE_ORDER.includes(normalized) ? normalized : null;
};

const resultIndex = (feature, fallback = null) => {
  const normalized = Number(feature?.resultIndex ?? feature?.result_index ?? fallback);
  return Number.isSafeInteger(normalized) && normalized >= 0 ? normalized : null;
};

const renderedId = (feature) => text(
  feature?.renderedFeatureSvgId
  ?? feature?.rendered_feature_svg_id
  ?? feature?.renderedSvgId
  ?? feature?.rendered_svg_id
  ?? feature?.svgId
  ?? feature?.svg_id
);

const legendCaption = (feature) => text(
  feature?.effectiveLegendCaption
  ?? feature?.effective_legend_caption
  ?? feature?.legendCaption
  ?? feature?.legend_caption
  ?? feature?.caption
  ?? feature?.type
);

const normalizeFeature = ({
  feature,
  resultKey,
  resultName,
  fallbackResultIndex,
  rendered,
  featureStrokeOverrides,
  legendStrokeOverrides,
  originalSvgStroke
}) => {
  const recordKey = text(feature?.recordKey ?? feature?.record_key);
  const biologicalFeatureId = text(
    feature?.biologicalFeatureId ?? feature?.biological_feature_id
  );
  const featureKey = stableFeatureTargetKey({
    ...feature,
    resultKey,
    recordKey,
    biologicalFeatureId
  });
  const svgIds = [...new Set([
    ...(Array.isArray(rendered) ? rendered.map(renderedId) : []),
    renderedId(feature)
  ].filter(Boolean))];
  if (!resultKey || !recordKey || !biologicalFeatureId || !featureKey || svgIds.length === 0) {
    return null;
  }
  const caption = legendCaption(feature);
  const overrideKey = biologicalFeatureKey(recordKey, biologicalFeatureId);
  const exactOverride = featureStrokeOverrides?.[overrideKey] || null;
  const captionOverride = caption ? legendStrokeOverrides?.[caption] || null : null;
  const viewModel = resolveFeatureStrokeViewModel({
    exactOverride,
    captionOverride,
    legendCaption: caption,
    svgDefaultStroke: originalSvgStroke
  });
  const aliases = [...new Set([
    featureKey,
    overrideKey,
    text(feature?.instanceHash),
    text(feature?.targetFeatureKey ?? feature?.target_feature_key),
    ...svgIds
  ].filter(Boolean))];
  return {
    featureKey,
    aliases,
    resultKey,
    resultName,
    resultIndex: resultIndex(feature, fallbackResultIndex),
    recordKey,
    biologicalFeatureId,
    overrideKey,
    instanceHash: isFeatureInstanceHash(feature?.instanceHash)
      ? text(feature.instanceHash)
      : '',
    type: text(feature?.type),
    legendCaption: caption,
    renderedSvgIds: svgIds,
    exactOverride: cloneJson(exactOverride),
    viewModel,
    raw: cloneJson(feature)
  };
};

const catalogFeatures = (snapshot) => {
  if (Array.isArray(snapshot?.features)) {
    return snapshot.features.map((feature) => ({
      feature,
      resultKey: text(feature?.resultKey ?? feature?.result_key),
      resultName: text(feature?.resultName ?? feature?.result_name),
      fallbackResultIndex: resultIndex(feature),
      rendered: [feature]
    }));
  }
  return (Array.isArray(snapshot?.catalog?.items) ? snapshot.catalog.items : [])
    .flatMap((item, itemIndex) => {
      const resultKey = catalogResultKey(item);
      const resultName = text(item?.resultName);
      const renderedByPair = new Map();
      (Array.isArray(item?.features) ? item.features : []).forEach((entry) => {
        const key = biologicalFeatureKey(entry?.recordKey, entry?.biologicalFeatureId);
        if (!key) return;
        if (!renderedByPair.has(key)) renderedByPair.set(key, []);
        renderedByPair.get(key).push(entry);
      });
      return (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : [])
        .map((feature) => ({
          feature,
          resultKey,
          resultName,
          fallbackResultIndex: itemIndex,
          rendered: renderedByPair.get(biologicalFeatureKey(
            feature?.recordKey,
            feature?.biologicalFeatureId
          )) || []
        }));
    });
};

const normalizeSnapshot = (snapshot) => {
  const featureStrokeOverrides = cloneJson(snapshot?.featureStrokeOverrides || {});
  const legendStrokeOverrides = cloneJson(snapshot?.legendStrokeOverrides || {});
  const originalSvgStroke = cloneJson(snapshot?.originalSvgStroke || {});
  const byKey = new Map();
  catalogFeatures(snapshot).forEach((entry) => {
    const normalized = normalizeFeature({
      ...entry,
      featureStrokeOverrides,
      legendStrokeOverrides,
      originalSvgStroke
    });
    if (normalized && !byKey.has(normalized.featureKey)) {
      byKey.set(normalized.featureKey, normalized);
    }
  });
  const features = [...byKey.values()].sort((left, right) => (
    left.resultKey.localeCompare(right.resultKey)
    || left.featureKey.localeCompare(right.featureKey)
  ));
  return {
    features,
    featureStrokeOverrides,
    legendStrokeOverrides,
    originalSvgStroke,
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

export const normalizeFeatureStrokeIntent = (intent = {}) => ({
  targetFeatureKey: text(intent?.targetFeatureKey) || text(intent?.targetFeatureKeys?.[0]),
  targetFeatureKeys: Array.isArray(intent?.targetFeatureKeys)
    ? intent.targetFeatureKeys.map(text).filter(Boolean)
    : [],
  value: normalizeFeatureStrokeValue(intent?.value),
  source: text(intent?.source)
    ? (VALID_SOURCES.has(text(intent.source)) ? text(intent.source) : null)
    : 'popup',
  originResultKey: text(intent?.originResultKey),
  originResultIndex: resultIndex(intent, intent?.originResultIndex),
  resultGenerationKey: Number(intent?.resultGenerationKey ?? 0),
  documentEpoch: Number(intent?.documentEpoch ?? 0),
  styleRevision: Number(intent?.styleRevision ?? intent?.semanticStyleRevision ?? 0),
  styleFingerprint: text(intent?.styleFingerprint),
  requestedSemanticScope: normalizeFeatureStrokeSemanticScope(
    intent?.semanticScope ?? intent?.requestedSemanticScope
  ),
  token: text(intent?.token)
});

const findFeature = (features, key) => features.find(
  (feature) => feature.featureKey === key || feature.aliases.includes(key)
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
    return `Update only this feature in ${target.resultName || 'this output'}`;
  }
  const prefix = scope === 'legend-caption'
    ? `Update legend caption “${target.legendCaption}”`
    : 'Update selected features';
  return `${prefix} — ${count} features in ${resultCount} outputs`;
};

const buildCandidate = (scope, target, rawTargets) => {
  const targets = [...rawTargets].sort((left, right) => (
    left.resultKey.localeCompare(right.resultKey)
    || left.featureKey.localeCompare(right.featureKey)
  ));
  const grouped = targetsByResult(targets);
  const signature = targetSetSignature(targets.map((feature) => feature.featureKey));
  return {
    id: `${scope}:${signatureHash(signature)}`,
    semanticScope: scope,
    resultExtent: scope === 'single' ? 'current-result' : 'all-affected-results',
    label: candidateLabel(scope, target, targets.length, grouped.length),
    targetCount: targets.length,
    affectedResultCount: grouped.length,
    targetsByResult: grouped,
    targetFeatureKeys: targets.map((feature) => feature.featureKey),
    targetSetSignature: signature,
    durableStrokeIntent: scope === 'legend-caption'
      ? { kind: 'legend-caption', caption: target.legendCaption }
      : {
          kind: 'exact-feature-overrides',
          selectors: targets.map((feature) => ({
            featureKey: feature.featureKey,
            instanceHash: feature.instanceHash,
            legendCaption: feature.legendCaption
          }))
        }
  };
};

const invalidPlan = (intent, diagnostics) => immutable({
  token: intent.token || `stroke-plan:${signatureHash({ intent, diagnostics })}`,
  status: 'invalid',
  intent,
  candidates: [],
  diagnostics
});

export const planFeatureStrokeChange = (rawSnapshot = {}, rawIntent = {}) => {
  const snapshot = normalizeSnapshot(rawSnapshot);
  const intent = normalizeFeatureStrokeIntent(rawIntent);
  const rawScope = text(rawIntent?.semanticScope ?? rawIntent?.requestedSemanticScope);
  if (text(rawIntent?.source) && !intent.source) {
    return invalidPlan(intent, [{ code: 'unsupported-source', message: 'Feature stroke source is unavailable' }]);
  }
  if (rawScope && !intent.requestedSemanticScope) {
    return invalidPlan(intent, [{
      code: rawScope === 'matching-rule' ? 'unsupported-matching-rule' : 'unsupported-scope',
      message: `Stroke scope ${rawScope} is unavailable`
    }]);
  }
  if (!intent.value) {
    return invalidPlan(intent, [{ code: 'invalid-value', message: 'Invalid Feature stroke value' }]);
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
    return invalidPlan(intent, [{ code: 'missing-target', message: 'Feature target is unavailable' }]);
  }
  if (
    (intent.originResultKey && intent.originResultKey !== target.resultKey)
    || (
      intent.originResultIndex !== null
      && intent.originResultIndex !== target.resultIndex
    )
  ) {
    return invalidPlan(intent, [{ code: 'stale-result', message: 'Origin Result changed' }]);
  }

  const candidates = [];
  if (target.legendCaption) {
    const captionTargets = snapshot.features.filter(
      (feature) => captionKey(feature.legendCaption) === captionKey(target.legendCaption)
    );
    if (captionTargets.length > 1) {
      candidates.push(buildCandidate('legend-caption', target, captionTargets));
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
    FEATURE_STROKE_SCOPE_ORDER.indexOf(left.semanticScope)
      - FEATURE_STROKE_SCOPE_ORDER.indexOf(right.semanticScope)
  ));

  const diagnostics = [];
  if (!exactAvailable || (selected.length > 1 && !selectedExact)) {
    diagnostics.push({ code: 'exact-scope-unavailable', message: EXACT_SCOPE_DIAGNOSTIC });
  }
  if (candidates.length === 0) {
    return invalidPlan(intent, diagnostics.length ? diagnostics : [{
      code: 'no-targets',
      message: 'No durable Feature stroke target is available'
    }]);
  }
  const requested = intent.requestedSemanticScope
    ? candidates.find((candidate) => candidate.semanticScope === intent.requestedSemanticScope)
    : null;
  if (intent.requestedSemanticScope && !requested) {
    return invalidPlan(intent, [...diagnostics, {
      code: 'unsupported-scope',
      message: `Stroke scope ${intent.requestedSemanticScope} is unavailable`
    }]);
  }
  const selectedCandidate = requested
    || (intent.source === 'selection-toolbar'
      ? candidates.find((candidate) => candidate.semanticScope === 'selected-features')
      : null)
    || (candidates.every((candidate) => candidate.semanticScope === 'single')
      ? candidates[0]
      : null);
  const token = intent.token || `stroke-plan:${signatureHash({
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
    strokeStateFingerprint: featureStrokeStateFingerprint(snapshot),
    ...(selectedCandidate ? { selectedCandidateId: selectedCandidate.id } : {})
  });
};

export const resolveFeatureStrokePlan = (
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
        code: 'missing-scope',
        message: 'The selected Feature stroke scope is unavailable'
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
      strokeStateFingerprint: plan.strokeStateFingerprint
    },
    semanticAfter: {
      value: plan.intent.value,
      durableStrokeIntent: candidate.durableStrokeIntent
    }
  });
};

export { stableFeatureTargetKey, stableFeatureFillStringify };
