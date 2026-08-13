import { normalizeSpecificRule } from '../specific-color-rules.js';
import {
  biologicalFeatureKey,
  catalogResultKey
} from '../../services/feature-catalog.js';
import {
  FEATURE_INSTANCE_HASH_QUALIFIER,
  normalizeFeatureFillValue,
  resolveFeatureFillViewModel,
  resolveOrderedSpecificRule
} from './fill-view-model.js';
import {
  specificRuleIdentitySignature,
  stableFeatureTargetKey
} from './fill-scope-plan.js';
import {
  normalizeStyleSnapshot,
  styleFingerprint
} from '../../services/style-revision.js';
import {
  buildStyleSnapshotCommand,
  exactJsonByteLength,
  immutableStyleCommandSnapshot
} from './style-snapshot-command.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
import { isFeatureInstanceHash } from '../../services/feature-instance-identity.js';
import {
  FEATURE_SEMANTIC_SCOPE_QUALIFIER,
  encodeFeatureSemanticSelector,
  parseFeatureSemanticSelector
} from './semantic-fill-selectors.js';

const text = (value) => String(value ?? '').trim();

const cloneRules = (rules) => (Array.isArray(rules) ? rules : []).map((rule) => ({
  ...normalizeSpecificRule(rule, { fromFile: Boolean(rule?.fromFile) })
}));

const cloneOverrides = (overrides) => Object.fromEntries(
  Object.entries(overrides || {}).map(([key, value]) => [key, cloneJson(value)])
);

const readRef = (value) => (
  value && typeof value === 'object' && 'value' in value
    ? value.value
    : value
);

const writeRef = (target, value) => {
  if (target && typeof target === 'object' && 'value' in target) {
    target.value = value;
    return;
  }
  throw new Error('Feature fill command requires a writable ref.');
};

const isWritableRef = (value) => (
  value
  && typeof value === 'object'
  && 'value' in value
  && value.__v_isReadonly !== true
);

const replaceObject = (target, source) => {
  Object.keys(target || {}).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = cloneJson(value);
  });
};

const replaceRules = (target, rules) => {
  if (!Array.isArray(target)) throw new Error('Feature fill rules are unavailable.');
  target.splice(0, target.length, ...cloneRules(rules));
};

const restoreRuleReferences = (target, rules) => {
  if (!Array.isArray(target) || !Array.isArray(rules)) {
    throw new Error('Feature fill rules cannot be restored exactly.');
  }
  target.splice(0, target.length, ...rules);
};

const restoreObjectReferences = (target, source) => {
  Object.keys(target || {}).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = value;
  });
};

const catalogBindings = (catalog) => {
  const byResultKey = new Map();
  const byTargetKey = new Map();
  (Array.isArray(catalog?.items) ? catalog.items : []).forEach((item, resultIndex) => {
    const resultKey = catalogResultKey(item);
    if (!resultKey || byResultKey.has(resultKey)) {
      throw new Error('Feature fill catalogue has missing or duplicate Result identity.');
    }
    const featureByPair = new Map();
    const renderedByPair = new Map();
    (Array.isArray(item?.features) ? item.features : []).forEach((feature) => {
      const pairKey = biologicalFeatureKey(feature?.recordKey, feature?.biologicalFeatureId);
      if (!pairKey) return;
      if (!renderedByPair.has(pairKey)) renderedByPair.set(pairKey, []);
      renderedByPair.get(pairKey).push(feature);
    });
    (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).forEach((feature) => {
      const pairKey = biologicalFeatureKey(feature?.recordKey, feature?.biologicalFeatureId);
      const targetKey = stableFeatureTargetKey({ ...feature, resultKey });
      if (!pairKey || !targetKey || featureByPair.has(pairKey) || byTargetKey.has(targetKey)) {
        throw new Error('Feature fill catalogue has ambiguous target identity.');
      }
      const rendered = renderedByPair.get(pairKey) || [];
      const binding = { resultKey, resultIndex, item, feature, rendered, pairKey, targetKey };
      featureByPair.set(pairKey, binding);
      byTargetKey.set(targetKey, binding);
    });
    byResultKey.set(resultKey, { resultKey, resultIndex, item, featureByPair });
  });
  return { byResultKey, byTargetKey };
};

const validatePlanTargets = (plan, bindings) => {
  if (!plan || plan.status !== 'ready') {
    throw new Error('A resolved ready Feature fill plan is required.');
  }
  const affected = [...new Set((plan.affectedResultKeys || []).map(text).filter(Boolean))];
  if (affected.length === 0 || affected.length !== (plan.affectedResultKeys || []).length) {
    throw new Error('Feature fill plan has invalid affected Result identities.');
  }
  const groups = Array.isArray(plan.targetsByResult) ? plan.targetsByResult : [];
  const planned = new Set();
  const plannedResults = new Set();
  groups.forEach((group) => {
    const resultKey = text(group?.resultKey);
    if (
      !affected.includes(resultKey)
      || !bindings.byResultKey.has(resultKey)
      || plannedResults.has(resultKey)
    ) {
      throw new Error(`Feature fill target Result is unavailable: ${resultKey}`);
    }
    plannedResults.add(resultKey);
    const featureKeys = Array.isArray(group?.featureKeys) ? group.featureKeys : [];
    if (featureKeys.length === 0) {
      throw new Error(`Feature fill target Result has no targets: ${resultKey}`);
    }
    featureKeys.forEach((featureKey) => {
      const binding = bindings.byTargetKey.get(text(featureKey));
      if (!binding || binding.resultKey !== resultKey || planned.has(binding.targetKey)) {
        throw new Error(`Feature fill target is stale or cross-Result: ${text(featureKey)}`);
      }
      planned.add(binding.targetKey);
    });
  });
  if (
    planned.size === 0
    || plannedResults.size !== affected.length
    || affected.some((resultKey) => !plannedResults.has(resultKey))
  ) {
    throw new Error('Feature fill plan does not cover every affected Result exactly once.');
  }
  return { affected, targetKeys: [...planned].sort() };
};

const valueColor = (value) => {
  const normalized = normalizeFeatureFillValue(value);
  if (!normalized) throw new Error('Feature fill plan has an invalid value.');
  if (normalized.kind === 'inherit') return null;
  return normalized.kind === 'none' ? 'none' : normalized.color;
};

const exactRule = ({ selector, value, caption }) => ({
  feat: text(selector?.feat),
  qual: FEATURE_INSTANCE_HASH_QUALIFIER,
  val: text(selector?.val),
  color: valueColor(value),
  cap: text(caption),
  match: 'literal'
});

const targetRuleMatchesSelector = (rule, selector) => (
  text(rule?.feat) === text(selector?.feat)
  && text(rule?.qual) === FEATURE_INSTANCE_HASH_QUALIFIER
  && text(rule?.val) === text(selector?.val)
  && text(rule?.match).toLowerCase() !== 'regex'
);

const targetRuleMatchesSemanticSelector = (rule, selector, featureType = null) => (
  (featureType === null || text(rule?.feat) === text(featureType))
  && text(rule?.qual) === FEATURE_SEMANTIC_SCOPE_QUALIFIER
  && text(rule?.val) === text(selector)
  && text(rule?.match).toLowerCase() !== 'regex'
);

const applyExactRuleIntent = ({ rules, intent, value, caption }) => {
  const selectors = Array.isArray(intent?.selectors) ? intent.selectors : [];
  if (selectors.length === 0) throw new Error('Exact Feature fill scope has no selectors.');
  const next = cloneRules(rules);
  selectors.forEach((selector) => {
    const ownedIndexes = [];
    next.forEach((rule, index) => {
      if (targetRuleMatchesSelector(rule, selector)) ownedIndexes.push(index);
    });
    if (value.kind === 'inherit') {
      if (ownedIndexes.length > 0) next.splice(ownedIndexes[0], 1);
      return;
    }
    const current = ownedIndexes.length > 0 ? next[ownedIndexes[0]] : null;
    const replacement = {
      ...exactRule({
        selector,
        value,
        caption: text(caption) || text(current?.cap) || text(selector?.caption)
      }),
      ...(current?.fromFile ? { fromFile: true } : {})
    };
    if (!replacement.feat || !replacement.val) {
      throw new Error('Exact Feature fill selector is incomplete.');
    }
    if (ownedIndexes.length > 0) {
      next.splice(ownedIndexes[0], 1, replacement);
    } else {
      next.unshift(replacement);
    }
  });
  return next;
};

const updateOwnedRule = ({ rules, intent, value, caption }) => {
  const next = cloneRules(rules);
  const index = Number(intent?.ruleIndex);
  const current = next[index];
  if (!Number.isInteger(index) || !current) {
    throw new Error('The selected specific-color rule is unavailable.');
  }
  if (specificRuleIdentitySignature(current, index) !== text(intent?.ruleIdentity)) {
    throw new Error('The selected specific-color rule changed before commit.');
  }
  if (value.kind === 'inherit') {
    next.splice(index, 1);
    return next;
  }
  next.splice(index, 1, {
    ...current,
    color: valueColor(value),
    cap: text(caption) || text(current.cap)
  });
  return next;
};

const semanticSelectorForIntent = (intent) => {
  const kind = text(intent?.kind);
  if (kind === 'semantic-selector') {
    const parsed = parseFeatureSemanticSelector(intent?.selector);
    if (!parsed) return null;
    return {
      kind: parsed.kind,
      value: parsed.value,
      ...(parsed.kind === 'feature-type'
        ? { primaryFeatureType: parsed.value }
        : {})
    };
  }
  if (kind === 'feature-type') {
    return {
      kind,
      value: text(intent?.featureType),
      primaryFeatureType: text(intent?.featureType)
    };
  }
  if (kind === 'legend-caption') {
    return { kind: 'base-legend-caption', value: text(intent?.caption) };
  }
  if (kind === 'rendered-label') {
    return { kind, value: text(intent?.label) };
  }
  if (kind === 'source-annotation-label') {
    return { kind, value: text(intent?.label) };
  }
  if (kind === 'similarity-group') {
    return { kind, value: text(intent?.groupId) };
  }
  return null;
};

const semanticRuleFeatureTypes = ({ selector, targetKeys, bindings, rules }) => {
  if (selector.primaryFeatureType) return [selector.primaryFeatureType];

  // A wildcard semantic row is the durable owner for cross-type scopes. The
  // renderer intentionally resolves an exact Feature-type rule set before the
  // wildcard set, so add one bounded exact-type companion only when a current
  // target is already captured by an exact-type rule. Rule growth therefore
  // follows distinct selector precedence barriers, never Feature count.
  const blockedFeatureTypes = new Set();
  targetKeys.forEach((featureKey) => {
    const feature = bindings.byTargetKey.get(featureKey)?.feature;
    const featureType = text(feature?.type);
    const resolved = resolveOrderedSpecificRule(feature, rules);
    if (featureType && text(resolved?.rule?.feat) === featureType) {
      blockedFeatureTypes.add(featureType);
    }
  });
  return ['*', ...blockedFeatureTypes].filter(Boolean);
};

const applySemanticRuleIntent = ({
  rules,
  intent,
  value,
  caption,
  targetKeys,
  bindings
}) => {
  const selector = semanticSelectorForIntent(intent);
  if (!selector?.value) {
    throw new Error(`Feature fill semantic scope ${text(intent?.kind)} is incomplete.`);
  }
  const literal = encodeFeatureSemanticSelector(selector.kind, selector.value);
  const next = cloneRules(rules);
  const owned = next.filter((rule) => targetRuleMatchesSemanticSelector(rule, literal));
  for (let index = next.length - 1; index >= 0; index -= 1) {
    if (targetRuleMatchesSemanticSelector(next[index], literal)) next.splice(index, 1);
  }
  if (value.kind === 'inherit') return next;

  const featureTypes = semanticRuleFeatureTypes({
    selector,
    targetKeys,
    bindings,
    rules: next
  });
  const effectiveCaption = text(caption) || text(owned[0]?.cap) || selector.value;
  const replacements = featureTypes.map((featureType) => {
    const previous = owned.find(
      (rule) => targetRuleMatchesSemanticSelector(rule, literal, featureType)
    );
    return {
      feat: featureType,
      qual: FEATURE_SEMANTIC_SCOPE_QUALIFIER,
      val: literal,
      color: valueColor(value),
      cap: effectiveCaption,
      match: 'literal',
      ...(previous?.fromFile ? { fromFile: true } : {})
    };
  });
  next.unshift(...replacements);
  return next;
};

const buildNextRules = ({ plan, beforeRules, targetKeys, bindings }) => {
  const value = normalizeFeatureFillValue(plan?.semanticAfter?.value ?? plan?.intent?.value);
  if (!value) throw new Error('Feature fill plan has an invalid value.');
  const caption = text(
    plan?.semanticAfter?.requestedCaption
    ?? plan?.intent?.requestedCaption
  );
  const intent = plan?.semanticAfter?.durableRuleIntent;
  if (intent?.kind === 'update-specific-rule') {
    return updateOwnedRule({ rules: beforeRules, intent, value, caption });
  }
  if (intent?.kind === 'exact-feature-rules') {
    const plannedSelectors = Array.isArray(intent.selectors) ? intent.selectors : [];
    const plannedKeys = new Set(targetKeys);
    if (
      plannedSelectors.length !== plannedKeys.size
      || plannedSelectors.some((selector) => (
        !plannedKeys.has(text(selector?.featureKey))
        || !isFeatureInstanceHash(selector?.val)
      ))
    ) {
      throw new Error('Exact Feature fill selectors do not match the resolved plan.');
    }
    return applyExactRuleIntent({ rules: beforeRules, intent, value, caption });
  }
  return applySemanticRuleIntent({
    rules: beforeRules,
    intent,
    value,
    caption,
    targetKeys,
    bindings
  });
};

const paletteSnapshot = (state, explicit = null) => ({
  appliedPaletteName: text(explicit?.appliedPaletteName ?? readRef(state?.appliedPaletteName)) || 'default',
  appliedPaletteColors: cloneJson(
    explicit?.appliedPaletteColors ?? readRef(state?.appliedPaletteColors) ?? {}
  )
});

const semanticSide = ({
  rules,
  palette,
  overrides,
  resultFingerprintByKey = {},
  expectedFingerprint = ''
}) => {
  const semantic = normalizeStyleSnapshot({ rules, ...palette });
  return immutableStyleCommandSnapshot({
    rules: cloneRules(rules),
    ...palette,
    overrides: cloneOverrides(overrides),
    fingerprint: styleFingerprint(semantic),
    expectedFingerprint: text(expectedFingerprint),
    resultFingerprintByKey: cloneJson(resultFingerprintByKey)
  });
};

export const deriveFeatureFillOverrides = ({
  catalog,
  rules = [],
  paletteColors = {},
  svgDefaultColor = '#cccccc'
} = {}) => {
  const overrides = {};
  const effectiveByTargetKey = {};
  const counters = { features: 0, ruleResolutionUpperBound: 0, explicitOverrides: 0 };
  const normalizedRules = cloneRules(rules);
  const bindings = catalogBindings(catalog);
  [...bindings.byTargetKey.values()]
    .sort((left, right) => left.targetKey.localeCompare(right.targetKey))
    .forEach((binding) => {
      counters.features += 1;
      counters.ruleResolutionUpperBound += normalizedRules.length;
      const resolved = resolveOrderedSpecificRule(binding.feature, normalizedRules);
      const model = resolveFeatureFillViewModel({
        feature: binding.feature,
        specificRules: normalizedRules,
        paletteColors,
        svgDefaultColor
      });
      const caption = text(
        resolved?.rule?.cap
        || binding.feature?.type
        || binding.feature?.effectiveLegendCaption
        || binding.feature?.effective_legend_caption
        || binding.feature?.legendCaption
        || binding.feature?.legend_caption
      );
      const value = { color: model.effectiveColor, caption };
      overrides[binding.pairKey] = value;
      effectiveByTargetKey[binding.targetKey] = value;
      if (model.explicitValue !== null) counters.explicitOverrides += 1;
    });
  return immutableStyleCommandSnapshot({ overrides, effectiveByTargetKey, counters });
};

export const buildFeatureFillSemanticSnapshots = ({
  plan,
  catalog,
  rules = [],
  appliedPaletteName = 'default',
  appliedPaletteColors = {},
  resultFingerprintByKey = {},
  svgDefaultColor = '#cccccc'
} = {}) => {
  const bindings = catalogBindings(catalog);
  const targets = validatePlanTargets(plan, bindings);
  const beforeRules = cloneRules(rules);
  const nextRules = buildNextRules({
    plan,
    beforeRules,
    targetKeys: targets.targetKeys,
    bindings
  });
  const palette = { appliedPaletteName, appliedPaletteColors };
  const beforeDerived = deriveFeatureFillOverrides({
    catalog,
    rules: beforeRules,
    paletteColors: appliedPaletteColors,
    svgDefaultColor
  });
  const afterDerived = deriveFeatureFillOverrides({
    catalog,
    rules: nextRules,
    paletteColors: appliedPaletteColors,
    svgDefaultColor
  });
  const plannedTargets = new Set(targets.targetKeys);
  Object.keys(afterDerived.effectiveByTargetKey).forEach((targetKey) => {
    const beforeValue = beforeDerived.effectiveByTargetKey[targetKey];
    const afterValue = afterDerived.effectiveByTargetKey[targetKey];
    if (
      !plannedTargets.has(targetKey)
      && (
        beforeValue?.color !== afterValue?.color
        || beforeValue?.caption !== afterValue?.caption
      )
    ) {
      throw new Error(`Feature fill rule intent changed an unplanned target: ${targetKey}`);
    }
  });
  const before = semanticSide({
    rules: beforeRules,
    palette,
    overrides: beforeDerived.overrides,
    resultFingerprintByKey,
    expectedFingerprint: plan?.semanticBefore?.styleFingerprint
  });
  const afterFingerprint = styleFingerprint({ rules: nextRules, ...palette });
  const nextLedger = Object.fromEntries(
    [...bindings.byResultKey.keys()].map((resultKey) => [resultKey, afterFingerprint])
  );
  const after = semanticSide({
    rules: nextRules,
    palette,
    overrides: afterDerived.overrides,
    resultFingerprintByKey: nextLedger
  });
  return immutableStyleCommandSnapshot({
    before,
    after,
    affectedResultKeys: targets.affected,
    targetFeatureKeys: targets.targetKeys,
    fillsBeforeByTargetKey: beforeDerived.effectiveByTargetKey,
    fillsAfterByTargetKey: afterDerived.effectiveByTargetKey,
    counters: {
      targetFeatures: targets.targetKeys.length,
      affectedResults: targets.affected.length,
      beforeRuleResolutionUpperBound: beforeDerived.counters.ruleResolutionUpperBound,
      afterRuleResolutionUpperBound: afterDerived.counters.ruleResolutionUpperBound
    }
  });
};

const stateFingerprint = (state) => styleFingerprint({
  rules: state?.manualSpecificRules,
  ...paletteSnapshot(state)
});

const resultKeySetIsCurrent = (results, catalog, affectedResultKeys) => {
  if (!Array.isArray(results) || !Array.isArray(catalog?.items) || results.length !== catalog.items.length) {
    return false;
  }
  const keys = new Set();
  for (let resultIndex = 0; resultIndex < catalog.items.length; resultIndex += 1) {
    const item = catalog.items[resultIndex];
    const result = results[resultIndex];
    const resultKey = catalogResultKey(item);
    if (
      !resultKey
      || keys.has(resultKey)
      || !result
      || (item?.resultIndex !== undefined && Number(item.resultIndex) !== resultIndex)
      || (text(item?.resultName) && text(item.resultName) !== text(result?.name))
    ) return false;
    keys.add(resultKey);
  }
  return affectedResultKeys.every((key) => keys.has(key));
};

const resultLedgerIsCurrent = (ledger, catalog, fingerprint) => (
  (Array.isArray(catalog?.items) ? catalog.items : []).every((item) => {
    const resultKey = catalogResultKey(item);
    return resultKey && text(ledger?.[resultKey]) === fingerprint;
  })
);

const targetBindingIdentity = (binding) => JSON.stringify({
  type: text(binding?.feature?.type),
  instanceHash: text(binding?.feature?.instanceHash),
  qualifiers: binding?.feature?.qualifiers || null,
  effectiveLegendCaption: text(
    binding?.feature?.effectiveLegendCaption
    ?? binding?.feature?.effective_legend_caption
    ?? binding?.feature?.legendCaption
    ?? binding?.feature?.legend_caption
  ),
  renderedLabel: text(binding?.feature?.renderedLabel ?? binding?.feature?.rendered_label),
  sourceAnnotationLabel: text(
    binding?.feature?.sourceAnnotationLabel ?? binding?.feature?.source_annotation_label
  ),
  hidden: binding?.feature?.hidden === true,
  renderedFlag: binding?.feature?.rendered !== false,
  renderedSvgIds: (Array.isArray(binding?.rendered) ? binding.rendered : [])
    .map((feature) => text(feature?.svgId ?? feature?.svg_id))
    .filter(Boolean)
    .sort()
});

const targetSetIsCurrent = (plan, catalog, expectedIdentityByTargetKey = null) => {
  try {
    const bindings = catalogBindings(catalog);
    const targets = validatePlanTargets(plan, bindings);
    if (expectedIdentityByTargetKey) {
      for (const targetKey of targets.targetKeys) {
        if (
          text(expectedIdentityByTargetKey[targetKey])
          !== targetBindingIdentity(bindings.byTargetKey.get(targetKey))
        ) return false;
      }
    }
    return true;
  } catch {
    return false;
  }
};

const currentMounted = (getMountedContext) => {
  const context = typeof getMountedContext === 'function' ? getMountedContext() : null;
  if (!context) return { resultIndex: null, resultKey: '', svg: null };
  return {
    ...context,
    resultIndex: context.resultIndex !== null
      && context.resultIndex !== undefined
      && Number.isInteger(Number(context.resultIndex))
      && Number(context.resultIndex) >= 0
      ? Number(context.resultIndex)
      : null,
    resultKey: text(context.resultKey),
    svg: context.svg || null
  };
};

const sameMountedOwnership = (left, right) => (
  (left?.resultIndex ?? null) === (right?.resultIndex ?? null)
  && text(left?.resultKey) === text(right?.resultKey)
  && (left?.svg || null) === (right?.svg || null)
);

const mountedContextIsCurrent = (mounted, catalog) => {
  if (mounted?.resultIndex === null || mounted?.resultIndex === undefined) {
    return !mounted?.resultKey && !mounted?.svg;
  }
  return (
    mounted.resultIndex < (catalog?.items?.length || 0)
    && catalogResultKey(catalog.items[mounted.resultIndex]) === mounted.resultKey
  );
};

const resultContentSignature = (result) => (
  String(result?.content ?? '')
);

const projectedResultsMatch = (published, projected) => (
  Array.isArray(published)
  && Array.isArray(projected)
  && published.length === projected.length
  && published.every((result, index) => (
    text(result?.name) === text(projected[index]?.name)
    && resultContentSignature(result) === resultContentSignature(projected[index])
  ))
);

const exactStateSnapshot = ({ state, captureRuntimeState }) => {
  const runtime = typeof captureRuntimeState === 'function' ? captureRuntimeState() : null;
  if (runtime && typeof runtime.then === 'function') {
    throw new Error('Feature fill runtime capture must be synchronous.');
  }
  return {
    rules: [...state.manualSpecificRules],
    appliedPaletteName: cloneJson(readRef(state.appliedPaletteName)),
    appliedPaletteColors: readRef(state.appliedPaletteColors),
    overrides: { ...state.featureColorOverrides },
    results: readRef(state.results),
    legendEntries: readRef(state.legendEntries) || [],
    ledger: readRef(state.validatedStyleFingerprintByResultKey) || {},
    revision: Number(readRef(state.semanticStyleRevision) || 0),
    fingerprint: text(readRef(state.semanticStyleFingerprint)),
    presentation: {
      ...(isWritableRef(state.clickedFeature)
        ? { clickedFeature: readRef(state.clickedFeature) }
        : {}),
      ...(isWritableRef(state.selectedFeatures)
        ? { selectedFeatures: readRef(state.selectedFeatures) }
        : {}),
      ...(isWritableRef(state.selectedResultIndex)
        ? { selectedResultIndex: readRef(state.selectedResultIndex) }
        : {})
    },
    runtime
  };
};

const restoreStateSnapshot = ({ state, snapshot, restoreRuntimeState }) => {
  restoreRuleReferences(state.manualSpecificRules, snapshot.rules);
  writeRef(state.appliedPaletteName, cloneJson(snapshot.appliedPaletteName));
  writeRef(state.appliedPaletteColors, snapshot.appliedPaletteColors);
  restoreObjectReferences(state.featureColorOverrides, snapshot.overrides);
  writeRef(state.results, snapshot.results);
  writeRef(state.legendEntries, snapshot.legendEntries);
  writeRef(state.validatedStyleFingerprintByResultKey, snapshot.ledger);
  writeRef(state.semanticStyleRevision, snapshot.revision);
  writeRef(state.semanticStyleFingerprint, snapshot.fingerprint);
  if (
    Object.prototype.hasOwnProperty.call(snapshot.presentation || {}, 'clickedFeature')
    && isWritableRef(state.clickedFeature)
  ) {
    writeRef(state.clickedFeature, snapshot.presentation.clickedFeature);
  }
  if (
    Object.prototype.hasOwnProperty.call(snapshot.presentation || {}, 'selectedFeatures')
    && isWritableRef(state.selectedFeatures)
  ) {
    writeRef(state.selectedFeatures, snapshot.presentation.selectedFeatures);
  }
  if (
    Object.prototype.hasOwnProperty.call(snapshot.presentation || {}, 'selectedResultIndex')
    && isWritableRef(state.selectedResultIndex)
  ) {
    writeRef(state.selectedResultIndex, snapshot.presentation.selectedResultIndex);
  }
  if (typeof restoreRuntimeState === 'function') {
    const restored = restoreRuntimeState(snapshot.runtime);
    if (restored && typeof restored.then === 'function') {
      throw new Error('Feature fill runtime rollback must be synchronous.');
    }
    if (restored === false) return false;
  }
  return true;
};

/**
 * Build one fill-specific History command. Result and legend preparation are
 * injected so the command owns atomicity without owning DOM parsing/layout.
 */
export const buildFeatureFillCommand = async ({
  plan,
  state,
  catalog,
  prepareResultProjection,
  prepareLegendProjection = null,
  getMountedContext = null,
  commitMountedProjection = null,
  restoreMountedProjection = null,
  captureRuntimeState = null,
  restoreRuntimeState = null,
  reconcile = null,
  refreshPresentation = null,
  label = 'Change feature fill',
  svgDefaultColor = '#cccccc'
} = {}) => {
  if (!state || typeof prepareResultProjection !== 'function') {
    throw new Error('Feature fill command requires state and a Result projection owner.');
  }
  const expectedEpoch = Number(plan?.intent?.documentEpoch);
  const expectedGeneration = Number(plan?.intent?.resultGenerationKey);
  const expectedFingerprint = text(plan?.semanticBefore?.styleFingerprint);
  if (
    Number(readRef(state.documentEpoch) || 0) !== expectedEpoch
    || Number(readRef(state.resultGenerationKey) || 0) !== expectedGeneration
    || (expectedFingerprint && stateFingerprint(state) !== expectedFingerprint)
    || (readRef(state.featureCatalog) && readRef(state.featureCatalog) !== catalog)
  ) {
    throw new Error('The Feature fill plan became stale before command preparation.');
  }
  const currentLedger = readRef(state.validatedStyleFingerprintByResultKey) || {};
  const semantic = buildFeatureFillSemanticSnapshots({
    plan,
    catalog,
    rules: state.manualSpecificRules,
    ...paletteSnapshot(state),
    resultFingerprintByKey: currentLedger,
    svgDefaultColor
  });
  const initialBindings = catalogBindings(catalog);
  const targetIdentityByTargetKey = Object.fromEntries(
    semantic.targetFeatureKeys.map((targetKey) => [
      targetKey,
      targetBindingIdentity(initialBindings.byTargetKey.get(targetKey))
    ])
  );
  const validateCurrent = ({ from, prepared = null }) => {
    const ledger = readRef(state.validatedStyleFingerprintByResultKey) || {};
    if (
      Number(readRef(state.documentEpoch) || 0) !== expectedEpoch
      || Number(readRef(state.resultGenerationKey) || 0) !== expectedGeneration
      || (from.expectedFingerprint && from.expectedFingerprint !== from.fingerprint)
      || stateFingerprint(state) !== from.fingerprint
      || JSON.stringify(cloneRules(state.manualSpecificRules)) !== JSON.stringify(from.rules)
      || text(readRef(state.semanticStyleFingerprint)) !== from.fingerprint
      || !resultLedgerIsCurrent(ledger, catalog, from.fingerprint)
      || (readRef(state.featureCatalog) && readRef(state.featureCatalog) !== catalog)
      || !resultKeySetIsCurrent(readRef(state.results), catalog, semantic.affectedResultKeys)
      || !targetSetIsCurrent(plan, catalog, targetIdentityByTargetKey)
    ) return false;
    const mounted = currentMounted(getMountedContext);
    if (!mountedContextIsCurrent(mounted, catalog)) return false;
    if (
      prepared
      && (
        readRef(state.results) !== prepared.sourceResults
        || !sameMountedOwnership(mounted, prepared.preflightMounted)
        || prepared.sourceResultRecords.some(({ resultIndex, result, content }) => (
          readRef(state.results)?.[resultIndex] !== result
          || resultContentSignature(result) !== content
        ))
        || readRef(state.legendEntries) !== prepared.sourceLegendEntries
        || JSON.stringify(readRef(state.legendEntries) || []) !== prepared.sourceLegendSignature
      )
    ) return false;
    if (mounted.resultKey && semantic.affectedResultKeys.includes(mounted.resultKey) && !mounted.svg) {
      return false;
    }
    return true;
  };

  const prepareTransition = async ({ from, to, direction }) => {
    const mounted = currentMounted(getMountedContext);
    const sourceResults = readRef(state.results);
    const sourceLegendEntries = readRef(state.legendEntries);
    const sourceLegendSignature = JSON.stringify(sourceLegendEntries || []);
    const sourceResultRecords = semantic.affectedResultKeys.map((resultKey) => {
      const resultIndex = catalog.items.findIndex((item) => catalogResultKey(item) === resultKey);
      const result = sourceResults?.[resultIndex];
      return { resultIndex, result, content: resultContentSignature(result) };
    });
    const fillsByTargetKey = direction === 'undo'
      ? semantic.fillsBeforeByTargetKey
      : semantic.fillsAfterByTargetKey;
    const legend = typeof prepareLegendProjection === 'function'
      ? await prepareLegendProjection({
          direction,
          from,
          to,
          plan,
          catalog,
          rules: to.rules,
          paletteColors: to.appliedPaletteColors,
          affectedResultKeys: semantic.affectedResultKeys,
          mounted
        })
      : null;
    const preparedSvgByResultKey = legend?.preparedSvgByResultKey || (
      legend?.staged instanceof Map
        ? new Map([...legend.staged].map(([resultKey, entry]) => [resultKey, entry?.svg || entry]))
        : null
    );
    if (typeof prepareLegendProjection === 'function') {
      for (const resultKey of semantic.affectedResultKeys) {
        const svg = preparedSvgByResultKey instanceof Map
          ? preparedSvgByResultKey.get(resultKey)
          : preparedSvgByResultKey?.[resultKey];
        if (!svg?.cloneNode) {
          throw new Error(`Feature fill legend projection is incomplete for Result ${resultKey}.`);
        }
      }
    }
    const targetFeatureKeysByResult = new Map(
      (Array.isArray(plan?.targetsByResult) ? plan.targetsByResult : []).map((entry) => [
        text(entry?.resultKey),
        cloneJson(entry?.featureKeys || [])
      ])
    );
    const projection = await prepareResultProjection({
      direction,
      from,
      to,
      plan,
      catalog,
      results: sourceResults,
      fillsByTargetKey,
      affectedResultKeys: semantic.affectedResultKeys,
      mounted,
      legend,
      preparedSvgByResultKey,
      targetFeatureKeysByResult
    });
    if (!projection || !Array.isArray(projection.nextResults)) {
      throw new Error('Feature fill Result projection did not prepare a complete Result array.');
    }
    if (readRef(state.results) !== sourceResults) {
      throw new Error('Feature fill Results changed during projection preparation.');
    }
    if (projection.nextResults.length !== sourceResults.length) {
      throw new Error('Feature fill Result projection changed Result topology.');
    }
    const affectedIndexes = new Set(semantic.affectedResultKeys.map((resultKey) => (
      catalog.items.findIndex((item) => catalogResultKey(item) === resultKey)
    )));
    projection.nextResults.forEach((result, resultIndex) => {
      const previous = sourceResults[resultIndex];
      if (!result || text(result?.name) !== text(previous?.name)) {
        throw new Error(`Feature fill Result identity changed at index ${resultIndex}.`);
      }
      if (!affectedIndexes.has(resultIndex) && result !== previous) {
        throw new Error(`Feature fill projection replaced unaffected Result ${resultIndex}.`);
      }
    });
    const resultIndexByKey = new Map(
      catalog.items.map((item, resultIndex) => [catalogResultKey(item), resultIndex])
    );
    const affectedResults = semantic.affectedResultKeys.map((resultKey) => {
      const resultIndex = resultIndexByKey.get(resultKey);
      return {
        resultKey,
        before: sourceResults[resultIndex],
        beforeAdmission: projection.admissionMetadataByResultKey?.[resultKey]?.before || null,
        after: projection.nextResults[resultIndex],
        afterAdmission: projection.admissionMetadataByResultKey?.[resultKey]?.after || null
      };
    });
    const retainedForHistory = {
      affectedResults,
      projection: projection.retainedForHistory || null,
      legend: legend?.retainedForHistory || null
    };
    return {
      projection,
      legend,
      sourceResults,
      sourceResultRecords,
      sourceLegendEntries,
      sourceLegendSignature,
      preflightMounted: mounted,
      counters: {
        projectedResults: Number(projection.counters?.affectedResults || semantic.affectedResultKeys.length),
        mountedSerializations: Number(projection.counters?.mountedSerializations || 0),
        detachedPasses: Number(projection.counters?.detachedPasses || 0),
        changedResults: Number(projection.counters?.changedResults || 0),
        resultArraySwaps: 1,
        legendMeasurements: Number(legend?.counters?.measurements || 0),
        legendAdditions: Number(legend?.counters?.additions || 0),
        legendRenames: Number(legend?.counters?.renames || 0),
        legendRemovals: Number(legend?.counters?.removals || 0),
        legendColorUpdates: Number(legend?.counters?.colorUpdates || 0)
      },
      selectedLegendEntries: cloneJson(
        legend?.selectedEntries
        ?? projection.selectedLegendEntries
        ?? readRef(state.legendEntries)
        ?? []
      ),
      retainedForHistory,
      retainedBytes: exactJsonByteLength(retainedForHistory)
    };
  };

  const commitTransition = ({ to, direction, prepared }) => {
    replaceRules(state.manualSpecificRules, to.rules);
    writeRef(state.appliedPaletteName, to.appliedPaletteName);
    writeRef(state.appliedPaletteColors, cloneJson(to.appliedPaletteColors));
    replaceObject(state.featureColorOverrides, to.overrides);
    if (typeof commitMountedProjection === 'function') {
      const mountedCommitted = commitMountedProjection({
        direction,
        prepared: prepared.projection,
        legend: prepared.legend,
        plan
      });
      if (mountedCommitted && typeof mountedCommitted.then === 'function') {
        throw new Error('Mounted Feature fill commit must be synchronous.');
      }
      if (mountedCommitted === false) return false;
    }
    if (semantic.affectedResultKeys.includes(prepared.preflightMounted.resultKey)) {
      writeRef(state.legendEntries, cloneJson(prepared.selectedLegendEntries));
    }
    writeRef(state.results, prepared.projection.nextResults);
    const publishedResults = readRef(state.results);
    if (!projectedResultsMatch(publishedResults, prepared.projection.nextResults)) {
      throw new Error('Feature fill Result publication did not preserve the staged projection.');
    }
    prepared.publishedResults = publishedResults;
    writeRef(
      state.validatedStyleFingerprintByResultKey,
      Object.freeze(cloneJson(to.resultFingerprintByKey))
    );
    const reconciled = typeof reconcile === 'function'
      ? reconcile({ direction, plan, prepared, to })
      : undefined;
    if (reconciled && typeof reconciled.then === 'function') {
      throw new Error('Feature fill reconciliation must be synchronous.');
    }
    if (reconciled === false) return false;
    const refreshed = typeof refreshPresentation === 'function'
      ? refreshPresentation({ direction, plan, prepared, to })
      : undefined;
    if (refreshed && typeof refreshed.then === 'function') {
      throw new Error('Feature fill presentation refresh must be synchronous.');
    }
    if (refreshed === false) return false;
    writeRef(state.semanticStyleRevision, Number(readRef(state.semanticStyleRevision) || 0) + 1);
    writeRef(state.semanticStyleFingerprint, to.fingerprint);
    return true;
  };

  const captureExactState = () => exactStateSnapshot({ state, captureRuntimeState });
  const restoreExactState = (snapshot, context) => {
    if (typeof restoreMountedProjection === 'function') {
      const restored = restoreMountedProjection({ snapshot, context, plan });
      if (restored && typeof restored.then === 'function') {
        throw new Error('Mounted Feature fill rollback must be synchronous.');
      }
      if (restored === false) return false;
    }
    return restoreStateSnapshot({ state, snapshot, restoreRuntimeState });
  };

  const command = await buildStyleSnapshotCommand({
    label,
    before: semantic.before,
    after: semantic.after,
    metadata: {
      planToken: text(plan?.token),
      documentEpoch: expectedEpoch,
      resultGenerationKey: expectedGeneration,
      affectedResultKeys: semantic.affectedResultKeys,
      targetFeatureKeys: semantic.targetFeatureKeys,
      semanticScope: text(plan?.semanticScope),
      resultExtent: text(plan?.resultExtent)
    },
    counters: semantic.counters,
    isNoop: (before, after) => (
      before.fingerprint === after.fingerprint
      && JSON.stringify(before.overrides) === JSON.stringify(after.overrides)
    ),
    validateCurrent,
    prepareTransition,
    commitTransition,
    captureExactState,
    restoreExactState,
    assertCommitted: ({ to, prepared }) => (
      stateFingerprint(state) === to.fingerprint
      && text(readRef(state.semanticStyleFingerprint)) === to.fingerprint
      && resultLedgerIsCurrent(
        readRef(state.validatedStyleFingerprintByResultKey) || {},
        catalog,
        to.fingerprint
      )
      && JSON.stringify(cloneRules(state.manualSpecificRules)) === JSON.stringify(to.rules)
      && JSON.stringify(cloneOverrides(state.featureColorOverrides)) === JSON.stringify(to.overrides)
      && readRef(state.results) === prepared.publishedResults
      && projectedResultsMatch(readRef(state.results), prepared.projection.nextResults)
    )
  });
  return command;
};
