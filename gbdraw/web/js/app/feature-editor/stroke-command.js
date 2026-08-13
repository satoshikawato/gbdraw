import {
  biologicalFeatureKey,
  catalogResultKey
} from '../../services/feature-catalog.js';
import { isFeatureInstanceHash } from '../../services/feature-instance-identity.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
import { styleFingerprint } from '../../services/style-revision.js';
import {
  buildStyleSnapshotCommand,
  exactJsonByteLength,
  immutableStyleCommandSnapshot
} from './style-snapshot-command.js';
import {
  featureStrokeStateFingerprint,
  stableFeatureStrokeStringify,
  stableFeatureTargetKey
} from './stroke-scope-plan.js';
import {
  normalizeFeatureStrokeColor,
  normalizeFeatureStrokeWidth,
  normalizeFeatureStrokeValue,
  resolveFeatureStrokeViewModel
} from './stroke-view-model.js';
import { prepareFeatureStrokeResultProjection } from './stroke-result-projection.js';

const text = (value) => String(value ?? '').trim();
const hasOwn = (value, key) => Object.prototype.hasOwnProperty.call(value || {}, key);
const cloneObject = (value) => cloneJson(value || {});

const readRef = (value) => (
  value && typeof value === 'object' && 'value' in value ? value.value : value
);

const writeRef = (target, value) => {
  if (target && typeof target === 'object' && 'value' in target) {
    target.value = value;
    return;
  }
  throw new Error('Feature stroke command requires a writable ref.');
};

const isWritableRef = (value) => Boolean(
  value && typeof value === 'object' && 'value' in value && value.__v_isReadonly !== true
);

const replaceObject = (target, source) => {
  Object.keys(target || {}).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = cloneJson(value);
  });
};

const restoreObjectReferences = (target, source) => {
  Object.keys(target || {}).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = value;
  });
};

const semanticStrokeProperties = (override) => Boolean(
  normalizeFeatureStrokeColor(override?.strokeColor)
  || normalizeFeatureStrokeWidth(override?.strokeWidth) !== null
);

const validateStrokeOverrideMap = (map, label) => {
  Object.entries(map || {}).forEach(([key, override]) => {
    if (!override || typeof override !== 'object' || Array.isArray(override)) {
      throw new Error(`${label} has an invalid entry: ${key}`);
    }
    if (
      (
        hasOwn(override, 'strokeColor')
        && text(override.strokeColor)
        && !normalizeFeatureStrokeColor(override.strokeColor)
      )
      || (
        hasOwn(override, 'strokeWidth')
        && override.strokeWidth !== null
        && override.strokeWidth !== undefined
        && override.strokeWidth !== ''
        && normalizeFeatureStrokeWidth(override.strokeWidth) === null
      )
    ) {
      throw new Error(`${label} has an invalid stroke value: ${key}`);
    }
  });
};

const removeStrokeProperties = (map, key, properties = ['strokeColor', 'strokeWidth']) => {
  const current = map[key];
  if (!current || typeof current !== 'object') return;
  properties.forEach((property) => delete current[property]);
  if (!semanticStrokeProperties(current)) delete map[key];
};

const applyStrokeProperties = (map, key, value) => {
  const current = map[key] && typeof map[key] === 'object'
    ? { ...map[key] }
    : {};
  if (hasOwn(value, 'strokeColor')) current.strokeColor = value.strokeColor;
  if (hasOwn(value, 'strokeWidth')) current.strokeWidth = value.strokeWidth;
  map[key] = current;
};

const catalogBindings = (catalog) => {
  const byResultKey = new Map();
  const byTargetKey = new Map();
  (Array.isArray(catalog?.items) ? catalog.items : []).forEach((item, resultIndex) => {
    const resultKey = catalogResultKey(item);
    if (!resultKey || byResultKey.has(resultKey)) {
      throw new Error('Feature stroke catalogue has missing or duplicate Result identity.');
    }
    const renderedByPair = new Map();
    (Array.isArray(item?.features) ? item.features : []).forEach((rendered) => {
      const pairKey = biologicalFeatureKey(rendered?.recordKey, rendered?.biologicalFeatureId);
      if (!pairKey) return;
      if (!renderedByPair.has(pairKey)) renderedByPair.set(pairKey, []);
      renderedByPair.get(pairKey).push(rendered);
    });
    const featureByPair = new Map();
    (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).forEach((feature) => {
      const pairKey = biologicalFeatureKey(feature?.recordKey, feature?.biologicalFeatureId);
      const targetKey = stableFeatureTargetKey({ ...feature, resultKey });
      const rendered = renderedByPair.get(pairKey) || [];
      if (
        !pairKey
        || !targetKey
        || featureByPair.has(pairKey)
        || byTargetKey.has(targetKey)
      ) {
        throw new Error('Feature stroke catalogue has ambiguous target identity.');
      }
      const binding = {
        resultKey,
        resultIndex,
        item,
        feature,
        pairKey,
        targetKey,
        rendered
      };
      featureByPair.set(pairKey, binding);
      byTargetKey.set(targetKey, binding);
    });
    byResultKey.set(resultKey, { resultKey, resultIndex, item, featureByPair });
  });
  return { byResultKey, byTargetKey };
};

const validatePlanTargets = (plan, bindings) => {
  if (!plan || plan.status !== 'ready') {
    throw new Error('A resolved ready Feature stroke plan is required.');
  }
  if (!['legend-caption', 'selected-features', 'single'].includes(text(plan.semanticScope))) {
    throw new Error('Feature stroke plan has an unsupported semantic scope.');
  }
  const affected = [...new Set((plan.affectedResultKeys || []).map(text).filter(Boolean))];
  if (affected.length === 0 || affected.length !== (plan.affectedResultKeys || []).length) {
    throw new Error('Feature stroke plan has invalid affected Result identities.');
  }
  const groups = Array.isArray(plan.targetsByResult) ? plan.targetsByResult : [];
  const plannedTargets = new Set();
  const plannedResults = new Set();
  groups.forEach((group) => {
    const resultKey = text(group?.resultKey);
    if (
      !affected.includes(resultKey)
      || !bindings.byResultKey.has(resultKey)
      || plannedResults.has(resultKey)
    ) {
      throw new Error(`Feature stroke target Result is unavailable: ${resultKey}`);
    }
    plannedResults.add(resultKey);
    const featureKeys = Array.isArray(group?.featureKeys) ? group.featureKeys : [];
    if (featureKeys.length === 0) {
      throw new Error(`Feature stroke target Result has no targets: ${resultKey}`);
    }
    featureKeys.forEach((rawFeatureKey) => {
      const featureKey = text(rawFeatureKey);
      const binding = bindings.byTargetKey.get(featureKey);
      if (!binding || binding.resultKey !== resultKey || plannedTargets.has(featureKey)) {
        throw new Error(`Feature stroke target is stale or cross-Result: ${featureKey}`);
      }
      if (
        ['selected-features', 'single'].includes(text(plan.semanticScope))
        && !isFeatureInstanceHash(binding.feature?.instanceHash)
      ) {
        throw new Error(`Feature stroke target has no valid exact identity: ${featureKey}`);
      }
      plannedTargets.add(featureKey);
    });
  });
  if (
    plannedTargets.size === 0
    || plannedResults.size !== affected.length
    || affected.some((resultKey) => !plannedResults.has(resultKey))
  ) {
    throw new Error('Feature stroke plan does not cover every affected Result exactly once.');
  }
  return { affected, targetKeys: [...plannedTargets].sort() };
};

const legendCaption = (feature) => text(
  feature?.effectiveLegendCaption
  ?? feature?.effective_legend_caption
  ?? feature?.legendCaption
  ?? feature?.legend_caption
  ?? feature?.caption
  ?? feature?.type
);

const strokeForBinding = ({
  binding,
  plannedCaption = '',
  featureStrokeOverrides,
  legendStrokeOverrides,
  originalSvgStroke
}) => {
  const caption = text(plannedCaption) || legendCaption(binding.feature);
  const model = resolveFeatureStrokeViewModel({
    exactOverride: featureStrokeOverrides?.[binding.pairKey] || null,
    captionOverride: caption ? legendStrokeOverrides?.[caption] || null : null,
    legendCaption: caption,
    svgDefaultStroke: originalSvgStroke
  });
  return Object.freeze({
    strokeColor: model.effectiveColor,
    strokeWidth: model.effectiveWidth,
    caption
  });
};

const buildNextMaps = ({ plan, targets, bindings, beforeFeature, beforeLegend }) => {
  const value = normalizeFeatureStrokeValue(plan?.semanticAfter?.value ?? plan?.intent?.value);
  if (!value) throw new Error('Feature stroke plan has an invalid value.');
  const nextFeature = cloneObject(beforeFeature);
  const nextLegend = cloneObject(beforeLegend);
  const isInherit = value.kind === 'inherit';
  const properties = isInherit
    ? ['strokeColor', 'strokeWidth']
    : [
        ...(hasOwn(value, 'strokeColor') ? ['strokeColor'] : []),
        ...(hasOwn(value, 'strokeWidth') ? ['strokeWidth'] : [])
      ];

  if (plan.semanticScope === 'legend-caption') {
    const caption = text(plan?.semanticAfter?.durableStrokeIntent?.caption);
    if (!caption) throw new Error('Feature stroke legend-caption scope has no caption.');
    targets.targetKeys.forEach((targetKey) => {
      removeStrokeProperties(nextFeature, bindings.byTargetKey.get(targetKey).pairKey, properties);
    });
    if (isInherit) removeStrokeProperties(nextLegend, caption, properties);
    else applyStrokeProperties(nextLegend, caption, value);
    return { featureStrokeOverrides: nextFeature, legendStrokeOverrides: nextLegend };
  }

  targets.targetKeys.forEach((targetKey) => {
    const key = bindings.byTargetKey.get(targetKey).pairKey;
    if (isInherit) removeStrokeProperties(nextFeature, key, properties);
    else applyStrokeProperties(nextFeature, key, value);
  });
  return { featureStrokeOverrides: nextFeature, legendStrokeOverrides: nextLegend };
};

const plannedCaptions = (plan, targetKeys) => {
  const intent = plan?.semanticAfter?.durableStrokeIntent || {};
  if (intent.kind === 'legend-caption') {
    return Object.fromEntries(targetKeys.map((key) => [key, text(intent.caption)]));
  }
  return Object.fromEntries(
    (Array.isArray(intent.selectors) ? intent.selectors : [])
      .map((selector) => [text(selector?.featureKey), text(selector?.legendCaption)])
      .filter(([key]) => key)
  );
};

const effectiveTargets = ({ plan, targets, bindings, maps, originalSvgStroke }) => {
  const captions = plannedCaptions(plan, targets.targetKeys);
  return Object.fromEntries(
  targets.targetKeys.map((targetKey) => [
    targetKey,
    strokeForBinding({
      binding: bindings.byTargetKey.get(targetKey),
      plannedCaption: captions[targetKey],
      featureStrokeOverrides: maps.featureStrokeOverrides,
      legendStrokeOverrides: maps.legendStrokeOverrides,
      originalSvgStroke
    })
  ])
  );
};

export const deriveFeatureStrokeTargets = ({
  plan,
  catalog,
  featureStrokeOverrides = {},
  legendStrokeOverrides = {},
  originalSvgStroke = {}
} = {}) => {
  const bindings = catalogBindings(catalog);
  const targets = validatePlanTargets(plan, bindings);
  return immutableStyleCommandSnapshot(effectiveTargets({
    plan,
    targets,
    bindings,
    maps: { featureStrokeOverrides, legendStrokeOverrides },
    originalSvgStroke
  }));
};

export const buildFeatureStrokeSemanticSnapshots = ({
  plan,
  catalog,
  featureStrokeOverrides = {},
  legendStrokeOverrides = {},
  originalSvgStroke = {}
} = {}) => {
  validateStrokeOverrideMap(featureStrokeOverrides, 'Feature stroke overrides');
  validateStrokeOverrideMap(legendStrokeOverrides, 'Legend stroke overrides');
  const bindings = catalogBindings(catalog);
  const targets = validatePlanTargets(plan, bindings);
  const beforeMaps = {
    featureStrokeOverrides: cloneObject(featureStrokeOverrides),
    legendStrokeOverrides: cloneObject(legendStrokeOverrides)
  };
  const afterMaps = buildNextMaps({
    plan,
    targets,
    bindings,
    beforeFeature: beforeMaps.featureStrokeOverrides,
    beforeLegend: beforeMaps.legendStrokeOverrides
  });
  const beforeStrokes = effectiveTargets({
    plan,
    targets,
    bindings,
    maps: beforeMaps,
    originalSvgStroke
  });
  const afterStrokes = effectiveTargets({
    plan,
    targets,
    bindings,
    maps: afterMaps,
    originalSvgStroke
  });
  const semanticSide = (maps) => immutableStyleCommandSnapshot({
    ...maps,
    strokeFingerprint: featureStrokeStateFingerprint({ ...maps, originalSvgStroke })
  });
  return immutableStyleCommandSnapshot({
    before: semanticSide(beforeMaps),
    after: semanticSide(afterMaps),
    affectedResultKeys: targets.affected,
    targetFeatureKeys: targets.targetKeys,
    strokesBeforeByTargetKey: beforeStrokes,
    strokesAfterByTargetKey: afterStrokes,
    legendCaption: plan.semanticScope === 'legend-caption'
      ? text(plan?.semanticAfter?.durableStrokeIntent?.caption)
      : '',
    counters: {
      targetFeatures: targets.targetKeys.length,
      affectedResults: targets.affected.length
    }
  });
};

const currentStyleFingerprint = (state) => styleFingerprint({
  rules: state?.manualSpecificRules,
  appliedPaletteName: readRef(state?.appliedPaletteName),
  appliedPaletteColors: readRef(state?.appliedPaletteColors)
});

const currentStrokeFingerprint = (state) => featureStrokeStateFingerprint({
  featureStrokeOverrides: state?.featureStrokeOverrides,
  legendStrokeOverrides: state?.legendStrokeOverrides,
  originalSvgStroke: readRef(state?.originalSvgStroke) || {}
});

const strokeMapsMatch = (state, side) => (
  stableFeatureStrokeStringify(state?.featureStrokeOverrides || {})
    === stableFeatureStrokeStringify(side?.featureStrokeOverrides || {})
  && stableFeatureStrokeStringify(state?.legendStrokeOverrides || {})
    === stableFeatureStrokeStringify(side?.legendStrokeOverrides || {})
);

const resultContent = (result) => String(result?.content ?? '');
const resultSignature = (result) => JSON.stringify({
  name: text(result?.name),
  content: resultContent(result)
});

const resultTopologyIsCurrent = (results, catalog, affectedResultKeys) => {
  if (!Array.isArray(results) || !Array.isArray(catalog?.items) || results.length !== catalog.items.length) {
    return false;
  }
  const resultKeys = new Set();
  for (let index = 0; index < catalog.items.length; index += 1) {
    const item = catalog.items[index];
    const resultKey = catalogResultKey(item);
    if (
      !resultKey
      || resultKeys.has(resultKey)
      || !results[index]
      || (item?.resultIndex !== undefined && Number(item.resultIndex) !== index)
      || (text(item?.resultName) && text(item.resultName) !== text(results[index]?.name))
    ) return false;
    resultKeys.add(resultKey);
  }
  return affectedResultKeys.every((key) => resultKeys.has(key));
};

const targetIdentity = (binding) => JSON.stringify({
  instanceHash: text(binding?.feature?.instanceHash),
  caption: legendCaption(binding?.feature),
  renderedSvgIds: binding?.rendered.map((entry) => text(entry?.svgId ?? entry?.svg_id)).sort()
});

const targetsAreCurrent = (plan, catalog, expectedIdentity) => {
  try {
    const bindings = catalogBindings(catalog);
    const targets = validatePlanTargets(plan, bindings);
    return targets.targetKeys.every(
      (targetKey) => expectedIdentity[targetKey] === targetIdentity(bindings.byTargetKey.get(targetKey))
    );
  } catch {
    return false;
  }
};

const mountedContext = (getMountedContext) => {
  const raw = typeof getMountedContext === 'function' ? getMountedContext() : null;
  if (!raw) return { resultIndex: null, resultKey: '', svg: null };
  const index = Number(raw.resultIndex);
  return {
    ...raw,
    resultIndex: Number.isInteger(index) && index >= 0 ? index : null,
    resultKey: text(raw.resultKey),
    svg: raw.svg || null
  };
};

const mountIsCurrent = (mounted, catalog) => {
  if (mounted.resultIndex === null) return !mounted.resultKey && !mounted.svg;
  return mounted.resultIndex < (catalog?.items?.length || 0)
    && catalogResultKey(catalog.items[mounted.resultIndex]) === mounted.resultKey;
};

const sameMount = (left, right) => (
  left.resultIndex === right.resultIndex
  && left.resultKey === right.resultKey
  && left.svg === right.svg
);

const projectedResultsMatch = (published, projected) => (
  Array.isArray(published)
  && Array.isArray(projected)
  && published.length === projected.length
  && published.every((result, index) => resultSignature(result) === resultSignature(projected[index]))
);

/** Build one History-compatible, all-Result Feature stroke command. */
export const buildFeatureStrokeCommand = async ({
  plan,
  state,
  catalog,
  prepareResultProjection = prepareFeatureStrokeResultProjection,
  getMountedContext = null,
  commitMountedProjection = null,
  restoreMountedProjection = null,
  captureRuntimeState = null,
  restoreRuntimeState = null,
  reconcile = null,
  refreshPresentation = null,
  label = 'Change feature stroke'
} = {}) => {
  if (
    !state
    || !state.featureStrokeOverrides
    || !state.legendStrokeOverrides
    || typeof prepareResultProjection !== 'function'
  ) {
    throw new Error('Feature stroke command requires state and a Result projection owner.');
  }
  const expectedEpoch = Number(plan?.intent?.documentEpoch ?? 0);
  const expectedGeneration = Number(plan?.intent?.resultGenerationKey ?? 0);
  const expectedStyleFingerprint = text(plan?.semanticBefore?.styleFingerprint);
  let expectedRevision = Number(plan?.semanticBefore?.styleRevision ?? plan?.intent?.styleRevision ?? 0);
  const originalSvgStroke = cloneObject(readRef(state.originalSvgStroke) || {});
  if (
    Number(readRef(state.documentEpoch) || 0) !== expectedEpoch
    || Number(readRef(state.resultGenerationKey) || 0) !== expectedGeneration
    || Number(readRef(state.semanticStyleRevision) || 0) !== expectedRevision
    || (expectedStyleFingerprint && currentStyleFingerprint(state) !== expectedStyleFingerprint)
    || (readRef(state.featureCatalog) && readRef(state.featureCatalog) !== catalog)
  ) {
    throw new Error('The Feature stroke plan became stale before command preparation.');
  }
  const semantic = buildFeatureStrokeSemanticSnapshots({
    plan,
    catalog,
    featureStrokeOverrides: state.featureStrokeOverrides,
    legendStrokeOverrides: state.legendStrokeOverrides,
    originalSvgStroke
  });
  if (
    text(plan?.semanticBefore?.strokeStateFingerprint)
    && text(plan.semanticBefore.strokeStateFingerprint) !== semantic.before.strokeFingerprint
  ) {
    throw new Error('The Feature stroke state changed before command preparation.');
  }
  const initialBindings = catalogBindings(catalog);
  const expectedTargetIdentity = Object.fromEntries(
    semantic.targetFeatureKeys.map((key) => [key, targetIdentity(initialBindings.byTargetKey.get(key))])
  );
  const resultIndexByKey = new Map(
    catalog.items.map((item, index) => [catalogResultKey(item), index])
  );
  let expectedResults = readRef(state.results);
  let expectedResultSignatures = Object.fromEntries(
    semantic.affectedResultKeys.map((key) => [key, resultSignature(expectedResults?.[resultIndexByKey.get(key)])])
  );

  const validateCurrent = ({ from, prepared = null }) => {
    const results = readRef(state.results);
    const mounted = mountedContext(getMountedContext);
    if (
      Number(readRef(state.documentEpoch) || 0) !== expectedEpoch
      || Number(readRef(state.resultGenerationKey) || 0) !== expectedGeneration
      || Number(readRef(state.semanticStyleRevision) || 0) !== expectedRevision
      || (expectedStyleFingerprint && currentStyleFingerprint(state) !== expectedStyleFingerprint)
      || currentStrokeFingerprint(state) !== from.strokeFingerprint
      || !strokeMapsMatch(state, from)
      || results !== expectedResults
      || !resultTopologyIsCurrent(results, catalog, semantic.affectedResultKeys)
      || semantic.affectedResultKeys.some((key) => (
        resultSignature(results?.[resultIndexByKey.get(key)]) !== expectedResultSignatures[key]
      ))
      || !targetsAreCurrent(plan, catalog, expectedTargetIdentity)
      || (readRef(state.featureCatalog) && readRef(state.featureCatalog) !== catalog)
      || !mountIsCurrent(mounted, catalog)
      || (mounted.resultKey && semantic.affectedResultKeys.includes(mounted.resultKey) && !mounted.svg)
    ) return false;
    if (prepared && (
      results !== prepared.sourceResults
      || !sameMount(mounted, prepared.preflightMounted)
      || prepared.sourceResultRecords.some(({ index, result, signature }) => (
        results[index] !== result || resultSignature(result) !== signature
      ))
    )) return false;
    return true;
  };

  const prepareTransition = async ({ from, to, direction }) => {
    const mounted = mountedContext(getMountedContext);
    const sourceResults = readRef(state.results);
    const sourceResultRecords = semantic.affectedResultKeys.map((resultKey) => {
      const index = resultIndexByKey.get(resultKey);
      const result = sourceResults[index];
      return { index, result, signature: resultSignature(result) };
    });
    const strokesByTargetKey = direction === 'undo'
      ? semantic.strokesBeforeByTargetKey
      : semantic.strokesAfterByTargetKey;
    const targetFeatureKeysByResult = new Map(
      plan.targetsByResult.map((entry) => [entry.resultKey, cloneJson(entry.featureKeys)])
    );
    let legendStroke = null;
    if (semantic.legendCaption) {
      const firstTarget = semantic.targetFeatureKeys[0];
      const binding = initialBindings.byTargetKey.get(firstTarget);
      legendStroke = strokeForBinding({
        binding,
        featureStrokeOverrides: {},
        legendStrokeOverrides: to.legendStrokeOverrides,
        originalSvgStroke
      });
    }
    const projection = await prepareResultProjection({
      direction,
      from,
      to,
      plan,
      catalog,
      results: sourceResults,
      strokesByTargetKey,
      affectedResultKeys: semantic.affectedResultKeys,
      targetFeatureKeysByResult,
      legendCaption: semantic.legendCaption,
      legendStroke,
      mounted
    });
    if (!projection || !Array.isArray(projection.nextResults)) {
      throw new Error('Feature stroke Result projection did not prepare a complete Result array.');
    }
    if (readRef(state.results) !== sourceResults) {
      throw new Error('Feature stroke Results changed during projection preparation.');
    }
    if (projection.nextResults.length !== sourceResults.length) {
      throw new Error('Feature stroke Result projection changed Result topology.');
    }
    const affectedIndexes = new Set(
      semantic.affectedResultKeys.map((key) => resultIndexByKey.get(key))
    );
    projection.nextResults.forEach((result, index) => {
      if (!result || text(result.name) !== text(sourceResults[index]?.name)) {
        throw new Error(`Feature stroke Result identity changed at index ${index}.`);
      }
      if (!affectedIndexes.has(index) && result !== sourceResults[index]) {
        throw new Error(`Feature stroke projection replaced unaffected Result ${index}.`);
      }
    });
    const retainedForHistory = {
      affectedResults: semantic.affectedResultKeys.map((resultKey) => {
        const index = resultIndexByKey.get(resultKey);
        return {
          resultKey,
          before: sourceResults[index],
          after: projection.nextResults[index],
          admission: projection.admissionMetadataByResultKey?.[resultKey] || null
        };
      })
    };
    return {
      projection,
      sourceResults,
      sourceResultRecords,
      preflightMounted: mounted,
      retainedForHistory,
      retainedBytes: exactJsonByteLength(retainedForHistory),
      counters: {
        projectedResults: Number(projection.counters?.affectedResults || semantic.affectedResultKeys.length),
        mountedSerializations: Number(projection.counters?.mountedSerializations || 0),
        detachedPasses: Number(projection.counters?.detachedPasses || 0),
        changedResults: Number(projection.counters?.changedResults || 0),
        resultArraySwaps: 1
      }
    };
  };

  const captureExactState = () => {
    const runtime = typeof captureRuntimeState === 'function' ? captureRuntimeState() : null;
    if (runtime && typeof runtime.then === 'function') {
      throw new Error('Feature stroke runtime capture must be synchronous.');
    }
    return {
      featureOverrides: { ...state.featureStrokeOverrides },
      legendOverrides: { ...state.legendStrokeOverrides },
      results: readRef(state.results),
      revision: Number(readRef(state.semanticStyleRevision) || 0),
      fingerprint: text(readRef(state.semanticStyleFingerprint)),
      expectedRevision,
      expectedResults,
      expectedResultSignatures: { ...expectedResultSignatures },
      presentation: {
        ...(isWritableRef(state.clickedFeature) ? { clickedFeature: readRef(state.clickedFeature) } : {}),
        ...(isWritableRef(state.selectedFeatures) ? { selectedFeatures: readRef(state.selectedFeatures) } : {}),
        ...(isWritableRef(state.legendEntries) ? { legendEntries: readRef(state.legendEntries) } : {}),
        ...(isWritableRef(state.selectedResultIndex)
          ? { selectedResultIndex: readRef(state.selectedResultIndex) }
          : {})
      },
      runtime
    };
  };

  const restoreExactState = (snapshot, context) => {
    if (typeof restoreMountedProjection === 'function') {
      const restored = restoreMountedProjection({ snapshot, context, plan });
      if (restored && typeof restored.then === 'function') {
        throw new Error('Mounted Feature stroke rollback must be synchronous.');
      }
      if (restored === false) return false;
    }
    restoreObjectReferences(state.featureStrokeOverrides, snapshot.featureOverrides);
    restoreObjectReferences(state.legendStrokeOverrides, snapshot.legendOverrides);
    writeRef(state.results, snapshot.results);
    writeRef(state.semanticStyleRevision, snapshot.revision);
    if (isWritableRef(state.semanticStyleFingerprint)) {
      writeRef(state.semanticStyleFingerprint, snapshot.fingerprint);
    }
    Object.entries(snapshot.presentation).forEach(([key, value]) => {
      if (isWritableRef(state[key])) writeRef(state[key], value);
    });
    expectedRevision = snapshot.expectedRevision;
    expectedResults = snapshot.expectedResults;
    expectedResultSignatures = { ...snapshot.expectedResultSignatures };
    if (typeof restoreRuntimeState === 'function') {
      const restored = restoreRuntimeState(snapshot.runtime);
      if (restored && typeof restored.then === 'function') {
        throw new Error('Feature stroke runtime rollback must be synchronous.');
      }
      if (restored === false) return false;
    }
    return true;
  };

  const commitTransition = ({ to, direction, prepared }) => {
    replaceObject(state.featureStrokeOverrides, to.featureStrokeOverrides);
    replaceObject(state.legendStrokeOverrides, to.legendStrokeOverrides);
    if (semantic.affectedResultKeys.includes(prepared.preflightMounted.resultKey)) {
      if (typeof commitMountedProjection !== 'function' || typeof restoreMountedProjection !== 'function') {
        throw new Error('Affected mounted Feature stroke projection has no atomic commit/restore owner.');
      }
      const committed = commitMountedProjection({
        direction,
        prepared: prepared.projection,
        plan
      });
      if (committed && typeof committed.then === 'function') {
        throw new Error('Mounted Feature stroke commit must be synchronous.');
      }
      if (committed === false) return false;
    }
    writeRef(state.results, prepared.projection.nextResults);
    const published = readRef(state.results);
    if (!projectedResultsMatch(published, prepared.projection.nextResults)) {
      throw new Error('Feature stroke Result publication did not preserve the staged projection.');
    }
    const reconciled = typeof reconcile === 'function'
      ? reconcile({ direction, plan, prepared, to })
      : undefined;
    if (reconciled && typeof reconciled.then === 'function') {
      throw new Error('Feature stroke reconciliation must be synchronous.');
    }
    if (reconciled === false) return false;
    const refreshed = typeof refreshPresentation === 'function'
      ? refreshPresentation({ direction, plan, prepared, to })
      : undefined;
    if (refreshed && typeof refreshed.then === 'function') {
      throw new Error('Feature stroke presentation refresh must be synchronous.');
    }
    if (refreshed === false) return false;
    expectedRevision += 1;
    writeRef(state.semanticStyleRevision, expectedRevision);
    expectedResults = published;
    expectedResultSignatures = Object.fromEntries(
      semantic.affectedResultKeys.map((key) => [key, resultSignature(published[resultIndexByKey.get(key)])])
    );
    prepared.publishedResults = published;
    return true;
  };

  return buildStyleSnapshotCommand({
    label,
    before: semantic.before,
    after: semantic.after,
    metadata: {
      planToken: text(plan.token),
      documentEpoch: expectedEpoch,
      resultGenerationKey: expectedGeneration,
      semanticScope: text(plan.semanticScope),
      resultExtent: text(plan.resultExtent),
      affectedResultKeys: semantic.affectedResultKeys,
      targetFeatureKeys: semantic.targetFeatureKeys
    },
    counters: semantic.counters,
    isNoop: (before, after) => before.strokeFingerprint === after.strokeFingerprint,
    validateCurrent,
    prepareTransition,
    commitTransition,
    captureExactState,
    restoreExactState,
    assertCommitted: ({ to, prepared }) => (
      currentStrokeFingerprint(state) === to.strokeFingerprint
      && strokeMapsMatch(state, to)
      && Number(readRef(state.semanticStyleRevision) || 0) === expectedRevision
      && readRef(state.results) === prepared.publishedResults
      && projectedResultsMatch(readRef(state.results), prepared.projection.nextResults)
    )
  });
};
