import { biologicalFeatureKey, catalogResultKey } from '../../services/feature-catalog.js';
import { isFeatureInstanceHash } from '../../services/feature-instance-identity.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
import { styleFingerprint } from '../../services/style-revision.js';
import {
  buildStyleSnapshotCommand,
  exactJsonByteLength,
  immutableStyleCommandSnapshot
} from './style-snapshot-command.js';
import {
  featureLabelStateFingerprint,
  stableFeatureLabelStringify,
  stableFeatureTargetKey
} from './label-scope-plan.js';
import { resolveFeatureLabelViewModel } from './label-view-model.js';
import { prepareFeatureLabelResultProjection } from './label-result-projection.js';

const text = (value) => String(value ?? '').trim();
const readRef = (value) => (
  value && typeof value === 'object' && 'value' in value ? value.value : value
);
const writeRef = (target, value) => {
  if (!target || typeof target !== 'object' || !('value' in target)) {
    throw new Error('Feature label command requires a writable ref.');
  }
  target.value = value;
};
const isWritableRef = (value) => Boolean(
  value && typeof value === 'object' && 'value' in value && value.__v_isReadonly !== true
);
const cloneObject = (value) => cloneJson(value || {});

const replaceObject = (target, source) => {
  Object.keys(target || {}).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => { target[key] = cloneJson(value); });
};

const restoreObjectReferences = (target, source) => {
  Object.keys(target || {}).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => { target[key] = value; });
};

const catalogBindings = (catalog) => {
  const byResultKey = new Map();
  const byTargetKey = new Map();
  (Array.isArray(catalog?.items) ? catalog.items : []).forEach((item, resultIndex) => {
    const resultKey = catalogResultKey(item);
    if (!resultKey || byResultKey.has(resultKey)) {
      throw new Error('Feature label catalogue has missing or duplicate Result identity.');
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
        !pairKey || !targetKey
        || featureByPair.has(pairKey) || byTargetKey.has(targetKey)
      ) {
        throw new Error('Feature label catalogue has ambiguous target identity.');
      }
      const binding = { resultKey, resultIndex, item, feature, pairKey, targetKey, rendered };
      featureByPair.set(pairKey, binding);
      byTargetKey.set(targetKey, binding);
    });
    byResultKey.set(resultKey, { resultKey, resultIndex, item, featureByPair });
  });
  return { byResultKey, byTargetKey };
};

const selectorsByTargetKey = (plan) => new Map(
  (Array.isArray(plan?.semanticAfter?.durableLabelIntent?.selectors)
    ? plan.semanticAfter.durableLabelIntent.selectors
    : []).map((selector) => [text(selector?.featureKey), selector])
    .filter(([key]) => key)
);

const validatePlanTargets = (plan, bindings) => {
  if (!plan || plan.status !== 'ready') {
    throw new Error('A resolved ready Feature label plan is required.');
  }
  if (!['rendered-label', 'source-annotation-label', 'selected-features', 'single']
    .includes(text(plan.semanticScope))) {
    throw new Error('Feature label plan has an unsupported semantic scope.');
  }
  const affected = [...new Set((plan.affectedResultKeys || []).map(text).filter(Boolean))];
  if (affected.length === 0 || affected.length !== (plan.affectedResultKeys || []).length) {
    throw new Error('Feature label plan has invalid affected Result identities.');
  }
  const selectors = selectorsByTargetKey(plan);
  const targetKeys = new Set();
  const resultKeys = new Set();
  (Array.isArray(plan.targetsByResult) ? plan.targetsByResult : []).forEach((group) => {
    const resultKey = text(group?.resultKey);
    if (!affected.includes(resultKey) || !bindings.byResultKey.has(resultKey) || resultKeys.has(resultKey)) {
      throw new Error(`Feature label target Result is unavailable: ${resultKey}`);
    }
    resultKeys.add(resultKey);
    const keys = Array.isArray(group?.featureKeys) ? group.featureKeys : [];
    if (keys.length === 0) throw new Error(`Feature label target Result has no targets: ${resultKey}`);
    keys.forEach((rawKey) => {
      const key = text(rawKey);
      const binding = bindings.byTargetKey.get(key);
      const selector = selectors.get(key);
      if (!binding || binding.resultKey !== resultKey || targetKeys.has(key) || !selector) {
        throw new Error(`Feature label target is stale or cross-Result: ${key}`);
      }
      const catalogSvgIds = binding.rendered
        .map((entry) => text(entry?.svgId ?? entry?.svg_id)).filter(Boolean).sort();
      const selectorSvgIds = [...new Set((selector?.renderedSvgIds || []).map(text).filter(Boolean))].sort();
      if (stableFeatureLabelStringify(catalogSvgIds) !== stableFeatureLabelStringify(selectorSvgIds)) {
        throw new Error(`Feature label SVG identity changed: ${key}`);
      }
      if (
        ['selected-features', 'single'].includes(text(plan.semanticScope))
        && !isFeatureInstanceHash(binding.feature?.instanceHash)
      ) {
        throw new Error(`Feature label target has no valid exact identity: ${key}`);
      }
      targetKeys.add(key);
    });
  });
  if (
    targetKeys.size === 0 || resultKeys.size !== affected.length
    || affected.some((key) => !resultKeys.has(key)) || selectors.size !== targetKeys.size
  ) {
    throw new Error('Feature label plan does not cover every affected Result exactly once.');
  }
  return { affected, targetKeys: [...targetKeys].sort(), selectors };
};

const resolveLabelForTarget = ({ binding, selector, maps }) => {
  const renderedSvgIds = binding.rendered
    .map((entry) => text(entry?.svgId ?? entry?.svg_id)).filter(Boolean);
  const model = resolveFeatureLabelViewModel({
    feature: {
      ...binding.feature,
      sourceAnnotationLabel: String(selector?.sourceText ?? ''),
      sourceAnnotationQualifier: String(selector?.qualifier ?? '')
    },
    renderedSvgIds,
    featureOverrides: maps.featureOverrides,
    bulkOverrides: maps.bulkOverrides,
    featureOverrideSources: maps.featureOverrideSources,
    presentationText: String(selector?.fromText ?? ''),
    presentationSourceText: String(selector?.sourceText ?? '')
  });
  return Object.freeze({
    text: model.effectiveText,
    sourceText: model.sourceText,
    fromText: String(selector?.fromText ?? model.effectiveText),
    renderedSvgIds
  });
};

const effectiveLabels = ({ targets, bindings, maps }) => Object.fromEntries(
  targets.targetKeys.map((key) => [
    key,
    resolveLabelForTarget({
      binding: bindings.byTargetKey.get(key),
      selector: targets.selectors.get(key),
      maps
    })
  ])
);

const nextMaps = ({ plan, targets, before }) => {
  const next = {
    featureOverrides: cloneObject(before.featureOverrides),
    bulkOverrides: cloneObject(before.bulkOverrides),
    featureOverrideSources: cloneObject(before.featureOverrideSources)
  };
  const newText = plan?.semanticAfter?.newText;
  if (typeof newText !== 'string') throw new Error('Feature label plan has invalid label text.');
  if (plan.semanticScope === 'source-annotation-label') {
    const sourceText = String(plan?.semanticAfter?.durableLabelIntent?.sourceText ?? '');
    if (!sourceText) throw new Error('Feature label source scope has no source text.');
    next.bulkOverrides[sourceText] = newText;
    targets.targetKeys.forEach((key) => {
      const selector = targets.selectors.get(key);
      (selector?.renderedSvgIds || []).forEach((svgId) => {
        delete next.featureOverrides[svgId];
        delete next.featureOverrideSources[svgId];
      });
    });
    return next;
  }
  targets.targetKeys.forEach((key) => {
    const selector = targets.selectors.get(key);
    const sourceText = String(selector?.sourceText ?? '');
    (selector?.renderedSvgIds || []).forEach((svgId) => {
      next.featureOverrides[svgId] = newText;
      if (sourceText) next.featureOverrideSources[svgId] = sourceText;
      else delete next.featureOverrideSources[svgId];
    });
  });
  return next;
};

export const buildFeatureLabelSemanticSnapshots = ({
  plan,
  catalog,
  labelTextFeatureOverrides = {},
  labelTextBulkOverrides = {},
  labelTextFeatureOverrideSources = {}
} = {}) => {
  const bindings = catalogBindings(catalog);
  const targets = validatePlanTargets(plan, bindings);
  const beforeMaps = {
    featureOverrides: cloneObject(labelTextFeatureOverrides),
    bulkOverrides: cloneObject(labelTextBulkOverrides),
    featureOverrideSources: cloneObject(labelTextFeatureOverrideSources)
  };
  const afterMaps = nextMaps({ plan, targets, before: beforeMaps });
  const side = (maps) => immutableStyleCommandSnapshot({
    ...maps,
    labelFingerprint: featureLabelStateFingerprint({
      labelTextFeatureOverrides: maps.featureOverrides,
      labelTextBulkOverrides: maps.bulkOverrides,
      labelTextFeatureOverrideSources: maps.featureOverrideSources
    })
  });
  return immutableStyleCommandSnapshot({
    before: side(beforeMaps),
    after: side(afterMaps),
    affectedResultKeys: targets.affected,
    targetFeatureKeys: targets.targetKeys,
    labelsBeforeByTargetKey: effectiveLabels({ targets, bindings, maps: beforeMaps }),
    labelsAfterByTargetKey: effectiveLabels({ targets, bindings, maps: afterMaps }),
    projectionMatch: {
      semanticScope: text(plan.semanticScope),
      matchText: String(plan?.semanticAfter?.durableLabelIntent?.matchText ?? ''),
      sourceText: String(plan?.semanticAfter?.durableLabelIntent?.sourceText ?? '')
    },
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

const currentLabelFingerprint = (state) => featureLabelStateFingerprint({
  labelTextFeatureOverrides: state?.labelTextFeatureOverrides,
  labelTextBulkOverrides: state?.labelTextBulkOverrides,
  labelTextFeatureOverrideSources: state?.labelTextFeatureOverrideSources
});

const mapsMatch = (state, side) => (
  stableFeatureLabelStringify(state?.labelTextFeatureOverrides || {})
    === stableFeatureLabelStringify(side?.featureOverrides || {})
  && stableFeatureLabelStringify(state?.labelTextBulkOverrides || {})
    === stableFeatureLabelStringify(side?.bulkOverrides || {})
  && stableFeatureLabelStringify(state?.labelTextFeatureOverrideSources || {})
    === stableFeatureLabelStringify(side?.featureOverrideSources || {})
);

const resultSignature = (result) => JSON.stringify({
  name: text(result?.name), content: String(result?.content ?? '')
});

const topologyIsCurrent = (results, catalog, affected) => {
  if (!Array.isArray(results) || results.length !== (catalog?.items?.length || 0)) return false;
  const keys = new Set();
  for (let index = 0; index < catalog.items.length; index += 1) {
    const key = catalogResultKey(catalog.items[index]);
    if (!key || keys.has(key) || !results[index]
      || (text(catalog.items[index]?.resultName) && text(catalog.items[index].resultName) !== text(results[index]?.name))) {
      return false;
    }
    keys.add(key);
  }
  return affected.every((key) => keys.has(key));
};

const ledgerIsCurrent = (state, catalog, fingerprint) => {
  if (!state?.validatedStyleFingerprintByResultKey) return true;
  const ledger = readRef(state.validatedStyleFingerprintByResultKey) || {};
  return catalog.items.every((item) => text(ledger[catalogResultKey(item)]) === fingerprint);
};

const targetIdentity = (binding) => stableFeatureLabelStringify({
  instanceHash: text(binding?.feature?.instanceHash),
  renderedSvgIds: binding?.rendered.map((entry) => text(entry?.svgId ?? entry?.svg_id)).sort(),
  product: binding?.feature?.product,
  gene: binding?.feature?.gene,
  locusTag: binding?.feature?.locus_tag,
  qualifiers: binding?.feature?.qualifiers || null
});

const targetsAreCurrent = (plan, catalog, expected) => {
  try {
    const bindings = catalogBindings(catalog);
    const targets = validatePlanTargets(plan, bindings);
    return targets.targetKeys.every(
      (key) => expected[key] === targetIdentity(bindings.byTargetKey.get(key))
    );
  } catch {
    return false;
  }
};

const mountedContext = (getter) => {
  const raw = typeof getter === 'function' ? getter() : null;
  if (!raw) return { resultIndex: null, resultKey: '', svg: null };
  const index = Number(raw.resultIndex);
  return {
    ...raw,
    resultIndex: Number.isInteger(index) && index >= 0 ? index : null,
    resultKey: text(raw.resultKey),
    svg: raw.svg || null
  };
};

const mountIsCurrent = (mounted, catalog) => (
  mounted.resultIndex === null
    ? !mounted.resultKey && !mounted.svg
    : mounted.resultIndex < (catalog?.items?.length || 0)
      && catalogResultKey(catalog.items[mounted.resultIndex]) === mounted.resultKey
);

const sameMount = (left, right) => (
  left.resultIndex === right.resultIndex && left.resultKey === right.resultKey && left.svg === right.svg
);

const projectedResultsMatch = (left, right) => (
  Array.isArray(left) && Array.isArray(right) && left.length === right.length
  && left.every((result, index) => resultSignature(result) === resultSignature(right[index]))
);

/** Build one all-before-mutate, History-compatible Feature label command. */
export const buildFeatureLabelCommand = async ({
  plan,
  state,
  catalog,
  prepareResultProjection = prepareFeatureLabelResultProjection,
  prepareLabelReflow = null,
  getMountedContext = null,
  commitMountedProjection = null,
  restoreMountedProjection = null,
  captureRuntimeState = null,
  restoreRuntimeState = null,
  reconcile = null,
  refreshPresentation = null,
  label = 'Change feature label'
} = {}) => {
  if (
    !state || !state.labelTextFeatureOverrides || !state.labelTextBulkOverrides
    || !state.labelTextFeatureOverrideSources || typeof prepareResultProjection !== 'function'
  ) {
    throw new Error('Feature label command requires state and a Result projection owner.');
  }
  const expectedEpoch = Number(plan?.intent?.documentEpoch ?? 0);
  const expectedGeneration = Number(plan?.intent?.resultGenerationKey ?? 0);
  const expectedStyleFingerprint = text(plan?.semanticBefore?.styleFingerprint);
  let expectedRevision = Number(plan?.semanticBefore?.styleRevision ?? 0);
  if (
    Number(readRef(state.documentEpoch) || 0) !== expectedEpoch
    || Number(readRef(state.resultGenerationKey) || 0) !== expectedGeneration
    || Number(readRef(state.semanticStyleRevision) || 0) !== expectedRevision
    || (expectedStyleFingerprint && currentStyleFingerprint(state) !== expectedStyleFingerprint)
    || !ledgerIsCurrent(state, catalog, expectedStyleFingerprint)
    || (readRef(state.featureCatalog) && readRef(state.featureCatalog) !== catalog)
  ) {
    throw new Error('The Feature label plan became stale before command preparation.');
  }
  const semantic = buildFeatureLabelSemanticSnapshots({
    plan,
    catalog,
    labelTextFeatureOverrides: state.labelTextFeatureOverrides,
    labelTextBulkOverrides: state.labelTextBulkOverrides,
    labelTextFeatureOverrideSources: state.labelTextFeatureOverrideSources
  });
  if (
    text(plan?.semanticBefore?.labelStateFingerprint)
    && text(plan.semanticBefore.labelStateFingerprint) !== semantic.before.labelFingerprint
  ) {
    throw new Error('The Feature label state changed before command preparation.');
  }
  const bindings = catalogBindings(catalog);
  const expectedTargetIdentity = Object.fromEntries(
    semantic.targetFeatureKeys.map((key) => [key, targetIdentity(bindings.byTargetKey.get(key))])
  );
  const resultIndexByKey = new Map(catalog.items.map((item, index) => [catalogResultKey(item), index]));
  let expectedResults = readRef(state.results);
  let expectedResultSignatures = Object.fromEntries(
    semantic.affectedResultKeys.map((key) => [key, resultSignature(expectedResults?.[resultIndexByKey.get(key)])])
  );
  const expectedLedger = readRef(state.validatedStyleFingerprintByResultKey);
  const expectedLedgerSignature = stableFeatureLabelStringify(expectedLedger || {});

  const validateCurrent = ({ from, prepared = null }) => {
    const results = readRef(state.results);
    const mounted = mountedContext(getMountedContext);
    if (
      Number(readRef(state.documentEpoch) || 0) !== expectedEpoch
      || Number(readRef(state.resultGenerationKey) || 0) !== expectedGeneration
      || Number(readRef(state.semanticStyleRevision) || 0) !== expectedRevision
      || (expectedStyleFingerprint && currentStyleFingerprint(state) !== expectedStyleFingerprint)
      || currentLabelFingerprint(state) !== from.labelFingerprint
      || !mapsMatch(state, from)
      || results !== expectedResults
      || !topologyIsCurrent(results, catalog, semantic.affectedResultKeys)
      || semantic.affectedResultKeys.some((key) => (
        resultSignature(results?.[resultIndexByKey.get(key)]) !== expectedResultSignatures[key]
      ))
      || !targetsAreCurrent(plan, catalog, expectedTargetIdentity)
      || (readRef(state.featureCatalog) && readRef(state.featureCatalog) !== catalog)
      || !ledgerIsCurrent(state, catalog, expectedStyleFingerprint)
      || readRef(state.validatedStyleFingerprintByResultKey) !== expectedLedger
      || stableFeatureLabelStringify(readRef(state.validatedStyleFingerprintByResultKey) || {})
        !== expectedLedgerSignature
      || !mountIsCurrent(mounted, catalog)
      || (mounted.resultKey && semantic.affectedResultKeys.includes(mounted.resultKey) && !mounted.svg)
    ) return false;
    if (prepared && (
      results !== prepared.sourceResults || !sameMount(mounted, prepared.preflightMounted)
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
    const labelsByTargetKey = direction === 'undo'
      ? semantic.labelsBeforeByTargetKey
      : semantic.labelsAfterByTargetKey;
    const reflow = typeof prepareLabelReflow === 'function'
      ? await prepareLabelReflow({
          direction, from, to, plan, catalog, results: sourceResults,
          labelsByTargetKey,
          affectedResultKeys: semantic.affectedResultKeys,
          mounted
        })
      : null;
    if (readRef(state.results) !== sourceResults) {
      throw new Error('Feature label state changed during reflow preparation.');
    }
    const preparedSvgByResultKey = reflow?.preparedSvgByResultKey || null;
    if (reflow) {
      for (const resultKey of semantic.affectedResultKeys) {
        const svg = preparedSvgByResultKey instanceof Map
          ? preparedSvgByResultKey.get(resultKey)
          : preparedSvgByResultKey?.[resultKey];
        if (!svg?.cloneNode) {
          throw new Error(`Feature label reflow projection is incomplete for Result ${resultKey}.`);
        }
      }
    }
    const targetFeatureKeysByResult = new Map(
      plan.targetsByResult.map((entry) => [entry.resultKey, cloneJson(entry.featureKeys)])
    );
    const projection = await prepareResultProjection({
      direction, from, to, plan, catalog, results: sourceResults,
      labelsByTargetKey,
      affectedResultKeys: semantic.affectedResultKeys,
      targetFeatureKeysByResult,
      semanticScope: semantic.projectionMatch.semanticScope,
      matchText: semantic.projectionMatch.matchText,
      sourceText: semantic.projectionMatch.sourceText,
      mounted,
      preparedSvgByResultKey
    });
    if (!projection || !Array.isArray(projection.nextResults)) {
      throw new Error('Feature label Result projection did not prepare a complete Result array.');
    }
    if (readRef(state.results) !== sourceResults || projection.nextResults.length !== sourceResults.length) {
      throw new Error('Feature label Result projection changed state or Result topology during preflight.');
    }
    const affectedIndexes = new Set(semantic.affectedResultKeys.map((key) => resultIndexByKey.get(key)));
    projection.nextResults.forEach((result, index) => {
      if (!result || text(result.name) !== text(sourceResults[index]?.name)) {
        throw new Error(`Feature label Result identity changed at index ${index}.`);
      }
      if (!affectedIndexes.has(index) && result !== sourceResults[index]) {
        throw new Error(`Feature label projection replaced unaffected Result ${index}.`);
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
      }),
      reflow: reflow?.retainedForHistory || null
    };
    return {
      projection, reflow, sourceResults, sourceResultRecords, preflightMounted: mounted,
      retainedForHistory,
      retainedBytes: exactJsonByteLength(retainedForHistory),
      counters: {
        projectedResults: Number(projection.counters?.affectedResults || semantic.affectedResultKeys.length),
        mountedSerializations: Number(projection.counters?.mountedSerializations || 0),
        detachedPasses: Number(projection.counters?.detachedPasses || 0),
        changedResults: Number(projection.counters?.changedResults || 0),
        changedLabels: Number(projection.counters?.changedLabels || 0),
        reflowedResults: Number(reflow?.counters?.affectedResults || 0),
        resultArraySwaps: 1
      }
    };
  };

  const captureExactState = () => {
    const runtime = typeof captureRuntimeState === 'function' ? captureRuntimeState() : null;
    if (runtime && typeof runtime.then === 'function') {
      throw new Error('Feature label runtime capture must be synchronous.');
    }
    return {
      featureOverrides: { ...state.labelTextFeatureOverrides },
      bulkOverrides: { ...state.labelTextBulkOverrides },
      featureOverrideSources: { ...state.labelTextFeatureOverrideSources },
      results: readRef(state.results),
      revision: Number(readRef(state.semanticStyleRevision) || 0),
      fingerprint: text(readRef(state.semanticStyleFingerprint)),
      ledger: readRef(state.validatedStyleFingerprintByResultKey),
      expectedRevision,
      expectedResults,
      expectedResultSignatures: { ...expectedResultSignatures },
      presentation: {
        ...(isWritableRef(state.clickedFeature) ? { clickedFeature: readRef(state.clickedFeature) } : {}),
        ...(isWritableRef(state.editableLabels) ? { editableLabels: readRef(state.editableLabels) } : {}),
        ...(isWritableRef(state.selectedFeatures) ? { selectedFeatures: readRef(state.selectedFeatures) } : {}),
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
        throw new Error('Mounted Feature label rollback must be synchronous.');
      }
      if (restored === false) return false;
    }
    restoreObjectReferences(state.labelTextFeatureOverrides, snapshot.featureOverrides);
    restoreObjectReferences(state.labelTextBulkOverrides, snapshot.bulkOverrides);
    restoreObjectReferences(state.labelTextFeatureOverrideSources, snapshot.featureOverrideSources);
    writeRef(state.results, snapshot.results);
    writeRef(state.semanticStyleRevision, snapshot.revision);
    if (isWritableRef(state.semanticStyleFingerprint)) writeRef(state.semanticStyleFingerprint, snapshot.fingerprint);
    if (isWritableRef(state.validatedStyleFingerprintByResultKey)) {
      writeRef(state.validatedStyleFingerprintByResultKey, snapshot.ledger);
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
        throw new Error('Feature label runtime rollback must be synchronous.');
      }
      if (restored === false) return false;
    }
    return true;
  };

  const commitTransition = ({ to, direction, prepared }) => {
    replaceObject(state.labelTextFeatureOverrides, to.featureOverrides);
    replaceObject(state.labelTextBulkOverrides, to.bulkOverrides);
    replaceObject(state.labelTextFeatureOverrideSources, to.featureOverrideSources);
    if (semantic.affectedResultKeys.includes(prepared.preflightMounted.resultKey)) {
      if (typeof commitMountedProjection !== 'function' || typeof restoreMountedProjection !== 'function') {
        throw new Error('Affected mounted Feature label projection has no atomic commit/restore owner.');
      }
      const committed = commitMountedProjection({ direction, prepared: prepared.projection, plan });
      if (committed && typeof committed.then === 'function') {
        throw new Error('Mounted Feature label commit must be synchronous.');
      }
      if (committed === false) return false;
    }
    writeRef(state.results, prepared.projection.nextResults);
    const published = readRef(state.results);
    if (!projectedResultsMatch(published, prepared.projection.nextResults)) {
      throw new Error('Feature label Result publication did not preserve the staged projection.');
    }
    const reconciled = typeof reconcile === 'function'
      ? reconcile({ direction, plan, prepared, to })
      : undefined;
    if (reconciled && typeof reconciled.then === 'function') {
      throw new Error('Feature label reconciliation must be synchronous.');
    }
    if (reconciled === false) return false;
    const refreshed = typeof refreshPresentation === 'function'
      ? refreshPresentation({ direction, plan, prepared, to })
      : undefined;
    if (refreshed && typeof refreshed.then === 'function') {
      throw new Error('Feature label presentation refresh must be synchronous.');
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
    isNoop: (before, after) => before.labelFingerprint === after.labelFingerprint,
    validateCurrent,
    prepareTransition,
    commitTransition,
    captureExactState,
    restoreExactState,
    assertCommitted: ({ to, prepared }) => (
      currentLabelFingerprint(state) === to.labelFingerprint
      && mapsMatch(state, to)
      && Number(readRef(state.semanticStyleRevision) || 0) === expectedRevision
      && currentStyleFingerprint(state) === expectedStyleFingerprint
      && ledgerIsCurrent(state, catalog, expectedStyleFingerprint)
      && readRef(state.validatedStyleFingerprintByResultKey) === expectedLedger
      && readRef(state.results) === prepared.publishedResults
      && projectedResultsMatch(readRef(state.results), prepared.projection.nextResults)
    )
  });
};
