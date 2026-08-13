import { catalogResultKey } from '../../services/feature-catalog.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
import {
  normalizeStyleSnapshot,
  styleFingerprint
} from '../../services/style-revision.js';
import { deriveUsedFeatureFillGroupsByResult } from '../legend/feature-fill-projection.js';
import { normalizeSpecificRule } from '../specific-color-rules.js';
import { deriveFeatureFillOverrides } from './fill-command.js';
import { resolveOrderedSpecificRule } from './fill-view-model.js';
import { stableFeatureTargetKey } from './fill-scope-plan.js';
import {
  buildStyleSnapshotCommand,
  exactJsonByteLength,
  immutableStyleCommandSnapshot
} from './style-snapshot-command.js';

const text = (value) => String(value ?? '').trim();

const readRef = (value) => (
  value && typeof value === 'object' && 'value' in value
    ? value.value
    : value
);

const isWritableRef = (value) => Boolean(
  value
  && typeof value === 'object'
  && 'value' in value
  && value.__v_isReadonly !== true
);

const writeRef = (target, value) => {
  if (!isWritableRef(target)) throw new Error('Bulk Feature style command requires a writable ref.');
  target.value = value;
};

const normalizedRules = (rules) => (Array.isArray(rules) ? rules : []).map((rule, index) => {
  const normalized = normalizeSpecificRule(rule, { fromFile: Boolean(rule?.fromFile) });
  if (!normalized.feat || !normalized.qual || !normalized.val || !normalized.color) {
    throw new Error(`Specific-color rule ${index + 1} is incomplete.`);
  }
  return normalized;
});

const stableJson = (value) => {
  if (Array.isArray(value)) return value.map(stableJson);
  if (!value || typeof value !== 'object') return value;
  return Object.fromEntries(
    Object.keys(value).sort().map((key) => [key, stableJson(value[key])])
  );
};

const sameJson = (left, right) => (
  JSON.stringify(stableJson(left)) === JSON.stringify(stableJson(right))
);

const cloneRules = (rules) => normalizedRules(rules).map((rule) => ({ ...rule }));

const replaceRules = (target, rules) => {
  if (!Array.isArray(target)) throw new Error('Specific-color rule state is unavailable.');
  target.splice(0, target.length, ...cloneRules(rules));
};

const restoreRuleReferences = (target, rules) => {
  if (!Array.isArray(target) || !Array.isArray(rules)) {
    throw new Error('Specific-color rules cannot be restored exactly.');
  }
  target.splice(0, target.length, ...rules);
};

const replaceObject = (target, source) => {
  if (!target || typeof target !== 'object' || Array.isArray(target)) {
    throw new Error('Feature override cache is unavailable.');
  }
  Object.keys(target).forEach((key) => delete target[key]);
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

/** Normalize one complete ordered rule and applied-palette semantic side. */
export const normalizeFeatureBulkStyleSnapshot = (snapshot = {}) => {
  const normalized = normalizeStyleSnapshot({
    rules: snapshot.rules || [],
    appliedPaletteName: snapshot.appliedPaletteName ?? snapshot.paletteName ?? 'default',
    appliedPaletteColors: snapshot.appliedPaletteColors ?? snapshot.paletteColors ?? {}
  });
  return immutableStyleCommandSnapshot({
    rules: cloneRules(snapshot.rules || []),
    appliedPaletteName: normalized.appliedPaletteName,
    appliedPaletteColors: cloneJson(normalized.appliedPaletteColors)
  });
};

/** Copy one semantic side while replacing any complete writer-owned field. */
export const replaceFeatureBulkStyleSnapshot = (snapshot, replacement = {}) => (
  normalizeFeatureBulkStyleSnapshot({
    rules: replacement.rules ?? snapshot?.rules ?? [],
    appliedPaletteName: replacement.appliedPaletteName
      ?? replacement.paletteName
      ?? snapshot?.appliedPaletteName
      ?? 'default',
    appliedPaletteColors: replacement.appliedPaletteColors
      ?? replacement.paletteColors
      ?? snapshot?.appliedPaletteColors
      ?? {}
  })
);

const resultBindings = (catalog) => {
  const bindings = [];
  const seen = new Set();
  (Array.isArray(catalog?.items) ? catalog.items : []).forEach((item, resultIndex) => {
    const resultKey = catalogResultKey(item);
    if (!resultKey || seen.has(resultKey)) {
      throw new Error(`Bulk Feature style catalogue has missing or duplicate Result identity at ${resultIndex}.`);
    }
    if (item?.resultIndex !== undefined && Number(item.resultIndex) !== resultIndex) {
      throw new Error(`Bulk Feature style catalogue index is stale for Result ${resultKey}.`);
    }
    seen.add(resultKey);
    const renderedPairs = new Set(
      (Array.isArray(item?.features) ? item.features : []).map((feature) => JSON.stringify([
        text(feature?.recordKey),
        text(feature?.biologicalFeatureId)
      ]))
    );
    const targetKeys = [];
    const seenTargetKeys = new Set();
    const renderedTargetKeys = new Set();
    (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).forEach((feature) => {
      const targetKey = stableFeatureTargetKey({ ...feature, resultKey });
      if (!targetKey || seenTargetKeys.has(targetKey)) {
        throw new Error(`Bulk Feature style target identity is missing or duplicate in Result ${resultKey}.`);
      }
      seenTargetKeys.add(targetKey);
      targetKeys.push(targetKey);
      const pair = JSON.stringify([
        text(feature?.recordKey),
        text(feature?.biologicalFeatureId)
      ]);
      if (
        renderedPairs.has(pair)
        && feature?.hidden !== true
        && feature?.rendered !== false
      ) renderedTargetKeys.add(targetKey);
    });
    bindings.push({ item, resultIndex, resultKey, targetKeys, renderedTargetKeys });
  });
  return bindings;
};

const canonicalLegendEntries = (entries) => (Array.isArray(entries) ? entries : []).map((entry) => ({
  caption: text(entry?.caption),
  color: text(entry?.color).toLowerCase(),
  owner: text(entry?.owner),
  featureIds: [...new Set((Array.isArray(entry?.featureIds) ? entry.featureIds : [])
    .map(text)
    .filter(Boolean))].sort()
}));

const legendProjectionMap = (derived) => new Map(
  (Array.isArray(derived?.projections) ? derived.projections : []).map((projection) => [
    text(projection?.resultKey),
    projection
  ])
);

const styleIdentity = (style) => ({
  rules: style.rules,
  appliedPaletteName: style.appliedPaletteName,
  appliedPaletteColors: style.appliedPaletteColors
});

/**
 * Compare complete bulk semantic snapshots over every catalogue Feature.
 * The returned target map includes only rendered Features whose effective
 * fill/caption changed; hidden Features still participate in the override cache.
 */
export const deriveFeatureBulkStyleChange = ({
  catalog,
  before,
  after,
  manualLegendEntries = [],
  existingEntriesByResult = {},
  nonFeatureResultKeys = [],
  selectedFeatureTypes = [],
  resolveFeatureVisibility = null,
  svgDefaultColor = '#cccccc'
} = {}) => {
  const beforeStyle = normalizeFeatureBulkStyleSnapshot(before);
  const afterStyle = normalizeFeatureBulkStyleSnapshot(after);
  const bindings = resultBindings(catalog);
  const beforeDerived = deriveFeatureFillOverrides({
    catalog,
    rules: beforeStyle.rules,
    paletteColors: beforeStyle.appliedPaletteColors,
    svgDefaultColor
  });
  const afterDerived = deriveFeatureFillOverrides({
    catalog,
    rules: afterStyle.rules,
    paletteColors: afterStyle.appliedPaletteColors,
    svgDefaultColor
  });
  const beforeLegend = deriveUsedFeatureFillGroupsByResult({
    catalog,
    rules: beforeStyle.rules,
    paletteColors: beforeStyle.appliedPaletteColors,
    manualLegendEntries,
    existingEntriesByResult,
    svgDefaultColor
  });
  const afterLegend = deriveUsedFeatureFillGroupsByResult({
    catalog,
    rules: afterStyle.rules,
    paletteColors: afterStyle.appliedPaletteColors,
    manualLegendEntries,
    existingEntriesByResult,
    svgDefaultColor
  });
  const beforeLegendByResult = legendProjectionMap(beforeLegend);
  const afterLegendByResult = legendProjectionMap(afterLegend);
  const featureAffectedResultKeys = [];
  const geometryChangedResultKeys = [];
  const targetFeatureKeysByResult = {};
  let changedRenderedFeatures = 0;
  let changedGeometryFeatures = 0;
  const selectedTypes = new Set(
    Array.from(selectedFeatureTypes || []).map(text).filter(Boolean)
  );

  bindings.forEach((binding) => {
    const changedTargets = binding.targetKeys.filter((targetKey) => {
      const beforeValue = beforeDerived.effectiveByTargetKey[targetKey];
      const afterValue = afterDerived.effectiveByTargetKey[targetKey];
      return binding.renderedTargetKeys.has(targetKey) && !sameJson(beforeValue, afterValue);
    });
    const legendChanged = !sameJson(
      canonicalLegendEntries(beforeLegendByResult.get(binding.resultKey)?.entries),
      canonicalLegendEntries(afterLegendByResult.get(binding.resultKey)?.entries)
    );
    const geometryTargets = binding.targetKeys.filter((targetKey, featureIndex) => {
      const feature = binding.item.biologicalFeatures[featureIndex];
      const currentlyRendered = binding.renderedTargetKeys.has(targetKey);
      const beforeRule = resolveOrderedSpecificRule(feature, beforeStyle.rules)?.rule || null;
      const afterRule = resolveOrderedSpecificRule(feature, afterStyle.rules)?.rule || null;
      let beforeVisible;
      let afterVisible;
      if (typeof resolveFeatureVisibility === 'function') {
        beforeVisible = resolveFeatureVisibility({
          phase: 'before', feature, item: binding.item, resultKey: binding.resultKey,
          style: beforeStyle, matchedRule: beforeRule, currentlyRendered
        });
        afterVisible = resolveFeatureVisibility({
          phase: 'after', feature, item: binding.item, resultKey: binding.resultKey,
          style: afterStyle, matchedRule: afterRule, currentlyRendered
        });
        if (typeof beforeVisible !== 'boolean' || typeof afterVisible !== 'boolean') {
          throw new Error('Bulk Feature style visibility resolver must return booleans.');
        }
        return beforeVisible !== afterVisible;
      }
      const typeSelected = selectedTypes.size === 0 || selectedTypes.has(text(feature?.type));
      const beforePotential = typeSelected || Boolean(beforeRule);
      const afterPotential = typeSelected || Boolean(afterRule);
      if (beforePotential === afterPotential) return false;
      // Visibility rules can override the base selection. Treat every base
      // transition as geometry-sensitive unless an exact resolver is supplied.
      return true;
    });
    const allTargets = [...new Set([...changedTargets, ...geometryTargets])];
    if (changedTargets.length > 0 || legendChanged || geometryTargets.length > 0) {
      if (allTargets.length === 0) {
        throw new Error(`Bulk Feature style legend changed without a rendered Feature target in ${binding.resultKey}.`);
      }
      featureAffectedResultKeys.push(binding.resultKey);
      targetFeatureKeysByResult[binding.resultKey] = Object.freeze(allTargets);
      changedRenderedFeatures += changedTargets.length;
      changedGeometryFeatures += geometryTargets.length;
      if (geometryTargets.length > 0) geometryChangedResultKeys.push(binding.resultKey);
    }
  });

  const beforeFingerprint = styleFingerprint(beforeStyle);
  const afterFingerprint = styleFingerprint(afterStyle);
  const knownResultKeys = new Set(bindings.map(({ resultKey }) => resultKey));
  const nonFeatureAffectedResultKeys = [...new Set(
    Array.from(nonFeatureResultKeys || []).map(text).filter(Boolean)
  )];
  nonFeatureAffectedResultKeys.forEach((resultKey) => {
    if (!knownResultKeys.has(resultKey)) {
      throw new Error(`Non-Feature palette target Result is unavailable: ${resultKey}.`);
    }
  });
  const featureAffected = new Set(featureAffectedResultKeys);
  const nonFeatureAffected = new Set(nonFeatureAffectedResultKeys);
  const affectedResultKeys = bindings
    .map(({ resultKey }) => resultKey)
    .filter((resultKey) => featureAffected.has(resultKey) || nonFeatureAffected.has(resultKey));
  const resultFingerprintByKey = Object.freeze(Object.fromEntries(
    bindings.map(({ resultKey }) => [resultKey, afterFingerprint])
  ));
  return immutableStyleCommandSnapshot({
    before: {
      ...beforeStyle,
      overrides: beforeDerived.overrides,
      effectiveByTargetKey: beforeDerived.effectiveByTargetKey,
      legendProjections: beforeLegend.projections,
      fingerprint: beforeFingerprint,
      identity: styleIdentity(beforeStyle)
    },
    after: {
      ...afterStyle,
      overrides: afterDerived.overrides,
      effectiveByTargetKey: afterDerived.effectiveByTargetKey,
      legendProjections: afterLegend.projections,
      fingerprint: afterFingerprint,
      identity: styleIdentity(afterStyle),
      resultFingerprintByKey
    },
    affectedResultKeys,
    featureAffectedResultKeys,
    nonFeatureAffectedResultKeys,
    geometryChangedResultKeys,
    unaffectedResultKeys: bindings
      .map(({ resultKey }) => resultKey)
      .filter((resultKey) => !affectedResultKeys.includes(resultKey)),
    targetFeatureKeysByResult,
    resultKeys: bindings.map(({ resultKey }) => resultKey),
    counters: {
      results: bindings.length,
      affectedResults: affectedResultKeys.length,
      featureAffectedResults: featureAffectedResultKeys.length,
      nonFeatureAffectedResults: nonFeatureAffectedResultKeys.length,
      geometryChangedResults: geometryChangedResultKeys.length,
      changedRenderedFeatures,
      changedGeometryFeatures,
      beforeRuleResolutionUpperBound:
        beforeDerived.counters.ruleResolutionUpperBound + beforeLegend.counters.ruleResolutionUpperBound,
      afterRuleResolutionUpperBound:
        afterDerived.counters.ruleResolutionUpperBound + afterLegend.counters.ruleResolutionUpperBound
    }
  });
};

const stateStyle = (state) => normalizeFeatureBulkStyleSnapshot({
  rules: state?.manualSpecificRules,
  appliedPaletteName: readRef(state?.appliedPaletteName),
  appliedPaletteColors: readRef(state?.appliedPaletteColors)
});

const stateFingerprint = (state) => styleFingerprint(stateStyle(state));

const currentMounted = (getMountedContext) => {
  const raw = typeof getMountedContext === 'function' ? getMountedContext() : null;
  if (!raw) return { resultIndex: null, resultKey: '', svg: null };
  const number = Number(raw.resultIndex);
  return {
    ...raw,
    resultIndex: Number.isInteger(number) && number >= 0 ? number : null,
    resultKey: text(raw.resultKey),
    svg: raw.svg || null
  };
};

const sameMount = (left, right) => (
  (left?.resultIndex ?? null) === (right?.resultIndex ?? null)
  && text(left?.resultKey) === text(right?.resultKey)
  && (left?.svg || null) === (right?.svg || null)
);

const mountIsCurrent = (mounted, catalog) => {
  if (mounted.resultIndex === null) return !mounted.resultKey && !mounted.svg;
  return (
    mounted.resultIndex < (catalog?.items?.length || 0)
    && catalogResultKey(catalog.items[mounted.resultIndex]) === mounted.resultKey
  );
};

const selectedResult = (state, catalog) => {
  const index = Number(readRef(state?.selectedResultIndex));
  if (!Number.isInteger(index) || index < 0 || index >= (catalog?.items?.length || 0)) {
    return { resultIndex: null, resultKey: '' };
  }
  return { resultIndex: index, resultKey: catalogResultKey(catalog.items[index]) };
};

const sameSelectedResult = (left, right) => (
  (left?.resultIndex ?? null) === (right?.resultIndex ?? null)
  && text(left?.resultKey) === text(right?.resultKey)
);

const resultContent = (result) => String(result?.content ?? '');

const resultsMatchCatalog = (results, catalog) => (
  Array.isArray(results)
  && results.length === (catalog?.items?.length || 0)
  && results.every((result, resultIndex) => {
    const expectedName = text(catalog.items[resultIndex]?.resultName);
    return result && (!expectedName || text(result.name) === expectedName);
  })
);

const resultGuard = (results) => (Array.isArray(results) ? results : []).map((result) => ({
  result,
  name: text(result?.name),
  content: resultContent(result)
}));

const resultsAreCurrent = (results, guard) => (
  Array.isArray(results)
  && results.length === guard.length
  && guard.every((entry, index) => (
    results[index] === entry.result
    && text(results[index]?.name) === entry.name
    && resultContent(results[index]) === entry.content
  ))
);

const publishedResultsMatch = (published, projected) => (
  Array.isArray(published)
  && Array.isArray(projected)
  && published.length === projected.length
  && published.every((result, index) => (
    text(result?.name) === text(projected[index]?.name)
    && resultContent(result) === resultContent(projected[index])
  ))
);

const ledgerMatches = (ledger, resultKeys, fingerprint) => (
  resultKeys.every((resultKey) => text(ledger?.[resultKey]) === fingerprint)
);

const catalogSignature = (catalog) => JSON.stringify(catalog);

const preparedSvgMap = (legend) => {
  if (legend?.preparedSvgByResultKey instanceof Map) return legend.preparedSvgByResultKey;
  if (legend?.preparedSvgByResultKey && typeof legend.preparedSvgByResultKey === 'object') {
    return legend.preparedSvgByResultKey;
  }
  if (legend?.staged instanceof Map) {
    return new Map([...legend.staged].map(([resultKey, entry]) => [resultKey, entry?.svg || entry]));
  }
  return null;
};

const rootForResult = (roots, resultKey) => (
  roots instanceof Map ? roots.get(resultKey) : roots?.[resultKey]
);

const legendForResult = (side, resultKey) => (
  (Array.isArray(side?.legendProjections) ? side.legendProjections : [])
    .find((projection) => text(projection?.resultKey) === resultKey)
);

const paletteUiState = (state) => {
  const values = {};
  if (isWritableRef(state?.selectedPalette)) {
    values.selectedPalette = readRef(state.selectedPalette);
  }
  if (isWritableRef(state?.currentColors)) {
    values.currentColors = cloneJson(readRef(state.currentColors) || {});
  }
  if (isWritableRef(state?.pendingPaletteName)) {
    values.pendingPaletteName = readRef(state.pendingPaletteName);
  }
  if (isWritableRef(state?.pendingPaletteColors)) {
    values.pendingPaletteColors = cloneJson(readRef(state.pendingPaletteColors) || {});
  }
  return values;
};

const acceptedPaletteUiState = (state, style) => {
  const values = {};
  if (isWritableRef(state?.selectedPalette)) values.selectedPalette = style.appliedPaletteName;
  if (isWritableRef(state?.currentColors)) {
    values.currentColors = cloneJson(style.appliedPaletteColors);
  }
  if (isWritableRef(state?.pendingPaletteName)) values.pendingPaletteName = '';
  if (isWritableRef(state?.pendingPaletteColors)) values.pendingPaletteColors = {};
  return values;
};

const paletteUiHasPendingDraft = (values) => (
  text(values?.pendingPaletteName) !== ''
);

const auxiliaryState = ({ palette, legend = null }) => ({
  palette: cloneJson(palette || {}),
  ...(legend ? { legend: cloneJson(legend) } : {})
});

const legendPresentationState = (state, patch = null) => {
  if (!patch || !state?.adv || typeof state.adv !== 'object') return null;
  const values = {};
  if (Object.prototype.hasOwnProperty.call(patch, 'legendBoxSize')) {
    values.legendBoxSize = state.adv.legend_box_size;
  }
  if (Object.prototype.hasOwnProperty.call(patch, 'legendFontSize')) {
    values.legendFontSize = state.adv.legend_font_size;
  }
  return values;
};

const replaceAuxiliaryState = (state, values = {}) => {
  const palette = values.palette || {};
  Object.entries({
    selectedPalette: palette.selectedPalette,
    currentColors: palette.currentColors,
    pendingPaletteName: palette.pendingPaletteName,
    pendingPaletteColors: palette.pendingPaletteColors
  }).forEach(([key, value]) => {
    if (!Object.prototype.hasOwnProperty.call(palette, key)) return;
    writeRef(state[key], cloneJson(value));
  });
  if (values.legend && state?.adv && typeof state.adv === 'object') {
    if (Object.prototype.hasOwnProperty.call(values.legend, 'legendBoxSize')) {
      state.adv.legend_box_size = values.legend.legendBoxSize;
    }
    if (Object.prototype.hasOwnProperty.call(values.legend, 'legendFontSize')) {
      state.adv.legend_font_size = values.legend.legendFontSize;
    }
  }
};

const auxiliaryStateMatches = (state, values = {}) => {
  const current = auxiliaryState({
    palette: paletteUiState(state),
    legend: values.legend ? legendPresentationState(state, {
      ...(Object.prototype.hasOwnProperty.call(values.legend, 'legendBoxSize')
        ? { legendBoxSize: true }
        : {}),
      ...(Object.prototype.hasOwnProperty.call(values.legend, 'legendFontSize')
        ? { legendFontSize: true }
        : {})
    }) : null
  });
  return sameJson(current, values);
};

const captureExactState = ({ state, captureRuntimeState }) => {
  const runtime = typeof captureRuntimeState === 'function' ? captureRuntimeState() : null;
  if (runtime && typeof runtime.then === 'function') {
    throw new Error('Bulk Feature style runtime capture must be synchronous.');
  }
  return {
    rules: [...state.manualSpecificRules],
    appliedPaletteName: readRef(state.appliedPaletteName),
    appliedPaletteColors: readRef(state.appliedPaletteColors),
    selectedPalette: isWritableRef(state.selectedPalette)
      ? readRef(state.selectedPalette)
      : undefined,
    currentColors: isWritableRef(state.currentColors)
      ? readRef(state.currentColors)
      : undefined,
    pendingPaletteName: isWritableRef(state.pendingPaletteName)
      ? readRef(state.pendingPaletteName)
      : undefined,
    pendingPaletteColors: isWritableRef(state.pendingPaletteColors)
      ? readRef(state.pendingPaletteColors)
      : undefined,
    legendPresentation: state?.adv && typeof state.adv === 'object'
      ? {
          hasBoxSize: Object.prototype.hasOwnProperty.call(state.adv, 'legend_box_size'),
          boxSize: state.adv.legend_box_size,
          hasFontSize: Object.prototype.hasOwnProperty.call(state.adv, 'legend_font_size'),
          fontSize: state.adv.legend_font_size
        }
      : null,
    overrides: { ...state.featureColorOverrides },
    results: readRef(state.results),
    legendEntries: readRef(state.legendEntries),
    manualLegendEntries: isWritableRef(state.manualLegendEntries)
      ? readRef(state.manualLegendEntries)
      : undefined,
    ledger: readRef(state.validatedStyleFingerprintByResultKey),
    revision: readRef(state.semanticStyleRevision),
    fingerprint: readRef(state.semanticStyleFingerprint),
    documentEpoch: readRef(state.documentEpoch),
    resultGenerationKey: readRef(state.resultGenerationKey),
    presentation: Object.fromEntries(
      ['clickedFeature', 'selectedFeatures', 'selectedResultIndex']
        .filter((key) => isWritableRef(state[key]))
        .map((key) => [key, readRef(state[key])])
    ),
    runtime
  };
};

const restoreExactState = ({ state, snapshot, restoreRuntimeState }) => {
  restoreRuleReferences(state.manualSpecificRules, snapshot.rules);
  writeRef(state.appliedPaletteName, snapshot.appliedPaletteName);
  writeRef(state.appliedPaletteColors, snapshot.appliedPaletteColors);
  if (snapshot.selectedPalette !== undefined && isWritableRef(state.selectedPalette)) {
    writeRef(state.selectedPalette, snapshot.selectedPalette);
  }
  if (snapshot.currentColors !== undefined && isWritableRef(state.currentColors)) {
    writeRef(state.currentColors, snapshot.currentColors);
  }
  if (snapshot.pendingPaletteName !== undefined && isWritableRef(state.pendingPaletteName)) {
    writeRef(state.pendingPaletteName, snapshot.pendingPaletteName);
  }
  if (snapshot.pendingPaletteColors !== undefined && isWritableRef(state.pendingPaletteColors)) {
    writeRef(state.pendingPaletteColors, snapshot.pendingPaletteColors);
  }
  if (snapshot.legendPresentation && state?.adv && typeof state.adv === 'object') {
    if (snapshot.legendPresentation.hasBoxSize) {
      state.adv.legend_box_size = snapshot.legendPresentation.boxSize;
    } else {
      delete state.adv.legend_box_size;
    }
    if (snapshot.legendPresentation.hasFontSize) {
      state.adv.legend_font_size = snapshot.legendPresentation.fontSize;
    } else {
      delete state.adv.legend_font_size;
    }
  }
  restoreObjectReferences(state.featureColorOverrides, snapshot.overrides);
  writeRef(state.results, snapshot.results);
  writeRef(state.legendEntries, snapshot.legendEntries);
  if (snapshot.manualLegendEntries !== undefined && isWritableRef(state.manualLegendEntries)) {
    writeRef(state.manualLegendEntries, snapshot.manualLegendEntries);
  }
  writeRef(state.validatedStyleFingerprintByResultKey, snapshot.ledger);
  writeRef(state.semanticStyleRevision, snapshot.revision);
  writeRef(state.semanticStyleFingerprint, snapshot.fingerprint);
  writeRef(state.documentEpoch, snapshot.documentEpoch);
  writeRef(state.resultGenerationKey, snapshot.resultGenerationKey);
  Object.entries(snapshot.presentation).forEach(([key, value]) => writeRef(state[key], value));
  if (typeof restoreRuntimeState === 'function') {
    const restored = restoreRuntimeState(snapshot.runtime);
    if (restored && typeof restored.then === 'function') {
      throw new Error('Bulk Feature style runtime rollback must be synchronous.');
    }
    if (restored === false) return false;
  }
  return true;
};

/**
 * Build one History command for any accepted complete rule/preset/import or
 * applied-palette/default-color snapshot.
 */
export const buildFeatureBulkStyleCommand = async ({
  state,
  catalog,
  before,
  after,
  writerKind = 'bulk-style',
  label = 'Change feature styles',
  prepareLegendProjection,
  prepareResultProjection,
  prepareGeometryProjection = null,
  getMountedContext = null,
  commitMountedProjection = null,
  restoreMountedProjection = null,
  captureRuntimeState = null,
  restoreRuntimeState = null,
  reconcile = null,
  refreshPresentation = null,
  existingEntriesByResult = {},
  nonFeatureAffectedResultKeys = [],
  presentationPatch = null,
  selectedFeatureTypes = null,
  resolveFeatureVisibility = null,
  svgDefaultColor = '#cccccc'
} = {}) => {
  if (!state || !catalog || (
    typeof prepareResultProjection !== 'function'
    && typeof prepareGeometryProjection !== 'function'
  )) {
    throw new Error('Bulk Feature style command requires state, catalogue, and Result preparation.');
  }
  const manualLegendEntries = cloneJson(readRef(state.manualLegendEntries) || []);
  const change = deriveFeatureBulkStyleChange({
    catalog,
    before,
    after,
    manualLegendEntries,
    existingEntriesByResult,
    nonFeatureResultKeys: nonFeatureAffectedResultKeys,
    selectedFeatureTypes: selectedFeatureTypes
      ?? Array.from(state?.adv?.features || []),
    resolveFeatureVisibility,
    svgDefaultColor
  });
  const requestPaletteUi = paletteUiState(state);
  const paletteChanged = !sameJson(
    {
      name: change.before.appliedPaletteName,
      colors: change.before.appliedPaletteColors
    },
    {
      name: change.after.appliedPaletteName,
      colors: change.after.appliedPaletteColors
    }
  );
  const beforePaletteUi = paletteChanged && !paletteUiHasPendingDraft(requestPaletteUi)
    ? acceptedPaletteUiState(state, change.before)
    : requestPaletteUi;
  const afterPaletteUi = paletteChanged
    ? acceptedPaletteUiState(state, change.after)
    : requestPaletteUi;
  const beforeLegendPresentation = legendPresentationState(state, presentationPatch);
  const afterLegendPresentation = beforeLegendPresentation
    ? {
        ...beforeLegendPresentation,
        ...(Object.prototype.hasOwnProperty.call(presentationPatch, 'legendBoxSize')
          ? { legendBoxSize: presentationPatch.legendBoxSize }
          : {}),
        ...(Object.prototype.hasOwnProperty.call(presentationPatch, 'legendFontSize')
          ? { legendFontSize: presentationPatch.legendFontSize }
          : {})
      }
    : null;
  const commandBefore = {
    ...change.before,
    auxiliary: auxiliaryState({
      palette: beforePaletteUi,
      legend: beforeLegendPresentation
    })
  };
  const commandAfter = {
    ...change.after,
    auxiliary: auxiliaryState({
      palette: afterPaletteUi,
      legend: afterLegendPresentation
    })
  };
  const expected = Object.freeze({
    documentEpoch: Number(readRef(state.documentEpoch) || 0),
    resultGenerationKey: Number(readRef(state.resultGenerationKey) || 0),
    revision: Number(readRef(state.semanticStyleRevision) || 0),
    catalogReference: readRef(state.featureCatalog) || catalog,
    catalogSignature: catalogSignature(catalog)
  });
  const initialResults = readRef(state.results);
  if (
    expected.catalogReference !== catalog
    || !resultsMatchCatalog(initialResults, catalog)
    || !sameJson(stateStyle(state), change.before.identity)
    || stateFingerprint(state) !== change.before.fingerprint
    || text(readRef(state.semanticStyleFingerprint)) !== change.before.fingerprint
    || !sameJson(state.featureColorOverrides, change.before.overrides)
    || !ledgerMatches(
      readRef(state.validatedStyleFingerprintByResultKey) || {},
      change.resultKeys,
      change.before.fingerprint
    )
  ) {
    throw new Error('The bulk Feature style snapshot is stale before command preparation.');
  }

  const expectedRevision = { apply: expected.revision, undo: null, redo: null };
  const validateCurrent = ({ from, direction, prepared = null }) => {
    const mounted = currentMounted(getMountedContext);
    const results = readRef(state.results);
    const revision = Number(readRef(state.semanticStyleRevision) || 0);
    const expectedAuxiliary = direction === 'apply'
      ? auxiliaryState({
          palette: requestPaletteUi,
          legend: beforeLegendPresentation
        })
      : from.auxiliary;
    if (
      Number(readRef(state.documentEpoch) || 0) !== expected.documentEpoch
      || Number(readRef(state.resultGenerationKey) || 0) !== expected.resultGenerationKey
      || revision !== expectedRevision[direction]
      || (readRef(state.featureCatalog) || catalog) !== catalog
      || catalogSignature(catalog) !== expected.catalogSignature
      || !sameJson(stateStyle(state), from.identity)
      || stateFingerprint(state) !== from.fingerprint
      || text(readRef(state.semanticStyleFingerprint)) !== from.fingerprint
      || !sameJson(state.featureColorOverrides, from.overrides)
      || !ledgerMatches(
        readRef(state.validatedStyleFingerprintByResultKey) || {},
        change.resultKeys,
        from.fingerprint
      )
      || !Array.isArray(results)
      || !resultsMatchCatalog(results, catalog)
      || !mountIsCurrent(mounted, catalog)
      || !auxiliaryStateMatches(state, expectedAuxiliary)
    ) return false;
    if (
      mounted.resultKey
      && change.affectedResultKeys.includes(mounted.resultKey)
      && !mounted.svg
    ) return false;
    if (prepared && (
      results !== prepared.sourceResults
      || !resultsAreCurrent(results, prepared.sourceResultGuard)
      || !sameMount(mounted, prepared.preflightMounted)
      || !sameSelectedResult(selectedResult(state, catalog), prepared.preflightSelectedResult)
      || readRef(state.legendEntries) !== prepared.sourceLegendEntries
      || !sameJson(readRef(state.legendEntries), prepared.sourceLegendValue)
      || !sameJson(readRef(state.manualLegendEntries) || [], manualLegendEntries)
    )) return false;
    return true;
  };

  const prepareTransition = async ({ from, to, direction }) => {
    const sourceResults = readRef(state.results);
    const sourceResultGuard = resultGuard(sourceResults);
    const sourceLegendEntries = readRef(state.legendEntries);
    const sourceLegendValue = cloneJson(sourceLegendEntries || []);
    const preflightMounted = currentMounted(getMountedContext);
    const preflightSelectedResult = selectedResult(state, catalog);
    const projections = change.affectedResultKeys.map((resultKey) => {
      const projection = legendForResult(to, resultKey);
      if (!projection) throw new Error(`Bulk Feature style legend projection is missing for ${resultKey}.`);
      return projection;
    });
    let legend = null;
    let roots = null;
    const needsGeometry = change.geometryChangedResultKeys.length > 0;
    if (change.affectedResultKeys.length > 0 && !needsGeometry) {
      if (typeof prepareLegendProjection !== 'function') {
        throw new Error('Bulk Feature style changes require Result-local legend preparation.');
      }
      legend = await prepareLegendProjection({
        direction,
        from,
        to,
        catalog,
        results: sourceResults,
        affectedResultKeys: change.affectedResultKeys,
        projections,
        mounted: preflightMounted,
        manualLegendEntries,
        existingEntriesByResult,
        change
      });
      roots = preparedSvgMap(legend);
      change.affectedResultKeys.forEach((resultKey) => {
        if (!rootForResult(roots, resultKey)?.cloneNode) {
          throw new Error(`Bulk Feature style legend preparation is incomplete for ${resultKey}.`);
        }
      });
    }
    if (
      readRef(state.results) !== sourceResults
      || !resultsAreCurrent(readRef(state.results), sourceResultGuard)
      || !sameMount(currentMounted(getMountedContext), preflightMounted)
      || !sameSelectedResult(selectedResult(state, catalog), preflightSelectedResult)
    ) {
      throw new Error('Bulk Feature style state changed during legend preparation.');
    }
    let projection;
    if (needsGeometry) {
      if (typeof prepareGeometryProjection !== 'function') {
        throw new Error(
          `Bulk Feature style geometry changed in ${change.geometryChangedResultKeys.length} Result(s); complete worker preflight is required.`
        );
      }
      projection = await prepareGeometryProjection({
        direction,
        from,
        to,
        catalog,
        results: sourceResults,
        fillsByTargetKey: to.effectiveByTargetKey,
        affectedResultKeys: change.affectedResultKeys,
        geometryChangedResultKeys: change.geometryChangedResultKeys,
        projections,
        mounted: preflightMounted,
        targetFeatureKeysByResult: change.targetFeatureKeysByResult
      });
      legend = projection?.legend || null;
    } else if (change.affectedResultKeys.length > 0) {
      projection = await prepareResultProjection({
          direction,
          from,
          to,
          catalog,
          results: sourceResults,
          fillsByTargetKey: to.effectiveByTargetKey,
          affectedResultKeys: change.affectedResultKeys,
          mounted: preflightMounted,
          legend,
          preparedSvgByResultKey: roots,
          targetFeatureKeysByResult: change.targetFeatureKeysByResult
        });
    } else {
      projection = {
        previousResults: sourceResults,
        nextResults: [...sourceResults],
        admissionMetadataByResultKey: {},
        preparedMountedSvg: null,
        counters: {
          affectedResults: 0,
          mountedSerializations: 0,
          detachedPasses: 0,
          changedResults: 0
        }
      };
    }
    if (
      readRef(state.results) !== sourceResults
      || !resultsAreCurrent(readRef(state.results), sourceResultGuard)
      || !sameMount(currentMounted(getMountedContext), preflightMounted)
      || !sameSelectedResult(selectedResult(state, catalog), preflightSelectedResult)
    ) {
      throw new Error('Bulk Feature style state changed during Result preparation.');
    }
    if (!Array.isArray(projection?.nextResults) || projection.nextResults.length !== sourceResults.length) {
      throw new Error('Bulk Feature style Result preparation changed Result topology.');
    }
    if (projection.previousResults && projection.previousResults !== sourceResults) {
      throw new Error('Bulk Feature style Result preparation retained a stale source array.');
    }
    const affected = new Set(change.affectedResultKeys);
    const admission = projection.admissionMetadataByResultKey || {};
    projection.nextResults.forEach((result, resultIndex) => {
      const resultKey = change.resultKeys[resultIndex];
      const previous = sourceResults[resultIndex];
      if (!result || text(result.name) !== text(previous?.name)) {
        throw new Error(`Bulk Feature style Result identity changed at index ${resultIndex}.`);
      }
      if (!affected.has(resultKey) && result !== previous) {
        throw new Error(`Bulk Feature style projection replaced unaffected Result ${resultKey}.`);
      }
      if (affected.has(resultKey)) {
        if (result === previous || resultContent(result) === resultContent(previous)) {
          throw new Error(`Bulk Feature style projection did not change affected Result ${resultKey}.`);
        }
        if (!admission[resultKey]?.before || !admission[resultKey]?.after) {
          throw new Error(`Bulk Feature style admission metadata is incomplete for ${resultKey}.`);
        }
      }
    });
    const selectedProjection = legendForResult(to, preflightSelectedResult.resultKey);
    const selectedLegendEntries = affected.has(preflightSelectedResult.resultKey)
      ? cloneJson(
          projection?.selectedLegendEntries
          ?? (preflightSelectedResult.resultKey === preflightMounted.resultKey
            ? legend?.selectedEntries
            : null)
          ?? selectedProjection?.entries
          ?? []
        )
      : sourceLegendEntries;
    const affectedResults = change.affectedResultKeys.map((resultKey) => {
      const resultIndex = change.resultKeys.indexOf(resultKey);
      return {
        resultKey,
        before: sourceResults[resultIndex],
        beforeAdmission: admission[resultKey]?.before || null,
        after: projection.nextResults[resultIndex],
        afterAdmission: admission[resultKey]?.after || null
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
      sourceResultGuard,
      sourceLegendEntries,
      sourceLegendValue,
      selectedLegendEntries,
      preflightMounted,
      preflightSelectedResult,
      retainedForHistory,
      retainedBytes: exactJsonByteLength(retainedForHistory),
      counters: {
        projectedResults: Number(projection.counters?.affectedResults || 0),
        mountedSerializations: Number(projection.counters?.mountedSerializations || 0),
        detachedPasses: Number(projection.counters?.detachedPasses || 0),
        changedResults: Number(projection.counters?.changedResults || 0),
        resultArraySwaps: 1,
        legendMeasurements: Number(legend?.counters?.measurements || 0),
        geometryRenders: needsGeometry ? change.affectedResultKeys.length : 0
      }
    };
  };

  const commitTransition = ({ to, direction, prepared }) => {
    replaceRules(state.manualSpecificRules, to.rules);
    writeRef(state.appliedPaletteName, to.appliedPaletteName);
    writeRef(state.appliedPaletteColors, cloneJson(to.appliedPaletteColors));
    replaceAuxiliaryState(state, to.auxiliary);
    replaceObject(state.featureColorOverrides, to.overrides);
    if (
      change.affectedResultKeys.includes(prepared.preflightMounted.resultKey)
      && typeof commitMountedProjection === 'function'
    ) {
      const committed = commitMountedProjection({
        direction,
        prepared: prepared.projection,
        legend: prepared.legend,
        change
      });
      if (committed && typeof committed.then === 'function') {
        throw new Error('Bulk Feature style mounted commit must be synchronous.');
      }
      if (committed === false) return false;
    }
    if (change.affectedResultKeys.includes(prepared.preflightSelectedResult.resultKey)) {
      writeRef(state.legendEntries, cloneJson(prepared.selectedLegendEntries));
    }
    writeRef(state.results, prepared.projection.nextResults);
    if (!publishedResultsMatch(readRef(state.results), prepared.projection.nextResults)) {
      throw new Error('Bulk Feature style Result publication differs from its staged projection.');
    }
    prepared.publishedResults = readRef(state.results);
    writeRef(
      state.validatedStyleFingerprintByResultKey,
      Object.freeze(cloneJson(to.resultFingerprintByKey))
    );
    const reconciled = typeof reconcile === 'function'
      ? reconcile({ direction, prepared, change, to })
      : undefined;
    if (reconciled && typeof reconciled.then === 'function') {
      throw new Error('Bulk Feature style reconciliation must be synchronous.');
    }
    if (reconciled === false) return false;
    const refreshed = typeof refreshPresentation === 'function'
      ? refreshPresentation({ direction, prepared, change, to })
      : undefined;
    if (refreshed && typeof refreshed.then === 'function') {
      throw new Error('Bulk Feature style presentation refresh must be synchronous.');
    }
    if (refreshed === false) return false;
    writeRef(state.semanticStyleRevision, Number(readRef(state.semanticStyleRevision) || 0) + 1);
    writeRef(state.semanticStyleFingerprint, to.fingerprint);
    return true;
  };

  const restore = (snapshot, context) => {
    if (typeof restoreMountedProjection === 'function') {
      const restored = restoreMountedProjection({ snapshot, context, change });
      if (restored && typeof restored.then === 'function') {
        throw new Error('Bulk Feature style mounted rollback must be synchronous.');
      }
      if (restored === false) return false;
    }
    return restoreExactState({ state, snapshot, restoreRuntimeState });
  };

  return buildStyleSnapshotCommand({
    label,
    before: {
      ...commandBefore,
      resultFingerprintByKey: cloneJson(readRef(state.validatedStyleFingerprintByResultKey) || {})
    },
    after: commandAfter,
    metadata: {
      writerKind: text(writerKind) || 'bulk-style',
      documentEpoch: expected.documentEpoch,
      resultGenerationKey: expected.resultGenerationKey,
      affectedResultKeys: change.affectedResultKeys,
      geometryChangedResultKeys: change.geometryChangedResultKeys
    },
    counters: change.counters,
    isNoop: (from, to) => (
      sameJson(from.identity, to.identity)
      && sameJson(from.auxiliary, to.auxiliary)
    ),
    validateCurrent,
    prepareTransition,
    commitTransition,
    captureExactState: () => captureExactState({ state, captureRuntimeState }),
    restoreExactState: restore,
    assertCommitted: ({ to, direction, prepared }) => {
      const revision = Number(readRef(state.semanticStyleRevision) || 0);
      const valid = (
        sameJson(stateStyle(state), to.identity)
        && stateFingerprint(state) === to.fingerprint
        && text(readRef(state.semanticStyleFingerprint)) === to.fingerprint
        && sameJson(state.featureColorOverrides, to.overrides)
        && auxiliaryStateMatches(state, to.auxiliary)
        && ledgerMatches(
          readRef(state.validatedStyleFingerprintByResultKey) || {},
          change.resultKeys,
          to.fingerprint
        )
        && readRef(state.results) === prepared.publishedResults
        && publishedResultsMatch(readRef(state.results), prepared.projection.nextResults)
      );
      if (valid) {
        if (direction === 'undo') expectedRevision.redo = revision;
        else expectedRevision.undo = revision;
      }
      return valid;
    }
  });
};
