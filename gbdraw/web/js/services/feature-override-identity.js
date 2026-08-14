import { stableFeatureOverrideKey } from './feature-catalog.js';

const text = (value) => String(value ?? '').trim();
const hasOwn = (object, key) => Object.prototype.hasOwnProperty.call(object || {}, key);

const stableFeatureId = (feature) => text(
  feature?.stable_feature_id
  ?? feature?.stableFeatureId
  ?? feature?.stable_svg_id
  ?? feature?.stableSvgId
  ?? feature?.selector?.hash
);

const recordKeyForFeature = (feature) => text(
  feature?.record_key ?? feature?.recordKey
);

const recordIdForFeature = (feature) => text(
  feature?.record_id ?? feature?.recordId
);

const recordIndexForFeature = (feature) => {
  const value = feature?.fileIdx
    ?? feature?.file_idx
    ?? feature?.record_idx
    ?? feature?.recordIndex
    ?? feature?.record_index;
  const index = Number(value);
  return Number.isInteger(index) ? index : null;
};

const legacyFeatureAliases = (
  feature,
  index,
  { includePositional = false } = {}
) => {
  const candidates = [
    feature?.id,
    feature?.legacy_id,
    feature?.legacyId,
    feature?.svg_id,
    feature?.svgId,
    feature?.rendered_svg_id,
    feature?.renderedSvgId,
    feature?.rendered_feature_svg_id,
    feature?.renderedFeatureSvgId,
    feature?.stable_svg_id,
    feature?.stableSvgId,
    feature?.stable_feature_id,
    feature?.stableFeatureId,
    feature?.selector?.hash,
    includePositional && Number.isInteger(index) ? `f${index}` : ''
  ];
  return Array.from(new Set(candidates.map(text).filter(Boolean)));
};

export const featureOverrideKey = (feature) => (
  stableFeatureOverrideKey(feature)
  || text(feature?.stable_override_key)
  || text(feature?.id)
  || text(feature?.svg_id ?? feature?.svgId)
);

export const getFeatureOverride = (overrides, feature) => {
  const key = featureOverrideKey(feature);
  return key && hasOwn(overrides, key) ? overrides[key] : undefined;
};

export const migrateLegacyFeatureOverrides = (
  overrides,
  features,
  {
    warn = console.warn,
    legacyFeatures = [],
    legacyRecordKeys = [],
    onDiagnostic = null
  } = {}
) => {
  const result = {
    migrated: 0,
    unresolved: 0,
    ambiguous: 0,
    collisions: 0
  };
  if (!overrides || typeof overrides !== 'object') return result;

  const targets = new Set();
  const targetsByAlias = new Map();
  const currentFeatures = Array.isArray(features) ? features : [];
  const currentDescriptorsByStableId = new Map();
  const diagnostics = {
    currentDescriptorCount: 0,
    legacyFeatureCount: Array.isArray(legacyFeatures) ? legacyFeatures.length : 0,
    legacyFeaturesVisited: 0,
    legacyKeysNeedingMigration: 0,
    fullDescriptorComparisons: 0,
    indexedDescriptorComparisons: 0,
    skippedLegacyFeatureScan: false
  };
  const reportDiagnostic = () => {
    if (typeof onDiagnostic === 'function') onDiagnostic({ ...diagnostics });
  };
  const addAliasTarget = (alias, target) => {
    if (!alias || !target || alias === target) return;
    if (!targetsByAlias.has(alias)) targetsByAlias.set(alias, new Set());
    targetsByAlias.get(alias).add(target);
  };
  currentFeatures.forEach((feature, index) => {
    const target = stableFeatureOverrideKey(feature);
    if (!target) return;
    targets.add(target);
    legacyFeatureAliases(feature, index).forEach((alias) => {
      addAliasTarget(alias, target);
    });
    const descriptor = {
      target,
      stableId: stableFeatureId(feature),
      recordKey: recordKeyForFeature(feature),
      recordId: recordIdForFeature(feature)
    };
    diagnostics.currentDescriptorCount += 1;
    if (descriptor.stableId) {
      if (!currentDescriptorsByStableId.has(descriptor.stableId)) {
        currentDescriptorsByStableId.set(descriptor.stableId, []);
      }
      currentDescriptorsByStableId.get(descriptor.stableId).push(descriptor);
    }
  });
  if (targets.size === 0) {
    reportDiagnostic();
    return result;
  }

  const legacyKeys = Object.keys(overrides).filter((key) => !targets.has(key));
  const legacyKeySet = new Set(legacyKeys);
  diagnostics.legacyKeysNeedingMigration = legacyKeys.length;
  if (legacyKeys.length === 0) {
    diagnostics.skippedLegacyFeatureScan = true;
    reportDiagnostic();
    return result;
  }

  (Array.isArray(legacyFeatures) ? legacyFeatures : []).forEach((feature, index) => {
    const relevantAliases = legacyFeatureAliases(
      feature,
      index,
      { includePositional: true }
    ).filter((alias) => legacyKeySet.has(alias));
    if (relevantAliases.length === 0) return;
    diagnostics.legacyFeaturesVisited += 1;
    const sourceStableId = stableFeatureId(feature);
    if (!sourceStableId) return;
    let sourceRecordKey = recordKeyForFeature(feature);
    if (!sourceRecordKey) {
      const recordIndex = recordIndexForFeature(feature);
      if (recordIndex !== null) {
        if (Array.isArray(legacyRecordKeys)) {
          sourceRecordKey = text(legacyRecordKeys[recordIndex]);
        } else if (legacyRecordKeys && typeof legacyRecordKeys === 'object') {
          sourceRecordKey = text(legacyRecordKeys[recordIndex]);
        }
      }
    }
    const sourceRecordId = recordIdForFeature(feature);
    if (!sourceRecordKey && !sourceRecordId) return;
    let matches = currentDescriptorsByStableId.get(sourceStableId) || [];
    diagnostics.indexedDescriptorComparisons += matches.length;
    if (sourceRecordKey) {
      matches = matches.filter((candidate) => candidate.recordKey === sourceRecordKey);
    } else if (sourceRecordId) {
      matches = matches.filter((candidate) => candidate.recordId === sourceRecordId);
    }
    relevantAliases.forEach((alias) => {
      matches.forEach((candidate) => addAliasTarget(alias, candidate.target));
    });
  });

  Object.keys(overrides).forEach((legacyKey) => {
    if (targets.has(legacyKey)) return;
    const matches = targetsByAlias.get(legacyKey);
    if (!matches || matches.size === 0) {
      result.unresolved += 1;
      return;
    }
    if (matches.size > 1) {
      result.ambiguous += 1;
      return;
    }
    const [target] = matches;
    if (hasOwn(overrides, target)) {
      result.collisions += 1;
      return;
    }
    overrides[target] = overrides[legacyKey];
    delete overrides[legacyKey];
    result.migrated += 1;
  });

  const unresolvedCount = result.unresolved + result.ambiguous + result.collisions;
  if (unresolvedCount > 0 && typeof warn === 'function') {
    warn(
      `${unresolvedCount} legacy feature overrides remain unresolved `
      + `(${result.unresolved} unmatched, ${result.ambiguous} ambiguous, `
      + `${result.collisions} key collisions).`
    );
  }
  reportDiagnostic();
  return result;
};
