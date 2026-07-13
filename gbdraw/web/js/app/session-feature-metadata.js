import {
  buildLinearRegionExtractionContext,
  extractFeatureMetadataForPreview
} from './feature-metadata-extraction.js';
import { cloneJsonValue } from '../services/json-clone.js';

const MISSING_INPUTS_WARNING =
  'Feature metadata could not be recovered because the session does not include embedded GenBank inputs. Generate the diagram again or save a session with embedded inputs.';
const UNSUPPORTED_INPUT_WARNING =
  'Feature metadata recovery currently supports GenBank-backed sessions only. Generate the diagram again to enable feature editing.';
const EXTRACTION_FAILED_WARNING =
  'Feature metadata recovery failed. The SVG preview and pairwise popups remain available.';
const MIGRATION_SKIPPED_WARNING =
  'Feature metadata was recovered, but some old feature edits could not be matched to recovered feature IDs.';

export const normalizeRenderedFeatureId = (value) =>
  String(value || '').trim().replace(/__part\d+$/, '');

const normalizeKey = (value) => String(value || '').trim();

export const normalizeRecordIndex = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  return Number.isInteger(numeric) && numeric >= 0 ? numeric : null;
};

const addNormalizedId = (ids, value) => {
  const normalized = normalizeRenderedFeatureId(value);
  if (normalized) ids.add(normalized);
};

const renderedFeatureIdRegex = () => /data-gbdraw-feature-id\s*=\s*["']([^"']+)["']/g;

const createRenderedIdentityCollection = () => ({
  byRenderedId: new Map(),
  renderedIds: new Set(),
  byStableId: new Map(),
  byStableRecordKey: new Map(),
  totalRenderedCount: 0
});

export const stableRecordKey = (stableId, recordIndex) => {
  const normalizedStableId = normalizeRenderedFeatureId(stableId);
  const normalizedRecordIndex = normalizeRecordIndex(recordIndex);
  return normalizedStableId && normalizedRecordIndex !== null
    ? `${normalizedStableId}\u001f${normalizedRecordIndex}`
    : '';
};

const addToIdentityListMap = (target, key, identity) => {
  if (!key) return;
  const items = target.get(key) || [];
  items.push(identity);
  target.set(key, items);
};

const addRenderedIdentity = (collection, identity) => {
  const renderedId = normalizeRenderedFeatureId(identity?.renderedId);
  if (!renderedId || collection.byRenderedId.has(renderedId)) return;
  const stableId = normalizeRenderedFeatureId(identity?.stableId) || renderedId;
  const recordIndex = normalizeRecordIndex(identity?.recordIndex);
  const normalized = {
    renderedId,
    stableId,
    recordIndex,
    recordId: normalizeKey(identity?.recordId),
    elementId: normalizeKey(identity?.elementId) || renderedId
  };
  collection.byRenderedId.set(renderedId, normalized);
  collection.renderedIds.add(renderedId);
  addToIdentityListMap(collection.byStableId, stableId, normalized);
  addToIdentityListMap(collection.byStableRecordKey, stableRecordKey(stableId, recordIndex), normalized);
  collection.totalRenderedCount = collection.renderedIds.size;
};

const elementRenderedId = (element) =>
  normalizeRenderedFeatureId(element?.getAttribute?.('data-gbdraw-feature-id')) ||
  normalizeRenderedFeatureId(element?.getAttribute?.('id'));

const elementStableId = (element, renderedId) =>
  normalizeRenderedFeatureId(element?.getAttribute?.('data-gbdraw-stable-feature-id')) ||
  renderedId;

const collectRenderedFeatureElements = (doc) => {
  const selector = [
    '[data-gbdraw-feature-id]',
    'path[id^="f"]',
    'polygon[id^="f"]',
    'rect[id^="f"]'
  ].join(', ');
  return Array.from(doc.querySelectorAll(selector) || []);
};

export const collectRenderedFeatureIdentitiesFromSvg = (svgText) => {
  const collection = createRenderedIdentityCollection();
  const text = String(svgText || '');
  if (!text) return collection;

  if (typeof DOMParser !== 'undefined') {
    try {
      const doc = new DOMParser().parseFromString(text, 'image/svg+xml');
      if (doc.querySelector('parsererror')) throw new Error('SVG parse error');
      collectRenderedFeatureElements(doc).forEach((element) => {
        const renderedId = elementRenderedId(element);
        if (!renderedId) return;
        addRenderedIdentity(collection, {
          renderedId,
          stableId: elementStableId(element, renderedId),
          recordIndex: element.getAttribute('data-gbdraw-record-index'),
          recordId: element.getAttribute('data-gbdraw-record-id'),
          elementId: element.getAttribute('id') || renderedId
        });
      });
      return collection;
    } catch (_err) {
      // Regex compatibility is intentionally ID-only for malformed SVG.
    }
  }

  for (const match of text.matchAll(renderedFeatureIdRegex())) {
    const renderedId = normalizeRenderedFeatureId(match[1]);
    addRenderedIdentity(collection, {
      renderedId,
      stableId: renderedId,
      recordIndex: null,
      recordId: '',
      elementId: renderedId
    });
  }
  return collection;
};

export const collectRenderedFeatureIdsFromSvg = (svgText) => {
  const identities = collectRenderedFeatureIdentitiesFromSvg(svgText);
  return identities.renderedIds;
};

export const collectRenderedFeatureIdsFromResults = (results) => {
  const byResultIndex = [];
  const identitiesByResultIndex = [];
  const allIds = new Set();
  const allIdentities = createRenderedIdentityCollection();
  (Array.isArray(results) ? results : []).forEach((result, index) => {
    const identities = collectRenderedFeatureIdentitiesFromSvg(result?.content || '');
    identitiesByResultIndex[index] = identities;
    byResultIndex[index] = identities.renderedIds;
    identities.byRenderedId.forEach((identity) => {
      allIds.add(identity.renderedId);
      addRenderedIdentity(allIdentities, identity);
    });
  });
  return {
    byResultIndex,
    identitiesByResultIndex,
    renderedIdentitiesByResultIndex: identitiesByResultIndex,
    allIds,
    allIdentities,
    totalRenderedCount: allIds.size
  };
};

const collectExactMetadataFeatureIds = (features) => {
  const ids = new Set();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    addNormalizedId(ids, feature?.svg_id);
  });
  return ids;
};

export const featureStableCandidate = (feature) =>
  normalizeRenderedFeatureId(
    feature?.stable_svg_id ||
    feature?.stableFeatureSvgId ||
    feature?.stable_feature_id ||
    feature?.svg_id
  );

export const featureRecordIndexCandidate = (feature) => {
  for (const value of [feature?.fileIdx, feature?.record_idx, feature?.recordIndex, feature?.record_index]) {
    const normalized = normalizeRecordIndex(value);
    if (normalized !== null) return normalized;
  }
  return null;
};

const collectAliasMetadataFeatureIds = (features) => {
  const ids = new Set();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    addNormalizedId(ids, featureStableCandidate(feature));
  });
  return ids;
};

export const classifyFeatureMetadataState = ({
  results,
  selectedResultIndex,
  extractedFeatures
}) => {
  const rendered = collectRenderedFeatureIdsFromResults(results);
  const resultCount = Array.isArray(results) ? results.length : 0;
  const safeIndex = resultCount > 0
    ? Math.max(0, Math.min(Number(selectedResultIndex) || 0, resultCount - 1))
    : 0;
  const selectedIdentities = rendered.identitiesByResultIndex[safeIndex] || createRenderedIdentityCollection();
  const selectedIds = selectedIdentities.renderedIds || new Set();
  const exactMetadataIds = collectExactMetadataFeatureIds(extractedFeatures);
  const aliasMetadataIds = collectAliasMetadataFeatureIds(extractedFeatures);
  const renderedCount = selectedIds.size;
  const metadataCount = Array.isArray(extractedFeatures) ? extractedFeatures.length : 0;
  let exactMatchingCount = 0;
  let aliasMatchingCount = 0;
  selectedIds.forEach((id) => {
    if (exactMetadataIds.has(id)) {
      exactMatchingCount += 1;
      return;
    }
    const identity = selectedIdentities.byRenderedId?.get(id);
    if (identity?.stableId && aliasMetadataIds.has(identity.stableId)) aliasMatchingCount += 1;
  });
  const missingExactCount = Math.max(0, renderedCount - exactMatchingCount);
  const matchRatio = renderedCount > 0 ? exactMatchingCount / renderedCount : 1;

  let state = 'stale';
  if (rendered.totalRenderedCount === 0) {
    state = 'not-needed';
  } else if (renderedCount === 0) {
    state = 'ready';
  } else if (metadataCount === 0) {
    state = 'missing';
  } else if (missingExactCount === 0) {
    state = 'ready';
  } else if (aliasMatchingCount > 0) {
    state = 'alignable';
  }

  return {
    state,
    renderedCount,
    metadataCount,
    exactMatchingCount,
    aliasMatchingCount,
    missingExactCount,
    matchingCount: exactMatchingCount,
    missingRenderedCount: missingExactCount,
    matchRatio,
    selectedResultIndex: safeIndex,
    totalRenderedCount: rendered.totalRenderedCount
  };
};

const featureId = (feature) => normalizeKey(feature?.id);
const featureSvgId = (feature) => normalizeRenderedFeatureId(feature?.svg_id);

const firstPresentValue = (...values) => {
  for (const value of values) {
    if (Array.isArray(value)) {
      const nested = firstPresentValue(...value);
      if (nested) return nested;
      continue;
    }
    const normalized = normalizeKey(value);
    if (normalized) return normalized;
  }
  return '';
};

const qualifierValues = (feature, key) => {
  const qualifiers = feature?.qualifiers && typeof feature.qualifiers === 'object' ? feature.qualifiers : {};
  const selectorQualifiers =
    feature?.selector?.qualifiers && typeof feature.selector.qualifiers === 'object'
      ? feature.selector.qualifiers
      : {};
  const raw = qualifiers[key] ?? selectorQualifiers[key] ?? feature?.[key];
  if (Array.isArray(raw)) {
    return raw.map((value) => normalizeKey(value)).filter(Boolean).join('\u001f');
  }
  return normalizeKey(raw);
};

const stableFeatureKey = (feature) => {
  const key = featureStableCandidate(feature);
  return key || '';
};

const biologicalFeatureKey = (feature) => {
  const type = normalizeKey(feature?.type);
  const start = feature?.start;
  const end = feature?.end;
  const strand = normalizeKey(feature?.strand);
  if (!type || start === null || start === undefined || end === null || end === undefined || !strand) return '';
  return JSON.stringify({
    fileIdx: firstPresentValue(feature?.fileIdx, feature?.file_idx),
    record_idx: firstPresentValue(feature?.record_idx, feature?.recordIndex),
    record_id: normalizeKey(feature?.record_id),
    type,
    start: String(start),
    end: String(end),
    strand,
    locus_tag: qualifierValues(feature, 'locus_tag'),
    gene: qualifierValues(feature, 'gene'),
    protein_id: qualifierValues(feature, 'protein_id')
  });
};

const uniqueFeatureMap = (features, keyFn) => {
  const counts = new Map();
  const items = new Map();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const key = keyFn(feature);
    if (!key) return;
    counts.set(key, (counts.get(key) || 0) + 1);
    if (!items.has(key)) items.set(key, feature);
  });
  const unique = new Map();
  counts.forEach((count, key) => {
    if (count === 1) unique.set(key, items.get(key));
  });
  return unique;
};

const countFeatureStableCandidates = (features) => {
  const counts = new Map();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const stableId = featureStableCandidate(feature);
    if (!stableId) return;
    counts.set(stableId, (counts.get(stableId) || 0) + 1);
  });
  return counts;
};

const identityCollectionOrEmpty = (renderedIdentities) => {
  if (renderedIdentities?.byRenderedId instanceof Map && renderedIdentities?.renderedIds instanceof Set) {
    return renderedIdentities;
  }
  return createRenderedIdentityCollection();
};

const plainMapWithEntry = (target, fromKey, toKey) => {
  const from = normalizeRenderedFeatureId(fromKey);
  const to = normalizeRenderedFeatureId(toKey);
  if (from && to && from !== to) target[from] = to;
};

export const alignRecoveredFeatureIdsToRenderedSvg = ({
  features,
  renderedIdentities
}) => {
  const sourceFeatures = Array.isArray(features) ? features : [];
  const rendered = identityCollectionOrEmpty(renderedIdentities);
  const stableFeatureCounts = countFeatureStableCandidates(sourceFeatures);
  const nextFeatures = [];
  const svgIdMap = {};
  const featureIdMap = {};
  const resolutions = [];
  let exactCount = 0;
  let alignedCount = 0;
  let ambiguousCount = 0;
  let unresolvedCount = 0;
  let changedCount = 0;

  sourceFeatures.forEach((feature, index) => {
    const oldSvgId = featureSvgId(feature);
    if (oldSvgId && rendered.renderedIds.has(oldSvgId)) {
      exactCount += 1;
      nextFeatures.push(feature);
      resolutions.push({
        index,
        status: 'exact',
        fromSvgId: oldSvgId,
        renderedId: oldSvgId,
        stableId: featureStableCandidate(feature),
        method: 'svg_id'
      });
      return;
    }

    const stableId = featureStableCandidate(feature);
    const recordIndex = featureRecordIndexCandidate(feature);
    const recordKey = stableRecordKey(stableId, recordIndex);
    const recordMatches = recordKey ? rendered.byStableRecordKey.get(recordKey) || [] : [];
    const stableMatches = stableId ? rendered.byStableId.get(stableId) || [] : [];
    let match = null;
    let method = '';

    if (recordMatches.length === 1) {
      match = recordMatches[0];
      method = 'stable-record';
    } else if (
      stableId &&
      (stableFeatureCounts.get(stableId) || 0) === 1 &&
      stableMatches.length === 1
    ) {
      match = stableMatches[0];
      method = 'unique-stable';
    }

    if (match?.renderedId) {
      const oldId = featureId(feature);
      const newFeature = {
        ...feature,
        svg_id: match.renderedId,
        stable_svg_id: stableId || match.stableId
      };
      nextFeatures.push(newFeature);
      plainMapWithEntry(svgIdMap, oldSvgId || stableId, match.renderedId);
      const newId = featureId(newFeature);
      if (oldId && newId && oldId !== newId) featureIdMap[oldId] = newId;
      alignedCount += 1;
      changedCount += 1;
      resolutions.push({
        index,
        status: 'aligned',
        fromSvgId: oldSvgId,
        renderedId: match.renderedId,
        stableId: stableId || match.stableId,
        recordIndex,
        method
      });
      return;
    }

    const isAmbiguous = Boolean(
      stableId &&
      (
        recordMatches.length > 1 ||
        stableMatches.length > 1 ||
        (stableMatches.length === 1 && (stableFeatureCounts.get(stableId) || 0) > 1)
      )
    );
    if (isAmbiguous) ambiguousCount += 1;
    else unresolvedCount += 1;
    nextFeatures.push(feature);
    resolutions.push({
      index,
      status: isAmbiguous ? 'ambiguous' : 'unresolved',
      fromSvgId: oldSvgId,
      renderedId: '',
      stableId,
      recordIndex,
      method: ''
    });
  });

  return {
    features: nextFeatures,
    changedCount,
    exactCount,
    alignedCount,
    ambiguousCount,
    unresolvedCount,
    svgIdMap,
    featureIdMap,
    resolutions
  };
};

const normalizeMigrationMap = (source) => {
  const normalized = {};
  if (source instanceof Map) {
    source.forEach((value, key) => {
      const from = normalizeKey(key);
      const to = normalizeKey(value);
      if (from && to) normalized[from] = to;
    });
    return normalized;
  }
  if (!source || typeof source !== 'object' || Array.isArray(source)) return normalized;
  Object.entries(source).forEach(([key, value]) => {
    const from = normalizeKey(key);
    const to = normalizeKey(value);
    if (from && to) normalized[from] = to;
  });
  return normalized;
};

export const buildFeatureOverrideMigration = ({
  previousFeatures,
  recoveredFeatures,
  featureIdMap: explicitFeatureIdMap = {},
  svgIdMap: explicitSvgIdMap = {}
}) => {
  const featureIdMap = normalizeMigrationMap(explicitFeatureIdMap);
  const svgIdMap = normalizeMigrationMap(explicitSvgIdMap);
  const mappedOldFeatures = new Set();
  const mappedNewFeatures = new Set();

  const addFeatureMapping = (oldFeature, newFeature) => {
    if (!oldFeature || !newFeature || mappedOldFeatures.has(oldFeature) || mappedNewFeatures.has(newFeature)) return;
    mappedOldFeatures.add(oldFeature);
    mappedNewFeatures.add(newFeature);

    const oldId = featureId(oldFeature);
    const newId = featureId(newFeature);
    const oldSvgId = featureSvgId(oldFeature);
    const newSvgId = featureSvgId(newFeature);
    if (oldId && newId && !Object.prototype.hasOwnProperty.call(featureIdMap, oldId)) featureIdMap[oldId] = newId;
    if (oldSvgId && newSvgId && !Object.prototype.hasOwnProperty.call(svgIdMap, oldSvgId)) svgIdMap[oldSvgId] = newSvgId;
  };

  const addUniqueMatches = (keyFn) => {
    const oldByKey = uniqueFeatureMap(previousFeatures, keyFn);
    const newByKey = uniqueFeatureMap(recoveredFeatures, keyFn);
    oldByKey.forEach((oldFeature, key) => {
      const newFeature = newByKey.get(key);
      if (newFeature) addFeatureMapping(oldFeature, newFeature);
    });
  };

  addUniqueMatches(stableFeatureKey);
  addUniqueMatches(biologicalFeatureKey);

  return {
    featureIdMap,
    svgIdMap,
    migratedCount: mappedOldFeatures.size,
    skippedCount: Math.max(0, (Array.isArray(previousFeatures) ? previousFeatures.length : 0) - mappedOldFeatures.size)
  };
};

const recoveredKeySets = (features) => {
  const featureIds = new Set();
  const svgIds = new Set();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const id = featureId(feature);
    const svgId = featureSvgId(feature);
    if (id) featureIds.add(id);
    if (svgId) svgIds.add(svgId);
  });
  const all = new Set([...featureIds, ...svgIds]);
  return { featureIds, svgIds, all };
};

const lookupMapValue = (mapLike, key) => {
  const normalized = normalizeKey(key);
  if (!normalized) return '';
  if (mapLike instanceof Map) return normalizeKey(mapLike.get(normalized));
  if (mapLike && typeof mapLike === 'object') return normalizeKey(mapLike[normalized]);
  return '';
};

const rewriteOverrideMap = (source, preserveKeys, mappingLookups, skipped) => {
  const rewritten = {};
  if (!source || typeof source !== 'object' || Array.isArray(source)) return rewritten;

  Object.entries(source).forEach(([rawKey, rawValue]) => {
    const key = normalizeKey(rawKey);
    if (!key) return;
    const value = cloneJsonValue(rawValue, rawValue);

    if (preserveKeys.has(key)) {
      rewritten[key] = value;
      return;
    }

    const mappedKey = mappingLookups.map((mapping) => lookupMapValue(mapping, key)).find(Boolean);
    if (mappedKey && preserveKeys.has(mappedKey) && !Object.prototype.hasOwnProperty.call(rewritten, mappedKey)) {
      rewritten[mappedKey] = value;
      return;
    }

    skipped.count += 1;
  });

  return rewritten;
};

export const migrateFeatureOverrideState = ({
  featureState,
  editorState,
  migration
}) => {
  const nextFeatureState = cloneJson(featureState || {}, {});
  const nextEditorState = cloneJson(editorState || {}, {});
  const keySets = recoveredKeySets(nextFeatureState.extractedFeatures);
  const skipped = { count: 0 };
  const featureIdMaps = [migration?.featureIdMap || {}];
  const svgIdMaps = [migration?.svgIdMap || {}];
  const allMaps = [migration?.featureIdMap || {}, migration?.svgIdMap || {}];

  nextFeatureState.featureColorOverrides = rewriteOverrideMap(
    nextFeatureState.featureColorOverrides,
    keySets.featureIds,
    featureIdMaps,
    skipped
  );
  nextFeatureState.featureVisibilityOverrides = rewriteOverrideMap(
    nextFeatureState.featureVisibilityOverrides,
    keySets.svgIds,
    svgIdMaps,
    skipped
  );
  nextFeatureState.labelTextFeatureOverrides = rewriteOverrideMap(
    nextFeatureState.labelTextFeatureOverrides,
    keySets.all,
    allMaps,
    skipped
  );
  nextFeatureState.labelTextFeatureOverrideSources = rewriteOverrideMap(
    nextFeatureState.labelTextFeatureOverrideSources,
    keySets.all,
    allMaps,
    skipped
  );
  nextFeatureState.labelVisibilityOverrides = rewriteOverrideMap(
    nextFeatureState.labelVisibilityOverrides,
    keySets.all,
    allMaps,
    skipped
  );

  if (!nextEditorState.featureStrokes || typeof nextEditorState.featureStrokes !== 'object') {
    nextEditorState.featureStrokes = { overrides: {} };
  }
  nextEditorState.featureStrokes.overrides = rewriteOverrideMap(
    nextEditorState.featureStrokes.overrides,
    keySets.all,
    allMaps,
    skipped
  );

  const warnings = skipped.count > 0 ? [MIGRATION_SKIPPED_WARNING] : [];
  return {
    featureState: nextFeatureState,
    editorState: nextEditorState,
    warnings,
    skippedOverrideCount: skipped.count
  };
};

const orthogroupIndexKey = (recordIndex, svgId) => `${Number(recordIndex)}:${String(svgId || '').trim()}`;

const enrichFeatureWithOrthogroup = (orthogroupIndex, feature, recordIndex) => {
  const svgId = String(feature?.svg_id || '').trim();
  if (!svgId) return feature;
  const stableSvgId = featureStableCandidate(feature);
  const entry = orthogroupIndex instanceof Map
    ? (
        orthogroupIndex.get(orthogroupIndexKey(recordIndex, svgId)) ||
        orthogroupIndex.get(svgId) ||
        orthogroupIndex.get(orthogroupIndexKey(recordIndex, stableSvgId)) ||
        orthogroupIndex.get(stableSvgId)
      )
    : null;
  if (!entry) return feature;
  return {
    ...feature,
    proteinId: entry.proteinId,
    sourceProteinId: entry.sourceProteinId,
    orthogroupId: entry.orthogroupId,
    orthogroupMemberCount: entry.orthogroupMemberCount,
    orthogroupRecordCoverage: entry.orthogroupRecordCoverage,
    orthogroupRepresentative: entry.orthogroupRepresentative,
    orthogroupScope: entry.orthogroupScope,
    orthogroupSourceRecordIndex: entry.orthogroupSourceRecordIndex,
    orthogroupMember: entry.orthogroupMember
  };
};

const recoverability = (snapshot = {}) => {
  const { mode, cInputType, lInputType } = snapshot;
  if (mode === 'circular' && cInputType === 'gb') {
    return snapshot.files?.c_gb
      ? { recoverable: true, reason: 'circular-gb' }
      : { recoverable: false, reason: 'missing-inputs', warning: MISSING_INPUTS_WARNING };
  }
  if (mode === 'linear' && lInputType === 'gb') {
    const seqs = Array.isArray(snapshot.linearSeqs) ? snapshot.linearSeqs : [];
    const hasAllGenbankInputs = seqs.length > 0 && seqs.every((seq) => Boolean(seq?.gb));
    return hasAllGenbankInputs
      ? { recoverable: true, reason: 'linear-gb' }
      : { recoverable: false, reason: 'missing-inputs', warning: MISSING_INPUTS_WARNING };
  }
  return { recoverable: false, reason: 'unsupported-input-type', warning: UNSUPPORTED_INPUT_WARNING };
};

const buildRecoveredFeatureState = (snapshot, payload) => ({
  ...(snapshot.featureState || {}),
  extractedFeatures: Array.isArray(payload.extractedFeatures) ? payload.extractedFeatures : [],
  featureSelectorSafetyScope: Array.isArray(payload.featureSelectorSafetyScope)
    ? payload.featureSelectorSafetyScope
    : [],
  featureRecordIds: Array.isArray(payload.featureRecordIds) ? payload.featureRecordIds : [],
  selectedFeatureRecordIdx: Number.isInteger(payload.selectedFeatureRecordIdx)
    ? payload.selectedFeatureRecordIdx
    : 0
});

const selectedRenderedIdentities = (snapshot, selectedResultIndex) => {
  const rendered = collectRenderedFeatureIdsFromResults(snapshot?.results);
  return rendered.identitiesByResultIndex[selectedResultIndex] || createRenderedIdentityCollection();
};

const alignFeatureStateToRenderedSvg = ({ snapshot, featureState, selectedResultIndex }) => {
  const alignment = alignRecoveredFeatureIdsToRenderedSvg({
    features: featureState?.extractedFeatures || [],
    renderedIdentities: selectedRenderedIdentities(snapshot, selectedResultIndex)
  });
  const alignedFeatureState = {
    ...(featureState || {}),
    extractedFeatures: alignment.features
  };
  const alignedValidation = classifyFeatureMetadataState({
    results: snapshot?.results,
    selectedResultIndex: snapshot?.selectedResultIndex,
    extractedFeatures: alignedFeatureState.extractedFeatures
  });
  return { alignment, alignedFeatureState, alignedValidation };
};

const isClickReadyValidation = (validation) =>
  validation?.state === 'ready' || validation?.state === 'not-needed';

const hasMigrationEntries = (migration) =>
  Object.keys(migration?.featureIdMap || {}).length > 0 ||
  Object.keys(migration?.svgIdMap || {}).length > 0;

export const buildSessionFeatureRecoveryPlan = async ({
  snapshot,
  featureVisibilityTsv,
  readFeatureExtractionDataImpl,
  extractFeatureMetadataForPreviewImpl = extractFeatureMetadataForPreview
}) => {
  const validation = classifyFeatureMetadataState({
    results: snapshot?.results,
    selectedResultIndex: snapshot?.selectedResultIndex,
    extractedFeatures: snapshot?.featureState?.extractedFeatures
  });

  if (isClickReadyValidation(validation)) {
    return {
      status: 'ready',
      reason: validation.state,
      validation,
      errors: []
    };
  }

  if (validation.state === 'alignable') {
    const existingAlignment = alignFeatureStateToRenderedSvg({
      snapshot,
      featureState: snapshot?.featureState || {},
      selectedResultIndex: validation.selectedResultIndex
    });
    if (isClickReadyValidation(existingAlignment.alignedValidation)) {
      const migrated = migrateFeatureOverrideState({
        featureState: existingAlignment.alignedFeatureState,
        editorState: snapshot?.editorState,
        migration: existingAlignment.alignment
      });
      return {
        status: 'aligned',
        reason: 'stale-metadata-aligned',
        validation,
        recoveredValidation: existingAlignment.alignedValidation,
        alignment: existingAlignment.alignment,
        recoveredFeatureState: migrated.featureState,
        migratedEditorState: migrated.editorState,
        warning: migrated.warnings[0] || null,
        errors: []
      };
    }
  }

  const recovery = recoverability(snapshot);
  if (!recovery.recoverable) {
    return {
      status: 'unrecoverable',
      reason: recovery.reason,
      validation,
      warning: recovery.warning,
      errors: []
    };
  }

  try {
    const linearContext = snapshot.mode === 'linear' && snapshot.lInputType === 'gb'
      ? buildLinearRegionExtractionContext(snapshot.linearSeqs, snapshot.lInputType)
      : { regionSpecs: [], recordSelectors: [], reverseFlags: [] };
    const payload = await extractFeatureMetadataForPreviewImpl({
      mode: snapshot.mode,
      cInputType: snapshot.cInputType,
      lInputType: snapshot.lInputType,
      circularFile: snapshot.files?.c_gb || null,
      linearSeqs: snapshot.linearSeqs || [],
      regionSpecs: linearContext.regionSpecs,
      recordSelectors: linearContext.recordSelectors,
      reverseFlags: linearContext.reverseFlags,
      featureVisibilityTablePath: featureVisibilityTsv ? '/web_feature_visibility_table.tsv' : null,
      featureVisibilityTsv,
      enrichFeature: (feature, recordIndex) => enrichFeatureWithOrthogroup(
        snapshot.orthogroupIndex,
        feature,
        recordIndex
      ),
      readFeatureExtractionDataImpl
    });

    const recoveredFeatureState = buildRecoveredFeatureState(snapshot, payload);
    const recoveredAlignment = alignFeatureStateToRenderedSvg({
      snapshot,
      featureState: recoveredFeatureState,
      selectedResultIndex: validation.selectedResultIndex
    });
    const alignedRecoveredFeatureState = recoveredAlignment.alignedFeatureState;
    const recoveredValidation = recoveredAlignment.alignedValidation;

    if (!isClickReadyValidation(recoveredValidation)) {
      return {
        status: 'failed',
        reason: 'validation-not-ready-after-recovery',
        validation,
        recoveredValidation,
        alignment: recoveredAlignment.alignment,
        warning: EXTRACTION_FAILED_WARNING,
        errors: payload.errors || []
      };
    }

    const migration = buildFeatureOverrideMigration({
      previousFeatures: snapshot.featureState?.extractedFeatures || [],
      recoveredFeatures: alignedRecoveredFeatureState.extractedFeatures,
      featureIdMap: recoveredAlignment.alignment.featureIdMap,
      svgIdMap: recoveredAlignment.alignment.svgIdMap
    });

    if (validation.state === 'stale' || validation.state === 'alignable' || hasMigrationEntries(migration)) {
      const migrated = migrateFeatureOverrideState({
        featureState: alignedRecoveredFeatureState,
        editorState: snapshot.editorState,
        migration
      });
      return {
        status: 'recovered',
        reason: 'stale-metadata-recovered',
        validation,
        recoveredValidation,
        recoveredFeatureState: migrated.featureState,
        migratedEditorState: migrated.editorState,
        warning: migrated.warnings[0] || null,
        errors: payload.errors || []
      };
    }

    return {
      status: 'recovered',
      reason: 'missing-metadata-recovered',
      validation,
      recoveredValidation,
      alignment: recoveredAlignment.alignment,
      recoveredFeatureState: alignedRecoveredFeatureState,
      migratedEditorState: null,
      warning: null,
      errors: payload.errors || []
    };
  } catch (error) {
    return {
      status: 'failed',
      reason: 'extraction-failed',
      validation,
      warning: EXTRACTION_FAILED_WARNING,
      errors: [error]
    };
  }
};
