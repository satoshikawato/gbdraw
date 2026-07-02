import {
  buildLinearRegionExtractionContext,
  extractFeatureMetadataForPreview
} from './feature-metadata-extraction.js';

const READY_SMALL_RENDERED_COUNT = 20;
const READY_MATCH_RATIO = 0.98;
const READY_MAX_MISSING_COUNT = 20;

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

const addNormalizedId = (ids, value) => {
  const normalized = normalizeRenderedFeatureId(value);
  if (normalized) ids.add(normalized);
};

const renderedFeatureIdRegex = () => /data-gbdraw-feature-id\s*=\s*["']([^"']+)["']/g;

export const collectRenderedFeatureIdsFromSvg = (svgText) => {
  const ids = new Set();
  const text = String(svgText || '');
  if (!text) return ids;

  if (typeof DOMParser !== 'undefined') {
    try {
      const doc = new DOMParser().parseFromString(text, 'image/svg+xml');
      if (doc.querySelector('parsererror')) throw new Error('SVG parse error');
      doc.querySelectorAll('[data-gbdraw-feature-id]').forEach((element) => {
        addNormalizedId(ids, element.getAttribute('data-gbdraw-feature-id'));
      });
      return ids;
    } catch (_err) {
      // Fall back to the same attribute pattern when DOMParser rejects malformed SVG.
    }
  }

  for (const match of text.matchAll(renderedFeatureIdRegex())) {
    addNormalizedId(ids, match[1]);
  }
  return ids;
};

export const collectRenderedFeatureIdsFromResults = (results) => {
  const byResultIndex = [];
  const allIds = new Set();
  (Array.isArray(results) ? results : []).forEach((result, index) => {
    const ids = collectRenderedFeatureIdsFromSvg(result?.content || '');
    byResultIndex[index] = ids;
    ids.forEach((id) => allIds.add(id));
  });
  return { byResultIndex, allIds, totalRenderedCount: allIds.size };
};

const collectMetadataFeatureIds = (features) => {
  const ids = new Set();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    addNormalizedId(ids, feature?.svg_id);
    addNormalizedId(ids, feature?.stable_svg_id);
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
  const selectedIds = rendered.byResultIndex[safeIndex] || new Set();
  const metadataIds = collectMetadataFeatureIds(extractedFeatures);
  const renderedCount = selectedIds.size;
  const metadataCount = Array.isArray(extractedFeatures) ? extractedFeatures.length : 0;
  let matchingCount = 0;
  selectedIds.forEach((id) => {
    if (metadataIds.has(id)) matchingCount += 1;
  });
  const missingRenderedCount = Math.max(0, renderedCount - matchingCount);
  const matchRatio = renderedCount > 0 ? matchingCount / renderedCount : 1;

  let state = 'stale';
  if (rendered.totalRenderedCount === 0) {
    state = 'not-needed';
  } else if (renderedCount === 0) {
    state = 'ready';
  } else if (metadataCount === 0) {
    state = 'missing';
  } else if (
    (renderedCount <= READY_SMALL_RENDERED_COUNT && missingRenderedCount === 0) ||
    (
      renderedCount > READY_SMALL_RENDERED_COUNT &&
      matchRatio >= READY_MATCH_RATIO &&
      missingRenderedCount <= READY_MAX_MISSING_COUNT
    )
  ) {
    state = 'ready';
  }

  return {
    state,
    renderedCount,
    metadataCount,
    matchingCount,
    missingRenderedCount,
    matchRatio,
    selectedResultIndex: safeIndex,
    totalRenderedCount: rendered.totalRenderedCount
  };
};

const normalizeKey = (value) => String(value || '').trim();

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
  const key = normalizeKey(feature?.stable_svg_id);
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

export const buildFeatureOverrideMigration = ({
  previousFeatures,
  recoveredFeatures
}) => {
  const featureIdMap = {};
  const svgIdMap = {};
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
    if (oldId && newId) featureIdMap[oldId] = newId;
    if (oldSvgId && newSvgId) svgIdMap[oldSvgId] = newSvgId;
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

const cloneJson = (value, fallback) => {
  try {
    return JSON.parse(JSON.stringify(value));
  } catch (_err) {
    return fallback;
  }
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
    const value = cloneJson(rawValue, rawValue);

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
  const entry = orthogroupIndex instanceof Map
    ? (orthogroupIndex.get(orthogroupIndexKey(recordIndex, svgId)) || orthogroupIndex.get(svgId))
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

  if (validation.state === 'ready' || validation.state === 'not-needed') {
    return {
      status: 'ready',
      reason: validation.state,
      validation,
      errors: []
    };
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
    const recoveredValidation = classifyFeatureMetadataState({
      results: snapshot.results,
      selectedResultIndex: snapshot.selectedResultIndex,
      extractedFeatures: recoveredFeatureState.extractedFeatures
    });

    if (recoveredValidation.state !== 'ready' && recoveredValidation.state !== 'not-needed') {
      return {
        status: 'failed',
        reason: 'validation-not-ready-after-recovery',
        validation,
        recoveredValidation,
        warning: EXTRACTION_FAILED_WARNING,
        errors: payload.errors || []
      };
    }

    if (validation.state === 'stale') {
      const migration = buildFeatureOverrideMigration({
        previousFeatures: snapshot.featureState?.extractedFeatures || [],
        recoveredFeatures: recoveredFeatureState.extractedFeatures
      });
      const migrated = migrateFeatureOverrideState({
        featureState: recoveredFeatureState,
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
      recoveredFeatureState,
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
