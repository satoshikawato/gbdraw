// @ts-check
import {
  buildLinearRegionExtractionContext,
  extractFeatureMetadataForPreview
} from './feature-metadata-extraction.js';
import { normalizeCaughtError, retryCanSucceed } from '../utils/error-normalization.js';
import { cloneJsonValue } from './json-clone.js';
import { enrichFeatureWithOrthogroup } from './orthogroup-feature-metadata.js';
import {
  RECORD_INDEX_KEYS,
  RENDERED_FEATURE_ID_KEYS,
  STABLE_FEATURE_ID_KEYS,
  nonnegativeIntegerAliasStatus,
  textAliasStatus
} from './feature-identity.js';
import { getCommittedSvgResultMetadata } from './svg-result-ingestion.js';
import {
  addRenderedIdentity,
  createRenderedIdentityCollection,
  markRenderedIdentityAmbiguous,
  normalizeFeatureRecordIndex,
  normalizeRenderedFeatureId,
  stableRenderedFeatureRecordKey
} from './session-feature-metadata.js';

export { normalizeRenderedFeatureId } from './session-feature-metadata.js';

/** @import { ExtractedFeature, ExtractFeatureMetadataOptions, FeatureRecordInput, LinearRegionSeq } from './feature-metadata-extraction.js' */
/** @import { IntegerAliasStatus, TextAliasStatus } from './feature-identity.js' */
/** @import { OrthogroupFeatureIndex } from './orthogroup-feature-metadata.js' */
/** @import { RenderedFeatureIdentity, RenderedIdentityCollection } from './session-feature-metadata.js' */
/** @import { SvgResultEntry } from './svg-result-ingestion.js' */

/** @typedef {readonly (SvgResultEntry | null | undefined)[]} RecoveryResults The Results a Session holds. */
/**
 * What the recovery reads of a feature besides its identity keys.
 * @typedef {{ qualifiers?: Record<string, unknown> | null, selector?: { qualifiers?: Record<string, unknown> | null } | null, [key: string]: unknown }} QualifiedFeatureSource
 */
/** @typedef {Map<string, unknown> | Record<string, unknown>} MigrationMaps Old key to new key, as a Map or a plain object. */
/**
 * The feature state of a Session: the recovered metadata and the edits keyed by feature id.
 * @typedef {{
 *   extractedFeatures?: ExtractedFeature[],
 *   biologicalFeatures?: ExtractedFeature[],
 *   featureRecordIds?: string[],
 *   selectedFeatureRecordIdx?: number,
 *   featureColorOverrides?: unknown,
 *   [key: string]: unknown
 * }} RecoveryFeatureState
 */
/** @typedef {{ featureStrokes?: { overrides?: unknown }, [key: string]: unknown }} RecoveryEditorState */
/**
 * The Session data the recovery reads; the loader builds it from the candidate drawing and the saved Results.
 * @typedef {object} FeatureRecoverySnapshot
 * @property {string} mode
 * @property {string} cInputType
 * @property {string} lInputType
 * @property {{ c_gb?: File | null }} [files]
 * @property {(LinearRegionSeq & { gb?: unknown })[]} [linearSeqs]
 * @property {RecoveryResults} [results]
 * @property {unknown} [selectedResultIndex]
 * @property {RecoveryFeatureState} [featureState]
 * @property {RecoveryEditorState} [editorState]
 * @property {OrthogroupFeatureIndex | null} [orthogroupIndex]
 */
// Texts shown later (a Session Load's feature status): this module is copied alone by a
// test, so it declares the key marker of utils/i18n.js itself; the shower translates with `tm`.
/** @param {string} text @returns {string} */
const msg = (text) => text;
const MISSING_INPUTS_WARNING =
  msg('Feature metadata could not be recovered because the session does not include embedded GenBank inputs. Generate the diagram again or save a session with embedded inputs.');
const UNSUPPORTED_INPUT_WARNING =
  msg('Feature metadata recovery currently supports GenBank-backed sessions only. Generate the diagram again to enable feature editing.');
const EXTRACTION_FAILED_WARNING =
  msg('Feature metadata recovery failed. The SVG preview and pairwise popups remain available.');
const MIGRATION_SKIPPED_WARNING =
  msg('Feature metadata was recovered, but some old feature edits could not be matched to recovered feature IDs.');

/** @param {unknown} value */
const normalizeKey = (value) => String(value || '').trim();

export const normalizeRecordIndex = normalizeFeatureRecordIndex;

/**
 * @param {Set<string>} ids
 * @param {unknown} value
 */
const addNormalizedId = (ids, value) => {
  const normalized = normalizeRenderedFeatureId(value);
  if (normalized) ids.add(normalized);
};

export const stableRecordKey = stableRenderedFeatureRecordKey;

/** @param {RecoveryResults | null | undefined} results */
export const collectRenderedFeatureIdsFromResults = (results) => {
  /** @type {Set<string>[]} */
  const byResultIndex = [];
  /** @type {RenderedIdentityCollection[]} */
  const identitiesByResultIndex = [];
  const allIds = new Set();
  const allIdentities = createRenderedIdentityCollection();
  (Array.isArray(results) ? results : []).forEach((result, index) => {
    const identities = getCommittedSvgResultMetadata(result)?.renderedFeatureIdentities
      || createRenderedIdentityCollection();
    identitiesByResultIndex[index] = identities;
    byResultIndex[index] = identities.renderedIds;
    identities.renderedIds.forEach((renderedId) => {
      allIds.add(renderedId);
      if (identities.ambiguousRenderedIds.has(renderedId)) {
        markRenderedIdentityAmbiguous(allIdentities, renderedId);
        return;
      }
      const identity = identities.byRenderedId.get(renderedId);
      if (identity) addRenderedIdentity(allIdentities, identity);
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

/** @param {ExtractedFeature[] | null | undefined} features */
const collectExactMetadataFeatureIds = (features) => {
  /** @type {Set<string>} */
  const ids = new Set();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    addNormalizedId(ids, feature?.rendered_svg_id);
    addNormalizedId(ids, feature?.renderedSvgId);
    addNormalizedId(ids, feature?.rendered_feature_svg_id);
    addNormalizedId(ids, feature?.renderedFeatureSvgId);
    // Legacy sessions stored the rendered identity in svg_id.
    addNormalizedId(ids, feature?.svg_id);
  });
  return ids;
};

/** @param {ExtractedFeature | null | undefined} feature */
export const featureStableCandidate = (feature) =>
  normalizeRenderedFeatureId(
    feature?.stable_svg_id ||
    feature?.stableFeatureSvgId ||
    feature?.stable_feature_id ||
    feature?.stableFeatureId ||
    feature?.stable_feature_svg_id ||
    feature?.svg_id
  );

/** @param {ExtractedFeature | null | undefined} feature */
export const featureRenderedCandidate = (feature) =>
  normalizeRenderedFeatureId(
    feature?.rendered_svg_id ||
    feature?.renderedSvgId ||
    feature?.rendered_feature_svg_id ||
    feature?.renderedFeatureSvgId
  );

/**
 * A feature as the alignment reads it: `stableKey` and `recordIndex` find its
 * drawn element, `stable` and `recordIndex` must agree with that element, and
 * `drawnId` is the ID it is drawn under when the reader knows it.
 * @typedef {{ stableKey: string, stable: TextAliasStatus, recordIndex: IntegerAliasStatus, drawnId: string }} AlignmentReading
 */
/** @typedef {(feature: ExtractedFeature | null | undefined) => AlignmentReading} AlignmentReader */

/**
 * Saved feature metadata names its drawn element by its own keys.
 * @type {AlignmentReader}
 */
const readSavedFeature = (feature) => ({
  stableKey: featureStableCandidate(feature),
  stable: textAliasStatus(feature, STABLE_FEATURE_ID_KEYS),
  recordIndex: nonnegativeIntegerAliasStatus(feature, RECORD_INDEX_KEYS),
  drawnId: ''
});

const INPUT_INDEX_KEYS = Object.freeze(['fileIdx', 'file_idx']);
const RECORD_WITHIN_INPUT_KEYS = Object.freeze(RECORD_INDEX_KEYS.filter((key) => !INPUT_INDEX_KEYS.includes(key)));

/** @param {ExtractedFeature | null | undefined} feature */
const drawnHashOf = (feature) => {
  const selector = feature?.drawn_selector;
  return selector && typeof selector === 'object' && 'hash' in selector
    ? normalizeRenderedFeatureId(selector.hash)
    : '';
};

/**
 * How the features read again from the sources name their elements in the SVG
 * of a Session below 40 (OV-751). The SVG names a feature by its hash in the
 * drawn record, which the extraction gives as `drawn_selector.hash`: a cropped,
 * reverse-complemented or rotated record draws other hashes than its source.
 * Only drawn features carry it; another feature takes it from the drawn
 * feature of its source hash and record. The SVG numbers records among all
 * drawn records; the extraction names a Linear feature by its input
 * (`fileIdx`) and its record within that input (`record_idx`), and lists the
 * drawn records in order (`recordInputs`). Without that list (Circular), the
 * record index is the drawn one.
 * @param {{ recordInputs?: readonly FeatureRecordInput[] | null, drawnFeatures?: readonly ExtractedFeature[] | null }} payload
 * @returns {AlignmentReader}
 */
const extractedFeatureReader = ({ recordInputs, drawnFeatures }) => {
  /** @type {Map<string, number>} */
  const drawnIndexByInputRecord = new Map(
    (Array.isArray(recordInputs) ? recordInputs : [])
      .map(({ fileIdx, recordIdx }, drawnIndex) => /** @type {[string, number]} */ ([`${fileIdx}:${recordIdx}`, drawnIndex]))
  );
  /** @param {ExtractedFeature | null | undefined} feature @returns {IntegerAliasStatus} */
  const recordIndexOf = (feature) => {
    const input = nonnegativeIntegerAliasStatus(feature, INPUT_INDEX_KEYS);
    const record = nonnegativeIntegerAliasStatus(feature, RECORD_WITHIN_INPUT_KEYS);
    if (drawnIndexByInputRecord.size === 0 || !input.supplied || !record.supplied) {
      return nonnegativeIntegerAliasStatus(feature, RECORD_INDEX_KEYS);
    }
    const drawnIndex = input.valid && record.valid
      ? drawnIndexByInputRecord.get(`${input.value}:${record.value}`)
      : undefined;
    return drawnIndex === undefined
      ? { valid: false, supplied: true, value: null }
      : { valid: true, supplied: true, value: drawnIndex };
  };
  /** @param {IntegerAliasStatus} recordIndex @param {ExtractedFeature | null | undefined} feature */
  const sourceKey = (recordIndex, feature) => `${recordIndex.value}\u0000${featureStableCandidate(feature)}`;
  /** @type {Map<string, string>} */
  const drawnHashBySource = new Map();
  (Array.isArray(drawnFeatures) ? drawnFeatures : []).forEach((feature) => {
    const drawnHash = drawnHashOf(feature);
    const recordIndex = recordIndexOf(feature);
    if (drawnHash && recordIndex.valid) drawnHashBySource.set(sourceKey(recordIndex, feature), drawnHash);
  });
  return (feature) => {
    const recordIndex = recordIndexOf(feature);
    const drawnHash = drawnHashOf(feature)
      || (recordIndex.valid ? drawnHashBySource.get(sourceKey(recordIndex, feature)) : '')
      || '';
    if (!drawnHash) return { ...readSavedFeature(feature), recordIndex };
    return {
      stableKey: drawnHash,
      stable: { valid: true, supplied: true, value: drawnHash },
      recordIndex,
      drawnId: drawnHash
    };
  };
};

/**
 * Whether a feature can be the drawn element `identity`. An ID of the feature
 * that names no drawn element (an extraction's prediction, an older render)
 * says nothing about which element it is.
 * @param {ExtractedFeature | null | undefined} feature
 * @param {AlignmentReading} reading
 * @param {RenderedFeatureIdentity | null | undefined} identity
 * @param {ReadonlySet<string>} drawnIds
 */
const featureAgreesWithRenderedIdentity = (feature, reading, identity, drawnIds) => {
  if (!identity?.renderedId) return false;
  const { stable, recordIndex } = reading;
  const rendered = textAliasStatus(feature, RENDERED_FEATURE_ID_KEYS);
  const legacySvg = textAliasStatus(feature, ['svg_id', 'svgId']);
  const recordId = textAliasStatus(feature, ['record_id', 'recordId']);
  if (![stable, rendered, legacySvg, recordIndex, recordId].every((status) => status.valid)) {
    return false;
  }
  /** @param {TextAliasStatus} status */
  const namesDrawnElement = (status) => status.supplied && drawnIds.has(normalizeRenderedFeatureId(status.value));
  if (stable.supplied && stable.value !== identity.stableId) return false;
  if (namesDrawnElement(rendered) && rendered.value !== identity.renderedId) return false;
  if (
    namesDrawnElement(legacySvg) &&
    legacySvg.value !== identity.renderedId &&
    legacySvg.value !== identity.stableId
  ) return false;
  // An element drawn without its record (every Circular figure; a Linear
  // figure of a Session before per-element record attributes) contradicts no
  // record: its feature hash already covers the record ID (OV-751).
  if (
    recordIndex.supplied &&
    identity.recordIndex !== null &&
    recordIndex.value !== identity.recordIndex
  ) return false;
  if (recordId.supplied && identity.recordId && recordId.value !== identity.recordId) return false;
  return true;
};

/** @param {ExtractedFeature[] | null | undefined} features */
const collectAliasMetadataFeatureIds = (features) => {
  /** @type {Set<string>} */
  const ids = new Set();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    addNormalizedId(ids, featureStableCandidate(feature));
  });
  return ids;
};

/**
 * @param {{
 *   results?: RecoveryResults | null,
 *   selectedResultIndex?: unknown,
 *   extractedFeatures?: ExtractedFeature[] | null
 * }} input
 */
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
    if (selectedIdentities.ambiguousRenderedIds?.has(id)) return;
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

/** @param {ExtractedFeature | null | undefined} feature */
const featureId = (feature) => normalizeKey(feature?.id);
/** @param {ExtractedFeature | null | undefined} feature */
const featureSvgId = (feature) => normalizeRenderedFeatureId(feature?.svg_id);

/**
 * @param {ExtractedFeature} feature
 * @param {RenderedFeatureIdentity | null | undefined} identity
 * @param {unknown} stableId
 * @param {{ writeRenderedSvgIdToSvgId?: boolean }} [options]
 * @returns {ExtractedFeature}
 */
const withRenderedFeatureIdentity = (
  feature,
  identity,
  stableId,
  { writeRenderedSvgIdToSvgId = false } = {}
) => {
  const renderedId = normalizeRenderedFeatureId(identity?.renderedId);
  if (!renderedId) return feature;
  const resolvedStableId = normalizeRenderedFeatureId(stableId || identity?.stableId);
  /** @type {Record<string, unknown>} */
  const updates = {
    rendered_svg_id: renderedId,
    renderedSvgId: renderedId
  };
  // An extraction's predicted rendered ID is replaced by the drawn one.
  ['rendered_feature_svg_id', 'renderedFeatureSvgId'].forEach((key) => {
    if (Object.prototype.hasOwnProperty.call(feature || {}, key)) updates[key] = renderedId;
  });
  if (writeRenderedSvgIdToSvgId) updates.svg_id = renderedId;
  if (resolvedStableId) {
    updates.stable_svg_id = feature?.stable_svg_id || resolvedStableId;
    updates.stable_feature_id = feature?.stable_feature_id || resolvedStableId;
  }
  const changed = Object.entries(updates).some(([key, value]) => feature?.[key] !== value);
  return changed ? { ...feature, ...updates } : feature;
};

/** @param {ExtractedFeature} feature */
const withoutRenderedFeatureIdentity = (feature) => {
  const renderedKeys = [
    'rendered_svg_id',
    'renderedSvgId',
    'rendered_feature_svg_id',
    'renderedFeatureSvgId'
  ];
  if (!renderedKeys.some((key) => Object.prototype.hasOwnProperty.call(feature || {}, key))) return feature;
  const nextFeature = { ...feature };
  renderedKeys.forEach((key) => delete nextFeature[key]);
  return nextFeature;
};

/**
 * @param {...unknown} values
 * @returns {string}
 */
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

/**
 * @param {QualifiedFeatureSource | null | undefined} feature
 * @param {string} key
 */
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

/** @param {ExtractedFeature | null | undefined} feature */
const stableFeatureKey = (feature) => {
  const key = featureStableCandidate(feature);
  return key || '';
};

/** @param {ExtractedFeature | null | undefined} feature */
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

/**
 * @param {ExtractedFeature[] | null | undefined} features
 * @param {(feature: ExtractedFeature) => string} keyFn
 */
const uniqueFeatureMap = (features, keyFn) => {
  /** @type {Map<string, number>} */
  const counts = new Map();
  /** @type {Map<string, ExtractedFeature>} */
  const items = new Map();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const key = keyFn(feature);
    if (!key) return;
    counts.set(key, (counts.get(key) || 0) + 1);
    if (!items.has(key)) items.set(key, feature);
  });
  /** @type {Map<string, ExtractedFeature | undefined>} */
  const unique = new Map();
  counts.forEach((count, key) => {
    if (count === 1) unique.set(key, items.get(key));
  });
  return unique;
};

/** @param {readonly AlignmentReading[]} readings */
const countFeatureStableKeys = (readings) => {
  /** @type {Map<string, number>} */
  const counts = new Map();
  readings.forEach(({ stableKey }) => {
    if (!stableKey) return;
    counts.set(stableKey, (counts.get(stableKey) || 0) + 1);
  });
  return counts;
};

/** @param {RenderedIdentityCollection | null | undefined} renderedIdentities */
const identityCollectionOrEmpty = (renderedIdentities) => {
  if (renderedIdentities?.byRenderedId instanceof Map && renderedIdentities?.renderedIds instanceof Set) {
    return renderedIdentities;
  }
  return createRenderedIdentityCollection();
};

/**
 * @param {Record<string, string>} target
 * @param {unknown} fromKey
 * @param {unknown} toKey
 */
const plainMapWithEntry = (target, fromKey, toKey) => {
  const from = normalizeRenderedFeatureId(fromKey);
  const to = normalizeRenderedFeatureId(toKey);
  if (from && to && from !== to) target[from] = to;
};

/**
 * @param {{
 *   features?: ExtractedFeature[] | null,
 *   renderedIdentities?: RenderedIdentityCollection | null,
 *   writeRenderedSvgIdToSvgId?: boolean,
 *   readFeature?: AlignmentReader
 * }} input `readFeature` reads how a feature names its drawn element (default: saved metadata).
 */
export const alignRecoveredFeatureIdsToRenderedSvg = ({
  features,
  renderedIdentities,
  writeRenderedSvgIdToSvgId = false,
  readFeature = readSavedFeature
}) => {
  const sourceFeatures = Array.isArray(features) ? features : [];
  const rendered = identityCollectionOrEmpty(renderedIdentities);
  const readings = sourceFeatures.map((feature) => readFeature(feature));
  const stableFeatureCounts = countFeatureStableKeys(readings);
  /** @type {ExtractedFeature[]} */
  const nextFeatures = [];
  /** @type {Record<string, string>} */
  const svgIdMap = {};
  /** @type {Record<string, string>} */
  const featureIdMap = {};
  /** @type {{ index: number, status: string, fromSvgId: string, renderedId: string, stableId: string, recordIndex?: number | null, method: string }[]} */
  const resolutions = [];
  let exactCount = 0;
  let alignedCount = 0;
  let ambiguousCount = 0;
  let unresolvedCount = 0;
  let changedCount = 0;

  sourceFeatures.forEach((feature, index) => {
    const reading = readings[index];
    const oldSvgId = featureSvgId(feature);
    const explicitRenderedStatus = textAliasStatus(feature, RENDERED_FEATURE_ID_KEYS);
    const legacySvgStatus = textAliasStatus(feature, ['svg_id', 'svgId']);
    const explicitRenderedId = explicitRenderedStatus.valid ? explicitRenderedStatus.value : '';
    /** @param {string} id */
    const drawnOnce = (id) => Boolean(id) && rendered.renderedIds.has(id) && !rendered.ambiguousRenderedIds?.has(id);
    const exactRenderedId = (
      explicitRenderedStatus.valid && legacySvgStatus.valid && drawnOnce(explicitRenderedId)
        ? explicitRenderedId
        : (
          explicitRenderedStatus.valid && legacySvgStatus.valid &&
          !explicitRenderedStatus.supplied && drawnOnce(oldSvgId)
            ? oldSvgId
            : (drawnOnce(reading.drawnId) ? reading.drawnId : '')
        )
    );
    const exactIdentity = exactRenderedId ? rendered.byRenderedId.get(exactRenderedId) : null;
    if (exactIdentity && featureAgreesWithRenderedIdentity(feature, reading, exactIdentity, rendered.renderedIds)) {
      const identity = exactIdentity;
      const nextFeature = withRenderedFeatureIdentity(
        feature,
        identity,
        featureStableCandidate(feature),
        { writeRenderedSvgIdToSvgId }
      );
      exactCount += 1;
      if (writeRenderedSvgIdToSvgId) {
        plainMapWithEntry(svgIdMap, oldSvgId, exactRenderedId);
      }
      if (nextFeature !== feature) changedCount += 1;
      nextFeatures.push(nextFeature);
      resolutions.push({
        index,
        status: 'exact',
        fromSvgId: oldSvgId,
        renderedId: exactRenderedId,
        stableId: featureStableCandidate(feature),
        method: explicitRenderedId === exactRenderedId
          ? 'rendered_svg_id'
          : (oldSvgId === exactRenderedId ? 'legacy-svg_id' : 'drawn-hash')
      });
      return;
    }

    const stableId = featureStableCandidate(feature);
    const { stableKey, recordIndex: recordStatus } = reading;
    const recordIndex = recordStatus.valid && recordStatus.supplied ? recordStatus.value : null;
    const recordKey = stableRecordKey(stableKey, recordIndex);
    const recordMatches = recordKey ? rendered.byStableRecordKey.get(recordKey) || [] : [];
    const stableMatches = stableKey ? rendered.byStableId.get(stableKey) || [] : [];
    /** @type {RenderedFeatureIdentity | null} */
    let match = null;
    let method = '';

    if (recordMatches.length === 1) {
      match = recordMatches[0];
      method = 'stable-record';
    }

    if (match?.renderedId && featureAgreesWithRenderedIdentity(feature, reading, match, rendered.renderedIds)) {
      const oldId = featureId(feature);
      const newFeature = withRenderedFeatureIdentity(
        feature,
        match,
        stableId || match.stableId,
        { writeRenderedSvgIdToSvgId }
      );
      nextFeatures.push(newFeature);
      plainMapWithEntry(svgIdMap, oldSvgId || stableId, match.renderedId);
      const newId = featureId(newFeature);
      if (oldId && newId && oldId !== newId) featureIdMap[oldId] = newId;
      alignedCount += 1;
      if (newFeature !== feature) changedCount += 1;
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
      stableKey &&
      (
        recordMatches.length > 1 ||
        stableMatches.length > 1 ||
        (stableMatches.length === 1 && (stableFeatureCounts.get(stableKey) || 0) > 1)
      )
    );
    if (isAmbiguous) ambiguousCount += 1;
    else unresolvedCount += 1;
    const unresolvedFeature = withoutRenderedFeatureIdentity(feature);
    if (unresolvedFeature !== feature) changedCount += 1;
    nextFeatures.push(unresolvedFeature);
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

  // Features that differ only by a record the drawn element does not carry
  // all agree with it (OV-751): when the features claiming one element differ
  // by record, or one's record is unknown, none is aligned. Twins of one record
  // (one hash drawn twice) stay aligned.
  /** @param {number} index */
  const claimRecord = (index) => {
    const { valid, supplied, value } = readings[index].recordIndex;
    return valid && supplied && value !== null ? String(value) : '';
  };
  /** @type {Map<string, string[]>} */
  const claims = new Map();
  resolutions.forEach(({ renderedId }, index) => {
    if (renderedId) claims.set(renderedId, [...(claims.get(renderedId) || []), claimRecord(index)]);
  });
  /** @param {string} renderedId */
  const claimedAcrossRecords = (renderedId) => {
    const records = claims.get(renderedId) || [];
    return records.length > 1 && (new Set(records).size > 1 || records.includes(''));
  };
  resolutions.forEach((resolution, index) => {
    const { renderedId, status } = resolution;
    if (!renderedId || !claimedAcrossRecords(renderedId)) return;
    const feature = sourceFeatures[index];
    if (status === 'exact') exactCount -= 1;
    else alignedCount -= 1;
    ambiguousCount += 1;
    if (nextFeatures[index] !== feature) changedCount -= 1;
    nextFeatures[index] = withoutRenderedFeatureIdentity(feature);
    if (nextFeatures[index] !== feature) changedCount += 1;
    Object.keys(svgIdMap).forEach((from) => { if (svgIdMap[from] === renderedId) delete svgIdMap[from]; });
    delete featureIdMap[featureId(feature)];
    resolutions[index] = { ...resolution, status: 'ambiguous', renderedId: '', method: '' };
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

/** @param {unknown} source */
const normalizeMigrationMap = (source) => {
  /** @type {Record<string, string>} */
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

/**
 * @param {{
 *   previousFeatures?: ExtractedFeature[] | null,
 *   recoveredFeatures?: ExtractedFeature[] | null,
 *   featureIdMap?: MigrationMaps | null,
 *   svgIdMap?: MigrationMaps | null
 * }} input
 */
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

  /**
   * @param {ExtractedFeature | undefined} oldFeature
   * @param {ExtractedFeature | undefined} newFeature
   */
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

  /** @param {(feature: ExtractedFeature) => string} keyFn */
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

/** @param {ExtractedFeature[] | null | undefined} features */
const recoveredKeySets = (features) => {
  /** @type {Set<string>} */
  const featureIds = new Set();
  /** @type {Set<string>} */
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

/**
 * @param {MigrationMaps | null | undefined} mapLike
 * @param {unknown} key
 */
const lookupMapValue = (mapLike, key) => {
  const normalized = normalizeKey(key);
  if (!normalized) return '';
  if (mapLike instanceof Map) return normalizeKey(mapLike.get(normalized));
  if (mapLike && typeof mapLike === 'object') return normalizeKey(mapLike[normalized]);
  return '';
};

/**
 * @param {unknown} source
 * @param {ReadonlySet<string>} preserveKeys
 * @param {MigrationMaps[]} mappingLookups
 * @param {{ count: number }} skipped
 */
const rewriteOverrideMap = (source, preserveKeys, mappingLookups, skipped) => {
  /** @type {Record<string, unknown>} */
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

/**
 * @param {{
 *   featureState?: RecoveryFeatureState | null,
 *   editorState?: RecoveryEditorState | null,
 *   migration?: { featureIdMap?: MigrationMaps, svgIdMap?: MigrationMaps } | null
 * }} input
 */
export const migrateFeatureOverrideState = ({
  featureState,
  editorState,
  migration
}) => {
  /** @type {RecoveryFeatureState} */
  const nextFeatureState = cloneJsonValue(featureState || {}, {});
  /** @type {RecoveryEditorState} */
  const nextEditorState = cloneJsonValue(editorState || {}, {});
  const keySets = recoveredKeySets(nextFeatureState.extractedFeatures);
  const skipped = { count: 0 };
  const featureIdMaps = [migration?.featureIdMap || {}];
  const allMaps = [migration?.featureIdMap || {}, migration?.svgIdMap || {}];

  nextFeatureState.featureColorOverrides = rewriteOverrideMap(
    nextFeatureState.featureColorOverrides,
    keySets.featureIds,
    featureIdMaps,
    skipped
  );
  // Rendered-ID Feature visibility and label edits are not realigned here:
  // services/feature-edit-migration.js moves the saved maps onto source
  // identities (design Q4 4.3).

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

/** @param {Partial<FeatureRecoverySnapshot>} [snapshot] */
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

/**
 * @param {FeatureRecoverySnapshot} snapshot
 * @param {{ extractedFeatures?: ExtractedFeature[], biologicalFeatures?: ExtractedFeature[],
 *   featureRecordIds?: string[], selectedFeatureRecordIdx?: number }} payload What the extraction returned.
 */
const buildRecoveredFeatureState = (snapshot, payload) => ({
  ...(snapshot.featureState || {}),
  extractedFeatures: Array.isArray(payload.extractedFeatures) ? payload.extractedFeatures : [],
  biologicalFeatures: Array.isArray(payload.biologicalFeatures)
    ? payload.biologicalFeatures
    : (Array.isArray(payload.extractedFeatures) ? payload.extractedFeatures : []),
  featureRecordIds: Array.isArray(payload.featureRecordIds) ? payload.featureRecordIds : [],
  selectedFeatureRecordIdx: Number.isInteger(payload.selectedFeatureRecordIdx)
    ? payload.selectedFeatureRecordIdx
    : 0
});

/**
 * @param {FeatureRecoverySnapshot | null | undefined} snapshot
 * @param {number} selectedResultIndex
 */
const selectedRenderedIdentities = (snapshot, selectedResultIndex) => {
  const rendered = collectRenderedFeatureIdsFromResults(snapshot?.results);
  return rendered.identitiesByResultIndex[selectedResultIndex] || createRenderedIdentityCollection();
};

/**
 * @param {{
 *   snapshot: FeatureRecoverySnapshot | null | undefined,
 *   featureState: RecoveryFeatureState | null | undefined,
 *   selectedResultIndex: number,
 *   readFeature?: AlignmentReader
 * }} input
 */
const alignFeatureStateToRenderedSvg = ({ snapshot, featureState, selectedResultIndex, readFeature }) => {
  const renderedIdentities = selectedRenderedIdentities(snapshot, selectedResultIndex);
  const alignment = alignRecoveredFeatureIdsToRenderedSvg({
    features: featureState?.extractedFeatures || [],
    renderedIdentities,
    writeRenderedSvgIdToSvgId: true,
    readFeature
  });
  const biologicalAlignment = alignRecoveredFeatureIdsToRenderedSvg({
    features: featureState?.biologicalFeatures || featureState?.extractedFeatures || [],
    renderedIdentities,
    readFeature
  });
  const renderedFeatures = alignment.features
    .filter((feature) => Boolean(featureRenderedCandidate(feature)));
  const alignedFeatureState = {
    ...(featureState || {}),
    extractedFeatures: renderedFeatures,
    biologicalFeatures: biologicalAlignment.features
  };
  const alignedValidation = classifyFeatureMetadataState({
    results: snapshot?.results,
    selectedResultIndex: snapshot?.selectedResultIndex,
    extractedFeatures: alignedFeatureState.extractedFeatures
  });
  return {
    alignment,
    biologicalAlignment,
    alignedFeatureState,
    alignedValidation,
    featureStateChanged: (
      alignment.changedCount > 0 ||
      biologicalAlignment.changedCount > 0 ||
      renderedFeatures.length !== (Array.isArray(featureState?.extractedFeatures)
        ? featureState.extractedFeatures.length
        : 0)
    )
  };
};

/** @param {{ state?: string } | null | undefined} validation */
const isClickReadyValidation = (validation) =>
  validation?.state === 'ready' || validation?.state === 'not-needed';

/** @param {RecoveryFeatureState | null | undefined} featureState */
export const hasUsableBiologicalFeatureCatalog = (featureState) => {
  const features = featureState?.biologicalFeatures;
  if (!Array.isArray(features) || features.length === 0) return false;
  return features.some((feature) => Boolean(
    normalizeKey(feature?.nucleotide_sequence || feature?.nucleotideSequence) ||
    normalizeKey(feature?.amino_acid_sequence || feature?.aminoAcidSequence)
  ));
};

/** @param {{ featureIdMap?: MigrationMaps, svgIdMap?: MigrationMaps } | null | undefined} migration */
const hasMigrationEntries = (migration) =>
  Object.keys(migration?.featureIdMap || {}).length > 0 ||
  Object.keys(migration?.svgIdMap || {}).length > 0;

// Reads the snapshot's GenBank sources with its crops, record selectors, and
// orientations, as the renderer drew them.
/**
 * @param {FeatureRecoverySnapshot} snapshot
 * @param {{
 *   featureVisibilityTsv?: string,
 *   enrichFeature?: ExtractFeatureMetadataOptions['enrichFeature'],
 *   readFeatureExtractionDataImpl?: ExtractFeatureMetadataOptions['readFeatureExtractionDataImpl'],
 *   extractFeatureMetadataForPreviewImpl?: typeof extractFeatureMetadataForPreview
 * }} options
 */
const extractSnapshotFeatures = (snapshot, {
  featureVisibilityTsv = '',
  enrichFeature,
  readFeatureExtractionDataImpl,
  extractFeatureMetadataForPreviewImpl = extractFeatureMetadataForPreview
}) => {
  const linearContext = snapshot.mode === 'linear' && snapshot.lInputType === 'gb'
    ? buildLinearRegionExtractionContext(snapshot.linearSeqs, snapshot.lInputType)
    : { regionSpecs: [], recordSelectors: [], reverseFlags: [] };
  return extractFeatureMetadataForPreviewImpl({
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
    ...(enrichFeature ? { enrichFeature } : {}),
    readFeatureExtractionDataImpl
  });
};

/**
 * Every source feature of a Session without a feature catalog, read again with
 * its crops and orientations and without visibility rules, so each feature
 * carries its drawn hash beside its source identity. The rendered-ID edit
 * migration (services/feature-edit-migration.js) names features through it.
 * Returns null when the Session has no GenBank sources to read or a source
 * cannot be read as saved; the migration then uses the saved feature metadata.
 * A runtime failure (the Worker does not start, transport, cancellation) is
 * thrown: the Load fails and keeps the previous Session instead of dropping
 * the edits as unmatched (OV-39).
 */
/**
 * @param {{
 *   snapshot: FeatureRecoverySnapshot,
 *   readFeatureExtractionDataImpl?: ExtractFeatureMetadataOptions['readFeatureExtractionDataImpl'],
 *   extractFeatureMetadataForPreviewImpl?: typeof extractFeatureMetadataForPreview
 * }} input
 */
export const extractSessionSourceFeatures = async ({
  snapshot,
  readFeatureExtractionDataImpl,
  extractFeatureMetadataForPreviewImpl = extractFeatureMetadataForPreview
}) => {
  if (!recoverability(snapshot).recoverable) return null;
  try {
    const payload = await extractSnapshotFeatures(snapshot, {
      readFeatureExtractionDataImpl,
      extractFeatureMetadataForPreviewImpl
    });
    return payload?.errors?.length || !Array.isArray(payload?.extractedFeatures) ? null : payload;
  } catch (error) {
    if (retryCanSucceed(normalizeCaughtError(error))) throw error;
    return null;
  }
};

/**
 * @param {{
 *   snapshot: FeatureRecoverySnapshot,
 *   featureVisibilityTsv?: string,
 *   readFeatureExtractionDataImpl?: ExtractFeatureMetadataOptions['readFeatureExtractionDataImpl'],
 *   extractFeatureMetadataForPreviewImpl?: typeof extractFeatureMetadataForPreview
 * }} input
 */
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
  const biologicalCatalogReady = hasUsableBiologicalFeatureCatalog(snapshot?.featureState);

  if (validation.state === 'not-needed') {
    return {
      status: 'ready',
      reason: 'not-needed',
      validation,
      errors: []
    };
  }

  if (isClickReadyValidation(validation) && validation.renderedCount > 0) {
    const existingAlignment = alignFeatureStateToRenderedSvg({
      snapshot,
      featureState: snapshot?.featureState || {},
      selectedResultIndex: validation.selectedResultIndex
    });
    // Saved metadata that matched every drawn feature is not replaced by an
    // alignment that drops some of them (OV-752).
    if (existingAlignment.featureStateChanged && isClickReadyValidation(existingAlignment.alignedValidation)) {
      const migrated = migrateFeatureOverrideState({
        featureState: existingAlignment.alignedFeatureState,
        editorState: snapshot?.editorState,
        migration: existingAlignment.alignment
      });
      return {
        status: 'aligned',
        reason: 'rendered-identities-indexed',
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

  if (biologicalCatalogReady && isClickReadyValidation(validation)) {
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
    if (biologicalCatalogReady && isClickReadyValidation(existingAlignment.alignedValidation)) {
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
    if (isClickReadyValidation(validation)) {
      return {
        status: 'ready',
        reason: 'biological-catalog-unavailable',
        validation,
        warning: recovery.warning,
        errors: []
      };
    }
    return {
      status: 'unrecoverable',
      reason: recovery.reason,
      validation,
      warning: recovery.warning,
      errors: []
    };
  }

  try {
    const payload = await extractSnapshotFeatures(snapshot, {
      featureVisibilityTsv,
      enrichFeature: (feature, recordIndex) => enrichFeatureWithOrthogroup(
        snapshot.orthogroupIndex,
        feature,
        recordIndex
      ),
      readFeatureExtractionDataImpl,
      extractFeatureMetadataForPreviewImpl
    });

    const recoveredFeatureState = buildRecoveredFeatureState(snapshot, payload);
    const recoveredAlignment = alignFeatureStateToRenderedSvg({
      snapshot,
      featureState: recoveredFeatureState,
      selectedResultIndex: validation.selectedResultIndex,
      readFeature: extractedFeatureReader({
        recordInputs: payload.featureRecordInputs,
        drawnFeatures: recoveredFeatureState.extractedFeatures
      })
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
