// @ts-check
/** @import { FeatureCatalogItem, FeatureCatalogRow } from './feature-catalog.js' */

/**
 * One annotation warning of the Result metadata, as Python writes it.
 * @typedef {object} AnnotationWarning
 * @property {string} code
 * @property {string} setId
 * @property {string} annotationId
 * @property {string} recordId
 * @property {number} recordIndex
 * @property {number} missingCount
 * @property {string} message
 * @property {number} resultIndex
 * @property {string} resultName
 */

/**
 * A rendered feature as the SVG draws it: the rendered id, the stable
 * (source) id it belongs to, and the record it is drawn in.
 * @typedef {object} RenderedFeatureIdentity
 * @property {string} renderedId
 * @property {string} stableId
 * @property {number | null} recordIndex
 * @property {string} recordId
 * @property {string} elementId
 */

/**
 * The rendered identities of one Result, indexed for lookup.
 * @typedef {object} RenderedIdentityCollection
 * @property {Map<string, RenderedFeatureIdentity>} byRenderedId
 * @property {Set<string>} renderedIds
 * @property {Set<string>} ambiguousRenderedIds
 * @property {Map<string, RenderedFeatureIdentity[]>} byStableId
 * @property {Map<string, RenderedFeatureIdentity[]>} byStableRecordKey
 * @property {number} totalRenderedCount
 */

/**
 * An identity before normalization; every field may be missing or a raw attribute value.
 * @typedef {{ renderedId?: unknown, stableId?: unknown, recordIndex?: unknown, recordId?: unknown, elementId?: unknown }} RenderedIdentityInput
 */

/**
 * @param {unknown} warnings
 * @param {ReadonlyArray<{ name?: unknown } | undefined> | null | undefined} results The Results the warnings name.
 * @returns {AnnotationWarning[]}
 */
export const validateAnnotationWarnings = (warnings, results) => {
  if (warnings === undefined) return [];
  const fields = ['code', 'setId', 'annotationId', 'recordId', 'recordIndex',
    'missingCount', 'message', 'resultIndex', 'resultName'];
  const codes = new Set(['feature_selector_unmatched', 'empty_span',
    'out_of_bounds_skipped', 'out_of_bounds_clipped']);
  if (!Array.isArray(warnings) || warnings.some((warning) => (
    (warning === null || typeof warning !== 'object' || Array.isArray(warning))
    || Object.keys(warning).length !== fields.length
    || fields.some((field) => !Object.prototype.hasOwnProperty.call(warning, field))
    || !codes.has(warning.code)
    || ['setId', 'annotationId', 'recordId', 'message', 'resultName']
      .some((field) => typeof warning[field] !== 'string')
    || !warning.setId || !warning.annotationId
    || ['recordIndex', 'missingCount', 'resultIndex']
      .some((field) => !Number.isSafeInteger(warning[field]) || warning[field] < 0)
    || (warning.code === 'feature_selector_unmatched') !== (warning.missingCount > 0)
    || results?.[warning.resultIndex]?.name !== warning.resultName
  ))) throw new Error('Annotation warnings do not match the successful Result metadata schema.');
  return warnings;
};

const FEATURE_ID_SELECTOR = [
  '[data-gbdraw-feature-id]',
  'path[id^="f"]',
  'polygon[id^="f"]',
  'rect[id^="f"]'
].join(', ');

/** @param {unknown} value */
export const normalizeRenderedFeatureId = (value) =>
  String(value || '').trim().replace(/__part\d+$/, '');

/**
 * @param {unknown} value
 * @returns {number | null}
 */
export const normalizeFeatureRecordIndex = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const text = String(value).trim();
  if (typeof value === 'boolean' || !/^\d+$/.test(text)) return null;
  const numeric = Number(text);
  return Number.isSafeInteger(numeric) ? numeric : null;
};

/**
 * @param {unknown} stableId
 * @param {unknown} recordIndex
 */
export const stableRenderedFeatureRecordKey = (stableId, recordIndex) => {
  const normalizedStableId = normalizeRenderedFeatureId(stableId);
  const normalizedRecordIndex = normalizeFeatureRecordIndex(recordIndex);
  return normalizedStableId && normalizedRecordIndex !== null
    ? `${normalizedStableId}\u001f${normalizedRecordIndex}`
    : '';
};

/** @returns {RenderedIdentityCollection} */
export const createRenderedIdentityCollection = () => ({
  byRenderedId: new Map(),
  renderedIds: new Set(),
  ambiguousRenderedIds: new Set(),
  byStableId: new Map(),
  byStableRecordKey: new Map(),
  totalRenderedCount: 0
});

/**
 * @param {Map<string, RenderedFeatureIdentity[]>} target
 * @param {string} key
 * @param {RenderedFeatureIdentity} identity
 */
const addToIdentityListMap = (target, key, identity) => {
  if (!key) return;
  const items = target.get(key) || [];
  items.push(identity);
  target.set(key, items);
};

/**
 * @param {Map<string, RenderedFeatureIdentity[]>} target
 * @param {string} key
 * @param {RenderedFeatureIdentity} identity
 */
const removeFromIdentityListMap = (target, key, identity) => {
  if (!key) return;
  const remaining = (target.get(key) || []).filter((item) => item !== identity);
  if (remaining.length > 0) target.set(key, remaining);
  else target.delete(key);
};

/**
 * @param {RenderedIdentityCollection} collection
 * @param {unknown} renderedId
 */
export const markRenderedIdentityAmbiguous = (collection, renderedId) => {
  const normalizedRenderedId = normalizeRenderedFeatureId(renderedId);
  if (!normalizedRenderedId) return;
  const existing = collection.byRenderedId.get(normalizedRenderedId);
  if (existing) {
    collection.byRenderedId.delete(normalizedRenderedId);
    removeFromIdentityListMap(collection.byStableId, existing.stableId, existing);
    removeFromIdentityListMap(
      collection.byStableRecordKey,
      stableRenderedFeatureRecordKey(existing.stableId, existing.recordIndex),
      existing
    );
  }
  collection.renderedIds.add(normalizedRenderedId);
  collection.ambiguousRenderedIds.add(normalizedRenderedId);
  collection.totalRenderedCount = collection.renderedIds.size;
};

/**
 * @param {RenderedIdentityCollection} collection
 * @param {RenderedIdentityInput | null | undefined} identity
 */
export const addRenderedIdentity = (collection, identity) => {
  const renderedId = normalizeRenderedFeatureId(identity?.renderedId);
  if (!renderedId || collection.ambiguousRenderedIds.has(renderedId)) return;
  const stableId = normalizeRenderedFeatureId(identity?.stableId) || renderedId;
  const recordIndex = normalizeFeatureRecordIndex(identity?.recordIndex);
  const normalized = {
    renderedId,
    stableId,
    recordIndex,
    recordId: String(identity?.recordId || '').trim(),
    elementId: String(identity?.elementId || '').trim() || renderedId
  };
  const existing = collection.byRenderedId.get(renderedId);
  if (existing) {
    const agrees = (
      existing.stableId === normalized.stableId
      && existing.recordIndex === normalized.recordIndex
      && existing.recordId === normalized.recordId
    );
    if (!agrees) markRenderedIdentityAmbiguous(collection, renderedId);
    return;
  }
  collection.byRenderedId.set(renderedId, normalized);
  collection.renderedIds.add(renderedId);
  addToIdentityListMap(collection.byStableId, stableId, normalized);
  addToIdentityListMap(
    collection.byStableRecordKey,
    stableRenderedFeatureRecordKey(stableId, recordIndex),
    normalized
  );
  collection.totalRenderedCount = collection.renderedIds.size;
};

/**
 * @param {FeatureCatalogItem | null | undefined} item
 * @returns {RenderedIdentityCollection}
 */
export const collectRenderedFeatureIdentitiesFromCatalogItem = (item) => {
  const collection = createRenderedIdentityCollection();
  const recordIndexByKey = new Map(
    (Array.isArray(item?.recordKeys) ? item.recordKeys : []).map(
      (recordKey, recordIndex) => /** @type {[string, number]} */ ([String(recordKey || '').trim(), recordIndex])
    )
  );
  /** @type {Map<string, FeatureCatalogRow>} */
  const biologicalByKey = new Map();
  (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).forEach(
    (feature) => {
      const recordKey = String(feature?.recordKey || '').trim();
      const featureId = String(feature?.biologicalFeatureId || '').trim();
      if (recordKey && featureId) biologicalByKey.set(`${recordKey}\u0000${featureId}`, feature);
    }
  );
  (Array.isArray(item?.features) ? item.features : []).forEach((feature) => {
    const recordKey = String(feature?.recordKey || '').trim();
    const featureId = String(feature?.biologicalFeatureId || '').trim();
    /** @type {FeatureCatalogRow} */
    const biological = biologicalByKey.get(`${recordKey}\u0000${featureId}`) || {};
    const renderedId = String(feature?.svgId || '').trim();
    addRenderedIdentity(collection, {
      renderedId,
      stableId: biological.stableFeatureId || featureId || renderedId,
      recordIndex: recordIndexByKey.get(recordKey),
      recordId: biological.record_id || biological.recordId || recordKey,
      elementId: renderedId
    });
  });
  return collection;
};

/**
 * Build feature identity metadata from the sanitized, detached ingress root.
 * This function never parses SVG text; the ingestion owner supplies the one root.
 * @param {ParentNode | null | undefined} svg
 * @returns {RenderedIdentityCollection}
 */
export const collectRenderedFeatureIdentitiesFromSvgRoot = (svg) => {
  const collection = createRenderedIdentityCollection();
  Array.from(svg?.querySelectorAll?.(FEATURE_ID_SELECTOR) || []).forEach((element) => {
    const renderedId = normalizeRenderedFeatureId(
      element?.getAttribute?.('data-gbdraw-feature-id')
      || element?.getAttribute?.('id')
    );
    if (!renderedId) return;
    addRenderedIdentity(collection, {
      renderedId,
      stableId: element.getAttribute('data-gbdraw-stable-feature-id') || renderedId,
      recordIndex: element.getAttribute('data-gbdraw-record-index'),
      recordId: element.getAttribute('data-gbdraw-record-id'),
      elementId: element.getAttribute('id') || renderedId
    });
  });
  return collection;
};
