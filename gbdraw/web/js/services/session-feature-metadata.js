const FEATURE_ID_SELECTOR = [
  '[data-gbdraw-feature-id]',
  'path[id^="f"]',
  'polygon[id^="f"]',
  'rect[id^="f"]'
].join(', ');

export const normalizeRenderedFeatureId = (value) =>
  String(value || '').trim().replace(/__part\d+$/, '');

export const normalizeFeatureRecordIndex = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const text = String(value).trim();
  if (typeof value === 'boolean' || !/^\d+$/.test(text)) return null;
  const numeric = Number(text);
  return Number.isSafeInteger(numeric) ? numeric : null;
};

export const stableRenderedFeatureRecordKey = (stableId, recordIndex) => {
  const normalizedStableId = normalizeRenderedFeatureId(stableId);
  const normalizedRecordIndex = normalizeFeatureRecordIndex(recordIndex);
  return normalizedStableId && normalizedRecordIndex !== null
    ? `${normalizedStableId}\u001f${normalizedRecordIndex}`
    : '';
};

export const createRenderedIdentityCollection = () => ({
  byRenderedId: new Map(),
  renderedIds: new Set(),
  ambiguousRenderedIds: new Set(),
  byStableId: new Map(),
  byStableRecordKey: new Map(),
  totalRenderedCount: 0
});

const addToIdentityListMap = (target, key, identity) => {
  if (!key) return;
  const items = target.get(key) || [];
  items.push(identity);
  target.set(key, items);
};

const removeFromIdentityListMap = (target, key, identity) => {
  if (!key) return;
  const remaining = (target.get(key) || []).filter((item) => item !== identity);
  if (remaining.length > 0) target.set(key, remaining);
  else target.delete(key);
};

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

export const collectRenderedFeatureIdentitiesFromCatalogItem = (item) => {
  const collection = createRenderedIdentityCollection();
  const recordIndexByKey = new Map(
    (Array.isArray(item?.recordKeys) ? item.recordKeys : []).map(
      (recordKey, recordIndex) => [String(recordKey || '').trim(), recordIndex]
    )
  );
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
