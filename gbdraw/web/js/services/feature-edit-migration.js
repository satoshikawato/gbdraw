import { featureIdentityKey, normalizeFeatureOverrideLabelText } from './feature-placement.js';

// Session 44 and older kept per-feature edits in four maps keyed by rendered
// SVG ID. Session 45 keys them by original-source feature identity
// (`features.featureOverrides`, design Q4 4.3). This reader maps each old key
// through the Session's saved feature catalog once, then the current model
// applies.
export const RENDERED_ID_FEATURE_EDIT_FIELDS = Object.freeze([
  'featureVisibilityOverrides',
  'labelVisibilityOverrides',
  'labelTextFeatureOverrides',
  'labelTextFeatureOverrideSources'
]);

export const FEATURE_EDIT_MIGRATION_WARNING = (count) => (
  `${count} feature edit(s) from an older Session could not be matched to a feature of its saved diagram and were dropped.`
);

const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);
const text = (value) => String(value ?? '').trim();
const RENDERED_PART_SUFFIX = /__(?:part|line)\d+$/;
const INSTANCE_SUFFIX = /__instance_([A-Za-z0-9_.-]+?)_[0-9a-f]{16}$/;

// Reads `<h>[_record_<n>][__instance_<s>_<digest>]...` (gbdraw/features/ids.py,
// gbdraw/svg/ids.py::instance_svg_id): the stable hash, the one-based record
// position in its Result, and the source feature index of a duplicate.
const parseRenderedId = (renderedId) => {
  let rest = text(renderedId).replace(RENDERED_PART_SUFFIX, '');
  let recordOrdinal = null;
  let sourceIndex = null;
  for (let match = rest.match(INSTANCE_SUFFIX); match; match = rest.match(INSTANCE_SUFFIX)) {
    const instance = match[1];
    const record = instance.match(/^record_([1-9]\d*)$/);
    if (record) recordOrdinal ??= Number(record[1]);
    else if (/^(0|[1-9]\d*)$/.test(instance)) sourceIndex ??= Number(instance);
    rest = rest.slice(0, match.index);
  }
  const linearRecord = rest.match(/^(.*)_record_([1-9]\d*)$/);
  if (linearRecord) {
    rest = linearRecord[1];
    recordOrdinal ??= Number(linearRecord[2]);
  }
  return { stableId: rest, recordOrdinal, sourceIndex };
};

const catalogIndex = (catalog) => {
  const renderedById = new Map();
  const biological = [];
  (isObject(catalog) && Array.isArray(catalog.items) ? catalog.items : []).forEach((item) => {
    const recordKeys = Array.isArray(item?.recordKeys) ? item.recordKeys.map(text) : [];
    (Array.isArray(item?.features) ? item.features : []).forEach((feature) => {
      const key = featureIdentityKey(text(feature?.recordKey), text(feature?.biologicalFeatureId));
      const svgId = text(feature?.svgId).replace(RENDERED_PART_SUFFIX, '');
      if (!key || !svgId) return;
      if (!renderedById.has(svgId)) renderedById.set(svgId, new Set());
      renderedById.get(svgId).add(key);
    });
    (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).forEach((feature) => {
      const recordKey = text(feature?.recordKey);
      const biologicalFeatureId = text(feature?.biologicalFeatureId);
      const key = featureIdentityKey(recordKey, biologicalFeatureId);
      if (!key) return;
      biological.push({
        key,
        stableId: text(feature?.stableFeatureId) || biologicalFeatureId,
        sourceIndex: Number.isSafeInteger(feature?.sourceFeatureIndex) ? feature.sourceFeatureIndex : null,
        recordOrdinal: recordKeys.indexOf(recordKey) + 1
      });
    });
  });
  return { renderedById, biological };
};

// A Session before feature catalogs (31-33, 39 without one) keeps only feature
// metadata recovered from its sources or saved with it. Its identities are the
// renderer's: the record key of the request record the feature belongs to, and
// the stable hash, with `~<source index>` when the record has the hash twice
// (gbdraw/features/ids.py::disambiguate_feature_ids). Without one record key
// per exactly-one record input, no key resolves.
const legacyIndex = (features, records) => {
  const renderedById = new Map();
  const biological = [];
  const recordKeys = (Array.isArray(records) ? records : []).map((record) => (
    record?.cardinality === 'all' ? '' : text(record?.recordKey)
  ));
  if (recordKeys.length === 0 || recordKeys.some((key) => !key)) return { renderedById, biological };
  const entries = (Array.isArray(features) ? features : []).map((feature) => {
    const recordIndex = Number(feature?.record_idx ?? feature?.recordIndex);
    const sourceIndex = [feature?.source_feature_index, feature?.sourceFeatureIndex, feature?.feature_index]
      .find((value) => Number.isSafeInteger(value) && value >= 0);
    return {
      recordKey: Number.isInteger(recordIndex) ? recordKeys[recordIndex] || '' : '',
      recordOrdinal: Number.isInteger(recordIndex) ? recordIndex + 1 : 0,
      stableId: text(feature?.stable_feature_id ?? feature?.stableFeatureId ?? feature?.stable_svg_id),
      sourceIndex: sourceIndex ?? null,
      svgId: text(feature?.rendered_feature_svg_id ?? feature?.svg_id).replace(RENDERED_PART_SUFFIX, '')
    };
  }).filter((entry) => entry.recordKey && entry.stableId);
  const counts = new Map();
  entries.forEach((entry) => {
    const id = JSON.stringify([entry.recordKey, entry.stableId]);
    counts.set(id, (counts.get(id) || 0) + 1);
  });
  entries.forEach((entry) => {
    const duplicated = counts.get(JSON.stringify([entry.recordKey, entry.stableId])) > 1;
    if (duplicated && entry.sourceIndex === null) return;
    const key = featureIdentityKey(entry.recordKey, duplicated ? `${entry.stableId}~${entry.sourceIndex}` : entry.stableId);
    if (!key) return;
    if (entry.svgId) {
      if (!renderedById.has(entry.svgId)) renderedById.set(entry.svgId, new Set());
      renderedById.get(entry.svgId).add(key);
    }
    biological.push({ key, stableId: entry.stableId, sourceIndex: entry.sourceIndex, recordOrdinal: entry.recordOrdinal });
  });
  return { renderedById, biological };
};

// Rule 1: the old key is a rendered ID of the saved catalog, so it names every
// identity drawn with that ID (the live projection reached all of them).
// Rule 2: otherwise it names the one biological feature with its stable hash in
// the record position (and source index) its suffixes give.
const resolveOldKey = (renderedId, index) => {
  const drawn = index.renderedById.get(text(renderedId).replace(RENDERED_PART_SUFFIX, ''));
  if (drawn?.size) return [...drawn];
  const { stableId, recordOrdinal, sourceIndex } = parseRenderedId(renderedId);
  if (!stableId) return [];
  const candidates = index.biological.filter((feature) => feature.stableId === stableId
    && (recordOrdinal === null || feature.recordOrdinal === recordOrdinal)
    && (sourceIndex === null || feature.sourceIndex === sourceIndex));
  return candidates.length === 1 ? [candidates[0].key] : [];
};

const FEATURE_VISIBILITY_MODES = Object.freeze({
  on: 'on', off: 'off', exclude_matching: 'exclude_matching', suppress: 'exclude_matching'
});

export const hasRenderedIdFeatureEdits = (features) => isObject(features)
  && RENDERED_ID_FEATURE_EDIT_FIELDS.some((field) => isObject(features[field]) && Object.keys(features[field]).length > 0);

/**
 * Maps the rendered-ID edit maps of a Session older than 45 to draft rows
 * keyed by source identity. `catalog` is the Session's saved feature catalog
 * (schema 3, 4, or 5); a Session without one gives its saved or recovered
 * feature metadata and request records (`legacy`). Returns the rows, the
 * number of dropped edits, and whether the Session had label edits.
 */
export const migrateRenderedIdFeatureEdits = ({ features, catalog = null, legacy = null }) => {
  const rows = {};
  let droppedCount = 0;
  const migratedLabelEdits = ['labelVisibilityOverrides', 'labelTextFeatureOverrides']
    .some((field) => isObject(features?.[field]) && Object.keys(features[field]).length > 0);
  if (!hasRenderedIdFeatureEdits(features)) return { featureOverrides: rows, droppedCount, migratedLabelEdits };
  const index = catalog ? catalogIndex(catalog) : legacyIndex(legacy?.features, legacy?.records);
  const rowFor = (key) => {
    if (!rows[key]) {
      const [recordKey, biologicalFeatureId] = JSON.parse(key);
      rows[key] = {
        recordKey,
        biologicalFeatureId,
        featureVisibility: null,
        labelVisibility: null,
        labelText: null,
        labelSourceText: null
      };
    }
    return rows[key];
  };
  const migrateMap = (field, assign) => {
    Object.entries(isObject(features[field]) ? features[field] : {}).forEach(([oldKey, value]) => {
      const keys = resolveOldKey(oldKey, index);
      if (keys.length === 0 || !keys.every((key) => assign(rowFor(key), value) !== false)) {
        droppedCount += 1;
      }
    });
  };
  migrateMap('featureVisibilityOverrides', (row, value) => {
    const mode = FEATURE_VISIBILITY_MODES[text(value).toLowerCase()];
    if (!mode) return false;
    row.featureVisibility ??= mode;
    return true;
  });
  migrateMap('labelVisibilityOverrides', (row, value) => {
    const mode = text(value).toLowerCase();
    if (mode !== 'on' && mode !== 'off') return false;
    row.labelVisibility ??= mode;
    return true;
  });
  migrateMap('labelTextFeatureOverrides', (row, value) => {
    const labelText = normalizeFeatureOverrideLabelText(String(value ?? ''));
    // A blank text hid the label (its table row drew an empty label).
    if (labelText) row.labelText ??= labelText;
    else if (row.labelVisibility !== 'on') row.labelVisibility = 'off';
    return true;
  });
  migrateMap('labelTextFeatureOverrideSources', (row, value) => {
    const sourceText = String(value ?? '');
    if (!sourceText || sourceText.includes('\0')) return false;
    row.labelSourceText ??= sourceText;
    return true;
  });
  // A source text alone is kept only for a feature a bulk label edit can reach.
  Object.entries(rows).forEach(([key, row]) => {
    if (['featureVisibility', 'labelVisibility', 'labelText', 'labelSourceText'].every((field) => row[field] === null)) {
      delete rows[key];
    }
  });
  return { featureOverrides: rows, droppedCount, migratedLabelEdits };
};

/**
 * The Session 45 `features` of an older Session's `features`: the edit maps
 * become `featureOverrides`. When a label map is migrated, the saved label
 * table is cleared: it was the copy of those maps built at the last Generate
 * (the same invariant the editor keeps when label edits change).
 */
export const migrateSessionFeatureEdits = ({ features, catalog = null, legacy = null }) => {
  const source = isObject(features) ? features : {};
  const migration = migrateRenderedIdFeatureEdits({ features: source, catalog, legacy });
  const migrated = { ...source };
  RENDERED_ID_FEATURE_EDIT_FIELDS.forEach((field) => delete migrated[field]);
  migrated.featureOverrides = migration.featureOverrides;
  if (migration.migratedLabelEdits) migrated.labelOverrideRows = [];
  return { features: migrated, droppedCount: migration.droppedCount };
};
