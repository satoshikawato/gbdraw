// @ts-check
/** @import { FeatureOverrideDraft, FeatureOverrideDraftRow, FeaturePlacementDraft, FeatureRequestRecord } from './feature-placement.js' */
import {
  canonicalFeaturePlacements,
  featureIdentityKey,
  normalizeFeatureOverrideLabelText,
  recordKeyBelongsToRequest
} from './feature-placement.js';
import { diagnosticError } from '../utils/error-normalization.js';

// Session 44 and older kept per-feature edits in four maps keyed by rendered
// SVG ID. Session 46 keys them by original-source feature identity
// (`features.featureOverrides`, design Q4 4.3). This reader maps each old key
// through the Session's saved feature catalog once, then the current model
// applies. Each migrated row names the mode of the Session's diagram, whose
// records the catalog names (R2): the Session 46 split (mode-scoped-migration.js)
// moves it into that mode's drawing and drops the mode. The Python twins of
// these migrations share their vectors, so the rows keep that shape.

// A migrated row's key: [scope, recordKey, biologicalFeatureId].
/**
 * @param {unknown} scope
 * @param {string} recordKey
 * @param {string} biologicalFeatureId
 * @returns {string}
 */
const scopedIdentityKey = (scope, recordKey, biologicalFeatureId) => (
  (scope === 'circular' || scope === 'linear') && featureIdentityKey(recordKey, biologicalFeatureId)
    ? JSON.stringify([scope, recordKey, biologicalFeatureId])
    : ''
);
/** @param {Record<string, any>} row */
const scopedIdentityKeyOf = (row) => scopedIdentityKey(row?.scope, row?.recordKey, row?.biologicalFeatureId);
export const RENDERED_ID_FEATURE_EDIT_FIELDS = Object.freeze([
  'featureVisibilityOverrides',
  'labelVisibilityOverrides',
  'labelTextFeatureOverrides',
  'labelTextFeatureOverrideSources'
]);

export const FEATURE_EDIT_MIGRATION_WARNING = (count) => (
  `${count} feature edit(s) from an older Session could not be matched to a feature of its saved diagram and were dropped.`
);

export const FEATURE_VISIBILITY_NARROWED_NOTICE = (count) => (
  `${count} Feature visibility edit(s) from an older Session hid every feature with the same hash, such as each copy`
  + ' of a duplicated record. Each now applies only to the feature that was edited, so the next Generate draws the others.'
);

/**
 * The input of a Session feature-edit migration. `catalog` is the Session's
 * saved catalog (schema 3, 4, or 5), unvalidated.
 * @typedef {object} FeatureEditMigrationInput
 * @property {Record<string, any>} features
 * @property {string} mode
 * @property {Record<string, any> | null} [catalog]
 * @property {Record<string, any> | null} [legacy]
 */

/**
 * @param {unknown} value
 * @returns {value is Record<string, any>}
 */
const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);
const text = (value) => String(value ?? '').trim();
const RENDERED_PART_SUFFIX = /__(?:part|line)\d+$/;
const INSTANCE_SUFFIX = /__instance_([A-Za-z0-9_.-]+?)_[0-9a-f]{16}$/;

// Reads `<h>[_record_<n>][__instance_<s>_<digest>]...` (gbdraw/features/ids.py,
// gbdraw/svg/ids.py::instance_svg_id): the stable hash, the one-based record
// position in its Result, and the source feature index of a duplicate.
const parseRenderedId = (renderedId) => {
  let rest = text(renderedId).replace(RENDERED_PART_SUFFIX, '');
  /** @type {number | null} */
  let recordOrdinal = null;
  /** @type {number | null} */
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

const catalogIndex = (catalog, mode) => {
  const renderedById = new Map();
  const biological = [];
  (isObject(catalog) && Array.isArray(catalog.items) ? catalog.items : []).forEach((item) => {
    const recordKeys = Array.isArray(item?.recordKeys) ? item.recordKeys.map(text) : [];
    (Array.isArray(item?.features) ? item.features : []).forEach((feature) => {
      const key = scopedIdentityKey(mode, text(feature?.recordKey), text(feature?.biologicalFeatureId));
      const svgId = text(feature?.svgId).replace(RENDERED_PART_SUFFIX, '');
      if (!key || !svgId) return;
      if (!renderedById.has(svgId)) renderedById.set(svgId, new Set());
      renderedById.get(svgId).add(key);
    });
    (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).forEach((feature) => {
      const recordKey = text(feature?.recordKey);
      const biologicalFeatureId = text(feature?.biologicalFeatureId);
      const key = scopedIdentityKey(mode, recordKey, biologicalFeatureId);
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

// A Session before feature catalogs (31-33, 39) keys an edit by the rendered ID
// `<drawn hash>[_record_<n>][__instance_<s>_<digest>]`: the hash of the
// feature in the drawn (cropped, reverse-complemented) record and the record's
// position in the Result. `features` are the Session's source features read
// again with its crops and orientations (each with its drawn hash beside its
// source hash), or else its saved feature metadata, which has the drawn hash
// only and so serves records drawn untransformed, where both hashes are equal.
// A feature's input is its request record (Linear: `fileIdx`; Circular: one
// file) and `record_idx` its record within that input. Identities are the
// renderer's: the record key, `<recordKey>:<n>` for each record of an ALL
// input with several records (gbdraw/api/record_planning.py), and the source
// hash, `~<source index>` when the record has it twice
// (gbdraw/features/ids.py::disambiguate_feature_ids).
const LINEAR_RECORD_SUFFIX = /_record_([1-9]\d*)(?=__|$)/;

// A request record drawn cropped, reverse-complemented, or rotated: its drawn
// hashes are not its source hashes.
const drawnTransformed = (record) => Boolean(record?.region || record?.presentation?.reverseComplement
  || Number.isSafeInteger(record?.display?.startCoordinate));

const nonnegativeInteger = (value) => (
  value !== null && value !== '' && Number.isSafeInteger(Number(value)) && Number(value) >= 0 ? Number(value) : null
);

const legacyIndex = ({ features = [], biologicalFeatures = [], records = [], mode = 'circular' } = {}) => {
  const linear = mode === 'linear';
  const requestRecords = Array.isArray(records) ? records : [];
  const describe = (feature) => {
    const svgOrdinal = Number(text(feature?.svg_id ?? feature?.svgId).match(LINEAR_RECORD_SUFFIX)?.[1]) || null;
    const fileIdx = nonnegativeInteger(feature?.fileIdx);
    return {
      input: linear ? fileIdx ?? (svgOrdinal ? svgOrdinal - 1 : null) : 0,
      recordIndex: linear && fileIdx === null ? 0 : nonnegativeInteger(feature?.record_idx ?? feature?.recordIndex),
      sourceIndex: [feature?.source_feature_index, feature?.sourceFeatureIndex, feature?.feature_index]
        .map(nonnegativeInteger).find((value) => value !== null) ?? null,
      sourceHash: text(feature?.stable_feature_id ?? feature?.stableFeatureId ?? feature?.stable_svg_id),
      drawnHash: text(feature?.drawn_selector?.hash ?? feature?.drawnSelector?.hash)
    };
  };
  const listed = (Array.isArray(features) ? features : []).map(describe)
    .filter((entry) => entry.input !== null && entry.recordIndex !== null && entry.sourceHash);
  const sources = (Array.isArray(biologicalFeatures) && biologicalFeatures.length > 0
    ? biologicalFeatures.map(describe).filter((entry) => entry.input !== null && entry.recordIndex !== null && entry.sourceHash)
    : listed);
  const recordCounts = new Map();
  [...listed, ...sources].forEach(({ input, recordIndex }) => {
    // `listed` and `sources` were filtered above to entries whose `recordIndex` is not null.
    recordCounts.set(input, Math.max(recordCounts.get(input) || 1, /** @type {number} */ (recordIndex) + 1));
  });
  const recordOf = (input, recordIndex) => (
    linear ? requestRecords[input] : requestRecords[requestRecords.length > 1 ? recordIndex : 0]
  );
  const recordKeyOf = (input, recordIndex) => {
    const record = recordOf(input, recordIndex);
    const recordKey = text(record?.recordKey);
    if (!recordKey) return '';
    return record.cardinality === 'all' && (linear || requestRecords.length === 1) && recordCounts.get(input) > 1
      ? `${recordKey}:${recordIndex + 1}`
      : recordKey;
  };
  const offsets = [];
  for (let input = 0, offset = 0; input < requestRecords.length; input += 1) {
    offsets[input] = offset;
    offset += recordCounts.get(input) || 1;
  }
  const ordinalOf = (input, recordIndex) => (linear ? (offsets[input] ?? 0) + recordIndex + 1 : recordIndex + 1);
  const transformed = (input, recordIndex) => drawnTransformed(recordOf(input, recordIndex));
  const hashCounts = new Map();
  sources.forEach((entry) => {
    const id = JSON.stringify([recordKeyOf(entry.input, entry.recordIndex), entry.sourceHash]);
    hashCounts.set(id, (hashCounts.get(id) || 0) + 1);
  });
  const biological = [];
  const seen = new Set();
  listed.forEach((entry) => {
    const drawnHash = entry.drawnHash || (transformed(entry.input, entry.recordIndex) ? '' : entry.sourceHash);
    const recordKey = recordKeyOf(entry.input, entry.recordIndex);
    if (!drawnHash || !recordKey) return;
    const duplicated = (hashCounts.get(JSON.stringify([recordKey, entry.sourceHash])) || 0) > 1;
    if (duplicated && entry.sourceIndex === null) return;
    const key = scopedIdentityKey(mode, recordKey, duplicated ? `${entry.sourceHash}~${entry.sourceIndex}` : entry.sourceHash);
    const recordOrdinal = ordinalOf(entry.input, entry.recordIndex);
    const once = JSON.stringify([key, drawnHash, recordOrdinal]);
    if (!key || seen.has(once)) return;
    seen.add(once);
    biological.push({ key, stableId: drawnHash, sourceIndex: entry.sourceIndex, recordOrdinal });
  });
  return { renderedById: new Map(), biological };
};

// Each feature a saved catalog (schema 3-5) drew, by the hash it was drawn
// with: the hash in its rendered ID (`svgId`), which schema 5 also keeps as
// `drawnSelector.hash`, and the record position of that ID or of the row in
// its item's records. The twin of `_catalog_drawn_features` in
// gbdraw/session_io.py.
const catalogDrawnFeatures = (catalog, mode) => {
  const drawn = [];
  (isObject(catalog) && Array.isArray(catalog.items) ? catalog.items : []).forEach((item) => {
    const recordKeys = Array.isArray(item?.recordKeys) ? item.recordKeys.map(text) : [];
    (Array.isArray(item?.features) ? item.features : []).forEach((feature) => {
      const recordKey = text(feature?.recordKey);
      const key = scopedIdentityKey(mode, recordKey, text(feature?.biologicalFeatureId));
      const rendered = parseRenderedId(feature?.svgId);
      const stableId = text(feature?.drawnSelector?.hash) || rendered.stableId;
      if (!key || !stableId) return;
      drawn.push({
        key,
        stableId,
        sourceIndex: rendered.sourceIndex,
        recordOrdinal: rendered.recordOrdinal ?? recordKeys.indexOf(recordKey) + 1
      });
    });
  });
  return drawn;
};

// OV-401, release D-39: a Session 44 or older matched `hash=` against the hash
// of the drawn feature, so on a cropped or reverse-complemented record its
// value is not the source hash `hash=` names now. The value, written as
// `<hash>` or `^<hash>$` and maybe with the rendered ID's record and instance
// suffixes (release 0.13.0), names the features of `drawn` drawn with that
// hash, in `recordKey` when given. Returns the value naming their one source
// hash, in the same form without the suffixes, and their identity keys; null
// when it names no drawn feature or features of two source hashes. The twin of
// `_source_hash_selector_value` in gbdraw/session_io.py.
const sourceHashSelectorValue = (value, drawn, recordKey = '') => {
  const saved = text(value);
  const anchored = saved.match(/^\^([\s\S]*)\$$/);
  const { stableId, recordOrdinal, sourceIndex } = parseRenderedId(anchored ? anchored[1] : saved);
  const keys = [...new Set(drawn.filter((feature) => stableId && feature.stableId === stableId
    && (recordOrdinal === null || feature.recordOrdinal === recordOrdinal)
    && (sourceIndex === null || feature.sourceIndex === sourceIndex)
    && (!recordKey || JSON.parse(feature.key)[1] === recordKey)).map((feature) => feature.key))];
  const sources = new Set(keys.map((key) => JSON.parse(key)[2].replace(/~\d+$/, '')));
  if (sources.size !== 1) return null;
  const [source] = sources;
  return { value: anchored ? `^${source}$` : source, keys };
};

/**
 * The request records of the Linear cards of a Session before 31, which has no
 * request (release 0.13.0): each card's record key, crop, and orientation, as
 * the hash readers read them. The twin of `_legacy_linear_request_records` in
 * gbdraw/session_io.py.
 * @param {unknown} linearSeqs
 */
export const legacyLinearRequestRecords = (linearSeqs) => (Array.isArray(linearSeqs) ? linearSeqs : [])
  .map((seq) => {
    const cropped = nonnegativeInteger(seq?.region_start) !== null || nonnegativeInteger(seq?.region_end) !== null;
    return {
      recordKey: text(seq?.uid),
      cardinality: cropped || text(seq?.region_record_id) ? 'exactly_one' : 'all',
      region: cropped ? { reverseComplement: Boolean(seq?.region_reverse) } : null,
      presentation: { reverseComplement: !cropped && Boolean(seq?.region_reverse) }
    };
  });

/**
 * Whether color rules (`qual`) or Feature visibility rules (`qualifier`) hold a `hash` rule.
 * @param {...unknown} ruleLists
 */
export const hasHashSelectorRules = (...ruleLists) => ruleLists.some((rules) => Array.isArray(rules)
  && rules.some((rule) => text(rule?.qual ?? rule?.qualifier).toLowerCase() === 'hash'));

export const HASH_SELECTOR_UNMAPPED_NOTICE = (count) => (
  `${count} hash= rule(s) or annotation target(s) from an older Session could not be matched to a feature of its saved diagram,`
  + ' which crops or reverse-complements a record. hash= now names a feature by its hash in the source record,'
  + ' so each may now match another feature or none.'
);

/**
 * S6 of the one feature-hash builder: the `hash` color rules (`rules`, `qual`
 * and `val`) and Feature visibility rules (`qualifier` and `value`) of a
 * Session 44 or older get the source hash of the features drawn with their
 * hash, through the saved `catalog` or else `legacy` (as
 * `migrateRenderedIdFeatureEdits` reads them), so each matches the features it
 * matched. When a request record is drawn transformed, every other `hash` value
 * is counted as unmapped and kept. Lists without a change are returned as
 * saved. The twin of `migrate_session_hash_rules` in gbdraw/session_io.py.
 * @param {{
 *   rules?: unknown,
 *   featureVisibilityManualRules?: unknown,
 *   mode: string,
 *   catalog?: Record<string, any> | null,
 *   legacy?: Record<string, any> | null,
 *   records?: unknown
 * }} input
 * @returns {{ rules: any, featureVisibilityManualRules: any, unmappedCount: number }}
 */
export const migrateSessionHashRules = ({
  rules, featureVisibilityManualRules, mode, catalog = null, legacy = null, records = []
}) => {
  const transformed = (Array.isArray(records) ? records : []).some((record) => isObject(record) && drawnTransformed(record));
  /** @type {any[] | null} */
  let drawn = null;
  let unmappedCount = 0;
  const migrate = (entries, qualifierField, valueField) => {
    if (!Array.isArray(entries)) return entries;
    let changed = false;
    const migrated = entries.map((entry) => {
      if (!isObject(entry) || text(entry[qualifierField]).toLowerCase() !== 'hash' || !text(entry[valueField])) return entry;
      drawn ??= catalog ? catalogDrawnFeatures(catalog, mode) : legacyIndex({ ...legacy, mode }).biological;
      const resolved = sourceHashSelectorValue(entry[valueField], drawn);
      if (!resolved) {
        if (transformed) unmappedCount += 1;
        return entry;
      }
      if (resolved.value === entry[valueField]) return entry;
      changed = true;
      return { ...entry, [valueField]: resolved.value };
    });
    return changed ? migrated : entries;
  };
  return {
    rules: migrate(rules, 'qual', 'val'),
    featureVisibilityManualRules: migrate(featureVisibilityManualRules, 'qualifier', 'value'),
    unmappedCount
  };
};

// Rule 1: the old key is a rendered ID of the saved catalog, so it names every
// identity drawn with that ID (the live projection reached all of them).
// Rule 2: otherwise it names the one feature with its hash in the record
// position (and source index) its suffixes give: the source hash of a catalog
// feature, the drawn hash of a feature of a Session without a catalog.
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

// The identities a Session before 46 reached with the `hash` row it sent for
// a Feature visibility edit: every feature drawn with, or whose source has,
// that hash (each copy of a duplicated record, each feature at the same
// coordinates). Owner decision Q1 = A keeps the edit on the one it named.
const identitiesWithHash = (index, hash) => {
  if (!index.byHash) {
    index.byHash = new Map();
    const add = (stableId, key) => {
      if (!index.byHash.has(stableId)) index.byHash.set(stableId, new Set());
      index.byHash.get(stableId).add(key);
    };
    index.renderedById.forEach((keys, svgId) => keys.forEach((key) => add(parseRenderedId(svgId).stableId, key)));
    index.biological.forEach(({ stableId, key }) => add(stableId, key));
  }
  return index.byHash.get(hash) || new Set();
};

// Maps, so a saved value such as `constructor` names no mode (OV-134).
const FEATURE_VISIBILITY_MODES = new Map([
  ['on', 'on'], ['off', 'off'], ['exclude_matching', 'exclude_matching'], ['suppress', 'exclude_matching']
]);

/**
 * @param {unknown} features
 * @returns {boolean}
 */
export const hasRenderedIdFeatureEdits = (features) => isObject(features)
  && RENDERED_ID_FEATURE_EDIT_FIELDS.some((field) => isObject(features[field]) && Object.keys(features[field]).length > 0);

/**
 * Maps the rendered-ID edit maps of a Session older than 46 to draft rows
 * keyed by source identity in the mode of its diagram (`mode`). `catalog` is
 * the Session's saved feature catalog (schema 3, 4, or 5); a Session without
 * one gives `legacy`: its request records and source features read again
 * (`features` with their drawn hashes, `biologicalFeatures`) or its saved
 * feature metadata. Returns the rows, the
 * number of dropped edits (Feature visibility, Label visibility, and label
 * text entries; a label's source text is not an edit of its own), the number
 * of Feature visibility edits that now reach fewer features than the `hash`
 * row the Session sent, and whether the Session had label edits.
 * @param {FeatureEditMigrationInput} input
 * @returns {{
 *   featureOverrides: FeatureOverrideDraft,
 *   droppedCount: number,
 *   narrowedVisibilityCount: number,
 *   migratedLabelEdits: boolean
 * }}
 */
export const migrateRenderedIdFeatureEdits = ({ features, mode, catalog = null, legacy = null }) => {
  /** @type {FeatureOverrideDraft} */
  const rows = {};
  let droppedCount = 0;
  let narrowedVisibilityCount = 0;
  const migratedLabelEdits = ['labelVisibilityOverrides', 'labelTextFeatureOverrides']
    .some((field) => isObject(features?.[field]) && Object.keys(features[field]).length > 0);
  if (!hasRenderedIdFeatureEdits(features)) {
    return { featureOverrides: rows, droppedCount, narrowedVisibilityCount, migratedLabelEdits };
  }
  const index = catalog ? catalogIndex(catalog, mode) : legacyIndex({ ...legacy, mode });
  const rowFor = (key) => {
    if (!rows[key]) {
      // `key` comes from `scopedIdentityKey` through the index.
      const [scope, recordKey, biologicalFeatureId] = JSON.parse(key);
      rows[key] = /** @type {FeatureOverrideDraftRow} */ ({
        scope,
        recordKey,
        biologicalFeatureId,
        featureVisibility: null,
        labelVisibility: null,
        labelText: null,
        labelSourceText: null
      });
    }
    return rows[key];
  };
  const migrateMap = (field, assign) => {
    Object.entries(isObject(features[field]) ? features[field] : {}).forEach(([oldKey, value]) => {
      const keys = resolveOldKey(oldKey, index);
      if (keys.length === 0 || !keys.every((key) => assign(rowFor(key), value) !== false)) {
        if (field !== 'labelTextFeatureOverrideSources') droppedCount += 1;
      } else if (field === 'featureVisibilityOverrides'
        && [...identitiesWithHash(index, parseRenderedId(oldKey).stableId)].some((key) => !keys.includes(key))) {
        narrowedVisibilityCount += 1;
      }
    });
  };
  migrateMap('featureVisibilityOverrides', (row, value) => {
    const mode = FEATURE_VISIBILITY_MODES.get(text(value).toLowerCase());
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
  return { featureOverrides: rows, droppedCount, narrowedVisibilityCount, migratedLabelEdits };
};

/**
 * The Session 46 `features` of an older Session's `features`: the edit maps
 * become `featureOverrides`. When a label map is migrated, the saved label
 * table is cleared: it was the copy of those maps built at the last Generate
 * (the same invariant the editor keeps when label edits change).
 * @param {FeatureEditMigrationInput} input
 * @returns {{
 *   features: Record<string, any>,
 *   droppedCount: number,
 *   narrowedVisibilityCount: number
 * }}
 */
export const migrateSessionFeatureEdits = ({ features, mode, catalog = null, legacy = null }) => {
  const source = isObject(features) ? features : {};
  const migration = migrateRenderedIdFeatureEdits({ features: source, mode, catalog, legacy });
  const migrated = { ...source };
  RENDERED_ID_FEATURE_EDIT_FIELDS.forEach((field) => delete migrated[field]);
  migrated.featureOverrides = migration.featureOverrides;
  if (migration.migratedLabelEdits) migrated.labelOverrideRows = [];
  return {
    features: migrated,
    droppedCount: migration.droppedCount,
    narrowedVisibilityCount: migration.narrowedVisibilityCount
  };
};

export const ANNOTATION_TARGET_MIGRATION_NOTICE = (count) => (
  `${count} annotation(s) from an older Session named a feature by hash=. Each now names that feature by its source,`
  + ' so it stays on the feature when the crop or orientation changes.'
);

// The records the renderer binds an annotation target in (the request's drawn
// records in order, gbdraw/annotations/resolve.py::_bind_record) and the
// source features by source hash, from the Session's saved feature catalog.
const annotationCatalogIndex = (catalog) => {
  const recordIds = new Map();
  const featuresByHash = new Map();
  (isObject(catalog) && Array.isArray(catalog.items) ? catalog.items : []).forEach((item) => {
    (Array.isArray(item?.recordKeys) ? item.recordKeys : []).forEach((recordKey) => {
      if (!recordIds.has(text(recordKey))) recordIds.set(text(recordKey), new Set());
    });
    (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).forEach((feature) => {
      const recordKey = text(feature?.recordKey);
      const biologicalFeatureId = text(feature?.biologicalFeatureId);
      if (!recordKey || !biologicalFeatureId) return;
      recordIds.get(recordKey)?.add(text(feature?.record_id ?? feature?.recordId));
      const hash = text(feature?.stableFeatureId) || biologicalFeatureId.replace(/~\d+$/, '');
      if (!featuresByHash.has(hash)) featuresByHash.set(hash, new Map());
      featuresByHash.get(hash).set(JSON.stringify([recordKey, biologicalFeatureId]), { recordKey, biologicalFeatureId });
    });
  });
  return { recordKeys: [...recordIds.keys()], recordIds, featuresByHash };
};

const boundRecordKey = (selector, { recordKeys, recordIds }) => {
  if (selector == null) return recordKeys.length === 1 ? recordKeys[0] : '';
  if (selector.kind === 'recordIndex') return recordKeys[selector.index] || '';
  if (selector.kind !== 'recordId') return '';
  // A record without catalog features has no known ID, so the binding is not certain.
  if (recordKeys.some((recordKey) => recordIds.get(recordKey).size !== 1)) return '';
  const matches = recordKeys.filter((recordKey) => recordIds.get(recordKey).has(text(selector.value)));
  return matches.length === 1 ? matches[0] : '';
};

/**
 * R-7 (Owner decision 2026-10-05): a Session before 46 named a selected
 * feature in an annotation by `hash=<hash>` (a featureSpan target), which the
 * renderer matched in the drawn record, and a hand-written target looks the
 * same. Load moves such a target to the feature's source identity, in the mode
 * of the Session's diagram (`mode`), only when the figure cannot change. On a
 * record drawn without a crop, reverse complement, or rotation, the hash must
 * name exactly one feature of the saved catalog, in the record the target
 * binds. On a record drawn transformed, it must name the features the catalog
 * drew with that hash there (S7, `sourceHashSelectorValue`): one feature gives
 * its identity, features of one source hash keep the featureSpan with that
 * hash. Every other target stays as saved, and one that may sit on a
 * transformed record is counted as unmapped. Returns the annotation sets, the
 * number of moved targets, and the number of unmapped targets.
 * @param {{
 *   annotationSets: unknown,
 *   mode: string,
 *   catalog?: Record<string, any> | null,
 *   records?: Array<FeatureRequestRecord & Record<string, any>>
 * }} input
 * @returns {{ annotationSets: any, migratedCount: number, unmappedCount: number }}
 */
export const migrateSessionAnnotationTargets = ({ annotationSets, mode, catalog = null, records = [] }) => {
  const index = annotationCatalogIndex(catalog);
  const requestRecords = Array.isArray(records) ? records : [];
  const drawn = catalogDrawnFeatures(catalog, mode);
  const anyTransformed = requestRecords.some((record) => isObject(record) && drawnTransformed(record));
  let migratedCount = 0;
  let unmappedCount = 0;
  const identityTarget = (target) => {
    const selector = target?.kind === 'featureSpan' && Array.isArray(target.selectors) && target.selectors.length === 1
      ? target.selectors[0] : null;
    if (selector?.key !== 'hash') return null;
    const recordKey = boundRecordKey(target.record, index);
    const request = recordKey ? requestRecords.find((record) => recordKeyBelongsToRequest(recordKey, [record])) : null;
    let matches;
    if (!request || drawnTransformed(request)) {
      // S7: on a transformed record the hash names the features the saved
      // catalog drew with it there.
      const resolved = request ? sourceHashSelectorValue(selector.value, drawn, recordKey) : null;
      if (!resolved) {
        if (anyTransformed) unmappedCount += 1;
        return null;
      }
      if (resolved.keys.length !== 1) return { ...target, selectors: [{ ...selector, value: resolved.value }] };
      const [, boundKey, biologicalFeatureId] = JSON.parse(resolved.keys[0]);
      matches = [{ recordKey: boundKey, biologicalFeatureId }];
    } else {
      matches = [...(index.featuresByHash.get(text(selector.value))?.values() || [])];
      if (matches.length !== 1 || matches[0].recordKey !== recordKey) return null;
    }
    const migrated = { kind: 'featureIdentity', scope: mode, ...matches[0],
      envelope: target.envelope, circularPath: target.circularPath };
    return scopedIdentityKeyOf(migrated) ? migrated : null;
  };
  const migratedSets = (Array.isArray(annotationSets) ? annotationSets : []).map((set) => (
    !Array.isArray(set?.annotations) ? set : {
      ...set,
      annotations: set.annotations.map((annotation) => {
        const target = identityTarget(annotation?.target);
        if (!target) return annotation;
        migratedCount += 1;
        return { ...annotation, target };
      })
    }
  ));
  return { annotationSets: migratedCount > 0 ? migratedSets : annotationSets, migratedCount, unmappedCount };
};

const LANE_MODES = new Map([['outward', 'circular'], ['inward', 'circular'], ['above', 'linear'], ['below', 'linear']]);

/**
 * The mode-scoped Feature placement drafts of a Session 41-44 draft keyed by
 * [recordKey, biologicalFeatureId]. Such a row reached every request with its
 * record key, so a Main row is kept for both modes and a lane row for the mode
 * of its side. A row whose key does not encode it is kept as is, and the draft
 * check rejects it.
 * @param {any} placements
 * @returns {any} The Feature placement draft, or the unvalidated input when it is not a map.
 */
export const migrateSessionFeaturePlacements = (placements) => {
  if (!isObject(placements)) return placements;
  // Entries, so a saved key such as `__proto__` stays a row of the draft.
  /** @type {Array<[string, any]>} */
  const migrated = [];
  Object.entries(placements).forEach(([key, row]) => {
    const target = row?.placement;
    const modes = key === JSON.stringify([row?.recordKey, row?.biologicalFeatureId])
      ? target?.kind === 'main' ? ['circular', 'linear'] : [LANE_MODES.get(target?.side)].filter(Boolean) : [];
    if (modes.length === 0) migrated.push([key, row]);
    modes.forEach((scope) => {
      const scoped = { scope, ...row };
      migrated.push([scopedIdentityKeyOf(scoped) || key, scoped]);
    });
  });
  return Object.fromEntries(migrated);
};

// A Session 41-44 draft's placement rows as `migrateSessionFeaturePlacements`
// names their mode (`[scope, recordKey, biologicalFeatureId]`); each is checked
// as a row of its mode. The Session 46 split moves each row into its mode's
// drawing.
/** @param {unknown} placements */
export const validateScopedFeaturePlacements = (placements) => {
  if (!isObject(placements)) throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
  for (const [key, row] of Object.entries(/** @type {Record<string, any>} */ (placements))) {
    const { scope, ...rest } = isObject(row) ? row : {};
    if (!['circular', 'linear'].includes(scope) || key !== scopedIdentityKeyOf(row)) {
      throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
    }
    canonicalFeaturePlacements(/** @type {FeaturePlacementDraft} */ ({
      [featureIdentityKey(rest.recordKey, rest.biologicalFeatureId)]: rest
    }), scope);
  }
};
