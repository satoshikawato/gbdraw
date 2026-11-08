// @ts-check
import {
  annotationRecordSelectorValue,
  parseAnnotationRecordSelectorValue
} from './target-actions.js';
import { ANNOTATION_RECORD_BINDING_KEY } from '../../services/annotation-state.js';
import { bindAnnotationRecords } from '../../services/drawing/annotation-sets.js';
import { bindingKeyOfRequestRecord } from '../../services/record-draw-selection.js';
import { withRecordSelectors } from './record-catalog.js';
import { t } from '../../utils/i18n.js';

/**
 * @import { AnnotationCatalogRecord, AnnotationRecordCatalog } from './record-catalog.js'
 * @import { AnnotationRecordSelector } from '../../services/annotation-state.js'
 * @import { AnnotationRecordBinding, AnnotationRow, AnnotationSetsSlice, AnnotationSetsView } from '../../services/drawing/annotation-sets.js'
 * @import { RecordSettingsOwner } from '../../services/record-draw-selection.js'
 * @import { RegionAxisRecord } from '../../services/region-axis.js'
 * @import { ResultRecordBindings } from '../../services/sequence-search.js'
 */
/**
 * The annotation sets of a drawing, and the handle their writers take.
 * @typedef {AnnotationSetsView & { readonly slices: { readonly annotationSets: AnnotationSetsSlice } }} AnnotationSetsDrawing
 */

/** @param {unknown} value */
const cleanText = (value) => String(value ?? '').trim();

/** @param {AnnotationRow | null | undefined} annotation */
const bindingFor = (annotation) => cleanText(annotation?.metadata?.[ANNOTATION_RECORD_BINDING_KEY]);
export const annotationRecordBinding = bindingFor;

/**
 * @param {AnnotationRow | null | undefined} annotation
 * @returns {{ selector: AnnotationRecordSelector | null, error: string }}
 */
const selectorFromTarget = (annotation) => {
  const raw = annotation?.target?.record;
  if (raw == null) return { selector: null, error: '' };
  if (raw?.kind === 'recordIndex') {
    const index = Number(raw.index);
    return Number.isInteger(index) && index >= 0
      ? { selector: { kind: 'recordIndex', index }, error: '' }
      : { selector: null, error: 'Record index must be a non-negative integer.' };
  }
  if (raw?.kind === 'recordId') {
    const parsed = parseAnnotationRecordSelectorValue(raw.value);
    if (parsed.error || parsed.selector?.kind !== 'recordId') {
      return { selector: null, error: parsed.error || 'Record ID is reserved as an unspecified selector.' };
    }
    return parsed;
  }
  return { selector: null, error: 'Unsupported record selector.' };
};
export const annotationRecordSelectorFromTarget = selectorFromTarget;

/**
 * @param {AnnotationCatalogRecord[]} records
 * @param {AnnotationRecordSelector | null} selector
 */
const recordForSelector = (records, selector) => {
  if (!selector) return null;
  if (selector.kind === 'recordIndex') return records[selector.index] || null;
  const matches = records.filter((record) => record.recordId === selector.value);
  return matches.length === 1 ? matches[0] : null;
};

/**
 * @param {AnnotationRecordCatalog | null | undefined} catalog
 * @param {AnnotationRow | null | undefined} annotation
 */
const resolvedRecord = (catalog, annotation) => {
  const records = Array.isArray(catalog?.records) ? catalog.records : [];
  const binding = bindingFor(annotation);
  if (binding) return records.find((record) => record.key === binding) || null;
  const parsed = selectorFromTarget(annotation);
  return parsed.error ? null : recordForSelector(records, parsed.selector);
};
export const resolveAnnotationRecord = resolvedRecord;

/**
 * The annotation's metadata with the record binding `key` (none when empty).
 * @param {AnnotationRow} annotation
 * @param {string} key
 * @returns {Record<string, unknown>}
 */
const metadataWithBinding = (annotation, key) => {
  const { metadata } = annotation;
  /** @type {Record<string, unknown>} */
  const next = metadata && typeof metadata === 'object' && !Array.isArray(metadata) ? { ...metadata } : {};
  if (key) next[ANNOTATION_RECORD_BINDING_KEY] = key;
  else delete next[ANNOTATION_RECORD_BINDING_KEY];
  return next;
};

/**
 * @param {AnnotationRecordCatalog | null | undefined} catalog
 * @param {AnnotationRow | null | undefined} annotation
 */
export const annotationRecordValue = (catalog, annotation) => {
  const binding = bindingFor(annotation);
  if (binding) return binding;
  const record = resolvedRecord(catalog, annotation);
  if (record) return record.key;
  const saved = annotationRecordSelectorValue(annotation?.target?.record);
  return saved ? `saved:${saved}` : '';
};

/**
 * How a row binds a record of a Result (its request key), as Add to annotation set and Sequence search bind one: by
 * the catalog key the record select writes, that of the record settings owner that holds it (`bindingKeyOf`; an OFF
 * record's too), and by the selector that names it among the records the Result draws, those of the committed request
 * it was made from (`withRecordSelectors`, as the catalog names the draft's; review r2 B1). The draft's catalog is not
 * read: a record change that waits for Generate changes neither. A record no owner holds (its card was deleted after
 * Generate) binds by a key no catalog record has (a catalog key ends in a JSON array; this one is a JSON object), so
 * its row's target is unavailable as any row's whose record is gone, and Generate refuses it until the row is deleted
 * or rebound, never binding it by selector to another record (OV-2631, D-57); the live draw still names the record.
 * A `#n` selector counts records, not drawn units: the windows of a record drawn in several are one record (I-6), so
 * a record's number is its unit's index less the later windows drawn before it.
 * @param {readonly Pick<RegionAxisRecord, 'recordKey' | 'recordId' | 'recordIndex' | 'windowIndex'>[]} axisRecords The
 *   Result's records.
 * @param {readonly RecordSettingsOwner[]} owners The drawing's record settings owners.
 * @returns {ResultRecordBindings}
 */
export const resultRecordBindings = (axisRecords, owners) => {
  // Each record once, by its first drawn unit.
  const drawn = [...axisRecords].sort((a, b) => a.recordIndex - b.recordIndex)
    .filter((record, index, sorted) => sorted.findIndex((other) => other.recordKey === record.recordKey) === index);
  /** @param {number} recordIndex */
  const recordNumber = (recordIndex) => recordIndex + 1
    - axisRecords.filter((unit) => Boolean(unit.windowIndex) && unit.recordIndex < recordIndex).length;
  const named = withRecordSelectors(drawn.map(({ recordKey, recordId, recordIndex }) => ({
    recordKey, recordId, selector: `#${recordNumber(recordIndex)}`
  })));
  /** @param {string} recordKey */
  const bindingKeyOf = (recordKey) => bindingKeyOfRequestRecord(owners, recordKey) || JSON.stringify({ resultRecord: recordKey });
  return {
    bindingKeyOf,
    selectorOf: (recordKey) => named.find((record) => record.recordKey === recordKey)?.backendSelector ?? null,
    selectorOfBinding: (bindingKey) => (bindingKey
      ? named.find(({ recordKey }) => bindingKeyOf(recordKey) === bindingKey)?.backendSelector ?? null : null)
  };
};

/**
 * A new annotation bound to the catalog `record` as the record select binds one: its target names the record's
 * selector and its metadata the record's key.
 * @template {AnnotationRow} T
 * @param {T} annotation
 * @param {Pick<AnnotationCatalogRecord, 'key' | 'backendSelector'>} record
 * @returns {T}
 */
export const boundAnnotation = (annotation, record) => ({
  ...annotation,
  target: { ...annotation.target, record: { ...record.backendSelector } },
  metadata: metadataWithBinding(annotation, record.key)
});

/**
 * Binds an annotation of the drawing to the record `value` names (none when
 * empty); a key of no record changes nothing.
 * @param {AnnotationRecordCatalog | null | undefined} catalog
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationRow | null | undefined} annotation
 * @param {unknown} value
 */
export const setAnnotationRecordValue = (catalog, slice, annotation, value) => {
  if (!annotation?.target) return;
  const key = cleanText(value);
  if (!key) {
    bindAnnotationRecords(slice, [{ annotation, record: null, metadata: metadataWithBinding(annotation, '') }]);
    return;
  }
  const record = (Array.isArray(catalog?.records) ? catalog.records : [])
    .find((candidate) => candidate.key === key);
  if (!record) return;
  bindAnnotationRecords(slice, [{ annotation, record: record.backendSelector, metadata: metadataWithBinding(annotation, record.key) }]);
};

/**
 * Binds each annotation of the drawing whose record the catalog resolves to
 * that record (a selected-feature target names its record by identity).
 * @param {AnnotationSetsDrawing} drawing
 * @param {AnnotationRecordCatalog | null | undefined} catalog
 */
export const reconcileAnnotationRecordBindings = (drawing, catalog) => {
  const records = Array.isArray(catalog?.records) ? catalog.records : [];
  /** @type {AnnotationRecordBinding[]} */
  const bindings = [];
  (Array.isArray(drawing.annotationSets) ? drawing.annotationSets : []).forEach((set) => {
    (Array.isArray(set?.annotations) ? set.annotations : []).forEach((/** @type {AnnotationRow} */ annotation) => {
      if (!annotation?.target || annotation.target.kind === 'featureIdentity') return;
      const binding = bindingFor(annotation);
      const record = binding
        ? records.find((candidate) => candidate.key === binding)
        : resolvedRecord(catalog, annotation);
      if (record) bindings.push({ annotation, record: record.backendSelector, metadata: metadataWithBinding(annotation, record.key) });
    });
  });
  bindAnnotationRecords(drawing.slices.annotationSets, bindings);
};

// The select shows only for a ready catalog (`showsSelect`), so its options
// need no wording for a catalog that is not ready.
/**
 * @param {AnnotationRecordCatalog | null | undefined} catalog
 * @param {AnnotationRow | null | undefined} annotation
 */
export const annotationRecordOptions = (catalog, annotation) => {
  const records = Array.isArray(catalog?.records) ? catalog.records : [];
  const currentValue = annotationRecordValue(catalog, annotation);
  const options = records.map((record) => ({
    value: record.key,
    label: record.label,
    synthetic: false
  }));
  const selectedIsMissing = Boolean(currentValue) && !options.some((option) => option.value === currentValue);
  const savedSelector = annotationRecordSelectorValue(annotation?.target?.record) || t('saved target');
  return [
    { value: '', label: records.length > 1 ? t('Select target record') : t('Automatic (single record)'), synthetic: false },
    ...(selectedIsMissing
      ? [{ value: currentValue, label: t('{selector} (target record unavailable)', { selector: savedSelector }), synthetic: true }]
      : []),
    ...options
  ];
};

// OV-373: a loaded Session reads no record bytes, so a catalog can wait for
// discovery; a record an annotation names is then not inspected, not unavailable.
/** @param {AnnotationRecordCatalog} catalog */
const discoveryPending = (catalog) => (Array.isArray(catalog?.issues) ? catalog.issues : [])
  .some((issue) => issue?.context?.reason === 'DISCOVERY_PENDING');

/**
 * @typedef {object} AnnotationRecordSelectorOptions
 * @property {() => AnnotationRecordCatalog | null | undefined} getCatalog The current catalog of the mode.
 */

/**
 * @param {AnnotationRecordSelectorOptions} options
 */
export const createAnnotationRecordSelector = ({ getCatalog }) => {
  const catalog = () => getCatalog?.() || {
    status: 'error', records: [], issues: [], requiresSelection: false
  };
  /** @param {AnnotationRow | null | undefined} annotation */
  const isMissing = (annotation) => {
    const current = catalog();
    if (current.status !== 'ready') return false;
    const parsed = selectorFromTarget(annotation);
    if (parsed.error) return true;
    if (!parsed.selector) return Boolean(current.requiresSelection);
    return !resolvedRecord(current, annotation);
  };
  // A single-record input needs no record select: the row names
  // its record in a caption, and only several records or a missing target
  // show the select.
  /** @param {AnnotationRow} annotation */
  const showsSelect = (annotation) => {
    const current = catalog();
    return current.status === 'ready' && (current.records.length > 1 || isMissing(annotation));
  };
  /** @param {AnnotationRow} annotation */
  const captionFor = (annotation) => {
    if (showsSelect(annotation)) return '';
    const current = catalog();
    const named = bindingFor(annotation) || annotation?.target?.record;
    if (!named) return '';
    const record = current.status === 'ready' ? resolvedRecord(current, annotation) : null;
    const name = record?.recordId || annotationRecordSelectorValue(annotation?.target?.record);
    if (!name) return '';
    const unread = current.status !== 'ready' && !discoveryPending(current);
    return unread ? t('{name} (record list unavailable)', { name }) : name;
  };
  return {
    showsSelect,
    captionFor,
    /** @param {AnnotationRow | null | undefined} annotation */
    optionsFor: (annotation) => annotationRecordOptions(catalog(), annotation),
    /** @param {AnnotationRow | null | undefined} annotation */
    valueFor: (annotation) => annotationRecordValue(catalog(), annotation),
    /**
     * @param {AnnotationSetsSlice} slice The handle of the annotation's drawing.
     * @param {AnnotationRow | null | undefined} annotation
     * @param {unknown} value
     */
    setValue: (slice, annotation, value) => setAnnotationRecordValue(catalog(), slice, annotation, value),
    isRequired: () => Boolean(catalog().requiresSelection),
    isMissing,
    missingMessage: () => t('Choose the record that this annotation targets.')
  };
};
