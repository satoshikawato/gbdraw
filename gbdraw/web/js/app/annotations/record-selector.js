import {
  annotationRecordSelectorValue,
  parseAnnotationRecordSelectorValue
} from './target-actions.js';

export const ANNOTATION_RECORD_BINDING_KEY = '_gbdraw_web_target_record_key';

const cleanText = (value) => String(value ?? '').trim();

const bindingFor = (annotation) => cleanText(annotation?.metadata?.[ANNOTATION_RECORD_BINDING_KEY]);
export const annotationRecordBinding = bindingFor;

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

const recordForSelector = (records, selector) => {
  if (!selector) return null;
  if (selector.kind === 'recordIndex') return records[selector.index] || null;
  const matches = records.filter((record) => record.recordId === selector.value);
  return matches.length === 1 ? matches[0] : null;
};

const resolvedRecord = (catalog, annotation) => {
  const records = Array.isArray(catalog?.records) ? catalog.records : [];
  const binding = bindingFor(annotation);
  if (binding) return records.find((record) => record.key === binding) || null;
  const parsed = selectorFromTarget(annotation);
  return parsed.error ? null : recordForSelector(records, parsed.selector);
};
export const resolveAnnotationRecord = resolvedRecord;

const setBinding = (annotation, key) => {
  if (!annotation || typeof annotation !== 'object') return;
  if (!annotation.metadata || typeof annotation.metadata !== 'object' || Array.isArray(annotation.metadata)) {
    annotation.metadata = {};
  }
  if (key) annotation.metadata[ANNOTATION_RECORD_BINDING_KEY] = key;
  else delete annotation.metadata[ANNOTATION_RECORD_BINDING_KEY];
};

export const annotationRecordValue = (catalog, annotation) => {
  const binding = bindingFor(annotation);
  if (binding) return binding;
  const record = resolvedRecord(catalog, annotation);
  if (record) return record.key;
  const saved = annotationRecordSelectorValue(annotation?.target?.record);
  return saved ? `saved:${saved}` : '';
};

export const setAnnotationRecordValue = (catalog, annotation, value) => {
  if (!annotation?.target) return;
  const key = cleanText(value);
  if (!key) {
    annotation.target.record = null;
    setBinding(annotation, '');
    return;
  }
  const record = (Array.isArray(catalog?.records) ? catalog.records : [])
    .find((candidate) => candidate.key === key);
  if (!record) return;
  annotation.target.record = { ...record.backendSelector };
  setBinding(annotation, record.key);
};

export const reconcileAnnotationRecordBindings = (sets, catalog) => {
  const records = Array.isArray(catalog?.records) ? catalog.records : [];
  (Array.isArray(sets) ? sets : []).forEach((set) => {
    (Array.isArray(set?.annotations) ? set.annotations : []).forEach((annotation) => {
      if (!annotation?.target) return;
      const binding = bindingFor(annotation);
      const record = binding
        ? records.find((candidate) => candidate.key === binding)
        : resolvedRecord(catalog, annotation);
      if (!record) return;
      annotation.target.record = { ...record.backendSelector };
      setBinding(annotation, record.key);
    });
  });
};

export const annotationRecordOptions = (catalog, annotation) => {
  const records = catalog?.allowExplicitSelectors === false
    ? []
    : (Array.isArray(catalog?.records) ? catalog.records : []);
  const currentValue = annotationRecordValue(catalog, annotation);
  const options = records.map((record) => ({
    value: record.key,
    label: record.label,
    synthetic: false
  }));
  const selectedIsMissing = Boolean(currentValue) && !options.some((option) => option.value === currentValue);
  const savedSelector = annotationRecordSelectorValue(annotation?.target?.record) || 'saved target';
  return [
    {
      value: '',
      label: catalog?.allowExplicitSelectors === false
        ? 'Automatic (current output record)'
        : (catalog?.status !== 'ready'
            ? 'Record list unavailable'
            : (records.length > 1 ? 'Select target record' : 'Automatic (single record)')),
      synthetic: false
    },
    ...(selectedIsMissing
      ? [{ value: currentValue, label: `${savedSelector} (target record unavailable)`, synthetic: true }]
      : []),
    ...options
  ];
};

export const createAnnotationRecordSelector = ({ getCatalog }) => {
  const catalog = () => getCatalog?.() || {
    status: 'error', records: [], issues: [], requiresSelection: false, allowExplicitSelectors: true
  };
  const isMissing = (annotation) => {
    const current = catalog();
    if (current.status !== 'ready') return false;
    const parsed = selectorFromTarget(annotation);
    if (parsed.error) return true;
    if (!parsed.selector) return Boolean(current.requiresSelection);
    if (current.allowExplicitSelectors === false) return true;
    return !resolvedRecord(current, annotation);
  };
  return {
    optionsFor: (annotation) => annotationRecordOptions(catalog(), annotation),
    valueFor: (annotation) => annotationRecordValue(catalog(), annotation),
    setValue: (annotation, value) => setAnnotationRecordValue(catalog(), annotation, value),
    isRequired: () => Boolean(catalog().requiresSelection),
    isMissing,
    isDisabled: (annotation) => {
      const current = catalog();
      if (current.status !== 'ready') return true;
      return current.allowExplicitSelectors === false && annotation?.target?.record == null;
    },
    missingMessage: (annotation) => (
      catalog().allowExplicitSelectors === false && selectorFromTarget(annotation).selector
        ? 'Clear the target record or enable Multi-record canvas.'
        : 'Choose the record that this annotation targets.'
    )
  };
};
