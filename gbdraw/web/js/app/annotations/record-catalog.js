import {
  annotationRecordSelectorFromValue,
  parseAnnotationRecordSelectorValue
} from './target-actions.js';
import { buildDisambiguatedRecordEntries, formatRecordLength } from '../record-options.js';

const cleanText = (value) => String(value ?? '').trim();
const fileFingerprint = (file) => file
  ? [String(file.name || ''), Number(file.size || 0)]
  : null;

export const annotationSourceKey = ({
  scope,
  uid = '',
  inputType = 'gb',
  primaryFile = null,
  pairedFile = null
}) => JSON.stringify([
  cleanText(scope),
  cleanText(uid),
  cleanText(inputType),
  fileFingerprint(primaryFile),
  fileFingerprint(pairedFile)
]);

const normalizeRecord = (record, fallbackIndex, sourceKey) => {
  const length = Number(record?.recordLength ?? record?.record_length);
  const selector = cleanText(record?.selector) || `#${fallbackIndex + 1}`;
  const parsedSelector = parseAnnotationRecordSelectorValue(selector).selector;
  const localIndex = parsedSelector?.kind === 'recordIndex' ? parsedSelector.index : fallbackIndex;
  const recordId = cleanText(record?.recordId ?? record?.record_id) || `Record_${fallbackIndex + 1}`;
  const recordLength = Number.isInteger(length) && length > 0 ? length : null;
  return {
    sourceKey,
    localIndex,
    key: `${sourceKey}::${JSON.stringify([localIndex, recordId, recordLength])}`,
    recordId,
    recordLength
  };
};

const sourceRecords = (source, fallbackSourceKey) => (
  (Array.isArray(source?.records) ? source.records : [])
    .map((record, index) => normalizeRecord(
      record,
      index,
      cleanText(source?.sourceKey) || fallbackSourceKey
    ))
);

const selectSourceRecords = (records, selectorValue) => {
  const parsed = parseAnnotationRecordSelectorValue(selectorValue);
  if (parsed.error) return { records: [], error: parsed.error, explicit: true };
  const selector = parsed.selector;
  if (!selector) return { records, error: '', explicit: false };
  if (selector.kind === 'recordIndex') {
    const selected = records.find((record) => record.localIndex === selector.index);
    return {
      records: selected ? [selected] : [],
      error: selected ? '' : `record selector #${selector.index + 1} is out of range`,
      explicit: true
    };
  }
  const matches = records.filter((record) => record.recordId === selector.value);
  if (matches.length !== 1) {
    return {
      records: [],
      error: matches.length > 1
        ? `record ID ${JSON.stringify(selector.value)} is duplicated; choose a #index entry`
        : `record ID ${JSON.stringify(selector.value)} is not available`,
      explicit: true
    };
  }
  return { records: matches, error: '', explicit: true };
};

const finalizeRecords = (records) => {
  const indexed = records.map((record, index) => ({ ...record, selector: `#${index + 1}` }));
  return buildDisambiguatedRecordEntries(indexed).map((record, index) => ({
    ...record,
    index,
    backendSelector: annotationRecordSelectorFromValue(record.value),
    label: `#${index + 1} · ${record.recordId}${record.recordLength ? ` · ${formatRecordLength(record.recordLength)}` : ''}`
  }));
};

const catalogStatus = (issues, sources) => {
  if (issues.length === 0) return 'ready';
  if (sources.some((source) => source?.status === 'error')) return 'error';
  return sources.some((source) => source?.status === 'loading') ? 'loading' : 'error';
};

const buildLinearCatalog = (sources, inputType, loadComparison) => {
  const normalizedSources = Array.isArray(sources) ? sources : [];
  const records = [];
  const issues = [];
  normalizedSources.forEach((source, sourceIndex) => {
    const sourceLabel = `Sequence #${sourceIndex + 1}`;
    if (!source?.hasInput) {
      issues.push(`${sourceLabel}: upload the required sequence file(s).`);
      return;
    }
    if (source?.status !== 'ready') {
      issues.push(source?.status === 'error'
        ? `${sourceLabel}: ${cleanText(source?.error) || 'record discovery failed.'}`
        : `${sourceLabel}: wait for the record list to finish loading.`);
      return;
    }
    const selected = selectSourceRecords(sourceRecords(source, `linear-source-${sourceIndex + 1}`), source?.selector);
    if (selected.error) {
      issues.push(`${sourceLabel}: ${selected.error}.`);
      return;
    }
    let materialized = selected.records;
    if (!selected.explicit && loadComparison) {
      const truncate = cleanText(inputType) === 'gff' || normalizedSources.length > 1;
      if (truncate) materialized = materialized.slice(0, 1);
    }
    records.push(...materialized);
  });
  const finalized = finalizeRecords(records);
  return {
    mode: 'linear',
    status: catalogStatus(issues, normalizedSources),
    records: finalized,
    issues,
    requiresSelection: finalized.length > 1,
    allowExplicitSelectors: true,
    signature: finalized.map((record) => record.key).join('|')
  };
};

const buildCircularCatalog = (source, multiRecordCanvas) => {
  const issues = [];
  if (!source?.hasInput) {
    issues.push('Upload the required circular sequence file(s).');
  } else if (source?.status !== 'ready') {
    issues.push(source?.status === 'error'
      ? cleanText(source?.error) || 'Circular record discovery failed.'
      : 'Wait for the circular record list to finish loading.');
  }
  const records = issues.length === 0 ? finalizeRecords(sourceRecords(source, 'circular-source')) : [];
  const allowExplicitSelectors = Boolean(multiRecordCanvas) || records.length <= 1;
  return {
    mode: 'circular',
    status: issues.length === 0 ? 'ready' : (source?.status === 'loading' ? 'loading' : 'error'),
    records,
    issues,
    requiresSelection: Boolean(multiRecordCanvas) && records.length > 1,
    allowExplicitSelectors,
    signature: records.map((record) => record.key).join('|')
  };
};

export const buildAnnotationRecordCatalog = ({
  mode,
  inputType = 'gb',
  loadComparison = false,
  multiRecordCanvas = false,
  circularSource = null,
  linearSources = []
} = {}) => (
  cleanText(mode) === 'linear'
    ? buildLinearCatalog(linearSources, inputType, Boolean(loadComparison))
    : buildCircularCatalog(circularSource, Boolean(multiRecordCanvas))
);
