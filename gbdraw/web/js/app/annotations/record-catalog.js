// @ts-check
/**
 * @import { AnnotationRecordSelector } from '../../services/annotation-state.js'
 */
import {
  annotationRecordSelectorFromValue,
  parseAnnotationRecordSelectorValue
} from './target-actions.js';
import { buildDisambiguatedRecordEntries, formatRecordLength } from '../../services/record-options.js';
import { normalizeUserFacingError } from '../../utils/error-normalization.js';

/**
 * @typedef {{
 *   key: string, value: string, label: string, recordId: string, recordLength: number | null,
 *   index: number, localIndex: number, sourceIndex?: number, sourceKey: string,
 *   backendSelector: AnnotationRecordSelector | null
 * }} AnnotationCatalogRecord A record a target may name; `key` survives reorders, `value` is the record ID or `#<n>`.
 * @typedef {{ code: string, context: Record<string, any> }} AnnotationCatalogIssue A producer diagnostic.
 * @typedef {{
 *   mode?: string, status: string, records: AnnotationCatalogRecord[],
 *   issues: AnnotationCatalogIssue[], requiresSelection: boolean, signature?: string
 * }} AnnotationRecordCatalog `status` is 'ready', 'loading', or 'error'.
 * @typedef {{
 *   sourceKey?: string, hasInput: boolean, status: string, error?: any, selector?: string,
 *   records?: Record<string, any>[]
 * }} AnnotationCatalogSource A discovered input: the Circular discovery state or a Linear source.
 *   Each record carries `recordId` or `record_id`, `recordLength` or `record_length`, and `selector`.
 * @typedef {object} AnnotationRecordCatalogOptions
 * @property {string} [mode] 'linear' selects `linearSources`; anything else `circularSource`.
 * @property {AnnotationCatalogSource | null} [circularSource]
 * @property {AnnotationCatalogSource[]} [linearSources]
 */

const cleanText = (value) => String(value ?? '').trim();
// Catalog issues are producer diagnostics ({ code, context }); wording stays
// with the error normalizer.
const catalogIssue = (code, context = {}) => ({ code, context });
const discoveryIssue = (error, context = {}) => {
  const model = normalizeUserFacingError(error || null);
  return model && model.code !== 'UNKNOWN'
    ? catalogIssue(model.code, { ...model.context, ...context })
    : catalogIssue('INPUT_UNREADABLE', context);
};
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
  if (parsed.error) return { records: [], error: 'SELECTOR_FORMAT', explicit: true };
  const selector = parsed.selector;
  if (!selector) return { records, error: '', explicit: false };
  if (selector.kind === 'recordIndex') {
    const selected = records.find((record) => record.localIndex === selector.index);
    return {
      records: selected ? [selected] : [],
      error: selected ? '' : 'OUT_OF_RANGE',
      explicit: true
    };
  }
  const matches = records.filter((record) => record.recordId === selector.value);
  if (matches.length !== 1) {
    return {
      records: [],
      error: matches.length > 1 ? 'AMBIGUOUS' : 'NO_MATCH',
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

const buildLinearCatalog = (sources) => {
  const normalizedSources = Array.isArray(sources) ? sources : [];
  const records = [];
  const issues = [];
  normalizedSources.forEach((source, sourceIndex) => {
    const inputOrdinal = sourceIndex + 1;
    if (!source?.hasInput) {
      issues.push(catalogIssue('INPUT_REQUIRED', { inputOrdinal }));
      return;
    }
    if (source?.status !== 'ready') {
      issues.push(source?.status === 'error'
        ? discoveryIssue(source?.error, { inputOrdinal })
        : catalogIssue('RECORD_SELECTION', { inputOrdinal, reason: 'DISCOVERY_PENDING' }));
      return;
    }
    const availableRecords = sourceRecords(source, `linear-source-${sourceIndex + 1}`);
    if (availableRecords.length === 0) {
      issues.push(catalogIssue('NO_RECORDS', { inputOrdinal }));
      return;
    }
    const selected = selectSourceRecords(availableRecords, source?.selector);
    if (selected.error) {
      issues.push(catalogIssue('RECORD_SELECTION', { inputOrdinal, reason: selected.error }));
      return;
    }
    records.push(...selected.records.map((record) => ({ ...record, sourceIndex })));
  });
  const finalized = finalizeRecords(records);
  return {
    mode: 'linear',
    status: catalogStatus(issues, normalizedSources),
    records: finalized,
    issues,
    requiresSelection: finalized.length > 1,
    signature: finalized.map((record) => record.key).join('|')
  };
};

// circularSource.records are the records the Circular request draws
// (resolveCircularRequestRecordSet); several records need an explicit target.
const buildCircularCatalog = (source) => {
  const issues = [];
  if (!source?.hasInput) {
    issues.push(catalogIssue('INPUT_REQUIRED'));
  } else if (source?.status !== 'ready') {
    issues.push(source?.status === 'error'
      ? discoveryIssue(source?.error)
      : catalogIssue('RECORD_SELECTION', { reason: 'DISCOVERY_PENDING' }));
  }
  const records = issues.length === 0 ? finalizeRecords(sourceRecords(source, 'circular-source')) : [];
  return {
    mode: 'circular',
    status: issues.length === 0 ? 'ready' : (source?.status === 'loading' ? 'loading' : 'error'),
    records,
    issues,
    requiresSelection: records.length > 1,
    signature: records.map((record) => record.key).join('|')
  };
};

/**
 * @param {AnnotationRecordCatalogOptions} [options]
 * @returns {AnnotationRecordCatalog}
 */
export const buildAnnotationRecordCatalog = ({
  mode,
  circularSource = null,
  linearSources = []
} = {}) => (
  cleanText(mode) === 'linear'
    ? buildLinearCatalog(linearSources)
    : buildCircularCatalog(circularSource)
);
