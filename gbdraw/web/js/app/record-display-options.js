// Surface preparation for the shared display/placement writer integration.
// These adapters do not enable controls or write a schema-6 request/session.
import { resolveDisambiguatedRecordSelection } from './record-options.js';

export const recordDisplayKey = ({ scope, sourceUid, selector }) => {
  if (!['circular', 'linear'].includes(scope) || !sourceUid || !/^#[1-9]\d*$/.test(selector)) {
    throw new Error('Record display requires a scope, source UID, and exact record selector.');
  }
  return JSON.stringify([scope, sourceUid, selector]);
};

export const buildRecordDisplayRows = ({ scope, sourceUid, source, records, selector = '' }) => {
  const selection = resolveDisambiguatedRecordSelection(records, selector);
  if (!['unspecified', 'resolved'].includes(selection.status)) {
    throw new Error(`Record selector is ${selection.status}: ${selector}.`);
  }
  return (selection.record ? [selection.record] : selection.entries).map((record) => {
    const row = { ...record, scope, sourceUid, source };
    return { ...row, key: recordDisplayKey(row) };
  });
};

export const reconcileRecordDisplayDrafts = (drafts, discoveredRows, replacedSourceUids = []) => {
  // Pass all discovered rows, including inactive selectors. A source replacement
  // purges its drafts even when the existing source card keeps its UID.
  const keys = new Set(discoveredRows.map(recordDisplayKey));
  return drafts.filter((draft) => !replacedSourceUids.includes(draft.sourceUid)
    && keys.has(recordDisplayKey(draft)));
};

export const parseRecordDisplayStart = (value) => {
  if (value === null || (typeof value === 'string' && !value.trim())) return null;
  if (!['string', 'number'].includes(typeof value)) throw new Error('Display start must be an integer.');
  const number = Number(value);
  if (!Number.isSafeInteger(number) || number < 1) {
    throw new Error('Display start must be a positive integer.');
  }
  return number;
};

export const recordDisplaySurface = (row, draft = {}, { cropped = false, reverse = false } = {}) => {
  const override = draft.topologyOverride ?? null;
  if (override !== null && typeof override !== 'boolean') throw new Error('Topology override must be boolean or null.');
  const effectiveCircular = override ?? row.detectedTopology === 'circular';
  const lengthKnown = Number.isSafeInteger(row.recordLength) && row.recordLength > 0;
  return {
    effectiveCircular,
    startEnabled: effectiveCircular && lengthKnown && !cropped,
    disabledReason: cropped ? 'Display start is unavailable for a cropped record.'
      : !lengthKnown ? 'Record length is unavailable.'
        : !effectiveCircular ? 'Display start requires a circular record.' : '',
    currentStart: reverse ? (lengthKnown ? row.recordLength : null) : 1
  };
};

export const selectedFeatureDisplayStart = ({ row, committedRow, selectedFeatures, shortcut }) => {
  if (!committedRow || recordDisplayKey(row) !== recordDisplayKey(committedRow)
    || !row.source || row.source !== committedRow.source
    || row.recordLength !== committedRow.recordLength || !committedRow.recordKey
    || selectedFeatures.length !== 1) {
    throw new Error('Select one feature bound to the current source record.');
  }
  const feature = selectedFeatures[0];
  if (feature?.record_key !== committedRow.recordKey) throw new Error('Selected feature belongs to another record.');
  const parts = feature.location_parts;
  const strand = feature.strand;
  if (!['+', '-'].includes(strand) || !Array.isArray(parts) || !parts.length
    || !Number.isSafeInteger(row.recordLength) || row.recordLength <= 0
    || parts.some((part) => part.strand !== strand
      || !Number.isSafeInteger(part.start) || !Number.isSafeInteger(part.end)
      || part.start < 0 || part.end <= part.start || part.end > row.recordLength)) {
    throw new Error('Feature shortcut requires nonempty source parts with one known strand.');
  }
  if (!['five-prime', 'midpoint'].includes(shortcut)) throw new Error('Unknown feature shortcut.');
  const length = parts.reduce((sum, part) => sum + part.end - part.start, 0);
  let offset = shortcut === 'five-prime' ? 0 : Math.floor((length - 1) / 2);
  for (const part of parts) {
    if (offset < part.end - part.start) {
      return strand === '+' ? part.start + 1 + offset : part.end - offset;
    }
    offset -= part.end - part.start;
  }
};
