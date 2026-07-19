const cleanText = (value) => String(value ?? '').trim();
const NULL_RECORD_SELECTOR_TOKENS = new Set([
  'none', 'null', 'jsnull', 'undefined', 'jsundefined', '-'
]);

export const isUnspecifiedRecordSelectorValue = (value) => {
  const text = cleanText(value);
  return !text || NULL_RECORD_SELECTOR_TOKENS.has(text.toLowerCase());
};

const isSafeRecordId = (value) => {
  const recordId = cleanText(value);
  return !isUnspecifiedRecordSelectorValue(recordId) && !recordId.startsWith('#');
};

export const formatRecordLength = (value) => {
  const numeric = Number(value);
  if (!Number.isInteger(numeric) || numeric <= 0) return 'length unavailable';
  return `${numeric.toLocaleString('en-US')} bp`;
};

export const buildDisambiguatedRecordEntries = (records) => {
  const normalized = (Array.isArray(records) ? records : []).map((record, index) => ({
    ...record,
    selector: cleanText(record?.selector) || `#${index + 1}`,
    recordId: cleanText(record?.recordId) || `Record_${index + 1}`
  }));
  const idCounts = new Map();
  normalized.forEach((record) => {
    idCounts.set(record.recordId, (idCounts.get(record.recordId) || 0) + 1);
  });
  return normalized.map((record) => {
    const duplicate = (idCounts.get(record.recordId) || 0) > 1;
    const usesIndex = duplicate || !isSafeRecordId(record.recordId);
    return {
      ...record,
      duplicate,
      usesIndex,
      value: usesIndex ? record.selector : record.recordId
    };
  });
};
