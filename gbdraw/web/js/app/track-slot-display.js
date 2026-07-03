const normalizeOptionalText = (value) => {
  const text = String(value ?? '').trim();
  return text.length > 0 ? text : null;
};

const roundDisplayNumber = (value, digits = 1) => {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return '';
  const fixed = numeric.toFixed(digits);
  return fixed.replace(/\.0+$/, '').replace(/(\.\d*?)0+$/, '$1');
};

export const isManualSlotValue = (value) => normalizeOptionalText(value) !== null;

export const formatPxAuto = (value, source = 'estimated') => {
  void source;
  const text = roundDisplayNumber(value, 1);
  return text ? `${text} px (auto)` : '';
};

export const formatRadiusFactorAuto = (value, source = 'estimated') => {
  void source;
  const text = roundDisplayNumber(value, 2);
  return text ? `${text} R (auto)` : '';
};

export const parseCircularScalarDisplay = (value) => {
  const text = normalizeOptionalText(value);
  if (text === null) return null;
  if (/px$/i.test(text)) {
    const numeric = Number(text.slice(0, -2));
    return Number.isFinite(numeric) ? { value: numeric, unit: 'px', text } : null;
  }
  const numeric = Number(text);
  return Number.isFinite(numeric) ? { value: numeric, unit: 'factor', text } : null;
};

export const formatCircularWidthValue = (value) => {
  const parsed = parseCircularScalarDisplay(value);
  if (!parsed) return '';
  return parsed.unit === 'px' ? `${roundDisplayNumber(parsed.value, 1)} px` : `${roundDisplayNumber(parsed.value, 3)} R`;
};

export const manualOrAutoDisplay = ({ manualValue, autoValue, formatter }) => {
  const manualText = normalizeOptionalText(manualValue);
  if (manualText !== null) return manualText;
  return typeof formatter === 'function' ? formatter(autoValue) : String(autoValue ?? '');
};

export const slotGeometryInstanceKey = (resultIndex, recordIndex, slotIndex) =>
  `${Number(resultIndex) || 0}:${Number(recordIndex) || 0}:${Number(slotIndex) || 0}`;

export const findTrackSlotGeometry = ({
  geometry,
  resultIndex = 0,
  recordIndex = 0,
  slotIndex
} = {}) => {
  if (!geometry || typeof geometry !== 'object' || !Array.isArray(geometry.records)) return null;
  const wantedResult = Number(resultIndex) || 0;
  const wantedRecord = Number(recordIndex) || 0;
  const wantedSlot = Number(slotIndex);
  if (!Number.isInteger(wantedSlot)) return null;
  const matchingRecord = geometry.records.find((record) => (
    Number(record?.resultIndex ?? 0) === wantedResult &&
    Number(record?.recordIndex ?? 0) === wantedRecord
  )) || geometry.records.find((record) => Number(record?.resultIndex ?? 0) === wantedResult);
  if (!matchingRecord || !Array.isArray(matchingRecord.slots)) return null;
  return matchingRecord.slots.find((slot) => Number(slot?.slotIndex) === wantedSlot) || null;
};
