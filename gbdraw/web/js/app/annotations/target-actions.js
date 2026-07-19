import { isUnspecifiedRecordSelectorValue } from '../record-options.js';

const cleanNullable = (value) => {
  const text = String(value ?? '').trim();
  return text || null;
};

export const parseAnnotationRecordSelectorValue = (value) => {
  const text = cleanNullable(value);
  if (isUnspecifiedRecordSelectorValue(text)) {
    return { selector: null, error: '' };
  }
  if (!text.startsWith('#')) {
    return { selector: { kind: 'recordId', value: text }, error: '' };
  }
  const indexText = text.slice(1).trim();
  if (!/^[0-9]+$/.test(indexText)) {
    return {
      selector: null,
      error: `Invalid record selector ${JSON.stringify(text)}. Use #<number> or record_id.`
    };
  }
  const index = Number(indexText) - 1;
  if (!Number.isSafeInteger(index) || index < 0) {
    return {
      selector: null,
      error: `Record index must be >= 1 in selector ${JSON.stringify(text)}.`
    };
  }
  return { selector: { kind: 'recordIndex', index }, error: '' };
};

export const isSafeRecordIdSelector = (value) => {
  const parsed = parseAnnotationRecordSelectorValue(value);
  return !parsed.error && parsed.selector?.kind === 'recordId';
};

export const annotationRecordSelector = (recordId, recordIndex) => {
  const id = cleanNullable(recordId);
  if (id && isSafeRecordIdSelector(id)) return { kind: 'recordId', value: id };
  if (recordIndex == null || recordIndex === '') return null;
  const index = Number(recordIndex);
  if (Number.isInteger(index) && index >= 0) return { kind: 'recordIndex', index };
  return null;
};

export const annotationRecordSelectorValue = (selector) => {
  if (selector?.kind === 'recordId') return cleanNullable(selector.value) || '';
  if (selector?.kind === 'recordIndex') {
    const index = Number(selector.index);
    return Number.isInteger(index) && index >= 0 ? `#${index + 1}` : '';
  }
  return '';
};

export const annotationRecordSelectorFromValue = (value) => {
  const parsed = parseAnnotationRecordSelectorValue(value);
  if (parsed.error) throw new Error(parsed.error);
  return parsed.selector;
};

export const coordinateTarget = ({ start, end, recordId = null, recordIndex = null, coordinateSpace = 'source' }) => ({
  kind: 'coordinateSpan',
  record: annotationRecordSelector(recordId, recordIndex),
  start: Number(start),
  end: Number(end),
  coordinateSpace: coordinateSpace === 'local' ? 'local' : 'source',
  wrapsOrigin: Number(start) > Number(end),
  outOfBounds: 'clip'
});

const parseFeatureSelector = (value) => {
  if (value && typeof value === 'object') {
    return { key: value.key == null ? null : String(value.key), value: String(value.value || '') };
  }
  const text = String(value || '').trim();
  const split = text.indexOf('=');
  return split > 0
    ? { key: text.slice(0, split), value: text.slice(split + 1) }
    : { key: null, value: text };
};

export const featureTarget = ({ selector, selectors = null, recordId = null, recordIndex = null, extent = 'outer_bounds', circularPath = 'shortest' }) => ({
  kind: 'featureSpan',
  record: annotationRecordSelector(recordId, recordIndex),
  selectors: (Array.isArray(selectors) ? selectors : String(selector || '').split(';')).filter((value) => String(value?.value ?? value).trim()).map(parseFeatureSelector),
  envelope: extent === 'segments' ? 'segments' : 'outer_bounds',
  circularPath: ['forward', 'reverse'].includes(circularPath) ? circularPath : 'shortest'
});

export const featureTargetsFromSelection = (features, options = {}) => (
  (Array.isArray(features) ? features : []).map((feature, index) => {
    const qualifier = feature?.qualifiers || {};
    const first = (value) => Array.isArray(value) ? value[0] : value;
    const locusTag = first(qualifier.locus_tag ?? feature?.locus_tag);
    const gene = first(qualifier.gene ?? feature?.gene);
    const selector = locusTag ? `locus_tag=${locusTag}` : (gene ? `gene=${gene}` : `#${index + 1}`);
    return featureTarget({
      ...options,
      selector,
      recordId: options.recordId ?? feature?.record_id ?? null,
      recordIndex: options.recordIndex ?? feature?.record_index ?? feature?.record_idx ?? feature?.fileIdx ?? null
    });
  })
);
