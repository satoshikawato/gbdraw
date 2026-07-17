const cleanNullable = (value) => {
  const text = String(value ?? '').trim();
  return text || null;
};

const recordSelector = (recordId, recordIndex) => {
  const id = cleanNullable(recordId);
  if (id) return { kind: 'recordId', value: id };
  if (recordIndex != null && recordIndex !== '') return { kind: 'recordIndex', index: Math.max(0, Number(recordIndex) || 0) };
  return null;
};

export const coordinateTarget = ({ start, end, recordId = null, recordIndex = null, coordinateSpace = 'source' }) => ({
  kind: 'coordinateSpan',
  record: recordSelector(recordId, recordIndex),
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
  record: recordSelector(recordId, recordIndex),
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
      recordIndex: options.recordIndex ?? feature?.record_index ?? null
    });
  })
);
