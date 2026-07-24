const DEFAULT_STYLE = Object.freeze({
  stroke: '#404040',
  strokeWidth: 1.5,
  strokeDasharray: [],
  lineCap: 'tick',
  fill: '#94a3b8',
  fillOpacity: 0.2,
  hatch: null,
  labelColor: '#202020',
  labelFontSize: null,
  labelOrientation: 'auto',
  labelPosition: 'center',
  labelOffset: 4
});

const cleanId = (value, fallback) => String(value || '').trim() || fallback;
const clone = (value) => JSON.parse(JSON.stringify(value));
const normalizeTarget = (target) => {
  const source = target && typeof target === 'object' ? clone(target) : {};
  if (source.kind === 'featureSpan') {
    return {
      kind: 'featureSpan',
      record: source.record ?? null,
      selectors: Array.isArray(source.selectors) ? source.selectors.map((selector) => ({
        key: selector?.key == null || selector.key === '' ? null : String(selector.key),
        value: String(selector?.value || '')
      })) : [],
      envelope: source.envelope === 'segments' ? 'segments' : 'outer_bounds',
      circularPath: ['forward', 'reverse'].includes(source.circularPath) ? source.circularPath : 'shortest'
    };
  }
  const start = Math.max(1, Number(source.start) || 1);
  const end = Math.max(1, Number(source.end) || 1);
  return {
    kind: 'coordinateSpan',
    record: source.record ?? null,
    start,
    end,
    coordinateSpace: source.coordinateSpace === 'local' ? 'local' : 'source',
    wrapsOrigin: start > end,
    outOfBounds: ['skip', 'error'].includes(source.outOfBounds) ? source.outOfBounds : 'clip'
  };
};

export const createDefaultAnnotationStyle = (overrides = {}) => ({
  ...DEFAULT_STYLE,
  ...(overrides && typeof overrides === 'object' ? clone(overrides) : {})
});

export const createAnnotationSet = (overrides = {}) => ({
  id: cleanId(overrides.id, 'annotations'),
  annotations: Array.isArray(overrides.annotations) ? clone(overrides.annotations) : [],
  defaultStyle: createDefaultAnnotationStyle(overrides.defaultStyle),
  legendLabel: overrides.legendLabel == null ? null : String(overrides.legendLabel)
});

export const normalizeAnnotationSets = (sets) => {
  const usedSetIds = new Set();
  return (Array.isArray(sets) ? sets : []).map((rawSet, setIndex) => {
    const set = createAnnotationSet(rawSet);
    let setId = cleanId(set.id, `annotations_${setIndex + 1}`);
    while (usedSetIds.has(setId)) setId = `${setId}_${setIndex + 1}`;
    usedSetIds.add(setId);
    const usedItemIds = new Set();
    set.id = setId;
    set.annotations = set.annotations.map((rawItem, itemIndex) => {
      const item = rawItem && typeof rawItem === 'object' ? clone(rawItem) : {};
      let id = cleanId(item.id, `region_${itemIndex + 1}`);
      while (usedItemIds.has(id)) id = `${id}_${itemIndex + 1}`;
      usedItemIds.add(id);
      return {
        id,
        target: normalizeTarget(item.target),
        label: String(item.label || ''),
        mark: ['line', 'bracket', 'band', 'highlight'].includes(item.mark) ? item.mark : 'bracket',
        lane: item.lane == null || item.lane === '' ? null : Math.max(0, Number(item.lane) || 0),
        style: item.style == null ? null : createDefaultAnnotationStyle(item.style),
        legendLabel: item.legendLabel == null ? null : String(item.legendLabel),
        metadata: item.metadata && typeof item.metadata === 'object' ? clone(item.metadata) : {}
      };
    });
    return set;
  });
};

export const uniqueAnnotationSetId = (sets, base = 'annotations') => {
  const ids = new Set((Array.isArray(sets) ? sets : []).map((set) => String(set?.id || '')));
  const stem = cleanId(base, 'annotations');
  let id = stem;
  let suffix = 2;
  while (ids.has(id)) {
    id = `${stem}_${suffix}`;
    suffix += 1;
  }
  return id;
};

export const annotationOptionsPayload = (sets) => ({ sets: normalizeAnnotationSets(sets), table: null, tableFile: null });
