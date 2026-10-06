// @ts-check
import { diagnosticError } from '../../services/error-normalization.js';
import { featureIdentityKeyOf, rowBelongsToRequest } from '../../services/feature-placement.js';

/**
 * @import { AnnotationRecordSelector } from './target-actions.js'
 */

/**
 * @typedef {{ angle: number, spacing: number, color: string, width: number, cross: boolean }} AnnotationHatch
 * @typedef {{
 *   stroke: string, strokeWidth: number, strokeDasharray: number[], lineCap: string,
 *   fill: string | null, fillOpacity: number, hatch: AnnotationHatch | null, labelColor: string,
 *   labelFontSize: number | null, labelOrientation: string, labelPosition: string, labelOffset: number
 * }} AnnotationStyle
 * @typedef {{
 *   kind: string, record?: AnnotationRecordSelector | null, envelope?: string, circularPath?: string,
 *   start?: number, end?: number, coordinateSpace?: string, wrapsOrigin?: boolean, outOfBounds?: string,
 *   selectors?: { key: string | null, value: string }[],
 *   scope?: string, recordKey?: string, biologicalFeatureId?: string
 * }} AnnotationTarget One shape for `coordinateSpan`, `featureSpan`, and `featureIdentity`; each kind reads its own fields.
 * @typedef {{
 *   id: string, target: AnnotationTarget, label: string, mark: string, lane: number | null,
 *   style: AnnotationStyle | null, legendLabel: string | null, metadata: Record<string, any>
 * }} AnnotationItem
 * @typedef {{
 *   id: string, annotations: AnnotationItem[], defaultStyle: AnnotationStyle, legendLabel: string | null
 * }} AnnotationSet
 */

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
const normalizeEnvelope = (value) => value === 'segments' ? 'segments' : 'outer_bounds';
const normalizeCircularPath = (value) => ['forward', 'reverse'].includes(value) ? value : 'shortest';
const normalizeTarget = (target) => {
  const source = target && typeof target === 'object' ? clone(target) : {};
  // One feature named by its original-source identity (request schema 9,
  // design Q4); Python resolves it after crop and reverse complement. In the
  // draft it also names the mode it was selected in (`scope`), as per-feature
  // edits do (R2). The editor writes exact targets, so a malformed one can come
  // only from a Session file.
  if (source.kind === 'featureIdentity') {
    if (!featureIdentityKeyOf(source)) throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
    return {
      kind: 'featureIdentity',
      scope: source.scope,
      recordKey: source.recordKey,
      biologicalFeatureId: source.biologicalFeatureId,
      envelope: normalizeEnvelope(source.envelope),
      circularPath: normalizeCircularPath(source.circularPath)
    };
  }
  if (source.kind === 'featureSpan') {
    return {
      kind: 'featureSpan',
      record: source.record ?? null,
      selectors: Array.isArray(source.selectors) ? source.selectors.map((selector) => ({
        key: selector?.key == null || selector.key === '' ? null : String(selector.key),
        value: String(selector?.value || '')
      })) : [],
      envelope: normalizeEnvelope(source.envelope),
      circularPath: normalizeCircularPath(source.circularPath)
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

/**
 * @param {Partial<AnnotationStyle> | null} [overrides]
 * @returns {AnnotationStyle}
 */
export const createDefaultAnnotationStyle = (overrides = {}) => ({
  ...DEFAULT_STYLE,
  ...(overrides && typeof overrides === 'object' ? clone(overrides) : {})
});

/**
 * @param {Partial<AnnotationSet>} [overrides] A set from a Session, a table, or the editor.
 * @returns {AnnotationSet}
 */
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

// A selected-feature target names its mode and record by key: a request
// carries only the targets of its own mode and records, without the draft-only
// `scope`, and the others stay in the draft for the records and mode that draw
// them (design Q4 3.2, R2).
export const annotationOptionsPayload = (sets, mode, records = []) => ({
  sets: normalizeAnnotationSets(sets).map((set) => ({
    ...set,
    annotations: set.annotations.flatMap((item) => {
      if (item.target.kind !== 'featureIdentity') return [item];
      const { scope: _scope, ...target } = item.target;
      return rowBelongsToRequest(item.target, mode, records) ? [{ ...item, target }] : [];
    })
  })),
  table: null,
  tableFile: null
});

// The draft sets of a request's annotation sets, whose selected-feature targets
// are of the request's mode.
export const draftAnnotationSetsOfRequest = (sets, mode) => normalizeAnnotationSets(
  (Array.isArray(sets) ? sets : []).map((set) => ({
    ...set,
    annotations: (Array.isArray(set?.annotations) ? set.annotations : []).map((item) => (
      item?.target?.kind === 'featureIdentity' ? { ...item, target: { scope: mode, ...item.target } } : item
    ))
  }))
);
