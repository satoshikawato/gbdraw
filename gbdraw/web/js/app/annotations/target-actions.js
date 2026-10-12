// @ts-check
import { isUnspecifiedRecordSelectorValue } from '../../services/record-options.js';
import { featureIdentityKeyOf } from '../../services/feature-placement.js';

/**
 * @import { AnnotationRecordSelector } from '../../services/annotation-state.js'
 * @typedef {string | { key?: unknown, value?: unknown }} FeatureSelectorInput One `key=value` selector of a feature target.
 * @typedef {Record<string, unknown>} SelectedFeatureRow A selected feature, by its snake_case or camelCase identity fields.
 */

/** @param {unknown} value */
const cleanNullable = (value) => {
  const text = String(value ?? '').trim();
  return text || null;
};

/**
 * @param {unknown} value A table cell, a record ID, or `#<n>`.
 * @returns {{ selector: AnnotationRecordSelector | null, error: string }}
 */
export const parseAnnotationRecordSelectorValue = (value) => {
  const text = cleanNullable(value);
  if (text === null || isUnspecifiedRecordSelectorValue(text)) {
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

/** @param {unknown} value */
export const isSafeRecordIdSelector = (value) => {
  const parsed = parseAnnotationRecordSelectorValue(value);
  return !parsed.error && parsed.selector?.kind === 'recordId';
};

/**
 * @param {unknown} recordId
 * @param {unknown} recordIndex
 * @returns {AnnotationRecordSelector | null}
 */
export const annotationRecordSelector = (recordId, recordIndex) => {
  const id = cleanNullable(recordId);
  if (id && isSafeRecordIdSelector(id)) return { kind: 'recordId', value: id };
  if (recordIndex == null || recordIndex === '') return null;
  const index = Number(recordIndex);
  if (Number.isInteger(index) && index >= 0) return { kind: 'recordIndex', index };
  return null;
};

/** @param {AnnotationRecordSelector | null | undefined} selector */
export const annotationRecordSelectorValue = (selector) => {
  if (selector?.kind === 'recordId') return cleanNullable(selector.value) || '';
  if (selector?.kind === 'recordIndex') {
    const index = Number(selector.index);
    return Number.isInteger(index) && index >= 0 ? `#${index + 1}` : '';
  }
  return '';
};

/** @param {unknown} value */
export const annotationRecordSelectorFromValue = (value) => {
  const parsed = parseAnnotationRecordSelectorValue(value);
  if (parsed.error) throw new Error(parsed.error);
  return parsed.selector;
};

/**
 * @param {{ start: unknown, end: unknown, recordId?: string | null, recordIndex?: number | null, coordinateSpace?: string }} params
 */
export const coordinateTarget = ({ start, end, recordId = null, recordIndex = null, coordinateSpace = 'source' }) => ({
  kind: 'coordinateSpan',
  record: annotationRecordSelector(recordId, recordIndex),
  start: Number(start),
  end: Number(end),
  coordinateSpace: coordinateSpace === 'local' ? 'local' : 'source',
  wrapsOrigin: Number(start) > Number(end),
  outOfBounds: 'clip'
});

/** @param {FeatureSelectorInput} value */
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

/**
 * The `.filter` cast reads `.value` of a string selector as undefined, which keeps the string.
 * @param {{ selector?: unknown, selectors?: FeatureSelectorInput[] | null, recordId?: string | null, recordIndex?: number | null, extent?: string, circularPath?: string }} params
 */
export const featureTarget = ({ selector, selectors = null, recordId = null, recordIndex = null, extent = 'outer_bounds', circularPath = 'shortest' }) => ({
  kind: 'featureSpan',
  record: annotationRecordSelector(recordId, recordIndex),
  selectors: (Array.isArray(selectors) ? selectors : String(selector || '').split(';')).filter((value) => String(/** @type {{ value?: unknown }} */ (value)?.value ?? value).trim()).map(parseFeatureSelector),
  envelope: extent === 'segments' ? 'segments' : 'outer_bounds',
  circularPath: ['forward', 'reverse'].includes(circularPath) ? circularPath : 'shortest'
});

// A selected feature is named by its original-source identity, which Python
// resolves after crop, reverse complement, reordering, and duplication (design
// Q4, OV-03); the target is of the drawing of the Result's mode (R2). Null when
// a selected feature has none: a Result without a feature catalog (a Session
// before 40 until Generate).
/** @param {unknown} features Selected feature rows; anything but an array selects none. */
export const featureTargetsFromSelection = (features) => {
  /** @type {SelectedFeatureRow[]} */
  const selected = Array.isArray(features) ? features : [];
  if (!selected.every((feature) => featureIdentityKeyOf(feature))) return null;
  return selected.map((feature) => ({
    kind: 'featureIdentity',
    recordKey: feature.record_key ?? feature.recordKey,
    biologicalFeatureId: feature.biological_feature_id ?? feature.biologicalFeatureId,
    envelope: 'outer_bounds',
    circularPath: 'shortest'
  }));
};
