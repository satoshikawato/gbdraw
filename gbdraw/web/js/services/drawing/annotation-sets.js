// @ts-check
// The annotation-sets slice of a drawing: its annotation sets and their rows
// (Session 46 `modes.<m>.config.annotationSets`). The drawing shows the sets
// read-only to the type checker, rows and fields included; only this module
// writes them, through `drawing.slices.annotationSets`: the annotation
// editor's writers below (app/annotations.js), the Delete settings of a
// record, and Load, a History restore, and Reset (`install`, `reset`).
// A set is the set object of the list, an annotation the row object of its
// set (D-15-44): a writer given another object changes nothing.
import { reactive } from '../../utils/reactivity.js';
import { cloneJsonData } from '../json-clone.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';
/** @import { DeepReadonly } from '../../utils/reactivity.js' */
/** @import { SliceHandle } from './slice-handle.js' */
/** @import { AnnotationItem, AnnotationRecordSelector, AnnotationSet, AnnotationStyle, AnnotationTarget } from '../annotation-state.js' */

/** @typedef {{ annotationSets: AnnotationSet[] }} AnnotationSetsData What `capture` returns and `install` writes. */
/** @typedef {DeepReadonly<AnnotationSet>} AnnotationSetRow A set as the drawing shows it. */
/** @typedef {DeepReadonly<AnnotationItem>} AnnotationRow An annotation as the drawing shows it. */
/** @typedef {{ readonly annotationSets: readonly AnnotationSetRow[] }} AnnotationSetsView */
/** @typedef {SliceHandle<'annotationSets'>} AnnotationSetsSlice The handle a writer takes. */
/**
 * A record binding of an annotation (app/annotations/record-selector.js): the
 * record its target names, and its metadata with the binding.
 * @typedef {{ annotation: AnnotationRow, record: AnnotationRecordSelector | null, metadata: Record<string, unknown> }} AnnotationRecordBinding
 */
/**
 * A value a number input gives (Vue's `v-model.number`): a number, or the
 * text when it reads as none.
 * @typedef {number | string} NumberInputValue
 */

// An annotation-sets handle exists only through `create`, which gives it these members.
/** @param {AnnotationSetsSlice} slice */
const membersOf = (slice) => /** @type {AnnotationSetsData} */ (sliceMembers(slice, 'annotationSets'));

/**
 * The set object of the list that `set` is.
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSetRow | null | undefined} set
 */
const setOf = (slice, set) => membersOf(slice).annotationSets.find((candidate) => candidate === set);

/**
 * The row object of `set` that `item` is.
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSetRow | null | undefined} set
 * @param {AnnotationRow | null | undefined} item
 */
const rowOf = (slice, set, item) => setOf(slice, set)?.annotations.find((candidate) => candidate === item);

/**
 * A new slice, or one over a drawing's existing members (a test fixture).
 * @param {Partial<AnnotationSetsData>} [members]
 * @returns {{ handle: AnnotationSetsSlice, view: AnnotationSetsView }}
 */
const create = (members = {}) => {
  const data = { annotationSets: members.annotationSets ?? reactive(/** @type {AnnotationSet[]} */ ([])) };
  return { handle: createSliceHandle('annotationSets', data), view: Object.freeze({ annotationSets: data.annotationSets }) };
};

/**
 * @param {AnnotationSetsSlice} slice
 * @returns {AnnotationSetsData}
 */
const capture = (slice) => ({ annotationSets: cloneJsonData([...membersOf(slice).annotationSets]) });

// The list becomes a copy of the sets `data` holds; without them it is kept.
/**
 * @param {AnnotationSetsSlice} slice
 * @param {{ annotationSets?: readonly AnnotationSet[] }} data
 */
const install = (slice, { annotationSets }) => {
  if (!annotationSets) return;
  const sets = cloneJsonData([...annotationSets]);
  const list = membersOf(slice).annotationSets;
  list.splice(0, list.length, ...sets);
};

/** @param {AnnotationSetsSlice} slice */
const reset = (slice) => install(slice, { annotationSets: [] });

/**
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSet} set A new set (the slice keeps the object).
 */
export const addAnnotationSet = (slice, set) => {
  membersOf(slice).annotationSets.push(set);
};

/**
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSetRow} set
 */
export const removeAnnotationSet = (slice, set) => {
  const list = membersOf(slice).annotationSets;
  const index = list.findIndex((candidate) => candidate === set);
  if (index >= 0) list.splice(index, 1);
};

/**
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSetRow} set
 * @param {string} id An id no other set of the list has.
 */
export const renameAnnotationSet = (slice, set, id) => {
  const row = setOf(slice, set);
  if (row) row.id = id;
};

/**
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSetRow} set
 * @param {string} legendLabel
 */
export const setAnnotationSetLegendLabel = (slice, set, legendLabel) => {
  const row = setOf(slice, set);
  if (row) row.legendLabel = legendLabel;
};

/**
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSetRow} set
 * @param {readonly AnnotationItem[]} items New rows (the slice keeps the objects).
 */
export const addAnnotations = (slice, set, items) => {
  setOf(slice, set)?.annotations.push(...items);
};

/**
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSetRow} set
 * @param {AnnotationRow} item
 */
export const removeAnnotation = (slice, set, item) => {
  const rows = setOf(slice, set)?.annotations;
  const index = rows ? rows.findIndex((candidate) => candidate === item) : -1;
  if (rows && index >= 0) rows.splice(index, 1);
};

/**
 * Removes the rows of the list that `items` holds, from their sets; a value
 * that is no row changes nothing.
 * @param {AnnotationSetsSlice} slice
 * @param {readonly object[]} items Rows, by identity.
 */
export const removeAnnotations = (slice, items) => {
  const removed = new Set(items);
  membersOf(slice).annotationSets.forEach(({ annotations }) => {
    const kept = annotations.filter((row) => !removed.has(row));
    if (kept.length < annotations.length) annotations.splice(0, annotations.length, ...kept);
  });
};

/**
 * An annotation's own fields.
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSetRow} set
 * @param {AnnotationRow} item
 * @param {{ id?: string, label?: string, mark?: string, lane?: NumberInputValue | null }} fields
 */
export const setAnnotationFields = (slice, set, item, fields) => {
  const row = rowOf(slice, set, item);
  // A lane is the number input's value as Vue's `v-model.number` gave it.
  if (row) Object.assign(row, /** @type {Partial<AnnotationItem>} */ (fields));
};

/**
 * One style field; a row without a style of its own starts from `defaults`.
 * @template {keyof AnnotationStyle} K
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSetRow} set
 * @param {AnnotationRow} item
 * @param {K} field
 * @param {AnnotationStyle[K]} value
 * @param {AnnotationStyle} defaults A new style object.
 */
export const setAnnotationStyleField = (slice, set, item, field, value, defaults) => {
  const row = rowOf(slice, set, item);
  if (!row) return;
  row.style ??= defaults;
  row.style[field] = value;
};

/**
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSetRow} set
 * @param {AnnotationRow} item
 * @param {AnnotationTarget} target A new target object.
 */
export const replaceAnnotationTarget = (slice, set, item, target) => {
  const row = rowOf(slice, set, item);
  if (row) row.target = target;
};

/**
 * The start or end of a coordinate target, as the number input gives it.
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSetRow} set
 * @param {AnnotationRow} item
 * @param {'start' | 'end'} field
 * @param {NumberInputValue} value
 */
export const setAnnotationTargetBound = (slice, set, item, field, value) => {
  const row = rowOf(slice, set, item);
  // As Vue's `v-model.number` wrote it.
  if (row) row.target[field] = /** @type {number} */ (value);
};

/**
 * The qualifier or value of a feature target's first selector.
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSetRow} set
 * @param {AnnotationRow} item
 * @param {'key' | 'value'} field
 * @param {string} value
 */
export const setAnnotationSelectorField = (slice, set, item, field, value) => {
  const selector = rowOf(slice, set, item)?.target.selectors?.[0];
  if (selector) selector[field] = value;
};

/**
 * Writes record bindings: each annotation's target names `record`, and its
 * metadata becomes `metadata` (a new object). An annotation of no set of the
 * list is skipped.
 * @param {AnnotationSetsSlice} slice
 * @param {readonly AnnotationRecordBinding[]} bindings
 */
export const bindAnnotationRecords = (slice, bindings) => {
  const sets = membersOf(slice).annotationSets;
  bindings.forEach(({ annotation, record, metadata }) => {
    const row = sets.flatMap((set) => set.annotations).find((candidate) => candidate === annotation);
    if (!row) return;
    row.target.record = record && { ...record };
    row.metadata = metadata;
  });
};

export const ANNOTATION_SETS_SLICE = Object.freeze({ members: Object.freeze(['annotationSets']), create, capture, install, reset });
