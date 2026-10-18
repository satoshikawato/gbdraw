// @ts-check
// The annotation-sets slice of a drawing: its list of annotation sets
// (Session 46 `modes.<m>.config.annotationSets`). The drawing shows the list
// read-only to the type checker; only this module adds, removes, or replaces
// sets, through `drawing.slices.annotationSets`. A set's own fields and rows
// are the annotation editor's (app/annotations.js), as the template edits
// them in place.
import { reactive } from '../../utils/reactivity.js';
import { cloneJsonData } from '../json-clone.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';
/** @import { SliceHandle } from './slice-handle.js' */
/** @import { AnnotationSet } from '../annotation-state.js' */

/** @typedef {{ annotationSets: AnnotationSet[] }} AnnotationSetsData What `capture` returns and `install` writes. */
/** @typedef {{ readonly annotationSets: readonly AnnotationSet[] }} AnnotationSetsView */
/** @typedef {SliceHandle<'annotationSets'>} AnnotationSetsSlice The handle a writer takes. */

// An annotation-sets handle exists only through `create`, which gives it these members.
/** @param {AnnotationSetsSlice} slice */
const membersOf = (slice) => /** @type {AnnotationSetsData} */ (sliceMembers(slice));

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

// The list becomes the sets `data` holds (the objects themselves: a caller
// passes sets it read or reconciled).
/**
 * @param {AnnotationSetsSlice} slice
 * @param {{ annotationSets: readonly AnnotationSet[] }} data
 */
const install = (slice, { annotationSets }) => {
  const list = membersOf(slice).annotationSets;
  list.splice(0, list.length, ...annotationSets);
};

/** @param {AnnotationSetsSlice} slice */
const reset = (slice) => install(slice, { annotationSets: [] });

/**
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSet} set
 */
export const addAnnotationSet = (slice, set) => {
  membersOf(slice).annotationSets.push(set);
};

/**
 * @param {AnnotationSetsSlice} slice
 * @param {AnnotationSet} set A set of the list; another value changes nothing.
 */
export const removeAnnotationSet = (slice, set) => {
  const list = membersOf(slice).annotationSets;
  const index = list.indexOf(set);
  if (index >= 0) list.splice(index, 1);
};

export const ANNOTATION_SETS_SLICE = Object.freeze({ members: Object.freeze(['annotationSets']), create, capture, install, reset });
