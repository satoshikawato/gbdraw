// @ts-check
// The feature-edits slice of a drawing: its per-feature edits
// (`featureOverrides`, Session 46 `modes.<m>.features.featureOverrides`) and
// its Feature placement rows (`featurePlacementOverrides`, `config`). Both are
// draft maps keyed by original-source feature identity (the draft key
// `featureIdentityKeyOf` builds, services/feature-placement.js, which owns the
// row shape and the request projection). The drawing shows them read-only to
// the type checker; only this module writes them, through
// `drawing.slices.featureEdits`: a row set or removed by its draft key (the
// feature editor's writers), and Load, a History restore, and Reset
// (`install`, `reset`). A caller passes rows it checked.
import { reactive } from '../../utils/reactivity.js';
import { cloneJsonData } from '../json-clone.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';
/** @import { DeepReadonly } from '../../utils/reactivity.js' */
/** @import { SliceHandle } from './slice-handle.js' */
/** @import { FeatureOverrideDraft, FeaturePlacementDraft } from '../feature-placement.js' */
/** @import { FeatureOverrideDraftRow } from '../../contracts/session.generated.js' */
/** @import { FeaturePlacementDraftRow } from '../feature-placement.js' */

/** @typedef {{ featureOverrides: FeatureOverrideDraft, featurePlacementOverrides: FeaturePlacementDraft }} FeatureEditsMembers */
/** @typedef {FeatureEditsMembers} FeatureEditsData What `capture` returns. */
/**
 * @typedef {{
 *   readonly featureOverrides: DeepReadonly<FeatureOverrideDraft>,
 *   readonly featurePlacementOverrides: DeepReadonly<FeaturePlacementDraft>
 * }} FeatureEditsView
 */
/** @typedef {SliceHandle<'featureEdits'>} FeatureEditsSlice The handle a writer takes. */

// A feature-edits handle exists only through `create`, which gives it these members.
/** @param {FeatureEditsSlice} slice */
const membersOf = (slice) => /** @type {FeatureEditsMembers} */ (sliceMembers(slice, 'featureEdits'));

/**
 * @template T
 * @param {Record<string, T>} target
 * @param {Readonly<Record<string, T>>} entries Copied first, so `target` may be `entries`.
 */
const replaceEntries = (target, entries) => {
  const copy = cloneJsonData({ ...entries });
  Object.keys(target).forEach((key) => { delete target[key]; });
  Object.assign(target, copy);
};

/**
 * A new slice, or one over a drawing's existing members (a test fixture).
 * @param {Partial<FeatureEditsMembers>} [members]
 * @returns {{ handle: FeatureEditsSlice, view: FeatureEditsView }}
 */
const create = (members = {}) => {
  /** @type {FeatureEditsMembers} */
  const data = {
    featureOverrides: members.featureOverrides ?? reactive(/** @type {FeatureOverrideDraft} */ ({})),
    featurePlacementOverrides: members.featurePlacementOverrides ?? reactive(/** @type {FeaturePlacementDraft} */ ({}))
  };
  return { handle: createSliceHandle('featureEdits', data), view: Object.freeze({ ...data }) };
};

/**
 * The drafts, read-only, for a writer's caller that holds only the handle.
 * @param {FeatureEditsSlice} slice
 * @returns {FeatureEditsView}
 */
export const featureEditsOf = (slice) => membersOf(slice);

/**
 * @param {FeatureEditsSlice} slice
 * @returns {FeatureEditsData}
 */
const capture = (slice) => {
  const { featureOverrides, featurePlacementOverrides } = membersOf(slice);
  return { featureOverrides: cloneJsonData({ ...featureOverrides }), featurePlacementOverrides: cloneJsonData({ ...featurePlacementOverrides }) };
};

// Each member `data` holds becomes a copy of it; the others are kept (a
// Generate's checkpoint restores the placements on their own).
/**
 * @param {FeatureEditsSlice} slice
 * @param {Partial<DeepReadonly<FeatureEditsData>>} data
 */
const install = (slice, { featureOverrides, featurePlacementOverrides }) => {
  const members = membersOf(slice);
  // A copy of the rows of a draft map is a draft map (the casts: the copy is writable).
  if (featureOverrides) replaceEntries(members.featureOverrides, /** @type {FeatureOverrideDraft} */ (featureOverrides));
  if (featurePlacementOverrides) {
    replaceEntries(members.featurePlacementOverrides, /** @type {FeaturePlacementDraft} */ (featurePlacementOverrides));
  }
};

/** @param {FeatureEditsSlice} slice */
const reset = (slice) => install(slice, { featureOverrides: {}, featurePlacementOverrides: {} });

/**
 * One feature's per-feature edit row, by its draft key; null removes it.
 * @param {FeatureEditsSlice} slice
 * @param {string} key The row's draft key (`featureIdentityKeyOf(row)`).
 * @param {Readonly<FeatureOverrideDraftRow> | null} row Copied.
 */
export const setFeatureOverrideRow = (slice, key, row) => {
  const { featureOverrides } = membersOf(slice);
  if (row) featureOverrides[key] = { ...row };
  else delete featureOverrides[key];
};

/**
 * Clears the per-feature edits (a new Session).
 * @param {FeatureEditsSlice} slice
 */
export const clearFeatureOverrides = (slice) => install(slice, { featureOverrides: {} });

/**
 * One feature's placement row, by its draft key; null removes it (Auto).
 * @param {FeatureEditsSlice} slice
 * @param {string} key The row's draft key (`featureIdentityKeyOf(row)`).
 * @param {DeepReadonly<FeaturePlacementDraftRow> | null} row Copied.
 */
export const setFeaturePlacementRow = (slice, key, row) => {
  const { featurePlacementOverrides } = membersOf(slice);
  // A copy of a placement row is a placement row (the cast: the copy is writable).
  if (row) featurePlacementOverrides[key] = /** @type {FeaturePlacementDraftRow} */ (cloneJsonData(row));
  else delete featurePlacementOverrides[key];
};

export const FEATURE_EDITS_SLICE = Object.freeze({
  members: Object.freeze(['featureOverrides', 'featurePlacementOverrides']), create, capture, install, reset
});
