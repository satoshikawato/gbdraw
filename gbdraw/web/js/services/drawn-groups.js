// @ts-check
// The groups of a Result that the Web app finds by what they draw: the scale,
// the dinucleotide and depth tracks, the plot title. This module is the one
// owner of their slot renderer names and ids, which Python writes
// (gbdraw/svg/ids.py, gbdraw/features/ids.py, and the track renderers).

/**
 * The groups that draw one kind of track. A Result marks each with its slot
 * renderer (`data-gbdraw-slot-renderer`), whatever its id (a custom track slot
 * is `track_slot_<hash>_<slot>`). A Result stored without the marks is read
 * by id: `<id>`, or `<id>_…` for a record of a grid, a Linear record, or a
 * child group.
 * @typedef {{ readonly renderers: readonly string[], readonly ids: readonly string[] }} TrackGroups
 */

/** @param {readonly string[]} renderers @param {readonly string[]} ids @returns {TrackGroups} */
const trackGroups = (renderers, ids) => Object.freeze({ renderers, ids });

export const TICK_GROUPS = trackGroups(['ticks'], ['tick']);
export const DINUCLEOTIDE_CONTENT_GROUPS = trackGroups(['dinucleotide_content'], ['gc_content']);
// Circular draws GC skew as `skew`, Linear as `gc_skew`.
export const DINUCLEOTIDE_SKEW_GROUPS = trackGroups(['dinucleotide_skew'], ['skew', 'gc_skew']);
export const DEPTH_GROUPS = trackGroups(['depth'], ['depth']);

// The Linear scale (bar or ruler style) and the plot title are one group each.
export const LENGTH_BAR_ID = 'length_bar';
export const PLOT_TITLE_ID = 'plot_title';

/** @param {string} renderer */
const rendererSelector = (renderer) => `g[data-gbdraw-slot-renderer="${renderer}"]`;
/** @param {string} id */
const idSelector = (id) => `g[id="${id}"], g[id^="${id}_"]`;

// Every group of a kind, marked or not, as one selector (for `closest` and for
// a list of groups to leave out).
/** @param {TrackGroups} groups */
export const trackGroupSelector = ({ renderers, ids }) => [...renderers.map(rendererSelector), ...ids.map(idSelector)].join(', ');

/**
 * The groups of a kind in a Result: the marked ones, or without any, those
 * found by id.
 * @param {ParentNode | null | undefined} svg
 * @param {TrackGroups} groups
 * @returns {Element[]}
 */
export const findTrackGroups = (svg, { renderers, ids }) => {
  if (!svg) return [];
  const marked = renderers.flatMap((renderer) => Array.from(svg.querySelectorAll(rendererSelector(renderer))));
  return [...new Set(marked.length > 0 ? marked : ids.flatMap((id) => Array.from(svg.querySelectorAll(idSelector(id)))))];
};
