// @ts-check
import { getDepthTrackFallbackLabel } from '../depth-tracks.js';
import { retireLegendPaint } from '../../services/drawing/paint.js';
/** @import { PaintSlice } from '../../services/drawing/paint.js' */

/** @param {unknown} value */
const text = (value) => String(value ?? '').trim();

/**
 * The Legend captions that a draft's track data names: the legend labels of
 * annotation sets and annotations, and the rows of Depth series that have a
 * source (the label of the series, or of a Depth row that names it). Python
 * draws a row from such data only, so a Legend style on a caption the data no
 * longer names can no longer be admitted (OV-65).
 * @param {{
 *   annotationSets?: readonly { legendLabel?: unknown, annotations?: unknown }[],
 *   depthTracks?: any[],
 *   depthSlots?: any[],
 *   sourcedDepthTrackIndexes?: number[]
 * }} [data]
 * @returns {Set<string>}
 */
export const trackDataLegendCaptions = ({
  annotationSets = [],
  depthTracks = [],
  depthSlots = [],
  sourcedDepthTrackIndexes = []
} = {}) => {
  const captions = new Set();
  /** @param {unknown} value */
  const add = (value) => {
    const caption = text(value);
    if (caption) captions.add(caption);
  };
  (Array.isArray(annotationSets) ? annotationSets : []).forEach((set) => {
    add(set?.legendLabel);
    (Array.isArray(set?.annotations) ? set.annotations : []).forEach((/** @type {{ legendLabel?: unknown } | null | undefined} */ annotation) => add(annotation?.legendLabel));
  });
  const sourced = new Set(sourcedDepthTrackIndexes);
  sourced.forEach((index) => add(depthTracks?.[index]?.label || getDepthTrackFallbackLabel(index)));
  (Array.isArray(depthSlots) ? depthSlots : []).forEach((slot) => {
    if (slot?.renderer === 'depth' && sourced.has(Number(slot?.params?.track_index ?? 0))) {
      add(slot.params?.legend_label);
    }
  });
  return captions;
};

/**
 * The one rule for Legend edits of track data: styles and names follow the
 * caption. Returns the port that runs a data change and then, for every
 * caption the track data named before the change and no longer names, retires
 * the Legend color and stroke stored under it (the row's key) and the Legend
 * rename of its row (OV-87). A caption
 * the new data still names keeps its rename and styles; a deleted row stays
 * deleted. The change and the retirement happen in the History step of the
 * caller, so Undo restores them together. `projectLegendEntries` shows the
 * entries with a retired rename on the displayed Result, as Undo and Redo do.
 * The renamed rows a Generate hid (OV-120) follow the same rule: a row whose
 * data is removed leaves `dormantLegendEntries` with its styles.
 * @param {{
 *   drawing: { readonly slices: { readonly paint: PaintSlice } },
 *   namedCaptions: () => Set<string>,
 *   projectLegendEntries: () => void
 * }} options
 * @returns {<T>(change: () => T) => T}
 */
export const buildLegendStyleRetirement = ({ drawing, namedCaptions, projectLegendEntries }) => (
  (change) => {
    const before = namedCaptions();
    const result = change();
    const after = namedCaptions();
    const retired = new Set([...before].filter((caption) => !after.has(caption)));
    if (retired.size > 0 && retireLegendPaint(drawing.slices.paint, retired)) projectLegendEntries();
    return result;
  }
);
