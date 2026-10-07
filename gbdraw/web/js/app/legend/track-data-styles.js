// @ts-check
import { getDepthTrackFallbackLabel } from '../depth-tracks.js';

const text = (value) => String(value ?? '').trim();

/**
 * The Legend captions that a draft's track data names: the legend labels of
 * annotation sets and annotations, and the rows of Depth series that have a
 * source (the label of the series, or of a Depth row that names it). Python
 * draws a row from such data only, so a Legend style on a caption the data no
 * longer names can no longer be admitted (OV-65).
 * @param {{
 *   annotationSets?: any[],
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
  const add = (value) => {
    const caption = text(value);
    if (caption) captions.add(caption);
  };
  (Array.isArray(annotationSets) ? annotationSets : []).forEach((set) => {
    add(set?.legendLabel);
    (Array.isArray(set?.annotations) ? set.annotations : []).forEach((annotation) => add(annotation?.legendLabel));
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
 * The one rule for Legend styles of track data: styles follow the caption.
 * Returns the port that runs a data change and then retires the Legend color
 * and stroke of every caption the track data named before the change and no
 * longer names. A caption the new data still names keeps its style. The
 * change and the retirement happen in the History step of the caller, so Undo
 * restores both.
 * @param {{
 *   legendColorOverrides: Record<string, any>,
 *   legendStrokeOverrides: Record<string, any>,
 *   namedCaptions: () => Set<string>
 * }} options
 * @returns {<T>(change: () => T) => T}
 */
export const buildLegendStyleRetirement = ({ legendColorOverrides, legendStrokeOverrides, namedCaptions }) => (
  (change) => {
    const before = namedCaptions();
    const result = change();
    const after = namedCaptions();
    before.forEach((caption) => {
      if (after.has(caption)) return;
      delete legendColorOverrides[caption];
      delete legendStrokeOverrides[caption];
    });
    return result;
  }
);
