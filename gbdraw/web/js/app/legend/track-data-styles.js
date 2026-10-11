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
 * The one rule for Legend edits of track data: styles and names follow the
 * caption. Returns the port that runs a data change and then, for every
 * caption the track data named before the change and no longer names, retires
 * the Legend color and stroke stored under it, and the Legend rename of its
 * row with the color and stroke stored under the new name (OV-87). A caption
 * the new data still names keeps its rename and styles; a deleted row stays
 * deleted. The change and the retirement happen in the History step of the
 * caller, so Undo restores them together. `projectLegendEntries` shows the
 * entries with a retired rename on the displayed Result, as Undo and Redo do.
 * The renamed rows a Generate hid (OV-120) follow the same rule: a row whose
 * data is removed leaves `dormantLegendEntries` with its styles.
 * @param {{
 *   legendColorOverrides: Record<string, any>,
 *   legendStrokeOverrides: Record<string, any>,
 *   legendEntries: { value: any[] },
 *   dormantLegendEntries?: { value: any[] },
 *   namedCaptions: () => Set<string>,
 *   projectLegendEntries: () => void
 * }} options
 * @returns {<T>(change: () => T) => T}
 */
export const buildLegendStyleRetirement = ({
  legendColorOverrides, legendStrokeOverrides, legendEntries, dormantLegendEntries, namedCaptions, projectLegendEntries
}) => (
  (change) => {
    const before = namedCaptions();
    const result = change();
    const after = namedCaptions();
    /** @param {string} caption */
    const retireStyles = (caption) => {
      delete legendColorOverrides[caption];
      delete legendStrokeOverrides[caption];
    };
    const retired = new Set([...before].filter((caption) => !after.has(caption)));
    retired.forEach(retireStyles);
    let renameRetired = false;
    const entries = (Array.isArray(legendEntries.value) ? legendEntries.value : []).map((entry) => {
      const original = text(entry?.originalCaption);
      const caption = text(entry?.caption);
      if (!retired.has(original) || !caption || caption === original) return entry;
      if (!after.has(caption)) retireStyles(caption);
      renameRetired = true;
      return { ...entry, caption: original };
    });
    if (dormantLegendEntries && Array.isArray(dormantLegendEntries.value)) {
      const dormant = dormantLegendEntries.value.filter((entry) => {
        const original = text(entry?.originalCaption);
        if (!retired.has(original)) return true;
        const caption = text(entry?.caption);
        if (caption && !after.has(caption)) retireStyles(caption);
        return false;
      });
      if (dormant.length !== dormantLegendEntries.value.length) dormantLegendEntries.value = dormant;
    }
    if (renameRetired) {
      legendEntries.value = entries;
      projectLegendEntries();
    }
    return result;
  }
);
