// @ts-check
import { getDepthTrackFallbackLabel } from '../depth-tracks.js';

/** @param {unknown} value */
const text = (value) => String(value ?? '').trim();

/** @typedef {{ caption?: unknown, originalCaption?: unknown }} KeyedLegendRow A Legend list row: its name and its key. */

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
 * The rows of a drawn track the draft hides (`hiddenRows`: Show Depth off,
 * Hide GC Content, Hide GC Skew) leave the list with the track, and a renamed
 * one waits in `dormantLegendEntries` as after a Generate (OV-1198, OV-1215);
 * once the track shows again, they return. The Legend list owner lists them
 * (`relistLegendEntries`, which leaves out the hidden rows), and a returned
 * row leaves `dormantLegendEntries`.
 * @param {{
 *   legendColorOverrides: Record<string, any>,
 *   legendStrokeOverrides: Record<string, any>,
 *   legendEntries: { value: any[] },
 *   dormantLegendEntries?: { value: any[] },
 *   namedCaptions: () => Set<string>,
 *   hiddenRows?: () => readonly string[],
 *   relistLegendEntries?: () => void,
 *   projectLegendEntries: () => void
 * }} options
 * @returns {<T>(change: () => T) => T}
 */
export const buildLegendStyleRetirement = ({
  legendColorOverrides, legendStrokeOverrides, legendEntries, dormantLegendEntries, namedCaptions,
  hiddenRows = () => [], relistLegendEntries = () => {}, projectLegendEntries
}) => (
  (change) => {
    const before = namedCaptions();
    const hiddenBefore = new Set(hiddenRows());
    const result = change();
    const after = namedCaptions();
    /** @param {string} caption */
    const retireStyles = (caption) => {
      delete legendColorOverrides[caption];
      delete legendStrokeOverrides[caption];
    };
    const retired = new Set([...before].filter((caption) => !after.has(caption)));
    retired.forEach(retireStyles);
    let listChanged = false;
    /** @param {KeyedLegendRow | null | undefined} entry */
    const keyOf = (entry) => text(entry?.originalCaption || entry?.caption);
    const entries = (Array.isArray(legendEntries.value) ? legendEntries.value : []).map((entry) => {
      const original = text(entry?.originalCaption);
      const caption = text(entry?.caption);
      if (!retired.has(original) || !caption || caption === original) return entry;
      listChanged = true;
      return { ...entry, caption: original };
    });
    const waiting = dormantLegendEntries && Array.isArray(dormantLegendEntries.value) ? dormantLegendEntries.value : null;
    let dormant = waiting ? waiting.filter((entry) => !retired.has(text(entry?.originalCaption))) : null;
    const hiddenAfter = new Set(hiddenRows());
    const hiding = new Set([...hiddenAfter].filter((key) => !hiddenBefore.has(key)));
    const showing = new Set([...hiddenBefore].filter((key) => !hiddenAfter.has(key)));
    const trackShownOrHidden = hiding.size > 0 || showing.size > 0;
    if (trackShownOrHidden) {
      const renamed = entries.filter((entry) => hiding.has(keyOf(entry)) && text(entry?.caption) !== keyOf(entry));
      if (dormant && renamed.length > 0) dormant = [...dormant.filter((entry) => !hiding.has(keyOf(entry))), ...renamed];
      // The list owner names a returning row from the list, else from the
      // waiting rename; once listed, the row no longer waits (review H1).
      if (dormantLegendEntries && dormant) dormantLegendEntries.value = dormant;
      legendEntries.value = entries;
      relistLegendEntries();
      const listed = new Set((Array.isArray(legendEntries.value) ? legendEntries.value : []).map(keyOf));
      dormant = dormant && dormant.filter((entry) => !(showing.has(keyOf(entry)) && listed.has(keyOf(entry))));
    }
    if (dormantLegendEntries && dormant && waiting && (dormant.length !== waiting.length
      || dormant.some((entry, index) => entry !== waiting[index]))) {
      dormantLegendEntries.value = dormant;
    }
    if (trackShownOrHidden) {
      projectLegendEntries();
      return result;
    }
    if (listChanged) {
      legendEntries.value = entries;
      projectLegendEntries();
    }
    return result;
  }
);
