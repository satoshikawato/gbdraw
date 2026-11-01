// @ts-check
import { defaultDepthSeriesLabel } from '../../services/depth-track-state.js';
import {
  legendEntryKey, paintOf, replaceDormantLegendRows, retireLegendPaint
} from '../../services/drawing/paint.js';
/** @import { PaintSlice } from '../../services/drawing/paint.js' */

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
  // A series without a label draws the name its index has among all the series.
  const seriesCount = Math.max(Array.isArray(depthTracks) ? depthTracks.length : 0, ...[...sourced].map((index) => index + 1));
  sourced.forEach((index) => add(depthTracks?.[index]?.label || defaultDepthSeriesLabel(index, seriesCount)));
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
 *   drawing: { readonly slices: { readonly paint: PaintSlice } },
 *   namedCaptions: () => Set<string>,
 *   hiddenRows?: () => readonly string[],
 *   relistLegendEntries?: () => void,
 *   projectLegendEntries: () => void
 * }} options
 * @returns {<T>(change: () => T) => T}
 */
export const buildLegendStyleRetirement = ({ drawing, namedCaptions,
  hiddenRows = () => [], relistLegendEntries = () => {}, projectLegendEntries }) => (
  (change) => {
    const before = namedCaptions();
    const hiddenBefore = new Set(hiddenRows());
    const result = change();
    const after = namedCaptions();
    const retired = new Set([...before].filter((caption) => !after.has(caption)));
    const paint = paintOf(drawing.slices.paint);
    const renameRetired = retired.size > 0 && retireLegendPaint(drawing.slices.paint, retired);
    const hiddenAfter = new Set(hiddenRows());
    const hiding = new Set([...hiddenAfter].filter((key) => !hiddenBefore.has(key)));
    const showing = new Set([...hiddenBefore].filter((key) => !hiddenAfter.has(key)));
    if (hiding.size === 0 && showing.size === 0) {
      if (renameRetired) projectLegendEntries();
      return result;
    }
    // A renamed row of a track the change hides waits with its rename. The list
    // owner names a returning row from the list, else from the waiting rename;
    // once listed, the row no longer waits (review H1).
    const renamed = paint.legendEntries.value
      .filter((entry) => hiding.has(legendEntryKey(entry)) && String(entry.caption ?? '').trim() !== legendEntryKey(entry));
    if (renamed.length > 0) {
      replaceDormantLegendRows(drawing.slices.paint, [
        ...paint.dormantLegendEntries.value.filter((entry) => !hiding.has(legendEntryKey(entry))), ...renamed
      ]);
    }
    relistLegendEntries();
    const listed = new Set(paint.legendEntries.value.map(legendEntryKey));
    const waiting = paint.dormantLegendEntries.value;
    const dormant = waiting.filter((entry) => !(showing.has(legendEntryKey(entry)) && listed.has(legendEntryKey(entry))));
    if (dormant.length !== waiting.length) replaceDormantLegendRows(drawing.slices.paint, dormant);
    projectLegendEntries();
    return result;
  }
);

/**
 * @typedef {{ readonly id: string, readonly slices: { readonly paint: PaintSlice } }} RetirementDrawing
 *   A drawing the Legend port changes: found by its ID, styled by its paint.
 */

/**
 * The Legend port of every drawing: `retire(drawing, change)` runs a track data
 * change of that drawing and retires the Legend styles of the captions its
 * data stops naming: `namedCaptions(drawing)` reads that drawing's form, and
 * `hiddenRows(drawing)` reads that drawing's form flags (Show Depth and the
 * like) over the rows displayed in the shown Result, so for a drawing not
 * shown only the flags are its own. A port write to a drawing not shown never
 * changes them (D-15-183), so its hidden set stays unchanged. The relist and
 * the projection read the displayed Result, so only the shown drawing, found
 * by its ID, runs them; a replacement
 * that reaches a drawing not shown retires its paint the same way (review S2).
 * @template {RetirementDrawing} D
 * @param {{
 *   shownDrawingId: () => string,
 *   namedCaptions: (drawing: D) => Set<string>,
 *   hiddenRows: (drawing: D) => readonly string[],
 *   relistLegendEntries: () => void,
 *   projectLegendEntries: () => void
 * }} options
 * @returns {<T>(drawing: D, change: () => T) => T}
 */
export const createDrawingLegendRetirement = ({
  shownDrawingId, namedCaptions, hiddenRows, relistLegendEntries, projectLegendEntries
}) => (drawing, change) => {
  const shown = drawing.id === shownDrawingId();
  return buildLegendStyleRetirement({
    drawing,
    namedCaptions: () => namedCaptions(drawing),
    hiddenRows: () => hiddenRows(drawing),
    relistLegendEntries: shown ? relistLegendEntries : () => {},
    projectLegendEntries: shown ? projectLegendEntries : () => {}
  })(change);
};
