// @ts-check
/** @import { DrawingState } from '../state.js' */
import { isLegendOrderEdited } from '../services/legend-svg.js';

// Which editor state each Result's bytes show (design §5.4). A Result
// displayed again receives a domain only when that state changed since it
// was last shown, so a Result without new edits gets no projection work and
// an Undo reaches a Result that is displayed again. The displayed Result
// follows every live edit; app-setup records the states, and a label
// rerender names the state it drew from (`editorPaintState`).

/**
 * @typedef {{
 *   colors: [unknown, string],
 *   visibility: string,
 *   strokes: string,
 *   labels: string,
 *   legendOrder: string
 * }} EditorPaintState
 */

const toRaw = globalThis.window?.Vue?.toRaw ?? (/** @type {<T>(value: T) => T} */ (value) => value);

/**
 * What the paint state reads of the application state besides the drawing.
 * @typedef {{ appliedPaletteColors: { value: unknown }, originalLegendOrder: { value: string[] } }} PaintStateSource
 */

/**
 * The paint state of a drawing's editor intent.
 * @param {PaintStateSource} state
 * @param {DrawingState} drawing
 * @returns {EditorPaintState}
 */
export const editorPaintState = (state, drawing) => ({
  colors: [
    toRaw(state.appliedPaletteColors.value),
    JSON.stringify([drawing.manualSpecificRules, drawing.featureColorOverrides, drawing.legendColorOverrides])
  ],
  visibility: JSON.stringify([
    Object.values(drawing.featureOverrides).map((row) => [row.recordKey, row.biologicalFeatureId, row.featureVisibility]),
    drawing.featureVisibilityManualRules
  ]),
  // A deleted row strokes no feature (OV-293).
  strokes: JSON.stringify([
    drawing.featureStrokeOverrides, drawing.legendStrokeOverrides,
    drawing.deletedLegendEntries.value.map((entry) => entry.originalCaption || entry.caption)
  ]),
  labels: JSON.stringify([
    Object.values(drawing.featureOverrides).map((row) => [
      row.recordKey, row.biologicalFeatureId, row.labelVisibility, row.labelText, row.labelSourceText
    ]),
    drawing.labelTextBulkOverrides
  ]),
  // An edited Legend order, or '' for the default order (D-08).
  legendOrder: isLegendOrderEdited(drawing.legendEntries.value, state.originalLegendOrder.value)
    ? JSON.stringify(drawing.legendEntries.value.map((entry) => entry.caption))
    : ''
});

// A paint state no editor state equals: a Result recorded with it shows
// every paint domain and the labels on its next display. The Legend order is
// kept: replaying the default order needs the Result's generated inventory.
/** @param {EditorPaintState} current @returns {EditorPaintState} */
const unknownPaint = (current) => ({ ...current, colors: [null, ''], visibility: '', strokes: '', labels: '' });

export const createResultPaintRecord = () => {
  /** @type {Map<string, EditorPaintState>} */
  const shownByResult = new Map();
  let bound = '';
  // The paint domains the displayed Result was not shown (a declined or
  // failed display); recorded with it when it leaves.
  /** @type {Partial<EditorPaintState> | null} */
  let behind = null;
  /** @param {EditorPaintState} state */
  const recordBound = (state) => {
    if (bound && shownByResult.has(bound)) shownByResult.set(bound, { ...state, ...behind });
    behind = null;
  };
  return {
    /**
     * A committed artifact: each Result first seen shows `drawn` (`current`
     * unless given), and a Result no longer live is forgotten. `displayed` is
     * the Result mounted by the commit. A `restored` artifact (Session Load,
     * a History restore) keeps each Result's bytes as they were saved or
     * kept: a batch Result not displayed since the last edits lacks them, so
     * every Result but the mounted one shows every paint domain on its first
     * display. A label rerender draws the fills and visibility of the state
     * it started from (`drawn`, OV-346): the displayed Result lacks those of
     * an edit made meanwhile (`lacking`) until they are shown.
     * @param {readonly string[]} identities
     * @param {string} displayed
     * @param {EditorPaintState} current
     * @param {{ restored?: boolean, drawn?: EditorPaintState }} [options]
     */
    commit(identities, displayed, current, { restored = false, drawn = current } = {}) {
      const live = new Set(identities.filter(Boolean));
      live.forEach((identity) => {
        if (shownByResult.has(identity)) return;
        shownByResult.set(identity, restored && identity !== displayed ? unknownPaint(current) : drawn);
      });
      [...shownByResult.keys()].forEach((identity) => {
        if (!live.has(identity)) shownByResult.delete(identity);
      });
      bound = displayed;
      const { colors, visibility } = drawn;
      behind = colors[0] !== current.colors[0] || colors[1] !== current.colors[1] || visibility !== current.visibility
        ? { colors, visibility } : null;
    },
    /**
     * The displayed Result leaves (a mode switch) having followed every live
     * edit. Returns its identity, or '' when none is recorded.
     * @param {EditorPaintState} current
     */
    depart(current) {
      const departed = bound && shownByResult.has(bound) ? bound : '';
      recordBound(current);
      bound = '';
      return departed;
    },
    /**
     * A Result is displayed: the Result shown until now followed every live
     * edit. Returns the state the arriving Result's bytes show.
     * @param {string} identity
     * @param {EditorPaintState} current
     * @returns {EditorPaintState}
     */
    display(identity, current) {
      if (bound !== identity) recordBound(current);
      bound = identity;
      return shownByResult.get(identity) || current;
    },
    /**
     * The displayed Result now shows `current`, but the fills and visibility
     * of `previous` when its palette and rules projection `declined`, and
     * every paint domain of `previous` when its projection or compile
     * `failed` (its labels follow in the binder's label step).
     * @param {string} identity
     * @param {EditorPaintState} current
     * @param {EditorPaintState} previous
     * @param {{ declined?: boolean, failed?: boolean }} [options]
     */
    shown(identity, current, previous, { declined = false, failed = false } = {}) {
      const { colors, visibility, strokes, legendOrder } = previous;
      behind = failed ? { colors, visibility, strokes, legendOrder } : (declined ? { colors, visibility } : null);
      shownByResult.set(identity, { ...current, ...behind });
    },
    /**
     * What the displayed Result `identity` still shows of the state before a
     * declined or failed display, or null (shown in full, or not displayed).
     * @param {string} identity
     * @returns {Partial<EditorPaintState> | null}
     */
    lacking(identity) {
      return bound === identity ? behind : null;
    },
    /** @param {string} identity The Result shown again with the intent it left with (E1). */
    rebind(identity) {
      bound = identity;
      behind = null;
    }
  };
};
