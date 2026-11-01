// @ts-check
/** @import { DrawingState } from '../state.js' */
/** @import { KeptLegendRow } from '../services/drawing/paint.js' */
/** @import { SpecificColorRule } from '../services/specific-color-rules.js' */
import { isLegendOrderEdited } from '../services/legend-svg.js';
import { legendEntryKey } from '../services/drawing/paint.js';
import { toRaw } from '../utils/vue-raw.js';
import { featureVisibilityRows } from '../services/feature-visibility.js';

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

/**
 * What the paint state reads of the application state besides the drawing.
 * @typedef {{ appliedPaletteColors: { value: unknown }, originalLegendOrder: { value: string[] } }} PaintStateSource
 */

/**
 * The Feature visibility a drawing's editor intent draws: its per-feature rows and its rules.
 * @param {DrawingState} drawing
 * @returns {string}
 */
export const visibilityPaintState = (drawing) => JSON.stringify([
  featureVisibilityRows(Object.values(drawing.featureOverrides)),
  drawing.featureVisibilityManualRules
]);

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
  visibility: visibilityPaintState(drawing),
  // A deleted row strokes no feature (OV-293).
  strokes: JSON.stringify([
    drawing.featureStrokeOverrides, drawing.legendStrokeOverrides,
    drawing.deletedLegendEntries.value.map((/** @type {KeptLegendRow} */ entry) => entry.originalCaption || entry.caption)
  ]),
  labels: JSON.stringify([
    Object.values(drawing.featureOverrides).map((row) => [
      row.recordKey, row.biologicalFeatureId, row.labelVisibility, row.labelText, row.labelSourceText
    ]),
    drawing.labelTextBulkOverrides
  ]),
  // An edited Legend order of row keys, or '' for the default order (D-08).
  legendOrder: isLegendOrderEdited(drawing.legendEntries.value, state.originalLegendOrder.value)
    ? JSON.stringify(drawing.legendEntries.value.map(legendEntryKey))
    : ''
});

// The rule table a paint state shows, or null when it is unknown (a Result
// restored without the state it shows, `unknownPaint`).
/** @param {EditorPaintState} paint @returns {SpecificColorRule[] | null} */
export const paintRuleTable = (paint) => (paint.colors[1]
  ? /** @type {[SpecificColorRule[], unknown, unknown]} */ (JSON.parse(paint.colors[1]))[0]
  : null);

// Whether two paint states show the same palette and color edits.
/** @param {EditorPaintState['colors']} left @param {EditorPaintState['colors']} right */
export const sameColors = (left, right) => left[0] === right[0] && left[1] === right[1];

// Whether two paint states show the same fills: palette and color edits, and visibility.
/** @param {Pick<EditorPaintState, 'colors' | 'visibility'>} left @param {Pick<EditorPaintState, 'colors' | 'visibility'>} right */
const sameFills = (left, right) => sameColors(left.colors, right.colors) && left.visibility === right.visibility;

// A paint state no editor state equals: a Result recorded with it shows
// every paint domain and the labels on its next display. The Legend order is
// kept: replaying the default order needs the Result's generated inventory.
/** @param {EditorPaintState} current @returns {EditorPaintState} */
const unknownPaint = (current) => ({ ...current, colors: [null, ''], visibility: '', strokes: '', labels: '' });

export const createResultPaintRecord = () => {
  /** @type {Map<string, EditorPaintState>} */
  const shownByResult = new Map();
  // The state of its drawing at the restore, for a restored Result not displayed:
  // what its bytes are taken to show (OV-2326).
  /** @type {Map<string, EditorPaintState>} */
  const stateAtRestore = new Map();
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
     * an edit made meanwhile (`lacking`) until they are shown. `drawingStateOf`
     * reads the state of the drawing that holds a restored Result not displayed.
     * @param {readonly string[]} identities
     * @param {string} displayed
     * @param {EditorPaintState} current
     * @param {{ restored?: boolean, drawn?: EditorPaintState, drawingStateOf?: (identity: string) => EditorPaintState }} [options]
     */
    commit(identities, displayed, current, { restored = false, drawn = current, drawingStateOf = () => current } = {}) {
      const live = new Set(identities.filter(Boolean));
      live.forEach((identity) => {
        if (shownByResult.has(identity)) return;
        const notDisplayed = restored && identity !== displayed;
        shownByResult.set(identity, notDisplayed ? unknownPaint(current) : drawn);
        if (notDisplayed) stateAtRestore.set(identity, drawingStateOf(identity));
      });
      [...shownByResult.keys()].forEach((identity) => {
        if (live.has(identity)) return;
        shownByResult.delete(identity);
        stateAtRestore.delete(identity);
      });
      bound = displayed;
      const { colors, visibility } = drawn;
      behind = sameFills(drawn, current) ? null : { colors, visibility };
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
    /**
     * Whether a label rerender drawn from `drawn` shows the `current` editor
     * state as Python draws it: the same fills and strokes (D-36). The labels
     * and the Legend order are projected on any figure. Read when its Result
     * is ready, so an edit made after its paint was committed counts too.
     * @param {EditorPaintState | null | undefined} drawn
     * @param {EditorPaintState} current
     */
    drewPaint(drawn, current) {
      return drawn ? sameFills(drawn, current) && drawn.strokes === current.strokes : false;
    },
    /**
     * The rule table a Result was restored with, when it is not known to show another.
     * @param {string} identity
     * @returns {SpecificColorRule[] | null}
     */
    restoredRules(identity) {
      const restored = stateAtRestore.get(identity);
      return restored ? paintRuleTable(restored) : null;
    },
    /**
     * Whether the label or Feature visibility intent changed since the Result
     * was last shown (`previous`, from `display`), or restored when it has not
     * been shown since: only a rerender places its labels for what its drawing
     * holds now. The display records what it shows, so a Result answers once
     * for a change.
     * @param {string} identity
     * @param {EditorPaintState} previous
     * @param {EditorPaintState} current
     */
    placedForOtherIntent(identity, previous, current) {
      const restored = stateAtRestore.get(identity) ?? current;
      return (previous.labels || restored.labels) !== current.labels
        || (previous.visibility || restored.visibility) !== current.visibility;
    },
    /** @param {string} identity The Result shown again with the intent it left with (E1). */
    rebind(identity) {
      bound = identity;
      behind = null;
    }
  };
};
