// @ts-check
/** @import { DrawingState } from '../state.js' */
/** @import { KeptLegendRow } from '../services/drawing/paint.js' */
/** @import { SpecificColorRule } from '../services/specific-color-rules.js' */
/** @import { SvgResultEntry } from '../services/svg-result-ingestion.js' */
import { isLegendOrderEdited } from '../services/legend-svg.js';
import { legendEntryKey } from '../services/drawing/paint.js';
import { toRaw } from '../utils/vue-raw.js';
import { featureVisibilityRows } from '../services/feature-visibility.js';
import { featureIdentityKey } from '../services/feature-placement.js';
import { paletteColorsEqual } from './svg-styles.js';

// Which editor state each Result's bytes were drawn with and show: one
// record per Result. A Result displayed again receives a
// domain only when that state changed since it was last shown, so a Result
// without new edits gets no projection work and an Undo reaches a Result that
// is displayed again. The displayed Result follows every live edit;
// app/editor-intent.js records the states, and a label rerender names the
// state it drew from (`editorPaintState`).

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

/**
 * The paint a Generate request draws, read as the request is built (its one producer), so that an edit made while it runs
 * is an edit since: `palette` is the one its commit applies, `bound` the drawing's paint with it as the request
 * starts, and `of` reads a drawing whose rows the request moved (OV-1949) with the same palette.
 * @param {DrawingState} drawing
 * @param {Readonly<{ name: string, colors: Record<string, string> }>} palette
 * @param {PaintStateSource['originalLegendOrder']} originalLegendOrder
 */
export const requestPaintOf = (drawing, palette, originalLegendOrder) => {
  /** @param {DrawingState} read */
  const of = (read) => editorPaintState({ appliedPaletteColors: { value: palette.colors }, originalLegendOrder }, read);
  return { palette, bound: of(drawing), of };
};

// A Generate keeps the fill edits that hold a color on a feature its Results
// draw (`replaceFeatureFills`): the paint it drew names the fills it keeps, so
// an edit it drops is no edit made since its request.
/** @param {EditorPaintState} paint @param {unknown} fills The fills the Generate keeps. @returns {EditorPaintState} */
export const withKeptFills = (paint, fills) => {
  const [rules, , legendFills] = /** @type {[unknown, unknown, unknown]} */ (JSON.parse(paint.colors[1]));
  return { ...paint, colors: [paint.colors[0], JSON.stringify([rules, fills, legendFills])] };
};

// The rule table a paint state shows, or null when it is unknown (a Result
// restored without the state it shows, `unknownPaint`).
/** @param {EditorPaintState} paint @returns {SpecificColorRule[] | null} */
export const paintRuleTable = (paint) => (paint.colors[1]
  ? /** @type {[SpecificColorRule[], unknown, unknown]} */ (JSON.parse(paint.colors[1]))[0]
  : null);

// Whether two paint states show the same palette and color edits.
/** @param {EditorPaintState['colors']} left @param {EditorPaintState['colors']} right */
export const sameColors = (left, right) => left[0] === right[0] && left[1] === right[1];

// A held paint whose palette gives every key the color of `palette`, the
// applied one, shows `palette` (a restore clones the palette it restores).
/** @param {EditorPaintState} paint @param {unknown} palette @returns {EditorPaintState} */
const withAppliedPalette = (paint, palette) => {
  const [held, edits] = paint.colors;
  return held !== palette && paletteColorsEqual(
    /** @type {Record<string, unknown> | null} */ (held), /** @type {Record<string, unknown> | null} */ (palette)
  ) ? { ...paint, colors: [palette, edits] } : paint;
};

// Whether two paint states show the same fills: palette and color edits, and visibility.
/** @param {Pick<EditorPaintState, 'colors' | 'visibility'>} left @param {Pick<EditorPaintState, 'colors' | 'visibility'>} right */
const sameFills = (left, right) => sameColors(left.colors, right.colors) && left.visibility === right.visibility;

/**
 * The label rows of `EditorPaintState.labels`: [recordKey, biologicalFeatureId,
 * labelVisibility, labelText, labelSourceText], and the bulk label texts by source text.
 * @typedef {[[string, string, string | null, string | null, string | null][], Record<string, string>]} LabelPaintRows
 */

// Whether a Result drawn with the label intent `drawn` lacks what `current`
// draws: a label On or Off that returned to Default, or an edited text that
// left. Only Python decides whether Default draws a label and which text it
// started from, so no projection of the Result's bytes can show them; another
// value (Off to On, one text to another) the label projection shows or asks for.
// An edit that had already left when the Result was last shown (`shown`, '' when
// unknown) is not asked about: the Result followed it, or asked once for it. A
// row of a feature the Result does not draw (`drawsFeature`) says nothing of it.
// The bytes hold no source text, so a bulk text counts for every Result.
/**
 * @param {string} drawn
 * @param {string} shown
 * @param {string} current
 * @param {(recordKey: string, featureId: string) => boolean} drawsFeature
 */
export const drawnLabelIntentGone = (drawn, shown, current, drawsFeature) => {
  const [drawnRows, drawnBulk] = /** @type {LabelPaintRows} */ (JSON.parse(drawn));
  const rowKey = (/** @type {LabelPaintRows[0][number]} */ [recordKey, featureId]) => JSON.stringify([recordKey, featureId]);
  const read = (/** @type {string} */ labels) => {
    const [rows, bulk] = /** @type {LabelPaintRows} */ (JSON.parse(labels));
    return { held: new Map(rows.map((row) => [rowKey(row), row])), bulk };
  };
  const now = read(current);
  const before = shown ? read(shown) : null;
  const shownByDefault = (/** @type {string | null | undefined} */ visibility) => visibility !== 'on' && visibility !== 'off';
  const rowGone = (/** @type {LabelPaintRows[0][number]} */ row, /** @type {ReturnType<typeof read>} */ { held }) => {
    const heldRow = held.get(rowKey(row));
    return (!shownByDefault(row[2]) && shownByDefault(heldRow?.[2])) || (row[3] !== null && (heldRow?.[3] ?? null) === null);
  };
  const bulkGone = (/** @type {string} */ source, /** @type {string} */ text, /** @type {ReturnType<typeof read>} */ { bulk }) => bulk[source] !== text;
  return drawnRows.some((row) => rowGone(row, now) && !(before && rowGone(row, before)) && drawsFeature(row[0], row[1]))
    || Object.entries(drawnBulk).some(([source, text]) => bulkGone(source, text, now) && !(before && bulkGone(source, text, before)));
};

// A paint state no editor state equals: a Result recorded with it shows
// every paint domain and the labels on its next display. The Legend order is
// kept: replaying the default order needs the Result's generated inventory.
/** @param {EditorPaintState} current @returns {EditorPaintState} */
const unknownPaint = (current) => ({ ...current, colors: [null, ''], visibility: '', strokes: '', labels: '' });

// The phases whose Results Python drew from the editor state of their bind:
// a Generate, a label rerender, and a committed candidate (a rotation, an
// alignment Apply). Every other commit (a Session Load, a History restore, a
// rollback) mounts bytes that were saved or kept.
const DRAW_PHASES = new Set(['generate', 'label-reflow', 'target-record-transform']);
/** @param {string} phase The phase of a mounted Result. */
export const drawsResult = (phase) => DRAW_PHASES.has(phase);

/**
 * What one Result's bytes hold of the editor state: `drawn`, the state Python
 * drew them with, which a live edit does not change (a rerender or a Generate
 * draws a new Result); `shown`, the state they show, which for the displayed
 * Result is the current state and the domains it lacks (`behind`), written
 * when it leaves; `left`, the intent signature it left a drawing with (E1),
 * or null.
 * @typedef {{ drawn: EditorPaintState, shown: EditorPaintState, left: string | null }} ResultEditorRecord
 */

/**
 * @typedef {object} ResultPaintRecordPorts
 * @property {(result: SvgResultEntry) => object} resultRuntimeKey
 *   The preview owner's key of a Result's runtime state, which every content
 *   version of the Result shares.
 */

/**
 * The editor state of each Result, held as long as the Result is referenced
 * (by a drawing, a History step, a rollback): one record per Result, keyed by
 * its runtime state, so a Result an Undo brings back keeps it. Never saved.
 * @param {ResultPaintRecordPorts} ports
 */
export const createResultPaintRecord = ({ resultRuntimeKey }) => {
  /** @type {WeakMap<object, ResultEditorRecord>} */
  const records = new WeakMap();
  /** @param {SvgResultEntry} result */
  const keyOf = (result) => toRaw(resultRuntimeKey(result));
  /** @param {SvgResultEntry} result */
  const recordOf = (result) => records.get(keyOf(result));
  /** @param {SvgResultEntry} result @param {EditorPaintState} drawn @param {EditorPaintState} [shown] */
  const recordNew = (result, drawn, shown = drawn) => {
    /** @type {ResultEditorRecord} */
    const record = { drawn, shown, left: null };
    records.set(keyOf(result), record);
    return record;
  };
  /** @type {ResultEditorRecord | null} */
  let bound = null;
  // The paint domains the displayed Result was not shown (a declined or
  // failed display); recorded with it when it leaves.
  /** @type {Partial<EditorPaintState> | null} */
  let behind = null;
  /** @param {EditorPaintState} state */
  const recordBound = (state) => {
    if (bound) bound.shown = { ...state, ...behind };
    behind = null;
  };
  return {
    /**
     * A committed artifact: each Result without a record gets one, and
     * `displayed` is the Result mounted by the commit. A commit that `draws`
     * (`drawsResult`) has Results Python drew from the paint it names, read
     * as its request was built: a label rerender draws the fills and visibility of the state it
     * started from (OV-346), so the displayed Result lacks those of an edit
     * made meanwhile (`lacking`) until they are shown. Any other commit (a
     * Session Load, a History restore, a rollback) keeps each Result's bytes
     * as they were saved or kept: a batch Result not displayed since the last
     * edits lacks them, so every Result but the mounted one shows every paint
     * domain on its first display, and was drawn with the state of the drawing
     * that holds it (`drawingStateOf`), which is all a Session or an Undo says
     * of its bytes (OV-3730). The Result displayed until now followed every
     * live edit until a draw; a restore may have restored the editor state
     * before the commit (an Undo of a Generate), so what that Result shows is
     * unknown then and its next display shows every domain.
     * @param {readonly SvgResultEntry[]} results The Results of all drawings.
     * @param {SvgResultEntry} displayed
     * @param {EditorPaintState} current
     * @param {{ draws?: EditorPaintState | null, drawingStateOf?: (result: SvgResultEntry) => EditorPaintState }} [options]
     *   `draws`: the paint the Results were drawn from, or null when the commit draws nothing.
     */
    commit(results, displayed, current, { draws = null, drawingStateOf = () => current } = {}) {
      const drawn = draws ?? current;
      const displayedKey = keyOf(displayed);
      if (bound && bound !== records.get(displayedKey)) {
        if (draws) recordBound(current);
        else bound.shown = unknownPaint(current);
      }
      [...results, displayed].forEach((result) => {
        if (recordOf(result)) return;
        if (draws) recordNew(result, drawn);
        else if (keyOf(result) === displayedKey) recordNew(result, current);
        else recordNew(result, drawingStateOf(result), unknownPaint(current));
      });
      bound = /** @type {ResultEditorRecord} */ (records.get(displayedKey)); // Recorded above.
      bound.left = null;
      const { colors, visibility } = drawn;
      behind = sameFills(drawn, current) ? null : { colors, visibility };
    },
    /**
     * Duplicate (OV-3732): each copied Result's bytes hold what its source's
     * hold, so its record starts as a copy of the source's; the displayed
     * Result followed every live edit (`current`, the copied drawing's state)
     * but what it lacks.
     * @param {readonly SvgResultEntry[]} sources
     * @param {readonly SvgResultEntry[]} copies The copies, in the order of `sources`.
     * @param {EditorPaintState} current
     */
    copy(sources, copies, current) {
      sources.forEach((source, index) => {
        const record = recordOf(source);
        if (record) records.set(keyOf(copies[index]), { ...record, shown: record === bound ? { ...current, ...behind } : record.shown, left: null });
      });
    },
    /**
     * The displayed Result leaves a drawing (a mode switch) having followed
     * every live edit: it keeps the intent it left with (`signature`, E1)
     * unless a declined display left it lacking a domain, which its next
     * display shows (OV-3731).
     * @param {EditorPaintState} current
     * @param {string} signature
     */
    depart(current, signature) {
      if (bound && !behind) bound.left = signature;
      recordBound(current);
      bound = null;
    },
    /**
     * The intent signature a Result left its drawing with, once: its next
     * display reads it (E1).
     * @param {SvgResultEntry} result
     */
    takeLeft(result) {
      const record = recordOf(result);
      const left = record?.left ?? null;
      if (record) record.left = null;
      return left;
    },
    /**
     * A Result is displayed: the Result shown until now followed every live
     * edit. Returns the state the arriving Result's bytes show.
     * @param {SvgResultEntry} result
     * @param {EditorPaintState} current
     * @returns {EditorPaintState}
     */
    display(result, current) {
      const record = recordOf(result) ?? recordNew(result, current);
      if (bound !== record) recordBound(current);
      bound = record;
      return record.shown;
    },
    /**
     * The displayed Result now shows `current`, but the fills and visibility
     * of `previous` when its palette and rules projection `declined`, and
     * every paint domain of `previous` when its projection or compile
     * `failed` (its labels follow in the binder's label step).
     * @param {SvgResultEntry} result
     * @param {EditorPaintState} current
     * @param {EditorPaintState} previous
     * @param {{ declined?: boolean, failed?: boolean }} [options]
     */
    shown(result, current, previous, { declined = false, failed = false } = {}) {
      const { colors, visibility, strokes, legendOrder } = previous;
      behind = failed ? { colors, visibility, strokes, legendOrder } : (declined ? { colors, visibility } : null);
      (recordOf(result) ?? recordNew(result, current)).shown = { ...current, ...behind };
    },
    /**
     * What the displayed Result still shows of the state before a declined
     * or failed display, or null (shown in full, or not displayed).
     * @param {SvgResultEntry} result
     * @returns {Partial<EditorPaintState> | null}
     */
    lacking(result) {
      return bound && bound === recordOf(result) ? behind : null;
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
     * The editor state the Result's bytes were drawn with.
     * @param {SvgResultEntry} result
     * @returns {EditorPaintState | undefined}
     */
    drawn(result) {
      return recordOf(result)?.drawn;
    },
    /**
     * Whether the displayed Result was drawn with a Feature visibility other than `visibility`: a live edit no
     * render drew since, so its input summary describes the render, not the figure (false while no Result is
     * displayed).
     * @param {string} visibility `visibilityPaintState` of the drawing that shows it.
     */
    displayedDrewOtherVisibility(visibility) {
      return bound !== null && bound.drawn.visibility !== visibility;
    },
    /**
     * The Label visibility On and Off the displayed Result was drawn with, by
     * feature identity key, or null when no Result is displayed. Only Python
     * decides whether Default draws a label, so the label owner asks for the
     * rerender when one of them returned to Default (UJ-01, R3), or when an
     * On row lacks its label on bytes not drawn with it (OV-2413).
     * @returns {ReadonlyMap<string, 'on' | 'off'> | null}
     */
    displayedDrawnLabelVisibility() {
      if (!bound) return null;
      const [rows] = /** @type {LabelPaintRows} */ (JSON.parse(bound.drawn.labels));
      return new Map(rows.flatMap(([recordKey, featureId, visibility]) => (
        visibility === 'on' || visibility === 'off' ? [/** @type {[string, 'on' | 'off']} */ ([featureIdentityKey(recordKey, featureId), visibility])] : []
      )));
    },
    /**
     * A Generate that moves a record's rows to its new key (OV-1949) moves them
     * after its Results were committed; the request drew the moved rows
     * (`moved`), which the Results' record takes then, so the moved rows ask
     * for no rerender (OV-1983).
     * @param {readonly SvgResultEntry[]} results The Results the Generate committed.
     * @param {EditorPaintState} moved
     */
    moveDrawn(results, moved) {
      results.forEach((result) => {
        const record = recordOf(result);
        if (record) record.drawn = moved;
      });
    },
    /**
     * What the bytes of `results` show now, in their order, for a figure
     * History holds (OV-2436): a Result not displayed shows its record's
     * `shown` (one without a record is recorded by its next display, as
     * before); the displayed Result shows the current state but what it lacks, or,
     * when `drawing` (a rerender's figure as its Legend is read, before an
     * edit made while it ran is shown on it), the paint it was drawn with.
     * @param {readonly SvgResultEntry[]} results
     * @param {EditorPaintState} current
     * @param {boolean} drawing
     * @returns {EditorPaintState[]}
     */
    heldPaint(results, current, drawing) {
      return results.map((result) => {
        const record = recordOf(result);
        if (!record) return unknownPaint(current);
        if (record !== bound) return record.shown;
        return drawing ? record.drawn : { ...current, ...behind };
      });
    },
    /**
     * A figure History holds shows again (OV-2436): its Results' bytes show
     * `paints` (`heldPaint`, in their order), whatever their records said
     * since, and left no drawing, so their next display shows the edits
     * since. A held palette that gives every key the color of `applied`, the
     * applied palette the restored intent left, is that palette: the restore
     * clones it, and the same colors ask for no fills. Called once the
     * displayed Result has left (`depart`).
     * @param {readonly SvgResultEntry[]} results
     * @param {readonly EditorPaintState[]} paints
     * @param {unknown} applied
     */
    arrive(results, paints, applied) {
      results.forEach((result, index) => {
        const record = recordOf(result);
        if (record) Object.assign(record, { shown: withAppliedPalette(paints[index], toRaw(applied)), left: null });
      });
    },
    /** @param {SvgResultEntry} result The Result shown again with the intent it left with (E1). */
    rebind(result) {
      bound = recordOf(result) ?? null;
      behind = null;
    }
  };
};
