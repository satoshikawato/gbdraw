// @ts-check
// A drawing's Legend draft as Python's rows: the listed and dormant entries, Python's keys of the generated and the
// deleted rows, and the row a style addresses by its key. The Generate compiler (`compilePlanBundle` in
// app/candidate-render.js) and the feature popup read the rows through `draftLegendRows`, and the request projection
// reads the edits Python draws from the same rows (`legendEditConfigLeaves`): the renames, deletions, and order of
// the rows Python draws, and the canvas padding, as the config leaves `objects.legends.row_names`,
// `objects.legends.hidden_rows`, `objects.legends.row_order`, and `canvas.padding` (gbdraw/config/models/objects.py
// and canvas.py read them; issue #1427, PD-OI-114 scenario revision 2).
import {
  getVisibleFeatureLegendGroup, isLegendOrderEdited, pythonLegendKey, pythonOrderedLegendRows, resultLegendRows
} from './legend-svg.js';
import { ruleLegendCaptions } from './specific-color-rules.js';
import { normalizeCanvasPadding } from './legend-layout.js';
import { refitFieldsBeforeStep } from './setting-application.js';
/** @import { PythonLegendKey } from './legend-svg.js' */
/** @import { SpecificColorRule } from './specific-color-rules.js' */
/** @import { CanvasPadding } from './drawing/canvas.js' */
/** @import { HistoryChangeLike, RefitField } from './setting-application.js' */

/**
 * A Legend entry as the compile reads it: a listed, deleted, or dormant row
 * (`LegendEntry`), or one a Session saved. Every field is unvalidated.
 * @typedef {{ caption?: unknown, originalCaption?: unknown, color?: unknown, xPos?: unknown, yPos?: unknown, featureIds?: unknown }} LegendEntryInput
 */

/**
 * A Legend entry after `normalizedLegendEntries`.
 * @typedef {{ caption: string, originalCaption: string, color: string, xPos: number | null, yPos: number | null, featureIds: readonly string[] }} NormalizedLegendEntry
 */

/**
 * The Legend lists of a draft, as the compile and the request projection read them.
 * @typedef {object} DraftLegendLists
 * @property {readonly (LegendEntryInput | null | undefined)[]} [legendEntries]
 * @property {readonly (LegendEntryInput | null | undefined)[]} [deletedLegendEntries]
 * @property {readonly (LegendEntryInput | null | undefined)[]} [dormantLegendEntries]
 * @property {readonly string[]} [originalLegendOrder] The shown Result's generated inventory (Python's keys).
 */

/** @param {unknown} value */
const text = (value) => String(value ?? '').trim();

/**
 * @param {readonly (LegendEntryInput | null | undefined)[] | null | undefined} entries
 * @returns {NormalizedLegendEntry[]}
 */
export const normalizedLegendEntries = (entries) => (
  Array.isArray(entries)
    ? entries.map((entry) => ({
        caption: text(entry?.caption),
        originalCaption: text(entry?.originalCaption || entry?.caption),
        color: text(entry?.color),
        xPos: Number.isFinite(Number(entry?.xPos)) ? Number(entry.xPos) : null,
        yPos: Number.isFinite(Number(entry?.yPos)) ? Number(entry.yPos) : null,
        featureIds: Object.freeze([
          ...new Set((Array.isArray(entry?.featureIds) ? entry.featureIds : []).map(text).filter(Boolean))
        ])
      })).filter((entry) => entry.caption)
    : []
);

// The draft's Legend rows: the listed and the dormant entries, Python's keys
// of the generated and the deleted rows, and the row a style addresses by its
// row key (a key of `legendColorOverrides` or `legendStrokeOverrides`, never
// the name the row shows: two rows may show one name); null for a deleted row.
// The compile and the feature popup read the rows through it, so a deleted
// row's stroke reaches no feature in either.
/** @param {DraftLegendLists} legend */
export const draftLegendRows = ({
  legendEntries = [], deletedLegendEntries = [], dormantLegendEntries = [], originalLegendOrder = []
}) => {
  const currentEntries = normalizedLegendEntries(legendEntries);
  const originalCaptions = new Set(
    (Array.isArray(originalLegendOrder) ? originalLegendOrder : []).map(text).filter(Boolean)
  );
  const deletedCaptions = new Set(
    (Array.isArray(deletedLegendEntries) ? deletedLegendEntries : [])
      .map((entry) => text(entry?.originalCaption || entry?.caption))
      .filter((caption) => originalCaptions.has(caption))
  );
  // OV-120: a renamed row an earlier Generate hid (GC off, Show Depth off)
  // waits in the drawing; this Result may draw it again or still leave it out.
  const shownOriginals = new Set(currentEntries.map((entry) => entry.originalCaption));
  const dormantEntries = normalizedLegendEntries(dormantLegendEntries)
    .filter((entry) => entry.caption !== entry.originalCaption && !shownOriginals.has(entry.originalCaption)
      && !deletedCaptions.has(entry.originalCaption));
  const dormantOriginals = new Set(dormantEntries.map((entry) => entry.originalCaption));
  const entriesByKey = new Map([...dormantEntries, ...currentEntries].map(entry => [entry.originalCaption, entry]));
  /** @param {string} key */
  const styledRow = (key) => {
    if (deletedCaptions.has(key)) return null;
    const dormant = dormantOriginals.has(key);
    return {
      entry: entriesByKey.get(key), dormant, isOriginal: dormant || originalCaptions.has(key), targetCaption: /** @type {PythonLegendKey} */ (key)
    };
  };
  return { currentEntries, originalCaptions, deletedCaptions, dormantEntries, styledRow };
};

/** The config leaves of the Legend row edits and the canvas padding (gbdraw/config/models). */
export const LEGEND_EDIT_CONFIG_PATHS = Object.freeze(/** @type {const} */ ({
  rowNames: 'objects.legends.row_names',
  hiddenRows: 'objects.legends.hidden_rows',
  rowOrder: 'objects.legends.row_order',
  canvasPadding: 'canvas.padding'
}));

/**
 * The config leaves a request carries for a drawing's Legend row edits and canvas padding. A leaf is present only
 * when it is not empty, so a drawing without these edits sends the request it sent before them: the renamed rows
 * Python draws (`[key, name]`), the deleted rows Python draws, the rows Python draws in their edited order (only
 * when the order is edited), all by key, and the padding (top, right, bottom, left, px; only when a side is not 0).
 * @typedef {{
 *   'objects.legends.row_names'?: Array<[PythonLegendKey, string]>,
 *   'objects.legends.hidden_rows'?: PythonLegendKey[],
 *   'objects.legends.row_order'?: PythonLegendKey[],
 *   'canvas.padding'?: [number, number, number, number]
 * }} LegendEditConfigLeaves
 */

/**
 * The draft a request projects into `LegendEditConfigLeaves`.
 * @typedef {DraftLegendLists & { rules: readonly Partial<SpecificColorRule>[], canvasPadding: Readonly<Partial<Record<keyof CanvasPadding, unknown>>> }} LegendEditDraft
 */

/**
 * The one projection of a drawing's Legend row edits and canvas padding into the config leaves Python draws them
 * from (Generate's request and the automatic rerender's). A row is named by its Python key (`originalCaption`),
 * never by its caption. Only rows the shown Result's inventory holds are Python's: an editor row stays the Web's.
 * A rule's row is renamed through its rule's caption, which the color table carries, so its name is not sent. The
 * canvas padding is the padding the Web draws (`normalizeCanvasPadding`).
 * @param {LegendEditDraft} draft
 * @returns {LegendEditConfigLeaves}
 */
export const legendEditConfigLeaves = ({ rules, canvasPadding, ...legend }) => {
  const rows = draftLegendRows(legend);
  const ruleRowOf = ruleLegendCaptions({
    rules: [...rules], pythonRows: new Map(), features: [], originalLegendOrder: [...rows.originalCaptions]
  });
  const ruleRowKeys = new Set(rules.filter((rule) => text(rule?.cap)).map(ruleRowOf));
  const pythonRows = rows.currentEntries.filter((entry) => rows.originalCaptions.has(entry.originalCaption));
  const key = (/** @type {NormalizedLegendEntry} */ entry) => /** @type {PythonLegendKey} */ (entry.originalCaption);
  const names = pythonRows
    .filter((entry) => entry.caption !== entry.originalCaption && !ruleRowKeys.has(key(entry)))
    .map((entry) => /** @type {[PythonLegendKey, string]} */ ([key(entry), entry.caption]));
  const hidden = /** @type {PythonLegendKey[]} */ ([...rows.deletedCaptions]);
  const order = isLegendOrderEdited(rows.currentEntries, legend.originalLegendOrder) ? pythonRows.map(key) : [];
  const padding = normalizeCanvasPadding(canvasPadding);
  const sides = /** @type {[number, number, number, number]} */ ([padding.top, padding.right, padding.bottom, padding.left]);
  /** @type {LegendEditConfigLeaves} */
  const leaves = {};
  if (names.length) leaves[LEGEND_EDIT_CONFIG_PATHS.rowNames] = names;
  if (hidden.length) leaves[LEGEND_EDIT_CONFIG_PATHS.hiddenRows] = hidden;
  if (order.length) leaves[LEGEND_EDIT_CONFIG_PATHS.rowOrder] = order;
  if (sides.some((side) => side !== 0)) leaves[LEGEND_EDIT_CONFIG_PATHS.canvasPadding] = sides;
  return leaves;
};

/**
 * What Python draws and sizes of a Legend draft and its canvas padding, as the request's config leaves carry it
 * (`legendEditConfigLeaves`), with a rule row's name too: without rules the projection names every renamed row
 * Python draws, and a rule row's name reaches Python through its rule. The color rules' captions and what they match
 * are part of it (OV-2556: a Legend rename of a feature type row can write rules of the new name, which Python draws as
 * a row the request's leaves do not name). Colors, strokes, and the rows only the Web draws are not part of it.
 * @param {Readonly<Record<RefitField, unknown>> & { originalLegendOrder: readonly string[] }} draft
 * @returns {string}
 */
const drawnLegendKey = ({ legendEntries, deletedLegendEntries, rules, originalLegendOrder, canvasPadding }) => {
  // Read as lists here, not handed to a helper: the paint intent writer guard counts a paint member passed to a
  // function it does not know as a write.
  /** @type {readonly (LegendEntryInput | null | undefined)[]} */
  const listed = Array.isArray(legendEntries) ? legendEntries : [];
  /** @type {readonly (LegendEntryInput | null | undefined)[]} */
  const deleted = Array.isArray(deletedLegendEntries) ? deletedLegendEntries : [];
  const padding = /** @type {Readonly<Partial<Record<keyof CanvasPadding, unknown>>>} */ (
    canvasPadding && typeof canvasPadding === 'object' ? canvasPadding : {});
  /** @type {readonly (Partial<SpecificColorRule> | null | undefined)[]} */
  const ruleRows = Array.isArray(rules) ? rules : [];
  return JSON.stringify([legendEditConfigLeaves({
    rules: [], canvasPadding: padding, legendEntries: listed, deletedLegendEntries: deleted, originalLegendOrder
  }), ruleRows.map((rule) => [rule?.feat, rule?.qual, rule?.val, rule?.cap])]);
};

/**
 * Whether a History step of the shown drawing changed what Python draws and sizes of its Legend rows or its canvas
 * padding (`drawnLegendKey`): a rename (a rule row's too), a deletion, a Restore, a Move, a Sort, or a padding change,
 * and the Undo or Redo of one; not a color or stroke edit. On a figure drawn to a requested size such a step asks the
 * automatic rerender to fit the figure again (PD-OI-114 rev 2, Q-D1 A).
 * @param {readonly HistoryChangeLike[] | null | undefined} changes The step's changes.
 * @param {'undo' | 'redo' | undefined} direction Undefined for the step's commit.
 * @param {Readonly<Record<RefitField, unknown>> & { originalLegendOrder: readonly string[] }} current The drawing's
 *   values as the step left them, and the shown Result's inventory.
 */
export const stepRefitsDrawnLegend = (changes, direction, current) => {
  const other = refitFieldsBeforeStep(changes, current, direction === 'undo');
  return other !== null
    && drawnLegendKey(current) !== drawnLegendKey({ ...other, originalLegendOrder: current.originalLegendOrder });
};

/** @param {unknown} value @returns {value is string} */
const isRowText = (value) => typeof value === 'string' && value.trim() !== '';

/**
 * The Legend row edit and canvas padding leaves a request carries, by path: its `configOverrides` leaf (a Web
 * request), else its full `config`'s (a CLI or Python Session's request holds the whole config). The one reader of
 * the leaves in a stored request. A leaf of another shape is no edit: Python refuses it, so no Result drew it.
 * @param {unknown} renderRequest
 * @returns {LegendEditConfigLeaves}
 */
export const readLegendEditLeaves = (renderRequest) => {
  const options = /** @type {{ diagramOptions?: { configOverrides?: unknown, config?: unknown } } | null} */ (
    renderRequest && typeof renderRequest === 'object' ? renderRequest : null
  )?.diagramOptions;
  const overrides = options?.configOverrides && typeof options.configOverrides === 'object' ? options.configOverrides : {};
  /** @param {string} path @returns {unknown} */
  const leaf = (path) => {
    if (Object.hasOwn(overrides, path)) return /** @type {Record<string, unknown>} */ (overrides)[path];
    /** @type {unknown} */
    let value = options?.config;
    for (const segment of path.split('.')) {
      value = value && typeof value === 'object' && !Array.isArray(value) && Object.hasOwn(value, segment)
        ? /** @type {Record<string, unknown>} */ (value)[segment] : undefined;
    }
    return value;
  };
  /** @param {unknown} value @returns {PythonLegendKey[]} */
  const keys = (value) => (Array.isArray(value) && value.every(isRowText) ? /** @type {PythonLegendKey[]} */ ([...value]) : []);
  const names = leaf(LEGEND_EDIT_CONFIG_PATHS.rowNames);
  const padding = leaf(LEGEND_EDIT_CONFIG_PATHS.canvasPadding);
  /** @type {LegendEditConfigLeaves} */
  const leaves = {};
  if (Array.isArray(names) && names.length
    && names.every((pair) => Array.isArray(pair) && pair.length === 2 && isRowText(pair[0]) && isRowText(pair[1]))) {
    leaves[LEGEND_EDIT_CONFIG_PATHS.rowNames] = names.map(([key, name]) => /** @type {[PythonLegendKey, string]} */ ([key, name]));
  }
  const hidden = keys(leaf(LEGEND_EDIT_CONFIG_PATHS.hiddenRows));
  if (hidden.length) leaves[LEGEND_EDIT_CONFIG_PATHS.hiddenRows] = hidden;
  const order = keys(leaf(LEGEND_EDIT_CONFIG_PATHS.rowOrder));
  if (order.length) leaves[LEGEND_EDIT_CONFIG_PATHS.rowOrder] = order;
  if (Array.isArray(padding) && padding.length === 4
    && padding.every((side) => typeof side === 'number' && Number.isFinite(side) && side >= 0) && padding.some((side) => side !== 0)) {
    leaves[LEGEND_EDIT_CONFIG_PATHS.canvasPadding] = /** @type {[number, number, number, number]} */ ([...padding]);
  }
  return leaves;
};

/**
 * A Legend row as a Session's `editorState.legend` lists it.
 * @typedef {{ caption: string, originalCaption: PythonLegendKey, color: string, xPos: number, yPos: number, featureIds: string[] }} LeafLegendEntry
 */

/**
 * The Web editor state of a drawing whose Session holds its Legend row edits and canvas padding only as request
 * leaves (a CLI or Python Session): the Legend lists (`editorState.legend`: the rows Python shows, in their order,
 * with their names; the hidden rows; Python's order of its rows) when a Legend leaf is set and the shown Result has
 * a Legend, and the padding (`ui.canvasPadding`) when it is set. Python drew the leaves into the Result, so the
 * rows, swatches, and order are read from it, and the names and hidden rows from the leaves.
 * @typedef {object} LeafEditorState
 * @property {{ entries: LeafLegendEntry[], deletedEntries: LeafLegendEntry[], originalOrder: PythonLegendKey[] } | null} legend
 * @property {CanvasPadding | null} canvasPadding
 */

/**
 * @param {LegendEditConfigLeaves} leaves `readLegendEditLeaves` of the drawing's request.
 * @param {SVGSVGElement | null} svg The drawing's shown Result.
 * @returns {LeafEditorState}
 */
export const editorStateOfLegendEditLeaves = (leaves, svg) => {
  const sides = leaves[LEGEND_EDIT_CONFIG_PATHS.canvasPadding];
  const canvasPadding = sides ? { top: sides[0], right: sides[1], bottom: sides[2], left: sides[3] } : null;
  const names = new Map(leaves[LEGEND_EDIT_CONFIG_PATHS.rowNames] ?? []);
  const hidden = leaves[LEGEND_EDIT_CONFIG_PATHS.hiddenRows] ?? [];
  const group = getVisibleFeatureLegendGroup(svg);
  const rows = (resultLegendRows(svg) ?? []).filter((row) => !row.editor);
  if (!group || rows.length === 0 || (names.size === 0 && hidden.length === 0 && !leaves[LEGEND_EDIT_CONFIG_PATHS.rowOrder])) {
    return { legend: null, canvasPadding };
  }
  const hiddenKeys = new Set(hidden);
  /** @param {typeof rows[number]} row @returns {LeafLegendEntry} */
  const entryOf = (row) => ({
    caption: names.get(row.key) ?? row.key, originalCaption: row.key, color: row.color, xPos: row.xPos, yPos: row.yPos, featureIds: []
  });
  const rowOf = new Map(rows.map((row) => [row.key, row]));
  return {
    legend: {
      entries: rows.filter((row) => row.pythonShown && !hiddenKeys.has(row.key)).map(entryOf),
      deletedEntries: [...hiddenKeys].flatMap((key) => { const row = rowOf.get(key); return row ? [entryOf(row)] : []; }),
      originalOrder: [...new Set(pythonOrderedLegendRows(group).map((row) => /** @type {PythonLegendKey} */ (pythonLegendKey(row))))]
    },
    canvasPadding
  };
};
