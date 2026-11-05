// @ts-check
// A drawing's Legend draft as Python's rows: the listed and dormant entries, Python's keys of the generated and the
// deleted rows, and the row a style addresses by its key. The Generate compiler (`compilePlanBundle` in
// app/candidate-render.js) and the feature popup read the rows through `draftLegendRows`, and the request projection
// reads the edits Python draws from the same rows (`legendEditConfigLeaves`): the renames, deletions, and order of
// the rows Python draws, and the canvas padding, as the config leaves `objects.legends.row_names`,
// `objects.legends.hidden_rows`, `objects.legends.row_order`, and `canvas.padding` (gbdraw/config/models/objects.py
// and canvas.py read them; issue #1427, PD-OI-114 scenario revision 2).
import { isLegendOrderEdited } from './legend-svg.js';
import { ruleLegendCaptions } from './specific-color-rules.js';
import { normalizeCanvasPadding } from './legend-layout.js';
/** @import { PythonLegendKey } from './legend-svg.js' */
/** @import { SpecificColorRule } from './specific-color-rules.js' */
/** @import { CanvasPadding } from './drawing/canvas.js' */

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
