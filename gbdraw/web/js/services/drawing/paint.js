// @ts-check
// The paint slice of a drawing (R2): its Legend rows (`legendEntries`), the
// rows a delete keeps for Restore and Generate (`deletedLegendEntries`), the
// renamed rows a Generate did not draw (`dormantLegendEntries`, OV-120), the
// captions of the rows the color rules add (`addedLegendCaptions`), and the
// Legend and per-feature color and stroke edits (Session 46
// `modes.<m>.editorState.legend`, `editorState.featureStrokes`, and
// `features.featureColorOverrides`). The drawing shows them read-only to the
// type checker, rows and fields included; only this module writes them,
// through `drawing.slices.paint`: the Legend and feature editors' writers
// below, which `PAINT_INTENT_WRITERS` names by route
// (tests/web/feature-visibility-rule-writers.test.mjs), and Load, a History
// restore, and Reset (`install`, `reset`). The rows belong to the slice
// (D-15-44): a row is addressed by its key, Python's key of the row
// (`legendEntryKey`), never by its caption or color, which are display values
// a rename and a Legend color write. A caller passes rows and edits it checked.
import { ref, reactive } from '../../utils/reactivity.js';
import { cloneJsonData } from '../json-clone.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';
/** @import { DeepReadonly, ReadonlyRef, Ref } from '../../utils/reactivity.js' */
/** @import { SliceHandle } from './slice-handle.js' */
/** @import { PythonLegendKey, RenderedFeatureId } from '../legend-svg.js' */
/** @import { StrokeOverride } from '../style-overrides.js' */

/** @typedef {string & { readonly __brand: 'LegendCaption' }} LegendCaption The caption a Legend row shows, which a rename writes. */
/**
 * The color a Legend row records, never its identity (OV-288): one the user
 * set (a Legend color, an editor row's color, the color of a rule's row), or,
 * for a row the palette colors, the swatch the displayed Result showed when
 * its Legend list was last extracted (a Generate, a Result display, Settings
 * Reset). That can be a draft palette color, not Python's swatch, and a
 * palette, Default colors or Reset edit does not refresh it by itself. A row
 * shows `displayedLegendRowColors` (app/candidate-render.js; OV-282, D-26).
 * @typedef {string & { readonly __brand: 'RecordedLegendColor' }} RecordedLegendColor
 */
/**
 * A row of a drawing's Legend list. Its identity is Python's key
 * (`originalCaption`; an editor row's own caption); the caption and the
 * recorded color are display values their writers own (rename, Legend color),
 * never a row's identity or the features it reaches (OV-288).
 * @typedef {object} LegendEntry
 * @property {LegendCaption} caption
 * @property {PythonLegendKey} originalCaption
 * @property {RecordedLegendColor} color
 * @property {number} xPos
 * @property {number} yPos
 * @property {RenderedFeatureId[]} featureIds The features the row lists (older Sessions); empty: Python's row decides.
 */
/**
 * A row the drawing keeps off its list: a deleted row (its listed row as it
 * was) or a dormant rename; neither has a place in the drawn Legend.
 * @typedef {Omit<LegendEntry, 'xPos' | 'yPos'> & Partial<Pick<LegendEntry, 'xPos' | 'yPos'>>} KeptLegendEntry
 */
/** @typedef {string | { color: string, caption?: string }} FeatureFill A feature's color edit (an older Session's: a bare color). */
/**
 * @typedef {{
 *   legendEntries: Ref<LegendEntry[]>, deletedLegendEntries: Ref<KeptLegendEntry[]>,
 *   dormantLegendEntries: Ref<KeptLegendEntry[]>, addedLegendCaptions: Ref<Set<string>>,
 *   legendColorOverrides: Record<string, string>, legendStrokeOverrides: Record<string, StrokeOverride>,
 *   featureColorOverrides: Record<string, FeatureFill>, featureStrokeOverrides: Record<string, StrokeOverride>
 * }} PaintMembers
 */
/**
 * What `capture` returns and `install` writes (the added captions as a list).
 * @typedef {{
 *   legendEntries: LegendEntry[], deletedLegendEntries: KeptLegendEntry[], dormantLegendEntries: KeptLegendEntry[],
 *   addedLegendCaptions: string[], legendColorOverrides: Record<string, string>,
 *   legendStrokeOverrides: Record<string, StrokeOverride>, featureColorOverrides: Record<string, FeatureFill>,
 *   featureStrokeOverrides: Record<string, StrokeOverride>
 * }} PaintData
 */
/**
 * @typedef {{
 *   readonly legendEntries: ReadonlyRef<LegendEntry[]>, readonly deletedLegendEntries: ReadonlyRef<KeptLegendEntry[]>,
 *   readonly dormantLegendEntries: ReadonlyRef<KeptLegendEntry[]>, readonly addedLegendCaptions: ReadonlyRef<Set<string>>,
 *   readonly legendColorOverrides: Readonly<Record<string, string>>,
 *   readonly legendStrokeOverrides: DeepReadonly<Record<string, StrokeOverride>>,
 *   readonly featureColorOverrides: DeepReadonly<Record<string, FeatureFill>>,
 *   readonly featureStrokeOverrides: DeepReadonly<Record<string, StrokeOverride>>
 * }} PaintView
 */
/** @typedef {SliceHandle<'paint'>} PaintSlice The handle a writer takes. */
/** @typedef {DeepReadonly<LegendEntry>} LegendRow A listed row as the drawing shows it. */
/** @typedef {DeepReadonly<KeptLegendEntry>} KeptLegendRow A deleted or dormant row as the drawing shows it. */
/** @typedef {DeepReadonly<{ color: string, caption?: string }>} FeatureFillRecord A feature color edit with its caption. */

// The key of a listed Legend row: its key record (`originalCaption`), else
// its caption (a row saved before U3a has no record).
/** @param {{ readonly caption?: unknown, readonly originalCaption?: unknown }} entry @returns {PythonLegendKey} */
export const legendEntryKey = (entry) => /** @type {PythonLegendKey} */ (String(entry.originalCaption || entry.caption || '').trim());

// A feature color edit as a record: an older Session's bare color has no
// caption and reads as no record, as the editors read it.
/** @param {DeepReadonly<FeatureFill> | null | undefined} fill @returns {FeatureFillRecord | null} */
export const fillRecordOf = (fill) => (fill && typeof fill === 'object' ? fill : null);

// A paint handle exists only through `create`, which gives it these members.
/** @param {PaintSlice} slice */
const membersOf = (slice) => /** @type {PaintMembers} */ (sliceMembers(slice, 'paint'));

/**
 * @template T
 * @param {Record<string, T>} target
 * @param {Readonly<Record<string, T>>} entries Copied first, so `target` may be `entries`.
 */
const replaceEntries = (target, entries) => {
  const copy = cloneJsonData({ ...entries });
  Object.keys(target).forEach((key) => { delete target[key]; });
  Object.assign(target, copy);
};

// Writable copies of rows a caller passes (the casts: a copy of a row is one).
/** @param {readonly LegendRow[]} rows */
const listedCopy = (rows) => /** @type {LegendEntry[]} */ (/** @type {unknown} */ (cloneJsonData([...rows])));
/** @param {readonly KeptLegendRow[]} rows */
const keptCopy = (rows) => /** @type {KeptLegendEntry[]} */ (/** @type {unknown} */ (cloneJsonData([...rows])));

/**
 * A new slice, or one over a drawing's existing members (a test fixture).
 * @param {Partial<PaintMembers>} [members]
 * @returns {{ handle: PaintSlice, view: PaintView }}
 */
const create = (members = {}) => {
  /** @type {PaintMembers} */
  const data = {
    featureColorOverrides: members.featureColorOverrides ?? reactive(/** @type {Record<string, FeatureFill>} */ ({})),
    featureStrokeOverrides: members.featureStrokeOverrides ?? reactive(/** @type {Record<string, StrokeOverride>} */ ({})),
    legendEntries: members.legendEntries ?? ref(/** @type {LegendEntry[]} */ ([])),
    deletedLegendEntries: members.deletedLegendEntries ?? ref(/** @type {KeptLegendEntry[]} */ ([])),
    dormantLegendEntries: members.dormantLegendEntries ?? ref(/** @type {KeptLegendEntry[]} */ ([])),
    legendColorOverrides: members.legendColorOverrides ?? reactive(/** @type {Record<string, string>} */ ({})),
    legendStrokeOverrides: members.legendStrokeOverrides ?? reactive(/** @type {Record<string, StrokeOverride>} */ ({})),
    addedLegendCaptions: members.addedLegendCaptions ?? ref(/** @type {Set<string>} */ (new Set()))
  };
  return { handle: createSliceHandle('paint', data), view: Object.freeze({ ...data }) };
};

/**
 * The paint, read-only, for a writer's caller that holds only the handle.
 * @param {PaintSlice} slice
 * @returns {PaintView}
 */
export const paintOf = (slice) => membersOf(slice);

// What Save and History capture: the Legend paint (`editorState`) and the
// feature color edits (`features`) apart, so each copies only its part.
/**
 * @param {PaintSlice} slice
 * @returns {Omit<PaintData, 'featureColorOverrides'>}
 */
export const captureLegendPaint = (slice) => {
  const members = membersOf(slice);
  return {
    featureStrokeOverrides: cloneJsonData({ ...members.featureStrokeOverrides }),
    legendEntries: cloneJsonData([...members.legendEntries.value]),
    deletedLegendEntries: cloneJsonData([...members.deletedLegendEntries.value]),
    dormantLegendEntries: cloneJsonData([...members.dormantLegendEntries.value]),
    legendColorOverrides: { ...members.legendColorOverrides },
    legendStrokeOverrides: cloneJsonData({ ...members.legendStrokeOverrides }),
    addedLegendCaptions: [...members.addedLegendCaptions.value]
  };
};

/**
 * @param {PaintSlice} slice
 * @returns {Record<string, FeatureFill>}
 */
export const captureFeatureFills = (slice) => cloneJsonData({ ...membersOf(slice).featureColorOverrides });

/**
 * @param {PaintSlice} slice
 * @returns {PaintData}
 */
const capture = (slice) => ({ featureColorOverrides: captureFeatureFills(slice), ...captureLegendPaint(slice) });

// Each member `data` holds becomes a copy of it; the others are kept (a
// Session's feature colors and its Legend edits are read apart).
/**
 * @param {PaintSlice} slice
 * @param {Partial<DeepReadonly<PaintData>>} data
 */
const install = (slice, data) => {
  const members = membersOf(slice);
  if (data.legendEntries) members.legendEntries.value = listedCopy(data.legendEntries);
  if (data.deletedLegendEntries) members.deletedLegendEntries.value = keptCopy(data.deletedLegendEntries);
  if (data.dormantLegendEntries) members.dormantLegendEntries.value = keptCopy(data.dormantLegendEntries);
  if (data.addedLegendCaptions) members.addedLegendCaptions.value = new Set(data.addedLegendCaptions);
  if (data.legendColorOverrides) replaceEntries(members.legendColorOverrides, data.legendColorOverrides);
  // A copy of a stroke or fill map is one (the casts: the copy is writable).
  if (data.legendStrokeOverrides) {
    replaceEntries(members.legendStrokeOverrides, /** @type {Record<string, StrokeOverride>} */ (data.legendStrokeOverrides));
  }
  if (data.featureColorOverrides) {
    replaceEntries(members.featureColorOverrides, /** @type {Record<string, FeatureFill>} */ (data.featureColorOverrides));
  }
  if (data.featureStrokeOverrides) {
    replaceEntries(members.featureStrokeOverrides, /** @type {Record<string, StrokeOverride>} */ (data.featureStrokeOverrides));
  }
};

/** @param {PaintSlice} slice */
const reset = (slice) => install(slice, {
  legendEntries: [], deletedLegendEntries: [], dormantLegendEntries: [], addedLegendCaptions: [],
  legendColorOverrides: {}, legendStrokeOverrides: {}, featureColorOverrides: {}, featureStrokeOverrides: {}
});

// ---- Legend rows ----

/**
 * The Legend list, in its order: the rows a Result lists (a Generate, a
 * Result display, a History step on another Result), the rows a rule commit
 * leaves, or the rows a Sort or Move ordered. A row of the list stays the
 * row object it is; a new row is the caller's.
 * @param {PaintSlice} slice
 * @param {readonly LegendRow[]} rows
 */
export const replaceLegendRows = (slice, rows) => {
  // A row of the list is a list row (the cast: the slice owns the rows).
  membersOf(slice).legendEntries.value = /** @type {LegendEntry[]} */ ([...rows]);
};

/**
 * Deletes the listed row at `index`: the deleted list keeps it.
 * @param {PaintSlice} slice
 * @param {number} index
 * @returns {boolean} Whether a row was deleted.
 */
export const deleteLegendRow = (slice, index) => {
  const { legendEntries, deletedLegendEntries } = membersOf(slice);
  const entry = legendEntries.value[index];
  if (!entry) return false;
  deletedLegendEntries.value = [...deletedLegendEntries.value, { ...entry }];
  legendEntries.value = legendEntries.value.filter((_, at) => at !== index);
  return true;
};

/**
 * Restores the deleted rows at `indexes`: the list becomes `rows`, which
 * holds them at their places.
 * @param {PaintSlice} slice
 * @param {ReadonlySet<number>} indexes Indexes into the deleted list.
 * @param {readonly LegendRow[]} rows
 */
export const restoreLegendRows = (slice, indexes, rows) => {
  const { deletedLegendEntries } = membersOf(slice);
  deletedLegendEntries.value = deletedLegendEntries.value.filter((_, at) => !indexes.has(at));
  replaceLegendRows(slice, rows);
};

/**
 * The renamed rows a draw left out (OV-120).
 * @param {PaintSlice} slice
 * @param {readonly KeptLegendRow[]} rows
 */
export const replaceDormantLegendRows = (slice, rows) => {
  membersOf(slice).dormantLegendEntries.value = keptCopy(rows);
};

/**
 * The captions of the rows the color rules add (a rule commit).
 * @param {PaintSlice} slice
 * @param {Iterable<string>} captions
 */
export const replaceAddedLegendCaptions = (slice, captions) => {
  membersOf(slice).addedLegendCaptions.value = new Set(captions);
};

/**
 * Renames the listed row `row` (the row object of the list): its caption
 * and its key record. A row the list does not hold is left alone.
 * @param {PaintSlice} slice
 * @param {LegendRow} row
 * @param {{ caption: LegendCaption, originalCaption: PythonLegendKey }} names
 */
export const renameListedLegendRow = (slice, row, { caption, originalCaption }) => {
  const entry = membersOf(slice).legendEntries.value.find((candidate) => candidate === row);
  if (!entry) return;
  entry.originalCaption = originalCaption;
  entry.caption = caption;
};

/**
 * The Legend color of the rows of `keys`: each key's color edit and the
 * recorded color of each listed row of them.
 * @param {PaintSlice} slice
 * @param {readonly PythonLegendKey[]} keys
 * @param {string} color The color the user set.
 */
export const recolorLegendRows = (slice, keys, color) => {
  const { legendEntries, legendColorOverrides } = membersOf(slice);
  const recolored = new Set(keys);
  // A color the user set is what the row records (`RecordedLegendColor`).
  const recorded = /** @type {RecordedLegendColor} */ (color);
  legendEntries.value.forEach((entry) => { if (recolored.has(legendEntryKey(entry))) entry.color = recorded; });
  keys.forEach((key) => { legendColorOverrides[key] = color; });
};

/**
 * Moves the records of re-keyed rows to their new keys (OV-660): the deleted
 * and dormant rows, the color and stroke edits, and the added captions, as
 * `rekeyLegendRecords` (services/legend-svg.js) gives them.
 * @param {PaintSlice} slice
 * @param {{ deletedEntries: readonly KeptLegendRow[], dormantEntries: readonly KeptLegendRow[],
 *   colorOverrides: Readonly<Record<string, string>>, strokeOverrides: DeepReadonly<Record<string, StrokeOverride>>,
 *   addedCaptions: readonly string[] }} moved
 */
export const rekeyLegendPaint = (slice, moved) => install(slice, {
  deletedLegendEntries: moved.deletedEntries, dormantLegendEntries: moved.dormantEntries,
  legendColorOverrides: moved.colorOverrides, legendStrokeOverrides: moved.strokeOverrides,
  addedLegendCaptions: moved.addedCaptions
});

/**
 * Retires the paint of the Legend rows of `keys`, whose data another owner's
 * edit removed (OV-87): their color and stroke edits, a listed row's rename
 * (the row shows its key again), and their dormant renames.
 * @param {PaintSlice} slice
 * @param {ReadonlySet<string>} keys
 * @returns {boolean} Whether a listed row's rename was retired.
 */
export const retireLegendPaint = (slice, keys) => {
  const { legendEntries, dormantLegendEntries, legendColorOverrides, legendStrokeOverrides } = membersOf(slice);
  keys.forEach((key) => {
    delete legendColorOverrides[key];
    delete legendStrokeOverrides[key];
  });
  /** @param {unknown} value */
  const text = (value) => String(value ?? '').trim();
  let renameRetired = false;
  const entries = legendEntries.value.map((entry) => {
    const original = text(entry.originalCaption);
    const caption = text(entry.caption);
    if (!keys.has(original) || !caption || caption === original) return entry;
    renameRetired = true;
    return { ...entry, caption: /** @type {LegendCaption} */ (original) };
  });
  const dormant = dormantLegendEntries.value.filter((entry) => !keys.has(text(entry.originalCaption)));
  if (dormant.length !== dormantLegendEntries.value.length) dormantLegendEntries.value = dormant;
  if (renameRetired) legendEntries.value = entries;
  return renameRetired;
};

/**
 * Settings Reset: the Legend and feature paint edits go, every renamed row
 * shows the caption Python draws, its key (R15-14), and the rows stay listed.
 * @param {PaintSlice} slice
 */
export const resetPaintEdits = (slice) => {
  const { legendEntries } = membersOf(slice);
  const rows = legendEntries.value.map((entry) => {
    const drawn = /** @type {LegendCaption} */ (/** @type {string} */ (entry.originalCaption));
    return entry.caption === drawn ? entry : { ...entry, caption: drawn };
  });
  install(slice, {
    deletedLegendEntries: [], dormantLegendEntries: [], addedLegendCaptions: [],
    legendColorOverrides: {}, legendStrokeOverrides: {}, featureColorOverrides: {}, featureStrokeOverrides: {}
  });
  legendEntries.value = rows;
};

/**
 * A new Session: no Legend list, and no feature paint edits.
 * @param {PaintSlice} slice
 */
export const clearSessionPaint = (slice) => {
  install(slice, { legendEntries: [], dormantLegendEntries: [], featureColorOverrides: {}, featureStrokeOverrides: {} });
};

// ---- Color and stroke edits ----

/**
 * A Legend row's color edit, by the row's key; null removes it.
 * @param {PaintSlice} slice
 * @param {PythonLegendKey} key
 * @param {string | null} color
 */
export const setLegendColor = (slice, key, color) => {
  const { legendColorOverrides } = membersOf(slice);
  if (color === null) delete legendColorOverrides[key];
  else legendColorOverrides[key] = color;
};

/**
 * A Legend row's stroke edit, by the row's key; null removes it.
 * @param {PaintSlice} slice
 * @param {PythonLegendKey} key
 * @param {Readonly<StrokeOverride> | null} stroke Copied.
 */
export const setLegendStroke = (slice, key, stroke) => {
  const { legendStrokeOverrides } = membersOf(slice);
  if (stroke) legendStrokeOverrides[key] = { ...stroke };
  else delete legendStrokeOverrides[key];
};

/**
 * A feature's color edit, by its feature key; null removes it.
 * @param {PaintSlice} slice
 * @param {string} key
 * @param {FeatureFill | null} fill Copied.
 */
export const setFeatureFill = (slice, key, fill) => {
  const { featureColorOverrides } = membersOf(slice);
  if (fill === null) delete featureColorOverrides[key];
  else featureColorOverrides[key] = typeof fill === 'string' ? fill : { ...fill };
};

/**
 * The feature color edits, whole (the rule commit's refresh, a
 * source-replacing Generate).
 * @param {PaintSlice} slice
 * @param {Readonly<Record<string, FeatureFill>>} fills
 */
export const replaceFeatureFills = (slice, fills) => install(slice, { featureColorOverrides: fills });

/**
 * A feature's stroke edit, by its feature key; null removes it.
 * @param {PaintSlice} slice
 * @param {string} key
 * @param {Readonly<StrokeOverride> | null} stroke Copied.
 */
export const setFeatureStroke = (slice, key, stroke) => {
  const { featureStrokeOverrides } = membersOf(slice);
  if (stroke) featureStrokeOverrides[key] = { ...stroke };
  else delete featureStrokeOverrides[key];
};

/**
 * Removes every Legend and feature stroke edit.
 * @param {PaintSlice} slice
 * @returns {boolean} Whether one was removed.
 */
export const clearStrokes = (slice) => {
  const { legendStrokeOverrides, featureStrokeOverrides } = membersOf(slice);
  const removed = Object.keys(legendStrokeOverrides).length > 0 || Object.keys(featureStrokeOverrides).length > 0;
  install(slice, { legendStrokeOverrides: {}, featureStrokeOverrides: {} });
  return removed;
};

export const PAINT_SLICE = Object.freeze({
  members: Object.freeze([
    'featureColorOverrides', 'featureStrokeOverrides', 'legendEntries', 'deletedLegendEntries', 'dormantLegendEntries',
    'legendColorOverrides', 'legendStrokeOverrides', 'addedLegendCaptions'
  ]),
  create, capture, install, reset
});
