// @ts-check
// The Web's move of a Session 46-48 to the drawings of a Session 49 (the
// 46-48 fast path, DESIGN §4.2; Python twin `drawings_document_from_mode_scoped`
// in gbdraw/session_drawings.py, pinned by
// tests/fixtures/sessions/drawings-vectors.json). A state-free module: every
// function reads only the document it is given.
//
// It owns the Session 46 Legend rekey (Session 47 keys a Legend's edits by
// row key): `legendEditsByRowKey` and `sessionLegendEditsByRowKey`, twins of
// `legend_edits_by_row_key` and `session_legend_edits_by_row_key`
// (gbdraw/session_io.py), pinned by
// tests/fixtures/sessions/legend-row-key-vectors.json.
import { ROW_KEYED_LEGEND_SESSION_VERSION } from '../contracts/schema.generated.js';
import { legendEntryKey } from './drawing/paint.js';
import { normalizeDefaultColor } from '../utils/color-utils.js';
import { normalizeLegendColorOverrides, normalizeStrokeOverrideMap } from './style-overrides.js';
/** @import { LegendEditsUnread } from './style-overrides.js' */

/** @param {unknown} value */
const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);

/** @typedef {'colorOverrides' | 'strokeOverrides'} LegendEditStore */
/**
 * One Session 46 Legend edit that did not stay as it was (`legendEditsByRowKey`):
 * `moved` to the key of the row whose caption it was stored under; kept under
 * its key as `unmatched` (no row shows or holds the key) or as `unshown` (the
 * key of a renamed row that no row shows: Session 46 applied it to that row
 * through Python's key); or a `collision`: a later edit of the store reached
 * the same row (`rowKey`) and Session 46 applied it over this one, dropped
 * (`toKey: null`) unless a stroke field of it is kept (`toKey`: `rowKey`). `slice`: the mode slice, or `draft` (`sessionLegendEditsByRowKey`).
 * @typedef {{ store: LegendEditStore, fromKey: string, toKey: string, outcome: 'moved' | 'unmatched' | 'unshown', slice?: string }
 *   | { store: LegendEditStore, fromKey: string, toKey: string | null, outcome: 'collision', rowKey: string, slice?: string }} LegendEditRekey
 */

// The rows a Session 46 Legend's edits reach (Session 46 Load,
// `normalizeSessionLegendEntries` and `normalizeEditorStateData` at
// 9839de8e: an entry with a caption and a usable color, the first of each
// caption in its list), by the caption each row shows (`draftLegendRows` at
// 9839de8e): a listed row before a dormant one (the later dormant entry of a
// caption wins), and a deleted row as its Restore shows it. A dormant entry
// that is not renamed, or whose key a listed or deleted row holds, is no row.
// `renamedKeys`: the keys of the renamed rows.
/** @param {Record<string, unknown>} legend */
const session46LegendRows = (legend) => {
  /** @param {string} list @returns {{ caption: string, key: string, dormant: boolean }[]} */
  const loaded = (list) => {
    const entries = legend[list];
    /** @type {Set<string>} */
    const captions = new Set();
    return (Array.isArray(entries) ? entries : []).flatMap((/** @type {unknown} */ entry) => {
      if (!isObject(entry)) return [];
      const row = /** @type {{ caption?: unknown, originalCaption?: unknown, color?: unknown, dormant?: unknown }} */ (entry);
      const caption = String(row.caption || '').trim();
      if (!caption || !normalizeDefaultColor(row.color) || captions.has(caption)) return [];
      captions.add(caption);
      return [{ caption, key: legendEntryKey(row), dormant: row.dormant === true }];
    });
  };
  const entries = loaded('entries');
  const listed = entries.filter(({ dormant }) => !dormant);
  const deleted = loaded('deletedEntries');
  const listedKeys = new Set(listed.map(({ key }) => key));
  const deletedKeys = new Set(deleted.map(({ key }) => key));
  const dormant = [...loaded('dormantEntries'), ...entries.filter((entry) => entry.dormant)]
    .filter(({ caption, key }) => caption !== key && !listedKeys.has(key) && !deletedKeys.has(key));
  /** @type {Map<string, string>} */
  const rowOf = new Map([...dormant, ...listed].map(({ caption, key }) => [caption, key]));
  deleted.forEach(({ caption, key }) => { if (!rowOf.has(caption)) rowOf.set(caption, key); });
  const renamedKeys = new Set([...listed, ...dormant, ...deleted].filter(({ caption, key }) => caption !== key).map(({ key }) => key));
  return { rowOf, renamedKeys };
};

// A Session 46 Legend (a slice's, or a flat draft's, `editorState.legend`)
// read as Session 47, and what the reading did to each edit that did not stay
// as it was. Session 46 stored the Legend color and stroke edits under the
// caption a row shows, Session 47 under the row key (`legendEntryKey`). Each
// row key gets the value Session 46 applied to that row (`compileDirectEditorMutationPlan`
// at 9839de8e, app/candidate-render.js): the stores as Session 46 Load read
// them (style-overrides.js), each edit reaching the row whose caption it is
// stored under, else the row of Python's key of that name; the color edits in
// their order, then the strokes in the order of the color edits' keys and then
// their own; a later edit painted over an earlier one, a stroke field by field
// (`strokeColor`, `strokeWidth`). Other fields are kept as they are. The
// outcomes are data for Load's repair report, computed before Load commits;
// the twin is `legend_edits_by_row_key` (gbdraw/session_io.py).
/**
 * @template {Record<string, unknown>} Legend
 * @param {Legend} legend
 * @param {LegendEditsUnread | null} [unread] Written: the edits the stores' reading drops (Load lists them).
 * @returns {{ legend: Legend, outcomes: LegendEditRekey[] }}
 */
export const legendEditsByRowKey = (legend, unread = null) => {
  const { rowOf, renamedKeys } = session46LegendRows(legend);
  /** @param {string} key */
  const rowKeyOf = (key) => rowOf.get(key) || key;
  const colors = normalizeLegendColorOverrides(legend.colorOverrides, unread?.colors);
  const strokes = normalizeStrokeOverrideMap(legend.strokeOverrides, { requireOverride: true, losses: unread?.strokes });
  /** @type {Record<string, unknown>} */
  const read = { ...legend };
  /** @type {LegendEditRekey[]} */
  const outcomes = [];
  /**
   * @template T
   * @param {LegendEditStore} store
   * @param {Record<string, T>} edits
   * @param {string[]} order
   * @param {(earlier: T | undefined, later: T) => T} paint
   * @param {(edit: T) => string[]} fieldsOf The fields an edit sets; a row keeps each from its last edit that sets it.
   */
  const keyByRow = (store, edits, order, paint, fieldsOf) => {
    if (!isObject(legend[store])) return;
    const lastOf = new Map(order.map((key) => [rowKeyOf(key), key]));
    /** @param {string} key @param {string} field */
    const fieldOf = (key, field) => JSON.stringify([rowKeyOf(key), field]);
    const fieldLastOf = new Map(order.flatMap((key) => fieldsOf(edits[key]).map((field) => [fieldOf(key, field), key])));
    /** @type {Record<string, T>} */
    const keyed = {};
    order.forEach((key) => {
      const rowKey = rowKeyOf(key);
      keyed[rowKey] = paint(keyed[rowKey], edits[key]);
      if (lastOf.get(rowKey) !== key) {
        const kept = fieldsOf(edits[key]).some((field) => fieldLastOf.get(fieldOf(key, field)) === key);
        outcomes.push({ store, fromKey: key, toKey: kept ? rowKey : null, outcome: 'collision', rowKey });
      } else if (rowKey !== key) outcomes.push({ store, fromKey: key, toKey: rowKey, outcome: 'moved' });
      else if (!rowOf.has(key)) outcomes.push({ store, fromKey: key, toKey: key, outcome: renamedKeys.has(key) ? 'unshown' : 'unmatched' });
    });
    read[store] = keyed;
  };
  const colorKeys = Object.keys(colors);
  keyByRow('colorOverrides', colors, colorKeys, (_earlier, later) => later, () => ['color']);
  keyByRow('strokeOverrides', strokes, [...new Set([...colorKeys, ...Object.keys(strokes)])].filter((key) => Object.hasOwn(strokes, key)),
    (earlier, later) => ({ ...earlier, ...later }), (edit) => Object.keys(edit));
  return { legend: /** @type {Legend} */ (read), outcomes };
};

// A Session document of `version` with each Legend's edits keyed by row key:
// every mode slice's (Session 46), and the flat draft's of a Session 27-44;
// with the outcomes of `legendEditsByRowKey`, each naming its slice. A
// Session 47 or later is returned as it is.
/**
 * @template {Record<string, unknown>} Session
 * @param {Session} session
 * @param {number} version
 * @param {{ slice: string, unread: LegendEditsUnread }[] | null} [unread] Written: what each slice's reading drops.
 * @returns {{ session: Session, outcomes: LegendEditRekey[] }}
 */
export const sessionLegendEditsByRowKey = (session, version, unread = null) => {
  /** @type {LegendEditRekey[]} */
  const outcomes = [];
  if (version >= ROW_KEYED_LEGEND_SESSION_VERSION) return { session, outcomes };
  /** @template {Record<string, unknown>} Holder @param {Holder} holder @param {string} slice @returns {Holder} */
  const withRowKeyedLegend = (holder, slice) => {
    const editor = holder.editorState;
    if (!isObject(editor) || !isObject(/** @type {Record<string, unknown>} */ (editor).legend)) return holder;
    const editorState = /** @type {Record<string, unknown>} */ (editor);
    /** @type {LegendEditsUnread} */
    const sliceUnread = { colors: [], strokes: { leaves: [], dropped: 0 } };
    const read = legendEditsByRowKey(/** @type {Record<string, unknown>} */ (editorState.legend), unread ? sliceUnread : null);
    read.outcomes.forEach((outcome) => outcomes.push({ ...outcome, slice }));
    unread?.push({ slice, unread: sliceUnread });
    return { ...holder, editorState: { ...editorState, legend: read.legend } };
  };
  const read = withRowKeyedLegend(session, 'draft');
  if (!isObject(read.modes)) return { session: read, outcomes };
  return {
    session: {
      ...read,
      modes: Object.fromEntries(Object.entries(/** @type {Record<string, unknown>} */ (read.modes)).map(([mode, slice]) => [
        mode, isObject(slice) ? withRowKeyedLegend(/** @type {Record<string, unknown>} */ (slice), mode) : slice
      ]))
    },
    outcomes
  };
};
