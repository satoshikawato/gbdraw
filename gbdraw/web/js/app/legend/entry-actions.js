// @ts-check
/** @import { DrawingState } from '../../state.js' */
/** @import { LegendRowName } from '../../services/legend-row-names.js' */
/** @import { PythonLegendKey } from '../../services/legend-svg.js' */
/** @import { RuleRowColor } from '../../services/specific-color-rules.js' */
/**
 * A Legend row record as the editor lists, restores, and compares it.
 * @typedef {{ readonly caption?: string, readonly originalCaption?: string, readonly color?: string,
 *   readonly featureIds?: readonly unknown[], readonly xPos?: number, readonly yPos?: number }} EntryRecord
 */
/** @typedef {EntryRecord & { caption: string, originalCaption: string }} ListedEntry A row a Result lists, under its caption and its key. */
/**
 * Where a rule commit puts a row it splits off: right after the row of key
 * `after` (D-29), or right before the row of key `before` (OV-2194).
 * @typedef {{ caption: PythonLegendKey, after: PythonLegendKey }} LegendRowPlacement
 */
/**
 * The rows a Result lists (`listLegendRows`) and the keys of those that are Python's.
 * @typedef {{ entries: ListedEntry[], generatedCaptions: Set<string> }} LegendRowListing
 */
import { normalizeOptionalHexColor, resolveColorToHex, toNativeColorInputValue } from '../../utils/color-utils.js';
import { diagnosticError } from '../../utils/error-normalization.js';
import { keepBothNames, legendSwatch } from '../../services/legend-row-names.js';
import {
  defaultLegendKeyOrder,
  drawnBlockStroke,
  getAllFeatureLegendGroups,
  isLegendOrderEdited,
  legendRowKeys,
  recordedLegendOrder,
  rekeyLegendRecords,
  resultLegendRows,
  shownLegendRowNames
} from '../../services/legend-svg.js';
import {
  deleteLegendRow,
  legendEntryKey,
  recolorLegendRows,
  rekeyLegendPaint,
  replaceDormantLegendRows,
  replaceLegendRows,
  restoreLegendRows,
  setLegendColor
} from '../../services/drawing/paint.js';
/** @import { KeptLegendRow, LegendRow } from '../../services/drawing/paint.js' */
import { diffLegendIntents, legendColorCopies, legendRowRules } from '../../services/specific-color-rules.js';
import { displayedLegendRowContext } from '../../services/feature-visibility.js';

/** @param {unknown} value */
const normalizedColor = (value) => {
  const resolved = String(resolveColorToHex(String(value || '').trim()) || value || '').trim().toLowerCase();
  return resolved.startsWith('#') ? toNativeColorInputValue(resolved) : resolved;
};

/** @param {EntryRecord | null | undefined} entry */
const legendCaption = (entry) => String(entry?.caption || '').trim();
/** @param {EntryRecord | null | undefined} entry */
const generatedCaption = (entry) => String(entry?.originalCaption || entry?.caption || '').trim();
// The rows this owner builds from a Result's rows (`listLegendRows`), the
// rows of the drawing's lists, and rows made from them are Legend rows (the
// casts: `EntryRecord` is their shape as this owner reads it).
/** @param {readonly EntryRecord[]} rows @returns {readonly LegendRow[]} */
const listedRowsOf = (rows) => /** @type {readonly LegendRow[]} */ (rows);
/** @param {readonly EntryRecord[]} rows @returns {readonly KeptLegendRow[]} */
const keptRowsOf = (rows) => /** @type {readonly KeptLegendRow[]} */ (rows);

// The entries in an order of row keys: entries the order does not list follow
// in their order; with `keepFollowed`, entries whose listed ones already
// follow the order keep it (B18), as `orderLegendEntries` orders a Legend group.
/**
 * @param {EntryRecord[]} entries
 * @param {readonly string[]} keys
 * @param {{ keepFollowed?: boolean }} [options]
 */
const entriesInOrder = (entries, keys, { keepFollowed = false } = {}) => {
  /** @type {Map<string, number>} */
  const rank = new Map();
  keys.forEach((key) => { if (key && !rank.has(key)) rank.set(key, rank.size); });
  const listed = entries.map((entry) => rank.get(legendEntryKey(entry))).filter((value) => value !== undefined);
  if (keepFollowed && listed.every((value, index) => index === 0 || listed[index - 1] < value)) return entries;
  /** @param {EntryRecord} entry */
  const at = (entry) => rank.get(legendEntryKey(entry)) ?? Infinity;
  return entries.map((entry, index) => ({ entry, index }))
    .sort((left, right) => (at(left.entry) - at(right.entry)) || (left.index - right.index))
    .map(({ entry }) => entry);
};

// A row of Python's that the Result does not show has no place on it: it goes
// before the first placed row that follows it in the Result's inventory, else
// last, so a returning row takes the place Python gives it (OV-154, OV-1198).
/**
 * @template T
 * @param {readonly T[]} placed In order.
 * @param {readonly T[]} unplaced In inventory order.
 * @param {readonly string[]} inventory
 * @param {(item: T) => string} keyOf
 * @returns {T[]}
 */
const placeByInventory = (placed, unplaced, inventory, keyOf) => {
  const order = [...placed];
  unplaced.forEach((item) => {
    const keys = order.map(keyOf);
    const next = inventory.slice(inventory.indexOf(keyOf(item)) + 1).find((later) => keys.includes(later));
    order.splice(next === undefined ? order.length : keys.indexOf(next), 0, item);
  });
  return order;
};

/**
 * @typedef {object} LegendEntryActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(() => string | undefined) | null} [readActiveResultIdentity]
 *   The preview owner's runtime identity of the mounted Result.
 * @property {import('../candidate-render.js').ReadLegendRowColorsPort} readLegendRowColors
 *   The colors a drawing's listed and deleted rows show (`displayedLegendRowColors`, OV-282), read
 *   from the composition root's derivation the Legend panel shows; never the color a row records.
 * @property {{ value: ReadonlyMap<PythonLegendKey, import('../../services/legend-svg.js').ShownLegendRowNames> }} [legendRowNames]
 *   The root's reactive holder of the names the displayed Result shows for its rows;
 *   this owner is its one writer.
 * @property {(drawing: DrawingState) => readonly string[]} [readHiddenTrackRows]
 *   The keys of the displayed Result's rows whose drawn track the drawing hides (OV-1198, OV-1215),
 *   which the Legend list leaves out as the live compile deletes them.
 */

// OV-2194, OV-2247, OV-2248 (Owner 2026-10-29): `items` (listed rows, or the
// keys of an inventory) where Generate draws them from the next rule table,
// by Python's answer (`legendOrder`, the keys of its feature rows in Legend
// order). An added item goes right before the first item the answer puts
// after it, else right after the last one it puts before it; one with neither
// stays out (`placed` names those placed). With `reorder` (the default order,
// which follows Python's), the items of the answer also take their places in
// its order, among the places they hold (an old row Python moves); the other
// items keep theirs. In an edited order, which Generate replays (PD-OI-063),
// the old rows stay.
/**
 * @template T
 * @param {readonly T[]} items
 * @param {ReadonlyMap<string, T>} added The items the commit adds, by key.
 * @param {readonly PythonLegendKey[]} legendOrder
 * @param {(item: T) => string} keyOf
 * @param {boolean} reorder
 * @returns {{ items: T[], placed: PythonLegendKey[] }}
 */
const inAnsweredLegendOrder = (items, added, legendOrder, keyOf, reorder) => {
  /** @type {Map<string, number>} */
  const rank = new Map();
  legendOrder.forEach((key, index) => { if (!rank.has(key)) rank.set(key, index); });
  /** @param {T} item */
  const rankOf = (item) => rank.get(keyOf(item)) ?? -1;
  let rows = [...items];
  /** @type {PythonLegendKey[]} */
  const placed = [];
  for (const [key, row] of [...added].filter(([key]) => rank.has(key)).sort(([a], [b]) => (rank.get(a) ?? 0) - (rank.get(b) ?? 0))) {
    const at = rank.get(key) ?? 0;
    const next = rows.findIndex((item) => rankOf(item) > at);
    const previous = rows.findLastIndex((item) => rankOf(item) >= 0 && rankOf(item) < at);
    const index = next >= 0 ? next : previous >= 0 ? previous + 1 : -1;
    if (index < 0) continue;
    rows = [...rows.slice(0, index), row, ...rows.slice(index)];
    placed.push(/** @type {PythonLegendKey} */ (key));
  }
  if (!reorder) return { items: rows, placed };
  const slots = rows.flatMap((item, index) => (rankOf(item) >= 0 ? [index] : []));
  const ordered = slots.map((index) => rows[index]).sort((a, b) => rankOf(a) - rankOf(b));
  slots.forEach((slot, index) => { rows[slot] = ordered[index]; });
  return { items: rows, placed };
};

// The Legend entry owner (U3a, R1): its writers edit the drawing's Legend
// intent (`legendEntries`, `deletedLegendEntries`, the override stores) and
// write no row; the composition root shows the intent on the displayed Result
// through the port (`showEditorIntent`), whose executor writes the rows.
/** @param {LegendEntryActionsOptions} options */
export const createLegendEntryActions = ({
  state,
  readActiveResultIdentity = null,
  readLegendRowColors,
  legendRowNames = { value: new Map() },
  readHiddenTrackRows = () => []
}) => {
  const {
    results,
    svgContainer,
    originalLegendOrder,
    originalLegendColors
  } = state;

  /** @type {((options?: { commit?: boolean }) => unknown) | null} */
  let legendGeometryChangedHandler = null;

  // The generated inventory of Legend captions in Python's order, per Result
  // (OV-47): each Result has its own order, so the displayed Result's order is
  // `originalLegendOrder`, and the others wait here under their runtime
  // identity until they are displayed. `replayedInventoryResults` are the
  // Results a draw replayed an edited order onto: their drawing no longer shows
  // Python's order.
  /** @type {Map<string, string[]>} */
  const inventoryByResult = new Map();
  /** @type {Set<string>} */
  const replayedInventoryResults = new Set();

  // The names the displayed Result shows for its Legend rows, by each row's
  // key: read when the rows are listed and after the executor named them. A
  // row's key is never shown; every surface shows `shownLegendName` of it.
  const readLegendRowNames = () => {
    legendRowNames.value = shownLegendRowNames(svgContainer.value?.querySelector?.('svg') || null);
  };
  // The name a caption stands for (D-15): the name of the row it is the key
  // of, as the displayed Result shows it, else as the color rules that draw the
  // row are captioned (a deleted or dormant row the Result does not hold),
  // else the caption itself (a row whose key is its name, or a typed name).
  // Every write of a caption that was read as a row takes this name (OV-640).
  /** @param {unknown} caption @returns {LegendRowName} */
  const legendRowName = (caption) => {
    const key = /** @type {PythonLegendKey} */ (String(caption ?? ''));
    return /** @type {LegendRowName} */ (legendRowNames.value.get(key)?.name
      ?? legendRowRules(key, displayedLegendRowContext(state, state.activeDrawing()))[0]?.cap
      ?? key);
  };
  /** @param {unknown} caption A row's key. @returns {string} The name the Legend draws for it. */
  const shownLegendName = (caption) => legendRowNames.value.get(/** @type {PythonLegendKey} */ (String(caption ?? '')))?.drawnName
    ?? legendRowName(caption);
  // The name of a row of the drawing's lists: the displayed Result's name of
  // its key, else its caption where a rename gave it one, else the name of
  // its key (a deleted editor row, or a row the Result does not hold).
  /** @param {EntryRecord} entry @returns {LegendRowName} */
  const legendEntryName = (entry) => {
    const key = legendEntryKey(entry);
    const caption = legendCaption(entry);
    return legendRowNames.value.get(key)?.name ?? (caption && caption !== key ? /** @type {LegendRowName} */ (caption) : legendRowName(key));
  };
  // Whether the Legend tab leaves out a row of the drawing's lists: the
  // displayed Result shows it as one row with an earlier row of its name and
  // swatch (D-28). Its entry stays, so its rename and edits stay.
  /** @param {EntryRecord} entry */
  const legendEntryMerged = (entry) => legendRowNames.value.get(legendEntryKey(entry))?.merged === true;
  // The name the Legend tab shows for a row of the drawing's lists.
  /** @param {EntryRecord} entry @returns {string} */
  const shownLegendEntryName = (entry) => legendRowNames.value.get(legendEntryKey(entry))?.drawnName ?? legendEntryName(entry);
  // Whether `typed` names the row of the key `caption`: its key, its name, or
  // its drawn name. Typing one of them into the row's name field is no rename.
  /** @param {unknown} caption @param {unknown} typed */
  const legendRowIsNamed = (caption, typed) => {
    const key = /** @type {PythonLegendKey} */ (String(caption ?? '').trim());
    const text = String(typed ?? '').trim();
    const names = legendRowNames.value.get(key);
    return text === key || (Boolean(names) && (text === names?.name || text === names?.drawnName));
  };
  // PD-OI-061 scenario revision 5: the row a typed name names, by its key. A
  // listed row the Legend draws whose drawn name is the text, exactly, then
  // ignoring case (a shown numbered name such as `CDS (2)` names that row);
  // else a listed or deleted row whose name is the text ignoring case, the one
  // in `color` first, else the first in legend order. Null when no row has
  // the name (a color rule's name that no row shows is the rule owner's).
  /** @param {unknown} typed @param {unknown} color @returns {PythonLegendKey | null} */
  const findLegendRowNamed = (typed, color) => {
    const text = String(typed ?? '').trim();
    if (!text) return null;
    const lower = text.toLowerCase();
    const drawing = state.activeDrawing();
    /** @type {EntryRecord[]} */
    const listed = drawing.legendEntries?.value || [];
    /** @type {EntryRecord[]} */
    const deleted = drawing.deletedLegendEntries?.value || [];
    const drawn = listed.flatMap((entry) => {
      const names = legendRowNames.value.get(legendEntryKey(entry));
      return names?.shown ? [{ key: legendEntryKey(entry), drawnName: String(names.drawnName) }] : [];
    });
    const byDrawn = drawn.find((row) => row.drawnName === text) || drawn.find((row) => row.drawnName.toLowerCase() === lower);
    if (byDrawn) return byDrawn.key;
    const wanted = normalizedColor(color);
    const shown = readLegendRowColors(drawing);
    const named = [
      ...listed.map((entry) => ({ entry, color: shown.listed(entry) })),
      ...deleted.map((entry) => ({ entry, color: shown.deleted(entry) }))
    ].filter(({ entry }) => String(legendEntryName(entry)).toLowerCase() === lower);
    const found = named.find((row) => normalizedColor(row.color) === wanted) || named[0];
    return found ? legendEntryKey(found.entry) : null;
  };

  // PD-OI-061 scenario revision 5: the names "Keep both" gives, for the rows
  // the displayed Result draws. The row renamed whole keeps its place and its
  // swatch; features split off their row (`renamedKey` null) are a row placed
  // right after the target row (D-29), drawn in `color` with the Result's
  // block stroke, as Python draws a color rule's row. A deleted target
  // (`deletedKey`) is restored first when the two would be one row (D-28):
  // the swatch its hidden row holds is the renamed row's (review L-1).
  /** @type {import('../feature-editor/color-actions.js').KeepBothNamesPort} */
  const keepBothNamesOf = ({ name, renamedKey, targetKey, color, deletedKey = null }) => {
    const rows = [...legendRowNames.value].filter(([, names]) => names.shown)
      .map(([key, names]) => ({ key, name: names.name, swatch: names.swatch }));
    const renamed = rows.find((row) => row.key === renamedKey);
    const block = renamed ? null : drawnBlockStroke(svgContainer.value?.querySelector?.('svg') || null);
    const swatch = renamed?.swatch ?? legendSwatch({ fill: color, stroke: block?.color, width: block?.width });
    const names = keepBothNames(rows, { name: /** @type {LegendRowName} */ (name), swatch, renamedKey: renamed ? renamed.key : null, targetKey });
    return { ...names, restores: deletedKey !== null && legendRowNames.value.get(deletedKey)?.swatch === swatch };
  };

  /** @param {(options?: { commit?: boolean }) => unknown} handler The root's refresh of the Legend geometry. */
  const setLegendGeometryChangedHandler = (handler) => {
    legendGeometryChangedHandler = typeof handler === 'function' ? handler : null;
  };

  /** @param {{ commit?: boolean }} [options] */
  const onLegendGeometryChanged = (options) => legendGeometryChangedHandler?.(options);

  // B19 (D-07, D-08): a History step made on another batch Result. Its
  // restored list describes that Result, so the displayed Result's list is
  // written from the rows it shows and the step's shared Legend intent, from
  // its other side `from` to the restored list: captions, colors, deletions,
  // returning rows, and the order. An entry only another Result draws is never
  // listed here and an entry only this Result draws stays; the port then shows
  // the list. Outside a batch, or when `from` lists the rows this Result shows,
  // the restored list describes this Result and stays. Returns whether the
  // list was rewritten.
  /** @param {{ from?: EntryRecord[] | null }} [options] */
  const adoptRestoredLegend = ({ from = null } = {}) => {
    const drawing = state.activeDrawing();
    const svg = svgContainer.value?.querySelector?.('svg');
    /** @type {EntryRecord[]} */
    const restored = drawing.legendEntries.value || [];
    if (!svg || results.value.length < 2 || !Array.isArray(from)) return false;
    const mounted = listLegendRows(svg, { entries: restored, dormant: drawing.dormantLegendEntries.value || [], asShown: true });
    if (!mounted) return false;
    // A row is its key and the name it shows: two rows may show one name.
    /** @param {EntryRecord} entry */
    const shown = (entry) => `${legendEntryKey(entry)}\u0000${legendCaption(entry)}`;
    const fromRows = new Set(from.map(shown));
    if (mounted.entries.length === fromRows.size && mounted.entries.every((entry) => fromRows.has(shown(entry)))) {
      return false;
    }
    /** @param {EntryRecord[]} list */
    const byKey = (list) => new Map(list.filter(legendCaption).map((entry) => [/** @type {string} */ (legendEntryKey(entry)), entry]));
    const before = byKey(from);
    const after = byKey(restored);
    // A row no Result's inventory lists is an editor row; a returning row is
    // this Result's own or an editor row.
    const inventories = [...inventoryByResult.values(), originalLegendOrder.value];
    /** @param {EntryRecord} entry */
    const editorRow = (entry) => !inventories.some((inventory) => inventory.includes(generatedCaption(entry)));
    const mountedKeys = new Set(mounted.entries.map(legendEntryKey));
    const returning = [...after.values()].filter((entry) => (
      !before.has(legendEntryKey(entry)) && !mountedKeys.has(legendEntryKey(entry))
      && (originalLegendOrder.value.includes(generatedCaption(entry)) || editorRow(entry))
    ));
    /** @param {EntryRecord} entry @param {EntryRecord} next */
    const stepColor = (entry, next) => {
      const previous = before.get(legendEntryKey(next));
      return !previous || normalizedColor(previous.color) !== normalizedColor(next.color) ? next.color : entry.color;
    };
    /** @type {EntryRecord[]} */
    let entries = mounted.entries.flatMap((entry) => {
      const next = after.get(legendEntryKey(entry));
      if (!next) return before.has(legendEntryKey(entry)) ? [] : [entry];
      return [{ ...entry, caption: legendCaption(next), color: stepColor(entry, next) }];
    });
    entries = [...entries, ...returning];
    // A Result already in the order keeps its own entries' places (B18); a
    // returned entry takes its place in the order, and own entries follow. A
    // step that leaves the default order of the Result it was made on gives
    // this Result its own default order (OV-47).
    const orderChanged = [...before.keys()].join('\u0000') !== [...after.keys()].join('\u0000');
    const generatedEntries = restored.filter((entry) => !editorRow(entry));
    const restoresDefault = [...inventoryByResult.values()].some((inventory) => (
      generatedEntries.every((entry) => inventory.includes(generatedCaption(entry)))
      && !isLegendOrderEdited(restored, inventory)
    ));
    if (orderChanged || returning.length > 0) {
      const order = originalLegendOrder.value.length > 0 && restoresDefault
        ? defaultLegendKeyOrder(entries, [...originalLegendOrder.value])
        : [...after.keys()];
      entries = entriesInOrder(entries, order, { keepFollowed: returning.length === 0 });
    }
    replaceLegendRows(drawing.slices.paint, listedRowsOf(entries));
    return true;
  };

  // The displayed Result's inventory (`originalLegendOrder`, Python's order
  // of its rows) and its stored copy, which `adoptResultInventory` adopts
  // when the Result is displayed again.
  /** @param {readonly string[]} inventory */
  const writeDisplayedInventory = (inventory) => {
    originalLegendOrder.value = [...inventory];
    const identity = readActiveResultIdentity?.() || '';
    if (identity) inventoryByResult.set(identity, [...inventory]);
  };

  // The Legend rows the specific-color rules draw, computed on the listed rows
  // (R1, R13): the rows the rules own (a caption and color of the previous rule
  // intents) and the listed rows of a wanted caption and color are diffed with
  // the wanted rows; a listed row the rules do not own with a wanted caption
  // in another color stops the commit. The rule owner runs its transition and
  // then `apply` in one History step, while `isCurrent` holds; `apply` writes
  // the rows into the intent, a row the commit splits off right after the row
  // `placement.after` (D-29), the rows of `legendOrder` where it puts them
  // (`inAnsweredLegendOrder`), and returns the rows the displayed Result shows
  // at once until Python draws them (gaps 1 and 2; null when the rows only
  // change color). A row of `rekeyed` (old key to new key) is one row whose
  // key the commit changes: a rule row the rule table change re-allocates
  // (OV-645), or a whole row a rename gives a new key (OV-158). The row the
  // commit adds under the new key takes the place of the row it removes under
  // the old one, with every record stored under the old key, in the list and
  // in the generated inventory, so the order stays the default order or the
  // user's. Resolves to false when the Result or the caller became stale.
  /**
   * @param {(EntryRecord & { name?: string })[]} intents
   * @param {{ drawing?: DrawingState, previousFileIntents?: EntryRecord[], isCurrent?: () => boolean,
   *   placement?: LegendRowPlacement | null, rekeyed?: ReadonlyMap<PythonLegendKey, PythonLegendKey>,
   *   legendOrder?: readonly PythonLegendKey[] | null }} [options]
   */
  const prepareFileLegendEntries = async (intents, {
    drawing = state.activeDrawing(), previousFileIntents = [], isCurrent = () => true, placement = null, rekeyed = new Map(),
    legendOrder = null
  } = {}) => {
    const mountedSvg = svgContainer.value?.querySelector('svg');
    if (getAllFeatureLegendGroups(mountedSvg).length === 0) {
      return { diff: { add: [], update: [], remove: [], unchanged: [] }, reorders: false, isCurrent: () => true, apply: () => null };
    }
    /** @type {Map<string, Set<string>>} */
    const provenance = new Map();
    for (const entry of previousFileIntents) {
      const caption = legendCaption(entry);
      if (!provenance.has(caption)) provenance.set(caption, new Set());
      provenance.get(caption)?.add(normalizedColor(entry?.color));
    }
    // A listed row is its key: an intent's `caption` is its row's key, and
    // its color the one it shows (OV-282: the palette's for a palette row),
    // as the previous intents name it.
    const shown = readLegendRowColors(drawing);
    /** @type {{ caption: string, color: string }[]} */
    const listed = (drawing.legendEntries.value || [])
      .map((/** @type {EntryRecord} */ entry) => ({ caption: legendEntryKey(entry), color: normalizedColor(shown.listed(entry)) }));
    const owned = listed.filter((entry) => provenance.get(entry.caption)?.has(entry.color));
    const desiredByCaption = new Map(intents.map((intent) => [intent.caption, normalizedColor(intent.color)]));
    for (const intent of intents) {
      const existing = listed.find((entry) => entry.caption === intent.caption);
      if (existing && !owned.includes(existing) && existing.color !== normalizedColor(intent.color)) {
        // OV-644: a typed refusal; its context names no row key.
        throw diagnosticError('LEGEND_ROW_CONFLICT', {}, { stage: 'rule-validation', operation: 'evaluateRules' });
      }
    }
    const reusable = listed.filter((entry) => !owned.includes(entry) && desiredByCaption.get(entry.caption) === entry.color);
    const diff = diffLegendIntents([...owned, ...reusable], intents);
    if (!isCurrent() || svgContainer.value?.querySelector('svg') !== mountedSvg) return false;
    const removedKeys = new Set(diff.remove.map(({ caption }) => caption));
    const inDefaultOrder = () => !isLegendOrderEdited(drawing.legendEntries.value || [], originalLegendOrder.value);
    /** @param {EntryRecord} entry */
    const keyOf = (entry) => legendEntryKey(entry);
    // Whether Python's answer, in the default order, moves a listed row the
    // commit keeps: the displayed Result shows the answer's order at once
    // (`order`), which no later live compile keeps, so the commit asks for
    // the automatic rerender (`reorders`).
    /** @param {readonly EntryRecord[]} kept @param {readonly EntryRecord[]} answered */
    const movesKeptRow = (kept, answered) => answered.filter((entry) => kept.includes(entry)).some((entry, index) => entry !== kept[index]);
    const listedKept = (drawing.legendEntries.value || []).filter((/** @type {EntryRecord} */ entry) => !removedKeys.has(keyOf(entry)));
    const reorders = Boolean(legendOrder) && inDefaultOrder()
      && movesKeptRow(listedKept, inAnsweredLegendOrder(listedKept, new Map(), legendOrder || [], keyOf, true).items);
    return {
      diff,
      reorders,
      isCurrent: () => svgContainer.value?.querySelector('svg') === mountedSvg,
      apply: () => {
        const inDefault = inDefaultOrder();
        const removed = removedKeys;
        // An intent's `caption` is its row's key; its `name` is the rule's caption.
        const names = new Map(intents.map((intent) => [legendCaption(intent), intent.name]));
        // A row the rules draw records their color, also one that already showed it (a palette row).
        const colors = new Map([...diff.update, ...diff.unchanged].map(({ caption, color }) => [caption, color]));
        const added = new Map(diff.add.map(({ caption, color }) => [caption, { caption, originalCaption: caption, color, featureIds: [] }]));
        /** @type {EntryRecord[]} */
        let entries = (drawing.legendEntries.value || []).map((/** @type {EntryRecord} */ entry) => (
          colors.has(legendEntryKey(entry)) ? { ...entry, color: colors.get(legendEntryKey(entry)) } : entry
        ));
        // A re-keyed row (OV-645) is the row the commit removes, added under
        // its new key: the row and every record stored under its key move to
        // the new key in place (OV-660), so it keeps its place, its rank in
        // the inventory, and its edits; the Result shows it before the row it
        // replaces.
        /** @type {Map<PythonLegendKey, PythonLegendKey>} */
        const moves = new Map([...rekeyed].filter(([from, to]) => removed.has(from) && added.has(to)));
        /** @type {Map<string, string>} A row the commit adds, and the row it is shown before. */
        const before = new Map([...moves].map(([from, to]) => [to, from]));
        if (moves.size > 0) {
          // The moved row is the added one: it shows the name of its rules.
          const movedTo = new Set(moves.values());
          entries = moveLegendRecords(drawing, entries, moves).map((entry) => (
            movedTo.has(legendEntryKey(entry)) ? { ...entry, caption: legendEntryKey(entry) } : entry
          ));
          moves.forEach((to) => added.delete(to));
        }
        // D-29: a row a commit splits off a row goes right after the row of
        // its name, so the editor list records the edited order Generate
        // replays (PD-OI-063), and the Result shows it there.
        /** @type {PythonLegendKey[]} The added rows placed in the list. */
        const placed = [];
        const split = placement && added.get(placement.caption);
        const target = split ? entries.findIndex((entry) => legendEntryKey(entry) === placement.after) : -1;
        if (placement && split && target >= 0) {
          entries = [...entries.slice(0, target + 1), split, ...entries.slice(target + 1)];
          added.delete(placement.caption);
          placed.push(placement.caption);
        }
        entries = entries.filter((entry) => !removed.has(legendEntryKey(entry)));
        /** @type {PythonLegendKey[] | null} */
        let order = null;
        if (legendOrder) {
          const kept = entries;
          const answered = inAnsweredLegendOrder(entries, added, legendOrder, keyOf, inDefault);
          entries = answered.items;
          answered.placed.forEach((key) => { added.delete(key); placed.push(key); });
          // Rows in Python's order are no edit (D-08): in the default order,
          // the inventory takes the answer's order and the rows placed, so
          // the list stays the default order, Generate replays no order on
          // another Result, and a Session saves none.
          if (inDefault) {
            const inventory = originalLegendOrder.value;
            const keys = new Map(answered.placed.filter((key) => !inventory.includes(key)).map((key) => [key, key]));
            writeDisplayedInventory(inAnsweredLegendOrder(inventory, keys, legendOrder, String, true).items);
            if (movesKeptRow(kept, entries)) order = [...entries, ...added.values()].map(keyOf);
          }
        }
        for (const key of placed) {
          const next = entries[entries.findIndex((entry) => legendEntryKey(entry) === key) + 1];
          if (next) before.set(key, legendEntryKey(next));
        }
        replaceLegendRows(drawing.slices.paint, listedRowsOf([...entries, ...added.values()]));
        const rows = {
          add: diff.add.map(({ caption, color }) => ({
            caption, color, ...(names.get(caption) ? { name: names.get(caption) } : {}),
            ...(before.has(caption) ? { before: before.get(caption) } : {})
          })),
          retire: [...removed],
          ...(order ? { order } : {})
        };
        return rows.add.length > 0 || rows.retire.length > 0 || order ? rows : null;
      }
    };
  };

  // Moves the drawing's Legend records, the generated inventory and colors,
  // and the open Stroke options of each re-keyed row to its new key
  // (`rekeyLegendRecords`, OV-660); returns `entries`, the list being
  // written, moved. For a drawing that is not shown (`hidden`: the Specific
  // rules file card's Remove, OV-1280), the records of the `removed` rows go
  // too, with the Legend colors that copy them (`legendColorCopies`), the
  // moved list is written here, and the inventory, colors, and Stroke
  // options, which are the shown Result's, stay.
  /**
   * @param {DrawingState} drawing
   * @param {readonly EntryRecord[]} entries
   * @param {ReadonlyMap<PythonLegendKey, PythonLegendKey>} moves
   * @param {{ removed: readonly RuleRowColor[] } | null} [hidden] The removed rule rows.
   * @returns {EntryRecord[]}
   */
  const moveLegendRecords = (drawing, entries, moves, hidden = null) => {
    /** @type {Set<string>} */
    const removed = new Set(hidden ? hidden.removed.map(({ caption }) => caption) : []);
    /** @param {readonly EntryRecord[]} list @returns {EntryRecord[]} */
    const kept = (list) => (hidden ? list.filter((entry) => !removed.has(legendEntryKey(entry))) : [...list]);
    const colorOverrides = { ...drawing.legendColorOverrides };
    if (hidden) legendColorCopies(colorOverrides, hidden.removed).forEach((key) => { delete colorOverrides[key]; });
    /** @type {Set<string>} */
    const strokeOptionsOpen = state.legendStrokeOptionsOpen || new Set();
    const moved = rekeyLegendRecords({
      entries: kept(entries),
      deletedEntries: kept(drawing.deletedLegendEntries.value || []),
      dormantEntries: kept(drawing.dormantLegendEntries.value || []),
      colorOverrides,
      strokeOverrides: { ...drawing.legendStrokeOverrides },
      addedCaptions: [...drawing.addedLegendCaptions.value].filter((caption) => !removed.has(caption)),
      originalOrder: hidden ? [] : /** @type {string[]} */ (originalLegendOrder.value),
      originalColors: hidden ? {} : /** @type {Record<string, string>} */ (originalLegendColors.value),
      strokeOptionsOpen: hidden ? [] : [...strokeOptionsOpen]
    }, moves);
    rekeyLegendPaint(drawing.slices.paint, {
      ...moved, deletedEntries: keptRowsOf(moved.deletedEntries), dormantEntries: keptRowsOf(moved.dormantEntries)
    });
    if (hidden) {
      replaceLegendRows(drawing.slices.paint, listedRowsOf(moved.entries));
      return moved.entries;
    }
    originalLegendOrder.value = moved.originalOrder;
    originalLegendColors.value = moved.originalColors;
    strokeOptionsOpen.clear();
    moved.strokeOptionsOpen.forEach((caption) => strokeOptionsOpen.add(caption));
    return moved.entries;
  };

  // The drawing's Legend intent that names a Result's rows.
  /** @param {DrawingState} drawing */
  const legendIntentOf = (drawing) => ({
    entries: drawing.legendEntries.value || [],
    dormant: drawing.dormantLegendEntries.value || [],
    deleted: drawing.deletedLegendEntries.value || [],
    inventory: originalLegendOrder.value || [],
    hidden: readHiddenTrackRows(drawing)
  });

  /**
   * The drawing's Legend list on the rows of a Result (U3b: one direction,
   * intent to shown rows). A row of Python's is Python's key: the record a
   * rename keeps, else the key of the intent entry that names the row (a
   * Result saved before U3a has no record). It is listed unless the intent
   * deletes that key, named by the intent's rename of it (else a dormant one,
   * OV-120), and keeps that entry's feature ids. An editor row is listed when
   * the intent lists it. The Result gives the rows' order and the paint the
   * executor showed from the intent; a row of Python's it does not show takes
   * its place in the inventory (`placeByInventory`). `asShown` lists the rows the Result shows
   * under the keys they show (B19 reads the list a Result showed). Also the
   * generated captions: the keys of the listed rows of Python's. A row of a
   * drawn track the intent hides (`hidden`) is not listed, as the live
   * compile deletes it (OV-1198).
   * @param {SVGSVGElement} svg
   * @param {{ entries: readonly EntryRecord[], dormant?: readonly EntryRecord[], deleted?: readonly EntryRecord[],
   *   inventory?: string[], hidden?: readonly string[], asShown?: boolean }} intent
   */
  const listLegendRows = (svg, { entries, dormant = [], deleted = [], inventory = [], hidden = [], asShown = false }) => {
    const rows = resultLegendRows(svg);
    if (!rows) return null;
    const known = new Set(inventory);
    const deletedKeys = new Set(deleted.map(generatedCaption).filter((key) => known.has(key)));
    const hiddenKeys = new Set(hidden);
    /** @param {(entry: EntryRecord) => boolean} test */
    const intentEntry = (test) => entries.find(test) || dormant.find(test);
    // The compile renames a key the inventory lists (else the entry is an
    // editor row), and a dormant row wherever it is drawn again unless a
    // listed entry names that key.
    /** @param {string} key */
    const renameOf = (key) => {
      const keyed = entries.filter((entry) => generatedCaption(entry) === key);
      if (keyed.length === 0) return dormant.find((entry) => generatedCaption(entry) === key);
      return known.has(key) ? keyed.find((entry) => legendCaption(entry) !== key) : undefined;
    };
    /** @type {ListedEntry[]} */
    const listed = [];
    // Python's rows the Result does not show, e.g. those of a track shown again.
    /** @type {Set<ListedEntry>} */
    const unshown = new Set();
    /** @type {Set<string>} */
    const generatedCaptions = new Set();
    // A row without a key record (a rename saved before U3a) is the intent
    // entry with its caption, unless another row of the Result is Python's
    // row of that entry's key: then it is Python's own row, e.g. a deleted
    // row hidden under the caption a rename took (review M1).
    /** @param {typeof rows[number]} row */
    const namedBy = (row) => {
      // An editor row shows its key; a row of Python's that shows the key of
      // an intent entry is that entry's row.
      const keyed = intentEntry((each) => generatedCaption(each) === row.shownKey);
      if (row.editor || keyed) return row.editor ? keyed || null : null;
      const entry = intentEntry((each) => legendCaption(each) === row.shownKey);
      return entry && (row.editor || !rows.some((other) => (
        other !== row && (other.recordedKey ?? other.shownKey) === generatedCaption(entry)
      ))) ? entry : null;
    };
    rows.forEach((row) => {
      if (asShown ? !row.shown : !row.pythonShown) return;
      const named = row.recordedKey === null ? namedBy(row) : null;
      if (row.editor) {
        if (!named && !asShown) return;
      } else {
        const key = row.recordedKey ?? (named ? generatedCaption(named) : row.shownKey);
        if (!asShown && deletedKeys.has(key)) return;
        // Python drew the row of a hidden track: it stays in the inventory.
        if (!asShown && hiddenKeys.has(key)) {
          generatedCaptions.add(key);
          return;
        }
        const caption = asShown ? row.shownKey : (legendCaption(renameOf(key)) || legendCaption(named) || key);
        const entry = intentEntry((each) => generatedCaption(each) === key && legendCaption(each) === caption)
          || intentEntry((each) => generatedCaption(each) === key);
        generatedCaptions.add(key);
        const item = {
          caption, originalCaption: key, color: row.color, xPos: row.xPos, yPos: row.yPos, featureIds: entry?.featureIds || []
        };
        listed.push(item);
        if (!asShown && !row.shown && !row.merged && known.has(key)) unshown.add(item);
        return;
      }
      listed.push({
        caption: asShown ? row.name : legendCaption(named) || row.name,
        originalCaption: row.shownKey,
        color: row.color,
        xPos: row.xPos,
        yPos: row.yPos,
        featureIds: named?.featureIds || []
      });
    });
    /** @param {ListedEntry} item */
    const keyOf = (item) => item.originalCaption;
    const placed = placeByInventory(listed.filter((item) => !unshown.has(item)),
      [...unshown].sort((a, b) => inventory.indexOf(keyOf(a)) - inventory.indexOf(keyOf(b))), inventory, keyOf);
    return { entries: placed, generatedCaptions };
  };

  // The inventory after a draw: while an edited order is replayed, the
  // surviving categories keep their default order and new ones follow it;
  // otherwise the Result shows the renderer's order, which becomes the default
  // order. Deleted entries keep their place.
  /**
   * @param {DrawingState} drawing
   * @param {{ inventory: string[], orderEdited: boolean, generatedCaptions: Set<string>, rendered: string[] }} drawn
   */
  const drawnInventory = (drawing, { inventory, orderEdited, generatedCaptions, rendered }) => {
    const retainedCaptions = new Set([
      ...generatedCaptions,
      ...drawing.deletedLegendEntries.value.map(entry => entry.originalCaption || entry.caption)
    ]);
    const surviving = inventory.filter(caption => retainedCaptions.has(caption));
    if (orderEdited) return [...new Set([...surviving, ...rendered])];
    // A deleted entry the renderer no longer shows keeps its place, so a
    // Restore returns it there (OV-154).
    const order = [...new Set(rendered)];
    return placeByInventory(order, surviving.filter(caption => !order.includes(caption)), inventory, (caption) => caption);
  };
  /** @param {LegendRowListing} listing */
  const renderedCaptions = ({ entries, generatedCaptions }) => entries
    .map(entry => entry.originalCaption)
    .filter(caption => generatedCaptions.has(caption));

  /**
   * The renamed rows a draw leaves out (OV-120): each keeps its generated
   * caption and its rename, unless the row is deleted.
   * @param {DrawingState} drawing
   * @param {{ previous: EntryRecord[], drawn: EntryRecord[] }} rows
   * @returns {EntryRecord[]}
   */
  const dormantAfterDraw = (drawing, { previous, drawn }) => {
    const drawnOriginals = new Set(drawn.map((entry) => String(entry?.originalCaption || entry?.caption || '').trim()));
    const deleted = new Set(drawing.deletedLegendEntries.value
      .map((entry) => String(entry?.originalCaption || entry?.caption || '').trim()));
    /** @type {Set<string>} */
    const seen = new Set();
    return previous.flatMap((entry) => {
      const caption = String(entry?.caption || '').trim();
      const originalCaption = String(entry?.originalCaption || '').trim();
      if (!caption || !originalCaption || caption === originalCaption || drawnOriginals.has(originalCaption)
        || deleted.has(originalCaption) || seen.has(originalCaption)) return [];
      seen.add(originalCaption);
      return [{ caption, originalCaption, color: entry.color, featureIds: [] }];
    });
  };

  // A drawn Result's inventory (OV-1486). Python's order of its rows is the
  // order the executor recorded when it replayed an edited order over them
  // (`recordedLegendOrder`), else the order the Result shows. The rows a
  // reorder moved are never Python's order: a row Python placed first that an
  // edited order put last would become a default the edited order then equals.
  /**
   * @param {DrawingState} drawing
   * @param {SVGSVGElement} svg
   * @param {LegendRowListing} listing The rows the Result lists.
   * @param {{ inventory: string[], orderEdited: boolean }} previous The inventory before
   *   the draw, and whether the draw replayed an edited order (read when nothing is recorded).
   * @returns {string[]}
   */
  const resultInventory = (drawing, svg, listing, { inventory, orderEdited }) => {
    const recorded = recordedLegendOrder(svg);
    return drawnInventory(drawing, {
      inventory,
      orderEdited: !recorded && orderEdited,
      generatedCaptions: recorded ? new Set(recorded) : listing.generatedCaptions,
      rendered: recorded || renderedCaptions(listing)
    });
  };

  /** @param {string[]} liveResultIdentities */
  const pruneResultInventories = (liveResultIdentities) => {
    const live = new Set(liveResultIdentities);
    [...inventoryByResult.keys()].forEach((key) => {
      if (!live.has(key)) inventoryByResult.delete(key);
    });
    [...replayedInventoryResults].forEach((key) => {
      if (!live.has(key)) replayedInventoryResults.delete(key);
    });
  };

  // A Result about to be displayed has no inventory yet: the Python keys of
  // the rows the drawing lists on it, as the renderer drew them. The stored SVG of a Result changes
  // only through edits made while it is displayed, so this runs before the
  // editor intent is projected onto it. The displayed Result keeps its own
  // inventory in `originalLegendOrder` until `adoptResultInventory`.
  /**
   * @param {SVGSVGElement} svg
   * @param {{ resultIdentity?: string, liveResultIdentities?: string[] }} [options]
   * @returns {string[]} The Result's inventory.
   */
  const captureResultInventory = (svg, { resultIdentity: identity = '', liveResultIdentities = [] } = {}) => {
    const drawing = state.activeDrawing();
    pruneResultInventories(liveResultIdentities);
    const stored = inventoryByResult.get(identity);
    if (!identity || stored) return stored || [];
    const mounted = listLegendRows(svg, legendIntentOf(drawing));
    const inventory = mounted
      ? resultInventory(drawing, svg, mounted, {
        inventory: originalLegendOrder.value, orderEdited: replayedInventoryResults.has(identity)
      })
      : [];
    inventoryByResult.set(identity, inventory);
    replayedInventoryResults.delete(identity);
    return inventory;
  };

  // E1: a Result of the other diagram mode keeps the inventory it was drawn
  // with while that mode is not shown; a Result loaded from a Session's
  // `otherModeResult` arrives with its saved one. The inventory a Result
  // arrives with is the one it left with, which followed every live edit made
  // while it was shown (a Legend rename takes the row's place in it, OV-348),
  // so it replaces a stored copy; `adoptResultInventory` adopts it when the
  // Result is displayed.
  /** @param {string | undefined} identity @param {string[]} inventory */
  const rememberResultInventory = (identity, inventory) => {
    if (!identity || !Array.isArray(inventory)) return;
    inventoryByResult.set(identity, [...inventory]);
  };

  // `originalLegendOrder` is the displayed Result's inventory: after a Result
  // is displayed and after a History restore wrote another Result's value.
  // `restored` marks a restore (Session Load, checkpoint) that installed the
  // displayed Result's own inventory.
  /** @param {string | undefined} identity @param {{ restored?: boolean }} [options] */
  const adoptResultInventory = (identity, { restored = false } = {}) => {
    if (!identity) return;
    const stored = inventoryByResult.get(identity);
    if (stored && !restored) originalLegendOrder.value = [...stored];
    else inventoryByResult.set(identity, [...originalLegendOrder.value]);
  };

  // The drawing's Legend list on the displayed Result (`listLegendRows`), and
  // after a draw the inventory and the dormant renames it leaves.
  /** @param {{ replaceGeneratedInventory?: boolean, liveResultIdentities?: string[] }} [options] */
  const extractLegendEntries = ({ replaceGeneratedInventory = false, liveResultIdentities = [] } = {}) => {
    const drawing = state.activeDrawing();
    readLegendRowNames();
    if (!svgContainer.value) {
      replaceLegendRows(drawing.slices.paint, []);
      return;
    }

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) {
      replaceLegendRows(drawing.slices.paint, []);
      return;
    }

    const previousEntries = drawing.legendEntries.value || [];
    const previousDormant = drawing.dormantLegendEntries.value || [];
    const mounted = listLegendRows(svg, legendIntentOf(drawing));
    if (!mounted) {
      replaceLegendRows(drawing.slices.paint, []);
      return;
    }
    const visuallySortedEntries = mounted.entries;

    // The entries shown before this extraction decide, as in Generate, whether
    // the displayed Result replays an edited order (D-08).
    const orderEdited = isLegendOrderEdited(drawing.legendEntries.value, originalLegendOrder.value);
    replaceLegendRows(drawing.slices.paint, listedRowsOf(visuallySortedEntries));
    // OV-120: a draw keeps the rename of each row it does not draw; the rename
    // applies again when a later Generate draws the row.
    if (replaceGeneratedInventory) {
      replaceDormantLegendRows(drawing.slices.paint, keptRowsOf(dormantAfterDraw(drawing, {
        previous: [...previousEntries, ...previousDormant],
        drawn: visuallySortedEntries
      })));
    }

    const identity = readActiveResultIdentity?.() || '';
    if (replaceGeneratedInventory || originalLegendOrder.value.length === 0) {
      // Keep deletion intent; live editor extraction alone must not advance
      // the accepted generated inventory.
      originalLegendOrder.value = resultInventory(drawing, svg, mounted, { inventory: originalLegendOrder.value, orderEdited });
      if (identity) {
        // A draw that replays an edited order shows it on every Result.
        if (replaceGeneratedInventory && orderEdited) {
          liveResultIdentities.filter((key) => key !== identity).forEach((key) => {
            inventoryByResult.delete(key);
            replayedInventoryResults.add(key);
          });
        }
        inventoryByResult.set(identity, [...originalLegendOrder.value]);
      }
    } else if (identity && !inventoryByResult.has(identity)) {
      inventoryByResult.set(identity, [...originalLegendOrder.value]);
    }

    // Kept as a Session holds them, so a Save and Load round trip keeps them (B24).
    if (Object.keys(originalLegendColors.value).length === 0 && visuallySortedEntries.length > 0) {
      visuallySortedEntries.forEach((entry) => {
        const color = normalizeOptionalHexColor(entry.color);
        if (color) originalLegendColors.value[entry.originalCaption] = color;
      });
    }
  };

  // The Legend color of the rows of `keys`: the drawing's override of each key
  // and the color of each listed row of them (D-35: the keys a merged row
  // stands for, `legendRowKeys`).
  /** @param {DrawingState} drawing @param {readonly PythonLegendKey[]} keys @param {string} color @param {{ inherit?: readonly PythonLegendKey[] }} [options] */
  const writeLegendRowColor = (drawing, keys, color, { inherit = [] } = {}) => {
    recolorLegendRows(drawing.slices.paint, keys, color);
    // The rows of `inherit` follow their type's default color again: no Legend color of their own.
    inherit.forEach((key) => setLegendColor(drawing.slices.paint, key, null));
  };

  // A Legend row color writes the editor intent only; the composition root
  // shows it on the Result through the executor (`editEditorIntent`). A color
  // the row already shows is no edit (OV-282: the palette's, for a row the
  // palette colors). It colors every row the shown row stands for (D-35).
  /** @param {number} idx @param {string} newColor */
  const updateLegendEntryColor = (idx, newColor) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = drawing.legendEntries.value[idx];
    if (!entry?.caption) return false;
    if (normalizedColor(readLegendRowColors(drawing).listed(entry)) === normalizedColor(newColor)) return false;
    writeLegendRowColor(drawing, legendRowKeys(svgContainer.value?.querySelector?.('svg'), legendEntryKey(entry)), newColor);
    return true;
  };

  // Delete a Legend row: the deleted list keeps it for Restore and Generate.
  /** @param {number} idx */
  const deleteLegendEntry = (idx) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    return deleteLegendRow(drawing.slices.paint, idx);
  };

  // The editor list with `entry` at its place in the default order: before the
  // first entry the generated inventory ranks after it, else after the last
  // ranked entry; an entry the inventory does not rank goes last, as the
  // default order puts the editor's own rows (`defaultLegendEntryOrder`).
  /**
   * @param {EntryRecord[]} entries
   * @param {EntryRecord} entry
   * @param {string[]} inventory
   */
  const withEntryAtDefaultPlace = (entries, entry, inventory) => {
    /** @param {EntryRecord} item */
    const rank = (item) => inventory.indexOf(generatedCaption(item));
    const own = rank(entry);
    if (own < 0) return [...entries, entry];
    const next = entries.findIndex((item) => rank(item) > own);
    const at = next >= 0 ? next : entries.reduce((last, item, index) => (rank(item) >= 0 ? index + 1 : last), 0);
    return [...entries.slice(0, at), entry, ...entries.slice(at)];
  };

  // OV-154: Restore and Restore all in the Legend editor. A restored row
  // leaves the deleted list, so Generate draws it again, and the editor list
  // takes it at its place in the default order (PD-OI-063), in the color
  // Generate gives it. Returns the Python keys of the restored rows, which the
  // composition root shows through the port (Python's row returns as Python
  // drew it, its stroke included) and asks Python for when the Result lacks
  // them (O-2); false when no row returned.
  /** @param {number[] | null} [indexes] Indexes into the deleted list; every row when omitted. */
  const restoreDeletedLegendEntries = (indexes = null) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    /** @type {EntryRecord[]} */
    const deleted = drawing.deletedLegendEntries.value || [];
    const picked = new Set(indexes ?? deleted.map((/** @type {unknown} */ _, /** @type {number} */ index) => index));
    const returning = deleted.filter((/** @type {unknown} */ _, /** @type {number} */ index) => picked.has(index));
    if (!svgContainer.value?.querySelector?.('svg') || returning.length === 0) return false;
    const inventory = (originalLegendOrder.value || []).map((/** @type {unknown} */ caption) => String(caption || '').trim());
    const context = displayedLegendRowContext(state, drawing);
    const shown = readLegendRowColors(drawing);
    let entries = [...(drawing.legendEntries.value || [])];
    for (const entry of returning) {
      const key = legendEntryKey(entry);
      const recorded = Object.hasOwn(drawing.legendColorOverrides, key) ? drawing.legendColorOverrides[key] : undefined;
      // OV-282: a palette row returns in the palette's color it shows, never the swatch it recorded.
      const color = String(recorded || legendRowRules(legendEntryKey(entry), context)[0]?.color || shown.deleted(entry));
      entries = withEntryAtDefaultPlace(entries, { ...entry, color }, inventory);
    }
    restoreLegendRows(drawing.slices.paint, picked, listedRowsOf(entries));
    return returning.map(generatedCaption).filter((key) => inventory.includes(key));
  };

  return {
    adoptRestoredLegend,
    adoptResultInventory,
    captureResultInventory,
    deleteLegendEntry,
    extractLegendEntries,
    findLegendRowNamed,
    keepBothNamesOf,
    legendEntryMerged,
    legendEntryName,
    legendRowIsNamed,
    legendRowName,
    moveLegendRecords,
    onLegendGeometryChanged,
    prepareFileLegendEntries,
    readLegendRowNames,
    rememberResultInventory,
    restoreDeletedLegendEntries,
    setLegendGeometryChangedHandler,
    shownLegendEntryName,
    shownLegendName,
    updateLegendEntryColor,
    writeLegendRowColor
  };
};
