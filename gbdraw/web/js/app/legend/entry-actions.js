// @ts-check
/** @import { DrawingState } from '../../state.js' */
import { normalizeOptionalHexColor, resolveColorToHex, toNativeColorInputValue } from '../../utils/color-utils.js';
import {
  defaultLegendCaptionOrder,
  getAllFeatureLegendGroups,
  isLegendOrderEdited,
  recordedLegendOrder,
  resultLegendRows
} from '../../services/legend-svg.js';
import { diffLegendIntents, legendRowRules } from '../../services/specific-color-rules.js';

const normalizedColor = (value) => {
  const resolved = String(resolveColorToHex(String(value || '').trim()) || value || '').trim().toLowerCase();
  return resolved.startsWith('#') ? toNativeColorInputValue(resolved) : resolved;
};

const legendCaption = (entry) => String(entry?.caption || '').trim();
const generatedCaption = (entry) => String(entry?.originalCaption || entry?.caption || '').trim();

// The entries in a caption order: entries the order does not list follow in
// their order; with `keepFollowed`, entries whose listed ones already follow
// the order keep it (B18), as `orderLegendEntries` orders a Legend group.
/**
 * @param {Record<string, any>[]} entries
 * @param {string[]} captions
 * @param {{ keepFollowed?: boolean }} [options]
 */
const entriesInOrder = (entries, captions, { keepFollowed = false } = {}) => {
  /** @type {Map<string, number>} */
  const rank = new Map();
  captions.forEach((caption) => { if (caption && !rank.has(caption)) rank.set(caption, rank.size); });
  const listed = entries.map((entry) => rank.get(legendCaption(entry))).filter((value) => value !== undefined);
  if (keepFollowed && listed.every((value, index) => index === 0 || listed[index - 1] < value)) return entries;
  /** @param {Record<string, any>} entry */
  const at = (entry) => rank.get(legendCaption(entry)) ?? Infinity;
  return entries.map((entry, index) => ({ entry, index }))
    .sort((left, right) => (at(left.entry) - at(right.entry)) || (left.index - right.index))
    .map(({ entry }) => entry);
};

/**
 * @typedef {object} LegendEntryActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(() => string | undefined) | null} [readActiveResultIdentity]
 *   The preview owner's runtime identity of the mounted Result.
 */

// The Legend entry owner (U3a, R1): its writers edit the drawing's Legend
// intent (`legendEntries`, `deletedLegendEntries`, the override stores) and
// write no row; the composition root shows the intent on the displayed Result
// through the port (`showEditorIntent`), whose executor writes the rows.
/** @param {LegendEntryActionsOptions} options */
export const createLegendEntryActions = ({
  state,
  readActiveResultIdentity = null
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
  /** @param {{ from?: Record<string, any>[] | null }} [options] */
  const adoptRestoredLegend = ({ from = null } = {}) => {
    const drawing = state.activeDrawing();
    const svg = svgContainer.value?.querySelector?.('svg');
    const restored = drawing.legendEntries.value || [];
    if (!svg || results.value.length < 2 || !Array.isArray(from)) return false;
    const mounted = listLegendRows(svg, { entries: restored, dormant: drawing.dormantLegendEntries.value || [], asShown: true });
    if (!mounted) return false;
    const fromCaptions = new Set(from.map(legendCaption));
    if (mounted.entries.length === fromCaptions.size && mounted.entries.every((entry) => fromCaptions.has(entry.caption))) {
      return false;
    }
    /** @param {Record<string, any>[]} list @param {(entry: Record<string, any>) => string} key */
    const captionMap = (list, key) => new Map(list.filter(legendCaption).map((entry) => [key(entry), entry]));
    const before = captionMap(from, legendCaption);
    const beforeByGenerated = captionMap(from, generatedCaption);
    const after = captionMap(restored, legendCaption);
    const afterByGenerated = captionMap(restored, generatedCaption);
    /** @type {Map<string, Record<string, any>>} */
    const renames = new Map();
    before.forEach((entry, caption) => {
      const next = afterByGenerated.get(generatedCaption(entry));
      if (next && !after.has(caption) && !before.has(legendCaption(next))) renames.set(caption, next);
    });
    const renamed = new Set([...renames.values()].map(legendCaption));
    // A row no Result's inventory lists is an editor row; a returning row is
    // this Result's own or an editor row.
    const inventories = [...inventoryByResult.values(), originalLegendOrder.value];
    /** @param {Record<string, any>} entry */
    const editorRow = (entry) => !inventories.some((inventory) => inventory.includes(generatedCaption(entry)));
    const returning = [...after.values()].filter((entry) => (
      !before.has(legendCaption(entry)) && !renamed.has(legendCaption(entry))
      && (originalLegendOrder.value.includes(generatedCaption(entry)) || editorRow(entry))
    ));
    /** @param {Record<string, any>} entry @param {Record<string, any>} next */
    const stepColor = (entry, next) => {
      const previous = before.get(legendCaption(next)) || beforeByGenerated.get(generatedCaption(next));
      return !previous || normalizedColor(previous.color) !== normalizedColor(next.color) ? next.color : entry.color;
    };
    let entries = mounted.entries.flatMap((entry) => {
      const next = renames.get(entry.caption) || after.get(entry.caption);
      if (!next) return before.has(entry.caption) ? [] : [entry];
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
        ? defaultLegendCaptionOrder(entries, [...originalLegendOrder.value])
        : [...after.keys()];
      entries = entriesInOrder(entries, order, { keepFollowed: returning.length === 0 });
    }
    drawing.legendEntries.value = entries;
    return true;
  };

  // The Legend rows the specific-color rules draw, computed on the listed rows
  // (R1, R13): the rows the rules own (a caption and color of the previous rule
  // intents) and the listed rows of a wanted caption and color are diffed with
  // the wanted rows; a listed row the rules do not own with a wanted caption
  // in another color stops the commit. The rule owner runs its transition and
  // then `apply` in one History step, while `isCurrent` holds; `apply` writes
  // the rows into the intent, the row a rename draws in the place of the row
  // it renames (OV-158), and returns the rows the displayed Result shows at
  // once until Python draws them (gaps 1 and 2; null when the rows only change
  // color). Resolves to false when the Result or the caller became stale.
  /**
   * @param {Record<string, any>[]} intents
   * @param {{ drawing?: DrawingState, previousFileIntents?: Record<string, any>[], isCurrent?: () => boolean, placement?: { caption: string, at: string } | null }} [options]
   */
  const prepareFileLegendEntries = async (intents, {
    drawing = state.activeDrawing(), previousFileIntents = [], isCurrent = () => true, placement = null
  } = {}) => {
    const mountedSvg = svgContainer.value?.querySelector('svg');
    if (getAllFeatureLegendGroups(mountedSvg).length === 0) {
      return { diff: { add: [], update: [], remove: [], unchanged: [] }, isCurrent: () => true, apply: () => null };
    }
    /** @type {Map<string, Set<string>>} */
    const provenance = new Map();
    for (const entry of previousFileIntents) {
      const caption = legendCaption(entry);
      if (!provenance.has(caption)) provenance.set(caption, new Set());
      provenance.get(caption)?.add(normalizedColor(entry?.color));
    }
    const listed = (drawing.legendEntries.value || [])
      .map((/** @type {Record<string, any>} */ entry) => ({ caption: legendCaption(entry), color: normalizedColor(entry.color) }));
    const owned = listed.filter((entry) => provenance.get(entry.caption)?.has(entry.color));
    const desiredByCaption = new Map(intents.map((intent) => [intent.caption, normalizedColor(intent.color)]));
    for (const intent of intents) {
      const existing = listed.find((entry) => entry.caption === intent.caption);
      if (existing && !owned.includes(existing) && existing.color !== normalizedColor(intent.color)) {
        throw new Error(`Legend entry "${intent.caption}" already exists with a different color.`);
      }
    }
    const reusable = listed.filter((entry) => !owned.includes(entry) && desiredByCaption.get(entry.caption) === entry.color);
    const diff = diffLegendIntents([...owned, ...reusable], intents);
    if (!isCurrent() || svgContainer.value?.querySelector('svg') !== mountedSvg) return false;
    return {
      diff,
      isCurrent: () => svgContainer.value?.querySelector('svg') === mountedSvg,
      apply: () => {
        const removed = new Set(diff.remove.map(({ caption }) => caption));
        const colors = new Map(diff.update.map(({ caption, color }) => [caption, color]));
        const added = new Map(diff.add.map(({ caption, color }) => [caption, { caption, originalCaption: caption, color, featureIds: [] }]));
        /** @type {Record<string, any>[]} */
        let entries = (drawing.legendEntries.value || []).map((/** @type {Record<string, any>} */ entry) => (
          colors.has(legendCaption(entry)) ? { ...entry, color: colors.get(legendCaption(entry)) } : entry
        ));
        // OV-158 (Owner decision 2026-10-07): the row a rename draws takes the
        // place of the row it renames, so the editor list records the edited
        // order Generate replays (PD-OI-063).
        let before = '';
        if (placement && intents.some((intent) => intent.caption === placement.caption)) {
          const row = added.get(placement.caption) || entries.find((entry) => legendCaption(entry) === placement.caption);
          const others = entries.filter((entry) => legendCaption(entry) !== placement.caption);
          const at = others.findIndex((entry) => legendCaption(entry) === placement.at);
          if (row && at >= 0) {
            entries = [...others.slice(0, at), row, ...others.slice(at)];
            if (added.delete(placement.caption)) before = placement.at;
          }
        }
        drawing.legendEntries.value = [...entries.filter((entry) => !removed.has(legendCaption(entry))), ...added.values()];
        const rows = {
          add: diff.add.map(({ caption, color }) => ({ caption, color, ...(before && caption === placement?.caption ? { before } : {}) })),
          retire: [...removed]
        };
        return rows.add.length > 0 || rows.retire.length > 0 ? rows : null;
      }
    };
  };

  // The drawing's Legend intent that names a Result's rows.
  /** @param {DrawingState} drawing */
  const legendIntentOf = (drawing) => ({
    entries: drawing.legendEntries.value || [],
    dormant: drawing.dormantLegendEntries.value || [],
    deleted: drawing.deletedLegendEntries.value || [],
    inventory: originalLegendOrder.value || []
  });

  /**
   * The drawing's Legend list on the rows of a Result (U3b: one direction,
   * intent to shown rows). A row of Python's is Python's key: the record a
   * rename keeps, else the key of the intent entry that names the row (a
   * Result saved before U3a has no record). It is listed unless the intent
   * deletes that key, named by the intent's rename of it (else a dormant one,
   * OV-120), and keeps that entry's feature ids. An editor row is listed when
   * the intent lists it. The Result gives the rows' order and the paint the
   * executor showed from the intent. `asShown` lists the rows the Result shows
   * under the keys they show (B19 reads the list a Result showed). Also the
   * generated captions: the keys of the listed rows of Python's.
   * @param {SVGSVGElement} svg
   * @param {{ entries: Record<string, any>[], dormant?: Record<string, any>[], deleted?: Record<string, any>[],
   *   inventory?: string[], asShown?: boolean }} intent
   */
  const listLegendRows = (svg, { entries, dormant = [], deleted = [], inventory = [], asShown = false }) => {
    const rows = resultLegendRows(svg);
    if (!rows) return null;
    const known = new Set(inventory);
    const deletedKeys = new Set(deleted.map(generatedCaption).filter((key) => known.has(key)));
    /** @param {(entry: Record<string, any>) => boolean} test */
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
    /** @type {Record<string, any>[]} */
    const listed = [];
    /** @type {Set<string>} */
    const generatedCaptions = new Set();
    // A row without a key record (a rename saved before U3a) is the intent
    // entry with its caption, unless another row of the Result is Python's
    // row of that entry's key: then it is Python's own row, e.g. a deleted
    // row hidden under the caption a rename took (review M1).
    /** @param {typeof rows[number]} row */
    const namedBy = (row) => {
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
        const caption = asShown ? row.shownKey : (legendCaption(renameOf(key)) || legendCaption(named) || key);
        const entry = intentEntry((each) => generatedCaption(each) === key && legendCaption(each) === caption)
          || intentEntry((each) => generatedCaption(each) === key);
        generatedCaptions.add(key);
        listed.push({
          caption, originalCaption: key, color: row.color, xPos: row.xPos, yPos: row.yPos, featureIds: entry?.featureIds || []
        });
        return;
      }
      listed.push({
        caption: row.shownKey,
        originalCaption: named ? generatedCaption(named) : row.shownKey,
        color: row.color,
        xPos: row.xPos,
        yPos: row.yPos,
        featureIds: named?.featureIds || []
      });
    });
    return { entries: listed, generatedCaptions };
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
    // A deleted entry the renderer no longer shows goes before the first
    // rendered entry that followed it, so a Restore returns it there (OV-154).
    const order = [...new Set(rendered)];
    surviving.filter(caption => !order.includes(caption)).forEach((caption) => {
      const next = inventory.slice(inventory.indexOf(caption) + 1).find(later => order.includes(later));
      order.splice(next === undefined ? order.length : order.indexOf(next), 0, caption);
    });
    return order;
  };
  const renderedCaptions = ({ entries, generatedCaptions }) => entries
    .map(entry => entry.originalCaption)
    .filter(caption => generatedCaptions.has(caption));

  /**
   * The renamed rows a draw leaves out (OV-120): each keeps its generated
   * caption and its rename, unless the row is deleted.
   * @param {DrawingState} drawing
   * @param {{ previous: Record<string, any>[], drawn: Record<string, any>[] }} rows
   * @returns {Record<string, any>[]}
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
    // A Result whose rows were reordered keeps Python's order as a record (L1).
    const recorded = recordedLegendOrder(svg);
    const inventory = mounted
      ? drawnInventory(drawing, {
        inventory: originalLegendOrder.value,
        orderEdited: !recorded && replayedInventoryResults.has(identity),
        generatedCaptions: recorded ? new Set(recorded) : mounted.generatedCaptions,
        rendered: recorded || renderedCaptions(mounted)
      })
      : [];
    inventoryByResult.set(identity, inventory);
    replayedInventoryResults.delete(identity);
    return inventory;
  };

  // E1: a Result of the other diagram mode keeps the inventory it was drawn
  // with while that mode is not shown; a Result loaded from a Session's
  // `otherModeResult` arrives with its saved one. Kept only for a Result with
  // no stored inventory; `adoptResultInventory` adopts it when it is displayed.
  /** @param {string | undefined} identity @param {string[]} inventory */
  const rememberResultInventory = (identity, inventory) => {
    if (!identity || inventoryByResult.has(identity) || !Array.isArray(inventory)) return;
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
    if (!svgContainer.value) {
      drawing.legendEntries.value = [];
      return;
    }

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) {
      drawing.legendEntries.value = [];
      return;
    }

    const previousEntries = drawing.legendEntries.value || [];
    const previousDormant = drawing.dormantLegendEntries.value || [];
    const mounted = listLegendRows(svg, legendIntentOf(drawing));
    if (!mounted) {
      drawing.legendEntries.value = [];
      return;
    }
    const { entries: visuallySortedEntries, generatedCaptions } = mounted;

    // The entries shown before this extraction decide, as in Generate, whether
    // the displayed Result replays an edited order (D-08).
    const orderEdited = isLegendOrderEdited(drawing.legendEntries.value, originalLegendOrder.value);
    drawing.legendEntries.value = visuallySortedEntries;
    // OV-120: a draw keeps the rename of each row it does not draw; the rename
    // applies again when a later Generate draws the row.
    if (replaceGeneratedInventory) {
      drawing.dormantLegendEntries.value = dormantAfterDraw(drawing, {
        previous: [...previousEntries, ...previousDormant],
        drawn: visuallySortedEntries
      });
    }

    const identity = readActiveResultIdentity?.() || '';
    if (replaceGeneratedInventory || originalLegendOrder.value.length === 0) {
      // Keep deletion intent; live editor extraction alone must not advance
      // the accepted generated inventory.
      originalLegendOrder.value = drawnInventory(drawing, {
        inventory: originalLegendOrder.value,
        orderEdited,
        generatedCaptions,
        rendered: renderedCaptions(mounted)
      });
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
        if (color) originalLegendColors.value[entry.caption] = color;
      });
    }
  };

  // A Legend row color writes the editor intent only; the composition root
  // shows it on the Result through the executor (`editEditorIntent`).
  /** @param {number} idx @param {string} newColor */
  const updateLegendEntryColor = (idx, newColor) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = drawing.legendEntries.value[idx];
    if (!entry?.caption) return false;
    if (normalizedColor(entry.color) === normalizedColor(newColor)) return false;
    entry.color = newColor;
    drawing.legendColorOverrides[entry.caption] = newColor;
    return true;
  };

  // Delete a Legend row: the deleted list keeps it for Restore and Generate.
  /** @param {number} idx */
  const deleteLegendEntry = (idx) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = drawing.legendEntries.value[idx];
    if (!entry) return false;
    drawing.deletedLegendEntries.value.push({ ...entry });
    drawing.legendEntries.value = drawing.legendEntries.value.filter((/** @type {unknown} */ _, /** @type {number} */ index) => index !== idx);
    return true;
  };

  // The editor list with `entry` at its place in the default order: before the
  // first entry the generated inventory ranks after it, else after the last
  // ranked entry; an entry the inventory does not rank goes last, as the
  // default order puts the editor's own rows (`defaultLegendEntryOrder`).
  /**
   * @param {Record<string, any>[]} entries
   * @param {Record<string, any>} entry
   * @param {string[]} inventory
   */
  const withEntryAtDefaultPlace = (entries, entry, inventory) => {
    /** @param {Record<string, any>} item */
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
    const deleted = drawing.deletedLegendEntries.value || [];
    const picked = new Set(indexes ?? deleted.map((/** @type {unknown} */ _, /** @type {number} */ index) => index));
    const returning = deleted.filter((/** @type {unknown} */ _, /** @type {number} */ index) => picked.has(index));
    if (!svgContainer.value?.querySelector?.('svg') || returning.length === 0) return false;
    const inventory = (originalLegendOrder.value || []).map((/** @type {unknown} */ caption) => String(caption || '').trim());
    const context = {
      rules: drawing.manualSpecificRules,
      legendEntries: drawing.legendEntries.value || [],
      originalLegendOrder: inventory
    };
    let entries = [...(drawing.legendEntries.value || [])];
    for (const entry of returning) {
      const caption = legendCaption(entry);
      const color = String(drawing.legendColorOverrides[caption] || legendRowRules(caption, context)[0]?.color || entry.color);
      entries = withEntryAtDefaultPlace(entries, { ...entry, color }, inventory);
    }
    drawing.deletedLegendEntries.value = deleted.filter((/** @type {unknown} */ _, /** @type {number} */ index) => !picked.has(index));
    drawing.legendEntries.value = entries;
    return returning.map(generatedCaption).filter((key) => inventory.includes(key));
  };

  return {
    adoptRestoredLegend,
    adoptResultInventory,
    captureResultInventory,
    deleteLegendEntry,
    extractLegendEntries,
    onLegendGeometryChanged,
    prepareFileLegendEntries,
    rememberResultInventory,
    restoreDeletedLegendEntries,
    setLegendGeometryChangedHandler,
    updateLegendEntryColor
  };
};
