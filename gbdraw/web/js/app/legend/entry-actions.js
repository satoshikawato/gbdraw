// @ts-check
/** @import { DrawingState } from '../../state.js' */
import { normalizeOptionalHexColor, resolveColorToHex, toNativeColorInputValue } from '../../utils/color-utils.js';
import {
  defaultLegendCaptionOrder,
  getAllFeatureLegendGroups,
  getVisibleFeatureLegendGroup,
  isLegendOrderEdited,
  legendRowShown,
  recordedLegendOrder,
  parseTransformXY
} from '../../services/legend-svg.js';
import { parseCompositionMetadata } from '../legend-layout/composition-actions.js';
import { diffLegendIntents, legendRowRules } from '../../services/specific-color-rules.js';

const normalizedColor = (value) => {
  const resolved = String(resolveColorToHex(String(value || '').trim()) || value || '').trim().toLowerCase();
  return resolved.startsWith('#') ? toNativeColorInputValue(resolved) : resolved;
};

const directLegendEntryGroups = (targetGroup) => {
  const direct = Array.from(targetGroup?.children || []).filter(
    (child) => child.tagName?.toLowerCase() === 'g' && child.hasAttribute?.('data-legend-key')
  );
  return direct.length > 0
    ? direct
    : Array.from(targetGroup?.querySelectorAll?.('g[data-legend-key]') || []);
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
    originalLegendColors,
    newLegendCaption,
    newLegendColor
  } = state;

  /** @type {((options?: { commit?: boolean }) => unknown) | null} */
  let legendGeometryChangedHandler = null;

  const targetGroupKey = (targetGroup, index) => (
    String(targetGroup?.id || targetGroup?.parentElement?.id || `legend-target-${index}`)
  );

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

  const captureLegendEntryOwners = () => {
    const svg = svgContainer.value?.querySelector?.('svg');
    return svg ? getAllFeatureLegendGroups(svg).map((group, index) => ({
      target: targetGroupKey(group, index),
      entries: directLegendEntryGroups(group).map(entry => ({
        caption: String(entry.getAttribute('data-legend-key') || '').trim(),
        owner: String(entry.getAttribute('data-legend-owner') || '')
      }))
    })) : [];
  };

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
    const mounted = readMountedLegend(svg, restored, drawing.dormantLegendEntries.value || []);
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

  /**
   * Read the Legend of a mounted Result: its entries in visual order and the
   * captions the renderer generated (not the editor's direct entries). The one
   * reader of a mounted Legend; `previousEntries` keep stroke, feature ids,
   * and the generated caption of a renamed entry; a renamed row an earlier
   * Generate hid (`dormantEntries`, OV-120) keeps its generated caption too.
   * @param {SVGSVGElement} svg
   * @param {Record<string, any>[]} previousEntries
   * @param {Record<string, any>[]} [dormantEntries]
   */
  const readMountedLegend = (svg, previousEntries, dormantEntries = []) => {
    const targetGroup = getVisibleFeatureLegendGroup(svg);
    if (!targetGroup) return null;

    const entries = [];
    const generatedCaptions = new Set();

    const entryGroups = targetGroup.querySelectorAll('g[data-legend-key]');
    entryGroups.forEach((entryGroup) => {
      const caption = entryGroup.getAttribute('data-legend-key');
      if (!caption || !legendRowShown(entryGroup)) return;

      let color = '#cccccc';
      const paths = entryGroup.querySelectorAll('path');
      for (const path of paths) {
        const fill = path.getAttribute('fill');
        if (fill && fill !== 'none' && !fill.startsWith('url(')) {
          color = fill;
          break;
        }
      }

      let xPos = 0,
        yPos = 0;
      const groupTransform = parseTransformXY(entryGroup.getAttribute('transform'));
      const textEl = entryGroup.querySelector('text');
      if (textEl) {
        const textTransform = parseTransformXY(textEl.getAttribute('transform'));
        xPos = groupTransform.x + textTransform.x;
        yPos = groupTransform.y + textTransform.y;
      } else if (groupTransform.x !== 0 || groupTransform.y !== 0) {
        xPos = groupTransform.x;
        yPos = groupTransform.y;
      }

      const existingEntry = previousEntries.find((entry) => (
        entry.caption === caption
        && normalizedColor(entry.color) === normalizedColor(color)
      )) || dormantEntries.find((entry) => entry.caption === caption);
      const existingFeatureIds = existingEntry?.featureIds || [];
      const originalCaption = existingEntry?.originalCaption || caption;
      if (entryGroup.getAttribute('data-legend-owner') !== 'direct-editor') {
        generatedCaptions.add(originalCaption);
      }

      entries.push({
        caption,
        originalCaption,
        color,
        xPos,
        yPos,
        featureIds: existingFeatureIds
      });
    });

    const visuallySortedEntries = [...entries].sort((a, b) => {
      const yDelta = a.yPos - b.yPos;
      if (Math.abs(yDelta) < 1) {
        const xDelta = a.xPos - b.xPos;
        if (Math.abs(xDelta) < 1) {
          return a.caption.localeCompare(b.caption, undefined, { sensitivity: 'base' });
        }
        return xDelta;
      }
      return yDelta;
    });

    return { entries: visuallySortedEntries, generatedCaptions };
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

  // A Result about to be displayed has no inventory yet: read it from the
  // mounted Legend as the renderer drew it. The stored SVG of a Result changes
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
    const mounted = readMountedLegend(svg, drawing.legendEntries.value || [], drawing.dormantLegendEntries.value || []);
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
    const mounted = readMountedLegend(svg, previousEntries, previousDormant);
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

  // Add a Legend row of the editor's own: listed last, under a caption no
  // listed, deleted, or generated row has (a ` (n)` suffix), so Generate adds
  // it (OV-86). A listed row of the same caption and color is no new row.
  const addNewLegendEntry = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const caption = newLegendCaption.value.trim();
    if (!caption) return false;
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg || getAllFeatureLegendGroups(svg).length === 0) return false;
    if (!parseCompositionMetadata(svg).legendReflow) {
      throw new Error('This diagram has no legend reflow metadata. Regenerate it before editing the legend.');
    }
    const color = newLegendColor.value;
    const entries = drawing.legendEntries.value || [];
    newLegendCaption.value = '';
    newLegendColor.value = '#808080';
    if (entries.some((/** @type {Record<string, any>} */ entry) => (
      legendCaption(entry) === caption && normalizedColor(entry.color) === normalizedColor(color)
    ))) return false;
    const taken = new Set([
      ...[...entries, ...(drawing.deletedLegendEntries.value || []), ...(drawing.dormantLegendEntries.value || [])]
        .flatMap((/** @type {Record<string, any>} */ entry) => [legendCaption(entry), generatedCaption(entry)]),
      ...originalLegendOrder.value
    ]);
    const base = caption.replace(/\s*\(\d+\)$/, '');
    let finalCaption = caption;
    for (let counter = 1; taken.has(finalCaption); counter += 1) finalCaption = `${base} (${counter})`;
    drawing.legendEntries.value = [...entries, { caption: finalCaption, originalCaption: finalCaption, color, featureIds: [] }];
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
    addNewLegendEntry,
    captureLegendEntryOwners,
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
