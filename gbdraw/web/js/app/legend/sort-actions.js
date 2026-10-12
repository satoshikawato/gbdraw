// @ts-check
/** @import { DrawingState, LegendCaption, LegendEntry } from '../../state.js' */
// Sort and Move compute the requested caption order, write it into the
// drawing's Legend intent, and show it through the root's port once (R1, R3,
// R13): the executor orders the displayed Result's rows.
/**
 * @typedef {object} LegendSortActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {() => unknown} showLegendStructure The root's show of the Legend structure intent.
 */

/** @param {LegendSortActionsOptions} options */
export const createLegendSortActions = ({ state, showLegendStructure }) => {
  const { originalLegendOrder } = state;

  /** @param {DrawingState} drawing @param {LegendCaption[]} captionOrder */
  const applyLegendEntryOrder = (drawing, captionOrder) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    /** @type {LegendEntry[]} */
    const entries = drawing.legendEntries.value || [];
    const listed = new Set(captionOrder);
    const ordered = [
      ...[...listed].flatMap((caption) => entries.filter((entry) => entry.caption === caption)),
      ...entries.filter((entry) => !listed.has(entry.caption))
    ];
    if (ordered.every((entry, index) => entry === entries[index])) return;
    drawing.legendEntries.value = ordered;
    showLegendStructure();
  };

  /** @param {DrawingState} drawing */
  const getVisibleLegendOrder = (drawing) => /** @type {LegendEntry[]} */ (drawing.legendEntries.value).map((entry) => entry.caption).filter(Boolean);

  // Moves the row at `from` to `to` (the Legend panel's drag handle and the
  // handle's arrow keys).
  /** @param {number} from @param {number} to */
  const moveLegendEntry = (from, to) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const currentOrder = getVisibleLegendOrder(drawing);
    if (from === to || from < 0 || to < 0 || from >= currentOrder.length || to >= currentOrder.length) return;
    const nextOrder = [...currentOrder];
    nextOrder.splice(to, 0, ...nextOrder.splice(from, 1));
    applyLegendEntryOrder(drawing, nextOrder);
  };

  const sortLegendEntries = (direction = 'asc') => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const currentOrder = getVisibleLegendOrder(drawing);
    if (currentOrder.length < 2) return;

    const currentIndex = new Map(currentOrder.map((caption, idx) => [caption, idx]));
    const sortedOrder = [...currentOrder].sort((a, b) => {
      const cmp = a.localeCompare(b, undefined, { sensitivity: 'base' });
      if (cmp === 0) {
        return (currentIndex.get(a) ?? 0) - (currentIndex.get(b) ?? 0);
      }
      return direction === 'asc' ? cmp : -cmp;
    });

    applyLegendEntryOrder(drawing, sortedOrder);
  };

  const sortLegendEntriesByDefault = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const currentOrder = getVisibleLegendOrder(drawing);
    if (currentOrder.length < 2) return;
    if (originalLegendOrder.value.length === 0) return;

    const currentIndex = new Map(currentOrder.map((caption, idx) => [caption, idx]));
    // A renamed generated entry keeps its generated caption as its default slot.
    const generatedCaption = new Map(/** @type {LegendEntry[]} */ (drawing.legendEntries.value).map((entry) => [entry.caption, entry.originalCaption || entry.caption]));
    const sortedOrder = [...currentOrder].sort((a, b) => {
      const aOrigIdx = originalLegendOrder.value.indexOf(generatedCaption.get(a) ?? a);
      const bOrigIdx = originalLegendOrder.value.indexOf(generatedCaption.get(b) ?? b);

      if (aOrigIdx !== -1 && bOrigIdx !== -1) {
        return aOrigIdx - bOrigIdx;
      }
      if (aOrigIdx !== -1) return -1;
      if (bOrigIdx !== -1) return 1;
      const cmp = a.localeCompare(b, undefined, { sensitivity: 'base' });
      if (cmp !== 0) return cmp;
      return (currentIndex.get(a) ?? 0) - (currentIndex.get(b) ?? 0);
    });

    applyLegendEntryOrder(drawing, sortedOrder);
  };

  return {
    moveLegendEntry,
    sortLegendEntries,
    sortLegendEntriesByDefault
  };
};
