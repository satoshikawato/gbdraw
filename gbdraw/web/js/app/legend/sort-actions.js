// @ts-check
/** @import { DrawingState } from '../../state.js' */
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

  /** @param {DrawingState} drawing @param {string[]} captionOrder */
  const applyLegendEntryOrder = (drawing, captionOrder) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    /** @type {Record<string, any>[]} */
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
  const getVisibleLegendOrder = (drawing) => drawing.legendEntries.value.map((entry) => entry.caption).filter(Boolean);

  const moveLegendEntryUp = (idx) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (idx <= 0) return;
    swapLegendEntries(idx, idx - 1);
  };

  const moveLegendEntryDown = (idx) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (idx >= drawing.legendEntries.value.length - 1) return;
    swapLegendEntries(idx, idx + 1);
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
    const generatedCaption = new Map(drawing.legendEntries.value.map((entry) => [entry.caption, entry.originalCaption || entry.caption]));
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

  const swapLegendEntries = (idx1, idx2) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const currentOrder = getVisibleLegendOrder(drawing);
    if (idx1 < 0 || idx2 < 0 || idx1 >= currentOrder.length || idx2 >= currentOrder.length) return;

    const nextOrder = [...currentOrder];
    [nextOrder[idx1], nextOrder[idx2]] = [nextOrder[idx2], nextOrder[idx1]];
    applyLegendEntryOrder(drawing, nextOrder);
  };

  return {
    moveLegendEntryDown,
    moveLegendEntryUp,
    sortLegendEntries,
    sortLegendEntriesByDefault,
    swapLegendEntries
  };
};
