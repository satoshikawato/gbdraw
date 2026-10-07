// @ts-check
/** @import { DrawingState } from '../../state.js' */
// Sort and Move compute the requested caption order; the Legend entry owner
// orders the mounted Legend through the `orderMountedLegend` port (R3, R13).
/**
 * @typedef {object} LegendSortActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(options?: { replaceGeneratedInventory?: boolean }) => void} extractLegendEntries
 *   The Legend entry owner's re-read of the mounted Legend.
 * @property {(captionOrder: string[], options?: { keepFollowed?: boolean }) => (boolean | null)} orderMountedLegend
 *   The Legend entry owner's one ordering of the mounted Legend; null when no Legend is mounted.
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 */

/** @param {LegendSortActionsOptions} options */
export const createLegendSortActions = ({ state, extractLegendEntries, orderMountedLegend, commitActiveResultEdit = null }) => {
  const { originalLegendOrder } = state;

  const persistLegendOrder = () => {
    commitActiveResultEdit?.('legend-order');
    extractLegendEntries();
  };

  /** @param {DrawingState} drawing */
  const applyLegendEntryOrder = (drawing, captionOrder) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const changed = orderMountedLegend(captionOrder);
    if (changed === null) return;

    if (!changed) {
      const currentOrder = drawing.legendEntries.value.map((entry) => entry.caption);
      const normalizedRequested = captionOrder.filter(Boolean);
      if (
        currentOrder.length === normalizedRequested.length &&
        currentOrder.every((caption, idx) => caption === normalizedRequested[idx])
      ) {
        return;
      }
    }

    persistLegendOrder();
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
