import { normalizeLegendOrderIntent } from './manual-intent-command.js';

const captionKey = (value) => String(value ?? '').trim().toLocaleLowerCase();

export const createLegendSortActions = ({ state, requestLegendOrderIntent }) => {
  const { legendEntries, originalLegendOrder } = state;

  const getVisibleLegendOrder = () => normalizeLegendOrderIntent(
    (Array.isArray(legendEntries.value) ? legendEntries.value : []).map((entry) => entry?.caption)
  );

  const applyLegendEntryOrder = (captionOrder, label = 'Reorder legend') => {
    const normalized = normalizeLegendOrderIntent(captionOrder);
    if (normalized.length < 2) return false;
    if (typeof requestLegendOrderIntent !== 'function') {
      throw new Error('Legend ordering requires the document-global History command.');
    }
    return requestLegendOrderIntent({ kind: 'order', order: normalized }, label);
  };

  const swapLegendEntries = (idx1, idx2) => {
    const currentOrder = getVisibleLegendOrder();
    if (idx1 < 0 || idx2 < 0 || idx1 >= currentOrder.length || idx2 >= currentOrder.length) {
      return false;
    }
    const nextOrder = [...currentOrder];
    [nextOrder[idx1], nextOrder[idx2]] = [nextOrder[idx2], nextOrder[idx1]];
    return applyLegendEntryOrder(nextOrder);
  };

  const moveLegendEntryUp = (idx) => (
    idx <= 0 ? false : swapLegendEntries(idx, idx - 1)
  );

  const moveLegendEntryDown = (idx) => (
    idx >= legendEntries.value.length - 1 ? false : swapLegendEntries(idx, idx + 1)
  );

  const sortLegendEntries = (direction = 'asc') => {
    const currentOrder = getVisibleLegendOrder();
    if (currentOrder.length < 2) return false;
    const currentIndex = new Map(currentOrder.map((caption, idx) => [captionKey(caption), idx]));
    const sortedOrder = [...currentOrder].sort((left, right) => {
      const compared = left.localeCompare(right, undefined, { sensitivity: 'base' });
      if (compared === 0) {
        return (currentIndex.get(captionKey(left)) ?? 0) - (currentIndex.get(captionKey(right)) ?? 0);
      }
      return direction === 'asc' ? compared : -compared;
    });
    return applyLegendEntryOrder(
      sortedOrder,
      direction === 'asc' ? 'Sort legend A–Z' : 'Sort legend Z–A'
    );
  };

  const sortLegendEntriesByDefault = () => {
    const currentOrder = getVisibleLegendOrder();
    const originalOrder = normalizeLegendOrderIntent(originalLegendOrder.value || []);
    if (currentOrder.length < 2 || originalOrder.length === 0) return false;
    const originalIndex = new Map(
      originalOrder.map((caption, index) => [captionKey(caption), index])
    );
    const currentIndex = new Map(
      currentOrder.map((caption, index) => [captionKey(caption), index])
    );
    const sortedOrder = [...currentOrder].sort((left, right) => {
      const leftOriginal = originalIndex.get(captionKey(left));
      const rightOriginal = originalIndex.get(captionKey(right));
      if (leftOriginal !== undefined || rightOriginal !== undefined) {
        if (leftOriginal === undefined) return 1;
        if (rightOriginal === undefined) return -1;
        if (leftOriginal !== rightOriginal) return leftOriginal - rightOriginal;
      }
      const compared = left.localeCompare(right, undefined, { sensitivity: 'base' });
      if (compared !== 0) return compared;
      return (currentIndex.get(captionKey(left)) ?? 0) - (currentIndex.get(captionKey(right)) ?? 0);
    });
    return applyLegendEntryOrder(sortedOrder, 'Restore default legend order');
  };

  return {
    moveLegendEntryDown,
    moveLegendEntryUp,
    sortLegendEntries,
    sortLegendEntriesByDefault,
    swapLegendEntries
  };
};
