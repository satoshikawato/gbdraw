import { getAllFeatureLegendGroups, orderLegendEntries } from './utils.js';

export const createLegendSortActions = ({ state, extractLegendEntries, previewRuntime = null }) => {
  const { svgContainer, legendEntries, originalLegendOrder } = state;

  const getCurrentSvg = () => {
    if (!svgContainer.value) return null;
    return svgContainer.value.querySelector('svg');
  };

  const persistLegendOrder = () => {
    previewRuntime?.commitActiveResultEdit('legend-order');
    extractLegendEntries();
  };

  const applyLegendEntryOrder = (captionOrder) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const svg = getCurrentSvg();
    if (!svg) return;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return;

    let changed = false;
    for (const targetGroup of targetGroups) {
      changed = orderLegendEntries(targetGroup, captionOrder) || changed;
    }

    if (!changed) {
      const currentOrder = legendEntries.value.map((entry) => entry.caption);
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

  const getVisibleLegendOrder = () => legendEntries.value.map((entry) => entry.caption).filter(Boolean);

  const moveLegendEntryUp = (idx) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (idx <= 0) return;
    swapLegendEntries(idx, idx - 1);
  };

  const moveLegendEntryDown = (idx) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (idx >= legendEntries.value.length - 1) return;
    swapLegendEntries(idx, idx + 1);
  };

  const sortLegendEntries = (direction = 'asc') => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const currentOrder = getVisibleLegendOrder();
    if (currentOrder.length < 2) return;

    const currentIndex = new Map(currentOrder.map((caption, idx) => [caption, idx]));
    const sortedOrder = [...currentOrder].sort((a, b) => {
      const cmp = a.localeCompare(b, undefined, { sensitivity: 'base' });
      if (cmp === 0) {
        return (currentIndex.get(a) ?? 0) - (currentIndex.get(b) ?? 0);
      }
      return direction === 'asc' ? cmp : -cmp;
    });

    applyLegendEntryOrder(sortedOrder);
  };

  const sortLegendEntriesByDefault = () => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const currentOrder = getVisibleLegendOrder();
    if (currentOrder.length < 2) return;
    if (originalLegendOrder.value.length === 0) return;

    const currentIndex = new Map(currentOrder.map((caption, idx) => [caption, idx]));
    // A renamed generated entry keeps its generated caption as its default slot.
    const generatedCaption = new Map(legendEntries.value.map((entry) => [entry.caption, entry.originalCaption || entry.caption]));
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

    applyLegendEntryOrder(sortedOrder);
  };

  const swapLegendEntries = (idx1, idx2) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const currentOrder = getVisibleLegendOrder();
    if (idx1 < 0 || idx2 < 0 || idx1 >= currentOrder.length || idx2 >= currentOrder.length) return;

    const nextOrder = [...currentOrder];
    [nextOrder[idx1], nextOrder[idx2]] = [nextOrder[idx2], nextOrder[idx1]];
    applyLegendEntryOrder(nextOrder);
  };

  return {
    moveLegendEntryDown,
    moveLegendEntryUp,
    sortLegendEntries,
    sortLegendEntriesByDefault,
    swapLegendEntries
  };
};
