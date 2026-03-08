import { getAllFeatureLegendGroups, parseTransform, parseTransformXY } from './utils.js';

export const createLegendSortActions = ({ state, extractLegendEntries }) => {
  const { svgContainer, legendEntries, originalLegendOrder, selectedResultIndex, results, skipCaptureBaseConfig } =
    state;

  const getCurrentSvg = () => {
    if (!svgContainer.value) return null;
    return svgContainer.value.querySelector('svg');
  };

  const getEntryGroups = (targetGroup) => {
    const directGroups = Array.from(targetGroup?.children || []).filter(
      (child) => child.tagName?.toLowerCase() === 'g' && child.hasAttribute('data-legend-key')
    );
    if (directGroups.length > 0) return directGroups;
    return Array.from(targetGroup?.querySelectorAll('g[data-legend-key]') || []);
  };

  const getEntryCaption = (entryGroup) => String(entryGroup?.getAttribute('data-legend-key') || '').trim();

  const getEntryAnchor = (entryGroup) => {
    if (!entryGroup) return { x: 0, y: 0 };

    const groupOffset = parseTransform(entryGroup.getAttribute('transform'));
    const textEl = entryGroup.querySelector('text');
    if (textEl) {
      const textPos = parseTransformXY(textEl.getAttribute('transform'));
      return { x: groupOffset.x + textPos.x, y: groupOffset.y + textPos.y };
    }

    const colorPath = Array.from(entryGroup.querySelectorAll('path')).find((path) => {
      const fill = path.getAttribute('fill');
      return fill && fill !== 'none' && !fill.startsWith('url(');
    });
    if (colorPath) {
      const pathPos = parseTransformXY(colorPath.getAttribute('transform'));
      return { x: groupOffset.x + pathPos.x, y: groupOffset.y + pathPos.y };
    }

    return groupOffset;
  };

  const translateEntryGroupChildren = (entryGroup, deltaX, deltaY) => {
    if (!entryGroup) return;
    if (Math.abs(deltaX) < 1e-6 && Math.abs(deltaY) < 1e-6) return;

    const transformedNodes = entryGroup.querySelectorAll('[transform]');
    transformedNodes.forEach((node) => {
      const { x, y } = parseTransform(node.getAttribute('transform'));
      node.setAttribute('transform', `translate(${x + deltaX}, ${y + deltaY})`);
    });
  };

  const persistLegendSvg = (svg) => {
    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }

    extractLegendEntries();
  };

  const applyLegendEntryOrder = (captionOrder) => {
    const svg = getCurrentSvg();
    if (!svg) return;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return;

    let changed = false;

    for (const targetGroup of targetGroups) {
      const entryGroups = getEntryGroups(targetGroup);
      if (entryGroups.length < 2) continue;

      const groupByCaption = new Map(entryGroups.map((entryGroup) => [getEntryCaption(entryGroup), entryGroup]));
      const orderedGroups = captionOrder.map((caption) => groupByCaption.get(caption)).filter(Boolean);
      if (orderedGroups.length < 2) continue;

      const positionedEntries = entryGroups
        .map((entryGroup) => ({ entryGroup, anchor: getEntryAnchor(entryGroup) }))
        .sort((a, b) => {
          const yDelta = a.anchor.y - b.anchor.y;
          if (Math.abs(yDelta) < 1) return a.anchor.x - b.anchor.x;
          return yDelta;
        });

      const slots = positionedEntries.map(({ anchor }) => anchor);
      orderedGroups.forEach((entryGroup, idx) => {
        const targetAnchor = slots[idx];
        if (!targetAnchor) return;

        const currentAnchor = getEntryAnchor(entryGroup);
        const deltaX = targetAnchor.x - currentAnchor.x;
        const deltaY = targetAnchor.y - currentAnchor.y;
        if (Math.abs(deltaX) >= 1e-6 || Math.abs(deltaY) >= 1e-6) {
          translateEntryGroupChildren(entryGroup, deltaX, deltaY);
          changed = true;
        }
      });

      const remainingGroups = entryGroups.filter((entryGroup) => !orderedGroups.includes(entryGroup));
      [...orderedGroups, ...remainingGroups].forEach((entryGroup) => targetGroup.appendChild(entryGroup));
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

    persistLegendSvg(svg);
  };

  const getVisibleLegendOrder = () => legendEntries.value.map((entry) => entry.caption).filter(Boolean);

  const moveLegendEntryUp = (idx) => {
    if (idx <= 0) return;
    swapLegendEntries(idx, idx - 1);
  };

  const moveLegendEntryDown = (idx) => {
    if (idx >= legendEntries.value.length - 1) return;
    swapLegendEntries(idx, idx + 1);
  };

  const sortLegendEntries = (direction = 'asc') => {
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
    const currentOrder = getVisibleLegendOrder();
    if (currentOrder.length < 2) return;
    if (originalLegendOrder.value.length === 0) return;

    const currentIndex = new Map(currentOrder.map((caption, idx) => [caption, idx]));
    const sortedOrder = [...currentOrder].sort((a, b) => {
      const aOrigIdx = originalLegendOrder.value.indexOf(a);
      const bOrigIdx = originalLegendOrder.value.indexOf(b);

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
