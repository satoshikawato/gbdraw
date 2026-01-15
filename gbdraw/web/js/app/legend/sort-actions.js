import { getAllFeatureLegendGroups } from './utils.js';

export const createLegendSortActions = ({ state, extractLegendEntries }) => {
  const { svgContainer, legendEntries, originalLegendOrder, selectedResultIndex, results, skipCaptureBaseConfig } =
    state;

  const moveLegendEntryUp = (idx) => {
    if (idx <= 0) return;
    swapLegendEntries(idx, idx - 1);
  };

  const moveLegendEntryDown = (idx) => {
    if (idx >= legendEntries.value.length - 1) return;
    swapLegendEntries(idx, idx + 1);
  };

  const sortLegendEntries = (direction = 'asc') => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const entries = [...legendEntries.value];
    if (entries.length < 2) return;

    const yPositions = entries.map((e) => e.yPos).sort((a, b) => a - b);

    const sortedEntries = [...entries].sort((a, b) => {
      const cmp = a.caption.localeCompare(b.caption, undefined, { sensitivity: 'base' });
      return direction === 'asc' ? cmp : -cmp;
    });

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return;

    const yMapping = new Map();
    entries.forEach((entry) => {
      const sortedIdx = sortedEntries.findIndex((e) => e.caption === entry.caption && e.color === entry.color);
      if (sortedIdx !== -1) {
        yMapping.set(entry.yPos, yPositions[sortedIdx]);
      }
    });

    for (const targetGroup of targetGroups) {
      const allElements = targetGroup.querySelectorAll('path, text');
      allElements.forEach((el) => {
        const transform = el.getAttribute('transform');
        if (!transform) return;

        const match = transform.match(/translate\(([^,]+),\s*([\d.-]+)\)/);
        if (!match) return;

        const x = match[1];
        const y = parseFloat(match[2]);

        for (const [oldY, newY] of yMapping) {
          if (Math.abs(y - oldY) < 1) {
            el.setAttribute('transform', `translate(${x}, ${newY})`);
            break;
          }
        }
      });
    }

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }

    extractLegendEntries();
  };

  const sortLegendEntriesByDefault = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const entries = [...legendEntries.value];
    if (entries.length < 2) return;
    if (originalLegendOrder.value.length === 0) return;

    const yPositions = entries.map((e) => e.yPos).sort((a, b) => a - b);

    const sortedEntries = [...entries].sort((a, b) => {
      const aOrigIdx = originalLegendOrder.value.indexOf(a.caption);
      const bOrigIdx = originalLegendOrder.value.indexOf(b.caption);

      if (aOrigIdx !== -1 && bOrigIdx !== -1) {
        return aOrigIdx - bOrigIdx;
      }
      if (aOrigIdx !== -1) return -1;
      if (bOrigIdx !== -1) return 1;
      return a.caption.localeCompare(b.caption, undefined, { sensitivity: 'base' });
    });

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return;

    const yMapping = new Map();
    entries.forEach((entry) => {
      const sortedIdx = sortedEntries.findIndex((e) => e.caption === entry.caption && e.color === entry.color);
      if (sortedIdx !== -1) {
        yMapping.set(entry.yPos, yPositions[sortedIdx]);
      }
    });

    for (const targetGroup of targetGroups) {
      const allElements = targetGroup.querySelectorAll('path, text');
      allElements.forEach((el) => {
        const transform = el.getAttribute('transform');
        if (!transform) return;

        const match = transform.match(/translate\(([^,]+),\s*([\d.-]+)\)/);
        if (!match) return;

        const x = match[1];
        const y = parseFloat(match[2]);

        for (const [oldY, newY] of yMapping) {
          if (Math.abs(y - oldY) < 1) {
            el.setAttribute('transform', `translate(${x}, ${newY})`);
            break;
          }
        }
      });
    }

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }

    extractLegendEntries();
  };

  const swapLegendEntries = (idx1, idx2) => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return;

    const entry1 = legendEntries.value[idx1];
    const entry2 = legendEntries.value[idx2];
    if (!entry1 || !entry2) return;

    const y1 = entry1.yPos;
    const y2 = entry2.yPos;

    for (const targetGroup of targetGroups) {
      const allElements = targetGroup.querySelectorAll('path, text');
      allElements.forEach((el) => {
        const transform = el.getAttribute('transform');
        if (!transform) return;

        const match = transform.match(/translate\(([^,]+),\s*([\d.-]+)\)/);
        if (!match) return;

        const x = match[1];
        const y = parseFloat(match[2]);

        if (Math.abs(y - y1) < 1) {
          el.setAttribute('transform', `translate(${x}, ${y2})`);
        } else if (Math.abs(y - y2) < 1) {
          el.setAttribute('transform', `translate(${x}, ${y1})`);
        }
      });
    }

    legendEntries.value[idx1].yPos = y2;
    legendEntries.value[idx2].yPos = y1;

    const temp = { ...legendEntries.value[idx1] };
    legendEntries.value[idx1] = { ...legendEntries.value[idx2], yPos: y1 };
    legendEntries.value[idx2] = { ...temp, yPos: y2 };

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }

    extractLegendEntries();
  };

  return {
    moveLegendEntryDown,
    moveLegendEntryUp,
    sortLegendEntries,
    sortLegendEntriesByDefault,
    swapLegendEntries
  };
};
