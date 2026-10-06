// @ts-check
import { parseTransform } from '../legend-layout/transform-utils.js';

export { parseTransform };

export const getLegendChildById = (parent, id) => {
  if (!parent) return null;
  for (const child of parent.children || []) {
    if (child.id === id) return child;
    if (
      id === 'pairwise_legend' &&
      (
        child.getAttribute?.('data-gbdraw-role') === 'comparison-legend' ||
        /^pairwise_legend_(?:h|v)$/.test(String(child.id || ''))
      )
    ) {
      return child;
    }
  }
  return null;
};

export const PAIRWISE_LEGEND_SELECTOR =
  '[data-gbdraw-role="comparison-legend"][data-gbdraw-orientation="h"], [data-gbdraw-role="comparison-legend"][data-gbdraw-orientation="v"], #pairwise_legend, #pairwise_legend_h, #pairwise_legend_v';

export const COMPARISON_LEGEND_SELECTOR =
  `[data-gbdraw-role="comparison-legend"], ${PAIRWISE_LEGEND_SELECTOR}, #conservation_identity_legend`;

export const getComparisonLegendGroup = (parent) =>
  getLegendChildById(parent, 'pairwise_legend') ||
  getLegendChildById(parent, 'conservation_identity_legend');

export const isInsideComparisonLegend = (el) =>
  Boolean(el?.closest?.(COMPARISON_LEGEND_SELECTOR));

export const parseTransformXY = (transform) => {
  if (!transform) return { x: 0, y: 0 };
  const number = '([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)';
  const match = String(transform).match(
    new RegExp(`translate\\(\\s*${number}(?:\\s*,\\s*|\\s+)${number}\\s*\\)`)
  );
  return match ? { x: parseFloat(match[1]), y: parseFloat(match[2]) } : { x: 0, y: 0 };
};

export const getAllFeatureLegendGroups = (svg) => {
  if (!svg) return [];

  const legendGroup = svg.getElementById('legend');
  if (!legendGroup) return [];

  const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
  const verticalLegend = legendGroup.querySelector('#legend_vertical');

  const groups = [];

  if (horizontalLegend && verticalLegend) {
    const hFeatureLegend = horizontalLegend.querySelector('#feature_legend_h');
    const vFeatureLegend = verticalLegend.querySelector('#feature_legend_v');
    if (hFeatureLegend) groups.push(hFeatureLegend);
    if (vFeatureLegend) groups.push(vFeatureLegend);
    return groups;
  }

  const featureLegendGroup = legendGroup.querySelector('#feature_legend');
  return featureLegendGroup ? [featureLegendGroup] : [legendGroup];
};

export const getVisibleFeatureLegendGroup = (svg) => {
  if (!svg) return null;

  const legendGroup = svg.getElementById('legend');
  if (!legendGroup) return null;

  const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
  const verticalLegend = legendGroup.querySelector('#legend_vertical');

  if (horizontalLegend && verticalLegend) {
    const hVisible =
      horizontalLegend.style.display !== 'none' &&
      horizontalLegend.getAttribute('display') !== 'none';
    const targetLegend = hVisible ? horizontalLegend : verticalLegend;
    const featureGroup = hVisible
      ? targetLegend.querySelector('#feature_legend_h')
      : targetLegend.querySelector('#feature_legend_v');
    return featureGroup || targetLegend;
  }

  const featureLegendGroup = legendGroup.querySelector('#feature_legend');
  return featureLegendGroup || legendGroup;
};

export const isCurrentLegendHorizontal = (svg) => {
  const legendGroup = svg?.getElementById('legend');
  if (!legendGroup) return false;

  const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
  const verticalLegend = legendGroup.querySelector('#legend_vertical');

  if (horizontalLegend && verticalLegend) {
    const hVisible =
      horizontalLegend.style.display !== 'none' &&
      horizontalLegend.getAttribute('display') !== 'none';
    return hVisible;
  }

  return false;
};

export const getLegendEntrySwatch = (entryGroup) => Array.from(
  entryGroup?.querySelectorAll?.('path') || []
).find((path) => {
  const fill = path.getAttribute('fill');
  return fill && fill !== 'none' && !fill.startsWith('url(');
}) || null;

/** The anchor of a Legend entry: its text position, else its swatch position. */
export const legendEntryAnchor = (entryGroup) => {
  const groupOffset = parseTransformXY(entryGroup?.getAttribute?.('transform'));
  const target = entryGroup?.querySelector?.('text') || getLegendEntrySwatch(entryGroup);
  const targetOffset = parseTransformXY(target?.getAttribute?.('transform'));
  return { x: groupOffset.x + targetOffset.x, y: groupOffset.y + targetOffset.y };
};

export const moveLegendEntryToAnchor = (entryGroup, xPos, yPos) => {
  if (!Number.isFinite(xPos) || !Number.isFinite(yPos)) return false;
  const current = legendEntryAnchor(entryGroup);
  const deltaX = xPos - current.x;
  const deltaY = yPos - current.y;
  if (Math.abs(deltaX) < 1e-6 && Math.abs(deltaY) < 1e-6) return false;
  if (entryGroup.hasAttribute?.('transform')) {
    const groupOffset = parseTransformXY(entryGroup.getAttribute('transform'));
    entryGroup.setAttribute(
      'transform',
      `translate(${groupOffset.x + deltaX}, ${groupOffset.y + deltaY})`
    );
    return true;
  }
  const transformedChildren = Array.from(entryGroup.querySelectorAll?.('[transform]') || []);
  if (transformedChildren.length === 0) {
    entryGroup.setAttribute('transform', `translate(${deltaX}, ${deltaY})`);
    return true;
  }
  transformedChildren.forEach((node) => {
    const position = parseTransformXY(node.getAttribute('transform'));
    node.setAttribute('transform', `translate(${position.x + deltaX}, ${position.y + deltaY})`);
  });
  return true;
};

/**
 * Legend entries in the default order (D-08): the generated entries ranked by
 * the generated order (a renamed entry ranks by its generated caption),
 * followed by the other entries in their current order. Without a generated
 * order, null.
 */
export const defaultLegendEntryOrder = (entries, generatedOrder) => {
  const rank = new Map();
  (Array.isArray(generatedOrder) ? generatedOrder : []).forEach((caption) => {
    const key = String(caption ?? '').trim();
    if (key && !rank.has(key)) rank.set(key, rank.size);
  });
  if (rank.size === 0) return null;
  const list = Array.isArray(entries) ? entries : [];
  const generated = (entry) => String(entry?.originalCaption || entry?.caption || '').trim();
  return [
    ...list.filter((entry) => rank.has(generated(entry)))
      .sort((left, right) => rank.get(generated(left)) - rank.get(generated(right))),
    ...list.filter((entry) => !rank.has(generated(entry)))
  ];
};

/**
 * Whether Legend entries are in an edited order, not the default order (D-08).
 * Generate replays an edited order; the extraction of a Generate keeps the
 * generated order only while an edited order is replayed.
 */
export const isLegendOrderEdited = (entries, generatedOrder) => {
  const defaultOrder = defaultLegendEntryOrder(entries, generatedOrder);
  return Boolean(defaultOrder) && defaultOrder.some((entry, index) => entry !== entries[index]);
};

const legendEntryGroups = (targetGroup) => {
  const direct = Array.from(targetGroup?.children || []).filter(
    (child) => String(child.tagName || '').toLowerCase() === 'g' && child.hasAttribute?.('data-legend-key')
  );
  return direct.length > 0
    ? direct
    : Array.from(targetGroup?.querySelectorAll?.('g[data-legend-key]') || []);
};

const legendEntryCaption = (entryGroup) => String(entryGroup.getAttribute('data-legend-key') || '').trim();

// One Legend group's entries in slot order (row-major anchor order).
const legendEntriesInSlotOrder = (targetGroup) => legendEntryGroups(targetGroup)
  .map((entryGroup, index) => ({ entryGroup, index, anchor: legendEntryAnchor(entryGroup) }))
  .sort((left, right) => {
    const yDelta = left.anchor.y - right.anchor.y;
    if (Math.abs(yDelta) >= 1) return yDelta;
    return left.anchor.x - right.anchor.x || left.index - right.index;
  });

// Whether the group's entries listed in the order already follow it.
const followsOrder = (positioned, rank) => {
  const listed = positioned
    .map(({ entryGroup }) => rank.get(legendEntryCaption(entryGroup)))
    .filter((value) => value !== undefined);
  return listed.every((value, index) => index === 0 || listed[index - 1] < value);
};

/**
 * Put one Legend group's entries in the requested caption order. Entries keep
 * the group's current slots (row-major anchor order) and receive them in the
 * new order; captions missing from the order follow in their current order.
 * Live Sort and Move and the Generate replay of the order (D-08) share this.
 * With `keepFollowed`, a group whose listed entries already follow the order
 * keeps its order, so the entries the order does not list keep their places:
 * a displayed batch Result draws entries that the order, read from another
 * Result, does not list (D-07, B18).
 * Returns whether any entry moved or changed document order.
 */
export const orderLegendEntries = (targetGroup, captionOrder, { keepFollowed = false } = {}) => {
  const entryGroups = legendEntryGroups(targetGroup);
  if (entryGroups.length < 2) return false;
  const rank = new Map();
  (Array.isArray(captionOrder) ? captionOrder : []).forEach((caption) => {
    const key = String(caption ?? '').trim();
    if (key && !rank.has(key)) rank.set(key, rank.size);
  });
  const caption = legendEntryCaption;
  const positioned = legendEntriesInSlotOrder(targetGroup);
  if (keepFollowed && followsOrder(positioned, rank)) return false;
  const ordered = positioned
    .map((item, slot) => ({ ...item, slot }))
    .sort((left, right) => {
      const leftRank = rank.has(caption(left.entryGroup)) ? rank.get(caption(left.entryGroup)) : Infinity;
      const rightRank = rank.has(caption(right.entryGroup)) ? rank.get(caption(right.entryGroup)) : Infinity;
      return leftRank === rightRank ? left.slot - right.slot : leftRank - rightRank;
    });
  let changed = false;
  ordered.forEach(({ entryGroup }, slot) => {
    const { x, y } = positioned[slot].anchor;
    changed = moveLegendEntryToAnchor(entryGroup, x, y) || changed;
  });
  ordered.forEach(({ entryGroup }, index) => {
    if (entryGroups[index] !== entryGroup) changed = true;
    targetGroup.appendChild(entryGroup);
  });
  return changed;
};
