const PAIRWISE_LEGEND_SELECTOR =
  '[data-gbdraw-role="comparison-legend"][data-gbdraw-orientation="h"], '
  + '[data-gbdraw-role="comparison-legend"][data-gbdraw-orientation="v"], '
  + '#pairwise_legend, #pairwise_legend_h, #pairwise_legend_v';

export const getGroupsByBaseIds = (svg, baseIds, slotRenderers = []) => {
  if (!svg) return [];
  const seen = new Set();
  const groups = [];
  const addGroups = (selector) => {
    svg.querySelectorAll(selector).forEach((group) => {
      if (seen.has(group)) return;
      seen.add(group);
      groups.push(group);
    });
  };
  (Array.isArray(slotRenderers) ? slotRenderers : []).forEach((renderer) => {
    if (renderer) addGroups(`g[data-gbdraw-slot-renderer="${renderer}"]`);
  });
  if (groups.length > 0) return groups;
  (Array.isArray(baseIds) ? baseIds : []).forEach((baseId) => {
    if (baseId) addGroups(`g[id="${baseId}"], g[id^="${baseId}_"]`);
  });
  return groups;
};

const parseTranslate = (value) => {
  const number = '([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)';
  const match = String(value || '').match(
    new RegExp(`translate\\(\\s*${number}(?:\\s*,\\s*|\\s+)${number}\\s*\\)`)
  );
  return match ? { x: Number(match[1]), y: Number(match[2]) } : { x: 0, y: 0 };
};

const featureLegendGroups = (svg) => {
  const legend = svg?.getElementById?.('legend');
  if (!legend) return [];
  const horizontal = legend.querySelector('#legend_horizontal');
  const vertical = legend.querySelector('#legend_vertical');
  if (horizontal && vertical) {
    return [
      horizontal.querySelector('#feature_legend_h') || horizontal,
      vertical.querySelector('#feature_legend_v') || vertical
    ];
  }
  return [legend.querySelector('#feature_legend') || legend];
};

export const normalizeLegacyLegendEntryGroups = (svg) => {
  const documentOwner = svg?.ownerDocument || globalThis.document;
  if (!documentOwner?.createElementNS) return false;
  let changed = false;
  featureLegendGroups(svg).forEach((targetGroup) => {
    if (targetGroup.querySelectorAll('g[data-legend-key]').length > 0) return;
    const texts = Array.from(targetGroup.querySelectorAll('text'));
    const paths = Array.from(targetGroup.querySelectorAll('path'));
    const additions = [];
    texts.forEach((textElement) => {
      const caption = textElement.textContent?.trim();
      if (!caption) return;
      const textPosition = parseTranslate(textElement.getAttribute('transform'));
      let bestPath = null;
      let bestX = -Infinity;
      paths.forEach((path) => {
        const fill = path.getAttribute('fill');
        if (!fill || fill === 'none' || fill.startsWith('url(')) return;
        const position = parseTranslate(path.getAttribute('transform'));
        if (
          Math.abs(position.y - textPosition.y) < 2
          && position.x < textPosition.x
          && position.x > bestX
        ) {
          bestPath = path;
          bestX = position.x;
        }
      });
      if (!bestPath) return;
      const group = documentOwner.createElementNS('http://www.w3.org/2000/svg', 'g');
      group.setAttribute('data-legend-key', caption);
      group.appendChild(bestPath);
      group.appendChild(textElement);
      additions.push(group);
    });
    additions.forEach((group) => targetGroup.appendChild(group));
    changed = additions.length > 0 || changed;
  });
  return changed;
};

export const ensureUniqueSkewClipPathIds = (svg) => {
  if (!svg) return false;
  const skewGroups = getGroupsByBaseIds(
    svg,
    ['skew', 'gc_skew'],
    ['dinucleotide_skew']
  );
  if (skewGroups.length === 0) return false;
  let changed = false;
  const idCounts = new Map();
  const usedIds = new Set();
  svg.querySelectorAll('[id]').forEach((element) => {
    const id = element.getAttribute('id');
    if (!id) return;
    idCounts.set(id, (idCounts.get(id) || 0) + 1);
    usedIds.add(id);
  });

  skewGroups.forEach((skewGroup, groupIndex) => {
    const groupId = skewGroup.getAttribute('id') || `skew_group_${groupIndex}`;
    skewGroup.querySelectorAll('clipPath[id]').forEach((clipPath, clipIndex) => {
      const oldId = clipPath.getAttribute('id');
      if (!oldId || (idCounts.get(oldId) || 0) <= 1) return;
      const baseNewId = `${oldId}_${groupId}_${clipIndex}`;
      let newId = baseNewId;
      let suffix = 1;
      while (usedIds.has(newId) && newId !== oldId) {
        newId = `${baseNewId}_${suffix}`;
        suffix += 1;
      }
      if (newId === oldId) return;
      clipPath.setAttribute('id', newId);
      usedIds.add(newId);
      changed = true;
      const oldReference = `url(#${oldId})`;
      const newReference = `url(#${newId})`;
      skewGroup.querySelectorAll('[clip-path], [clipPath]').forEach((element) => {
        if (element.getAttribute('clip-path') === oldReference) {
          element.setAttribute('clip-path', newReference);
          changed = true;
        }
        if (element.getAttribute('clipPath') === oldReference) {
          element.setAttribute('clipPath', newReference);
          changed = true;
        }
      });
    });
  });
  return changed;
};

export const ensureUniquePairwiseGradientIds = (svg) => {
  if (!svg) return false;
  const legendGroup = svg.getElementById?.('legend');
  const horizontalLegend = legendGroup?.querySelector('#legend_horizontal');
  const verticalLegend = legendGroup?.querySelector('#legend_vertical');
  if (!horizontalLegend || !verticalLegend) return false;
  let changed = false;
  const fixGradientId = (legend, suffix) => {
    const pairwiseLegend = legend.querySelector(PAIRWISE_LEGEND_SELECTOR);
    if (!pairwiseLegend) return;
    pairwiseLegend.querySelectorAll('linearGradient').forEach((gradient) => {
      const currentId = gradient.id;
      if (!currentId || currentId.endsWith(`_${suffix}`)) return;
      const newId = `${currentId.replace(/_[hv]$/, '')}_${suffix}`;
      gradient.setAttribute('id', newId);
      changed = true;
      pairwiseLegend.querySelectorAll('path').forEach((path) => {
        if (path.getAttribute('fill') === `url(#${currentId})`) {
          path.setAttribute('fill', `url(#${newId})`);
          changed = true;
        }
      });
    });
  };
  fixGradientId(horizontalLegend, 'h');
  fixGradientId(verticalLegend, 'v');
  return changed;
};

export const normalizeSvgResultIds = (svg) => {
  const skewChanged = ensureUniqueSkewClipPathIds(svg);
  const gradientChanged = ensureUniquePairwiseGradientIds(svg);
  return Boolean(skewChanged || gradientChanged);
};
