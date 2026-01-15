export const getLegendChildById = (parent, id) => {
  if (!parent) return null;
  for (const child of parent.children || []) {
    if (child.id === id) return child;
  }
  return null;
};

export const parseTransformXY = (transform) => {
  if (!transform) return { x: 0, y: 0 };
  const match = transform.match(/translate\(\s*([\d.-]+)\s*,\s*([\d.-]+)\s*\)/);
  return match ? { x: parseFloat(match[1]), y: parseFloat(match[2]) } : { x: 0, y: 0 };
};

export const parseTransform = (transformStr) => {
  if (!transformStr) return { x: 0, y: 0 };
  const match = transformStr.match(/translate\(\s*([-\d.]+)\s*,?\s*([-\d.]+)?\s*\)/);
  if (match) {
    return { x: parseFloat(match[1]) || 0, y: parseFloat(match[2]) || 0 };
  }
  return { x: 0, y: 0 };
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
