export const parseTransform = (transformStr) => {
  if (!transformStr) return { x: 0, y: 0 };
  const match = transformStr.match(/translate\(\s*([-\d.]+)\s*,?\s*([-\d.]+)?\s*\)/);
  if (match) {
    return { x: parseFloat(match[1]) || 0, y: parseFloat(match[2]) || 0 };
  }
  return { x: 0, y: 0 };
};

export const getTransformedBBox = (el) => {
  if (!el) return null;
  const bbox = el.getBBox();
  const offset = parseTransform(el.getAttribute('transform'));
  return { x: bbox.x + offset.x, y: bbox.y + offset.y, width: bbox.width, height: bbox.height };
};

export const getElementsBounds = (elements) => {
  let minX = Infinity;
  let minY = Infinity;
  let maxX = -Infinity;
  let maxY = -Infinity;
  elements.forEach((el) => {
    if (!el) return;
    const bounds = getTransformedBBox(el);
    if (!bounds) return;
    minX = Math.min(minX, bounds.x);
    minY = Math.min(minY, bounds.y);
    maxX = Math.max(maxX, bounds.x + bounds.width);
    maxY = Math.max(maxY, bounds.y + bounds.height);
  });
  if (minX === Infinity) return null;
  return { x: minX, y: minY, width: maxX - minX, height: maxY - minY };
};

export const getLocalVerticalBounds = (group) => {
  let minY = Number.POSITIVE_INFINITY;
  let maxY = Number.NEGATIVE_INFINITY;
  const texts = Array.from(group?.children || [])
    .filter((child) => String(child?.tagName || '').toLowerCase() === 'text');
  texts.forEach((textEl) => {
    const yValue = Number(String(textEl.getAttribute('y') || '0').replace('px', ''));
    const fontSize = Number(String(textEl.getAttribute('font-size') || '0').replace('px', ''));
    const halfHeight = Number.isFinite(fontSize) && fontSize > 0 ? 0.5 * fontSize : 0;
    minY = Math.min(minY, yValue - halfHeight);
    maxY = Math.max(maxY, yValue + halfHeight);
  });
  if (!Number.isFinite(minY) || !Number.isFinite(maxY)) {
    return { minY: 0, maxY: 0, height: 0 };
  }
  return { minY, maxY, height: maxY - minY };
};

export const getDefinitionGroupTranslate = (group, canvasWidth, canvasHeight, position, edgeMargin = 24) => {
  const normalizedPosition = ['center', 'top', 'bottom'].includes(position) ? position : 'center';
  const bounds = getLocalVerticalBounds(group);
  let y = 0.5 * canvasHeight;
  if (normalizedPosition === 'top') {
    y = edgeMargin - bounds.minY;
  } else if (normalizedPosition === 'bottom') {
    y = canvasHeight - edgeMargin - bounds.maxY;
  }
  return { x: 0.5 * canvasWidth, y };
};
