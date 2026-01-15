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
