export const hexToRgb = (hex) => {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  return result
    ? {
        r: parseInt(result[1], 16),
        g: parseInt(result[2], 16),
        b: parseInt(result[3], 16)
      }
    : { r: 128, g: 128, b: 128 };
};

export const rgbToHex = (r, g, b) => {
  return (
    '#' +
    [r, g, b]
      .map((x) => {
        const hex = Math.round(Math.max(0, Math.min(255, x))).toString(16);
        return hex.length === 1 ? '0' + hex : hex;
      })
      .join('')
  );
};

export const interpolateColor = (color1, color2, factor) => {
  const c1 = hexToRgb(color1);
  const c2 = hexToRgb(color2);
  return rgbToHex(
    c1.r + (c2.r - c1.r) * factor,
    c1.g + (c2.g - c1.g) * factor,
    c1.b + (c2.b - c1.b) * factor
  );
};

export const estimateColorFactor = (currentColor, minColor, maxColor) => {
  const current = hexToRgb(currentColor);
  const min = hexToRgb(minColor);
  const max = hexToRgb(maxColor);

  const totalDist = Math.sqrt(
    Math.pow(max.r - min.r, 2) +
      Math.pow(max.g - min.g, 2) +
      Math.pow(max.b - min.b, 2)
  );
  if (totalDist < 1) return 0.5;

  const currentDist = Math.sqrt(
    Math.pow(current.r - min.r, 2) +
      Math.pow(current.g - min.g, 2) +
      Math.pow(current.b - min.b, 2)
  );
  return Math.max(0, Math.min(1, currentDist / totalDist));
};
