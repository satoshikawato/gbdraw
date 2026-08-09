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

export const COLLINEAR_ORIENTATION_COLOR_KEYS = Object.freeze({
  plus: 'collinear_block_plus',
  minus: 'collinear_block_minus'
});

export const COLLINEAR_ORIENTATION_MIN_COLOR_KEYS = Object.freeze({
  plus: 'collinear_block_plus_min',
  minus: 'collinear_block_minus_min'
});

export const DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS = Object.freeze({
  plus: '#f0f1f5',
  minus: '#FFE7E7'
});

export const DEFAULT_COLLINEAR_ORIENTATION_COLORS = Object.freeze({
  plus: '#8b9cc1',
  minus: '#E15759'
});

export const DEFAULT_COMPARISON_COLORS = Object.freeze({
  pairwise_match: '#d3d3d3',
  pairwise_match_min: '#FFE7E7',
  pairwise_match_max: '#FF7272',
  collinear_block_plus_min: DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS.plus,
  collinear_block_plus: DEFAULT_COLLINEAR_ORIENTATION_COLORS.plus,
  collinear_block_minus_min: DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS.minus,
  collinear_block_minus: DEFAULT_COLLINEAR_ORIENTATION_COLORS.minus
});

export const COMPARISON_COLOR_KEYS = Object.freeze(Object.keys(DEFAULT_COMPARISON_COLORS));

export const resolvePairwiseLegendGradientColorKeys = (legendKey) => {
  const normalizedKey = String(legendKey || '').trim();
  if (
    normalizedKey === 'Collinear' ||
    normalizedKey === 'Same direction' ||
    normalizedKey === 'Collinear identity'
  ) {
    return {
      minKey: COLLINEAR_ORIENTATION_MIN_COLOR_KEYS.plus,
      maxKey: COLLINEAR_ORIENTATION_COLOR_KEYS.plus
    };
  }
  if (normalizedKey === 'Inverted' || normalizedKey === 'Inverted identity') {
    return {
      minKey: COLLINEAR_ORIENTATION_MIN_COLOR_KEYS.minus,
      maxKey: COLLINEAR_ORIENTATION_COLOR_KEYS.minus
    };
  }
  return { minKey: 'pairwise_match_min', maxKey: 'pairwise_match_max' };
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

let namedColorContext;

const resolveBrowserNamedColor = (value) => {
  if (!/^[a-z]+$/i.test(value) || !globalThis.document?.createElement) return null;
  if (namedColorContext === undefined) {
    namedColorContext = globalThis.document.createElement('canvas').getContext?.('2d') || null;
  }
  if (!namedColorContext) return null;

  namedColorContext.fillStyle = '#010203';
  namedColorContext.fillStyle = value;
  const resolved = String(namedColorContext.fillStyle || '');
  if (resolved.toLowerCase() === '#010203') return null;
  if (/^#[0-9a-f]{6}$/i.test(resolved)) return resolved.toUpperCase();

  const rgb = resolved.match(/^rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)$/i);
  return rgb ? rgbToHex(Number(rgb[1]), Number(rgb[2]), Number(rgb[3])).toUpperCase() : null;
};

export const resolveColorToHex = (colorValue) => {
  if (!colorValue || typeof colorValue !== 'string') return colorValue;
  const trimmed = colorValue.trim();
  if (!trimmed) return trimmed;
  if (trimmed.startsWith('#')) return trimmed;
  return resolveBrowserNamedColor(trimmed) || trimmed;
};

export const colorValueMode = (colorValue) => {
  if (colorValue === null || colorValue === undefined || String(colorValue).trim() === '') {
    return 'auto';
  }
  return String(colorValue).trim().toLowerCase() === 'none' ? 'none' : 'color';
};

export const toNativeColorInputValue = (colorValue, fallback = '#000000') => {
  const resolved = String(resolveColorToHex(colorValue) || '').trim();
  const shortHex = resolved.match(/^#([0-9a-f]{3})$/i);
  if (shortHex) {
    return `#${shortHex[1].split('').map((value) => `${value}${value}`).join('')}`.toLowerCase();
  }
  if (/^#[0-9a-f]{6}$/i.test(resolved)) return resolved.toLowerCase();
  const normalizedFallback = String(resolveColorToHex(fallback) || '#000000').trim();
  return /^#[0-9a-f]{6}$/i.test(normalizedFallback)
    ? normalizedFallback.toLowerCase()
    : '#000000';
};

export const colorValueForMode = (mode, currentColor = null, fallback = '#000000') => {
  const normalizedMode = String(mode || '').trim().toLowerCase();
  if (normalizedMode === 'auto') return null;
  if (normalizedMode === 'none') return 'none';
  return toNativeColorInputValue(currentColor, fallback);
};

export const normalizePaletteColors = (colors = {}) => {
  const normalized = { ...(colors || {}) };
  Object.keys(normalized).forEach((key) => {
    if (/^collinear_block_\d+$/.test(key)) delete normalized[key];
  });
  if (normalized.collinear_block_plus_max && !normalized.collinear_block_plus) {
    normalized.collinear_block_plus = normalized.collinear_block_plus_max;
  }
  Object.entries(DEFAULT_COMPARISON_COLORS).forEach(([key, value]) => {
    if (!normalized[key]) normalized[key] = value;
  });
  return normalized;
};

export const normalizePaletteDefinitions = (palettes = {}) => {
  const normalized = {};
  Object.entries(palettes || {}).forEach(([name, colors]) => {
    if (name === 'title') return;
    normalized[name] = normalizePaletteColors(colors || {});
  });
  return normalized;
};

const normalizeComparableColor = (colorValue) => String(resolveColorToHex(String(colorValue || '').trim()) || '')
  .trim()
  .toLowerCase();

export const buildPaletteColorOverrideRows = ({
  colors = {},
  paletteColors = {}
} = {}) => {
  const rows = [];
  Object.entries(colors || {}).forEach(([key, color]) => {
    const normalizedKey = String(key || '').trim();
    const normalizedColor = String(color || '').trim();
    if (!normalizedKey || !normalizedColor) return;
    const paletteColor = paletteColors?.[normalizedKey];
    if (normalizeComparableColor(normalizedColor) === normalizeComparableColor(paletteColor)) return;
    rows.push([normalizedKey, normalizedColor]);
  });
  return rows;
};

export const buildDefaultColorOverrideTsv = ({
  colors = {},
  paletteColors = {}
} = {}) => buildPaletteColorOverrideRows({ colors, paletteColors })
  .map(([key, color]) => `${key}\t${color}`)
  .join('\n');

export const resolveCollinearMatchColor = ({
  blockId,
  colorMode,
  orientation,
  identityFactor,
  colors = {}
}) => {
  const normalizedBlockId = String(blockId || '').trim();
  if (!normalizedBlockId) return null;

  const normalizedMode = String(colorMode || '').trim().toLowerCase().replace(/-/g, '_');
  const normalizedOrientation = String(orientation || '').trim().toLowerCase();
  const colorKey = COLLINEAR_ORIENTATION_COLOR_KEYS[normalizedOrientation];
  const minColorKey = COLLINEAR_ORIENTATION_MIN_COLOR_KEYS[normalizedOrientation];
  if (!colorKey) return null;
  const orientationColor = colors[colorKey] || DEFAULT_COLLINEAR_ORIENTATION_COLORS[normalizedOrientation];

  if (normalizedMode === 'orientation') return orientationColor;
  if (normalizedMode === 'orientation_identity') {
    const factor = Number(identityFactor);
    if (!Number.isFinite(factor)) return null;
    const minColor = colors[minColorKey] || DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS[normalizedOrientation];
    return interpolateColor(
      minColor,
      orientationColor,
      factor
    );
  }

  return null;
};
