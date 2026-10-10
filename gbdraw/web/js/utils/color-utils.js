// @ts-check
import { namedColorHex } from './named-colors.js';
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

// A CSS color name resolves through the table Python shares (OV-160), so
// Load, the Session split and Python agree with or without a browser; a name
// outside that table, such as a system color, stays as written (OV-272).
export const resolveColorToHex = (colorValue) => {
  if (!colorValue || typeof colorValue !== 'string') return colorValue;
  const trimmed = colorValue.trim();
  if (!trimmed) return trimmed;
  if (trimmed.startsWith('#')) return trimmed;
  return namedColorHex(trimmed) || trimmed;
};

// Specific-color table domain, shared with Python's read_color_table:
// `none`, #RGB, #RRGGBB, or a color name of the table Python uses, resolved
// to hex. Names only a browser knows are rejected as in Python (OV-271).
export const normalizeSpecificRuleColor = (colorValue) => {
  const color = String(colorValue ?? '').trim().toLowerCase();
  if (color === 'none' || /^#(?:[0-9a-f]{3}|[0-9a-f]{6})$/.test(color)) return color;
  return namedColorHex(color)?.toLowerCase() || null;
};

// Default colors (-d) domain, shared with Python's is_user_color minus
// svgwrite's other paint values (OV-272, OV-302): `none` or `transparent`, a
// color name of the shared table resolved to hex, #RGB, #RGBA, #RRGGBB,
// #RRGGBBAA, or rgb()/rgba()/hsl()/hsla() with the arguments
// `_is_css_color_function` reads, both as written. Anything else is null.
const CSS_NUMBER = String.raw`[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?`;
const CSS_NUMBER_OR_PERCENT = new RegExp(`^${CSS_NUMBER}%?$`);
const CSS_HUE = new RegExp(`^${CSS_NUMBER}(?:deg|grad|rad|turn)?$`, 'i');
/**
 * @param {unknown} colorValue
 * @returns {string | null}
 */
export const normalizeDefaultColor = (colorValue) => {
  const color = String(colorValue ?? '').trim();
  const keyword = color.toLowerCase();
  if (keyword === 'none' || keyword === 'transparent') return keyword;
  if (/^#(?:[0-9a-f]{3,4}|[0-9a-f]{6}|[0-9a-f]{8})$/i.test(color)) return color;
  const fn = /^(rgba?|hsla?)\(([\s\S]*)\)$/i.exec(color);
  if (!fn) return namedColorHex(color) || null;
  const args = fn[2].trim();
  let channels;
  let alpha;
  if (args.includes(',')) {
    const parts = args.split(',').map((part) => part.trim());
    if (args.includes('/') || (parts.length !== 3 && parts.length !== 4)) return null;
    [channels, alpha] = [parts.slice(0, 3), parts.slice(3)];
  } else {
    const slash = args.indexOf('/');
    const words = (/** @type {string} */ text) => text.split(/\s+/).filter(Boolean);
    channels = words(slash < 0 ? args : args.slice(0, slash));
    alpha = slash < 0 ? [] : words(args.slice(slash + 1));
    if (channels.length !== 3 || (slash >= 0 && alpha.length !== 1)) return null;
  }
  const first = fn[1].toLowerCase().startsWith('hsl') ? CSS_HUE : CSS_NUMBER_OR_PERCENT;
  return first.test(channels[0]) && [...channels.slice(1), ...alpha].every((part) => CSS_NUMBER_OR_PERCENT.test(part))
    ? color
    : null;
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

/**
 * @param {Record<string, string>} [colors]
 * @returns {Record<string, string>}
 */
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

// OV-262, OV-263: the colors the shown Result draws, for every reader of a
// feature or Legend color: the applied palette's colors under the non-empty
// applied colors. A key that is Auto (empty), as a live edit, the Generate
// commit, Load, and a History restore leave it, shows the applied palette's
// color, as Generate draws a key without a `-d` row.
/**
 * @param {{
 *   paletteDefinitions?: { value: Record<string, Record<string, string>> | null },
 *   appliedPaletteName?: { value: string },
 *   appliedPaletteColors?: { value: Record<string, string | null> | null }
 * }} state
 * @returns {Record<string, string>}
 */
export const appliedFeatureColors = ({ paletteDefinitions, appliedPaletteName, appliedPaletteColors }) => {
  /** @type {Record<string, string>} */
  const colors = { ...(paletteDefinitions?.value?.[appliedPaletteName?.value ?? ''] || {}) };
  Object.entries(appliedPaletteColors?.value || {}).forEach(([key, color]) => {
    if (color && String(color).trim() !== '') colors[key] = color;
  });
  return colors;
};

// One color, compared: case, a color name, and `#rgb` against `#rrggbb` (D-15).
const normalizeComparableColor = (colorValue) => {
  const color = String(resolveColorToHex(String(colorValue || '').trim()) || '').trim().toLowerCase();
  return /^#[0-9a-f]{3}$/.test(color) ? color.replace(/[0-9a-f]/g, '$&$&') : color;
};

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
