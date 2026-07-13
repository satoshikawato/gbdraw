import { resolveColorToHex } from './color-utils.js';
import { normalizeOptionalText } from './track-slot-display.js';

const SKEW_COLOR_PARAM_TO_PALETTE_KEY = Object.freeze({
  positive_color: 'skew_high',
  negative_color: 'skew_low'
});

const FINAL_SKEW_COLOR_FALLBACKS = Object.freeze({
  skew_high: '#6dded3',
  skew_low: '#ad72e3'
});

const unwrapRef = (value) => {
  if (value && typeof value === 'object' && Object.prototype.hasOwnProperty.call(value, 'value')) {
    return value.value;
  }
  return value;
};

export const normalizeColorInputValue = (value) => {
  const text = normalizeOptionalText(value);
  if (text === null) return null;
  const resolved = String(resolveColorToHex(text) || '').trim();
  const fullHex = resolved.match(/^#([0-9a-f]{6})$/i);
  if (fullHex) return `#${fullHex[1].toLowerCase()}`;
  const shortHex = resolved.match(/^#([0-9a-f])([0-9a-f])([0-9a-f])$/i);
  if (shortHex) {
    return `#${shortHex[1]}${shortHex[1]}${shortHex[2]}${shortHex[2]}${shortHex[3]}${shortHex[3]}`.toLowerCase();
  }
  return null;
};

const resolvePaletteName = (selectedPalette) => {
  const text = String(unwrapRef(selectedPalette) || 'default').trim();
  return text || 'default';
};

export const resolveInheritedSkewSlotColor = ({
  key,
  currentColors = {},
  paletteDefinitions = {},
  selectedPalette = 'default'
} = {}) => {
  const paletteKey = SKEW_COLOR_PARAM_TO_PALETTE_KEY[key];
  if (!paletteKey) return '#777777';

  const colors = unwrapRef(currentColors) || {};
  const palettes = unwrapRef(paletteDefinitions) || {};
  const selectedName = resolvePaletteName(selectedPalette);
  const candidates = [
    colors?.[paletteKey],
    palettes?.[selectedName]?.[paletteKey],
    palettes?.default?.[paletteKey],
    FINAL_SKEW_COLOR_FALLBACKS[paletteKey]
  ];

  for (const candidate of candidates) {
    const normalized = normalizeColorInputValue(candidate);
    if (normalized !== null) return normalized;
  }
  return '#777777';
};

export const resolveTrackSlotSkewColorValue = ({
  slot,
  key,
  currentColors = {},
  paletteDefinitions = {},
  selectedPalette = 'default'
} = {}) => {
  const inherited = resolveInheritedSkewSlotColor({
    key,
    currentColors,
    paletteDefinitions,
    selectedPalette
  });
  const explicit = normalizeColorInputValue(slot?.params?.[key]);
  return explicit || inherited;
};
