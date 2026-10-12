// @ts-check
import { resolveColorToHex } from '../utils/color-utils.js';
import { normalizeOptionalText } from '../services/track-slot-display.js';

/**
 * Where an inherited skew color comes from. The color and palette sources may
 * be Vue refs or plain values (`unwrapRef`).
 * @typedef {object} InheritedSkewSlotColorOptions
 * @property {string} [key] The slot param: `positive_color` or `negative_color`.
 * @property {MaybeRef<Record<string, unknown> | null | undefined>} [currentColors]
 * @property {MaybeRef<Record<string, Record<string, unknown> | undefined> | null | undefined>} [paletteDefinitions]
 * @property {unknown} [selectedPalette]
 */

/**
 * @typedef {InheritedSkewSlotColorOptions & { slot?: { params?: Record<string, unknown> } }} TrackSlotSkewColorOptions
 */

/**
 * @template T
 * @typedef {T | { value: T }} MaybeRef A plain value or a Vue ref holding it.
 */

const SKEW_COLOR_PARAM_TO_PALETTE_KEY = Object.freeze({
  positive_color: 'skew_high',
  negative_color: 'skew_low'
});

const FINAL_SKEW_COLOR_FALLBACKS = Object.freeze({
  skew_high: '#6dded3',
  skew_low: '#ad72e3'
});

/**
 * @template T
 * @param {MaybeRef<T>} value
 * @returns {T}
 */
const unwrapRef = (value) => {
  if (value && typeof value === 'object' && 'value' in value) {
    return value.value;
  }
  return value;
};

/** @param {unknown} value */
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

/** @param {unknown} selectedPalette */
const resolvePaletteName = (selectedPalette) => {
  const text = String(unwrapRef(selectedPalette) || 'default').trim();
  return text || 'default';
};

/** @param {InheritedSkewSlotColorOptions} [options] */
export const resolveInheritedSkewSlotColor = ({
  key,
  currentColors = {},
  paletteDefinitions = {},
  selectedPalette = 'default'
} = {}) => {
  if (!key) return '#777777';
  // The callers pass the two slot params named in the table; an unnamed key reads undefined below.
  const paletteKey = SKEW_COLOR_PARAM_TO_PALETTE_KEY[/** @type {keyof typeof SKEW_COLOR_PARAM_TO_PALETTE_KEY} */ (key)];
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

/** @param {TrackSlotSkewColorOptions} [options] */
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
  const explicit = normalizeColorInputValue(key === undefined ? undefined : slot?.params?.[key]);
  return explicit || inherited;
};
