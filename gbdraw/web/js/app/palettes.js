// @ts-check
import {
  COMPARISON_COLOR_KEYS,
  normalizePaletteColors,
  normalizePaletteDefinitions
} from '../utils/color-utils.js';

/**
 * The palette fields of the Web state (`state.js`) the loader reads and writes.
 * Each `{ value }` is a Vue ref.
 * @typedef {object} PaletteLoaderState
 * @property {{ value: any }} paletteDefinitions
 * @property {{ value: string[] }} paletteNames
 * @property {{ value: string }} appliedPaletteName
 * @property {{ value: Record<string, string> }} appliedPaletteColors
 * @property {() => PaletteLoaderDrawing} activeDrawing The shown mode's drawing.
 * @property {{ circular: PaletteLoaderDrawing, linear: PaletteLoaderDrawing }} drawings
 */

/**
 * The palette members of a drawing (`DrawingState` of state.js).
 * @typedef {object} PaletteLoaderDrawing
 * @property {{ value: string }} selectedPalette
 * @property {{ value: Record<string, string> }} currentColors
 * @property {{ value: string }} pendingPaletteName
 * @property {{ value: Record<string, string> }} pendingPaletteColors
 */

/**
 * @typedef {object} PaletteLoaderOptions
 * @property {PaletteLoaderState} state
 */

/** @param {PaletteLoaderOptions} options */
export const createPaletteLoader = ({ state }) => {
  const {
    paletteDefinitions,
    paletteNames,
    appliedPaletteName,
    appliedPaletteColors
  } = state;

  const hasColorEntries = (colors) => (
    Boolean(colors) &&
    typeof colors === 'object' &&
    !Array.isArray(colors) &&
    Object.keys(colors).length > 0
  );
  const comparisonColorKeys = new Set(COMPARISON_COLOR_KEYS);
  const hasPaletteColorEntries = (colors) => (
    hasColorEntries(colors) &&
    Object.keys(colors).some((key) => !comparisonColorKeys.has(key))
  );

  /**
   * Gives one drawing its palette's colors where it has none yet.
   * @param {PaletteLoaderDrawing} drawing
   * @param {Record<string, any>} normalizedPalettes
   */
  const initializeDrawingPalette = (drawing, normalizedPalettes) => {
    const requestedPalette = String(drawing.selectedPalette.value || 'default').trim() || 'default';
    const resolvedPalette = normalizedPalettes[requestedPalette] ? requestedPalette : 'default';
    const resolvedColors = normalizePaletteColors(
      normalizedPalettes[resolvedPalette] || normalizedPalettes.default || {}
    );
    const currentHasPaletteColors = hasPaletteColorEntries(drawing.currentColors.value);
    drawing.selectedPalette.value = resolvedPalette;
    if (!currentHasPaletteColors) drawing.currentColors.value = resolvedColors;
    if (String(drawing.pendingPaletteName.value || '').trim() && !hasPaletteColorEntries(drawing.pendingPaletteColors.value)) {
      drawing.pendingPaletteColors.value = { ...drawing.currentColors.value };
    }
  };

  // Every drawing gets its palette (PR-1: each mode keeps its own colors); the
  // applied palette, which the committed Result's Legend shows, starts from the
  // shown drawing's.
  const applyPalettes = (allPalettes) => {
    if (!allPalettes || typeof allPalettes !== 'object') return false;
    const normalizedPalettes = normalizePaletteDefinitions(allPalettes);
    if (Object.keys(normalizedPalettes).length === 0) return false;
    paletteDefinitions.value = normalizedPalettes;
    paletteNames.value = Object.keys(normalizedPalettes).sort();
    initializeDrawingPalette(state.drawings.circular, normalizedPalettes);
    initializeDrawingPalette(state.drawings.linear, normalizedPalettes);
    if (!hasPaletteColorEntries(appliedPaletteColors.value)) {
      const shown = state.activeDrawing();
      appliedPaletteName.value = shown.selectedPalette.value;
      appliedPaletteColors.value = { ...shown.currentColors.value };
    }
    return true;
  };

  const loadPaletteAsset = async () => {
    const url = new URL('../../gallery/palettes/palettes.json', import.meta.url);
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`Palette asset request failed (${response.status}).`);
    }
    const payload = await response.json();
    if (!applyPalettes(payload?.palettes || payload)) {
      throw new Error('Palette asset is empty or malformed.');
    }
  };

  return { applyPalettes, loadPaletteAsset };
};
