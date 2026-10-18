// @ts-check
import {
  COMPARISON_COLOR_KEYS,
  normalizePaletteColors,
  normalizePaletteDefinitions
} from '../utils/color-utils.js';
import { queuePalette, selectPalette } from '../services/drawing/palette.js';
/** @import { PaletteDrawing } from '../services/drawing/palette.js' */

/**
 * The palette fields of the Web state (`state.js`) the loader reads and writes.
 * Each `{ value }` is a Vue ref.
 * @typedef {object} PaletteLoaderState
 * @property {{ value: Record<string, Record<string, string>> | null }} paletteDefinitions
 * @property {{ value: string[] }} paletteNames
 * @property {{ value: string }} appliedPaletteName
 * @property {{ value: Record<string, string> }} appliedPaletteColors
 * @property {() => PaletteLoaderDrawing} activeDrawing The shown mode's drawing.
 * @property {() => readonly PaletteLoaderDrawing[]} drawingList Every drawing, in creation order.
 */

/** @typedef {PaletteDrawing} PaletteLoaderDrawing The palette of a drawing (`DrawingState` of state.js). */

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

  /**
   * @param {Record<string, string> | null | undefined} colors
   * @returns {colors is Record<string, string>} True for a non-empty plain object.
   */
  const hasColorEntries = (colors) => (
    colors != null &&
    typeof colors === 'object' &&
    !Array.isArray(colors) &&
    Object.keys(colors).length > 0
  );
  const comparisonColorKeys = new Set(COMPARISON_COLOR_KEYS);
  /** @param {Record<string, string> | null | undefined} colors */
  const hasPaletteColorEntries = (colors) => (
    hasColorEntries(colors) &&
    Object.keys(colors).some((key) => !comparisonColorKeys.has(key))
  );

  /**
   * Gives one drawing its palette's colors where it has none yet.
   * @param {PaletteLoaderDrawing} drawing
   * @param {Record<string, Record<string, string>>} normalizedPalettes
   */
  const initializeDrawingPalette = (drawing, normalizedPalettes) => {
    const requestedPalette = String(drawing.selectedPalette.value || 'default').trim() || 'default';
    const resolvedPalette = Object.hasOwn(normalizedPalettes, requestedPalette) && normalizedPalettes[requestedPalette]
      ? requestedPalette : 'default';
    const resolvedColors = normalizePaletteColors(
      normalizedPalettes[resolvedPalette] || normalizedPalettes.default || {}
    );
    selectPalette(drawing.slices.palette, resolvedPalette, hasPaletteColorEntries(drawing.currentColors.value) ? undefined : resolvedColors);
    if (String(drawing.pendingPaletteName.value || '').trim() && !hasPaletteColorEntries(drawing.pendingPaletteColors.value)) {
      queuePalette(drawing.slices.palette, drawing.pendingPaletteName.value, drawing.currentColors.value);
    }
  };

  // Every drawing gets its palette (PR-1: each mode keeps its own colors); the
  // applied palette, which the committed Result's Legend shows, starts from the
  // shown drawing's.
  /** @param {Record<string, Record<string, string> | string> | null | undefined} allPalettes The palette file: palette tables and a `title`. */
  const applyPalettes = (allPalettes) => {
    if (!allPalettes || typeof allPalettes !== 'object') return false;
    const normalizedPalettes = normalizePaletteDefinitions(allPalettes);
    if (Object.keys(normalizedPalettes).length === 0) return false;
    paletteDefinitions.value = normalizedPalettes;
    paletteNames.value = Object.keys(normalizedPalettes).sort();
    state.drawingList().forEach((drawing) => initializeDrawingPalette(drawing, normalizedPalettes));
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
