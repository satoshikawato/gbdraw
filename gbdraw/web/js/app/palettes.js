import { COMPARISON_COLOR_KEYS } from './color-utils.js';

export const createPaletteLoader = ({ state }) => {
  const {
    paletteDefinitions,
    paletteNames,
    selectedPalette,
    currentColors,
    appliedPaletteName,
    appliedPaletteColors,
    pendingPaletteName,
    pendingPaletteColors,
    normalizePaletteColors,
    normalizePaletteDefinitions
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

  const applyPalettes = (allPalettes) => {
    if (!allPalettes || typeof allPalettes !== 'object') return false;
    const normalizedPalettes = normalizePaletteDefinitions(allPalettes);
    if (Object.keys(normalizedPalettes).length === 0) return false;
    paletteDefinitions.value = normalizedPalettes;
    paletteNames.value = Object.keys(normalizedPalettes).sort();
    const requestedPalette = String(selectedPalette.value || 'default').trim() || 'default';
    const resolvedPalette = normalizedPalettes[requestedPalette] ? requestedPalette : 'default';
    const resolvedColors = normalizePaletteColors(
      normalizedPalettes[resolvedPalette] || normalizedPalettes.default || {}
    );
    const currentHasPaletteColors = hasPaletteColorEntries(currentColors.value);
    const appliedHasColors = hasPaletteColorEntries(appliedPaletteColors.value);
    const pendingHasColors = hasPaletteColorEntries(pendingPaletteColors.value);

    selectedPalette.value = resolvedPalette;
    if (!currentHasPaletteColors) currentColors.value = resolvedColors;
    if (!appliedHasColors) {
      appliedPaletteName.value = resolvedPalette;
      appliedPaletteColors.value = {
        ...(currentHasPaletteColors ? currentColors.value : resolvedColors)
      };
    }
    if (String(pendingPaletteName.value || '').trim() && !pendingHasColors) {
      pendingPaletteColors.value = {
        ...(currentHasPaletteColors ? currentColors.value : resolvedColors)
      };
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
