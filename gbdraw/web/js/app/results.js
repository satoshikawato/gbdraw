// @ts-check
import { normalizePaletteColors } from '../utils/color-utils.js';

/**
 * The palette refs the manager reads and writes (state.js owns them).
 * @typedef {object} ResultsManagerState
 * @property {{ value: Record<string, Record<string, string>> | null }} paletteDefinitions
 * @property {{ value: string }} selectedPalette
 * @property {{ value: Record<string, string> }} currentColors
 * @property {{ value: boolean }} paletteInstantPreviewEnabled
 * @property {{ value: string }} appliedPaletteName
 * @property {{ value: Record<string, string> }} appliedPaletteColors
 * @property {{ value: string }} pendingPaletteName
 * @property {{ value: Record<string, string> }} pendingPaletteColors
 * @property {() => any} [sessionOperationAvailability] The busy outcome of a Session operation, if any.
 */

/** @param {{ state: ResultsManagerState }} options */
export const createResultsManager = ({ state }) => {
  const {
    paletteDefinitions,
    selectedPalette,
    currentColors,
    paletteInstantPreviewEnabled,
    appliedPaletteName,
    appliedPaletteColors,
    pendingPaletteName,
    pendingPaletteColors
  } = state;
  const cloneColors = (colors) => ({ ...(colors || {}) });
  const getPaletteMap = () => {
    if (paletteDefinitions.value && Object.keys(paletteDefinitions.value).length > 0) {
      return paletteDefinitions.value;
    }
    return {};
  };
  const getPaletteBaseColors = (paletteName) => {
    const allPalettes = getPaletteMap();
    return normalizePaletteColors(cloneColors(allPalettes[paletteName] || {}));
  };
  const setAppliedPaletteState = (paletteName, colors = currentColors.value) => {
    appliedPaletteName.value = String(paletteName || selectedPalette.value || 'default');
    appliedPaletteColors.value = cloneColors(colors);
  };
  const setPendingPaletteState = (paletteName, colors = currentColors.value) => {
    pendingPaletteName.value = String(paletteName || selectedPalette.value || '');
    pendingPaletteColors.value = cloneColors(colors);
  };
  const clearPendingPaletteDraft = () => {
    pendingPaletteName.value = '';
    pendingPaletteColors.value = {};
  };
  const applyPaletteDraftToPreview = () => {
    setAppliedPaletteState(selectedPalette.value, currentColors.value);
    clearPendingPaletteDraft();
  };
  const syncPaletteDraftState = () => {
    if (paletteInstantPreviewEnabled.value) {
      applyPaletteDraftToPreview();
      return;
    }

    if (String(pendingPaletteName.value || '').trim() !== '') {
      setPendingPaletteState(selectedPalette.value, currentColors.value);
      return;
    }

    setAppliedPaletteState(selectedPalette.value, currentColors.value);
  };

  const updatePalette = () => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const selectedName = String(selectedPalette.value || '').trim() || 'default';

    if (!paletteInstantPreviewEnabled.value && selectedName === appliedPaletteName.value) {
      currentColors.value = cloneColors(appliedPaletteColors.value);
      clearPendingPaletteDraft();
      return;
    }

    currentColors.value = getPaletteBaseColors(selectedName);
    if (paletteInstantPreviewEnabled.value) {
      applyPaletteDraftToPreview();
      return;
    }

    setPendingPaletteState(selectedName, currentColors.value);
  };

  const resetColors = () => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const selectedName = String(selectedPalette.value || '').trim() || 'default';
    currentColors.value = getPaletteBaseColors(selectedName);
    if (paletteInstantPreviewEnabled.value) {
      applyPaletteDraftToPreview();
      return;
    }

    if (String(pendingPaletteName.value || '').trim() !== '') {
      setPendingPaletteState(selectedName, currentColors.value);
      return;
    }

    setAppliedPaletteState(selectedName, currentColors.value);
  };

  return {
    updatePalette,
    resetColors,
    applyPaletteDraftToPreview,
    clearPendingPaletteDraft,
    setAppliedPaletteState,
    setPendingPaletteState,
    syncPaletteDraftState
  };
};
