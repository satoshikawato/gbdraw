// @ts-check
import { normalizePaletteColors } from '../utils/color-utils.js';

/**
 * The palette refs the manager reads and writes (state.js owns them).
 * @typedef {object} ResultsManagerState
 * @property {{ value: Record<string, Record<string, string>> | null }} paletteDefinitions
 * @property {{ value: boolean }} paletteInstantPreviewEnabled
 * @property {{ value: string }} appliedPaletteName
 * @property {{ value: Record<string, string> }} appliedPaletteColors
 * @property {() => ResultsManagerDrawing} activeDrawing The shown mode's drawing.
 * @property {() => any} [sessionOperationAvailability] The busy outcome of a Session operation, if any.
 */

/**
 * The palette members of a drawing (`DrawingState` of state.js).
 * @typedef {object} ResultsManagerDrawing
 * @property {{ value: string }} selectedPalette
 * @property {{ value: Record<string, string> }} currentColors
 * @property {{ value: string }} pendingPaletteName
 * @property {{ value: Record<string, string> }} pendingPaletteColors
 */

/** @param {{ state: ResultsManagerState }} options */
export const createResultsManager = ({ state }) => {
  const {
    paletteDefinitions,
    paletteInstantPreviewEnabled,
    appliedPaletteName,
    appliedPaletteColors
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
  /** @param {ResultsManagerDrawing} drawing */
  const setAppliedPaletteState = (drawing, paletteName, colors = drawing.currentColors.value) => {
    appliedPaletteName.value = String(paletteName || drawing.selectedPalette.value || 'default');
    appliedPaletteColors.value = cloneColors(colors);
  };
  /** @param {ResultsManagerDrawing} drawing */
  const setPendingPaletteState = (drawing, paletteName, colors = drawing.currentColors.value) => {
    drawing.pendingPaletteName.value = String(paletteName || drawing.selectedPalette.value || '');
    drawing.pendingPaletteColors.value = cloneColors(colors);
  };
  /** @param {ResultsManagerDrawing} drawing */
  const clearPendingPaletteDraft = (drawing) => {
    drawing.pendingPaletteName.value = '';
    drawing.pendingPaletteColors.value = {};
  };
  const applyPaletteDraftToPreview = () => {
    const drawing = state.activeDrawing();
    setAppliedPaletteState(drawing, drawing.selectedPalette.value, drawing.currentColors.value);
    clearPendingPaletteDraft(drawing);
  };
  const syncPaletteDraftState = () => {
    const drawing = state.activeDrawing();
    if (paletteInstantPreviewEnabled.value) {
      applyPaletteDraftToPreview();
      return;
    }

    if (String(drawing.pendingPaletteName.value || '').trim() !== '') {
      setPendingPaletteState(drawing, drawing.selectedPalette.value, drawing.currentColors.value);
      return;
    }

    setAppliedPaletteState(drawing, drawing.selectedPalette.value, drawing.currentColors.value);
  };

  const updatePalette = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const selectedName = String(drawing.selectedPalette.value || '').trim() || 'default';

    if (!paletteInstantPreviewEnabled.value && selectedName === appliedPaletteName.value) {
      drawing.currentColors.value = cloneColors(appliedPaletteColors.value);
      clearPendingPaletteDraft(drawing);
      return;
    }

    drawing.currentColors.value = getPaletteBaseColors(selectedName);
    if (paletteInstantPreviewEnabled.value) {
      applyPaletteDraftToPreview();
      return;
    }

    setPendingPaletteState(drawing, selectedName, drawing.currentColors.value);
  };

  const resetColors = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const selectedName = String(drawing.selectedPalette.value || '').trim() || 'default';
    drawing.currentColors.value = getPaletteBaseColors(selectedName);
    if (paletteInstantPreviewEnabled.value) {
      applyPaletteDraftToPreview();
      return;
    }

    if (String(drawing.pendingPaletteName.value || '').trim() !== '') {
      setPendingPaletteState(drawing, selectedName, drawing.currentColors.value);
      return;
    }

    setAppliedPaletteState(drawing, selectedName, drawing.currentColors.value);
  };

  return {
    updatePalette,
    resetColors,
    applyPaletteDraftToPreview,
    syncPaletteDraftState
  };
};
