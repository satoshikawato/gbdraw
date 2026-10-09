// @ts-check
/** @import { PaletteColorsDialog } from '../state.js' */
import { buildPaletteColorOverrideRows, normalizePaletteColors } from '../utils/color-utils.js';

// The palette dialog body names up to three keys (Q6 A, Owner 2026-10-09).
/** @param {string[]} keys */
const describeKeys = (keys) => {
  const rest = keys.length - 3;
  const parts = rest > 0 ? [...keys.slice(0, 3), `${rest} more`] : keys;
  if (parts.length <= 2) return parts.join(' and ');
  return `${parts.slice(0, -1).join(', ')}, and ${parts[parts.length - 1]}`;
};

/**
 * The palette refs the manager reads and writes (state.js owns them).
 * @typedef {object} ResultsManagerState
 * @property {{ value: Record<string, Record<string, string>> | null }} paletteDefinitions
 * @property {{ value: boolean }} paletteInstantPreviewEnabled
 * @property {{ value: string }} appliedPaletteName
 * @property {{ value: Record<string, string> }} appliedPaletteColors
 * @property {() => ResultsManagerDrawing} activeDrawing The shown mode's drawing.
 * @property {() => any} [sessionOperationAvailability] The busy outcome of a Session operation, if any.
 * @property {PaletteColorsDialog} paletteColorsDialog
 */

/**
 * The palette members of a drawing (`DrawingState` of state.js).
 * @typedef {object} ResultsManagerDrawing
 * @property {{ value: string }} selectedPalette
 * @property {{ value: Record<string, string> }} currentColors
 * @property {{ value: string }} pendingPaletteName
 * @property {{ value: Record<string, string> }} pendingPaletteColors
 */

/**
 * @param {{
 *   state: ResultsManagerState,
 *   closeAfterDialogChoice?: (close: () => unknown) => void
 * }} options `closeAfterDialogChoice` runs the palette dialog's close at once,
 *   or once the History step of its choice in flight ends (D-12, OIC-028).
 */
export const createResultsManager = ({ state, closeAfterDialogChoice = (close) => { close(); } }) => {
  const {
    paletteDefinitions,
    paletteInstantPreviewEnabled,
    appliedPaletteName,
    appliedPaletteColors,
    paletteColorsDialog
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
  const paletteNameOf = (drawing) => String(drawing.selectedPalette.value || '').trim() || 'default';
  // D-15: a user default color is a value that differs from the selected
  // palette's color for its key, by the comparator of the `-d` table; Auto
  // (empty) is none.
  /**
   * @param {ResultsManagerDrawing} drawing
   * @param {Record<string, string | null>} [colors]
   * @returns {[string, string][]}
   */
  const userDefaultColorRows = (drawing, colors = drawing.currentColors.value) => buildPaletteColorOverrideRows({
    colors, paletteColors: getPaletteBaseColors(paletteNameOf(drawing))
  });
  /**
   * @param {ResultsManagerDrawing} drawing
   * @param {string} key
   * @returns {string | null}
   */
  const readUserDefaultColor = (drawing, key) => (
    userDefaultColorRows(drawing, { [key]: drawing.currentColors.value?.[key] })[0]?.[1] ?? null
  );
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

  /**
   * @param {PaletteColorsDialog['kind']} kind
   * @param {ResultsManagerDrawing} drawing
   * @param {[string, string][]} rows
   * @param {string} toPalette
   */
  const openPaletteColorsDialog = (kind, drawing, rows, toPalette) => {
    // Switching back to the applied palette drops the queued one (Q5 A, D-24).
    const queued = String(drawing.pendingPaletteName.value || '').trim();
    Object.assign(paletteColorsDialog, {
      show: true, kind, fromPalette: paletteNameOf(drawing), toPalette,
      droppedPalette: kind === 'switch' && queued && toPalette === appliedPaletteName.value ? queued : '',
      count: rows.length, keysText: describeKeys(rows.map(([key]) => key))
    });
  };
  const closePaletteColorsDialog = () => closeAfterDialogChoice(() => {
    Object.assign(paletteColorsDialog, {
      show: false, kind: 'switch', fromPalette: '', toPalette: '', droppedPalette: '', count: 0, keysText: ''
    });
  });

  // D-15: a palette switch keeps the user default colors (`keep`) or takes the
  // palette's colors (`palette`); without a choice it asks while any exist.
  /**
   * @param {string} name
   * @param {'' | 'keep' | 'palette'} [choice]
   */
  const selectPalette = (name, choice = '') => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const paletteName = String(name || '').trim() || 'default';
    const rows = userDefaultColorRows(drawing);
    if (!choice && rows.length > 0) {
      openPaletteColorsDialog('switch', drawing, rows, paletteName);
      return undefined;
    }
    const colors = getPaletteBaseColors(paletteName);
    if (choice === 'keep') rows.forEach(([key, color]) => { colors[key] = color; });
    drawing.selectedPalette.value = paletteName;
    drawing.currentColors.value = colors;
    // The applied palette queues nothing: switching back to it while another
    // is queued applies its colors live (Q5 A, Owner 2026-10-09).
    if (paletteInstantPreviewEnabled.value || paletteName === appliedPaletteName.value) {
      applyPaletteDraftToPreview();
      return undefined;
    }

    setPendingPaletteState(drawing, paletteName, colors);
    return undefined;
  };

  // The Palette select shows the drawing's palette until a switch applies, so
  // a dialog's Cancel leaves it unchanged.
  /** @param {{ target?: { value: string } | null } | null} [event] */
  const requestPaletteChange = (event) => {
    const target = event?.target;
    const name = String(target?.value || '');
    if (target) target.value = state.activeDrawing().selectedPalette.value;
    return selectPalette(name);
  };

  /** @param {ResultsManagerDrawing} drawing */
  const resetColors = (drawing) => {
    drawing.currentColors.value = getPaletteBaseColors(paletteNameOf(drawing));
    syncPaletteDraftState();
  };

  // D-15: Default colors Reset asks before it discards user default colors.
  const requestResetColors = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const rows = userDefaultColorRows(drawing);
    if (rows.length > 0) {
      openPaletteColorsDialog('reset', drawing, rows, paletteNameOf(drawing));
      return undefined;
    }
    resetColors(drawing);
    return undefined;
  };

  // A choice of the palette dialog; app-setup.js makes it one History step
  // and answers Cancel with `closePaletteColorsDialog` (PD-OI-088).
  /** @param {string} choice */
  const handlePaletteColorsChoice = (choice) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (paletteColorsDialog.kind === 'reset') {
      if (choice === 'reset') resetColors(drawing);
    } else if (choice === 'keep' || choice === 'palette') {
      selectPalette(paletteColorsDialog.toPalette, choice);
    }
    closePaletteColorsDialog();
    return undefined;
  };

  // D-15: the popup's "Apply to all" on a palette row sets the type's default
  // color, as an edit in the Default colors list does. While a palette is
  // queued, the shown Result takes the color now too: a user color wins over
  // any palette, so the next Generate draws it (Q1 B, Owner 2026-10-09).
  /**
   * @param {ResultsManagerDrawing} drawing
   * @param {string} key
   * @param {string} color
   */
  const setDefaultColor = (drawing, key, color) => {
    drawing.currentColors.value = { ...drawing.currentColors.value, [key]: color };
    if (drawing !== state.activeDrawing()) return;
    syncPaletteDraftState();
    if (String(drawing.pendingPaletteName.value || '').trim() !== '') {
      appliedPaletteColors.value = { ...appliedPaletteColors.value, [key]: color };
    }
  };

  return {
    requestPaletteChange,
    selectPalette,
    requestResetColors,
    handlePaletteColorsChoice,
    cancelPaletteColorsDialog: closePaletteColorsDialog,
    readUserDefaultColor,
    setDefaultColor,
    applyPaletteDraftToPreview,
    syncPaletteDraftState
  };
};
