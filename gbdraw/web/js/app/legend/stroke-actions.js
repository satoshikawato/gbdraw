// @ts-check
import { drawnLegendRowStroke } from '../../services/legend-svg.js';

/**
 * @typedef {object} LegendStrokeActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 */

// The Legend row stroke edits write the editor intent only. The composition
// root shows them on the Result through the executor as one History step
// (`editEditorIntent` in app/app-setup.js), which records Python's strokes,
// so a reset returns each feature part and swatch to the stroke Python drew.
/** @param {LegendStrokeActionsOptions} options */
export const createLegendStrokeActions = ({ state }) => {
  const {
    svgContainer,
    legendStrokeOptionsOpen
  } = state;

  // OV-157: a row's Stroke options button shows or hides its stroke controls.
  // That is view state, kept out of the Legend entries, so the click records
  // no History step and the Session does not save it.
  /** @param {string} caption */
  const isLegendStrokeOptionsOpen = (caption) => Boolean(legendStrokeOptionsOpen?.has(caption));
  /** @param {string} caption */
  const toggleLegendStrokeOptions = (caption) => {
    if (!legendStrokeOptionsOpen?.delete(caption)) legendStrokeOptionsOpen?.add(caption);
  };
  const closeLegendStrokeOptions = () => legendStrokeOptionsOpen?.clear();

  /** @param {string} caption */
  const drawnSwatchStroke = (caption) => drawnLegendRowStroke(svgContainer.value?.querySelector?.('svg'), caption);

  const getLegendEntryStrokeColor = (idx) => {
    const drawing = state.activeDrawing();
    const entry = drawing.legendEntries.value[idx];
    if (!entry) return '';
    const override = drawing.legendStrokeOverrides[entry.caption];
    if (override && override.strokeColor !== undefined) return override.strokeColor;
    return '';
  };

  const getLegendEntryStrokeWidth = (idx) => {
    const drawing = state.activeDrawing();
    const entry = drawing.legendEntries.value[idx];
    if (!entry) return '';
    const override = drawing.legendStrokeOverrides[entry.caption];
    if (override && override.strokeWidth !== undefined) return override.strokeWidth;
    return '';
  };

  /** @param {number} idx @param {string} color */
  const updateLegendEntryStrokeColor = (idx, color) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = drawing.legendEntries.value[idx];
    if (!entry) return false;
    const normalized = String(color || '').trim();
    if (String(drawing.legendStrokeOverrides[entry.caption]?.strokeColor || '').trim() === normalized) {
      return false;
    }
    drawing.legendStrokeOverrides[entry.caption] = {
      ...(drawing.legendStrokeOverrides[entry.caption] || {
        ...drawnSwatchStroke(entry.caption),
        strokeWidth: getLegendEntryStrokeWidth(idx)
      }),
      strokeColor: normalized
    };
    return true;
  };

  /** @param {number} idx @param {string | number} width */
  const updateLegendEntryStrokeWidth = (idx, width) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = drawing.legendEntries.value[idx];
    if (!entry) return false;
    const widthVal = parseFloat(String(width));
    if (isNaN(widthVal)) return false;
    if (Number(drawing.legendStrokeOverrides[entry.caption]?.strokeWidth) === widthVal) return false;
    drawing.legendStrokeOverrides[entry.caption] = {
      ...(drawing.legendStrokeOverrides[entry.caption] || {
        ...drawnSwatchStroke(entry.caption),
        strokeColor: getLegendEntryStrokeColor(idx)
      }),
      strokeWidth: widthVal
    };
    return true;
  };

  /** @param {number} idx @param {string | null} value null removes the row's stroke color. */
  const setLegendEntryStrokeColorValue = (idx, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = drawing.legendEntries.value[idx];
    if (!entry) return false;
    if (value !== null) return updateLegendEntryStrokeColor(idx, String(value || '').trim());
    const override = drawing.legendStrokeOverrides[entry.caption];
    if (!override || !Object.prototype.hasOwnProperty.call(override, 'strokeColor')) return false;
    delete override.strokeColor;
    if (override.strokeWidth === undefined || override.strokeWidth === '') {
      delete drawing.legendStrokeOverrides[entry.caption];
    }
    return true;
  };

  /** @param {number} idx */
  const resetLegendEntryStroke = (idx) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = drawing.legendEntries.value[idx];
    if (!entry || !Object.prototype.hasOwnProperty.call(drawing.legendStrokeOverrides, entry.caption)) return false;
    delete drawing.legendStrokeOverrides[entry.caption];
    return true;
  };

  const resetAllStrokes = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const overridesRemoved =
      Object.keys(drawing.legendStrokeOverrides).length > 0 ||
      Object.keys(drawing.featureStrokeOverrides).length > 0;
    Object.keys(drawing.legendStrokeOverrides).forEach((key) => delete drawing.legendStrokeOverrides[key]);
    Object.keys(drawing.featureStrokeOverrides).forEach((key) => delete drawing.featureStrokeOverrides[key]);
    return overridesRemoved;
  };

  return {
    closeLegendStrokeOptions,
    getLegendEntryStrokeColor,
    getLegendEntryStrokeWidth,
    isLegendStrokeOptionsOpen,
    resetAllStrokes,
    resetLegendEntryStroke,
    setLegendEntryStrokeColorValue,
    toggleLegendStrokeOptions,
    updateLegendEntryStrokeColor,
    updateLegendEntryStrokeWidth
  };
};
