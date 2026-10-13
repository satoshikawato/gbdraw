// @ts-check
import { drawnLegendRowStroke, legendEntryKey } from '../../services/legend-svg.js';
/** @import { PythonLegendKey } from '../../services/legend-svg.js' */
/** @typedef {{ caption?: unknown, originalCaption?: unknown }} StrokedLegendEntry A listed Legend row. */

/**
 * @typedef {object} LegendStrokeActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 */

// The Legend row stroke edits write the editor intent only. The composition
// root shows them on the Result through the executor as one History step
// (`editEditorIntent` in app/app-setup.js), which records Python's strokes,
// so a reset returns each feature part and swatch to the stroke Python drew.
// A row's stroke and its open Stroke options are stored under the row's key
// (`legendEntryKey`), so two rows of one name keep their own.
/** @param {LegendStrokeActionsOptions} options */
export const createLegendStrokeActions = ({ state }) => {
  const {
    svgContainer,
    legendStrokeOptionsOpen
  } = state;

  // OV-157: a row's Stroke options button shows or hides its stroke controls.
  // That is view state, kept out of the Legend entries, so the click records
  // no History step and the Session does not save it.
  /** @param {StrokedLegendEntry} entry */
  const isLegendStrokeOptionsOpen = (entry) => Boolean(legendStrokeOptionsOpen?.has(legendEntryKey(entry)));
  /** @param {StrokedLegendEntry} entry */
  const toggleLegendStrokeOptions = (entry) => {
    const key = legendEntryKey(entry);
    if (!legendStrokeOptionsOpen?.delete(key)) legendStrokeOptionsOpen?.add(key);
  };
  const closeLegendStrokeOptions = () => legendStrokeOptionsOpen?.clear();

  /** @param {PythonLegendKey} key */
  const drawnSwatchStroke = (key) => drawnLegendRowStroke(svgContainer.value?.querySelector?.('svg'), key);

  /** @param {number} idx */
  const getLegendEntryStrokeColor = (idx) => {
    const drawing = state.activeDrawing();
    const entry = drawing.legendEntries.value[idx];
    if (!entry) return '';
    const override = drawing.legendStrokeOverrides[legendEntryKey(entry)];
    if (override && override.strokeColor !== undefined) return override.strokeColor;
    return '';
  };

  /** @param {number} idx */
  const getLegendEntryStrokeWidth = (idx) => {
    const drawing = state.activeDrawing();
    const entry = drawing.legendEntries.value[idx];
    if (!entry) return '';
    const override = drawing.legendStrokeOverrides[legendEntryKey(entry)];
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
    const key = legendEntryKey(entry);
    if (String(drawing.legendStrokeOverrides[key]?.strokeColor || '').trim() === normalized) {
      return false;
    }
    // A key is free text: only an own entry is an override (OV-532).
    const existing = Object.hasOwn(drawing.legendStrokeOverrides, key)
      ? drawing.legendStrokeOverrides[key]
      : undefined;
    drawing.legendStrokeOverrides[key] = {
      ...(existing || {
        ...drawnSwatchStroke(key),
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
    const key = legendEntryKey(entry);
    if (Number(drawing.legendStrokeOverrides[key]?.strokeWidth) === widthVal) return false;
    // A key is free text: only an own entry is an override (OV-532).
    const existing = Object.hasOwn(drawing.legendStrokeOverrides, key)
      ? drawing.legendStrokeOverrides[key]
      : undefined;
    drawing.legendStrokeOverrides[key] = {
      ...(existing || {
        ...drawnSwatchStroke(key),
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
    const key = legendEntryKey(entry);
    const override = drawing.legendStrokeOverrides[key];
    if (!override || !Object.prototype.hasOwnProperty.call(override, 'strokeColor')) return false;
    delete override.strokeColor;
    if (override.strokeWidth === undefined || override.strokeWidth === '') {
      delete drawing.legendStrokeOverrides[key];
    }
    return true;
  };

  /** @param {number} idx */
  const resetLegendEntryStroke = (idx) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = drawing.legendEntries.value[idx];
    if (!entry || !Object.prototype.hasOwnProperty.call(drawing.legendStrokeOverrides, legendEntryKey(entry))) return false;
    delete drawing.legendStrokeOverrides[legendEntryKey(entry)];
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
