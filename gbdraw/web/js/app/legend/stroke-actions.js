// @ts-check
import { drawnLegendRowStroke, legendRowKeys } from '../../services/legend-svg.js';
import { clearStrokes, legendEntryKey, setLegendStroke } from '../../services/drawing/paint.js';
/** @import { PythonLegendKey } from '../../services/legend-svg.js' */
/** @import { StrokeOverride } from '../../services/style-overrides.js' */
/** @import { WebState } from '../../state.js' */
/** @typedef {{ caption?: unknown, originalCaption?: unknown }} StrokedLegendEntry A listed Legend row. */

/**
 * @typedef {object} LegendStrokeActionsOptions
 * @property {Pick<WebState, 'activeDrawing' | 'sessionOperationAvailability' | 'svgContainer' | 'legendStrokeOptionsOpen'>} state
 */

// The Legend row stroke edits write the editor intent only. The editor intent
// shows them on the Result through the executor as one History step
// (`editEditorIntent` in app/editor-intent.js), which records Python's strokes,
// so a reset returns each feature part and swatch to the stroke Python drew.
// A row's stroke and its open Stroke options are stored under the row's key
// (`legendEntryKey`), so two rows of one name keep their own; an edit of a row
// that the Legend draws for merged rows writes each of their keys (D-35).
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

  // The listed row at `idx`, and the keys of the rows it stands for (D-35:
  // a merged row's edit reaches every row it merges, `legendRowKeys`).
  /** @param {number} idx @returns {{ entry: StrokedLegendEntry, keys: PythonLegendKey[] } | null} */
  const rowAt = (idx) => {
    const entry = state.activeDrawing().legendEntries.value[idx];
    return entry ? { entry, keys: legendRowKeys(svgContainer.value?.querySelector?.('svg'), legendEntryKey(entry)) } : null;
  };

  /** @param {number} idx @param {'strokeColor' | 'strokeWidth'} field */
  const strokeFieldAt = (idx, field) => {
    const entry = state.activeDrawing().legendEntries.value[idx];
    const override = entry ? state.activeDrawing().legendStrokeOverrides[legendEntryKey(entry)] : undefined;
    return override && override[field] !== undefined ? override[field] : '';
  };
  /** @param {number} idx */
  const getLegendEntryStrokeColor = (idx) => strokeFieldAt(idx, 'strokeColor');
  /** @param {number} idx */
  const getLegendEntryStrokeWidth = (idx) => strokeFieldAt(idx, 'strokeWidth');

  // Writes one field of the stroke of every row `idx` stands for; a row
  // without its own override starts from the stroke Python drew on its swatch
  // and holds only the field written. No edit when every row has the value.
  /**
   * @template {'strokeColor' | 'strokeWidth'} F
   * @param {number} idx
   * @param {F} field
   * @param {Required<StrokeOverride>[F]} value
   * @param {(current: unknown) => boolean} same Whether a row's current value of `field` is `value`.
   */
  const writeStrokeField = (idx, field, value, same) => {
    const drawing = state.activeDrawing();
    const row = rowAt(idx);
    if (!row) return false;
    // A key is free text: only an own entry is an override (OV-532).
    /** @param {PythonLegendKey} key */
    const ownOverride = (key) => (Object.hasOwn(drawing.legendStrokeOverrides, key) ? drawing.legendStrokeOverrides[key] : undefined);
    if (row.keys.every((key) => same(ownOverride(key)?.[field]))) return false;
    row.keys.forEach((key) => setLegendStroke(drawing.slices.paint, key, { ...(ownOverride(key) || drawnSwatchStroke(key)), [field]: value }));
    return true;
  };

  /** @param {number} idx @param {string} color */
  const updateLegendEntryStrokeColor = (idx, color) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const normalized = String(color || '').trim();
    return writeStrokeField(idx, 'strokeColor', normalized, (current) => String(current || '').trim() === normalized);
  };

  /** @param {number} idx @param {string | number} width */
  const updateLegendEntryStrokeWidth = (idx, width) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const widthVal = parseFloat(String(width));
    if (isNaN(widthVal)) return false;
    return writeStrokeField(idx, 'strokeWidth', widthVal, (current) => Number(current) === widthVal);
  };

  /** @param {number} idx @param {string | null} value null removes the row's stroke color. */
  const setLegendEntryStrokeColorValue = (idx, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const row = rowAt(idx);
    if (!row) return false;
    if (value !== null) return updateLegendEntryStrokeColor(idx, String(value || '').trim());
    const stroked = row.keys.filter((key) => Object.hasOwn(drawing.legendStrokeOverrides[key] || {}, 'strokeColor'));
    stroked.forEach((key) => {
      const { strokeColor: _removed, ...override } = drawing.legendStrokeOverrides[key];
      setLegendStroke(drawing.slices.paint, key, override.strokeWidth === undefined ? null : override);
    });
    return stroked.length > 0;
  };

  /** @param {number} idx */
  const resetLegendEntryStroke = (idx) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const stroked = (rowAt(idx)?.keys || []).filter((key) => Object.hasOwn(drawing.legendStrokeOverrides, key));
    stroked.forEach((key) => setLegendStroke(drawing.slices.paint, key, null));
    return stroked.length > 0;
  };

  const resetAllStrokes = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    return clearStrokes(drawing.slices.paint);
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
