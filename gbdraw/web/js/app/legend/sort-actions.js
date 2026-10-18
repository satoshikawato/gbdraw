// @ts-check
/** @import { DrawingState } from '../../state.js' */
/** @import { LegendRow } from '../../services/drawing/paint.js' */
import { legendEntryKey, replaceLegendRows } from '../../services/drawing/paint.js';
// Sort and Move compute the requested row order, write it into the
// drawing's Legend intent, and show it through the root's port once (R1, R3,
// R13): the executor orders the displayed Result's rows.
/**
 * @typedef {object} LegendSortActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {() => unknown} showLegendStructure The root's show of the Legend structure intent.
 */

/** @param {LegendSortActionsOptions} options */
export const createLegendSortActions = ({ state, showLegendStructure }) => {
  const { originalLegendOrder } = state;

  // Each row is moved as itself, never by its name: two rows may show one name.
  /** @param {DrawingState} drawing @param {readonly LegendRow[]} ordered The listed rows in their new order. */
  const applyLegendEntryOrder = (drawing, ordered) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entries = drawing.legendEntries.value || [];
    if (ordered.every((entry, index) => entry === entries[index])) return;
    replaceLegendRows(drawing.slices.paint, ordered);
    showLegendStructure();
  };

  /** @param {DrawingState} drawing @returns {LegendRow[]} */
  const listedRows = (drawing) => (drawing.legendEntries.value || []).filter((entry) => entry.caption);

  // Moves the row at `from` to `to` (the Legend panel's drag handle and the
  // handle's arrow keys).
  /** @param {number} from @param {number} to */
  const moveLegendEntry = (from, to) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const current = listedRows(drawing);
    if (from === to || from < 0 || to < 0 || from >= current.length || to >= current.length) return;
    const next = [...current];
    next.splice(to, 0, ...next.splice(from, 1));
    applyLegendEntryOrder(drawing, next);
  };

  /** @param {LegendRow[]} rows @param {(left: LegendRow, right: LegendRow) => number} compare */
  const stableSorted = (rows, compare) => rows.map((entry, index) => ({ entry, index }))
    .sort((left, right) => compare(left.entry, right.entry) || left.index - right.index)
    .map(({ entry }) => entry);

  const sortLegendEntries = (direction = 'asc') => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const current = listedRows(drawing);
    if (current.length < 2) return;
    const sign = direction === 'asc' ? 1 : -1;
    applyLegendEntryOrder(drawing, stableSorted(current, (left, right) => (
      sign * String(left.caption).localeCompare(String(right.caption), undefined, { sensitivity: 'base' })
    )));
  };

  const sortLegendEntriesByDefault = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const current = listedRows(drawing);
    if (current.length < 2) return;
    if (originalLegendOrder.value.length === 0) return;
    // A renamed generated row keeps its key's place in the generated order.
    /** @param {LegendRow} entry */
    const rank = (entry) => originalLegendOrder.value.indexOf(legendEntryKey(entry));
    applyLegendEntryOrder(drawing, stableSorted(current, (left, right) => {
      const [leftRank, rightRank] = [rank(left), rank(right)];
      if (leftRank !== -1 && rightRank !== -1) return leftRank - rightRank;
      if (leftRank !== -1) return -1;
      if (rightRank !== -1) return 1;
      return String(left.caption).localeCompare(String(right.caption), undefined, { sensitivity: 'base' });
    }));
  };

  return {
    moveLegendEntry,
    sortLegendEntries,
    sortLegendEntriesByDefault
  };
};
