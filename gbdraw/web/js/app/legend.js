// @ts-check
/** @import { DrawingState } from '../state.js' */
import { createLegendDragActions } from './legend/drag-actions.js';
import { createLegendEntryActions } from './legend/entry-actions.js';
import { createLegendLayoutActions } from './legend/layout-actions.js';
import { createLegendSortActions } from './legend/sort-actions.js';
import { createLegendStrokeActions } from './legend/stroke-actions.js';
import {
  getAllFeatureLegendGroups,
  getVisibleFeatureLegendGroup,
  isCurrentLegendHorizontal,
  pythonLegendRows
} from '../services/legend-svg.js';
import { legendRowRules } from '../services/specific-color-rules.js';

/**
 * @typedef {object} LegendManagerOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(rows: Array<Record<string, any>>, label?: string) => Promise<any>} commitLegendRowRules
 *   The rule owner's commit of the specific-color rules that a Legend row edit changes (R13).
 * @property {((label?: string, options?: { source?: string, owner?: unknown }) => Promise<any>) | null} [beginHistoryTransaction]
 *   History's begin of one step (R11); resolves to the transaction, or null when History is busy.
 * @property {((transaction: any, options?: Record<string, any>) => Promise<any>) | null} [commitHistoryTransaction]
 *   History's commit of the step that `beginHistoryTransaction` opened.
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 * @property {(() => string | undefined) | null} [readActiveResultIdentity]
 *   The preview owner's runtime identity of the mounted Result.
 * @property {() => unknown} showLegendStructure The root's show of the Legend structure intent (R1).
 */

/** @param {LegendManagerOptions} options */
export const createLegendManager = ({
  state,
  commitLegendRowRules,
  beginHistoryTransaction = null,
  commitHistoryTransaction = null,
  // R13: the preview owner's ports; the Legend owners never hold it.
  commitActiveResultEdit = null,
  readActiveResultIdentity = null,
  showLegendStructure
}) => {
  const layoutActions = createLegendLayoutActions();
  const entryActions = createLegendEntryActions({ state, readActiveResultIdentity });
  const sortActions = createLegendSortActions({ state, showLegendStructure });
  const strokeActions = createLegendStrokeActions({ state });
  // The rules a row draws, by Python's rows of the displayed Result (OV-294).
  /** @param {DrawingState} drawing @param {number} index */
  const rowRulesAt = (drawing, index) => legendRowRules(drawing.legendEntries.value[index]?.caption, {
    rules: drawing.manualSpecificRules,
    pythonRows: pythonLegendRows(state.svgContainer?.value?.querySelector?.('svg')),
    features: state.extractedFeatures?.value || [],
    originalLegendOrder: state.originalLegendOrder?.value || []
  });
  const dragActions = createLegendDragActions({
    state,
    beginHistoryTransaction,
    commitHistoryTransaction,
    commitActiveResultEdit
  });

  return {
    ...entryActions,
    // A row a rule draws, including its N-06 "<caption> [<hex>]" row, edits
    // that rule through the rule owner's port (R13); any other row is a
    // legend-only edit, which the root shows on the Result.
    /** @param {number} index */
    legendRowHasRules: (index) => rowRulesAt(state.activeDrawing(), index).length > 0,
    updateLegendEntryColor: (index, color) => {
      const drawing = state.activeDrawing();
      const rowRules = rowRulesAt(drawing, index);
      if (rowRules.length) {
        return commitLegendRowRules(drawing.manualSpecificRules.map(rule => rowRules.includes(rule) ? { ...rule, color } : { ...rule }), 'Change legend color');
      }
      return entryActions.updateLegendEntryColor(index, color);
    },
    ...layoutActions,
    ...sortActions,
    ...strokeActions,
    ...dragActions,
    getAllFeatureLegendGroups,
    getVisibleFeatureLegendGroup,
    isCurrentLegendHorizontal
  };
};
