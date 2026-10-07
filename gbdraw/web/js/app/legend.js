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
  isCurrentLegendHorizontal
} from './legend/utils.js';
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
 * @property {() => ({ diagramOptions?: Record<string, any> } | null)} [getCommittedRequest]
 *   The committed canonical request (Python owns the option fields, R7).
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
  getCommittedRequest = () => null
}) => {
  const layoutActions = createLegendLayoutActions();
  const entryActions = createLegendEntryActions({
    state,
    updatePairwiseLegendPositions: layoutActions.updatePairwiseLegendPositions,
    reflowDualLegendLayout: layoutActions.reflowDualLegendLayout,
    compactLegendEntries: layoutActions.compactLegendEntries,
    commitActiveResultEdit,
    readActiveResultIdentity,
    getCommittedRequest
  });
  const sortActions = createLegendSortActions({
    state,
    extractLegendEntries: entryActions.extractLegendEntries,
    orderMountedLegend: entryActions.orderMountedLegend,
    commitActiveResultEdit
  });
  const strokeActions = createLegendStrokeActions({ state, commitActiveResultEdit });
  /** @param {DrawingState} drawing */
  const rowRulesAt = (drawing, index) => legendRowRules(drawing.legendEntries.value[index]?.caption, {
    rules: drawing.manualSpecificRules,
    legendEntries: drawing.legendEntries.value,
    originalLegendOrder: state.originalLegendOrder?.value || []
  });
  const dragActions = createLegendDragActions({
    state,
    extractLegendEntries: entryActions.extractLegendEntries,
    beginHistoryTransaction,
    commitHistoryTransaction,
    commitActiveResultEdit
  });

  return {
    ...entryActions,
    // A row a rule draws, including its N-06 "<caption> [<hex>]" row, edits
    // that rule through the rule owner's port (R13); any other row is a
    // legend-only edit.
    updateLegendEntryColor: (index, color) => {
      const drawing = state.activeDrawing();
      const rowRules = rowRulesAt(drawing, index);
      if (rowRules.length) {
        return commitLegendRowRules(drawing.manualSpecificRules.map(rule => rowRules.includes(rule) ? { ...rule, color } : { ...rule }), 'Change legend color');
      }
      return entryActions.updateLegendEntryColor(index, color);
    },
    updateLegendEntryCaption: (index, caption) => {
      const drawing = state.activeDrawing();
      const rowRules = rowRulesAt(drawing, index);
      if (rowRules.length) {
        return commitLegendRowRules(drawing.manualSpecificRules.map(rule => rowRules.includes(rule) ? { ...rule, cap: caption } : { ...rule }), 'Rename legend item');
      }
      return entryActions.updateLegendEntryCaption(index, caption);
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
