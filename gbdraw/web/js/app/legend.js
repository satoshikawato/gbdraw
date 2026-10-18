// @ts-check
/** @import { DrawingState } from '../state.js' */
/** @import { LegendRow } from '../services/drawing/paint.js' */
/** @import { SpecificRuleBuild } from '../services/specific-color-rules.js' */
import { createLegendDragActions } from './legend/drag-actions.js';
import { createLegendEntryActions } from './legend/entry-actions.js';
import { createLegendLayoutActions } from './legend/layout-actions.js';
import { createLegendSortActions } from './legend/sort-actions.js';
import { createLegendStrokeActions } from './legend/stroke-actions.js';
import {
  getAllFeatureLegendGroups,
  getVisibleFeatureLegendGroup,
  isCurrentLegendHorizontal,
  legendRowKeys
} from '../services/legend-svg.js';
import { legendEntryKey } from '../services/drawing/paint.js';
import { legendRowRules } from '../services/specific-color-rules.js';
import { displayedLegendRowContext } from '../services/feature-visibility.js';
import { registeredLegendFaces } from '../services/font-store.js';

/**
 * @typedef {object} LegendManagerOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(build: SpecificRuleBuild, label?: string, options?: { drawing?: DrawingState, afterCommit?: () => void }) => Promise<any>} commitLegendRowRules
 *   The rule owner's write of the rule table (`commitSpecificRules`).
 *   The rule owner's commit of the specific-color rules that a Legend row edit changes (R13).
 * @property {((label?: string, options?: { source?: string, owner?: unknown }) => Promise<any>) | null} [beginHistoryTransaction]
 *   History's begin of one step (R11); resolves to the transaction, or null when History is busy.
 * @property {((transaction: any, options?: Record<string, any>) => Promise<any>) | null} [commitHistoryTransaction]
 *   History's commit of the step that `beginHistoryTransaction` opened.
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 * @property {(() => string | undefined) | null} [readActiveResultIdentity]
 *   The preview owner's runtime identity of the mounted Result.
 * @property {import('./candidate-render.js').LegendRowColorOf} readShownLegendColor
 *   The root's color of a listed row as the Legend panel shows it (OV-282).
 * @property {import('./legend/entry-actions.js').LegendEntryActionsOptions['legendRowNames']} [legendRowNames]
 *   The root's reactive holder of the names the displayed Result shows for its rows; the entry owner writes it.
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
  readShownLegendColor,
  legendRowNames,
  showLegendStructure
}) => {
  const layoutActions = createLegendLayoutActions({
    readRegisteredFontFaces: () => registeredLegendFaces(state.loadedFonts.value)
  });
  const entryActions = createLegendEntryActions({ state, readActiveResultIdentity, readShownLegendColor, legendRowNames });
  const sortActions = createLegendSortActions({ state, showLegendStructure });
  const strokeActions = createLegendStrokeActions({ state });
  // The keys a listed row stands for on the displayed Result (D-35: its own
  // and those of the rows merged into it), each with the rules it draws, by its
  // key among Python's rows (OV-294).
  /** @param {DrawingState} drawing @param {number} index */
  const rowRulesAt = (drawing, index) => {
    /** @type {LegendRow | undefined} */
    const entry = drawing.legendEntries.value[index];
    if (!entry) return [];
    const context = displayedLegendRowContext(state, drawing);
    return legendRowKeys(state.svgContainer.value?.querySelector?.('svg'), legendEntryKey(entry))
      .map((key) => ({ key, rules: legendRowRules(key, context) }));
  };
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
    // legend-only edit, which the root shows on the Result. A row that stands
    // for rows of both kinds (D-35) edits the rules, and the Legend color of
    // the others in the rule commit's step.
    /** @param {number} index */
    legendRowHasRules: (index) => rowRulesAt(state.activeDrawing(), index).some((row) => row.rules.length > 0),
    /** @param {number} index @param {string} color */
    updateLegendEntryColor: (index, color) => {
      const drawing = state.activeDrawing();
      const rows = rowRulesAt(drawing, index);
      const rowRules = rows.flatMap((row) => row.rules);
      if (rowRules.length) {
        const legendOnly = rows.filter((row) => row.rules.length === 0).map((row) => row.key);
        return commitLegendRowRules((current, revise) => current.map(rule => (rowRules.includes(rule) ? revise(rule, { color }) : rule)), 'Change legend color',
          { drawing, afterCommit: () => entryActions.writeLegendRowColor(drawing, legendOnly, color) });
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
