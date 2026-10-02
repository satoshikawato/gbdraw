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
import { legendRowRules } from './specific-color-rules.js';

export const createLegendManager = ({
  state,
  rulePreparation,
  commitSpecificRules,
  history = null,
  previewRuntime = null,
  getCommittedRequest = () => null
}) => {
  const layoutActions = createLegendLayoutActions({ state });
  const entryActions = createLegendEntryActions({
    state,
    layoutActions,
    previewRuntime,
    getCommittedRequest
  });
  const sortActions = createLegendSortActions({ state, extractLegendEntries: entryActions.extractLegendEntries });
  const strokeActions = createLegendStrokeActions({ state, previewRuntime });
  const rowRulesAt = (index) => legendRowRules(state.legendEntries.value[index]?.caption, {
    rules: state.manualSpecificRules,
    legendEntries: state.legendEntries.value,
    originalLegendOrder: state.originalLegendOrder?.value || []
  });
  const dragActions = createLegendDragActions({
    state,
    extractLegendEntries: entryActions.extractLegendEntries,
    history
  });

  return {
    ...entryActions,
    // A row a rule draws, including its N-06 "<caption> [<hex>]" row, edits
    // that rule; any other row is a legend-only edit.
    updateLegendEntryColor: (index, color) => {
      const rowRules = rowRulesAt(index);
      if (rowRules.length) {
        return commitSpecificRules(state.manualSpecificRules.map(rule => rowRules.includes(rule) ? { ...rule, color } : { ...rule }), 'Change legend color');
      }
      return entryActions.updateLegendEntryColor(index, color);
    },
    updateLegendEntryCaption: (index, caption) => {
      const rowRules = rowRulesAt(index);
      if (rowRules.length) {
        return commitSpecificRules(state.manualSpecificRules.map(rule => rowRules.includes(rule) ? { ...rule, cap: caption } : { ...rule }), 'Rename legend item');
      }
      return entryActions.updateLegendEntryCaption(index, caption);
    },
    ...layoutActions,
    ...sortActions,
    ...strokeActions,
    captureOriginalStrokeValues: (...args) => rulePreparation.run(state.manualSpecificRules, () => strokeActions.captureOriginalStrokeValues(...args)),
    ...dragActions,
    getAllFeatureLegendGroups,
    getVisibleFeatureLegendGroup,
    isCurrentLegendHorizontal
  };
};
