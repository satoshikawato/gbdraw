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

export const createLegendManager = ({
  state,
  rulePreparation,
  history = null,
  previewRuntime = null
}) => {
  const layoutActions = createLegendLayoutActions({ state });
  const entryActions = createLegendEntryActions({
    state,
    layoutActions,
    previewRuntime
  });
  const sortActions = createLegendSortActions({ state, extractLegendEntries: entryActions.extractLegendEntries });
  const strokeActions = createLegendStrokeActions({ state, previewRuntime });
  const dragActions = createLegendDragActions({
    state,
    extractLegendEntries: entryActions.extractLegendEntries,
    history
  });

  return {
    ...entryActions,
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
