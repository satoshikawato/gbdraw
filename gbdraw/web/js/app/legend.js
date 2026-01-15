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

export const createLegendManager = ({ state, getPyodide, debugLog }) => {
  const layoutActions = createLegendLayoutActions({ state });
  const entryActions = createLegendEntryActions({ state, getPyodide, layoutActions });
  const sortActions = createLegendSortActions({ state, extractLegendEntries: entryActions.extractLegendEntries });
  const strokeActions = createLegendStrokeActions({ state, debugLog });
  const dragActions = createLegendDragActions({ state, extractLegendEntries: entryActions.extractLegendEntries });

  return {
    ...entryActions,
    ...layoutActions,
    ...sortActions,
    ...strokeActions,
    ...dragActions,
    getAllFeatureLegendGroups,
    getVisibleFeatureLegendGroup,
    isCurrentLegendHorizontal
  };
};
