import { createLegendDragActions } from './legend/drag-actions.js';
import { createLegendEntryActions } from './legend/entry-actions.js';
import { createLegendLayoutActions } from './legend/layout-actions.js';
import { createLegendSortActions } from './legend/sort-actions.js';
import {
  getAllFeatureLegendGroups,
  getVisibleFeatureLegendGroup,
  isCurrentLegendHorizontal
} from './legend/utils.js';

export const createLegendManager = ({
  state,
  getPyodide,
  ensurePyodide = null,
  history = null,
  previewRuntime = null,
  requestFeatureLegendIntent = null,
  requestBuiltInLegendIntent = null,
  buildManualLegendCommand = undefined,
  manualCommandAdapters = undefined
}) => {
  const layoutActions = createLegendLayoutActions({ state });
  const entryActions = createLegendEntryActions({
    state,
    getPyodide,
    ensurePyodide,
    layoutActions,
    previewRuntime,
    history,
    requestFeatureLegendIntent,
    requestBuiltInLegendIntent,
    ...(buildManualLegendCommand ? { buildManualLegendCommand } : {}),
    ...(manualCommandAdapters ? { manualCommandAdapters } : {})
  });
  const sortActions = createLegendSortActions({
    state,
    requestLegendOrderIntent: entryActions.requestManualLegendIntent
  });
  const dragActions = createLegendDragActions({
    state,
    extractLegendEntries: entryActions.extractLegendEntries,
    history
  });

  return {
    ...entryActions,
    ...layoutActions,
    ...sortActions,
    ...dragActions,
    getAllFeatureLegendGroups,
    getVisibleFeatureLegendGroup,
    isCurrentLegendHorizontal
  };
};
