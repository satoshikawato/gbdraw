import { createInputFileStore } from '../../../gbdraw/web/js/services/history-files.js';
import { sessionResourceSource } from '../../../gbdraw/web/js/services/session-resource-backing.js';
import { bindLoadedInputs } from '../../../gbdraw/web/js/services/drawing/inputs.js';

/** An input-file store as state.js builds it: refs of one Session resource share a `sourceId`. */
export const createTestInputStore = () => createInputFileStore({
  sourceOf: (file) => /** @type {{ descriptor?: object } | null} */ (sessionResourceSource(file))?.descriptor
});

/**
 * A state as the app binds it: each input file of `source` is registered in
 * the store and named by its ref (registering reads no bytes).
 */
export const boundState = (source, store = createTestInputStore()) => {
  const bound = bindLoadedInputs(store, {
    files: source.files || {}, linearSeqs: source.linearSeqs || [],
    linearComparisonPlan: source.linearComparisonPlan || { edges: [] }, circularRecordDiscovery: {}
  });
  return {
    ...source, inputFiles: store, files: bound.files, linearSeqs: bound.linearSeqs,
    ...(source.linearComparisonPlan ? { linearComparisonPlan: bound.linearComparisonPlan } : {})
  };
};
