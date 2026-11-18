// The Load rollback port a node test passes to `importSession`: the composition
// root's own `createSessionImportRollbackState`, built over the real History
// snapshot service and artifact slot store of the given `state`. The modules
// load on call, after the test has installed its window and document.

/** @import { SvgResultEntry } from '../../../gbdraw/web/js/services/svg-result-ingestion.js' */
/** @typedef {ReturnType<typeof import('../../../gbdraw/web/js/services/history.js').createHistoryRouter>} HistoryRouter */
/**
 * @param {{ drawingOf: (id: string) => any, activeDrawing: () => any, results: { value: SvgResultEntry[] } }} state The state.js `state` the Load runs on.
 * @param {Pick<HistoryRouter, 'captureMembers' | 'replaceMembers' | 'restoreMembers' | 'releaseCapturedMembers'> | null} [history] A History
 *   router whose members the Load replaces; without one the members are an empty list nobody reads.
 * @param {{ endReviewOutside: (drawingIds: readonly string[]) => void }} [alignment] The Similarity alignment owner's
 *   port; by default no review is open.
 * @param {{ keepDrawingSearch?: (drawing: any) => void, returnToDrawingSearch?: (drawing: any) => void,
 *   drawingCopyFollowers?: readonly ((drawing: any, copy: any) => void)[],
 *   slots?: ReturnType<typeof import('../../../gbdraw/web/js/services/artifact-slot.js').createArtifactSlotStore> }} [ports]
 *   The Preview search owner's port, by default no search typed; the followers of a drawing's same-ID copy, by default
 *   the CLI options' own; and the slot store of the drawings not shown, by default a new, empty one.
 */
export const sessionRollbackPort = async (
  state, history = null, alignment = { endReviewOutside: () => {} }, { slots = null, ...ports } = {}
) => {
  const { createSessionImportRollbackState } = await import('../../../gbdraw/web/js/app/app-setup.js');
  const { followDrawingCopyCliOptions } = await import('../../../gbdraw/web/js/services/config.js');
  const { createArtifactSlotStore } = await import('../../../gbdraw/web/js/services/artifact-slot.js');
  const { createHistorySnapshotService } = await import('../../../gbdraw/web/js/services/history-snapshot.js');
  const { createInputFileStore } = await import('../../../gbdraw/web/js/services/history-files.js');
  const artifactPort = createHistorySnapshotService({
    state: /** @type {any} */ (state), fileStore: createInputFileStore(), currentResults: () => state.results.value
  });
  const historyMembers = Object.freeze({ histories: [] });
  return createSessionImportRollbackState({
    artifactSlots: slots ?? createArtifactSlotStore(),
    captureHistoryMembers: () => (history ? history.captureMembers() : /** @type {any} */ (historyMembers)),
    replaceHistoryMembers: (drawingIds) => { history?.replaceMembers(drawingIds); },
    restoreHistoryMembers: (members) => { history?.restoreMembers(members); },
    releaseHistoryMembers: (members) => { history?.releaseCapturedMembers(members); },
    captureDisplayedArtifact: () => artifactPort.captureArtifactSlot(),
    installDisplayedArtifact: (slot) => artifactPort.installArtifactSlot(slot, { drawing: state.drawingOf(slot.drawingId) }),
    featureListScrollTop: { value: 0 },
    featureListScrollRef: { value: null },
    selectedPairwiseBlockOrthogroupId: { value: '' },
    alignment,
    shownDrawing: state.activeDrawing,
    keepDrawingSearch: () => {},
    returnToDrawingSearch: () => {},
    drawingCopyFollowers: [followDrawingCopyCliOptions],
    ...ports
  });
};
