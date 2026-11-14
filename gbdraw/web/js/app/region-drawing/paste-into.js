// @ts-check
// The existing drawings Copy to another drawing pastes into (PD-OI-090 〔作成〕(4)): their draft requests, and the
// paste as one History step of the drawing it goes into. The dialog (app/region-drawing.js) reads them as its
// `pasteTargets` port; the app passes the owners (the drawing store, the record read, the projection, the drawing
// switch, History, and the paste's reader and install in services/config.js).
import { msg } from '../../utils/i18n.js';

/** @import { DrawingTab, DrawingState } from '../../state.js' */
/** @import { DrawingId } from '../../services/drawing-view.js' */
/** @import { PasteDerivation, PreparedPaste, ProjectedDraftRequest } from '../../services/config.js' */
/** @import { PasteTargetPorts } from './paste.js' */
/** @import { DrawingListOutcome } from '../drawing-list.js' */

/**
 * @typedef {object} PasteIntoOwners
 * @property {() => readonly DrawingTab[]} tabs
 * @property {(drawingId: DrawingId) => DrawingState} drawingOf
 * @property {(drawing: DrawingState) => Promise<unknown>} readRecords Reads a drawing's records when a Load left them
 *   unread (a Circular drawing's projection names its records by them).
 * @property {(drawing: DrawingState) => Promise<ProjectedDraftRequest>} project A drawing's draft request, shown or not.
 * @property {(drawingId: DrawingId) => unknown} show Shows a drawing; a refused switch answers `{ status }` not `ok`.
 * @property {(label: string, run: () => void) => unknown} runUndoable One History step of the shown drawing.
 * @property {(paste: PasteDerivation) => Promise<PreparedPaste>} prepare `preparePastedRecords`
 * @property {(drawing: DrawingState, paste: PreparedPaste) => void} install `installPastedRecords`
 */

/**
 * The dialog's port to the existing drawings. A paste is read first, then installed, with its drawing shown, as one
 * step of that drawing's History; its figure is drawn by the next Generate (the drawing keeps its Result until then).
 * A switch the app refuses (a Generate runs) pastes nothing and answers `busy` with its reason.
 * @param {PasteIntoOwners} owners
 * @returns {PasteTargetPorts}
 */
export const createPasteTargets = ({ tabs, drawingOf, readRecords, project, show, runUndoable, prepare, install }) => ({
  tabs,
  drawingOf,
  projectDrawing: async (drawingId) => {
    await readRecords(drawingOf(drawingId));
    return project(drawingOf(drawingId));
  },
  pasteInto: async (drawingId, paste) => {
    const prepared = await prepare(paste);
    const shown = show(drawingId);
    if (shown && typeof shown === 'object' && 'status' in shown && shown.status !== 'ok') {
      return /** @type {DrawingListOutcome} */ ({ status: 'busy', reason: 'reason' in shown ? shown.reason : null });
    }
    await runUndoable(msg('Copy to another drawing'), () => install(drawingOf(drawingId), prepared));
    return /** @type {DrawingListOutcome} */ ({ status: 'ok', drawingId });
  }
});
