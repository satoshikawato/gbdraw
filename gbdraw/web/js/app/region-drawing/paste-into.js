// @ts-check
// The existing drawings Copy to another drawing pastes into (PD-OI-090 〔作成〕(4)): their draft requests, and the
// paste as one History step of the drawing it goes into. The dialog (app/region-drawing.js) reads them as its
// `pasteTargets` port; the app passes the owners (the drawing store, the record read, the projection, the drawing
// switch, History, and the paste's reader and install in services/config.js).
import { msg } from '../../utils/i18n.js';

/** @import { DrawingTab, DrawingState } from '../../state.js' */
/** @import { DrawingId } from '../../services/drawing-view.js' */
/** @import { DrawingAlignment, PasteDerivation, PreparedPaste, ProjectedDraftRequest } from '../../services/config.js' */
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
 * @property {(label: string, run: () => void) => unknown} runUndoable One History step of the shown drawing; while a
 *   Session operation blocks edits it runs nothing and answers `{ status: 'busy', reason }`.
 * @property {(paste: PasteDerivation) => Promise<PreparedPaste>} prepare `preparePastedRecords`
 * @property {(drawing: DrawingState, paste: PreparedPaste) => DrawingAlignment} install `installPastedRecords`: it
 *   answers the drawing's alignment after the paste.
 * @property {(drawingId: DrawingId, values: DrawingAlignment) => void} writeArtifactValues Writes a drawing's own
 *   artifact values (History's writer: the page's while it is shown).
 */

/**
 * The refusal in an answer of the drawing switch or of History (`{ status }` other than `ok`), as the drawing list
 * answers one; null for any other answer.
 * @param {unknown} answer
 * @returns {DrawingListOutcome | null}
 */
const refusalOf = (answer) => (answer && typeof answer === 'object' && 'status' in answer && answer.status !== 'ok'
  ? { status: 'busy', reason: 'reason' in answer ? answer.reason : null } : null);

/**
 * The dialog's port to the existing drawings. A paste is read first, then installed, with its drawing shown, as one
 * step of that drawing's History; its figure is drawn by the next Generate (the drawing keeps its Result until then).
 * A switch or a History step the app refuses (a Generate runs) pastes nothing and answers `busy` with its reason.
 * @param {PasteIntoOwners} owners
 * @returns {PasteTargetPorts}
 */
export const createPasteTargets = ({ tabs, drawingOf, readRecords, project, show, runUndoable, prepare, install, writeArtifactValues }) => ({
  tabs,
  drawingOf,
  projectDrawing: async (drawingId) => {
    await readRecords(drawingOf(drawingId));
    return project(drawingOf(drawingId));
  },
  pasteInto: async (drawingId, paste) => {
    const prepared = await prepare(paste);
    const refused = refusalOf(show(drawingId)) ?? refusalOf(await runUndoable(msg('Copy to another drawing'),
      () => writeArtifactValues(drawingId, install(drawingOf(drawingId), prepared))));
    return refused ?? { status: 'ok', drawingId };
  }
});
