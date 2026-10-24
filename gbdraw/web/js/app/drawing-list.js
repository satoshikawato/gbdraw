// @ts-check
// The operations on the drawing list (PD-OI-090, R11): New, Rename, and Move.
// Each is one project step of the History router (`runListStep`), and these
// are the only callers of the drawing store's list writers. A step that adds
// a drawing shows it, and its Undo shows the drawing shown before (R16-2).
import { allocateDrawingId, defaultDrawingName, drawingNamesProblem } from '../services/drawing-view.js';
import { diagnosticError } from '../utils/error-normalization.js';

/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { ModeProfileName } from '../mode-profiles.js' */
/** @import { DrawingState, DrawingTab } from '../state.js' */
/** @import { ProjectListStep } from '../services/history.js' */

/**
 * @typedef {() => readonly DrawingTab[]} ReadDrawingTabs The drawings' IDs and names, in tab order.
 * @typedef {() => DrawingId} ReadShownDrawingId
 * @typedef {() => unknown} ReadSwitchBusy Why a drawing switch waits now, or a falsy value (R15).
 * @typedef {(drawingId: DrawingId, mode: ModeProfileName) => DrawingState} CreateDrawing A new drawing at its defaults.
 * @typedef {(drawing: DrawingState, name: string, index: number) => void} AddDrawingToStore
 * @typedef {(drawingId: DrawingId) => void} RemoveDrawingFromStore Removes a drawing that is not shown, with its slot.
 * @typedef {(drawingId: DrawingId, name: string) => void} RenameDrawingInStore
 * @typedef {(drawingId: DrawingId, index: number) => void} MoveDrawingInStore
 * @typedef {(drawingId: DrawingId) => void} AddDrawingHistory An empty History for a drawing the store holds (R11).
 * @typedef {(label: string, step: ProjectListStep) => void} RunListStep One project step: applies the step now.
 * @typedef {(drawingId: DrawingId) => void} ShowDrawingNow The drawing switch's work, without its gate.
 * @typedef {{ status: 'ok', drawingId: DrawingId } | { status: 'busy', reason: unknown }} DrawingListOutcome
 *   A list operation that ran, or the reason it waits (as a drawing switch waits).
 */

/**
 * @typedef {object} DrawingListPorts
 * @property {ReadDrawingTabs} drawingTabs
 * @property {ReadShownDrawingId} shownDrawingId
 * @property {ReadSwitchBusy} switchBusy
 * @property {CreateDrawing} createDrawing
 * @property {AddDrawingToStore} addDrawing
 * @property {RemoveDrawingFromStore} removeDrawing
 * @property {RenameDrawingInStore} renameDrawing
 * @property {MoveDrawingInStore} moveDrawing
 * @property {AddDrawingHistory} addDrawingHistory
 * @property {RunListStep} runListStep
 * @property {ShowDrawingNow} showDrawingNow
 */

/** @param {DrawingListPorts} ports */
export const createDrawingList = ({
  drawingTabs, shownDrawingId, switchBusy, createDrawing, addDrawing, removeDrawing, renameDrawing: renameInStore,
  moveDrawing: moveInStore, addDrawingHistory, runListStep, showDrawingNow
}) => {
  /** @param {string} operation @param {string} reason @param {DrawingId} drawingId */
  const listError = (operation, reason, drawingId) => diagnosticError(
    'INPUT_INVALID', { field: 'drawing', reason, drawingId }, { stage: 'drawing-list', operation }
  );
  /** @param {string} operation @param {DrawingId} drawingId @returns {number} */
  const tabIndexOf = (operation, drawingId) => {
    const index = drawingTabs().findIndex((tab) => tab.id === drawingId);
    if (index < 0) throw listError(operation, 'UNKNOWN_DRAWING', drawingId);
    return index;
  };
  /**
   * @param {() => DrawingId} run The operation's first application, through `runListStep`.
   * @returns {DrawingListOutcome}
   */
  const unlessBusy = (run) => {
    const reason = switchBusy();
    if (reason) return { status: 'busy', reason };
    return { status: 'ok', drawingId: run() };
  };

  /**
   * A new drawing of `mode` at the end of the tabs, with the ID and name of
   * its kind (`linear-2` is "Linear 2") and an empty History, shown (R11).
   * Its Undo shows the drawing shown before; its Redo shows it again (R16-2).
   * @param {ModeProfileName} mode
   */
  const newDrawing = (mode) => unlessBusy(() => {
    const before = shownDrawingId();
    const drawingId = allocateDrawingId(mode, drawingTabs().map((tab) => tab.id));
    const drawing = createDrawing(drawingId, mode);
    const name = defaultDrawingName(mode, drawingId);
    // Redo gets the drawing's History back from the router, which holds it.
    runListStep('New drawing', {
      apply: () => {
        addDrawing(drawing, name, drawingTabs().length);
        addDrawingHistory(drawingId);
        showDrawingNow(drawingId);
      },
      redo: () => {
        addDrawing(drawing, name, drawingTabs().length);
        showDrawingNow(drawingId);
      },
      undo: () => {
        showDrawingNow(before);
        removeDrawing(drawingId);
      }
    });
    return drawingId;
  });

  /**
   * A drawing's new name: refused when blank or when another drawing has it,
   * whatever its letter case (Python's rule for a Session's drawings).
   * @param {DrawingId} drawingId @param {string} name
   * @returns {DrawingListOutcome}
   */
  const renameDrawing = (drawingId, name) => {
    const tabs = drawingTabs();
    const index = tabIndexOf('renameDrawing', drawingId);
    const previous = tabs[index].name;
    if (name === previous) return { status: 'ok', drawingId };
    const problem = drawingNamesProblem(tabs.map((tab, i) => (i === index ? name : tab.name)));
    if (problem) throw listError('renameDrawing', problem, drawingId);
    return unlessBusy(() => {
      runListStep('Rename drawing', {
        redo: () => renameInStore(drawingId, name),
        undo: () => renameInStore(drawingId, previous)
      });
      return drawingId;
    });
  };

  /**
   * A drawing's tab at `index` of the tab order.
   * @param {DrawingId} drawingId @param {number} index @returns {DrawingListOutcome}
   */
  const moveDrawing = (drawingId, index) => {
    const from = tabIndexOf('moveDrawing', drawingId);
    if (index === from) return { status: 'ok', drawingId };
    return unlessBusy(() => {
      runListStep('Move drawing', {
        redo: () => moveInStore(drawingId, index),
        undo: () => moveInStore(drawingId, from)
      });
      return drawingId;
    });
  };

  return { newDrawing, renameDrawing, moveDrawing };
};
