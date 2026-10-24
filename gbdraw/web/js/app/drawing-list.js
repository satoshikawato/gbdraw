// @ts-check
// The operations on the drawing list (PD-OI-090, R11): New, Rename, Move, and
// Delete. Each is one project step of the History router (`runListStep`), and
// these are the only callers of the drawing store's list writers. A step that
// adds or deletes a drawing shows the drawing shown after it, and its Undo the
// drawing shown before (R16-2). A step that takes a drawing out holds it, its
// Result slot, and (through the router) its History, for its Undo (R16-3).
import { allocateDrawingId, defaultDrawingName, drawingNamesProblem } from '../services/drawing-view.js';
import { diagnosticError } from '../utils/error-normalization.js';

/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { ModeProfileName } from '../mode-profiles.js' */
/** @import { DrawingState, DrawingTab } from '../state.js' */
/** @import { ProjectListStep } from '../services/history.js' */
/** @import { ArtifactSlot } from '../services/artifact-slot.js' */

/**
 * @typedef {() => readonly DrawingTab[]} ReadDrawingTabs The drawings' IDs and names, in tab order.
 * @typedef {() => DrawingId} ReadShownDrawingId
 * @typedef {() => unknown} ReadSwitchBusy Why a drawing switch waits now, or a falsy value (R15).
 * @typedef {(drawingId: DrawingId, mode: ModeProfileName) => DrawingState} CreateDrawing A new drawing at its defaults.
 * @typedef {(drawingId: DrawingId) => DrawingState} ReadDrawing A drawing the store holds.
 * @typedef {Readonly<ArtifactSlot> | null} HeldSlot A removed drawing's Result slot (null: it had none).
 * @typedef {(drawing: DrawingState, name: string, index: number) => void} AddDrawingToStore
 * @typedef {(drawingId: DrawingId) => HeldSlot} RemoveDrawingFromStore
 *   Removes a drawing that is not shown, and hands over its slot.
 * @typedef {(slot: Readonly<ArtifactSlot>) => void} PutSlot A held slot back, with its drawing.
 * @typedef {(drawingId: DrawingId) => boolean} ReadDrawingHas Whether a drawing has a Result, or edits.
 * @typedef {(drawingId: DrawingId) => Iterable<string>} ReadBoundFileIds The files a drawing binds.
 * @typedef {(drawingId: DrawingId, name: string) => void} RenameDrawingInStore
 * @typedef {(drawingId: DrawingId, index: number) => void} MoveDrawingInStore
 * @typedef {(drawingId: DrawingId) => void} AddDrawingHistory An empty History for a drawing the store holds (R11).
 * @typedef {(drawingId: DrawingId) => void} RemoveDrawingHistory
 *   Hands a drawing's History to the open list step, which holds it (R16-3).
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
 * @property {ReadDrawing} drawingOf
 * @property {ReadDrawingHas} drawingHasResult
 * @property {ReadDrawingHas} drawingHasEdits
 * @property {ReadBoundFileIds} boundFileIdsOf
 * @property {CreateDrawing} createDrawing
 * @property {AddDrawingToStore} addDrawing
 * @property {RemoveDrawingFromStore} removeDrawing
 * @property {PutSlot} putSlot
 * @property {RenameDrawingInStore} renameDrawing
 * @property {MoveDrawingInStore} moveDrawing
 * @property {AddDrawingHistory} addDrawingHistory
 * @property {RemoveDrawingHistory} removeDrawingHistory
 * @property {RunListStep} runListStep
 * @property {ShowDrawingNow} showDrawingNow
 */

/** @param {DrawingListPorts} ports */
export const createDrawingList = ({
  drawingTabs, shownDrawingId, switchBusy, drawingOf, drawingHasResult, drawingHasEdits, boundFileIdsOf, createDrawing, addDrawing,
  removeDrawing, putSlot, renameDrawing: renameInStore, moveDrawing: moveInStore, addDrawingHistory,
  removeDrawingHistory, runListStep, showDrawingNow
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
    // Redo gets the drawing's History back from the router, which holds it,
    // and its slot (a Result it got before the Undo) from this step.
    /** @type {HeldSlot} */
    let slot = null;
    runListStep('New drawing', {
      apply: () => {
        addDrawing(drawing, name, drawingTabs().length);
        addDrawingHistory(drawingId);
        showDrawingNow(drawingId);
      },
      redo: () => {
        addDrawing(drawing, name, drawingTabs().length);
        if (slot) putSlot(slot);
        showDrawingNow(drawingId);
      },
      undo: () => {
        showDrawingNow(before);
        slot = removeDrawing(drawingId);
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

  /**
   * Deletes a drawing, as one project step that holds it with its settings,
   * edits, Result slot, and History until the step is evicted (R16-3). The
   * last drawing cannot be deleted (R11). Deleting the shown drawing first
   * shows the next tab, else the previous (R16-2).
   * @param {DrawingId} drawingId @returns {DrawingListOutcome}
   */
  const deleteDrawing = (drawingId) => {
    const tabs = drawingTabs();
    const index = tabIndexOf('deleteDrawing', drawingId);
    if (tabs.length <= 1) throw listError('deleteDrawing', 'LAST_DRAWING', drawingId);
    return unlessBusy(() => {
      const before = shownDrawingId();
      const after = before === drawingId ? tabs[index + 1 < tabs.length ? index + 1 : index - 1].id : before;
      const drawing = drawingOf(drawingId);
      const { name } = tabs[index];
      /** @type {HeldSlot} */
      let slot = null;
      // The step keeps the files only this drawing binds until it is evicted.
      const fileIds = [...boundFileIdsOf(drawingId)];
      // The switch settles the drawing's open step while its History is still a member.
      runListStep('Delete drawing', {
        fileIds: () => fileIds,
        apply: () => {
          showDrawingNow(after);
          removeDrawingHistory(drawingId);
          slot = removeDrawing(drawingId);
        },
        redo: () => {
          showDrawingNow(after);
          slot = removeDrawing(drawingId);
        },
        undo: () => {
          addDrawing(drawing, name, index);
          if (slot) putSlot(slot);
          showDrawingNow(before);
        }
      });
      return drawingId;
    });
  };

  /**
   * Whether Delete asks first: only for a drawing with a Result or edits, as
   * its Undo restores any deleted drawing (R17-6).
   * @param {DrawingId} drawingId
   */
  const needsDeleteConfirmation = (drawingId) => {
    tabIndexOf('needsDeleteConfirmation', drawingId);
    return drawingHasResult(drawingId) || drawingHasEdits(drawingId);
  };

  return { newDrawing, renameDrawing, moveDrawing, deleteDrawing, needsDeleteConfirmation };
};
