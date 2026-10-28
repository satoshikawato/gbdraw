// @ts-check
// The operations on the drawing list (PD-OI-090, R11): New, Rename, Move, and
// Delete. Each is one project step of the History router (`runListStep`), and
// these are the only callers of the drawing store's list writers. A step that
// adds or deletes a drawing shows the drawing shown after it, and its Undo the
// drawing shown before (R16-2). A step that takes a drawing out holds it, its
// Result slot, and (through the router) its History, for its Undo (R16-3).
import { allocateDrawingId, defaultDrawingName, drawingNamesProblem } from '../services/drawing-view.js';
import { diagnosticError, normalizeUserFacingError } from '../utils/error-normalization.js';
import { errorText } from '../utils/error-text.js';
import { msg, t, tm } from '../utils/i18n.js';

/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { ModeProfileName } from '../mode-profiles.js' */
/** @import { DrawingState, DrawingTab } from '../state.js' */
/** @import { ProjectListStep } from '../services/history.js' */
/** @import { ArtifactSlot } from '../services/artifact-slot.js' */
/** @import { AskChoice, AskText } from './app-messages.js' */
/** @import { DiagnosticError } from '../utils/error-normalization.js' */

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

/**
 * Whether `error` is the list's refusal of a name (blank or repeated), which
 * the Rename dialog shows and asks again; any other error is not its to show.
 * @param {unknown} error @returns {error is DiagnosticError}
 */
const isNameRefusal = (error) => error instanceof Error && 'stage' in error && error.stage === 'drawing-list'
  && 'context' in error && typeof error.context === 'object' && error.context !== null && 'reason' in error.context
  && (error.context.reason === 'BLANK_NAME' || error.context.reason === 'REPEATED_NAME');

/** The UI language's word for each mode, for a default name in the bar. */
const MODE_LABELS = Object.freeze({ circular: msg('Circular'), linear: msg('Linear') });

/**
 * One tab of the drawing bar. `label` is the name as the tab shows it: a
 * default name ("Circular", "Linear 2") in the UI language, any other name as
 * stored. The tab's accessible name is its label.
 * @typedef {object} DrawingBarTab
 * @property {DrawingId} id
 * @property {string} label
 * @property {ModeProfileName} mode
 * @property {boolean} selected the shown drawing
 * @property {boolean} canMoveLeft
 * @property {boolean} canMoveRight
 * @property {boolean} canDelete false for the last drawing, which a Session always keeps (R11)
 */
/**
 * @typedef {(drawingId: DrawingId) => ModeProfileName} ReadDrawingMode
 * @typedef {(drawingId: DrawingId, name: string) => DrawingListOutcome} RenameDrawingStep
 *   The list's Rename: one History step, or a typed refusal (`stage: 'drawing-list'`).
 * @typedef {(drawingId: DrawingId, index: number) => DrawingListOutcome} MoveDrawingStep The list's Move.
 * @typedef {(drawingId: DrawingId) => DrawingListOutcome} DeleteDrawingStep The list's Delete.
 * @typedef {(drawingId: DrawingId) => boolean} NeedsDeleteConfirmation
 * @typedef {object} DrawingBarPorts
 * @property {RenameDrawingStep} renameDrawing
 * @property {MoveDrawingStep} moveDrawing
 * @property {DeleteDrawingStep} deleteDrawing
 * @property {NeedsDeleteConfirmation} needsDeleteConfirmation
 * @property {ReadDrawingTabs} drawingTabs
 * @property {ReadShownDrawingId} shownDrawingId
 * @property {ReadDrawingMode} modeOf
 * @property {AskText} askText
 * @property {AskChoice} askChoice
 */

/**
 * The tab the arrow keys focus next in a tablist of `count` tabs, wrapping at
 * the ends, or `null` for a key the tab pattern leaves to the button: Enter
 * and Space select the focused tab (manual activation), so moving the focus
 * never switches a drawing.
 * @param {string} key @param {number} index @param {number} count
 * @returns {number | null}
 */
export const focusedTabIndex = (key, index, count) => {
  if (key === 'ArrowRight') return (index + 1) % count;
  if (key === 'ArrowLeft') return (index - 1 + count) % count;
  if (key === 'Home') return 0;
  if (key === 'End') return count - 1;
  return null;
};

/**
 * The tablist's keydown: the arrow keys, Home and End move the focus among
 * the enabled tabs (the ARIA tab pattern with manual activation).
 * @param {KeyboardEvent} event
 */
export const moveTabFocus = (event) => {
  const tablist = event.currentTarget;
  if (!(tablist instanceof HTMLElement)) return;
  const tabs = /** @type {HTMLButtonElement[]} */ (Array.from(tablist.querySelectorAll('[role="tab"]:not(:disabled)')));
  const at = tabs.findIndex((tab) => tab === document.activeElement);
  const next = at < 0 ? null : focusedTabIndex(event.key, at, tabs.length);
  if (next === null) return;
  event.preventDefault();
  tabs[next].focus();
};

/**
 * The drawing bar's view (PD-OI-090 〔drawing bar と切り替え〕): its tabs, and
 * the tab menu's Rename, Move, and Delete, each one list operation (one History step).
 * @param {DrawingBarPorts} ports
 */
export const createDrawingBar = ({
  renameDrawing: renameStep, moveDrawing: moveStep, deleteDrawing: deleteStep, needsDeleteConfirmation,
  drawingTabs, shownDrawingId, modeOf, askText, askChoice
}) => {
  /** @returns {DrawingBarTab[]} */
  const tabs = () => {
    const all = drawingTabs();
    const shown = shownDrawingId();
    return all.map(({ id, name }, index) => {
      const mode = modeOf(id);
      return {
        id,
        mode,
        label: name === defaultDrawingName(mode, id) ? defaultDrawingName(mode, id, tm(MODE_LABELS[mode])) : name,
        selected: id === shown,
        canMoveLeft: index > 0,
        canMoveRight: index < all.length - 1,
        canDelete: all.length > 1
      };
    });
  };

  /**
   * Rename through the app's text dialog. A refused name (blank, or another
   * drawing's) asks again with the error owner's text; Cancel renames nothing.
   * @param {DrawingId} drawingId
   * @returns {Promise<DrawingListOutcome | null>}
   */
  const renameDrawing = async (drawingId) => {
    let value = drawingTabs().find((tab) => tab.id === drawingId)?.name ?? '';
    let description = '';
    for (;;) {
      const input = await askText({
        title: t('Rename drawing'), description, label: t('Drawing name'), value, actionLabel: t('Rename')
      });
      if (input === null) return null;
      try {
        return renameStep(drawingId, input);
      } catch (error) {
        if (!isNameRefusal(error)) throw error;
        value = input;
        // A refusal always normalizes (null only for a null input).
        const refusal = normalizeUserFacingError(error);
        description = refusal ? errorText(refusal.summaryParts) : '';
      }
    }
  };

  /**
   * One place left (-1) or right (+1) in the tab order.
   * @param {DrawingId} drawingId @param {-1 | 1} step
   * @returns {DrawingListOutcome}
   */
  const moveDrawing = (drawingId, step) => moveStep(
    drawingId, drawingTabs().findIndex((tab) => tab.id === drawingId) + step
  );

  /**
   * Delete, asking first only for a drawing with a Result or edits (R17-6);
   * Undo brings any deleted drawing back (R16-3). Cancel deletes nothing and
   * records no step.
   * @param {DrawingId} drawingId
   * @returns {Promise<DrawingListOutcome | null>}
   */
  const deleteDrawing = async (drawingId) => {
    if (needsDeleteConfirmation(drawingId)) {
      const name = tabs().find((tab) => tab.id === drawingId)?.label ?? drawingId;
      const choice = await askChoice({
        title: t('Delete {name}?', { name }),
        description: t('Its settings, edits, and Result are removed from this Session. Undo brings the drawing back.'),
        actions: [{ id: 'delete', label: t('Delete drawing'), tone: 'danger' }],
        cancelLabel: t('Keep drawing')
      });
      if (choice !== 'delete') return null;
    }
    return deleteStep(drawingId);
  };

  return { tabs, renameDrawing, moveDrawing, deleteDrawing };
};
