// @ts-check
// The operations on the drawing list (PD-OI-090, R11): New, New with a style, Duplicate, Rename,
// Move, and Delete. Each is one project step of the History router (`runListStep`), and
// these are the only callers of the drawing store's list writers. A step that
// adds or deletes a drawing shows the drawing shown after it, and its Undo the
// drawing shown before (R16-2). A step that takes a drawing out holds it, its
// Result slot, and (through the router) its History, for its Undo (R16-3).
import { allocateDrawingId, copyDrawingName, defaultDrawingName, drawingNamesProblem } from '../services/drawing-view.js';
import { projectInputReplacement } from '../services/drawing/inputs.js';
import { diagnosticError, normalizeUserFacingError } from '../utils/error-normalization.js';
import { errorText } from '../utils/error-text.js';
import { msg, t, tm } from '../utils/i18n.js';

/** @import { CopyNamePhrase, DrawingId } from '../services/drawing-view.js' */
/** @import { ModeProfileName } from '../mode-profiles.js' */
/** @import { DrawingState, DrawingTab } from '../state.js' */
/** @import { ProjectListStep } from '../services/history.js' */
/** @import { ArtifactSlot } from '../services/artifact-slot.js' */
/** @import { AskChoice, AskText } from './app-messages.js' */
/** @import { DiagnosticError } from '../utils/error-normalization.js' */
/** @import { InputFileCandidate, InputFileRef } from '../services/history-files.js' */

/**
 * @typedef {() => readonly DrawingTab[]} ReadDrawingTabs The drawings' IDs and names, in tab order.
 * @typedef {() => DrawingId} ReadShownDrawingId
 * @typedef {{ readonly reason: string }} SwitchBusy What a busy drawing switch says; its owners add the fields the app uses.
 * @typedef {() => (SwitchBusy | null)} ReadSwitchBusy Why a drawing switch waits now, or null (R15).
 * @typedef {(drawingId: DrawingId, mode: ModeProfileName) => DrawingState} CreateDrawing A new drawing at its defaults.
 * @typedef {(drawingId: DrawingId) => DrawingState} ReadDrawing A drawing the store holds.
 * @typedef {Readonly<ArtifactSlot> | null} HeldSlot A removed drawing's Result slot (null: it had none).
 * @typedef {(sourceId: DrawingId, drawingId: DrawingId, mode: ModeProfileName) => DrawingState} StyleDrawing
 *   A new drawing of `mode` at its defaults with the style of a drawing the store holds (PD-OI-090 〔作成〕(3)).
 * @typedef {{ slot: HeldSlot }} HeldDrawing What a drawing out of the store keeps for its return: its Result slot,
 *   which holds its LOSAT evidence too (D-221), so it comes back with the values it left with.
 * @typedef {{ drawing: DrawingState, held: HeldDrawing }} MadeDrawing A new drawing and what it starts with.
 * @typedef {(sourceId: DrawingId, copyId: DrawingId) => MadeDrawing} CopyDrawing
 *   A copy of a drawing the store holds, under `copyId`, with its Result slot (null: it has no Result), which shares
 *   its source's LOSAT evidence (OV-3205).
 * @typedef {import('../services/config.js').RegionDerivation} RegionDerivation A derived region drawing and the request it was derived from.
 * @typedef {(derivation: RegionDerivation, create: { runSavedComparison: boolean }) => Promise<(drawingId: DrawingId) => DrawingState>} RegionDrawing
 *   Reads a derived region drawing; the result makes its state under the new drawing's ID. Find in other files' Create
 *   runs its saved comparison at Generate (PD-OI-123); New drawing from a region leaves the choice to the user.
 * @typedef {(drawing: DrawingState, name: string, index: number) => void} AddDrawingToStore
 * @typedef {(drawingId: DrawingId) => HeldDrawing} RemoveDrawingFromStore
 *   Removes a drawing that is not shown and hands over its slot.
 * @typedef {(held: HeldDrawing) => void} PutHeld A held slot back, with its drawing.
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
 * @property {StyleDrawing} styleDrawing
 * @property {CopyDrawing} copyDrawing
 * @property {RegionDrawing} regionDrawing
 * @property {AddDrawingToStore} addDrawing
 * @property {RemoveDrawingFromStore} removeDrawing
 * @property {PutHeld} putHeld
 * @property {RenameDrawingInStore} renameDrawing
 * @property {MoveDrawingInStore} moveDrawing
 * @property {AddDrawingHistory} addDrawingHistory
 * @property {RemoveDrawingHistory} removeDrawingHistory
 * @property {RunListStep} runListStep
 * @property {ShowDrawingNow} showDrawingNow
 */

/** @param {string} operation @param {string} reason @param {DrawingId} drawingId */
const listError = (operation, reason, drawingId) => diagnosticError(
  'INPUT_INVALID', { field: 'drawing', reason, drawingId }, { stage: 'drawing-list', operation }
);

/** The UI language's word for each mode, for a default name in the bar. */
const MODE_LABELS = Object.freeze({ circular: msg('Circular'), linear: msg('Linear') });

/**
 * A drawing's name as its tab shows it: a default name ("Circular", "Linear
 * 2") in the UI language, any other name as stored. A default name stays
 * stored in English. The tab's label and a copy's name read it (OV-2477).
 * @param {DrawingTab} tab @param {ModeProfileName} mode
 * @returns {string}
 */
const shownDrawingName = ({ id, name }, mode) => (
  name === defaultDrawingName(mode, id) ? defaultDrawingName(mode, id, tm(MODE_LABELS[mode])) : name
);

/**
 * A copy's name phrase in the UI language, stored as made (Owner D-35):
 * "<name> copy", then "<name> copy 2" in English.
 * @type {CopyNamePhrase}
 */
const copyNamePhrase = (name, number) => (
  number > 1 ? t('{name} copy {number}', { name, number }) : t('{name} copy', { name })
);

/** @type {HeldDrawing} What a drawing with no slot holds. */
const NOTHING_HELD = Object.freeze({ slot: null });
/**
 * Where a new drawing's tab goes: right after the tab of `after`, or at the end of the tabs without it.
 * @param {readonly DrawingTab[]} tabs @param {DrawingId | null} after
 */
const newTabIndex = (tabs, after) => (after === null ? tabs.length : tabs.findIndex((tab) => tab.id === after) + 1);

/** @param {DrawingListPorts} ports */
export const createDrawingList = ({
  drawingTabs, shownDrawingId, switchBusy, drawingOf, drawingHasResult, drawingHasEdits, boundFileIdsOf, createDrawing, styleDrawing, copyDrawing,
  regionDrawing, addDrawing, removeDrawing, putHeld, renameDrawing: renameInStore, moveDrawing: moveInStore, addDrawingHistory,
  removeDrawingHistory, runListStep, showDrawingNow
}) => {
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
   * A new drawing of `mode` at the end of the tabs (or after the tab of `after`), with the ID and name of
   * its kind (`linear-2` is "Linear 2") and an empty History, shown (R11).
   * Its Undo shows the drawing shown before; its Redo shows it again (R16-2).
   * @param {string} label The project step's.
   * @param {ModeProfileName} mode
   * @param {(drawingId: DrawingId) => MadeDrawing} make The new drawing's state and what it starts with.
   * @param {DrawingId | null} [after]
   * @param {(drawingId: DrawingId) => string} [nameOf]
   */
  const addNewDrawing = (label, mode, make, after = null, nameOf = (drawingId) => defaultDrawingName(mode, drawingId)) => unlessBusy(() => {
    const before = shownDrawingId();
    const drawingId = allocateDrawingId(mode, drawingTabs().map((tab) => tab.id));
    const { drawing, held: made } = make(drawingId);
    let held = made;
    const name = nameOf(drawingId);
    // Redo gets the drawing's History back from the router, which holds it, and its slot (a Result and LOSAT
    // evidence it got before the Undo) from this step.
    const add = (first = false) => {
      addDrawing(drawing, name, newTabIndex(drawingTabs(), after));
      if (first) addDrawingHistory(drawingId);
      putHeld(held);
      showDrawingNow(drawingId);
    };
    runListStep(label, {
      apply: () => add(true),
      redo: () => add(),
      undo: () => {
        showDrawingNow(before);
        held = removeDrawing(drawingId);
      }
    });
    return drawingId;
  });
  /** @param {DrawingState} drawing @returns {MadeDrawing} */
  const startingEmpty = (drawing) => ({ drawing, held: NOTHING_HELD });

  /** @param {ModeProfileName} mode */
  const newDrawing = (mode) => addNewDrawing('New drawing', mode, (drawingId) => startingEmpty(createDrawing(drawingId, mode)));

  /**
   * A new drawing of `mode` with the style of a drawing (PD-OI-090 〔作成〕(3),
   * R17-4: its settings that `mode` has and the table rows the user entered),
   * as New adds one.
   * @param {DrawingId} sourceId
   * @param {ModeProfileName} mode
   * @returns {DrawingListOutcome}
   */
  const newDrawingWithStyle = (sourceId, mode) => {
    tabIndexOf('newDrawingWithStyle', sourceId);
    return addNewDrawing('New drawing with this style', mode, (drawingId) => startingEmpty(styleDrawing(sourceId, drawingId, mode)));
  };

  /**
   * A copy of a drawing at the end of the tabs (PD-OI-090 〔作成〕): the ID of
   * its kind, the source tab's shown name with the copy phrase in the UI
   * language ("<name> copy", `copyDrawingName`), the source's files, settings,
   * edits, Result and LOSAT results (OV-3205), and an empty History, shown (R11). Its Undo
   * shows the drawing shown before; its Redo shows the copy again (R16-2).
   * @param {DrawingId} sourceId
   * @returns {DrawingListOutcome}
   */
  const duplicateDrawing = (sourceId) => {
    const tabs = drawingTabs();
    const index = tabIndexOf('duplicateDrawing', sourceId);
    const { mode } = drawingOf(sourceId);
    return addNewDrawing('Duplicate drawing', mode, (drawingId) => copyDrawing(sourceId, drawingId), null,
      () => copyDrawingName(shownDrawingName(tabs[index], mode), tabs.map((tab) => tab.name), copyNamePhrase));
  };

  /**
   * A new drawing from a region of a drawing (PD-OI-090 〔作成〕(4)), or of Find in other files (`search`, PD-OI-123:
   * right after the source's tab): the derived drawing read as Load reads a stored one, added as New adds one. The
   * source drawing (where the region was chosen) does not change.
   * @param {DrawingId} sourceId @param {RegionDerivation} derivation @param {'region' | 'search'} [kind]
   * @returns {Promise<DrawingListOutcome>}
   */
  const newDrawingFromRegion = async (sourceId, derivation, kind = 'region') => {
    tabIndexOf('newDrawingFromRegion', sourceId);
    const make = await regionDrawing(derivation, { runSavedComparison: kind === 'search' });
    return addNewDrawing(kind === 'region' ? 'New drawing from a region' : 'Find in other files', derivation.derived.mode,
      (drawingId) => startingEmpty(make(drawingId)), kind === 'region' ? null : sourceId);
  };

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
   * edits, Result slot (with its LOSAT evidence), and History until the step is evicted (R16-3). The
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
      /** @type {HeldDrawing} */
      let held = NOTHING_HELD;
      // The step keeps the files only this drawing binds until it is evicted.
      const fileIds = [...boundFileIdsOf(drawingId)];
      // The switch settles the drawing's open step while its History is still a member.
      runListStep('Delete drawing', {
        fileIds: () => fileIds,
        apply: () => {
          showDrawingNow(after);
          removeDrawingHistory(drawingId);
          held = removeDrawing(drawingId);
        },
        redo: () => {
          showDrawingNow(after);
          held = removeDrawing(drawingId);
        },
        undo: () => {
          addDrawing(drawing, name, index);
          putHeld(held);
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

  return { newDrawing, newDrawingWithStyle, newDrawingFromRegion, duplicateDrawing, renameDrawing, moveDrawing, deleteDrawing, needsDeleteConfirmation };
};

/**
 * @typedef {import('../services/drawing/inputs.js').InputReplacementSite<DrawingId>} DrawingInputSite
 * @typedef {(sites: readonly DrawingInputSite[], fileId: string, file: InputFileCandidate | InputFileRef) => void} ApplyInputReplacementToDrawings
 *   Writes `file` in place of the file `fileId` at each site, in the site's drawing, with that
 *   drawing's reactions to its new source; every slot a site names has a writer.
 * @typedef {(label: string, apply: () => void) => unknown} RunReplacementStep One History step that runs `apply`.
 * @typedef {'drawing' | 'all'} InputReplacementScope The shown drawing only, or every drawing that uses the file.
 * @typedef {object} InputReplacementPorts
 * @property {ReadDrawingTabs} drawingTabs
 * @property {ReadShownDrawingId} shownDrawingId
 * @property {ReadSwitchBusy} switchBusy
 * @property {ReadDrawing} drawingOf
 * @property {RunReplacementStep} runProjectStep A project step: one Undo from any drawing (R11).
 * @property {RunReplacementStep} runDrawingStep A step of the shown drawing.
 * @property {ApplyInputReplacementToDrawings} applyInputReplacementToDrawings
 * @property {(previous: InputFileRef, otherDrawingIds: readonly DrawingId[]) => Promise<InputReplacementScope | null>} askScope
 *   Asks "This drawing only" / "Every drawing that uses this file"; null is Cancel.
 * @property {(busy: SwitchBusy) => void} onBusy Shows why "Every drawing" wrote nothing: the drawing switch waits (D-15-186).
 */

/**
 * The replacement of an input file that drawings bind (PD-OI-090 〔入力ファイル〕,
 * R16-5, R17-3, R17-30): one project step, after which no drawing binds the
 * replaced file. The sites come from `projectInputReplacement` by drawing ID
 * and slot, row uid, or comparison id; the five cards are never sites.
 * @param {InputReplacementPorts} ports
 */
export const createInputReplacement = ({
  drawingTabs, shownDrawingId, switchBusy, drawingOf, runProjectStep, runDrawingStep, applyInputReplacementToDrawings, askScope, onBusy
}) => {
  /** @param {string} fileId @returns {DrawingInputSite[]} every site of `fileId`, in every drawing */
  const sitesOf = (fileId) => projectInputReplacement(
    drawingTabs().map((tab) => ({ drawingId: tab.id, drawing: drawingOf(tab.id) })), fileId
  );
  /**
   * Replaces `fileId` with `file` in every drawing that uses it. `writeShown`,
   * the uploader's own write of the control the user chose the file in, runs
   * first in the same step (R17-30); the replacement port swaps the file at
   * every site still bound after it (R17-28: only the file changes). While the
   * drawing switch waits it writes nothing, not even `writeShown`, and `onBusy`
   * shows the reason, since the answer to the question has no control to
   * disable (D-15-186).
   * @param {string} fileId the bound file being replaced
   * @param {InputFileCandidate | InputFileRef} file the new file, or a file in this Session
   * @param {() => unknown} writeShown
   * @returns {DrawingListOutcome}
   */
  const replaceInputFile = (fileId, file, writeShown) => {
    const reason = switchBusy();
    if (reason) {
      onBusy(reason);
      return { status: 'busy', reason };
    }
    if (sitesOf(fileId).length) {
      runProjectStep(msg('Replace input file'), () => {
        writeShown();
        applyInputReplacementToDrawings(sitesOf(fileId), fileId, file);
      });
    }
    return { status: 'ok', drawingId: shownDrawingId() };
  };
  /**
   * The drawings other than the shown one that bind `fileId` where a
   * replacement writes: the replacement asks only when there is one.
   * @param {string} fileId @returns {DrawingId[]}
   */
  const otherDrawingsUsing = (fileId) => {
    const shown = shownDrawingId();
    return [...new Set(sitesOf(fileId).map((site) => site.drawingId).filter((drawingId) => drawingId !== shown))];
  };
  /**
   * An uploader of the shown drawing writes a new file in place of `previous`
   * (PD-OI-090 〔入力ファイル〕). Its control always gets the uploader's own
   * write `writeShown` (R17-30), in the uploader's step when nothing is asked.
   * When another drawing uses `previous`, the scope is asked first, after
   * which the uploader's step has ended: Cancel records nothing; "This drawing
   * only" is that write alone, a step of the shown drawing; "Every drawing that
   * uses this file" is `replaceInputFile` with it, which swaps every other
   * site. A Remove (`file` null) is not a replacement and never asks. The only
   * place that asks.
   * @template T
   * @param {InputFileRef | null} previous the file the uploader binds now
   * @param {InputFileCandidate | InputFileRef | null} file
   * @param {() => T} writeShown
   * @returns {T | Promise<unknown>}
   */
  const replaceFromShownDrawing = (previous, file, writeShown) => {
    const others = previous && file && previous !== file ? otherDrawingsUsing(previous.fileId) : [];
    if (!previous || !file || !others.length) return writeShown();
    return askScope(previous, others).then((scope) => {
      if (scope === 'drawing') return runDrawingStep(msg('Replace input file'), writeShown);
      return scope === 'all' ? replaceInputFile(previous.fileId, file, writeShown) : null;
    });
  };
  return { replaceInputFile, otherDrawingsUsing, replaceFromShownDrawing };
};

/**
 * Whether `error` is the list's refusal of a name (blank or repeated), which
 * the Rename dialog shows and asks again; any other error is not its to show.
 * @param {unknown} error @returns {error is DiagnosticError}
 */
const isNameRefusal = (error) => error instanceof Error && 'stage' in error && error.stage === 'drawing-list'
  && 'context' in error && typeof error.context === 'object' && error.context !== null && 'reason' in error.context
  && (error.context.reason === 'BLANK_NAME' || error.context.reason === 'REPEATED_NAME');

/**
 * One tab of the drawing bar. `label` is the name as the tab shows it: a
 * default name ("Circular", "Linear 2") in the UI language, any other name as
 * stored. The tab's accessible name is its label.
 * @typedef {object} DrawingBarTab
 * @property {'drawing'} kind
 * @property {DrawingId} id
 * @property {string} label
 * @property {ModeProfileName} mode
 * @property {boolean} selected the shown drawing
 * @property {boolean} canMoveLeft
 * @property {boolean} canMoveRight
 * @property {boolean} canDelete false for the last drawing, which a Session always keeps (R11)
 */
/**
 * The Linear offer tab (PD-OI-090: a new Session starts with one Circular
 * drawing, and Linear is made from a tab that offers an empty Linear drawing).
 * It follows the drawing tabs while the project has no Linear drawing (T-B4-5);
 * choosing it is New Linear.
 * @typedef {object} DrawingBarOffer
 * @property {'offer'} kind
 * @property {ModeProfileName} mode
 * @property {string} label the mode's name in the UI language
 */
/**
 * @typedef {(drawingId: DrawingId) => ModeProfileName} ReadDrawingMode
 * @typedef {(drawingId: DrawingId, name: string) => DrawingListOutcome} RenameDrawingStep
 *   The list's Rename: one History step, or a typed refusal (`stage: 'drawing-list'`).
 * @typedef {(drawingId: DrawingId, index: number) => DrawingListOutcome} MoveDrawingStep The list's Move.
 * @typedef {(drawingId: DrawingId) => DrawingListOutcome} DeleteDrawingStep The list's Delete.
 * @typedef {(drawingId: DrawingId) => boolean} NeedsDeleteConfirmation
 * @typedef {(mode: ModeProfileName) => DrawingListOutcome} NewDrawingStep The list's New.
 * @typedef {(drawingId: DrawingId) => DrawingListOutcome} DuplicateDrawingStep The list's Duplicate.
 * @typedef {(drawingId: DrawingId, mode: ModeProfileName) => DrawingListOutcome} NewDrawingWithStyleStep
 *   The list's New with the style of a drawing, as a drawing of the target mode.
 * @typedef {object} DrawingBarPorts
 * @property {NewDrawingStep} newDrawing
 * @property {DuplicateDrawingStep} duplicateDrawing
 * @property {NewDrawingWithStyleStep} newDrawingWithStyle
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
 * The drawing bar's view (PD-OI-090 〔drawing bar と切り替え〕): its tabs with
 * the Linear offer tab, and the tab menu's New Circular, New Linear, Duplicate,
 * New Circular or Linear drawing with this style (R17-27), Rename, Move, and
 * Delete, each one list operation (one History step).
 * @param {DrawingBarPorts} ports
 */
export const createDrawingBar = ({
  newDrawing, duplicateDrawing, newDrawingWithStyle, renameDrawing: renameStep, moveDrawing: moveStep, deleteDrawing: deleteStep, needsDeleteConfirmation,
  drawingTabs, shownDrawingId, modeOf, askText, askChoice
}) => {
  /** @returns {Array<DrawingBarTab | DrawingBarOffer>} */
  const tabs = () => {
    const all = drawingTabs();
    const shown = shownDrawingId();
    /** @type {Array<DrawingBarTab | DrawingBarOffer>} */
    const bar = all.map((tab, index) => {
      const { id } = tab;
      const mode = modeOf(id);
      return {
        kind: /** @type {const} */ ('drawing'),
        id,
        mode,
        label: shownDrawingName(tab, mode),
        selected: id === shown,
        canMoveLeft: index > 0,
        canMoveRight: index < all.length - 1,
        canDelete: all.length > 1
      };
    });
    if (!bar.some((tab) => tab.mode === 'linear')) bar.push({ kind: 'offer', mode: 'linear', label: tm(MODE_LABELS.linear) });
    return bar;
  };

  /**
   * The bar's tab of a drawing, by its ID; an unknown ID is a bug.
   * @param {string} operation @param {DrawingId} drawingId
   * @returns {DrawingBarTab}
   */
  const drawingTabOf = (operation, drawingId) => {
    const tab = tabs().find((each) => each.kind === 'drawing' && each.id === drawingId);
    if (!tab || tab.kind !== 'drawing') throw listError(operation, 'UNKNOWN_DRAWING', drawingId);
    return tab;
  };

  /**
   * Rename through the app's text dialog, starting from the tab's shown
   * name (OV-2477). A refused name (blank, or another
   * drawing's) asks again with the error owner's text; Cancel renames nothing.
   * @param {DrawingId} drawingId
   * @returns {Promise<DrawingListOutcome | null>}
   */
  const renameDrawing = async (drawingId) => {
    let value = drawingTabOf('renameDrawing', drawingId).label;
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
      const name = drawingTabOf('deleteDrawing', drawingId).label;
      const choice = await askChoice({
        title: t('Delete {name}?', { name }),
        lines: [
          { icon: 'ph-warning', text: t('Its edits and Result are removed') },
          { icon: 'ph-arrow-counter-clockwise', text: t('Undo restores the drawing') }
        ],
        actions: [{ id: 'delete', label: t('Delete drawing'), tone: 'danger' }],
        cancelLabel: t('Keep drawing')
      });
      if (choice !== 'delete') return null;
    }
    return deleteStep(drawingId);
  };

  return { tabs, newDrawing, duplicateDrawing, newDrawingWithStyle, renameDrawing, moveDrawing, deleteDrawing };
};
