// @ts-check
// The edits of a Linear drawing's Files and record cards: the one transition
// of its record rows (`applyLinearSeqMutation`), Add File, the removal and
// clearing of a File with its dialog and focus, a File's source file
// replacement, moving a File, and a record card's selector, crop, definition,
// and subtitle. The composition root passes the state, the shown-drawing
// readers and admission, History's undoable step, the alignment owner's clear
// and stable reorder, the Depth stack's source change and loss, the record
// discovery, and the comparison plan's replacement and invalidation as ports.
// The File cards are in app/linear-file-cards.js.
import { createLinearSeq, normalizeLinearSeqList, replaceLinearSeqs } from '../services/drawing/inputs.js';
import { padDepthFileSlots } from '../services/depth-track-state.js';
import { inputSourceFiles } from '../services/input-types.js';
import { reconcileLinearComparisonPlan, setWindowLosatFilenames, windowLosatFilenamesOf } from '../services/linear-comparisons.js';
import { windowLosatFilenamesAfterRemoval } from '../services/linear-units.js';
import { placeNewLayoutWindow, planLinearSourceRowMove, reconcileLinearRecordLayout } from '../services/linear-record-layout.js';
import { addCardWindow, cardWindowEditKind, cardWindowList, editCardWindow, removeCardWindow } from '../services/linear-card-windows.js';
import {
  groupLinearSourceRecords,
  inferredDefinitionForRecord,
  isPristineLinearSource,
  linearSourceHasPrimaryInput,
  moveLinearSourceGroup,
  planLinearSourceRemoval
} from '../services/linear-sources.js';
import { linearSeqRecordKey, pruneRecordsOff, wholeRecordCardWindows } from '../services/record-draw-selection.js';
import { msg, t } from '../utils/i18n.js';
import { linearDepthLogicalWidth, linearDepthSources } from './depth-tracks.js';
import { askTrackSourceRemoval } from './track-slot-edits.js';
/** @import { DiagramMode, DrawingState, LinearSeq, SessionBusy, WebState } from '../state.js' */
/** @import { AskingStackPorts, SwitchedTrackRow } from './track-slot-edits.js' */
/** @import { FileLike } from '../services/session-resource-backing.js' */
/** @import { InputFileCandidate, InputFileRef } from '../services/history-files.js' */
/** @import { LinearSourceGroup, LinearSourceRemovalPlan } from '../services/linear-sources.js' */
/** @import { LinearRecordRow } from '../services/linear-record-layout.js' */
/** @import { DepthFileSlots, DepthSourceLoss } from '../services/depth-track-state.js' */
/** @import { SequenceRecord } from '../services/record-discovery.js' */
/** @import { RecordSelectorStatus } from './linear-record-selector.js' */
/** @import { LinearUnitPlace } from '../services/linear-units.js' */
/** @import { LinearComparisonPlan } from '../services/linear-comparisons.js' */

const { computed, nextTick } = window.Vue;

/**
 * A mode panel's mutating handler: the session's busy answer, nothing when the
 * shown drawing is of another mode, else `action` on the shown drawing.
 * @typedef {<T>(drawingMode: DiagramMode, action: (drawing: DrawingState) => T) => T | SessionBusy | undefined} ShownAction
 */
/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */
/** @typedef {() => SessionBusy | null} SessionOperationAvailability The session's busy answer, or null when an edit may run. */
/** @typedef {(label: string, fn: () => unknown) => Promise<unknown>} RunUndoable History's undoable step. */
/** @typedef {(reason: string) => unknown} ClearAlignmentForMutation The alignment owner clears its plan; `reason` is a `msg` clause. */
/** @typedef {(drawing: DrawingState) => unknown} DiscoverDrawingRecords Reads the records of the drawing's inputs. */
/**
 * The edit of a drawing's record rows that every File and record card change
 * goes through: the session's busy answer, or nothing once it is written.
 * @typedef {(drawing: DrawingState, items: readonly Partial<LinearSeq>[], options?: LinearSeqMutationOptions) => SessionBusy | undefined} ApplyLinearSeqMutation
 */
/**
 * @typedef {object} LinearSeqMutationOptions
 * @property {boolean} [preserveLosatCacheInfo] Keeps the LOSAT cache rows of the comparisons.
 * @property {LinearRecordRow[]} [layoutEntries] The record rows the new rows are laid out from.
 * @property {string} [alignmentMutation] The `msg` clause that clears the alignment plan, or 'stable-reorder'.
 * @property {string[]} [retiredRecordKeys] Cards that stay but now read another file.
 * @property {() => void} [turnOffReaching] Turns off the track rows that reach a removed Depth source.
 */
/**
 * A Linear File removal's second question (`askTrackSourceRemoval`): null when nothing asks.
 * @typedef {(drawing: DrawingState, kept: readonly Partial<LinearSeq>[], label: string, apply: (turnOffReaching: () => void) => unknown) => Promise<unknown> | null} AskLinearFileRemoval
 */
/**
 * @typedef {(drawing: DrawingState, plan: Extract<LinearSourceRemovalPlan, { allowed: true }>, intent: 'clear' | 'delete', turnOffReaching?: () => void) => SessionBusy | undefined} RemoveLinearSource
 */
/**
 * The element a closed dialog returns the focus to.
 * @typedef {EventTarget & { readonly isConnected?: boolean, focus?: () => void }} ReturnFocusTarget
 */

/**
 * @typedef {object} LinearSourceEditsPorts
 * @property {SessionOperationAvailability} sessionOperationAvailability
 * @property {ShownAction} shownAction
 * @property {<T>(drawing: DrawingState, change: () => T) => T} retireLegendStylesOfUnnamedCaptions
 *   Runs `change` and retires the Legend styles of the captions it leaves unnamed.
 * @property {(drawing: DrawingState, mutate: () => unknown) => unknown} changeLinearDepthSources
 *   Runs `mutate`, a change of the drawing's Depth sources, and keeps the Depth rows in step.
 * @property {(drawing: DrawingState, nextDepth: (DepthFileSlots | undefined)[]) => DepthSourceLoss} linearDepthSourceLoss
 *   What a Linear Depth file removal does to the stack.
 * @property {<T>(drawing: DrawingState, mutate: () => T) => T} changeRecordSources
 *   The transition of every write of a drawing's row source.
 * @property {(slotsKey: 'circular_track_slots' | 'linear_track_slots') => AskingStackPorts<SwitchedTrackRow>} askingStackPorts
 *   The ports of the question a removal asks when it leaves track rows that Generate refuses.
 * @property {(drawing: DrawingState) => { expansions: Map<string, boolean>, uploads: Set<string> }} pendingDiscoveryOf
 *   The record discoveries a drawing's rows wait for.
 * @property {(drawing: DrawingState, nextPlan: unknown, options?: { invalidate?: boolean }) => void} replaceLinearComparisonPlan
 * @property {(drawing: DrawingState, options?: { preserveLosatCacheInfo?: boolean }) => void} invalidateLinearComparisonArtifacts
 * @property {{ value: string }} linearReorderNotice
 * @property {(recordKeys: string[]) => unknown} retainAlignmentForStableReorder
 *   The alignment owner keeps its plan for the records in their new order.
 * @property {ClearAlignmentForMutation} clearAlignmentForMutation
 * @property {() => void} rebindAnnotationRecords
 * @property {DiscoverDrawingRecords} discoverDrawingRecords
 */

/**
 * The one transition of a Linear drawing's record rows, Add File, and the
 * removal of a File.
 * @param {LinearSourceEditsPorts} ports
 */
export const createLinearSourceEdits = ({
  sessionOperationAvailability, shownAction, retireLegendStylesOfUnnamedCaptions, changeLinearDepthSources,
  linearDepthSourceLoss, changeRecordSources, askingStackPorts, pendingDiscoveryOf, replaceLinearComparisonPlan,
  invalidateLinearComparisonArtifacts, linearReorderNotice, retainAlignmentForStableReorder, clearAlignmentForMutation,
  rebindAnnotationRecords, discoverDrawingRecords
}) => {
  /** @param {DrawingState} drawing @param {readonly Partial<LinearSeq>[]} items */
  const applyLinearSeqMutation = (
    drawing, items,
    {
      preserveLosatCacheInfo = false,
      layoutEntries = drawing.linearRecordRows,
      alignmentMutation = msg('source set changed.'),
      retiredRecordKeys = /** @type {string[]} */ ([]),
      turnOffReaching = /** @type {(() => void) | undefined} */ (undefined)
    } = {}
  ) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const depthWidth = linearDepthLogicalWidth(drawing);
    const next = normalizeLinearSeqList(items);
    if (depthWidth > 0) {
      next.forEach((seq) => {
        seq.depth = padDepthFileSlots(seq.depth, depthWidth);
      });
    }
    // Removing a File can remove the last source of a Depth series.
    retireLegendStylesOfUnnamedCaptions(drawing, () => changeLinearDepthSources(drawing, () => {
      changeRecordSources(drawing, () => replaceLinearSeqs(drawing.slices.inputs, next));
      turnOffReaching?.();
    }));
    const activeUids = new Set(next.map((seq) => seq.uid));
    // A removed row of this drawing no longer waits; another drawing's rows are its own.
    const pending = pendingDiscoveryOf(drawing);
    pending.expansions.forEach((_, uid) => {
      if (!activeUids.has(uid)) pending.expansions.delete(uid);
    });
    pending.uploads.forEach((uid) => {
      if (!activeUids.has(uid)) pending.uploads.delete(uid);
    });
    const nextRows = reconcileLinearRecordLayout(drawing.linearSeqs, layoutEntries);
    drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...nextRows);
    // A card that stays but now reads another file (`retiredRecordKeys`) is ON.
    const retired = new Set(retiredRecordKeys);
    const recordsOff = pruneRecordsOff(drawing.recordsOff, new Set([...activeUids].filter((uid) => !retired.has(uid))));
    if (recordsOff.length !== drawing.recordsOff.length) drawing.recordsOff.splice(0, Infinity, ...recordsOff);
    replaceLinearComparisonPlan(
      drawing,
      reconcileLinearComparisonPlan(drawing.linearComparisonPlan, drawing.linearSeqs),
      { invalidate: false }
    );
    invalidateLinearComparisonArtifacts(drawing, { preserveLosatCacheInfo });
    linearReorderNotice.value = '';
    if (alignmentMutation === 'stable-reorder') {
      retainAlignmentForStableReorder(
        drawing.linearSeqs.map((sequence) => sequence.uid)
      );
    } else {
      clearAlignmentForMutation(alignmentMutation);
    }
    rebindAnnotationRecords();
    void discoverDrawingRecords(drawing);
  };

  const addLinearSeq = () => shownAction('linear', (drawing) => {
    return applyLinearSeqMutation(drawing, [...drawing.linearSeqs, createLinearSeq()]);
  });
  /**
   * Deletes a Linear File, or clears it to one blank File (`clear`), inside
   * the caller's History step: the removal dialog and Remove File (D-06).
   * @param {DrawingState} drawing
   * @param {Extract<LinearSourceRemovalPlan, { allowed: true }>} plan
   * @param {'clear' | 'delete'} intent
   * @param {() => void} [turnOffReaching]
   */
  const removeLinearSource = (drawing, plan, intent, turnOffReaching) => {
    const next = [...plan.retainedSequences];
    if (intent === 'clear') next.splice(plan.insertionIndex, 0, createLinearSeq());
    return applyLinearSeqMutation(drawing, next, { turnOffReaching });
  };
  /**
   * A Linear File removal asks, after its own question, when it takes the last
   * Depth source of a series that another row reaches (`askTrackSourceRemoval`).
   * @param {DrawingState} drawing The shown Linear drawing the File leaves.
   * @param {readonly Partial<LinearSeq>[]} kept The record cards with a Depth file that stay.
   * @param {string} label
   * @param {(turnOffReaching: () => void) => unknown} apply
   */
  const askLinearFileRemoval = (drawing, kept, label, apply) => askTrackSourceRemoval(askingStackPorts('linear_track_slots'), {
    kind: 'linear-file',
    loss: linearDepthSourceLoss(drawing, kept.map((seq) => seq.depth)),
    sources: () => linearDepthSources(drawing),
    label,
    apply
  });
  return { applyLinearSeqMutation, addLinearSeq, removeLinearSource, askLinearFileRemoval };
};

/**
 * @typedef {object} LinearSourceRemovalDialogPorts
 * @property {{ readonly value: LinearSourceGroup[] }} linearSourceGroups The Files of the shown Linear drawing.
 * @property {{ open: boolean, sourceUid: string, origin: string }} linearSourceRemovalDialog
 * @property {{ value: ReturnFocusTarget | null }} linearSourceRemovalReturnFocus
 * @property {SessionOperationAvailability} sessionOperationAvailability
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {RunUndoable} runUndoable
 * @property {RemoveLinearSource} removeLinearSource
 * @property {AskLinearFileRemoval} askLinearFileRemoval
 */

/**
 * The File removal dialog of the shown Linear drawing (Clear File, Delete
 * File, and the last File's removal) and where the focus goes after it.
 * @param {LinearSourceRemovalDialogPorts} ports
 */
export const createLinearSourceRemovalDialog = ({
  linearSourceGroups, linearSourceRemovalDialog, linearSourceRemovalReturnFocus, sessionOperationAvailability,
  shownDrawingOfMode, runUndoable, removeLinearSource, askLinearFileRemoval
}) => {
  const restoreLinearSourceRemovalFocus = async () => {
    await nextTick();
    const target = linearSourceRemovalReturnFocus.value;
    linearSourceRemovalReturnFocus.value = null;
    if (target?.isConnected && typeof target.focus === 'function') target.focus();
  };
  const closeLinearSourceRemovalDialog = ({ restoreFocus = true } = {}) => {
    linearSourceRemovalDialog.open = false;
    linearSourceRemovalDialog.sourceUid = '';
    linearSourceRemovalDialog.origin = '';
    if (restoreFocus) void restoreLinearSourceRemovalFocus();
    else linearSourceRemovalReturnFocus.value = null;
  };
  const focusLinearSourceRemovalDialog = async () => {
    await nextTick();
    /** @type {HTMLElement | null} */ (document.querySelector('[data-linear-source-removal-primary]'))?.focus();
  };
  // A folded card takes the focus on its header line and stays folded.
  /** @param {number} sourceIndex */
  const focusLinearSourceAfterRemoval = async (sourceIndex) => {
    await nextTick();
    const cards = document.querySelectorAll('[data-linear-source-card]');
    const targetIndex = Math.max(0, Math.min(Number(sourceIndex) || 0, cards.length - 1));
    const target = /** @type {HTMLElement | null} */ (cards[targetIndex]?.querySelector(
      '[data-linear-source-fold]:not([open]) > summary, .upload-zone[tabindex="0"]'
    ) || document.querySelector('[data-linear-file-add]'));
    if (typeof target?.focus === 'function') target.focus();
  };
  /** @param {LinearSourceGroup} source @param {string} origin @param {EventTarget | null} [returnFocus] */
  const openLinearSourceRemovalDialog = (source, origin, returnFocus = null) => {
    if (!source?.uid) return false;
    linearSourceRemovalDialog.sourceUid = source.uid;
    linearSourceRemovalDialog.origin = origin;
    linearSourceRemovalDialog.open = true;
    linearSourceRemovalReturnFocus.value = returnFocus;
    void focusLinearSourceRemovalDialog();
    return true;
  };
  const cancelLinearSourceRemoval = () => closeLinearSourceRemovalDialog();
  /** @param {'clear' | 'delete'} intent */
  const applyLinearSourceRemoval = async (intent) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const drawing = shownDrawingOfMode('linear');
    if (!drawing) return false;
    if (linearSourceRemovalDialog.origin === 'global' && intent !== 'delete') return false;
    const plan = planLinearSourceRemoval({
      sequences: drawing.linearSeqs,
      sourceUid: linearSourceRemovalDialog.sourceUid,
      intent
    });
    if (!plan.allowed) {
      closeLinearSourceRemovalDialog();
      return false;
    }
    const label = intent === 'clear' ? msg('Clear Linear File') : msg('Delete Linear File');
    const asked = askLinearFileRemoval(drawing, plan.retainedSequences, label, (turnOffReaching) => removeLinearSource(drawing, plan, intent, turnOffReaching));
    if (asked) {
      // The File's dialog gives way to the second question; its Cancel leaves the File and returns the focus.
      const returnFocus = linearSourceRemovalReturnFocus.value;
      closeLinearSourceRemovalDialog({ restoreFocus: false });
      const answer = await asked;
      if (answer === false) {
        linearSourceRemovalReturnFocus.value = returnFocus;
        void restoreLinearSourceRemovalFocus();
      } else await focusLinearSourceAfterRemoval(plan.sourceIndex);
      return answer;
    }
    const operation = await runUndoable(label, () => removeLinearSource(drawing, plan, intent));
    closeLinearSourceRemovalDialog({ restoreFocus: false });
    await focusLinearSourceAfterRemoval(plan.sourceIndex);
    return operation;
  };
  /** @param {LinearSourceGroup} source @param {EventTarget | null} [returnFocus] */
  const requestLinearSourceRemoval = (source, returnFocus = null) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!source || !linearSourceHasPrimaryInput(source)) return false;
    return openLinearSourceRemovalDialog(source, 'card', returnFocus);
  };
  /** @param {{ currentTarget?: EventTarget | null } | null} [event] */
  const removeLastLinearSeq = (event = null) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const source = linearSourceGroups.value.at(-1);
    if (!source || linearSourceGroups.value.length <= 1) return false;
    if (isPristineLinearSource(source)) {
      linearSourceRemovalDialog.sourceUid = source.uid;
      linearSourceRemovalDialog.origin = 'global';
      return applyLinearSourceRemoval('delete');
    }
    return openLinearSourceRemovalDialog(source, 'global', event?.currentTarget || null);
  };
  return { cancelLinearSourceRemoval, applyLinearSourceRemoval, requestLinearSourceRemoval, removeLastLinearSeq };
};

// CI-03: the per-record presentation a Linear source file carries. Replacing
// or removing the file resets it in the replacement's History step.
const RETIRED_LINEAR_RECORD_PRESENTATION = Object.freeze({ region_record_id: '', ...wholeRecordCardWindows(), definition: '', record_subtitle: '' });
// The record reference names a record of the old file, so any replacement retires it.
const RETIRED_LINEAR_RECORD_REFERENCE = Object.freeze({ region_record_id: '' });

/**
 * @typedef {object} LinearSourceFileReplacementPorts
 * @property {WebState} state
 * @property {ShownAction} shownAction
 * @property {(drawing: DrawingState) => { expansions: Map<string, boolean>, uploads: Set<string> }} pendingDiscoveryOf
 *   The record discoveries a drawing's rows wait for.
 * @property {ApplyLinearSeqMutation} applyLinearSeqMutation
 * @property {AskLinearFileRemoval} askLinearFileRemoval
 * @property {(replacement: { removed: boolean, previousRecordCount: number }) => boolean} sourceReplacementRetiresPresentation
 *   Whether a source file replacement retires the record crop and titles.
 */

/**
 * A Linear record row's source file: the replacement every writer of it
 * (upload, replace, Remove) goes through, and the shown drawing's control.
 * @param {LinearSourceFileReplacementPorts} ports
 */
export const createLinearSourceFileReplacement = ({
  state, shownAction, pendingDiscoveryOf, applyLinearSeqMutation, askLinearFileRemoval, sourceReplacementRetiresPresentation
}) => {
  /**
   * A Linear record row's source file (upload, replace, Remove), found by the
   * row's uid, with the drawing's reactions to its new source: the shown
   * drawing's File card and the replacement for every drawing write through it.
   * @param {DrawingState} drawing
   * @param {string} rowUid
   * @param {'gb' | 'gff' | 'fasta'} field
   * @param {InputFileCandidate | InputFileRef} value
   */
  const replaceLinearSourceFile = (drawing, rowUid, field, value) => {
    const seq = drawing.linearSeqs.find((entry) => entry.uid === rowUid);
    if (!seq) return;
    const nextValue = state.inputFiles.registerFile(value);
    if (seq[field] === nextValue) return;
    // `groupLinearSourceRecords` puts every record card of the drawing in a group, and `seq` is one of them.
    const group = /** @type {LinearSourceGroup} */ (groupLinearSourceRecords(drawing.linearSeqs).find((entry) => (
      entry.records.some(({ sequence }) => sequence.uid === seq.uid)
    )));
    const members = new Set(group.records.map(({ sequence }) => sequence.uid));
    const replacement = createLinearSeq({
      ...group.sequence,
      [field]: nextValue,
      ...(field === 'gb' && nextValue ? { file_definition: '', file_subtitle: '' } : {}),
      ...(!seq[field] ? {} : sourceReplacementRetiresPresentation({
        removed: !nextValue, previousRecordCount: group.records.length
      }) ? RETIRED_LINEAR_RECORD_PRESENTATION : RETIRED_LINEAR_RECORD_REFERENCE)
    });
    // The card stays when it gets a file, or, on Remove, while it still holds a file its input type reads: a file
    // left in another type's slot keeps no card.
    const keepSource = Boolean(nextValue) || inputSourceFiles(drawing.inputType.value, replacement).holdsFile;
    const nextSeqs = drawing.linearSeqs.flatMap((entry) => (
      entry.uid === group.uid ? (keepSource ? [replacement] : [])
        : members.has(entry.uid) ? [] : [entry]
    ));
    /** @param {() => void} [turnOffReaching] */
    const apply = (turnOffReaching) => applyLinearSeqMutation(drawing, nextSeqs, {
      alignmentMutation: msg('source replaced.'), retiredRecordKeys: [...members], turnOffReaching
    });
    const asked = keepSource ? null : askLinearFileRemoval(drawing, nextSeqs, msg('Delete Linear File'), apply);
    if (asked) return asked;
    apply();
    if (keepSource) {
      const pending = pendingDiscoveryOf(drawing);
      pending.expansions.set(replacement.uid, !seq[field]);
      pending.uploads.add(replacement.uid);
    }
  };
  /** @param {number} index @param {'gb' | 'gff' | 'fasta'} field @param {FileLike | null | undefined} value */
  const setLinearSeqPrimaryFile = (index, field, value) => shownAction('linear', (drawing) => {
    const idx = Number(index);
    if (!Number.isInteger(idx) || idx < 0 || idx >= drawing.linearSeqs.length) return;
    if (!['gb', 'gff', 'fasta'].includes(field)) return;
    replaceLinearSourceFile(drawing, drawing.linearSeqs[idx].uid, field, value);
  });
  return { replaceLinearSourceFile, setLinearSeqPrimaryFile };
};

/**
 * @typedef {object} LinearSourceMovesPorts
 * @property {WebState} state
 * @property {{ readonly value: LinearSourceGroup[] }} linearSourceGroups The Files of the shown Linear drawing.
 * @property {ShownAction} shownAction
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {RunUndoable} runUndoable
 * @property {ApplyLinearSeqMutation} applyLinearSeqMutation
 */

/**
 * Moving a File of the shown Linear drawing, and why it cannot move.
 * @param {LinearSourceMovesPorts} ports
 */
export const createLinearSourceMoves = ({
  state, linearSourceGroups, shownAction, shownDrawingOfMode, runUndoable, applyLinearSeqMutation
}) => {
  /** @param {DrawingState} drawing @param {number} sourceIndex @param {number} direction */
  const linearSourceMovePlan = (drawing, sourceIndex, direction) => planLinearSourceRowMove({
    sourceGroups: linearSourceGroups.value,
    entries: drawing.linearRecordRows,
    sourceIndex,
    direction
  });
  const linearSourceMoveBlockedReason = computed(() => {
    const drawing = state.activeDrawing();
    return (
      linearSourceMovePlan(drawing, 0, 1).reason === 'custom-layout'
        ? t('Files cannot be moved while Record Layout does not give each File its own consecutive rows. Change the rows in Record Layout under Advanced comparison and layout.')
        : ''
    );
  });
  /** @param {number} sourceIndex @param {number} direction */
  const canMoveLinearSource = (sourceIndex, direction) => {
    const drawing = shownDrawingOfMode('linear');
    return drawing !== null && linearSourceMovePlan(drawing, sourceIndex, direction).allowed;
  };

  /** @param {number} sourceIndex @param {number} direction */
  const moveLinearSource = (sourceIndex, direction) => shownAction('linear', (drawing) => {
    const plan = linearSourceMovePlan(drawing, sourceIndex, direction);
    if (!plan.allowed) return;
    const next = moveLinearSourceGroup(drawing.linearSeqs, sourceIndex, direction);
    return runUndoable(msg('Move File'), () => {
      applyLinearSeqMutation(drawing, next, {
        preserveLosatCacheInfo: true,
        layoutEntries: plan.rows,
        alignmentMutation: 'stable-reorder'
      });
    });
  });
  return { linearSourceMoveBlockedReason, canMoveLinearSource, moveLinearSource };
};

/**
 * @typedef {object} LinearRecordEditsPorts
 * @property {SessionOperationAvailability} sessionOperationAvailability
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {RunUndoable} runUndoable
 * @property {ClearAlignmentForMutation} clearAlignmentForMutation
 * @property {(drawing: DrawingState, sequence: LinearSeq) => SequenceRecord[]} linearRecordsFor
 *   A card's records, read in the drawing that holds it.
 * @property {(drawing: DrawingState, sequence: LinearSeq | undefined) => RecordSelectorStatus} linearRecordStatusFor
 * @property {() => void} rebindAnnotationRecords
 * @property {ApplyLinearSeqMutation} applyLinearSeqMutation
 * @property {() => ReadonlyMap<string, LinearUnitPlace>} shownUnitPlaces The drawn units of the shown Linear Result.
 */

/**
 * A Linear record card's record selector, its windows (the Region block's crop, or the Windows list; PD-OI-124),
 * and the resets of its definition and subtitle.
 * @param {LinearRecordEditsPorts} ports
 */
export const createLinearRecordEdits = ({
  sessionOperationAvailability, shownDrawingOfMode, runUndoable, clearAlignmentForMutation, linearRecordsFor,
  linearRecordStatusFor, rebindAnnotationRecords, applyLinearSeqMutation, shownUnitPlaces
}) => {
  /** @param {LinearSeq} sequence @param {unknown} value */
  const setLinearRecordSelector = (sequence, value) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!sequence) return false;
    const next = String(value || '');
    if (String(sequence.region_record_id || '') === next) return false;
    clearAlignmentForMutation(msg('record selector changed.'));
    sequence.region_record_id = next;
    // The inferred definition follows the record the row now selects (D-12);
    // the row is a card of the shown Linear drawing.
    const drawing = shownDrawingOfMode('linear');
    if (drawing !== null && drawing.inputType.value === 'gb' && linearRecordStatusFor(drawing, sequence) === 'ready') {
      sequence.inferred_definition = inferredDefinitionForRecord(linearRecordsFor(drawing, sequence), next);
    }
    rebindAnnotationRecords();
    return true;
  };

  /**
   * A field of the card's window `index` (the Region block is window 0): Start or End (blank is empty), Reverse
   * complement, or Across the origin. One History step: the control's own. A change of the span clears the alignment;
   * Reverse complement alone is a direction the alignment keeps (`cardWindowEditKind`).
   * @param {LinearSeq} sequence @param {number} index
   * @param {'start' | 'end' | 'reverseComplement' | 'acrossOrigin'} field @param {unknown} value
   */
  const setLinearCardWindow = (sequence, index, field, value) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!sequence) return false;
    const next = field === 'reverseComplement' || field === 'acrossOrigin' ? Boolean(value)
      : value === '' || value === null || value === undefined ? null : Number(value);
    const edit = { [field]: next };
    const kind = cardWindowEditKind(sequence, index, edit);
    if (!kind) return false;
    if (kind === 'crop') clearAlignmentForMutation(msg('record crop changed.'));
    editCardWindow(sequence, index, edit);
    return true;
  };

  /**
   * Replaces a card's windows (Add window, Remove): one History step through the record rows' transition, which also
   * drops the stored comparisons whose unit indexes the change shifts (the next Generate plans them again).
   * @param {string} label @param {LinearSeq} sequence @param {(card: LinearSeq, drawing: DrawingState) => boolean} change
   */
  const changeLinearCardWindows = (label, sequence, change) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const drawing = shownDrawingOfMode('linear');
    if (!drawing || !sequence) return false;
    return runUndoable(label, () => {
      const card = { ...sequence, windows: sequence.windows.map((window) => ({ ...window })) };
      if (!change(card, drawing)) return;
      applyLinearSeqMutation(drawing, drawing.linearSeqs.map((seq) => (seq.uid === card.uid ? card : seq)), {
        preserveLosatCacheInfo: true, alignmentMutation: msg('record crop changed.')
      });
    });
  };
  /** "+ Add window": a window after the card's last one, in its row's next place. @param {LinearSeq} sequence */
  const addLinearCardWindow = (sequence) => changeLinearCardWindows(msg('Add window'), sequence, (card, drawing) => {
    const placed = drawing.linearRecordLayoutEnabled.value
      ? placeNewLayoutWindow(drawing.linearSeqs, drawing.linearRecordRows, card.uid) : null;
    const single = card.windows.length === 0;
    // The placement wrote the places of the card's own windows on the drawing's card.
    if (placed && !single) card.windows = sequence.windows.map((window) => ({ ...window }));
    addCardWindow(card, placed ? placed.column : null);
    if (placed && single) card.windows[0].column = placed.lastColumn;
    if (placed) drawing.linearRecordRows.splice(0, drawing.linearRecordRows.length, ...placed.entries);
    return true;
  });
  /**
   * Remove of window `index`; the last window stays. A name typed on a window pair stays with its window or goes with
   * it, in the same History step (OV-2742).
   * @param {LinearSeq} sequence @param {number} index
   */
  const removeLinearCardWindow = (sequence, index) => changeLinearCardWindows(
    msg('Remove window'), sequence, (card, drawing) => {
      if (!removeCardWindow(card, index)) return false;
      const recordKey = linearSeqRecordKey(card, drawing.linearSeqs.indexOf(sequence));
      const places = shownUnitPlaces();
      /** @type {LinearComparisonPlan} */ (drawing.linearComparisonPlan).edges.forEach((edge) => setWindowLosatFilenames(edge, windowLosatFilenamesAfterRemoval(
        windowLosatFilenamesOf(edge), places, recordKey, index + 1, cardWindowList(card).length)));
      return true;
    }
  );

  /** @param {LinearSeq} seq */
  const resetLinearRecordDefinition = (seq) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!seq) return;
    runUndoable(msg('Reset record definition'), () => {
      seq.definition = '';
    });
  };

  /** @param {LinearSeq} seq */
  const resetLinearRecordSubtitle = (seq) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!seq) return;
    runUndoable(msg('Reset record subtitle'), () => {
      seq.record_subtitle = '';
    });
  };
  return {
    setLinearRecordSelector, setLinearCardWindow, addLinearCardWindow, removeLinearCardWindow, resetLinearRecordDefinition,
    resetLinearRecordSubtitle
  };
};
