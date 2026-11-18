// @ts-check
// The Circular drawing's source file and record controls: the Multi-Record
// Canvas checkbox, the source file upload, replace and Remove with the record
// presentation they retire, the record presentation selector, the links that
// reveal the Multi-Record Canvas setting and the Ticks row, and the scroll
// anchoring when the record presentation panel appears. The composition root
// passes the state, the shown-drawing readers and admission, the record
// display's source change, the record discovery, the rebind of annotations,
// the settings pane's reveal, the panels' refs, and the record form defaults
// as ports.
import {
  addCircularSource, bindCircularSourceFile, circularSlotFileKind, createCircularSource, markCircularDepthRows, moveCircularSource,
  removeCircularSource, sessionFilesOfKind
} from '../services/drawing/inputs.js';
import { isInputFileRef } from '../services/history-files.js';
import { isRecordMajorDepthFileMatrix } from '../services/depth-track-state.js';
import { circularSourceInputs, inputSourceFiles, inputTypeSpec, pickedFileCards } from '../services/input-types.js';
import { circularMarkRecords, settledCircularCardMarks } from '../services/record-discovery.js';
import { cardHoldsRecordKey, cardRowsLeaving } from '../services/record-draw-selection.js';
import { cardCropWindow } from '../services/drawing/inputs.js';
import { circularRangesWithoutCard, circularRecordEntries, circularRecordRange } from '../services/record-options.js';
import { msg, t } from '../utils/i18n.js';
import { circularMultiRecordDefaults, updateActiveLayoutPreference } from '../services/layout-preferences.js';
import { settingsScrollOwner } from '../services/settings-sections.js';
import { chooseRecordLayout } from '../services/one-circle-layout.js';
import { drawingTabLabels } from './input-replacement-ports.js';
/** @import { DiagramMode, DrawingState, SessionBusy, WebState } from '../state.js' */
/** @import { InputFileCandidate, InputFileRef } from '../services/history-files.js' */
/** @import { CircularTrackSlot } from '../services/circular-track-slot-model.js' */
/** @import { LinearSeqWindow } from '../services/drawing/inputs.js' */
/** @import { CircularListRecord } from '../services/record-discovery.js' */
/** @import { FileLike } from '../services/session-resource-backing.js' */
/** @import { AskChoice, ShowNotice } from './app-messages.js' */
/** @import { DrawingBarOffer, DrawingBarTab, DrawingInputSite } from './drawing-list.js' */
/** @import { RevealSetting } from './settings-pane.js' */
/** @import { CircularRecordLayoutChoice } from '../services/one-circle-layout.js' */

const { computed, nextTick, watch } = window.Vue;

/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */

/**
 * @typedef {object} CircularMultiRecordCanvasPorts
 * @property {WebState} state
 * @property {() => void} rebindAnnotationRecords Rebinds the shown drawing's annotations to its record catalog (R10).
 * @property {(drawing: DrawingState, hidden: boolean) => void} projectCircularSkewHidden Writes the drawing's
 *   Hide GC Skew and projects it onto its track rows, without the Hide GC Skew checkbox's question: a layout's
 *   own value hides no row for good (its marked rows come back when the layout is left).
 */

/**
 * Registers nothing; the root calls it after its `rebindAnnotationRecords`.
 * @param {CircularMultiRecordCanvasPorts} ports
 */
export const createCircularMultiRecordCanvas = ({ state, rebindAnnotationRecords, projectCircularSkewHidden }) => {
  // A Circular drawing whose Multi-Record Canvas is on and has no multi-record
  // sides stored gets the toggle's smart defaults. Load, its rollback, and a
  // History restore install the sides they read or captured (R11, D-15-48).
  /** @param {DrawingState} drawing */
  const storeCircularMultiRecordDefaults = (drawing) => {
    const defaults = drawing.mode === 'circular' && drawing.form.multi_record_canvas
      ? circularMultiRecordDefaults(drawing.layoutPreferences) : null;
    if (defaults) updateActiveLayoutPreference(drawing.layoutPreferences, 'circular', true, defaults);
  };
  // After the record layout choice wrote the shown drawing: in the control's History step.
  const changeCircularMultiRecordCanvas = () => {
    storeCircularMultiRecordDefaults(state.activeDrawing());
    rebindAnnotationRecords();
  };
  // The record layout choice (PD-OI-019, `chooseRecordLayout`): separate figures, the Multi-Record Canvas grid,
  // one circle, or (`'left'`) the layout one circle left, with its GC values; one control, one History step.
  /** @param {CircularRecordLayoutChoice | 'left'} choice */
  const chooseCircularRecordLayout = (choice) => {
    const drawing = state.activeDrawing();
    const skewHidden = chooseRecordLayout(drawing, choice);
    if (skewHidden !== null && skewHidden !== Boolean(drawing.form.suppress_skew)) {
      projectCircularSkewHidden(drawing, skewHidden);
    }
    changeCircularMultiRecordCanvas();
  };
  return { chooseCircularRecordLayout };
};

/**
 * A mode panel's mutating handler: the session's busy answer, nothing when the
 * shown drawing is of another mode, else `action` on the shown drawing.
 * @typedef {<T>(drawingMode: DiagramMode, action: (drawing: DrawingState) => T) => T | SessionBusy | undefined} ShownAction
 */
/**
 * The record form fields of a Circular drawing as a new drawing has them.
 * @typedef {Readonly<{ circular_record_selector: string, circular_region_start: number | null,
 *   circular_region_end: number | null, circular_reverse: boolean, circular_windows: readonly LinearSeqWindow[],
 *   circular_record_label: string, circular_record_subtitle: string }>} CircularRecordFormDefaults
 */

/**
 * @typedef {object} CircularSourceEditsPorts
 * @property {WebState} state
 * @property {CircularRecordFormDefaults} CIRCULAR_RECORD_FORM_DEFAULTS
 * @property {<T>(drawing: DrawingState, mutate: () => T) => T} changeRecordSources
 *   Runs a write of the drawing's record sources; the replaced source's display drafts go.
 * @property {(drawing: DrawingState) => { circularUpload: InputFileRef | null }} pendingDiscoveryOf
 *   What the drawing's next record discovery does with a Circular upload.
 * @property {(drawing: DrawingState) => unknown} discoverDrawingRecords Reads the records of the drawing's inputs.
 * @property {ShownAction} shownAction
 * @property {(label: string, fn: () => unknown) => Promise<unknown>} runUndoable History's undoable step.
 */

/**
 * Before a card edit, in its History step (OV-2868): the rows of a
 * record-major Depth matrix stand for the cards as last read, so the read
 * after the edit moves them with their records (`circularDepthRowsAfterRead`
 * in the record reader). Marks a Load or an earlier edit whose read has not
 * settled stay; a series list or no Depth needs none.
 * @param {DrawingState} drawing
 */
const markDepthRowsBeforeEdit = (drawing) => {
  if (!isRecordMajorDepthFileMatrix(drawing.files.c_depth) || drawing.files.c_depth_rows) return;
  markCircularDepthRows(drawing.slices.inputs, settledCircularCardMarks(drawing));
};

/**
 * The moves of a Circular drawing's source cards: one place up or down, one
 * History step, the Depth matrix following its records once they are read.
 * @param {Pick<CircularSourceEditsPorts, 'state' | 'discoverDrawingRecords' | 'shownAction' | 'runUndoable'>} ports
 */
const createCircularSourceMoves = ({ state, discoverDrawingRecords, shownAction, runUndoable }) => {
  /** @param {DrawingState} drawing @param {number} index @param {number} direction */
  const movedIndex = (drawing, index, direction) => {
    const target = index + Math.sign(direction);
    return Number.isInteger(index) && index >= 0 && index < drawing.circularSources.length
      && target >= 0 && target < drawing.circularSources.length ? target : null;
  };
  /** @param {number} index @param {number} direction */
  const canMoveCircularSource = (index, direction) => {
    const drawing = state.activeDrawing();
    return drawing.mode === 'circular' && movedIndex(drawing, index, direction) !== null;
  };
  /**
   * Moves a source card one place up (`direction` -1) or down (1), one History step.
   * @param {number} index @param {number} direction
   */
  const moveCircularSourceCard = (index, direction) => shownAction('circular', (drawing) => {
    const target = movedIndex(drawing, index, direction);
    if (target === null) return undefined;
    const { uid } = drawing.circularSources[index];
    return runUndoable(msg('Move File'), () => {
      markDepthRowsBeforeEdit(drawing);
      moveCircularSource(drawing.slices.inputs, uid, target);
      void discoverDrawingRecords(drawing);
    });
  });
  return { canMoveCircularSource, moveCircularSourceCard };
};

/**
 * The slot of the input type that a file of this Session fills: the first
 * slot whose kind the Session binds it as (a ref chosen from "Files in this
 * Session" keeps its kind, whatever its name says); null for a file from disk.
 * @param {Parameters<typeof sessionFilesOfKind>[0]} bindings
 * @param {unknown} inputType
 * @returns {(file: InputFileCandidate | InputFileRef) => 'gb' | 'gff' | 'fasta' | null}
 */
const sessionSlotOf = (bindings, inputType) => (file) => {
  if (!isInputFileRef(file)) return null;
  const { primary, paired } = inputTypeSpec(inputType);
  for (const slot of paired ? [primary, paired] : [primary]) {
    if (sessionFilesOfKind(bindings, circularSlotFileKind(slot)).some((ref) => ref.fileId === file.fileId)) return slot;
  }
  return null;
};

// D-32: one record replacing one record (a new version of the same genome)
// keeps its crop and titles; a removed file or a replaced multi-record file
// retires them. The record selector names a record of the old file, so any
// replacement retires it.
/** @param {{ removed: boolean, previousRecordCount: number }} replacement */
const sourceReplacementRetiresPresentation = ({ removed, previousRecordCount }) => removed || previousRecordCount > 1;

/**
 * Keeps the items of a reactive list that `keep` accepts, in place; a list that loses none is not written.
 * @template T @param {T[]} list @param {(item: T) => boolean} keep
 */
const keepInPlace = (list, keep) => {
  const kept = list.filter(keep);
  if (kept.length !== list.length) list.splice(0, Infinity, ...kept);
};
/** @typedef {{ uid: string, field: 'gb' | 'gff' | 'fasta' }} ReplacedCardFile One file of a source card a Replace writes. */
/** @typedef {InputFileRef | NonNullable<InputFileCandidate>} ReplacementFile A card's file: bound, or the one a Replace writes. */
/** @typedef {import('../services/record-discovery.js').CardRecordIds<ReplacementFile>} ReplacementCardRecords */
/**
 * What a Replace was picked on, compared by identity: the drawing's input type and its cards, in order, with their
 * files. A change ends the Replace (a later pick of the card, the card removed or moved, an Undo).
 * @param {DrawingState} drawing @returns {unknown[]}
 */
const replacedCards = (drawing) => [
  drawing.inputType.value, ...drawing.circularSources.flatMap((card) => [card.uid, card.gb, card.gff, card.fasta])
];
/**
 * The Depth files a drawing's Replace question names, compared by identity: the matrix and the cards its rows stand
 * for. A change before the dialog opens reads them again; a change while it is open makes the answer do nothing.
 * @param {DrawingState} drawing @returns {unknown[]}
 */
const replacementDepthFiles = (drawing) => [
  drawing.files.c_depth_rows,
  ...(isRecordMajorDepthFileMatrix(drawing.files.c_depth)
    ? /** @type {(InputFileRef | null)[][]} */ (drawing.files.c_depth).flatMap((row) => [row.length, ...row]) : [drawing.files.c_depth])
];
/**
 * The rows of the drawing's record-major Depth matrix when one of them holds a file; else null.
 * @param {DrawingState} drawing
 */
const depthRowsWithFiles = (drawing) => {
  const matrix = drawing.files.c_depth;
  return isRecordMajorDepthFileMatrix(matrix) && /** @type {(InputFileRef | null)[][]} */ (matrix).some((row) => row.some(Boolean))
    ? /** @type {(InputFileRef | null)[][]} */ (matrix) : null;
};
/** @param {readonly unknown[]} asked @param {readonly unknown[]} now */
const sameEntries = (asked, now) => now.length === asked.length && now.every((entry, index) => entry === asked[index]);

/**
 * The Depth files of `drawing` that leave when `value` replaces the card files `replaced` (OV-2318), each with the
 * record it stands for. The rows stand for the cards of `c_depth_rows` (else the cards as last read, which the edit
 * marks), and the read after the edit keeps the row of each record with the record it continues as
 * (`circularDepthRowsAfterRead` with `resliceCardRows`): both sets of cards are read by `circularMarkRecords`, as that
 * read reads them, and the rows of the records the new file has no record for leave (`cardRowsLeaving`: by record ID,
 * else by place, D-05). None when the rows cannot be placed (a file that cannot be read, rows that are not those of
 * their cards): that read reports the matrix, as before.
 * @param {DrawingState} drawing
 * @param {readonly ReplacedCardFile[]} replaced
 * @param {ReplacementFile} value
 * @param {(file: ReplacementFile | null) => FileLike | null | undefined} fileOf
 * @returns {Promise<{ k: number, record: string | null, files: string[] }[]>}
 */
const circularDepthFilesLeaving = async (drawing, replaced, value, fileOf) => {
  const settledMarks = settledCircularCardMarks(drawing);
  /** @type {readonly ReplacementCardRecords[]} */
  const settled = settledMarks ?? [];
  /** @type {readonly ReplacementCardRecords[] | null} */
  const marks = drawing.files.c_depth_rows ?? settledMarks;
  const rows = depthRowsWithFiles(drawing);
  if (!marks || !rows) return [];
  const inputType = drawing.inputType.value;
  /** @type {import('../services/input-types.js').CircularSourceInput<ReplacementFile>[]} */
  const inputs = circularSourceInputs(inputType, drawing.circularSources);
  const marked = await circularMarkRecords({ marks, now: settled, inputs, fileOf });
  if (!marked || marked.reduce((sum, mark) => sum + mark.recordIds.length, 0) !== rows.length) return [];
  /** @type {import('../services/input-types.js').CircularSourceInput<ReplacementFile>[]} */
  const after = circularSourceInputs(inputType, drawing.circularSources.map((card) => replaced.reduce(
    /** @param {Readonly<Partial<Record<'gb' | 'gff' | 'fasta', ReplacementFile | null>>> & { uid: string }} next */
    (next, { uid, field }) => (uid === card.uid ? { ...next, [field]: value } : next), card
  )));
  // The cards as marked are read already; a replaced card's new files are read.
  const known = [...settled, ...marks.map((mark, index) => ({ ...mark, recordIds: marked[index].recordIds }))];
  const now = await circularMarkRecords({
    marks: after.map(({ uid, primaryFile, pairedFile }) => ({ uid, primaryFile, pairedFile, recordIds: null })), now: known, inputs: after, fileOf
  });
  if (!now) return [];
  return cardRowsLeaving(rows, marked, now).flatMap(({ uid, k, row }) => {
    const files = row.flatMap((ref) => (ref ? [ref.name] : []));
    const record = drawing.circularRecordList.value.find((entry) => entry.uid === uid && entry.k === k);
    return files.length === 0 ? [] : [{ k, record: record ? record.record_id : null, files }];
  });
};

/**
 * The question a card's file control asks before a Replace that removes Depth files (OV-2318). Registers nothing.
 * @param {{ state: WebState, askChoice: AskChoice, showNotice: ShowNotice, barTabs: () => Array<DrawingBarTab | DrawingBarOffer> }} ports
 *   `barTabs`: the drawing bar's tabs, named as they show.
 * @returns {(sites: readonly DrawingInputSite[], value: ReplacementFile, pickHolds: () => boolean) => Promise<boolean> | null}
 */
export const createCircularSourceReplacementQuestion = ({ state, askChoice, showNotice, barTabs }) => {
  /** @param {ReplacementFile | null} file */
  const fileOf = (file) => /** @type {FileLike | null} */ (isInputFileRef(file) ? state.inputFiles.fileOf(file) : file);
  /**
   * Before a Replace writes `value` at the source card `sites` (the shown
   * card, or every card of every drawing the chosen scope replaces the file
   * in): when Depth files of some drawing would leave
   * (`circularDepthFilesLeaving`), one dialog names them, with the drawing
   * when the sites are in more than one (Replace file / Cancel). The files are
   * read first, with no dialog open: a change of the replaced drawings' cards
   * (a later pick, the card removed or moved, an Undo or Redo of the card, a
   * replaced drawing deleted), or a pick that no longer holds (`pickHolds`:
   * its drawing not shown once the file is read), ends the Replace (false).
   * The cards are compared first: their change is the user's own later act
   * and ends it silently. Otherwise another drawing shown does not supersede
   * the pick, so a notice says the file was not put on the card, naming its
   * drawing (D-09). A History step of another setting changes nothing the
   * pick holds. A change of their Depth files reads
   * them again, so the dialog names the files as they are. An answer given
   * after either change does nothing (as `askTrackSourceRemoval`). Null when
   * no Depth file can leave (the control's write runs in its own step); else
   * whether to replace, after which the write needs a step of its own.
   * @param {readonly DrawingInputSite[]} sites
   * @param {ReplacementFile} value
   * @param {() => boolean} pickHolds
   * @returns {Promise<boolean> | null}
   */
  const askCircularSourceReplacement = (sites, value, pickHolds) => {
    /** @type {Map<DrawingState['id'], ReplacedCardFile[]>} */
    const replacedIn = new Map();
    for (const { drawingId, slot, rowUid } of sites) {
      const field = slot === 'source.gb' ? 'gb' : slot === 'source.gff' ? 'gff' : slot === 'source.fasta' ? 'fasta' : null;
      if (field === null || rowUid === null) continue;
      replacedIn.set(drawingId, [...(replacedIn.get(drawingId) ?? []), { uid: rowUid, field }]);
    }
    const inScope = [...replacedIn].map(([drawingId, replaced]) => ({ drawingId, drawing: state.drawingOf(drawingId), replaced }));
    // The drawings where the Replace changes a card holding Depth files.
    const asking = () => inScope.filter(({ drawing, replaced }) => depthRowsWithFiles(drawing)
      && replaced.some(({ uid, field }) => drawing.circularSources.some((card) => card.uid === uid && card[field] !== value)));
    if (asking().length === 0) return null;
    const cardsNow = () => inScope.flatMap(({ drawing }) => replacedCards(drawing));
    const depthNow = () => inScope.flatMap(({ drawing }) => replacementDepthFiles(drawing));
    const picked = cardsNow();
    // The replaced drawings, by identity, and their cards as at the pick.
    const cardsHold = () => inScope.every(({ drawing }) => state.drawingList().includes(drawing)) && sameEntries(picked, cardsNow());
    const [shown] = sites;
    const cardNumber = state.drawingOf(shown.drawingId).circularSources.findIndex((card) => card.uid === shown.rowUid) + 1;
    const stillShown = () => {
      if (pickHolds()) return true;
      showNotice({
        text: () => t('{file} was not put on File {n} of the drawing “{drawing}”: another drawing was shown while it was read. Pick it again.', {
          file: value.name ?? '', n: cardNumber, drawing: drawingTabLabels(barTabs(), [shown.drawingId]).join('')
        }),
        tone: 'warning'
      });
      return false;
    };
    return (async () => {
      /** @type {unknown[]} */
      let asked;
      /** @type {{ drawingId: DrawingState['id'], rows: Awaited<ReturnType<typeof circularDepthFilesLeaving>> }[]} */
      let leaving;
      do {
        asked = depthNow();
        leaving = await Promise.all(asking().map(async ({ drawingId, drawing, replaced }) => ({
          drawingId, rows: await circularDepthFilesLeaving(drawing, replaced, value, fileOf)
        })));
        if (!cardsHold() || !stillShown()) return false;
      } while (!sameEntries(asked, depthNow()));
      const named = inScope.length > 1;
      const items = leaving.flatMap(({ drawingId, rows }) => rows.map(({ k, record, files }) => {
        const texts = { n: k, files: files.join(', '), drawing: drawingTabLabels(barTabs(), [drawingId]).join('') };
        if (named) return record ? t('#{n} {record} in {drawing}: {files}', { ...texts, record }) : t('#{n} in {drawing}: {files}', texts);
        return record ? t('#{n} {record}: {files}', { ...texts, record }) : t('#{n}: {files}', texts);
      }));
      if (items.length === 0) return true;
      const answer = await askChoice({
        title: t('Replace the file and remove Depth files?'),
        description: t('The new file has no record for these records of the File card. Replacing the file removes their Depth files.'),
        items,
        actions: [{ id: 'replace', label: t('Replace file'), tone: 'danger' }]
      });
      return answer === 'replace' && cardsHold() && stillShown() && sameEntries(asked, depthNow());
    })();
  };
  return askCircularSourceReplacement;
};

/**
 * The edits of a Circular drawing's source cards (circular-multi-file D-01):
 * add cards with their files, replace or remove one file of a card, remove a
 * card, and move a card, each one History step (R11), as the Linear File edits
 * (app/linear-source-edits.js). A card's records keep their identities
 * (`<uid>:<k>`) through every edit of the other cards, so their rows, OFF
 * state, Record select and Record Order stay; the edits retire only what names
 * the records of the card they replace or remove. A record-major Depth matrix
 * follows its records once the cards are read again (`markDepthRowsBeforeEdit`). Registers nothing.
 * @param {CircularSourceEditsPorts} ports
 */
export const createCircularSourceEdits = ({
  state, CIRCULAR_RECORD_FORM_DEFAULTS, changeRecordSources, pendingDiscoveryOf, discoverDrawingRecords, shownAction,
  runUndoable
}) => {
  // CI-03: the per-record presentation a source file carries in each mode.
  // Replacing or removing the file resets it in the replacement's History step.
  /** @param {DrawingState} drawing */
  const retireCircularRecordSelector = ({ form, adv }) => {
    form.circular_record_selector = CIRCULAR_RECORD_FORM_DEFAULTS.circular_record_selector;
    if (adv.circular_grouping_intent === 'single') adv.circular_grouping_intent = 'auto';
  };
  // The whole record presentation: the record form fields return to their
  // defaults (the one list a style copy leaves out too). Each field is written
  // by name, so the track-layout writer guard reads no whole-form write; a
  // node case pins that every field of the list is written.
  /** The single-record crop fields. @param {DrawingState} drawing */
  const retireCircularCropFields = ({ form }) => {
    form.circular_region_start = CIRCULAR_RECORD_FORM_DEFAULTS.circular_region_start;
    form.circular_region_end = CIRCULAR_RECORD_FORM_DEFAULTS.circular_region_end;
    form.circular_reverse = CIRCULAR_RECORD_FORM_DEFAULTS.circular_reverse;
    form.circular_windows = [...CIRCULAR_RECORD_FORM_DEFAULTS.circular_windows];
  };
  /** @param {DrawingState} drawing */
  const retireCircularRecordPresentation = (drawing) => {
    const { form } = drawing;
    retireCircularRecordSelector(drawing);
    retireCircularCropFields(drawing);
    form.circular_record_label = CIRCULAR_RECORD_FORM_DEFAULTS.circular_record_label;
    form.circular_record_subtitle = CIRCULAR_RECORD_FORM_DEFAULTS.circular_record_subtitle;
  };
  /**
   * The card's records leave the drawing's references: its OFF records and the
   * records it names of its file (`recordsNamed`, one value with them), and the
   * Record select (with the presentation when `retirePresentation`) when it
   * names one of them.
   * @param {DrawingState} drawing @param {string} uid @param {boolean} retirePresentation
   */
  const retireCardRecords = (drawing, uid, retirePresentation) => {
    keepInPlace(drawing.recordsOff, (key) => !cardHoldsRecordKey(uid, key));
    keepInPlace(drawing.recordsNamed, (card) => card.uid !== uid);
    // The grid or batch crops of the card's records go with their presentation (PD-OI-044 scenario revision 5).
    const ranges = circularRangesWithoutCard(drawing.form.circular_record_ranges, uid);
    if (retirePresentation && ranges.length !== drawing.form.circular_record_ranges.length) {
      drawing.form.circular_record_ranges = ranges;
    }
    // The crop fields filled from one of those crops go with it (OV-3140), so no other record draws it.
    if (retirePresentation && cardHoldsRecordKey(uid, drawing.form.circular_range_fields_record)) {
      retireCircularCropFields(drawing);
      drawing.form.circular_range_fields_record = '';
    }
    if (!cardHoldsRecordKey(uid, String(drawing.form.circular_record_selector || '').trim())) return;
    if (retirePresentation) retireCircularRecordPresentation(drawing);
    else retireCircularRecordSelector(drawing);
  };
  /** @param {DrawingState} drawing @param {InputFileRef | null} ref */
  const readUpload = (drawing, ref) => {
    if (ref) pendingDiscoveryOf(drawing).circularUpload = ref;
    void discoverDrawingRecords(drawing);
  };
  /**
   * One file of a Circular source card (upload, replace, Remove), found by the
   * card's uid, with the drawing's reactions to its new source: the shown
   * drawing's card and the replacement for every drawing that uses a file
   * write through it (`INPUT_REPLACEMENT_WRITERS`). A card left with no file
   * its input type reads is removed.
   * @param {DrawingState} drawing
   * @param {string} uid
   * @param {'gb' | 'gff' | 'fasta'} field
   * @param {InputFileCandidate | InputFileRef | null} value
   */
  const replaceCircularSourceFile = (drawing, uid, field, value) => {
    const card = drawing.circularSources.find((source) => source.uid === uid);
    if (!card) return;
    const previous = card[field];
    const nextRef = state.inputFiles.registerFile(value ?? null);
    if (previous === nextRef) return;
    const keepCard = Boolean(nextRef) || inputSourceFiles(drawing.inputType.value, { ...card, [field]: null }).holdsFile;
    const previousRecordCount = drawing.circularRecordList.value.filter((/** @type {CircularListRecord} */ record) => record.uid === uid).length;
    markDepthRowsBeforeEdit(drawing);
    changeRecordSources(drawing, () => {
      if (keepCard) bindCircularSourceFile(state.inputFiles, drawing.slices.inputs, uid, field, nextRef);
      else removeCircularSource(drawing.slices.inputs, uid);
    });
    readUpload(drawing, nextRef);
    // An upload into an empty slot starts reading another input (after an
    // Input switch), whose records the OFF list does not name.
    retireCardRecords(drawing, uid, !previous ? false : sourceReplacementRetiresPresentation({
      removed: !nextRef, previousRecordCount
    }));
  };
  /**
   * New source cards at the end of the drawing's list from the files of one
   * pick or drop, in their order, as the drawing's input type makes cards of
   * them (`pickedFileCards`: one per file; GFF3 + FASTA pairs by file-name
   * stem; a file of this Session goes to the slot of the kind it is bound as,
   * `sessionSlotOf`): one History step for every file.
   * @param {readonly (InputFileCandidate | InputFileRef)[]} files
   */
  const addCircularSourceFiles = (files) => shownAction('circular', (drawing) => (files.length === 0 ? undefined
    : runUndoable(msg('Add File'), () => {
      /** @type {InputFileRef | null} */
      let first = null;
      markDepthRowsBeforeEdit(drawing);
      changeRecordSources(drawing, () => pickedFileCards(drawing.inputType.value, files, sessionSlotOf(state, drawing.inputType.value)).forEach((slots) => {
        const card = createCircularSource();
        addCircularSource(drawing.slices.inputs, card);
        for (const [field, file] of /** @type {['gb' | 'gff' | 'fasta', InputFileCandidate | InputFileRef][]} */ (Object.entries(slots))) {
          const ref = bindCircularSourceFile(state.inputFiles, drawing.slices.inputs, card.uid, field, file);
          first ??= ref;
        }
      }));
      readUpload(drawing, first);
    })));
  /**
   * Removes a source card and every reference to its records, inside the
   * caller's History step (the record list's Remove File dialog).
   * @param {DrawingState} drawing @param {string} uid
   * @returns {boolean} whether the drawing held the card
   */
  const removeCircularSourceOf = (drawing, uid) => {
    if (!drawing.circularSources.some((source) => source.uid === uid)) return false;
    markDepthRowsBeforeEdit(drawing);
    changeRecordSources(drawing, () => removeCircularSource(drawing.slices.inputs, uid));
    retireCardRecords(drawing, uid, true);
    void discoverDrawingRecords(drawing);
    return true;
  };
  /**
   * Removes a source card of the shown drawing, one History step.
   * @param {string} uid
   */
  const removeCircularSourceCard = (uid) => shownAction('circular', (drawing) => (
    drawing.circularSources.some((source) => source.uid === uid)
      ? runUndoable(msg('Remove File'), () => removeCircularSourceOf(drawing, uid)) : undefined));
  /**
   * The file controls of one card of the shown drawing.
   * @param {string} uid
   * @param {'gb' | 'gff' | 'fasta'} field
   * @param {File | null} value
   */
  const setCircularSourceFile = (uid, field, value) => shownAction('circular', (drawing) => replaceCircularSourceFile(drawing, uid, field, value));
  return {
    sourceReplacementRetiresPresentation, replaceCircularSourceFile, setCircularSourceFile, addCircularSourceFiles,
    removeCircularSourceOf, removeCircularSourceCard,
    ...createCircularSourceMoves({ state, discoverDrawingRecords, shownAction, runUndoable })
  };
};

/**
 * @typedef {object} CircularRecordDisplayPorts
 * @property {WebState} state
 * @property {{ readonly value: boolean }} semanticMutationAvailable No session operation holds the inputs.
 * @property {{ readonly value: { status: string, inputs: readonly { uid: string, primaryFile: InputFileRef | null, pairedFile: InputFileRef | null }[] } }} circularRecordDiscoveryState
 *   The shown drawing's Circular record discovery.
 * @property {{ readonly value: boolean }} circularSingleRecordPresentationEnabled
 *   Whether the shown Circular drawing presents one record.
 * @property {() => SessionBusy | null} sessionOperationAvailability
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {() => void} rebindAnnotationRecords Rebinds the shown drawing's annotations to its record catalog (R10).
 * @property {RevealSetting} revealSetting Brings a setting into view.
 * @property {{ value: boolean }} circularTrackSlotsPanelOpen Whether the Custom Track Slots panel is open.
 * @property {{ readonly value: HTMLDetailsElement | null }} circularRecordPresentationPanel The record presentation panel.
 * @property {() => Promise<unknown>} afterPaint Settles after the next paint.
 */

/**
 * @typedef {object} RecordPresentationScrollPorts
 * @property {() => unknown[]} presentation The mode, whether one record is presented, the cards' files
 *   (`inputKey`, third) and the record selector: a change of any opens the panel.
 * @property {() => string} inputKey The shown drawing's cards and their files.
 * @property {{ readonly value: boolean }} circularSingleRecordPresentationEnabled
 * @property {{ readonly value: HTMLDetailsElement | null }} circularRecordPresentationPanel
 * @property {() => Promise<unknown>} afterPaint Settles after the next paint.
 */

/**
 * The scroll anchoring when the record presentation panel appears: the Settings pane keeps the
 * focused control where it was (or its scroll position) once the panel has opened and painted.
 * @param {RecordPresentationScrollPorts} ports
 * @returns {() => Promise<unknown>} Settles once the pane is back where every change so far found
 *   it; a capture tool awaits it after an operation instead of counting the frames it takes.
 */
const anchorScrollOnRecordPresentation = ({
  presentation, inputKey, circularSingleRecordPresentationEnabled, circularRecordPresentationPanel, afterPaint
}) => {
  /** @type {Promise<unknown>} */
  let anchoring = Promise.resolve();
  /** @param {unknown[]} current */
  const anchor = async (current) => {
    const origin = document.activeElement;
    const pane = document.querySelector('.settings-scroll');
    const scrollOwner = pane ? settingsScrollOwner(pane) : document.scrollingElement;
    const anchored = origin && pane?.contains(origin) && origin.getClientRects().length;
    const top = anchored ? origin.getBoundingClientRect().top : null;
    const scrollTop = scrollOwner?.scrollTop;
    await nextTick();
    if (!circularSingleRecordPresentationEnabled.value || inputKey() !== current[2]) return;
    if (circularRecordPresentationPanel.value) circularRecordPresentationPanel.value.open = true;
    await afterPaint();
    if (scrollOwner && document.activeElement === origin) {
      // `top` is a number whenever `anchored` is truthy; `scrollTop` was read from
      // the same non-null `scrollOwner`.
      scrollOwner.scrollTop = anchored && origin.isConnected
        ? scrollOwner.scrollTop + origin.getBoundingClientRect().top - /** @type {number} */ (top)
        : /** @type {number} */ (scrollTop);
    }
  };
  watch(presentation, (/** @type {unknown[]} */ current, /** @type {unknown[]} */ previous = []) => {
    if (!current[1] || current.every((value, index) => Object.is(value, previous[index]))) return undefined;
    const run = anchor(current);
    // Vue reports a failed run, as it did when the watch itself was async; the wait still settles.
    anchoring = Promise.all([anchoring, run.catch(() => undefined)]).then(() => undefined);
    return run;
  }, { flush: 'sync' });
  return () => anchoring;
};

/**
 * Registers the `circularRecordSelectionEnabled` computed and the scroll
 * anchoring watch; the root calls it where they stood.
 * @param {CircularRecordDisplayPorts} ports
 */
export const createCircularRecordDisplay = ({
  state, semanticMutationAvailable, circularRecordDiscoveryState, circularSingleRecordPresentationEnabled,
  sessionOperationAvailability, shownDrawingOfMode, rebindAnnotationRecords, revealSetting, circularTrackSlotsPanelOpen,
  circularRecordPresentationPanel, afterPaint
}) => {
  const { mode } = state;
  /**
   * A record chosen alone in Record select takes the crop its drawing keeps for it in a grid or batch into the
   * single-record crop fields, where it can be changed (PD-OI-044 scenario revision 5): its one window in Start/End
   * (or the list, across the origin), and none for a record without one. The fields keep a change while they hold
   * that record's crop (`circular_range_fields_record`), through All records too, so the changed crop draws when the
   * record is chosen alone again; another record refills them. A drawing that keeps no such crop keeps its fields.
   * @param {DrawingState} drawing @param {string} identity
   */
  const fillRecordRange = (drawing, identity) => {
    const { form } = drawing;
    if (!form.circular_record_ranges.length || identity === form.circular_range_fields_record) return;
    const entries = circularRecordEntries(drawing.circularRecordList.value);
    const entry = entries.find((candidate) => candidate.value === identity);
    const windows = (entry ? circularRecordRange(form.circular_record_ranges, entries, entry)?.windows : null) ?? [];
    const only = cardCropWindow(windows);
    form.circular_region_start = only ? only.start : null;
    form.circular_region_end = only ? only.end : null;
    form.circular_reverse = only ? only.reverseComplement : false;
    form.circular_windows = only ? [] : windows.map((window) => ({ ...window }));
    form.circular_range_fields_record = identity;
  };
  const circularRecordSelectionEnabled = computed(() => semanticMutationAvailable.value
    && mode.value === 'circular'
    && !state.activeDrawing().form.multi_record_canvas && circularRecordDiscoveryState.value.status === 'ready');
  /** @param {unknown} value */
  const setCircularRecordPresentationSelector = (value) => {
    const busy = sessionOperationAvailability();
    if (busy) return busy;
    const drawing = shownDrawingOfMode('circular');
    if (!drawing) return { status: 'unavailable' };
    if (!circularRecordSelectionEnabled.value) return { status: 'unavailable' };
    const normalized = String(value || '').trim();
    if (normalized) fillRecordRange(drawing, normalized);
    drawing.form.circular_record_selector = normalized;
    drawing.adv.circular_grouping_intent = normalized ? 'single' : 'auto';
    rebindAnnotationRecords();
    return { status: 'ok' };
  };
  const showCircularCanvasSetting = async () => {
    const drawing = shownDrawingOfMode('circular');
    const control = /** @type {HTMLElement | null} */ (document.querySelector('[data-circular-canvas-setting]'));
    if (!control || !drawing || !drawing.form.multi_record_canvas) {
      return { status: 'unavailable' };
    }
    await revealSetting(control);
    return { status: 'ok' };
  };
  // With Custom Track Slots on, a Ticks row owns the coordinate scale: the link beside Show Coordinate Scale
  // reveals that row (the first enabled one, else the first), or the Custom Track Slots panel when the stack
  // has none. It changes no setting.
  const showCircularTicksSlot = async () => {
    const drawing = state.activeDrawing();
    if (mode.value !== 'circular' || !drawing.adv.circular_track_slots_enabled) return { status: 'unavailable' };
    circularTrackSlotsPanelOpen.value = true;
    await nextTick();
    /** @type {CircularTrackSlot[]} */
    const slots = Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots : [];
    const ticks = slots.filter((slot) => slot.renderer === 'ticks');
    const slot = ticks.find((candidate) => candidate.enabled) || ticks[0];
    const row = slot ? document.querySelector(`[data-capture="circular-track-slot-${CSS.escape(slot.id)}"]`) : null;
    const panel = document.getElementById('circular-custom-track-slots-panel');
    const toggle = /** @type {HTMLElement | null} */ (document.querySelector('[data-track-slots-toggle]'));
    const control = /** @type {HTMLElement | null} */ (row?.querySelector('input, select, button') ?? null);
    if (row instanceof HTMLElement && control) {
      await revealSetting(row, { focus: control });
    } else if (panel && toggle) {
      await revealSetting(panel, { align: 'start', focus: toggle });
    } else {
      return { status: 'unavailable' };
    }
    return { status: 'ok' };
  };
  // The shown drawing's cards and their files, as one value a watch compares.
  const inputKey = () => circularRecordDiscoveryState.value.inputs
    .map(({ uid, primaryFile, pairedFile }) => `${uid}/${primaryFile?.fileId ?? ''}/${pairedFile?.fileId ?? ''}`).join('|');
  const waitForCircularRecordScrollAnchor = anchorScrollOnRecordPresentation({
    presentation: () => [mode.value, circularSingleRecordPresentationEnabled.value, inputKey(),
      state.activeDrawing().form.circular_record_selector],
    inputKey, circularSingleRecordPresentationEnabled, circularRecordPresentationPanel, afterPaint
  });
  return {
    circularRecordSelectionEnabled, setCircularRecordPresentationSelector, showCircularCanvasSetting,
    showCircularTicksSlot, waitForCircularRecordScrollAnchor
  };
};
