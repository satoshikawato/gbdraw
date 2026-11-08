// @ts-check
// The Circular File cards of the shown drawing (one card per source file, or
// per GFF3 + FASTA pair, in drawing order): each card's header line, state and
// records line, its fold, the place a record's Topology card names, Add file,
// and the File removal dialog's target. The cards show in the shell the Linear
// File cards use (app/file-card.js). The edits are the source card writers'
// (app/circular-record-controls.js), passed as ports; this owner writes no
// drawing.
import { circularSourceInputs } from '../services/input-types.js';
import { linearCardFileFields } from '../services/linear-sources.js';
import { formatRecordLength } from '../services/record-options.js';
import { tn } from '../utils/i18n.js';
import { fileCardLine } from './file-card.js';
/** @import { DiagramMode, DrawingState, SessionBusy } from '../state.js' */
/** @import { CircularSource } from '../services/drawing/inputs.js' */
/** @import { InputFileCandidate, InputFileRef } from '../services/history-files.js' */
/** @import { InputSlot } from '../services/input-types.js' */
/** @import { CircularListRecord } from '../services/record-discovery.js' */

const { computed, reactive } = window.Vue;

/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */

/**
 * @typedef {object} CircularFileCardsPorts
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {{ readonly value: { status: string, failedSources: readonly string[] } }} circularRecordDiscoveryState
 *   The shown drawing's record catalog: its status, and the cards whose read failed.
 * @property {(field: InputSlot, files: readonly (InputFileCandidate | InputFileRef)[]) => unknown} addCircularSourceFiles
 *   New cards at the end of the list, one per file, each with its file in `field`: one History step.
 * @property {(uid: string) => unknown} removeCircularSourceCard Removes a card, one History step.
 */

/**
 * One card as the list shows it.
 * @typedef {object} CircularFileCardView
 * @property {string} uid
 * @property {number} index
 * @property {CircularSource} source The card's file slots.
 * @property {string} line The header line (`fileCardLine`).
 * @property {boolean} ready Every slot the input type reads holds a file.
 * @property {boolean} failed Its read failed.
 * @property {CircularListRecord[]} records Its records, once read.
 * @property {string} recordsLine "1 record · 92,721 bp" once its records are read, else ''.
 */

/** @param {readonly CircularListRecord[]} records */
const recordsLineOf = (records) => {
  if (records.length === 0) return '';
  const count = tn(records.length, '{count} record', '{count} records');
  const lengths = records.map((record) => record.record_length);
  return lengths.every((length) => length !== null)
    ? `${count} · ${formatRecordLength(lengths.reduce((sum, length) => sum + Number(length), 0))}` : count;
};

/**
 * The File cards of the shown Circular drawing.
 * @param {CircularFileCardsPorts} ports
 */
export const createCircularFileCards = ({
  shownDrawingOfMode, circularRecordDiscoveryState, addCircularSourceFiles, removeCircularSourceCard
}) => {
  /** The file slots a card shows for the drawing's input type, in the order of its uploaders. */
  const fields = computed(() => {
    const drawing = shownDrawingOfMode('circular');
    return linearCardFileFields(drawing ? drawing.inputType.value : 'gb');
  });
  /** @type {{ readonly value: CircularFileCardView[] }} */
  const cards = computed(() => {
    const drawing = shownDrawingOfMode('circular');
    if (drawing === null) return [];
    const inputType = drawing.inputType.value;
    const discovery = circularRecordDiscoveryState.value;
    const read = discovery.status === 'ready' || discovery.status === 'error';
    /** @type {CircularListRecord[]} */
    const listed = read ? drawing.circularRecordList.value : [];
    const inputs = circularSourceInputs(inputType, drawing.circularSources);
    return drawing.circularSources.map((source, index) => {
      const failed = discovery.failedSources.includes(source.uid);
      const records = failed ? [] : listed.filter((record) => record.uid === source.uid);
      return {
        uid: source.uid, index, source, line: fileCardLine(inputType, source), ready: inputs[index].complete, failed,
        records, recordsLine: recordsLineOf(records)
      };
    });
  });

  // A card folds to its header line; it shows open only after the user opened
  // it, or when it was added still waiting for a file (a GFF3 without its
  // FASTA). Fold state is view state (no Session, History step, or Reset holds
  // it), by card uid so it follows a move.
  /** @type {Set<string>} */
  const openCards = reactive(new Set());
  /** @param {string} uid */
  const isFolded = (uid) => !openCards.has(uid);
  /** @param {string} uid @param {boolean} folded */
  const setFolded = (uid, folded) => {
    if (folded) openCards.delete(uid);
    else openCards.add(uid);
  };

  /**
   * New cards, one per file, each with its file in `field` (one History step);
   * a new card that still waits for a file shows open.
   * @param {InputSlot} field
   * @param {readonly (InputFileCandidate | InputFileRef)[]} files
   */
  const addFiles = async (field, files) => {
    const before = new Set(cards.value.map((card) => card.uid));
    const outcome = await addCircularSourceFiles(field, files);
    cards.value.filter((card) => !before.has(card.uid) && !card.ready).forEach((card) => openCards.add(card.uid));
    return outcome;
  };
  /**
   * Add file: the picked files become new cards in the slot the input type
   * reads first (GenBank, GFF3, or FASTA).
   * @param {Event} event the hidden file input's change
   */
  const addPickedFiles = (event) => {
    const input = /** @type {HTMLInputElement} */ (event.target);
    const files = Array.from(input.files || []);
    input.value = '';
    return files.length ? addFiles(fields.value[0], files) : undefined;
  };

  /**
   * The place a record's Topology card names it by: its file, and its `#k`
   * when the file holds more than one record (circular-multi-file D-10).
   * @param {{ sourceUid: string, selector: string }} row
   */
  const recordPlace = (row) => {
    const card = cards.value.find((entry) => entry.uid === row.sourceUid);
    if (!card) return row.selector;
    return card.records.length > 1 ? `${card.line} ${row.selector}` : card.line;
  };

  // The File removal dialog: the card it asks about, by uid.
  const removal = reactive({ open: false, uid: '' });
  const removalTarget = computed(() => (removal.open ? cards.value.find((card) => card.uid === removal.uid) || null : null));
  /** @param {string} uid */
  const requestRemoval = (uid) => Object.assign(removal, { open: true, uid });
  const cancelRemoval = () => Object.assign(removal, { open: false, uid: '' });
  /** @returns {Promise<SessionBusy | undefined | unknown>} */
  const applyRemoval = async () => {
    const target = removalTarget.value;
    if (!target) return undefined;
    const outcome = await removeCircularSourceCard(target.uid);
    cancelRemoval();
    return outcome;
  };

  return {
    fields, cards, isFolded, setFolded, addFiles, addPickedFiles, recordPlace, removal, removalTarget, requestRemoval,
    cancelRemoval, applyRemoval
  };
};
