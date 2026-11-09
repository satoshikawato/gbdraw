// @ts-check
// The Circular File cards of the shown drawing (one card per source file, or
// per GFF3 + FASTA pair, in drawing order): each card's header line, state and
// records line, its fold, the place a record's Topology card names, the Add
// files zone (several files of one pick or drop become cards), the read error
// that names its card, and the File removal dialog's target. The cards show in the shell the Linear
// File cards use (app/file-card.js). The edits are the source card writers'
// (app/circular-record-controls.js), passed as ports; this owner writes no
// drawing.
import { circularSlotFileKind } from '../services/drawing/inputs.js';
import { circularSourceInputs, inputTypeSpec } from '../services/input-types.js';
import { formatRecordLength } from '../services/record-options.js';
import { errorText } from '../utils/error-text.js';
import { t, tn } from '../utils/i18n.js';
import { fileCardLine } from './file-card.js';
/** @import { DiagramMode, DrawingState, SessionBusy } from '../state.js' */
/** @import { CircularSource } from '../services/drawing/inputs.js' */
/** @import { InputFileCandidate, InputFileRef } from '../services/history-files.js' */
/** @import { InputType } from '../services/input-types.js' */
/** @import { CircularListRecord } from '../services/record-discovery.js' */
/** @import { InputFileKind } from '../services/drawing/inputs.js' */
/** @import { UserFacingError } from '../utils/error-normalization.js' */

const { computed, reactive } = window.Vue;

/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */

/**
 * @typedef {object} CircularFileCardsPorts
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {{ readonly value: { status: string, error: UserFacingError | '', failedSources: readonly string[] } }} circularRecordDiscoveryState
 *   The shown drawing's record catalog: its status, the first failed card's error, and the cards whose read failed.
 * @property {(files: readonly (InputFileCandidate | InputFileRef)[]) => unknown} addCircularSourceFiles
 *   New cards at the end of the list from the files of one pick or drop: one History step.
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

/**
 * What the Add files zone names, takes and offers from this Session.
 * @typedef {object} CircularAddZone
 * @property {string} label
 * @property {readonly InputFileKind[]} fileKinds The kinds of the Session files it offers while the list is empty:
 *   those of every slot the input type reads (a GFF3 + FASTA zone offers both, and a chosen one goes to the slot of
 *   its kind, as the empty list's two uploaders offered them); none once a card shows, whose own uploaders offer
 *   them, so the zone lists no file already in the list.
 * @property {string} accept
 */
const COMMON_ACCEPT = '.txt,text/plain,text/*,.gz,application/gzip';
/** @type {Readonly<Record<InputType, { label: () => string, accept: string }>>} */
const ADD_ZONES = Object.freeze({
  gb: { label: () => t('GenBank/DDBJ files'), accept: `.gb,.gbk,.gbff,.genbank,.ddbj,.embl,.emb,${COMMON_ACCEPT}` },
  gff: { label: () => t('GFF3 and FASTA files'), accept: `.gff,.gff3,.fa,.fas,.fasta,.fna,.ffn,.faa,${COMMON_ACCEPT}` },
  fasta: { label: () => t('FASTA files'), accept: `.fa,.fas,.fasta,.fna,.ffn,${COMMON_ACCEPT}` }
});
/** @param {InputType} inputType @returns {InputFileKind[]} the kinds of the slots the input type reads */
const fileKindsOf = (inputType) => {
  const { primary, paired } = inputTypeSpec(inputType);
  return (paired ? [primary, paired] : [primary]).map(circularSlotFileKind);
};

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
  /** @type {{ readonly value: CircularAddZone }} */
  const addZone = computed(() => {
    const drawing = shownDrawingOfMode('circular');
    const inputType = drawing ? drawing.inputType.value : 'gb';
    const zone = ADD_ZONES[inputType];
    return { label: zone.label(), accept: zone.accept, fileKinds: drawing && drawing.circularSources.length ? [] : fileKindsOf(inputType) };
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
   * The files of one pick or drop on the Add files zone become new cards (one
   * History step); a new card that still waits for a file shows open. The zone
   * stays where it is, above the cards, so the focus stays on it.
   * @param {readonly (InputFileCandidate | InputFileRef)[]} files
   */
  const addFiles = async (files) => {
    if (files.length === 0) return undefined;
    const before = new Set(cards.value.map((card) => card.uid));
    const outcome = await addCircularSourceFiles(files);
    cards.value.filter((card) => !before.has(card.uid) && !card.ready).forEach((card) => openCards.add(card.uid));
    return outcome;
  };
  /**
   * The read error the Source records line shows: the first failed card's,
   * named by that card's files when the list holds several cards, so a wrong
   * kind of file says which file it is.
   */
  const readError = computed(() => {
    const discovery = circularRecordDiscoveryState.value;
    if (discovery.status !== 'error') return '';
    const message = (discovery.error && errorText(discovery.error.summaryParts)) || t('Source inspection failed. Retry or replace the input.');
    const failed = cards.value.find((card) => card.failed);
    return failed && cards.value.length > 1 ? t('{file}: {message}', { file: failed.line, message }) : message;
  });

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
    addZone, cards, isFolded, setFolded, addFiles, readError, recordPlace, removal, removalTarget, requestRemoval,
    cancelRemoval, applyRemoval
  };
};
