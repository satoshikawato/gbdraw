// @ts-check
// The Linear File cards of the shown drawing: the cards (one per source), the
// header line of a card, whether a card holds the files its input type reads,
// a card's fold, and the File removal dialog's target. The composition root
// passes the shown-drawing readers and the folded cards as ports; the removal
// itself stays in app-setup.js.
import { circularInputSlots, inputSourceFiles } from '../services/input-types.js';
import {
  groupLinearSourceRecords,
  linearCardFileFields,
  linearCardInputComplete,
  setLinearSourceFolded
} from '../services/linear-sources.js';
import { t } from '../utils/i18n.js';
/** @import { DiagramMode, DrawingState, LinearSeq } from '../state.js' */
/** @import { InputSlot } from '../services/input-types.js' */
/** @import { FoldedLinearSources, LinearSourceGroup } from '../services/linear-sources.js' */

const { computed, reactive, ref } = window.Vue;

/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */

/** @param {{ name?: string } | null} file */
const chosenFileName = (file) => String(file?.name || t('Unnamed file'));
const LINEAR_FILE_SLOT_NAMES = { gb: 'GenBank', gff: 'GFF3', fasta: 'FASTA' };

/**
 * @typedef {object} LinearFileCardsPorts
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {() => readonly LinearSeq[]} shownLinearSeqs The record rows of the shown Linear drawing.
 * @property {FoldedLinearSources} linearFoldedSources The File cards folded to their header line.
 */

/**
 * The File cards of the shown Linear drawing and the File removal dialog's target.
 * @param {LinearFileCardsPorts} ports
 */
export const createLinearFileCards = ({ shownDrawingOfMode, shownLinearSeqs, linearFoldedSources }) => {
  /** @type {{ readonly value: LinearSourceGroup[] }} */
  const linearSourceGroups = computed(() => groupLinearSourceRecords(shownLinearSeqs()));
  const linearSourceRemovalDialog = reactive({ open: false, sourceUid: '', origin: '' });
  const linearSourceRemovalReturnFocus = ref(null);
  const linearSourceRemovalTarget = computed(() => (
    linearSourceGroups.value.find((source) => source.uid === linearSourceRemovalDialog.sourceUid) || null
  ));
  const linearSourceRemovalCanDelete = computed(() => linearSourceGroups.value.length > 1);
  /**
   * The file names of a Linear File, as its dialogs name it: every file the
   * row holds, which a removal takes, also one of the other input type.
   * @param {LinearSourceGroup} source
   */
  const linearSourceFileName = (source) => {
    const names = [source.sequence.gb, source.sequence.gff, source.sequence.fasta].filter(Boolean).map(chosenFileName);
    return names.length ? names.join(' + ') : t('File {n}', { n: linearSourceGroups.value.indexOf(source) + 1 });
  };
  /**
   * What a File card's header line says of the slots the card shows: each
   * file's name, "no FASTA" for a slot still empty, "No file" for none.
   * @param {LinearSourceGroup} source
   */
  const linearCardFileLine = (source) => {
    const drawing = shownDrawingOfMode('linear');
    if (drawing === null) return t('No file');
    const fields = linearCardFileFields(drawing.inputType.value);
    if (!fields.some((field) => source.sequence[field])) return t('No file');
    return fields.map((field) => (
      source.sequence[field] ? chosenFileName(source.sequence[field]) : t('no {file}', { file: LINEAR_FILE_SLOT_NAMES[field] })
    )).join(' + ');
  };
  /**
   * Whether the shown drawing of a mode has every file its input type reads.
   * @param {'circular' | 'linear'} drawingMode
   * @param {Readonly<Partial<Record<InputSlot, unknown>>>} [linearSlots] a Linear card's slots; Circular reads the drawing's own files
   */
  const shownInputComplete = (drawingMode, linearSlots = {}) => {
    const drawing = shownDrawingOfMode(drawingMode);
    if (drawing === null) return false;
    const slots = drawingMode === 'circular' ? circularInputSlots(drawing.files) : linearSlots;
    return inputSourceFiles(drawing.inputType.value, slots).complete;
  };
  /** @param {LinearSourceGroup} source */
  const linearCardHoldsInput = (source) => {
    const drawing = shownDrawingOfMode('linear');
    return drawing !== null && linearCardInputComplete(source, drawing.inputType.value);
  };
  // A File card folds to its header line. The disclosure reports every change
  // of its own state: the user's press, and a path that opened the card to
  // show a setting inside it.
  /** @param {LinearSourceGroup} source */
  const isLinearSourceFolded = (source) => linearFoldedSources.has(source.uid);
  /** @param {LinearSourceGroup} source @param {Event} event */
  const syncLinearSourceFold = (source, event) => setLinearSourceFolded(
    linearFoldedSources, source.uid, !(/** @type {HTMLDetailsElement} */ (event.currentTarget)).open
  );
  const linearSourceRemovalTargetName = computed(() => {
    const source = linearSourceRemovalTarget.value;
    return source ? linearSourceFileName(source) : t('Unavailable File');
  });
  return {
    linearSourceGroups, linearSourceRemovalDialog, linearSourceRemovalReturnFocus, linearSourceRemovalTarget,
    linearSourceRemovalCanDelete, linearSourceFileName, linearCardFileLine, shownInputComplete, linearCardHoldsInput,
    isLinearSourceFolded, syncLinearSourceFold, linearSourceRemovalTargetName
  };
};
