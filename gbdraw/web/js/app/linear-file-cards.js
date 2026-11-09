// @ts-check
// The Linear File cards of the shown drawing: the cards (one per source), the
// header line of a card, whether a card holds the files its input type reads,
// a card's fold, and the File removal dialog's target. The composition root
// passes the shown-drawing readers and the folded cards as ports; the removal
// itself stays in app-setup.js.
import { circularInputComplete, circularSourceInputs, inputSourceFiles } from '../services/input-types.js';
import { groupLinearSourceRecords, linearCardInputComplete, setLinearSourceFolded } from '../services/linear-sources.js';
import { t } from '../utils/i18n.js';
import { chosenFileName, fileCardLine } from './file-card.js';
/** @import { DiagramMode, DrawingState, LinearSeq } from '../state.js' */
/** @import { InputSlot } from '../services/input-types.js' */
/** @import { FoldedLinearSources, LinearSourceGroup } from '../services/linear-sources.js' */

const { computed, reactive, ref } = window.Vue;

/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */

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
  /** A File card's header line (`fileCardLine`). @param {LinearSourceGroup} source */
  const linearCardFileLine = (source) => {
    const drawing = shownDrawingOfMode('linear');
    return drawing === null ? t('No file') : fileCardLine(drawing.inputType.value, source.sequence);
  };
  /**
   * Whether the shown drawing of a mode has every file its input type reads.
   * @param {'circular' | 'linear'} drawingMode
   * @param {Readonly<Partial<Record<InputSlot, unknown>>>} [linearSlots] a Linear card's slots; Circular reads the drawing's own cards
   */
  const shownInputComplete = (drawingMode, linearSlots = {}) => {
    const drawing = shownDrawingOfMode(drawingMode);
    if (drawing === null) return false;
    return drawingMode === 'circular'
      ? circularInputComplete(circularSourceInputs(drawing.inputType.value, drawing.circularSources))
      : inputSourceFiles(drawing.inputType.value, linearSlots).complete;
  };
  /** @param {LinearSourceGroup} source */
  const linearCardHoldsInput = (source) => {
    const drawing = shownDrawingOfMode('linear');
    return drawing !== null && linearCardInputComplete(source, drawing.inputType.value);
  };
  // A File card folds to its header line. The disclosure reports every change
  // of its own state (the shell's `fold`): the user's press, and a path that
  // opened the card to show a setting inside it.
  /** @param {LinearSourceGroup} source */
  const isLinearSourceFolded = (source) => linearFoldedSources.has(source.uid);
  /** @param {LinearSourceGroup} source @param {boolean} folded */
  const syncLinearSourceFold = (source, folded) => setLinearSourceFolded(linearFoldedSources, source.uid, folded);
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
