// @ts-check
// The ports the replacement of an input file takes (its owner is
// `createInputReplacement` in app/drawing-list.js): the writer of each slot a
// replacement names, applied to the sites in their drawings
// (`applyInputReplacementToDrawings`), and the question whether the shown
// drawing only or every drawing that uses the file takes the new file
// (`askScope`). The composition root passes the state, the slot writers of
// the source, Depth, comparison-ring and comparison edits, the drawing bar's
// tabs, and its choice dialog as ports.
import { circularDepthRows } from './depth-tracks.js';
import { normalizeFileList } from '../services/conservation-series.js';
import { depthFileSlotsFromValue } from '../services/depth-track-state.js';
import { listText, t, tn } from '../utils/i18n.js';
/** @import { DrawingState, LinearSeq, WebState } from '../state.js' */
/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { InputFileCandidate, InputFileRef } from '../services/history-files.js' */
/** @import { ReplaceableInputSlot } from '../services/drawing/inputs.js' */
/** @import { AskChoice } from './app-messages.js' */
/** @import { ApplyInputReplacementToDrawings, DrawingBarOffer, DrawingBarTab, DrawingInputSite, InputReplacementScope } from './drawing-list.js' */

/**
 * @typedef {InputFileCandidate | InputFileRef} ReplacementFile
 * @typedef {object} InputReplacementPortsPorts
 * @property {WebState} state
 * @property {(drawing: DrawingState, field: 'c_gb' | 'c_gff' | 'c_fasta', value: ReplacementFile) => unknown} replaceCircularSource
 * @property {(drawing: DrawingState, index: number, file: ReplacementFile) => unknown} replaceCircularDepthFile
 * @property {(drawing: DrawingState, nextFiles: readonly (InputFileRef | null)[]) => unknown} replaceCircularConservationBlasts
 * @property {(drawing: DrawingState, rowUid: string, field: 'gb' | 'gff' | 'fasta', value: ReplacementFile) => unknown} replaceLinearSourceFile
 * @property {(drawing: DrawingState, sequences: LinearSeq[], index: number, file: ReplacementFile) => unknown} setLinearDepthFiles
 * @property {(drawing: DrawingState, id: string, file: ReplacementFile) => unknown} swapLinearComparisonFile
 * @property {() => Array<DrawingBarTab | DrawingBarOffer>} barTabs The drawing bar's tabs, named as they show.
 * @property {AskChoice} askChoice
 */

/**
 * Registers nothing. The root creates it right before the input replacement,
 * once every writer it reaches exists.
 * @param {InputReplacementPortsPorts} ports
 * @returns {{ applyInputReplacementToDrawings: ApplyInputReplacementToDrawings,
 *   askScope: (previous: InputFileRef, otherDrawingIds: readonly DrawingId[]) => Promise<InputReplacementScope | null> }}
 */
export const createInputReplacementPorts = ({
  state, replaceCircularSource, replaceCircularDepthFile, replaceCircularConservationBlasts, replaceLinearSourceFile,
  setLinearDepthFiles, swapLinearComparisonFile, barTabs, askChoice
}) => {
  /** @param {readonly (InputFileRef | null)[]} slots @param {string} fileId @returns {number[]} */
  const slotIndexesOfFile = (slots, fileId) => slots.flatMap((ref, idx) => (ref && ref.fileId === fileId ? [idx] : []));
  /**
   * The writer of each slot an input file replacement names: every
   * `ReplaceableInputSlot` has one (a missing one fails the type check), and
   * each is the body the shown drawing's control of that slot runs.
   * @type {Readonly<Record<ReplaceableInputSlot, (drawing: DrawingState, site: DrawingInputSite, fileId: string, file: InputFileCandidate | InputFileRef) => void>>}
   */
  const INPUT_REPLACEMENT_WRITERS = Object.freeze({
    c_gb: (drawing, _site, _fileId, file) => replaceCircularSource(drawing, 'c_gb', file),
    c_gff: (drawing, _site, _fileId, file) => replaceCircularSource(drawing, 'c_gff', file),
    c_fasta: (drawing, _site, _fileId, file) => replaceCircularSource(drawing, 'c_fasta', file),
    c_depth: (drawing, _site, fileId, file) => {
      const indexes = new Set(circularDepthRows(drawing).flatMap((row) => slotIndexesOfFile(row, fileId)));
      indexes.forEach((idx) => replaceCircularDepthFile(drawing, idx, file));
    },
    c_conservation_blasts: (drawing, _site, fileId, file) => {
      const ref = state.inputFiles.registerFile(file);
      replaceCircularConservationBlasts(drawing, normalizeFileList(drawing.files.c_conservation_blasts)
        .map((entry) => (entry.fileId === fileId ? ref : entry)));
    },
    'row.gb': (drawing, site, _fileId, file) => replaceLinearSourceFile(drawing, String(site.rowUid), 'gb', file),
    'row.gff': (drawing, site, _fileId, file) => replaceLinearSourceFile(drawing, String(site.rowUid), 'gff', file),
    'row.fasta': (drawing, site, _fileId, file) => replaceLinearSourceFile(drawing, String(site.rowUid), 'fasta', file),
    'row.depth': (drawing, site, fileId, file) => {
      const seq = drawing.linearSeqs.find((entry) => entry.uid === site.rowUid);
      if (!seq) return;
      slotIndexesOfFile(depthFileSlotsFromValue(seq.depth), fileId).forEach((idx) => setLinearDepthFiles(drawing, [seq], idx, file));
    },
    edge: (drawing, site, _fileId, file) => swapLinearComparisonFile(drawing, String(site.edgeId), file)
  });
  /** @type {ApplyInputReplacementToDrawings} */
  const applyInputReplacementToDrawings = (sites, fileId, file) => sites.forEach((site) => {
    INPUT_REPLACEMENT_WRITERS[site.slot](state.drawingOf(site.drawingId), site, fileId, file);
  });
  /**
   * @param {InputFileRef} previous
   * @param {readonly DrawingId[]} otherDrawingIds
   * @returns {Promise<InputReplacementScope | null>}
   */
  const askScope = async (previous, otherDrawingIds) => {
    // Named as their tabs show them: a default name in the UI language.
    const others = barTabs().flatMap((tab) => (
      tab.kind === 'drawing' && otherDrawingIds.includes(tab.id) ? [tab.label] : []));
    const choice = await askChoice({
      title: t('Replace input file?'),
      description: tn(others.length, '{name} is also used in {drawings}.', '{name} is also used in {count} drawings: {drawings}.', {
        name: previous.name, drawings: listText(others)
      }),
      actions: [
        { id: 'drawing', label: t('This drawing only'), tone: 'secondary' },
        { id: 'all', label: t('Every drawing that uses this file'), tone: 'primary' }
      ]
    });
    return choice === 'drawing' || choice === 'all' ? choice : null;
  };

  return { applyInputReplacementToDrawings, askScope };
};
