// @ts-check
// The shown drawing's input files and input type: the uploaders of its bound
// input files (each names the file its control binds and its own write, and
// the input replacement decides whether to ask), and the Circular and Linear
// input type. The composition root passes the shown-drawing readers and
// admission, the input replacement, the source file writers, the record
// discovery, and the alignment owner's clear as ports.
import { normalizeInputType } from '../services/input-types.js';
import { msg } from '../utils/i18n.js';
import { circularDepthRepresentatives, getLinearDepthFile } from './depth-tracks.js';
/** @import { DiagramMode, DrawingState, LinearSeq, SessionBusy } from '../state.js' */
/** @import { InputFileCandidate, InputFileRef } from '../services/history-files.js' */
/** @import { LinearSourceGroup } from '../services/linear-sources.js' */
/** @import { LinearComparisonPair } from '../services/linear-comparisons.js' */

/**
 * A mode panel's mutating handler: the session's busy answer, nothing when the
 * shown drawing is of another mode, else `action` on the shown drawing.
 * @typedef {<T>(drawingMode: DiagramMode, action: (drawing: DrawingState) => T) => T | SessionBusy | undefined} ShownAction
 */
/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */
/** @typedef {() => SessionBusy | null} SessionOperationAvailability The session's busy answer, or null when an edit may run. */
/** @typedef {(reason: string) => unknown} ClearAlignmentForMutation The alignment owner clears its plan; `reason` is a `msg` clause. */
/** @typedef {(drawing: DrawingState) => unknown} DiscoverDrawingRecords Reads the records of the drawing's inputs. */

/**
 * @typedef {object} InputUploadsPorts
 * @property {ShownAction} shownAction
 * @property {<T>(previous: InputFileRef | null, file: InputFileCandidate | InputFileRef | null, writeShown: () => T) => T | Promise<unknown>} replaceFromShownDrawing
 *   Writes `file` in place of `previous` with the uploader's own write, after asking when another drawing uses `previous`.
 * @property {(drawing: DrawingState, field: 'c_gb' | 'c_gff' | 'c_fasta', value: InputFileCandidate | InputFileRef) => void} replaceCircularSource
 * @property {(drawing: DrawingState, index: number, file: InputFileCandidate | InputFileRef) => unknown} replaceCircularDepthFile
 * @property {(drawing: DrawingState, rowUid: string, field: 'gb' | 'gff' | 'fasta', value: InputFileCandidate | InputFileRef) => unknown} replaceLinearSourceFile
 * @property {(drawing: DrawingState, sequences: LinearSeq[], index: number, file: InputFileCandidate | InputFileRef | null) => unknown} setLinearDepthFiles
 * @property {(source: LinearSourceGroup, index: number, file: File | null) => unknown} setLinearSourceDepthFile
 * @property {(edgeKey: string) => LinearComparisonPair | null} linearComparisonPairForEdgeKey
 * @property {(edgeKey: string, file: File | null) => unknown} setLinearComparisonCardFile
 */

/**
 * The uploaders of the shown drawing's bound input files.
 * @param {InputUploadsPorts} ports
 */
export const createInputUploads = ({
  shownAction, replaceFromShownDrawing, replaceCircularSource, replaceCircularDepthFile, replaceLinearSourceFile,
  setLinearDepthFiles, setLinearSourceDepthFile, linearComparisonPairForEdgeKey, setLinearComparisonCardFile
}) => {
  /** @param {'c_gb' | 'c_gff' | 'c_fasta'} field @param {File | null} file */
  const uploadCircularSourceFile = (field, file) => shownAction('circular', (drawing) => (
    replaceFromShownDrawing(drawing.files[field], file, () => replaceCircularSource(drawing, field, file))
  ));
  /** @param {number} index @param {File | null} file */
  const uploadCircularDepthFile = (index, file) => shownAction('circular', (drawing) => (
    replaceFromShownDrawing(circularDepthRepresentatives(drawing)[Number(index)] || null, file, () => replaceCircularDepthFile(drawing, index, file))
  ));
  /** @param {number} index @param {'gb' | 'gff' | 'fasta'} field @param {File | null} file */
  const uploadLinearSeqPrimaryFile = (index, field, file) => shownAction('linear', (drawing) => {
    const seq = drawing.linearSeqs[Number(index)];
    return seq ? replaceFromShownDrawing(seq[field], file, () => replaceLinearSourceFile(drawing, seq.uid, field, file)) : undefined;
  });
  /** @param {LinearSeq} seq @param {number} index @param {File | null} file */
  const uploadLinearDepthFile = (seq, index, file) => shownAction('linear', (drawing) => (
    replaceFromShownDrawing(getLinearDepthFile(seq, index), file, () => setLinearDepthFiles(drawing, [seq], index, file))
  ));
  /** @param {LinearSourceGroup} source @param {number} index @param {File | null} file */
  const uploadLinearSourceDepthFile = (source, index, file) => shownAction('linear', () => (
    replaceFromShownDrawing(getLinearDepthFile(source.sequence, index), file, () => setLinearSourceDepthFile(source, index, file))
  ));
  /** @param {string} edgeKey @param {File | null} file */
  const uploadLinearComparisonCardFile = (edgeKey, file) => shownAction('linear', () => {
    const pair = linearComparisonPairForEdgeKey(edgeKey);
    return replaceFromShownDrawing(pair && pair.draft ? pair.draft.file : null, file, () => setLinearComparisonCardFile(edgeKey, file));
  });
  return {
    uploadCircularSourceFile, uploadCircularDepthFile, uploadLinearSeqPrimaryFile, uploadLinearDepthFile,
    uploadLinearSourceDepthFile, uploadLinearComparisonCardFile
  };
};

/**
 * @typedef {object} InputTypeEditsPorts
 * @property {SessionOperationAvailability} sessionOperationAvailability
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {<T>(drawing: DrawingState, mutate: () => T) => T} changeRecordSources
 *   The transition of every write of a drawing's row source.
 * @property {DiscoverDrawingRecords} discoverDrawingRecords
 * @property {ClearAlignmentForMutation} clearAlignmentForMutation
 */

/**
 * The input type of the shown Circular or Linear drawing.
 * @param {InputTypeEditsPorts} ports
 */
export const createInputTypeEdits = ({
  sessionOperationAvailability, shownDrawingOfMode, changeRecordSources, discoverDrawingRecords, clearAlignmentForMutation
}) => {
  /** @param {unknown} value */
  const setCircularInputType = (value) => {
    const busy = sessionOperationAvailability();
    if (busy) return busy;
    const drawing = shownDrawingOfMode('circular');
    if (!drawing) return { status: 'unavailable' };
    const next = normalizeInputType(value);
    if (drawing.inputType.value === next) return { status: 'ok' };
    changeRecordSources(drawing, () => { drawing.inputType.value = next; });
    // The other input's records are not the ones the OFF list names.
    if (drawing.recordsOff.length) drawing.recordsOff.splice(0);
    void discoverDrawingRecords(drawing);
    return { status: 'ok' };
  };

  /** @param {unknown} value */
  const setLinearInputType = (value) => {
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const drawing = shownDrawingOfMode('linear');
    if (!drawing) return false;
    const next = normalizeInputType(value);
    if (drawing.inputType.value === next) return false;
    clearAlignmentForMutation(msg('source type changed.'));
    changeRecordSources(drawing, () => { drawing.inputType.value = next; });
    void discoverDrawingRecords(drawing);
    return true;
  };
  return { setCircularInputType, setLinearInputType };
};
