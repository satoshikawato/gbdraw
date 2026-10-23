// @ts-check
// The input bindings: which stored file each input field names. A field holds
// an `InputFileRef` of the project's file store (`state.inputFiles`), never a
// File. The writers here take the object they write (the project's bindings
// now; a drawing's inputs slice later), so the store registers every bound
// file through one owner. The bound IDs keep their files stored: the store
// releases a file only when no History step and no binding names it.

/** @import { HistoryFile, InputFileCandidate, InputFileRef, InputFileStore } from '../history-files.js' */
import { isInputFileRef } from '../history-files.js';

/**
 * The bindings that name stored files: the input fields, the Linear record
 * rows, and each drawing's uploaded comparison files.
 * @typedef {object} InputBindings
 * @property {Readonly<Record<string, unknown>>} files
 * @property {readonly Readonly<Record<string, unknown>>[]} linearSeqs
 * @property {() => readonly { readonly linearComparisonPlan?: { readonly edges?: readonly { readonly file?: unknown }[] } | null }[]} drawingList
 */


/**
 * A comparison file a Linear diagram uploaded for one edge, by ref.
 * @typedef {{ file?: InputFileRef | null }} InputComparisonFile
 */

/**
 * The input file fields of the project (`state.files`): each names its stored
 * files by `InputFileRef` of `state.inputFiles`.
 * @typedef {object} InputFiles
 * @property {InputFileRef | null} c_gb
 * @property {InputFileRef | null} c_gff
 * @property {InputFileRef | null} c_fasta
 * @property {InputFileRef | (InputFileRef | null)[] | (InputFileRef | null)[][] | null} c_depth one record's series, or a record-major matrix
 * @property {InputFileRef[]} c_conservation_blasts
 * @property {string | null} c_conservation_blasts_source a flag (`'losat-cache'`), not a file
 * @property {InputFileRef[]} c_conservation_fastas
 * @property {(InputFileRef | null)[]} c_conservation_sequence_sources
 * @property {InputComparisonFile[]} linearCanonicalComparisons
 * @property {InputFileRef | null} d_color
 * @property {InputFileRef | null} t_color
 * @property {InputFileRef | null} blacklist
 * @property {InputFileRef | null} whitelist
 * @property {InputFileRef | null} qualifier_priority
 */

/** The `InputFiles` fields that hold one ref or null. @typedef {'c_gb' | 'c_gff' | 'c_fasta' | 'd_color' | 't_color' | 'blacklist' | 'whitelist' | 'qualifier_priority'} InputFileSlot */

/**
 * Bind a file to one input field: the store registers it and the field holds
 * its ref. This is the one writer of a single-file field.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputFiles} files
 * @param {InputFileSlot} slot
 * @param {InputFileCandidate} file
 * @returns {InputFileRef | null}
 */
export const bindInputFile = (store, files, slot, file) => {
  const ref = store.registerFile(file);
  files[slot] = ref;
  return ref;
};

/** The `InputFiles` fields that hold a list of refs. @typedef {'c_conservation_blasts' | 'c_conservation_fastas'} InputFileListSlot */

/**
 * Bind the files of one list field; a place that holds no file is dropped.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputFiles} files
 * @param {InputFileListSlot} slot
 * @param {InputFileCandidate | readonly InputFileCandidate[]} list
 * @returns {InputFileRef[]}
 */
export const bindInputFileList = (store, files, slot, list) => {
  const refs = registerInputFiles(store, list).filter((ref) => ref !== null);
  files[slot] = refs;
  return refs;
};

/**
 * The refs of a file or a list of files, each place kept (null for no file).
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputFileCandidate | readonly InputFileCandidate[]} list
 * @returns {(InputFileRef | null)[]}
 */
export const registerInputFiles = (store, list) => (
  (Array.isArray(list) ? list : (list ? [list] : [])).map((file) => store.registerFile(file))
);

// The fields of `files` that bind input files; the other members are flags
// (`c_conservation_blasts_source`) or generated data.
const BOUND_FILE_KEYS = Object.freeze([
  'c_gb', 'c_gff', 'c_fasta', 'c_depth', 'c_conservation_blasts', 'c_conservation_fastas',
  'c_conservation_sequence_sources', 'd_color', 't_color', 'blacklist', 'whitelist', 'qualifier_priority'
]);
const LINEAR_SEQ_FILE_KEYS = Object.freeze(['gb', 'gff', 'fasta', 'depth']);

/**
 * The IDs of every file a binding names. A field that still holds a File is
 * registered, so its ID is the one History records.
 * @param {InputBindings} bindings
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @returns {Set<string>}
 */
export const boundInputFileIds = (bindings, store) => {
  /** @type {Set<string>} */
  const fileIds = new Set();
  /** @param {unknown} value */
  const add = (value) => {
    if (Array.isArray(value)) {
      value.forEach(add);
      return;
    }
    const ref = isInputFileRef(value) ? value : store.registerFile(/** @type {HistoryFile | null} */ (value));
    if (ref) fileIds.add(ref.fileId);
  };
  const files = bindings.files;
  BOUND_FILE_KEYS.forEach((key) => add(files[key]));
  (Array.isArray(files.linearCanonicalComparisons) ? files.linearCanonicalComparisons : [])
    .forEach((/** @type {{ file?: unknown } | null | undefined} */ comparison) => add(comparison?.file));
  bindings.linearSeqs.forEach((sequence) => LINEAR_SEQ_FILE_KEYS.forEach((key) => add(sequence[key])));
  bindings.drawingList().forEach((drawing) => {
    (drawing.linearComparisonPlan?.edges || []).forEach((edge) => add(edge.file));
  });
  return fileIds;
};

/**
 * The store's release for History's budget: it keeps the files History names
 * and every bound file.
 * @param {InputBindings} bindings
 * @param {Pick<InputFileStore, 'registerFile' | 'retainOnly'>} store
 * @returns {(historyFileIds: Iterable<string>) => void}
 */
export const retainHistoryAndBoundFiles = (bindings, store) => (historyFileIds) => {
  const retained = boundInputFileIds(bindings, store);
  for (const fileId of historyFileIds) retained.add(fileId);
  store.retainOnly(retained);
};
