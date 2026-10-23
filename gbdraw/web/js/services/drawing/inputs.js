// @ts-check
// The input bindings: which stored file each input field names. A field holds
// an `InputFileRef` of the project's file store (`state.inputFiles`), never a
// File. The writers here take the object they write (the project's bindings
// now; a drawing's inputs slice later), so the store registers every bound
// file through one owner. The bound IDs keep their files stored: the store
// releases a file only when no History step and no binding names it.

/** @import { HistoryFile, InputFileRef, InputFileStore } from '../history-files.js' */

/**
 * The bindings that name stored files: the input fields, the Linear record
 * rows, and each drawing's uploaded comparison files.
 * @typedef {object} InputBindings
 * @property {Readonly<Record<string, unknown>>} files
 * @property {readonly Readonly<Record<string, unknown>>[]} linearSeqs
 * @property {() => readonly { readonly linearComparisonPlan?: { readonly edges?: readonly { readonly file?: unknown }[] } | null }[]} drawingList
 */

/**
 * @param {unknown} value
 * @returns {value is InputFileRef}
 */
export const isInputFileRef = (value) => Boolean(value) && typeof value === 'object'
  && typeof /** @type {{ fileId?: unknown }} */ (value).fileId === 'string';

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
