// @ts-check
// The input bindings: which stored file each input field names. A field holds
// an `InputFileRef` of the project's file store (`state.inputFiles`), never a
// File. The writers here take the object they write (the project's bindings
// now; a drawing's inputs slice later), so the store registers every bound
// file through one owner. The bound IDs keep their files stored: the store
// releases a file only when no History step and no binding names it.

/** @import { InputFileCandidate, InputFileRef, InputFileStore } from '../history-files.js' */
/** @import { RequestCanonicalComparison } from '../session-request.js' */
/** @import { LinearComparisonEdge, LinearComparisonPlan } from '../linear-comparisons.js' */
/** @import { LinearSeq } from '../../state.js' */
/** @import { DepthFileSlots } from '../depth-track-state.js' */
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
 * A Linear comparison a Session committed, as the bindings hold it: the
 * request entry's fields, with a resource-backed entry's file by ref.
 * @typedef {Omit<RequestCanonicalComparison, 'file'> & { file?: InputFileRef | null }} BoundCanonicalComparison
 */
/**
 * A binding's value: one ref, or a list of them (a Depth matrix nests lists).
 * @typedef {InputFileRef | null | undefined | readonly BoundFileValue[]} BoundFileValue
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
 * @property {(InputFileRef | null)[]} c_conservation_fastas
 * @property {(InputFileRef | null)[]} c_conservation_sequence_sources
 * @property {BoundCanonicalComparison[]} linearCanonicalComparisons
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
 * The IDs of every file a binding names.
 * @param {InputBindings} bindings
 * @returns {Set<string>}
 */
export const boundInputFileIds = (bindings) => {
  /** @type {Set<string>} */
  const fileIds = new Set();
  /** @param {unknown} value */
  const add = (value) => {
    if (Array.isArray(value)) {
      value.forEach(add);
      return;
    }
    if (isInputFileRef(value)) fileIds.add(value.fileId);
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
 * @param {Pick<InputFileStore, 'retainOnly'>} store
 * @returns {(historyFileIds: Iterable<string>) => void}
 */
export const retainHistoryAndBoundFiles = (bindings, store) => (historyFileIds) => {
  const retained = boundInputFileIds(bindings);
  for (const fileId of historyFileIds) retained.add(fileId);
  store.retainOnly(retained);
};

/**
 * A loaded value of a binding: a File or Session resource view (a Load
 * candidate), a ref (a rolled-back snapshot), or a list of them.
 * @typedef {InputFileCandidate | InputFileRef | readonly LoadedFileValue[]} LoadedFileValue
 */

/**
 * The refs of a loaded value, each place and each nested list kept.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {LoadedFileValue} value
 * @returns {BoundFileValue}
 */
export const registerInputFileTree = (store, value) => (
  Array.isArray(value)
    ? value.map((entry) => registerInputFileTree(store, entry))
    : store.registerFile(/** @type {InputFileCandidate} */ (value))
);

/**
 * The comparisons a Session committed, as the bindings hold them: a
 * resource-backed entry's loaded file is registered and named by its ref.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {readonly (Omit<BoundCanonicalComparison, 'file'> & { file?: unknown })[]} comparisons
 * @returns {BoundCanonicalComparison[]}
 */
export const bindCanonicalComparisons = (store, comparisons) => comparisons.map((comparison) => (
  Object.hasOwn(comparison, 'file')
    ? { ...comparison, file: store.registerFile(/** @type {InputFileCandidate} */ (comparison.file)) }
    : /** @type {BoundCanonicalComparison} */ (comparison)
));

/**
 * A Linear record row as a Load candidate or a rolled-back snapshot holds it.
 * @typedef {Omit<LinearSeq, 'gb' | 'gff' | 'fasta' | 'depth'> & {
 *   gb: LoadedFileValue, gff: LoadedFileValue, fasta: LoadedFileValue, depth: LoadedFileValue
 * }} LoadedLinearSeq
 */
/**
 * A comparison plan as a Load candidate or a rolled-back snapshot holds it.
 * @typedef {Omit<LinearComparisonPlan, 'edges'> & {
 *   edges: (Omit<LinearComparisonEdge, 'file'> & { file: LoadedFileValue })[]
 * }} LoadedComparisonPlan
 */
/**
 * The project inputs a Load installs (its candidate, which holds the loaded
 * Files) or a rollback restores (refs).
 * @typedef {object} LoadedInputs
 * @property {Readonly<Record<string, unknown>>} files
 * @property {readonly LoadedLinearSeq[]} linearSeqs
 * @property {LoadedComparisonPlan} linearComparisonPlan
 * @property {{ primaryFile?: LoadedFileValue, pairedFile?: LoadedFileValue }} circularRecordDiscovery
 */
/**
 * The bindings a Load installs: every file it holds is registered here, at
 * install, and each field holds its ref (a ref of a rollback stays bound).
 * The Load candidate's own reads before the commit read its Files directly.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {LoadedInputs} loaded
 */
export const bindLoadedInputs = (store, loaded) => {
  /** @param {unknown} value */
  const tree = (value) => registerInputFileTree(store, /** @type {LoadedFileValue} */ (value));
  /** @param {unknown} value */
  const one = (value) => store.registerFile(/** @type {InputFileCandidate} */ (value));
  /** @type {Record<string, unknown>} */
  const files = { ...loaded.files };
  BOUND_FILE_KEYS.forEach((key) => {
    if (Object.hasOwn(files, key)) files[key] = tree(files[key]);
  });
  if (Array.isArray(files.linearCanonicalComparisons)) {
    files.linearCanonicalComparisons = bindCanonicalComparisons(store, files.linearCanonicalComparisons);
  }
  /** @type {LinearSeq[]} */
  const linearSeqs = loaded.linearSeqs.map((seq) => ({
    ...seq,
    gb: one(seq.gb),
    gff: one(seq.gff),
    fasta: one(seq.fasta),
    depth: /** @type {DepthFileSlots<InputFileRef>} */ (tree(seq.depth))
  }));
  /** @type {LinearComparisonPlan} */
  const linearComparisonPlan = {
    ...loaded.linearComparisonPlan,
    edges: loaded.linearComparisonPlan.edges.map((edge) => ({ ...edge, file: one(edge.file) }))
  };
  const discovery = loaded.circularRecordDiscovery;
  return {
    files,
    linearSeqs,
    linearComparisonPlan,
    circularRecordDiscovery: {
      ...discovery,
      ...(Object.hasOwn(discovery, 'primaryFile') ? { primaryFile: one(discovery.primaryFile) } : {}),
      ...(Object.hasOwn(discovery, 'pairedFile') ? { pairedFile: one(discovery.pairedFile) } : {})
    }
  };
};
