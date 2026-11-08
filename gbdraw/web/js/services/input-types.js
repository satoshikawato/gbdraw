// @ts-check
// The input types of a drawing's source files (GenBank, GFF3 + FASTA, FASTA):
// which file slots each type reads and which request record source kind it
// projects to. The one table of input types: readers look a type up here
// instead of branching on it. It imports only the generated contracts, so
// services and app modules read it without the inputs slice.
import { RECORD_SOURCE_FILES } from '../contracts/schema.generated.js';
/** @import { CanonicalRecordSource } from '../contracts/request.generated.js' */

/** @typedef {'gb' | 'gff' | 'fasta'} InputType */
/** @typedef {'gb' | 'gff' | 'fasta'} InputSlot A source file slot of a Linear card or of a Circular source card. */
/**
 * What an input type reads: the slot of its primary file, the slot of the FASTA a
 * GFF3 file pairs with (or null), and the record source kind the request projects
 * (`CanonicalRecordSource.kind`). The one table of input types: readers look a
 * type up here instead of branching on it.
 * @typedef {object} InputTypeSpec
 * @property {InputSlot} primary
 * @property {'fasta' | null} paired
 * @property {RecordSourceKind} sourceKind
 */
/** @typedef {CanonicalRecordSource['kind']} RecordSourceKind */
/** @type {Readonly<Record<InputType, Readonly<InputTypeSpec>>>} */
export const INPUT_TYPES = Object.freeze({
  gb: Object.freeze({ primary: 'gb', paired: null, sourceKind: 'genbank' }),
  gff: Object.freeze({ primary: 'gff', paired: 'fasta', sourceKind: 'gffFasta' }),
  fasta: Object.freeze({ primary: 'fasta', paired: null, sourceKind: 'fasta' })
});
/** @param {unknown} value A stored or chosen input type; anything else is GenBank. @returns {InputType} */
export const normalizeInputType = (value) => (value === 'gff' || value === 'fasta' ? value : 'gb');
/** @param {unknown} inputType */
export const inputTypeSpec = (inputType) => INPUT_TYPES[normalizeInputType(inputType)];
/** @param {unknown} sourceKind A record source kind. @returns {InputType} the input type that projects it (GenBank when none does) */
export const inputTypeOfSourceKind = (sourceKind) => normalizeInputType(
  /** @type {InputType[]} */ (Object.keys(INPUT_TYPES)).find((inputType) => INPUT_TYPES[inputType].sourceKind === sourceKind)
);
/**
 * The source files of one input: a Linear card's slots, or a Circular source card's
 * (`circularSourceInputs`). `complete` says the input has every file its type reads; `holdsFile`, that it has one of
 * them (a file in a slot its type does not read counts for neither).
 * @template F
 * @param {unknown} inputType
 * @param {Readonly<Partial<Record<InputSlot, F | null>>>} slots
 * @returns {{ spec: Readonly<InputTypeSpec>, primaryFile: F | null, pairedFile: F | null, complete: boolean, holdsFile: boolean }}
 */
export const inputSourceFiles = (inputType, slots) => {
  const spec = inputTypeSpec(inputType);
  const primaryFile = slots[spec.primary] ?? null;
  const pairedFile = spec.paired ? slots[spec.paired] ?? null : null;
  return {
    spec, primaryFile, pairedFile,
    complete: Boolean(primaryFile && (!spec.paired || pairedFile)),
    holdsFile: Boolean(primaryFile || pairedFile)
  };
};
/**
 * The file that holds an input's sequences (a GFF3 input's FASTA; otherwise its
 * primary file) and the format a sequence reader reads it in.
 * @param {unknown} inputType
 * @returns {{ spec: Readonly<InputTypeSpec>, slot: InputSlot, format: 'genbank' | 'fasta' }}
 */
export const inputSequenceSource = (inputType) => {
  const spec = inputTypeSpec(inputType);
  return { spec, slot: spec.paired ?? spec.primary, format: spec.sourceKind === 'genbank' ? 'genbank' : 'fasta' };
};
/**
 * One file of a record source: the resource kind that holds it and the source
 * field that names its resource id (`RECORD_SOURCE_FILES`, generated from
 * `_SOURCE_FILE_RESOURCES` in gbdraw/session_request_codec.py), with the input
 * slot that holds the file.
 * @typedef {(typeof RECORD_SOURCE_FILES)[RecordSourceKind][number] & { slot: InputSlot }} SourceSlotFile
 */
/** @param {Readonly<InputTypeSpec>} spec @returns {readonly SourceSlotFile[]} */
const sourceSlotFiles = ({ primary, paired, sourceKind }) => {
  // A record source names its files in this order: the primary file, then a GFF3 file's FASTA.
  const slots = paired ? [primary, paired] : [primary];
  /** @type {readonly (typeof RECORD_SOURCE_FILES)[RecordSourceKind][number][]} */
  const files = RECORD_SOURCE_FILES[sourceKind];
  if (files.length !== slots.length) throw new TypeError(`Record source kind '${sourceKind}' names ${files.length} files.`);
  return Object.freeze(files.map((file, index) => Object.freeze({ ...file, slot: slots[index] })));
};
/** The files of each input type's record source, in source order. @type {Readonly<Record<InputType, readonly SourceSlotFile[]>>} */
const INPUT_SOURCE_FILES = Object.freeze({
  gb: sourceSlotFiles(INPUT_TYPES.gb), gff: sourceSlotFiles(INPUT_TYPES.gff), fasta: sourceSlotFiles(INPUT_TYPES.fasta)
});
/**
 * The files of the record source an input type projects, in the order the source names them.
 * @param {unknown} inputType
 */
export const inputSourceSlotFiles = (inputType) => INPUT_SOURCE_FILES[normalizeInputType(inputType)];
/**
 * Whether a value is a record source: a kind of `RECORD_SOURCE_FILES` and a non-empty resource id in each
 * field the kind names, and no other field.
 * @param {Readonly<Record<string, unknown>>} value
 * @returns {value is CanonicalRecordSource}
 */
export const isRecordSource = (value) => {
  const kind = String(value.kind);
  if (!Object.hasOwn(RECORD_SOURCE_FILES, kind)) return false;
  /** @type {readonly { field: string }[]} */
  const files = RECORD_SOURCE_FILES[/** @type {RecordSourceKind} */ (kind)];
  return Object.keys(value).length === files.length + 1
    && files.every(({ field }) => typeof value[field] === 'string' && value[field] !== '');
};
/**
 * The record source of `sourceKind` whose files have the resource ids `idOf` gives (each file, in source order),
 * or null when a file has none. The one builder of a record source.
 * @param {RecordSourceKind} sourceKind
 * @param {(file: SourceSlotFile, index: number) => string | null | undefined} idOf
 * @returns {CanonicalRecordSource | null}
 */
export const recordSourceOf = (sourceKind, idOf) => {
  /** @type {Record<string, string>} */
  const source = { kind: sourceKind };
  for (const [index, file] of inputSourceSlotFiles(inputTypeOfSourceKind(sourceKind)).entries()) {
    const resourceId = idOf(file, index);
    if (!resourceId) return null;
    source[file.field] = resourceId;
  }
  return isRecordSource(source) ? source : null;
};
/**
 * The resource ids of a record source's files, in the order the source names them.
 * @param {CanonicalRecordSource} source
 * @returns {string[]}
 */
export const recordSourceResourceIds = (source) => {
  /** @type {Readonly<Record<string, unknown>>} */
  const fields = source;
  /** @type {readonly { field: string }[]} */
  const files = RECORD_SOURCE_FILES[source.kind];
  return files.map(({ field }) => String(fields[field]));
};
/**
 * The files of a record source kind (`RECORD_SOURCE_FILES`), in source order; none for an unknown kind.
 * @param {unknown} kind
 * @returns {readonly (typeof RECORD_SOURCE_FILES)[RecordSourceKind][number][]}
 */
export const recordSourceFilesOf = (kind) => (
  Object.hasOwn(RECORD_SOURCE_FILES, String(kind)) ? RECORD_SOURCE_FILES[/** @type {RecordSourceKind} */ (kind)] : []
);
/**
 * The resource id of each file of a stored record source, in the order its kind names them ('' for a missing one);
 * none for a source whose kind is unknown.
 * @param {Readonly<Record<string, unknown>> | null | undefined} source
 * @returns {string[]}
 */
export const storedRecordSourceResourceIds = (source) => (
  source ? recordSourceFilesOf(source.kind).map(({ field }) => String(source[field] || '').trim()) : []
);
/**
 * The resource id of each input slot's file in a stored record source, or '' for a slot its kind does not read
 * (and for every slot of a source whose kind is unknown).
 * @param {Readonly<Record<string, unknown>> | null | undefined} source
 * @returns {Record<InputSlot, string>}
 */
export const recordSourceSlotIds = (source) => {
  /** @type {Record<InputSlot, string>} */
  const ids = { gb: '', gff: '', fasta: '' };
  const resourceIds = storedRecordSourceResourceIds(source);
  if (resourceIds.length === 0) return ids;
  inputSourceSlotFiles(inputTypeOfSourceKind(source?.kind)).forEach(({ slot }, index) => { ids[slot] = resourceIds[index]; });
  return ids;
};
/**
 * The source files of one Circular source card (`CircularSource` of the inputs
 * slice): its uid and its files by the drawing's input type (`inputSourceFiles`).
 * @template F
 * @typedef {{ uid: string, spec: Readonly<InputTypeSpec>, primaryFile: F | null, pairedFile: F | null,
 *   complete: boolean, holdsFile: boolean }} CircularSourceInput
 */
/**
 * The source files of each Circular source card, in drawing order: the one
 * reader of a Circular drawing's input files by its input type.
 * @template F
 * @param {unknown} inputType
 * @param {readonly (Readonly<Partial<Record<InputSlot, F | null>>> & { readonly uid: string })[]} sources
 * @returns {CircularSourceInput<F>[]}
 */
export const circularSourceInputs = (inputType, sources) => sources.map((source) => ({
  uid: source.uid, ...inputSourceFiles(inputType, source)
}));
/**
 * Whether a Circular drawing's cards make an input Generate reads: at least one
 * card, and every card with each file its input type reads.
 * @param {readonly { complete: boolean }[]} inputs `circularSourceInputs`
 */
export const circularInputComplete = (inputs) => inputs.length > 0 && inputs.every((input) => input.complete);
