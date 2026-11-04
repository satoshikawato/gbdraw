// @ts-check
// The input types of a drawing's source files (GenBank, GFF3 + FASTA, FASTA):
// which file slots each type reads and which request record source kind it
// projects to. The one table of input types: readers look a type up here
// instead of branching on it. No imports: services and app modules read it
// without the inputs slice.

/** @typedef {'gb' | 'gff' | 'fasta'} InputType */
/** @typedef {'gb' | 'gff' | 'fasta'} InputSlot A source file slot of a Linear card, or a Circular `c_` slot. */
/**
 * What an input type reads: the slot of its primary file, the slot of the FASTA a
 * GFF3 file pairs with (or null), and the record source kind the request projects
 * (`CanonicalRecordSource.kind`). The one table of input types: readers look a
 * type up here instead of branching on it.
 * @typedef {object} InputTypeSpec
 * @property {InputSlot} primary
 * @property {'fasta' | null} paired
 * @property {'genbank' | 'gffFasta' | 'fasta'} sourceKind
 */
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
 * The source files of one input: a Linear card's slots, or the Circular slots
 * (`circularInputSlots`). `complete` says the input has every file its type reads.
 * @template F
 * @param {unknown} inputType
 * @param {Readonly<Partial<Record<InputSlot, F | null>>>} slots
 * @returns {{ spec: Readonly<InputTypeSpec>, primaryFile: F | null, pairedFile: F | null, complete: boolean }}
 */
export const inputSourceFiles = (inputType, slots) => {
  const spec = inputTypeSpec(inputType);
  const primaryFile = slots[spec.primary] ?? null;
  const pairedFile = spec.paired ? slots[spec.paired] ?? null : null;
  return { spec, primaryFile, pairedFile, complete: Boolean(primaryFile && (!spec.paired || pairedFile)) };
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
 * The resource kind and record source field of each input slot's file, as the
 * Python codec names them (`_SOURCE_FILE_RESOURCES` in gbdraw/session_request_codec.py).
 * @type {Readonly<Record<InputSlot, readonly ['genbank' | 'gff3' | 'fasta', 'resourceId' | 'gffResourceId' | 'fastaResourceId']>>}
 */
export const INPUT_SLOT_RESOURCES = Object.freeze({
  gb: Object.freeze(/** @type {const} */ (['genbank', 'resourceId'])),
  gff: Object.freeze(/** @type {const} */ (['gff3', 'gffResourceId'])),
  fasta: Object.freeze(/** @type {const} */ (['fasta', 'fastaResourceId']))
});
/** The Circular file field of each source slot. @type {Readonly<Record<InputSlot, 'c_gb' | 'c_gff' | 'c_fasta'>>} */
export const CIRCULAR_INPUT_FIELDS = Object.freeze({ gb: 'c_gb', gff: 'c_gff', fasta: 'c_fasta' });
/**
 * The Circular source slots under the names a Linear card uses.
 * @template F
 * @param {{ readonly c_gb?: F | null, readonly c_gff?: F | null, readonly c_fasta?: F | null }} files
 * @returns {Record<InputSlot, F | null>}
 */
export const circularInputSlots = (files) => ({ gb: files.c_gb ?? null, gff: files.c_gff ?? null, fasta: files.c_fasta ?? null });
