// @ts-check
// The sequence hashes that key a Circular drawing's LOSAT conservation rings
// (DESIGN §1, K-3, D-15-123): one builder that Generate (to reuse or run a
// ring's search) and Save (to write into a drawing only the entries its rings
// name) both call. Every read is kept per File (`readOncePerFile`), so a File
// already hashed by a Generate or an earlier Save is not read again.
import { normalizeFileList, orderedConservationSources } from './conservation-series.js';
import { sha256Hex, textToBytes } from './byte-utils.js';
import { readAllSequenceRecords, readComparisonSequenceDigest, readOncePerFile } from './record-discovery.js';

/** @import { FileLike } from './session-resource-backing.js' */
/** @import { InputFileRef, InputFileStore } from './history-files.js' */
/** @import { CircularConservationLike, OrderedConservationSource } from './conservation-series.js' */

/**
 * What the key builder reads of a Circular drawing: its input type, the reference and ring files it binds, and its
 * conservation settings (the ring order and whether the rings read sequences or uploaded BLAST tables).
 * @typedef {{
 *   files: {
 *     readonly c_gb: InputFileRef | null, readonly c_fasta: InputFileRef | null,
 *     readonly c_conservation_fastas: (InputFileRef | null)[]
 *   },
 *   inputType: { readonly value: string },
 *   circularConservation: CircularConservationLike & { readonly source?: unknown }
 * }} ConservationDrawing
 */

/** @param {string} sequence */
const wrapFastaSequence = (sequence) => {
  const lines = [];
  for (let i = 0; i < sequence.length; i += 60) lines.push(sequence.slice(i, i + 60));
  return lines.join('\n');
};
/**
 * A record's LOSAT FASTA text, as the LOSAT search reads it.
 * @param {{ id: string, sequence: string }} record
 */
export const buildFastaText = (record) => `>${record.id}\n${wrapFastaSequence(record.sequence)}\n`;

/**
 * The conservation subject: the LOSAT FASTA of every record of the Circular reference (the records the CLI loader
 * reads, OV-1599; a record without a sequence is left out) and its SHA-256, the entry's `subjectCanonicalHash`.
 * @typedef {{ file: FileLike, fasta: string, hash: string }} ConservationSubject
 */
/**
 * One conservation ring in ring order: its source, its bound File and that File's query hash (the entry's
 * `queryCanonicalHash`, from the Python ring reader); `file` and the hash are null when the File is not bound.
 * @typedef {{ source: OrderedConservationSource<InputFileRef>, file: FileLike, queryCanonicalHash: string }
 *   | { source: OrderedConservationSource<InputFileRef>, file: null, queryCanonicalHash: null }} ConservationRingKey
 */
/**
 * @typedef {{ subject: ConservationSubject | null, rings: ConservationRingKey[] }} CircularConservationKeys
 */

/**
 * @param {FileLike} file
 * @param {'genbank' | 'fasta'} fmt
 * @returns {Promise<ConservationSubject>}
 */
const readConservationSubject = (file, fmt) => readOncePerFile([file], `ring-subject:${fmt}`, async () => {
  const records = (await readAllSequenceRecords(file, fmt))
    .flatMap(({ id, sequence }) => (sequence === null ? [] : [{ id, sequence }]));
  if (!records.length) throw new Error('No records found for circular conservation reference.');
  const fasta = records.map((record) => buildFastaText(record)).join('');
  return Object.freeze({ file, fasta, hash: await sha256Hex(textToBytes(fasta)) });
});

/**
 * The subject and the rings of a Circular drawing's sequence conservation (its LOSAT rings), from its bound files.
 * A drawing whose rings read uploaded BLAST tables has none. `beforeRing` runs before each ring is read (Generate's
 * cancel check).
 * @param {ConservationDrawing} drawing
 * @param {InputFileStore['fileOf']} fileOf
 * @param {() => void} [beforeRing]
 * @returns {Promise<CircularConservationKeys>}
 */
export const circularConservationKeys = async (drawing, fileOf, beforeRing = () => {}) => {
  if (String(drawing.circularConservation.source || '').trim().toLowerCase() === 'upload') return { subject: null, rings: [] };
  const sources = orderedConservationSources(normalizeFileList(drawing.files.c_conservation_fastas), drawing.circularConservation);
  if (!sources.length) return { subject: null, rings: [] };
  const genbank = drawing.inputType.value === 'gb';
  const subjectFile = fileOf(genbank ? drawing.files.c_gb : drawing.files.c_fasta);
  const subject = subjectFile ? await readConservationSubject(subjectFile, genbank ? 'genbank' : 'fasta') : null;
  /** @type {ConservationRingKey[]} */
  const rings = [];
  for (const source of sources) {
    beforeRing();
    const file = fileOf(source.file);
    rings.push(file
      ? { source, file, queryCanonicalHash: await readComparisonSequenceDigest({ file }) }
      : { source, file: null, queryCanonicalHash: null });
  }
  return { subject, rings };
};

/**
 * One ring's key in a set of pairs: its `subjectCanonicalHash` and `queryCanonicalHash`.
 * @param {unknown} subjectCanonicalHash
 * @param {unknown} queryCanonicalHash
 */
export const conservationPairKey = (subjectCanonicalHash, queryCanonicalHash) => JSON.stringify([subjectCanonicalHash, queryCanonicalHash]);

/**
 * The hash pairs (`conservationPairKey`) a drawing's bound conservation rings name: Save keeps a conservation entry
 * for the drawing only when its own two hashes are one of them. A ring or a subject that is not bound (after a
 * settings-only Load) names none.
 * @param {CircularConservationKeys} keys
 * @returns {Set<string>}
 */
export const conservationHashPairs = ({ subject, rings }) => new Set(subject
  ? rings.flatMap(({ queryCanonicalHash }) => (queryCanonicalHash === null ? [] : [conservationPairKey(subject.hash, queryCanonicalHash)]))
  : []);
