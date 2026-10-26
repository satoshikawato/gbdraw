// @ts-check
// The sequence hashes that key a Circular drawing's LOSAT conservation rings
// (DESIGN §1, K-3, D-15-123, D-15-124): one builder, which only Generate calls
// (to reuse or run a ring's search). Every read is kept per File
// (`readOncePerFile`), so a File a Generate already hashed is not read again.
// The hash pairs a request names (`ConservationPairs`) are recorded with the
// request when Generate or Load installs it; Save reads that record, no file.
import { normalizeFileList, orderedConservationSources } from './conservation-series.js';
import { sha256Hex, textToBytes } from './byte-utils.js';
import { CIRCULAR_CONSERVATION_FLOW } from './losat-cache.js';
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
 * The conservation hash pairs (`conservationPairKey`) a drawing's committed request names: runtime state of the
 * request, kept beside it (never written). Save writes the drawing a conservation entry only when the entry's own two
 * hashes are one of them.
 * @typedef {ReadonlySet<string>} ConservationPairs
 */
/**
 * The pairs Generate gives the request it installs (`adoptCanonicalRenderArtifacts`): its ring keys' pairs, or
 * `'committed'` when it re-renders the committed request without a conservation search, which keeps that request's.
 * @typedef {ConservationPairs | 'committed'} AdoptedConservationPairs
 */
/** A request that names no conservation ring, and a drawing without a request. @type {ConservationPairs} */
export const NO_CONSERVATION_PAIRS = Object.freeze(new Set());

/**
 * The pairs Generate records for the request it installs: those of the drawing's bound rings. A ring or a subject
 * that is not bound names none.
 * @param {CircularConservationKeys} keys
 * @returns {ConservationPairs}
 */
export const conservationHashPairs = ({ subject, rings }) => new Set(subject
  ? rings.flatMap(({ queryCanonicalHash }) => (queryCanonicalHash === null ? [] : [conservationPairKey(subject.hash, queryCanonicalHash)]))
  : []);

/**
 * The pairs Load records for a request it installs: those of the conservation entries the file wrote into the
 * request's drawing (unvalidated; an entry that is not a conservation entry names none).
 * @param {unknown} entries
 * @returns {ConservationPairs}
 */
export const conservationPairsOfEntries = (entries) => new Set((Array.isArray(entries) ? entries : [])
  .flatMap((/** @type {unknown} */ entry) => (
    entry !== null && typeof entry === 'object' && 'flow' in entry && entry.flow === CIRCULAR_CONSERVATION_FLOW
      && 'subjectCanonicalHash' in entry && typeof entry.subjectCanonicalHash === 'string'
      && 'queryCanonicalHash' in entry && typeof entry.queryCanonicalHash === 'string'
      ? [conservationPairKey(entry.subjectCanonicalHash, entry.queryCanonicalHash)] : []
  )));
