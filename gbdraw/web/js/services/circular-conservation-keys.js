// @ts-check
// The sequence hashes that key a Circular drawing's LOSAT conservation rings
// (DESIGN §1, K-3, D-15-123, D-15-124): one builder, which only Generate calls
// (to reuse or run a ring's search). Every read is kept per File
// (`readOncePerFile`), so a File a Generate already hashed is not read again.
// The hash pairs a request names (`ConservationPairs`) are recorded with the
// request when Generate or Load installs it; Save reads that record, no file.
import { inputSequenceSource } from './input-types.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { orderedConservationSources } from './conservation-series.js';
import { ringGenomeFileLists } from './drawing/inputs.js';
import { sha256Hex, textToBytes } from './byte-utils.js';
import { CIRCULAR_CONSERVATION_FLOW } from './losat-cache.js';
import { readAllSequenceRecords, readComparisonSequenceDigest, readOncePerFile } from './record-discovery.js';

/** @import { FileLike } from './session-resource-backing.js' */
/** @import { InputFileRef, InputFileStore } from './history-files.js' */
/** @import { CircularConservationLike, OrderedConservationSource } from './conservation-series.js' */
/** @import { CircularSource } from './drawing/inputs.js' */

/**
 * What the key builder reads of a Circular drawing: its input type, the reference and ring files it binds, and its
 * conservation settings (the ring order and whether the rings read sequences or uploaded BLAST tables).
 * @typedef {{
 *   files: { readonly c_conservation_fastas: readonly (readonly InputFileRef[])[] },
 *   circularSources: readonly Readonly<CircularSource>[],
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
 * The conservation subject: the LOSAT FASTA of every record of the Circular reference, each source card's file
 * whole, in card order (the records the CLI loader reads, OV-1599, as `ring_subject_records` in
 * gbdraw/comparisons/circular_losat.py takes them; a record without a sequence is left out) and its SHA-256, the
 * entry's `subjectCanonicalHash`, which the CLI's ring of the same files has. `format` is the files' format (a GFF3
 * card gives its FASTA file).
 * @typedef {{ files: readonly FileLike[], format: 'genbank' | 'fasta', fasta: string, hash: string }} ConservationSubject
 */
/**
 * One conservation ring in ring order: its source, the bound Files of its comparison genome (in the listed order) and
 * their query hash (the entry's `queryCanonicalHash`: of the genome's LOSAT query FASTA, from the Python ring reader,
 * as the CLI ring of the same files hashes it); `files` and the hash are null when a File is not bound.
 * @typedef {{ source: OrderedConservationSource<InputFileRef>, files: readonly [FileLike, ...FileLike[]], queryCanonicalHash: string }
 *   | { source: OrderedConservationSource<InputFileRef>, files: null, queryCanonicalHash: null }} ConservationRingKey
 */
/**
 * @typedef {{ subject: ConservationSubject | null, rings: ConservationRingKey[] }} CircularConservationKeys
 */

/**
 * One file's part of the subject: its records' IDs and LOSAT FASTA, with the FASTA's SHA-256.
 * @param {FileLike} file
 * @param {'genbank' | 'fasta'} fmt
 * @returns {Promise<{ ids: readonly string[], fasta: string, hash: string }>}
 */
const readSubjectFile = (file, fmt) => readOncePerFile([file], `ring-subject:${fmt}`, async () => {
  const records = (await readAllSequenceRecords(file, fmt))
    .flatMap(({ id, sequence }) => (sequence === null ? [] : [{ id, sequence }]));
  const fasta = records.map((record) => buildFastaText(record)).join('');
  return Object.freeze({ ids: Object.freeze(records.map(({ id }) => id)), fasta, hash: await sha256Hex(textToBytes(fasta)) });
});

/**
 * The subject of the source cards' files. Two cards that hold a record of one ID would put a ring's hits on either
 * record, so the subject refuses them before any search (Python refuses such a frame too, `_one_frame_per_id`).
 * @param {readonly FileLike[]} files one per card, in card order
 * @param {'genbank' | 'fasta'} fmt
 * @returns {Promise<ConservationSubject>}
 */
export const readConservationSubject = async (files, fmt) => {
  const parts = await Promise.all(files.map((file) => readSubjectFile(file, fmt)));
  /** @type {Map<string, number>} */
  const cardOfId = new Map();
  parts.forEach(({ ids }, card) => ids.forEach((recordId) => {
    const other = cardOfId.get(recordId);
    if (other !== undefined && other !== card) {
      throw diagnosticError('COMPARISON_INPUT', { reason: 'UNIQUE_IDS', field: 'records', recordId });
    }
    cardOfId.set(recordId, card);
  }));
  if (cardOfId.size === 0) throw new Error('No records found for circular conservation reference.');
  // One file's subject is its own part, whose hash is kept with it.
  const fasta = parts.map((part) => part.fasta).join('');
  const hash = parts.length === 1 ? parts[0].hash : await sha256Hex(textToBytes(fasta));
  return Object.freeze({ files: Object.freeze([...files]), format: fmt, fasta, hash });
};

/**
 * The rings a LOSAT Generate of the drawing searches, in ring order: each ring's genome files at its place, matched to
 * its series entry as Generate matches them (`orderedConservationSources`), with the label, colour and genetic code it
 * searches the ring with. A drawing whose rings read uploaded BLAST tables searches none.
 * @param {Pick<ConservationDrawing, 'files' | 'circularConservation'>} drawing
 * @returns {OrderedConservationSource<InputFileRef>[]}
 */
export const searchedRingSources = (drawing) => (
  String(drawing.circularConservation.source || '').trim().toLowerCase() === 'upload'
    ? [] : orderedConservationSources(ringGenomeFileLists(drawing.files.c_conservation_fastas), drawing.circularConservation)
);
/**
 * Generate's own rings, which its adoption of the result reads (`adoptCircularRingResult`): the rings it searched
 * (`searchedRingSources`, request ring k is `searched[k]`), and whether the draft's rings at the adoption are still
 * those (`builtFromDraft`; a change from outside the ring controls may have changed them while it ran).
 * @typedef {{ searched: readonly OrderedConservationSource<InputFileRef>[], builtFromDraft: boolean }} SearchedRingAdoption
 */
/**
 * One ring of a LOSAT Generate's request: its search result (`name`, `text`), the label and colour Generate searched
 * it with, its source index and its genome files as the request's serialized files hold them.
 * @template F
 * @typedef {{ name: string, text: string, sourceIndex: number, label: string, color: string, fastas: F[] }} LosatRequestRing
 */
/**
 * The request's rings of a LOSAT Generate, one per ring it searched, in ring order: ring k's table (`tables[k]`), and the
 * genome files at its place of the serialized ring genome lists (`ringFastas`), so the request names the files of the
 * rings it searched and no other.
 * @template F
 * @param {readonly OrderedConservationSource<InputFileRef>[]} searched
 * @param {readonly { name: string, text: string }[]} tables
 * @param {readonly (readonly F[])[]} ringFastas
 * @returns {LosatRequestRing<F>[]}
 */
export const losatRequestRings = (searched, tables, ringFastas) => searched.map((source, index) => ({
  name: tables[index].name,
  text: tables[index].text,
  sourceIndex: source.sourceIndex,
  label: source.label,
  color: source.color,
  fastas: [...(ringFastas[source.placeIndex] || [])]
}));

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
  const sources = searchedRingSources(drawing);
  if (!sources.length) return { subject: null, rings: [] };
  const { slot, format } = inputSequenceSource(drawing.inputType.value);
  const subjectFiles = drawing.circularSources.map((source) => fileOf(source[slot]));
  const subject = subjectFiles.length > 0 && subjectFiles.every(Boolean)
    ? await readConservationSubject(/** @type {FileLike[]} */ (subjectFiles), format) : null;
  /** @type {ConservationRingKey[]} */
  const rings = [];
  for (const source of sources) {
    beforeRing();
    const files = source.files.map((ref) => fileOf(ref));
    rings.push(files.every(Boolean)
      ? {
          source,
          files: /** @type {[FileLike, ...FileLike[]]} */ (files),
          queryCanonicalHash: await readComparisonSequenceDigest({ files: /** @type {[FileLike, ...FileLike[]]} */ (files) })
        }
      : { source, files: null, queryCanonicalHash: null });
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
