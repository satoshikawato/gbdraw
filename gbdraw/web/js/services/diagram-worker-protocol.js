// @ts-check
export const DIAGRAM_HELPER_OPERATIONS = Object.freeze({
  EVALUATE_RULES: 'evaluateRules',
  EMBED_FONTS: 'embedFonts',
  PDF_FONTS: 'pdfFonts',
  REGISTER_FONTS: 'registerFonts',
  BUILD_PROTEIN_LOSAT_CACHE_KEYS: 'buildProteinLosatCacheKeys',
  COMPUTE_LABEL_FIT_ROW: 'computeLabelFitRow',
  CONVERT_LOSATP_PAIRS_TO_GENOMIC_PAYLOAD: 'convertLosatpPairsToGenomicPayload',
  CONVERT_MAIN_SESSION_COMPARISON_FRAME: 'convertMainSessionComparisonFrame',
  CONVERT_STORED_PROTEIN_RESULT_FRAMES: 'convertStoredProteinResultFrames',
  EXTRACT_CDS_PROTEIN_FASTA: 'extractCdsProteinFasta',
  EXTRACT_FIRST_FASTA: 'extractFirstFasta',
  HYDRATE_PROTEIN_LOSAT_TSV: 'hydrateProteinLosatTsv',
  LIST_GFF_FASTA_RECORDS: 'listGffFastaRecords',
  LIST_SEQUENCE_RECORDS: 'listSequenceRecords',
  PROMOTE_LEGACY_LOSATP_CACHE: 'promoteLegacyLosatpCache',
  READ_COMPARISON_SEQUENCE: 'readComparisonSequence',
  READ_FEATURE_OVERRIDE_TABLE: 'readFeatureOverrideTable',
  RESOLVE_LEGACY_PROTEIN_REFERENCES: 'resolveLegacyProteinReferences',
  RESOLVE_SIMILARITY_ALIGNMENT: 'resolveSimilarityAlignment',
  VALIDATE_CONFIG_OVERRIDES: 'validateConfigOverrides'
});

export const DIAGRAM_HELPER_OPERATION_NAMES = Object.freeze(
  Object.values(DIAGRAM_HELPER_OPERATIONS)
);

/**
 * One feature catalog build in one Worker lifetime: the cache epoch and the
 * build numbers of its items (gbdraw/web/CLAUDE.md, Runtime data flow). Opaque;
 * compared for equality only.
 * @typedef {string} FeatureCatalogRevisionId
 */

/**
 * A catalog revision with the digest of the catalog bytes the Worker last sent
 * for it, so an "unchanged" reply gets the same artifact fingerprint without
 * hashing again.
 * @typedef {object} FeatureCatalogRevision
 * @property {FeatureCatalogRevisionId} revision
 * @property {string} sha256 64 lowercase hex digits.
 * @property {number} bytes
 */

/**
 * Read a catalog revision from a message. `null` and `undefined` give `null`.
 * @param {unknown} value
 * @param {() => unknown} invalid Builds the caller's classified error for a malformed value.
 * @returns {Readonly<FeatureCatalogRevision> | null}
 */
export const readFeatureCatalogRevision = (value, invalid) => {
  if (value === null || value === undefined) return null;
  const record = /** @type {{ revision?: unknown, sha256?: unknown, bytes?: unknown }} */ (
    typeof value === 'object' && !Array.isArray(value) ? value : {}
  );
  const { revision, sha256, bytes } = record;
  if (
    typeof revision !== 'string' || !revision
    || typeof sha256 !== 'string' || !/^[0-9a-f]{64}$/.test(sha256)
    || typeof bytes !== 'number' || !Number.isSafeInteger(bytes) || bytes < 0
  ) throw invalid();
  return Object.freeze({ revision, sha256, bytes });
};

const SHA256_HEX = /^[0-9a-f]{64}$/;

/** @param {unknown} value */
const isByteCount = (value) => typeof value === 'number' && Number.isSafeInteger(value) && value >= 0;

/**
 * A children text the page holds, as a Generate request names it: its
 * revision (the sha256 of its UTF-8 bytes) and its UTF-8 length.
 * @typedef {object} WorkerHeldSegment
 * @property {string} revision 64 lowercase hex digits.
 * @property {number} bytes
 */

/**
 * Read the held segments of a run message: sorted by revision, each once.
 * `null` and `undefined` give an empty list.
 * @param {unknown} value
 * @param {() => unknown} invalid Builds the caller's classified error for a malformed value.
 * @returns {readonly Readonly<WorkerHeldSegment>[]}
 */
export const readHeldSegments = (value, invalid) => {
  if (value === null || value === undefined) return Object.freeze([]);
  if (!Array.isArray(value)) throw invalid();
  /** @type {Readonly<WorkerHeldSegment>[]} */
  const segments = [];
  for (const entry of value) {
    const record = /** @type {{ revision?: unknown, bytes?: unknown }} */ (
      entry && typeof entry === 'object' && !Array.isArray(entry) ? entry : {}
    );
    const { revision, bytes } = record;
    if (typeof revision !== 'string' || !SHA256_HEX.test(revision) || !isByteCount(bytes)) throw invalid();
    const previous = segments[segments.length - 1];
    if (previous && previous.revision >= revision) throw invalid();
    segments.push(Object.freeze({ revision, bytes: /** @type {number} */ (bytes) }));
  }
  return Object.freeze(segments);
};

/**
 * A placed content group of a Generate reply: its start tag, and its children
 * text named by `revision` and `bytes`. A sent text lies at `offset` in the
 * Result's `segments` buffer with `length` (equal to `bytes`); a held one has
 * neither. Its end tag is always `</g>`.
 * @typedef {object} SvgContentPiece
 * @property {'content'} kind
 * @property {string} start
 * @property {string} revision
 * @property {number} bytes
 * @property {number} [offset]
 * @property {number} [length]
 */

/**
 * Any other top-level element of a Generate reply, as text.
 * @typedef {{ kind: 'raw', text: string }} SvgRawPiece
 */

/**
 * One Result's SVG document as Python wrote it: the root start tag, the
 * top-level pieces in order, and the root end tag (gbdraw/web/CLAUDE.md,
 * Runtime data flow).
 * @typedef {object} SvgDocumentParts
 * @property {string} head
 * @property {'</svg>'} tail
 * @property {readonly (Readonly<SvgContentPiece> | Readonly<SvgRawPiece>)[]} pieces
 */

/**
 * Read the document of a Generate reply's Result. The sent children texts must
 * fill the Result's `segments` buffer of `segmentBytes` bytes, in order.
 * @param {unknown} value
 * @param {number} segmentBytes The byte length of the Result's `segments`, 0 without one.
 * @param {() => unknown} invalid Builds the caller's classified error for a malformed value.
 * @returns {Readonly<SvgDocumentParts>}
 */
export const readSvgDocumentParts = (value, segmentBytes, invalid) => {
  const record = /** @type {{ head?: unknown, tail?: unknown, pieces?: unknown }} */ (
    value && typeof value === 'object' && !Array.isArray(value) ? value : {}
  );
  const { head, tail, pieces } = record;
  if (typeof head !== 'string' || !head.startsWith('<svg') || tail !== '</svg>' || !Array.isArray(pieces)) {
    throw invalid();
  }
  let offset = 0;
  /** @type {(Readonly<SvgContentPiece> | Readonly<SvgRawPiece>)[]} */
  const read = [];
  for (const entry of pieces) {
    const piece = /** @type {Record<string, unknown>} */ (
      entry && typeof entry === 'object' && !Array.isArray(entry) ? entry : {}
    );
    if (piece.kind === 'raw' && typeof piece.text === 'string') {
      read.push(Object.freeze({ kind: 'raw', text: piece.text }));
      continue;
    }
    const { start, revision, bytes } = piece;
    if (
      piece.kind !== 'content' || typeof start !== 'string'
      || typeof revision !== 'string' || !SHA256_HEX.test(revision) || !isByteCount(bytes)
    ) throw invalid();
    const pieceBytes = /** @type {number} */ (bytes);
    if (piece.offset === undefined && piece.length === undefined) {
      read.push(Object.freeze({ kind: 'content', start, revision, bytes: pieceBytes }));
      continue;
    }
    if (piece.offset !== offset || piece.length !== pieceBytes) throw invalid();
    read.push(Object.freeze({ kind: 'content', start, revision, bytes: pieceBytes, offset, length: pieceBytes }));
    offset += pieceBytes;
  }
  if (offset !== segmentBytes) throw invalid();
  return Object.freeze({ head, tail: '</svg>', pieces: Object.freeze(read) });
};
