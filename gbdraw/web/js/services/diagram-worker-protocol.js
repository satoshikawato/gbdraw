// @ts-check
/** @import { buildPyodideAssetManifest } from './pyodide-assets.js' */
/** @import { FontIdentity } from './font-store.js' */
/** @import { CommittedCanonicalSession } from './session-request.js' */
/** @import { CanonicalRecord, CanonicalRenderRequest } from '../contracts/request.generated.js' */
/** @import { MATCH_SPAN_REASONS } from '../contracts/schema.generated.js' */

/**
 * A font file the page sends a Worker (init `fonts`, and each font it
 * registers later): the identity of its bytes, its file name or installed
 * face name, and the bytes. The Worker's font registry reads it.
 * @typedef {object} WorkerFontFile
 * @property {FontIdentity} identity
 * @property {string} label
 * @property {ArrayBuffer | ArrayBufferView} bytes
 */

/**
 * A file a helper request or a feature extraction gives the Worker by its
 * `role` or its `path` (an input before a request names it, or a table the
 * page wrote): its content as a Blob, which postMessage passes as a handle,
 * so no input's bytes cross it (OV-2749). The Worker reads it into its
 * request workspace.
 * @typedef {object} WorkerInputFile
 * @property {string} [role]
 * @property {string} [path]
 * @property {string} [name]
 * @property {Blob} blob
 */

/**
 * The page's first message to a new diagram Worker: its boot manifest, the
 * page's fonts, and, on a cross-origin-isolated page, the interrupt buffer
 * through which a Cancel stops Python (D7); `null` elsewhere, where a Cancel
 * terminates the Worker.
 * @typedef {ReturnType<typeof buildPyodideAssetManifest> & {
 *   type: 'init', id: string, fonts: WorkerFontFile[], interruptBuffer: Uint8Array | null
 * }} DiagramWorkerInitMessage
 */

/**
 * A Cancel while the Worker is kept: every request up to `through` (the last
 * request ID the page sent) is canceled; the Worker drops their replies,
 * clears the interrupt and the Python cancel flag after the operation in
 * flight, and answers with the same `id`.
 * @typedef {{ type: 'cancel', id: string, through: number }} DiagramWorkerCancelMessage
 */

/**
 * The Worker's answer to a Cancel: it is ready for the next request. `error`: the out-of-memory error of a reply it
 * dropped (OV-3027), so the page still replaces that engine.
 * @typedef {{ type: 'cancel', id: string, ok: true, error?: unknown }} DiagramWorkerCancelReply
 */

/**
 * Where a request's reply carries its error: `error` of a failed reply (`ok: false`, any request type), or
 * `results.error` of a Generate whose engine returned a structured error. A Cancel's drain reads a reply that was
 * on its way when the Cancel settled its request (OV-3027).
 * @typedef {object} DiagramWorkerReplyError
 * @property {unknown} [ok]
 * @property {unknown} [error]
 * @property {{ error?: unknown }} [results]
 */

/** The signal a Cancel writes into the interrupt buffer: SIGINT, as Pyodide's `setInterruptBuffer` reads it. */
export const DIAGRAM_WORKER_INTERRUPT_SIGNAL = 2;

export const DIAGRAM_HELPER_OPERATIONS = Object.freeze({
  EVALUATE_RULES: 'evaluateRules',
  EMBED_FONTS: 'embedFonts',
  PDF_FONTS: 'pdfFonts',
  REGISTER_FONTS: 'registerFonts',
  BUILD_PROTEIN_LOSAT_CACHE_KEYS: 'buildProteinLosatCacheKeys',
  CHECK_EXPORT_LIMIT: 'checkExportLimit',
  COMPUTE_LABEL_FIT_ROW: 'computeLabelFitRow',
  CONVERT_LOSATP_PAIRS_TO_GENOMIC_PAYLOAD: 'convertLosatpPairsToGenomicPayload',
  CONVERT_RING_TABLES_TO_SEARCH_FRAME: 'convertRingTablesToSearchFrame',
  CONVERT_STORED_PROTEIN_RESULT_FRAMES: 'convertStoredProteinResultFrames',
  DECIDE_AUTO_ORIENTATION: 'decideAutoOrientation',
  DERIVE_REGION_DRAWING: 'deriveRegionDrawing',
  EXTRACT_CDS_PROTEIN_FASTA: 'extractCdsProteinFasta',
  EXTRACT_FIRST_FASTA: 'extractFirstFasta',
  FIND_IN_OTHER_FILES: 'findInOtherFiles',
  HYDRATE_PROTEIN_LOSAT_TSV: 'hydrateProteinLosatTsv',
  LIST_GFF_FASTA_RECORDS: 'listGffFastaRecords',
  LIST_SEQUENCE_RECORDS: 'listSequenceRecords',
  PLAN_CARD_WINDOWS: 'planCardWindows',
  PLAN_GENOME_SCALES: 'planGenomeScales',
  PLAN_LINEAR_UNITS: 'planLinearUnits',
  PLAN_RECORD_PLACES: 'planRecordPlaces',
  PLAN_RING_SCALES: 'planRingScales',
  PROMOTE_LEGACY_LOSATP_CACHE: 'promoteLegacyLosatpCache',
  READ_ALL_SEQUENCE_RECORDS: 'readAllSequenceRecords',
  READ_COMPARISON_SEQUENCE: 'readComparisonSequence',
  READ_EMBED_SEQUENCES: 'readEmbedSequences',
  READ_FEATURE_OVERRIDE_TABLE: 'readFeatureOverrideTable',
  READ_FEATURE_SEQUENCES: 'readFeatureSequences',
  READ_MATCH_SEQUENCES: 'readMatchSequences',
  READ_RECORD_SPAN_SEQUENCES: 'readRecordSpanSequences',
  READ_SEARCH_QUERIES: 'readSearchQueries',
  RESOLVE_FEATURE_ANCHOR: 'resolveFeatureAnchor',
  RESOLVE_LEGACY_PROTEIN_REFERENCES: 'resolveLegacyProteinReferences',
  RESOLVE_SIMILARITY_ALIGNMENT: 'resolveSimilarityAlignment',
  SEARCH_HIT_ANNOTATION_SET: 'searchHitAnnotationSet',
  STORE_SAVED_PROTEIN_RESULTS: 'storeSavedProteinResults',
  SYNTHESIS_QUALIFIERS: 'synthesisQualifiers',
  SYNTHESIZE_COLOR_RULES: 'synthesizeColorRules',
  UPGRADE_LEGACY_SESSION: 'upgradeLegacySession',
  VALIDATE_CONFIG_OVERRIDES: 'validateConfigOverrides'
});

export const DIAGRAM_HELPER_OPERATION_NAMES = Object.freeze(
  Object.values(DIAGRAM_HELPER_OPERATIONS)
);

/**
 * The payload fields a staging helper reads its resource references from.
 * @typedef {object} StagingHelperPayload
 * @property {unknown} [canonicalRequest]
 * @property {{ canonicalRequest?: unknown } | null} [projection]
 * @property {unknown} [recordSource]
 * @property {readonly CanonicalRecord['source'][]} [referenceSources]
 */

/**
 * The helpers that read records of the drawing, each with the references it
 * names: the page stages those resources as a Generate does and adds the
 * staged payload keys; the Worker accepts those keys from these helpers only.
 * A helper's references come from the committed request, from the request of
 * a projection, from the one record source it reads (blastp extraction, OV-1555), or from the
 * reference's record sources (the ring scales, which the next Generate's render reuses).
 * @type {Readonly<Record<string, (payload: StagingHelperPayload) => unknown>>}
 */
const STAGED_HELPER_REFERENCES = Object.freeze({
  [DIAGRAM_HELPER_OPERATIONS.READ_FEATURE_OVERRIDE_TABLE]: (payload) => payload.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.RESOLVE_FEATURE_ANCHOR]: (payload) => payload.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.CONVERT_STORED_PROTEIN_RESULT_FRAMES]: (payload) => payload.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.CONVERT_RING_TABLES_TO_SEARCH_FRAME]: (payload) => payload.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.DECIDE_AUTO_ORIENTATION]: (payload) => payload.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.DERIVE_REGION_DRAWING]: (payload) => payload.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.FIND_IN_OTHER_FILES]: (payload) => payload.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.READ_RECORD_SPAN_SEQUENCES]: (payload) => payload.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.READ_FEATURE_SEQUENCES]: (payload) => payload.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.READ_MATCH_SEQUENCES]: (payload) => payload.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.READ_EMBED_SEQUENCES]: (payload) => payload.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.RESOLVE_SIMILARITY_ALIGNMENT]: (payload) => payload.projection?.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.EXTRACT_CDS_PROTEIN_FASTA]: (payload) => payload.recordSource,
  [DIAGRAM_HELPER_OPERATIONS.PLAN_GENOME_SCALES]: (payload) => payload.canonicalRequest,
  [DIAGRAM_HELPER_OPERATIONS.PLAN_RING_SCALES]: (payload) => payload.referenceSources
});

/**
 * The payload of `decideAutoOrientation`: the committed Linear request of the
 * last Generate, whose resources the page stages (`resources`, its Session
 * resource table).
 * @typedef {object} AutoOrientationHelperPayload
 * @property {CanonicalRenderRequest} canonicalRequest
 * @property {CommittedCanonicalSession['resources']} resources
 */

/**
 * The payload of `planGenomeScales`: a Linear drawing's draft request (`projectDraftRequest`), whose resources the
 * page stages (`resources`, the request's resource table), and the LOSAT program of its searches (`losatn`,
 * `tlosatx` or `losatp`; a request holds no nucleotide search, which the page resolves before it renders).
 * @typedef {object} GenomeScalesHelperPayload
 * @property {CanonicalRenderRequest} canonicalRequest
 * @property {CommittedCanonicalSession['resources']} resources
 * @property {'losatn' | 'tlosatx' | 'losatp'} program
 */

/**
 * One record of the `decideAutoOrientation` reply (`RecordOrientation` in
 * gbdraw/web_support/auto_orientation.py): the record key the committed
 * request draws it under, and the orientation `--auto_orient` draws it with.
 * @typedef {object} RecordOrientationDecision
 * @property {string} recordKey
 * @property {boolean} reverseComplement
 */

/**
 * The `decideAutoOrientation` reply (`AutoOrientationReply`): every drawn record, in drawn order.
 * @typedef {object} AutoOrientationReply
 * @property {readonly RecordOrientationDecision[]} records
 */

/**
 * One region of `deriveRegionDrawing` (`RegionInput` in
 * gbdraw/web_support/region_drawing.py): the source drawing's record key and
 * source coordinates, 1-based and inclusive.
 * @typedef {{ recordKey: string, start: number, end: number } | { recordKey: string, featureId: string }} RegionDrawingRegion
 *   Coordinates (`CoordinateRegionInput`) or a feature by its biological feature ID (`FeatureRegionInput`).
 */

/**
 * The payload of `deriveRegionDrawing`: the source drawing's draft projection
 * (its canonical request as the next Generate would send it), whose resources
 * the page stages, the regions, and whether hand-set sizes return to Auto.
 * @typedef {object} RegionDrawingHelperPayload
 * @property {CanonicalRenderRequest} canonicalRequest
 * @property {CommittedCanonicalSession['resources']} resources
 * @property {readonly RegionDrawingRegion[]} regions
 * @property {boolean} adaptSizes
 */

/*
 * The derived drawing (`DerivedRegionDrawing`: its mode, its canonical request as a Session stores a committed one,
 * the source drawing's resource each resource ID reads, the sizes that returned to Auto, and what is not carried) is
 * declared once, in Python (gbdraw/web_support/region_drawing.py), and generated into contracts/locate.generated.js.
 */
/** @import { DerivedRegionDrawing } from '../contracts/locate.generated.js' */

/**
 * A refused region (`RefusedRegionDrawing`): a region of a linear record whose
 * start is after its end (`StartAfterEndReason`); on a circular record such a
 * region is one window across the origin.
 * @typedef {object} RefusedRegionDrawing
 * @property {'refused'} status
 * @property {'START_AFTER_END'} reason
 * @property {string} recordKey
 * @property {number} start
 * @property {number} end
 */

/** The `deriveRegionDrawing` reply (`RegionDrawingReply`). @typedef {DerivedRegionDrawing | RefusedRegionDrawing} RegionDrawingReply */

/**
 * One span of `readRecordSpanSequences` (`RecordSpanInput` in gbdraw/web_support/record_spans.py): a record key of the
 * request and source coordinates, 1-based and inclusive; a start after the end crosses the origin.
 * @typedef {object} RecordSpanInput
 * @property {string} recordKey
 * @property {number} start
 * @property {number} end
 */

/**
 * A record the request draws (`DrawnRecordInput`): Python reads the bases it draws and, when a span of it may cross
 * its origin, its origin window (its last and first bases, as many as `originWindow`, a sequence search's longest query,
 * and at most half the record each); a record asked twice is read once, and the windows follow every drawn span.
 * @typedef {object} DrawnRecordInput
 * @property {string} recordKey
 * @property {number} originWindow
 */

/**
 * What to read: a span by its coordinates, or a drawn record; one read asks for one kind.
 * @typedef {RecordSpanInput | DrawnRecordInput} RecordSpanRead
 */

/**
 * The payload of `readRecordSpanSequences`: the request the shown Result was made from (the committed request),
 * whose resources the page stages, and the spans to read.
 * @typedef {object} RecordSpanSequencesPayload
 * @property {CanonicalRenderRequest} canonicalRequest
 * @property {CommittedCanonicalSession['resources']} resources
 * @property {readonly RecordSpanRead[]} spans
 */

/**
 * One span's bases on the source strand (`RecordSpanSequence`), with the FASTA header that names them (unique in one
 * read), the length of its record, and whether it was asked as an origin window.
 * @typedef {object} RecordSpanSequence
 * @property {string} recordKey
 * @property {string} recordId
 * @property {number} start
 * @property {number} end
 * @property {number} length
 * @property {boolean} wrapsOrigin
 * @property {boolean} originWindow
 * @property {number} recordLength
 * @property {string} header
 * @property {string} sequence
 */

/**
 * The `readRecordSpanSequences` reply (`RecordSpanSequencesReply`): one span per span asked by coordinates, in order;
 * or, for drawn records, one span per record in the order first asked, then the origin windows; and their FASTA text.
 * @typedef {object} RecordSpanSequencesReply
 * @property {readonly RecordSpanSequence[]} spans
 * @property {string} fasta
 */

// The sequences the feature catalog does not store (D5), read by Python from the committed request
// (gbdraw/web_support/catalog_sequences.py; each typedef below is the twin of the TypedDict it names).

/**
 * One feature by identity (`FeatureSequenceInput`): the record key and the biological feature ID of its catalog row.
 * @typedef {object} FeatureSequenceInput
 * @property {string} recordKey
 * @property {string} biologicalFeatureId
 */

/**
 * One feature's sequences with the identity asked (`FeatureSequencesOf`).
 * @typedef {object} FeatureSequencesOf
 * @property {string} recordKey
 * @property {string} biologicalFeatureId
 * @property {string} nucleotideSequence
 * @property {string} aminoAcidSequence
 * @property {readonly string[]} sequenceWarnings
 */

/**
 * The payload of `readFeatureSequences`: the committed request, its resources, and the features (each once).
 * @typedef {object} FeatureSequencesPayload
 * @property {CanonicalRenderRequest} canonicalRequest
 * @property {CommittedCanonicalSession['resources']} resources
 * @property {readonly FeatureSequenceInput[]} features
 */

/** The `readFeatureSequences` reply (`FeatureSequencesReply`): the features in the order asked. @typedef {{ features: readonly FeatureSequencesOf[] }} FeatureSequencesReply */

/**
 * A span coordinate (`SpanCoordinate`): a whole number, or null when the drawing's value is not a whole number
 * (Python says so after it finds the record, as the match popup did before D5).
 * @typedef {number | null} SpanCoordinate
 */

/**
 * A span of a drawn record by the record key of its drawn unit (`RecordMatchSpan`), with the record ID or alias the
 * match names ('' when none), which the unit must carry: displayed coordinates, 1-based and inclusive; a start after
 * the end is the reverse strand.
 * @typedef {object} RecordMatchSpan
 * @property {'record'} kind
 * @property {string} recordKey
 * @property {string} recordId
 * @property {SpanCoordinate} start
 * @property {SpanCoordinate} end
 */

/**
 * A span of a drawn record named by its record ID or an alias (`ReferenceMatchSpan`): a Circular comparison band
 * names its reference by the ID its hit carries, and Python finds the record among the drawn records by ID, then by
 * alias (name, accession), among every drawn unit or only the unit of `recordKey` (a band on a cropped or reversed
 * reference carries its record index), which it reads when the span names no ID; coordinates as in `RecordMatchSpan`.
 * @typedef {object} ReferenceMatchSpan
 * @property {'reference'} kind
 * @property {string} recordId
 * @property {string | null} recordKey
 * @property {SpanCoordinate} start
 * @property {SpanCoordinate} end
 */

/**
 * A span of a comparison FASTA record (`ComparisonMatchSpan`), by the index of its comparison source (null: the match
 * names none, and every source is searched) and its record ID or alias; coordinates as in `RecordMatchSpan`.
 * @typedef {object} ComparisonMatchSpan
 * @property {'comparison'} kind
 * @property {number | null} sourceIndex
 * @property {string} recordId
 * @property {SpanCoordinate} start
 * @property {SpanCoordinate} end
 */

/** @typedef {RecordMatchSpan | ReferenceMatchSpan | ComparisonMatchSpan} MatchSpanInput */
/** Why a span's bases cannot be read (`MatchSpanReason`, Python's list). @typedef {typeof MATCH_SPAN_REASONS[number]} MatchSpanReason */

/**
 * The bases of one span on its strand, or why they cannot be read (`MatchSpanSequence`), with the length of the
 * record the span reads (0 when no record is found).
 * @typedef {object} MatchSpanSequence
 * @property {'ok' | 'unavailable'} status
 * @property {MatchSpanReason | null} reason
 * @property {string} recordId
 * @property {'+' | '-'} orientation
 * @property {number} length
 * @property {string} sequence
 * @property {number} sourceLength
 */

/**
 * The payload of `readMatchSequences`: the committed request, its resources, and the spans.
 * @typedef {object} MatchSequencesPayload
 * @property {CanonicalRenderRequest} canonicalRequest
 * @property {CommittedCanonicalSession['resources']} resources
 * @property {readonly MatchSpanInput[]} spans
 */

/** The `readMatchSequences` reply (`MatchSequencesReply`): the spans in the order asked. @typedef {{ spans: readonly MatchSpanSequence[] }} MatchSequencesReply */

/**
 * The part of the stored catalog item an export read sends (`readEmbedSequences`): the item's drawn records and
 * matches, and the biological rows of the part.
 * @typedef {object} EmbedItemPart
 * @property {readonly string[]} recordKeys
 * @property {Readonly<Record<string, string>>} [windowOf]
 * @property {readonly unknown[]} [comparisonMatches]
 * @property {readonly Readonly<Record<string, unknown>>[]} biologicalFeatures
 */

/**
 * Which whole-record sources a part returns (`EmbedPartInput`): those of the records named by record key, and, with
 * `comparisonSources`, those of the comparison FASTA records.
 * @typedef {{ sourcesOf: readonly string[], comparisonSources: boolean }} EmbedPartInput
 */

/**
 * The file of an export read: the JSON of its `EmbedItemPart` as a Blob, like every helper file the Worker stages,
 * which the Worker stages as a file for Python to read (a long JavaScript string handed to Python takes many times
 * its size of Worker heap).
 * @typedef {{ role: 'embedItem', blob: Blob }} EmbedItemFile
 */

/**
 * The payload of `readEmbedSequences`: the committed request, its resources, and one part of the stored catalog
 * item of the Result to export (`files`); the reply is that part in the export encoding as its text
 * (`EmbedPartText`, `read_embed_sequences_json`).
 * @typedef {object} EmbedSequencesPayload
 * @property {CanonicalRenderRequest} canonicalRequest
 * @property {CommittedCanonicalSession['resources']} resources
 * @property {[EmbedItemFile]} files
 * @property {EmbedPartInput} part
 */

/**
 * The `readEmbedSequences` reply (`EmbedPartText`): the part's rows in the order sent, with their sequences, the
 * whole-record sources the part owns, each with its place among the item's `sourceCount` sources, and, in the part
 * that holds the comparison sources, what the item's match spans read (`matchSpanReads`, Python's
 * `embedded_match_span_reads`, read once per item) and the sources only a span reads, after the others (so its
 * `sourceCount` may be larger than the other parts'), all as the text the file embeds (`embed_json_text`), so neither
 * the Worker nor the page parses them.
 * @typedef {object} EmbedPartText
 * @property {string} biologicalFeatures The rows' text joined by commas ('' for none).
 * @property {number} rowCount
 * @property {readonly { index: number, source: string }[]} sequenceSources
 * @property {number} sourceCount
 * @property {string | null} matchSpanReads Null in the other parts and when no match names a span.
 */

/** @type {readonly string[]} */
export const DIAGRAM_HELPER_STAGING_OPERATIONS = Object.freeze(Object.keys(STAGED_HELPER_REFERENCES));

/**
 * The references a helper call stages, or `null` when it stages nothing.
 * @param {string} operation
 * @param {StagingHelperPayload} payload
 * @returns {unknown}
 */
export const stagedHelperReferences = (operation, payload) => (
  Object.hasOwn(STAGED_HELPER_REFERENCES, operation) ? STAGED_HELPER_REFERENCES[operation](payload) ?? null : null
);

/**
 * The payload keys a resource staging adds (`PreparedResources` without `commit`).
 * @type {readonly ['resourceManifest', 'stagedResources', 'retainedResourceTokens']}
 */
export const STAGED_RESOURCE_PAYLOAD_KEYS = Object.freeze(
  /** @type {const} */ (['resourceManifest', 'stagedResources', 'retainedResourceTokens'])
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
