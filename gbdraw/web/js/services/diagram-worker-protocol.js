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
  EXTRACT_CDS_PROTEIN_FASTA: 'extractCdsProteinFasta',
  EXTRACT_FIRST_FASTA: 'extractFirstFasta',
  HYDRATE_PROTEIN_LOSAT_TSV: 'hydrateProteinLosatTsv',
  LIST_GFF_FASTA_RECORDS: 'listGffFastaRecords',
  LIST_SEQUENCE_RECORDS: 'listSequenceRecords',
  MEASURE_LEGEND_TEXT: 'measureLegendText',
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
