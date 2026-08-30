export const DIAGRAM_HELPER_OPERATIONS = Object.freeze({
  BUILD_PROTEIN_LOSAT_CACHE_KEY: 'buildProteinLosatCacheKey',
  CONVERT_LOSAT_NUCLEOTIDE_TO_DISPLAY_TSV: 'convertLosatNucleotideToDisplayTsv',
  CONVERT_LOSATP_PAIRS_TO_GENOMIC_PAYLOAD: 'convertLosatpPairsToGenomicPayload',
  EXTRACT_CDS_PROTEIN_FASTA: 'extractCdsProteinFasta',
  EXTRACT_FIRST_FASTA: 'extractFirstFasta',
  GENERATE_LEGEND_ENTRY_SVG: 'generateLegendEntrySvg',
  HYDRATE_PROTEIN_LOSAT_TSV: 'hydrateProteinLosatTsv',
  LIST_GFF_FASTA_RECORDS: 'listGffFastaRecords',
  LIST_SEQUENCE_RECORDS: 'listSequenceRecords',
  MEASURE_LEGEND_TEXT: 'measureLegendText',
  PROMOTE_LEGACY_LOSATP_CACHE: 'promoteLegacyLosatpCache',
  REGENERATE_DEFINITION_SVGS: 'regenerateDefinitionSvgs',
  RESOLVE_LEGACY_PROTEIN_REFERENCES: 'resolveLegacyProteinReferences',
  VALIDATE_CONFIG_OVERRIDES: 'validateConfigOverrides'
});

export const DIAGRAM_HELPER_OPERATION_NAMES = Object.freeze(
  Object.values(DIAGRAM_HELPER_OPERATIONS)
);
