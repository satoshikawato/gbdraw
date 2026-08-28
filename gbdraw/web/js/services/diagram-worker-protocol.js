export const DIAGRAM_HELPER_OPERATIONS = Object.freeze({
  ASSEMBLE_PROTEIN_ANALYSIS: 'assembleProteinAnalysis',
  CONVERT_LOSAT_NUCLEOTIDE_TO_DISPLAY_TSV: 'convertLosatNucleotideToDisplayTsv',
  EXTRACT_CDS_PROTEIN_FASTA: 'extractCdsProteinFasta',
  EXTRACT_FIRST_FASTA: 'extractFirstFasta',
  GENERATE_LEGEND_ENTRY_SVG: 'generateLegendEntrySvg',
  HYDRATE_PROTEIN_LOSAT_TSV: 'hydrateProteinLosatTsv',
  LIST_GFF_FASTA_RECORDS: 'listGffFastaRecords',
  LIST_SEQUENCE_RECORDS: 'listSequenceRecords',
  MEASURE_LEGEND_TEXT: 'measureLegendText',
  PLAN_PROTEIN_ANALYSIS: 'planProteinAnalysis',
  PROMOTE_LEGACY_LOSATP_CACHE: 'promoteLegacyLosatpCache',
  REGENERATE_DEFINITION_SVGS: 'regenerateDefinitionSvgs',
  RESOLVE_LEGACY_PROTEIN_REFERENCES: 'resolveLegacyProteinReferences'
});

export const DIAGRAM_HELPER_OPERATION_NAMES = Object.freeze(
  Object.values(DIAGRAM_HELPER_OPERATIONS)
);
