const RESOURCE_BACKED_CANONICAL_COMPARISON_KINDS = new Set([
  'precomputedProteinComparison',
  'orthogroupResult',
  'collinearityResult'
]);

export const isResourceBackedCanonicalComparison = (comparison) => (
  RESOURCE_BACKED_CANONICAL_COMPARISON_KINDS.has(comparison?.kind)
);

export const canonicalComparisonResourceKind = (comparison) => ({
  precomputedProteinComparison: 'canonical-tsv',
  orthogroupResult: 'orthogroup-result',
  collinearityResult: 'collinearity-result'
}[comparison?.kind] || null);

export const mapResourceBackedCanonicalComparison = (
  comparison,
  mapFile = (file) => file
) => {
  if (!isResourceBackedCanonicalComparison(comparison)) return null;
  const kind = comparison.kind;
  return {
    kind,
    encoding: String(
      comparison.encoding || (
        kind === 'precomputedProteinComparison'
          ? 'canonicalTsv'
          : 'canonicalJson'
      )
    ),
    ...(kind === 'precomputedProteinComparison'
      ? {
          queryRecordIndex: Number(comparison.queryRecordIndex),
          subjectRecordIndex: Number(comparison.subjectRecordIndex)
        }
      : {}),
    ...(kind === 'collinearityResult'
      ? { valueKind: String(comparison.valueKind || 'result') }
      : {}),
    ...(comparison.canonicalInput === true ? { canonicalInput: true } : {}),
    file: mapFile(comparison.file)
  };
};
