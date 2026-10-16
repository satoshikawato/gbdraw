// @ts-check
/**
 * What these helpers read of a canonical comparison entry.
 * @typedef {{
 *   kind?: string,
 *   resourceId?: string,
 *   encoding?: string,
 *   queryRecordIndex?: number,
 *   subjectRecordIndex?: number,
 *   valueKind?: string,
 *   canonicalInput?: boolean,
 *   file?: unknown
 * }} ResourceBackedComparison
 */

/** @type {Set<string | undefined>} */
const RESOURCE_BACKED_CANONICAL_COMPARISON_KINDS = new Set([
  'precomputedProteinComparison',
  'orthogroupResult',
  'collinearityResult'
]);

/**
 * @template {ResourceBackedComparison | null | undefined} C
 * @param {C} comparison
 * @returns {comparison is C & ResourceBackedComparison & { kind: 'precomputedProteinComparison' | 'orthogroupResult' | 'collinearityResult' }}
 *   The kind is one of the set's; the other members are as `C` and the view type them.
 */
export const isResourceBackedCanonicalComparison = (comparison) => (
  RESOURCE_BACKED_CANONICAL_COMPARISON_KINDS.has(comparison?.kind)
);

/** @param {ResourceBackedComparison | null | undefined} comparison */
export const canonicalComparisonResourceKind = (comparison) => (/** @type {Record<string, string>} */ ({
  precomputedProteinComparison: 'canonical-tsv',
  orthogroupResult: 'orthogroup-result',
  collinearityResult: 'collinearity-result'
}))[/** @type {string} */ (comparison?.kind)] || null;

/**
 * @param {ResourceBackedComparison | null | undefined} comparison
 * @param {(file: unknown) => unknown} [mapFile]
 */
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
