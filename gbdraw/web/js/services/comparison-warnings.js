// @ts-check
// Comparison table rows whose sequence IDs match no displayed record are drawn by
// position (PD-OI-074). The Worker carries the warning the CLI logs.
import { COMPARISON_WARNING_CODES, COMPARISON_WARNING_FIELDS } from '../contracts/schema.generated.js';
/** @import { SessionComparisonWarning } from '../contracts/session.generated.js' */

/** @type {ReadonlySet<unknown>} */
const CODES = new Set(COMPARISON_WARNING_CODES);

/**
 * @param {unknown} warnings
 * @param {{ readonly [index: number]: { readonly name?: unknown } | null | undefined } | null | undefined} results Read by index: an unvalidated value too.
 * @returns {SessionComparisonWarning[]}
 */
export const validateComparisonWarnings = (warnings, results) => {
  if (warnings === undefined) return [];
  const fields = COMPARISON_WARNING_FIELDS;
  if (!Array.isArray(warnings) || warnings.some((warning) => (
    (warning === null || typeof warning !== 'object' || Array.isArray(warning))
    || Object.keys(warning).length !== fields.length
    || fields.some((field) => !Object.prototype.hasOwnProperty.call(warning, field))
    || !CODES.has(warning.code)
    || ['queryRecordId', 'subjectRecordId', 'message', 'resultName']
      .some((field) => typeof warning[field] !== 'string')
    || !warning.message
    || !Array.isArray(warning.exampleIds)
    || warning.exampleIds.some((/** @type {unknown} */ id) => typeof id !== 'string')
    || ['queryRecordIndex', 'subjectRecordIndex', 'resultIndex']
      .some((field) => !Number.isSafeInteger(warning[field]) || warning[field] < 0)
    || !Number.isSafeInteger(warning.rowCount) || warning.rowCount < 1
    || results?.[warning.resultIndex]?.name !== warning.resultName
  ))) throw new Error('Comparison warnings do not match the successful Result metadata schema.');
  return warnings;
};
