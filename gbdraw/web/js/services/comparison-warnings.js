// Comparison table rows whose sequence IDs match no displayed record are drawn by
// position (PD-OI-074). The Worker carries the warning the CLI logs.
export const validateComparisonWarnings = (warnings, results) => {
  if (warnings === undefined) return [];
  const fields = ['code', 'queryRecordIndex', 'subjectRecordIndex', 'queryRecordId',
    'subjectRecordId', 'rowCount', 'exampleIds', 'message', 'resultIndex', 'resultName'];
  if (!Array.isArray(warnings) || warnings.some((warning) => (
    (warning === null || typeof warning !== 'object' || Array.isArray(warning))
    || Object.keys(warning).length !== fields.length
    || fields.some((field) => !Object.prototype.hasOwnProperty.call(warning, field))
    || warning.code !== 'comparison_record_id_unmatched'
    || ['queryRecordId', 'subjectRecordId', 'message', 'resultName']
      .some((field) => typeof warning[field] !== 'string')
    || !warning.message
    || !Array.isArray(warning.exampleIds)
    || warning.exampleIds.some((id) => typeof id !== 'string')
    || ['queryRecordIndex', 'subjectRecordIndex', 'resultIndex']
      .some((field) => !Number.isSafeInteger(warning[field]) || warning[field] < 0)
    || !Number.isSafeInteger(warning.rowCount) || warning.rowCount < 1
    || results?.[warning.resultIndex]?.name !== warning.resultName
  ))) throw new Error('Comparison warnings do not match the successful Result metadata schema.');
  return warnings;
};
