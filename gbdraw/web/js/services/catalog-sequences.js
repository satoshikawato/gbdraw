// @ts-check
// Skeleton (S0 probe).
import { diagnosticError } from '../utils/error-normalization.js';

/** @param {unknown} result */
export const readFeatureSequencesReply = (result) => {
  if (!result) throw diagnosticError('RESULT_INVALID', {}, { stage: 'helper', operation: 'readFeatureSequences' });
  return result;
};
