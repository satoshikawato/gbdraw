// Run from the repository root. Synthetic inputs only; no application mutation.
import { normalizeUserFacingError } from '../../../../gbdraw/web/js/services/error-normalization.js';
import { deserializeWorkerError } from '../../../../gbdraw/web/js/services/diagram-generation.js';
import { runFeatureSearch } from '../../../../gbdraw/web/js/app/feature-search/search-core.js';

// The Worker module installs a message handler at import time. Do not dispatch it.
globalThis.self = {};
const { serializeError } = await import('../../../../gbdraw/web/js/workers/diagram-generation-worker.js');

const known = normalizeUserFacingError({
  type: 'ValueError',
  message: 'Comparison source feature index conflicts with its view feature ID.'
});
const tracebackOnly = normalizeUserFacingError({
  stderr: 'Traceback (most recent call last):\n  File "synthetic.py", line 1\nValueError: invalid comparison binding'
});
const structuredTransport = serializeError(Object.assign(new Error('Synthetic failure'), {
  code: 'COMPARISON_FEATURE_IDENTITY_CONFLICT',
  stage: 'render',
  context: { row: 2 }
}));
const transported = deserializeWorkerError(structuredTransport);
const patterns = ['(?i)NADH', '(?P<enzyme>NADH)', '(?<enzyme>NADH)', 'NADH'];
const search = patterns.map(query => {
  const result = runFeatureSearch({ features: [], query, useRegex: true });
  return { query, error: result.error };
});
console.log(JSON.stringify({
  known,
  tracebackOnly,
  structuredTransport: {
    keys: Object.keys(structuredTransport),
    code: structuredTransport.code ?? null,
    stage: structuredTransport.stage ?? null,
    context: structuredTransport.context ?? null
  },
  deserialized: { code: transported.code ?? null, stage: transported.stage ?? null },
  search
}, null, 2));
