import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';

const runAnalysis = await readFile(
  new URL('../../gbdraw/web/js/app/run-analysis.js', import.meta.url),
  'utf8'
);
const pythonHelpers = await readFile(
  new URL('../../gbdraw/web/js/app/python-helpers.js', import.meta.url),
  'utf8'
);
const workerProtocol = await readFile(
  new URL('../../gbdraw/web/js/services/diagram-worker-protocol.js', import.meta.url),
  'utf8'
);
const sessionRequest = await readFile(
  new URL('../../gbdraw/web/js/services/session-request.js', import.meta.url),
  'utf8'
);

for (const retiredOwner of [
  'buildLosatDerivedPayloadCachePayload',
  'resolveProteinBlastpCandidateLimit',
  'BUILD_PROTEIN_LOSAT_CACHE_KEY',
  'CONVERT_LOSATP_PAIRS_TO_GENOMIC_PAYLOAD',
  '--max-hsps-per-subject',
  '--max-target-seqs'
]) {
  assert.equal(
    runAnalysis.includes(retiredOwner),
    false,
    `${retiredOwner} must not remain in the Web execution owner`
  );
}

assert.equal(
  pythonHelpers.includes('def convert_losatp_blastp_pairs_to_genomic_payload'),
  false,
  'the embedded Python duplicate derived-stage algorithm must be removed'
);
assert.match(runAnalysis, /DIAGRAM_HELPER_OPERATIONS\.PLAN_PROTEIN_ANALYSIS/);
assert.match(runAnalysis, /DIAGRAM_HELPER_OPERATIONS\.ASSEMBLE_PROTEIN_ANALYSIS/);
assert.match(runAnalysis, /error\.name = String\(failure\.type/);
assert.match(runAnalysis, /error\.status = failure\.status/);
assert.match(workerProtocol, /PLAN_PROTEIN_ANALYSIS: 'planProteinAnalysis'/);
assert.match(workerProtocol, /ASSEMBLE_PROTEIN_ANALYSIS: 'assembleProteinAnalysis'/);

for (const retiredNormalizer of [
  'normalizeCollinearAnchorMode',
  'normalizeCollinearSearchScope',
  'normalizeOrthogroupMembershipMode'
]) {
  assert.equal(
    sessionRequest.includes(retiredNormalizer),
    false,
    `${retiredNormalizer} must not alter canonical protein request fields`
  );
}

const resultCommit = runAnalysis.indexOf('results.value = candidateCommit.results;');
const cacheCommit = runAnalysis.lastIndexOf('commitLosatCaches?.();');
assert.ok(resultCommit >= 0);
assert.ok(cacheCommit > resultCommit, 'cache publication must follow Result admission');
