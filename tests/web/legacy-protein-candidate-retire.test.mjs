// OV-3561: a legacy schema-2 LOSATP candidate that verifies for a search this Generate answers with another entry (a
// current one, or another candidate it promoted) leaves the Session with the Generate's commit, as on the command line
// (gbdraw/api/session_compat.py `_promote_legacy_candidates`). The Python owner decides which candidates those are
// (`promote_legacy_protein_raw_cache_entries`, its `duplicates`); these tests run the production source of
// run-analysis.js that asks it and retires them.
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { join } from 'node:path';
import test from 'node:test';

const runSource = await readFile(join(process.cwd(), 'gbdraw/web/js/app/run-analysis.js'), 'utf8');

/** The production source between two markers of run-analysis.js. */
const slice = (/** @type {string} */ start, /** @type {string} */ end) => {
  const from = runSource.indexOf(start);
  const to = runSource.indexOf(end, from);
  assert.ok(from > 0 && to > from, `${start} .. ${end}`);
  return runSource.slice(from, to);
};

const tryPromoteSource = slice(
  '        const tryPromoteLegacyProteinEntry = async ({',
  '        /** @param {LosatJobSpec} spec */'
);
const commitSource = slice(
  '          if (legacyPromotionTransaction.length',
  '          structuredLosatTelemetry = {'
);

const legacyEntry = (/** @type {string} */ name) => ({ schema: 2, kind: 'raw-losat', key: name, text: '' });
const pendingEnvelope = () => ({
  schema: 1,
  entries: ['web-copy', 'cli-copy'].map((name) => ({
    state: 'pending', originalEntry: legacyEntry(name), rejectionReason: null
  }))
});

/** `tryPromoteLegacyProteinEntry` with a helper that answers `reply`, and what it asked and recorded. */
const promoter = (/** @type {Record<string, unknown>} */ reply) => {
  /** @type {any[]} */
  const payloads = [];
  const context = {
    workingLegacyProteinRawCandidates: pendingEnvelope(),
    workingProteinIdentityManifest: { schema: 1 },
    DIAGRAM_HELPER_OPERATIONS: { PROMOTE_LEGACY_LOSATP_CACHE: 'promoteLegacyLosatpCache' },
    async runDiagramHelperOperation(/** @type {string} */ _operation, /** @type {any} */ payload) {
      payloads.push(payload);
      return { result: structuredClone(reply) };
    },
    cloneJsonData: (/** @type {unknown} */ value) => structuredClone(value),
    normalizeLosatArgs: (/** @type {unknown[]} */ args) => [...args],
    transitionLegacyProteinCandidate() { throw new Error('no candidate is rejected here'); },
    LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA: 2,
    cacheMap: new Map(),
    getCurrentRawLosatCacheEntry: (/** @type {Map<string, any>} */ map, /** @type {string} */ key) => (
      map.has(key) ? { cacheKey: key, text: map.get(key).text } : null
    ),
    legacyPromotionTransaction: /** @type {any[]} */ ([]),
    legacyRetiredCandidateIndexes: /** @type {number[]} */ ([])
  };
  const tryPromote = new Function('context', `let {${Object.keys(context).join(',')}} = context;
    ${tryPromoteSource}
    return tryPromoteLegacyProteinEntry;`)(context);
  return { tryPromote, context, payloads };
};

const metadata = { identityKind: 'protein', program: 'blastp', outfmt: '6', args: ['--max-hsps-per-subject', '1'] };
const entries = { queryEntry: { fasta: '>q\nM\n', proteinMap: {} }, subjectEntry: { fasta: '>s\nM\n', proteinMap: {} } };

test('a search with a current entry promotes nothing and retires the candidates that verify for it', async () => {
  const { tryPromote, context, payloads } = promoter({ status: 'no-match', rejections: [], duplicates: [0, 1] });
  const found = await tryPromote({ cacheKey: 'search-key', metadata, ...entries, answered: true });
  assert.equal(found, null);
  assert.deepEqual(payloads[0].expectedOptions.answeredKeys, ['search-key']);
  assert.equal(payloads[0].candidates.length, 2, 'every pending candidate is asked about');
  assert.deepEqual(context.legacyRetiredCandidateIndexes, [0, 1]);
  assert.deepEqual(context.legacyPromotionTransaction, []);
  assert.equal(context.cacheMap.size, 0);
});

test('a promotion stores one candidate and retires the other copies of its search', async () => {
  const { tryPromote, context, payloads } = promoter({
    status: 'promoted', candidateIndex: 0, rejections: [], duplicates: [1],
    entry: { key: 'search-key', filename: 'pair.tsv', display: true, schema: 4, kind: 'raw-losat' },
    text: 'q\ts\t90\n'
  });
  const found = await tryPromote({ cacheKey: 'search-key', metadata, ...entries, answered: false });
  assert.deepEqual(found, { cacheKey: 'search-key', text: 'q\ts\t90\n' });
  assert.deepEqual(payloads[0].expectedOptions.answeredKeys, []);
  assert.deepEqual(context.legacyPromotionTransaction.map(({ candidateIndex }) => candidateIndex), [0]);
  assert.deepEqual(context.legacyRetiredCandidateIndexes, [1]);
  assert.deepEqual(context.cacheMap.get('search-key'), {
    schema: 4, kind: 'raw-losat', text: 'q\ts\t90\n', migratedFromSchema: 2
  }, 'the stored entry has no key, file name or display flag');
});

test('a Generate that only retires candidates commits them with its Result', () => {
  /** @type {any[]} */
  const writes = [];
  const context = {
    legacyPromotionTransaction: [{ candidateIndex: 2 }],
    legacyRetiredCandidateIndexes: [0, 1],
    writtenLegacyProteinArtifacts(/** @type {any} */ write) {
      writes.push(write);
      return { candidates: 'written', evidence: 'kept' };
    },
    commitProteinMigration: null,
    legacyPromotionCommitted: false
  };
  const commit = (/** @type {typeof context} */ state) => new Function('context', `let {${Object.keys(context).join(',')}} = context;
    ${commitSource}
    return commitProteinMigration;`)(state);
  const migrate = commit(context);
  assert.equal(typeof migrate, 'function');
  migrate({ legacyProteinRawCandidates: 'candidates', legacyProteinDerivedEvidence: 'evidence' });
  assert.deepEqual(writes[0].promoted, [2, 0, 1]);
  const retiredOnly = commit({ ...context, legacyPromotionTransaction: [], legacyRetiredCandidateIndexes: [1] });
  assert.equal(typeof retiredOnly, 'function', 'a Generate that answers every search from the cache still retires');
  retiredOnly({ legacyProteinRawCandidates: 'candidates', legacyProteinDerivedEvidence: 'evidence' });
  assert.deepEqual(writes[1].promoted, [1]);
  assert.equal(commit({ ...context, legacyPromotionTransaction: [], legacyRetiredCandidateIndexes: [] }), null);
});
