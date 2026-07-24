import assert from 'node:assert/strict';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-losat-cache-'));
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempRoot, 'losat-cache.js'),
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'losat-cache.js'), 'utf8'),
  'utf8'
);

const cache = await import(pathToFileURL(join(tempRoot, 'losat-cache.js')));
assert.equal(cache.SESSION_LOSAT_CACHE_BYTE_LIMIT, 109_051_904);

const manifest = {
  schema: 1,
  proteinSets: {
    'sha256:set-a': { schema: 1, proteins: [{ featureAnalysisId: 'f_a', aaSha256: 'a'.repeat(64) }] },
    'sha256:set-b': { schema: 1, proteins: [{ featureAnalysisId: 'f_b', aaSha256: 'b'.repeat(64) }] }
  },
  recordAnalyses: {
    'sha256:analysis-a': { schema: 1, recordSourceId: 'A', proteinSetHash: 'sha256:set-a' },
    'sha256:analysis-b': { schema: 1, recordSourceId: 'B', proteinSetHash: 'sha256:set-b' }
  },
  recordInstances: {
    'record-1': {
      schema: 1,
      recordAnalysisId: 'sha256:analysis-a',
      bindingHash: 'sha256:binding-a',
      transportIds: { f_a: 'A@record-1|protein-a~f_a' }
    },
    'record-2': {
      schema: 1,
      recordAnalysisId: 'sha256:analysis-b',
      bindingHash: 'sha256:binding-b',
      transportIds: { f_b: 'B@record-2|protein-b~f_b' }
    }
  }
};

const proteinEntry = {
  schema: 3,
  kind: 'raw-losat',
  identityKind: 'protein',
  key: 'protein-key',
  text: 'A@record-1|protein-a~f_a\tB@record-2|protein-b~f_b\t100\t1\t0\t0\t1\t1\t1\t1\t0\t50\n',
  program: 'blastp',
  outfmt: '6',
  args: [],
  queryProteinSetHash: 'sha256:set-a',
  subjectProteinSetHash: 'sha256:set-b',
  queryBindingHash: 'sha256:binding-a',
  subjectBindingHash: 'sha256:binding-b',
  queryRecordInstanceKey: 'record-1',
  subjectRecordInstanceKey: 'record-2'
};
const nucleotideEntry = {
  schema: 2,
  kind: 'raw-losat',
  identityKind: 'nucleotide',
  key: 'nucleotide-key',
  text: '',
  program: 'blastn',
  outfmt: '6',
  args: [],
  queryCanonicalHash: 'q',
  subjectCanonicalHash: 's'
};
const legacyProteinEntry = {
  schema: 2,
  kind: 'raw-losat',
  key: 'legacy-key',
  text: 'p_r_old_0_3_1_deadbeefdead\tp_r_other_0_3_1_deadbeefdead\n',
  program: 'blastp',
  outfmt: '6',
  args: [],
  queryCanonicalHash: 'old-q',
  subjectCanonicalHash: 'old-s'
};

assert.equal(cache.classifyRawLosatCacheEntry(proteinEntry), 'protein-current');
assert.equal(cache.classifyRawLosatCacheEntry(nucleotideEntry), 'nucleotide-current');
assert.equal(cache.classifyRawLosatCacheEntry(legacyProteinEntry), 'protein-legacy');
assert.equal(cache.validateProteinIdentityManifest(manifest), true);
assert.equal(cache.validateProteinRawEntryReferences(proteinEntry, manifest), true);
const queryTransportIds = new Set(['A@record-1|protein-a~f_a']);
const subjectTransportIds = new Set(['B@record-2|protein-b~f_b']);
assert.equal(
  cache.rawProteinTextMatchesBindings(proteinEntry.text, queryTransportIds, subjectTransportIds),
  true
);
assert.equal(
  cache.rawProteinTextMatchesBindings(
    'A@record-1|protein-a~f_a\tB@record-2|protein-b~f_b\n',
    queryTransportIds,
    subjectTransportIds
  ),
  false,
  'binding IDs alone are not a valid LOSAT outfmt 6 row'
);
assert.equal(
  cache.rawProteinTextMatchesBindings(
    `${proteinEntry.text.trimEnd()}\textra\n`,
    queryTransportIds,
    subjectTransportIds
  ),
  false,
  'LOSAT outfmt 6 rows must not contain unexpected columns'
);

const manifestWithQueryTransportIds = (transportIds) => ({
  ...manifest,
  recordInstances: {
    ...manifest.recordInstances,
    'record-1': {
      ...manifest.recordInstances['record-1'],
      transportIds
    }
  }
});
const mismatchedBindingManifest = manifestWithQueryTransportIds({
  f_missing: 'A@record-1|protein-a~f_missing'
});
assert.equal(
  cache.validateProteinIdentityManifest(mismatchedBindingManifest),
  false,
  'a binding must reference exactly the features in its protein set'
);
assert.equal(
  cache.validateProteinRawEntryReferences(proteinEntry, mismatchedBindingManifest),
  false
);
assert.equal(
  cache.validateProteinIdentityManifest(manifestWithQueryTransportIds({
    f_a: 'A@record-1|protein-a~f_missing'
  })),
  false,
  'a transport ID must identify the feature referenced by its binding key'
);

const rawMap = new Map([['protein-key', { ...proteinEntry, key: undefined }]]);
assert.equal(
  cache.getCurrentRawLosatCacheEntry(
    rawMap,
    'protein-key',
    { program: 'blastp', outfmt: '6', args: [], queryBindingHash: 'sha256:binding-a', subjectBindingHash: 'sha256:binding-b' },
    manifest
  )?.entry?.text,
  proteinEntry.text
);
assert.equal(
  cache.getCurrentRawLosatCacheEntry(
    rawMap,
    'protein-key',
    { program: 'blastp', outfmt: '6', args: [], queryBindingHash: 'sha256:binding-b', subjectBindingHash: 'sha256:binding-a' },
    manifest
  ),
  null,
  'reverse query/subject bindings must not be a direct hit'
);

const invalidTextMap = new Map([['protein-key', {
  ...proteinEntry,
  key: undefined,
  text: 'unknown\tB@record-2|protein-b~f_b\t100\t1\t0\t0\t1\t1\t1\t1\t0\t50\n'
}]]);
assert.equal(
  cache.getCurrentRawLosatCacheEntry(invalidTextMap, 'protein-key', {
    program: 'blastp', outfmt: '6', args: [],
    queryBindingHash: 'sha256:binding-a', subjectBindingHash: 'sha256:binding-b'
  }, manifest),
  null
);

const pending = cache.createLegacyProteinCandidateEnvelope([legacyProteinEntry, nucleotideEntry]);
assert.equal(pending.entries.length, 1);
assert.equal(pending.entries[0].state, 'pending');
const savedBeforeGenerate = cache.serializableLegacyProteinCandidateEnvelope(pending);
assert.deepEqual(savedBeforeGenerate, pending);
assert.notStrictEqual(savedBeforeGenerate.entries[0].originalEntry, legacyProteinEntry);

const rejected = cache.transitionLegacyProteinCandidate(pending, 0, 'rejected', 'ambiguous mapping');
assert.equal(rejected.entries[0].state, 'rejected');
assert.equal(rejected.entries[0].rejectionReason, 'ambiguous mapping');
assert.equal(pending.entries[0].state, 'pending', 'candidate transition must be copy-on-write');
const promoted = cache.transitionLegacyProteinCandidate(pending, 0, 'promoted');
assert.equal(cache.serializableLegacyProteinCandidateEnvelope(promoted).entries.length, 0);

assert.equal(cache.isLosatDerivedCacheEntry({
  schema: 1, kind: 'derived-losatp-payload', key: 'legacy-derived', payload: {}
}), true);
assert.equal(cache.isLosatDerivedCacheEntry({
  schema: 2, kind: 'derived-losatp-payload', key: 'current-derived', payload: {}
}, { allowLegacy: false }), true);

const merged = cache.mergeProteinIdentityManifests([
  { ...manifest, recordInstances: { 'record-1': manifest.recordInstances['record-1'] } },
  { ...manifest, recordInstances: { 'record-2': manifest.recordInstances['record-2'] } }
]);
assert.deepEqual(Object.keys(merged.recordInstances).sort(), ['record-1', 'record-2']);
assert.throws(
  () => cache.mergeProteinIdentityManifests([
    manifest,
    {
      ...manifest,
      recordInstances: {
        ...manifest.recordInstances,
        'record-1': { ...manifest.recordInstances['record-1'], bindingHash: 'different' }
      }
    }
  ]),
  /conflicting record instance/
);

const rawCacheEntries = [
  { key: 'raw-first-visible', display: true, text: 'a'.repeat(40) },
  { key: 'raw-middle-visible', display: true, text: 'b'.repeat(40) },
  { key: 'raw-last-hidden', display: false, text: 'c'.repeat(40) }
];
const derivedCacheEntries = [
  { key: 'derived-first', payload: 'd'.repeat(40) },
  { key: 'derived-last', payload: 'e'.repeat(40) }
];
const exactBudget = cache.serializedLosatArtifactByteLength({
  rawEntries: [rawCacheEntries[2]],
  derivedEntries: [derivedCacheEntries[1]]
});
const prunedArtifacts = cache.pruneSerializedLosatArtifacts({
  rawEntries: rawCacheEntries,
  derivedEntries: derivedCacheEntries,
  maxBytes: exactBudget
});
assert.deepEqual(
  prunedArtifacts.rawEntries.map((entry) => entry.key),
  ['raw-last-hidden']
);
assert.deepEqual(
  prunedArtifacts.derivedEntries.map((entry) => entry.key),
  ['derived-last']
);
assert.equal(
  cache.serializedLosatArtifactByteLength(prunedArtifacts),
  exactBudget
);
