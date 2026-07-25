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
const numericContract = JSON.parse(await readFile(
  join(repoRoot, 'tests', 'fixtures', 'losat-outfmt6-numeric-contract.json'),
  'utf8'
));

const featureA = `f_${'a'.repeat(64)}`;
const featureB = `f_${'b'.repeat(64)}`;
const runtimeA = `h_${'a'.repeat(26)}`;
const runtimeB = `h_${'b'.repeat(26)}`;
const manifest = {
  schema: 2,
  proteinSets: {
    'sha256:set-a': { schema: 1, proteins: [{ featureAnalysisId: featureA, aaSha256: 'a'.repeat(64) }] },
    'sha256:set-b': { schema: 1, proteins: [{ featureAnalysisId: featureB, aaSha256: 'b'.repeat(64) }] }
  },
  recordAnalyses: {
    'sha256:analysis-a': { schema: 1, recordSourceId: 'A', proteinSetHash: 'sha256:set-a' },
    'sha256:analysis-b': { schema: 1, recordSourceId: 'B', proteinSetHash: 'sha256:set-b' }
  },
  recordInstances: {
    'record-1': {
      schema: 2,
      recordAnalysisId: 'sha256:analysis-a',
      runtimeBindingHash: 'sha256:runtime-binding-a',
      displayBindingHash: 'sha256:display-binding-a',
      runtimeIds: { [featureA]: runtimeA },
      featureMetadata: { [featureA]: { displayAlias: 'protein-a', exportOrdinal: null } }
    },
    'record-2': {
      schema: 2,
      recordAnalysisId: 'sha256:analysis-b',
      runtimeBindingHash: 'sha256:runtime-binding-b',
      displayBindingHash: 'sha256:display-binding-b',
      runtimeIds: { [featureB]: runtimeB },
      featureMetadata: { [featureB]: { displayAlias: 'protein-b', exportOrdinal: null } }
    }
  }
};

const proteinEntry = {
  schema: 4,
  kind: 'raw-losat',
  identityKind: 'protein',
  idEncoding: 'runtime-handle-v1',
  key: 'protein-key',
  text: `${runtimeA}\t${runtimeB}\t100\t1\t0\t0\t1\t1\t1\t1\t0\t50\n`,
  program: 'blastp',
  outfmt: '6',
  args: [],
  queryProteinSetHash: 'sha256:set-a',
  subjectProteinSetHash: 'sha256:set-b',
  queryRuntimeBindingHash: 'sha256:runtime-binding-a',
  subjectRuntimeBindingHash: 'sha256:runtime-binding-b',
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
const queryRuntimeIds = new Set([runtimeA]);
const subjectRuntimeIds = new Set([runtimeB]);
assert.equal(
  cache.rawProteinTextMatchesBindings(proteinEntry.text, queryRuntimeIds, subjectRuntimeIds),
  true
);
assert.equal(
  cache.rawProteinTextMatchesBindings(
    `${runtimeA}\t${runtimeB}\n`,
    queryRuntimeIds,
    subjectRuntimeIds
  ),
  false,
  'binding IDs alone are not a valid LOSAT outfmt 6 row'
);
assert.equal(
  cache.rawProteinTextMatchesBindings(
    `${proteinEntry.text.trimEnd()}\textra\n`,
    queryRuntimeIds,
    subjectRuntimeIds
  ),
  false,
  'LOSAT outfmt 6 rows must not contain unexpected columns'
);
assert.equal(
  cache.rawProteinTextMatchesBindings(
    `${runtimeA}\t${runtimeB}\tnot-a-number\t1\t0\t0\t1\t1\t1\t1\t0\t50\n`,
    queryRuntimeIds,
    subjectRuntimeIds
  ),
  false,
  'LOSAT outfmt 6 numeric columns must be parseable'
);
assert.equal(
  cache.rawProteinTextMatchesBindings(
    `${runtimeA}\t${runtimeB}\t100\t1\t0\t0\t1\t1\t1\t1\t0\t\n`,
    queryRuntimeIds,
    subjectRuntimeIds
  ),
  false,
  'LOSAT outfmt 6 numeric columns must not be empty'
);

const manifestWithQueryRuntimeIds = (runtimeIds) => ({
  ...manifest,
  recordInstances: {
    ...manifest.recordInstances,
    'record-1': {
      ...manifest.recordInstances['record-1'],
      runtimeIds
    }
  }
});
const mismatchedBindingManifest = manifestWithQueryRuntimeIds({
  [`f_${'c'.repeat(64)}`]: runtimeA
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
  cache.validateProteinIdentityManifest(manifestWithQueryRuntimeIds({
    [featureA]: 'not-a-runtime-handle'
  })),
  false,
  'a runtime handle must use the compact handle grammar'
);

const rawMap = new Map([['protein-key', { ...proteinEntry, key: undefined }]]);
assert.equal(
  cache.getCurrentRawLosatCacheEntry(
    rawMap,
    'protein-key',
    {
      program: 'blastp',
      outfmt: '6',
      args: [],
      queryRuntimeBindingHash: 'sha256:runtime-binding-a',
      subjectRuntimeBindingHash: 'sha256:runtime-binding-b'
    },
    manifest
  )?.entry?.text,
  proteinEntry.text
);
assert.equal(
  cache.getCurrentRawLosatCacheEntry(
    rawMap,
    'protein-key',
    {
      program: 'blastp',
      outfmt: '6',
      args: [],
      queryRuntimeBindingHash: 'sha256:runtime-binding-b',
      subjectRuntimeBindingHash: 'sha256:runtime-binding-a'
    },
    manifest
  ),
  null,
  'reverse query/subject bindings must not be a direct hit'
);

const invalidTextMap = new Map([['protein-key', {
  ...proteinEntry,
  key: undefined,
  text: `unknown\t${runtimeB}\t100\t1\t0\t0\t1\t1\t1\t1\t0\t50\n`
}]]);
assert.equal(
  cache.getCurrentRawLosatCacheEntry(invalidTextMap, 'protein-key', {
    program: 'blastp', outfmt: '6', args: [],
    queryRuntimeBindingHash: 'sha256:runtime-binding-a',
    subjectRuntimeBindingHash: 'sha256:runtime-binding-b'
  }, manifest),
  null
);

const v35Manifest = {
  schema: 1,
  proteinSets: manifest.proteinSets,
  recordAnalyses: manifest.recordAnalyses,
  recordInstances: {
    'record-1': {
      schema: 1,
      recordAnalysisId: 'sha256:analysis-a',
      bindingHash: 'sha256:v35-binding-a',
      transportIds: { [featureA]: `A@record-1|protein-a~${featureA}` }
    },
    'record-2': {
      schema: 1,
      recordAnalysisId: 'sha256:analysis-b',
      bindingHash: 'sha256:v35-binding-b',
      transportIds: { [featureB]: `B@record-2|protein-b~${featureB}` }
    }
  }
};
const v35ProteinEntry = {
  schema: 3,
  kind: 'raw-losat',
  identityKind: 'protein',
  key: 'v35-protein-key',
  text: (
    `A@record-1|protein-a~${featureA}\t` +
    `B@record-2|protein-b~${featureB}\t100\t1\t0\t0\t1\t1\t1\t1\t0\t50\n`
  ),
  program: 'blastp',
  outfmt: '6',
  args: [],
  queryProteinSetHash: 'sha256:set-a',
  subjectProteinSetHash: 'sha256:set-b',
  queryBindingHash: 'sha256:v35-binding-a',
  subjectBindingHash: 'sha256:v35-binding-b',
  queryRecordInstanceKey: 'record-1',
  subjectRecordInstanceKey: 'record-2'
};
assert.equal(cache.classifyRawLosatCacheEntry(v35ProteinEntry), 'protein-v35');
assert.equal(cache.validateV35ProteinIdentityManifest(v35Manifest), true);
assert.equal(cache.validateV35ProteinRawEntryReferences(v35ProteinEntry, v35Manifest), true);

const outfmt6ColumnIndexes = new Map([
  'query',
  'subject',
  'identity',
  'alignment_length',
  'mismatches',
  'gap_opens',
  'qstart',
  'qend',
  'sstart',
  'send',
  'evalue',
  'bitscore'
].map((column, index) => [column, index]));
const entryWithNumericCase = (entry, numericCase) => {
  const columns = entry.text.trimEnd().split('\t');
  columns[outfmt6ColumnIndexes.get(numericCase.column)] = numericCase.value;
  return { ...entry, text: `${columns.join('\t')}\n` };
};
for (const numericCase of numericContract) {
  const currentCase = entryWithNumericCase(proteinEntry, numericCase);
  const v35Case = entryWithNumericCase(v35ProteinEntry, numericCase);
  assert.equal(
    cache.rawProteinTextMatchesBindings(
      currentCase.text,
      queryRuntimeIds,
      subjectRuntimeIds
    ),
    numericCase.valid,
    `${numericCase.name}: raw text contract`
  );
  assert.equal(
    cache.validateProteinRawEntryReferences(currentCase, manifest),
    numericCase.valid,
    `${numericCase.name}: schema-4 validation`
  );
  assert.equal(
    cache.validateV35ProteinRawEntryReferences(v35Case, v35Manifest),
    numericCase.valid,
    `${numericCase.name}: version-35 validation`
  );
}

assert.deepEqual(
  cache.buildV35ProteinReferenceMap(v35Manifest, manifest),
  {
    [`A@record-1|protein-a~${featureA}`]: runtimeA,
    [`B@record-2|protein-b~${featureB}`]: runtimeB
  },
  'version-35 UI references resolve from manifest authority without a raw entry'
);
assert.throws(
  () => cache.buildV35ProteinReferenceMap(
    v35Manifest,
    {
      ...manifest,
      recordInstances: {
        ...manifest.recordInstances,
        'record-2': {
          ...manifest.recordInstances['record-2'],
          runtimeIds: {}
        }
      }
    }
  ),
  /valid version-35 and current manifests|membership changed/,
  'reference migration must fail closed when current membership changes'
);
const v35ManifestOnly = cache.createV35ProteinCandidateEnvelope([], v35Manifest);
assert.deepEqual(v35ManifestOnly, {
  schema: 1,
  sourceManifest: v35Manifest,
  entries: []
});
assert.deepEqual(
  cache.serializableV35ProteinCandidateEnvelope(v35ManifestOnly),
  v35ManifestOnly,
  'source manifest authority must survive save-before-generate with zero raw candidates'
);
const v35Pending = cache.createV35ProteinCandidateEnvelope(
  [v35ProteinEntry],
  v35Manifest
);
assert.equal(v35Pending.entries.length, 1);
assert.deepEqual(
  cache.serializableV35ProteinCandidateEnvelope(v35Pending),
  v35Pending
);
const v35Promoted = cache.transitionV35ProteinCandidate(
  v35Pending,
  0,
  'promoted'
);
assert.equal(
  cache.serializableV35ProteinCandidateEnvelope(v35Promoted).entries.length,
  0
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
  schema: 2, kind: 'derived-losatp-payload', key: 'v35-derived', payload: {}
}), true);
assert.equal(cache.isLosatDerivedCacheEntry({
  schema: 3,
  kind: 'derived-losatp-payload',
  idEncoding: 'runtime-handle-v1',
  key: 'current-derived',
  payload: {}
}, { allowLegacy: false }), true);

const zeroHitDerivedEntry = (mode, { includeIdentity = true } = {}) => {
  const payload = {
    pairs: [{
      pair_index: 0,
      query_index: 0,
      subject_index: 1,
      tsv: '',
      rows: [],
      hit_count: 0
    }],
    orthogroups: []
  };
  if (includeIdentity) {
    payload.identity = {
      cacheSchema: 3,
      idEncoding: 'runtime-handle-v1',
      converter: 'convert_losatp_blastp_pairs_to_genomic_payload',
      mode,
      rawCacheKeys: ['raw-key']
    };
  }
  if (mode === 'collinear') {
    Object.assign(payload, {
      collinearGroups: [],
      collinearGroupScope: 'adjacent_local',
      collinearityBlocks: []
    });
  }
  return {
    schema: 3,
    kind: 'derived-losatp-payload',
    idEncoding: 'runtime-handle-v1',
    key: `zero-hit-${mode}`,
    mode,
    payload
  };
};

for (const entry of [
  zeroHitDerivedEntry('orthogroup'),
  zeroHitDerivedEntry('collinear'),
  zeroHitDerivedEntry('collinear', { includeIdentity: false })
]) {
  assert.equal(
    cache.validateDerivedProteinReferences(entry, manifest),
    true,
    `${entry.mode} zero-hit results must remain cacheable`
  );
}

const zeroHitNearMisses = [
  ['arbitrary field', (payload) => { payload.note = 'arbitrary'; }],
  ['nonempty rows', (payload) => { payload.pairs[0].rows = [{}]; }],
  ['nonempty TSV', (payload) => { payload.pairs[0].tsv = 'unexpected'; }],
  ['nonzero hit count', (payload) => { payload.pairs[0].hit_count = 1; }],
  ['noninteger pair index', (payload) => { payload.pairs[0].pair_index = '0'; }],
  ['partial record indices', (payload) => { delete payload.pairs[0].subject_index; }],
  ['missing pair results', (payload) => { payload.pairs = []; }],
  ['nonempty orthogroups', (payload) => { payload.orthogroups = [{}]; }],
  ['nonempty collinear groups', (payload) => { payload.collinearGroups = [{}]; }],
  ['nonempty collinearity blocks', (payload) => { payload.collinearityBlocks = [{}]; }],
  ['invalid collinear scope', (payload) => { payload.collinearGroupScope = 'invalid'; }],
  ['null identity', (payload) => { payload.identity = null; }],
  ['invalid raw cache binding', (payload) => { payload.identity.rawCacheKeys = [null]; }],
  ['identity mode mismatch', (payload) => { payload.identity.mode = 'orthogroup'; }]
];
for (const [label, mutate] of zeroHitNearMisses) {
  const entry = zeroHitDerivedEntry('collinear');
  mutate(entry.payload);
  assert.equal(
    cache.validateDerivedProteinReferences(entry, manifest),
    false,
    `${label} must not qualify as an empty derived result`
  );
}

const compoundDerivedEntry = {
  schema: 3,
  kind: 'derived-losatp-payload',
  idEncoding: 'runtime-handle-v1',
  key: 'compound-derived',
  mode: 'collinear',
  payload: {
    rows: [{
      query_protein_id: `${runtimeA};${runtimeB}`,
      subject_protein_id: runtimeB
    }]
  }
};
assert.equal(
  cache.validateDerivedProteinReferences(compoundDerivedEntry, manifest),
  true,
  'collinear block rows may contain semicolon-delimited runtime handles'
);
assert.equal(
  cache.validateDerivedProteinReferences({
    ...compoundDerivedEntry,
    payload: {
      rows: [{
        query_protein_id: `${runtimeA};h_${'c'.repeat(26)}`,
        subject_protein_id: runtimeB
      }]
    }
  }, manifest),
  false,
  'every handle in a compound collinear reference must resolve'
);
const runtimeC = `h_${'c'.repeat(26)}`;
for (const unitReferenceKey of [
  'queryUnitId',
  'subjectUnitId',
  'query_unit_id',
  'subject_unit_id'
]) {
  assert.equal(
    cache.validateDerivedProteinReferences({
      ...compoundDerivedEntry,
      payload: {
        proteinId: runtimeA,
        [unitReferenceKey]: `${runtimeA};${runtimeB}`
      }
    }, manifest),
    true,
    `${unitReferenceKey} must accept manifest-owned runtime handles`
  );
  assert.equal(
    cache.validateDerivedProteinReferences({
      ...compoundDerivedEntry,
      payload: {
        proteinId: runtimeA,
        [unitReferenceKey]: `${runtimeA};${runtimeC}`
      }
    }, manifest),
    false,
    `${unitReferenceKey} must reject an unresolved runtime handle`
  );
}
assert.equal(
  cache.validateDerivedProteinReferences({
    ...compoundDerivedEntry,
    payload: {
      proteinId: runtimeA,
      query_unit_id: 'gbd_r0001_unit000001',
      subject_unit_id: 'gbd_r0002_unit000002'
    }
  }, manifest),
  true,
  'synthetic non-protein collinearity unit IDs must remain valid'
);
const supportingEdge = `${runtimeA}->${runtimeB}:rbh`;
const pathEdge = `og_1:0:${runtimeA}->1:${runtimeB}:rbh`;
const compoundEdgesByKey = {
  supportingEdge,
  supportingEdges: [supportingEdge],
  supporting_edge: supportingEdge,
  supporting_edges: [supportingEdge],
  edgeId: pathEdge,
  edgeIds: [pathEdge],
  edge_id: pathEdge,
  edge_ids: [pathEdge]
};
for (const [key, value] of Object.entries(compoundEdgesByKey)) {
  assert.equal(
    cache.validateDerivedProteinReferences({
      ...compoundDerivedEntry,
      payload: { [key]: value }
    }, manifest),
    true,
    `${key} must resolve both runtime-handle endpoints`
  );
  const unknownValue = Array.isArray(value)
    ? value.map((item) => item.replace(runtimeB, runtimeC))
    : value.replace(runtimeB, runtimeC);
  assert.equal(
    cache.validateDerivedProteinReferences({
      ...compoundDerivedEntry,
      payload: { [key]: unknownValue }
    }, manifest),
    false,
    `${key} must reject an unresolved runtime-handle endpoint`
  );
}
assert.equal(
  cache.validateDerivedProteinReferences({
    ...compoundDerivedEntry,
    payload: { supportingEdges: [`${runtimeA}->${runtimeB}`] }
  }, manifest),
  false,
  'compound edge references must use a recognized complete grammar'
);
for (const invalidCompoundValue of [null, 1, { edge: supportingEdge }]) {
  assert.equal(
    cache.validateDerivedProteinReferences({
      ...compoundDerivedEntry,
      payload: {
        proteinId: runtimeA,
        supportingEdges: invalidCompoundValue
      }
    }, manifest),
    false,
    'compound edge reference fields must be strings or arrays of strings'
  );
}

const derivedPayloadWithNote = (note) => ({
  ...compoundDerivedEntry,
  payload: { proteinId: runtimeA, note }
});
assert.equal(
  cache.validateDerivedProteinReferences(
    derivedPayloadWithNote('A plain p_r_ fragment is descriptive text, not a legacy ID.'),
    manifest
  ),
  true
);
for (const forbiddenReference of [
  'p_r_old_0_3_1_deadbeefdead',
  `A@record-1|protein-a~${featureA}`,
  featureA
]) {
  assert.equal(
    cache.validateDerivedProteinReferences(
      derivedPayloadWithNote(`embedded legacy reference: ${forbiddenReference}`),
      manifest
    ),
    false,
    `embedded legacy reference ${forbiddenReference} must be rejected`
  );
  assert.equal(
    cache.validateDerivedProteinReferences({
      ...compoundDerivedEntry,
      payload: {
        proteinId: runtimeA,
        [`embedded legacy reference: ${forbiddenReference}`]: true
      }
    }, manifest),
    false,
    `embedded legacy reference ${forbiddenReference} in a key must be rejected`
  );
}

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
        'record-1': {
          ...manifest.recordInstances['record-1'],
          runtimeBindingHash: 'different'
        }
      }
    }
  ]),
  /conflicting record instance/
);
