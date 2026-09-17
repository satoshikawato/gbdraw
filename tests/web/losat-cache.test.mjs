import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
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
assert.equal(
  cache.classifyRawLosatCacheEntry({ ...proteinEntry, schema: 3 }),
  'invalid',
  'branch-internal protein raw schema 3 must not be accepted'
);
assert.equal(cache.validateProteinIdentityManifest(manifest), true);
const identityIndex = cache.buildValidatedProteinIdentityIndex(manifest);
assert.ok(identityIndex, 'a valid manifest must produce a reusable identity index');
const canonicalJson = (value) => {
  if (Array.isArray(value)) return `[${value.map(canonicalJson).join(',')}]`;
  if (value && typeof value === 'object') {
    return `{${Object.keys(value).sort().map((key) => (
      `${JSON.stringify(key)}:${canonicalJson(value[key])}`
    )).join(',')}}`;
  }
  return JSON.stringify(value);
};
for (const forbiddenKey of [
  'viewFeatureSvgId',
  'queryViewFeatureSvgId',
  'processed_view_feature_svg_id',
  'subjectViewFeatureHashParts',
  'renderedFeatureSvgId',
  'queryRenderedFeatureSvgId',
  'processed_rendered_svg_id'
]) {
  const maliciousManifest = JSON.parse(JSON.stringify(manifest));
  const instance = maliciousManifest.recordInstances['record-1'];
  instance.featureMetadata[featureA].extension = {
    [forbiddenKey]: 'presentation-only'
  };
  const displayPayload = {
    recordAnalysisId: instance.recordAnalysisId,
    recordSourceId: maliciousManifest.recordAnalyses[instance.recordAnalysisId].recordSourceId,
    recordInstanceKey: 'record-1',
    featureMetadata: instance.featureMetadata
  };
  instance.displayBindingHash = `sha256:${createHash('sha256')
    .update(canonicalJson(displayPayload), 'utf8')
    .digest('hex')}`;
  assert.equal(
    cache.validateProteinIdentityManifest(maliciousManifest),
    false,
    `${forbiddenKey} must not cross the presentation-independent manifest boundary`
  );
}
assert.equal(
  cache.validateProteinIdentityManifest({ ...manifest, schema: 1 }),
  false,
  'branch-internal top-level manifest schema 1 must not be accepted'
);
assert.equal(cache.validateProteinRawEntryReferences(proteinEntry, manifest), true);
assert.equal(
  cache.validateProteinRawEntryReferences(proteinEntry, manifest, { identityIndex }),
  true,
  'raw cache validation must accept the reusable identity index'
);
assert.equal(cache.releaseValidatedProteinIdentityIndex(identityIndex), true);
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
    Boolean(cache.getCurrentRawLosatCacheEntry(
      new Map([['protein-key', currentCase]]), 'protein-key', proteinEntry, manifest
    )),
    numericCase.valid,
    `${numericCase.name}: raw getter`
  );
}

const getProteinEntry = (entry = proteinEntry, identity = manifest, metadata = proteinEntry) => (
  cache.getCurrentRawLosatCacheEntry(new Map([['protein-key', entry]]), 'protein-key', metadata, identity)
);
assert.ok(getProteinEntry({ ...proteinEntry, text: '' }), 'empty raw is a successful hit');
assert.equal(getProteinEntry(proteinEntry, mismatchedBindingManifest), null);
for (const patch of [
  { schema: 3 }, { program: 'blastn' }, { outfmt: '7' }, { args: ['--max-target-seqs', '1'] },
  { searchContext: 'a'.repeat(64) }, { queryProteinSetHash: 'wrong' },
  { subjectProteinSetHash: 'wrong' }, { queryRecordInstanceKey: 'missing' },
  { subjectRecordInstanceKey: 'missing' }, { queryRuntimeBindingHash: 'wrong' },
  { subjectRuntimeBindingHash: 'wrong' },
  { text: `${proteinEntry.text.trimEnd()}\textra\n` },
  { text: `${runtimeA}\t${runtimeB}\n` },
  { text: proteinEntry.text.replace(runtimeB, 'unknown') }
]) assert.equal(getProteinEntry({ ...proteinEntry, ...patch }), null, JSON.stringify(patch));
assert.ok(getProteinEntry(
  { ...proteinEntry, searchContext: 'a'.repeat(64) }, manifest,
  { ...proteinEntry, searchContext: 'a'.repeat(64) }
));
const reverseEntry = {
  ...proteinEntry,
  text: proteinEntry.text.replace(`${runtimeA}\t${runtimeB}`, `${runtimeB}\t${runtimeA}`),
  queryProteinSetHash: proteinEntry.subjectProteinSetHash,
  subjectProteinSetHash: proteinEntry.queryProteinSetHash,
  queryRuntimeBindingHash: proteinEntry.subjectRuntimeBindingHash,
  subjectRuntimeBindingHash: proteinEntry.queryRuntimeBindingHash,
  queryRecordInstanceKey: proteinEntry.subjectRecordInstanceKey,
  subjectRecordInstanceKey: proteinEntry.queryRecordInstanceKey
};
assert.ok(getProteinEntry(reverseEntry, manifest, reverseEntry), 'matching reverse direction hits');

// Count actual validation entry points, without a production telemetry API.
const probePath = join(tempRoot, 'counted-cache.js');
let probeSource = await readFile(join(tempRoot, 'losat-cache.js'), 'utf8');
for (const [signature, counter] of [
  ['export const validateProteinIdentityManifest = (manifest) => {', 'manifest'],
  ['export const rawProteinTextMatchesBindings = (text, queryIds, subjectIds) => {', 'tsv']
]) {
  assert.ok(probeSource.includes(signature));
  probeSource = probeSource.replace(signature, `${signature}\n  counts.${counter}++;${counter === 'manifest' ? ' manifestCalls.push(manifest);' : ''}`);
}
probeSource = probeSource.replace(
  /new Set\(Object.values\((\w+)\.runtimeIds\)\)/g,
  'new Set((counts.runtimeSets++, Object.values($1.runtimeIds)))'
);
probeSource = probeSource.replaceAll(
  "metadata.displayAlias.normalize('NFC').trim()",
  "(aliasNormalizations++, metadata.displayAlias.normalize('NFC').trim())"
);
await writeFile(probePath, 'export let aliasNormalizations = 0;\nexport const manifestCalls = [];\nexport const counts = {manifest: 0, tsv: 0, runtimeSets: 0};\n' + probeSource);
const counted = await import(pathToFileURL(probePath));
assert.ok(counted.getCurrentRawLosatCacheEntry(rawMap, 'protein-key', proteinEntry, manifest));
assert.deepEqual(counted.counts, {manifest: 1, tsv: 1, runtimeSets: 2}, 'one validation per getter hit');
const resetCounts = () => Object.assign(counted.counts, {manifest: 0, tsv: 0, runtimeSets: 0});
resetCounts();
const privateManifest = counted.mergeProteinIdentityManifests([manifest]);
// The merge owns deep copies, including nested runtime IDs and metadata.
assert.notStrictEqual(privateManifest.recordInstances['record-1'].runtimeIds,
  manifest.recordInstances['record-1'].runtimeIds);
resetCounts();
const batchIndex = counted.buildValidatedProteinIdentityIndex(privateManifest);
assert.deepEqual(counted.counts, {manifest: 1, tsv: 0, runtimeSets: 2});
resetCounts();
for (const entry of [proteinEntry, reverseEntry, {...proteinEntry, text: ''}]) {
  assert.ok(counted.getCurrentRawLosatCacheEntry(
    new Map([['protein-key', entry]]), 'protein-key', entry, privateManifest, {identityIndex: batchIndex}
  ));
}
assert.deepEqual(counted.counts, {manifest: 0, tsv: 3, runtimeSets: 0}, 'batch reuses only manifest work');
const changedManifest = structuredClone(privateManifest);
changedManifest.recordInstances['record-1'].runtimeIds[featureA] = `h_${'c'.repeat(26)}`;
assert.equal(counted.getCurrentRawLosatCacheEntry(
  rawMap, 'protein-key', proteinEntry, changedManifest, {identityIndex: batchIndex}
), null, 'an index for another manifest cannot approve stale IDs');
for (const numericCase of numericContract) {
  assert.equal(Boolean(counted.getCurrentRawLosatCacheEntry(
    new Map([['protein-key', entryWithNumericCase(proteinEntry, numericCase)]]),
    'protein-key', proteinEntry, privateManifest, {identityIndex: batchIndex}
  )), numericCase.valid, `${numericCase.name}: indexed getter`);
}
assert.equal(counted.releaseValidatedProteinIdentityIndex(batchIndex), true);
assert.equal(counted.releaseValidatedProteinIdentityIndex(batchIndex), false);
privateManifest.recordInstances['record-1'].runtimeIds[featureA] = `h_${'c'.repeat(26)}`;
assert.equal(counted.getCurrentRawLosatCacheEntry(
  rawMap, 'protein-key', proteinEntry, privateManifest, {identityIndex: batchIndex}
), null, 'a released index cannot approve changed contents of the same object');
privateManifest.schema = 1;
assert.equal(counted.getCurrentRawLosatCacheEntry(
  rawMap, 'protein-key', proteinEntry, privateManifest, {identityIndex: batchIndex}
), null, 'released index falls back to full manifest validation');


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
  schema: 2, kind: 'derived-losatp-payload', key: 'unsupported-derived', payload: {}
}), false);
assert.equal(cache.isLosatDerivedCacheEntry({
  schema: 3,
  kind: 'derived-losatp-payload',
  idEncoding: 'runtime-handle-v1',
  key: 'current-derived',
  payload: {
    queryViewFeatureSvgId: 'processed-view-id',
    viewTransform: { length: 100, reverse: true }
  }
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

// Execute the production pair-loop body with real validators. The surrounding
// asynchronous I/O is controlled so an exception/cancel occurs while the loan
// is active, rather than after the loop has already released its index.
const runSource = await readFile(join(repoRoot, 'gbdraw/web/js/app/run-analysis.js'), 'utf8');
// Exercise the actual Generate merge boundary, including its user-facing error.
const mergeStart = runSource.indexOf('          const manifests = proteinEntries.map(');
const mergeEnd = runSource.indexOf('          const legacyReferenceIds =', mergeStart);
assert.ok(mergeStart > 0 && mergeEnd > mergeStart);
const generateMerge = new Function(
  'proteinEntries', 'mergeProteinIdentityManifests', 'validateProteinIdentityManifest',
  `let workingProteinIdentityManifest;\n${runSource.slice(mergeStart, mergeEnd)}
   return workingProteinIdentityManifest;`
);
const mergeForGenerate = (inputs) => generateMerge(
  inputs.map((identityManifest) => ({ identityManifest })),
  counted.mergeProteinIdentityManifests, counted.validateProteinIdentityManifest
);
const firstRecord = structuredClone({
  ...manifest, recordInstances: { 'record-1': manifest.recordInstances['record-1'] }
});
const secondRecord = structuredClone({
  ...manifest, recordInstances: { 'record-2': manifest.recordInstances['record-2'] }
});
const mergeInputs = [firstRecord, secondRecord, cache.emptyProteinIdentityManifest()];
const originalInputs = structuredClone(mergeInputs);
counted.manifestCalls.length = 0;
const generatedManifest = mergeForGenerate(mergeInputs);
assert.deepEqual(generatedManifest, manifest);
assert.deepEqual(mergeInputs, originalInputs);
assert.deepEqual(counted.manifestCalls, [...mergeInputs, generatedManifest],
  'Generate validates R inputs once each, then the distinct merged manifest once');
for (let index = 0; index < mergeInputs.length; index++) {
  assert.equal(counted.manifestCalls[index], mergeInputs[index]);
}
assert.equal(counted.manifestCalls.at(-1), generatedManifest);
for (const map of ['proteinSets', 'recordAnalyses', 'recordInstances']) {
  for (const [key, value] of Object.entries(firstRecord[map])) {
    assert.notEqual(generatedManifest[map][key], value, `merged ${map} values are copied`);
  }
}
generatedManifest.proteinSets['sha256:set-a'].proteins[0].aaSha256 = 'changed';
generatedManifest.recordInstances['record-1'].runtimeIds[featureA] = runtimeB;
generatedManifest.recordInstances['record-1'].featureMetadata[featureA].displayAlias = 'changed';
assert.deepEqual(mergeInputs, originalInputs, 'nested merged edits cannot alter source inputs');

const reloadMessage = 'Protein comparison metadata could not be validated. Reload the page and try again.';
const invalidInputMessage = 'Cannot merge an invalid protein identity manifest.';
const mergedMessage = 'Merged protein identity manifest is invalid.';
for (const invalid of [null, {}, { ...firstRecord, schema: 1 }]) {
  for (let position = 0; position < mergeInputs.length; position++) {
    const inputs = mergeInputs.with(position, invalid);
    assert.throws(() => mergeForGenerate(inputs), { name: 'Error', message: reloadMessage });
    assert.throws(() => cache.mergeProteinIdentityManifests(inputs), { name: 'Error', message: invalidInputMessage });
  }
}
for (const [map, key, field, value, label] of [
  ['proteinSets', 'sha256:set-a', 'extra', true, 'protein set'],
  ['recordAnalyses', 'sha256:analysis-a', 'recordSourceId', 'different', 'record analysis'],
  ['recordInstances', 'record-1', 'runtimeBindingHash', 'different', 'record instance']
]) {
  const conflict = structuredClone(firstRecord);
  conflict[map][key][field] = value;
  assert.equal(cache.validateProteinIdentityManifest(conflict), true);
  const conflictPattern = new RegExp(`conflicting ${label}`);
  assert.throws(() => mergeForGenerate([firstRecord, conflict]), conflictPattern);
  assert.throws(() => cache.mergeProteinIdentityManifests([firstRecord, conflict]), conflictPattern);
  for (let position = 0; position <= 2; position++) {
    const inputs = [firstRecord, conflict].toSpliced(position, 0, null);
    assert.throws(() => mergeForGenerate(inputs), { name: 'Error', message: reloadMessage },
      'Generate prioritizes invalid input even after an earlier merge conflict');
    assert.throws(() => cache.mergeProteinIdentityManifests(inputs),
      position === 2 ? conflictPattern : { name: 'Error', message: invalidInputMessage },
      'the default helper keeps its existing interleaved failure precedence');
  }
}
const collidingRuntime = structuredClone(secondRecord);
collidingRuntime.recordInstances['record-2'].runtimeIds[featureB] = runtimeA;
assert.equal(cache.validateProteinIdentityManifest(collidingRuntime), true);
for (const merge of [cache.mergeProteinIdentityManifests, mergeForGenerate]) {
  assert.throws(() => merge([firstRecord, collidingRuntime]), { name: 'Error', message: mergedMessage });
  assert.deepEqual(merge([]), cache.emptyProteinIdentityManifest());
  assert.deepEqual(merge([firstRecord, firstRecord]), firstRecord, 'identical identities deduplicate');
}
for (const input of [undefined, null, {}, 'not an array']) {
  assert.deepEqual(cache.mergeProteinIdentityManifests(input), cache.emptyProteinIdentityManifest());
}
console.log('RW-03: 3 input checks + 1 merged check; default/Generate failure precedence and deep copy passed');

const loopStart = runSource.indexOf('          const identityIndex = useProteinBlastp && preparedJobs.length > 0');
const loopEnd = runSource.indexOf('          const sourceJobs = [];', loopStart);
assert.ok(loopStart > 0 && loopEnd > loopStart);
const loopSource = runSource.slice(loopStart, loopEnd);
const exercisePairLoop = async ({failure = null, promote = false, mutateSource = false} = {}) => {
  const externalManifest = structuredClone(manifest);
  const workingManifest = counted.mergeProteinIdentityManifests([externalManifest]);
  const entries = [proteinEntry, reverseEntry, {...proteinEntry, text: ''}];
  const cacheMap = new Map(promote ? [] : entries.map((entry, index) => [String(index), entry]));
  let loan = null, released = false, reads = 0;
  const cancellation = new Error('canceled inside pair preparation');
  const context = {
    useProteinBlastp: true,
    preparedJobs: entries.map((entry, index) => ({
      spec: {ordinal: index, queryIndex: 0, subjectIndex: 1, edgeKey: String(index)},
      losatArgs: [], cacheMetadata: entry, batch: {}
    })),
    workingProteinIdentityManifest: workingManifest,
    buildValidatedProteinIdentityIndex(value) {
      loan = counted.buildValidatedProteinIdentityIndex(value);
      return loan;
    },
    releaseValidatedProteinIdentityIndex(value) {
      assert.strictEqual(value, loan);
      released = counted.releaseValidatedProteinIdentityIndex(value);
    },
    throwIfGenerationCanceled() { if (failure === 'cancel' && reads >= 3) throw cancellation; },
    async getSeqEntry() {
      reads++;
      await Promise.resolve();
      if (failure === 'error' && reads >= 3) throw new Error('sequence preparation failed');
      if (mutateSource) {
        externalManifest.recordInstances['record-1'].runtimeIds[featureA] = `h_${'c'.repeat(26)}`;
        externalManifest.schema = 1;
      }
      return {sequenceKey: 'sequence', fasta: 'data'};
    },
    proteinCacheKeys: ['0', '1', '2'],
    async getSeqHash() { return 'hash'; },
    sequenceEntriesByKey: new Map(), cacheMap,
    getReusableLosatCacheEntry(map, key, metadata, identity, index) {
      assert.strictEqual(index, loan);
      return counted.getCurrentRawLosatCacheEntry(map, key, metadata, identity, {identityIndex: index});
    },
    async tryPromoteLegacyProteinEntry({cacheKey, metadata, identityIndex: index}) {
      await Promise.resolve();
      assert.strictEqual(index, loan);
      cacheMap.set(cacheKey, entries[Number(cacheKey)]);
      return counted.getCurrentRawLosatCacheEntry(cacheMap, cacheKey, metadata, workingManifest, {identityIndex: index});
    },
    promoteRawLosatCacheEntry() {},
    losatTiming: {totalPairs: 0, cacheHits: 0, cacheMisses: 0},
    comparisonResolution: {edges: []}, losatPairs: [], cacheInfo: [],
    buildCacheFilename() { return 'pair.tsv'; },
    pendingJobKeys: new Set(), losatJobs: [], losatProgram: {value: 'blastp'}, losat: {outfmt: '6'}
  };
  const execute = new Function('context', `return (async () => {
    const {${Object.keys(context).join(',')}} = context;
    ${loopSource}
  })()`);
  resetCounts();
  if (failure) await assert.rejects(execute(context), failure === 'cancel' ? /canceled inside/ : /sequence preparation failed/);
  else await execute(context);
  assert.equal(released, true, 'every exit releases the loop index');
  assert.equal(counted.releaseValidatedProteinIdentityIndex(loan), false);
  assert.equal(counted.counts.manifest, 1, 'one manifest validation in the pair loop');
  assert.equal(counted.counts.runtimeSets, 2, 'one runtime ID Set per record');
  assert.equal(counted.counts.tsv, failure ? 1 : 3, 'one TSV validation for every hit, including promotion');
  if (!failure) assert.equal(context.losatTiming.cacheHits, 3);
  // The very same manifest object may change once the loan has ended. Its old
  // index must no longer bypass ID or manifest validation on a later lookup.
  workingManifest.schema = 1;
  assert.equal(counted.getCurrentRawLosatCacheEntry(
    new Map([['protein-key', proteinEntry]]), 'protein-key', proteinEntry, workingManifest, {identityIndex: loan}
  ), null);
};
await exercisePairLoop({mutateSource: true});
await exercisePairLoop({promote: true, mutateSource: true});
await exercisePairLoop({failure: 'error'});
await exercisePairLoop({failure: 'cancel'});
await exercisePairLoop(); // Retry must acquire a new index.

// Alias validity and ordinal grouping must use the same NFC-trimmed value.
for (const [alias, valid] of [
  ['protein-a', true], [' e\u0301 ', true], ['\u00e9', true], ['\u200b', true],
  ['', false], [' \t\r\n\u00a0\u3000', false],
  [null, false], [undefined, false], [42, false], [false, false], [[], false], [{}, false]
]) {
  const candidate = structuredClone(manifest);
  candidate.recordInstances['record-1'].featureMetadata[featureA].displayAlias = alias;
  const before = structuredClone(candidate);
  assert.equal(cache.validateProteinIdentityManifest(candidate), valid);
  assert.equal(Boolean(getProteinEntry(proteinEntry, candidate)), valid);
  assert.deepEqual(candidate, before, 'validation must not rewrite the stored alias');
}
for (const metadata of [null, [], {}, { displayAlias: 'valid', exportOrdinal: 1 }]) {
  const candidate = structuredClone(manifest);
  candidate.recordInstances['record-1'].featureMetadata[featureA] = metadata;
  assert.equal(cache.validateProteinIdentityManifest(candidate), false);
}
for (const [left, right, ordinals, valid] of [
  [' e\u0301 ', '\u00e9', [1, 2], true],
  [' e\u0301 ', '\u00e9', [2, 1], false],
  [' e\u0301 ', '\u00e9', [null, null], false],
  [' e\u0301 ', '\u00e9', ['1', '2'], false],
  ['A', 'a', [null, null], true],
  ['\uff21', 'A', [null, undefined], true]
]) {
  const candidate = structuredClone(firstRecord);
  candidate.proteinSets['sha256:set-a'].proteins.push({ featureAnalysisId: featureB });
  const instance = candidate.recordInstances['record-1'];
  instance.runtimeIds = { [featureB]: runtimeB, [featureA]: runtimeA };
  instance.featureMetadata = {
    [featureB]: { displayAlias: right, exportOrdinal: ordinals[1] },
    [featureA]: { displayAlias: left, exportOrdinal: ordinals[0] }
  };
  const before = structuredClone(candidate);
  assert.equal(cache.validateProteinIdentityManifest(candidate), valid,
    'NFC collisions use feature-ID order, preserving case and compatibility distinctions');
  assert.deepEqual(candidate, before);
}
const blankAlias = structuredClone(firstRecord);
blankAlias.recordInstances['record-1'].featureMetadata[featureA].displayAlias = ' \t\u3000';
const aliasConflict = structuredClone(firstRecord);
aliasConflict.recordInstances['record-1'].featureMetadata[featureA].displayAlias = 'different';
assert.throws(() => mergeForGenerate([firstRecord, aliasConflict, blankAlias]),
  { name: 'Error', message: reloadMessage });
assert.throws(() => cache.mergeProteinIdentityManifests([firstRecord, aliasConflict, blankAlias]),
  /conflicting record instance/);
assert.throws(() => cache.mergeProteinIdentityManifests([blankAlias]),
  { name: 'Error', message: invalidInputMessage });
console.log('RW-04: alias validity, NFC collisions, ordinals, immutability and failure precedence passed');
const aliasCountBefore = counted.aliasNormalizations;
assert.equal(counted.validateProteinIdentityManifest(manifest), true);
const aliasOperations = counted.aliasNormalizations - aliasCountBefore;
console.log(`RW-04: 2 valid feature aliases, ${aliasOperations} NFC/trim operations`);
assert.equal(aliasOperations, 2, 'one normalization per valid feature alias in each manifest validation');
