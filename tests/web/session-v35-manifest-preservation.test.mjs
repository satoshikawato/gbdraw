import assert from 'node:assert/strict';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {};

const {
  buildSessionLegacyArtifacts,
  quarantineV35ProteinArtifacts
} = await import(
  '../../gbdraw/web/js/services/config.js'
);

const featureId = `f_${'a'.repeat(64)}`;
const transportId = `record@record-1|protein-a~${featureId}`;
const sourceManifest = {
  schema: 1,
  proteinSets: {
    'sha256:set': {
      schema: 1,
      proteins: [{ featureAnalysisId: featureId, aaSha256: 'a'.repeat(64) }]
    }
  },
  recordAnalyses: {
    'sha256:analysis': {
      schema: 1,
      recordSourceId: 'record',
      proteinSetHash: 'sha256:set'
    }
  },
  recordInstances: {
    'record-1': {
      schema: 1,
      recordAnalysisId: 'sha256:analysis',
      bindingHash: 'sha256:binding',
      transportIds: { [featureId]: transportId }
    }
  }
};
const derivedEntry = {
  schema: 2,
  kind: 'derived-losatp-payload',
  key: 'derived-v35',
  payload: {
    orthogroups: [{
      members: [{ proteinId: transportId }]
    }]
  }
};
const source = {
  losatCache: { entries: [] },
  losatDerivedCache: { entries: [derivedEntry] },
  proteinIdentityManifest: sourceManifest
};

const migrated = quarantineV35ProteinArtifacts(source, 35);
assert.deepEqual(migrated.losatCache, { entries: [] });
assert.deepEqual(migrated.losatDerivedCache, { entries: [] });
assert.deepEqual(migrated.proteinIdentityManifest, {
  schema: 2,
  proteinSets: {},
  recordAnalyses: {},
  recordInstances: {}
});
assert.deepEqual(
  migrated.legacyArtifacts.proteinRawV35Candidates,
  {
    schema: 1,
    sourceManifest,
    entries: []
  },
  'version-35 source manifest must be quarantined even without raw candidates'
);
assert.deepEqual(
  migrated.legacyArtifacts.proteinDerivedV35Evidence,
  { schema: 1, entries: [derivedEntry] }
);
assert.deepEqual(source.proteinIdentityManifest, sourceManifest, 'migration must not mutate input');
assert.deepEqual(
  buildSessionLegacyArtifacts({
    legacyRawCandidates: { schema: 1, entries: [] },
    legacyV35RawCandidates: migrated.legacyArtifacts.proteinRawV35Candidates,
    legacyDerivedEvidence: { schema: 1, entries: [] },
    legacyV35DerivedEvidence: { schema: 1, entries: [derivedEntry] }
  }),
  migrated.legacyArtifacts,
  'save-before-generate must serialize v35 source authority with zero raw candidates'
);

console.log('version-35 manifest preservation tests passed');
