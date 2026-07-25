import assert from 'node:assert/strict';
import test from 'node:test';

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

const { validateSessionLosatArtifacts } = await import(
  '../../gbdraw/web/js/services/config.js'
);

const rawEntry = (key) => ({
  schema: 2,
  kind: 'raw-losat',
  key,
  program: 'blastn',
  identityKind: 'nucleotide',
  text: ''
});

const proteinRawEntry = (key, version) => ({
  schema: version === 35 ? 3 : 4,
  kind: 'raw-losat',
  ...(version === 36 ? { idEncoding: 'runtime-handle-v1' } : {}),
  key,
  program: 'blastp',
  identityKind: 'protein',
  queryProteinSetHash: 'sha256:query',
  subjectProteinSetHash: 'sha256:subject',
  ...(version === 35
    ? {
        queryBindingHash: 'sha256:query-binding',
        subjectBindingHash: 'sha256:subject-binding'
      }
    : {
        queryRuntimeBindingHash: 'sha256:query-binding',
        subjectRuntimeBindingHash: 'sha256:subject-binding'
      }),
  queryRecordInstanceKey: 'query',
  subjectRecordInstanceKey: 'subject',
  text: ''
});

const derivedEntry = (key, version) => ({
  schema: version === 35 ? 2 : 3,
  kind: 'derived-losatp-payload',
  ...(version === 36 ? { idEncoding: 'runtime-handle-v1' } : {}),
  key,
  payload: {}
});

const session = (version, rawEntries, derivedEntries) => ({
  losatCache: { entries: rawEntries },
  losatDerivedCache: { entries: derivedEntries },
  proteinIdentityManifest: version === 36
    ? {
        schema: 2,
        proteinSets: {},
        recordAnalyses: {},
        recordInstances: {}
      }
    : undefined
});

for (const version of [35, 36]) {
  test(`version-${version} import rejects duplicate raw LOSAT cache keys`, () => {
    assert.throws(
      () => validateSessionLosatArtifacts(
        session(version, [rawEntry('raw-key'), rawEntry('raw-key')], []),
        version
      ),
      /Duplicate LOSAT cache key/
    );
  });

  test(`version-${version} import rejects duplicate derived LOSATP cache keys`, () => {
    assert.throws(
      () => validateSessionLosatArtifacts(
        session(
          version,
          [],
          [derivedEntry('derived-key', version), derivedEntry('derived-key', version)]
        ),
        version
      ),
      /Duplicate derived LOSATP cache key/
    );
  });

  for (const invalidKey of [undefined, '', null, 1]) {
    test(`version-${version} import rejects a nucleotide raw cache key of ${String(invalidKey)}`, () => {
      assert.throws(
        () => validateSessionLosatArtifacts(
          session(version, [rawEntry(invalidKey)], []),
          version
        ),
        /LOSAT cache entry at losatCache\.entries\[0\] requires a key/
      );
    });
  }

  test(`version-${version} import rejects an empty protein raw cache key`, () => {
    assert.throws(
      () => validateSessionLosatArtifacts(
        session(version, [proteinRawEntry('', version)], []),
        version
      ),
      /LOSAT cache entry at losatCache\.entries\[0\] requires a key/
    );
  });

  for (const field of ['losatCache', 'losatDerivedCache']) {
    test(`version-${version} import rejects a non-object ${field}`, () => {
      const malformed = session(version, [], []);
      malformed[field] = [];
      assert.throws(
        () => validateSessionLosatArtifacts(malformed, version),
        new RegExp(`Session ${field} must be an object when present`)
      );
    });

    test(`version-${version} import rejects a non-array ${field}.entries`, () => {
      const malformed = session(version, [], []);
      malformed[field] = { entries: null };
      assert.throws(
        () => validateSessionLosatArtifacts(malformed, version),
        new RegExp(`Session ${field}\\.entries must be an array`)
      );
    });
  }
}

test('version-35 import rejects a present-invalid manifest without protein cache entries', () => {
  const malformed = session(35, [], []);
  malformed.proteinIdentityManifest = {};
  assert.throws(
    () => validateSessionLosatArtifacts(malformed, 35),
    /Session version 35 contains invalid protein LOSAT artifacts/
  );
});
