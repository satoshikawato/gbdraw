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
const alerts = [];
globalThis.alert = (message) => alerts.push(String(message));

const { importSession, validateSessionLosatArtifacts } = await import(
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

const proteinRawEntry = (key) => ({
  schema: 4,
  kind: 'raw-losat',
  idEncoding: 'runtime-handle-v1',
  key,
  program: 'blastp',
  identityKind: 'protein',
  queryProteinSetHash: 'sha256:query',
  subjectProteinSetHash: 'sha256:subject',
  queryRuntimeBindingHash: 'sha256:query-binding',
  subjectRuntimeBindingHash: 'sha256:subject-binding',
  queryRecordInstanceKey: 'query',
  subjectRecordInstanceKey: 'subject',
  text: ''
});

const derivedEntry = (key) => ({
  schema: 3,
  kind: 'derived-losatp-payload',
  idEncoding: 'runtime-handle-v1',
  key,
  payload: {}
});

const session = (rawEntries, derivedEntries) => ({
  losatCache: { entries: rawEntries },
  losatDerivedCache: { entries: derivedEntries },
  proteinIdentityManifest: {
    schema: 2,
    proteinSets: {},
    recordAnalyses: {},
    recordInstances: {}
  }
});

for (const version of [36, 37]) {
  test(`version-${version} import rejects duplicate raw LOSAT cache keys`, () => {
    assert.throws(
      () => validateSessionLosatArtifacts(
        session([rawEntry('raw-key'), rawEntry('raw-key')], []),
        version
      ),
      /Duplicate LOSAT cache key/
    );
  });

  test(`version-${version} import rejects duplicate derived LOSATP cache keys`, () => {
    assert.throws(
      () => validateSessionLosatArtifacts(
          session(
            [],
            [derivedEntry('derived-key'), derivedEntry('derived-key')]
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
          session([rawEntry(invalidKey)], []),
          version
        ),
        /LOSAT cache entry at losatCache\.entries\[0\] requires a key/
      );
    });
  }

  test(`version-${version} import rejects an empty protein raw cache key`, () => {
    assert.throws(
      () => validateSessionLosatArtifacts(
        session([proteinRawEntry('')], []),
        version
      ),
      /LOSAT cache entry at losatCache\.entries\[0\] requires a key/
    );
  });

  for (const field of ['losatCache', 'losatDerivedCache']) {
    test(`version-${version} import rejects a non-object ${field}`, () => {
      const malformed = session([], []);
      malformed[field] = [];
      assert.throws(
        () => validateSessionLosatArtifacts(malformed, version),
        new RegExp(`Session ${field} must be an object when present`)
      );
    });

    test(`version-${version} import rejects a non-array ${field}.entries`, () => {
      const malformed = session([], []);
      malformed[field] = { entries: null };
      assert.throws(
        () => validateSessionLosatArtifacts(malformed, version),
        new RegExp(`Session ${field}\\.entries must be an array`)
      );
    });
  }
  test(`version-${version} import rejects branch-internal raw schema 3`, () => {
    const malformed = session([{ ...proteinRawEntry('raw-key'), schema: 3 }], []);
    assert.throws(
      () => validateSessionLosatArtifacts(malformed, version),
      new RegExp(`Session version ${version} contains a non-current raw LOSAT entry`)
    );
  });

  test(`version-${version} import rejects branch-internal derived schema 2`, () => {
    const malformed = session([], [{ ...derivedEntry('derived-key'), schema: 2 }]);
    assert.throws(
      () => validateSessionLosatArtifacts(malformed, version),
      new RegExp(`Session version ${version} contains an invalid derived LOSATP entry`)
    );
  });
}

for (const unsupportedVersion of [34, 35]) {
  test(`session import rejects branch-internal version ${unsupportedVersion}`, async () => {
    alerts.length = 0;
    const file = new Blob([JSON.stringify({
      format: 'gbdraw-session',
      version: unsupportedVersion
    })], { type: 'application/json' });
    const event = { target: { files: [file], value: 'selected' } };
    const result = await importSession(event);
    assert.equal(result.status, 'error');
    assert.match(
      String(result.error?.message || ''),
      new RegExp(`Unsupported session version: ${unsupportedVersion}`)
    );
    assert.deepEqual(alerts, [
      `Failed to load session: Unsupported session version: ${unsupportedVersion}.`
    ]);
    assert.equal(event.target.value, '');
  });
}
