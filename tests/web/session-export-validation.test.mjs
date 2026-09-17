import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';
import { gunzipSync } from 'node:zlib';

if (!globalThis.crypto?.subtle) {
  Object.defineProperty(globalThis, 'crypto', {
    configurable: true,
    enumerable: true,
    value: webcrypto,
    writable: true
  });
}

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {
  body: { appendChild: () => {} },
  createElement: () => ({
    addEventListener: () => {},
    click: () => {},
    parentNode: null
  })
};

const {
  adoptCanonicalRenderArtifacts,
  exportSession,
  validateSessionLosatArtifacts
} = await import(
  '../../gbdraw/web/js/services/config.js'
);
const { buildCanonicalRenderRequest } = await import(
  '../../gbdraw/web/js/services/session-request.js'
);
const { state } = await import('../../gbdraw/web/js/state.js');

let inputReads = 0;
state.mode.value = 'circular';
state.cInputType.value = 'gb';
state.files.c_gb = {
  name: 'legacy.gbk',
  type: 'application/genbank',
  size: 12,
  lastModified: 1,
  async arrayBuffer() {
    inputReads += 1;
    return new TextEncoder().encode('LOCUS legacy').buffer;
  }
};
state.results.value = [{
  name: 'legacy.svg',
  content: '<svg xmlns="http://www.w3.org/2000/svg"></svg>'
}];
state.featureCatalog.value = null;

const saveError = /Generate again before using Save Session\./;
await assert.rejects(exportSession('legacy-session'), saveError);
assert.equal(inputReads, 0);

state.featureCatalog.value = { schema: 3, items: [] };
await assert.rejects(exportSession('invalid-catalog'), saveError);
assert.equal(inputReads, 0);

state.results.value = [];
state.featureCatalog.value = null;
state.adv.circular_track_slots_enabled = true;
state.adv.circular_track_slots_axis_index = 0;
state.adv.circular_track_slots.splice(
  0,
  state.adv.circular_track_slots.length,
  {
    id: 'features',
    renderer: 'features',
    enabled: true,
    side: 'inside',
    width: null,
    radius: null,
    inner_gap_px: null,
    outer_gap_px: null,
    z: 0,
    params: { lane_direction: 'inside' }
  }
);
const committed = buildCanonicalRenderRequest({
  state,
  filesData: {
    c_gb: {
      name: 'legacy.gbk',
      type: 'application/genbank',
      size: 12,
      lastModified: 1,
      encoding: 'base64',
      data: btoa('LOCUS legacy')
    }
  }
});
adoptCanonicalRenderArtifacts(committed, { adoptOwnedRequest: true });
state.adv.circular_track_slots[0].width = '16px';

const retainedCircularFiles = {
  blasts: [{ name: 'retained-blast.tsv' }],
  fastas: [{ name: 'retained-subject.fa' }],
  sources: [{ name: 'retained-source.fa' }]
};
state.files.c_conservation_blasts = retainedCircularFiles.blasts;
state.files.c_conservation_fastas = retainedCircularFiles.fastas;
state.files.c_conservation_sequence_sources = retainedCircularFiles.sources;
state.files.c_conservation_blasts_source = 'upload';
Object.assign(state.circularConservation, {
  enabled: false,
  source: 'upload',
  reference: 'query',
  labels: 'Retained',
  ring_width: 0.12,
  ring_gap: 0.03
});
state.circularConservation.series.splice(
  0,
  state.circularConservation.series.length,
  { label: 'Retained', sourceIndex: 0 }
);
const failedAdoption = structuredClone(committed);
failedAdoption.renderRequest.diagramOptions.tracks.circularTrackSlots[0].width = {
  value: 16,
  unit: 'px'
};
failedAdoption.resources['invalid-conservation'] = {
  kind: 'blast',
  name: 'invalid.tsv',
  type: 'text/tab-separated-values',
  size: 3,
  lastModified: 1,
  encoding: 'base64',
  data: '%%%'
};
failedAdoption.renderRequest.diagramOptions.conservationBlastFiles = [{
  resourceId: 'invalid-conservation'
}];
failedAdoption.webFiles = {
  conservationBlastSource: 'losat-cache'
};
assert.throws(
  () => adoptCanonicalRenderArtifacts(failedAdoption),
  /Invalid character/
);
assert.equal(state.files.c_conservation_blasts, retainedCircularFiles.blasts);
assert.equal(state.files.c_conservation_fastas, retainedCircularFiles.fastas);
assert.equal(
  state.files.c_conservation_sequence_sources,
  retainedCircularFiles.sources
);
assert.equal(state.files.c_conservation_blasts_source, 'upload');
assert.deepEqual(state.circularConservation.series, [{
  label: 'Retained',
  sourceIndex: 0
}]);
state.files.c_conservation_blasts = [];
state.files.c_conservation_fastas = [];
state.files.c_conservation_sequence_sources = [];
state.files.c_conservation_blasts_source = null;

let downloadedBlobs = 0;
let downloadedBlob = null;
let compressionAttempts = 0;
const originalCreateObjectUrl = URL.createObjectURL;
const OriginalCompressionStream = globalThis.CompressionStream;
URL.createObjectURL = (blob) => {
  downloadedBlobs += 1;
  downloadedBlob = blob;
  return 'blob:active-draft-session';
};
globalThis.CompressionStream = function CountingCompressionStream(...args) {
  compressionAttempts += 1;
  return new OriginalCompressionStream(...args);
};
try {
  const saved = await exportSession('divergent-active-draft');
  assert.equal(saved.status, 'saved');
  const session = JSON.parse(gunzipSync(
    Buffer.from(await saved.blob.arrayBuffer())
  ).toString('utf8'));
  assert.equal(session.config.adv.circular_track_slots[0].width, '16px');
  assert.equal(
    session.renderRequest.diagramOptions.tracks.circularTrackSlots[0].width,
    null
  );

  // Replacing a source leaves older raw entries in the live cache, but their
  // bindings cannot be published with the replacement's identity manifest.
  const feature = `f_${'a'.repeat(64)}`;
  const handle = `h_${'a'.repeat(26)}`;
  state.proteinIdentityManifest.value = {
    schema: 2,
    proteinSets: { current: { schema: 1, proteins: [{ featureAnalysisId: feature }] } },
    recordAnalyses: { current: { schema: 1, proteinSetHash: 'current' } },
    recordInstances: { record: {
      schema: 2, recordAnalysisId: 'current', runtimeBindingHash: 'binding',
      displayBindingHash: 'display', runtimeIds: { [feature]: handle },
      featureMetadata: { [feature]: { displayAlias: 'protein', exportOrdinal: null } }
    } }
  };
  const currentRaw = {
    schema: 4, kind: 'raw-losat', identityKind: 'protein', idEncoding: 'runtime-handle-v1',
    program: 'blastp', outfmt: '6', args: [],
    queryProteinSetHash: 'current', subjectProteinSetHash: 'current',
    queryRuntimeBindingHash: 'binding', subjectRuntimeBindingHash: 'binding',
    queryRecordInstanceKey: 'record', subjectRecordInstanceKey: 'record',
    text: `${handle}\t${handle}\t100\t1\t0\t0\t1\t1\t1\t1\t0\t50\n`
  };
  const entries = [
    ['old-displayed', { ...currentRaw, queryRuntimeBindingHash: 'old-binding' }],
    ['current', currentRaw],
    ['old-dormant', { ...currentRaw, subjectProteinSetHash: 'former-source' }],
    ['removed-record', { ...currentRaw, queryRecordInstanceKey: 'removed' }],
    ['valid-dormant', { ...currentRaw, args: ['--max-target-seqs', '5'], text: '' }],
    ['nucleotide', { schema: 2, kind: 'raw-losat', program: 'blastn', text: '' }]
  ];
  state.losatCache.value = new Map(entries);
  state.losatCacheInfo.value = [{ key: 'old-displayed' }, { key: 'current' }];
  const regenerated = await exportSession('replacement-source');
  const replacementSession = JSON.parse(gunzipSync(
    Buffer.from(await regenerated.blob.arrayBuffer())
  ).toString('utf8'));
  assert.deepEqual(replacementSession.losatCache.entries.map(entry => entry.key),
    ['current', 'valid-dormant', 'nucleotide']);
  assert.equal(replacementSession.losatCache.entries[0].text, currentRaw.text);
  assert.doesNotThrow(() => validateSessionLosatArtifacts(replacementSession, replacementSession.version));
  assert.deepEqual([...state.losatCache.value], entries, 'Save must not mutate live cache/History state');
  state.proteinIdentityManifest.value = null;
  await assert.rejects(exportSession('invalid-protein-manifest'), /valid protein identity manifest/);
  // A missing protein raw entry must not turn an invalid manifest into a
  // successful save with silently discarded identity metadata.
  for (const rawEntries of [[], [entries.at(-1)]]) {
    state.losatCache.value = new Map(rawEntries);
    state.losatCacheInfo.value = [];
    for (const invalidManifest of [null, { schema: -1 }, {
      ...structuredClone(replacementSession.proteinIdentityManifest),
      recordAnalyses: {}
    }]) {
      state.proteinIdentityManifest.value = invalidManifest;
      await assert.rejects(exportSession('invalid-manifest-without-protein-raw'),
        /valid protein identity manifest/);
      assert.equal(state.proteinIdentityManifest.value, invalidManifest);
      assert.deepEqual([...state.losatCache.value], rawEntries);
      assert.equal(downloadedBlobs, 2);
      assert.equal(compressionAttempts, 2);
    }
  }
  for (const rawEntries of [[], [entries.at(-1)]]) {
    state.losatCache.value = new Map(rawEntries);
    state.proteinIdentityManifest.value = structuredClone(session.proteinIdentityManifest);
    const valid = await exportSession(`valid-empty-manifest-${rawEntries.length}`);
    const document = JSON.parse(gunzipSync(Buffer.from(await valid.blob.arrayBuffer())).toString('utf8'));
    assert.equal(valid.status, 'saved');
    assert.deepEqual(document.proteinIdentityManifest, session.proteinIdentityManifest);
    assert.equal(document.losatCache.entries.length, rawEntries.length);
    assert.doesNotThrow(() => validateSessionLosatArtifacts(document, document.version));
  }
} finally {
  URL.createObjectURL = originalCreateObjectUrl;
  globalThis.CompressionStream = OriginalCompressionStream;
}
assert.ok(inputReads > 0, 'saving must bind the active input file');
assert.equal(compressionAttempts, 4);
assert.equal(downloadedBlobs, 4);
assert.ok(downloadedBlob instanceof Blob);
