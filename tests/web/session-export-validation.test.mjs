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
  adoptCanonicalRenderArtifacts,
  exportSession
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
adoptCanonicalRenderArtifacts(committed);
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

let downloadedBlobs = 0;
let compressionAttempts = 0;
const originalCreateObjectUrl = URL.createObjectURL;
const OriginalCompressionStream = globalThis.CompressionStream;
URL.createObjectURL = () => {
  downloadedBlobs += 1;
  return 'blob:unexpected-session';
};
globalThis.CompressionStream = function CountingCompressionStream(...args) {
  compressionAttempts += 1;
  return new OriginalCompressionStream(...args);
};
try {
  await assert.rejects(
    exportSession('stale-active-draft'),
    /Generate again before using Save Session\. The active Custom Track draft has changed/
  );
} finally {
  URL.createObjectURL = originalCreateObjectUrl;
  globalThis.CompressionStream = OriginalCompressionStream;
}
assert.equal(inputReads, 0, 'active-draft rejection must happen before resource reads');
assert.equal(compressionAttempts, 0, 'active-draft rejection must happen before gzip Blob creation');
assert.equal(downloadedBlobs, 0, 'active-draft rejection must not create a download Blob URL');
