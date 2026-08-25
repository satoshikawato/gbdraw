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
  exportSession
} = await import(
  '../../gbdraw/web/js/services/config.js'
);
const { buildCanonicalRenderRequest } = await import(
  '../../gbdraw/web/js/services/session-request.js'
);
const { resolveLinearComparisonPlan } = await import(
  '../../gbdraw/web/js/app/linear-comparisons.js'
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
} finally {
  URL.createObjectURL = originalCreateObjectUrl;
  globalThis.CompressionStream = OriginalCompressionStream;
}
assert.ok(inputReads > 0, 'saving must bind the active input file');
assert.equal(compressionAttempts, 1);
assert.equal(downloadedBlobs, 1);
assert.ok(downloadedBlob instanceof Blob);

const linearFile = (name, recordId) => {
  const text = [
    `LOCUS       ${recordId.padEnd(16)} 12 bp    DNA     linear`,
    `DEFINITION  ${recordId}.`,
    `ACCESSION   ${recordId}`,
    'ORIGIN',
    '        1 atgcatgcatgc',
    '//',
    ''
  ].join('\n');
  const bytes = new TextEncoder().encode(text);
  return {
    live: {
      name,
      type: 'application/genbank',
      size: bytes.byteLength,
      lastModified: 1,
      async arrayBuffer() { return bytes.buffer.slice(0); }
    },
    serialized: {
      name,
      type: 'application/genbank',
      size: bytes.byteLength,
      lastModified: 1,
      encoding: 'base64',
      data: btoa(text)
    }
  };
};
const linearA = linearFile('record-a.gbk', 'RecordA');
const linearB = linearFile('record-b.gbk', 'RecordB');
state.mode.value = 'linear';
state.lInputType.value = 'gb';
state.results.value = [];
state.featureCatalog.value = null;
state.form.normalize_length = false;
state.linearSeqs.splice(0, state.linearSeqs.length,
  {
    uid: 'record-a',
    gb: linearA.live,
    gff: null,
    fasta: null,
    depth: null,
    losat_gencode: 1,
    definition: '',
    record_subtitle: '',
    region_record_id: '',
    region_start: null,
    region_end: null,
    region_reverse: false
  },
  {
    uid: 'record-b',
    gb: linearB.live,
    gff: null,
    fasta: null,
    depth: null,
    losat_gencode: 1,
    definition: '',
    record_subtitle: '',
    region_record_id: '',
    region_start: null,
    region_end: null,
    region_reverse: false
  }
);
state.linearRecordLayoutEnabled.value = true;
state.linearRecordRows.splice(0, state.linearRecordRows.length,
  { uid: 'record-a', row: 2 },
  { uid: 'record-b', row: 1 }
);
state.linearComparisonPlan.mode = 'none';
state.linearComparisonPlan.defaultSource = 'losat';
state.linearComparisonPlan.edges.splice(0, state.linearComparisonPlan.edges.length);
const linearCommitted = buildCanonicalRenderRequest({
  state,
  filesData: {
    linearSeqs: [
      { ...state.linearSeqs[0], gb: linearA.serialized },
      { ...state.linearSeqs[1], gb: linearB.serialized }
    ],
    linearComparisons: []
  },
  comparisonPlanSnapshot: resolveLinearComparisonPlan({
    plan: state.linearComparisonPlan,
    sequences: state.linearSeqs,
    layout: state.linearRecordRows,
    losatProgram: state.losatProgram.value,
    blastpMode: state.losat.blastp.mode
  })
});
adoptCanonicalRenderArtifacts(linearCommitted, { adoptOwnedRequest: true });
linearCommitted.renderRequest.records.forEach((record) => {
  record.selector = null;
  record.region = null;
});
linearCommitted.renderRequest.layout.multiRecordPositions = [
  'RecordB@1',
  'RecordA@2'
];

const linearSaved = await exportSession('source-record-id-layout');
assert.equal(linearSaved.status, 'saved');
const linearSession = JSON.parse(gunzipSync(
  Buffer.from(await linearSaved.blob.arrayBuffer())
).toString('utf8'));
assert.deepEqual(linearSession.renderRequest.layout.multiRecordPositions, [
  'RecordB@1',
  'RecordA@2'
]);
assert.deepEqual(linearSession.config.linearRecordLayout.rows, [
  { uid: 'record-a', row: 2 },
  { uid: 'record-b', row: 1 }
]);
