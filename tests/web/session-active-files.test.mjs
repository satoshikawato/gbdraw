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

const { serializeActiveRenderFiles } = await import(
  '../../gbdraw/web/js/services/config.js'
);

const readableFile = (name, text) => {
  let reads = 0;
  const bytes = new TextEncoder().encode(text);
  return {
    name,
    type: 'text/plain',
    size: bytes.byteLength,
    lastModified: 1,
    arrayBuffer: async () => {
      reads += 1;
      return bytes.slice().buffer;
    },
    reads: () => reads
  };
};
const unreadableFile = (name) => ({
  name,
  type: 'text/plain',
  size: 1,
  lastModified: 1,
  arrayBuffer: async () => {
    throw new Error(`inactive file was read: ${name}`);
  }
});

const circularInput = readableFile('active.gb', 'LOCUS active');
const retainedDepth = unreadableFile('retained-depth.tsv');
const inactiveLinear = unreadableFile('inactive-linear.gb');
const sharedDraftFile = unreadableFile('inactive-colors.tsv');
const sourceState = {
  form: { show_depth: false },
  adv: {
    circular_track_slots_enabled: false,
    circular_track_slots: [],
    linear_track_slots_enabled: false,
    linear_track_slots: []
  },
  circularConservation: { enabled: false, source: 'upload' },
  files: {
    c_gb: circularInput,
    c_gff: null,
    c_fasta: null,
    c_depth: retainedDepth,
    c_conservation_blasts: [],
    c_conservation_blasts_source: null,
    c_conservation_fastas: [],
    c_conservation_sequence_sources: [],
    d_color: sharedDraftFile,
    t_color: null,
    blacklist: null,
    whitelist: null,
    qualifier_priority: null,
    linearCanonicalComparisons: []
  },
  linearSeqs: [{
    uid: 'inactive',
    gb: inactiveLinear,
    gff: null,
    fasta: null,
    depth: null,
    blast: null
  }],
  linearComparisons: []
};

const circular = await serializeActiveRenderFiles('circular', sourceState);
assert.equal(circularInput.reads(), 1);
assert.equal(circular.c_gb.name, 'active.gb');
assert.equal(circular.c_depth, null);
assert.deepEqual(circular.linearSeqs, []);
assert.equal(circular.d_color, null);

const linearInput = readableFile('linear.gb', 'LOCUS linear');
const linearState = {
  ...sourceState,
  files: {
    ...sourceState.files,
    c_gb: unreadableFile('inactive-circular.gb')
  },
  linearSeqs: [{
    uid: 'linear-uid',
    gb: linearInput,
    gff: null,
    fasta: null,
    depth: retainedDepth,
    blast: null
  }]
};
const linear = await serializeActiveRenderFiles('linear', linearState);
assert.equal(linearInput.reads(), 1);
assert.equal(linear.linearSeqs[0].gb.name, 'linear.gb');
assert.equal(linear.linearSeqs[0].depth, null);
assert.equal(linear.c_gb, null);
