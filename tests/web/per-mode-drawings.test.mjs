// Per-mode drawings (plan OV-80 §8 "Per-mode values", §5 Depth series count):
// a value set in one mode's drawing stays out of the other drawing, and each
// mode's request carries its own drawing's values.
import assert from 'node:assert/strict';
import { test } from 'node:test';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};

const { state, createDefaultAdv, createDefaultForm } = await import('../../gbdraw/web/js/state.js');
const { buildCanonicalRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
const { resolveLinearComparisonPlan } = await import('../../gbdraw/web/js/services/linear-comparisons.js');
const { MODE_SCOPED_SETTINGS } = await import('../../gbdraw/web/js/mode-scoped-settings.generated.js');
const MODE_SCOPED_ROWS = MODE_SCOPED_SETTINGS.rows;

const MODES = ['circular', 'linear'];
const genbankText = `LOCUS       PERMODE                    40 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Per-mode drawing fixture.
ACCESSION   PERMODE
VERSION     PERMODE
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             1..30
                     /product="test protein"
ORIGIN
        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//
`;
const file = (name, text) => ({
  name, type: 'text/plain', size: new TextEncoder().encode(text).byteLength, lastModified: 0, data: btoa(text)
});
const genbank = () => file('per-mode.gb', genbankText);
const depthFile = (name) => file(name, 'position\tdepth\n1\t10\n2\t12\n');

const resetDrawing = (mode) => {
  const drawing = state.drawings[mode];
  Object.assign(drawing.form, createDefaultForm());
  Object.assign(drawing.adv, createDefaultAdv(mode));
  drawing.linearRecordLayoutEnabled.value = false;
  drawing.linearComparisonPlan.mode = 'none';
  drawing.linearComparisonPlan.edges.splice(0);
};
// Generate builds the request of the shown mode from that mode's drawing.
const requestOf = (mode, { linearDepth = null } = {}) => {
  const drawing = state.drawings[mode];
  state.mode.value = mode;
  state.cInputType.value = 'gb';
  state.lInputType.value = 'gb';
  state.circularRecordList.value = [];
  const filesData = mode === 'circular'
    ? { c_gb: genbank(), linearSeqs: [] }
    : {
        linearSeqs: [{
          uid: 'record-1', gb: genbank(), depth: linearDepth, losat_gencode: 1, region_record_id: '',
          region_start: null, region_end: null, region_reverse: false
        }],
        linearComparisons: []
      };
  const comparisonPlanSnapshot = mode === 'linear'
    ? resolveLinearComparisonPlan({
        plan: drawing.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
        losatProgram: drawing.losatProgram.value, blastpMode: drawing.losat.blastp.mode
      })
    : null;
  return buildCanonicalRenderRequest({ state, drawing, filesData, comparisonPlanSnapshot });
};
const SAMPLE_PATHS = {
  circular: ['labels.font_size.long', 'objects.axis.circular.stroke_width.long', 'canvas.strandedness',
    'objects.gc_content.mode', 'objects.depth.min_depth'],
  linear: ['labels.font_size.linear.long', 'objects.axis.linear.stroke_width.long', 'canvas.strandedness',
    'objects.gc_content.mode', 'objects.depth.min_depth']
};
const samples = (mode, canonical) => SAMPLE_PATHS[mode]
  .map((path) => canonical.renderRequest.diagramOptions.configOverrides[path]);

test('no object of one drawing is shared with the other drawing', () => {
  // A shared nested object (a track slot list, a Legend map) would carry an
  // edit across modes; every registry row then holds a value per mode.
  const seen = new Map();
  const walk = (value, path, mode) => {
    if (!value || typeof value !== 'object') return;
    const raw = Object.hasOwn(value, 'value') && Object.keys(value).length === 1 ? value.value : value;
    if (!raw || typeof raw !== 'object' || Object.isFrozen(raw)) return;
    const other = seen.get(raw);
    assert.ok(!other || other.mode === mode, `${path} (${mode}) is ${other?.path} (${other?.mode})`);
    if (other) return;
    seen.set(raw, { path, mode });
    for (const [key, child] of Object.entries(raw)) walk(child, `${path}.${key}`, mode);
  };
  for (const mode of MODES) {
    for (const [key, member] of Object.entries(state.drawings[mode])) walk(member, key, mode);
  }
  assert.ok(MODE_SCOPED_ROWS.length > 100, 'the registry names the per-mode values');
});

test('each mode request carries its own drawing values', () => {
  MODES.forEach(resetDrawing);
  const linearDefaults = samples('linear', requestOf('linear'));
  const circular = state.drawings.circular;
  circular.adv.label_font_size = 13;
  circular.adv.axis_stroke_width = 3;
  circular.form.separate_strands = false;
  circular.adv.gc_content_mode = 'percent';
  circular.adv.depth_min = 4;
  assert.deepEqual(samples('circular', requestOf('circular')), [13, 3, false, 'percent', 4]);
  assert.deepEqual(samples('linear', requestOf('linear')), linearDefaults, 'the Linear request keeps the Linear values');
  const linear = state.drawings.linear;
  linear.adv.label_font_size = 7;
  linear.adv.depth_min = 2;
  linear.adv.gc_content_mode = 'deviation';
  const fromLinear = samples('linear', requestOf('linear'));
  assert.deepEqual([fromLinear[0], fromLinear[3], fromLinear[4]], [7, 'deviation', 2]);
  assert.deepEqual(samples('circular', requestOf('circular')), [13, 3, false, 'percent', 4]);
  state.mode.value = 'circular';
  MODES.forEach(resetDrawing);
});

test('two Circular Depth series and one Linear Depth file both build their requests (Depth series count)', () => {
  MODES.forEach(resetDrawing);
  // Two series configured while Circular had two Depth files.
  state.drawings.circular.form.show_depth = true;
  state.drawings.circular.adv.depth_tracks.splice(0, state.drawings.circular.adv.depth_tracks.length,
    { label: 'Illumina', color: '#112233' }, { label: 'Nanopore', color: '#445566' });
  state.drawings.linear.form.show_depth = true;
  state.mode.value = 'linear';
  const linear = requestOf('linear', { linearDepth: [depthFile('a.depth.tsv')] });
  assert.equal(linear.renderRequest.diagramOptions.depthTracks?.length, 1);
  assert.equal(state.drawings.circular.adv.depth_tracks.length, 2, 'the Linear request leaves the Circular series');
  state.mode.value = 'circular';
  MODES.forEach(resetDrawing);
});
