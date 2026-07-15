import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourceRoot = join(repoRoot, 'gbdraw', 'web', 'js');
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-session-request-'));
await cp(sourceRoot, join(tempRoot, 'js'), { recursive: true });
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');

const { buildCanonicalSessionRequest, projectCanonicalSessionRequest } = await import(
  pathToFileURL(join(tempRoot, 'js', 'services', 'session-request.js'))
);

const ref = (value) => ({ value });
const state = {
  mode: ref('circular'),
  cInputType: ref('gb'),
  lInputType: ref('gb'),
  circularRecordList: ref([]),
  form: {
    prefix: 'web-session', species: '', strain: '', plot_title: '', legend: 'right',
    multi_record_canvas: false, suppress_gc: false, suppress_skew: false,
    show_gc: false, show_skew: false, show_depth: false, separate_strands: true,
    labels_mode: 'none', show_labels_linear: 'none', track_type: 'tuckin',
    linear_track_layout: 'middle', linear_ruler_on_axis: false, align_center: false,
    keep_definition_left_aligned: false, scale_style: 'bar', normalize_length: false
  },
  adv: {
    features: ['CDS'], feature_shapes: { CDS: 'arrow' }, nt: 'GC', evalue: '1e-5',
    min_bitscore: 50, identity: 70, alignment_length: 0, plot_title_position: 'none',
    gc_content_mode: 'deviation', gc_content_show_axis: true, gc_content_show_ticks: true,
    depth_color: '#4A90E2', depth_show_axis: true, depth_show_ticks: true,
    pairwise_match_style: 'ribbon', multi_record_size_mode: 'auto',
    multi_record_min_radius_ratio: 0.55, multi_record_column_gap_ratio: 0.10,
    multi_record_row_gap_ratio: 0.05, multi_record_positions: [], depth_tracks: [],
    circular_track_slots_enabled: false, circular_track_slots_axis_index: null,
    circular_track_slots: [], linear_track_slots_enabled: false,
    linear_track_slots_axis_index: null, linear_track_slots: []
  },
  normalizePaletteColors: (value) => value,
  paletteDefinitions: ref({ default: {} }),
  currentColors: ref({}),
  selectedPalette: ref('default'),
  manualSpecificRules: [],
  featureVisibilityRules: ref([]),
  filterMode: ref('None'),
  manualBlacklist: ref(''),
  manualWhitelist: [],
  manualPriorityRules: [],
  labelTextFeatureOverrides: {},
  labelTextBulkOverrides: {},
  labelTextFeatureOverrideSources: {},
  labelVisibilityOverrides: {},
  editableLabels: ref([]),
  extractedFeatures: ref([]),
  circularConservation: { reference: 'auto', labels: '', series: [] },
  blastSource: ref('files'),
  losatProgram: ref('blastn'),
  losat: { blastp: {} },
  selectedOrthogroupAlignmentFeature: ref('')
};
const genbankText = `LOCUS       WEBTEST                    4 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Web canonical session fixture.
ACCESSION   WEBTEST
VERSION     WEBTEST
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
ORIGIN
        1 atgc
//
`;
const genbank = {
  name: 'input.gb',
  type: 'text/plain',
  size: new TextEncoder().encode(genbankText).byteLength,
  lastModified: 0,
  data: btoa(genbankText)
};
const filesData = { c_gb: genbank, linearSeqs: [] };

const canonical = buildCanonicalSessionRequest({ state, filesData });
assert.equal(canonical.renderRequest.schema, 1);
assert.equal(canonical.renderRequest.mode, 'circular');
assert.equal(canonical.renderRequest.records[0].source.resourceId, 'record-1-genbank');
assert.equal(canonical.resources['record-1-genbank'].kind, 'genbank');
assert.equal(canonical.resources['record-1-genbank'].encoding, 'base64');
assert.equal(canonical.renderRequest.output.prefix, 'web-session');

const projection = projectCanonicalSessionRequest(canonical);
assert.equal(projection.mode, 'circular');
assert.equal(projection.inputType, 'gb');
assert.equal(projection.files.c_gb.data, genbank.data);
assert.equal(projection.config.form.prefix, 'web-session');

if (process.argv.includes('--print')) console.log(JSON.stringify(canonical));
