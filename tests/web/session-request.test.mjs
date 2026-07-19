import assert from 'node:assert/strict';
import { cp, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { gunzipSync } from 'node:zlib';

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
  selectedOrthogroupAlignmentFeature: ref(''),
  linearRecordLayoutEnabled: ref(false),
  linearRecordGap: ref(24),
  linearRecordRows: [],
  linearComparisons: [],
  annotationSets: []
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
assert.equal(canonical.renderRequest.schema, 2);
assert.equal(canonical.renderRequest.mode, 'circular');
assert.equal(canonical.renderRequest.records[0].source.resourceId, 'record-1-genbank');
assert.equal(canonical.resources['record-1-genbank'].kind, 'genbank');
assert.equal(canonical.resources['record-1-genbank'].encoding, 'base64');
assert.equal(canonical.renderRequest.output.prefix, 'web-session');
assert.equal(canonical.renderRequest.diagramOptions.configOverrides.circular_label_placement, 'horizontal');

const blastTable = { ...genbank, name: 'hits.tsv', data: btoa('ref\tcmp\t99\t4\t0\t0\t1\t4\t4\t1\t1e-20\t50\n') };
const comparisonFasta = { ...genbank, name: 'comparison.fna', data: btoa('>cmp\nAACCGG\n') };
state.circularConservation.series = [{ label: 'Comparison', color: '#E15759' }];
const conservationCanonical = buildCanonicalSessionRequest({
  state,
  filesData: {
    ...filesData,
    c_conservation_blasts: [blastTable],
    c_conservation_sequence_sources: [comparisonFasta],
  },
});
assert.equal(conservationCanonical.renderRequest.diagramOptions.conservationSequenceSources, undefined);
assert.deepEqual(conservationCanonical.webFiles.conservationSequenceSources, ['conservation-sequence-sources-1']);
assert.equal(
  projectCanonicalSessionRequest(conservationCanonical).files.c_conservation_sequence_sources[0].data,
  comparisonFasta.data
);
state.circularConservation.series = [];

state.annotationSets = [{
  id: 'review',
  annotations: [{
    id: 'window',
    target: { kind: 'coordinateSpan', record: null, start: 1, end: 4, coordinateSpace: 'source', wrapsOrigin: false, outOfBounds: 'clip' },
    label: 'Window', mark: 'band', lane: null, style: null, legendLabel: null, metadata: {}
  }],
  defaultStyle: {
    stroke: '#404040', strokeWidth: 1.5, strokeDasharray: [], lineCap: 'tick', fill: null,
    fillOpacity: 0.2, hatch: null, labelColor: '#202020', labelFontSize: null,
    labelOrientation: 'auto', labelPosition: 'center', labelOffset: 4
  },
  legendLabel: null
}];
const annotationCanonical = buildCanonicalSessionRequest({ state, filesData });
assert.equal(annotationCanonical.renderRequest.diagramOptions.annotations.sets[0].id, 'review');
assert.equal(projectCanonicalSessionRequest(annotationCanonical).config.annotationSets[0].annotations[0].id, 'window');
state.annotationSets = [];

state.adv.circular_label_placement = 'radial';
const radialCanonical = buildCanonicalSessionRequest({ state, filesData });
assert.equal(radialCanonical.renderRequest.diagramOptions.configOverrides.circular_label_placement, 'radial');

const projection = projectCanonicalSessionRequest(canonical);
assert.equal(projection.mode, 'circular');
assert.equal(projection.inputType, 'gb');
assert.equal(projection.files.c_gb.data, genbank.data);
assert.equal(projection.config.form.prefix, 'web-session');
assert.equal(projection.config.adv.circular_label_placement, 'horizontal');
const radialProjection = projectCanonicalSessionRequest(radialCanonical);
assert.equal(radialProjection.config.adv.circular_label_placement, 'radial');

const customTrackProjection = projectCanonicalSessionRequest({
  renderRequest: {
    ...canonical.renderRequest,
    diagramOptions: {
      ...canonical.renderRequest.diagramOptions,
      tracks: {
        circularTrackSlots: [
          'features:features@side=overlay,lane_direction=split',
          'plastome_regions:annotations@set_id=plastome_regions,side=inside,r=0.65,w=20px,show_labels=true,padding_px=1,overflow=compress,inner_gap_px=1,outer_gap_px=1',
          'gc_content:dinucleotide_content@side=inside,r=0.56,w=0.08'
        ],
        circularTrackAxisIndex: 0,
        centerReservedRadius: 140
      }
    }
  },
  resources: canonical.resources
});
assert.equal(customTrackProjection.config.adv.circular_track_slots_enabled, true);
assert.equal(customTrackProjection.config.adv.circular_track_slots_axis_index, 0);
assert.equal(customTrackProjection.config.adv.center_reserved_radius, 140);
assert.equal(customTrackProjection.config.adv.circular_track_slots.length, 3);
assert.deepEqual(customTrackProjection.config.adv.circular_track_slots[1], {
  id: 'plastome_regions',
  renderer: 'annotations',
  enabled: true,
  width: '20px',
  radius: '0.65',
  spacing: null,
  inner_gap_px: '1',
  outer_gap_px: '1',
  side: 'inside',
  z: 0,
  params: {
    set_id: 'plastome_regions',
    show_labels: true,
    padding_px: '1',
    overflow: 'compress',
    layer: 'foreground'
  }
});

const secondGenbank = {
  ...genbank,
  name: 'second.gb',
  data: btoa(genbankText.replaceAll('WEBTEST', 'WEBTWO'))
};
const multiCircularProjection = projectCanonicalSessionRequest({
  renderRequest: {
    ...canonical.renderRequest,
    records: [
      canonical.renderRequest.records[0],
      { ...canonical.renderRequest.records[0], source: { kind: 'genbank', resourceId: 'record-2-genbank' } }
    ]
  },
  resources: {
    ...canonical.resources,
    'record-2-genbank': { kind: 'genbank', ...secondGenbank }
  }
});
const combinedCircularGenbank = atob(multiCircularProjection.files.c_gb.data);
assert.equal(combinedCircularGenbank.match(/^LOCUS/gm)?.length, 2);
assert.match(combinedCircularGenbank, /WEBTEST/);
assert.match(combinedCircularGenbank, /WEBTWO/);

state.mode.value = 'linear';
state.lInputType.value = 'gb';
const linearFilesData = {
  linearSeqs: [
    { uid: 'first', gb: genbank, region_record_id: '', region_start: null, region_end: null, region_reverse: false },
    { uid: 'second', gb: genbank, region_record_id: 'RecA', region_start: null, region_end: null, region_reverse: false },
    { uid: 'third', gb: genbank, region_record_id: '#2', region_start: 10, region_end: 20, region_reverse: true }
  ],
  linearComparisons: []
};
const linearCanonical = buildCanonicalSessionRequest({ state, filesData: linearFilesData });
assert.equal(linearCanonical.renderRequest.records[0].selector, null);
assert.deepEqual(linearCanonical.renderRequest.records[1].selector, { kind: 'recordId', value: 'RecA' });
assert.equal(linearCanonical.renderRequest.records[2].selector, null);
assert.deepEqual(linearCanonical.renderRequest.records[2].region, {
  selector: { kind: 'recordIndex', index: 1 },
  start: 10,
  end: 20,
  reverseComplement: true
});
assert.deepEqual(
  linearCanonical.renderRequest.records.map((record) => record.recordKey),
  ['first', 'second', 'third']
);

state.linearRecordLayoutEnabled.value = true;
state.linearRecordGap.value = 30;
state.linearRecordRows.splice(0, state.linearRecordRows.length,
  { uid: 'first', row: 1 }, { uid: 'second', row: 1 }, { uid: 'third', row: 2 });
const arrangedCanonical = buildCanonicalSessionRequest({ state, filesData: linearFilesData });
assert.deepEqual(arrangedCanonical.renderRequest.layout, {
  recordGapPx: 30,
  multiRecordPositions: ['#1@1', '#2@1', '#3@2']
});

state.losatProgram.value = 'blastp';
linearFilesData.linearComparisons = [{
  id: 'selected-losat-pair', queryUid: 'first', subjectUid: 'third', source: 'losat', file: null
}];
const losatPairCanonical = buildCanonicalSessionRequest({ state, filesData: linearFilesData });
const generatedProtein = losatPairCanonical.renderRequest.comparisons.find(
  (comparison) => comparison.kind === 'generatedProteinComparison'
);
assert.deepEqual(generatedProtein.pairs, [{ queryRecordIndex: 0, subjectRecordIndex: 2 }]);
assert.equal(
  projectCanonicalSessionRequest(losatPairCanonical).files.linearComparisons[0].source,
  'losat'
);

const linearProjection = projectCanonicalSessionRequest(linearCanonical);
assert.deepEqual(
  linearProjection.files.linearSeqs.map((seq) => seq.region_record_id),
  ['', 'RecA', '#2']
);
assert.equal(linearProjection.files.linearSeqs[2].region_start, 10);
assert.equal(linearProjection.files.linearSeqs[2].region_end, 20);
assert.equal(linearProjection.files.linearSeqs[2].region_reverse, true);

const projectSessionIndex = process.argv.indexOf('--project-session');
if (projectSessionIndex >= 0) {
  const sessionPath = process.argv[projectSessionIndex + 1];
  assert.ok(sessionPath, '--project-session requires a path');
  const sessionBytes = await readFile(sessionPath);
  const sessionText = sessionBytes[0] === 0x1f && sessionBytes[1] === 0x8b
    ? gunzipSync(sessionBytes).toString('utf8')
    : sessionBytes.toString('utf8');
  const session = JSON.parse(sessionText);
  const projectedSession = projectCanonicalSessionRequest({
    renderRequest: session.renderRequest,
    resources: session.resources
  });
  assert.ok(['circular', 'linear'].includes(projectedSession.mode));
  assert.ok(
    projectedSession.files.c_gb || projectedSession.files.linearSeqs.length > 0,
    'projected session must restore at least one input record'
  );
}

if (process.argv.includes('--print')) console.log(JSON.stringify(canonical));
