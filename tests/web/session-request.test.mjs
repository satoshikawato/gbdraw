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
const {
  buildLinearTrackSlotSpec,
  linearTrackAxisIndexForEnabledSlots,
  migrateLinearTrackSlotsToCurrentSchema,
  parseLinearTrackSlotSpec,
  parseLinearTrackSlotSpecs
} = await import(
  pathToFileURL(join(tempRoot, 'js', 'app', 'linear-track-slots.js'))
);
const { validateTrackSlotBindingInvariants } = await import(
  pathToFileURL(join(tempRoot, 'js', 'app', 'track-slot-validation.js'))
);

assert.equal(linearTrackAxisIndexForEnabledSlots([
  { id: 'above', renderer: 'depth', enabled: true, params: { track_index: 0 } },
  { id: 'removed', renderer: 'depth', enabled: false, params: {} },
  { id: 'features', renderer: 'features', enabled: true, params: {} }
], 2), 1);

const legacyLinearSlots = [
  {
    id: 'features', renderer: 'features', enabled: false, side: 'overlay',
    height: '48px', spacing: '7px', z: 3, params: { legend_label: 'Genes' }
  },
  {
    id: 'gc', renderer: 'dinucleotide_content', enabled: true, side: 'below',
    height: '22px', spacing: '4px', z: 1, params: { nt: 'AT' }
  }
];
const migratedLinearSlots = migrateLinearTrackSlotsToCurrentSchema(legacyLinearSlots, 1);
assert.equal(migratedLinearSlots[0].height, undefined);
assert.equal(migratedLinearSlots[0].spacing, undefined);
assert.equal(migratedLinearSlots[0].enabled, false);
assert.equal(migratedLinearSlots[0].side, 'overlay');
assert.equal(migratedLinearSlots[0].z, 3);
assert.deepEqual(migratedLinearSlots[0].params, { legend_label: 'Genes' });
assert.equal(migratedLinearSlots[1].height, '22px');
assert.equal(migratedLinearSlots[1].spacing, '4px');
assert.equal(legacyLinearSlots[0].height, '48px', 'migration must not mutate imported session data');
const currentLinearSlots = migrateLinearTrackSlotsToCurrentSchema(legacyLinearSlots, 2);
assert.equal(currentLinearSlots[0].height, '48px');
assert.equal(currentLinearSlots[0].spacing, '7px');
assert.throws(
  () => migrateLinearTrackSlotsToCurrentSchema(legacyLinearSlots, '1'),
  /Unsupported linear track slot schema version/
);

const parsedLinearSlot = parseLinearTrackSlotSpec(
  'genes:features@side=above,h=31px,spacing=5px,z=2,legend_label=Genes'
);
assert.equal(parsedLinearSlot.id, 'genes');
assert.equal(parsedLinearSlot.renderer, 'features');
assert.equal(parsedLinearSlot.height, '31px');
assert.equal(parsedLinearSlot.spacing, '5px');
assert.equal(parsedLinearSlot.side, 'above');
assert.equal(parsedLinearSlot.z, 2);
assert.equal(parsedLinearSlot.params.legend_label, 'Genes');
assert.equal(buildLinearTrackSlotSpec(parsedLinearSlot), 'genes:features@side=above,h=31px,spacing=5px,z=2,legend_label=Genes');

const aliasedLinearSlot = parseLinearTrackSlotSpec(
  'depth:depth@ID=coverage,TYPE=features,SHOW=no,HEIGHT=40px,spacing=5,Z_INDEX=4'
);
assert.equal(aliasedLinearSlot.id, 'coverage');
assert.equal(aliasedLinearSlot.renderer, 'features');
assert.equal(aliasedLinearSlot.enabled, false);
assert.equal(aliasedLinearSlot.height, '40px');
assert.equal(aliasedLinearSlot.spacing, '5px');
assert.equal(aliasedLinearSlot.z, 4);

const structuredLinearSlot = parseLinearTrackSlotSpec({
  kind: 'linearTrackSlot',
  id: 'depth_2',
  renderer: 'depth',
  enabled: true,
  side: 'below',
  height: { value: 12.5, unit: 'px' },
  spacing: { value: 3, unit: 'px' },
  z: 2,
  params: { track_index: 1, legend_label: 'Depth 2' }
});
assert.deepEqual(structuredLinearSlot, {
  id: 'depth_2',
  renderer: 'depth',
  enabled: true,
  side: 'below',
  height: '12.5px',
  spacing: '3px',
  z: 2,
  params: { track_index: 1, legend_label: 'Depth 2' }
});
assert.throws(() => parseLinearTrackSlotSpec('missing-renderer'), /requires '<slot_id>:<renderer>'/);
assert.throws(() => parseLinearTrackSlotSpec('mystery:not_a_renderer'), /Unsupported linear track renderer/);
assert.throws(() => parseLinearTrackSlotSpec('features:features@spacing='), /Invalid linear track slot spacing/);
assert.throws(
  () => parseLinearTrackSlotSpec({
    kind: 'linearTrackSlot', id: 'bad', renderer: 'features', enabled: true,
    side: 'overlay', height: { value: 1, unit: 'factor' }, spacing: null, z: 0, params: {}
  }),
  /only accepts px/
);
assert.throws(
  () => parseLinearTrackSlotSpec({
    kind: 'linearTrackSlot', id: 'bad', renderer: 'features', enabled: true,
    side: 'overlay', height: null, spacing: null, z: 0, params: {}, unknown: true
  }),
  /unsupported field unknown/
);
assert.throws(
  () => parseLinearTrackSlotSpecs(['same:features', 'same:depth']),
  /Duplicate canonical linear track slot id/
);
assert.throws(() => parseLinearTrackSlotSpecs([]), /cannot be empty/);
assert.throws(
  () => parseLinearTrackSlotSpec('depth:depth@track_index=1.0'),
  /non-negative integer/
);
assert.throws(
  () => validateTrackSlotBindingInvariants([
    {
      id: 'depth', renderer: 'depth', enabled: true, side: 'below', z: 0,
      params: { track_index: 1.5 }
    }
  ], {
    modeLabel: 'Linear',
    layoutKind: 'linear',
    supportedRenderers: ['features', 'depth', 'annotations', 'spacer'],
    supportedSides: ['above', 'below', 'overlay'],
    anchorlessRenderers: ['spacer'],
    depthTrackCount: 1
  }),
  /non-negative integer/
);
assert.throws(
  () => validateTrackSlotBindingInvariants([
    {
      id: 'depth', renderer: 'depth', enabled: true, side: 'below', z: 0,
      params: { track_index: 2 }
    }
  ], {
    modeLabel: 'Linear',
    layoutKind: 'linear',
    supportedRenderers: ['features', 'depth', 'annotations', 'spacer'],
    supportedSides: ['above', 'below', 'overlay'],
    anchorlessRenderers: ['spacer'],
    depthTrackCount: 1
  }),
  /available range is 0\.\.0/
);
assert.throws(
  () => validateTrackSlotBindingInvariants([
    {
      id: 'annotation', renderer: 'annotations', enabled: true,
      side: 'overlay', z: 2, params: { anchor_slot: 'missing', layer: 'foreground' }
    }
  ], {
    modeLabel: 'Linear',
    layoutKind: 'linear',
    supportedRenderers: ['features', 'depth', 'annotations', 'spacer'],
    supportedSides: ['above', 'below', 'overlay'],
    anchorlessRenderers: ['spacer']
  }),
  /unknown anchor_slot/
);
assert.throws(
  () => validateTrackSlotBindingInvariants([
    {
      id: 'depth', renderer: 'depth', enabled: true, side: 'overlay', z: 0,
      params: { track_index: 0 }
    }
  ], {
    modeLabel: 'Linear',
    layoutKind: 'linear',
    supportedRenderers: ['features', 'depth', 'annotations', 'spacer'],
    supportedSides: ['above', 'below', 'overlay'],
    anchorlessRenderers: ['spacer'],
    depthTrackCount: 1
  }),
  /only supported for features and annotations/
);
assert.throws(
  () => validateTrackSlotBindingInvariants([
    { id: 'features_a', renderer: 'features', enabled: true, side: 'overlay', z: 0, params: {} },
    { id: 'features_b', renderer: 'features', enabled: true, side: 'above', z: 0, params: {} }
  ], {
    modeLabel: 'Linear',
    layoutKind: 'linear',
    supportedRenderers: ['features', 'depth', 'annotations', 'spacer'],
    supportedSides: ['above', 'below', 'overlay'],
    anchorlessRenderers: ['spacer']
  }),
  /only one enabled features slot/
);
assert.throws(
  () => validateTrackSlotBindingInvariants([
    {
      id: 'annotation', renderer: 'annotations', enabled: true,
      side: 'above', z: 2, params: { anchor_slot: 'features' }
    },
    { id: 'features', renderer: 'features', enabled: true, side: 'overlay', z: 0, params: {} }
  ], {
    modeLabel: 'Linear',
    layoutKind: 'linear',
    supportedRenderers: ['features', 'depth', 'annotations', 'spacer'],
    supportedSides: ['above', 'below', 'overlay'],
    anchorlessRenderers: ['spacer']
  }),
  /uses anchor_slot without side=overlay/
);
assert.throws(
  () => validateTrackSlotBindingInvariants([
    {
      id: 'features', renderer: 'features', enabled: true, side: 'overlay',
      height: '-1px', z: 0, params: {}
    }
  ], {
    modeLabel: 'Linear',
    layoutKind: 'linear',
    supportedRenderers: ['features', 'depth', 'annotations', 'spacer'],
    supportedSides: ['above', 'below', 'overlay'],
    anchorlessRenderers: ['spacer']
  }),
  /height must be positive/
);
assert.throws(
  () => validateTrackSlotBindingInvariants([
    {
      id: 'depth', renderer: 'depth', enabled: true, side: 'inside',
      spacing: '0.1', inner_gap_px: '2', z: 0, params: { track_index: 0 }
    }
  ], {
    modeLabel: 'Circular',
    layoutKind: 'circular',
    supportedRenderers: ['features', 'ticks', 'depth', 'annotations', 'spacer'],
    supportedSides: ['inside', 'outside', 'overlay'],
    anchorlessRenderers: ['ticks', 'spacer'],
    depthTrackCount: 1
  }),
  /cannot combine spacing/
);
assert.throws(
  () => validateTrackSlotBindingInvariants([
    { id: 'features', renderer: 'features', enabled: true, side: 'overlay', z: 0, params: null }
  ], {
    modeLabel: 'Linear',
    layoutKind: 'linear',
    supportedRenderers: ['features', 'depth', 'annotations', 'spacer'],
    supportedSides: ['above', 'below', 'overlay'],
    anchorlessRenderers: ['spacer']
  }),
  /params must be an object/
);
assert.throws(
  () => parseLinearTrackSlotSpec({
    kind: 'linearTrackSlot', id: 'depth', renderer: 'depth', enabled: true,
    side: 'below', height: null, spacing: null, z: 0, params: { track_index: 1.5 }
  }),
  /non-negative integer/
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

state.paletteDefinitions.value = { default: { CDS: '#cccccc' } };
state.currentColors.value = { CDS: '#112233' };
state.manualSpecificRules.push(
  { feat: 'CDS', qual: 'gene', val: '^alpha$', color: '#445566', cap: 'Alpha', fromFile: false },
  { feat: 'CDS', qual: 'product', val: '^beta$', color: '#445566', cap: 'Alpha', fromFile: true }
);
state.filterMode.value = 'Whitelist';
state.manualWhitelist.push({ feat: 'CDS', qual: 'gene', key: 'alpha' });
state.manualPriorityRules.push({ feat: 'CDS', order: 'gene,product' });
state.featureVisibilityRules.value = [{
  id: 'visibility-1', source: 'manual', recordId: '*', featureType: 'CDS',
  qualifier: 'gene', value: '^alpha$', action: 'off'
}];
state.labelTextFeatureOverrides.f1 = 'Renamed alpha';
const semanticCanonical = buildCanonicalSessionRequest({ state, filesData });
const semanticProjection = projectCanonicalSessionRequest(semanticCanonical);
assert.deepEqual(semanticProjection.config.colors, { CDS: '#112233' });
assert.equal(semanticProjection.config.colorsAreOverrides, true);
assert.equal(semanticProjection.config.palette, 'default');
assert.equal(semanticProjection.config.rules.length, 2);
assert.equal(semanticProjection.config.rules.some((rule) => rule.fromFile), false);
assert.deepEqual(semanticProjection.config.whitelist, [
  { feat: 'CDS', qual: 'gene', key: 'alpha' }
]);
assert.deepEqual(semanticProjection.config.qualifierPriorityRules, [
  { feat: 'CDS', order: 'gene,product' }
]);
assert.equal(semanticProjection.config.filterMode, 'Whitelist');
assert.equal(semanticProjection.semanticFeatureState.featureVisibilityManualRules.length, 1);
assert.equal(semanticProjection.semanticFeatureState.labelOverrideRows.length, 1);
assert.equal(semanticProjection.semanticFeatureState.labelOverrideRows[0].labelText, 'Renamed alpha');
state.paletteDefinitions.value = { default: {} };
state.currentColors.value = {};
state.manualSpecificRules.splice(0);
state.filterMode.value = 'None';
state.manualWhitelist.splice(0);
state.manualPriorityRules.splice(0);
state.featureVisibilityRules.value = [];
delete state.labelTextFeatureOverrides.f1;

const invalidCircularDepthIndex = structuredClone(canonical);
invalidCircularDepthIndex.renderRequest.diagramOptions.tracks = {
  circularTrackSlots: ['depth:depth@track_index=1.5']
};
assert.throws(
  () => projectCanonicalSessionRequest(invalidCircularDepthIndex),
  /non-negative integer/
);

const invalidCircularOverlay = structuredClone(canonical);
invalidCircularOverlay.renderRequest.diagramOptions.tracks = {
  circularTrackSlots: [
    'review:annotations@side=overlay,set_id=review,anchor_slot=missing,layer=foreground,z=1'
  ]
};
assert.throws(
  () => projectCanonicalSessionRequest(invalidCircularOverlay),
  /references unknown anchor_slot/
);

const emptyDepthColumn = structuredClone(canonical);
emptyDepthColumn.renderRequest.diagramOptions.depthTrackFiles = [[null]];
emptyDepthColumn.renderRequest.diagramOptions.depthTrackLabels = ['Empty'];
assert.throws(
  () => projectCanonicalSessionRequest(emptyDepthColumn),
  /logical track index 0.*no source/
);

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
    label: 'Window', mark: 'band', lane: null, style: null, legendLabel: null,
    metadata: { _gbdraw_web_target_record_key: 'linear-source::#2' }
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
assert.equal(
  projectCanonicalSessionRequest(annotationCanonical).config.annotationSets[0].annotations[0]
    .metadata._gbdraw_web_target_record_key,
  'linear-source::#2'
);
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
state.adv.comparison_height = 42.5;
const linearFilesData = {
  linearSeqs: [
    { uid: 'first', gb: genbank, region_record_id: '', region_start: null, region_end: null, region_reverse: false },
    { uid: 'second', gb: genbank, region_record_id: 'RecA', region_start: null, region_end: null, region_reverse: false },
    { uid: 'third', gb: genbank, region_record_id: '#2', region_start: 10, region_end: 20, region_reverse: true }
  ],
  linearComparisons: []
};
const linearCanonical = buildCanonicalSessionRequest({ state, filesData: linearFilesData });
assert.equal(linearCanonical.renderRequest.diagramOptions.configOverrides.comparison_height, 42.5);
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
assert.equal(linearProjection.config.adv.comparison_height, 42.5);
assert.deepEqual(
  linearProjection.files.linearSeqs.map((seq) => seq.region_record_id),
  ['', 'RecA', '#2']
);
assert.equal(linearProjection.files.linearSeqs[2].region_start, 10);
assert.equal(linearProjection.files.linearSeqs[2].region_end, 20);
assert.equal(linearProjection.files.linearSeqs[2].region_reverse, true);

state.adv.comparison_height = null;
const autoHeightCanonical = buildCanonicalSessionRequest({ state, filesData: linearFilesData });
assert.equal(autoHeightCanonical.renderRequest.diagramOptions.configOverrides.comparison_height, null);
assert.equal(projectCanonicalSessionRequest(autoHeightCanonical).config.adv.comparison_height, null);
state.adv.comparison_height = -2;
assert.throws(
  () => buildCanonicalSessionRequest({ state, filesData: linearFilesData }),
  /Pairwise Match Height must be Auto or a positive finite number/
);
const historicalInvalidHeight = structuredClone(linearCanonical);
historicalInvalidHeight.renderRequest.diagramOptions.configOverrides.comparison_height = -2;
assert.throws(
  () => projectCanonicalSessionRequest(historicalInvalidHeight),
  /Pairwise Match Height must be Auto or a positive finite number/
);
assert.equal(
  projectCanonicalSessionRequest({
    ...historicalInvalidHeight,
    repairInvalidComparisonHeight: true
  }).config.adv.comparison_height,
  null
);
state.adv.comparison_height = 42.5;

const depthA = {
  ...genbank,
  name: 'sample-a.depth.tsv',
  type: 'text/tab-separated-values',
  data: btoa('position\tdepth\n1\t10\n')
};
const depthB = {
  ...genbank,
  name: 'sample-b.depth.tsv',
  type: 'text/tab-separated-values',
  data: btoa('position\tdepth\n1\t20\n')
};
state.form.show_depth = true;
state.adv.resolve_overlaps = true;
state.adv.feature_height = 9;
state.adv.depth_tracks = [
  { label: 'Sample A', color: '#112233', height: 12 },
  { label: 'Sample B', color: '#445566', height: 18 }
];

const combinedGenbankText = `${genbankText}${genbankText.replaceAll('WEBTEST', 'WEBTWO')}`;
const combinedGenbank = {
  ...genbank,
  name: 'combined.gb',
  size: new TextEncoder().encode(combinedGenbankText).byteLength,
  data: btoa(combinedGenbankText)
};
state.mode.value = 'circular';
state.form.multi_record_canvas = true;
state.circularRecordList.value = [
  { selector: '#1' },
  { selector: '#2' }
];
const circularSparseFilesData = {
  c_gb: combinedGenbank,
  c_depth: [
    [depthA, null],
    [null, depthB]
  ],
  linearSeqs: []
};
const circularSparseCanonical = buildCanonicalSessionRequest({
  state,
  filesData: circularSparseFilesData
});
assert.equal(circularSparseCanonical.renderRequest.records.length, 2);
assert.deepEqual(
  circularSparseCanonical.renderRequest.diagramOptions.depthTrackFiles.map((row) => (
    row.map((entry) => Boolean(entry?.resourceId))
  )),
  [[true, false], [false, true]]
);
const circularSparseProjection = projectCanonicalSessionRequest(circularSparseCanonical);
assert.deepEqual(
  circularSparseProjection.files.c_depth.map((row) => row.map(Boolean)),
  [[true, false], [false, true]]
);
assert.equal(circularSparseProjection.files.c_depth[0][0].data, depthA.data);
assert.equal(circularSparseProjection.files.c_depth[0][1], null);
assert.equal(circularSparseProjection.files.c_depth[1][0], null);
assert.equal(circularSparseProjection.files.c_depth[1][1].data, depthB.data);
assert.deepEqual(circularSparseProjection.config.adv.depth_tracks.map((track) => track.label), [
  'Sample A',
  'Sample B'
]);
const circularSparseRebuilt = buildCanonicalSessionRequest({
  state,
  filesData: circularSparseProjection.files
});
assert.deepEqual(
  circularSparseRebuilt.renderRequest.diagramOptions.depthTrackFiles.map((row) => (
    row.map((entry) => Boolean(entry?.resourceId))
  )),
  [[true, false], [false, true]]
);

const circularLegacyFlatCanonical = buildCanonicalSessionRequest({
  state,
  filesData: {
    c_gb: combinedGenbank,
    c_depth: [depthA, depthB],
    linearSeqs: []
  }
});
assert.deepEqual(
  circularLegacyFlatCanonical.renderRequest.diagramOptions.depthTrackFiles.map((row) => (
    row.map((entry) => Boolean(entry?.resourceId))
  )),
  [[true, true], [true, true]]
);
assert.throws(
  () => buildCanonicalSessionRequest({
    state,
    filesData: {
      c_gb: combinedGenbank,
      c_depth: [[depthA, depthB]],
      linearSeqs: []
    }
  }),
  /Circular Depth matrix has 1 record rows; expected 2/
);

state.mode.value = 'linear';
state.form.multi_record_canvas = false;
state.circularRecordList.value = [];
state.adv.linear_track_slots_enabled = true;
state.adv.linear_track_slots_axis_index = 1;
state.adv.linear_track_slots = [
  {
    id: 'depth_b', renderer: 'depth', enabled: true, side: 'above',
    height: '18px', spacing: '2px', params: { track_index: 1 }
  },
  {
    id: 'features', renderer: 'features', enabled: true, side: 'overlay',
    height: '44px', spacing: '6px', params: {}
  },
  { id: 'depth_a', renderer: 'depth', enabled: true, side: 'below', params: { track_index: 0 } }
];
const sparseDepthFilesData = {
  ...linearFilesData,
  linearSeqs: linearFilesData.linearSeqs.slice(0, 2).map((seq, index) => ({
    ...seq,
    depth: index === 0 ? [depthA, null] : [null, depthB]
  }))
};
const sparseDepthCanonical = buildCanonicalSessionRequest({ state, filesData: sparseDepthFilesData });
assert.equal(sparseDepthCanonical.renderRequest.diagramOptions.depthTrackFiles.length, 2);
assert.equal(sparseDepthCanonical.renderRequest.diagramOptions.depthTrackFiles[0].length, 2);
assert.equal(sparseDepthCanonical.renderRequest.diagramOptions.depthTrackFiles[0][1], null);
assert.equal(sparseDepthCanonical.renderRequest.diagramOptions.depthTrackFiles[1][0], null);
assert.deepEqual(
  sparseDepthCanonical.renderRequest.diagramOptions.depthTrackLabels,
  ['Sample A', 'Sample B']
);
assert.equal(sparseDepthCanonical.renderRequest.diagramOptions.tracks.linearTrackAxisIndex, 1);
assert.deepEqual(
  sparseDepthCanonical.renderRequest.diagramOptions.tracks.linearTrackSlots.map((slot) => slot.split(':')[0]),
  ['depth_b', 'features', 'depth_a']
);
const sparseDepthProjection = projectCanonicalSessionRequest(sparseDepthCanonical);
assert.equal(sparseDepthProjection.files.linearSeqs[0].depth.length, 2);
assert.equal(sparseDepthProjection.files.linearSeqs[0].depth[0].data, depthA.data);
assert.equal(sparseDepthProjection.files.linearSeqs[0].depth[1], null);
assert.equal(sparseDepthProjection.files.linearSeqs[1].depth[0], null);
assert.equal(sparseDepthProjection.files.linearSeqs[1].depth[1].data, depthB.data);
assert.equal(sparseDepthProjection.config.adv.linear_track_slots_enabled, true);
assert.equal(sparseDepthProjection.config.adv.linear_track_slots_schema_version, 2);
assert.equal(sparseDepthProjection.config.adv.linear_track_slots_axis_index, 1);
assert.deepEqual(
  sparseDepthProjection.config.adv.linear_track_slots.map((slot) => slot.id),
  ['depth_b', 'features', 'depth_a']
);
assert.equal(sparseDepthProjection.config.adv.linear_track_slots[0].height, '18px');
assert.equal(sparseDepthProjection.config.adv.linear_track_slots[0].spacing, '2px');
assert.equal(sparseDepthProjection.config.adv.linear_track_slots[1].height, '44px');
assert.equal(sparseDepthProjection.config.adv.linear_track_slots[1].spacing, '6px');
assert.equal(sparseDepthProjection.config.adv.resolve_overlaps, true);
assert.equal(sparseDepthProjection.config.adv.feature_height, 9);

const structuredSparseDepthCanonical = {
  ...sparseDepthCanonical,
  renderRequest: {
    ...sparseDepthCanonical.renderRequest,
    diagramOptions: {
      ...sparseDepthCanonical.renderRequest.diagramOptions,
      tracks: {
        ...sparseDepthCanonical.renderRequest.diagramOptions.tracks,
        linearTrackSlots: [
          {
            kind: 'linearTrackSlot', id: 'depth_b', renderer: 'depth', enabled: true,
            side: 'above', height: { value: 18, unit: 'px' }, spacing: { value: 2, unit: 'px' },
            z: 3, params: { track_index: 1, legend_label: 'Sample B' }
          },
          {
            kind: 'linearTrackSlot', id: 'features', renderer: 'features', enabled: true,
            side: 'overlay', height: { value: 44, unit: 'px' }, spacing: { value: 6, unit: 'px' },
            z: 1, params: {}
          },
          {
            kind: 'linearTrackSlot', id: 'depth_a', renderer: 'depth', enabled: true,
            side: 'below', height: null, spacing: null, z: 0,
            params: { track_index: 0, legend_label: 'Sample A' }
          }
        ]
      }
    }
  }
};
const structuredV33Projection = projectCanonicalSessionRequest(structuredSparseDepthCanonical);
assert.deepEqual(
  structuredV33Projection.config.adv.linear_track_slots.map((slot) => slot.id),
  ['depth_b', 'features', 'depth_a']
);
assert.equal(structuredV33Projection.config.adv.linear_track_slots[0].renderer, 'depth');
assert.equal(structuredV33Projection.config.adv.linear_track_slots[0].height, '18px');
assert.equal(structuredV33Projection.config.adv.linear_track_slots[0].spacing, '2px');
assert.equal(structuredV33Projection.config.adv.linear_track_slots[0].z, 3);
assert.deepEqual(
  structuredV33Projection.config.adv.linear_track_slots[0].params,
  { track_index: 1, legend_label: 'Sample B' }
);
assert.equal(structuredV33Projection.config.adv.linear_track_slots[1].height, '44px');
assert.equal(structuredV33Projection.config.adv.linear_track_slots[1].spacing, '6px');

const legacySparseDepthProjection = projectCanonicalSessionRequest({
  ...structuredSparseDepthCanonical,
  linearTrackSlotSchemaVersion: 1
});
assert.equal(legacySparseDepthProjection.config.adv.linear_track_slots[0].height, '18px');
assert.equal(legacySparseDepthProjection.config.adv.linear_track_slots[0].spacing, '2px');
assert.equal(legacySparseDepthProjection.config.adv.linear_track_slots[1].height, undefined);
assert.equal(legacySparseDepthProjection.config.adv.linear_track_slots[1].spacing, undefined);

state.adv.linear_track_slots_axis_index = 1;
state.adv.linear_track_slots = [
  {
    id: 'custom_removed', renderer: 'depth', enabled: false, side: 'above',
    depth_binding_error: 'Depth logical track index 0 is no longer available.', params: {}
  },
  { id: 'features', renderer: 'features', enabled: true, side: 'overlay', params: {} },
  { id: 'depth_b', renderer: 'depth', enabled: true, side: 'below', params: { track_index: 1 } }
];
const filteredAxisCanonical = buildCanonicalSessionRequest({ state, filesData: sparseDepthFilesData });
assert.equal(filteredAxisCanonical.renderRequest.diagramOptions.tracks.linearTrackAxisIndex, 0);
assert.deepEqual(
  filteredAxisCanonical.renderRequest.diagramOptions.tracks.linearTrackSlots.map((slot) => slot.split(':')[0]),
  ['features', 'depth_b']
);

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
    resources: session.resources,
    linearTrackSlotSchemaVersion: Number(session.version) <= 32 ? 1 : 2
  });
  assert.ok(['circular', 'linear'].includes(projectedSession.mode));
  assert.ok(
    projectedSession.files.c_gb || projectedSession.files.linearSeqs.length > 0,
    'projected session must restore at least one input record'
  );
}

if (process.argv.includes('--print')) console.log(JSON.stringify(canonical));
