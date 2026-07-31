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

const { buildCanonicalRenderRequest, projectCanonicalSessionRequest } = await import(
  pathToFileURL(join(tempRoot, 'js', 'services', 'session-request.js'))
);
const {
  buildLinearTrackSlotPayload,
  buildLinearTrackSlotSpec,
  linearTrackAxisIndexForEnabledSlots,
  migrateLinearTrackSlotsToCurrentSchema,
  parseLinearTrackSlotSpec,
  parseLinearTrackSlotSpecs
} = await import(
  pathToFileURL(join(tempRoot, 'js', 'app', 'linear-track-slots.js'))
);
const {
  buildCircularTrackSlotPayload,
  buildCircularTrackSlotSpec,
  migrateLegacyCircularTrackSlot,
  migrateLegacyCircularTrackSlotSpec,
  parseCircularTrackSlotSpec,
  resolveCircularTrackFeaturePlacement
} = await import(
  pathToFileURL(join(tempRoot, 'js', 'app', 'circular-track-slots.js'))
);
const { validateTrackSlotBindingInvariants } = await import(
  pathToFileURL(join(tempRoot, 'js', 'app', 'track-slot-validation.js'))
);
const { orderedConservationSources } = await import(
  pathToFileURL(join(tempRoot, 'js', 'app', 'conservation-series.js'))
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
  /obsolete field 'spacing'/
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
    keep_definition_left_aligned: false, show_scale: true, scale_style: 'bar',
    normalize_length: false
  },
  adv: {
    features: ['CDS'], feature_shapes: { CDS: 'arrow' }, nt: 'GC', evalue: '1e-5',
    min_bitscore: 50, identity: 70, alignment_length: 0, plot_title_position: 'none',
    gc_content_mode: 'deviation', gc_content_show_axis: true, gc_content_show_ticks: true,
    depth_color: '#4A90E2', depth_show_axis: true, depth_show_ticks: true,
    depth_large_tick_interval: 25,
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
  losat: { blastp: { collinearMaxUnitGap: 2 } },
  selectedOrthogroupAlignmentFeature: ref(''),
  linearRecordLayoutEnabled: ref(false),
  linearRecordGap: ref(24),
  linearRecordRows: [],
  linearComparisons: [],
  annotationSets: []
};

const stateForCanonicalProjection = (projection) => {
  const config = projection.config || {};
  const linearLayout = config.linearRecordLayout || {};
  const palette = String(config.palette || 'default');
  return {
    ...state,
    mode: ref(projection.mode),
    cInputType: ref(projection.inputType),
    lInputType: ref(projection.inputType),
    form: { ...state.form, ...structuredClone(config.form || {}) },
    adv: { ...state.adv, ...structuredClone(config.adv || {}) },
    losat: structuredClone(config.losat || { threadsPerJob: 'auto', blastp: {} }),
    blastSource: ref(config.blastSource || 'files'),
    losatProgram: ref(config.losatProgram || 'blastn'),
    selectedOrthogroupAlignmentFeature: ref(
      projection.pipelineState?.selectedOrthogroupAlignmentFeature || ''
    ),
    currentColors: ref(structuredClone(config.colors || {})),
    selectedPalette: ref(palette),
    paletteDefinitions: ref({ [palette]: {} }),
    manualSpecificRules: structuredClone(config.rules || []),
    featureVisibilityRules: ref(structuredClone(
      projection.semanticFeatureState?.featureVisibilityManualRules || []
    )),
    filterMode: ref(config.filterMode || 'None'),
    manualBlacklist: ref(String(config.blacklistText || '')),
    manualWhitelist: structuredClone(config.whitelist || []),
    manualPriorityRules: structuredClone(config.qualifierPriorityRules || []),
    linearRecordLayoutEnabled: ref(Boolean(linearLayout.enabled)),
    linearRecordGap: ref(linearLayout.recordGap ?? 24),
    linearRecordRows: structuredClone(linearLayout.rows || []),
    linearComparisons: structuredClone(linearLayout.comparisons || []),
    annotationSets: structuredClone(config.annotationSets || [])
  };
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

state.form.multi_record_canvas = true;
const canonical = buildCanonicalRenderRequest({ state, filesData });
state.form.multi_record_canvas = false;
assert.equal(canonical.renderRequest.schema, 5);
assert.equal(canonical.renderRequest.mode, 'circular');
assert.equal(canonical.renderRequest.grouping, 'grid');
assert.equal(canonical.renderRequest.records[0].source.resourceId, 'record-1-genbank');
assert.equal(canonical.resources['record-1-genbank'].kind, 'genbank');
assert.equal(canonical.resources['record-1-genbank'].name, 'record-1-genbank-input.gb');
assert.equal(canonical.resources['record-1-genbank'].encoding, 'base64');
assert.equal(canonical.webFiles.resourceOriginalNames['record-1-genbank'], 'input.gb');
assert.equal(canonical.webFiles.circularInputOriginalName, 'input.gb');
assert.equal(canonical.renderRequest.output.prefix, 'web-session');
assert.deepEqual(canonical.renderRequest.output.formats, ['svg']);
assert.equal(canonical.renderRequest.output.overwrite, false);
assert.equal(canonical.renderRequest.output.interactiveMetadataPolicy, 'omit');
assert.deepEqual(canonical.renderRequest.diagramOptions.output, {
  legend: 'right',
  plotTitlePosition: 'none'
});
assert.equal(
  Object.prototype.hasOwnProperty.call(
    canonical.renderRequest.diagramOptions.output,
    'outputPrefix'
  ),
  false
);
const canonicalWithWebBindings = structuredClone(canonical);
canonicalWithWebBindings.resources['inactive-linear-source'] = {
  kind: 'web-file',
  name: 'inactive-linear-source.gb',
  type: 'application/octet-stream',
  size: genbank.size,
  lastModified: 0,
  encoding: 'base64',
  data: genbank.data
};
canonicalWithWebBindings.resources['inactive-comparison'] = {
  kind: 'web-file',
  name: 'inactive-comparison.tsv',
  type: 'application/octet-stream',
  size: 5,
  lastModified: 0,
  encoding: 'base64',
  data: btoa('match')
};
canonicalWithWebBindings.webFiles.bindings = {
  schema: 1,
  c_gb: {
    resourceId: 'record-1-genbank',
    name: 'bound-circular.gb',
    type: 'text/x-genbank',
    lastModified: 123
  },
  linearSeqs: [{
    uid: 'inactive-linear-uid',
    gb: {
      resourceId: 'inactive-linear-source',
      name: 'bound-linear.gb',
      type: 'application/genbank',
      lastModified: 456
    },
    gff: null,
    fasta: null,
    depth: null,
    blast: null,
    losat_gencode: 11
  }],
  linearComparisons: [{
    id: 'inactive-comparison-uid',
    queryUid: 'inactive-linear-uid',
    subjectUid: 'inactive-linear-uid-2',
    source: 'upload',
    file: {
      resourceId: 'inactive-comparison',
      name: 'bound-comparison.tsv',
      type: 'text/tab-separated-values',
      lastModified: 789
    }
  }],
  linearCanonicalComparisons: []
};
const webBindingProjection = projectCanonicalSessionRequest(
  canonicalWithWebBindings
);
assert.equal(webBindingProjection.files.c_gb.name, 'bound-circular.gb');
assert.equal(webBindingProjection.files.c_gb.type, 'text/x-genbank');
assert.equal(webBindingProjection.files.c_gb.lastModified, 123);
assert.equal(webBindingProjection.files.linearSeqs[0].uid, 'inactive-linear-uid');
assert.equal(webBindingProjection.files.linearSeqs[0].gb.name, 'bound-linear.gb');
assert.equal(
  webBindingProjection.files.linearComparisons[0].file.name,
  'bound-comparison.tsv'
);
const circularConfigOverrides = canonical.renderRequest.diagramOptions.configOverrides;
assert.ok(Object.keys(circularConfigOverrides).every((path) => path.includes('.')));
assert.ok(Object.values(circularConfigOverrides).every((value) => value !== null));
assert.equal(circularConfigOverrides['objects.scale.show'], true);
assert.equal(projectCanonicalSessionRequest(canonical).config.form.show_scale, true);
assert.equal(circularConfigOverrides['labels.circular.scope'], 'none');
assert.equal(circularConfigOverrides['labels.circular.placement'], 'horizontal');
assert.equal(canonical.renderRequest.diagramOptions.featureShapes.repeat_region, 'underlay');
assert.equal(canonical.renderRequest.diagramOptions.evalue, 1e-5);
assert.equal(canonical.renderRequest.diagramOptions.bitscore, 50);
assert.equal(canonical.renderRequest.diagramOptions.identity, 70);
assert.equal(canonical.renderRequest.diagramOptions.alignmentLength, 0);
assert.equal(
  circularConfigOverrides['objects.depth.large_tick_interval'],
  25
);
assert.equal(projectCanonicalSessionRequest(canonical).config.adv.depth_large_tick_interval, 25);
const sparseCircularScaleCanonical = structuredClone(canonical);
delete sparseCircularScaleCanonical.renderRequest.diagramOptions.configOverrides[
  'objects.scale.show'
];
assert.equal(
  projectCanonicalSessionRequest(sparseCircularScaleCanonical).config.form.show_scale,
  true
);
state.form.show_scale = false;
const hiddenCircularScaleCanonical = buildCanonicalRenderRequest({ state, filesData });
assert.equal(
  hiddenCircularScaleCanonical.renderRequest.diagramOptions.configOverrides[
    'objects.scale.show'
  ],
  false
);
assert.equal(
  projectCanonicalSessionRequest(hiddenCircularScaleCanonical).config.form.show_scale,
  false
);
state.form.show_scale = true;
state.form.labels_mode = 'out';
const outerLabelsCanonical = buildCanonicalRenderRequest({ state, filesData });
assert.equal(
  outerLabelsCanonical.renderRequest.diagramOptions.configOverrides[
    'labels.circular.scope'
  ],
  'outer'
);
assert.equal(
  projectCanonicalSessionRequest(outerLabelsCanonical).config.form.labels_mode,
  'out'
);
state.form.labels_mode = 'both';
const bothLabelsCanonical = buildCanonicalRenderRequest({ state, filesData });
assert.equal(
  bothLabelsCanonical.renderRequest.diagramOptions.configOverrides[
    'labels.circular.scope'
  ],
  'both'
);
assert.equal(
  projectCanonicalSessionRequest(bothLabelsCanonical).config.form.labels_mode,
  'both'
);
state.form.labels_mode = 'none';

state.circularRecordList.value = [{ selector: '#1', record_id: 'single.id' }];
state.form.prefix = '';
const implicitSingleCanonical = buildCanonicalRenderRequest({ state, filesData });
assert.equal(implicitSingleCanonical.renderRequest.grouping, 'single');
assert.equal(implicitSingleCanonical.webFiles.circularOutputPrefixExplicit, false);
assert.deepEqual(
  implicitSingleCanonical.renderRequest.records[0].selector,
  { kind: 'recordId', value: 'single.id' }
);
assert.equal(implicitSingleCanonical.renderRequest.output.prefix, 'single.id');
assert.equal(projectCanonicalSessionRequest(implicitSingleCanonical).config.form.prefix, '');

state.form.prefix = 'release.v1';
const explicitSingleCanonical = buildCanonicalRenderRequest({ state, filesData });
assert.equal(explicitSingleCanonical.renderRequest.output.prefix, 'release.v1');
assert.equal(explicitSingleCanonical.webFiles.circularOutputPrefixExplicit, true);
assert.equal(
  projectCanonicalSessionRequest(explicitSingleCanonical).config.form.prefix,
  'release.v1'
);

const explicitOneRecordBatch = structuredClone(explicitSingleCanonical);
explicitOneRecordBatch.renderRequest.grouping = 'batch';
explicitOneRecordBatch.renderRequest.output = [
  explicitOneRecordBatch.renderRequest.output
];
const oneRecordBatchProjection = projectCanonicalSessionRequest(
  explicitOneRecordBatch
);
assert.equal(
  oneRecordBatchProjection.config.adv.circular_grouping_intent,
  'batch'
);
Object.assign(state.form, oneRecordBatchProjection.config.form);
Object.assign(state.adv, oneRecordBatchProjection.config.adv);
const resavedOneRecordBatch = buildCanonicalRenderRequest({
  state,
  filesData: oneRecordBatchProjection.files
});
assert.equal(resavedOneRecordBatch.renderRequest.grouping, 'batch');
assert.equal(Array.isArray(resavedOneRecordBatch.renderRequest.output), true);
assert.equal(
  resavedOneRecordBatch.renderRequest.output.every(
    (entry) => entry.overwrite === false
  ),
  true
);
state.adv.circular_grouping_intent = 'auto';

const namingRecords = [
  { selector: '#1', record_id: 'dup' },
  { selector: '#2', record_id: 'dup_2' },
  { selector: '#3', record_id: 'dup' }
];
state.circularRecordList.value = namingRecords;
state.form.prefix = '';
const implicitBatchCanonical = buildCanonicalRenderRequest({ state, filesData });
assert.equal(implicitBatchCanonical.renderRequest.grouping, 'batch');
assert.deepEqual(
  implicitBatchCanonical.renderRequest.output.map((output) => output.prefix),
  ['dup', 'dup_2', 'dup_3']
);
assert.equal(
  projectCanonicalSessionRequest(implicitBatchCanonical).config.form.prefix,
  ''
);
assert.equal(
  projectCanonicalSessionRequest(implicitBatchCanonical).config.form.multi_record_canvas,
  false
);

state.form.prefix = 'release.v1';
const explicitBatchCanonical = buildCanonicalRenderRequest({ state, filesData });
assert.deepEqual(
  explicitBatchCanonical.renderRequest.output.map((output) => output.prefix),
  ['release.v1_1', 'release.v1_2', 'release.v1_3']
);
assert.equal(
  projectCanonicalSessionRequest(explicitBatchCanonical).config.form.prefix,
  'release.v1'
);

state.form.prefix = '';
state.form.multi_record_canvas = true;
const gridCanonical = buildCanonicalRenderRequest({ state, filesData });
assert.equal(gridCanonical.renderRequest.grouping, 'grid');
assert.equal(gridCanonical.renderRequest.output.prefix, 'dup');
assert.equal(projectCanonicalSessionRequest(gridCanonical).config.form.multi_record_canvas, true);

for (const schema of [3, 4]) {
  const branchOnlyCanonical = structuredClone(implicitBatchCanonical);
  branchOnlyCanonical.renderRequest.schema = schema;
  assert.throws(
    () => projectCanonicalSessionRequest(branchOnlyCanonical),
    /Unsupported canonical renderRequest schema/
  );
}
const missingCurrentGrouping = structuredClone(canonical);
delete missingCurrentGrouping.renderRequest.grouping;
assert.throws(
  () => projectCanonicalSessionRequest(missingCurrentGrouping),
  /grouping/
);

state.form.prefix = 'web-session';
state.form.multi_record_canvas = false;
state.circularRecordList.value = [];

state.adv.multi_record_size_mode = 'sqrt';
assert.throws(
  () => buildCanonicalRenderRequest({ state, filesData }),
  /Circular multi-record size mode/
);
state.adv.multi_record_size_mode = 'auto';
state.form.linear_track_layout = 'tuckin';
assert.throws(
  () => buildCanonicalRenderRequest({ state, filesData }),
  /Linear track layout/
);
state.form.linear_track_layout = 'middle';
state.adv.label_placement = 'on_feature';
assert.throws(
  () => buildCanonicalRenderRequest({ state, filesData }),
  /Linear label placement/
);
state.adv.label_placement = 'auto';
state.adv.circular_track_slots = [
  { id: 'stale', renderer: 'dinucleotide_skew', spacing: null, params: {} }
];
assert.doesNotThrow(
  () => buildCanonicalRenderRequest({ state, filesData }),
  'an inactive custom-stack draft must not block a simple-mode request'
);
state.adv.circular_track_slots_enabled = true;
state.adv.circular_track_slots_axis_index = 1;
state.adv.circular_track_slots = [
  {
    id: 'features',
    renderer: 'features',
    enabled: true,
    side: 'inside',
    params: { lane_direction: 'inside' }
  },
  {
    id: 'stale',
    renderer: 'dinucleotide_skew',
    enabled: false,
    spacing: null,
    params: {}
  }
];
assert.doesNotThrow(
  () => buildCanonicalRenderRequest({ state, filesData }),
  'a disabled invalid row must stay in the draft without entering the request'
);
state.adv.circular_track_slots[1].enabled = true;
assert.throws(
  () => buildCanonicalRenderRequest({ state, filesData }),
  /obsolete/
);
state.adv.circular_track_slots_enabled = false;
state.adv.circular_track_slots_axis_index = null;
state.adv.circular_track_slots = [];

state.adv.depth_tick_interval = 10;
assert.throws(
  () => buildCanonicalRenderRequest({ state, filesData }),
  /depth_tick_interval is obsolete/
);
delete state.adv.depth_tick_interval;
state.adv.depth_tracks = [{ tick_interval: 10 }];
assert.throws(
  () => buildCanonicalRenderRequest({ state, filesData }),
  /tick_interval is obsolete/
);
state.adv.depth_tracks = [];
state.losat.blastp.collinearMaxGeneGap = 3;
assert.throws(
  () => buildCanonicalRenderRequest({ state, filesData }),
  /collinearMaxGeneGap is obsolete/
);
delete state.losat.blastp.collinearMaxGeneGap;

for (const obsoleteKey of ['spacing', 'strict', 'compress', 'reserve']) {
  assert.throws(
    () => buildCircularTrackSlotSpec({
      id: `obsolete_${obsoleteKey}`,
      renderer: 'dinucleotide_skew',
      side: 'inside',
      [obsoleteKey]: obsoleteKey === 'spacing' ? '5' : true,
      params: { nt: 'GC' }
    }),
    /obsolete/
  );
  assert.throws(
    () => parseCircularTrackSlotSpec(
      `obsolete_${obsoleteKey}:dinucleotide_skew@${obsoleteKey}=true`
    ),
    /obsolete/
  );
}
const migratedLegacyCircularObject = migrateLegacyCircularTrackSlot({
  id: 'legacy_object',
  renderer: 'dinucleotide_skew',
  spacing: '5px',
  strict: true,
  params: { nt: 'GC', compress: true, reserve: true }
});
assert.equal(migratedLegacyCircularObject.spacing, undefined);
assert.equal(migratedLegacyCircularObject.strict, undefined);
assert.equal(migratedLegacyCircularObject.params.compress, undefined);
assert.equal(migratedLegacyCircularObject.params.reserve, undefined);
assert.equal(migratedLegacyCircularObject.inner_gap_px, '5px');
assert.equal(migratedLegacyCircularObject.outer_gap_px, '5px');
const migratedLegacyCircularSpec = migrateLegacyCircularTrackSlotSpec(
  'legacy_spec:dinucleotide_skew@spacing=4,strict=true,compress=true,reserve=true'
);
assert.equal(
  migratedLegacyCircularSpec,
  'legacy_spec:dinucleotide_skew@inner_gap_px=4,outer_gap_px=4'
);

state.adv.circular_track_slots_axis_index = 2;
const circularSlotsDisabled = buildCanonicalRenderRequest({ state, filesData });
assert.equal(circularSlotsDisabled.renderRequest.diagramOptions.tracks.circularTrackSlots, null);
assert.equal(circularSlotsDisabled.renderRequest.diagramOptions.tracks.circularTrackAxisIndex, null);
state.adv.circular_track_slots_axis_index = null;

const circularFeaturePresetCases = [
  { preset: 'middle', laneDirection: 'split', side: 'overlay', axisIndex: 0 },
  { preset: 'tuckin', laneDirection: 'inside', side: 'inside', axisIndex: 0 },
  { preset: 'spreadout', laneDirection: 'outside', side: 'outside', axisIndex: 1 }
];
for (const { preset, laneDirection, side, axisIndex } of circularFeaturePresetCases) {
  const featureSlot = {
    id: 'features',
    renderer: 'features',
    enabled: true,
    width: null,
    radius: null,
    inner_gap_px: null,
    outer_gap_px: null,
    side,
    z: 0,
    params: { lane_direction: laneDirection }
  };
  assert.deepEqual(
    resolveCircularTrackFeaturePlacement(featureSlot, preset),
    { laneDirection, side }
  );
  const liveSpec = buildCircularTrackSlotSpec(featureSlot, state.adv.nt, preset);
  assert.match(liveSpec, new RegExp(`(?:@|,)lane_direction=${laneDirection}(?:,|$)`));

  state.form.track_type = preset;
  state.adv.circular_track_slots_enabled = true;
  state.adv.circular_track_slots_axis_index = axisIndex;
  state.adv.circular_track_slots = [featureSlot];
  const presetCanonical = buildCanonicalRenderRequest({ state, filesData });
  const canonicalSpec = presetCanonical.renderRequest.diagramOptions.tracks.circularTrackSlots[0];
  assert.deepEqual(
    canonicalSpec,
    buildCircularTrackSlotPayload(featureSlot, state.adv.nt, preset)
  );

  const presetProjection = projectCanonicalSessionRequest(presetCanonical);
  const projectedFeature = presetProjection.config.adv.circular_track_slots[0];
  assert.deepEqual(
    resolveCircularTrackFeaturePlacement(projectedFeature, preset),
    { laneDirection, side }
  );
  assert.deepEqual(
    buildCircularTrackSlotPayload(projectedFeature, state.adv.nt, preset),
    canonicalSpec
  );
}

for (const explicitLaneDirection of ['inside', 'outside', 'split']) {
  const explicitSide = explicitLaneDirection === 'split' ? 'overlay' : explicitLaneDirection;
  const explicitSlot = {
    id: 'features',
    renderer: 'features',
    enabled: true,
    side: explicitSide,
    params: { lane_direction: explicitLaneDirection }
  };
  for (const preset of ['middle', 'tuckin', 'spreadout']) {
    const spec = buildCircularTrackSlotSpec(explicitSlot, 'GC', preset);
    assert.match(spec, new RegExp(`(?:@|,)lane_direction=${explicitLaneDirection}(?:,|$)`));
    const parsed = parseCircularTrackSlotSpec(spec, 0, 'GC', preset);
    assert.deepEqual(
      resolveCircularTrackFeaturePlacement(parsed, preset),
      { laneDirection: explicitLaneDirection, side: explicitSide }
    );
  }
}

state.form.track_type = 'tuckin';
state.adv.circular_track_slots_enabled = false;
state.adv.circular_track_slots_axis_index = null;
state.adv.circular_track_slots = [];

const legacyRepeatCanonical = structuredClone(canonical);
legacyRepeatCanonical.renderRequest.schema = 2;
legacyRepeatCanonical.renderRequest.diagramOptions.output.outputPrefix = 'ignored';
legacyRepeatCanonical.renderRequest.diagramOptions.selectedFeaturesSet = ['repeat_region'];
delete legacyRepeatCanonical.renderRequest.diagramOptions.featureShapes.repeat_region;
assert.equal(
  projectCanonicalSessionRequest(legacyRepeatCanonical).config.adv.feature_shapes.repeat_region,
  'rectangle'
);

const currentRepeatCanonical = structuredClone(canonical);
currentRepeatCanonical.renderRequest.diagramOptions.selectedFeaturesSet = ['repeat_region'];
delete currentRepeatCanonical.renderRequest.diagramOptions.featureShapes.repeat_region;
assert.equal(
  projectCanonicalSessionRequest(currentRepeatCanonical).config.adv.feature_shapes.repeat_region,
  'underlay'
);

const explicitRepeatCanonical = structuredClone(currentRepeatCanonical);
explicitRepeatCanonical.renderRequest.diagramOptions.featureShapes.repeat_region = 'underlay';
assert.equal(
  projectCanonicalSessionRequest(explicitRepeatCanonical).config.adv.feature_shapes.repeat_region,
  'underlay'
);

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
const semanticCanonical = buildCanonicalRenderRequest({
  state,
  filesData: {
    ...filesData,
    d_color: { name: 'original-default-colors.tsv' },
    t_color: { name: 'original-specific-colors.tsv' },
    whitelist: { name: 'original-whitelist.tsv' },
    qualifier_priority: { name: 'original-priority.tsv' }
  }
});
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
assert.deepEqual(
  {
    dColor: semanticProjection.files.d_color.name,
    tColor: semanticProjection.files.t_color.name,
    whitelist: semanticProjection.files.whitelist.name,
    priority: semanticProjection.files.qualifier_priority.name
  },
  {
    dColor: 'original-default-colors.tsv',
    tColor: 'original-specific-colors.tsv',
    whitelist: 'original-whitelist.tsv',
    priority: 'original-priority.tsv'
  }
);
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
const blastTableB = { ...genbank, name: 'hits-b.tsv', data: btoa('ref\tcmp_b\t98\t4\t0\t0\t1\t4\t4\t1\t1e-18\t45\n') };
const comparisonFasta = { ...genbank, name: 'comparison.fna', data: btoa('>cmp\nAACCGG\n') };
const comparisonFastaB = {
  ...comparisonFasta,
  name: 'comparison-b.fna',
  data: btoa('>cmp_b\nTTGGCC\n')
};
state.circularConservation.source = 'upload';
state.circularConservation.reference = 'subject';
state.circularConservation.labels = 'Comparison';
state.circularConservation.series = [{ label: 'Comparison', color: '#E15759' }];
state.circularConservation.ring_width = 14;
state.circularConservation.ring_gap = 3;
const conservationCanonical = buildCanonicalRenderRequest({
  state,
  filesData: {
    ...filesData,
    c_conservation_blasts: [blastTable],
    c_conservation_sequence_sources: [comparisonFasta],
  },
});
assert.deepEqual(
  conservationCanonical.renderRequest.diagramOptions.conservationFastaFiles,
  [{ resourceId: 'conservation-fasta-files-1', representation: 'file' }]
);
assert.equal(conservationCanonical.webFiles.conservationSequenceSources, undefined);
assert.equal(
  projectCanonicalSessionRequest(conservationCanonical).files.c_conservation_sequence_sources[0].data,
  comparisonFasta.data
);
const conservationProjection = projectCanonicalSessionRequest(conservationCanonical);
assert.equal(conservationProjection.config.circularConservation.enabled, true);
assert.equal(conservationProjection.config.circularConservation.source, 'upload');
assert.equal(conservationProjection.config.circularConservation.reference, 'subject');
assert.equal(conservationProjection.config.circularConservation.labels, 'Comparison');
assert.equal(conservationProjection.config.circularConservation.ring_width, 14);
assert.equal(conservationProjection.config.circularConservation.ring_gap, 3);
assert.deepEqual(
  conservationProjection.config.circularConservation.series.map((entry) => ({
    sourceIndex: entry.sourceIndex,
    label: entry.label,
    color: entry.color
  })),
  [{ sourceIndex: 0, label: 'Comparison', color: '#e15759' }]
);
assert.equal(
  conservationProjection.config.circularConservation.series[0].fileName,
  conservationProjection.files.c_conservation_blasts[0].name
);

state.circularConservation.labels = 'Stale A,Stale B';
state.circularConservation.series = [
  { sourceIndex: 1, label: 'Edited B', color: '#4E79A7' },
  { sourceIndex: 0, label: 'Edited A', color: '#E15759' }
];
const reorderedUploadCanonical = buildCanonicalRenderRequest({
  state,
  filesData: {
    ...filesData,
    c_conservation_blasts: [blastTable, blastTableB],
    c_conservation_sequence_sources: [comparisonFasta, comparisonFastaB]
  }
});
assert.deepEqual(
  reorderedUploadCanonical.renderRequest.diagramOptions.conservationLabels,
  ['Edited B', 'Edited A']
);
assert.deepEqual(
  reorderedUploadCanonical.renderRequest.diagramOptions.conservationColors,
  ['#4e79a7', '#e15759']
);
assert.deepEqual(
  reorderedUploadCanonical.renderRequest.diagramOptions.conservationBlastFiles
    .map((ref) => reorderedUploadCanonical.resources[ref.resourceId].data),
  [blastTableB.data, blastTable.data]
);
assert.deepEqual(
  reorderedUploadCanonical.renderRequest.diagramOptions.conservationFastaFiles
    .map((ref) => reorderedUploadCanonical.resources[ref.resourceId].data),
  [comparisonFastaB.data, comparisonFasta.data]
);
const reorderedUploadProjection = projectCanonicalSessionRequest(reorderedUploadCanonical);
assert.deepEqual(
  reorderedUploadProjection.config.circularConservation.series.map((entry) => entry.label),
  ['Edited B', 'Edited A']
);

state.circularConservation.source = 'losat';
state.circularConservation.series = [
  { sourceIndex: 1, label: 'Comparison B', color: '#4E79A7' },
  { sourceIndex: 0, label: 'Comparison A', color: '#E15759' }
];
const losatConservationCanonical = buildCanonicalRenderRequest({
  state,
  filesData: {
    ...filesData,
    c_conservation_blasts: [blastTable],
    c_conservation_fastas: [comparisonFasta, comparisonFastaB]
  }
});
assert.equal(
  losatConservationCanonical.renderRequest.diagramOptions.conservationBlastFiles,
  undefined
);
assert.deepEqual(
  losatConservationCanonical.webFiles.conservationLosatFastaSources,
  [
    'conservation-losat-fasta-files-1',
    'conservation-losat-fasta-files-2'
  ]
);
assert.equal(
  losatConservationCanonical.resources['conservation-losat-fasta-files-1'].data,
  comparisonFastaB.data
);
assert.deepEqual(
  projectCanonicalSessionRequest(losatConservationCanonical)
    .files.c_conservation_fastas.map((file) => file.data),
  [comparisonFastaB.data, comparisonFasta.data]
);
const losatReplayCanonical = buildCanonicalRenderRequest({
  state,
  filesData: {
    ...filesData,
    c_conservation_blasts: [blastTable, blastTableB],
    c_conservation_blasts_source: 'losat-cache',
    c_conservation_fastas: [comparisonFasta, comparisonFastaB]
  }
});
assert.equal(
  losatReplayCanonical.renderRequest.diagramOptions.conservationBlastFiles.length,
  2
);
assert.deepEqual(
  losatReplayCanonical.renderRequest.diagramOptions.conservationBlastFiles
    .map((ref) => losatReplayCanonical.resources[ref.resourceId].data),
  [blastTableB.data, blastTable.data]
);
assert.deepEqual(
  losatReplayCanonical.renderRequest.diagramOptions.conservationLabels,
  ['Comparison B', 'Comparison A']
);
assert.equal(losatReplayCanonical.webFiles.conservationBlastSource, 'losat-cache');
assert.deepEqual(
  losatReplayCanonical.webFiles.conservationLosatFastaSources
    .map((resourceId) => losatReplayCanonical.resources[resourceId].data),
  [comparisonFastaB.data, comparisonFasta.data]
);
const losatReplayProjection = projectCanonicalSessionRequest(losatReplayCanonical);
assert.equal(losatReplayProjection.files.c_conservation_blasts_source, 'losat-cache');
const mixedFastaCanonical = buildCanonicalRenderRequest({
  state,
  filesData: {
    ...filesData,
    c_conservation_blasts: [blastTable, blastTableB],
    c_conservation_blasts_source: 'losat-cache',
    c_conservation_fastas: [comparisonFasta, null]
  }
});
assert.deepEqual(
  mixedFastaCanonical.webFiles.conservationLosatFastaSources,
  [null, 'conservation-losat-fasta-files-2']
);
const mixedFastaProjection = projectCanonicalSessionRequest(mixedFastaCanonical);
assert.equal(mixedFastaProjection.files.c_conservation_fastas.length, 2);
assert.equal(mixedFastaProjection.files.c_conservation_fastas[0], null);
assert.equal(
  mixedFastaProjection.files.c_conservation_fastas[1].data,
  comparisonFasta.data
);
const mixedFastaRebuilt = buildCanonicalRenderRequest({
  state: {
    ...stateForCanonicalProjection(mixedFastaProjection),
    circularConservation: {
      ...structuredClone(mixedFastaProjection.config.circularConservation),
      source: 'losat'
    }
  },
  filesData: mixedFastaProjection.files
});
assert.deepEqual(
  mixedFastaRebuilt.webFiles.conservationLosatFastaSources,
  [null, 'conservation-losat-fasta-files-2']
);
const missingFastasCanonical = buildCanonicalRenderRequest({
  state,
  filesData: {
    ...filesData,
    c_conservation_blasts: [blastTable, blastTableB],
    c_conservation_blasts_source: 'losat-cache',
    c_conservation_fastas: []
  }
});
assert.deepEqual(
  missingFastasCanonical.webFiles.conservationLosatFastaSources,
  [null, null]
);
const missingFastasProjection = projectCanonicalSessionRequest(
  missingFastasCanonical
);
assert.deepEqual(
  missingFastasProjection.files.c_conservation_fastas,
  [null, null]
);
const missingFastasRebuilt = buildCanonicalRenderRequest({
  state: {
    ...stateForCanonicalProjection(missingFastasProjection),
    circularConservation: {
      ...structuredClone(missingFastasProjection.config.circularConservation),
      source: 'losat'
    }
  },
  filesData: missingFastasProjection.files
});
assert.deepEqual(
  missingFastasRebuilt.webFiles.conservationLosatFastaSources,
  [null, null]
);
state.circularConservation = {
  ...losatReplayProjection.config.circularConservation,
  source: 'losat'
};
const losatResavedCanonical = buildCanonicalRenderRequest({
  state,
  filesData: losatReplayProjection.files
});
assert.equal(
  losatResavedCanonical.renderRequest.diagramOptions.conservationBlastFiles.length,
  2
);
const oldDerivedSession = structuredClone(losatReplayCanonical);
delete oldDerivedSession.webFiles.conservationBlastSource;
assert.equal(
  projectCanonicalSessionRequest({
    ...oldDerivedSession,
    storedConfig: { circularConservation: { source: 'losat' } }
  }).files.c_conservation_blasts_source,
  'losat-cache'
);
state.circularConservation.source = 'upload';
state.circularConservation.reference = 'auto';
state.circularConservation.labels = '';
state.circularConservation.series = [];
state.circularConservation.ring_width = null;
state.circularConservation.ring_gap = null;

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
const annotationCanonical = buildCanonicalRenderRequest({ state, filesData });
assert.equal(annotationCanonical.renderRequest.diagramOptions.annotations.sets[0].id, 'review');
assert.equal(projectCanonicalSessionRequest(annotationCanonical).config.annotationSets[0].annotations[0].id, 'window');
assert.equal(
  projectCanonicalSessionRequest(annotationCanonical).config.annotationSets[0].annotations[0]
    .metadata._gbdraw_web_target_record_key,
  'linear-source::#2'
);
state.annotationSets = [];

state.adv.circular_label_placement = 'radial';
state.adv.outer_label_x_offset = 0.9;
state.adv.outer_label_y_offset = 0.91;
state.adv.inner_label_x_offset = 0.97;
state.adv.inner_label_y_offset = 0.98;
const radialCanonical = buildCanonicalRenderRequest({ state, filesData });
assert.equal(
  radialCanonical.renderRequest.diagramOptions.configOverrides['labels.circular.placement'],
  'radial'
);
assert.equal(
  radialCanonical.renderRequest.diagramOptions.configOverrides[
    'labels.unified_adjustment.outer_labels.x_radius_offset'
  ],
  0.9
);

const projection = projectCanonicalSessionRequest(canonical);
assert.equal(projection.mode, 'circular');
assert.equal(projection.inputType, 'gb');
assert.equal(projection.files.c_gb.data, genbank.data);
assert.equal(projection.files.c_gb.name, 'input.gb');
assert.equal(projection.config.form.prefix, 'web-session');
assert.equal(projection.config.form.labels_mode, 'none');
assert.equal(projection.config.adv.circular_label_placement, 'horizontal');
assert.equal(projection.config.adv.evalue, '0.00001');
assert.equal(projection.config.adv.min_bitscore, 50);
assert.equal(projection.config.adv.identity, 70);
assert.equal(projection.config.adv.alignment_length, 0);
const legacyCircularLabelScope = structuredClone(canonical);
legacyCircularLabelScope.renderRequest.schema = 2;
legacyCircularLabelScope.renderRequest.diagramOptions.output.outputPrefix = 'ignored';
legacyCircularLabelScope.renderRequest.diagramOptions.configOverrides = {
  show_labels: true,
  allow_inner_labels: true
};
assert.equal(
  projectCanonicalSessionRequest(legacyCircularLabelScope).config.form.labels_mode,
  'both'
);
const intermediateCircularLabelScope = structuredClone(canonical);
intermediateCircularLabelScope.renderRequest.diagramOptions.configOverrides = {
  'canvas.circular.show_labels': true,
  'canvas.circular.allow_inner_labels': false
};
assert.equal(
  projectCanonicalSessionRequest(intermediateCircularLabelScope).config.form.labels_mode,
  'out'
);
const legacyOutputCanonical = structuredClone(canonical);
legacyOutputCanonical.renderRequest.schema = 2;
legacyOutputCanonical.renderRequest.diagramOptions.output.outputPrefix = 'ignored-nested';
assert.equal(
  projectCanonicalSessionRequest(legacyOutputCanonical).config.form.prefix,
  'web-session'
);
const missingLegacyOutputPrefix = structuredClone(legacyOutputCanonical);
delete missingLegacyOutputPrefix.renderRequest.diagramOptions.output.outputPrefix;
assert.throws(
  () => projectCanonicalSessionRequest(missingLegacyOutputPrefix),
  /outputPrefix/
);
const currentWithLegacyOutputPrefix = structuredClone(canonical);
currentWithLegacyOutputPrefix.renderRequest.diagramOptions.output.outputPrefix = 'stale';
assert.throws(
  () => projectCanonicalSessionRequest(currentWithLegacyOutputPrefix),
  /outputPrefix/
);
const legacyCircularOptions = structuredClone(legacyOutputCanonical);
legacyCircularOptions.renderRequest.layout = { multiRecordSizeMode: 'sqrt' };
assert.equal(
  projectCanonicalSessionRequest(legacyCircularOptions).config.adv.multi_record_size_mode,
  'auto'
);
const currentCircularOptions = structuredClone(canonical);
currentCircularOptions.renderRequest.layout = { multiRecordSizeMode: 'sqrt' };
assert.throws(
  () => projectCanonicalSessionRequest(currentCircularOptions),
  /Circular multi-record size mode/
);
for (const schema of [1, 2]) {
  const sparseCircularCanonical = structuredClone(canonical);
  sparseCircularCanonical.renderRequest.schema = schema;
  sparseCircularCanonical.renderRequest.diagramOptions.output.outputPrefix = 'ignored';
  delete sparseCircularCanonical.renderRequest.diagramOptions.evalue;
  delete sparseCircularCanonical.renderRequest.diagramOptions.bitscore;
  delete sparseCircularCanonical.renderRequest.diagramOptions.identity;
  delete sparseCircularCanonical.renderRequest.diagramOptions.alignmentLength;
  delete sparseCircularCanonical.renderRequest.diagramOptions.selectedFeaturesSet;
  delete sparseCircularCanonical.renderRequest.diagramOptions.config;
  delete sparseCircularCanonical.renderRequest.diagramOptions.configOverrides;
  const sparseCircularProjection = projectCanonicalSessionRequest(sparseCircularCanonical);
  assert.equal(sparseCircularProjection.config.adv.evalue, '0.00001');
  assert.equal(sparseCircularProjection.config.adv.min_bitscore, 50);
  assert.equal(sparseCircularProjection.config.adv.identity, 70);
  assert.equal(sparseCircularProjection.config.adv.alignment_length, 0);
  assert.deepEqual(
    sparseCircularProjection.config.adv.features,
    ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'misc_RNA', 'repeat_region']
  );
  assert.equal(sparseCircularProjection.config.form.suppress_gc, false);
  assert.equal(sparseCircularProjection.config.form.suppress_skew, false);
}
const repeatedResourcePrefix = structuredClone(canonical);
repeatedResourcePrefix.resources['record-1-genbank'].name =
  'record-1-genbank-record-1-genbank-original.gb';
delete repeatedResourcePrefix.webFiles.resourceOriginalNames;
delete repeatedResourcePrefix.webFiles.circularInputOriginalName;
assert.equal(projectCanonicalSessionRequest(repeatedResourcePrefix).files.c_gb.name, 'original.gb');
const radialProjection = projectCanonicalSessionRequest(radialCanonical);
assert.equal(radialProjection.config.adv.circular_label_placement, 'radial');
assert.deepEqual(
  {
    outerX: radialProjection.config.adv.outer_label_x_offset,
    outerY: radialProjection.config.adv.outer_label_y_offset,
    innerX: radialProjection.config.adv.inner_label_x_offset,
    innerY: radialProjection.config.adv.inner_label_y_offset
  },
  { outerX: 0.9, outerY: 0.91, innerX: 0.97, innerY: 0.98 }
);

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
  inner_gap_px: '1',
  outer_gap_px: '1',
  side: 'inside',
  z: 0,
  params: {
    set_id: 'plastome_regions',
    show_labels: true,
    padding_px: 1,
    overflow: 'compress',
    layer: 'foreground'
  }
});

const legacyCircularSlotsCanonical = structuredClone(canonical);
legacyCircularSlotsCanonical.renderRequest.schema = 2;
legacyCircularSlotsCanonical.renderRequest.diagramOptions.output.outputPrefix = 'ignored';
legacyCircularSlotsCanonical.renderRequest.diagramOptions.tracks = {
  circularTrackSlots: [
    {
      id: 'legacy_object',
      renderer: 'dinucleotide_skew',
      side: 'inside',
      spacing: { value: 5, unit: 'px' },
      strict: true,
      params: { nt: 'GC', reserve: true }
    },
    'legacy_spec:dinucleotide_content@side=inside,spacing=4,compress=true,nt=AT'
  ],
  circularTrackAxisIndex: 0,
  centerReservedRadius: null
};
const legacyCircularSlotsProjection = projectCanonicalSessionRequest(
  legacyCircularSlotsCanonical
);
assert.deepEqual(
  legacyCircularSlotsProjection.config.adv.circular_track_slots.map((slot) => ({
    id: slot.id,
    inner: slot.inner_gap_px,
    outer: slot.outer_gap_px
  })),
  [
    { id: 'legacy_object', inner: '5', outer: '5' },
    { id: 'legacy_spec', inner: '4', outer: '4' }
  ]
);
const currentCircularObjectSlot = structuredClone(canonical);
currentCircularObjectSlot.renderRequest.diagramOptions.tracks.circularTrackSlots = [
  { id: 'old', renderer: 'dinucleotide_skew', spacing: 5, params: {} }
];
assert.throws(
  () => projectCanonicalSessionRequest(currentCircularObjectSlot),
  /obsolete/
);
const currentStructuredCircularSlot = structuredClone(canonical);
currentStructuredCircularSlot.renderRequest.diagramOptions.tracks.circularTrackSlots = [
  {
    id: 'current',
    renderer: 'dinucleotide_skew',
    innerGapPx: { value: 2, unit: 'px' },
    outerGapPx: { value: 3, unit: 'px' },
    params: {}
  }
];
assert.deepEqual(
  projectCanonicalSessionRequest(
    currentStructuredCircularSlot
  ).config.adv.circular_track_slots[0],
  {
    id: 'current',
    renderer: 'dinucleotide_skew',
    enabled: true,
    width: null,
    radius: null,
    inner_gap_px: '2',
    outer_gap_px: '3',
    side: null,
    z: 0,
    params: {}
  }
);
const currentCircularStringSlot = structuredClone(canonical);
currentCircularStringSlot.renderRequest.diagramOptions.tracks.circularTrackSlots = [
  'old:dinucleotide_skew@reserve=true'
];
assert.throws(
  () => projectCanonicalSessionRequest(currentCircularStringSlot),
  /obsolete/
);

const secondGenbank = {
  ...genbank,
  name: 'second.gb',
  data: btoa(genbankText.replaceAll('WEBTEST', 'WEBTWO'))
};
const multiCircularProjection = projectCanonicalSessionRequest({
  renderRequest: {
    ...canonical.renderRequest,
    grouping: 'grid',
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
state.adv.linear_track_slots_axis_index = 2;
const linearFilesData = {
  linearSeqs: [
    {
      uid: 'first', gb: genbank, losat_gencode: 11, losat_filename: 'first-to-second.losat.tsv',
      region_record_id: '', region_start: null, region_end: null, region_reverse: false
    },
    {
      uid: 'second', gb: genbank, losat_gencode: 4, losat_filename: 'second-to-third.losat.tsv',
      region_record_id: 'RecA', region_start: null, region_end: null, region_reverse: false
    },
    {
      uid: 'third', gb: genbank, losat_gencode: 1, losat_filename: '',
      region_record_id: '#2', region_start: 10, region_end: 20, region_reverse: true
    }
  ],
  linearComparisons: []
};
state.form.prefix = '';
const linearDefaultCanonical = buildCanonicalRenderRequest({ state, filesData: linearFilesData });
assert.equal(linearDefaultCanonical.renderRequest.grouping, 'single');
assert.equal(linearDefaultCanonical.renderRequest.output.prefix, 'out');
state.form.prefix = 'web-session';
const linearCanonical = buildCanonicalRenderRequest({ state, filesData: linearFilesData });
assert.ok(
  Object.keys(linearCanonical.renderRequest.diagramOptions.configOverrides)
    .every((path) => path.includes('.'))
);
assert.ok(
  Object.values(linearCanonical.renderRequest.diagramOptions.configOverrides)
    .every((value) => value !== null)
);
assert.equal(
  linearCanonical.renderRequest.diagramOptions.configOverrides['labels.linear.scope'],
  'none'
);
assert.equal(
  linearCanonical.renderRequest.diagramOptions.configOverrides['objects.scale.show'],
  true
);
assert.equal(projectCanonicalSessionRequest(linearCanonical).config.form.show_scale, true);
state.adv.block_stroke_width = 2;
state.adv.def_font_size = 16;
state.adv.label_font_size = 14;
state.adv.linear_definition_line_styles = {
  name: { font_size: 13, font_weight: 'bold', fill: '#112233' }
};
const styledLinearCanonical = buildCanonicalRenderRequest({
  state,
  filesData: linearFilesData
});
const styledLinearOverrides = styledLinearCanonical.renderRequest.diagramOptions.configOverrides;
assert.equal(styledLinearOverrides['objects.features.block_stroke_width.short'], 2);
assert.equal(styledLinearOverrides['objects.features.block_stroke_width.long'], 2);
assert.equal(styledLinearOverrides['objects.definition.linear.font_size.short'], 16);
assert.equal(styledLinearOverrides['objects.definition.linear.font_size.long'], 16);
assert.equal(styledLinearOverrides['labels.font_size.linear.short'], 14);
assert.equal(styledLinearOverrides['labels.font_size.linear.long'], 14);
assert.equal(
  styledLinearOverrides['objects.definition.linear.line_styles.name.font_size'],
  13
);
assert.equal(
  styledLinearOverrides['objects.definition.linear.line_styles.name.font_weight'],
  'bold'
);
assert.equal(
  styledLinearOverrides['objects.definition.linear.line_styles.name.fill'],
  '#112233'
);
assert.equal(
  Object.prototype.hasOwnProperty.call(
    styledLinearOverrides,
    'objects.definition.linear.line_styles'
  ),
  false
);
delete state.adv.block_stroke_width;
delete state.adv.def_font_size;
delete state.adv.label_font_size;
delete state.adv.linear_definition_line_styles;
const legacyLinearOptions = structuredClone(linearCanonical);
legacyLinearOptions.renderRequest.schema = 2;
legacyLinearOptions.renderRequest.diagramOptions.output.outputPrefix = 'ignored';
legacyLinearOptions.renderRequest.diagramOptions.configOverrides = {
  label_placement: 'on_feature',
  linear_track_layout: 'spreadout',
  show_labels: 'first'
};
const legacyLinearOptionsProjection = projectCanonicalSessionRequest(
  legacyLinearOptions
);
assert.equal(
  legacyLinearOptionsProjection.config.adv.label_placement,
  'above_feature'
);
assert.equal(
  legacyLinearOptionsProjection.config.form.linear_track_layout,
  'above'
);
assert.equal(
  legacyLinearOptionsProjection.config.form.show_labels_linear,
  'first'
);
const intermediateLinearLabelScope = structuredClone(linearCanonical);
intermediateLinearLabelScope.renderRequest.diagramOptions.configOverrides = {
  'canvas.linear.show_labels': 'orthogroup_top'
};
assert.equal(
  projectCanonicalSessionRequest(intermediateLinearLabelScope)
    .config.form.show_labels_linear,
  'orthogroup_top'
);
const currentLinearOptions = structuredClone(linearCanonical);
currentLinearOptions.renderRequest.diagramOptions.configOverrides.label_placement = 'on_feature';
assert.throws(
  () => projectCanonicalSessionRequest(currentLinearOptions),
  /Linear label placement/
);
currentLinearOptions.renderRequest.diagramOptions.configOverrides.label_placement = 'auto';
currentLinearOptions.renderRequest.diagramOptions.configOverrides.linear_track_layout = 'tuckin';
assert.throws(
  () => projectCanonicalSessionRequest(currentLinearOptions),
  /Linear track layout/
);
assert.equal(
  linearCanonical.renderRequest.diagramOptions.configOverrides[
    'objects.axis.linear.stroke_color'
  ],
  'lightgray'
);
state.form.linear_ruler_on_axis = true;
const rulerAxisCanonical = buildCanonicalRenderRequest({ state, filesData: linearFilesData });
assert.equal(
  rulerAxisCanonical.renderRequest.diagramOptions.configOverrides[
    'objects.axis.linear.stroke_color'
  ],
  'dimgray'
);
state.form.show_scale = false;
const hiddenRulerAxisCanonical = buildCanonicalRenderRequest({
  state,
  filesData: linearFilesData
});
assert.equal(
  hiddenRulerAxisCanonical.renderRequest.diagramOptions.configOverrides[
    'objects.scale.show'
  ],
  false
);
assert.equal(
  hiddenRulerAxisCanonical.renderRequest.diagramOptions.configOverrides[
    'canvas.linear.ruler_on_axis'
  ],
  true,
  'ordinary Web state must preserve a dormant ruler-on-axis choice'
);
assert.equal(
  hiddenRulerAxisCanonical.renderRequest.diagramOptions.configOverrides[
    'objects.axis.linear.stroke_color'
  ],
  'lightgray',
  'a hidden dormant ruler must not apply the ruler-axis managed color'
);
const hiddenRulerProjection = projectCanonicalSessionRequest(hiddenRulerAxisCanonical);
assert.equal(hiddenRulerProjection.config.form.show_scale, false);
assert.equal(hiddenRulerProjection.config.form.linear_ruler_on_axis, true);
state.form.show_scale = true;
state.form.linear_ruler_on_axis = false;
assert.deepEqual(linearCanonical.webFiles.linearRecordMetadata, [
  { recordKey: 'first', losatGencode: 11, losatFilename: 'first-to-second.losat.tsv' },
  { recordKey: 'second', losatGencode: 4, losatFilename: 'second-to-third.losat.tsv' },
  { recordKey: 'third', losatGencode: 1, losatFilename: '' }
]);
assert.equal(linearCanonical.renderRequest.diagramOptions.tracks.linearTrackSlots, null);
assert.equal(linearCanonical.renderRequest.diagramOptions.tracks.linearTrackAxisIndex, null);
state.adv.linear_track_slots_axis_index = null;
assert.equal(
  linearCanonical.renderRequest.diagramOptions.configOverrides[
    'canvas.linear.comparison_height'
  ],
  42.5
);
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
const arrangedCanonical = buildCanonicalRenderRequest({ state, filesData: linearFilesData });
assert.deepEqual(arrangedCanonical.renderRequest.layout, {
  recordGapPx: 30,
  multiRecordPositions: ['#1@1', '#2@1', '#3@2']
});

state.losatProgram.value = 'blastp';
linearFilesData.linearComparisons = [{
  id: 'selected-losat-pair', queryUid: 'first', subjectUid: 'third', source: 'losat', file: null
}];
const losatPairCanonical = buildCanonicalRenderRequest({ state, filesData: linearFilesData });
const generatedProtein = losatPairCanonical.renderRequest.comparisons.find(
  (comparison) => comparison.kind === 'generatedProteinComparison'
);
assert.deepEqual(generatedProtein.pairs, [{ queryRecordIndex: 0, subjectRecordIndex: 2 }]);
assert.equal(
  generatedProtein.settings.collinearityParams.parameters.maxUnitGap,
  2
);
assert.equal(
  projectCanonicalSessionRequest(losatPairCanonical).files.linearComparisons[0].source,
  'losat'
);
state.selectedOrthogroupAlignmentFeature.value = 'resolved-feature-anchor';
const resolvedProteinPlotTitlePosition = state.adv.plot_title_position;
state.adv.plot_title_position = 'bottom';
const resolvedProteinCanonical = buildCanonicalRenderRequest({
  state,
  filesData: linearFilesData,
  resolvedComparisons: [{
    kind: 'precomputedProteinComparison',
    queryRecordIndex: 0,
    subjectRecordIndex: 2,
    rows: [{
      query: 'protein-a',
      subject: 'protein-b',
      identity: 95,
      alignment_length: 20,
      mismatches: 1,
      gap_opens: 0,
      qstart: 10,
      qend: 30,
      sstart: 50,
      send: 70,
      evalue: 1e-20,
      bitscore: 120,
      group_kind: 'orthogroup',
      group_scope: 'cross_record'
    }]
  }]
});
const resolvedProtein = resolvedProteinCanonical.renderRequest.comparisons.find(
  (comparison) => comparison.kind === 'precomputedProteinComparison'
);
assert.equal(resolvedProtein.queryRecordIndex, 0);
assert.equal(resolvedProtein.subjectRecordIndex, 2);
const resolvedProteinSettings = resolvedProteinCanonical.renderRequest.comparisons.find(
  (comparison) => comparison.kind === 'generatedProteinComparison'
);
assert.ok(resolvedProteinSettings);
assert.equal(resolvedProteinSettings.mode, 'none');
assert.deepEqual(resolvedProteinSettings.pairs, []);
assert.equal(
  resolvedProteinSettings.settings.alignOrthogroupFeature,
  'resolved-feature-anchor'
);
const resolvedProteinTsv = Buffer.from(
  resolvedProteinCanonical.resources[resolvedProtein.resourceId].data,
  'base64'
).toString('utf8');
assert.match(resolvedProteinTsv, /group_kind/);
assert.match(resolvedProteinTsv, /cross_record/);
const resolvedProteinProjection = projectCanonicalSessionRequest(
  resolvedProteinCanonical
);
assert.deepEqual(
  resolvedProteinProjection.files.linearComparisons,
  [],
  'precomputed protein tables must not be projected as nucleotide uploads'
);
assert.deepEqual(
  resolvedProteinProjection.files.linearCanonicalComparisons.map(
    (comparison) => comparison.kind
  ),
  ['precomputedProteinComparison', 'generatedProteinComparison']
);
assert.equal(
  resolvedProteinProjection.config.losat.blastp.collinearMaxUnitGap,
  2
);
assert.equal(
  resolvedProteinProjection.pipelineState.selectedOrthogroupAlignmentFeature,
  'resolved-feature-anchor'
);
state.blastSource.value = resolvedProteinProjection.config.blastSource;
state.losatProgram.value = resolvedProteinProjection.config.losatProgram;
state.losat = structuredClone(resolvedProteinProjection.config.losat);
state.selectedOrthogroupAlignmentFeature.value =
  resolvedProteinProjection.pipelineState.selectedOrthogroupAlignmentFeature;
const resolvedProteinRoundTripCanonical = buildCanonicalRenderRequest({
  state,
  filesData: resolvedProteinProjection.files
});
const roundTripPrecomputed = resolvedProteinRoundTripCanonical.renderRequest.comparisons.find(
  (comparison) => comparison.kind === 'precomputedProteinComparison'
);
const roundTripGenerated = resolvedProteinRoundTripCanonical.renderRequest.comparisons.find(
  (comparison) => comparison.kind === 'generatedProteinComparison'
);
assert.deepEqual(
  {
    encoding: roundTripPrecomputed.encoding,
    queryRecordIndex: roundTripPrecomputed.queryRecordIndex,
    subjectRecordIndex: roundTripPrecomputed.subjectRecordIndex
  },
  {
    encoding: resolvedProtein.encoding,
    queryRecordIndex: resolvedProtein.queryRecordIndex,
    subjectRecordIndex: resolvedProtein.subjectRecordIndex
  }
);
assert.equal(
  Buffer.from(
    resolvedProteinRoundTripCanonical.resources[roundTripPrecomputed.resourceId].data,
    'base64'
  ).toString('utf8'),
  resolvedProteinTsv
);
assert.deepEqual(roundTripGenerated, resolvedProteinSettings);
const orthogroupResourceText = '{"schema":1,"valueKind":"orthogroupResult","value":{}}\n';
const collinearityResourceText = '{"schema":1,"valueKind":"blocks","value":[]}\n';
const resolvedWithMetadataCanonical = structuredClone(resolvedProteinCanonical);
resolvedWithMetadataCanonical.resources['orthogroup-result-test'] = {
  kind: 'orthogroup-result',
  name: 'orthogroups.json',
  type: 'application/json',
  size: new TextEncoder().encode(orthogroupResourceText).byteLength,
  lastModified: 0,
  encoding: 'base64',
  data: btoa(orthogroupResourceText)
};
resolvedWithMetadataCanonical.resources['collinearity-result-test'] = {
  kind: 'collinearity-result',
  name: 'collinearity.json',
  type: 'application/json',
  size: new TextEncoder().encode(collinearityResourceText).byteLength,
  lastModified: 0,
  encoding: 'base64',
  data: btoa(collinearityResourceText)
};
const generatedMetadataIndex =
  resolvedWithMetadataCanonical.renderRequest.comparisons.findIndex(
    (comparison) => comparison.kind === 'generatedProteinComparison'
  );
resolvedWithMetadataCanonical.renderRequest.comparisons.splice(
  generatedMetadataIndex,
  0,
  {
    kind: 'orthogroupResult',
    resourceId: 'orthogroup-result-test',
    encoding: 'canonicalJson'
  },
  {
    kind: 'collinearityResult',
    resourceId: 'collinearity-result-test',
    encoding: 'canonicalJson',
    valueKind: 'blocks'
  }
);
const resolvedWithMetadataProjection = projectCanonicalSessionRequest(
  resolvedWithMetadataCanonical
);
const resolvedWithFreshTable = buildCanonicalRenderRequest({
  state: stateForCanonicalProjection(resolvedWithMetadataProjection),
  filesData: resolvedWithMetadataProjection.files,
  resolvedComparisons: [{
    kind: 'precomputedProteinComparison',
    queryRecordIndex: 0,
    subjectRecordIndex: 2,
    rows: [{
      query: 'fresh-a',
      subject: 'fresh-b',
      identity: 88,
      alignment_length: 10,
      mismatches: 1,
      gap_opens: 0,
      qstart: 1,
      qend: 10,
      sstart: 1,
      send: 10,
      evalue: 1e-10,
      bitscore: 80
    }]
  }]
});
assert.deepEqual(
  resolvedWithFreshTable.renderRequest.comparisons.map(
    (comparison) => comparison.kind
  ),
  [
    'orthogroupResult',
    'collinearityResult',
    'precomputedProteinComparison',
    'generatedProteinComparison'
  ]
);
const rebuiltOrthogroups = resolvedWithFreshTable.renderRequest.comparisons.find(
  (comparison) => comparison.kind === 'orthogroupResult'
);
const rebuiltCollinearity = resolvedWithFreshTable.renderRequest.comparisons.find(
  (comparison) => comparison.kind === 'collinearityResult'
);
assert.equal(
  Buffer.from(
    resolvedWithFreshTable.resources[rebuiltOrthogroups.resourceId].data,
    'base64'
  ).toString('utf8'),
  orthogroupResourceText
);
assert.equal(rebuiltCollinearity.valueKind, 'blocks');
assert.equal(
  Buffer.from(
    resolvedWithFreshTable.resources[rebuiltCollinearity.resourceId].data,
    'base64'
  ).toString('utf8'),
  collinearityResourceText
);
const typedCanonicalFiles = structuredClone(resolvedWithMetadataProjection.files);
typedCanonicalFiles.linearCanonicalComparisons =
  typedCanonicalFiles.linearCanonicalComparisons.filter(
    (comparison) => [
      'orthogroupResult',
      'collinearityResult'
    ].includes(comparison.kind)
  );
const typedCanonicalState = stateForCanonicalProjection(
  resolvedWithMetadataProjection
);
typedCanonicalState.blastSource.value = 'files';
typedCanonicalState.losatProgram.value = 'blastn';
const typedCanonicalResult = buildCanonicalRenderRequest({
  state: typedCanonicalState,
  filesData: typedCanonicalFiles
});
assert.deepEqual(
  typedCanonicalResult.renderRequest.comparisons.map(
    (comparison) => comparison.kind
  ),
  [
    'orthogroupResult',
    'collinearityResult'
  ],
  'direct typed canonical results must not require Web protein-pipeline state'
);
const inactiveProteinState = stateForCanonicalProjection(
  resolvedWithMetadataProjection
);
inactiveProteinState.blastSource.value = 'files';
inactiveProteinState.losatProgram.value = 'blastn';
const withoutStaleProteinArtifacts = buildCanonicalRenderRequest({
  state: inactiveProteinState,
  filesData: resolvedWithMetadataProjection.files
});
assert.equal(
  withoutStaleProteinArtifacts.renderRequest.comparisons.some(
    (comparison) => [
      'precomputedProteinComparison',
      'orthogroupResult',
      'collinearityResult',
      'generatedProteinComparison'
    ].includes(comparison.kind)
  ),
  false,
  'saved protein artifacts must not leak into an active nucleotide pipeline'
);
state.selectedOrthogroupAlignmentFeature.value = '';
state.adv.plot_title_position = resolvedProteinPlotTitlePosition;

const linearProjection = projectCanonicalSessionRequest(linearCanonical);
assert.equal(linearProjection.config.adv.comparison_height, 42.5);
assert.equal(linearProjection.config.form.show_labels_linear, 'none');
let sparseLinearCanonical;
for (const schema of [1, 2]) {
  sparseLinearCanonical = structuredClone(linearCanonical);
  sparseLinearCanonical.renderRequest.schema = schema;
  sparseLinearCanonical.renderRequest.diagramOptions.output.outputPrefix = 'ignored';
  delete sparseLinearCanonical.renderRequest.diagramOptions.evalue;
  delete sparseLinearCanonical.renderRequest.diagramOptions.bitscore;
  delete sparseLinearCanonical.renderRequest.diagramOptions.identity;
  delete sparseLinearCanonical.renderRequest.diagramOptions.alignmentLength;
  delete sparseLinearCanonical.renderRequest.diagramOptions.selectedFeaturesSet;
  delete sparseLinearCanonical.renderRequest.diagramOptions.config;
  delete sparseLinearCanonical.renderRequest.diagramOptions.configOverrides;
  const sparseLinearProjection = projectCanonicalSessionRequest(sparseLinearCanonical);
  assert.equal(sparseLinearProjection.config.adv.evalue, '0.00001');
  assert.equal(sparseLinearProjection.config.adv.min_bitscore, 50);
  assert.equal(sparseLinearProjection.config.adv.identity, 70);
  assert.equal(sparseLinearProjection.config.adv.alignment_length, 0);
  assert.deepEqual(
    sparseLinearProjection.config.adv.features,
    ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'misc_RNA', 'repeat_region']
  );
  assert.equal(sparseLinearProjection.config.form.show_gc, true);
  assert.equal(sparseLinearProjection.config.form.show_skew, true);
  assert.equal(sparseLinearProjection.config.adv.axis_stroke_color, 'gray');
}
const nullOverrideSparseLinear = structuredClone(sparseLinearCanonical);
nullOverrideSparseLinear.renderRequest.diagramOptions.configOverrides = {
  show_gc: null,
  show_skew: null,
  linear_axis_stroke_color: null
};
const nullOverrideSparseLinearProjection = projectCanonicalSessionRequest(
  nullOverrideSparseLinear
);
assert.equal(nullOverrideSparseLinearProjection.config.form.show_gc, true);
assert.equal(nullOverrideSparseLinearProjection.config.form.show_skew, true);
assert.equal(nullOverrideSparseLinearProjection.config.adv.axis_stroke_color, 'gray');

const sparseCurrentLinear = structuredClone(linearCanonical);
delete sparseCurrentLinear.renderRequest.diagramOptions.evalue;
delete sparseCurrentLinear.renderRequest.diagramOptions.bitscore;
delete sparseCurrentLinear.renderRequest.diagramOptions.identity;
delete sparseCurrentLinear.renderRequest.diagramOptions.alignmentLength;
delete sparseCurrentLinear.renderRequest.diagramOptions.selectedFeaturesSet;
delete sparseCurrentLinear.renderRequest.diagramOptions.config;
delete sparseCurrentLinear.renderRequest.diagramOptions.configOverrides;
const sparseCurrentLinearProjection = projectCanonicalSessionRequest(sparseCurrentLinear);
assert.equal(sparseCurrentLinearProjection.config.adv.evalue, '1e-2');
assert.equal(sparseCurrentLinearProjection.config.adv.min_bitscore, 50);
assert.equal(sparseCurrentLinearProjection.config.adv.identity, 0);
assert.equal(sparseCurrentLinearProjection.config.adv.alignment_length, 0);
assert.equal(sparseCurrentLinearProjection.config.form.show_gc, false);
assert.equal(sparseCurrentLinearProjection.config.form.show_skew, false);
assert.equal(sparseCurrentLinearProjection.config.adv.axis_stroke_color, 'lightgray');

const explicitLinearThresholds = structuredClone(sparseLinearCanonical);
explicitLinearThresholds.renderRequest.diagramOptions.evalue = 0;
explicitLinearThresholds.renderRequest.diagramOptions.bitscore = 0;
explicitLinearThresholds.renderRequest.diagramOptions.identity = 0;
explicitLinearThresholds.renderRequest.diagramOptions.alignmentLength = 1;
explicitLinearThresholds.renderRequest.diagramOptions.selectedFeaturesSet = [];
explicitLinearThresholds.renderRequest.diagramOptions.configOverrides = {
  show_gc: false,
  show_skew: false,
  linear_axis_stroke_color: ''
};
const explicitLinearProjection = projectCanonicalSessionRequest(explicitLinearThresholds);
assert.equal(explicitLinearProjection.config.adv.evalue, '0');
assert.equal(explicitLinearProjection.config.adv.min_bitscore, 0);
assert.equal(explicitLinearProjection.config.adv.identity, 0);
assert.equal(explicitLinearProjection.config.adv.alignment_length, 1);
assert.deepEqual(explicitLinearProjection.config.adv.features, []);
assert.equal(explicitLinearProjection.config.form.show_gc, false);
assert.equal(explicitLinearProjection.config.form.show_skew, false);
assert.equal(explicitLinearProjection.config.adv.axis_stroke_color, '');
assert.deepEqual(
  linearProjection.files.linearSeqs.map((seq) => seq.gb.name),
  ['input.gb', 'input.gb', 'input.gb']
);
assert.deepEqual(
  linearProjection.files.linearSeqs.map((seq) => seq.losat_gencode),
  [11, 4, 1]
);
assert.deepEqual(
  linearProjection.files.linearSeqs.map((seq) => seq.losat_filename),
  ['first-to-second.losat.tsv', 'second-to-third.losat.tsv', '']
);
assert.deepEqual(
  linearProjection.files.linearSeqs.map((seq) => seq.region_record_id),
  ['', 'RecA', '#2']
);
assert.equal(linearProjection.files.linearSeqs[2].region_start, 10);
assert.equal(linearProjection.files.linearSeqs[2].region_end, 20);
assert.equal(linearProjection.files.linearSeqs[2].region_reverse, true);

const pythonConfigCanonical = structuredClone(linearCanonical);
delete pythonConfigCanonical.renderRequest.diagramOptions.configOverrides;
pythonConfigCanonical.renderRequest.diagramOptions.config = {
  canvas: {
    show_gc: true,
    show_skew: true,
    show_depth: false,
    show_labels: 'orthogroup_top',
    strandedness: true,
    resolve_overlaps: true,
    circular: { track_type: 'spreadout', allow_inner_labels: false },
    linear: {
      align_center: true,
      keep_definition_left_aligned: true,
      track_layout: 'above',
      track_axis_gap: 7,
      ruler_on_axis: true,
      comparison_height: 31,
      default_cds_height: { short: 23, long: 23 },
      default_gc_height: 18,
      depth_height: 12,
      normalize_length: true
    }
  },
  labels: {
    rendering: 'auto',
    font_size: { short: 13, long: 13, linear: { short: 15, long: 15 } },
    spacing: { circular: 4, linear: 5 },
    circular: { placement: 'radial' },
    linear: { placement: 'above_feature', rotation: 30 },
    filtering: { blacklist_keywords: ['hypothetical'] },
    unified_adjustment: {
      outer_labels: { x_radius_offset: 0.9, y_radius_offset: 0.91 },
      inner_labels: { x_radius_offset: 0.97, y_radius_offset: 0.98 }
    }
  },
  objects: {
    features: {
      block_stroke_color: '#111111', block_stroke_width: { short: 1, long: 1 },
      line_stroke_color: '#222222', line_stroke_width: { short: 2, long: 2 }
    },
    axis: {
      circular: { stroke_color: '#333333', stroke_width: { short: 3, long: 3 } },
      linear: { stroke_color: '#444444', stroke_width: { short: 4, long: 4 } }
    },
    definition: {
      linear: {
        font_size: { short: 16, long: 16 },
        show_replicon: true,
        show_accession: false,
        show_length: false,
        line_styles: { name: { font_weight: 'bold', fill: '#112233' } }
      },
      circular: { font_size: 18, plot_title_font_size: 30 }
    },
    legends: { color_rect_size: { short: 17, long: 17 }, font_size: { short: 19, long: 19 } },
    scale: {
      show: false,
      style: 'bar',
      font_size: { short: 21, long: 21 },
      ruler_label_font_size: { short: 22, long: 22 }
    },
    blast_match: { style: 'curve' }
  }
};
const pythonConfigProjection = projectCanonicalSessionRequest(pythonConfigCanonical);
assert.equal(pythonConfigProjection.config.form.separate_strands, true);
assert.equal(pythonConfigProjection.config.form.align_center, true);
assert.equal(pythonConfigProjection.config.form.keep_definition_left_aligned, true);
assert.equal(pythonConfigProjection.config.form.linear_track_layout, 'above');
assert.equal(pythonConfigProjection.config.form.linear_ruler_on_axis, true);
assert.equal(pythonConfigProjection.config.form.show_scale, false);
assert.equal(pythonConfigProjection.config.form.normalize_length, true);
assert.equal(pythonConfigProjection.config.form.show_labels_linear, 'orthogroup_top');
assert.equal(pythonConfigProjection.config.adv.comparison_height, 31);
assert.equal(pythonConfigProjection.config.adv.track_axis_gap, 7);
assert.equal(pythonConfigProjection.config.adv.linear_show_replicon, true);
assert.equal(pythonConfigProjection.config.adv.linear_show_accession, false);
assert.equal(pythonConfigProjection.config.adv.linear_show_length, false);
assert.equal(pythonConfigProjection.config.adv.block_stroke_width, 1);
assert.equal(pythonConfigProjection.config.adv.block_stroke_color, '#111111');
assert.equal(pythonConfigProjection.config.adv.line_stroke_width, 2);
assert.equal(pythonConfigProjection.config.adv.line_stroke_color, '#222222');
assert.equal(pythonConfigProjection.config.adv.axis_stroke_width, 4);
assert.equal(pythonConfigProjection.config.adv.axis_stroke_color, '#444444');
assert.equal(pythonConfigProjection.config.adv.def_font_size, 16);
assert.equal(pythonConfigProjection.config.adv.feature_height, 23);
assert.equal(pythonConfigProjection.config.adv.label_font_size, 15);
assert.equal(pythonConfigProjection.config.adv.label_placement, 'above_feature');
assert.equal(pythonConfigProjection.config.adv.legend_box_size, 17);
assert.equal(pythonConfigProjection.config.adv.legend_font_size, 19);
assert.equal(pythonConfigProjection.config.adv.scale_font_size, 21);
assert.deepEqual(
  pythonConfigProjection.config.adv.linear_definition_line_styles,
  { name: { font_weight: 'bold', fill: '#112233' } }
);
assert.equal(pythonConfigProjection.config.adv.pairwise_match_style, 'curve');
assert.equal(pythonConfigProjection.config.blacklistText, 'hypothetical');

const rulerFontCanonical = structuredClone(pythonConfigCanonical);
rulerFontCanonical.renderRequest.diagramOptions.config.objects.scale.style = 'ruler';
assert.equal(projectCanonicalSessionRequest(rulerFontCanonical).config.adv.scale_font_size, 22);

const sparsePythonScaleCanonical = structuredClone(pythonConfigCanonical);
delete sparsePythonScaleCanonical.renderRequest.diagramOptions.config.objects.scale.show;
assert.equal(
  projectCanonicalSessionRequest(sparsePythonScaleCanonical).config.form.show_scale,
  true
);

const pythonCircularConfigCanonical = structuredClone(canonical);
delete pythonCircularConfigCanonical.renderRequest.diagramOptions.configOverrides;
pythonCircularConfigCanonical.renderRequest.diagramOptions.config = structuredClone(
  pythonConfigCanonical.renderRequest.diagramOptions.config
);
pythonCircularConfigCanonical.renderRequest.diagramOptions.config.canvas.show_labels = true;
const pythonCircularConfigProjection = projectCanonicalSessionRequest(pythonCircularConfigCanonical);
assert.equal(pythonCircularConfigProjection.config.adv.axis_stroke_width, 3);
assert.equal(pythonCircularConfigProjection.config.adv.axis_stroke_color, '#333333');
assert.equal(pythonCircularConfigProjection.config.adv.def_font_size, 18);
assert.equal(pythonCircularConfigProjection.config.adv.plot_title_font_size, 30);
assert.equal(pythonCircularConfigProjection.config.adv.label_font_size, 13);
assert.equal(pythonCircularConfigProjection.config.adv.outer_label_x_offset, 0.9);
assert.equal(pythonCircularConfigProjection.config.adv.outer_label_y_offset, 0.91);
assert.equal(pythonCircularConfigProjection.config.adv.inner_label_x_offset, 0.97);
assert.equal(pythonCircularConfigProjection.config.adv.inner_label_y_offset, 0.98);
assert.equal(pythonCircularConfigProjection.config.form.show_scale, false);

const emptyBlacklistCanonical = structuredClone(pythonConfigCanonical);
emptyBlacklistCanonical.renderRequest.diagramOptions.config.labels.filtering.blacklist_keywords = [];
const emptyBlacklistProjection = projectCanonicalSessionRequest(emptyBlacklistCanonical);
assert.equal(emptyBlacklistProjection.config.filterMode, 'None');
assert.equal(emptyBlacklistProjection.config.blacklistText, '');

const explicitOverrideCanonical = structuredClone(pythonConfigCanonical);
explicitOverrideCanonical.renderRequest.diagramOptions.configOverrides = {
  'canvas.strandedness': false,
  'canvas.linear.align_center': false,
  'canvas.linear.keep_definition_left_aligned': false
};
const explicitOverrideProjection = projectCanonicalSessionRequest(explicitOverrideCanonical);
assert.equal(explicitOverrideProjection.config.form.separate_strands, false);
assert.equal(explicitOverrideProjection.config.form.align_center, false);
assert.equal(explicitOverrideProjection.config.form.keep_definition_left_aligned, false);

state.adv.comparison_height = null;
const autoHeightCanonical = buildCanonicalRenderRequest({ state, filesData: linearFilesData });
assert.equal(
  Object.prototype.hasOwnProperty.call(
    autoHeightCanonical.renderRequest.diagramOptions.configOverrides,
    'canvas.linear.comparison_height'
  ),
  false
);
assert.equal(projectCanonicalSessionRequest(autoHeightCanonical).config.adv.comparison_height, null);
state.adv.comparison_height = -2;
assert.throws(
  () => buildCanonicalRenderRequest({ state, filesData: linearFilesData }),
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

const depthAText = 'position\tdepth\n1\t10\n';
const depthA = {
  ...genbank,
  name: 'sample-a.depth.tsv',
  type: 'text/tab-separated-values',
  size: new TextEncoder().encode(depthAText).byteLength,
  data: btoa(depthAText)
};
const depthBText = 'position\tdepth\n1\t20\n';
const depthB = {
  ...genbank,
  name: 'sample-b.depth.tsv',
  type: 'text/tab-separated-values',
  size: new TextEncoder().encode(depthBText).byteLength,
  data: btoa(depthBText)
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
const circularSparseCanonical = buildCanonicalRenderRequest({
  state,
  filesData: circularSparseFilesData
});
assert.equal(circularSparseCanonical.renderRequest.records.length, 2);
assert.deepEqual(
  circularSparseCanonical.renderRequest.diagramOptions.depthTracks.map((track) => (
    track.source.map((entry) => Boolean(entry?.resourceId))
  )),
  [[true, false], [false, true]]
);
assert.equal(
  circularSparseCanonical.renderRequest.diagramOptions.depthTracks.every(
    (track) => track.height === null
  ),
  true
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
const circularSparseRebuilt = buildCanonicalRenderRequest({
  state,
  filesData: circularSparseProjection.files
});
assert.deepEqual(
  circularSparseRebuilt.renderRequest.diagramOptions.depthTracks.map((track) => (
    track.source.map((entry) => Boolean(entry?.resourceId))
  )),
  [[true, false], [false, true]]
);

const circularLegacyFlatCanonical = buildCanonicalRenderRequest({
  state,
  filesData: {
    c_gb: combinedGenbank,
    c_depth: [depthA, depthB],
    linearSeqs: []
  }
});
assert.deepEqual(
  circularLegacyFlatCanonical.renderRequest.diagramOptions.depthTracks.map(
    (track) => Boolean(track.source?.resourceId)
  ),
  [true, true]
);
assert.deepEqual(
  projectCanonicalSessionRequest(circularLegacyFlatCanonical).files.c_depth.map(
    (row) => row.map((entry) => Boolean(entry))
  ),
  [[true, true], [true, true]]
);
assert.throws(
  () => buildCanonicalRenderRequest({
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
state.adv.plot_title_position = 'bottom';
state.form.multi_record_canvas = false;
state.circularRecordList.value = [];
const nullOnlyDepthCanonical = buildCanonicalRenderRequest({
  state,
  filesData: {
    ...linearFilesData,
    linearSeqs: linearFilesData.linearSeqs.map((seq) => ({ ...seq, depth: [null] }))
  }
});
assert.equal(nullOnlyDepthCanonical.renderRequest.diagramOptions.depthTrackFiles, undefined);
assert.equal(nullOnlyDepthCanonical.renderRequest.diagramOptions.depthTracks, undefined);
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
const sparseDepthCanonical = buildCanonicalRenderRequest({ state, filesData: sparseDepthFilesData });
assert.equal(sparseDepthCanonical.renderRequest.diagramOptions.depthTracks.length, 2);
assert.equal(sparseDepthCanonical.renderRequest.diagramOptions.depthTracks[0].source.length, 2);
assert.equal(sparseDepthCanonical.renderRequest.diagramOptions.depthTracks[0].source[1], null);
assert.equal(sparseDepthCanonical.renderRequest.diagramOptions.depthTracks[1].source[0], null);
assert.deepEqual(
  sparseDepthCanonical.renderRequest.diagramOptions.depthTracks.map((track) => track.label),
  ['Sample A', 'Sample B']
);
[
  'depthTrackFiles',
  'depthTrackLabels',
  'depthTrackColors',
  'depthTrackHeights',
  'depthTrackLargeTickIntervals',
  'depthTrackSmallTickIntervals',
  'depthTrackTickFontSizes'
].forEach((fieldName) => {
  assert.equal(
    Object.hasOwn(sparseDepthCanonical.renderRequest.diagramOptions, fieldName),
    false
  );
});
const mixedCanonicalDepth = structuredClone(sparseDepthCanonical);
mixedCanonicalDepth.renderRequest.diagramOptions.depthTrackFiles = [
  [sparseDepthCanonical.renderRequest.diagramOptions.depthTracks[0].source[0], null],
  [null, sparseDepthCanonical.renderRequest.diagramOptions.depthTracks[1].source[1]]
];
assert.throws(
  () => projectCanonicalSessionRequest(mixedCanonicalDepth),
  /depthTracks cannot be combined with legacy depth fields: depthTrackFiles/
);
const wrongCanonicalDepthCardinality = structuredClone(sparseDepthCanonical);
wrongCanonicalDepthCardinality.renderRequest.diagramOptions.depthTracks[0].source =
  wrongCanonicalDepthCardinality.renderRequest.diagramOptions.depthTracks[0].source.slice(0, 1);
assert.throws(
  () => projectCanonicalSessionRequest(wrongCanonicalDepthCardinality),
  /one source per displayed record \(2\)/
);
const emptyCanonicalDepth = structuredClone(sparseDepthCanonical);
emptyCanonicalDepth.renderRequest.diagramOptions.depthTracks[0].source = [null, null];
assert.throws(
  () => projectCanonicalSessionRequest(emptyCanonicalDepth),
  /logical track index 0.*no source/
);
const circularDepthHeight = structuredClone(circularSparseCanonical);
circularDepthHeight.renderRequest.diagramOptions.depthTracks[0].height = 12;
assert.throws(
  () => projectCanonicalSessionRequest(circularDepthHeight),
  /height must be null for Circular requests/
);
const missingCanonicalDepthField = structuredClone(sparseDepthCanonical);
delete missingCanonicalDepthField.renderRequest.diagramOptions.depthTracks[0].label;
assert.throws(
  () => projectCanonicalSessionRequest(missingCanonicalDepthField),
  /missing required field\(s\): label/
);
const blankCanonicalDepthLabel = structuredClone(sparseDepthCanonical);
blankCanonicalDepthLabel.renderRequest.diagramOptions.depthTracks[0].label = ' ';
assert.throws(
  () => projectCanonicalSessionRequest(blankCanonicalDepthLabel),
  /label must be null or a non-empty string/
);
const nonStringCanonicalDepthColor = structuredClone(sparseDepthCanonical);
nonStringCanonicalDepthColor.renderRequest.diagramOptions.depthTracks[0].color = 123;
assert.throws(
  () => projectCanonicalSessionRequest(nonStringCanonicalDepthColor),
  /color must be null or a non-empty string/
);
const unknownCanonicalDepthField = structuredClone(sparseDepthCanonical);
unknownCanonicalDepthField.renderRequest.diagramOptions.depthTracks[0].legacyLabel = 'old';
assert.throws(
  () => projectCanonicalSessionRequest(unknownCanonicalDepthField),
  /contains unknown field\(s\): legacyLabel/
);
assert.equal(sparseDepthCanonical.renderRequest.diagramOptions.tracks.linearTrackAxisIndex, 1);
assert.deepEqual(
  sparseDepthCanonical.renderRequest.diagramOptions.tracks.linearTrackSlots.map((slot) => slot.id),
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
const filteredAxisCanonical = buildCanonicalRenderRequest({ state, filesData: sparseDepthFilesData });
assert.equal(filteredAxisCanonical.renderRequest.diagramOptions.tracks.linearTrackAxisIndex, 0);
assert.deepEqual(
  filteredAxisCanonical.renderRequest.diagramOptions.tracks.linearTrackSlots.map((slot) => slot.id),
  ['features', 'depth_b']
);

const circularRequestState = ({
  form: formOverrides = {},
  adv: advOverrides = {},
  annotationSets: nextAnnotationSets = []
} = {}) => ({
  ...state,
  mode: ref('circular'),
  cInputType: ref('gb'),
  form: {
    ...state.form,
    show_depth: false,
    suppress_gc: false,
    suppress_skew: false,
    track_type: 'tuckin',
    ...formOverrides
  },
  adv: {
    ...state.adv,
    features: ['CDS'],
    feature_shapes: { CDS: 'arrow' },
    circular_track_slots_enabled: false,
    circular_track_slots_axis_index: null,
    circular_track_slots: [],
    linear_track_slots_enabled: false,
    depth_tracks: [],
    feature_width_circular: null,
    depth_width_circular: null,
    gc_content_width_circular: null,
    gc_content_radius_circular: null,
    gc_skew_width_circular: null,
    gc_skew_radius_circular: null,
    ...advOverrides
  },
  annotationSets: structuredClone(nextAnnotationSets)
});

const geometryCases = [
  ['feature_width_circular', 'features', 'width', { value: 17, unit: 'px' }],
  ['depth_width_circular', 'depth', 'width', { value: 18, unit: 'px' }],
  ['gc_content_width_circular', 'dinucleotide_content', 'width', { value: 19, unit: 'px' }],
  ['gc_content_radius_circular', 'dinucleotide_content', 'radius', { value: 0.61, unit: 'factor' }],
  ['gc_skew_width_circular', 'dinucleotide_skew', 'width', { value: 20, unit: 'px' }],
  ['gc_skew_radius_circular', 'dinucleotide_skew', 'radius', { value: 0.52, unit: 'factor' }]
];
const sparseScaleState = circularRequestState();
delete sparseScaleState.form.show_scale;
assert.equal(
  buildCanonicalRenderRequest({ state: sparseScaleState, filesData })
    .renderRequest.diagramOptions.configOverrides['objects.scale.show'],
  true,
  'sparse Web form state must retain the visible-scale default'
);
for (const [field, renderer, property, expected] of geometryCases) {
  const needsDepth = field === 'depth_width_circular';
  const geometryState = circularRequestState({
    form: { show_depth: needsDepth },
    adv: { [field]: expected.value }
  });
  const geometryCanonical = buildCanonicalRenderRequest({
    state: geometryState,
    filesData: {
      ...filesData,
      c_depth: needsDepth ? [depthA] : null
    }
  });
  const slots = geometryCanonical.renderRequest.diagramOptions.tracks.circularTrackSlots;
  assert.ok(Array.isArray(slots), `${field} should materialize the implicit stack`);
  const slot = slots.find((entry) => entry.renderer === renderer);
  assert.deepEqual(slot[property], expected, `${field} should patch ${renderer}.${property}`);
}

const hiddenScaleGeometryCanonical = buildCanonicalRenderRequest({
  state: circularRequestState({
    form: { show_scale: false },
    adv: { feature_width_circular: 17 }
  }),
  filesData
});
assert.equal(
  hiddenScaleGeometryCanonical.renderRequest.diagramOptions.configOverrides[
    'objects.scale.show'
  ],
  false
);
assert.equal(
  hiddenScaleGeometryCanonical.renderRequest.diagramOptions.tracks.circularTrackSlots
    .some((slot) => slot.renderer === 'ticks'),
  false,
  'simple geometry materialization must not restore hidden coordinate ticks'
);

const explicitTicksWithHiddenSimpleScale = buildCanonicalRenderRequest({
  state: circularRequestState({
    form: { show_scale: false },
    adv: {
      circular_track_slots_enabled: true,
      circular_track_slots_axis_index: 0,
      circular_track_slots: [
        {
          id: 'features',
          renderer: 'features',
          enabled: true,
          side: 'inside',
          params: { lane_direction: 'inside' }
        },
        {
          id: 'ticks',
          renderer: 'ticks',
          enabled: true,
          side: 'inside',
          params: { tick_label_layout: 'label_in_tick_out' }
        }
      ]
    }
  }),
  filesData
});
assert.equal(
  explicitTicksWithHiddenSimpleScale.renderRequest.diagramOptions.tracks.circularTrackSlots
    .some((slot) => slot.renderer === 'ticks'),
  true,
  'enabled custom Ticks slots must remain authoritative over the dormant simple setting'
);

const nullGeometryCanonical = buildCanonicalRenderRequest({
  state: circularRequestState(),
  filesData
});
assert.equal(
  nullGeometryCanonical.renderRequest.diagramOptions.tracks.circularTrackSlots,
  null,
  'all-Auto Circular geometry must preserve the implicit default request'
);

for (const invalidGeometry of [0, -1, Number.NaN, Number.POSITIVE_INFINITY]) {
  assert.throws(
    () => buildCanonicalRenderRequest({
      state: circularRequestState({
        adv: { feature_width_circular: invalidGeometry }
      }),
      filesData
    }),
    /Feature Width must be Auto or a positive finite number/
  );
}

const retainedDepthFile = { ...genbank, name: 'retained-depth.tsv', data: btoa('position\tdepth\n1\t5\n') };
const inactiveDepthCanonical = buildCanonicalRenderRequest({
  state: circularRequestState({ form: { show_depth: false } }),
  filesData: { ...filesData, c_depth: [retainedDepthFile] }
});
assert.equal(inactiveDepthCanonical.renderRequest.diagramOptions.depthTracks, undefined);
assert.equal(
  inactiveDepthCanonical.renderRequest.diagramOptions.configOverrides['canvas.show_depth'],
  false
);
assert.ok(
  !Object.values(inactiveDepthCanonical.resources)
    .some((resource) => resource.kind === 'depth-track-file')
);

const customDepthState = circularRequestState({
  form: { show_depth: false },
  adv: {
    circular_track_slots_enabled: true,
    circular_track_slots_axis_index: 0,
    circular_track_slots: [{
      id: 'depth',
      renderer: 'depth',
      enabled: true,
      side: 'inside',
      width: null,
      radius: null,
      inner_gap_px: null,
      outer_gap_px: null,
      z: 0,
      params: { track_index: 0 }
    }]
  }
});
const customDepthCanonical = buildCanonicalRenderRequest({
  state: customDepthState,
  filesData: { ...filesData, c_depth: [retainedDepthFile] }
});
assert.equal(
  customDepthCanonical.renderRequest.diagramOptions.configOverrides['canvas.show_depth'],
  true
);
assert.equal(customDepthCanonical.renderRequest.diagramOptions.depthTracks.length, 1);

const customNoDepthCanonical = buildCanonicalRenderRequest({
  state: circularRequestState({
    form: { show_depth: true },
    adv: {
      circular_track_slots_enabled: true,
      circular_track_slots_axis_index: 0,
      circular_track_slots: [{
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
      }]
    }
  }),
  filesData: { ...filesData, c_depth: [retainedDepthFile] }
});
assert.equal(customNoDepthCanonical.renderRequest.diagramOptions.depthTracks, undefined);
assert.equal(
  customNoDepthCanonical.renderRequest.diagramOptions.configOverrides['canvas.show_depth'],
  false
);

const circularAxisRebaseCanonical = buildCanonicalRenderRequest({
  state: circularRequestState({
    adv: {
      circular_track_slots_enabled: true,
      circular_track_slots_axis_index: 2,
      circular_track_slots: [
        {
          id: 'outside',
          renderer: 'spacer',
          enabled: true,
          side: 'outside',
          width: '6px',
          radius: null,
          inner_gap_px: 0,
          outer_gap_px: 0,
          z: 0,
          params: {}
        },
        {
          id: 'disabled',
          renderer: 'spacer',
          enabled: false,
          side: 'outside',
          width: '6px',
          radius: null,
          inner_gap_px: 0,
          outer_gap_px: 0,
          z: 0,
          params: {}
        },
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
      ]
    }
  }),
  filesData
});
assert.equal(
  circularAxisRebaseCanonical.renderRequest.diagramOptions.tracks.circularTrackAxisIndex,
  1
);
assert.deepEqual(
  circularAxisRebaseCanonical.renderRequest.diagramOptions.tracks.circularTrackSlots
    .map((slot) => slot.id),
  ['outside', 'features']
);

const nestedStyleOverride = {
  stroke: '#112233',
  strokeWidth: 2,
  hatch: { pattern: 'diagonal', spacing: 4 },
  label: { color: '#445566', fontSize: 9 }
};
const p3AnnotationSet = [{
  id: 'review',
  annotations: [{
    id: 'line',
    target: {
      kind: 'coordinateSpan',
      record: null,
      start: 1,
      end: 3,
      coordinateSpace: 'source',
      wrapsOrigin: false,
      outOfBounds: 'clip'
    },
    label: 'Line',
    mark: 'line',
    lane: null,
    style: null,
    legendLabel: null,
    metadata: {}
  }],
  defaultStyle: {
    stroke: '#404040',
    strokeWidth: 1.5,
    strokeDasharray: [],
    lineCap: 'tick',
    fill: null,
    fillOpacity: 0.2,
    hatch: null,
    labelColor: '#202020',
    labelFontSize: null,
    labelOrientation: 'auto',
    labelPosition: 'center',
    labelOffset: 4
  },
  legendLabel: null
}];

const circularP3State = circularRequestState({
  annotationSets: p3AnnotationSet,
  adv: {
    circular_track_slots_enabled: true,
    circular_track_slots_axis_index: 0,
    circular_track_slots: [
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
      },
      {
        id: 'space',
        renderer: 'spacer',
        enabled: true,
        side: 'inside',
        width: '8px',
        radius: null,
        inner_gap_px: 1,
        outer_gap_px: 2,
        z: 0,
        params: {}
      },
      {
        id: 'review_overlay',
        renderer: 'annotations',
        enabled: true,
        side: 'overlay',
        width: null,
        radius: null,
        inner_gap_px: null,
        outer_gap_px: null,
        z: 2,
        params: {
          set_id: 'review',
          marks: ['line', 'band'],
          lane_gap_px: 5,
          padding_px: 0,
          overflow: 'error',
          show_labels: true,
          anchor_slot: 'features',
          layer: 'foreground',
          cover_anchor: true,
          style_override: nestedStyleOverride
        }
      }
    ]
  }
});
const circularP3Canonical = buildCanonicalRenderRequest({
  state: circularP3State,
  filesData
});
const circularP3Slots = circularP3Canonical.renderRequest.diagramOptions.tracks.circularTrackSlots;
assert.equal(circularP3Slots.find((slot) => slot.id === 'space').renderer, 'spacer');
assert.deepEqual(
  circularP3Slots.find((slot) => slot.id === 'review_overlay').params,
  {
    set_id: 'review',
    marks: ['line', 'band'],
    lane_gap_px: 5,
    padding_px: 0,
    anchor_slot: 'features',
    cover_anchor: true,
    style_override: nestedStyleOverride
  }
);
assert.deepEqual(
  projectCanonicalSessionRequest(circularP3Canonical)
    .config.adv.circular_track_slots.find((slot) => slot.id === 'review_overlay')
    .params.style_override,
  nestedStyleOverride
);

const linearP3State = {
  ...state,
  mode: ref('linear'),
  lInputType: ref('gb'),
  form: {
    ...state.form,
    show_depth: false,
    linear_track_layout: 'middle'
  },
  adv: {
    ...state.adv,
    linear_track_slots_enabled: true,
    linear_track_slots_axis_index: 0,
    linear_track_slots: [
      {
        id: 'features',
        renderer: 'features',
        enabled: true,
        side: 'overlay',
        height: '',
        spacing: '',
        z: 0,
        params: {}
      },
      {
        id: 'space',
        renderer: 'spacer',
        enabled: true,
        side: 'below',
        height: '14px',
        spacing: '3px',
        z: 0,
        params: {}
      },
      {
        id: 'review_overlay',
        renderer: 'annotations',
        enabled: true,
        side: 'overlay',
        height: '',
        spacing: '',
        z: 2,
        params: {
          set_id: 'review',
          marks: ['bracket'],
          lane_gap_px: 6,
          padding_px: 1,
          anchor_slot: 'features',
          layer: 'foreground',
          cover_anchor: true,
          style_override: nestedStyleOverride
        }
      }
    ]
  },
  annotationSets: structuredClone(p3AnnotationSet)
};
const linearP3Canonical = buildCanonicalRenderRequest({
  state: linearP3State,
  filesData: linearFilesData
});
const linearP3Slots = linearP3Canonical.renderRequest.diagramOptions.tracks.linearTrackSlots;
assert.deepEqual(
  linearP3Slots.find((slot) => slot.id === 'space'),
  buildLinearTrackSlotPayload(linearP3State.adv.linear_track_slots[1])
);
assert.deepEqual(
  linearP3Slots.find((slot) => slot.id === 'review_overlay').params.style_override,
  nestedStyleOverride
);

const defaultAnnotationPayload = buildCircularTrackSlotPayload({
  id: 'default_annotation',
  renderer: 'annotations',
  enabled: true,
  side: 'outside',
  width: null,
  radius: null,
  inner_gap_px: null,
  outer_gap_px: null,
  z: 0,
  params: {
    set_id: 'review',
    marks: [],
    lane_gap_px: 3,
    padding_px: 2,
    overflow: 'error',
    show_labels: true,
    layer: 'foreground',
    cover_anchor: false
  }
});
assert.deepEqual(defaultAnnotationPayload.params, { set_id: 'review' });

const implicitAnnotationsWithGeometry = buildCanonicalRenderRequest({
  state: circularRequestState({
    annotationSets: [{
      ...p3AnnotationSet[0],
      annotations: [
        ...p3AnnotationSet[0].annotations,
        {
          ...p3AnnotationSet[0].annotations[0],
          id: 'highlight',
          mark: 'highlight'
        }
      ]
    }],
    adv: { feature_width_circular: 16 }
  }),
  filesData
});
assert.deepEqual(
  implicitAnnotationsWithGeometry.renderRequest.diagramOptions.tracks.circularTrackSlots
    .slice(0, 2)
    .map((slot) => [slot.id, slot.renderer, slot.params.marks]),
  [
    ['annotations_1', 'annotations', ['line']],
    ['annotations_1_highlight', 'annotations', ['highlight']]
  ]
);

const combinedGeometryState = circularRequestState({
  form: { show_depth: true },
  annotationSets: [{
    ...p3AnnotationSet[0],
    annotations: [
      ...p3AnnotationSet[0].annotations,
      {
        ...p3AnnotationSet[0].annotations[0],
        id: 'highlight',
        mark: 'highlight'
      }
    ]
  }],
  adv: {
    depth_tracks: [
      { label: 'Sample A', color: '#112233' },
      { label: 'Sample B', color: '#445566' }
    ],
    feature_width_circular: 16,
    depth_width_circular: 10,
    gc_content_width_circular: 18,
    gc_content_radius_circular: 0.66,
    gc_skew_width_circular: 20,
    gc_skew_radius_circular: 0.48
  }
});
combinedGeometryState.circularConservation = {
  enabled: true,
  source: 'upload',
  reference: 'subject',
  labels: 'Comparison',
  series: [{ sourceIndex: 0, label: 'Comparison', color: '#E15759' }],
  ring_width: 14,
  ring_gap: 3
};
const combinedGeometryCanonical = buildCanonicalRenderRequest({
  state: combinedGeometryState,
  filesData: {
    ...filesData,
    c_depth: [depthA, depthB],
    c_conservation_blasts: [blastTable],
    c_conservation_sequence_sources: [comparisonFasta]
  }
});
const combinedGeometrySlots =
  combinedGeometryCanonical.renderRequest.diagramOptions.tracks.circularTrackSlots;
assert.deepEqual(
  combinedGeometrySlots.map((slot) => [slot.id, slot.renderer]),
  [
    ['annotations_1', 'annotations'],
    ['annotations_1_highlight', 'annotations'],
    ['features', 'features'],
    ['ticks', 'ticks'],
    ['depth_1', 'depth'],
    ['depth_2', 'depth'],
    ['gc_content', 'dinucleotide_content'],
    ['gc_skew', 'dinucleotide_skew']
  ],
  'geometry materialization must retain the full implicit stack and automatic annotations'
);
assert.equal(
  combinedGeometryCanonical.renderRequest.diagramOptions.tracks.circularTrackAxisIndex,
  2
);
assert.equal(
  combinedGeometrySlots.filter((slot) => slot.renderer === 'depth').length,
  2
);
assert.deepEqual(
  combinedGeometrySlots
    .filter((slot) => slot.renderer === 'depth')
    .map((slot) => slot.width),
  [
    { value: 10, unit: 'px' },
    { value: 10, unit: 'px' }
  ]
);
assert.equal(
  combinedGeometryCanonical.renderRequest.diagramOptions.conservationBlastFiles.length,
  1,
  'materializing geometry must retain managed comparison input'
);

const [managedComparisonSource] = orderedConservationSources(
  [blastTable],
  combinedGeometryState.circularConservation
);
const managedComparisonState = circularRequestState({
  adv: {
    circular_track_slots_enabled: true,
    circular_track_slots_axis_index: 0,
    circular_track_slots: [
      {
        id: 'features',
        renderer: 'features',
        enabled: true,
        side: 'inside',
        z: 0,
        params: { lane_direction: 'inside' }
      },
      {
        id: 'comparison',
        renderer: 'sequence_conservation',
        enabled: true,
        side: 'inside',
        z: 0,
        params: {
          managed: 'circular_conservation',
          series_key: managedComparisonSource.sourceKey,
          source_index: managedComparisonSource.orderIndex,
          track_index: managedComparisonSource.orderIndex + 1
        }
      }
    ]
  }
});
managedComparisonState.circularConservation = structuredClone(
  combinedGeometryState.circularConservation
);
const managedComparisonCanonical = buildCanonicalRenderRequest({
  state: managedComparisonState,
  filesData: {
    ...filesData,
    c_conservation_blasts: [blastTable],
    c_conservation_sequence_sources: [comparisonFasta]
  }
});
assert.equal(
  managedComparisonCanonical.renderRequest.diagramOptions.tracks.circularTrackSlots
    .filter((slot) => slot.renderer === 'sequence_conservation').length,
  1
);
assert.equal(
  managedComparisonCanonical.renderRequest.diagramOptions.conservationBlastFiles.length,
  1
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
    webFiles: session.webFiles,
    legacyFiles: session.files,
    storedConfig: session.config,
    fileBindings: session.cliInvocation?.fileBindings,
    linearTrackSlotSchemaVersion: Number(session.version) <= 32 ? 1 : 2
  });
  assert.ok(['circular', 'linear'].includes(projectedSession.mode));
  assert.ok(
    projectedSession.files.c_gb || projectedSession.files.linearSeqs.length > 0,
    'projected session must restore at least one input record'
  );
  if (sessionPath.includes('hepatoplasmataceae_')) {
    assert.deepEqual(
      projectedSession.files.linearSeqs.map((sequence) => sequence.gb.name),
      ['AP027078.gb', 'AP027131.gb', 'AP027133.gb', 'AP027132.gb', 'NZ_CP006932.gb']
    );
    assert.equal(projectedSession.config.filterMode, 'None');
  }
  if (sessionPath.includes('majanivirus_orthogroup')) {
    assert.equal(projectedSession.files.d_color.name, 'modified_default_colors.tsv');
    assert.equal(projectedSession.files.t_color.name, 'majani_custom_color_table.tsv');
    assert.equal(projectedSession.config.adv.def_font_size, 18);
    assert.equal(projectedSession.config.adv.block_stroke_width, 1);
    assert.equal(projectedSession.config.adv.line_stroke_width, 2);
    assert.equal(projectedSession.config.adv.legend_box_size, 18);
    assert.equal(projectedSession.config.adv.legend_font_size, 18);
  }
  if (sessionPath.includes('tobacco-chloroplast')) {
    assert.equal(projectedSession.files.t_color.name, 'chloroplast_specific_table.tsv');
    assert.equal(projectedSession.files.qualifier_priority.name, 'qualifier_priority.tsv');
    assert.equal(projectedSession.config.adv.def_font_size, 28);
    assert.equal(projectedSession.config.adv.block_stroke_width, 1);
    assert.equal(projectedSession.config.adv.block_stroke_color, 'black');
    assert.equal(projectedSession.config.adv.line_stroke_width, 2);
    assert.equal(projectedSession.config.adv.axis_stroke_width, 3);
    assert.equal(projectedSession.config.adv.outer_label_x_offset, 0.9);
    assert.equal(projectedSession.config.adv.outer_label_y_offset, 0.9);
    assert.equal(projectedSession.config.adv.inner_label_x_offset, 0.975);
    assert.equal(projectedSession.config.adv.inner_label_y_offset, 0.975);
  }
  if (sessionPath.includes('vibrio-harveyi-group-collinear')) {
    assert.equal(projectedSession.files.linearSeqs.length, 11);
    assert.equal(
      projectedSession.files.linearSeqs[0].gb.name,
      'NZ_CP125875.1__GCF_030060435.1_ASM3006043v1_genomic.gbff'
    );
    assert.equal(projectedSession.config.adv.block_stroke_width, 0);
    assert.equal(projectedSession.config.adv.line_stroke_width, 1);
    assert.equal(projectedSession.config.adv.axis_stroke_width, 2);
    assert.equal(projectedSession.config.adv.def_font_size, 16);
    assert.equal(projectedSession.config.filterMode, 'None');
  }
  if (sessionPath.includes('python-canonical-depth')) {
    assert.deepEqual(
      projectedSession.files.linearSeqs.map((sequence) => sequence.depth.map(Boolean)),
      [[true, true], [true, false]]
    );
    assert.deepEqual(
      projectedSession.config.adv.depth_tracks,
      [
        {
          label: 'Shared',
          color: '#112233',
          height: 18,
          large_tick_interval: 10,
          small_tick_interval: null,
          tick_font_size: null
        },
        {
          label: 'Sparse',
          color: '#445566',
          height: 24,
          large_tick_interval: null,
          small_tick_interval: 5,
          tick_font_size: 9
        }
      ]
    );
  }
  if (sessionPath.includes('WSSV_genome_comparison')) {
    assert.equal(projectedSession.files.c_conservation_blasts.length, 20);
    assert.equal(projectedSession.files.c_conservation_blasts_source, 'losat-cache');
    assert.equal((projectedSession.files.c_conservation_fastas || []).length, 20);
  }
}

const roundTripSessionIndex = process.argv.indexOf('--round-trip-session');
if (roundTripSessionIndex >= 0) {
  const sessionPath = process.argv[roundTripSessionIndex + 1];
  if (!sessionPath) throw new Error('--round-trip-session requires a JSON path.');
  const session = JSON.parse(await readFile(sessionPath, 'utf8'));
  const sessionProjection = projectCanonicalSessionRequest({
    renderRequest: session.renderRequest,
    resources: session.resources,
    webFiles: session.webFiles || {},
    legacyFiles: session.files || null,
    storedConfig: session.config || null
  });
  console.log(JSON.stringify(buildCanonicalRenderRequest({
    state: stateForCanonicalProjection(sessionProjection),
    filesData: sessionProjection.files
  })));
} else if (process.argv.includes('--print-resolved-protein')) {
  console.log(JSON.stringify(resolvedProteinRoundTripCanonical));
} else if (process.argv.includes('--print-depth')) {
  console.log(JSON.stringify(sparseDepthCanonical));
} else if (process.argv.includes('--print')) {
  console.log(JSON.stringify(canonical));
}
