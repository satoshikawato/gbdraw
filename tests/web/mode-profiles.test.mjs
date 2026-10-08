import assert from 'node:assert/strict';
import { cp, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const semanticParity = JSON.parse(await readFile(
  join(repoRoot, 'tests', 'fixtures', 'mode_semantic_parity.json'),
  'utf8'
));
const expectedModes = semanticParity.modes;
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-mode-profiles-'));
await cp(join(repoRoot, 'gbdraw', 'web', 'js'), join(tempDir, 'js'), { recursive: true });
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}', 'utf8');

const {
  MODE_DEFAULT_FEATURE_TYPES,
  MODE_PROFILE_VERSION,
  comparisonFiltersForMode,
  comparisonProfileDefault,
  comparisonStateForMode,
  effectiveLinearAxisColor,
  managedAdvStateForMode,
  modeProfile,
  trackDefaultsForMode
} = await import(pathToFileURL(join(tempDir, 'js', 'mode-profiles.js')));
const {
  WEB_UX_PROFILE,
  WEB_UX_PROFILE_VERSION
} = await import(pathToFileURL(join(tempDir, 'js', 'web-ux-profile.js')));
const { resolveLinearComparisonPlan } = await import(
  pathToFileURL(join(tempDir, 'js', 'services', 'linear-comparisons.js'))
);

const normalizedComparison = (filters) => ({
  evalue: Number(filters.evalue),
  bitscore: Number(filters.bitscore),
  identity: Number(filters.identity),
  alignmentLength: Number(filters.alignment_length)
});

const differentLeafPaths = (left, right, prefix = '') => {
  const paths = new Set();
  const keys = new Set([...Object.keys(left), ...Object.keys(right)]);
  keys.forEach((key) => {
    const path = prefix ? `${prefix}.${key}` : key;
    const leftValue = left[key];
    const rightValue = right[key];
    if (
      leftValue && rightValue &&
      typeof leftValue === 'object' && typeof rightValue === 'object' &&
      !Array.isArray(leftValue) && !Array.isArray(rightValue)
    ) {
      differentLeafPaths(leftValue, rightValue, path).forEach((entry) => paths.add(entry));
    } else if (!Object.is(leftValue, rightValue)) {
      paths.add(path);
    }
  });
  return paths;
};

assert.equal(semanticParity.schema, 1);
assert.equal(MODE_PROFILE_VERSION, semanticParity.profileVersion);
assert.equal(WEB_UX_PROFILE_VERSION, 1);
assert.deepEqual(WEB_UX_PROFILE, {
  separateStrands: true,
  circular: {
    singleRecordGrouping: 'single',
    multiRecordGrouping: 'batch',
    gridByDefault: true,
    legend: 'left',
    plotTitlePosition: 'none'
  },
  linear: {
    arrangeInRowsByDefault: true,
    legend: 'bottom',
    plotTitlePosition: 'bottom'
  }
});
assert.deepEqual(
  [...MODE_DEFAULT_FEATURE_TYPES],
  semanticParity.featureTypes
);
assert.deepEqual(
  normalizedComparison(comparisonFiltersForMode('circular')),
  expectedModes.circular.comparison
);
assert.deepEqual(
  normalizedComparison({
    ...comparisonStateForMode('linear'),
    bitscore: comparisonStateForMode('linear').min_bitscore
  }),
  expectedModes.linear.comparison
);
assert.equal(
  comparisonProfileDefault('circular', 'identity'),
  expectedModes.circular.comparison.identity
);
assert.deepEqual(trackDefaultsForMode('circular'), expectedModes.circular.tracks);
assert.deepEqual(trackDefaultsForMode('linear'), expectedModes.linear.tracks);
assert.equal(
  modeProfile('linear').linearAxisColor,
  expectedModes.linear.linearAxisColor
);
assert.equal(
  modeProfile('circular').linearRulerAxisColor,
  expectedModes.circular.linearRulerAxisColor
);
assert.equal(
  modeProfile('linear').linearRulerAxisColor,
  expectedModes.linear.linearRulerAxisColor
);
assert.deepEqual(
  [...differentLeafPaths(expectedModes.circular, expectedModes.linear)].sort(),
  semanticParity.reviewedModeDifferences.map(({ path }) => path).sort()
);
assert.ok(semanticParity.reviewedModeDifferences.every(({ reason }) => reason.trim()));

assert.equal(
  effectiveLinearAxisColor(),
  expectedModes.linear.linearAxisColor
);
assert.equal(
  effectiveLinearAxisColor({ rulerOnAxis: true }),
  expectedModes.linear.linearRulerAxisColor
);
assert.equal(
  effectiveLinearAxisColor({
    axisColor: expectedModes.linear.linearAxisColor,
    rulerOnAxis: true,
    managed: true
  }),
  expectedModes.linear.linearRulerAxisColor
);
assert.equal(
  effectiveLinearAxisColor({ axisColor: '#123456', rulerOnAxis: true }),
  '#123456'
);

assert.throws(() => comparisonStateForMode('radial'), /Unsupported diagram mode/);
assert.throws(
  () => comparisonProfileDefault('linear', 'unknown'),
  /Unsupported comparison field/
);

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({
      get value() {
        return getter();
      }
    }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
const { createDefaultAdv, createDefaultForm, createDefaultLosat, state } = await import(
  pathToFileURL(join(tempDir, 'js', 'state.js'))
);
assert.equal(
  createDefaultLosat().blastp.collinearSearchScope,
  'adjacent',
  'fresh Collinear LOSATP must search adjacent record pairs by default'
);
// Each mode has its own drawing (PR-1). A drawing starts from its mode's
// profile, and its `form.legend` and `adv.plot_title_position` read and write
// its own mode's layout slot, whichever mode is shown.
const { circular: circularDrawing, linear: linearDrawing } = state.drawings;
assert.notEqual(circularDrawing, linearDrawing);
for (const drawing of [circularDrawing, linearDrawing]) {
  assert.equal(Object.keys(drawing.form).includes('legend'), false);
  assert.equal(Object.keys(drawing.adv).includes('plot_title_position'), false);
  assert.equal(drawing.form.multi_record_canvas, true);
  assert.equal(drawing.linearRecordLayoutEnabled.value, true);
}
for (const [drawing, modeName] of [[circularDrawing, 'circular'], [linearDrawing, 'linear']]) {
  const managed = managedAdvStateForMode(modeName);
  assert.deepEqual(
    Object.fromEntries(Object.keys(managed).map((key) => [key, drawing.adv[key]])),
    managed,
    `the ${modeName} drawing starts from its mode's profile`
  );
}
state.mode.value = 'linear';
circularDrawing.form.multi_record_canvas = false;
circularDrawing.form.legend = 'right';
circularDrawing.adv.plot_title_position = 'top';
assert.deepEqual(circularDrawing.layoutPreferences.circular.single, {
  legend: 'right',
  plotTitlePosition: 'top'
});
state.mode.value = 'circular';
linearDrawing.form.legend = 'left';
linearDrawing.adv.plot_title_position = 'center';
assert.deepEqual(linearDrawing.layoutPreferences.linear, {
  legend: 'left',
  plotTitlePosition: 'center'
});
assert.deepEqual(circularDrawing.layoutPreferences.circular.single, {
  legend: 'right',
  plotTitlePosition: 'top'
}, 'a Linear layout edit leaves the Circular drawing');
// A value set in one drawing stays out of the other, and a mode switch moves
// nothing between them.
circularDrawing.adv.identity = 88;
linearDrawing.adv.identity = 77;
state.mode.value = 'linear';
state.mode.value = 'circular';
assert.deepEqual([circularDrawing.adv.identity, linearDrawing.adv.identity], [88, 77]);
assert.equal(state.activeDrawing(), circularDrawing);
const formDefaults = createDefaultForm();
assert.deepEqual(
  {
    multi_record_canvas: formDefaults.multi_record_canvas,
    separate_strands: formDefaults.separate_strands,
    suppress_gc: formDefaults.suppress_gc,
    suppress_skew: formDefaults.suppress_skew,
    show_gc: formDefaults.show_gc,
    show_skew: formDefaults.show_skew,
    show_scale: formDefaults.show_scale
  },
  {
    multi_record_canvas: WEB_UX_PROFILE.circular.gridByDefault,
    separate_strands: WEB_UX_PROFILE.separateStrands,
    suppress_gc: false,
    suppress_skew: false,
    show_gc: false,
    show_skew: false,
    show_scale: true
  }
);
const circularAdv = createDefaultAdv();
assert.deepEqual(
  {
    evalue: circularAdv.evalue,
    identity: circularAdv.identity,
    features: circularAdv.features
  },
  {
    evalue: '1e-5',
    identity: 70,
    features: ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'misc_RNA', 'repeat_region']
  }
);
const linearAdv = createDefaultAdv('linear');
assert.deepEqual(
  {
    evalue: linearAdv.evalue,
    identity: linearAdv.identity,
    axis_stroke_color: linearAdv.axis_stroke_color
  },
  { evalue: '1e-2', identity: 0, axis_stroke_color: 'lightgray' }
);

{
  const { resetSettings } = await import(
    pathToFileURL(join(tempDir, 'js', 'services', 'reset.js'))
  );
  state.mode.value = 'circular';
  state.drawings.circular.adv.identity = 88;
  state.drawings.circular.form.plot_title = 'Circular title';
  state.mode.value = 'linear';
  state.activeDrawing().adv.identity = 77;
  state.activeDrawing().form.show_scale = false;
  state.activeDrawing().linearTypographyLinked.value = false;
  state.activeDrawing().adv.scale_font_size = 18;
  state.activeDrawing().adv.ruler_label_font_size = 11;
  state.activeDrawing().adv.circular_track_slots_enabled = true;
  state.activeDrawing().adv.circular_track_slots_axis_index = 1;
  state.activeDrawing().adv.circular_track_slots.splice(
    0,
    state.activeDrawing().adv.circular_track_slots.length,
    {
      id: 'custom_annotation',
      renderer: 'annotations',
      enabled: false,
      side: 'outside',
      params: {
        set_id: 'review',
        style_override: {
          stroke: '#123456',
          hatch: { angle: 45, spacing: 4 }
        }
      }
    }
  );
  state.activeDrawing().adv.linear_track_slots_enabled = true;
  state.activeDrawing().adv.linear_track_slots_axis_index = 1;
  state.activeDrawing().adv.linear_track_slots.splice(
    0,
    state.activeDrawing().adv.linear_track_slots.length,
    {
      id: 'custom_spacer',
      renderer: 'spacer',
      enabled: false,
      side: 'below',
      height: '19px',
      spacing: '4px',
      params: {}
    }
  );
  const retainedFile = { name: 'retained.gb' };
  state.files.c_gb = retainedFile;
  const retainedComparisonFile = { name: 'retained-comparison.tsv' };
  state.activeDrawing().linearComparisonPlan.mode = 'selected';
  state.activeDrawing().linearComparisonPlan.defaultSource = 'upload';
  state.activeDrawing().linearComparisonPlan.edges.splice(
    0,
    state.activeDrawing().linearComparisonPlan.edges.length,
    {
      id: 'generated-only',
      queryUid: 'a',
      subjectUid: 'b',
      included: true,
      fileActive: false,
      losatFilenameActive: false,
      source: 'losat',
      file: null,
      losatFilename: ''
    },
    {
      id: 'retained-file',
      queryUid: 'a',
      subjectUid: 'b',
      included: true,
      fileActive: true,
      losatFilenameActive: false,
      source: 'upload',
      file: retainedComparisonFile,
      losatFilename: ''
    },
    {
      id: 'retained-name',
      queryUid: 'b',
      subjectUid: 'c',
      included: true,
      fileActive: false,
      losatFilenameActive: true,
      source: 'losat',
      file: null,
      losatFilename: 'custom-subject.fna'
    }
  );
  state.activeDrawing().losat.blastp.collinearSearchScope = 'adjacent';
  state.activeDrawing().unmanagedConfigOverrides['objects.gc_content.percent_background_opacity'] = 0.42;

  resetSettings(state);
  assert.equal(state.activeDrawing().form.multi_record_canvas, true);
  assert.equal(state.activeDrawing().linearRecordLayoutEnabled.value, true);

  const resetAdvDefaults = createDefaultAdv('linear');
  // Reset returns both drawings to their own mode's defaults.
  for (const modeName of ['circular', 'linear']) {
    const drawing = state.drawings[modeName];
    const managed = managedAdvStateForMode(modeName);
    assert.deepEqual(Object.fromEntries(Object.keys(managed).map((key) => [key, drawing.adv[key]])), managed);
    assert.equal(drawing.form.plot_title, '');
  }
  assert.equal(state.activeDrawing().adv.circular_track_slots_enabled, false);
  assert.equal(state.activeDrawing().adv.linear_track_slots_enabled, false);
  assert.equal(state.activeDrawing().form.show_scale, true);
  assert.equal(state.activeDrawing().linearTypographyLinked.value, true);
  assert.equal(state.activeDrawing().adv.scale_font_size, null);
  assert.equal(state.activeDrawing().adv.ruler_label_font_size, null);
  assert.deepEqual(state.activeDrawing().unmanagedConfigOverrides, {});
  assert.equal(
    state.activeDrawing().losat.blastp.collinearSearchScope,
    'adjacent',
    'Reset must restore the fresh adjacent Collinear default'
  );
  assert.deepEqual(
    state.activeDrawing().adv.circular_track_slots,
    resetAdvDefaults.circular_track_slots
  );
  assert.deepEqual(
    state.activeDrawing().adv.linear_track_slots,
    resetAdvDefaults.linear_track_slots
  );
  assert.equal(
    state.activeDrawing().adv.circular_track_slots_axis_index,
    resetAdvDefaults.circular_track_slots_axis_index
  );
  assert.equal(
    state.activeDrawing().adv.linear_track_slots_axis_index,
    resetAdvDefaults.linear_track_slots_axis_index
  );
  assert.equal(state.files.c_gb, retainedFile);
  assert.equal(state.activeDrawing().linearComparisonPlan.mode, 'none');
  assert.equal(state.activeDrawing().linearComparisonPlan.defaultSource, 'losat');
  assert.deepEqual(
    state.activeDrawing().linearComparisonPlan.edges.map((edge) => edge.id),
    ['retained-file', 'retained-name']
  );
  assert.equal(state.activeDrawing().linearComparisonPlan.edges[0].file, retainedComparisonFile);
  assert.equal(state.activeDrawing().linearComparisonPlan.edges[0].source, 'upload');
  assert.equal(state.activeDrawing().linearComparisonPlan.edges[0].included, false);
  assert.equal(state.activeDrawing().linearComparisonPlan.edges[0].fileActive, false);
  assert.equal(state.activeDrawing().linearComparisonPlan.edges[0].losatFilenameActive, false);
  assert.equal(state.activeDrawing().linearComparisonPlan.edges[1].losatFilename, 'custom-subject.fna');
  assert.equal(state.activeDrawing().linearComparisonPlan.edges[1].source, 'losat');
  assert.equal(state.activeDrawing().linearComparisonPlan.edges[1].included, false);
  assert.equal(state.activeDrawing().linearComparisonPlan.edges[1].fileActive, false);
  assert.equal(state.activeDrawing().linearComparisonPlan.edges[1].losatFilenameActive, false);
}

{
  const { applyConfigData, buildConfigData } = await import(
    pathToFileURL(join(tempDir, 'js', 'services', 'config.js'))
  );
  state.mode.value = 'circular';
  Object.assign(state.activeDrawing().adv, createDefaultAdv('circular'));
  state.activeDrawing().adv.identity = 88;
  state.mode.value = 'linear';
  state.activeDrawing().adv.identity = 77;
  state.activeDrawing().form.show_scale = false;
  state.activeDrawing().losat.blastp.collinearSearchScope = 'adjacent';
  state.activeDrawing().unmanagedConfigOverrides['objects.blast_match.curve_tension'] = 0.25;
  const savedConfig = structuredClone(buildConfigData(state.activeDrawing()));

  // A drawing's configuration is its own mode's (Session 46 `modes.<mode>.config`).
  assert.equal(savedConfig.form.show_scale, false);
  assert.equal(Object.hasOwn(savedConfig, 'modeProfiles'), false);
  assert.equal(savedConfig.adv.identity, 77);
  assert.equal(buildConfigData(state.drawings.circular).adv.identity, 88);
  assert.equal(savedConfig.losat.blastp.collinearSearchScope, 'adjacent');
  assert.deepEqual(savedConfig.unmanagedConfigOverrides, {
    'objects.blast_match.curve_tension': 0.25
  });

  Object.assign(state.activeDrawing().adv, createDefaultAdv('linear'));
  state.activeDrawing().form.show_scale = true;
  state.activeDrawing().unmanagedConfigOverrides.stale = true;
  applyConfigData(state.activeDrawing(), savedConfig);
  assert.equal(state.activeDrawing().form.show_scale, false);
  assert.equal(state.activeDrawing().adv.identity, 77);
  assert.deepEqual(state.activeDrawing().unmanagedConfigOverrides, {
    'objects.blast_match.curve_tension': 0.25
  });
  assert.equal(
    state.activeDrawing().losat.blastp.collinearSearchScope,
    'adjacent',
    'an explicit saved-session adjacent scope must survive current-reader loading'
  );
  assert.equal(state.drawings.circular.adv.identity, 88, 'applying the Linear drawing leaves the Circular one');

  // Arrange in rows: omission takes the fresh default; explicit false is kept.
  for (const [layout, expected] of [[undefined, true], [{ rows: [] }, true], [{ enabled: false, rows: [] }, false]]) {
    const layoutConfig = structuredClone(savedConfig);
    if (layout === undefined) delete layoutConfig.linearRecordLayout;
    else layoutConfig.linearRecordLayout = layout;
    applyConfigData(state.activeDrawing(), layoutConfig);
    assert.equal(state.activeDrawing().linearRecordLayoutEnabled.value, expected, JSON.stringify(layout));
  }

  const cliProjectedNumericConfig = structuredClone(savedConfig);
  cliProjectedNumericConfig.adv.arrow_head_length_ratio = '1.25';
  cliProjectedNumericConfig.adv.arrow_shaft_width_ratio = '0.25';
  applyConfigData(state.activeDrawing(), cliProjectedNumericConfig);
  assert.equal(state.activeDrawing().adv.arrow_head_length_ratio, 1.25);
  assert.equal(state.activeDrawing().adv.arrow_shaft_width_ratio, 0.25);

  const cliProjectedAutoConfig = structuredClone(savedConfig);
  cliProjectedAutoConfig.adv.arrow_head_length_ratio = 'auto';
  cliProjectedAutoConfig.adv.arrow_shaft_width_ratio = '1';
  applyConfigData(state.activeDrawing(), cliProjectedAutoConfig);
  assert.equal(state.activeDrawing().adv.arrow_head_length_ratio, null);
  assert.equal(state.activeDrawing().adv.arrow_shaft_width_ratio, 1.0);

  // D-04: a Session saved before the min-1 field may hold 0 or less, which drew automatic.
  for (const [stored, expected] of [[0, null], [-5, null], [250, 250]]) {
    const scaleConfig = structuredClone(savedConfig);
    scaleConfig.adv.scale_interval = stored;
    applyConfigData(state.activeDrawing(), scaleConfig);
    assert.equal(state.activeDrawing().adv.scale_interval, expected, `scale_interval ${stored}`);
  }

  Object.keys(state.activeDrawing().unmanagedConfigOverrides).forEach((path) => {
    delete state.activeDrawing().unmanagedConfigOverrides[path];
  });

}

const { buildCanonicalRenderRequest } = await import(
  pathToFileURL(join(tempDir, 'js', 'services', 'session-request.js'))
);
const genbankText = `LOCUS       MODEPARITY                  4 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Mode semantic parity fixture.
ACCESSION   MODEPARITY
VERSION     MODEPARITY
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
  name: 'mode-parity.gb',
  type: 'text/plain',
  size: new TextEncoder().encode(genbankText).byteLength,
  lastModified: 0,
  data: btoa(genbankText)
};

for (const modeName of ['circular', 'linear']) {
  state.mode.value = modeName;
  Object.assign(state.activeDrawing().form, createDefaultForm());
  Object.assign(state.activeDrawing().adv, createDefaultAdv(modeName));
  state.cInputType.value = 'gb';
  state.lInputType.value = 'gb';
  state.circularRecordList.value = [];
  state.activeDrawing().linearRecordLayoutEnabled.value = false;
  if (modeName === 'linear') {
    state.activeDrawing().linearComparisonPlan.mode = 'adjacent';
    state.activeDrawing().linearComparisonPlan.defaultSource = 'losat';
    state.activeDrawing().linearComparisonPlan.edges.splice(0);
  }

  const filesData = modeName === 'circular'
    ? { c_gb: genbank, linearSeqs: [] }
    : {
        linearSeqs: ['record-1', 'record-2'].map((uid) => ({
          uid,
          gb: genbank,
          region_record_id: '',
          region_start: null,
          region_end: null,
          region_reverse: false
        })),
        linearComparisons: []
      };
  const comparisonPlanSnapshot = modeName === 'linear'
    ? resolveLinearComparisonPlan({
        plan: state.activeDrawing().linearComparisonPlan,
        sequences: filesData.linearSeqs,
        layout: [],
        losatProgram: state.activeDrawing().losatProgram.value,
        blastpMode: state.activeDrawing().losat.blastp.mode
      })
    : null;
  const canonical = buildCanonicalRenderRequest({
    state,
    drawing: state.activeDrawing(),
    filesData,
    comparisonPlanSnapshot
  });
  const options = canonical.renderRequest.diagramOptions;
  const expected = expectedModes[modeName];
  assert.equal(canonical.renderRequest.mode, modeName);
  assert.equal(
    canonical.renderRequest.grouping,
    modeName === 'circular' ? 'grid' : 'single'
  );
  assert.deepEqual({
    evalue: options.evalue,
    bitscore: options.bitscore,
    identity: options.identity,
    alignmentLength: options.alignmentLength
  }, expected.comparison);
  assert.deepEqual(options.selectedFeaturesSet, semanticParity.featureTypes);
  assert.equal(options.configOverrides['canvas.show_gc'], expected.tracks.gc);
  assert.equal(options.configOverrides['canvas.show_skew'], expected.tracks.skew);
  assert.equal(options.configOverrides['objects.scale.show'], true);
  if (modeName === 'linear') {
    assert.equal(
      options.configOverrides['objects.axis.linear.stroke_color'],
      expected.linearAxisColor
    );
  } else {
    assert.equal(
      Object.hasOwn(options.configOverrides, 'objects.axis.linear.stroke_color'),
      false
    );
  }
}

state.mode.value = 'linear';
Object.assign(state.activeDrawing().form, createDefaultForm());
Object.assign(state.activeDrawing().adv, createDefaultAdv('linear'));
Object.assign(state.activeDrawing().losat, createDefaultLosat());
state.activeDrawing().losatProgram.value = 'blastp';
state.activeDrawing().losat.blastp.mode = 'collinear';
state.activeDrawing().linearComparisonPlan.mode = 'adjacent';
state.activeDrawing().linearComparisonPlan.defaultSource = 'losat';
state.activeDrawing().linearComparisonPlan.edges.splice(0);
const defaultCollinearRecords = ['record-1', 'record-2', 'record-3'].map((uid) => ({
  uid,
  gb: genbank,
  region_record_id: '',
  region_start: null,
  region_end: null,
  region_reverse: false
}));
const defaultCollinearSnapshot = resolveLinearComparisonPlan({
  plan: state.activeDrawing().linearComparisonPlan,
  sequences: defaultCollinearRecords,
  layout: [],
  losatProgram: state.activeDrawing().losatProgram.value,
  blastpMode: state.activeDrawing().losat.blastp.mode
});
const defaultCollinearRequest = buildCanonicalRenderRequest({
  state,
  drawing: state.activeDrawing(),
  filesData: {
    linearSeqs: defaultCollinearRecords,
    linearComparisons: [],
    linearCanonicalComparisons: []
  },
  comparisonPlanSnapshot: defaultCollinearSnapshot
}).renderRequest;
const defaultCollinearComparison = defaultCollinearRequest.comparisons.find(
  (comparison) => comparison.kind === 'generatedProteinComparison'
);
assert.equal(defaultCollinearComparison.mode, 'collinear');
assert.equal(
  defaultCollinearComparison.settings.collinearitySearchScope,
  'adjacent',
  'the fresh adjacent scope must reach the canonical Collinear request'
);

{
  // D-04: a loaded Session whose interval is 0 or less regenerates with the
  // automatic interval, so the Source recipe never writes a value the CLI rejects.
  const { applyConfigData, buildConfigData } = await import(
    pathToFileURL(join(tempDir, 'js', 'services', 'config.js'))
  );
  const { buildSourceRecipe } = await import(pathToFileURL(join(tempDir, 'js', 'app', 'run-info.js')));
  state.mode.value = 'circular';
  Object.assign(state.activeDrawing().form, createDefaultForm());
  Object.assign(state.activeDrawing().adv, createDefaultAdv('circular'));
  state.cInputType.value = 'gb';
  for (const [stored, flag] of [[-5, null], [0, null], [250, '250']]) {
    const loaded = structuredClone(buildConfigData(state.activeDrawing()));
    loaded.adv.scale_interval = stored;
    applyConfigData(state.activeDrawing(), loaded);
    const canonical = buildCanonicalRenderRequest({
      state, drawing: state.activeDrawing(), filesData: { c_gb: genbank, linearSeqs: [] }, comparisonPlanSnapshot: null
    });
    const recipe = await buildSourceRecipe(canonical);
    assert.equal(recipe.available, true, recipe.unavailableReason);
    const index = recipe.args.indexOf('--scale_interval');
    assert.equal(index === -1 ? null : recipe.args[index + 1], flag, `stored scale_interval ${stored}`);
  }
}
