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
  MODE_PROFILE_STATE_SCHEMA,
  MODE_PROFILE_VERSION,
  comparisonFiltersForMode,
  comparisonProfileDefault,
  comparisonStateForMode,
  createModeProfileStateManager,
  effectiveLinearAxisColor,
  managedAdvStateForMode,
  modeProfile,
  trackDefaultsForMode
} = await import(pathToFileURL(join(tempDir, 'js', 'mode-profiles.js')));
const {
  WEB_UX_PROFILE,
  WEB_UX_PROFILE_VERSION
} = await import(pathToFileURL(join(tempDir, 'js', 'web-ux-profile.js')));

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
    gridByDefault: false,
    legend: 'left',
    plotTitlePosition: 'none'
  },
  linear: {
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

const managedState = managedAdvStateForMode('circular');
const manager = createModeProfileStateManager('circular', managedState);
manager.transition(managedState, 'circular', 'linear');
assert.deepEqual(managedState, managedAdvStateForMode('linear'));
managedState.identity = 70;
managedState.alignment_length = '';
manager.transition(managedState, 'linear', 'circular');
assert.deepEqual(managedState, managedAdvStateForMode('circular'));
manager.transition(managedState, 'circular', 'linear');
assert.equal(managedState.identity, 70);
assert.equal(
  managedState.alignment_length,
  '',
  'a blank customization must not be treated as the numeric default zero'
);
assert.equal(managedState.axis_stroke_color, 'lightgray');

const importedLinearState = {
  ...managedAdvStateForMode('linear'),
  identity: 91,
  axis_stroke_color: 'gray'
};
const importedManager = createModeProfileStateManager(
  'circular',
  managedAdvStateForMode('circular')
);
importedManager.invalidate('linear');
importedManager.transition(importedLinearState, 'linear', 'circular');
importedManager.transition(importedLinearState, 'circular', 'linear');
assert.equal(importedLinearState.identity, 91);
assert.equal(importedLinearState.axis_stroke_color, 'gray');

{
  const liveState = {
    ...managedAdvStateForMode('circular'),
    identity: 88,
    axis_stroke_color: '#123456'
  };
  const profileManager = createModeProfileStateManager('circular', liveState);
  profileManager.transition(liveState, 'circular', 'linear');
  liveState.evalue = '1e-9';
  liveState.identity = 77;

  const exported = profileManager.exportState();
  assert.equal(exported.schema, MODE_PROFILE_STATE_SCHEMA);
  assert.equal(exported.activeMode, 'linear');
  assert.deepEqual(Object.keys(exported.profiles).sort(), ['circular', 'linear']);
  assert.equal(exported.profiles.circular.values.identity, 88);
  assert.equal(exported.profiles.circular.values.axis_stroke_color, '#123456');
  assert.equal(exported.profiles.circular.managed.identity, false);
  assert.equal(exported.profiles.linear.values.evalue, '1e-9');
  assert.equal(exported.profiles.linear.values.identity, 77);
  assert.equal(exported.profiles.linear.managed.evalue, false);
  assert.equal(exported.profiles.linear.managed.identity, false);

  const restoredState = managedAdvStateForMode('circular');
  const restoredManager = createModeProfileStateManager('circular', restoredState);
  restoredManager.importState(exported, 'linear', restoredState);
  assert.equal(restoredState.evalue, '1e-9');
  assert.equal(restoredState.identity, 77);
  restoredState.identity = 66;
  restoredManager.transition(restoredState, 'linear', 'circular');
  assert.equal(restoredState.identity, 88);
  assert.equal(restoredState.axis_stroke_color, '#123456');
  restoredManager.transition(restoredState, 'circular', 'linear');
  assert.equal(restoredState.identity, 66);

  const roundTrip = restoredManager.exportState();
  assert.equal(roundTrip.profiles.linear.values.identity, 66);
  assert.equal(roundTrip.profiles.linear.managed.identity, false);
}

{
  const restoredState = {
    ...managedAdvStateForMode('linear'),
    identity: 91,
    axis_stroke_color: 'gray'
  };
  const restoredManager = createModeProfileStateManager('circular', {});
  restoredManager.importState(null, 'linear', restoredState);
  const exported = restoredManager.exportState();
  assert.equal(exported.activeMode, 'linear');
  assert.equal(exported.profiles.linear.values.identity, 91);
  assert.equal(exported.profiles.linear.managed.identity, false);
  assert.deepEqual(
    exported.profiles.circular.values,
    managedAdvStateForMode('circular')
  );
  assert.throws(
    () => restoredManager.importState({ schema: 999, profiles: {} }, 'linear', restoredState),
    /mode profile state schema/
  );
}

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
const { createDefaultAdv, createDefaultForm, state } = await import(
  pathToFileURL(join(tempDir, 'js', 'state.js'))
);
assert.equal(Object.keys(state.form).includes('legend'), false);
assert.equal(Object.keys(state.adv).includes('plot_title_position'), false);
state.form.legend = 'right';
state.adv.plot_title_position = 'top';
assert.deepEqual(state.layoutPreferences.circular.single, {
  legend: 'right',
  plotTitlePosition: 'top'
});
state.mode.value = 'linear';
state.form.legend = 'left';
state.adv.plot_title_position = 'center';
assert.deepEqual(state.layoutPreferences.linear, {
  legend: 'left',
  plotTitlePosition: 'center'
});
state.mode.value = 'circular';
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
  Object.assign(state.adv, createDefaultAdv('circular'));
  state.modeProfileStateManager.reset('circular', state.adv);
  state.adv.identity = 88;
  state.modeProfileStateManager.transition(state.adv, 'circular', 'linear');
  state.mode.value = 'linear';
  state.adv.identity = 77;
  state.form.show_scale = false;
  state.adv.circular_track_slots_enabled = true;
  state.adv.circular_track_slots_axis_index = 1;
  state.adv.circular_track_slots.splice(
    0,
    state.adv.circular_track_slots.length,
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
  state.adv.linear_track_slots_enabled = true;
  state.adv.linear_track_slots_axis_index = 1;
  state.adv.linear_track_slots.splice(
    0,
    state.adv.linear_track_slots.length,
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

  resetSettings(state);

  const resetAdvDefaults = createDefaultAdv('linear');
  const resetProfiles = state.modeProfileStateManager.exportState();
  assert.deepEqual(
    resetProfiles.profiles.circular.values,
    managedAdvStateForMode('circular')
  );
  assert.deepEqual(
    resetProfiles.profiles.linear.values,
    managedAdvStateForMode('linear')
  );
  assert.ok(Object.values(resetProfiles.profiles.circular.managed).every(Boolean));
  assert.ok(Object.values(resetProfiles.profiles.linear.managed).every(Boolean));
  assert.equal(state.adv.circular_track_slots_enabled, false);
  assert.equal(state.adv.linear_track_slots_enabled, false);
  assert.equal(state.form.show_scale, true);
  assert.deepEqual(
    state.adv.circular_track_slots,
    resetAdvDefaults.circular_track_slots
  );
  assert.deepEqual(
    state.adv.linear_track_slots,
    resetAdvDefaults.linear_track_slots
  );
  assert.equal(
    state.adv.circular_track_slots_axis_index,
    resetAdvDefaults.circular_track_slots_axis_index
  );
  assert.equal(
    state.adv.linear_track_slots_axis_index,
    resetAdvDefaults.linear_track_slots_axis_index
  );
  assert.equal(state.files.c_gb, retainedFile);
}

{
  const { applyConfigData, buildConfigData } = await import(
    pathToFileURL(join(tempDir, 'js', 'services', 'config.js'))
  );
  state.mode.value = 'circular';
  Object.assign(state.adv, createDefaultAdv('circular'));
  state.modeProfileStateManager.reset('circular', state.adv);
  state.adv.identity = 88;
  state.modeProfileStateManager.transition(state.adv, 'circular', 'linear');
  state.mode.value = 'linear';
  state.adv.identity = 77;
  state.form.show_scale = false;
  const savedConfig = structuredClone(buildConfigData());

  assert.equal(savedConfig.form.show_scale, false);
  assert.equal(savedConfig.modeProfiles.profiles.circular.values.identity, 88);
  assert.equal(savedConfig.modeProfiles.profiles.circular.managed.identity, false);
  assert.equal(savedConfig.modeProfiles.profiles.linear.values.identity, 77);
  assert.equal(savedConfig.modeProfiles.profiles.linear.managed.identity, false);

  Object.assign(state.adv, createDefaultAdv('linear'));
  state.modeProfileStateManager.reset('linear', state.adv);
  state.form.show_scale = true;
  applyConfigData(savedConfig);
  assert.equal(state.form.show_scale, false);
  assert.equal(state.adv.identity, 77);
  state.modeProfileStateManager.transition(state.adv, 'linear', 'circular');
  assert.equal(state.adv.identity, 88);

  const version39Config = structuredClone(savedConfig);
  delete version39Config.modeProfiles;
  state.mode.value = 'linear';
  Object.assign(state.adv, createDefaultAdv('linear'));
  state.modeProfileStateManager.reset('linear', state.adv);
  applyConfigData(version39Config);
  const migratedProfiles = state.modeProfileStateManager.exportState();
  assert.equal(migratedProfiles.profiles.linear.values.identity, 77);
  assert.equal(migratedProfiles.profiles.linear.managed.identity, false);
  assert.deepEqual(
    migratedProfiles.profiles.circular.values,
    managedAdvStateForMode('circular')
  );
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
  Object.assign(state.form, createDefaultForm());
  Object.assign(state.adv, createDefaultAdv(modeName));
  state.modeProfileStateManager.reset(modeName, state.adv);
  state.cInputType.value = 'gb';
  state.lInputType.value = 'gb';
  state.circularRecordList.value = [];
  state.linearRecordLayoutEnabled.value = false;

  const filesData = modeName === 'circular'
    ? { c_gb: genbank, linearSeqs: [] }
    : {
        linearSeqs: [{
          uid: 'record-1',
          gb: genbank,
          region_record_id: '',
          region_start: null,
          region_end: null,
          region_reverse: false
        }],
        linearComparisons: []
      };
  const canonical = buildCanonicalRenderRequest({ state, filesData });
  const options = canonical.renderRequest.diagramOptions;
  const expected = expectedModes[modeName];
  assert.equal(canonical.renderRequest.mode, modeName);
  assert.equal(
    canonical.renderRequest.grouping,
    modeName === 'circular'
      ? WEB_UX_PROFILE.circular.singleRecordGrouping
      : 'single'
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
