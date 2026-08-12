import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js', 'services');
const appSourceDir = join(repoRoot, 'gbdraw', 'web', 'js', 'app');
const webSourceDir = join(repoRoot, 'gbdraw', 'web', 'js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-history-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempDir, 'web-ux-profile.js'),
  await readFile(join(webSourceDir, 'web-ux-profile.js'), 'utf8'),
  'utf8'
);
await mkdir(join(tempDir, 'services'), { recursive: true });
await mkdir(join(tempDir, 'app'), { recursive: true });
for (const filename of [
  'canonical-comparisons.js',
  'history.js',
  'history-files.js',
  'history-snapshot.js',
  'json-clone.js',
  'svg-serialization.js'
]) {
  await writeFile(
    join(tempDir, 'services', filename),
    await readFile(join(sourceDir, filename), 'utf8'),
    'utf8'
  );
}
for (const filename of [
  'feature-selector.js',
  'feature-visibility.js',
  'layout-preferences.js',
  'plot-title-position.js'
]) {
  await writeFile(
    join(tempDir, 'app', filename),
    await readFile(join(appSourceDir, filename), 'utf8'),
    'utf8'
  );
}

const { createHistoryManager } = await import(pathToFileURL(join(tempDir, 'services', 'history.js')));
const { createHistoryFileStore } = await import(pathToFileURL(join(tempDir, 'services', 'history-files.js')));
const { createHistorySnapshotService } = await import(
  pathToFileURL(join(tempDir, 'services', 'history-snapshot.js'))
);

const ref = (value) => ({ value });
const makeFile = (name, size = 10) => ({ name, size, type: 'text/plain', lastModified: 1 });
const createLayoutPreferences = () => ({
  circular: {
    single: { legend: 'left', plotTitlePosition: 'none' },
    multi: { legend: null, plotTitlePosition: null }
  },
  linear: { legend: 'bottom', plotTitlePosition: 'bottom' }
});

{
  let value = 0;
  const applied = [];
  const buildState = async () => ({ value });
  const applyState = async (snapshot) => {
    value = snapshot.value;
    applied.push(value);
  };
  const history = createHistoryManager({
    buildIntent: buildState,
    applyIntent: applyState,
    buildCheckpoint: buildState,
    applyCheckpoint: applyState
  });

  await history.captureBaseline('initial');
  await history.runUndoable('Increment', () => {
    value = 1;
  });
  assert.equal(history.getUndoCount(), 1);
  assert.equal(history.canUndo(), true);
  assert.equal(history.undoLabel(), 'Increment');

  await history.undo();
  assert.equal(value, 0);
  assert.deepEqual(applied, [0]);
  assert.equal(history.canRedo(), true);

  await history.redo();
  assert.equal(value, 1);
  assert.equal(history.getUndoCount(), 1);
}

{
  let value = 0;
  const buildState = () => ({ value });
  const history = createHistoryManager({
    buildIntent: buildState,
    applyIntent: (snapshot) => {
      value = snapshot.value;
    },
    buildCheckpoint: buildState,
    applyCheckpoint: (snapshot) => {
      value = snapshot.value;
    }
  });

  await history.captureBaseline('initial');
  const reset = history.runUndoableCheckpoint('Synchronous reset', () => {
    value = 1;
  });
  assert.equal(value, 1);
  await reset;
  assert.equal(history.undoLabel(), 'Synchronous reset');
  await history.undo();
  assert.equal(value, 0);
}

{
  let setting = 0;
  let artifact = 'artifact-a';
  const history = createHistoryManager({
    buildIntent: () => ({ setting }),
    applyIntent: (intent) => {
      setting = intent.setting;
    },
    buildCheckpoint: () => ({ setting, artifact }),
    applyCheckpoint: (checkpoint) => {
      setting = checkpoint.setting;
      artifact = checkpoint.artifact;
    }
  });

  await history.captureBaseline('initial');
  const pendingEdit = await history.begin('Edit setting');
  setting = 1;
  await history.runUndoableCheckpoint('Generate', () => {
    artifact = 'artifact-b';
  });

  assert.equal(pendingEdit.closed, true);
  assert.equal(history.getUndoCount(), 2);
  assert.equal(history.undoLabel(), 'Generate');
  await history.undo();
  assert.equal(setting, 1);
  assert.equal(artifact, 'artifact-a');
  assert.equal(history.undoLabel(), 'Edit setting');
  await history.undo();
  assert.equal(setting, 0);
  assert.equal(artifact, 'artifact-a');
  await history.redo();
  await history.redo();
  assert.equal(setting, 1);
  assert.equal(artifact, 'artifact-b');
}

{
  let value = 0;
  const buildState = async () => ({ value });
  const applyState = async (snapshot) => {
    value = snapshot.value;
  };
  const history = createHistoryManager({
    buildIntent: buildState,
    applyIntent: applyState,
    buildCheckpoint: buildState,
    applyCheckpoint: applyState
  });

  await history.captureBaseline('initial');
  await history.runUndoable('No-op', () => {});
  assert.equal(history.getUndoCount(), 0);

  await history.runUndoable('Set one', () => {
    value = 1;
  });
  await history.undo();
  assert.equal(history.canRedo(), true);
  await history.runUndoable('Set two', () => {
    value = 2;
  });
  assert.equal(history.canRedo(), false);
}

{
  let value = 0;
  const buildState = async () => ({ value });
  const applyState = async (snapshot) => {
    value = snapshot.value;
  };
  const history = createHistoryManager({
    maxActions: 2,
    buildIntent: buildState,
    applyIntent: applyState,
    buildCheckpoint: buildState,
    applyCheckpoint: applyState
  });

  await history.captureBaseline('initial');
  for (const nextValue of [1, 2, 3]) {
    await history.runUndoable(`Set ${nextValue}`, () => {
      value = nextValue;
    });
  }
  assert.equal(history.getUndoCount(), 2);
  await history.undo();
  assert.equal(value, 2);
  await history.undo();
  assert.equal(value, 1);
  assert.equal(history.canUndo(), false);
}

{
  let value = 'small';
  const buildState = async () => ({ value });
  const applyState = async (snapshot) => {
    value = snapshot.value;
  };
  const history = createHistoryManager({
    maxBytes: 220,
    buildIntent: buildState,
    applyIntent: applyState,
    buildCheckpoint: buildState,
    applyCheckpoint: applyState
  });

  await history.captureBaseline('initial');
  await history.runUndoable('Large state', () => {
    value = 'x'.repeat(1000);
  });
  assert.equal(history.canUndo(), false);
  assert.equal(history.historyLimitMessage.value.length > 0, true);
}

{
  const fileStore = createHistoryFileStore();
  const file = makeFile('a.gb', 100);
  const first = fileStore.describeFile(file);
  const second = fileStore.describeFile(file);
  assert.equal(first.fileId, second.fileId);

  let currentFile = file;
  const buildState = async () => ({ file: fileStore.describeFile(currentFile) });
  const applyState = async (snapshot) => {
    currentFile = fileStore.restoreFile(snapshot.file);
  };
  const history = createHistoryManager({
    maxActions: 0,
    fileStore,
    buildIntent: buildState,
    applyIntent: applyState,
    buildCheckpoint: buildState,
    applyCheckpoint: applyState
  });
  await history.captureBaseline('initial');
  currentFile = null;
  await history.runUndoable('Remove file', () => {});
  assert.equal(fileStore.has(first.fileId), false);
}

{
  let value = 0;
  let history = null;
  const buildState = async () => ({ value });
  const applyState = async (snapshot) => {
    value = snapshot.value;
    await history.runUndoable('Nested restore edit', () => {
      value = snapshot.value;
    });
  };
  history = createHistoryManager({
    buildIntent: buildState,
    applyIntent: applyState,
    buildCheckpoint: buildState,
    applyCheckpoint: applyState
  });

  await history.captureBaseline('initial');
  await history.runUndoable('Set one', () => {
    value = 1;
  });
  await history.undo();
  assert.equal(value, 0);
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 1);
}

{
  let value = 0;
  let intentBuildCount = 0;
  let intentApplyCount = 0;
  let checkpointBuildCount = 0;
  let checkpointApplyCount = 0;
  let commandSizeComputations = 0;
  const history = createHistoryManager({
    buildIntent: async () => {
      intentBuildCount += 1;
      return { value };
    },
    applyIntent: async (snapshot) => {
      intentApplyCount += 1;
      value = snapshot.value;
    },
    buildCheckpoint: async () => {
      checkpointBuildCount += 1;
      return { value };
    },
    applyCheckpoint: async (snapshot) => {
      checkpointApplyCount += 1;
      value = snapshot.value;
    }
  });

  await history.captureBaseline('initial');
  await history.runUndoable('Snapshot set one', () => {
    value = 1;
  });
  const checkpointsAfterIntentEdit = checkpointBuildCount;

  await history.runUndoableCommand('Command set two', () => ({
    apply: () => {
      value = 2;
    },
    revert: () => {
      value = 1;
    },
    estimateBytes: () => {
      commandSizeComputations += 1;
      return 8;
    }
  }));

  assert.equal(value, 2);
  assert.equal(history.getUndoCount(), 2);
  assert.equal(history.undoLabel(), 'Command set two');
  assert.equal(checkpointBuildCount, checkpointsAfterIntentEdit);
  assert.equal(commandSizeComputations, 1);

  await history.undo();
  assert.equal(value, 1);
  assert.equal(checkpointApplyCount, 0);
  assert.equal(history.undoLabel(), 'Snapshot set one');

  await history.undo();
  assert.equal(value, 0);
  assert.equal(intentApplyCount, 1);

  await history.redo();
  assert.equal(value, 1);
  await history.redo();
  assert.equal(value, 2);
  assert.equal(intentApplyCount, 2);
  assert.equal(checkpointApplyCount, 0);
  assert.equal(intentBuildCount > 0, true);
  assert.equal(commandSizeComputations, 1);
}

{
  let value = 0;
  let allowApply = true;
  let allowRevert = true;
  const warnings = [];
  const originalWarn = console.warn;
  console.warn = (message) => warnings.push(String(message));
  const buildState = async () => ({ value });
  const applyState = async (snapshot) => {
    value = snapshot.value;
  };
  const history = createHistoryManager({
    buildIntent: buildState,
    applyIntent: applyState,
    buildCheckpoint: buildState,
    applyCheckpoint: applyState
  });

  try {
    await history.captureBaseline('initial');
    assert.equal(
      await history.runUndoableCommand('Rejected command', () => ({
        apply: () => false,
        revert: () => true
      })),
      false
    );
    assert.equal(history.getUndoCount(), 0);

    assert.equal(
      await history.runUndoableCommand('Conditional command', () => ({
        apply: () => {
          if (!allowApply) return false;
          value = 1;
          return true;
        },
        revert: () => {
          if (!allowRevert) return false;
          value = 0;
          return true;
        }
      })),
      true
    );
    assert.equal(history.getUndoCount(), 1);
    assert.equal(history.undoLabel(), 'Conditional command');

    allowRevert = false;
    assert.equal(await history.undo(), false);
    assert.equal(value, 1);
    assert.equal(history.getUndoCount(), 1);
    assert.equal(history.getRedoCount(), 0);

    allowRevert = true;
    assert.equal(await history.undo(), true);
    assert.equal(value, 0);
    assert.equal(history.getUndoCount(), 0);
    assert.equal(history.getRedoCount(), 1);
    assert.equal(history.redoLabel(), 'Conditional command');

    allowApply = false;
    assert.equal(await history.redo(), false);
    assert.equal(value, 0);
    assert.equal(history.getUndoCount(), 0);
    assert.equal(history.getRedoCount(), 1);

    allowApply = true;
    assert.equal(await history.redo(), true);
    assert.equal(value, 1);
    assert.equal(history.getUndoCount(), 1);
    assert.equal(history.getRedoCount(), 0);
  } finally {
    console.warn = originalWarn;
  }
  assert.equal(warnings.length, 3);
}

{
  const fileStore = createHistoryFileStore();
  const file = makeFile('restore.gb', 25);
  const proteinTable = makeFile('resolved-protein.tsv', 40);
  const orthogroupsJson = makeFile('orthogroups.json', 50);
  const collinearityJson = makeFile('collinearity.json', 60);
  const comparisonUpload = makeFile('selected-comparison.tsv', 70);
  const nestedStyleOverride = {
    stroke: '#123456',
    strokeWidth: 2,
    hatch: {
      angle: 45,
      spacing: 4,
      color: '#654321',
      width: 1,
      cross: true
    }
  };
  let modeProfiles = {
    schema: 1,
    activeMode: 'circular',
    profiles: {
      circular: { values: { identity: 88 }, managed: { identity: false } },
      linear: { values: { identity: 77 }, managed: { identity: false } }
    }
  };
  const state = {
    form: { prefix: 'before', show_scale: true },
    adv: {
      features: ['CDS'],
      feature_shapes: { CDS: 'arrow' },
      arrow_head_length_ratio: 1.5,
      arrow_shaft_width_ratio: 0.25,
      circular_track_slots_enabled: true,
      circular_track_slots_axis_index: 1,
      circular_track_slots: [
        {
          id: 'notes',
          renderer: 'annotations',
          enabled: true,
          side: 'outside',
          params: {
            set_id: 'review',
            lane_gap_px: 5,
            style_override: nestedStyleOverride
          }
        },
        {
          id: 'disabled_space',
          renderer: 'spacer',
          enabled: false,
          side: 'inside',
          width: '8px',
          params: {}
        }
      ],
      linear_track_slots_enabled: false,
      linear_track_slots_axis_index: 1,
      linear_track_slots: [
        {
          id: 'inactive_space',
          renderer: 'spacer',
          enabled: true,
          side: 'below',
          height: '12px',
          spacing: '3px',
          params: {}
        }
      ]
    },
    files: {
      c_gb: file,
      c_gff: null,
      c_fasta: null,
      c_depth: null,
      c_conservation_blasts: [],
      c_conservation_blasts_source: 'losat-cache',
      c_conservation_fastas: [null, makeFile('comparison-2.fa', 30)],
      linearCanonicalComparisons: [
        {
          kind: 'precomputedProteinComparison',
          encoding: 'canonicalTsv',
          queryRecordIndex: 0,
          subjectRecordIndex: 2,
          file: proteinTable
        },
        {
          kind: 'orthogroupResult',
          encoding: 'canonicalJson',
          file: orthogroupsJson
        },
        {
          kind: 'collinearityResult',
          encoding: 'canonicalJson',
          valueKind: 'blocks',
          file: collinearityJson
        },
        {
          kind: 'generatedProteinComparison',
          mode: 'none',
          pairs: [],
          settings: { alignOrthogroupFeature: 'feature-anchor' }
        }
      ],
      d_color: null,
      t_color: null,
      blacklist: null,
      whitelist: null,
      qualifier_priority: null
    },
    linearSeqs: [
      { uid: 'record-a', gb: null, gff: null, fasta: null, depth: null },
      { uid: 'record-b', gb: null, gff: null, fasta: null, depth: null }
    ],
    linearComparisonPlan: {
      mode: 'selected',
      defaultSource: 'losat',
      edges: [{
        id: 'edge-a-b',
        queryUid: 'record-a',
        subjectUid: 'record-b',
        included: true,
        fileActive: true,
        losatFilenameActive: true,
        source: 'upload',
        file: comparisonUpload,
        losatFilename: 'custom-query.fna'
      }]
    },
    results: ref([{ name: 'r1', content: '<svg id="a"></svg>' }]),
    selectedResultIndex: ref(0),
    mode: ref('circular'),
    cInputType: ref('gb'),
    lInputType: ref('gb'),
    downloadDpi: ref(300),
    canvasPadding: { top: 1, right: 2, bottom: 3, left: 4 },
    layoutPreferences: createLayoutPreferences(),
    extractedFeatures: ref([]),
    featureRecordIds: ref([]),
    selectedFeatureRecordIdx: ref(0),
    featureColorOverrides: { f1: { color: '#111111', caption: 'A' } },
    featureVisibilityManualRules: [],
    featureVisibilityOverrides: {},
    featureVisibilitySelectorCache: {},
    featureStrokeOverrides: {},
    labelTextFeatureOverrides: {},
    labelTextBulkOverrides: {},
    labelTextFeatureOverrideSources: {},
    labelVisibilityOverrides: {},
    labelOverrideContextKey: ref(''),
    orthogroups: ref([]),
    selectedOrthogroupId: ref(''),
    selectedOrthogroupAlignmentFeature: ref(''),
    orthogroupNameOverrides: {},
    orthogroupDescriptionOverrides: {},
    lastRunInfo: ref(null),
    pairwiseMatchFactors: ref({}),
    skipCaptureBaseConfig: ref(false),
    skipPositionReapply: ref(false),
    skipExtractOnSvgChange: ref(false)
  };
  const snapshots = createHistorySnapshotService({
    state,
    fileStore,
    nextTick: async () => {},
    buildConfigData: () => ({
      form: state.form,
      adv: state.adv,
      modeProfiles,
      linearComparisonPlan: {
        mode: state.linearComparisonPlan.mode,
        defaultSource: state.linearComparisonPlan.defaultSource,
        edges: state.linearComparisonPlan.edges.map(({ file: _file, ...edge }) => edge)
      }
    }),
    applyConfigData: (config) => {
      state.form = { ...config.form };
      state.adv = { ...config.adv };
      modeProfiles = structuredClone(config.modeProfiles);
      state.linearComparisonPlan.mode = config.linearComparisonPlan.mode;
      state.linearComparisonPlan.defaultSource = config.linearComparisonPlan.defaultSource;
      state.linearComparisonPlan.edges.splice(
        0,
        state.linearComparisonPlan.edges.length,
        ...structuredClone(config.linearComparisonPlan.edges).map((edge) => ({
          ...edge,
          file: null
        }))
      );
    }
  });

  const snapshot = await snapshots.buildArtifactCheckpoint();
  state.form.prefix = 'after';
  state.form.show_scale = false;
  state.adv.arrow_head_length_ratio = null;
  state.adv.arrow_shaft_width_ratio = 1.0;
  state.files.c_gb = null;
  state.files.c_conservation_blasts_source = null;
  state.files.linearCanonicalComparisons = [];
  state.linearComparisonPlan.mode = 'none';
  state.linearComparisonPlan.defaultSource = 'upload';
  state.linearComparisonPlan.edges.splice(0);
  state.results.value = [{ name: 'r2', content: '<svg id="b"></svg>' }];
  state.featureColorOverrides.f1.color = '#222222';
  modeProfiles.profiles.circular.values.identity = 70;

  assert.equal(snapshot.config.form.prefix, 'before');
  assert.equal(snapshot.files.c_gb.name, 'restore.gb');
  await snapshots.applyArtifactCheckpoint(snapshot);
  assert.equal(state.form.prefix, 'before');
  assert.equal(state.form.show_scale, true);
  assert.equal(state.adv.feature_shapes.CDS, 'arrow');
  assert.equal(state.adv.arrow_head_length_ratio, 1.5);
  assert.equal(state.adv.arrow_shaft_width_ratio, 0.25);
  assert.equal(state.files.c_gb.name, 'restore.gb');
  assert.equal(state.files.c_conservation_blasts_source, 'losat-cache');
  assert.equal(state.files.c_conservation_fastas.length, 2);
  assert.equal(state.files.c_conservation_fastas[0], null);
  assert.equal(state.files.c_conservation_fastas[1].name, 'comparison-2.fa');
  assert.equal(state.linearComparisonPlan.mode, 'selected');
  assert.equal(state.linearComparisonPlan.defaultSource, 'losat');
  assert.equal(state.linearComparisonPlan.edges[0].id, 'edge-a-b');
  assert.equal(state.linearComparisonPlan.edges[0].source, 'upload');
  assert.equal(state.linearComparisonPlan.edges[0].included, true);
  assert.equal(state.linearComparisonPlan.edges[0].fileActive, true);
  assert.equal(state.linearComparisonPlan.edges[0].losatFilenameActive, true);
  assert.equal(state.linearComparisonPlan.edges[0].losatFilename, 'custom-query.fna');
  assert.equal(state.linearComparisonPlan.edges[0].file.name, 'selected-comparison.tsv');
  assert.equal(
    state.files.linearCanonicalComparisons[0].file.name,
    'resolved-protein.tsv'
  );
  assert.equal(
    state.files.linearCanonicalComparisons[1].file.name,
    'orthogroups.json'
  );
  assert.equal(
    state.files.linearCanonicalComparisons[2].valueKind,
    'blocks'
  );
  assert.equal(
    state.files.linearCanonicalComparisons[2].file.name,
    'collinearity.json'
  );
  assert.equal(
    state.files.linearCanonicalComparisons[3].settings.alignOrthogroupFeature,
    'feature-anchor'
  );
  assert.equal(state.results.value[0].name, 'r1');
  assert.equal(state.featureColorOverrides.f1.color, '#111111');
  assert.equal(modeProfiles.profiles.circular.values.identity, 88);
  assert.equal(modeProfiles.profiles.circular.managed.identity, false);
  assert.equal(modeProfiles.profiles.linear.values.identity, 77);
  assert.deepEqual(state.featureVisibilityManualRules, []);
  assert.equal(state.canvasPadding.top, 1);

  const history = createHistoryManager({
    fileStore,
    buildIntent: snapshots.buildHistoryIntent,
    applyIntent: snapshots.applyHistoryIntent,
    buildCheckpoint: snapshots.buildArtifactCheckpoint,
    applyCheckpoint: snapshots.applyArtifactCheckpoint
  });
  await history.captureBaseline('P3 baseline');
  const stableResults = state.results.value;
  const stableExtractedFeatures = state.extractedFeatures.value;
  await history.runUndoable('Hide coordinate scale', () => {
    state.form.show_scale = false;
  });
  assert.equal(state.form.show_scale, false);
  await history.undo();
  assert.equal(state.form.show_scale, true);
  assert.equal(state.results.value, stableResults);
  assert.equal(state.extractedFeatures.value, stableExtractedFeatures);
  assert.equal(state.linearComparisonPlan.edges[0].file, comparisonUpload);
  const serializedAfterConfigUndo = await snapshots.buildArtifactCheckpoint();
  assert.equal(
    fileStore.restoreValue(serializedAfterConfigUndo.files.linearComparisons[0].file),
    comparisonUpload
  );
  await history.redo();
  assert.equal(state.form.show_scale, false);
  assert.equal(state.linearComparisonPlan.edges[0].file, comparisonUpload);
  await history.undo();
  assert.equal(state.form.show_scale, true);
  assert.equal(state.linearComparisonPlan.edges[0].file, comparisonUpload);

  const replacementComparisonUpload = makeFile('replacement-comparison.tsv', 39);
  await history.runUndoable('Replace comparison upload', () => {
    state.linearComparisonPlan.edges[0].file = replacementComparisonUpload;
  });
  assert.equal(state.linearComparisonPlan.edges[0].file, replacementComparisonUpload);
  await history.undo();
  assert.equal(state.linearComparisonPlan.edges[0].file, comparisonUpload);
  await history.redo();
  assert.equal(state.linearComparisonPlan.edges[0].file, replacementComparisonUpload);
  await history.undo();
  assert.equal(state.linearComparisonPlan.edges[0].file, comparisonUpload);

  await history.runUndoable('Edit arrow geometry', () => {
    state.adv.arrow_head_length_ratio = 2;
    state.adv.arrow_shaft_width_ratio = 0.75;
  });
  assert.equal(state.adv.arrow_head_length_ratio, 2);
  assert.equal(state.adv.arrow_shaft_width_ratio, 0.75);
  await history.undo();
  assert.equal(state.adv.arrow_head_length_ratio, 1.5);
  assert.equal(state.adv.arrow_shaft_width_ratio, 0.25);
  await history.redo();
  assert.equal(state.adv.arrow_head_length_ratio, 2);
  assert.equal(state.adv.arrow_shaft_width_ratio, 0.75);
  await history.undo();

  await history.runUndoable('Edit Annotation lane gap', () => {
    state.adv.circular_track_slots[0].params.lane_gap_px = 9;
  });
  assert.deepEqual(
    state.adv.circular_track_slots[0].params.style_override,
    nestedStyleOverride
  );
  await history.undo();
  assert.equal(state.adv.circular_track_slots[0].params.lane_gap_px, 5);
  assert.deepEqual(
    state.adv.circular_track_slots[0].params.style_override,
    nestedStyleOverride
  );
  await history.redo();
  assert.equal(state.adv.circular_track_slots[0].params.lane_gap_px, 9);
  assert.deepEqual(
    state.adv.circular_track_slots[0].params.style_override,
    nestedStyleOverride
  );

  const beforeCustomReset = structuredClone({
    slots: state.adv.circular_track_slots,
    axisIndex: state.adv.circular_track_slots_axis_index,
    linearSlots: state.adv.linear_track_slots,
    linearEnabled: state.adv.linear_track_slots_enabled
  });
  await history.runUndoable('Reset Circular Custom Tracks', () => {
    state.adv.circular_track_slots.splice(
      0,
      state.adv.circular_track_slots.length,
      {
        id: 'features',
        renderer: 'features',
        enabled: true,
        side: 'inside',
        params: { lane_direction: 'inside' }
      }
    );
    state.adv.circular_track_slots_axis_index = 0;
  });
  await history.undo();
  assert.deepEqual(state.adv.circular_track_slots, beforeCustomReset.slots);
  assert.equal(
    state.adv.circular_track_slots_axis_index,
    beforeCustomReset.axisIndex
  );
  assert.deepEqual(state.adv.linear_track_slots, beforeCustomReset.linearSlots);
  assert.equal(
    state.adv.linear_track_slots_enabled,
    beforeCustomReset.linearEnabled
  );
  assert.deepEqual(
    state.adv.circular_track_slots[0].params.style_override,
    nestedStyleOverride
  );
  await history.redo();
  assert.deepEqual(
    state.adv.circular_track_slots.map((slot) => slot.id),
    ['features']
  );
}

{
  const fileStore = createHistoryFileStore();
  let pendingModeReset = false;
  const mode = {
    _value: 'circular',
    get value() {
      return this._value;
    },
    set value(next) {
      if (next !== this._value) pendingModeReset = true;
      this._value = next;
    }
  };
  const state = {
    form: { prefix: '', legend: 'bottom' },
    adv: { features: ['CDS'], plot_title_position: 'bottom' },
    files: {
      c_gb: null,
      c_gff: null,
      c_fasta: null,
      c_depth: null,
      c_conservation_blasts: [],
      c_conservation_fastas: [],
      d_color: null,
      t_color: null,
      blacklist: null,
      whitelist: null,
      qualifier_priority: null
    },
    linearSeqs: [],
    results: ref([]),
    selectedResultIndex: ref(0),
    mode,
    cInputType: ref('gb'),
    lInputType: ref('gb'),
    downloadDpi: ref(300),
    canvasPadding: { top: 0, right: 0, bottom: 0, left: 0 },
    layoutPreferences: createLayoutPreferences(),
    extractedFeatures: ref([]),
    featureRecordIds: ref([]),
    selectedFeatureRecordIdx: ref(0),
    featureColorOverrides: {},
    featureVisibilityManualRules: [],
    featureVisibilityOverrides: {},
    featureVisibilitySelectorCache: {},
    featureStrokeOverrides: {},
    labelTextFeatureOverrides: {},
    labelTextBulkOverrides: {},
    labelTextFeatureOverrideSources: {},
    labelVisibilityOverrides: {},
    labelOverrideContextKey: ref(''),
    orthogroups: ref([]),
    selectedOrthogroupId: ref(''),
    selectedOrthogroupAlignmentFeature: ref(''),
    orthogroupNameOverrides: {},
    orthogroupDescriptionOverrides: {},
    lastRunInfo: ref(null),
    pairwiseMatchFactors: ref({}),
    skipCaptureBaseConfig: ref(false),
    skipPositionReapply: ref(false),
    skipExtractOnSvgChange: ref(false)
  };
  const snapshots = createHistorySnapshotService({
    state,
    fileStore,
    nextTick: async () => {
      if (!pendingModeReset) return;
      pendingModeReset = false;
      state.extractedFeatures.value = [];
      state.featureRecordIds.value = [];
      state.orthogroups.value = [];
      state.selectedOrthogroupId.value = '';
    },
    applyConfigData: (config) => {
      state.form = { ...state.form, ...config.form };
      state.adv = { ...state.adv, ...config.adv };
    }
  });

  await snapshots.applyArtifactCheckpoint({
    config: { form: { prefix: 'restored', legend: 'bottom' }, adv: { features: ['CDS'] } },
    ui: { mode: 'linear', lInputType: 'gb', selectedResultIndex: 0 },
    files: { linearSeqs: [] },
    results: [{ name: 'linear.svg', content: '<svg></svg>' }],
    features: {
      extractedFeatures: [{ id: 'feat-1', svg_id: 'f1', type: 'CDS' }],
      featureRecordIds: ['record-1'],
      selectedFeatureRecordIdx: 0
    },
    orthogroupState: {
      groups: [{ id: 'og_1', members: [] }],
      selectedOrthogroupId: 'og_1'
    },
    editorState: {},
    runState: {}
  });

  assert.equal(state.mode.value, 'linear');
  assert.equal(state.extractedFeatures.value.length, 1);
  assert.equal(state.featureRecordIds.value.length, 1);
  assert.equal(state.orthogroups.value.length, 1);
  assert.equal(state.selectedOrthogroupId.value, 'og_1');
}

{
  let setting = 0;
  let artifact = { id: 'artifact-a', svg: '<svg>' + 'x'.repeat(100_000) + '</svg>' };
  let checkpointBuilds = 0;
  const history = createHistoryManager({
    buildIntent: async () => ({ setting }),
    applyIntent: async (intent) => {
      setting = intent.setting;
    },
    buildCheckpoint: async () => {
      checkpointBuilds += 1;
      return { setting, artifact };
    },
    applyCheckpoint: async (checkpoint) => {
      setting = checkpoint.setting;
      artifact = checkpoint.artifact;
    }
  });

  await history.captureBaseline('large artifact baseline');
  const baseline = history.getDiagnostics();
  const originalArtifact = artifact;
  for (let index = 1; index <= 10; index += 1) {
    await history.runUndoable(`Setting ${index}`, () => {
      setting = index;
    });
  }
  const afterEdits = history.getDiagnostics();
  assert.equal(checkpointBuilds, 1);
  assert.equal(afterEdits.artifactCheckpointBuilds - baseline.artifactCheckpointBuilds, 0);
  assert.equal(afterEdits.signatureComputations - baseline.signatureComputations, 20);
  assert.equal(afterEdits.byteEstimateComputations - baseline.byteEstimateComputations, 10);
  assert.equal(afterEdits.historySvgBytes - baseline.historySvgBytes, 0);
  assert.equal(artifact, originalArtifact);

  const beforeNoop = history.getDiagnostics();
  await history.runUndoable('No-op setting', () => {});
  const afterNoop = history.getDiagnostics();
  assert.equal(history.getUndoCount(), 10);
  assert.equal(afterNoop.byteEstimateComputations, beforeNoop.byteEstimateComputations);

  await history.undo();
  assert.equal(setting, 9);
  assert.equal(artifact, originalArtifact);
  await history.redo();
  assert.equal(setting, 10);
  assert.equal(artifact, originalArtifact);

  const replacementArtifact = { id: 'artifact-b', svg: '<svg id="b" />' };
  await history.runUndoableCheckpoint('Generate', () => {
    artifact = replacementArtifact;
  });
  assert.equal(history.undoLabel(), 'Generate');
  await history.undo();
  assert.equal(artifact, originalArtifact);
  await history.redo();
  assert.equal(artifact, replacementArtifact);
}

{
  let setting = 0;
  let rejectRestore = false;
  const history = createHistoryManager({
    buildIntent: async () => ({ setting }),
    applyIntent: async (intent) => {
      if (rejectRestore) throw new Error('restore rejected');
      setting = intent.setting;
    },
    buildCheckpoint: async () => ({ setting }),
    applyCheckpoint: async (checkpoint) => {
      setting = checkpoint.setting;
    }
  });
  await history.captureBaseline();
  await history.runUndoable('Set one', () => {
    setting = 1;
  });
  rejectRestore = true;
  await assert.rejects(history.undo(), /restore rejected/);
  assert.equal(history.getUndoCount(), 1);
  assert.equal(history.getRedoCount(), 0);
}

{
  const fileStore = createHistoryFileStore();
  let forbiddenArtifactBuilds = 0;
  const state = {
    form: { setting: 1 },
    adv: {},
    files: {
      c_gb: null,
      c_gff: null,
      c_fasta: null,
      c_depth: null,
      c_conservation_blasts: [],
      c_conservation_blasts_source: 'losat-cache',
      c_conservation_fastas: [],
      c_conservation_sequence_sources: [makeFile('generated-source.fa', 9_000_000)],
      linearCanonicalComparisons: [{ kind: 'generated', pairs: ['x'.repeat(100_000)] }],
      d_color: null,
      t_color: null,
      blacklist: null,
      whitelist: null,
      qualifier_priority: null
    },
    linearSeqs: [],
    linearComparisonPlan: { mode: 'none', defaultSource: 'losat', edges: [] },
    results: ref([{ name: 'large.svg', content: '<svg>' + 'x'.repeat(100_000) + '</svg>' }]),
    selectedResultIndex: ref(0),
    mode: ref('circular'),
    cInputType: ref('gb'),
    lInputType: ref('gb'),
    downloadDpi: ref(300),
    canvasPadding: { top: 0, right: 0, bottom: 0, left: 0 },
    layoutPreferences: createLayoutPreferences(),
    extractedFeatures: ref([{ featureCatalog: 'x'.repeat(100_000) }]),
    selectedFeatureRecordIdx: ref(0),
    featureColorOverrides: {},
    featureVisibilityManualRules: [],
    featureVisibilityOverrides: {},
    featureStrokeOverrides: {},
    labelTextFeatureOverrides: {},
    labelTextBulkOverrides: {},
    labelTextFeatureOverrideSources: {},
    labelVisibilityOverrides: {},
    labelOverrideContextKey: ref(''),
    orthogroups: ref([{ members: ['x'.repeat(100_000)] }]),
    selectedOrthogroupId: ref(''),
    selectedOrthogroupAlignmentFeature: ref(''),
    orthogroupNameOverrides: {},
    orthogroupDescriptionOverrides: {},
    legendEntries: ref([]),
    deletedLegendEntries: ref([]),
    legendColorOverrides: {},
    legendStrokeOverrides: {},
    addedLegendCaptions: ref(new Set()),
    semanticFileWatchersSuppressed: ref(false)
  };
  const snapshots = createHistorySnapshotService({
    state,
    fileStore,
    buildConfigData: () => ({ form: state.form, adv: state.adv }),
    buildFeatureStateData: () => {
      forbiddenArtifactBuilds += 1;
      throw new Error('feature artifacts must not be built for intent');
    },
    buildEditorStateData: () => {
      forbiddenArtifactBuilds += 1;
      throw new Error('feature catalog must not be built for intent');
    },
    buildOrthogroupStateData: () => {
      forbiddenArtifactBuilds += 1;
      throw new Error('orthogroup artifacts must not be built for intent');
    },
    serializeResults: () => {
      forbiddenArtifactBuilds += 1;
      throw new Error('Results must not be serialized for intent');
    }
  });
  const intent = await snapshots.buildHistoryIntent();
  assert.equal(forbiddenArtifactBuilds, 0);
  assert.equal(Object.prototype.hasOwnProperty.call(intent, 'results'), false);
  assert.equal(Object.prototype.hasOwnProperty.call(intent, 'runState'), false);
  assert.equal(Object.prototype.hasOwnProperty.call(intent.features, 'extractedFeatures'), false);
  assert.equal(Object.prototype.hasOwnProperty.call(intent.editorState, 'featureCatalog'), false);
  assert.equal(Object.prototype.hasOwnProperty.call(intent.orthogroupState, 'groups'), false);
  assert.equal(Object.prototype.hasOwnProperty.call(intent.files, 'linearCanonicalComparisons'), false);
  assert.equal(Object.prototype.hasOwnProperty.call(intent.files, 'c_conservation_sequence_sources'), false);
  assert.equal(Object.prototype.hasOwnProperty.call(intent.files, 'c_conservation_blasts'), false);
}

{
  const fileStore = createHistoryFileStore();
  const state = {
    form: {},
    adv: {},
    files: {},
    linearSeqs: [],
    linearComparisonPlan: { mode: 'none', defaultSource: 'losat', edges: [] },
    results: ref([]),
    selectedResultIndex: ref(0),
    selectedAnnotation: ref({ setId: 'set-a', id: 'annotation-a' }),
    selectedSpecificPreset: ref('bakta'),
    newSpecRule: { feat: 'CDS', qual: 'product', val: 'before', color: '#112233', cap: 'Before' },
    newPriorityRule: { feat: 'CDS', order: 'product,gene' },
    newColorFeat: ref('gene'),
    newColorVal: ref('#123456'),
    newFeatureToAdd: ref('mobile_element'),
    newLegendCaption: ref('Draft legend'),
    newLegendColor: ref('#654321'),
    fileLegendCaptions: ref(new Set(['Imported legend'])),
    semanticFileWatchersSuppressed: ref(false)
  };
  const snapshots = createHistorySnapshotService({
    state,
    fileStore,
    buildUiStateData: () => ({}),
    applyUiStateData: () => {},
    serializeResults: () => [],
    applyResultsData: () => {}
  });
  const history = createHistoryManager({
    fileStore,
    buildIntent: snapshots.buildHistoryIntent,
    applyIntent: snapshots.applyHistoryIntent,
    buildCheckpoint: snapshots.buildArtifactCheckpoint,
    applyCheckpoint: snapshots.applyArtifactCheckpoint
  });
  const draftState = () => ({
    selectedAnnotation: structuredClone(state.selectedAnnotation.value),
    selectedSpecificPreset: state.selectedSpecificPreset.value,
    newSpecRule: structuredClone(state.newSpecRule),
    newPriorityRule: structuredClone(state.newPriorityRule),
    newColorFeat: state.newColorFeat.value,
    newColorVal: state.newColorVal.value,
    newFeatureToAdd: state.newFeatureToAdd.value,
    newLegendCaption: state.newLegendCaption.value,
    newLegendColor: state.newLegendColor.value,
    fileLegendCaptions: Array.from(state.fileLegendCaptions.value)
  });
  const before = draftState();

  await history.captureBaseline('draft baseline');
  await history.runUndoable('Edit rule drafts', () => {
    state.selectedAnnotation.value = { setId: 'set-b', id: 'annotation-b' };
    state.selectedSpecificPreset.value = 'pharokka';
    Object.assign(state.newSpecRule, { val: 'after', cap: 'After' });
    Object.assign(state.newPriorityRule, { feat: 'gene', order: 'gene,note' });
    state.newColorFeat.value = 'CDS';
    state.newColorVal.value = '#abcdef';
    state.newFeatureToAdd.value = 'repeat_region';
    state.newLegendCaption.value = 'Edited legend';
    state.newLegendColor.value = '#fedcba';
    state.fileLegendCaptions.value = new Set(['Edited imported legend']);
  });
  const after = draftState();

  assert.equal(history.getUndoCount(), 1);
  await history.undo();
  assert.deepEqual(draftState(), before);
  await history.redo();
  assert.deepEqual(draftState(), after);

  await history.runUndoableCheckpoint('Reset settings', () => {
    state.selectedAnnotation.value = null;
    state.selectedSpecificPreset.value = '';
    Object.assign(state.newSpecRule, { feat: 'CDS', qual: 'product', val: '', color: '#ff0000', cap: '' });
    Object.assign(state.newPriorityRule, { feat: 'CDS', order: 'product,gene,locus_tag' });
    state.newColorFeat.value = 'gene';
    state.newColorVal.value = '#d3d3d3';
    state.newFeatureToAdd.value = 'mobile_element';
    state.newLegendCaption.value = '';
    state.newLegendColor.value = '#808080';
    state.fileLegendCaptions.value = new Set();
  });
  await history.undo();
  assert.deepEqual(draftState(), after);
}

console.log('history tests passed');
