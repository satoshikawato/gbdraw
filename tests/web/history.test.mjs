import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js', 'services');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-history-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'services'), { recursive: true });
for (const filename of ['history.js', 'history-files.js', 'history-snapshot.js']) {
  await writeFile(
    join(tempDir, 'services', filename),
    await readFile(join(sourceDir, filename), 'utf8'),
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

{
  let value = 0;
  const applied = [];
  const history = createHistoryManager({
    buildSnapshot: async () => ({ value }),
    applySnapshot: async (snapshot) => {
      value = snapshot.value;
      applied.push(value);
    }
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
  const history = createHistoryManager({
    buildSnapshot: async () => ({ value }),
    applySnapshot: async (snapshot) => {
      value = snapshot.value;
    }
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
  const history = createHistoryManager({
    maxActions: 2,
    buildSnapshot: async () => ({ value }),
    applySnapshot: async (snapshot) => {
      value = snapshot.value;
    }
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
  const history = createHistoryManager({
    maxBytes: 220,
    buildSnapshot: async () => ({ value }),
    applySnapshot: async (snapshot) => {
      value = snapshot.value;
    }
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
  const history = createHistoryManager({
    maxActions: 0,
    fileStore,
    buildSnapshot: async () => ({ file: fileStore.describeFile(currentFile) }),
    applySnapshot: async (snapshot) => {
      currentFile = fileStore.restoreFile(snapshot.file);
    }
  });
  await history.captureBaseline('initial');
  currentFile = null;
  await history.runUndoable('Remove file', () => {});
  assert.equal(fileStore.has(first.fileId), false);
}

{
  let value = 0;
  let history = null;
  history = createHistoryManager({
    buildSnapshot: async () => ({ value }),
    applySnapshot: async (snapshot) => {
      value = snapshot.value;
      await history.runUndoable('Nested restore edit', () => {
        value = snapshot.value;
      });
    }
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
  const fileStore = createHistoryFileStore();
  const file = makeFile('restore.gb', 25);
  const state = {
    form: { prefix: 'before' },
    adv: { features: ['CDS'] },
    files: {
      c_gb: file,
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
    results: ref([{ name: 'r1', content: '<svg id="a"></svg>' }]),
    selectedResultIndex: ref(0),
    mode: ref('circular'),
    cInputType: ref('gb'),
    lInputType: ref('gb'),
    downloadDpi: ref(300),
    canvasPadding: { top: 1, right: 2, bottom: 3, left: 4 },
    extractedFeatures: ref([]),
    featureRecordIds: ref([]),
    selectedFeatureRecordIdx: ref(0),
    featureColorOverrides: { f1: { color: '#111111', caption: 'A' } },
    featureVisibilityOverrides: {},
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
    buildConfigData: () => ({ form: state.form, adv: state.adv }),
    applyConfigData: (config) => {
      state.form = { ...config.form };
      state.adv = { ...config.adv };
    }
  });

  const snapshot = await snapshots.buildHistorySnapshot();
  state.form.prefix = 'after';
  state.files.c_gb = null;
  state.results.value = [{ name: 'r2', content: '<svg id="b"></svg>' }];
  state.featureColorOverrides.f1.color = '#222222';

  assert.equal(snapshot.config.form.prefix, 'before');
  assert.equal(snapshot.files.c_gb.name, 'restore.gb');
  await snapshots.applyHistorySnapshot(snapshot);
  assert.equal(state.form.prefix, 'before');
  assert.equal(state.files.c_gb.name, 'restore.gb');
  assert.equal(state.results.value[0].name, 'r1');
  assert.equal(state.featureColorOverrides.f1.color, '#111111');
  assert.equal(state.canvasPadding.top, 1);
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
    extractedFeatures: ref([]),
    featureRecordIds: ref([]),
    selectedFeatureRecordIdx: ref(0),
    featureColorOverrides: {},
    featureVisibilityOverrides: {},
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

  await snapshots.applyHistorySnapshot({
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

console.log('history tests passed');
